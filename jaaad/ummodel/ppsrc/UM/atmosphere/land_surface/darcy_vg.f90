
! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2000, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************
!    SUBROUTINE DARCY_VG-----------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE DARCY_VG(NPNTS,SOIL_PTS,SOIL_INDEX                     &
     &,                 BEXP,KS,SATHH1                                  &
     &,                 STHU1,DZ1,STHU2,DZ2,WFLUX                       &
     &,                 DWFLUX_DSTHU1,DWFLUX_DSTHU2                     &
! LOGICAL LTIMER
     &,LTIMER                                                           &
     &)

      IMPLICIT NONE
!
! Description:
!     Calculates the Darcian fluxes between adjacent soil layers
!     using Van Genuchten formulation
!
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : Martin Best
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  7.1     23/05/08  New routine
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Global variables:

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Number of gridpoints.                  
     &,SOIL_PTS             ! IN Number of soil points.


!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL                                                              &
     & BEXP(NPNTS)                                                      &
                            ! IN Clapp-Hornberger exponent for K.       
     &,DZ1                                                              &
                            ! IN Thickness of the upper layer (m).     
     &,DZ2                                                              &
                            ! IN Thickness of the lower layer (m).      
     &,KS(NPNTS)                                                        &
                            ! IN Saturated hydraulic conductivity        
!                           !    (kg/m2/s).                             
     &,SATHH1(NPNTS)                                                    &
                            ! IN Saturated soil water pressure (m).      
     &,STHU1(NPNTS)                                                     &
                            ! IN Unfrozen soil moisture content of upper 
!                           !    layer as a fraction of saturation.
!
     &,STHU2(NPNTS)         ! IN Unfrozen soil moisture content of lower
!                           !    layer as a fraction of saturation.

      LOGICAL LTIMER        ! Logical switch for TIMER diags


!   Array arguments with intent(OUT) :
      REAL                                                              &
     & WFLUX(NPNTS)                                                     &
                            ! OUT The flux of water between layers       
!                           !     (kg/m2/s).
     &,DWFLUX_DSTHU1(NPNTS)                                             &
                            ! OUT The rate of change of the explicit     
!                           !     flux with STHU1 (kg/m2/s).
     &,DWFLUX_DSTHU2(NPNTS) ! OUT The rate of change of the explicit
!                           !     flux with STHU2 (kg/m2/s).

! Local scalars:
      INTEGER                                                           &
     & I,J,N                ! WORK Loop counters.

      REAL                                                              &
     & DTHK_DTH1,DTHK_DTH2                                              &
                            ! WORK DTHETAK/DTHETA(1:2).                  
     &,DK_DTH1,DK_DTH2                                                  &
                            ! WORK DK/DTHETA(1:2) (kg/m2/s).             
     &,PD                   ! WORK Hydraulic potential difference (m).

! Local arrays:
      REAL                                                              &
     & SDUM(NPNTS,2)                                                    &
                            ! WORK bounded values of THETAK.             
     &,THETA(NPNTS,2)                                                   &
                            ! WORK Fractional saturation of the upper    
!                           !      and lower layer respectively.        
     &,THETAK(NPNTS)                                                    &
                            ! WORK Fractional saturation at the layer    
!                           !      boundary.
     &,K(NPNTS)                                                         &
                            ! WORK The hydraulic conductivity between    
!                           !      layers (kg/m2/s).
     &,DK_DTHK(NPNTS)                                                   &
                            ! WORK The rate of change of K with THETAK   
!                           !      (kg/m2/s).
     &,PSI(NPNTS,2)                                                     &
                            ! WORK The soil water suction of the upper   
!                           !      and lower layer respectively (m).
     &,DPSI_DTH(NPNTS,2)                                                &
                            ! WORK The rate of change of PSI with        
!                           !      THETA(1:2) (m).
     &,BRACKET(NPNTS,2)                                                 &
                            ! WORK 1-S^(-b-1)                            
     &,DBR_DTH(NPNTS,2)                                                 &
                            ! WORK The rate of change of BRACKET with    
!                           !      THETA(1:2) (m).
     &,B(NPNTS,2)                                                       &
                            ! WORK Clapp-Hornberger exponent for PSI.    
     &,SATHH(NPNTS,2)       ! WORK Saturated soil water pressure (m).

      REAL                                                              &
     & THETA_MIN                                                        &
                        ! Minimum value of THETAK for K calculation.     
     &,THETA_MAX        ! Maximum value of THETAK for K calculation.
      PARAMETER (THETA_MIN=0.05, THETA_MAX=0.95)

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('DARCYVG ',103)
      ENDIF

!-----------------------------------------------------------------------
! Calculate the fractional saturation of the layers
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        THETA(I,1)=STHU1(I)
        THETA(I,2)=STHU2(I)
        B(I,1)=BEXP(I)
        B(I,2)=BEXP(I)
        SATHH(I,1)=SATHH1(I)
        SATHH(I,2)=SATHH1(I)

        SDUM(I,1)=MAX(THETA(I,1),THETA_MIN)
        SDUM(I,1)=MIN(SDUM(I,1),THETA_MAX)
        SDUM(I,2)=MAX(THETA(I,2),THETA_MIN)
        SDUM(I,2)=MIN(SDUM(I,2),THETA_MAX)
      ENDDO

!-----------------------------------------------------------------------
! Calculate the soil water suction of the layers.
!-----------------------------------------------------------------------
      DO N=1,2
!CDIR NODEP
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          BRACKET(I,N)=-1+SDUM(I,N)**(-B(I,N)-1)
          PSI(I,N)=SATHH(I,N)*BRACKET(I,N)**(B(I,N)/(B(I,N)+1))

!----------------------------------------------------------------------
! To avoid blow-up of implicit increments approximate by piecewise
! linear functions
! (a) for THETA>THETA_MAX  (ensuring that PSI=0 at THETA=1)
! (b) for THETA<THETA_MIN  (extrapolating using DPSI_DTH at THETA_MIN)
!----------------------------------------------------------------------
          IF (THETA(I,N).GT.THETA_MAX) THEN
            DPSI_DTH(I,N)=-PSI(I,N)/(1.0-THETA_MAX)
            PSI(I,N)=PSI(I,N)                                           &
     &              +DPSI_DTH(I,N)*(MIN(THETA(I,N),1.0)-THETA_MAX)
          ELSE
            DBR_DTH(I,N)=(-B(I,N)-1)*SDUM(I,N)**(-B(I,N)-2)
            DPSI_DTH(I,N)=SATHH(I,N)*B(I,N)/(B(I,N)+1)                  &
     &                 *BRACKET(I,N)**(-1.0/(B(I,N)+1))                 &
     &                 *DBR_DTH(I,N)
            IF (THETA(I,N).LT.THETA_MIN) THEN
              PSI(I,N)=PSI(I,N)                                         &
     &                +DPSI_DTH(I,N)*(MAX(THETA(I,N),0.0)-THETA_MIN)
            ENDIF
          ENDIF

          IF ((THETA(I,N).GT.1.0).OR.(THETA(I,N).LT.0.0)) THEN
            DPSI_DTH(I,N)=0.0
            PSI(I,N)=MAX(PSI(I,N),0.0)
          ENDIF

        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Estimate the fractional saturation at the layer boundary by
! interpolating the soil moisture.
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        THETAK(I)=(DZ2*THETA(I,1)+DZ1*THETA(I,2))/(DZ2+DZ1)
      ENDDO
      DTHK_DTH1=DZ2/(DZ1+DZ2)
      DTHK_DTH2=DZ1/(DZ2+DZ1)

!-----------------------------------------------------------------------
! Calculate the hydraulic conductivities for transport between layers.
!-----------------------------------------------------------------------
! DEPENDS ON: hyd_con_vg
      CALL HYD_CON_VG(NPNTS,SOIL_PTS,SOIL_INDEX,BEXP,KS,THETAK,K          &
     &,             DK_DTHK,LTIMER)

!-----------------------------------------------------------------------
! Calculate the Darcian flux from the upper to the lower layer.
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        PD=(2.0*(PSI(I,2)-PSI(I,1))/(DZ2+DZ1)+1)
        WFLUX(I)=K(I)*PD

!-----------------------------------------------------------------------
! Calculate the rate of change of WFLUX with respect to the STHU1 and
! STHU2.
!-----------------------------------------------------------------------
        DK_DTH1=DK_DTHK(I)*DTHK_DTH1
        DK_DTH2=DK_DTHK(I)*DTHK_DTH2
        DWFLUX_DSTHU1(I)=DK_DTH1*PD-2*K(I)*DPSI_DTH(I,1)/(DZ1+DZ2)
        DWFLUX_DSTHU2(I)=DK_DTH2*PD+2*K(I)*DPSI_DTH(I,2)/(DZ1+DZ2)
      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('DARCYVG ',104)
      ENDIF

      RETURN
      END
