#if defined(A08_7A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE DARCY--------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE DARCY (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,SATHH,           &
     &                  STHU1,DZ1,STHU2,DZ2,WFLUX                       &
     &,                 DWFLUX_DSTHU1,DWFLUX_DSTHU2                     &
! LOGICAL LTIMER
     &,LTIMER                                                           &
     &)

      IMPLICIT NONE
!
! Description:
!     Calculates the Darcian fluxes between adjacent soil layers.
!                                                     (Cox, 6/95)
!
!
! Documentation : UM Documentation Paper 25
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.2     15/11/00  New deck.  M. Best
!  6.1  08/12/03  Add !CDIR NODEP to force vectorisation. R Barnes
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
     & B(NPNTS)                                                         &
                            ! IN Clapp-Hornberger exponent.
     &,DZ1                                                              &
                            ! IN Thickness of the upper layer (m).
     &,DZ2                                                              &
                            ! IN Thickness of the lower layer (m).
     &,KS(NPNTS)                                                        &
                            ! IN Saturated hydraulic conductivity
!                           !    (kg/m2/s).
     &,SATHH(NPNTS)                                                     &
                            ! IN Saturated soil water pressure (m).
     &,STHU1(NPNTS)                                                     &
                            ! IN Unfrozen soil moisture content of upper
!                           !    layer as a fraction of saturation.
!
     &,STHU2(NPNTS)         ! IN Unfrozen soil moisture content of lower
!                           !    layer as a fraction of saturation.
!
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

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL HYD_CON
      EXTERNAL TIMER

!-----------------------------------------------------------------------

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
     & THETA(NPNTS,2)                                                   &
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
     &,DPSI_DTH(NPNTS,2)    ! WORK The rate of change of PSI with
!                           !      THETA(1:2) (m).
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('DARCY   ',103)
      ENDIF

!-----------------------------------------------------------------------
! Calculate the fractional saturation of the layers
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        THETA(I,1)=STHU1(I)
        THETA(I,2)=STHU2(I)
      ENDDO

!-----------------------------------------------------------------------
! Calculate the soil water suction of the layers.
!-----------------------------------------------------------------------
      DO N=1,2
!CDIR NODEP
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          IF (THETA(I,N) <= 0.01) THEN  ! Prevent blow up for dry soil.
            PSI(I,N)=SATHH(I)/(0.01**B(I))
            DPSI_DTH(I,N)=0.0
          ELSE
            PSI(I,N)=SATHH(I)/(THETA(I,N)**B(I))
            DPSI_DTH(I,N)=-B(I)*PSI(I,N)/THETA(I,N)
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
! DEPENDS ON: hyd_con
      CALL HYD_CON (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,THETAK,K,DK_DTHK     &
     &,             LTIMER)

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
        CALL TIMER('DARCY   ',104)
      ENDIF

      RETURN
      END SUBROUTINE DARCY
#endif
