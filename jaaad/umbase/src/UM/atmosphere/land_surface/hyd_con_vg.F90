#if defined(A08_7A)
! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2000, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HYD_CON_VG---------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE HYD_CON_VG(NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,THETAK,K     &
     &,                   DK_DTHK                                       &
! LOGICAL LTIMER
     &,LTIMER                                                           &
     &)

      IMPLICIT NONE
!
! Description:
!     Calculates the hydraulic conductivity using Van Genuchten curves
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
! Subroutine arguments:
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & NPNTS                                                            &
                        ! IN points in grid                               
     &,SOIL_PTS         ! IN Number of soil points.

!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & SOIL_INDEX(NPNTS)! IN Array of soil points.

      REAL                                                              &
     & B(NPNTS)                                                         &
                        ! IN Exponent in conductivity and soil water     
!                       !    suction fits.
     &,KS(NPNTS)                                                        &
                        ! IN The saturated hydraulic conductivity (kg/m2 
     &,THETAK(NPNTS)    ! IN Fractional saturation.

      LOGICAL LTIMER    ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :
      REAL                                                              &
     & K(NPNTS)                                                         &
                        ! OUT The hydraulic conductivity (kg/m2/s).      
     &,DK_DTHK(NPNTS)   ! OUT The rate of change of K with THETAK
!                       !     (kg/m2/s).

! Local scalars:
      REAL                                                              &
     & BRACKET(NPNTS)                                                   &
                        ! WORK 1-S^(b+1)                                 
     &,DBR_DTHK(NPNTS)                                                  &
                        ! WORK The rate of change of BRACKET             
!                       !       with THETAK.
     &,KRED(NPNTS)                                                      &
                        ! WORK KSAT*S^L_WAG (kg/m2/s).                   
     &,SDUM(NPNTS)      ! WORK Bounded THETAK value.

      INTEGER                                                           &
     & I,J              ! WORK Loop counter.

      REAL                                                              &
     & L_WAG            ! Exponent in the van Mualem / Van Genuchten
!                       ! fit to the hydraulic conduactivity curve.
      PARAMETER (L_WAG=0.5)

      REAL                                                              &
     & THETA_MIN                                                        &
                        ! Minimum value of THETAK for K calculation.     
     &,THETA_MAX        ! Maximum value of THETAK for K calculation.
      PARAMETER (THETA_MIN=0.05, THETA_MAX=0.95)


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('HYDCONVG',103)
      ENDIF

!CDIR NODEP
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)

        SDUM(I)=MAX(THETAK(I),THETA_MIN)
        SDUM(I)=MIN(SDUM(I),THETA_MAX)

        DK_DTHK(I)=0.0
        BRACKET(I)=1-SDUM(I)**(B(I)+1)
        KRED(I)=KS(I)*SDUM(I)**L_WAG

        K(I)=KRED(I)*(1-BRACKET(I)**(1.0/(B(I)+1)))**2

!----------------------------------------------------------------------
! To avoid blow-up of implicit increments approximate by piecewise
! linear functions
! (a) for THETA>THETA_MAX  (ensuring that K=KS at THETA=1)
! (b) for THETA<THETA_MIN  (ensuring that K=0 at THETA=THETA_MIN)
!----------------------------------------------------------------------
        IF (THETAK(I).LT.THETA_MIN) THEN
          DK_DTHK(I)=K(I)/THETA_MIN
          K(I)=K(I)+DK_DTHK(I)*(MAX(THETAK(I),0.0)-THETA_MIN)
        ELSEIF (THETAK(I).GT.THETA_MAX) THEN
          DK_DTHK(I)=(KS(I)-K(I))/(1.0-THETA_MAX)
          K(I)=K(I)+DK_DTHK(I)*(MIN(THETAK(I),1.0)-THETA_MAX)
        ELSE
          DBR_DTHK(I)=-(B(I)+1)*SDUM(I)**B(I)
          DK_DTHK(I)=L_WAG*K(I)/SDUM(I)                                 &
     &            -2*KRED(I)/(B(I)+1)                                   &
     &            *(1-BRACKET(I)**(1.0/(B(I)+1)))                       &
     &            *(BRACKET(I)**(-B(I)/(B(I)+1)))                       &
     &            *DBR_DTHK(I)
        ENDIF

        IF ((THETAK(I).GT.1.0).OR.(THETAK(I).LT.0.0)) THEN
          DK_DTHK(I)=0.0
          K(I)=MAX(K(I),0.0)
          K(I)=MIN(K(I),KS(I))
        ENDIF

      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('HYDCONVG',104)
      ENDIF

      RETURN
      END

#endif
