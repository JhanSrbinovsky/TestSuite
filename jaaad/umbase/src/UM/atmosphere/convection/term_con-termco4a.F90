#if defined(A05_4A)  || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE TERM_CON-----------------------------------------------
!LL
!LL  PURPOSE : RETURENS A MASK FOR POINTS AT WHICH CONVECTION
!LL            IS TERMINATING
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
!LL  VERSION NO. 1
!LL
!LL  LOGICAL COMPONENTS COVERED: P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE TERM_CON(NPNTS,NLEV,K,BTERM,BWKP1,FLXKP1,THEKP1,QEKP1, &
                          THPI,QPI,QSEKP1,DELTAK,EXPI,EKP14,EKP34,      &
                          NEW_TERMC,PSTAR,PK,PKP1,XSBMIN)

      Use cv_run_mod, Only:                                             &
          qstice

      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! MODEL CONSTANTS
!-----------------------------------------------------------------------
!
#include "c_r_cp.h"
#include "mparfl.h"
#include "c_lheat.h"
#include "c_epslon.h"
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
      INTEGER NPNTS          ! IN VECTOR LENGTH
!
      INTEGER NLEV           ! IN NUMBER OF MODEL LAYER
!
      INTEGER K              ! IN PRESENT MODEL LAYER
!
      INTEGER I              ! LOOP COUNTER
!
!
!-----------------------------------------------------------------------
! VARIABLES THAT ARE INPUT
!-----------------------------------------------------------------------
!
      REAL THEKP1(NPNTS)     ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
!
      REAL QEKP1(NPNTS)      ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL QSEKP1(NPNTS)     ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL THPI(NPNTS)       ! IN INITIAL PARCEL POTENTIAL TEMPERATURE
                             !    (K)
!
      REAL QPI(NPNTS)        ! IN INITIAL PARCEL MIXING RATIO (KG/KG)
!
      REAL FLXKP1(NPNTS)     ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
!
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
!
      REAL DELTAK(NPNTS)     ! IN FORCED DETRAINMENT IN LAYER K
                             !    MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!
      REAL EXPI(NPNTS)       ! IN INITIAL PARCEL EXNER PRESSURE
!
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT RATE FOR LEVEL K+1/4
                             !    MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!
      REAL EKP34(NPNTS)      ! IN ENTRAINMENT RATE FOR LEVEL K+3/4
                             !    MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!
      REAL PSTAR(NPNTS)      ! IN SURFACE PRESSURE (PA)
!
      REAL PK(NPNTS)         ! Pressure at mid-point of level K
!
      REAL PKP1(NPNTS)       ! Pressure at mid_point of K+1
!
      REAL XSBMIN(NPNTS)     ! IN THRESHOLD FOR TERMINAL DETRAINMENT
!                            !    Function of delta P
      INTEGER, INTENT(IN) :: NEW_TERMC      !FLAG FOR SIMPLIFIED
                                            !TERMINATION
!                                           !OF CONVECTION
!
!-----------------------------------------------------------------------
! VARIABLES THAT ARE OUTPUT
!-----------------------------------------------------------------------
!
      LOGICAL BTERM(NPNTS)   ! OUT MASK OF THOSE POINTS AT WHICH
                             !     CONVECTION IS ENDING
!
!-----------------------------------------------------------------------
! VARIABLES THAT ARE DEFINED LOCALLY
!-----------------------------------------------------------------------
!
      REAL EL                ! LATENT HEAT OF CONDENSATION OR
                             ! (CONDENSATION + FUSION) (J/KG)
!
      REAL FLXMIN            ! MINIMUM CONVECTIVE MASSFLUX BELOW
                             ! WHICH TERMINAL DETRAINMENT OCCURS
                             ! (PA/S)
!
      REAL THVUNDI           ! POTENTIAL TEMPERATURE OF AN
                             ! UNDILUTE PARCEL IN LAYER K+1
                             ! FROM THE STARTING LAYER OF
                             ! CONVECTION (K)
!
      REAL THVEKP1           ! VIRTUAL POTENTIAL TEMPERATURE
                             ! OF ENVIRONMENT IN LAYER K+1 (K)
!
!*---------------------------------------------------------------------
!
!----------------------------------------------------------------------
!  CALCULATE MINIMUM MASS FLUX BELOW WHICH CONVECTION IS TERMINATED
!----------------------------------------------------------------------
!
      DO 10 I=1,NPNTS
        FLXMIN = MPARFL*(1.+EKP14(I))*(1.+EKP34(I))*PSTAR(I)
!
!-----------------------------------------------------------------------
!   CREATE A VECTOR OF LATENT HEATS
!-----------------------------------------------------------------------
!
       IF (BWKP1(I)) THEN
          EL = LC
       ELSE
          EL = LC + LF
       ENDIF
!L
!L----------------------------------------------------------------------
!L  PARCELS ARE ONLY CONSIDERED FOR TERMINATION IF THEY ARE DETRAINING
!L  EXCEPT AT THE TOP MODEL LAYER, WHERE ALL CONVECTION TERMINATES
!L
!L  IF THE PARCEL HAS A POTENTIAL TEMPETURE GREATER THAN THE
!L  POTENTIAL TEMPERATURE OF AN UNDILUTE PARCEL FORM THE STARTING
!L  LAYER OF CONVECION IN LAYER K+1 THEN CONVECTION IS TERMINATED
!L
!L  UM DOCUMENTATION PAPER P27
!L  SECTION (7), EQUATION (32)
!L
!L  CONVECTION IS ALSO TERMINATED IF MASS FLUX IN LAYER K+1 IS LESS
!L  IS LESS THAN A MINIMUM VALUE
!L
!L  UM DOCUMENTATION PAPER P27
!L  SECTION (7), EQUATION (33)
!L----------------------------------------------------------------------
!L
       THVUNDI=( THPI(I) + (EL/(EXPI(I)*CP)) *(QPI(I) - QSEKP1(I))      &
     &         +((LC-EL)/(EXPI(I)*CP))*MAX(0.0,(QPI(I)-QSTICE))         &
     &         )*(1.+C_VIRTUAL*QSEKP1(I))
!

!
       THVEKP1 = (THEKP1(I)*(1.+C_VIRTUAL*QEKP1(I)) + XSBMIN(I))
!

       IF (.NOT. BTERM(I)) THEN
! Depending on whether option has been chosen in UMUI, either use
! original 4a termination condition or Martin Willett's simplified
! termination condition (new_termc=1)
         IF(NEW_TERMC  ==  1) THEN
           BTERM(I) = (FLXKP1(I)  <   FLXMIN) .OR. ((K+1)  ==  NLEV)
         ELSE
           BTERM(I) = (((FLXKP1(I)  <   FLXMIN) .OR.                    &
     &             (THVUNDI  <   THVEKP1))                              &
     &               .AND. (DELTAK(I) >  0.0)) .OR. (K+1)  ==  NLEV
         ENDIF
       ENDIF

!
  10  CONTINUE
!
      RETURN
      END SUBROUTINE TERM_CON
#endif
