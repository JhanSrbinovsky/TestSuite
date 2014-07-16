
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE THP_DET------------------------------------------------
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL  5.4     6/8/02   New deck created for version 4A convection
!LL                   scheme.            Gill Martin
!LL  6.2   03/02/05  Added section 5A. R A Stratton
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
!LL  VERSION NO. 1
!LL
!LL  LOGICAL COMPONENTS COVERED : P27
!LL
!LL  PURPOSE : CALCULATES POTENTIAL TEMPERATURE OF THE
!LL            PARCEL IN LAYER K+1 AFTER FORCED DETRAINMENT
!LL            IN LAYER K
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
!LL                  SECTION (6), EQUATION (28)
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,       &
     &                    DQSKP1,BGMKP1,BCALC,PK,PKP1,XSBMIN)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! MODEL CONSTANTS
!-----------------------------------------------------------------------
!
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
      INTEGER NPNTS           ! IN VECTOR LENGTH
!
      INTEGER I               ! LOOP COUNTER
!
!
!-----------------------------------------------------------------------
! VARAIBLES WHICH ARE INPUT
!-----------------------------------------------------------------------
!
      REAL THEKP1(NPNTS)      ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                              !    IN LAYER K+1 (K)
!
      REAL QPKP1(NPNTS)       ! IN PARCEL MIXING RATIO IN LAYER K+1
                              !    (KG/KG)
!
      REAL QSEKP1(NPNTS)      ! IN ENVIRONMENT SATURATED MIXING RATIO
                              !    IN LAYER K+1 (KG/KG)
!
      REAL DQSKP1(NPNTS)      ! IN GRADIENT OF SATURATION MIXING RATIO
                              !    POTENTIAL TEMPERATURE FOR THE
                              !    ENVIRONMENT IN LAYER K+1 (KG/KG/K)
!
      REAL QEKP1(NPNTS)       ! IN ENVIRONMENT MIXING RATIO IN
                              !    LAYER K+1 (KG/KG)
!
      REAL PK(NPNTS)         ! Pressure at mid-point of level K
!
      REAL PKP1(NPNTS)       ! Pressure at mid_point of K+1
!
      REAL XSBMIN(NPNTS)     ! IN THRESHOLD FOR FORCED DETRAINMENT
!                            !    Function of delta P
!
      LOGICAL BGMKP1(NPNTS)   ! IN MASK FOR PARCELS WHICH ARE SATURATED
                              !    IN LAYER K+1
!
      LOGICAL BCALC(NPNTS)    ! IN MASK FOR PARCELS AT WHICH
                              !    CALCULATIONS OF THIS SUBROUTINE ARE
                              !    TO BE CARRIED OUT
!
!
!-----------------------------------------------------------------------
! VARAIBLES WHICH ARE OUTPUT
!-----------------------------------------------------------------------
!
      REAL THPKP1(NPNTS)      ! OUT PARCEL POTENTIAL TEMPERATURE
                              !     IN LAYER K+1 AFTER FORCED
                              !     DETRAINMENT (K)
!
!*---------------------------------------------------------------------
!L
!L---------------------------------------------------------------------
!L  NO SIGNIFICANT STRUCTURE
!L---------------------------------------------------------------------
!L
!
      DO 10 I=1,NPNTS
        IF (BCALC(I))THEN
         IF (BGMKP1(I)) THEN
           THPKP1(I) = THEKP1(I) +                                      &
     &                (C_VIRTUAL*THEKP1(I)*                             &
     &                           (QEKP1(I)-QSEKP1(I)) + XSBMIN(I))      &
     &               /( 1. + C_VIRTUAL*THEKP1(I)*DQSKP1(I) )
!
         ELSE
           THPKP1(I) = (THEKP1(I)*(1. + C_VIRTUAL*QEKP1(I))             &
     &                                                + XSBMIN(I))      &
     &                    /(1. + C_VIRTUAL*QPKP1(I))
         END IF
        END IF
  10  CONTINUE
!
      RETURN
      END SUBROUTINE THP_DET
