#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE TERMDD-------------------------------------------------
!LL
!LL  PURPOSE : CALCULATE WHETHER DOWNDRAUGHT IS ABLE TO CONTINUE
!LL
!LL            CALCULATE BUOYANCY
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL  VERSION NO. 4  DATED 5/2/92
!LL
!LL  SYSTEM TASK : P27
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE TERMDD (NPNTS,BDD_START,THDD_K,QDD_K,THE_K,QE_K,K,     &
                         B_DD_END,BDD_ON,PPN_MIX_DD)

      Use cv_run_mod, Only:                                             &
          dd_opt

      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! MODEL CONSTANTS USED IN THIS ROUTINE
!-----------------------------------------------------------------------
!
#include "c_epslon.h"
#include "ddkmdet.h"
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
      INTEGER NPNTS                ! IN VECTOR LENGTH
!
      INTEGER I                    ! LOOP COUNTER
!
      INTEGER K                    ! IN PRESENT MODEL LAYER
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!-----------------------------------------------------------------------
!
!
      REAL THDD_K(NPNTS)           ! IN MODEL POTENTIAL TEMPERATURE
                                   !    OF DOWNDRAUGHT AT LAYER K (K)
!
      REAL QDD_K(NPNTS)            ! IN MODEL MIXING RATIO OF
                                   !    DOWNDRAUGHT AT LAYER K
!
      REAL THE_K(NPNTS)            ! IN POTENTIAL TEMPERATURE OF
                                   !    ENVIRONMENTAL AIR IN LAYER K
!
      REAL QE_K(NPNTS)             ! IN MODEL MIXING RATIO AT LAYER K
!
      REAL PPN_MIX_DD(NPNTS)       ! IN PRECIP MIXING RATIO      
!
      LOGICAL BDD_START(NPNTS)     ! IN MASK FOR THOSE POINTS WHERE
                                   !    DOWNDRAUGHT MAY OCCUR IN
                                   !    LAYER K-1
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!-----------------------------------------------------------------------
!
      LOGICAL B_DD_END(NPNTS)      ! OUT MASK FOR THOSE POINTS WHERE
                                   !     DOWNDRAUGHT IS TERMINATING
!
      LOGICAL BDD_ON(NPNTS)        ! OUT MASK FOR THOSE POINTS WHERE
                                   !     DOWNDRAUGHT CONTINUES TO LAYER
                                   !     K-1 (AS BDD_START HERE)
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!-----------------------------------------------------------------------
!
      REAL BUOY1                   ! BUOYANCY OF PARCEL
!
      REAL THDD_V                  ! USED IN CALCULATION OF BUOYANCY
!
      REAL THE_V                   ! USED IN CALCULATION OF BUOYANCY
!
!
      IF (DD_OPT == 1) THEN
!-----------------------------------------------------------------------
! CHECK IF PARCEL STILL NEGATIVELY BUOYANT SUCH THAT DOWNDRAUGHT
! CAN CONTINUE TO NEXT LAYER
!-----------------------------------------------------------------------
!
       DO I=1,NPNTS
         THDD_V = THDD_K(I)*(1.0+C_VIRTUAL*QDD_K(I)) &
          /(1.+PPN_MIX_DD(I))
         THE_V = THE_K(I)*(1.0+C_VIRTUAL*QE_K(I))
         BUOY1 = THDD_V - THE_V
!
!-----------------------------------------------------------------------
! CALCULATE STATE OF DOWNDRAUGHT
!-----------------------------------------------------------------------
!
         IF (BDD_START(I) .AND. BUOY1 >  0.5) THEN
            BDD_ON(I) = .FALSE.
         ELSE IF (BUOY1 >  0.5 .OR. K == 2) THEN
            B_DD_END(I) = .TRUE.
         END IF
       END DO
!
      ELSE
!-----------------------------------------------------------------------
! CHECK IF PARCEL STILL NEGATIVELY BUOYANT SUCH THAT DOWNDRAUGHT
! CAN CONTINUE TO NEXT LAYER
!-----------------------------------------------------------------------
!
       DO I=1,NPNTS
         THDD_V = THDD_K(I)*(1.0+C_VIRTUAL*QDD_K(I))
         THE_V = THE_K(I)*(1.0+C_VIRTUAL*QE_K(I))
         BUOY1 = THDD_V - THE_V
!
!-----------------------------------------------------------------------
! CALCULATE STATE OF DOWNDRAUGHT
!-----------------------------------------------------------------------
!
         IF (BDD_START(I) .AND. BUOY1 >  0.5) THEN
            BDD_ON(I) = .FALSE.
         ELSE IF (BUOY1 >  0.5 .OR. K == 2) THEN
            B_DD_END(I) = .TRUE.
         END IF
       END DO
      END IF
!
      RETURN
      END SUBROUTINE TERMDD
!
#endif
