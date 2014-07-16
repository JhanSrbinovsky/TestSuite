#if defined(A05_4A) || defined (A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE LAYER_DD--------------------------------------------
!LL
!LL  PURPOSE : CALCULATES LAYER DEPENDENT CONSTANTS FOR LAYER K
!LL            -PRESSURE
!LL            -LAYER THICKNESS
!LL            -ENTRAINMENT COEFFICIENTS
!LL            -DETRAINMENT COEFFICIENTS
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL  VERSION NO. 4  DATED 5/2/92
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE LAYER_DD(NPNTS,K,KCT,THE_K,THE_KM1,FLX_STRT,           &
     &                    P_LAYER_CENTRES,                              &
     &                    P_LAYER_BOUNDARIES,                           &
     &                    EXNER_LAYER_CENTRES,                          &
     &                    EXNER_KM12,                                   &
     &                    EXNER_KP12,EXNER_KM32,PSTAR,PK,PKM1,DELPK,    &
     &                    DELPKM1,EXK,EXKM1,AMDETK,EKM14,EKM34,KMIN,    &
     &                    BDDI, recip_pstar )

      Use cv_run_mod, Only:                                             &
          dd_opt

      IMPLICIT NONE
!
!----------------------------------------------------------------------
! MODEL CONSTANTS
!----------------------------------------------------------------------
#include "c_0_dg_c.h"
#include "c_r_cp.h"
#include "entcnst.h"
#include "entdd.h"
#include "ddkmdet.h"
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTER
!----------------------------------------------------------------------
!
      INTEGER NPNTS             ! IN VECTOR LENGTH
!
      INTEGER K                 ! IN PRESENT MODEL LAYER
!
      INTEGER I                 ! COUNTER FOR DO LOOPS
!
      INTEGER KCT               ! IN CONVECTIVE CLOUD TOP LAYER
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!

      REAL P_LAYER_CENTRES(NPNTS,0:KCT+2)
      REAL P_LAYER_BOUNDARIES(NPNTS,0:KCT+1)
      REAL EXNER_LAYER_CENTRES(NPNTS,0:KCT+1)
!
!
      REAL PSTAR(NPNTS)         ! IN SURFACE PRESSURE (PA)
!
      REAL EXNER_KM12(NPNTS)    ! IN EXNER FUNCTION AT LAYER K-1/2
!
      REAL EXNER_KP12(NPNTS)    ! IN EXNER FUNCTION AT LAYER K+1/2
!
      REAL EXNER_KM32(NPNTS)    ! IN EXNER FUNCTION AT LAYER K-3/2
!
      REAL FLX_STRT(NPNTS)      ! IN UPDRAUGHT MASSFLUX AT LEVEL WHERE
                                !    DOWNDRAUGHT STARTS (PA/S)
!
      REAL THE_K(NPNTS)         ! IN POTENTIAL TEMPERATURE OF
                                !    ENVIRONMENT IN LAYER K (K)
!
      REAL THE_KM1(NPNTS)       ! IN POTENTIAL TEMPERATURE OF
                                !    ENVIRONMENT IN LAYER K-1 (K)
!
      LOGICAL BDDI(NPNTS)       ! IN MASK FOR POINTS WHERE DOWNDRAUGHT
                                !    MAY INITIATE
      REAL recip_PSTAR(NPNTS)! Reciprocal of pstar array

!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!----------------------------------------------------------------------
!
      INTEGER KMIN(NPNTS)       ! INOUT
                                ! FREEZING LEVEL
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!----------------------------------------------------------------------
!
      REAL PK(NPNTS)            ! OUT PRESSURE AT LAYER K (PA)
!
      REAL PKM1(NPNTS)          ! OUT PRESSURE AT LAYER K-1 (PA)
!
      REAL DELPK(NPNTS)         ! OUT THICKNESS OF LAYER K (PA)
!
      REAL DELPKM1(NPNTS)       ! OUT THICHNESS OF LAYER K-1 (PA)
!
      REAL EKM14(NPNTS)         ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K-1/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
!
      REAL EKM34(NPNTS)         ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K-3/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
!
      REAL AMDETK(NPNTS)        ! OUT MIXING DETRAINMENT COEFFICIENT
                                !     AT LEVEL K MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
!
      REAL EXK(NPNTS)           ! OUT EXNER FUNCTION AT LEVEL K
!
      REAL EXKM1(NPNTS)         ! OUT EXNER FUNCTION AT LEVEL K-1
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
!
      REAL TTK                  ! TEMPERATURE STORE AT LAYER K
!
      REAL TTKM1                ! TEMPERATURE STORE AT LAYER K-1
!
      REAL THKM12               ! POTENTIAL TEMPERATURE STORE AT
                                ! LAYER K-1/2
!
      REAL TTKM12               ! TEMPERATURE STORE AT LAYER K-1/2
!
      REAL INCR_FAC             ! INCREMENT FACTOR FOR ENTRAINMENT
                                ! RATES AT FREEZING LEVEL
!
      REAL                                                              &
     &    PU,PL
      REAL DDCOEF2A             ! COEFFICIENT USED IN CALCULATION 
                                ! OF DOWNDRAUGHT ENTRAINMENT RATES
!                               
!----------------------------------------------------------------------
! SET KMIN TO INITIAL VALUE
!L CALCULATE PK, DELPK AND EXNER FUNCTION - IF K = KCT THEN
!L VALUES FOR PREVIOUS PASS THROUGH ROUTINE AT (K-1)+1 ARE TAKEN
!----------------------------------------------------------------------
!
      IF (K == KCT+1) THEN
       DO I=1,NPNTS
        KMIN(I) = KCT+2
        PK(I) = P_LAYER_CENTRES(I,K)
        DELPK(I) =  -(P_LAYER_BOUNDARIES(I,K) -                         &
     &              P_LAYER_BOUNDARIES(I,K-1))
        EXK(I) = EXNER_LAYER_CENTRES(I,K)
       END DO
      ELSE
       DO I=1,NPNTS
        PK(I) = PKM1(I)
        DELPK(I) = DELPKM1(I)
        EXK(I) = EXKM1(I)
       END DO
      END IF
!L
!L---------------------------------------------------------------------
!L CALCULATE PKM1, DELPKM1
!L CALCULATE EXNER FUNCTIONS AT MID-LAYES K AND K-1, AND
!L DIFFERENCE OF EXNER FUNCTION ACROSS LAYER K
!L---------------------------------------------------------------------
!L
      DO I=1,NPNTS
        PKM1(I) = P_LAYER_CENTRES(I,K-1)
        DELPKM1(I) =  -(P_LAYER_BOUNDARIES(I,K-1)                       &
     &                - P_LAYER_BOUNDARIES(I,K-2))
        EXKM1(I) = EXNER_LAYER_CENTRES(I,K-1)
      End Do
!
!---------------------------------------------------------------------
! SET DDCOEF2A DEPENDING UPON WHICH REVISION OF THE DD CODE IS USED.
!---------------------------------------------------------------------
!
      IF (DD_OPT == 1) THEN 
        DDCOEF2A=2.0
      ELSE
        DDCOEF2A=3.0
      END IF
!L
!L---------------------------------------------------------------------
!L CALCULATE FREEZING LEVEL : CHECK IF FREEZING LEVEL IN THIS LAYER
!L---------------------------------------------------------------------
!L
      DO I=1,NPNTS
       IF (KMIN(I) == KCT+2) THEN
        TTK = THE_K(I)*EXK(I)
        TTKM1 = THE_KM1(I)*EXKM1(I)
        THKM12 = (THE_KM1(I)+THE_K(I))*0.5
        TTKM12 = THKM12*EXNER_KM12(I)
        IF (TTKM12  >=  TM .AND. TTK  <   TM) THEN
           KMIN(I) = K
        ELSE IF (TTKM1  >=  TM .AND. TTKM12  <   TM) THEN
           KMIN(I) = K-1
        END IF
       END IF
!
!L
!L---------------------------------------------------------------------
!L CALCULATE ENTRAINMENT COEFFICIENTS MULTIPLIED BY
!L APPROPRIATE LAYER THICKNESS
!L
!L CALCULATE MIXING DETRAINMENT COEFFICIENT MULTIPLIED BY
!L APPROPRIATE LAYER THICKNESS
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (2C), EQUATION(14)
!L---------------------------------------------------------------------
!L
      IF (PK(I) <  PSTAR(I)-DET_LYR) THEN
       EKM14(I) = AE2 *                                                 &
     &         (P_LAYER_BOUNDARIES(I,K-1)-PK(I)) * recip_PSTAR(I)
       EKM34(I) = AE2 *                                                 &
     &         (PKM1(I)-P_LAYER_BOUNDARIES(I,K-1)) * recip_PSTAR(I)
       AMDETK(I) = (EKM14(I)+EKM34(I)) * (1.0-1.0/AE2)
      ELSE
       EKM14(I) = 0.0
       EKM34(I) = 0.0
       AMDETK(I) = DELPK(I) /(PSTAR(I)-                                 &
     &             P_LAYER_BOUNDARIES(I,K))
      END IF
!
      IF (BDDI(I)) THEN
!
      
      IF (K == KMIN(I) .AND. PK(I) <  PSTAR(I)-DET_LYR) THEN
        INCR_FAC = FLX_STRT(I)*DDCOEF1*recip_pstar(I)
        IF (INCR_FAC >  6.0) INCR_FAC=6.0
        EKM14(I) = EKM14(I)*INCR_FAC
        EKM34(I) = EKM34(I)*INCR_FAC
      ELSE
        EKM14(I) = EKM14(I)*DDCOEF2A
        EKM34(I) = EKM34(I)*DDCOEF2A
        IF ((DD_OPT == 1) .OR. (KMIN(I) /= KCT+2 .AND. K <  KMIN(I) .AND. PK(I) <           &
     &    PSTAR(I)-DET_LYR)) THEN
          AMDETK(I) = AMDETK(I)*DDCOEF2A
        END IF
      END IF
!
      END IF
      END DO
!
      RETURN
      END SUBROUTINE LAYER_DD
!
#endif
