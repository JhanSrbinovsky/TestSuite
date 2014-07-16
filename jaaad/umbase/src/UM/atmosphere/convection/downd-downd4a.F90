#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DOWND--------------------------------------------------
!LL
!LL  PURPOSE : CALL DOWNDRAUGHT CALCULATION
!LL
!LL            CHANGE OF PHASE CALCULATION WHERE NO DOWNDRAUGHT OCCURS
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
      SUBROUTINE DOWND (NPNTS,NP_FULL,K,KCT,THDD_K,QDD_K,THE_K,THE_KM1, &
     &                  QE_K,QE_KM1,DTHBYDT_K,DTHBYDT_KM1,DQBYDT_K,     &
     &                  DQBYDT_KM1,FLX_DD_K,P_KM1,DELPK,DELPKM1,EXK,    &
     &                  EXKM1,DELTD,DELQD,AMDETK,EKM14,EKM34,PRECIP_K,  &
     &                  RAIN,SNOW,ICCB,BWATER_K,BDD_START,              &
     &                  BDDWT_K,BDDWT_KM1,BDD_ON,RAIN_ENV,SNOW_ENV,     &
     &                  RAIN_DD,SNOW_DD,FLX_UD_K,TIMESTEP,CCA,NDDON_A,  &
     &                  LR_UD_REF,                                      &
     &                  L_TRACER,NTRA,TRADD_K,TRAE_K,                   &
     &                  TRAE_KM1,DTRABYDT_K,DTRABYDT_KM1,DELTRAD)


      Use cv_run_mod, Only:                                             &
          dd_opt

      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! MODEL CONSTANTS
!-----------------------------------------------------------------------
!
#include "c_0_dg_c.h"
#include "c_g.h"
#include "ddptef.h"
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
!
      INTEGER I,KTRA             ! LOOP COUNTERS
!
      INTEGER K                  ! IN PRESENT MODEL LAYER
!
      INTEGER NPNTS              ! IN NUMBER OF POINTS
!
      INTEGER NDDON,NDDON_A      ! NUMBER OF POINTS AT WHICH
                                 ! DOWNDRAUGHT DOES OCCUR
!
      INTEGER NP_FULL            ! IN FULL VECTOR LENGTH
!
      INTEGER NTRA               ! NUMBER OF TRACERS
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!-----------------------------------------------------------------------
!
      INTEGER KCT                ! IN CONVECTIVE CLOUD TOP LAYER
!
      REAL THDD_K(NPNTS)         ! IN MODEL POTENTIAL TEMPERATURE
                                 !    OF DOWNDRAUGHT IN LAYER K (K)
!
      REAL QDD_K(NPNTS)          ! IN MIXING RATIO OF DOWNDRAUGHT IN
                                 !    LAYER K (KG/KG)
!

      REAL TRADD_K(NP_FULL,NTRA) ! IN TRACER CONTENT OF DOWNDRAUGHT
                                 !    IN LAYER K (KG/KG)
!
      REAL THE_K(NPNTS)          ! IN POTENTIAL TEMPERATURE OF
                                 !    ENVIRONMENT IN LAYER K (K)
!
      REAL THE_KM1(NPNTS)        ! IN POTENTIAL TEMPERATURE OF
                                 !    ENVIRONMENT IN LAYER K-1 (K)
!
      REAL QE_K(NPNTS)           ! IN MIXING RATIO OF ENVIRONMENT IN
                                 !    LAYER K (KG/KG)
!
      REAL QE_KM1(NPNTS)         ! IN MIXING RATIO OF ENVIRONMENT IN
                                 !    LAYER K-1 (KG/KG)
!
      REAL TRAE_K(NP_FULL,NTRA)  ! IN TRACER CONTENT OF ENVIRONMENT
                                 !    IN LAYER K (KG/KG)
!
      REAL TRAE_KM1(NP_FULL,NTRA)! IN TRACER CONTENT OF ENVIRONMENT
                                 !    IN LAYER K-1 (KG/KG)
!
      REAL FLX_DD_K(NPNTS)       ! IN DOWNDRAUGHT MASS FLUX OF LAYER K
                                 !    (PA/S)
!
      REAL P_KM1(NPNTS)          ! IN PRESSURE OF LAYER K-1 (PA)
!
      REAL DELPK(NPNTS)          ! IN PRESSURE DIFFERENCE ACROSS
                                 !    LAYER K (PA)
!
      REAL DELPKM1(NPNTS)        ! IN PRESSURE DIFFERENCE ACROSS
                                 !    LAYER K-1 (PA)
!
!
      REAL EXK(NPNTS)            ! IN EXNER RATIO FOR LAYER K
!
      REAL EXKM1(NPNTS)          ! IN EXNER RATIO FOR LAYER K-1
!
      REAL PRECIP_K(NPNTS)       ! IN PRECIPITATION ADDED WHEN
                                 !    DESCENDING FROM LAYER K TO K-1
                                 !    (KG/M**2/S)
!
      REAL AMDETK(NPNTS)         ! IN MIXING DETRAINMENT AT LEVEL K
                                 !    MULTIPLIED BY APPROPRIATE LAYER
                                 !    THICKNESS
!
      REAL EKM14(NPNTS)          ! IN EXNER RATIO AT LAYER K-1/4
!
      REAL EKM34(NPNTS)          ! IN EXNER RATIO AT LAYER K-3/4
!
      REAL DELTD(NPNTS)          ! IN COOLING NECESSARY TO
                                 !    ACHIEVE SATURATION (K)
!
      REAL DELQD(NPNTS)          ! IN MOISTENING NECESSARY TO
                                 !    ACHIEVE SATURATION (KG/KG)
!
      REAL DELTRAD(NP_FULL,NTRA) ! IN DEPLETION OF ENVIRONMENT TRACER
                                 !    DUE TO DOWNDRAUGHT FORMATION
!
      INTEGER ICCB(NPNTS)           ! IN CLOUD BASE LEVEL
!
      LOGICAL BWATER_K(NPNTS)    ! IN MASK FOR THOSE POINTS AT WHICH
                                 !    CONDENSATE IS WATER IN LAYER K
!
      LOGICAL BDDWT_K(NPNTS)     ! IN MASK FOR THOSE POINTS IN
                                 !    DOWNDRAUGHT WHERE PRECIPITATION
                                 !    IS LIQUID IN LAYER K
!
      LOGICAL BDDWT_KM1(NPNTS)   ! IN MASK FOR THOSE POINTS IN
                                 !    DOWNDRAUGHT WHERE PRECIPITATION
                                 !    IS LIQUID IN LAYER K-1
!
      LOGICAL L_TRACER           ! IN SWITCH FOR INCLUSION OF TRACERS
!
      REAL RAIN_ENV(NPNTS)       ! IN AMOUNT OF RAIN FALLING THROUGH
                                 !    THE ENVIRONMENT
!
      REAL SNOW_ENV(NPNTS)       ! IN AMOUNT OF SNOW FALLING THROUGH
                                 !    THE ENVIRONMENT
!
      REAL RAIN_DD(NPNTS)        ! IN AMOUNT OF RAIN FALLING THROUGH
                                 !    THE DOWNDRAUGHT
!
      REAL SNOW_DD(NPNTS)        ! IN AMOUNT OF SNOW FALLING THROUGH
                                 !    THE DOWNDRAUGHT
!
      REAL FLX_UD_K(NPNTS)       ! IN UPDRAUGHT MASSFLUX AT LAYER K
!
      REAL TIMESTEP              ! IN MODEL TIMESTEP (S)
!
      REAL CCA(NPNTS)            ! IN CONVECTIVE CLOUD AMOUNT
!
      REAL LR_UD_REF(NPNTS)      ! IN UD PPN MIXING RATION IN LOWEST
                                 !    PRECIPITATING LAYER IN UD
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!-----------------------------------------------------------------------
!
      LOGICAL BDD_START(NPNTS)   ! INOUT
                                 ! IN  MASK FOR THOSE POINTS WHERE
                                 !     DOWNDRAUGHT MAY FORM IN LAYER K
                                 ! OUT MASK FOR THOSE POINTS WHERE
                                 !     DOWNDRAUGHT MAY FORM IN LAYER
                                 !     K-1
!
      REAL DTHBYDT_K(NPNTS)      ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE OF LAYER K (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE OF LAYER K
                                 !     (K/S)
!
      REAL DTHBYDT_KM1(NPNTS)    ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE OF LAYER K-1 (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE OF
                                 !     LAYER K-1 (K/S)
!
      REAL DQBYDT_K(NPNTS)       ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING
                                 !     RATIO OF LAYER K (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     MIXING RATIO OF LAYER K (KG/KG/S)
!
      REAL DQBYDT_KM1(NPNTS)     ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING
                                 !     RATIO OF LAYER K-1 (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     MIXING RATIO OF
                                 !     LAYER K-1 (KG/KG/S)
!
      REAL DTRABYDT_K(NP_FULL,                                          &
                                 ! INOUT
     &                NTRA)      ! IN  INCREMENT TO MODEL TRACER OF
                                 !     LAYER K (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL

!
      REAL DTRABYDT_KM1(NP_FULL,                                        &
                                 ! INOUT
     &                 NTRA)     ! IN  INCREMENT TO MODEL TRACER OF
                                 !     LAYER K-1 (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     TRACER IN LAYER K-1
                                 !     (KG/KG/S)
!
      REAL RAIN (NPNTS)          ! INOUT
                                 ! IN  INITIALISED RAINFALL (KG/M**2/S)
                                 ! OUT SURFACE RAINFALL (KG/M**2/S)
!
      REAL SNOW(NPNTS)           ! INOUT
                                 ! IN  INITIALISED SNOWFALL (KG/M**2/S)
                                 ! OUT SURFACE SNOWFALL (KG/M**2/S)
!
      LOGICAL BDD_ON(NPNTS)      ! INOUT
                                 ! IN  MASK FOR THOSE POINTS WHERE DD
                                 !     HAS CONTINUED FROM PREVIOUS LAYER
                                 ! OUT MASK FOR THOSE POINTS WHERE DD
                                 !     CONTINUES TO LAYER K-1
!
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!-----------------------------------------------------------------------
!
      REAL WORK(NDDON_A,38)      !  WORK SPACE
!
      LOGICAL BWORK(NDDON_A,5)   !  WORK SPACE FOR 'BIT' MASKS
!
      INTEGER INDEX1(NDDON_A)    !  INDEX FOR COMPRESS AND
!
      LOGICAL B_DD_END(NPNTS)    !  MASK FOR POINTS WHERE DOWNDRAUGHT
                                 ! HAS ENDED
!
      REAL TRADD_K_C(NDDON_A,                                           &
                                 ! TRACER CONTENT IN DOWNDRAUGHT AT
     &               NTRA)       ! LAYER K - COMPRESSED (KG/KG)
!
      REAL TRAE_K_C(NDDON_A,NTRA)! TRACER CONTENT OF ENVIRONMENT AT
                                 ! LAYER K - COMPRESSED (KG/KG)
!
      REAL TRAE_KM1_C(NDDON_A,                                          &
                                 ! TRACER CONTENT OF ENVIRONMENT IN
     &                NTRA)      ! LAYER K-1 - COMPRESSED (KG/KG)
!
      REAL DTRA_K_C(NDDON_A,NTRA)! INCREMENT TO MODEL TRACER IN LAYER
                                 ! K - COMPRESSED (KG/KG/S)
!
      REAL DTRA_KM1_C(NDDON_A,                                          &
                                 ! INCREMENT TO MODEL TRACER IN LAYER
     &                NTRA)      ! K-1 - COMPRESSED (KG/KG/S)
!
      REAL DELTRAD_C(NDDON_A,                                           &
                                 ! DEPLETION OF ENVIRONMENT TRACER
     &               NTRA)       ! DUE TO DOWNDRAUGHT FORMATION
                                 ! COMPRESSED
!
      REAL PPN_MIX_DD_C(NDDON_A) ! PRECIP MIXING RATIO IN DD
!
      REAL FACTOR                !  PROPORTION OF RAINFALL GOING INTO
                                 !  DOWNDRAUGHT FROM UD
!
      REAL FACTOR_ENV            !  PROPORTION OF RAINFALL GOING INTO
                                 !  DD FROM FALLING PPN
!
      REAL PPN_DD_REF            !  REFERENCE DD PPN MASS
!
!-----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!-----------------------------------------------------------------------
!
      EXTERNAL CHG_PHSE, PEVP_BCB, DDRAUGHT
!
!-----------------------------------------------------------------------
! START OF MAIN LOOP
!   UPDATE PRECIPITATION AND CALCULATE MASK FOR WHERE PRECIPITATION
!   IS LIQUID
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
        B_DD_END(I) = .FALSE.
      END DO
!
      IF (K == KCT+1) THEN
        DO I=1,NPNTS
         RAIN_DD(I) = 0.0
         RAIN_ENV(I) = 0.0
        END DO
        DO I=1,NPNTS
         SNOW_DD(I) = 0.0
         SNOW_ENV(I) = 0.0
        END DO
      END IF
!
!----------------------------------------------------------------------
! INJECTION OF PRECIPITATION FROM UD AT LEVEL K
!----------------------------------------------------------------------
!
      IF (DD_OPT == 1) THEN
!----------------------------------------------------------------------
! CORRECTED DD
! TRANSFER ALL THE PRECIPITATION CREATED IN THE UD INTO THE ENVIRONMENT
!----------------------------------------------------------------------
       DO I=1,NPNTS

        IF (BWATER_K(I)) THEN
         RAIN_ENV(I) = RAIN_ENV(I) + PRECIP_K(I)
        ELSE
         SNOW_ENV(I) = SNOW_ENV(I) + PRECIP_K(I)
        END IF
        PRECIP_K(I) = 0.0
       END DO

      ELSE
!----------------------------------------------------------------------
! ORIGINAL DD
!----------------------------------------------------------------------
       DO I=1,NPNTS
        FACTOR = 0.0
        IF (BDD_ON(I) .AND. FLX_UD_K(I) >  0.0) THEN
         FACTOR = G * FLX_DD_K(I)/FLX_UD_K(I)
         FACTOR = AMIN1(FACTOR, 1.0)
        END IF
!
        IF (BWATER_K(I)) THEN
         RAIN_DD(I)  = RAIN_DD(I) + PRECIP_K(I)*FACTOR
         RAIN_ENV(I) = RAIN_ENV(I) + PRECIP_K(I)*(1.0-FACTOR)
        ELSE
         SNOW_DD(I)  = SNOW_DD(I) + PRECIP_K(I)*FACTOR
         SNOW_ENV(I) = SNOW_ENV(I) + PRECIP_K(I)*(1.0-FACTOR)
        END IF
!
       END DO
      ENDIF
!
!----------------------------------------------------------------------
! INTERACTION OF DOWNDRAUGHT WITH RESERVE OF PRECIPITATION OUTSIDE
! DOWNDRAUGHT
!
! BASED UPON CONTINUITY OF PRECIPITATION MIXING RATIO WITHIN
! DOWNDRAUGHT - EITHER AFTER INJECTION OF RAIN FROM UD IN LEVEL
! K OR WITH PPN MIXING RATIO IN LOWEST PRECIPITATING LAYER
!
! IF DOWNDRAUGHT INCREASES IN MASS THEN WATER INJECTED
! IF DOWNDRAUGHT DECREASES IN MASS THEN WATER IS REMOVED
!
!----------------------------------------------------------------------
!

      IF (DD_OPT == 1) THEN
!----------------------------------------------------------------------
! CORRECTED DD
!----------------------------------------------------------------------
       DO I=1,NPNTS
!
        IF (BDD_ON(I)) THEN
   
         IF (FLX_UD_K(I) > 0.0) THEN
          FACTOR_ENV = DDPTEF*FLX_DD_K(I)/FLX_UD_K(I)*  &
           DELPK(I)/5000.0
         ELSE
          FACTOR_ENV = 1.0*DELPK(I)/5000.0
         END IF
         FACTOR_ENV = AMIN1(FACTOR_ENV,1.0)


         IF (FACTOR_ENV >  0.0) THEN
          RAIN_DD(I)  = RAIN_DD(I) + RAIN_ENV(I)*FACTOR_ENV
          RAIN_ENV(I) = RAIN_ENV(I) * (1.0-FACTOR_ENV)
          SNOW_DD(I)  = SNOW_DD(I) + SNOW_ENV(I)*FACTOR_ENV
          SNOW_ENV(I) = SNOW_ENV(I) * (1.0-FACTOR_ENV)
         ELSE
          RAIN_ENV(I) = RAIN_ENV(I) - RAIN_DD(I)*FACTOR_ENV
          RAIN_DD(I)  = RAIN_DD(I) * (1.0+FACTOR_ENV)
          SNOW_ENV(I) = SNOW_ENV(I) - SNOW_DD(I)*FACTOR_ENV
          SNOW_DD(I)  = SNOW_DD(I) * (1.0+FACTOR_ENV)
         END IF
!
        END IF
!
       END DO

      ELSE
!----------------------------------------------------------------------
! ORIGINAL DD
!----------------------------------------------------------------------
       DO I=1,NPNTS
!
        IF (BDD_ON(I)) THEN
!
         FACTOR_ENV = 0.0
         IF (PRECIP_K(I) >  0.0 .AND. FLX_UD_K(I) >  0.0) THEN
!
!---------------------------------------------------------------------
! CALCULATE NEW REFERENCE PPN MIXING RATIO
! DD PPN MIXING RATIO IN LAYER KM1 BASED ON CONTINUITY
! WITH THAT IN LAYER K
!---------------------------------------------------------------------
!
          LR_UD_REF(I) = G * PRECIP_K(I)/FLX_UD_K(I)
          PPN_DD_REF   = RAIN_DD(I)+SNOW_DD(I)
         ELSE
!
!---------------------------------------------------------------------
! DD PPN MIXING RATIO IN LAYER KM1 BASED ON CONTINUITY
! WITH THAT IN LAST PRECIPITATING UD LAYER
!---------------------------------------------------------------------
!
          PPN_DD_REF = LR_UD_REF(I) * FLX_DD_K(I)
         END IF
!
!--------------------------------------------------------------------
! INJECT PPN INTO DD FROM PPN FALLING OUTSIDE OF THE DD
!--------------------------------------------------------------------
!
         IF ((RAIN_ENV(I) + SNOW_ENV(I))  >   0.0) THEN
!-------Already inside IF ( BDD_ON(I)) block----------------------------
          FACTOR_ENV = ( (PPN_DD_REF * (1.0+EKM14(I))*                   &
     &                    (1.0+EKM34(I))*(1.0-AMDETK(I))) -             &
     &                         (RAIN_DD(I)+SNOW_DD(I)) ) /              &
     &                          (RAIN_ENV(I)+SNOW_ENV(I))
          FACTOR_ENV = AMIN1(FACTOR_ENV,1.0)
          FACTOR_ENV = AMAX1(FACTOR_ENV,-1.0)
         END IF
!
         IF (FACTOR_ENV >  0.0) THEN
          RAIN_DD(I)  = RAIN_DD(I) + RAIN_ENV(I)*FACTOR_ENV
          RAIN_ENV(I) = RAIN_ENV(I) * (1.0-FACTOR_ENV)
          SNOW_DD(I)  = SNOW_DD(I) + SNOW_ENV(I)*FACTOR_ENV
          SNOW_ENV(I) = SNOW_ENV(I) * (1.0-FACTOR_ENV)
         ELSE
          RAIN_ENV(I) = RAIN_ENV(I) - RAIN_DD(I)*FACTOR_ENV
          RAIN_DD(I)  = RAIN_DD(I) * (1.0+FACTOR_ENV)
          SNOW_ENV(I) = SNOW_ENV(I) - SNOW_DD(I)*FACTOR_ENV
          SNOW_DD(I)  = SNOW_DD(I) * (1.0+FACTOR_ENV)
         END IF
!
        END IF
!
!--------------------------------------------------------------------
! ZERO PRECIPITATION RATE IN LAYER K
!--------------------------------------------------------------------
!
        PRECIP_K(I) = 0.0
!
       END DO

      ENDIF   !(DD_OPT)
!
!
!-----------------------------------------------------------------------
! COMPRESS OUT ON BASIS OF BIT VECTOR BDDON - THOSE POINTS WITH A
! DOWNDRAUGHT
!-----------------------------------------------------------------------
!
      NDDON=0
!
      DO I=1,NPNTS
        IF (BDD_ON(I)) THEN
           NDDON = NDDON+1
           INDEX1(NDDON) = I
        END IF
      END DO
!
      IF (NDDON  /=  0) THEN
         DO I=1,NDDON
          WORK(I,1) = THDD_K(INDEX1(I))
          WORK(I,2) = QDD_K(INDEX1(I))
         END DO
         DO I=1,NDDON
          WORK(I,3) = THE_K(INDEX1(I))
          WORK(I,4) = THE_KM1(INDEX1(I))
         END DO
         DO I=1,NDDON
          WORK(I,5) = QE_K(INDEX1(I))
          WORK(I,6) = QE_KM1(INDEX1(I))
         END DO
         DO I=1,NDDON
          WORK(I,7) = DTHBYDT_K(INDEX1(I))
          WORK(I,8) = DTHBYDT_KM1(INDEX1(I))
         END DO
         DO I=1,NDDON
          WORK(I,9) = DQBYDT_K(INDEX1(I))
          WORK(I,10) = DQBYDT_KM1(INDEX1(I))
         END DO
         DO I=1,NDDON
          WORK(I,11) = FLX_DD_K(INDEX1(I))
          WORK(I,12) = P_KM1(INDEX1(I))
         END DO
         DO I=1,NDDON
          WORK(I,13) = DELPK(INDEX1(I))
          WORK(I,14) = DELPKM1(INDEX1(I))
         END DO
         DO I=1,NDDON
          WORK(I,15) = EXK(INDEX1(I))
          WORK(I,16) = EXKM1(INDEX1(I))
         END DO
         DO I=1,NDDON
          WORK(I,17) = DELTD(INDEX1(I))
          WORK(I,18) = DELQD(INDEX1(I))
         END DO
         DO I=1,NDDON
          WORK(I,19) = AMDETK(INDEX1(I))
          WORK(I,20) = EKM14(INDEX1(I))
         END DO
         DO I=1,NDDON
          WORK(I,21) = EKM34(INDEX1(I))
          WORK(I,22) = RAIN_DD(INDEX1(I))
         END DO
         DO I=1,NDDON
          WORK(I,23) = SNOW_DD(INDEX1(I))
          WORK(I,24) = CCA(INDEX1(I))
         END DO
         DO I=1,NDDON
          BWORK(I,1) = BDD_START(INDEX1(I))
          BWORK(I,2) = BDDWT_K(INDEX1(I))
         END DO
         DO I=1,NDDON
          BWORK(I,3) = BDDWT_KM1(INDEX1(I))
          BWORK(I,4) = BDD_ON(INDEX1(I))
          BWORK(I,5) = B_DD_END(INDEX1(I))
      END DO      
!
         DO I=1,NDDON
          PPN_MIX_DD_C(I) = G*(RAIN_DD(INDEX1(I)) + &
           SNOW_DD(INDEX1(I)))/FLX_DD_K(INDEX1(I))
         END DO

      IF(L_TRACER)THEN
!
      DO KTRA=1,NTRA
        DO I=1,NDDON
          TRADD_K_C(I,KTRA) = TRADD_K(INDEX1(I),KTRA)
          TRAE_K_C(I,KTRA) = TRAE_K(INDEX1(I),KTRA)
         END DO
         DO I=1,NDDON
          TRAE_KM1_C(I,KTRA) = TRAE_KM1(INDEX1(I),KTRA)
          DTRA_K_C(I,KTRA)  = DTRABYDT_K(INDEX1(I),KTRA)
         END DO
         DO I=1,NDDON
          DTRA_KM1_C(I,KTRA) = DTRABYDT_KM1(INDEX1(I),KTRA)
          DELTRAD_C(I,KTRA) = DELTRAD(INDEX1(I),KTRA)
        END DO
      END DO
!
      END IF
!
!
!-----------------------------------------------------------------------
! START DOWNDRAUGHT CALCULATION
!-----------------------------------------------------------------------
!
!
! DEPENDS ON: ddraught
         CALL DDRAUGHT (NDDON,NDDON_A,K,KCT,WORK(1,1),WORK(1,2),        &
     &                  WORK(1,3),WORK(1,4),WORK(1,5),WORK(1,6),        &
     &                  WORK(1,7),WORK(1,8),WORK(1,9),WORK(1,10),       &
     &                  WORK(1,11),WORK(1,12),WORK(1,13),WORK(1,14),    &
     &                  WORK(1,15),WORK(1,16),WORK(1,17),WORK(1,18),    &
     &                  WORK(1,19),WORK(1,20),WORK(1,21),WORK(1,22),    &
     &                  WORK(1,23),BWORK(1,1),BWORK(1,2),BWORK(1,3),    &
     &                  BWORK(1,4),BWORK(1,5),WORK(1,24),               &
                        PPN_MIX_DD_C,                                   &  
     &                  L_TRACER,NTRA,TRADD_K_C,TRAE_K_C,TRAE_KM1_C,    &
     &                  DTRA_K_C,DTRA_KM1_C,DELTRAD_C)
!
!-----------------------------------------------------------------------
! EXPAND REQUIRED VECTORS BACK TO FULL FIELDS
!-----------------------------------------------------------------------
!
      DO I=1,NDDON
       THDD_K(INDEX1(I)) = WORK(I,1)
       QDD_K(INDEX1(I)) = WORK(I,2)
      END DO
      DO I=1,NDDON
       DTHBYDT_K(INDEX1(I)) = WORK(I,7)
       DTHBYDT_KM1(INDEX1(I)) = WORK(I,8)
      END DO
      DO I=1,NDDON
       DQBYDT_K(INDEX1(I)) = WORK(I,9)
       DQBYDT_KM1(INDEX1(I)) = WORK(I,10)
      END DO
      DO I=1,NDDON
       FLX_DD_K(INDEX1(I)) = WORK(I,11)
       RAIN_DD(INDEX1(I)) = WORK(I,22)
      END DO
      DO I=1,NDDON
       SNOW_DD(INDEX1(I)) = WORK(I,23)
       BDD_START(INDEX1(I)) = BWORK(I,1)
      END DO
      DO I=1,NDDON
       BDDWT_K(INDEX1(I)) = BWORK(I,2)
       BDDWT_KM1(INDEX1(I)) = BWORK(I,3)
      END DO
      DO I=1,NDDON
       BDD_ON(INDEX1(I)) = BWORK(I,4)
       B_DD_END(INDEX1(I)) = BWORK(I,5)
      END DO
!
      IF(L_TRACER)THEN
!
      DO KTRA=1,NTRA
        DO I=1,NDDON
          TRADD_K(INDEX1(I),KTRA) = TRADD_K_C(I,KTRA)
          DTRABYDT_K(INDEX1(I),KTRA) = DTRA_K_C(I,KTRA)
          DTRABYDT_KM1(INDEX1(I),KTRA) = DTRA_KM1_C(I,KTRA)
        END DO
      END DO
!
      END IF
!
      END IF
!
!-----------------------------------------------------------------------
! RESET PRECIPITATION FALLING THROUGH ENVIRONMENT IF DOWNDRAUGHT
! DID NOT FORM
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
        IF (.NOT.BDD_ON(I).AND..NOT.B_DD_END(I)) THEN
          RAIN_ENV(I) = RAIN_ENV(I)+RAIN_DD(I)
          SNOW_ENV(I) = SNOW_ENV(I)+SNOW_DD(I)
          RAIN_DD(I) = 0.0
          SNOW_DD(I) = 0.0
        END IF
      END DO
!
!-----------------------------------------------------------------------
! CARRY OUT CHANGE OF PHASE CALCULATION FOR PRECIPITATION FALLING
! THROUGH ENVIRONMENT
!-----------------------------------------------------------------------
!
! DEPENDS ON: chg_phse
         CALL CHG_PHSE (NPNTS,K,RAIN_ENV,SNOW_ENV,DTHBYDT_KM1,          &
                        EXK,EXKM1,DELPKM1,THE_K,THE_KM1,TIMESTEP,CCA)
!
!-----------------------------------------------------------------------
! EVAPORATE RAIN FALLING THROUGH ENVIRONMENT IF LAYER K BELOW
! CLOUD BASE
!-----------------------------------------------------------------------
!
! DEPENDS ON: pevp_bcb
         CALL PEVP_BCB (NPNTS,K-1,ICCB,THE_KM1,P_KM1,QE_KM1,DELPKM1,    &
     &                  RAIN_ENV,SNOW_ENV,DTHBYDT_KM1,DQBYDT_KM1,       &
     &                  EXKM1,TIMESTEP,CCA)
!
!-----------------------------------------------------------------------
! RESET PRECIPITATION FALLING THROUGH ENVIRONMENT IF DOWNDRAUGHT
! TERMINATES
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
        IF (B_DD_END(I)) THEN
          RAIN_ENV(I) = RAIN_ENV(I)+RAIN_DD(I)
          SNOW_ENV(I) = SNOW_ENV(I)+SNOW_DD(I)
          RAIN_DD(I) = 0.0
          SNOW_DD(I) = 0.0
        END IF
      END DO
!
!-----------------------------------------------------------------------
! UPDATE RAIN AND SNOW
!-----------------------------------------------------------------------
!
       IF (K == 2) THEN
         DO I=1,NPNTS
           RAIN(I) = RAIN(I)+RAIN_DD(I)+RAIN_ENV(I)
           SNOW(I) = SNOW(I)+SNOW_DD(I)+SNOW_ENV(I)
         END DO
       END IF
!
      RETURN
      END SUBROUTINE DOWND
!
#endif
