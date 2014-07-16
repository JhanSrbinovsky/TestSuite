#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DDRAUGHT-----------------------------------------------
!LL
!LL  PURPOSE : DOWNDRAUGHT ROUTINE
!LL
!LL            CONVECTIVE DOWNDRAUGHT BASED ON PARCEL THEORY
!LL
!LL            CARRY OUT DRY DESCENT
!LL
!LL            CALCULATE SUBSATURATION
!LL
!LL            CALCULATE EFFECT ON THE ENVIRONMENT
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
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
      SUBROUTINE DDRAUGHT (NPNTS,NP_FULL,K,KCT,THDD_K,QDD_K,THE_K,      &
     &                     THE_KM1,QE_K,QE_KM1,DTHBYDT_K,DTHBYDT_KM1,   &
     &                     DQBYDT_K,DQBYDT_KM1,FLX_DD_K,P_KM1,DELPK,    &
     &                     DELPKM1,EXK,EXKM1,DELTD,DELQD,AMDETK,EKM14,  &
     &                     EKM34,RAIN,SNOW,BDD_START,BDDWT_K,BDDWT_KM1, &
     &                     BDD_ON,B_DD_END,CCA,PPN_MIX_DD,L_TRACER,     &
     &                     NTRA,TRADD_K,TRAE_K,TRAE_KM1,DTRABYDT_K,     &
     &                     DTRABYDT_KM1,DELTRAD)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! MODEL CONSTANTS
!-----------------------------------------------------------------------
!
#include "c_0_dg_c.h"
#include "ddkmdet.h"
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
!
      INTEGER I,KTRA                ! LOOP COUNTERS
!
      INTEGER NPNTS                 ! IN NUMBER OF POINTS
!
      INTEGER NP_FULL               ! IN FULL VECTOR LENGTH
!
      INTEGER NTRA                  ! NUMBER OF TRACERS
!
      INTEGER K                     ! IN PRESENT MODEL LAYER
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!-----------------------------------------------------------------------
!
      INTEGER KCT                   ! IN CONVECTIVE CLOUD TOP
!
      REAL THE_KM1(NPNTS)           ! IN POTENTIAL TEMPERATURE OF
                                    !    ENVIRONMENT IN LAYER K-1 (K)
!
      REAL QE_KM1(NPNTS)            ! IN MIXING RATIO OF ENVIRONMENT IN
                                    !    LAYER K-1 (KG/KG)
!
      REAL TRAE_KM1(NP_FULL,NTRA)   ! IN TRACER CONTENT OF ENVIRONMENT
                                    !    IN LAYER K-1 (KG/KG)
!
      REAL P_KM1(NPNTS)             ! IN PRESSURE OF LAYER K-1 (PA)
!
      REAL DELPK(NPNTS)             ! IN CHANGE IN PRESSURE ACROSS
                                    !    LAYER K (PA)
!
      REAL DELPKM1(NPNTS)           ! IN CHANGE IN PRESSURE ACROSS
                                    !    LAYER K-1 (PA)
!
      REAL EXK(NPNTS)               ! IN EXNER RATIO IN LAYER K
!
      REAL EXKM1(NPNTS)             ! IN EXNER RATIO IN LAYER K-1
!
      REAL AMDETK(NPNTS)            ! IN MIXING DETRAINMENT RATE
!
      REAL EKM14(NPNTS)             ! IN EXNER RATIO AT LAYER K-1/4
!
      REAL EKM34(NPNTS)             ! IN EXNER RATIO AT LAYER K-3/4
!
      REAL DELTD(NPNTS)             ! IN COOLING NECESSARY TO ACHIEVE
                                    !    SATURATION (K)
!
      REAL DELQD(NPNTS)             ! IN MOISTENING NECESSARY TO ACHIEVE
                                    !    SATURATION (KG/KG)
!
      REAL DELTRAD(NP_FULL,NTRA)    ! IN DEPLETION OF ENV. TRACER DUE
                                    !    TO DOWNDRAUGHT FORMATION
!
      LOGICAL BDDWT_K(NPNTS)        ! IN MASK FOR THOSE POINTS IN
                                    !    DOWNDRAUGHT WHERE PRECIPITATION
                                    !    IS LIQUID IN LAYER K
!
      LOGICAL BDDWT_KM1(NPNTS)      ! IN MASK FOR THOSE POINTS IN
                                    !    DOWNDRAUGHT WHERE PRECIPITATION
                                    !    IS LIQUID IN LAYER K-1
!
      LOGICAL L_TRACER              ! IN SWITCH FOR INCLUSION OF TRACERS
!
      REAL CCA(NPNTS)               ! IN CONVECTIVE CLOUD AMOUNT
!
      REAL PPN_MIX_DD(NPNTS)        ! IN PRECIP MIXING RATIO      
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!-----------------------------------------------------------------------
!
      REAL THDD_K(NPNTS)            ! INOUT
                                    ! IN  POTENTIAL TEMPERATURE OF
                                    !     DOWNDRAUGHT IN LAYER K (K)
                                    ! OUT POTENTIAL TEMPERATURE RESET
                                    !     FOR NEXT LAYER (K)
!
      REAL QDD_K(NPNTS)             ! INOUT
                                    ! IN  DOWNDRAUGHT MIXING RATIO OF
                                    !     LAYER K (KG/KG)
                                    ! OUT MIXING RATIO RESET FOR NEXT
                                    !     LAYER (KG/KG)
!
!
      REAL TRADD_K(NP_FULL,NTRA)    ! INOUT
                                    ! IN  DOWNDRAUGHT TRACER CONTENT OF
                                    !     LAYER K (KG/KG)
                                    ! OUT TRACER CONTENT RESET FOR NEXT
                                    !     LAYER (KG/KG)
!
      REAL THE_K(NPNTS)             ! INOUT
                                    ! IN  POTENTIAL TEMPERATURE OF
                                    !     ENVIRONMENT IN LAYER K (K)
                                    ! OUT ENVIRONMENT POTENTIAL
                                    !     TEMPERATURE RESET FOR NEXT
                                    !     LAYER (K)
!
      REAL QE_K(NPNTS)              ! INOUT
                                    ! IN  MIXING RATIO OF ENVIRONMENT
                                    !     LAYER K (KG/KG)
                                    ! OUT ENVIRONMENT MIXING RATIO
                                    !     RESET FOR NEXT LAYER (KG/KG)
!
      REAL TRAE_K(NP_FULL,NTRA)     ! INOUT
                                    ! IN  TRACER CONTENT OF ENVIRONMENT
                                    !     IN LAYER K (KG/KG)
                                    ! OUT ENVIRONMENT TRACER CONTENT
                                    !     RESET FOR NEXT LAYER (KG/KG)
!
      REAL FLX_DD_K(NPNTS)          ! INOUT
                                    ! IN  DOWNDRAUGHT MASS FLUX OF
                                    !     LAYER K (PA/S)
                                    ! OUT DOWNDRAUGHT MASS FLUX RESET
                                    !      FOR NEXT LAYER (PA/S)
!
      REAL RAIN(NPNTS)              ! INOUT
                                    ! IN  AMOUNT OF RAIN (KG/M**2/S)
                                    ! OUT UPDATED RAINFALL (KG/M**2/S)
!
      REAL SNOW(NPNTS)              ! INOUT
                                    ! IN  AMOUNT OF SNOW(KG/M**2/S)
                                    ! OUT UPDATED SNOWFALL (KG/M**2/S)
!
      REAL DTHBYDT_K(NPNTS)         ! INOUT
                                    ! IN  INCREMENT TO MODEL POTENTIAL
                                    !     TEMPERATURE OF LAYER K (K/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     POTENTIAL TEMPERATURE IN
                                    !     LAYER K (K/S)
!
      REAL DTHBYDT_KM1(NPNTS)       ! INOUT
                                    ! IN  INCREMENT TO MODEL POTENTIAL
                                    !     TEMPERATURE IN LAYER K-1 (K/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     POTENTIAL TEMPERATURE IN
                                    !     LAYER K-1 (K/S)
!
      REAL DQBYDT_K(NPNTS)          ! INOUT
                                    ! IN  INCREMENT TO MODEL MIXING
                                    !     RATIO IN LAYER K (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     MIXING RATIO IN LAYER K
                                    !     (KG/KG/S)
!
      REAL DQBYDT_KM1(NPNTS)        ! INOUT
                                    ! IN  INCREMENT TO MODEL MIXING
                                    !     RATIO IN LAYER K-1 (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     MIXING RATIO IN LAYER K-1
                                    !     (KG/KG/S)
!
      REAL DUBYDT_K(NPNTS)          ! INOUT
                                    ! IN  INCREMENT TO MODEL U IN
                                    !     UV LAYER K (M/S**2)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     U IN UV LAYER K (M/S**2)
!
      REAL DUBYDT_KM1(NPNTS)        ! INOUT
                                    ! IN  INCREMENT TO MODEL U IN
                                    !     UV LAYER K-1 (M/S**2)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     U IN UV LAYER K-1 (M/S**2)
!
      REAL DTRABYDT_K(NP_FULL,NTRA) ! INOUT
                                    ! IN  INCREMENT TO MODEL TRACER
                                    !     CONTENTOF LAYER K (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     TRACER CONTENT IN LAYER K
                                    !     (KG/KG/S)
!
      REAL DTRABYDT_KM1(NP_FULL,                                        &
                                    ! INOUT
     &                  NTRA)       ! IN  INCREMENT TO MODEL TRACER
                                    !     CONTENT IN LAYER K-1
                                    !     (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     TRACER CONTENT IN LAYER K-1
                                    !     (KG/KG/S)
!
      LOGICAL BDD_ON(NPNTS)         ! INOUT
                                    ! IN  MASK FOR THOSE POINTS WHERE DD
                                    !     HAS CONTINUED FROM LAYER K+1
                                    ! OUT MASK FOR THOSE POINTS WHERE DD
                                    !     CONTINUES TO LAYER K-1
!
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!-----------------------------------------------------------------------
!
      LOGICAL BDD_START(NPNTS)      ! OUT MASK FOR THOSE POINTS WHERE
                                    !     DOWNDRAUGHT MAY START IN
                                    !     LAYER K-1
!
      LOGICAL B_DD_END(NPNTS)       ! OUT MASK FOR THOSE POINTS WHERE
                                    !     DOWNDRAUGHT IS ENDING IN
                                    !     LAYER K-1
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!-----------------------------------------------------------------------
!
!
      REAL THDD_KM1(NPNTS)          ! POTENTIAL TEMPERATURE OF
                                    ! DOWNDRAUGHT IN LAYER K-1 (K)
!
      REAL QDD_KM1(NPNTS)           ! DOWNDRAUGHT MIXING RATIO OF
                                    ! LAYER K-1 (KG/KG)
!
      REAL TRADD_KM1(NPNTS,NTRA)    ! TRACER CONTENT OF DOWNDRAUGHT
                                    ! IN LAYER K-1 (KG/KG)
!
      REAL QSATDD(NPNTS)            ! SATURATED DOWNDRAUGHT MIXING
                                    ! RATIO (KG/KG)
!
      REAL TDD_KM1(NPNTS)           ! TEMPERATURE OF DOWNDRAUGHT
                                    ! IN LAYER K-1 (K)
!
      REAL THDDS(NPNTS)             ! POTENTIAL TEMPERATURE OF
                                    ! SATURATED DOWNDRAUGHT (K)
!
      REAL QDDS(NPNTS)              ! SATURATED DOWNDRAUGHT MIXING
                                    ! RATIO (KG/KG)
!
      REAL FLX_DD_KM1(NPNTS)        ! DOWNDRAUGHT MASS FLUX IN
                                    ! LAYER K-1 (PA/S)
!
      REAL RAIN_TMP(NPNTS)          ! LIQUID PRECIPITATION STORE
!
      REAL SNOW_TMP(NPNTS)          ! SNOW STORE
!
!
!-----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!-----------------------------------------------------------------------
!
      EXTERNAL SATCAL, CRS_FRZL, QSAT, DEVAP, TERMDD,                   &
     &         DD_ENV, EVP
!
!-----------------------------------------------------------------------
! CALCULATE MASK FOR THOSE POINTS IN DOWNDRAUGHT WHERE PRECIPITATION
! IS LIQUID
!
! STORE PRECIPITATION IN LAYER K IN TEMPORARY VARIABLES
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
        IF (K  ==  KCT+1 .OR. BDD_START(I)) THEN
           BDDWT_K(I) = THDD_K(I)*EXK(I)  >   TM
        ELSE
          BDDWT_K(I) = BDDWT_KM1(I)
        END IF
          RAIN_TMP(I) = RAIN(I)
          SNOW_TMP(I) = SNOW(I)
!
!-----------------------------------------------------------------------
! DRY DESCENT FROM LAYER K TO K-1
!
! ENTRAINMENT CALCULATION
!-----------------------------------------------------------------------
!
          THDD_KM1(I) = (THDD_K(I)+(EKM14(I)*THE_K(I)) +                &
     &                  (1.0+EKM14(I))*EKM34(I)*THE_KM1(I)) /           &
     &                  ((1.0+EKM14(I))*(1.0+EKM34(I)))
          QDD_KM1(I) = (QDD_K(I)+(EKM14(I)*QE_K(I)) +                   &
     &                 (1.0+EKM14(I))*EKM34(I)*QE_KM1(I))/              &
     &                 ((1.0+EKM14(I))*(1.0+EKM34(I)))
      END DO
!
!
!----------------------------------------------------------------------
!  DRY DESCENT FOR TRACERS
!----------------------------------------------------------------------
!
      IF(L_TRACER)THEN
!
        DO KTRA=1,NTRA
          DO I=1,NPNTS
!
          TRADD_KM1(I,KTRA)=(TRADD_K(I,KTRA)+(EKM14(I)*                 &
     &                      TRAE_K(I,KTRA)) +                           &
     &                      (1.0+EKM14(I))*EKM34(I)*TRAE_KM1(I,KTRA))/  &
     &                      ((1.0+EKM14(I))*(1.0+EKM34(I)))
!
          END DO
        END DO
!
      END IF
!-----------------------------------------------------------------------
! UPDATE MASS FLUX  AND CALCULATE TEMPERATURE OF LAYER K-1
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
          FLX_DD_KM1(I) = FLX_DD_K(I)*(1.0+EKM34(I))*(1.0+EKM14(I))*    &
     &                (1.0-AMDETK(I))
!
          TDD_KM1(I) = THDD_KM1(I)*EXKM1(I)
      END DO
!
!-----------------------------------------------------------------------
! CALCULATE SUBSATURATION
! CALCULATE TEMPERATURE IF BROUGHT TO SATURATION
!-----------------------------------------------------------------------
!
! DEPENDS ON: satcal
       CALL SATCAL(NPNTS,TDD_KM1,THDD_KM1,P_KM1,QDDS,THDDS,             &
     &             K,EXKM1,QDD_KM1,THE_KM1)
!
      DO I=1,NPNTS
        BDDWT_KM1(I) = THDDS(I)*EXKM1(I)  >   TM
      END DO
!
!-----------------------------------------------------------------------
! CALCULATE CHANGE OF PHASE DUE TO DOWNDRAUGHT SATURATION TEMPERATURE
!-----------------------------------------------------------------------
!
! DEPENDS ON: crs_frzl
       CALL CRS_FRZL (NPNTS,RAIN,SNOW,THDD_KM1,EXKM1,FLX_DD_KM1,        &
     &                BDDWT_KM1)
!
      DO I=1,NPNTS
        TDD_KM1(I) = THDD_KM1(I)*EXKM1(I)
      END DO
!
!-----------------------------------------------------------------------
! RECALCULATE SUBSATURATION TEMPERATURE
!-----------------------------------------------------------------------
!
! DEPENDS ON: satcal
       CALL SATCAL(NPNTS,TDD_KM1,THDD_KM1,P_KM1,QDDS,THDDS,             &
     &             K,EXKM1,QDD_KM1,THE_KM1)
!
!-----------------------------------------------------------------------
! CALCULATE MOISTURE SUBSATURATION
!-----------------------------------------------------------------------
!
! DEPENDS ON: qsat
       CALL QSAT(QSATDD,TDD_KM1,P_KM1,NPNTS)
!
!-----------------------------------------------------------------------
! EVAPORATION CALCULATION AND ADJUSTMENT OF DOWNDRAUGHT TEMPERATURE
! AND MOISTURE
!-----------------------------------------------------------------------
!
! DEPENDS ON: devap
       CALL DEVAP (NPNTS,THDD_K,THDD_KM1,QDD_KM1,THDDS,QDDS,            &
     &             FLX_DD_KM1,EXK,EXKM1,QSATDD,RAIN,SNOW,               &
     &             DELPKM1,BDDWT_KM1,CCA,P_KM1)
!
!-----------------------------------------------------------------------
! CHECK IF PARCEL STILL NEGATIVELY BUOYANT SUCH THAT DOWNDRAUGHT CAN
! CONTINUE TO K-1
!-----------------------------------------------------------------------
!
! DEPENDS ON: termdd
       CALL TERMDD (NPNTS,BDD_START,THDD_KM1,QDD_KM1,THE_KM1,           &
     &              QE_KM1,K,B_DD_END,BDD_ON,PPN_MIX_DD)
!
!-----------------------------------------------------------------------
! CALCULATE THE EFFECT ON THE ENVIRONMENT IN LAYER K
!-----------------------------------------------------------------------
!
! DEPENDS ON: dd_env
       CALL DD_ENV (NPNTS,NP_FULL,THDD_K,THDD_KM1,QDD_K,QDD_KM1,THE_K,  &
     &              THE_KM1,QE_K,QE_KM1,DTHBYDT_K,DTHBYDT_KM1,DQBYDT_K, &
     &              DQBYDT_KM1,FLX_DD_K,FLX_DD_KM1,DELPK,DELPKM1,       &
     &              DELTD,DELQD,AMDETK,EKM14,B_DD_END,BDD_START,BDD_ON, &
     &              L_TRACER,NTRA,TRADD_K,TRADD_KM1,TRAE_K,             &
     &              TRAE_KM1,DTRABYDT_K,DTRABYDT_KM1,DELTRAD)
!
!-----------------------------------------------------------------------
! RESET DOWNDRAUGHT BIT VECTORS
!
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
       BDD_START(I) = .FALSE.
       IF (.NOT. BDD_ON(I)) THEN
         RAIN(I) = RAIN_TMP(I)
         SNOW(I) = SNOW_TMP(I)
       END IF
       IF (B_DD_END(I)) BDD_ON(I) = .FALSE.
      END DO
!
!-----------------------------------------------------------------------
! SWITCH POTENTIAL TEMPERATURE, MIXING RATIO, MOMENTUM, MASS FLUX
! AND TRACER READY FOR CALCULATION AT NEXT MODEL LAYER
!-----------------------------------------------------------------------
!
      IF (K >  2) THEN
        DO I=1,NPNTS
         IF (BDD_ON(I)) THEN
          THDD_K(I) = THDD_KM1(I)
          QDD_K(I) = QDD_KM1(I)
          FLX_DD_K(I) = FLX_DD_KM1(I)
         END IF
        END DO
!
       IF(L_TRACER)THEN
!
        DO KTRA=1,NTRA
          DO I=1,NPNTS
            IF(BDD_ON(I))THEN
              TRADD_K(I,KTRA) = TRADD_KM1(I,KTRA)
            END IF
          END DO
        END DO
!
       END IF
!
      END IF

      RETURN
      END SUBROUTINE DDRAUGHT
!
#endif
