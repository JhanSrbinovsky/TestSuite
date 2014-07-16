#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DD_ENV-------------------------------------------------
!LL
!LL  PURPOSE : CALCULATE THE EFFECT OF THE DOWNDRAUGHT
!LL            ON THE LARGE_SCALE ATMOSPHERE
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL   5.4   6/8/02    New deck created for version 4A of convection
!LL                   scheme. No CMT by downdraughts in this version.
!LL                                                      Gill Martin
!     6.2   03/02/05  Added section 5A R.A.Stratton
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
      SUBROUTINE DD_ENV (NPNTS,NP_FULL,THDD_K,THDD_KM1,QDD_K,QDD_KM1,   &
     &                   THE_K,THE_KM1,QE_K,QE_KM1,DTHBYDT_K,           &
     &                   DTHBYDT_KM1,DQBYDT_K,DQBYDT_KM1,FLX_DD_K,      &
     &                   FLX_DD_KM1,DELPK,DELPKM1,DELTD,DELQD,AMDETK,   &
     &                   EKM14,B_DD_END,BDD_START,BDD_ON,               &
     &                   L_TRACER,NTRA,TRADD_K,TRADD_KM1,TRAE_K,        &
     &                   TRAE_KM1,DTRABYDT_K,DTRABYDT_KM1,DELTRAD)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
      INTEGER NPNTS                 ! IN VECTOR LENGTH
!
      INTEGER NP_FULL               ! IN FULL VECTOR LENGTH
!
      INTEGER I,KTRA                ! LOOP COUNTERS
!
      INTEGER NTRA                  ! IN NUMBER OF TRACERS
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!-----------------------------------------------------------------------
!
      REAL THDD_K(NPNTS)            ! IN DOWNDRAUGHT POTENTIAL
                                    !    TEMPERATURE IN LAYER K (K)
!
      REAL THDD_KM1(NPNTS)          ! IN DOWNDRAUGHT POTENTIAL
                                    !    TEMPERATURE IN LAYER K-1 (K)
!
      REAL QDD_K(NPNTS)             ! IN DOWNDRAUGHT MIXING RATIO
                                    !    AT LAYER K (KG/KG)
!
      REAL QDD_KM1(NPNTS)           ! IN DOWNDRAUGHT MIXING RATIO
                                    !    AT LAYER K-1 (KG/KG)
!
      REAL TRADD_K(NP_FULL,NTRA)    ! IN DOWNDRAUGHT TRACER CONTENT
                                    !    IN LAYER K (KG/KG)
!
      REAL TRADD_KM1(NPNTS,NTRA)    ! IN DOWNDRAUGHT TRACER CONTENT
                                    !    AT LAYER K-1 (KG/KG)
!
      REAL THE_K(NPNTS)             ! IN POTENTIAL TEMPERATURE OF
                                    !    ENVIRONMENT IN LAYER K (K)
!
      REAL THE_KM1(NPNTS)           ! IN POTENTIAL TEMPERATURE OF
                                    !    ENVIRONMENT IN LAYER K-1 (K)
!
      REAL QE_K(NPNTS)              ! IN MIXING RATIO AT LAYER K (KG/KG)
!
      REAL QE_KM1(NPNTS)            ! IN MIXING RATIO AT LAYER K-1
                                    !    (KG/KG)
!
      REAL TRAE_K(NP_FULL,NTRA)     ! IN TRACER CONTENT OF
                                    !    ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL TRAE_KM1(NP_FULL,NTRA)   ! IN TRACER CONTENT OF ENVIRONMENT
                                    !    IN LAYER K-1 (KG/KG)
!
      REAL FLX_DD_K(NPNTS)          ! IN MASS FLUX IN LAYER K (PA/S)
!
      REAL FLX_DD_KM1(NPNTS)        ! IN MASS FLUX IN LAYER K-1 (PA/S)
!
      REAL DELPK(NPNTS)             ! IN DIFFERENCE IN PRESSURE ACROSS
                                    !    LAYER K (PA)
!
      REAL DELPKM1(NPNTS)           ! IN DIFFERENCE IN PRESSURE ACROSS
                                    !    LAYER K-1 (PA)
!
      REAL DELTD(NPNTS)             ! IN COOLING NECESSARY TO ACHIEVE
                                    !    SATURATION (K)
!
      REAL DELQD(NPNTS)             ! IN MOISTENING NECESSARY TO ACHIEVE
                                    !    SATURATION (KG/KG)
!
      REAL DELUD(NPNTS)             ! IN CHANGE TO ENVIRONMENT U DUE TO
                                    !    DOWNDRAUGHT FORMATION (M/S)
!
      REAL DELVD(NPNTS)             ! IN CHANGE TO ENVIRONMENT V DUE TO
                                    !    DOWNDRAUGHT FORMATION (M/S)
!
      REAL DELTRAD(NP_FULL,NTRA)    ! IN DEPLETION OF ENVIRONMENT TRACER
                                    !    DUE TO DOWNDRAUGHT FORMATION
!
      REAL AMDETK(NPNTS)            ! IN MIXING DETRAINMENT AT LEVEL K
                                    !    MULTIPLIED BY APPROPRIATE LAYER
                                    !    THICKNESS
!
      REAL EKM14(NPNTS)             ! IN EXNER RATIO AT LAYER K-1/4
!
      LOGICAL B_DD_END(NPNTS)       ! IN MASK FOR THOSE POINTS WHERE
                                    !    DOWNDRAUGHT IS TERMINATING
!
      LOGICAL BDD_START(NPNTS)      ! IN MASK FOR THOSE POINTS WHERE
                                    !    DOWNDRAUGHT IS STARTING
!
      LOGICAL BDD_ON(NPNTS)         ! IN MASK FOR THOSE POINTS WHERE
                                    !    DOWNDRAUGHT IS ON
!
      LOGICAL L_TRACER              ! IN SWITCH FOR INCLUSION OF TRACERS
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!-----------------------------------------------------------------------
!
      REAL DTHBYDT_K(NPNTS)         ! INOUT
                                    ! IN  INCREMENT TO POTENTIAL
                                    !     TEMPERATURE IN LAYER K (K/S)
                                    ! OUT UPDATED INCREMENT TO POTENTIAL
                                    !     TEMPERATURE IN LAYER K (K/S)
!
      REAL DTHBYDT_KM1(NPNTS)       ! INOUT
                                    ! IN  INCREMENT TO POTENTIAL
                                    !     TEMPERATURE AT LAYER K-1 (K/S)
                                    ! OUT UPDATED INCREMENT TO POTENTIAL
                                    !     TEMPERATURE AT LAYER K-1 (K/S)
!
      REAL DQBYDT_K(NPNTS)          ! INOUT
                                    ! IN  INCREMENT TO MIXING RATIO
                                    !     AT LAYER K (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MIXING
                                    !     RATIO AT LAYER K (KG/KG/S)
!
      REAL DQBYDT_KM1(NPNTS)        ! INOUT
                                    ! IN  INCREMENT TO MIXING RATIO AT
                                    !     LAYER K-1 (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MIXING
                                    !     RATIO AT LAYER K-1 (KG/KG/S)
!
      REAL DTRABYDT_K(NP_FULL,NTRA) ! INOUT
                                    ! IN  INCREMENT TO TRACER CONTENT
                                    !     IN LAYER K (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO TRACER
                                    !     CONTENT IN LAYER K (KG/KG/S)
!
      REAL DTRABYDT_KM1(NP_FULL,                                        &
                                    ! INOUT
     &                  NTRA)       ! IN  INCREMENT TO TRACER CONTENT
                                    !     AT LAYER K-1 (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO TRACER
                                    !     CONTENT AT LAYER K-1
                                    !     (KG/KG/S)
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!-----------------------------------------------------------------------
!
      REAL TEMPRY                   ! USED IN CALCULATIONS OF THE
                                    ! EFFECT ON THE ENVIRONMENT
!
!
!-----------------------------------------------------------------------
! CALCULATE THE EFFECT ON THE ENVIRONMENT IN LAYER K
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
       IF (BDD_ON(I)) THEN
!
!-----------------------------------------------------------------------
! SUBTRACT ENERGY USED TO FORM DOWNDRAUGHT
!-----------------------------------------------------------------------
!
       TEMPRY = FLX_DD_K(I)/DELPK(I)
       IF (BDD_START(I)) THEN
         DTHBYDT_K(I) = DTHBYDT_K(I)-TEMPRY*DELTD(I)
         DQBYDT_K(I) = DQBYDT_K(I)-TEMPRY*DELQD(I)
       END IF
!
!-----------------------------------------------------------------------
! EFFECT OF CONVECTION AND DOWNDRAUGHT UPON POTENTIAL TEMPERATURE OF
! LAYER K
!-----------------------------------------------------------------------
!
       DTHBYDT_K(I) = DTHBYDT_K(I) + TEMPRY * (                         &
      
     &          (1.0+EKM14(I)) * (1.0-AMDETK(I)) *                      &
                                                        ! COMPENSATING
     &           (THE_KM1(I)-THE_K(I))                                  &
                                                        ! SUBSIDENCE
     &        +                                                         &
     &          AMDETK(I)* (THDD_K(I)-THE_K(I))                         &
                                                        ! MIXING
     &        )                                         ! DETRAINMENT
!
!-----------------------------------------------------------------------
! EFFECT OF CONVECTION AND DOWNDRAUGHT UPON MIXING RATIO OF
! LAYER K
!-----------------------------------------------------------------------
!
       DQBYDT_K(I) = DQBYDT_K(I) + TEMPRY * (                           &
      
     &      (1.0+EKM14(I)) * (1.0-AMDETK(I)) *                          &
                                                     ! COMPENSATING
     &      (QE_KM1(I)-QE_K(I))                                         &
                                                     ! SUBSIDENCE
     &    +                                                             &
     &      AMDETK(I)* (QDD_K(I)-QE_K(I))                               &
                                                     ! MIXING
     &    )                                          ! DETRAINMENT
!
!
!-----------------------------------------------------------------------
! TERMINAL DETRAINMENT AND SUBSIDENCE IN TERMINAL LAYER
!-----------------------------------------------------------------------
!
         IF (B_DD_END(I)) THEN
           TEMPRY = FLX_DD_KM1(I)/DELPKM1(I)
           DTHBYDT_KM1(I) = DTHBYDT_KM1(I)+TEMPRY*                      &
     &                      (THDD_KM1(I)-THE_KM1(I))
           DQBYDT_KM1(I) = DQBYDT_KM1(I)+TEMPRY*(QDD_KM1(I)-QE_KM1(I))
!
         END IF
!
       END IF
      END DO
!
!-----------------------------------------------------------------------
! EFFECT OF CONVECTION AND DOWNDRAUGHT UPON TRACER CONTENT OF
! LAYER K
!-----------------------------------------------------------------------
!
      IF(L_TRACER)THEN
!
      DO KTRA=1,NTRA
       DO I=1,NPNTS
        IF(BDD_ON(I))THEN
!
        TEMPRY = FLX_DD_K(I)/DELPK(I)
        IF(BDD_START(I))THEN
        DTRABYDT_K(I,KTRA) = DTRABYDT_K(I,KTRA)-TEMPRY*DELTRAD(I,KTRA)
        END IF
        DTRABYDT_K(I,KTRA) = DTRABYDT_K(I,KTRA) + TEMPRY * (            &
      
     &   (1.0+EKM14(I)) * (1.0-AMDETK(I)) *                             &
                                                        ! COMPENSATING
     &   (TRAE_KM1(I,KTRA)-TRAE_K(I,KTRA))                              &
                                                        ! SUBSIDENCE
     &        +                                                         &
     &   AMDETK(I)* (TRADD_K(I,KTRA)-TRAE_K(I,KTRA))                    &
                                                        ! MIXING
     &        )                                         ! DETRAINMENT
!
!-----------------------------------------------------------------------
! TERMINAL DETRAINMENT OF TRACER
!-----------------------------------------------------------------------
!
           IF(B_DD_END(I))THEN
             TEMPRY = FLX_DD_KM1(I)/DELPKM1(I)
             DTRABYDT_KM1(I,KTRA)=DTRABYDT_KM1(I,KTRA)+TEMPRY*          &
     &                            (TRADD_KM1(I,KTRA)-TRAE_KM1(I,KTRA))
           END IF
!
        END IF
       END DO
      END DO
!
      END IF
!
      RETURN
      END SUBROUTINE DD_ENV
!
#endif
