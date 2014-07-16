
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE LAYER_CN-----------------------------------------------
!LL
!LL  PURPOSE : CALCULATES LAYER DEPENDENT CONSTANTS FOR LAYER K
!LL            -PRESSURE
!LL            -LAYER THICKNESS
!LL            -ENTRAINMENT COEFFICIENTS
!LL            -DETRAINMENT COEFFICIENTS
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL C.W. , D.G. <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
!LL
!LL  MODEL            MODIFICATION HISTORY:
!LL VERSION  DATE
!LL
!LL   5.4    6/8/02   New deck created for version 4A of convection
!LL                   scheme.               Gill Martin
!     6.2   09/12/05  Add adaptive detrainment option to
!                     convection scheme    Anna Maidens
!LL
!     6.2   03/02/05  Added section 5A. R A Stratton
!LL  PROGRAMMING STANDARDS :
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE LAYER_CN(K,NP_FIELD,NPNTS,NLEV,                        &
     &                    EXNER_LAYER_BOUNDARIES,                       &
     &                    EXNER_LAYER_CENTRES,                          &
     &                    P_LAYER_BOUNDARIES,                           &
     &                    P_LAYER_CENTRES,                              &
     &                    PSTAR,PK,PKP1,DELPK,DELPKP1,                  &
     &                    DELPKP12,EKP14,EKP34,AMDETK,EXK,EXKP1,        &
     &                    DELEXKP1,DELP_UV_K, DELP_UV_KP1,recip_PSTAR,  &
     &                    RHUM,L_SHALLOW,NTML,NTPAR                     &
     &                    ,CUMULUS,MDET_ON, ENT_ON, BCONV               &
     &                    ,THEK, QEK,QSEK, THEKP1,QEKP1,QSEKP1          &
     &                    ,THPK,QPK                                     &
     &                    ,BWK, BWKP1,EKM14                             &
     &                    ,ZK, ZKP12, ZKP1                              &
     &                    )

      IMPLICIT NONE
!
!----------------------------------------------------------------------
! MODEL CONSTANTS
!----------------------------------------------------------------------
!
! ENTCNST start

!  6.2   28/02/05  Modifications for STPH_RP
      ! coefficients used in calculation of entrainment rate
      REAL,PARAMETER:: AE1     = 1.0
      REAL,PARAMETER:: AE2     = 1.5
      REAL:: ENTCOEF
      COMMON/STPH_RP_1/entcoef
      REAL,PARAMETER:: SH_FAC  = 1.0

! ENTCNST end
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTER
!----------------------------------------------------------------------
!
      INTEGER NP_FIELD          ! IN FULL LENGTH OF DATA
!
      INTEGER NPNTS             ! IN VECTOR LENGTH
!
      INTEGER NLEV              ! IN NUMBER OF MODEL LEVELS
!
      INTEGER K                 ! IN PRESENT MODEL LAYER
!
      INTEGER I                 ! COUNTER FOR DO LOOPS
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      REAL EXNER_LAYER_BOUNDARIES(NP_FIELD,0:NLEV)
                                  ! IN EXNER FUNCTION AT LAYER

      REAL EXNER_LAYER_CENTRES(NP_FIELD,0:NLEV)
                                  ! IN EXNER FUNCTION AT LAYER

                                  ! BOUNDARIES STARTING AT LEVEL K-1/2
      REAL P_LAYER_CENTRES(NP_FIELD,0:NLEV)     ! IN PRESSURE.
!
      REAL P_LAYER_BOUNDARIES(NP_FIELD,0:NLEV)  ! IN PRESSURE.

      REAL PSTAR(NP_FIELD)      ! IN SURFACE PRESSURE (PA)
!
!
      REAL RHUM(NPNTS)         ! IN Relative humidity at level K

      INTEGER NTML(NPNTS)      ! IN Number of levels in the surface-base
!                              !    turbulently mixed layer

      INTEGER NTPAR(NPNTS)     ! IN Top level of initial parcel
!                              !    ascent in BL scheme.

      LOGICAL L_SHALLOW(NPNTS) ! IN logical for shallow convection

      LOGICAL, INTENT(IN) :: CUMULUS(NPNTS)   !logical for cumulus

      INTEGER, INTENT(IN) :: MDET_ON          !flag for adaptive
                                              !mixing detrainment
                                              !  on = 1, off = 0

      INTEGER, INTENT(IN) :: ENT_ON           !flag for adaptive
                                              !entrainment

      LOGICAL, INTENT(IN) :: BCONV(NPNTS)     !mask for points at
                                              !which
                                              !convection takes place

      REAL, INTENT(IN) :: THEK(NPNTS)         !theta for environment
                                              !on k

      REAL, INTENT(IN) :: QEK(NPNTS)          !q for environment on k

      REAL, INTENT(IN) :: QSEK(NPNTS)         !q sat for environment
                                              !on k

      REAL, INTENT(IN) :: THEKP1(NPNTS)       !theta for environment
                                              !on k+1

      REAL, INTENT(IN) :: QEKP1(NPNTS)        !q for environment on k+1

      REAL, INTENT(IN) :: QSEKP1(NPNTS)       !q sat for environment
                                              !on k+1

      REAL, INTENT(IN) :: THPK(NPNTS)         !parcel theta on k

      REAL, INTENT(IN) :: QPK(NPNTS)          !parcel q on k

      LOGICAL, INTENT(IN) :: BWK(NPNTS)       !mask for points where
                                              !condensate is liquid

      LOGICAL, INTENT(IN) :: BWKP1(NPNTS)     !mask for points where
                                              ! condensate is liquid

      REAL, INTENT(IN) :: EKM14(NPNTS)        !ek14 from previous
                                              ! pass i.e.
                                              !ek for  k - 1 + 1/4

      REAL, INTENT(IN) :: ZK(NPNTS)           !height on k

      REAL, INTENT(IN) :: ZKP12(NPNTS)        !height on k+ 1/2

      REAL, INTENT(IN) :: ZKP1(NPNTS)         !height on k+1

!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!----------------------------------------------------------------------
!
      REAL PK(NPNTS)            ! OUT PRESSURE AT LAYER K (PA)
!
      REAL PKP1(NPNTS)          ! OUT PRESSURE AT LAYER K+1 (PA)
!
      REAL DELPK(NPNTS)         ! OUT THICKNESS OF LAYER K (PA)
!
      REAL DELPKP1(NPNTS)       ! OUT THICHNESS OF LAYER K+1 (PA)
!
      REAL DELPKP12(NPNTS)      ! OUT THICKNESS BETWEEN LAYER K AND K+1
                                !     (PA)
      REAL DELP_UV_K(NPNTS)      ! OUT THICKNESS OF UV LAYER K (PA)
!
      REAL DELP_UV_KP1(NPNTS)   ! OUT THICKNESS OF UV LAYER K + 1(PA)
!
!
      REAL EKP14(NPNTS)         ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K+1/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
!
      REAL EKP34(NPNTS)         ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K+3/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
!
      REAL AMDETK(NPNTS)        ! OUT MIXING DETRAINMENT COEFFICIENT
                                !     AT LEVEL K MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
!
      REAL EXK(NPNTS)           ! EXNER FUNCTION AT LEVEL K
!
      REAL EXKP1(NPNTS)         ! EXNER FUNCTION AT LEVEL K+1
!
      REAL DELEXKP1(NPNTS)      ! DIFFERENCE IN EXNER FUNCTION
                                ! BETWEEN K+3/2 AND K+1/2
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!----------------------------------------------------------------------
!
      REAL AEKP14,AEKP34        ! USED IN CALCULATION OF ENTRAINMENT
                                ! RATE
!

      REAL                                                              &
     &    PU,PL
      REAL recip_PSTAR(NPNTS)   ! Reciprocal of pstar array
      REAL HSAT_EK(NPNTS)
      REAL HSAT_EKP1(NPNTS)     !saturated moist static energies
      REAL HSAT_EKP12(NPNTS)
      REAL H_PK(NPNTS)          !moist static energy
      REAL H_PKP12(NPNTS)
      REAL H_EK(NPNTS)
      REAL H_EKP1(NPNTS)
      REAL H_EKP12(NPNTS)
      REAL DELTA_HSAT_ENV(NPNTS)
                                !HSAT_EKP1 - HSAT_EK
      REAL DHPBYDP(NPNTS)
      REAL DHEBYDP(NPNTS)
      REAL RENT                 !parameter controlling fractional
                                !entrainment
      REAL EL                   !latent heat of condensation
                                !(/freezing)
      REAL EPSILONP14(NPNTS)    !fractional entrainment coef calc from
                                !hsat
      REAL EPSILONP34(NPNTS)    !fractional entrainment coef calc from
                                !hsat
      REAL EKP14_AD(NPNTS)      !EPS * LAYER THICKNESS FOR 1/4 LAYER
      REAL EKP34_AD(NPNTS)      !EPS * LAYER THICKNESS FOR 3/4 LAYER
      REAL AMDET_FAC            !factor to scale down AMDETK
                                !'cos eps is bigger
!*---------------------------------------------------------------------
!
!----------------------------------------------------------------------
! SET CONSTANT AE USED IN CALCULATION OF ENTARINMENT AND
! DETRAINMENT RATES DEPENDING UPON LEVEL
!----------------------------------------------------------------------
!
        AEKP14 = AE2
        AEKP34 = AE2
!
       IF (ENT_ON  ==  1 ) THEN
       RENT=0.5

!bit compare -possible fix, initialise everything in sight!
!except ekm14, obviously
      DO I=1, NPNTS
        HSAT_EK(I)=0.0
        HSAT_EKP1(I)=0.0
        HSAT_EKP12(I)=0.0
        H_PK(I)=0.0
        H_PKP12(I)=0.0
        H_EK(I)=0.0
        H_EKP1(I)=0.0
        H_EKP12(I)=0.0
        DELTA_HSAT_ENV(I)=0.0
        DHPBYDP(I)=0.0
        DHEBYDP(I)=0.0
        EPSILONP14(I)=0.0
        EPSILONP34(I)=0.0
        EKP14_AD(I)=0.0
        EKP34_AD(I)=0.0
      END DO
         DO  I=1,NPNTS
!
! ---------------------------------------------------------------------
!  CALCULATE PK AND DELPK - IF K = 1 (LOWEST MODEL LAYER) THEN
!  VALUES FOR PREVIOUS PASS THROUGH ROUTINE AT (K-1)+1 ARE TAKEN
! ---------------------------------------------------------------------
!
           IF(K == 2)THEN
              PK(I) = P_LAYER_CENTRES(I,K)
! note minus sign so be careful using it in discretizing
!derivatives
              DELPK(I) = - (P_LAYER_BOUNDARIES(I,K) -                   &
     &                     P_LAYER_BOUNDARIES(I,K-1))
           ELSE
             PK(I) =  PKP1(I)
             DELPK(I) = DELPKP1(I)
           END IF

!
! ---------------------------------------------------------------------
!  CALCULATE PKP1, DELPKP1 AND DELPK+1/2
! ---------------------------------------------------------------------
           PKP1(I) = P_LAYER_CENTRES(I,K+1)
!note minus sign so be careful using it in discretizing
!derivatives
              DELPKP1(I) = - (P_LAYER_BOUNDARIES(I,K+1) -               &
     &                     P_LAYER_BOUNDARIES(I,K))
!note opposite sign convention to other DELPetc
           DELPKP12(I) = PK(I) - PKP1(I)
!
!Calculate hsat in case where adaptive entrainment is switched on
!Already inside loop over I
             IF(BWK(I)) THEN
               EL = LC
             ELSE
               EL = LC+LF
             ENDIF
             HSAT_EK(I)=CP * EXNER_LAYER_CENTRES(I,K) * THEK(I)         &
     &                   + EL * QSEK(I) + G*ZK(I)
             H_EK(I)=CP * EXNER_LAYER_CENTRES(I,K) * THEK(I)            &
     &                   + EL * QEK(I) + G*ZK(I)
             H_PK(I)=CP * EXNER_LAYER_CENTRES(I,K) * THPK(I)            &
     &                   + EL * QPK(I) + G*ZK(I)


             IF(BWKP1(I)) THEN
               EL = LC
             ELSE
               EL = LC+LF
             ENDIF
             HSAT_EKP1(I)=CP * EXNER_LAYER_CENTRES(I,K+1) * THEKP1(I)   &
     &                   + EL * QSEKP1(I) + G*ZKP1(I)
             H_EKP1(I)=CP * EXNER_LAYER_CENTRES(I,K+1) * THEKP1(I)      &
     &                + EL * QEKP1(I) + G*ZKP1(I)


             DELTA_HSAT_ENV(I)=(HSAT_EKP1(I)-HSAT_EK(I))


!version below has stability criterion, but doesn't really help
!          DELTA_HSAT_ENV(I)=MIN((HSAT_EKP1(I)-HSAT_EK(I)), 0.0)
!Not sure precisely what to do with this- if DELTA_HSAT_ENV
!before correction was positive, convection shouldn't be taking
!place at all, never mind the value of the entrainment coefficient.


               HSAT_EKP12(I) = HSAT_EK(I)+                              &
     &                  (P_LAYER_BOUNDARIES(I,K)- P_LAYER_CENTRES(I,K))/&
     &                  (P_LAYER_CENTRES(I,K+1)- P_LAYER_CENTRES(I,K)) *&
     &             DELTA_HSAT_ENV(I)



       !NB already inside loop over NPNTS

          DHPBYDP(I) = (1- RENT) * DELTA_HSAT_ENV(I)/                   &
     &               (P_LAYER_CENTRES(I,K+1)- P_LAYER_CENTRES(I,K))
          H_PKP12(i) = H_PK(I) + (P_LAYER_BOUNDARIES(I,K) -             &
     &                P_LAYER_CENTRES(I,K)) * DHPBYDP(I)

          DHEBYDP(I) = (H_EKP1(I) - H_EK(I))/                           &
     &                 (P_LAYER_CENTRES(I,K+1)- P_LAYER_CENTRES(I,K))
          H_EKP12(I) = H_EK(I) + (P_LAYER_BOUNDARIES(I,K) -             &
     &                P_LAYER_CENTRES(I,K)) * DHEBYDP(I)



! ---------------------------------------------------------------------
!  CALCULATE DELP_UV_K and DELP_UV_KP1
!  NB: Ensure positive by taking - of difference
! ---------------------------------------------------------------------
!
           DELP_UV_K(I) = P_LAYER_CENTRES(I,K-1) - P_LAYER_CENTRES(I,K)
           DELP_UV_KP1(I) = P_LAYER_CENTRES(I,K) -                      &
     &                                   P_LAYER_CENTRES(I,K+1)
!
! ---------------------------------------------------------------------
!  CALCULATE EXNER FUNCTIONS AT MID-LAYES K AND K+1, AND
!  DIFFERENCE OF EXNER FUNCTION ACROSS LAYER K
! ---------------------------------------------------------------------
!
           IF(K == 2)THEN
              EXK(I) = EXNER_LAYER_CENTRES(I,K)
           ELSE
             EXK(I) = EXKP1(I)
           END IF
             EXKP1(I) = EXNER_LAYER_CENTRES(I,K+1)
           DELEXKP1(I) = EXNER_LAYER_BOUNDARIES(I,K)-                   &
     &                   EXNER_LAYER_BOUNDARIES(I,K+1)
!
! ---------------------------------------------------------------------
!  CALCULATE ENTRAINMENT AND MIXING DETRAINMENT COEFFICIENTS
! ---------------------------------------------------------------------
!
!
! ---------------------------------------------------------------------
!  CALCULATE ENTRAINMENT COEFFICIENTS MULTIPLIED BY
!  APPROPRIATE LAYER THICKNESS
!
!  UM DOCUMENTATION PAPER 27
!  SECTION (2C), EQUATION(14)
! ---------------------------------------------------------------------
!
           IF (L_SHALLOW(I) .AND. K  >   NTML(I) ) THEN

              EKP14(I) = (PK(I) - P_LAYER_BOUNDARIES(I,K))              &
     &            * 0.03 * EXP( -1.0*( P_LAYER_BOUNDARIES(I,NTML(I))    &
     &                            - PK(I) )                             &
     &          / ( P_LAYER_BOUNDARIES(I,NTML(I))                       &
     &            - P_LAYER_BOUNDARIES(I,NTPAR(I)) ) )                  &
     &            / ( P_LAYER_BOUNDARIES(I,NTML(I))                     &
     &              - P_LAYER_BOUNDARIES(I,NTPAR(I)) )

           ELSE

            EKP14(I) = ENTCOEF * AEKP14 * PK(I) *                       &
     &                (PK(I) - P_LAYER_BOUNDARIES(I,K))*                &
     &                recip_PSTAR(I) * recip_PSTAR(I)

           ENDIF


           IF (L_SHALLOW(I) .AND. K  ==  NTML(I)                        &
     &                          ) EKP14(I) = 0.0

           IF (L_SHALLOW(I) .AND. K  >=  NTML(I) ) THEN

              EKP34(I) = (P_LAYER_BOUNDARIES(I,K) - PKP1(I))            &
     &             * 0.03 * EXP( -1.0*( P_LAYER_BOUNDARIES(I,NTML(I))   &
     &              - P_LAYER_BOUNDARIES(I,K) )                         &
     &          / ( P_LAYER_BOUNDARIES(I,NTML(I))                       &
     &            - P_LAYER_BOUNDARIES(I,NTPAR(I)) ) )                  &
     &            / ( P_LAYER_BOUNDARIES(I,NTML(I))                     &
     &              - P_LAYER_BOUNDARIES(I,NTPAR(I)) )

           ELSE

         EKP34(I) = ENTCOEF * AEKP34 * (P_LAYER_BOUNDARIES(I,K)) *      &
     &                (P_LAYER_BOUNDARIES(I,K)  - PKP1(I)) *            &
     &                recip_PSTAR(I) * recip_PSTAR(I)

           ENDIF

        !-------------------
        !calculate epsilon14 and epsilon 34
        !-------------------
           IF(BCONV(I)) THEN
             IF((H_PK(I) - H_EK(I))  >   0.0) THEN
               EPSILONP14(I)=  DHPBYDP(I)/                              &
     &            MAX((H_PK(I) - H_EK(I)),500.0)
          !500 chosen with view to not letting buoyancy drop too much
         !approx 0.2 cp + 0.2/1000 L
             ELSE
               EPSILONP14(I)= 0.0
             END IF
             IF((H_PKP12(I) - H_EKP12(I))  >   0.0) THEN
               EPSILONP34(I)=   DHPBYDP(I)/                             &
     &            MAX((H_PKP12(I) - H_EKP12(I)),500.0)
             ELSE
               EPSILONP34(I)= 0.0
             END IF


           IF(((H_PK(I) - H_EK(I))  >   500)                            &
     &       .AND.(DHPBYDP(I)  >   0.000000001)) THEN
              EKP14_AD(I) = MIN((EPSILONP14(I) *                        &
     &                (PK(I) - P_LAYER_BOUNDARIES(I,K))),1.0)
           ELSE
             EKP14_AD(I) = 0.0
           END IF

           IF(((H_PKP12(I) - H_EKP12(I))  >   0.0000000001)             &
     &       .AND. (DHPBYDP(I)  >   0.0000000001)) THEN
             EKP34_AD(I) = MIN((EPSILONP34(I) *                         &
     &                (P_LAYER_BOUNDARIES(I,K)  - PKP1(I))),1.0)


!Check on sign of difference
             IF((((EKP14_AD(I) - EKM14(I)) <   0.0) .AND.               &
     &           ((EKP34_AD(I) - EKP14_AD(I)) >   0.0))                 &
     &          .OR.                                                    &
     &          (((EKP14_AD(I) - EKM14(I)) >   0.0) .AND.               &
     &           ((EKP34_AD(I) - EKP14_AD(I)) <   0.0))) THEN

               EKP34_AD(I) = EKP14_AD(I)

             END IF

           ELSE
             EKP34_AD(I) = 0.0
           END IF



!       reset values of entrainment coefficients to adaptive versions
!        EKP14(I)=MAX(EKP14_AD(I),0.0)
!        EKP34(I)=MAX(EKP34_AD(I),0.0)
!Try an average of the adaptive and original versions
!but using GR as minimum value below which  entrainment doesn't fall
           IF(EKP14_AD(I)  >   0.0) THEN
             EKP14(I)=0.5*(EKP14_AD(I)+EKP14(I))
           ENDIF
           IF(EKP34(I)  >   0.0) THEN
             EKP34(I)=0.5*(EKP34_AD(I)+EKP34(I))
           ENDIF
         END IF

!
! ---------------------------------------------------------------------
!  CALCULATE MIXING DETRAINMENT COEFFICIENT MULTIPLIED BY
!  APPROPRIATE LAYER THICKNESS
!
!  UM DOCUMENTATION PAPER 27
!  SECTION (2C), EQUATION(15)
! ---------------------------------------------------------------------
!
           AMDET_FAC=1
           IF(K == 1)THEN
             AMDETK(I) = 0.0
           ELSE IF (L_SHALLOW(I) .AND. K  >=  NTML(I) ) THEN
              IF (RHUM(I)  <=  0.85) THEN
                AMDETK(I) = (1.0 + 0.3)*(EKP14(I) + EKP34(I))
              ELSE
                AMDETK(I) = (1.0 + (0.3/0.15)*(1.0-RHUM(I)))            &
     &                                    *(EKP14(I) + EKP34(I))
                IF (RHUM(I)  >   1.0)                                   &
     &                   AMDETK(I) = 1.0*(EKP14(I) + EKP34(I))
              ENDIF
           ELSE IF ((.NOT. L_SHALLOW(I)) .AND. CUMULUS(I)               &
     &            .AND. K  >=  NTML(I) .AND. MDET_ON  ==  1) THEN
             IF (RHUM(I)  <=  1.0) THEN
               AMDETK(I) = AMDET_FAC* (EKP14(I) + EKP34(I))*            &
     &                         (1-RHUM(I))
             ELSE
                  AMDETK(I) = 0
             ENDIF
           ELSE
             AMDETK(I) = (EKP14(I) + EKP34(I)) * (1.0-1.0/AEKP34)
           ENDIF

         END DO !npnts
       ELSE     ! end of adaptive entrainment calculation
! calculate coefficients without adaptive entrainment
      DO 10 I=1,NPNTS
!L
!L---------------------------------------------------------------------
!L CALCULATE PK AND DELPK - IF K = 1 (LOWEST MODEL LAYER) THEN
!L VALUES FOR PREVIOUS PASS THROUGH ROUTINE AT (K-1)+1 ARE TAKEN
!L---------------------------------------------------------------------
!L
        IF(K == 2)THEN
           PK(I) = P_LAYER_CENTRES(I,K)
           DELPK(I) = - (P_LAYER_BOUNDARIES(I,K) -                      &
     &                  P_LAYER_BOUNDARIES(I,K-1))
        ELSE
          PK(I) =  PKP1(I)
          DELPK(I) = DELPKP1(I)
        END IF
!L
!L---------------------------------------------------------------------
!L CALCULATE PKP1, DELPKP1 AND DELPK+1/2
!L---------------------------------------------------------------------
        PKP1(I) = P_LAYER_CENTRES(I,K+1)
           DELPKP1(I) = - (P_LAYER_BOUNDARIES(I,K+1) -                  &
     &                  P_LAYER_BOUNDARIES(I,K))
        DELPKP12(I) = PK(I) - PKP1(I)
!L
!L---------------------------------------------------------------------
!L CALCULATE DELP_UV_K and DELP_UV_KP1
!L NB: Ensure positive by taking - of difference
!L---------------------------------------------------------------------
!L
        DELP_UV_K(I) = P_LAYER_CENTRES(I,K-1) - P_LAYER_CENTRES(I,K)
        DELP_UV_KP1(I) = P_LAYER_CENTRES(I,K) - P_LAYER_CENTRES(I,K+1)
!L
!L---------------------------------------------------------------------
!L CALCULATE EXNER FUNCTIONS AT MID-LAYES K AND K+1, AND
!L DIFFERENCE OF EXNER FUNCTION ACROSS LAYER K
!L---------------------------------------------------------------------
!L
        IF(K == 2)THEN
          EXK(I) = EXNER_LAYER_CENTRES(I,K)
        ELSE
          EXK(I) = EXKP1(I)
        END IF
          EXKP1(I) = EXNER_LAYER_CENTRES(I,K+1)
        DELEXKP1(I) = EXNER_LAYER_BOUNDARIES(I,K)-                      &
     &                EXNER_LAYER_BOUNDARIES(I,K+1)
!L
!L---------------------------------------------------------------------
!L CALCULATE ENTRAINMENT AND MIXING DETRAINMENT COEFFICIENTS
!L---------------------------------------------------------------------
!L
!L
!L---------------------------------------------------------------------
!L CALCULATE ENTRAINMENT COEFFICIENTS MULTIPLIED BY
!L APPROPRIATE LAYER THICKNESS
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (2C), EQUATION(14)
!L---------------------------------------------------------------------
!L
        IF (L_SHALLOW(I) .AND. K  >   NTML(I) ) THEN

           EKP14(I) = (PK(I) - P_LAYER_BOUNDARIES(I,K))                 &
     &         * 0.03 * EXP( -1.0*( P_LAYER_BOUNDARIES(I,NTML(I))       &
     &                         - PK(I) )                                &
     &       / ( P_LAYER_BOUNDARIES(I,NTML(I))                          &
     &         - P_LAYER_BOUNDARIES(I,NTPAR(I)) ) )                     &
     &         / ( P_LAYER_BOUNDARIES(I,NTML(I))                        &
     &           - P_LAYER_BOUNDARIES(I,NTPAR(I)) )

        ELSE

         EKP14(I) = ENTCOEF * AEKP14 * PK(I) *                          &
     &             (PK(I) - P_LAYER_BOUNDARIES(I,K))*                   &
     &             recip_PSTAR(I) * recip_PSTAR(I)

        ENDIF


        IF (L_SHALLOW(I) .AND. K  ==  NTML(I)                           &
     &                       ) EKP14(I) = 0.0

        IF (L_SHALLOW(I) .AND. K  >=  NTML(I) ) THEN

           EKP34(I) = (P_LAYER_BOUNDARIES(I,K) - PKP1(I))               &
     &          * 0.03 * EXP( -1.0*( P_LAYER_BOUNDARIES(I,NTML(I))      &
     &           - P_LAYER_BOUNDARIES(I,K) )                            &
     &       / ( P_LAYER_BOUNDARIES(I,NTML(I))                          &
     &         - P_LAYER_BOUNDARIES(I,NTPAR(I)) ) )                     &
     &         / ( P_LAYER_BOUNDARIES(I,NTML(I))                        &
     &           - P_LAYER_BOUNDARIES(I,NTPAR(I)) )

        ELSE

      EKP34(I) = ENTCOEF * AEKP34 * (P_LAYER_BOUNDARIES(I,K)) *         &
     &             (P_LAYER_BOUNDARIES(I,K)  - PKP1(I)) *               &
     &             recip_PSTAR(I) * recip_PSTAR(I)

        ENDIF


!L
!L---------------------------------------------------------------------
!L CALCULATE MIXING DETRAINMENT COEFFICIENT MULTIPLIED BY
!L APPROPRIATE LAYER THICKNESS
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (2C), EQUATION(15)
!L---------------------------------------------------------------------
!L
        IF(K == 1)THEN
          AMDETK(I) = 0.0
        ELSE IF (L_SHALLOW(I) .AND. K  >=  NTML(I) ) THEN
           IF (RHUM(I)  <=  0.85) THEN
             AMDETK(I) = (1.0 + 0.3)*(EKP14(I) + EKP34(I))
           ELSE
             AMDETK(I) = (1.0 + (0.3/0.15)*(1.0-RHUM(I)))               &
     &                                 *(EKP14(I) + EKP34(I))
             IF (RHUM(I)  >   1.0)                                      &
     &                AMDETK(I) = 1.0*(EKP14(I) + EKP34(I))
           ENDIF
        ELSE IF ((.NOT. L_SHALLOW(I)) .AND. CUMULUS(I)                  &
     &         .AND. K  >=  NTML(I) .AND. MDET_ON  ==  1) THEN
          IF (RHUM(I)  <=  1.0) THEN
            AMDETK(I) = (EKP14(I) + EKP34(I))*                          &
     &                      (1-RHUM(I))
          ELSE
               AMDETK(I) = 0
          ENDIF
        ELSE
          AMDETK(I) = (EKP14(I) + EKP34(I)) * (1.0-1.0/AEKP34)
        ENDIF

 10   CONTINUE
       END IF ! entrainment off condition


!
      RETURN
      END SUBROUTINE LAYER_CN
