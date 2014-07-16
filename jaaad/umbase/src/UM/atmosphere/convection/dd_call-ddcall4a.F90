#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DD_CALL------------------------------------------------
!LL
!LL  PURPOSE : CALCULATE INITIAL DOWNDRAUGHT MASSFLUX
!LL
!LL            RESET EN/DETRAINMENT RATES FOR DOWNDRAUGHT
!LL
!LL            COMPRESS/EXPAND VARIABLES
!LL
!LL            INITIALISE DOWNDRAUGHT
!LL
!LL            CALL DOWNDRAUGHT ROUTINE
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
      SUBROUTINE DD_CALL (NP_FIELD,NPNTS,KCT,THP,QP,THE,QE,DTHBYDT,     &
     &                    DQBYDT,FLX,PSTAR,                             &
     &                    P_LAYER_BOUNDARIES,                           &
     &                    P_LAYER_CENTRES,                              &
     &                    EXNER_LAYER_BOUNDARIES,                       &
     &                    EXNER_LAYER_CENTRES,                          &
     &                    PRECIP,RAIN,SNOW,RAIN_3D,SNOW_3D,ICCB,ICCT,   &
     &                    BWATER,BTERM,BGMK,TIMESTEP,CCA,NTERM,         &
     &                    L_TRACER,NTRA,                                &
     &                    TRAP,TRAE,DTRABYDT,NLEV,TRLEV,recip_pstar,    &
     &                    DD_FLUX,FLG_DWN_FLX,ENTRAIN_DWN,FLG_ENTR_DWN, &
     &                    DETRAIN_DWN,FLG_DETR_DWN,L_SHALLOW,INDEX1     &
     &                    )

      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
!
      INTEGER I,KTRA,I2          ! LOOP COUNTERS
!
      INTEGER K                  ! PRESENT MODEL LAYER
!
      INTEGER NPNTS              ! IN NUMBER OF POINTS
!
      INTEGER NDD,NTERM          ! COMPRESSED VECTOR LENGTH FOR
                                 ! DOWNDRAUGHT CALCULATION
!
      INTEGER NP_FIELD           ! IN FULL VECTOR LENGTH
!
      INTEGER NDDON_TMP          ! NUMBER OF POINTS WITH ACTIVE
                                 ! DOWNDRAUGHT
!
      INTEGER NTRA               ! NUMBER OF TRACER VARIABLES
!
      INTEGER NLEV               ! NUMBER OF MODEL LEVELS
!
      INTEGER TRLEV              ! NUMBER OF TRACER LEVELS
!
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!-----------------------------------------------------------------------
!
      INTEGER KCT                 ! IN CONVECTIVE CLOUD TOP LAYER
!
      REAL P_LAYER_CENTRES(NP_FIELD,0:KCT+2)

      REAL P_LAYER_BOUNDARIES(NP_FIELD,0:KCT+1)
!
      REAL EXNER_LAYER_BOUNDARIES(NP_FIELD,0:KCT+1)
      REAL EXNER_LAYER_CENTRES(NP_FIELD,0:KCT+1)
                                 ! IN EXNER FUNCTION AT LAYER BOUNDARIES
                                 !    STARTING AT LEVEL K-1/2
!
      REAL THP(NPNTS,KCT+1)      ! IN POTENTIAL TEMPERATURE OF
                                 !    PARCEL (K)
!
      REAL QP(NPNTS,KCT+1)       ! IN MODEL MIXING RATIO (KG/KG)
!
      REAL TRAP(NPNTS,NLEV,                                             &
                                ! IN PARCEL TRACER (KG/KG)
     &          NTRA)
!
      REAL THE(NP_FIELD,KCT+1)   ! IN MODEL ENVIRONMENTAL POTENTIAL
                                 !    TEMPERATURE (K)
!
      REAL QE(NP_FIELD,KCT+1)    ! IN ENVIRONMENT MIXING RATIO
                                 !    (KG/KG)
!
      REAL TRAE(NP_FIELD,TRLEV,                                         &
                                 ! IN ENVIRONMENT TRACER (KG/KG)
     &          NTRA)
!
      REAL FLX(NPNTS,KCT+1)      ! IN CONVECTIVE MASSFLUX (PA/S)
!
      REAL PSTAR(NP_FIELD)       ! IN SURFACE PRESSURE (PA)
!
      REAL PRECIP(NPNTS,KCT+1)   ! IN PRECIPITATION ADDED WHEN
                                 !    DESCENDING FROM LAYER K TO K-1
                                 !    (KG/M**2/S)
!
      INTEGER ICCB(NP_FIELD)     ! IN CLOUD BASE LEVEL
!
      INTEGER ICCT(NP_FIELD)     ! IN CLOUD TOP LEVEL
!
      REAL CCA(NP_FIELD)         ! IN CONVECTIVE CLOUD AMOUNT
!
      LOGICAL L_SHALLOW(NPNTS)   ! IN Logical flag for shallow Cu

      LOGICAL BWATER(NPNTS,2:KCT+1)!IN  MASK FOR THOSE POINTS AT WHICH
                                   !     CONDENSATE IS WATER IN LAYER K
!
      LOGICAL BTERM(NPNTS)       ! IN MASK FOR THOSE POINTS WHERE
                                 !    UPDRAUGHT IS TERMINATING
!
      LOGICAL BGMK(NPNTS)        ! IN MASK FOR POINTS WHERE PARCEL IN
                                 !    LAYER K IS SATURATED
!
      LOGICAL L_TRACER           ! IN SWITCH FOR INCLUSION OF TRACERS
!
                                 !    MOMENTUM TRANSPORTS
      LOGICAL FLG_DWN_FLX        ! STASH FLAG FOR DOWNDRAUGHT MASS FLUX
!
      LOGICAL FLG_ENTR_DWN       ! STASH FLAG FOR DOWNDRAUGHT ENTRAINMNT
!
      LOGICAL FLG_DETR_DWN       ! STASH FLAG FOR DOWNDRAUGHT DETRANMNT

      REAL TIMESTEP
      REAL recip_PSTAR(NPNTS)    ! Reciprocal of pstar array
!
      INTEGER INDEX1(NTERM) ! IN INDEX OF DOWNDRAUGHTS POSSIBLE
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT AND OUTPUT
!-----------------------------------------------------------------------
!
      REAL DTHBYDT(NP_FIELD,KCT+1) ! INOUT
                                   ! IN  INCREMENT TO MODEL POTENTIAL
                                   !     TEMPERATURE (K/S)
                                   ! OUT UPDATED INCREMENT TO MODEL
                                   !     POTENTIAL TEMPERATURE (K/S)
!
      REAL DQBYDT(NP_FIELD,KCT+1)  ! INOUT
                                   ! IN  INCREMENT TO MODEL MIXING
                                   !     RATIO (KG/KG/S)
                                   ! OUT UPDATED INCREMENT TO MODEL
                                   !     MIXING RATIO (KG/KG/S)
!
      REAL DTRABYDT(NPNTS,NLEV,                                         &
                                  ! INOUT
     &              NTRA)          ! IN  INCREMENT TO MODEL
                                   !     TRACER (KG/KG/S)
                                   ! OUT UPDATED INCREMENT TO
                                   !     MODEL TRACER (KG/KG/S)
!
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!-----------------------------------------------------------------------
!
      REAL RAIN(NP_FIELD)         ! OUT RAINFALL AT SURFACE (KG/M**2/S)
!
      REAL SNOW(NP_FIELD)         ! OUT SNOWFALL AT SURFACE (KG/M**2/S)
!
      REAL RAIN_3D(NP_FIELD,NLEV) ! OUT RAINFALL FLUX PROFILE (KG/M**2/S)
!
      REAL SNOW_3D(NP_FIELD,NLEV) ! OUT SNOWFALL FLUX PROFILE (KG/M**2/S)

      REAL DD_FLUX(NP_FIELD,KCT+1) ! OUT DOWN DRAUGHT MASS FLUX
!
      REAL ENTRAIN_DWN(NP_FIELD,KCT+1) ! OUT FRACTIONAL ENTRAINMENT
                                       ! RATE FOR DOWN DRAUGHT
!
      REAL DETRAIN_DWN(NP_FIELD,KCT+1) ! OUT FRACTIONAL DETRAINMENT
!
!-----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!-----------------------------------------------------------------------
!
!
      REAL EXNER_KM12_C(NTERM)   ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K
!
      REAL EXNER_KP12_C(NTERM)   ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K+1
!
      REAL EXNER_KM32_C(NTERM)   ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K-1
!
      REAL PK(NTERM)             ! PRESSURE OF LAYER K (PA)
!
      REAL P_KM1(NTERM)          ! PRESSURE OF LAYER K-1 (PA)
!
      REAL EXK(NTERM)            ! EXNER RATIO FOR LAYER K
!
      REAL EXKM1(NTERM)          ! EXNER RATIO FOR LAYER K-1
!
      REAL DELPK(NTERM)          ! PRESSURE DIFFERENCE ACROSS LAYER K
                                 ! (PA)
!
      REAL DELPKM1(NTERM)        ! PRESSURE DIFFERENCE ACROSS
                                 ! LAYER K-1 (PA)
!
      REAL AMDETK(NTERM)         ! MIXING DETRAINMENT AT LEVEL K
                                 ! MULTIPLIED BY APPROPRIATE LAYER
                                 ! THICKNESS
!
      REAL EKM14(NTERM)          ! EXNER RATIO AT LAYER K-1/4
!
      REAL EKM34(NTERM)          ! EXNER RATIO AT LAYER K-3/4
!
      LOGICAL BWATER_K_C(NTERM)  ! COMPRESSED MASK FOR THOSE
                                 ! POINTS AT WHICH CONDENSATE
                                 ! IS WATER IN LAYER K
!
      REAL PRECIP_K_C(NTERM)     ! COMPRESSED PRECIPITATION
                                 ! ADDED WHEN DESCENDING FROM
                                 ! LAYER K TO K-1 (KG/M**2/S)
!
      REAL Q_K_C(NTERM)          ! COMPRESSED PARCEL MIXING RATIO
                                 ! OF LAYER K (KG/KG)
!
      REAL TH_K_C(NTERM)         ! COMPRESSED PARCEL POTENTIAL
                                 ! TEMPERATURE OF LAYER K (K)
!
      REAL TRA_K_C(NTERM,NTRA)   ! COMPRESSED PARCEL TRACER IN
                                 ! LAYER K (KG/KG)
!
      REAL PSTAR_C(NTERM)        ! COMPRESSED SURFACE PRESSURE (PA)
!
      REAL recip_PSTAR_C(NTERM)  ! Reciprocal of comp. pstar array
!
      REAL P_LAYER_CENTRES_C(NTERM,0:KCT+2)
      REAL P_LAYER_BOUNDARIES_C(NTERM,0:KCT+1)
      REAL EXNER_LAYER_CENTRES_C(NTERM,0:KCT+1)

!
      INTEGER ICCB_C(NTERM)      ! COMPRESSED CLOUD BASE LEVEL
!
      REAL DTHBYDT_K_C(NTERM)    ! COMPRESSED INCREMENT TO MODEL
                                 ! POTENTIAL TEMPERATURE OF LAYER K
                                 ! (K/S)
!
      REAL DTHBYDT_KM1_C(NTERM)  ! COMPRESSED INCREMENT TO MODEL
                                 ! POTENTIAL TEMPERATURE OF LAYER K-1
                                 ! (K/S)
!
      REAL DQBYDT_K_C(NTERM)     ! COMPRESSED INCREMENT TO MODEL
                                 ! MIXING RATIO OF LAYER K (KG/KG/S)
!
      REAL DQBYDT_KM1_C(NTERM)   ! COMPRESSED INCREMENT TO MODEL
                                 ! MIXING RATIO OF LAYER K-1 (KG/KG/S)
!
      REAL DTRA_K_C(NTERM,NTRA)  ! COMPRESSED INCREMENT TO MODEL
                                 ! TRACER OF LAYER K (KG/KG/S)
!
      REAL DTRA_KM1_C(NTERM,NTRA)! COMPRESSED INCREMENT TO MODEL
                                 ! TRACER OF LAYER K-1 (KG/KG/S)
!
      REAL DELTD(NTERM)          ! COOLING NECESSARY TO
                                 ! ACHIEVE SATURATION (K)
!
      REAL DELQD(NTERM)          ! MOISTENING NECESSARY TO
                                 ! ACHIEVE SATURATION (KG/KG)
!
      REAL DELTRAD(NTERM,NTRA)   ! DEPLETION OF ENVIRONMENT TRACER
                                 ! DUE TO DOWNDRAUGHT FORMATION (KG/KG)
!
      REAL QDD_K(NTERM)          ! MIXING RATIO OF DOWNDRAUGHT IN
                                 ! LAYER K (KG/KG)
!
      REAL THDD_K(NTERM)         ! MODEL POTENTIAL TEMPERATURE
                                 ! OF DOWNDRAUGHT IN LAYER K (K)
!
      REAL TRADD_K(NTERM,NTRA)   ! MODEL TRACER OF DOWNDRAUGHT
                                 ! IN LAYER K (KG/KG)
!
      REAL FLX_DD_K(NPNTS)       ! DOWNDRAUGHT INITIAL MASS FLUX
                                 ! (PA/S)
!
      REAL FLX_DD_K_C(NTERM)     ! COMPRESSED DOWNDRAUGHT INITIAL
                                 ! MASS FLUX (PA/S)
!
      LOGICAL BDDI(NPNTS)        ! MASK FOR POINTS WHERE DOWNDRAUGHT
                                 ! MIGHT OCCUR
!
      LOGICAL BDDI_C(NTERM)      ! COMPRESSED MASK FOR POINTS WHERE
                                 ! DOWNDRAUGHT MAY INITIATE
!
      REAL QE_K_C(NTERM)         ! COMPRESSED ENVIRONMENT MIXING
                                 ! RATIO OF LAYER K (KG/KG)
!
      REAL QE_KM1_C(NTERM)       ! COMPRESSED ENVIRONMENT MIXING
                                 ! RATIO OF LAYER K-1 (KG/KG)
!
      REAL THE_K_C(NTERM)        ! COMPRESSED POTENTIAL TEMPERATURE
                                 ! OF ENVIRONMENT IN LAYER K (K)
!
      REAL THE_KM1_C(NTERM)      ! COMPRESSED POTENTIAL TEMPERATURE
                                 ! OF ENVIRONMENT IN LAYER K-1 (K)
      REAL TRAE_K_C(NTERM,NTRA)  ! COMPRESSED TRACER OF ENVIRONMENT
                                 ! IN LAYER K (KG/KG)
!
      REAL TRAE_KM1_C(NTERM,NTRA)! COMPRESSED TRACER OF ENVIRONMENT
                                 ! IN LAYER K-1 (KG/KG)
!
      REAL RAIN_C(NTERM)         ! COMPRESSED SURFACE RAINFALL
                                 ! (KG/M**2/S)
!
      REAL SNOW_C(NTERM)         ! COMPRESSED SURFACE SNOWFALL
                                 ! (KG/M**2/S)
!
      REAL FLX_UD_K_C(NTERM)     ! UPDRAUGHT MASS FLUX AT LAYER K
!
      REAL RAIN_ENV(NTERM)       ! AMOUNT OF RAINFALL PASSING THROUGH
                                 ! ENVIRONMENT (KG/M**2/S)
!
      REAL SNOW_ENV(NTERM)       ! AMOUNT OF SNOWFALL PASSING THROUGH
                                 ! ENVIRONMENT (KG/M**2/S)
!
      REAL RAIN_DD(NTERM)        ! AMOUNT OF RAINFALL PASSING THROUGH
                                 ! DOWNDRAUGHT (KG/M**2/S)
!
      REAL SNOW_DD(NTERM)        ! AMOUNT OF SNOWFALL PASSING THROUGH
                                 ! DOWNDRAUGHT (KG/M**2/S)
!
      LOGICAL BDD_START(NPNTS)   ! MASK FOR THOSE POINT WHERE
                                 ! DOWNDRAUGHT IS ABLE TO START
                                 ! FROM LEVEL K
!
      LOGICAL BDD_START_C(NTERM) ! COMPRESSED MASK FOR THOSE POINT
                                 ! WHERE DOWNDRAUGHT IS ABLE TO START
                                 ! FROM LEVEL K
!
      LOGICAL BDDWT_K(NPNTS)     ! MASK FOR POINTS IN DOWNDRAUGHT
                                 ! WHERE PPT IN LAYER K IS LIQUID
!
      LOGICAL BDDWT_K_C(NTERM)   ! COMPRESSED MASK FOR POINTS IN DD
                                 ! WHERE PPT IN LAYER K IS LIQUID
!
      LOGICAL BDDWT_KM1(NPNTS)   ! MASK FOR POINTS IN DOWNDRAUGHT
                                 ! WHERE PPT IN LAYER K-1 IS LIQUID
!
      LOGICAL BDDWT_KM1_C(NTERM) ! COMPRESSED MASK FOR POINTS IN DD
                                 ! WHERE PPT IN LAYER K-1 IS LIQUID
!
      LOGICAL BDD_ON(NPNTS)      ! MASK FOR THOSE POINTS WHERE DD
                                 ! CONTINUES FROM LAYER K+1
!
      LOGICAL BDD_ON_C(NTERM)    ! COMPRESSED MASK FOR POINTS WHERE DD
                                 ! CONTINUES FROM LAYER K+1
!
      LOGICAL L_SHALLOW_C(NTERM) ! COMPRESSED FLAG FOR SHALLOW CU
!
      INTEGER KMIN(NTERM)        ! FREEZING LEVEL WHERE ENTRAINMENT
                                 ! RATES ARE INCREASED
!
      REAL FLX_STRT(NPNTS)       ! MASSFLUX AT LEVEL WHERE DOWNDRAUGHT
                                 ! STARTS (PA/S)
!
      REAL FLX_STRT_C(NTERM)     ! COMPRESSED VALUE OF FLX_STRT
!
      REAL CCA_C(NTERM)          ! COMPRESSED CONVECTIVE CLOUD AMOUNT
!
      INTEGER INDEX2(NTERM)      ! INDEX OF WHERE ACTICE DOWNDRAUGHT
                                 ! OCCURS
!
      REAL LR_UD_REF(NTERM)      ! PRECIPITATION MIXING RATIO AT LOWEST
                                 ! PRECIPITATING LEVEL OF UD
!
!
      REAL P_CLD_TOP             ! PRESSURE AT CLOUD TOP (PA)
!
      REAL P_CLD_BASE            ! PRESSURE AT CLOUD BASE (PA)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!-----------------------------------------------------------------------
!
      EXTERNAL FLX_INIT, LAYER_DD, DD_INIT, DOWND
!
!-----------------------------------------------------------------------
! CALCULATE INDEX FOR COMPRESS ON BASIS OF BTERM
!-----------------------------------------------------------------------
!
      NDD=NTERM

      DO K=0,KCT+1
        DO I=1,NDD
          P_LAYER_CENTRES_C(I,K) = P_LAYER_CENTRES(INDEX1(I),K)
          EXNER_LAYER_CENTRES_C(I,K) =                                  &
     &                           EXNER_LAYER_CENTRES(INDEX1(I),K)
          P_LAYER_BOUNDARIES_C(I,K) =                                   &
     &                           P_LAYER_BOUNDARIES(INDEX1(I),K)
        ENDDO
      ENDDO
!
!----------------------------------------------------------------------
! INITIALISE LOGICAL ARRAYS AS FALSE
!-----------------------------------------------------------------------
!
      DO I=1,NPNTS
       BDDI(I) = .FALSE.
       BDD_START(I) = .FALSE.
      ENDDO
      DO I=1,NPNTS
       BDDWT_K(I) = .FALSE.
       BDDWT_KM1(I) = .FALSE.
      ENDDO
      DO I=1,NPNTS
       BDD_ON(I) = .FALSE.
      ENDDO
      DO I2=1,NDD
        BDDI(INDEX1(I2)) = .TRUE.
      END DO
!
!----------------------------------------------------------------------
! CALCULATE INITIAL DOWNDRAUGHT MASS FLUX
!-----------------------------------------------------------------------
!
! DEPENDS ON: flx_init
        CALL FLX_INIT (NPNTS, KCT, ICCB, ICCT, FLX, FLX_DD_K, BDDI,     &
     &   FLX_STRT)
!
!-----------------------------------------------------------------------
! COMPRESS ALL INPUT ARRAYS FOR THE DOWNDRAUGHT CALCULATION
!-----------------------------------------------------------------------
!
      DO 10 K = KCT+1,2,-1
!
         DO I=1,NDD
            TH_K_C(I) = THP(INDEX1(I),K)
            Q_K_C(I) = QP(INDEX1(I),K)
         ENDDO
         DO I=1,NDD
            THE_K_C(I) = THE(INDEX1(I),K)
            THE_KM1_C(I) = THE(INDEX1(I),K-1)
         ENDDO
         DO I=1,NDD
            QE_K_C(I) = QE(INDEX1(I),K)
            QE_KM1_C(I) = QE(INDEX1(I),K-1)
         ENDDO
         DO I=1,NDD
            DTHBYDT_K_C(I) = DTHBYDT(INDEX1(I),K)
            DTHBYDT_KM1_C(I) = DTHBYDT(INDEX1(I),K-1)
         ENDDO
         DO I=1,NDD
            DQBYDT_K_C(I) = DQBYDT(INDEX1(I),K)
            DQBYDT_KM1_C(I) = DQBYDT(INDEX1(I),K-1)
            EXNER_KM12_C(I) = EXNER_layer_boundaries(INDEX1(I),K-1)
            EXNER_KP12_C(I) = EXNER_layer_boundaries(INDEX1(I),K)
         ENDDO
         DO I=1,NDD
            EXNER_KM32_C(I) = EXNER_layer_boundaries(INDEX1(I),K-2)
         ENDDO
         DO I=1,NDD
            PRECIP_K_C(I) = PRECIP(INDEX1(I),K)
         ENDDO
         DO I=1,NDD
            FLX_UD_K_C(I) = FLX(INDEX1(I),K)
            BWATER_K_C(I) = BWATER(INDEX1(I),K)
         END DO
!
         IF(L_TRACER)THEN
!
         DO KTRA=1,NTRA
           DO I=1,NDD
             TRA_K_C(I,KTRA) = TRAP(INDEX1(I),K,KTRA)
             TRAE_K_C(I,KTRA) = TRAE(INDEX1(I),K,KTRA)
         ENDDO
         DO I=1,NDD
             TRAE_KM1_C(I,KTRA) = TRAE(INDEX1(I),K-1,KTRA)
             DTRA_K_C(I,KTRA) = DTRABYDT(INDEX1(I),K,KTRA)
             DTRA_KM1_C(I,KTRA) = DTRABYDT(INDEX1(I),K-1,KTRA)
           END DO
         END DO
!
         END IF
!
         IF (K == KCT+1) THEN
          DO I=1,NDD
            FLX_DD_K_C(I) = FLX_DD_K(INDEX1(I))
            FLX_STRT_C(I) = FLX_STRT(INDEX1(I))
         ENDDO
         DO I=1,NDD
            PSTAR_C(I) = PSTAR(INDEX1(I))
            recip_pstar_c(I)=recip_pstar(index1(I))
         ENDDO
         DO I=1,NDD
            ICCB_C(I) = ICCB(INDEX1(I))
            BDDI_C(I) = BDDI(INDEX1(I))
         ENDDO
         DO I=1,NDD
            BDD_START_C(I) = BDD_START(INDEX1(I))
            RAIN_C(I) = RAIN(INDEX1(I))
         ENDDO
         DO I=1,NDD
            SNOW_C(I) = SNOW(INDEX1(I))
            BDDWT_K_C(I) = BDDWT_K(INDEX1(I))
         ENDDO
         DO I=1,NDD
            BDDWT_KM1_C(I) = BDDWT_KM1(INDEX1(I))
            BDD_ON_C(I) = BDD_ON(INDEX1(I))
         ENDDO
         DO I=1,NDD
            CCA_C(I) = CCA(INDEX1(I))
            LR_UD_REF(I) = 0.0
          END DO
         END IF
!
!----------------------------------------------------------------------
! IF BELOW CONVECTIVE CLOUD BASE DOWNDRAUGHT NOT ALLOWED TO FORM
!----------------------------------------------------------------------
!
      DO I=1,NDD
       IF (K <  ICCB_C(I)) BDDI_C(I)=.FALSE.
      END DO
!
!-----------------------------------------------------------------------
! RESET EN/DETRAINMENT RATES FOR DOWNDRAUGHT
!-----------------------------------------------------------------------
!
! DEPENDS ON: layer_dd
      CALL LAYER_DD (NDD,K,KCT,THE_K_C,THE_KM1_C,FLX_STRT_C,            &
     &               P_LAYER_CENTRES_C,                                 &
     &               P_LAYER_BOUNDARIES_C,                              &
     &               EXNER_LAYER_CENTRES_C,                             &
     &               EXNER_KM12_C,EXNER_KP12_C,                         &
     &               EXNER_KM32_C,PSTAR_C,PK,P_KM1,DELPK,DELPKM1,EXK,   &
     &               EXKM1,AMDETK,EKM14,EKM34,KMIN,BDDI_C,              &
     &               recip_pstar_c)
!----------------------------------------------------------------------
! IF LEVEL K WITHIN 150MB OF SURFACE THEN DOWNDRAUGHT NOT ALLOWED TO
! FORM
!----------------------------------------------------------------------
!
      DO I=1,NDD
       IF (PK(I) >  (PSTAR_C(I)-15000.0)) BDDI_C(I)=.FALSE.
      END DO
!
!
!-----------------------------------------------------------------------
! INITIALISE DOWNDRAUGHT
! DOWNDRAUGHT NOT ALLOWED TO FORM FROM CLOUD TOP LAYER (KCT+1)
! OR FROM BELOW CLOUD BASE
!-----------------------------------------------------------------------
!
! DEPENDS ON: dd_init
       CALL DD_INIT(NDD,NTERM,TH_K_C,Q_K_C,THE_K_C,QE_K_C,PK,EXK,       &
     &              THDD_K,QDD_K,DELTD,DELQD,BDD_START_C,K,BDDI_C,      &
     &              BDD_ON_C,L_TRACER,NTRA,TRA_K_C,                     &
     &              TRAE_K_C,TRADD_K,DELTRAD)
!
!-----------------------------------------------------------------------
! UPDATE MASK FOR WHERE DOWNDRAUGHT OCCURS
!-----------------------------------------------------------------------
!
      DO I=1,NDD
        IF (BDD_START_C(I).OR.BDD_ON_C(I)) BDD_ON_C(I)=.TRUE.
      END DO
!
      IF(FLG_DWN_FLX) THEN
!
! IF DOWNDRAUGHT INITIATED SET DIAGNOSTIC ARRAY
!
       DO I=1,NDD
       IF(BDD_START_C(I)) DD_FLUX(INDEX1(I),K)=FLX_DD_K(INDEX1(I))
       END DO
      END IF
!
      NDDON_TMP = 0
      DO I=1,NDD
        IF (BDD_ON_C(I)) THEN
          NDDON_TMP = NDDON_TMP+1
        END IF
      END DO
!
!-----------------------------------------------------------------------
! CALL DOWNDRAUGHT ROUTINE
!-----------------------------------------------------------------------
!

! DEPENDS ON: downd
      CALL DOWND(NDD,NTERM,K,KCT,THDD_K,QDD_K,THE_K_C,THE_KM1_C,        &
     &           QE_K_C,QE_KM1_C,DTHBYDT_K_C,DTHBYDT_KM1_C,DQBYDT_K_C,  &
     &           DQBYDT_KM1_C,FLX_DD_K_C,P_KM1,DELPK,DELPKM1,EXK,       &
     &           EXKM1,DELTD,DELQD,AMDETK,EKM14,EKM34,PRECIP_K_C,       &
     &           RAIN_C,SNOW_C,ICCB_C,BWATER_K_C,BDD_START_C,           &
     &           BDDWT_K_C,BDDWT_KM1_C,BDD_ON_C,RAIN_ENV,SNOW_ENV,      &
     &           RAIN_DD,SNOW_DD,FLX_UD_K_C,TIMESTEP,CCA_C,NDDON_TMP,   &
     &           LR_UD_REF,L_TRACER,NTRA,TRADD_K,                       &
     &           TRAE_K_C,TRAE_KM1_C,DTRA_K_C,DTRA_KM1_C,DELTRAD        &
     &           )
!
!-----------------------------------------------------------------------
! DECOMPRESS/EXPAND THOSE VARIABLES WHICH ARE TO BE OUTPUT
!-----------------------------------------------------------------------
!
        DO I=1,NDD
         DTHBYDT(INDEX1(I),K) = DTHBYDT_K_C(I)
         DTHBYDT(INDEX1(I),K-1) = DTHBYDT_KM1_C(I)
         ENDDO
         DO I=1,NDD
         DQBYDT(INDEX1(I),K) = DQBYDT_K_C(I)
         DQBYDT(INDEX1(I),K-1) = DQBYDT_KM1_C(I)
         ENDDO
         IF(FLG_DWN_FLX) THEN
          DO I=1,NDD
           IF(BDD_ON_C(I)) THEN
!
! NEED TO CHECK THAT POINT WOULD BE SELECTED IN S.R DOWND OR ELSE
! NOT SENSIBLE TO SET ENTRAINMENT AND DETRAINMENT RATES IN DIAGNOSTICS
!
            DD_FLUX(INDEX1(I),K-1) = FLX_DD_K_C(I)
           ENDIF
          ENDDO
         ENDIF
         IF(FLG_ENTR_DWN) THEN
          DO I=1,NDD
           IF(BDD_ON_C(I)) THEN
            ENTRAIN_DWN(INDEX1(I),K)=(1.0-AMDETK(I))*                   &
     &                           (EKM14(I)+EKM34(I)*(1.0+EKM14(I)))*    &
     &                             DD_FLUX(INDEX1(I),K)
           ENDIF
          ENDDO
         ENDIF
         IF(FLG_DETR_DWN) THEN
          DO I=1,NDD
           IF(BDD_ON_C(I)) THEN
            DETRAIN_DWN(INDEX1(I),K)=-AMDETK(I)*DD_FLUX(INDEX1(I),K)
           ENDIF
          ENDDO
         ENDIF
         DO I=1,NDD
         IF (K == 2) THEN
          RAIN(INDEX1(I)) = RAIN_C(I)
          SNOW(INDEX1(I)) = SNOW_C(I)
         END IF
         PRECIP(INDEX1(I),K) = PRECIP_K_C(I)
        END DO

        Do I=1,NDD
          RAIN_3D(INDEX1(I), K-1) = RAIN_3D(INDEX1(I), K-1)                   &
                                  + RAIN_DD(I) + RAIN_ENV(I)

          SNOW_3D(INDEX1(I), K-1) = SNOW_3D(INDEX1(I), K-1)                   &
                                  + SNOW_DD(I) + SNOW_ENV(I)
        End Do
!
       IF(L_TRACER)THEN
!
        DO KTRA=1,NTRA
          DO I=1,NDD
            DTRABYDT(INDEX1(I),K,KTRA) = DTRA_K_C(I,KTRA)
            DTRABYDT(INDEX1(I),K-1,KTRA) = DTRA_KM1_C(I,KTRA)
          END DO
        END DO
!
       END IF
!
!----------------------------------------------------------------------
!   END OF MAIN K LOOP
!----------------------------------------------------------------------
!
 10   CONTINUE
!
      RETURN
      END SUBROUTINE DD_CALL
!
#endif
