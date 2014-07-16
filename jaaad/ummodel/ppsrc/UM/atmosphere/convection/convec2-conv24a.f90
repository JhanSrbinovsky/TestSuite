
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE CONVEC2------------------------------------------------
!LL
!LL  PURPOSE : COMPLETES LIFTING OF THE PARCEL FROM LAYER K TO K+1
!LL
!LL            CALL SUBROUTINE PARCEL AND ENVIRON
!LL
!LL            SUBROUTINE PARCEL CALCULATES AN INITIAL MASS FLUX,
!LL            CARRIES OUT THE DETRAINMENT CALCULATION, TESTS
!LL            TO SEE IF CONVECTION IS TERMINATING AND CALCULATES THE
!LL            PRECIPITATION RATE FROM LAYER K+1
!LL
!LL            SUBROUTINE ENVIRON CALCULATES THE EFFECT OF CONVECTION
!LL            UPON THE LARGE-SCALE ATMOSPHERE
!LL
!LL  SUITABLE FOR SINGLE COLUMN MODEL USE
!LL
!LL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!LL
!LL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
!LL
!LLEND-----------------------------------------------------------------
!
!*L  ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE CONVEC2 (NPNTS,NP_FULL,NLEV,K,THEK,THEKP1,QEK,QEKP1,   &
     &                   QCLEK,QCLEKP1,QCFEK,QCFEKP1,                   &
     &                   CFLEK,CFLEKP1,CFFEK,CFFEKP1,BCFEK,             &
     &                   BCFEKP1,QSEKP1,DQSKP1,PSTAR,                   &
     &                   THPK,QPK,QCLPK,QCFPK,THPKP1,                   &
     &                   QPKP1,QCLPKP1,QCFPKP1,START_LEV,               &
     &         XSQKP1,RBUOY,QSEK,DQSK,THPI,QPI,EXPI,XPK,FLXK,BWK,BWKP1, &
     &                   BGMKP1,BGMK,BLOWST,BLAND,BTERM,DEPTH,PREKP1,   &
     &                   DTHEK,DQEK,DQCLEK,DQCFEK,DCFLEK,               &
     &                   DCFFEK,DBCFEK,DTHEKP1,DQEKP1,                  &
     &                   DQCLEKP1,DQCFEKP1,DCFLEKP1,DCFFEKP1,           &
     &                   DBCFEKP1,BINIT,CCA,ICCB,ICCT,                  &
     &                   TCW,EKP14,EKP34,AMDETK,PK,PKP1,                &
     &                   EXK,EXKP1,DELEXKP1,DELPK,DELPKP1,              &
     &                   CCLWP,CCW,LCCA,LCBASE,LCTOP,LCCLWP,T1_SD,      &
     &                   Q1_SD,L_MOM_GK,UEK,UEKP1,VEK,VEKP1,UPK,VPK,    &
     &                   UPKP1,VPKP1,DUEK,DUEKP1,DVEK,DVEKP1,           &
     &                   EFLUX_U_UD,EFLUX_V_UD,                         &
     &                   delp_uv_k, delp_uv_kp1,                        &
     &                   THPIXS_V,QPIXS_V,XSBMIN_V,                     &
     &                   L_SHALLOW,L_MID,L_TRACER,NTRA,TRAEK,TRAEKP1,   &
     &                   TRAPK,TRAPKP1,DTRAEK,DTRAEKP1,CAPE,DCPBYDT,    &
     &                   MAX_CFL_C, TIMESTEP,RBUOY_P_HERE,THE_HERE,     &
     &                   THP_HERE,QE_HERE,QP_HERE,RBUOY_P_OLD,          &
     &                   AD_ON,SDET_ON,NEW_TERMC,DELTAK,L_CALC_DXEK,    &
     &                   L_Q_INTERACT,FLXKP12,CUMULUS,RELH,DPTOT)
!
      Use cv_run_mod, Only:                                             &
          cape_opt

      IMPLICIT NONE
!
!----------------------------------------------------------------------
!  MODEL CONSTANTS
!----------------------------------------------------------------------
!
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
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
!
!----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!----------------------------------------------------------------------
!
      INTEGER NPNTS          ! IN VECTOR LENGTH
!
      INTEGER NP_FULL        ! IN FULL VECTOR LENGTH
!
      INTEGER NLEV           ! IN NUMBER OF MODEL LAYERS
!
      INTEGER NTRA           ! IN NUMBER OF TRACER VARIABLES
!
      INTEGER I,KTRA         ! LOOP COUNTER
!
      INTEGER K              ! PRESENT MODEL LAYER
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT
!----------------------------------------------------------------------
!
      REAL THEK(NPNTS)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K (K)
!
      REAL THEKP1(NPNTS)     ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
!
      REAL QEK(NPNTS)        ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL QEKP1(NPNTS)      ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL QCLEK(NPNTS)      ! IN LIQUID CONDENSATE MIXING RATIO OF
!                                 CLOUD ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL QCLEKP1(NPNTS)    ! IN LIQUID CONDENSATE MIXING RATIO OF
!                                 CLOUD ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL QCFEK(NPNTS)      ! IN FROZEN CONDENSATE MIXING RATIO OF
!                                 CLOUD ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL QCFEKP1(NPNTS)    ! IN FROZEN CONDENSATE MIXING RATIO OF
!                                 CLOUD ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL BCFEK(NPNTS)      ! IN TOTAL  CLOUD VOLUME FRACTION OF
!                                 CLOUD ENVIRONMENT IN LAYER K ( )
!
      REAL BCFEKP1(NPNTS)    ! IN TOTAL  CLOUD VOLUME FRACTION OF
!                                 CLOUD ENVIRONMENT IN LAYER K+1 ( )
!
      REAL CFLEK(NPNTS)      ! IN LIQUID CLOUD VOLUME FRACTION OF
!                                 CLOUD ENVIRONMENT IN LAYER K ( )
!
      REAL CFLEKP1(NPNTS)    ! IN LIQUID CLOUD VOLUME FRACTION OF
!                                 CLOUD ENVIRONMENT IN LAYER K+1 ( )
!
      REAL CFFEK(NPNTS)      ! IN FROZEN CLOUD VOLUME FRACTION OF
!                                 CLOUD ENVIRONMENT IN LAYER K ( )
!
      REAL CFFEKP1(NPNTS)    ! IN FROZEN CLOUD VOLUME FRACTION OF
!                                 CLOUD ENVIRONMENT IN LAYER K+1 ( )
!
      REAL QCLPK(NPNTS)      ! IN PARCEL LIQUID CONDENSATE MIXING RATIO
!                                 IN LAYER K (KG/KG)
!
      REAL QCFPK(NPNTS)      ! IN PARCEL FROZEN CONDENSATE MIXING RATIO
!                                 IN LAYER K (KG/KG)
!
      REAL UEK(NPNTS)        ! IN U IN ENVIRONMENT IN LAYER K (M/S)
!
      REAL UEKP1(NPNTS)      ! IN U IN ENVIRONMENT IN LAYER K+1 (M/S)
!
      REAL VEK(NPNTS)        ! IN V IN ENVIRONMENT IN LAYER K (M/S)
!
      REAL VEKP1(NPNTS)      ! IN V IN ENVIRONMENT IN LAYER K+1 (M/S)
!
      REAL TRAEK(NP_FULL,                                               &
                             ! IN TRACER CONTENT OF CLOUD
     &           NTRA)       !    ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL TRAEKP1(NP_FULL,                                             &
                             ! IN TRACER CONTENT OF CLOUD
     &             NTRA)     !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL QSEKP1(NPNTS)     ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
!
      REAL DQSKP1(NPNTS)     ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT IN LAYER K+1
                             !    (KG/KG)
!
      REAL PSTAR(NPNTS)      ! IN SURFACE PRESSURE (PA)
!
      REAL THPKP1(NPNTS)     ! IN PARCEL POTENTIAL TEMPERATURE
                             !    IN LAYER K+1 (K)
!
      REAL QPKP1(NPNTS)      ! IN PARCEL MIXING RATIO IN LAYER K+1
                             !    (KG/KG)
!
      REAL UPKP1(NPNTS)      ! IN PARCEL U IN LAYER K+1 (M/S)
!
      REAL VPKP1(NPNTS)      ! IN PARCEL V IN LAYER K+1 (M/S)
!
      REAL TRAPKP1(NP_FULL,                                             &
                             ! IN PARCEL TRACER CONTENT IN LAYER
     &             NTRA)     !    K+1 (KG/KG)
!
      REAL XSQKP1(NPNTS)     ! IN EXCESS WATER IN PARCEL AFTER LIFTING
                             !    LAYER K TO K+1 (KG/KG)
!
      REAL RBUOY(NPNTS)      ! IN PARCEL BUOYANCY IN LAYER K+1 (K)
!
      REAL QSEK(NPNTS)       ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
!
      REAL DQSK(NPNTS)       ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT OF LAYER K
                             !    (KG/KG/K)
!
      REAL THPI(NPNTS)       ! IN INITIAL PARCEL POTENTIAL TEMPERATURE
                             !    (K)
!
      REAL QPI(NPNTS)        ! IN INITIAL PARCEL MIXING RATIO
                             !    (KG/KG)
!
      REAL EXPI(NPNTS)       ! IN INITIAL PARCEL EXNER PRESSURE
!
      REAL TIMESTEP          ! IN CONVECTION TIMESTEP (SECS)
!                                (= model timestep if conv. called once)
!
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
      LOGICAL BWK(NPNTS)     ! IN  And layer K
!
      LOGICAL BGMKP1(NPNTS)  ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K+1
!
      LOGICAL BLAND(NPNTS)   ! IN LAND/SEA MASK
!
      LOGICAL BINIT(NPNTS)   ! IN MASK FOR THOSE POINTS AT WHICH
                             !    CONVECTION IS OCCURING
!
      LOGICAL BLOWST(NPNTS)  ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
!
      LOGICAL L_SHALLOW(NPNTS),                                         &
                                ! IN SWITCHES FOR TYPE OF CONVECTION
     &        L_MID(NPNTS)      !    LIKELY TO DEVELOP
!
      LOGICAL CUMULUS(NPNTS) ! IN CUMULUS CONVECTION (DEEP OR SHALLOW)
!
      LOGICAL L_TRACER       ! IN LOGICAL SWITCH FOR INCLUSION OF
                             !    TRACERS
!
      LOGICAL L_MOM_GK       ! IN LOGICAL SWITCH FOR INCLUSION OF
                             !    MOMENTUM TRANSPORT (Gregory-Kershaw)
!
      LOGICAL L_CALC_DXEK    ! IN Switch for calculating condensate inc.
!
      LOGICAL L_Q_INTERACT   ! IN Switch allows overwriting of parcel
!                                 variables when calculating condensate
!                                 increments (will alter results).
!
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!
      REAL EKP34(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!
      REAL AMDETK(NPNTS)     ! IN MIXING DETRAINMENT COEFFICIENT
                             !    AT LEVEL K MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!
      REAL DELPKP12(NPNTS)   ! IN PRESSURE DIFFERENCE BETWEEN
                             !    MID-POINTS OF LAYERS K AND K+1
                             !    (PA)
!
      REAL PK(NPNTS)         ! IN PRESSURE AT MID-POINT OF LAYER K
                             !    (PA)
!
      REAL PKP1(NPNTS)       ! IN PRESSURE AT MID-POINT OF LAYER K+1
                             !    (PA)
!
      REAL EXK(NPNTS)        ! IN EXNER RATIO AT MID-POINT OF LAYER K
!
      REAL EXKP1(NPNTS)      ! IN EXNER RATIO AT MID-POINT OF LAYER K+1
!
      REAL DELEXKP1(NPNTS)   ! IN DIFFERENCE IN EXNER RATIO BETWEEN
                             !    MID-POINTS OF LAYERS K AND K+1
!
      REAL DELPK(NPNTS)      ! IN DIFFERENCE IN PRESSURE ACROSS LAYER K
                             !    (PA)
!
      REAL DELPKP1(NPNTS)    ! IN DIFFERENCE IN PRESSURE ACROSS
                             !    LAYER K+1 (PA)
!
      REAL DELP_UV_K(NPNTS)   ! IN PRESSURE DIFFERENCE ACROSS UV LAYER K
                              !    (PA)
!
      REAL DELP_UV_KP1(NPNTS) ! IN PRESSURE DIFFERENCE ACROSS UV LAYER
                              !    K+1 (PA)
!
      REAL T1_SD(NPNTS)      ! IN Standard deviation of turbulent
                             !    fluctuations of layer 1
                             !    temperature (K).
      REAL Q1_SD(NPNTS)      ! IN Standard deviation of turbulent
                             !    fluctuations of layer 1
                             !    humidity (kg/kg).
      REAL THPIXS_V(NPNTS)   ! IN  PARCEL EXCESSES OF THETA
      REAL QPIXS_V(NPNTS)    !     AND Q
      REAL XSBMIN_V(NPNTS)   ! IN  THRESHOLD BUOYANCY FOR
                             !     FORCED DETRAINMENT

      INTEGER START_LEV(NPNTS)  ! IN LEVEL AT WHICH CONVECTION INITIATED
!
      INTEGER, INTENT(IN) :: AD_ON      ! Flag for adaptive detrainment
!
      INTEGER, INTENT(IN) :: SDET_ON    ! Flag smoothed forced detrainment
!
      INTEGER, INTENT(IN) :: NEW_TERMC  ! Flag for simplified
                                        ! termination of convection

      REAL, INTENT(IN) :: RBUOY_P_OLD(NPNTS)
                                        ! BUOYANCY FROM PREVIOUS LEVEL

!----------------------------------------------------------------------
! VARIABLES WHICH ARE INPUT BUT WHICH ARE ALSO UPDATED IN THIS ROUTINE
!----------------------------------------------------------------------
!
      REAL THPK(NPNTS)       ! INOUT
                             ! IN  PARCEL POTENTIAL TEMPERATURE
                             !     IN LAYER K (K)
                             ! OUT PARCEL POTENTIAL TEMPERATURE
                             !     IN LAYER K+1 (K)
!
      REAL QPK(NPNTS)        ! INOUT
                             ! IN  PARCEL MIXING RATIO IN LAYER K
                             !     (KG/KG)
                             ! OUT PARCEL MIXING RATIO IN LAYER K+1
                             !     (KG/KG)
      REAL QCLPKP1(NPNTS)    ! INOUT
!                              IN  PARCEL LIQUID CONDENSATE MIXING RATIO
!                                  IN LAYER K+1 AFTER DRY ASCENT ONLY
!                                  (KG/KG)
!                              OUT PARCEL LIQUID CONDENSATE IN LAYER K+1
!                                  (KG/KG)
!
      REAL QCFPKP1(NPNTS)    ! INOUT
!                              IN  PARCEL FROZEN CONDENSATE MIXING RATIO
!                                  IN LAYER K+1 AFTER DRY ASCENT ONLY
!                                  (KG/KG)
!                              OUT PARCEL FROZEN CONDENSATE IN LAYER K+1
!                                  (KG/KG)
!
      REAL UPK(NPNTS)        ! INOUT
                             ! IN  PARCEL U IN LAYER K (M/S)
                             ! OUT PARCEL U IN LAYER K+1 (M/S)
!
      REAL VPK(NPNTS)        ! INOUT
                             ! IN  PARCEL V IN LAYER K (M/S)
                             ! OUT PARCEL V IN LAYER K+1 (M/S)
!
      REAL TRAPK(NP_FULL,                                               &
                             ! INOUT
     &           NTRA)       ! IN  PARCEL TRACER CONTENT IN LAYER K
                             !     (KG/KG)
                             ! OUT PARCEL TRACER CONTENT IN LAYER K+1
                             !     (KG/KG)
!
      REAL XPK(NPNTS)        ! INOUT
                             ! IN  PARCEL CLOUD WATER IN LAYER K
                             !     (KG/KG)
                             ! OUT PARCEL CLOUD WATER IN LAYER K+1
                             !     (KG/KG)
! NOTE: XPK is identical with CONDENSATE MIXING RATIO OF CLOUD PARCEL
!       XPK will be redundant when L_calc_dxek == .True.  because
!       it can be calculated directly from (Qclpk + Qcfpk).
!
      REAL FLXK(NPNTS)       ! INOUT
                             ! IN  PARCEL MASSFLUX IN LAYER K (PA/S)
                             ! OUT PARCEL MASSFLUX IN LAYER K+1 (PA/S)
!
      LOGICAL BGMK(NPNTS)    ! INOUT
                             ! IN  MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K
                             ! OUT MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1
!
      REAL DTHEK(NPNTS)      ! INOUT
                             ! IN  INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (MAY BE NONE ZERO
                             !     DUE TO A PREVIOUS SPLIT FINAL
                             !     DETRAINEMNT CALCULATION) (K/S)
                             ! OUT UPDATED INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (K/S)
!
      REAL DQEK(NPNTS)       ! INOUT
                             ! IN  INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K DUE TO CONVECTION
                             !     (MAY BE NONE ZERO DUE TO A
                             !     PREVIOUS SPLIT FINAL DETRAINEMNT
                             !     CALCULATION) (KG/KG/S)
                             ! OUT UPDATED INCREMENT TO MODEL MIXING
                             !     RATIO IN LAYER K DUE TO CONVECTION
                             !     (KG/KG/S)
!
      REAL DQCLEK(NPNTS)     ! INOUT
!                              IN  INCREMENT TO MODEL LIQUID CONDENSATE
!                                  MIXING RATIO IN LAYER K DUE TO
!                                  CONVECTION (KG/KG/S)
!                                  (MAY BE NON-ZERO DUE TO A PREVIOUS
!                                  SPLIT FINAL DETRAINMENT CALCULATION)
!                              OUT UPDATED INCREMENT TO MODEL LIQUID
!                                  CONDENSATE MIXING RATIO IN LAYER K
!                                  DUE TO CONVECTION (KG/KG/S)
!
      REAL DQCFEK(NPNTS)     ! INOUT
!                              IN  INCREMENT TO MODEL FROZEN CONDENSATE
!                                  MIXING RATIO IN LAYER K DUE TO
!                                  CONVECTION (KG/KG/S)
!                                  (MAY BE NON-ZERO DUE TO A PREVIOUS
!                                  SPLIT FINAL DETRAINMENT CALCULATION)
!                              OUT UPDATED INCREMENT TO MODEL FROZEN
!                                  CONDENSATE MIXING RATIO IN LAYER K
!                                  DUE TO CONVECTION (KG/KG/S)
!
      REAL DBCFEK(NPNTS)     ! INOUT
!                              IN  INCREMENT TO MODEL TOTAL  CLOUD
!                                  VOLUME FRACTION IN LAYER K DUE TO
!                                  CONVECTION ( /S)
!                                  (MAY BE NON-ZERO DUE TO A PREVIOUS
!                                  SPLIT FINAL DETRAINMENT CALCULATION)
!                              OUT UPDATED INCREMENT TO MODEL TOTAL
!                                  CLOUD VOLUME FRACTION IN LAYER K
!                                  DUE TO CONVECTION ( /S)
!
      REAL DCFLEK(NPNTS)     ! INOUT
!                              IN  INCREMENT TO MODEL LIQUID CLOUD
!                                  VOLUME FRACTION IN LAYER K DUE TO
!                                  CONVECTION ( /S)
!                                  (MAY BE NON-ZERO DUE TO A PREVIOUS
!                                  SPLIT FINAL DETRAINMENT CALCULATION)
!                              OUT UPDATED INCREMENT TO MODEL LIQUID
!                                  CLOUD VOLUME FRACTION IN LAYER K
!                                  DUE TO CONVECTION ( /S)
!
      REAL DCFFEK(NPNTS)     ! INOUT
!                              IN  INCREMENT TO MODEL FROZEN CLOUD
!                                  VOLUME FRACTION IN LAYER K DUE TO
!                                  CONVECTION ( /S)
!                                  (MAY BE NON-ZERO DUE TO A PREVIOUS
!                                  SPLIT FINAL DETRAINMENT CALCULATION)
!                              OUT UPDATED INCREMENT TO MODEL FROZEN
!                                  CLOUD VOLUME FRACTION IN LAYER K
!                                  DUE TO CONVECTION ( /S)
!
      REAL DUEK(NPNTS)       ! INOUT
                             ! IN  INCREMENT TO MODEL U IN LAYER K
                             !     DUE TO CONVECTION (M/S**2)
                             ! OUT UPDATED INCREMENT TO U IN LAYER K
                             !     DUE TO CONVECTION (M/S**2)
!
      REAL DVEK(NPNTS)       ! INOUT
                             ! IN  INCREMENT TO MODEL V IN LAYER K
                             !     DUE TO CONVECTION (M/S**2)
                             ! OUT UPDATED INCREMENT TO V IN LAYER K
                             !     DUE TO CONVECTION (M/S**2)
!
      REAL DTRAEK(NP_FULL,                                              &
                             ! INOUT
     &            NTRA)      ! IN  INCREMENT TO MODEL TRACER IN LAYER
                             !     K DUE TO CONVECTION (MAY BE NON
                             !     ZERO DUE TO A PREVIOUS SPLIT
                             !     FINAL DETRAINMENT CALCULATION)
                             !     (KG/KG/S)
                             ! OUT UPDATED INCREMENT TO MODEL TRACER
                             !     IN LAYER K DUE TO CONVECTION
                             !     (KG/KG/S)
!
      REAL TCW(NPNTS)        ! INOUT
                             ! IN  TOTAL CONDENSED WATER SUMMED TO
                             !     LAYER K (KG/M**2/S)
                             ! OUT UPDATED TOTAL CONDENSED WATER
                             !     SUMMED TO LAYER K+1 (KG/M**2/S)
!
      REAL DEPTH(NPNTS)      ! INOUT
                             ! IN  DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K (M)
                             ! OUT UPDATED DEPTH OF CONVECTIVE CLOUD
                             !     TO LAYER K+1 (M)
!
      REAL CCLWP(NPNTS)      ! INOUT
                             ! IN  CONDENSED WATER PATH SUMMED TO
                             !     LAYER K (KG/M**2)
                             ! OUT UPDATED CONDENSED WATER PATH
                             !     SUMMED TO LAYER K+1 (KG/M**2)
!
      REAL CAPE(NPNTS)       ! IN  CONVECTIVE AVAILABLE POTENTIAL ENERGY
                             !     UP TO THE CURRENT CONVECTING
                             !     LAYER
                             ! OUT CONVECTIVE AVAILABLE POTENTIAL ENERGY
                             !     INCLUDING ADDITION DUE TO THE CAPE
                             !     WITHIN THE CURRENT LAYER
!
      REAL DCPBYDT(NPNTS)    ! IN  RATE OF CHANGE OF CAPE
                             ! OUT RATE OF CHANGE OF CAPE INCLUDING
                             !     CONTRIBUTION FROM CURRENT LAYER
!
      REAL EFLUX_U_UD(NPNTS),                                           &
                             ! IN  EDDY FLUX OF MOMENTUM DUE TO UD AT
     &     EFLUX_V_UD(NPNTS) !     BOTTOM OF LAYER
                             ! OUT EDDY FLUX OF MOMENTUM DUE TO UD AT
                             !     TOP OF LAYER
!
!
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE OUTPUT
!----------------------------------------------------------------------
!
      LOGICAL BTERM(NPNTS)   ! OUT MASK FOR PARCELS WHICH TERMINATE IN
                             !     LAYER K+1
!
      REAL PREKP1(NPNTS)     ! OUT PRECIPITATION FROM PARCEL AS IT
                             !     RISES FROM LAYER K TO K+1 (KG/M**2/S)
!
      REAL DTHEKP1(NPNTS)    ! OUT INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K+1 DUE TO
                             !     CONVECTION (K/S)
!
      REAL DQEKP1(NPNTS)     ! OUT INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K+1 DUE TO CONVECTION
                             !     (KG/KG/S)
!
      REAL DQCLEKP1(NPNTS)   ! OUT INCREMENT TO MODEL LIQUID CONDENSATE
!                                  MIXING RATIO IN LAYER K+1 DUE TO
!                                  CONVECTION (KG/KG/S)
!
      REAL DQCFEKP1(NPNTS)   ! OUT INCREMENT TO MODEL FROZEN CONDENSATE
!                                  MIXING RATIO IN LAYER K+1 DUE TO
!                                  CONVECTION (KG/KG/S)
!
      REAL DBCFEKP1(NPNTS)   ! OUT INCREMENT TO MODEL TOTAL  CLOUD
!                                  VOLUME FRACTION IN LAYER K+1 DUE TO
!                                  CONVECTION ( /S)
!
      REAL DCFLEKP1(NPNTS)   ! OUT INCREMENT TO MODEL LIQUID CLOUD
!                                  VOLUME FRACTION IN LAYER K+1 DUE TO
!                                  CONVECTION ( /S)
!
      REAL DCFFEKP1(NPNTS)   ! OUT INCREMENT TO MODEL FROZEN CLOUD
!                                  VOLUME FRACTION IN LAYER K+1 DUE TO
!                                  CONVECTION ( /S)
!
      REAL DUEKP1(NPNTS)     ! OUT INCREMENT TO MODEL U IN LAYER K+1
                             !     DUE TO CONVECTION (M/S**2)
!
      REAL DVEKP1(NPNTS)     ! OUT INCREMENT TO MODEL V IN LAYER K+1
                             !     DUE TO CONVECTION (M/S**2)
!
      REAL DTRAEKP1(NP_FULL,                                            &
                             ! OUT INCREMENT TO MODEL TRACER IN
     &              NTRA)    !     LAYER K+1 DUE TO CONVECTION
                             !     (KG/KG/S)
!
      REAL CCA(NPNTS)        ! OUT CONVECTIVE CLOUD AMOUNT (%)
!
      INTEGER ICCB(NPNTS)    ! OUT CONVECTIVE CLOUD BASE LEVEL
!
      INTEGER ICCT(NPNTS)    ! OUT CONVECTIVE CLOUD TOP LEVEL
!
      REAL CCW(NPNTS)        ! OUT CONVECTIVE CLOUD WATER(G/KG) ON
                             ! MODEL LEVELS
!
      REAL LCCA(NPNTS)       ! OUT LOWEST CONV.CLOUD AMOUNT (%)
!
      INTEGER LCBASE(NPNTS)  ! OUT LOWEST CONV.CLOUD BASE LEVEL
!
      INTEGER LCTOP(NPNTS)   ! OUT LOWEST CONV.CLOUD TOP LEVEL
!
      REAL LCCLWP(NPNTS)     ! OUT LOWEST CONV.CLOUD LIQ.WATER PATH
!
      REAL DELTAK(NPNTS)     ! PARCEL FORCED DETRAINMENT RATE
                             ! IN LAYER K MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
!
      REAL FLXKP12(NPNTS)    ! HALF LEVEL MASS FLUX
!
      REAL MAX_CFL_C(NPNTS)  ! OUT CFL RATIO
!
      REAL RELH(NPNTS)       !OUT RH integral (average when
!                            !        convection terminates)
      REAL DPTOT(NPNTS)      !OUT Delta P integral
!
!
      REAL, INTENT(OUT) :: RBUOY_P_HERE(NPNTS) ! BUOYANCY FOR SCM DIAGS

      REAL, INTENT(OUT) :: THE_HERE(NPNTS)   ! TH ENVIRON FOR SCM DIAGS

      REAL, INTENT(OUT) :: THP_HERE(NPNTS)   ! TH PARCEL FOR SCM DIAGS

      REAL, INTENT(OUT) :: QE_HERE(NPNTS)    ! Q ENVIRON FOR SCM DIAGS

      REAL, INTENT(OUT) :: QP_HERE(NPNTS)    ! Q PARCEL FOR SCM DIAGS
!
!----------------------------------------------------------------------
! VARIABLES WHICH ARE DEFINED LOCALLY
!
      REAL THRK(NPNTS)       ! PARCEL DETRAINMENT POTENTIAL
                             ! TEMPERATURE IN LAYER K (K)
!
      REAL QRK(NPNTS)        ! PARCEL DETRAINMENT MIXING RATIO
                             ! IN LAYER K (KG/KG)
!
      REAL XPKP1(NPNTS)      ! PARCEL CLOUD WATER IN LAYER K+1 (KG/KG)
!
      REAL FLXKP1(NPNTS )    ! PARCEL MASSFLUX IN LAYER K+1 (PA/S)
!
      REAL THVP,THVE,RHO     ! VIRTUAL TEMPERATURE OF PARCEL, VIRTUAL
                             ! TEMPERATURE OF ENVIRONMENT AND
                             ! DENSITY REQUIRED IN CAPE CALCULATIONS
!
      REAL DQEK_NONPC2(NPNTS)   ! Increment to qek if PC2 was not
!                                 in place
      REAL DTHEK_NONPC2(NPNTS)  ! Increment to thek if PC2 was not
!                                 in place
      REAL DQEKP1_NONPC2(NPNTS) ! Increment to thek if PC2 was not
!                                 in place
      REAL DTHEKP1_NONPC2(NPNTS)! Increment to thekp1 if PC2 was not
!                                 in place
      REAL TMP_DCPBYDT          ! Temporary dcpbydt
!
!----------------------------------------------------------------------
! EXTERNAL ROUTINES CALLED
!----------------------------------------------------------------------
!
      EXTERNAL PARCEL,ENVIRON
!
!*---------------------------------------------------------------------
!L
!L----------------------------------------------------------------------
!L COMPLETE LIFTING PARCELS TO LAYER K+1
!L
!L SUBROUTINE PARCEL
!L
!L UM DOCUMENTATION PAPER 27
!L SECTIONS (5),(6),(7),(8),(9)
!L----------------------------------------------------------------------
!L
! DEPENDS ON: parcel
       CALL PARCEL (K,NPNTS,NLEV,PSTAR,THEKP1,THEK,QEKP1,QEK,           &
     &              QSEK,QSEKP1,DQSK,DQSKP1,BLAND,BWKP1,                &
     &              DELTAK,FLXK,THPK,QPK,THRK,QRK,                      &
     &              BTERM,THPKP1,QPKP1,PREKP1,XPK,XPKP1,FLXKP1,         &
     &              QCLPK,QCFPK,QCLPKP1,QCFPKP1,                        &
     &              XSQKP1,THPI,QPI,EXPI,BGMK,BGMKP1,BLOWST,RBUOY,      &
     &              XSBMIN_V,                                           &
     &              CCA,ICCB,ICCT,TCW,DEPTH,                            &
     &              EKP14,EKP34,AMDETK,DELPKP1,PK,PKP1,                 &
     &              EXK,EXKP1,DELEXKP1,CCLWP,CCW,                       &
     &              LCCA,LCBASE,LCTOP,LCCLWP,L_SHALLOW,                 &
     &              RBUOY_P_HERE,THE_HERE,                              &
     &              THP_HERE,QE_HERE, QP_HERE ,RBUOY_P_OLD, AD_ON,      &
     &              NEW_TERMC,L_Q_INTERACT,START_LEV,FLXKP12)
!L
!L----------------------------------------------------------------------
!L CALCULATE THE EFFECT ON THE ENVIRONMENT (EXCEPT FOR THE
!L THE EVAPORATION OF PRECIPITATION AND CHANGE OF PHASE)
!L
!L SUBROUTINE ENVIRON
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (10)
!L----------------------------------------------------------------------
!L
! DEPENDS ON: environ
       CALL ENVIRON (K,NPNTS,NP_FULL,DTHEK,DQEK,DTHEKP1,DQEKP1,         &
     &               DQCLEK, DQCFEK, DQCLEKP1, DQCFEKP1,                &
     &               DBCFEK,DCFLEK,DCFFEK, DBCFEKP1,DCFLEKP1,DCFFEKP1,  &
     &               DQEK_NONPC2, DTHEK_NONPC2,                         &
     &               DQEKP1_NONPC2, DTHEKP1_NONPC2,                     &
     &               THEK,QEK,QCLEK,QCFEK,BCFEK,CFLEK,CFFEK,            &
     &               DELTAK, FLXK, THPK, QPK, QCLPK, QCFPK,             &
     &               THRK,QRK,THEKP1,QEKP1,QCLEKP1,                     &
     &               QCFEKP1,BCFEKP1,CFLEKP1,CFFEKP1,QCLPKP1,QCFPKP1,   &
     &               BTERM,THPKP1,QPKP1,XPK,XPKP1,BWK,BWKP1,FLXKP1,     &
     &               BLOWST,EKP14,EXK,EXKP1,DELPK,DELPKP1,              &
     &               DELP_UV_K, DELP_UV_KP1,                            &
     &               THPIXS_V,QPIXS_V,AMDETK,T1_SD,Q1_SD,               &
     &               L_MOM_GK,DUEK,DVEK,DUEKP1,DVEKP1,                  &
     &               UEK,VEK,UPK,VPK,UEKP1,VEKP1,UPKP1,VPKP1,EFLUX_U_UD,&
     &               EFLUX_V_UD,L_SHALLOW,CUMULUS,                      &
     &               L_MID,L_TRACER,NTRA,DTRAEK,DTRAEKP1,TRAEK,         &
     &               TRAPK,TRAEKP1,TRAPKP1,MAX_CFL_C,                   &
     &               TIMESTEP,L_CALC_DXEK,L_Q_INTERACT,                 &
     &               SDET_ON)
!
      DO 10 I=1,NPNTS
!
!-----------------------------------------------------------------------
! RESET BINIT WHERE CONVECTION HAS TERMINATED
!-----------------------------------------------------------------------
!
        BINIT(I) = .NOT.BTERM(I)
   10 CONTINUE
!
!L---------------------------------------------------------------------
!L CALCULATE CONTRIBUTION TO CAPE AND RATE OF CHANGE OF CAPE DUE TO
!L THE UPDRAUGHT
!L---------------------------------------------------------------------
!
      IF (CAPE_OPT  /=  0) THEN
       DO I=1,NPNTS
         THVP=THPK(I)*(1.0+C_VIRTUAL*QPK(I))
         THVE=THEK(I)*(1.0+C_VIRTUAL*QEK(I))
         RHO=PK(I)/(R*THEK(I)*EXK(I))
!
         CAPE(I)=CAPE(I)+(THVP-THVE)*DELPK(I)/(RHO*THVE)
         RELH(I)=RELH(I)+(QEK(I)/QSEK(I))*DELPK(I)
         DPTOT(I)=DPTOT(I)+DELPK(I)
!
         IF (L_Q_INTERACT) THEN
!         PC2 is interactive. This significantly alters the dthek and
!         dqek terms. In order to calculate reasonable mass
!         flux profiles from the CAPE closure (which is designed
!         without condensate in mind) we need to use values
!         of dthek and dqek that are similar to those calculated
!         assuming that PC2 is not interactive.
          TMP_DCPBYDT=(DTHEK_NONPC2(I)*(1.0+C_VIRTUAL*QEK(I))+          &
     &               C_VIRTUAL*THEK(I)*DQEK_NONPC2(I))*                 &
     &               (DELPK(I)/(RHO*THVE))
         ELSE
!         PC2 is switched off or not interactive. This is the original
!         code.
          TMP_DCPBYDT=(DTHEK(I)*(1.0+C_VIRTUAL*QEK(I))+                 &
     &               C_VIRTUAL*THEK(I)*DQEK(I))*                        &
     &               (DELPK(I)/(RHO*THVE))
         END IF
!
         IF (TMP_DCPBYDT  >   0.0) THEN
          DCPBYDT(I) = DCPBYDT(I) + TMP_DCPBYDT
         ENDIF
!
         IF(BTERM(I))THEN
          THVP=THPKP1(I)*(1.0+C_VIRTUAL*QPKP1(I))
          THVE=THEKP1(I)*(1.0+C_VIRTUAL*QEKP1(I))
          RHO=PKP1(I)/(R*THEKP1(I)*EXKP1(I))
!
          CAPE(I)=CAPE(I)+(THVP-THVE)*DELPKP1(I)/(RHO*THVE)
          RELH(I)=RELH(I)+(QEKP1(I)/QSEKP1(I))*DELPKP1(I)
          DPTOT(I)=DPTOT(I)+DELPKP1(I)
!
          IF (L_Q_INTERACT) THEN
!          PC2 is interactive. This significantly alters the dthek and
!          dqek terms. In order to calculate reasonable mass
!          flux profiles from the CAPE closure (which is designed
!          without condensate in mind) we need to use values
!          of dthek and dqek that are similar to those calculated
!          assuming that PC2 is not interactive.
           TMP_DCPBYDT=(DTHEKP1_NONPC2(I)                               &
     &               *(1.0+C_VIRTUAL*QEKP1(I))+                         &
     &               C_VIRTUAL*THEKP1(I)*DQEKP1_NONPC2(I))*             &
     &               (DELPKP1(I)/(RHO*THVE))
          ELSE
!          PC2 is switched off or not interactive. This is the original
!          code.
           TMP_DCPBYDT=(DTHEKP1(I)*(1.0+C_VIRTUAL*QEKP1(I))+            &
     &               C_VIRTUAL*THEKP1(I)*DQEKP1(I))*                    &
     &               (DELPKP1(I)/(RHO*THVE))
          END IF
!
          IF (TMP_DCPBYDT  >   0.0) THEN
           DCPBYDT(I) = DCPBYDT(I) + TMP_DCPBYDT
          END IF
         END IF
!
       END DO
      ELSE !CAPE_OPT
       DO I=1,NPNTS
         THVP=THPK(I)*(1.0+C_VIRTUAL*QPK(I))
         THVE=THEK(I)*(1.0+C_VIRTUAL*QEK(I))
         RHO=PK(I)/(R*THEK(I)*EXK(I))
!
         CAPE(I)=CAPE(I)+(THVP-THVE)*DELPK(I)/(RHO*THVE)
         RELH(I)=RELH(I)+(QEK(I)/QSEK(I))*DELPK(I)
         DPTOT(I)=DPTOT(I)+DELPK(I)
!
         IF (L_Q_INTERACT) THEN
!         PC2 is interactive. This significantly alters the dthek and
!         dqek terms. In order to calculate reasonable mass
!         flux profiles from the CAPE closure (which is designed
!         without condensate in mind) we need to use values
!         of dthek and dqek that are similar to those calculated
!         assuming that PC2 is not interactive.
          DCPBYDT(I)=DCPBYDT(I)+                                        &
     &               (DTHEK_NONPC2(I)*(1.0+C_VIRTUAL*QEK(I))+           &
     &               C_VIRTUAL*THEK(I)*DQEK_NONPC2(I))*                 &
     &               (DELPK(I)/(RHO*THVE))
         ELSE
!         PC2 is switched off or not interactive. This is the original
!         code.
          DCPBYDT(I)=DCPBYDT(I)+(DTHEK(I)*(1.0+C_VIRTUAL*QEK(I))+       &
     &               C_VIRTUAL*THEK(I)*DQEK(I))*                        &
     &               (DELPK(I)/(RHO*THVE))
         END IF
!
         IF(BTERM(I))THEN
          THVP=THPKP1(I)*(1.0+C_VIRTUAL*QPKP1(I))
          THVE=THEKP1(I)*(1.0+C_VIRTUAL*QEKP1(I))
          RHO=PKP1(I)/(R*THEKP1(I)*EXKP1(I))
!
          CAPE(I)=CAPE(I)+(THVP-THVE)*DELPKP1(I)/(RHO*THVE)
          RELH(I)=RELH(I)+(QEKP1(I)/QSEKP1(I))*DELPKP1(I)
          DPTOT(I)=DPTOT(I)+DELPKP1(I)
!
          IF (L_Q_INTERACT) THEN
!          PC2 is interactive. This significantly alters the dthek and
!          dqek terms. In order to calculate reasonable mass
!          flux profiles from the CAPE closure (which is designed
!          without condensate in mind) we need to use values
!          of dthek and dqek that are similar to those calculated
!          assuming that PC2 is not interactive.
           DCPBYDT(I)=DCPBYDT(I)+(DTHEKP1_NONPC2(I)                     &
     &               *(1.0+C_VIRTUAL*QEKP1(I))+                         &
     &               C_VIRTUAL*THEKP1(I)*DQEKP1_NONPC2(I))*             &
     &               (DELPKP1(I)/(RHO*THVE))
          ELSE
!          PC2 is switched off or not interactive. This is the original
!          code.
           DCPBYDT(I)=DCPBYDT(I)+(DTHEKP1(I)*(1.0+C_VIRTUAL*QEKP1(I))+  &
     &               C_VIRTUAL*THEKP1(I)*DQEKP1(I))*                    &
     &               (DELPKP1(I)/(RHO*THVE))
          END IF
         END IF
!
       END DO
      END IF ! CAPE_OPT
!
!L
!L---------------------------------------------------------------------
!L SWAP PARCEL VALUES READY FOR THE NEXT PART OF ASCENT
!L FROM LAYER K+1 TO K+2 (LOOPING OVER NO. OF TRACERS)
!L---------------------------------------------------------------------
!L
        DO 30 I=1,NPNTS
            THPK(I) = THPKP1(I)
            QPK(I) = QPKP1(I)
            XPK(I) = XPKP1(I)
            QCLPK(I) = QCLPKP1(I)
            QCFPK(I) = QCFPKP1(I)
            FLXK(I) = FLXKP1(I)
            BGMK(I) = BGMKP1(I)
 30     CONTINUE
!
        IF(L_MOM_GK)THEN
          DO I=1,NPNTS
            UPK(I) = UPKP1(I)
            VPK(I) = VPKP1(I)
          END DO
        END IF
!
          IF(L_TRACER)THEN
!
          DO KTRA = 1,NTRA
            DO I=1,NPNTS
              IF(BINIT(I))THEN
                TRAPK(I,KTRA) = TRAPKP1(I,KTRA)
              END IF
            END DO
          END DO
!
          END IF
!
      RETURN
      END SUBROUTINE CONVEC2
