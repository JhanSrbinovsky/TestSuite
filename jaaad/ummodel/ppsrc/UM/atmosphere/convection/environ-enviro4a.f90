
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE ENVIRON------------------------------------------------
!LL
!LL  PURPOSE : CALCULATE THE EFFECT OF CONVECTION UPON THE
!LL            LARGE-SCALE ATMOSPHERE
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
      SUBROUTINE ENVIRON (K,NPNTS,NP_FULL,DTHEK,DQEK,DTHEKP1,DQEKP1,    &
     &                    DQCLEK, DQCFEK, DQCLEKP1, DQCFEKP1, DBCFEK,   &
     &                    DCFLEK, DCFFEK, DBCFEKP1, DCFLEKP1, DCFFEKP1, &
     &                    DQEK_NONPC2, DTHEK_NONPC2,                    &
     &                    DQEKP1_NONPC2, DTHEKP1_NONPC2,                &
     &                    THEK,QEK,QCLEK,QCFEK,BCFEK,CFLEK,CFFEK,       &
     &                    DELTAK, FLXK, THPK, QPK, QCLPK, QCFPK,        &
     &                    THRK,QRK,THEKP1,QEKP1,                        &
     &                    QCLEKP1,QCFEKP1,BCFEKP1,CFLEKP1,CFFEKP1,      &
     &                    QCLPKP1,QCFPKP1,BTERM,THPKP1,                 &
     &                    QPKP1,XPK,XPKP1,BWK,BWKP1,FLXKP1,BLOWST,      &
     &                    EKP14,EXK,EXKP1,DELPK,DELPKP1,                &

! Add new layer pressures for uv layers on charney-phillips grid
     &                    DELP_UV_K, DELP_UV_KP1,                       &
     &                    THPIXS_V,QPIXS_V,                             &
     &                    AMDETK,T1_SD,Q1_SD,                           &
     &                    L_MOM_GK,DUEK,DVEK,DUEKP1,DVEKP1,UEK,VEK,     &
     &                    UPK,VPK,UEKP1,VEKP1,UPKP1,VPKP1,EFLUX_U_UD,   &
     &                    EFLUX_V_UD,L_SHALLOW,CUMULUS,                 &
     &                    L_MID,L_TRACER,NTRA,DTRAEK,DTRAEKP1,          &
     &                    TRAEK,TRAPK,TRAEKP1,TRAPKP1,                  &
     &                    MAX_CFL_C, TIMESTEP,L_CALC_DXEK,L_Q_INTERACT, &
     &                    SDET_ON                                       &
     &                   )

!
      Use cv_run_mod, Only:                                             &
          l_sdxs, l_xscomp, bl_cnv_mix

      IMPLICIT NONE
!
!-----------------------------------------------------------------------
! MODEL CONSTANTS
!-----------------------------------------------------------------------
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
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
! PARXS start
      ! initial excess potential temperature (k) and mixing ratio
      ! (kg/kg) for deep convection
      REAL, PARAMETER :: THPIXS_DEEP= 0.2
      REAL, PARAMETER :: QPIXS_DEEP =0.0

      ! initial excess potential temperature (k) and mixing ratio
      ! (kg/kg) for shallow convection
      REAL, PARAMETER :: THPIXS_SHALLOW = 0.2
      REAL, PARAMETER :: QPIXS_SHALLOW  = 0.0

      ! initial excess potential temperature (k) and mixing ratio
      ! (kg/kg) for mid-level convection
      REAL, PARAMETER :: THPIXS_MID= 0.2
      REAL, PARAMETER :: QPIXS_MID =0.0
! PARXS end
!
!-----------------------------------------------------------------------
! VECTOR LENGTHS AND LOOP COUNTERS
!-----------------------------------------------------------------------
!
      INTEGER NPNTS          ! IN VECTOR LENGTH
!
      INTEGER NP_FULL        ! IN FULL VECTOR LENGTH
!
      INTEGER NTRA           ! IN NUMBER OF TRACERS
!
      INTEGER I,KTRA         ! LOOP COUNTERS
!
      INTEGER K              ! IN NUMBER OF MODEL LEVELS
!
!-----------------------------------------------------------------------
! VARIABLES THAT ARE INPUT
!-----------------------------------------------------------------------
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
      REAL UEK(NPNTS)        ! IN ENVIRONMENT U IN LAYER K (M/S)
!                            !    (UV levels)
!
      REAL UEKP1(NPNTS)      ! IN ENVIRONMENT U IN LAYER K+1 (M/S)
!
      REAL VEK(NPNTS)        ! IN ENVIRONMENT V IN LAYER K (M/S)
!                            !    (UV levels)
!
      REAL VEKP1(NPNTS)      ! IN ENVIRONMENT V IN LAYER K+1 (M/S)
!
      REAL TRAEK(NP_FULL,                                               &
                             ! IN TRACER OF CLOUD ENVIRONMENT
     &           NTRA)       !    IN LAYER K (KG/KG)
!
      REAL TRAEKP1(NP_FULL,                                             &
                             ! IN TRACER OF CLOUD ENVIRONMENT
     &             NTRA)     !    IN LAYER K+1 (KG/KG)
!
      REAL THPK(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K (K)
!
      REAL QPK(NPNTS)        ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
!
      REAL UPK(NPNTS)        ! IN PARCEL U IN LAYER K (M/S)
!                            !    (UV levels)
!
      REAL VPK(NPNTS)        ! IN PARCEL V IN LAYER K (M/S)
!                            !    (UV levels)
!
      REAL TRAPK(NP_FULL,                                               &
                             ! IN PARCEL TRACER IN LAYER K (KG/KG)
     &           NTRA)
!
      REAL THPKP1(NPNTS)     ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K+1 (K)
!
      REAL QPKP1(NPNTS)      ! IN PARCEL MIXING RATIO IN LAYER K+1
                             !    (KG/KG)
!
      REAL UPKP1(NPNTS)      ! IN PARCEL U IN LAYER K+1 (M/S)
!
      REAL VPKP1(NPNTS)      ! IN PARCEL V IN LAYER K+1 (M/S)
!
      REAL TRAPKP1(NP_FULL,                                             &
                             ! IN PARCEL TRACER IN LAYER K+1
     &             NTRA)     !    (KG/KG)
!
!          Xpk should be redundant if we set L_q_interact .True
      REAL XPK(NPNTS)        ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
!
      REAL QCLPK(NPNTS)      ! IN PARCEL LIQUID CONDENSATE MIXING RATIO
!                                 IN LAYER K (KG/KG)
!
      REAL QCFPK(NPNTS)      ! IN PARCEL FROZEN CONDENSATE MIXING RATIO
!                                 IN LAYER K (KG/KG)
!
      REAL QCLPKP1(NPNTS)    ! IN PARCEL LIQUID CONDENSATE MIXING RATIO
!                                 IN LAYER K+1 (KG/KG)
!
      REAL QCFPKP1(NPNTS)    ! IN PARCEL FROZEN CONDENSATE MIXING RATIO
!                                 IN LAYER K+1 (KG/KG)
!
      REAL FLXK(NPNTS)       ! IN PARCEL MASSFLUX IN LAYER K (PA/S)
!
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
      LOGICAL BWK(NPNTS)     ! IN  And layer K
!
      LOGICAL BTERM(NPNTS)   ! IN MASK FOR PARCELS WHICH TERMINATE IN
                             !    LAYER K+1
!
      LOGICAL BLOWST(NPNTS)  ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
!
      LOGICAL L_SHALLOW(NPNTS),                                         &
                                !IN SWITCHES FOR TYPE OF CONVECTION
     &          L_MID(NPNTS)    !   LIKELY TO DEVELOP
!
      LOGICAL CUMULUS(NPNTS) ! IN CUMULUS CONVECTION FLAG (DEEP
                             !    AND SHALLOW)
!
      LOGICAL L_TRACER       ! IN SWITCH FOR INCLUSION OF TRACERS

      LOGICAL L_MOM_GK       ! IN LOGICAL SWITCH FOR INCLUSION OF
                             !    MOMENTUM TRANSPORT (Gregory-Kershaw)

      LOGICAL L_CALC_DXEK    ! IN Switch for calculating condensate inc.
!
      LOGICAL L_Q_INTERACT   ! IN Switch allows overwriting of parcel
!                                 variables when calculating condensate
!                                 increments (will alter results).
!
      INTEGER SDET_ON        ! IN Flag smoothed forced detrainment
!
      REAL TIMESTEP          ! IN CONVECTION TIMESTEP (SECS)
!                                (= model timestep if conv. called once)
!
      REAL THRK(NPNTS)       ! IN PARCEL DETRAINMENT POTENTIAL
                             !    TEMPERATURE IN LAYER K (K)
!
      REAL QRK(NPNTS)        ! IN PARCEL DETRAINMENT MIXING RATIO
                             !    IN LAYER K (KG/KG)
!
!          Xpkp1 should be redundant if we set L_q_interact .True.
      REAL XPKP1(NPNTS)      ! IN PARCEL CLOUD WATER IN LAYER K+1
                             !    (KG/KG)
!
      REAL FLXKP1(NPNTS)     ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
!
      REAL DELTAK(NPNTS)     ! IN PARCEL FORCED DETRAINMENT RATE
                             !    IN LAYER K MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
!
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT RATE FOR LEVEL K+1/4
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
!
      REAL EXK(NPNTS)        ! IN EXNER RATIO FOR MID-POINT OF LAYER K
!
      REAL EXKP1(NPNTS)      ! IN EXNER RATIO FOR MID-POINT OF
                             !    LAYER K+1
!
      REAL DELPK(NPNTS)      ! IN PRESSURE DIFFERENCE ACROSS LAYER K
                             !    (PA)
!
      REAL DELPKP1(NPNTS)    ! IN PRESSURE DIFFERENCE ACROSS LAYER K+1
                             !    (PA)

!
      REAL DELP_UV_K(NPNTS)  ! IN PRESSURE DIFFERENCE ACROSS UV LAYER K
                             !    (ie: model level k-1/2) (PA)
!
      REAL DELP_UV_KP1(NPNTS) ! IN PRESSURE DIFFERENCE ACROSS UV LAYER
                              !    K+1 (ie: model level k+1/2)   (PA)
!
      REAL THPIXS_V(NPNTS)    ! IN  PARCEL EXCESSES OF THETA
      REAL QPIXS_V(NPNTS)     ! AND Q

!
      REAL AMDETK(NPNTS)     ! IN MIXING DETRIANMENT AT LEVEL K
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
!
      REAL T1_SD(NPNTS)      ! IN Standard deviation of turbulent
!                            !    fluctuations of layer 1
!                            !    temperature (K).
      REAL Q1_SD(NPNTS)      ! IN Standard deviation of turbulent
!                            !    fluctuations of layer 1
!                            !    humidity (kg/kg).
!
!-----------------------------------------------------------------------
! VARIABLES THAT ARE INPUT AND OUTPUT
!-----------------------------------------------------------------------
!
      REAL DTHEK(NPNTS)      ! INOUT
                             ! IN  INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (MAY BE NONE ZERO
                             !     DUE TO A PREVIOUS SPLIT FINAL
                             !     DETRAINMENT CALCULATION) (K/S)
                             ! OUT UPDATED INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (K/S)
!
      REAL DQEK(NPNTS)       ! INOUT
                             ! IN  INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K DUE TO CONVECTION
                             !     (MAY BE NONE ZERO
                             !     DUE TO A PREVIOUS SPLIT FINAL
                             !     DETRAINMENT CALCULATION) (KG/KG/S)
                             ! OUT UPDATED INCREMENT TO MODEL MIXING
                             !     RATIO IN LAYER K DUE TO
                             !     CONVECTION (KG/KG/S)
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
                             ! IN  INCREMENT TO MODEL U DUE TO
                             !     CONVECTION (M/S)
                             ! OUT UPDATED INCREMENT TO MODEL U
                             !     DUE TO CONVECTION
!                            !     (at UV level K)
!
      REAL DVEK(NPNTS)       ! INOUT
                             ! IN  INCREMENT TO MODEL V DUE TO
                             !     CONVECTION (M/S)
                             ! OUT UPDATED INCREMENT TO MODEL V
                             !     DUE TO CONVECTION
!                            !     (at UV level K)
!
      REAL DTRAEK(NP_FULL,                                              &
                             ! INOUT
     &            NTRA)      ! IN INCREMENT TO MODEL TRACER IN
                             !    LAYER K DUE TO CONVECTION
                             !    (MAY BE NON ZERO DUE TO
                             !    A PREVIOUS SPLIT FINAL DETRAINMENT
                             !    CALCULATION (KG/KG/S)
                             ! OUT UPDATED INCREMENT TO MODEL TRACER
                             !    IN LAYER K DUE TO CONVECTION
                             !    (KG/KG/S)

!
!-----------------------------------------------------------------------
! VARIABLES THAT ARE OUTPUT
!-----------------------------------------------------------------------
!
      REAL DTHEKP1(NPNTS)    ! OUT INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K+1 DUE TO
                             !     CONVECTION (K/S)
!
      REAL DQEKP1(NPNTS)     ! OUT INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K+1 DUE TO CONVECTION
                             !     (KG/KG/S)
!
      REAL DQEK_NONPC2(NPNTS)   ! OUT Increment to qek if PC2 was not
!                                     in place
      REAL DTHEK_NONPC2(NPNTS)  ! OUT Increment to thek if PC2 was not
!                                     in place
      REAL DQEKP1_NONPC2(NPNTS) ! OUT Increment to qekp1 if PC2 was not
!                                     in place
      REAL DTHEKP1_NONPC2(NPNTS)! OUT Increment to thekp1 if PC2 was not
!                                     in place
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
      REAL DUEKP1(NPNTS)     ! OUT INCREMENT TO MODEL U IN UV LAYER K+1
                             !     DUE TO CONVECTION
!
      REAL DVEKP1(NPNTS)     ! OUT INCREMENT TO MODEL V IN UV LAYER K+1
                             !     DUE TO CONVECTION
!
      REAL DTRAEKP1(NP_FULL,                                            &
                             ! OUT INCREMENT TO MODEL TRACER
     &              NTRA)    !     IN LAYER K+1 DUE TO CONVECTION
                             !     (KG/KG)
!
      REAL EFLUX_U_UD(NPNTS),                                           &
                                ! INOUT
     &     EFLUX_V_UD(NPNTS)    ! IN  EDDY FLUX OF MOMENTUM AT BOTTOM
                                !     OF A LAYER DUE TO UD
                                ! OUT EDDY FLUX OF MOMENTUM AT CURRENT
                                !     LAYER DUE TO UD

      REAL MAX_CFL_C(NPNTS)     ! OUT CFL RATIO
!
!
!-----------------------------------------------------------------------
! VARIABLES THAT ARE DEFINED LOCALLY
!-----------------------------------------------------------------------
!
      REAL FD_DTHEK(NPNTS)   ! THE FORCED DETRAINMENT P. TEMP INC AT LEVEL
                             ! K USED WITH SMOOTHED FORCED DETRAINMENT.
!
      REAL TMP_FD_DTHEK      ! THE FORCED DETRAINMENT P. TEMP INC AT LEVEL
                                ! K USED WITH SMOOTHED FORCED DETRAINMENT.
!
      REAL FD_DTHEKP1(NPNTS) ! THE FORCED DETRAINMENT P. TEMP INC AT LEVEL
                             ! K+1 USED WITH SMOOTHED FORCED DETRAINMENT.
!
      REAL FD_DQEK(NPNTS)    ! THE FORCED DETRAINMENT HUMIDITY INC AT LEVEL
                             ! K USED WITH SMOOTHED FORCED DETRAINMENT.
!
      REAL TMP_FD_DQEK       ! THE FORCED DETRAINMENT HUMIDITY INC AT LEVEL
                             ! K USED WITH SMOOTHED FORCED DETRAINMENT.
!
      REAL FD_DQEKP1(NPNTS)  ! THE FORCED DETRAINMENT HUMIDITY INC AT LEVEL
                             ! K+1 USED WITH SMOOTHED FORCED DETRAINMENT.
!
      REAL EL                ! LATENT HEAT OF CONDENSATION OR
                             ! (CONDENSATION + FUSION) (J/KG)
!
      REAL TEMPRY            ! TEMPORARY ARRAY
!
      REAL TEMP_COM_SUB      ! TEMPORARY ARRAY scalar expands to vector.
!
      REAL TEMP_FOR_DET      ! TEMPORARY ARRAY scalar expands to vector.
!
      REAL THPIXS,QPIXS      ! PARCEL EXCESS POTENTIAL TEMP(K)
                             ! AND MOISTURE(KG/KG)
!
      REAL FLX_U_KP0P5       ! FLUX OF ZONAL MOMENTUM IN CLOUD AT
                             ! TOP OF CURRENT UV LAYER
!
      REAL FLX_V_KP0P5       ! FLUX OF MERIDIONAL MOM. IN CLOUD AT
                             ! TOP OF CURRENT UV LAYER
!
      REAL TEMPDCFXK   ! TEMPORARY WORKSPACE scalar expands to vector.
!
      REAL DBYDT_OF_QCLEK    ! Storage space for liquid condensate rate.
!
      REAL DBYDT_OF_QCFEK    ! Storage space for frozen condensate rate.
!
      REAL DBYDT_OF_CFEK     ! Storage space for total cloud increment.
!
      REAL CFMEK_I     ! Storage space for mixed phase cloud on K.
!
      REAL CFMEKP1_I   ! Storage space for mixed phase cloud on K+1.
!
      REAL QCMEK_I     ! Storage space for mixed phase condensate on K.
!
      REAL QCMEKP1_I   ! Storage space for mixed phase condensate K+1.
!
      REAL DQCMEK_I    ! Storage space for mixed phase condensate inc K.
!
      REAL DELTAXL     ! 1 IF CONVECTION PARCEL IS LIQUID, 0 ELSE.
!
      REAL DELTAXF     ! 1 IF CONVECTION PARCEL IS FROZEN, 0 ELSE.
!
      REAL CALC_DXEK   ! 0 IF CALCULATING CONDENSATE INCREMENT, 1 ELSE.
!
      REAL TOLERANCE   ! THRESHOLD TO AVOID DIVIDE BY ZERO ERRORS.
      PARAMETER ( TOLERANCE = 1.0E-10 )
!
      REAL LS0         ! Minimum value for ls - l
      PARAMETER ( LS0 = 5.0E-5)
!
      REAL A_SMTH      ! PARAMETER DETERMINING THE WEIGHTING BETWEEN THE 
                       ! INCREMENTS AT K AND K-1 USED WHEN
                       ! SMOOTHED FORCED DETRAINMENT IS SELECTED.
      PARAMETER (A_SMTH = 0.5)      
      
!
      LOGICAL L_CALC_QCM     ! Switch for mixed phase cloud calc. option
!
!
!*---------------------------------------------------------------------
!
!
      L_CALC_QCM = .TRUE.
!
      IF (L_Q_INTERACT) THEN
        CALC_DXEK = 0.0
!
!-----------------------------------------------------------------------
!     Xpk is probably undefined in this case and may never be used. At
!     some point it would be best to tidy up and remove it, but not if
!     the .false. option is still being retained. Meanwhile OVERWRITE!
!-----------------------------------------------------------------------
!
        DO I=1,NPNTS
!       NB can tell from bwk and bwkp1 which of qcl or qcf is zero.
          XPK(I)   = QCLPK(I)   + QCFPK(I)
          XPKP1(I) = QCLPKP1(I) + QCFPKP1(I)
        END DO
!
      ELSE
        CALC_DXEK = 1.0
      ENDIF
!
      DO I=1,NPNTS
!
!-----------------------------------------------------------------------
!   CREATE A VECTOR OF LATENT HEATS
!-----------------------------------------------------------------------
!
       IF (BWK(I)) THEN
          EL = LC
       ELSE
          EL = LC + LF
       ENDIF
!
!----------------------------------------------------------------------
! CALCULATE PARCEL MASSFLUX DIVIDED BY THE THICKNESS OF LAYER K
! THIS VALUE IS USED IN SEVERAL PLACES IN THE SUBROUTINE
!----------------------------------------------------------------------
!
       TEMPRY = FLXK(I)/DELPK(I)
!
! L_calc_dxek_if1:
        IF (L_CALC_DXEK)  THEN
!
          TEMP_COM_SUB = (1.0 + EKP14(I)) * (1.0 - DELTAK(I)) *         &
     &                                      (1.0 - AMDETK(I))
          TEMP_FOR_DET = DELTAK(I) * (1.0 - AMDETK(I))
!
! ----------------------------------------------------------------------
!         DELTAXL should be zero when QCLPK is zero and 1 otherwise,
!         DELTAXF should be zero when QCFPK is zero and 1 otherwise, as
!         currently, parcel has single phase only. It would be wise to
!         think through implications before allowing both = 1 at once!
! ----------------------------------------------------------------------
!
          IF (BWK(I)) THEN
            DELTAXL = 1.
            DELTAXF = 0.
          ELSE
            DELTAXL = 0.
            DELTAXF = 1.
          ENDIF
!
        END IF  ! L_calc_dxek_if1
!
!
       IF ( BLOWST(I) .AND. L_XSCOMP .AND. (BL_CNV_MIX == 0) ) THEN
!L
!L----------------------------------------------------------------------
!L AT THE LOWEST CONVECTIVE LAYER, THE PARCEL MASS FLUX IS A FLUX FROM
!L THE ENVIRONMENT. IE. THE INITIAL MASS FLUX IS ENTRAINED WITH EXCESS
!L POTENTIAL TEMPERATURE AND MIXING RATIO TPIXS, QPIXS
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (10), EQUATION (39)
!L----------------------------------------------------------------------
!L

         THPIXS= THPIXS_V(I)
         QPIXS = QPIXS_V(I)
!L
         IF ( L_SDXS .AND. K  ==  1 ) THEN
           DTHEK(I) = DTHEK(I) - TEMPRY*MAX(THPIXS , T1_SD(I)/EXK(I))
           DQEK(I) = DQEK(I) - TEMPRY*MAX(QPIXS , Q1_SD(I))
         ELSE
           DTHEK(I) = DTHEK(I) - TEMPRY*THPIXS
           DQEK(I) = DQEK(I) - TEMPRY*QPIXS
         ENDIF
       END IF

!-----------------------------------------------------------------------
! CALCULATE THE INCREMENTS DUE TO FORCED CONVECTION AND IF 
! SMOOTHED FORCED DETRAINMENT IS SELECTED THEN SPLIT ACROSS K AND K+1
!-----------------------------------------------------------------------

       IF (SDET_ON == 1) THEN
!         P. TEMP INCREMENT
          TMP_FD_DTHEK  = DELTAK(I) * (1.0-AMDETK(I)) *                 &
     &           (THRK(I)-THEK(I)-                                      &
     &           (CALC_DXEK*(EL/CP)*XPK(I)/EXK(I)))
!      
          FD_DTHEK(I)   = A_SMTH*TMP_FD_DTHEK
!
          FD_DTHEKP1(I) = (1.0-A_SMTH) * EXK(I)/EXKP1(I) *              &
     &           DELPK(I)/DELPKP1(I)*TMP_FD_DTHEK  
!
!         HUMIDITY INCREMENT
          TMP_FD_DQEK   = DELTAK(I) * (1.0-AMDETK(I)) *                 &
     &           (QRK(I)-QEK(I)+(CALC_DXEK*XPK(I)))                     
!
          FD_DQEK(I)    = A_SMTH*TMP_FD_DQEK
!
          FD_DQEKP1(I)  = (1.0-A_SMTH) *                                &
     &           DELPK(I)/DELPKP1(I)*TMP_FD_DQEK     
       ELSE
!         P. TEMP INCREMENT
          FD_DTHEK(I)   = DELTAK(I) * (1.0-AMDETK(I)) *                 &
     &           (THRK(I)-THEK(I)-                                      &
     &           (CALC_DXEK*(EL/CP)*XPK(I)/EXK(I)))
!
          FD_DTHEKP1(I) = 0.0
!
!         Q INCREMENT     
          FD_DQEK(I)    = DELTAK(I) * (1.0-AMDETK(I)) *                 &
     &           (QRK(I)-QEK(I)+(CALC_DXEK*XPK(I)))                     
!
          FD_DQEKP1(I) = 0.0
       END IF

!L
!L---------------------------------------------------------------------
!L EFFECT OF CONVECTION UPON POTENTIAL TEMPERATURE OF LAYER K
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (10), EQUATION (38A)
! (Modified to remove re-evaporation of condensate to vapour when the
!  effect of convection upon condensate is formally calculated.)
!L--------------------------------------------------------------------
!L
!
       IF (L_Q_INTERACT) THEN
!        This is the old calculation
         DTHEK_NONPC2(I) = DTHEK(I) + TEMPRY * (                        &
      
     &           (1+EKP14(I)) * (1.0-DELTAK(I)) *                       &
                                                         ! COMPENSATING
     &           (1-AMDETK(I)) * (THEKP1(I)-THEK(I))                    &
                                                         ! SUBSIDENCE
     &         +                                                        &
     &           DELTAK(I) * (1.0-AMDETK(I)) *                          &
                                                         ! FORCED
     &           (THRK(I)-THEK(I)-                                      &
                                                         ! DETRAINMENT
     &                    ((EL/CP)*XPK(I)/EXK(I)))                      &
     &         +                                                        &
     &           AMDETK(I) * (THPK(I)-THEK(I)-                          &
                                                         ! MIXING
     &                    ((EL/CP)*XPK(I)/EXK(I)))                      &
                                                         ! DETRAINMENT
     &         )
       END IF
!
!      This is the modified calculation suitable for PC2 and non-PC2
       DTHEK(I) = DTHEK(I) + TEMPRY * (                                 &
      
     &           (1+EKP14(I)) * (1.0-DELTAK(I)) *                       &
                                                         ! COMPENSATING
     &           (1-AMDETK(I)) * (THEKP1(I)-THEK(I))                    &
                                                         ! SUBSIDENCE
     &         +                                                        &
     &           FD_DTHEK(I)                                            &
                                                         ! FORCED
                                                         ! DETRAINMENT
     &         +                                                        &
     &           AMDETK(I) * (THPK(I)-THEK(I)-                          &
                                                         ! MIXING
     &              (CALC_DXEK*(EL/CP)*XPK(I)/EXK(I)))                  &
                                                         ! DETRAINMENT
     &         )
!
       IF (.NOT. L_Q_INTERACT) THEN
!        For safety, copy dthek to dthek_nonpc2 although it is not used
         DTHEK_NONPC2(I) = DTHEK(I)
       END IF

       IF (SDET_ON == 1) THEN
         DTHEKP1(I) = DTHEKP1(I)+TEMPRY*FD_DTHEKP1(I)
       END IF

!
!L
!L---------------------------------------------------------------------
!L EFFECT OF CONVECTION UPON MIXING RATIO OF LAYER K
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (10), EQUATION (38B)
! (Modified to remove re-evaporation of condensate to vapour when the
!  effect of convection upon condensate is formally calculated.)
!L--------------------------------------------------------------------
!L
       IF (L_Q_INTERACT) THEN
!        This is the old calculation
         DQEK_NONPC2(I) = DQEK(I) + TEMPRY * (                          &
      
     &           (1+EKP14(I)) * (1.0-DELTAK(I)) *                       &
                                                         ! COMPENSATING
     &           (1-AMDETK(I)) * (QEKP1(I)-QEK(I))                      &
                                                         ! SUBSIDENCE
     &         +                                                        &
     &           DELTAK(I) * (1.0-AMDETK(I)) *                          &
                                                         ! FORCED
     &           (QRK(I)-QEK(I)+XPK(I))                                 &
                                                         ! DETRAINMENT
     &         +                                                        &
     &           AMDETK(I) * (QPK(I)-QEK(I)+                            &
                                                         ! MIXING
     &                                XPK(I))                           &
                                                         ! DETRAINMENT
     &         )
       END IF
!
!      This is the modified calculation suitable for PC2 and non-PC2
       DQEK(I) = DQEK(I) + TEMPRY * (                                   &
      
     &           (1+EKP14(I)) * (1.0-DELTAK(I)) *                       &
                                                         ! COMPENSATING
     &           (1-AMDETK(I)) * (QEKP1(I)-QEK(I))                      &
                                                         ! SUBSIDENCE
     &         +                                                        &
     &           FD_DQEK(I)                                             &
                                                         ! FORCED
                                                         ! DETRAINMENT
     &         +                                                        &
     &           AMDETK(I) * (QPK(I)-QEK(I)+                            &
                                                         ! MIXING
     &                          (CALC_DXEK*XPK(I)))                     &
                                                         ! DETRAINMENT
     &         )
       IF (.NOT. L_Q_INTERACT) THEN
!        For safety, copy dqek to dqek_nonpc2, although it is not used
         DQEK_NONPC2(I) = DQEK(I)
       END IF


!----------------------------------------------------------------------
! ADD THE FORCED DETRAINMENT INCREMENT TO K+1
!---------------------------------------------------------------------- 
       IF (SDET_ON == 1) THEN
         DQEKP1(I)  = DQEKP1(I)+TEMPRY*FD_DQEKP1(I)
       END IF
      
!
! ----------------------------------------------------------------------
!  OPTIONAL EFFECT OF CONVECTION UPON CONDENSATE MIXING RATIO OF LAYER K
!  (Bushell et al, QJRMS 2003, 129, pp1435-1455)
! ----------------------------------------------------------------------
!
! L_calc_dxek_if2:
        IF (L_CALC_DXEK)  THEN
!
! ----------------------------------------------------------------------
!       Calculate the Q4 rates (liquid, frozen, mixed phase condensate).
! ----------------------------------------------------------------------
!
!     --Increment to Liquid Condensate Mixing Ratio--
          DBYDT_OF_QCLEK = TEMPRY * (                                   &
!                                                    COMPENSATING
     &       TEMP_COM_SUB * (QCLEKP1(I)-QCLEK(I))                       &
                                                   ! SUBSIDENCE
     &      +                                                           &
     &       TEMP_FOR_DET * (QCLPK(I) - QCLEK(I))                       &
                                                   ! FORCED DETRAINMENT
!   Forced detrainment of condensate preserves parcel condensate values.
     &      +                                                           &
     &       AMDETK(I) * (QCLPK(I) - QCLEK(I))                          &
                                                   ! MIXING DETRAINMENT
     &         )
!
!         --------------------------------------------------------------
!         Initial convecting layer - flux through cloud base is zero.
!         --------------------------------------------------------------
          IF (BLOWST(I))  DBYDT_OF_QCLEK = DBYDT_OF_QCLEK -             &
     &                    TEMPRY * (1+EKP14(I)) * (QCLEKP1(I)-QCLEK(I))
!
          DQCLEK(I) = DQCLEK(I) + DBYDT_OF_QCLEK
!
!     --Increment rate to Frozen Condensate Mixing Ratio--
          DBYDT_OF_QCFEK = TEMPRY * (                                   &
!                                                    COMPENSATING
     &       TEMP_COM_SUB * (QCFEKP1(I)-QCFEK(I))                       &
                                                   ! SUBSIDENCE
     &      +                                                           &
     &       TEMP_FOR_DET * (QCFPK(I) - QCFEK(I))                       &
                                                   ! FORCED DETRAINMENT
!   Forced detrainment of condensate preserves parcel condensate values.
     &      +                                                           &
     &       AMDETK(I) * (QCFPK(I) - QCFEK(I))                          &
                                                   ! MIXING DETRAINMENT
     &         )
!
!         --------------------------------------------------------------
!         Initial convecting layer - flux through cloud base is zero.
!         --------------------------------------------------------------
          IF (BLOWST(I))  DBYDT_OF_QCFEK = DBYDT_OF_QCFEK -             &
     &                    TEMPRY * (1+EKP14(I)) * (QCFEKP1(I)-QCFEK(I))
!
          DQCFEK(I) = DQCFEK(I) + DBYDT_OF_QCFEK
!
!     --Optional: Increment rate to Mixed Condensate Mixing Ratio-------
!
! L_calc_qcm_if1:
          IF (L_CALC_QCM) THEN
!
            IF (CFLEK(I)  >   TOLERANCE  .AND.  CFFEK(I)  >   TOLERANCE)&
     &      THEN
              CFMEK_I   = CFLEK(I)   + CFFEK(I)   - BCFEK(I)
              QCMEK_I   = CFMEK_I * ( (QCLEK(I) / CFLEK(I)) +           &
     &                                (QCFEK(I) / CFFEK(I)) )
            ELSE
              CFMEK_I   = 0.0
              QCMEK_I   = 0.0
            END IF
!
          IF (CFLEKP1(I)  >   TOLERANCE .AND. CFFEKP1(I)  >   TOLERANCE)&
     &      THEN
              CFMEKP1_I = CFLEKP1(I) + CFFEKP1(I) - BCFEKP1(I)
              QCMEKP1_I = CFMEKP1_I * ( (QCLEKP1(I) / CFLEKP1(I)) +     &
     &                                  (QCFEKP1(I) / CFFEKP1(I)) )
            ELSE
              CFMEKP1_I = 0.0
              QCMEKP1_I = 0.0
            END IF
!
!           Parcel condensate for mixed phase is zero
            DQCMEK_I = TEMPRY * (                                       &
!                                                    COMPENSATING
     &       TEMP_COM_SUB * (QCMEKP1_I - QCMEK_I)                       &
                                                   ! SUBSIDENCE
     &      -                                                           &
     &       TEMP_FOR_DET * (QCMEK_I)                                   &
                                                   ! FORCED DETRAINMENT
!   Forced detrainment of condensate preserves parcel condensate values.
     &      -                                                           &
     &       AMDETK(I) * (QCMEK_I)                                      &
                                                   ! MIXING DETRAINMENT
     &         )
!
          END IF  ! L_calc_qcm_if1
        END IF  ! L_calc_dxek_if2
!
! ----------------------------------------------------------------------
!  OPTIONAL EFFECT OF CONVECTION UPON CLOUD VOLUME FRACTIONS OF LAYER K
!  Block separate from condensate for clarity, but actually dependent on
!  calculations of Q4 condensation rates.
! ----------------------------------------------------------------------
!
! L_calc_dxek_if3:
        IF (L_CALC_DXEK)  THEN
!
! ----------------------------------------------------------------------
!       Calculate cloud volume rates (liquid, frozen and total).
! ----------------------------------------------------------------------
!
!     --Increment to Liquid Cloud Volume Fraction----
!         Using an explicit form of calculation with enforced limit
          TEMPDCFXK =      DELTAXL * MAX( (QCLPK(I)-QCLEK(I)  ) ,LS0)   &
     &              + (1.0-DELTAXL)*    (          -QCLEK(I)  )
!           DELTAXL should be zero when QCLPK is zero and 1 otherwise.
!
          IF (ABS(TEMPDCFXK)  >   TOLERANCE) THEN
            TEMPDCFXK = DBYDT_OF_QCLEK*(DELTAXL - CFLEK(I)) / TEMPDCFXK
            TEMPDCFXK = Max(TEMPDCFXK, (-1. * CFLEK(I) / TIMESTEP))
            TEMPDCFXK = Min(TEMPDCFXK, ((1. - CFLEK(I)) / TIMESTEP))
          ELSE
            TEMPDCFXK = 0.0                        ! Zero rate
          END IF
!
          DCFLEK(I) = DCFLEK(I) + TEMPDCFXK
!         Initialize total cloud increment
          DBYDT_OF_CFEK = TEMPDCFXK
!
!     --Increment to Frozen Cloud Volume Fraction----
!         Using an explicit form of calculation with enforced limit
          TEMPDCFXK =      DELTAXF * MAX( (QCFPK(I)-QCFEK(I)  ) ,LS0)   &
     &              + (1.0-DELTAXF)*    (          -QCFEK(I)  )
!           DELTAXF should be zero when QCFPK is zero and 1 otherwise.
!
          IF (ABS(TEMPDCFXK)  >   TOLERANCE) THEN
            TEMPDCFXK = DBYDT_OF_QCFEK*(DELTAXF - CFFEK(I)) / TEMPDCFXK
            TEMPDCFXK = Max(TEMPDCFXK, (-1. * CFFEK(I) / TIMESTEP))
            TEMPDCFXK = Min(TEMPDCFXK, ((1. - CFFEK(I)) / TIMESTEP))
          ELSE
            TEMPDCFXK = 0.0                        ! Zero rate
          END IF
!
          DCFFEK(I) = DCFFEK(I) + TEMPDCFXK
          DBYDT_OF_CFEK = DBYDT_OF_CFEK + TEMPDCFXK
!
!     --Increment to Mixed Phase Cloud Volume Fraction----
! L_calc_qcm_if2:
          IF (L_CALC_QCM) THEN
!           Using an explicit form of calculation with enforced limit
            TEMPDCFXK = QCMEK_I
!
              IF (ABS(TEMPDCFXK)   >    TOLERANCE)  THEN
              TEMPDCFXK = DQCMEK_I  * CFMEK_I / TEMPDCFXK
              TEMPDCFXK = Max(TEMPDCFXK, (-1. * CFMEK_I / TIMESTEP))
              TEMPDCFXK = Min(TEMPDCFXK, ((1. - CFMEK_I) / TIMESTEP))
            ELSE
              TEMPDCFXK = 0.0                      ! Zero rate
            END IF
!
! ----------------------------------------------------------------------
!       Total values.
! ----------------------------------------------------------------------
!
!     --Increment to Total  Cloud Volume Fraction----
!       METHOD 1: Summation.
            TEMPDCFXK = DBYDT_OF_CFEK - TEMPDCFXK
          ELSE  ! L_calc_qcm_if2
!
!     --Increment to Total  Cloud Volume Fraction----
!       METHOD 2: Direct specification.
!           Using an explicit form of calculation with enforced limit
            TEMPDCFXK = MAX((QCLPK(I)+QCFPK(I)-QCLEK(I)-QCFEK(I)),LS0)
              IF (ABS(TEMPDCFXK)   >    TOLERANCE)  THEN
              TEMPDCFXK = (DBYDT_OF_QCLEK + DBYDT_OF_QCFEK) *           &
     &                    (1.0 - BCFEK(I)) / TEMPDCFXK
              TEMPDCFXK = Max(TEMPDCFXK, (-1. * BCFEK(I) / TIMESTEP))
              TEMPDCFXK = Min(TEMPDCFXK, ((1. - BCFEK(I)) / TIMESTEP))
              ELSE
              TEMPDCFXK = 0.0                      ! Zero rate
            END IF
!
          END IF  ! L_calc_qcm_if2
!
          DBCFEK(I) = DBCFEK(I) + TEMPDCFXK
!
        END IF  ! L_calc_dxek_if3
!
!
! CALCULATE COURANT NUMBER USING MASS FLUX AT HALF LEVEL
!
          MAX_CFL_C(I)=TEMPRY*(1+EKP14(I))*(1.0-DELTAK(I))*             &
     &                        (1.0-AMDETK(I))
!L
!L----------------------------------------------------------------------
!L TERMINAL DETRAINMENT AND SUBSIDENCE IN TERMINAL LAYER
!L
!L UM DOCUMENTATION PAPER 27
!L SECTION (10), EQUATION (40)
!L--------------------------------------------------------------------
!L
       IF ( BTERM(I) ) THEN
          IF (BWKP1(I)) THEN
             EL = LC
           ELSE
             EL = LC + LF
          ENDIF
          TEMPRY = FLXKP1(I)/DELPKP1(I)
          IF (L_Q_INTERACT) THEN
!           This is the old calculation
            DTHEKP1_NONPC2(I) = DTHEKP1(I)                              &
     &                            + TEMPRY*((THPKP1(I)-THEKP1(I))       &
     &                            - EL*XPKP1(I)/(EXKP1(I)*CP))
            DQEKP1_NONPC2(I)  = DQEKP1(I) + TEMPRY*(QPKP1(I)-QEKP1(I)   &
     &                                                + XPKP1(I))
          END IF
          DTHEKP1(I) = DTHEKP1(I) + TEMPRY*((THPKP1(I)-THEKP1(I))       &
     &                    - CALC_DXEK * EL*XPKP1(I)/(EXKP1(I)*CP))
          DQEKP1(I)  = DQEKP1(I) + TEMPRY*(QPKP1(I)-QEKP1(I)            &
     &                                     + CALC_DXEK * XPKP1(I))
          IF (.NOT. L_Q_INTERACT) THEN
!           For safety, set dthekp1_nonpc2 and dqekp1_nonpc2, even
!           though they aren't used.
            DTHEKP1_NONPC2(I) = DTHEKP1(I)
            DQEKP1_NONPC2(I)  = DQEKP1(I)
          END IF
!
          IF (L_CALC_DXEK)  THEN
            DQCLEKP1(I)  = DQCLEKP1(I) + TEMPRY*(QCLPKP1(I)-QCLEKP1(I))
!
            DQCFEKP1(I)  = DQCFEKP1(I) + TEMPRY*(QCFPKP1(I)-QCFEKP1(I))
!
! Terminal detrainment of liquid
            IF ( (QCLPKP1(I)-QCLEKP1(I))  >   LS0 ) THEN
!         Go ahead with terminal detraiment
              DCFLEKP1(I)= DCFLEKP1(I) + TEMPRY *(DELTAXL - CFLEKP1(I))
            ELSE
!         Scale down terminal detrainment of liquid cloud fraction
              DCFLEKP1(I)= DCFLEKP1(I) + TEMPRY *(DELTAXL - CFLEKP1(I)) &
     &                   * MAX( (QCLPKP1(I)-QCLEKP1(I)),0.0)/LS0
            END IF
!
! Terminal detrainment of ice
            IF ( (QCFPKP1(I)-QCFEKP1(I))  >   LS0 ) THEN
!         Go ahead with terminal detraiment
              DCFFEKP1(I)= DCFFEKP1(I) + TEMPRY *(DELTAXF - CFFEKP1(I))
            ELSE
!         Scale down terminal detrainment of frozen cloud fraction
              DCFFEKP1(I)= DCFFEKP1(I) + TEMPRY *(DELTAXF - CFFEKP1(I)) &
     &                   * MAX( (QCFPKP1(I)-QCFEKP1(I)),0.0)/LS0
            END IF
!
! Terminal detrainment of bulk cloud fraction
            IF ( (QCFPKP1(I)-QCFEKP1(I)+QCLPKP1(I)-QCLEKP1(I))          &
     &             >   LS0 ) THEN
!         Go ahead with terminal detraiment
              DBCFEKP1(I)= DBCFEKP1(I) + TEMPRY * ( 1.0   - BCFEKP1(I))
            ELSE
!         Scale down terminal detrainment of bulk cloud fraction
              DBCFEKP1(I)= DBCFEKP1(I) + TEMPRY *( 1.0 - BCFEKP1(I))    &
     &                   * MAX( (QCFPKP1(I)-QCFEKP1(I)                  &
     &                          +QCLPKP1(I)-QCLEKP1(I)),0.0)/LS0
            END IF
!
! endif l_calc_dxek
          END IF
! endif bterm(i)
       END IF
!L
      END DO
!L
!L---------------------------------------------------------------------
!L CALCULATE EFFECT OF CONVECTION UPON MOMENTUM OF LAYER K
!L AND DO TERMINAL DETRAINMENT OF MOMENTUM
!L
!L RATE OF CHANGE OF WIND FIELD BY CONVECTION IS ESTIMATED USING A
!L DIVERGENCE OF VERTICAL EDDY MOMENTUM FLUX ACROSS THE LAYER
!L--------------------------------------------------------------------
!L
!
! ALL CONVECTIVE MOMENTUM TRANSPORT CALCULATIONS FOR CUMULUS CONVECTION
! (DEEP AND SHALLOW) DONE WHEN CONVECTION TERMINATES
!

      IF(L_MOM_GK) THEN      ! Gregory-Kershaw CMT
!
      DO I=1,NPNTS
!----------------------------------------------------------------------
! ESTIMATE EDDY FLUX AT TOP OF CURRENT UV LAYER DUE TO CONVECTION
!----------------------------------------------------------------------

       FLX_U_KP0P5 = FLXK(I) * (1.0-AMDETK(I)) * (1.0-DELTAK(I)) *      &
     &               (1.0+EKP14(I)) * (UPK(I)-UEKP1(I))
       FLX_V_KP0P5 = FLXK(I) * (1.0-AMDETK(I)) * (1.0-DELTAK(I)) *      &
     &               (1.0+EKP14(I)) * (VPK(I)-VEKP1(I))
!
       IF (BLOWST(I)) THEN
!----------------------------------------------------------------------
! INITIAL CONVECTING LAYER - NO FLUX AT BASE OF LAYER
!----------------------------------------------------------------------
       DUEK(I) = DUEK(I) - FLX_U_KP0P5 / DELP_UV_K(I)
       DVEK(I) = DVEK(I) - FLX_V_KP0P5 / DELP_UV_K(I)
!----------------------------------------------------------------------
! STORE EDDY FLUX AT TOP OF CURRENT UV LAYER READY FOR CALCULATION
! OF NEXT LAYER
!----------------------------------------------------------------------
       EFLUX_U_UD(I) = FLX_U_KP0P5
       EFLUX_V_UD(I) = FLX_V_KP0P5
!
       ELSE
!----------------------------------------------------------------------
! CONVECTING LAYER - TAKE EDDY FLUX DIVERGENCE ACROSS THE LAYER
!----------------------------------------------------------------------
       DUEK(I) = DUEK(I) - ( (FLX_U_KP0P5 - EFLUX_U_UD(I)) /            &
     &                                                DELP_UV_K(I) )
       DVEK(I) = DVEK(I) - ( (FLX_V_KP0P5 - EFLUX_V_UD(I)) /            &
     &                                                DELP_UV_K(I) )
!----------------------------------------------------------------------
! STORE EDDY FLUX AT TOP OF CURRENT UV LAYER READY FOR CALCULATION
! OF NEXT LAYER
!----------------------------------------------------------------------
       EFLUX_U_UD(I) = FLX_U_KP0P5
       EFLUX_V_UD(I) = FLX_V_KP0P5
!
      END IF
!
      IF(BTERM(I))THEN
!----------------------------------------------------------------------
! CONVECTION TERMINATES - CALCULATE INCREMENT DUE TO CONVECTION
! IN TOP LAYER - NO FLUX OUT OF TOP OF LAYER
!----------------------------------------------------------------------
          DUEKP1(I)  = EFLUX_U_UD(I) / DELP_UV_KP1(I)
          DVEKP1(I)  = EFLUX_V_UD(I) / DELP_UV_KP1(I)
!----------------------------------------------------------------------
! ZERO EDDY FLUX OUT OF TOP OF LAYER
!----------------------------------------------------------------------
          EFLUX_U_UD(I) = 0.0
          EFLUX_V_UD(I) = 0.0
!
      END IF

!
      END DO
!
      END IF
!
!L_____________________________________________________________________
!L
!L EFFECT OF CONVECTION ON TRACER CONTENT OF LAYER K
!L (LOOPING OVER NUMBER OF TRACER VARIABLES)
!L AND DO TERMINAL DETRAINMENT OF TRACER
!L_____________________________________________________________________
!L
      IF(L_TRACER)THEN
!L
      DO KTRA = 1,NTRA
      DO I = 1,NPNTS
!L
       TEMPRY = FLXK(I)/DELPK(I)
       DTRAEK(I,KTRA) = DTRAEK(I,KTRA) + TEMPRY * (                     &
      
     & (1+EKP14(I)) * (1.0-DELTAK(I)) *                                 &
                                                        ! COMPENSATING
     & (1-AMDETK(I)) * (TRAEKP1(I,KTRA)-TRAEK(I,KTRA))                  &
                                                        ! SUBSIDENCE
     &+                                                                 &
     & DELTAK(I) * (1.0-AMDETK(I)) *                                    &
                                                        ! FORCED
     & (TRAPK(I,KTRA)-TRAEK(I,KTRA))                                    &
                                                        ! DETRAINMENT
     &+                                                                 &
     & AMDETK(I) * (TRAPK(I,KTRA)-TRAEK(I,KTRA))                        &
                                                        ! MIXING
     &           )                                      ! DETRAINMENT
!L
!
        IF(BTERM(I))THEN
          TEMPRY = FLXKP1(I)/DELPKP1(I)
          DTRAEKP1(I,KTRA) = DTRAEKP1(I,KTRA) +TEMPRY*                  &
     &                       (TRAPKP1(I,KTRA)-TRAEKP1(I,KTRA))
        END IF
!
      END DO
!
      END DO
!
      END IF

      RETURN
      END SUBROUTINE ENVIRON
