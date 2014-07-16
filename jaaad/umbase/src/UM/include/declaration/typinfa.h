#if defined(ATMOS)
! TYPINFA start
! 10/11/00 5.2 Add LBC_ETA_THETA and LBC_ETA_RHO and rename
!              MAX_INTF_P_LEVELS to MAX_INTF_MODEL_LEVELS. D.Robinson
! 18/08/04 6.1 Add AV_INDEX_* and AV_WEIGHT_*. D Robinson
!
! This file needs files TYPESIZ, TYPESIZA and TYPSIZO to be included

      ! Headers for atmosphere interface data sets
      INTEGER::FIXHD_INTFA(LEN_FIXHD,N_INTF_A)        ! Fixed header
      INTEGER::INTHD_INTFA(PP_LEN_INTHD,N_INTF_A)     ! Integer header
      INTEGER::LOOKUP_INTFA(LEN1_LOOKUP,INTF_LOOKUPSA,N_INTF_A)! Lookups

      REAL :: REALHD_INTFA(PP_LEN_REALHD,N_INTF_A)   ! Real header
      REAL :: LEVDEPC_INTFA(MAX_INTF_MODEL_LEVELS,INTF_LEN2_LEVDEPC,    &
     &  N_INTF_A)
      REAL::ROWDEPC_INTFA(MAX_LBCROWS,INTF_LEN2_ROWDEPC,N_INTF_A)
      REAL::COLDEPC_INTFA(MAX_LBCROW_LENGTH,INTF_LEN2_COLDEPC,N_INTF_A)

      ! Interpolation constants for atmosphere interface data sets.
      !                              Index of corner in source grid box:
      INTEGER::AP_INDEX_B_L(TOT_LEN_INTFA_P) ! Bottom left  ( p grid)
      INTEGER::AP_INDEX_B_R(TOT_LEN_INTFA_P) ! Bottom right ( p grid)
      INTEGER::AU_INDEX_B_L(TOT_LEN_INTFA_U) ! Bottom left  ( u grid)
      INTEGER::AU_INDEX_B_R(TOT_LEN_INTFA_U) ! Bottom right ( u grid)
      INTEGER::AV_INDEX_B_R(TOT_LEN_INTFA_U) ! Bottom right ( v grid)
      INTEGER::AV_INDEX_B_L(TOT_LEN_INTFA_U) ! Bottom right ( v grid)
      !                                      Weight applied to value at:
      REAL :: AP_WEIGHT_T_R(TOT_LEN_INTFA_P) ! Top    right (p grid)
      REAL :: AP_WEIGHT_B_L(TOT_LEN_INTFA_P) ! Bottom left  (p grid)
      REAL :: AP_WEIGHT_B_R(TOT_LEN_INTFA_P) ! Bottom right (p grid)
      REAL :: AP_WEIGHT_T_L(TOT_LEN_INTFA_P) ! Top    left  (p grid)
      REAL :: AU_WEIGHT_T_R(TOT_LEN_INTFA_U) ! Top    right (u grid)
      REAL :: AU_WEIGHT_B_L(TOT_LEN_INTFA_U) ! Bottom left  (u grid)
      REAL :: AU_WEIGHT_B_R(TOT_LEN_INTFA_U) ! Bottom right (u grid)
      REAL :: AU_WEIGHT_T_L(TOT_LEN_INTFA_U) ! Top    left  (u grid)
      REAL :: AV_WEIGHT_T_R(TOT_LEN_INTFA_U) ! Top    right (v grid)
      REAL :: AV_WEIGHT_B_L(TOT_LEN_INTFA_U) ! Bottom left  (v grid)
      REAL :: AV_WEIGHT_B_R(TOT_LEN_INTFA_U) ! Bottom right (v grid)
      REAL :: AV_WEIGHT_T_L(TOT_LEN_INTFA_U) ! Top    left  (v grid)

      ! Rotation coefficients for atmosphere interface data sets
      REAL :: COEFF1(TOT_LEN_INTFA_U)
      REAL :: COEFF2(TOT_LEN_INTFA_U)
      REAL :: COEFF3(U_FIELD_INTFA)
      REAL :: COEFF4(U_FIELD_INTFA)

      ! Vertical levels for each area
      REAL :: INTF_AKH(MAX_INTF_MODEL_LEVELS+1,N_INTF_A)
      REAL :: INTF_BKH(MAX_INTF_MODEL_LEVELS+1,N_INTF_A)
      REAL :: INTF_AK (MAX_INTF_MODEL_LEVELS  ,N_INTF_A)
      REAL :: INTF_BK (MAX_INTF_MODEL_LEVELS  ,N_INTF_A)

      ! Eta Values for each area
      REAL :: LBC_ETA_THETA (MAX_INTF_MODEL_LEVELS+1,N_INTF_A)
      REAL :: LBC_ETA_RHO   (MAX_INTF_MODEL_LEVELS  ,N_INTF_A)
! TYPINFA end
#endif
