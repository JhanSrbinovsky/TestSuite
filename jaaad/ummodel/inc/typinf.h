! TYPINF start
      ! Headers for interface data sets
      INTEGER :: FIXHD_INTF(LEN_FIXHD,N_INTF)        ! Fixed header
      INTEGER :: INTHD_INTF(PP_LEN_INTHD,N_INTF)     ! Integer header
      INTEGER :: LOOKUP_INTF(LEN1_LOOKUP,INTF_LOOKUPS,N_INTF) ! Lookups

      REAL :: REALHD_INTF(PP_LEN_REALHD,N_INTF)   ! Real header
      REAL :: LEVDEPC_INTF(MAX_INTF_P_LEVELS,INTF_LEN2_LEVDEPC,N_INTF)
      REAL ::ROWDEPC_INTFA(MAX_LBCROWS,INTF_LEN2_ROWDEPC,N_INTF_A)
      REAL ::COLDEPC_INTFA(MAX_LBCROW_LENGTH,INTF_LEN2_COLDEPC,N_INTF_A)
 
      ! Interpolation constants for interface data sets.
      !                             Index of corner in source grid box:
      INTEGER :: P_INDEX_B_L(TOT_LEN_INTF_P) ! Bottom left  ( p grid)
      INTEGER :: P_INDEX_B_R(TOT_LEN_INTF_P) ! Bottom right ( p grid)
      INTEGER :: U_INDEX_B_L(TOT_LEN_INTF_U) ! Bottom left  ( u grid)
      INTEGER :: U_INDEX_B_R(TOT_LEN_INTF_U) ! Bottom right ( u grid)
      !                                      Weight applied to value at:
      INTEGER :: P_WEIGHT_T_R(TOT_LEN_INTF_P) ! Top    right (p grid)
      INTEGER :: P_WEIGHT_B_L(TOT_LEN_INTF_P) ! Bottom left  (p grid)
      INTEGER :: P_WEIGHT_B_R(TOT_LEN_INTF_P) ! Bottom right (p grid)
      INTEGER :: P_WEIGHT_T_L(TOT_LEN_INTF_P) ! Top    left  (p grid)
      INTEGER :: U_WEIGHT_T_R(TOT_LEN_INTF_U) ! Top    right (u grid)
      INTEGER :: U_WEIGHT_B_L(TOT_LEN_INTF_U) ! Bottom left  (u grid)
      INTEGER :: U_WEIGHT_B_R(TOT_LEN_INTF_U) ! Bottom right (u grid)
      INTEGER :: U_WEIGHT_T_L(TOT_LEN_INTF_U) ! Top    left  (u grid)

      ! Rotation coefficients for atmosphere interface data sets
      REAL :: COEFF1(TOT_LEN_INTF_U)
      REAL :: COEFF2(TOT_LEN_INTF_U)
      REAL :: COEFF3(NPTS_U_FIELD)
      REAL :: COEFF4(NPTS_U_FIELD)
      REAL :: COEFF5(TOT_LEN_INTF_U)
      REAL :: COEFF6(TOT_LEN_INTF_U)
! TYPINF end
