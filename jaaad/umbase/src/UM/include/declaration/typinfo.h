#if defined(OCEAN)
! TYPINFO needs TYPESIZ and TYPESIZO to be included
      ! Headers for atmosphere interface data sets
      INTEGER::FIXHD_INTFO(LEN_FIXHD,N_INTF_O)        ! Fixed header
      INTEGER::INTHD_INTFO(PP_LEN_INTHD,N_INTF_O)     ! Integer header
      INTEGER::LOOKUP_INTFO(LEN1_LOOKUP,INTF_LOOKUPSO,N_INTF_O)!Lookups

      REAL :: REALHD_INTFO(PP_LEN_REALHD,N_INTF_O)   ! Real header
      REAL :: LEVDEPC_INTFO(MAX_INTF_P_LEVELS_O,INTF_LEN2_LEVDEPC_O,    &
     &    N_INTF_O)

      ! Interpolation constants for atmosphere interface data sets.
      !                              Index of corner in source grid box:
      INTEGER::OP_INDEX_B_L(TOT_LEN_INTFO_P) ! Bottom left  ( p grid)
      INTEGER::OP_INDEX_B_R(TOT_LEN_INTFO_P) ! Bottom right ( p grid)
      INTEGER::OU_INDEX_B_L(TOT_LEN_INTFO_U) ! Bottom left  ( u grid)
      INTEGER::OU_INDEX_B_R(TOT_LEN_INTFO_U) ! Bottom right ( u grid)
      !                                      Weight applied to value at:
      REAL :: OP_WEIGHT_T_R(TOT_LEN_INTFO_P) ! Top    right (p grid)
      REAL :: OP_WEIGHT_B_L(TOT_LEN_INTFO_P) ! Bottom left  (p grid)
      REAL :: OP_WEIGHT_B_R(TOT_LEN_INTFO_P) ! Bottom right (p grid)
      REAL :: OP_WEIGHT_T_L(TOT_LEN_INTFO_P) ! Top    left  (p grid)
      REAL :: OU_WEIGHT_T_R(TOT_LEN_INTFO_U) ! Top    right (u grid)
      REAL :: OU_WEIGHT_B_L(TOT_LEN_INTFO_U) ! Bottom left  (u grid)
      REAL :: OU_WEIGHT_B_R(TOT_LEN_INTFO_U) ! Bottom right (u grid)
      REAL :: OU_WEIGHT_T_L(TOT_LEN_INTFO_U) ! Top    left  (u grid)

      ! Rotation coefficients for ocean interface data sets
      REAL :: O_COEFF1(TOT_LEN_INTFO_U)
      REAL :: O_COEFF2(TOT_LEN_INTFO_U)
      REAL :: O_COEFF3(NPTS_U_FIELD_O)
      REAL :: O_COEFF4(NPTS_U_FIELD_O)
      REAL :: O_COEFF5(TOT_LEN_INTFO_U)
      REAL :: O_COEFF6(TOT_LEN_INTFO_U)
! TYPINFO end
#endif
