#if defined(ATMOS)
! ARGINFA Headers for atmosphere interface data sets

! 14/06/94     DEF LBOUTA replaced by ATMOS for vn 3.4   S.J.Swarbrick
! 29/07/98     Add INTF_AKH, INTF_BKH, INTF_AK and INTF_BK. D. Robinson.
! 10/11/00 5.2 Add LBC_ETA_THETA and LBC_ETA_RHO. D.Robinson
! 18/08/04 6.1 Add AV_INDEX_* and AV_WEIGHT_* D Robinson.
     &  FIXHD_INTFA, INTHD_INTFA, LOOKUP_INTFA,                         &
     &  REALHD_INTFA,LEVDEPC_INTFA,                                     &
      ! Interpolation constants for atmosphere interface data sets.
     &  AP_INDEX_B_L, AP_INDEX_B_R, AU_INDEX_B_L, AU_INDEX_B_R,         &
     &  AP_WEIGHT_T_R, AP_WEIGHT_B_L, AP_WEIGHT_B_R, AP_WEIGHT_T_L,     &
     &  AU_WEIGHT_T_R, AU_WEIGHT_B_L, AU_WEIGHT_B_R, AU_WEIGHT_T_L,     &
      ! Rotation coefficients for atmosphere interface data sets
     &  COEFF1, COEFF2, COEFF3, COEFF4,                                 &
     &  INTF_AKH, INTF_BKH, INTF_AK, INTF_BK,                           &
      ! Eta values for LBC levels
     &  LBC_ETA_THETA, LBC_ETA_RHO,                                     &
     &  AV_INDEX_B_L,  AV_INDEX_B_R,                                    &
     &  AV_WEIGHT_T_R, AV_WEIGHT_B_L, AV_WEIGHT_B_R, AV_WEIGHT_T_L,     &
      ! Row/Col DEPC for variable resolution LBCs
     &  ROWDEPC_INTFA, COLDEPC_INTFA,                                   &  
! ARGINFA end
#endif
