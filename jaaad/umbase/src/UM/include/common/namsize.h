!
! Description:
! NLSIZES namelist. Set up by UMUI and put in SIZES member of Job.
!
      NAMELIST/NLSIZES/                                                 &
     &  global_ROW_LENGTH,global_ROWS,LAND_FIELD,                       &
     &  MODEL_LEVELS,WET_LEVELS,                                        &
     &  NTILES,                                                         &
     &  CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS,           &
     &  OZONE_LEVELS,tpps_ozone_levels,TR_VARS,TR_UKCA,                 &
     &  RIVER_ROWS, RIVER_ROW_LENGTH,                                   &
     &  A_LEN_INTHD,A_LEN_REALHD,                                       &
     &  A_LEN2_LEVDEPC,A_LEN2_ROWDEPC,A_LEN2_COLDEPC,                   &
     &  A_LEN2_FLDDEPC,A_LEN_EXTCNST,                                   &
     &  A_LEN_CFI1,A_LEN_CFI2,A_LEN_CFI3,                               &

     &  KM,IMT,JMT,NT,NISLE,LSEG,ISEGM,N_STRAIT,N_STRAIT_CLM,           &
     &  NICE,O_LEN_COMPRESSED,                                          &

     &  NANCIL_LOOKUPSA, N_INTF_A,                                      &

     &  RIMWIDTHA, NRIM_TIMESA,                                         &

     &  PP_LEN_INTHD,PP_LEN_REALHD,                                     &

     & THETA_PV_P_LEVS, N_AOBS
! NAMSIZE end
