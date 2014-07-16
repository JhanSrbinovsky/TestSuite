!     ------------------------------------------------------------------
!     Module to set properties of representations of clouds.
!
      DATA                                                              &
     &    NP_CLOUD_TYPE(IP_CLOUD_HOMOGEN)/1/                            &
     &  , NP_CLOUD_TYPE(IP_CLOUD_ICE_WATER)/2/                          &
     &  , NP_CLOUD_TYPE(IP_CLOUD_CONV_STRAT)/2/                         &
     &  , NP_CLOUD_TYPE(IP_CLOUD_CSIW)/4/                               &
#if defined(XICE)
     &  , NP_CLOUD_TYPE(IP_CLOUD_CSIW_CRYS)/4/
#endif
!
!
!     The array IP_CLOUD_TYPE_MAP indicates to which type of cloud
!     each component belongs in a particular representation. An
!     entry of 0 indicates that that component should not be
!     present in the representation.
!
      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_HOMOGEN)        &
     &       /IP_CLOUD_TYPE_HOMOGEN/                                    &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_HOMOGEN)          &
     &       /IP_CLOUD_TYPE_HOMOGEN/                                    &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_HOMOGEN)       &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_HOMOGEN)         &
     &       /0/                                                        &
#if defined(XICE)
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_COL, IP_CLOUD_HOMOGEN)      &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_COL, IP_CLOUD_HOMOGEN)     &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_ROS, IP_CLOUD_HOMOGEN)      &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_ROS, IP_CLOUD_HOMOGEN)     &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_PLT, IP_CLOUD_HOMOGEN)      &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_PLT, IP_CLOUD_HOMOGEN)     &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_PYC, IP_CLOUD_HOMOGEN)      &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_PYC, IP_CLOUD_HOMOGEN)     &
     &       /0/
#endif
      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_ICE_WATER)      &
     &       /IP_CLOUD_TYPE_WATER/                                      &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_ICE_WATER)        &
     &       /IP_CLOUD_TYPE_ICE/                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_ICE_WATER)     &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_ICE_WATER)       &
     &       /0/                                                        &
#if defined(XICE)
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_COL, IP_CLOUD_ICE_WATER)    &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_COL, IP_CLOUD_ICE_WATER)   &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_ROS, IP_CLOUD_ICE_WATER)    &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_ROS, IP_CLOUD_ICE_WATER)   &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_PLT, IP_CLOUD_ICE_WATER)    &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_PLT, IP_CLOUD_ICE_WATER)   &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_PYC, IP_CLOUD_ICE_WATER)    &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_PYC, IP_CLOUD_ICE_WATER)   &
     &       /0/
#endif
      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_CONV_STRAT)     &
     &       /IP_CLOUD_TYPE_STRAT/                                      &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_CONV_STRAT)       &
     &       /IP_CLOUD_TYPE_CONV/                                       &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_CONV_STRAT)    &
     &       /IP_CLOUD_TYPE_STRAT/                                      &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_CONV_STRAT)      &
     &       /IP_CLOUD_TYPE_CONV/                                       &
#if defined(XICE)
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_COL, IP_CLOUD_CONV_STRAT)   &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_COL, IP_CLOUD_CONV_STRAT)  &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_ROS, IP_CLOUD_CONV_STRAT)   &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_ROS, IP_CLOUD_CONV_STRAT)  &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_PLT, IP_CLOUD_CONV_STRAT)   &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_PLT, IP_CLOUD_CONV_STRAT)  &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_PYC, IP_CLOUD_CONV_STRAT)   &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_PYC, IP_CLOUD_CONV_STRAT)  &
     &       /0/
#endif
      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_CSIW)           &
     &       /IP_CLOUD_TYPE_SW/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_CSIW)             &
     &       /IP_CLOUD_TYPE_SI/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_CSIW)          &
     &       /IP_CLOUD_TYPE_CW/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_CSIW)            &
     &       /IP_CLOUD_TYPE_CI/                                         &
#if defined(XICE)
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_COL, IP_CLOUD_CSIW)         &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_COL, IP_CLOUD_CSIW)        &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_ROS, IP_CLOUD_CSIW)         &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_ROS, IP_CLOUD_CSIW)        &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_PLT, IP_CLOUD_CSIW)         &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_PLT, IP_CLOUD_CSIW)        &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_PYC, IP_CLOUD_CSIW)         &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_PYC, IP_CLOUD_CSIW)        &
     &       /0/
#endif
#if defined(XICE)
      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_CSIW_CRYS)      &
     &       /IP_CLOUD_TYPE_SW/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_CSIW_CRYS)        &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_CSIW_CRYS)     &
     &       /IP_CLOUD_TYPE_CW/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_CSIW_CRYS)       &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_COL, IP_CLOUD_CSIW_CRYS)    &
     &       /IP_CLOUD_TYPE_SI/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_COL, IP_CLOUD_CSIW_CRYS)   &
     &       /IP_CLOUD_TYPE_CI/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_ROS, IP_CLOUD_CSIW_CRYS)    &
     &       /IP_CLOUD_TYPE_SI/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_ROS, IP_CLOUD_CSIW_CRYS)   &
     &       /IP_CLOUD_TYPE_CI/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_PLT, IP_CLOUD_CSIW_CRYS)    &
     &       /IP_CLOUD_TYPE_SI/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_PLT, IP_CLOUD_CSIW_CRYS)   &
     &       /IP_CLOUD_TYPE_CI/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE_PYC, IP_CLOUD_CSIW_CRYS)    &
     &       /IP_CLOUD_TYPE_SI/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE_PYC, IP_CLOUD_CSIW_CRYS)   &
     &       /IP_CLOUD_TYPE_CI/
#endif
!
!     ------------------------------------------------------------------
