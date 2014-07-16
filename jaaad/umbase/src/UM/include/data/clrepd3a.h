! CLREPD3A defines representations of clouds in two-stream radiation
! code.

      DATA NP_CLOUD_TYPE(IP_CLOUD_HOMOGEN)/1/
      DATA NP_CLOUD_TYPE(IP_CLOUD_ICE_WATER)/2/
      DATA NP_CLOUD_TYPE(IP_CLOUD_CONV_STRAT)/2/
      DATA NP_CLOUD_TYPE(IP_CLOUD_CSIW)/4/

      ! the array ip_cloud_type_map indicates to which type of cloud
      ! each component belongs in a particular representation. an
      ! entry of 0 indicates that that component should not be
      ! present in the representation.

      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_HOMOGEN)        &
     &       /IP_CLOUD_TYPE_HOMOGEN/                                    &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_HOMOGEN)          &
     &       /IP_CLOUD_TYPE_HOMOGEN/                                    &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_HOMOGEN)       &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_HOMOGEN)         &
     &       /0/
      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_ICE_WATER)      &
     &       /IP_CLOUD_TYPE_WATER/                                      &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_ICE_WATER)        &
     &       /IP_CLOUD_TYPE_ICE/                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_ICE_WATER)     &
     &       /0/                                                        &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_ICE_WATER)       &
     &       /0/
      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_CONV_STRAT)     &
     &       /IP_CLOUD_TYPE_STRAT/                                      &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_CONV_STRAT)       &
     &       /IP_CLOUD_TYPE_STRAT/                                      &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_CONV_STRAT)    &
     &       /IP_CLOUD_TYPE_CONV/                                       &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_CONV_STRAT)      &
     &       /IP_CLOUD_TYPE_CONV/
      DATA                                                              &
     &    IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_WATER, IP_CLOUD_CSIW)           &
     &       /IP_CLOUD_TYPE_SW/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_ST_ICE, IP_CLOUD_CSIW)             &
     &       /IP_CLOUD_TYPE_SI/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_WATER, IP_CLOUD_CSIW)          &
     &       /IP_CLOUD_TYPE_CW/                                         &
     &  , IP_CLOUD_TYPE_MAP(IP_CLCMP_CNV_ICE, IP_CLOUD_CSIW)            &
     &       /IP_CLOUD_TYPE_CI/

! CLREPD3A end
