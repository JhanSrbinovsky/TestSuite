!     ------------------------------------------------------------------
!     Module to set numbers of parameters for ice cloud schemes
!
      DATA                                                              &
     &    ( N_ICE_CLOUD_FIT_PARAMETER(IP_SLINGO_SCHRECKER_ICE, INLIST)  &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/2, 2, 2, 0, 0, 0/       &
     &  , ( N_ICE_CLOUD_FIT_PARAMETER(IP_SUN_SHINE_VN2_VIS, INLIST)     &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/2, 2, 2, 0, 0, 0/       &
     &  , ( N_ICE_CLOUD_FIT_PARAMETER(IP_SUN_SHINE_VN2_IR, INLIST)      &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/2, 0, 0, 0, 0, 0/       &
     &  , ( N_ICE_CLOUD_FIT_PARAMETER(IP_ICE_ADT, INLIST)               &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/10, 10, 10, 0, 0, 0/    &
     &  , ( N_ICE_CLOUD_FIT_PARAMETER(IP_ICE_ADT_10, INLIST)            &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/12, 12, 12, 0, 0, 0/    &
     &  , ( N_ICE_CLOUD_FIT_PARAMETER(IP_ICE_PARAMETRIZATION_TEST       &
     &    , INLIST)                                                     &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/12, 12, 12, 0, 0, 0/    &
     &  , ( N_ICE_CLOUD_FIT_PARAMETER(IP_ICE_FU_SOLAR, INLIST)          &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/2, 3, 4, 0, 0, 4/       &
     &  , ( N_ICE_CLOUD_FIT_PARAMETER(IP_ICE_FU_IR, INLIST)             &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/3, 0, 4, 0, 4, 0/
      DATA                                                              &
     &    NAME_ICE_PARAMETRIZATION(IP_SLINGO_SCHRECKER_ICE)             &
     &      /'SLINGO-SCHRECKER:ICE'/                                    &
     &  , NAME_ICE_PARAMETRIZATION(IP_ICE_UNPARAMETRIZED)               &
     &      /'UNPARAMETRIZED DATA.'/                                    &
     &  , NAME_ICE_PARAMETRIZATION(IP_SUN_SHINE_VN2_VIS)                &
     &      /'SUN-SHINE (VIS. VN2)'/                                    &
     &  , NAME_ICE_PARAMETRIZATION(IP_SUN_SHINE_VN2_IR)                 &
     &      /'SUN-SHINE (IR VN2)  '/                                    &
     &  , NAME_ICE_PARAMETRIZATION(IP_ICE_ADT)                          &
     &      /'ADT FOR ICE CRYSTALS'/                                    &
     &  , NAME_ICE_PARAMETRIZATION(IP_ICE_ADT_10)                       &
     &      /'10TH ORDER ICE ADT  '/                                    &
     &  , NAME_ICE_PARAMETRIZATION(IP_ICE_PARAMETRIZATION_TEST)         &
     &      /'TEST SCHEME FOR ICE '/                                    &
     &  , NAME_ICE_PARAMETRIZATION(IP_ICE_FU_SOLAR)                     &
     &      /'Solar scheme of Fu  '/                                    &
     &  , NAME_ICE_PARAMETRIZATION(IP_ICE_FU_IR)                        &
     &      /'IR scheme of Fu     '/
!
!     ------------------------------------------------------------------
