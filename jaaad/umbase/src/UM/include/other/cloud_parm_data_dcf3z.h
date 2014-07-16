!     ------------------------------------------------------------------
!     Module to set numbers of parameters for cloud schemes
!
      DATA                                                              &
     &    ( N_CLOUD_FIT_PARAMETER(IP_SLINGO_SCHRECKER, INLIST)          &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/2, 2, 2, 0, 0, 0/       &
     &  , ( N_CLOUD_FIT_PARAMETER(IP_ACKERMAN_STEPHENS, INLIST)         &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/3, 3, 3, 0, 0, 0/       &
     &  , ( N_CLOUD_FIT_PARAMETER(IP_DROP_PARAMETRIZATION_TEST, INLIST) &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/6, 5, 5, 0, 0, 0/       &
     &  , ( N_CLOUD_FIT_PARAMETER(IP_DROP_PADE_2, INLIST)               &
     &  , INLIST=1, NPD_CLOUD_OPTICAL_VARIABLE)/6, 5, 5, 0, 0, 0/
      DATA                                                              &
     &     NAME_CLOUD_PARAMETRIZATION(IP_SLINGO_SCHRECKER)              &
     &   / 'SLINGO-SCHRECKER    '/                                      &
     &  , NAME_CLOUD_PARAMETRIZATION(IP_ACKERMAN_STEPHENS)              &
     &   / 'ACKERMAN-STEPHENS   '/                                      &
     &  , NAME_CLOUD_PARAMETRIZATION(IP_DROP_UNPARAMETRIZED)            &
     &   / 'UNPARAMETRIZED DATA '/                                      &
     &  , NAME_CLOUD_PARAMETRIZATION(IP_DROP_PARAMETRIZATION_TEST)      &
     &   / 'TEST PARAMETRIZATION'/                                      &
     &  , NAME_CLOUD_PARAMETRIZATION(IP_DROP_PADE_2)                    &
     &   / '2ND ORDER PADE APPX '/
!
!     ------------------------------------------------------------------
