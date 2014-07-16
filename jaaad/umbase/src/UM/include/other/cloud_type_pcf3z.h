!     ------------------------------------------------------------------
!     Module to set types of clouds.
!
      INTEGER                                                           &
     &    IP_CLOUD_TYPE_HOMOGEN                                         &
!           Cloud composed of mixed water and ice
     &  , IP_CLOUD_TYPE_WATER                                           &
!           Cloud composed only of water
     &  , IP_CLOUD_TYPE_ICE                                             &
!           Cloud composed only of ice
     &  , IP_CLOUD_TYPE_STRAT                                           &
!           Mixed-phase stratiform cloud
     &  , IP_CLOUD_TYPE_CONV                                            &
!           Mixed-phase convective cloud
     &  , IP_CLOUD_TYPE_SW                                              &
!           Stratiform water cloud
     &  , IP_CLOUD_TYPE_SI                                              &
!           Stratiform ice cloud
     &  , IP_CLOUD_TYPE_CW                                              &
!           Convective water cloud
     &  , IP_CLOUD_TYPE_CI
!           Convective ice cloud
!
      PARAMETER(                                                        &
     &    IP_CLOUD_TYPE_HOMOGEN=1                                       &
     &  , IP_CLOUD_TYPE_WATER=1                                         &
     &  , IP_CLOUD_TYPE_ICE=2                                           &
     &  , IP_CLOUD_TYPE_STRAT=1                                         &
     &  , IP_CLOUD_TYPE_CONV=2                                          &
     &  , IP_CLOUD_TYPE_SW=1                                            &
     &  , IP_CLOUD_TYPE_SI=2                                            &
     &  , IP_CLOUD_TYPE_CW=3                                            &
     &  , IP_CLOUD_TYPE_CI=4                                            &
     &  )
!
!     ------------------------------------------------------------------
