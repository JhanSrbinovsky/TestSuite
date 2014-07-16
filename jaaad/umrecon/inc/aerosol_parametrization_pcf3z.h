!     ------------------------------------------------------------------
!     Module to set the parametrizations available for aerosols.
!
      INTEGER                                                           &
     &    IP_AEROSOL_PARAM_DRY                                          &
!           Parametrization for dry aerosols
     &  , IP_AEROSOL_PARAM_MOIST                                        &
!           Parametrization for moist aerosols
     &  , IP_AEROSOL_UNPARAMETRIZED                                     &
!          Observational aerosol data
     &  , IP_AEROSOL_PARAM_PHF_DRY                                      &
!           Parametrization of the phase function for dry aerosols
     &  , IP_AEROSOL_PARAM_PHF_MOIST
!           Parametrization of the phase function for moist aerosols
!
      PARAMETER(                                                        &
     &    IP_AEROSOL_PARAM_DRY=1                                        &
     &  , IP_AEROSOL_PARAM_MOIST=2                                      &
     &  , IP_AEROSOL_UNPARAMETRIZED=3                                   &
     &  , IP_AEROSOL_PARAM_PHF_DRY=4                                    &
     &  , IP_AEROSOL_PARAM_PHF_MOIST=5                                  &
     &  )
!
!     ------------------------------------------------------------------
