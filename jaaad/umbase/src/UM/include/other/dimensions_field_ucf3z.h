!     ------------------------------------------------------------------
!     Module setting the dimensions of physical fields
!     ------------------------------------------------------------------
!
      INTEGER                                                           &
     &    NPD_LATITUDE                                                  &
!           Number of latitudes
     &  , NPD_LONGITUDE                                                 &
!           Number of longitudes
     &  , NPD_PROFILE                                                   &
!           Number of atmospheric profiles
     &  , NPD_LAYER                                                     &
!           Number of atmospheric layers
     &  , NPD_COLUMN                                                    &
!           Maximum number of cloudy subcolumns
     &  , NPD_DIRECTION                                                 &
!           Maximum number of directions for radiances
     &  , NPD_MAX_ORDER                                                 &
!           Maximum order of spherical harmonics used
#if defined(OBSERVED)
     &  , NPD_PROFILE_AEROSOL_PRSC                                      &
!           Size allocated for profiles of prescribed
!           cloudy optical properties
     &  , NPD_PROFILE_CLOUD_PRSC                                        &
!           Size allocated for profiles of prescribed
!           aerosol optical properties
     &  , NPD_OPT_LEVEL_AEROSOL_PRSC                                    &
!           Size allocated for levels of prescribed
!           cloudy optical properties
     &  , NPD_OPT_LEVEL_CLOUD_PRSC
!           Size allocated for levels of prescribed
!           aerosol optical properties
#endif
!
      PARAMETER(                                                        &
     &    NPD_LONGITUDE=1                                               &
     &  , NPD_LATITUDE=3                                                &
     &  , NPD_PROFILE=3                                                 &
     &  , NPD_LAYER=61                                                  &
     &  , NPD_COLUMN=2                                                  &
     &  , NPD_DIRECTION=11                                              &
     &  , NPD_MAX_ORDER=101                                             &
#if defined(OBSERVED)
     &  , NPD_PROFILE_AEROSOL_PRSC=1                                    &
     &  , NPD_PROFILE_CLOUD_PRSC=1                                      &
     &  , NPD_OPT_LEVEL_AEROSOL_PRSC=10                                 &
     &  , NPD_OPT_LEVEL_CLOUD_PRSC=2                                    &
#endif
     &  )
!
!     ------------------------------------------------------------------
