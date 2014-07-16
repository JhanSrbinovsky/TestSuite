!     ------------------------------------------------------------------
!     Module to set options for weightings.
!
      INTEGER                                                           &
     &    IP_WEIGHT_PLANCK                                              &
!           Planckian weighting
     &  , IP_WEIGHT_DPLANCK                                             &
!           Differential planckian weighting
     &  , IP_WEIGHT_SOLAR                                               &
!           Solar weighting
     &  , IP_WEIGHT_UNIFORM
!           Uniform weighting across the band
!
      PARAMETER(                                                        &
     &    IP_WEIGHT_PLANCK=1                                            &
     &  , IP_WEIGHT_DPLANCK=2                                           &
     &  , IP_WEIGHT_SOLAR=3                                             &
     &  , IP_WEIGHT_UNIFORM=4                                           &
     &  )
!
!     ------------------------------------------------------------------
