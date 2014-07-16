!     ------------------------------------------------------------------
!     Module to set methods of weighting optical properties.
!
      INTEGER                                                           &
     &    IP_WEIGHT_THIN                                                &
!           Weighting in thin limit
     &  , IP_WEIGHT_THICK                                               &
!           Weighting in thick limit
     &  , IP_WEIGHT_THIN_THICK
!           Combined thin and thick weighting
!
      PARAMETER(                                                        &
     &    IP_WEIGHT_THIN=1                                              &
     &  , IP_WEIGHT_THICK=2                                             &
     &  , IP_WEIGHT_THIN_THICK=3                                        &
     &  )
!
!     ------------------------------------------------------------------
