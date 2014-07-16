!     ------------------------------------------------------------------
!     Module to set types of scaling for absorber amounts
!
      INTEGER                                                           &
     &    IP_SCALE_FNC_NULL                                             &
!           Null scaling function
     &  , IP_SCALE_POWER_LAW                                            &
!           Power law scaling function
     &  , IP_SCALE_POWER_QUAD                                           &
!           Power law for p; quadratic for T
     &  , IP_SCALE_DOPPLER_QUAD
!           Power law for p; quadratic for T with implicit
!           Doppler correction
!
      PARAMETER(                                                        &
     &    IP_SCALE_FNC_NULL=0                                           &
     &  , IP_SCALE_POWER_LAW=1                                          &
     &  , IP_SCALE_POWER_QUAD=2                                         &
     &  , IP_SCALE_DOPPLER_QUAD=3                                       &
     &  )
!
!     -----------------------------------------------------------------
