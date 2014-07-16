!     ------------------------------------------------------------------
!     Module to set identification of distribution types.
!
      INTEGER                                                           &
     &    NPD_DISTRIBUTION
!           Number of supported distributions
      PARAMETER(                                                        &
     &    NPD_DISTRIBUTION=4                                            &
     &  )
!
      INTEGER                                                           &
     &    IP_EXTERNAL                                                   &
!           External file type
     &  , IP_LOG_NORMAL                                                 &
!           Log normal distribution
     &  , IP_MODIFIED_GAMMA                                             &
!           Modified gamma distribution
     &  , IP_HEYMSFIELD_PLATT                                           &
!           Heymsfield and Platt's distribution
     &  , IP_MITCHELL_97                                                &
!           Mitchell's size distribution
     &  , IP_MITCHELL_TROP_00                                           &
!           Mitchell's distribution for tropical Ci
     &  , IP_IVANOVA_MLAT_00
!           Ivanova's distribution for mid-latitude Ci
!
      PARAMETER(                                                        &
     &    IP_EXTERNAL=1                                                 &
     &  , IP_LOG_NORMAL=2                                               &
     &  , IP_MODIFIED_GAMMA=3                                           &
     &  , IP_HEYMSFIELD_PLATT=4                                         &
     &  , IP_MITCHELL_97=5                                              &
     &  , IP_MITCHELL_TROP_00=6                                         &
     &  , IP_IVANOVA_MLAT_00=7                                          &
     &  )
!
!     ------------------------------------------------------------------
