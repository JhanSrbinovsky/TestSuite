!     ------------------------------------------------------------------
!     Module to set the indexing numbers for the optical variables.
!
      INTEGER                                                           &
     &    NPD_CLOUD_OPTICAL_VARIABLE                                    &
!           Number of optical variables
     &  , IP_EXTINCTION                                                 &
!           Index number for extinction
     &  , IP_COALBEDO                                                   &
!           Index number for coalbedo
     &  , IP_ASYMMETRY                                                  &
!           Index number for asymmetry
     &  , IP_SCATTERING                                                 &
!           Index number for scattering
     &  , IP_ABSORPTION                                                 &
!           Index number for absorption
     &  , IP_FORWARD_FU
!           Index number for Fu's f_delta contribution to forward
!           scattering
!
      PARAMETER(                                                        &
     &    NPD_CLOUD_OPTICAL_VARIABLE=6                                  &
     &  , IP_EXTINCTION=1                                               &
     &  , IP_COALBEDO=2                                                 &
     &  , IP_ASYMMETRY=3                                                &
     &  , IP_SCATTERING=4                                               &
     &  , IP_ABSORPTION=5                                               &
     &  , IP_FORWARD_FU=6                                               &
     &  )
!
!     ------------------------------------------------------------------
