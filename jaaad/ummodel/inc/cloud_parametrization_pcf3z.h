!     ------------------------------------------------------------------
!     Module to set numbers for water cloud schemes.
!
      INTEGER                                                           &
     &    NPD_CLOUD_FIT                                                 &
!           Number of cloud fitting schemes
     &  , IP_SLINGO_SCHRECKER                                           &
!           Parametrization of Slingo & Schrecker
     &  , IP_ACKERMAN_STEPHENS                                          &
!           Parametrization of Ackerman & Stephens
#if ! defined(UM)
     &  , IP_DROP_UNPARAMETRIZED                                        &
!           Unparametrized droplet data
     &  , IP_DROP_PARAMETRIZATION_TEST                                  &
!           Test parametrization
#endif
     &  , IP_DROP_PADE_2                                                &
!           Pade approximation of the second order
!           (third order for the extinction)
     &  , IP_SLINGO_SCHR_PHF
!           Parameterization of Slingo & Schrecker +
!           Moments of Phase Functions
!
      PARAMETER(                                                        &
     &    NPD_CLOUD_FIT=5                                               &
     &  , IP_SLINGO_SCHRECKER=1                                         &
     &  , IP_ACKERMAN_STEPHENS=2                                        &
#if ! defined(UM)
     &  , IP_DROP_UNPARAMETRIZED=3                                      &
     &  , IP_DROP_PARAMETRIZATION_TEST=4                                &
#endif
     &  , IP_DROP_PADE_2=5                                              &
     &  , IP_SLINGO_SCHR_PHF=6                                          &
     &  )
!
!     ------------------------------------------------------------------
