!     ------------------------------------------------------------------
!     Module to set the defining numbers for the two-stream
!     schemes.
!
      INTEGER                                                           &
     &    IP_EDDINGTON                                                  &
!           Eddington approximation
     &  , IP_DISCRETE_ORD                                               &
!           Discrete ordinate method
     &  , IP_IFM                                                        &
!           Improved flux method
     &  , IP_PIFM85                                                     &
!           Practical improved flux method
!           (version of zdunkowski et al. 1985)
     &  , IP_ZDK_FLUX                                                   &
!           Zdunkowski's flux method
     &  , IP_KRSCHG_FLUX                                                &
!           Kerschgen's flux method
     &  , IP_COAKLEY_CHYLEK_1                                           &
!           Coakley & chylek's 1st method
     &  , IP_COAKLEY_CHYLEK_2                                           &
!           Coakley & chylek's 2nd method
     &  , IP_MEADOR_WEAVER                                              &
!           Meador & weaver's method
     &  , IP_ELSASSER                                                   &
!           Elsasser's diffusivity scheme
     &  , IP_2S_TEST                                                    &
!           User's defined test approximation.
     &  , IP_HEMI_MEAN                                                  &
!           Hemispheric mean approximation.
     &  , IP_PIFM80
!           Hemispheric mean approximation.
!             (version of zdunkowski et al. 1980)
!
      PARAMETER(                                                        &
     &    IP_EDDINGTON=2                                                &
     &  , IP_DISCRETE_ORD=4                                             &
     &  , IP_IFM=5                                                      &
     &  , IP_PIFM85=6                                                   &
     &  , IP_ZDK_FLUX=7                                                 &
     &  , IP_KRSCHG_FLUX=8                                              &
     &  , IP_COAKLEY_CHYLEK_1=9                                         &
     &  , IP_COAKLEY_CHYLEK_2=10                                        &
     &  , IP_MEADOR_WEAVER=11                                           &
     &  , IP_ELSASSER=12                                                &
     &  , IP_2S_TEST=14                                                 &
     &  , IP_HEMI_MEAN=15                                               &
     &  , IP_PIFM80=16                                                  &
     &  )
!
!     ------------------------------------------------------------------
