!     ------------------------------------------------------------------
!     Module to set indices of aerosol components.
!
!
      INTEGER                                                           &
     &    NPD_AEROSOL_COMPONENT
!           Maximum number of aerosol components
      INTEGER                                                           &
     &    IP_WATER_SOLUBLE                                              &
!           Water soluble aerosol
     &  , IP_DUST_LIKE                                                  &
!           Dust-like aerosol
     &  , IP_OCEANIC                                                    &
!           Oceanic aerosol
     &  , IP_SOOT                                                       &
!           Soot aerosol
     &  , IP_ASH                                                        &
!           Volcanic ash
     &  , IP_SULPHURIC                                                  &
!           Sulphuric acid
     &  , IP_AMMONIUM_SULPHATE                                          &
!           Ammonium sulphate aerosol
     &  , IP_AEROSOL_UNCHARACTERIZED                                    &
!           Uncharacterized aerosol (for observations)
     &  , IP_SAHARAN_DUST                                               &
!           Saharan dust
     &  , IP_ACCUM_SULPHATE                                             &
!           Accumulation mode sulphate
     &  , IP_AITKEN_SULPHATE                                            &
!           Aitken mode sulphate
     &  , IP_FRESH_SOOT                                                 &
!           Fresh soot
     &  , IP_AGED_SOOT                                                  &
!           Aged soot
     &  , IP_SODIUM_CHLORIDE                                            &
!           Sodium chloride (generic aerosol)
     &  , IP_SEASALT_FILM                                               &
!           Sodium chloride (film mode)
     &  , IP_SEASALT_JET                                                &
!           Sodium chloride (jet mode)
     &  , IP_DUST_1                                                  &
                                        !steph
!           Dust, ision1                  !steph
     &  , IP_DUST_2                                                  &
                                        !steph
!           Dust, ision2                  !steph
     &  , IP_DUST_3                                                  &
                                        !steph
!           Dust, ision3                  !steph
     &  , IP_DUST_4                                                  &
                                        !steph
!           Dust, ision4                  !steph
     &  , IP_DUST_5                                                  &
                                        !steph
!           Dust, ision5                  !steph
     &  , IP_DUST_6                                                  &
                                        !steph
!           Dust, ision6                  !steph
     &  , IP_BIOMASS_1                                                  &
!           Biomass (ision 1)
     &  , IP_BIOMASS_2                                                  &
!           Biomass (ision 2)
     &  , IP_BIOGENIC                                                   &
!           Biogenic
     &  , IP_OCFF_FRESH                                                 &
!           Organic Carbon Fossil Fuel (fresh)
     &  , IP_OCFF_AGED                                                  &
!           Organic Carbon Fossil Fuel (aged)
     &  , IP_DELTA
!           Delta aerosol
!
      PARAMETER(                                                        &
     &    NPD_AEROSOL_COMPONENT=28                                      &
     &   )
!
      PARAMETER(                                                        &
     &    IP_WATER_SOLUBLE=1                                            &
     &  , IP_DUST_LIKE=2                                                &
     &  , IP_OCEANIC=3                                                  &
     &  , IP_SOOT=4                                                     &
     &  , IP_ASH=5                                                      &
     &  , IP_SULPHURIC=6                                                &
     &  , IP_AMMONIUM_SULPHATE=7                                        &
     &  , IP_AEROSOL_UNCHARACTERIZED=8                                  &
     &  , IP_SAHARAN_DUST=9                                             &
     &  , IP_ACCUM_SULPHATE=10                                          &
     &  , IP_AITKEN_SULPHATE=11                                         &
     &  , IP_FRESH_SOOT=12                                              &
     &  , IP_AGED_SOOT=13                                               &
     &  , IP_SODIUM_CHLORIDE=14                                         &
     &  , IP_SEASALT_FILM=15                                            &
     &  , IP_SEASALT_JET=16                                             &
     &  , IP_DUST_1=17                                               &
     &  , IP_DUST_2=18                                               &
     &  , IP_DUST_3=19                                               &
     &  , IP_DUST_4=20                                               &
     &  , IP_DUST_5=21                                               &
     &  , IP_DUST_6=22                                               &
     &  , IP_BIOMASS_1=23                                               &
     &  , IP_BIOMASS_2=24                                               &
     &  , IP_BIOGENIC=25                                                &
     &  , IP_OCFF_FRESH=26                                              &
     &  , IP_OCFF_AGED=27                                               &
     &  , IP_DELTA=28                                                   &
     &  )
!
!     ------------------------------------------------------------------
