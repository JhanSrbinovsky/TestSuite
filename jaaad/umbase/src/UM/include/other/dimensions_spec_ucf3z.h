!     ------------------------------------------------------------------
!     Module setting sizes of spectral arrays.
!
      INTEGER                                                           &
     &    NPD_BAND                                                      &
!           Number of spectral bands
     &  , NPD_EXCLUDE                                                   &
!           Numer of excluded bands
     &  , NPD_ESFT_TERM                                                 &
!           Number of esft terms
     &  , NPD_TYPE                                                      &
!           Number of data types
     &  , NPD_SPECIES                                                   &
!           Number of gaseous species
     &  , NPD_SCALE_VARIABLE                                            &
!           Number of scaling variables
     &  , NPD_SURFACE                                                   &
!           Number of surface types
     &  , NPD_BRDF_BASIS_FNC                                            &
!           Number of BRDF basis functions
     &  , NPD_BRDF_TRUNC                                                &
!           Order of BRDF truncation
     &  , NPD_ALBEDO_PARM                                               &
!           Number of albedo parameters
     &  , NPD_CONTINUUM                                                 &
!           Number of continua
     &  , NPD_DROP_TYPE                                                 &
!           Number of drop types
     &  , NPD_ICE_TYPE                                                  &
!           Number of ice crystal types
     &  , NPD_AEROSOL_SPECIES                                           &
!           Number of aerosol species
     &  , NPD_THERMAL_COEFF                                             &
!           Number of thermal coefficients
     &  , NPD_CLOUD_PARAMETER                                           &
!           Number of cloud parameters
     &  , NPD_HUMIDITIES                                                &
!           Number of humidities
     &  , NPD_PHASE_TERM                                                &
!           Number of terms in the phase function
     &  , NPD_CHANNEL
!           Number of spectral channels permitted for output
!
      PARAMETER(                                                        &
     &    NPD_BAND=20                                                   &
     &  , NPD_EXCLUDE=2                                                 &
     &  , NPD_ESFT_TERM=25                                              &
     &  , NPD_TYPE=20                                                   &
     &  , NPD_SPECIES=11                                                &
     &  , NPD_ALBEDO_PARM=4                                             &
     &  , NPD_SCALE_VARIABLE=4                                          &
     &  , NPD_SURFACE=10                                                &
     &  , NPD_BRDF_BASIS_FNC=2                                          &
     &  , NPD_BRDF_TRUNC=5                                              &
     &  , NPD_CONTINUUM=3                                               &
     &  , NPD_DROP_TYPE=6                                               &
     &  , NPD_ICE_TYPE=10                                               &
     &  , NPD_AEROSOL_SPECIES=11                                        &
     &  , NPD_THERMAL_COEFF=9                                           &
     &  , NPD_CLOUD_PARAMETER=36                                        &
     &  , NPD_HUMIDITIES=21                                             &
     &  , NPD_PHASE_TERM=101                                            &
     &  , NPD_CHANNEL=2                                                 &
     &  )
!
!     ------------------------------------------------------------------
