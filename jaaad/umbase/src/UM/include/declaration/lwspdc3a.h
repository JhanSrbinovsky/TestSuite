!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE CONTAINING DECLARATIONS FOR REDUCED LW-SPECTRAL FILE.
!     (NOTE: LWSPDC3A, LWSPCM3A AND LWSARG3A MUST BE CONSISTENT)
!     ------------------------------------------------------------------
!
!
!     DIMENSIONS FOR THE REDUCED LW SPECTRAL FILE
!
      INTEGER                                                           &
     &     NPD_TYPE_LW                                                  &
!             NUMBER OF TYPES OF DATA IN LW SPECTRUM
     &   , NPD_BAND_LW                                                  &
!             NUMBER OF SPECTRAL BANDS IN LW SPECTRUM
     &   , NPD_EXCLUDE_LW                                               &
!             NUMBER OF EXCLUDED BANDS IN LW SPECTRUM
     &   , NPD_SPECIES_LW                                               &
!             NUMBER OF GASEOUS SPECIES IN LW SPECTRUM
     &   , NPD_ESFT_TERM_LW                                             &
!             NUMBER OF ESFT TERMS IN LW SPECTRUM
     &   , NPD_SCALE_FNC_LW                                             &
!             NUMBER OF SCALING FUNCTIONS IN LW SPECTRUM
     &   , NPD_SCALE_VARIABLE_LW                                        &
!             NUMBER OF SCALING VARIABLES IN LW SPECTRUM
     &   , NPD_SURFACE_LW                                               &
!             NUMBER OF SURFACE TYPES IN LW SPECTRUM
     &   , NPD_ALBEDO_PARM_LW                                           &
!             NUMBER OF ALBEDO PARAMETERS IN LW SPECTRUM
     &   , NPD_CONTINUUM_LW                                             &
!             NUMBER OF CONTINUA IN LW SPECTRUM
     &   , NPD_DROP_TYPE_LW                                             &
!             NUMBER OF DROP TYPES IN LW SPECTRUM
     &   , NPD_ICE_TYPE_LW                                              &
!             NUMBER OF ICE CRYSTAL TYPES IN LW SPECTRUM
     &   , NPD_AEROSOL_SPECIES_LW                                       &
!             NUMBER OF AEROSOL SPECIES IN LW SPECTRUM
     &   , NPD_CLOUD_PARAMETER_LW                                       &
!             MAX NUMBER OF CLOUD PARAMETERS IN LW SPECTRUM
     &   , NPD_HUMIDITIES_LW                                            &
!             MAXIMUM NUMBER OF HUMIDITIES IN LW SPECTRUM
     &   , NPD_THERMAL_COEFF_LW                                         &
!             NUMBER OF THERMAL COEFFICIENTS IN LW SPECTRUM
     &   , NPD_AOD_WAVEL_LW
!             NUMBER OF WAVELENGTHS FOR AOD IN LW SPECTRUM
!
!
!     GENERAL FIELDS:
!
      LOGICAL                                                           &
     &     L_PRESENT_LW(0: NPD_TYPE_LW)
!             FLAG FOR TYPES OF DATA PRESENT
!
!
!
!     PROPERTIES OF THE SPECTRAL BANDS:
!
      INTEGER                                                           &
     &     N_BAND_LW
!             NUMBER OF SPECTRAL BANDS
!
      REAL                                                              &
     &     WAVE_LENGTH_SHORT_LW(NPD_BAND_LW)                            &
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG_LW(NPD_BAND_LW)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
!
      INTEGER                                                           &
     &     N_BAND_EXCLUDE_LW(NPD_BAND_LW)                               &
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE_LW(NPD_EXCLUDE_LW, NPD_BAND_LW)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL                                                              &
     &     SOLAR_FLUX_BAND_LW(NPD_BAND_LW)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL                                                              &
     &     RAYLEIGH_COEFFICIENT_LW(NPD_BAND_LW)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER                                                           &
     &     N_ABSORB_LW                                                  &
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB_LW(NPD_BAND_LW)                                &
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB_LW(NPD_SPECIES_LW, NPD_BAND_LW)                 &
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB_LW(NPD_SPECIES_LW)                               &
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT_LW(NPD_BAND_LW, NPD_SPECIES_LW)                  &
!             NUMBER OF ESFT TERMS IN EACH BAND FOR EACH GAS
     &   , I_SCALE_ESFT_LW(NPD_BAND_LW, NPD_SPECIES_LW)                 &
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC_LW(NPD_BAND_LW, NPD_SPECIES_LW)
!             TYPE OF SCALING FUNCTION
!
      REAL                                                              &
     &     K_ESFT_LW(NPD_ESFT_TERM_LW, NPD_BAND_LW, NPD_SPECIES_LW)     &
!             ESFT EXPONENTS
     &   , W_ESFT_LW(NPD_ESFT_TERM_LW, NPD_BAND_LW, NPD_SPECIES_LW)     &
!             ESFT WEIGHTS
     &   , SCALE_VECTOR_LW(NPD_SCALE_VARIABLE_LW, NPD_ESFT_TERM_LW      &
     &        , NPD_BAND_LW, NPD_SPECIES_LW)                            &
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE_LW(NPD_SPECIES_LW, NPD_BAND_LW)                  &
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE_LW(NPD_SPECIES_LW, NPD_BAND_LW)
!             REFERENCE TEMPERATURE FOR SCALING FUNCTION
!
!
!
!     REPRESENTATION OF THE PLANCKIAN:
!
      INTEGER                                                           &
     &     N_DEG_FIT_LW
!             DEGREE OF THERMAL POLYNOMIAL
!
      REAL                                                              &
     &     THERMAL_COEFFICIENT_LW(0: NPD_THERMAL_COEFF_LW-1             &
     &   , NPD_BAND_LW)                                                 &
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK_LW
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER                                                           &
     &     I_SPEC_SURFACE_LW(NPD_SURFACE_LW)                            &
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT_LW(NPD_SURFACE_LW)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL                                                           &
     &     L_SURFACE_LW(NPD_SURFACE_LW)
!             SURFACE TYPES INCLUDED
!
      REAL                                                              &
     &     SURFACE_ALBEDO_LW(NPD_BAND_LW, NPD_SURFACE_LW)               &
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM_LW(0: NPD_ALBEDO_PARM_LW                  &
     &      , NPD_BAND_LW, NPD_SURFACE_LW)                              &
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND_LW(NPD_BAND_LW, NPD_SURFACE_LW)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER                                                           &
     &     N_BAND_CONTINUUM_LW(NPD_BAND_LW)                             &
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM_LW(NPD_BAND_LW, NPD_CONTINUUM_LW)            &
!             LIST OF CONTINUA IN EACH BAND
     &   , INDEX_WATER_LW                                               &
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT_LW(NPD_BAND_LW, NPD_CONTINUUM_LW)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL                                                              &
     &     K_CONTINUUM_LW(NPD_BAND_LW, NPD_CONTINUUM_LW)                &
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM_LW(NPD_SCALE_VARIABLE_LW                     &
     &      , NPD_BAND_LW, NPD_CONTINUUM_LW)                            &
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM_LW(NPD_CONTINUUM_LW, NPD_BAND_LW)            &
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM_LW(NPD_CONTINUUM_LW, NPD_BAND_LW)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER                                                           &
     &     I_DROP_PARAMETRIZATION_LW(NPD_DROP_TYPE_LW)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL                                                           &
     &     L_DROP_TYPE_LW(NPD_DROP_TYPE_LW)
!             TYPES OF DROPLET PRESENT
!
      REAL                                                              &
     &     DROP_PARAMETER_LIST_LW(NPD_CLOUD_PARAMETER_LW                &
     &        , NPD_BAND_LW, NPD_DROP_TYPE_LW)                          &
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM_LW(NPD_DROP_TYPE_LW)                       &
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM_LW(NPD_DROP_TYPE_LW)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER                                                           &
     &     N_AEROSOL_LW                                                 &
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL_LW(NPD_AEROSOL_SPECIES_LW)                      &
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION_LW(NPD_AEROSOL_SPECIES_LW)         &
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY_LW(NPD_AEROSOL_SPECIES_LW)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL                                                           &
     &     L_AEROSOL_SPECIES_LW(NPD_AEROSOL_SPECIES_LW)
!             AEROSOL SPECIES INCLUDED
!
      REAL                                                              &
     &     AEROSOL_ABSORPTION_LW(NPD_HUMIDITIES_LW                      &
     &        , NPD_AEROSOL_SPECIES_LW, NPD_BAND_LW)                    &
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING_LW(NPD_HUMIDITIES_LW                      &
     &        , NPD_AEROSOL_SPECIES_LW, NPD_BAND_LW)                    &
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY_LW(NPD_HUMIDITIES_LW                       &
     &        , NPD_AEROSOL_SPECIES_LW, NPD_BAND_LW)                    &
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES_LW(NPD_HUMIDITIES_LW, NPD_AEROSOL_SPECIES_LW)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER                                                           &
     &     I_ICE_PARAMETRIZATION_LW(NPD_ICE_TYPE_LW)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL                                                           &
     &     L_ICE_TYPE_LW(NPD_ICE_TYPE_LW)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL                                                              &
     &     ICE_PARAMETER_LIST_LW(NPD_CLOUD_PARAMETER_LW                 &
     &        , NPD_BAND_LW, NPD_ICE_TYPE_LW)                           &
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM_LW(NPD_ICE_TYPE_LW)                         &
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM_LW(NPD_ICE_TYPE_LW)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL                                                           &
     &     L_DOPPLER_PRESENT_LW(NPD_SPECIES_LW)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL                                                              &
     &     DOPPLER_CORRECTION_LW(NPD_SPECIES_LW)
!             OFFSET TO PRESSURE TO REPRESENT DOPPLER BROADENING
!
!
!
!     FIELDS FOR AEROSOL OPTICAL DEPTH
!
      INTEGER                                                           &
     &     N_AOD_WAVEL_LW
!             NUMBER OF WAVELENGTHS
      REAL                                                              &
     &     AOD_ABSORPTION_LW(NPD_HUMIDITIES_LW                          &
     &       , NPD_AEROSOL_SPECIES_LW, NPD_AOD_WAVEL_LW)
!             MONOCHROMATIC SPECIFIC ABSORPTION COEFFICIENT
      REAL                                                              &
     &     AOD_SCATTERING_LW(NPD_HUMIDITIES_LW                          &
     &       , NPD_AEROSOL_SPECIES_LW, NPD_AOD_WAVEL_LW)
!             MONOCHROMATIC SPECIFIC SCATTERING COEFFICIENT
      INTEGER                                                           &
     &     I_AOD_TYPE_LW(NPD_AEROSOL_SPECIES_LW)
!             RELATIONSHIP BETWEEN AEROSOL COMPONENT AND TYPE
!
!
!
!    ------------------------------------------------------------------
