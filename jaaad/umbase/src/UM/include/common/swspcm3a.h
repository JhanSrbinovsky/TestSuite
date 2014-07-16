! SWSPCM3A declares common block containing the reduced sw spectral
! file.(note: SWSPDC3A, SWSPCM3A and SWSARG3A must be consistent)

      COMMON/R2SWSPCM/                                                  &
      ! dimensions of arrays:
     &  NPD_TYPE_SW, NPD_BAND_SW, NPD_EXCLUDE_SW,                       &
     &  NPD_SPECIES_SW, NPD_ESFT_TERM_SW, NPD_SCALE_FNC_SW,             &
     &  NPD_SCALE_VARIABLE_SW,                                          &
     &  NPD_THERMAL_COEFF_SW,                                           &
     &  NPD_SURFACE_SW, NPD_ALBEDO_PARM_SW,                             &
     &  NPD_CONTINUUM_SW,                                               &
     &  NPD_DROP_TYPE_SW, NPD_ICE_TYPE_SW, NPD_CLOUD_PARAMETER_SW,      &
     &  NPD_AEROSOL_SPECIES_SW, NPD_HUMIDITIES_SW,                      &
     &  NPD_AOD_WAVEL_SW,                                               &
      ! general arrays:
     &  L_PRESENT_SW,                                                   &
      ! properties of bands:
     &  N_BAND_SW, WAVE_LENGTH_SHORT_SW, WAVE_LENGTH_LONG_SW,           &
      ! exclusions from bands:
     &  N_BAND_EXCLUDE_SW, INDEX_EXCLUDE_SW,                            &
      ! solar fields:
     &  SOLAR_FLUX_BAND_SW, RAYLEIGH_COEFFICIENT_SW,                    &
      ! gaseous absorption:
     &  N_ABSORB_SW, N_BAND_ABSORB_SW, INDEX_ABSORB_SW,                 &
     &  TYPE_ABSORB_SW,                                                 &
     &  I_BAND_ESFT_SW, I_SCALE_ESFT_SW, I_SCALE_FNC_SW,                &
     &  K_ESFT_SW, W_ESFT_SW,                                           &
     &  SCALE_VECTOR_SW, P_REFERENCE_SW, T_REFERENCE_SW,                &
      ! thermal source function:
     &  N_DEG_FIT_SW, THERMAL_COEFFICIENT_SW, T_REF_PLANCK_SW,          &
      ! surface properties:
     &  I_SPEC_SURFACE_SW, N_DIR_ALBEDO_FIT_SW, L_SURFACE_SW,           &
     &  SURFACE_ALBEDO_SW, DIRECT_ALBEDO_PARM_SW ,                      &
     &  EMISSIVITY_GROUND_SW,                                           &
      ! continua:
     &  N_BAND_CONTINUUM_SW, INDEX_CONTINUUM_SW, INDEX_WATER_SW,        &
     &  I_SCALE_FNC_CONT_SW, K_CONTINUUM_SW,                            &
     &  SCALE_CONTINUUM_SW, P_REF_CONTINUUM_SW, T_REF_CONTINUUM_SW,     &
      ! water droplets:
     &  I_DROP_PARAMETRIZATION_SW, L_DROP_TYPE_SW,                      &
     &  DROP_PARAMETER_LIST_SW,                                         &
     &  DROP_PARM_MIN_DIM_SW, DROP_PARM_MAX_DIM_SW,                     &
      ! aerosols:
     &  N_AEROSOL_SW, TYPE_AEROSOL_SW, I_AEROSOL_PARAMETRIZATION_SW,    &
     &  NHUMIDITY_SW, HUMIDITIES_SW, L_AEROSOL_SPECIES_SW,              &
     &  AEROSOL_ABSORPTION_SW, AEROSOL_SCATTERING_SW,                   &
     &  AEROSOL_ASYMMETRY_SW,                                           &
      ! ice crystals:
     &  I_ICE_PARAMETRIZATION_SW, L_ICE_TYPE_SW,                        &
     &  ICE_PARAMETER_LIST_SW,                                          &
     &  ICE_PARM_MIN_DIM_SW, ICE_PARM_MAX_DIM_SW,                       &
      ! doppler broadening:
     &  L_DOPPLER_PRESENT_SW, DOPPLER_CORRECTION_SW,                    &
      ! aerosol optical depth
     &  N_AOD_WAVEL_SW, AOD_ABSORPTION_SW, AOD_SCATTERING_SW,           &
     &  I_AOD_TYPE_SW

! SWSPCM3A end
