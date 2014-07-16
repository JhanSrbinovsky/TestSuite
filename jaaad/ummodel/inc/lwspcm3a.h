!     ------------------------------------------------------------------
!     MODULE DECLARING COMMON BLOCK CONTAINING THE REDUCED LW SPECTRAL
!     FILE.
!     (NOTE: LWSPDC3A, LWSPCM3A AND LWSARG3A MUST BE CONSISTENT)
!
!
      COMMON/R2LWSPCM/                                                  &
!
!     DIMENSIONS OF ARRAYS:
     &     NPD_TYPE_LW, NPD_BAND_LW, NPD_EXCLUDE_LW                     &
     &   , NPD_SPECIES_LW, NPD_ESFT_TERM_LW, NPD_SCALE_FNC_LW           &
     &   , NPD_SCALE_VARIABLE_LW                                        &
     &   , NPD_THERMAL_COEFF_LW                                         &
     &   , NPD_SURFACE_LW, NPD_ALBEDO_PARM_LW                           &
     &   , NPD_CONTINUUM_LW                                             &
     &   , NPD_DROP_TYPE_LW, NPD_ICE_TYPE_LW, NPD_CLOUD_PARAMETER_LW    &
     &   , NPD_AEROSOL_SPECIES_LW, NPD_HUMIDITIES_LW, NPD_AOD_WAVEL_LW  &
!
!     GENERAL ARRAYS:
     &   , L_PRESENT_LW                                                 &
!
!     PROPERTIES OF BANDS:
     &   , N_BAND_LW, WAVE_LENGTH_SHORT_LW, WAVE_LENGTH_LONG_LW         &
!
!     EXCLUSIONS FROM BANDS:
     &   , N_BAND_EXCLUDE_LW, INDEX_EXCLUDE_LW                          &
!
!     SOLAR FIELDS:
     &   , SOLAR_FLUX_BAND_LW, RAYLEIGH_COEFFICIENT_LW                  &
!
!     GASEOUS ABSORPTION:
     &   , N_ABSORB_LW, N_BAND_ABSORB_LW, INDEX_ABSORB_LW               &
     &   , TYPE_ABSORB_LW                                               &
     &   , I_BAND_ESFT_LW, I_SCALE_ESFT_LW, I_SCALE_FNC_LW              &
     &   , K_ESFT_LW, W_ESFT_LW                                         &
     &   , SCALE_VECTOR_LW, P_REFERENCE_LW, T_REFERENCE_LW              &
!
!     THERMAL SOURCE FUNCTION:
     &   , N_DEG_FIT_LW, THERMAL_COEFFICIENT_LW, T_REF_PLANCK_LW        &
!
!     SURFACE PROPERTIES:
     &   , I_SPEC_SURFACE_LW, N_DIR_ALBEDO_FIT_LW, L_SURFACE_LW         &
     &   , SURFACE_ALBEDO_LW, DIRECT_ALBEDO_PARM_LW                     &
     &   , EMISSIVITY_GROUND_LW                                         &
!
!     CONTINUA:
     &   , N_BAND_CONTINUUM_LW, INDEX_CONTINUUM_LW, INDEX_WATER_LW      &
     &   , I_SCALE_FNC_CONT_LW, K_CONTINUUM_LW                          &
     &   , SCALE_CONTINUUM_LW, P_REF_CONTINUUM_LW, T_REF_CONTINUUM_LW   &
!
!     WATER DROPLETS:
     &   , I_DROP_PARAMETRIZATION_LW, L_DROP_TYPE_LW                    &
     &   , DROP_PARAMETER_LIST_LW                                       &
     &   , DROP_PARM_MIN_DIM_LW, DROP_PARM_MAX_DIM_LW                   &
!
!     AEROSOLS:
     &   , N_AEROSOL_LW, TYPE_AEROSOL_LW, I_AEROSOL_PARAMETRIZATION_LW  &
     &   , NHUMIDITY_LW, HUMIDITIES_LW, L_AEROSOL_SPECIES_LW            &
     &   , AEROSOL_ABSORPTION_LW, AEROSOL_SCATTERING_LW                 &
     &   , AEROSOL_ASYMMETRY_LW                                         &
!
!     ICE CRYSTALS:
     &   , I_ICE_PARAMETRIZATION_LW, L_ICE_TYPE_LW                      &
     &   , ICE_PARAMETER_LIST_LW                                        &
     &   , ICE_PARM_MIN_DIM_LW, ICE_PARM_MAX_DIM_LW                     &
!
!     DOPPLER BROADENING:
     &   , L_DOPPLER_PRESENT_LW, DOPPLER_CORRECTION_LW                  &
!
!     AEROSOL OPTICAL DEPTH
     &   , N_AOD_WAVEL_LW, AOD_ABSORPTION_LW, AOD_SCATTERING_LW         &
     &   , I_AOD_TYPE_LW
!
!     ------------------------------------------------------------------
