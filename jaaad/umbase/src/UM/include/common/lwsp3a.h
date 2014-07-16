! LWSP3A declares elements of a spectral file as a namelist.
!
! 6.2  24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!

      NAMELIST/R2LWSP/                                                  &
        ! Blocks Present
     &  L_PRESENT,                                                      &
        ! Block 0
     &  N_BAND, N_ABSORB, N_AEROSOL, TYPE_ABSORB, TYPE_AEROSOL,         &
        ! Block 1
     &  WAVE_LENGTH_SHORT, WAVE_LENGTH_LONG,                            &
        ! Block 2
     &  SOLAR_FLUX_BAND,                                                &
        ! Block 3
     &  RAYLEIGH_COEFFICIENT,                                           &
        ! Block 4
     &  N_BAND_ABSORB, INDEX_ABSORB,                                    &
        ! Block 5
     &  I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC,                         &
     &  P_REFERENCE, T_REFERENCE, K_ESFT, W_ESFT, SCALE_VECTOR,         &
        ! Block 6
     &  N_DEG_FIT, T_REF_PLANCK, THERMAL_COEFFICIENT,                   &
        ! Block 7
     &  I_SPEC_SURFACE, N_DIR_ALBEDO_FIT, L_SURFACE,                    &
     &  SURFACE_ALBEDO, DIRECT_ALBEDO_PARM, EMISSIVITY_GROUND,          &
        ! Block 8
     &  N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER,                 &
        ! Block 9
     &  I_SCALE_FNC_CONT, P_REF_CONTINUUM, T_REF_CONTINUUM,             &
     &  K_CONTINUUM, SCALE_CONTINUUM,                                   &
        ! Block 10
     &  I_DROP_PARAMETRIZATION, L_DROP_TYPE, DROP_PARAMETER_LIST,       &
#if defined(A02_3A)
     &  DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM,                           &
#endif
#if defined(A02_3C) || defined(A02_3Z)
     &  N_DROP_PHF_TERM, DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM,          &
#endif
        ! Block 11
     &  I_AEROSOL_PARAMETRIZATION, NHUMIDITY, L_AEROSOL_SPECIES,        &
#if defined(A02_3A)
     &  AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_ASYMMETRY,      &
     &  HUMIDITIES,                                                     &
#endif
#if defined(A02_3C) || defined(A02_3Z)
     &  AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_ASYMMETRY,      &
     &  N_AEROSOL_PHF_TERM, AEROSOL_PHASE_FNC,                          &
     &  HUMIDITIES,                                                     &
#endif
        ! Block 12
     &  I_ICE_PARAMETRIZATION, L_ICE_TYPE, ICE_PARAMETER_LIST,          &
#if defined(A02_3A)
     &  ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM,                             &
#endif
#if defined(A02_3C) || defined(A02_3Z)
     &  N_ICE_PHF_TERM,ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM,              &
#endif
        ! Block 13
     &  L_DOPPLER_PRESENT, DOPPLER_CORRECTION,                          &
        ! Block 14
     &  N_BAND_EXCLUDE, INDEX_EXCLUDE,                                  &
        ! Block 15
     &  N_AOD_WAVEL, AOD_ABSORPTION, AOD_SCATTERING, I_AOD_TYPE

! LWSP3A end
