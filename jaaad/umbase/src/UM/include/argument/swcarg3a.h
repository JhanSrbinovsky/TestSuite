!     ------------------------------------------------------------------
!     ARGUMENT LIST OF CONTROLLING OPTIONS FOR SHORTWAVE RADIATION.
!
     &   ,                                                              &
     &  I_2STREAM_SW, I_GAS_OVERLAP_SW, I_CLOUD_SW,                     &
     &  I_CLOUD_REPRESENTATION_SW, I_SOLVER_SW,                         &
     &  L_O2_SW,                                                        &
     &  I_ST_WATER_SW, I_CNV_WATER_SW, I_ST_ICE_SW, I_CNV_ICE_SW,       &
     &  L_LOCAL_CNV_PARTITION_SW, L_EXTRA_TOP_SW                        &
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
     &  ,SPECTRAL_FILE_SW, FIRST_BAND_SW, LAST_BAND_SW,                 &
     &  I_ANGULAR_INTEGRATION_SW, I_TRUNCATION_SW,LS_GLOBAL_TRUNC_SW,   &
     &  MS_MIN_SW, MS_MAX_SW, LS_BRDF_TRUNC_SW, I_SPH_ALGORITHM_SW,     &
     &  N_ORDER_PHASE_SOLAR_SW, I_SPH_MODE_SW,  ACCURACY_ADAPTIVE_SW,   &
     &  L_HENYEY_GREENSTEIN_PF_SW, L_EULER_TRNF_SW                      &
#endif
!
!     ------------------------------------------------------------------
