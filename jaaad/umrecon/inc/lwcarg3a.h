!     ------------------------------------------------------------------
!     ARGUMENT LIST OF CONTROLLING OPTIONS FOR THE LONGWAVE RADIATION.
!
     &   ,                                                              &
     &  I_2STREAM_LW, L_IR_SOURCE_QUAD_LW, I_GAS_OVERLAP_LW,            &
     &  I_CLOUD_LW, I_CLOUD_REPRESENTATION_LW, I_SOLVER_LW,             &
     &  I_SCATTER_METHOD_LW,                                            &
     &  L_N2O_LW, L_CH4_LW, L_CFC11_LW , L_CFC12_LW, L_CFC113_LW,       &
     &  L_HCFC22_LW, L_HFC125_LW, L_HFC134A_LW,                         &
     &  I_ST_WATER_LW, I_CNV_WATER_LW, I_ST_ICE_LW,                     &
     &  I_CNV_ICE_LW, L_MICROPHYSICS_LW,                                &
     &  L_LOCAL_CNV_PARTITION_LW, L_EXTRA_TOP_LW                        &
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
     &  ,SPECTRAL_FILE_LW, FIRST_BAND_LW, LAST_BAND_LW,                 &
     &  I_ANGULAR_INTEGRATION_LW, I_TRUNCATION_LW,LS_GLOBAL_TRUNC_LW,   &
     &  MS_MIN_LW, MS_MAX_LW, LS_BRDF_TRUNC_LW, I_SPH_ALGORITHM_LW,     &
     &  I_SPH_MODE_LW, ACCURACY_ADAPTIVE_LW,                            &
     &  L_HENYEY_GREENSTEIN_PF_LW, L_EULER_TRNF_LW                      &
#endif
!
!     ------------------------------------------------------------------
