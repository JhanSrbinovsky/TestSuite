!     ------------------------------------------------------------------
!     Module defining namelists required for input of options to
!     the radiation code.
!
!     Current Owner of Code: J. M. Edwards
!
!     History:
!     Version  Date      Comment.
!     5.2      05/12/00  Logical for extra top level added.
!                                J. M. Edwards
!     5.3      04/10/01  Option for LW scattering added.
!                                J. M. Edwards
! 6.2  24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!
!
!     Elements of namelists for options in the radiation code.
#if defined(A01_3C) || defined(A01_3Z)
      NAMELIST/R2SWNCAL/                                                &
     &     n_swcall
#endif
!
!
      NAMELIST/R2SWCLNL/                                                &
#if defined(A01_3C) || defined(A01_3Z)
     &     spectral_file_sw, first_band_sw, last_band_sw                &
!                       Algorithmic Options
     &   , I_2STREAM_SW, I_GAS_OVERLAP_SW                               &
     &   , I_CLOUD_SW, I_CLOUD_REPRESENTATION_SW                        &
     &   , I_SOLVER_SW                                                  &
     &   , L_O2_SW                                                      &
     &   , I_ST_WATER_SW, I_CNV_WATER_SW, I_ST_ICE_SW, I_CNV_ICE_SW     &
     &   , L_LOCAL_CNV_PARTITION_SW, L_EXTRA_TOP_SW                     &
     &   , i_angular_integration_sw                                     &
     &   , i_truncation_sw, n_order_phase_solar_sw                      &
     &   , ls_global_trunc_sw, ls_brdf_trunc_sw                         &
     &   , ms_min_sw, ms_max_sw, i_sph_algorithm_sw                     &
     &   , accuracy_adaptive_sw, l_euler_trnf_sw                        &
     &   , l_henyey_greenstein_pf_sw, i_sph_mode_sw                     &
     &   , l_subsample                                                  &
     &   , l_geostationary, sat_desc, sat_hgt, sat_lon, sat_lat         &
     &   , max_view_lon, min_view_lon, max_view_lat, min_view_lat
#else
!                       Algorithmic Options
     &     I_2STREAM_SW, I_GAS_OVERLAP_SW                               &
     &   , I_CLOUD_SW, I_CLOUD_REPRESENTATION_SW                        &
     &   , I_SOLVER_SW                                                  &
     &   , L_O2_SW                                                      &
     &   , I_ST_WATER_SW, I_CNV_WATER_SW, I_ST_ICE_SW, I_CNV_ICE_SW     &
     &   , L_LOCAL_CNV_PARTITION_SW, L_EXTRA_TOP_SW
#endif
!
!
#if defined(A02_3C) || defined(A02_3Z)
!
      NAMELIST/R2LWNCAL/                                                &
     &     n_lwcall
#endif
!
      NAMELIST/R2LWCLNL/                                                &
#if defined(A02_3C) || defined(A02_3Z)
     &     spectral_file_lw, first_band_lw, last_band_lw                &
!                       Algorithmic Options
     &   , I_2STREAM_LW, L_IR_SOURCE_QUAD_LW, I_GAS_OVERLAP_LW          &
     &   , I_CLOUD_LW, I_CLOUD_REPRESENTATION_LW                        &
     &   , I_SOLVER_LW, I_SCATTER_METHOD_LW                             &
     &   , L_N2O_LW, L_CH4_LW, L_CFC11_LW, L_CFC12_LW                   &
     &   , L_CFC113_LW, L_HCFC22_LW, L_HFC125_LW, L_HFC134A_LW          &
     &   , I_ST_WATER_LW, I_CNV_WATER_LW, I_ST_ICE_LW, I_CNV_ICE_LW     &
     &   , L_MICROPHYSICS_LW, L_LOCAL_CNV_PARTITION_LW                  &
     &   , L_EXTRA_TOP_LW                                               &
     &   , i_angular_integration_lw, i_truncation_lw                    &
     &   , ls_global_trunc_lw, ls_brdf_trunc_lw                         &
     &   , ms_min_lw, ms_max_lw, i_sph_algorithm_lw                     &
     &   , accuracy_adaptive_lw, l_euler_trnf_lw                        &
     &   , l_henyey_greenstein_pf_lw, i_sph_mode_lw                     &
     &   , l_subsample                                                  &
     &   , l_geostationary, sat_desc, sat_hgt, sat_lon, sat_lat         &
     &   , max_view_lon, min_view_lon, max_view_lat, min_view_lat
#else
!                       Algorithmic Options
     &     I_2STREAM_LW, L_IR_SOURCE_QUAD_LW, I_GAS_OVERLAP_LW          &
     &   , I_CLOUD_LW, I_CLOUD_REPRESENTATION_LW                        &
     &   , I_SOLVER_LW, I_SCATTER_METHOD_LW                             &
     &   , L_N2O_LW, L_CH4_LW, L_CFC11_LW, L_CFC12_LW                   &
     &   , L_CFC113_LW, L_HCFC22_LW, L_HFC125_LW, L_HFC134A_LW          &
     &   , I_ST_WATER_LW, I_CNV_WATER_LW, I_ST_ICE_LW, I_CNV_ICE_LW     &
     &   , L_MICROPHYSICS_LW, L_LOCAL_CNV_PARTITION_LW                  &
     &   , L_EXTRA_TOP_LW
#endif
!
!
!     ------------------------------------------------------------------
