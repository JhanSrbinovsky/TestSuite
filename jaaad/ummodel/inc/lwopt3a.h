!     ------------------------------------------------------------------
!     Module defining LW options for the radiation code.
!
!     NOTE: LWOPT3A and LWCAVR3A must be consistent.
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
#if defined(A02_3C) || defined(A02_3Z)
      CHARACTER  (LEN=132) :: spectral_file_lw
!       Spectral file
      INTEGER :: first_band_lw
!       first band in lw calculations
      INTEGER :: last_band_lw
!       last band in lw calculation
#endif
      INTEGER                                                           &
     &     I_2STREAM_LW                                                 &
!             TWO-STREAM SCHEME
     &   , I_GAS_OVERLAP_LW                                             &
!             TREATMENT OF GASEOUS OVERLAPS
     &   , I_CLOUD_LW                                                   &
!             TREATMENT OF CLOUDY OVERLAPS
     &   , I_CLOUD_REPRESENTATION_LW                                    &
!             REPRESENTATION OF CLOUDS
     &   , I_SOLVER_LW                                                  &
!             SOLVER SELECTED
     &   , I_SCATTER_METHOD_LW
!             TREATMENT OF SCATTERING (THIS CONTROLS THE CALCULATION
!             OF OPTICAL PROPERTIES)
#if defined(A02_3C) || defined(A02_3Z)
      INTEGER :: i_angular_integration_lw
!       Method of angular integration
      INTEGER :: i_truncation_lw
!       Type of spherical truncation
      INTEGER :: ls_global_trunc_lw
!       Order of truncation
      INTEGER :: ms_min_lw
!       Minimum azimuthal order
      INTEGER :: ms_max_lw
!       Maximum azimuthal order
      INTEGER :: ls_brdf_trunc_lw
!       Order of truncation applied to brdfs
      INTEGER :: i_sph_algorithm_lw
!       Algorithm employed in calculating spherical harmonics
      INTEGER :: i_sph_mode_lw
!       Mode of application of spherical harmonic code
      REAL  :: accuracy_adaptive_lw
!       Accuracy of the adaptive truncation in the lw
#endif
      LOGICAL                                                           &
     &     L_IR_SOURCE_QUAD_LW                                          &
!             REPRESENTATION OF THE IR-SOURCE TERM
     &   , L_MICROPHYSICS_LW                                            &
!             FLAG FOR MICROPHYSICS IN LONG WAVE
     &   , L_LOCAL_CNV_PARTITION_LW                                     &
!             FLAG TO PARTITION CONVECTIVE CLOUD USING THE
!             LOCAL TEMPERATURE.
     &   , L_EXTRA_TOP_LW
!             Flag to include an extra layer at the top of the
!             atmosphere in radiative calculations
#if defined(A02_3C) || defined(A02_3Z)
      LOGICAL :: l_henyey_greenstein_pf_lw
!       Flag to use henyey_greenstein phase functions
      LOGICAL :: l_euler_trnf_lw
!       Flag to apply euler's transformation to series
!       of spherical harmonic terms
#endif
!     OPTIONS FOR TRACE GASES:
      LOGICAL                                                           &
     &     L_N2O_LW                                                     &
!             FLAG FOR NITROUS OXIDE
     &   , L_CH4_LW                                                     &
!             FLAG FOR METHANE
     &   , L_CFC11_LW                                                   &
!             FLAG FOR CFC11
     &   , L_CFC12_LW                                                   &
!             FLAG FOR CFC12
     &   , L_CFC113_LW                                                  &
!             FLAG FOR CFC113
     &   , L_HCFC22_LW                                                  &
!             FLAG FOR HCFC22
     &   , L_HFC125_LW                                                  &
!             FLAG FOR HFC125
     &   , L_HFC134A_LW
!             FLAG FOR HFC134A
!
!     TYPES OF DROPLETS OR ICE CRYSTALS USED FOR PARAMETRIZATIONS
      INTEGER                                                           &
     &     I_ST_WATER_LW                                                &
!             TYPE FOR STRATIFORM WATER
     &   , I_CNV_WATER_LW                                               &
!             TYPE FOR CONVECTIVE WATER
     &   , I_ST_ICE_LW                                                  &
!             TYPE FOR STRATIFORM ICE
     &   , I_CNV_ICE_LW
!             TYPE FOR CONVECTIVE ICE
!
!     ------------------------------------------------------------------
