

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interface between Atmos_physics1 and the radiation code
!
! Purpose:
!   This glue routine has been inserted in order to
!   manage the versions 3A,3C and 3Z of the radiation code
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      Subroutine NI_RAD_CTL (                                           &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc            &
     &, n_procx, n_procy, neighbour, g_rows, g_row_length               &
     &, g_datastart, me                                                 &

! model dimensions.
     &, row_length, rows, n_rows                                        &
     &, model_levels, wet_model_levels, bl_levels                       &
     &, Ozone_levels, cloud_levels, N_cca_levels                        &
     &, NTILES, LAND_FIELD, DUST_DIM1, DUST_DIM2, biogenic_dim1         &
     &, biogenic_dim2, sulp_dim1, sulp_dim2, soot_dim1, soot_dim2       &
     &, bmass_dim1, bmass_dim2, salt_dim1, salt_dim2, salt_dim3         &
     &, co2_dim_len, co2_dim_row, co2_dim2, arcl_dim1, arcl_dim2        &
     &, n_arcl_species, n_arcl_compnts, i_arcl_compnts                  &
     &, ocff_dim1, ocff_dim2                                            &

! Model switches
     &, model_domain                                                    &
     &, L_Rad_Step, L_Rad_Step_diag, L_Rad_Step_prog                    &
     &, L_Forcing, L_Timestep, L_Radiance, L_Wenyi                      &
     &, L_CAL360, L_SEC_VAR, L_EqT                                      &
     &, L_MICROPHY, L_INHOM_CLOUD, L_USE_DUST, L_USE_BIOGENIC           &
     &, L_use_sulpc_direct, L_use_soot_direct, L_use_soot_indirect      &
     &, L_use_bmass_direct, L_use_bmass_indirect                        &
     &, L_use_ocff_direct, L_use_ocff_indirect                          &
     &, L_use_sulpc_indirect_SW, L_use_sulpc_indirect_LW                &
     &, L_emcorr,L_climat_aerosol,L_clim_aero_hgt,L_HadGEM1_Clim_Aero   &
     &, Ltimer, L_ssice_albedo                                          &
     &, L_cable,L_snow_albedo                                           &
     &, L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a              &
     &, L_use_seasalt_indirect                                          &
     &, L_use_seasalt_direct                                            &
     &, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                           &
     &, L_pc2, l_mixing_ratio                                           &
     &, L_MURK_RAD                                                      &
     &, L_rad_deg                                                       &
     &, L_co2_interactive                                               &
! kdcorbin, 06/10 - added L_CO2_RADIATION
     &, L_CO2_RADIATION                                                 &
     &, L_use_stochem_CH4, L_use_clearrh, L_ukca, L_USE_ARCL            &
     &, L_USE_SPEC_SEA, L_MOD_BARKER_ALBEDO, L_MOD_K_FLUX               &
     &, L_SCVARY, L_VOLCTS                                              &

! model Parameters
     &, A_SW_segments, A_SW_seg_size, A_LW_segments, A_LW_seg_size      &
     &, INHOM_CLOUD_SW, INHOM_CLOUD_LW, DP_CORR_STRAT, DP_CORR_CONV     &
     &, CO2_MMR, O2MMR, N2O_mix_ratio, CH4_mix_ratio                    &
     &, CFC11_mix_ratio, CFC12_mix_ratio                                &
     &, C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR                       &
     &, SW_alpham, SW_alphac, SW_alphab, SW_dtice                       &
     &, dt_bare,dalb_bare_wet,pen_rad_frac,SW_beta                      &
     &, min_trop_level, max_trop_level                                  &
     &, Ntot_land, Ntot_sea, aero_bl_levels                             &

! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &

! in coordinate information
     &, rho_r2, r_rho_levels, r_theta_levels                            &
     &, eta_theta_levels, eta_rho_levels                                &
     &, delta_lambda, delta_phi                                         &
     &, lat_rot_NP, long_rot_NP                                         &

! in time stepping information.
     &, timestep, radiation_timestep                                    &
     &, radiation_tstep_diag, radiation_tstep_prog                      &
     &, val_year, val_day_number, val_hour, val_minute                  &
     &, val_second, timestep_number                                     &
     &,         PREVIOUS_TIME                                           &

! trig arrays
     &, sin_theta_longitude, cos_theta_longitude                        &
     &, FV_cos_theta_latitude                                           &

! grid-dependent arrays
     &, f3_at_u, true_longitude, true_latitude                          &

! diagnostic info
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork1,                                                      &
     & STASHwork2                                                       &

! SCM diagnostics switches (dummy in full UM)
     &, nSCMDpkgs,L_SCMDiags                                            &

! in data fields.
     &, p_star                                                          &
     &, p_layer_boundaries, p_layer_centres                             &
     &, p, p_theta_levels                                               &
     &, exner_rho_levels, exner_theta_levels                            &
     &, land_sea_mask, fland, land0p5, l_ctile                          &
     &, T_SURF,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE,AREA_CLOUD_FRACTION      &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, biogenic, so4_aitken, so4_accu, so4_diss                        &
     &, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld   &
     &, ocff_new, ocff_agd, ocff_cld, aerosol, arcl                     &
     &, co2_3D                                                          &
     &, ch4_stochem                                                     &
     &, cos_zenith_angle                                                &
     &, can_rad_mod                                                     &
     &, frac_control                                                    &
! chemical greenhouse gas fields
     &, ngrgas, grgas_field                                             &

! ancillary fields and fields needed to be kept from timestep to
! timestep
     &, snow_depth, snow_depth_sea, ice_fract, rgrain, soot             &
     &, cca, ccb, cct, cclwp, ccw, lcbase                               &
     &, ozone, SW_incs, LW_incs, dirpar_inc                             &
     &, O3_trop_level, O3_trop_height                                   &
     &, T_trop_level, T_trop_height, zh                                 &
     &, land_index, albsoil, lai, snow_tile                             &
     &, frac, tstar_tile, z0_tile                                       &
     &, dOLR_rts, LW_down, SW_tile_rts                                  &
     &, u_1, v_1                                                        &
     &, land_alb,sice_alb                                               &
     &, ES_SPACE_INTERP, A_SW_radstep, A_LW_radstep                     &
     &, A_SW_radstep_diag, A_SW_radstep_prog                            &
     &, A_LW_radstep_diag, A_LW_radstep_prog, RAD_MASK                  &

! in/out
     &, T_n, q_n, qcl_n, qcf_n, cf_n, cfl_n, cff_n                      &
     &, qcf2_n, qrain_n, qgraup_n                                       &
     &, T_inc, q_inc, qcl_inc, cf_inc, cfl_inc                          &
     &, sum_eng_fluxes                                                  &

! out.
     &, photosynth_act_rad, surf_radflux, rad_hr, dOLR, SW_tile         &
!    EAK
     &, istep_cur                                                       &
     &, TILE_PTS,SM_LEVELS,TILE_INDEX                                   &
     &, surf_down_sw,alb_tile                                           &
     &, SNOW_TMP3L,SNOW_RHO1L,TSOIL_TILE,ISNOW_FLG3L,LAND_ALBEDO        &
! Section Information
     &, maxsects, h_sect                                                &
! error information
     &, Error_code  )

       Implicit None

! No choice as we are using new microphysics
      Logical                                                           &
     & L_cloud_water_partition

      Parameter (L_cloud_water_partition = .true.)


! arguments with intent in. ie: input variables.

! Parallel setup variables
      Integer                                                           &
     &  halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, global_row_length                                               &
                           ! number of points on a row
     &, global_rows                                                     &
                           ! NUMBER OF global rows
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbours
                             ! in the horizontal plane
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)                                      &
     &, me         ! My processor number

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north, sout
                         ! east or west of the processor grid
!
! Model dimensions
!
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, model_levels                                                    &
     &, wet_model_levels                                                &
     &, bl_levels                                                       &
     &, Ozone_levels                                                    &
     &, cloud_levels                                                    &
     &, N_cca_levels                                                    &
     &, ntiles                                                          &
     &, land_field                                                      &
     &, DUST_DIM1                                                       &
                  !dimensions of mineral dust arrays
     &, DUST_DIM2                                                       &
                  !dimensions of biogenic aerosol arrays
     &, biogenic_dim1                                                   &
     &, biogenic_dim2                                                   &
     &, sulp_dim1                                                       &
                                ! dimensions of S Cyc arrays
     &, sulp_dim2                                                       &
     &, soot_dim1                                                       &
                                ! dimensions of soot arrays
     &, soot_dim2                                                       &
     &, bmass_dim1                                                      &
                                ! dimensions of biomass arrays
     &, bmass_dim2                                                      &
     &, ocff_dim1                                                       &
     &, ocff_dim2                                                       &
                                ! dimensions of fossil-fuel arrays
     &, salt_dim1                                                       &
                                ! dimensions of sea-salt arrays
     &, salt_dim2                                                       &
     &, salt_dim3                                                       &
     &, co2_dim2                                                        &
                               ! dimensions of CO2 array
     &, co2_dim_len                                                     &
     &, co2_dim_row                                                     &
                               ! dimensions of aerosol clim for NWP
     &, arcl_dim1                                                       &
     &, arcl_dim2
!
! Model switches
!
      Integer                                                           &
     &  model_domain

      Logical                                                           &
     &  L_Rad_Step                                                      &
                             ! true on radiation timestep (3A)
     &, L_Rad_Step_diag                                                 &
                             ! true if fast radiation timestep (3C)
     &, L_Rad_Step_prog                                                 &
                             ! true if slow radiation timestep (3C)
     &, L_Timestep                                                      &
                             ! true if improved timestepping is used
     &, L_Forcing                                                       &
                             ! true if radiative frocing required
     &, L_Radiance                                                      &
                             ! true if radiances are to be calculated
     &, L_Wenyi                                                         &
                        ! true if Wenyi scaling required  (3A/3C)
     &, L_CAL360                                                        &
                        ! true if using 360 day calender
     &, L_SEC_VAR                                                       &
                        ! True if using time varying astronomy
     &, L_EqT                                                           &
                        ! true if including the equation of time
     &, L_microphy                                                      &
                        ! true if using cloud microphysics (set to true)
     &, L_INHOM_CLOUD                                                   &
                        ! switch to simulate inhomogeneous cloud
     &, L_climat_aerosol                                                &
                        ! true
     &, L_clim_aero_hgt                                                 &
                          ! True if using the prognostic BL depth to
!                         ! specify the BL component of the aerosol
!                         ! climatology.
     &, L_HadGEM1_Clim_Aero                                             &
                             ! True if using HadGEM1 setting for
!                            ! climatological aerosols in FILL3A/3Z
     &, L_USE_DUST                                                      &
                   ! include mineral dust direct effect
     &, L_USE_BIOGENIC                                                  &
                   ! include biogenic aerosol direct effect
     &, L_use_sulpc_direct                                              &
                                ! Include direct effect of sulphate
     &, L_use_soot_direct                                               &
                                ! Include direct effect of soot
     &, L_use_bmass_direct                                              &
                                ! Include direct effect of biomass smoke
     &, L_use_ocff_direct                                               &
                                ! Include direct effect of OCFF
     &, L_use_sulpc_indirect_SW                                         &
                                ! Use sulphate aerosol to determine
     &, L_use_sulpc_indirect_LW                                         &
                                ! cloud droplet number in SW and LW
!                               ! radiation respectively.
     &, L_use_soot_indirect                                             &
                                ! Include indirect effect of soot
     &, L_use_bmass_indirect                                            &
                                ! Include indirect effect of biomass
     &, L_use_ocff_indirect                                             &
                                ! Include indirect effect of OCFF
     &, L_emcorr                                                        &
                    ! true if energy correction scheme is to be used.
     &, Ltimer                                                          &
                 ! true then output some timing information
     &, L_ssice_albedo                                                  &
                       ! Switch on the effect of snow on sea-ice albedo
     &, L_sice_meltponds                                                &
                         ! true if sea ice meltponds albedo scheme
     &, L_sice_scattering                                               &
                          ! true if sea ice internal scattering scheme
     &, L_sice_hadgem1a                                                 &
                        ! true if HadGEM1 albedo bug corrected
     &, L_MOSES_II                                                      &
                       ! Always true - should really be removed
     &, L_snow_albedo                                                   &
                       ! True if spectral albedo scheme selected
     &, L_use_seasalt_indirect                                          &
                              ! Switch for indirect effect of sea-salt.
     &, L_use_seasalt_direct                                            &
                              ! Switch for direct effect of sea-salt.
     &, l_mcr_qcf2                                                      &
                       ! Use second ice category
     &, l_mcr_qrain                                                     &
                       ! Use prognostic rain
     &, l_mcr_qgraup                                                    &
                       ! Use graupel
     &, L_pc2                                                           &
                       ! True if using the PC2 cloud scheme
     &, L_mixing_ratio                                                  &
                           ! True if using mixing ratios
     &, L_MURK_RAD                                                      &
                           ! True if using radiative effects of 'murk'.
     &, l_ctile                                                         &
                 !coastal tiling switch
     &, L_rad_deg                                                       &
                   ! Controls the spatial degradation of E-S code
     &, L_co2_interactive                                               &
                           ! Controls the use of 3D CO2 field
     &, L_CO2_RADIATION                                                 &
          !Controls the use of CO2 in radiation - kdcorbin, 06/10
     &, L_use_stochem_CH4                                               &
                           ! Control the use of STOCHEM CH4 field
     &, L_ukca                                                          &
                           ! Switch for UKCA sub-model
     &, L_use_clearrh                                                   &
                           ! RH used for hygroscopic aerosols
     &, L_USE_SPEC_SEA                                                  &
                       ! Switch for spectrally dependent sea albedos
     &, L_MOD_BARKER_ALBEDO                                             &
!                      ! Use modified Barker albedo (open sea)
     &, L_SCVARY, L_VOLCTS
!                      ! switches for solar/volcanic forcing

!     Use modulus of fluxes to remove negative effective extinctions
      LOGICAL, INTENT(IN) :: L_MOD_K_FLUX

! physical constants
      Real                                                              &
     &  lc, lf, cp                                                      &
     &, two_Omega                                                       &
                        ! twice Earth's rotation rate
     &, p_zero                                                          &
     &, kappa                                                           &
     &, R, g, Lapse_Rate, earth_radius, Pi

! model parameters
      Real                                                              &
     &  timestep                                                        &
     &, radiation_timestep                                              &
     &, radiation_tstep_diag                                            &
     &, radiation_tstep_prog                                            &
     &, CO2_MMR                                                         &
                        ! set equal to co2_start
     &, O2MMR                                                           &
     &, N2O_mix_ratio                                                   &
     &, CH4_mix_ratio                                                   &
     &, CFC11_mix_ratio                                                 &
     &, CFC12_mix_ratio                                                 &
     &, SW_alphab                                                       &
     &, SW_alphac                                                       &
     &, SW_alpham                                                       &
     &, SW_dtice                                                        &
     &, dt_bare                                                         &
     &, dalb_bare_wet                                                   &
     &, pen_rad_frac                                                    &
     &, SW_beta
      REAL C113MMR      ! CFC113 mmr
      REAL HCFC22MMR    ! HCFC22 mmr
      REAL HFC125MMR    ! HFC125 mmr
      REAL HFC134AMMR   ! HFC134A mmr

      Integer                                                           &
     &  min_trop_level                                                  &
                        ! Lowest permitted level for the tropopause
!                       ! used for radiative purposes.
     &, max_trop_level                                                  &
                        ! Highest permitted level for the tropopause
!                       ! used for radiative purposes.
     &, can_rad_mod     ! which canopy light mod are we using using Q

      Integer                                                           &
     &  A_SW_segments                                                   &
     &, A_SW_seg_size                                                   &
     &, A_LW_segments                                                   &
     &, A_LW_seg_size                                                   &
     &, A_SW_radstep                                                    &
     &, A_LW_radstep                                                    &
     &, A_SW_radstep_diag                                               &
     &, A_SW_radstep_prog                                               &
     &, A_LW_radstep_diag                                               &
     &, A_LW_radstep_prog                                               &

     &, aero_bl_levels
!             Common number of layers taken to be occupied by the
!             boundary-layer aerosol if the boundary layer
!             depth is not used to determine the number separately
!             at each grid-point

! DIMFIX3A defines internal dimensions tied to algorithms for
! two-stream radiation code, mostly for clouds

      ! number of components of clouds
      INTEGER,PARAMETER:: NPD_CLOUD_COMPONENT=4

      ! number of permitted types of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_TYPE=4

      ! number of permitted representations of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_REPRESENTATION=4

      ! number of overlap coefficients for clouds
      INTEGER,PARAMETER:: NPD_OVERLAP_COEFF=18

      ! number of coefficients for two-stream sources
      INTEGER,PARAMETER:: NPD_SOURCE_COEFF=2

      ! number of regions in a layer
      INTEGER,PARAMETER:: NPD_REGION=3

! DIMFIX3A end

      ! Scaling factors to simulate inhomogeneous cloud.
      Real, Dimension(NPD_CLOUD_COMPONENT) :: INHOM_CLOUD_SW
      Real, Dimension(NPD_CLOUD_COMPONENT) :: INHOM_CLOUD_LW

      ! Decorrelation pressure scale for large scale cloud
      Real :: DP_CORR_STRAT
      ! Decorrelation pressure scale for convective cloud
      Real :: DP_CORR_CONV

! CSUBMODL start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.5    07/04/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!  1. Internal model and submodel dump partition identifiers - fixed
!     for all experiments.
! CSMID start
!
! Description:
!    Hold parameters defining internal model identifiers and submodel
!    data partition (ie main D1 data array and consequent dump), both
!    short and long form.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.3    26/10/93   M. Carter. Part of an extensive mod that:
!                    1.Removes the limit on primary STASH item numbers.
!                    2.Removes the assumption that (section,item)
!                      defines the sub-model.
!                    3.Thus allows for user-prognostics.
!                    Add index to submodel home dump.
! 3.5    13/03/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
! 6.0    02/07/03   Add X_IM and X_SM for small exec.      E.Leung
!
! Declarations:
!
!   Hold parameters defining internal model identifiers and submodel
!   data partition (ie main D1 data array and consequent dump), both
!   short and long form
      ! Internal models
      INTEGER,PARAMETER:: A_IM      = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: ATMOS_IM  = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: O_IM      = 2 ! Ocean internal model
      INTEGER,PARAMETER:: OCEAN_IM  = 2 ! Ocean internalmodel
      INTEGER,PARAMETER:: S_IM      = 3 ! Slab internal model
      INTEGER,PARAMETER:: SLAB_IM   = 3 ! Slab internal model
      INTEGER,PARAMETER:: W_IM      = 4 ! Wave internal model
      INTEGER,PARAMETER:: WAVE_IM   = 4 ! Wave internal model
      INTEGER,PARAMETER:: I_IM      = 5 ! Sea=ice internal model
      INTEGER,PARAMETER:: SEAICE_IM = 5 ! Sea=ice internal model
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_IM      = 6 ! ND internal model
      INTEGER,PARAMETER:: NATMOS_IM = 6 ! ND internal model
      ! Small Executables
      INTEGER,PARAMETER:: X_IM      = 7 ! SX indicator

      ! Submodels
      INTEGER,PARAMETER:: A_SM      = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: ATMOS_SM  = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: O_SM      = 2 ! Ocean submodel
      INTEGER,PARAMETER:: OCEAN_SM  = 2 ! Ocean submodel
      INTEGER,PARAMETER:: W_SM      = 4 ! Wave submodel
      INTEGER,PARAMETER:: WAVE_SM   = 4 ! Wave submodel
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_SM      = 6 ! ND submodel
      INTEGER,PARAMETER:: NATMOS_SM = 6 ! ND submodel
      ! Small Executables
      INTEGER,PARAMETER:: X_SM      = 7 ! SX indicator

! CSMID end

!
!  2. Maximum internal model/submodel array sizes for this version.
!
! CSUBMAX start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.5    13/07/95   Original code. D.M. Goddard
! 4.0     3/11/95   Reduce max internal model, submodel from 10 to 4
!                   to save space in model. At 4.0 the max no of
!                   supported models is 3, 1 slot is reserved for
!                   expansion. Rick Rawlins.
!  4.1  21/02/96  Wave model introduced as 4th sub-model.  RTHBarnes
!
! Declarations:
!
!
!  1. Maximum internal model/submodel array sizes for this version.
!
      ! Max no. of internal models
      INTEGER,PARAMETER:: N_INTERNAL_MODEL_MAX=4

      ! Max no. of submodel dump partitions
      INTEGER,PARAMETER:: N_SUBMODEL_PARTITION_MAX=4

      ! Max value of internal model id
      INTEGER,PARAMETER:: INTERNAL_ID_MAX=N_INTERNAL_MODEL_MAX

      ! Max value of submodel dump id
      INTEGER,PARAMETER:: SUBMODEL_ID_MAX=N_SUBMODEL_PARTITION_MAX

! CSUBMAX end
!
!  3. Lists of internal models and their submodel dump partitions -
!     initialised by the user interface - experiment specific.
      INTEGER :: N_INTERNAL_MODEL          ! No. of internal models
      INTEGER :: N_SUBMODEL_PARTITION      ! No. of submodel partitions

      ! Internal models
      INTEGER :: INTERNAL_MODEL_LIST(N_INTERNAL_MODEL_MAX)

      ! Submodel identifier for each internal model in list
      INTEGER :: SUBMODEL_FOR_IM    (N_INTERNAL_MODEL_MAX)

      ! Submodel number for each submodel id
      INTEGER :: SUBMODEL_FOR_SM(N_INTERNAL_MODEL_MAX)

      ! Namelist for information in 3.
      NAMELIST/NSUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,          &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM

      ! 4. Lists calculated in model from user interface supplied arrays
      ! experiment specific.

      ! No of internal models in each submodel partition indexed by sm
      !  identifier
      INTEGER :: N_INTERNAL_FOR_SM(SUBMODEL_ID_MAX)

      ! List of  submodel partition identifiers
      INTEGER :: SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION_MAX)

      ! Submodel partition identifier indexed by internal model identifie
      INTEGER :: SUBMODEL_PARTITION_INDEX(INTERNAL_ID_MAX)

      ! Sequence number of internal model indexed by internal model
      ! identifier: required to map from id to STASH internal model
      ! sequence
      INTEGER :: INTERNAL_MODEL_INDEX(INTERNAL_ID_MAX)


      ! Last internal model within a submodel partition if .TRUE.,
      ! indexed by internal model id.
      LOGICAL :: LAST_IM_IN_SM(INTERNAL_ID_MAX)

      ! Common block for information in 3. and 4.
      COMMON/SUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,             &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM,SUBMODEL_FOR_SM,            &
     &  N_INTERNAL_FOR_SM,SUBMODEL_PARTITION_LIST,                      &
     &  SUBMODEL_PARTITION_INDEX,                                       &
     &  INTERNAL_MODEL_INDEX,                                           &
     &  LAST_IM_IN_SM

!
!  5. Time information specifying coupling frequencies between internal
!     models and submodels, and multipliers, indexed by sequence of
!     internal models and submodels (ie left to right along node tree).
!     {Not required at this release}.
!
! Namelists for information in 5. {Not required at this release}
!
!
!  6. Lists of coupling nodes defining coupling frequencies between
!     internal models and between submodel partitions. (Not defined
!     yet at this release).
!CALL CNODE
!
!  7. Variables dealing with general coupling switches at the control
!     level. {These will require revision at the next release when
!     coupling between internal models is dealt with more generally.
!     Logicals below are set in routine SETGRCTL.}

      ! new internal model next group of timesteps if .true.
      LOGICAL :: new_im

      ! new submodel dump  next group of timesteps if .true.
      LOGICAL :: new_sm

      COMMON/CSUBMGRP/new_im,new_sm

      INTEGER SUBMODEL_IDENT
      COMMON/SUBMODID/SUBMODEL_IDENT
! CSUBMODL end
! TYPSTS starts
! CSUBMODL must be included before this file
!Applicable to all configurations (except MOS variables)
!STASH related variables for describing output requests and space
!management.
!LL
!LL   AUTHOR            Rick Rawlins
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL   3.2             Code creation for Dynamic allocation
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL   3.5  Apr. 95   Sub-Models project.
!LL                  Dimensioning of various STASH arrays altered in
!LL                  accordance with internal model separation scheme.
!LL                  Arrays PPXREF, INDEX_PPXREF deleted as they are no
!LL                  longer required.
!LL                  S.J.Swarbrick
!LL
!
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS

      INTEGER :: MOS_MASK_LEN         ! Size of bit mask for MOS

      COMMON/DSIZE_AO/  MOS_MASK_LEN

! TYPSTSZ end
!LL  Comdeck: CPPXREF --------------------------------------------------
!LL
!LL  Purpose: Holds PARAMETER definitions to describe the structure of
!LL           each STASHmaster file record plus some valid entries.
!LL
!LL  Author    Dr T Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                  Add a PPXREF record for model number.
!LL  4.0   26/07/95  T.Johns.  Add codes for real/int/log data types.
!LL  3.5   10/3/94   Sub-Models project:
!LL                 List of PPXREF addressing codes augmented, in order
!LL                 to include all of the pre_STASH master information
!LL                 in the new PPXREF file.
!LL                 PPXREF_CODELEN increased to 38.
!LL                 PPXREF_IDLEN deleted - no longer relevant.
!LL                   S.J.Swarbrick
!LL  4.1   June 96  Wave model parameters included.
!LL                 ppx_ address parameters adjusted to allow for
!LL                  reading option code as 4x5 digit groups.
!LL                   S.J.Swarbrick
!LL  5.0   29/06/99  Add halo type parameter for new dynamics.
!LL                  New grid codes for LAM boundary conditions
!LL                  D.M. Goddard
!LL  5.1   07/03/00  Fixed/Free format conversion
!LL  5.2   19/09/00  Added ppx_atm_lbc_orog descriptor   P.Burton
!LL  5.3   21/08/01  Added ocean lbc descriptors.   M. J. Bell
!LL  5.3   23/07/01  Add valid pp_lbvc codes referenced in UM. R Rawlins
!LL  5.5   30/01/03  Option code increase from 20 to 30 digits thus
!LL                  requiring option code address range increase by
!LL                  2 so all subsequent addressing codes need to be
!LL                  increased by 2 to make a gap.
!LL                  W Roseblade
!LL
!LL  Logical components covered: C40
!LL
!-----------------------------------------------------------------------
! Primary file record definition
      ! length of ID in a record
      Integer, Parameter :: PPXREF_IDLEN      = 2

      ! total length of characters *WARNING* must be multiple of 4
      ! to avoid overwriting
      Integer, Parameter :: PPXREF_CHARLEN    = 36

      ! number of packing profiles
      Integer, Parameter :: PPXREF_PACK_PROFS = 10

      ! total length of codes = no. of codes (excluding profs)
      ! + pack_profs
      Integer, Parameter :: PPXREF_CODELEN    = 33 + PPXREF_PACK_PROFS

! Derived file record sizes
      ! Assume that an integer is at least 4 bytes long. Wastes some
      ! space on an 8 byte machine.
      ! ppx_charword = 9.
      Integer, Parameter :: PPX_CHARWORD      = ((PPXREF_CHARLEN+3)/4)

      ! read buffer record length
      Integer, Parameter :: PPX_RECORDLEN = PPX_CHARWORD+PPXREF_CODELEN
!
!-----------------------------------------------------------------------
! Addressing codes within PPXREF
      Integer, Parameter ::  ppx_model_number   = 1  ! Model number
                                                     ! address
      Integer, Parameter ::  ppx_section_number = 2  ! Section number
                                                     ! address
      Integer, Parameter ::  ppx_item_number    = 3  ! Item number
                                                     ! address
      Integer, Parameter ::  ppx_version_mask   = 4  ! Version mask
                                                     ! address
      Integer, Parameter ::  ppx_space_code     = 5  ! Space code
                                                     ! address
      Integer, Parameter ::  ppx_timavail_code  = 6  ! Time availability
                                                     !  code  address
      Integer, Parameter ::  ppx_grid_type      = 7  ! Grid type code
                                                     ! address
      Integer, Parameter ::  ppx_lv_code        = 8  ! Level type code
                                                     ! address
      Integer, Parameter ::  ppx_lb_code        = 9  ! First level code
                                                     !  address
      Integer, Parameter ::  ppx_lt_code        =10  ! Last level code
                                                     ! address
      Integer, Parameter ::  ppx_lev_flag       =11  ! Level compression
                                                     !  flag  address
      Integer, Parameter ::  ppx_opt_code       =12  ! Sectional option
                                                     ! code  address
      Integer, Parameter ::  ppx_pt_code        =18  ! Pseudo dimension
                                                     ! type  address
      Integer, Parameter ::  ppx_pf_code        =19  ! First pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_pl_code        =20  ! Last pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_ptr_code       =21  ! Section 0 point-
                                                     ! back code address
      Integer, Parameter ::  ppx_dump_packing   =22  ! Dump packing code
                                                     ! address
      Integer, Parameter ::  ppx_lbvc_code      =23  ! PP LBVC code
                                                     ! address
      Integer, Parameter ::  ppx_rotate_code    =24  ! Rotation code
                                                     ! address
      Integer, Parameter ::  ppx_field_code     =25  ! PP field code
                                                     ! address
      Integer, Parameter ::  ppx_user_code      =26  ! User code address
      Integer, Parameter ::  ppx_meto8_levelcode=27  ! CF level code
                                                     ! address
      Integer, Parameter ::  ppx_meto8_fieldcode=28  ! CF field code
                                                     ! address
      Integer, Parameter ::  ppx_cf_levelcode   =27
      Integer, Parameter ::  ppx_cf_fieldcode   =28
      Integer, Parameter ::  ppx_base_level     =29  ! Base level code
                                                     ! address
      Integer, Parameter ::  ppx_top_level      =30  ! Top level code
                                                     ! address
      Integer, Parameter ::  ppx_ref_lbvc_code  =31  ! Ref level LBVC
                                                     ! code address
      Integer, Parameter ::  ppx_data_type      =32  ! Data type code
                                                     ! address
      Integer, Parameter ::  ppx_halo_type      =33
      Integer, Parameter ::  ppx_packing_acc    =34  ! Packing accuracy
                                                     ! code  address
      Integer, Parameter ::  ppx_pack_acc       =34  ! Must be last:


                                                 ! multiple pack_acc to
                                                 ! fill up remaining
                                                 ! array elements


!-------------------------------------------------------------------
! Valid grid type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_atm_nonstd=0      ! Non-standard atmos
                                                  ! grid
      Integer, Parameter :: ppx_atm_tall=1        ! All T points (atmos)
      Integer, Parameter :: ppx_atm_tland=2       ! Land-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tsea=3        ! Sea-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tzonal=4      ! Zonal field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_tmerid=5      ! Merid field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_uall=11       ! All u points (atmos)
      Integer, Parameter :: ppx_atm_uland=12      ! Land-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_usea=13       ! Sea-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_uzonal=14     ! Zonal field at u
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_umerid=15     ! Merid field at u
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_scalar=17     ! Scalar (atmos)
      Integer, Parameter :: ppx_atm_cuall=18      ! All C-grid (u)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_cvall=19      ! All C-grid (v)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_compressed=21 ! Compressed land
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_ozone=22      ! Field on ozone
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_river=23      ! River routing
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_rim=25        ! Rim type field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_theta=26  ! All T points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_u=27      ! All u points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_v=28      ! All v points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_orog=29   ! Orography field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_ocn_nonstd=30     ! Non-standard ocean
                                                  ! grid
      Integer, Parameter :: ppx_ocn_tcomp=31      ! Compressed T points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_ucomp=32      ! Compressed u points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_tall=36       ! All T points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_uall=37       ! All u points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_cuall=38      ! All C-grid (u)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_cvall=39      ! All C-grid (v)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_tfield=41     ! All non-cyclic T
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_ufield=42     ! All non-cyclic u
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_tzonal=43     ! Zonal n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_uzonal=44     ! Zonal n-c field at
                                                  ! u points (ocean)
      Integer, Parameter :: ppx_ocn_tmerid=45     ! Merid n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_umerid=46     ! Merid n-c field at
                                                  ! u points  (ocean)
      Integer, Parameter :: ppx_ocn_scalar=47     ! Scalar (ocean)
      Integer, Parameter :: ppx_ocn_rim=51        ! Rim type field
                                                  ! (LAM BCs ocean)
      Integer, Parameter :: ppx_ocn_lbc_theta=52  ! Ocean rim fields
      Integer, Parameter :: ppx_ocn_lbc_u=53      ! on T & U grids
      Integer, Parameter :: ppx_wam_all=60        ! All points (wave
                                                  ! model)
      Integer, Parameter :: ppx_wam_sea=62        ! Sea points only
                                                  ! (wave model)
      Integer, Parameter :: ppx_wam_rim=65        ! Rim type field
                                                  ! (LAM BCs wave)

!--------------------------------------------------------------------
! Valid rotation type codes
!--------------------------------------------------------------------
      Integer, Parameter :: ppx_unrotated=0       ! Unrotated output
                                                  ! field
      Integer, Parameter :: ppx_elf_rotated=1     ! Rotated ELF field

!-------------------------------------------------------------------
! Valid level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_full_level=1      ! Model full level
      Integer, Parameter :: ppx_half_level=2      ! Model half level
      Integer, Parameter :: ppx_rho_level=1       ! Model rho level
      Integer, Parameter :: ppx_theta_level=2     ! Model theta level
      Integer, Parameter :: ppx_single_level=5    ! Model single level
      Integer, Parameter :: ppx_soil_level=6      ! Deep Soil level

!-------------------------------------------------------------------
! Valid data type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_type_real=1       ! Real data type
      Integer, Parameter :: ppx_type_int=2        ! Integer data type
      Integer, Parameter :: ppx_type_log=3        ! Logical data type

!-------------------------------------------------------------------
! Valid meto8 level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_meto8_surf=9999   ! MetO8 surface type
                                                  ! code

!-------------------------------------------------------------------
! Valid dump packing codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_pack_off=0        ! Field not packed
                                                  ! (ie. 64 bit)
      Integer, Parameter :: ppx_pack_32=-1        ! Field packed to
                                                  ! 32 bit in  dump
      Integer, Parameter :: ppx_pack_wgdos=1      ! Field packed by
                                                  ! WGDOS method
      Integer, Parameter :: ppx_pack_cfi1=11      ! Field packed using
                                                  ! CFI1  (ocean)

!-------------------------------------------------------------------
! Add valid lbvc codes referenced in model (pp header output labels)
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_lbvc_height  =  1 ! height
      Integer, Parameter :: ppx_lbvc_depth   =  2 ! depth (ocean)
      Integer, Parameter :: ppx_lbvc_pressure=  8 ! pressure
      Integer, Parameter :: ppx_lbvc_theta   = 19 ! potential T
      Integer, Parameter :: ppx_lbvc_hybrid  = 65 ! hybrid height(atmos)
      Integer, Parameter :: ppx_lbvc_PV      = 82 ! potential vorticity
      Integer, Parameter :: ppx_lbvc_surface =129 ! surface
! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      INTEGER :: MOS_OUTPUT_LENGTH
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP,MOS_OUTPUT_LENGTH

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
      INTEGER:: MOS_MASK(MOS_MASK_LEN)
! TYPSTS end
!------------------------ nstypes.h ----------------------------------
!jhan:further renovation of ths file may be necessary params are dependent on dataset
!jhan: ALSO nstypes_cable.h should be unecessary nsoil/soil is only difference
      !--- Number of non-vegetation surface types
      Integer, Parameter :: NNVG  = 4

      !--- Number of plant functional types.
      Integer, Parameter :: NPFT  = 13
      
      !--- Number of surface types.
      Integer, Parameter :: NTYPE =17 
      
      !--- Index of the surface type 'Soil'
      !Integer, Parameter :: SOIL  = 16 
      !dhb599, 20110615: change made as per Peter Vohralik, item 1:
      Integer, Parameter :: SOIL  = 14

!--- Land surface types :
!--- original veg. tiles 
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!--- for testing these tiles are set = 1:5 
!     6 - Broadleaf Tree
!     7 - Needleleaf Tree
!     8 - C3 Grass
!     9 - C4 Grass
!    10 - Shrub
!--- for testing these tiles are set = 0
!    11 - 0 
!    11 - 0
!    11 - 0
!--- original non-veg tiles moved to these indices
!     14 - Urban
!     15 - Water
!     16 - Soil
!     17 - Ice


!     Missing number indicators
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Diagnostics info
      Real                                                              &
     & STASHwork1(*)                                                    &
                         ! STASH workspace
     &,STASHwork2(*)     ! STASH workspace


! Data arrays
      Real                                                              &
     &  p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, p_layer_centres(row_length, rows, 0:model_levels)               &
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = pstar, and at top = 0.
     &, p_star(row_length, rows)                                        &
     &, p_theta_levels(1-off_x:row_length+off_x,                        &
     &                 1-off_y:rows+off_y, model_levels)                &
     &, p(1-off_x:row_length+off_x,                                     &
     &    1-off_y:rows+off_y, model_levels)                             &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                     1-off_y:rows+off_y, model_levels)

      Real                                                              &
     &  cca (row_length, rows, N_cca_levels)                            &
     &, cclwp(row_length, rows)                                         &
                                ! condensed water path (KG/M**2)
     &, area_cloud_fraction(row_length, rows, wet_model_levels)

      Integer                                                           &
     &  ccb (row_length, rows)                                          &
     &, cct (row_length, rows)

! ancillary arrays and fields required to be saved from timestep to
! timestep.

      Real                                                              &
     &  T_surf(row_length, rows)                                        &
     &, TSTAR_LAND(ROW_LENGTH,ROWS)                                     &
                                    ! IN Land mean sfc temperature (K)
     &, TSTAR_SEA(ROW_LENGTH,ROWS)                                      &
                                    ! IN Open sea sfc temperature (K).
     &, TSTAR_SICE(ROW_LENGTH,ROWS) ! IN Sea-ice sfc temperature (K).


      logical                                                           &
     &  land_sea_mask(row_length, rows)                                 &
     &, land0p5(row_length, rows)                                       &
!         A mask set to .TRUE. if the fraction of land in the grid-box
!         exceeds 0.5.
     &, RAD_MASK(row_length, rows)
!         A mask which ensures a chequerboard pattern of radiation
!         calculations over the whole domain (not just one PE)


      Real                                                              &
     &  fland(land_field)
!         Fractional amount of land at each land point

      Real                                                              &
     &  Ntot_land                                                       &
                             ! Number of droplets over land / m-3
     &, Ntot_sea             ! Number of droplets over sea / m-3

      Real                                                              &
     &  snow_depth (row_length, rows)                                   &
                                      ! snow/qrclim.snow.(month)
     &, snow_depth_sea (row_length, rows) !  snow depth on sea ice

      Real                                                              &
     &  ice_fract (row_length, rows)                                    &
                                     ! ice/qrclim.ice.(month)
     &, soot (row_length, rows)

      Real                                                              &
     &  albsoil(land_field)                                             &
     &, lai(land_field, npft)                                           &
     &, rgrain(land_field, ntiles)                                      &
     &, snow_tile(land_field, ntiles)                                   &
     &, frac(land_field, ntype)                                         &
     &, frac_control(land_field, ntype)                                 &
     &, tstar_tile(land_field, ntiles)                                  &
     &, z0_tile(land_field, ntiles)                                     &
     &, dOLR_rts(row_length, rows)                                      &
                                        ! TOA - surface upward LW
     &, LW_down(row_length, rows)                                       &
                                        ! Surface downward LW
     &, SW_tile_rts(land_field, ntiles)                                 &
                                        ! Surface net SW on land tiles
     &, land_alb(row_length, rows)                                      &
                                        ! Mean land albedo
     &, sice_alb(row_length, rows)                                      &
                                        ! Mean sea-ice albedo
     &, ES_SPACE_INTERP(4, row_length, rows)                            &
!              Coeffs for spatial interpolation of radiation quantities
!     EAK
     &, alb_tile(land_field,ntiles,4)                                   &
     &, surf_down_sw(row_length,rows,4)

     LOGICAL                                &
     & l_cable                              
                       ! true if cable albedo is used
!     & ,L_TILE_PTS(land_field,ntiles)                              

     Integer istep_cur

     INTEGER       &
!       LAND_PTS    & ! IN No of land points being processed.
       SM_LEVELS   & ! IN No. of soil moisture levels
      ,TILE_INDEX(land_field,NTILES)  &  ! OUT Index of tile points
      ,TILE_PTS(NTILES)               &  ! OUT Number of tile points
      ,ISNOW_FLG3L(land_field,NTILES) &
      ,day
      REAL                                  &
       LAND_ALBEDO(ROW_LENGTH,ROWS,4)         &! Land albedo calculated by Cable
      ,SNOW_RHO1L(land_field,NTILES)          &! snow cover in the ice tile.
      ,TSOIL_TILE(land_field,NTILES,SM_LEVELS)&! Mean snow density  (or 1 layer)
      ,SNOW_TMP3L(land_field,NTILES,3)         ! Snow temperature (3 layer)



! Aerosol climatology for NWP

! arcl_dim.h
!
! Maximum dimensions for the aerosol climatology for NWP
!

      integer, parameter :: NPD_ARCL_SPECIES = 7
      integer, parameter :: NPD_ARCL_COMPNTS = 20

! end of arcl_dim.h
      
      ! Number of requested species within the climatology
      Integer n_arcl_species
      
      ! Corresponding number of requested components
      Integer n_arcl_compnts
      
      ! Model switch for each species
      Logical L_USE_ARCL(NPD_ARCL_SPECIES)
      
      ! Array indices of components
      Integer i_arcl_compnts(NPD_ARCL_COMPNTS)
      
      ! Mass-mixing ratios
      Real                                                              &
     &  arcl(row_length, rows, model_levels, n_arcl_compnts)

      Integer                                                           &
     &  land_index(land_field)

      Real                                                              &
     &  ozone(row_length, rows, ozone_levels)                           &
     &, O3_trop_level(row_length,rows)                                  &
     &, O3_trop_height(row_length,rows)                                 &
     &, T_trop_level(row_length,rows)                                   &
     &, T_trop_height(row_length,rows)                                  &
     &, u_1(salt_dim1, salt_dim2)                                       &
     &, v_1(salt_dim1, salt_dim2)                                       &
     &, SW_incs(row_length, rows, 0:model_levels+1)                     &
     &, LW_incs(row_length, rows, 0:model_levels)                       &
     &, dirpar_inc(row_length, rows)                                    &
     &, zh(row_length, rows)                                            &
                             ! boundary layer height
     &, co2_3D(row_length, rows, model_levels)                          &
     &, ch4_stochem(row_length, rows, model_levels)                     &
     &, cos_zenith_angle(row_length, rows)
! chemical greenhouse gas fields
      INTEGER, INTENT(IN) :: ngrgas
      REAL, INTENT(IN) ::                                               &
     &        grgas_field(row_length,rows,model_levels,ngrgas)

! Co-ordinate arrays
      Real                                                              &
           ! local vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &          model_levels)                                           &
                               ! Density*radius^2
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, delta_lambda                                                    &
     &, delta_phi

! Trig arrays
      real                                                              &
     &  cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)

! Grid-dependent arrays
      Real                                                              &
     &  f3_at_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y)          &
     &, true_longitude(row_length, rows)                                &
     &, true_latitude(row_length,rows)                                  &
     &, FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)

! arguments with intent in. ie: input variables.


      Real   , intent(in) :: ccw   (row_length, rows, wet_model_levels)
      Integer, intent(in) :: lcbase(row_length, rows)


! time information for current timestep
      Integer                                                           &
     &  val_year                                                        &
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second                                                      &
     &, timestep_number                                                 &
     &,         PREVIOUS_TIME(7)

      Real                                                              &
     &  recip_timestep

! Diagnostic variables
      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP

! Additional variables for SCM diagnostics which are dummy in full UM
      Integer                                                           &
     &  nSCMDpkgs          ! No of diagnostics packages

      Logical                                                           &
                           ! Logicals for diagnostics packages
     &  L_SCMDiags(nSCMDpkgs)

      Integer                                                           &
     &   maxsects

      Character*3                                                       &
     &   h_sect(0:maxsects)

      Integer                                                           &
     &  Error_code

! Variables with intent (in/out)

      Real                                                              &
     &  T_n(row_length, rows, model_levels)                             &
     &, q_n(row_length, rows, wet_model_levels)                         &
     &, qcl_n(row_length, rows, wet_model_levels)                       &
     &, qcf_n(row_length, rows, wet_model_levels)                       &
     &, qcf2_n(row_length, rows, wet_model_levels)                      &
                                                     ! 2nd ice prog
     &, qrain_n(row_length, rows, wet_model_levels)                     &
                                                     ! Rain prognostic
     &, qgraup_n(row_length, rows, wet_model_levels)                    &
                                                     ! Graupel
     &, cf_n(row_length, rows, wet_model_levels)                        &
     &, cfl_n(row_length, rows, wet_model_levels)                       &
     &, cff_n(row_length, rows, wet_model_levels)                       &
     &, T_inc(row_length, rows, model_levels)                           &
     &, q_inc(row_length, rows, wet_model_levels)                       &
     &, qcl_inc(row_length, rows, wet_model_levels)                     &
     &, cf_inc(row_length, rows, wet_model_levels)                      &
     &, cfl_inc(row_length, rows, wet_model_levels)                     &
     &, sea_salt_film(salt_dim1, salt_dim2, salt_dim3)                  &
     &, sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)                   &
     &, height(salt_dim1, salt_dim2, salt_dim3)                         &
     &, sum_eng_fluxes(row_length, rows)

      Real, Intent(InOut) ::                                            &
     &  so4_aitken(1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &                                               model_levels)      &
     &, so4_accu(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, so4_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, soot_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, soot_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, soot_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, bmass_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &                                               model_levels)      &
     &, bmass_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &                                               model_levels)      &
     &, bmass_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &                                               model_levels)      &
     &, ocff_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, ocff_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, ocff_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                               model_levels)      &
     &, DUST_DIV1(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, DUST_DIV2(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, DUST_DIV3(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, DUST_DIV4(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, DUST_DIV5(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, DUST_DIV6(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &                                              MODEL_LEVELS)       &
     &, biogenic(row_length, rows, model_levels)                        &
     &, aerosol(row_length, rows, model_levels)                         
!
! Variables with intent out
      Real                                                              &
     &  photosynth_act_rad(row_length, rows)                            &
                                             ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
     &, surf_radflux(row_length, rows)                                  &
                                       ! Net downward
     &, rad_hr(row_length, rows, bl_levels, 2)                          &
                                               !
!                               ! BL radiative (LW,SW) heating rates
     &, dOLR(row_length, rows)                                          &
                                    ! TOA - surface upward LW
     &, SW_tile(land_field, ntiles) ! Surface net SW on land tiles

!
! MOSES_II always true from UM6.3 onwards - should really be removed
! completely but non-trivial due to L_DUMMY arguments
!
       L_MOSES_II = .TRUE.

!
! Call the radiation code
!

! DEPENDS ON: glue_rad
      Call Glue_Rad (                                                   &

! Parallel variables
     &    halo_i, halo_j, off_x, off_y, global_row_length, global_rows  &
     &,   proc_row_group, proc_col_group, at_extremity, n_proc, n_procx &
     &,   n_procy, neighbour, g_rows, g_row_length, g_datastart, me     &

! model dimensions.
     &,   row_length, rows, n_rows                                      &
     &,   model_levels, wet_model_levels, bl_levels                     &
     &,   Ozone_levels, cloud_levels, N_cca_levels                      &
     &,   NTILES, LAND_FIELD, DUST_DIM1, DUST_DIM2, biogenic_dim1       &
     &,   biogenic_dim2, sulp_dim1, sulp_dim2, soot_dim1, soot_dim2     &
     &,   bmass_dim1, bmass_dim2, salt_dim1, salt_dim2, salt_dim3       &
     &,   co2_dim_len, co2_dim_row, co2_dim2, arcl_dim1, arcl_dim2      &
     &,   n_arcl_species, n_arcl_compnts, i_arcl_compnts                &
     &,   ocff_dim1, ocff_dim2                                          &

! Model switches
     &,   model_domain                                                  &
     &,   L_Rad_Step, L_Rad_Step_diag, L_Rad_Step_prog                  &
     &,   L_Forcing, L_Timestep, L_Radiance, L_Wenyi                    &
     &,   L_CAL360, L_SEC_VAR, L_EqT                                    &
     &,   L_MICROPHY, L_INHOM_CLOUD, L_USE_DUST, L_USE_BIOGENIC         &
     &,   L_use_sulpc_direct, L_use_soot_direct, L_use_soot_indirect    &
     &,   L_use_bmass_direct, L_use_bmass_indirect                      &
     &,   L_use_ocff_direct, L_use_ocff_indirect                        &
     &,   L_use_sulpc_indirect_SW, L_use_sulpc_indirect_LW              &
     &,   L_emcorr, L_climat_aerosol, L_clim_aero_hgt                   &
     &,   L_HadGEM1_Clim_Aero, Ltimer, L_ssice_albedo                   &
     !kdcorbin, 08/10 - changed l_cable flag
     &,   L_moses_ii, L_snow_albedo                                        &
     &,   L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a            &
     &,   L_use_seasalt_indirect                                        &
     &,   L_use_seasalt_direct                                          &
     &,   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                         &
     &,   L_pc2, l_mixing_ratio                                         &
     &,   L_MURK_RAD                                                    &
     &,   L_rad_deg                                                     &
!kdcorbin, 06/10 - added l_co2_radiation flag
     &,   L_co2_interactive, l_co2_radiation                            &
     &,   L_use_stochem_CH4, L_use_clearrh, L_ukca, L_USE_ARCL          &
     &,   L_USE_SPEC_SEA, L_MOD_BARKER_ALBEDO, L_MOD_K_FLUX             &
     &,   L_SCVARY, L_VOLCTS                                            &
     &,   CAN_RAD_MOD                                                   &

! model Parameters
     &,   A_SW_segments, A_SW_seg_size, A_LW_segments, A_LW_seg_size    &
     &,   INHOM_CLOUD_SW, INHOM_CLOUD_LW, DP_CORR_STRAT, DP_CORR_CONV   &
     &,   CO2_MMR, O2MMR, N2O_mix_ratio, CH4_mix_ratio                  &
     &,   CFC11_mix_ratio, CFC12_mix_ratio                              &
     &,   C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR                     &
     &,   SW_alpham, SW_alphac, SW_alphab, SW_dtice                     &
     &,   dt_bare,dalb_bare_wet,pen_rad_frac,SW_beta                    &
     &,   min_trop_level, max_trop_level                                &
     &,   Ntot_land, Ntot_sea, aero_bl_levels                           &

! Physical constants
     &,   lc, lf, cp, two_Omega, p_zero, kappa                          &
     &,   R, g, Lapse_Rate, earth_radius, Pi                            &

! in coordinate information
     &,   rho_r2, r_rho_levels, r_theta_levels                          &
     &,   eta_theta_levels, eta_rho_levels                              &
     &,   delta_lambda, delta_phi                                       &
     &,   lat_rot_NP, long_rot_NP                                       &

! in time stepping information.

     &,   timestep, radiation_timestep                                  &
     &,   radiation_tstep_diag, radiation_tstep_prog                    &
     &,   val_year, val_day_number, val_hour, val_minute                &
     &,   val_second, timestep_number                                   &
     &,   PREVIOUS_TIME                                                 &

! trig arrays
     &,   sin_theta_longitude, cos_theta_longitude                      &
     &,   FV_cos_theta_latitude                                         &

! grid-dependent arrays
     &,   f3_at_u, true_longitude, true_latitude                        &

! diagnostic info
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
! ARGSTS end
     &   STASHwork1                                                     &
     &,  STASHwork2                                                     &
!
! SCM diagnostics switches (dummy in full um)
     &,   nSCMDpkgs,L_SCMDiags                                          &

! in data fields.
     &,   p_star                                                        &
     &,   p_layer_boundaries, p_layer_centres                           &
     &,   p, p_theta_levels                                             &
     &,   exner_rho_levels, exner_theta_levels                          &
     &,   land_sea_mask, fland,land0p5,l_ctile                          &
     &,   T_SURF,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE , AREA_CLOUD_FRACTION  &
     &,   DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6   &
     &,   biogenic, so4_aitken, so4_accu, so4_diss                      &
     &,   soot_new, soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld &
     &,   ocff_new, ocff_agd, ocff_cld, aerosol, arcl                   &
     &,   co2_3D                                                        &
     &,   ch4_stochem                                                   &
     &,   frac_control                                                  &
! chemical greenhouse gas fields
     &, ngrgas, grgas_field                                             &

! ancillary fields and fields needed to be kept from timestep to
! timestep
     &,   snow_depth, snow_depth_sea, ice_fract, rgrain, soot           &
     &,   cca, ccb, cct, cclwp, ccw, lcbase                             &
     &,   ozone, SW_incs, LW_incs,dirpar_inc                            &
     &,   O3_trop_level, O3_trop_height                                 &
     &,   T_trop_level, T_trop_height, zh                               &
     &,   land_index, albsoil, lai, snow_tile                           &
     &,   frac, tstar_tile, z0_tile                                     &
     &,   dOLR_rts, LW_down, SW_tile_rts                                &
     &,   u_1, v_1                                                      &
     &,   land_alb,sice_alb                                             &
     &,   ES_SPACE_INTERP, A_SW_radstep, A_LW_radstep                   &
     &,   A_SW_radstep_diag, A_SW_radstep_prog                          &
     &,   A_LW_radstep_diag, A_LW_radstep_prog, RAD_MASK                &

! in/out
     &,   T_n, q_n, qcl_n, qcf_n, cf_n, cfl_n, cff_n                    &
     &,   qcf2_n, qrain_n, qgraup_n                                     &
     &,   T_inc, q_inc, qcl_inc, cf_inc, cfl_inc                        &
     &,   sum_eng_fluxes                                                &
     &,   cos_zenith_angle                                              &

! out.
     &,   photosynth_act_rad, surf_radflux, rad_hr, dOLR, SW_tile       &
!EAK
     &, l_cable                                                         &
     &, istep_cur                                                       &
!     &, L_TILE_PTS,TILE_PTS,SM_LEVELS,TILE_INDEX           &
     &, TILE_PTS,SM_LEVELS,TILE_INDEX           &
     &,   surf_down_sw,alb_tile                                         &
     &, SNOW_TMP3L,SNOW_RHO1L,TSOIL_TILE,ISNOW_FLG3L,LAND_ALBEDO        &
! Section Information
     &, maxsects, h_sect                                                &
! error information
     &,   Error_code  )

!


      return
      END SUBROUTINE NI_RAD_CTL
