
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interface to Atmospheric Physics radiation code.
!
! Purpose:
!   This is the top-level radiation control routine. Major book-keeping
!   is carried out here. Separate calls are made to calculate SW and LW 
!   radiative fluxes and heating rates on radiation time-steps.
!   Radiation increments are applied on all physics time-steps.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      Subroutine Glue_Rad (                                             &

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
     &, model_domain, L_Rad_Step                                        &
     &, L_Rad_Step_Diag, L_Rad_Step_prog                                &
     &, L_Forcing, L_Timestep, L_Radiance, L_Wenyi                      &
     &, L_CAL360, L_SEC_VAR, L_Eqt                                      &
     &, L_MICROPHY, L_INHOM_CLOUD, L_USE_DUST, L_USE_BIOGENIC           &
     &, L_use_sulpc_direct, L_use_soot_direct, L_use_soot_indirect      &
     &, L_use_bmass_direct, L_use_bmass_indirect                        &
     &, L_use_ocff_direct, L_use_ocff_indirect                          &
     &, L_use_sulpc_indirect_SW, L_use_sulpc_indirect_LW                &
     &, L_emcorr, L_climat_aerosol,L_clim_aero_hgt,L_HadGEM1_Clim_Aero  &
     &, Ltimer, L_ssice_albedo, L_MOSES_II, L_snow_albedo               &
     &, L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a              &
     &, L_use_seasalt_indirect                                          &
     &, L_use_seasalt_direct                                            &
     &, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                           &
     &, L_pc2, l_mixing_ratio                                           &
     &, L_MURK_RAD                                                      &
     &, L_rad_deg                                                       &
! kdcorbin, 06/10 - added l_co2_radiation flag
     &, L_co2_interactive, L_co2_radiation                              &
     &, L_use_stochem_CH4, L_use_clearrh, L_ukca, L_USE_ARCL            &
     &, L_USE_SPEC_SEA, L_MOD_BARKER_ALBEDO, L_MOD_K_FLUX               &
     &, L_SCVARY, L_VOLCTS                                              &
     &, can_rad_mod                                                     &

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
     &, radiation_tstep_diag                                            &
     &, radiation_tstep_prog                                            &
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
!
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
     &, frac_control                                                    &
! chemical greenhouse gas fields
     &, ngrgas, grgas_field                                             &

! ancillary fields and fields needed to be kept from timestep to
! timestep
     &, snow_depth, snow_depth_sea, ice_fract, rgrain, soot             &
     &, cca, ccb, cct, cclwp, ccw, lcbase                               &
     &, ozone, SW_incs, LW_incs, dirpar_inc                             &
     &, O3_trop_level, O3_trop_height                                   &
     &, T_trop_level, T_trop_height, zh, land_index, albsoil, lai       &
     &, snow_tile, frac, tstar_tile, z0_tile                            &
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
     &, cos_zenith_angle                                                &

! out.
     &, photosynth_act_rad, surf_radflux, rad_hr, dOLR, SW_tile         &
! EAK
     &, l_cable                                                         &
     &, istep_cur                                                       &
!     &, L_TILE_PTS,TILE_PTS,SM_LEVELS,TILE_INDEX           &
     &, TILE_PTS,SM_LEVELS,TILE_INDEX           &
     &, surf_down_sw,alb_tile                                           &
     &, SNOW_TMP3L,SNOW_RHO1L,TSOIL_TILE,ISNOW_FLG3L,LAND_ALBEDO        &
! Section Information
     &, maxsects, h_sect                                                &
! error information
     &, Error_code  )


! Modules required for the orography scheme

      Use solinc_data, Only: sol_bearing, orog_corr, f_orog,            &
     &                       slope_aspect, slope_angle, L_orog

      Use swrdiag_mod, Only: StrSWDiag
      Use lwrdiag_mod, Only: StrLWDiag
      Use volcts_Mod
      Use csenario_mod
!dhb599 20110812: get SC
!      REAL SC                           !  Solar constant
!      PARAMETER ( SC = 1365. )
!20110812: define SC in coupling control
!
! Now get SC from namelist of coupling control
       use auscom_cpl_data_mod, only: SC,VOLCTS_val
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

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

! Parameters
!dhb599 20110812: now SC is passed in by "swsc.h" (see above)
!      Real SC                           !  solar constant
!      Parameter ( SC = 1365. )

      Integer icode   
      Parameter (icode = 0 )
      Character*80 cmessage

! Model dimensions
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
     &, BIOGENIC_DIM1                                                   &
                      !dimensions of biogenic aerosol arrays
     &, BIOGENIC_DIM2                                                   &
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
                                ! dimensions of OCFF arrays
     &, salt_dim1                                                       &
                                ! dimensions of sea-salt arrays
     &, salt_dim2                                                       &
     &, salt_dim3                                                       &
     &, co2_dim1                                                        &
                               ! dimensions of CO2 array
     &, co2_dim2                                                        &
     &, co2_dim_len                                                     &
     &, co2_dim_row                                                     &
                               ! dimensions of aerosol clim for NWP
     &, arcl_dim1                                                       &
     &, arcl_dim2

! Model switches
      INTEGER                                                           &
     & can_rad_mod                                                      
                                         ! Which canopy radiation model
                                         ! we're using
      Integer                                                           &
     & model_domain

      Logical                                                           &
     &  L_Rad_Step                                                      &
                        ! true if a radiation timestep
     &, L_rad_step_diag                                                 &
                             ! Dummy variable in version 3A
     &, L_rad_step_prog                                                 &
                             ! Dummy variable in version 3A
     &, L_Forcing                                                       &
                             ! Dummy variable in version 3A
     &, L_Timestep                                                      &
                             ! Dummy variable in version 3A
     &, L_Radiance                                                      &
                             ! Dummy variable in version 3A
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
!                            ! climatological aerosols in FILL3A
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
     & ,L_sice_scattering                                               &
                          ! true if sea ice internal scattering scheme
     & ,L_sice_hadgem1a                                                 &
                        ! true if HadGEM1 albedo bug corrected
     &, L_MOSES_II                                                      &
                       ! True if 8A boundary layer selected
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
     &, L_co2_radiation                                                 &
          ! Controls the CO2/radiation interaction - kdcorbin, 06/10
     &, L_use_stochem_CH4                                               &
                           ! Control the use of STOCHEM CH4 field
     &, L_use_clearrh                                                   &
                      ! RH used for hygroscopic aerosols
     &, L_ukca                                                          &
                ! Switch for UKCA sub-model
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
     &, radiation_tstep_prog                                            &
                             ! Dummy variable in version 3A
     &, radiation_tstep_diag                                            &
                             ! Dummy variable in version 3A
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
     &, max_trop_level  ! Highest permitted level for the tropopause
!                       ! used for radiative purposes.

      Integer                                                           &
     &  A_SW_segments                                                   &
     &, A_SW_seg_size                                                   &
     &, A_LW_segments                                                   &
     &, A_LW_seg_size                                                   &
     &, A_SW_radstep                                                    &
     &, A_LW_radstep                                                    &
     &, A_SW_radstep_diag                                               &
                             ! Dummy variable in version 3A
     &, A_LW_radstep_diag                                               &
                             ! Dummy variable in version 3A
     &, A_SW_radstep_prog                                               &
                             ! Dummy variable in version 3A
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


      Real    :: ccw   (row_length, rows, wet_model_levels)
      Integer :: lcbase(row_length, rows)

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
     &, frac_control(land_field,ntype)                                  &
                                       ! Dummy variable in version 3A
     &, tstar_tile(land_field, ntiles)                                  &
     &, z0_tile(land_field, ntiles)                                     &
     &, dOLR_rts(row_length, rows)                                      &
                                        ! TOA - surface upward LW
     &, LW_down(row_length, rows)                                       &
                                        ! Surface downward LW
     &, SW_tile_rts(land_field, ntiles)                                 &
                                        ! Surface net SW on land tiles
     &, SW_tile_down_rts(land_field, ntiles)                            &
     &, land_alb(row_length, rows)                                      &
                                        ! Mean land albedo
     &, sice_alb(row_length, rows)                                      &
                                        ! Mean sea-ice albedo
     &, ES_SPACE_INTERP(4, row_length, rows)
!              Coeffs for spatial interpolation of radiation quantities

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
!
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
     &, ch4_stochem(row_length, rows, model_levels)

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
     &, true_latitude(row_length, rows)                                 &
     &, FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)


! Allocatable arrays for diagnostic variables - required to save memory
! use during radiation routines
      Real,Dimension(:,:,:),Allocatable::                               &
     & T_incr_diagnostic     ! temperature increment for STASH


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

!  Additional variables for SCM diagnostics which are dummy in full UM
      Integer                                                           &
     &  nSCMDpkgs          ! No of SCM diagnostics packages

      Logical                                                           &
                           ! Logicals for SCM diagnostics packages
     &  L_SCMDiags(nSCMDpkgs)

      Integer                                                           &
     &  maxsects

      Character*3                                                       &
     &  h_sect(0:maxsects)

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
     &, aerosol(row_length, rows, model_levels)       

      Real, Intent(In) ::                                               &
     &  biogenic(row_length, rows, model_levels)
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
     &, SW_tile(land_field, ntiles)                                     &
                                    ! Surface net SW on land tiles
     &, SW_tile_down(land_field, ntiles)                                &
                                    ! Surface down SW on land tiles
     &, cos_zenith_angle(row_length,rows)                               

! local variables.

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end

      Character*(*), Parameter ::  RoutineName = 'glue_rad'

! loop counters
      Integer                                                           &
     &  i, j, k, info                                                   &
     &, l, n, point

      Integer JS 
      Integer IJ

! local variables
      Integer                                                           &
     &  daylight_points                                                 &
     &, N_SW_segments                                                   &
     &, step                                                            &
     &, step_sw                                                         &
     &, seg_start                                                       &
     &, lit_points                                                      &
     &, start_point                                                     &
     &, seg_points                                                      &
     &, first_point                                                     &
     &, ptr_local                                                       &
                      ! Pointer for LW arrays used to reduce length
!                     ! of lines, repeating current element of
!                     ! first_point_local
     &, List_LW_points(row_length*rows)                                 &
     &, LW_POINTS, num_lw_points                                        &
!         Variables to enable scatter/gather of LW (like SW): creates
!         a LIST of points to do a calculation, and also
!              facilitates segmenting.
     &, INTERP_INDEX                                                    &
!              Variable which alternates between 0 and 1 every other
!              radiation timestep
     &, FIRST_DATA_INTERP                                               &
!              The first data point, in data co-ords, which needs to be
!              interpolated on a PE (for use in interpolation routine)
     &, first_data_interp_sw                                            &
!              Saved value to use with ISCCP diagnostics
     &, first_row, last_row                                             &
     &, tot_daylight_points
!              Total number of daylight points in whole domain

    
! LW segmentation allocatable arrays
      Integer, Allocatable :: first_point_local(:)
      Integer, Allocatable :: seg_points_local(:)

! SW segmentation allocatable arrays
      Integer, Allocatable :: first_point_temp(:)
      Integer, Allocatable :: first_point_temp_sw(:)
      Integer, Allocatable :: seg_points_temp(:)
      Integer, Allocatable :: seg_points_temp_sw(:)

! Local arrays holding information to be passed between physics
! routines.
      Real                                                              &
     &  T_large(1-off_x:row_length+off_x,                               &
     &          1-off_y:rows+off_y, model_levels)                       &
!
!
! The following arrays are required by the PC2 cloud scheme
!
     &,T_latest(row_length,rows,model_levels)                           &
     &,q_latest(row_length,rows,wet_model_levels)                       &
     &,qcl_latest(row_length,rows,wet_model_levels)                     &
     &,cf_latest(row_length,rows,wet_model_levels)                      &
     &,cfl_latest(row_length,rows,wet_model_levels)                     &
     &,cff_latest(row_length,rows,wet_model_levels)                     &
     &,delta_T(row_length,rows,model_levels)                            &
     &,zeros(row_length,rows,wet_model_levels)                          &
     &,work_3d(row_length,rows,model_levels)
!
! For calculation of masses of levels - these arrays will not have
! halo information
!
      Real                                                              &
     & r_rho_levels_nh(row_length*rows,model_levels)                    &
     &,r_theta_levels_nh(row_length*rows,0:model_levels)                &
     &,rho_r2_nh(row_length*rows,model_levels)

! Diagnostics controlled by Diagnostic switches

! Radiation fields 1. SW & common with LW.
!
!
      Real                                                              &
     &  day_fraction(row_length,rows)                                   &
     &, volcmass(row_length,rows)      					&
! VOLCTS expanded to full fields,
!     and converted from optical depth to mass loading.
     &, mean_cos_zenith_angle(row_length,rows)                          &
     &, flandg(row_length, rows)                                        &
                                  ! Land fractions at all grid-points
     &, land_albedo(row_length, rows, 4)                                &
     &, land_and_ice_albedo(row_length, rows)                           &
     &, open_sea_albedo(row_length, rows, 2)                            &
     &, netSW(row_length, rows)                                         &
                                  ! Net short-wave absorbed by planet
     &, SWsea(row_length, rows)                                         &
                                  ! Net short-wave absorbed by planet
     &, DUST_1(DUST_DIM1,DUST_DIM2)                                     &
     &, DUST_2(DUST_DIM1,DUST_DIM2)                                     &
     &, DUST_3(DUST_DIM1,DUST_DIM2)                                     &
     &, DUST_4(DUST_DIM1,DUST_DIM2)                                     &
     &, DUST_5(DUST_DIM1,DUST_DIM2)                                     &
     &, DUST_6(DUST_DIM1,DUST_DIM2)                                     &
     &, local_biogenic(BIOGENIC_DIM1,BIOGENIC_DIM2)                     &
     &, accum_sulphate(sulp_dim1,sulp_dim2)                             &
     &, aitken_sulphate(sulp_dim1,sulp_dim2)                            &
     &, diss_sulphate(sulp_dim1,sulp_dim2)                              &
     &, fresh_soot(soot_dim1,soot_dim2),aged_soot(soot_dim1,soot_dim2)  &
     &, fresh_bmass(bmass_dim1,bmass_dim2)                              &
     &, aged_bmass(bmass_dim1,bmass_dim2)                               &
     &, cloud_bmass(bmass_dim1,bmass_dim2)                              &
     &, fresh_ocff(ocff_dim1, ocff_dim2)                                &
     &, aged_ocff(ocff_dim1, ocff_dim2)                                 &
     &, cloud_ocff(ocff_dim1, ocff_dim2)                                &
     &, sin_true_latitude(row_length, rows)                             &
     &, height_rho(row_length, rows, model_levels)                      &
     &, height_theta(row_length, rows, 0:model_levels)                  &
     &, sw_net_land(row_length, rows)                                   &
                                       !SW net local flux over land
     &, sw_net_sice(row_length, rows)  !SW net local flux over sea-ice
!
      ! Aerosol climatology for NWP
      Real                                                              &
     &  local_arcl(arcl_dim1, arcl_dim2, n_arcl_compnts)
      Real                                                              &
     &  arcl_multf ! Mass-mixing ratio conversion factor
!
      Integer                                                           &
     &  List_daylight_points(row_length*rows)                           &
     &, diag_row_list(row_length*rows)                                  &
                                       ! List of row indices of points
!                                      ! where diagnostics are
!                                      ! calculated
     &, diag_row_list_sw(row_length*rows)                               &
                                           ! Saved version for ISCCP
     &, diag_col_list(row_length*rows)                                  &
                                       ! List of column indices of
!                                      ! points where diagnostics are
!                                      ! calculated
     &, diag_col_list_sw(row_length*rows)                               &
                                           ! Saved version for ISCCP
     &, trindx(row_length, rows)

      Logical                                                           &
     &  switch(row_length,rows)

      Real                                                              &
     &  Sindec                                                          &
                         !  sin(solar declination)
     &, hour_angle(row_length,rows)                                     &
                         !  Hour Angle
     &, SCS                                                             &
                         !  solar constant scaling factor
     &, seconds_since_midnight                                          &
     &, Eq_Time         ! The Equation of time

      Integer                                                           &
     &  global_cloud_top                                                &
     &, FIRST_POINT_DUST                                                &
     &, FIRST_POINT_BIOGENIC                                            &
     &, first_point_sulpc                                               &
     &, first_point_soot                                                &
     &, first_point_biomass                                             &
     &, first_point_ocff                                                &
     &, first_point_seasalt                                             &
     &, first_point_arcl                                                &
     &, nd_rad_pts                                                      &
     &, n_rad_layers
!
      Logical :: L_complete_North
!                  Flag to complete field on the northern polar row
      Logical :: L_complete_South
!                  Flag to complete field on the southern polar row
      Logical :: L_complete_deg
!                  Flag to complete field because of degradation
!
      Logical                                                           &
     &  L_CO2_3D    !  controls use of 3D co2 field
!
!     Declaration of SW diagnostics.
!
      Type (StrSWDiag) :: SW_diag
!
      Logical :: L_T_incr_sw    ! Flag for SW temperature increments
!
! SW diagnostics not on stash flag
      Real                                                              &
     &  itoasw(row_length,rows)                                         &
     &, surfsw(row_length,rows)

! SW diagnostics on always true stash flag
      Logical                                                           &
     &  l_flux_below_690nm_surf

      Parameter (l_flux_below_690nm_surf = .true.)

      Real                                                              &
     &  flux_below_690nm_surf(row_length,rows)





! Radiation fields 2. LW

      Real                                                              &
     &  LWsea(row_length, rows)                                         &
                                  !
     &, OLR(row_length, rows)                                           &
                                !
     &, surflw(row_length,rows)                                         &
                                !
     &, top_absorption(row_length, rows)
! needed by energy correction
      Real                                                              &
     &  net_atm_flux (row_length, rows)
!
      Type (StrLWDiag) :: LW_diag
!
      Logical :: L_T_incr_lw    ! Flag for LW temperature increments
!
! needed for 8A boundary layer
      Real                                                              &
     &  alb_tile(land_field,ntiles,4)                                   &
     &, surf_down_sw(row_length,rows,4)                                 &
     &, T_sol_rad(row_length, rows)                                     &
                                     ! Effective radiative temperature
!                                    ! over solid surface of gridbox.
     &, tile_frac(land_field,ntiles)
      Integer                                                           &
     &  land_index_i(land_field)                                        &
     &, land_index_j(land_field)                                        &
     &, tile_pts(ntype)                                                 &
     &, tile_index(land_field,ntype)
!   Workspace
      Real                                                              &
     &  fracsolid(row_length, rows)                                     &
!                           !Solid surface fraction in gridbox
     &, sw_net_rts(row_length, rows)                                    &
!                           !Net surface SW on radiation timesteps.
     &, sw_down_rts(row_length, rows)
!                           !Downward surface SW on radiation timesteps.

      Real                                                              &
     &  albsolid
!                           !Mean solid surface albedo

      real Recip_Pi_Over_180, alat, alon
!            Latitude & longitude in degree
! EAK
     Integer istep_cur
     INTEGER       &
!       LAND_PTS    & !  land_filed used; (IN) No of land points being processed
       SM_LEVELS   & ! IN No. of soil moisture levels
!      ,TILE_INDEX(land_field,ntiles)  &  ! OUT Index of tile points
!      ,TILE_PTS(ntype)             &  ! OUT Number of tile points
      ,ISNOW_FLG3L(land_field,ntiles)
       Real                                                              &
!     &  alb_tile(land_field,ntiles,4)                                   &
!     &  surf_down_sw(row_length,rows,4)                                 &
!     &, cos_zenith_angle(row_length,rows)                               &
     & SNOW_RHO1L(land_field,ntiles)          &! snow cover in the ice tile.
!     &,TSOIL_TILE(land_field,ntiles,1)&! Mean snow density (or 1 layer)
     &,TSOIL_TILE(land_field,ntiles,SM_LEVELS)&! Mean snow density (or 1 layer)
     &,SNOW_TMP3L(land_field,ntiles,3)         ! Snow temperature (3 layer)

     LOGICAL                                &
     & l_cable                              
!     & ,L_TILE_PTS(land_field,ntiles)                              

! CSIGMA start
      ! Stefan-Boltzmann constant (W/m**2/K**4).
      REAL, PARAMETER ::  SBCON=5.67E-8
! CSIGMA end
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------

      Real rowvol         ! Variable to pass aerosol mass to VOLCMASS
      Real mascon         ! Conversion factor from optical depth to mass
      Parameter ( mascon =     						&
     &    .0001 * 1.16898 / ( 6.09715E+03 + 4.60178E-04 ))

!     .0001 converts the numbers in the DATA statement to optical depth.
!     1.16898 converts .55 micron optical depth to the Slingo (1989) 
!     band-1 optical depth (from the HadCM2 mod) & then divide by the
!     total k (scattering+absorption) for bands 2 & 3, originally from
!     spec3a_sw_3_tknlaero but still the same in spec3a_sw_hadgem1_2.

      Real delslt         ! "error" to allow when checking sin(latitude)
      Parameter ( delslt = 1.E-6 )


! Set max sizes for E-S radiation.
!     ------------------------------------------------------------------
!     MODULE SETTING MAXIMUM DIMENSIONS OF ARRAYS IN THE RADIATION CODE.
!   4.5   Aug 1998     Increment by 2 the no. of aerosol species
!                      affecting the radiation.    Luke Robinson
!
!
!   5.2   Dec 2000     Increment by 2 the no. of aerosol species
!                      affecting the radiation.    Andy Jones
!   5.3   Oct 2001     Allow six more types of aerosol in preparation
!                      for runs with dust.      (J. M. Edwards)
!   5.5   Feb 2003     Increase aerosol species to include dust
!                      and biomass    (S. Woodward)
!   5.5   Feb 2003     Increase the number of types of ice crystal.
!                                                 John Edwards
!   6.2   Nov 2005     Increase number of types of data.
!                      Add number of wavelengths for aerosol optical
!                      depth.
!                                                 N. Bellouin
!   6.2   Jan 2006     Increase values of the cloud parameters and
!                      of the number of terms in the phase function
!                      under versions 3C and 3Z of the radiation
!                      code.                         (J.-C. Thelen)
!
      INTEGER,PARAMETER :: NPD_TYPE      = 16 ! number of types of data
      INTEGER,PARAMETER :: NPD_BAND      = 20 ! number of spectral bands
      INTEGER,PARAMETER :: NPD_EXCLUDE   = 2  ! nubmer of excluded bands
      INTEGER,PARAMETER :: NPD_SPECIES   = 11 ! number of gaseous species
      INTEGER,PARAMETER :: NPD_ESFT_TERM = 16 ! number of esft terms
      INTEGER,PARAMETER :: NPD_SCALE_FNC = 4  ! number of scaling funcs

      ! number of scaling variables
      INTEGER,PARAMETER :: NPD_SCALE_VARIABLE=4

      INTEGER,PARAMETER :: NPD_SURFACE     = 1 ! no of surface types
      INTEGER,PARAMETER :: NPD_ALBEDO_PARM = 4 ! no of albedo parameters
      INTEGER,PARAMETER :: NPD_CONTINUUM   = 2 ! no of continua
      INTEGER,PARAMETER :: NPD_DROP_TYPE   = 6 ! no of drop types
      INTEGER,PARAMETER :: NPD_ICE_TYPE    = 8 ! no of ice crystal types


      ! max no of cloud parameters
      INTEGER,PARAMETER :: NPD_CLOUD_PARAMETER=30

      INTEGER, PARAMETER :: NPD_AEROSOL_SPECIES=24

      ! maximum no of humidities
      INTEGER,PARAMETER :: NPD_HUMIDITIES=21

      ! no of thermal coefficients
      INTEGER,PARAMETER :: NPD_THERMAL_COEFF=9

      ! number of wavelengths for aerosol optical depth
      INTEGER,PARAMETER :: NPD_AOD_WAVEL=6
! MXSIZE3A end
!     Declaration of spectral files for 3A-radiation.
! SWSPDL3A contains declarations for reduced sw-spectral file in
! two-stream radiation code.

!     note: SWSPDC3A, SWSPCM3A and SWSARG3A must be consistent
!     note: since the arrays here will be passed in a common block
!     their sizes must be fixed, even though variable sizes are used
!     lower in the code. they are accordingly defined as 1-dimensional
!     arrays with fixed maximum sizes at this level.

! dimensions for the SW spectrum

      INTEGER::  NPD_TYPE_SW            ! number of types of data
      INTEGER::  NPD_BAND_SW            ! number of spectral bands
      INTEGER::  NPD_EXCLUDE_SW         ! number of excluded bands
      INTEGER::  NPD_SPECIES_SW         ! number of gaseous species
      INTEGER::  NPD_ESFT_TERM_SW       ! number of esft terms
      INTEGER::  NPD_SCALE_FNC_SW       ! number of scaling functions
      INTEGER::  NPD_SCALE_VARIABLE_SW  ! number of scaling variables
      INTEGER::  NPD_SURFACE_SW         ! number of surface types
      INTEGER::  NPD_ALBEDO_PARM_SW     ! number of albedo parameters
      INTEGER::  NPD_CONTINUUM_SW       ! number of continua
      INTEGER::  NPD_DROP_TYPE_SW       ! number of drop types
      INTEGER::  NPD_ICE_TYPE_SW        ! number of ice crystal types
      INTEGER::  NPD_AEROSOL_SPECIES_SW ! number of aerosol species
      INTEGER::  NPD_CLOUD_PARAMETER_SW ! max number of cloud parameters
      INTEGER::  NPD_HUMIDITIES_SW      ! maximum number of humidities
      INTEGER::  NPD_THERMAL_COEFF_SW   ! number of thermal coefficients
      INTEGER::  NPD_AOD_WAVEL_SW       ! number of wavelengths for AOD

! general fields:

      ! flag for types of data present

      LOGICAL :: L_PRESENT_SW(0: NPD_TYPE)

! properties of the spectral bands:

      INTEGER :: N_BAND_SW  ! number of spectral bands

      REAL :: WAVE_LENGTH_SHORT_SW(NPD_BAND) ! shorter wavelength limits
      REAL :: WAVE_LENGTH_LONG_SW(NPD_BAND)  ! longer wavelength limits

! exclusion of specific bands from parts of the spectrum:

      ! number of excluded bands within each spectral band
      INTEGER :: N_BAND_EXCLUDE_SW(NPD_BAND)

      ! indices of excluded bands

      INTEGER :: INDEX_EXCLUDE_SW(NPD_EXCLUDE, NPD_BAND)

! fields for the solar flux:

      ! fraction of the incident solar flux in each band
      REAL :: SOLAR_FLUX_BAND_SW(NPD_BAND)

! fields for Rayleigh scattering:

      REAL :: RAYLEIGH_COEFFICIENT_SW(NPD_BAND) ! Rayleigh coefficients

! fields for gaseous absorption:

      INTEGER :: N_ABSORB_SW ! number of absorbers

      ! number of absorbers in each band
      INTEGER :: N_BAND_ABSORB_SW(NPD_BAND)

      ! list of absorbers in each band
      INTEGER :: INDEX_ABSORB_SW(NPD_SPECIES, NPD_BAND)

      ! types of each gas in the spectral file
      INTEGER :: TYPE_ABSORB_SW(NPD_SPECIES)

      ! number of esft terms in each band for each gas
      INTEGER :: I_BAND_ESFT_SW(NPD_BAND, NPD_SPECIES)

      ! type of esft scaling
      INTEGER :: I_SCALE_ESFT_SW(NPD_BAND, NPD_SPECIES)

      ! type of scaling function

      INTEGER :: I_SCALE_FNC_SW(NPD_BAND, NPD_SPECIES)

      !  esft exponents
      REAL :: K_ESFT_SW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)

      ! esft weights
      REAL :: W_ESFT_SW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)

      ! scaling parameters for each absorber and term

      REAL :: SCALE_VECTOR_SW(NPD_SCALE_VARIABLE, NPD_ESFT_TERM,        &
     &  NPD_BAND, NPD_SPECIES)

      ! reference pressure for scaling function
      REAL :: P_REFERENCE_SW(NPD_SPECIES, NPD_BAND)

      ! reference temperature for scaling function
      REAL :: T_REFERENCE_SW(NPD_SPECIES, NPD_BAND)

! representation of the Planckian:

      INTEGER :: N_DEG_FIT_SW ! degree of thermal polynomial

      ! coefficients in polynomial fit to source function
      REAL :: THERMAL_COEFFICIENT_SW(0: NPD_THERMAL_COEFF-1,NPD_BAND)

      REAL :: T_REF_PLANCK_SW ! Planckian reference temperature

! surface properties:

      ! method of specifying properties of surface

      INTEGER :: I_SPEC_SURFACE_SW(NPD_SURFACE)

      ! number of parameters fitting the direct albedo
      INTEGER :: N_DIR_ALBEDO_FIT_SW(NPD_SURFACE)


      LOGICAL :: L_SURFACE_SW(NPD_SURFACE) ! surface types included

      ! surface albedos
      REAL :: SURFACE_ALBEDO_SW(NPD_BAND, NPD_SURFACE)

      ! coefficients for fitting direct albedo

      REAL :: DIRECT_ALBEDO_PARM_SW(0: NPD_ALBEDO_PARM, NPD_BAND,       &
     &  NPD_SURFACE)

      ! surface emissivities
      REAL :: EMISSIVITY_GROUND_SW(NPD_BAND, NPD_SURFACE)

! fields for continua:

      ! number of continua in each band
      INTEGER :: N_BAND_CONTINUUM_SW(NPD_BAND)

      ! list of continua in each band
      INTEGER :: INDEX_CONTINUUM_SW(NPD_BAND, NPD_CONTINUUM)

      INTEGER :: INDEX_WATER_SW ! index of water vapour

      ! type of scaling function for continuum
      INTEGER :: I_SCALE_FNC_CONT_SW(NPD_BAND, NPD_CONTINUUM)

      ! grey extinction coefficients for continuum
      REAL :: K_CONTINUUM_SW(NPD_BAND, NPD_CONTINUUM)

      ! scaling parameters for continuum
      REAL :: SCALE_CONTINUUM_SW(NPD_SCALE_VARIABLE,NPD_BAND,           &
     &  NPD_CONTINUUM)

      ! reference pressure for scaling of continuum
      REAL :: P_REF_CONTINUUM_SW(NPD_CONTINUUM, NPD_BAND)

      ! reference temperature for scaling of continuum
      REAL :: T_REF_CONTINUUM_SW(NPD_CONTINUUM, NPD_BAND)

! fields for water droplets:

      ! parametrization type of droplets
      INTEGER :: I_DROP_PARAMETRIZATION_SW(NPD_DROP_TYPE)

      ! types of droplet present
      LOGICAL :: L_DROP_TYPE_SW(NPD_DROP_TYPE)


      ! parameters used to fit optical properties of clouds
      REAL :: DROP_PARAMETER_LIST_SW(NPD_CLOUD_PARAMETER, NPD_BAND,     &
     &  NPD_DROP_TYPE)

      ! minimum dimension permissible in the parametrization
      REAL ::DROP_PARM_MIN_DIM_SW(NPD_DROP_TYPE)

      ! maximum dimension permissible in the parametrization
      REAL :: DROP_PARM_MAX_DIM_SW(NPD_DROP_TYPE)

! fields for aerosols:

      INTEGER :: N_AEROSOL_SW ! number of species of aerosol

      ! types of aerosols
      INTEGER :: TYPE_AEROSOL_SW(NPD_AEROSOL_SPECIES)

      ! parametrization of aerosols
      INTEGER :: I_AEROSOL_PARAMETRIZATION_SW(NPD_AEROSOL_SPECIES)

      ! numbers of humidities
      INTEGER :: NHUMIDITY_SW(NPD_AEROSOL_SPECIES)

      ! aerosol species included
      LOGICAL :: L_AEROSOL_SPECIES_SW(NPD_AEROSOL_SPECIES)

      ! absorption by aerosols
      REAL :: AEROSOL_ABSORPTION_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,&
     &  NPD_BAND)

      ! scattering by aerosols
      REAL :: AEROSOL_SCATTERING_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,&
     &  NPD_BAND)

      ! asymmetry of aerosols
      REAL :: AEROSOL_ASYMMETRY_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES, &
     &  NPD_BAND)

      ! humidities for components
      REAL :: HUMIDITIES_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)

! fields for ice crystals:

      ! types of parametrization of ice crystals
      INTEGER :: I_ICE_PARAMETRIZATION_SW(NPD_ICE_TYPE)

      ! types of ice crystal present
      LOGICAL :: L_ICE_TYPE_SW(NPD_ICE_TYPE)

      ! parameters used to fit single scattering of ice crystals
      REAL :: ICE_PARAMETER_LIST_SW(NPD_CLOUD_PARAMETER, NPD_BAND,      &
     &  NPD_ICE_TYPE)

      ! minimum dimension permissible in the parametrization

      REAL :: ICE_PARM_MIN_DIM_SW(NPD_ICE_TYPE)

      ! maximum dimension permissible in the parametrization
      REAL :: ICE_PARM_MAX_DIM_SW(NPD_ICE_TYPE)

! fields for doppler broadening:

      ! flag for doppler broadening for each species
      LOGICAL :: L_DOPPLER_PRESENT_SW(NPD_SPECIES)

      ! offset to pressure to represent doppler broadening

      REAL :: DOPPLER_CORRECTION_SW(NPD_SPECIES)

! fields for aerosol optical depth

      ! number of wavelengths
      INTEGER :: N_AOD_WAVEL_SW

      ! monochromatic specific coefficient of absorption
      REAL                                                              &
     &     AOD_ABSORPTION_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &                       NPD_AOD_WAVEL)

      ! monochromatic specific coefficient of scattering
      REAL                                                              &
     &     AOD_SCATTERING_SW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &                       NPD_AOD_WAVEL)

      ! relationship between aerosol component and type
      INTEGER                                                           &
     &     I_AOD_TYPE_SW(NPD_AEROSOL_SPECIES)

! SWSPDL3A end
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE CONTAINING DECLARATIONS FOR REDUCED LW-SPECTRAL FILE.
!     NOTE: LWSPDC3A, LWSPCM3A AND LWSARG3A MUST BE CONSISTENT.
!     NOTE: SINCE THE ARRAYS HERE WILL BE PASSED IN A COMMON BLOCK
!     THEIR SIZES MUST BE FIXED, EVEN THOUGH VARIABLE SIZES ARE USED
!     LOWER IN THE CODE. THEY ARE ACCORDINGLY DEFINED AS 1-DIMENSIONAL
!     ARRAYS WITH FIXED MAXIMUM SIZES AT THIS LEVEL.
!
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
     &     L_PRESENT_LW(0: NPD_TYPE)
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
     &     WAVE_LENGTH_SHORT_LW(NPD_BAND)                               &
!             SHORTER WAVELENGTH LIMITS
     &   , WAVE_LENGTH_LONG_LW(NPD_BAND)
!             LONGER WAVELENGTH LIMITS
!
!
!
!     EXCLUSION OF SPECIFIC BANDS FROM PARTS OF THE SPECTRUM:
!
      INTEGER                                                           &
     &     N_BAND_EXCLUDE_LW(NPD_BAND)                                  &
!             NUMBER OF EXCLUDED BANDS WITHIN EACH SPECTRAL BAND
     &   , INDEX_EXCLUDE_LW(NPD_EXCLUDE, NPD_BAND)
!             INDICES OF EXCLUDED BANDS
!
!
!
!     FIELDS FOR THE SOLAR FLUX:
!
      REAL                                                              &
     &     SOLAR_FLUX_BAND_LW(NPD_BAND)
!             FRACTION OF THE INCIDENT SOLAR FLUX IN EACH BAND
!
!
!
!     FIELDS FOR RAYLEIGH SCATTERING:
!
      REAL                                                              &
     &     RAYLEIGH_COEFFICIENT_LW(NPD_BAND)
!             RAYLEIGH COEFFICIENTS
!
!
!
!     FIELDS FOR GASEOUS ABSORPTION:
!
      INTEGER                                                           &
     &     N_ABSORB_LW                                                  &
!             NUMBER OF ABSORBERS
     &   , N_BAND_ABSORB_LW(NPD_BAND)                                   &
!             NUMBER OF ABSORBERS IN EACH BAND
     &   , INDEX_ABSORB_LW(NPD_SPECIES, NPD_BAND)                       &
!             LIST OF ABSORBERS IN EACH BAND
     &   , TYPE_ABSORB_LW(NPD_SPECIES)                                  &
!             TYPES OF EACH GAS IN THE SPECTRAL FILE
     &   , I_BAND_ESFT_LW(NPD_BAND, NPD_SPECIES)                        &
!             NUMBER OF ESFT TERMS IN EACH BAND FOR EACH GAS
     &   , I_SCALE_ESFT_LW(NPD_BAND, NPD_SPECIES)                       &
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC_LW(NPD_BAND, NPD_SPECIES)
!             TYPE OF SCALING FUNCTION
!
      REAL                                                              &
     &     K_ESFT_LW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)              &
!             ESFT EXPONENTS
     &   , W_ESFT_LW(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)              &
!             ESFT WEIGHTS
     &   , SCALE_VECTOR_LW(NPD_SCALE_VARIABLE, NPD_ESFT_TERM            &
     &        , NPD_BAND, NPD_SPECIES)                                  &
!             SCALING PARAMETERS FOR EACH ABSORBER AND TERM
     &   , P_REFERENCE_LW(NPD_SPECIES, NPD_BAND)                        &
!             REFERENCE PRESSURE FOR SCALING FUNCTION
     &   , T_REFERENCE_LW(NPD_SPECIES, NPD_BAND)
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
     &     THERMAL_COEFFICIENT_LW(0: NPD_THERMAL_COEFF-1                &
     &        , NPD_BAND)                                               &
!             COEFFICIENTS IN POLYNOMIAL FIT TO SOURCE FUNCTION
     &   , T_REF_PLANCK_LW
!             PLANCKIAN REFERENCE TEMPERATURE
!
!
!
!     SURFACE PROPERTIES:
!
      INTEGER                                                           &
     &     I_SPEC_SURFACE_LW(NPD_SURFACE)                               &
!             METHOD OF SPECIFYING PROPERTIES OF SURFACE
     &   , N_DIR_ALBEDO_FIT_LW(NPD_SURFACE)
!             NUMBER OF PARAMETERS FITTING THE DIRECT ALBEDO
!
      LOGICAL                                                           &
     &     L_SURFACE_LW(NPD_SURFACE)
!             SURFACE TYPES INCLUDED
!
      REAL                                                              &
     &     SURFACE_ALBEDO_LW(NPD_BAND, NPD_SURFACE)                     &
!             SURFACE ALBEDOS
     &   , DIRECT_ALBEDO_PARM_LW(0: NPD_ALBEDO_PARM                     &
     &      , NPD_BAND, NPD_SURFACE)                                    &
!             COEFFICIENTS FOR FITTING DIRECT ALBEDO
     &   , EMISSIVITY_GROUND_LW(NPD_BAND, NPD_SURFACE)
!             SURFACE EMISSIVITIES
!
!
!
!     FIELDS FOR CONTINUA:
!
      INTEGER                                                           &
     &     N_BAND_CONTINUUM_LW(NPD_BAND)                                &
!             NUMBER OF CONTINUA IN EACH BAND
     &   , INDEX_CONTINUUM_LW(NPD_BAND, NPD_CONTINUUM)                  &
!             LIST OF CONTINUA IN EACH BAND
     &   , INDEX_WATER_LW                                               &
!             INDEX OF WATER VAPOUR
     &   , I_SCALE_FNC_CONT_LW(NPD_BAND, NPD_CONTINUUM)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
!
      REAL                                                              &
     &     K_CONTINUUM_LW(NPD_BAND, NPD_CONTINUUM)                      &
!             GREY EXTINCTION COEFFICIENTS FOR CONTINUUM
     &   , SCALE_CONTINUUM_LW(NPD_SCALE_VARIABLE                        &
     &      , NPD_BAND, NPD_CONTINUUM)                                  &
!             SCALING PARAMETERS FOR CONTINUUM
     &   , P_REF_CONTINUUM_LW(NPD_CONTINUUM, NPD_BAND)                  &
!             REFERENCE PRESSURE FOR SCALING OF CONTINUUM
     &   , T_REF_CONTINUUM_LW(NPD_CONTINUUM, NPD_BAND)
!             REFERENCE TEMPERATURE FOR SCALING OF CONTINUUM
!
!
!
!     FIELDS FOR WATER DROPLETS:
!
      INTEGER                                                           &
     &     I_DROP_PARAMETRIZATION_LW(NPD_DROP_TYPE)
!             PARAMETRIZATION TYPE OF DROPLETS
!
      LOGICAL                                                           &
     &     L_DROP_TYPE_LW(NPD_DROP_TYPE)
!             TYPES OF DROPLET PRESENT
!
      REAL                                                              &
     &     DROP_PARAMETER_LIST_LW(NPD_CLOUD_PARAMETER                   &
     &        , NPD_BAND, NPD_DROP_TYPE)                                &
!             PARAMETERS USED TO FIT OPTICAL PROPERTIES OF CLOUDS
     &   , DROP_PARM_MIN_DIM_LW(NPD_DROP_TYPE)                          &
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , DROP_PARM_MAX_DIM_LW(NPD_DROP_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR AEROSOLS:
!
      INTEGER                                                           &
     &     N_AEROSOL_LW                                                 &
!             NUMBER OF SPECIES OF AEROSOL
     &   , TYPE_AEROSOL_LW(NPD_AEROSOL_SPECIES)                         &
!             TYPES OF AEROSOLS
     &   , I_AEROSOL_PARAMETRIZATION_LW(NPD_AEROSOL_SPECIES)            &
!             PARAMETRIZATION OF AEROSOLS
     &   , NHUMIDITY_LW(NPD_AEROSOL_SPECIES)
!             NUMBERS OF HUMIDITIES
!
      LOGICAL                                                           &
     &     L_AEROSOL_SPECIES_LW(NPD_AEROSOL_SPECIES)
!             AEROSOL SPECIES INCLUDED
!
      REAL                                                              &
     &     AEROSOL_ABSORPTION_LW(NPD_HUMIDITIES                         &
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)                          &
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING_LW(NPD_HUMIDITIES                         &
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)                          &
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY_LW(NPD_HUMIDITIES                          &
     &        , NPD_AEROSOL_SPECIES, NPD_BAND)                          &
!             ASYMMETRY OF AEROSOLS
     &   , HUMIDITIES_LW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             HUMIDITIES FOR COMPONENTS
!
!
!
!     FIELDS FOR ICE CRYSTALS:
!
      INTEGER                                                           &
     &     I_ICE_PARAMETRIZATION_LW(NPD_ICE_TYPE)
!             TYPES OF PARAMETRIZATION OF ICE CRYSTALS
!
      LOGICAL                                                           &
     &     L_ICE_TYPE_LW(NPD_ICE_TYPE)
!             TYPES OF ICE CRYSTAL PRESENT
!
      REAL                                                              &
     &     ICE_PARAMETER_LIST_LW(NPD_CLOUD_PARAMETER                    &
     &        , NPD_BAND, NPD_ICE_TYPE)                                 &
!             PARAMETERS USED TO FIT SINGLE SCATTERING OF ICE CRYSTALS
     &   , ICE_PARM_MIN_DIM_LW(NPD_ICE_TYPE)                            &
!             MINIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
     &   , ICE_PARM_MAX_DIM_LW(NPD_ICE_TYPE)
!             MAXIMUM DIMENSION PERMISSIBLE IN THE PARAMETRIZATION
!
!
!
!     FIELDS FOR DOPPLER BROADENING:
!
      LOGICAL                                                           &
     &     L_DOPPLER_PRESENT_LW(NPD_SPECIES)
!             FLAG FOR DOPPLER BROADENING FOR EACH SPECIES
!
      REAL                                                              &
     &     DOPPLER_CORRECTION_LW(NPD_SPECIES)
!             OFFSET TO PRESSURE TO REPRESENT DOPPLER BROADENING
!
!
!
!     FIELDS FOR AEROSOL OPTICAL DEPTH
!
      INTEGER                                                           &
     &     N_AOD_WAVEL_LW
!             NUMBER OF WAVELENGTHS
!
      REAL                                                              &
     &     AOD_ABSORPTION_LW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &                       NPD_AOD_WAVEL)
!             MONOCHROMATIC SPECIFIC COEFFICIENT OF ABSORPTION
!
      REAL                                                              &
     &     AOD_SCATTERING_LW(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &                       NPD_AOD_WAVEL)
!             MONOCHROMATIC SPECIFIC COEFFICIENT OF SCATTERING
!
      INTEGER                                                           &
     &     I_AOD_TYPE_LW(NPD_AEROSOL_SPECIES)
!             RELATIONSHIP BETWEEN AEROSOL COMPONENT AND TYPE
!
!
!
!
!    ------------------------------------------------------------------
!
!     Common blocks of spectral data for 3A-radiation.
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
!     Declaration of algorithmic options for 3A-radiation.
! SWOPT3A defines SW options for the radiation code.
!
!     NOTE: SWOPT3A and SWCAVR3A must be consistent.
!
!     Current Owner of Code: J. M. Edwards
!
!     History:
!     Version  Date      Comment.
!     5.2      05/12/00  Logical for extra top level added.
!                                J. M. Edwards
! 6.2  24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!
!
      INTEGER :: I_2STREAM_SW     ! two_stream scheme selected
      INTEGER :: I_GAS_OVERLAP_SW ! treatment of gaseous overlaps
      INTEGER :: I_CLOUD_SW       ! treatment of cloud overlaps
      INTEGER :: I_CLOUD_REPRESENTATION_SW  ! representation of clouds
      INTEGER :: I_SOLVER_SW                ! solver selected

      LOGICAL :: L_O2_SW ! flag for oxygen

      ! types of droplets or ice crystals used for parametrizations
      INTEGER :: I_ST_WATER_SW  ! type for stratiform water
      INTEGER :: I_CNV_WATER_SW ! type for convective water
      INTEGER :: I_ST_ICE_SW    ! type for stratiform ice
      INTEGER :: I_CNV_ICE_SW   ! type for convective ice

      ! flag to partition convective clouds using the local temperature
      LOGICAL :: L_LOCAL_CNV_PARTITION_SW


      ! Flag to include an extra layer at the top of the atmosphere in
      ! radiative calculations
      LOGICAL :: L_EXTRA_TOP_SW

! SWOPT3A end
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
!     Common blocks of algorithmic options for 3A-radiation.
!     ------------------------------------------------------------------
!     COMMON BLOCK CONTAINING OPTIONS FOR 3A-RADIATION CODE.
!     VERSION FOR SHORTWAVE CALCULATIONS.
!
      COMMON/R2SWOPT/                                                   &
     &  I_2STREAM_SW, I_GAS_OVERLAP_SW, I_CLOUD_SW,                     &
     &  I_CLOUD_REPRESENTATION_SW, I_SOLVER_SW,                         &
     &  L_O2_SW,                                                        &
     &  I_ST_WATER_SW, I_CNV_WATER_SW, I_ST_ICE_SW, I_CNV_ICE_SW,       &
     &  L_LOCAL_CNV_PARTITION_SW, L_EXTRA_TOP_SW
!
!     ------------------------------------------------------------------
! LWCOPT3A contains options for 3a-radiation code version for
! longwave calculations.

      COMMON/R2LWOPT/                                                   &
     &  I_2STREAM_LW, L_IR_SOURCE_QUAD_LW, I_GAS_OVERLAP_LW,            &
     &  I_CLOUD_LW, I_CLOUD_REPRESENTATION_LW, I_SOLVER_LW,             &
     &  I_SCATTER_METHOD_LW,                                            &
     &  L_N2O_LW, L_CH4_LW, L_CFC11_LW , L_CFC12_LW, L_CFC113_LW,       &
     &  L_HCFC22_LW, L_HFC125_LW, L_HFC134A_LW,                         &
     &  I_ST_WATER_LW, I_CNV_WATER_LW, I_ST_ICE_LW,                     &
     &  I_CNV_ICE_LW, L_MICROPHYSICS_LW,                                &
     &  L_LOCAL_CNV_PARTITION_LW, L_EXTRA_TOP_LW
! LWCOPT3A end

!     Component IDs for the aerosol climatology for NWP
! arcl_ids.h
!
! Constants used by the aerosol climatology for NWP
!
      !
      ! Available species:
      !
      !   1. Sulphate
      !   2. Mineral dust
      !   3. Sea-salt
      !   4. Black-carbon (also named soot)
      !   5. Biomass-burning
      !   6. Fossil-fuel organic carbon
      !   7. Delta aerosol
      !
      ! Note: When adding a species, increase NPD_ARCL_SPECIES in
      !       arcl_dim.h. NPD_ARCL_SPECIES must be equal to the
      !       largest IP_ARCL_????

      integer, parameter :: IP_ARCL_SULP = 1
      integer, parameter :: IP_ARCL_DUST = 2
      integer, parameter :: IP_ARCL_SSLT = 3
      integer, parameter :: IP_ARCL_BLCK = 4
      integer, parameter :: IP_ARCL_BIOM = 5
      integer, parameter :: IP_ARCL_OCFF = 6
      integer, parameter :: IP_ARCL_DLTA = 7
      
      !
      ! List of components.
      !
      !   The number of components depends on the species:
      !
      !   1. Sulphate: 3 (accumulation, Aitken, dissolved)
      !   2. Mineral dust: 6 size bins
      !   3. Sea-salt: 2 (film and jet)
      !   4. Black-carbon: 2 (fresh and aged)
      !   5. Biomass-burning: 3 (fresh, aged, in-cloud)
      !   6. Fossil-fuel organic carbon: 3 (fresh, aged, in-cloud)
      !   7. Delta aerosol: 1
      !
      ! Note: when adding a component, increase NPD_ARCL_COMPNTS
      !       in arcl_dim.h. NPD_ARCL_COMPNTS must be equal to the
      !       largest IP_ARCL_????_??.
      !       Array i_arcl_compnts_per_species in set_arcl_dimensions()
      !       will also have to be updated.
      !
      
      integer, parameter :: IP_ARCL_SULP_AC = 1
      integer, parameter :: IP_ARCL_SULP_AK = 2
      integer, parameter :: IP_ARCL_SULP_DI = 3
      integer, parameter :: IP_ARCL_DUST_B1 = 4
      integer, parameter :: IP_ARCL_DUST_B2 = 5
      integer, parameter :: IP_ARCL_DUST_B3 = 6
      integer, parameter :: IP_ARCL_DUST_B4 = 7
      integer, parameter :: IP_ARCL_DUST_B5 = 8
      integer, parameter :: IP_ARCL_DUST_B6 = 9
      integer, parameter :: IP_ARCL_SSLT_FI = 10
      integer, parameter :: IP_ARCL_SSLT_JT = 11
      integer, parameter :: IP_ARCL_BLCK_FR = 12
      integer, parameter :: IP_ARCL_BLCK_AG = 13
      integer, parameter :: IP_ARCL_BIOM_FR = 14
      integer, parameter :: IP_ARCL_BIOM_AG = 15
      integer, parameter :: IP_ARCL_BIOM_IC = 16
      integer, parameter :: IP_ARCL_OCFF_FR = 17
      integer, parameter :: IP_ARCL_OCFF_AG = 18
      integer, parameter :: IP_ARCL_OCFF_IC = 19
      integer, parameter :: IP_ARCL_DLTA_DL = 20

! end of arcl_ids.h

            Real :: tmp_abso

! External Routines:
      External Swap_Bounds
      External r2_global_cloud_top, solpos, ftsa, solang
      External solinc
      External r2_swrad, tropin
      External gc_imax
      External r2_lwrad, timer
      External flux_diag
      External tile_albedo
      External tilepts
      External set_seasalt
      External isccp
      External diagnostics_lw
      External diagnostics_sw
      External rad3d_inp

      ! prevent tripping -check uninit when l_rad_deg==.false.
      tot_daylight_points=0

! ----------------------------------------------------------------------
! Section INI. Initialisation of stash output variables.
! ----------------------------------------------------------------------
!kdcorbin, 06/10 - added switch for radiation interaction
     L_CO2_3D = L_CO2_radiation


!     Map the stashflags to the internal controlling logicals.
      L_T_incr_sw                 = (sf(181,1) .or. sf(161,1))
      SW_diag%L_solar_out_toa     = sf(208,1)
      SW_diag%L_solar_out_clear   = sf(209,1)
      SW_diag%L_surface_down_flux = (sf(235,1) .OR. L_MOSES_II .OR. L_CABLE)
      SW_diag%L_surf_down_clr     = sf(210,1)
      SW_diag%L_surf_up_clr       = sf(211,1)
      SW_diag%L_flux_direct       = sf(230,1)
      SW_diag%L_flux_diffuse      = sf(231,1)     
      SW_diag%L_clear_hr          = sf(233,1)
      SW_diag%L_net_flux_trop     = sf(237,1)
      SW_diag%L_up_flux_trop      = sf(238,1)
! Microphysical diagnostics
      SW_diag%re_strat_flag       = sf(221,1)
      SW_diag%wgt_strat_flag      = sf(223,1)
      SW_diag%lwp_strat_flag      = sf(224,1)
      SW_diag%re_conv_flag        = sf(225,1)
      SW_diag%wgt_conv_flag       = sf(226,1)
      SW_diag%ntot_diag_flag      = sf(241,1)
      SW_diag%strat_lwc_diag_flag = sf(242,1)
      SW_diag%so4_ccn_diag_flag   = sf(243,1)
      SW_diag%cond_samp_wgt_flag  = sf(244,1)
      SW_diag%weighted_re_flag    = sf(245,1)
      SW_diag%sum_weight_re_flag  = sf(246,1)
      SW_diag%wgtd_warm_re_flag   = sf(254,1)
      SW_diag%sum_wgt_warm_re_flag= sf(255,1)
      SW_diag%seasalt_film_flag   = sf(247,1)
      SW_diag%seasalt_jet_flag    = sf(248,1)
      SW_diag%Nc_diag_flag        = sf(280,1)
      SW_diag%Nc_weight_flag      = sf(281,1)
! PAR diagnostics
      SW_diag%L_direct_par        = sf(291,1)
! Diagnostics for MOSES II
      SW_diag%L_FlxSolBelow690nmSurf = (sf(259,1) .OR.                  &
     &  (L_flux_below_690nm_surf .AND. L_ctile) )
      SW_diag%L_FlxSeaBelow690nmSurf = (sf(260,1) .OR.                  &
     &  (L_flux_below_690nm_surf .AND. L_ctile) )
! Diagnostics for orography correction
      If (.not.L_orog) Then
         sf(292:296,1) = .False.
      End If
      SW_diag%L_sol_bearing = sf(292,1)
      SW_diag%L_orog_corr   = sf(295,1)
!  Isccp diagnostics
      LW_diag%L_isccp_weights                  = sf(269,2)
      LW_diag%L_isccp_cf                       = sf(270,2)
      LW_diag%L_isccp_cf_tau_0_to_p3           = sf(271,2)
      LW_diag%L_isccp_cf_tau_p3_to_1p3         = sf(272,2)
      LW_diag%L_isccp_cf_tau_1p3_to_3p6        = sf(273,2)
      LW_diag%L_isccp_cf_tau_3p6_to_9p4        = sf(274,2)
      LW_diag%L_isccp_cf_tau_9p4_to_23         = sf(275,2)
      LW_diag%L_isccp_cf_tau_23_to_60          = sf(276,2)
      LW_diag%L_isccp_cf_tau_ge_60             = sf(277,2)
      LW_diag%L_meanalbedocld                  = sf(290,2)
      LW_diag%L_meantaucld                     = sf(291,2)
      LW_diag%L_meanptop                       = sf(292,2)
      LW_diag%L_totalcldarea                   = sf(293,2)
! Extinction and absorptivity diagnostics
      SW_diag%L_cloud_extinction               = sf(262,1)
      SW_diag%L_cloud_weight_extinction        = sf(263,1)
      SW_diag%L_ls_cloud_extinction            = (sf(264,1) .OR.        &
     &  (LW_diag%L_isccp_weights))
      SW_diag%L_ls_cloud_weight_extinction     = (sf(265,1) .OR.        &
     &  (LW_diag%L_isccp_weights))
      SW_diag%L_cnv_cloud_extinction           = (sf(266,1) .OR.        &
     &  (LW_diag%L_isccp_weights))
      SW_diag%L_cnv_cloud_weight_extinction    = (sf(267,1) .OR.        &
     &  (LW_diag%L_isccp_weights))
      LW_diag%L_cloud_absorptivity             = sf(262,2)
      LW_diag%L_cloud_weight_absorptivity      = sf(263,2)
      LW_diag%L_ls_cloud_absorptivity          = (sf(264,2) .OR.        &
     &  (LW_diag%L_isccp_weights))
      LW_diag%L_ls_cloud_weight_absorptivity   = (sf(265,2) .OR.        &
     &  (LW_diag%L_isccp_weights))
      LW_diag%L_cnv_cloud_absorptivity         = (sf(266,2) .OR.        &
     &  (LW_diag%L_isccp_weights))
      LW_diag%L_cnv_cloud_weight_absorptivity  = (sf(267,2) .OR.        &
     &  (LW_diag%L_isccp_weights))
      LW_diag%L_total_cloud_on_levels          = sf(261,2)
! Grid-box mean cloud diagnostics as seen by radiation:
      LW_diag%L_ls_qcl_rad                     = sf(308,2)
      LW_diag%L_ls_qcf_rad                     = sf(309,2)
      LW_diag%L_cc_qcl_rad                     = sf(310,2)
      LW_diag%L_cc_qcf_rad                     = sf(311,2)
      LW_diag%L_ls_cl_rad                      = sf(312,2)
      LW_diag%L_ls_cf_rad                      = sf(313,2)
      LW_diag%L_cc_cl_rad                      = sf(314,2)
      LW_diag%L_cc_cf_rad                      = sf(315,2)

! Switches enabled as STASHflags: LW
      L_T_incr_lw                    = (sf(181,2) .or. sf(161,2))
      LW_diag%L_total_cloud_cover    = sf(204,2)
      LW_diag%L_clear_olr            = sf(206,2)
      LW_diag%L_surface_down_flux    = (L_MOSES_II .OR. L_CABLE .OR. sf(207,2) )
      LW_diag%L_surf_down_clr        = sf(208,2)
      LW_diag%L_clear_hr             = sf(233,2)
      LW_diag%L_net_flux_trop        = sf(237,2)
      LW_diag%L_down_flux_trop       = sf(238,2)

! Aerosol optical depth diagnostics
      LW_diag%L_aod_sulphate         = sf(284, 2)
      LW_diag%L_aod_dust             = sf(285, 2)
      LW_diag%L_aod_seasalt          = sf(286, 2)
      LW_diag%L_aod_soot             = sf(287, 2)
      LW_diag%L_aod_biomass          = sf(288, 2)
      LW_diag%L_aod_biogenic         = sf(289, 2)
      LW_diag%L_aod_ocff             = sf(295, 2)
      LW_diag%L_aod_delta            = sf(296, 2)

! ----------------------------------------------------------------------
! Set tile_pts, tile_index if MOSES II selected.
! Set land_index.
! ----------------------------------------------------------------------
      if (L_MOSES_II .OR. L_CABLE)                                      &
! DEPENDS ON: tilepts
     &  call tilepts(land_field,frac,tile_pts,tile_index)
        do l = 1, land_field
          j = (land_index(l)-1)/row_length + 1
          land_index_i(l) = land_index(l) - (j-1)*row_length
          land_index_j(l) = j
        enddo

! ----------------------------------------------------------------------
! Set global land fraction
! ----------------------------------------------------------------------
      do j = 1, rows
        do i = 1, row_length
          flandg(i,j)=0.0
        enddo
      enddo
      do l=1,land_field
        i = land_index_i(l)
        j = land_index_j(l)
        flandg(i,j)=fland(l)
      enddo

! ----------------------------------------------------------------------
! Check that Coastal tiling switch is off if not running with MOSES2:
! ----------------------------------------------------------------------
      IF(.NOT.(L_MOSES_II.OR.L_CABLE).AND.L_CTILE)THEN
        write(6,*)'****ERROR***Coastal tiling cannot be '
        write(6,*)'run with pre-MOSES2 sfc/boundary layer code'
! DEPENDS ON: ereport
        CALL EREPORT('NI_RAD_CTL ', 1000,                               &
     &  'Coastal tiling cannot be run with pre-MOSES2 sfc/boundary layer&
     & code')

      ENDIF
! set CO2 array dimensions for passing to SW & LW schemes
! If used:
! CO2_DIM1 = row_length*rows (calc'd from CO2_DIM_LEN & CO2_DIM_ROW)
! CO2_DIM2 = model_levels (already passed in from CO2_DIM_LEV)
      CO2_DIM1 = CO2_DIM_LEN * CO2_DIM_ROW

! ----------------------------------------------------------------------
! Section RAD.0 Radiation scheme.
!               Assumes no sulphur cycle.
!               This incorporates radiation code for non-radiation
!               timesteps which is in CLDCTL1.dk in the UM.
!               Set up required fields and constants.
!-----------------------------------------------------------------------

      If (error_code  ==  0) Then
! DEPENDS ON: timer
        If (Ltimer) Call timer ('SW Rad  ',3)

! Calculate number of seconds since midnight to the beginning of
! the timetsep.

      seconds_since_midnight = Real( PREVIOUS_TIME(4) * 3600            &
     &         + PREVIOUS_TIME(5) * 60  + PREVIOUS_TIME(6))

! Calculate true latitude at each grid-point.
        Do j = 1, rows
          Do i = 1, row_length
            sin_true_latitude(i,j) = f3_at_u(i,j) / two_Omega
          End Do
        End Do

! Set radiation fields to zero
        where (ice_fract > 0. ) surf_radflux = 0. ! EAK
        Do j = 1, rows
          Do i = 1, row_length
            photosynth_act_rad(i,j) = 0.
!            surf_radflux(i,j) = 0.
            flux_below_690nm_surf(i,j)=0.
            sw_net_land(i,j) = 0.0
            sw_net_sice(i,j) = 0.0
            dOLR(i,j) = 0.
            !Alison McLaren advised bug fix here: commented the folloing line
            !which wrong SW radiation over coastal sea-ice cells, namely, it
            !applys SW Only at the beginning of a radiation step on to these cells
            ! fixed: 20111201.
            !!! land_alb(i,j)=0.0
            !Note this change is also made to the same place in another file:
            !     src/atmosphere/short_wave_radiation/glue_rad-rad_ctl3c.F90
            !     BUT "glue_rad" subroutine HERE is the one that is ACTUALLY USED
            !     in the access hg3-C configuration! (2011-12-21).
             FRACSOLID(I,J)=FLANDG(I,J)                                  &
     &        +(1.-FLANDG(I,J))*ICE_FRACT(I,J)

          End Do
        End Do

        If (L_orog) Then
           Allocate(f_orog(row_length,rows))
           f_orog=0.0
        Else
           Allocate(f_orog(1,1))
           Allocate(slope_aspect(1,1))
           Allocate(slope_angle(1,1))
        End If

        Do n=1,ntiles
          Do l = 1, land_field
            SW_tile(l,n) = 0.
            SW_tile_down(l,n) = 0.
          End Do
        End Do
!
! Code for the mixing ratio calculation of air density.
! Copy coordinate information into arrays that do not contain
! halo information
      If (l_mixing_ratio) then

        ! Level 0 theta level information
        Do j=1,rows
          Do i=1,row_length
            ij=i+(j-1)*row_length
            r_theta_levels_nh(ij,0) = r_theta_levels(i,j,0)
          End do
        End do

        Do k=1,model_levels
          ! Other level information
          Do j=1,rows
            Do i=1,row_length
              ij=i+(j-1)*row_length
              rho_r2_nh(ij,k)   = rho_r2(i,j,k)
              r_theta_levels_nh(ij,k) = r_theta_levels(i,j,k)
              r_rho_levels_nh(ij,k)   = r_rho_levels(i,j,k)
            End do
          End do
        End do

      End If  ! l_mixing_ratio
!
!     Code for the mineral dust aerosol scheme. We copy
!     dust aerosol into local arrays if the direct effect is
!     switched on.
!
      IF (L_USE_DUST) THEN
        IF (DUST_DIM1  ==  ROWS*ROW_LENGTH .AND.                        &
     &                     DUST_DIM2  ==  MODEL_LEVELS) THEN
          DO K=1,MODEL_LEVELS
            DO J=1,ROWS
              DO I=1,ROW_LENGTH
                IJ=I+(J-1)*ROW_LENGTH
                DUST_1(IJ,K) = DUST_DIV1(I,J,K)
                DUST_2(IJ,K) = DUST_DIV2(I,J,K)
                DUST_3(IJ,K) = DUST_DIV3(I,J,K)
                DUST_4(IJ,K) = DUST_DIV4(I,J,K)
                DUST_5(IJ,K) = DUST_DIV5(I,J,K)
                DUST_6(IJ,K) = DUST_DIV6(I,J,K)
             END DO
            END DO
          END DO
        ELSE
          WRITE(6,*)                                                    &
     &       'RAD_CTL: DUST_DIM INCONSISTENT WITH L_USE_DUST'
          ERROR_CODE=1
          RETURN
        END IF
      END IF
!
!     Code for the biogenic aerosol. Copy into local arrays if
!     this aerosol was requested.
!
      IF (L_USE_BIOGENIC) THEN
        IF (BIOGENIC_DIM1  ==  rows*row_length .AND.                    &
     &      BIOGENIC_DIM2  ==  model_levels) THEN
          DO k=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                IJ=i+(j-1)*row_length
                local_biogenic(IJ,k) = biogenic(i,j,k)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          WRITE(6,*)                                                    &
     &        'BIOGENIC_DIM INCONSISTENT WITH L_USE_BIOGENIC'
          ERROR_CODE=1
          RETURN
        END IF
      END IF
!
!     Code for the Sulphur Cycle. We multiply by 4.125 to convert from
!     mass mixing ratio of sulphur atoms to mass mixing ratio of
!     ammonium sulphate.
      IF (L_USE_SULPC_DIRECT .OR. L_USE_SULPC_INDIRECT_SW               &
     &                       .OR. L_USE_SULPC_INDIRECT_LW) THEN
        IF (SULP_DIM1  ==  rows*row_length .AND.                        &
     &                       SULP_DIM2  ==  model_levels)  THEN
          Do k=1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                IJ=i+(j-1)*row_length
                accum_sulphate(IJ,k)=so4_accu(i,j,k)*4.125
                aitken_sulphate(IJ,k)=so4_aitken(i,j,k)*4.125
                diss_sulphate(IJ,k)=so4_diss(i,j,k)*4.125
              End Do
            End Do
          End Do
        ELSE
         WRITE(6,*)                                                     &
     &       'SULP_DIM INCONSISTENT WITH L_USE_SULPC, EXIT RAD_CTL'
         Error_Code = 1
         RETURN
        ENDIF
      ENDIF
!
!     Code for the Soot Scheme. As for the Sulphur Cycle (above), we
!     copy soot into local arrays if the direct effect is switched on,
!     but no multiplication is required.
!
      If (L_use_soot_direct) Then
        If (soot_dim1  ==  rows*row_length .and.                        &
     &                     soot_dim2  ==  model_levels) Then
          Do k=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                ij=i+(j-1)*row_length
                fresh_soot(ij,k) = soot_new(i,j,k)
                aged_soot(ij,k) = soot_agd(i,j,k)
              End Do
            End Do
          End Do
        Else
          Write(6,*)                                                    &
     &       'RAD_CTL: SOOT_DIM INCONSISTENT WITH L_USE_SOOT_DIRECT'
          Error_Code=1
          Return
        End If
      End If
!
!     Code for the biomass aerosol scheme. As for soot, we copy biomass
!     smoke aerosol into local arrays if the direct/indirect effect is
!     switched on.
!
      If (L_use_bmass_direct .or. L_use_bmass_indirect) Then
        If (bmass_dim1  ==  rows*row_length .and.                       &
     &                     bmass_dim2  ==  model_levels) Then
          Do k=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                ij=i+(j-1)*row_length
                fresh_bmass(ij,k) = bmass_new(i,j,k)
                aged_bmass(ij,k)  = bmass_agd(i,j,k)
                cloud_bmass(ij,k) = bmass_cld(i,j,k)
              End Do
            End Do
          End Do
        Else
          Write(6,*)                                                    &
     &       'RAD_CTL: BMASS_DIM INCONSISTENT WITH L_USE_BMASS'
          Error_Code=1
          Return
        End If
      End If

!
!     Code for the fossil-fuel organic carbon aerosol scheme. As for 
!     biomass, we copy fossil-fuel organic carbon aerosol into local 
!     arrays if the direct/indirect effect is switched on.
!
      If (L_use_ocff_direct .or. L_use_ocff_indirect) Then
        If (ocff_dim1  ==  rows*row_length .and.                       &
     &      ocff_dim2  ==  model_levels) Then
          Do k=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                ij=i+(j-1)*row_length
                fresh_ocff(ij,k) = ocff_new(i,j,k)
                aged_ocff(ij,k)  = ocff_agd(i,j,k)
                cloud_ocff(ij,k) = ocff_cld(i,j,k)
              End Do
            End Do
          End Do
        Else
          Write(6,*)                                                    &
     &       'RAD_CTL: OCFF_DIM INCONSISTENT WITH L_USE_OCFF'
          Error_Code=1
          Return
        End If
      End If
!
!
!     Code for the NWP aerosol climatology. As above, we copy the input
!     into local arrays if the climatology was requested.
!
      If (n_arcl_species > 0) Then
        
        Do l = 1, NPD_ARCL_COMPNTS
        
          If (i_arcl_compnts(l) /= -1) Then
            If (arcl_dim1 == rows * row_length .and.                    &
     &          arcl_dim2 == model_levels) Then
              
              !
              ! Sulphur mass-mixing ratios have to be converted
              ! into ammonium sulphate.
              !
              If ((l == IP_ARCL_SULP_AC).or.(l == IP_ARCL_SULP_AK)) Then
                arcl_multf = 4.125
              Else
                arcl_multf = 1.0
              End If
              
              Do k=1, model_levels
                Do j=1, rows
                  Do i=1, row_length
                    ij=i+(j-1)*row_length
                    local_arcl(ij,k,i_arcl_compnts(l)) =                &
     &                    arcl_multf *                                  &
     &                    arcl(i,j,k,i_arcl_compnts(l))
                  End Do ! i
                End Do ! j
              End Do ! k
            Else
              Write(6,*)                                                &
     &           'RAD_CTL: ARCL_DIM IS INCONSISTENT WITH L_USE_ARCL'
              Error_Code=1
              Return
            End If
          End If
        
        End Do
      End If ! n_arcl_species
!
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                T_latest(i,j,k)= T_n(i,j,k)
              End Do
            End Do
          End Do
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_latest(i,j,k)= q_n(i,j,k)
                qcl_latest(i,j,k)=qcl_n(i,j,k)
                cf_latest(i,j,k)=cf_n(i,j,k)
                cfl_latest(i,j,k)=cfl_n(i,j,k)
                cff_latest(i,j,k)=cff_n(i,j,k)
                zeros(i,j,k)=0.0
              End Do
            End Do
          End Do
!
        If (L_Rad_Step) Then
          height_rho(:,:,:)=r_rho_levels(1:row_length,1:rows,:)
          height_theta(:,:,:)=r_theta_levels(1:row_length,1:rows,:)
        End If
!
! If we use this then more code is required from RADCTL
        If ( L_Rad_step .and. (l_climat_aerosol         .or.            &
     &                         LW_diag%l_net_flux_trop  .or.            &
     &                         LW_diag%l_down_flux_trop .or.            &
     &                         SW_diag%l_net_flux_trop  .or.            &
     &                         SW_diag%l_up_flux_trop )) Then

! need to copy t_n into an array with haloes and swopbound
! Only required for tropin call.

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                T_large(i,j,k) = T_n(i,j,k)
              End Do
            End Do
          End Do

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &         T_large, row_length, rows,                               &
     &         model_levels, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
          Call FILL_EXTERNAL_HALOS(T_Large,row_length,rows,             &
     &                             model_levels,off_x,off_y)

!         The variables min_trop_level and max_trop_level
!         index rho-levels in the control routines,
!         which start from a zeroth level at the surface.
!         Layer boundaries in physical routines are indexed with
!         the convention that the lowest is indexed by 1.
!         Since the first rho-level is omitted from the physics
!         grid, the two methods of indexing refer to the same
!         horizontal level from the second rho-level upwards,
!         so there is actually no need to adjust these variables
!         for the change of indexing convention.
! DEPENDS ON: tropin
          Call tropin (T_large, exner_rho_levels, exner_theta_levels,   &
     &                 row_length, rows, model_levels, off_x, off_y,    &
     &                 at_extremity,asin(sin_true_latitude(1,1)),       &
     &                 height_theta(1,1,0:model_levels),                &
     &                 min_trop_level, max_trop_level, trindx )

        End If

        If (L_Rad_step .AND.                                            &
     &        (L_use_seasalt_indirect .OR. L_use_seasalt_direct)) Then
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                height(i,j,k) = r_theta_levels(i,j,k)                   &
     &                                          -r_theta_levels(i,j,0)
              End Do
            End Do
          End Do

! DEPENDS ON: set_seasalt
          Call set_seasalt(u_1, v_1, height, flandg, ice_fract          &
     &                    , row_length, rows, model_levels              &
     &                    , bl_levels, sea_salt_film, sea_salt_jet)

        Else
          sea_salt_film(1, 1, 1)=0.0
          sea_salt_jet(1, 1, 1)=0.0
        End if
!
! ----------------------------------------------------------------------
! Section RAD.1 Short Wave Radiation Code.
!-----------------------------------------------------------------------

!
!       Allocate space for those SW diagnostic arrays required and zero
!       the elements explicitly; we will overwrite at lit points.
!
        If (SW_diag%L_solar_out_toa) then
          allocate(SW_diag%solar_out_toa(row_length,rows))
          SW_diag%solar_out_toa(:,:) = 0.0
        End if
!
        If (SW_diag%L_solar_out_clear) then
          allocate(SW_diag%solar_out_clear(row_length,rows))
          SW_diag%solar_out_clear(:,:) = 0.0
        End if
!
        If (SW_diag%L_surface_down_flux) then
          allocate(SW_diag%surface_down_flux(row_length,rows))
          SW_diag%surface_down_flux(:,:) = 0.0
        End if
!
        If (SW_diag%L_surf_down_clr) then
          allocate(SW_diag%surf_down_clr(row_length,rows))
          SW_diag%surf_down_clr(:,:) = 0.0
        End if
!
        If (SW_diag%L_surf_up_clr) then
          allocate(SW_diag%surf_up_clr(row_length,rows))
          SW_diag%surf_up_clr(:,:) = 0.0
        End if
!
        If (SW_diag%L_net_flux_trop) then
          allocate(SW_diag%net_flux_trop(row_length, rows))
          SW_diag%net_flux_trop(:,:) = 0.0
        End if
!
        If (SW_diag%L_up_flux_trop) then
          allocate(SW_diag%up_flux_trop(row_length, rows))
          SW_diag%up_flux_trop(:,:) = 0.0
        End if
!
        If (SW_diag%L_clear_hr) then
          allocate(SW_diag%clear_hr(row_length, rows, model_levels))
          SW_diag%clear_hr(:,:,:) = 0.0
        End if
!
! The diffuse and direct downward flux
        
        If (SW_diag%L_flux_direct) then
          allocate(SW_diag%flux_direct(row_length,rows,model_levels+1))
          SW_diag%flux_direct(:,:,:)=0.0
        End if  
        
        If (SW_diag%L_flux_diffuse) then
          allocate(SW_diag%flux_diffuse(row_length,rows,model_levels+1))
          SW_diag%flux_diffuse(:,:,:)=0.0
        End if  
!
! Microphysical diagnostics
!
        If (SW_diag%re_strat_flag) then
          allocate(SW_diag%re_strat(row_length, rows, cloud_levels))
          SW_diag%re_strat(:,:,:) = 0.0
        End if
        If (SW_diag%wgt_strat_flag) then
          allocate(SW_diag%wgt_strat(row_length, rows, cloud_levels))
          SW_diag%wgt_strat(:,:,:) = 0.0
        End if
        If (SW_diag%lwp_strat_flag) then
          allocate(SW_diag%lwp_strat(row_length, rows, cloud_levels))
          SW_diag%lwp_strat(:,:,:) = 0.0
        End if
        If (SW_diag%re_conv_flag) then
          allocate(SW_diag%re_conv(row_length, rows, cloud_levels))
          SW_diag%re_conv(:,:,:) = 0.0
        End if
        If (SW_diag%wgt_conv_flag) then
          allocate(SW_diag%wgt_conv(row_length, rows, cloud_levels))
          SW_diag%wgt_conv(:,:,:) = 0.0
        End if
        If (SW_diag%ntot_diag_flag) then
          allocate(SW_diag%ntot_diag(row_length, rows, cloud_levels))
          SW_diag%ntot_diag(:,:,:) = 0.0
        End if
        If (SW_diag%strat_lwc_diag_flag) then
          allocate(SW_diag%strat_lwc_diag(row_length, rows,             &
     &                                    cloud_levels))
          SW_diag%strat_lwc_diag(:,:,:) = 0.0
        End if
        If (SW_diag%so4_ccn_diag_flag) then
          allocate(SW_diag%so4_ccn_diag(row_length, rows,               &
     &                                  cloud_levels))
          SW_diag%so4_ccn_diag(:,:,:) = 0.0
        End if
!
        If (SW_diag%cond_samp_wgt_flag) then
          allocate(SW_diag%cond_samp_wgt(row_length, rows,              &
     &                                   cloud_levels))
          SW_diag%cond_samp_wgt(:,:,:) = 0.0
        End if
!
        If (SW_diag%weighted_re_flag) then
          allocate(SW_diag%weighted_re(row_length, rows))
          SW_diag%weighted_re(:,:) = 0.0
        End if
!
        If (SW_diag%sum_weight_re_flag) then
          allocate(SW_diag%sum_weight_re(row_length, rows))
          SW_diag%sum_weight_re(:,:) = 0.0
        End if
!
        If (SW_diag%wgtd_warm_re_flag) then
          allocate(SW_diag%weighted_warm_re(row_length, rows))
          SW_diag%weighted_warm_re(:,:) = 0.0
        End if
!
        If (SW_diag%sum_wgt_warm_re_flag) then
          allocate(SW_diag%sum_weight_warm_re(row_length, rows))
          SW_diag%sum_weight_warm_re(:,:) = 0.0
        End if
!
        If (SW_diag%Nc_diag_flag) then
          allocate(SW_diag%Nc_diag(row_length, rows))
          SW_diag%Nc_diag(:,:) = 0.0
        End if
!
        If (SW_diag%Nc_weight_flag) then
          allocate(SW_diag%Nc_weight(row_length, rows))
          SW_diag%Nc_weight(:,:) = 0.0
        End if
!
!       Diagnostics for MOSES II:
!
        If (SW_diag%L_FlxSolBelow690nmSurf) Then
          allocate(SW_diag%FlxSolBelow690nmSurf(row_length, rows))
          SW_diag%FlxSolBelow690nmSurf(:,:) = 0.0
        End if
!
        If (SW_diag%L_FlxSeaBelow690nmSurf) Then
          allocate(SW_diag%FlxSeaBelow690nmSurf(row_length, rows))
          SW_diag%FlxSeaBelow690nmSurf(:,:) = 0.0
        End if
!
!
!       Direct component of PAR flux for STOCHEM
!
        If (SW_diag%L_direct_par) Then
          allocate(SW_diag%flxdirparsurf(row_length, rows))
          SW_diag%flxdirparsurf(:,:) = 0.0
        End if
!
!
!       Orography correction diagnostics:

        If (SW_diag%L_orog_corr) Then
          Allocate(SW_diag%orog_corr(row_length, rows))
          SW_diag%orog_corr(:,:) = 1.0
        End If

        If (SW_diag%L_sol_bearing) Then
          Allocate(SW_diag%sol_bearing(row_length, rows))
          SW_diag%sol_bearing(:,:) = 0.0
        End If
!
!
!        Extinction diagnostics:
!
        If (SW_diag%L_cloud_extinction) then
          allocate(SW_diag%cloud_extinction(row_length, rows,           &
     &                                    cloud_levels))
          SW_diag%cloud_extinction(:,:,:) = 0.0
        End if
!
        If (SW_diag%L_cloud_weight_extinction) then
          allocate(SW_diag%cloud_weight_extinction(row_length, rows,    &
     &                                    cloud_levels))
          SW_diag%cloud_weight_extinction(:,:,:) = 0.0
        End if
!
        If (SW_diag%L_ls_cloud_extinction) then
          allocate(SW_diag%ls_cloud_extinction(row_length, rows,        &
     &                                    cloud_levels))
          SW_diag%ls_cloud_extinction(:,:,:) = 0.0
        End if
!
        If (SW_diag%L_ls_cloud_weight_extinction) then
          allocate(SW_diag%ls_cloud_weight_extinction(row_length, rows, &
     &                                    cloud_levels))
          SW_diag%ls_cloud_weight_extinction(:,:,:) = 0.0
        End if
!
        If (SW_diag%L_cnv_cloud_extinction) then
          allocate(SW_diag%cnv_cloud_extinction(row_length, rows,       &
     &                                    cloud_levels))
          SW_diag%cnv_cloud_extinction(:,:,:) = 0.0
        End if
!
        If (SW_diag%L_cnv_cloud_weight_extinction) then
          allocate(SW_diag%cnv_cloud_weight_extinction(row_length,      &
     &                                    rows, cloud_levels))
          SW_diag%cnv_cloud_weight_extinction(:,:,:) = 0.0
        End if
!
!
!       Repeat for the case when diagnostics are not selected, but
!       allocate only minimal sizes. This should not be required,
!       but seems to be necessary on certain compilers.
        If (.NOT.SW_diag%L_solar_out_toa) Then
          allocate(SW_diag%solar_out_toa(1,1))
        End if
!
        If (.NOT.SW_diag%L_solar_out_clear) Then
          allocate(SW_diag%solar_out_clear(1,1))
        End if
!
        If (.NOT.SW_diag%L_surface_down_flux) Then
          allocate(SW_diag%surface_down_flux(1,1))
        End if
!
        If (.NOT.SW_diag%L_surf_down_clr) Then
          allocate(SW_diag%surf_down_clr(1,1))
        End if
!
        If (.NOT.SW_diag%L_surf_up_clr) Then
          allocate(SW_diag%surf_up_clr(1,1))
        End if
!
        If (.NOT.SW_diag%L_net_flux_trop) Then
          allocate(SW_diag%net_flux_trop(1, 1))
        End if
!
        If (.NOT.SW_diag%L_up_flux_trop) Then
          allocate(SW_diag%up_flux_trop(1, 1))
        End if
!
        If (.NOT.SW_diag%L_clear_hr) Then
          allocate(SW_diag%clear_hr(1, 1, 1))
        End if
!
! Direct and Diffuse Downward Flux
!       
        If (.NOT.SW_diag%L_flux_direct) Then
          allocate(SW_diag%flux_direct(1,1,1))
        End if  
        
        If (.NOT.SW_diag%L_flux_diffuse) Then
          allocate(SW_diag%flux_diffuse(1,1,1))
        End if  
!
! Microphysical diagnostics
!
        If (.NOT.SW_diag%re_strat_flag) Then
          allocate(SW_diag%re_strat(1, 1, 1))
        End if
        If (.NOT.SW_diag%wgt_strat_flag) Then
          allocate(SW_diag%wgt_strat(1, 1, 1))
        End if
        If (.NOT.SW_diag%lwp_strat_flag) Then
          allocate(SW_diag%lwp_strat(1, 1, 1))
        End if
        If (.NOT.SW_diag%re_conv_flag) Then
          allocate(SW_diag%re_conv(1, 1, 1))
        End if
        If (.NOT.SW_diag%wgt_conv_flag) Then
          allocate(SW_diag%wgt_conv(1, 1, 1))
        End if
        If (.NOT.SW_diag%ntot_diag_flag) Then
          allocate(SW_diag%ntot_diag(1, 1, 1))
        End if
        If (.NOT.SW_diag%strat_lwc_diag_flag) Then
          allocate(SW_diag%strat_lwc_diag(1, 1, 1))
        End if
        If (.NOT.SW_diag%so4_ccn_diag_flag) Then
          allocate(SW_diag%so4_ccn_diag(1, 1, 1))
        End if
!
        If (.NOT.SW_diag%cond_samp_wgt_flag) Then
          allocate(SW_diag%cond_samp_wgt(1, 1, 1))
        End if
!
        If (.NOT.SW_diag%weighted_re_flag) Then
          allocate(SW_diag%weighted_re(1, 1))
        End if
!
        If (.NOT.SW_diag%sum_weight_re_flag) Then
          allocate(SW_diag%sum_weight_re(1, 1))
        End if
!
        If (.NOT.SW_diag%wgtd_warm_re_flag) Then
          allocate(SW_diag%weighted_warm_re(1, 1))
        End if
!
        If (.NOT.SW_diag%sum_wgt_warm_re_flag) Then
          allocate(SW_diag%sum_weight_warm_re(1, 1))
        End if
!
        If (.NOT.SW_diag%Nc_diag_flag) Then
          allocate(SW_diag%Nc_diag(1, 1))
        End if
!
        If (.NOT.SW_diag%Nc_weight_flag) Then
          allocate(SW_diag%Nc_weight(1, 1))
        End if
!
!       Diagnostics for MOSES II
!
        If (.NOT.SW_diag%L_FlxSolBelow690nmSurf) Then
          allocate(SW_diag%FlxSolBelow690nmSurf(1, 1))
        End if
!
        If (.NOT.SW_diag%L_FlxSeaBelow690nmSurf) Then
          allocate(SW_diag%FlxSeaBelow690nmSurf(1, 1))
        End if
!
!
!       Direct component of PAR flux for STOCHEM
!
        If (.NOT.SW_diag%L_direct_par) Then
          allocate(SW_diag%flxdirparsurf(1, 1))
        End if
!
!
!       Orography correction diagnostics:

        If (.Not.SW_diag%L_orog_corr) Then
          Allocate(SW_diag%orog_corr(1,1))
        End If

        If (.Not.SW_diag%L_sol_bearing) Then
          Allocate(SW_diag%sol_bearing(1,1))
        End If
!
!
!         Extinction diagnostics:
!
        If (.NOT.SW_diag%L_cloud_extinction) Then
          allocate(SW_diag%cloud_extinction(1, 1, 1))
        End if
        If (.NOT.SW_diag%L_cloud_weight_extinction) Then
          allocate(SW_diag%cloud_weight_extinction(1, 1, 1))
        End if
        If (.NOT.SW_diag%L_ls_cloud_extinction) Then
          allocate(SW_diag%ls_cloud_extinction(1, 1, 1))
        End if
        If (.NOT.SW_diag%L_ls_cloud_weight_extinction) Then
          allocate(SW_diag%ls_cloud_weight_extinction(1, 1, 1))
        End if
        If (.NOT.SW_diag%L_cnv_cloud_extinction) Then
          allocate(SW_diag%cnv_cloud_extinction(1, 1, 1))
        End if
        If (.NOT.SW_diag%L_cnv_cloud_weight_extinction) Then
          allocate(SW_diag%cnv_cloud_weight_extinction(1, 1, 1))
        End if
!
! Calculates sine of the solar declination and the scaling
! factor for solar intensity from the day number and year.

! DEPENDS ON: solpos
        Call solpos (PREVIOUS_TIME(7), PREVIOUS_TIME(1),                &
     &               L_CAL360, L_SEC_VAR, L_EqT, Eq_Time, Sindec, SCS)

!  If L_SCVARY set TRUE then
!  Tweak SCS and SW spectral file to allow for solar variability -
!    the total output & its spectral distribution:

      if (L_SCVARY) then
! DEPENDS ON: solvar
      Call solvar (PREVIOUS_TIME(1), SCS, NPD_BAND_SW,			&
     &     SOLAR_FLUX_BAND_SW,  RAYLEIGH_COEFFICIENT_SW,		&
     &     W_ESFT_SW(1+2*NPD_ESFT_TERM_SW*NPD_BAND_SW,1,1),		&
!     !  2* to skip two species & get to ozone, the third - can't put
!     !  W_ESFT_SW(1,1,3) as here it is dimensioned by NPD_ESFT_TERM
!     !  & NPD_BAND, but its contents are arranged as if dimensioned 
!     !  by NPD_ESFT_TERM_SW & NPD_BAND_SW.
     &     .FALSE., cmessage, icode)

      If ( icode .GT. 0 ) Then
        write(6,*) ' icode GT 0 in solvar'
        Return
      Endif
      endif  ! L_SCVARY

          If (L_orog) Allocate(sol_bearing(row_length,rows))

! DEPENDS ON: solang
          Call solang(                                                  &
! input constants
     &                Sindec, seconds_since_midnight,                   &
     &                Radiation_timestep, Eq_Time,                      &
! row and column dependent constants
     &                sin_true_latitude,                                &
     &                true_longitude,                                   &
! size variables
     &                row_length*rows,                                  &
! output fields
     &                day_fraction, cos_zenith_angle, hour_angle )

! EAK
!    Cable needs cos_zenith_angle every time step
        If (L_Rad_Step) Then

! Set rounding-error size values to zero - the criterion depends
! on the frequency of full SW calculations because on the physics
! timesteps which are not SW timesteps a test has to be done to
! avoid using the unset data for such points.
! NB: Copied from UM, including comment.
          Do j = 1, rows
            Do i= 1, row_length
              If ( cos_zenith_angle(i,j) * day_fraction(i,j)  <         &
     &            ( 1.e-10 * timestep / radiation_timestep ) ) Then
                cos_zenith_angle(i,j) = 0.0
                day_fraction(i,j) = 0.0
              End If
            End Do
          End Do

! Transfer solar bearing diagnostic before it is deallocated:

          If (SW_diag%L_sol_bearing) Then
            SW_diag%sol_bearing = sol_bearing
          End If

! Calculate orography correction:

          If (L_orog) Then
            Allocate(orog_corr(row_length,rows))
! DEPENDS ON: solinc
            Call solinc(row_length, rows, cos_zenith_angle)
            Deallocate(sol_bearing)
          End If

! DEPENDS ON: ftsa
          Call ftsa (                                                   &
! input fields
     &               land_sea_mask, flandg, ice_fract,                  &
     &               T_surf, tstar_sice,                                &
     &               cos_zenith_angle, snow_depth, snow_depth_sea,      &
! max and min sea ice albedo specifications
     &               SW_alpham, SW_alphac, SW_alphab, SW_dtice,         &
     &               L_MOSES_II,L_ssice_albedo,                         &
     &               L_MOD_BARKER_ALBEDO,                               &
     &               L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a,&
     &               dt_bare,dalb_bare_wet,pen_rad_frac,SW_beta,        &
! Version of Shortwave
     &               '03A',                                             &
! size and control variables
     &               row_length*rows, row_length*rows,                  &
! output arguments
     &               land_alb,sice_alb,                                 &
     &               open_sea_albedo )

!-----------------------------------------------------------------------
! Calculate MOSES II tile albedos, then reset TILE_PTS and TILE_INDEX
! and set tile fractions to 1 if aggregate tiles are used (NTILES=1).
!-----------------------------------------------------------------------
     if ( L_MOSES_II .OR. L_CABLE ) then
! Code required to initialize and update snow soot content
        if ( L_snow_albedo ) then
          Do j = 1, rows
            Do i = 1, row_length
              soot(i,j) = 0.
            End Do
          End Do
        endif

! DEPENDS ON: tile_albedo
            call tile_albedo ( row_length*rows,                         &
     &        land_field,land_index,ntiles,tile_pts,tile_index,         &
     &        L_snow_albedo,albsoil,cos_zenith_angle,frac,lai,          &
     &        rgrain,snow_tile,soot,tstar_tile,z0_tile,                 &
     &        alb_tile,land_albedo,can_rad_mod )
            tile_frac = 0.
            if (ntiles == 1) then
              tile_pts(1) = land_field
              do l=1,land_field
                tile_frac(l,1) = 1.
                tile_index(l,1) = l
              enddo
            else
              do n=1,ntype
                do j=1,tile_pts(n)
                  l = tile_index(j,n)
                  tile_frac(l,n) = frac(l,n)
                enddo
              enddo
            endif

       ENDIF ! ( L_MOSES_II .or. l_cable) then

       IF ( l_cable ) then

! land_alb n cable_RADNum is only set on points where there is incoming
! SW radiation at the previous timestep, therefore it will be zero over some
! land points
          if(istep_cur .gt. 1) then
! DEPENDS ON: cable_rad_driver.o
                  CALL cable_rad_driver(                &
! IN atmospheric forcing
!     &           surf_down_sw, cos_zenith_angle, L_TILE_PTS,   &
     &           surf_down_sw, cos_zenith_angle,                       &
! IN soil/vegetation/land surface data :
     &           SNOW_TILE,SNOW_TMP3L,SNOW_RHO1L,TSOIL_TILE,            &
     &           ISNOW_FLG3L,                                           &
! IN  Time stepping infomation:
     &           albsoil,                            &
     &           LAND_ALBEDO,ALB_TILE,LAND_ALB                          &
     &  )
         endif  ! if if(istep_cur .gt. 1)
      endif  ! if l_cable


!  If L_VOLCTS set TRUE then
!L  Expand aerosol optical depth time series to global field, ignoring 
!L    haloes, which aren't passed down:
!   (The Old Dynamics code used decomposition variables like GLSIZE,
!    DATASTART, OFFY, ROW & P_ROWS to work out where points were - but
!    it's easier & more general (thanks, JCT) to use SIN_TRUE_LATITUDE)
! There are 4 values for each month, for the 4 equal-area 
!  quarterspheres bounded at 30' N & S, & the Equator, i.e. where
!  sin(latitude)=+-.5 & 0.  But we normally have rows right at these 
!  latitudes, where we want to set the volcanic optical depth to the
!  average - & of course the value may not be identically correct.

      if (L_VOLCTS) then
      
      write (6,*) '------------------------- '
      write (6,*) 'Warning: Volcano mod code '
      write (6,*) 'has been included but the '
      write (6,*) 'check that the run is global'
      write (6,*) 'has been disabled.  Only use'
      write (6,*) 'this FCM branch for global'
      write (6,*) 'simulations'
      write (6,*) '------------------------- '

!     If date is outside range of given data then set VOLCMASS
!     to control value.
        IF ( (PREVIOUS_TIME(1) .LT. 1850) .OR. 			&
     &       (PREVIOUS_TIME(1) .GT. 2300) ) THEN
         DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            VOLCMASS(I,J) = MASCON * VOLCTS_val ! 1.86604E-6
          ENDDO
         ENDDO
         write (6,*) 'toplevel: control year:',PREVIOUS_TIME(1), &
                  ' VOLCTS=',VOLCTS_val,' VOLCMASS=', VOLCMASS(1,1)

        ELSE
         DO J = 1, ROWS
! We assume rows are latitude rows, which could be checked & trapped
! A check is made that global is defined. So rows should be latitude
! rows.
          IF ( SIN_TRUE_LATITUDE(1,J) .LT. -.5-DELSLT ) THEN
             ROWVOL = VOLCTS(4,PREVIOUS_TIME(2),PREVIOUS_TIME(1))
            ELSE IF ( ABS(SIN_TRUE_LATITUDE(1,J)+.5) .LT. DELSLT ) THEN
             ROWVOL = .5 * ( VOLCTS(4,PREVIOUS_TIME(2),PREVIOUS_TIME(1))  &
     &                 + VOLCTS(3,PREVIOUS_TIME(2),PREVIOUS_TIME(1)) )
            ELSE IF ( SIN_TRUE_LATITUDE(1,J) .LT. -DELSLT ) THEN
             ROWVOL = VOLCTS(3,PREVIOUS_TIME(2),PREVIOUS_TIME(1))
            ELSE IF ( ABS(SIN_TRUE_LATITUDE(1,J)) .LT. DELSLT ) THEN
             ROWVOL = .5 * ( VOLCTS(3,PREVIOUS_TIME(2),PREVIOUS_TIME(1))  &
     &                 + VOLCTS(2,PREVIOUS_TIME(2),PREVIOUS_TIME(1)) )
            ELSE IF ( SIN_TRUE_LATITUDE(1,J) .LT. .5-DELSLT ) THEN
             ROWVOL = VOLCTS(2,PREVIOUS_TIME(2),PREVIOUS_TIME(1))
            ELSE IF ( ABS(SIN_TRUE_LATITUDE(1,J)-.5) .LT. DELSLT ) THEN
             ROWVOL = .5 * ( VOLCTS(2,PREVIOUS_TIME(2),PREVIOUS_TIME(1))  &
     &                 + VOLCTS(1,PREVIOUS_TIME(2),PREVIOUS_TIME(1)) )
            ELSE !   ( SIN_TRUE_LATITUDE well over .5)
             ROWVOL = VOLCTS(1,PREVIOUS_TIME(2),PREVIOUS_TIME(1))
          ENDIF
          ROWVOL = MASCON * ROWVOL
          DO I = 1, ROW_LENGTH
            VOLCMASS(I,J) = ROWVOL
          ENDDO
         ENDDO
        ENDIF

        endif   ! L_VOLCTS

!
! To obtain reproducible results independent of the
! decomposition of the domain used on an mpp machine a global
! value for the topmost cloudy layer is used.
! The original faster code can be restored by omitting this code.
! The results will then not be independent of the number of
! segments or the configuration of processors used.
!
! DEPENDS ON: r2_global_cloud_top
          Call r2_global_cloud_top(row_length*rows, model_levels,       &
     &                              cloud_levels, cca, cct,             &
     &                              area_cloud_fraction,                &
!                       Calculated top of cloud fields.
     &                              global_cloud_top,                   &
     &                              row_length*rows )

!        GLOBAL_CLOUD_TOP returned from R2_GLOBAL_CLOUD_TOP is
!        the cloud top for the local domain. Derive a Global
!        value from the local values.
          Call gc_imax(1, n_proc, info, global_cloud_top)

          If (L_rad_deg) Then
!
            INTERP_INDEX=MOD((timestep_number/A_SW_RADSTEP),2)
!
            If ( (interp_index == 0) .EQV. rad_mask(1,1)) Then
              first_data_interp = 1
            Else
              first_data_interp = 0
            Endif
! Save interpolation parameter for later use with isccp diagnostics
            first_data_interp_sw = first_data_interp
!
            Do j=1, rows
              Do i=1, row_length
                switch(i,j) = ( (day_fraction(i,j) > 0.0) .AND.         &
     &                          (rad_mask(i,j) .EQV.                    &
     &                          (interp_index == 0) ) )
              Enddo
            Enddo
!
          Else   ! Spatial degradation is switched off
!
            Do j = 1, rows
              Do i = 1, row_length
                switch(i,j) = day_fraction(i,j) >  0.
              Enddo
            Enddo
!
          Endif   ! If l_rad_deg
!
          If (model_domain  ==  mt_global) Then
! If at a pole then only calculate first point of polar row if all
! are lit. ie set switch to false for all other points.
            If (at_extremity(PNorth) .AND.                              &
     &           ( day_fraction(1,rows) > 0.0 ) ) Then
              switch(1,rows) = .TRUE.
              Do i = 2, row_length
                switch(i,rows) = .false.
              End Do
            End If
            If (at_extremity(PSouth) .AND.                              &
     &           ( day_fraction(1,1) > 0.0 ) ) Then
              switch(1,1) = .TRUE.
              Do i = 2, row_length
                switch(i,1) = .false.
              End Do
            End If
          End If

          daylight_points = 0
          Do j = 1, rows
            Do i = 1, row_length
              If (switch(i,j)) Then
                daylight_points = daylight_points + 1
                list_daylight_points(daylight_points) =                 &
     &                                           (j-1)*row_length + i
!               The following arrays are 2-D analogues of the above
!               used for diagnostic purposes from 5.3.
                diag_row_list(daylight_points) = j
                diag_col_list(daylight_points) = i
                diag_row_list_sw(daylight_points) = j
                diag_col_list_sw(daylight_points) = i
              End If
            End Do
          End Do

! Add up the total number of daylight points in the whole domain using
! the GC routine ISUM (integer sum).
          IF (L_rad_deg) THEN
            tot_daylight_points=daylight_points
            CALL GC_ISUM(1,n_proc,Error_code,tot_daylight_points)
            If ( error_code /= 0 ) Then
! DEPENDS ON: ereport
              Call Ereport('NI_rad_ctl', Error_code,                    &
     &          'Unable to determine total lit points.')
            Endif
          ENDIF   ! If l_rad_deg

!         Zero the calculated outputs: this will later simplify the
!         treatment of the case when there are no lit points.
          SW_incs(:,:,:)=0.0
          If (SW_diag%L_direct_par) dirpar_inc(:,:) = 0.0
          netSW (:,:) = 0.0
          SWsea (:,:) = 0.0
          If (L_EXTRA_TOP_SW) Then
            top_absorption(:,:) = 0.0
          Endif
          surf_down_sw (:,:,:) = 0.0
!
!
          If ( daylight_points  >   0 ) Then
! Calculate length of segments into which to split
! short wave calculations.
            If ( A_SW_seg_size > 0 ) Then
              ! We are specifying size of segments
              step = A_SW_seg_size
              A_SW_segments = Ceiling ( Real(daylight_points) /          &
                                        Real(step) )
              N_SW_segments = A_SW_segments
            Else
              ! The original method
              N_SW_segments = min(A_SW_segments, daylight_points)
              step = daylight_points/N_SW_segments
            End If

            ! Allocate space for segmentation arrays
            Allocate(    first_point_temp( A_SW_segments ) )
            Allocate( first_point_temp_sw( A_SW_segments) )
            Allocate(     seg_points_temp( A_SW_segments) )
            Allocate(  seg_points_temp_sw( A_SW_segments) )

            seg_start=1

! Short wave radiation calculations called in segments
! Code left over from parallel version.
            first_point_temp(1) = 1
            Do i = 1, N_SW_segments
              lit_points = step
              start_point = 1 + (i-1) *step
              seg_points = list_daylight_points(i*step) - seg_start + 1
              If ( i  ==  N_SW_segments ) Then
                lit_points = daylight_points - step* (N_SW_segments - 1)
                seg_points = row_length*rows-seg_start+1
              End If
              Do j = start_point, start_point+lit_points-1
                list_daylight_points(j) = list_daylight_points(j)       &
     &                                  - seg_start + 1
              End Do
              seg_points_temp(i) = seg_points
              If ( i  <   N_SW_segments ) Then
                first_point_temp(i+1) = first_point_temp(i)+seg_points
              End If
              seg_start = seg_start + seg_points
            End Do
            Do i = 1, N_SW_segments

! Save values for later use in call to ISCCP code
              first_point_temp_sw(i)=first_point_temp(i)
              seg_points_temp_sw(i)=seg_points_temp(i)
              step_sw=step

              lit_points = step
              start_point = 1+(i-1)*step
              If ( i  ==  N_SW_segments ) Then
                lit_points = daylight_points - step* (N_SW_segments - 1)
              End If
              first_point = first_point_temp(i)

!           Set the first point of the dust arrays to be used.
            IF (L_USE_DUST) THEN
              FIRST_POINT_DUST=FIRST_POINT
            ELSE
              FIRST_POINT_DUST=1
            ENDIF

!           Set the first point of the biogenic array.
            IF (L_USE_BIOGENIC) THEN
              FIRST_POINT_BIOGENIC=FIRST_POINT
            ELSE
              FIRST_POINT_BIOGENIC=1
            ENDIF

!           Set the first point of the array of sulphate to be used.
!           A separate assignment is necessary since this array will
!           not be of the full size unless the sulphur cycle is on.
            IF (L_USE_SULPC_DIRECT .OR. L_USE_SULPC_INDIRECT_SW) THEN
               first_point_sulpc=first_point
            ELSE
               first_point_sulpc=1
            ENDIF
!
            If (L_use_soot_direct) then
              first_point_soot=first_point
            Else
              first_point_soot=1
            Endif
!
            If (L_use_bmass_direct .or. L_use_bmass_indirect) then
              first_point_biomass=first_point
            Else
              first_point_biomass=1
            EndIf
!
            If (L_use_ocff_direct .or. L_use_ocff_indirect) then
              first_point_ocff=first_point
            Else
              first_point_ocff=1
            End If
!
              If (L_use_seasalt_indirect .OR. L_use_seasalt_direct) Then
                first_point_seasalt=first_point
              Else
                first_point_seasalt=1
              Endif
!
              If (n_arcl_species > 0) Then
                first_point_arcl = first_point
              Else
                first_point_arcl = 1
              End If
!
!             Set the actual size of arrays in the radiation code:
!             for some architectures (such as that of Cray vector
!             machines) on odd size is preferred to avoid memory
!             bank conflicts.
              nd_rad_pts=2*(lit_points/2)+1
!
!             Set the number of layers seen in the radiation code.
!             This may optionally be 1 greater than the number used
!             in the rest of the model to avoid spurious effects
!             resulting from the upper boundary (principally in
!             stratospheric configurations).
              if (L_EXTRA_TOP_SW) then
                n_rad_layers=model_levels+1
              else
                n_rad_layers=model_levels
              endif
!
! DEPENDS ON: r2_swrad
              Call r2_swrad(error_code,                                 &
!                       Mixing Ratios
     &        q_n(first_point,1,1), CO2_MMR,ozone(first_point,1,1),     &
     &        O2MMR,CO2_DIM1,CO2_DIM2,CO2_3D(first_point,1,1),L_CO2_3D, &
     &        L_use_stochem_CH4, CH4_stochem(first_point,1,1),          &
!                       Pressures and Temperatures
     &        p_star(first_point,1),p_layer_boundaries(first_point,1,0),&
     &        p_layer_centres(first_point,1,0),T_n(first_point,1,1),    &
!                       Options for treating clouds
     &        .true., global_cloud_top, L_microphy,                     &
     & L_INHOM_CLOUD, INHOM_CLOUD_SW, DP_CORR_STRAT, DP_CORR_CONV,      &
     &        sin_true_latitude(first_point,1),                         &
!            Stratiform Cloud Fields
     &        L_cloud_water_partition, L_PC2,                           &
     & area_cloud_fraction(first_point,1,1),cf_n(first_point,1,1),      &
     &        qcl_n(first_point,1,1), qcf_n(first_point,1,1),           &

!                       Convective Cloud Fields
     &        cca(first_point,1,1), cclwp(first_point,1),               &
     &        ccw(first_point,1,1), lcbase(first_point,1),              &
     &        ccb(first_point,1), cct(first_point,1),                   &

!                       Surface Fields
     &        land_albedo(first_point,1,1),L_MOSES_II,L_CABLE, L_CTILE, &
     &        L_USE_SPEC_SEA,  &
     &        land_alb(first_point,1), sice_alb(first_point,1),         &
     &        flandg(first_point,1), open_sea_albedo(first_point,1,1),  &
     &        ice_fract(first_point,1), land_sea_mask(first_point,1),   &
     &        land0p5(first_point,1), snow_depth(first_point,1),        &
!                       Solar Fields
     & cos_zenith_angle(first_point,1),day_fraction(first_point,1),     &
     & list_daylight_points(start_point), SCS, l_climat_aerosol,        &
     & l_clim_aero_hgt, L_HadGEM1_Clim_Aero, zh, aero_bl_levels,        &
     &        L_use_clearrh,L_USE_DUST,DUST_DIM1,DUST_DIM2,             &
     &        DUST_1(FIRST_POINT_DUST,1),DUST_2(FIRST_POINT_DUST,1),    &
     &        DUST_3(FIRST_POINT_DUST,1),DUST_4(FIRST_POINT_DUST,1),    &
     &        DUST_5(FIRST_POINT_DUST,1),DUST_6(FIRST_POINT_DUST,1),    &
     &        L_USE_BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2,             &
     &        local_biogenic(first_point_biogenic, 1),                  &
     &        L_use_sulpc_direct, L_use_sulpc_indirect_SW,              &
     &        sulp_dim1, sulp_dim2,accum_sulphate(first_point_sulpc, 1),&
     &        aitken_sulphate(first_point_sulpc, 1),                    &
     &        diss_sulphate(first_point_sulpc, 1),			&
     &        L_volcts, volcmass(first_point,1),		        &
     &        sea_salt_film(first_point_seasalt,1,1),                   &
     &        sea_salt_jet(first_point_seasalt,1,1),                    &
     &        L_use_seasalt_indirect, L_use_seasalt_direct,             &
     &        salt_dim1*salt_dim2, salt_dim3, L_use_soot_direct,        &
     &        soot_dim1, soot_dim2, fresh_soot(first_point_soot, 1),    &
     &        aged_soot(first_point_soot, 1), L_use_bmass_direct,       &
     &        bmass_dim1,bmass_dim2,fresh_bmass(first_point_biomass,1), &
     &        aged_bmass(first_point_biomass, 1),                       &
     &        cloud_bmass(first_point_biomass, 1), L_use_bmass_indirect,&
     &        L_use_ocff_direct, ocff_dim1, ocff_dim2,                  &
     & fresh_ocff(first_point_ocff, 1), aged_ocff(first_point_ocff, 1), &
     &        cloud_ocff(first_point_ocff, 1), L_use_ocff_indirect,     &
     &        L_USE_ARCL, arcl_dim1, arcl_dim2, n_arcl_species,         &
     &        n_arcl_compnts, i_arcl_compnts,                           &
     &        local_arcl(first_point_arcl,1,1),                         &
     &        aerosol(first_point,1,1), L_MURK_RAD,                     &
     &        Ntot_land, Ntot_sea, trindx(first_point,1)                &
!                       Spectrum
!     ------------------------------------------------------------------
!     ARGUMENT LIST FOR THE REDUCED SW SPECTRAL FILE.
!     (NOTE: SWSPDC3A, SWSPCM3A AND SWSARG3A MUST BE CONSISTENT)
!
     &   , NPD_TYPE_SW, NPD_BAND_SW, NPD_EXCLUDE_SW                     &
     &   , NPD_SPECIES_SW, NPD_ESFT_TERM_SW, NPD_SCALE_FNC_SW           &
     &   , NPD_SCALE_VARIABLE_SW, NPD_THERMAL_COEFF_SW                  &
     &   , NPD_SURFACE_SW, NPD_ALBEDO_PARM_SW, NPD_CONTINUUM_SW         &
     &   , NPD_DROP_TYPE_SW, NPD_ICE_TYPE_SW, NPD_CLOUD_PARAMETER_SW    &
     &   , NPD_AEROSOL_SPECIES_SW, NPD_HUMIDITIES_SW, NPD_AOD_WAVEL_SW  &
     &   , L_PRESENT_SW, N_BAND_SW, WAVE_LENGTH_SHORT_SW                &
     &   , WAVE_LENGTH_LONG_SW, N_BAND_EXCLUDE_SW, INDEX_EXCLUDE_SW     &
     &   , SOLAR_FLUX_BAND_SW, RAYLEIGH_COEFFICIENT_SW                  &
     &   , N_ABSORB_SW, N_BAND_ABSORB_SW, INDEX_ABSORB_SW               &
     &   , TYPE_ABSORB_SW, I_BAND_ESFT_SW, I_SCALE_ESFT_SW              &
     &   , I_SCALE_FNC_SW, K_ESFT_SW, W_ESFT_SW                         &
     &   , SCALE_VECTOR_SW, P_REFERENCE_SW, T_REFERENCE_SW              &
     &   , N_DEG_FIT_SW, THERMAL_COEFFICIENT_SW, T_REF_PLANCK_SW        &
     &   , I_SPEC_SURFACE_SW, N_DIR_ALBEDO_FIT_SW, L_SURFACE_SW         &
     &   , SURFACE_ALBEDO_SW, DIRECT_ALBEDO_PARM_SW                     &
     &   , EMISSIVITY_GROUND_SW, N_BAND_CONTINUUM_SW, INDEX_CONTINUUM_SW&
     &   , INDEX_WATER_SW, I_SCALE_FNC_CONT_SW, K_CONTINUUM_SW          &
     &   , SCALE_CONTINUUM_SW, P_REF_CONTINUUM_SW, T_REF_CONTINUUM_SW   &
     &   , I_DROP_PARAMETRIZATION_SW, L_DROP_TYPE_SW                    &
     &   , DROP_PARAMETER_LIST_SW                                       &
     &   , DROP_PARM_MIN_DIM_SW, DROP_PARM_MAX_DIM_SW                   &
     &   , N_AEROSOL_SW, TYPE_AEROSOL_SW, I_AEROSOL_PARAMETRIZATION_SW  &
     &   , NHUMIDITY_SW, HUMIDITIES_SW, L_AEROSOL_SPECIES_SW            &
     &   , AEROSOL_ABSORPTION_SW, AEROSOL_SCATTERING_SW                 &
     &   , AEROSOL_ASYMMETRY_SW                                         &
     &   , I_ICE_PARAMETRIZATION_SW, L_ICE_TYPE_SW                      &
     &   , ICE_PARAMETER_LIST_SW                                        &
     &   , ICE_PARM_MIN_DIM_SW, ICE_PARM_MAX_DIM_SW                     &
     &   , L_DOPPLER_PRESENT_SW, DOPPLER_CORRECTION_SW                  &
     &   , N_AOD_WAVEL_SW, AOD_ABSORPTION_SW, AOD_SCATTERING_SW         &
     &   , I_AOD_TYPE_SW                                                &
!
!     ------------------------------------------------------------------
!                       Algorithmic options
!     ------------------------------------------------------------------
!     ARGUMENT LIST OF CONTROLLING OPTIONS FOR SHORTWAVE RADIATION.
!
     &   ,                                                              &
     &  I_2STREAM_SW, I_GAS_OVERLAP_SW, I_CLOUD_SW,                     &
     &  I_CLOUD_REPRESENTATION_SW, I_SOLVER_SW,                         &
     &  L_O2_SW,                                                        &
     &  I_ST_WATER_SW, I_CNV_WATER_SW, I_ST_ICE_SW, I_CNV_ICE_SW,       &
     &  L_LOCAL_CNV_PARTITION_SW, L_EXTRA_TOP_SW                        &
!
!     ------------------------------------------------------------------
!
     &   ,    timestep, l_mod_k_flux, L_Wenyi,                          &
!
!                       All diagnostics and associated arrays
     &        SW_diag,                                                  &
     &        diag_row_list(start_point), diag_col_list(start_point),   &
!
! Dimensions
     &        lit_points,seg_points_temp(i),model_levels,n_rad_layers,  &
     &        cloud_levels,wet_model_levels,ozone_levels,row_length,    &
     &        rows,rows*row_length,nd_rad_pts,n_rad_layers,1,           &
     &        n_cca_levels,                                             &
! Output data
!
     &        surf_down_sw(first_point,1,1),                            &
     &        flux_below_690nm_surf(first_point,1),                     &
     &        L_flux_below_690nm_surf,                                  &
     &        netsw(first_point,1), top_absorption(first_point,1),      &
     &        swsea(first_point,1), SW_incs(first_point,1,0),           &
     &        dirpar_inc(first_point,1), SW_diag%L_direct_par,          &
!
! Variables needed to calculate layer masses
     &        rho_r2_nh(first_point,1),                                 &
     &        r_rho_levels_nh(first_point,1),                           &
     &        r_theta_levels_nh(first_point,0),                         &
     &        q_n(first_point,1,1), qcl_n(first_point,1,1),             &
     &        qcf_n(first_point,1,1), qcf2_n(first_point,1,1),          &
     &        qrain_n(first_point,1,1), qgraup_n(first_point,1,1),      &
     &        l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio     &
!
     &        )


            End Do
          End If ! End If on daylight_points > 0

          If (L_orog) Deallocate(orog_corr)

!         Radiative fluxes may not have been calculated at all
!         points: we now fill in as required. In the case of
!         spatial degradation calls to SWAPBOUNDS are made and
!         these must be made on all PEs, which requires that
!         this code shall executed even when the PE contains
!         no daylit points.
!
!         At the North Pole in a global domain calculations
!         are performed only at the first point if lit, so
!         the rest of the row must be filled in if this point
!         is lit.
          L_complete_North = ( model_domain == mt_global ) .AND.        &
     &                       ( at_extremity(PNorth) ) .AND.             &
     &                       ( switch(1,rows) )

!         At the South Pole in a global domain calculations
!         are performed only at the first point if lit, so
!         the rest of the row must be filled in if this point
!         is lit.
          L_complete_South = ( model_domain == mt_global ) .AND.        &
     &                       ( at_extremity(PSouth) ) .AND.             &
     &                       ( switch(1,1) )
!
!         When spatial degradation is performed fields must be
!         filled in at alternate points.
          L_complete_deg = ( L_rad_deg ) .AND.                          &
     &                     ( tot_daylight_points > 0 )
!
!         Set addressing limits for spatial degradation.
          If ( L_complete_deg ) Then
            first_row=1
            last_row=rows
            If (model_domain  ==  mt_global) Then
              If (at_extremity(PNorth)) Then
                last_row=rows-1
              End If
              If (at_extremity(PSouth)) Then
                first_row=2
              End If
            End If
          Endif
!
!         Call appropriate subroutines to fill in missing data
!         as required.
!
          If ( L_complete_North .OR. L_complete_South .OR.              &
     &         L_complete_deg ) Then
!
!           Primary Fields:
!
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp,                       &
     &        model_levels+2,                                           &
     &        SW_incs                                                   &
     &        )
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, 1,                    &
     &        netSW                                                     &
     &        )
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, 1,                    &
     &        SWsea                                                     &
     &        )
            If ( L_flux_below_690nm_surf ) Then
!             This field is purely diagnostic in some versions of
!             the model, but is always required with MOSES.
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          flux_below_690nm_surf                                   &
     &          )
              If (SW_diag%L_direct_par) Then
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 1,                &
     &            dirpar_inc                                            &
     &            )
              Endif
            Endif
            If (L_extra_top_SW) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          top_absorption                                          &
     &          )
            Endif
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp,                       &
     &        4,                                                        &
     &        surf_down_SW                                              &
     &        )
            If (L_orog) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, 1,                    &
     &        f_orog                                                    &
     &        )
            Endif
!
!           Complete the diagnostic fields as required.
            If ( SW_diag%L_solar_out_toa ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%solar_out_toa                                   &
     &          )
            Endif
            If ( SW_diag%L_solar_out_clear ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%solar_out_clear                                 &
     &          )
            Endif
            If ( SW_diag%L_surface_down_flux ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%surface_down_flux                               &
     &          )
            Endif
            If ( SW_diag%L_surf_down_clr ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%surf_down_clr                                   &
     &          )
            Endif
            If ( SW_diag%L_surf_up_clr ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%surf_up_clr                                     &
     &          )
            Endif
            If ( SW_diag%L_clear_hr ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          model_levels,                                           &
     &          SW_diag%clear_hr                                        &
     &          )
            Endif
            If ( SW_diag%L_net_flux_trop ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%net_flux_trop                                   &
     &          )
            Endif
            If ( SW_diag%L_up_flux_trop ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%up_flux_trop                                    &
     &          )
            Endif
!
! Direct and Diffuse SW Downward Flux       
                    
            If ( SW_diag%L_flux_direct ) Then
! DEPENDS ON: rad3d_inp     
              Call rad3d_inp(                                           & 
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          model_levels+1,                                         &
     &          SW_diag%flux_direct                                     &
     &          )
            Endif
            
            If ( SW_diag%L_flux_diffuse ) Then
              Call rad3d_inp(                                           &  
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          model_levels+1,                                         &
     &          SW_diag%flux_diffuse                                    &
     &          )
            Endif           
            
!
!           Microphysical diagnostics
!
            If ( SW_diag%re_conv_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          cloud_levels,                                           &
     &          SW_diag%re_conv                                         &
     &          )
            Endif
            If ( SW_diag%re_strat_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          cloud_levels,                                           &
     &          SW_diag%re_strat                                        &
     &          )
            Endif
            If ( SW_diag%wgt_conv_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          cloud_levels,                                           &
     &          SW_diag%wgt_conv                                        &
     &          )
            Endif
            If ( SW_diag%wgt_strat_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          cloud_levels,                                           &
     &          SW_diag%wgt_strat                                       &
     &          )
            Endif
            If ( SW_diag%lwp_strat_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          cloud_levels,                                           &
     &          SW_diag%lwp_strat                                       &
     &          )
            Endif
            If ( SW_diag%weighted_re_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%weighted_re                                     &
     &          )
            Endif
            If ( SW_diag%sum_weight_re_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%sum_weight_re                                   &
     &          )
            Endif
            If ( SW_diag%wgtd_warm_re_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%weighted_warm_re                                &
     &          )
            Endif
            If ( SW_diag%sum_wgt_warm_re_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%sum_weight_warm_re                              &
     &          )
            Endif
            If ( SW_diag%Nc_diag_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%Nc_diag                                         &
     &          )
            Endif
            If ( SW_diag%Nc_weight_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%Nc_weight                                       &
     &          )
            Endif
            If ( SW_diag%ntot_diag_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          cloud_levels,                                           &
     &          SW_diag%ntot_diag                                       &
     &          )
            Endif
            If ( SW_diag%strat_lwc_diag_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          cloud_levels,                                           &
     &          SW_diag%strat_lwc_diag                                  &
     &          )
            Endif
            If ( SW_diag%so4_ccn_diag_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          cloud_levels,                                           &
     &          SW_diag%so4_ccn_diag                                    &
     &          )
            Endif
            If ( SW_diag%cond_samp_wgt_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          cloud_levels,                                           &
     &          SW_diag%cond_samp_wgt                                   &
     &          )
            Endif
            If ( SW_diag%seasalt_film_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          cloud_levels,                                           &
     &          sea_salt_film                                           &
     &          )
            Endif
            If ( SW_diag%seasalt_jet_flag ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp,                     &
     &          cloud_levels,                                           &
     &          sea_salt_jet                                            &
     &          )
            Endif
!
!     Diagnostics for MOSES
!
            If ( SW_diag%L_FlxSolBelow690nmSurf ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%FlxSolBelow690nmSurf                            &
     &          )
            Endif
            If ( SW_diag%L_FlxSeaBelow690nmSurf ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%FlxSeaBelow690nmSurf                            &
     &          )
            Endif
!
!     Diagnostics for orography correction

            If ( SW_diag%L_orog_corr ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, 1,                  &
     &          SW_diag%orog_corr                                       &
     &          )
            Endif
!
!        Extinction diagnostics:
!
            If ( SW_diag%L_cloud_extinction ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, cloud_levels,       &
     &          SW_diag%cloud_extinction                                &
     &          )
            Endif
           If ( SW_diag%L_cloud_weight_extinction ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, cloud_levels,       &
     &          SW_diag%cloud_weight_extinction                         &
     &          )
            Endif
            If ( SW_diag%L_ls_cloud_extinction ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, cloud_levels,       &
     &          SW_diag%ls_cloud_extinction                             &
     &          )
            Endif
           If ( SW_diag%L_ls_cloud_weight_extinction ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, cloud_levels,       &
     &          SW_diag%ls_cloud_weight_extinction                      &
     &          )
            Endif
            If ( SW_diag%L_cnv_cloud_extinction ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, cloud_levels,       &
     &          SW_diag%cnv_cloud_extinction                            &
     &          )
            Endif
           If ( SW_diag%L_cnv_cloud_weight_extinction ) Then
! DEPENDS ON: rad3d_inp
              Call rad3d_inp(                                           &
     &          L_complete_North, L_complete_South, L_complete_deg,     &
     &          row_length, rows, off_x, off_y, first_row, last_row,    &
     &          first_data_interp, ES_space_interp, cloud_levels,       &
     &          SW_diag%cnv_cloud_weight_extinction                     &
     &          )
            Endif
!



          End If ! on number of daylight points.
!
          Do j = 1, rows
            Do i= 1, row_length
              mean_cos_zenith_angle(i,j) = cos_zenith_angle(i,j) *      &
     &                                     day_fraction(i,j)
            End Do
          End Do
! calculate net surface SW for diagnostic
          If (L_CTILE) then
            Do j = 1, rows
              Do i= 1, row_length
                surfsw(i,j) = SW_incs(i,j,0)                            &
     &                      * mean_cos_zenith_angle(i,j)                &
     &                      + (1.-flandg(i,j))*SWsea(i,j)
!                land_alb(i,j)=0.0
                if(flandg(i,j) == 1.0)SWsea(i,j)=rmdi
              End Do
            End Do
          Else
            Do j = 1, rows
              Do i= 1, row_length
                surfsw(i,j) = SW_incs(i,j,0)                            &
     &                      * mean_cos_zenith_angle(i,j)                &
     &                      + SWsea(i,j)
                if(land_sea_mask(i,j))SWsea(i,j)=rmdi
              End Do
            End Do
          Endif
          If (L_MOSES_II) then

! calculate net surface SW on tiles if required for MOSES II
            Do l=1,land_field
              i = land_index_i(l)
              j = land_index_j(l)
              sw_down_rts(i,j) = 0.0
              sw_net_rts(i,j) = 0.0
            End Do
            Do n=1,ntiles
              Do l=1,land_field
                SW_tile_rts(l,n) = 0.
                SW_tile_down_rts(l,n) = 0.
              End Do
            End Do
            Do n=1,ntiles
              Do point=1,tile_pts(n)
                l = tile_index(point,n)
                i = land_index_i(l)
                j = land_index_j(l)
                Do k=1,4
                  SW_tile_rts(l,n) = SW_tile_rts(l,n) +                 &
     &                        (1. - alb_tile(l,n,k))*surf_down_sw(i,j,k)
                  SW_tile_down_rts(l,n) = SW_tile_down_rts(l,n) +       &
     &                                               surf_down_sw(i,j,k)
                End Do
                sw_net_rts(i,j)=sw_net_rts(i,j) +                       &
     &            SW_tile_rts(l,n)*tile_frac(l,n)
              End Do
            End Do
            IF(L_CTILE)THEN
              Do l=1,land_field
                i = land_index_i(l)
                j = land_index_j(l)
                Do k=1,4
                  sw_down_rts(i,j) = sw_down_rts(i,j) +                 &
     &                  surf_down_sw(i,j,k)
                Enddo
                ! calculate land_alb only on tiles with
                ! non-zero absorption
                ! this fixes crash below when calculating
                ! sw_net_land for coastal tiles case
                if(sw_down_rts(i,j) >  0.0) then
                   tmp_abso=sw_net_rts(i,j)/sw_down_rts(i,j)
                   if(tmp_abso > 0.0)                                   &
     &                land_alb(i,j)=1.0-tmp_abso
                end if
              End Do
            Endif
          End If

        End If ! on a radiation timestep

! Calculate day fraction and mean cos(solar zenith angle while
! the sun is up) for each grid point for this physics timestep:
! (if in fact full SW calculations are being done every timestep, this
! is of course unnecessary, as are various calculations later on)

        If ( Radiation_timestep  >   timestep ) Then

! DEPENDS ON: solang
          Call solang(                                                  &
! input constants
     &                Sindec, seconds_since_midnight,                   &
     &                timestep, Eq_Time,                                &
! row and column dependent constants
     &                sin_true_latitude,                                &
     &                true_longitude,                                   &
! size variables
     &                row_length*rows,                                  &
! output fields
     &                day_fraction, cos_zenith_angle, hour_angle )

        End If

!  write(6,*) 'XXX glue_rad-rad_ctl2 (GLUE_RAD):  SC = ', SC

! Combine the two terms to give the mean cos zenith angle over the
! whole of the physics timestep.
! calculate incoming SW at top of atmosphere
        Do j = 1, rows
          Do i= 1, row_length
            cos_zenith_angle(i,j) = cos_zenith_angle(i,j) *             &
     &                              day_fraction(i,j)
            itoasw(i,j) = cos_zenith_angle(i,j) * SCS * SC
          End Do
        End Do
!

! Is the PC2 cloud scheme being used?
!
        If (L_pc2) then
!
! Reset _latest values to _n values
!
          Do k = 1,model_levels
            Do j = 1, rows
              Do i = 1, row_length
                delta_T(i,j,k) = SW_incs(i,j,k) * cos_zenith_angle(i,j)
                T_latest(i,j,k)   = T_n(i,j,k)
              End Do
            End Do
          End Do
          Do k = 1,wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_latest(i,j,k)   = q_n(i,j,k)
                qcl_latest(i,j,k) = qcl_n(i,j,k)
                cf_latest(i,j,k)  = cf_n(i,j,k)
                cfl_latest(i,j,k) = cfl_n(i,j,k)
              End Do
            End Do
          End Do
!
! ----------------------------------------------------------------------
! Homogeneous forcing. Note the temperature increment from the shortwave
! is added in this routine
! ----------------------------------------------------------------------
!
! DEPENDS ON: pc2_homog_plus_turb
          CALL PC2_HOMOG_PLUS_TURB(p_layer_centres(1,1,1),              &
     &      wet_model_levels,                                           &
     &      row_length, rows, timestep, T_latest, cf_latest, cfl_latest,&
     &      cff_latest, q_latest, qcl_latest, delta_t(1,1,1),           &
     &      zeros, zeros, zeros, 0.0, 0.0, l_mixing_ratio)
!
! Add increments from the homogeneous forcing to the increment variables
!
          Do k = 1,model_levels
            Do j = 1, rows
              Do i = 1, row_length
                T_inc(i,j,k) = T_inc(i,j,k) + T_latest(i,j,k)-T_n(i,j,k)
              End Do
            End Do
          End Do
          Do k = 1,wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_inc(i,j,k) = q_inc(i,j,k) + q_latest(i,j,k)-q_n(i,j,k)
                qcl_inc(i,j,k) = qcl_inc(i,j,k)                         &
     &                           + qcl_latest(i,j,k)-qcl_n(i,j,k)
                cf_inc(i,j,k) = cf_inc(i,j,k)                           &
     &                           + cf_latest(i,j,k)-cf_n(i,j,k)
                cfl_inc(i,j,k) = cfl_inc(i,j,k)                         &
     &                           + cfl_latest(i,j,k)-cfl_n(i,j,k)
              End Do
            End Do
          End Do
!
        Else  ! L_pc2
!
! add SW radiative heating to temperatures
          Do k = 1,model_levels
            Do j = 1, rows
              Do i = 1, row_length
                delta_T(i,j,k) = SW_incs(i,j,k) * cos_zenith_angle(i,j)
                T_inc(i,j,k) = T_inc(i,j,k) + delta_T(i,j,k)
                T_latest(i,j,k) = T_n(i,j,k) + delta_T(i,j,k)
              End Do
            End Do
          End Do
!
        End If  ! L_pc2
!
! Get T_incr for output as STASH diagnostic
      If ( ( L_T_incr_sw ).or.                                          &
     &     ( sf(232,1) ) ) Then
!       Increments will be calculated if to be diagnosed directly
!       or to be used to determine heating rates.
        Allocate ( T_incr_diagnostic(row_length,rows,model_levels) )
        Do k=1,model_levels
         Do j=1,rows
          Do i=1,row_length
            T_incr_diagnostic(i,j,k) =  delta_T(i,j,k)
          Enddo ! i
         Enddo ! j
        Enddo ! k
!
      Else
        Allocate ( T_incr_diagnostic(1,1,1) )
      Endif                    ! on STASHflag

! Set up net down surface SW radiation flux in SURF_RADFLUX

! Only set on sea-ice points if MOSES2:
       If (L_CTILE) then
        Do j = 1, rows
          Do i= 1, row_length

! land_alb is only set on points where there is incoming sw radiation
! at the previous timestep, therefore it will be zero over some
! land points
            If (fracsolid(i,j) >  0.0)then

              If (flandg(i,j) >  0.0.and.land_alb(i,j) <= 0.0)then
                sw_net_land(i,j) = 0.0
                sw_net_sice(i,j) = 0.0
              Else
                albsolid = ( flandg(i,j) * land_alb(i,j) +              &
     &            (1.0-flandg(i,j)) * ice_fract(i,j) * sice_alb(i,j) )  &
     &            /fracsolid(i,j)

                if(flandg(i,j) >  0.0)sw_net_land(i,j) = SW_incs(i,j,0) &
     &            * cos_zenith_angle(i,j) / fracsolid(i,j)              &
     &            * (1.0-land_alb(i,j))/(1.0-albsolid)

                if(ice_fract(i,j) >  0.0)then
                  sw_net_sice(i,j) = SW_incs(i,j,0)                     &
     &            * cos_zenith_angle(i,j) / fracsolid(i,j)              &
     &            * (1.0-sice_alb(i,j))/(1.0-albsolid)

                  surf_radflux(i,j) = sw_net_sice(i,j) * ice_fract(i,j)
                Endif
              Endif

            Endif
          End Do
        End Do
      Else
        Do j = 1, rows
          Do i= 1, row_length
            surf_radflux(i,j) = SW_incs(i,j,0) * cos_zenith_angle(i,j)
          End Do
        End Do
      Endif


! Set up photosynthetically active surface radiation, if calculated

        If (l_flux_below_690nm_surf) Then
          Do j = 1, rows
            Do i= 1, row_length
              photosynth_act_rad(i,j) = SW_incs(i,j,model_levels+1)     &
     &                                  * cos_zenith_angle(i,j)
            End Do
          End Do
        If (SW_diag%L_direct_par) Then
          Do j = 1, rows
            Do i = 1, row_length
              SW_diag%flxdirparsurf(i,j) = dirpar_inc(i,j) *            &
     &          cos_zenith_angle(i,j)
            End Do
          End Do
        End If
        End If

! Set up net surface SW on tiles if required for MOSES II
          if (L_MOSES_II) then
            If (ntiles == 1) then
              Do l=1,land_field
                i = land_index_i(l)
                j = land_index_j(l)
                SW_tile(l,1) = SW_tile_rts(l,1)                         &
     &                         * cos_zenith_angle(i,j)
                SW_tile_down(l,1) = SW_tile_down_rts(l,1)               &
     &                         * cos_zenith_angle(i,j)
              End Do
            Else
              do n=1,ntiles
                do point=1,tile_pts(n)
                  l = tile_index(point,n)
                  i = land_index_i(l)
                  j = land_index_j(l)
                  SW_tile(l,n) = SW_tile_rts(l,n)                       &
     &                           * cos_zenith_angle(i,j)
                  SW_tile_down(l,n) = SW_tile_down_rts(l,n)             &
     &                           * cos_zenith_angle(i,j)
                End Do
              End Do
            End If
          End If
! DEPENDS ON: timer
        If (Ltimer) Call timer ('SW Rad  ',4)

! ----------------------------------------------------------------------
! Section RAD.1.1 Short Wave Radiation Energy correction code
!-----------------------------------------------------------------------

! DEPENDS ON: timer
        If (Ltimer) Call timer ('Eng Corr',3)

        If (L_Rad_Step .and. L_emcorr) Then

! Sum short wave fluxes into the atmosphere and
! add into the net diabatic fluxes into the
! atmosphere for use in the energy correction
! procedure

          If (L_EXTRA_TOP_SW) Then
!           The energy absorbed above the top of the model in
!           the radiation scheme does not contribute to the
!           energy absorbed, but the diagnostics are calculated
!           at the top of the atmosphere, so the net atmospheric
!           flux must be adjusted.
            Do j = 1, rows
              Do i = 1, row_length
                net_atm_flux(i,j) = netSW(i,j) - surfsw(i,j)            &
     &                            - top_absorption(i,j)
              End Do
            End Do
          Else
            Do j = 1, rows
              Do i = 1, row_length
                net_atm_flux(i,j) = netSW(i,j) - surfsw(i,j)
              End Do
            End Do
          End If

          If (L_orog) net_atm_flux = net_atm_flux + f_orog

! DEPENDS ON: flux_diag
          Call flux_diag(net_atm_flux, FV_cos_theta_latitude,           &
     &                   row_length, rows ,off_x,off_y, 1.0,            &
     &                   sum_eng_fluxes,radiation_timestep)

        End If

! DEPENDS ON: timer
        If (Ltimer) Call timer ('Eng Corr',4)

! ----------------------------------------------------------------------
! Section RAD.1.2 Short Wave Radiation diagnostics
!-----------------------------------------------------------------------

! Check that sw diagnostics requested this timestep
        If (error_code  ==  0 .and. sf(0,1)) Then
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',3)

! DEPENDS ON: diagnostics_sw
          Call diagnostics_sw(                                          &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels, cloud_levels             &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      timestep                                   &
     &,                      T_n, T_inc                                 &
     &,                      q_n, qcl_n, cf_n, cfl_n                    &
     &,                      T_latest, q_latest, qcl_latest             &
     &,                      cf_latest, cfl_latest                      &
     &,                      surfsw, itoasw, SW_diag%solar_out_toa      &
     &,                      SW_diag%solar_out_clear                    &
     &,                      SW_diag%surface_down_flux                  &
     &,                      SW_diag%surf_down_clr, SW_diag%surf_up_clr &
     &,                      SWsea, flux_below_690nm_surf               &
     &,                      photosynth_act_rad                         &
     &,                      SW_diag%flxdirparsurf                      &
     &,                      SW_diag%FlxSolBelow690nmSurf               &
     &,                      SW_diag%FlxSeaBelow690nmSurf               &
     &,                      SW_diag%orog_corr, SW_diag%sol_bearing     &
     &,                      f_orog                                     &
     &,                      slope_aspect, slope_angle                  &
     &,                      sw_net_land,sw_net_sice                    &
     &,                      T_incr_diagnostic                          &
     &,                      SW_diag%clear_hr                           &
     &,                      SW_diag%flux_direct,SW_diag%flux_diffuse   &
     &,                      SW_diag%net_flux_trop, SW_diag%up_flux_trop&
     &,                      SW_diag%re_strat, SW_diag%wgt_strat        &
     &,                      SW_diag%lwp_strat                          &
     &,                      SW_diag%re_conv, SW_diag%wgt_conv          &
     &,                      SW_diag%ntot_diag, SW_diag%strat_lwc_diag  &
     &,                      SW_diag%so4_ccn_diag, SW_diag%cond_samp_wgt&
     &,                      SW_diag%weighted_re, SW_diag%sum_weight_re &
     &,                      SW_diag%weighted_warm_re                   &
     &,                      SW_diag%sum_weight_warm_re                 &
     &,                      SW_diag%Nc_diag, SW_diag%Nc_weight         &
     &,                      sea_salt_film, sea_salt_jet                &
     &,                      salt_dim1, salt_dim2, salt_dim3            &
     &,                      SW_diag%cloud_extinction                   &
     &,                      SW_diag%cloud_weight_extinction            &
     &,                      SW_diag%ls_cloud_extinction                &
     &,                      SW_diag%ls_cloud_weight_extinction         &
     &,                      SW_diag%cnv_cloud_extinction               &
     &,                      SW_diag%cnv_cloud_weight_extinction        &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork1                                                       &
     & )

! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',4)
        Endif ! on error_code .and. sf(0,1)


      Deallocate ( T_incr_diagnostic )
      If (L_orog) Then
        Deallocate(f_orog)
      Else
        Deallocate(f_orog,slope_aspect,slope_angle)
      Endif

!     Deallocate the SW diagnostic space which is no longer required.
      If (SW_diag%L_solar_out_toa) Deallocate(SW_diag%solar_out_toa)
      If (SW_diag%L_solar_out_clear)                                    &
     &  Deallocate(SW_diag%solar_out_clear)
      If (SW_diag%L_surface_down_flux)                                  &
     &  Deallocate(SW_diag%surface_down_flux)
      If (SW_diag%L_surf_down_clr) Deallocate(SW_diag%surf_down_clr)
      If (SW_diag%L_surf_up_clr) Deallocate(SW_diag%surf_up_clr)
      If (SW_diag%L_net_flux_trop) Deallocate(SW_diag%net_flux_trop)
      If (SW_diag%L_up_flux_trop) Deallocate(SW_diag%up_flux_trop)
      If (SW_diag%L_flux_direct) Deallocate(SW_diag%flux_direct)
      If (SW_diag%L_flux_diffuse) Deallocate(SW_diag%flux_diffuse)
      If (SW_diag%L_clear_hr) Deallocate(SW_diag%clear_hr)
! Microphysical diagnostics
      If (SW_diag%re_strat_flag) Deallocate(SW_diag%re_strat)
      If (SW_diag%wgt_strat_flag) Deallocate(SW_diag%wgt_strat)
      If (SW_diag%lwp_strat_flag) Deallocate(SW_diag%lwp_strat)
      If (SW_diag%re_conv_flag) Deallocate(SW_diag%re_conv)
      If (SW_diag%wgt_conv_flag) Deallocate(SW_diag%wgt_conv)
      If (SW_diag%ntot_diag_flag) Deallocate(SW_diag%ntot_diag)
      If (SW_diag%strat_lwc_diag_flag)                                  &
     &  Deallocate(SW_diag%strat_lwc_diag)
      If (SW_diag%so4_ccn_diag_flag) Deallocate(SW_diag%so4_ccn_diag)
      If (SW_diag%cond_samp_wgt_flag)                                   &
     &  Deallocate(SW_diag%cond_samp_wgt)
      If (SW_diag%weighted_re_flag) Deallocate(SW_diag%weighted_re)
      If (SW_diag%sum_weight_re_flag)                                   &
     &  Deallocate(SW_diag%sum_weight_re)
      If (SW_diag%wgtd_warm_re_flag)                                    &
     &  Deallocate(SW_diag%weighted_warm_re)
      If (SW_diag%sum_wgt_warm_re_flag)                                 &
     &  Deallocate(SW_diag%sum_weight_warm_re)
      If (SW_diag%Nc_diag_flag) Deallocate(SW_diag%Nc_diag)
      If (SW_diag%Nc_weight_flag) Deallocate(SW_diag%Nc_weight)
! Diagnostics for MOSES II
      If (SW_diag%L_FlxSolBelow690nmSurf)                               &
     &  Deallocate(SW_diag%FlxSolBelow690nmSurf)
      If (SW_diag%L_FlxSeaBelow690nmSurf)                               &
     &  Deallocate(SW_diag%FlxSeaBelow690nmSurf)
! STOCHEM diagnostics
      If (SW_diag%L_direct_par)                                         &
     &  Deallocate(SW_diag%flxdirparsurf)
! Diagnostics for orography correction
      If (SW_diag%L_orog_corr) Deallocate(SW_diag%orog_corr)
      If (SW_diag%L_sol_bearing) Deallocate(SW_diag%sol_bearing)
!
! Note: the extinction diagnostics are deallocated later, as they
! may be needed to calculate the ISCCP diagnostics.
!
!
!     Carry out deallocation in case when diagnostics were not
!     selected. This is required only because minimal space was
!     allocated initially for running on certain compilers.
      If (.NOT.SW_diag%L_solar_out_toa)                                 &
     &  Deallocate(SW_diag%solar_out_toa)
      If (.NOT.SW_diag%L_solar_out_clear)                               &
     &  Deallocate(SW_diag%solar_out_clear)
      If (.NOT.SW_diag%L_surface_down_flux)                             &
     &  Deallocate(SW_diag%surface_down_flux)
      If (.NOT.SW_diag%L_surf_down_clr)                                 &
     &  Deallocate(SW_diag%surf_down_clr)
      If (.NOT.SW_diag%L_surf_up_clr)                                   &
     &  Deallocate(SW_diag%surf_up_clr)
      If (.NOT.SW_diag%L_net_flux_trop)                                 &
     &  Deallocate(SW_diag%net_flux_trop)
      If (.NOT.SW_diag%L_up_flux_trop)                                  &
     &  Deallocate(SW_diag%up_flux_trop)
      If (.NOT.SW_diag%L_flux_direct) Deallocate(SW_diag%flux_direct)
      If (.NOT.SW_diag%L_flux_diffuse) Deallocate(SW_diag%flux_diffuse)
      If (.NOT.SW_diag%L_clear_hr) Deallocate(SW_diag%clear_hr)
! Microphysical diagnostics
      If (.NOT.SW_diag%re_strat_flag) Deallocate(SW_diag%re_strat)
      If (.NOT.SW_diag%wgt_strat_flag)                                  &
     &  Deallocate(SW_diag%wgt_strat)
      If (.NOT.SW_diag%lwp_strat_flag)                                  &
     &  Deallocate(SW_diag%lwp_strat)
      If (.NOT.SW_diag%re_conv_flag) Deallocate(SW_diag%re_conv)
      If (.NOT.SW_diag%wgt_conv_flag) Deallocate(SW_diag%wgt_conv)
      If (.NOT.SW_diag%ntot_diag_flag)                                  &
     &  Deallocate(SW_diag%ntot_diag)
      If (.NOT.SW_diag%strat_lwc_diag_flag)                             &
     &  Deallocate(SW_diag%strat_lwc_diag)
      If (.NOT.SW_diag%so4_ccn_diag_flag)                               &
     &  Deallocate(SW_diag%so4_ccn_diag)
      If (.NOT.SW_diag%cond_samp_wgt_flag)                              &
     &  Deallocate(SW_diag%cond_samp_wgt)
      If (.NOT.SW_diag%weighted_re_flag)                                &
     &  Deallocate(SW_diag%weighted_re)
      If (.NOT.SW_diag%sum_weight_re_flag)                              &
     &  Deallocate(SW_diag%sum_weight_re)
      If (.NOT.SW_diag%wgtd_warm_re_flag)                               &
     &  Deallocate(SW_diag%weighted_warm_re)
      If (.NOT.SW_diag%sum_wgt_warm_re_flag)                            &
     &  Deallocate(SW_diag%sum_weight_warm_re)
      If (.NOT.SW_diag%Nc_diag_flag) Deallocate(SW_diag%Nc_diag)
      If (.NOT.SW_diag%Nc_weight_flag) Deallocate(SW_diag%Nc_weight)
! Diagnostics for MOSES II
      If (.NOT.SW_diag%L_FlxSolBelow690nmSurf)                          &
     &  Deallocate(SW_diag%FlxSolBelow690nmSurf)
      If (.NOT.SW_diag%L_FlxSeaBelow690nmSurf)                          &
     &  Deallocate(SW_diag%FlxSeaBelow690nmSurf)
! STOCHEM diagnostics
      If (.NOT.SW_diag%L_direct_par)                                    &
     &  Deallocate(SW_diag%flxdirparsurf)
! Diagnostics for orography correction
      If (.Not.SW_diag%L_orog_corr) Deallocate(SW_diag%orog_corr)
      If (.Not.SW_diag%L_sol_bearing) Deallocate(SW_diag%sol_bearing)


! ----------------------------------------------------------------------
! Section RAD.2 Long Wave Radiation Code.
!-----------------------------------------------------------------------

!
!     Allocate space for the diagnostic arrays required.
!
      If (LW_diag%L_total_cloud_cover)                                  &
     &  allocate(LW_diag%total_cloud_cover(row_length, rows))
      If (LW_diag%L_clear_olr)                                          &
     &  allocate(LW_diag%clear_olr(row_length, rows))
      If (LW_diag%L_surface_down_flux)                                  &
     &  allocate(LW_diag%surface_down_flux(row_length, rows))
      If (LW_diag%L_surf_down_clr)                                      &
     &  allocate(LW_diag%surf_down_clr(row_length, rows))
      If (LW_diag%L_clear_hr)                                           &
     &  allocate(LW_diag%clear_hr(row_length, rows, model_levels))
      If (LW_diag%L_net_flux_trop)                                      &
     &  allocate(LW_diag%net_flux_trop(row_length, rows))
      If (LW_diag%L_down_flux_trop)                                     &
     &  allocate(LW_diag%down_flux_trop(row_length, rows))
!
! Absorptivity diagnostics:
      If (LW_diag%L_cloud_absorptivity)                                 &
     &  allocate(LW_diag%cloud_absorptivity(row_length,                 &
     &                                        rows, cloud_levels))
      If (LW_diag%L_cloud_weight_absorptivity)                          &
     &  allocate(LW_diag%cloud_weight_absorptivity(row_length,          &
     &                                        rows, cloud_levels))
      If (LW_diag%L_ls_cloud_absorptivity)                              &
     &  allocate(LW_diag%ls_cloud_absorptivity(row_length,              &
     &                                        rows, cloud_levels))
      If (LW_diag%L_ls_cloud_weight_absorptivity)                       &
     &  allocate(LW_diag%ls_cloud_weight_absorptivity(row_length,       &
     &                                        rows, cloud_levels))
      If (LW_diag%L_cnv_cloud_absorptivity)                             &
     &  allocate(LW_diag%cnv_cloud_absorptivity(row_length,             &
     &                                        rows, cloud_levels))
      If (LW_diag%L_cnv_cloud_weight_absorptivity)                      &
     &  allocate(LW_diag%cnv_cloud_weight_absorptivity(row_length,      &
     &                                        rows, cloud_levels))
      If (LW_diag%L_total_cloud_on_levels)                              &
     &  allocate(LW_diag%total_cloud_on_levels(row_length,              &
     &                                        rows, cloud_levels))
!
! Grid-box mean cloud diagnostics as seen by radiation:
      If (LW_diag%L_ls_qcl_rad)                                         &
     &  allocate(LW_diag%ls_qcl_rad(row_length, rows, model_levels))
      If (LW_diag%L_ls_qcf_rad)                                         &
     &  allocate(LW_diag%ls_qcf_rad(row_length, rows, model_levels))
      If (LW_diag%L_cc_qcl_rad)                                         &
     &  allocate(LW_diag%cc_qcl_rad(row_length, rows, model_levels))
      If (LW_diag%L_cc_qcf_rad)                                         &
     &  allocate(LW_diag%cc_qcf_rad(row_length, rows, model_levels))
      If (LW_diag%L_ls_cl_rad)                                          &
     &  allocate(LW_diag%ls_cl_rad(row_length, rows, model_levels))
      If (LW_diag%L_ls_cf_rad)                                          &
     &  allocate(LW_diag%ls_cf_rad(row_length, rows, model_levels))
      If (LW_diag%L_cc_cl_rad)                                          &
     &  allocate(LW_diag%cc_cl_rad(row_length, rows, model_levels))
      If (LW_diag%L_cc_cf_rad)                                          &
     &  allocate(LW_diag%cc_cf_rad(row_length, rows, model_levels))
!
! Aerosol optical depth diagnostics:
      If (LW_diag%L_aod_sulphate) then
        allocate(LW_diag%aod_sulphate(row_length, rows,                 &
     &                                N_AOD_WAVEL_LW))
        LW_diag%aod_sulphate(:,:,:) = 0.0
      Endif
      If (LW_diag%L_aod_dust) then
        allocate(LW_diag%aod_dust(row_length, rows,                     &
     &                            N_AOD_WAVEL_LW))
        LW_diag%aod_dust(:,:,:) = 0.0
      Endif
      If (LW_diag%L_aod_seasalt) then
        allocate(LW_diag%aod_seasalt(row_length, rows,                  &
     &                               N_AOD_WAVEL_LW))
        LW_diag%aod_seasalt(:,:,:) = 0.0
      Endif
      If (LW_diag%L_aod_soot) then
        allocate(LW_diag%aod_soot(row_length, rows,                     &
     &                            N_AOD_WAVEL_LW))
        LW_diag%aod_soot(:,:,:) = 0.0
      Endif
      If (LW_diag%L_aod_biomass) then
        allocate(LW_diag%aod_biomass(row_length, rows,                  &
     &                               N_AOD_WAVEL_LW))
        LW_diag%aod_biomass(:,:,:) = 0.0
      Endif
      If (LW_diag%L_aod_biogenic) then
        allocate(LW_diag%aod_biogenic(row_length, rows,                 &
     &                                N_AOD_WAVEL_LW))
        LW_diag%aod_biogenic(:,:,:) = 0.0
      Endif
      If (LW_diag%L_aod_ocff) then
        allocate(LW_diag%aod_ocff(row_length, rows,                     &
     &                            N_AOD_WAVEL_LW))
        LW_diag%aod_ocff(:, :, :) = 0.0
      Endif
      If (LW_diag%L_aod_delta) then
        allocate(LW_diag%aod_delta(row_length, rows,                    &
     &                             N_AOD_WAVEL_LW))
        LW_diag%aod_delta(:, :, :) = 0.0
      Endif

!
!     With some compilers, it seems that an allocation must
!     always be made, so a minimal allocation is made in the case
!     when the diagnostic is not selected.
!
      If (.NOT.LW_diag%L_total_cloud_cover)                             &
     &  allocate(LW_diag%total_cloud_cover(1, 1))
      If (.NOT.LW_diag%L_clear_olr)                                     &
     &  allocate(LW_diag%clear_olr(1, 1))
      If (.NOT.LW_diag%L_surface_down_flux)                             &
     &  allocate(LW_diag%surface_down_flux(1, 1))
      If (.NOT.LW_diag%L_surf_down_clr)                                 &
     &  allocate(LW_diag%surf_down_clr(1, 1))
      If (.NOT.LW_diag%L_clear_hr)                                      &
     &  allocate(LW_diag%clear_hr(1, 1, 1))
      If (.NOT.LW_diag%L_net_flux_trop)                                 &
     &  allocate(LW_diag%net_flux_trop(1, 1))
      If (.NOT.LW_diag%L_down_flux_trop)                                &
     &  allocate(LW_diag%down_flux_trop(1, 1))
!
!
      If (.NOT.LW_diag%L_cloud_absorptivity)                            &
     &  allocate(LW_diag%cloud_absorptivity(1, 1, 1))
      If (.NOT.LW_diag%L_cloud_weight_absorptivity)                     &
     &  allocate(LW_diag%cloud_weight_absorptivity(1, 1, 1))
      If (.NOT.LW_diag%L_ls_cloud_absorptivity)                         &
     &  allocate(LW_diag%ls_cloud_absorptivity(1, 1, 1))
      If (.NOT.LW_diag%L_ls_cloud_weight_absorptivity)                  &
     &  allocate(LW_diag%ls_cloud_weight_absorptivity(1, 1, 1))
      If (.NOT.LW_diag%L_cnv_cloud_absorptivity)                        &
     &  allocate(LW_diag%cnv_cloud_absorptivity(1, 1, 1))
      If (.NOT.LW_diag%L_cnv_cloud_weight_absorptivity)                 &
     &  allocate(LW_diag%cnv_cloud_weight_absorptivity(1, 1, 1))
      If (.NOT.LW_diag%L_total_cloud_on_levels)                         &
     &  allocate(LW_diag%total_cloud_on_levels(1, 1, 1))
      If (.NOT.LW_diag%L_ls_qcl_rad)                                    &
     &  allocate(LW_diag%ls_qcl_rad(1, 1, 1))
      If (.NOT.LW_diag%L_ls_qcf_rad)                                    &
     &  allocate(LW_diag%ls_qcf_rad(1, 1, 1))
      If (.NOT.LW_diag%L_cc_qcl_rad)                                    &
     &  allocate(LW_diag%cc_qcl_rad(1, 1, 1))
      If (.NOT.LW_diag%L_cc_qcf_rad)                                    &
     &  allocate(LW_diag%cc_qcf_rad(1, 1, 1))
      If (.NOT.LW_diag%L_ls_cl_rad)                                     &
     &  allocate(LW_diag%ls_cl_rad(1, 1, 1))
      If (.NOT.LW_diag%L_ls_cf_rad)                                     &
     &  allocate(LW_diag%ls_cf_rad(1, 1, 1))
      If (.NOT.LW_diag%L_cc_cl_rad)                                     &
     &  allocate(LW_diag%cc_cl_rad(1, 1, 1))
      If (.NOT.LW_diag%L_cc_cf_rad)                                     &
     &  allocate(LW_diag%cc_cf_rad(1, 1, 1))
!
      If (.NOT.LW_diag%L_aod_sulphate)                                  &
     &  allocate(LW_diag%aod_sulphate(1,1,1))
      If (.NOT.LW_diag%L_aod_dust)                                      &
     &  allocate(LW_diag%aod_dust(1,1,1))
      If (.NOT.LW_diag%L_aod_seasalt)                                   &
     &  allocate(LW_diag%aod_seasalt(1,1,1))
      If (.NOT.LW_diag%L_aod_soot)                                      &
     &  allocate(LW_diag%aod_soot(1,1,1))
      If (.NOT.LW_diag%L_aod_biomass)                                   &
     &  allocate(LW_diag%aod_biomass(1,1,1))
      If (.NOT.LW_diag%L_aod_biogenic)                                  &
     &  allocate(LW_diag%aod_biogenic(1,1,1))
      If (.NOT.LW_diag%L_aod_ocff)                                      &
     &  allocate(LW_diag%aod_ocff(1,1,1))
      If (.NOT.LW_diag%L_aod_delta)                                     &
     &  allocate(LW_diag%aod_delta(1,1,1))
!
! DEPENDS ON: timer
        If (Ltimer) Call timer ('LW Rad  ',3)

        If (L_Rad_Step) Then

! Effective solid surface surface radiative temperature
          do j = 1, rows
            do i = 1, row_length
              T_sol_rad(i,j) = 0.0
            end do
          end do
          if ( L_CTILE ) then
            do j = 1, rows
              do i = 1, row_length
                T_sol_rad(i,j) = (1.-flandg(i,j))*                      &
     &            ice_fract(i,j)*Tstar_sice(i,j)**4
              end do
            end do
            do n=1,ntiles
              do point=1,tile_pts(n)
                l = tile_index(point,n)
                i = land_index_i(l)
                j = land_index_j(l)
                T_sol_rad(i,j) = T_sol_rad(i,j) + flandg(i,j)*          &
     &                      tile_frac(l,n) * tstar_tile(l,n)**4
              enddo
            enddo
            do j = 1, rows
              do i = 1, row_length
                if(fracsolid(i,j) >  0.0)T_sol_rad(i,j) =               &
     &             (T_sol_rad(i,j)/fracsolid(i,j))**0.25
              enddo
            enddo
          elseif ( L_MOSES_II ) then
            do n=1,ntiles
              do point=1,tile_pts(n)
                l = tile_index(point,n)
                i = land_index_i(l)
                j = land_index_j(l)
                T_sol_rad(i,j) = T_sol_rad(i,j) + tile_frac(l,n) *      &
     &                                              tstar_tile(l,n)**4
              enddo
            enddo
            do l=1,land_field
              i = land_index_i(l)
              j = land_index_j(l)
              T_sol_rad(i,j) = T_sol_rad(i,j)**0.25
            enddo

          endif

        lw_points = 0
        IF (L_rad_deg) THEN
          interp_index=MOD((timestep_number/A_LW_RADSTEP),2)
!
          If ( ( interp_index == 0 ) .EQV. rad_mask(1,1) ) Then
            first_data_interp = 1
          Else
            first_data_interp = 0
          Endif
!
!         Calculations need only be performed at one point on polar
!         rows, so special treatment is required in these cases.
          first_row = 1
          last_row  = rows
          If (model_domain == mt_global) Then
            If (at_extremity(PSouth)) Then
              LW_POINTS = LW_POINTS + 1
              LIST_LW_points(LW_POINTS)=1
              first_row = 2
            Endif
            If (at_extremity(PNorth)) Then
              last_row = rows-1
            Endif
          Endif
!
          Do j=first_row, last_row
            Do i=1, row_length
              If ( rad_mask(i,j) .EQV. ( interp_index == 0 ) ) Then
                lw_points                 = lw_points + 1
                list_lw_points(lw_points) = i+(j-1)*row_length
              Endif
            Enddo
          Enddo
!
!         Include one point on the northern polar row of a global
!         domain.
          If (model_domain == mt_global.and.                            &
     &                                 at_extremity(PNorth)) Then
              lw_points                 = lw_points + 1
              LIST_LW_points(lw_points) = (rows-1)*row_length + 1
          Endif
!
! Initialize all output fields.
!
          OLR(:,:)       = 0.0
          LWsea(:,:)     = 0.0
          LW_incs(:,:,:) = 0.0
          IF (L_extra_top_LW) Then
            top_absorption(:,:)=0.0
          Endif
!
!
        ELSE  ! SPATIAL DEGRADATION IS SWITCHED OFF
!
!         Calculations need only be performed at one point on polar
!         rows, so special treatment is required in these cases.
          first_row = 1
          last_row  = rows
          If (model_domain == mt_global) Then
            If (at_extremity(PSouth)) Then
              LW_POINTS = LW_POINTS + 1
              LIST_LW_points(LW_POINTS)=1
              first_row = 2
            Endif
            If (at_extremity(PNorth)) Then
              last_row = rows-1
            Endif
          Endif
!
          Do j=first_row, last_row
            Do i=1, row_length
              lw_points                 = lw_points + 1
              list_lw_points(lw_points) = i+(j-1)*row_length
            Enddo
          Enddo
!
!         Include one point on the northern polar row of a global
!         domain.
          If (model_domain == mt_global.and.                            &
     &                                 at_extremity(PNorth)) Then
              lw_points                 = lw_points + 1
              LIST_LW_points(lw_points) = (rows-1)*row_length + 1
          Endif
!
        ENDIF  ! If L_rad_deg block
!
!       Infer the row and column indices of the points where
!       calculations are required. This is very mildly inefficient,
!       but is cleaner than repeating code throughout the preceding
!       block.
        Do j=1, lw_points
          diag_row_list(j) = (LIST_LW_points(j) - 1)/row_length + 1
          diag_col_list(j) = LIST_LW_points(j) - row_length             &
     &      * (diag_row_list(j) - 1)
        Enddo
!
!

! Now, segment all the LW calculation points
        If (a_lw_seg_size > 0) Then
          ! We are specifying size of segments rather than number
          step = a_lw_seg_size
          a_lw_segments = Ceiling( Real(lw_points) / Real(step) )
        Else
          step = lw_points/a_lw_segments
        End If

        ! Alloate space for segmentation arrays
        Allocate( first_point_local( a_lw_segments ) )
        Allocate(  seg_points_local( a_lw_segments ) )

        seg_start = LIST_LW_points(1)
        first_point=1
        Do I = 1, a_lw_segments
          num_lw_points=step
          first_point=1 + (i-1) *step
          seg_points = LIST_LW_points(I*STEP)-seg_start+1
          If (i  ==  a_lw_segments) Then
            num_lw_points = lw_points - step* (a_lw_segments - 1)
            seg_points = row_length*rows-seg_start+1
          End If
          Do j = first_point, first_point+num_lw_points-1
            LIST_LW_points(j) = LIST_LW_points(j) - seg_start + 1
          End Do
          first_point_local(i) = seg_start
          seg_points_local(i) = num_lw_points
          seg_start = seg_start+seg_points
        End Do

          Do i = 1, A_LW_segments

             first_point = 1 + (i-1) * step
!
!             Set the pointer to the beginning of the current segment.
!             This is done solely to reduce the number of continuation
!             lines required by the call.
              ptr_local=first_point_local(i)
!
!             Set the actual size of arrays in the radiation code:
!             for some architectures (such as that of Cray vector
!             machines) on odd size is preferred to avoid memory
!             bank conflicts.
              nd_rad_pts=2*(seg_points_local(i)/2)+1
!
!             Set the number of layers seen in the radiation code.
!             This may optionally be 1 greater than the number used
!             in the rest of the model to avoid spurious effects
!             resulting from the upper boundary (principally in
!             stratospheric configurations).
              if (L_EXTRA_TOP_LW) then
                n_rad_layers=model_levels+1
              else
                n_rad_layers=model_levels
              endif

!           Set the first point of the dust arrays to be used.
            IF (L_USE_DUST) THEN
              FIRST_POINT_DUST=FIRST_POINT_LOCAL(i)
            ELSE
              FIRST_POINT_DUST=1
            ENDIF

!           Set the first point of the biogenic array
            IF (L_USE_BIOGENIC) THEN
              FIRST_POINT_BIOGENIC=FIRST_POINT_LOCAL(i)
            ELSE
              FIRST_POINT_BIOGENIC=1
            ENDIF
!           Set the first points of the arrays of sulphates to be used.
!           A separate assignment is necessary since
!           not be of the full size unless the sulphur cycle is on.
            IF (L_USE_SULPC_DIRECT .OR. L_USE_SULPC_INDIRECT_LW) THEN
               first_point_sulpc=first_point_local(i)
            ELSE
               first_point_sulpc=1
            ENDIF
!
            If (L_use_soot_direct) then
              first_point_soot=first_point_local(i)
            Else
              first_point_soot=1
            Endif
!
            If (L_use_bmass_direct .or. L_use_bmass_indirect) then
              first_point_biomass=first_point_local(i)
            Else
              first_point_biomass=1
            Endif
!
            If (L_use_ocff_direct .or. L_use_ocff_indirect) then
              first_point_ocff=first_point_local(i)
            Else
              first_point_ocff=1
            Endif
!
            If (L_use_seasalt_indirect .OR. L_use_seasalt_direct) Then
              first_point_seasalt=first_point_local(i)
            Else
              first_point_seasalt=1
            Endif
!
            If (n_arcl_species > 0) Then
              first_point_arcl = first_point_local(i)
            Else
              first_point_arcl = 1
            Endif
!


! DEPENDS ON: r2_lwrad
            Call r2_lwrad(error_code,                                   &
! Input data
     &        q_n(ptr_local,1,1), CO2_MMR, ozone(ptr_local,1,1),        &
     &        CO2_DIM1, CO2_DIM2, CO2_3D(ptr_local,1,1), L_CO2_3D,      &
     &        L_use_stochem_CH4, CH4_stochem(ptr_local,1,1),            &
! chemical greenhouse gas fields
     &        ngrgas, grgas_field(ptr_local,1,1,1),                     &
     &        n2o_mix_ratio, ch4_mix_ratio, cfc11_mix_ratio,            &
     &        cfc12_mix_ratio, C113MMR,HCFC22MMR,HFC125MMR,HFC134AMMR,  &
     &        T_n(ptr_local,1,1),T_surf(ptr_local,1),                   &
     &        T_sol_rad(ptr_local,1),Tstar_sea(ptr_local,1),L_CTILE,    &
     &        p_star(ptr_local,1),p_layer_boundaries(ptr_local,1,0),    &
     &        p_layer_centres(ptr_local,1,0),                           &
     &        height_theta(ptr_local,1,0), height_rho(ptr_local,1,1),   &
! Options for treating clouds
     &        .true., global_cloud_top, L_INHOM_CLOUD, INHOM_CLOUD_LW,  &
     &        DP_CORR_STRAT, DP_CORR_CONV,                              &
     &        sin_true_latitude(ptr_local,1),                           &
! Stratiform Cloud Fields
     &        L_cloud_water_partition, L_PC2,                           &
     &        area_cloud_fraction(ptr_local,1,1),                       &
     &        cf_n(ptr_local,1,1),                                      &
     &        qcl_n(ptr_local,1,1),qcf_n(ptr_local,1,1),                &

! Convective Cloud Fields
     &        cca(ptr_local,1,1), cclwp(ptr_local,1),                   &
     &        ccw(ptr_local,1,1), lcbase(ptr_local,1),                  &
     &        ccb(ptr_local,1),   cct(ptr_local,1),                     &

!                       Surface Fields
! Only want the 0.5 threshold LAND mask and fractional land:
     &        land0p5(ptr_local,1),flandg(ptr_local,1),                 &
     &        ice_fract(ptr_local,1),snow_depth(ptr_local,1),           &
! Aerosol Fields
     &        l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero,   &
     &        zh, aero_bl_levels,                                       &
     &        L_use_clearrh,L_USE_DUST, DUST_DIM1, DUST_DIM2,           &
     &        DUST_1(FIRST_POINT_DUST,1),DUST_2(FIRST_POINT_DUST,1),    &
     &        DUST_3(FIRST_POINT_DUST,1),DUST_4(FIRST_POINT_DUST,1),    &
     &        DUST_5(FIRST_POINT_DUST,1),DUST_6(FIRST_POINT_DUST,1),    &
     &        L_USE_BIOGENIC, biogenic_dim1, biogenic_dim2,             &
     &        local_biogenic(first_point_biogenic, 1),                  &
     &        L_use_sulpc_direct, L_use_sulpc_indirect_LW, sulp_dim1,   &
     &        sulp_dim2, accum_sulphate(first_point_sulpc, 1),          &
     &        aitken_sulphate(first_point_sulpc, 1),                    &
     &        diss_sulphate(first_point_sulpc, 1),                      &
     &        L_volcts, volcmass(ptr_local, 1),                         &
     &        sea_salt_film(first_point_seasalt,1,1),                   &
     &        sea_salt_jet(first_point_seasalt,1,1),                    &
     &        L_use_seasalt_indirect, L_use_seasalt_direct,             &
     &        salt_dim1*salt_dim2, salt_dim3, L_use_soot_direct,        &
     &        soot_dim1, soot_dim2, fresh_soot(first_point_soot, 1),    &
     &        aged_soot(first_point_soot, 1), L_use_bmass_direct,       &
     &        bmass_dim1,bmass_dim2,fresh_bmass(first_point_biomass,1), &
     &        aged_bmass(first_point_biomass, 1),                       &
     &        cloud_bmass(first_point_biomass, 1), L_use_bmass_indirect,&
     &        L_use_ocff_direct, ocff_dim1, ocff_dim2,                  &
     & fresh_ocff(first_point_ocff, 1), aged_ocff(first_point_ocff, 1), &
     &        cloud_ocff(first_point_ocff, 1), L_use_ocff_indirect,     &
     &        L_USE_ARCL, arcl_dim1, arcl_dim2, n_arcl_species,         &
     &        n_arcl_compnts, i_arcl_compnts,                           &
     &        local_arcl(first_point_arcl,1,1),                         &
     &        aerosol(first_point_local(i),1,1), L_MURK_RAD,            &
     &        Ntot_land, Ntot_sea,                                      &
! Level of tropopause
     &        trindx(ptr_local, 1)                                      &
! Size and control variables
!
!                       Spectral data
!     ------------------------------------------------------------------
!     ARGUMENT LIST OF LW SPECTRAL DATA.
!     (NOTE: LWSPDC3A, LWSPCM3A AND LWSARG3A MUST BE CONSISTENT)
!
     &   , NPD_TYPE_LW, NPD_BAND_LW, NPD_EXCLUDE_LW                     &
     &   , NPD_SPECIES_LW, NPD_ESFT_TERM_LW, NPD_SCALE_FNC_LW           &
     &   , NPD_SCALE_VARIABLE_LW, NPD_THERMAL_COEFF_LW                  &
     &   , NPD_SURFACE_LW, NPD_ALBEDO_PARM_LW, NPD_CONTINUUM_LW         &
     &   , NPD_DROP_TYPE_LW, NPD_ICE_TYPE_LW, NPD_CLOUD_PARAMETER_LW    &
     &   , NPD_AEROSOL_SPECIES_LW, NPD_HUMIDITIES_LW, NPD_AOD_WAVEL_LW  &
     &   , L_PRESENT_LW, N_BAND_LW, WAVE_LENGTH_SHORT_LW                &
     &   , WAVE_LENGTH_LONG_LW, N_BAND_EXCLUDE_LW, INDEX_EXCLUDE_LW     &
     &   , SOLAR_FLUX_BAND_LW, RAYLEIGH_COEFFICIENT_LW, N_ABSORB_LW     &
     &   , N_BAND_ABSORB_LW, INDEX_ABSORB_LW, TYPE_ABSORB_LW            &
     &   , I_BAND_ESFT_LW, I_SCALE_ESFT_LW, I_SCALE_FNC_LW, K_ESFT_LW   &
     &   , W_ESFT_LW, SCALE_VECTOR_LW, P_REFERENCE_LW, T_REFERENCE_LW   &
     &   , N_DEG_FIT_LW, THERMAL_COEFFICIENT_LW, T_REF_PLANCK_LW        &
     &   , I_SPEC_SURFACE_LW, N_DIR_ALBEDO_FIT_LW, L_SURFACE_LW         &
     &   , SURFACE_ALBEDO_LW, DIRECT_ALBEDO_PARM_LW                     &
     &   , EMISSIVITY_GROUND_LW, N_BAND_CONTINUUM_LW                    &
     &   , INDEX_CONTINUUM_LW, INDEX_WATER_LW, I_SCALE_FNC_CONT_LW      &
     &   , K_CONTINUUM_LW, SCALE_CONTINUUM_LW, P_REF_CONTINUUM_LW       &
     &   , T_REF_CONTINUUM_LW, I_DROP_PARAMETRIZATION_LW                &
     &   , L_DROP_TYPE_LW, DROP_PARAMETER_LIST_LW                       &
     &   , DROP_PARM_MIN_DIM_LW, DROP_PARM_MAX_DIM_LW                   &
     &   , N_AEROSOL_LW, TYPE_AEROSOL_LW, I_AEROSOL_PARAMETRIZATION_LW  &
     &   , NHUMIDITY_LW, HUMIDITIES_LW, L_AEROSOL_SPECIES_LW            &
     &   , AEROSOL_ABSORPTION_LW, AEROSOL_SCATTERING_LW                 &
     &   , AEROSOL_ASYMMETRY_LW, I_ICE_PARAMETRIZATION_LW               &
     &   , L_ICE_TYPE_LW, ICE_PARAMETER_LIST_LW                         &
     &   , ICE_PARM_MIN_DIM_LW, ICE_PARM_MAX_DIM_LW                     &
     &   , L_DOPPLER_PRESENT_LW, DOPPLER_CORRECTION_LW                  &
     &   , N_AOD_WAVEL_LW, AOD_ABSORPTION_LW, AOD_SCATTERING_LW         &
     &   , I_AOD_TYPE_LW                                                &
!
!     ------------------------------------------------------------------
!
!                       Algorithmic options
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
!
!     ------------------------------------------------------------------
!
     & ,      timestep,l_mod_k_flux,L_Wenyi,LIST_LW_points(first_point),&

! All diagnostics
     &        LW_diag,                                                  &
     &        diag_row_list(first_point), diag_col_list(first_point),   &

! Physical Dimensions
     &        seg_points_local(i),model_levels,n_rad_layers,            &
     &        cloud_levels,wet_model_levels,ozone_levels,row_length,    &
     &        rows,row_length*rows, nd_rad_pts, n_rad_layers, 1,        &
     &        n_cca_levels,                                             &

! Output data
     &        OLR(ptr_local, 1), top_absorption(ptr_local, 1),          &
     &        LWsea(ptr_local, 1), LW_incs(ptr_local,1,0),              &
     &        rho_r2_nh(ptr_local,1),                                   &
     & r_rho_levels_nh(ptr_local,1), r_theta_levels_nh(ptr_local,0),    &
     & q_n(ptr_local,1,1), qcl_n(ptr_local,1,1),qcf_n(ptr_local,1,1),   &
     & qcf2_n(ptr_local,1,1),qrain_n(ptr_local,1,1),                    &
     & qgraup_n(ptr_local,1,1),                                         &
     & l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio            &
     & )

          End Do ! end loop over long-wave segments

        ! Deallocate the segmentation arrays
        Deallocate( first_point_local )
        Deallocate(  seg_points_local )


!         Radiative fluxes may not have been calculated at all
!         points. We now fill in as required.
!
!         At the North Pole in a global domain calculations
!         are performed only at the first point, so the rest of
!         the row must be filled in.
          L_complete_North = ( model_domain == mt_global ) .AND.        &
     &                       ( at_extremity(PNorth) )

!         At the South Pole in a global domain calculations
!         are performed only at the first point, so the rest of
!         the row must be filled in.
          L_complete_South = ( model_domain == mt_global ) .AND.        &
     &                       ( at_extremity(PSouth) )
!
!         When spatial degradation is performed fields must be
!         filled in at alternate points.
          L_complete_deg = ( L_rad_deg )
!
!         Set addressing limits for spatial degradation.
          If ( L_complete_deg ) Then
            first_row=1
            last_row=rows
            If (model_domain  ==  mt_global) Then
              If (at_extremity(PNorth)) Then
                last_row=rows-1
              End If
              If (at_extremity(PSouth)) Then
                first_row=2
              End If
            End If
          Endif
!
!
!         Primary Fields:
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp,                         &
     &      model_levels+1,                                             &
     &      LW_incs                                                     &
     &      )
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      OLR                                                         &
     &      )
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      LWsea                                                       &
     &      )
          If ( L_extra_top_LW ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, 1,                    &
     &        top_absorption                                            &
     &        )
          Endif
!
!         LW Diagnostics:
!
          If ( LW_diag%L_total_cloud_cover ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, 1,                    &
     &        LW_diag%total_cloud_cover                                 &
     &        )
          Endif
          If ( LW_diag%L_clear_olr ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, 1,                    &
     &        LW_diag%clear_olr                                         &
     &        )
          Endif
          If ( LW_diag%L_surface_down_flux ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, 1,                    &
     &        LW_diag%surface_down_flux                                 &
     &        )
          Endif
          If ( LW_diag%L_surf_down_clr ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, 1,                    &
     &        LW_diag%surf_down_clr                                     &
     &        )
          Endif
          If ( LW_diag%L_clear_hr ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp,                       &
     &        model_levels,                                             &
     &        LW_diag%clear_hr                                          &
     &        )
          Endif
          If ( LW_diag%L_net_flux_trop ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, 1,                    &
     &        LW_diag%net_flux_trop                                     &
     &        )
          Endif
          If ( LW_diag%L_down_flux_trop ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, 1,                    &
     &        LW_diag%down_flux_trop                                    &
     &        )
          Endif
          If ( LW_diag%L_total_cloud_on_levels ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, cloud_levels,         &
     &        LW_diag%total_cloud_on_levels                             &
     &        )
          Endif
! Grid-box mean cloud diagnostics as seen by radiation:
          If ( LW_diag%L_ls_qcl_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag%ls_qcl_rad                                        &
     &        )
          Endif
          If ( LW_diag%L_ls_qcf_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag%ls_qcf_rad                                        &
     &        )
          Endif
          If ( LW_diag%L_cc_qcl_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag%cc_qcl_rad                                        &
     &        )
          Endif
          If ( LW_diag%L_cc_qcf_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag%cc_qcf_rad                                        &
     &        )
          Endif
          If ( LW_diag%L_ls_cl_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag%ls_cl_rad                                         &
     &        )
          Endif
          If ( LW_diag%L_ls_cf_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag%ls_cf_rad                                         &
     &        )
          Endif
          If ( LW_diag%L_cc_cl_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag%cc_cl_rad                                         &
     &        )
          Endif
          If ( LW_diag%L_cc_cf_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag%cc_cf_rad                                         &
     &        )
          Endif
!   Absorptivity diagnostics:
          If ( LW_diag%L_cloud_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, cloud_levels,         &
     &        LW_diag%cloud_absorptivity                                &
     &        )
          Endif
          If ( LW_diag%L_cloud_weight_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, cloud_levels,         &
     &        LW_diag%cloud_weight_absorptivity                         &
     &        )
          Endif
          If ( LW_diag%L_ls_cloud_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, cloud_levels,         &
     &        LW_diag%ls_cloud_absorptivity                             &
     &        )
          Endif
          If ( LW_diag%L_ls_cloud_weight_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, cloud_levels,         &
     &        LW_diag%ls_cloud_weight_absorptivity                      &
     &        )
          Endif
          If ( LW_diag%L_cnv_cloud_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, cloud_levels,         &
     &        LW_diag%cnv_cloud_absorptivity                            &
     &        )
          Endif
          If ( LW_diag%L_cnv_cloud_weight_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, cloud_levels,         &
     &        LW_diag%cnv_cloud_weight_absorptivity                     &
     &        )
          Endif
! Aerosol optical depth diagnostics
          If ( LW_diag%L_aod_sulphate) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, N_AOD_WAVEL_LW,       &
     &        LW_diag%aod_sulphate                                      &
     &        )
          Endif
          If ( LW_diag%L_aod_dust) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, N_AOD_WAVEL_LW,       &
     &        LW_diag%aod_dust                                          &
     &        )
          Endif
          If ( LW_diag%L_aod_seasalt) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, N_AOD_WAVEL_LW,       &
     &        LW_diag%aod_seasalt                                       &
     &        )
          Endif
          If ( LW_diag%L_aod_soot) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, N_AOD_WAVEL_LW,       &
     &        LW_diag%aod_soot                                          &
     &        )
          Endif
          If ( LW_diag%L_aod_biomass) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, N_AOD_WAVEL_LW,       &
     &        LW_diag%aod_biomass                                       &
     &        )
          Endif
          If ( LW_diag%L_aod_biogenic) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, N_AOD_WAVEL_LW,       &
     &        LW_diag%aod_biogenic                                      &
     &        )
          Endif
          If ( LW_diag%L_aod_ocff) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, N_AOD_WAVEL_LW,       &
     &        LW_diag%aod_ocff                                          &
     &        )
          Endif
          If ( LW_diag%L_aod_delta) THEN
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, N_AOD_WAVEL_LW,       &
     &        LW_diag%aod_delta                                         &
     &        )
          Endif
!
!
!
          If ( L_CTILE ) then
            Do j = 1, rows
              Do i= 1, row_length
                surflw(i,j) = LW_incs(i,j,0)                            &
     &                      + (1.-flandg(i,j))*LWsea(i,j)
                if(flandg(i,j) == 1.0)LWsea(i,j)=rmdi
              End Do
            End Do
          Else
            Do j = 1, rows
              Do i= 1, row_length
                surflw(i,j) = LW_incs(i,j,0) + LWsea(i,j)
                if(land_sea_mask(i,j))LWsea(i,j)=rmdi
              End Do
            End Do
          Endif

!
!         The logic in the following IF-block is rather subtle
!         and is best clarified. Note first that:
!           * L_CTILE and L_MOSES_II are independent flags
!           * T_sol_rad is the mean over all parts of the grid-box,
!               excluding open sea
!           * dOLR is a grid-box mean quantity
!           * ice_fract is the fraction of sea-ice in the portion
!             of the grid-box which is not land
!
!         Now consider the possible cases:
!           * If there is no land in the grid-box, ice_fract will
!             be equal to the total fraction of sea-ice in the
!             grid-boxes, so the change to dOLR in the IF-block under
!             L_MOSES_II will indeed be calculated using grid-box
!             mean quantities.
!           * If there is any land in the grid-box, dOLR should
!             be adjusted using fracsolid, the total fraction of
!             solid surface in the grid-box, but since this loop
!             is defined only over grid-boxes containing land, it 
!             catches cases not addressed by the foregoing loop.
!             T_sol_rad is calculated earlier under logic involving
!             L_CTILE. If this is true, T_sol_rad is clearly the
!             correct mean over the solid fraction of the grid-box;
!             but if it's false the mean is defined only over the
!             land fraction of the grid-box, which appears at first
!             sight to be incorrect. However, if L_CTILE is false,
!             a grid-box must be either all land or all sea, so
!             T_sol_rad provides the correct mean in this case as
!             well.  
!

          If ( L_MOSES_II ) then
! Store downward surface LW and dOLR on radiation timestep
! dOLR is TOA outward LW - surface upward LW for land and sea-ice
            Do j = 1, rows
              Do i= 1, row_length
                dOLR(i,j) = OLR(i,j)
                If(.not.land_sea_mask(i,j).and.ice_fract(i,j) >  0.)then
                  dOLR(i,j) = dOLR(i,j) -                               &
     &                        ice_fract(i,j)*sbcon*tstar_sice(i,j)**4
                End If
              End Do
            End Do
            Do l=1,land_field
              i = land_index_i(l)
              j = land_index_j(l)
              dOLR(i,j) = dOLR(i,j) -                                   &
     &          fracsolid(i,j)*sbcon*T_sol_rad(i,j)**4
            End Do
            Do j = 1, rows
              Do i= 1, row_length
                LW_down(i,j) = LW_Diag%surface_down_flux(i,j)
                dOLR_rts(i,j) = dOLR(i,j)
              End Do
            End Do
          End If

        End If ! end conditional on being a radiation timestep
!
! Is the PC2 cloud scheme being used?
!
        If (L_pc2) then
!
! Reset _latest values to _n values
!
          Do k = 1,model_levels
            Do j = 1, rows
              Do i = 1, row_length
                T_latest(i,j,k)   = T_n(i,j,k)
              End Do
            End Do
          End Do
          Do k = 1,wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_latest(i,j,k)   = q_n(i,j,k)
                qcl_latest(i,j,k) = qcl_n(i,j,k)
                cf_latest(i,j,k)  = cf_n(i,j,k)
                cfl_latest(i,j,k) = cfl_n(i,j,k)
              End Do
            End Do
          End Do
!
! ----------------------------------------------------------------------
! Homogeneous forcing. Note the temperature increment from longwave
! is added in this routine
! ----------------------------------------------------------------------
! DEPENDS ON: pc2_homog_plus_turb
          CALL PC2_HOMOG_PLUS_TURB(p_layer_centres(1,1,1),              &
     &      wet_model_levels,                                           &
     &      row_length, rows, timestep, T_latest, cf_latest, cfl_latest,&
     &      cff_latest, q_latest, qcl_latest, LW_incs(1,1,1),           &
     &      zeros, zeros, zeros, 0.0, 0.0, l_mixing_ratio)

! Add increments from the homogeneous forcing to the increment variables
!
          Do k = 1,model_levels
            Do j = 1, rows
              Do i = 1, row_length
                T_inc(i,j,k) = T_inc(i,j,k) + T_latest(i,j,k)-T_n(i,j,k)
              End Do
            End Do
          End Do
          Do k = 1,wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_inc(i,j,k) = q_inc(i,j,k) + q_latest(i,j,k)-q_n(i,j,k)
                qcl_inc(i,j,k) = qcl_inc(i,j,k)                         &
     &                           + qcl_latest(i,j,k)-qcl_n(i,j,k)
                cf_inc(i,j,k) = cf_inc(i,j,k)                           &
     &                           + cf_latest(i,j,k)-cf_n(i,j,k)
                cfl_inc(i,j,k) = cfl_inc(i,j,k)                         &
     &                           + cfl_latest(i,j,k)-cfl_n(i,j,k)
              End Do
            End Do
          End Do
!
        Else  ! L_pc2
!
          Do k = 1,model_levels
            Do j = 1, rows
              Do i = 1, row_length
                T_inc(i,j,k) = T_inc(i,j,k) + LW_incs(i,j,k)
                T_latest(i,j,k) = T_n(i,j,k) + LW_incs(i,j,k)

              End Do
            End Do
          End Do
!
        End If  ! L_pc2
!
! Get T_incr for output as STASH diagnostic
      If ( L_T_incr_lw ) Then  ! STASHflag set
        Allocate ( T_incr_diagnostic(row_length,rows,model_levels) )
        Do k=1,model_levels
         Do j=1,rows
          Do i=1,row_length
            T_incr_diagnostic(i,j,k) =  LW_incs(i,j,k)
          Enddo ! i
         Enddo ! j
        Enddo ! k
      Else
        Allocate ( T_incr_diagnostic(1,1,1) )
      Endif                    ! on STASHflag

! Set up total net down surface radiation flux in SURF_RADFLUX

        If (L_MOSES_II) then
          Do j = 1, rows
            Do i= 1, row_length
              If (ice_fract(i,j) >  0.)                                 &
     &          surf_radflux(i,j) = surf_radflux(i,j) +                 &
     &                              ice_fract(i,j)*lw_down(i,j)
            End Do
          End Do
        Else
          Do j = 1, rows
            Do i= 1, row_length
              surf_radflux(i,j) = LW_incs(i,j,0) + surf_radflux(i,j)
            End Do
          End Do
        End If

! Set up radiative heating rates for 6A boundary layer code

        recip_timestep = 1./ timestep
        Do k = 1, bl_levels
          Do j = 1, rows
            Do i= 1, row_length
              rad_hr(i,j,k,1) = LW_incs(i,j,k) * recip_timestep
              rad_hr(i,j,k,2) = SW_incs(i,j,k)                          &
     &                 * cos_zenith_angle(i,j) * recip_timestep
            End Do
          End Do
        End Do

! Copy dOLR from last LW timestep for MOSES II
        If (L_MOSES_II) then
          Do j = 1, rows
            Do i= 1, row_length
              dOLR(i,j) = dOLR_rts(i,j)
            End Do
          End Do
        End If

! DEPENDS ON: timer
        If (Ltimer) Call timer ('LW Rad  ',4)

! ----------------------------------------------------------------------
! Section RAD.2.1 Isccp diagnostics
!-----------------------------------------------------------------------

! Allocate space for the diagnostic arrays required, and zero the
! elements explicitly: we will overwrite at lit points.

      If (LW_diag%L_isccp_weights) then
       allocate(LW_diag%isccp_weights(row_length,rows))
       LW_diag%isccp_weights(:,:) = 0.0
      End if
      If (LW_diag%L_isccp_cf) then
       allocate(LW_diag%isccp_cf(row_length,rows,7))
       LW_diag%isccp_cf(:,:,:) = 0.0
      End if
      If (LW_diag%L_isccp_cf_tau_0_to_p3) then
       allocate(LW_diag%isccp_cf_tau_0_to_p3(row_length,rows,7))
       LW_diag%isccp_cf_tau_0_to_p3(:,:,:) = 0.0
      End if
      If (LW_diag%L_isccp_cf_tau_p3_to_1p3) then
       allocate(LW_diag%isccp_cf_tau_p3_to_1p3(row_length,rows,7))
       LW_diag%isccp_cf_tau_p3_to_1p3(:,:,:) = 0.0
      End if
      If (LW_diag%L_isccp_cf_tau_1p3_to_3p6) then
       allocate(LW_diag%isccp_cf_tau_1p3_to_3p6(row_length,rows,7))
       LW_diag%isccp_cf_tau_1p3_to_3p6(:,:,:) = 0.0
      End if
      If (LW_diag%L_isccp_cf_tau_3p6_to_9p4) then
       allocate(LW_diag%isccp_cf_tau_3p6_to_9p4(row_length,rows,7))
       LW_diag%isccp_cf_tau_3p6_to_9p4(:,:,:) = 0.0
      End if
      If (LW_diag%L_isccp_cf_tau_9p4_to_23) then
       allocate(LW_diag%isccp_cf_tau_9p4_to_23(row_length,rows,7))
       LW_diag%isccp_cf_tau_9p4_to_23(:,:,:) = 0.0
      End if
      If (LW_diag%L_isccp_cf_tau_23_to_60) then
       allocate(LW_diag%isccp_cf_tau_23_to_60(row_length,rows,7))
       LW_diag%isccp_cf_tau_23_to_60(:,:,:) = 0.0
      End if
      If (LW_diag%L_isccp_cf_tau_ge_60) then
       allocate(LW_diag%isccp_cf_tau_ge_60(row_length,rows,7))
       LW_diag%isccp_cf_tau_ge_60(:,:,:) = 0.0
      End if
      If (LW_diag%L_meanalbedocld) then
       allocate(LW_diag%meanalbedocld(row_length,rows))
       LW_diag%meanalbedocld(:,:) = 0.0
      End if
      If (LW_diag%L_meantaucld) then
       allocate(LW_diag%meantaucld(row_length,rows))
       LW_diag%meantaucld(:,:) = 0.0
      End if
      If (LW_diag%L_meanptop) then
       allocate(LW_diag%meanptop(row_length,rows))
       LW_diag%meanptop(:,:) = 0.0
      End if
      If (LW_diag%L_totalcldarea) then
       allocate(LW_diag%totalcldarea(row_length,rows))
       LW_diag%totalcldarea(:,:) = 0.0
      End if

!     With some compilers, it seems that an allocation must
!     always be made, so a minimal allocation is made in the case
!     when the diagnostic is not selected.

      If (.NOT.LW_diag%L_isccp_weights)                                 &
     &  allocate(LW_diag%isccp_weights(1, 1))
      If (.NOT.LW_diag%L_isccp_cf)                                      &
     &  allocate(LW_diag%isccp_cf(1, 1, 1))
      If (.NOT.LW_diag%L_isccp_cf_tau_0_to_p3)                          &
     &  allocate(LW_diag%isccp_cf_tau_0_to_p3(1, 1, 1))
      If (.NOT.LW_diag%L_isccp_cf_tau_p3_to_1p3)                        &
     &  allocate(LW_diag%isccp_cf_tau_p3_to_1p3(1, 1, 1))
      If (.NOT.LW_diag%L_isccp_cf_tau_1p3_to_3p6)                       &
     &  allocate(LW_diag%isccp_cf_tau_1p3_to_3p6(1, 1, 1))
      If (.NOT.LW_diag%L_isccp_cf_tau_3p6_to_9p4)                       &
     &  allocate(LW_diag%isccp_cf_tau_3p6_to_9p4(1, 1, 1))
      If (.NOT.LW_diag%L_isccp_cf_tau_9p4_to_23)                        &
     &  allocate(LW_diag%isccp_cf_tau_9p4_to_23(1, 1, 1))
      If (.NOT.LW_diag%L_isccp_cf_tau_23_to_60)                         &
     &  allocate(LW_diag%isccp_cf_tau_23_to_60(1, 1, 1))
      If (.NOT.LW_diag%L_isccp_cf_tau_ge_60)                            &
     &  allocate(LW_diag%isccp_cf_tau_ge_60(1, 1, 1))
      If (.NOT.LW_diag%L_meanalbedocld)                                 &
     &  allocate(LW_diag%meanalbedocld(1, 1))
      If (.NOT.LW_diag%L_meantaucld)                                    &
     &  allocate(LW_diag%meantaucld(1, 1))
      If (.NOT.LW_diag%L_meanptop)                                      &
     &  allocate(LW_diag%meanptop(1, 1))
      If (.NOT.LW_diag%L_totalcldarea)                                  &
     &  allocate(LW_diag%totalcldarea(1, 1))

      If ( L_Rad_Step .AND. (LW_diag%L_isccp_cf .OR.                    &
     &                       LW_diag%L_isccp_cf_tau_0_to_p3 .OR.        &
     &                       LW_diag%L_isccp_cf_tau_p3_to_1p3 .OR.      &
     &                       LW_diag%L_isccp_cf_tau_1p3_to_3p6 .OR.     &
     &                       LW_diag%L_isccp_cf_tau_3p6_to_9p4 .OR.     &
     &                       LW_diag%L_isccp_cf_tau_9p4_to_23 .OR.      &
     &                       LW_diag%L_isccp_cf_tau_23_to_60 .OR.       &
     &                       LW_diag%L_isccp_cf_tau_ge_60) .AND.        &
     &                       error_code  ==  0) Then


! This code is repeated from the earlier call to sw_rad (the
! diagnostics are calculated at lit points only)

          If ( daylight_points  >   0 ) Then

            Do i = 1, N_SW_segments

! Restore values saved from shortwave
              first_point_temp(i)=first_point_temp_sw(i)
              seg_points_temp(i)=seg_points_temp_sw(i)
              step=step_sw

              lit_points = step
              start_point = 1+(i-1)*step
              If ( i  ==  N_SW_segments ) Then
                lit_points = daylight_points - step* (N_SW_segments - 1)
              End If
              first_point = first_point_temp(i)

!             Set the actual size of arrays in the radiation code:
!             for some architectures (such as that of Cray vector
!             machines) on odd size is preferred to avoid memory
!             bank conflicts.
              nd_rad_pts=2*(lit_points/2)+1
!
!             Set the number of layers seen in the radiation code.
!             This may optionally be 1 greater than the number used
!             in the rest of the model to avoid spurious effects
!             resulting from the upper boundary (principally in
!             stratospheric configurations).
              if (L_EXTRA_TOP_SW) then
                n_rad_layers=model_levels+1
              else
                n_rad_layers=model_levels
              endif

! DEPENDS ON: isccp
              Call ISCCP(error_code, i,                                 &
!                       Mixing Ratios
     &        q_n(first_point,1,1),                                     &
!                       Pressures and Temperatures
     &        T_surf(first_point,1),p_star(first_point,1),              &
     &        p_layer_boundaries(first_point,1,0),                      &
     &        p_layer_centres(first_point,1,0),T_n(first_point,1,1),    &

!                       Stratiform Cloud Fields
     &        area_cloud_fraction(first_point,1,1),                     &
     &        cf_n(first_point,1,1),                                    &
!                       Convective Cloud Fields
     &        cca(first_point,1,1), cclwp(first_point,1),               &
     &        ccb(first_point,1), cct(first_point,1),                   &
!                       Solar Fields
     &        day_fraction(first_point,1),                              &
     &        list_daylight_points(start_point), SCS,                   &
!                       Level of tropopause
     &        trindx(first_point,1),                                    &
!                       All diagnostics and associated arrays
     &        LW_diag, SW_diag,                                         &
     &        diag_row_list_sw(start_point),                            &
     &        diag_col_list_sw(start_point),                            &
!                       Dimensions
     &        lit_points,seg_points_temp(i),model_levels,n_rad_layers,  &
     &        cloud_levels,wet_model_levels,ozone_levels,               &
     &        rows*row_length,nd_rad_pts,n_rad_layers,1,                &
     &        n_cca_levels,                                             &
!
! Variables needed to calculate layer masses
     &        rho_r2_nh(first_point,1),                                 &
     &        r_rho_levels_nh(first_point,1),                           &
     &        r_theta_levels_nh(first_point,0),                         &
     &        q_n(first_point,1,1), qcl_n(first_point,1,1),             &
     &        qcf_n(first_point,1,1), qcf2_n(first_point,1,1),          &
     &        qrain_n(first_point,1,1), qgraup_n(first_point,1,1),      &
     &        l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio     &
     &        )

        End Do


!           Radiative fluxes may not have been calculated at all
!           points. We now fill in as required.
!
!           At the North Pole in a global domain calculations
!           are performed only at the first point if lit, so
!           the rest of the row must be filled in if this point
!           is lit.
            L_complete_North = ( model_domain == mt_global ) .AND.      &
     &                         ( at_extremity(PNorth) ) .AND.           &
     &                         ( switch(1,rows) )

!           At the South Pole in a global domain calculations
!           are performed only at the first point if lit, so
!           the rest of the row must be filled in if this point
!           is lit.
            L_complete_South = ( model_domain == mt_global ) .AND.      &
     &                         ( at_extremity(PSouth) ) .AND.           &
     &                         ( switch(1,1) )
!
!           When spatial degradation is performed fields must be
!           filled in at alternate points.
            L_complete_deg = ( L_rad_deg ) .AND.                        &
     &                       ( tot_daylight_points > 0 )
!
!           Set addressing limits for spatial degradation.
            If ( L_complete_deg ) Then
              If (model_domain  ==  mt_global) Then
                If (at_extremity(PNorth)) Then
                  first_row=1
                  last_row=rows-1
                Else If (at_extremity(PSouth)) Then
                  first_row=2
                  last_row=rows
                Else
                  first_row=1
                  last_row=rows
                End If
              Else
                first_row=1
                last_row=rows
              End If
            Endif
!
            first_data_interp = first_data_interp_sw
!
!           Call appropriate subroutines to fill in missing data
!           as required.
!
            If ( L_complete_North .OR. L_complete_South .OR.            &
     &           L_complete_deg ) Then

!   Isccp diagnostics

              If ( LW_diag%L_isccp_weights ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 1,                &
     &            LW_diag%isccp_weights                                 &
     &            )
              Endif

              If ( LW_diag%L_isccp_cf ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 7,                &
     &            LW_diag%isccp_cf                                      &
     &            )
              Endif

              If ( LW_diag%L_isccp_cf_tau_0_to_p3 ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 7,                &
     &            LW_diag%isccp_cf_tau_0_to_p3                          &
     &            )
              Endif

              If ( LW_diag%L_isccp_cf_tau_p3_to_1p3 ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 7,                &
     &            LW_diag%isccp_cf_tau_p3_to_1p3                        &
     &            )
              Endif

              If ( LW_diag%L_isccp_cf_tau_1p3_to_3p6 ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 7,                &
     &            LW_diag%isccp_cf_tau_1p3_to_3p6                       &
     &            )
              Endif

              If ( LW_diag%L_isccp_cf_tau_3p6_to_9p4 ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 7,                &
     &            LW_diag%isccp_cf_tau_3p6_to_9p4                       &
     &            )
              Endif

              If ( LW_diag%L_isccp_cf_tau_9p4_to_23 ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 7,                &
     &            LW_diag%isccp_cf_tau_9p4_to_23                        &
     &            )
              Endif

              If ( LW_diag%L_isccp_cf_tau_23_to_60 ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 7,                &
     &            LW_diag%isccp_cf_tau_23_to_60                         &
     &            )
              Endif

              If ( LW_diag%L_isccp_cf_tau_ge_60 ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 7,                &
     &            LW_diag%isccp_cf_tau_ge_60                            &
     &            )
              Endif

              If ( LW_diag%L_meanalbedocld ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 1,                &
     &            LW_diag%meanalbedocld                                 &
     &            )
              Endif

              If ( LW_diag%L_meantaucld ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 1,                &
     &            LW_diag%meantaucld                                    &
     &            )
              Endif

              If ( LW_diag%L_meanptop ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 1,                &
     &            LW_diag%meanptop                                      &
     &            )
              Endif

              If ( LW_diag%L_totalcldarea ) THEN
! DEPENDS ON: rad3d_inp
                Call rad3d_inp(                                         &
     &            L_complete_North, L_complete_South, L_complete_deg,   &
     &            row_length, rows, off_x, off_y, first_row, last_row,  &
     &            first_data_interp, ES_space_interp, 1,                &
     &            LW_diag%totalcldarea                                  &
     &            )
              Endif

            Endif

          Endif  ! On daylit points

        Endif  ! If a radiation timestep and want ISCCP diagnostics

        ! Deallocate the SW segmentation arrays
        If ( Allocated( first_point_temp ) )  Then
          Deallocate( first_point_temp )
        End If

        If ( Allocated( first_point_temp_sw ) )  Then
          Deallocate( first_point_temp_sw )
        End If

        If ( Allocated( seg_points_temp ) )  Then
          Deallocate( seg_points_temp )
        End If

        If ( Allocated( seg_points_temp_sw ) )  Then
          Deallocate( seg_points_temp_sw )
        End If

! ----------------------------------------------------------------------
! Section RAD.2.2 Long Wave Radiation Energy correction code
!-----------------------------------------------------------------------

! DEPENDS ON: timer
        If (Ltimer) Call timer ('Eng Corr',3)

        If (L_Rad_Step .and. L_emcorr) Then

! Sum long wave fluxes into the atmosphere and
! add into the net diabatic fluxes into the
! atmosphere for use in the energy correction
! procedure

          If (L_EXTRA_TOP_LW) Then
!           The energy absorbed above the top of the model in
!           the radiation scheme does not contribute to the
!           energy absorbed, but the diagnostics are calculated
!           at the top of the atmosphere, so the net atmospheric
!           flux must be adjusted.
            Do j = 1, rows
              Do i = 1, row_length
                net_atm_flux(i,j) = -OLR(i,j) - surflw(i,j)             &
     &                              -top_absorption(i,j)
              End Do
            End Do
          Else
            Do j = 1, rows
              Do i = 1, row_length
                net_atm_flux(i,j) = -OLR(i,j) - surflw(i,j)
              End Do
            End Do
          End if

! DEPENDS ON: flux_diag
          Call flux_diag(net_atm_flux, FV_cos_theta_latitude,           &
     &                   row_length, rows ,off_x,off_y,1.0,             &
     &                    sum_eng_fluxes,radiation_timestep)

        End If

! DEPENDS ON: timer
        If (Ltimer) Call timer ('Eng Corr',4)

! ----------------------------------------------------------------------
! Section RAD.2.3 Long Wave Radiation diagnostics
!-----------------------------------------------------------------------

! Check that lw diagnostics requested this timestep
        If (error_code  ==  0 .and. sf(0,2)) Then
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',3)

! DEPENDS ON: diagnostics_lw
          Call diagnostics_lw(                                          &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels, ozone_levels             &
     &,                      cloud_levels                               &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      timestep                                   &
     &,                      T_n, T_inc                                 &
     &,                      q_n, qcl_n, cf_n, cfl_n                    &
     &,                      T_latest, q_latest, qcl_latest             &
     &,                      cf_latest, cfl_latest                      &
     &,                      surflw, OLR, LW_diag%clear_olr             &
     &,                      LW_diag%total_cloud_cover                  &
     &,                      LW_diag%surface_down_flux                  &
     &,                      LW_diag%surf_down_clr                      &
     &,                      T_incr_diagnostic                          &
     &,                      LW_diag%clear_hr                           &
     &,                      LW_diag%net_flux_trop                      &
     &,                      LW_diag%down_flux_trop                     &
     &,                      LW_diag%total_cloud_on_levels              &
     &,                      LW_diag%cloud_absorptivity                 &
     &,                      LW_diag%cloud_weight_absorptivity          &
     &,                      LW_diag%ls_cloud_absorptivity              &
     &,                      LW_diag%ls_cloud_weight_absorptivity       &
     &,                      LW_diag%cnv_cloud_absorptivity             &
     &,                      LW_diag%cnv_cloud_weight_absorptivity      &
     &,                      LW_diag%isccp_weights                      &
     &,                      LW_diag%isccp_cf                           &
     &,                      LW_diag%isccp_cf_tau_0_to_p3               &
     &,                      LW_diag%isccp_cf_tau_p3_to_1p3             &
     &,                      LW_diag%isccp_cf_tau_1p3_to_3p6            &
     &,                      LW_diag%isccp_cf_tau_3p6_to_9p4            &
     &,                      LW_diag%isccp_cf_tau_9p4_to_23             &
     &,                      LW_diag%isccp_cf_tau_23_to_60              &
     &,                      LW_diag%isccp_cf_tau_ge_60                 &
     &,                      LW_diag%meanalbedocld                      &
     &,                      LW_diag%meantaucld                         &
     &,                      LW_diag%meanptop                           &
     &,                      LW_diag%totalcldarea                       &
     &,                      LW_diag%ls_qcl_rad                         &
     &,                      LW_diag%ls_qcf_rad                         &
     &,                      LW_diag%cc_qcl_rad                         &
     &,                      LW_diag%cc_qcf_rad                         &
     &,                      LW_diag%ls_cl_rad                          &
     &,                      LW_diag%ls_cf_rad                          &
     &,                      LW_diag%cc_cl_rad                          &
     &,                      LW_diag%cc_cf_rad                          &
     &,                      ozone                                      &
     &,                      O3_trop_level                              &
     &,                      O3_trop_height                             &
     &,                      T_trop_level                               &
     &,                      T_trop_height                              &
     &,                      LW_incs, LWsea                             &
     &,                      N_AOD_WAVEL_LW                             &
     &,                      LW_diag%aod_sulphate                       &
     &,                      LW_diag%aod_dust                           &
     &,                      LW_diag%aod_seasalt                        &
     &,                      LW_diag%aod_soot                           &
     &,                      LW_diag%aod_biomass                        &
     &,                      LW_diag%aod_biogenic                       &
     &,                      LW_diag%aod_ocff                           &
     &,                      LW_diag%aod_delta                          &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork2                                                       &
     & )

! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',4)
        Endif ! on error_code .and. sf(0,2)


      Deallocate ( T_incr_diagnostic )
!
!     Deallocate the diagnostic space that is no longer required.
!
      If (LW_diag%L_total_cloud_cover)                                  &
     &  Deallocate(LW_diag%total_cloud_cover)
      If (LW_diag%L_clear_olr) Deallocate(LW_diag%clear_olr)
      If (LW_diag%L_surface_down_flux)                                  &
     &  Deallocate(LW_diag%surface_down_flux)
      If (LW_diag%L_surf_down_clr)                                      &
     &  Deallocate(LW_diag%surf_down_clr)
      If (LW_diag%L_clear_hr) Deallocate(LW_diag%clear_hr)
      If (LW_diag%L_net_flux_trop)                                      &
     &  Deallocate(LW_diag%net_flux_trop)
      If (LW_diag%L_down_flux_trop)                                     &
     &  Deallocate(LW_diag%down_flux_trop)
      If (LW_diag%L_total_cloud_on_levels)                              &
     &  Deallocate(LW_diag%total_cloud_on_levels)
!
! Deallocate the extinction and absorptivity diagnostics
!  Note that the deallocation of the (SW) extinction diagnostics was
!  delayed as they are required to calculate the ISCCP diagnostics.
      If (LW_diag%L_cloud_absorptivity)                                 &
     &  Deallocate(LW_diag%cloud_absorptivity)
      If (LW_diag%L_cloud_weight_absorptivity)                          &
     &  Deallocate(LW_diag%cloud_weight_absorptivity)
      If (LW_diag%L_ls_cloud_absorptivity)                              &
     &  Deallocate(LW_diag%ls_cloud_absorptivity)
      If (LW_diag%L_ls_cloud_weight_absorptivity)                       &
     &  Deallocate(LW_diag%ls_cloud_weight_absorptivity)
      If (LW_diag%L_cnv_cloud_absorptivity)                             &
     &  Deallocate(LW_diag%cnv_cloud_absorptivity)
      If (LW_diag%L_cnv_cloud_weight_absorptivity)                      &
     &  Deallocate(LW_diag%cnv_cloud_weight_absorptivity)
      If (SW_diag%L_cloud_extinction)                                   &
     &  Deallocate(SW_diag%cloud_extinction)
      If (SW_diag%L_cloud_weight_extinction)                            &
     &  Deallocate(SW_diag%cloud_weight_extinction)
      If (SW_diag%L_ls_cloud_extinction)                                &
     &  Deallocate(SW_diag%ls_cloud_extinction)
      If (SW_diag%L_ls_cloud_weight_extinction)                         &
     &  Deallocate(SW_diag%ls_cloud_weight_extinction)
      If (SW_diag%L_cnv_cloud_extinction)                               &
     &  Deallocate(SW_diag%cnv_cloud_extinction)
      If (SW_diag%L_cnv_cloud_weight_extinction)                        &
     &  Deallocate(SW_diag%cnv_cloud_weight_extinction)

! Deallocate the isccp diagnostics
      If (LW_diag%L_isccp_weights)                                      &
     &  Deallocate(LW_diag%isccp_weights)
      If (LW_diag%L_isccp_cf)                                           &
     &  Deallocate(LW_diag%isccp_cf)
      If (LW_diag%L_isccp_cf_tau_0_to_p3)                               &
     &  Deallocate(LW_diag%isccp_cf_tau_0_to_p3)
      If (LW_diag%L_isccp_cf_tau_p3_to_1p3)                             &
     &  Deallocate(LW_diag%isccp_cf_tau_p3_to_1p3)
      If (LW_diag%L_isccp_cf_tau_1p3_to_3p6)                            &
     &  Deallocate(LW_diag%isccp_cf_tau_1p3_to_3p6)
      If (LW_diag%L_isccp_cf_tau_3p6_to_9p4)                            &
     &  Deallocate(LW_diag%isccp_cf_tau_3p6_to_9p4)
      If (LW_diag%L_isccp_cf_tau_9p4_to_23)                             &
     &  Deallocate(LW_diag%isccp_cf_tau_9p4_to_23)
      If (LW_diag%L_isccp_cf_tau_23_to_60)                              &
     &  Deallocate(LW_diag%isccp_cf_tau_23_to_60)
      If (LW_diag%L_isccp_cf_tau_ge_60)                                 &
     &  Deallocate(LW_diag%isccp_cf_tau_ge_60)
      If (LW_diag%L_meanalbedocld)                                      &
     &  Deallocate(LW_diag%meanalbedocld)
      If (LW_diag%L_meantaucld)                                         &
     &  Deallocate(LW_diag%meantaucld)
      If (LW_diag%L_meanptop)                                           &
     &  Deallocate(LW_diag%meanptop)
      If (LW_diag%L_totalcldarea)                                       &
     &  Deallocate(LW_diag%totalcldarea)
! Deallocate the grid-box mean cloud diagnostics
      If (LW_diag%L_ls_qcl_rad)                                         &
     &  Deallocate(LW_diag%ls_qcl_rad)
      If (LW_diag%L_ls_qcf_rad)                                         &
     &  Deallocate(LW_diag%ls_qcf_rad)
      If (LW_diag%L_cc_qcl_rad)                                         &
     &  Deallocate(LW_diag%cc_qcl_rad)
      If (LW_diag%L_cc_qcf_rad)                                         &
     &  Deallocate(LW_diag%cc_qcf_rad)
      If (LW_diag%L_ls_cl_rad)                                          &
     &  Deallocate(LW_diag%ls_cl_rad)
      If (LW_diag%L_ls_cf_rad)                                          &
     &  Deallocate(LW_diag%ls_cf_rad)
      If (LW_diag%L_cc_cl_rad)                                          &
     &  Deallocate(LW_diag%cc_cl_rad)
      If (LW_diag%L_cc_cf_rad)                                          &
     &  Deallocate(LW_diag%cc_cf_rad)
! Deallocate the aerosol optical depth diagnostics
      If (LW_diag%L_aod_sulphate)                                       &
     &  Deallocate(LW_diag%aod_sulphate)
      If (LW_diag%L_aod_dust)                                           &
     &  Deallocate(LW_diag%aod_dust)
      If (LW_diag%L_aod_seasalt)                                        &
     &  Deallocate(LW_diag%aod_seasalt)
      If (LW_diag%L_aod_soot)                                           &
     &  Deallocate(LW_diag%aod_soot)
      If (LW_diag%L_aod_biomass)                                        &
     &  Deallocate(LW_diag%aod_biomass)
      If (LW_diag%L_aod_biogenic)                                       &
     &  Deallocate(LW_diag%aod_biogenic)
      If (LW_diag%L_aod_ocff)                                           &
     &  Deallocate(LW_diag%aod_ocff)
      If (LW_diag%L_aod_delta)                                          &
     &  Deallocate(LW_diag%aod_delta)
!
!     Deallocate the minimal space for the diagnostics which have
!     not been selected.
      If (.NOT.LW_diag%L_total_cloud_cover)                             &
     &  Deallocate(LW_diag%total_cloud_cover)
      If (.NOT.LW_diag%L_clear_olr) Deallocate(LW_diag%clear_olr)
      If (.NOT.LW_diag%L_surface_down_flux)                             &
     &  Deallocate(LW_diag%surface_down_flux)
      If (.NOT.LW_diag%L_surf_down_clr)                                 &
     &  Deallocate(LW_diag%surf_down_clr)
      If (.NOT.LW_diag%L_clear_hr) Deallocate(LW_diag%clear_hr)
      If (.NOT.LW_diag%L_net_flux_trop)                                 &
     &  Deallocate(LW_diag%net_flux_trop)
      If (.NOT.LW_diag%L_down_flux_trop)                                &
     &  Deallocate(LW_diag%down_flux_trop)
      If (.NOT.LW_diag%L_total_cloud_on_levels)                         &
     &  Deallocate(LW_diag%total_cloud_on_levels)
      If (.NOT.LW_diag%L_cloud_absorptivity)                            &
     &  Deallocate(LW_diag%cloud_absorptivity)
      If (.NOT.LW_diag%L_cloud_weight_absorptivity)                     &
     &  Deallocate(LW_diag%cloud_weight_absorptivity)
      If (.NOT.LW_diag%L_ls_cloud_absorptivity)                         &
     &  Deallocate(LW_diag%ls_cloud_absorptivity)
      If (.NOT.LW_diag%L_ls_cloud_weight_absorptivity)                  &
     &  Deallocate(LW_diag%ls_cloud_weight_absorptivity)
      If (.NOT.LW_diag%L_cnv_cloud_absorptivity)                        &
     &  Deallocate(LW_diag%cnv_cloud_absorptivity)
      If (.NOT.LW_diag%L_cnv_cloud_weight_absorptivity)                 &
     &  Deallocate(LW_diag%cnv_cloud_weight_absorptivity)
      If (.NOT.SW_diag%L_cloud_extinction)                              &
     &  Deallocate(SW_diag%cloud_extinction)
      If (.NOT.SW_diag%L_cloud_weight_extinction)                       &
     &  Deallocate(SW_diag%cloud_weight_extinction)
      If (.NOT.SW_diag%L_ls_cloud_extinction)                           &
     &  Deallocate(SW_diag%ls_cloud_extinction)
      If (.NOT.SW_diag%L_ls_cloud_weight_extinction)                    &
     &  Deallocate(SW_diag%ls_cloud_weight_extinction)
      If (.NOT.SW_diag%L_cnv_cloud_extinction)                          &
     &  Deallocate(SW_diag%cnv_cloud_extinction)
      If (.NOT.SW_diag%L_cnv_cloud_weight_extinction)                   &
     &  Deallocate(SW_diag%cnv_cloud_weight_extinction)
           If (.NOT.LW_diag%L_isccp_weights)                            &
     &  Deallocate(LW_diag%isccp_weights)
      If (.NOT.LW_diag%L_isccp_cf)                                      &
     &  Deallocate(LW_diag%isccp_cf)
      If (.NOT.LW_diag%L_isccp_cf_tau_0_to_p3)                          &
     &  Deallocate(LW_diag%isccp_cf_tau_0_to_p3)
      If (.NOT.LW_diag%L_isccp_cf_tau_p3_to_1p3)                        &
     &  Deallocate(LW_diag%isccp_cf_tau_p3_to_1p3)
      If (.NOT.LW_diag%L_isccp_cf_tau_1p3_to_3p6)                       &
     &  Deallocate(LW_diag%isccp_cf_tau_1p3_to_3p6)
      If (.NOT.LW_diag%L_isccp_cf_tau_3p6_to_9p4)                       &
     &  Deallocate(LW_diag%isccp_cf_tau_3p6_to_9p4)
      If (.NOT.LW_diag%L_isccp_cf_tau_9p4_to_23)                        &
     &  Deallocate(LW_diag%isccp_cf_tau_9p4_to_23)
      If (.NOT.LW_diag%L_isccp_cf_tau_23_to_60)                         &
     &  Deallocate(LW_diag%isccp_cf_tau_23_to_60)
      If (.NOT.LW_diag%L_isccp_cf_tau_ge_60)                            &
     &  Deallocate(LW_diag%isccp_cf_tau_ge_60)
      If (.NOT.LW_diag%L_meanalbedocld)                                 &
     &  Deallocate(LW_diag%meanalbedocld)
      If (.NOT.LW_diag%L_meantaucld)                                    &
     &  Deallocate(LW_diag%meantaucld)
      If (.NOT.LW_diag%L_meanptop)                                      &
     &  Deallocate(LW_diag%meanptop)
      If (.NOT.LW_diag%L_totalcldarea)                                  &
     &  Deallocate(LW_diag%totalcldarea)
      If (.NOT.LW_diag%L_aod_sulphate)                                  &
     &  Deallocate(LW_diag%aod_sulphate)
      If (.NOT.LW_diag%L_aod_dust)                                      &
     &  Deallocate(LW_diag%aod_dust)
      If (.NOT.LW_diag%L_aod_seasalt)                                   &
     &  Deallocate(LW_diag%aod_seasalt)
      If (.NOT.LW_diag%L_aod_soot)                                      &
     &  Deallocate(LW_diag%aod_soot)
      If (.NOT.LW_diag%L_aod_biomass)                                   &
     &  Deallocate(LW_diag%aod_biomass)
      If (.NOT.LW_diag%L_aod_biogenic)                                  &
     &  Deallocate(LW_diag%aod_biogenic)
      If (.NOT.LW_diag%L_aod_ocff)                                      &
     &  Deallocate(LW_diag%aod_ocff)
      If (.NOT.LW_diag%L_aod_delta)                                     &
     &  Deallocate(LW_diag%aod_delta)
!

!

      End If ! on error_code

! end of routine NI_rad_ctl

      Return
      END SUBROUTINE Glue_Rad
