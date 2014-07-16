#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
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
     &, L_Rad_Step_diag, L_Rad_Step_prog                                &
     &, L_Forcing, L_Timestep, L_Radiance, L_Wenyi                      &
     &, L_CAL360, L_SEC_VAR, L_EqT                                      &
     &, L_MICROPHY, L_INHOM_CLOUD, L_USE_DUST, L_USE_BIOGENIC           &
     &, L_use_sulpc_direct, L_use_soot_direct, L_use_soot_indirect      &
     &, L_use_bmass_direct, L_use_bmass_indirect                        &
     &, L_use_ocff_direct, L_use_ocff_indirect                          &
     &, L_use_sulpc_indirect_SW, L_use_sulpc_indirect_LW                &
     &, L_emcorr,L_climat_aerosol,L_clim_aero_hgt,L_HadGEM1_Clim_Aero   &
     &, Ltimer, L_ssice_albedo                                          &
     &, L_MOSES_II, L_snow_albedo                                       &
     &, L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a              &
     &, L_use_seasalt_indirect                                          &
     &, L_use_seasalt_direct                                            &
     &, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                           &
     &, L_pc2, l_mixing_ratio                                           &
     &, L_MURK_RAD                                                      &
     &, L_rad_deg                                                       &
     &, L_co2_interactive,                                              &
!kdcorbin, 05/10 - added L_co2_radiation
     &, L_co2_radiation                                                 &
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
     &, timestep,radiation_timestep                                     &
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
#include "argsts.h"
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
     &, frac_control                                                    &
! chemical greenhouse gas fields
     &, ngrgas, grgas_field                                             &

! ancillary fields and fields needed to be kept from timestep to
! timestep
     &, snow_depth, snow_depth_sea, ice_fract, rgrain, soot             &
     &, cca, ccb, cct, cclwp,ccw,lcbase                                 &
     &, ozone, SW_incs, LW_incs, dirpar_inc                             &
     &, O3_trop_level, O3_trop_height                                   &
     &, T_trop_level, T_trop_height, zh                                 &
     &, land_index, albsoil, lai, snow_tile                             &
     &, frac, tstar_tile, z0_tile                                       &
     &, dOLR_rts, LW_down, SW_tile_rts                                  &
     &, u_1, v_1                                                        &
     &, land_alb,sice_alb                                               &
     &, ES_SPACE_INTERP, A_SW_radstep, A_LW_radstep                     &
     &, A_SW_radstep_diag,A_SW_radstep_prog                             &
     &, A_LW_radstep_diag, A_LW_radstep_prog, RAD_MASK                  &

! in/out
     &, T_n, q_n, qcl_n, qcf_n, cf_n, cfl_n, cff_n                      &
     &, qcf2_n, qrain_n, qgraup_n                                       &
     &, T_inc, q_inc, qcl_inc, cf_inc, cfl_inc                          &
     &, sum_eng_fluxes                                                  &
     &, cos_zenith_angle                                                &

! out.
     &, photosynth_act_rad, surf_radflux, rad_hr, dOLR, SW_tile         &
! Section Information
     &, maxsects, h_sect                                                &
! error information
     &, Error_code  )


! Modules required for multiple calls to radiation

      USE max_calls
      USE diagoffset

! Modules required for reading in spectral files and control options
! as structures

      USE spec_sw_lw
      USE sw_control_struct
      USE lw_control_struct

! Modules required for radiative forcing

      USE mm_ratios
      USE coradoca

! Modules required for the orography scheme

      USE solinc_data, Only:                                            &
        sol_bearing, orog_corr, f_orog,                                 &
        slope_aspect, slope_angle, L_orog

! Modules for diagnostics

      USE sw_diag_mod
      USE lw_diag_mod
      Use volcts_Mod
      USE csenario_mod
! get SC from namelist of coupling control      
!      use auscom_cpl_data_mod, only: SC 
!dhb599 20110812: get SC
#include "swsc.h"

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
!dhb599 20110812: SC is passed in by 'swsc.h' (see above)
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
     &  can_rad_mod                                                      
                                         ! Which canopy radiation model
                                         ! we're using
      Integer                                                           &
     &  model_domain

!       LOGICAL, Parameter :: l_surfalb_forcing=.true.

      Logical                                                           &
     &  L_Rad_Step                                                      &
                             ! dummy variable in versions 3C/3Z
     &, L_Rad_Step_diag                                                 &
                             ! true if fast radiation timestep
     &, L_Rad_Step_prog                                                 &
                             ! true if slow radiation timestep
     &, L_Timestep                                                      &
                             ! true if improved timestepping is used
     &, L_Forcing                                                       &
                             ! true if radiative frocing required
     &, L_Radiance                                                      &
                             ! true if radiances are to be calculated
     &, L_Wenyi                                                         &
                        ! true if Wenyi P & T scaling required 
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
              ! Controls radiation interaction with CO2 - kdcorbin, 06/10
     &, L_use_stochem_CH4                                               &
                           ! Control the use of STOCHEM CH4 field
     &, L_use_clearrh                                                   &
                           ! RH used for hygroscopic aerosols
     &, L_ukca                                                          &
                           ! Switch for UKCA sub-model
     &, L_USE_SPEC_SEA                                                  &
                       ! Switch for spectrally dependent sea albedos
     &, L_MOD_BARKER_ALBEDO						& 
 			   ! Use modified Barker albedo (open sea)
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
                                     !Dummy variable in version 3C/3Z
     &, radiation_tstep                                                 &
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
     &, max_trop_level  ! Highest permitted level for the tropopause
!                       ! used for radiative purposes.

      Integer                                                           &
     &  A_SW_segments                                                   &
     &, A_SW_seg_size                                                   &
     &, A_LW_segments                                                   &
     &, A_LW_seg_size                                                   &
     &, A_SW_radstep                                                    &
                           ! Dummy variable in version 3C/3Z
     &, A_LW_radstep                                                    &
                           ! Dummy variable in version 3C/3Z
     &, A_SW_radstep_diag                                               &
     &, A_SW_radstep_prog                                               &
     &, A_LW_radstep_diag                                               &
     &, A_LW_radstep_prog                                               &

     &, aero_bl_levels
!             Common number of layers taken to be occupied by the
!             boundary-layer aerosol if the boundary layer
!             depth is not used to determine the number separately
!             at each grid-point

#include "dimfix3a.h"
      ! Scaling factors to simulate inhomogeneous cloud.
      Real, Dimension(NPD_CLOUD_COMPONENT) :: INHOM_CLOUD_SW
      Real, Dimension(NPD_CLOUD_COMPONENT) :: INHOM_CLOUD_LW

      ! Decorrelation pressure scale for large scale cloud
      Real :: DP_CORR_STRAT
      ! Decorrelation pressure scale for convective cloud
      Real :: DP_CORR_CONV

#include "csubmodl.h"
#include "typsts.h"
#include "nstypes.h"
!     Missing number indicators
#include "c_mdi.h"

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
     &, ES_SPACE_INTERP(4, row_length, rows)
!              Coeffs for spatial interpolation of radiation quantities

! Aerosol climatology for NWP

#include "arcl_dim.h"
      
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
     &, true_latitude(row_length,rows)                                  &
     &, FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)

#if defined(SCMA)
! Include parameters necessary for calls to SCMoutput...
#include "s_scmop.h"
#endif

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

! Additional variables for SCM diagnostics which are dummy in full UM
      Integer                                                           &
     &  nSCMDpkgs          ! No of diagnostics packages

      Logical                                                           &
                           ! Logicals for diagnostics packages
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
     &, cos_zenith_angle(row_length,rows)                               

! local variables.

#include "fldtype.h"
#include "domtyp.h"

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

#if defined(SCMA)

! The following arrays are required for calculating increments
! and heating rates for SCM diagnostic output

     &,work_3dw(row_length,rows,wet_model_levels)                       &
     &,heating_rate(row_length,rows,model_levels)                       &
#endif

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


! Radiation fields 1. SW & common with LW.
!
!
      Real                                                              &
     &  day_fraction(row_length,rows,n_swcall)                          &
     &, cos_zenith_angle_2(row_length,rows,n_swcall)                    &
     &, volcmass(row_length,rows)      					&
! VOLCTS expanded to full fields,
!     and converted from optical depth to mass loading.
     &, mean_cos_zenith_angle(row_length,rows)                          &
     &, flandg(row_length, rows)                                        &
     &, land_albedo(row_length, rows, 4)                                &
     &, land_albedo2(row_length,rows,4,n_swcall)                        &
     &, land_and_ice_albedo(row_length, rows)                           &
     &, open_sea_albedo(row_length, rows, 2,n_swcall)                   &
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
     &  List_daylight_points(row_length*rows,n_swcall)                  &
     &, diag_row_list(row_length*rows,n_swcall)                         &
!           List of row indices of points where diagnostics
!           are calculated.
     &, diag_row_list_sw(row_length*rows)                               &
                                           ! Saved version for ISCCP
     &, diag_col_list(row_length*rows,n_swcall)                         &
!           List of column indices of points where diagnostics
!           are calculated.
     &, diag_col_list_sw(row_length*rows)                               &
                                           ! Saved version for ISCCP
     &, trindx(row_length, rows)

      Logical                                                           &
     &  switch(row_length,rows,n_swcall)

      Real                                                              &
     &  Sindec                                                          &
                         !  sin(solar declination)
     &, hour_angle(row_length,rows)                                     &
                                          ! Hour Angle
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
      Logical :: L_CO2_3D       ! Controls use of 3D co2 field
!
      Logical :: L_T_incr_sw    ! Flag for SW temperature increments
!
! SW diagnostics not on stash flag
      Real                                                              &
     &  itoasw(row_length,rows)                                         &
     &, surfsw(row_length,rows)                                         &
     &, land_alb2(row_length,rows,n_swcall)                             &
     &, sice_alb2(row_length,rows,n_swcall)

! SW diagnostics on always true stash flag
      Logical                                                           &
     &  l_flux_below_690nm_surf

      Parameter (l_flux_below_690nm_surf = .true.)

      Real                                                              &
     &  flux_below_690nm_surf(row_length,rows)


! Radiation fields 2. LW

      Real                                                              &
     &  LWsea(row_length, rows)                                         &
     &, OLR(row_length, rows)                                           &
     &, surflw(row_length,rows)                                         &
     &, top_absorption(row_length, rows)
! needed by energy correction
      Real                                                              &
     &  net_atm_flux (row_length, rows)
!
      Logical :: L_T_incr_lw    ! Flag for LW temperature increments
!
! needed for 8A boundary layer
      Real, Allocatable, Save :: alb_tile(:,:,:,:)

      Real                                                              &
     &  frac_tile_alb(land_field, ntype)                                &
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

      Real :: albsolid      !Mean solid surface albedo

!
! Local variables required for multiple time-stepping
!

      LOGICAL :: l_call_swrad
      LOGICAL :: l_call_lwrad

! Local variables required for advancing the model. These are required
! for the improved timestepping algorithm.

      Real, allocatable, save :: SW_incs_local(:,:,:,:)
      Real, allocatable, save :: LW_incs_local(:,:,:,:)
      Real, allocatable, save :: netSW_local(:,:,:)
      Real, allocatable, save :: SWsea_local(:,:,:)
      Real, allocatable, save :: LWsea_local(:,:,:)
      Real, allocatable, save :: top_abs_sw(:,:,:)
      Real, allocatable, save :: top_abs_lw(:,:,:)
      Real, allocatable, save :: olr_local(:,:,:)
      Real, allocatable, save :: flux_b690nm_local(:,:,:)
      Real, allocatable, save :: surf_down_sw_local(:,:,:,:)
      Real, allocatable, save :: dirpar_inc_local(:,:,:)

! Dimension of SW_diag or LW_diag required for forcing,
! improved timestepping or radiances

      INTEGER :: n_swdiag
      INTEGER :: n_lwdiag

! Looping variables for extra calls to radiation

      INTEGER :: j_sw
      INTEGER :: j_lw

! Offset in STASH item number for different sets of diagnostics

      INTEGER :: i_off
!
! Local variables relating to the radiance code
!
! Dynamically defined dimensions used by the radiance code
      INTEGER :: nd_profile
      INTEGER :: nd_field_flux_diag
      INTEGER :: nd_field_rad_diag
      INTEGER :: nd_flux_profile
      INTEGER :: nd_radiance_profile
      INTEGER :: nd_viewing_level
      INTEGER :: nd_direction
      INTEGER :: nd_column
      INTEGER :: nd_cloud_component
      INTEGER :: nd_cloud_type
      INTEGER :: nd_brdf_basis_fnc
      INTEGER :: nd_brdf_trunc
      INTEGER :: nd_channel
      INTEGER :: nd_point_tile
      INTEGER :: nd_tile
      INTEGER :: id_ct
!
!     Generic logicals
!
      LOGICAL :: l_scale_inc
!
!     Mapping of spectral bands to satellite channels
      INTEGER :: n_channel
!       Number of viewing channels
      INTEGER, Pointer, Dimension(:) :: map_channel
!       Mapping of spectral bands to satellite channels
!
!     Satellite viewing geometry
      LOGICAL :: l_viewed
!       Flag set to .TRUE. if point can be viewed by current satellite
      INTEGER :: n_viewing_direction
      INTEGER :: n_viewing_level
      REAL :: viewing_direction(row_length*rows, 1, 2)
      REAL :: viewing_level(model_levels+1)
!
!
! Local Variables required to calculate radiative forcings.
!
! Variables holding the mass mixing ratios of either the
! reference call or the time advancing call.

      REAL                                                              &
     &  j_co2_mmr                                                       &
     &, j_o2_mmr                                                        &
     &, j_n2o_mmr                                                       &
     &, j_ch4_mmr                                                       &
     &, j_cfc11_mmr                                                     &
     &, j_cfc12_mmr                                                     &
     &, j_C113_mmr                                                      &
     &, j_hcfc22_mmr                                                    &
     &, j_hfc125_mmr                                                    &
     &, j_hfc134_mmr

! The variable j_ozone contains the reference ozone field for the
! diagnostic call and the 'normal' ozone field for the time-advancing
! call.

       REAL :: j_ozone(row_length, rows, ozone_levels)

! Same as above for logicals defining the use of aerosols

      LOGICAL                                                           &
     &  j_l_use_sulpc_direct                                            &
     &, j_l_use_seasalt_direct                                          &
     &, j_l_use_soot_direct                                             &
     &, j_l_use_bmass_direct                                            &
     &, j_l_use_ocff_direct

! Are we required to find fluxes at tropopause.

      LOGICAL :: l_find_trop


! SCM Dummy variables to keep call to tropin consistent.
      REAL                                                              &
     & scm_dummy_1d(1,1)                                                &
     &,scm_dummy_2d(1,1,0:model_levels)

#include "csigma.h"
#include "c_0_dg_c.h"

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

! Component IDs for the aerosol climatology for NWP
#include "arcl_ids.h"

! Initialise n_channel as it is not set in the 3C scheme.
      n_channel=1

!---------------------------------------------------------------------
! Section: CHECKS. Check that flags and timesteps are in agreement
!---------------------------------------------------------------------
      if (l_forcing) then
        if ( (A_SW_Radstep_prog >  A_SW_Radstep_diag).or.               &
     &       (A_LW_Radstep_prog >  A_LW_Radstep_diag) ) then

! DEPENDS ON: ereport
           Call Ereport("NI_RAD_CTL",20,                                &
     &       "If radiative forcing is required the diagnostic" //       &
     &       "timestep has to be larger or equal to the "   //          &
     &       "prognostic timestep.")
        endif
        if ((mod(A_SW_Radstep_diag,A_SW_Radstep_prog) /= 0).or.         &
     &       (mod(A_SW_Radstep_diag,A_SW_Radstep_prog) /= 0)) then

! DEPENDS ON: ereport
          Call Ereport("NI_RAD_CTL",20,                                 &
     &      "The diagnostic(fast) timestep needs to be a" //            &
     &      "a multiple of the prognostic (slow) timestep.")
        endif

      endif

      if (l_timestep) then
        if ( (A_SW_Radstep_prog <  A_SW_Radstep_diag).or.               &
     &       (A_LW_Radstep_prog <  A_LW_Radstep_diag) ) then

! DEPENDS ON: ereport
          Call Ereport("NI_RAD_CTL",20,                                 &
     &      "If improved timestepping  is required the diagnostic" //   &
     &      "timestep has to be smaller than the prognostic"   //       &
     &      "timestep.")
        endif
        if ((mod(A_SW_Radstep_prog,A_SW_Radstep_diag) /= 0).or.         &
     &       (mod(A_SW_Radstep_prog,A_SW_Radstep_diag) /= 0)) then

! DEPENDS ON: ereport
          Call Ereport("NI_RAD_CTL",20,                                 &
     &      "The prognostic(slow) timestep needs to be a" //            &
     &      "a multiple of the diagnostic (fast) timestep.")
        endif
      endif

      if ((n_swcall >  npd_swcall).or.(n_lwcall >  npd_lwcall))         &
! DEPENDS ON: ereport
     &    call Ereport("NI_RAD_CTL",19,"Too may calls to radiation")


! ----------------------------------------------------------------------
! Section INI. Initialisation of stash output variables.
! ----------------------------------------------------------------------

! kdcorbin, 06/10 - added switch to turn on/off radiation interactions
      !L_CO2_3D = L_CO2_interactive
      L_CO2_3D = L_CO2_RADIATION

! Map the stashflags to the internal controlling logicals.
! Initialise all the internal controlling logicals for the diagnostics
!
#if !defined(SCMA)
      L_T_incr_sw                 = (sf(181,1) .or. sf(161,1))
      L_T_incr_lw                 = (sf(181,2) .or. sf(161,2))
#else
      L_T_incr_sw                      = .true.
      L_T_incr_lw                      = .true.
#endif
!
! Set up the logicals for the diagnostics.
!
! Set STASH flags to false if orography scheme is off.
!
#if !defined(SCMA)
      if (.not.L_orog) sf(292:296,1) = .false.
#endif

! If the forcing or the improved timestepping is not required then
! the diagnostic timestep is equal to the prognostic timestep.

      if ((.not.l_forcing).and.(.not.l_timestep)) then
         L_rad_step_diag=L_rad_step_prog
      endif
      if (l_forcing.or.l_timestep.or.l_radiance) then
         n_swdiag=n_swcall
         n_lwdiag=n_lwcall
      else
         n_swdiag=1
         n_lwdiag=1
      endif

! DEPENDS ON: init_swdiag_logic
      call init_swdiag_logic(1)
! DEPENDS ON: init_lwdiag_logic
      call init_lwdiag_logic(1)

      if (l_rad_step_prog) then

         i_off=0

! The LW logicals need to be set before the SW ones since
! the latter ones depends on the LW ones.

! DEPENDS ON: set_lwdiag_logic
         call set_lwdiag_logic(sf                                       &
     &,      nitems,nsects,l_moses_II,l_ctile                           &
     &,      1, i_off)

! DEPENDS ON: set_swdiag_logic
         call set_swdiag_logic(sf                                       &
     &,      nitems,nsects,l_moses_II,l_ctile                           &
     &,      L_flux_below_690nm_surf                                    &
     &,      1, i_off)

      endif

      if (l_radiance) then

         do j_lw=2,n_lwcall

            i_off=90+(diagnostic_offset/20)*(j_lw-1)

! DEPENDS ON: init_lwdiag_logic
            call init_lwdiag_logic(j_lw)

! DEPENDS ON: set_lwdiag_logic
            call set_lwdiag_logic(sf                                    &
     &,         nitems,nsects,l_moses_II,l_ctile                        &
     &,         j_lw, i_off)

         enddo
         do j_sw=2,n_swcall

            i_off=90+(diagnostic_offset/20)*(j_sw-1)

! DEPENDS ON: init_swdiag_logic
            call init_swdiag_logic(j_sw)

! DEPENDS ON: set_swdiag_logic
            call set_swdiag_logic(sf                                    &
     &,         nitems,nsects,l_moses_II,l_ctile                        &
     &,         L_flux_below_690nm_surf                                 &
     &,         j_sw, i_off)
         enddo
      endif

      if (l_forcing.or.l_timestep) then

         if (l_timestep) i_off=0
         if (l_forcing)  i_off=diagnostic_offset

! DEPENDS ON: init_swdiag_logic
         call init_swdiag_logic(2)
! DEPENDS ON: init_lwdiag_logic
         call init_lwdiag_logic(2)

         if (l_rad_step_diag) then

! DEPENDS ON: set_lwdiag_logic
             call set_lwdiag_logic(sf                                   &
     &,         nitems,nsects,l_moses_II,l_ctile                        &
     &,         2, i_off)

! DEPENDS ON: set_swdiag_logic
             call set_swdiag_logic(sf                                   &
     &,         nitems,nsects,l_moses_II,l_ctile                        &
     &,         L_flux_below_690nm_surf                                 &
     &,         2, i_off)

         endif
      endif


! ----------------------------------------------------------------------
! Set tile_pts, tile_index if MOSES II selected.
! Set land_index.
! ----------------------------------------------------------------------
      if (L_MOSES_II)                                                   &
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
      IF(.NOT.L_MOSES_II.AND.L_CTILE)THEN
        write(6,*)'****ERROR***Coastal tiling cannot be '
        write(6,*)'run with pre-MOSES2 sfc/boundary layer code'
        CALL ABORT
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
        Do j = 1, rows
          Do i = 1, row_length
            photosynth_act_rad(i,j) = 0.
            surf_radflux(i,j) = 0.
            flux_below_690nm_surf(i,j)=0.
            sw_net_land(i,j) = 0.0
            sw_net_sice(i,j) = 0.0
            dOLR(i,j) = 0.
            !Alison McLaren advised bug fix here: commented the folloing line
            !which wrong SW radiation over coastal sea-ice cells, namely, it 
            !applys SW Only at the beginning of a radiation step on to these cells
            ! fixed: 20111201.  
            !!! land_alb(i,j)=0.0
            !Note for the ACCESS hg3-C configuration this change won't take effect because
            !     it's the file "src/control/top_level/glue_rad-rad_ctl2.F90"
            !     that is actually used !!!  (2011-12-21)   
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
      If (L_Rad_Step_diag.or.L_rad_step_prog) Then
         height_rho(:,:,:)=r_rho_levels(1:row_length,1:rows,:)
         height_theta(:,:,:)=r_theta_levels(1:row_length,1:rows,:)
      End If
!
! If we use this then more code is required from RADCTL. Namely, if
! climatological aerosols are present or if any radiation calls require
! fluxes at the tropopause then it is necessary to find the tropopause.
!
      l_find_trop = l_climat_aerosol

      do j_sw=1,n_swcall

         l_find_trop=l_find_trop.or.                                    &
     &                  sw_diag(j_sw)%l_net_flux_trop.or.               &
     &                  sw_diag(j_sw)%l_up_flux_trop
      enddo
      do j_lw=1,n_lwcall

         l_find_trop=l_find_trop.or.                                    &
     &                  lw_diag(j_lw)%l_net_flux_trop.or.               &
     &                  lw_diag(j_lw)%l_down_flux_trop
      enddo

      if ((l_rad_step_diag.or.l_rad_step_prog).and.l_find_trop) then

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
!
! Set SCM dummy values to zero
!
          scm_dummy_1d(:,:)   = 0.0
          scm_dummy_2d(:,:,:) = 0.0

! DEPENDS ON: tropin
          Call tropin (T_large, exner_rho_levels, exner_theta_levels,   &
     &                 row_length, rows, model_levels, off_x, off_y,    &
     &                 at_extremity,scm_dummy_1d,scm_dummy_2d,          &
     &                 min_trop_level, max_trop_level, trindx )

        End If

        If ((L_Rad_step_diag.or.L_rad_step_prog) .AND.                  &
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
!-----------------------------------------------------------------------
! Section Rad.0.1 Allocate SW variables and diagnostics
!-----------------------------------------------------------------------

      if (l_rad_step_prog) then

! DEPENDS ON: allocate_sw_diag
        call allocate_sw_diag( row_length, rows, model_levels,          &
     &                 cloud_levels, 1 )

        allocate(alb_tile(land_field,ntiles,4,n_swcall))

        if (l_timestep) then
          allocate(SW_incs_local(row_length, rows,                      &
     &                                0:model_levels+1,n_swcall))
          allocate(surf_down_sw_local(row_length,rows,4,n_swcall))
          allocate(netSW_local(row_length, rows, n_swcall))
          allocate(SWsea_local(row_length, rows ,n_swcall))
          allocate(top_abs_sw(row_length, rows, n_swcall))
          allocate(flux_b690nm_local(row_length,rows,n_swcall))
          allocate(dirpar_inc_local(row_length,rows,n_swcall))
        endif

         if (l_radiance) then
           do j_sw=2,n_swcall

! DEPENDS ON: allocate_sw_diag
             call allocate_sw_diag( row_length, rows,                   &
     &                model_levels, cloud_levels, j_sw )
           enddo
         endif

      endif

      if (l_rad_step_diag) then
        if (l_timestep.or.l_forcing) then
! DEPENDS ON: allocate_sw_diag
            call allocate_sw_diag( row_length, rows, model_levels,      &
     &                    cloud_levels, 2 )
        endif
      endif
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
      Call solvar (PREVIOUS_TIME(1), SCS,                               &
        sw_spectrum(1)%npd_band,                                        &
        sw_spectrum(1)%solar_flux_band,                                 &
        sw_spectrum(1)%rayleigh_coefficient,                            &
        sw_spectrum(1)%w_esft(1+2*sw_spectrum(1)%npd_esft_term          &
          *sw_spectrum(1)%npd_band,1,1),		                &
        .FALSE., cmessage, icode)

      If ( icode .GT. 0 ) Then
        write(6,*) ' icode GT 0 in solvar'
        Return
      Endif
      endif  ! L_SCVARY

!------------------------------------------------------------
! Section Rad.0.2 Do preliminary timestep calculations
!------------------------------------------------------------
       If (L_rad_step_prog.or.L_rad_step_diag) Then
          do j_sw=1,n_swcall

             l_call_swrad=.false.

             if (l_timestep) then

!
! Here we need to be careful with the timestepping. In the
! case of the improved timestepping the prognostic (slow)
! call is done with the 'astronomy' calculated on the
! prognostic timestep and the diagnostic (fast) call is done
! with the astronomy calculated on the diganostic timestep.
! This ensures that all the averaging over different timesteps
! is done accurately.
!
                if ((j_sw==1).and.(l_rad_step_prog)) then
                   radiation_tstep = radiation_tstep_prog
                   frac_tile_alb=frac
                   l_call_swrad=.true.
                endif
                if ((j_sw==2).and.(l_rad_step_diag)) then
                   radiation_tstep = radiation_tstep_diag
                   frac_tile_alb=frac
                   l_call_swrad=.true.
                endif

             elseif (l_forcing) then

!
! In the case of the forcing the astronomy is always calculated
! on the prognostic timescale, i.e. we always use the prognostic
! timestep to calculate the astronomy.
!
                if ((j_sw==1).and.(l_rad_step_prog)) then
                   radiation_tstep = radiation_tstep_prog
                   frac_tile_alb=frac
                   l_call_swrad=.true.
                endif
                if ((j_sw==2).and.(l_rad_step_diag)) then
                   radiation_tstep = radiation_tstep_prog
                   if (C2C_LAND_S(j_sw-1)) then
                      frac_tile_alb=frac_control
                   else
                      frac_tile_alb=frac
                   endif
                   l_call_swrad=.true.
                endif

             elseif (l_radiance) then
                radiation_tstep = radiation_tstep_prog
                frac_tile_alb=frac
                l_call_swrad=.true.
             else
                radiation_tstep = radiation_tstep_prog
                frac_tile_alb=frac
                l_call_swrad=.true.
             endif

             if (l_call_swrad) then

                If (L_orog) Allocate(sol_bearing(row_length,rows))

! DEPENDS ON: solang
                Call solang(                                            &
! input constants
     &             Sindec, seconds_since_midnight,                      &
     &             Radiation_tstep, Eq_Time,                            &
! row and column dependent constants
     &             sin_true_latitude,                                   &
     &             true_longitude,                                      &
! size variables
     &             row_length*rows,                                     &
! output fields
     &             day_fraction(1,1,j_sw) ,                             &
     &             cos_zenith_angle_2(1,1,j_sw),                        &
     &             hour_angle )
!
! Set rounding-error size values to zero - the criterion depends
! on the frequency of full SW calculations because on the physics
! timesteps which are not SW timesteps a test has to be done to
! avoid using the unset data for such points.
! NB: Copied from UM, including comment.
!
                Do j = 1, rows
                   Do i= 1, row_length
                      If ( cos_zenith_angle_2(i,j,j_sw)*                &
     &                                   day_fraction(i,j,j_sw) <       &
     &                     (1.e-10*timestep/radiation_tstep) ) Then
                         cos_zenith_angle_2(i,j,j_sw) = 0.0
                         day_fraction(i,j,j_sw) = 0.0
                      End If
                   End Do
                End Do

! Transfer solar bearing diagnostic before it is deallocated:

                If (SW_diag(j_sw)%L_sol_bearing) Then
                  SW_diag(j_sw)%sol_bearing = sol_bearing
                End If

! Calculate orography correction:

                If (L_orog) Then
                  Allocate(orog_corr(row_length,rows))
! DEPENDS ON: solinc
                  Call solinc(row_length, rows,                         &
     &                        cos_zenith_angle_2(:,:,j_sw))
                  Deallocate(sol_bearing)
                End If

! DEPENDS ON: ftsa
                Call ftsa (                                             &
! input fields
     &            land_sea_mask, flandg, ice_fract,                     &
     &            T_surf, tstar_sice,                                   &
     &            cos_zenith_angle_2(1,1,j_sw),                         &
     &            snow_depth, snow_depth_sea,                           &
! max and min sea ice albedo specifications
     &            SW_alpham, SW_alphac, SW_alphab, SW_dtice,            &
     &            L_MOSES_II,L_ssice_albedo,                            &
     &            L_MOD_BARKER_ALBEDO,                                  &
     &            L_sice_meltponds, L_sice_scattering, L_sice_hadgem1a, &
     &            dt_bare, dalb_bare_wet, pen_rad_frac, SW_beta,        &
! Version of Shortwave
     &            h_sect(1),                                            &
! size and control variables
     &            row_length*rows, row_length*rows,                     &
! output arguments
     &            land_alb2(1,1,j_sw),sice_alb2(1,1,j_sw),              &
     &            open_sea_albedo(1,1,1,j_sw) )
!-----------------------------------------------------------------------
! Calculate MOSES II tile albedos, then reset TILE_PTS and TILE_INDEX
! and set tile fractions to 1 if aggregate tiles are used (NTILES=1).
!-----------------------------------------------------------------------
                if ( L_MOSES_II ) then
!
! Code required to initialize and update snow soot content
!
                   if ( L_snow_albedo ) then
                      Do j = 1, rows
                         Do i = 1, row_length
                            soot(i,j) = 0.
                         End Do
                      End Do
                   endif


! DEPENDS ON: tile_albedo
                   call tile_albedo ( row_length*rows,                  &
     &                land_field,land_index,ntiles,tile_pts,            &
     &                tile_index,L_snow_albedo,albsoil,                 &
     &                cos_zenith_angle_2(1,1,j_sw),frac_tile_alb,       &
     &                lai,rgrain, snow_tile,soot,tstar_tile,            &
     &                z0_tile, alb_tile(1,1,1,j_sw),                    &
     &                land_albedo2(1,1,1,j_sw),can_rad_mod )


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
                endif           ! if MOSES II


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
        cmessage = 'WARNING: VOLCANIC FORCING'
	WRITE(6,*) 'Time varying volcanic forcing has been selected but'
	WRITE(6,*) 'the check that the run is global has been disabled.'
        WRITE(6,*) 'Only use this option for global simulations.'
!DEPENDS ON: ereport
        CALL EREPORT('NI_rad_ctl3c', -1, cmessage) 

!     If date is outside range of given data then set VOLCMASS
!     to control value.
        IF ( (PREVIOUS_TIME(1) .LT. 1850) .OR. 			&
     &       (PREVIOUS_TIME(1) .GT. 2300) ) THEN
         DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            VOLCMASS(I,J) = 1.86604E-6
          ENDDO
         ENDDO
        ELSE
         DO J = 1, ROWS
! We assume rows are latitude rows, which could be checked & trapped
! A check is made that GLOBAL is defined. So rows should be latitude
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


             endif              ! if sw radiation call required
          enddo                 ! loop over radiation calls
       endif                    ! if radiation timestep
!
! To obtain reproducible results independent of the
! decomposition of the domain used on an MPP machine a global
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

       If (L_rad_step_prog.or.L_rad_step_diag) Then
          do j_sw=n_swcall,1,-1

! For each SW call set the correct mass mixing ratios and switch the
! aerosols on/off.
! Note that in the case of the aerosols the reference state is
! 'no aerosols' and hence they are switched off in the diagnostic call.

! Since the "default" for each variable is to have it set in the
! diagnostic call to the same as in the prognostic, we set them all
! unconditionally to the prognostic & then re-set as needed.

             j_co2_mmr = CO2_MMR
             j_o2_mmr  = O2MMR
             j_ozone = ozone  ! Note this is an array assignment
             j_l_use_sulpc_direct   = L_use_sulpc_direct
             j_l_use_seasalt_direct = L_use_seasalt_direct
             j_l_use_soot_direct    = L_use_soot_direct
             j_l_use_bmass_direct   = L_use_bmass_direct
             j_l_use_ocff_direct    = L_use_ocff_direct

             l_call_swrad=.false.

!
! if Improved timestepping required
!
             if (L_timestep) then
                if ((j_sw==1).and.(l_rad_step_prog)) then
                   a_sw_radstep=a_sw_radstep_prog
                   land_alb=land_alb2(:,:,j_sw)
                   sice_alb=sice_alb2(:,:,j_sw)
                   land_albedo=land_albedo2(:,:,:,j_sw)
                   l_call_swrad=.true.
                endif
                if ((j_sw==2).and.(l_rad_step_diag)) then
                   a_sw_radstep=a_sw_radstep_diag
                   land_alb=land_alb2(:,:,j_sw)
                   sice_alb=sice_alb2(:,:,j_sw)
                   land_albedo=land_albedo2(:,:,:,j_sw)
                   l_call_swrad=.true.
                endif
!
! if Forcing required
!
             elseif (L_Forcing) then
!
! Decide if it is necessary to call radiation
!
                if ((j_sw==1).and.(l_rad_step_prog)) then
                   a_sw_radstep=a_sw_radstep_prog
                   land_alb=land_alb2(:,:,j_sw)
                   sice_alb=sice_alb2(:,:,j_sw)
                   land_albedo = land_albedo2(:,:,:,j_sw)
                   l_call_swrad=.true.
                endif
                if ((j_sw==2).and.(l_rad_step_diag)) then
                   a_sw_radstep=a_sw_radstep_prog
                   land_alb=land_alb2(:,:,j_sw)
                   sice_alb=sice_alb2(:,:,j_sw)
                   land_albedo = land_albedo2(:,:,:,j_sw)
                   l_call_swrad=.true.
                endif
!
                IF ( J_SW  ==  1 ) THEN        ! Prognostic Call
!
! Switch back to time varying land use for prognostic call
!
!                   if (l_surfalb_forcing) then
!                      land_albedo = land_albedo2(:,:,:,1)
!                   endif

                ELSE                ! J_SW > 1: Diagnostic Call
!
! Forcing due to land_surface.  Set dummy variable to current land
! albedo for use later in prognostic call.
!
!                   IF (l_surfalb_forcing) then
!                      land_albedo = land_albedo2(:,:,:,2)
!                   ENDIF
!
! Gases
!
                   IF ( C2C_CO2(J_SW-1) ) j_co2_mmr = CO2_MMR_D
                   IF ( C2C_O2(J_SW-1) )  j_o2_mmr  = O2_MMR_D
!
! Ozone: Note that aerosol contains the MURK array.  It is in fact what
! is called d1(MURK) in Simon Tett's old vn4.5 mod.  Care needs to be
! taken so that MURK has the correct dimensions before it arrives in
! RAD_CTL2.
!
                   IF ( C2C_O3(J_SW-1) ) j_ozone = aerosol
!
! Aerosols
!
                   IF ( C2C_SULPC_D(J_SW-1) )                           &
     &                j_l_use_sulpc_direct   = .false.
                   IF ( C2C_SEAS_D(J_SW-1) )                            &
     &                j_l_use_seasalt_direct = .false.
                   IF ( C2C_SOOT_D(J_SW-1) )                            &
     &                j_l_use_soot_direct    = .false.
                   IF ( C2C_BMB_D(J_SW-1) )                             &
     &                j_l_use_bmass_direct   = .false.
                   IF ( C2C_OCFF_D(J_SW-1) )                            &
     &                j_l_use_ocff_direct    = .false.

                ENDIF
             elseif (l_radiance) then
                a_sw_radstep=a_sw_radstep_prog
                land_alb=land_alb2(:,:,j_sw)
                sice_alb=sice_alb2(:,:,j_sw)
                land_albedo=land_albedo2(:,:,:,j_sw)
                l_call_swrad=.true.
             else
                a_sw_radstep=a_sw_radstep_prog
                land_alb=land_alb2(:,:,j_sw)
                sice_alb=sice_alb2(:,:,j_sw)
                land_albedo=land_albedo2(:,:,:,j_sw)
                l_call_swrad=.true.
             endif
!
! If call to SW Radiation required
!
             if (l_call_swrad) then

! DEPENDS ON: prelim_swrad
                call prelim_swrad(error_code                            &
     &,            at_extremity, n_proc                                 &
! Model Dimensions
     &,            row_length,rows,model_levels                         &
! Model Switches
     &,            model_domain,mt_global,L_rad_deg                     &
     &,            SW_Control(j_sw)%L_subsample                         &
     &,            SW_Control(j_sw)%L_geostationary                     &
! Time stepping Information
     &,            timestep_number,a_sw_radstep                         &
! Satellite Geometry
     &,            SW_Control(j_sw)%min_view_lon                        &
     &,            SW_Control(j_sw)%max_view_lon                        &
     &,            SW_Control(j_sw)%min_view_lat                        &
     &,            SW_Control(j_sw)%max_view_lat                        &
! Number of Call
     &,            j_sw                                                 &
! Other variables
     &,            sin_true_latitude,true_longitude                     &
     &,            seconds_since_midnight                               &
     &,            tot_daylight_points,daylight_points                  &
     &,            day_fraction(1,1,j_sw)                               &
     &,            list_daylight_points(1,j_sw), rad_mask               &
     &,            switch(1,1,j_sw), first_data_interp                  &
     &,            first_data_interp_sw                                 &
     &,            diag_row_list(1,j_sw), diag_row_list_sw              &
     &,            diag_col_list(1,j_sw), diag_col_list_sw )

! Zero the calculated outputs: this will later simplify the
! treatment of the case when there are no lit points.

                SW_incs(:,:,:)=0.0
                If (SW_diag(j_sw)%L_direct_par) dirpar_inc(:,:) = 0.0
                netSW (:,:) = 0.0
                SWsea (:,:) = 0.0
                If (SW_Control(1)%l_extra_top) Then
                   top_absorption(:,:) = 0.0
                Endif
                surf_down_sw (:,:,:) = 0.0


                If ( daylight_points  >   0 ) Then
!
! Calculate length of segments into which to split
! short wave calculations.
!
                   If ( A_SW_seg_size > 0 ) Then
                     ! We are specifying size of segments
                     step = A_SW_seg_size
                     A_SW_segments = Ceiling ( Real(daylight_points) /  &
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
!
! Short wave radiation calculations called in segments
! Code left over from parallel version.
!
                   first_point_temp(1) = 1
                   Do i = 1, N_SW_segments
                      lit_points = step
                      start_point = 1 + (i-1) *step
                      seg_points = list_daylight_points(i*step,j_sw)    &
     &                                              - seg_start + 1
                      If ( i  ==  N_SW_segments ) Then
                         lit_points = daylight_points                   &
     &                                   - step* (N_SW_segments - 1)
                         seg_points = row_length*rows-seg_start+1
                      End If
                      Do j = start_point, start_point+lit_points-1
                         list_daylight_points(j,j_sw)                   &
     &                     =list_daylight_points(j,j_sw)-seg_start+1
                      End Do
                      seg_points_temp(i) = seg_points
                      If ( i  <   N_SW_segments ) Then
                         first_point_temp(i+1) = first_point_temp(i)    &
     &                                                   +seg_points
                      End If
                      seg_start = seg_start + seg_points
                   End Do

                   Do i = 1, N_SW_segments
!
! Save values for later use in call to ISCCP code
!
                      first_point_temp_sw(i)=first_point_temp(i)
                      seg_points_temp_sw(i)=seg_points_temp(i)
                      step_sw=step

                      lit_points = step
                      start_point = 1+(i-1)*step
                      If ( i  ==  N_SW_segments ) Then
                         lit_points = daylight_points - step            &
     &                                          * (N_SW_segments - 1)
                      End If
                      first_point = first_point_temp(i)
!
! Set the first point of the dust arrays to be used.
!
                      IF (L_USE_DUST) THEN
                         FIRST_POINT_DUST=FIRST_POINT
                      ELSE
                        FIRST_POINT_DUST=1
                      ENDIF

! Set the first point of the biogenic array.
                      IF (L_USE_BIOGENIC) THEN
                        FIRST_POINT_BIOGENIC=FIRST_POINT
                      ELSE
                        FIRST_POINT_BIOGENIC=1
                      ENDIF
!
! Set the first point of the array of sulphate to be used.
! A separate assignment is necessary since this array will
! not be of the full size unless the sulphur cycle is on.
!
                      IF (L_USE_SULPC_DIRECT .OR.                       &
     &                                 L_USE_SULPC_INDIRECT_SW) THEN
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
                      If (L_use_bmass_direct .or.                       &
     &                                      L_use_bmass_indirect) then
                         first_point_biomass=first_point
                      Else
                         first_point_biomass=1
                      EndIf
!
                      If (L_use_ocff_direct .or.                        &
     &                    L_use_ocff_indirect) then
                         first_point_ocff=first_point
                      Else
                         first_point_ocff=1
                      End If
!
                      If (L_use_seasalt_indirect .OR.                   &
     &                                      L_use_seasalt_direct) Then
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
! Set the actual size of arrays in the radiation code:
! for some architectures (such as that of Cray vector
! machines) on odd size is preferred to avoid memory
! bank conflicts.
!
                      nd_rad_pts=2*(lit_points/2)+1
!
! Set the number of layers seen in the radiation code.
! This may optionally be 1 greater than the number used
! in the rest of the model to avoid spurious effects
! resulting from the upper boundary (principally in
! stratospheric configurations).
!
                      If (SW_Control(1)%l_extra_top) Then
                         n_rad_layers=model_levels+1
                      else
                         n_rad_layers=model_levels
                      endif


                      If (h_sect(1) == '03Z') Then

! Set dynamically determined array sizes required for
! radiation calculations. There is redundancy here
! because the viewing direction is calculated for the
! whole PE within each segement.

                        Allocate(map_channel(SW_Spectrum(j_sw)%n_band))

! DEPENDS ON: r2_set_rad_dim
                        Call r2_set_rad_dim(row_length, rows            &
     & ,                   sin_true_latitude, true_longitude            &
     & ,                   .TRUE., sindec, hour_angle                   &
     & ,                   SW_Control(j_sw)%l_geostationary             &
     & ,                   SW_Control(j_sw)%sat_hgt                     &
     & ,                   SW_Control(j_sw)%sat_lat                     &
     & ,                   SW_Control(j_sw)%sat_lon                     &
     & ,                   SW_Control(j_sw)%i_cloud                     &
     & ,                   SW_Control(j_sw)%i_angular_integration       &
     & ,                   SW_Control(j_sw)%i_sph_mode                  &
     & ,                   row_length*rows, seg_points_temp(i)          &
     & ,                   model_levels, cloud_levels                   &
     & ,                   nd_cloud_type, nd_cloud_component            &
     & ,                   SW_Control(j_sw)%l_extra_top                 &
     & ,                   id_ct, n_rad_layers, nd_column               &
     & ,                   SW_Spectrum(j_sw)%npd_band                   &
     & ,                   SW_Spectrum(j_sw)%n_band                     &
     & ,                   map_channel, n_channel, nd_channel           &
     & ,                   nd_viewing_level, n_viewing_level            &
     & ,                   viewing_level                                &
     & ,                   nd_direction, n_viewing_direction            &
     & ,                   viewing_direction                            &
     & ,                   nd_brdf_basis_fnc, nd_brdf_trunc             &
     & ,                   nd_profile, nd_flux_profile                  &
     & ,                   nd_radiance_profile                          &
     & ,                   nd_field_flux_diag, nd_field_rad_diag        &
     & ,                   L_MOSES_II, L_CTILE, nd_tile                 &
     & ,                   nd_point_tile)

! Heating rates are converted to increments across
! a timestep to advance the model, but are more
! conveniently processed as heating rates on other
! calls.
                        l_scale_inc = (j_sw == 1)
!
#if defined(A01_3Z)
! DEPENDS ON: r2_swrad3z
                        Call r2_swrad3z(error_code,                     &
! Mixing Ratios
     &   q_n(first_point,1,1), J_CO2_MMR,j_ozone(first_point,1,1),      &
     &   J_O2_MMR, CO2_DIM1, CO2_DIM2,                                  &
     &   CO2_3D(first_point,1,1), L_CO2_3D,                             &
! Pressures and Temperatures
     &   p_star(first_point,1), p_layer_boundaries(first_point,1,0),    &
     &   p_layer_centres(first_point,1,0),T_n(first_point,1,1),         &
! Options for treating clouds
     &   global_cloud_top, L_microphy,                                  &
     &   L_INHOM_CLOUD, INHOM_CLOUD_SW, DP_CORR_STRAT, DP_CORR_CONV,    &
! Stratiform Cloud Fields
     &   L_cloud_water_partition,                                       &
     &   area_cloud_fraction(first_point,1,1), cf_n(first_point,1,1),   &
     &   qcl_n(first_point,1,1), qcf_n(first_point,1,1),                &
! Convective Cloud Fields
     &   cca(first_point,1,1), cclwp(first_point,1),                    &
     &   ccb(first_point,1), cct(first_point,1),                        &
! Surface Fields
     &   land_albedo(first_point,1,1), L_MOSES_II, L_CTILE,             &
     &   land_alb(first_point,1),                                       &
     &   sice_alb(first_point,1), flandg(first_point,1),                &
     &   open_sea_albedo(first_point,1,1,j_sw),                         &
     &   ice_fract(first_point,1), land_sea_mask(first_point,1),        &
     &   land0p5(first_point,1), snow_depth(first_point,1),             &
! Solar Fields
     &   cos_zenith_angle_2(first_point,1,j_sw),                        &
     &   day_fraction(first_point,1,j_sw),                              &
     &   list_daylight_points(start_point,j_sw), SCS,                   &
! Aerosol Fields
     &   l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero,        &
     &   zh, aero_bl_levels,                                            &
     &   L_use_sulpc_direct, L_use_sulpc_indirect_SW,                   &
     &   sulp_dim1, sulp_dim2,                                          &
     &   accum_sulphate(first_point_sulpc, 1),                          &
     &   aitken_sulphate(first_point_sulpc, 1),                         &
     &   diss_sulphate(first_point_sulpc, 1),                           &
     &   sea_salt_film(first_point_seasalt,1,1),                        &
     &   sea_salt_jet(first_point_seasalt,1,1),                         &
     &   L_use_seasalt_indirect, L_use_seasalt_direct,                  &
     &   salt_dim1*salt_dim2, salt_dim3,                                &
! set l_use_soot_direct to false, dimensions to 1
     &   .false., 1, 1, fresh_soot, aged_soot,                          &
     &   L_USE_BIOGENIC, biogenic_dim1, biogenic_dim2,                  &
     &   local_biogenic(first_point_biogenic, 1),                       &
     &   l_use_arcl, arcl_dim1, arcl_dim2, n_arcl_species,              &
     &   n_arcl_compnts,i_arcl_compnts,local_arcl(first_point_arcl,1,1),&
     &   aerosol(first_point,1,1), L_MURK_RAD, Ntot_land, Ntot_sea,     &
! Level of tropopause
     &   trindx(first_point,1),                                         &
! Spectrum
     &   SW_Spectrum(j_sw),                                             &
! Algorithmic options
     &   SW_Control(j_sw), timestep, l_mod_k_flux, l_scale_inc, j_sw,   &
#if defined(RAD_DBG)
     &   i,                                                             &
#endif
! Satellite viewing geometry
     &n_viewing_direction, viewing_direction(first_point,1, 1),         &
     &viewing_direction(first_point,1,2),n_viewing_level,viewing_level, &
!
! All diagnostics and associated arrays
     & n_channel, map_channel, SW_diag(j_sw),                           &
     & diag_row_list(start_point,j_sw),diag_col_list(start_point,j_sw), &
!
! Dimensions
     &   lit_points,seg_points_temp(i),model_levels,n_rad_layers,       &
     &   cloud_levels,wet_model_levels,ozone_levels,row_length,         &
     &   rows,rows*row_length,nd_field_flux_diag,                       &
     &   nd_field_rad_diag,nd_profile,n_rad_layers,nd_column,           &
     &   n_cca_levels,nd_channel, nd_flux_profile,                      &
     &   nd_radiance_profile, nd_viewing_level, nd_direction,           &
     &   nd_cloud_component, nd_cloud_type,                             &
     &   nd_brdf_basis_fnc, nd_brdf_trunc,                              &
     &   nd_point_tile, nd_tile, id_ct,                                 &
! Output data
!
     &   surf_down_sw(first_point,1,1),                                 &
     &   flux_below_690nm_surf(first_point,1),                          &
     &   L_flux_below_690nm_surf,                                       &
     &   netsw(first_point,1),                                          &
     &   top_absorption(first_point,1),                                 &
     &   swsea(first_point,1),                                          &
     &   SW_incs(first_point,1,0),                                      &
!     &   dirpar_inc(first_point,1),
!
! Variables needed to calculate layer masses
     & rho_r2_nh(first_point,1),r_rho_levels_nh(first_point,1),         &
     &        r_theta_levels_nh(first_point,0),                         &
     &        q_n(first_point,1,1), qcl_n(first_point,1,1),             &
     &        qcf_n(first_point,1,1), qcf2_n(first_point,1,1),          &
     &        qrain_n(first_point,1,1), qgraup_n(first_point,1,1),      &
     &        l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio     &
!
     &        )
#endif

                        Deallocate(map_channel)

                      Else

#if defined(A01_3C)
! DEPENDS ON: r2_swrad3c
                        Call r2_swrad3c(error_code,                     &
! Mixing Ratios
     &   q_n(first_point,1,1), j_co2_mmr, j_ozone(first_point,1,1),     &
     &   j_o2_mmr, CO2_DIM1, CO2_DIM2, CO2_3D(first_point,1,1),         &
     &   L_CO2_3D, L_use_stochem_CH4, CH4_stochem(first_point,1,1),     &

! Pressures and Temperatures
     &   p_star(first_point,1), p_layer_boundaries(first_point,1,0),    &
     &   p_layer_centres(first_point,1,0), T_n(first_point,1,1),        &

! Options for treating clouds
     &   .true., global_cloud_top, L_microphy,                          &
     &   L_INHOM_CLOUD, INHOM_CLOUD_SW, DP_CORR_STRAT, DP_CORR_CONV,    &

! Stratiform Cloud Fields
     &   L_cloud_water_partition, L_PC2,                                &
     &   area_cloud_fraction(first_point,1,1), cf_n(first_point,1,1),   &
     &   qcl_n(first_point,1,1), qcf_n(first_point,1,1),                &

! Convective Cloud Fields
     &   cca(first_point,1,1), cclwp(first_point,1),                    &
     &   ccw(first_point,1,1), lcbase(first_point,1),                   &    
     &   ccb(first_point,1), cct(first_point,1),                        &

! Surface Fields
     &   land_albedo(first_point,1,1), L_MOSES_II, L_CTILE,             &
     &   L_USE_SPEC_SEA, land_alb(first_point,1),                       &
     &   sice_alb(first_point,1), flandg(first_point,1),                &
     &   open_sea_albedo(first_point,1,1,j_sw),                         &
     &   ice_fract(first_point,1), land_sea_mask(first_point,1),        &
     &   land0p5(first_point,1), snow_depth(first_point,1),             &

! Solar Fields
     &   cos_zenith_angle_2(first_point,1,j_sw),                        &
     &   day_fraction(first_point,1,j_sw),                              &
     &   list_daylight_points(start_point,j_sw), SCS,                   &

! Aerosol Fields
     & l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero, zh,      &
     & aero_bl_levels, L_use_clearrh, L_USE_DUST, DUST_DIM1, DUST_DIM2, &
     &   DUST_1(FIRST_POINT_DUST,1), DUST_2(FIRST_POINT_DUST,1),        &
     &   DUST_3(FIRST_POINT_DUST,1), DUST_4(FIRST_POINT_DUST,1),        &
     &   DUST_5(FIRST_POINT_DUST,1), DUST_6(FIRST_POINT_DUST,1),        &
     &   L_USE_BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2,                  &
     &   local_biogenic(first_point_biogenic, 1),                       &
     &   j_l_use_sulpc_direct, L_use_sulpc_indirect_SW, sulp_dim1,      &
     &   sulp_dim2, accum_sulphate(first_point_sulpc, 1),               &
     &   aitken_sulphate(first_point_sulpc, 1),                         &
     &   diss_sulphate(first_point_sulpc, 1),                           &
     &   L_volcts, volcmass(first_point,1),			        &
     &   sea_salt_film(first_point_seasalt,1,1),                        &
     &   sea_salt_jet(first_point_seasalt,1,1),                         &
     &   L_use_seasalt_indirect, j_l_use_seasalt_direct,                &
     &   salt_dim1*salt_dim2, salt_dim3, j_l_use_soot_direct,           &
     &   soot_dim1, soot_dim2, fresh_soot(first_point_soot, 1),         &
     &   aged_soot(first_point_soot, 1), j_l_use_bmass_direct,          &
     &   bmass_dim1, bmass_dim2, fresh_bmass(first_point_biomass, 1),   &
     &   aged_bmass(first_point_biomass, 1),                            &
     &   cloud_bmass(first_point_biomass, 1), L_use_bmass_indirect,     &
     &   j_l_use_ocff_direct, ocff_dim1, ocff_dim2,                     &
     & fresh_ocff(first_point_ocff, 1), aged_ocff(first_point_ocff, 1), &
     &   cloud_ocff(first_point_ocff, 1), L_use_ocff_indirect,          &
     &   L_USE_ARCL, arcl_dim1, arcl_dim2, n_arcl_species,              &
     &   n_arcl_compnts,i_arcl_compnts,local_arcl(first_point_arcl,1,1),&
     &   aerosol(first_point,1,1), L_MURK_RAD, Ntot_land, Ntot_sea,     &

! Level of tropopause
     &   trindx(first_point,1),                                         &

! Spectrum
     &   SW_Spectrum(j_sw),                                             &

! Algorithmic options
     &   SW_Control(j_sw), timestep, l_mod_k_flux, L_Wenyi,             &

! All diagnostics and associated arrays
     &   SW_diag(j_sw), diag_row_list(start_point,j_sw),                &
     &   diag_col_list(start_point,j_sw),                               &

! Dimensions
     &   lit_points, seg_points_temp(i), model_levels, n_rad_layers,    &
     &   cloud_levels, wet_model_levels, ozone_levels, row_length,      &
     &   rows, rows*row_length, nd_rad_pts, n_rad_layers, 1,            &
     &   n_cca_levels,                                                  &

! Output data
     &   surf_down_sw(first_point,1,1),                                 &
     &   flux_below_690nm_surf(first_point,1),                          &
     &   L_flux_below_690nm_surf,                                       &
     &   netsw(first_point,1),                                          &
     &   top_absorption(first_point,1),                                 &
     &   swsea(first_point,1),                                          &
     &   SW_incs(first_point,1,0),                                      &
     &   dirpar_inc(first_point,1),                                     &
!
! Variables needed to calculate layer masses
     &   rho_r2_nh(first_point,1),r_rho_levels_nh(first_point,1),       &
     &   r_theta_levels_nh(first_point,0),                              &
     &   q_n(first_point,1,1), qcl_n(first_point,1,1),                  &
     &   qcf_n(first_point,1,1), qcf2_n(first_point,1,1),               &
     &   qrain_n(first_point,1,1), qgraup_n(first_point,1,1),           &
     &   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio          &
!
     &        )
#endif

                      End If

                   enddo              ! Loop over segments
                endif                 ! If daylight_points>0

                If (L_orog) Deallocate(orog_corr)

#if !defined(SCMA)
!
! Radiative fluxes may not have been calculated at all
! points: we now fill in as required. In the case of
! spatial degradation calls to SWAPBOUNDS are made and
! these must be made on all PEs, which requires that
! this code shall executed even when the PE contains
! no daylit points.
!
! At the North Pole in a global domain calculations
! are performed only at the first point if lit, so
! the rest of the row must be filled in if this point
! is lit.
!
                L_complete_North=(model_domain==mt_global).AND.         &
     &                       ( at_extremity(PNorth) ) .AND.             &
     &                       ( switch(1,rows,j_sw) )
!
! At the South Pole in a global domain calculations
! are performed only at the first point if lit, so
! the rest of the row must be filled in if this point
! is lit.
!
                L_complete_South=(model_domain==mt_global).AND.         &
     &                       ( at_extremity(PSouth) ) .AND.             &
     &                       ( switch(1,1,j_sw) )
!
! When spatial degradation is performed fields must be
! filled in at alternate points.
!
                L_complete_deg = ( L_rad_deg ) .AND.                    &
     &                     ( tot_daylight_points > 0 )
!
! Set addressing limits for spatial degradation.
!
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
! Call appropriate subroutines to fill in missing data
! as required.
!
               If ( L_complete_North .OR. L_complete_South .OR.         &
     &                                   L_complete_deg ) Then

! DEPENDS ON: fill_missing_data_sw
                  call fill_missing_data_sw(                            &
     &               off_x, off_y, row_length, rows,model_levels,       &
     &               salt_dim1, salt_dim2, salt_dim3,                   &
     &               cloud_levels,first_row,last_row,                   &
     &               first_data_interp, ES_space_interp,                &
     &               L_complete_North, L_complete_South,                &
     &               L_complete_deg,L_flux_below_690nm_surf,            &
     &               n_channel, j_sw,                                   &
     &               SW_Control(1)%l_extra_top,                         &
     &               SW_incs,                                           &
     &               netsw,                                             &
     &               SWsea,                                             &
     &               flux_below_690nm_surf,                             &
     &               top_absorption,                                    &
     &               surf_down_sw,                                      &
     &               dirpar_inc,                                        &
     &               sea_salt_film, sea_salt_jet)
               endif
#endif
              if (l_timestep) then
                 sw_incs_local(:,:,:,j_sw)=sw_incs
                 netsw_local(:,:,j_sw)=netsw
                 swsea_local(:,:,j_sw)=swsea
                 surf_down_sw_local(:,:,:,j_sw)=surf_down_sw
                 flux_b690nm_local(:,:,j_sw)=flux_below_690nm_surf
                 if (SW_Control(1)%l_extra_top) then
                    top_abs_sw(:,:,j_sw)=top_absorption
                 endif
                 if (SW_diag(j_sw)%L_direct_par) then
                    dirpar_inc_local(:,:,j_sw)=dirpar_inc
                 endif
              endif

            endif                                   ! l_call_swrad
         enddo                                      ! radiation calls

         if (l_timestep) then
            if (l_rad_step_diag) then
               sw_incs(:,:,:) = sw_incs_local(:,:,:,1)                  &
     &                        + sw_incs_local(:,:,:,2)
               netSW(:,:)     = netSW_local(:,:,1)                      &
     &                        + netSW_local(:,:,2)
               SWsea(:,:)     = SWsea_local(:,:,1)                      &
     &                        + SWsea_local(:,:,2)
               surf_down_sw(:,:,:) = surf_down_sw_local(:,:,:,1)        &
     &                             + surf_down_sw_local(:,:,:,2)
               flux_below_690nm_surf(:,:) = flux_b690nm_local(:,:,1)    &
     &                                    + flux_b690nm_local(:,:,2)
               if (SW_Control(1)%l_extra_top) then
                  top_absorption(:,:) = top_abs_sw(:,:,1)               &
     &                                + top_abs_sw(:,:,2)
               endif
               if (SW_diag(1)%L_direct_par .or.                         &
     &             SW_diag(2)%L_direct_par ) then
                  dirpar_inc(:,:)     = dirpar_inc_local(:,:,1)         &
     &                                + dirpar_inc_local(:,:,2)
               endif
            endif
         endif

         if (l_timestep) then
            if (l_rad_step_diag) then
               Do j = 1, rows
                  Do i= 1, row_length
                     mean_cos_zenith_angle(i,j) =                       &
     &                   cos_zenith_angle_2(i,j,2)*day_fraction(i,j,2)
                  End Do
               End Do
            endif
         else
            if (l_rad_step_prog) then
               Do j = 1, rows
                  Do i= 1, row_length
                     mean_cos_zenith_angle(i,j) =                       &
     &                   cos_zenith_angle_2(i,j,1)*day_fraction(i,j,1)
                  End Do
               End Do
            endif
         endif
!
! calculate net surface SW for diagnostic
!
         If (L_CTILE) then
            Do j = 1, rows
               Do i= 1, row_length
                  surfsw(i,j) = SW_incs(i,j,0)                          &
     &                        * mean_cos_zenith_angle(i,j)              &
     &                        + (1.-flandg(i,j))*SWsea(i,j)
                  land_alb(i,j)=0.0
                  if(flandg(i,j) == 1.0)SWsea(i,j)=rmdi
               End Do
            End Do
         Else
            Do j = 1, rows
               Do i= 1, row_length
                  surfsw(i,j) = SW_incs(i,j,0)                          &
     &                        * mean_cos_zenith_angle(i,j)              &
     &                        + SWsea(i,j)
                  if(land_sea_mask(i,j))SWsea(i,j)=rmdi
               End Do
            End Do
         Endif
!
! calculate net surface SW on tiles if required for MOSES II
!
         If (L_MOSES_II) then

            Do l=1,land_field
               i = land_index_i(l)
               j = land_index_j(l)
               sw_down_rts(i,j) = 0.0
               sw_net_rts(i,j) = 0.0
            End Do
            Do n=1,ntiles
               Do l=1,land_field
                  SW_tile_rts(l,n) = 0.
               End Do
            End Do

            if (l_timestep) then
               if (l_rad_step_diag) then
                  Do n=1,ntiles
                     Do point=1,tile_pts(n)
                        l = tile_index(point,n)
                        i = land_index_i(l)
                        j = land_index_j(l)
                        Do k=1,4
                           SW_tile_rts(l,n)=SW_tile_rts(l,n)            &
     &                        + (1. - alb_tile(l,n,k,1))                &
     &                        * surf_down_sw_local(i,j,k,1)             &
     &                        + (1. - alb_tile(l,n,k,2))                &
     &                        * surf_down_sw_local(i,j,k,2)
                        End Do
                        sw_net_rts(i,j)=sw_net_rts(i,j) +               &
     &                            SW_tile_rts(l,n)*tile_frac(l,n)
                     End Do
                  End Do
               endif
            else
               if (l_rad_step_prog) then
                  Do n=1,ntiles
                     Do point=1,tile_pts(n)
                        l = tile_index(point,n)
                        i = land_index_i(l)
                        j = land_index_j(l)
                        Do k=1,4
                           SW_tile_rts(l,n) = SW_tile_rts(l,n)          &
     &                        + (1. - alb_tile(l,n,k,1))                &
     &                        *surf_down_sw(i,j,k)
                        End Do
                        sw_net_rts(i,j)=sw_net_rts(i,j) +               &
     &                           SW_tile_rts(l,n)*tile_frac(l,n)
                     End Do
                  End Do
               endif
            endif
            IF (L_CTILE) THEN
               Do l=1,land_field
                  i = land_index_i(l)
                  j = land_index_j(l)
                  Do k=1,4
                     sw_down_rts(i,j) = sw_down_rts(i,j) +              &
     &                       surf_down_sw(i,j,k)
                  Enddo
                  if (sw_down_rts(i,j) >  0.0)                          &
     &               land_alb(i,j)=1.0-sw_net_rts(i,j)                  &
     &                                     /sw_down_rts(i,j)
               End Do
            Endif
         End If
      End If ! on a radiation timestep

! Calculate day fraction and mean cos(solar zenith angle while
! the sun is up) for each grid point for this physics timestep:
! (if in fact full SW calculations are being done every timestep, this
! is of course unnecessary, as are various calculations later on)

!        If ( Radiation_tstep_prog  >   timestep ) Then

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
     &                day_fraction(1,1,1), cos_zenith_angle,            &
     &                hour_angle)

!        else

!           do j=1,rows
!              do i=1,row_length
!                 cos_zenith_angle(i,j)=cos_zenith_angle_2(i,j,1)
!              enddo
!           enddo

!        End If

! Combine the two terms to give the mean cos zenith angle over the
! whole of the physics timestep.
! calculate incoming SW at top of atmosphere

!        write(6,*) 'XXX glue_rad-rad_ctl3c (Glue_Rad): SC = ', SC

        Do j = 1, rows
          Do i= 1, row_length
            cos_zenith_angle(i,j) = cos_zenith_angle(i,j) *             &
     &                              day_fraction(i,j,1)
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
#if !defined(SCMA)
      If ( ( L_T_incr_sw ).or.                                          &
     &     ( sf(232,1) ) ) Then
#else
      If (L_T_incr_sw ) Then
#endif
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

#if defined(SCMA)
        If (L_SCMDiags(SCMDiag_pc2)) Then
!         Calculate equivalent to stash code 1,181 - the total
!         temperature increment including the pc2 scheme
          Do k = 1,model_levels
            Do j = 1,rows
              Do i = 1,row_length
                  work_3d(i,j,k) = T_latest(i,j,k) - T_n(i,j,k)
              Enddo ! i
            Enddo   ! j
          Enddo     ! k
          ! stash item = 1,181  ! temperature increment minus pc2
! DEPENDS ON: scmoutput
          call SCMoutput(work_3d,                                       &
               'sw1pc2','SW heating rate incl PC2','K/timestep',        &
               t_acc,d_all,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
          call SCMoutput(work_3d,                                       &
               'sw2pc2','SW heating rate incl PC2','K/day',             &
               t_mult,d_all,default_streams,'ntspday',RoutineName)
        Endif
#endif
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
              if (SW_diag(1)%L_direct_par) Then
                 SW_diag(1)%flxdirparsurf(i,j) = dirpar_inc(i,j)        &
     &                                  * cos_zenith_angle(i,j)
              endif
              if (SW_diag(2)%L_direct_par) Then
                 SW_diag(2)%flxdirparsurf(i,j) = dirpar_inc(i,j)        &
     &                                  * cos_zenith_angle(i,j)
              endif
            End Do
          End Do
        End If

! Set up net surface SW on tiles if required for MOSES II
          if (L_MOSES_II) then
            If (ntiles == 1) then
              Do l=1,land_field
                i = land_index_i(l)
                j = land_index_j(l)
                SW_tile(l,1) = SW_tile_rts(l,1)                         &
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

        If ((L_rad_step_prog.or.L_Rad_Step_diag) .and. L_emcorr) Then

! Sum short wave fluxes into the atmosphere and
! add into the net diabatic fluxes into the
! atmosphere for use in the energy correction
! procedure

          If (SW_Control(1)%l_extra_top) Then
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
     &                   sum_eng_fluxes,radiation_tstep_prog)

        End If

! DEPENDS ON: timer
        If (Ltimer) Call timer ('Eng Corr',4)

#if !defined(SCMA)
! ----------------------------------------------------------------------
! Section RAD.1.2 Short Wave Radiation diagnostics
!-----------------------------------------------------------------------

      do j_sw=1,n_swdiag

         l_call_swrad=.false.
!
! Check that sw diagnostics requested this timestep
!
         If (error_code  ==  0 .and. sf(0,1)) Then
! DEPENDS ON: timer
            If (Ltimer) Call timer ('Diags   ',3)

            if (l_timestep) then

               i_off=0

               if ((j_sw >  1).and.(l_rad_step_diag)) then

                  if (SW_diag(2)%L_solar_out_toa)                       &
     &                SW_diag(2)%solar_out_toa =                        &
     &                SW_diag(2)%solar_out_toa +                        &
     &                SW_diag(1)%solar_out_toa

                  if (SW_diag(2)%L_solar_out_clear)                     &
     &                SW_diag(2)%solar_out_clear =                      &
     &                SW_diag(2)%solar_out_clear +                      &
     &                SW_diag(1)%solar_out_clear

                  if (SW_diag(2)%L_surface_down_flux)                   &
     &                SW_diag(2)%surface_down_flux =                    &
     &                SW_diag(2)%surface_down_flux +                    &
     &                SW_diag(1)%surface_down_flux

                  if (SW_diag(2)%L_surf_down_clr)                       &
     &                SW_diag(2)%surf_down_clr =                        &
     &                SW_diag(2)%surf_down_clr +                        &
     &                SW_diag(1)%surf_down_clr

                  if (SW_diag(2)%L_surf_up_clr)                         &
     &                SW_diag(2)%surf_up_clr =                          &
     &                SW_diag(2)%surf_up_clr +                          &
     &                SW_diag(1)%surf_up_clr

                  if (SW_diag(2)%L_clear_hr)                            &
     &                SW_diag(2)%clear_hr =                             &
     &                SW_diag(2)%clear_hr +                             &
     &                SW_diag(1)%clear_hr

                  if (SW_diag(2)%L_net_flux_trop)                       &
     &                SW_diag(2)%net_flux_trop =                        &
     &                SW_diag(2)%net_flux_trop +                        &
     &                SW_diag(1)%net_flux_trop

                  if (SW_diag(2)%L_up_flux_trop)                        &
     &                SW_diag(2)%up_flux_trop =                         &
     &                SW_diag(2)%up_flux_trop +                         &
     &                SW_diag(1)%up_flux_trop

                  if (SW_diag(2)%L_flux_direct)                         &
     &                SW_diag(2)%flux_direct =                          &
     &                SW_diag(2)%flux_direct +                          &
     &                SW_diag(1)%flux_direct

                  if (SW_diag(2)%L_flux_diffuse)                        &
     &                SW_diag(2)%flux_diffuse =                         &
     &                SW_diag(2)%flux_diffuse +                         &
     &                SW_diag(1)%flux_diffuse

                  if (SW_diag(2)%L_uvflux_direct)                       &
     &                SW_diag(2)%uvflux_direct =                        &
     &                SW_diag(2)%uvflux_direct +                        &
     &                SW_diag(1)%uvflux_direct

                  if (SW_diag(2)%L_uvflux_up)                           &
     &                SW_diag(2)%uvflux_up =                            &
     &                SW_diag(2)%uvflux_up +                            &
     &                SW_diag(1)%uvflux_up

                  if (SW_diag(2)%L_uvflux_net)                          &
     &                SW_diag(2)%uvflux_net =                           &
     &                SW_diag(2)%uvflux_net +                           &
     &                SW_diag(1)%uvflux_net

                  if (SW_diag(2)%L_FlxSolBelow690nmSurf)                &
     &                SW_diag(2)%FlxSolBelow690nmSurf =                 &
     &                SW_diag(2)%FlxSolBelow690nmSurf +                 &
     &                SW_diag(1)%FlxSolBelow690nmSurf

                  if (SW_diag(2)%L_FlxSeaBelow690nmSurf)                &
     &                SW_diag(2)%FlxSeaBelow690nmSurf =                 &
     &                SW_diag(2)%FlxSeaBelow690nmSurf +                 &
     &                SW_diag(1)%FlxSeaBelow690nmSurf

                  if (SW_diag(2)%L_cloud_extinction)                    &
     &                SW_diag(2)%cloud_extinction =                     &
     &                SW_diag(2)%cloud_extinction +                     &
     &                SW_diag(1)%cloud_extinction

                  if (SW_diag(2)%L_cloud_weight_extinction)             &
     &                SW_diag(2)%cloud_weight_extinction =              &
     &                SW_diag(2)%cloud_weight_extinction +              &
     &                SW_diag(1)%cloud_weight_extinction

                  if (SW_diag(2)%L_ls_cloud_extinction)                 &
     &                SW_diag(2)%ls_cloud_extinction =                  &
     &                SW_diag(2)%ls_cloud_extinction +                  &
     &                SW_diag(1)%ls_cloud_extinction

                  if (SW_diag(2)%L_ls_cloud_weight_extinction)          &
     &                SW_diag(2)%ls_cloud_weight_extinction =           &
     &                SW_diag(2)%ls_cloud_weight_extinction +           &
     &                SW_diag(1)%ls_cloud_weight_extinction

                  if (SW_diag(2)%L_cnv_cloud_extinction)                &
     &                SW_diag(2)%cnv_cloud_extinction =                 &
     &                SW_diag(2)%cnv_cloud_extinction +                 &
     &                SW_diag(1)%cnv_cloud_extinction

                  if (SW_diag(2)%L_cnv_cloud_weight_extinction)         &
     &                SW_diag(2)%cnv_cloud_weight_extinction =          &
     &                SW_diag(2)%cnv_cloud_weight_extinction +          &
     &                SW_diag(1)%cnv_cloud_weight_extinction

                  l_call_swrad=.true.

               endif


            elseif (l_forcing.and.l_rad_step_diag) then

! Calculate SW forcings. Note that forcings are only calculated for
! fluxes. Difference between j_sw=1 (call with forcing in place)
! and j_sw= 2 (call with reference values). N.B. Call to radiation
! schemes done first with j_sw=2 then with j_sw=1

               i_off=diagnostic_offset*(j_sw-1)

               if (j_sw == n_swdiag) then

                  if (SW_diag(j_sw)%L_solar_out_toa)                    &
     &                SW_diag(j_sw)%solar_out_toa =                     &
     &                SW_diag(1)%solar_out_toa -                        &
     &                SW_diag(j_sw)%solar_out_toa

                  if (SW_diag(j_sw)%L_solar_out_clear)                  &
     &                SW_diag(j_sw)%solar_out_clear =                   &
     &                SW_diag(1)%solar_out_clear -                      &
     &                SW_diag(j_sw)%solar_out_clear

                  if (SW_diag(j_sw)%L_surface_down_flux)                &
     &                SW_diag(j_sw)%surface_down_flux =                 &
     &                SW_diag(1)%surface_down_flux -                    &
     &                SW_diag(j_sw)%surface_down_flux

                  if (SW_diag(j_sw)%L_surf_down_clr)                    &
     &                SW_diag(j_sw)%surf_down_clr =                     &
     &                SW_diag(1)%surf_down_clr -                        &
     &                SW_diag(j_sw)%surf_down_clr

                  if (SW_diag(j_sw)%L_surf_up_clr)                      &
     &                SW_diag(j_sw)%surf_up_clr =                       &
     &                SW_diag(1)%surf_up_clr -                          &
     &                SW_diag(j_sw)%surf_up_clr

                  if (SW_diag(j_sw)%L_clear_hr)                         &
     &                SW_diag(j_sw)%clear_hr =                          &
     &                SW_diag(1)%clear_hr -                             &
     &                SW_diag(j_sw)%clear_hr

                  if (SW_diag(j_sw)%L_net_flux_trop)                    &
     &                SW_diag(j_sw)%net_flux_trop =                     &
     &                SW_diag(1)%net_flux_trop -                        &
     &                SW_diag(j_sw)%net_flux_trop

                  if (SW_diag(j_sw)%L_up_flux_trop)                     &
     &                SW_diag(j_sw)%up_flux_trop =                      &
     &                SW_diag(1)%up_flux_trop -                         &
     &                SW_diag(j_sw)%up_flux_trop

                  if (SW_diag(j_sw)%L_flux_direct)                      &
     &                SW_diag(j_sw)%flux_direct =                       &
     &                SW_diag(1)%flux_direct -                          &
     &                SW_diag(j_sw)%flux_direct

                  if (SW_diag(j_sw)%L_flux_diffuse)                     &
     &                SW_diag(j_sw)%flux_diffuse =                      &
     &                SW_diag(1)%flux_diffuse -                         &
     &                SW_diag(j_sw)%flux_diffuse

                  if (SW_diag(j_sw)%L_uvflux_direct)                    &
     &                SW_diag(j_sw)%uvflux_direct =                     &
     &                SW_diag(1)%uvflux_direct -                        &
     &                SW_diag(j_sw)%uvflux_direct

                  if (SW_diag(j_sw)%L_uvflux_up)                        &
     &                SW_diag(j_sw)%uvflux_up =                         &
     &                SW_diag(1)%uvflux_up -                            &
     &                SW_diag(j_sw)%uvflux_up

                  if (SW_diag(j_sw)%L_uvflux_net)                       &
     &                SW_diag(j_sw)%uvflux_net =                        &
     &                SW_diag(1)%uvflux_net -                           &
     &                SW_diag(j_sw)%uvflux_net

               endif
               l_call_swrad=.true.
            elseif (l_radiance) then
               if (j_sw == 1) then
                  i_off = 0
               else
                  i_off=90+(diagnostic_offset/20)*(j_sw-1)
               endif
               l_call_swrad=.true.
            else
               i_off=0
               if (l_rad_step_prog) then
                  l_call_swrad=.true.
               endif
            endif

! DEPENDS ON: diagnostics_sw
            Call diagnostics_sw(                                        &
     &                     row_length, rows, model_levels               &
     &,                    wet_model_levels, cloud_levels               &
     &,                    n_rows, global_row_length, global_rows       &
     &,                    halo_i, halo_j, off_x, off_y, me             &
     &,                    n_proc, n_procx, n_procy                     &
     &,                    g_rows, g_row_length, g_datastart            &
     &,                    at_extremity                                 &
     &,                    timestep,i_off                               &
     &,                    T_n, T_inc                                   &
     &,                    q_n, qcl_n, cf_n, cfl_n                      &
     &,                    T_latest, q_latest, qcl_latest               &
     &,                    cf_latest, cfl_latest                        &
     &,                    surfsw, itoasw                               &
     &,                    SWsea, flux_below_690nm_surf                 &
     &,                    photosynth_act_rad                           &
     &,                    f_orog, slope_aspect, slope_angle            &
     &,                    sw_net_land,sw_net_sice                      &
     &,                    T_incr_diagnostic                            &
     &,                    SW_spectrum(1)%n_band                        &
     &,                    sea_salt_film, sea_salt_jet                  &
     &,                    salt_dim1, salt_dim2, salt_dim3              &
     &,                    l_call_swrad, j_sw                           &
     &,                                                                 &
#include "argsts.h"
     & STASHwork1)

! DEPENDS ON: timer
           If (Ltimer) Call timer ('Diags   ',4)
         endif
      enddo
#else
      If (L_SCMDiags(SCMDiag_rad)) Then
!       ! Output some SCM diagnostics for SW radiation
!       ! Equiv to stash 1,161
! DEPENDS ON: scmoutput
        call SCMoutput(T_incr_diagnostic,                               &
             'sw1','SW heating rate','K/timestep',                      &
             t_acc,d_all,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
        call SCMoutput(T_incr_diagnostic,                               &
             'sw2','SW heating rate','K/day',                           &
             t_mult,d_all,default_streams,'ntspday',RoutineName)
! DEPENDS ON: scmoutput
        call SCMoutput(surfsw,                                          &
             'surfsw','Net surface SW flux','W/m2',                     &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
        call SCMoutput(itoasw,                                          &
             'is_toa','Incoming solar radn (TOA)','W/m2',               &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
        call SCMoutput(flux_below_690nm_surf,                           &
             'surf_sw_b1','Net SW surface flux in band 1','W/m2',       &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 1,208
! DEPENDS ON: scmoutput
        call SCMoutput(SW_diag(j_sw)%solar_out_toa,                     &
             'os_toa','Outgoing solar radn (TOA)','W/m2',               &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 1,209
! DEPENDS ON: scmoutput
        call SCMoutput(SW_diag(j_sw)%solar_out_clear,                   &
             'cs_os','Clear-sky outgoing SW','W/m2',                    &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 1,210
! DEPENDS ON: scmoutput
        call SCMoutput(SW_diag(j_sw)%surf_down_clr,                     &
             'cs_surf_dnsw','Clear-sky down SW flux','W/m2',            &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 1,211
! DEPENDS ON: scmoutput
        call SCMoutput(SW_diag(j_sw)%surf_up_clr,                       &
             'cs_surf_upsw','Clear-sky up SW flux','W/m2',              &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 1,221
! DEPENDS ON: scmoutput
        call SCMoutput(SW_diag(j_sw)%re_strat,                          &
             're_strat',                                                &
             'Layer cld liq effective radius * layer cld weight','-',   &
             t_avg+only_radsteps,d_cloud,default_streams,'',            &
             RoutineName)
!       ! stash 1,223
! DEPENDS ON: scmoutput
        call SCMoutput(SW_diag(j_sw)%wgt_strat,                         &
             'wgt_strat','Layer cloud weight for microphysics','-',     &
             t_avg+only_radsteps,d_cloud,default_streams,'',            &
             RoutineName)
!       ! stash 1,224
! DEPENDS ON: scmoutput
        call SCMoutput(SW_diag(j_sw)%lwp_strat,                         &
             'lwp_strat','Layer cld liquid water path * weight',        &
             'kg/m2',t_avg+only_radsteps,d_cloud,default_streams,'',    &
              RoutineName)
!       ! stash 1,225
! DEPENDS ON: scmoutput
        call SCMoutput(SW_diag(j_sw)%re_conv,                           &
             're_conv','Conv cloud liq re * conv cld weight','-',       &
             t_avg+only_radsteps,d_cloud,default_streams,'',            &
             RoutineName)
!       ! stash 1,226
! DEPENDS ON: scmoutput
        call SCMoutput(SW_diag(j_sw)%wgt_conv,                          &
             'wgt_conv','Conv cloud weight for microphysics','-',       &
             t_avg+only_radsteps,d_cloud,default_streams,'',            &
             RoutineName)
!   ! stash 1,182  Vapour increment
        Do k = 1,wet_model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3dw(i,j,k) = q_latest(i,j,k) - q_n(i,j,k)
            Enddo ! i
          Enddo   ! j
        Enddo     ! k
! DEPENDS ON: scmoutput
         call SCMoutput(work_3dw,                                       &
             'dq_sw','Specific humidity increment swrad','kg/kg',       &
             t_avg+only_radsteps,d_wet,default_streams,'',RoutineName)
!       ! stash 1,183  liquid water content increment
        Do k = 1,wet_model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3dw(i,j,k) = qcl_latest(i,j,k) - qcl_n(i,j,k)
            Enddo ! i
          Enddo   ! j
        Enddo     ! k
! DEPENDS ON: scmoutput
        call SCMoutput(work_3dw,                                        &
             'dqcl_sw','qcl increment swrad','kg/kg',                   &
             t_avg+only_radsteps,d_wet,default_streams,'',RoutineName)
!       ! stash 1,192  total cloud fraction increment
        Do k = 1,wet_model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3dw(i,j,k) = cf_latest(i,j,k) - cf_n(i,j,k)
            Enddo ! i
          Enddo   ! j
        Enddo     ! k
! DEPENDS ON: scmoutput
         call SCMoutput(work_3dw,                                       &
             'dbcf_sw',                                                 &
             'bulk cloud fraction increment swrad','fraction',          &
             t_avg+only_radsteps,d_wet,default_streams,'',RoutineName)
!       ! stash 1,193  liquid cloud fraction increment
        Do k = 1,wet_model_levels
          Do j = 1,rows
            Do i = 1,row_length
              work_3dw(i,j,k) = cfl_latest(i,j,k) - cfl_n(i,j,k)
            Enddo ! i
          Enddo   ! j
        Enddo     ! k
! DEPENDS ON: scmoutput
         call SCMoutput(work_3dw,                                       &
             'dcfl_sw',                                                 &
             'liquid cloud fraction increment swrad','fraction',        &
             t_avg+only_radsteps,d_wet,default_streams,'',RoutineName)
!       ! stash 1,232
!       ! sW heating = sw rad temp inc /timestep/timestep
        Do k = 1,model_levels
          Do j = 1,rows
            Do i = 1,row_length
               heating_rate(i,j,k) = T_incr_diagnostic(i,j,k) /         &
     &                                  timestep
            Enddo ! i
          Enddo   ! j
        Enddo     ! k
! DEPENDS ON: scmoutput
         call SCMoutput(heating_rate,                                   &
             'sw1rate','SW heating rate all timesteps','K/ts/ts',       &
             t_avg,d_all,default_streams,'',RoutineName)
!       ! stash 1,233
! DEPENDS ON: scmoutput
        call SCMoutput(SW_diag(j_sw)%clear_hr,                          &
             'dt_cssw','Clear-sky SW heating rates','K/s',              &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 1,235
! DEPENDS ON: scmoutput
        call SCMoutput(SW_diag(j_sw)%surface_down_flux,                 &
             'surfsw_dn','Total down SW flux','W/m2',                   &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!
      Endif
#endif

      Deallocate ( T_incr_diagnostic )
      If (L_orog) Then
        Deallocate(f_orog)
      Else
        Deallocate(f_orog,slope_aspect,slope_angle)
      Endif

      if (l_timestep) then
         if (l_rad_step_diag) then
! DEPENDS ON: deallocate_swdiag
            call deallocate_swdiag(2)
         endif
         if (timestep_number >  1) then
            if (mod(timestep_number,a_sw_radstep_prog) == 0) then

               deallocate(alb_tile)
               deallocate(sw_incs_local)
               deallocate(swsea_local)
               deallocate(netsw_local)
               deallocate(top_abs_sw)
               deallocate(flux_b690nm_local)
               deallocate(surf_down_sw_local)
               deallocate(dirpar_inc_local)

! DEPENDS ON: deallocate_swdiag
               call deallocate_swdiag(1)

            endif
         endif
      elseif (l_forcing) then
         if (l_rad_step_prog) then

            deallocate(alb_tile)

! DEPENDS ON: deallocate_swdiag
            call deallocate_swdiag(1)
         endif
         if (l_rad_step_diag) then
! DEPENDS ON: deallocate_swdiag
            call deallocate_swdiag(2)
         endif
      elseif (l_radiance) then
         if (l_rad_step_prog) then

            deallocate(alb_tile)

            do j_sw=1,n_swcall
! DEPENDS ON: deallocate_swdiag
               call deallocate_swdiag(j_sw)
            enddo
         endif
      else
         if (l_rad_step_prog) then

            deallocate(alb_tile)

! DEPENDS ON: deallocate_swdiag
            call deallocate_swdiag(1)
         endif
      endif

! ----------------------------------------------------------------------
! Section RAD.2 Long Wave Radiation Code.
!-----------------------------------------------------------------------

      if (l_rad_step_prog) then

! DEPENDS ON: allocate_lw_diag
        call allocate_lw_diag(row_length, rows, model_levels,           &
     &                 cloud_levels, 1 )

        if (l_timestep) then
          allocate(LW_incs_local(row_length, rows, 0:model_levels,      &
     &                                                    n_lwcall))
          allocate(LWsea_local(row_length, rows, n_lwcall))
          allocate(OLR_local(row_length, rows, n_lwcall))
          allocate(top_abs_lw(row_length, rows, n_lwcall))
        endif

         if (l_radiance) then
           do j_lw=2,n_lwcall

! DEPENDS ON: allocate_lw_diag
             call allocate_lw_diag( row_length, rows,                   &
     &                model_levels, cloud_levels, j_lw )
           enddo
         endif

      endif

      if (l_rad_step_diag) then
        if (l_timestep.or.l_forcing) then

! DEPENDS ON: allocate_lw_diag
            call allocate_lw_diag(row_length, rows, model_levels,       &
     &                    cloud_levels, 2 )
        endif
      endif

! DEPENDS ON: timer
        If (Ltimer) Call timer ('LW Rad  ',3)

        If (L_rad_step_prog.or.L_rad_step_diag) Then

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

           do j_lw=n_lwcall,1,-1

! For each LW call set the correct mass mixing ratios and switch the
! aerosols on/off, &c.  Note that as for the SW, in the case of the
! aerosols the reference  state is 'no aerosols' and hence they are
! switched off in the diagnostic call.
!
! As in the SW, since the "default" for each variable is to have it
! set in the diagnostic call to the same as in the prognostic, we set
! them all unconditionally to the prognostic & then re-set as needed.


!
! Set default values for greenhouse gases
!
              j_co2_mmr    = CO2_MMR
              j_n2o_mmr    = n2o_mix_ratio
              j_ch4_mmr    = ch4_mix_ratio
              j_cfc11_mmr  = cfc11_mix_ratio
              j_cfc12_mmr  = cfc12_mix_ratio
              j_c113_mmr   = C113MMR
              j_hcfc22_mmr = HCFC22MMR
              j_hfc125_mmr = HFC125MMR
              j_hfc134_mmr = HFC134AMMR
              j_ozone      = ozone
!
! Set default flags for aerosols
!
              j_l_use_sulpc_direct   = L_use_sulpc_direct
              j_l_use_seasalt_direct = L_use_seasalt_direct
              j_l_use_soot_direct    = L_use_soot_direct
              j_l_use_bmass_direct   = L_use_bmass_direct
              j_l_use_ocff_direct    = L_use_ocff_direct

              l_call_lwrad=.false.
!
! if timestepping required
!
              if (l_timestep) then

                 if ((j_lw==1).and.(l_rad_step_prog)) then
                    a_lw_radstep=a_lw_radstep_prog
                    l_call_lwrad=.true.
                endif
                if ((j_lw==2).and.(l_rad_step_diag)) then
                   a_lw_radstep=a_lw_radstep_diag
                   l_call_lwrad=.true.
                endif
! If Forcing required
!
              elseif (l_forcing) then
!
! Decide whether it is necessary to call radiation
!
                if ((j_lw==1).and.(l_rad_step_prog)) then
                    a_lw_radstep=a_lw_radstep_prog
                    l_call_lwrad=.true.
                endif
                if ((j_lw==2).and.(l_rad_step_diag)) then
                   a_lw_radstep=a_lw_radstep_prog
                   l_call_lwrad=.true.
                endif
!
                IF ( J_LW  ==  1 ) THEN        ! Prognostic Call


                ELSE                ! J_LW > 1: Diagnostic Call
!
! Greenhouse Gases
!
                    IF ( C2C_CO2(J_LW-1) ) j_co2_mmr = CO2_MMR_D
                    IF ( C2C_N2O(J_LW-1) ) j_n2o_mmr = n2o_mix_ratio_d
                    IF ( C2C_CH4(J_LW-1) ) j_ch4_mmr = ch4_mix_ratio_d
                    IF ( C2C_CFC11(J_LW-1)) j_cfc11_mmr=                &
     &                                               cfc11_mix_ratio_d
                    IF ( C2C_CFC12(J_LW-1)) j_cfc12_mmr=                &
     &                                               cfc12_mix_ratio_d
                    IF ( C2C_C113(J_LW-1) ) j_c113_mmr = C113MMR_D
                    IF ( C2C_HCFC22(J_LW-1) ) j_hcfc22_mmr = HCFC22MMR_D
                    IF ( C2C_HFC125(J_LW-1) ) j_hfc125_mmr = HFC125MMR_D
                    IF ( C2C_HFC134(J_LW-1) ) j_hfc134_mmr =            &
     &                                                      HFC134AMMR_D
!
! Ozone: Note that aerosol contains the MURK array.  It is in fact what
! is called d1(MURK) in Simon Tett's old vn4.5 mod.  Care needs to be
! taken so that MURK has the correct dimensions before it arrives in
! RAD_CTL2.
!
                    IF ( C2C_O3(J_LW-1) ) j_ozone = aerosol
!
! Aerosols
!
                    IF ( C2C_SULPC_D(J_LW-1) )                          &
     &                 j_l_use_sulpc_direct   = .false.
                    IF ( C2C_SEAS_D(J_LW-1) )                           &
     &                 j_l_use_seasalt_direct = .false.
                    IF ( C2C_SOOT_D(J_LW-1) )                           &
     &                 j_l_use_soot_direct    = .false.
                    IF ( C2C_BMB_D(J_LW-1) )                            &
     &                 j_l_use_bmass_direct   = .false.
                    IF ( C2C_OCFF_D(J_LW-1) )                           &
     &                 j_l_use_ocff_direct    = .false.

                 endif
              elseif (l_radiance) then
                 a_lw_radstep=a_lw_radstep_prog
                 l_call_lwrad=.true.
              else
                 a_lw_radstep=a_lw_radstep_prog
                 l_call_lwrad=.true.
              endif

              if (l_call_lwrad) then

! DEPENDS ON: prelim_lwrad
                 call prelim_lwrad(                                     &
! Parallel variables
     &                at_extremity, n_proc                              &
! Model Dimensions
     &,               row_length,rows,model_levels                      &
! Model Switches
     &,               model_domain,mt_global                            &
     &,               L_rad_deg, LW_Control(1)%l_extra_top              &
     &,               LW_Control(j_lw)%L_subsample                      &
     &,               LW_Control(j_lw)%L_geostationary                  &
! Time stepping Information
     &,               timestep_number,a_lw_radstep                      &
! ancillary fields and fields needed to be kept from timestep to
! timestep
     &,               lw_incs                                           &
! Satellite Geometry
     &,               LW_Control(j_lw)%min_view_lon                     &
     &,               LW_Control(j_lw)%max_view_lon                     &
     &,               LW_Control(j_lw)%min_view_lat                     &
     &,               LW_Control(j_lw)%max_view_lat                     &
! Number of Call
     &,               j_lw                                              &
! Other variables
     &,               sin_true_latitude,true_longitude                  &
     &,               seconds_since_midnight                            &
     &,               rad_mask, list_lw_points, first_data_interp       &
     &,               first_row,last_row,diag_row_list(1,j_lw)          &
     &,               diag_col_list(1,j_lw)                             &
     &,               lw_points                                         &
     &,               OLR, Lwsea, top_absorption )

!
! Now, segment all the LW calculation points
!
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
                       num_lw_points=lw_points-step                     &
     &                                         *(a_lw_segments-1)
                       seg_points = row_length*rows-seg_start+1
                    End If
                    Do j = first_point, first_point+num_lw_points-1
                       LIST_LW_points(j)= LIST_LW_points(j)             &
     &                                              -seg_start+1
                    End Do
                    first_point_local(i) = seg_start
                    seg_points_local(i) = num_lw_points
                    seg_start = seg_start+seg_points
                 End Do


                 Do i = 1, A_LW_segments

                    first_point = 1 + (i-1) * step
!
! Set the pointer to the beginning of the current segment.
! This is done solely to reduce the number of continuation
! lines required by the call.
!
                    ptr_local=first_point_local(i)

! Set the actual size of arrays in the radiation code:
! for some architectures (such as that of Cray vector
! machines) on odd size is preferred to avoid memory
! bank conflicts.
!
                    nd_rad_pts=2*(seg_points_local(i)/2)+1
!
! Set the number of layers seen in the radiation code.
! This may optionally be 1 greater than the number used
! in the rest of the model to avoid spurious effects
! resulting from the upper boundary (principally in
! stratospheric configurations).
!
                    IF (LW_Control(1)%l_extra_top) Then
                       n_rad_layers=model_levels+1
                    else
                       n_rad_layers=model_levels
                    endif
!
! Set the first point of the dust arrays to be used.
!
                    IF (L_USE_DUST) THEN
                       FIRST_POINT_DUST=FIRST_POINT_LOCAL(i)
                    ELSE
                       FIRST_POINT_DUST=1
                    ENDIF

! Set the first point of the biogenic array
                    IF (L_USE_BIOGENIC) THEN
                       FIRST_POINT_BIOGENIC=FIRST_POINT_LOCAL(i)
                    ELSE
                       FIRST_POINT_BIOGENIC=1
                    ENDIF
!
! Set the first points of the arrays of sulphates to be used.
!
! A separate assignment is necessary since
! not be of the full size unless the sulphur cycle is on.
!
                    IF (L_USE_SULPC_DIRECT .OR.                         &
     &                        L_USE_SULPC_INDIRECT_LW) THEN
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
                    If (L_use_bmass_direct .or.                         &
     &                          L_use_bmass_indirect) then
                       first_point_biomass=first_point_local(i)
                    Else
                       first_point_biomass=1
                    Endif
!
                    If (L_use_ocff_direct .or.                         &
     &                  L_use_ocff_indirect) then
                       first_point_ocff=first_point_local(i)
                    Else
                       first_point_ocff=1
                    Endif
!
                    If (L_use_seasalt_indirect .OR.                     &
     &                            L_use_seasalt_direct) Then
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

                    If (h_sect(2) == '03Z') Then

! Set dynamically determined array sizes required for
! radiation calculations. There is redundancy here in
! that the viewing directions are calculated for the
! whole PE in each segemnt: ideally, this will be
! tidied.
                      Allocate(map_channel(LW_Spectrum(j_lw)%n_band))

! DEPENDS ON: r2_set_rad_dim
                      Call r2_set_rad_dim(row_length, rows              &
     & ,               sin_true_latitude, true_longitude                &
     & ,               .FALSE., sindec, hour_angle                      &
     & ,               LW_Control(j_lw)%l_geostationary                 &
     & ,               LW_Control(j_lw)%sat_hgt                         &
     & ,               LW_Control(j_lw)%sat_lat                         &
     & ,               LW_Control(j_lw)%sat_lon                         &
     & ,               LW_Control(j_lw)%i_cloud                         &
     & ,               LW_Control(j_lw)%i_angular_integration           &
     & ,               LW_Control(j_lw)%i_sph_mode                      &
     & ,               row_length*rows, seg_points_local(i)             &
     & ,               model_levels, cloud_levels                       &
     & ,               nd_cloud_type, nd_cloud_component                &
     & ,               LW_control(j_lw)%l_extra_top                     &
     & ,               id_ct, n_rad_layers, nd_column                   &
     & ,               LW_Spectrum(j_lw)%npd_band                       &
     & ,               LW_Spectrum(j_lw)%n_band                         &
     & ,               map_channel, n_channel, nd_channel               &
     & ,               nd_viewing_level, n_viewing_level                &
     & ,               viewing_level                                    &
     & ,               nd_direction, n_viewing_direction                &
     & ,               viewing_direction                                &
     & ,               nd_brdf_basis_fnc, nd_brdf_trunc                 &
     & ,               nd_profile, nd_flux_profile                      &
     & ,               nd_radiance_profile                              &
     & ,               nd_field_flux_diag, nd_field_rad_diag            &
     & ,               L_MOSES_II, L_CTILE, nd_tile                     &
     & ,               nd_point_tile)

! Heating rates are converted to increments across
! a timestep to advance the model, but are more
! conveniently processed as heating rates on other
! calls.
                      l_scale_inc = (j_lw == 1)

#if defined(A02_3Z)
! DEPENDS ON: r2_lwrad3z
                      Call r2_lwrad3z(error_code,                       &
! Input data
     &  q_n(ptr_local,1,1), J_CO2_MMR, j_ozone(ptr_local,1,1),          &
     &  1, 1, CO2_3D, L_CO2_3D,                                         &
     &  j_n2o_mmr, j_ch4_mmr, j_cfc11_mmr, j_cfc12_mmr,                 &
     &  J_C113_MMR, J_HCFC22_MMR, J_HFC125_MMR, J_HFC134_MMR,           &
     &  T_n(ptr_local,1,1),T_surf(ptr_local,1),                         &
     &  T_sol_rad(ptr_local,1),Tstar_sea(ptr_local,1),L_CTILE,          &
     &  p_star(ptr_local,1),p_layer_boundaries(ptr_local,1,0),          &
     &  p_layer_centres(ptr_local,1,0),                                 &
     &  height_theta(ptr_local,1,0), height_rho(ptr_local,1,1),         &
! Options for treating clouds
     &  global_cloud_top, L_INHOM_CLOUD, INHOM_CLOUD_LW,                &
     &  DP_CORR_STRAT, DP_CORR_CONV,                                    &
! Stratiform Cloud Fields
     &  L_cloud_water_partition,                                        &
     &  area_cloud_fraction(ptr_local,1,1), cf_n(ptr_local,1,1),        &
     &  qcl_n(ptr_local,1,1),qcf_n(ptr_local,1,1),                      &
! Convective Cloud Fields
     &  cca(ptr_local,1,1),cclwp(ptr_local,1),                          &
     &  ccb(ptr_local,1), cct(ptr_local,1),                             &
! Surface Fields. Only want the 0.5 threshold LAND mask 
! and fractional land:
     &  land0p5(ptr_local,1),flandg(ptr_local,1),                       &
     &  ice_fract(ptr_local,1),snow_depth(ptr_local,1),                 &
! Aerosol Fields
     &  l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero,         &
     &  zh, aero_bl_levels,                                             &
     &  L_use_sulpc_direct, L_use_sulpc_indirect_LW,                    &
     &  sulp_dim1, sulp_dim2,                                           &
     &  accum_sulphate(first_point_sulpc, 1),                           &
     &  aitken_sulphate(first_point_sulpc, 1),                          &
     &  diss_sulphate(first_point_sulpc, 1),                            &
     &  sea_salt_film(first_point_seasalt,1,1),                         &
     &  sea_salt_jet(first_point_seasalt,1,1),                          &
     &  L_use_seasalt_indirect, L_use_seasalt_direct,                   &
     &  salt_dim1*salt_dim2, salt_dim3,                                 &
! set l_use_soot_direct to false, dimensions to 1
     &  .false., 1, 1, fresh_soot, aged_soot,                           &
     &   L_USE_BIOGENIC, biogenic_dim1, biogenic_dim2,                  &
     &   local_biogenic(first_point_biogenic, 1),                       &
     &  L_USE_ARCL, arcl_dim1, arcl_dim2, n_arcl_species,               &
     &  n_arcl_compnts, i_arcl_compnts,local_arcl(first_point_arcl,1,1),&
     &  aerosol(first_point_local(i),1,1), L_MURK_RAD,                  &
     &  Ntot_land, Ntot_sea,                                            &
! Level of tropopause
     &   trindx(ptr_local, 1),                                          &
! Spectral data
     &   LW_Spectrum(j_lw),                                             &
!
! Algorithmic options
     &   LW_control(j_lw), timestep, l_mod_k_flux, l_scale_inc,         &
     &   LIST_LW_points(first_point), j_lw,                             &
#if defined(RAD_DBG)
     &   i,                                                             &
#endif
!
! Viewing geometry for satellite simulations
     &   n_viewing_direction,                                           &
     &   viewing_direction(ptr_local, 1, 1),                            &
     &   viewing_direction(ptr_local, 1, 2),                            &
     &   n_viewing_level, viewing_level,                                &

! All diagnostics
     &   n_channel, map_channel,                                        &
     &   LW_diag(j_lw), diag_row_list(first_point,j_lw),                &
     &   diag_col_list(first_point,j_lw),                               &

! Physical Dimensions
     &   seg_points_local(i),model_levels,n_rad_layers,                 &
     &   cloud_levels,wet_model_levels,ozone_levels,row_length,         &
     &   rows,row_length*rows, nd_field_flux_diag,                      &
     &   nd_field_rad_diag, nd_profile, n_rad_layers, nd_column,        &
     &   n_cca_levels,                                                  &
     &   nd_channel,                                                    &
     &   nd_flux_profile, nd_radiance_profile,                          &
     &   nd_viewing_level, nd_direction,                                &
     &   nd_cloud_component, nd_cloud_type,                             &
     &   nd_brdf_basis_fnc, nd_brdf_trunc,                              &
     &   nd_point_tile, nd_tile, id_ct,                                 &

! Output data
     &   OLR(ptr_local, 1), top_absorption(ptr_local, 1),               &
     &   LWsea(ptr_local, 1), LW_incs(ptr_local,1,0),                   &

! Data for calculation of layer masses
     & rho_r2_nh(ptr_local,1),                                          &
     & r_rho_levels_nh(ptr_local,1), r_theta_levels_nh(ptr_local,0),    &
     & q_n(ptr_local,1,1), qcl_n(ptr_local,1,1),qcf_n(ptr_local,1,1),   &
     & qcf2_n(ptr_local,1,1),qrain_n(ptr_local,1,1),                    &
     & qgraup_n(ptr_local,1,1),                                         &
     & l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio            &
     & )
#endif

                      Deallocate(map_channel)

                    Else

#if defined(A02_3C)
! DEPENDS ON: r2_lwrad3c
                      Call r2_lwrad3c(error_code,                       &
! Input data
     &   q_n(ptr_local,1,1), j_co2_mmr, j_ozone(ptr_local,1,1),         &
     &   CO2_DIM1, CO2_DIM2, CO2_3D(ptr_local,1,1), L_CO2_3D,           &
     &   L_use_stochem_CH4, CH4_stochem(ptr_local,1,1),                 &
! chemical greenhouse gas fields
     &   ngrgas, grgas_field(ptr_local,1,1,1),                          &
     &   j_n2o_mmr, j_ch4_mmr, j_cfc11_mmr, j_cfc12_mmr,                &
     &   j_c113_mmr, j_hcfc22_mmr, j_hfc125_mmr, j_hfc134_mmr,          &
     &   T_n(ptr_local,1,1), T_surf(ptr_local,1),                       &
     &   T_sol_rad(ptr_local,1), Tstar_sea(ptr_local,1), L_CTILE,       &
     &   p_star(ptr_local,1), p_layer_boundaries(ptr_local,1,0),        &
     &   p_layer_centres(ptr_local,1,0),                                &
     &   height_theta(ptr_local,1,0), height_rho(ptr_local,1,1),        &

! Options for treating clouds
     &   .true., global_cloud_top, L_INHOM_CLOUD, INHOM_CLOUD_LW,       &
     &   DP_CORR_STRAT, DP_CORR_CONV,                                   &

! Stratiform Cloud Fields
     &   L_cloud_water_partition, L_PC2,                                &
     &   area_cloud_fraction(ptr_local,1,1), cf_n(ptr_local,1,1),       &
     &   qcl_n(ptr_local,1,1), qcf_n(ptr_local,1,1),                    &

! Convective Cloud Fields
     &   cca(ptr_local,1,1), cclwp(ptr_local,1),                        &
     &   ccw(first_point,1,1), lcbase(first_point,1),                   &
     &   ccb(ptr_local,1), cct(ptr_local,1),                            &

! Surface Fields. Only want the 0.5 threshold LAND mask
! and fractional land:
     &   land0p5(ptr_local,1), flandg(ptr_local,1),                     &
     &   ice_fract(ptr_local,1), snow_depth(ptr_local,1),               &

! Aerosol Fields
     &   l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero,        &
     &   zh, aero_bl_levels,                                            &
     &   L_use_clearrh, L_USE_DUST, DUST_DIM1, DUST_DIM2,               &
     &   DUST_1(FIRST_POINT_DUST,1), DUST_2(FIRST_POINT_DUST,1),        &
     &   DUST_3(FIRST_POINT_DUST,1), DUST_4(FIRST_POINT_DUST,1),        &
     &   DUST_5(FIRST_POINT_DUST,1), DUST_6(FIRST_POINT_DUST,1),        &
     &   L_USE_BIOGENIC, biogenic_dim1, biogenic_dim2,                  &
     &   local_biogenic(first_point_biogenic, 1),                       &
     &   j_l_use_sulpc_direct, L_use_sulpc_indirect_LW, sulp_dim1,      &
     &   sulp_dim2, accum_sulphate(first_point_sulpc, 1),               &
     &   aitken_sulphate(first_point_sulpc, 1),                         &
     &   diss_sulphate(first_point_sulpc, 1),                           &
     &   L_volcts, volcmass(ptr_local, 1),   	                        &
     &   sea_salt_film(first_point_seasalt,1,1),                        &
     &   sea_salt_jet(first_point_seasalt,1,1),                         &
     &   L_use_seasalt_indirect, j_l_use_seasalt_direct,                &
     &   salt_dim1*salt_dim2, salt_dim3, j_l_use_soot_direct,           &
     &   soot_dim1, soot_dim2, fresh_soot(first_point_soot, 1),         &
     &   aged_soot(first_point_soot, 1), j_l_use_bmass_direct,          &
     &   bmass_dim1, bmass_dim2, fresh_bmass(first_point_biomass, 1),   &
     &   aged_bmass(first_point_biomass, 1),                            &
     &   cloud_bmass(first_point_biomass, 1), L_use_bmass_indirect,     &
     &   j_l_use_ocff_direct, ocff_dim1, ocff_dim2,                     &
     & fresh_ocff(first_point_ocff, 1), aged_ocff(first_point_ocff, 1), &
     &   cloud_ocff(first_point_ocff, 1), L_use_ocff_indirect,          &
     &   L_USE_ARCL, arcl_dim1, arcl_dim2, n_arcl_species,              &
     &   n_arcl_compnts,i_arcl_compnts,local_arcl(first_point_arcl,1,1),&
     &   aerosol(first_point_local(i),1,1), L_MURK_RAD,                 &
     &   Ntot_land, Ntot_sea,                                           &

! Level of tropopause
     &   trindx(ptr_local, 1),                                          &

! Spectral Data
     &   LW_Spectrum(j_lw),                                             &

! Algotithmic Options
     &   LW_Control(j_lw),timestep, l_mod_k_flux, L_Wenyi,              &
     &   LIST_LW_points(first_point),                                   &

! All diagnostics
     &   LW_diag(j_lw), diag_row_list(first_point,j_lw),                &
     &   diag_col_list(first_point,j_lw),                               &

! Physical Dimensions
     &   seg_points_local(i), model_levels, n_rad_layers,               &
     &   cloud_levels, wet_model_levels, ozone_levels,                  &
     &   row_length, rows, row_length*rows, nd_rad_pts,                 &
     &   n_rad_layers, 1, n_cca_levels,                                 &

! Output data
     &   OLR(ptr_local, 1), top_absorption(ptr_local, 1),               &
     &   LWsea(ptr_local, 1), LW_incs(ptr_local,1,0),                   &
     &   rho_r2_nh(ptr_local,1),                                        &
     & r_rho_levels_nh(ptr_local,1), r_theta_levels_nh(ptr_local,0),    &
     & q_n(ptr_local,1,1), qcl_n(ptr_local,1,1),qcf_n(ptr_local,1,1),   &
     & qcf2_n(ptr_local,1,1),qrain_n(ptr_local,1,1),                    &
     & qgraup_n(ptr_local,1,1),                                         &
     & l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio            &
     & )
#endif

                    End If

                 End Do               ! end loop over long-wave segments

                 ! Deallocate the segmentation arrays
                 Deallocate( first_point_local )
                 Deallocate(  seg_points_local )


#if !defined(SCMA)

! Radiative fluxes may not have been calculated at all
! points. We now fill in as required.
!
! At the North Pole in a global domain calculations
! are performed only at the first point, so the rest of
! the row must be filled in.
!
                 L_complete_North=(model_domain==mt_global) .AND.       &
     &                       ( at_extremity(PNorth) )
!
! At the South Pole in a global domain calculations
! are performed only at the first point, so the rest of
! the row must be filled in.
!
                 L_complete_South=(model_domain==mt_global) .AND.       &
     &                       ( at_extremity(PSouth) )
!
! When spatial degradation is performed fields must be
! filled in at alternate points.
!
                 L_complete_deg = ( L_rad_deg )
!
!         Set addressing limits for spatial degradation.
!
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

! DEPENDS ON: fill_missing_data_lw
                 call fill_missing_data_lw(                             &
     &              off_x, off_y, row_length, rows, model_levels,       &
     &              cloud_levels, lw_spectrum(j_lw)%n_aod_wavel,        &
     &              first_row,last_row,                                 &
     &              first_data_interp, ES_space_interp,                 &
     &              L_complete_North, L_complete_South,                 &
     &              L_complete_deg, n_channel, j_lw,                    &
     &              LW_Control(1)%l_extra_top,                          &
     &              LW_incs, OLR, LWsea, top_absorption )

#endif
                 if (l_timestep) then
                    lw_incs_local(:,:,:,j_lw)=lw_incs
                    olr_local(:,:,j_lw)=olr
                    lwsea_local(:,:,j_lw)=lwsea
                    if (LW_Control(1)%l_extra_top) then
                       top_abs_lw(:,:,j_lw)=top_absorption
                    endif
                 endif

              endif                ! if radiation call required
           enddo                   ! loop over radiation calls

           if (l_timestep.and.l_rad_step_diag) then
              lw_incs(:,:,:) = lw_incs_local(:,:,:,1)                   &
     &                       + lw_incs_local(:,:,:,2)
              olr(:,:) = olr_local(:,:,1)                               &
     &                 + olr_local(:,:,2)
              lwsea(:,:) = lwsea_local(:,:,1)                           &
     &                   + lwsea_local(:,:,2)

              if (LW_Control(1)%l_extra_top) then
                 top_absorption(:,:) = top_abs_lw(:,:,1)                &
     &                               + top_abs_lw(:,:,2)
              endif
           endif
!
!
!
           If ( L_CTILE ) then
              Do j = 1, rows
                 Do i= 1, row_length
                    surflw(i,j) = LW_incs(i,j,0)                        &
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
!
! Store downward surface LW and dOLR on radiation timestep
! dOLR is TOA outward LW - surface upward LW for land and sea-ice
!
              Do j = 1, rows
                 Do i= 1, row_length
                    dOLR(i,j) = OLR(i,j)
                    If (.not.land_sea_mask(i,j).and.                    &
     &                              ice_fract(i,j) >  0.)then
                       dOLR(i,j) = dOLR(i,j) -                          &
     &                        ice_fract(i,j)*sbcon*tstar_sice(i,j)**4
                    End If
                 End Do
              End Do
              Do l=1,land_field
                 i = land_index_i(l)
                 j = land_index_j(l)
                 dOLR(i,j) = dOLR(i,j) -                                &
     &              fracsolid(i,j)*sbcon*T_sol_rad(i,j)**4
              End Do
              if (l_timestep) then
                 Do j = 1, rows
                    Do i= 1, row_length
                       LW_down(i,j) = LW_Diag(1)%surface_down_flux(i,j) &
     &                              + LW_Diag(2)%surface_down_flux(i,j)
                       dOLR_rts(i,j) = dOLR(i,j)
                    End Do
                 End Do
              else
                 Do j = 1, rows
                    Do i= 1, row_length
                       LW_down(i,j) = LW_Diag(1)%surface_down_flux(i,j)
                       dOLR_rts(i,j) = dOLR(i,j)
                    End Do
                 enddo
              endif
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
!
! DEPENDS ON: pc2_homog_plus_turb
          CALL PC2_HOMOG_PLUS_TURB(p_layer_centres(1,1,1),              &
     &      wet_model_levels,                                           &
     &      row_length, rows, timestep, T_latest, cf_latest, cfl_latest,&
     &      cff_latest, q_latest, qcl_latest, LW_incs(1,1,1),           &
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
#if defined(SCMA)
        If (L_SCMDiags(SCMDiag_pc2)) Then
!         Calculate equivalent to stash code 2,181 - the total
!         temperature increment including the pc2 scheme
          Do k = 1,model_levels
            Do j = 1,rows
              Do i = 1,row_length
                  work_3d(i,j,k) = T_latest(i,j,k) - T_n(i,j,k)
              Enddo ! i
            Enddo   ! j
          Enddo     ! k
          ! stash item = 2,181  ! temperature increment minus pc2
! DEPENDS ON: scmoutput
          call SCMoutput(work_3d,                                       &
               'lw1pc2','LW heating rate incl PC2','K/timestep',        &
               t_acc,d_all,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
          call SCMoutput(work_3d,                                       &
               'lw2pc2','LW heating rate incl PC2','K/day',             &
               t_mult,d_all,default_streams,'ntspday',RoutineName)
        Endif
#endif
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


      If ( L_Rad_Step_prog .AND. (LW_diag(1)%L_isccp_cf .OR.            &
     &               LW_diag(1)%L_isccp_cf_tau_0_to_p3 .OR.             &
     &               LW_diag(1)%L_isccp_cf_tau_p3_to_1p3 .OR.           &
     &               LW_diag(1)%L_isccp_cf_tau_1p3_to_3p6 .OR.          &
     &               LW_diag(1)%L_isccp_cf_tau_3p6_to_9p4 .OR.          &
     &               LW_diag(1)%L_isccp_cf_tau_9p4_to_23 .OR.           &
     &               LW_diag(1)%L_isccp_cf_tau_23_to_60 .OR.            &
     &               LW_diag(1)%L_isccp_cf_tau_ge_60) .AND.             &
     &               error_code  ==  0) Then


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
              If (SW_Control(1)%l_extra_top) Then
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
     &        day_fraction(first_point,1,1),                            &
     &        list_daylight_points(start_point,1), SCS,                 &
!                       Level of tropopause
     &        trindx(first_point,1),                                    &
!                       All diagnostics and associated arrays
     &        LW_diag(1), SW_diag(1),                                   &
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

#if !defined(SCMA)

!           Radiative fluxes may not have been calculated at all
!           points. We now fill in as required.
!
!           At the North Pole in a global domain calculations
!           are performed only at the first point if lit, so
!           the rest of the row must be filled in if this point
!           is lit.
            L_complete_North = ( model_domain == mt_global ) .AND.      &
     &                         ( at_extremity(PNorth) ) .AND.           &
     &                         ( switch(1,rows,1) )

!           At the South Pole in a global domain calculations
!           are performed only at the first point if lit, so
!           the rest of the row must be filled in if this point
!           is lit.
            L_complete_South = ( model_domain == mt_global ) .AND.      &
     &                         ( at_extremity(PSouth) ) .AND.           &
     &                         ( switch(1,1,1) )
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

! DEPENDS ON: fill_missing_data_isccp
          Call fill_missing_data_isccp(                                 &
     &      off_x, off_y, row_length, rows,                             &
     &      first_row,last_row,                                         &
     &      first_data_interp, 1, ES_space_interp,                      &
     &      L_complete_North, L_complete_South, L_complete_deg)
            Endif
#endif

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

        If ((L_Rad_step_prog.or.L_Rad_Step_diag) .and. L_emcorr) Then

! Sum long wave fluxes into the atmosphere and
! add into the net diabatic fluxes into the
! atmosphere for use in the energy correction
! procedure

        IF (LW_Control(1)%l_extra_top) Then
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
     &                    sum_eng_fluxes,radiation_tstep_prog)

        End If

! DEPENDS ON: timer
        If (Ltimer) Call timer ('Eng Corr',4)

#if !defined(SCMA)
! ----------------------------------------------------------------------
! Section RAD.2.3 Long Wave Radiation diagnostics
!-----------------------------------------------------------------------

      do j_lw=1,n_lwdiag

         l_call_lwrad=.false.
!
! Check that lw diagnostics requested this timestep
!
         If (error_code  ==  0 .and. sf(0,2)) Then
! DEPENDS ON: timer
            If (Ltimer) Call timer ('Diags   ',3)

            if (l_timestep) then

               i_off=0

               if ((j_lw == n_lwdiag).and.(l_rad_step_diag)) then

                  if (LW_diag(2)%L_clear_olr)                           &
     &               LW_diag(2)%clear_olr =                             &
     &               LW_diag(2)%clear_olr +                             &
     &               LW_diag(1)%clear_olr

                  if (LW_diag(2)%L_surface_down_flux)                   &
     &               LW_diag(2)%surface_down_flux =                     &
     &               LW_diag(2)%surface_down_flux +                     &
     &               LW_diag(1)%surface_down_flux

                  if (LW_diag(2)%L_surf_down_clr)                       &
     &               LW_diag(2)%surf_down_clr =                         &
     &               LW_diag(2)%surf_down_clr +                         &
     &               LW_diag(1)%surf_down_clr

                  if (LW_diag(2)%L_clear_hr)                            &
     &               LW_diag(2)%clear_hr =                              &
     &               LW_diag(2)%clear_hr +                              &
     &               LW_diag(1)%clear_hr

                  if (LW_diag(2)%L_net_flux_trop)                       &
     &               LW_diag(2)%net_flux_trop =                         &
     &               LW_diag(2)%net_flux_trop +                         &
     &               LW_diag(1)%net_flux_trop

                  if (LW_diag(2)%L_down_flux_trop)                      &
     &               LW_diag(2)%down_flux_trop =                        &
     &               LW_diag(2)%down_flux_trop +                        &
     &               LW_diag(1)%down_flux_trop

                  if (LW_diag(2)%L_cloud_absorptivity)                  &
     &                LW_diag(2)%cloud_absorptivity =                   &
     &                LW_diag(2)%cloud_absorptivity +                   &
     &                LW_diag(1)%cloud_absorptivity

                  if (LW_diag(2)%L_cloud_weight_absorptivity)           &
     &                LW_diag(2)%cloud_weight_absorptivity =            &
     &                LW_diag(2)%cloud_weight_absorptivity +            &
     &                LW_diag(1)%cloud_weight_absorptivity

                  if (LW_diag(2)%L_ls_cloud_absorptivity)               &
     &                LW_diag(2)%ls_cloud_absorptivity =                &
     &                LW_diag(2)%ls_cloud_absorptivity +                &
     &                LW_diag(1)%ls_cloud_absorptivity

                  if (LW_diag(2)%L_ls_cloud_weight_absorptivity)        &
     &                LW_diag(2)%ls_cloud_weight_absorptivity =         &
     &                LW_diag(2)%ls_cloud_weight_absorptivity +         &
     &                LW_diag(1)%ls_cloud_weight_absorptivity          

                  if (LW_diag(2)%L_cnv_cloud_absorptivity)              & 
     &                LW_diag(2)%cnv_cloud_absorptivity =               &
     &                LW_diag(2)%cnv_cloud_absorptivity +               & 
     &                LW_diag(1)%cnv_cloud_absorptivity

                  if (LW_diag(2)%L_cnv_cloud_weight_absorptivity)       &
     &                LW_diag(2)%cnv_cloud_weight_absorptivity =        &
     &                LW_diag(2)%cnv_cloud_weight_absorptivity +        &
     &                LW_diag(1)%cnv_cloud_weight_absorptivity

                  l_call_lwrad=.true.

               endif

            else if (l_forcing.and.l_rad_step_diag) then
!
! Calculate LW Forcings. Note that forcings are only calculated
! for fluxes. Difference between j_lw=1 (call with forcing in
! place) and j_lw= 2 (call with reference values). N.B. Call to
! radiation schemes done first with j_lw=2 then with j_lw =1
!

               i_off=diagnostic_offset*(j_lw-1)

               If (j_lw >  1) then

                  if (LW_diag(j_lw)%L_clear_olr)                        &
     &                LW_diag(j_lw)%clear_olr =                         &
     &                LW_diag(1)%clear_olr -                            &
     &                LW_diag(j_lw)%clear_olr

                  if (LW_diag(j_lw)%L_surface_down_flux)                &
     &                LW_diag(j_lw)%surface_down_flux =                 &
     &                LW_diag(1)%surface_down_flux -                    &
     &                LW_diag(j_lw)%surface_down_flux

                  if (LW_diag(j_lw)%L_surf_down_clr)                    &
     &                LW_diag(j_lw)%surf_down_clr =                     &
     &                LW_diag(1)%surf_down_clr -                        &
     &                LW_diag(j_lw)%surf_down_clr

                  if (LW_diag(j_lw)%L_clear_hr)                         &
     &                LW_diag(j_lw)%clear_hr =                          &
     &                LW_diag(1)%clear_hr -                             &
     &                LW_diag(j_lw)%clear_hr

                  if (LW_diag(j_lw)%L_net_flux_trop)                    &
     &                LW_diag(j_lw)%net_flux_trop =                     &
     &                LW_diag(1)%net_flux_trop -                        &
     &                LW_diag(j_lw)%net_flux_trop

                  if (LW_diag(j_lw)%L_down_flux_trop)                   &
     &                LW_diag(j_lw)%down_flux_trop =                    &
     &                LW_diag(1)%down_flux_trop -                       &
     &                LW_diag(j_lw)%down_flux_trop

               endif
               l_call_lwrad=.true.
            elseif (l_radiance) then
               if (j_lw == 1) then
                  i_off=0
               else
                  i_off=90+(diagnostic_offset/20)*(j_lw-1)
               endif
                  l_call_lwrad=.true.
            else
               i_off=0
               if (l_rad_step_prog) then
                  l_call_lwrad=.true.
               endif
            endif

! DEPENDS ON: diagnostics_lw
            Call diagnostics_lw(                                        &
     &                     row_length, rows, model_levels               &
     &,                    wet_model_levels, ozone_levels               &
     &,                    cloud_levels                                 &
     &,                    n_rows, global_row_length, global_rows       &
     &,                    halo_i, halo_j, off_x, off_y, me             &
     &,                    n_proc, n_procx, n_procy                     &
     &,                    g_rows, g_row_length, g_datastart            &
     &,                    at_extremity                                 &
     &,                    timestep,i_off                               &
     &,                    T_n, T_inc                                   &
     &,                    q_n, qcl_n, cf_n, cfl_n                      &
     &,                    T_latest, q_latest, qcl_latest               &
     &,                    cf_latest, cfl_latest                        &
     &,                    surflw, OLR                                  &
     &,                    T_incr_diagnostic                            &
     &,                    LW_spectrum(1)%n_band                        &
     &,                    ozone                                        &
     &,                    O3_trop_level                                &
     &,                    O3_trop_height                               &
     &,                    T_trop_level                                 &
     &,                    T_trop_height                                &
     &,                    LW_incs, LWsea                               &
     &,                    lw_spectrum(j_lw)%n_aod_wavel                &
     &,                    l_call_lwrad, j_lw                           &
     &,                                                                 &
#include "argsts.h"
     &                     STASHwork2 )

! DEPENDS ON: timer
            If (Ltimer) Call timer ('Diags   ',4)
         endif
      enddo
#else
      If (L_SCMDiags(SCMDiag_rad)) Then
!       ! Output some SCM diagnostics for LW radiation
!       ! stash 2,161
! DEPENDS ON: scmoutput
        call SCMoutput(T_incr_diagnostic,                               &
             'lw1','LW heating rate minus PC2','K/timestep',            &
             t_acc,d_all,default_streams,'',RoutineName)
! DEPENDS ON: scmoutput
        call SCMoutput(T_incr_diagnostic,                               &
             'lw2','LW heating rate minus PC2','K/day',                 &
             t_mult,d_all,default_streams,'ntspday',RoutineName)
! DEPENDS ON: scmoutput
        call SCMoutput(surflw,                                          &
             'surf_lw','Net surface LW flux','W/m2',                    &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 2,204
! DEPENDS ON: scmoutput
        call SCMoutput(LW_diag(j_lw)%total_cloud_cover,                 &
             'tca_lw','Total cloud amount in LW rad','fraction',        &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 2,205
! DEPENDS ON: scmoutput
        call SCMoutput(OLR,                                             &
             'olr_toa','Outgoing LW','W/m2',                            &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 2,206
! DEPENDS ON: scmoutput
        call SCMoutput(LW_diag(j_lw)%clear_olr,                         &
             'cs_olr','Clear-sky outgoing LW','W/m2',                   &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 2,208
! DEPENDS ON: scmoutput
        call SCMoutput(LW_diag(j_lw)%surf_down_clr,                     &
             'cs_surf_dnlw','Clear-sky down LW flux','W/m2',            &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 2,233
! DEPENDS ON: scmoutput
        call SCMoutput(LW_diag(j_lw)%clear_hr,                          &
             'dt_cslw','Clear-sky LW heating rates','K/s',              &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)
!       ! stash 2,182  Vapour increment
        Do k = 1,wet_model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3dw(i,j,k) = q_latest(i,j,k) - q_n(i,j,k)
            Enddo ! i
          Enddo   ! j
        Enddo     ! k
! DEPENDS ON: scmoutput
        call SCMoutput(work_3dw,                                        &
             'dq_lw','Specific humidity increment lwrad','kg/kg',       &
             t_avg+only_radsteps,d_wet,default_streams,'',RoutineName)
!       ! stash 2,183  liquid water content increment
        Do k = 1,wet_model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3dw(i,j,k) = qcl_latest(i,j,k) - qcl_n(i,j,k)
            Enddo ! i
          Enddo   ! j
        Enddo     ! k
! DEPENDS ON: scmoutput
        call SCMoutput(work_3dw,                                        &
             'dqcl_lw','qcl increment lwrad','kg/kg',                   &
             t_avg+only_radsteps,d_wet,default_streams,'',RoutineName)
!       ! stash 2,192  total cloud fraction increment
        Do k = 1,wet_model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3dw(i,j,k) = cf_latest(i,j,k) - cf_n(i,j,k)
            Enddo ! i
          Enddo   ! j
        Enddo     ! k
! DEPENDS ON: scmoutput
        call SCMoutput(work_3dw,                                        &
             'dcf_lw','bulk cloud fraction increment lwrad',            &
             'fraction',                                                &
             t_avg+only_radsteps,d_wet,default_streams,'',RoutineName)
!       ! stash 2,193  liquid cloud fraction increment
        Do k = 1,wet_model_levels
          Do j = 1,rows
            Do i = 1,row_length
               work_3dw(i,j,k) = cfl_latest(i,j,k) - cfl_n(i,j,k)
            Enddo ! i
          Enddo   ! j
        Enddo     ! k
! DEPENDS ON: scmoutput
        call SCMoutput(work_3dw,                                        &
             'dcfl_lw','liquid cloud fraction increment lwrad',         &
             'fraction',                                                &
             t_avg+only_radsteps,d_wet,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        call SCMoutput(LW_diag(j_lw)%surface_down_flux,                 &
             'surf_dnlw','Downward LW surface flux','W/m2',             &
             t_avg+only_radsteps,d_sl,default_streams,'',RoutineName)

!
! In-radiation cloud diagnostics --------------------------------------
!
! DEPENDS ON: scmoutput
        Call SCMoutput(LW_diag(j_lw)%ls_qcl_rad,                        &
             'ls_qcl_rad','Stratiform cloud liquid water','kg/kg',      &
             t_avg+only_radsteps,d_all,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(LW_diag(j_lw)%ls_qcf_rad,                        &
             'ls_qcf_rad','Stratiform cloud ice water','kg/kg',         &
             t_avg+only_radsteps,d_all,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(LW_diag(j_lw)%cc_qcl_rad,                        &
             'cc_qcl_rad','Convective cloud liquid water','kg/kg',      &
             t_avg+only_radsteps,d_all,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(LW_diag(j_lw)%cc_qcf_rad,                        &
             'cc_qcf_rad','Convective cloud ice water','kg/kg',         &
             t_avg+only_radsteps,d_all,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(LW_diag(j_lw)%ls_cl_rad,                         &
             'ls_cl_rad','Stratiform liquid cloud fraction','',         &
             t_avg+only_radsteps,d_all,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(LW_diag(j_lw)%ls_cf_rad,                         &
             'ls_cf_rad','Stratiform ice cloud fraction','',            &
             t_avg+only_radsteps,d_all,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(LW_diag(j_lw)%cc_cl_rad,                         &
             'ccal_rad','Convective liquid cloud fraction','',          &
             t_avg+only_radsteps,d_all,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(LW_diag(j_lw)%cc_cf_rad,                         &
             'ccaf_rad','Convective ice cloud fraction','',             &
             t_avg+only_radsteps,d_all,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_rad)

#endif

      Deallocate ( T_incr_diagnostic )
!
!     Deallocate the diagnostic space that is no longer required.
!
      if (l_timestep) then
         if (l_rad_step_diag) then
! DEPENDS ON: deallocate_swdiag_ext
            call deallocate_swdiag_ext(2)
! DEPENDS ON: deallocate_lwdiag
            call deallocate_lwdiag(2)
         endif
         if (timestep_number >  1) then
            if (mod(timestep_number,a_lw_radstep_prog) == 0) then

               deallocate(lw_incs_local)
               deallocate(lwsea_local)
               deallocate(olr_local)
               deallocate(top_abs_lw)

! DEPENDS ON: deallocate_swdiag_ext
               call deallocate_swdiag_ext(1)
! DEPENDS ON: deallocate_lwdiag
               call deallocate_lwdiag(1)

            endif
         endif
      elseif (l_forcing) then
         if (l_rad_step_prog) then

! DEPENDS ON: deallocate_swdiag_ext
            call deallocate_swdiag_ext(1)
! DEPENDS ON: deallocate_lwdiag
            call deallocate_lwdiag(1)
         endif
         if (l_rad_step_diag) then
! DEPENDS ON: deallocate_swdiag_ext
            call deallocate_swdiag_ext(2)
! DEPENDS ON: deallocate_lwdiag
            call deallocate_lwdiag(2)
         endif
      elseif (l_radiance) then
         if (l_rad_step_prog) then

            do j_sw=1,n_swcall
! DEPENDS ON: deallocate_swdiag_ext
               call deallocate_swdiag_ext(j_sw)
            enddo
            do j_lw=1,n_lwcall
! DEPENDS ON: deallocate_lwdiag
               call deallocate_lwdiag(j_lw)
            enddo
         endif
      else
         if (l_rad_step_prog) then

! DEPENDS ON: deallocate_swdiag_ext
            call deallocate_swdiag_ext(1)
! DEPENDS ON: deallocate_lwdiag
            call deallocate_lwdiag(1)
         endif
      endif

      End If ! on error_code

      Return
      END SUBROUTINE Glue_Rad
#endif
