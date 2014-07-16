#if defined(ATMOS)
#if defined(CONTROL) || defined(SCMA)
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

#if defined(SCMA)
! Include parameters necessary for calls to SCMoutput...
#include "s_scmop.h"
#endif

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
#include "argsts.h"
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
#endif
#endif
