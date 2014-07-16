#if defined(ATMOS)
#if defined(CONTROL) || defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interface to Atmos Physics parametrizations before S-L advection.
!
! Subroutine Interface:
      SUBROUTINE Atmos_Physics1(                                        &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &

! model dimensions.
     &, row_length, rows, n_rows, land_points, model_levels, wet_levels &
     &, bl_levels, dst_levels, dsm_levels, Ozone_levels, cloud_levels   &
     &, land_ice_points, soil_points, n_cca_levels, ntiles              &
     &, salt_dim1, salt_dim2, salt_dim3, tr_levels, tr_vars             &
     &, co2_dim_len, co2_dim_row, co2_dim_lev                           &
     &, n_arcl_species, n_arcl_compnts, i_arcl_compnts                  &

! Model switches
     &, model_domain, L_regular, L_SEC_VAR, L_EqT                       &
     &, L_Rad_Step, L_Rad_Step_diag,L_Rad_Step_prog                     &
     &, L_Forcing, L_Timestep, L_Radiance, L_Wenyi                      &
     &, L_CAL360, L_microphy, L_emcorr, L_climat_aerosol, Ltimer        &
     &, L_gwd, L_use_ussp, l_taus_scale, l_fix_gwsatn                   &
     &, l_gwd_40km, l_ussp_opaque, sat_scheme                           &
     &, l_use_clearrh, L_ssice_albedo, L_RHCPT, L_murk                  &
     &, L_murk_source, l_murk_bdry, L_MURK_RAD, L_DUST, l_sulpc_so2     &
     &, l_sulpc_nh3, l_soot, l_biomass, l_ocff, l_co2_interactive       &
! Added switch for co2/radiation interaction - kdcorbin, 06/10
     &, l_co2_radiation                                                 &
     &, Lflux_reset,L_clim_aero_hgt,L_HadGEM1_Clim_Aero,L_use_dust      &
     &, L_use_sulphate_autoconv, L_auto_debias, L_use_seasalt_autoconv  &
     &, L_use_seasalt_indirect, L_use_seasalt_direct                    &
     &, L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a              &
     &, L_snow_albedo, L_ctile, L_radiation, L_rain                     &
     &, L_INHOM_CLOUD  , L_USE_BIOGENIC                                 &
     &, L_use_sulpc_direct, L_use_soot_direct, L_use_soot_indirect      &
     &, L_use_soot_autoconv, L_use_bmass_direct                         &
     &, L_use_bmass_indirect, L_use_bmass_autoconv                      &
     &, L_use_ocff_direct, L_use_ocff_indirect, L_use_ocff_autoconv     &
     &, L_use_sulpc_indirect_SW, L_use_sulpc_indirect_LW                &
     &, L_pc2, L_eacf, L_mixing_ratio, l_cry_agg_dep, l_droplet_settle  &
     &, L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup                           &

     &, L_USE_METHOX, L_rad_deg, L_TRIFFID                              &
     &, l_use_stochem_ch4, L_it_melting                                 &
     &, L_ukca, L_USE_ARCL                                              &
     &, L_USE_SPEC_SEA, L_MOD_BARKER_ALBEDO, L_MOD_K_FLUX               &
     &, L_SCVARY, L_VOLCTS                                              &

! model Parameters
     &, rhcrit, cw_sea, cw_land                                         &
     &, A_SW_segments, A_SW_seg_size, A_LW_segments, A_LW_seg_size      &
     &, A_SW_radstep, A_LW_radstep, A_SW_radstep_diag, A_SW_radstep_prog&
     &, A_LW_radstep_diag, A_LW_radstep_prog, aero_bl_levels            &
     &, INHOM_CLOUD_SW, INHOM_CLOUD_LW, DP_CORR_STRAT, DP_CORR_CONV     &
     &, CO2_MMR, SW_alpham, SW_alphac, SW_alphab, SW_dtice              &
     &, dt_bare,dalb_bare_wet,pen_rad_frac,SW_beta                      &
     &, min_trop_level, max_trop_level, GW_kay,GWD_FRC,gwd_fsat, O2MMR  &
     &, N2O_mix_ratio, CH4_mix_ratio, CFC11_mix_ratio, CFC12_mix_ratio  &
     &, C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR                       &
     &, ngrgas , grgas_addr                                             &
     &, cloud_fraction_method,overlap_ice_liquid                        &
     &, ice_fraction_method,ctt_weight,t_weight,qsat_fixed,sub_cld      &
     &, dbsdtbs_turb_0                                                  &
     &, L_seq_mcr,L_autoc_3b,L_autolim_3b,L_autoconv_murk               &
     &, ec_auto,N_drop_land                                             &
     &, N_drop_sea,N_drop_land_cr,N_drop_sea_cr,Ntot_land, Ntot_sea     &
     &, x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic                        &
     &, lsp_ei,lsp_fi,lsp_eic,lsp_fic                                   &
! parameter for stochastic physics random parameters
     &, M_CI                                                            &

! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &

! in coordinate information
     &, r_rho_levels, r_theta_levels, r_at_u, r_at_v, eta_theta_levels  &
     &, eta_rho_levels, delta_lambda, delta_phi, lat_rot_NP, long_rot_NP&

! in time stepping information.
     &, timestep, radiation_timestep                                    &
     &, radiation_tstep_diag, radiation_tstep_prog                      &
     &, val_year, val_day_number, val_hour, val_minute                  &
     &, val_second, timestep_number, PREVIOUS_TIME,istep_cur            &

! trig arrays
     &, sin_theta_longitude, cos_theta_longitude                        &
     &, FV_cos_theta_latitude, sin_theta_latitude                       &

! grid-dependent arrays
     &, f3_at_u, true_longitude, true_latitude,                         &

! diagnostic info
#include "argsts.h"
     & STASHwork1,STASHwork2,STASHwork4,STASHwork6,STASHwork14          &
!
! Additional variables for SCM diagnostics
     &, nSCMDpkgs, L_SCMDiags                                           &
!
! in data fields.
     &, theta, q, qcl, qcf, qcf2, qrain, qgraup, rho, u, v, p, p_star   &
     &, exner_rho_levels, exner_theta_levels                            &
     &, land_sea_mask, p_theta_levels,fland_ctile                       &

     &     ,frac_control                                                &
! ancillary fields and fields needed to be kept from timestep to
! timestep

     &, land_index,rgrain,soot,ntml,cumulus,ice_fract,cca,ccb,cct,cclwp &
     &, ccw,lcbase                                                      &
     &, t_surf, tstar_land_ctile, tstar_sea_ctile, tstar_sice_ctile     &
     &, sice_alb_ctile,land_alb_ctile,snow_depth,snow_depth_sea         &
     &, ozone, SW_incs, LW_incs, dirpar_inc                             &
     &, O3_trop_level, O3_trop_height, T_trop_level, T_trop_height      &
     &, zh, sd_orog_land , orog_grad_xx_land, orog_grad_xy_land         &
     &, orog_grad_yy_land, area_cloud_fraction, cf, cfl, cff            &
     &, aerosol_em, arcl                                                &
     &, albsoil, lai, snow_tile, tile_frac, tstar_tile, z0_tile         &
     &, dOLR_rts, LW_down, SW_tile_rts, ES_SPACE_INTERP, RAD_MASK       &
     &, ch4_stochem, cos_zenith_angle, can_rad_mod                      & 

! in/out
     &, theta_star, q_star, qcl_star, qcf_star, qcf2_star, qrain_star   &
     &, qgraup_star, cf_star, cfl_star, cff_star, u_inc, v_inc          &
     &, energy_correction, sum_eng_fluxes, sum_moist_flux, AEROSOL      &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, SO2, SO4_AITKEN, SO4_ACCU, SO4_DISS, NH3                        &
     &, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld   &
     &, ocff_new, ocff_agd, ocff_cld                                    &
     &, co2, free_tracers, BIOGENIC, ASTEPS_SINCE_TRIFFID               &
!    EAK
!    IN
     &, surf_down_sw,alb_tile                                           &
!sxy
!     &, day,L_TILE_PTS,TILE_PTS,SM_LEVELS,TILE_INDEX                   &
     &, day,TILE_PTS,SM_LEVELS,TILE_INDEX                               &
     &, SNOW_TMP3L,SNOW_RHO1L,TSOIL_TILE,ISNOW_FLG3L,LAND_ALBEDO        & 
     &, l_cable                                                         &
! out fields
     &, ls_rain, ls_snow, micro_tends, unscaled_dry_rho                 &
     &, photosynth_act_rad, rad_hr, surf_radflux, dOLR, SW_tile         &
! Section Information
     &, maxsects, h_sect                                                &

! error information
     &, Error_code  )

      USE CSENARIO_MOD
      Use cv_cntl_mod, Only:                                            &
          lcv_pc2_diag_sh

#if defined(ACCESS)
      USE auscom_cpl_data_mod,                                          &
     &    Only : auscom_salinity, access_tfs, ocn_sss
#endif

      IMPLICIT NONE
!
! Description: This version interfaces to physics schemes in sequence:
!    energy correction              (optional)
!    microphysics (cloud and large scale precipitation schemes)
!    radiation
!    gravity wave drag
!
!          CALLed before Semi-Lagrangian in atmosphere timestep.
! Method:
!
!
! Current Code Owner: Rick Rawlins
!
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component: Control Atmos
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
! for N_INTERNAL_MODEL in typsts.h
#include "csubmodl.h"
#include "typsts.h"
#include "parparm.h"
#include "domtyp.h"
#include "nstypes.h"
#include "c_0_dg_c.h"
#include "c_mdi.h"
#include "c_dust_ndiv.h"
#include "feedback.h"
#if defined(SCMA)
! INOUT SCMop is declared in here
#include "s_scmop.h"
#endif

! Subroutine arguments

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
     &, g_datastart(3,0:n_proc-1)                                       &
     &, me         ! My processor number

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north, sout
                         ! east or west of the processor grid


! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, salt_dim1                                                       &
                      !
     &, salt_dim2                                                       &
                      ! Dimensions of sea-salt aerosol arrays
     &, salt_dim3                                                       &
                      !
     &, land_points                                                     &
                    ! IN No.of land points being processed, can be 0.
     &, model_levels                                                    &
     &, wet_levels                                                      &
     &, bl_levels                                                       &
     &, dst_levels                                                      &
                    ! number of deep soil temperature levels
     &, dsm_levels                                                      &
                    ! number of deep soil moisture levels
     &, Ozone_levels                                                    &
     &, cloud_levels                                                    &
     &, land_ice_points                                                 &
                        ! number of land ice points
     &, soil_points                                                     &
                        ! number of soil points
     &, n_cca_levels                                                    &
                      ! Number of levels for conv cloud
     &, tr_levels                                                       &
                      ! number of free tracer levels
     &, tr_vars                                                         &
                      ! number of free tracers
                      ! amount: 1 for 2D, nlevs for 3D.
     &, ntiles                                                          &
     &,CO2_DIM_LEN                                                      &
                     !\ For dimension 3-D CO2 field to be passed
     &,CO2_DIM_ROW                                                      &
                     !/ to NI_rad_ctl
     &,CO2_DIM_LEV   !/

! Model switches
      Integer                                                           &
     &  model_domain                                                    &
     &, can_rad_mod                                                     &
     &, sat_scheme



      Logical                                                           &
     &  Ltimer   ! true then output some timing information

      Logical                                                           &
     &  L_Rad_Step                                                      &
                        ! true if a radiation timestep
     &, L_Rad_Step_diag                                                 &
                         ! true if fast radiation timestep    (3C)
     &, L_Rad_Step_prog                                                 &
                         ! true if slow radiation timestep    (3C)
     &, L_Forcing                                                       &
                         ! true if radiative forcing required (3C)
     &, L_Timestep                                                      &
                         ! true if new time-stepping is used  (3C)
     &, L_Radiance                                                      &
                         ! true if radiances required         (3C)
     &, L_Wenyi                                                         &
                         ! true if Wenyi scaling required  (3A/3C)
     &, L_CAL360                                                        &
                        ! true if using 360 day calender
     &, L_regular                                                       &
                     ! True if NOT variable resolution
     &, L_SEC_VAR                                                       &
                        ! True if using time varying astronomy
     &, L_EqT                                                           &
                        ! True if including the equation of time
     &, L_microphy                                                      &
                        ! true if using cloud microphysics (set to true)
     &, L_climat_aerosol                                                &
                         ! True if using the climatological aerosol
     &, L_clim_aero_hgt                                                 &
                         ! True if using the prognostic BL depth to
!                        ! specify the BL component of the aerosol
!                        ! climatology.
     &, L_HadGEM1_Clim_Aero                                             &
                             ! True if using HadGEM1 setting for
!                            ! climatological aerosols in FILL3A
     &, L_USE_DUST                                                      &
                    ! Include direct radiative effect of mieral dust
     &, L_USE_BIOGENIC                                                  &
                       ! Include direct radiative effect of biogenic aer
     &, L_use_sulpc_direct                                              &
                                ! Include direct effect of sulphate
     &, L_use_soot_direct                                               &
                                ! Include direct effect of soot
     &, L_use_bmass_direct                                              &
                                ! Include direct effect of biomass smoke
     &, L_use_ocff_direct                                               &
                                ! Include direct effect of fossil-fuel
                                ! organic carbon aerosol
     &, L_use_sulpc_indirect_SW                                         &
                                ! Use sulphate aerosol to determine
     &, L_use_sulpc_indirect_LW                                         &
                                ! cloud droplet number in SW and LW
!                               ! radiation respectively.
     &, L_auto_debias                                                   &
                                ! True if autoconversion debiasing on
     &, L_use_sulphate_autoconv                                         &
                                ! True if sulphate used in microphysics
     &, L_use_seasalt_autoconv                                          &
                                ! True if sea-salt used in microphysics
     &, L_use_seasalt_indirect                                          &
                                ! True if using sea-salt indirect effect
     &, L_use_seasalt_direct                                            &
                                ! True if using sea-salt direct effect
     &, L_use_soot_indirect                                             &
                                ! True if using soot indirect effect
     &, L_use_soot_autoconv                                             &
                                ! True if soot used in microphysics
     &, L_use_bmass_indirect                                            &
                                ! True if using biomass indirect effect
     &, L_use_bmass_autoconv                                            &
                                ! True if biomass used in microphysics
     &, L_use_ocff_indirect                                             &
                                ! True if using OCFF indirect effect
     &, L_use_ocff_autoconv                                             &
                                ! True if OCFF used in microphysics
     &, L_emcorr                                                        &
                        ! true if energy correction scheme used
     &, L_GWD                                                           &
                     ! if true then enables orographic GWD scheme
     &, L_use_ussp                                                      &
                     ! IN if true non-orographic USSP scheme invoked
     &, l_taus_scale                                                    &
                      ! if true then orog surf stress depends on Fr
     &, l_fix_gwsatn                                                    &
                      ! if true then invoke minor bug fixes to gwsatn4a
     &, l_gwd_40km                                                      &
                      ! if true then turn off orographic GWD above 40km
     &, l_ussp_opaque                                                   &
                      ! if true then use opaque lid option in gw_ussp
     &, L_ssice_albedo                                                  &
                       ! Switch on the effect of snow on sea-ice albedo
     &, L_sice_meltponds                                                &
                         ! true if sea ice meltponds albedo scheme
     &, L_sice_scattering                                               &
                          ! true if sea ice internal scattering scheme
     &, L_sice_hadgem1a                                                 &
                        ! true if HadGEM1 albedo bug corrected
     &, l_cable                                                         &
                       ! true if cable albedo is used
     &,   L_USE_METHOX                                                  &
                           ! true if methane oxidation is switched on
     &, L_snow_albedo                                                   &
                       ! True if spectral albedo scheme selected
     &,  l_ctile       ! True if coastal tiling selected.

      Logical                                                           &


     &  L_RHCPT                                                         &
                    ! Switch for 3D RHcrit Diagnostic, not 1D parameter
     &, L_murk                                                          &
                    ! Switch for aerosol
     &, L_murk_source                                                   &
                      ! Switch for source/scavenging of aerosol
     &, L_murk_bdry                                                     &
                      ! switch for bdry for UK MES
     &, L_DUST                                                          &
                    ! mineral dust switch
     &, l_sulpc_so2                                                     &
                    ! Sulphur cycle main switch
     &, l_sulpc_nh3                                                     &
                    ! Sulphur cycle - nh3 switch
     &, l_soot                                                          &
                    ! Soot cycle switch
     &, l_biomass                                                       &
                    ! Biomass aerosol switch
     &, l_ocff                                                          &
                    ! Fossil-fuel organic carbon aerosol switch
     &, l_co2_interactive                                               &
                           ! carbon cycle switch
     &, l_co2_radiation                                                 &
          !carbon/radiation interaction switch - kdcorbin, 06/10
     &, Lflux_reset                                                     &
                    ! true if timestep to reset flux to zero
     &, L_mcr_qcf2                                                      &
                     ! true if using second prognostic cloud ice
     &, L_mcr_qrain                                                     &
                     ! true if using prognostic rain
     &, L_mcr_qgraup                                                    &
                     ! true if using prognostic graupel
     &, L_pc2                                                           &
                    ! true if PC2 cloud scheme is used
     &, L_mixing_ratio                                                  &
                       ! true if mixing ratios used
     &, L_cry_agg_dep                                                   &
                    ! Limit the supersaturation that can be
                    ! removed by deposition depending on
                    ! amount of ice in each category
     &, l_droplet_settle                                                &
                    ! Allow cloud droplets to settle
     &, L_eacf                                                          &
                    ! true if using empirically adjusted cloud fraction
     &, L_it_melting                                                    &
                     ! Use iterative melting
     &, L_MURK_RAD                                                      &
     &, L_rad_deg                                                       &
                   ! Controls the spatial degradation of E-S code
     &,     L_TRIFFID                                                   &
                            ! IN Indicates whether TRIFFID in use.
     &, L_radiation                                                     &
                       ! switch to turn off radiation
     &, L_rain                                                          &
                       ! switch to turn off precip scheme
     &, L_INHOM_CLOUD                                                   &
                       ! switch to simulate inhomogeneous cloud
     &, l_use_stochem_ch4                                               &
     &, l_use_clearrh                                                   &
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
     &, radiation_tstep_diag                                            &
     &, radiation_tstep_prog                                            &
     &, cw_sea                                                          &
                        ! IN threshold cloud liquid water content
                        !    over sea for conversion to ppn
                        !   (kg water per m**3)
     &, cw_land                                                         &
                        ! IN threshold cloud liquid water content
                        !    over land for conversion to ppn
                        !    (kg water per m**3)
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
     &, SW_beta                                                         &
     &, GW_kay                                                          &
     &, GWD_frc                                                         &
     &, gwd_fsat                                                        &
     &, C113MMR                                                         &
                       ! CFC113 mmr
     &, HCFC22MMR                                                       &
                       ! HCFC22 mmr
     &, HFC125MMR                                                       &
                       ! HFC125 mmr
     &, HFC134AMMR     ! HFC134A mmr

      Logical                                                           &
                        !, Intent(IN)
     &    L_seq_mcr                                                     &
                               ! Use sequential updating of mphys
     &   ,L_autoc_3b                                                    &
                               ! Use 3B autoconversion method
     &   ,L_autolim_3b                                                  &
                               ! Use fixed 3B values for the
     &   ,L_autoconv_murk                                               &
                               ! Use murk aerosol to calc. drop number
     &   ,L_psd
                               ! Use generic ice particle size distn.

      Integer                                                           &
                        !, Intent(IN)
     &    cloud_fraction_method                                         &
                                 ! Method for calculating
                               ! total cloud fraction
     &   ,ice_fraction_method  ! Method for calculating ice cloud frac.

      Real                                                              &
                        !, Intent(IN)
     &    overlap_ice_liquid                                            &
                               ! Overlap between ice and liquid phases
     &   ,ctt_weight                                                    &
                               ! Weighting of cloud top temperature
     &   ,t_weight                                                      &
                               ! Weighting of local temperature
     &   ,qsat_fixed                                                    &
                               ! Fixed value of saturation humidity
     &   ,sub_cld                                                       &
                               ! Scaling parameter
     &   ,dbsdtbs_turb_0
                               ! PC2 erosion parameter / s-1

      Real                                                              &
                        !, Intent(IN)
     &    ec_auto                                                       &
                               ! Collision coalescence efficiency
     &   ,N_drop_land                                                   &
                               ! Number of droplets over land / m-3
     &   ,N_drop_sea                                                    &
                               ! Number of droplets over sea / m-3
     &   ,N_drop_land_cr                                                &
                               ! N_drop_land ^ (-1/3) / m
     &   ,N_drop_sea_cr                                                 &
                               ! N_drop_sea ^ (-1/3) / m
     &   ,Ntot_land                                                     &
                               ! Number of droplets over land / m-3
     &   ,Ntot_sea             ! Number of droplets over sea / m-3

      Real                                                              &
                        !, Intent(IN)
     &    x1i                                                           &
                    ! Intercept of aggregate size distribution
     &   ,x1ic                                                          &
                    ! Intercept of crystal size distribution
     &   ,x1r                                                           &
                    ! Intercept of raindrop size distribution
     &   ,x2r                                                           &
                    ! Scaling parameter of raindrop size distribn
     &   ,x4r                                                           &
                    ! Shape parameter of raindrop size distribution
     &   ,ai,  bi                                                       &
                    ! Ice aggregate mass-size relationship m(D)=ai D^bi
     &   ,aic, bic                                                      &
                    ! Ice crystal mass-size relationship m(D)=aic D^bic
     &   ,lsp_ei,  lsp_fi                                               &
                    ! Ice aggregate Best number and Reynolds number 
                    ! relationship: Re(D) =LSP_EI Be^LSP_FI
     &   ,lsp_eic, lsp_fic                                              &
                    ! Ice crystal Best number and Reynolds number 
                    ! relationship: Re(D) =LSP_EIC Be^LSP_FIC

     &   ,M_CI      ! variable to modify ice fall speed for LSPCON3C
                    !  for stochastic physics random parameters     

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

      Real                                                              &
     &  rhcrit(wet_levels)   ! IN Critical relative humidity.
                             ! the values need to be tuned
                             ! for the given set of levels.

! Diagnostics info
      REAL                                                              &
     & STASHWORK1(*)                                                    &
                         ! STASH workspace for section 1 (SW rad)
     &,STASHWORK2(*)                                                    &
                         ! STASH workspace for section 2 (LW rad)
     &,STASHWORK4(*)                                                    &
                         ! STASH workspace for section 4 (LS precip)
     &,STASHWORK6(*)                                                    &
                         ! STASH workspace for section 6 (gw drag)
     &,STASHWORK14(*)    ! STASH workspace for section 14 (Energy cor)

! Data arrays
      Real, Intent (InOut) ::                                           &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels + 1)                                           &
     &, p_theta_levels(1-off_x:row_length+off_x,                        &
     &                   1-off_y:rows+off_y, model_levels)              &
     &, p_star(row_length, rows)                                        &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                   1-off_y:rows+off_y, model_levels + 1)          &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &        wet_levels)                                               &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, qcf2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &        wet_levels)                                               &
     &, qrain(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        wet_levels)                                               &
     &, qgraup(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &        wet_levels)                                               &
     &, cf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &        wet_levels)                                               &
     &, cfl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, cff(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)

      Real                                                              &
     &  energy_correction

! ancillary arrays and fields required to be saved from timestep to
! timestep.

      Real                                                              &
     &  T_surf(row_length, rows)
! The following are only used if coastal tiling is switched on:
      REAL                                                              &
     & FLAND_CTILE(LAND_POINTS)                                         &
                                   ! IN Land fraction on land points.
     &,TSTAR_LAND_CTILE(ROW_LENGTH,ROWS)                                &
!                                  ! IN Land mean sfc temperature (K)
     &,TSTAR_SEA_CTILE(ROW_LENGTH,ROWS)                                 &
!                                  ! IN Open sea sfc temperature (K).
     &,TSTAR_SICE_CTILE(ROW_LENGTH,ROWS)                                &
!                                  ! IN Sea-ice sfc temperature (K).
     &,LAND_ALB_CTILE(ROW_LENGTH,ROWS)                                  &
!                                  ! INOUT Mean land albedo.
     &,SICE_ALB_CTILE(ROW_LENGTH,ROWS)
!                                  ! INOUT Sea-ice albedo.

      logical                                                           &
     &  land_sea_mask(row_length, rows)                                 &
     &, RAD_MASK(row_length, rows)
!  A mask which ensures a chequerboard pattern of radiation calculations
!  over the whole domain (not just one PE)

      Integer                                                           &
     &  land_index (land_points)                                        &
                                      ! set from land_sea_mask
     &, ntml(row_length, rows)

      Logical                                                           &
     &  cumulus (row_length, rows) ! bl convection flag

      Real                                                              &
     &  snow_depth (row_length, rows)                                   &
                                      ! snow/qrclim.snow.(month)
     &, snow_depth_sea (row_length, rows)                               &
                                          ! snow depth on sea ice
     &, sd_orog_land (land_points)                                      &
                                   ! orog/qrparm.orog.stdev
     &, orog_grad_xx_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmaxx
     &, orog_grad_xy_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmaxy
     &, orog_grad_yy_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmayy
     &, zh(row_length, rows)                                            &
                             ! boundary layer height
     &, aerosol_em( 1 : row_length, 1 : rows, 1 :model_levels )

      Real                                                              &
     &  ozone(row_length, rows, ozone_levels)                           &
     &, O3_trop_level(row_length,rows)                                  &
     &, O3_trop_height(row_length,rows)                                 &
     &, ch4_stochem(row_length, rows, model_levels)                     &
     &, cos_zenith_angle(row_length, rows)                              &
     &, T_trop_level(row_length,rows)                                   &
     &, T_trop_height(row_length,rows)                                  &
     &, SW_incs(row_length, rows, 0:model_levels+1)                     &
     &, LW_incs(row_length, rows, 0:model_levels)                       &
     &, dirpar_inc(row_length, rows)                                    &
     &, soot(row_length, rows)                                          &
                                       ! Snow soot
     &, ice_fract (row_length, rows) ! ice/qrclim.ice.(month)

      Real                                                              &
     &  albsoil(land_points)                                            &
     &, lai(land_points, npft)                                          &
     &, rgrain(land_points, ntiles)                                     &
     &, snow_tile(land_points, ntiles)                                  &
     &, tile_frac(land_points, ntype)                                   &
     &, tstar_tile(land_points, ntiles)                                 &
     &, z0_tile(land_points, ntiles)                                    &
     &, dOLR_rts(row_length, rows)                                      &
                                         ! TOA - surface upward LW
     &, LW_down(row_length, rows)                                       &
                                         ! Surface downward LW
     &, SW_tile_rts(land_points, ntiles)                                &
                                         ! Surface net SW on land tiles
     &, ES_SPACE_INTERP(4, row_length, rows)                            &
!              Coeffs for spatial interpolation of radiation quantities
!
!     EAK
     &, alb_tile(land_points,ntiles,4)                                  &
     &, surf_down_sw(row_length,rows,4)


! Aerosol climatology for NWP

#include "arcl_dim.h"
      
      ! Number of requested species within the climatology
      Integer n_arcl_species
      
      ! Corresponding number of requested components
      Integer n_arcl_compnts
      
      ! Model switch for each species
      Logical L_USE_ARCL(NPD_ARCL_SPECIES)
      
      ! Mass-mixing ratios
      Real                                                              &
     &   arcl(row_length, rows, model_levels, n_arcl_compnts)
     
      ! Array index of each component
       Integer i_arcl_compnts(NPD_ARCL_COMPNTS)

! Convection Scheme

      Real    :: ccw(row_length, rows, wet_levels)
      Integer :: lcbase(row_length, rows)
 
      Real                                                              &
     &  cca (row_length, rows, n_cca_levels)                            &
     &, cclwp(row_length, rows) ! condensed water path (KG/M**2)

      Integer                                                           &
     &  ccb (row_length, rows)                                          &
     &, cct (row_length, rows)

! Co-ordinate arrays
      Real                                                              &
           ! local vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
     &  r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)

! Trig arrays
      real                                                              &
     &  cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)                          &
     &, sin_theta_latitude  (row_length, rows)                          &
     &, FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)

! Grid-dependent arrays
      Real                                                              &
     &  f3_at_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y)          &
     &, true_longitude(row_length, rows)                                &
     &, true_latitude(row_length, rows)

! time information for current timestep
      Integer                                                           &
     &  val_year                                                        &
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second                                                      &
     &, timestep_number                                                 &
     &,         PREVIOUS_TIME(7)                                        &
     &,     ASTEPS_SINCE_TRIFFID   ! INOUT  Number of atmospheric
!                                  !        timesteps since last call
!                                  !        to TRIFFID.

!===sxy arguments with call cable_RADNum_interface
     Integer istep_cur
     INTEGER       &
!       LAND_PTS    & ! IN No of land points being processed.
       SM_LEVELS   & ! IN No. of soil moisture levels
      ,TILE_INDEX(land_points,NTILES)  &  ! OUT Index of tile points
      ,TILE_PTS(NTILES)             &  ! OUT Number of tile points
      ,ISNOW_FLG3L(land_points,NTILES) &
      ,day 
      REAL                                  &
       LAND_ALBEDO(ROW_LENGTH,ROWS,4) &! Land albedo calculated by Cable
      ,SNOW_RHO1L(land_points,NTILES)          &! snow cover in the ice tile.
      ,TSOIL_TILE(land_points,NTILES,SM_LEVELS)&! Mean snow density  (or 1 layer)
      ,SNOW_TMP3L(land_points,NTILES,3)         ! Snow temperature (3 layer)
!     LOGICAL                                &
!       L_TILE_PTS(land_points,ntiles)
!===sxy

! Diagnostic variables
      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP

      Real                                                              &
     &  area_cloud_fraction(row_length, rows, wet_levels)

! arguments with intent in/out. ie: input variables changed on output.

      Real, Intent (InOut) ::                                           &
     &  theta_star(1-off_x:row_length+off_x,                            &
     &               1-off_y:rows+off_y, model_levels)                  &
     &, q_star(1-off_x:row_length+off_x,                                &
     &           1-off_y:rows+off_y, wet_levels)                        &
     &, qcl_star(1-off_x:row_length+off_x,                              &
     &             1-off_y:rows+off_y, wet_levels)                      &
     &, qcf_star(1-off_x:row_length+off_x,                              &
     &             1-off_y:rows+off_y, wet_levels)                      &
     &, qcf2_star(1-off_x:row_length+off_x,                             &
     &             1-off_y:rows+off_y, wet_levels)                      &
     &, qrain_star(1-off_x:row_length+off_x,                            &
     &             1-off_y:rows+off_y, wet_levels)                      &
     &, qgraup_star(1-off_x:row_length+off_x,                           &
     &             1-off_y:rows+off_y, wet_levels)                      &
     &, cf_star(1-off_x:row_length+off_x,                               &
     &           1-off_y:rows+off_y, wet_levels)                        &
     &, cfl_star(1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, wet_levels)                        &
     &, cff_star(1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, wet_levels)                        &
     &, u_inc(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &        model_levels)                                             &
     &, v_inc(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,           &
     &        model_levels)                                             &
     &, sum_eng_fluxes(row_length,rows)                                 &
                                         ! sum atmosphere fluxes
     &, sum_moist_flux(row_length,rows)  ! sum moist fluxes

      Real, Intent(InOut) ::                                            &
     &  aerosol( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,      &
     &           model_levels )                                         &
     &, DUST_DIV1( 1-OFF_X : ROW_LENGTH+OFF_X, 1-OFF_Y : ROWS+OFF_Y,    &
     &           MODEL_LEVELS )                                         &
     &, DUST_DIV2( 1-OFF_X : ROW_LENGTH+OFF_X, 1-OFF_Y : ROWS+OFF_Y,    &
     &           MODEL_LEVELS )                                         &
     &, DUST_DIV3( 1-OFF_X : ROW_LENGTH+OFF_X, 1-OFF_Y : ROWS+OFF_Y,    &
     &           MODEL_LEVELS )                                         &
     &, DUST_DIV4( 1-OFF_X : ROW_LENGTH+OFF_X, 1-OFF_Y : ROWS+OFF_Y,    &
     &           MODEL_LEVELS )                                         &
     &, DUST_DIV5( 1-OFF_X : ROW_LENGTH+OFF_X, 1-OFF_Y : ROWS+OFF_Y,    &
     &           MODEL_LEVELS )                                         &
     &, DUST_DIV6( 1-OFF_X : ROW_LENGTH+OFF_X, 1-OFF_Y : ROWS+OFF_Y,    &
     &           MODEL_LEVELS )                                         &
     &, so2( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,          &
     &       model_levels )                                             &
     &, so4_aitken( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,   &
     &              model_levels )                                      &
     &, so4_accu( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,     &
     &            model_levels )                                        &
     &, so4_diss( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,     &
     &            model_levels )                                        &
     &, nh3( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,          &
     &       model_levels )                                             &
     &, soot_new( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,     &
     &             model_levels )                                       &
     &, soot_agd( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,     &
     &             model_levels )                                       &
     &, soot_cld( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,     &
     &             model_levels )                                       &
     &, bmass_new( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,    &
     &             model_levels )                                       &
     &, bmass_agd( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,    &
     &             model_levels )                                       &
     &, bmass_cld( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,    &
     &             model_levels )                                       &
     &, ocff_new( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,     &
     &             model_levels )                                       &
     &, ocff_agd( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,     &
     &             model_levels )                                       &
     &, ocff_cld( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,     &
     &             model_levels )                                       &
     &, co2( 1-off_x : row_length+off_x, 1-off_y : rows+off_y,          &
     &       model_levels )                                             &
     &, free_tracers( 1-off_x : row_length+off_x, 1-off_y : rows+off_y, &
     &                tr_levels, tr_vars )                              
     
      Real, Intent(In) ::                                               &
     &  biogenic(row_length, rows, model_levels)

! arguments with intent out. ie: output variables.

! arrays holding information to be passed between physics
! routines.

      Real                                                              &
     &  ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, micro_tends(row_length, rows, bl_levels, 2)
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)

! Radiation fields 1. SW & common with LW.
      Real                                                              &
     &  photosynth_act_rad(row_length, rows)                            &
                                             ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
     &, rad_hr(row_length, rows, bl_levels, 2)                          &
                                               !
!                               ! BL radiative (LW,SW) heating rates
     &, surf_radflux(row_length, rows)                                  &
     &, dOLR(row_length, rows)                                          &
                                    ! TOA - surface upward LW
     &, SW_tile(land_points, ntiles)! Surface net SW on land tiles

! Position of greenhouse gases in free_tracers array
      INTEGER, INTENT(IN) :: ngrgas
      INTEGER, INTENT(IN) :: grgas_addr(ngrgas)
!  Additional variables for SCM diagnostics (dummy in full UM)
      Integer                                                           &
     &  nSCMDpkgs             ! No of SCM diagnostics packages

      Logical                                                           &
     &  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

      Integer                                                           &
     &   maxsects

      Character*3                                                       &
     &   h_sect(0:maxsects)
      Integer                                                           &
     &  Error_code

! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='Atm_Physics1')
      Real                                                              &
     &  amp                                                             &
                        ! amplitude of diurnal variation in tracers
     &, tau_decay                                                       &
                        ! time constant for decay of tracer to clim
     &, clim_murk_land                                                  &
                        ! climatological murk value over land points
     &, clim_murk_sea   ! climatological murk value over sea points
      Parameter (                                                       &
     &            amp=0.7                                               &
     &,           tau_decay=1.728E5                                     &
     &,           clim_murk_land=25.0                                   &
     &,           clim_murk_sea=12.5                                    &
     &           )

! Local scalars:

! loop counters
      Integer                                                           &
     &  i, j, k                                                         &
     &, l

! local variables
      Integer                                                           &
     &  rhc_row_length                                                  &
                        ! Row length for RHcrit array
     &, rhc_rows                                                        &
                        ! Row number for RHcrit array
     &, lspice_dim1,lspice_dim2,lspice_dim3                             &
!                       ! Required for array dimensions MCR_CTL2
     &, sulp_dim1                                                       &
                                ! dimensions for sulphate arrays in
     &, DUST_DIM1                                                       &
                                ! dimensions for mineral dust arrays i
     &, DUST_DIM2                                                       &
                                ! in rad_ctl
     &, BIOGENIC_DIM1                                                   &
                                ! dimensions of biogenic arrays in
     &, BIOGENIC_DIM2                                                   &
                                !
     &, sulp_dim2                                                       &
                                !   RAD_CTL
     &, soot_dim1, soot_dim2                                            &
                                ! dimensions of soot arrays in RAD_CTL
     &, bmass_dim1, bmass_dim2                                          &
                                ! dimensions of biomass arrays in radn
     &, ocff_dim1, ocff_dim2                                            &
                                ! dimensions of OCFF arrays in radiation
     &, arcl_dim1, arcl_dim2                                            &
                     ! dimensions of aerosol clim for NWP arrays in radn
     &, i_start                                                         &
                        ! Row start point for polar row tidying
     &, istat           ! Status (error code) indicator

! Local data arrays

      Real                                                              &
     &  T_n(row_length, rows, model_levels)                             &
     &, q_n(row_length, rows, wet_levels)                               &
     &, qcl_n(row_length, rows, wet_levels)                             &
     &, qcf_n(row_length, rows, wet_levels)                             &
     &, cf_n(row_length, rows, wet_levels)                              &
     &, cfl_n(row_length, rows, wet_levels)                             &
     &, cff_n(row_length, rows, wet_levels)                             &
     &, dirpar_local(row_length, rows)                                  &
     &, u_1(salt_dim1, salt_dim2)                                       &
     &, v_1(salt_dim1, salt_dim2)                                       &
     &, co2_3D(co2_dim_len, co2_dim_row, co2_dim_lev)                   &
     &, u_1_mean                                                        &
     &, v_1_mean                                                        &
     &, frac_control(land_points,ntype)

      Real                                                              &
     &  theta_inc(row_length,rows, model_levels)                        &
     &, q_inc(row_length,rows, wet_levels)                              &
     &, qcl_inc(row_length,rows, wet_levels)                            &
     &, qcf_inc(row_length, rows, wet_levels)                           &
     &, cf_inc(row_length, rows, wet_levels)                            &
     &, cfl_inc(row_length, rows, wet_levels)                           &
     &, cff_inc(row_length, rows, wet_levels)

      ! Local Arrays to store additional microphysics fields if in use
      Real, Dimension (:,:,:), Allocatable ::                           &
     &  qcf2_n,   qrain_n,   qgraup_n                                   &
     &, qcf2_inc, qrain_inc, qgraup_inc                                 &
     &, T_inc_diag,  q_inc_diag,   qcl_inc_diag, qcf_inc_diag           &
     &, cf_inc_diag, cfl_inc_diag, cff_inc_diag

      Real                                                              &
     &  p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, p_layer_centres(row_length, rows, 0:model_levels)               &
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = pstar, and at top = 0.
     &, exner_layer_boundaries(row_length, rows, 0:model_levels)        &
     &, exner_layer_centres(row_length, rows, 0:model_levels)

      Real :: moist(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,   &
              wet_levels)
      ! holds total moisture for conversion from wet to dry density
      Real :: unscaled_dry_rho(1-off_x:row_length+off_x,                &
                              1-off_y:rows+off_y, model_levels)
      ! unscaled dry density
      Real :: weight1
      Real :: weight2
      Real :: weight3
      Real :: temp


      REAL                                                              &
     & FLAND(LAND_POINTS)                                               &
                                   ! Land fraction on land points.
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
                                   ! Land mean sfc temperature (K)
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! Open sea sfc temperature (K).
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! Sea-ice sfc temperature (K).
     &,LAND_ALB(ROW_LENGTH,ROWS)                                        &
                                   ! Mean land albedo.
     &,SICE_ALB(ROW_LENGTH,ROWS)   ! Sea-ice albedo.

      LOGICAL                                                           &
     & land0p5(row_length, rows)

! Diagnostics controlled by Diagnostic switches

! fields output by ls_cld not using stash flag

! Energy correction work variables

      Real                                                              &
     &  tot_precip_scaled(row_length, rows)                             &
     &, lclf

! Temporary logicals
      Logical                                                           &
     &  L_zero_boundaries                                               &
     &, L_zero_halos                                                    &
     &, L_use_dirpar

      INTEGER                                                           &
     &  level

      REAL, ALLOCATABLE :: grgas_field(:,:,:,:)

      REAL :: ltfs

#if defined(ACCESS)
      ltfs = access_tfs
#else
      ltfs = tfs
#endif

! ----------------------------------------------------------------------
! Section INI. Initialisation of variables.
! ----------------------------------------------------------------------

! Set convective outputs to zero in case convection scheme not called.
! Set radiation outputs to zero in case convection scheme not called.
! This code can be removed later.
! set temporary logicals to disabled un-called physics.
      Do j = 1, rows
        Do i = 1, row_length
          photosynth_act_rad(i,j) = 0.
          surf_radflux(i,j) = 0.
        End Do
      End Do

! set p at layer boundaries.
! NB: some arrays have haloes but are unset, if never used they will
!     be removed.
      Do j = 1, rows
        Do i = 1, row_length
          p_layer_boundaries(i,j,0) = p_star(i,j)
          p_layer_centres(i,j,0) = p_star(i,j)
          exner_layer_boundaries(i,j,0) = (p_layer_boundaries(i,j,0)/   &
     &                                     p_zero)**kappa
          exner_layer_centres(i,j,0) = exner_layer_boundaries(i,j,0)

        End Do
      End Do
      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
            p_layer_boundaries(i,j,k) = p(i,j,k+1)
            p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
            exner_layer_boundaries(i,j,k) = exner_rho_levels(i,j,k+1)
            exner_layer_centres(i,j,k) = exner_theta_levels(i,j,k)
          End Do
        End Do
      End Do
      k = model_levels
      Do j = 1, rows
        Do i = 1, row_length
          p_layer_boundaries(i,j,model_levels) = 0.0
          p_layer_centres(i,j,model_levels) =                           &
     &                           p_theta_levels(i,j,model_levels)
          exner_layer_boundaries(i,j,k) = 0.
          exner_layer_centres(i,j,k) = exner_theta_levels(i,j,k)
        End Do
      End Do

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_inc(i,j,k) = 0.0
          End Do
        End Do
      End Do
      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            q_inc(i,j,k) = 0.0
            qcl_inc(i,j,k) = 0.0
            qcf_inc(i,j,k) = 0.0
            cf_inc(i,j,k) = 0.0
            cfl_inc(i,j,k) = 0.0
            cff_inc(i,j,k) = 0.0
          End Do
        End Do
      End Do
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            u_inc(i,j,k) = 0.0
          End Do
        End Do
      End Do
      Do k = 1, model_levels
        Do j = 1, n_rows
          Do i = 1, row_length
            v_inc(i,j,k) = 0.0
          End Do
        End Do
      End Do
      IF (L_TRIFFID) THEN
!-----------------------------------------------------------------------
! Increment counter for number of atmosphere timesteps since last
! call to TRIFFID vegetation model
!-----------------------------------------------------------------------
        ASTEPS_SINCE_TRIFFID = ASTEPS_SINCE_TRIFFID + 1
      ENDIF

      If (L_mcr_qcf2) Then  ! Second cloud ice variable in use
        Allocate ( qcf2_inc(row_length, rows, wet_levels) )
        qcf2_inc(:,:,:) = 0.0
      Else
        Allocate ( qcf2_inc(1,1,1) )
      End If

      If (L_mcr_qrain) Then  ! Prognostic rain in use
        Allocate ( qrain_inc(row_length, rows, wet_levels) )
        qrain_inc(:,:,:) = 0.0
      Else
        Allocate ( qrain_inc(1,1,1) )
      End If

      If (L_mcr_qgraup) Then  ! Prognostic graupel in use
        Allocate ( qgraup_inc(row_length, rows, wet_levels) )
        qgraup_inc(:,:,:) = 0.0
      Else
        Allocate ( qgraup_inc(1,1,1) )
      End If

! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics:
! ----------------------------------------------------------------------
      IF(L_CTILE)THEN

        DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            TSTAR_LAND(I,J)=TSTAR_LAND_CTILE(I,J)
            TSTAR_SEA(I,J)=TSTAR_SEA_CTILE(I,J)
            TSTAR_SICE(I,J)=TSTAR_SICE_CTILE(I,J)
            LAND_ALB(I,J)=LAND_ALB_CTILE(I,J)
            SICE_ALB(I,J)=SICE_ALB_CTILE(I,J)
            LAND0P5(I,J)=.FALSE.
          ENDDO
        ENDDO

        DO L=1,LAND_POINTS
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          FLAND(L)=FLAND_CTILE(L)
          IF(FLAND(L) >= 0.5)LAND0P5(I,J)=.TRUE.
        ENDDO

      ELSE

        DO L=1,LAND_POINTS
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          IF(LAND_SEA_MASK(I,J))THEN
            FLAND(L)=1.0
          ELSE
            FLAND(L)=0.0
          ENDIF
        ENDDO

        DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            TSTAR_LAND(I,J)=T_SURF(I,J)
            IF(.NOT.LAND_SEA_MASK(I,J))THEN
              IF(ICE_FRACT(I,J) <= 0.0)THEN
                TSTAR_SEA(I,J)=T_SURF(I,J)
                TSTAR_SICE(I,J)=T_SURF(I,J)
              ELSE
#if defined(ACCESS)
                IF (ocn_sss) THEN
                   ltfs =  ZeroDegC - 0.054 * auscom_salinity(I,J)
                END IF
#endif
                TSTAR_SEA(I,J)=LTFS
                TSTAR_SICE(I,J)=(T_SURF(I,J)                            &
     &            -(1.-ICE_FRACT(I,J))*TSTAR_SEA(I,J))/ICE_FRACT(I,J)
              ENDIF
            ELSE
              TSTAR_SEA(I,J)=T_SURF(I,J)
              TSTAR_SICE(I,J)=T_SURF(I,J)
            ENDIF
            LAND0P5(I,J)=LAND_SEA_MASK(I,J)

            LAND_ALB(I,J)=RMDI
            SICE_ALB(I,J)=RMDI
          ENDDO
        ENDDO
      ENDIF

! ----------------------------------------------------------------------
! Create an unscaled dry density variable
! ----------------------------------------------------------------------
      moist= q + qcl + qcf
      If(L_mcr_qcf2)Then
        moist      = moist + qcf2
      End If
      If(L_mcr_qrain)Then
        moist      = moist + qrain
      End If
      If(L_mcr_qgraup)Then
        moist      = moist + qgraup
      End If

      k = 1
        Do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            unscaled_dry_rho(i,j,k) = rho(i,j,k) *                      &
                  (1. - moist(i,j,k) )/(r_rho_levels(i,j,k)**2)
          End Do
        End Do

      Do k = 2, wet_levels
        Do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            weight2 = r_rho_levels(i,j,k)                               &
                      - r_theta_levels(i,j,k-1)
            weight1 = r_theta_levels(i,j,k)                             &
                      - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k)                             &
                      - r_theta_levels(i,j,k-1)
            temp = ( weight2 * (1. - moist(i,j,k)  ) +                  &
                     weight1 * (1. - moist(i,j,k-1)) )                  &
                   / weight3
            unscaled_dry_rho(i,j,k) = rho(i,j,k) *                      &
                            temp/(r_rho_levels(i,j,k)**2)
          End Do
        End Do
      End Do

      k = wet_levels + 1
      If ( k <= model_levels ) Then
        Do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            weight2 = r_rho_levels(i,j,k)                               &
                      - r_theta_levels(i,j,k-1)
            weight1 = r_theta_levels(i,j,k)                             &
                      - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k)                             &
                      - r_theta_levels(i,j,k-1)
            temp = ( weight2  +                                         &
                     weight1 * (1. - moist(i,j,k-1)) )                  &
                   / weight3
            unscaled_dry_rho(i,j,k) = rho(i,j,k) *                      &
                                     temp/(r_rho_levels(i,j,k)**2)
          End Do
        End Do
      End If !  k <= model_levels

      Do k = wet_levels+2, model_levels
        Do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            unscaled_dry_rho(i,j,k) = rho(i,j,k)                        &
                                      /(r_rho_levels(i,j,k)**2)
          End Do
        End Do
      End Do

! ----------------------------------------------------------------------
! Section ENG.1  Add energy correction increments to temperature
! ----------------------------------------------------------------------
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Energy Correct.',5)

      If (L_emcorr) then

        If (Lflux_reset) then
! reinitialise net flux field at beginning of energy correction period
            Do j = 1, rows
              Do i = 1, row_length
              sum_eng_fluxes(i,j)=0.0
              sum_moist_flux(i,j)=0.0
            Enddo
          Enddo
        Endif

! Add energy correction increments every timestep.
! This is a temperature increment

! DEPENDS ON: timer
        If (Ltimer) Call timer ('add_eng_corr',3)

! DEPENDS ON: add_eng_corr
        Call add_eng_corr (energy_correction,Theta_inc,row_length,rows, &
     &                     model_levels,timestep,at_extremity,          &
#include "argsts.h"
     &                         STASHwork14)

! DEPENDS ON: timer
        If (Ltimer) Call timer ('add_eng_corr',4)

      End If    ! (L_emcorr)

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Energy Correct.',6)

! ----------------------------------------------------------------------
! Section Set-up time-level n
! ----------------------------------------------------------------------

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            T_n(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
          End Do
        End Do
      End Do
      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            q_n(i,j,k) = q(i,j,k)
            qcl_n(i,j,k) = qcl(i,j,k)
            qcf_n(i,j,k) = qcf(i,j,k)
            cf_n(i,j,k)  = cf(i,j,k)
            cfl_n(i,j,k) = cfl(i,j,k)
            cff_n(i,j,k) = cff(i,j,k)
          End Do
        End Do
      End Do
      If (L_CO2_INTERACTIVE) then
       Do k = 1, co2_dim_lev
         Do j = 1, co2_dim_row
           Do i = 1, co2_dim_len
             co2_3D(i,j,k) = co2(i,j,k)
           End Do
         End Do
       End Do
      Else
       co2_3D(1,1,1) = 0.0
      End If

      If (L_mcr_qcf2) Then  ! Second cloud ice variable in use
        Allocate ( qcf2_n(row_length, rows, wet_levels) )
        qcf2_n(:,:,:) = qcf2(1:row_length, 1:rows, :)
      Else
        Allocate ( qcf2_n(1,1,1) )
      End If

      If (L_mcr_qrain) Then  ! Prognostic rain in use
        Allocate ( qrain_n(row_length, rows, wet_levels) )
        qrain_n(:,:,:) = qrain(1:row_length, 1:rows, :)
      Else
        Allocate ( qrain_n(1,1,1) )
      End If

      If (L_mcr_qgraup) Then  ! Prognostic graupel in use
        Allocate ( qgraup_n(row_length, rows, wet_levels) )
        qgraup_n(:,:,:) = qgraup(1:row_length, 1:rows, :)
      Else
        Allocate ( qgraup_n(1,1,1) )
      End If

      If (L_use_seasalt_direct .OR. L_use_seasalt_indirect              &
     &                         .OR. L_use_seasalt_autoconv) Then
        Do j = 1, rows
          Do i = 1, row_length
            u_1(i,j) = u(i,j,1)
            v_1(i,j) = v(i,j,1)
          End Do
        End Do

! Tidy up at North Pole; not done at South Pole as it's land.
! Polar values are mean of 1st interior row:
        If (at_extremity(PNorth)) Then
! Start point of first interior (i.e. non-polar) row:
          i_start  = (rows - 2) * row_length + 1
! Sum over points on PEs in order along first interior row:
#if defined(REPROD)
          Call gcg_rvecsumr(row_length*rows, row_length, i_start, 1,    &
     &                      u_1, proc_row_group, istat, u_1_mean)
          Call gcg_rvecsumr(row_length*rows, row_length, i_start, 1,    &
     &                      v_1, proc_row_group, istat, v_1_mean)
#else
          Call gcg_rvecsumf(row_length*rows, row_length, i_start, 1,    &
     &                      u_1, proc_row_group, istat, u_1_mean)
          Call gcg_rvecsumf(row_length*rows, row_length, i_start, 1,    &
     &                      v_1, proc_row_group, istat, v_1_mean)
#endif
          u_1_mean=u_1_mean/global_row_length
          v_1_mean=v_1_mean/global_row_length
          Do i = 1, row_length
            u_1(i, rows) = u_1_mean
            v_1(i, rows) = v_1_mean
          End Do
        Endif

      Endif

! In the LAM set qcl_n, qcf_n and cloud fractions to zero on model
! boundaries. This avoids failures do to inconsistences in these fields.
! As the increments on the boundary are purely from the boundary
! conditions we are free to do anything sensible to these values.
! Note that this only applies to the variables in the physics.
      If (model_domain  ==  mt_LAM) Then

        L_zero_boundaries=.TRUE.
        L_zero_halos=.FALSE.

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,QCL_N,               &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,QCF_N,               &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

        If (L_mcr_qcf2)                                                 &
                        ! prognostic second cloud ice in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qcf2_n,             &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

        If (L_mcr_qrain)                                                &
                         ! prognostic rain in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qrain_n,            &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

        If (L_mcr_qgraup)                                               &
                          ! prognostic graupel in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qgraup_n,           &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,                     &
     &   AREA_CLOUD_FRACTION,                                           &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,                     &
     &   cf_n,                                                          &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,                     &
     &   cfl_n,                                                         &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,                     &
     &   cff_n,                                                         &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)


      End If

! ----------------------------------------------------------------------
! Section Microphysics. Call microphys_ctl routine
! ----------------------------------------------------------------------

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Microphysics',5)
! DEPENDS ON: timer
      If (Ltimer) Call timer ('MICROPHYS_CTL',3)
!
      If (L_RHCPT) Then
!     Dimension diagnostic 3D RHcrit array
        rhc_row_length = row_length
        rhc_rows = rows
      Else
!     RHcrit will be a 1D parametrized array input from user interface
        rhc_row_length = 1
        rhc_rows = 1
      End if
      IF(L_rain)then
!
! Only allow dynamic allocation of full space for arrays for 3D
! precipitation diagnostics if they are being used. Otherwise save space
! and give them a minimum size of 1 by 1.
#if !defined(SCMA)
      IF ( SF(222,4) .OR. SF(223,4) .OR. SF(224,4) .OR. SF(225,4)       &
     &               .OR. L_DUST                                        &
     &               .OR. SF(227,4) .OR. L_SULPC_SO2                    &
     &               .OR. L_SOOT .OR. L_BIOMASS) THEN
        LSPICE_DIM1 = row_length
        LSPICE_DIM2 = rows
        LSPICE_DIM3 = wet_levels
      ELSE
        LSPICE_DIM1 = 1
        LSPICE_DIM2 = 1
        LSPICE_DIM3 = 1
      END IF
#else
        LSPICE_DIM1 = row_length
        LSPICE_DIM2 = rows
        LSPICE_DIM3 = wet_levels
#endif
!
! DEPENDS ON: microphys_ctl
      Call microphys_ctl (                                              &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &

! model dimensions.
     &, row_length, rows, rhc_row_length, rhc_rows, n_rows, land_points &
     &, model_levels, wet_levels, bl_levels                             &
     &, lspice_dim1,lspice_dim2,lspice_dim3                             &
     &, dst_levels, dsm_levels, Ozone_levels, cloud_levels              &
     &, land_ice_points, soil_points                                    &
     &, n_cca_levels                                                    &
     &, salt_dim1, salt_dim2, salt_dim3                                 &

! Model switches
     &, model_domain, L_CAL360, Ltimer, L_RHCPT, L_Murk, L_DUST         &
     &, l_sulpc_so2, l_sulpc_nh3, l_soot, l_biomass, l_ocff             &
     &, L_use_sulphate_autoconv, L_auto_debias, L_use_seasalt_autoconv  &
     &, L_use_soot_autoconv, L_use_bmass_autoconv, l_use_ocff_autoconv  &
     &, L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_it_melting             &
     &, L_eacf, L_mixing_ratio, l_cry_agg_dep, l_droplet_settle         &
     &, L_pc2, L_use_biogenic                                           &

! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &

! in coordinate information
     &, r_rho_levels, r_theta_levels, r_at_u, r_at_v                    &
     &, eta_theta_levels, eta_rho_levels                                &
     &, delta_lambda, delta_phi                                         &
     &, lat_rot_NP, long_rot_NP                                         &

! in time stepping information.
     &, timestep, radiation_timestep                                    &
     &, radiation_tstep_diag                                            &
     &, radiation_tstep_prog                                            &
     &, val_year, val_day_number, val_hour, val_minute                  &
     &, val_second, timestep_number                                     &

! Model parameters
     &, rhcrit, cw_sea, cw_land                                         &
     &, cloud_fraction_method,overlap_ice_liquid                        &
     &, ice_fraction_method,ctt_weight,t_weight,qsat_fixed,sub_cld      &
     &, L_seq_mcr,L_autoc_3b,L_autolim_3b,L_autoconv_murk               &
     &, ec_auto,N_drop_land                                             &
     &, N_drop_sea,N_drop_land_cr,N_drop_sea_cr,Ntot_land, Ntot_sea     &
     &, x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic                        & 
     &, lsp_ei,lsp_fi,lsp_eic,lsp_fic                                   &
! Primary fields passed in
     &, T_n, q_n, qcl_n, qcf_n, qcf2_n, qrain_n, qgraup_n               &
     &, cf_n, cfl_n, cff_n                                              &
     &, snow_depth                                                      &
     &, land_sea_mask, ice_fract                                        &
     &, p_layer_centres, p_layer_boundaries, exner_theta_levels         &
     &, rho                                                             &
     &, AEROSOL                                                         &
     &, DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6&
     &, SO2, NH3, SO4_AITKEN, SO4_ACCU, SO4_DISS                        &
     &, soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld             &
     &, ocff_agd, ocff_cld, biogenic                                    &

! Other fields passed in
     &, energy_correction, cca, ccb, cct, p, p_star, p_theta_levels     &
     &, exner_rho_levels                                                &
     &, ntml, cumulus                                                   &
     &, fland, land_index                                               &
     &, u_1, v_1                                                        &
! Variables for stochastic physics random parameters
     &, M_CI                                                            &
 

! diagnostic info
     &     ,                                                            &
#include "argsts.h"
     &  STASHwork4                                                      &
!
! SCM diagnostics switches (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! Increment fields passed in/out
     &, Theta_inc, q_inc, qcl_inc, qcf_inc                              &
     &, qcf2_inc, qrain_inc, qgraup_inc                                 &
     &, cf_inc, cfl_inc, cff_inc                                        &


! Fields required elsewhere
     &, ls_rain, ls_snow, micro_tends                                   &

! Section information
     &, maxsects, h_sect                                                &
! error information
     &, Error_code  )

      else
        ls_rain(:,:) = 0.0
        ls_snow(:,:) = 0.0
           micro_tends(:,:,:,:) = 0.0
      endif
! DEPENDS ON: timer
      If (Ltimer) Call timer ('MICROPHYS_CTL',4)
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Microphysics',6)

! ----------------------------------------------------------------------
! Section Turbulence. Call pc2_turbulence_ctl routine
! ----------------------------------------------------------------------
! Earlier versions of PC2 had the erosion term included here,
! parallel with the rest of the slow physics. We have since decided to
! move the erosion to be parallel with the PC2 response to convection
! since this results in improved numerical balances in shallow
! convection. 
! In the case of PC2 using diagnostic shallow cloud it is better
! to do the erosion term here so the call is now under a switch for this.

      If (l_pc2 .and. lcv_pc2_diag_sh) then

! DEPENDS ON: pc2_turbulence_ctl
        Call pc2_turbulence_ctl (                                       &

! Parallel variables
     &    halo_i, halo_j, off_x, off_y, global_row_length, global_rows  &
     &,   proc_row_group, proc_col_group, at_extremity, n_proc, n_procx &
     &,   n_procy, neighbour, g_rows, g_row_length, g_datastart, me     &

! model dimensions.
     &,   row_length, rows, n_rows                                      &
     &,   model_levels, wet_levels                                      &

! timing
     &,   Ltimer                                                        &

! logical control
     &,   L_mixing_ratio                                                &

! Model parameters
     &,   dbsdtbs_turb_0                                                &

! in time stepping information.
     &,   timestep                                                      &

! Primary fields passed in
     &,   T_n, q_n, qcl_n, qcf_n, cf_n, cfl_n, cff_n                    &
     &,   p_layer_centres                                               &

! Other fields passed in

! diagnostic info
     &,                                                                 &
#include "argsts.h"
     &    STASHwork4                                                    &
!
! SCM diagnostics switches (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! Increment fields passed in/out
     &,   Theta_inc, q_inc, qcl_inc, qcf_inc, cf_inc, cfl_inc, cff_inc  &

! error information
     &,   Error_code  )

      End If  ! L_pc2 and lcv_pc2_diag_sh
!
! ----------------------------------------------------------------------
! Section RAD   Radiation scheme.
!               This incorporates radiation code for non-radiation
!               timesteps which is in CLDCTL1.dk in the UM.
!-----------------------------------------------------------------------

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Radiation',5)
!  Set dimensions of mineral dust arrays for passing to rad_ctl
      IF (L_USE_DUST) THEN
        DUST_DIM1 = ROW_LENGTH*ROWS
        DUST_DIM2 = MODEL_LEVELS
      ELSE
        DUST_DIM1 = 1
        DUST_DIM2 = 1
      END IF
!  Set dimensions of _SULPHATE arrays for passing to RAD_CTL
!   (avoids wasting space if aerosol not required)
      IF (L_USE_SULPC_DIRECT .OR. L_USE_SULPC_INDIRECT_SW               &
     &                       .OR. L_USE_SULPC_INDIRECT_LW) THEN
        SULP_DIM1 = rows*row_length
        SULP_DIM2 = model_levels
      ELSE
        SULP_DIM1 = 1
        SULP_DIM2 = 1
      END IF
!  Set dimensions of soot arrays for passing to RAD_CTL
      If (L_use_soot_direct) then
        soot_dim1 = rows*row_length
        soot_dim2 = model_levels
      Else
        soot_dim1 = 1
        soot_dim2 = 1
      End If
!  Set dimensions of array for aerosol climatology for NWP
      If (n_arcl_species > 0) then
        arcl_dim1 = rows*row_length
        arcl_dim2 = model_levels
      Else
        arcl_dim1 = 1
        arcl_dim2 = 1
      End If
!  Set dimensions of biomass arrays for passing to RAD_CTL
      If (L_use_bmass_direct .or. L_use_bmass_indirect) then
        bmass_dim1 = rows*row_length
        bmass_dim2 = model_levels
      Else
        bmass_dim1 = 1
        bmass_dim2 = 1
      End If
!  Set dimensions of biogenic array for passing to RAD_CTL
      IF (L_USE_BIOGENIC) THEN
        biogenic_dim1 = rows*row_length
        biogenic_dim2 = model_levels
      ELSE
        biogenic_dim1 = 1
        biogenic_dim2 = 1
      ENDIF
!  Set dimensions of OCFF array for passing to RAD_CTL
      If (L_use_ocff_direct .or. l_use_ocff_indirect) then
        ocff_dim1 = rows*row_length
        ocff_dim2 = model_levels
      Else
        ocff_dim1 = 1
        ocff_dim2 = 1
      End If

!     Code to calculate mixing ratio of well-mixed greenhouse gases
!       if scenarios for their time variation have been prescribed.
! DEPENDS ON: timer
      If ( Ltimer ) Call timer ('gas_calc', 3)

        If ( clim_fcg_nyears(s_co2) > 0 ) Then  !  CO2 level calculated
! depends on: gas_calc
          Call gas_calc ( co2_mmr,                                      &
            clim_fcg_nyears(s_co2),  clim_fcg_years(1,s_co2),           &
            clim_fcg_levls(1,s_co2), clim_fcg_rates(1,s_co2),           &
            lenscen, error_code)
          If ( error_code /= 0 ) Return
        End if

        If ( clim_fcg_nyears(s_n2o) > 0 ) Then  !  same for N2O
! depends on: gas_calc
          Call gas_calc ( n2o_mix_ratio,                                &
             clim_fcg_nyears(s_n2o),  clim_fcg_years(1,s_n2o),          &
             clim_fcg_levls(1,s_n2o), clim_fcg_rates(1,s_n2o),          &
             lenscen, error_code)
          If ( error_code /= 0 ) Return
        End if

        If ( clim_fcg_nyears(s_ch4) > 0 ) Then  !  CH4
! depends on: gas_calc
          Call gas_calc ( ch4_mix_ratio,                                &
             clim_fcg_nyears(s_ch4),  clim_fcg_years(1,s_ch4),          &
             clim_fcg_levls(1,s_ch4), clim_fcg_rates(1,s_ch4),          &
             lenscen, error_code)
          If ( error_code /= 0 ) Return
        End if

        If ( clim_fcg_nyears(s_cfc11) > 0 ) Then   !  "CFC11"
! depends on: gas_calc
          Call gas_calc ( cfc11_mix_ratio,                              &
             clim_fcg_nyears(s_cfc11),  clim_fcg_years(1,s_cfc11),      &
             clim_fcg_levls(1,s_cfc11), clim_fcg_rates(1,s_cfc11),      &
             lenscen, error_code)
          If ( error_code /= 0 ) Return
        End if

        If ( clim_fcg_nyears(s_cfc12) > 0 ) Then   !  "CFC12"
! depends on: gas_calc
          Call gas_calc ( cfc12_mix_ratio,                              &
             clim_fcg_nyears(s_cfc12),  clim_fcg_years(1,s_cfc12),      &
             clim_fcg_levls(1,s_cfc12), clim_fcg_rates(1,s_cfc12),      &
             lenscen, error_code)
          If ( error_code /= 0 ) Return
        End if

        If ( clim_fcg_nyears(s_cfc113) > 0 ) Then  !  "CFC113"
! depends on: gas_calc
          Call gas_calc ( c113mmr,                                      &
             clim_fcg_nyears(s_cfc113),  clim_fcg_years(1,s_cfc113),    &
             clim_fcg_levls(1,s_cfc113), clim_fcg_rates(1,s_cfc113),    &
             lenscen, error_code)
          If ( error_code /= 0 ) Return
        End if

        If ( clim_fcg_nyears(s_hcfc22) > 0 ) Then  !  "HCFC22"
! depends on: gas_calc
          Call gas_calc ( hcfc22mmr,                                    &
             clim_fcg_nyears(s_hcfc22),  clim_fcg_years(1,s_hcfc22),    &
             clim_fcg_levls(1,s_hcfc22), clim_fcg_rates(1,s_hcfc22),    &
             lenscen, error_code)
          If ( error_code /= 0 ) Return
        End if

        If ( clim_fcg_nyears(s_hfc125) > 0 ) Then  !  "HFC125"
! depends on: gas_calc
          Call gas_calc ( hfc125mmr,                                    &
             clim_fcg_nyears(s_hfc125),  clim_fcg_years(1,s_hfc125),    &
             clim_fcg_levls(1,s_hfc125), clim_fcg_rates(1,s_hfc125),    &
             lenscen, error_code)
          If ( error_code /= 0 ) Return
        Endif

        If ( clim_fcg_nyears(s_hfc134a) > 0 ) Then  ! "HFC134a"
! depends on: gas_calc
          Call gas_calc ( hfc134ammr,                                   &
             clim_fcg_nyears(s_hfc134a),  clim_fcg_years(1,s_hfc134a),  &
             clim_fcg_levls(1,s_hfc134a), clim_fcg_rates(1,s_hfc134a),  &
             lenscen, error_code)
          If ( error_code /= 0 ) Return
        End if
! DEPENDS ON: timer
      If ( Ltimer ) Call timer ('gas_calc', 4)

! Copy greenhouse gases into greenhouse gas array
      ALLOCATE(grgas_field(row_length, rows, model_levels, ngrgas))

      DO i=1,ngrgas
        IF (grgas_addr(i) > 0)                                          &
     &    grgas_field(:,:,:,i) =                                        &
     &      free_tracers(1:row_length,1:rows,:,grgas_addr(i))
      END DO

      IF (grgas_addr(p_o3) <= 0)                                        &
     &  grgas_field(:,:,1:ozone_levels,p_o3) = ozone

      IF (grgas_addr(p_ch4) <= 0)                                       &
     &  grgas_field(:,:,:,p_ch4) = ch4_mix_ratio ! global constant

      IF (grgas_addr(p_n2o) <= 0)                                       &
     &  grgas_field(:,:,:,p_n2o) = n2o_mix_ratio ! global constant

      IF (grgas_addr(p_f11) <= 0) THEN
        grgas_field(:,:,:,p_f11) = cfc11_mix_ratio
! Commented out rescaling option - l_ukca_stratcfc not available
!      ELSEIF (.NOT. L_ukca_stratcfc) THEN
! rescale field so that surface mixing ratio equals cfc11_mix_ratio.
!        DO k=model_levels,1,-1
!          grgas_field(:,:,k,p_f11) = grgas_field(:,:,k,p_f11) *         &
!     &      cfc11_mix_ratio / grgas_field(:,:,1,p_f11)
!        END DO
      END IF

      IF (grgas_addr(p_f12) <= 0) THEN
        grgas_field(:,:,:,p_f12) = cfc12_mix_ratio
! Commented out rescaling option - l_ukca_stratcfc not available
!      ELSEIF (.NOT. L_ukca_stratcfc) THEN
! rescale field so that surface mixing ratio equals cfc11_mix_ratio.
!        DO k=model_levels,1,-1
!          grgas_field(:,:,k,p_f12) = grgas_field(:,:,k,p_f12) *         &
!     &      cfc12_mix_ratio / grgas_field(:,:,1,p_f12)
!        END DO
      END IF

      IF (grgas_addr(p_f113) <= 0)                                      &
     &  grgas_field(:,:,:,p_f113) = c113mmr  ! global constant

      IF (grgas_addr(p_f22) <= 0)                                       &
     &  grgas_field(:,:,:,p_f22) = hcfc22mmr  ! global constant

      grgas_field(:,:,1:wet_levels,p_h2os) =                            & 
     &  q_n(1:row_length, 1:rows, :)


! DEPENDS ON: timer
      If (Ltimer) Call timer ('RAD_CTL',3)

      If (L_radiation) then
! DEPENDS ON: ni_rad_ctl
       Call NI_rad_ctl (                                                &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &

! model dimensions.
     &, row_length, rows, n_rows                                        &
     &, model_levels, wet_levels, bl_levels                             &
     &, Ozone_levels, cloud_levels, N_cca_levels                        &
     &, NTILES, LAND_POINTS, DUST_DIM1, DUST_DIM2, biogenic_dim1        &
     &, biogenic_dim2, sulp_dim1, sulp_dim2, soot_dim1, soot_dim2       &
     &, bmass_dim1, bmass_dim2, salt_dim1, salt_dim2, salt_dim3         &
     &, co2_dim_len, co2_dim_row, co2_dim_lev, arcl_dim1, arcl_dim2     &
     &, n_arcl_species, n_arcl_compnts, i_arcl_compnts                  &
     &, ocff_dim1, ocff_dim2                                            &

! Model switches
     &, model_domain                                                    &
     &, L_Rad_Step,L_Rad_Step_diag, L_Rad_Step_prog                     &
     &, L_Forcing, L_Timestep, L_Radiance , L_Wenyi                     &
     &, L_CAL360, L_SEC_VAR, L_EqT                                      &
     &, L_MICROPHY, L_INHOM_CLOUD, L_USE_DUST, L_USE_BIOGENIC           &
     &, L_use_sulpc_direct, L_use_soot_direct, L_use_soot_indirect      &
     &, L_use_bmass_direct, L_use_bmass_indirect                        &
     &, L_use_ocff_direct, L_use_ocff_indirect                          &
     &, L_use_sulpc_indirect_SW, L_use_sulpc_indirect_LW                &
     &, L_emcorr, L_climat_aerosol, L_clim_aero_hgt, L_HadGEM1_Clim_Aero&
     &, Ltimer, L_ssice_albedo, l_cable, L_snow_albedo                  &
     &, L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a              &
     &, L_use_seasalt_indirect                                          &
     &, L_use_seasalt_direct                                            &
     &, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                           &
     &, L_pc2, l_mixing_ratio                                           &
     &, L_MURK_RAD                                                      &
     &, L_rad_deg                                                       &
! kdcorbin, 06/10 - added l_co2_radiation switch 
     &, L_co2_interactive, l_co2_radiation                              &
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
     &, rho, r_rho_levels, r_theta_levels                               &
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
     &, sin_theta_longitude, cos_theta_longitude, FV_cos_theta_latitude &

! grid-dependent arrays
     &, f3_at_u, true_longitude, true_latitude                          &

! diagnostic info
     &     ,                                                            &
#include "argsts.h"
     & STASHwork1,                                                      &
     & STASHwork2                                                       &
!
! SCM diagnostics switches (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! in data fields.
     &, p_star                                                          &
     &, p_layer_boundaries, p_layer_centres                             &
     &, p, p_theta_levels                                               &
     &, exner_rho_levels, exner_theta_levels                            &
     &, land_sea_mask, fland,land0p5,l_ctile                            &
     &, T_SURF,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE , AREA_CLOUD_FRACTION    &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, biogenic, so4_aitken, so4_accu, so4_diss                        &
     &, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld   &
     &, ocff_new, ocff_agd, ocff_cld                                    &
     &, aerosol(1:row_length, 1:rows, 1:model_levels), arcl             &
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
! chemical ozone; replaces climatological ozone
     &, grgas_field(:,:,1:ozone_levels,p_o3)                            &
     &, SW_incs, LW_incs, dirpar_local                                  &
     &, O3_trop_level, O3_trop_height                                   &
     &, T_trop_level, T_trop_height, zh, land_index, albsoil, lai       &
     &, snow_tile, tile_frac, tstar_tile, z0_tile                       &
     &, dOLR_rts, LW_down, SW_tile_rts                                  &
     &, u_1, v_1                                                        &
     &, land_alb,sice_alb                                               &
     &, ES_SPACE_INTERP,A_SW_radstep, A_LW_radstep                      &
     &, A_SW_radstep_diag, A_SW_radstep_prog                            &
     &, A_LW_radstep_diag, A_LW_radstep_prog, RAD_MASK                  &

! in/out
! chemical water replaces hydrological water
     &, T_n, grgas_field(:,:,1:wet_levels,p_h2os)                       &
     &, qcl_n, qcf_n, cf_n, cfl_n, cff_n                                &
     &, qcf2_n, qrain_n, qgraup_n                                       &
     &, Theta_inc, q_inc, qcl_inc, cf_inc, cfl_inc                      &
     &, sum_eng_fluxes                                                  &

! out.
     &, photosynth_act_rad, surf_radflux, rad_hr, dOLR, SW_tile         &
!    EAK
!    IN
     &, istep_cur                                                       &
!     &, L_TILE_PTS,TILE_PTS,SM_LEVELS,TILE_INDEX           &
     &, TILE_PTS,SM_LEVELS,TILE_INDEX                                   &
     &,   surf_down_sw,alb_tile                                         &
     &, SNOW_TMP3L,SNOW_RHO1L,TSOIL_TILE,ISNOW_FLG3L,LAND_ALBEDO        &
! Section Information
     &, maxsects, h_sect                                                &

! error information
     &, Error_code  )
      else
        rad_hr(:,:,:,:)=0.0
        dolr(:,:)=0.0
        sw_tile(:,:)=0.0
      endif

! DEPENDS ON: timer
      If (Ltimer) Call timer ('RAD_CTL',4)
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Radiation',6)
      DEALLOCATE(grgas_field)

!-----------------------------------------------------------------------
! Set dirpar_inc dependent prognostic for output to D1 array:
!-----------------------------------------------------------------------

! Use prognostic dirpar only if Stochem is switched on.

      IF ( H_SECT(25) == "01A" ) THEN
         L_USE_DIRPAR=.TRUE.
      ELSE
         L_USE_DIRPAR=.FALSE.
      ENDIF              
         
      IF (L_USE_DIRPAR) THEN     
         DO J=1, ROWS
           DO I=1, ROW_LENGTH
              DIRPAR_INC(I,J)=DIRPAR_LOCAL(I,J)
           ENDDO
         ENDDO  
      ENDIF      

! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics for output to D1 array:
! ----------------------------------------------------------------------
      IF(L_CTILE)THEN
        DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            LAND_ALB_CTILE(I,J)=LAND_ALB(I,J)
            SICE_ALB_CTILE(I,J)=SICE_ALB(I,J)
          ENDDO
        ENDDO
      ENDIF

      If (Error_code  ==  0) Then

! Convert temperature held in t_n to potential temperature.
! Convert temperature increment to theta increment

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              T_n(i,j,k) = Theta(i,j,k)
              Theta_inc(i,j,k) = Theta_inc(i,j,k) /                     &
     &                           exner_theta_levels(i,j,k)
            End Do
          End Do
        End Do

      End If ! on error code equal to zero

! ----------------------------------------------------------------------
! Section CNV.2 Energy correction code
! ----------------------------------------------------------------------

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Conv Eng Corr',5)
      If (L_emcorr .and. Error_code  ==  0) Then

! Add convective + large scale, rain and snow, at the surface to the
! diabatic heating for use in the energy correction
! procedure.
! Scale variables by conversion factor so that only one call is required

        lclf = lc + lf
        Do j = 1, rows
          Do i = 1, row_length
            tot_precip_scaled(i,j) = ls_rain(i,j) * lc +                &
     &                               ls_snow(i,j) * lclf
          End Do
        End Do

! DEPENDS ON: timer
      If (Ltimer) Call timer ('flux_diag',3)

! DEPENDS ON: flux_diag
        Call flux_diag(tot_precip_scaled, FV_cos_theta_latitude,        &
     &                 row_length, rows ,off_x, off_y, 1.0,             &
     &                  sum_eng_fluxes,timestep)
! DEPENDS ON: timer
      If (Ltimer) Call timer ('flux_diag',4)
! moist fluxes
        Do j = 1, rows
          Do i = 1, row_length
            tot_precip_scaled(i,j) = -ls_rain(i,j) -ls_snow(i,j)
          End Do
        End Do

! DEPENDS ON: timer
      If (Ltimer) Call timer ('flux_diag',3)
! DEPENDS ON: flux_diag
        Call flux_diag(tot_precip_scaled, FV_cos_theta_latitude,        &
     &                 row_length, rows ,off_x, off_y, 1.0,             &
     &                  sum_moist_flux,timestep)
! DEPENDS ON: timer
      If (Ltimer) Call timer ('flux_diag',4)


      End If   ! L_emcorr

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Conv Eng Corr',6)

! ----------------------------------------------------------------------
! Section GWD.1 Call gravity wave drag
! l_gwd      = orographic GWD
! l_use_ussp = middle atmosphere non-orographic GWD scheme
! ----------------------------------------------------------------------

      If ( l_gwd .or. l_use_ussp ) Then

        If (error_code  ==  0 ) Then
! DEPENDS ON: timer
          If (Ltimer) Call timer ('AP1 G-wave drag',5)

! reset p at layer boundaries for spectral GWD code
          Do k = 1, model_levels - 1
            Do j = 1, rows
              Do i = 1, row_length
                p_layer_boundaries(i,j,k) = p_theta_levels(i,j,k)
              End Do
            End Do
          End Do
          Do j = 1, rows
            Do i = 1, row_length
              p_layer_boundaries(i,j,model_levels) = 0.0
            End Do
          End Do


! DEPENDS ON: timer
        If (Ltimer) Call timer ('NI_gwd_ctl',3)

! DEPENDS ON: ni_gwd_ctl
        Call NI_gwd_ctl(                                                &
     &                   halo_i, halo_j, off_x, off_y                   &
     &,                  global_row_length, global_rows                 &
     &,                  proc_row_group, at_extremity, n_proc, n_procx  &
     &,                  n_procy, neighbour, g_rows, g_row_length       &
     &,                  g_datastart, me                                &

! model dimensions.
     &,                  row_length, rows, n_rows, land_points          &
     &,                  model_levels                                   &

! Model switches
     &,                  model_domain, gw_kay, gwd_frc                  &
     &,                  Ltimer, l_gwd, L_use_ussp                      &
     &,                  L_taus_scale, L_fix_gwsatn, l_gwd_40km         &
     &,                  L_ussp_opaque                                  &
     &,                  sat_scheme, gwd_fsat                           &

! trig arrays
     &,                  sin_theta_longitude, sin_theta_latitude        &

! in coordinate information
     &,                  r_rho_levels, r_theta_levels                   &
     &,                  eta_theta_levels, eta_rho_levels               &

! in time stepping information.
     &,                  timestep, timestep_number                      &

! diagnostic info
     &     ,                                                            &
#include "argsts.h"
     & STASHwork6                                                       &

! in data fields.
     &,                  u, v, land_sea_mask                            &
     &,                  p_layer_boundaries                             &
     &,                  rho, t_n, sd_orog_land                         &
     &,                  orog_grad_xx_land, orog_grad_xy_land           &
     &,                  orog_grad_yy_land, land_index                  &

! in/out
     &,                  u_inc, v_inc                                   &

! error information
     &, Error_code  )

! DEPENDS ON: timer
          If (Ltimer) Call timer ('NI_gwd_ctl',4)
! DEPENDS ON: timer
          If (Ltimer) Call timer ('AP1 G-wave drag',6)
        End If

      End If

!-----------------------------------------------------------------------
! Tracer Source and Boundary updating where appliable
!-----------------------------------------------------------------------
! Aerosol
      If ( L_murk_source ) Then

        Do level = 1, model_levels
! DEPENDS ON: trsrce
          Call trsrce(                                                  &
     &       rows, row_length, off_x, off_y                             &
     &,      halo_i, halo_j, model_levels, wet_levels                   &
     &, 0, 0                                                            &
     &,      r_rho_levels, r_theta_levels                               &
     &,      theta, q_n, qcl_n, qcf_n, exner_rho_levels, rho            &
     &,      aerosol(:,:,level), aerosol_em (:,:,level )                &
     &,      level, timestep, val_hour, val_minute, amp                 &
     &)
        End Do
      End If

      If ( L_murk_bdry ) Then
! DEPENDS ON: trbdry
        Call trbdry(                                                    &
     &       row_length, rows, n_rows, model_levels                     &
     &,      off_x, off_y, at_extremity                                 &
     &,      p, u, v                                                    &
     &,      aerosol, timestep                                          &
     &)
      End If

      If ( L_murk ) Then
        If (.NOT. L_murk_bdry) Then
          ! Set the external halos if not specified by the
          ! formula for UK mes.

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &         aerosol, row_length, rows,                               &
     &         model_levels, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
          Call Fill_External_Halos(                                     &
     &               aerosol,      row_length, rows,                    &
     &               model_levels, off_x,      off_y )

        End If

      End If

!---------------------------------------------------------------------
! Section METHOX Call methane oxidation
!---------------------------------------------------------------------

      If (L_USE_METHOX) Then

        If (error_code == 0) Then
! DEPENDS ON: timer
          If (Ltimer) Call timer ('NI_methox',5)
!
!    call methane oxidation directly (no interface routine needed)
!
! DEPENDS ON: ni_methox
          call NI_methox(                                               &
! Parallel variables
     &  halo_i, halo_j                                                  &

! model dimensions.
     &, row_length, rows                                                &
     &, model_levels, wet_levels                                        &

! model levels
     &, eta_theta_levels                                                &

! Model switches
     &, Ltimer                                                          &

! in time stepping information.
     &, timestep                                                        &

! in/out
     &,  q_n,q_inc                                                      &

! error information
     &, Error_code  )

! DEPENDS ON: timer
          If (Ltimer) Call timer ('NI_methox',6)
        End If
      End If

! ----------------------------------------------------------------------
! Copy increment arrays into star locations
! ----------------------------------------------------------------------

! In the LAM set physics increments on boundary to zero.
      If (model_domain  ==  mt_LAM) Then

        L_zero_boundaries=.TRUE.
        L_zero_halos=.FALSE.

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,MODEL_LEVELS,fld_type_p,THETA_INC,         &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,Q_INC,               &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,QCL_INC,             &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,QCF_INC,             &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

        If (L_mcr_qcf2)                                                 &
                        ! prognostic second cloud ice in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qcf2_inc,           &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

        If (L_mcr_qrain)                                                &
                         ! prognostic rain in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qrain_inc,          &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

        If (L_mcr_qgraup)                                               &
                          ! prognostic graupel in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qgraup_inc,         &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,CF_INC,              &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,CFL_INC,             &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,CFF_INC,             &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,n_ROWS,off_x,off_y,MODEL_LEVELS,fld_type_v,V_INC,   &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,off_x,off_y,MODEL_LEVELS,fld_type_u,U_INC,     &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

      End If

      If (L_pc2) then
!
#if defined(SCMA)
      Allocate( T_inc_diag  (row_length,rows,model_levels) )
      Allocate( q_inc_diag  (row_length,rows,wet_levels)   )
      Allocate( qcl_inc_diag(row_length,rows,wet_levels)   )
      Allocate( qcf_inc_diag(row_length,rows,wet_levels)   )
      Allocate( cf_inc_diag (row_length,rows,wet_levels)   )
      Allocate( cfl_inc_diag(row_length,rows,wet_levels)   )
      Allocate( cff_inc_diag(row_length,rows,wet_levels)   )
#endif
!
! Now call a consistency check for the moisture and cloud fields.
! The catch is that _inc variables hold the increments and not the
! full variables yet, so temporarily form them as _inc, check them,
! then rewrite them as increments.
!
! 1. Create full variables (not increments). Only need to do this
!    on wet_levels for theta_inc (theta_inc and T_n hold potential
!    temperature) since the checking is only over wet levels.
!
      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_inc(i,j,k) = (t_n(i,j,k)   + theta_inc(i,j,k))        &
     &                         * exner_theta_levels(i,j,k)
! theta_inc now holds temperature, t_n still holds potential temp
            q_inc(i,j,k)     = q_n(i,j,k)   + q_inc(i,j,k)
            qcl_inc(i,j,k)   = qcl_n(i,j,k) + qcl_inc(i,j,k)
            qcf_inc(i,j,k)   = qcf_n(i,j,k) + qcf_inc(i,j,k)
            cf_inc(i,j,k)    = cf_n(i,j,k)  + cf_inc(i,j,k)
            cfl_inc(i,j,k)   = cfl_n(i,j,k) + cfl_inc(i,j,k)
            cff_inc(i,j,k)   = cff_n(i,j,k) + cff_inc(i,j,k)
#if defined(SCMA)
            t_inc_diag(i,j,k)   = theta_inc(i,j,k)
            q_inc_diag(i,j,k)   = q_inc(i,j,k)
            qcl_inc_diag(i,j,k) = qcl_inc(i,j,k)
            qcf_inc_diag(i,j,k) = qcf_inc(i,j,k)
            cf_inc_diag(i,j,k)  = cf_inc(i,j,k)
            cfl_inc_diag(i,j,k) = cfl_inc(i,j,k)
            cff_inc_diag(i,j,k) = cff_inc(i,j,k)
#endif
          End Do
        End Do
      End Do
!
! 2. Call the checking routine (needs to use temperature, not
!    potential temperature)
!
! DEPENDS ON: pc2_checks
      call pc2_checks(p_layer_centres(1,1,1),wet_levels                 &
     &,               row_length, rows, theta_inc, cf_inc, cfl_inc      &
     &,               cff_inc, q_inc, qcl_inc, qcf_inc                  &
     &,               l_mixing_ratio)
!
#if defined(SCMA)
! Now form the increment diagnostics in t_inc_diag etc.
!
      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            t_inc_diag(i,j,k)   = theta_inc(i,j,k) - t_inc_diag(i,j,k)
            q_inc_diag(i,j,k)   = q_inc(i,j,k)     - q_inc_diag(i,j,k)
            qcl_inc_diag(i,j,k) = qcl_inc(i,j,k)   - qcl_inc_diag(i,j,k)
            qcf_inc_diag(i,j,k) = qcf_inc(i,j,k)   - qcf_inc_diag(i,j,k)
            cf_inc_diag(i,j,k)  = cf_inc(i,j,k)    - cf_inc_diag(i,j,k)
            cfl_inc_diag(i,j,k) = cfl_inc(i,j,k)   - cfl_inc_diag(i,j,k)
            cff_inc_diag(i,j,k) = cff_inc(i,j,k)   - cff_inc_diag(i,j,k)
          End Do ! i
        End Do ! j
      End Do ! k
#endif
!
! 3. Update fields to produce the net increments which are written to
!    the _star variables. Note that theta_inc
!    currently holds the full temperature value up to wet_levels, and
!    the increment in potential temp from there to model_levels
!
      Do k = wet_levels+1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_star(i,j,k) = theta_inc(i,j,k)
#if defined(SCMA)
            ! No increment occurs above wet_levels
            t_inc_diag(i,j,k) = 0.0
#endif
          End Do
        End Do
      End Do
!
#if defined(SCMA)

!-----------------------------------------------------------------------
!     SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_PC2)) Then

!         Stash 4,141
! DEPENDS ON: scmoutput
          Call SCMoutput(T_inc_diag,                                    &
               'dt_pc2ck','Temperature increment PC2 checks','K',       &
               t_avg,d_all,default_streams,'',RoutineName)

!         Stash 4,142
! DEPENDS ON: scmoutput
          Call SCMoutput(q_inc_diag,                                    &
               'dq_pc2ck','Specific humidity increment PC2 checks',     &
               'kg/kg',                                                 &
               t_avg,d_wet,default_streams,'',RoutineName)

!         Stash 4,143
! DEPENDS ON: scmoutput
          Call SCMoutput(qcl_inc_diag,                                  &
               'dqcl_pc2ck','QCL increment PC2 checks','kg/kg',         &
               t_avg,d_wet,default_streams,'',RoutineName)

!         Stash 4,144
! DEPENDS ON: scmoutput
          Call SCMoutput(qcf_inc_diag,                                  &
               'dqcf_pc2ck','QCF increment PC2 checks','kg/kg',         &
               t_avg,d_wet,default_streams,'',RoutineName)

!         Stash 4,152
! DEPENDS ON: scmoutput
          Call SCMoutput(cf_inc_diag,                                   &
               'dbcf_pc2ck','Bulk cloud fraction increment PC2 checks', &
               'Fraction',                                              &
               t_avg,d_wet,default_streams,'',RoutineName)

!         Stash 4,153
! DEPENDS ON: scmoutput
          Call SCMoutput(cfl_inc_diag,                                  &
               'dcfl_pc2ck',                                            &
               'Liquid cloud fraction increment PC2 checks','Fraction', &
               t_avg,d_wet,default_streams,'',RoutineName)

!         Stash 4,154
! DEPENDS ON: scmoutput
          Call SCMoutput(cff_inc_diag,                                  &
               'dcff_pc2ck',                                            &
               'Ice cloud fraction increment PC2 checks','Fraction',    &
               t_avg,d_wet,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_PC2)

      Deallocate( T_inc_diag   )
      Deallocate( q_inc_diag   )
      Deallocate( qcl_inc_diag )
      Deallocate( qcf_inc_diag )
      Deallocate( cf_inc_diag  )
      Deallocate( cfl_inc_diag )
      Deallocate( cff_inc_diag )

#endif
!
! 6. Update fields to produce the net increments which are written to
!    the _star variables.

      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_star(i,j,k) = theta_inc(i,j,k)                        &
     &                        /exner_theta_levels(i,j,k) - t_n(i,j,k)
! theta_star now holds the potential temperature increment
            q_star(i,j,k)     = q_inc(i,j,k)     - q_n(i,j,k)
            qcl_star(i,j,k)   = qcl_inc(i,j,k)   - qcl_n(i,j,k)
            qcf_star(i,j,k)   = qcf_inc(i,j,k)   - qcf_n(i,j,k)
            cf_star(i,j,k)    = cf_inc(i,j,k)    - cf_n(i,j,k)
            cfl_star(i,j,k)   = cfl_inc(i,j,k)   - cfl_n(i,j,k)
            cff_star(i,j,k)   = cff_inc(i,j,k)   - cff_n(i,j,k)
          End Do
        End Do
      End Do

      Else  ! L_pc2

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_star(i,j,k) = theta_inc(i,j,k)
          End Do
        End Do
      End Do

      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            q_star  (i,j,k) = q_inc  (i,j,k)
            qcl_star(i,j,k) = qcl_inc(i,j,k)
            qcf_star(i,j,k) = qcf_inc(i,j,k)
            cf_star (i,j,k) = cf_inc (i,j,k)
            cfl_star(i,j,k) = cfl_inc(i,j,k)
            cff_star(i,j,k) = cff_inc(i,j,k)
          End Do
        End Do
      End Do

      End If  ! L_pc2

      If (L_mcr_qcf2)                                                   &
                      ! prognostic second cloud ice in use
     &  qcf2_star(1:row_length, 1:rows, :) = qcf2_inc(:,:,:)

      If (L_mcr_qrain)                                                  &
                       ! prognostic rain in use
     &  qrain_star(1:row_length, 1:rows, :) = qrain_inc(:,:,:)

      If (L_mcr_qgraup)                                                 &
                        ! prognostic graupel in use
     &  qgraup_star(1:row_length, 1:rows, :) = qgraup_inc(:,:,:)

      ! Deallocate additional microphysics variables
      Deallocate ( qcf2_inc )
      Deallocate ( qrain_inc )
      Deallocate ( qgraup_inc )
      Deallocate ( qcf2_n )
      Deallocate ( qrain_n )
      Deallocate ( qgraup_n )

! end of routine Atmos_physics1
      RETURN
      END SUBROUTINE Atmos_Physics1
#endif
#endif
