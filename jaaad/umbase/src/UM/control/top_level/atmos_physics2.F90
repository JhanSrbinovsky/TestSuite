#if defined(ATMOS)
#if defined(CONTROL) || defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interface to Atmos Physics parametrizations after S-L advection.
!
! Subroutine Interface:
      SUBROUTINE Atmos_Physics2(                                        &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &
     &, NumCycles, CycleNo                                              &

! model dimensions.
     &, row_length, rows, n_rows, land_points, model_levels, nice       &
     &, wet_levels, bl_levels, dst_levels, dsm_levels, cloud_levels     &
     &, land_ice_points, soil_points, n_cca_levels, ntiles, tr_levels   &
     &, first_constant_r_rho_level, DIM_CS1, DIM_CS2                    &

! IN Substepping information and switches
     &, Num_Substeps, L_phys2_substep                                   &
! Model switches
     &, L_regular, l_mixing_ratio, model_domain, L_dry                  &
     &, FORMDRAG, OROG_DRAG_PARAM, L_CAL360, L_emcorr, Ltimer           &
     &, L_us_blsol, BL_OPTIONS, L_area_cloud, L_ACF_Cusack              &
     &, L_ACF_Brooks, L_RHCPT, L_hydrology, L_bl                        &
     &, L_Murk, L_murk_advect, L_murk_source, L_bl_tracer_mix, L_DUST   &
     &, L_CAM_DUST,L_sulpc_so2,L_sulpc_nh3,L_SBLeq,L_SBLco,L_sulpc_dms  &
     &, L_soot, L_biomass, L_ocff, L_co2_interactive, L_ctile           &
     &, L_co2_emits, L_use_bl_diag_term, L_pc2, L_pc2_reset, L_eacf     &
     &, L_LAMBDAM2, L_FULL_LAMBDAS, L_UKCA, L_sice_heatflux, L_OHADGEM1 &
     &, L_USE_CARIOLLE, l_anthrop_heat_src                              &

! Model Parameters
     &, alpha_cd, Puns, Pstb, rhcrit, CO2_MMR, tr_vars, tr_ukca         &
     &, Muw_SBL, Mwt_SBL, cloud_fraction_method, overlap_ice_liquid     &
     &, ice_fraction_method,ctt_weight,t_weight,qsat_fixed,sub_cld      &
     &, x1i, x1ic, x1r, x2r, x4r, l_psd, ai, bi, aic, bic, lsp_ei       & 
     &, lsp_fi, lsp_eic, lsp_fic                                        &  
     &, dbsdtbs_turb_0, Charnock, SeaSalinityFactor                     &

! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &

! in coordinate information
     &, r_rho_levels, r_theta_levels, r_at_u, r_at_v, unscaled_dry_rho  &
     &, eta_theta_levels, eta_rho_levels, delta_lambda, delta_phi       &
     &, dlambda_p, dphi_p, wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v &
     &, lat_rot_NP, long_rot_NP , f3_at_u                               &
! Variables required by STPH_SCV
     &, mix_v_phys2, mix_cl_phys2, mix_cf_phys2                         &

! in time stepping information.
     &, timestep, val_year, val_day_number, val_hour, val_minute        &
     &, val_second, timestep_number                                     &

! trig arrays
     &, sin_theta_longitude, cos_theta_longitude, sin_theta_latitude,   &
!     &, sin_theta_longitude, cos_theta_longitude, FV_cos_theta_latitude,&
     &  FV_cos_theta_latitude, sec_theta_latitude,                       &
!     &  sec_theta_latitude,                                             &

! River routing
     & AOCPL_ROW_LENGTH,AOCPL_P_ROWS,L_RIVERS,XPA,XUA,XVA,YPA,YUA,YVA,  &
     & G_P_FIELD, G_R_FIELD, A_STEPS_SINCE_RIV, RIVER_ROW_LENGTH,       &
     & RIVER_ROWS, global_river_row_length, global_river_rows,          &
     & RIVER_STEP, RIVER_VEL, RIVER_MCOEF,                              &
! Add inland basin outflow to arguments
     & TRIVDIR, TRIVSEQ, TWATSTOR,INLANDOUT_ATM,                        &
     & L_INLAND,                                                        &


! Grid-to-grid river routing
     &  R_AREA, SLOPE, FLOWOBS1, R_INEXT, R_JNEXT, R_LAND,              &
     &  SUBSTORE, SURFSTORE, FLOWIN, BFLOWIN,                           &
!
! diagnostic info
#include "argsts.h"
     & STASHwork3,STASHwork5,STASHwork8,STASHwork9,STASHwork19,         &
     & STASHwork26                                                      &
!
! SCM Diagnostics (dummy values in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &
!
! in data fields.
     &, theta, q, qcl, qcf, rho_rsq, u, v, w, w_adv, p, p_star          &
     &, exner_rho_levels, exner_theta_levels                            &
     &, land_sea_mask, p_theta_levels                                   &
! variables for subgrid turbulence scheme
     &, visc_BL_m, FM_3D, FH_3D, L_subfilter_vert, L_subfilter_horiz    &
     &, L_subfilter_blend, max_diff,turb_startlev_vert,turb_endlev_vert &
     &, BL_COEF_KM, BL_COEF_KH                                          &

! ancillary fields and fields needed to be kept from timestep to
! timestep

     &, land_index, land_ice_index, soil_index, canopy_gb, snow_depth   &
     &, hcon, hcap, smvccl, smvcwt, smvcst, sthf, sthu                  &
     &, sil_orog_land, ho2r2_orog, di, ice_fract, u_0, v_0, u_0_p       &
     &, v_0_p, cca, ccb, cct, cclwp, ccw_out, lcbase_out, t_soil, ti    &
     &, ti_gb                                                           &
     &, t_surf, z0msea, ice_fract_ncat,di_ncat,satcon,sathh,clapp       &
     &, soil_layer_moisture, t1_sd, q1_sd, zh                           &
     &, area_cloud_fraction, bulk_cloud_fraction_halos                  &
     &, cloud_fraction_liquid_halos, cloud_fraction_frozen_halos        &
     &, ls_rain, ls_snow, micro_tends                                   &
     &, photosynth_act_rad, rad_hr, surf_radflux                        &
     &, SOIL_CLAY,SOIL_SILT,SOIL_SAND,DUST_MREL1,DUST_MREL2,DUST_MREL3  &
     &, DUST_MREL4,DUST_MREL5,DUST_MREL6                                &
     &, so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em          &
     &, ocff_hilem, ocff_em, co2_emits, co2flux                         &

! Tracer Fluxes - kdcorbin, 05/10
     &, tracer_flux1, tracer_flux2, tracer_flux3, tracer_flux4          &
     &, tracer_flux5, tracer_flux6, tracer_flux7, tracer_flux8          &
     &, tracer_flux9, tracer_flux10, tracer_flux11, tracer_flux12       &
     &, tracer_flux13, tracer_flux14, tracer_flux15, tracer_flux16      &
     &, tracer_flux17, tracer_flux18, tracer_flux19, tracer_flux20      &
! CO2 Global Emissions - kdcorbin, 05/10
     &, co2emitmass                                                     &   
! rml 1/7/13 flag for co2 flux into passive tracer
     &, L_CO2_TRACER                                                    &

! in/out
     &, theta_star, q_star, qcl_star, qcf_star, cf_star, cfl_star       &
     &, cff_star, R_u, R_v, R_w, sum_eng_fluxes, sum_moist_flux         &

! In/Out tracer fields
     &, aerosol, free_tracers, ukca_tracers                             &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &

     &, so2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new         &
     &, soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld           &
     &, ocff_new, ocff_aged, ocff_cld, co2                              &
! IN/OUT STPH_RP
     &,G0_RP,par_mezcla                                                 &

! IN/OUT RIVERS
     &, TOT_SURF_RUNOFF, TOT_SUB_RUNOFF                                 &
!
! out fields
     &, rhokm, cH_term, ntml, cumulus, nbdsc, ntdsc ,rhcpt              &
     &, rhc_row_length, rhc_rows                                        &

! Additional variables for MOSES II
     &, frac, frac_disturb, canht_ft, lai_ft, canopy, catch, catch_snow &
     &, snow_grnd, snow_tile, z0_tile, tstar_tile, infil_tile, rgrain  &
     &, cs, gs, co2_dim_row, co2_dim_len, l_neg_tstar, l_snow_albedo    &
     &, l_phenol, l_triffid, l_trif_eq, l_q10, asteps_since_triffid     &
     &, a_step, phenol_period, triffid_period, can_model                &
     &, g_leaf_acc, g_leaf_phen_acc, npp_ft_acc, resp_w_ft_acc          &
     &, resp_s_acc, land_pts_trif, npft_trif, olr, lw_down, sw_tile     &

     &, FLAND_CTILE,TSTAR_LAND_CTILE,TSTAR_SEA_CTILE,TSTAR_SICE_CTILE   &
     &, ALBSOIL,COS_ZENITH_ANGLE,CAN_RAD_MOD,ILAYERS                    &
     &,RADNET_TILE                                                      &

!    EAK
!    IN
     &, surf_down_sw,alb_tile,l_tile_pts                                &
!     &, surf_down_sw,alb_tile,cos_zenith_angle               &
     &, lat,long,day,time_sec,SW_DOWN                                   &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L                   &
     &, SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE             &
     &, TSOIL_TILE,T_SURF_TILE,HCONS                                    &
     &, SOIL_TYPE,VEG_TYPE                                              &
     &, ISNOW_FLG3L,total_nsteps                                        &
     &, FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB                         &
     &, TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB                     &
     &, USTAR_CAB,SURF_HTF_CAB                                          &
! Lestevens - needed for wblake fix
     &, l_cable,wblake_ratio,WB_LAKE                                    &
     &, TOT_WBLAKE,TOT_SUBRUN                                           &
!sxy
     &, TOT_ALB                                                         &      
     &, U_S_CAB,CH_CAB,CD_CAB                                           &
     &, CD,CH                                                           &
     &, TILE_PTS,TILE_INDEX                                             &
     &, SNAGE_TILE,RTSOIL_TILE                                          &
     &, GFLUX_TILE,SGFLUX_TILE                                          &
     &, TRANSP_TILE                                                     &
     &, CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER                     &
     &, NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE                           &

! Additional variables required for large-scale hydrology:
     &, L_TOP,L_PDM,FEXP,GAMTOT,TI_MEAN,TI_SIG,FSAT,FWETL,ZW            &
     &, STHZW,A_FSAT,C_FSAT,A_FWET,C_FWET,L_SOIL_SAT_DOWN               &
! Cariolle ozone and associated ancillaries
     &, OZONE_TRACER                                                    &

! Variables for SCM and idealised UM
     &, L_flux_bc,flux_e,flux_h,L_spec_z0,z0m_scm,z0h_scm               &
!
#if defined(SCMA)
! Additional variable for the SCM diagnostics
     &, nfor, ntrop, time_string, daycount, layer_cloud                 &
     &, tstar, nconvars, ntrad, ntrad1, conv_mode                       &
#endif

! error information
     &, Error_code,                                                     &
      ! end step of experiment
     & endstep, mype )

      Use cv_cntl_mod, Only:                                            &
          lcv_3d_cca, lcv_ccrad, lcv_pc2_diag_sh

      Use cv_run_mod,  Only:                                            &
          cape_opt, cape_bottom, cape_top, convection_option

      Use bl_diags_mod, Only:                                           &
          strnewbldiag

#if defined(ACCESS)
      USE auscom_cpl_data_mod,                                          &
     &    Only : auscom_salinity, access_tfs, ocn_sss
#endif

      IMPLICIT NONE
!
! Description: This version interfaces to physics schemes in sequence:
!    convection                     (optional)
!    boundary layer
!    convection                     (optional)
!    hydrology
!    river routing                  (optional)
!
!          CALLed after Semi-Lagrangian in atmosphere timestep.
! Method:
!
!
! Current Code Owner: 
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
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
#if defined(SCMA)
! Include parameters necessary for calls to SCMoutput...
#include "s_scmop.h"
#endif
! Need INVERT_OCEAN (hence N->S configuration of ATMOS grid)
! for river routing
#include "typocpar.h"

#include "c_dust_ndiv.h"


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
     &, me                                                              &
                   ! My processor number
     &, NumCycles                                                       &
                   ! Number of cycles (iterations) for iterative SISL
     &, CycleNo                                                         &
                   ! Sweep number
     &, Num_Substeps ! Substepping Info

      Logical                                                           &
     &  L_phys2_substep                                                 &
                        ! false: default BL scheme without substepping.
                        ! true: (i) BL scheme exchange coefficients
                        !        and fluxes are calculated using
                        !        latest (modified by phys1 + dyn)
                        !        values (ii) phys2 substepping can
                        !        be used.
     &, L_us_blsol      ! true: use uncond stable numerical BL solver
                        ! false: use standard implicit solver

      Real :: Puns, Pstb ! parameters for uncond stable numerical solver
                         ! Puns : used in an unstable BL column
                         ! Pstb : used in an stable BL column
      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north, sout
                         ! east or west of the processor grid

! Parameters

! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, land_points                                                     &
                    ! IN No.of land points being processed, can be 0.
     &, model_levels                                                    &
     &, wet_levels                                                      &
     &, bl_levels                                                       &
     &, dst_levels                                                      &
                    ! number of deep soil temperature levels
     &, dsm_levels                                                      &
                    ! number of deep soil moisture levels
     &, cloud_levels                                                    &
     &, land_ice_points                                                 &
                        ! number of land ice points
     &, soil_points                                                     &
                        ! number of soil points
     &, n_cca_levels                                                    &
                      ! Number of levels for conv cloud
                      ! amount: 1 for 2D, nlevs for 3D.
     &, ntiles                                                          &
                      ! No. of land-surface tiles ( MOSES II )
     &, tr_levels                                                       &
                      ! No. of free tracer levels
     &, first_constant_r_rho_level                                      &
                                   ! 1st rho level on which r constant
     &, nice                                                            &
                      ! No. of sea ice catagories
     &, DIM_CS1, DIM_CS2  ! soil carbon dimensions

! Model switches
      Integer                                                           &
     &  model_domain

      INTEGER                                                           &
     & FORMDRAG     ! IN switch for orographic form drag
      Logical                                                           &
     &  L_SBLeq                                                         &
                    ! IN Switch for Equilibrium SBL model

     &, L_SBLco                                                         &
                    ! IN Switch for coupled gradient method in
!                   !    Equilibrium SBL model
     &, L_regular                                                       &
                     ! True if NOT variable resolution
     &, Ltimer                                                          &
                 ! true then output some timing information
     &, L_dry       ! true if model to be run with no moisture

      INTEGER, DIMENSION(20) :: BL_OPTIONS   ! IN BL switches

      INTEGER trweights1
                           ! local switch to use implicit weights
                           ! of 1 for tracers


      Logical                                                           &
     &  L_CAL360                                                        &
                        ! true if using 360 day calender
     &, L_emcorr                                                        &
                        ! true if energy correction scheme used
     &, L_area_cloud                                                    &
                          ! true if using area cloud fraction (ACF)
     &, L_ACF_Cusack                                                    &
                          ! ... and selected Cusack and PC2 off
     &, L_ACF_Brooks                                                    &
                          ! ... and selected Brooks
     &, L_eacf                                                          &
     &, L_RHCPT                                                         &
                    ! Switch for 3D RHcrit Diagnostic, not 1D parameter
     &, L_Murk                                                          &
                    ! aerosol field
     &, L_murk_advect                                                   &
                          ! Switch for advecting aerosol
     &, L_murk_source                                                   &
                          ! Switch for source/scavenging aerosol
     &, L_bl_tracer_mix                                                 &
                          ! Switch for BL mixing of free tracers
     &, L_DUST                                                          &
                          ! switch for mineral dust
     &, L_CAM_DUST                                                      &
               !Old version of dust_uplift scheme used in CAM NWP models
     &, L_sulpc_so2                                                     &
                          ! Switch for Sulphur Cycle
     &, L_sulpc_nh3                                                     &
                          ! NH3 included in Sulphur Cycle
     &, L_sulpc_dms                                                     &
                          ! DMS included in Sulphur Cycle
     &, L_soot                                                          &
                          ! Switch for Soot Cycle
     &, L_biomass                                                       &
                          ! Switch for Biomass aerosol
     &, L_ocff                                                          &
                          ! Switch for fossil-fuel organic carbon aer.
     &, L_co2_interactive                                               &
                          ! Switch for interactive CO2
     &, L_pc2                                                           &
                          ! Activate PC2 cloud scheme
     &, L_pc2_reset                                                     &
                          ! Run PC2 cloud scheme in diagnostic mode
     &, L_psd                                                           &
                          ! Use generic ice particle size distribution
     &, L_mixing_ratio                                                  &
                          ! Use mixing ratios (if available)
     &, L_ctile                                                         &
                          ! True if coastal tiling selected.
     &, L_co2_emits                                                     &
                          ! Switch to include surface emissions
     &, L_CO2_TRACER                                                    &
              ! rml, 1/7/13 Switch to put co2 flux into passive tracer 
     &, L_use_bl_diag_term                                              &
                           ! Switch to calc ch_term in ni_imp_ctl()
     &, L_hydrology                                                     &
                        ! switch to turn off hydrology
     &, L_bl                                                            &
                        ! switch to turn off boundary layer scheme
     &,L_LAMBDAM2                                                       &
                   ! IN LambdaM=2*LambdaH (operational setting).
     &,L_FULL_LAMBDAS                                                   &
                      ! Lambdas NOT reduced above NTML_LOCAL+1
     &, L_ukca                                                          &
                          ! switch to turn on UKCA
     &, L_sice_heatflux                                                 &
                           ! Switch for semi-impl sea-ice temperature
     &, L_OHADGEM1                                                      &
                           ! HadGEM1 ocean sci. Here for coupled runs.
     &, L_USE_CARIOLLE     ! switch to turn off cariolle scheme for ozone

!
! Boundary Layer
      REAL, Intent(IN) :: SeaSalinityFactor
!                         ! Factor allowing for the effect
!                         ! of the salinity of sea water on
!                         ! the evaporative flux.
!
      REAL, Intent(IN) :: OROG_DRAG_PARAM
!                         ! Drag coefficient for orographic form drag

                          
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
     &  timestep

      Real                                                              &
     &  alpha_cd(bl_levels)                                             &
     &, Muw_SBL,Mwt_SBL                                                 &
                        ! IN Powers to use in prescription of equilib
!                       !    profiles of stress and buoyancy flux in
!                       !    Equilibrium SBL model
     &, rhcrit(wet_levels)                                              &
                            ! IN Critical relative humidity.
                                  ! the values need to be tuned
                                  ! for the given set of levels.
     &, CO2_MMR                                                         &
                        ! set equal to co2_start
     &, dbsdtbs_turb_0
                        ! IN PC2 erosion rate / s-1

      Real                                                              &
                  !, INTENT(IN)
     &  x1i                                                             &
                  ! Intercept of aggregate size distribution
     &, x1ic                                                            &
                  ! Intercept of crystal size distribution
     &, x1r                                                             &
                  ! Intercept of raindrop size distribution
     &, x2r                                                             &
                  ! Scaling parameter of raindrop size distribution
     &, x4r                                                             &
                  ! Shape parameter of raindrop size distribution
     &, ai, bi, aic, bic                                                &
                  ! Ice particle mass-size relationships m(D)=ai D^bi
     &, lsp_ei, lsp_fi, lsp_eic, lsp_fic
                  ! Ice particle Best-Reynolds number relationships
                  ! Re(D) =LSP_EI(C) Be^LSP_FI(C)

      Integer                                                           &
                               !, INTENT(IN)
     &  cloud_fraction_method                                           &
                               ! Method for calculating
                               ! total cloud fraction
     &, ice_fraction_method    ! Method for calculating ice cloud frac.

      Real                                                              &
                               !, INTENT(IN)
     &  overlap_ice_liquid                                              &
                               ! Overlap between ice and liquid phases
     &, ctt_weight                                                      &
                               ! Weighting of cloud top temperature
     &, t_weight                                                        &
                               ! Weighting of local temperature
     &, qsat_fixed                                                      &
                               ! Fixed value of saturation humidity
     &, sub_cld                ! Scaling parameter
! Convection

      Integer                                                           &
     &  tr_vars                                                         &
                          ! number of free tracer variables
     &, tr_ukca                                             
                          ! number of ukca tracer variables

      Real                                                              &
     &  Charnock   ! Charnock parameter for sea surface
! RIVER routing

      INTEGER                                                           &
     & AOCPL_P_ROWS, AOCPL_ROW_LENGTH                                   &
     &, G_P_FIELD                                                       &
                                  ! IN size of global ATMOS field
     &, G_R_FIELD                                                       &
                                  ! IN Size of global river field
     &, RIVER_ROW_LENGTH                                                &
                                  ! IN local river row length
     &, RIVER_ROWS                                                      &
                                  ! IN local river rows
     &, GLOBAL_RIVER_ROW_LENGTH                                         &
                                  ! IN global river row length
     &, GLOBAL_RIVER_ROWS                                               &
                                        ! IN GLOBAL river rows
     &, A_STEPS_SINCE_RIV         ! IN No. Physics timsteps since last
!                                 ! call to river routing
      LOGICAL                                                           &
     & L_RIVERS               ! IN: TRUE if River routing needed

! Local parameters:
      INTEGER                                                           &
     &       swap_levels                                                &
                                  ! no. levels for SWAPBOUNDS
     &, gather_pe_trip                                                  &
                                  ! pe River routing to be run on
     &, info                                                            &
                              ! Return code from MPP
     &, icode                 ! Return code : 0 Normal Exit : >0 Error
      PARAMETER(swap_levels=1)              ! by definition for A- T

! data to regrid runoff from ATMOS to river routing grid
      REAL                                                              &
     & XPA(AOCPL_ROW_LENGTH+1)                                          &
                               ! IN Atmosphere TP long coordinates
     &,XUA(0:AOCPL_ROW_LENGTH)                                          &
                               ! IN Atmosphere U long coordinates
     &,XVA(AOCPL_ROW_LENGTH+1)                                          &
                               ! IN Atmosphere V long coordinates
     &,YPA(AOCPL_P_ROWS)                                                &
                               ! IN Atmosphere TP lat coordinates
     &,YUA(AOCPL_P_ROWS)                                                &
                               ! IN Atmosphere U lat coordinates
     &,YVA(0:AOCPL_P_ROWS)     ! IN Atmosphere V lat coordinates
! Data to run river routing
      REAL                                                              &
     & TRIVDIR(RIVER_ROW_LENGTH, RIVER_ROWS)                            &
                                              ! IN River direction file
     &,TRIVSEQ(RIVER_ROW_LENGTH, RIVER_ROWS)                            &
                                              ! IN River sequence file
     &,TWATSTOR(RIVER_ROW_LENGTH, RIVER_ROWS)                           &
                                              ! IN/OUT Water storage
!                                             ! file (Kg)
     &,RIVER_VEL                                                        &
                                              ! IN river velocity (m/s)
     &,RIVER_MCOEF                                                      &
                                              ! IN meander coefficient
     &,RIVER_STEP                    ! IN River routing step (secs)


#include "c_mdi.h"
!
       REAL                                                             &
     & r_area(ROW_LENGTH,ROWS),                                         &
!               accumulated areas file
     & r_inext(ROW_LENGTH,ROWS),                                        &
!               x-coordinate of downstream grid point
     & r_jnext(ROW_LENGTH,ROWS),                                        &
!               y-coordinate of downstream grid point
     &  slope(ROW_LENGTH,ROWS),                                         &
!             slopes (not used yet)
     &  flowobs1(ROW_LENGTH,ROWS),                                      &
!             initialisation for flows
     & r_land(ROW_LENGTH,ROWS),                                         &
!             land/river/sea
     & substore(ROW_LENGTH,ROWS),                                       &
!          routing sub_surface store (mm)
     & surfstore(ROW_LENGTH,ROWS),                                      &
!          routing surface store (mm)
     & flowin(ROW_LENGTH,ROWS),                                         &
!          surface lateral inflow (mm)
     & bflowin(ROW_LENGTH,ROWS)
!          sub-surface lateral inflow (mm)
! Diagnostics info
      REAL                                                              &
     & STASHWORK3(*)                                                    &
                         ! STASH workspace for section 3 (b layer)
     &,STASHWORK5(*)                                                    &
                         ! STASH workspace for section 5 (convection)
     &,STASHWORK8(*)                                                    &
                         ! STASH workspace for section 8 (hydrology)
     &,STASHwork9(*)                                                    &
                         ! STASH workspace for section 9 (LS Cloud)
     &,STASHwork19(*)                                                   &
                         ! STASH workspace for section 19 (Veg)
     &,STASHwork26(*)    ! STASH workspace for sect. 26 (River routing)
!

! Data arrays
      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &     0:model_levels)                                              &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          0:model_levels)                                         &
     &, rho_rsq(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &        model_levels)                                             &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, p_theta_levels(1-off_x:row_length+off_x,                        &
     &                   1-off_y:rows+off_y, model_levels)              &
     &, p_star(row_length, rows)                                        &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                   1-off_y:rows+off_y, model_levels)              &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      wet_levels)                                                 &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, visc_BL_m(1:row_length, 1:rows, bl_levels)                      &
!            ! visc_m only on BL levels
     &, FM_3D(row_length,rows,BL_LEVELS)                                &
!            ! stability function for momentum transport.
!            ! level 1 value is dummy for use in diagnostics
     &, FH_3D(row_length,rows,BL_LEVELS)                                &
!             ! stability function for heat and moisture.
!             ! level 1 value is dummy for use in diagnostics
     &, max_diff                                                        &
                  ! max diffusion coeff for run
     &, BL_COEF_KM(1:row_length, 1:rows, bl_levels-1)                   &
!            ! RHOKM from BL scheme
     &, BL_COEF_KH(1:row_length, 1:rows, bl_levels-1)
!            ! RHOKH from BL scheme

      Real :: mix_v_phys2 (1-off_x:row_length+off_x,                    &
     &                   1-off_y:rows+off_y, wet_levels)
             ! q_star converted to a mixing ratio
      Real :: mix_cl_phys2(1-off_x:row_length+off_x,                    &
     &                   1-off_y:rows+off_y, wet_levels)
             ! qcl_star converted to a mixing ratio
      Real :: mix_cf_phys2(1-off_x:row_length+off_x,                    &
     &                   1-off_y:rows+off_y, wet_levels)
             ! qcf_star converted to a mixing ratio
      Real :: unscaled_dry_rho (1-off_x:row_length+off_x,               &
     &                   1-off_y:rows+off_y, model_levels)
             ! unscaled dry density

      Logical                                                           &
     &  L_subfilter_vert                                                &
                            ! subgrid turbulence scheme in vertical     
     &, L_subfilter_horiz                                               &
                            ! subgrid turbulence scheme in horizontal   
     &, L_subfilter_blend                                               
                            ! Blending of diffusion coefficients in     
                            ! vertical

      Integer                                                           &
     & turb_startlev_vert                                               &
                            ! start and end vertical levels
     &,turb_endlev_vert     
                            ! for 3D turbulence scheme

! ancillary arrays and fields required to be saved from timestep to
! timestep.

      Real                                                              &
     &  T_surf(row_length, rows)

      logical                                                           &
     &  land_sea_mask(row_length, rows)

      Integer                                                           &
     &  land_index (land_points)                                        &
                                      ! set from land_sea_mask
     &, land_ice_index (land_points)                                    &
                                      ! Array of land ice points.
     &, soil_index(land_points)       ! Array of soil points.

      Real                                                              &
     &  u_0(row_length, rows)                                           &
                                ! set to zero
     &, v_0(row_length, n_rows)                                         &
                                ! set to zero
     &, u_0_p(row_length, rows)                                         &
                                  ! set to zero
     &, v_0_p(row_length, rows)                                         &
                                ! set to zero
     &, hcon (land_points)                                              &
                             ! soil/qrparm.soil.hcond
     &, hcap (land_points)                                              &
                             ! soil/qrparm.soil.hcap
     &, smvccl (land_points)                                            &
                             ! soil/qrparm.soil.crit
     &, smvcwt (land_points)                                            &
                             ! soil/qrparm.soil.wilt
     &, smvcst (land_points)                                            &
                             ! soil/qrparm.soil.satn
     &, sthf(land_points,dsm_levels)                                    &
                                      ! IN Frozen soil moisture content
!                of each layer as a fraction of saturation.
     &, sthu(land_points,dsm_levels)                                    &
                                      ! IN Unfrozen soil moisture
!                content of each layer as a fraction of saturation.
     &, canopy_gb (land_points)                                         &
                                ! set to zero.
     &, snow_depth (row_length, rows) ! snow/qrclim.snow.(month)

!SURFACE CURRENT INTERPOLATION VARIABLES 
      Real UHALO(1-off_x:row_length+off_x,1-off_y:rows+off_y) 
      Real VHALO(1-off_x:row_length+off_x,1-off_y:n_rows+off_y) 

      Real                                                              &
     &  ice_fract (row_length, rows)                                    &
                                     ! ice/qrclim.ice.(month)
     &, di(row_length, rows)                                            &
                              ! ice/qrclim.ice_thick.(month)
     &, ice_fract_ncat(row_length, rows, nice)                          &
     &, di_ncat(row_length, rows, nice)                                 &
     &, z0msea(row_length, rows)                                        &
                                   ! Sea surface roughness
     &, z0m_scm(row_length, rows)                                       &
                                   ! Fixed sea surface roughness
                                   ! length(m) for MOMENTUM,
                                   ! used in SCM
     &, z0h_scm(row_length, rows)  ! Fixed sea surface roughness
                                   ! length(m) for HEAT,
                                   ! used in SCM

      Real                                                              &
     &  sil_orog_land (land_points)                                     &
                                    ! orog/qrparm.orog.as
     &, ho2r2_orog (land_points)                                        &
                                 ! orog/qrparm.orog.h2root2
     &, t_soil(land_points,dsm_levels)                                  &
                                       ! slt/qrclim.slt_pm(lev).(month)
     &, ti_gb(row_length, rows)                                         &
                                   ! ice temp as grib box mean
     &, ti(row_length, rows, nice)                                      &
                                   ! ice temp on catagories
     &, t1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, q1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, zh(row_length, rows)                                            &
                             ! boundary layer height
     &, zh_prev(row_length, rows) ! boundary layer height from previous
!                                 ! timestep

      Real                                                              &
     &  clapp(land_points)                                              &
                             !  qrparm.soil.bwag ?
!                               Clapp-Hornberger exponent.
     &, satcon(land_points)                                             &
                                 !  qrparm.soil.satcon
     &, sathh(land_points)                                              &
                             !  soil water suction
     &, soil_layer_moisture(land_points, dsm_levels)
                             !  qrclim.smc_pm(lev).(month)


! CLoud fields

! local variables.
      Integer                                                           &
     &  rhc_row_length                                                  &
                        ! Row length for RHcrit array
     &, rhc_rows        ! Row number for RHcrit array
      Real                                                              &
     &  area_cloud_fraction(row_length, rows, wet_levels)               &
     &, bulk_cloud_fraction_halos(1-halo_i:row_length+halo_i,           &
     &                      1-halo_j:rows+halo_j, wet_levels)           &
     &, cloud_fraction_liquid_halos(1-halo_i:row_length+halo_i,         &
     &                      1-halo_j:rows+halo_j, wet_levels)           &
     &, cloud_fraction_frozen_halos(1-halo_i:row_length+halo_i,         &
     &                      1-halo_j:rows+halo_j, wet_levels)           &
     &, rhcpt(rhc_row_length, rhc_rows, wet_levels)

! Rain fields
      Real                                                              &
     &  ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, micro_tends(row_length, rows, bl_levels, 2)                     &
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)
     &, conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)

! Radiation fields
      Real                                                              &
     &  photosynth_act_rad(row_length, rows)                            &
                                             ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
     &, rad_hr(row_length, rows, bl_levels, 2)                          &
                                               !
!                               ! BL radiative (LW,SW) heating rates
     &, surf_radflux(row_length, rows)

! Fields for mineral dust source flux calculations
      REAL, INTENT(IN) ::                                               &
     &  SOIL_CLAY ( ROW_LENGTH, ROWS )                                  &
     &, SOIL_SILT ( ROW_LENGTH, ROWS )                                  &
     &, SOIL_SAND ( ROW_LENGTH, ROWS )                                  &
     &, DUST_MREL1 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL2 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL3 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL4 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL5 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL6 ( ROW_LENGTH, ROWS )

! Tracer emissions
      Real, Intent(In) ::                                               &
     &  so2_hilem ( row_length, rows )                                  &
     &, so2_em    ( row_length, rows )                                  &
     &, nh3_em    ( row_length, rows )                                  &
     &, dms_em    ( row_length, rows )                                  &
     &, soot_hilem( row_length, rows )                                  &
     &, soot_em   ( row_length, rows )                                  &
     &, ocff_hilem( row_length, rows )                                  &
     &, ocff_em   ( row_length, rows )

! CO2 fields
      Real                                                              &
     &  co2_emits ( row_length, rows )                                  &
     &, co2flux ( row_length, rows )

! Tracer fluxes - kdcorbin, 05/10
      Real                                                              &
     &  tracer_flux1 ( row_length, rows ), tracer_flux2( row_length, rows)  &
     &, tracer_flux3 ( row_length, rows ), tracer_flux4( row_length, rows)  &
     &, tracer_flux5 ( row_length, rows ), tracer_flux6( row_length, rows)  &
     &, tracer_flux7 ( row_length, rows ), tracer_flux8( row_length, rows)  &
     &, tracer_flux9 ( row_length, rows ), tracer_flux10(row_length, rows)  &
     &, tracer_flux11( row_length, rows ), tracer_flux12(row_length, rows)  &
     &, tracer_flux13( row_length, rows ), tracer_flux14(row_length, rows)  &
     &, tracer_flux15( row_length, rows ), tracer_flux16(row_length, rows)  &
     &, tracer_flux17( row_length, rows ), tracer_flux18(row_length, rows)  &
     &, tracer_flux19( row_length, rows ), tracer_flux20(row_length, rows)  

! CO2 Global Emissions - kdcorbin, 05/10
    REAL co2emitmass

! Additional Heat Sources for anthropogenic urban heat source. 
      REAL ANTHROP_HEAT(NTILES)                           ! (w/m2)
      LOGICAL, INTENT(INOUT) :: L_ANTHROP_HEAT_SRC 

! Convection
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
           !VarRes horizontal co-ordinate information
     &  dlambda_p(1-halo_i:row_length+halo_i)                           &
     &, dphi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)        &
     &, wt_lambda_p(1-halo_i:row_length+halo_i)                         &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)      &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

      Real                                                              &
     &  r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)

! Trig arrays
      real                                                              &
     &  cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)                          &
     &, sin_theta_latitude(row_length, rows)                            &
     &, FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)                      &
     &, sec_theta_latitude(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y)
! Information for STPH_SCV
      REAL,INTENT(In) ::                                                &
     & f3_at_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

! time information for current timestep
      Integer                                                           &
     &  val_year                                                        &
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second                                                      &
     &, timestep_number

! Diagnostic variables
      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP

! arguments with intent in/out. ie: input variables changed on output.

      Real                                                              &
     &  theta_star(1-off_x:row_length+off_x,                            &
     &               1-off_y:rows+off_y, model_levels)                  &
     &, q_star(1-off_x:row_length+off_x,                                &
     &           1-off_y:rows+off_y, wet_levels)                        &
     &, qcl_star(1-off_x:row_length+off_x,                              &
     &             1-off_y:rows+off_y, wet_levels)                      &
     &, qcf_star(1-off_x:row_length+off_x,                              &
     &             1-off_y:rows+off_y, wet_levels)                      &
     &, cf_star(1-off_x:row_length+off_x,                               &
     &           1-off_y:rows+off_y, wet_levels)                        &
     &, cfl_star(1-off_x:row_length+off_x,                              &
     &             1-off_y:rows+off_y, wet_levels)                      &
     &, cff_star(1-off_x:row_length+off_x,                              &
     &             1-off_y:rows+off_y, wet_levels)                      &
     &, R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels)                             &
     &, sum_eng_fluxes(row_length, rows)                                &
     &, sum_moist_flux(row_length, rows)

! Used in SCM for prescribed surface flux forcing
      Real                                                              &
     &  flux_e(row_length, rows)                                        &
                                 ! Surface latent heat flux (W/m^2)
     &, flux_h(row_length, rows) ! Surface sensible heat flux (W/m^2)

!    IN logicals for surface forcing
      LOGICAL                                                           &
     & L_flux_bc                                                        &
                    ! T if prescribed surface fluxes to be used
     &,L_spec_z0    ! T if roughness lengths have been specified

! Tracer variables
      Real, Intent(InOut) ::                                            &
     &  aerosol     ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,free_tracers( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y,                         &
     &                tr_levels, tr_vars)                               &
     & ,ukca_tracers( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y,                         &
     &                tr_levels, tr_ukca)
      REAL, INTENT(INOUT) ::                                            &
     &  DUST_DIV1   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV2   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV3   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV4   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV5   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV6   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )


      Real, Intent(InOut) ::                                            &
     &  so2         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,dms         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_aitken  ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_accu    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_diss    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,nh3         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )

      Real, Intent(InOut) ::                                            &
     &  soot_new    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,soot_aged   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,soot_cld    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_new   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_aged  ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_cld   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_new    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_aged   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_cld    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,co2         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )

! Add definitions for the cariolle scheme.

      Real, Intent(InOut) ::                                            &
     &   OZONE_TRACER( 1 - off_x : row_length + off_x,                  &
     &                 1 - off_y : rows + off_y, model_levels )         

! Variables with intent InOut for STPH_RP
      REAL,INTENT(InOut) :: par_mezcla ! Used to modify LAMBDAH,LAMBDAM
                                       !neutral mixing length in EXCOEF
      REAL,INTENT(InOut) :: G0_RP ! Stability function parameter in EXCOEF
! arguments with intent out. ie: output variables.

      Real                                                              &
     &  rhokm(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &           0:bl_levels-1)                                         &
     &, cH_term(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels-1)                                         &
     &, rhokm_u(row_length, rows, bl_levels)                            &
     &, rhokm_v(row_length, n_rows, bl_levels)

      Integer                                                           &
     &  ntml (row_length, rows)                                         &
     &, nbdsc(row_length, rows)                                         &
     &, ntdsc(row_length, rows)

      Logical                                                           &
     &  cumulus (row_length, rows) ! bl convection flag

      Integer                                                           &
     &  Error_code


! Additional variables for MOSES II
#include "nstypes.h"
#include "c_0_dg_c.h"

      LOGICAL                                                           &
     & L_NEG_TSTAR                                                      &
                                   ! IN Switch for -ve TSTAR error check
     &,L_SNOW_ALBEDO                                                    &
                                   ! IN Flag for prognostic snow albedo
     &,L_PHENOL                                                         &
                                   ! IN Indicates whether phenology
!                                  !    in use
     &,L_TRIFFID                                                        &
                                   ! IN Indicates whether TRIFFID
!                                  !    in use.
     &,L_TRIF_EQ                                                        &
                                   ! IN .T. for vegetation equilibrium.
     &,L_Q10                       ! IN .T. to use Q10 for soil resp

      Integer                                                           &
     & CO2_DIM_LEN                                                      &
                                   ! IN Length of a CO2 field row.
     &,CO2_DIM_ROW                                                      &
                                   ! IN Number of CO2 field rows.
     &,LAND_PTS_TRIF                                                    &
                                   ! IN For dimensioning land fields
     &,NPFT_TRIF                                                        &
                                   ! IN For dimensioning PFT fields
!                                  !    available only with TRIFFID.
!                                  !    Set to NPFT when TRIFFID on,
!                                  !    set to 1 when TRIFFID off.
     &,ASTEPS_SINCE_TRIFFID                                             &
                                   ! IN Number of atmospheric
!                                  !    timesteps since last call
!                                  !    to TRIFFID.
     &,A_STEP                                                           &
                                   ! IN Atmospheric timestep number.
     &,PHENOL_PERIOD                                                    &
                                   ! IN Phenology period (days).
     &,TRIFFID_PERIOD                                                   &
                                   ! IN TRIFFID period (days).
     &,CAN_MODEL                   ! IN Swith for thermal vegetation
!                                  !    canopy

      Real                                                              &
     & FRAC(LAND_POINTS,NTYPE)                                          &
                                      ! IN Fractions of surface types.
     &,FRAC_DISTURB(LAND_POINTS)                                        &
                                    ! IN Fraction of gridbox in which
!                                   !    vegetation is disturbed.
     &,CANHT_FT(LAND_POINTS,NPFT)                                       &
                                      ! IN Canopy height (m)
     &,LAI_FT(LAND_POINTS,NPFT)                                         &
                                      ! IN Leaf area index
     &,CANOPY(LAND_POINTS,NTILES)                                       &
                                      ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
     &,CATCH(LAND_POINTS,NTILES)                                        &
                                      ! IN Surface/canopy water capacity
!                                  !    of snow-free land tiles (kg/m2).
     &,CATCH_SNOW(LAND_POINTS)                                          &
                                   ! IN Coniferous canopy snow capacity
!                                  !    (kg/m2).
     &,SNOW_TILE(LAND_POINTS,NTILES)                                    &
                                      ! IN Lying snow on tiles (kg/m2)
     &,Z0_TILE(LAND_POINTS,NTILES)                                      &
                                      ! IN Tile roughness lengths (m).
     &,TSTAR_TILE(LAND_POINTS,NTILES)                                   &
                                      ! IN Surface tile temperatures
     &,INFIL_TILE(LAND_POINTS,NTILES)                                   &
!                          ! IN Maximum surface infiltration
     &,RGRAIN(LAND_POINTS,NTILES)                                       &
                                 ! INOUT Snow grain size (microns).
     &,SNOW_GRND(LAND_POINTS)                                           &
                                   ! INOUT Snow below canopy (kg/m2).
     &,CS(LAND_POINTS,DIM_CS1)                                          &
                                        ! IN Soil carbon (kg C/m2).
     &,GS(LAND_POINTS)                                                  &
                                      ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
     &,G_LEAF_ACC(LAND_POINTS,NPFT)                                     &
                                      ! INOUT Accumulated G_LEAF
     &,G_LEAF_PHEN_ACC(LAND_POINTS,NPFT)                                &
!                                  ! INOUT Accumulated leaf turnover
!                                  !       rate including phenology.
! Lestevens 17may13 - change npft to ntiles
     &,NPP_FT_ACC(LAND_PTS_TRIF,NTILES)                                 &
!    &,NPP_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                              &
!                                  ! INOUT Accumulated NPP_FT
     &,RESP_W_FT_ACC(LAND_PTS_TRIF,NTILES)                              &
!    &,RESP_W_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                           &
!                                  ! INOUT Accum RESP_W_FT
     &,RESP_S_ACC(LAND_PTS_TRIF,DIM_CS1)                                &
                                        ! INOUT Accumulated RESP_S
     &,OLR(ROW_LENGTH,ROWS)                                             &
                                ! IN    TOA - surface upward LW on
!                               !       last radiation timestep
!                               ! OUT   Corrected TOA outward LW
     &,LW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface downward LW radiation
!                                  !    (W/m2).
     &,SW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! Surface downward SW radiation (W/m2).
     &,SW_TILE(LAND_POINTS,NTILES)    ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
!-----------------------------------------------------------------------
!     EAK
      Real                                                          &
     &  alb_tile(land_points,ntiles,4)                              &
     &, SNOW_DEPTH3L(land_points,NTILES,3)                          &
     &, SNOW_MASS3L(land_points,NTILES,3)                           &
     &, SNOW_COND(land_points,NTILES,3)                             &
     &, SNOW_TMP3L(land_points,NTILES,3)                            &
     &, SNOW_RHO3L(land_points,NTILES,3)                            &
     &, SNOW_RHO1L(land_points,NTILES)                              &
     &, SNAGE_TILE(land_points,NTILES)                              &
     &, SMCL_TILE(land_points,NTILES,dsm_levels)                    &
     &, STHU_TILE(land_points,NTILES,dsm_levels)                    &
     &, STHF_TILE(land_points,NTILES,dsm_levels)                    &
     &, TSOIL_TILE(land_points,NTILES,dsm_levels)                   &
     &, T_SURF_TILE(land_points,NTILES)                             &
     &, RTSOIL_TILE(land_points,NTILES)                             &
     &, GFLUX_TILE(land_points,NTILES)                              &
     &, SGFLUX_TILE(land_points,NTILES)                             &
     &, HCONS(land_points)                                          &
     &, TRANSP_TILE(land_points,NTILES)                             &
     &, CPOOL_TILE(land_points,NTILES,10)                           &  
     &, NPOOL_TILE(land_points,NTILES,10)                           &
     &, PPOOL_TILE(land_points,NTILES,12)                           &
     &, SOIL_ORDER(land_points)                                     &
     &, NIDEP(land_points)                                          &
     &, NIFIX(land_points)                                          &
     &, PWEA(land_points)                                           &
     &, PDUST(land_points)                                          &
     &, GLAI(land_points,NTILES)                                    &
     &, PHENPHASE(land_points,NTILES)                               &

     &, surf_down_sw(row_length,rows,4)                                 &
     &, lat(row_length, rows)                                           &
     &, long(row_length, rows)                                          &
     &, time_sec
                                ! actual time of day in secs.
!       diagnostic variables for CABLE output
      Real                                                              &
     & FTL_TILE_CAB(LAND_POINTS,NTILES)                                 &
     &,FTL_CAB(LAND_POINTS)                                             &
     &,LE_TILE_CAB(LAND_POINTS,NTILES)                                  &
     &,LE_CAB(LAND_POINTS)                                              &
     &,TSTAR_TILE_CAB(LAND_POINTS,NTILES)                               &
     &,TSTAR_CAB(LAND_POINTS)                                           &
     &,SMCL_CAB(LAND_POINTS,DSM_LEVELS)                                 &
     &,TSOIL_CAB(LAND_POINTS,DSM_LEVELS)                                &
     &,USTAR_CAB(LAND_POINTS)                                           &
     &,SURF_HTF_CAB(LAND_POINTS)                                        &
!sxy
     &,TOT_ALB(LAND_POINTS,NTILES)                                      &
     &,CH_CAB(LAND_POINTS)                                              &
     &,CD_CAB(LAND_POINTS)                                              &
     &,U_S_CAB(LAND_POINTS)                                             &

     &,T1P5M(row_length, rows)                                          &
     &,Q1P5M(row_length, rows)

      Integer                                                           &
     &  day                                                             &
     &, total_nsteps                                                    &
                                ! Total number of steps in run
     &, SOIL_TYPE(row_length,rows)                                      &
     &, VEG_TYPE(row_length,rows)                                       &
     &, ISNOW_FLG3L(land_points,NTILES)

      Logical                                                           &
     & l_cable         !
      LOGICAL, DIMENSION(LAND_POINTS,NTILES) :: L_TILE_PTS

!-----------------------------------------------------------------------

! The following are only used if coastal tiling is switched on:
      REAL                                                              &
     & FLAND_CTILE(LAND_POINTS)                                         &
                                   ! IN Land fraction on land tiles.
     &,TSTAR_LAND_CTILE(ROW_LENGTH,ROWS)                                &
!                                  ! INOUT Land mean sfc temperature (K)
     &,TSTAR_SEA_CTILE(ROW_LENGTH,ROWS)                                 &
!                                  ! IN    Open sea sfc temperature (K).
     &,TSTAR_SICE_CTILE(ROW_LENGTH,ROWS)                                &
!                                  ! INOUT Sea-ice sfc temperature (K).
     &,ALBSOIL(LAND_POINTS)                                             &
                                   ! Soil albedo.
     &,COS_ZENITH_ANGLE(ROW_LENGTH,ROWS)
!                                  ! Cosine of the zenith angle
      INTEGER                                                           &
     & CAN_RAD_MOD                                                      &
                                   !Switch for canopy radiation model 
     &,ILAYERS                     !No of layers in canopy radiation model



! Additional variables required for large-scale hydrology:
      LOGICAL                                                           &
     & L_TOP                                                            &
                                  ! IN True if running TOPMODEL-based
                                  !    hydrology.
     &,L_PDM                                                            &
                                  ! IN True if running PDM.
     &,L_SOIL_SAT_DOWN            ! IN Direction of super-saturated 
                                  !    soil moisture

      REAL                                                              &
     & FEXP(LAND_POINTS)                                                &
                            ! IN Decay factor in Sat. Conductivity
!                           !    in water table layer.
     &,GAMTOT(LAND_POINTS)                                              &
                            ! IN Integrated complete Gamma function.
     &,TI_MEAN(LAND_POINTS)                                             &
                            ! IN Mean topographic index.
     &,TI_SIG(LAND_POINTS)  ! IN Standard dev. of topographic index.
!                           !    in water table layer.
      REAL                                                              &
     & FSAT(LAND_POINTS)                                                &
                            ! INOUT Surface saturation fraction.
     &,FWETL(LAND_POINTS)                                               &
                            ! INOUT Wetland fraction.
     &,ZW(LAND_POINTS)                                                  &
                            ! INOUT Water table depth (m).
     &,STHZW(LAND_POINTS)                                               &
                            ! INOUT soil moist fract. in deep-zw layer.
     &,A_FSAT(LAND_POINTS)                                              &
                            ! IN Fitting parameter for Fsat in LSH model
     &,C_FSAT(LAND_POINTS)                                              &
                            ! IN Fitting parameter for Fsat in LSH model
     &,A_FWET(LAND_POINTS)                                              &
                            ! IN Fitting parameter for Fwet in LSH model
     &,C_FWET(LAND_POINTS)  ! IN Fitting parameter for Fwet in LSH model

      REAL                                                              &
     & DUN_ROFF(LAND_POINTS)                                            &
                            ! OUT Dunne part of sfc runoff (kg/m2/s).
     &,QBASE(LAND_POINTS)                                               &
                            ! OUT Base flow (kg/m2/s).
     &,QBASE_ZW(LAND_POINTS)                                            &
                            ! OUT Base flow from ZW layer (kg/m2/s).
     &,DRAIN(LAND_POINTS)                                               &
                         ! Drainage out of nshyd'th level (kg/m2/s).
     &,FCH4_WETL(LAND_POINTS)
!                           ! OUT Wetland methane flux. (kg C/m2/s).
#if defined(SCMA)
!  Additional variables for SCM diagnostics

      Integer                                                           &
     &  nfor                                                            &
     &, ntrop                                                           &
     &, daycount                                                        &
     &, nconvars

      Character*8                                                       &
     &  time_string             ! String containing actual time

      Real                                                              &
     &  layer_cloud(row_length, rows,wet_levels)                        &
                                ! layer cloud amount
                                !  (decimal fraction)
     &  ,tstar(row_length, rows)! Surface temperature (K)

#endif

! Additional variables for SCM diagnostics which are dummy in full UM
      Integer                                                           &
     &  nSCMDpkgs               ! No of SCM diagnostics packages

      Logical                                                           &
     &  L_SCMDiags(nSCMDpkgs)   ! Logicals for SCM diagnostics packages

! loop counters
      Integer                                                           &
     &  i, j, k, l                                                      &
     &, j_begin, j_end                                                  &
     &, Substep_Number

      Logical                                                           &
     &  L_poles   !  include poles in etadot calc (false for LAMs)

! Diagnostic switches
! a) hydrology
      Logical                                                           &
     &  stf_sub_surf_roff                                               &
     &, smlt

! Sub timestep for ATMPHYS2
      Real                                                              &
     & sub_timestep        ! sub_timestep = timestep/Num_Substeps
                           ! i.e. the phys2 integration timestep.
! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='Atm_Physics2')

! local variables

! Local data arrays

#if defined(SCMA)
      Real                                                              &
     &  resdump(row_length, rows, nconvars)                             &
     &, t(row_length, rows, model_levels)
#endif
      Real                                                              &
     &  T_latest(row_length, rows, model_levels)                        &
     &, q_latest(row_length, rows, wet_levels)                          &
     &, qcl_latest(row_length, rows, wet_levels)                        &
     &, qcf_latest(row_length, rows, wet_levels)                        &
     &, cf_latest(row_length, rows, wet_levels)                         &
     &, cfl_latest(row_length, rows, wet_levels)                        &
     &, cff_latest(row_length, rows, wet_levels)                        &
     &, bulk_cloud_fraction(row_length, rows, wet_levels)               &
     &, cloud_fraction_liquid(row_length, rows, wet_levels)             &
     &, cloud_fraction_frozen(row_length, rows, wet_levels)             &
     &, theta_conv(row_length, rows, model_levels)                      &
     &, u_conv(1-off_x:row_length+off_x,1-off_y:rows+off_y,bl_levels)   &
     &, v_conv(1-off_x:row_length+off_x,1-off_y:n_rows+off_y,bl_levels) &
! The above two arrays are activated when L_phys2_substep=T and are
! equal to the horiz wind comp before the application of the expl BL.
! They are used to form the non turbulent terms in the BL solver.

     &, q_conv(row_length, rows, wet_levels)                            &
     &, qcl_conv(row_length, rows, wet_levels)                          &
     &, qcf_conv(row_length, rows, wet_levels)                          &
     &, cf_liquid_conv(row_length, rows, wet_levels)                    &
     &, cf_frozen_conv(row_length, rows, wet_levels)                    &
     &, bulk_cf_conv(row_length, rows, wet_levels)                      &
     &, theta_inc(row_length, rows, model_levels)                       &
     &, q_inc(row_length, rows, wet_levels)                             &
     &, qcl_inc(row_length, rows, wet_levels)                           &
     &, qcf_inc(row_length, rows, wet_levels)                           &
     &, cf_liquid_inc(row_length, rows, wet_levels)                     &
     &, cf_frozen_inc(row_length, rows, wet_levels)                     &
     &, bulk_cf_inc(row_length, rows, wet_levels)                       &
     &, l_s_poles(row_length, first_constant_r_rho_level-1)             &
     &, l_n_poles(row_length, first_constant_r_rho_level-1)
!
      Real                                                              &
     &  p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, p_layer_centres(row_length, rows, 0:model_levels)               &
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = pstar, and at top = 0.
     &, exner_layer_boundaries(row_length, rows, 0:model_levels)        &
     &, exner_layer_centres(row_length, rows, 0:model_levels)

! Local arrays used by both conv_diag and convection. Calculations moved
! from conv_diag & convection to save CPU.

      Real ::                                                           &
     &  z_theta(row_length, rows, model_levels)                         &
                        ! height of theta levels above surface (m)
     &, z_rho(row_length, rows, model_levels)                           &
                        ! height of rho levels above surface (m)
     &, rho_wet(row_length, rows, model_levels)                         &
                        ! wet density on rho levels (kg/m3)
     &, rho_wet_theta(row_length, rows, model_levels-1)                 &
                        ! wet density on theta levels (not top) (kg/m3)
     &, rho_dry(row_length, rows, model_levels)
                        ! dry density on rho levels (kg/m3)

      Real                                                              &
     &  etadot_copy(row_length,rows,0:model_levels)                     &
     &, Z_TOP_OF_MODEL     ! eta = (r-r_surf)/(r_top-r_surf)
                           !     = z/z_top

! Local arrays holding information to be passed between physics
! routines.

      Real                                                              &
     &  ls_rain_land(land_points)                                       &
     &, ls_snow_land(land_points)                                       &
     &, conv_rain_land(land_points)                                     &
     &, conv_snow_land(land_points)
! Local Arrays for convection
      Integer                                                           &
     &  lcbase (row_length, rows)                                       &
                                     ! lowest convective cloud base leve
     &, lctop (row_length, rows)     ! lowest convective cloud top level

!                                        ! shallow convection
      Real                                                              &
     &  lcca(row_length, rows)                                          &
                                ! lowest convective cloud amount (%)
     &, lcclwp(row_length, rows)                                        &
                                 ! lowest cloud condensed water
                                  ! path (KG/M**2)
     &, ccw(row_length, rows, wet_levels)
                                  ! convective cloud liquid water
                                  ! (G/KG) on model levels

      !=================================================================
      ! DO NOT REMOVE OR WRITE TO if lcv_ccrad=.FALSE.
      !
      ! Intermediary versions of ccw lcbase. These are used when
      ! lcv_ccrad=.TRUE.. They exist because no space is allocated in the
      ! D1 array in atm_step when lcv_ccrad=.FALSE., without them
      ! overwriting will occur. They should not be removed unless an
      ! alternative in atm_step can be found.

      Real    :: ccw_out    (row_length, rows, wet_levels)
      Integer :: lcbase_out (row_length, rows)            

      ! DO NOT REMOVE OR WRITE TO if lcv_ccrad=.FALSE.
      !=================================================================

      Integer                                                           &
     &  ntpar(row_length, rows)                                         &
                                     ! top of diagnostic parcel ascent
     &, nlcl(row_length, rows)       ! lifting condensation level

      Logical ::                      &
      l_shallow(row_length, rows)     & ! Logical indicator of shallow 
                                        ! shallow 
     , l_congestus(row_length, rows)  & ! Logical indicator of congestus 
     , l_congestus2(row_length, rows) & ! Logical indicator of congestus
     , l_pc2_diag_sh_pts(row_length, rows)  ! PC2 is using diagnostic
                                            ! shallow convection

      real                                                              &
     & CIN_undilute(row_length, rows)                                   &
                                       ! undilute CIN from parcel ascent
                                       ! (m2/s2)
     &,CAPE_undilute(row_length, rows) ! undilute CAPE from parcel
                                       ! ascent (m2/s2)

      Real ::                                                           &
     &  delthvu(row_length, rows)                                       &
                                     ! buoyancy integral
     &, zhpar(row_length, rows)                                         &
                                     ! height of ntpar
     &, zlcl(row_length, rows)                                          &
                                     ! height of lcl accurate value rather
                                     ! than a model level (m)
     &, zlcl_uv(row_length,rows)                                        &
                                     ! height of nlcl for uv grid
     &, ql_ad(row_length,rows)       ! adiabatic liquid water content
                                     ! at inversion (kg/kg)

      Real                                                              &
     &  uw0(row_length, rows)                                           &
!                       ! U-component of surface wind stress (P-grid)
     &, vw0(row_length, rows)
!                       ! V-component of surface wind stress (P-grid)
!     Declaration of new BL diagnostics.
      Type (Strnewbldiag) :: BL_diag

      Real                                                              &
     &  wstar(row_length, rows)                                         &
                                     ! surface-based mixed layer
!                                    ! velocity scale
     &, wthvs(row_length, rows)                                         &
                                     ! surface buoyancy flux
     &, w_max(row_length, rows) ! max w in column

! Diagnostics controlled by Diagnostic switches

! Convection diagnostics. Moved from CONV_CTL2 to enable substepping.

      Real                                                              &
     &  uw_dp(row_length, rows, model_levels)                           &
                                ! x-stress for deep convection
                                ! (kg m-1 s-2)
     &, vw_dp(row_length, rows, model_levels)                           &
                                ! y-stress for deep convection
                                ! (kg m-1 s-2)
     &, uw_shall(row_length, rows, model_levels)                        &
                                ! x-stress for shallow convection
                                ! (kg m-1 s-2)
     &, vw_shall(row_length, rows, model_levels)                        &
                                ! y-stress for shallow convection
                                ! (kg m-1 s-2)
                                ! from downdraughts
     &, up_flux(row_length, rows, model_levels)                         &
                                                !  updraught mass flux
     &, dwn_flux(row_length, rows, model_levels)                        &
                                                 ! downdraught mass flux
     &, entrain_up(row_length, rows, model_levels)                      &
                                ! fractional entrainment rate
     &, detrain_up(row_length, rows, model_levels)                      &
                                ! fractional detrainment rate
     &, entrain_dwn(row_length, rows,model_levels)                      &
                                ! fractional entrainment rate
                                ! into downdraughts
     &, detrain_dwn(row_length, rows,model_levels)                      &
                                ! fractional detrainment rate
                                ! from downdraughts
     &, precip_deep(row_length,rows)                                    &
                                      ! deep precipitation
     &, precip_shall(row_length,rows)                                   &
                                      ! shallow precipitation
     &, precip_mid(row_length,rows)                                     &
                                      ! mid precipitation
     &, precip_cong(row_length,rows)                                    &
                                      ! congestus precipitation
     &, cape(row_length,rows)                                           &
                                      ! CAPE
     &, shallow_ind(row_length,rows)                                    &
                                      ! shallow indicator
     &, congestus_ind(row_length,rows)                                  &
                                        ! congestus indicator
     &, congestus_ind2(row_length,rows)                                 &
                                         ! congestus indicator2
     &, mid_ind(row_length,rows)                                        &
                                      ! mid indicator
     &, ntml_diag(row_length,rows)                                      &
                                      ! NTML output diag
     &, ntpar_diag(row_length,rows)                                     &
                                      ! NTPAR output diag
     &, freeze_diag(row_length,rows)                                    &
                                       ! freezing level output diag
     &, kterm_diag(row_length,rows)                                     &
                                      ! deep termination level diag
     &, wstar_up_diag(row_length,rows)                                  &
                                       ! cumulus layer convective vel
     &, wstar_dn_diag(row_length,rows)                                  &
                                       ! subcloud layer convective vel
     &, mb1_diag(row_length,rows)                                       &
                                        ! cloud base mass flux 1
     &, mb2_diag(row_length,rows)                                       &
                                        ! cloud base mass flux 2
     &, wqt_cb(row_length, rows)                                        &
                                     ! w'qt cloud base flux
     &, wthetal_cb(row_length, rows)                                    &
                                     ! w'thetal' cb flux
     &, wqt_inv(row_length, rows)                                       &
                                      ! w'qt inversion base flux
     &, wthetal_inv(row_length, rows)                                   &
                                      ! w'thetal' inversion flux
     &, sh_top(row_length, rows)                                        &
                                  ! height of top of shallow conv (m)
     &, sh_base(row_length, rows)                                       &
                                  ! height of base of shallow conv (m)
     &, cg_top(row_length, rows)                                        &
                                  ! height of top of congestus conv (m)
     &, cg_base(row_length, rows)                                       &
                                  ! height of base of congestus conv (m)
     &, cg_term(row_length, rows)                                       &
                                  ! congestus conv termination level
     &,  CONSCAV_DUST(ROW_LENGTH, ROWS, 6)                              &
                                          !col total scvngd for
                                          !DIV1,DIV2,...,DIV6 dust
     &,  conscav_so4ait(row_length, rows)                               &
                                          !column total scvngd so4ait
     &,  conscav_so4acc(row_length, rows)                               &
                                          !column total scvngd so4acc
     &,  conscav_so4dis(row_length, rows)                               &
                                          !column total scvngd so4dis
     &,  conscav_agedsoot(row_length, rows)                             &
                                            !colm total scvngd agd soot
     &,  conscav_agedbmass(row_length, rows)                            &
                                            !colm total scvngd bmass
     &,  conscav_agedocff(row_length, rows)
                                            !colm total scvngd ocff

      Real                                                              &
     &  sub_surf_roff(land_points)                                      &
                                    ! sub-surface runoff
     &, snomlt_sub_htf(land_points)&! subsurface snowmelt heat flux
     &, wblake_ratio               &! ratio of wblake/subroff
     &, WBLAKE(land_points)        &! 
     &, SUBROFF(land_points)       &! 
     &, SUBROFF0(land_points)       &! 
     &, SUBROFF1(land_points)       &! 
     &, TOT_WBLAKE(land_points)    &! 
     &, TOT_SUBRUN(land_points)   &! 
     &, WB_LAKE(land_points,ntiles) ! 

       Real TOTWBLAKE,TOTWBLAKE0,TOTWBLAKE1,TOTRN,TOTWBX,TOTWBX0,TOTWBX1, &
            TOTAREA

      INTEGER   N

      Real                                                              &
     &  infil(land_points)                                              &
                            ! max surface infiltration rate (kg/m2/s)
     &, smc(land_points)                                                &
                            ! out available soil moisture in the
!                                 rootzone (kg/m2).
     &, snow_melt(land_points)                                          &
                                ! snowmelt (kg/m2/s).
     &, surf_roff(land_points)                                          &
                                ! surface runoff (kg/m2/s).
     &, tot_tfall(land_points)  ! total throughfall (kg/m2/s).

! Fields passed from BDY_LAYR to IMP_SOLVER
      Real                                                              &
     &  alpha1_sice(row_length,rows)                                    &
     &, ashtf(row_length,rows)                                          &
     &, dtrdz_charney_grid(row_length,rows,BL_LEVELS)                   &
     &, rdz_charney_grid(row_length,rows,BL_LEVELS)                     &
     &, dtrdz_u(row_length,rows,BL_LEVELS)                              &
     &, dtrdz_v(row_length,n_rows,BL_LEVELS)                            &
     &, rdz_u(row_length,rows,2:BL_LEVELS)                              &
     &, rdz_v(row_length,n_rows,2:BL_LEVELS)                            &
     &, cdr10m_u(row_length,rows)                                       &
     &, cdr10m_v(row_length,n_rows)                                     &
     &, z1_tq(row_length,rows)                                          &
     &, rhokh (row_length, rows, bl_levels)


! OUT variables which need to be maintained between substeps. They will
!     be used by BDY_LAYR(). Moved in atmos_physics2 as they don't need
!     to be recalculated at each substep. Calculate at 1st substep and
!     maintain their value in memory.
!
      REAL                                                              &
     & RHO_UV(row_length,rows,BL_LEVELS+1)                              &
!                               ! RHO_UV(*,K) is the density at half
!                               ! level k-1/2.
     &,RHO_TQ(row_length,rows,BL_LEVELS)                                &
!                               ! RHO_TQ(*,K) is the density at half
!                               ! level k+1/2.
     &,DZL_charney(row_length,rows,BL_LEVELS)                           &
                                ! DZL(,K) is depth in m of layer
!                                 K, i.e. distance from boundary
!                                 K-1/2 to boundary K+1/2.
     &,RDZ(row_length,rows,BL_LEVELS)                                   &
                                ! RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
     &,Z1_UV(row_length,rows)                                           &
                                ! Height of lowest u,v level.
     &,Z_FULL(row_length,rows,model_levels)                             &
                                ! Z_FULL(*,K) is height of full level k.
     &,Z_UV(row_length,rows,BL_LEVELS)                                  &
                                ! Z_UV(*,K) is height of half level
!                               ! k-1/2.
     &,Z_TQ(row_length,rows,BL_LEVELS)                                  &
                                ! Z_TQ(*,K) is height of half level
!                               ! k+1/2.
     &,P_HALF(row_length,rows,model_levels)                             &
                                ! P_HALF(*,K) is pressure at half
!                               ! level k-1/2.
     &,DELTAP(row_length,rows,BL_LEVELS)
                                ! Difference in pressure between levels

! Additional fields
      REAL                                                              &
     & FLAND(LAND_POINTS)                                               &
                                   ! Land fraction on land tiles.
     &,FLANDG(ROW_LENGTH,ROWS)                                          &
                                   ! Land fraction on global field.
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
                                   ! Land mean sfc temperature (K)
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! Open sea sfc temperature (K).
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! Sea-ice sfc temperature (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS)  ! Sea mean sfc temperature (K).
! diagnostics
      Real                                                              &
            ! output from bdy_layr.
     &  cd(row_length, rows)                                            &
     &, ch(row_length, rows)                                            &
     &, e_sea(row_length, rows)                                         &
     &, fqT(row_length, rows, bl_levels)                                &
     &, ftl(row_length, rows, bl_levels)                                &
     &, h_sea(row_length, rows)                                         &
     &, rib_gb(row_length,rows)                                         &
     &, taux(row_length, rows, bl_levels)                               &
     &, tauy(row_length, n_rows, bl_levels)                             &
     &, vshr(row_length, rows)                                          &
     &, zht(row_length, rows)                                           &
     &, shallowc(row_length,rows)                                       &
                                   !  Indicator set to 1.0 if shallow,
!                                  !   0.0 if not shallow or not cumulus
     &, cu_over_orog(row_length,rows)                                   &
!                                  !  Indicator for cumulus
!                                  !     over steep orography
!                                  !   Indicator set to 1.0 if true,
!                                  !   0.0 if false. Exclusive.
     &, bl_type_1(row_length,rows)                                      &
                                   !  Indicator set to 1.0 if stable
!                                  !     b.l. diagnosed, 0.0 otherwise.
     &, bl_type_2(row_length,rows)                                      &
                                   !  Indicator set to 1.0 if Sc over
!                                  !     stable surface layer diagnosed,
!                                  !     0.0 otherwise.
     &, bl_type_3(row_length,rows)                                      &
                                   !  Indicator set to 1.0 if well
!                                  !     mixed b.l. diagnosed,
!                                  !     0.0 otherwise.
     &, bl_type_4(row_length,rows)                                      &
                                   !  Indicator set to 1.0 if
!                                  !     decoupled Sc layer (not over
!                                  !     cumulus) diagnosed,
!                                  !     0.0 otherwise.
     &, bl_type_5(row_length,rows)                                      &
                                   !  Indicator set to 1.0 if
!                                  !     decoupled Sc layer over cumulus
!                                  !     diagnosed, 0.0 otherwise.
     &, bl_type_6(row_length,rows)                                      &
                                   !  Indicator set to 1.0 if a
!                                  !     cumulus capped b.l. diagnosed,
!                                  !     0.0 otherwise.
     &, bl_type_7(row_length,rows)                                      &
                                   !  Indicator set to 1.0 if a
!                                  !     shear-dominated b.l.
!                                  !      diagnosed, 0.0 otherwise.
     &, z0m_eff_gb(row_length,rows)                                     &
     &, z0h_eff_gb(row_length,rows)                                     &
     &, fme(row_length, rows)                                           &
     &, WT_EXT(LAND_POINTS,dsm_levels)                                  &
                                       ! cumulative fract of transp'n
     &, RA(LAND_POINTS) ! aerodynamic resiatance (s/m)

#if defined(SCMA)
      integer ntrad                                                     &
                       ! No. of timesteps between
                       ! calls to radiation
     &     ,ntrad1                                                      &
                       ! 1st timestep on which radiation called
     &     ,conv_mode  ! Determines actions for convection scheme
#endif

! data required for tracer mixing :
      Real                                                              &
     &  rho_aresist(row_length,rows)                                    &
     &, aresist(row_length,rows)                                        &
     &, resist_b(row_length,rows)                                       &
     &, R_B_DUST(ROW_LENGTH,ROWS,NDIV)                                  &
                                      ! surface layer resist for dust
     &, DUST_FLUX(ROW_LENGTH,ROWS,NDIV)                                 &
                                       ! dust emissions (kg m-2 s-1)
     &, DRYDEP2(ROW_LENGTH,ROWS,NDIV)                                   &
                                      !dry dep though grav. set.
     &, U_S_T_TILE(LAND_POINTS,NTILES,NDIV)                             &
                                           !OUT threshold frict. vel
     &, U_S_T_DRY_TILE(LAND_POINTS,NTILES,NDIV)                         &
                                               !OUT dry soil value
     &, U_S_STD_TILE(LAND_POINTS,NTILES)                                &
                                        !OUT friction velocity
     &, WE_LIM(row_length,rows,3)                                       &
                                    !  rho*entrainment rate implied by
!                                   !     placing of subsidence
     &, ZRZI(row_length,rows,3)                                         &
                                    !  (z-z_base)/(z_i-z_base)
     &, T_FRAC(row_length,rows,3)                                       &
                                    !  a fraction of the timestep
     &, WE_LIM_DSC(row_length,rows,3)                                   &
!                                   !  rho*entrainment rate implied by
!                                   !     placing of subsidence
     &, ZRZI_DSC(row_length,rows,3)                                     &
                                    !  (z-z_base)/(z_i-z_base)
     &, T_FRAC_DSC(row_length,rows,3)                                   &
!                                   !  a fraction of the timestep
     &, Z_HALF(row_length,rows,bl_levels)                               &

                                    !  Z_HALF(*,K) is height of half
!                                   ! level k-1/2.
     &, ZHSC(row_length,rows)       !  Top of decoupled layer

      Integer                                                           &
     &  KENT(row_length,rows)                                           &
                                    !  grid-level of SML inversion
     &, KENT_DSC(row_length,rows)   !  grid-level of DSC inversion

! data required for 4D-VAR :
      Real                                                              &
             ! (output from sf_exch)
     &  rho_cd_modv1(row_length, rows)                                  &
     &, rho_km_var(row_length,rows,2:bl_levels)

! Fields passed out of boundary layer into hydrology
      Real                                                              &
             !
     &  ecan(row_length, rows)                                          &
                                        !output from sf_evap.
     &, ei(row_length, rows)                                            &
                                        !output from sf_evap.
     &, snowmelt(row_length, rows)                                      &
                                        !output from sf_evap.
     &, ext(land_points,dsm_levels) ! Extraction of water from each
!                                    soil layer (kg/m2/s).

      Real                                                              &
             ! ( needed for soil_hft )
     &  surf_ht_flux_land(row_length, rows)                             &
                                             !
     &, surf_ht_flux_ld(land_points)                                    &
                                             !
     &, snomlt_surf_htf(row_length, rows)

      Real                                                              &
             ! (STASH diagnostic)
     &  LYING_SNOW(LAND_POINTS)  ! Gridbox snowmass (kg/m2).
      Real                                                              &
     &  cca_2d(row_length, rows)                                        &
     &, lclf

      Real                                                              &
     &  tot_precip_scaled(row_length, rows)

! logicals
      Logical                                                           &
     &  L_zero_boundaries                                               &
     &, L_zero_halos
!
      Logical                                                           &
     & L_scrn                                                           &
                                 ! Logical to control output
                                 !    of screen level T,Q,QCL,QCF
     &,L_plsp                    ! Logical to control output
                                 !    of Probability of LS Precip
!
      Logical                                                           &
     &  L_calc_dxek                                                     &
                     ! Switch for calculation of condensate increment
     &, L_q_interact ! Switch allows overwriting of parcel variables
!                      when calculating condensate increments.
!

! Additional variables for MOSES II
      INTEGER                                                           &
     & TILE_INDEX(LAND_POINTS,NTYPE)                                    &
                                      ! Index of tile points
     &,TILE_PTS(NTYPE)             ! Number of tile points

      Real                                                              &
     & FTL_TILE(LAND_POINTS,NTILES)                                     &
                                      ! Surface FTL for land tiles
     &,LE_TILE(LAND_POINTS,NTILES)                                      &
                                      ! Surface latent heat flux for
!                                  !     land tiles
     &,RADNET_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! Surface net radiation on
!                                  !     sea-ice (W/m2)
     &,RADNET_TILE(LAND_POINTS,NTILES)                                  &
                                      ! Surface net radiation on
!                                  !     land tiles (W/m2)
     &,RIB_TILE(LAND_POINTS,NTILES)                                     &
                                      ! RIB for land tiles.
     &,RHO_ARESIST_TILE(LAND_POINTS,NTILES)                             &
!                                  ! RHOSTAR*CD_STD*VSHR on land
!                                  !     tiles
     &,ARESIST_TILE(LAND_POINTS,NTILES)                                 &
!                                  ! 1/(CD_STD*VSHR) on land tiles
     &,RESIST_B_TILE(LAND_POINTS,NTILES)                                &
!                                  ! (1/CH-1/CD_STD)/VSHR on land
!                                  !     tiles
     &,ALPHA1(LAND_POINTS,NTILES)                                       &
                                   ! Mean gradient of saturated
!                                  !     specific humidity with respect
!                                  !     to temperature between the
!                                  !     bottom model layer and tile
!                                  !     surfaces
     &,ASHTF_TILE(LAND_POINTS,NTILES)                                   &
                                      !Coefficient to calculate
!                                  !     surface heat flux into land
!                                  !     tiles.
     &,FQT_TILE(LAND_POINTS,NTILES)                                     &
                                      ! Surface FQT for land tiles
     &,EPOT_TILE(LAND_POINTS,NTILES)                                    &
!                                  ! OUT Local EPOT for land tiles.
     &,FQT_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! Surface FQT for sea-ice
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! Surface FTL for sea-ice
     &,FRACA(LAND_POINTS,NTILES)                                        &
                                   ! Fraction of surface moisture
!                                  !     flux with only aerodynamic
!                                  !     resistance for snow-free land
!                                  !     tiles.
     &,RESFS(LAND_POINTS,NTILES)                                        &
                                   ! Combined soil, stomatal
!                                  !     and aerodynamic resistance
!                                  !     factor for fraction (1-FRACA)
!                                  !     of snow-free land tiles.
     &,RESFT(LAND_POINTS,NTILES)                                        &
                                   ! Total resistance factor.
!                                  !     FRACA+(1-FRACA)*RESFS for
!                                  !     snow-free land, 1 for snow.
     &,RHOKH_TILE(LAND_POINTS,NTILES)                                   &
                                      ! Surface exchange coefficients
!                                  !     for land tiles
     &,RHOKH_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! Surface exchange coefficients
!                                  !     for sea and sea-ice
     &,RHOKPM(LAND_POINTS,NTILES)                                       &
                                   ! Surface exchange coefficient.
     &,RHOKPM_POT(LAND_POINTS,NTILES)                                   &
!                                  ! Potential evaporation
!                                  !     exchange coefficient.
     &,RHOKPM_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! Sea-ice surface exchange coeff.
     &,Z0HSSI(ROW_LENGTH,ROWS)                                          &
     &,Z0MSSI(ROW_LENGTH,ROWS)                                          &
                                   ! Roughness lengths over sea (m).
     &,Z0H_TILE(LAND_POINTS,NTILES)                                     &
                                   ! Tile roughness lengths for heat
!                                  !     and moisture (m).
     &,Z0M_GB(ROW_LENGTH,ROWS)                                          &
                                   ! Gridbox mean roughness length for
!                                  !     momentum (m).
     &,Z0M_TILE(LAND_POINTS,NTILES)                                     &
                                   ! Tile roughness lengths for
!                                  !     momentum.
     &,CHR1P5M(LAND_POINTS,NTILES)                                      &
                                   ! Ratio of coefffs for
!                                  !     calculation of 1.5m temp for
!                                  !     land tiles.
     &,CHR1P5M_SICE(ROW_LENGTH,ROWS)                                    &
!                                  ! CHR1P5M for sea and sea-ice
!                                  !     (leads ignored).
     &,GPP(LAND_POINTS)                                                 &
                                   ! Gross primary productivity
!                                  !     (kg C/m2/s).
     &,NPP(LAND_POINTS)                                                 &
                                   ! Net primary productivity
!                                  !     (kg C/m2/s).
     &,RESP_P(LAND_POINTS)                                              &
                                   ! Plant respiration (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,G_LEAF(LAND_POINTS,NTILES)                                         &
                                   ! Leaf turnover rate (/360days).
     !kdcorbin, 11/10 - changed from NPFT
     &,GPP_FT(LAND_POINTS,NTILES)                                       &
                                   ! Gross primary productivity
                                   !     on PFTs (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,NPP_FT(LAND_POINTS,NTILES)                                       &
                                   ! Net primary productivity
                                   !     (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,RESP_P_FT(LAND_POINTS,NTILES)                                    &
                                   ! Plant respiration on PFTs
                                   !     (kg C/m2/s).
     &,RESP_S(LAND_POINTS,DIM_CS1)                                      &
                                        ! Soil respiration (kg C/m2/s).
     &,RESP_S_TOT(DIM_CS2)                                              &
                                        ! OUT total RESP_S over pools
     &,RESP_S_TILE(LAND_POINTS,NTILES)                                  &
                                  ! Soil respiration on tiles (kg C/m2/s).
                                  ! kdcorbin, 10/10
     &,RESP_W_FT(LAND_POINTS,NPFT)                                      &
                                   ! Wood maintenance respiration
!                                  !     (kg C/m2/s).
     &,GC(LAND_POINTS,NTILES)                                           &
                                   ! "Stomatal" conductance to
!                                  !      evaporation for land tiles
!                                  !      (m/s).
     &,CANHC_TILE(LAND_POINTS,NTILES)                                   &
                                      ! Areal heat capacity of canopy
!                                  !    for land tiles (J/K/m2).
     &,WT_EXT_TILE(LAND_POINTS,DSM_LEVELS,NTILES)                       &
!                                  ! Fraction of evapotranspiration
!                                  !    which is extracted from each
!                                  !    soil layer by each tile.
     &,FLAKE(LAND_POINTS,NTILES)                                        &
                                   ! Lake fraction.
     &,TILE_FRAC(LAND_POINTS,NTILES)                                    &
                                      ! Tile fractions including
!                                  !     snow cover in the ice tile.
     &,FSMC(LAND_POINTS,NPFT)                                           &
                                   ! Soil moisture availability factor.
     &,FLANDG_U(ROW_LENGTH,ROWS)                                        &
                                   ! Land frac (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,FLANDG_V(ROW_LENGTH,N_ROWS)                                      &
                                   ! Land frac (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,rib_ssi(row_length,rows)                                         &
                                   ! Rib over sea part of grid box.
     &,vshr_land(row_length, rows)                                      &
                                   ! Vshr over land part of grid box.
     &,vshr_ssi(row_length, rows)                                       &
                                   ! Vshr over sea part of grid box.
     &,TAUX_LAND(ROW_LENGTH,ROWS)                                       &
                                   ! Taux over land part of grid box.
     &,TAUX_SSI(ROW_LENGTH,ROWS)                                        &
                                   ! Taux over sea part of grid box.
     &,TAUY_LAND(ROW_LENGTH,N_ROWS)                                     &
                                   ! Tauy over land part of grid box.
     &,TAUY_SSI(ROW_LENGTH,N_ROWS) ! Tauy over sea part of grid box.

      REAL                                                              &
     & ESOIL_TILE(LAND_POINTS,NTILES)                                   &
                                ! Evaporation from bare soil (kg/m2)
     &,ES(row_length,rows)                                              &
                                ! Surface evapotranspiration from
!                               !     soil moisture store (kg/m2/s).
     &,EI_TILE(LAND_POINTS,NTILES)                                      &
                                   ! EI for land tiles
     &,Q1P5M_TILE(LAND_POINTS,NTILES)                                   &
!                               ! Q1P5M over land tiles.
     &,T1P5M_TILE(LAND_POINTS,NTILES)                                   &
!                               ! T1P5M over land tiles.
     &,ECAN_TILE(LAND_POINTS,NTILES)                                    &
                                ! ECAN for land tiles
     &,MELT_TILE(LAND_POINTS,NTILES)
!                               ! Snowmelt on tiles (kg/m2/s).


      INTEGER                                                           &
     & PHENOL_CALL                                                      &
                    ! indicates whether phenology is to be called
     &,TRIFFID_CALL                                                     &
                    ! indicates whether TRIFFID is to be called
     &,NSTEP_TRIF   ! Number of atmospheric timesteps between calls to
!                   ! TRIFFID vegetation model.

! River routing:
      INTEGER :: nstep_trip

      REAL :: TOT_SURF_RUNOFF(LAND_POINTS)! Accumulated runoff over
      REAL :: TOT_SUB_RUNOFF(LAND_POINTS) ! river routing timestep
                                         ! (Kg/m2/s)
      REAL :: A_BOXAREAS(ROW_LENGTH, ROWS) ! Atmos.gridbox areas (m2)

      REAL :: RIVEROUT(ROW_LENGTH, ROWS) ! river outflow at
                                         ! seapoints on the ATMOS grid
! Declare local inland basin variables
       REAL                                                             &
     & INLANDOUT_ATMOS(row_length,rows)                                 &
                                                ! INLAND
!                                    BASIN FLOW  ATMOS GRID  kg/m2/s
     &,INLANDOUT_ATM(land_points)         ! INLAND BASIN FLOW
!                                  land points only kg/m2/s

      LOGICAL                                                           &
     & L_INLAND                   ! IN True if re-routing inland
                                  !   basin flow to soil moisture

      REAL :: INLANDOUT_RIV(RIVER_ROW_LENGTH,RIVER_ROWS)
      ! inland basin OUTFLOW on the trip grid

      REAL :: BOX_OUTFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)
                                         ! gridbox outflow on river
                                         ! routing grid (Kg/s)
      REAL :: BOX_INFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)
                                         ! gridbox runoff on river
                                         ! routing grid (Kg/s)

      LOGICAL :: INVERT_ATMOS            ! marker whether ATMOS
!                                        ! grid is inverted N/S or not
      LOGICAL :: TRIP_CALL               !If River routing called
      LOGICAL :: FIRST_ROUTING      !.T. on first call river routing

      real ::                                                           &
     & w_copy(row_length,rows,0:model_levels)  ! copy of w to pass to
                                               ! conv_diag, BL?

! Convection and BL INC diagnostics that
! need to be kept between substeps
! Moved from CONV_CTL2, IMP_CTL2 to enable substepping.
#if !defined(SCMA)
      Real,Dimension(:,:,:),Allocatable ::                              &
     & up_flux_half                                                     &
                                 !up flux on half levels
     &,T_incr_diag_conv                                                 &
                                 ! temperature increment for conv
     &,q_incr_diag_conv                                                 &
                                 ! humidity increment for conv
     &,qcl_incr_diag_conv                                               &
                                 ! qCL   increment for conv
     &,qcf_incr_diag_conv                                               &
                                 ! qCF   increment for conv
     &,cf_liquid_incr_diag_conv                                         &
                                 ! cf_l  increment for conv
     &,cf_frozen_incr_diag_conv                                         &
                                 ! cf_f  increment for conv
     &,bulk_cf_incr_diag_conv                                           &
                                 ! bcf   increment for conv
     &,u_incr_diag_conv                                                 &
                                 ! u wind  increment for conv
     &,v_incr_diag_conv                                                 &
                                 ! v wind  increment for conv
     &,u_incr_diag_bl                                                   &
                                 ! u wind      increment for BL
     &,v_incr_diag_bl                                                   &
                                 ! v wind      increment for BL
     &,T_incr_diag_bl                                                   &
                                 ! temperature increment for BL
     &,q_incr_diag_bl                                                   &
                                 ! humidity    increment for BL
     &,qcl_incr_diag_bl                                                 &
                                 ! cl liq frac increment for BL
     &,qcf_incr_diag_bl                                                 &
                                 ! cl fro frac increment for BL
     &,bulk_cf_incr_diag_bl                                             &
                                 ! cf_l  increment for BL
     &,cf_liquid_incr_diag_bl                                           &
                                 ! cf_f  increment for BL
     &,cf_frozen_incr_diag_bl    ! bcf   increment for BL

      Real,Dimension(:,:,:),Allocatable ::                              &
     &  mf_deep                                                         &
                                 ! mass flux for deep convection
     &, mf_congest                                                      &
                                 ! mass flux for congestus convection
     &, mf_shall                                                        &
                                 ! mass flux for shallow convection
     &, mf_midlev                                                       &
                                 ! mass flux for mid-level convection
     &, dt_deep                                                         &
                                 ! dT for deep convection
     &, dt_congest                                                      &
                                 ! dT for congestus convection
     &, dt_shall                                                        &
                                 ! dT for shallow convection
     &, dt_midlev                                                       &
                                 ! dT for mid-level convection
     &, dq_deep                                                         &
                                 ! dq for deep convection
     &, dq_congest                                                      &
                                 ! dq for congestus convection
     &, dq_shall                                                        &
                                 ! dq for shallow convection
     &, dq_midlev                                                       &
                                 ! dq for mid-level convection
     &, du_deep                                                         &
                                 ! du for deep convection
     &, du_congest                                                      &
                                 ! du for congestus convection
     &, du_shall                                                        &
                                 ! du for shallow convection
     &, du_midlev                                                       &
                                 ! du for mid-level convection
     &, dv_deep                                                         &
                                 ! dv for deep convection
     &, dv_congest                                                      &
                                 ! dv for congestus convection
     &, dv_shall                                                        &
                                 ! dv for shallow convection
     &, dv_midlev                                                       &
                                 ! dv for mid-level convection
     &, wqt_flux_sh                                                     &
                                    ! w'qt flux  shallow
     &, wql_flux_sh                                                     &
                                    ! w'qt flux shallow
     &, wthetal_flux_sh                                                 &
                                    ! w'thetal' flux shallow
     &, wthetav_flux_sh                                                 &
                                    ! w'thetav' flux shallow
     &, dubydt_pout                                                     &
                                 ! du on p grid
     &, dvbydt_pout                                                     &
                                 ! dv on p grid
     &, conv_rain_3d                                                    &
                                 ! 3d conv rainfall rate
     &, conv_snow_3d             ! 3d conv snowfall rate
#else
! Diagnostic output on model levels
! Boundary Layer output arrays
      Real ::                                                           &
     & u_incr_diag_bl(row_length,rows,model_levels)                     &
     &,v_incr_diag_bl(row_length,rows,model_levels)                     &
     &,T_incr_diag_bl(row_length,rows,model_levels)                     &
     &,q_incr_diag_bl(row_length,rows,wet_levels)                       &
     &,qcl_incr_diag_bl(row_length,rows,wet_levels)                     &
     &,qcf_incr_diag_bl(row_length,rows,wet_levels)                     &
     &,bulk_cf_incr_diag_bl(row_length,rows,wet_levels)                 &
     &,cf_liquid_incr_diag_bl(row_length,rows,wet_levels)               &
     &,cf_frozen_incr_diag_bl(row_length,rows,wet_levels)

! Convection output arrays
      Real ::                                                           &
     & T_incr_diag_conv(row_length,rows,model_levels)                   &
     &,u_incr_diag_conv(row_length,rows,model_levels)                   &
     &,v_incr_diag_conv(row_length,rows,model_levels)                   &
     &,q_incr_diag_conv(row_length,rows,wet_levels)                     &
     &,qcl_incr_diag_conv(row_length,rows,wet_levels)                   &
     &,qcf_incr_diag_conv(row_length,rows,wet_levels)                   &
     &,cf_liquid_incr_diag_conv(row_length,rows,wet_levels)             &
     &,cf_frozen_incr_diag_conv(row_length,rows,wet_levels)             &
     &,up_flux_half(row_length,rows,model_levels)                       &
     &,bulk_cf_incr_diag_conv(row_length,rows,wet_levels)

      Real ::                                                           &
     &  mf_deep(1,1,wet_levels)                                         &
     &, mf_congest(1,1,wet_levels)                                      &
     &, mf_shall(1,1,wet_levels)                                        &
     &, mf_midlev(1,1,wet_levels)                                       &
     &, dt_deep(1,1,wet_levels)                                         &
     &, dt_congest(1,1,wet_levels)                                      &
     &, dt_shall(1,1,wet_levels)                                        &
     &, dt_midlev(1,1,wet_levels)                                       &
     &, dq_deep(1,1,wet_levels)                                         &
     &, dq_congest(1,1,wet_levels)                                      &
     &, dq_shall(1,1,wet_levels)                                        &
     &, dq_midlev(1,1,wet_levels)                                       &
     &, du_deep (1,1,wet_levels)                                        &
     &, du_congest(1,1,wet_levels)                                      &
     &, du_shall(1,1,wet_levels)                                        &
     &, du_midlev(1,1,wet_levels)                                       &
     &, dv_deep(1,1,wet_levels)                                         &
     &, dv_congest(1,1,wet_levels)                                      &
     &, dv_shall(1,1,wet_levels)                                        &
     &, dv_midlev(1,1,wet_levels)                                       &
     &, wqt_flux_sh(1,1,wet_levels)                                     &
                                            ! w'qt flux  shallow
     &, wql_flux_sh(1,1,wet_levels)                                     &
                                            ! w'qt flux shallow
     &, wthetal_flux_sh(1,1,wet_levels)                                 &
                                            ! w'thetal' flux shallow
     &, wthetav_flux_sh(1,1,wet_levels)                                 &
                                            ! w'thetav' flux shallow
     &, dubydt_pout(1,1,wet_levels)                                     &
     &, dvbydt_pout(1,1,wet_levels)                                     &
     &, conv_rain_3d(1,1,wet_levels)                                    &
     &, conv_snow_3d(1,1,wet_levels)
#endif

! STASHflag switches for increment diagnostics:
      Logical                                                           &
     & l_u_incr_bl                                                      &
                             ! u wind
     &,l_v_incr_bl                                                      &
                             ! v wind
     &,L_T_incr_bl_lsc                                                  &
                             ! T across BL and LS CLD
     &,L_Tl_incr_bl_lsc                                                 &
                             ! Tl across BL (and LS CLD)
     &,L_q_incr_bl_lsc                                                  &
                             ! Q across BL and LS CLD
     &,L_qtl_incr_bl_lsc                                                &
                             ! QT (q+qCL) across BL (and LS CLD)
     &,L_qcl_incr_bl_lsc                                                &
                             ! qCL across BL and LS CLD
     &,L_qcl_incr_bl                                                    &
                             ! qcl across BL
     &,L_q_incr_bl                                                      &
                             ! q across BL
     &,L_T_incr_bl                                                      &
                             ! T across BL
     &,L_qcf_incr_bl_lsc                                                &
                             ! qCF across BL (and LS CLD)
     &,L_cf_incr_bl                                                     &
                             ! tot cl frac increment for BL
     &,L_cfl_incr_bl                                                    &
                             ! liq cl frac increment for BL
     &,L_cff_incr_bl                                                    &
                             ! ice cl frac increment for BL
     &,L_apply_diag          ! flag to determine when to apply
                             ! diagnostics when iterating

!  Workspace :-

!  Local scalars :-

      REAL                                                              &
     &  MAG_VECTOR_NP (bl_levels)                                       &
     &, DIR_VECTOR_NP (bl_levels)                                       &
     &, MAG_VECTOR_SP (bl_levels)                                       &
     &, DIR_VECTOR_SP (bl_levels)

      Real                                                              &
     & pptn_rate(row_length,rows)                                       &
                                    ! Total precipitation
                                    ! (convective+large scale) (kg/m2/s)
     &,accum_pptn(row_length,rows)                                      &
                                    ! Accumulated total precip (kg/m2)
     &,tot_rain(row_length,rows)                                        &
                                    ! Total rain (kg/m2/s)
     &,tot_snow(row_length,rows)    ! Total snow (kg/m2/s)

      ! end step of experiment
      integer :: endstep, mype

      REAL ltfs

#if defined(ACCESS)
      ltfs = access_tfs
#else
      ltfs = TFS
#endif

! ----------------------------------------------------------------------
! Section INI. Initialisation of variables.
! ----------------------------------------------------------------------
       L_apply_diag = CycleNo == NumCycles

!      ! Copy flag for tracer implicit weights into local variable
       trweights1 = BL_OPTIONS(5)

! Hydrology
#if !defined(SCMA)
      stf_sub_surf_roff = sf(205,8) .or. sf(235,8) .or. l_rivers
      smlt = sf(202,8)
#else
      stf_sub_surf_roff=.true.
      smlt=.true.
#endif

! Set convective outputs to zero in case convection scheme not called.
      conv_rain(:,:) = 0.0
      conv_snow(:,:) = 0.0
! Set radiation outputs to zero in case convection scheme not called.
! This code can be removed later.
! set temporary logicals to disabled un-called physics.
      if(.Not. L_bl .and. convection_option == 3) then
        write(*,*)' convection on and boundary layer off,'
        write(*,*)' not a sensible choice'
        write(*,*)' will try and run but results may be garbage'
      endif

! set p at layer boundaries.
! NB: some arrays have haloes but are unset, if never used they will
!     be removed.
      Do j = 1, rows
        Do i = 1, row_length
          zh_prev(i,j) = zh(i,j)  ! make a copy of zh, as passed in
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
!
! Variables required by conv_diag and convection_control
! Heights of model levels, density; all without halos
!
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            z_theta(i,j,k) = r_theta_levels(i,j,k)                      &
     &                                   - r_theta_levels(i,j,0)
            z_rho(i,j,k)   = r_rho_levels(i,j,k)                        &
     &                                   - r_theta_levels(i,j,0)
            rho_wet(i,j,k) = rho_rsq(i,j,k) /(r_rho_levels(i,j,k)       &
     &                                    *r_rho_levels(i,j,k))
            rho_dry(i,j,k) = unscaled_dry_rho(i,j,k)

          End Do
        End Do
      End Do

! density on theta levels (not top level).

! DEPENDS ON: p_to_t
      call p_to_t(row_length,rows, 0,0,0,0, model_levels-1,             &
     &            r_theta_levels(1:row_length,1:rows,0:model_levels),   &
     &            r_rho_levels(1:row_length,1:rows,1:model_levels),     &
     &            rho_wet, rho_wet_theta)


!
      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            bulk_cloud_fraction(i,j,k) =                                &
     &        bulk_cloud_fraction_halos(i,j,k)
            cloud_fraction_liquid(i,j,k) =                              &
     &        cloud_fraction_liquid_halos(i,j,k)
            cloud_fraction_frozen(i,j,k) =                              &
     &        cloud_fraction_frozen_halos(i,j,k)
          End Do
        End Do
      End Do

      L_poles = .FALSE.
      j_begin = 1
      j_end = rows
      If (model_domain  /=  mt_bi_cyclic_lam) Then
        If (at_extremity(PSouth)) j_begin = 2
        If (at_extremity(PNorth)) j_end = rows-1
      Endif
      If (model_domain  ==  mt_Global)  L_poles = .TRUE.
! zero qcl, qcf and bulk_cloud_fraction on LAM boundaries to avoid
! inconsistencies between cloud fields and cloud amounts.
! As results on boundaries do not affect model results this can be done.
      If (model_domain  ==  mt_LAM) Then

        L_zero_boundaries=.TRUE.
        L_zero_halos=.FALSE.

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QCL,       &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QCF,       &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
         CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,                     &
     &   BULK_CLOUD_FRACTION,                                           &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)
      End If

!     Bug fix to inculde surface currents in coupling. 
!     Do not include if HadGEM1. 
      If (.NOT.L_OHADGEM1) Then 
! SURFACE CURRENT INTERPOLATION START- 
! interpolate U_0,V_0 onto U_0_P,V_0_P 
! as nobody else seems to have thought of it yet... 
 
        UHALO(1:ROW_LENGTH,1:ROWS) = U_0(1:ROW_LENGTH,1:ROWS) 
        VHALO(1:ROW_LENGTH,1:N_ROWS) = V_0(1:ROW_LENGTH,1:N_ROWS) 
 
! Need to call swap bounds as halo points not set 
! DEPENDS ON: swap_bounds 
        CALL swap_bounds(uhalo,row_length,rows,1,                       &
     &                   off_x,off_y,fld_type_u,.true.) 
! DEPENDS ON: swap_bounds 
        CALL swap_bounds(vhalo,row_length,n_rows,1,                     &
     &                   off_x,off_y,fld_type_v,.true.) 
 
! DEPENDS ON: fill_external_halos 
        CALL fill_external_halos(uhalo,ROW_LENGTH,ROWS,1,off_x,off_y) 
! DEPENDS ON: fill_external_halos 
        CALL fill_external_halos(vhalo,ROW_LENGTH,N_ROWS,1,off_x,off_y) 
 
! DEPENDS ON: v_to_p 
        CALL v_to_p(vhalo,row_length,rows,n_rows,1,                     &
     &              off_x,off_y,model_domain,at_extremity,v_0_p) 
 
! DEPENDS ON: u_to_p 
        CALL u_to_p(uhalo,row_length,rows,1,                            &
     &              off_x,off_y,model_domain,at_extremity,u_0_p) 
 
      Else 
        UHALO(1:ROW_LENGTH,1:ROWS) = 0.0 
        VHALO(1:ROW_LENGTH,1:N_ROWS) = 0.0 
      End If 
 
! copy of w or w_adv for use in conv_diag

#if defined(A05_4A) 
! to give bit comparison for HadGEM1 etc currently must continue
! to use w_adv

      Do k=0,model_levels
        Do j=1,rows
          Do i=1,row_length
            w_copy(i,j,k) = w_adv(i,j,k)
          End do
        End do
      End do
#else
! All new convection or no convection scheme switch to using w

      Do k=0,model_levels
        Do j=1,rows
          Do i=1,row_length
            w_copy(i,j,k) = w(i,j,k)
          End do
        End do
      End do
#endif

!
! Calculate etadot
!
#if defined(SCMA)
      do j=1,rows
      do i=1,row_length
      do k=0,model_levels
         Z_TOP_OF_MODEL = z_theta(i,j,model_levels)                     &
     &                  / eta_theta_levels(model_levels)
         ETADOT_COPY(I,J,K) = W(I,J,K)/Z_TOP_OF_MODEL
      ENDDO
      ENDDO
      ENDDO
#else
! Copied from PEHELMEU2A
! DEPENDS ON: etadot_calc
      Call Etadot_Calc (                                                &
     &                  r_theta_levels, r_rho_levels,                   &
     &                  eta_theta_levels, eta_rho_levels,               &
     &                  u, v, w,                                        &
     &                  sec_theta_latitude,                             &
     &                  row_length, rows, n_rows, model_levels,         &
     &                  delta_lambda, delta_phi,                        &
     &                  dlambda_p, dphi_p,                              &
     &                  wt_lambda_p, wt_lambda_u,                       &
     &                  wt_phi_p, wt_phi_v,                             &
     &                  model_domain, first_constant_r_rho_level,       &
     &                  proc_row_group, at_extremity, global_row_length,&
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  off_x, off_y, off_x, off_y, 0, 0,               &
     &                  1, row_length, 1, rows, j_begin, j_end,         &
     &                  L_regular, L_poles,                             &
     &                  L_s_poles, l_n_poles, etadot_copy)

#endif

! ---------------------------------------------
! Call CONV_DIAG to diagnose convection
! ---------------------------------------------


! Set parameters for sub-stepping and begin the loop.
      sub_timestep = timestep/(1.*Num_Substeps)

! ----------------------------------------------------------------------
! Invoke BL/CONV/HYDROL Num_Substep times
! ----------------------------------------------------------------------
      Do Substep_Number=1, Num_Substeps

#if defined(SCMA)
        ! The SCM output diagnostic system needs to know what sub-step
        ! it's currently on, if at all.

! DEPENDS ON: scm_substep_start
        Call SCM_substep_start(Substep_Number)
#endif
!
! latest values needed for substepping/fully sequential BL.
        If ( L_phys2_substep ) Then

          Do k = 1, bl_levels
            Do j = 1, rows
              Do i =1, row_length
                u_conv(i,j,k) = U(i,j,k) + R_u(i,j,k)
              Enddo
            Enddo
          Enddo
! DEPENDS ON: swap_bounds
          Call Swap_Bounds( u_conv, row_length, rows, bl_levels,        &
     &                      off_x, off_y, fld_type_u, .true.  )
! DEPENDS ON: fill_external_halos
          Call FILL_EXTERNAL_HALOS( u_conv,row_length, rows,            &
     &                              bl_levels,off_x,off_y )

          Do k = 1, bl_levels
            Do j = 1, n_rows
              Do i =1, row_length
                V_conv(i,j,k) = V(i,j,k) + R_v(i,j,k)
              Enddo
            Enddo
          Enddo
! DEPENDS ON: swap_bounds
          Call Swap_Bounds( V_conv, row_length, n_rows, bl_levels,      &
     &                      off_x, off_y, fld_type_v, .true.  )
! DEPENDS ON: fill_external_halos
          Call FILL_EXTERNAL_HALOS( V_conv,row_length, n_rows,          &
     &                              bl_levels,off_x,off_y)

          IF(MODEL_DOMAIN  ==  mt_global) THEN

! Overwrite values of U_P, V_P at the poles with the magnitude of
! the vector wind.

! DEPENDS ON: polar_vector_wind_n
            Call Polar_vector_wind_n(                                   &
     &                       v_conv,                                    &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, bl_levels, mag_vector_np,          &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       off_x, off_y, global_row_length,           &
     &                       proc_row_group, at_extremity)

            If (at_extremity(PSouth) ) Then
              DO K=1,BL_LEVELS
                DO I=1,ROW_LENGTH
                  v_conv(I,1,K) = MAG_VECTOR_SP(k)
                  u_conv(I,1,K) = 0.0
                END DO
              End Do
            End If
            If (at_extremity(PNorth) ) Then
              DO K=1,BL_LEVELS
                DO I=1,ROW_LENGTH
                  v_conv(I,rows,K) = MAG_VECTOR_NP(k)
                  u_conv(I,rows,K) = 0.0
                END DO
              End Do
            End If
          ENDIF
! Reset zh_prev when substepping BL
          If ( Substep_Number  >   1 ) Then
            Do j = 1, rows
              Do i = 1, row_length
                zh_prev(i,j) = zh(i,j)
              Enddo
            Enddo
          End If

        Endif ! If ( L_phys2_substep )

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              T_latest(i,j,k) = theta_star(i,j,k)                       &
     &                        * exner_theta_levels(i,j,k)
            End Do
          End Do
        End Do
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              q_latest(i,j,k) = q_star(i,j,k)
              qcl_latest(i,j,k) = qcl_star(i,j,k)
              qcf_latest(i,j,k) = qcf_star(i,j,k)
              cf_latest(i,j,k)  = cf_star(i,j,k)
              cfl_latest(i,j,k) = cfl_star(i,j,k)
              cff_latest(i,j,k) = cff_star(i,j,k)
            End Do
          End Do
        End Do

!  Initialise output arrays
          Do j = 1, rows
            Do i = 1, row_length
              NTML(i,j)=1
              NTPAR(i,j)=1
              NLCL(i,j)=1
              CUMULUS(i,j)=.FALSE.
              L_SHALLOW(i,j)=.FALSE.
              DELTHVU(i,j)=0.0
              ql_ad(i,j) = 0.0
              ZHPAR(i,j)=0.0
              ZLCL(i,j)=0.0
              ZLCL_UV(i,j)=0.0
            End do
          End do


! DEPENDS ON: timer
          If (Ltimer) Call timer ('CONV_DIAG',3)

          If ( .NOT. L_phys2_substep ) Then

! DEPENDS ON: conv_diag
            CALL CONV_DIAG(                                             &

! IN Parallel variables
     &      row_length, rows, n_rows, halo_i, halo_j, off_x, off_y      &
     &,     global_row_length, proc_row_group, at_extremity             &

! IN model dimensions.
     &,     model_domain, bl_levels, model_levels, wet_levels           &
     &,     p, p_layer_centres(1,1,1),exner_rho_levels                  &
     &,     r_rho_levels, r_theta_levels                                &
     &,     rho_wet, rho_wet_theta, z_theta, z_rho                      &
     &,     SIN_THETA_LONGITUDE,COS_THETA_LONGITUDE                     &

! IN Model switches
     &,     l_mixing_ratio                                              &

! IN cloud data
     &,     qcf(1:row_length,1:rows,1:wet_levels)                       &
     &,     qcl(1:row_length,1:rows,1:wet_levels), bulk_cloud_fraction  &

! IN everything not covered so far :

     &,     p_star, q(1:row_length,1:rows,1:wet_levels)                 &
     &,     theta, exner_theta_levels, u, v, u_0_p, v_0_p               &
     &,     L_flux_bc, flux_e, flux_h, L_spec_z0, z0h_scm               &
     &, t_surf, land_sea_mask, timestep                                 &
     &, w_copy                                                          &
!
! SCM Diagnostics (dummy values in full UM)
     &, nSCMDpkgs,L_SCMDiags                                            &

! OUT data required elsewhere in UM system :
     &,     ZH,ZHPAR,ZLCL,ZLCL_UV,DELTHVU,ql_ad,NTML,NTPAR,NLCL,CUMULUS &
     &,     L_SHALLOW,l_congestus,l_congestus2                          &
     &,     CIN_undilute,CAPE_undilute,Error_code                       &
     &     )

          Else

! DEPENDS ON: conv_diag
            CALL CONV_DIAG(                                             &

! IN Parallel variables
     &     row_length, rows, n_rows, halo_i, halo_j, off_x, off_y       &
     &,    global_row_length, proc_row_group, at_extremity              &

! IN model dimensions.
     &,    model_domain, bl_levels, model_levels, wet_levels            &
     &,    p, p_layer_centres(1,1,1), exner_rho_levels                  &
     &,    r_rho_levels, r_theta_levels                                 &
     &,    rho_wet, rho_wet_theta, z_theta, z_rho                       &
     &,    SIN_THETA_LONGITUDE,COS_THETA_LONGITUDE                      &

! IN Model switches
     &,    l_mixing_ratio                                               &

! IN cloud data
     &,    qcf_latest, qcl_latest, bulk_cloud_fraction                  &

! IN everything not covered so far :
     &,    p_star, q_latest, theta_star, exner_theta_levels             &
     &,    u_conv, v_conv, u_0_p, v_0_p                                 &
     &,    L_flux_bc,flux_e, flux_h, L_spec_z0, z0h_scm                 &
     &,    t_surf, land_sea_mask, timestep                              &
     &,    w_copy                                                       &
!
! SCM Diagnostics (dummy values in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! OUT data required elsewhere in UM system :
     &,    ZH,ZHPAR,ZLCL,ZLCL_UV,DELTHVU,ql_ad,NTML,NTPAR,NLCL,CUMULUS  &
     &,    L_SHALLOW,l_congestus,l_congestus2                           &
     &,    CIN_undilute,CAPE_undilute,Error_code                        &
     &      )

          Endif ! if NOT L_phys2_substep

! DEPENDS ON: timer
          If (Ltimer) Call timer ('CONV_DIAG',4)


!        Endif
! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics:
! ----------------------------------------------------------------------
      IF(L_CTILE)THEN

        DO L=1,LAND_POINTS
          FLAND(L)=FLAND_CTILE(L)
        ENDDO
        DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            TSTAR_LAND(I,J)=TSTAR_LAND_CTILE(I,J)
            TSTAR_SEA(I,J)=TSTAR_SEA_CTILE(I,J)
            TSTAR_SICE(I,J)=TSTAR_SICE_CTILE(I,J)
            IF(ICE_FRACT(I,J) <= 0.0)THEN
              TSTAR_SSI(I,J)=TSTAR_SEA(I,J)
            ELSE
              TSTAR_SSI(I,J)=ICE_FRACT(I,J)*TSTAR_SICE(I,J)             &
     &          +(1.0-ICE_FRACT(I,J))*TSTAR_SEA(I,J)
            ENDIF
          ENDDO
        ENDDO

      ELSE

        DO L=1,LAND_POINTS
           FLAND(L)=1.0
        ENDDO
        DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            TSTAR_LAND(I,J)=T_SURF(I,J)
            TSTAR_SSI(I,J)=T_SURF(I,J)

            IF(.NOT.LAND_SEA_MASK(I,J))THEN
              IF(ICE_FRACT(I,J) <= 0.0)THEN
                TSTAR_SEA(I,J)=T_SURF(I,J)
                TSTAR_SICE(I,J)=T_SURF(I,J)
              ELSE
#if defined(ACCESS)
                if (ocn_sss) then
                    ltfs =  ZeroDegC - 0.054 * auscom_salinity(I,J)
                end if
#endif
                TSTAR_SEA(I,J)=LTFS
                TSTAR_SICE(I,J)=(T_SURF(I,J)                            &
     &            -(1.-ICE_FRACT(I,J))*TSTAR_SEA(I,J))/ICE_FRACT(I,J)
              ENDIF
            ELSE
              TSTAR_SEA(I,J)=T_SURF(I,J)
              TSTAR_SICE(I,J)=T_SURF(I,J)
            ENDIF

          ENDDO
        ENDDO
      ENDIF
!----------------------------------------------------------------------
! Section BL Call Explicit part of Boundary Layer scheme.
! ---------------------------------------------------------------------

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Boundary Layer',5)

! create a 2d convective cloud array to pass to the boundary layer
! scheme
      If (lcv_3d_cca) Then
        Do j = 1, rows
          Do i = 1, row_length
            If (ccb(i,j)  >   0) Then
              cca_2d(i,j)= cca(i,j,ccb(i,j))
            Else
              cca_2d(i,j)=0.0
            End If
          End Do
        End Do
      Else
        Do j = 1, rows
          Do i = 1, row_length
            cca_2d(i,j)=cca(i,j,1)
          End Do
        End Do
      End If

      Do j = 1, rows
        Do i = 1, row_length
! Initialise arrays which will be passed from BL to convection
           wstar(i,j) = 0.0
           wthvs(i,j) = 0.0
        End do
      End do

!
! DEPENDS ON: Generate_Anthropogenic_Heat
      Call Generate_Anthropogenic_Heat(                                 &
     & val_year, val_day_number, val_hour, val_minute, val_second       &
     &, ntiles, anthrop_heat, l_anthrop_heat_src                        &
     & )
!
#if !defined(SCMA)
      L_plsp=sf(281,3) .or. sf(282,3) .or. sf(283,3)
      L_scrn=L_plsp
#else
      L_plsp=.true.
      L_scrn=L_plsp
#endif
!
      IF(L_bl)then
! DEPENDS ON: timer
      If (Ltimer) Call timer ('NI_BL_CTL',3)

! DEPENDS ON: ni_bl_ctl
      Call NI_bl_ctl (                                                  &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &

! IN Substep number for ATMPHYS2 and logical switches
     &, Substep_Number, Num_Substeps, L_phys2_substep,NumCycles,CycleNo &

! model dimensions.
     &, row_length, rows, n_rows, land_points, ntiles, model_levels     &
     &, wet_levels, bl_levels, DIM_CS1, DIM_CS2                         &
     &, dst_levels, dsm_levels                                          &

! Model switches
     &, model_domain, formdrag, orog_drag_param, L_CAL360               &
     &, L_area_Cloud, L_mixing_ratio, L_emcorr, Ltimer, L_dry           &
     &, L_scrn, L_plsp, BL_OPTIONS, L_SBLeq, L_SBLco, L_ctile           &
     &, L_DUST, L_CAM_DUST, L_ukca, L_LAMBDAM2, L_FULL_LAMBDAS          &

! model Parameters
     &, ALPHA_CD, CO2_MMR, MUW_SBL, MWT_SBL                             &
     &, Charnock, SeaSalinityFactor                                     &

! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &

! in coordinate information
     &, r_rho_levels, r_theta_levels, r_at_u, r_at_v                    &
     &, eta_theta_levels, eta_rho_levels, delta_lambda                  &
     &, delta_phi, lat_rot_NP, long_rot_NP                              &

! in time stepping information.
     &, sub_timestep, val_year, val_day_number, val_hour, val_minute    &
     &, val_second, timestep_number                                     &

! trig arrays
     &, sin_theta_longitude, cos_theta_longitude, FV_cos_theta_latitude &
     &, sin_theta_latitude                                              &

! in data fields.
     &, p, p_layer_centres, p_layer_boundaries,rho_rsq,rho_wet,rho_dry  &
     &, u, v, w_copy, etadot_copy                                       &
     &, u_conv, v_conv, land_sea_mask, q, qcl, qcf, p_star, theta       &
     &, EXNER_THETA_LEVELS, RAD_HR, MICRO_TENDS, SOIL_LAYER_MOISTURE    &
! variables for subgrid turbulence scheme
     &, visc_BL_m, FM_3D, FH_3D, L_subfilter_vert,L_subfilter_horiz     &
     &, L_subfilter_blend, max_diff,turb_startlev_vert,turb_endlev_vert &
     &, BL_COEF_KM, BL_COEF_KH                                          &

! ancillary fields and fields needed to be kept from timestep to
! timestep
     &, hcon, smvccl, smvcwt, smvcst, sthf, sthu                        &
     &, sil_orog_land, ho2r2_orog, ice_fract                            &
     &, u_0, v_0, u_0_p, v_0_p, land_index                              &
     &, cca_2d, ccb, cct, photosynth_act_rad, surf_radflux              &
     &, SOIL_CLAY,SOIL_SILT,SOIL_SAND,DUST_MREL1,DUST_MREL2             &
     &, DUST_MREL3,DUST_MREL4,DUST_MREL5,DUST_MREL6                     &

! IN additional variables for MOSES II
     &,LAND_PTS_TRIF,NPFT_TRIF                                          &
     &,CANOPY,CATCH,CATCH_SNOW                                          &
     &,SNOW_TILE,Z0_TILE,LW_DOWN                                        &
     &,SW_TILE,TSTAR_TILE                                               &
     &,CO2( 1 : CO2_DIM_LEN, 1 : CO2_DIM_ROW, 1 )                       &
     &,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE                        &
     &,L_PHENOL,L_TRIFFID,L_Q10,ASTEPS_SINCE_TRIFFID,CAN_MODEL          &
     &,CS,FRAC,CANHT_FT,LAI_FT                                          &
     &, FLAND,FLANDG,TSTAR_SEA,ANTHROP_HEAT                             &
     &,ALBSOIL,COS_ZENITH_ANGLE,CAN_RAD_MOD,ILAYERS                     &

!    EAK
!    IN
     &, surf_down_sw,alb_tile,l_tile_pts                                &
!     &, surf_down_sw,alb_tile,cos_zenith_angle               &
     &, ls_rain,ls_snow,SW_DOWN                                         &
     &, lat,long,day,time_sec                                           &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L                   &
     &, SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE             &
     &, TSOIL_TILE,T_SURF_TILE,HCONS,clapp                              &
     &, SATHH,SATCON,HCAP                                               &
     &, SOIL_TYPE,VEG_TYPE                                              &
     &, ISNOW_FLG3L,total_nsteps                                        &
     &, FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB                         &
     &, TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB                     &
     &, USTAR_CAB,SURF_HTF_CAB                                          &
     &, l_cable                                                         & 
!sxy
     &, U_S_CAB,CH_CAB,CD_CAB                                           &
     &, SNAGE_TILE,RTSOIL_TILE                                          &
     &, GFLUX_TILE,SGFLUX_TILE                                          &
     &, CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER                     &
     &, NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE                           &

! in/out
     &, t_soil, ti_gb, t_surf, z0msea, L_spec_z0, z0m_scm, z0h_scm      &
     &, area_cloud_fraction, bulk_cloud_fraction                        &
     &, T_latest, q_latest, qcl_latest, qcf_latest                      &
     &, R_u, R_v, cloud_fraction_liquid, cloud_fraction_frozen          &
     &, zh, zh_prev, flux_e, flux_h, L_flux_bc                          &

! INOUT additional variables for MOSES II
     &,GS,G_LEAF_ACC,NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC                &
     &, TSTAR_LAND,TSTAR_SICE,TSTAR_SSI                                 &

! INOUT Variables for STPH_RP
     &,G0_RP,par_mezcla                                                 &
! diagnostic info
! SCM Diagnostics (dummy values in full UM)
     &, nSCMDpkgs, L_SCMDiags,                                          &
!
! STASH Diagnostics
#include "argsts.h"

! out variables required in IMP_SOLVER
     &  alpha1_sice, ashtf, dtrdz_charney_grid, rdz_charney_grid        &
     &, dtrdz_u, dtrdz_v, rdz_u, rdz_v, cdr10m_u, cdr10m_v, z1_tq       &

! OUT variables which need to be maintained between substeps
     &, RHO_UV,RHO_TQ,DZL_CHARNEY,RDZ                                   &
     &, Z1_UV,Z_FULL,Z_UV,Z_TQ,P_HALF,DELTAP                            &


!START: new to atmos_phys
! out diagnostics (done after implicit solver)
     &, cd, ch, e_sea, fqT                                              &
     &, ftl, h_sea                                                      &
     &, rib_gb                                                          &
     &, taux, tauy, vshr                                                &
     &, zht                                                             &
     &, shallowc, cu_over_orog                                          &
     &, bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6     &
     &, bl_type_7, z0m_eff_gb, z0h_eff_gb, fme                          &

! OUT diagnostics required for soil moisture nudging scheme :
     &, WT_EXT,RA                                                       &

! out data required for tracer mixing :
     &, rho_aresist,aresist,resist_b                                    &

!OUT variables required for mineral dust scheme
     &, R_B_DUST,DUST_FLUX,U_S_T_TILE,U_S_T_DRY_TILE,U_S_STD_TILE       &
     &, KENT, WE_LIM, T_FRAC, ZRZI                                      &
     &, KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC                      &
     &, ZHSC, Z_HALF                                                    &

! OUT data required for 4D-VAR :
     &, rho_cd_modv1, rho_km_var                                        &
!END: new to atmos_phys

! OUT additional variables for MOSES II
     &,FTL_TILE,LE_TILE,RADNET_SICE                                     &
     &,RADNET_TILE,RIB_TILE,RHO_ARESIST_TILE,ARESIST_TILE               &
     &,RESIST_B_TILE,ALPHA1,ASHTF_TILE,FQT_TILE,EPOT_TILE,FQT_ICE       &
     &,FTL_ICE,FRACA,RESFS,RESFT,RHOKH_TILE,RHOKH_SICE                  &
     &,RHOKPM,RHOKPM_POT                                                &
     &,RHOKPM_SICE,Z0HSSI,Z0H_TILE,Z0M_GB,Z0MSSI,Z0M_TILE,CHR1P5M       &
     &,CHR1P5M_SICE,SMC,GPP,NPP,RESP_P,G_LEAF,GPP_FT,NPP_FT             &
     &,RESP_P_FT,RESP_S,RESP_S_TOT,RESP_W_FT                            &
     &,GC,CANHC_TILE,WT_EXT_TILE,FLAKE                                  &
     &,TILE_INDEX,TILE_PTS,TILE_FRAC,FSMC                               &
     &, FLANDG_U,FLANDG_V,rib_ssi                                       &
     &, taux_land,taux_ssi,tauy_land,tauy_ssi,vshr_land,vshr_ssi        &

! now in imp_solv: es_gb, ext, snowmelt
     &, t1_sd, q1_sd, ntml, cumulus, nbdsc, ntdsc                       &
     &, ntpar, nlcl, zhpar, zlcl, l_shallow, wstar, wthvs, delthvu      &
     &, uw0,vw0                                                         &
! now in imp_solv: surf_ht_flux_gb,  snomlt_surf_htf, cH_term
     &, rhokm, rhokm_u, rhokm_v                                         &
!new to atmos_phys
     &, rhokh                                                           &

! error information
     &, Error_code,BL_diag, &
      ! end step of experiment, this step, step width, processor num
     & endstep, mype )

! DEPENDS ON: timer
      If (Ltimer) Call timer ('NI_BL_CTL',4)
      endif

! PC2 required mask
#if defined(A05_5A)
          ! Calculate for each point whether we want to carry
          ! convective cloud information for shallow conv in PC2
      If ( lcv_pc2_diag_sh) then
        ! Carry convective cloud information
        Do j = 1, rows
          Do i = 1, row_length
            l_pc2_diag_sh_pts(i,j) = l_shallow(i,j)
          End Do
        End Do
      Else
        ! Do not carry convective cloud information
        Do j = 1, rows
          Do i = 1, row_length
            l_pc2_diag_sh_pts(i,j) = .false.
          End Do
        End Do
      Endif

#else
          ! Calculate for each point whether we want to carry
          ! convective cloud information for shallow conv in PC2
      Do j = 1, rows
        Do i = 1, row_length
          If (lcv_pc2_diag_sh .and. l_shallow(i,j) .and.                &
                                       bl_type_6(i,j) == 1.0) Then
            ! Carry convective cloud information
            l_pc2_diag_sh_pts(i,j) = .true.
          Else
            ! Do not carry convective cloud information
            l_pc2_diag_sh_pts(i,j) = .false.
          End if
        End Do
      End Do
#endif      
! ----------------------------------------------------------------------
! Section CNV.1 Call Convection scheme.
! ----------------------------------------------------------------------

#if !defined(SCMA)
      If (Convection_option  ==  3) Then
#else
      If (Convection_option  ==  3 .and. conv_mode  <   2) Then
#endif

      If (Error_code  ==  0) Then
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Convection',5)

!
! Set up logicals for condensate increment calculation.
! The PC2 Convection code is implemented in three phases:
!  1. OFF-OFF No extra diagnostic space reserved, no calculations.
!  2.  ON-OFF NEW diagnostics calculated, no existing fields touched.
!  3.  ON-ON  Full prognostic interactions and accompanying diagnostics.
!
! Other parts of the PC2 code use L_PC2 to perform overwriting calcs,
! then L_PC2_RESET to write back the original fields for non-interacting
! mode if required. Hence use of different logicals here.
!
          L_calc_dxek  = ( L_pc2 )
          L_q_interact = ( L_pc2  .AND.  ( .NOT. L_pc2_reset ) )
!
!  Initialise w_max array
          If((cape_opt  ==  3) .or. (cape_opt  ==  4) .or.              &
     &        (cape_opt  ==  5)) then
!w_max only used in these cases
            Do j = 1, rows
              Do i = 1, row_length
                w_max(i,j) = 0.0
              Enddo
            Enddo
!   and find w_max for each column
            Do k =  cape_bottom, cape_top
              Do j = 1, rows
                Do i = 1, row_length
                  if( w_max(i,j) < w_adv(i,j,k) ) then
                      w_max(i,j) = w_adv(i,j,k)
                  endif ! w_max(i,j) < w_adv(i,j,k)
                Enddo
              Enddo
            Enddo !  k =  cape_bottom, cape_top
          End if  !test on cape_opt
!
! Store desired potential temperature in theta_conv array

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                theta_conv(i,j,k) = Theta_star(i,j,k)
                theta_inc(i,j,k) = 0.0
              End Do
            End Do
          End Do
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_conv(i,j,k) = q_star(i,j,k)
                q_inc(i,j,k) = 0.0
                   qcl_conv(i,j,k) = qcl_star(i,j,k)
                qcl_inc(i,j,k)        = 0.0
                   qcf_conv(i,j,k) = qcf_star(i,j,k)
                qcf_inc(i,j,k)        = 0.0
                cf_liquid_conv(i,j,k) = cfl_star(i,j,k)
                cf_liquid_inc(i,j,k)  = 0.0
                cf_frozen_conv(i,j,k) = cff_star(i,j,k)
                cf_frozen_inc(i,j,k)  = 0.0
                bulk_cf_conv(i,j,k)   = cf_star(i,j,k)
                bulk_cf_inc(i,j,k)    = 0.0
              End Do
            End Do
          End Do

#if defined(SCMA)
! If conv_mode = 1 then we want to call convection purely to get
! diagnostics out - initial values will then need to be restored.

      If (conv_mode  ==  1 .and. a_step  >=  ntrad1                     &
     &  .and. mod(a_step-ntrad1,ntrad)  ==  0) then

! DEPENDS ON: save_conv
        Call Save_Conv(row_length, rows, model_levels, wet_levels,      &
     &                 tr_levels, tr_vars,                              &
     &  OFF_X, OFF_Y,                                                   &
     &                 theta_conv, q_conv, qcl_conv, qcf_conv,          &
     &                 cf_liquid_conv, cf_frozen_conv, bulk_cf_conv,    &
     &                 theta_inc, q_inc, qcl_inc, qcf_inc,              &
     &                 cf_liquid_inc, cf_frozen_inc, bulk_cf_inc,       &
     &                 R_U, R_V, AEROSOL,                               &
     &                 DUST_DIV1,DUST_DIV2,DUST_DIV3,                   &
     &                 DUST_DIV4,DUST_DIV5,DUST_DIV6,                   &
     &                 SO2, SO4_AITKEN, SO4_ACCU,                       &
     &                 so4_diss, dms, nh3, soot_new, soot_aged,         &
     &                 soot_cld, ocff_new, ocff_aged, ocff_cld,         &
     &                 co2, free_tracers, OZONE_TRACER,                 &
     &                 cclwp, conv_rain, conv_snow, nconvars,           &
     &                 resdump)

      End If
#endif

! DEPENDS ON: timer
        If (Ltimer) Call timer ('NI_conv_ctl',3)

!---------------------------------------------------------------------
! Intercept values of physics increments before being updated by
! convection for optional output of convective increments
!---------------------------------------------------------------------
#if !defined(SCMA)
      If ( Substep_Number  ==  1 ) Then

        If (((sf(181,5).OR.sf(187,5)) .AND. L_apply_diag) .OR. L_calc_dxek) Then
          Allocate ( T_incr_diag_conv(row_length,rows,model_levels) )
          Do k=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                T_incr_diag_conv(i,j,k) = theta_inc(i,j,k)
              Enddo ! i
            Enddo ! j
          Enddo ! k
        Else
          Allocate ( T_incr_diag_conv(1,1,1) )
        Endif

        If (((sf(182,5).OR.sf(188,5)) .AND. L_apply_diag) .OR. L_calc_dxek) Then
          Allocate (q_incr_diag_conv(row_length,rows,wet_levels))
! Hold input q increment
          Do k=1,wet_levels
            Do j=1,rows
              Do i=1,row_length
                q_incr_diag_conv(i,j,k) = q_inc(i,j,k)
              Enddo ! i
            Enddo ! j
          Enddo ! k
        Else
          Allocate ( q_incr_diag_conv(1,1,1) )
        Endif

        If ( ( sf(183,5) .AND. L_apply_diag ) .OR. L_calc_dxek ) Then
          Allocate (qcl_incr_diag_conv(row_length,rows,wet_levels))
!
! Hold input qcl increment
          Do k=1,wet_levels
            Do j=1,rows
              Do i=1,row_length
                qcl_incr_diag_conv(i,j,k) = qcl_inc(i,j,k)
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate (qcl_incr_diag_conv(1,1,1))
        Endif                   ! on STASHflag
!
        If ( ( sf(184,5) .AND. L_apply_diag ) .OR. L_calc_dxek ) Then
          Allocate (qcf_incr_diag_conv(row_length,rows,wet_levels))
!
! Hold input qcf increment
          Do k=1,wet_levels
            Do j=1,rows
              Do i=1,row_length
                qcf_incr_diag_conv(i,j,k) = qcf_inc(i,j,k)
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate (qcf_incr_diag_conv(1,1,1))
        Endif                   ! on STASHflag
!
        If ( ( sf(193,5) .AND. L_apply_diag ) .OR. L_calc_dxek ) Then
          Allocate(cf_liquid_incr_diag_conv(row_length,rows,            &
     &             wet_levels))
!
! Hold input cf_liquid increment
          Do k=1,wet_levels
            Do j=1,rows
              Do i=1,row_length
                cf_liquid_incr_diag_conv(i,j,k) = cf_liquid_inc(i,j,k)
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate(cf_liquid_incr_diag_conv(1,1,1))
        Endif                   ! on STASHflag
!
        If ( ( sf(194,5) .AND. L_apply_diag ) .OR. L_calc_dxek ) Then
          Allocate(cf_frozen_incr_diag_conv(row_length,rows,            &
     &             wet_levels))
!
! Hold input cf_frozen increment
          Do k=1,wet_levels
            Do j=1,rows
              Do i=1,row_length
                cf_frozen_incr_diag_conv(i,j,k) = cf_frozen_inc(i,j,k)
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate(cf_frozen_incr_diag_conv(1,1,1))
        Endif                   ! on STASHflag
!
        If ( ( sf(195,5) .AND. L_apply_diag ) .OR. L_calc_dxek ) Then
          Allocate(bulk_cf_incr_diag_conv(row_length,rows,wet_levels))
!
! Hold input bulk_cf increment
          Do k=1,wet_levels
            Do j=1,rows
              Do i=1,row_length
                bulk_cf_incr_diag_conv(i,j,k) = bulk_cf_inc(i,j,k)
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate(bulk_cf_incr_diag_conv(1,1,1))
        Endif                   ! on STASHflag

        If ( sf(185,5) .AND. L_apply_diag ) Then
          Allocate ( u_incr_diag_conv(row_length,rows,model_levels) )
          Do k=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                u_incr_diag_conv(i,j,k) = 0.0
              Enddo ! i
            Enddo ! j
          Enddo ! k
        Else
          Allocate ( u_incr_diag_conv(1,1,1) )
        Endif

        If ( sf(186,5) .AND. L_apply_diag ) Then
          Allocate ( v_incr_diag_conv(row_length,n_rows,model_levels) )
          Do k=1,model_levels
            Do j=1,n_rows
              Do i=1,row_length
                v_incr_diag_conv(i,j,k) = 0.0
              Enddo ! i
            Enddo ! j
          Enddo ! k
        Else
          Allocate ( v_incr_diag_conv(1,1,1) )
        Endif


! Other model level diagnostics :
! setup allocatable arrays to limit memory over head if diagnostics not
! required.

        If ( sf(227,5) ) then
          Allocate ( conv_rain_3d(row_length,rows,wet_levels) )
        Else
          Allocate ( conv_rain_3d(1,1,1) )
        Endif
        If ( sf(228,5) ) then
          Allocate ( conv_snow_3d(row_length,rows,wet_levels) )
        Else
          Allocate ( conv_snow_3d(1,1,1) )
        Endif
        If ( sf(290,5) ) then
          Allocate ( wqt_flux_sh(row_length,rows,wet_levels) )
        Else
          Allocate ( wqt_flux_sh(1,1,1) )
        Endif
        If ( sf(291,5) ) then
          Allocate ( wql_flux_sh(row_length,rows,wet_levels) )
        Else
          Allocate ( wql_flux_sh(1,1,1) )
        Endif
        If ( sf(292,5) ) then
          Allocate ( wthetal_flux_sh(row_length,rows,wet_levels) )
        Else
          Allocate ( wthetal_flux_sh(1,1,1) )
        Endif
        If ( sf(293,5) ) then
          Allocate ( wthetav_flux_sh(row_length,rows,wet_levels) )
        Else
          Allocate ( wthetav_flux_sh(1,1,1) )
        Endif

        If ( sf(256,5) ) then
          Allocate ( dubydt_pout(row_length,rows,wet_levels) )
        Else
          Allocate ( dubydt_pout(1,1,1) )
        Endif
        If ( sf(257,5) ) then
          Allocate ( dvbydt_pout(row_length,rows,wet_levels) )
        Else
          Allocate ( dvbydt_pout(1,1,1) )
        Endif

        If ( sf(320,5) ) then
          Allocate ( mf_deep(row_length,rows,wet_levels) )
        Else
          Allocate ( mf_deep(1,1,1) )
        Endif
        If ( sf(321,5) ) then
          Allocate ( mf_congest(row_length,rows,wet_levels) )
        Else
          Allocate ( mf_congest(1,1,1) )
        Endif
        If ( sf(322,5) ) then
          Allocate ( mf_shall(row_length,rows,wet_levels) )
        Else
          Allocate ( mf_shall(1,1,1) )
        Endif
        If ( sf(323,5) ) then
          Allocate ( mf_midlev(row_length,rows,wet_levels) )
        Else
          Allocate ( mf_midlev(1,1,1) )
        Endif
        If ( sf(324,5) ) then
          Allocate ( dt_deep(row_length,rows,wet_levels) )
        Else
          Allocate ( dt_deep(1,1,1) )
        Endif
        If ( sf(325,5) ) then
          Allocate ( dt_congest(row_length,rows,wet_levels) )
        Else
          Allocate ( dt_congest(1,1,1) )
        Endif
        If ( sf(326,5) ) then
          Allocate ( dt_shall(row_length,rows,wet_levels) )
        Else
          Allocate ( dt_shall(1,1,1) )
        Endif
        If ( sf(327,5) ) then
          Allocate ( dt_midlev(row_length,rows,wet_levels) )
        Else
          Allocate ( dt_midlev(1,1,1) )
        Endif
        If ( sf(328,5) ) then
          Allocate ( dq_deep(row_length,rows,wet_levels) )
        Else
          Allocate ( dq_deep(1,1,1) )
        Endif
        If ( sf(329,5) ) then
          Allocate ( dq_congest(row_length,rows,wet_levels) )
        Else
          Allocate ( dq_congest(1,1,1) )
        Endif
        If ( sf(330,5) ) then
          Allocate ( dq_shall(row_length,rows,wet_levels) )
        Else
          Allocate ( dq_shall(1,1,1) )
        Endif
        If ( sf(331,5) ) then
          Allocate ( dq_midlev(row_length,rows,wet_levels) )
        Else
          Allocate ( dq_midlev(1,1,1) )
        Endif
        If ( sf(332,5) ) then
          Allocate ( du_deep(row_length,rows,wet_levels) )
        Else
          Allocate ( du_deep(1,1,1) )
        Endif
        If ( sf(333,5) ) then
          Allocate ( du_congest(row_length,rows,wet_levels) )
        Else
          Allocate ( du_congest(1,1,1) )
        Endif
        If ( sf(334,5) ) then
          Allocate ( du_shall(row_length,rows,wet_levels) )
        Else
          Allocate ( du_shall(1,1,1) )
        Endif
        If ( sf(335,5) ) then
          Allocate ( du_midlev(row_length,rows,wet_levels) )
        Else
          Allocate ( du_midlev(1,1,1) )
        Endif
        If ( sf(336,5) ) then
          Allocate ( dv_deep(row_length,rows,wet_levels) )
        Else
          Allocate ( dv_deep(1,1,1) )
        Endif
        If ( sf(337,5) ) then
          Allocate ( dv_congest(row_length,rows,wet_levels) )
        Else
          Allocate ( dv_congest(1,1,1) )
        Endif
        If ( sf(338,5) ) then
          Allocate ( dv_shall(row_length,rows,wet_levels) )
        Else
          Allocate ( dv_shall(1,1,1) )
        Endif
        If ( sf(339,5) ) then
          Allocate ( dv_midlev(row_length,rows,wet_levels) )
        Else
          Allocate ( dv_midlev(1,1,1) )
        Endif


        If ( sf(249,5) .OR. sf(246,5) ) Then
          Allocate ( up_flux_half(row_length,rows,model_levels) )
          Do k=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                up_flux_half(i,j,k) = 0.0
              Enddo ! i
            Enddo ! j
          Enddo ! k
        Else
          Allocate ( up_flux_half(1,1,1) )
        Endif

      Endif  ! if ( Substep_Number  ==  1 )
#else
!
! SCM need to initialise diagnostic space
          Do k=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                T_incr_diag_conv(i,j,k) = theta_inc(i,j,k)
                u_incr_diag_conv(i,j,k) = 0.0
                v_incr_diag_conv(i,j,k) = 0.0
              End Do ! i
            End Do ! j
          End Do ! k

          Do k=1,wet_levels
            Do j=1,rows
              Do i=1,row_length
                q_incr_diag_conv(i,j,k)   = q_inc(i,j,k)
                qcl_incr_diag_conv(i,j,k) = qcl_inc(i,j,k)
                qcf_incr_diag_conv(i,j,k) = qcf_inc(i,j,k)
                cf_liquid_incr_diag_conv(i,j,k) = cf_liquid_inc(i,j,k)
                cf_frozen_incr_diag_conv(i,j,k) = cf_frozen_inc(i,j,k)
                bulk_cf_incr_diag_conv(i,j,k)   = bulk_cf_inc(i,j,k)
              End Do ! i
            End Do ! j
          End Do ! k
#endif

        If (lcv_ccrad) Then

          ! Copy ccw/lcbase values for Radiation if lcv_ccrad in use.
          ! See comments at declaration of ccw_out/lcbase_out

          Do k=1,wet_levels
            Do j=1,rows
              Do i=1,row_length
                ccw(i,j,k) = ccw_out(i,j,k)
              End Do
            End Do
          End Do

          Do j=1,rows
            Do i=1,row_length
              lcbase(i,j)  = lcbase_out(i,j)
            End Do
          End Do

        End If ! lcv_ccrad

! DEPENDS ON: ni_conv_ctl
        Call NI_conv_ctl (                                              &
! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, at_extremity, n_proc, n_procx, n_procy          &
     &, neighbour, g_rows, g_row_length, g_datastart, me                &
! parameters for cycling physics-dynamics
     &, NumCycles, CycleNo                                              &
     &, Substep_Number, Num_Substeps                                    &
! model dimensions.
     &, row_length, rows, n_rows                                        &
     &, rows*row_length                                                 &
     &, model_levels, wet_levels, bl_levels, n_cca_levels               &
     &, tr_levels, tr_vars, tr_ukca                                     &

! Model switches
     &, model_domain, L_CAL360                                          & 
     &, L_calc_dxek, L_q_interact, l_mixing_ratio                       &
     &, LTIMER, L_MURK_SOURCE, L_DUST, L_SULPC_SO2                      &
     &, L_sulpc_dms, L_sulpc_nh3, L_soot, L_ocff, L_biomass             &
     &, L_co2_interactive, L_USE_CARIOLLE                               &

! Model parameters
     &, dbsdtbs_turb_0                                                  &

! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &
! in coordinate information
     &, r_rho_levels, r_theta_levels ,z_rho, z_theta                    &
     &, eta_theta_levels, eta_rho_levels                                &
     &, delta_lambda, delta_phi                                         &
     &, lat_rot_NP, long_rot_NP                                         &
! in time stepping information.
     &, sub_timestep                                                    &
     &, val_year, val_day_number, val_hour, val_minute                  &
     &, val_second, timestep_number,                                    &
! diagnostic info
#include "argsts.h"
     & STASHwork5                                                       &
     &, ls_rain, ls_snow                                                &
!
! SCM diagnostics (dummy in full UM)
     &, nSCMDpkgs,L_SCMDiags                                            &

! in data fields.
     &, rho_rsq, rho_wet, rho_wet_theta                                 &
     &, u, v, p, p_star                                                 &
     &, land_sea_mask                                                   &
     &, p_layer_boundaries, p_layer_centres                             &
     &, exner_layer_boundaries, exner_layer_centres                     &
     &, t1_sd, q1_sd, exner_theta_levels                                &
     &, uw0, vw0, w_max, zlcl, zlcl_uv, zhpar                           &
     &, cumulus, l_shallow, l_congestus, l_congestus2,l_pc2_diag_sh_pts &
     &, ntml, ntpar, ntdsc, nbdsc                                       &
     &, wstar, wthvs, delthvu, ql_ad, ftl ,fqt                          &
     &, shallowc, cu_over_orog, cape_undilute, cin_undilute             &

! IN for STPH_SCV
     &, f3_at_u, FV_cos_theta_latitude, cos_theta_longitude             &
! in/out stash diagnostics
     &, up_flux, up_flux_half, dwn_flux, entrain_up, entrain_dwn        &
     &, detrain_up, detrain_dwn, uw_dp, vw_dp, uw_shall, vw_shall       &
     &, wqt_flux_sh,wthetal_flux_sh,wthetav_flux_sh,wql_flux_sh         &
     &, wqt_cb,wthetal_cb,wqt_inv,wthetal_inv,sh_top,sh_base            &
     &, T_incr_diag_conv, q_incr_diag_conv                              &
     &, qcl_incr_diag_conv, qcf_incr_diag_conv                          &
     &, u_incr_diag_conv, v_incr_diag_conv                              &
     &, cf_liquid_incr_diag_conv, cf_frozen_incr_diag_conv              &
     &, bulk_cf_incr_diag_conv                                          &
     &, dubydt_pout,dvbydt_pout                                         &
     &, precip_deep, precip_shall, precip_mid,precip_cong, cape         &
     &, shallow_ind, congestus_ind,congestus_ind2                       &
     &, mid_ind, ntml_diag, ntpar_diag, freeze_diag                     &
     &, kterm_diag, wstar_up_diag, wstar_dn_diag, mb1_diag, mb2_diag    &
     &, cg_term,cg_top,cg_base                                          &
     &, mf_deep,mf_congest,mf_shall,mf_midlev                           &
     &, dt_deep,dt_congest,dt_shall,dt_midlev                           &
     &, dq_deep,dq_congest,dq_shall,dq_midlev                           &
     &, du_deep,du_congest,du_shall,du_midlev                           &
     &, dv_deep,dv_congest,dv_shall,dv_midlev                           &
! in/out
     &, theta_conv,q_conv,qcl_conv,qcf_conv,cf_liquid_conv              &
     &, cf_frozen_conv, bulk_cf_conv                                    &
     &, theta_inc, q_inc, qcl_inc, qcf_inc, cf_liquid_inc, cf_frozen_inc&
     &, bulk_cf_inc                                                     &
     &, R_U, R_V, AEROSOL                                               &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, CONSCAV_DUST, SO2, SO4_AITKEN, SO4_ACCU, SO4_DISS               &
     &, dms, nh3, soot_new, soot_aged, soot_cld, bmass_new, bmass_aged  &
     &, bmass_cld, ocff_new, ocff_aged, ocff_cld, co2, conscav_so4ait   &
     &, conscav_so4acc, conscav_so4dis, conscav_agedsoot                &
     &, conscav_agedbmass, conscav_agedocff, free_tracers, ukca_tracers &
     &, OZONE_TRACER                                                    &

! out fields
     &, cca, ccb, cct, cclwp, conv_rain, conv_snow                      &
     &, conv_rain_3d, conv_snow_3d                                      &
     &, ccw, lcca, lcclwp, lctop, lcbase                                &

! error information
     &, Error_code  )

      If (lcv_ccrad) Then

        ! Copy ccw/lcbase values for Radiation if lcv_ccrad in use.
        ! See comments at declaration of ccw_out/lcbase_out

        Do k=1,wet_levels
          Do j=1,rows
            Do i=1,row_length
              ccw_out(i,j,k) = ccw(i,j,k)
            End Do
          End Do
        End Do

        Do j=1,rows
          Do i=1,row_length
            lcbase_out(i,j)  = lcbase(i,j)
          End Do
        End Do

      End If ! lcv_ccrad

#if !defined(SCMA)
      If ( Substep_Number  ==  Num_Substeps ) Then
        Deallocate (dubydt_pout)
        Deallocate (dvbydt_pout)
        Deallocate (wqt_flux_sh)
        Deallocate (wql_flux_sh)
        Deallocate (wthetal_flux_sh)
        Deallocate (wthetav_flux_sh)
        Deallocate (mf_deep)
        Deallocate (mf_congest)
        Deallocate (mf_shall)
        Deallocate (mf_midlev)
        Deallocate (dt_deep)
        Deallocate (dt_congest)
        Deallocate (dt_shall)
        Deallocate (dt_midlev)
        Deallocate (dq_deep)
        Deallocate (dq_congest)
        Deallocate (dq_shall)
        Deallocate (dq_midlev)
        Deallocate (du_deep)
        Deallocate (du_congest)
        Deallocate (du_shall)
        Deallocate (du_midlev)
        Deallocate (dv_deep)
        Deallocate (dv_congest)
        Deallocate (dv_shall)
        Deallocate (dv_midlev)
        Deallocate(up_flux_half)
        Deallocate(conv_rain_3d)
        Deallocate(conv_snow_3d)

      Endif
      If ( Substep_Number  ==  Num_Substeps ) Then
! Clear up allocatable arrays
        Deallocate (T_incr_diag_conv)
        Deallocate (q_incr_diag_conv)
        Deallocate(qcl_incr_diag_conv)
        Deallocate(qcf_incr_diag_conv)
        Deallocate(u_incr_diag_conv)
        Deallocate(v_incr_diag_conv)
        Deallocate(cf_liquid_incr_diag_conv)
        Deallocate(cf_frozen_incr_diag_conv)
        Deallocate(bulk_cf_incr_diag_conv)
      Endif
#endif

! DEPENDS ON: timer
        If (Ltimer) Call timer ('NI_conv_ctl',4)

#if defined(SCMA)
! If conv_mode = 1 then we want to restore all the values before
! convection except for cca which we want to use in the diagnostics

      If (conv_mode  ==  1 .and. a_step  >=  ntrad1                     &
     &  .and. mod(a_step-ntrad1,ntrad)  ==  0) then

! DEPENDS ON: restore_conv
        Call Restore_Conv(resdump, nconvars,                            &
     &                 row_length, rows, model_levels, wet_levels,      &
     &                 tr_levels, tr_vars,                              &
     &  OFF_X, OFF_Y,                                                   &
     &                 theta_conv, q_conv, qcl_conv, qcf_conv,          &
     &                 cf_liquid_conv, cf_frozen_conv, bulk_cf_conv,    &
     &                 theta_inc, q_inc, qcl_inc, qcf_inc,              &
     &                 cf_liquid_inc, cf_frozen_inc, bulk_cf_inc,       &
     &                 R_U,R_V,AEROSOL,DUST_DIV1,DUST_DIV2,             &
     &                 DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6,         &
     &                 SO2,SO4_AITKEN,SO4_ACCU,                         &
     &                 so4_diss, dms, nh3, soot_new, soot_aged,         &
     &                 soot_cld, ocff_new, ocff_aged, ocff_cld,         &
     &                 co2, free_tracers, OZONE_TRACER,                 &
     &                 cclwp, conv_rain, conv_snow)

      End If

      If(conv_mode  <   2 .and. a_step  >=  ntrad1                      &
     &  .and. mod(a_step-ntrad1,ntrad)  ==  0) then

! convert theta to t

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              t(i,j,k) = theta(i,j,k)*exner_theta_levels(i,j,k)
            End Do
          End Do
        End Do

        If (conv_mode == 1) then
! DEPENDS ON: sub_data
          Call SUB_DATA(                                                &
     &      row_length, rows, model_levels, wet_levels,                 &
     &      nfor, bl_levels, dst_levels, dsm_levels, ntrop,             &
     &      ' After convect, before hydr        ', a_step,              &
     &      val_year, val_day_number, time_string, daycount,            &
     &      u, v, t, theta, q, qcl, qcf, layer_cloud, p, rho_rsq,       &
     &      exner_rho_levels, t_soil, smc,canopy, snow_depth,           &
     &      tstar, zh, z0msea, cca, ccb, cct, soil_layer_moisture)
        End if ! conv_mode

      End If
#endif
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Convection',6)

      End If ! on error code equal to zero

      End If ! on convection option

! The following needs to be executed regardless if the conv scheme
! is used.
! Restore thermo vars to their value before convection as they are
! needed in the BL for the calculation of NT terms and increments
! when L_phys2_substep=T
      If ( L_phys2_substep ) Then
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              q_conv(i,j,k)     = q_star(i,j,k)
              qcl_conv(i,j,k)   = qcl_star(i,j,k)
              qcf_conv(i,j,k)   = qcf_star(i,j,k)
! restore theta value and convert to temperature as BL works with T
              theta_conv(i,j,k) = theta_star(i,j,k)                     &
     &                          * exner_theta_levels(i,j,k)
            Enddo
          Enddo
        Enddo
        Do k = wet_levels+1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              theta_conv(i,j,k) = theta_star(i,j,k)                     &
     &                          * exner_theta_levels(i,j,k)
            Enddo
          Enddo
        Enddo
      Endif

#if !defined(SCMA)
      If (Convection_option  ==  3) Then
#else
      If (Convection_option  ==  3 .and. conv_mode  <   2) Then
#endif

      If ( Error_code  ==  0 ) Then

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Convection',5)

! update star variables by adding on increments
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_star(i,j,k) = theta_star(i,j,k)                       &
     &                        + theta_inc(i,j,k)
          End Do
        End Do
      End Do
      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            q_star(i,j,k) = q_star(i,j,k) + q_inc(i,j,k)
          End Do
        End Do
      End Do
! ----------------------------------------------------------------------
! Protected loop. Update increments only when PC2 scheme is fully ON-ON.
! ----------------------------------------------------------------------
!
! L_calc_dxek_if1:
      If (L_calc_dxek .AND. L_q_interact) Then
!
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              qcl_star(i,j,k) = qcl_star(i,j,k) + qcl_inc(i,j,k)
              qcf_star(i,j,k) = qcf_star(i,j,k) + qcf_inc(i,j,k)
              cfl_star(i,j,k) = cfl_star(i,j,k) + cf_liquid_inc(i,j,k)
              cff_star(i,j,k) = cff_star(i,j,k) + cf_frozen_inc(i,j,k)
              cf_star(i,j,k)  = cf_star(i,j,k)  + bulk_cf_inc(i,j,k)
            End Do
          End Do
        End Do
!
      End If  ! L_calc_dxek_if1
!
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Convection',6)

      Endif ! Error_code

! ----------------------------------------------------------------------
! Section CNV.2 Energy correction code
! ----------------------------------------------------------------------

      If ( CycleNo == NumCycles .and.                                   &
     &     Substep_Number == Num_Substeps .and.                         &
     &     L_emcorr .and. Error_code == 0 )  Then

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Conv Eng Corr',5)

! Add convective rain and snow, at the surface to the
! diabatic heating for use in the energy correction
! procedure.
! Scale variables by conversion factor so that only one call is required

        lclf = lc + lf
        Do j = 1, rows
          Do i = 1, row_length
            tot_precip_scaled(i,j) =  conv_rain(i,j)                    &
     &                               * lc +                             &
     &                                conv_snow(i,j)                    &
     &                               * lclf
          End Do
        End Do

! DEPENDS ON: timer
      If (Ltimer) Call timer ('flux_diag',3)

! DEPENDS ON: flux_diag
        Call flux_diag(tot_precip_scaled, FV_cos_theta_latitude,        &
     &                 row_length, rows ,off_x,off_y, 1.0,              &
     &                  sum_eng_fluxes,    timestep)

! DEPENDS ON: timer
      If (Ltimer) Call timer ('flux_diag',4)
        Do j = 1, rows
          Do i = 1, row_length
            tot_precip_scaled(i,j) = -conv_rain(i,j)-conv_snow(i,j)
          End Do
        End Do

! DEPENDS ON: timer
      If (Ltimer) Call timer ('flux_diag',3)

! DEPENDS ON: flux_diag
        Call flux_diag(tot_precip_scaled, FV_cos_theta_latitude,        &
     &                 row_length, rows ,off_x,off_y, 1.0,              &
     &                  sum_moist_flux,    timestep)

! DEPENDS ON: timer
      If (Ltimer) Call timer ('flux_diag',4)

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Conv Eng Corr',6)

      End If   ! L_emcorr .and. Substep_Number  ==  Num_Substeps

      End If ! on convection option

! ----------------------------------------------------------------------
! Section BL Call implicit solver
! ---------------------------------------------------------------------

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            T_latest(i,j,k) = theta_star(i,j,k)                         &
     &                      * exner_theta_levels(i,j,k)
          End Do
        End Do
      End Do
      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            q_latest(i,j,k) = q_star(i,j,k)
            qcl_latest(i,j,k) = qcl_star(i,j,k)
            qcf_latest(i,j,k) = qcf_star(i,j,k)
            cf_latest(i,j,k)  = cf_star(i,j,k)
            cfl_latest(i,j,k) = cfl_star(i,j,k)
            cff_latest(i,j,k) = cff_star(i,j,k)
          End Do
        End Do
      End Do
!
! create a 2d convective cloud array to pass to the boundary layer
! scheme, for the purposes of calculating the visibility in precip.
      If (lcv_3d_cca) Then
        Do j = 1, rows
          Do i = 1, row_length
            If (ccb(i,j)  >   0) Then
              cca_2d(i,j)= cca(i,j,ccb(i,j))
            Else
              cca_2d(i,j)=0.0
            End If
          End Do
        End Do
      Else
        Do j = 1, rows
          Do i = 1, row_length
            cca_2d(i,j)=cca(i,j,1)
          End Do
        End Do
      End If
!
      If ( Error_code  ==  0 ) Then

      IF (L_bl) Then

! DEPENDS ON: timer
      If (Ltimer) Call timer ('NI_IMP_CTL',3)
!
!
#if !defined(SCMA)
!---------------------------------------------------------------------
! Intercept values of physics increments before being updated by
! implicit solver for optional output of bl increments
!---------------------------------------------------------------------
! [ Note: u/v_incr_diag_bl has been redefined over bl_levels to
!   saved space, since upper level increments should be all zero. ]

      If ( Substep_Number  ==  1 ) Then
! Set diagnostic flags required for boundary layer diagnostics from
! STASHflags.
        L_u_incr_bl = sf(185,3) .and. L_apply_diag
        L_v_incr_bl = sf(186,3) .and. L_apply_diag
        L_T_incr_bl_lsc = sf(181,9) .and. L_apply_diag
        L_Tl_incr_bl_lsc = sf(189,3) .and. L_apply_diag
        L_q_incr_bl_lsc = sf(182,9) .and. L_apply_diag
        L_qtl_incr_bl_lsc = sf(190,3) .and. L_apply_diag
        L_qcl_incr_bl_lsc = sf(183,9) .and. L_apply_diag
        L_qcf_incr_bl_lsc = sf(184,3) .and. L_apply_diag
        L_cf_incr_bl  = sf(192,3) .and. L_apply_diag
        L_cfl_incr_bl = sf(193,3) .and. L_apply_diag
        L_cff_incr_bl = sf(194,3) .and. L_apply_diag
        L_qcl_incr_bl = sf(183,3) .and. L_apply_diag
        L_q_incr_bl   = sf(182,3) .and. L_apply_diag
        L_T_incr_bl   = sf(181,3) .and. L_apply_diag
!
! Allocate and initialize u,v,T,Q,Qcl,Qcf increments.
! Diagnostics for u,v,T,Q,Qcl,Qcf increments need to retain their values
! between calls, to allow their summation over the number of substeps,
! and thus had to be moved out of ni_imp_ctl().
!
        If ( L_u_incr_bl ) Then
          Allocate ( u_incr_diag_bl(row_length,rows,model_levels) )
          Do k=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                u_incr_diag_bl(i,j,k) = 0.0
              Enddo ! i
            Enddo ! j
          Enddo ! k
        Else
          Allocate ( u_incr_diag_bl(1,1,1) )
        Endif                   ! on STASHflag

        If ( L_v_incr_bl ) Then
          Allocate ( v_incr_diag_bl(row_length,n_rows,model_levels) )
          Do k=1,model_levels
            Do j=1,n_rows
              Do i=1,row_length
                v_incr_diag_bl(i,j,k) = 0.0
              Enddo ! i
            Enddo ! j
          Enddo ! k
        Else
          Allocate ( v_incr_diag_bl(1,1,1) )
        Endif                   ! on STASHflag

        If ( L_T_incr_bl_lsc .OR. L_Tl_incr_bl_lsc                      &
     &      .or. L_T_incr_bl ) Then
          Allocate ( T_incr_diag_bl(row_length,rows,model_levels) )
          Do k=1, model_levels
            Do j=1, rows
              Do i=1, row_length
                T_incr_diag_bl(i,j,k) = 0.0
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate ( T_incr_diag_bl(1,1,1) )
        End if                   ! on STASHflags

        If ( L_q_incr_bl_lsc .OR. L_qtl_incr_bl_lsc                     &
     &      .or. L_q_incr_bl ) Then
          Allocate ( q_incr_diag_bl(row_length,rows,wet_levels) )
          Do k=1, wet_levels
            Do j=1, rows
              Do i=1, row_length
                q_incr_diag_bl(i,j,k) = 0.0
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate ( q_incr_diag_bl(1,1,1) )
        End if                  ! on STASHflags

        If ( L_qcl_incr_bl_lsc .OR. L_Tl_incr_bl_lsc .OR.               &
     &       L_qtl_incr_bl_lsc .or. L_qcl_incr_bl) Then
          Allocate ( qcl_incr_diag_bl(row_length,rows,wet_levels) )
          Do k=1, wet_levels
            Do j=1, rows
              Do i=1, row_length
                qcl_incr_diag_bl(i,j,k) = 0.0
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate ( qcl_incr_diag_bl(1,1,1) )
        End if                  ! on STASHflag

        If ( L_qcf_incr_bl_lsc ) Then
          Allocate (qcf_incr_diag_bl(row_length,rows,wet_levels))
          Do k=1, wet_levels
            Do j=1, rows
              Do i=1, row_length
                qcf_incr_diag_bl(i,j,k) = 0.0
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate ( qcf_incr_diag_bl(1,1,1) )
        End if                  ! on STASHflag!

        If ( L_cf_incr_bl ) Then
          Allocate (bulk_cf_incr_diag_bl(row_length,rows,wet_levels))
          Do k=1, wet_levels
            Do j=1, rows
              Do i=1, row_length
                bulk_cf_incr_diag_bl(i,j,k) = 0.0
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate (bulk_cf_incr_diag_bl(1,1,1))
        Endif

        If ( L_cfl_incr_bl ) Then
          Allocate (cf_liquid_incr_diag_bl(row_length,rows,wet_levels))
          Do k=1, wet_levels
            Do j=1, rows
              Do i=1, row_length
                cf_liquid_incr_diag_bl(i,j,k) = 0.0
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate (cf_liquid_incr_diag_bl(1,1,1))
        Endif

        If ( L_cff_incr_bl ) Then
          Allocate (cf_frozen_incr_diag_bl(row_length,rows,wet_levels))
          Do k=1, wet_levels
            Do j=1, rows
              Do i=1, row_length
                cf_frozen_incr_diag_bl(i,j,k) = 0.0
              End Do ! i
            End Do ! j
          End Do ! k
        Else
          Allocate (cf_frozen_incr_diag_bl(1,1,1))
        Endif

      Endif ! if ( Substep_Number  ==  1 )
#else
!
! SCM set logicals for diagnostics
        L_u_incr_bl       = .TRUE.
        L_v_incr_bl       = .TRUE.
        L_T_incr_bl_lsc   = .TRUE.
        L_Tl_incr_bl_lsc  = .TRUE.
        L_q_incr_bl_lsc   = .TRUE.
        L_qtl_incr_bl_lsc = .TRUE.
        L_qcl_incr_bl_lsc = .TRUE.
        L_qcf_incr_bl_lsc = .TRUE.
        L_cf_incr_bl      = .TRUE.
        L_cfl_incr_bl     = .TRUE.
        L_cff_incr_bl     = .TRUE.

! SCM initialise arrays for diagnostics
      If ( Substep_Number == 1 ) Then
          Do k=1,model_levels
            Do j=1, rows
            Do i=1, row_length
              u_incr_diag_bl(i,j,k) = 0.0
              v_incr_diag_bl(i,j,k) = 0.0
             T_incr_diag_bl(i,j,k) = 0.0
            Enddo
            Enddo
          End Do ! k

          Do k=1,wet_levels
            Do j=1, rows
            Do i=1, row_length
              q_incr_diag_bl(i,j,k)   = 0.0
              qcl_incr_diag_bl(i,j,k) = 0.0
              qcf_incr_diag_bl(i,j,k) = 0.0
              bulk_cf_incr_diag_bl(i,j,k)   = 0.0
              cf_liquid_incr_diag_bl(i,j,k) = 0.0
              cf_frozen_incr_diag_bl(i,j,k) = 0.0
            Enddo
            Enddo
          End Do ! k

      End If ! If ( Substep_Number == 1 )
#endif

! DEPENDS ON: ni_imp_ctl
      Call NI_imp_ctl (                                                 &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &

! IN Substeping information
     &, Substep_Number, Num_Substeps, L_phys2_substep                   &

! model dimensions.
     &, row_length, rows, rhc_row_length, rhc_rows, n_rows, land_points &
     &, ntiles, model_levels, wet_levels, bl_levels, dst_levels         &
     &, dsm_levels, cloud_levels, n_cca_levels, tr_levels, nice         &
     &, DIM_CS1, DIM_CS2                                                &
!
! Model switches
     &, model_domain, L_CAL360, L_area_cloud, L_ACF_Cusack              &
     &, L_ACF_Brooks, L_RHCPT, L_emcorr, Ltimer, L_DRY,L_MURK           &
     &, L_MURK_ADVECT,L_BL_TRACER_MIX,L_DUST,L_CAM_DUST,L_SULPC_SO2     &
     &, L_sulpc_nh3, L_sulpc_dms, L_soot, L_biomass, L_ocff             &
     &, L_co2_interactive                                               &
! rml 1/7/13
     &, L_CO2_TRACER                                                    &
     &, L_co2_emits, L_us_blsol, L_eacf, L_pc2, L_use_bl_diag_term      &
     &, NumCycles, CycleNo, L_mixing_ratio                              &
     &, L_ukca, L_sice_heatflux                                         &
     &, L_subfilter_vert, L_subfilter_horiz, L_USE_CARIOLLE             &

! model Parameters
     &, alpha_cd, Puns, Pstb, trweights1, rhcrit, tr_vars               &
     &, cloud_fraction_method, overlap_ice_liquid,ice_fraction_method   &
     &, ctt_weight,t_weight,qsat_fixed,sub_cld,x1i,x1ic,x1r,x2r,x4r     &
     &, l_psd,ai,bi,aic,bic,lsp_ei,lsp_fi,lsp_eic,lsp_fic               &

! Physical constants
     &, lc,lf,cp,two_Omega,p_zero,kappa,R,g,Lapse_Rate,earth_radius,Pi  &

! in coordinate information
     &, r_rho_levels, r_theta_levels, r_at_u, r_at_v, eta_theta_levels  &
     &, eta_rho_levels, delta_lambda, delta_phi,lat_rot_NP,long_rot_NP  &

! in time stepping information.
     &, sub_timestep, val_year, val_day_number, val_hour, val_minute    &
     &, val_second, timestep_number                                     &

! trig arrays
     &, sin_theta_longitude, cos_theta_longitude, FV_cos_theta_latitude &

! diagnostic info
     &     ,                                                            &
#include "argsts.h"
     & STASHwork3, STASHwork9                                           &
!
! SCM Diagnostics (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! in data fields.
     &, p, p_layer_centres, p_layer_boundaries, rho_rsq, u, v, w        &
     &, land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels   &
     &, u_conv,v_conv,theta_conv,q_conv,qcl_conv,qcf_conv               &
! ancillary fields and fields needed to be kept from timestep to
! timestep
     &, smvccl, smvcwt, smvcst, sthf, sthu,sil_orog_land, ho2r2_orog    &
     &, di, ice_fract, di_ncat, ice_fract_ncat                          &
     &, u_0, v_0, land_index, cca, ccb, cct, ccw, surf_radflux          &
     &, ls_rain, ls_snow, conv_rain, conv_snow, cca_2d, L_scrn, L_plsp  &

! in variables required from BDY_LAYR
     &, alpha1_sice, ashtf, dtrdz_charney_grid, rdz_charney_grid        &
     &, dtrdz_u, dtrdz_v, rdz_u, rdz_v, cdr10m_u, cdr10m_v, z1_tq       &
!ajm   extra variable
     &, rhokm_u, rhokm_v                                                &

! in diagnostics (started or from) BDY_LAYR
! pass diagn arrays and reduce number of argument lines
     &, T_incr_diag_bl, q_incr_diag_bl, qcl_incr_diag_bl                &
     &, qcf_incr_diag_bl, u_incr_diag_bl, v_incr_diag_bl                &
     &, cf_liquid_incr_diag_bl, cf_frozen_incr_diag_bl                  &
     &, bulk_cf_incr_diag_bl                                            &
     &, CD,CH,E_SEA,FQT,FTL,H_SEA,RIB_GB,TAUX,TAUY,VSHR,ZHT             &
     &, bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6     &
     &, bl_type_7, z0m_gb, z0m_eff_gb, z0h_eff_gb, fme, rhokh           &

! in data required to calculate increment diagnostics
     &, theta_star,q_star                                               &

! IN logical for scm surface forcing
     &, L_flux_bc                                                       &

! in data required for tracer mixing :
     &, RHO_ARESIST,ARESIST,RESIST_B,R_B_DUST                           &
     &, KENT, WE_LIM, T_FRAC, ZRZI                                      &
     &, KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC                      &
     &, ZHSC,Z_HALF,DUST_FLUX,U_S_T_TILE,U_S_T_DRY_TILE,U_S_STD_TILE    &
     &, so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em          &
     &, ocff_hilem, ocff_em                                             &

! in data required for 4D-VAR :
     &, rho_cd_modv1, rho_km_var                                        &

! IN additional variables for MOSES II. Now includes lai_ft, canht_ft.
     &, TILE_PTS,TILE_INDEX,TILE_FRAC,CANOPY,L_NEG_TSTAR                &
     &, ALPHA1,FRACA,RHOKH_TILE,SMC,CHR1P5M,RESFS,Z0HSSI,Z0MSSI         &
     &, CANHC_TILE,FLAKE,WT_EXT_TILE,LW_DOWN,lai_ft,canht_ft            &
     &, SW_TILE,ASHTF_TILE,gc,aresist_tile,resist_b_tile                &
     &, FQT_ICE,FTL_ICE,RESFT,RHOKH_SICE,RHOKPM,RHOKPM_POT              &
     &, RHOKPM_SICE,Z0H_TILE,Z0M_TILE,CHR1P5M_SICE                      &
     &, FLAND,FLANDG,FLANDG_U,FLANDG_V,TSTAR_SEA,VSHR_LAND,VSHR_SSI     &

!    EAK
!    IN
     &, surf_down_sw,alb_tile,cos_zenith_angle,l_tile_pts               &
!     &, surf_down_sw,alb_tile,cos_zenith_angle               &
     &, lat,long,day,time_sec                                           &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L                   &
     &, SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE             &
     &, TSOIL_TILE,T_SURF_TILE,HCONS,clapp                              &
     &, SATHH,SATCON,HCAP,HCON,Z1_UV                                    &
     &, SOIL_LAYER_MOISTURE                                             &
     &, SOIL_TYPE,VEG_TYPE                                              &
     &, ISNOW_FLG3L,total_nsteps                                        &
     &, FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB                         &
     &, TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB                     &
     &, USTAR_CAB,SURF_HTF_CAB                                          &
     &, l_cable                                                         &

!sxy
     &, TOT_ALB                                                         &
     &, SNAGE_TILE,RTSOIL_TILE                                          &
     &, GFLUX_TILE,SGFLUX_TILE                                          &
     &, T1P5M,Q1P5M,CANOPY_GB                                           &
     &, TRANSP_TILE                                                     &
     &, CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER                     &
     &, GLAI,PHENPHASE                                                  &
     ! Lestevens 23apr13
     &, NPP_FT_ACC,RESP_W_FT_ACC                                        &
 
! IN MOSES II variables for STASH
! kdcorbin, 10/10 - added resp_s_tile
     &, GS,GPP,NPP,RESP_P,GPP_FT,NPP_FT,RESP_P_FT,RESP_S                &
     &, RESP_S_TOT,RESP_S_TILE,CS                                       &
     &, RIB_TILE,FSMC,CATCH,G_LEAF                                      &
     &, CO2_EMITS, CO2FLUX                                              &

! IN tracer fluxes - kdcorbin, 05/10
     &, TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4          &
     &, TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8          &
     &, TRACER_FLUX9, TRACER_FLUX10,TRACER_FLUX11,TRACER_FLUX12         &
     &, TRACER_FLUX13,TRACER_FLUX14,TRACER_FLUX15,TRACER_FLUX16         &
     &, TRACER_FLUX17,TRACER_FLUX18,TRACER_FLUX19,TRACER_FLUX20         &

! IN co2 global emissions - kdcorbin, 05/10
     &, co2emitmass                                                     &         

!IN additional variables for soil moisture nudging scheme
     &, WT_EXT,RA                                                       &

!IN additional variables for anthropogenic heat
     &, ANTHROP_HEAT                                                    &
     
! in/out
     &, t_soil, ti, t_surf, ti_gb                                       &
     &, area_cloud_fraction, bulk_cloud_fraction                        &
     &, T_latest, q_latest, qcl_latest, qcf_latest                      &
     &, cf_latest, cfl_latest, cff_latest                               &
     &, R_u, R_v, R_w, cloud_fraction_liquid, cloud_fraction_frozen     &
     &, zh, sum_eng_fluxes,sum_moist_flux, rhcpt                        &

! In/Out tracer fields
     &, aerosol, free_tracers                                           &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, drydep2, so2, dms, so4_aitken, so4_accu, so4_diss, nh3          &
     &, soot_new, soot_aged, soot_cld, bmass_new, bmass_aged            &
     &, bmass_cld, ocff_new, ocff_aged, ocff_cld, co2, OZONE_TRACER     &

! in/out fields
     &, ecan, ei, ext, snowmelt                                         &
     &, t1_sd, q1_sd, ntml, cumulus                                     &
     &, l_pc2_diag_sh_pts, nbdsc, ntdsc                                 &
     &, surf_ht_flux_land, snomlt_surf_htf, cH_term                     &

! INOUT additional variables for MOSES II
     &, TSTAR_TILE,FQT_TILE,EPOT_TILE,FTL_TILE                          &
     &, SNOW_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,OLR                   &
     &, TSTAR_LAND,TSTAR_SICE,TSTAR_SSI                                 &

! Additional diagnostics for MOSES II
     &, rib_ssi,taux_land,taux_ssi,tauy_land,tauy_ssi                   &

! OUT additional variables for MOSES II
     &, ESOIL_TILE,ES,EI_TILE                                           &
     &, Q1P5M_TILE,T1P5M_TILE,ECAN_TILE,MELT_TILE                       &

! error information
     &, Error_code,BL_diag  )




! DEPENDS ON: timer
      If (Ltimer) Call timer ('NI_IMP_CTL',4)

      Else

          Do k = 1, wet_levels
            RHCPT(1,1,k) = rhcrit(k)
          End Do
      End If

      End If ! error_code == 0

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Boundary Layer',6)

! ----------------------------------------------------------------------
! Section RESTORE: Copy _latest fields into _star locations
! Area cloud fraction has been done in NI_IMP_CTL
! ----------------------------------------------------------------------
      If (error_code  ==  0) Then
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              theta_star(i,j,k) = T_latest(i,j,k) /                     &
     &                            exner_theta_levels(i,j,k)
            End Do
          End Do
        End Do
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              q_star(i,j,k) = q_latest(i,j,k)
              qcl_star(i,j,k) = qcl_latest(i,j,k)
              qcf_star(i,j,k) = qcf_latest(i,j,k)
              cf_star(i,j,k)  = cf_latest(i,j,k)
              cfl_star(i,j,k) = cfl_latest(i,j,k)
              cff_star(i,j,k) = cff_latest(i,j,k)
            End Do
          End Do
        End Do
! Reset halo variables for writing back into D1 array on exit from
! this subroutine.
        If (.not. L_pc2 .and. (Substep_number  ==  Num_Substeps) ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                  bulk_cloud_fraction_halos(i,j,k) =                    &
     &              bulk_cloud_fraction(i,j,k)
                  cloud_fraction_liquid_halos(i,j,k) =                  &
     &              cloud_fraction_liquid(i,j,k)
                  cloud_fraction_frozen_halos(i,j,k) =                  &
     &              cloud_fraction_frozen(i,j,k)
              End Do
            End Do
          End Do
        End If  ! .not. L_pc2

      End If ! on error code equal to zero

#if !defined(SCMA)
      If ( L_bl ) Then
        If ( Substep_Number == Num_Substeps ) Then
          Deallocate ( u_incr_diag_bl )
          Deallocate ( v_incr_diag_bl )
          Deallocate ( T_incr_diag_bl )
          Deallocate ( q_incr_diag_bl )
          Deallocate ( qcl_incr_diag_bl )
          Deallocate ( qcf_incr_diag_bl )
          Deallocate ( bulk_cf_incr_diag_bl )
          Deallocate ( cf_liquid_incr_diag_bl )
          Deallocate ( cf_frozen_incr_diag_bl )
        End If
      End If
#endif

! ----------------------------------------------------------------------
! Section HYD.1 Compress fields to land points, then call hydrology.
! ----------------------------------------------------------------------

! Call hydrology only at last cycle
      If ( CycleNo == NumCycles) Then
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Hydrology',5)

      If (error_code  ==  0 .and. L_hydrology                           &
     &    .and. land_points  /=  0) Then

! Set inland basin outflow to zero unless river routing
! called in prev timestep

      IF(L_RIVERS)THEN
! If inland river routing basin required
       IF(L_INLAND)THEN
        IF (A_STEPS_SINCE_RIV /= 0) THEN
! If river routing called in previous timestep
          inlandout_atmos=0.0
          inlandout_atm=0.0
          INLANDOUT_RIV=0.0

        ENDIF

        ELSE
! If inland river routing basin not required
          inlandout_atmos=0.0
          INLANDOUT_RIV=0.0
       ENDIF
      ENDIF


! Compress fields to land points

        Do l = 1, land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          ls_rain_land(l)=ls_rain(i,j)
          conv_rain_land(l)=conv_rain(i,j)
          conv_snow_land(l)=conv_snow(i,j)
          ls_snow_land(l)=ls_snow(i,j)
          surf_ht_flux_ld(l) = surf_ht_flux_land(i,j)
          LYING_SNOW(L) = snow_depth(i,j)
          snow_melt(l) = snowmelt(i,j)
        End Do


#if defined(SCMA)
!-----------------------------------------------------------------------
!     SCM Land Points Diagnostics Package
!-----------------------------------------------------------------------
        If (L_SCMDiags(SCMDiag_land)) Then

          ! snomlt_surf_htf is overwritten by HYD_INTCTL, so need
          ! to output here
! DEPENDS ON: scmoutput
          Call SCMoutput(snomlt_surf_htf,'snomlt_surf_htf',             &
               'Snowmelt heatflux (boundary layer scheme)','W/m2',      &
               t_avg,d_sl,default_streams,'',RoutineName)

        End If ! L_SCMDiags(SCMDiag_land)
#endif

! DEPENDS ON: hyd_intctl
        Call HYD_INTCTL (                                               &
     &               land_ice_points, land_ice_index,                   &
     &               soil_points, soil_index,                           &
     &               land_points, dsm_levels, clapp,                    &
     &               conv_rain_land, conv_snow_land,                    &
     &               ext, hcap, hcon,                                   &
     &               ls_rain_land, ls_snow_land,                        &
     &               satcon, sathh,                                     &
     &               surf_ht_flux_ld, sub_timestep,                     &
     &               smvcst, smvcwt,                                    &
     &               canopy_gb, snomlt_surf_htf, smlt,                  &
     &               soil_layer_moisture,                               &
     &               sthf, sthu, t_soil, lying_snow,                    &
! output
! Beware snow_melt is actually in/out
     &               infil, smc, snow_melt, snomlt_sub_htf,             &
     &               stf_sub_surf_roff, sub_surf_roff, surf_roff,       &
     &               tot_tfall,                                         &
! add inland basin outflow in call to hydctl
     &               inlandout_atm,L_INLAND,                            &


! Additional variables for MOSES II
     &               NTILES,TILE_PTS,TILE_INDEX,                        &
     &               L_SNOW_ALBEDO,CAN_MODEL,                           &
     &               CATCH,ECAN_TILE,INFIL_TILE,                        &
     &               MELT_TILE,TSTAR_TILE,TILE_FRAC,CANOPY,             &
     &               RGRAIN,SNOW_TILE,                                  &
     &               CATCH_SNOW,SNOW_GRND,                              &

! Additional variables required for large-scale hydrology:
     &               L_TOP,L_PDM,FEXP,GAMTOT,TI_MEAN,TI_SIG,CS,         &
     &               DUN_ROFF,FSAT,FWETL,QBASE,QBASE_ZW,ZW,DRAIN,       &
     &               STHZW,A_FSAT,C_FSAT,A_FWET,C_FWET,FCH4_WETL,       &
     &               L_SOIL_SAT_DOWN,                                   &
! sxy
       L_TILE_PTS,ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS,                     &
!       ROW_LENGTH,ROWS,N_ROWS,BL_LEVELS,                     &
       LAND_INDEX,SMVCCL,SMVCST,SMVCWT,T_SURF,                          &
       FQT,FTL,FTL_TILE,LE_TILE,RADNET_TILE,FQT_TILE,                   &
       SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,                   &
       SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,             &
       TSOIL_TILE,T_SURF_TILE,HCONS,                                    &
       SOIL_TYPE,VEG_TYPE,ISNOW_FLG3L,total_nsteps,                     &
       FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,                         &
       TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,                     &
       SURF_HTF_CAB,SURF_HT_FLUX_LAND,                                  &
       ESOIL_TILE,EI_TILE,                                              &
       SNAGE_TILE,GFLUX_TILE,SGFLUX_TILE,                               &
! EAK
     &   l_cable, WB_LAKE,                                              &

! Timer diagnostics disabled
     &               Ltimer )


! Copy land points output back to full fields array.
        Do l = 1, land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          snow_depth(i,j) = LYING_SNOW(L)
        End Do

      End If

  WBLAKE  = 0.0
  SUBROFF = 0.0

if (L_CABLE) then
! Lestevens - moved for lake fix
! Calculate the Atmosphere gridbox areas
          DO J = 1, rows
            DO I = 1, row_length
              A_BOXAREAS(I,J) = r_theta_levels(i,j,0)                   &
     &                      * r_theta_levels(i,j,0)                     &
     &                      * delta_lambda * delta_phi                  &
     &                      * FV_cos_theta_latitude(i,j)
            ENDDO
          ENDDO

   !WBLAKE  = 0.0
   !SUBROFF = 0.0

       DO N=1,NTILES
       DO K=1,TILE_PTS(N)
         L = TILE_INDEX(K,N)
         J = (LAND_INDEX(L)-1)/ROW_LENGTH + 1
         I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
         WBLAKE(L)  = WBLAKE(L)  + FLAND(L)*TILE_FRAC(L,N)           &
                                 *WB_LAKE(L,N)*A_BOXAREAS(I,J)
         SUBROFF(L) = FLAND(L)*TIMESTEP*sub_surf_roff(L)*A_BOXAREAS(I,J)
         !SUBROFF0(L) = FLAND(L)*TIMESTEP*sub_surf_roff(L)
         !SUBROFF1(L) = TIMESTEP*sub_surf_roff(L)
       ENDDO
       ENDDO

END IF
!EAK
!       TOTWBLAKE =  sum ( WBLAKE )
!       TOTWBX =  sum ( SUBROFF )
!       TOTWBX0 =  sum ( SUBROFF0 )
!       TOTWBX1 =  sum ( SUBROFF1 )
!        print 88,TOTWBX,TOTWBX0,TOTWBX1,  &
!         TOTWBLAKE
! 88     format(1x,'atm,TOT2',10e11.4)


! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Hydrology',6)

!-------------------------------------------------------------------
! RIVER ROUTING

      IF ( L_RIVERS .AND. SUBSTEP_NUMBER == NUM_SUBSTEPS ) THEN

! Set the river routing to run on the 'last' PE as PE0 is very busy
        gather_pe_trip = n_proc-1

! Initialise diagnostics on non-river routing timesteps
          RIVEROUT = 0.0
          BOX_OUTFLOW = 0.0
          BOX_INFLOW = 0.0
! Initialise the gather fields at the beginning of TRIP timestep

        IF(A_STEPS_SINCE_RIV == 0)THEN
          TOT_SURF_RUNOFF=0.0
          TOT_SUB_RUNOFF=0.0
          TOT_WBLAKE  = 0.0
          TOT_SUBRUN = 0.0
        ENDIF

        A_STEPS_SINCE_RIV = A_STEPS_SINCE_RIV + 1
        NSTEP_TRIP=INT(RIVER_STEP/TIMESTEP)
        IF (A_STEPS_SINCE_RIV == NSTEP_TRIP) THEN
          TRIP_CALL=.TRUE.
        ELSE
          TRIP_CALL=.FALSE.
        ENDIF

! Accumulate the runoff as Kg/m2/s over the TRIP period
        DO l = 1, land_points
          IF(surf_roff(L) <  0.0)THEN
            WRITE(6,*)'surf_roff(',L,')= ',surf_roff(L)
          ELSE
            TOT_SURF_RUNOFF(L) = TOT_SURF_RUNOFF(L) +                   &
     &                   (surf_roff(L)/ REAL(NSTEP_TRIP))
          ENDIF
        ENDDO
        DO l = 1, land_points
          IF(sub_surf_roff(L) <  0.0)THEN
            WRITE(6,*)'sub_surf_roff(',L,')= ',sub_surf_roff(L)
          ELSE
            TOT_SUB_RUNOFF(L) = TOT_SUB_RUNOFF(L) +                     &
     &                   (sub_surf_roff(L)/ REAL(NSTEP_TRIP))
          ENDIF
        ENDDO
        IF (L_CABLE) THEN
        DO l = 1, land_points
          IF(wblake(L) <  0.0)THEN
            WRITE(6,*)'wblake(',L,')= ',wblake(L)
          ELSE
            TOT_WBLAKE(L) = TOT_WBLAKE(L) +                     &
     &                   (WBLAKE(L)/ REAL(NSTEP_TRIP))
          ENDIF
        ENDDO

        DO l = 1, land_points
          IF(subroff(L) <  0.0)THEN
            WRITE(6,*)'subroff(',L,')= ',subroff(L)
          ELSE
            TOT_SUBRUN(L) = TOT_SUBRUN(L) +                   &
     &                   (SUBROFF(L)/ REAL(NSTEP_TRIP))
          ENDIF
        ENDDO
        ENDIF

! detect first entry into river routing
      FIRST_ROUTING = .FALSE.
      IF (TIMESTEP_NUMBER == NSTEP_TRIP) FIRST_ROUTING = .TRUE.

        IF(TRIP_CALL)THEN

! DEPENDS ON: timer
        If (Ltimer) Call timer ('AP2 River Routing',5)

! If ATMOS fields are as Ocean (i.e. inverted NS) set INVERT_ATMOS
           INVERT_ATMOS = .FALSE.
           IF(.NOT.INVERT_OCEAN)INVERT_ATMOS = .TRUE.

if (.not.l_cable) then
! Calculate the Atmosphere gridbox areas
          DO J = 1, rows
            DO I = 1, row_length
              A_BOXAREAS(I,J) = r_theta_levels(i,j,0)                   &
     &                      * r_theta_levels(i,j,0)                     &
     &                      * delta_lambda * delta_phi                  &
     &                      * FV_cos_theta_latitude(i,j)
            ENDDO
          ENDDO
end if
! DEPENDS ON: riv_intctl
          CALL RIV_INTCTL(                                              &
     & XPA, XUA, XVA, YPA, YUA, YVA,                                    &
     & G_P_FIELD, G_R_FIELD, N_PROC, ME, RMDI,                          &
     & GATHER_PE_TRIP,LAND_POINTS,LAND_INDEX,                           &
     & INVERT_ATMOS, ROW_LENGTH, ROWS,                                  &
     & GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                                  &
     & RIVER_ROW_LENGTH, RIVER_ROWS,                                    &
     & GLOBAL_RIVER_ROW_LENGTH, GLOBAL_RIVER_ROWS,                      &
     & FLANDG, RIVER_STEP, RIVER_VEL, RIVER_MCOEF,                      &
     & TRIVDIR, TRIVSEQ, TWATSTOR, A_BOXAREAS,                          &
     &  DELTA_PHI,EARTH_RADIUS,FIRST_ROUTING,                           &
     &  r_area, slope, flowobs1,r_inext,r_jnext,r_land,                 &
     &  substore,surfstore,flowin,bflowin,                              &
! IN/OUT accumulated runoff
     & TOT_SURF_RUNOFF, TOT_SUB_RUNOFF,                                 &
! OUT
     &  BOX_OUTFLOW, BOX_INFLOW, RIVEROUT                               &
! Add inland basin arguments in call to rivctl
     &  ,INLANDOUT_ATMOS,INLANDOUT_RIV                                  &
     &  ,TOT_WBLAKE,TOT_SUBRUN,timestep,l_cable                        &
     &  ,wblake_ratio )

! compress inland basin outputs to land points only
          if (L_INLAND) then
            do l = 1,land_points
              j = (land_index(l)-1)/row_length +1
              i = land_index(l) - (j-1)*row_length
              inlandout_atm(l) = inlandout_atmos(i,j)
            end do
          end if

! Mult RIVEROUT by the number of physics timesteps per River routing
! timestep as DAGHYD stores RIVEROUT every timestep. Non-routing
! timestep vals are passed in as 0.0

          A_STEPS_SINCE_RIV = 0

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 River Routing',6)


! -------------------------------------------------------------------
! Section RIVER Output diagnostics
! ------------------------------------------------------------------

! Synchronize cpu's to avoid shmem time-outs due to load imbalance
      Call gc_ssync(n_proc, i)

! Check that river diagnostics requested this timestep
#if !defined(SCMA)
      If (error_code  ==  0  .and. sf(0,26) ) Then

! DEPENDS ON: timer
        If (Ltimer) Call timer ('diagnostics_riv',3)
! DEPENDS ON: diagnostics_riv
        Call diagnostics_riv(                                           &
     &                       row_length, rows                           &
     &,                      river_row_length, river_rows               &
     &,                      at_extremity                               &
     &,                      at_extremity                               &
     &,                      RIVEROUT                                   &
     &,                      BOX_OUTFLOW, BOX_INFLOW                    &
! Put inland basin outflow in call to dagriv
     &,                      TWATSTOR,INLANDOUT_RIV                     &

     &,                                                                 &
#include "argsts.h"
     & STASHwork26                                                      &
     & )


! DEPENDS ON: timer
        If (Ltimer) Call timer ('diagnostics_riv',4)

      End If ! on error code .and. L_rivers .and. sf(0,26)
#endif
        ENDIF                        !TRIP_CALL

      ENDIF                        ! L_RIVERS (ATMOS)

! ------------------------------------------------------------------

! Synchronize cpu's to avoid shmem time-outs due to load imbalance  ???
      Call gc_ssync(n_proc, i)

! ----------------------------------------------------------------------
! Section HYD.2 Output diagnostics
! ----------------------------------------------------------------------

!Lestevens 30.10.13 - fix sub-runoff for diagnostics
if (l_cable) then
  DO L = 1,LAND_POINTS
    SUB_SURF_ROFF(L) = max(0.,(1.-wblake_ratio)*(SUB_SURF_ROFF(L)))
  ENDDO
endif


! Check that hydrology diagnostics requested this timestep
#if !defined(SCMA)
      If (error_code  ==  0 .and. L_hydrology .and. sf(0,8)             &
     &     .and. Substep_Number  ==  Num_Substeps ) Then

! DEPENDS ON: timer
        If (Ltimer) Call timer ('diagnostics_hyd',3)
! DEPENDS ON: diagnostics_hyd
        Call diagnostics_hyd(                                           &
     &                       row_length, rows, model_levels             &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      land_points, dsm_levels                    &
 ! Put inland basin outflow in call to dagriv
     &,                      land_index,inlandout_atm                   &

     &,                      smc, surf_roff, sub_surf_roff              &
     &,                      lying_snow, snow_melt                      &
     &,                      canopy_gb,t_soil                           &
     &,                      soil_layer_moisture                        &
     &,                      ntiles, snomlt_surf_htf, sthu, sthf        &
     &,                      tot_tfall, snow_tile, melt_tile            &
     &,                      rgrain, land_sea_mask                      &
     &,                      dun_roff, drain, qbase, qbase_zw           &
     &,                      fch4_wetl                                  &
     &,                      fexp,gamtot,ti_mean,ti_sig                 &
     &,                      fsat,fwetl,zw,sthzw                        &
     &,                      timestep                                   &
     &,                                                                 &
#include "argsts.h"
     & STASHwork8                                                       &
     & )

! DEPENDS ON: timer
        If (Ltimer) Call timer ('diagnostics_hyd',4)

      End If ! on error code .and. L_hydrology .and. sf(0,8)
#endif

! ----------------------------------------------------------------------
! Section 19 -- VEGETATION DYNAMICS
! ----------------------------------------------------------------------

      IF ( SUBSTEP_NUMBER  ==  NUM_SUBSTEPS ) THEN
!-----------------------------------------------------------------------
! If leaf phenology is activated, check whether the atmosphere model
! has run an integer number of phenology calling periods.
!-----------------------------------------------------------------------
      PHENOL_CALL=1
      TRIFFID_CALL=1
      IF (L_PHENOL) THEN
        PHENOL_CALL = MOD ( FLOAT(A_STEP),(FLOAT(PHENOL_PERIOD)*        &
     &  (86400.0/TIMESTEP)) )
      ENDIF

      IF (L_TRIFFID) THEN


        NSTEP_TRIF=INT(86400.0*TRIFFID_PERIOD/TIMESTEP)
        IF (ASTEPS_SINCE_TRIFFID == NSTEP_TRIF) THEN
          TRIFFID_CALL=0
        ENDIF
      ENDIF

      IF ((PHENOL_CALL == 0).OR.(TRIFFID_CALL == 0)) THEN
! DEPENDS ON: veg_ctl
        CALL VEG_CTL(                                                   &
     &               row_length, rows, n_rows                           &
     &,              global_row_length, global_rows                     &
     &,              DIM_CS1, DIM_CS2                                   &
     &,              halo_i, halo_j, off_x, off_y, me                   &
     &,              n_proc, n_procx, n_procy                           &
     &,              g_rows, g_row_length, g_datastart                  &
     &,              at_extremity                                       &
     &,              LAND_POINTS,LAND_INDEX,NTILES,CAN_MODEL            &
     &,              A_STEP,ASTEPS_SINCE_TRIFFID                        &
     &,              PHENOL_PERIOD,TRIFFID_PERIOD                       &
     &,              L_PHENOL,L_TRIFFID,L_TRIF_EQ                       &
     &,              TIMESTEP,FRAC_DISTURB,SATCON                       &
     &,              G_LEAF_ACC,G_LEAF_PHEN_ACC,NPP_FT_ACC              &
     &,              RESP_S_ACC,RESP_W_FT_ACC                           &
     &,              CS,FRAC,LAI_FT,SOIL_CLAY,CANHT_FT                  &
     &,              CATCH_SNOW,CATCH,INFIL_TILE,Z0_TILE                &

     &,                                                                 &
#include "argsts.h"
     &               STASHwork19                                        &
     &                   )
      ENDIF

      ENDIF

! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics for output to D1 array:
! ----------------------------------------------------------------------
      IF(L_CTILE)THEN
        DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            TSTAR_LAND_CTILE(I,J)=TSTAR_LAND(I,J)
            TSTAR_SICE_CTILE(I,J)=TSTAR_SICE(I,J)
          ENDDO
        ENDDO
      ENDIF
      End If ! If ( CycleNo == NumCycles ) Then
      Enddo ! terminate substepping loop

#if defined(SCMA)
      ! The SCM output diagnostic system needs to know that
      ! sub-stepping has ended

! DEPENDS ON: scm_substepping_end
      Call SCM_substepping_end

! ----------------------------------------------------------------------
! Output some SCM diagnostics
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_bl)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(rib_gb,                                          &
             'rib','Richardson number (lowest layer)','-',              &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(taux,                                            &
             'taux','Wind stress x','N/m2',                             &
             t_avg,d_bl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(tauy,                                            &
             'tauy','Wind stress y','M/m2',                             &
             t_avg,d_bl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(fqT,                                             &
             'fqT','Turbulent moisture flux','kg/m2/s',                 &
             t_avg,d_bl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(ftl,                                             &
             'ftl','Turbulent sensible heat flux','W/m2',               &
             t_avg,d_bl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_bl)
!
!-----------------------------------------------------------------------
!     SCM Surface Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_surf)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(ei,                                              &
             'sublim','Sublimation from lying snow or sea ice',         &
             'kg/m2/day',                                               &
             t_mult,d_sl,default_streams,'sec_day',RoutineName)

      End If ! L_SCMDiags(SCMDiag_surf)
!
!-----------------------------------------------------------------------
!     SCM Land Points Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_land)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(ftl_tile,                                        &
             'shf_tile','Tile sensible heat flux','W/m2',               &
             t_avg,d_tile,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(le_tile,                                         &
             'lhf_tile','Tile latent heat flux','W/m2',                 &
             t_avg,d_tile,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(t1p5m_tile,                                      &
             't1p5m_tile','Tile 1.5m temperature','K',                  &
             t_avg,d_tile,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(q1p5m_tile,                                      &
     &       'q1p5m_tile','Tile 1.5m specific humidity','kg/kg',         &
     &       t_avg,d_tile,default_streams,'',RoutineName)    

! DEPENDS ON: scmoutput
        Call SCMoutput(surf_ht_flux_land,'surf_ht_flux_ld',             &
             'Net downward heat flux into land fraction','W/m2',        &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(ecan,                                            &
             'can_evap','Canopy evaporation','kg/m2/day',               &
             t_mult,d_sl,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(es,                                              &
             'soil_evap','Soil evapotranspiration','kg/m2/day',         &
             t_mult,d_sl,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(snomlt_sub_htf,                                  &
             'snomlt_sub_htf','Snow melt heat flux into sub-surface',   &
             'W/m2',                                                    &
             t_avg,d_land,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(snow_melt,                                       &
             'snomlt_hyd','Snow melt (hydrology scheme)','kg/m2/day',   &
             t_mult,d_land,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(snowmelt,                                        &
             'snomlt_bl','Snow melt (boundary layer scheme)',           &
             'kg/m2/day',                                               &
             t_mult,d_land,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(tot_tfall,                                       &
             'thro_fall','Throughfall','kg/m2/dat',                     &
             t_mult,d_land,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(surf_roff,                                       &
             'surf_roff','Surface runoff','kg/m2/day',                  &
             t_mult,d_land,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(sub_surf_roff,                                   &
             'Sub_roff','Sub-surface runoff','kg/m2/day',               &
             t_mult,d_land,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(gpp,                                             &
             'gpp','Gross primary productivity','kg C/m2/s',            &
             t_mult,d_land,default_streams,'oneKsecday',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(npp,                                             &
             'npp','Net primary productivity','kg C/m2/s',              &
             t_mult,d_land,default_streams,'oneKsecday',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(resp_p,                                          &
             'resp_p','Plant respiration','kg/m2/s',                    &
             t_mult,d_land,default_streams,'oneKsecday',RoutineName)

      End If ! L_SCMDiags(SCMDiag_land)
!
!-----------------------------------------------------------------------
!     SCM General OR Surface Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_gen)                                       &
     &    .OR. L_SCMDiags(SCMDiag_surf)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(ftl(1,1,1),                                      &
             'sens_ht','Surface sensible heat flux','W/m2',             &
             t_avg,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_surf)
!
!-----------------------------------------------------------------------
!     SCM Convection OR Increments Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_conv)                                      &
     &    .OR. L_SCMDiags(SCMDiag_incs)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(theta_inc,                                       &
             'dth_conv','Convective increment - theta','K',             &
             t_mult,d_all,default_streams,'sec_day',RoutineName)

      End If ! L_SCMDiags(SCMDiag_conv)
!
!-----------------------------------------------------------------------
!     SCM General OR Large Scale Precipitation Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_lsp)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(ls_rain,                                         &
             'ls_rain','Large scale rainfall rate','kg/m2/day',         &
             t_mult,d_sl,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(ls_snow,                                         &
             'ls_snow','Large scale snowfall rate','kg/m2/day',         &
             t_mult,d_sl,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(ls_rain,                                         &
             'ls_rain_inst','Large scale rainfall rate','kg/m2/s',      &
             t_inst,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(ls_snow,                                         &
             'ls_snow_inst','Large scale snowfall rate','kg/m2/s',      &
             t_inst,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_lsp)
!
!-----------------------------------------------------------------------
!     SCM General OR Convection Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_conv)) Then

! DEPENDS ON: scmoutput
        Call SCMoutput(conv_rain,                                       &
             'conv_rain','Convective rainfall','kg/m2/day',             &
             t_mult,d_sl,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(conv_snow,                                       &
             'conv_snow','Convective snowfall','kg/m2/day',             &
             t_mult,d_sl,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(conv_rain,                                       &
             'conv_rain_inst','Convective rainfall','kg/m2/s',          &
             t_inst,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(conv_snow,                                       &
             'conv_snow_inst','Convective snowfall','kg/m2/s',          &
             t_inst,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_conv)
!
!-----------------------------------------------------------------------
!     SCM Large Scale Precipitation OR Convection OR General
!     Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_lsp) .OR. L_SCMDiags(SCMDiag_conv)         &
     &    .OR. L_SCMDiags(SCMDiag_gen)) Then

         ! Calculate precipitation rates
         Do j=1,rows
           Do i=1,row_length
             tot_rain(i,j)   = ls_rain(i,j)   + conv_rain(i,j)
             tot_snow(i,j)   = ls_snow(i,j)   + conv_snow(i,j)
             pptn_rate(i,j)  = ls_rain(i,j)   + ls_snow(i,j)            &
     &                       + conv_rain(i,j) + conv_snow(i,j)
             accum_pptn(i,j) = timestep * pptn_rate(i,j)
           End Do ! i
         End Do ! j

! DEPENDS ON: scmoutput
         Call SCMoutput(pptn_rate,                                      &
              'tot_precip_avg','Precipitation rate','kg/m2/s',          &
              t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
         Call SCMoutput(accum_pptn,                                     &
              'tot_precip_acc','Accumulated Precipitation','kg/m2/s',   &
              t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 5,216
! DEPENDS ON: scmoutput
        Call SCMoutput(pptn_rate,                                       &
             'tot_precip','Total precipitation rate','kg/m2/day',       &
             t_mult,d_sl,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(pptn_rate,                                       &
             'tot_precip_inst','Total precipitation rate','kg/m2/s',    &
             t_inst,d_sl,default_streams,'',RoutineName)

!       Stash 5,214
! DEPENDS ON: scmoutput
        Call SCMoutput(tot_rain,                                        &
             'tot_rain','Total rainfall rate','kg/m2/day',              &
             t_mult,d_sl,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(tot_rain,                                        &
             'tot_rain_inst','Total rainfall rate','kg/m2/s',           &
             t_inst,d_sl,default_streams,'',RoutineName)

!       Stash 5,215
! DEPENDS ON: scmoutput
        Call SCMoutput(tot_snow,                                        &
             'tot_snow','Total snowfall rate','kg/m2/day',              &
             t_mult,d_sl,default_streams,'sec_day',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(tot_snow,                                        &
             'tot_snow_inst','Total snowfall rate','kg/m2/s',           &
             t_inst,d_sl,default_streams,'',RoutineName)

      End If ! L_SCMDiags(SCMDiag_lsp) .OR. L_SCMDiags(SCMDiag_conv)
             ! .OR. L_SCMDiags(SCMDiag_gen)

#endif

      RETURN      ! end of routine Atmos_physics2
      END SUBROUTINE Atmos_Physics2
#endif
#endif
