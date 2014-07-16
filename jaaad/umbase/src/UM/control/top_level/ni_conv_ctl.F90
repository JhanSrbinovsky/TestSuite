#if defined(ATMOS)
#if defined(CONTROL) || defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine NI_conv_ctl

      Subroutine NI_conv_ctl (                                          &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, at_extremity, n_proc, n_procx, n_procy          &
     &, neighbour, g_rows, g_row_length, g_datastart, me                &

! Parameters for iterative SISL
     &, NumCycles, CycleNo                                              &
     &, Substep_Number, Num_Substeps                                    &

! Model dimensions.
     &, row_length, rows, n_rows                                        &
     &, rowsrowlength                                                   &
     &, model_levels, wet_model_levels, bl_levels, n_cca_levels         &
     &, tr_levels, tr_vars, tr_ukca                                     &

! Model switches
     &, model_domain, L_CAL360, L_calc_dxek, L_q_interact               &
     &, l_mixing_ratio, Ltimer                                          &
     &, L_MURK_SOURCE, L_DUST, L_SULPC_SO2, L_SULPC_DMS, L_SULPC_NH3    &
     &, L_soot, L_ocff, L_biomass, L_co2_interactive                    &
     &, L_USE_CARIOLLE                                                  &

! Model parameters
     &, dbsdtbs_turb_0                                                  &

! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &

! in coordinate information
     &, r_rho_levels, r_theta_levels, z_rho, z_theta                    &
     &, eta_theta_levels, eta_rho_levels                                &
     &, delta_lambda, delta_phi                                         &
     &, lat_rot_NP, long_rot_NP                                         &

! in time stepping information.
     &, timestep                                                        &
     &, val_year, val_day_number, val_hour, val_minute                  &
     &, val_second, timestep_number                                     &

! diagnostic info
     &     ,                                                            &
#include "argsts.h"
     & STASHwork                                                        &
     &, ls_rain, ls_snow                                                &
!
! SCM diagnostics (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! in data fields.
     &, rho,  rho_only, rho_theta                                       &
     &, u, v, p, p_star                                                 &
     &, land_sea_mask                                                   &
     &, p_layer_boundaries, p_layer_centres                             &
     &, exner_layer_boundaries, exner_layer_centres                     &
     &, t1_sd, q1_sd, exner_theta_levels                                &
     &, uw0_p, vw0_p, w_max, zlcl,  zlcl_uv,ztop                        &
     &, cumulus, l_shallow, l_congestus,l_congestus2,l_pc2_diag_sh_pts  &
     &,  ntml, ntpar, ntdsc, nbdsc                                      &
     &, wstar, wthvs, delthvu, ql_ad, ftl, fqt                          &
     &, shallowc, cu_over_orog, cape_undilute, cin_undilute             &

! IN for STPH_SCV
     &, f3_at_u ,FV_cos_theta_latitude,cos_theta_longitude              &
! in/out stash diagnostics
     &, up_flux, up_flux_half, dwn_flux, entrain_up, entrain_dwn        &
     &, detrain_up, detrain_dwn, uw_dp, vw_dp, uw_shall, vw_shall       &
     &, wqt_flux,wthetal_flux,wthetav_flux,wql_flux                     &
     &, wqt_cb,wthetal_cb,wqt_inv,wthetal_inv,sh_top,sh_base            &
     &, T_incr_diagnostic, q_incr_diagnostic                            &
     &, qcl_incr_diagnostic, qcf_incr_diagnostic                        &
     &, u_incr_diagnostic, v_incr_diagnostic                            &
     &, cf_liquid_incr_diagnostic, cf_frozen_incr_diagnostic            &
     &, bulk_cf_incr_diagnostic                                         &
     &, dubydt_pout, dvbydt_pout                                        &
     &, precip_deep, precip_shall, precip_mid, precip_cong, cape_out    &
     &, shallow_ind, congestus_ind, congestus_ind2, mid_ind             &
     &, ntml_diag, ntpar_diag, freeze_diag                              &
     &, kterm_deep, wstar_dn_diag, wstar_up_diag, mb1_diag, mb2_diag    &
     &, cg_term,cg_top,cg_base                                          &
     &, mf_deep,mf_congest,mf_shall,mf_midlev                           &
     &, dt_deep,dt_congest,dt_shall,dt_midlev                           &
     &, dq_deep,dq_congest,dq_shall,dq_midlev                           &
     &, du_deep,du_congest,du_shall,du_midlev                           &
     &, dv_deep,dv_congest,dv_shall,dv_midlev                           &

! in/out
     &, theta_n, q_n, qcl_n, qcf_n, cf_liquid_n, cf_frozen_n, bulk_cf_n &
     &, theta_inc, q_inc, qcl_inc, qcf_inc, cf_liquid_inc, cf_frozen_inc&
     &, bulk_cf_inc                                                     &
     &, R_U, R_V, AEROSOL                                               &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, CONSCAV_DUST, SO2, SO4_AITKEN, SO4_ACCU, SO4_DISS               &
     &, dms, nh3, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd    &
     &, bmass_cld, ocff_new, ocff_agd, ocff_cld, co2, conscav_so4ait    &
     &, conscav_so4acc, conscav_so4dis, conscav_agedsoot                &
     &, conscav_agedbmass, conscav_agedocff, free_tracers, ukca_tracers &
     &, OZONE_TRACER                                                    &

! out fields
     &, cca_rad, ccb_rad, cct_rad, cclwp_rad, conv_rain, conv_snow      &
     &, conv_rain_3d, conv_snow_3d                                      &
     &, ccw_rad, lcca, lcclwp, lctop_rad, lcbase_rad                    &

! error information
     &, Error_code  )

! purpose: Interface to Atmospheric Physics convection code.
!
! method:
!
! code description:
!   language: fortran 77 + cray extensions
!   this code is written to umdp3 programming standards.

      Use cv_cntl_mod, Only:                                            &
          lcv_3d_cca, lcv_ccrad, lcv_pc2_diag_sh

      Use cv_run_mod,  Only:                                            &
          l_convcld_hadgem1,       l_fix_udfactor,        l_mom,        &
          l_cape,                  l_scv,                 ud_factor,    &
          iconv_shallow,           iconv_congestus,       iconv_deep,   &
          convection_option,       a_convect_segments,    n_conv_calls, &
          rad_cloud_decay_opt,     cld_life_opt,          cca_min,      &
          fixed_cld_life,          a_convect_seg_size

      Implicit None

! arguments with intent in. ie: input variables.

! Parallel setup variables
      Integer , intent(in) ::                                           &
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
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, global_rows                                                     &
                           ! NUMBER OF global rows
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
     &, me                                                              &
                   ! My processor number

     &, NumCycles                                                       &
     &, CycleNo

! Substep number for ATMPHYS2
      Integer, intent(in) ::                                            &
     &  Substep_Number                                                  &
     &, Num_Substeps
      Logical, intent(in) ::                                            &
     &  at_extremity(4)  ! Indicates if this processor is at north, sout
                         ! east or west of the processor grid

! Model dimensions
      Integer, intent(in) ::                                            &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, rowsrowlength                                                   &
                        ! rows*row_length
     &, model_levels                                                    &
     &, wet_model_levels                                                &
     &, bl_levels                                                       &
     &, n_cca_levels                                                    &
                      ! Number of levels for conv cloud
                      ! amount: 1 for 2D, nlevs for 3D.
     &, tr_levels                                                       &
                      ! free tracer levels
     &, tr_vars                                                         &
                      ! number of free tracers
     &, tr_ukca       ! number of ukca tracers

! Model switches
      Integer, intent(in) ::                                            &
     &  model_domain

      Logical, intent(in) ::                                            &
     &  L_CAL360                                                        &
                   ! true if using 360 day calender
     &, Ltimer   ! true then output some timing information



      Logical, intent(in) ::                                            &
     &  L_calc_dxek                                                     &
                     ! IN Switch for calculation of condensate increment
     &, L_q_interact                                                    &
                     ! IN Switch allows overwriting of parcel variables
!                         when calculating condensate increments.
     &, L_mixing_ratio                                                  &
                       ! IN Use mixing ratio formulation
     &, l_murk_source                                                   &
                     ! IN Switch for murk source/scavanging
     &, L_DUST                                                          &
                     ! IN Switch for mineral dust
     &, L_sulpc_so2                                                     &
                     ! IN Switch for sulphur cycle
     &, L_sulpc_dms                                                     &
                     ! IN Switch for DMS in S-cycle
     &, L_sulpc_nh3                                                     &
                     ! IN Switch for NH3 in S-cycle
     &, L_soot                                                          &
                     ! IN Switch for soot cycle
     &, L_biomass                                                       &
                     ! IN Switch for biomass aerosol scheme
     &, L_ocff                                                          &
                     ! IN Switch for fossil-fuel organic carbon scheme
     &, L_co2_interactive                                               &
                     ! IN Switch for CO2 cycle
     &, L_USE_CARIOLLE 
                     ! IN true if Cariolle ozone tracer scheme being used

! physical constants
      Real                                                              &
     &  lc, lf, cp                                                      &
     &, two_Omega                                                       & 
                        ! twice Earth's rotation rate
     &, p_zero                                                          &
     &, kappa                                                           &
     &, R, g, Lapse_Rate, earth_radius, Pi

! model parameters
      Real,  intent(in) ::                                              &
     &  timestep                                                        &

     &, dbsdtbs_turb_0
                       ! IN PC2 erosion rate / s-1

      Real                                                              &
     &  Cloud_lifetime                                                  
                       ! Timescale for cloud decay (seconds), used cloud
                       ! decay on full convection timestep.


      !----------------------------------------------
      ! CCRad - Convection related variables
      !----------------------------------------------
      Integer, intent(out) :: lcbase_rad (row_length, rows)
      Integer, intent(out) :: lctop_rad  (row_length, rows)
      Real   , intent(out) :: ccw_rad    (row_length, rows, wet_model_levels)

      ! CCRad - Local variables and parameters
      Real                 :: cld_life_3d(row_length, rows, n_cca_levels)
      Integer, Parameter   :: cld_life_constant       = 0
      Integer, Parameter   :: cld_life_func_hgt       = 1
      Integer, Parameter   :: rad_decay_off           = 0
      Integer, Parameter   :: rad_decay_full_timestep = 1
      Integer, Parameter   :: rad_decay_conv_substep  = 2


#include "csubmodl.h"
#include "typsts.h"
! Parameters for soot cycle tracer scavenging
#include "c_st_con.h"
! Parameters for biomass aerosol scavenging
#include "c_bm_con.h"
! Parameters for fossil-fuel organic carbon aerosol scavenging
#include "c_ocff_con.h"
!

! Diagnostics info
      Real                                                              &
     &  STASHwork(*)  ! STASH workspace for section 4 (LS precip)


      Real , intent(inout) ::                                           &
     &  ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)

      Logical, intent(inout) ::         &
       cumulus(row_length, rows)        & ! Logical switch from boundary
                                          !    layer for presence of Cu
     , l_shallow(row_length, rows)      & ! Logical switch for shallow Cu
     , l_congestus(row_length, rows)    & ! Logical switch for congestus
     , l_congestus2(row_length, rows)   & ! congestus in descending air
     , l_pc2_diag_sh_pts(row_length,rows) ! Carry diagnostic shallow convective
                                          ! information for PC2

      Integer, intent(in) ::                                            &
     &       ntml(row_length, rows)                                     &
                                ! IN top level of surface mixed layer
     &,      ntpar(row_length, rows)                                    &
                                ! IN Top level of initial parcel
                                ! ascent  
     &,      ntdsc(row_length, rows)                                    &
                                ! IN top level of decoupled sc layer
     &,      nbdsc(row_length, rows)
                                ! IN bottom level of decoupled Sc layer

! Data arrays coming in from conv_diag & BL

      Real, intent(in) ::             &
       wstar(row_length, rows)        & ! Mixed layer convective velocity scale
     , wthvs(row_length, rows)        & ! Surface flux of THV
     , delthvu(row_length, rows)      & ! Integral of undilute parcel
                                        ! buoyancy over convective cloud layer
     , ql_ad(row_length, rows)        & ! adiabatic liquid water content
                                        ! at inversion (kg/kg)
     , ftl(row_length, rows)          & ! Surface sensible heat flux from BL
                                        ! (W/m2) i.e. cp*rho*w'tl'
     , fqt(row_length, rows)          & ! Total water flux from surface
                                        ! (kg/m2/s) i.e. rho*w'qT'
     , shallowc(row_length,rows)      & ! Indicator set to 1.0 if shallow,
                                        !  0.0 if not shallow or not cumulus
     , cu_over_orog(row_length,rows)  & ! Indicator for cumulus over steep
                                        ! orography. Indicator set to 1.0 if
                                        ! true, 0.0 if false. Exclusive.
     , cape_undilute(row_length,rows) & ! Undilute CAPE from parcel ascent
                                        ! m2/s2
     , cin_undilute(row_length,rows)    ! Undilute CIN from parcel ascent
                                        ! m2/s2

      Real , intent(in) ::                                              &
     & rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &        model_levels),                                            &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, p_star(row_length, rows)                                        &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &

     &, uw0_p(row_length, rows)                                         &
                                ! X-COMP SURFACE STRESS ON P GRID
     &, vw0_p(row_length, rows)                                         &
                                ! Y-COMP SURFACE STRESS ON P GRID
     &, w_max(row_length, rows)                                         &
                                !   max w in column
     &, ZLCL(row_length, rows)                                          &
                                ! LIFTING CONDENSATION LEVEL actual (m)
     &, ZLCL_UV(row_length, rows)                                       &
                                ! LIFTING CONDENSATION LEVEL for GWD
     &, ZTOP(row_length, rows)  ! TOP OF CLOUD LAYER

      Real, intent(in) ::                                               &
     &  p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, p_layer_centres(row_length, rows, 0:model_levels)               &
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = pstar, and at top = 0.
     &, exner_layer_boundaries(row_length, rows, 0:model_levels)        &
     &, exner_layer_centres(row_length, rows, 0:model_levels)

! heights above model surface in metres
       real, intent(in) ::                                              &
     &  z_theta(row_length, rows, model_levels)                         &
                                                    ! theta levels
     &, z_rho(row_length, rows, model_levels)                           &
                                                    ! rho levels
     &, rho_only(row_length, rows, model_levels)                        &
                                                    ! density (kg/m3)
     &, rho_theta(row_length, rows, model_levels-1) ! on th lev (kg/m3)

      logical, intent(in) ::                                            &
     &  land_sea_mask(row_length, rows)

      Real , intent(in) ::                                              &
     &  t1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, q1_sd(row_length, rows) ! set to zero initially

! Co-ordinate arrays
      Real, intent(in) ::                                               &
                          ! local vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, delta_lambda                                                    &
     &, delta_phi
! Information needed for STPH_SCV
      REAL,INTENT(IN) :: cos_theta_longitude(row_length,rows)
      REAL,INTENT(IN) :: FV_cos_theta_latitude(1-off_x:row_length+off_x &
     &                                         ,1-off_y:rows+off_y)
      REAL,INTENT(IN) :: f3_at_u(1-off_x:row_length+off_x, 1-off_y:     &
     &                           rows+off_y)

! time information for current timestep
      Integer, intent(in) ::                                            &
     &  val_year                                                        &
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second                                                      &
     &, timestep_number

! Diagnostic variables
      Real, intent(in) ::                                               &
     &  lat_rot_NP                                                      &
     &, long_rot_NP

! arguments with intent in/out. ie: input variables changed on output.
      Real, intent(inout) ::                                            &
     &  Theta_n(row_length, rows, model_levels)                         &
     &, q_n(row_length, rows, wet_model_levels)                         &
     &, qcl_n(row_length, rows, wet_model_levels)                       &
     &, qcf_n(row_length, rows, wet_model_levels)                       &
     &, cf_liquid_n(row_length, rows, wet_model_levels)                 &
     &, cf_frozen_n(row_length, rows, wet_model_levels)                 &
     &, bulk_cf_n(row_length, rows, wet_model_levels)                   &
     &, Theta_inc(row_length, rows, model_levels)                       &
     &, q_inc(row_length, rows, wet_model_levels)                       &
     &, qcl_inc(row_length, rows, wet_model_levels)                     &
     &, qcf_inc(row_length, rows, wet_model_levels)                     &
     &, cf_liquid_inc(row_length, rows, wet_model_levels)               &
     &, cf_frozen_inc(row_length, rows, wet_model_levels)               &
     &, bulk_cf_inc(row_length, rows, wet_model_levels)                 &
     &, R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        model_levels)

      Real, Intent(InOut) ::                                            &
     &  DUST_DIV1(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div1
     &, DUST_DIV2(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div2
     &, DUST_DIV3(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div3
     &, DUST_DIV4(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div4
     &, DUST_DIV5(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div5
     &, DUST_DIV6(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS) !dust mmr in div6

      Real, Intent(InOut) ::                                            &
     &  aerosol( 1 - off_x : row_length + off_x,                        &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, so2( 1 - off_x : row_length + off_x,                            &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, so4_aitken( 1 - off_x : row_length + off_x,                     &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, so4_accu( 1 - off_x : row_length + off_x,                       &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, so4_diss( 1 - off_x : row_length + off_x,                       &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, dms( 1 - off_x : row_length + off_x,                            &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, nh3( 1 - off_x : row_length + off_x,                            &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, soot_new( 1 - off_x : row_length + off_x,                       &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, soot_agd( 1 - off_x : row_length + off_x,                       &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, soot_cld( 1 - off_x : row_length + off_x,                       &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, bmass_new( 1 - off_x : row_length + off_x,                      &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, bmass_agd( 1 - off_x : row_length + off_x,                      &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, bmass_cld( 1 - off_x : row_length + off_x,                      &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, ocff_new( 1 - off_x : row_length + off_x,                       &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, ocff_agd( 1 - off_x : row_length + off_x,                       &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, ocff_cld( 1 - off_x : row_length + off_x,                       &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, co2( 1 - off_x : row_length + off_x,                            &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, OZONE_TRACER( 1 - off_x : row_length + off_x,                   &
     &           1 - off_y : rows + off_y, model_levels )               &
     &, free_tracers( 1 - off_x : row_length + off_x,                   &
     &           1 - off_y : rows + off_y, tr_levels, tr_vars)          &
     &, ukca_tracers( 1 - off_x : row_length + off_x,                   &
     &           1 - off_y : rows + off_y, tr_levels, tr_ukca)


! arguments with intent out. ie: output variables.
      Real, intent(inout) ::                                            &
     &  conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)                                     &
     &, conv_rain_3d(row_length, rows, wet_model_levels)                &
     &, conv_snow_3d(row_length, rows, wet_model_levels)                &
     &, cca_rad(row_length, rows, n_cca_levels)                         &
     &, cclwp_rad(row_length, rows) ! condensed water path (KG/M**2)    

! diagnostic arrays required for substepping to hold between calls

      Real, intent(inout) ::                                            &
     &  kterm_deep(row_length,rows)                                     &
                                    ! terminating level for deep
     &, precip_deep(row_length,rows)                                    &
                                      ! deep precipitation
     &, precip_shall(row_length,rows)                                   &
                                      ! shallow precipitation
     &, precip_mid(row_length,rows)                                     &
                                      ! mid-level precipitation
     &, precip_cong(row_length,rows)                                    &
                                      ! congestus precipitation
     &, shallow_ind(row_length,rows)                                    &
                                      ! shallow indicator
     &, congestus_ind(row_length,rows)                                  &
                                        ! congestus indicator
     &, congestus_ind2(row_length,rows)                                 &
                                         ! congestus indicator2
     &, mid_ind(row_length,rows)                                        &
                                        ! mid indicator
     &, ntml_diag(row_length,rows)                                      &
                                        ! ntml output diagnostic
     &, ntpar_diag(row_length,rows)                                     &
                                        ! ntpar output diagnostic
     &, freeze_diag(row_length,rows)                                    &
                                      ! freezing level output diagnostic
     &, wstar_up_diag(row_length,rows)                                  &
                                        ! cumulus layer conv vel scale
     &, wstar_dn_diag(row_length,rows)                                  &
                                        ! subcloud layer conv vel scale
     &, mb1_diag(row_length,rows)                                       &
                                        ! cloud base mass flux 1
     &, mb2_diag(row_length,rows)                                       &
                                        ! cloud base mass flux 2
     &, cape_out (row_length, rows)                                     &
                                       ! CAPE for diagnostic purposes
     &, cg_term(row_length, rows)                                       &
                                       ! congestus termination level
     &, cg_top(row_length, rows)                                        &
                                       ! height of congestus top (m)
     &, cg_base(row_length, rows)      ! height of congestus base (m)

! 5A extra diagnostics on model levels
      Real, intent(inout) ::                                            &
     &  mf_deep(row_length, rows, wet_model_levels)                     &
                                                    ! mass flux for deep
     &, mf_congest(row_length, rows, wet_model_levels)                  &
                                                         ! for congestus
     &, mf_shall (row_length, rows, wet_model_levels)                   &
                                                         ! for shallow
     &, mf_midlev(row_length, rows, wet_model_levels)                   &
                                                         ! for mid-level
     &, dt_deep(row_length, rows, wet_model_levels)                     &
                                                        ! dt for deep
     &, dt_congest(row_length, rows, wet_model_levels)                  &
                                                        ! dt for congest
     &, dt_shall(row_length, rows, wet_model_levels)                    &
                                                        ! dt for shallow
     &, dt_midlev(row_length, rows, wet_model_levels)                   &
                                                        ! dt for mid-lev
     &, dq_deep(row_length, rows, wet_model_levels)                     &
                                                        ! dq for deep
     &, dq_congest(row_length, rows, wet_model_levels)                  &
                                                        ! dq for congest
     &, dq_shall(row_length, rows, wet_model_levels)                    &
                                                        ! dq for shallow
     &, dq_midlev(row_length, rows, wet_model_levels)                   &
                                                        ! dq for mid-lev
     &, du_deep(row_length, rows, wet_model_levels)                     &
                                                        ! du for deep
     &, du_congest(row_length, rows, wet_model_levels)                  &
                                                        ! du for congest
     &, du_shall(row_length, rows, wet_model_levels)                    &
                                                        ! du for shallow
     &, du_midlev (row_length, rows, wet_model_levels)                  &
                                                        ! du for mid-lev
     &, dv_deep(row_length, rows, wet_model_levels)                     &
                                                        ! dv for deep
     &, dv_congest(row_length, rows, wet_model_levels)                  &
                                                        ! dv for congest
     &, dv_shall(row_length, rows, wet_model_levels)                    &
                                                        ! dv for shallow
     &, dv_midlev(row_length, rows, wet_model_levels)   ! dv for mid-lev

      Integer, intent(inout) ::                                         &
     &  ccb_rad(row_length, rows)                                       &
     &, cct_rad(row_length, rows)

      Integer ::                                                        &
     &  freeze_lev(row_length,rows)                                     &
                                     ! freezing level no
     &, it_kterm_deep(row_length,rows)                                  &
                                       !lev no for terminating deep con
     &, it_cg_term(row_length,rows)  !lev no for terminating congestus

      logical                                                           &
     &  it_mid_level(row_length,rows)  ! true if mid level convection

! Additional variables for SCM diagnostics which are dummy in full UM
      Integer, Intent(IN) ::                                            &
     &  nSCMDpkgs             ! No of SCM diagnostics packages

      Logical, Intent(IN) ::                                            &
     &  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

      Integer                                                           &
     &  Error_code

! local variables.

#include "parparm.h"
#include "domtyp.h"
#include "c_dust_ndiv.h"
#include "c_dustscav.h"
#include "cprintst.h"
#include "c_pc2pdf.h"
#if defined(SCMA)
! INOUT SCMop is declared in here
#include "s_scmop.h"
#endif

! loop counters
      Integer                                                           &
     &  i, j, k, ii, jj                                                 &
     &, call_number                                                     &
     &, loop_number

! local variables
      Integer                                                           &
     &  levels_h                                                        &
     &, levels_v                                                        &
     &, n_conv_levels                                                   &
     &, n_conv_points                                                   &
     &, this_point                                                      &
     &, info                                                            &
     &, step                                                            &
     &, seg_points                                                      &
     &, first_point

! DECLARE LOCAL VARIABLES FOR STPH_SCV
!
      INTEGER :: minutos
!
      INTEGER :: NMCSMAX    ! MAX NUMBER OF VORTICITY SOURCES ON EACH PE
      PARAMETER(NMCSMAX = 50)
!
      INTEGER :: ISCBFREQ ! FREQUENCY OF UPDATING SCB IN HOURS
      PARAMETER(ISCBFREQ = 3)
!
      REAL :: SCB_MINRAD ! MINIMUM RADIUS FOR VORTICITY SOURCE
      PARAMETER(SCB_MINRAD = 50000.)
!
      INTEGER :: NMCSPROP
      PARAMETER(NMCSPROP=7) ! NO. OF PROPERTIES IN SCB_MCS_PROP
!
      LOGICAL :: L_SCV_UP   ! LOGICAL TO DETERMINE WHETHER TO UPDATE SCB
!
      REAL :: SCB_MCS_PROP(NMCSMAX*N_PROC,NMCSPROP)
! CONTAINS THE PROPERTIES OF THE VORTICITY SOURCES AS FOLLOWS:
!   1 - LATITUDE
!   2 - LONGITUDE
!   3 - RADIUS (M)
!   4 - MCV LEVEL (FREEZING LEVEL) (PA)
!   5 - ANTICYCLONE LEVEL (CCT) (PA)
!   6 - MCV DEPTH (PA)
!   7 - ANTICYCLONE DEPTH (PA)

!
      INTEGER :: ISCBSEED    ! RANDOM SEED
!
      REAL :: SCB_RAD     ! TEMPORARY SCALAR FOR VORTICITY SOURCE
!                      ! RADIUS CALCULATION
      REAL :: coriolis(row_length,rows)
      REAL :: SCB_FCORIOL ! TEMPORARY SCALAR FOR LOCAL CORIOLIS FORCE
!
      INTEGER :: NMCS,NMCSTOT ! VARIABLE COUNTING NUMBER OF VORT.SOURCES
!
      INTEGER :: ISCBCHANGE ! FLAG USED IN THE SORTING ALGORITHM
!
      REAL :: LAT1,LAT2,LONG1,LONG2,DIST    ! SCALARS USED CALCULATING
!                                           ! DISTANCE BETWEEN 2 POINTS
      INTEGER :: N,M  ! ADDITIONAL LOOP COUNTERS
!
! VARIABLES DEFINING THE STRUCTURE OF THE VORTICITY SOURCES
!
      REAL :: SCB_TOR ! TIMESCALE OVER VORTICITY INCREMENTS APPLIED
!
      REAL :: TSC             ! RECIPROCAL OF SCB_TOR
!
      REAL :: ZDEPTH_A,ZDEPTH_C ! HALF DEPTH OF (ANTI)CYCLONE IN M
!
      REAL :: PDEPTH_A,PDEPTH_C ! HALF DEPTH OF (ANTI)CYCLONE IN PA
!
      REAL :: C_ETA,RCETA,RCETA2 ! C_ETA IS THE RATIO OF THE CYCLONE RAD
                               ! TO THE ANTICYCLONE RADIUS. RCETA AND
                               ! RCETA2 ARE THE RECIPROCAL AND
                               ! RECIPROCAL SQUARED RESPECTIVELY.
!
      REAL :: DADV    ! RATIO OF ANTICYLONE AND CYCLONE HALF DEPTHS
      REAL :: XDIST           ! DISTANCE IN LONGITUDE ONLY
      REAL :: YDIST           ! DISTANCE IN LATITUDE ONLY
      REAL :: PRESSURE        ! LOCAL PRESSURE (PA)
      REAL :: PFAC            ! HEIGHT DEPENDENT FACTOR CONTROLLING
                              ! RADIUS OF VORTICITY SOURCE.
      REAL :: PMID_C,PMID_A   ! PRESSURE AT CENTRE OF MCV/ANTICYCLONE
      REAL :: TEMPERATURE     ! A TEMP TEMPERATURE VARIABLE (K)
      REAL :: SCB_TEMP        ! TEMPORARY SCALAR IN SHUFFLE ALGORITHM
      REAL :: SCB_DIST        ! declare real function scb_dist
!
      REAL :: SCB_DUDT(row_length,rows,model_levels)  ! SCV U INCREMENT
      REAL :: SCB_DVDT(row_length,rows,model_levels)  ! SCV V INCREMENT
!
      REAL :: azar ! Variables for setting up the random seed
! Variable to avoid conflict with k (loop)
      INTEGER :: L
!
! More variables needed for STPH_SCV: latitude, longitude and coriolis
       REAL :: latitude(row_length,rows)
       REAL :: longitude(row_length,rows)
!

      Logical ::                                                        &
     &  cumulus_1d(rowsrowlength)                                       &
                                       ! 1-d version of CUMULUS
     &, l_shallow_1d(rowsrowlength)                                     &
                                       ! 1-d version of SHALLOW_BL
     &, l_congestus_1d(rowsrowlength)  ! 1-d version of congestus


      Logical                                                           &
     &  L_copy_top                                                      &
     &, L_tracer                                                        &
                         ! Switch for tracer variables used
     &, L_full_zero      ! True if a dummy zero full field is required

! Local data arrays
!    Check Variables
      Integer                                                           &
     &  n_qcx(1+wet_model_levels,2) ! Record  out-of-range input events
!
!
! Allocatable arrays for PC2 scheme increment calculation - used to save
! memory when increment calculations are not requested.
      Real,Dimension(:,:,:),Allocatable::                               &
     &  T_earliest                                                      &
                         !  temperature at start of convection
     &, T_inc_latest                                                    &
                         !  temperature increment on input to PC2 homog
     &, q_earliest                                                      &
                         !  humidity    at start of convection
     &, qcl_earliest                                                    &
                         !  qCL         at start of convection
     &, cfl_earliest                                                    &
                         !  cf_liquid   at start of convection
     &, cff_earliest                                                    &
                         !  cf_frozen   at start of convection
     &, bcf_earliest                                                    &
                         !  bulk cloud  at start of convection
     &, theta_inc_PC2                                                   &
                         !  pot temperature increment due to PC2 homog
     &, q_inc_PC2                                                       &
                         !  humidity        increment due to PC2 homog
     &, qcl_inc_PC2                                                     &
                         !  qCL             increment due to PC2 homog
     &, cfl_inc_PC2                                                     &
                         !  cf_liquid       increment due to PC2 homog
     &, bcf_inc_PC2      !  bulk cloud      increment due to PC2 homog
!
      Real,Dimension(:,:,:),Allocatable::                               &
     &  full_zero        !  a dummy array for a zero field

! Allocatable arrays for diagnostic variables - required to save memory
! when diagnostic not requested
      Real,Dimension(:,:,:),Allocatable::                               &
     & qcl_incr_inhom_diag                                              &
                             ! qCL         increment for STASH
     &,qcf_incr_inhom_diag                                              &
                             ! qCF         increment for STASH
     &,bulk_cf_incr_inhom_diag                                          &
                                   ! bcf   increment for STASH
     &,cf_liquid_incr_inhom_diag                                        &
                                   ! cf_l  increment for STASH
     &,cf_frozen_incr_inhom_diag   ! cf_f  increment for STASH


! Diagnostics that need to be kept between substeps
      Real                                                              &
     & T_incr_diagnostic(row_length,rows,model_levels)                  &
                            ! temperature increment for STASH
     &,q_incr_diagnostic(row_length,rows,wet_model_levels)              &
                            ! humidity    increment for STASH
     &,qcl_incr_diagnostic(row_length,rows,wet_model_levels)            &
     &,qcf_incr_diagnostic(row_length,rows,wet_model_levels)            &
     &,u_incr_diagnostic(row_length,rows,model_levels)                  &
                            ! u wind      increment for STASH
     &,v_incr_diagnostic(row_length,n_rows,model_levels)                &
                            ! v wind      increment for STASH
     &,cf_liquid_incr_diagnostic(row_length,rows,wet_model_levels)      &
     &,cf_frozen_incr_diagnostic(row_length,rows,wet_model_levels)      &
     &,bulk_cf_incr_diagnostic(row_length,rows,wet_model_levels)
! Local arrays holding information to be passed between physics
! routines.

! Diagnostics controlled by Diagnostic switches

! Convection
      Real                                                              &
     &  lcca(row_length, rows)                                          &
                                ! lowest convective cloud amount (%)
     &, lcclwp(row_length, rows)                                        &
                                 ! lowest cloud condensed water
                                  ! path (KG/M**2)
     &, ccw(row_length, rows, wet_model_levels)                         &
                                  ! convective cloud liquid water
                                  ! (G/KG) on model levels
     &, dthbydt(row_length, rows, wet_model_levels)                     &
     &, dqbydt(row_length, rows, wet_model_levels)                      &
     &, dqclbydt(row_length, rows, wet_model_levels)                    &
                                                     ! Q4 Increment qCL
     &, dqcfbydt(row_length, rows, wet_model_levels)                    &
                                                     ! Q4 Increment qCF
     &, dcflbydt(row_length, rows, wet_model_levels)                    &
                                                     ! Cloud Increment
     &, dcffbydt(row_length, rows, wet_model_levels)                    &
                                                     ! Cloud Increment
     &, dbcfbydt(row_length, rows, wet_model_levels)                    &
                                                     ! Cloud Increment
     &, dubydt_p(row_length, rows, model_levels)                        &
     &, dvbydt_p(row_length, rows, model_levels)                        &
     &, dubydt_u(row_length, rows, model_levels)                        &
     &, dvbydt_v(row_length, n_rows, model_levels)                      &
     &, u_p(row_length, rows, model_levels)                             &
     &, v_p(row_length, rows, model_levels)                             &
     &, work_p_halo(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &              model_levels)                                       &
     &, work_v_halo(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,     &
     &              model_levels)                                       &
     &, Theta_diag(row_length, rows, model_levels)                      &
     &, q_diag(row_length, rows, wet_model_levels)



! mineral dust diagnostics
      REAL                                                              &
     &  CONSCAV_DUST(ROW_LENGTH, ROWS, 6) !COL TOTAL SCVNGD FOR
                                          !DIV1,DIV2,...,DIV6 DUST

! Diagnostics for S Cycle
      Real                                                              &
     &   conwash_so2(row_length, rows)                                  &
                                          !column total scvngd so2
     &,  conscav_so4ait(row_length, rows)                               &
                                          !column total scvngd so4ait
     &,  conscav_so4acc(row_length, rows)                               &
                                          !column total scvngd so4acc
     &,  conscav_so4dis(row_length, rows)                               &
                                          !column total scvngd so4dis
     &,  conwash_nh3(row_length, rows)    !column total scvngd nh3
!
! Diagnostics for soot cycle
      Real                                                              &
     &   conscav_agedsoot(row_length, rows) !colm total scvngd agd soot
! Diagnostics for biomass aerosol scheme
      Real                                                              &
     &   conscav_agedbmass(row_length, rows) !colm total scvngd bmass
! Diagnostics for fossil-fuel organic carbon scheme
      Real                                                              &
     &   conscav_agedocff(row_length, rows) !colm total scvngd ocff

      Integer                                                           &
     &  lcbase (row_length, rows)                                       &
                                   ! lowest convective cloud base level
     &, lctop (row_length, rows)  ! lowest convective cloud top level

!  allocatable array holding copies of all tracers
      Real, Dimension(:,:,:,:), Allocatable ::                          &
     &  tot_tracer

#if defined(SCMA)
!    Work array for calculating temperature at end of step
      Real                                                              &
     &  work_3d(row_length,rows,model_levels)                           &
     &, work_2d(row_length,rows)
#endif

!  sizes for temp. tracer calculations
      Integer ::                                                        &
     &  ntra_fld                                                        &
     &, ntra_tmp                                                        &
     &, ntra_lev
!
! Parameters for S Cycle tracer scavenging
      Real                                                              &
     &    krain_so4ait                                                  &
     &,   ksnow_so4ait                                                  &
     &,   krain_so4acc                                                  &
     &,   ksnow_so4acc                                                  &
     &,   krain_so4dis                                                  &
     &,   ksnow_so4dis
!
      Parameter(                                                        &
     &    krain_so4ait=0.3E-4                                           &
     &,   ksnow_so4ait=0.3E-4                                           &
     &,   krain_so4acc=0.3E-4                                           &
     &,   ksnow_so4acc=0.3E-4                                           &
     &,   krain_so4dis=0.3E-4                                           &
     &,   ksnow_so4dis=0.3E-4                                           &
     & )
!

      Real                                                              &
     &  it_lcca(row_length, rows)                                       &
                                   ! lowest convective cloud amount (%)
     &, it_lcclwp(row_length, rows)                                     &
                                    ! lowest cloud condensed water
                                  ! path (KG/M**2)
     &, it_ccw(row_length, rows, wet_model_levels)
                                  ! convective cloud liquid water
                                  ! (G/KG) on model levels
      Real                                                           :: &
     &  it_conv_rain    (row_length, rows)                              &
     &, it_conv_snow    (row_length, rows)                              &
     &, it_conv_rain_3d (row_length, rows, wet_model_levels)            &
     &, it_conv_snow_3d (row_length, rows, wet_model_levels)            &
     &, it_cca          (row_length, rows, n_cca_levels)                &
     &, cca             (row_length, rows, n_cca_levels)                &
     &, it_cclwp        (row_length, rows)                              &
     &, cclwp           (row_length, rows)                              &
                                   ! condensed water path (KG/M**2)
     &, it_cape_out(row_length, rows)                                   &
                                      ! CAPE
     &, it_precip_dp(row_length, rows)                                  &
                                       ! deep precip
     &, it_precip_sh(row_length, rows)                                  &
                                       ! shallow precip
     &, it_precip_md(row_length, rows)                                  &
                                       ! mid-level precip
     &, it_precip_cg(row_length, rows)                                  &
                                       ! congestus precip
     &, it_wstar_dn(row_length, rows)                                   &
     &, it_wstar_up(row_length, rows)                                   &
     &, it_mb1(row_length, rows)                                        &
     &, it_mb2(row_length, rows)                                        &
     &, it_up_flux_half(row_length, rows, model_levels)
                                       !up flux on half levs.

      Integer                                                           &
     &  it_ccb (row_length, rows)                                       &
     &, it_cct (row_length, rows)                                       &
     &, ccb (row_length, rows)                                          &
     &, cct (row_length, rows)

      Integer                                                           &
     &  it_lcbase (row_length, rows)                                    &
                                ! lowest convective cloud base level
     &, it_lctop (row_length, rows)
                                ! lowest convective cloud top level

      Logical                                                           &
     & flg_up_flx,                                                      &
                                ! stash flag for updraught mass flux
     & flg_up_flx_half,                                                 &
                                ! stash flag for updraught mass flux
     & flg_dwn_flx,                                                     &
                                ! stash flag for downdraught mass flux
     & flg_entr_up,                                                     &
                                ! stash flag for updraught entrainment
     & flg_entr_dwn,                                                    &
                                ! stash flag for downdraught entrainmn
     & flg_detr_up,                                                     &
                                ! stash flag for updraught detrainment
     & flg_detr_dwn,                                                    &
                                ! stash flag for downdraught detrainmn
     & flg_conv_rain_3d,                                                &
                                ! stash flag for 3d conv rainfall
     & flg_conv_snow_3d,                                                &
                                ! stash flag for 3d conv snowfall
! For A05_4A:
     & flg_uw_dp,                                                       &
                                ! stash flag for deep x-comp stress
     & flg_vw_dp,                                                       &
                                ! stash flag for deep y-comp stress
     & flg_uw_shall,                                                    &
                                ! stash flag for shallow x stress
     & flg_vw_shall             ! stash flag for shallow y stress

! For 5A turbulence based schemes
       logical ::                                                       &
     & flg_wqt_flux                                                     & 
                               ! stash flag for w'qt'
     &,flg_wql_flux                                                     & 
                               ! stash flag for w'ql'
     &,flg_wthetal_flux                                                 & 
                               ! stash flag for w'thetal'
     &,flg_wthetav_flux                                                 & 
                               ! stash flag for w'thetav'
     &,flg_mf_deep                                                      &
                               ! stash flag for deep mass flux
     &,flg_mf_congest                                                   &
                               ! stash flag for congestus mass flux
     &,flg_mf_shall                                                     &
                               ! stash flag for shallow mass flux
     &,flg_mf_midlev                                                    &
                               ! stash flag for mid_level mass flux
     &,flg_dt_deep                                                      &
                               ! stash flag for deep dT
     &,flg_dt_congest                                                   &
                               ! stash flag for congestus dT
     &,flg_dt_shall                                                     &
                               ! stash flag for shallow dT
     &,flg_dt_midlev                                                    &
                               ! stash flag for mid_level dT
     &,flg_dq_deep                                                      &
                               ! stash flag for deep dq
     &,flg_dq_congest                                                   &
                               ! stash flag for congestus dq
     &,flg_dq_shall                                                     &
                               ! stash flag for shallow dq
     &,flg_dq_midlev                                                    &
                               ! stash flag for mid_level dq
     &,flg_du_deep                                                      &
                               ! stash flag for deep du
     &,flg_du_congest                                                   &
                               ! stash flag for congestus du
     &,flg_du_shall                                                     &
                               ! stash flag for shallow du
     &,flg_du_midlev                                                    &
                               ! stash flag for mid_level du
     &,flg_dv_deep                                                      &
                               ! stash flag for deep dv
     &,flg_dv_congest                                                   &
                               ! stash flag for congestus dv
     &,flg_dv_shall                                                     &
                               ! stash flag for shallow dv
     &,flg_dv_midlev           ! stash flag for mid_level dv

! STASHflag switches for increment diagnostics:
      Logical                                                           &
     & l_qcl_incr_cinh                                                  &
                             ! liquid cloud condensate qCL (inhom)
     &,l_qcf_incr_cinh                                                  &
                             ! frozen cloud condensate qCF (inhom)
     &,l_cfl_incr_cinh                                                  &
                             ! liquid cloud amount cf_liquid (inhom)
     &,l_cff_incr_cinh                                                  &
                             ! frozen cloud amount cf_frozen (inhom)
     &,l_bcf_incr_cinh                                                  &
                             ! total (bulk) cloud amount bulk_cf (inhom)
     &,l_T_incr_conv                                                    &
                             ! temperature
     &,l_q_incr_conv                                                    &
                             ! humidity
     &,l_qcl_incr_conv                                                  &
                             ! liquid cloud condensate qCL
     &,l_qcf_incr_conv                                                  &
                             ! frozen cloud condensate qCF
     &,l_cfl_incr_conv                                                  &
                             ! liquid cloud amount cf_liquid
     &,l_cff_incr_conv                                                  &
                             ! frozen cloud amount cf_frozen
     &,l_bcf_incr_conv                                                  &
                             ! total (bulk) cloud amount bulk_cf
     &,l_u_incr_conv                                                    &
                             ! u wind
     &,l_v_incr_conv                                                    &
                             ! v wind
     &,l_up_incr_conv                                                   &
                             ! u wind
     &,l_vp_incr_conv        ! v wind

!    Dummy arrays for convect diagnostics : will have to
!    be passed to the diagnostics routine for stash processing
!    (Their above respective flags will then have to be turned on)
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
     &, it_uw_dp(row_length, rows, model_levels)                        &
     &, it_vw_dp(row_length, rows, model_levels)                        &
     &, it_uw_shall(row_length, rows, model_levels)                     &
     &, it_vw_shall(row_length, rows, model_levels)
      real, intent(inout) ::                                            &
     &  dubydt_pout(row_length, rows, model_levels)                     &
                                                    ! store mean inc
     &, dvbydt_pout(row_length, rows, model_levels)                     &
                                                    ! store mean inc
     &, up_flux(row_length, rows, model_levels)                         &
                                                !  updraught mass flux
     &, up_flux_half(row_length, rows, model_levels)                    &
                                !  updraught mass flux on half levs
     &, dwn_flux(row_length, rows, model_levels)                        &
                                                 ! downdraught mass flux
     &, entrain_up(row_length, rows, model_levels)                      &
                                ! fractional entrainment rate
     &, detrain_up(row_length, rows, model_levels)                      &
                                ! fractional detrainment rate
     &, entrain_dwn(row_length, rows,model_levels)                      &
                                ! fractional entrainment rate
                                ! into downdraughts
     &, detrain_dwn(row_length, rows,model_levels)
                                ! fractional detrainment rate
                                ! from downdraughts
       real ::                                                          &
     &  it_up_flux(row_length, rows, model_levels)                      &
     &, it_dwn_flux(row_length, rows, model_levels)                     &
     &, it_entrain_up(row_length, rows, model_levels)                   &
     &, it_detrain_up(row_length, rows, model_levels)                   &
     &, it_entrain_dwn(row_length, rows, model_levels)                  &
     &, it_detrain_dwn(row_length, rows, model_levels)                  &
     &, r_rho_lev(row_length, rows, model_levels)                       &
     &, r_theta_lev(row_length, rows, 0:model_levels)                   &
     &, one_over_conv_calls                                             &
     &, timestep_conv

!-----------------------------------------------------------------------
!  variables required by HadGEM1 cloud calculation only

      Real ::                                                           &
     &  it_mf_deep(row_length, rows, wet_model_levels)                  &
     &, it_mf_congest (row_length, rows, wet_model_levels)              &
     &, it_mf_shall(row_length, rows, wet_model_levels)                 &
     &, it_mf_midlev(row_length, rows, wet_model_levels)                &
     &, it_dt_deep(row_length, rows, wet_model_levels)                  &
     &, it_dt_congest(row_length, rows, wet_model_levels)               &
     &, it_dt_shall(row_length, rows, wet_model_levels)                 &
     &, it_dt_midlev(row_length, rows, wet_model_levels)                &
     &, it_dq_deep(row_length, rows, wet_model_levels)                  &
     &, it_dq_congest(row_length, rows, wet_model_levels)               &
     &, it_dq_shall(row_length, rows, wet_model_levels)                 &
     &, it_dq_midlev (row_length, rows, wet_model_levels)               &
     &, it_du_deep (row_length, rows, wet_model_levels)                 &
     &, it_du_congest(row_length, rows, wet_model_levels)               &
     &, it_du_shall(row_length, rows, wet_model_levels)                 &
     &, it_du_midlev (row_length, rows, wet_model_levels)               &
     &, it_dv_deep (row_length, rows, wet_model_levels)                 &
     &, it_dv_congest(row_length, rows, wet_model_levels)               &
     &, it_dv_shall (row_length, rows, wet_model_levels)                &
     &, it_dv_midlev (row_length, rows, wet_model_levels)

      real, intent(out) ::                                              &
     & wqt_flux(row_length, rows, wet_model_levels)                     & 
                                                     ! w'qt' flux
     &,wql_flux(row_length, rows, wet_model_levels)                     & 
                                                     ! w'ql' flux
     &,wthetal_flux(row_length, rows, wet_model_levels)                 & 
                                                        ! w'thetal' flux
     &,wthetav_flux(row_length, rows, wet_model_levels)                 & 
                                                        ! w'thetav' flux
     &,wqt_cb(row_length, rows)                                         & 
                                    ! w'qt' cloud base
     &,wthetal_cb(row_length, rows)                                     & 
                                    ! w'thetal' cloud base
     &,wqt_inv(row_length, rows)                                        & 
                                     ! w'qt' base of inversion
     &,wthetal_inv(row_length, rows)                                    & 
                                     ! w'thetal' base of inversion
     &,sh_top(row_length, rows)                                         &
                                ! height of top of shallow convection(m)
     &,sh_base(row_length, rows) ! height of base of shall convection(m)


      real ::                                                           &
     & it_wqt_flux(row_length, rows, wet_model_levels)                  &
     &,it_wql_flux(row_length, rows, wet_model_levels)                  &
     &,it_wthetal_flux(row_length, rows, wet_model_levels)              &
     &,it_wthetav_flux(row_length, rows, wet_model_levels)
!
      real ::                                                           &
     &  mb_temp                                                         &
                        ! cloud base mass flux
     &, wsc_temp                                                        &
                        ! convective velocity scale
     &, cca_2d(row_length,rows)

      logical ::                                                        &
     &  cct_change_map(row_length,rows) ! where we need to correct cct

! Segmentation variables
      Integer, Allocatable  :: sp_local(:)
      Integer, Allocatable  :: fp_local(:)
      Integer, Allocatable  :: fyp(:)
      Integer, Allocatable  :: fxp(:)
      Integer, Allocatable  :: n_cumulus(:)   ! Number of CUMULUS points
      Integer, Allocatable  :: n_deep(:)      ! Number of DEEP points
      Integer, Allocatable  :: n_shallow(:)   ! Number of SHALLOW points
      Integer, Allocatable  :: n_congestus(:) ! Number of congestus points

      real, parameter :: a_land=0.3
      real, parameter :: a_sea=0.3
      real, parameter :: b_land=0.025
      real, parameter :: b_sea=0.025
      real, parameter :: reduction_factor=0.75

!-----------------------------------------------------------------------


!     Check variables used to ensure that moist variables are sensible.
      Integer ::                                                        &
     &  n_test
!
! Error reporting variables
      Integer ::                                                        &
     &  ErrorStatus

      Logical                                                           &
     &  L_apply_diag

      Character (Len=80) ::                                             &
     &  Cmessage

      Character (Len=*), Parameter ::                                   &
     &  RoutineName='NI_conv_ctl'
!
!-----------------------------------------------------------------------

! Sensible names for STASHflag switches on increment diagnostics
#if !defined(SCMA)
! Apply diags at last cycle only
      L_apply_diag = CycleNo == NumCycles
      flg_conv_rain_3d = sf(227,5)
      flg_conv_snow_3d = sf(228,5)
      flg_dwn_flx = sf(251,5) .and. L_apply_diag
      flg_entr_up = sf(252,5) .and. L_apply_diag
      flg_detr_up = sf(253,5) .and. L_apply_diag
      flg_entr_dwn = sf(254,5) .and. L_apply_diag
      flg_detr_dwn = sf(255,5) .and. L_apply_diag
      flg_uw_dp = sf(258,5) .and. L_apply_diag
      flg_vw_dp = sf(259,5) .and. L_apply_diag
      flg_uw_shall = sf(260,5) .and. L_apply_diag
      flg_vw_shall = sf(261,5) .and. L_apply_diag

        flg_wqt_flux = sf(290,5).or.sf(304,5).or.sf(306,5)
        flg_wql_flux = sf(291,5)
        flg_wthetal_flux = sf(292,5).or.sf(305,5).or.sf(307,5)
        flg_wthetav_flux = sf(293,5)
        flg_mf_deep    = sf(320,5)
        flg_mf_congest = sf(321,5)
        flg_mf_shall   = sf(322,5)
        flg_mf_midlev  = sf(323,5)
        flg_dt_deep    = sf(324,5)
        flg_dt_congest = sf(325,5)
        flg_dt_shall   = sf(326,5)
        flg_dt_midlev  = sf(327,5)
        flg_dq_deep    = sf(328,5)
        flg_dq_congest = sf(329,5)
        flg_dq_shall   = sf(330,5)
        flg_dq_midlev  = sf(331,5)
        flg_du_deep    = sf(332,5)
        flg_du_congest = sf(333,5)
        flg_du_shall   = sf(334,5)
        flg_du_midlev  = sf(335,5)
        flg_dv_deep    = sf(336,5)
        flg_dv_congest = sf(337,5)
        flg_dv_shall   = sf(338,5)
        flg_dv_midlev  = sf(339,5)
      flg_up_flx = sf(250,5)
      flg_up_flx_half = ( sf(249,5) .OR. sf(246,5) ) .and. L_apply_diag
      l_qcl_incr_cinh = sf(163,5) .and. L_apply_diag
      l_qcf_incr_cinh = sf(164,5) .and. L_apply_diag
      l_bcf_incr_cinh = sf(172,5) .and. L_apply_diag
      l_cfl_incr_cinh = sf(173,5) .and. L_apply_diag
      l_cff_incr_cinh = sf(174,5) .and. L_apply_diag
      l_T_incr_conv = ( ( sf(181,5) .or. sf(187,5) ) .and. L_apply_diag ) .or. L_calc_dxek
      l_q_incr_conv = ( sf(182,5) .or. sf(188,5) ) .and. L_apply_diag
      l_qcl_incr_conv = sf(183,5) .and. L_apply_diag
      l_qcf_incr_conv = sf(184,5) .and. L_apply_diag
      l_u_incr_conv = sf(185,5) .and. L_apply_diag
      l_v_incr_conv = sf(186,5) .and. L_apply_diag
      l_up_incr_conv = sf(256,5) .and. L_apply_diag
      l_vp_incr_conv = sf(257,5) .and. L_apply_diag
      l_bcf_incr_conv = sf(192,5) .and. L_apply_diag
      l_cfl_incr_conv = sf(193,5) .and. L_apply_diag
      l_cff_incr_conv = sf(194,5) .and. L_apply_diag
#else
! Always apply diags in SCM when stash flags are on.
      L_apply_diag     = .TRUE.
      flg_conv_rain_3d = .TRUE.
      flg_conv_snow_3d = .TRUE.
      flg_up_flx_half  = .TRUE.
      flg_dwn_flx      = .TRUE.
      flg_entr_up      = .TRUE.
      flg_detr_up      = .TRUE.
      flg_entr_dwn     = .TRUE.
      flg_detr_dwn     = .TRUE.
      flg_uw_dp        = .TRUE.
      flg_vw_dp        = .TRUE.
      flg_uw_shall     = .TRUE.
      flg_vw_shall     = .TRUE.
      flg_up_flx       = .TRUE.
      l_qcl_incr_cinh  = .TRUE.
      l_qcf_incr_cinh  = .TRUE.
      l_bcf_incr_cinh  = .TRUE.
      l_cfl_incr_cinh  = .TRUE.
      l_cff_incr_cinh  = .TRUE.
      l_T_incr_conv    = .TRUE.
      l_q_incr_conv    = .TRUE.
      l_qcl_incr_conv  = .TRUE.
      l_qcf_incr_conv  = .TRUE.
      l_u_incr_conv    = .TRUE.
      l_v_incr_conv    = .TRUE.
      l_up_incr_conv   = .FALSE.
      l_vp_incr_conv   = .FALSE.
      l_bcf_incr_conv  = .TRUE.
      l_cfl_incr_conv  = .TRUE.
      l_cff_incr_conv  = .TRUE.

#if defined(A05_5A)
      ! Calculate 5A diagnostics so available if
      ! user wishes to add lines for their output
      flg_mf_deep    = .true.
      flg_mf_congest = .true.
      flg_mf_shall   = .true.
      flg_mf_midlev  = .true.
      flg_dt_deep    = .true.
      flg_dt_congest = .true.
      flg_dt_shall   = .true.
      flg_dt_midlev  = .true.
      flg_dq_deep    = .true.
      flg_dq_congest = .true.
      flg_dq_shall   = .true.
      flg_dq_midlev  = .true.
      flg_du_deep    = .true.
      flg_du_congest = .true.
      flg_du_shall   = .true.
      flg_du_midlev  = .true.
      flg_dv_deep    = .true.
      flg_dv_congest = .true.
      flg_dv_shall   = .true.
      flg_dv_midlev  = .true.
      flg_wqt_flux   = .true.
      flg_wql_flux   = .true.
      flg_wthetal_flux = .true.
      flg_wthetav_flux = .true.
#else
      ! Diagnostics for 5A scheme
      flg_mf_deep    = .false.
      flg_mf_congest = .false.
      flg_mf_shall   = .false.
      flg_mf_midlev  = .false.
      flg_dt_deep    = .false.
      flg_dt_congest = .false.
      flg_dt_shall   = .false.
      flg_dt_midlev  = .false.
      flg_dq_deep    = .false.
      flg_dq_congest = .false.
      flg_dq_shall   = .false.
      flg_dq_midlev  = .false.
      flg_du_deep    = .false.
      flg_du_congest = .false.
      flg_du_shall   = .false.
      flg_du_midlev  = .false.
      flg_dv_deep    = .false.
      flg_dv_congest = .false.
      flg_dv_shall   = .false.
      flg_dv_midlev  = .false.
      flg_wqt_flux   = .false.
      flg_wql_flux   = .false.
      flg_wthetal_flux = .false.
      flg_wthetav_flux = .false.
#endif
#endif

        timestep_conv=timestep/(n_conv_calls*1.0)

      If (error_code  ==  0 ) Then
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Convect ',3)
!
      L_full_zero = L_calc_dxek
!
! L_full_zero_if1:
      If (L_full_zero) Then
!     Set up a field of zero values if dummy field is required
!
        Allocate ( full_zero(row_length,rows,1) )
!
        Do J=1, rows
          Do I=1, row_length
            full_zero(i,j,1)  = 0.0
          End Do  ! I
        End Do  ! J
!
      End If  ! L_full_zero_if1
!
! ----------------------------------------------------------------------
! Prevent negative condensate problems by condensing vapour if needed.
! Input fields are updated without updating increments so change will
! only affect PC2 increments and not the top-level QX_STAR fields.
! ----------------------------------------------------------------------
! L_calc_dxek_if0:
      If (L_calc_dxek) Then
!
        n_qcx(1+wet_model_levels,1) = 0
        n_qcx(1+wet_model_levels,2) = 0
!
        Do K=1, wet_model_levels
!
          n_test = 0
          Do J=1, rows
            Do I=1, row_length
              If(qcf_n(I,J,K)  <   0.0) Then
!             Freeze some liquid to zero ice content
                qcl_n(I,J,K) = qcl_n(I,J,K) + qcf_n(I,J,K)
                Theta_n(I,J,K) = Theta_n(I,J,K) -                       &
     &        ( (qcf_n(I,J,K) * LF) / (CP * exner_theta_levels(I,J,K)) )
                qcf_n(I,J,K) = 0.0
                cf_frozen_n(I,J,K) = 0.0
                bulk_cf_n(I,J,K) = cf_liquid_n(I,J,K)
                n_test = n_test + 1
              End If
            End Do  ! I
          End Do  ! J
          n_qcx(K,1) = n_test
          n_qcx(1+wet_model_levels,1) = n_qcx(1+wet_model_levels,1) +   &
     &                                  n_qcx(K,1)
!
          n_test = 0
          Do J=1, rows
            Do I=1, row_length
              If(qcl_n(I,J,K)  <   0.0) Then
!             Condense some vapour to zero liquid content
                q_n(I,J,K) = q_n(I,J,K) + qcl_n(I,J,K)
                Theta_n(I,J,K) = Theta_n(I,J,K) -                       &
     &        ( (qcl_n(I,J,K) * LC) / (CP * exner_theta_levels(I,J,K)) )
                qcl_n(I,J,K) = 0.0
                cf_liquid_n(I,J,K) = 0.0
                bulk_cf_n(I,J,K) = cf_frozen_n(I,J,K)
                n_test = n_test + 1
              End If
            End Do  ! I
          End Do  ! J
          n_qcx(K,2) = n_test
          n_qcx(1+wet_model_levels,2) = n_qcx(1+wet_model_levels,2) +   &
     &                                  n_qcx(K,2)
!
!     Might also be necessary to place a limit on supersaturation.
!     Convection copes but cloud scheme response is less predictable.
!     Would need eg. LS_CLD_C to reduce supersaturation consistently.
!
          Do J=1, rows
            Do I=1, row_length
!
!           Ensure that input values of cloud fields lie within the
!           bounds of physical possibility (should do nothing).
!
              bulk_cf_n(I,J,K) = max(0., ( min(1., bulk_cf_n(I,J,K)) ))
              cf_liquid_n(I,J,K) =                                      &
     &                         max(0., ( min(1., cf_liquid_n(I,J,K)) ))
              cf_frozen_n(I,J,K) =                                      &
     &                         max(0., ( min(1., cf_frozen_n(I,J,K)) ))
            End Do  ! I
          End Do  ! J
!
        End Do  ! K
!
        If(n_qcx(1+wet_model_levels,1)  >   0)                          &
     &    WRITE(6,*) 'Qcf < 0 fixed by PC2 ',n_qcx(1+wet_model_levels,1)
        If(n_qcx(1+wet_model_levels,2)  >   0)                          &
     &    WRITE(6,*) 'Qcl < 0 fixed by PC2 ',n_qcx(1+wet_model_levels,2)
!
      End If  ! L_calc_dxek_if0
!
!---------------------------------------------------------------------
! Intercept values of physics increments before being updated by
! convection for optional output of convective increments
!---------------------------------------------------------------------
!
      If(l_qcl_incr_cinh) Then
!
        Allocate (qcl_incr_inhom_diag(row_length,rows,wet_model_levels))
!
! Hold input qcl increment (almost certainly zero)
        Do k=1,wet_model_levels
          Do j=1,rows
            Do i=1,row_length
              qcl_incr_inhom_diag(i,j,k) = qcl_inc(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k
!
      Endif                   ! on STASHflag
!
      If(l_qcf_incr_cinh) Then
!
        Allocate (qcf_incr_inhom_diag(row_length,rows,wet_model_levels))
!
! Hold input qcf increment (almost certainly zero)
        Do k=1,wet_model_levels
          Do j=1,rows
            Do i=1,row_length
              qcf_incr_inhom_diag(i,j,k) = qcf_inc(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k
!
      Endif                   ! on STASHflag
!
      If(l_bcf_incr_cinh) Then
!
        Allocate                                                        &
     & (bulk_cf_incr_inhom_diag(row_length,rows,wet_model_levels))
!
! Hold input bulk_cf increment (almost certainly zero)
        Do k=1,wet_model_levels
          Do j=1,rows
            Do i=1,row_length
              bulk_cf_incr_inhom_diag(i,j,k) = bulk_cf_inc(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k
!
      Endif                   ! on STASHflag
!
      If(l_cfl_incr_cinh .or. L_calc_dxek) Then
!
        Allocate                                                        &
     & (cf_liquid_incr_inhom_diag(row_length,rows,wet_model_levels))
!
! Hold input bulk_cf increment (almost certainly zero)
        Do k=1,wet_model_levels
          Do j=1,rows
            Do i=1,row_length
              cf_liquid_incr_inhom_diag(i,j,k) = cf_liquid_inc(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k
!
      Endif                   ! on STASHflag
!
      If(l_cff_incr_cinh) Then
!
        Allocate                                                        &
     & (cf_frozen_incr_inhom_diag(row_length,rows,wet_model_levels))
!
! Hold input bulk_cf increment (almost certainly zero)
        Do k=1,wet_model_levels
          Do j=1,rows
            Do i=1,row_length
              cf_frozen_incr_inhom_diag(i,j,k) = cf_frozen_inc(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k
!
      Endif                   ! on STASHflag
!
!
! ----------------------------------------------------------------------
! Section CNV.1 Set-up u and v on lorenz grid.
! Vertical interpolation to theta levels ONLY for 3C convection
! ----------------------------------------------------------------------

        If (l_mom) Then

! set number of levels to interpolate horizontally and vertically
          If (model_levels  >   wet_model_levels ) Then
! horizontally interpolate one more level than required and then
! vertically interpolate to all required levels.
            levels_h = wet_model_levels + 1
            levels_v = wet_model_levels
            L_copy_top = .false.
          Else
! horizontally interpolate all levels but can only vertically
! interpolate to levels - 1. Top level values are given by copying
! top level u and v values.
            levels_h = model_levels
            levels_v = model_levels - 1
            L_copy_top = .true.
          End If

      If(convection_option  <   3)then
         Do k = 1, levels_h
          Do j = 1 - off_y, rows + off_y
              Do i = 1 - off_x, row_length + off_x
                work_p_halo(i,j,k) = u(i,j,k)
              End Do
            End Do
          End Do

      Else         !convection_option 3
#if !defined(SCMA)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(                                           &
     &                       R_u, row_length, rows, model_levels,       &
     &                       off_x, off_y, fld_type_u, .true.  )
! DEPENDS ON: fill_external_halos
            CALL FILL_EXTERNAL_HALOS(R_u,row_length, rows,              &
     &                               model_levels,off_x,off_y)

        Do k = 1, levels_h
          Do j = 1 - off_y, rows + off_y
              Do i = 1 - off_x, row_length + off_x
                work_p_halo(i,j,k) = u(i,j,k) + R_u(i,j,k)
             End Do
            End Do
          End Do
#else
! SCM u on p grid
      Do k = 1,levels_h
        Do j = 1, rows
          Do i = 1, row_length
            work_p_halo(i,j,k) = u(i,j,k) + R_u(i,j,k)
          End Do ! i
        End Do ! j
      End Do ! k
#endif
      EndIf          !(convection_option  <   3)


! DEPENDS ON: u_to_p
           Call u_to_p (work_p_halo, row_length, rows, levels_h,        &
     &                 off_x, off_y, model_domain,                      &
     &                 at_extremity, u_p)

! set polar values of u to zero.
          If (model_domain  ==  mt_global .and.                         &
     &        at_extremity(PSouth) ) Then
            Do k = 1, levels_h
              Do i = 1, row_length
                u_p(i,1,k) = 0.
              End Do
            End Do
          End If
          If (model_domain  ==  mt_global .and.                         &
     &        at_extremity(PNorth) ) Then
            Do k = 1, levels_h
              Do i = 1, row_length
                u_p(i,rows,k) = 0.
              End Do
            End Do
          End If


      If(convection_option  <   3)then
         Do k = 1, levels_h
          Do j = 1 - off_y, n_rows + off_y
              Do i = 1 - off_x, row_length + off_x
                work_v_halo(i,j,k) = v(i,j,k)
              End Do
            End Do
          End Do

      Else         !convection_option 3

#if !defined(SCMA)
! DEPENDS ON: swap_bounds
            Call Swap_Bounds(                                           &
     &           R_v, row_length, n_rows, model_levels,                 &
     &            off_x, off_y, fld_type_v, .true.  )

! DEPENDS ON: fill_external_halos
            CALL FILL_EXTERNAL_HALOS(R_v,row_length, n_rows,            &
     &                               model_levels,off_x,off_y)

         Do k = 1, levels_h
          Do j = 1 - off_y, n_rows + off_y
              Do i = 1 - off_x, row_length + off_x
                work_v_halo(i,j,k) = v(i,j,k) + R_v(i,j,k)
              End Do
            End Do
          End Do
#else
! SCM v on p grid
      Do k = 1,levels_h
        Do j = 1, rows
          Do i = 1, row_length
            work_v_halo(i,j,k) = v(i,j,k) + R_v(i,j,k)
          End Do ! i
        End Do ! j
      End Do ! k
#endif
      EndIf          !(convection_option  <   3)


! DEPENDS ON: v_to_p
           Call v_to_p (work_v_halo, row_length, rows, n_rows,          &
     &                 levels_h, off_x, off_y, model_domain,            &
     &                 at_extremity, v_p)


! set polar values of v to zero.
          If (model_domain  ==  mt_global .and.                         &
     &        at_extremity(PSouth) ) Then
            Do k = 1, levels_h
              Do i = 1, row_length
                v_p(i,1,k) = 0.
              End Do
            End Do
          End If
          If (model_domain  ==  mt_global .and.                         &
     &        at_extremity(PNorth) ) Then
            Do k = 1, levels_h
              Do i = 1, row_length
                v_p(i,rows,k) = 0.
              End Do
            End Do
          End If

! need to set cumulus diagnosis to false at poles to ensure that
! not running cmt on points without valid winds and stress
!
        If (model_domain  ==  mt_global .and.                           &
     &        at_extremity(PSouth) ) Then
            Do i = 1, row_length
                cumulus(i,1) = .FALSE.
                l_shallow(i,1) = .FALSE.
                l_congestus(i,1) = .FALSE.
                l_congestus2(i,1) = .FALSE.
              End Do
        End If
        If (model_domain  ==  mt_global .and.                           &
     &        at_extremity(PNorth) ) Then
            Do i = 1, row_length
                cumulus(i,rows) = .FALSE.
                l_shallow(i,rows) = .FALSE.
                l_congestus(i,rows) = .FALSE.
                l_congestus2(i,rows) = .FALSE.
            End Do
          End If


        End If   !(l_mom)

!-----------------------------------------------------------------------
! Section 1.1  Setup total tracer variables
!-----------------------------------------------------------------------
      L_Tracer = ( ( tr_vars > 0 .OR. L_soot .OR. L_CO2_interactive .OR.&
     &               L_dust .OR. L_biomass .OR. L_sulpc_SO2 .OR.        &
     &               L_USE_CARIOLLE .OR. L_ocff .OR.                    &
     &               tr_ukca > 0 ) .AND. ( CycleNo == NumCycles ) )

! work with tracers only in final cycle
      If ( CycleNo == NumCycles ) Then
! Set up array dimensions for total tracer array (free + sulphur cycle
! tracers) so that convective transport of all tracers is done

        If ( (L_sulpc_SO2 .OR. L_soot .OR. L_CO2_interactive .OR.       &
     &        L_dust .OR. L_biomass .OR. L_USE_CARIOLLE .OR. L_ocff)    &
     &        .AND. (tr_vars > 0) .AND.                                 &
     &        (tr_levels /= model_levels ) )  Then ! exit
           Write(cmessage,*) 'Cannot call convect for tracer expts with'&
     &                      ,' tr_levels /= model_levels'
           ErrorStatus = -10
! DEPENDS ON: ereport
           Call Ereport(RoutineName, ErrorStatus, Cmessage )
       End If

       If (L_sulpc_SO2 .OR. L_soot .OR. L_dust .OR. L_biomass .OR.      &
     &     L_CO2_interactive .OR. L_USE_CARIOLLE .OR. L_ocff) Then
         ntra_lev = model_levels
         ntra_fld = 0        !Initialise to zero


         IF (L_DUST) THEN
           NTRA_FLD = NTRA_FLD + 6    !ADD 6 DUST SIZE CLASSES
         ENDIF

         If (L_SULPC_SO2) Then
           ntra_fld = ntra_fld + 4    !Add SO2 + 3 SO4 modes
           If (L_SULPC_NH3) Then
             ntra_fld = ntra_fld + 1  !Add NH3 field
           End If
           If (L_SULPC_DMS) Then
             ntra_fld = ntra_fld + 1  !Add DMS field
           End If
         End If

         If (L_SOOT) Then
           ntra_fld = ntra_fld + 3    !Add 3 modes of soot
         End If

         If (L_BIOMASS) Then
           ntra_fld = ntra_fld + 3    !Add 3 modes of biomass aerosol
         End If

         If (L_CO2_INTERACTIVE) Then
           ntra_fld = ntra_fld + 1    !Add CO2 field
         End If

         
         If (L_ocff) Then
           ntra_fld = ntra_fld + 3    !Add 3 modes of fossil-fuel OC
         End If
         
         
         If (L_USE_CARIOLLE) Then
           ntra_fld = ntra_fld + 1    !Add Cariolle Ozone tracer field
         End If

         If (tr_vars > 0) Then
           ntra_fld = ntra_fld + tr_vars
         End If

         If (tr_ukca > 0) Then
           ntra_fld = ntra_fld + tr_ukca
         End If

       Else
           If (tr_vars == 0 .AND. tr_ukca == 0) Then
             ntra_fld = 1             ! can't have 0 sized arrays
             ntra_lev = 1
           Else
             ntra_fld = tr_vars + tr_ukca
             ntra_lev = tr_levels
           End If
       End If

! Allocate the space in tot_tracer
       Allocate( tot_tracer(row_length, rows, ntra_lev, ntra_fld) )

! copy arrays into tot_tracer

       If (l_soot .OR. l_sulpc_SO2 .OR. l_CO2_interactive .OR.          &
     &     l_dust .OR. l_biomass .OR. L_USE_CARIOLLE .OR. L_ocff) Then
         ntra_tmp = 0

         IF (L_DUST) THEN
           NTRA_TMP = NTRA_TMP + 1
           TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP) = &
     &            DUST_DIV1(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS)

           NTRA_TMP = NTRA_TMP + 1
           TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP) = &
     &            DUST_DIV2(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS)

           NTRA_TMP = NTRA_TMP + 1
           TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP) = &
     &            DUST_DIV3(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS)

           NTRA_TMP = NTRA_TMP + 1
           TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP) = &
     &            DUST_DIV4(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS)

           NTRA_TMP = NTRA_TMP + 1
           TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP) = &
     &            DUST_DIV5(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS)

           NTRA_TMP = NTRA_TMP + 1
           TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP) = &
     &            DUST_DIV6(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS)

         ENDIF

         If (l_sulpc_SO2) Then
           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &            so2(1:row_length, 1:rows, 1:model_levels)

           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &     so4_aitken(1:row_length, 1:rows, 1:model_levels)

           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &       so4_accu(1:row_length, 1:rows, 1:model_levels)

           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &       so4_diss(1:row_length, 1:rows, 1:model_levels)

           If (l_sulpc_dms) Then
             ntra_tmp = ntra_tmp + 1
             tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)=&
     &              dms(1:row_length, 1:rows, 1:model_levels)
           End If

           If (l_sulpc_nh3) Then
             ntra_tmp = ntra_tmp + 1
             tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)=&
     &              nh3(1:row_length, 1:rows, 1:model_levels)
           End If
         End If

         If (l_soot) Then
           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &       soot_new(1:row_length, 1:rows, 1:model_levels)

           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &       soot_agd(1:row_length, 1:rows, 1:model_levels)

           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &       soot_cld(1:row_length, 1:rows, 1:model_levels)
         End If

         If (l_biomass) Then
           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &       bmass_new(1:row_length, 1:rows, 1:model_levels)

           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &       bmass_agd(1:row_length, 1:rows, 1:model_levels)

           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &       bmass_cld(1:row_length, 1:rows, 1:model_levels)
         End If

         If (l_CO2_interactive) Then
           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &            co2(1:row_length, 1:rows, 1:model_levels)
         End If
         
         If (l_ocff) Then
           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &       ocff_new(1:row_length, 1:rows, 1:model_levels)
 
           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &       ocff_agd(1:row_length, 1:rows, 1:model_levels)

           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &       ocff_cld(1:row_length, 1:rows, 1:model_levels)
         End If
         
         
         If (L_USE_CARIOLLE) Then
           ntra_tmp = ntra_tmp + 1
           tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
     &            OZONE_TRACER(1:row_length, 1:rows, 1:model_levels)
         End If

         If ( tr_vars > 0 ) Then
           tot_tracer(1:row_length, 1:rows, 1:tr_levels,                &
     &                ntra_tmp+1:ntra_tmp+tr_vars) =                    &
     &     free_tracers(1:row_length, 1:rows, 1:tr_levels,              &
     &                  1:tr_vars)
           ntra_tmp = ntra_tmp + tr_vars
         End If

         If ( tr_ukca > 0 ) Then
           tot_tracer(1:row_length, 1:rows, 1:tr_levels,                &
     &                ntra_tmp+1:ntra_tmp+tr_ukca) =                    &
     &     ukca_tracers(1:row_length, 1:rows, 1:tr_levels,              &
     &                  1:tr_ukca)
           ntra_tmp = ntra_tmp + tr_ukca
         End If

       Else    ! no soot, co2, sulphur, dust, biomass or OCFF
         If ( tr_vars > 0 ) Then
           tot_tracer  (1:row_length, 1:rows, 1:tr_levels, 1:tr_vars) = &
     &     free_tracers(1:row_length, 1:rows, 1:tr_levels, 1:tr_vars)
           ntra_tmp = tr_ukca
         End If

         If ( tr_ukca > 0 ) Then
           tot_tracer(1:row_length, 1:rows, 1:tr_levels,                &
     &                ntra_tmp+1:ntra_tmp+tr_ukca) =                    &
     &     ukca_tracers(1:row_length, 1:rows, 1:tr_levels,              &
     &                  1:tr_ukca)
           ntra_tmp = ntra_tmp + tr_ukca
         End If

         If ( tr_vars + tr_ukca == 0) Then ! make sure the value is set.
           tot_tracer(:,:,:,:) = 0.0
         End If
       End If !soot or co2, sulphur, dust, biomass, cariolle o3, or ocff

      Else
       ntra_fld = 1             ! can't have 0 sized arrays
       ntra_lev = 1
! Allocate dummy space in tot_tracer
       Allocate( tot_tracer(row_length, rows, ntra_fld, ntra_lev) )
      End If

! ----------------------------------------------------------------------
! Section CNV.2 Call Convection scheme.
! ----------------------------------------------------------------------
! Information required every call to convection

! Set number of levels to call convection for
      n_conv_levels = wet_model_levels

      If (l_mom) Then
! Limit convection calling levels to maxiumum of model_levels - 1
! This is because CMT increments to u and v exist on n_levels + 1
        If (n_conv_levels  >   model_levels - 1 )Then
          n_conv_levels = model_levels - 1
        End If
      End If

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
             r_rho_lev(i,j,k) = r_rho_levels(i,j,k)
          End Do
        End Do
      End Do

      Do k = 0, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            r_theta_lev(i,j,k) = r_theta_levels(i,j,k)
          End Do
        End Do
      End Do
!
! Sub-timestep scheme
!
        one_over_conv_calls = 1.0/(Num_Substeps*n_conv_calls*1.0)


      If ( Substep_Number  ==  1 ) Then

! 3D rain diagnostics

        If (flg_conv_rain_3d) Then
          Do k=1, wet_model_levels
            Do j=1, rows
              Do i=1, row_length
                conv_rain_3d(i,j,k) = 0.0
              End Do
            End Do
          End Do
        End If

! 3D snow diagnostics

        If (flg_conv_snow_3d) Then
          Do k=1, wet_model_levels
            Do j=1, rows
              Do i=1, row_length
                conv_snow_3d(i,j,k) = 0.0
              End Do
            End Do
          End Do
        End If

        Do j = 1, rows
          Do i = 1, row_length
            conv_rain(i,j) = 0.0
            conv_snow(i,j) = 0.0
            precip_deep(i,j) = 0.0
            precip_shall(i,j) = 0.0
            precip_mid(i,j) = 0.0
            precip_cong(i,j) = 0.0
            cape_out(i,j) = 0.0
            kterm_deep(i,j) = 0.0
            freeze_diag(i,j) = 0.0
            shallow_ind(i,j) = 0.0
            congestus_ind(i,j) = 0.0
            congestus_ind2(i,j) = 0.0
            mid_ind(i,j)   =0.0
! 5A diagnostics
            wstar_up_diag(i,j) = 0.0
            wstar_dn_diag(i,j) = 0.0
            mb1_diag(i,j) = 0.0
            mb2_diag(i,j) = 0.0
            cg_term(i,j) = 0.0
            cg_top(i,j) = 0.0
            cg_base(i,j) = 0.0

          End Do
        End Do
        Do j = 1, rows
          Do i = 1, row_length
            ntml_diag(i,j)  = float(ntml(i,j)*n_conv_calls)* one_over_conv_calls
            ntpar_diag(i,j) = float(ntpar(i,j)*n_conv_calls)*one_over_conv_calls
            wqt_cb(i,j)     = 0.0
            wthetal_cb(i,j) = 0.0
            wqt_inv(i,j)     = 0.0
            wthetal_inv(i,j) = 0.0
            sh_top(i,j)  = 0.0
            sh_base(i,j) = 0.0
          End Do
        End Do

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                up_flux(i,j,k)=0.0
                dwn_flux(i,j,k) = 0.0
                entrain_up(i,j,k) = 0.0
                entrain_dwn(i,j,k) = 0.0
                detrain_up(i,j,k) = 0.0
                detrain_dwn(i,j,k) = 0.0
                uw_dp(i,j,k) = 0.0
                vw_dp(i,j,k) = 0.0
                uw_shall(i,j,k) = 0.0
                vw_shall(i,j,k) = 0.0
              End Do
            End Do
          End Do

          If (flg_up_flx_half) then
            Do k=1,model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  up_flux_half(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif


          If (l_mom) Then
            If (l_up_incr_conv) Then
             Do k=1,model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dubydt_pout(i,j,k) = 0.0
                End Do
              End Do
             End Do
            Endif
            If (l_vp_incr_conv) Then
             Do k=1,model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dvbydt_pout(i,j,k) = 0.0
                End Do
              End Do
             End Do
            Endif
          Endif
! 5A diagnostics
          If (flg_wqt_flux) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  wqt_flux(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_wql_flux) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  wql_flux(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_wthetal_flux) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  wthetal_flux(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_wthetav_flux) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  wthetav_flux(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif

          If (flg_mf_deep) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mf_deep(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_mf_congest) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mf_congest(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_mf_shall) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mf_shall(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_mf_midlev) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mf_midlev(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dt_deep) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dt_deep(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dt_congest) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dt_congest(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dt_shall) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dt_shall(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dt_midlev) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dt_midlev(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dq_deep) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dq_deep(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dq_congest) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dq_congest(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dq_shall) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dq_shall(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dq_midlev) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dq_midlev(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_du_deep) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  du_deep(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_du_congest) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  du_congest(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_du_shall) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  du_shall(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_du_midlev) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  du_midlev(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dv_deep) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dv_deep(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dv_congest) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dv_congest(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dv_shall) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dv_shall(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif
          If (flg_dv_midlev) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dv_midlev(i,j,k) = 0.0
                End Do
              End Do
            End Do
          Endif


      Else    ! not first substep

! Set level diagnostics - required here because of substepping

        Do j = 1, rows
          Do i = 1, row_length
            ntml_diag(i,j)  = ntml_diag(i,j)                            &
                            + float(ntml(i,j)*n_conv_calls)*one_over_conv_calls
            ntpar_diag(i,j) = ntpar_diag(i,j)                           &
                            + float(ntpar(i,j)*n_conv_calls)*one_over_conv_calls
          End Do
        End Do

      Endif   ! end test on substep number


! At the last substep keep theta, q values as they
! will be used by the diagnostics subroutine
        If ( Substep_Number  ==  Num_Substeps ) Then
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                theta_diag(i,j,k) = theta_n(i,j,k)
              End Do
            End Do
          End Do
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_diag(i,j,k) = q_n(i,j,k)
              End Do
            End Do
          End Do
        Endif

        !---------------------------------------------------------------------
        ! Set up cloud decay lifetime
        !---------------------------------------------------------------------

        If (Rad_cloud_decay_opt /= rad_decay_off) Then
          Select Case (cld_life_opt)
            Case(cld_life_func_hgt)
              !---------------------------------------------------------------
              ! Make lifetime a function of cloud size by making it a function
              ! of height: 30 minutes for shallow, 2 hours for high cloud
              !---------------------------------------------------------------
              Do k=1, n_cca_levels
                Do j=1, rows
                  Do i=1, row_length
                    cld_life_3d(i,j,k) = 1800.0                               &
                                       + (fixed_cld_life*0.5 - 900.0)         &
                                       * (tanh((z_theta(i,j,k)/1000.0)-5.0)   &
                                       + 1.0)
                  End Do      ! i
                End Do      ! j
              End Do      ! k (n_cca_levels)

            Case(cld_life_constant)
              !---------------------------------------------------------------
              ! Set Cloud_lifetime to a constant
              !---------------------------------------------------------------
              Do k=1, n_cca_levels
                Do j=1, rows
                  Do i=1, row_length
                    cld_life_3d(i,j,k) = fixed_cld_life
                  End Do      ! i
                End Do      ! j
              End Do      ! k (n_cca_levels)
          End Select    ! cld_life_opt
        End If        !  Rad_cloud_decay_opt


        Do call_number = 1, n_conv_calls

          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                it_ccw          (i,j,k) = 0.0
                it_conv_rain_3d (i,j,k) = 0.0
                it_conv_snow_3d (i,j,k) = 0.0
              End Do
            End Do
          End Do

          Do k = 1, n_cca_levels
            Do j = 1, rows
              Do i = 1, row_length
                it_cca(i,j,k) = 0.0
              End Do
            End Do
          End Do
          Do j = 1, rows
            Do i = 1, row_length
              it_lcca(i,j) = 0.0
              it_lcbase(i,j) = 0
              it_lctop(i,j) = 0
              it_ccb(i,j) = 0
              it_cct(i,j) = 0
              it_conv_rain(i,j) = 0.0
              it_conv_snow(i,j) = 0.0
              it_precip_dp(i,j) = 0.0
              it_precip_sh(i,j) = 0.0
              it_precip_md(i,j) = 0.0
              it_precip_cg(i,j) = 0.0
              it_cclwp(i,j) = 0.0
              it_lcclwp(i,j) = 0.0
              it_cape_out(i,j) = 0.0
              it_kterm_deep(i,j) = 0
              it_mid_level(i,j) = .false.
              it_wstar_up(i,j) = 0.0
              it_wstar_dn(i,j) = 0.0
              it_mb1(i,j) = 0.0
              it_mb2(i,j) = 0.0
              it_cg_term(i,j) = 0
            End Do
          End Do

! NB: increments to t and q are added on inside routine.
! Increments to qCL, qCF, CFl, CFf are calculated but only added at the
! control level (see below).
! segment the call to convection, to reduce the memory required in
! the subroutine convect.

          If (a_convect_seg_size > 0) Then
            ! We are specifying size of segments rather than number
            step = a_convect_seg_size
            a_convect_segments = Ceiling( Real(rows*row_length)/        &
                                          Real(step) )
          Else
            step = rows*row_length/a_convect_segments
          End If

          ! Allocate space for segmentation arrays
          Allocate(    sp_local( a_convect_segments ) )
          Allocate(    fp_local( a_convect_segments ) )
          Allocate(         fyp( a_convect_segments ) )
          Allocate(         fxp( a_convect_segments ) )
          Allocate(   n_cumulus( a_convect_segments ) )
          Allocate(      n_deep( a_convect_segments ) )
          Allocate(   n_shallow( a_convect_segments ) )
          Allocate( n_congestus( a_convect_segments ) )

          first_point = 1
          Do i = 1, a_convect_segments
            seg_points = step
            If (i == a_convect_segments) Then
              seg_points = rows*row_length-step*(a_convect_segments-1)
            End If

            fp_local(i) = first_point
            fyp(i) = (first_point-1)/row_length + 1
            fxp(i) = first_point-(fyp(i)-1)*row_length
            sp_local(i) = seg_points
            first_point = first_point+step
          End Do
!
! WORK OUT NUMBER OF CUMULUS POINTS IN EACH CONVECTION SEGMENT
! Also calculate number of shallow and deep points in each convection
! segment
!
          ii=1
          Do j = 1, rows
            Do i = 1, row_length
              cumulus_1d(ii) = cumulus(i,j)
              l_shallow_1d(ii) = l_shallow(i,j)
              l_congestus_1d(ii) = l_congestus(i,j)
              ii=ii+1
            Enddo
          Enddo

! Nec compilier directive - small loop don't try to vectorise as often
!  a_convective_segments may be 1.
!CDIR NOVECTOR
          Do i=1,a_convect_segments
            n_cumulus(i)=0
            n_deep(i)   =0
            n_shallow(i)   = 0
            n_congestus(i) = 0

            Do j=fp_local(i),fp_local(i)+sp_local(i)-1
              If(cumulus_1d(j)) then
                n_cumulus(i)=n_cumulus(i)+1

                If (iconv_deep >  0.and..not. l_shallow_1d(j).and.      &
     &                      .not. l_congestus_1d(j)) then

                  n_deep(i) = n_deep(i)+1
                End If
                If (iconv_shallow >  0.and.l_shallow_1d(j)) then
                  n_shallow(i) = n_shallow(i)+1
                End If
                If (iconv_congestus >  0) then
                  If (l_congestus_1d(j)) then
                    n_congestus(i) = n_congestus(i)+1
                  End If
                Else
                  n_congestus(i) = 1     ! may be required for dim
                End If
              Endif
            End do

          End do

          Do i = 1, a_convect_segments

! DEPENDS ON: glue_conv
              Call glue_conv(rows*row_length, sp_local(i),              &
     &        n_conv_levels, bl_levels,                                 &
     &        Theta_n(fxp(i),fyp(i),1), q_n(fxp(i),fyp(i),1),           &
     &        qcl_n(fxp(i),fyp(i),1), qcf_n(fxp(i),fyp(i),1),           &
     &        cf_liquid_n(fxp(i),fyp(i),1),cf_frozen_n(fxp(i),fyp(i),1),&
     &        bulk_cf_n(fxp(i),fyp(i),1),                               &
     &        p_star(fxp(i),fyp(i)),land_sea_mask(fxp(i),fyp(i)),       &
     &        u_p(fxp(i),fyp(i),1), v_p(fxp(i),fyp(i),1),               &
     &        tot_tracer(fxp(i),fyp(i),1,1), dthbydt(fxp(i),fyp(i),1),  &
     &        dqbydt(fxp(i),fyp(i),1), dqclbydt(fxp(i),fyp(i),1),       &
     &        dqcfbydt(fxp(i),fyp(i),1), dcflbydt(fxp(i),fyp(i),1),     &
     &        dcffbydt(fxp(i),fyp(i),1), dbcfbydt(fxp(i),fyp(i),1),     &
     &        dubydt_p(fxp(i),fyp(i),1),                                &
     &        dvbydt_p(fxp(i),fyp(i),1),                                &
     &        it_conv_rain(fxp(i),fyp(i)),                              &
     &        it_conv_snow(fxp(i),fyp(i)),                              &
     &        it_conv_rain_3d(fxp(i),fyp(i),1),                         &
     &        it_conv_snow_3d(fxp(i),fyp(i),1),                         &
     &        it_cca(fxp(i),fyp(i),1),it_ccb(fxp(i),fyp(i)),            &
     &        it_cct(fxp(i),fyp(i)),                                    &
     &        it_cclwp(fxp(i),fyp(i)), it_ccw(fxp(i),fyp(i),1),         &
     &        it_lcbase(fxp(i),fyp(i)), it_lctop(fxp(i),fyp(i)),        &
     &        it_lcca(fxp(i),fyp(i)),freeze_lev(fxp(i),fyp(i)),         &
     &        it_mid_level(fxp(i),fyp(i)),it_kterm_deep(fxp(i),fyp(i)), &
     &        it_precip_dp(fxp(i),fyp(i)),it_precip_sh(fxp(i),fyp(i)),  &
     &        it_precip_md(fxp(i),fyp(i)),it_precip_cg(fxp(i),fyp(i)),  &
     &        it_wstar_dn(fxp(i),fyp(i)),                               &
     &        it_wstar_up(fxp(i),fyp(i)),it_mb1(fxp(i),fyp(i)),         &
     &        it_mb2(fxp(i),fyp(i)),it_cg_term(fxp(i),fyp(i)),          &
     &        n_cumulus(i),                                             &
     &        uw0_p(fxp(i),fyp(i)),vw0_p(fxp(i),fyp(i)),                &
     &        w_max(fxp(i),fyp(i)), zlcl(fxp(i),fyp(i)),                &
     &        zlcl_uv(fxp(i),fyp(i)),ztop(fxp(i),fyp(i)),               &
     &        it_lcclwp(fxp(i),fyp(i)),it_cape_out(fxp(i),fyp(i)),      &
     &        n_deep(i),n_congestus(i),n_shallow(i),                    &
     &        r_rho_lev(fxp(i),fyp(i),1),r_theta_lev(fxp(i),fyp(i),0)   &
     &       ,rho_only(fxp(i),fyp(i),1),rho_theta(fxp(i),fyp(i),1),     &
     &        exner_layer_boundaries(fxp(i),fyp(i),0),                  &
     &        exner_layer_centres(fxp(i),fyp(i),0),                     &
     &        p_layer_boundaries(fxp(i),fyp(i),0),                      &
     &        p_layer_centres(fxp(i),fyp(i),0),                         &
     &        z_theta(fxp(i),fyp(i),1),z_rho(fxp(i),fyp(i),1),          &
     &        timestep_conv, t1_sd(fxp(i),fyp(i)), q1_sd(fxp(i),fyp(i)),&
     &        NTML(fxp(i),fyp(i)),NTPAR(fxp(i),fyp(i)),                 &
     &        NBDSC(fxp(i),fyp(i)),NTDSC(fxp(i),fyp(i)),                &
     &        L_SHALLOW(fxp(i),fyp(i)),l_pc2_diag_sh_pts(fxp(i),fyp(i)),&
     &        l_congestus(fxp(i),fyp(i)),                               &
     &        CUMULUS(fxp(i),fyp(i)),                                   &
     &        WSTAR(fxp(i),fyp(i)),WTHVS(fxp(i),fyp(i)),                &
     &        DELTHVU(fxp(i),fyp(i)),ql_ad(fxp(i),fyp(i)),              &
     &        ftl(fxp(i),fyp(i)),fqt(fxp(i),fyp(i)),                    &
     &        L_tracer, ntra_fld, ntra_lev, n_cca_levels,               &
     &        l_mixing_ratio, L_calc_dxek, L_q_interact,                &
     &        it_up_flux_half(fxp(i),fyp(i),1),flg_up_flx_half,         &
     &        it_up_flux(fxp(i),fyp(i),1),flg_up_flx,                   &
     &        it_dwn_flux(fxp(i),fyp(i),1),flg_dwn_flx,                 &
     &        it_entrain_up(fxp(i),fyp(i),1),flg_entr_up,               &
     &        it_detrain_up(fxp(i),fyp(i),1),flg_detr_up,               &
     &        it_entrain_dwn(fxp(i),fyp(i),1),flg_entr_dwn,             &
     &        it_detrain_dwn(fxp(i),fyp(i),1),flg_detr_dwn              &
     &       ,it_uw_dp(fxp(i),fyp(i),1),flg_uw_dp                       &
     &       ,it_vw_dp(fxp(i),fyp(i),1),flg_vw_dp                       &
     &       ,it_uw_shall(fxp(i),fyp(i),1),flg_uw_shall                 &
     &       ,it_vw_shall(fxp(i),fyp(i),1),flg_vw_shall                 &
     &,  flg_wqt_flux,flg_wthetal_flux,flg_wthetav_flux,flg_wql_flux    &
     &,      flg_mf_deep,flg_mf_congest,flg_mf_shall,flg_mf_midlev      &
     &,      flg_dt_deep,flg_dt_congest,flg_dt_shall,flg_dt_midlev      &
     &,      flg_dq_deep,flg_dq_congest,flg_dq_shall,flg_dq_midlev      &
     &,      flg_du_deep,flg_du_congest,flg_du_shall,flg_du_midlev      &
     &,      flg_dv_deep,flg_dv_congest,flg_dv_shall,flg_dv_midlev      &
     &, it_wqt_flux(fxp(i),fyp(i),1),it_wthetal_flux(fxp(i),fyp(i),1)   &
     &, it_wthetav_flux(fxp(i),fyp(i),1),it_wql_flux(fxp(i),fyp(i),1)   &
     &, it_mf_deep(fxp(i),fyp(i),1),it_mf_congest(fxp(i),fyp(i),1)      &
     &, it_mf_shall(fxp(i),fyp(i),1),it_mf_midlev(fxp(i),fyp(i),1)      &
     &, it_dt_deep(fxp(i),fyp(i),1),it_dt_congest(fxp(i),fyp(i),1)      &
     &, it_dt_shall(fxp(i),fyp(i),1),it_dt_midlev(fxp(i),fyp(i),1)      &
     &, it_dq_deep(fxp(i),fyp(i),1),it_dq_congest(fxp(i),fyp(i),1)      &
     &, it_dq_shall(fxp(i),fyp(i),1),it_dq_midlev(fxp(i),fyp(i),1)      &
     &, it_du_deep(fxp(i),fyp(i),1),it_du_congest(fxp(i),fyp(i),1)      &
     &, it_du_shall(fxp(i),fyp(i),1),it_du_midlev(fxp(i),fyp(i),1)      &
     &, it_dv_deep(fxp(i),fyp(i),1),it_dv_congest(fxp(i),fyp(i),1)      &
     &, it_dv_shall(fxp(i),fyp(i),1),it_dv_midlev(fxp(i),fyp(i),1)      &
     &         )

          End Do  ! loop over number of segments

          ! Deallocate segmentation variables
          Deallocate( sp_local )
          Deallocate( fp_local )
          Deallocate( fyp )
          Deallocate( fxp )
          Deallocate( n_cumulus )
          Deallocate( n_deep )
          Deallocate( n_shallow )
          Deallocate( n_congestus )

         loop_number = call_number * Substep_Number

         If (L_fix_udfactor) Then 
!          ! Use ud_factor of 1 in convection to return CCW and cclwp
!          ! unscaled.  Then the true ud_factor can be applied here to
!          ! the whole profile, rather than just to precipitating levels
           Do j = 1, rows
             Do i = 1, row_length
                it_cclwp(i,j)  = ud_factor * it_cclwp(i,j)
                it_lcclwp(i,j) = ud_factor * it_lcclwp(i,j)
             End Do
           End Do
           Do k = 1, wet_model_levels
             Do j = 1, rows
               Do i = 1, row_length
                  it_ccw(i,j,k) = ud_factor * it_ccw(i,j,k)
               End Do
             End Do
           End Do

         End If


         If (l_convcld_hadgem1) Then

           ! Use values from first call only
           If (loop_number == 1) Then
             Do j=1, rows
               Do i=1, row_length
                 lcca   (i,j) = it_lcca   (i,j)
                 lcbase (i,j) = it_lcbase (i,j)
                 lctop  (i,j) = it_lctop  (i,j)
                 cclwp  (i,j) = it_cclwp  (i,j)
                 lcclwp (i,j) = it_lcclwp (i,j)
                 ccb    (i,j) = it_ccb    (i,j)
                 cct    (i,j) = it_cct    (i,j)
               End Do
             End Do

             Do k=1, wet_model_levels
               Do j=1, rows
                 Do i=1, row_length
                   ccw(i,j,k) = it_ccw(i,j,k)
                 End Do
               End Do
             End Do
           End If      ! Test loop_number

         Else     ! Test l_convcld_hadgem1

           ! Original code

           If (loop_number == 1) Then

             !----------------------------------------------------------------
             ! Code in this loop is performed on 1ST call to convection
             ! DURING A MODEL PHYSICS TIMESTEP
             !----------------------------------------------------------------

             Do j=1, rows
               Do i=1, row_length
                 lcca  (i,j) = it_lcca  (i,j)
                 lcbase(i,j) = it_lcbase(i,j)
                 lctop (i,j) = it_lctop (i,j)
                 lcclwp(i,j) = it_lcclwp(i,j) * one_over_conv_calls

                 ccb   (i,j) = it_ccb   (i,j)
                 cct   (i,j) = it_cct   (i,j)
                 cclwp (i,j) = it_cclwp (i,j) * one_over_conv_calls
               End Do      ! i
             End Do      ! j

             Do k=1, wet_model_levels
               Do j=1, rows
                 Do i=1, row_length
                   ccw(i,j,k) = it_ccw(i,j,k) * one_over_conv_calls
                 End Do      ! i
               End Do      ! j
             End Do      ! k

             Do k=1, n_cca_levels
               Do j=1, rows
                 Do i=1, row_length
                   cca(i,j,k) = it_cca(i,j,k) * one_over_conv_calls
                 End Do      ! i
               End Do      ! j
             End Do      ! k

             If (lcv_ccrad) Then
               Do k=1, n_cca_levels
                 Do j=1, rows
                   Do i=1, row_length
                     cca(i,j,k) = MIN(cca(i,j,k), 1.0)
                   End Do      ! i
                 End Do      ! j
               End Do      ! k
             End If      ! lcv_ccrad

           Else       ! loop_number /= 1

             !----------------------------------------------------------------
             ! Code in this loop is performed on every call to convection
             ! EXCEPT THE 1st in the MODEL PHYSICS TIMESTEP
             !----------------------------------------------------------------

             Do j=1, rows
               Do i=1, row_length

                 lctop (i,j) = MAX(lctop(i,j), it_lctop(i,j))

                 lcclwp(i,j) = lcclwp(i,j)                                    &
                             + (it_lcclwp(i,j)*one_over_conv_calls)

                 lcca  (i,j) = lcca(i,j)                                      &
                             + (it_lcca(i,j)*one_over_conv_calls)

                 If (lcbase(i,j) > 0 .AND. it_lcbase(i,j) > 0) Then
                   lcbase(i,j) = MIN(lcbase(i,j), it_lcbase(i,j))
                 Else
                   lcbase(i,j) = MAX(lcbase(i,j), it_lcbase(i,j))
                 End If

               End Do      ! i (row_length)
             End Do      ! j (rows)

             Do k=1, wet_model_levels
               Do j=1, rows
                 Do i=1, row_length
                   ccw(i,j,k) = ccw(i,j,k)                                    &
                              + (it_ccw(i,j,k)*one_over_conv_calls)
                 End Do      ! i
               End Do      ! j
             End Do      ! k


             If ( Rad_cloud_decay_opt /= rad_decay_off .OR.                   &
                  n_conv_calls > 1 ) Then

               Do j=1, rows
                 Do i=1, row_length

                   cct   (i,j) = MAX(cct(i,j), it_cct(i,j))
                   cclwp (i,j) = cclwp (i,j)                                  &
                               + (it_cclwp (i,j)*one_over_conv_calls)

                   If (ccb(i,j) > 0 .AND. it_ccb(i,j) > 0) Then
                     ccb(i,j) = MIN(ccb(i,j), it_ccb(i,j))
                   Else
                     ccb(i,j) = MAX(ccb(i,j), it_ccb(i,j))
                   End If

                 End Do      ! i (row_length)
               End Do      ! j (rows)

               Do k=1, n_cca_levels
                 Do j=1, rows
                   Do i=1, row_length
                     cca(i,j,k) = cca(i,j,k)                                  &
                                + (it_cca(i,j,k)*one_over_conv_calls)
                   End Do      ! i 
                 End Do      ! j
               End Do      ! k

             Else 

               ! Radiative cloud decay disabled. .AND.
               ! 1 call to convection per BL step.

               ! ccb_rad holds ccb values from previous steps

               Do j=1, rows
                 Do i=1, row_length

                   cct   (i,j) = MAX(cct_rad(i,j), it_cct(i,j))
                   cclwp (i,j) = cclwp_rad(i,j)                               &
                               + (it_cclwp (i,j)*one_over_conv_calls)

                   If (ccb_rad(i,j) > 0 .AND. it_ccb(i,j) > 0) Then
                     ccb(i,j) = MIN(ccb_rad(i,j), it_ccb(i,j))
                   Else
                     ccb(i,j) = MAX(ccb_rad(i,j), it_ccb(i,j))
                   End If

                 End Do      ! i (row_length)
               End Do      ! j (rows)

               Do k=1, n_cca_levels
                 Do j=1, rows
                   Do i=1, row_length
                     cca(i,j,k) = cca_rad(i,j,k)                              &
                                + (it_cca(i,j,k)*one_over_conv_calls)
                   End Do      ! i 
                 End Do      ! j
               End Do      ! k

             End If  ! Test on rad_cloud_decay_opt / n_conv_calls

             If (lcv_ccrad) Then
               Do k=1, n_cca_levels
                 Do j=1, rows
                   Do i=1, row_length
                     cca(i,j,k) = MIN(cca(i,j,k), 1.0)
                   End Do      ! i
                 End Do      ! j
               End Do      ! k
             End If      ! lcv_ccrad

           End If      ! Test on Loop Number
         End If      ! Test on Hadgem1


         If (rad_cloud_decay_opt == rad_decay_conv_substep) Then
           !------------------------------------------------------------------
           ! Rad_cloud_decay_opt 2, to use CCRad cloud decay on convection    
           ! substeps                                                         
           !------------------------------------------------------------------

           !------------------------------------------------------------------
           ! Convection Cloud Decay Code - Affects convective cloud seen by
           ! radiation scheme. Lock/Wong
           !------------------------------------------------------------------
           ! USE CODE TO ALLOW DECAY OF CONVECTIVE CLOUD FOR RADIATIVE IMPACTS
           ! ON CONVECTION SUBSTEP.
           !------------------------------------------------------------------
           ! NOTE: CCA_rad, CCW_rad profiles are the resultant profile seen by
           !       radiation, they may contain elements of cloud from previous
           !       timesteps/conv_timesteps.
           !
           !       With decayed cloud the lctop and ccb become meaningless.
           !
           !       The CCB_RAD and CCT_RAD exiting from this routine will be
           !       the EXTENTS of the resultant profile. LCBASE, LCTOP, CCB
           !       and CCT sent to the stash will continue to refer to
           !       convection on a given timestep.
           !------------------------------------------------------------------

           !------------------------------------------------------------------
           ! Reinitialise all parameters sent to radiation based on resultant
           ! cloud profiles sent to radiation. lcbase_rad, lctop_rad, cct_rad, 
           ! ccb_rad are all calculated from resultant cca profile.
           !------------------------------------------------------------------

           Do j=1, rows
             Do i=1, row_length
               lcbase_rad (i,j) = 0
               lctop_rad  (i,j) = 0
               ccb_rad    (i,j) = 0
               cct_rad    (i,j) = 0
               cclwp_rad  (i,j) = 0.0 
             End Do      ! i
           End Do      ! j


           !------------------------------------------------------------------
           ! DECAY OF PREVIOUS CCA PROFILE:
           ! Decays across convection substep, model levels with convective
           ! cloud on THIS substep overwrite the decaying profile on THOSE
           ! models levels. Calculate new cloud extents based on resultant
           ! profile.
           !------------------------------------------------------------------

           Do k=2, n_cca_levels
             Do j=1, rows
               Do i=1, row_length

                 !------------------------------------------------------------
                 ! Determine the higher value of convective cloud amount, CCA
                 ! from call to convection OR decayed CCA from previous
                 ! timesteps
                 !------------------------------------------------------------
                 cca_rad(i,j,k) = MAX( cca_rad(i,j,k)                         &
                                      *(1.0-timestep_conv/Cld_life_3d(i,j,k)) &
                                      ,it_cca(i,j,k))

                 cca_rad(i,j,k) = MIN( cca_rad(i,j,k),1.0 )



                 !------------------------------------------------------------
                 ! Apply threshold resultant profile for minimum CCA
                 !------------------------------------------------------------
                 If (cca_rad(i,j,k) < CCA_min ) Then 
                   cca_rad(i,j,k) = 0.0
                 End If



                 !------------------------------------------------------------
                 ! Check for resultant cloud bases for radiation
                 !------------------------------------------------------------
                 If (cca_rad(i,j,k)   >  0.0 .AND.                            &
                     cca_rad(i,j,k-1) == 0.0) Then

                   ccb_rad(i,j) = k

                   If (lcbase_rad(i,j) == 0) Then
                     lcbase_rad(i,j) = k
                   End If

                 End If


                 !------------------------------------------------------------
                 ! Check for resultant cloud tops for radiation
                 !------------------------------------------------------------
                 If (cca_rad(i,j,k)   == 0.0 .AND.                            &
                     cca_rad(i,j,k-1) >  0.0) Then

                   cct_rad(i,j) = k-1

                   If (lctop_rad(i,j) == 0) Then
                     lctop_rad(i,j) = k-1
                   End If
                 End If

               End Do      ! i,  (row_length)
             End Do      ! j,  (rows)
           End Do      ! k,  (n_cca_levels)


           !------------------------------------------------------------------
           ! DECAY OF PREVIOUS CCW PROFILE:
           ! Decay across convection substep, model levels with convective
           ! cloud on THIS substep overwrite the the decaying profile on THOSE
           ! models levels.
           !------------------------------------------------------------------
            
           Do k=1, wet_model_levels
             Do j=1, rows
               Do i=1, row_length
                  
                 If (cca_rad(i,j,k) > 0.0) Then
                     ccw_rad(i,j,k) = MAX( ccw_rad(i,j,k)                     &
                                    *(1.0-timestep_conv/cld_life_3d(i,j,k))   &
                                    , it_ccw(i,j,k) )

                     cclwp_rad(i,j) = cclwp_rad(i,j) + ccw_rad(i,j,k)         &
                                    * (  p_layer_boundaries(i,j,k-1)          &
                                       - p_layer_boundaries(i,j,k) )/g
                 Else
                   ccw_rad(i,j,k) = 0.0
                 End If

               End Do      ! i, (row_length)
             End Do      ! j, (rows)
           End Do      ! k, (wet_model_levels)

         End If      ! rad_cloud_decay_opt =
                     !            rad_decay_conv_substep (lcv_ccrad)


! All loops

          If (flg_conv_rain_3d) Then
            Do k=1, wet_model_levels
              Do j=1, rows
                Do i=1, row_length
                conv_rain_3d(i,j,k) =    conv_rain_3d(i,j,k)                  &
                                    + it_conv_rain_3d(i,j,k)                  &
                                    * one_over_conv_calls
                End Do      ! i (row_length)
              End Do      ! j (rows)
            End Do      ! k (wet_model_levels)
          End If

          If (flg_conv_snow_3d) Then
            Do k=1, wet_model_levels
              Do j=1, rows
                Do i=1, row_length
                  conv_snow_3d(i,j,k) =    conv_snow_3d(i,j,k)                &
                                      + it_conv_snow_3d(i,j,k)                &
                                      * one_over_conv_calls
                End Do      ! i (row_length)
              End Do      ! j (rows)
            End Do      ! k (wet_model_levels)
          End If

          Do j = 1, rows
            Do i = 1, row_length

              conv_rain(i,j) = conv_rain(i,j) + it_conv_rain(i,j)       &
     &                                           * one_over_conv_calls
              conv_snow(i,j) = conv_snow(i,j) + it_conv_snow(i,j)       &
     &                                           * one_over_conv_calls
              precip_deep(i,j) = precip_deep(i,j) + it_precip_dp(i,j)   &
     &                                           * one_over_conv_calls
              precip_shall(i,j)= precip_shall(i,j) + it_precip_sh(i,j)  &
     &                                           * one_over_conv_calls
              precip_mid(i,j) = precip_mid(i,j) + it_precip_md(i,j)     &
     &                                           * one_over_conv_calls
              precip_cong(i,j) = precip_cong(i,j) + it_precip_cg(i,j)   &
     &                                           * one_over_conv_calls
              cape_out(i,j)   = cape_out(i,j)                           &
     &                           +it_cape_out(i,j)*one_over_conv_calls
              if (it_mid_level(i,j)) then
                mid_ind(i,j) = mid_ind(i,j)+one_over_conv_calls
              endif
              if (l_shallow(i,j)) then
                shallow_ind(i,j) = shallow_ind(i,j)                     &
     &                                         +one_over_conv_calls
                 k=ntpar(i,j)+1
                 sh_top(i,j)=sh_top(i,j)+ z_rho(i,j,k)                  &
     &                                         *one_over_conv_calls
                 k=ntml(i,j)+1
                 sh_base(i,j)=sh_base(i,j)+ z_rho(i,j,k)                &
     &                                         *one_over_conv_calls
              endif
              if (l_congestus(i,j)) then
                congestus_ind(i,j) = congestus_ind(i,j)                 &
     &                                         +one_over_conv_calls
                 k=it_cg_term(i,j)+1
                 cg_top(i,j)=cg_top(i,j)+ z_rho(i,j,k)                  &
     &                                         *one_over_conv_calls
                 k=ntml(i,j)+1
                 cg_base(i,j)=cg_base(i,j)+ z_rho(i,j,k)                &
     &                                         *one_over_conv_calls

              endif
              if (l_congestus2(i,j)) then
                congestus_ind2(i,j) = congestus_ind2(i,j)               &
     &                                         +one_over_conv_calls
              endif
              kterm_deep(i,j) = kterm_deep(i,j)+                        &
     &               float(it_kterm_deep(i,j)) *one_over_conv_calls
              freeze_diag(i,j) = freeze_diag(i,j)+                      &
     &               float(freeze_lev(i,j)) *one_over_conv_calls
! 5A only
              wstar_up_diag(i,j) = wstar_up_diag(i,j)                   &
     &                           +it_wstar_up(i,j)*one_over_conv_calls
              wstar_dn_diag(i,j) = wstar_dn_diag(i,j)                   &
     &                           +it_wstar_dn(i,j)*one_over_conv_calls
              mb1_diag(i,j) = mb1_diag(i,j)                             &
     &                           +it_mb1(i,j)*one_over_conv_calls
              mb2_diag(i,j) = mb2_diag(i,j)                             &
     &                           +it_mb2(i,j)*one_over_conv_calls
              cg_term(i,j) = cg_term(i,j)+                              &
     &               float(it_cg_term(i,j)) *one_over_conv_calls

            End Do
          End Do
          If (flg_mf_deep) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mf_deep(i,j,k) = mf_deep(i,j,k)                       &
     &                           +it_mf_deep(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_mf_congest) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mf_congest(i,j,k) = mf_congest(i,j,k)                 &
     &                        +it_mf_congest(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_mf_shall) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mf_shall(i,j,k) = mf_shall(i,j,k)                     &
     &                         +it_mf_shall(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_mf_midlev) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mf_midlev(i,j,k) = mf_midlev(i,j,k)                   &
     &                       +it_mf_midlev(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dt_deep) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dt_deep(i,j,k) = dt_deep(i,j,k)                       &
     &                       +it_dt_deep(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dt_congest) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dt_congest(i,j,k) = dt_congest(i,j,k)                 &
     &                     +it_dt_congest(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dt_shall) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dt_shall(i,j,k) = dt_shall(i,j,k)                     &
     &                     +it_dt_shall(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dt_midlev) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dt_midlev(i,j,k) = dt_midlev(i,j,k)                   &
     &                     +it_dt_midlev(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dq_deep) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dq_deep(i,j,k) = dq_deep(i,j,k)                       &
     &                     +it_dq_deep(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dq_congest) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dq_congest(i,j,k) = dq_congest(i,j,k)                 &
     &                    +it_dq_congest(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dq_shall) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dq_shall(i,j,k) = dq_shall(i,j,k)                     &
     &                     +it_dq_shall(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dq_midlev) then
            Do k=1,n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dq_midlev(i,j,k) = dq_midlev(i,j,k)                   &
     &                     +it_dq_midlev(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_du_deep) then
            Do k=1,n_conv_levels+1
              Do j = 1, rows
                Do i = 1, row_length
                  du_deep(i,j,k) =  du_deep(i,j,k)                      &
     &                        +it_du_deep(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_du_congest) then
            Do k=1,n_conv_levels+1
              Do j = 1, rows
                Do i = 1, row_length
                  du_congest(i,j,k) = du_congest(i,j,k)                 &
     &                     +it_du_congest(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_du_shall) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  du_shall(i,j,k) = du_shall(i,j,k)                     &
     &                     +it_du_shall(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_du_midlev) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  du_midlev(i,j,k) = du_midlev(i,j,k)                   &
     &                     +it_du_midlev(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dv_deep) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dv_deep(i,j,k) = dv_deep(i,j,k)                       &
     &                      +it_dv_deep(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dv_congest) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dv_congest(i,j,k) = dv_congest(i,j,k)                 &
     &                     +it_dv_congest(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dv_shall) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dv_shall(i,j,k) = dv_shall(i,j,k)                     &
     &                      +it_dv_shall(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif
          If (flg_dv_midlev) then
            Do k=1,wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dv_midlev(i,j,k) = dv_midlev(i,j,k)                   &
     &                    +it_dv_midlev(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          Endif


! ----------------------------------------------------------------------
! Section CNV.3 Add on theta and q increments, qCL and qCF increments.
! ----------------------------------------------------------------------

! add on increments to theta and q for next convection call
          Do k = 1, n_conv_levels
            Do j = 1, rows
              Do i = 1, row_length
                theta_inc(i,j,k) = theta_inc(i,j,k)                     &
     &                           + dthbydt(i,j,k) * timestep_conv
                q_inc(i,j,k) = q_inc(i,j,k)                             &
     &                          + dqbydt(i,j,k) * timestep_conv
                qcl_inc(i,j,k) = qcl_inc(i,j,k) +                       &
     &                         (dqclbydt(i,j,k) * timestep_conv)
                qcf_inc(i,j,k) = qcf_inc(i,j,k) +                       &
     &                         (dqcfbydt(i,j,k) * timestep_conv)
                cf_liquid_inc(i,j,k) = cf_liquid_inc(i,j,k) +           &
     &                         (dcflbydt(i,j,k) * timestep_conv)
                cf_frozen_inc(i,j,k) = cf_frozen_inc(i,j,k) +           &
     &                         (dcffbydt(i,j,k) * timestep_conv)
                bulk_cf_inc(i,j,k)   = bulk_cf_inc(i,j,k) +             &
     &                         (dbcfbydt(i,j,k) * timestep_conv)
                theta_n(i,j,k) = theta_n(i,j,k)                         &
     &                           + dthbydt(i,j,k) * timestep_conv
                q_n(i,j,k) = q_n(i,j,k)                                 &
     &                          + dqbydt(i,j,k) * timestep_conv
                qcl_n(i,j,k) = qcl_n(i,j,k) +                           &
     &                         (dqclbydt(i,j,k) * timestep_conv)
                qcf_n(i,j,k) = qcf_n(i,j,k) +                           &
     &                         (dqcfbydt(i,j,k) * timestep_conv)
                cf_liquid_n(i,j,k) = cf_liquid_n(i,j,k) +               &
     &                         (dcflbydt(i,j,k) * timestep_conv)
                cf_frozen_n(i,j,k) = cf_frozen_n(i,j,k) +               &
     &                         (dcffbydt(i,j,k) * timestep_conv)
                bulk_cf_n(i,j,k)   = bulk_cf_n(i,j,k) +                 &
     &                         (dbcfbydt(i,j,k) * timestep_conv)
              End Do
            End Do
          End Do
!

! Update T, q diagnostic incerements
          If ( l_T_incr_conv ) Then
            Do k=1,n_conv_levels
              Do j=1,rows
                Do i=1,row_length
                  T_incr_diagnostic(i,j,k) = T_incr_diagnostic(i,j,k)   &
     &                           + dthbydt(i,j,k) * timestep_conv       &
     &                                     * exner_theta_levels(i,j,k)
                Enddo ! i
              Enddo ! j
            Enddo ! k
          Endif                   ! on STASHflag

          If ( l_q_incr_conv ) Then
            Do k=1,n_conv_levels
              Do j=1,rows
                Do i=1,row_length
                  q_incr_diagnostic(i,j,k) = q_incr_diagnostic(i,j,k)   &
     &                          + dqbydt(i,j,k) * timestep_conv
                Enddo ! i
              Enddo ! j
            Enddo ! k
          Endif                   ! on STASHflag

          If ( l_qcl_incr_conv ) Then
            Do k=1,n_conv_levels
              Do j=1,rows
                Do i=1,row_length
                  qcl_incr_diagnostic(i,j,k) =                          &
     &                     qcl_incr_diagnostic(i,j,k) +                 &
     &                         dqclbydt(i,j,k) * timestep_conv
                End Do ! i
              End Do ! j
            End Do
          Endif

          If ( l_qcf_incr_conv ) Then
            Do k=1,n_conv_levels
              Do j=1,rows
                Do i=1,row_length
                  qcf_incr_diagnostic(i,j,k) =                          &
     &                     qcf_incr_diagnostic(i,j,k) +                 &
     &                            dqcfbydt(i,j,k) * timestep_conv
                End Do ! i
              End Do ! j
            End Do
          Endif

          If ( l_cfl_incr_conv ) Then
            Do k=1,n_conv_levels
              Do j=1,rows
                Do i=1,row_length
                  cf_liquid_incr_diagnostic(i,j,k) =                    &
     &                          cf_liquid_incr_diagnostic(i,j,k) +      &
     &                          dcflbydt(i,j,k) * timestep_conv
                End Do ! i
              End Do ! j
            End Do
          Endif

          If ( l_cff_incr_conv ) Then
            Do k=1,n_conv_levels
              Do j=1,rows
                Do i=1,row_length
                  cf_frozen_incr_diagnostic(i,j,k) =                    &
     &                          cf_frozen_incr_diagnostic(i,j,k) +      &
     &                          dcffbydt(i,j,k) * timestep_conv
                End Do ! i
              End Do ! j
            End Do ! k
          Endif

          If ( l_bcf_incr_conv ) Then
            Do k=1,n_conv_levels
              Do j=1,rows
                Do i=1,row_length
                  bulk_cf_incr_diagnostic(i,j,k) =                      &
     &                       bulk_cf_incr_diagnostic(i,j,k)   +         &
     &                         (dbcfbydt(i,j,k) * timestep_conv)
                End Do ! i
              End Do ! j
            End Do ! k
          Endif

          If (l_mom) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  u_p(i,j,k) =u_p(i,j,k)+dubydt_p(i,j,k)*timestep_conv
                  v_p(i,j,k) =v_p(i,j,k)+dvbydt_p(i,j,k)*timestep_conv
                End Do
              End Do
            End Do
            If (l_up_incr_conv) then
              Do k = 1, n_conv_levels
                Do j = 1, rows
                  Do i = 1, row_length
                  dubydt_pout(i,j,k) =dubydt_pout(i,j,k)                &
     &                        +dubydt_p(i,j,k)*one_over_conv_calls
                  End Do
                End Do
              End Do
            Endif
            If (l_vp_incr_conv) then
              Do k = 1, n_conv_levels
                Do j = 1, rows
                  Do i = 1, row_length
                  dvbydt_pout(i,j,k) =dvbydt_pout(i,j,k)                &
     &                        +dvbydt_p(i,j,k)*one_over_conv_calls
                  End Do
                End Do
              End Do
            Endif
          End If

          If (flg_up_flx) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  up_flux(i,j,k)=up_flux(i,j,k)+it_up_flux(i,j,k)*      &
     &                                       one_over_conv_calls
                End Do
              End Do
            End Do
          End If
!
          If (flg_up_flx_half) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  up_flux_half(i,j,k)=up_flux_half(i,j,k)               &
     &                +it_up_flux_half(i,j,k)* one_over_conv_calls
                End Do
              End Do
            End Do
          End If
!
          If (flg_dwn_flx) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  dwn_flux(i,j,k)=dwn_flux(i,j,k)+it_dwn_flux(i,j,k)*   &
     &                                       one_over_conv_calls
                End Do
              End Do
            End Do
          End If
!
          If (flg_entr_up) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  entrain_up(i,j,k)=entrain_up(i,j,k)+                  &
     &                   it_entrain_up(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          End If
!
          If (flg_entr_dwn) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  entrain_dwn(i,j,k)=entrain_dwn(i,j,k)+                &
     &                   it_entrain_dwn(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          End If
!
          If (flg_detr_up) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  detrain_up(i,j,k)=detrain_up(i,j,k)+                  &
     &                   it_detrain_up(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          End If
!
          If (flg_detr_dwn) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  detrain_dwn(i,j,k)=detrain_dwn(i,j,k)+                &
     &                   it_detrain_dwn(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          End If
!
          If (flg_uw_dp) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  uw_dp(i,j,k)=uw_dp(i,j,k)+                            &
     &                   it_uw_dp(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          End If
!
          If (flg_vw_dp) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  vw_dp(i,j,k)=vw_dp(i,j,k)+                            &
     &                   it_vw_dp(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          End If
!
          If (flg_uw_shall) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  uw_shall(i,j,k)=uw_shall(i,j,k)+                      &
     &                   it_uw_shall(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          End If
!
          If (flg_vw_shall) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  vw_shall(i,j,k)=vw_shall(i,j,k)+                      &
     &                   it_vw_shall(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          End If

          If (flg_wqt_flux) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  wqt_flux(i,j,k)=wqt_flux(i,j,k)+                      &
     &                   it_wqt_flux(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
              Do j = 1, rows
                Do i = 1, row_length
                  k=ntml(i,j)+1
                  wqt_cb(i,j)=wqt_cb(i,j)+                              &
     &                   it_wqt_flux(i,j,k)*one_over_conv_calls
                End Do
              End Do
              Do j = 1, rows
                Do i = 1, row_length
                  k=ntpar(i,j)+1
                  wqt_inv(i,j)=wqt_inv(i,j)+                            &
     &                   it_wqt_flux(i,j,k)*one_over_conv_calls
                End Do
              End Do
          End If
          If (flg_wql_flux) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  wql_flux(i,j,k)=wql_flux(i,j,k)+                      &
     &                   it_wql_flux(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          End If
          If (flg_wthetal_flux) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  wthetal_flux(i,j,k)=wthetal_flux(i,j,k)+              &
     &                   it_wthetal_flux(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
              Do j = 1, rows
                Do i = 1, row_length
                  k=ntml(i,j)+1
                  wthetal_cb(i,j)=wthetal_cb(i,j)+                      &
     &                   it_wthetal_flux(i,j,k)*one_over_conv_calls
                End Do
              End Do
              Do j = 1, rows
                Do i = 1, row_length
                  k=ntpar(i,j)+1
                  wthetal_inv(i,j)=wthetal_inv(i,j)+                    &
     &                   it_wthetal_flux(i,j,k)*one_over_conv_calls
                End Do
              End Do
          End If
          If (flg_wthetav_flux) Then
            Do k = 1, n_conv_levels
              Do j = 1, rows
                Do i = 1, row_length
                  wthetav_flux(i,j,k)=wthetav_flux(i,j,k)+              &
     &                   it_wthetav_flux(i,j,k)*one_over_conv_calls
                End Do
              End Do
            End Do
          End If


!
! diagnose number of convecting points
          If (PrintStatus == PrStatus_Diag) Then
            If (n_conv_calls  >   1) Then
              n_conv_points = 0
              Do j = 1, rows
                Do i = 1, row_length
                  this_point = 0
                  Do k = 1, n_conv_levels
                    If (abs(dthbydt(i,j,k) * timestep)  >   0.0001) Then
                      this_point = 1
                    End If
                  End Do
                  n_conv_points = n_conv_points + this_point
                End Do
              End Do
              If (n_proc  >   1) Then
                Call gc_isum(1, n_proc, info, n_conv_points )
              End If
              If (me  ==  0) Then
                write(6,*) ' conv call ',call_number,' has ',           &
     &                     n_conv_points,' convecting points '
              End If

            End If
          End If


! ----------------------------------------------------------------------
! Section CNV.4 Move u, v, increments to C grid.
! ----------------------------------------------------------------------
          If (l_mom) Then
!------------------------------
! PLACE THE MAIN SCV CODE HERE
!------------------------------

!--------------------------------------------------------------
! FIRST OF ALL ZERO ALL THE INCREMENTS JUST IN CASE WE MUCK UP
!--------------------------------------------------------------
      IF( L_SCV)then
        DO L=1,model_levels
          DO J=1,rows
            DO I=1,row_length
              SCB_DUDT(I,J,L) = 0.0
              SCB_DVDT(I,J,L) = 0.0
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF (call_number == n_conv_calls) THEN
      IF( L_SCV .AND. L_CAPE )THEN
!
!------------------------------
! CALCULATE LONG,LAT
!------------------------------
       DO J=1,rows
        DO I=1,row_length
         latitude(i,j)=acos(FV_cos_theta_latitude(i,j))
         longitude(i,j)=acos(cos_theta_longitude(i,j))
         coriolis(i,j)=abs(0.5*(f3_at_u(i,j)+f3_at_u(i-1,j)))
        ENDDO
       ENDDO
!
! LET US DEFINE THE VALUES OF THE VARIABLES THAT CONTROL THE
! SIZE AND SHAPE OF THE VORTICITY SOURCES
!
      SCB_TOR=3.*3600.  ! THREE HOUR TIMESCALE
      TSC=1./SCB_TOR
      ZDEPTH_A=1500.    ! HALF DEPTHS
      ZDEPTH_C=2000.
      C_ETA=0.5         ! RADIUS OF MCV/RADIUS OF ANVIL
      RCETA=1./C_ETA
      RCETA2=1./(C_ETA**2)
      DADV=ZDEPTH_A/ZDEPTH_C  ! RATIO OF ANVIL AND MCV DEPTHS
!
!-------------------------------------------
! RUN THE STOCHASTIC CONVECTIVE BACKSCATTER
!-------------------------------------------
!
!-------------------------------------------------------------
! DETERMINE WHETHER TO UPDATE THE VORTICITY SOURCE INCREMENTS
!-------------------------------------------------------------

        minutos=timestep/60
        IF (MOD(val_HOUR,ISCBFREQ) == 0                                 &
     &      .AND. val_MINUTE == minutos) THEN
         L_SCV_UP=.true.
        ELSE
         L_SCV_UP=.false.
        ENDIF
!
        IF (L_SCV_UP) THEN
         IF (PrintStatus  >=  PrStatus_Normal) THEN
          WRITE(6,*) 'Inside L_SCV_UP ..',val_hour,val_minute
         ENDIF
!

!------------------------------
! UPDATE THE VORTICITY SOURCES
!--------------------------------------------------------------
! AND THE MCS PROPERTIES
          DO I=1,NMCSMAX*N_PROC
            DO N=1,NMCSPROP
              SCB_MCS_PROP(I,N)=0.0
            ENDDO
          ENDDO
!
!----------------------------------------------------------------
! IDENTIFY POTENTIAL VORTICITY SOURCES AND CALCULATE THEIR RADII
!----------------------------------------------------------------
          NMCS = 0
          DO J=1,rows
           DO I=1,row_length
            IF(conv_rain(I,J) >  0.0)THEN
              call random_number(azar)
              IF(cape_out(I,J) <  0.01) THEN
               cape_out(I,J)=0.01
              ENDIF
              SCB_FCORIOL = MAX(2.0E-5,coriolis(I,J))
              SCB_RAD = (0.1* 3600.*conv_rain(I,J)*azar                 &
     &                  *SQRT(cape_out(I,J)))/SCB_FCORIOL
!
! MULTIPLICATION BY 3600 ABOVE CONVERTS RAIN INTO MM/H
! ONLY CONSIDER POINTS WITH CCT ABOVE 350mb
!
              L=cct(I,J)  ! LEVEL OF CONVECTIVE CLOUD TOP
              PRESSURE = P_LAYER_CENTRES(I,J,L)
!
              IF(SCB_RAD >= SCB_MINRAD .AND. NMCS <  NMCSMAX            &
     &           .AND. PRESSURE <  35000.0 .AND. L /= 0 )THEN
                NMCS = NMCS + 1
                SCB_MCS_PROP(NMCS+ME*NMCSMAX,1) = LATITUDE(I,J)
                SCB_MCS_PROP(NMCS+ME*NMCSMAX,2) = LONGITUDE(I,J)
                SCB_MCS_PROP(NMCS+ME*NMCSMAX,3) = SCB_RAD
!
! FOR THE CENTRE OF THE MID-LEVEL CYCLONE TAKE 50MB
! ABOVE THE FREEZING LEVEL
!
                L=1
                PRESSURE = P_LAYER_CENTRES(I,J,L)
                TEMPERATURE=theta_n(I,J,L)*(PRESSURE/p_zero)            &
     &                       **(Kappa)
                DO WHILE(TEMPERATURE  >   273.16 )
                 L=L+1
                 PRESSURE = P_LAYER_CENTRES(I,J,L)
                 TEMPERATURE=theta_n(I,J,L)*(PRESSURE/p_zero)           &
     &                       **(Kappa)
                ENDDO
                SCB_MCS_PROP(NMCS+ME*NMCSMAX,4)=(PRESSURE-5000.)
! CONVERT ZDEPTH_C TO PRESSURE
                PRESSURE = SCB_MCS_PROP(NMCS+ME*NMCSMAX,4)
                TEMPERATURE=theta_n(I,J,L)*(PRESSURE/p_zero)            &
     &                      **(Kappa)
                SCB_MCS_PROP(NMCS+ME*NMCSMAX,6) =                       &
     &               (PRESSURE*G/(R*TEMPERATURE))*ZDEPTH_C
!
! FOR THE CENTRE OF THE UPPER LEVEL ANTICYCLONE TAKE 50 MB BELOW THE
! DIAGNOSED CONVECTIVE CLOUD TOP HEIGHT
!
                L=cct(I,J)  ! LEVEL OF CONVECTIVE CLOUD TOP
                PRESSURE = P_LAYER_CENTRES(I,J,L)
                SCB_MCS_PROP(NMCS+ME*NMCSMAX,5)=(PRESSURE+5000.)
! CONVERT ZDEPTH_A TO PRESSURE
                PRESSURE = SCB_MCS_PROP(NMCS+ME*NMCSMAX,5)
                TEMPERATURE=theta_n(I,J,L)*(PRESSURE/p_zero)            &
     &                      **(Kappa)
                SCB_MCS_PROP(NMCS+ME*NMCSMAX,7) =                       &
     &               (PRESSURE*G/(R*TEMPERATURE))*ZDEPTH_A
!
              ENDIF ! (if PRESSURE >  35000.0 .AND. L /= 0
            ENDIF  ! (if conv. rain gt 0)
           ENDDO ! end i
          ENDDO ! end j
!
!-----------------------------------------------------------------
! NOW FOR EACH INDIVIDUAL PE RANK THE VORTICITY SOURCES IN ORDER
! OF SIZE USING A SHUFFLE ALGORITHM.
!-----------------------------------------------------------------
!
         IF (PrintStatus  >   PrStatus_Normal) THEN
          WRITE(6,*) 'NMCS: ',nmcs
         ENDIF
!
          IF(NMCS >  0)THEN
            ISCBCHANGE=1
            DO WHILE(ISCBCHANGE == 1)
              ISCBCHANGE=0
              DO I=ME*NMCSMAX+2,ME*NMCSMAX+NMCS
                IF(SCB_MCS_PROP(I-1,3) <  SCB_MCS_PROP(I,3))THEN
                  ISCBCHANGE=1
! SWAP OVER THE SCB_MCS_PROP VALUES
                  DO L=1,NMCSPROP
                    SCB_TEMP=SCB_MCS_PROP(I-1,L)
                    SCB_MCS_PROP(I-1,L)=SCB_MCS_PROP(I,L)
                    SCB_MCS_PROP(I,L)=SCB_TEMP
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
!
!----------------------------------------------------
! NOW EXCLUDE ALL VORTICITY SOURCES WITHIN 2*RADIUS
! OF A LARGER VORTICITY SOURCE.
!-----------------------------------------------------------------------
! THE DISTANCE BETWEEN ANY TWO POINTS ON A SPHERE IS CALCULATED FROM
!   COS(D) = SIN(A1)SIN(A2) + COS(A1)COS(A2)COS(B1-B2)
! WHERE A IS THE LATITUDE, B IS THE LONGITUDE AND D IS THE ANGLE BETWEEN
! THE TWO POINTS.
! THIS WILL BE USED TO WORK OUT THE DISTANCE BETWEEN TWO POINTS.
!-----------------------------------------------------------------------
           IF(NMCS >  1)THEN
             DO M=1+ME*NMCSMAX,NMCS+ME*NMCSMAX
                DO I=M+1,NMCS+ME*NMCSMAX
                  LAT1=SCB_MCS_PROP(M,1)
                  LAT2=SCB_MCS_PROP(I,1)
                  LONG1=SCB_MCS_PROP(M,2)
                  LONG2=SCB_MCS_PROP(I,2)
! DEPENDS ON: scb_dist
                  DIST=SCB_DIST(LAT1,LAT2,LONG1,LONG2)
!------------------------------------------------
! COMPARE THE DISTANCE BETWEEN VORTICITY SOURCES
! WITH THE RADIUS OF THE LARGEST SOURCE.
!------------------------------------------------
                  IF( DIST <= 2.*SCB_MCS_PROP(M,3) )THEN
                    DO L=1,NMCSPROP
                      SCB_MCS_PROP(I,L)=0.0
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
!
         ENDIF   ! NMCS >  0
!
!----------------------------------------------------------------------
! NOW COLLECT THE SCB_MCS_PROP TOGETHER ON ALL PES USING GCOM ROUTINES
! BEWARE! GC_RSUM is used here in a non-standard way: SCB_MCS_PROP is
! defined as a (400,7) matrix but each processor only writes on 50
! lines. Thus, PE0 1-50; PE1 51-100, etc. Thus GC_RSUM can be used
! to group and distribute all data from all processors.
!----------------------------------------------------------------------
          CALL GC_RSUM(NMCSMAX*N_PROC*NMCSPROP,                         &
     &                 N_PROC,INFO,SCB_MCS_PROP)
!-----------------------------------------------------------------
! NOW SCB_MCS_PROP SHOULD BE ONE BIG ARRAY WITH ALL THE VORTICITY
! SOURCE INFORMATION IN. FIRST REMOVE ALL THE ZEROS AND THEN SORT
! INTO RANK ORDER WITH LARGEST MCS IN (1,*) ETC.
!-----------------------------------------------------------------
          NMCS=0
          DO I=1,NMCSMAX*N_PROC
            IF( SCB_MCS_PROP(I,3) /= 0.0 )THEN
              NMCS=NMCS+1
              IF(NMCS <  I)THEN
                DO L=1,NMCSPROP
                  SCB_MCS_PROP(NMCS,L)=SCB_MCS_PROP(I,L)
                  SCB_MCS_PROP(I,L)=0.0
                ENDDO
              ENDIF
            ENDIF
          ENDDO
! PRINT OUT ALL THE DETAILS OF THE SCB DATA
         IF (PrintStatus  >   PrStatus_Normal) THEN
          WRITE(6,*) 'GLOBAL SCB DATA FOLLOWS: NMCS = ',NMCS
         ENDIF
!----------------------------------------------------
! NOW RANK THEM INTO ORDER USING A SHUFFLE ALGORITHM
!----------------------------------------------------
          ISCBCHANGE=1
          DO WHILE(ISCBCHANGE == 1)
            ISCBCHANGE=0
            DO I=2,NMCS
              IF(SCB_MCS_PROP(I-1,3) <  SCB_MCS_PROP(I,3))THEN
                ISCBCHANGE=1
! SWAP OVER THE SCB_MCS_PROP VALUES
                DO L=1,NMCSPROP
                  SCB_TEMP=SCB_MCS_PROP(I-1,L)
                  SCB_MCS_PROP(I-1,L)=SCB_MCS_PROP(I,L)
                  SCB_MCS_PROP(I,L)=SCB_TEMP
                ENDDO
              ENDIF
            ENDDO
          ENDDO
!
!
!----------------------------------------------------
! NOW EXCLUDE ALL VORTICITY SOURCES WITHIN 2*RADIUS
! OF A LARGER VORTICITY SOURCE.
!----------------------------------------------------------------------
! THE DISTANCE BETWEEN ANY TWO POINTS ON A SPHERE IS CALCULATED FROM
!   COS(D) = SIN(A1)SIN(A2) + COS(A1)COS(A2)COS(B1-B2)
! WHERE A IS THE LATITUDE, B IS THE LONGITUDE AND D IS THE ANGLE BETWEEN
! THE TWO POINTS.
! THIS WILL BE USED TO WORK OUT THE DISTANCE BETWEEN TWO POINTS.
!----------------------------------------------------------------------
          DO M=1,NMCS
            IF(SCB_MCS_PROP(M,3) >  0.0)THEN
! MAKES SURE DON'T DO MCSs ALREADY EXCLUDED
              DO I=M+1,NMCS
                IF(SCB_MCS_PROP(I,3) >  0.0)THEN
                  LAT1=SCB_MCS_PROP(M,1)
                  LAT2=SCB_MCS_PROP(I,1)
                  LONG1=SCB_MCS_PROP(M,2)
                  LONG2=SCB_MCS_PROP(I,2)
! DEPENDS ON: scb_dist
                  DIST=SCB_DIST(LAT1,LAT2,LONG1,LONG2)
!------------------------------------------------
! COMPARE THE DISTANCE BETWEEN VORTICITY SOURCES
! WITH THE RADIUS OF THE LARGEST SOURCE.
!------------------------------------------------
                  IF( DIST <= 2.*SCB_MCS_PROP(M,3) )THEN
                    DO L=1,NMCSPROP
                      SCB_MCS_PROP(I,L)=0.0
                    ENDDO
                  ENDIF
                ENDIF    ! TEST MCS NOT EXCLUDED ALREADY
              ENDDO
            ENDIF        ! TESTING MCS NOT EXCLUDED ALREADY
          ENDDO
!
!------------------------------------------------
! FINALLY REMOVE ALL THE ZEROS AGAIN AND RECOUNT
!------------------------------------------------
          NMCSTOT=NMCS
          NMCS=0
          DO I=1,NMCSTOT
            IF( SCB_MCS_PROP(I,3) /= 0.0 )THEN
              NMCS=NMCS+1
              IF(NMCS <  I)THEN
                DO L=1,NMCSPROP
                  SCB_MCS_PROP(NMCS,L)=SCB_MCS_PROP(I,L)
                  SCB_MCS_PROP(I,L)=0.0
                ENDDO
              ENDIF
            ENDIF
          ENDDO
! PRINT OUT ALL THE DETAILS OF THE SCB DATA
         IF (PrintStatus  >   PrStatus_Normal) THEN
            WRITE(6,*)'GLOBAL SCB DATA FOLLOWS:'
            WRITE(6,*)'NMCS after exclusion = ',NMCS
            DO I=1,NMCS
              WRITE(6,*)I,(SCB_MCS_PROP(I,L),L=1,NMCSPROP)
            ENDDO
         ENDIF
!
!----------------------------------------------------------------
! SO AT THIS POINT WE HAVE ALL THE RADII, AND MCV AND
! ANTICYCLONE LEVELS FOR ALL THE VORTICITY SOURCES THAT
! WE ARE GOING TO CONSIDER IN THIS PERIOD. WITH THIS INFORMATION
! WE CAN NOW DETERMINE THE FIELD OF MOMENTUM INCREMENTS THAT
! CORRESPOND TO THE STOCHASTIC CONVECTIVE BACKSCATTER.
!----------------------------------------------------------------
! LOOP OVER ALL POINTS, THEN LOOP OVER ALL MCS VORTICITY SOURCES
!
          DO J=1,rows
           DO I=1,row_length
            DO N=1,NMCS
! CHECK IF GRID POINT IS WITHIN RANGE OF VORTICITY SOURCE
              LAT1=SCB_MCS_PROP(N,1)
              LAT2=LATITUDE(I,J)
              LONG1=SCB_MCS_PROP(N,2)
              LONG2=LONGITUDE(I,J)
! DEPENDS ON: scb_dist
              DIST=SCB_DIST(LAT1,LAT2,LONG1,LONG2)
              IF(DIST <= (2.*SCB_MCS_PROP(N,3)))THEN
! CALCULATE THE MID LEVELS AND DEPTH OF ANTICYCLONE
                PMID_A = SCB_MCS_PROP(N,5)
                PDEPTH_A = SCB_MCS_PROP(N,7)
! AND NOW FOR THE CYCLONE
                PMID_C = SCB_MCS_PROP(N,4)
                PDEPTH_C = SCB_MCS_PROP(N,6)
                SCB_FCORIOL = coriolis(I,J)  ! LOCAL CORIOLIS
!
! IFS TO AVOID DIVISION BY 0
               IF(dist <= 0.01) THEN
                dist=0.01
               ENDIF
               IF(PDEPTH_A <= 0.01) THEN
                PDEPTH_A=0.01
               ENDIF
               IF(PDEPTH_C <= 0.01) THEN
                PDEPTH_C=0.01
               ENDIF
!
! NOW TO WORK OUT THE GREAT CIRCLE DISTANCES FROM CENTRE TO POINT
! IN X AND Y DIRECTION
                LAT1=SCB_MCS_PROP(N,1)
                LAT2=LAT1
                LONG1=SCB_MCS_PROP(N,2)
                LONG2=LONGITUDE(I,J)
! DEPENDS ON: scb_dist
                XDIST = SCB_DIST(LAT1,LAT2,LONG1,LONG2)
                IF( LONG1 >  LONG2 )XDIST=-XDIST
!
                LAT1=SCB_MCS_PROP(N,1)
                LAT2=LATITUDE(I,J)
                LONG1=SCB_MCS_PROP(N,2)
                LONG2=LONG1
! DEPENDS ON: scb_dist
                YDIST = SCB_DIST(LAT1,LAT2,LONG1,LONG2)
                IF( LAT1 >  LAT2 )YDIST=-YDIST
! GRID POINT WITHIN VORTICITY SOURCE FOR MCS N
                DO L=1,wet_model_levels
                  PRESSURE = P_LAYER_CENTRES(I,J,L)
! FIRST DO UPPER-LEVEL ANTICYCLONE
                  IF( PRESSURE >  (PMID_A-PDEPTH_A) .AND.               &
     &                PRESSURE <  (PMID_A+PDEPTH_A) )THEN
                    PFAC = COS((0.5*PI*(PRESSURE - PMID_A))/PDEPTH_A)
                    IF( DIST <  (PFAC*SCB_MCS_PROP(N,3)))THEN
                      SCB_DUDT(I,J,L) = SCB_DUDT(I,J,L) +               &
     &                  (0.5*SCB_FCORIOL*YDIST*TSC)
                      SCB_DVDT(I,J,L) = SCB_DVDT(I,J,L) -               &
     &                  (0.5*SCB_FCORIOL*XDIST*TSC)
                    ELSEIF( DIST <= (2.*PFAC*SCB_MCS_PROP(N,3)))THEN
                      SCB_DUDT(I,J,L) = SCB_DUDT(I,J,L) +               &
     &           (0.5*SCB_FCORIOL*SCB_MCS_PROP(N,3)*PFAC*(YDIST/DIST)   &
     &                *((PFAC*SCB_MCS_PROP(N,3)/DIST)**2.)*TSC)
                      SCB_DVDT(I,J,L) = SCB_DVDT(I,J,L) -               &
     &           (0.5*SCB_FCORIOL*SCB_MCS_PROP(N,3)*PFAC*(XDIST/DIST)   &
     &                *((PFAC*SCB_MCS_PROP(N,3)/DIST)**2.)*TSC)
                    ENDIF
                  ENDIF
! AND NOW THE MID-LEVEL CYCLONE
                  IF( (PRESSURE >  (PMID_C-PDEPTH_C)) .AND.             &
     &                (PRESSURE <  (PMID_C+PDEPTH_C)) )THEN
                    PFAC = COS((0.5*PI*(PRESSURE - PMID_C))/PDEPTH_C)
                    IF( DIST <  (C_ETA*SCB_MCS_PROP(N,3)*PFAC))THEN
                      SCB_DUDT(I,J,L) = SCB_DUDT(I,J,L) -               &
     &                 (0.5*SCB_FCORIOL*YDIST*DADV*RCETA2*TSC)
                      SCB_DVDT(I,J,L) = SCB_DVDT(I,J,L) +               &
     &                 (0.5*SCB_FCORIOL*XDIST*DADV*RCETA2*TSC)
                    ELSEIF                                              &
     &               ( DIST <= (2.*C_ETA*SCB_MCS_PROP(N,3)*PFAC))THEN
                      SCB_DUDT(I,J,L) = SCB_DUDT(I,J,L) -               &
     &                 (0.5*SCB_FCORIOL*SCB_MCS_PROP(N,3)*DADV*RCETA    &
     &                *PFAC*(YDIST/DIST)*                               &
     &                 ((PFAC*SCB_MCS_PROP(N,3)*C_ETA/DIST)**2.)*TSC)
                      SCB_DVDT(I,J,L) = SCB_DVDT(I,J,L) +               &
     &                 (0.5*SCB_FCORIOL*SCB_MCS_PROP(N,3)*DADV*RCETA    &
     &                *PFAC*(XDIST/DIST)*                               &
     &                 ((PFAC*SCB_MCS_PROP(N,3)*C_ETA/DIST)**2.)*TSC)
                    ENDIF
                  ENDIF
                ENDDO  ! L=1,wet_model_levels
              ENDIF    ! IF WITHIN VORTICITY RADIUS
            ENDDO      ! N - NUMBER OF VORTICITY SOURCES
           ENDDO       ! I - NUMBER OF X GRID POINTS
          ENDDO        ! J - NUMBER OF Y GRID POINTS
!
!-----------------------------------------------------------------
! NOW ADD THE SCB U,V INCREMENTS TO THE CONVECTIVE MOMENTUM TRANSPORTS
!-----------------------------------------------------------------
!
        ENDIF  ! L_SCV_UP (Substitued by val_hour IF ....)
!
      ENDIF    ! L_SCV .AND. L_CAPE
      ENDIF    ! call_number=n_conv_calls
!

! U increments
! first need to copy increments into arrays with one point haloes
! and swop boundaries

           IF(L_SCV) THEN
            DO k = 1, n_conv_levels + 1
              DO j = 1, rows
                DO i = 1, row_length
                 work_p_halo(i,j,k) = dubydt_p(i,j,k)+scb_dudt(i,j,k)
                ENDDO
              ENDDO
            ENDDO
           ELSE
            DO k = 1, n_conv_levels + 1
              DO j = 1, rows
                DO i = 1, row_length
                  work_p_halo(i,j,k) = dubydt_p(i,j,k)
                ENDDO
              ENDDO
            ENDDO
           ENDIF

#if !defined(SCMA)
! DEPENDS ON: swap_bounds
            Call Swap_Bounds(                                           &
     &           work_p_halo, row_length, rows,                         &
     &           n_conv_levels+1, off_x, off_y, fld_type_p,             &
     &           .false.  )
! DEPENDS ON: fill_external_halos
            Call FILL_EXTERNAL_HALOS(WORK_P_HALO,ROW_LENGTH,ROWS,       &
     &                               N_CONV_LEVELS+1,off_x,off_y)
#endif

! interpolate to u grid
! DEPENDS ON: p_to_u
            Call p_to_u (work_p_halo, row_length, rows, n_conv_levels+1,&
     &                 off_x, off_y, dubydt_u)

! add on to increment field
            If ( l_u_incr_conv ) Then
! update R_u and diagnostics
              Do k=1,n_conv_levels+1
                Do j=1,rows
                  Do i=1,row_length
                    R_u(i,j,k) = R_u(i,j,k)                             &
     &                         + dubydt_u(i,j,k)*timestep_conv
                    u_incr_diagnostic(i,j,k) = u_incr_diagnostic(i,j,k) &
     &                                      + dubydt_u(i,j,k)           &
     &                                      * timestep_conv
                  Enddo ! i
                Enddo ! j
              Enddo ! k
            Else
              Do k = 1, n_conv_levels+1
                Do j = 1, rows
                  Do i = 1, row_length
                    R_u(i,j,k)=R_u(i,j,k)+dubydt_u(i,j,k)*timestep_conv
                  End Do
                End Do
              End Do
            Endif                   ! on STASHflag

! V increments
! first need to copy increments into arrays with one point haloes
! and swop boundaries

           IF(L_SCV) THEN
            DO k = 1, n_conv_levels+1
              DO j = 1, rows
                DO i = 1, row_length
                  work_p_halo(i,j,k) = dvbydt_p(i,j,k)+scb_dvdt(i,j,k)
                ENDDO
              ENDDO
            ENDDO
           ELSE
            DO k = 1, n_conv_levels+1
              DO j = 1, rows
                DO i = 1, row_length
                  work_p_halo(i,j,k) = dvbydt_p(i,j,k)
                ENDDO
              ENDDO
            ENDDO
           ENDIF

#if !defined(SCMA)
! DEPENDS ON: swap_bounds
            Call Swap_Bounds(                                           &
     &           work_p_halo, row_length, rows,                         &
     &           n_conv_levels+1, off_x, off_y, fld_type_p,             &
     &           .false.  )
! DEPENDS ON: fill_external_halos
            Call FILL_EXTERNAL_HALOS(WORK_P_HALO,ROW_LENGTH,ROWS,       &
     &                               N_CONV_LEVELS+1,off_x,off_y)
#endif

! interpolate to v grid
! DEPENDS ON: p_to_v
            Call p_to_v (work_p_halo, row_length, rows, n_rows,         &
     &                   n_conv_levels+1, off_x, off_y, dvbydt_v)

! add on to increment field
            If ( l_v_incr_conv ) Then
              Do k=1, n_conv_levels+1
                Do j=1,n_rows
                  Do i=1,row_length
                     R_v(i,j,k) = R_v(i,j,k)                            &
     &                          + dvbydt_v(i,j,k)*timestep_conv
                     v_incr_diagnostic(i,j,k) = v_incr_diagnostic(i,j,k)&
     &                                       + dvbydt_v(i,j,k)          &
     &                                       * timestep_conv
                 Enddo ! i
               Enddo ! j
             Enddo ! k
            Else
              Do k = 1, n_conv_levels+1
                Do j = 1, n_rows
                  Do i = 1, row_length
                   R_v(i,j,k)=R_v(i,j,k)+dvbydt_v(i,j,k)*timestep_conv
                  End Do
                End Do
              End Do
            Endif                   ! on STASHflag

          End If
        End Do ! loop over number of convection calls

!-----------------------------------------------------------------------
! HadGEM1 cloud changes - the location of this is not ideal, it should
! have been done in glue but due to what was used in the original
! modset it has to remain here to give the same results.
! The problem comes from the fact that HadGEM1 uses the cloud top
! and bottoms from the first call to convection, the freezing level
! from the last call to convection and rain and snow from all substeps.
!-----------------------------------------------------------------------

       If (l_convcld_hadgem1) then
         Do k = 1, n_cca_levels
           Do j = 1, rows
             Do i = 1, row_length
               cca(i,j,k) = 0.0
             End Do
           End Do
         End Do
          Do j = 1, rows
            Do i = 1, row_length
                cca_2d(i,j)=0.0
            End Do
          End Do
!
! modify shallow cloud amount and depth
!
         cct_change_map(:,:) = .FALSE.

         Do j = 1, rows
           Do i = 1, row_length
             If (l_shallow(i,j).and.ccb(i,j) /= 0) then
              mb_temp=0.03*wstar(i,j)
              wsc_temp=(delthvu(i,j)*mb_temp*g/                         &
     & (theta_n(i,j,ntml(i,j))*(1.0+0.61*q_n(i,j,ntml(i,j)))))**0.3333
               If ( lcbase(i,j)  ==  ccb(i,j)) then
!
! Record where we make these changes as will need to correct cct later
!
                 cct_change_map(i,j) = .TRUE.
!
! reset cloud base and top, calculate new cloud fraction
!
                 ccb(i,j)=ntml(i,j)+1
!
! cloud top is defined as level above cloud in 3d cal (cloud top level
! ntpar+1)
!
                 cct(i,j)=ntpar(i,j)+2

                 cca_2d(i,j)= 2.0*mb_temp/wsc_temp
                 lcbase(i,j)=ccb(i,j)
                 lctop(i,j)=cct(i,j)
                 lcca(i,j)=2.0*mb_temp/wsc_temp
               Else
!
! shallow cumulus not only cloud type, only change lowest cloud
!
                 lcbase(i,j)=ntml(i,j)+1
                 lctop(i,j)=ntpar(i,j)+1
                 lcca(i,j)=2.0*mb_temp/wsc_temp
               Endif
              Endif
            End Do
          End Do
!
! Deep cloud changes
!
          Do j = 1, rows
            Do i = 1, row_length

             If ( .not. l_shallow(i,j) .and. ccb(i,j)  /=  0.0) then
               If((conv_rain(i,j)+conv_snow(i,j))  >   0.0) Then

! Reduce cloud liquid water for HadGEM1

                 cclwp(i,j) = reduction_factor * cclwp(i,j)

                 If (land_sea_mask(i,j)) then
                   cca_2d(i,j)=a_land+b_land*alog(24.*3600.*            &
     &                        (conv_rain(i,j)+conv_snow(i,j)))
                 Else
                   cca_2d(i,j)=a_sea+b_sea*alog(24.*3600.*              &
     &                      (conv_rain(i,j)+conv_snow(i,j)))
                 Endif
                 If(cumulus(i,j)) then
!
! set cloud base and top to ntml and ntpar for diagnosed cumulus
!
                   If (cca_2d(i,j)  <=  0.0) then
                    ccb(i,j)=0
                    cct(i,j)=0
                   Else
                    ccb(i,j)=ntml(i,j)+1
                    cct(i,j)=ntpar(i,j)+1
                    lcbase(i,j)=ccb(i,j)
                    lctop(i,j)=cct(i,j)
                   Endif
                 Endif
               Endif
             Endif
           End Do
         End Do
         If (lcv_3d_cca) then
! DEPENDS ON: calc_3d_cca
           Call calc_3d_cca(rowsrowlength,rowsrowlength,n_conv_levels,  &
                         bl_levels, ccb,cct                             &
                        ,p_layer_boundaries,freeze_lev                  &
                        ,cca_2d,cca,z_theta,z_rho                       &
                        ,l_q_interact, .true. ,l_pc2_diag_sh_pts)

         Else
           Do j = 1, rows
             Do i = 1, row_length
              If (cca_2d(i,j)  /=  0.0) Then
                cca(i,j,1)=cca_2d(i,j)
              End If
             End Do
           End Do
         End If

!
! Correct cct
!
         Do j=1, rows
           Do i=1, row_length
             If (cct_change_map(i,j)) Then
               cct(i,j)=ntpar(i,j)+1
             End If
           End Do
         End Do

        Else      ! original code

          ! Check that CCA doesn't exceed 1.0
          Do k=1, n_cca_levels
            Do j=1, rows
              Do i=1, row_length
                 cca(i,j,k) = MIN(cca(i,j,k), 1.0)
              End Do
            End Do
          End Do

        End if      ! l_convcld_hadgem1

! ----------------------------------------------------------------------
! Section CNV 4.2 Copy increment fields for PC2 inhomog diagnostics.
! ----------------------------------------------------------------------
!
        If(l_qcl_incr_cinh) Then
!
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                qcl_incr_inhom_diag(i,j,k) =                            &
     &                      qcl_inc(i,j,k) - qcl_incr_inhom_diag(i,j,k)
              End Do
            End Do
          End Do
!
        Endif                   ! on STASHflag
!
        If(l_qcf_incr_cinh) Then
!
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                qcf_incr_inhom_diag(i,j,k) =                            &
     &                      qcf_inc(i,j,k) - qcf_incr_inhom_diag(i,j,k)
              End Do
            End Do
          End Do
!
        Endif                   ! on STASHflag
!
        If(l_bcf_incr_cinh) Then
!
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                bulk_cf_incr_inhom_diag(i,j,k) =                        &
     &              bulk_cf_inc(i,j,k) - bulk_cf_incr_inhom_diag(i,j,k)
              End Do
            End Do
          End Do
!
        Endif                   ! on STASHflag
!
        If(l_cfl_incr_cinh .or. L_calc_dxek) Then
!
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                cf_liquid_incr_inhom_diag(i,j,k) =                      &
     &          cf_liquid_inc(i,j,k) - cf_liquid_incr_inhom_diag(i,j,k)
              End Do
            End Do
          End Do
!
        Endif                   ! on STASHflag
!
        If(l_cff_incr_cinh) Then
!
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                cf_frozen_incr_inhom_diag(i,j,k) =                      &
     &          cf_frozen_inc(i,j,k) - cf_frozen_incr_inhom_diag(i,j,k)
              End Do
            End Do
          End Do
!
        Endif                   ! on STASHflag
!
!-----------------------------------------------------------------------
! Section CNV.5 Treat tracers
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
! Section CNV 5.1 Copy tracers back into variables from tot_tracers
!----------------------------------------------------------------------
       If ( L_Tracer ) Then
         ntra_tmp = 0


         IF (L_DUST) THEN
           NTRA_TMP = NTRA_TMP + 1
           DUST_DIV1(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS) =            &
     &     TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP)

           NTRA_TMP = NTRA_TMP + 1
           DUST_DIV2(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS) =            &
     &     TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP)

           NTRA_TMP = NTRA_TMP + 1
           DUST_DIV3(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS) =            &
     &     TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP)

           NTRA_TMP = NTRA_TMP + 1
           DUST_DIV4(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS) =            &
     &     TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP)

           NTRA_TMP = NTRA_TMP + 1
           DUST_DIV5(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS) =            &
     &     TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP)

           NTRA_TMP = NTRA_TMP + 1
           DUST_DIV6(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS) =            &
     &     TOT_TRACER(1:ROW_LENGTH, 1:ROWS, 1:MODEL_LEVELS, NTRA_TMP)
         ENDIF

         If (l_sulpc_SO2) Then
           ntra_tmp = ntra_tmp + 1
           so2       (1:row_length, 1:rows, 1:model_levels) =           &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)

           ntra_tmp = ntra_tmp + 1
           so4_aitken(1:row_length, 1:rows, 1:model_levels) =           &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)

           ntra_tmp = ntra_tmp + 1
           so4_accu  (1:row_length, 1:rows, 1:model_levels) =           &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)

           ntra_tmp = ntra_tmp + 1
           so4_diss  (1:row_length, 1:rows, 1:model_levels) =           &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)

           If (l_sulpc_dms) Then
             ntra_tmp = ntra_tmp + 1
             dms       (1:row_length, 1:rows, 1:model_levels) =         &
     &       tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)
           End If

           If (l_sulpc_nh3) Then
             ntra_tmp = ntra_tmp + 1
             nh3       (1:row_length, 1:rows, 1:model_levels) =         &
     &       tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)
           End If
         End If

         If (l_soot) Then
           ntra_tmp = ntra_tmp + 1
           soot_new  (1:row_length, 1:rows, 1:model_levels) =           &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)

           ntra_tmp = ntra_tmp + 1
           soot_agd  (1:row_length, 1:rows, 1:model_levels) =           &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)

           ntra_tmp = ntra_tmp + 1
           soot_cld  (1:row_length, 1:rows, 1:model_levels) =           &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)
         End If

         If (l_biomass) Then
           ntra_tmp = ntra_tmp + 1
           bmass_new  (1:row_length, 1:rows, 1:model_levels) =          &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)

           ntra_tmp = ntra_tmp + 1
           bmass_agd  (1:row_length, 1:rows, 1:model_levels) =          &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)

           ntra_tmp = ntra_tmp + 1
           bmass_cld  (1:row_length, 1:rows, 1:model_levels) =          &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)
         End If

         If (l_CO2_interactive) Then
           ntra_tmp = ntra_tmp + 1
           co2       (1:row_length, 1:rows, 1:model_levels) =           &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)
         End If
         
         If (l_ocff) Then
           ntra_tmp = ntra_tmp + 1
           ocff_new  (1:row_length, 1:rows, 1:model_levels) =           &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)
 
           ntra_tmp = ntra_tmp + 1
           ocff_agd  (1:row_length, 1:rows, 1:model_levels) =           &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)

           ntra_tmp = ntra_tmp + 1
           ocff_cld  (1:row_length, 1:rows, 1:model_levels) =           &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)
         End If
         
         If (L_USE_CARIOLLE) Then
           ntra_tmp = ntra_tmp + 1
           OZONE_TRACER (1:row_length, 1:rows, 1:model_levels) =        &
     &     tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)
         End If

         If (tr_vars > 0 ) Then
           free_tracers(1:row_length, 1:rows, 1:model_levels,1:tr_vars)=&
     &     tot_tracer  (1:row_length, 1:rows, 1:model_levels,           &
     &                  ntra_tmp+1:ntra_tmp + tr_vars )
           ntra_tmp = ntra_tmp + tr_vars
         End If

         If (tr_ukca > 0 ) Then
           ukca_tracers(1:row_length, 1:rows, 1:model_levels,1:tr_ukca)=&
     &     tot_tracer  (1:row_length, 1:rows, 1:model_levels,           &
     &                  ntra_tmp+1:ntra_tmp + tr_ukca )
           ntra_tmp = ntra_tmp + tr_ukca
         End If
       End If

! Free up the memory again
       Deallocate( tot_tracer )

!----------------------------------------------------------------------
! Section CNV 5.2 Scavenging for tracers
!----------------------------------------------------------------------

! work with tracers only in final cycle
      If ( CycleNo == NumCycles ) Then
! Scavenging for aerosol
        If (L_murk_source) Then
! DEPENDS ON: con_scav
          Call con_scav(                                                &
     &     timestep, row_length, rows, model_levels                     &
     &,    ccb, cct, conv_rain, conv_snow                               &
     &,    aerosol(1:row_length, 1:rows, 1:model_levels)                &
     &)
        End If

! Scavenging of other tracers should go here....
! Scavenge Mineral Dust tracers
!
        IF (L_DUST) THEN
!
! DEPENDS ON: scnscv2
          CALL SCNSCV2(                                                 &
     &           ROW_LENGTH, ROWS,                                      &
     &           OFF_X, OFF_Y, HALO_I, HALO_J,                          &
     &           MODEL_LEVELS, WET_MODEL_LEVELS,                        &
     &           R_RHO_LEVELS, R_THETA_LEVELS,                          &
     &           TIMESTEP, SUBSTEP_NUMBER,                              &
     &           RHO,                                                   &
     &           Q_N, QCL_N, QCF_N, DUST_DIV1,                          &
     &           CCB, CCT,                                              &
     &           CONV_RAIN, CONV_SNOW,                                  &
     &           .TRUE.,                                                &
     &           KRAIN_DUST(1), KSNOW_DUST(1),                          &
     &           CONSCAV_DUST(1,1,1)                                    &
     &           )
!
! DEPENDS ON: scnscv2
          CALL SCNSCV2(                                                 &
     &           ROW_LENGTH, ROWS,                                      &
     &           OFF_X, OFF_Y, HALO_I, HALO_J,                          &
     &           MODEL_LEVELS, WET_MODEL_LEVELS,                        &
     &           R_RHO_LEVELS, R_THETA_LEVELS,                          &
     &           TIMESTEP,SUBSTEP_NUMBER,                               &
     &           RHO,                                                   &
     &           Q_N, QCL_N, QCF_N, DUST_DIV2,                          &
     &           CCB, CCT,                                              &
     &           CONV_RAIN, CONV_SNOW,                                  &
     &           .TRUE.,                                                &
     &           KRAIN_DUST(2), KSNOW_DUST(2),                          &
     &           CONSCAV_DUST(1,1,2)                                    &
     &           )
!
! DEPENDS ON: scnscv2
          CALL SCNSCV2(                                                 &
     &           ROW_LENGTH, ROWS,                                      &
     &           OFF_X, OFF_Y, HALO_I, HALO_J,                          &
     &           MODEL_LEVELS, WET_MODEL_LEVELS,                        &
     &           R_RHO_LEVELS, R_THETA_LEVELS,                          &
     &           TIMESTEP,SUBSTEP_NUMBER,                               &
     &           RHO,                                                   &
     &           Q_N, QCL_N, QCF_N, DUST_DIV3,                          &
     &           CCB, CCT,                                              &
     &           CONV_RAIN, CONV_SNOW,                                  &
     &           .TRUE.,                                                &
     &           KRAIN_DUST(3), KSNOW_DUST(3),                          &
     &           CONSCAV_DUST(1,1,3)                                    &
     &           )
!
! DEPENDS ON: scnscv2
          CALL SCNSCV2(                                                 &
     &           ROW_LENGTH, ROWS,                                      &
     &           OFF_X, OFF_Y, HALO_I, HALO_J,                          &
     &           MODEL_LEVELS, WET_MODEL_LEVELS,                        &
     &           R_RHO_LEVELS, R_THETA_LEVELS,                          &
     &           TIMESTEP,SUBSTEP_NUMBER,                               &
     &           RHO,                                                   &
     &           Q_N, QCL_N, QCF_N, DUST_DIV4,                          &
     &           CCB, CCT,                                              &
     &           CONV_RAIN, CONV_SNOW,                                  &
     &           .TRUE.,                                                &
     &           KRAIN_DUST(4), KSNOW_DUST(4),                          &
     &           CONSCAV_DUST(1,1,4)                                    &
     &           )
!
! DEPENDS ON: scnscv2
          CALL SCNSCV2(                                                 &
     &           ROW_LENGTH, ROWS,                                      &
     &           OFF_X, OFF_Y, HALO_I, HALO_J,                          &
     &           MODEL_LEVELS, WET_MODEL_LEVELS,                        &
     &           R_RHO_LEVELS, R_THETA_LEVELS,                          &
     &           TIMESTEP,SUBSTEP_NUMBER,                               &
     &           RHO,                                                   &
     &           Q_N, QCL_N, QCF_N, DUST_DIV5,                          &
     &           CCB, CCT,                                              &
     &           CONV_RAIN, CONV_SNOW,                                  &
     &           .TRUE.,                                                &
     &           KRAIN_DUST(5), KSNOW_DUST(5),                          &
     &           CONSCAV_DUST(1,1,5)                                    &
     &           )
!
! DEPENDS ON: scnscv2
          CALL SCNSCV2(                                                 &
     &           ROW_LENGTH, ROWS,                                      &
     &           OFF_X, OFF_Y, HALO_I, HALO_J,                          &
     &           MODEL_LEVELS, WET_MODEL_LEVELS,                        &
     &           R_RHO_LEVELS, R_THETA_LEVELS,                          &
     &           TIMESTEP,SUBSTEP_NUMBER,                               &
     &           RHO,                                                   &
     &           Q_N, QCL_N, QCF_N, DUST_DIV6,                          &
     &           CCB, CCT,                                              &
     &           CONV_RAIN, CONV_SNOW,                                  &
     &           .TRUE.,                                                &
     &           KRAIN_DUST(6), KSNOW_DUST(6),                          &
     &           CONSCAV_DUST(1,1,6)                                    &
     &           )
!
        ENDIF !L_DUST

! Scavenge Sulphur Cycle tracers
!
        If (L_sulpc_so2) Then
!
! Scavenge SO2
! DEPENDS ON: scnwsh2
          Call scnwsh2(                                                 &
     &           row_length, rows,                                      &
     &           off_x, off_y, halo_i, halo_j,                          &
     &           model_levels, wet_model_levels,                        &
     &           r_rho_levels, r_theta_levels,                          &
     &           timestep, substep_number,                              &
     &           rho,                                                   &
     &           q_n, qcl_n, qcf_n, so2,                                &
     &           ccb, cct,                                              &
     &           conv_rain, conwash_so2                                 &
     &           )
!
! Scavenge NH3 if present
!
            If (L_sulpc_nh3) Then
!
! Scavenge NH3
! DEPENDS ON: ncnwsh2
            Call ncnwsh2(                                               &
     &           row_length, rows,                                      &
     &           off_x, off_y, halo_i, halo_j,                          &
     &           model_levels, wet_model_levels,                        &
     &           r_rho_levels, r_theta_levels,                          &
     &           timestep, substep_number,                              &
     &           rho,                                                   &
     &           q_n, qcl_n, qcf_n, nh3,                                &
     &           ccb, cct,                                              &
     &           conv_rain, conwash_nh3                                 &
     &           )
!
!
            End If              !End L_sulpc_nh3 condition
!
! Scavenge SO4_AIT
! DEPENDS ON: scnscv2
            Call scnscv2(                                               &
     &           row_length, rows,                                      &
     &           off_x, off_y, halo_i, halo_j,                          &
     &           model_levels, wet_model_levels,                        &
     &           r_rho_levels, r_theta_levels,                          &
     &           timestep, substep_number,                              &
     &           rho,                                                   &
     &           q_n, qcl_n, qcf_n, so4_aitken,                         &
     &           ccb, cct,                                              &
     &           conv_rain, conv_snow,                                  &
     &           .false.,                                               &
     &           krain_so4ait, ksnow_so4ait,                            &
     &           conscav_so4ait                                         &
     &           )
!
! Scavenge SO4_ACC
! DEPENDS ON: scnscv2
            Call scnscv2(                                               &
     &           row_length, rows,                                      &
     &           off_x, off_y, halo_i, halo_j,                          &
     &           model_levels, wet_model_levels,                        &
     &           r_rho_levels, r_theta_levels,                          &
     &           timestep, substep_number,                              &
     &           rho,                                                   &
     &           q_n, qcl_n, qcf_n, so4_accu,                           &
     &           ccb, cct,                                              &
     &           conv_rain, conv_snow,                                  &
     &           .false.,                                               &
     &           krain_so4acc, ksnow_so4acc,                            &
     &           conscav_so4acc                                         &
     &           )
!
! Scavenge SO4_DIS
! DEPENDS ON: scnscv2
            Call scnscv2(                                               &
     &           row_length, rows,                                      &
     &           off_x, off_y, halo_i, halo_j,                          &
     &           model_levels, wet_model_levels,                        &
     &           r_rho_levels, r_theta_levels,                          &
     &           timestep, substep_number,                              &
     &           rho,                                                   &
     &           q_n, qcl_n, qcf_n, so4_diss,                           &
     &           ccb, cct,                                              &
     &           conv_rain, conv_snow,                                  &
     &           .false.,                                               &
     &           krain_so4dis, ksnow_so4dis,                            &
     &           conscav_so4dis                                         &
     &           )
!
        End If              !End L_sulpc_so2 condition
!
!
        If (l_soot) then  ! If soot modelling is included
!
! Scavenge aged soot
! DEPENDS ON: scnscv2
          Call scnscv2(                                                 &
     &         row_length, rows,                                        &
     &         off_x, off_y, halo_i, halo_j,                            &
     &         model_levels, wet_model_levels,                          &
     &         r_rho_levels, r_theta_levels,                            &
     &         timestep, substep_number,                                &
     &         rho,                                                     &
     &         q_n, qcl_n, qcf_n, soot_agd,                             &
     &         ccb, cct,                                                &
     &         conv_rain, conv_snow,                                    &
     &         .false.,                                                 &
     &         krain_agedsoot, ksnow_agedsoot,                          &
     &         conscav_agedsoot                                         &
     &         )
!
        End If  ! l_soot
!
        If (l_biomass) then  ! If biomass aerosol is included
!
! Scavenge aged biomass aerosol
! DEPENDS ON: scnscv2
          Call scnscv2(                                                 &
     &         row_length, rows,                                        &
     &         off_x, off_y, halo_i, halo_j,                            &
     &         model_levels, wet_model_levels,                          &
     &         r_rho_levels, r_theta_levels,                            &
     &         timestep, substep_number,                                &
     &         rho,                                                     &
     &         q_n, qcl_n, qcf_n, bmass_agd,                            &
     &         ccb, cct,                                                &
     &         conv_rain, conv_snow,                                    &
     &         .false.,                                                 &
     &         krain_agedbmass, ksnow_agedbmass,                        &
     &         conscav_agedbmass                                        &
     &         )
!
        End If  ! l_biomass
!
        If (l_ocff) then ! If fossil-fuel OC aerosol is included
!
! Scavenge aged fossil-fuel organic carbon aerosol
! DEPENDS ON: scnscv2
          Call scnscv2(                                                 &
     &         row_length, rows,                                        &
     &         off_x, off_y, halo_i, halo_j,                            &
     &         model_levels, wet_model_levels,                          &
     &         r_rho_levels, r_theta_levels,                            &
     &         timestep, substep_number,                                &
     &         rho,                                                     &
     &         q_n, qcl_n, qcf_n, ocff_agd,                             &
     &         ccb, cct,                                                &
     &         conv_rain, conv_snow,                                    &
     &         .false.,                                                 &
     &         krain_agedocff, ksnow_agedocff,                          &
     &         conscav_agedocff                                         &
     &         )
!
        End If  ! l_ocff
!

      End If ! CycleNo == NumCycles
!
! ----------------------------------------------------------------------
! Section CNV 5.3 Calculate PC2 scheme increments to the increments.
! ----------------------------------------------------------------------
!
! L_calc_dxek_if2:
      If (L_calc_dxek) Then
!
          Allocate ( T_inc_latest(row_length,rows,1) )
!
          Allocate ( theta_inc_PC2(row_length,rows,1) )
          Allocate ( q_inc_PC2(row_length,rows,1) )
          Allocate ( qcl_inc_PC2(row_length,rows,1) )
          Allocate ( cfl_inc_PC2(row_length,rows,1) )
          Allocate ( bcf_inc_PC2(row_length,rows,1) )
!
          Allocate ( T_earliest(row_length,rows,1) )
          Allocate ( q_earliest(row_length,rows,1) )
          Allocate ( qcl_earliest(row_length,rows,1) )
          Allocate ( cfl_earliest(row_length,rows,1) )
          Allocate ( cff_earliest(row_length,rows,1) )
          Allocate ( bcf_earliest(row_length,rows,1) )
!
! Convert potential temperature increment to temperature increment
!
        Do k=1, n_conv_levels
          Do j=1, rows
            Do i=1, row_length
              T_inc_latest(i,j,1) =                                     &
     &          theta_inc(i,j,k) * exner_theta_levels(i,j,k)
              T_earliest(i,j,1)                                         &
     &          = theta_n(i,j,k)*exner_theta_levels(i,j,k)              &
     &          - T_inc_latest(i,j,1)
              q_earliest(i,j,1) = q_n(i,j,k) - q_inc(i,j,k)
              qcl_earliest(i,j,1) = qcl_n(i,j,k) - qcl_inc(i,j,k)
              bcf_earliest(i,j,1) = bulk_cf_n(i,j,k)                    &
     &                            - bulk_cf_inc(i,j,k)
              cfl_earliest(i,j,1) = cf_liquid_n(i,j,k)                  &
     &                            - cf_liquid_inc(i,j,k)
              cff_earliest(i,j,1) = cf_frozen_n(i,j,k)                  &
     &                            - cf_frozen_inc(i,j,k)
            End Do ! i
          End Do ! j

          If (lcv_pc2_diag_sh) then
         
           ! call pc2_hom_conv without erosion term i.e. dbsdtbs_turb_0=0.0
           ! and dbsdtbs_turb_1=0.0 
!
! DEPENDS ON: pc2_hom_conv
            Call PC2_HOM_CONV(                                          &
! Input variables
             p_layer_centres(1,1,k), 1, row_length, rows                &
       ,     timestep                                                   &
! INput variables
       ,     T_earliest(1,1,1), q_earliest(1,1,1), qcl_earliest(1,1,1)  &
       ,     bcf_earliest(1,1,1),cfl_earliest(1,1,1),cff_earliest(1,1,1)&
       ,     T_inc_latest,q_inc(1,1,k),full_zero(1,1,1),full_zero(1,1,1)&
       ,     cf_liquid_incr_inhom_diag(1,1,k)                           &
! OUTput variables
       ,     theta_inc_PC2,q_inc_PC2,qcl_inc_PC2,bcf_inc_PC2,cfl_inc_PC2&
! INput variables (other quantities)
       ,             0.0,              0.0 , l_mixing_ratio             &
                                 )
          Else
           ! call pc2_hom_conv with erosion term   
!
! DEPENDS ON: pc2_hom_conv
            Call PC2_HOM_CONV(                                          &
! Input variables
             p_layer_centres(1,1,k), 1, row_length, rows                &
       ,     timestep                                                   &
! INput variables
       ,     T_earliest(1,1,1), q_earliest(1,1,1), qcl_earliest(1,1,1)  &
       ,     bcf_earliest(1,1,1),cfl_earliest(1,1,1),cff_earliest(1,1,1)&
       ,     T_inc_latest,q_inc(1,1,k),full_zero(1,1,1),full_zero(1,1,1)&
       ,     cf_liquid_incr_inhom_diag(1,1,k)                           &
! OUTput variables
       ,     theta_inc_PC2,q_inc_PC2,qcl_inc_PC2,bcf_inc_PC2,cfl_inc_PC2&
! INput variables (other quantities)
       ,     dbsdtbs_turb_0, dbsdtbs_turb_1, l_mixing_ratio             &
                                 )
          End If 

!
! Calculate potential temperature increment (on convect levels only)
! from temperature increment output by PC2_Homog.
!
          Do j=1,rows
            Do i=1,row_length
              theta_inc_PC2(i,j,1) = theta_inc_PC2(i,j,1) /             &
     &                               exner_theta_levels(i,j,k)
            End Do ! i
          End Do ! j
!
! L_q_interact_if1:
          If (L_q_interact)  Then
!
            Do j = 1, rows
              Do i = 1, row_length
!
! Update increments to theta, moisture and cloud fields with additional
! increments from the response to environment changes (homogenous),
!
                theta_inc(i,j,k) = theta_inc(i,j,k)                     &
     &                           + theta_inc_PC2(i,j,1)
                q_inc(i,j,k)   =   q_inc(i,j,k) +   q_inc_PC2(i,j,1)
                qcl_inc(i,j,k) = qcl_inc(i,j,k) + qcl_inc_PC2(i,j,1)
! Not updated   qcf_inc(i,j,k) = qcf_inc(i,j,k)
                cf_liquid_inc(i,j,k) = cf_liquid_inc(i,j,k) +           &
     &                                 cfl_inc_PC2(i,j,1)
! Not updated   cf_frozen_inc(i,j,k) = cf_frozen_inc(i,j,k)
                bulk_cf_inc(i,j,k)   = bulk_cf_inc(i,j,k) +             &
     &                                 bcf_inc_PC2(i,j,1)
!
! ... and update working version of theta, moisture and cloud fields.
!
                theta_n(i,j,k) = theta_n(i,j,k) + theta_inc_PC2(i,j,1)
                q_n(i,j,k)     =     q_n(i,j,k) +     q_inc_PC2(i,j,1)
                qcl_n(i,j,k)   =   qcl_n(i,j,k) +   qcl_inc_PC2(i,j,1)
! Not updated   qcf_n(i,j,k)   =   qcf_n(i,j,k)
                cf_liquid_n(i,j,k) = cf_liquid_n(i,j,k) +               &
     &                                 cfl_inc_PC2(i,j,1)
! Not updated   cf_frozen_n(i,j,k) = cf_frozen_n(i,j,k)
                bulk_cf_n(i,j,k)   = bulk_cf_n(i,j,k) +                 &
     &                                 bcf_inc_PC2(i,j,1)
!
! ... and the diagnostics.
!
                T_incr_diagnostic(i,j,k) = T_incr_diagnostic(i,j,k)     &
     &           + exner_theta_levels(i,j,k) * theta_inc_PC2(i,j,1)
                q_incr_diagnostic(i,j,k) = q_incr_diagnostic(i,j,k)     &
     &                                   + q_inc_PC2(i,j,1)
                qcl_incr_diagnostic(i,j,k)= qcl_incr_diagnostic(i,j,k)  &
     &                                    + qcl_inc_PC2(i,j,1)
! Not updated   qcf_incr_diagnostic(i,j,k)= qcf_incr_diagnostic(i,j,k)
                cf_liquid_incr_diagnostic(i,j,k)                        &
     &        = cf_liquid_incr_diagnostic(i,j,k) + cfl_inc_PC2(i,j,1)
! Not updated   cf_frozen_incr_diagnostic(i,j,k)
!    &        = cf_frozen_incr_diagnostic(i,j,k)
                bulk_cf_incr_diagnostic(i,j,k)                          &
     &        = bulk_cf_incr_diagnostic(i,j,k)   + bcf_inc_PC2(i,j,1)
!
              End Do  ! i loop
            End Do  ! j
!
          End If  ! L_q_interact_if1

        End Do  ! k

      End If  ! L_calc_dxek_if2
!

      Select Case (Rad_cloud_decay_opt)

        Case(rad_decay_full_timestep)
          !-------------------------------------------------------------------
          ! Use Cloud decay on convective cloud. Take maximum cloud parameters
          ! when comparing this timestep with those of previous full timstep
          ! (allowing for cloud decay) 
          !-------------------------------------------------------------------

          Cloud_lifetime = fixed_cld_life 

          Do j=1, rows
            Do i=1, row_length
              ccb_rad(i,j) = 0
              cct_rad(i,j) = 0
            End Do
          End Do

          Do k=1, n_cca_levels
            Do j=1, rows
              Do i=1, row_length
                cca_rad(i,j,k) = MAX( cca_rad(i,j,k)                          &
                                     *(1.0-timestep/cloud_lifetime)           &
                                    , cca(i,j,k))

                ! Apply threshold for minimum CCA
                If (cca_rad(i,j,k) <  CCA_min) Then 
                  cca_rad(i,j,k) = 0.0
                End If

                If ((ccb_rad(i,j)   == 0  ) .AND.                             &
                    (cca_rad(i,j,k) >  0.0)) Then
                  ccb_rad(i,j) = k
                End If

                If (cca_rad(i,j,k) > 0.0) Then
                  cct_rad(i,j) = k+1
                End If   

              End Do
            End Do
          End Do

          Do j=1, rows
            Do i=1, row_length
              If (ccb_rad(i,j) == 0) Then
                cclwp_rad(i,j) = 0.0
              Else
                cclwp_rad(i,j) = MAX( cclwp_rad(i,j)                          &
                                     *(1.0-timestep/Cloud_lifetime)           &
                                    , cclwp(i,j))
              End If
            End Do
          End Do

        Case(rad_decay_conv_substep)
          !-------------------------------------------------------------------
          ! Use Cloud decay on convection substeps(ccrad). See earlier(2.0);
          ! at end of convection substepping loop 
          !-------------------------------------------------------------------

        Case(rad_decay_off)
          !-------------------------------------------------------------------
          ! Rad cloud decay disabled (neither == 1 or 2):
          ! No time-decay of convective cloud so simply pass out this
          ! timesteps parameters
          !-------------------------------------------------------------------

          Do j=1, rows
            Do i=1, row_length
              ccb_rad   (i,j) = ccb   (i,j)
              cct_rad   (i,j) = cct   (i,j)
              cclwp_rad (i,j) = cclwp (i,j)
              lcbase_rad(i,j) = lcbase(i,j)
              lctop_rad (i,j) = lctop (i,j)
            End Do
          End Do

          Do k=1, n_cca_levels
            Do j=1, rows
              Do i=1, row_length
                cca_rad(i,j,k) = cca(i,j,k)
              End Do
            End Do
          End Do

          Do k=1, wet_model_levels
            Do j=1, rows
              Do i=1, row_length
                ccw_rad(i,j,k) = ccw(i,j,k)
              End Do
            End Do
          End Do

      End Select  ! Rad_cloud_decay_opt


! DEPENDS ON: timer
        If (Ltimer) Call timer ('Convect ',4)
! ----------------------------------------------------------------------
! Section CNV.6 Call Convection diagnostics, also returns total
!               precipitation diagnostics.
! ----------------------------------------------------------------------

! Check that convection diagnostics requested this timestep
#if !defined(SCMA)
! diagnostics requested this timestep
        If ( sf(0,5) .and. Substep_Number  ==  Num_Substeps ) Then
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',3)

      If ( L_apply_diag ) Then

! DEPENDS ON: diagnostics_conv
        Call diagnostics_conv(                                          &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels                           &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      u, v, p, R_u, R_v                          &
     &,                      qcl_incr_inhom_diag, qcf_incr_inhom_diag   &
     &,                      bulk_cf_incr_inhom_diag                    &
     &,          cf_liquid_incr_inhom_diag, cf_frozen_incr_inhom_diag   &
     &,                      Theta_diag, theta_inc, q_diag, q_inc       &
     &,                      qcl_inc, qcf_inc, cf_liquid_inc            &
     &,                      cf_frozen_inc, bulk_cf_inc                 &
     &,                      T_incr_diagnostic,q_incr_diagnostic        &
     &,                      qcl_incr_diagnostic, qcf_incr_diagnostic   &
     &,                      u_incr_diagnostic,v_incr_diagnostic        &
     &,          cf_liquid_incr_diagnostic, cf_frozen_incr_diagnostic   &
     &,                      bulk_cf_incr_diagnostic                    &
     &,                      exner_theta_levels                         &
     &,                      ls_rain, ls_snow                           &
     &,                      ccw, conv_rain, conv_snow                  &
     &,                      conv_rain_3d, conv_snow_3d                 &
     &,                      cca, ccb, cct                              &
     &,                      dubydt_pout,dvbydt_pout                    &
     &,                      up_flux_half                               &
     &,                      up_flux,dwn_flux,entrain_up,detrain_up     &
     &,                      entrain_dwn,detrain_dwn,uw_dp,vw_dp        &
     &,                      uw_shall,vw_shall                          &
     &,              wqt_flux,wthetal_flux,wthetav_flux,wql_flux        &
     &,              wqt_cb,wthetal_cb,wqt_inv,wthetal_inv,sh_top       &
     &,              sh_base,shallow_ind,congestus_ind,congestus_ind2   &
     &,              cu_over_orog,cape_out,cape_undilute,cin_undilute   &
     &,              mid_ind                                            &
     &,              ntml_diag,ntpar_diag,freeze_diag,kterm_deep        &
     &,              precip_deep,precip_shall,precip_mid,precip_cong    &
     &,              wstar_dn_diag,wstar_up_diag,mb1_diag,mb2_diag      &
     &,              cg_term, cg_top,cg_base                            &
     &,              mf_deep,mf_congest,mf_shall,mf_midlev              &
     &,              dt_deep,dt_congest,dt_shall,dt_midlev              &
     &,              dq_deep,dq_congest,dq_shall,dq_midlev              &
     &,              du_deep,du_congest,du_shall,du_midlev              &
     &,              dv_deep,dv_congest,dv_shall,dv_midlev              &
     &,                      lcbase, lctop, lcca, n_cca_levels          &
     &,                      CONSCAV_DUST(1,1,1),CONSCAV_DUST(1,1,2)    &
     &,                      CONSCAV_DUST(1,1,3),CONSCAV_DUST(1,1,4)    &
     &,                      CONSCAV_DUST(1,1,5),CONSCAV_DUST(1,1,6)    &
     &,                      L_DUST                                     &
     &,                      conwash_nh3                                &
     &,                      conwash_so2, conscav_so4ait                &
     &,                      conscav_so4acc, conscav_so4dis             &
     &,                      conscav_agedsoot, conscav_agedbmass        &
     &,                      conscav_agedocff                           &
     &,                      timestep, Num_Substeps                     &
     &       ,                                                          &
#include "argsts.h"
     & STASHwork                                                        &
     & )

      End If !  L_apply_diag
!
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',4)
        End If            ! on sf(0,5)
#else
! SCM diagnostics
!
!-----------------------------------------------------------------------
!     SCM Convection OR Increments Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_conv)                                      &
     &      .OR. L_SCMDiags(SCMDiag_incs)) Then

!       Stash 5,163
! DEPENDS ON: scmoutput
        Call SCMoutput(qcl_incr_inhom_diag,                             &
             'dqcl_cinh','QCL inc: conv inhom','kg/kg',                 &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 5,164
! DEPENDS ON: scmoutput
        Call SCMoutput(qcf_incr_inhom_diag,                             &
             'dqcf_cinh','QCF inc: conv inhom','kg/kg',                 &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 5,172
! DEPENDS ON: scmoutput
        Call SCMoutput(bulk_cf_incr_inhom_diag,                         &
             'dbcf_cinh','Bulk cloud frac inc: conv inhom',             &
             'Fraction',                                                &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 5,173
! DEPENDS ON: scmoutput
        Call SCMoutput(cf_liquid_incr_inhom_diag,                       &
             'dcfl_cinh','Liquid cloud frac inc: conv inhom',           &
             'Fraction',                                                &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 5,174
! DEPENDS ON: scmoutput
        Call SCMoutput(cf_frozen_incr_inhom_diag,                       &
             'dcff_cinh','Frozen cloud frac inc: conv inhom',           &
             'Fraction',                                                &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 5,181
! DEPENDS ON: scmoutput
        Call SCMoutput(T_incr_diagnostic,                               &
             'dt_conv','T increment convection','K',                    &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,182
! DEPENDS ON: scmoutput
        Call SCMoutput(q_incr_diagnostic,                               &
             'dq_conv','q increment convection','kg/kg',                &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 5,183
! DEPENDS ON: scmoutput
        Call SCMoutput(qcl_incr_diagnostic,                             &
             'dqcl_conv','QCL increment convection','kg/kg',            &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 5,184
! DEPENDS ON: scmoutput
        Call SCMoutput(qcf_incr_diagnostic,                             &
             'dqcf_conv','QCF increment convection','kg/kg',            &
             t_avg,d_wet,default_streams,'',RoutineName)

!      Stash 5,185
! DEPENDS ON: scmoutput
        Call SCMoutput(u_incr_diagnostic,                               &
             'du_conv','U increment convection','m/s',                  &
             t_avg,d_all,default_streams,'',RoutineName)

!      Stash 5,186
! DEPENDS ON: scmoutput
        Call SCMoutput(v_incr_diagnostic,                               &
             'dv_conv','V increment convection','m/s',                  &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,192
! DEPENDS ON: scmoutput
        Call SCMoutput(bulk_cf_incr_diagnostic,                         &
             'dbcf_conv',                                               &
             'Bulk cloud frac increment convection','Fraction',         &
             t_avg,d_wet,default_streams,'',RoutineName)

!       Stash 5,193
! DEPENDS ON: scmoutput
        Call SCMoutput(cf_liquid_incr_diagnostic,                       &
             'dcfl_conv','Liquid cloud frac increment convection',      &
             'Fraction',                                                &
             t_avg,d_wet,default_streams,'',RoutineName)

!      Stash 5,194
! DEPENDS ON: scmoutput
        Call SCMoutput(cf_frozen_incr_diagnostic,                       &
             'dcff_conv','Frozen cloud frac increment convection',      &
             'Fraction',                                                &
             t_avg,d_wet,default_streams,'',RoutineName)
!
      End If ! L_SCMDiags(SCMDiag_conv) .OR. L_SCMDiags(SCMDiag_incs)
!
!-----------------------------------------------------------------------
!     SCM Convection Diagnostics Package
!-----------------------------------------------------------------------
      If (L_SCMDiags(SCMDiag_conv)) Then

!       Stash 5,209
        Do k = 1, model_levels
          Do j = 1, rows
             Do i = 1, row_length
               work_3d(i,j,k) = (theta_diag(i,j,k) +                    &
     &               theta_inc(i,j,k) )* exner_theta_levels(i,j,k)
             End Do ! i
          End Do ! j
        End Do ! k

! DEPENDS ON: scmoutput
        Call SCMoutput(work_3d,                                         &
             'T_afterconv','Temperature after convection','K',          &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,273
        work_2d=float(ntml)
! DEPENDS ON: scmoutput
        Call SCMoutput(work_2d,                                         &
             'NTML','Top level of surface mixed layer',                 &
             'Model level',                                             &
             t_inst,d_point,default_streams,'',RoutineName)

!       Stash 5,274
        work_2d=float(ntpar)
! DEPENDS ON: scmoutput
        Call SCMoutput(work_2d,                                         &
             'NTPAR','Top level of initial parcel ascent',              &
             'Model level',                                             &
             t_inst,d_point,default_streams,'',RoutineName)

!       Stash 5,275
        work_2d=float(freeze_lev)
! DEPENDS ON: scmoutput
        Call SCMoutput(work_2d,                                         &
             'freeze_lev','Freezing level','Model level',               &
             t_inst,d_point,default_streams,'',RoutineName)
!
        Do j=1,rows
          Do i=1,row_length
            If (cumulus(i,j)) Then
              work_2d(i,j)=1.0
            Else
              work_2d(i,j)=0.0
            End If
          End Do ! i
        End Do ! j

! DEPENDS ON: scmoutput
        Call SCMoutput(work_2d,                                         &
             'ind_cumulus','Indicator for cumulus convection',          &
             'Indicator',                                               &
             t_inst,d_point,default_streams,'',RoutineName)

!       Stash 5,270
        Do j=1,rows
          Do i=1,row_length
            If (l_shallow(i,j)) Then
              work_2d(i,j)=1.0
            Else
              work_2d(i,j)=0.0
            End If
          End Do ! i
        End Do ! j

! DEPENDS ON: scmoutput
        Call SCMoutput(work_2d,                                         &
             'ind_shallow','Indicator for shallow convection',          &
             'Indicator',                                               &
             t_inst,d_point,default_streams,'',RoutineName)

!       Stash 5,272
        Do j=1,rows
          Do i=1,row_length
            If (it_mid_level(i,j)) Then
              work_2d(i,j)=1.0
            Else
              work_2d(i,j)=0.0
            End If
          End Do ! i
        End Do ! j

! DEPENDS ON: scmoutput
        Call SCMoutput(work_2d,                                         &
             'ind_midconv','Indicator for mid-level convection',        &
             'Indicator',                                               &
             t_inst,d_point,default_streams,'',RoutineName)

        If (iconv_congestus == 1) Then  
!        Stash 5,310
          Do j=1,rows
            Do i=1,row_length
              If (l_congestus(i,j)) Then
                work_2d(i,j)=1.0
              Else
                work_2d(i,j)=0.0
              End If
            End Do ! i
          End Do ! j

! DEPENDS ON: scmoutput
          Call SCMoutput(work_2d,                                       & 
              'ind_congest','Indicator for congestus convection',       &
              'Indicator',                                              &
              t_inst,d_point,default_streams,'',RoutineName)
        End If

!       Stash 5,276
        work_2d=kterm_deep
! DEPENDS ON: scmoutput
        Call SCMoutput(work_2d,                                         &
             'kterm_deep','Deep convection termination level',          &
             'Model level',                                             &
             t_inst,d_point,default_streams,'',RoutineName)

!       Stash 5,277
! DEPENDS ON: scmoutput
        Call SCMoutput(precip_deep,                                     &
             'precip_deep','Deep convective precipitation',             &
             'kg/m2/s',                                                 &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 5,278
! DEPENDS ON: scmoutput
        Call SCMoutput(precip_shall,                                    &
             'precip_shall','Shallow convective precipitation',         &
             'kg/m2/s',                                                 &
             t_avg,d_sl,default_streams,'',RoutineName)

!       Stash 5,279
! DEPENDS ON: scmoutput
        Call SCMoutput(precip_mid,                                      &
             'precip_mid','Mid level convective precipitation',         &
             'kg/m2/s',                                                 &
             t_avg,d_sl,default_streams,'',RoutineName)

        If (iconv_congestus == 1) Then  
!        Stash 5,280
! DEPENDS ON: scmoutput
          Call SCMoutput(precip_cong,                                   &
             'precip_cong','Congestus convective precipitation',        &
             'kg/m2/s',                                                 &
             t_avg,d_sl,default_streams,'',RoutineName) 
        End If
!       Stash 5,249
! DEPENDS ON: scmoutput
        Call SCMoutput(up_flux_half,                                    &
             'up_halfmassflux','updraught mass flux on half levs',      &
             'Pa/s',                                                    &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,250
! DEPENDS ON: scmoutput
        Call SCMoutput(up_flux,                                         &
             'up_massflux','updraught mass flux','Pa/s',                &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,251
! DEPENDS ON: scmoutput
        Call SCMoutput(dwn_flux,                                        &
             'down_massflux','downdraught mass flux','Pa/s',            &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,252
! DEPENDS ON: scmoutput
        Call SCMoutput(entrain_up,                                      &
             'entrain_up','updraught entrainment rate','s-1',           &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,254
! DEPENDS ON: scmoutput
        Call SCMoutput(entrain_dwn,                                     &
             'entrain_dw','downdraught entrainment rate','s-1',         &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,253
! DEPENDS ON: scmoutput
        Call SCMoutput(detrain_up,                                      &
             'detrain_up','updraught detrainment rate','s-1',           &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,255
! DEPENDS ON: scmoutput
        Call SCMoutput(detrain_dwn,                                     &
             'detrain_dw','downdraught detrainment rate','s-1',         &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,258
! DEPENDS ON: scmoutput
        Call SCMoutput(uw_dp,                                           &
             'uw_dp','x_comp of stress from deep convection',           &
             'kg/m/s2',                                                 &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,259
! DEPENDS ON: scmoutput
        Call SCMoutput(vw_dp,                                           &
             'vw_dp','y_comp of stress from deep convection',           &
             'kg/m/s2',                                                 &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,260
! DEPENDS ON: scmoutput
        Call SCMoutput(uw_shall,                                        &
             'uw_shall','x_comp of stress from shallow convection',     &
             'kg/m/s2',                                                 &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,261
! DEPENDS ON: scmoutput
        Call SCMoutput(vw_shall,                                        &
             'vw_shall','y_comp of stress from shallow convection',     &
             'kg/m/s2',                                                 &
             t_avg,d_all,default_streams,'',RoutineName)

!       Stash 5,227 
! DEPENDS ON: scmoutput
        Call SCMoutput(conv_rain_3d, 'conv_rain_3d',                    &
             'Convective Rainfall Flux',                                &
             'kg/m2/s',t_avg,d_wet,default_streams,'',RoutineName)  

!       Stash 5,228
! DEPENDS ON: scmoutput
        Call SCMoutput(conv_snow_3d, 'conv_snow_3d',                    &
             'Convective Snowfall Flux',                                &
             'kg/m2/s',t_avg,d_wet,default_streams,'',RoutineName)  

!       Stash 5,213
! DEPENDS ON: scmoutput
        Call SCMoutput(ccw,                                         &
     &       'ccw','Convective cloud water','kg/kg',                &
     &       t_avg,d_all,default_streams,'')


! Extra added for testing

! DEPENDS ON: scmoutput
        Call SCMoutput(wstar,                                           &
             'wstar','Sub cloud convective velocity scale',             &
             'm/s',                                                     &
             t_avg,d_sl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
        Call SCMoutput(wthvs,                                           &
             'wthvs','wthetav flux at surface',                         &
             'K m/s',                                                   &
             t_avg,d_sl,default_streams,'',RoutineName)

        ! Stash 5,212
        ! Convective Cloud Amount
! DEPENDS ON: scmoutput
        Call SCMoutput(cca,                                             &
             'cca','Convective Cloud Amount', 'Fraction',               &
             t_avg,d_all,default_streams,'',RoutineName)

        ! Stash 5,213
        ! Convective Cloud Water
! DEPENDS ON: scmoutput
        Call SCMoutput(ccw,                                             &
             'ccw','Convective Cloud Water', 'kg/kg',                   &
             t_avg,d_wet,default_streams,'',RoutineName)

        ! Condensed water path convective cloud
! DEPENDS ON: scmoutput
        Call SCMoutput(cclwp,                                           &
             'cclwp','Condensed Cloud Water Path','Kg/m2',              &
             t_avg,d_sl,default_streams,'',RoutineName)



      End If ! L_SCMDiags(SCMDiag_conv)
#endif


! Clear up allocatable arrays
!
      If (L_full_zero) Then
        Deallocate ( full_zero )
      End If
!
      If (L_calc_dxek) Then
        Deallocate ( T_earliest )
        Deallocate ( q_earliest )
        Deallocate ( qcl_earliest )
        Deallocate ( cfl_earliest )
        Deallocate ( cff_earliest )
        Deallocate ( bcf_earliest )
        Deallocate ( T_inc_latest )
!
        Deallocate ( theta_inc_PC2 )
        Deallocate ( q_inc_PC2 )
        Deallocate ( qcl_inc_PC2 )
        Deallocate ( cfl_inc_PC2 )
        Deallocate ( bcf_inc_PC2 )
      End If
!
      If(l_qcl_incr_cinh) Then
         Deallocate (qcl_incr_inhom_diag)
      Endif
!
      If(l_qcf_incr_cinh) Then
         Deallocate (qcf_incr_inhom_diag)
      Endif
!
      If(l_bcf_incr_cinh) Then
         Deallocate (bulk_cf_incr_inhom_diag)
      Endif
!
      If(l_cfl_incr_cinh .or. L_calc_dxek) Then
         Deallocate (cf_liquid_incr_inhom_diag)
      Endif
!
      If(l_cff_incr_cinh) Then
         Deallocate (cf_frozen_incr_inhom_diag)
      Endif
!
      End If ! on error code equal to zero

! end of routine NI_conv_ctl
      Return
      END SUBROUTINE NI_conv_ctl
#endif
#endif
