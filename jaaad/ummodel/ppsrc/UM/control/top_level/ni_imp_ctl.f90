

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine NI_imp_ctl
! *********************************************************************
!
!+ Boundary layer Implicit solver and Large-scale (Area) Cloud Scheme.
! Subroutine Interface:
      Subroutine NI_imp_ctl (                                           &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &

! IN Substepping information
     &, Substep_Number, Num_Substeps, L_phys2_substep                   &

! model dimensions.
     &, row_length, rows, rhc_row_length, rhc_rows, n_rows, land_points &
     &, ntiles, model_levels, wet_model_levels, bl_levels, dst_levels   &
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
     &, NumCycles, CycleNo, lq_mix_bl, L_ukca, L_sice_heatflux          &
     &, L_subfilter_vert, L_subfilter_horiz, l_use_cariolle             &
!
! model Parameters
     &, alpha_cd, Puns, Pstb, trweights1, rhcrit, tr_vars               &
     &, cloud_fraction_method,overlap_ice_liquid,ice_fraction_method    &
     &, ctt_weight,t_weight,qsat_fixed,sub_cld,x1i,x1ic,x1r,x2r,x4r     &
     &, l_psd, ai, bi, aic, bic, lsp_ei, lsp_fi, lsp_eic, lsp_fic       &

! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &

! in coordinate information
     &, r_rho_levels, r_theta_levels, r_at_u, r_at_v, eta_theta_levels  &
     &, eta_rho_levels, delta_lambda, delta_phi,lat_rot_NP,long_rot_NP  &

! in time stepping information.
     &, timestep, val_year, val_day_number, val_hour, val_minute        &
     &, val_second, timestep_number                                     &

! trig arrays
     &, sin_theta_longitude, cos_theta_longitude, FV_cos_theta_latitude &

! diagnostic info
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork3, STASHwork9                                           &
!
! SCM diagnostics (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! in data fields.
     &, p, p_layer_centres, p_layer_boundaries, rho, u, v, w            &
     &, land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels   &
     &, u_conv,v_conv,T_conv,q_conv,qcl_conv,qcf_conv                   &
! ancillary fields and fields needed to be kept from timestep to
! timestep
     &, smvccl, smvcwt, smvcst, sthf, sthu, sil_orog_land, ho2r2_orog   &
     &, di, ice_fract, di_ncat, ice_fract_ncat                          &
     &, u_0, v_0,land_index, cca, ccb, cct, ccw, surf_radflux           &
     &, ls_rain, ls_snow, conv_rain, conv_snow, cca_2d, L_scrn, L_plsp  &

! in variables required from BDY_LAYR in IMP_SOLVER
! in variables required from BDY_LAYR in IMP_SOLVER
     &, alpha1_sice, ashtf, dtrdz_charney_grid, rdz_charney_grid        &
     &, dtrdz_u, dtrdz_v, rdz_u, rdz_v, cdr10m_u, cdr10m_v, z1_tq       &
!ajm    extra variable added
     &, rhokm_u,rhokm_v                                                 &

! in diagnostics (started or from) BDY_LAYR
     &, T_incr_diag_bl, q_incr_diag_bl, qcl_incr_diag_bl                &
     &, qcf_incr_diag_bl, u_incr_diag_bl, v_incr_diag_bl                &
     &, cfl_incr_diag_bl, cff_incr_diag_bl, cf_incr_diag_bl             &
     &, cd, ch, e_sea, fqT, ftl, h_sea, rib_gb                          &
     &, taux, tauy, vshr, zht                                           &

     &, bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6     &
     &, bl_type_7                                                       &
     &, z0m_gb, z0m_eff_gb, z0h_eff_gb                                  &
     &, fme, rhokh                                                      &

! in data required to calculate increment diagnostics
     &, theta_star,q_star                                               &
! in logical for scm surface forcing
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

! IN additional variables for MOSES II
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
     &, slm                                                             &
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

! IN Tracer fluxes - kdcorbin, 05/10
     &, TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4          &
     &, TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8          &
     &, TRACER_FLUX9, TRACER_FLUX10,TRACER_FLUX11,TRACER_FLUX12         &
     &, TRACER_FLUX13,TRACER_FLUX14,TRACER_FLUX15,TRACER_FLUX16         &
     &, TRACER_FLUX17,TRACER_FLUX18,TRACER_FLUX19,TRACER_FLUX20         &

! IN CO2 Global Emissions - kdcorbin, 05/10
     &, co2emitmass                                                     &

!IN additional variables for soil moisture nudging scheme
     &, WT_EXT,RA                                                       &

!IN additional variables for anthropogenic heat source 
     &, ANTHROP_HEAT                                                    &

! in/out
     &, t_soil, ti, t_surf, ti_gb                                       &
     &, area_cloud_fraction, bulk_cloud_fraction                        &
     &, T_latest, q_latest, qcl_latest, qcf_latest                      &
     &, cf_latest, cfl_latest, cff_latest                               &
     &, R_u, R_v, R_w, cloud_fraction_liquid, cloud_fraction_frozen     &
     &, zh,sum_eng_fluxes,sum_moist_flux,rhcpt                          &

! In/Out tracer fields
     &, aerosol, free_tracers                                           &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, DRYDEP2, so2, dms, so4_aitken, so4_accu, so4_diss, nh3          &
     &, soot_new, soot_aged, soot_cld, bmass_new, bmass_agd             &
     &, bmass_cld, ocff_new, ocff_aged, ocff_cld, co2, ozone_tracer     &

! in/out fields
     &, ecan, ei, ext, snowmelt                                         &
     &, t1_sd, q1_sd, ntml, cumulus, l_pc2_diag_sh_pts                  &
     &, nbdsc, ntdsc                                                    &
     &, surf_ht_flux_land, snomlt_surf_htf, cH_term                     &

! INOUT additional variables for MOSES II
     &, TSTAR_TILE,FQT_TILE,EPOT_TILE,FTL_TILE                          &
     &, SNOW_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,OLR                   &
     &, TSTAR_LAND,TSTAR_SICE,TSTAR_SSI                                 &
     &, rib_ssi,taux_land,taux_ssi,tauy_land,tauy_ssi                   &

! OUT additional variables for MOSES II
     &, ESOIL_TILE,ES,EI_TILE                                           &
     &, Q1P5M_TILE,T1P5M_TILE,ECAN_TILE,MELT_TILE                       &

! error information
     &, Error_code,BL_diag)

! purpose: Interface to boundary layer implicit solver and cloud scheme
!
!          IN/OUT etc intents to be added later.
!
! method:
!
! code description:
!   language: fortran 77 + cray extensions
!   this code is written to umdp3 programming standards.

      Use cv_cntl_mod, Only:                                            &
          lcv_3d_cca, lcv_pc2_diag_sh
      Use cv_run_mod, Only:                                             &
          tice
      Use bl_diags_mod, Only:                                           &
          strnewbldiag

      Implicit None

!C_DUST_NDIV.............................................................
! Description: Contains parameters for mineral dust code
! Current Code Owner: Stephanie Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  5.5      12/02/03  Original Code.   Stephanie Woodward
!
! Declarations:
!
      INTEGER NDIV        ! number of particle size divisions
      PARAMETER (NDIV = 6)
!.....................................................................
!C_DUSTGEN..............................................................
! Description: Contains parameters for mineral dust generation
! Current Code Owner: Stephanie Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  5.5      12/02/03  Original Code.   Stephanie Woodward
!  6.2      14/02/06  Alternative emissions terms for HadGEM1A.
!                                                     Stephanie Woodward
!  6.4      03/01/07  Alternative U*t terms.   Stephanie Woodward
! Parameters for mineral dust generation
!
      REAL, PARAMETER :: HORIZ_C = 2.61 ! C in horizontal flux calc.     
      REAL, PARAMETER :: VERT_A = 13.4 ! A in vertical flux calc
      REAL, PARAMETER :: VERT_B = -6. ! B in vertical flux calc

      REAL,PARAMETER :: RHOP = 2.65E+3  ! density of a dust particle
      REAL,PARAMETER :: Z0B = 0.0003    !roughness length for bare soil
      REAL Z0S(NDIV)  ! smooth roughness len (calc.d from part size)
      REAL DREP(NDIV) ! representative particle diameter
      REAL DMAX(NDIV) ! max diameter of particles in each div.
      REAL DMIN(NDIV) ! min diameter of particles in each div.
                         ! note that by using two arrays here we can set
                         ! up overlapping divisions, however this means
                         ! we have to be careful to make them consistent
!
      DATA Z0S/ .374894E-08, .118552E-07, .374894E-07, .118552E-06,     &
     &           .374894E-06, .118552E-05/
      DATA DREP/ .112468E-06, .355656E-06, .112468E-05, .355656E-05,    &
     &           .112468E-04, .355656E-04/
      DATA DMAX/2.0E-7,6.32456E-7,2.0E-6,                               &
     &          6.32456E-6,2.0E-5,6.32456E-5/
      DATA DMIN/6.32456E-8,2.0E-7,6.32456E-7,                           &
     &          2.0E-6,6.32456E-6,2.0E-5/
!.......................................................................


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

! Substep number for ATMPHYS2
      Integer                                                           &
     &  Substep_Number                                                  &
     &, Num_Substeps

! Switch for calculating exchange coeffs from latest values.
      Logical                                                           &
     &  L_phys2_substep

      Logical :: L_us_blsol   ! switch for BL numerical solver

      Real :: Puns, Pstb ! parameters for uncond stable numerical solver
                         ! Puns : used in an unstable BL column
                         ! Pstb : used in an stable BL column

      LOGICAL                                                           &
     & lq_mix_bl              ! TRUE if mixing ratios used in
!                             ! boundary layer code

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
!                          south, east or west of the processor grid
!                          (array index PNorth etc. from parparm.h)
! Switch for  surface flux forcing in SCM
     &, L_flux_bc

! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, rhc_row_length                                                  &
     &, rhc_rows                                                        &
     &, n_rows                                                          &
     &, land_points                                                     &
                    ! IN No.of land points being processed, can be 0.
     &, ntiles                                                          &
                    ! IN No. of land-surface tiles ( MOSES II )
     &, nice                                                            &
                    ! IN No. of sea ice catagories
     &, model_levels                                                    &
     &, wet_model_levels                                                &
     &, bl_levels                                                       &
     &, dst_levels                                                      &
                    ! number of deep soil temperature levels
     &, dsm_levels                                                      &
                    ! number of deep soil moisture levels
     &, cloud_levels                                                    &
     &, n_cca_levels                                                    &
                      ! Number of levels for conv cloud amount :
!                       1 for 2D, nlevs for 3D.
     &, tr_levels     ! number of free tracer levels
       Integer                                                          &
     &  NumCycles                                                       &
                   ! Number of cycles
     &, CycleNo                                                         &
                   ! Sweep number

     &, DIM_CS1, DIM_CS2
!
! Model switches
      Integer                                                           &
     &  model_domain

      Logical                                                           &
     &  L_CAL360                                                        &
                    ! true if using 360 day calender
     &, L_emcorr                                                        &
                    ! true if energy correction scheme is to be used.
     &, L_dry                                                           &
                    ! true if model to be run with no moisture
     &, L_eacf      ! true if using empirically adjusted cloud fraction

      Logical                                                           &
     &  L_area_cloud                                                    &
                           ! Switch for area cloud fraction param
     &, L_ACF_Cusack                                                    &
                           ! ... to select Cusack (if PC2 off)
     &, L_ACF_Brooks                                                    &
                           ! ... to select Brooks
     &, L_RHCPT                                                         &
                 ! Switch for 3D diagnosed RHcrit not 1D parameter
     &, Ltimer                                                          &
                 ! true then output some timing information
     &, L_Murk                                                          &
                          ! Switch for (visibility) aerosol
     &, L_murk_advect                                                   &
                          ! Switch for advecting aerosol
     &, L_bl_tracer_mix                                                 &
                          ! Switch for BL mixing of free tracers
     &, L_DUST                                                          &
                          ! Switch for mineral dust
     &, L_CAM_DUST                                                      &
                          ! Old dust uplift scheme used in CAM NWP models 
     &, L_sulpc_so2                                                     &
                          ! Switch for Sulphur Cycle
     &, L_sulpc_nh3                                                     &
                          ! NH3 included in Sulphur Cycle
     &, L_sulpc_dms                                                     &
                          ! DMS included in Sulphur Cycle
     &, L_soot                                                          &
                          ! Switch for Soot Cycle
     &, L_biomass                                                       &
                          ! Switch for biomass aerosol
     &, L_ocff                                                          &
                          ! Switch for fossil-fuel OC aerosol
     &, L_co2_interactive                                               &
                          ! Switch for interactive CO2
     &, L_pc2                                                           &
                          ! Switch for PC2 cloud scheme
     &, L_psd                                                           &
                          ! Use generic ice particle size distribution
     &, L_co2_emits                                                     &
                          ! Switch to include surface emissions
     &, L_CO2_TRACER                                                    &
            ! rml, 1/7/13   Switch to put co2 fluxes into passive tracer
     &, L_use_bl_diag_term                                              &
                          ! Switch to calc ch_term in ni_imp_ctl()
     &, L_ukca                                                          &
                          ! Switch for UKCA sub-model
     &, L_sice_heatflux                                                 &
                          ! Switch for semi-impl sea-ice temperature
     &, L_subfilter_vert                                                &
     &, L_subfilter_horiz                                               &
                          ! Switches for subgrid turbulence scheme
     &, l_use_cariolle   
                          ! Switch for cariolle ozone tracer scheme
      Integer                                                           &
     &  trweights1         ! IN switch to use implicit weights
!                          !    of 1 for tracers (if =1)
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
     &, rhcrit(wet_model_levels)                                        &
                                  ! IN Critical relative humidity.
                                  ! the values need to be tuned
                                  ! for the given set of levels.
     &, alpha_tr(bl_levels)  ! Implicit weights for tracers
!                            !  = alpha_cd from RUN_BL namelist
!                            ! or = 1 if TRWEIGHTS1 = ON

      Integer                                                           &
     &  tr_vars              ! IN number of free tracer variables

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
! Start blopt8a

! Description:
!   Permissible settings for BL options.
!
! Current Code Owner: J. M. Edwards
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 6.2      27/01/06 Original code.  J. M. Edwards
!
      INTEGER, PARAMETER :: Off = 0  ! Switch disabled
      INTEGER, PARAMETER :: On  = 1  ! Switch enabled
!
      INTEGER, PARAMETER :: BrownGrant97 = 1
      INTEGER, PARAMETER :: BrownGrant97_limited = 2
!     Options for non-gradient stress following
!     Brown and Grant (1997), version 2 including a limit on its size
!
!     Options for flux gradient formulation
      INTEGER, PARAMETER :: Locketal2000   = 0
!       Flux gradients as in Lock et al. (2000)
      INTEGER, PARAMETER :: HoltBov1993 = 1
!       Flux gradients as in Lock et al (2000) but using
!       coefficients from Holtslag and Boville (1993)
      INTEGER, PARAMETER :: LockWhelan2006 = 2
!       Flux gradients as in Lock and Whelan (2006)
!
!     Options for form drag
      INTEGER, PARAMETER :: No_drag         = 0
      INTEGER, PARAMETER :: Effective_z0    = 1
      INTEGER, PARAMETER :: Explicit_stress = 2
!
!     Options for marine boundary layers
      INTEGER, PARAMETER :: Fixed_Z0T = 0
!       Stanard flixed value of thermal roughness length over sea
      INTEGER, PARAMETER :: SurfDivZ0T = 1
!       Thermal roughness length over sea defined from surface
!       divergence theory
      INTEGER, PARAMETER :: DynDiag_ZL = 1
!       The ratio of the height of the inversion to the surface
!       Obukhov length is used a dynamic criterion in the
!       diagnosis of BL types
!
      INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
!
!     Options for stable boundary layers
      INTEGER, PARAMETER ::  Long_tails           = 0
      INTEGER, PARAMETER ::  Sharpest             = 1
      INTEGER, PARAMETER ::  Sharp_sea_long_land  = 2
      INTEGER, PARAMETER ::  Mes_tails            = 3
      INTEGER, PARAMETER ::  Louis_tails          = 4
      INTEGER, PARAMETER ::  Depth_based          = 5
      INTEGER, PARAMETER ::  Sharp_sea_mes_land   = 6
      INTEGER, PARAMETER ::  LEM_stability        = 7
      INTEGER, PARAMETER ::  Sharp_sea_Louis_land = 8

!     Options for Prandtl number (in local Ri scheme)
      INTEGER, PARAMETER ::  Constant_SBL = 0
      INTEGER, PARAMETER ::  LockMailhot2004 = 1

! End blopt8a

! Local variables
      Integer, Parameter :: sect = 3               ! BL section
      Integer            :: item                   ! stash item code
      Integer            :: im_index               ! internal model

      Character (Len=*), Parameter :: RoutineName='bl_w_mixing'
      Character (Len=80)           :: Cmessage

! Diagnostics info
       Real                                                             &
     & STASHwork3(*)                                                    &
                     ! STASH workspace for section 3 (Boundary Layer)
     &,STASHwork9(*) ! STASH workspace for section 9 (LS Cloud)


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
     &, FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)

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

!      Real
!     &  p_theta_levels(1-off_x:row_length+off_x,
!     &                 1-off_y:rows+off_y, model_levels)
!     &, exner_rho_levels(1-off_x:row_length+off_x,
!     &                   1-off_y:rows+off_y, model_levels)

! Data arrays
      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      0:model_levels)                                             &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, p_layer_centres(row_length, rows, 0:model_levels)               &
     &, p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, p_star(row_length, rows)                                        &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      wet_model_levels)                                           &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_model_levels)                                         &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_model_levels)                                         &
     &, T_conv(row_length, rows, model_levels)                          &
     &, q_conv(row_length, rows, wet_model_levels)                      &
     &, qcl_conv(row_length, rows, wet_model_levels)                    &
     &, qcf_conv(row_length, rows, wet_model_levels)                    &
     &, u_conv(1-off_x:row_length+off_x,1-off_y:rows+off_y,             &
     &              bl_levels)                                          &
     &, v_conv(1-off_x:row_length+off_x,1-off_y:n_rows+off_y,           &
     &              bl_levels)                                          &
     &, ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)                                     &
     &, cca_2d(row_length, rows)

      logical                                                           &
     &  land_sea_mask(row_length, rows)                                 &
     &, L_scrn                                                          &
                                 ! Logical to control output
                                 !    of screen level T,Q,QCL,QCF
     &, L_plsp                   ! Logical to control output
                                 !    of Probability of LS Precip

! ancillary arrays and fields required to be saved from timestep to
! timestep.
      Integer                                                           &
     &  land_index (land_points)      ! set from land_sea_mask

      Real                                                              &
     &  u_0(row_length, rows)                                           &
                                ! set to zero
     &, v_0(row_length, n_rows)                                         &
                                ! set to zero
     &, sil_orog_land (land_points)                                     &
                                    ! orog/qrparm.orog.as
     &, ho2r2_orog (land_points)                                        &
                                 ! orog/qrparm.orog.h2root2
     &, smvccl (land_points)                                            &
                             ! soil/qrparm.soil.crit
     &, smvcwt (land_points)                                            &
                             ! soil/qrparm.soil.wilt
     &, smvcst (land_points)                                            &
                             ! soil/qrparm.soil.satn
     &, sthf(land_points,dsm_levels)                                    &
                                ! IN Frozen soil moisture content of
                                !     each layer as a fraction of
                                !     saturation.
     &, sthu(land_points,dsm_levels)                                    &
                                ! IN Unfrozen soil moisture content
                                !    of each layer as a fraction of
                                !    saturation.
     &, ice_fract (row_length, rows)                                    &
                                     ! ice/qrclim.ice.(month)
     &, di(row_length, rows)                                            &
                              ! ice/qrclim.ice_thick.(month)
     &, ice_fract_ncat (row_length, rows, nice)                         &
                                                !ice fract on catagories
     &, di_ncat(row_length, rows, nice)  ! ice thickness on catagories

      Real                                                              &
     & SURF_RADFLUX(row_length,rows)! IN Surface net radiation
!                                    (W/sq m,  positive downwards).

      Real ::                                  &
       cca (row_length, rows, n_cca_levels)    &
     , ccw(row_length, rows, wet_model_levels)   ! Conv water kg/kg

      Integer                                                           &
     &  ccb (row_length, rows)                                          &
     &, cct (row_length, rows)

! in data required to calculate increment diagnostics
      Real                                                              &
     & theta_star(1-off_x:row_length+off_x,1-off_y:rows+off_y,          &
     &                                                model_levels)     &
     &,    q_star(1-off_x:row_length+off_x,1-off_y:rows+off_y,          &
     &                                            wet_model_levels)

! in variables passed from BDY_LAYR to IMP_SOLVER
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
     &, z1_tq(row_length,rows)
!ajm variable added
      Real                                                              &
     &  rhokm_u(row_length, rows, bl_levels)                            &
     &, rhokm_v(row_length, n_rows, bl_levels)

! in diagnostics (started or from) BDY_LAYR
      Real                                                              &
            ! output from bdy_layr.
     &  u_incr_diag_bl(row_length,rows,model_levels)                    &
                              ! u wind      increment for STASH
     &, v_incr_diag_bl(row_length,n_rows,model_levels)                  &
                              ! v wind      increment for STASH
     &, T_incr_diag_bl(row_length,rows,model_levels)                    &
     &, q_incr_diag_bl(row_length,rows,wet_model_levels)                &
     &, qcl_incr_diag_bl(row_length,rows,wet_model_levels)              &
     &, qcf_incr_diag_bl(row_length,rows,wet_model_levels)              &
     &, cf_incr_diag_bl(row_length,rows,wet_model_levels)               &
     &, cfl_incr_diag_bl(row_length,rows,wet_model_levels)              &
     &, cff_incr_diag_bl(row_length,rows,wet_model_levels)              &
     &, cd(row_length, rows)                                            &
                              ! needed as diagnostic ?
     &, ch(row_length, rows)                                            &
                              ! needed as diagnostic ?
     &, e_sea(row_length, rows)                                         &
                                 ! needed as diagnostic ?
     &, fqT(row_length, rows, bl_levels)                                &
                                          ! needed as diagnostic ?
     &, ftl(row_length, rows, bl_levels)                                &
                                          ! needed as diagnostic
     &, h_sea(row_length, rows)                                         &
                                 ! needed as diagnostic ?
     &, rib_gb(row_length,rows)                                         &
                                     ! Mean bulk Richardson number for
!                                     lowest layer.
     &, taux(row_length, rows, bl_levels)                               &
                                           ! needed as diagnostic
     &, tauy(row_length, n_rows, bl_levels)                             &
                                             ! needed as diagnostic
     &, vshr(row_length, rows)                                          &
     &, zht(row_length, rows)                                           &
                                      ! Max height of turb mixing
     &, bl_type_1(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if stable
!                                 !     b.l. diagnosed, 0.0 otherwise.
     &, bl_type_2(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if Sc over
!                                 !     stable surface layer diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_3(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if well
!                                 !     mixed b.l. diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_4(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer (not over
!                                 !     cumulus) diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_5(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer over cumulus
!                                 !     diagnosed, 0.0 otherwise.
     &, bl_type_6(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if a
!                                 !     cumulus capped b.l. diagnosed,
!                                 !     0.0 otherwise.
     &, bl_type_7(row_length,rows)                                      &
                                  ! IN Indicator set to 1.0 if a
!                                 !     shear-dominated b.l. diagnosed,
!                                 !     0.0 otherwise.
     &, z0m_eff_gb(row_length,rows)                                     &
     &, z0h_eff_gb(row_length,rows)                                     &
                                  ! IN Effective grid-box roughness
                                  !     lengths for diagnostics
     &, fme(row_length, rows)                                           &
     &, rhokh (row_length, rows, bl_levels)

      REAL :: stress,stress_1

      Real                                                              &
     &  rho_aresist(row_length,rows)                                    &
                                ! OUT RHOSTAR*CD_STD*VSHR for SULPHUR
                                !     cycle
     &, aresist(row_length,rows)                                        &
                                ! OUT 1/(CD_STD*VSHR) for Sulphur cycle
     &, resist_b(row_length,rows)                                       &
                                ! OUT (1/CH-1/(CD_STD)/VSHR for Sulphur
                                !     cycle
     &, R_B_DUST(ROW_LENGTH,ROWS,NDIV)                                  &
                                        ! IN surf layer res for dust
     &, WE_LIM(row_length,rows,3)                                       &
                                    ! IN rho*entrainment rate implied by
!                                   !     placing of subsidence
     &, ZRZI(row_length,rows,3)                                         &
                                    ! IN (z-z_base)/(z_i-z_base)
     &, T_FRAC(row_length,rows,3)                                       &
                                    ! IN a fraction of the timestep
     &, WE_LIM_DSC(row_length,rows,3)                                   &
!                                   ! IN rho*entrainment rate implied by
!                                   !     placing of subsidence
     &, ZRZI_DSC(row_length,rows,3)                                     &
                                    ! IN (z-z_base)/(z_i-z_base)
     &, T_FRAC_DSC(row_length,rows,3)                                   &
!                                   ! IN a fraction of the timestep
     &, Z_HALF(row_length,rows,BL_LEVELS)                               &
                                    ! IN Z_HALF(*,K) is height of half
!                                   !    level k-1/2.
     &, ZHSC(row_length,rows)                                           &
                                    ! IN Top of decoupled layer
     &, RA(land_points)                                                 &
                                    ! IN aerodynamic resiatnce (s/m)
     &, WT_EXT(land_points,dsm_levels) ! IN cumulative fract of trans'n

!     Declaration of new BL diagnostics.
      Type (Strnewbldiag) :: BL_diag


      Integer                                                           &
     &  KENT(row_length,rows)                                           &
                                    ! IN grid-level of SML inversion
     &, KENT_DSC(row_length,rows)   ! IN grid-level of DSC inversion


! Mineral dust source flux for tracer mixing
      REAL, INTENT(IN) ::                                               &
     &  DUST_FLUX ( ROW_LENGTH, ROWS, NDIV )                            &
     &, U_S_T_TILE(LAND_POINTS,NTILES,NDIV)                             &
                                           !OUT threshold frict. vel
     &, U_S_T_DRY_TILE(LAND_POINTS,NTILES,NDIV)                         &
                                               !OUT dry soil value
     &, U_S_STD_TILE(LAND_POINTS,NTILES)!OUT friction velocity

! Emissions for tracer mixing
      Real, Intent(In) ::                                               &
     &  so2_hilem ( row_length, rows )                                  &
     &, so2_em    ( row_length, rows )                                  &
     &, nh3_em    ( row_length, rows )                                  &
     &, dms_em    ( row_length, rows )                                  &
     &, soot_hilem( row_length, rows )                                  &
     &, soot_em   ( row_length, rows )                                  &
     &, ocff_hilem( row_length, rows )                                  &
     &, ocff_em   ( row_length, rows )

      Real                                                              &
             ! (output from sf_exch)
     &  rho_cd_modv1(row_length, rows)                                  &
     &, rho_km_var(row_length,rows,2:bl_levels)

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

! arguments with intent in/out. ie: input variables changed on output.
      Real                                                              &
     &  T_surf(row_length, rows)                                        &
     &, t_soil(land_points,dsm_levels)                                  &
                                       ! slt/qrclim.slt_pm(lev).(month)
     &, ti(row_length, rows, nice)                                      &
                                       ! ice temperatures on catagories
     &, R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels)                             &
     &, T_latest(row_length, rows, model_levels)                        &
     &, q_latest(row_length, rows, wet_model_levels)                    &
     &, qcl_latest(row_length, rows, wet_model_levels)                  &
     &, qcf_latest(row_length, rows, wet_model_levels)                  &
     &, cf_latest(row_length, rows, wet_model_levels)                   &
     &, cfl_latest(row_length, rows, wet_model_levels)                  &
     &, cff_latest(row_length, rows, wet_model_levels)                  &
     &, area_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, bulk_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, cloud_fraction_liquid(row_length, rows, wet_model_levels)       &
     &, cloud_fraction_frozen(row_length, rows, wet_model_levels)       &
     &, zh(row_length, rows)                                            &
     &, sum_eng_fluxes(row_length, rows)                                &
     &, sum_moist_flux(row_length, rows)

! Tracer variables
      Real, Intent(InOut) ::                                            &
     &  aerosol     ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,free_tracers( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y,                         &
     &                tr_levels, tr_vars)

      Real, Intent(InOut) ::                                            &
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
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DRYDEP2(ROW_LENGTH,ROWS,NDIV) !dry dep though grav. set.

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
     & ,bmass_agd   ( 1 - off_x : row_length + off_x,                   &
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
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ozone_tracer( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )         

      Real                                                              &
             !
     &  ecan(row_length, rows)                                          &
                                        !output from sf_evap.
     &, ei(row_length, rows)                                            &
                                        !output from sf_evap.
     &, ext(land_points,dsm_levels)                                     &
                                    ! Extraction of water from each
!                                    soil layer (kg/m2/s).
     &, snowmelt(row_length, rows)                                      &
                                        !output from sf_evap.
     &, t1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, q1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, surf_ht_flux_gb(row_length, rows)                               &
                                           !
     &, snomlt_surf_htf(row_length, rows)

      Logical ::                     &
       cumulus (row_length, rows)    & ! *APL bl convection flag
     , l_pc2_diag_sh_pts(row_length, rows) ! Carry diagnostic shallow
                                           ! convective information for PC2

      Integer                                                           &
     &  ntml (row_length, rows)                                         &
     &, nbdsc(row_length, rows)                                         &
     &, ntdsc(row_length, rows)

      Real                                                              &
     & cH_term(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &          model_levels-1)

! IN additional variables for MOSES II
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



      LOGICAL                                                           &
     & L_NEG_TSTAR                 ! IN Switch for -ve TSTAR error check

      INTEGER                                                           &
     & TILE_PTS(NTYPE)                                                  &
                                 ! IN Number of tile points.
     &,TILE_INDEX(LAND_POINTS,NTYPE)
!                                ! IN Index of tile points.

      REAL                                                              &
     & TILE_FRAC(land_points,ntiles)                                    &
                                ! IN fractional coverage for each
                                !    surface tile
     &,CANOPY(land_points,ntiles)                                       &
                                   ! IN Surface/canopy water (kg/m2)
     &,ALPHA1(land_points,ntiles)                                       &
                                  ! IN Mean gradient of saturated
!                                 specific humidity with
!                                 respect to temperature between
!                                 the bottom model layer and the
!                                 tile surfaces.
     &,FRACA(land_points,ntiles)                                        &
                                   ! IN Fraction of surface
                                !            moisture flux with only
                                !            aerodynamic resistance.
     &,RHOKH_TILE(land_points,ntiles)                                   &
                                      ! IN
!                                 Tile surface exchange coefficients
!                                 for heat
     &,SMC(LAND_POINTS,NTILES)                                          &
                                ! IN Soil moisture content in root depth
!                                  (kg/m2).
     &,CHR1P5M(LAND_POINTS,NTILES)                                      &
                                   ! IN Ratio of coefficients reqd for
!                                 calculation of 1.5 m T.
     &,RESFS(land_points,ntiles)                                        &
                                 ! IN Combined soil, stomatal
!                                 and aerodynamicresistance
!                                 factor = PSIS/(1+RS/RA) for
!                                 fraction (1-FRACA)
     &,Z0HSSI(row_length,rows)                                          &
     &,Z0MSSI(row_length,rows)                                          &
                                ! IN Roughness lengths over sea (m).
     &,Z0M_GB(row_length,rows)                                          &
                                ! IN Gridbox mean roughness length 
!                               !    for momentum (m).
     &,CANHC_TILE(LAND_POINTS,NTILES)                                   &
!                               ! IN Areal heat capacity of canopy
!                               !    for land tiles (J/K/m2).
     &,FLAKE(LAND_POINTS,NTILES)                                        &
                                   ! IN Lake fraction.
     &,WT_EXT_TILE(LAND_POINTS,DSM_LEVELS,NTILES)                       &
!                               ! IN Fraction of evapotranspiration
!                               !    which is extracted from each
!                               !    soil layer by each tile.
     &,LW_DOWN(ROW_LENGTH,ROWS)                                         &
                                ! IN Surface downward LW radiation
!                               !    (W/m2).
     &,lai_ft(land_points,npft)                                         &
                                   ! IN LAI on vegetated tiles
     &,canht_ft(land_points,npft)                                       &
                                   ! IN CANHT on vegetated tiles
     &,SW_TILE(LAND_POINTS,NTILES)                                      &
                                   ! IN Surface net SW radiation on land
!                               !    tiles (W/m2).
     &,ASHTF_TILE(LAND_POINTS,NTILES)                                   &
!                               ! IN Coefficient to calculate
!                               !    surface heat flux into land
!                               !    tiles.
     &,FQT_ICE(ROW_LENGTH,ROWS)                                         &
                                ! IN Surface FQT for sea-ice
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
                                ! IN Surface FTL for sea-ice
     &,RESFT(LAND_POINTS,NTILES)                                        &
                                   ! IN Total resistance factor.
!                               !    FRACA+(1-FRACA)*RESFS for
!                               !    snow-free land, 1 for snow.
     &,RHOKH_SICE(ROW_LENGTH,ROWS)                                      &
!                               ! IN Surface exchange coefficients
!                               !    for sea and sea-ice
     &,RHOKPM(LAND_POINTS,NTILES)                                       &
                                ! IN Land surface exchange coeff.
     &,RHOKPM_POT(LAND_POINTS,NTILES)                                   &
!                               ! IN Land surface exchange coeff.
!                                    for potential evaporation.
     &,EPOT_TILE(land_points,ntiles)                                    &
!                               ! INOUT surface tile potential
!                               !       evaporation
     &,RHOKPM_SICE(ROW_LENGTH,ROWS)                                     &
!                               ! IN Sea-ice surface exchange coeff.
     &,Z0H_TILE(LAND_POINTS,NTILES)                                     &
                                ! IN Tile roughness lengths for heat
!                               !    and moisture (m).
     &,Z0M_TILE(LAND_POINTS,NTILES)                                     &
                                ! IN Tile roughness lengths for
!                               !    momentum.
     &,CHR1P5M_SICE(ROW_LENGTH,ROWS)                                    &
!                               ! IN CHR1P5M for sea and sea-ice
!                               !    (leads ignored).
     &,FLAND(LAND_POINTS)                                               &
                                ! IN Land fraction on land tiles.
     &,FLANDG(ROW_LENGTH,ROWS)                                          &
                                ! IN Land fraction on all tiles.
     &,FLANDG_U(ROW_LENGTH,ROWS)                                        &
!                               ! IN Land frac (on U-grid, with 1st
!                               !    and last rows undefined or, at
!                               !    present, set to "missing data")
     &,FLANDG_V(ROW_LENGTH,N_ROWS)                                      &
!                               ! IN Land frac (on V-grid, with 1st
!                               !    and last rows undefined or, at
!                               !    present, set to "missing data")
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
!                               ! IN Open sea sfc temperature (K).
     &,VSHR_LAND(ROW_LENGTH,ROWS)                                       &
!                               ! IN VSHR over land part of gridbox.
     &,VSHR_SSI(ROW_LENGTH,ROWS)                                        &
                                ! IN VSHR over sea part of gridbox.
     &, gc(land_points,ntiles)                                          &
                                  !Stomatal conductance to evapn
!                                 !    for land tiles (m/s).
     &, aresist_tile(land_points,ntiles)                                &
!                                  !1/(CD_STD*VSHR) on land tiles
     &, resist_b_tile(land_points,ntiles)                               &
!                                  !(1/CH-1/CD_STD)/VSHR on land tiles
     &,ANTHROP_HEAT(NTILES)  
                                ! IN Additional heat source on tiles
                                !    for anthropogenic urban heat
                                !    source (W/m2) 

! IN MOSES II additional STASH variables
      REAL                                                              &
     & GS(LAND_POINTS)                                                  &
                                ! IN "Stomatal" conductance to
!                               !    evaporation (m/s).
     &,GPP(LAND_POINTS)                                                 &
                                ! IN Gross primary productivity
!                               !    (kg C/m2/s).
     &,NPP(LAND_POINTS)                                                 &
                                ! IN Net primary productivity
!                               !    (kg C/m2/s).
     &,RESP_P(LAND_POINTS)                                              &
                                ! IN Plant respiration (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,GPP_FT(LAND_POINTS,NTILES)                                       &
                                ! IN Gross primary productivity
                                !    on PFTs (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,NPP_FT(LAND_POINTS,NTILES)                                       &
                                ! IN Net primary productivity
                                !    on PFTs (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,RESP_P_FT(LAND_POINTS,NTILES)                                    &
                                  !IN Plant respiration on PFTs
                                  !     (kg C/m2/s).
     &,RESP_S(LAND_POINTS,DIM_CS1)                                      &
                                   ! IN Soil respiration (kg C/m2/s).
     &,RESP_S_TOT(DIM_CS2)                                              & 
                                   ! IN Total soil resp'n (kg C/m2/s).
     &,RESP_S_TILE(LAND_POINTS,NTILES)                                  &
                                   ! IN Soil respiration on tiles (kg C/m2/s).
                                   !   kdcorbin, 10/10
     &,CS(LAND_POINTS,DIM_CS1)                                          &
                                   ! IN Soil carbon
     &,RIB_TILE(LAND_POINTS,NTILES)                                     &
!                               ! IN RIB for land tiles.
     &,FSMC(LAND_POINTS,NPFT)                                           &
                                ! IN Moisture availability factor.
     &,CATCH(LAND_POINTS,NTILES)                                        &
                                ! IN Surface/canopy water capacity
                                !    of snow-free land tiles (kg/m2).
     !kdcorbin,11/10 - changed from NPFT
     &,G_LEAF(LAND_POINTS,NTILES)                                       &
                                ! IN Leaf turnover rate (/360days).
     &,CO2_EMITS(ROW_LENGTH,ROWS)                                       &
                                 !IN CO2 Emissions
     &,CO2FLUX(ROW_LENGTH,ROWS)                                         &
                                      ! IN ocean CO2 flux
     &,co2_flux_tot(row_length,rows)                                    &
                                      !  total CO2 flux
     &,land_co2(land_points)          !  terrestrial CO2 flux

! INOUT tracer fluxes - kdcorbin, 05/10
      REAL                                                              &
     & TRACER_FLUX1(ROW_LENGTH,ROWS),      TRACER_FLUX2(ROW_LENGTH,ROWS)   &
     &, TRACER_FLUX3(ROW_LENGTH,ROWS),     TRACER_FLUX4(ROW_LENGTH,ROWS)   &
     &, TRACER_FLUX5(ROW_LENGTH,ROWS),     TRACER_FLUX6(ROW_LENGTH,ROWS)   &
     &, TRACER_FLUX7(ROW_LENGTH,ROWS),     TRACER_FLUX8(ROW_LENGTH,ROWS)   &
     &, TRACER_FLUX9(ROW_LENGTH,ROWS),     TRACER_FLUX10(ROW_LENGTH,ROWS)  &
     &, TRACER_FLUX11(ROW_LENGTH,ROWS),    TRACER_FLUX12(ROW_LENGTH,ROWS)  &
     &, TRACER_FLUX13(ROW_LENGTH,ROWS),    TRACER_FLUX14(ROW_LENGTH,ROWS)  &
     &, TRACER_FLUX15(ROW_LENGTH,ROWS),    TRACER_FLUX16(ROW_LENGTH,ROWS)  &
     &, TRACER_FLUX17(ROW_LENGTH,ROWS),    TRACER_FLUX18(ROW_LENGTH,ROWS)  &
     &, TRACER_FLUX19(ROW_LENGTH,ROWS),    TRACER_FLUX20(ROW_LENGTH,ROWS)

! INOUT CO2 Global Emissions - kdcorbin, 05/10
      REAL co2emitmass


! INOUT additional variables for MOSES II
      REAL                                                              &
     & TSTAR_TILE(land_points,ntiles)                                   &
                                ! INOUT Surface tile temperature
     &,FQT_TILE(land_points,ntiles)                                     &
                                ! INOUT surface tile moisture flux
     &,FTL_TILE(land_points,ntiles)                                     &
!                               ! INOUT surface tile heat flux
     &,SNOW_TILE(LAND_POINTS,NTILES)                                    &
!                               ! INOUT Snow on tiles (kg/m2).
     &,LE_TILE(LAND_POINTS,NTILES)                                      &
                                   ! INOUT Surface latent heat flux for
!                               !       land tiles (W/m2).
     &,RADNET_SICE(ROW_LENGTH,ROWS)                                     &
!                               ! INOUT Sea-ice surface net radiation.
     &,RADNET_TILE(LAND_POINTS,NTILES)                                  &
!                               ! INOUT Tile surface net radiation.
     &,OLR(ROW_LENGTH,ROWS)                                             &
                                ! IN    TOA - surface upward LW on
!                               !       last radiation timestep
!                               ! OUT   Corrected TOA outward LW
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
                                   ! INOUT Land mean sfc temperature (K)
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! INOUT Sea-ice sfc temperature (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS)                                       &
                                   ! INOUT Sea mean sfc temperature (K).
     &,RIB_SSI(ROW_LENGTH,ROWS)                                         &
                                   ! INOUT Sea mean bulk Richardson
!                                          number for lowest layer.
     &,TAUX_LAND(ROW_LENGTH,ROWS)                                       & 
                                   ! INOUT W'ly compt of land sfc wind
!                                  !       stress (N/sq m). (On U-grid
!                                  !       with first and last rows
!                                  !       undefined or, at present,
!                                  !       set to missing data
     &,TAUX_SSI(ROW_LENGTH,ROWS)                                        & 
                                   ! INOUT W'ly compt of sea sfc wind
!                                  !       stress (N/sq m). (On U-grid
!                                  !       with first and last rows
!                                  !       undefined or, at present,
!                                  !       set to missing data
     &,TAUY_LAND(ROW_LENGTH,N_ROWS)                                     & 
                                   ! INOUT S'ly compt of land sfc wind
!                                  !       stress (N/sq m).  On V-grid;
!                                  !       comments as per TAUX.
     &,TAUY_SSI(ROW_LENGTH,N_ROWS) ! INOUT S'ly compt of sea sfc wind
!                                  !       stress (N/sq m).  On V-grid;
!                                  !       comments as per TAUX.

! OUT additional variables for MOSES II
      REAL                                                              &
     & ESOIL_TILE(land_points,ntiles)                                   &
                                ! OUT Evaporation from bare soil (kg/m2)
     &,ES(row_length,rows)                                              &
                                ! OUT Surface evapotranspiration from
!                               !     soil moisture store (kg/m2/s).
     &,EI_TILE(LAND_POINTS,NTILES)                                      &
                                   ! OUT EI for land tiles
     &,Q1P5M_TILE(LAND_POINTS,NTILES)                                   &
!                               ! OUT Q1P5M over land tiles.
     &,T1P5M_TILE(LAND_POINTS,NTILES)                                   &
!                               ! OUT T1P5M over land tiles.
     &,ECAN_TILE(LAND_POINTS,NTILES)                                    &
                                ! OUT ECAN for land tiles
     &,MELT_TILE(LAND_POINTS,NTILES)
!                               ! OUT Snowmelt on tiles (kg/m2/s).

! Local variables for MOSES II only
      REAL                                                              &
     &  E_SSI(ROW_LENGTH,ROWS)                                          &
                                    ! Evaporation from mean sea
     &, EI_SICE(ROW_LENGTH,ROWS)                                        &
                                    ! Output from sf_evap.
     &, SURF_HT_FLUX_LAND(ROW_LENGTH,ROWS)                              &
     &, SURF_HT_FLUX_SICE(ROW_LENGTH,ROWS)                              &
     &, FTL_SSI(ROW_LENGTH,ROWS)    ! Surface FTL for mean sea
! local variables for mineral dust

      REAL DUST_ALL(ROW_LENGTH,ROWS,MODEL_LEVELS,NDIV) !dust mmr
      REAL T_MODELLEVS(ROW_LENGTH,ROWS,MODEL_LEVELS) !T on model levs


      Integer                                                           &
     &  Error_code


! local variables.
      REAL                                                              &
     &  DENOM                                                           &
                   ! Denominator in PC2 inhomogeneous ice forcing calc.
     &, Q4                                                              &
                   ! QCF increment in PC2 inhomog.    ice forcing calc.
     &, w_int(row_length, rows, bl_levels)
                   ! w on interior points     
!
! loop counters
      Integer                                                           &
     &  i, j, k                                                         &
     &, kinvert                                                         &
                                ! vertical index for inverted arrays.
     &, IDIV

! Diagnostic switches
! a) boundary layer
      Logical                                                           &
     &  su10, sv10, slh, sq1p5, sT1p5, sq_T1p5                          &
     &, simlt, smlt, l_ftl, l_fqw, l_taux, l_tauy

! local variables

! CLSCHM3A defines reference numbers for cloud schemes in two-stream
! radiation code.

      ! maximum/random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_MAX=2

      ! random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_RANDOM=4

      ! maximum overlap in a column model
      INTEGER,PARAMETER:: IP_CLOUD_COLUMN_MAX=3

      ! clear column
      INTEGER,PARAMETER:: IP_CLOUD_CLEAR=5

      ! mixed column with split between  convective and layer cloud.
      INTEGER,PARAMETER:: IP_CLOUD_TRIPLE=6

      ! Coupled overlap with partial correlation of cloud
      INTEGER,Parameter:: IP_cloud_part_corr=7

      ! Coupled overlap with partial correlation of cloud
      ! with a separate treatment of convective cloud
      INTEGER,Parameter:: IP_cloud_part_corr_cnv=8

! CLSCHM3A end
!         clschm3a defines reference numbers for cloud overlap choices.
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
!
!       Specified in-plume value of ice content when there is not enough
!       information to calculate this parameter another way (kg kg-1)
      REAL,PARAMETER:: LS_BL0=1.0E-4
!
!
! Diagnostics controlled by Diagnostic switches

      Real                                                              &
     &  q1p5m(row_length, rows)                                         &
     &, t1p5m(row_length, rows)                                         &
     &, rho1(row_length, rows)                                          &
                                               ! Density at level 1
     &, u10m(row_length, rows)                                          &
     &, v10m(row_length, n_rows)                                        &
     &, latent_heat(row_length, rows)                                   &
     &, sice_mlt_htf(row_length, rows, nice)                            &
                                               !output seaice topmelt
     &, sea_ice_htf(row_length, rows, nice)                             &
                                               !output seaice botmelt
     &, ti_gb(row_length, rows)                !output seaice temp.

      Real                                                              &
     &  rhokh_mix (row_length, rows, bl_levels)

! Additional variables for SCM diagnostics which are dummy in full UM
      Integer                                                           &
     &  nSCMDpkgs              ! No of SCM diagnostics packages

      Logical                                                           &
     &  L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages

! local variables
      Integer                                                           &
     &  nclds      ! Number of radiation cloud levels ( <=  wet levels)
!
!-------Needed for area_cloud formulation-------------------------------
      Integer                                                           &
     &  levels_per_level                                                &
                               ! 3 is hardwired inside ls_arcld
     &, large_levels           ! depends on above and wet_model_levels
!
!-----------------------------------------------------------------------

!     EAK
      REAL                                                              &
     &  alb_tile(land_points,ntiles,4)                                  &
     &, surf_down_sw(row_length,rows,4)                                 &
     &, cos_zenith_angle(row_length,rows)                               &
     &, lat(row_length, rows)                                           &
                                          ! Lat. of gridpoint chosen
     &, long(row_length, rows)                                          &
                                          ! Long. of gridpoint chosen
     &, time_sec                                                        &
                                ! actual time of day in secs.
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
     &,Z1_UV(row_length,rows)                                       &
                                ! Height of lowest u,v level.
     &, HCONS(land_points)                                          &
     &, CPOOL_TILE(land_points,NTILES,10)                           &
     &, NPOOL_TILE(land_points,NTILES,10)                           &
     &, PPOOL_TILE(land_points,NTILES,12)                           &
     &, SOIL_ORDER(land_points)                                     &
     &, GLAI(land_points,NTILES)                                    &
     &, PHENPHASE(land_points,NTILES)                               &
     ! Lestevens 23apr13
     &, NPP_FT_ACC(land_points,NTILES)                              &
     &, RESP_W_FT_ACC(land_points,NTILES)

      Real                                                              &
     &  clapp(land_points)                                              &
                             !  qrparm.soil.bwag ?
!                               Clapp-Hornberger exponent.
     &, satcon(land_points)                                             &
                                 !  qrparm.soil.satcon
     &, sathh(land_points)                                              &
                             !  soil water suction
     &, hcon (land_points)                                              &
                             ! soil/qrparm.soil.hcond
     &, hcap (land_points)                                              &
                             ! soil/qrparm.soil.hcap
     &, slm(land_points, dsm_levels)
                             !  qrclim.smc_pm(lev).(month)

!  EAK    diagnostic variables for CABLE output
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
     &,CANOPY_GB(LAND_POINTS)                                           &
     &,TRANSP_TILE(LAND_POINTS,NTILES)

      Integer                                                           &
     &  day                                                             &
     &, total_nsteps                                                    &
                                ! Total number of steps in run
     &, SOIL_TYPE(row_length,rows)                                      &
     &, VEG_TYPE(row_length,rows)                                       &
     &, ISNOW_FLG3L(land_points,NTILES)

      Logical                                                           &
     &  l_tile_pts(land_points,ntiles)

      LOGICAL                                                           &
     & l_cable

! Local data arrays
      Logical                                                           &
     &  ext_LAND(0:rhc_row_length+1,0:rhc_rows+1)
!
      Real                                                              &
     &  T(row_length, rows, bl_levels)                                  &
     &, RHCPT(rhc_row_length, rhc_rows, wet_model_levels)               &
     &, ext_FLANDG(0:rhc_row_length+1,0:rhc_rows+1)

      Real                                                              &
     &  u_inc_bl(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           bl_levels)                                             &
     &, v_inc_bl(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,        &
     &           bl_levels)

      Real                                                              &
     &  interp(row_length, rows)                                        & 
                                 ! Workspace in calculation of Ch_term's
     &, temp1, temp2                                                    &
                                 ! Temporary variables in calculation of
                                 ! Ch_term's.
     &, drydep_str(row_length, rows)                                    &
     &, ext_p_layer_centres(0:rhc_row_length+1,0:rhc_rows+1,            &
     &                                         0:wet_model_levels)      &
     &, ext_TL(0:rhc_row_length+1, 0:rhc_rows+1, wet_model_levels)      &
     &, ext_QL(0:rhc_row_length+1, 0:rhc_rows+1, wet_model_levels)      &
     &, ext_QCF(0:rhc_row_length+1,0:rhc_rows+1, wet_model_levels)      &
     &, ext_ICE_FRAC(0:rhc_row_length+1,0:rhc_rows+1)                   &
     &, work2d_1(row_length, rows)                                      &
                                    ! Single-level work array (cloud)
     &, plsp(row_length, rows)      ! Probability of large-scale precip

!
! Allocatable arrays for diagnostic variables - required to save memory
! when diagnostic not requested
      Real,Dimension(:,:,:),Allocatable::                               &
     & combined_cloud                                                   &
                             ! Mixed CCA and CF per gridbox
     &,T_earliest                                                       &
     &,q_earliest                                                       &
     &,qcl_earliest                                                     &
     &,qcf_earliest                                                     &
     &,cf_earliest                                                      &
     &,cfl_earliest                                                     &
     &,cff_earliest                                                     &
     &,T_inc_PC2                                                        &
                        !  temperature     increment due to PC2 homog
     &,q_inc_PC2                                                        &
                        !  humidity        increment due to PC2 homog
     &,qcl_inc_PC2                                                      &
                        !  qCL             increment due to PC2 homog
     &,cfl_inc_PC2                                                      &
                        !  cf_liquid       increment due to PC2 homog
     &,bcf_inc_PC2      !  bulk cloud      increment due to PC2 homog
!
      Real, Dimension (:,:,:), Allocatable ::    &
       zeros              & ! Array of zero values
     , TL_force           & ! Forcing of TL by homogenous processes
     , QT_force           & ! Forcing of QT by homogenous processes
     , ccw_cca            & ! Convective cloud water * frac (i.e. gridbox
                            ! mean ccw)
     , cca_3d               ! 3D array of convective cloud frac
!
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
     &,L_qcf_incr_bl_lsc                                                &
                             ! qCF across BL (and LS CLD)
     &,L_apply_diag                                                     &
                             ! flag to determine when to apply
                             ! diagnostics when iterating
     &,L_cf_incr_bl                                                     &
                            ! total cloud fraction
     &,L_cfl_incr_bl                                                    &
                            ! liquid cloud fraction
     &,L_cff_incr_bl        ! frozen cloud fraction
!
! Switch for field calculations to support STASH diagnostics
      Logical                                                           &
     & L_combi_cld           ! combined cloud amount
      Logical                                                           &
     & last_substep           ! true when last phys2 substep

!
! External Routines:
      External ls_calc_rhcrit
      External ls_arcld, swap_bounds
      External flux_diag, tr_mix, bl_lsp
      External R2_calc_total_cloud_cover
      External timer
      External BL_TRACER_INTCTL
      EXTERNAL GRAVSETT
      EXTERNAL LS_ACF_Brooks
      External diagnostics_bl, diagnostics_lscld
!- End of Header
! ----------------------------------------------------------------------
!
      If ( error_code  ==  0) Then
! ----------------------------------------------------------------------
! Section BL.0 Initialisation of variables.
! ----------------------------------------------------------------------

        If (Substep_Number == Num_Substeps) last_substep = .TRUE.

        IF ( TRWEIGHTS1  ==  ON ) THEN
!         ! Set all implict weights used by tracers to one
!         !  - overweighting not necessary since tracers
!         !    have no feedback on the diffusion coefficients
          DO K=1,BL_LEVELS
            alpha_tr(k) = 1.0
          ENDDO
        ELSE
!        ! Set implict weights used by tracers to those input
          DO K=1,BL_LEVELS
            alpha_tr(k) = alpha_cd(k)
          ENDDO
        ENDIF

! Apply diags at last cycle only
        L_apply_diag = CycleNo == NumCycles

! Set diagnostic flags required for boundary layer diagnostics from
! STASHflags.
!        !--------------------------------------------------------
!        ! Note that an equivalent block of code exists in routine
!        ! ni_bl_ctl, and needs to be kept consistent.
!        !--------------------------------------------------------
!        ! Windspeed (227) and u, v at 10m on 'B' or 'C' grid
         su10 = (sf(227,3) .or. sf(225,3) .or. sf(209,3) .or. sf(463,3))&
     &           .and.  L_apply_diag
         sv10 = (sf(227,3) .or. sf(226,3) .or. sf(210,3) .or. sf(463,3))&
     &           .and.  L_apply_diag
         slh = sf(234,3) .and. L_apply_diag
         sq_T1p5 = ( sf(236,3) .or. sf(237,3) .or. sf(245,3)            &
     &        .or. sf(247,3) .or. sf(248,3) .or. sf(250,3) .or. L_scrn  &
     &        .or. sf(341,3) .or. sf(342,3)                             &
     &        .or. sf(253,3) .or. sf(328,3) .or. sf(329,3)              &
     &              ) .and. L_apply_diag
         sq1p5 = sq_T1p5 .and. L_apply_diag
         sT1p5 = sq_T1p5 .and. L_apply_diag
         simlt = sf(235,3) .and. L_apply_diag
         smlt = sf(258,3) .and. L_apply_diag
         smlt = sf(258,3) .and. last_substep .and. L_apply_diag
         l_ftl = sf(216,3) .and. last_substep .and. L_apply_diag
         l_fqw = sf(222,3) .and. last_substep .and. L_apply_diag
         l_taux = ( sf(219,3) .or. sf(221,3) ) .and. last_substep       &
     &           .and. L_apply_diag
         l_tauy = ( sf(220,3) .or. sf(221,3) ) .and. last_substep       &
     &           .and. L_apply_diag
         L_u_incr_bl = sf(185,3) .and. L_apply_diag
         L_v_incr_bl = sf(186,3) .and. L_apply_diag
         L_T_incr_bl_lsc = (sf(181,9).or.sf(181,3)) .and. L_apply_diag
         L_Tl_incr_bl_lsc = sf(189,3) .and. L_apply_diag
         L_q_incr_bl_lsc = (sf(182,9).or.sf(182,3)) .and. L_apply_diag
         L_qtl_incr_bl_lsc = sf(190,3) .and. L_apply_diag
         L_qcl_incr_bl_lsc =(sf(183,9).or.sf(183,3)).and. L_apply_diag
         L_qcf_incr_bl_lsc = sf(184,3) .and. L_apply_diag
         L_cf_incr_bl  = sf(192,3) .and. L_apply_diag
         L_cfl_incr_bl = sf(193,3) .and. L_apply_diag
         L_cff_incr_bl = sf(194,3) .and. L_apply_diag
!
! Flag required for pre-calculation of cloud-related fields needed for
! cloud and vis diagnostics. Not duplicated in ni_bl_ctl.
         L_combi_cld = ( L_plsp .OR.                                    &
     &           sf(208,9) .OR. sf(209,9) .OR. sf(210,9) .OR. sf(211,9) &
     &      .OR. sf(212,9) .OR. sf(213,9) .OR. sf(214,9) .OR. sf(215,9) &
     &      .OR. sf(216,9) .OR. sf(217,9) .OR. sf(223,9) .OR. sf(231,9) &
     &      ) .AND. L_apply_diag

!---------------------------------------------------------------------
! Intercept values of physics increments before being updated by
! implicit solver for optional output of bl wind increments
!---------------------------------------------------------------------
! We need to store information about the increments of the temperature
! and moisture variables, so copy these to the _earliest variables.
      If (L_T_incr_bl_lsc .OR. L_Tl_incr_bl_lsc .OR. L_pc2) Then
!
        Allocate ( T_earliest(row_length,rows,model_levels) )

! Hold initial value of Temperature
        Do k=1, model_levels
          Do j=1, rows
            Do i=1, row_length
              T_earliest(i,j,k) = T_latest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

      End if                   ! on STASHflags or PC2
!
      If(L_q_incr_bl_lsc .OR. L_qtl_incr_bl_lsc                         &
     &   .OR. L_qcl_incr_bl_lsc .OR. L_Tl_incr_bl_lsc                   &
     &   .OR. L_qcf_incr_bl_lsc .OR. L_pc2) Then
!
        Allocate ( q_earliest(row_length,rows,wet_model_levels) )
        Allocate ( qcl_earliest(row_length,rows,wet_model_levels) )
        Allocate ( qcf_earliest(row_length,rows,wet_model_levels) )
        Allocate ( cf_earliest(row_length,rows,wet_model_levels) )
        Allocate ( cfl_earliest(row_length,rows,wet_model_levels) )
        Allocate ( cff_earliest(row_length,rows,wet_model_levels) )
!
! Hold initial values of wet parameters
        Do k=1, wet_model_levels
          Do j=1, rows
            Do i=1, row_length
              q_earliest(i,j,k)   = q_latest(i,j,k)
              qcl_earliest(i,j,k) = qcl_latest(i,j,k)
              qcf_earliest(i,j,k) = qcf_latest(i,j,k)
              cf_earliest(i,j,k)  = cf_latest(i,j,k)
              cfl_earliest(i,j,k) = cfl_latest(i,j,k)
              cff_earliest(i,j,k) = cff_latest(i,j,k)
            End Do ! i
          End Do ! j
        End Do ! k

      End if                  ! on STASHflags or PC2
!
! ----------------------------------------------------------------------
! Section BL.1 Calculate T at old time level.
! Modified to use latest values to avoid time-level inconsistencies
! with cloud data.
! ---------------------------------------------------------------------

        Do k = 1, bl_levels
          Do j = 1, rows
            Do i = 1, row_length
              T(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
            End Do
          End Do
        End Do

!  Initialise output arrays to zero

        Do j= 1, rows
          Do i= 1, row_length
            ecan(i,j)=0.0
            ei(i,j) = 0.0
            snowmelt(i,j)=0.0
          End Do
        End Do

! ----------------------------------------------------------------------
! Section BL.2  Call Implicit solver
!               Call tracer mixing for qcf.
!               Call tracer mixing for other tracers.
! ----------------------------------------------------------------------

! DEPENDS ON: imps_intct
        Call IMPS_INTCT( halo_i, halo_j, off_x, off_y, row_length,      &
     &                  rows, n_rows, global_row_length, proc_row_group,&
     &                  at_extremity, n_proc, n_procx, n_procy,         &
     &                  neighbour,                                      &

     &                  ntiles, land_points, nice,                      &

! in values defining model domain, vertical grid of model atmosphere
!    and implicit weights. :
     &                  model_domain, bl_levels,                        &
     &                  r_rho_levels, r_theta_levels,                   &
     &                  alpha_cd,                                       &

! IN Substepping Information
     &              Substep_Number, Num_Substeps, L_phys2_substep,      &

! IN New BL solver parameters
     &              L_us_blsol, Puns, Pstb,                             &

! IN U and V momentum fields.
     &              u, v, u_conv, v_conv,                               &
! IN Non turbulent increments to momentum and thermodynamic variables.
!  (New dynamics only).
     &                  R_u, R_v,                                       &
! in soil/vegetation/land surface data :
     &                  land_sea_mask, land_index,                      &
     &                  dst_levels, dsm_levels,                         &

! in sea/sea-ice data :
     &                  di,ice_fract, di_ncat, ice_fract_ncat, u_0, v_0,&

! in cloud data :
! cloud_fraction passed in
     &              q, qcf, qcl, q_conv, qcf_conv, qcl_conv,            &
     &                  qcf_latest, qcl_latest, T, T_conv,              &
!
! in everything not covered so far :
     &                  p_star, surf_radflux, timestep, L_sice_heatflux,&

! in variables required in IMP_SOLVER
     &                  alpha1_sice, ashtf,                             &
     &                  dtrdz_charney_grid, rdz_charney_grid,           &
     &                  dtrdz_u, dtrdz_v, rdz_u, rdz_v,                 &
     &                  cdr10m_u, cdr10m_v, z1_tq,                      &
     &                  rhokm_u,rhokm_v,                                &
! in variables for new BL solver                                
     &                  bl_type_1,bl_type_2,bl_type_3,bl_type_4,        &
     &                  bl_type_5,bl_type_6,bl_type_7,                  &

! IN additional variables for MOSES II
     &                  TILE_PTS,TILE_INDEX,TILE_FRAC,CANOPY,           &
     &                  L_NEG_TSTAR,ALPHA1,FRACA,RHOKH_TILE,SMC,CHR1P5M,&
     &                  RESFS,Z0HSSI,Z0MSSI,CANHC_TILE,FLAKE,           &
     &                  WT_EXT_TILE,LW_DOWN,SW_TILE,ASHTF_TILE,         &
     &                  FQT_ICE,FTL_ICE,RESFT,RHOKH_SICE,               &
     &                  RHOKPM,RHOKPM_POT,RHOKPM_SICE,                  &
     &                  Z0H_TILE,Z0M_TILE,CHR1P5M_SICE,                 &
     &                  FLAND,FLANDG,FLANDG_U,FLANDG_V,TSTAR_SEA,       &
     &                  ANTHROP_HEAT,                                   &

!    EAK
!    IN
     &  l_cable                                                         &
     &, surf_down_sw,alb_tile,cos_zenith_angle,l_tile_pts               &
!     &, surf_down_sw,alb_tile,cos_zenith_angle               &
     &, lat,long,day,time_sec                                           &
     &, ls_rain, ls_snow, conv_rain, conv_snow                          &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L                   &
     &, SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE             &
     &, TSOIL_TILE,T_SURF_TILE,HCONS,clapp                              &
     &, SATHH,SATCON,HCAP,HCON,Z1_UV                                    &
     &, smvccl, smvcwt, smvcst, sthf, sthu                              &
     &, slm                                                             &
     &, SOIL_TYPE,VEG_TYPE                                              &
     &, ISNOW_FLG3L,total_nsteps                                        &
     &, FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB                         &
     &, TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB                     &
     &, USTAR_CAB,SURF_HTF_CAB                                          &
     &, TOT_ALB                                                         &
     &, SNAGE_TILE,RTSOIL_TILE                                          &
     &, GFLUX_TILE,SGFLUX_TILE                                          &
     &, GC,GS,CANOPY_GB                                                 & ! Added GC, pfv 24oct13
!    &, GS,CANOPY_GB                                                    & ! orig,     pfv 24oct13
! Lestevens March 2010: CO2 fluxes and dimensions to pass to cable_IMum
     &, DIM_CS1, DIM_CS2    &
     &, NPP, NPP_FT         &
     &, GPP, GPP_FT         &
     &, RESP_S, RESP_S_TOT  &
     &, RESP_S_TILE         &  !kdcorbin, 10/10 
     &, RESP_P, RESP_P_FT   &
     &, G_LEAF              &  !kdcorbin, 10/10
     &, TRANSP_TILE         &  !Lestevens, sep12
     &, CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER                     &
     &, GLAI,PHENPHASE                                                  &
     ! Lestevens 23apr13
     &, NPP_FT_ACC,RESP_W_FT_ACC,                                       &

! INOUT data :
     &                  t_soil, ti, ti_gb, t_surf,                      &
     &                  T_latest, Q_latest,                             &

! inout  diagnostic started in bdy_layr not requiring stash flags :
     &                  e_sea, fqT, ftl, h_sea, rhokh,                  &
     &                  taux, tauy,                                     &

! INOUT additional variables for MOSES II
     &                  TSTAR_TILE,FQT_TILE,EPOT_TILE,FTL_TILE,         &
     &                  SNOW_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,OLR,  &
     & TSTAR_SICE,TSTAR_SSI,                                            &
     & TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI,                           &

! out u and v increments
     &                  u_inc_bl, v_inc_bl,                             &

! out  diagnostic not requiring stash flags :
     &                  rhokh_mix,                                      &
     &                  sea_ice_htf, surf_ht_flux_gb,                   &
     &                  SURF_HT_FLUX_LAND,SURF_HT_FLUX_SICE,            &

! out diagnostic requiring stash flags :
     &                  sice_mlt_htf, snomlt_surf_htf,                  &
     &                  latent_heat, q1p5m, t1p5m, u10m, v10m,          &

! (in) stash flags :-
     &              simlt, smlt, slh, sq1p5, st1p5, su10, sv10,         &
     &              l_ftl, l_fqw, l_taux, l_tauy,                       &

! OUT additional variables for MOSES II
     &                  ESOIL_TILE,ES,EI_TILE,                          &
     &                  Q1P5M_TILE,T1P5M_TILE,ECAN_TILE,MELT_TILE,      &
     &                  TSTAR_LAND,E_SSI,EI_SICE,FTL_SSI,               &
     &                  ERROR_CODE,                                     &

! out data required elsewhere in um system :
     &                  ecan, ei, ext, snowmelt,                        &
     &                  lq_mix_bl,                                      &
     &                  L_flux_bc, Ltimer )

! add boundary layer increment to R_u and R_v
        If ( l_u_incr_bl ) Then
! add boundary layer increment to R_u and R_v
        Do k = 1, bl_levels
          Do j = 1, rows
            Do i = 1, row_length
              u_incr_diag_bl(i,j,k) = u_incr_diag_bl(i,j,k)             &
     &                              + ( u_inc_bl(i,j,k) - R_u(i,j,k) )
              R_u(i,j,k) = u_inc_bl(i,j,k)
            End Do
          End Do
        Enddo
        Else
        Do k = 1, bl_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_u(i,j,k) = u_inc_bl(i,j,k)
            End Do
          End Do
        Enddo
        Endif

        If ( l_v_incr_bl ) Then
        Do k = 1, bl_levels
          Do j = 1, n_rows
            Do i = 1, row_length
              v_incr_diag_bl(i,j,k) = v_incr_diag_bl(i,j,k)             &
     &                              + ( v_inc_bl(i,j,k) - R_v(i,j,k) )
              R_v(i,j,k) = v_inc_bl(i,j,k)
            End Do
          End Do
        End Do
        Else
        Do k = 1, bl_levels
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = v_inc_bl(i,j,k)
            End Do
          End Do
        End Do
        Endif

! Call tr_mix to mix qcf
! output qcf_flux in T to save workspace
! Pass in a zero field for source terms.
        Do j = 1, rows
          Do i = 1, row_length
            interp(i,j) = 0.
          End Do
        End Do

! DEPENDS ON: tr_mix
        Call  tr_mix (                                                  &
     &        halo_i, halo_j, row_length, rows, bl_levels, off_x, off_y &
     &        , alpha_tr                                                &
     &       ,rhokh_mix(1,1,2), rhokh_mix(1,1,1)                        &
     &       ,dtrdz_charney_grid, r_rho_levels, r_theta_levels          &
     &       ,timestep                                                  &
     &       ,T, qcf_latest                                             &
     &       ,interp, interp, drydep_str                                &
     &       ,KENT, WE_LIM, T_FRAC, ZRZI                                &
     &       ,KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC                &
     &       ,ZH ,ZHSC, Z_HALF                                          &
     &       ,error_code, .false.                                       &
     &        )
!
!Call tr_mix to mix w in the vertical when the subgrid turbulence
!scheme is activated
!
        If (L_subfilter_vert .or. L_subfilter_horiz) then

          Do k = 1, bl_levels
            Do j = 1, rows
              Do i = 1, row_length
                w_int(i,j,k) = w(i,j,k)
              End Do
            End Do
          End Do
!
!rhokm should be used in mixing of w.  rhokh_mix is used instead to avoid
!passing around rhokm on p-points.
!
! DEPENDS ON: tr_mix
          Call  tr_mix (                                                &
     &        halo_i, halo_j, row_length, rows, bl_levels, off_x, off_y &
     &       ,alpha_cd                                                  &
     &       ,rhokh_mix(1,1,2), rhokh_mix(1,1,1)                        &
     &       ,dtrdz_charney_grid, r_rho_levels, r_theta_levels          &
     &       ,timestep                                                  &
     &       ,T, w_int                                                  &
     &       ,interp, interp, drydep_str                                &
     &       ,KENT, WE_LIM, T_FRAC, ZRZI                                &
     &       ,KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC                &
     &       ,ZH ,ZHSC, Z_HALF                                          &
     &       ,error_code, .false.                                       &
     &        )

          Do k = 1, bl_levels
            Do j = 1, rows
              Do i = 1, row_length
                R_w(i,j,k) = R_w(i,j,k) + (w_int(i,j,k) - w(i,j,k))
              End Do
            End Do
          End Do

        End If     !L_subfilter_vert or L_subfilter_horiz

! apply BL tracer mixing and gravitational settling of dust
! on the last cycle only

! pfv, 29oct13a, access 1.4, put GRAVSETT before BL_TRACER_INTCTL
! (via Stephanie Woodward)

!     ---------------------------------
      If ( CycleNo == NumCycles ) Then
!     ---------------------------------
!
! Gravitational settling of mineral dust
!
!      ----------------
       IF (L_DUST) THEN
!      ----------------

         DO K = 1, MODEL_LEVELS
           DO J = 1,ROWS
             DO I = 1,ROW_LENGTH
               T_MODELLEVS(I,J,K)=THETA(I,J,K)*EXNER_THETA_LEVELS(I,J,K)
               DUST_ALL(I,J,K,1)=DUST_DIV1(I,J,K)
               DUST_ALL(I,J,K,2)=DUST_DIV2(I,J,K)
               DUST_ALL(I,J,K,3)=DUST_DIV3(I,J,K)
               DUST_ALL(I,J,K,4)=DUST_DIV4(I,J,K)
               DUST_ALL(I,J,K,5)=DUST_DIV5(I,J,K)
               DUST_ALL(I,J,K,6)=DUST_DIV6(I,J,K)
             ENDDO !ROW_LENGTH
           ENDDO !ROWS
         ENDDO !MODEL_LEVELS

         DO IDIV=1,NDIV
!
! DEPENDS ON: gravsett
           CALL GRAVSETT(                                               &
     & ROW_LENGTH,ROWS,MODEL_LEVELS,                                    &
     & DUST_ALL(1:ROW_LENGTH,1:ROWS,1:MODEL_LEVELS,IDIV),               &
     & DREP(IDIV),RHOP,P_LAYER_CENTRES,P_LAYER_BOUNDARIES,T_MODELLEVS,  &
     & TIMESTEP,NUM_SUBSTEPS,SUBSTEP_NUMBER,                            &
     & DRYDEP2(1:ROW_LENGTH,1:ROWS,IDIV),L_CAM_DUST                     &
     & )
!
         ENDDO
!
         DO K = 1, MODEL_LEVELS
           DO J = 1,ROWS
              DO I = 1,ROW_LENGTH
                DUST_DIV1(I,J,K)=DUST_ALL(I,J,K,1)
                DUST_DIV2(I,J,K)=DUST_ALL(I,J,K,2)
                DUST_DIV3(I,J,K)=DUST_ALL(I,J,K,3)
                DUST_DIV4(I,J,K)=DUST_ALL(I,J,K,4)
                DUST_DIV5(I,J,K)=DUST_ALL(I,J,K,5)
                DUST_DIV6(I,J,K)=DUST_ALL(I,J,K,6)
              ENDDO !ROW_LENGTH
            ENDDO !ROWS
          ENDDO !MODEL_LEVELS

!      ----------------
       ENDIF !L_DUST
!      ----------------

! Mixing for all non-qcf tracers done in subroutine

!       ------------------------------------------------------
        If ( l_bl_tracer_mix .OR. l_sulpc_so2 .OR. l_soot .OR.          &
     &       l_co2_interactive .OR. l_murk .OR. l_biomass .OR.          &
     &       l_dust .OR. l_ocff) Then
!       ------------------------------------------------------

! DEPENDS ON: bl_tracer_intctl
          Call BL_TRACER_INTCTL(                                        &
     &         row_length, rows, bl_levels, off_x, off_y, at_extremity  &
     &        ,dtrdz_charney_grid,r_rho_levels,r_theta_levels           &
     &        ,halo_i,halo_j                                            &
     &        ,tr_levels, tr_vars, model_levels                         &
     &        ,DIM_CS2                                                  &
     &        ,alpha_tr, rhokh_mix                                      &
     &        ,p_star, p, timestep                                      &
! Control logicals
     &        ,L_MURK,L_MURK_ADVECT,L_BL_TRACER_MIX,L_DUST,L_CAM_DUST   &
     &        ,L_SULPC_SO2, l_sulpc_nh3, l_sulpc_dms, l_soot, l_biomass &
     &        ,l_ocff, l_co2_interactive                                &
     &        ,l_co2_emits                                              &
! rml 1/7/13
     &        ,L_CO2_TRACER                                             &
     &        ,L_ukca, l_use_cariolle                                   &
! Fields to mix
     &        ,aerosol, free_tracers, ozone_tracer                      &
! Mineral dust
     &        ,DUST_DIV1,DUST_DIV2,DUST_DIV3                            &
     &        ,DUST_DIV4,DUST_DIV5,DUST_DIV6                            &
! Sulphur cycle
     &        ,so2, dms, so4_aitken, so4_accu, so4_diss, nh3            &
! Soot cycle
     &        ,soot_new, soot_aged, soot_cld                            &
! Biomass aerosol
     &        ,bmass_new, bmass_agd, bmass_cld                          &
! Fossil-fuel organic carbon aerosol
     &        ,ocff_new, ocff_aged, ocff_cld                            &
! Carbon cycle
     &        ,co2, co2_emits, co2flux, npp, resp_s_tot                 &
     &        ,co2_flux_tot, land_co2                                   &
! Tracer fluxes - kdcorbin, 05/10
     &        ,tracer_flux1, tracer_flux2, tracer_flux3, tracer_flux4   &
     &        ,tracer_flux5, tracer_flux6, tracer_flux7, tracer_flux8   &
     &        ,tracer_flux9, tracer_flux10,tracer_flux11,tracer_flux12  &
     &        ,tracer_flux13,tracer_flux14,tracer_flux15,tracer_flux16  &
     &        ,tracer_flux17,tracer_flux18,tracer_flux19,tracer_flux20  &
! Emissions fields
     &        ,DUST_FLUX                                                &
     &        ,so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em   &
     &        ,ocff_hilem, ocff_em, drydep_str                          &
     &        ,kent, we_lim, t_frac, zrzi                               &
     &        ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc               &
     &        ,zh, zhsc, z_half                                         &
! IN for dry deposition of tracers
     &,       rho_aresist, aresist, resist_b                            &
     &,        R_B_DUST,T_SURF                                          &
     &,       land_points, land_index, ice_fract                        &
! IN variables for MOSES II only
     &,       ntype, ntiles, tile_pts, tile_index, tile_frac            &
     &,       canopy, catch, snow_tile, gc                              &
     &,       aresist_tile, resist_b_tile, FLANDG,                      &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &        STASHwork3                                                &
     &        )

!       ------------------------------------------------------
        End If  ! if tracer mixing required
!       ------------------------------------------------------

!     ---------------------------------
      End If ! (CycleNo == NumCycles)
!     ---------------------------------

! ----------------------------------------------------------------------
! Section BL.3 Form Ch terms.
! ----------------------------------------------------------------------

! NB: Code assumes that the number of boundary layer levels is at
!     least 2.
        If ( L_use_bl_diag_term .and. last_substep ) Then
! Remove factor of r**2 from rho, store in interp

        Do k = 1, bl_levels

          Do j = 1, rows
            Do i =1, row_length
              interp(i,j) = ( ( r_theta_levels(i,j,k) -                 &
     &                          r_rho_levels(i,j,k) )                   &
     &                         * rho(i,j,k+1) /                         &
     &                          (r_rho_levels(i,j,k+1) *                &
     &                           r_rho_levels(i,j,k+1))                 &
     &                       + ( r_rho_levels(i,j,k+1) -                &
     &                           r_theta_levels(i,j,k) )                &
     &                         * rho(i,j,k) /                           &
     &                          (r_rho_levels(i,j,k) *                  &
     &                           r_rho_levels(i,j,k))                   &
     &                      ) / (r_rho_levels(i,j,k+1) -                &
     &                           r_rho_levels(i,j,k) )
            End Do
          End Do

          If (k  ==  1) Then
            Do j = 1, rows
              Do i =1, row_length

                temp1 = rhokh(i,j,k)

                temp2 = rhokh(i,j,k+1)/ (r_theta_levels(i,j,k+1)        &
     &                                   - r_theta_levels(i,j,k))

                Ch_term(i,j,k) = exner_theta_levels(i,j,k) * alpha_Cd(k)&
     &                           * ( temp1 + temp2 )                    &
     &                           / (interp(i,j) *                       &
     &                           (r_rho_levels(i,j,k+1) -               &
     &                            r_rho_levels(i,j,k)))

              End Do
            End Do


          Else If (k  <   bl_levels) Then
            Do j = 1, rows
              Do i =1, row_length

                temp1 = rhokh(i,j,k)/ (r_theta_levels(i,j,k)            &
     &                                 - r_theta_levels(i,j,k-1))

                temp2 = rhokh(i,j,k+1)/ (r_theta_levels(i,j,k+1)        &
     &                                   - r_theta_levels(i,j,k))

                Ch_term(i,j,k) = exner_theta_levels(i,j,k) * alpha_Cd(k)&
     &                           * ( temp1 + temp2 )                    &
     &                         / (interp(i,j) * (r_rho_levels(i,j,k+1)  &
     &                            - r_rho_levels(i,j,k)))

              End Do
            End Do

          Else ! top boundary layer level

            Do j = 1, rows
              Do i =1, row_length

                temp1 = rhokh(i,j,k)/ (r_theta_levels(i,j,k)            &
     &                                 - r_theta_levels(i,j,k-1))

                temp2 = 0.0

                Ch_term(i,j,k) = exner_theta_levels(i,j,k) * alpha_Cd(k)&
     &                           * ( temp1 + temp2 )                    &
     &                         / (interp(i,j) * (r_rho_levels(i,j,k+1)  &
     &                             - r_rho_levels(i,j,k)))

              End Do
            End Do

          End If

        End Do
! swop halo values for Ch_term in boundary layer
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   cH_term, row_length, rows,                     &
     &                   bl_levels, off_x, off_y, fld_type_p,  .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(CH_term,row_length, rows,              &
     &                     bl_levels,off_x,off_y)

! set zero value above Boundary layer
        Do k = bl_levels + 1, model_levels - 1
          Do j = 1-off_y, rows+off_y
            Do i =1-off_x, row_length+off_x
              Ch_term(i,j,k) = 0.0
            End Do
          End Do
        End Do
        Endif ! If  L_use_bl_diag_term


! ----------------------------------------------------------------------
! Section BL.4 Convert and calculate theta and q fields from qT and Tl.
! ----------------------------------------------------------------------

        If (.not. L_dry) Then
! If the mixed phase precipitation scheme is used then T and Q are
! required to contain T liquid and Q(vapour+liquid) but at this stage
! will actually contain T liquid ice and Q(vapour+liquid+ice)

! DEPENDS ON: bl_lsp
          Call bl_lsp( row_length, rows, bl_levels,                     &
     &                 qcf_latest, q_latest, t_latest )

        End If

        If (L_dry) Then
          Do k = bl_levels+1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_latest(i,j,k) = 0.0
                qcl_latest(i,j,k) = 0.0
                qcf_latest(i,j,k) = 0.0
              End Do
            End Do
          End Do
        End If

! Create Tl and qT outside boundary layer
        Do k = bl_levels+1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              T_latest(i,j,k) = t_latest(i,j,k) -                       &
     &                        (lc * qcl_latest(i,j,k))                  &
     &                         / cp
              q_latest(i,j,k) = q_latest(i,j,k) + qcl_latest(i,j,k)
            End Do
          End Do
        End Do
!
! Prepare for cloud scheme. Are we using PC2 or not?
!
! zero any negative q_latests
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              If (q_latest(i,j,k)  <   0.) Then
!              print*,' neg qT before BL cld call set to zero',i,j,k
                q_latest(i,j,k) = 0.
              End If
              If (qcf_latest(i,j,k)  <   0.) Then
!              print*,' neg qcf before BL cld call set to zero',i,j,k
                qcf_latest(i,j,k) = 0.
              End If
            End Do
          End Do
        End Do
!
! Calculate diagnostic RHcrit or read as parameter in from namelist
!
! Lrhcpt_if1:
        If (L_RHCPT) Then
!       RHCRIT is 3D diagnosed variable
! Wet_mlev_do1:
          Do k = 1, wet_model_levels
! Rhc_rows_do1:
            Do j = 1, rhc_rows
! Rhc_rowlen_do1:
              Do i = 1, rhc_row_length
                ext_p_layer_centres(i,j,k) = p_layer_centres(i,j,k)
                ext_TL(i,j,k) = T_latest(i,j,k)
                ext_QL(i,j,k) = q_latest(i,j,k)
                ext_QCF(i,j,k) = qcf_latest(i,j,k)
              End Do ! Rhc_rowlen_do1
            End Do ! Rhc_rows_do1
          End Do ! Wet_mlev_do1
!
! Rhc_rows_do2:
          Do j = 1, rhc_rows
! Rhc_rowlen_do2:
            Do i = 1, rhc_row_length
              ext_p_layer_centres(i,j,0) = p_layer_centres(i,j,0)
              ext_LAND(i,j) = land_sea_mask(i,j)
              ext_ICE_FRAC(i,j) = ice_fract(i,j)
              ext_FLANDG(i,j) = FLANDG(i,j)
            End Do ! Rhc_rowlen_do2
          End Do ! Rhc_rows_do2
!
! Synchronize haloes.
!
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &       ext_p_layer_centres, rhc_row_length, rhc_rows,             &
     &           wet_model_levels+1, 1, 1, fld_type_p,  .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_p_layer_centres,rhc_row_length,    &
     &                     rhc_rows,                                    &
     &                     wet_model_levels+1,1,1)
!
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     ext_TL, rhc_row_length, rhc_rows,            &
     &                     wet_model_levels, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_TL,rhc_row_length,                 &
     &                     rhc_rows,                                    &
     &                     wet_model_levels,1,1)
!
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     ext_QL, rhc_row_length, rhc_rows,            &
     &                     wet_model_levels, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_QL,rhc_row_length,                 &
     &                     rhc_rows,                                    &
     &                     wet_model_levels,1,1)
!
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     ext_QCF, rhc_row_length, rhc_rows,           &
     &                     wet_model_levels, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_QCF,rhc_row_length,                &
     &                     rhc_rows,                                    &
     &                     wet_model_levels,1,1)
!
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     ext_LAND, rhc_row_length, rhc_rows,          &
     &                     1, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_land,rhc_row_length,               &
     &                     rhc_rows,                                    &
     &                     1,1,1)
!
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     ext_FLANDG, rhc_row_length, rhc_rows,        &
     &                     1, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
          CALL FILL_EXTERNAL_HALOS(ext_FLANDG,rhc_row_length,           &
     &                     rhc_rows,                                    &
     &                     1,1,1)
!
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     ext_ICE_FRAC, rhc_row_length, rhc_rows,      &
     &                     1, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_ICE_FRAC,rhc_row_length,           &
     &                     rhc_rows,                                    &
     &                     1,1,1)
!
! DEPENDS ON: ls_calc_rhcrit
          Call ls_calc_rhcrit( ext_p_layer_centres                      &
!              Array dimensions
     &     , wet_model_levels, rhc_row_length, rhc_rows, bl_levels      &
     &     , global_row_length                                          &
!              Prognostic Fields
     &     , ext_TL, ext_QL, ext_QCF, ext_LAND, ext_FLANDG, ext_ICE_FRAC&
!              Logical control
     &     , lq_mix_bl                                                  &
                           ! Use mixing_ratio
!              Output
     &     , RHCPT)
!
        Else
!       RHCRIT is 1D Parameter read in from namelist
          Do k = 1, wet_model_levels
            RHCPT(1,1,k) = rhcrit(k)
          End Do
        End if  ! Lrhcpt_if1
!
! Which cloud scheme are we using?
!
        If (L_pc2) then
!
! ----------------------------------------------------------------------
! PC2 cloud scheme
! ----------------------------------------------------------------------
          allocate ( zeros   (row_length,rows,wet_model_levels) )
!
! ----------------------------------------------------------------------
! Inhomogenous forcing of ice
! ----------------------------------------------------------------------
!
! Calculate in-plume ice content (LS) by assuming that they are equal to
! the current values, except when the current value is not defined.
! Also calculate the forcing of ice content Q4F
!
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
!
! Calculate Q4. Only perform the calculation if the Q4 is non-zero.
! Since ice is mixed by tracer mixing regardless of whether cumulus
! is present we still need to provide an ice cloud fraction increment
! below cumulus base.
!
                Q4 = qcf_latest(i,j,k) - qcf_earliest(i,j,k)
                IF (Q4  /=  0.0) Then
!
! Calculate the change in total cloud fraction.
! Use a weighted (by ice cloud fraction) average of in-cloud ice
! content and a fixed value to specify the plume ice content.
! The denominator in the deltaCf calculation is then
! (qcf_earliest/cff_earliest + ls_bl0*(1-cff_earliest) - qcf_earliest
! and the qcf_earliest terms cancel.
!
                  DENOM = LS_bl0 * (1.0 - cff_earliest(i,j,k))
!
                  IF ( ABS(DENOM)  >   1.0E-10 ) THEN
!
                    DENOM = Q4 / DENOM
                    cf_latest(i,j,k)  = cf_latest(i,j,k)  +             &
     &                           (1.0 - cf_latest(i,j,k))  * DENOM
                    cff_latest(i,j,k) = cff_latest(i,j,k) +             &
     &                           (1.0 - cff_latest(i,j,k)) * DENOM
!
! Otherwise cloud fraction will go to one. In theory, cloud fraction
! can never be reduced by this process.
!
                  ELSE
                    cf_latest(i,j,k)  = 1.0
                    cff_latest(i,j,k) = 1.0
                  END IF
!
                END IF
!
! Set zero array
                zeros(i,j,k)        = 0.0
              end do
            end do
          end do
!
! ----------------------------------------------------------------------
! Homogenous forcing of the liquid cloud
! ----------------------------------------------------------------------
!
! Calculate forcing in qT and TL. Currently q_latest contains the vapour
! plus liquid content and q_earliest just the initial vapour content
!
          allocate ( QT_force(row_length,rows,wet_model_levels) )
          allocate ( TL_force(row_length,rows,wet_model_levels) )
          allocate ( T_inc_PC2(row_length,rows,wet_model_levels) )
          allocate ( q_inc_PC2(row_length,rows,wet_model_levels) )
          allocate ( qcl_inc_PC2(row_length,rows,wet_model_levels) )
          allocate ( cfl_inc_PC2(row_length,rows,wet_model_levels) )
          allocate ( bcf_inc_PC2(row_length,rows,wet_model_levels) )
!
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                qT_force(i,j,k) = ( q_latest(i,j,k)                     &
     &            - (q_earliest(i,j,k) + qcl_earliest(i,j,k)) )
                TL_force(i,j,k) = ( T_latest(i,j,k)                     &
     &            - (t_earliest(i,j,k)- lc * qcl_earliest(i,j,k) / cp) )
              End Do
            End Do
          End Do
!
! Call homogenous forcing routine
!
! DEPENDS ON: pc2_delta_hom_turb
          Call PC2_DELTA_HOM_TURB(                                      &
! INput variables
     &      p_layer_centres(1,1,1)                                      &
     & ,    wet_model_levels, row_length, rows, timestep                &
! INput variables
     & ,    T_earliest(1,1,1), q_earliest(1,1,1), qcl_earliest(1,1,1)   &
     & ,    cf_latest(1,1,1), cfl_latest(1,1,1), cff_latest(1,1,1)      &
     & ,    TL_force, qT_force, zeros(1,1,1), zeros(1,1,1)              &
! OUTput variables
     & ,    T_inc_PC2, q_inc_PC2, qcl_inc_PC2, bcf_inc_PC2, cfl_inc_PC2 &
! INput variables (other quantities)
     & ,     0.0, 0.0, lq_mix_bl                                        &
     &                            )
!
          If (lcv_pc2_diag_sh) Then
            allocate ( ccw_cca(row_length,rows,wet_model_levels) )
            allocate ( cca_3d (row_length,rows,wet_model_levels) )
           
!
            Do k = 1, wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
!
!
! If below cumulus cloud base or we are using a diagnostic
! shallow convective cloud then we simply zero the liquid water
! content instead of using the homogenous BL response. The forcings
! themselves are still applied but to q and T.
! Note that we find, for PC2, that zeroing up to ntml,
! rather than ntml+1, is more physically justified and gives slightly
! better cloud results.
!
                  If ( (l_pc2_diag_sh_pts(i,j) .and.                     &
                       k <= cct(i,j) ).or. (cumulus(i,j)                 &
                                   .AND. (k  <=  ntml(i,j))) ) Then
                    T_inc_PC2(i,j,k)   =  (-LC/CP) * qcl_earliest(i,j,k)
                    q_inc_PC2(i,j,k)   =  qcl_earliest(i,j,k)
                    qcl_inc_PC2(i,j,k) =  (-qcl_earliest(i,j,k))
                    cfl_inc_PC2(i,j,k) =  (-cfl_earliest(i,j,k))
                    bcf_inc_PC2(i,j,k) =  cff_latest(i,j,k)-cf_earliest(i,j,k)

                  End If

                ! Set convective cloud properties if we are
                ! using a diagnostic shallow convection for PC2
                ! Calculate convective properties

                  If (lcv_3d_cca) then
                  ! Three dimensional cloud field
                    ccw_cca(i,j,k) = ccw(i,j,k) * cca(i,j,k)
                    cca_3d(i,j,k)  = cca(i,j,k)
                  Else
                  ! Two dimensional cloud field
                  ! Will not be correct if mid level convection above shallow
                  ! convection 
                    If (k <= cct(i,j) .and. k >= ccb(i,j)) Then
                    ! Inside convective cloud
                      ccw_cca(i,j,k) = ccw(i,j,k) * cca_2d(i,j)
                      cca_3d(i,j,k)  = cca_2d(i,j)
                    Else
                    ! Outside of convective cloud
                      ccw_cca(i,j,k) = 0.0
                      cca_3d(i,j,k)  = 0.0
                    End if
                  End if  ! lcv_3d_cca

                  If ( l_pc2_diag_sh_pts(i,j) .and. k <= cct(i,j)        &
                                .and. t_earliest(i,j,k) >= tice) Then

                    ! Hand over convective attributes to the large scale 
                    T_inc_PC2(i,j,k)   =  T_inc_PC2(i,j,k)            &
                                        + (LC/CP) *  ccw_cca(i,j,k)
                    q_inc_PC2(i,j,k)   =  q_inc_PC2(i,j,k) -  ccw_cca(i,j,k)  
                    qcl_inc_PC2(i,j,k) =  qcl_inc_PC2(i,j,k)          &
                                        + ccw_cca(i,j,k)
                    cfl_inc_PC2(i,j,k) =  cfl_inc_PC2(i,j,k)          &
                                         + cca_3d(i,j,k)
                    bcf_inc_PC2(i,j,k) =  bcf_inc_PC2(i,j,k)          &
                                    + cca_3d(i,j,k) *(1.0-cff_latest(i,j,k))

                  End if

                  ! Reset the diagnostic convective cloud
                  ! back to zero
                  cca(i,j,k) = 0.0
                  cca_2d(i,j) = 0.0
                  ccw(i,j,k) = 0.0

                End Do  ! i loop
              End Do  ! j
            End Do  ! k
!
            deallocate (ccw_cca)
            deallocate (cca_3d)

          Else     ! original code

            Do k = 1, wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
!
!
! If below cumulus cloud base then we simply zero the liquid water
! content instead of using the homogenous BL response. The forcings
! themselves are still applied but to q and T.
                  If ( cumulus(i,j)                                      &
                                   .AND. (k  <=  ntml(i,j)) ) Then
                    T_inc_PC2(i,j,k)   =  (-LC/CP) * qcl_earliest(i,j,k)
                    q_inc_PC2(i,j,k)   =  qcl_earliest(i,j,k)
                    qcl_inc_PC2(i,j,k) =  (-qcl_earliest(i,j,k))
                    cfl_inc_PC2(i,j,k) =  (-cfl_earliest(i,j,k))
                    bcf_inc_PC2(i,j,k) =  cff_latest(i,j,k)              &
                                            -cf_earliest(i,j,k)
                  End If
                End Do  ! i loop
              End Do  ! j
            End Do  ! k


          End If    ! test on lcv_pc2_diag_sh

          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
!
! Update working version of temperature, moisture and cloud fields with
! increments from the PC2 homogenous response.
!
                T_latest(i,j,k)   = T_earliest(i,j,k) + TL_force(i,j,k) &
     &                             + T_inc_PC2(i,j,k)
                q_latest(i,j,k)   = q_earliest(i,j,k) + qT_force(i,j,k) &
     &                             + q_inc_PC2(i,j,k)
                qcl_latest(i,j,k) = qcl_earliest(i,j,k)                 &
     &                             + qcl_inc_PC2(i,j,k)
                cfl_latest(i,j,k) = cfl_earliest(i,j,k)                 &
     &                             + cfl_inc_PC2(i,j,k)
                cf_latest(i,j,k)  =  cf_earliest(i,j,k)                 &
     &                             + bcf_inc_PC2(i,j,k)
!               qcf_latest and cff_latest are not updated
!
              End Do  ! i loop
            End Do  ! j
          End Do  ! k
!
          deallocate (QT_force)
          deallocate (TL_force)
          deallocate (zeros   )
!
! ----------------------------------------------------------------------
! Copy updated cloud fractions to the in/out variables
! ----------------------------------------------------------------------
!
          If (.not. L_area_cloud) Then

          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
! For the moment set area cloud fraction to the bulk cloud fraction
                area_cloud_fraction(i,j,k)=cf_latest(i,j,k)
              End Do
            End Do
          End Do

          Else If (L_area_cloud) Then

            If (L_ACF_Brooks) Then

! DEPENDS ON: ls_acf_brooks
              Call LS_ACF_Brooks (                                      &
     &             halo_i, halo_j, off_x, off_y                         &
     &,            row_length, rows, model_levels, wet_model_levels     &
     &,            r_theta_levels, delta_lambda, delta_phi              &
     &,            FV_cos_theta_latitude                                &
     &,            bulk_cloud_fraction, cloud_fraction_liquid           &
     &,            cloud_fraction_frozen, cumulus                       &
     &,            area_cloud_fraction )

            End If ! L_ACF_Brooks

          End If ! L_area_cloud
!
! call sub-level interpolation parameterisation of cloud area
          If (L_acf_cusack) then
!
! Determine number of sublevels for vertical gradient area cloud
! Want an odd number of sublevels per level: 3 is hardwired in do loops
            levels_per_level = 3
            large_levels = ((wet_model_levels - 2)*levels_per_level) + 2
!
! depends on: pc2_hom_arcld
            CALL pc2_hom_arcld(p_layer_centres,p_layer_boundaries,      &
     &       wet_model_levels,row_length,rows,                          &
     &       large_levels,levels_per_level,                             &
     &       area_cloud_fraction,T_latest,bulk_cloud_fraction,          &
     &       cloud_fraction_liquid,cloud_fraction_frozen,q_latest,      &
     &       qcl_latest,qcf_latest,                                     &
     &       lq_mix_bl)
          End if
!
! ----------------------------------------------------------------------
! Provide an estimate of convective cloud fraction for visibility
! ----------------------------------------------------------------------
!
          Do j = 1, rows
            Do i = 1, row_length
              If (cca_2d(i,j)  ==  0.0 .AND.                            &
     &             (conv_rain(i,j)+conv_snow(i,j))  >   0.0) then
! Convective precipitation exists but no estimate for its cloud fraction
! Set it to a constant value of 0.2 as an estimate.
                cca_2d(i,j) = 0.2
              End If
            End Do
          End Do
!
! ----------------------------------------------------------------------
!  PC2: End of cloud section
! ----------------------------------------------------------------------
!
        Else   !  L_pc2
!
! ----------------------------------------------------------------------
! Section BL.4b Call cloud scheme to convert Tl and qT to T, q and qcl
! in boundary layer, calculate bulk_cloud fields from qT and qcf
! and calculate area_cloud fields.
! ----------------------------------------------------------------------
!
! Determine number of sublevels for vertical gradient area cloud
! Want an odd number of sublevels per level: 3 is hardwired in do loops
        levels_per_level = 3
        large_levels = ((wet_model_levels - 2)*levels_per_level) + 2
!
! DEPENDS ON: ls_arcld
        CALL ls_arcld( p_layer_centres, RHCPT, p_layer_boundaries,      &
     &                 model_levels, wet_model_levels, row_length, rows,&
     &                 rhc_row_length, rhc_rows, bl_levels,             &
     &                 cloud_fraction_method,overlap_ice_liquid,        &
     &                 ice_fraction_method,ctt_weight,t_weight,         &
     &                 qsat_fixed,sub_cld,                              &
     &                 levels_per_level, large_levels,                  &
     &                 L_area_cloud,L_ACF_Cusack,L_ACF_Brooks,L_eacf,   &
     &                 halo_i, halo_j, off_x, off_y,                    &
     &                 delta_lambda, delta_phi,                         &
     &                 r_theta_levels, FV_cos_theta_latitude,           &
     &                 ntml, cumulus, lq_mix_bl, qcf_latest,            &
     &                 T_latest, q_latest, qcl_latest,                  &
     &                 area_cloud_fraction, bulk_cloud_fraction,        &
     &                 cloud_fraction_liquid, cloud_fraction_frozen,    &
     &                 error_code, me )
!
        End If   ! L_pc2

        If ( L_T_incr_bl_lsc .OR. L_Tl_incr_bl_lsc ) Then
          Do k=1, model_levels
            Do j=1, rows
              Do i=1, row_length
                T_incr_diag_bl(i,j,k) = T_incr_diag_bl(i,j,k)           &
     &            + ( T_latest(i,j,k) - T_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        End if                   ! on STASHflags

        If ( L_q_incr_bl_lsc .OR. L_qtl_incr_bl_lsc ) Then
          Do k=1, wet_model_levels
            Do j=1, rows
              Do i=1, row_length
                q_incr_diag_bl(i,j,k) = q_incr_diag_bl(i,j,k)           &
     &            + ( q_latest(i,j,k) - q_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        End if                  ! on STASHflags

        If ( L_qcl_incr_bl_lsc .OR. L_Tl_incr_bl_lsc .OR.               &
     &       L_qtl_incr_bl_lsc ) Then
          Do k=1, wet_model_levels
            Do j=1, rows
              Do i=1, row_length
                qcl_incr_diag_bl(i,j,k) = qcl_incr_diag_bl(i,j,k)       &
     &            + ( qcl_latest(i,j,k) - qcl_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        End if                  ! on STASHflag

        If ( L_qcf_incr_bl_lsc ) Then
          Do k=1, wet_model_levels
            Do j=1, rows
              Do i=1, row_length
                qcf_incr_diag_bl(i,j,k) = qcf_incr_diag_bl(i,j,k)       &
     &            + ( qcf_latest(i,j,k) - qcf_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        End if                  ! on STASHflag!

        If ( L_cf_incr_bl ) Then
          Do k=1, wet_model_levels
            Do j=1, rows
              Do i=1, row_length
                cf_incr_diag_bl(i,j,k) = cf_incr_diag_bl(i,j,k)         &
     &            + ( cf_latest(i,j,k) - cf_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        Endif

        If ( L_cfl_incr_bl ) Then
          Do k=1, wet_model_levels
            Do j=1, rows
              Do i=1, row_length
                cfl_incr_diag_bl(i,j,k) = cfl_incr_diag_bl(i,j,k)       &
     &            + ( cfl_latest(i,j,k) - cfl_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        Endif

        If ( L_cff_incr_bl ) Then
          Do k=1, wet_model_levels
            Do j=1, rows
              Do i=1, row_length
                cff_incr_diag_bl(i,j,k) = cff_incr_diag_bl(i,j,k)       &
     &            + ( cff_latest(i,j,k) - cff_earliest(i,j,k) )
              End Do ! i
            End Do ! j
          End Do ! k
        Endif
!
      End If ! on error code zero

!
! ----------------------------------------------------------------------
! Section BL.4c Combined cloud field calculation for use by visibility
!               (section 3) and cloud scheme (section 9) diagnostics
!               09208 - 09217 and 09223.
! ----------------------------------------------------------------------
!
! L_combi_cld_if1:
      If (error_code  <=  0  .AND.  L_combi_cld)  Then
!     Set the combined cloud area fractions in each gridbox.
!     Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
!
        Allocate ( combined_cloud(row_length, rows, wet_model_levels) )
        nclds = MIN(cloud_levels, wet_model_levels)
!
! ***** Code adapted from R2_SET_CLOUD_FIELD. *****
!
!       Zero cloud amounts in the upper layers (if necessary).
! Nclds_if1:
        If (wet_model_levels  >   nclds) Then
! Rad_k_do1:
          Do k = 1, wet_model_levels-nclds
            Do j = 1, rows
              Do i = 1, row_length
                combined_cloud(i, j, k) = 0.0E+00
              End Do
            End Do
          End Do  ! Rad_k_do1
        End If  ! Nclds_if1
!
! Rad_k_do2:
        Do k = wet_model_levels+1-nclds, wet_model_levels
!
          kinvert = wet_model_levels+1-k
! lcv_3d_cca_if1:
          If (lcv_3d_cca) Then
            Do j = 1, rows
              Do i = 1, row_length
                combined_cloud(i, j, k) = cca(i,j,kinvert) +            &
     &          ( (1.0E+00 - cca(i, j, kinvert)) *                      &
     &             area_cloud_fraction(i, j, kinvert) )
              End Do
            End Do
!
          Else
            Do j = 1, rows
              Do i = 1, row_length
                If ( (cct(i,j)  >=  kinvert+1) .AND.                    &
     &               (ccb(i,j)  <=  kinvert  ) ) Then
                  combined_cloud(i, j, k) = cca(i, j, 1) +              &
     &           ((1.0E+00-cca(i,j,1))*area_cloud_fraction(i,j,kinvert))
                Else
                  combined_cloud(i,j,k)=area_cloud_fraction(i,j,kinvert)
                End If
              End Do
            End Do
!
          End If  ! lcv_3d_cca_if1
!
        End Do  ! Rad_k_do2
!
      End if  ! L_combi_cld_if1
!
! NB: Combined cloud area fractions in each gridbox set up above.
!     Convention in Sect 70 (Radiation) is to invert levels, 1 at top.
!
! L_plsp_if1:
      If ( error_code <= 0  .AND.  L_plsp .AND. L_apply_diag ) Then
!
! DEPENDS ON: r2_calc_total_cloud_cover
        Call R2_calc_total_cloud_cover(                                 &
     &         row_length*rows, wet_model_levels, nclds                 &
     &       , IP_CLOUD_MIX_MAX, combined_cloud(1,1,1), work2d_1        &
     &       , row_length*rows, wet_model_levels                        &
     &        )
!
        Do j = 1, rows
          Do i = 1, row_length
            If (cca_2d(i,j)  <   1.0) Then
              plsp(i,j)=max(0., (work2d_1(i,j)-cca_2d(i,j))/            &
     &                                         (1.-cca_2d(i,j)))
            Else
              plsp(i,j) = 0.0
            End If
          End Do
        End Do
!
      End If  ! L_plsp_if1:
      If ( L_apply_diag ) Then

! ----------------------------------------------------------------------
! Section BL.4d Cloud scheme (section 9) diagnostics.
! ----------------------------------------------------------------------
!
! Check that cloud diagnostics are requested this timestep
! DiagSect09_if1:
      If (error_code  ==  0 .and. sf(0,9)) Then
!
! Call timer for diagnostics code
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',3)
!
! DEPENDS ON: diagnostics_lscld
        Call Diagnostics_Lscld(                                         &
     &                       row_length, rows, model_levels             &
     &,                      rhc_row_length, rhc_rows                   &
     &,                      wet_model_levels, bl_levels                &
     &,                      cloud_levels                               &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity, p_layer_centres(1, 1, 1)     &
     &,                      p, r_theta_levels, r_rho_levels            &
     &,                      T_incr_diag_bl, q_incr_diag_bl             &
     &,                      qcl_incr_diag_bl                           &
     &,                      T_latest, q_latest, qcl_latest, qcf_latest &
     &,                      area_cloud_fraction, bulk_cloud_fraction   &
     &,                      p_star, rhcpt                              &
     &,                      combined_cloud, cca, ccb, cct              &
     &,                      n_cca_levels, L_murk, Aerosol, RHcrit      &
     &,                      lq_mix_bl                                  &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &                       STASHwork9                                 &
     &                        )
!
! Call timer for diagnostics code
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',4)
      End If  ! DiagSect09_if1
!
! ----------------------------------------------------------------------
! Section BL.5 Energy correction
! ----------------------------------------------------------------------

! DEPENDS ON: timer
      If (Ltimer) Call timer ('Eng Corr',3)
      If (L_emcorr .and. Error_code  ==  0) Then

! Add surface sensible heat flux into diabatic heating
! for use in energy correction procedure.

! DEPENDS ON: flux_diag
        Call flux_diag(ftl, FV_cos_theta_latitude,                      &
     &                 row_length, rows ,off_x,off_y,                   &
     &                 1.0, sum_eng_fluxes,timestep)
! moisture flux level 1 held in fqt
! Should be total moisture flux from surface to layer 1 ie evaporation
! DEPENDS ON: flux_diag
        Call flux_diag(fqt, FV_cos_theta_latitude,                      &
     &                 row_length, rows ,off_x,off_y,                   &
     &                 1.0, sum_moist_flux,timestep)

      End If   ! L_emcorr
! DEPENDS ON: timer
      If (Ltimer) Call timer ('Eng Corr',4)

! ----------------------------------------------------------------------
! Section BL.6 Output Diagnostics requested.
! -----------------__---------------------------------------------------
      ! Take rho values on the lowest level for visibility diagnostics
      Do j = 1, rows
        Do i = 1, row_length
          rho1(i,j)=rho(i,j,1)/(r_rho_levels(i,j,1)*r_rho_levels(i,j,1))
        End Do
      End Do

      ! Calculate the total global emissions put in CO2 array - kdcorbin, 05/10
      !DEPENDS ON: tracer_fluxemit
      call tracer_fluxemit(row_length, rows, 1, timestep,  &
             co2_flux_tot,off_x,off_y, &
             FV_cos_theta_latitude, delta_lambda, delta_phi, co2emitmass)

! diagnostics requested this timestep
        If ( sf(0,3) .and. last_substep ) Then
      If ( error_code  ==  0) Then
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',3)
! DEPENDS ON: diagnostics_bl

        Call diagnostics_bl(                                            &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels, bl_levels                &
     &,                      land_points, dsm_levels                    &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      DIM_CS1, DIM_CS2                           &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      sq_T1p5, L_murk, L_pc2                     &
     &,                      Aerosol(1:row_length,1:rows,1), RHcrit     &
     &,                      L_plsp, plsp, cca_2d, lq_mix_bl            &
     &,                      ls_rain, ls_snow, conv_rain, conv_snow     &
     &,                      at_extremity                               &
     &,                      timestep                                   &
     &,                      p_star                                     &
     &,                      rhcpt, rhc_row_length, rhc_rows            &
     &,                      cloud_fraction_method,overlap_ice_liquid   &
     &,                      ice_fraction_method,ctt_weight,t_weight    &
     &,                      qsat_fixed,sub_cld                         &
     &,                      x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic   &
     &,                      lsp_ei,lsp_fi,lsp_eic,lsp_fic              &
     &,                      ntml, cumulus, L_eacf, u, v, R_u, R_v      &
     &,                      T_latest, rho1, q_latest, qcl_latest       &
     &,                      qcf_latest, cf_latest, cfl_latest          &
     &,                      cff_latest                                 &
     &,                      cf_earliest, cfl_earliest, cff_earliest    &
     &,                      T_incr_diag_bl, q_incr_diag_bl             &
     &,                      qcl_incr_diag_bl, qcf_incr_diag_bl         &
     &,                      cf_incr_diag_bl                            &
     &,                      cfl_incr_diag_bl, cff_incr_diag_bl         &
     &,                      u_incr_diag_bl,v_incr_diag_bl              &
     &,                      exner_theta_levels                         &
     &,                      t1p5m, ftl, zh                             &
     &,                      u10m, v10m, q1p5m                          &
     &,                      e_sea, h_sea,ei                            &
     &,                      sea_ice_htf, sice_mlt_htf                  &
     &,                      snomlt_surf_htf                            &
     &,                      zht                                        &
     &,                      bl_type_1,bl_type_2,bl_type_3,bl_type_4    &
     &,                      bl_type_5,bl_type_6,bl_type_7              &
     &,                      fqt, ftl, z0m_gb, z0m_eff_gb, z0h_eff_gb   &
     &,                      rib_gb, latent_heat, taux, tauy, fme       &
     &,                      t_soil, land_index                         &
     &,                      surf_ht_flux_gb                            &
     &,                      surf_ht_flux_land,surf_ht_flux_sice        &
     &,                      rib_ssi,ftl_ssi,e_ssi,ei_sice              &
     &,                      vshr_land,vshr_ssi                         &
     &,                      taux_land,taux_ssi,tauy_land,tauy_ssi      &
     &,                      radnet_sice,flandg                         &
     &,                      ntiles,npft,land_sea_mask,nice             &
     &,                      sil_orog_land,ho2r2_orog,gs,gpp,npp,resp_p &
     &,                      ecan_tile,esoil_tile,gpp_ft,ftl_tile       &
     &,                      npp_ft,resp_p_ft,resp_s,resp_s_tot         &
     &,                      resp_s_tile,cs             & !kdcorbin, 10/10
     &,                      rib_tile,es,ecan,fsmc,radnet_tile          &
     &,                      tstar_tile,canopy,catch,z0m_tile,g_leaf    &
     &,                      t1p5m_tile,q1p5m_tile,le_tile,ei_tile,olr  &
     &,                      epot_tile,tile_frac                        &
     &,                      co2_flux_tot, land_co2                     &
     &,                      l_co2_interactive                          &
     &, L_DUST,DUST_FLUX,U_S_T_TILE,U_S_T_DRY_TILE,U_S_STD_TILE,DRYDEP2 &
! CABLE
     &, TSOIL_TILE,SMCL_TILE,STHF_TILE,SNOW_DEPTH3L,SNOW_MASS3L         &
     &, SNOW_TMP3L,SNOW_RHO3L,SNOW_RHO1L,SNAGE_TILE                     &
! variables required for soil moisture nudging scheme macro
     &,                      rhokh,resfs,chr1p5m,alpha1,ra,wt_ext       &
     &,                      lai_ft,canht_ft,gc                         &
! MGS extra bl vars for UKCA
     &, rhokh_mix, rho_aresist, aresist, resist_b, r_b_dust             &
     &, dtrdz_charney_grid, kent, we_lim, t_frac, zrzi, kent_dsc        &
     &, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork3                                                       &
     & ,model_domain                                                    &
     & ,sin_theta_longitude, cos_theta_longitude                        &
     & ,proc_row_group                                                  &
     & ,BL_diag,l_cable,surf_radflux, CPOOL_TILE, NPOOL_TILE            &
     & ,PPOOL_TILE, SOIL_ORDER, TRANSP_TILE, GLAI, PHENPHASE)

! Deallocate BL_diags
      If (BL_diag%L_oblen) Deallocate(BL_diag%oblen)
      If (BL_diag%L_ustar) Deallocate(BL_diag%ustar)
      If (BL_diag%L_wbsurf) Deallocate(BL_diag%wbsurf)
      If (BL_diag%L_gradrich) Deallocate(BL_diag%gradrich)
      If (BL_diag%L_wstar) Deallocate(BL_diag%wstar)
      If (BL_diag%L_dbdz) Deallocate(BL_diag%dbdz)
      If (BL_diag%L_dvdzm) Deallocate(BL_diag%dvdzm)
      If (BL_diag%L_rhokm) Deallocate(BL_diag%rhokm)
      If (BL_diag%L_rhokh) Deallocate(BL_diag%rhokh)
      If (BL_diag%L_tke) Deallocate(BL_diag%tke)
      If (BL_diag%L_ostressx) Deallocate(BL_diag%ostressx)
      If (BL_diag%L_ostressy) Deallocate(BL_diag%ostressy)
      If (BL_diag%L_dscbase) Deallocate(BL_diag%dscbase)
      If (BL_diag%L_cldbase) Deallocate(BL_diag%cldbase)
      If (BL_diag%L_weparm) Deallocate(BL_diag%weparm)
      If (BL_diag%L_weparm_dsc) Deallocate(BL_diag%weparm_dsc)

! Deallocate BL_diags for unity arrays
      If (.NOT. BL_diag%L_oblen) Deallocate(BL_diag%oblen)
      If (.NOT. BL_diag%L_ustar) Deallocate(BL_diag%ustar)
      If (.NOT. BL_diag%L_wbsurf) Deallocate(BL_diag%wbsurf)
      If (.NOT. BL_diag%L_gradrich) Deallocate(BL_diag%gradrich)
      If (.NOT. BL_diag%L_wstar) Deallocate(BL_diag%wstar)
      If (.NOT. BL_diag%L_dbdz) Deallocate(BL_diag%dbdz)
      If (.NOT. BL_diag%L_dvdzm) Deallocate(BL_diag%dvdzm)
      If (.NOT. BL_diag%L_rhokm) Deallocate(BL_diag%rhokm)
      If (.NOT. BL_diag%L_rhokh) Deallocate(BL_diag%rhokh)
      If (.NOT. BL_diag%L_tke) Deallocate(BL_diag%tke)
      If (.NOT. BL_diag%L_ostressx) Deallocate(BL_diag%ostressx)
      If (.NOT. BL_diag%L_ostressy) Deallocate(BL_diag%ostressy)
      If (.NOT. BL_diag%L_dscbase) Deallocate(BL_diag%dscbase)
      If (.NOT. BL_diag%L_cldbase) Deallocate(BL_diag%cldbase)
      If (.NOT. BL_diag%L_weparm) Deallocate(BL_diag%weparm)
      If (.NOT. BL_diag%L_weparm_dsc) Deallocate(BL_diag%weparm_dsc)

! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',4)
      End If ! on error code zero
        Endif            ! on sf(0,3)

!
! Clear up allocatable arrays
      If (error_code  <=  0  .AND.  L_combi_cld)  Then
        Deallocate ( combined_cloud )
      End If
!
      If (L_T_incr_bl_lsc .OR. L_Tl_incr_bl_lsc .OR. L_PC2) Then
        Deallocate (T_earliest)
      End If
!
      If(L_q_incr_bl_lsc .OR. L_qtl_incr_bl_lsc                         &
     &   .OR. L_qcl_incr_bl_lsc .OR. L_Tl_incr_bl_lsc                   &
     &   .OR. L_qcf_incr_bl_lsc .OR. L_PC2) Then
        Deallocate (q_earliest)
        Deallocate (qcl_earliest)
        Deallocate (qcf_earliest)
        Deallocate (cf_earliest)
        Deallocate (cfl_earliest)
        Deallocate (cff_earliest)
      End If
!
      If (L_PC2) Then
        Deallocate ( T_inc_PC2 )
        Deallocate ( q_inc_PC2 )
        Deallocate ( qcl_inc_PC2 )
        Deallocate ( cfl_inc_PC2 )
        Deallocate ( bcf_inc_PC2 )
      End If

      End If ! L_apply_diag
! end of routine NI_imp_ctl
      Return
      END SUBROUTINE NI_imp_ctl

