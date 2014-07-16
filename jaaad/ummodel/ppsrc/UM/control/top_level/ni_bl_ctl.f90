

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine NI_bl_ctl
! *********************************************************************


      Subroutine NI_bl_ctl (                                            &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &

! IN Substep number for ATMPHYS2
     &,Substep_Number, Num_Substeps, L_phys2_substep, NumCycles,CycleNo &

! model dimensions.
     &, row_length, rows, n_rows, land_points, ntiles, model_levels     &
     &, wet_model_levels, bl_levels, DIM_CS1, DIM_CS2                   &
     &, dst_levels, dsm_levels                                          &

! Model switches
     &, model_domain, formdrag, orog_drag_param, L_CAL360               &
     &, L_area_Cloud, L_mixing_ratio, L_emcorr, Ltimer, L_dry           &
     &, L_scrn, L_plsp, BL_OPTIONS, L_SBLeq, L_SBLco, L_ctile           &
     &, L_DUST, L_CAM_DUST, L_ukca, L_LAMBDAM2, L_FULL_LAMBDAS          &

! model Parameters
     &, alpha_cd, CO2_MMR, MUW_SBL, MWT_SBL                             &
     &, Charnock, SeaSalinityFactor                                     &

! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &

! in coordinate information
     &, r_rho_levels, r_theta_levels, r_at_u, r_at_v                    &
     &, eta_theta_levels, eta_rho_levels, delta_lambda, delta_phi       &
     &, lat_rot_NP, long_rot_NP                                         &

! in time stepping information.
     &, timestep, val_year, val_day_number, val_hour, val_minute        &
     &, val_second, timestep_number                                     &

! trig arrays
     &, sin_theta_longitude, cos_theta_longitude, FV_cos_theta_latitude &
     &, sin_theta_latitude                                               &

! in data fields.
     &, p, p_layer_centres, p_layer_boundaries, rho_rsq, rho_wet, rho_dry&
     &, u, v, w, etadot                                                 &
     &, u_conv, v_conv, land_sea_mask, q, qcl, qcf, p_star, theta       &
     &, EXNER_THETA_LEVELS, RAD_HR, MICRO_TENDS, SOIL_LAYER_MOISTURE    &
! variables for subgrid turbulence scheme
     &, visc_BL_m, FM_3D, FH_3D, L_subfilter_vert, L_subfilter_horiz    &
     &, L_subfilter_blend, max_diff,turb_startlev_vert,turb_endlev_vert &
     &, BL_COEF_KM, BL_COEF_KH                                          &

! ancillary fields and fields needed to be kept from timestep to
! timestep
     &, hcon, smvccl, smvcwt, smvcst, sthf, sthu, sil_orog_land         &
     &, ho2r2_orog, ice_fract, u_0, v_0, u_0_p, v_0_p                   &
     &, land_index, cca, ccb, cct, photosynth_act_rad, surf_radflux     &

! IN variables required for mineral dust scheme
     &, SOIL_CLAY,SOIL_SILT,SOIL_SAND,DUST_MREL1,DUST_MREL2,DUST_MREL3  &
     &, DUST_MREL4,DUST_MREL5,DUST_MREL6                                &

! IN additional variables for MOSES II
     &,LAND_PTS_TRIF,NPFT_TRIF,CANOPY,CATCH,CATCH_SNOW                  &
     &,SNOW_TILE,Z0_TILE,LW_DOWN,SW_TILE,TSTAR_TILE                     &
     &,CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE                 &
     &,L_PHENOL,L_TRIFFID,L_Q10,ASTEPS_SINCE_TRIFFID,CAN_MODEL          &
     &,CS,FRAC,CANHT_FT,LAI_FT,FLAND,FLANDG,TSTAR_SEA,ANTHROP_HEAT      &
     &,ALBSOIL,COS_ZENITH_ANGLE,CAN_RAD_MOD,ILAYERS                     &

!    EAK
!    IN
     &, surf_down_sw,alb_tile,l_tile_pts                                &
!     &, surf_down_sw,alb_tile,cos_zenith_angle               &
     &, ls_rain,ls_snow,SW_DOWN                                         &
     &, lat,long,day,time_sec                                           &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L                   &
     &, SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE             &
     &, TSOIL_TILE,T_SURF_TILE,HCONS,BEXP                               &
     &, SATHH,SATCON,HCAP                                               &
     &, SOIL_TYPE,VEG_TYPE                                              &
     &, ISNOW_FLG3L,total_nsteps                                        &
     &           ,FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB               &
     &           ,TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB           &
     &           ,USTAR_CAB,SURF_HTF_CAB                                &
     &, l_cable                                                         &
!
     &           ,U_S_CAB,CH_CAB,CD_CAB                                 &
     &           ,SNAGE_TILE,RTSOIL_TILE                                &
     &           ,GFLUX_TILE,SGFLUX_TILE                                &
     &, CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER                     &
     &, NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE                           &

! in/out
     &, t_soil, ti, t_surf, z0msea                                      &
! in SCM namelist data
     &, L_spec_z0, z0m_scm, z0h_scm                                     &

! in/out
     &, area_cloud_fraction, bulk_cloud_fraction                        &
     &, T_latest, q_latest, qcl_latest, qcf_latest                      &
     &, R_u, R_v, cloud_fraction_liquid, cloud_fraction_frozen          &
     &, zh, zh_prev, flux_e, flux_h, L_flux_bc                          &

! INOUT additional variables for MOSES II
     &,GS,G_LEAF_ACC,NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC                &
     &, TSTAR_LAND,TSTAR_SICE,TSTAR_SSI                                 &

! INOUT Variables for STPH_RP
     &,G0_RP,par_mezcla                                                 &
!
! diagnostic info
! SCM Diagnostics (dummy values in full UM)
     &, nSCMDpkgs,L_SCMDiags,                                           &
!
! STASH Diagnostics
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end

! out variables required in IMP_SOLVER
     & alpha1_sice, ashtf, dtrdz_charney_grid, rdz_charney_grid         &
     &, dtrdz_u, dtrdz_v, rdz_u, rdz_v, cdr10m_u, cdr10m_v, z1_tq       &

! OUT variables which need to be maintained between substeps
     &, RHO_UV,RHO_TQ,DZL_CHARNEY,RDZ                                   &
     &, Z1_UV,Z_FULL,Z_UV,Z_TQ,P_HALF,DELTAP                            &

! out diagnostics (done after implicit solver)
     &, cd, ch, e_sea, fqT, ftl, h_sea, rib_gb                          &
     &, taux, tauy, vshr, zht, shallowc, cu_over_orog                   &
     &, bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6     &
     &, bl_type_7, z0m_eff_gb, z0h_eff_gb, fme                          &


! OUT diagnostics required for soil moisture nudging scheme :
     &, WT_EXT,RA                                                       &

! out data required for tracer mixing :
     &, rho_aresist,aresist,resist_b                                    &

!OUT variables required for mineral dust scheme
     &, R_B_DUST,DUST_FLUX,U_S_T_TILE,U_S_T_DRY_TILE,U_S_STD_TILE       &
     &, KENT, WE_LIM, T_FRAC, ZRZI                                      &
     &, KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC, ZHSC, Z_HALF        &

! OUT data required for 4D-VAR :
     &, rho_cd_modv1, rho_km_var                                        &

! OUT additional variables for MOSES II
     &,FTL_TILE,LE_TILE,RADNET_SICE                                     &
     &,RADNET_TILE,RIB_TILE,RHO_ARESIST_TILE,ARESIST_TILE               &
     &,RESIST_B_TILE,ALPHA1,ASHTF_TILE,FQT_TILE,EPOT_TILE,FQT_ICE       &
     &,FTL_ICE,FRACA,RESFS,RESFT,RHOKH_TILE,RHOKH_SICE                  &
     &,RHOKPM,RHOKPM_POT                                                &
     &,RHOKPM_SICE,Z0HSSI,Z0H_TILE,Z0M_GB,Z0MSSI,Z0M_TILE,CHR1P5M       &
     &,CHR1P5M_SICE,SMC,GPP,NPP,RESP_P,G_LEAF,GPP_FT,NPP_FT             &
     &,RESP_P_FT,RESP_S,RESP_S_TOT,RESP_W_FT                            &
     &,GC,CANHC_TILE,WT_EXT_TILE,FLAKE,TILE_INDEX                       &
     &,TILE_PTS,TILE_FRAC,FSMC,FLANDG_U,FLANDG_V,rib_ssi                &
     &,taux_land,taux_ssi,tauy_land,tauy_ssi,vshr_land,vshr_ssi         &

! out fields
     &, t1_sd, q1_sd, ntml, cumulus, nbdsc, ntdsc                       &
     &, ntpar, nlcl, zhpar, zlcl, l_shallow, wstar, wthvs, delthvu      &
     &, uw0,vw0, rhokm,rhokm_u,rhokm_v, rhokh                           &

! error information
     &, Error_code,BL_diag,                                             &
      ! end step of experiment, processor num
     & endstep, mype )

! purpose: Interface to boundary layer mixing coefficients calculation
!
!          IN/OUT etc intents to be added later.
!
! method:
!
!   language: fortran 90 + cray extensions
!   this code is written to umdp3 programming standards.

      Use bl_diags_mod, Only :                                          &
          strnewbldiag

      Implicit None

      ! end step of experiment, processor num
      integer :: endstep, mype

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
     &, me                                                              &
                   ! My processor number
     &,NumCycles                                                        &
                  ! Number of phys-dyn iterations per tstep
     &,CycleNo    ! Iteration no

! Substep number for ATMPHYS2
      Integer                                                           &
     & Substep_Number                                                   &
     &,Num_Substeps

! Switch for calculating exchange coeffs from latest values.
      Logical                                                           &
     & L_phys2_substep
      LOGICAL                                                           &
     & L_mixing_ratio         ! TRUE if mixing ratios used in
!                             ! boundary layer code

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north, sout
                         ! east or west of the processor grid
     &, L_DUST                                                          &
               !switch for mineral dust
     &, L_CAM_DUST                                                      &
               !Old version of dust_uplift scheme used in CAM NWP models
     &, L_ukca !   switch for UKCA scheme

! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, land_points                                                     &
                    ! IN No.of land points being processed, can be 0.
     &, ntiles                                                          &
                    ! IN No. of land-surface tiles ( MOSES II )
     &, model_levels                                                    &
     &, wet_model_levels                                                &
     &, bl_levels                                                       &
     &, DIM_CS1, DIM_CS2                                                &
                              ! soil carbon dimensions
     &, dst_levels                                                      &
                    ! number of deep soil temperature levels
     &, dsm_levels  ! number of deep soil moisture levels

! Model switches
      Integer                                                           &
     &  model_domain                                                    &
     &, formdrag   ! No_drag if no orographic form drag
!                  ! Effective_z0 if effective roughness length used
!                  ! Explicit_stress if distributed form drag used

      Logical                                                           &
     &  L_SBLeq                                                         &
                    ! IN Switch for Equilibrium SBL model

     &, L_SBLco     ! IN Switch for coupled gradient method in
!                   !    Equilibrium SBL model
!
! Boundary Layer
      REAL, Intent(IN) :: SeaSalinityFactor
!                   ! Factor allowing for the effect of
!                   ! the salinity of sea water on the
!                   ! evaporative flux.
!
      REAL, Intent(IN) :: orog_drag_param
!                   ! Drag coefficient for orographic form drag

      Logical                                                           &
     &  L_CAL360                                                        &
                    ! true if using 360 day calender
     &, L_emcorr                                                        &
                    ! true if energy correction scheme is to be used.
     &, L_dry                                                           &
                    ! true if model to be run with no moisture
     &, L_ctile     ! true if model to be run with coastal tiling

      Logical                                                           &
     &  L_area_Cloud                                                    &

     &, Ltimer                                                          &
                 ! true then output some timing information
     &, L_scrn                                                          &
                                 ! Logical to control output
                                 !    of screen level T,Q,QCL,QCF
     &, L_plsp                                                          &
                                 ! Logical to control output
                                 !    of Probability of LS Precip
     &,L_LAMBDAM2                                                       &
                   ! IN LambdaM=2*LambdaH (operational setting).
     &,L_FULL_LAMBDAS ! Lambdas NOT reduced above NTML_LOCAL+1

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
     &, CO2_MMR         ! set equal to co2_start

      Real                                                              &
     &  alpha_cd(bl_levels)                                             &
     &, Muw_SBL,Mwt_SBL                                                 &
                        ! IN Powers to use in prescription of equilib
!                       !    profiles of stress and buoyancy flux in
!                       !    Equilibrium SBL model
     &, Charnock   ! Charnock parameter for sea surface

      INTEGER, DIMENSION(20) :: BL_OPTIONS   ! IN BL switches

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
     &, sin_theta_latitude (row_length, rows)                          &
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

! Diagnostic variables
      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP
!
! Additional variables for SCM diagnostics which are dummy in full UM
      Integer                                                           &
     & nSCMDpkgs                ! No of SCM diagnostics packages
!
      Logical                                                           &
     & L_SCMDiags(nSCMDpkgs)    ! Logicals for SCM diagnostics packages

! Data arrays
      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &,  u_conv(1-off_x:row_length+off_x,1-off_y:rows+off_y,            &
     &              bl_levels)                                          &
     &,  v_conv(1-off_x:row_length+off_x,1-off_y:n_rows+off_y,          &
     &              bl_levels)                                          &
     &, w(row_length, rows,0:model_levels)                              &
     &, etadot(row_length, rows,0:model_levels)                         &
     &, rho_rsq(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         model_levels)                                            &
!                       ! wet density times r^2 on rho levels (kg/m3)
     &, rho_wet(row_length, rows, model_levels)                         &
!                       ! wet density on rho levels (kg/m3)
     &, rho_dry(row_length, rows, model_levels)                         &
!                       ! dry density on rho levels (kg/m3)
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
     &, rad_hr(row_length,rows,bl_levels,2)                             &
                              ! IN (LW,SW) radiative heating rate (K/s)
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
     &, BL_COEF_KH(1:row_length, 1:rows, bl_levels-1)                   &
!            ! RHOKH from BL scheme
     &, micro_tends(row_length, rows, bl_levels, 2)
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)

!     Declaration of new BL diagnostics.
      Type (Strnewbldiag) :: BL_diag
      logical                                                           &
     &  L_subfilter_vert                                                &
                            ! subgrid turbulence scheme in vertical     
     &, L_subfilter_horiz                                               &
                            ! subgrid turbulence scheme in horizontal   
     &, L_subfilter_blend                                           
                            ! Blending of BL and Smag coefficients

      Integer                                                           &
     & turb_startlev_vert                                               &
                            ! start and end vertical levels for         
     &,turb_endlev_vert     
                            ! 3D turbulence scheme
     

       REAL                                                             &
     & SOIL_LAYER_MOISTURE(LAND_POINTS,DSM_LEVELS)!IN soil moisture
!                 ! per layer (kg m-2)
      logical                                                           &
     &  land_sea_mask(row_length, rows)

! ancillary arrays and fields required to be saved from timestep to
! timestep.
      Integer                                                           &
     &  land_index (land_points)      ! set from land_sea_mask

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
     &, smvccl (land_points)                                            &
                             ! soil/qrparm.soil.crit
     &, smvcwt (land_points)                                            &
                             ! soil/qrparm.soil.wilt
     &, smvcst (land_points)                                            &
                             ! soil/qrparm.soil.satn
     &, sthf(land_points,dsm_levels)                                    &
                                ! IN Frozen soil moisture content of
                                !    each layer as a fraction of
                                !    saturation.
     &, sthu(land_points,dsm_levels)
                                ! IN Unfrozen soil moisture content
                                !    of each layer as a fraction of
                                !    saturation.

      Real                                                              &
     &  sil_orog_land (land_points)                                     &
                                    ! orog/qrparm.orog.as
     &, ho2r2_orog (land_points) ! orog/qrparm.orog.h2root2

      Real                                                              &
     &  ice_fract (row_length, rows) ! ice/qrclim.ice.(month)

      Real                                                              &
     &  photosynth_act_rad(row_length, rows)                            &
                                             ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
     &, surf_radflux(row_length, rows) ! Net downward
      REAL                                                              &
           ! mineral dust fields
     &  SOIL_CLAY ( ROW_LENGTH, ROWS )                                  &
     &, SOIL_SILT ( ROW_LENGTH, ROWS )                                  &
     &, SOIL_SAND ( ROW_LENGTH, ROWS )                                  &
     &, DUST_MREL1 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL2 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL3 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL4 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL5 ( ROW_LENGTH, ROWS )                                 &
     &, DUST_MREL6 ( ROW_LENGTH, ROWS )

      Real                                                              &
     &  cca (row_length, rows)

      Integer                                                           &
     &  ccb (row_length, rows)                                          &
     &, cct (row_length, rows)


! arguments with intent in/out. ie: input variables changed on output.
      Real                                                              &
     &  T_surf(row_length, rows)                                        &
     &, z0msea(row_length, rows)                                        &
                                 ! Sea surface roughness length(m)
                                 ! for momentum
     &, z0m_scm(row_length, rows)                                       &
                                 ! Fixed Sea surface roughness
                                 ! length(m) for momentum (SCM)
     &, z0h_scm(row_length, rows)                                       &
                                 ! Fixed Sea surface roughness
                                 ! length(m) for heat (SCM)
     &, t_soil(land_points,dsm_levels)                                  &
                                       ! slt/qrclim.slt_pm(lev).(month)
     &, ti(row_length, rows)                                            &
                              ! set equal to tstar
     &, R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        model_levels)                                             &
     &, T_latest(row_length, rows, model_levels)                        &
     &, q_latest(row_length, rows, wet_model_levels)                    &
     &, qcl_latest(row_length, rows, wet_model_levels)                  &
     &, qcf_latest(row_length, rows, wet_model_levels)                  &
     &, area_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, bulk_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, cloud_fraction_liquid(row_length, rows, wet_model_levels)       &
     &, cloud_fraction_frozen(row_length, rows, wet_model_levels)       &
     &, zh(row_length, rows)                                            &
     &, zh_prev(row_length, rows)                                       &
                                  ! IN boundary layer height from
!                                 !    previous timestep
     &, flux_e(row_length,rows)                                         &
     &, flux_h(row_length,rows)
      LOGICAL                                                           &
     & L_flux_bc                                                        &
                    ! T if prescribed surface fluxes to be used
     &,L_spec_z0    ! T is roughness lengths have been specified

      Logical                                                           &
     &  cumulus (row_length, rows)                                      &
                                   ! *APL bl convection flag
     &, l_shallow(row_length, rows)      ! Logical indicator of
!                                        ! shallow convection
      Real                                                              &
     &  delthvu(row_length, rows)                                       &
                                     ! buoyancy integral
     &, zhpar(row_length, rows)                                         &
                                     ! height of ntpar
     &, zlcl(row_length, rows)                                          &
                                     ! height of lcl accurate value not 
                                     ! a model level height (m)
     &, wstar(row_length, rows)                                         &
                                     ! surface-based mixed layer
!                                    ! velocity scale
     &, wthvs(row_length, rows)      ! surface buoyancy flux

      Real                                                              &
     &  uw0(row_length, rows)                                           &
!                       ! U-component of surface wind stress (P-grid)
     &, vw0(row_length, rows)
!                       ! V-component of surface wind stress (P-grid)

! Definition of variables for STPH_RP
      REAL,INTENT(InOut) :: par_mezcla ! Used to modify
                                       ! LAMBDAH,LAMBDAM in the
                                       ! Boundary Layer 8A
      REAL,INTENT(InOut) :: G0_RP ! Flux profile param. in BL8A
! arguments with intent out. ie: output variables.
      Real                                                              &
     &  t1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, q1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, z0m_eff_gb(row_length,rows)                                     &
     &, z0h_eff_gb(row_length,rows)
                                ! Effective grid-box roughness
!                                 lengths for momentum and for 
!                                 heat, moisture

      Integer                                                           &
     &  ntml (row_length, rows)                                         &
                                    ! IN/OUT
     &, ntpar(row_length, rows)                                         &
                                    ! IN/OUT top of diagnostic parcel
!                                   !        ascent
     &, nlcl(row_length, rows)                                          &
                                    ! IN/OUT lifting condensation level
     &, nbdsc(row_length, rows)                                         &
     &, ntdsc(row_length, rows)                                         &
     &, KENT(row_length,rows)                                           &
                                    ! OUT grid-level of SML inversion
     &, KENT_DSC(row_length,rows)   ! OUT grid-level of DSC inversion

      Real                                                              &
     &  rhokm(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &           0:bl_levels-1)                                         &
     &  ,rhokm_u(row_length,rows,bl_levels)                             &
     &  ,rhokm_v(row_length,n_rows,bl_levels)

! variables passed from BDY_LAYR to IMP_SOLVER
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
! OUT variables which need to be maintained between substeps. They will
!     be used by BDY_LAYR()
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
     &,Z_FULL(row_length,rows,BL_LEVELS)                                &
                                ! Z_FULL(*,K) is height of full level k.
     &,Z_UV(row_length,rows,BL_LEVELS)                                  &
                                ! Z_UV(*,K) is height of half level
!                               ! k-1/2.
     &,Z_TQ(row_length,rows,BL_LEVELS)                                  &
                                ! Z_TQ(*,K) is height of half level
!                               ! k+1/2.
     &,P_HALF(row_length,rows,BL_LEVELS)                                &
                                ! P_HALF(*,K) is pressure at half
!                               ! level k-1/2.
     &,DELTAP(row_length,rows,BL_LEVELS)
                                ! Difference in pressure between levels

! Diagnostics needed in NI_imp_ctl

      Real                                                              &
     &  fme(row_length, rows)                                           &
     &, rho_aresist(row_length,rows)                                    &
                                ! OUT RHOSTAR*CD_STD*VSHR for SULPHUR
                                !     cycle
     &, aresist(row_length,rows)                                        &
                                ! OUT 1/(CD_STD*VSHR) for Sulphur cycle
     &, resist_b(row_length,rows)                                       &
                                ! OUT (1/CH-1/(CD_STD)/VSHR for
                                !     Sulphur cycle
     &, R_B_DUST(ROW_LENGTH,ROWS,ndiv)                                  &
                                      !OUT surface layer resist for dust
     &, DUST_FLUX(ROW_LENGTH,ROWS,NDIV)                                 &
                                       !OUT dust emissions (kg m-2 s-1)
     &, U_S_T_TILE(LAND_POINTS,NTILES,NDIV)                             &
                                           !OUT threshold frict. vel
     &, U_S_T_DRY_TILE(LAND_POINTS,NTILES,NDIV)                         &
                                               !OUT dry soil value
     &, U_S_STD_TILE(LAND_POINTS,NTILES)                                &
                                        !OUT friction velocity
     &, WE_LIM(row_length,rows,3)                                       &
                                    ! OUT rho*entrainment rate implied b
!                                   !     placing of subsidence
     &, ZRZI(row_length,rows,3)                                         &
                                    ! OUT (z-z_base)/(z_i-z_base)
     &, T_FRAC(row_length,rows,3)                                       &
                                    ! OUT a fraction of the timestep
     &, WE_LIM_DSC(row_length,rows,3)                                   &
!                                   ! OUT rho*entrainment rate implied b
!                                   !     placing of subsidence
     &, ZRZI_DSC(row_length,rows,3)                                     &
                                    ! OUT (z-z_base)/(z_i-z_base)
     &, T_FRAC_DSC(row_length,rows,3)                                   &
!                                   ! OUT a fraction of the timestep
     &, Z_HALF(row_length,rows,BL_LEVELS)                               &
                                    ! OUT Z_HALF(*,K) is height of half
!                                   ! level k-1/2.
     &, ZHSC(row_length,rows)                                           &
                                    ! OUT Top of decoupled layer
     &, zht(row_length, rows)                                           &
                                      ! Max height of turb mixing
     &, shallowc(row_length,rows)                                       &
                                   !OUT Indicator set to 1.0 if shallow,
!                                  !   0.0 if not shallow or not cumulus
     &, cu_over_orog(row_length,rows)                                   &
!                                  ! OUT Indicator for cumulus
!                                  !     over steep orography
!                                  !   Indicator set to 1.0 if true,
!                                  !   0.0 if false. Exclusive.
     &, bl_type_1(row_length,rows)                                      &
                                  ! OUT Indicator set to 1.0 if stable
!                                  !     b.l. diagnosed, 0.0 otherwise.
     &, bl_type_2(row_length,rows)                                      &
                                  ! OUT Indicator set to 1.0 if Sc over
!                                  !     stable surface layer diagnosed,
!                                  !     0.0 otherwise.
     &, bl_type_3(row_length,rows)                                      &
                                  ! OUT Indicator set to 1.0 if well
!                                  !     mixed b.l. diagnosed,
!                                  !     0.0 otherwise.
     &, bl_type_4(row_length,rows)                                      &
                                  ! OUT Indicator set to 1.0 if
!                                  !     decoupled Sc layer (not over
!                                  !     cumulus) diagnosed,
!                                  !     0.0 otherwise.
     &, bl_type_5(row_length,rows)                                      &
                                  ! OUT Indicator set to 1.0 if
!                                  !     decoupled Sc layer over cumulus
!                                  !     diagnosed, 0.0 otherwise.
     &, bl_type_6(row_length,rows)                                      &
                                  ! OUT Indicator set to 1.0 if a
!                                  !     cumulus capped b.l. diagnosed,
!                                  !     0.0 otherwise.
     &, bl_type_7(row_length,rows)! OUT Indicator set to 1.0 if a
!                                  !     shear-dominated  b.l.
!                                  !     diagnosed, 0.0 otherwise.

      Real                                                              &
             ! (output from sf_exch)
     &  rho_cd_modv1(row_length, rows)                                  &
     &, rho_km_var(row_length,rows,2:bl_levels)                         &
     &, vshr(row_length, rows)

      Real                                                              &
            ! output from bdy_layr.
     &  taux(row_length, rows, bl_levels)                               &
                                             ! needed as diagnostic
     &, tauy(row_length, n_rows, bl_levels)                             &
                                             ! needed as diagnostic
     &, ftl(row_length, rows, bl_levels)                                &
                                             ! needed as diagnostic
     &, fqT(row_length, rows, bl_levels)                                &
                                             ! needed as diagnostic ?
     &, h_sea(row_length, rows)                                         &
                                             ! needed as diagnostic ?
     &, e_sea(row_length, rows)                                         &
                                             ! needed as diagnostic ?
     &, ch(row_length, rows)                                            &
                                             ! needed as diagnostic ?
     &, cd(row_length, rows)                                            &
                                             ! needed as diagnostic ?
     &, rib_gb(row_length,rows)                                         &
                                ! Mean bulk Richardson number for
                                     !  lowest layer.
     &, rhokh (row_length, rows, bl_levels)                             &
     &, WT_EXT(LAND_POINTS,dsm_levels)                                  &
                                      !OUT cumulative fract of transp'n
     &, RA(LAND_POINTS)               !OUT Aerodynamic resistance (s/m)

      Integer                                                           &
     &  Error_code

! Additional MOSES II variables

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



      INTEGER                                                           &
     & LAND_PTS_TRIF                                                    &
                                   ! IN For dimensioning land fields
     &,NPFT_TRIF                                                        &
                                   ! IN For dimensioning PFT fields
!                                  !    available only with TRIFFID.
!                                  !    Set to NPFT when TRIFFID on,
!                                  !    set to 1 when TRIFFID off.
     &,CO2_DIM_LEN                                                      &
                                   ! IN Length of a CO2 field row.
     &,CO2_DIM_ROW                                                      &
                                   ! IN Number of CO2 field rows.
     &,ASTEPS_SINCE_TRIFFID                                             &
                                   ! IN Number of atmospheric
!                                  !    timesteps since last call
!                                  !    to TRIFFID.
     &,CAN_MODEL                   ! IN Swith for thermal vegetation
!                                  !    canopy

      REAL                                                              &
     & CANOPY(LAND_POINTS,NTILES)                                       &
                                      ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
     &,CATCH(LAND_POINTS,NTILES)                                        &
                                      ! IN Surface/canopy water capacity
!                                  !    of snow-free land tiles (kg/m2).
     &,CATCH_SNOW(LAND_POINTS)                                          &
                                   ! IN Snow interception capacity of
!                                  !    NLT tile (kg/m2).
     &,SNOW_TILE(LAND_POINTS,NTILES)                                    &
                                      ! IN Lying snow on tiles (kg/m2)
     &,Z0_TILE(LAND_POINTS,NTILES)                                      &
                                      ! IN Tile roughness lengths (m).
     &,LW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface downward LW radiation
!                                  !    (W/m2).
     &,SW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! Surface downward SW radiation (W/m2).
     &,SW_TILE(LAND_POINTS,NTILES)                                      
                                      ! IN Surface net SW radiation on
!------------------------------------------------------
!     EAK
      Real                                                          &
     &  alb_tile(land_points,ntiles,4)                              &
     &, SNOW_DEPTH3L(LAND_POINTS,NTILES,3)                          &
     &, SNOW_MASS3L(LAND_POINTS,NTILES,3)                           &
     &, SNOW_COND(LAND_POINTS,NTILES,3)                             &
     &, SNOW_TMP3L(LAND_POINTS,NTILES,3)                            &
     &, SNOW_RHO3L(LAND_POINTS,NTILES,3)                            &
     &, SNOW_RHO1L(LAND_POINTS,NTILES)                              &
     &, SNAGE_TILE(LAND_POINTS,NTILES)                              &
     &, SMCL_TILE(LAND_POINTS,NTILES,dsm_levels)                    &
     &, STHU_TILE(LAND_POINTS,NTILES,dsm_levels)                    &
     &, STHF_TILE(LAND_POINTS,NTILES,dsm_levels)                    &
     &, TSOIL_TILE(LAND_POINTS,NTILES,dsm_levels)                   &
     &, T_SURF_TILE(LAND_POINTS,NTILES)                             &
     &, RTSOIL_TILE(LAND_POINTS,NTILES)                             &
     &, GFLUX_TILE(LAND_POINTS,NTILES)                              &
     &, SGFLUX_TILE(LAND_POINTS,NTILES)                             &
     &, BEXP(LAND_POINTS)                                           &
     &, SATCON(LAND_POINTS)                                         &
     &, SATHH(LAND_POINTS)                                          &
     &, HCAP(LAND_POINTS)                                           &
     &, HCONS(LAND_POINTS)                                          &
     &, CPOOL_TILE(LAND_POINTS,NTILES,10)                           &
     &, NPOOL_TILE(LAND_POINTS,NTILES,10)                           &
     &, PPOOL_TILE(LAND_POINTS,NTILES,12)                           &
     &, SOIL_ORDER(LAND_POINTS)                                     &
     &, NIDEP(land_points)                                          &
     &, NIFIX(land_points)                                          &
     &, PWEA(land_points)                                           &
     &, PDUST(land_points)                                          &
     &, GLAI(land_points,NTILES)                                    &
     &, PHENPHASE(land_points,NTILES)                               &

     &, surf_down_sw(row_length,rows,4)                                 &
     &, ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, lat(row_length, rows)                                           &
     & ,long(row_length, rows)                                         &
     & ,time_sec                                                       

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
     &,CH_CAB(LAND_POINTS)                                              &
     &,CD_CAB(LAND_POINTS)                                              &
     &,U_S_CAB(LAND_POINTS)

      Integer                                                           &
     &  day                                                             &
     &, total_nsteps                                                    &
                                ! Total number of steps in run          
     &, SOIL_TYPE(row_length,rows)                                      &
     &, VEG_TYPE(row_length,rows)                                       &
     &, ISNOW_FLG3L(LAND_POINTS,NTILES)

      Logical l_cable
      LOGICAL, DIMENSION(LAND_POINTS,NTILES) :: L_TILE_PTS

!------------------------------------------------------

      REAL                                                              &
     & TSTAR_TILE(LAND_POINTS,NTILES)                                   &
                                      ! IN Surface tile temperatures
     &,CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)                                  &
!                                  ! IN 3D CO2 field if required.
     &,CS(LAND_POINTS,DIM_CS1)                                          &
                                 ! IN Soil carbon (kg C/m2).
     &,FRAC(LAND_POINTS,NTYPE)                                          &
                                      ! IN Fractions of surface types.
     &,CANHT_FT(LAND_POINTS,NPFT)                                       &
                                      ! IN Canopy height (m)
     &,LAI_FT(LAND_POINTS,NPFT)                                         &
                                      ! IN Leaf area index
     &,FLAND(LAND_POINTS)                                               &
                                   ! IN Land fraction on land points.
     &,FLANDG(ROW_LENGTH,ROWS)                                          &
                                   ! IN Land fraction on all points.
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! IN Open sea sfc temperature (K).
     &,ALBSOIL(LAND_POINTS)                                             &
                                   ! Soil albedo.
     &,COS_ZENITH_ANGLE(ROW_LENGTH, ROWS)
!                                  ! Cosine of the zenith angle
      INTEGER                                                           &
     &CAN_RAD_MOD                                                       &
                                   !Switch for canopy radiation model
     &,ILAYERS                     !No of layers in canopy radiation
                                   ! model


      LOGICAL                                                           &
     & L_CO2_INTERACTIVE                                                &
                                   ! IN Switch for 3D CO2 field
     &,L_PHENOL                                                         &
                                   ! IN Indicates whether phenology
!                                  !    in use
     &,L_TRIFFID                                                        &
                                   ! IN Indicates whether TRIFFID
!                                  !    in use.
     &,L_Q10                       ! IN Indicates Q10 for soil resp'n



      REAL                                                              &
     & GS(LAND_POINTS)                                                  &
                                      ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
     &,G_LEAF_ACC(LAND_POINTS,NPFT)                                     &
                                      ! INOUT Accumulated G_LEAF
! Lestevens 17may13 - change npft to ntiles
     &,NPP_FT_ACC(LAND_PTS_TRIF,NTILES)                                 &
!    &,NPP_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                              &
!                                  ! INOUT Accumulated NPP_FT
     &,RESP_W_FT_ACC(LAND_PTS_TRIF,NTILES)                              &
!    &,RESP_W_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                           &
!                                  ! INOUT Accum RESP_W_FT
     &,RESP_S_ACC(LAND_PTS_TRIF,DIM_CS2)                                &
                                         ! INOUT Accumulated RESP_S
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
                                   ! INOUT Land mean sfc temperature (K)
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! INOUT Sea-ice sfc temperature (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS)  ! INOUT Sea mean sfc temperature (K).

      INTEGER                                                           &
     & TILE_INDEX(LAND_POINTS,NTYPE)                                    &
                                      ! OUT Index of tile points
     &,TILE_PTS(NTYPE)             ! OUT Number of tile points


      REAL                                                              &
     & FTL_TILE(LAND_POINTS,NTILES)                                     &
                                      ! OUT Surface FTL for land tiles
     &,LE_TILE(LAND_POINTS,NTILES)                                      &
                                      ! OUT Surface latent heat flux for
!                                  !     land tiles
     &,RADNET_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! OUT Surface net radiation on
!                                  !     sea-ice (W/m2)
     &,RADNET_TILE(LAND_POINTS,NTILES)                                  &
                                      ! OUT Surface net radiation on
!                                  !     land tiles (W/m2)
     &,RIB_TILE(LAND_POINTS,NTILES)                                     &
                                      ! OUT RIB for land tiles.
     &,RHO_ARESIST_TILE(LAND_POINTS,NTILES)                             &
!                                  ! OUT RHOSTAR*CD_STD*VSHR on land
!                                  !     tiles
     &,ARESIST_TILE(LAND_POINTS,NTILES)                                 &
!                                  ! OUT 1/(CD_STD*VSHR) on land tiles
     &,RESIST_B_TILE(LAND_POINTS,NTILES)                                &
!                                  ! OUT (1/CH-1/CD_STD)/VSHR on land
!                                  !     tiles
     &,ALPHA1(LAND_POINTS,NTILES)                                       &
                                      ! OUT Mean gradient of saturated
!                                  !     specific humidity with respect
!                                  !     to temperature between the
!                                  !     bottom model layer and tile
!                                  !     surfaces
     &,ASHTF_TILE(LAND_POINTS,NTILES)                                   &
                                      !OUT Coefficient to calculate
!                                  !     surface heat flux into land
!                                  !     tiles.
     &,FQT_TILE(LAND_POINTS,NTILES)                                     &
                                      ! OUT Surface FQT for land tiles
     &,EPOT_TILE(LAND_POINTS,NTILES)                                    &
                                      ! OUT Local EPOT for land tiles.
     &,FQT_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface FQT for sea-ice
     &,FTL_ICE(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Surface FTL for sea-ice
     &,FRACA(LAND_POINTS,NTILES)                                        &
                                      ! OUT Fraction of surface moisture
!                                  !     flux with only aerodynamic
!                                  !     resistance for snow-free land
!                                  !     tiles.
     &,RESFS(LAND_POINTS,NTILES)                                        &
                                      ! OUT Combined soil, stomatal
!                                  !     and aerodynamic resistance
!                                  !     factor for fraction (1-FRACA)
!                                  !     of snow-free land tiles.
     &,RESFT(LAND_POINTS,NTILES)                                        &
                                      ! OUT Total resistance factor.
!                                  !     FRACA+(1-FRACA)*RESFS for
!                                  !     snow-free land, 1 for snow.
     &,RHOKH_TILE(LAND_POINTS,NTILES)                                   &
                                      ! OUT Surface exchange coefficient
!                                  !     for land tiles
     &,RHOKH_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! OUT Surface exchange coefficients
!                                  !     for sea and sea-ice
     &,RHOKPM(LAND_POINTS,NTILES)                                       &
                                      ! OUT Surface exchange coefficient
     &,RHOKPM_POT(LAND_POINTS,NTILES)                                   &
                                      ! OUT Potential evaporation
!                                     !     exchange coefficient.
     &,RHOKPM_SICE(ROW_LENGTH,ROWS)                                     &
                                   ! OUT Sea-ice surface exchange coeff.
     &,Z0HSSI(ROW_LENGTH,ROWS)                                          &
     &,Z0MSSI(ROW_LENGTH,ROWS)                                          &
                                   ! OUT Roughness lengths over sea (m).
     &,Z0H_TILE(LAND_POINTS,NTILES)                                     &
                                   ! OUT Tile roughness lengths for h
!                                  !     and moisture (m).
     &,Z0M_GB(ROW_LENGTH,ROWS)                                          &
                                   ! OUT Grid-box mean roughness length
!                                  !      for momentum (m).
     &,Z0M_TILE(LAND_POINTS,NTILES)                                     &
                                   ! OUT Tile roughness lengths for
!                                  !     momentum.
     &,CHR1P5M(LAND_POINTS,NTILES)                                      &
                                   ! OUT Ratio of coefffs for
!                                  !     calculation of 1.5m temp for
!                                  !     land tiles.
     &,CHR1P5M_SICE(ROW_LENGTH,ROWS)                                    &
!                                  ! OUT CHR1P5M for sea and sea-ice
!                                  !     (leads ignored).
     &,SMC(LAND_POINTS)                                                 &
                                      ! OUT Available moisture in the
!                                  !     soil profile (mm).
     &,GPP(LAND_POINTS)                                                 &
                                      ! OUT Gross primary productivity
!                                  !     (kg C/m2/s).
     &,NPP(LAND_POINTS)                                                 &
                                      ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
     &,RESP_P(LAND_POINTS)                                              &
                                      ! OUT Plant respiration (kg C/m2/s
     !kdcorbin, 11/10 - changed from NPFT, currently setup for Leaf Resp
     &,G_LEAF(LAND_POINTS,NTILES)                                       &
                                      ! OUT Leaf turnover rate (/360days)
     !kdcorbin, 11/10 - changed from NPFT
     &,GPP_FT(LAND_POINTS,NTILES)                                       &
                                      ! OUT Gross primary productivity
                                      !     on PFTs (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,NPP_FT(LAND_POINTS,NTILES)                                       &
                                      ! OUT Net primary productivity
                                      !     (kg C/m2/s).
     !kdcorbin, 11/10 - changed from NPFT
     &,RESP_P_FT(LAND_POINTS,NTILES)                                    &
                                      ! OUT Plant respiration on PFTs
                                      !     (kg C/m2/s).
     &,RESP_S(LAND_POINTS,DIM_CS1)                                      &
                                  ! OUT Soil respiration (kg C/m2/s)
     &,RESP_S_TOT(DIM_CS2)                                              &
                                 ! OUT total soil respiration
     &,RESP_W_FT(LAND_POINTS,NPFT)                                      &
                                      ! OUT Wood maintenance respiration
!                                  !     (kg C/m2/s).
     &,GC(LAND_POINTS,NTILES)                                           &
                                      ! OUT "Stomatal" conductance to
!                                  !      evaporation for land tiles
!                                  !      (m/s).
     &,CANHC_TILE(LAND_POINTS,NTILES)                                   &
                                      ! IN Areal heat capacity of canopy
!                                  !    for land tiles (J/K/m2).
     &,WT_EXT_TILE(LAND_POINTS,DSM_LEVELS,NTILES)                       &
!                                  ! IN Fraction of evapotranspiration
!                                  !    which is extracted from each
!                                  !    soil layer by each tile.
     &,FLAKE(LAND_POINTS,NTILES)                                        &
                                      ! IN Lake fraction.
     &,TILE_FRAC(LAND_POINTS,NTILES)                                    &
                                      ! OUT Tile fractions including
!                                  !     snow cover in the ice tile.
     &,FSMC(LAND_POINTS,NPFT)                                           &
                                      ! OUT Moisture availability
!                                  !     factor.
     &,FLANDG_U(ROW_LENGTH,ROWS)                                        &
                                   ! OUT Land frac (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,FLANDG_V(ROW_LENGTH,N_ROWS)                                      &
                                   ! OUT Land frac (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
     &,RIB_SSI(ROW_LENGTH,ROWS)                                         &
                                   ! OUT Sea mean bulk Richardson number
                                   !     for lowest layer.
     &,TAUX_LAND(ROW_LENGTH,ROWS)                                       &
                                   ! OUT W'ly component of land
!                                  !     surface wind stress (N/sq m).
!                                  !     (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
     &,TAUX_SSI(ROW_LENGTH,ROWS)                                        &
                                   ! OUT W'ly component of mean sea
!                                  !     surface wind stress (N/sq m).
!                                  !     (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
     &,TAUY_LAND(ROW_LENGTH,N_ROWS)                                     &
                                   ! OUT S'ly component of land
!                                  !     surface wind stress (N/sq m).
!                                  !     On V-grid;
!                                  !     comments as per TAUX.
     &,TAUY_SSI(ROW_LENGTH,N_ROWS)                                      &
                                   ! OUT S'ly component of mean sea
!                                  !     surface wind stress (N/sq m).
!                                  !     On V-grid;
!                                  !     comments as per TAUX.
     &, vshr_land(row_length, rows)                                     &
     &, vshr_ssi(row_length, rows)                                      &
     &,ANTHROP_HEAT(NTILES)        ! IN  Additional heat source on 
                                   !     tiles for anthropogenic urban
                                   !     heat source  (W/m2)



! local variables.

! loop counters
      Integer                                                           &
     &  i, j, k                                                         &
     &,l


! Diagnostic switches
      Logical                                                           &
     &  su10, sv10, slh, sq1p5, sT1p5, sq_T1p5                          &
     &, sfme, simlt, smlt, sz0heff

! local variables

! Local data arrays

      Real                                                              &
     &  T(row_length, rows, bl_levels)

      Real                                                              &
     & Z_land(row_length, rows)  !land height over
!                                !fractional land points
! External Routines:
      External bdy_layr
      External timer


! error information

!      print *,' ni_bl_ctl hcon',hcon,soil_type,veg_type,t_soil
      If ( error_code  ==  0) Then
! ----------------------------------------------------------------------
! Section BL.0 Initialisation of variables.
! ----------------------------------------------------------------------


! Set diagnostic flags required for boundary layer diagnostics from
! STASHflags.
!        ! --------------------------------------------------------
!        ! Note that an equivalent block of code exists in routine
!        ! ni_imp_ctl, and needs to be kept consistent.
!        ! --------------------------------------------------------
!        ! Windspeed (227) and u, v at 10m on 'B' or 'C' grid
         su10 = sf(227,3) .or. sf(225,3) .or. sf(209,3) .or. sf(463,3)
         sv10 = sf(227,3) .or. sf(226,3) .or. sf(210,3) .or. sf(463,3)
         slh = sf(234,3)
         sq_T1p5 = sf(236,3) .or. sf(237,3) .or. sf(245,3)              &
     &        .or. sf(247,3) .or. sf(248,3) .or. sf(250,3)              &
     &        .or. L_scrn                                               &
     &        .or. sf(341,3) .or. sf(342,3)                             &
     &        .or. sf(253,3) .or. sf(328,3) .or. sf(329,3)
         sq1p5 = sq_T1p5
         sT1p5 = sq_T1p5
         sfme  = sf(224,3)
         simlt = sf(235,3)
         smlt  = sf(258,3)
         sz0heff = sf(027,3)


!
!       Allocate space for those BL diagnostic arrays required and zero
!       the elements explicitly
!
! Obukhov length
! Also required for gustiness diagnostic (463)

        BL_diag%L_oblen    = sf(464,3) .or. sf(463,3)
        If (BL_diag%L_oblen) then
          allocate(BL_diag%oblen(row_length,rows))
          BL_diag%oblen(:,:) = 0.0
        End if

! Friction velocity
! Also required for gustiness diagnostic (463)

        BL_diag%L_ustar   = sf(465,3) .or. sf(463,3)
        If (BL_diag%L_ustar) then
          allocate(BL_diag%ustar(row_length,rows))
          BL_diag%ustar(:,:) = 0.0
        End if
! Surface buoyancy flux
        BL_diag%L_wbsurf    = sf(467,3)
        If (BL_diag%L_wbsurf) then
          allocate(BL_diag%wbsurf(row_length,rows))
          BL_diag%wbsurf(:,:) = 0.0
        End if
! Gradient Richardson number
        BL_diag%L_gradrich    = sf(468,3)
        If (BL_diag%L_gradrich) then
          allocate(BL_diag%gradrich(row_length,rows,bl_levels))
          BL_diag%gradrich(:,:,:) = 0.0
        End if
! Convective velocity scale
        BL_diag%L_wstar   = sf(466,3)
        If (BL_diag%L_wstar) then
          allocate(BL_diag%wstar(row_length,rows))
          BL_diag%wstar(:,:) = 0.0
        End if
! Stratification
        BL_diag%L_dbdz    = sf(469,3)
        If (BL_diag%L_dbdz) then
          allocate(BL_diag%dbdz(row_length,rows,bl_levels))
          BL_diag%dbdz(:,:,:) = 0.0
        End if
! Modulus of shear
        BL_diag%L_dvdzm    = sf(470,3)
        If (BL_diag%L_dvdzm) then
          allocate(BL_diag%dvdzm(row_length,rows,bl_levels))
          BL_diag%dvdzm(:,:,:) = 0.0
        End if
! Momentum diffusivity
        BL_diag%L_rhokm    = sf(471,3)
        If (BL_diag%L_rhokm) then
          allocate(BL_diag%rhokm(row_length,rows,bl_levels))
          BL_diag%rhokm(:,:,:) = 0.0
        End if
! Thermal diffusivity
        BL_diag%L_rhokh    = sf(472,3)
        If (BL_diag%L_rhokh) then
          allocate(BL_diag%rhokh(row_length,rows,bl_levels))
          BL_diag%rhokh(:,:,:) = 0.0
        End if
! Turbulent kinetic energy
        BL_diag%L_tke    = sf(473,3)
        If (BL_diag%L_tke) then
          allocate(BL_diag%tke(row_length,rows,bl_levels))
          BL_diag%tke(:,:,:) = 0.0
        End if
! x component of orographic stress
        BL_diag%L_ostressx    = sf(474,3)
        If (BL_diag%L_ostressx) then
          allocate(BL_diag%ostressx(row_length,rows,bl_levels))
          BL_diag%ostressx(:,:,:) = 0.0
        End if
! y component of orographic stress
        BL_diag%L_ostressy    = sf(475,3)
        If (BL_diag%L_ostressy) then
          allocate(BL_diag%ostressy(row_length,rows,bl_levels))
          BL_diag%ostressy(:,:,:) = 0.0
        End if
! 
        BL_diag%L_dscbase    = sf(360,3)
        If (BL_diag%L_dscbase) then
          allocate(BL_diag%dscbase(row_length,rows))
          BL_diag%dscbase(:,:) = 0.0
        End if
! 
        BL_diag%L_cldbase    = sf(361,3)
        If (BL_diag%L_cldbase) then
          allocate(BL_diag%cldbase(row_length,rows))
          BL_diag%cldbase(:,:) = 0.0
        End if
! 
        BL_diag%L_weparm    = sf(362,3)
        If (BL_diag%L_weparm) then
          allocate(BL_diag%weparm(row_length,rows))
          BL_diag%weparm(:,:) = 0.0
        End if
! 
        BL_diag%L_weparm_dsc    = sf(363,3)
        If (BL_diag%L_weparm_dsc) then
          allocate(BL_diag%weparm_dsc(row_length,rows))
          BL_diag%weparm_dsc(:,:) = 0.0
        End if

!
! Code to allocate unity arrays when not
! used (for portability)
!

        If (.NOT. BL_diag%L_oblen) then
          allocate(BL_diag%oblen(1,1))
          BL_diag%oblen(:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_ustar) then
          allocate(BL_diag%ustar(1,1))
          BL_diag%ustar(:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_wbsurf) then
          allocate(BL_diag%wbsurf(1,1))
          BL_diag%wbsurf(:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_gradrich) then
          allocate(BL_diag%gradrich(1,1,1))
          BL_diag%gradrich(:,:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_wstar) then
          allocate(BL_diag%wstar(1,1))
          BL_diag%wstar(:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_dbdz) then
          allocate(BL_diag%dbdz(1,1,1))
          BL_diag%dbdz(:,:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_dvdzm) then
          allocate(BL_diag%dvdzm(1,1,1))
          BL_diag%dvdzm(:,:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_rhokm) then
          allocate(BL_diag%rhokm(1,1,1))
          BL_diag%rhokm(:,:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_rhokh) then
          allocate(BL_diag%rhokh(1,1,1))
          BL_diag%rhokh(:,:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_tke) then
          allocate(BL_diag%tke(1,1,1))
          BL_diag%tke(:,:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_ostressx) then
          allocate(BL_diag%ostressx(1,1,1))
          BL_diag%ostressx(:,:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_ostressy) then
          allocate(BL_diag%ostressy(1,1,1))
          BL_diag%ostressy(:,:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_dscbase) then
          allocate(BL_diag%dscbase(1,1))
          BL_diag%dscbase(:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_cldbase) then
          allocate(BL_diag%cldbase(1,1))
          BL_diag%cldbase(:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_weparm) then
          allocate(BL_diag%weparm(1,1))
          BL_diag%weparm(:,:) = 0.0
        End if

        If (.NOT. BL_diag%L_weparm_dsc) then
          allocate(BL_diag%weparm_dsc(1,1))
          BL_diag%weparm_dsc(:,:) = 0.0
        End if


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

! ----------------------------------------------------------------------
! Calculate the land height over fractional land points:
! ----------------------------------------------------------------------
        l=0
        Do j = 1, rows
           Do i = 1, row_length
             Z_land(i,j) = 0.0

             If(land_sea_mask(i,j))Then
               l=l+1

               If(l_ctile.and.Fland(l) >  0.0.and.Fland(l) <  1.0)Then
                 Z_land(i,j) = r_theta_levels(i,j,0) - earth_radius
                 If(Z_land(i,j) <  0.0)Z_land(i,j) = 0.0
               Endif

             Endif
           End Do
        End Do

! ----------------------------------------------------------------------
! Section BL.2a Call boundary_layer scheme.
! ----------------------------------------------------------------------

!      print *,' ni_bl_ctl bef BL_INTCT'
!      print *,'atmos physic2 lat,lon',lat,long
!      print *,'surf_d 0',surf_down_sw
!      print *,'surf_d 1',alb_tile
!      print *,'surf_d 2',cos_zenith_angle
!      print *,'ls_rain', ls_rain,ls_snow
!      print *,'lat,long,day,time_sec',lat,long,day,time_sec
!      print *,'SN_D',  SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L
!      print *,'SN_R', SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE
!      print *,'TSOIL_T', TSOIL_TILE,T_SURF_TILE,HCONS
!      print *,'ISNOW+F', ISNOW_FLG3L
!       print *,'bef bl_intct snow_tile',SNOW_TILE


! DEPENDS ON: bl_intct
         Call BL_INTCT( halo_i, halo_j, off_x, off_y, row_length, rows, &
     &  n_rows,global_row_length,proc_row_group,at_extremity,n_proc,    &
     &  n_procx,n_procy,neighbour,Substep_Number,Num_Substeps,          &
     &  L_phys2_substep,NumCycles,CycleNo,ntiles,land_points,           &

! in values defining model domain, vertical grid of model atmosphere
!    and implicit weights. :
     &        model_domain, bl_levels,                                  &
     &        r_rho_levels, r_theta_levels, eta_theta_levels,           &
     &        p, p_layer_centres(1,1,1), rho_rsq, rho_wet, rho_dry,     &
     &        alpha_cd, sin_theta_longitude, cos_theta_longitude,       &
     &        sin_theta_latitude,                                       &

! IN U and V momentum fields.
     &        u, v, w, etadot, u_conv, v_conv,                          &
! IN Non turbulent increments to momentum and thermodynamic variables.
!  (New dynamics only).
     &        R_u, R_v,                                                 &
! variables for subgrid turbulence scheme
     & visc_BL_m, delta_lambda, delta_phi, FM_3D,FH_3D,L_subfilter_vert,&
     & L_subfilter_horiz, L_subfilter_blend,max_diff,turb_startlev_vert,&
     & turb_endlev_vert, BL_COEF_KM, BL_COEF_KH,                        &
! in soil/vegetation/land surface data :
     &        land_sea_mask, land_index,                                &
     &        dst_levels, dsm_levels,                                   &
     &        DIM_CS1, DIM_CS2,                                         &
     &        hcon,                                                     &
     &        smvccl, smvcst, smvcwt, sthf, sthu,                       &
     &        sil_orog_land,                                            &
     &        formdrag, orog_drag_param, ho2r2_orog,                    &
     & SOIL_LAYER_MOISTURE,                                             &

! in sea/sea-ice data :
     & ice_fract,u_0,v_0,u_0_p,v_0_p,Charnock,SeaSalinityFactor,        &

! in cloud data :
! cloud_fraction passed in
     &        bulk_cloud_fraction, q, qcf, qcl,                         &
     &        qcf_latest, qcl_latest, cca, ccb, cct, T,                 &

! in everything not covered so far :
     &        L_LAMBDAM2, L_FULL_LAMBDAS,                               &
     &        CO2_MMR, photosynth_act_rad, p_star, surf_radflux,        &
     &        rad_hr, micro_tends, timestep, L_mixing_ratio, zh_prev,   &
     &        L_ctile,BL_OPTIONS,                                       &
     &        L_SBLeq,L_SBLco,Muw_SBL,Mwt_SBL,                          &

! IN Variables for: prescribed surface flux forcing
!                   prescribed sea-surface roughness lengths
     &        flux_e, flux_h, L_flux_bc, L_spec_z0, Z0M_SCM, Z0H_SCM,   &

! IN variables required for mineral dust scheme
     & L_DUST, L_CAM_DUST,                                              &
     & SOIL_CLAY,SOIL_SILT,SOIL_SAND,                                   &
     & DUST_MREL1,DUST_MREL2,DUST_MREL3,                                &
     & DUST_MREL4,DUST_MREL5,DUST_MREL6,                                &

! IN additional variables for MOSES II
     &        LAND_PTS_TRIF,NPFT_TRIF,                                  &
     &        CANOPY,CATCH,CATCH_SNOW,                                  &
     &        SNOW_TILE,Z0_TILE,LW_DOWN,                                &
     &        SW_TILE,TSTAR_TILE,                                       &
     &        CO2_3D,CO2_DIM_LEN,CO2_DIM_ROW,L_CO2_INTERACTIVE,         &
     &        L_PHENOL,L_TRIFFID,L_Q10,ASTEPS_SINCE_TRIFFID,CAN_MODEL,  &
     &        CS,FRAC,CANHT_FT,LAI_FT,                                  &
     &        FLAND,FLANDG,TSTAR_SEA,Z_LAND,ANTHROP_HEAT,               &
     &        ALBSOIL,COS_ZENITH_ANGLE,CAN_RAD_MOD,ILAYERS,             &
    
!    EAK      
!    IN       
     &        surf_down_sw,alb_tile,l_tile_pts,                         &
!     &        surf_down_sw,alb_tile,cos_zenith_angle,        &
     &        ls_rain,ls_snow,SW_DOWN,                                  &
     &        lat,long,day,time_sec,                                    &
     &  SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,                  &
     &  SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,            &
     &  TSOIL_TILE,T_SURF_TILE,HCONS,BEXP,                              &
     &  SATHH,SATCON,HCAP,                                              &
     & SOIL_TYPE,VEG_TYPE,                                              &
     &  ISNOW_FLG3L,total_nsteps,                                       &
     &            FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB,              &
     &            TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB,          &
     &            USTAR_CAB,SURF_HTF_CAB,                               &
     &            U_S_CAB,CH_CAB,CD_CAB,                                &
     &            l_cable,                                              &
!EAK
     &            SNAGE_TILE,RTSOIL_TILE,                               &
     &            GFLUX_TILE,SGFLUX_TILE ,                              &
     &  CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                    &
     &  NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                          &

! INOUT data :
     &        t_soil, ti, t_surf, z0msea,                               &
! Variables for STPH_RP (G0_RP,par_mezcla)
     &        T_latest, Q_latest,G0_RP, par_mezcla,                     &

! INOUT additional variables for MOSES II
     &        GS,G_LEAF_ACC,NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC,        &
     &        TSTAR_SICE,TSTAR_SSI,                                     &

! out  diagnostic not requiring stash flags :
     &        cd, ch, e_sea, fqT,                                       &
     &        ftl, h_sea, rhokh, rhokm,                                 &
     &        rhokm_u,rhokm_v,                                          &
     &        rib_gb,                                                   &
     &        taux, tauy, vshr,                                         &
     &        zht,                                                      &
     &        shallowc,cu_over_orog,                                    &
     &     bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6, &
     &     bl_type_7,                                                   &
     &        z0m_eff_gb,z0h_eff_gb,                                    &

! out diagnostic requiring stash flags :
     &        fme,WT_EXT,RA,                                            &

! (in) stash flags :-
     &        sfme, simlt, smlt, slh,                                   &
     &        sq1p5, st1p5, su10, sv10, sz0heff,                        &
!
! SCM Diagnostics (dummy values in full UM)
     &        nSCMDpkgs,L_SCMDiags,                                     &

! out data required for tracer mixing :
     &        rho_aresist,aresist,resist_b,                             &

!OUT variables required for mineral dust scheme
     & R_B_DUST,DUST_FLUX,                                              &
     & U_S_T_TILE,U_S_T_DRY_TILE,U_S_STD_TILE,                          &
     &        ntml,                                                     &
     &        KENT, WE_LIM, T_FRAC, ZRZI,                               &
     &        KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC,               &
     &        ZHSC, Z_HALF,                                             &


! OUT data required for 4D-VAR :
     &        rho_cd_modv1, rho_km_var,                                 &

! OUT variables required in IMP_SOLVER
     &        alpha1_sice, ashtf,                                       &
     &        dtrdz_charney_grid, rdz_charney_grid,                     &
     &        dtrdz_u, dtrdz_v, rdz_u, rdz_v,                           &
     &        cdr10m_u, cdr10m_v,                                       &
     &        z1_tq,                                                    &

! OUT variables which need to be maintained between substeps
     &        RHO_UV,RHO_TQ,DZL_CHARNEY,RDZ,                            &
     &        Z1_UV,Z_FULL,Z_UV,Z_TQ,P_HALF,DELTAP,                     &

! OUT additional variables for MOSES II
     &        FTL_TILE,LE_TILE,RADNET_SICE,                             &
     &        RADNET_TILE,RIB_TILE,RHO_ARESIST_TILE,ARESIST_TILE,       &
     &RESIST_B_TILE,ALPHA1,ASHTF_TILE,FQT_TILE,EPOT_TILE,FQT_ICE,       &
     &FTL_ICE,FRACA,RESFS,RESFT,RHOKH_TILE,RHOKH_SICE,RHOKPM,RHOKPM_POT,&
     &RHOKPM_SICE,Z0HSSI,Z0H_TILE,Z0M_GB,Z0MSSI,Z0M_TILE,CHR1P5M,       &
     &        CHR1P5M_SICE,SMC,GPP,NPP,RESP_P,G_LEAF,GPP_FT,NPP_FT,     &
     &        RESP_P_FT,RESP_S,RESP_S_TOT,RESP_W_FT,                    &
     &        GC,CANHC_TILE,WT_EXT_TILE,FLAKE,                          &
     &        TILE_INDEX,TILE_PTS,TILE_FRAC,FSMC,                       &
     &        TSTAR_LAND,RIB_SSI,TAUX_LAND,TAUX_SSI,TAUY_LAND,TAUY_SSI, &
     &        VSHR_LAND,VSHR_SSI,                                       &
     &        FLANDG_U,FLANDG_V,                                        &

! out data required elsewhere in um system :
     &        zh, t1_sd, q1_sd,NTPAR,NLCL,ZHPAR,ZLCL,L_SHALLOW,         &
     &        ntdsc, nbdsc, cumulus,WSTAR,WTHVS,DELTHVU,uw0,vw0,        &
     & L_ukca,                                                          &
     & Error_code,Ltimer,BL_diag,                                       &
      ! end step of experiment, this step, step width, processor num
     & endstep, timestep_number, mype )
!print *,'aft bl_intct snow_tile',SNOW_TILE


      End If                    ! on error code = 0

!      print *,'after ni_bl_ctl'
! end of routine NI_bl_ctl
      Return
      END SUBROUTINE NI_bl_ctl
