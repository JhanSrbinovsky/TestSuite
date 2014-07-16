
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_conv

      Subroutine diagnostics_conv(                                      &
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
     &,                      Theta, theta_inc, q, q_inc                 &
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
     &,                      dubydt_p,dvbydt_p,up_flux_half,up_flux     &
     &,                      dwn_flux                                   &
     &,                      entrain_up,detrain_up,entrain_dwn          &
     &,                      detrain_dwn,uw_dp,vw_dp,uw_shall,vw_shall  &
     &,                 wqt_flux,wthetal_flux,wthetav_flux,wql_flux     &
     &,                 wqt_cb,wthetal_cb,wqt_inv,wthetal_inv           &
     &,                 sh_top,sh_base                                  &
     &,               shallowc,congestus_ind,congestus_ind2,cu_over_orog&
     &,               cape, cape_undilute, cin_undilute                 &
     &,               mid_level,ntml_diag,ntpar_diag,freeze_diag        &
     &,               kterm_deep,precip_deep,precip_shall,precip_mid    &
     &,               precip_cong, wstar_dn,wstar_up,mb1,mb2            &
     &,               cg_term,cg_top,cg_base                            &
     &,              mf_deep,mf_congest,mf_shall,mf_midlev              &
     &,              dt_deep,dt_congest,dt_shall,dt_midlev              &
     &,              dq_deep,dq_congest,dq_shall,dq_midlev              &
     &,              du_deep,du_congest,du_shall,du_midlev              &
     &,              dv_deep,dv_congest,dv_shall,dv_midlev              &
     &,                      lcbase, lctop, lcca, n_cca_levels          &
     &,                      CONSCAV_DUST1,CONSCAV_DUST2,CONSCAV_DUST3  &
     &,                      CONSCAV_DUST4,CONSCAV_DUST5,CONSCAV_DUST6  &
     &,                      L_DUST                                     &
     &,                      conwash_nh3                                &
     &,                      conwash_so2, conscav_so4ait                &
     &,                      conscav_so4acc, conscav_so4dis             &
     &,                      conscav_agedsoot, conscav_agedbmass        &
     &,                      conscav_agedocff                           &
     &,                      timestep, Num_Substeps                     &
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &     STASHwork                                                    &
     &     )

     Use ac_diagnostics_mod, Only :                                     &
         CVRR, CVSR, CONVCC, TINC_CVN

     Use cv_cntl_mod, Only:                                             &
         lcv_3d_cca

     Use cv_run_mod,  Only:                                             &
         tower_factor
!
! Description:
!   Calculates convection-related diagnostics (held in STASH section 5).
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Intercepted arrays and diagnostic arrays are input from the
!   convection routines called previously. Each diagnostic is simply
!   copied into the STASHwork array to be passed on to STASH for
!   output processing.
!
!   Diagnostics currently available (in order calculated):
!   Item  Description
!    163  Liquid Cloud Condensate increment (inhomog)
!    164  Frozen Cloud Condensate increment (inhomog)
!    172  Total  Cloud Volume increment (inhomog)
!    173  Liquid Cloud Volume increment (inhomog)
!    174  Frozen Cloud Volume increment (inhomog)
!    181  Temperature increment
!    182  Humidity increment
!    183  Liquid Cloud Condensate increment
!    184  Frozen Cloud Condensate increment
!    185  u wind increment
!    186  v wind increment
!    187  Temperature increment in no shallow convection columns
!    188  Humidity increment in no shallow convection columns
!    192  Total  Cloud Volume increment
!    193  Liquid Cloud Volume increment
!    194  Frozen Cloud Volume increment
!    201  Convective rainfall
!    202  Convective snowfall
!    205  Convective rainfall rates
!    206  Convective snowfall rates
!    212  Convective cloud amount
!    213  Convective cloud liquid water
!    226  Total precipitation
!    214  Total rain rate
!    215  Total snow rate
!    216  Total precipitation rate
!    217  CAPE
!    207  Pressure at convective cloud base
!    208  Pressure at convective cloud top
!    210  Convective cloud base (ICAO standard atmosphere heights)
!    211  Convective cloud top (ICAO standard atmosphere heights)
!    218  Lowest convective cloud base
!    219  Lowest convective cloud top
!    220  Lowest convective cloud amount
!    222  Pressure at lowest convective cloud base
!    223  Pressure at lowest convective cloud top
!    224  Lowest convective cloud base(ICAO standard atmosphere heights)
!    225  Lowest convective cloud top (ICAO standard atmosphere heights)
!    227  3D Convective rainfall rate 
!    228  3D Convective snowfall rate 
!    262  2D Convective cloud amount
!    209  Temperature at end of timestep
!     10  q at end of timestep
!    233  undilute parcel CAPE
!    234  undilute parcel CIN
!    235  u latest
!    236  v latest
!    237  NH3 SCAVENGED BY CONV PPN KG/M2/SEC
!    238  SO2 SCAVENGED BY CONV PPN KG/M2/SEC
!    239  SO4 AIT SCAVNGD BY CONV PPN KG/M2/S
!    240  SO4 ACC SCAVNGD BY CONV PPN KG/M2/S
!    241  SO4 DIS SCAVNGD BY CONV PPN KG/M2/S
!    242  Aged soot scavenged by convective precipitation kg/m2/s
!    243  Aged biomass scavenged by convective precipitation kg/m2/s
!    244  Aged fossil-fuel OC scavenged by convective precipitation kg/m2/s
!    246  Component of UPDRAUGHT MASS FLUX ON HALF (RHO) LEVELS (Pa/s) needed for 4D-VAR
!    249  UPDRAUGHT MASS FLUX ON HALF (RHO) LEVELS (Pa/s)
!    250  UPDRAUGHT MASS FLUX (Pa/s)
!    251  DOWNDRAUGHT MASS FLUX
!    252  UPDRAUGHT ENTRAINMENT RATE
!    253  UPDRAUGHT DETRAINMENT RATE
!    254  DOWNDRAUGHT ENTRAINMENT RATE
!    255  DOWNDRAUGHT DETRAINMENT RATE
! NEW DIAGNOSTICS FOR 4A CONVECTION SCHEME ONLY
!    256  U INCREMENT (P GRID)
!    257  V INCREMENT (P GRID)
!    258  DEEP UW STRESS (P GRID)
!    259  DEEP VW STRESS (P GRID)
!    260  SHALLOW UW STRESS (P GRID)
!    261  SHALLOW UW STRESS (P GRID)
!    270  SHALLOW CUMULUS INDICATOR
!    271  CU OVER OROGRAPHY INDICATOR
!    272  Mid level onvection indicator
!    273  NTML
!    274  NTPAR
!    275  Freezing level number
!    276  KTERM deep
!    277  total precipitation from deep convection
!    278  total precipitation from shallow convection
!    279  total precipitation from mid-level convection
!    280  total precipitation from congestus convection
! NEW DIAGNOSTICS FOR 5A turbulence part of convection scheme only
!    290  w'qt' flux
!    291  w'ql' flux
!    292  w'thetal' flux
!    293  w'thetav' flux
!    300  wstar_dn
!    301  wstar_up
!    302  cloud base mass flux  from wstar_dn
!    303  alternative cloud base mass flux
!    304  cloud base w'qt'
!    305  cloud base w'thetal'
!    306  inversion w'qt'
!    307  inversion w'thetal'
!    308  base of shallow convection (m)
!    309  top of shallow convection (m)
!
!    310  Congestus indicator
!    311  Congestus indicator 2
!    312  Termination level for congestus
!    313  height of top of congestus
!    314  height of base of congestus
!
!    320  mass flux deep convection
!    321  mass flux congestus convection
!    322  mass flux shallow convection
!    323  mass flux mid_lev convection
!    324  dT deep convection
!    325  dT congestus convection
!    326  dT shallow convection
!    327  dT mid_lev convection
!    328  dq deep convection
!    329  dq congestus convection
!    330  dq shallow convection
!    331  dq mid_lev convection
!    332  du deep convection
!    333  du congestus convection
!    334  du shallow convection
!    335  du mid_lev convection
!    336  dv deep convection
!    337  dv congestus convection
!    338  dv shallow convection
!    339  dv mid_lev convection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      Implicit None

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN. ie: Input variables.
      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held
     &, n_cca_levels                                                    &
     &, Num_Substeps

      Logical                                                           &
     & L_DUST

      Integer                                                           &
     &  global_row_length                                               &
                            !IN. NUMBER OF points on a global row
     &, global_rows                                                     &
                            !IN. NUMBER OF global rows
     &, me                                                              &
                            !IN. Processor number
     &, halo_i                                                          &
                            !IN. size of large halo in x direction
     &, halo_j                                                          &
                            !IN. size of large halo in y direction
     &, off_x                                                           &
                            !IN. size of small halo in x direction
     &, off_y                                                           &
                            !IN. size of small halo in y direction
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)

      Real                                                              &
     &  timestep


! Primary Arrays used in all models
      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &      model_levels)                                               &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)

      Real                                                              &
     &  Theta(row_length, rows, model_levels)                           &
     &, Theta_inc(row_length, rows, model_levels)                       &
     &, q(row_length,rows, wet_model_levels)                            &
     &, q_inc(row_length,rows, wet_model_levels)                        &
     &, qcl_inc(row_length,rows, wet_model_levels)                      &
     &, qcf_inc(row_length,rows, wet_model_levels)                      &
     &, cf_liquid_inc(row_length,rows, wet_model_levels)                &
     &, cf_frozen_inc(row_length,rows, wet_model_levels)                &
     &, bulk_cf_inc(row_length, rows, wet_model_levels)                 &
     &, ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, ccw(row_length, rows, wet_model_levels)                         &
     &, conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)                                     &
     &, conv_rain_3d(row_length, rows, wet_model_levels)                & 
     &, conv_snow_3d(row_length, rows, wet_model_levels)                & 
     &, cca(row_length, rows, n_cca_levels)                             &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
! wind increments on p grid (a05_4a)
     &, dubydt_p(row_length,rows,model_levels)                          &
     &, dvbydt_p(row_length,rows,model_levels)                          &
! mass flux and entrainment/detrainment diagnostics
     &, up_flux(row_length,rows,model_levels)                           &
     &, up_flux_half(row_length,rows,model_levels)                      &
     &, dwn_flux(row_length,rows,model_levels)                          &
     &, entrain_up(row_length,rows,model_levels)                        &
     &, detrain_up(row_length,rows,model_levels)                        &
     &, entrain_dwn(row_length,rows,model_levels)                       &
     &, detrain_dwn(row_length,rows,model_levels)                       &
! convective momentum transport diagnostics (a05_4a)
     &, uw_dp(row_length,rows,model_levels)                             &
     &, vw_dp(row_length,rows,model_levels)                             &
     &, uw_shall(row_length,rows,model_levels)                          &
     &, vw_shall(row_length,rows,model_levels)                          &
! Type diagnostics (a05_4a)
     &, shallowc(row_length,rows)                                       &
                                  ! IN Indicator set to 1.0 if shallow,
!                                 !   0.0 if not shallow or not cumulus
     &, congestus_ind(row_length,rows)                                  &
                                       ! Indicator, 1.0 if congestus
     &, congestus_ind2(row_length,rows)                                 &
                                        ! Indicator 2, 1.0 if congestus
     &, cu_over_orog(row_length,rows)                                   &
!                                  ! IN Indicator for cumulus
!                                  !     over steep orography
!                                  !   Indicator set to 1.0 if true,
!                                  !   0.0 if false. Exclusive.
! diagnostic increment fields
!     Intent INOUT , IN initial value in convection, OUT increment.
     &, qcl_incr_inhom_diag(row_length,rows,wet_model_levels)           &
     &, qcf_incr_inhom_diag(row_length,rows,wet_model_levels)           &
     &, bulk_cf_incr_inhom_diag(row_length, rows, wet_model_levels)     &
     &, cf_liquid_incr_inhom_diag(row_length,rows,wet_model_levels)     &
     &, cf_frozen_incr_inhom_diag(row_length,rows,wet_model_levels)     &
     &, T_incr_diagnostic(row_length,  rows,    model_levels)           &
     &, q_incr_diagnostic(row_length,  rows,wet_model_levels)           &
     &, qcl_incr_diagnostic(row_length,rows,wet_model_levels)           &
     &, qcf_incr_diagnostic(row_length,rows,wet_model_levels)           &
     &, u_incr_diagnostic(row_length,  rows,    model_levels)           &
     &, v_incr_diagnostic(row_length,n_rows,    model_levels)           &
     &, cf_liquid_incr_diagnostic(row_length,rows,wet_model_levels)     &
     &, cf_frozen_incr_diagnostic(row_length,rows,wet_model_levels)     &
     &, bulk_cf_incr_diagnostic(row_length, rows, wet_model_levels)     &
! mineral dust diagnostics
     &, CONSCAV_DUST1(ROW_LENGTH, ROWS)                                 &
                                           !COL TOTAL SCVNGD DIV1 DUST
     &, CONSCAV_DUST2(ROW_LENGTH, ROWS)                                 &
                                           !COL TOTAL SCVNGD DIV2 DUST
     &, CONSCAV_DUST3(ROW_LENGTH, ROWS)                                 &
                                           !COL TOTAL SCVNGD DIV3 DUST
     &, CONSCAV_DUST4(ROW_LENGTH, ROWS)                                 &
                                           !COL TOTAL SCVNGD DIV4 DUST
     &, CONSCAV_DUST5(ROW_LENGTH, ROWS)                                 &
                                           !COL TOTAL SCVNGD DIV5 DUST
     &, CONSCAV_DUST6(ROW_LENGTH, ROWS)                                 &
                                           !COL TOTAL SCVNGD DIV6 DUST

! S Cycle diagnostics
     &,  conwash_nh3(row_length, rows)                                  &
                                          !column total scvngd n in nh3
     &,  conwash_so2(row_length, rows)                                  &
                                          !column total scvngd so2
     &,  conscav_so4ait(row_length, rows)                               &
                                          !column total scvngd so4ait
     &,  conscav_so4acc(row_length, rows)                               &
                                          !column total scvngd so4acc
     &,  conscav_so4dis(row_length, rows)                               &
                                          !column total scvngd so4dis
! Soot cycle diagnostics
     &,  conscav_agedsoot(row_length, rows)                             &
                                            !col total scvngd aged soot
! Biomass aerosol scheme diagnostics
     &,  conscav_agedbmass(row_length, rows)                            & 
                                            !col total scvngd biomass
! Fossil-fuel organic carbon scheme diagnostics
     &,  conscav_agedocff(row_length, rows)
                                            !col total scvngd OCFF

      Real, intent(in) ::                                               &
     &  mid_level(row_length,rows)                                      &
                                       ! true mid-lev
     &, kterm_deep(row_length,rows)                                     &
                                       ! term lev for deep convection
     &, cape(row_length,rows)                                           &
                                       ! CAPE
     &, cape_undilute(row_length,rows)                                  &
                                       ! CAPE from undilute parcel ascent
     &, cin_undilute(row_length,rows)                                   &
                                       ! CIN from undilute parcel ascent
     &, precip_deep(row_length,rows)                                    &
                                       ! deep precipitation
     &, precip_shall(row_length,rows)                                   &
                                       ! shallow precipitation
     &, precip_mid(row_length,rows)                                     &
                                       ! mid-level precipitation
     &, precip_cong(row_length,rows)                                    &
                                       ! congestus precipitation
     &, wstar_up(row_length,rows)                                       &
                                     ! cumulus con vel scale
     &, wstar_dn(row_length,rows)                                       &
                                     ! cumulus con vel scale
     &, mb1(row_length,rows)                                            &
                                ! cloud base mass flux method 1
     &, mb2(row_length,rows)                                            &
                                ! cloud base mass flux method 2
     &, wqt_flux(row_length,rows,wet_model_levels)                      & 
                                                        ! w'qt'
     &, wql_flux(row_length,rows,wet_model_levels)                      & 
                                                        ! w'ql'
     &, wthetal_flux(row_length,rows,wet_model_levels)                  & 
                                                        ! w'thetal'
     &, wthetav_flux(row_length,rows,wet_model_levels)                  & 
                                                        ! w'thetav'
     &, wqt_cb(row_length,rows)                                         & 
                                     ! w'qt' cloud base
     &, wthetal_cb(row_length,rows)                                     & 
                                     ! w'thetal' cloud base
     &, wqt_inv(row_length,rows)                                        & 
                                      ! w'qt' inversion base
     &, wthetal_inv(row_length,rows)                                    & 
                                      ! w'thetal' inversion base
     &, sh_top(row_length,rows)                                         &
                                 ! height of top of shallow convection
     &, sh_base(row_length,rows)                                        &
                                  ! base f shallow convection
     &, ntml_diag(row_length,rows)                                      &
                                       ! top mixing level
     &, ntpar_diag(row_length,rows)                                     &
                                       ! parcel top
     &, freeze_diag(row_length,rows)                                    &
                                       ! freezing lev
     &, cg_term(row_length,rows)                                        &
                                       ! termination level congestus
     &, cg_top(row_length,rows)                                         &
                                       ! top congestus    (m)
     &, cg_base(row_length,rows)                                        &
                                       ! base congestus  (m)
     &, mf_deep(row_length,rows,wet_model_levels)                       &
                                                   ! mass flux deep
     &, mf_congest(row_length,rows,wet_model_levels)                    &
                                                     ! mass flux congest
     &, mf_shall(row_length,rows,wet_model_levels)                      &
                                                    ! mass flux shallow
     &, mf_midlev(row_length,rows,wet_model_levels)                     &
                                                    ! mass flux mid lev
     &, dt_deep(row_length,rows,wet_model_levels)                       &
                                                     ! dT deep
     &, dt_congest(row_length,rows,wet_model_levels)                    &
                                                     ! dT congest
     &, dt_shall(row_length,rows,wet_model_levels)                      &
                                                     ! dT shallow
     &, dt_midlev(row_length,rows,wet_model_levels)                     &
                                                     ! dT mid lev
     &, dq_deep(row_length,rows,wet_model_levels)                       &
                                                     ! dq deep
     &, dq_congest(row_length,rows,wet_model_levels)                    &
                                                     ! dq congest
     &, dq_shall(row_length,rows,wet_model_levels)                      &
                                                     ! dq shallow
     &, dq_midlev(row_length,rows,wet_model_levels)                     &
                                                     ! dq mid lev
     &, du_deep(row_length,rows,wet_model_levels)                       &
                                                     ! du deep
     &, du_congest(row_length,rows,wet_model_levels)                    &
                                                     ! du congest
     &, du_shall(row_length,rows,wet_model_levels)                      &
                                                     ! du shallow
     &, du_midlev(row_length,rows,wet_model_levels)                     &
                                                     ! du mid lev
     &, dv_deep(row_length,rows,wet_model_levels)                       &
                                                     ! dv deep
     &, dv_congest(row_length,rows,wet_model_levels)                    &
                                                     ! dv congest
     &, dv_shall(row_length,rows,wet_model_levels)                      &
                                                     ! dv shallow
     &, dv_midlev(row_length,rows,wet_model_levels)  ! dv mid lev

      Integer                                                           &
     &  ccb(row_length, rows)                                           &
     &, cct(row_length, rows)                                           &
     &, lcbase(row_length, rows)                                        &
     &, lctop(row_length, rows)
      REAL                                                              &
     &  lcca(row_length, rows)

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
     & STASHwork(*)     ! STASH workspace

! Local variables
      Integer                                                           &
     &  i, j, k, ji                                                     &
     &,    icode                                                        &
                                ! Return code  =0 Normal exit  >1 Error
     & ,item                                                            &
                        ! STASH item
     & ,sect            ! STASH section
      Parameter( sect = 5 ) ! for convection

      Character*80  cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_conv')

! max_2d_cca is the maximum aount of 2d cca permitted
      REAL, PARAMETER :: max_2d_cca = 0.5
! min_val is the amount of difference in convective cloud which
! constitites an anvil
      Integer                                                           &
     &  im_index        ! internal model index


      Real                                                              &
     &  interp_data(row_length,rows)                                    &
     &, interp_data_3(row_length,rows,model_levels)                     &
     &, interp_data_3_n(row_length,n_rows,model_levels)                 &
     &, T(row_length, rows, model_levels)                               &
     &, icao_height(row_length,rows)                                    &
     &, scf        ! scaling factor for LS precipitation and scavenging
                   ! diagnostics to allow correct calculation when
                   ! substepping is used
! External routines
      External                                                          &
     &     copydiag, copydiag_3d                                        &
     &  ,Ereport                                                        &
     &, icao_ht


      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output.
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! DIAG.05163 Copy qCL INC: conv inhom to stashwork
! ----------------------------------------------------------------------
      item = 163
! Diag05163_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        qcl_incr_inhom_diag,                                      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
        cmessage="conv_ctl  : error in copydiag_3d(Qcl INC: conv inhom)"
        End if
      End if  ! Diag05163_if1
!
! ----------------------------------------------------------------------
! DIAG.05164 Copy qCF INC: conv inhom to stashwork
! ----------------------------------------------------------------------
      item = 164
! Diag05164_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        qcf_incr_inhom_diag,                                      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
        cmessage="conv_ctl  : error in copydiag_3d(Qcf INC: conv inhom)"
        End if
      End if  ! Diag05164_if1
!
! ----------------------------------------------------------------------
! DIAG.05172 Copy BULK_CF INC: conv inhom to stashwork
! ----------------------------------------------------------------------
      item = 172
! Diag05172_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        bulk_cf_incr_inhom_diag,                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
        cmessage="conv_ctl  : error in copydiag_3d(BCF INC: conv inhom)"
        End if
      End if  ! Diag05172_if1
!
! ----------------------------------------------------------------------
! DIAG.05173 Copy CF_LIQUID INC: conv inhom to stashwork
! ----------------------------------------------------------------------
      item = 173
! Diag05173_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        cf_liquid_incr_inhom_diag,                                &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
        cmessage="conv_ctl  : error in copydiag_3d(CFl INC: conv inhom)"
        End if
      End if  ! Diag05173_if1
!
! ----------------------------------------------------------------------
! DIAG.05174 Copy CF_FROZEN INC: conv inhom to stashwork
! ----------------------------------------------------------------------
      item = 174
! Diag05174_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        cf_frozen_incr_inhom_diag,                                &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
        cmessage="conv_ctl  : error in copydiag_3d(CFf INC: conv inhom)"
        End if
      End if  ! Diag05174_if1
!
! increment diagnostics= modified - previous

      item = 181  ! temperature increment
      If (icode <= 0 .and. sf(item,sect)) Then

      If (.not.allocated(TINC_CVN)) Then
        Allocate ( TINC_CVN(row_length*rows,model_levels) )
      End If

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        T_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 181)"//cmessage
         End if

         Do k = 1,model_levels
           Do j = 1,rows
             Do i = 1,row_length
               ji = (j-1)*row_length+i
               TINC_CVN(ji,k) = T_incr_diagnostic(i,j,k)
             End Do
           End Do
         End Do
 
      End if  !  sf(item,sect)

      item = 182  ! humidity increment
      If (icode <= 0 .and. sf(item,sect)) Then


! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        q_incr_diagnostic,                                        &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 182)"//cmessage
         End if

      End if  !  sf(item,sect)
!
! ----------------------------------------------------------------------
! DIAG.05183 Copy qCL INC: convection to stashwork
! ----------------------------------------------------------------------
      item = 183
! Diag05183_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        qcl_incr_diagnostic,                                      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
        cmessage="conv_ctl  : error in copydiag_3d(Qcl INC: convection)"
        End if
      End if  ! Diag05183_if1
!
! ----------------------------------------------------------------------
! DIAG.05184 Copy qCF INC: convection to stashwork
! ----------------------------------------------------------------------
      item = 184
! Diag05184_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        qcf_incr_diagnostic,                                      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
        cmessage="conv_ctl  : error in copydiag_3d(Qcf INC: convection)"
        End if
      End if  ! Diag05184_if1
!
! ----------------------------------------------------------------------
! DIAG.05185 Copy U Wind INC: convection to stashwork
! ----------------------------------------------------------------------
      item = 185  ! u wind increment
      If (icode <= 0 .and. sf(item,sect)) Then


! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        u_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 185)"//cmessage
         End if

      End if  !  sf(item,sect)

      item = 186  ! v wind increment
      If (icode <= 0 .and. sf(item,sect)) Then


! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        v_incr_diagnostic,                                        &
     &        row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 186)"//cmessage
         End if

      End if  !  sf(item,sect)
! ----------------------------------------------------------------------
! convection source terms in non shallow regions
! ----------------------------------------------------------------------
      item = 187  ! temperature increment
      If (icode <= 0 .and. sf(item,sect)) Then

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              interp_data_3(i,j,k) = T_incr_diagnostic(i,j,k)           &
                                      *(1.0-shallowc(i,j))
            End Do
          End Do
        End Do

         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              interp_data_3,                                            &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)

         If (icode > 0) Then
            cmessage=": error in copydiag_3d(item 187)"//cmessage
         End if

      Endif  !  sf(item,sect)

      item = 188  ! humidity increment
      If (icode <= 0 .and. sf(item,sect)) Then

        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              interp_data_3(i,j,k) = q_incr_diagnostic(i,j,k)           &
                                      *(1.0-shallowc(i,j))
            End Do
          End Do
        End Do

         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
              interp_data_3,                                            &
              row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)

         If (icode > 0) Then
            cmessage=": error in copydiag_3d(item 188)"//cmessage
         End if

      Endif  !  sf(item,sect)
!
! ----------------------------------------------------------------------
! DIAG.05192 Copy BULK_CF INC: convection to stashwork
! ----------------------------------------------------------------------
      item = 192
! Diag05192_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        bulk_cf_incr_diagnostic,                                  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
        cmessage="conv_ctl  : error in copydiag_3d(BCF INC: convection)"
        End if
      End if  ! Diag05192_if1
!
! ----------------------------------------------------------------------
! DIAG.05193 Copy CF_LIQUID INC: convection to stashwork
! ----------------------------------------------------------------------
      item = 193
! Diag05193_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        cf_liquid_incr_diagnostic,                                &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
        cmessage="conv_ctl  : error in copydiag_3d(CFl INC: convection)"
        End if
      End if  ! Diag05193_if1
!
! ----------------------------------------------------------------------
! DIAG.05194 Copy CF_FROZEN INC: convection to stashwork
! ----------------------------------------------------------------------
      item = 194
! Diag05194_if1:
      If (icode  <=  0  .AND.  sf(item,sect)) Then
!
!
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        cf_frozen_incr_diagnostic,                                &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
        cmessage="conv_ctl  : error in copydiag_3d(CFf INC: convection)"
        End if
      End if  ! Diag05194_if1
!
! ----------------------------------------------------------------------
! Section  Convective Rain
! ----------------------------------------------------------------------
! Item 201 Convective rainfall,resolve to accumulate over timestep

      item = 201 ! convective rainfall accumulation
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(stashwork(si(201,5,im_index)),conv_rain,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,201,                                           &
     &        icode,cmessage)

! scaling factor to convert accum
! over model timestep as the following timestep may be a phys2 substep
         scf = 1.*Num_Substeps


         Do i=1,row_length*rows
            stashwork(si(201,5,im_index)+i-1)=                          &
     &           stashwork(si(201,5,im_index)+i-1)*scf*timestep
         End do

      End if


! ----------------------------------------------------------------------
! Section  Convective Snow
! ----------------------------------------------------------------------
! Item 202 Convective Snowfall, resolve to accumulate overtimestep.

      item = 202 ! convective snowfall accumulation
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag
         Call copydiag(stashwork(si(202,5,im_index)),conv_snow,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,202,                                           &
     &        icode,cmessage)


! scaling factor to convert accum
! over model timestep as the following timestep may be a phys2 substep
         scf = 1.*Num_Substeps

         Do i=1,row_length*rows
            stashwork(si(202,5,im_index)+i-1)=                          &
     &           stashwork(si(202,5,im_index)+i-1)*scf*timestep
         End do

      End if

! ----------------------------------------------------------------------
! Section  Convective rainfall rates
! ----------------------------------------------------------------------
! C Item 205 Convective rainfall rates

      item = 205 ! convective rainfall rate
      If (icode <= 0 .and. sf(item,sect)) Then

      If (.not.allocated(CVRR)) Then
        Allocate ( CVRR(row_length*rows) )
      End If

! DEPENDS ON: copydiag
         Call copydiag(stashwork(si(205,5,im_index)),CONV_RAIN,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,205,                                           &
     &        icode,cmessage)

        Do j = 1,rows
          Do i = 1,row_length
            ji = (j-1)*row_length+i
            CVRR(ji) = CONV_RAIN(i,j)
          End Do
        End Do

      End if


! ----------------------------------------------------------------------
! Section  206 Convective snowfall rates
! ----------------------------------------------------------------------
! C Item 206 Convective snowfall rates

      item = 206 ! convective snowfall rate
      If (icode <= 0 .and. sf(item,sect)) Then

      If (.not.allocated(CVSR)) Then
        Allocate ( CVSR(row_length*rows) )
      End If

! DEPENDS ON: copydiag
         Call copydiag(stashwork(si(206,5,im_index)),CONV_SNOW,         &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,206,                                           &
     &        icode,cmessage)

        Do j = 1,rows
          Do i = 1,row_length
            ji = (j-1)*row_length+i
            CVSR(ji) = CONV_SNOW(i,j)
          End Do
        End Do

      End if


! ----------------------------------------------------------------------
! Section Convective Cloud Amount
! ----------------------------------------------------------------------
! Item 212

      item = 212 ! convective cloud amount
      If (icode <= 0 .and. sf(item,sect)) Then

         If (lcv_3d_cca) Then
            If (.not.allocated(CONVCC)) Then
              Allocate ( CONVCC(row_length*rows,n_cca_levels) )
            End If
            Do k = 1, n_cca_levels
               Do j = 1, rows
                  Do i = 1, row_length
                     ji = (j-1)*row_length+i
                     interp_data_3(i,j,k) = cca(i,j,k)
                     CONVCC(ji,k) = interp_data_3(i,j,k)
                  End Do
               End Do
            End Do
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (stashwork(si(item,sect,im_index)),        &
     &           interp_data_3,                                         &
     &           row_length,rows,n_cca_levels,0,0,0,0, at_extremity,    &
     &           stlist(1,stindex(1,item,sect,im_index)),len_stlist,    &
     &           stash_levels,num_stash_levels+1,                       &
     &           atmos_im,sect,item,                                    &
     &           icode,cmessage)
         Else
            If (.not.allocated(CONVCC)) Then
              Allocate ( CONVCC(row_length*rows,wet_model_levels) )
            End If
            Do k = 1, wet_model_levels
               Do j = 1, rows
                  Do i = 1, row_length
                     If (k  >=  ccb(i,j) .and.                          &
     &                   k  <   cct(i,j)) Then
                        interp_data_3(i,j,k) = cca(i,j,1)
                     Else
                        interp_data_3(i,j,k) = 0.0
                     End If
                     ji = (j-1)*row_length+i
                     CONVCC(ji,k) = interp_data_3(i,j,k)
                  End Do
               End Do
            End Do
! DEPENDS ON: copydiag_3d
            Call copydiag_3d (stashwork(si(item,sect,im_index)),        &
     &           interp_data_3,                                         &
     &           row_length,rows,wet_model_levels,0,0,0,0, at_extremity,&
     &           stlist(1,stindex(1,item,sect,im_index)),len_stlist,    &
     &           stash_levels,num_stash_levels+1,                       &
     &           atmos_im,sect,item,                                    &
     &           icode,cmessage)
         End If


         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 212)"//cmessage
         End if

      End if

! ----------------------------------------------------------------------
! Section Convective Cloud Liquid Water
! ----------------------------------------------------------------------
! Item 213

      item = 213 ! convective cloud liquid water
      If (icode <= 0 .and. sf(item,sect)) Then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d (stashwork(si(item,sect,im_index)), ccw,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

         If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 213)"//cmessage
         End if

      End if

! ----------------------------------------------------------------------
! Section Total precipitation
! ----------------------------------------------------------------------
! Item 226

      item = 226 ! total precipitation accumulation
      If (icode <= 0 .and. sf(item,sect)) Then


! scaling factor to convert accum
! over model timestep as the following timestep may be a phys2 substep
         scf = 1.*Num_Substeps

         Do j = 1, rows
            Do i = 1, row_length

               interp_data(i,j) = ( ls_rain(i,j) + ls_snow(i,j)         &
     &                              + conv_rain(i,j) + conv_snow(i,j) ) &
     &                          * scf * timestep
            End Do
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(226,5,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,226,                                           &
     &        icode,cmessage)

      End if

! ----------------------------------------------------------------------
! Section Total Rain rate
! ----------------------------------------------------------------------
! Item 214 ls_rain +conv_rain

      item = 214 ! total rain rate
      If (icode <= 0 .and. sf(item,sect)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = (ls_rain(i,j) +conv_rain(i,j) )
            End Do
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(214,5,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,214,                                           &
     &        icode,cmessage)

      End if


! ----------------------------------------------------------------------
! Section Total Snow rate
! ----------------------------------------------------------------------
! Item 215 ls_snow + conv_snow

      item = 215 ! total snow rate
      If (icode <= 0 .and. sf(item,sect)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = (ls_snow(i,j) + conv_snow(i,j) )
            End Do
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(215,5,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,215,                                           &
     &        icode,cmessage)

      End if

! ----------------------------------------------------------------------
! Section Total Precipitation rate
! ----------------------------------------------------------------------
! Item 216 ls_rain + conv_rain + ls_snow + conv_snow

      item = 216 ! total precipitation rate
      If (icode <= 0 .and. sf(item,sect)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = (ls_rain(i,j) + conv_rain(i,j) +      &
     &                             ls_snow(i,j) + conv_snow(i,j))
            End Do
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(216,5,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,216,                                           &
     &        icode,cmessage)

      End if

! ----------------------------------------------------------------------
! Calculate pressure at convective cloud base
! ----------------------------------------------------------------------

      IF (icode <= 0 .and. (sf(207,sect) .or. sf(210,sect))) THEN
         DO j = 1, rows
            DO i = 1, row_length
               IF (ccb(i,j)  /=  0 )THEN
                 interp_data(i,j) = p(i,j,ccb(i,j))
               ELSE
                 interp_data(i,j) =  RMDI
               END IF
            END DO
         END DO
       END IF

! ----------------------------------------------------------------------
! Section pressure at convective cloud base
! ----------------------------------------------------------------------
! Item 207 Pressure at convective cloud base

      item = 207 ! Pressure at convective cloud base
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(207,5,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,207,                                           &
     &        icode,cmessage)

      END IF

! ----------------------------------------------------------------------
! Calculate height at convective cloud base
!  (ICAO standard atmosphere heights)
! ----------------------------------------------------------------------

      IF (icode <= 0 .and. sf(210,sect)) THEN
! DEPENDS ON: icao_ht
         CALL Icao_ht(interp_data,row_length,rows,icao_height)
      END IF

! ----------------------------------------------------------------------
! Section ICAO height at convective cloud base
! ----------------------------------------------------------------------
! Item 210 ICAO height at convective cloud base

      item = 210 ! ICAO height at convective cloud base
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(210,5,im_index)),icao_height,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,210,                                           &
     &        icode,cmessage)

      END IF

! ----------------------------------------------------------------------
! Calculate pressure at convective cloud top
! ----------------------------------------------------------------------

      IF (icode <= 0 .and. (sf(208,sect) .or. sf(211,sect))) THEN
         DO j = 1, rows
            DO i = 1, row_length
               IF (cct(i,j)  /=  0 )THEN
                 interp_data(i,j) = p(i,j,cct(i,j))
               ELSE
                 interp_data(i,j) =  RMDI
               END IF
            END DO
         END DO
       END IF

! ----------------------------------------------------------------------
! Section pressure at convective cloud top
! ----------------------------------------------------------------------
! Item 208 Pressure at convective cloud top

      item = 208 ! Pressure at convective cloud top
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(208,5,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,208,                                           &
     &        icode,cmessage)

      END IF

! ----------------------------------------------------------------------
! Calculate height at convective cloud top
!  (ICAO standard atmosphere heights)
! ----------------------------------------------------------------------

      IF (icode <= 0 .and. sf(211,sect)) THEN
! DEPENDS ON: icao_ht
         CALL Icao_ht(interp_data,row_length,rows,icao_height)
      END IF

! ----------------------------------------------------------------------
! Section ICAO height at convective cloud top
! ----------------------------------------------------------------------
! Item 210 ICAO height at convective cloud top

      item = 211 ! ICAO height at convective cloud top
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(211,5,im_index)),icao_height,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,211,                                           &
     &        icode,cmessage)

      END IF

! ----------------------------------------------------------------------
!  CAPE
! ----------------------------------------------------------------------
! Item 217 CAPE

      item=217
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(217,5,im_index)),cape,               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,217,                                            &
     &       ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Section Lowest Convective Cloud Base
! ----------------------------------------------------------------------
! Item 218 Lowest convective cloud base

      item = 218 ! Lowest convective cloud base
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(218,5,im_index)),lcbase,            &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,218,                                           &
     &        icode,cmessage)

      END IF

! ----------------------------------------------------------------------
! Section Lowest Convective Cloud Top
! ----------------------------------------------------------------------
! Item 219 Lowest convective cloud top

      item = 219 ! Lowest convective cloud top
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(219,5,im_index)),lctop,             &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,219,                                           &
     &        icode,cmessage)

      END IF

! ----------------------------------------------------------------------
! Section Lowest Convective Cloud Amount
! ----------------------------------------------------------------------
! Item 220 Lowest convective cloud amount

      item = 220 ! Lowest convective cloud amount
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(220,5,im_index)),lcca,              &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,220,                                           &
     &        icode,cmessage)

      END IF

! ----------------------------------------------------------------------
! Calculate pressure at lowest convective cloud base
! ----------------------------------------------------------------------

      IF (icode <= 0 .and. (sf(222,sect) .or. sf(224,sect))) THEN
         DO j = 1, rows
            DO i = 1, row_length
              IF (lcbase(i,j)  /=  0)THEN
                interp_data(i,j) = p(i,j,lcbase(i,j))
              ELSE
                interp_data(i,j) = RMDI
              END IF
            END DO
         END DO
       END IF

! ----------------------------------------------------------------------
! Section pressure at lowest convective cloud base
! ----------------------------------------------------------------------
! Item 222 Pressure at lowest convective cloud base

      item = 222 ! Pressure at lowest convective cloud base
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(222,5,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,222,                                           &
     &        icode,cmessage)

      END IF

! ----------------------------------------------------------------------
! Calculate height at convective cloud base
!  (ICAO standard atmosphere heights)
! ----------------------------------------------------------------------

      IF (icode <= 0 .and. sf(224,sect)) THEN
! DEPENDS ON: icao_ht
         CALL Icao_ht(interp_data,row_length,rows,icao_height)
      END IF

! ----------------------------------------------------------------------
! Section ICAO height at lowest convective cloud base
! ----------------------------------------------------------------------
! Item 224 ICAO height at convective cloud base

      item = 224 ! ICAO height at convective cloud base
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(224,5,im_index)),icao_height,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,224,                                           &
     &        icode,cmessage)

      END IF

! ----------------------------------------------------------------------
! Calculate pressure at convective cloud top
! ----------------------------------------------------------------------

      IF (icode <= 0 .and. (sf(223,sect) .or. sf(225,sect))) THEN
         DO j = 1, rows
            DO i = 1, row_length
              IF (lctop(i,j)  /=  0)THEN
                interp_data(i,j) = p(i,j,lctop(i,j))
              ELSE
                interp_data(i,j) = RMDI
              END IF
            END DO
         END DO
       END IF

! ----------------------------------------------------------------------
! Section pressure at lowest convective cloud top
! ----------------------------------------------------------------------
! Item 223 Pressure at lowest convective cloud top

      item = 223 ! Pressure at lowest convective cloud top
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(223,5,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,223,                                           &
     &        icode,cmessage)

      END IF

! ----------------------------------------------------------------------
! Calculate ICAO height at convective cloud top
! ----------------------------------------------------------------------

      IF (icode <= 0 .and. sf(225,sect)) THEN
! DEPENDS ON: icao_ht
         CALL Icao_ht(interp_data,row_length,rows,icao_height)
      END IF

! ----------------------------------------------------------------------
! Section ICAO height at lowest convective cloud top
! ----------------------------------------------------------------------
! Item 225 ICAO height at convective cloud top

      item = 225 ! ICAO height at convective cloud top
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(225,5,im_index)),icao_height,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,225,                                           &
     &        icode,cmessage)

      END IF

! ----------------------------------------------------------------------
! Calculate 2D convective cloud amount
! ----------------------------------------------------------------------

      IF (icode <= 0 .and. sf(262,sect)) THEN

        DO j = 1, rows
          DO i = 1, row_length
           IF( lcv_3d_cca)THEN
! Convert 3d cca to 2d cca
            IF( ccb(i,j)  ==  0 .AND. cct(i,j)  ==  0 )THEN
              interp_data(i,j) = 0.0
            ELSE IF ( cca(i,j,cct(i,j)-1)  >   cca(i,j,ccb(i,j))) THEN
            If (tower_factor  >   0.0) then
              ! Back out the original 2D convective cloud amount
              ! before the anvil parametrization was applied
              interp_data(i,j) = cca(i,j,ccb(i,j))  / tower_factor
            Else
              ! Tower_factor is zero (likely to be a PC2 run)
              ! and hence assume that no convective cloud is
              ! required.
              interp_data(i,j) = 0.0
            End if  ! tower_factor  >   0.0
            ELSE
              interp_data(i,j)=cca(i,j,ccb(i,j))
            END IF
!   Ensure that the maximum 2d cca isn't breached
            IF ( interp_data(i,j)  >   max_2d_cca )THEN
              interp_data(i,j) = max_2d_cca
            END IF
           ELSE
! Copy 2d cca
            interp_data(i,j) = cca(i,j,1)
           END IF
          END DO
        END DO
      END IF

!-----------------------------------------------------------------------
! Section 3D Convective rainfall rate
! Item 227
! ---------------------------------------------------------------------- 

      item = 227                         ! 3D convective rainfall rate 

      If (icode <= 0 .and. sf(item,sect)) Then 
 
! DEPENDS ON: copydiag_3d 
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        conv_rain_3d,                                             & 
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   & 
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       & 
     &        stash_levels,num_stash_levels+1,                          & 
     &        atmos_im,sect,item,                                       & 
     &        icode,cmessage) 
         
         If (icode >  0) Then 
           cmessage=": error in copydiag_3d(item 227)"//cmessage 
         End if 
 
      End if 

! ---------------------------------------------------------------------- 
! Section 3D Convective snowfall rate
! Item 228
! ---------------------------------------------------------------------- 
                         
      item = 228                        ! 3D convective snowfall rate 

      If (icode <= 0 .and. sf(item,sect)) Then 
         
! DEPENDS ON: copydiag_3d 
         Call copydiag_3d (stashwork(si(item,sect,im_index)),           &
     &        conv_snow_3d,                                             & 
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   & 
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       & 
     &        stash_levels,num_stash_levels+1,                          & 
     &        atmos_im,sect,item,                                       & 
     &        icode,cmessage) 
         
         If (icode >  0) Then 
            cmessage=": error in copydiag_3d(item 228)"//cmessage 
         End if 

       End if 

! ----------------------------------------------------------------------
! Section 2D convective cloud amount
! ----------------------------------------------------------------------
! Item 262 2D convective cloud amount

      item = 262 ! 2D convective cloud amount
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(262,5,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,5,262,                                           &
     &        icode,cmessage)
      END IF

! ----------------------------------------------------------------------
! Section 1.2  Tendency diagnostics
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section Temperature at end of timestep.
! ----------------------------------------------------------------------
! Item 209 Temperature

      item = 209 ! temperature
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1, model_levels
            Do j = 1, rows
               Do i = 1, row_length
                  T(i,j,k) = (theta(i,j,k) + theta_inc(i,j,k) )*        &
     &                 exner_theta_levels(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(209,5,im_index)),T,              &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,209,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,209,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag_3d(item 209)"//cmessage
         End if

      End if


! ----------------------------------------------------------------------
! Section q at end of timestep.
! ----------------------------------------------------------------------
! Item 010 q

      item =  10 ! temperature
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1, wet_model_levels
            Do j = 1, rows
               Do i = 1, row_length
                  interp_data_3(i,j,k) = q(i,j,k) + q_inc(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(010,5,im_index)),interp_data_3,  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,010,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,010,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag_3d(item 010)"//cmessage
         End if

      End if


! ----------------------------------------------------------------------
!  233 undilute CAPE
! ----------------------------------------------------------------------
! Item 233 undilute CAPE

      item=233
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),cape_undilute,     &
             row_length,rows,0,0,0,0, at_extremity,                     &
             atmos_im,5,item,                                           &
             ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
!  234 undilute CIN
! ----------------------------------------------------------------------
! Item 234 undilute CIN

      item=234
      IF (icode <= 0 .and. sf(item,sect)) THEN

       ! Change sign of CIN to be positive.
        Do j = 1, rows
          Do i = 1, row_length
            interp_data(i,j) = -1.0*cin_undilute(i,j)
          End Do
        End Do

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),interp_data,       &
             row_length,rows,0,0,0,0, at_extremity,                     &
             atmos_im,5,item,                                           &
             ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Section u latest.
! ----------------------------------------------------------------------
! Item 235 u

      item = 235 ! u
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1, model_levels
            Do j = 1, rows
               Do i = 1, row_length
                  interp_data_3(i,j,k) = u(i,j,k) + R_u(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(235,5,im_index)),interp_data_3,  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,235,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,235,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag_3d(item 235)"//cmessage
         End if

      End if


! ----------------------------------------------------------------------
! Section v latest
! ----------------------------------------------------------------------
! Item 236 v

      item = 236 ! u
      If (icode <= 0 .and. sf(item,sect)) Then

         Do k = 1, model_levels
            Do j = 1, n_rows
               Do i = 1, row_length
                  interp_data_3_n(i,j,k) = v(i,j,k) + R_v(i,j,k)
               End Do
            End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(236,5,im_index)),interp_data_3_n,&
     &        row_length,n_rows,model_levels,0,0,0,0, at_extremity,     &
     &        stlist(1,stindex(1,236,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,236,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage=": error in copydiag_3d(item 236)"//cmessage
         End if

      End if
! ----------------------------------------------------------------------
! deep convec. massflux components (for 4D-Var)
! ----------------------------------------------------------------------
! item 246 rho-level massflux in a) no shallow regions and 
!                                b) predominantly entraining regions 
!                                   (mass flux is growing with heigth)

      item=246 ! component of half level massflux
      If (icode == 0 .and. sf(item,sect)) Then

        Do k = 1, model_levels-1
          Do j = 1, rows
            Do i = 1, row_length
              If (up_flux_half(i,j,k).lt.up_flux_half(i,j,k+1)) Then
                interp_data_3(i,j,k) = up_flux_half(i,j,k)              &
                                      *(1.0-shallowc(i,j))
              Else
                interp_data_3(i,j,k) =0.0
              Endif
            End Do
          End Do
        End Do
        Do j = 1, rows
          Do i = 1, row_length
            interp_data_3(i,j,model_levels) = 0.0
          End Do
        End Do

         Call copydiag_3d(STASHwork(si(item,5,im_index)),               &
              interp_data_3,                                            &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,5,item,                                          &
              icode,cmessage)

         If (icode >  0) Then
            cmessage="conv_ctl  :                                       &
     & error in copydiag_3d(mass flux component 246)"
            goto 9999
         Endif

      End if


! item 249 updraught mass flux on half levels

      item=249 ! updraught mass flux on half levels
      If (icode == 0 .and. sf(item,sect)) then


! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(249,5,im_index)),up_flux_half,   &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,249,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,249,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="conv_ctl  :                                       &
     & error in copydiag_3d(up mass flux half levs)"
            go to 9999
         End if

      End if

! item 250 updraught mass flux

      item=250 ! updraught mass flux
      If (icode == 0 .and. sf(item,sect)) then


! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(250,5,im_index)),up_flux,        &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,250,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,250,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="conv_ctl  : error in copydiag_3d(up mass flux)"
            go to 9999
         End if

      End if

! item 251 downdraught mass flux

      item=251
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(251,5,im_index)),dwn_flux,       &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,251,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,251,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="conv_ctl  : error in copydiag_3d(dwn mass flux)"
            go to 9999
         End if

      End if

! item 252 updraught entrainment rate

      item=252
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(252,5,im_index)),entrain_up,     &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,252,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,252,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="conv_ctl  : error in copydiag_3d(up entrainment)"
            go to 9999
         End if

      End if

! item 253 updraught detrainment rate

      item=253
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(253,5,im_index)),detrain_up,     &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,253,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,253,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
            cmessage="conv_ctl  : error in copydiag_3d(up detrainment)"
            go to 9999
         End if

      End if

! item 254 downdraught entrainment rate

      item=254
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(254,5,im_index)),entrain_dwn,    &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,254,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,254,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dwn entrainment)"
          go to 9999
         End if

      End if

! item 255 downdraught detrainment rate

      item=255
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(255,5,im_index)),detrain_dwn,    &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,255,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,255,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dwn detrainment)"
            go to 9999
         End if

      End if

! item 256 time rate of change of u on the p-grid

      item=256
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(256,5,im_index)),dubydt_p,       &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,256,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,256,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dudt (p grid)"
            go to 9999
         End if

      End if

! item 257 time rate of change of v on the p-grid

      item=257
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(257,5,im_index)),dvbydt_p,       &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,257,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,257,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dvdt (p grid)"
            go to 9999
         End if
      End if

! item 258 u-stress for deep convection (p-grid)

      item=258
      If (icode == 0.and.sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(258,5,im_index)),uw_dp,          &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,258,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,258,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(uw deep)"
            go to 9999
         End if
      End if

! item 259 v-stress for deep convection (p-grid)

      item=259
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(259,5,im_index)),vw_dp,          &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,259,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,259,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(vw deep)"
            go to 9999
         End if
      End if

! item 260 u stress for shallow convection (p-grid)

      item=260
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(260,5,im_index)),uw_shall,       &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,260,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,260,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(uw shallow)"
            go to 9999
         End if
      End if

! item 261 v-stress for shallow convection (p-grid)

      item=261
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(261,5,im_index)),vw_shall,       &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,261,5,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,261,                                           &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(vw shallow)"
            go to 9999
         End if

      End if

! ----------------------------------------------------------------------
!  Shallow convection indicator
! ----------------------------------------------------------------------
! Item 270 Shallow convection diagnostic

      item=270
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(270,5,im_index)),SHALLOWC,           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,270,                                            &
     &       ICODE,CMESSAGE)

      END IF

! ----------------------------------------------------------------------
!  Convection over orography indicator
! ----------------------------------------------------------------------
! Item 271 Convection over orography indicator

      item=271
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(271,5,im_index)),cu_over_orog,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,271,                                            &
     &       ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
!  Item 272 - indicator for Mid level convection
! ----------------------------------------------------------------------

      item=272
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(272,5,im_index)),mid_level,          &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,272,                                            &
     &       ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 273 - NTML level number for top of mixing level
! ----------------------------------------------------------------------

      item=273
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(273,5,im_index)),ntml_diag,          &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,273,                                            &
     &       ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 274 -  NTPAR model level number for top of initial parcel ascent
! ----------------------------------------------------------------------

      item=274
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(274,5,im_index)),ntpar_diag,         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,                                           &
     &       ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
!  Item 275 - Freeze_lev
! ----------------------------------------------------------------------

      item=275
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(275,5,im_index)),freeze_diag,        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,275,                                            &
     &       ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 276 - kterm deep convection
! ----------------------------------------------------------------------

      item=276
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(276,5,im_index)),kterm_deep,         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,276,                                            &
     &       ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 277 -  Total precipitation from deep convection
! ----------------------------------------------------------------------

      item=277
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(277,5,im_index)),precip_deep,        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,277,                                            &
     &       ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 278 -  Total precipitation from shallow convection
! ----------------------------------------------------------------------

      item=278
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(278,5,im_index)),precip_shall,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,278,                                            &
     &       ICODE,CMESSAGE)

      END IF

! ----------------------------------------------------------------------
! Item 279 -  Total precipitation from mid_level convection
! ----------------------------------------------------------------------

      item=279
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(279,5,im_index)),precip_mid,         &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,279,                                            &
     &       ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 280 -  Total precipitation from congestus convection
! ----------------------------------------------------------------------

      item=280
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),precip_cong,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,                                           &
     &       ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
!  Diagnostics only available from turbulence convection scheme
! ----------------------------------------------------------------------
! Item 290 -  w'qt' flux from turbulent convection
! ----------------------------------------------------------------------

      item=290
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),wqt_flux,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(wqt_flux)"
            go to 9999
         End if

      End if
! ----------------------------------------------------------------------
! Item 291 -  w'ql' flux from turbulent convection
! ----------------------------------------------------------------------

      item=291
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),wql_flux,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(wql_flux)"
            go to 9999
         End if

      End if
! ----------------------------------------------------------------------
! Item 292 -  w'thetal' flux from turbulent convection
! ----------------------------------------------------------------------

      item=292
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),wthetal_flux,  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(wthetal_flux)"
            go to 9999
         End if

      End if
! ----------------------------------------------------------------------
! Item 293 -  w'thetav' flux from turbulent convection
! ----------------------------------------------------------------------

      item=293
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),wthetav_flux,  &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(wthetav_flux)"
            go to 9999
         End if

      End if
! ----------------------------------------------------------------------
! Item 300 -  subcloud layer convective velocity scale
! ----------------------------------------------------------------------

      item=300
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),wstar_dn,          &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 301 -  cumulus layer convective velocity scale
! ----------------------------------------------------------------------

      item=301
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),wstar_up,          &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF

! ----------------------------------------------------------------------
! Item 302 -   cloud base mass flux 1
! ----------------------------------------------------------------------

      item=302
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),mb1,               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 303 -   cloud base mass flux 2
! ----------------------------------------------------------------------

      item=303
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),mb2,               &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 304 - wqt at cloud base
! ----------------------------------------------------------------------

      item=304
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),wqt_cb,            &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 305 - wthetal at cloud base
! ----------------------------------------------------------------------

      item=305
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),wthetal_cb,        &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 306 - wqt at inversion
! ----------------------------------------------------------------------

      item=306
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),wqt_inv,           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 307 - wthetal at inversion
! ----------------------------------------------------------------------

      item=307
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),wthetal_inv,       &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 308 - height of top of shallow convection
!            Diagnostic = height * 1 if shallow convection otherwise
!            zero
! ----------------------------------------------------------------------

      item=308
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),sh_top,            &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
! Item 309 - height of base of shallow convection
!            Diagnostic = height * 1 if shallow convection otherwise
!            zero
! ----------------------------------------------------------------------

      item=309
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),sh_base,           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
!  Item 310 - Congestus convection indicator
! ----------------------------------------------------------------------

      item=310
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),congestus_ind,     &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF
! ----------------------------------------------------------------------
!  Item 311 - Congestus convection indicator  2
! ----------------------------------------------------------------------

      item=311
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),congestus_ind2,    &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF

! ----------------------------------------------------------------------
!  Item 312 - Congestus termination level
! ----------------------------------------------------------------------

      item=312
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),cg_term,           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF

! ----------------------------------------------------------------------
!  Item 313 - Congestus top height
! ----------------------------------------------------------------------

      item=313
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),cg_top,            &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF

! ----------------------------------------------------------------------
!  Item 314 - Congestus base height
! ----------------------------------------------------------------------

      item=314
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag
        CALL COPYDIAG(STASHWORK(SI(item,5,im_index)),cg_base,           &
     &       row_length,rows,0,0,0,0, at_extremity,                     &
     &       atmos_im,5,item,ICODE,CMESSAGE)

      END IF

! ----------------------------------------------------------------------
! Item 320 - mass flux deep
! ----------------------------------------------------------------------

      item=320
      If (icode == 0 .and. sf(item,sect)) then

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),mf_deep,       &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(mf_deep)"
            go to 9999
         End if

      End if
! ----------------------------------------------------------------------
! Item 321 - mass flux congestus
! ----------------------------------------------------------------------

      item=321
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),mf_congest,    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(mf_congest)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 322 - mass flux shallow
! ----------------------------------------------------------------------

      item=322
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),mf_shall,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(mf_shall)"
            go to 9999
         End if
      END IF

! ----------------------------------------------------------------------
! Item 323 - mass flux mid -level
! ----------------------------------------------------------------------

      item=323
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),mf_midlev,     &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(mf_midlev)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 324 - DT deep
! ----------------------------------------------------------------------

      item=324
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dt_deep,       &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dt_deep)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 325 - dt congestus
! ----------------------------------------------------------------------

      item=325
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dt_congest,    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dt_congest)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 326 - dT shallow
! ----------------------------------------------------------------------

      item=326
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dt_shall,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dt_shall)"
            go to 9999
         End if
      END IF

! ----------------------------------------------------------------------
! Item 327 - dT mid -level
! ----------------------------------------------------------------------

      item=327
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dt_midlev,     &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dt_midlev)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 328 - dq deep
! ----------------------------------------------------------------------

      item=328
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dq_deep,       &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dq_deep)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 329 - dq congestus
! ----------------------------------------------------------------------

      item=329
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dq_congest,    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dq_congest)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 330 - dq shallow
! ----------------------------------------------------------------------

      item=330
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dq_shall,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dq_shall)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 331 - dq mid -level
! ----------------------------------------------------------------------

      item=331
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dq_midlev,     &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dq_midlev)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 332 - du deep
! ----------------------------------------------------------------------

      item=332
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),du_deep,       &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(du_deep)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 333 - du congestus
! ----------------------------------------------------------------------

      item=333
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),du_congest,    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(du_congest)"
            go to 9999
         End if
      END IF

! ----------------------------------------------------------------------
! Item 334 - du shallow
! ----------------------------------------------------------------------

      item=334
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),du_shall,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(du_shall)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 335 - du mid -level
! ----------------------------------------------------------------------

      item=335
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),du_midlev,     &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(du_midlev)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 336 - dv deep
! ----------------------------------------------------------------------

      item=336
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dv_deep,       &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dv_deep)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 337 - dv congestus
! ----------------------------------------------------------------------

      item=337
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dv_congest,    &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dv_congest)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 338 - dv shallow
! ----------------------------------------------------------------------

      item=338
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dv_shall,      &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dv_shall)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Item 339 - dv mid -level
! ----------------------------------------------------------------------

      item=339
      IF (icode <= 0 .and. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(item,5,im_index)),dv_midlev,     &
     &        row_length,rows,wet_model_levels,0,0,0,0, at_extremity,   &
     &        stlist(1,stindex(1,item,5,im_index)),len_stlist,          &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,5,item,                                          &
     &        icode,cmessage)

         If (icode  >   0) Then
          cmessage="conv_ctl  : error in copydiag_3d(dv_midlev)"
            go to 9999
         End if

      END IF

! ----------------------------------------------------------------------
! Section Mineral Dust scavenging diagnostics
! ----------------------------------------------------------------------
!
      IF (L_DUST) THEN

        scf = 1.*Num_Substeps
        ITEM = 281
        IF (ICODE <= 0 .AND. SF(ITEM,SECT)) THEN
!
!  Convert scavenged dust per timestep to flux per sec
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              CONSCAV_DUST1(I,J) = CONSCAV_DUST1(I,J)/(scf*TIMESTEP)
            END DO
          END DO
!
! DEPENDS ON: copydiag
          CALL COPYDIAG(STASHWORK(SI(ITEM,SECT,IM_INDEX)),CONSCAV_DUST1,&
     &        ROW_LENGTH,ROWS,0,0,0,0,AT_EXTREMITY,                     &
     &        ATMOS_IM,SECT,ITEM,                                       &
     &        ICODE,CMESSAGE)
!
          IF (ICODE >  0) THEN
            CMESSAGE=": ERROR IN COPYDIAG(ITEM 281)"//CMESSAGE
          END IF
!
        END IF
!
!
        ITEM = 282
        IF (ICODE <= 0 .AND. SF(ITEM,SECT)) THEN
!
!  Convert scavenged dust per timestep to flux per sec
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              CONSCAV_DUST2(I,J) = CONSCAV_DUST2(I,J)/(scf*TIMESTEP)
            END DO
          END DO
!
! DEPENDS ON: copydiag
          CALL COPYDIAG(STASHWORK(SI(ITEM,SECT,IM_INDEX)),CONSCAV_DUST2,&
     &        ROW_LENGTH,ROWS,0,0,0,0,AT_EXTREMITY,                     &
     &        ATMOS_IM,SECT,ITEM,                                       &
     &        ICODE,CMESSAGE)
!
          IF (ICODE >  0) THEN
            CMESSAGE=": ERROR IN COPYDIAG(ITEM 282)"//CMESSAGE
          END IF
!
        END IF
!
!
        ITEM = 283
        IF (ICODE <= 0 .AND. SF(ITEM,SECT)) THEN
!
!  Convert scavenged dust per timestep to flux per sec
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              CONSCAV_DUST3(I,J) = CONSCAV_DUST3(I,J)/(scf*TIMESTEP)
            END DO
          END DO
!
! DEPENDS ON: copydiag
          CALL COPYDIAG(STASHWORK(SI(ITEM,SECT,IM_INDEX)),CONSCAV_DUST3,&
     &        ROW_LENGTH,ROWS,0,0,0,0,AT_EXTREMITY,                     &
     &        ATMOS_IM,SECT,ITEM,                                       &
     &        ICODE,CMESSAGE)
!
          IF (ICODE >  0) THEN
            CMESSAGE=": ERROR IN COPYDIAG(ITEM 283)"//CMESSAGE
          END IF
!
        END IF
!
!
        ITEM = 284
        IF (ICODE <= 0 .AND. SF(ITEM,SECT)) THEN
!
!  Convert scavenged dust per timestep to flux per sec
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              CONSCAV_DUST4(I,J) = CONSCAV_DUST4(I,J)/(scf*TIMESTEP)
            END DO
          END DO
!
! DEPENDS ON: copydiag
          CALL COPYDIAG(STASHWORK(SI(ITEM,SECT,IM_INDEX)),CONSCAV_DUST4,&
     &        ROW_LENGTH,ROWS,0,0,0,0,AT_EXTREMITY,                     &
     &        ATMOS_IM,SECT,ITEM,                                       &
     &        ICODE,CMESSAGE)
!
          IF (ICODE >  0) THEN
            CMESSAGE=": ERROR IN COPYDIAG(ITEM 284)"//CMESSAGE
          END IF
!
        END IF
!
!
        ITEM = 285
        IF (ICODE <= 0 .AND. SF(ITEM,SECT)) THEN
!
!  Convert scavenged dust per timestep to flux per sec
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              CONSCAV_DUST5(I,J) = CONSCAV_DUST5(I,J)/(scf*TIMESTEP)
            END DO
          END DO
!
! DEPENDS ON: copydiag
          CALL COPYDIAG(STASHWORK(SI(ITEM,SECT,IM_INDEX)),CONSCAV_DUST5,&
     &        ROW_LENGTH,ROWS,0,0,0,0,AT_EXTREMITY,                     &
     &        ATMOS_IM,SECT,ITEM,                                       &
     &        ICODE,CMESSAGE)
!
          IF (ICODE >  0) THEN
            CMESSAGE=": ERROR IN COPYDIAG(ITEM 285)"//CMESSAGE
          END IF
!
        END IF
!
!
        ITEM = 286
        IF (ICODE <= 0 .AND. SF(ITEM,SECT)) THEN
!
!  Convert scavenged dust per timestep to flux per sec
          DO J=1,ROWS
            DO I=1,ROW_LENGTH
              CONSCAV_DUST6(I,J) = CONSCAV_DUST6(I,J)/(scf*TIMESTEP)
            END DO
          END DO
!
! DEPENDS ON: copydiag
          CALL COPYDIAG(STASHWORK(SI(ITEM,SECT,IM_INDEX)),CONSCAV_DUST6,&
     &        ROW_LENGTH,ROWS,0,0,0,0,AT_EXTREMITY,                     &
     &        ATMOS_IM,SECT,ITEM,                                       &
     &        ICODE,CMESSAGE)
!
          IF (ICODE >  0) THEN
            CMESSAGE=": ERROR IN COPYDIAG(ITEM 286)"//CMESSAGE
          END IF
!
      END IF
!
      END IF !L_DUST
!
! ----------------------------------------------------------------------
! Section Sulphur Cycle scavenging diagnostics
! ----------------------------------------------------------------------
!
      item = 237                      !wet scav flux N in NH3
      IF ( sf(item,sect) ) THEN
!
!  Convert scavenged nh3 per timestep to flux per sec
        Do j=1,rows
          Do i=1,row_length
            conwash_nh3(i,j) = conwash_nh3(i,j)/timestep
          End Do
        End Do
!
! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(item,sect,im_index)),conwash_nh3,    &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage=": error in copydiag(item 237)"//cmessage
        End if
!
      END IF
!
      item = 238                      !wet scav flux SO2
      IF (icode <= 0 .and. sf(item,sect)) THEN
!
!  Convert scavenged so2 per timestep to flux per sec
        Do j=1,rows
          Do i=1,row_length
            conwash_so2(i,j) = conwash_so2(i,j)/timestep
          End Do
        End Do
!
! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(item,sect,im_index)),conwash_so2,    &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage=": error in copydiag(item 238)"//cmessage
        End if
!
      END IF
!
      item = 239                      !wet scav flux SO4_AIT
      IF (icode <= 0 .and. sf(item,sect)) THEN
!
!  Convert scavenged so4_ait per timestep to flux per sec
        Do j=1,rows
          Do i=1,row_length
            conscav_so4ait(i,j) = conscav_so4ait(i,j)/timestep
          End Do
        End Do
!
! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(item,sect,im_index)),conscav_so4ait, &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage=": error in copydiag(item 239)"//cmessage
        End if
!
      END IF
!
      item = 240                      !wet scav flux SO4_ACC
      IF (icode <= 0 .and. sf(item,sect)) THEN
!
!  Convert scavenged so4_acc per timestep to flux per sec
        Do j=1,rows
          Do i=1,row_length
            conscav_so4acc(i,j) = conscav_so4acc(i,j)/timestep
          End Do
        End Do
!
! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(item,sect,im_index)),conscav_so4acc, &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage=": error in copydiag(item 240)"//cmessage
        End if
!
      END IF
!
      item = 241                      !wet scav flux SO4_DIS
      IF (icode <= 0 .and. sf(item,sect)) THEN
!
!  Convert scavenged so4_dis per timestep to flux per sec
        Do j=1,rows
          Do i=1,row_length
            conscav_so4dis(i,j) = conscav_so4dis(i,j)/timestep
          End Do
        End Do
!
! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(item,sect,im_index)),conscav_so4dis, &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage=": error in copydiag(item 241)"//cmessage
        End if
 !
      END IF
!
! ----------------------------------------------------------------------
! Section Soot scheme scavenging diagnostics
! ----------------------------------------------------------------------
!
      item = 242                      !wet scav flux soot
      IF (icode <= 0 .and. sf(item,sect)) THEN
!
!  Convert scavenged soot per timestep to flux per sec
        Do j=1,rows
          Do i=1,row_length
            conscav_agedsoot(i,j) = conscav_agedsoot(i,j)/timestep
          End Do
        End Do
!
! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(item,sect,im_index)),                &
     &        conscav_agedsoot,                                         &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage=": error in copydiag(item 242)"//cmessage
        End if
!
      END IF
!
! ----------------------------------------------------------------------
! Section Biomass aerosol scheme scavenging diagnostics
! ----------------------------------------------------------------------
!
      item = 243                      !wet scav flux biomass
      IF (icode <= 0 .and. sf(item,sect)) THEN
!
!  Convert scavenged biomass per timestep to flux per sec
        Do j=1,rows
          Do i=1,row_length
            conscav_agedbmass(i,j) = conscav_agedbmass(i,j)/timestep
          End Do
        End Do
!
! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(item,sect,im_index)),                &
     &        conscav_agedbmass,                                        &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode >  0) Then
          cmessage=": error in copydiag(item 243)"//cmessage
        End if
!
      END IF
!
! ----------------------------------------------------------------------
! Section Fossil-fuel OC aerosol scheme scavenging diagnostics
! ----------------------------------------------------------------------
!
      item = 244                     !wet scav flux ocff
      IF (icode.le.0 .and. sf(item,sect)) THEN
!
!  Convert scavenged OCFF per timestep to flux per sec
        Do j=1,rows
          Do i=1,row_length
            conscav_agedocff(i,j) = conscav_agedocff(i,j)/timestep
          End Do
        End Do
!
! DEPENDS ON: copydiag
        Call copydiag(STASHwork(si(item,sect,im_index)),                &
     &        conscav_agedocff,                                         &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
!
        If (icode.gt.0) Then
          cmessage=": error in copydiag(item 244)"//cmessage
        End if
!
      END IF
!

 9999 continue

      If (icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      End if

      Return
      END SUBROUTINE diagnostics_conv