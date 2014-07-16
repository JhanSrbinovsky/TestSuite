#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Perform a 1-timestep integration of the Atmosphere Model
!
! Subroutine Interface:
      SUBROUTINE Atm_Step(                                              &
#include "argd1.h"
#include "argduma.h"
#include "arglndm.h"
#include "arg_atm_fields.h"
#include "argbnd.h"
#include "argsts.h"
#include "argppx.h"
! Constant arrays needed for Atmosphere-routing coupling
#include "argatcpl.h"
     & G_P_FIELD,                                                       &
     & G_R_FIELD,                                                       &
     & OBS_FLAG,OBS,obs_flag_len,obs_len,                               &
     &      ngrgas,grgas_addr,                                          &
     &                     IAU_lookup,                                  &
                                       ! inout
     &                     D1_IAU_k4,                                   &
                                       ! inout
     &                     D1_TFilt,   &! inout
     &                     endstep)
      USE rad_switches_mod
      Use cv_cntl_mod, Only:                                            &
          lcv_phase_lim, lcv_3d_cca, lcv_3d_ccw, lcv_ccrad,             &
          lcv_pc2_diag_sh             
      Use cv_run_mod,  Only:                                            &
          l_rp, l_rp2, convection_option
      Use level_heights_mod
      Use trignometric_mod
      Use dyn_coriolis_mod
      Use dyn_var_res_mod
      Use diff_coeff_mod
      Use rad_mask_trop_mod
      Use rot_coeff_mod
      use cable_iday_mod
      IMPLICIT NONE
!
! Description: Perform a 1-timestep integration of the Atmosphere Model,
!   including assimilation, physics and dynamics processing.
!
! Method: Sequential execution of code sections below, with control and
!  logical switches determining whether each section is accessed.
!
!!ND ND code section and UM added, tidy and re-number required.
! Section 0.  Initialisation.
! Section 0.1  Filter near polar winds and theta
! Section 1.0  Call Atmospheric Physics1
!     STASH diagnostics for Physics1.
! Section 2. Calculate new time-level estimates to Theta, q, qcl, qcf
! Section 2.1 Calculate advective Momentum increments.
! Section 2.2 Calculate Semi-Lagrangian part of continuity equation.
! Section 2.3 Diagnostics at end of advection.
!     STASH diagnostics for Dynamics advection.
! Section 3.0  Call Atmospheric Physics2
!     STASH diagnostics for Physics2.
! Section 4.  Form and solve Helmholtz equation and update variables.
! Section 4.1 Form and solve Helmholtz equation, return all corrections
!     Set lateral boundaries.
! Section 5.0 Calculate rho at new time level, using flux form.
! Section 6.0 Calculate u, v, w, theta, q, qcl, qcf at new time level.
! Section 7.0 Mean all polar variables on a level to remove deviation.
! Section 8.0 Update pressure to value at new-time level.
!     VAR or AC assimilation.
! Section 9.0 Diagnostics at end of timestep
!     STASH diagnostics - dynamics based (section 15).
!     STASH diagnostics - physics  based (section 16).
!     STASH diagnostics - climate end of timestep (section 30).
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
#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"

#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typcona.h"
#include "typlndm.h"
#include "ctracera.h"
#include "typ_atm_fields.h"
#include "typbnd.h"
#include "typsts.h"
#include "nstypes.h"
! River routing
#include "typatcpl.h"
      INTEGER                                                           &
     &  G_P_FIELD                                                       &
                              ! IN : global horiz domain for atmos
     &, G_R_FIELD             ! IN : global horiz domain for rivers
      INTEGER :: obs_flag_len,obs_len
      INTEGER :: OBS_FLAG(obs_flag_len)
      REAL    :: OBS(obs_len)
!
#include "chsunits.h"
#include "ccontrol.h"
#include "cruntimc.h"
#include "problem.h"
#include "surface.h"
#include "vgrid.h"
#include "o3intp.h"
#include "ctime.h"
#include "c_global.h"
#include "cphyscon.h"
#include "c_writd.h"
#include "ppxlook.h"
#include "cprintst.h"
#include "c_dust_ndiv.h"
#include "ctfilt.h"
#include "c_kinds.h"

! Variables for FLUME
      integer :: tr_size ! size of a single tracer

! Subroutine arguments:
      INTEGER :: IAU_lookup(IAU_Len1Lookup, IAU_Len2Lookup)

      REAL(KIND=real32) :: D1_IAU_k4(D1_IAU_k4_len) ! Packed IAU fields
      REAL         :: D1_TFilt (D1_TFilt_len ) ! Unpacked IAU/TDF fields
      
! 3-D fields of species to be passed down to radiation       
      INTEGER, INTENT(IN) :: ngrgas 
      INTEGER, INTENT(IN) :: grgas_addr(ngrgas)

! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='Atm_Step')

      LOGICAL gather   ! Convert to sea_ice points within BL code
      PARAMETER ( gather = .True. ) ! (was L_COMPRESS_SEAICE)

! Local scalars:
      Integer                                                           &
     &  i, j, k, l                                                      &
                          ! loop counters
     &, CycleNo                                                         &
                   ! Number of cycles (iterations) for iterative SISL.
     &, ij                                                              &
                          ! 2D array index for offx/y variables
     &, lbc_pt                                                          &
                          ! pointer to LBC array
     &, ij_u                                                            &
                          ! 2D array index for U offx/y variables
     &, ij_v                                                            &
                          ! 2D array index for V offx/y variables
     &, ji                                                              &
                          ! 2D array index for halo_i/j variables
     &, gi                                                              &
                          ! derived from loop counters
     &, lambda_start                                                    &
                    ! pointer for start of lambda_p/lambda_u on this pe

     &, first_constant_r_rho_level_m1                                   &
                                      ! value used to dimension
!                                     ! arrays, max of (1 and
!                                     ! first_constant_r_rho_level)
     &, icount                                                          &
                ! diagnostic count
     &, i_start, i_stop                                                 &
     &, j_start, j_stop                                                 &
     &, j_begin, j_end                                                  &
     &, j_ti, j_ice_fract, j_ice_thick

      Integer                                                           &
     &  n_Y_arrays                                                      &
                         ! = 1 for global, 3 for LAM
     &, n_Yw_arrays                                                     &
                         ! = 1 for global, 2 for LAM
     &, n_Yd_arrays                                                     &
                         ! = 1 for global, 3 for LAM
     &, n_Ydw_arrays     ! = 1 for global, 2 for LAM
!
      Integer :: nd_o3  ! Total size of ozone array supplied
      Integer :: nd_stochem  ! size of stochem ozone array

      Integer                                                           &
     &  check_bottom_levels ! used in interpolation code, and is
!                           ! the number of levels to check to see
!                           ! if the departure point lies inside the
!                           ! orography.

      INTEGER                                                           &
     &  timestep_number                                                 &
                             ! no. of atmos timesteps since basis time
     &,first_constant_r_rho_level                                       &
                                  ! 1st rho level on which r is constant
     &, info       ! icode return from UM_FORT_FLUSH

      INTEGER                                                           &
     &  DIM_CS1                                                         &
                         ! soil C dimension: 1 for single, 4 for RothC
     &, DIM_CS2          ! soil C dimension: 1 for single,
                         !                        LAND_FIELD for RothC
      REAL                                                              &
     & constant                                                         &
     &,h_print                                                          &
                         ! needed for printing idealised orography
     &,timestep                                                         &
                            ! atmosphere model timestep
     &,radiation_timestep                                               &
                            ! timestep of radiation scheme
     &,radiation_tstep_diag                                             &
                              ! timestep of fast radiation (3C)
     &,radiation_tstep_prog                                             &
                              ! timestep of slow radiation (3C)
     &,chemistry_timestep             ! must be  <=  model timestep

      REAL                                                              &
     & pos_timestep                                                     &
                    ! = +timestep.
     &,neg_timestep ! = -timestep.

      LOGICAL                                                           &
     & l_rad_step                                                       &
                            ! :T :activate radiation this timestep
     &,L_Rad_step_diag                                                  &
                         ! T :activate fast radiation this timestep (3C)
     &,L_Rad_step_prog                                                  &
                         ! T :activate slow radiation this timestep (3C)
     &,L_print_L2norms                                                  &
                         ! diagnostic printing of l2norms
     &,L_print_L2helm                                                   &
                         ! diagnostic printing of l2norms in solver
     &,L_Tracer          ! T if *any* tracer variables present

! Code to do with tropopause diagnostics from O3_to_3D
! Declare logicals for Stash Diagnostics used in call to O3_to_3D
      Logical :: L_O3_trop_level   ! STASH code 2,280
      Logical :: L_O3_trop_height  ! STASH code 2,281
      Logical :: L_T_trop_level    ! STASH code 2,282
      Logical :: L_T_trop_height   ! STASH code 2,283

! Declare the tropopause variables output for O3_to_3D as allocatable
      Real, DIMENSION (:,:), ALLOCATABLE::                              &
     & O3_trop_level, T_trop_level, O3_trop_height, T_trop_height

! END OF Code to do with tropopause diagnostics from O3_to_3D

      INTEGER                                                           &
     & LAND_PTS_TRIF                                                    &
                     !\ For dimensioning variables in NI_bl_ctl
     &,NPFT_TRIF                                                        &
                     !/ depending on whether TRIFFID is in use
     &,CO2_DIM_LEN                                                      &
                     !\ For dimension 3-D CO2 field to be passed
     &,CO2_DIM_ROW                                                      &
                     !/ to NI_bl_ctl
     &,CO2_DIM_LEV   !/ and NI_rad_ctl

      INTEGER :: Sec         ! Second
      LOGICAL :: L_CallTFilt ! Call TFilt_cntl?

! Variables required for call to SET_LATERAL_BOUNDARIES
      INTEGER                                                           &
     &  lbc_size           ! size of a single level of LBC

      LOGICAL                                                           &
     &  L_do_halos                                                      &
                          ! update the halos?
     &, L_do_boundaries
                          ! update the boundaries?


! Variables required for ice category selection
      REAL, POINTER :: p_ti(:,:,:)
      REAL, POINTER :: p_ice_fract(:,:,:)
      REAL, POINTER :: p_ice_thick(:,:,:)

! soil carbon content & accumulated soil respiration 
! (target varies according to L_TRIFFID)
      REAL, POINTER :: CS(:), RSA(:)

! Array dimensions for sea-salt aerosol
      INTEGER                                                           &
     &  salt_dim1                                                       &
     &, salt_dim2                                                       &
     &, salt_dim3

! Array dimensions for Aero_Ctl
      INTEGER                                                           &
     & aero_dim1                                                        &
     &,aero_dim2                                                        &
     &,aero_dim3

! CABLE declarations
#include "c_soilh.h"
      INTEGER                               &
     &  TILE_INDEX(land_points,NTYPE)       &! OUT Index of tile points     
     &, TILE_PTS(NTYPE)                     &! OUT Number of tile points
     &, day                                 &
     &, total_nsteps                        &! Total number of steps in run
!     &, istep_cur                           &
     &, SOIL_TYPE(row_length,rows)          &
     &, VEG_TYPE(row_length,rows)
      INTEGER, save :: istep_cur=0          ! current time step  
          
      LOGICAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: L_TILE_PTS
!      LOGICAL l_cable
!*** note in ni-rad_ctl & glue_rad alb_tile(land_field,ntiles) dim ***

! CABLE prognostic variables
! Declared as allocable here so they can be saved. Eventually they'll all
! move to the D1 array
      Real, DIMENSION(:,:,:), ALLOCATABLE, SAVE ::                      &
     &  SNOW_COND                                                       &
     &, STHU_TILE
      Real, DIMENSION(:,:), ALLOCATABLE, SAVE ::                        &
     &  T_SURF_TILE

!code differences:
      ! Some of these should eventually be promoted to proper prognostic
      ! variables. Just save them for now - MRD
      ! Save doesn't work because they're automatic variables.
!!$      save SNAGE_TILE, RTSOIL_TILE, GFLUX_TILE, SGFLUX_TILE, SNOW_DEPTH3L    &
!!$     , SNOW_MASS3L, SNOW_COND, SNOW_TMP3L, SNOW_RHO3L, SNOW_RHO1L, SMCL_TILE &
!!$     , STHU_TILE, STHF_TILE, TSOIL_TILE, T_SURF_TILE, HCONS, alb_tile        &
!!$     , TOT_ALB, CH_CAB, CD_CAB, U_S_CAB, LAND_ALBEDO_CABLE
!code differences:
!      Call RUN_INIT(                                ...                 &
!     &  ,l_use_seasalt_direct, l_use_biogenic, l_murk_rad, l_use_aod    &
!atmos_physics2 needs:

      Real, Dimension(:,:,:), ALLOCATABLE, SAVE ::                            & 
      ALB_TILE,SURF_DOWN_SW,LAND_ALBEDO_CABLE
 
      Real, Dimension(:,:), ALLOCATABLE, SAVE ::                              & 
      FTL_TILE_CAB,LE_TILE_CAB,TSTAR_TILE_CAB,SMCL_CAB,TSOIL_CAB,RADNET_TILE, & 
      SW_DOWN,RTSOIL_TILE,GFLUX_TILE,SGFLUX_TILE,TOT_ALB,CD,CHX,TRANSP_TILE,  &
      WB_LAKE
 
      Real, Dimension(:), ALLOCATABLE, SAVE ::                                & 
      FTL_CAB,LE_CAB,TSTAR_CAB,USTAR_CAB,SURF_HTF_CAB,CH_CAB,CD_CAB,U_S_CAB,  & 
      HCONS,ALBSOIL,TOT_WBLAKE,TOT_SUBRUN

      Real, allocatable, save ::                                        &
     & wblake_ratio    

      Real                                                              &
                lat(row_length, rows),long(row_length, rows),time_sec
                                          ! Lat. of gridpoint chosen
                                          ! Long. of gridpoint chosen
                                          ! actual time of day in secs.
! CABLE end ----------------------------------------------------

! Local Arrays
      Real                                                              &
     &  theta_star(1-offx:row_length+offx,                              &
     &               1-offy:rows+offy, model_levels)                    &
     &, q_star(1-offx:row_length+offx,                                  &
     &           1-offy:rows+offy, wet_levels)                          &
     &, qcl_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, qcf_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, cf_star(1-offx:row_length+offx,                                 &
     &           1-offy:rows+offy, wet_levels)                          &
     &, cfl_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, cff_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, exner_star(1-offx:row_length+offx,                              &
     &             1-offy:rows+offy, wet_levels)                        &
     &, temp_ozone(1-offx:row_length+offx,                              &
     &             1-offy:rows+offy, model_levels)                      &
     &, frac_control(land_field,ntype)   !Forcing for land surface (3C)

      Real, DIMENSION (:,:,:), ALLOCATABLE ::                           &
        ! Local Arrays to store microphysics fields
     &  qcf2_star, qrain_star, qgraup_star

! Declare allocatable arrays for passing cloud fractions
! to LS_ACF_Brooks
      Real, DIMENSION (:,:,:), ALLOCATABLE::                            &
     & cf_bulk_nohalo, cf_liquid_nohalo, cf_frozen_nohalo

      Real                                                              &
     &  R_u(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &        model_levels)                                             &
     &, R_v(1-offx:row_length+offx, 1-offy:n_rows+offy,                 &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels)

      Real                                                              &
     &  rho_n (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)

      Real                                                              &
     &   biogenic(row_length, rows, model_levels)

! Local arrays for using the aerosol climatology for NWP

#include "arcl_dim.h"
      
      ! Internal model switches
      Logical L_USE_ARCL(NPD_ARCL_SPECIES)
      
      ! Internal array of mass-mixing ratios
      Real, dimension(:,:,:,:), allocatable :: arcl
      
      ! Number of requested species within the climatology
      Integer n_arcl_species
      
      ! Corresponding number of requested components
      Integer n_arcl_compnts
      
      ! Array index of each component
      Integer i_arcl_compnts(NPD_ARCL_COMPNTS)
      
      Integer cloud_tol

      Real                                                              &
     &  mag_vector_np (model_levels)                                    &
     &, dir_vector_np (model_levels)                                    &
     &, mag_vector_sp (model_levels)                                    &
     &, dir_vector_sp (model_levels)                                    &
     &, lambda_a (row_length) ! delta_lambda term for polar wind

      REAL, DIMENSION (:,:,:), ALLOCATABLE::                            &
     &  inc_u, inc_v, inc_w, inc_t, inc_rho                             &
     &, inc_q, inc_qcl, inc_qcf, inc_cf, inc_cfl, inc_cff

      Logical                                                           &
     &  L_do_inc_vels                                                   &
     &, L_do_inc_rho               ! flag for rho incr diagnostic

! LAM LBC tendency
      REAL                                                              &
     &  U_LBC_REAL_TEND(LENRIMA(fld_type_u,halo_type_extended,          &
     &                    rima_type_norm),MODEL_LEVELS)                 &
     &, V_LBC_REAL_TEND(LENRIMA(fld_type_v,halo_type_extended,          &
     &                    rima_type_norm),MODEL_LEVELS)                 &
     &, W_LBC_REAL_TEND(LENRIMA(fld_type_p,halo_type_extended,          &
     &                    rima_type_norm),MODEL_LEVELS)                 &
     &, EXNER_LBC_REAL_TEND(LENRIMA(fld_type_p,halo_type_extended,      &
     &                       rima_type_norm),MODEL_LEVELS+1)

! Physics arrays needed by dynamics
      Real                                                              &
     &  rho_km(1-offx:row_length+offx, 1-offy:rows+offy,                &
     &           0:bl_levels-1)                                         &
     &, cH(1-offx:row_length+offx, 1-offy:rows+offy,                    &
     &       model_levels-1)                                            &
     &, wet_to_dry_n (1-offx:row_length+offx, 1-offy:rows+offy,         &
     &          model_levels)                                           &
     &, wet_to_dry_np1 (1-offx:row_length+offx, 1-offy:rows+offy,       &
     &          model_levels)

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
!                                 BL (LW,SW) rad heating rates
     &, surf_radflux(row_length, rows)                                  &
     &, dolr(row_length,rows)                                           &
!       local field "dolr" is distinguished from "dolr_field" (in atm_fields_mod)
                                   ! TOA - surface upward LW
     &, SW_tile(land_field,ntiles)                                      &
                                   ! Surface net SW on land tiles
     &, cos_zenith_angle(row_length, rows)


! MPP-related Arrays
      INTEGER                                                           &
     & g_row_length(0:nproc-1)                                          &
                               ! Table of number of points on a row
     &,g_rows(0:nproc-1)                                                &
                               ! Table number of rows in theta field
     &,g_i_pe(1-halo_i:global_row_length+halo_i) ! processor on my
!               processor-row holding a given value in i direction

! Workspace defined as allocatable arrays, since they each communicate
! fields between near adjacent calls and only need to use memory for
! a subset of the total routine.
      REAL, DIMENSION (:,:,:), ALLOCATABLE::                            &
     &  exner_prime                                                     &
                        ! solution to helmholtz solver
     &, ozone3D                                                         &
                        ! 3d ozone (expanded from zonal) for radiation
     &, dtheta_dr_term                                                  &
     &, depart_lambda, depart_phi, depart_r_theta, depart_r_w

! Allocatable arrays for use in AC_CTL call
      REAL, DIMENSION (:,:,:), ALLOCATABLE::                            &
     &  work_q, work_qcl, work_qcf

      INTEGER                                                           &
     & ErrorStatus      ! Return code : 0 Normal Exit : >0 Error

      CHARACTER*256                                                     &
     & CMessage         ! Error message if return code >0

! Local dynamic arrays:
      REAL, DIMENSION (:), ALLOCATABLE::                                &
!  STASHworki = STASHwork for section i
     & STASHwork1,STASHwork2,STASHwork3,STASHwork4,STASHwork5           &
     &,STASHwork6,STASHwork8,STASHwork9,STASHwork12,STASHwork13         &
     &,STASHwork14,STASHwork17,STASHwork19,STASHwork26,STASHwork30      &
     &,STASHwork10,STASHwork18

! increment diagnostics:
      REAL, DIMENSION (:,:,:), ALLOCATABLE::                            &
     & u_incr_diagnostic, v_incr_diagnostic, T_incr_diagnostic          &
     &,q_incr_diagnostic ,qcl_incr_diagnostic, qcf_incr_diagnostic      &
     &,cf_incr_diagnostic, cfl_incr_diagnostic, cff_incr_diagnostic     &
     &,w_incr_diagnostic

      REAL, DIMENSION (:,:), ALLOCATABLE::                              &
     & w_local_mask
      REAL STASHwork0_dummy(1) ! STASHwork not defined for section 0,
                               !  but required as dummy argument.

! Local arrays for phys1 and phys2 increments for tracers:
      Real :: super_tracer_phys1(1-halo_i:row_length+halo_i,            &
                                 1-halo_j:rows+halo_j,                  &
                                 model_levels, super_array_size)
      Real :: super_tracer_phys2(row_length,                            &
                                 rows,                                  &
                                 model_levels, super_array_size)
      Real :: tracer_phys1(1-halo_i:row_length+halo_i,                  &
                           1-halo_j:rows+halo_j,                        &
                           tr_levels, tr_vars)
      Real :: tracer_phys2(row_length,                                  &
                           rows,                                        &
                           tr_levels, tr_vars)

      Real :: unscaled_dry_rho(1-offx:row_length+offx,                  &
                               1-offy:rows+offy, model_levels)

! Local dynamic arrays for phys1 and phys2 increments for moisture:
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     &  q_phys1, qcl_phys1, qcf_phys1, q_phys2, qcl_phys2, qcf_phys2    &
     &, qcf2_phys1, qrain_phys1, qgraup_phys1                           &
     &, qcf2_phys2, qrain_phys2, qgraup_phys2                           &
     &, mix_v_phys1, mix_cl_phys1, mix_cf_phys1                         &
     &, mix_v_phys2, mix_cl_phys2, mix_cf_phys2                         &
     &, mix_cf2_phys1, mix_rain_phys1, mix_graup_phys1                  &
     &, mix_cf2_phys2, mix_rain_phys2, mix_graup_phys2                  &
     &, cf_phys1, cfl_phys1, cff_phys1, cf_phys2, cfl_phys2, cff_phys2

! local dynamic arrays for PC2
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     &  t_inc_pres, q_inc_pres, qcl_inc_pres, qcf_inc_pres, cf_inc_pres &
     &, cfl_inc_pres, cff_inc_pres                                      &
     &, t_dini, q_dini, qcl_dini, qcf_dini, cf_dini                     &
     &, cfl_dini, cff_dini, rhts, qtts, tlts, ptts

!
! Extra variables needed for cycling
! Vars ending in _phys1 are copies holding the value the original
! variable had after exiting phys1.
! Vars ending in _np1 are tn+1 estimates holding the value the original
! variable had at the end of the last cycle (provided that CycleNo>1).
! obtained from the last
! cycle when CycleNo>1.
!
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     &  R_u_phys1, R_v_phys1, thetastar_phys1, qstar_phys1              &
     &, qclstar_phys1, qcfstar_phys1, qcf2_star_phys1, qrain_star_phys1 &
     &, qgraup_star_phys1, ti_phys1, cca_phys1, area_cld_frac_phys1     &
     &, bulk_cld_frac_phys1, bulk_cld_liq_phys1, bulk_cld_fr_phys1      &
     &, u_np1, v_np1, w_np1, theta_np1, rho_np1, q_np1, qcl_np1         &
     &, qcf_np1, qcf2_np1, qrain_np1, qgraup_np1

      REAL, DIMENSION (:,:), ALLOCATABLE::                              &
     &  z0msea_phys1, zh_phys1, t_land_ctile_phys1, t_sice_ctile_phys1, &
     &  t_surf_phys1, t_sf_tile_phys1, snow_tile_phys1, dolr_phys1,     &
     &  RHO_LBC_REAL_TEND

      INTEGER, DIMENSION (:,:), ALLOCATABLE :: ccb_phys1, cct_phys1

      REAL GS1(land_field)

! local variable
      REAL                                                              &
     &  tot_dry_mass_final                                              &
                            ! mass at end of energy correction period
     &, tot_energy_final                                                &
                            ! energy at end of energy correction period
     &, tot_moist_final                                                 &
                            ! moist at end of energy correction period
     &, energy_corr_now                                                 &
                            ! instanteous energy correction
     &, increment_factor    ! For calculating value of LBC at next TL

! Monsoon variables
      Real                                                              &
     &  lambda_half_width                                               &
     &, phi_half_width                                                  &
     &, p_max                                                           &
     &, p_top                                                           &
     &, p_bottom                                                        &
     &, p_arbitrary                                                     &
     &, lambda_heat_centre                                              &
     &, phi_heat_centre                                                 &
     &, max_heat_per_day                                                &
     &, Mons_newtonian_timescale

!    Local arrays for when using mixing ratios
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     &  mix_v, mix_cl, mix_cf                                           &
     &, mix_v_star, mix_cl_star, mix_cf_star                            &
     &, mix_v_inter, mix_cl_inter, mix_cf_inter                         &
     &, mix_cf2, mix_rain, mix_graup                                    &
     &, mix_cf2_star, mix_rain_star, mix_graup_star                     &
     &, mix_cf2_inter, mix_rain_inter, mix_graup_inter                  &
     &, q_store, qcl_store, qcf_store                                   &
     &, qcf2_store, qrain_store, qgraup_store
!
! Additional variables needed for cycling when mixing ratios are used
!
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     &  mix_v_star_phys1, mix_cl_star_phys1, mix_cf_star_phys1          &
     &, mix_cf2_star_phys1, mix_rain_star_phys1, mix_graup_star_phys1   &
     &, mix_v_np1, mix_cl_np1, mix_cf_np1                               &
     &, mix_cf2_np1, mix_rain_np1, mix_graup_np1

!  Suarez-Held variables now declared in CRUNTIMC
      Real, ALLOCATABLE :: RHCPT(:,:,:)
      Integer rhc_row_length
      Integer rhc_rows
      INTEGER    :: minutos ! LOCAL Store value of timestep in min.
!
!  Dummy variables for SCM Diagnostics,
!  for passing down to atmos_physics1 and 2 routines
      INTEGER, PARAMETER :: nSCMDpkgs = 12
      LOGICAL L_SCMDiags(nSCMDpkgs)
      LOGICAL                                                           &
     & L_flux_bc                 ! T if prescribed surface fluxes to be used

      REAL                                                              &
     &  flux_e(row_length, rows)                                        &
                                 ! Surface latent heat flux (W/m^2)
     &, flux_h(row_length, rows)                                        &
                                 ! Surface sensible heat flux (W/m^2)
     &, z0m_scm(row_length, rows)                                       &
                                 ! SCM specified z0m (m)
     &, z0h_scm(row_length, rows)! SCM specified z0m (m)
      REAL, DIMENSION(:,:,:), ALLOCATABLE ::                            &
     &  q_inc_subs, th_inc_subs                                         &
                                  ! subsidence increments
     &, q_inc_ls, th_inc_ls                                             &
                                  ! large scale increments
     &, u_inc_dmp, q_inc_dmp, th_inc_dmp                                &   
                                  !Damping incs
     &, v_inc_dmp
! Tolerance for CycleNo >1
      REAL                                                              &
     &  GCR_run_tol_abs                                                 &
     &, GCR_run_tol_res

! Subgrid turbulence scheme variables.
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     & visc_m                                                           &
!            ! diffusion coefficient for momentum from subgrid
!            ! turbulence scheme
     &,visc_h                                                           &
!            ! diffusion coefficient for heat and moisture from
!            ! subgrid turbulence scheme
     &,FM_3D                                                            &
!            ! stability function for momentum transport.
!            ! level 1 value is dummy for use in diagnostics
     &,FH_3D                                                            &
!            ! stability function for heat and moisture.
!            ! level 1 value is dummy for use in diagnostics
     &,visc_BL_m                                                        &
                     ! visc_m only on BL levels
     &,RNEUTML                                                          &
                     ! mixing length scale (m) (lambda)
     &,shear                                                            &
                     ! S from subgrid turbulence scheme
     &,BL_COEF_KM                                                       &
                     ! RHOKM from BL scheme
     &,BL_COEF_KH    ! RHOKH from BL scheme

       Real                                                             &
     & max_diff  ! max diffusion coeff for run
! Function & Subroutine calls:
      External                                                          &
     & Swap_Bounds                                                      &
     &,Polar_filter                                                     &
     &,Atmos_Physics1                                                   &
     &,STASH                                                            &
     &,SL_Thermo                                                        &
     &,SL_Full_wind                                                     &
     &,Atmos_Physics2                                                   &
     &,dry_static_adj                                                   &
     &,Diff_Divdamp_Ctl                                                 &
     &,pg_update                                                        &
     &,Polar_filter_incs                                                &
     &,Polar_vector_wind_n                                              &
     &,PE_Helmholtz_eul                                                 &
     &,Flux_Rho                                                         &
     &,Polar_Reset_Mean                                                 &
     &,Calc_P_star                                                      &
     &,Calc_Exner_at_theta                                              &
     &,Calc_P_from_Exner                                                &
     &,BOUNDVAL                                                         &
     &,AC_CTL                                                           &
     &,St_diag1                                                         &
     &,St_diag2                                                         &
     &,St_diag3                                                         &
     &,Ereport                                                          &
     &,bottom_w_Calc                                                    &
     &,TFilt_cntl                                                       &
     &,LS_ACF_Brooks                                                    &
     &,set_arcl_dimensions                                              &
     &,set_arcl_clim                                                    &
     &,exppxi               ! Function to extract ppxref info


      Integer                                                           &
     &     exppxi               ! Function to extract ppxref info

      LOGICAL                                                           &
     &  first_atmstep_call                                              &
     &, L_update_lbcs                                                   &
     &, L_apply_lbcs                                                    &
     &, GCR_zero_guess_it

      Integer :: itemp

! time-stepping weights. Values may be different
! at different cycles (UMUI controlled).
      REAL                                                              &
     & alpha1, alpha2, alpha3, alpha4

      ! Local parameters for mixing ratio physics
      Logical, Parameter ::                                             &
                                ! Mixing ratios for atmos_physics_1 and 
                                ! _2 are defined through namelist
     & l_mr_acctl    = .false.                                          &
                                ! Use mixing ratios for ac_ctl
     &,l_mr_qtbalcld = .false.                                          &
                                ! Use mixing ratios for qt_bal_cld
     &,l_mr_tfiltctl = .false.                                          &
                                ! Use mixing ratios for tfilt_ctl
     &,l_mr_pc2      = .false.  ! Use mixing ratios for PC2 routines

      SAVE first_atmstep_call

      DATA first_atmstep_call /.TRUE./

      Logical :: L_physics_store
!
! Oxidant mass-mixing ratios and concentrations, for use in sulphur
! cycle.
      REAL, DIMENSION(:,:,:), ALLOCATABLE ::                            &
     &  O3_MMR                                                          &
     &, H2O2_MMR                                                        &
     &, OH_conc                                                         &
     &, HO2_conc 

! Array to hold atmospheric loss for tracers - kdcorbin, 05/10
   REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: lossrate

   !end step of experiment
   integer :: endstep
!
!- End of header

! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

!! FLUME 
      tr_size = (tr_levels * theta_off_size) ! size of a single tracer

!glr cable 2 lines ----------------------------------------------------
!      L_CABLE = .false. !.true. ! to call CABLE paralell to MOSES 
!      l_cable = .true.  ! to call CABLE instead of MOSES

!      print *,'atm_step before allocation',istep_cur,land_points,sm_levels,NTILES
!      print *,'atm_step: t_i, mype, landpoints',istep_cur, mype,land_points

!     Temporary CABLE initialisation - MRD
      if ( first_atmstep_call ) then
         allocate(SNOW_COND(land_points,NTILES,3))
         allocate(STHU_TILE(land_points,NTILES,sm_levels))
!!$         allocate(STHF_TILE(land_points,NTILES,sm_levels))
         allocate(T_SURF_TILE(land_points,NTILES))
!pfu MOD
         allocate(L_TILE_PTS(land_points,NTILES))
      ALLOCATE(FTL_TILE_CAB(LAND_POINTS,NTILES)   )
      ALLOCATE(FTL_CAB(LAND_POINTS)               )
      ALLOCATE(LE_TILE_CAB(LAND_POINTS,NTILES)    )
      ALLOCATE(LE_CAB(LAND_POINTS)                )
      ALLOCATE(TSTAR_TILE_CAB(LAND_POINTS,NTILES) )
      ALLOCATE(TSTAR_CAB(LAND_POINTS)             )
      ALLOCATE(SMCL_CAB(LAND_POINTS,SM_LEVELS)    )
      ALLOCATE(TSOIL_CAB(LAND_POINTS,SM_LEVELS)   )
      ALLOCATE(USTAR_CAB(LAND_POINTS)                         )
      ALLOCATE(SURF_HTF_CAB(LAND_POINTS)          )
      ALLOCATE(RADNET_TILE(LAND_POINTS,NTILES)    )
      ALLOCATE(CH_CAB(land_points)    )
      ALLOCATE(CD_CAB(land_points)     )
      ALLOCATE(U_S_CAB(land_points)     )
      ALLOCATE(SW_DOWN(row_length,rows)           )
      ALLOCATE(RTSOIL_TILE(land_points,NTILES)    )
      ALLOCATE(GFLUX_TILE(land_points,NTILES)     )
      ALLOCATE(SGFLUX_TILE(land_points,NTILES)    )
      ALLOCATE(HCONS(land_points)                 )
      ALLOCATE(alb_tile(land_points,ntiles,4)     )
      ALLOCATE(surf_down_sw(row_length,rows,4)    )
      ALLOCATE(TOT_ALB(land_points,ntiles)        )
      ALLOCATE(albsoil(land_points)                )
      ALLOCATE(CD(row_length,rows)                )
      ALLOCATE(CHX(row_length,rows)               )
      ALLOCATE(LAND_ALBEDO_CABLE(row_length,rows,4))
      ALLOCATE(TRANSP_TILE(LAND_POINTS,NTILES))
      ALLOCATE(WB_LAKE(LAND_POINTS,NTILES))
      ALLOCATE(TOT_WBLAKE(LAND_POINTS))
      ALLOCATE(TOT_SUBRUN(LAND_POINTS))
      ALLOCATE(wblake_ratio)
         SNOW_COND = -huge(1.)
         STHU_TILE = -huge(1.)
         T_SURF_TILE = -huge(1.)
         L_TILE_PTS = .false.
         if (l_cable) then
          wblake_ratio = 0.0
          ! need init ?
          !WB_LAKE     = 0.0
          !TOT_WBLAKE  = 0.0
          !TOT_SUBRUN  = 0.0
         end if
      end if

       R_w = 0.


! Set non-initialised variables to one

       gi             = 1
       lbc_size       = 1
       ij_v           = 1
       ij_u           = 1
! set Error code to zero
      ErrorStatus = 0

! set the L_Tracer flag
      L_Tracer = ( tr_vars > 0 .OR. L_soot .OR. L_CO2_interactive .OR.  &
     &         L_sulpc_SO2 .OR. L_murk .OR. L_dust .OR. L_biomass .OR.  &
     &         L_USE_CARIOLLE .OR. L_ocff .OR. tr_ukca > 0 )
! timestep information
      timestep = SECS_PER_STEPim(atmos_im) ! timestep in seconds
      radiation_timestep = timestep * A_SW_RADSTEP
      radiation_tstep_diag = timestep * A_SW_RADSTEP_DIAG
      radiation_tstep_prog = timestep * A_SW_RADSTEP_PROG
      timestep_number = STEPim(atmos_im)   ! no. of steps since basis
      l_rad_step = L_SW_RADIATE .OR. L_LW_RADIATE ! Activate radiation
! Activate radiation
      l_rad_step_diag = L_SW_RADIATE_DIAG .OR. L_LW_RADIATE_DIAG
      l_rad_step_prog = L_SW_RADIATE_PROG .OR. L_LW_RADIATE_PROG

      pos_timestep =  timestep
      neg_timestep = -timestep

      IF (L_Backwards) timestep = neg_timestep

      IF (PrintStatus  >=  PrStatus_Normal) THEN

        IF (timestep_number  ==  1 .AND. L_Backwards) THEN
          WRITE(6,*) ''
          WRITE(6,*) '          *************************'
          WRITE(6,*) 'Atm_Step: * INTEGRATING BACKWARDS *'
          WRITE(6,*) '          *************************'
          WRITE(6,*) ''
        ENDIF

        IF (timestep_number  ==  1 .AND. L_RHCPT) THEN
          WRITE(6,*) ''
          WRITE(6,*) 'Atm_Step: Running with diagnostic RHcrit'         &
     &                      //' option active'
          WRITE(6,*) ''
        ENDIF

        IF( L_print_pe .or. mype ==0 ) then
          WRITE(6,*) 'Atm_Step: Timestep ', timestep_number
        ENDIF ! L_print_pe .or. mype ==0

      ENDIF     ! PrintStatus

      IF ((model_domain == mt_global) .and. (L_trivial_trigs)) THEN
        ErrorStatus=123
! DEPENDS ON: ereport
        Call Ereport("ATM_STEP", ErrorStatus,                           &
     &        "Unable to run global model with L_trivial_trigs.")
      ENDIF

      ! Lest 16/12/13 - moved up to allocate loop
      !if (l_cable) then
      ! if (first_atmstep_call) wblake_ratio = 0.0
      !end if

! Set radiation switches in rad_switches_mod
      lrad_ctile_fix    = l_ctile_fix
      lrad_cldtop_t_fix = l_cldtop_t_fix
      lrad_quad_src_fix = l_quad_src_fix
      lrad_ccrad        = l_ccrad
      lrad_3d_ccw       = l_3d_ccw
      lrad_ovrlap       = l_ovrlap
      lrad_ccw_scav     = l_ccw_scav
      lrad_emis_land_gen= l_emis_land_gen
      rad_emis_land_gen = emis_land_gen
! sza
! set cloud inhomegenous option
      if ( .not. L_INHOM_CLOUD ) then
         LRAD_TRIPLECLOUDS = .TRUE.
      endif
! sza end

! Set convection switches in cv_cntl_mod
      lcv_phase_lim     = l_phase_lim
      lcv_3d_cca        = l_3d_cca
      lcv_3d_ccw        = l_3d_ccw
      lcv_ccrad         = l_ccrad


      lcv_pc2_diag_sh   = l_pc2_diag_sh
! set up logical switch for diagnostic printing of l2norms
      L_print_L2norms = .false.
      L_print_L2helm = .false.
      itemp = timestep_number - first_norm_print
      If( print_step > 0 .and. itemp >=0) then
        If( mod( itemp , print_step ) == 0 ) then
          If ( L_diag_L2norms ) then
            if( itemp / print_step  > 11 ) then
         write(6,*)'l2norms printing too often, limited to 12 occasions'
            else
            L_print_L2norms = .true.
            endif  ! itemp / print_step > 11
          endif ! L_diag_L2norms
          If(  L_diag_L2helm ) then
            if( itemp / print_step  > 11 ) then
         write(6,*)'l2norms printing too often, limited to 12 occasions'
            else
            L_print_L2helm = .true.
            endif  ! itemp / print_step > 11
          endif !  L_diag_L2helm
        endif !  mod( itemp , print_step ) == 0)
      endIf ! print_step > 0

! grid information
      If (L_initialise_data .AND.                                       &
     &    first_constant_r_rho_level_new /= -1) Then
        first_constant_r_rho_level = first_constant_r_rho_level_new
      else   !   /=  L_initialise_data
        first_constant_r_rho_level = a_inthd(ih_1_c_rho_level)
      endif  !  L_initialise_data

! Set mpp arrays from parvars.h information

      DO i= 0,nproc-1
        g_row_length(i) = g_lasize(1,fld_type_p,halo_type_no_halo,i)
        g_rows      (i) = g_lasize(2,fld_type_p,halo_type_no_halo,i)
      ENDDO ! i processors

      DO i= 1-halo_i, 0
        g_i_pe(i) = 0
      ENDDO ! i

      DO i= 1,global_row_length
        g_i_pe(i) = g_pe_index_EW(i)
      ENDDO ! i row points

      DO i= global_row_length+1, global_row_length+halo_i
        g_i_pe(i) = nproc_x-1
      ENDDO ! i

! End of Set mpp arrays from parvars.h information

      If (first_constant_r_rho_level  ==  1 ) Then
        first_constant_r_rho_level_m1 = first_constant_r_rho_level
      Else
        first_constant_r_rho_level_m1 = first_constant_r_rho_level - 1
      End If

      lambda_start = datastart(1) - halo_i
      If ( L_regular ) Then
        Do i = 1, row_length
          gi = datastart(1) + i - 1
          lambda_a(i) = (gi - .5) * delta_lambda
        endDo
      else  ! variable resolution
        Do i = 1, row_length
          gi = datastart(1) + i - 1
          lambda_a(i) = glambda_u(gi) - base_lambda
        endDo
      endIf ! L_regular

! set number of levels to check for trajectory inside orography at
! bottom of model, the same parameter is also used inside the
! interpolation routine.

      check_bottom_levels = max(bl_levels,                              &
     &                          interp_vertical_search_tol )

! Define number of work arrays needed in sl_full_wind routine depending
! on model domain choice.
      If (model_domain  ==  mt_global .or.                              &
     &    model_domain  ==  mt_cyclic_lam  .or.                         &
     &    model_domain  ==  mt_bi_cyclic_lam ) Then
        n_Y_arrays    = 1
        n_Yw_arrays   = 1
        n_Yd_arrays   = 1
        n_Ydw_arrays  = 1
      Else If (model_domain  ==  mt_lam) Then
        n_Y_arrays    = 3
        n_Yw_arrays   = 2
        n_Yd_arrays   = 3
        n_Ydw_arrays  = 2
      End If

!  i_start, i_stop, j_start, j_stop define the solution domain
!  j_begin, j_end exclude the poles
      i_start = 1
      i_stop = row_length
      j_start = 1
      j_stop = rows
      j_begin = 1
      j_end = rows
!  mt_bi_cyclic_lam keeps defaults; other domains change as below
      If (model_domain == mt_global) Then
        If (at_extremity(PSouth)) j_begin = 2
        If (at_extremity(PNorth)) j_end = rows - 1
      EndIf  ! model_domain  ==  mt_global
      If (model_domain == mt_lam) Then
        If (at_extremity(PSouth)) then
          j_start = 2
        EndIf ! at_extremity(PSouth)
        If (at_extremity(PNorth)) then
          j_stop = rows - 1
        EndIf ! at_extremity(PNorth)        
        If(at_extremity(PWest)) i_start = 2
        If(at_extremity(PEast)) i_stop = row_length - 1
      End If  ! model_domain  ==  mt_LAM
      If (model_domain == mt_cyclic_lam) Then
        If (at_extremity(PSouth)) then
          j_start = 2
          j_begin = j_start
        EndIf ! at_extremity(PSouth)
        If (at_extremity(PNorth)) then
          j_stop = rows - 1
          j_end = j_stop
        EndIf ! at_extremity(PNorth)        
      End If  ! model_domain ==  mt_cyclic_lam

!  Set variables for surface forcing to false/0.0 for standard UM
      L_flux_bc     = .FALSE.
      flux_e(:,:)  = 0.0
      flux_h(:,:)  = 0.0
! Set variables for roughness length to false/0.0 if not initialising      
      If (.NOT. L_initialise_data) Then
        L_spec_z0     = .FALSE.
        z0m_scm(:,:) = 0.0
        z0h_scm(:,:) = 0.0
      End If
!  Set dummy variables for SCM Diagnostics to false for full UM
      L_SCMDiags(1:nSCMDpkgs)  = .FALSE.
!
! ---------------------------------------------------------------------
! Section 0.1  Initialisation for idealised test problems
!              For standard runs go to end section 0.1
! ---------------------------------------------------------------------

      h_print=0.0

      If (L_subfilter_vert) Then

        ALLOCATE (BL_COEF_KM(1:row_length, 1:rows, bl_levels-1))
        ALLOCATE (BL_COEF_KH(1:row_length, 1:rows, bl_levels-1))

        Do k = 1, BL_LEVELS-1
          Do j = 1, rows
            Do i = 1, row_length
              BL_COEF_KM(i,j,k) = 0.0
              BL_COEF_KH(i,j,k) = 0.0
            End Do
          End Do
        End Do

      Else

        ALLOCATE (BL_COEF_KM(1,1,1))
        ALLOCATE (BL_COEF_KH(1,1,1))

      Endif
      if(L_initialise_data)then

! DEPENDS ON: idl_ni_init
        Call IDL_NI_Init(                                               &
     & R, g, kappa, epsilon, Cp, p_zero, Earth_radius, Pi, two_omega    &
     &,model_domain, row_length, rows, n_rows, model_levels,wet_levels  &
     &,TR_VARS, TR_LEVELS, bl_levels, first_constant_r_rho_level        &
     &,cos_theta_latitude, sec_theta_latitude, f3_at_u, f3_at_v         &
     &,timestep, first_atmstep_call, L_regular                          &
     &,delta_x, delta_y, delta_lambda, delta_phi, base_phi, base_lambda &
     &,A_REALHD(rh_rotlat), A_REALHD(rh_rotlong)                        &
     &, glambda_p(lambda_start), phi_p, glambda_u(lambda_start)         &
     &, phi_v, lambda_p_end, phi_p_end                                  &
     &,r_theta_levels, r_rho_levels, r_at_u, r_at_v, z_orog_print       &
     &,eta_theta_levels, eta_rho_levels                                 &
! Multi-processor
     &,offx,offy,halo_i,halo_j, mype, nproc, at_extremity, datastart    &
     &,gc_all_proc_group, global_row_length, global_rows                &
     &,g_rows, g_row_length, g_datastart, nproc_x                       &
! Primary fields
     &,THETA,RHO, EXNER_THETA_LEVELS                                    &
     &,EXNER_RHO_LEVELS,P                                               &
     &,P_THETA_LEVELS, PSTAR                                            &
     &,Q, QCL, QCF,QCF2,QRAIN                                           &
     &,QGRAUP,CF_BULK, CF_LIQUID                                        &
     &,CF_FROZEN,U, V, W                                                &
     &,U_ADV, V_ADV, W_ADV                                              &
! Lateral boundaries
     &,RIMWIDTHA(rima_type_norm), RIMWEIGHTSA                           &
     &,LENRIMA(1,1,rima_type_norm), LBC_SIZEA(1,1,1,rima_type_norm)     &
     &,LBC_STARTA(1,1,1,rima_type_norm)                                 &
     &,THETA_LBC, THETA_LBC_TEND, EXNER_LBC                             &
     &,EXNER_LBC_TEND, RHO_LBC, RHO_LBC_TEND                            &
     &,Q_LBC, Q_LBC_TEND, QCL_LBC, QCL_LBC_TEND                         &
     &,QCF_LBC, QCF_LBC_TEND, QCF2_LBC                                  &
     &,QCF2_LBC_TEND, QRAIN_LBC, QRAIN_LBC_TEND                         &
     &,QGRAUP_LBC, QGRAUP_LBC_TEND, CF_BULK_LBC                         &
     &,CF_BULK_LBC_TEND, CF_LIQUID_LBC                                  &
     &,CF_LIQUID_LBC_TEND, CF_FROZEN_LBC                                &
     &,CF_FROZEN_LBC_TEND, U_LBC, U_LBC_TEND                            &
     &,V_LBC, V_LBC_TEND, W_LBC, W_LBC_TEND                             &
     &,U_ADV_LBC, U_ADV_LBC_TEND, V_ADV_LBC                             &
     &,V_ADV_LBC_TEND, W_ADV_LBC, W_ADV_LBC_TEND                        &
! Grid info for idealised
     &,A_REALHD(rh_z_top_theta), height_domain, big_layers              &
     &,transit_layers, mod_layers, surface_type, p_surface              &
! Profile settings
     &,tprofile_number, qprofile_number, uvprofile_number, Brunt_Vaisala&
     &,theta_surface, dtheta_dz1, height_dz1, u_in, v_in, height_u_in   &
     &,ujet_lat, ujet_width, u_ramp_start, u_ramp_end, f_plane, r_plane &
     &,q1, num_profile_data, zprofile_data, tprofile_data, qprofile_data&
     &,                  num_uvprofile_data, z_uvprofile_data           &
     &,                  uprofile_data, vprofile_data                   &
     &,max_model_levels, max_num_profile_data, max_num_force_times      &
     &,tforce_option, qforce_option, uvforce_option, num_tforce_levels  &
     &,num_tforce_times, num_qforce_levels, num_qforce_times            &
     &,num_uvforce_levels, num_uvforce_times, z_tforce_data, tforce_data&
     &,z_qforce_data, qforce_data, z_uvforce_data, uforce_data          &
     &,vforce_data, tforce_data_modlev, qforce_data_modlev              &
     &,uforce_data_modlev, vforce_data_modlev                           &
! Dynamical core settings
     &,SuHe_pole_equ_deltaT, SuHe_static_stab                           &
     &,base_frictional_timescale, frictional_timescale                  &
     &,SuHe_sigma_cutoff, SuHe_level_weight, L_SH_Williamson            &
!  Horizontal function parameters
     &,                  t_horizfn_number, uv_horizfn_number            &
     &, t_horizfn_data, L_perturb_t, perturb_magnitude_t                &
     &, L_perturb_q, perturb_magnitude_q, L_perturb_correlate_tq        &
     &, L_perturb_correlate_vert, L_perturb_correlate_time              &
     &, perturb_type, perturb_height                                    &
!  Profiles for fixed lbcs and sponge zones
     &,                  u_ref, v_ref, theta_ref, exner_ref, rho_ref    &
     &,                  q_ref                                          &
     &,L_fix_orog_hgt_lbc, orog_hgt_lbc                                 &
     &,zprofile_orog, idl_interp_option, hf                             &
!  Options
     &,L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2                     &
     &,L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc     &
     &,L_constant_dz, L_rotating, L_fixed_lbcs, L_polar_wind_zero       &
     &,L_wind_balance, L_rotate_winds, L_pressure_balance, L_physics    &
     &,L_dry, L_sponge                                                  &
     &,                  L_trivial_trigs, L_perturb, L_code_test        &
     &,L_cyclone, L_baroclinic                                          &
     &,h_print, timestep_number, h_o_actual, grow_steps, h_o_per_step   &
     &,h_o, grid_number, grid_flat, first_theta_height                  &
     &,thin_theta_height, big_factor, mag, lambda_fraction, phi_fraction&
     &,half_width_x, half_width_y, plat_size_x, plat_size_y, Witch_power&
     &,idl_max_num_bubbles, idl_bubble_option, idl_bubble_max           &
     &,idl_bubble_height, idl_bubble_xoffset, idl_bubble_yoffset        &
     &,idl_bubble_width,  idl_bubble_depth, L_idl_bubble_saturate       &
     &,                  a_ixsts,len_a_ixsts, a_spsts,len_a_spsts       &
     &,                  nproc_y, gc_proc_row_group,n_cca_lev           &
     &,IdlSurfFluxSeaOption,IdlSurfFluxSeaParams,L_flux_bc,flux_h,flux_e&
     &,L_spec_z0, z0m_scm, z0h_scm, roughlen_z0m, roughlen_z0h          &
     &,i_hour, i_minute, i_second                                       &
     &,                  problem_number, rad_hr                         &
     &,OROGRAPHY, TSTAR_TILE, ntiles, land_field, land_index            &
     &,CUMULUS, NBDSC, NTDSC, CCA, CCB                                  &
     &,CCT, CCLWP, TSTAR, LAND, SW_INCS                                 &
     &,LW_INCS, T1_SD, Q1_SD, ZH, CF_AREA                               &
     &,TI, Z0, NTML, U_SEA, V_SEA, U_0_P                                &
     &,                  V_0_P )

        elseif(Problem_number  /=   standard                            &
     &                      .and. timestep_number  ==  1 ) Then
          if (.not. L_Physics .and. mype  ==  0)then
           print*,'Data from dump being used without initialising. '
           print*,'If source of dump is a full model run with orography'
           print*,'then there is no guarantee that run will work'
           print*,'since physics is OFF'
          endif !.not. L_Physics.and. mype  ==  0

        end if   ! L_initialise_data
! ---------------------------------------------------------------------
! End Section 0.1  End of Initialisation for idealised test problems
! ---------------------------------------------------------------------
! DEPENDS ON: timer
      If (Ltimer) Call timer('AS Swap_Bounds',5)

! Update domain halos for time-dependent fields
! DEPENDS ON: set_halos
      call set_halos(                                                   &
     &               U, V, W,                                           &
     &               U_ADV, V_ADV, W_ADV,                               &
     &               THETA, Q, QCL, QCF,                                &
     &               QCF2, QRAIN, QGRAUP,                               &
     &               CF_BULK, CF_LIQUID,                                &
     &               CF_FROZEN,                                         &
     &               RHO, P, P_THETA_LEVELS,                            &
     &               EXNER_RHO_LEVELS,                                  &
     &               EXNER_THETA_LEVELS,                                &
     &               MURK,                                              &
     &               DUST_DIV1,DUST_DIV2,                               &
     &               DUST_DIV3,DUST_DIV4,                               &
     &               DUST_DIV5,DUST_DIV6,                               &
     &               SO2, SO4_AITKEN,                                   &
     &               SO4_ACCU, SO4_DISS, DMS,                           &
     &               NH3, SOOT_NEW, SOOT_AGD,                           &
     &               SOOT_CLD, BMASS_NEW,                               &
     &               BMASS_AGD, BMASS_CLD,                              &
     &               OCFF_NEW, OCFF_AGD, OCFF_CLD,                      &
     &               CO2, TRACER, tracer_ukca,                          &
     &               row_length, rows, n_rows, model_levels, wet_levels,&
     &               offx, offy, halo_i, halo_j, tr_levels, tr_vars,    &
     &               tr_ukca,                                           &
     &               L_MURK,                                            &
     &               L_DUST,                                            &
     &               L_SULPC_SO2, L_SULPC_NH3, L_SULPC_DMS,             &
     &               l_soot, l_biomass, l_ocff, l_co2_interactive,      &
     &               L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,             &
     &               OZONE_TRACER,                                      &
     &               L_USE_CARIOLLE                                     &
     &               )


! DEPENDS ON: timer
      If (Ltimer) Call timer('AS Swap_Bounds',6)

! ---------------------------------------------------------------------
!   Section 0.3  Update lbcs for LAMs
! ---------------------------------------------------------------------
      If ((ErrorStatus == 0) .and. (model_domain == mt_lam)) Then
      
       If ( L_lbc_new ) THEN
         L_update_lbcs = .false.
       Else If ( first_atmstep_call ) THEN
         L_update_lbcs = .false.
       Else If ( RIM_STEPSA == 0 ) THEN
         L_update_lbcs = .false.
       Else If ( L_Fixed_lbcs ) THEN
         L_update_lbcs = .false.
       Else !  Old lbcs and NOT first_atmstep_call 
         L_update_lbcs = .true.
       end if ! L_lbc_new

       If ( L_update_lbcs ) Then

        If (MOD (bndary_offsetim(atmos_im) + stepim(atmos_im),          &
     &           RIM_STEPSA) /= 1 ) Then

          If (RIM_STEPSA  /=  1) Then

! DEPENDS ON: BOUNDVAL
              Call BOUNDVAL(LENRIMA(1,1,rima_type_norm),                &
     &                      L_mcr_qcf2_lbc, L_mcr_qrain_lbc,            &
     &                      L_mcr_qgraup_lbc, L_pc2_lbc,                &
     &                      L_murk_lbc, L_int_uvw_lbc,                  &
     &                      U_LBC, U_LBC_TEND,                          &
     &                      V_LBC, V_LBC_TEND,                          &
     &                      W_LBC, W_LBC_TEND,                          &
     &                      RHO_LBC, RHO_LBC_TEND,                      &
     &                      THETA_LBC, THETA_LBC_TEND,                  &
     &                      Q_LBC, Q_LBC_TEND,                          &
     &                      QCL_LBC, QCL_LBC_TEND,                      &
     &                      QCF_LBC, QCF_LBC_TEND,                      &
     &                      QCF2_LBC, QCF2_LBC_TEND,                    &
     &                      QRAIN_LBC, QRAIN_LBC_TEND,                  &
     &                      QGRAUP_LBC, QGRAUP_LBC_TEND,                &
     &                      CF_BULK_LBC, CF_BULK_LBC_TEND,              &
     &                      CF_LIQUID_LBC, CF_LIQUID_LBC_TEND,          &
     &                      CF_FROZEN_LBC, CF_FROZEN_LBC_TEND,          &
     &                      EXNER_LBC, EXNER_LBC_TEND,                  &
     &                      U_ADV_LBC, U_ADV_LBC_TEND,                  &
     &                      V_ADV_LBC, V_ADV_LBC_TEND,                  &
     &                      W_ADV_LBC, W_ADV_LBC_TEND,                  &
     &                      MURK_LBC, MURK_LBC_TEND,                    &
     &                      TRACER_LBC, TRACER_LBC_TEND,                &
     &                      0, 1, ErrorStatus, CMESSAGE)

           End If       ! RIM_STEPSA  /=  1
         End If ! MOD(bndary_offsetim+stepim,RIM_STEPSA) /= 1 )

        End If !  L_update_lbcs 


        !--------------------------------------------------------------
        ! Idealised UM LBC forcing
        !  If active, update lateral boundary arrays to contain
        !  idealised namelist profile data interpolated in time.
        !--------------------------------------------------------------
        If (L_initialise_data .and. L_force_lbc) Then

! DEPENDS ON: idl_force_lbc
          Call IDL_Force_LBC (                                          &
     &             R, g, Cp, kappa, epsilon, p_zero                     &
     &,            row_length, rows, offx, offy                         &
     &,            halo_i, halo_j, Earth_radius                         &
     &,            LENRIMA(1,1,rima_type_norm)                          &
     &,            timestep, timestep_number                            &
     &,            model_levels, wet_levels                             &
     &,            max_model_levels, max_num_force_times                &
     &,            U_LBC, V_LBC                                         &
     &,            THETA_LBC,Q_LBC                                      &
     &,            U_ADV_LBC,V_ADV_LBC                                  &
     &,            EXNER_LBC                                            &
     &,            r_theta_levels, r_rho_levels                         &
     &,            eta_theta_levels, eta_rho_levels                     &
     &,            height_domain, theta_surface                         &
     &,            pforce_option                                        &
     &,            tforce_option, qforce_option, uvforce_option         &
     &,            num_pforce_times                                     &
     &,            num_tforce_times, num_qforce_times                   &
     &,            num_uvforce_times                                    &
     &,            pforce_time_interval                                 &
     &,            tforce_time_interval, qforce_time_interval           &
     &,            uvforce_time_interval                                &
     &,            p_surface_data                                       &
     &,            tforce_data_modlev, qforce_data_modlev               &
     &,            uforce_data_modlev, vforce_data_modlev               &
     &,            newtonian_timescale )

        End If ! on (L_initialise_data .and. L_force_lbc)

       If ( first_atmstep_call ) THEN
         L_apply_lbcs = .true.
       Else If ( L_lbc_new ) THEN
         L_apply_lbcs = .false.
       Else !  Old lbcs  
         L_apply_lbcs = .true.
       end if ! L_lbc_new

       If ( L_apply_lbcs ) THEN

        !--------------------------------------------------------------
        !           Update primary fields with LAM LBC data
        !--------------------------------------------------------------

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UPDATE_LAM_LBCS',3)

! DEPENDS ON: UPDATE_LAM_LBCS
        CALL UPDATE_LAM_LBCS(                                           &
     &    r_rho_levels, r_theta_levels,                                 &
     &    ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,               &
     &    OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY,                         &
     &    L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,                 &
     &    L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc, &
     &    L_murk, L_murk_lbc,                                           &
     &    L_LBC_balance, L_int_uvw_lbc,                                 &
     &    RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                        &
     &    LENRIMA(1,1,rima_type_norm),                                  &
     &    LBC_SIZEA(1,1,1,rima_type_norm),                              &
     &    LBC_STARTA(1,1,1,rima_type_norm),                             &
     &    THETA_LBC,Q_LBC,QCL_LBC,                                      &
     &    QCF_LBC,QCF2_LBC,QRAIN_LBC,                                   &
     &    QGRAUP_LBC, CF_BULK_LBC,CF_LIQUID_LBC,                        &
     &    CF_FROZEN_LBC, RHO_LBC,EXNER_LBC,                             &
     &    U_LBC,V_LBC,W_LBC,                                            &
     &    U_ADV_LBC,V_ADV_LBC,W_ADV_LBC,                                &
     &    MURK_LBC,                                                     &
     &    THETA,Q,QCL,QCF,                                              &
     &    QCF2,QRAIN,QGRAUP,                                            &
     &    CF_BULK,CF_LIQUID,CF_FROZEN,                                  &
     &    RHO,EXNER_RHO_LEVELS,                                         &
     &    U,V,W,                                                        &
     &    U_ADV,V_ADV,W_ADV,                                            &
     &    MURK,                                                         &
     &    DELTA_PHI, DELTA_LAMBDA,                                      &
     &    BASE_PHI, BASE_LAMBDA,                                        &
     &    TWO_OMEGA, DATASTART,                                         &
     &    lat_rot_NP,                                                   &
     &    GLOBAL_ROW_LENGTH, GLOBAL_ROWS                                &
     &     )

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UPDATE_LAM_LBCS',4)
!
! Must now re-calculate the pressure-based variables, namely pressure
! on both rho and theta levels, exner on theta levels and pstar so that
! they are all consistent with the new LBC-updated values of exner on
! rho levels.
!
! DEPENDS ON: consistent_pressure
        call     Consistent_Pressure (                                  &
     &           exner_rho_levels,                                      &
     &           offx,offy,halo_i,halo_J,                               &
     &           row_length,rows,model_levels,                          &
     &           kappa, g, r_theta_levels, r_rho_levels, rho, p_zero,   &
     &           p, pstar, p_theta_levels,exner_theta_levels)
! ----------------------------------------------------------------------
! Now check that cloud is consistent with moisture fields
! Only really needed in lateral boundary zone but done everywhere
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
!  Check liquid/ice cloud fraction is zero if qcl/qcf = 0
! ----------------------------------------------------------------------
!

! DEPENDS ON: cloud_check
        Call Cloud_check(                                               &
     &                 j_start, j_stop, i_start, i_stop                 &
     &,                rows, row_length, wet_levels                     &
     &,                halo_i, halo_j                                   &
     &,                QCL, QCF                                         &
     &,                CF_LIQUID, CF_FROZEN                             &
     &,                CF_AREA, CF_BULK)

        End If ! L_apply_lbcs

        if ( L_int_uvw_lbc ) then

! Obtain tendencies in the boundary zone for u, v, w and use to make.
! lbcs for u_adv, v_adv, w_adv,
! This means that  u_adv, v_adv and w_adv
! can be removed from the lateral boundary files.
! Similar code is needed later for tendencies for the solver

          L_do_halos=.TRUE.
          L_do_boundaries=.TRUE.

          If (RIM_STEPSA  ==  0) Then
            increment_factor=0.0
          Else
            increment_factor=1.0/                                       &
     &      (RIM_STEPSA-MOD(Timestep_Number-1,RIM_STEPSA))
          End If

          lbc_size=LENRIMA(fld_type_p,halo_type_extended,               &
     &             rima_type_norm)

          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              W_LBC_REAL_TEND(i,k) = W_LBC(i,k) +                       &
     &                                   0.5 * increment_factor *       &
     &                           ( W_LBC_TEND(i,k) - W_LBC(i,k) )
            END DO
          END DO

! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH, ROWS, HALO_I, HALO_J,                             &
     &    MODEL_LEVELS, fld_type_p, W_ADV,                              &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    W_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm), RIMWIDTHA(rima_type_norm),         &
     &    RIMWEIGHTSA, AT_EXTREMITY,                                    &
     &    L_do_boundaries, L_do_halos)

          lbc_size=LENRIMA(fld_type_u,halo_type_extended,               &
     &             rima_type_norm)
 
          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              U_LBC_REAL_TEND(i,k) = U_LBC(i,k) +                       &
     &                                      0.5 * increment_factor *    &
     &                          ( U_LBC_TEND(i,k) - U_LBC(i,k) )
            END DO
          END DO

! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH, ROWS, HALO_I, HALO_J,                             &
     &    MODEL_LEVELS, fld_type_u, U_ADV,                              &
     &    LENRIMA(fld_type_u,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_u,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_u,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    U_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm), RIMWIDTHA(rima_type_norm),         &
     &    RIMWEIGHTSA, AT_EXTREMITY,                                    &
     &    L_do_boundaries, L_do_halos)
     
          lbc_size=LENRIMA(fld_type_v,halo_type_extended,               &
     &             rima_type_norm)
 
          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              V_LBC_REAL_TEND(i,k) = V_LBC(i,k) +                       &
     &                                    0.5 * increment_factor *      &
     &                            ( V_LBC_TEND(i,k) - V_LBC(i,k) )
            END DO
          END DO

! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH, N_ROWS, HALO_I, HALO_J,                           &
     &    MODEL_LEVELS, fld_type_v, V_ADV,                              &
     &    LENRIMA(fld_type_v,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_v,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_v,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    V_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm), RIMWIDTHA(rima_type_norm),         &
     &    RIMWEIGHTSA, AT_EXTREMITY,                                    &
     &    L_do_boundaries, L_do_halos)

        endif !  L_int_uvw_lbc
          

      ENDIF     !   model_domain  ==  mt_lam
!QAN - No code for mt_cyclic_lam yet
       
! ----------------------------------------------------------------------
! Section 0.1  Filter winds and theta near poles if active
!              Do horizontal diffusion as a filter if active
! ----------------------------------------------------------------------
!    diagnostic printing of l2norms
      If( L_print_L2norms ) then
        If( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms start of Timestep ***'
        End If ! L_print_pe .or. mype ==0
#include "prfldnorm.h"
      End If !  L_print_L2norms

! Call timer for diffusion code
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('NI_filter_Ctl',3)

      If ( L_filter ) then

! section 13:
        IF( SF(0,13) ) THEN    ! Diagnostics required for this section

! Allocate diagnostic space for STASH
          ALLOCATE (STASHwork13(STASH_maxlen(13,A_im)))
        
        ENDIF

! DEPENDS ON: ni_filter_ctl
       Call NI_filter_Ctl(                                              &
     &                      THETA,                                      &
     &                      U, V, W,                                    &
     &                      EXNER_RHO_LEVELS, RHO,                      &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      r_theta_levels, r_rho_levels,               &
     &                      r_at_u, r_at_v, delta_lambda, delta_phi,    &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      cos_v_latitude, sec_v_latitude,             &
     &                      polar_filter_north_lat_limit,               &
     &                      polar_filter_south_lat_limit,               &
     &                      polar_filter_coefficient,                   &
     &                      polar_filter_n_sweeps,                      &
     &                      polar_filter_step_per_sweep,                &
     &                      polar_filter_lat_limit,                     &
     &                      max_121_rows, u_sweeps, v_sweeps,           &
     &                      global_u_filter, global_v_filter,           &
     &                      u_begin, u_end, v_begin, v_end,             &
     &                      diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
     &                      diff_coeff_thermo, diff_coeff_wind,         &
     &                      diff_order_thermo, diff_order_wind,         &
     &                      first_constant_r_rho_level,                 &
     &                      first_constant_r_rho_level_m1,              &
     &                      top_filt_start, top_filt_end,               &
     &                      up_diff, max_updiff_levels,                 &
     &                      horizontal_level, mype,                     &
     &                      global_row_length, global_rows,             &
     &                      offx, offy, halo_i, halo_j,                 &
     &                      nproc, nproc_x, nproc_y, datastart,         &
     &                      neighbour, at_extremity, model_domain,      &
     &                      gc_proc_row_group, gc_proc_col_group,       &
     &                      L_polar_filter, L_pofil_new,                &
     &                      L_pfcomb, L_pftheta, L_pfuv,                &
     &                      L_pfw, L_pfexner, L_diff_exner,             &
     &                      L_diff_thermo, L_diff_wind, L_diff_w,       &
     &                      L_pofil_hadgem2, Ltimer,exner_theta_levels, &
#include "argsts.h"
     &                      STASHwork13)

      if(L_pfexner .and. L_pofil_new)then
! DEPENDS ON: consistent_pressure
        call     Consistent_Pressure (                                  &
     &           exner_rho_levels,                                      &
     &           offx,offy,halo_i,halo_J,                               &
     &           row_length,rows,model_levels,                          &
     &           kappa, g, r_theta_levels, r_rho_levels, rho, p_zero,   &
     &           p, pstar, p_theta_levels,exner_theta_levels)
      endif  !  (L_pfexner .and. l_pofil_new)

      End If    !  L_filter

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('NI_filter_Ctl',4)

! ----------------------------------------------------------------------
! Section 0.2  If free-slip then put non-zero w field on bottom boundary
! ----------------------------------------------------------------------
      if ( L_free_slip ) then
! DEPENDS ON: bottom_w_calc
        Call Bottom_w_Calc(                                             &
     &                     r_theta_levels,                              &
     &                     u, v,                                        &
     &                     w, w_adv,                                    &
     &                     sec_theta_latitude, delta_lambda, delta_phi, &
     &                     glambda_p(lambda_start), phi_p,              &
     &                     glambda_u(lambda_start), phi_v,              &
     &                     rows, n_rows, row_length, model_levels,      &
     &                     model_domain, L_regular,                     &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     gc_proc_row_group, at_extremity,             &
     &                     global_row_length)

! DEPENDS ON: swap_bounds
        call Swap_Bounds                                                &
     &                  (W,                                             &
     &                   row_length, rows, model_levels+1,              &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
        call Swap_Bounds                                                &
     &                  (W_ADV,                                         &
     &                   row_length, rows, model_levels+1,              &
     &                   halo_i, halo_j, fld_type_p, .false.)

      endif !( L_free_slip )
! Call to the RANDOM PARAMETERS (STPH_RP) subroutine
       IF (L_Physics .and. L_RP) THEN
        minutos=timestep/60
        IF(mod(i_hour,3) == 0 .and. i_minute == minutos) THEN
         IF (PrintStatus  >=   PrStatus_Normal) THEN
            WRITE(6,*) 'CALLING RANDOM PARAMETERS'
         ENDIF
! DEPENDS ON: stph_rp
         CALL stph_rp(max_model_levels,                                 &
     &            RHCRIT,RHCRIT_max,RHCRIT_min,                         &
     &            CI_max,CI_min,                                        &
     &            GWD_FRC,GWD_FRC_max,GWD_FRC_min,                      &
     &            KAY_GWAVE,KAY_GWAVE_max,KAY_GWAVE_min,                &
     &            par_mezcla_max,par_mezcla_min,                        &
     &            G0_max,G0_min,                                        &
     &            G0_RP,par_mezcla)
        ELSE
         IF (PrintStatus  >=  PrStatus_Normal) THEN
             WRITE(6,*) 'NOT CALLING RANDOM PARAMETERS'
             WRITE(6,*) 'This routine is only called every 3hrs'
         ENDIF
        ENDIF
      ENDIF
! Call to the RANDOM PARAMETERS2 (STPH_RP2) subroutine
       IF (L_Physics .and. L_RP2) THEN
        minutos=timestep/60
        IF(mod(i_hour,3) == 0 .and. i_minute == minutos) THEN
         IF (PrintStatus  >=   PrStatus_Normal) THEN
            WRITE(6,*) 'CALLING RANDOM PARAMETERS2'
         ENDIF
! DEPENDS ON: stph_rp2
          CALL stph_rp2(max_model_levels,                               &
     &            RHCRIT,RHCRIT_max,RHCRIT_min,                         &
     &            CI_max,CI_min,                                        &
     &            GWD_FRC,GWD_FRC_max,GWD_FRC_min,                      &
     &            KAY_GWAVE,KAY_GWAVE_max,KAY_GWAVE_min,                &
     &            par_mezcla_max,par_mezcla_min,                        &
     &            G0_max,G0_min,G0_RP,par_mezcla,M_CI,                  &
     &            M_CI_max,M_CI_min,                                    &
     &            Charnock,Charnock_max,Charnock_min)
        ELSE
         IF (PrintStatus  >=  PrStatus_Normal) THEN
             WRITE(6,*) 'NOT CALLING RANDOM PARAMETERS2'
             WRITE(6,*) 'This routine is only called every 3hrs'
         ENDIF
        ENDIF
      ENDIF
! ----------------------------------------------------------------------
!    diagnostic printing of l2norms
      If( L_print_L2norms ) then
        If( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after Filter_ctl ***'
        End If ! L_print_pe .or. mype ==0
#include "prfldnorm.h"
      End If !  L_print_L2norms
! ---------------------------------------------------------------

! ---------------------------------------------------------------
!    diagnostic printing at timestep 1
! ---------------------------------------------------------------
       if( L_diag_print .and. timestep_number ==1 ) then
! DEPENDS ON: print_diag
          Call Print_diag(                                              &
     &                  U, V, THETA,                                    &
     &                  RHO, W,                                         &
     &                  Q, QCL, QCF,                                    &
     &                  rows, n_rows, row_length,                       &
     &                  model_levels, wet_levels, model_domain,         &
     &                  global_row_length, global_rows,                 &
     &                  r_theta_levels, r_rho_levels,                   &
     &                  r_at_u, r_at_v , Pi,                            &
     &                  FV_sec_theta_latitude, FV_cos_theta_latitude,   &
     &                  cos_v_latitude,                                 &
     &                  cos_theta_longitude, sin_theta_longitude,       &
     &                  offx, offy, halo_i, halo_j,                     &
     &                  mype, nproc, at_extremity, datastart,           &
     &                  gc_proc_row_group, delta_lambda, delta_phi,     &
     &                  timestep_number, 1, diag_interval,              &
     &                  rpemax, rpemin, ipesum, rpesum, w_print_limit,  &
     &                  L_print_pe, L_print_w,                          &
     &                  L_print_wmax, L_print_max_wind,                 &
     &                  L_print_div, L_print_lapse, L_print_theta1,     &
     &                  L_print_shear, L_diag_wind, L_diag_noise,       &
     &                  max_w_run, max_wind_run, min_theta1_run,        &
     &                  dtheta1_run, max_div_run, min_div_run,          &
     &                  min_lapse_run, max_shear_run, time_max_shear,   &
     &                  time_div_max, time_div_min, time_lapse_min,     &
     &                  time_w_max, time_max_wind, time_theta1_min,     &
     &                  max_KE_run, min_KE_run, max_noise_run,          &
     &                  time_KE_max, time_KE_min, time_noise_max )
! DEPENDS ON: um_fort_flush
         if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
       endif     !  timestep_number ==1
! ---------------------------------------------------------------
! Section 1.0  Call Atmospheric Physics1
! ----------------------------------------------------------------------
      
      !
      ! Biogenic aerosol climatology for the climate and NWP models
      !
      IF (L_USE_BIOGENIC) THEN 
        do k=1, model_levels
          do j=1, rows
            do i=1, row_length
              biogenic(i,j,k) = ARCLBIOG_BG(i,j,k)
            end do
          end do
        end do
      END IF
            
      !
      ! Aerosol climatologies for the Numerical Weather Prediction
      ! model. Model switches and climatologies are gathered into 
      ! bigger arrays.
      !
      ! First, set the internal model switches according to the
      ! value of the CNTLATM switches, and determine how many 
      ! components we need.
      !
!DEPENDS ON: set_arcl_dimensions
      CALL set_arcl_dimensions(                                         &
     &                         L_USE_ARCLBIOM,                          &
     &                         L_USE_ARCLBLCK,                          &
     &                         L_USE_ARCLSSLT,                          &
     &                         L_USE_ARCLSULP,                          &
     &                         L_USE_ARCLDUST,                          &
     &                         L_USE_ARCLOCFF,                          &
     &                         L_USE_ARCLDLTA,                          &
     &                         n_arcl_species,                          &
     &                         n_arcl_compnts,                          &
     &                         L_USE_ARCL                               &
     &                        )
      
      !
      ! If the aerosol climatology for NWP is used, n_arcl_species
      ! is larger than 0. In that case, allocate the array gathering
      ! component mass-mixing ratio and take the values from the
      ! arrays in arg_atm_fields.h.
      !
      if (n_arcl_species > 0) Then
      
        allocate(arcl(row_length, rows, model_levels, n_arcl_compnts))
      
!DEPENDS ON: set_arcl_clim
        call set_arcl_clim(                                             &
                           ! Array dimensions
     &                     row_length, rows, model_levels,              &
     &                     n_arcl_compnts,                              &
                           ! Internal model switches
     &                     L_USE_ARCL,                                  &
                           ! Climatologies from ancillary files
     &                     arclbiom_fr, arclbiom_ag, arclbiom_ic,       &
     &                     arclblck_fr, arclblck_ag,                    &
     &                     arclsslt_fi, arclsslt_jt,                    &
     &                     arclsulp_ac, arclsulp_ak, arclsulp_di,       &
     &                     arcldust_b1, arcldust_b2, arcldust_b3,       &
     &                     arcldust_b4, arcldust_b5, arcldust_b6,       &
     &                     arclocff_fr, arclocff_ag, arclocff_ic,       &
     &                     arcldlta_dl,                                 &
                           ! Internal climatology array
     &                     arcl,                                        &
                           ! Component array indices
     &                     i_arcl_compnts                               &
     &       )
        
      else
         allocate ( arcl(1,1,1,1) )
      end if
            
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS1 Atmos_Phys1',5)

      IF (L_Physics .AND. ErrorStatus == 0) THEN

! Set logicals for tropopause diagnostics

        L_O3_trop_level  = ( sf(280,2) ) .AND.                          &
     &                     ( ( i_ozone_int  ==  IO3_TROP_MAP) .OR.      &
     &                       ( i_ozone_int  ==  IO3_TROP_MAP_MASSCON) )
        L_O3_trop_height = ( sf(281,2) ) .AND.                          &
     &                     ( ( i_ozone_int  ==  IO3_TROP_MAP) .OR.      &
     &                       ( i_ozone_int  ==  IO3_TROP_MAP_MASSCON) )
        L_T_trop_level   = ( sf(282,2) ) .AND.                          &
     &                     ( ( i_ozone_int  ==  IO3_TROP_MAP) .OR.      &
     &                       ( i_ozone_int  ==  IO3_TROP_MAP_MASSCON) )
        L_T_trop_height  = ( sf(283,2) ) .AND.                          &
     &                     ( ( i_ozone_int  ==  IO3_TROP_MAP) .OR.      &
     &                       ( i_ozone_int  ==  IO3_TROP_MAP_MASSCON) )

! Allocate space for tropopause diagnostics
        IF( L_O3_trop_level ) THEN
          ALLOCATE ( O3_trop_level(row_length,rows) )
        ELSE
          ALLOCATE ( O3_trop_level(1,1) )
        ENDIF
        IF( L_O3_trop_height ) THEN
          ALLOCATE ( O3_trop_height(row_length,rows) )
        ELSE
          ALLOCATE ( O3_trop_height(1,1) )
        ENDIF
        IF( L_T_trop_level ) THEN
          ALLOCATE ( T_trop_level(row_length,rows) )
        ELSE
          ALLOCATE ( T_trop_level(1,1) )
        ENDIF
        IF( L_T_trop_height ) THEN
          ALLOCATE ( T_trop_height(row_length,rows) )
        ELSE
          ALLOCATE ( T_trop_height(1,1) )
        ENDIF

! Check whether ozone option is applicable to the specified
! ozone ancillary file. If not output error message.
        IF (PrintStatus  >=   PrStatus_Diag) THEN   
            WRITE(6,*) 'Atm Step: Lexpand_ozone', Lexpand_ozone
        END IF
  
        IF ((Lexpand_ozone).and.(I_ozone_int == 1)) THEN
           ErrorStatus=123
! DEPENDS ON: ereport
           Call Ereport("ATM_STEP", ErrorStatus,                        &
     &        "A 2D ozone ancillary has been specified with a" //       &
     &        "3D ozone option.")
        ENDIF

        IF (.not.Lexpand_ozone) THEN
           IF (I_ozone_int == 2) THEN
              ErrorStatus=123
! DEPENDS ON: ereport
              Call Ereport("ATM_STEP", ErrorStatus,                     &
     &           "A 3D ozone ancillary has been specified with a" //    &
     &           "2D ozone option.")
           ENDIF

           IF (I_ozone_int == 5) THEN
              ErrorStatus=123
! DEPENDS ON: ereport
              Call Ereport("ATM_STEP", ErrorStatus,                     &
     &           "A 3D ozone ancillary has been specified with a" //    &
     &           "2D ozone option.")
           ENDIF
        ENDIF


!     Convert the ozone field supplied from the dump to a full 3-D
!     array. At 5.2 we require the ozone field to be specified on
!     the same number of levels as the model's grid: this will later
!     be generalized to allow for a different number of levels.
        If (Lexpand_ozone) then
          nd_o3=rows*ozone_levels
        Else
          nd_o3=row_length*rows*ozone_levels
        Endif
!  STOCHEM size
        nd_stochem=row_length*rows*model_levels

! Allocate workspace, just required in atmos_physics1
        ALLOCATE (ozone3D (row_length, rows, ozone_levels) )
! DEPENDS ON: o3_to_3d
        Call O3_to_3D(                                                  &
     &                lexpand_ozone, i_ozone_int,                       &
     &                rows, row_length, model_levels, ozone_levels,     &
     &                halo_i, halo_j, offx, offy, at_extremity,         &
     &                a_realhd(rh_z_top_theta),                         &
     &                THETA,                                            &
     &                r_theta_levels, eta_theta_levels,                 &
     &                EXNER_THETA_LEVELS,                               &
     &                RHO,                                              &
     &                r_rho_levels, eta_rho_levels,                     &
     &                EXNER_RHO_LEVELS,                                 &
     &                nd_o3, O3,                                        &
     &                l_use_stochem_o3,                                 &
     &                nd_stochem,O3_STOCH,                              &
     &  min_trop_level, max_trop_level,                                 &
     &  L_O3_trop_level,L_O3_trop_height,                               &
     &  L_T_trop_level,L_T_trop_height,                                 &
     &  O3_trop_level,O3_trop_height,                                   &
     &  T_trop_level,T_trop_height,                                     &
     &  gc_proc_row_group,                                              &
     &  global_row_length,                                              &
     &                ozone3D,                                          &
     &                ErrorStatus, cmessage                             &
     &                )
!
!     If the ozone_tracer is initialised to 0 then it should be reset to something realistic
!     so set it to climatology expanded to 3D
     
       IF (L_USE_CARIOLLE) THEN
          IF (OZONE_TRACER(1,1,1) == 0.0) THEN
             WRITE(6,*)'O3 tracer must not be set to 0 reset to clim'
        
             DO k=1,model_levels
               DO j=1,rows
                 DO i=1,row_length        
                   OZONE_TRACER(i,j,k)=ozone3D(i,j,k)      
                 END DO
               END DO
             END DO
         
! Halos updated
! DEPENDS ON: swap_bounds
             call Swap_Bounds(OZONE_TRACER,                             &
                        row_length, rows, model_levels,                 &
                        offx, offy, fld_type_p, .false.)
!DEPENDS ON: fill_external_halos
             CALL FILL_EXTERNAL_HALOS(OZONE_TRACER,row_length, rows,    &
                          model_levels,offx,offy)
          
          ELSE  
            WRITE(6,*) 'At least first row of Ozone tracer is not zero'
          END IF ! End of check for zero ozone
       END IF
!     After a successful call ozone3D will contain the mixing ratios of
!     ozone on the model's grid. Uses ozone tracer from the last timestep
!     in the radiation scheme. Radiation is done first so all variables used 
!     are from the last timestep.
!    

       IF (L_USE_OZONEINRAD) THEN                           
             DO k = 1, model_levels                           
               DO j = 1, rows                             
                 DO i = 1, row_length                                  
                   temp_ozone(i,j,k)=OZONE_TRACER(i,j,k)

!  Use the 3D ozone in radiation but check for zero or -ve ozone as some 
!  tropospheric analysed ozone may be -ve; do not use -ve values.
                     IF (temp_ozone(i,j,k) > 0.0) THEN   
                        ozone3D(i,j,k) = temp_ozone(i,j,k)
                     ELSE   
                       IF (j > 1 .and. j < rows .and.     &
     &                    i > 1 .and. i < row_length)       THEN
                         IF (OZONE_TRACER(i,j-1,k) > 0.0)   THEN
                           ozone3D(i,j,k) = OZONE_TRACER(i,j-1,k)
                         ELSE
                           IF (OZONE_TRACER(i,j+1,k) > 0.0) THEN
                             ozone3D(i,j,k) = OZONE_TRACER(i,j+1,k)
                           ELSE
                             IF (OZONE_TRACER(i-1,j,k) > 0.0)    THEN
                              ozone3D(i,j,k) = OZONE_TRACER(i-1,j,k)
                             ELSE
                                IF (OZONE_TRACER(i+1,j,k) > 0.0) THEN
                                   ozone3D(i,j,k) = OZONE_TRACER(i+1,j,k)
                                END IF   ! Check on i+1>0
                             END IF    ! Check on i-1>0
                           END IF    ! Check on j+1>0
                         END IF    ! Check on j-1>0
                       END IF    ! Check on j (1:rows) and i (1:row_length)
                     END IF    ! Check on temp_ozone > 0                                          
                 END DO                                    
               END DO                              
             END DO

       END IF                                                    
      
      END IF !  L_Physics
!
! Set LAND_PTS_TRIF and NPFT_TRIF according to TRIFFID on/off
! and set pointers to single or multi pool soil carbon, and
! soil carbon dimensions to be used in subroutines.
      IF (L_TRIFFID) THEN
        LAND_PTS_TRIF = LAND_FIELD
        NPFT_TRIF = NPFT
        CS => SOIL_CARB1
        RSA => RSP_S_ACC1
        DIM_CS1 = 4
        DIM_CS2 = LAND_FIELD
      ELSE
        LAND_PTS_TRIF = 1
        NPFT_TRIF = 1
        CS => SOIL_CARB
        RSA => RSP_S_ACC
        DIM_CS1 = 1
        DIM_CS2 = 1
      ENDIF
!     dim_cs2 needs to be modified for Carbon fluxes
      IF ( l_cable ) DIM_CS2 = LAND_FIELD
     ! Lestevens 23apr13
      IF ( l_cable ) LAND_PTS_TRIF = LAND_FIELD
      IF ( l_cable ) NPFT_TRIF = NPFT
!
!  set up CO2 field to be passed down
!
      IF (L_CO2_INTERACTIVE) THEN
        CO2_DIM_LEN = row_length
        CO2_DIM_ROW = rows
        CO2_DIM_LEV = model_levels
      ELSE
        CO2_DIM_LEN = 1
        CO2_DIM_ROW = 1
        CO2_DIM_LEV = 1
      ENDIF

      if (L_use_seasalt_autoconv .OR. L_use_seasalt_sulpc .OR.          &
     &         L_use_seasalt_indirect .OR. L_use_seasalt_direct) then
        salt_dim1=row_length
        salt_dim2=rows
        salt_dim3=model_levels
      else
        salt_dim1=1
        salt_dim2=1
        salt_dim3=1
      endif

      if (L_use_seasalt_sulpc .OR. L_DMS_em_inter                       &
     &                        .OR. L_DMS_Ointer) then
        aero_dim1=row_length
        aero_dim2=rows
      else
        aero_dim1=1
        aero_dim2=1
      endif

      if (L_use_seasalt_sulpc) then
        aero_dim3=model_levels
      else
        aero_dim3=1
      endif

! Allocate additional microphysics variables to full size
! if in use, otherwise allocate minimum amount of space

      If (L_mcr_qcf2) Then  ! Allocate second cloud ice
        Allocate ( qcf2_star(1-offx:row_length+offx,                    &
     &                       1-offy:rows+offy, wet_levels) )
      Else
        Allocate ( qcf2_star(1,1,1) )
      End If

      If (L_mcr_qrain) Then  ! Allocate rain
        Allocate ( qrain_star(1-offx:row_length+offx,                   &
     &                        1-offy:rows+offy, wet_levels) )
      Else
        Allocate ( qrain_star(1,1,1) )
      End If

      If (L_mcr_qgraup) Then  ! Allocate graupel
        Allocate ( qgraup_star(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
      Else
        Allocate ( qgraup_star(1,1,1) )
      End If

! The _star fields are used to store the increments to theta, q, qcl,
! and qcf
      IF (L_Physics .AND. ErrorStatus == 0) THEN

! Allocate diagnostic space for STASH
      If (L_radiation) then
        ALLOCATE (STASHwork1(STASH_maxlen(1,A_im)))
        ALLOCATE (STASHwork2(STASH_maxlen(2,A_im)))
      endif
      If (L_rain) then
        ALLOCATE (STASHwork4(STASH_maxlen(4,A_im)))
      endif
      IF(L_gwd .or. l_use_ussp)then
        ALLOCATE (STASHwork6(STASH_maxlen(6,A_im)))
      endif
        IF (L_EMCORR) THEN   ! only if energy correction required
          ALLOCATE (STASHwork14(STASH_maxlen(14,A_im)))
        ENDIF   ! l_emcorr

      ! -----------------------------------------------------------
      ! Convert to mixing ratios from specific humidities if needed
      ! -----------------------------------------------------------
      If (l_mr_physics1) Then

        ! Allocate q_store variables

        allocate ( q_store  (1-halo_i:row_length+halo_i,                &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        allocate ( qcl_store(1-halo_i:row_length+halo_i,                &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        allocate ( qcf_store(1-halo_i:row_length+halo_i,                &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        If (L_mcr_qcf2) then
          allocate ( qcf2_store(1-halo_i:row_length+halo_i,             &
     &                       1-halo_j:rows+halo_j, wet_levels) )
        Else
          allocate ( qcf2_store(1,1,1) )
        End if
        If (L_mcr_qrain) then
          allocate ( qrain_store(1-halo_i:row_length+halo_i,            &
     &                        1-halo_j:rows+halo_j, wet_levels) )
        Else
          allocate ( qrain_store(1,1,1) )
        End if
        If (L_mcr_qgraup) then
          allocate ( qgraup_store(1-halo_i:row_length+halo_i,           &
     &                         1-halo_j:rows+halo_j, wet_levels) )
        Else
          allocate ( qgraup_store(1,1,1) )
        End if

        ! Hold a copy of d1(q) etc. in q_store to place back in d1
        ! after the conversions are finished

        Do k = 1, wet_levels
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              q_store(i,j,k)   = q(i,j,k)   ! Vapour
              qcl_store(i,j,k) = qcl(i,j,k) ! Liquid
              qcf_store(i,j,k) = qcf(i,j,k) ! Ice
            End do
          End do
        End do

        If (L_mcr_qcf2) then  ! Ice2
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qcf2_store(i,j,k) = qcf2(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qrain) then  ! Rain
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qrain_store(i,j,k) = qrain(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qgraup) then  ! Graupel
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qgraup_store(i,j,k) = qgraup(i,j,k)
              End do
            End do
          End do
        End if

        ! Convert d1 values to mixing ratios, including the halos,
        ! using the _store values as inputs (d1s are outputs)

! DEPENDS ON: q_to_mix_halo
        call q_to_mix_halo (row_length, rows, wet_levels,               &
     &           halo_i, halo_j,                                        &
     &           q_store, qcl_store, qcf_store,                         &
     &           qcf2_store, qrain_store, qgraup_store,                 &
     &           L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                 &
     &           Q, QCL, QCF,                                           &
     &           QCF2, QRAIN, QGRAUP                                    &
     &          )

        ! d1 values now contain mixing ratios (including halos).
        ! _store values are specific humidities.

      ! -----------------------------------------------------------
      ! End of conversion to mixing ratios for atmos_physics1 call
      ! -----------------------------------------------------------
      End if  ! l_mr_physics1

        If (L_pc2) Then
!         Calculate the intital Relative humidity wrt TL (rhts)
          allocate(rhts(row_length,rows,wet_levels))
! DEPENDS ON: pc2_rhtl
          Call pc2_rhtl(halo_i, halo_j, offx, offy                      &
     &,     wet_levels,row_length,rows                                  &
     &,     THETA, EXNER_THETA_LEVELS                                   &
     &,     Q,QCL,P_THETA_LEVELS,rhts                                   &
     &,     l_mr_physics1)
! Cusack interpolation needs RH on large_levels from start of timestep
          allocate(tlts(row_length,rows,wet_levels))
          allocate(qtts(row_length,rows,wet_levels))
          allocate(ptts(row_length,rows,wet_levels))
	  do k=1,wet_levels
            do j=1,rows
              do i=1,row_length
                 tlts(i,j,k)=theta(i,j,k)*exner_theta_levels(i,j,k)     &
     &              -QCL(i,j,k)*LC/CP
                 ptts(i,j,k)=p_theta_levels(i,j,k)
                 qtts(i,j,k)=Q(i,j,k)+QCL(i,j,k)
              end do
            end do
          end do
        EndIf

! NB: the star variables and R_u and R_v have not been set in the
!     halo region yet.

! DEPENDS ON: timer
        If (Ltimer) Call timer('Atmos_Physics1',3)

           do l=1,land_field
              frac_control(l,1)=FRAC_CON1(l)
              frac_control(l,2)=FRAC_CON2(l)
              frac_control(l,3)=FRAC_CON3(l)
              frac_control(l,4)=FRAC_CON4(l)
              frac_control(l,5)=FRAC_CON5(l)
              frac_control(l,6)=FRAC_CON6(l)
              frac_control(l,7)=FRAC_CON7(l)
              frac_control(l,8)=FRAC_CON8(l)
              frac_control(l,9)=FRAC_CON9(l)
           enddo
          
      istep_cur = istep_cur + 1  ! For CABLE
! NB if you are changing the argument list to atmos_physics1, please
! do an equivalent change in routine scm_main to keep the single column
! model consistent.

   CALL iday_kick(i_day_number, cable_lai)

! DEPENDS ON: atmos_physics1
        Call Atmos_Physics1(                                            &
! Parallel variables
     & halo_i, halo_j, offx, offy, global_row_length, global_rows       &
     &,gc_proc_row_group,gc_proc_col_group, at_extremity, nproc,nproc_x &
     &,nproc_y, neighbour, g_rows, g_row_length, g_datastart, mype      &
! field dimensions etc.
     &,row_length, rows, n_rows, land_field, model_levels, wet_levels   &
     &,bl_levels, st_levels, sm_levels, Ozone_levels, cloud_levels      &
     &,land_ice_points,soil_points,n_cca_lev,ntiles,salt_dim1,salt_dim2 &
     &,salt_dim3,tr_levels,tr_vars,co2_dim_len,co2_dim_row,co2_dim_lev  &
     &,n_arcl_species, n_arcl_compnts, i_arcl_compnts                   &
! model switches
     &,     model_domain, L_regular, L_SEC_VAR, L_EqT                   &
     &,     L_Rad_Step, L_Rad_Step_diag, L_Rad_Step_prog                &
     &,     L_Forcing, L_Timestep, L_Radiance, L_Wenyi                  &
     &,     LCAL360, L_MICROPHY, L_emcorr, L_climat_aerosol, Ltimer     &
     &,     L_gwd, L_use_ussp, l_taus_scale, l_fix_gwsatn, l_gwd_40km   &
     &,     l_ussp_opaque, sat_scheme, l_use_clearrh, L_ssice_albedo    &
     &,     L_RHCPT, L_Murk, l_murk_source, l_murk_bdry, L_MURK_RAD     &
     &,     L_DUST, L_SULPC_SO2, L_SULPC_NH3, L_SOOT, L_BIOMASS, L_OCFF &
! kdcorbin, 06/10 - added L_CO2_RADIATION flag
     &,l_co2_interactive,l_co2_radiation                                &
     &,Lflux_reset,L_clim_aero_hgt,L_HadGEM1_Clim_Aero                  &
     &,L_use_dust, L_use_sulphate_autoconv, L_auto_debias,              &
     &L_use_seasalt_autoconv,L_use_seasalt_indirect,L_use_seasalt_direct&
     &,     L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a          &
     &,     L_snow_albedo, L_ctile, L_radiation, L_rain                 &
     &,     L_INHOM_CLOUD, L_USE_BIOGENIC                               &
     &,     L_USE_SULPC_DIRECT, l_use_soot_direct, l_use_soot_indirect  &
     &,     l_use_soot_autoconv, l_use_bmass_direct                     &
     &,     l_use_bmass_indirect, l_use_bmass_autoconv                  &
     &,     l_use_ocff_direct, l_use_ocff_indirect, l_use_ocff_autoconv &
     &,     L_USE_SULPC_INDIRECT_SW, L_USE_SULPC_INDIRECT_LW            &
     &,     L_pc2, L_eacf, L_mr_physics1,l_cry_agg_dep,l_droplet_settle &
     &,     L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup                       &
     &,     L_USE_METHOX, L_rad_deg, L_TRIFFID                          &
     &,     l_use_stochem_ch4, L_it_melting, L_ukca, L_USE_ARCL         &
     &,     L_USE_SPEC_SEA, L_MOD_BARKER_ALBEDO, L_MOD_K_FLUX           &
     &,     L_SCVARY, L_VOLCTS                                          &
! model Parameters
     &,     RHcrit, cw_sea, cw_land                                     &
     &,     A_SW_segments, A_SW_seg_size, A_LW_segments, A_LW_seg_size  &
     &,     A_SW_radstep, A_LW_radstep, A_SW_radstep_diag               &
     &,     A_SW_radstep_prog, A_LW_radstep_diag, A_LW_radstep_prog     &
     &,     aero_bl_levels,INHOM_CLOUD_SW, INHOM_CLOUD_LW               &
     &,     DP_CORR_STRAT, DP_CORR_CONV, CO2_MMR, alpham,alphac,alphab  &
     &,     dtice,dt_bare,dalb_bare_wet,pen_rad_frac,SW_beta            &
     &,min_trop_level, max_trop_level, KAY_GWAVE, GWD_FRC, gwd_fsat     &
     &,O2MMR, N2OMMR                                                    &
     &,CH4MMR, C11MMR, C12MMR, C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR&
! Pass down position of greenhouse gases in free_tracers
!                array, for chemical coupling
     &,     ngrgas ,grgas_addr                                          &
     &,     cloud_fraction_method,overlap_ice_liquid                    &
     &,     ice_fraction_method,ctt_weight,t_weight                     &
     &,     qsat_fixed,sub_cld,dbsdtbs_turb_0                           &
     &,     L_seq_mcr,L_autoc_3b,L_autolim_3b,L_autoconv_murk           &
     &,     ec_auto,N_drop_land                                         &
     &,     N_drop_sea,N_drop_land_cr,N_drop_sea_cr,Ntot_land, Ntot_sea &
     &,     x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic                    &
     &,     lsp_ei,lsp_fi,lsp_eic,lsp_fic                               &
! parameter for stochastic physics random parameters2
     &,     M_CI                                                        &
! Physical constants
     &,Lc,Lf,Cp, two_Omega, p_zero, kappa,R, g, Lapse, earth_radius, Pi &
! Vertical coordinate levels.
     &, r_rho_levels, r_theta_levels, r_at_u, r_at_v, eta_theta_levels  &
     &, eta_rho_levels, delta_lambda, delta_phi, lat_rot_NP, long_rot_NP&
! Time stepping information
     &,     timestep, radiation_timestep, radiation_tstep_diag          &
     &,     radiation_tstep_prog, I_year, I_day_number, I_hour          &
     &,     I_minute, I_second, timestep_number, PREVIOUS_TIME          &
!glr cable 1 line
     &, istep_cur                                                       &
! trig arrays
     &,     sin_theta_longitude, cos_theta_longitude                    &
     &,     FV_cos_theta_latitude, sin_theta_latitude                   &
! grid-dependent arrays
     &,     f3_at_u, true_longitude, true_latitude,                     &
! diagnostic info
#include "argsts.h"
     &      STASHwork1,STASHwork2,STASHwork4,STASHwork6,STASHwork14     &
!
! Additional variables for SCM diagnostics
     &,     nSCMDpkgs, L_SCMDiags                                       &
!
! Data Fields.
     &,     THETA, Q, QCL, QCF, QCF2                                    &
     &,     QRAIN, QGRAUP, RHO, U, V                                    &
     &,     P, PSTAR, EXNER_RHO_LEVELS                                  &
     &,     EXNER_THETA_LEVELS, LAND                                    &
     &,     P_THETA_LEVELS, FRAC_LAND,frac_control                      &

! ancillary fields and fields needed to be kept from timestep to
! timestep
     &,land_index, RGRAIN_TILE,SNSOOT,NTML, CUMULUS                     &
     &,ICE_FRACTION, CCA, CCB, CCT, CCLWP, CCW_RAD, LCBASE              &
     &,TSTAR,TSTAR_LAND,TSTAR_SEA, TSTAR_SICE                           &
     &,SICE_ALB, LAND_ALB, SNODEP,SNODEP_SEA                            &
     &,ozone3D,SW_INCS                                                  &
     &, LW_INCS, DIRPAR                                                 &
     &,     O3_trop_level, O3_trop_height, T_trop_level, T_trop_height  &
     &,ZH, OROG_SD, OROG_GRAD_XX, OROG_GRAD_XY                          &
     &,OROG_GRAD_YY, CF_AREA, CF_BULK                                   &
     &,CF_LIQUID,CF_FROZEN,MURK_SOURCE,arcl                             &
     &,SOIL_ALB,LAI_PFT,SNODEP_TILE, FRAC_TYP                           &
     &,TSTAR_TILE,Z0_TILE,DOLR_FIELD,LW_DOWN,SW_TILE_RTS                &
     &,ES_SPACE_INTERP, RAD_MASK, CH4_STOCH                             &
     &, cos_zenith_angle,can_rad_mod                                    &
! IN/OUT
     &,theta_star, q_star, qcl_star, qcf_star, qcf2_star, qrain_star    &
     &,qgraup_star, cf_star, cfl_star, cff_star, R_u, R_v               &
     &,a_realhd(rh_energy_corr), NET_FLUX, NET_MFLUX                    &
     &,MURK, DUST_DIV1, DUST_DIV2                                       &
     &,DUST_DIV3, DUST_DIV4, DUST_DIV5                                  &
     &,DUST_DIV6, SO2, SO4_AITKEN                                       &
     &,SO4_ACCU, SO4_DISS, NH3, SOOT_NEW                                &
     &,SOOT_AGD, SOOT_CLD, BMASS_NEW                                    &
     &,BMASS_AGD,BMASS_CLD, OCFF_NEW, OCFF_AGD, OCFF_CLD, CO2, TRACER   &
     &,biogenic, A_INTHD(23)                                            &
!glr cable 
     &, surf_down_sw,alb_tile                                           &
     &, day,TILE_PTS,SM_LEVELS,TILE_INDEX                               &
     &, SNOW_TMP3L,SNOW_RHO1L,TSOIL_TILE,SNOW_FLG3L,LAND_ALBEDO_CABLE   &
     &, l_cable                                                         &
! OUT Fields
     &,     ls_rain, ls_snow, micro_tends, unscaled_dry_rho             &
     &,     photosynth_act_rad, rad_hr, surf_radflux, dolr, SW_tile     &
! Section information
     &,     maxsects,h_sect                                             &
! error information
     &,     ErrorStatus  )

! DEPENDS ON: timer
        If (Ltimer) Call timer('Atmos_Physics1',4)
      
      ! Deallocate the array of the aerosol climatology for NWP
      deallocate(arcl)
      
      ! --------------------------------------------------------------
      ! Convert back from  mixing ratios to spec. humidities if needed
      ! --------------------------------------------------------------
      If (l_mr_physics1) Then

        ! Allocate mix_star variables to hold mixing ratio increments
        ! from atmos_physics1

        allocate ( mix_v_star (1-offx:row_length+offx,                  &
     &                   1-offy:rows+offy, wet_levels) )
        allocate ( mix_cl_star(1-offx:row_length+offx,                  &
     &                   1-offy:rows+offy, wet_levels) )
        allocate ( mix_cf_star(1-offx:row_length+offx,                  &
     &                   1-offy:rows+offy, wet_levels) )

        If (L_mcr_qcf2) then
          allocate ( mix_cf2_star(1-offx:row_length+offx,               &
     &                   1-offy:rows+offy, wet_levels) )
        Else
          allocate ( mix_cf2_star(1,1,1) )
        End if

        If (L_mcr_qrain) then
          allocate ( mix_rain_star(1-offx:row_length+offx,              &
     &                   1-offy:rows+offy, wet_levels) )
        Else
          allocate ( mix_rain_star(1,1,1) )
        End if

        If (L_mcr_qgraup) then
          allocate ( mix_graup_star(1-offx:row_length+offx,             &
     &                   1-offy:rows+offy, wet_levels) )
        Else
          allocate ( mix_graup_star(1,1,1) )
        End if

        ! q_star currently holds mixing ratio increments from physics1
        ! but we wish it to hold specific humidity increments.
        ! Start by copying q_star into mix_v_star etc.

        Do k = 1, wet_levels
          Do j = 1-offy, rows+offy
            Do i = 1-offx, row_length+offx
              mix_v_star(i,j,k)  = q_star(i,j,k)     ! Vapour
              mix_cl_star(i,j,k) = qcl_star(i,j,k)   ! Liquid
              mix_cf_star(i,j,k) = qcf_star(i,j,k)   ! Ice
            End do
          End do
        End do

        If (L_mcr_qcf2) then  ! Ice2
          Do k = 1, wet_levels
            Do j = 1-offy, rows+offy
              Do i = 1-offx, row_length+offx
                mix_cf2_star(i,j,k) = qcf2_star(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qrain) then  ! Rain
          Do k = 1, wet_levels
            Do j = 1-offy, rows+offy
              Do i = 1-offx, row_length+offx
                mix_rain_star(i,j,k) = qrain_star(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qgraup) then  ! Graupel
          Do k = 1, wet_levels
            Do j = 1-offy, rows+offy
              Do i = 1-offx, row_length+offx
                mix_graup_star = qgraup_star(i,j,k)
              End do
            End do
          End do
        End if

        ! Now convert the mixing ratio increments (stored in
        ! mix_v_star etc.) to specific humidity increments
        ! and write these s.h. incs. back into the q_star variables.
        ! d1 variables currently contain the mixing ratios before
        ! atmos_physics1.
        ! q_store variables currently contain specific humidities
        ! before atmos_physics1.

! DEPENDS ON: calc_q_star
        call calc_q_star (row_length, rows, wet_levels,                 &
     &           halo_i, halo_j, offx, offy,                            &
     &           Q, QCL, QCF,                                           &
     &           QCF2, QRAIN, QGRAUP,                                   &
     &           mix_v_star, mix_cl_star, mix_cf_star,                  &
     &           mix_cf2_star, mix_rain_star, mix_graup_star,           &
     &           L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                 &
     &           q_store, qcl_store, qcf_store,                         &
     &           qcf2_store, qrain_store, qgraup_store,                 &
     &           q_star, qcl_star, qcf_star,                            &
     &           qcf2_star, qrain_star, qgraup_star                     &
     &                   )

        ! Copy contents of q_store (specific humidities before
        ! atmos_physics1) back into the d1(q) variables etc. Include
        ! the halo points.

        Do k = 1, wet_levels
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              q(i,j,k)   = q_store(i,j,k)    ! Vapour
              qcl(i,j,k) = q_store(i,j,k)    ! Liquid
              qcf(i,j,k) = q_store(i,j,k)    ! Ice
            End do
          End do
        End do

        If (L_mcr_qcf2) then  ! Ice2
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qcf2(i,j,k) = qcf2_store(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qrain) then  ! Rain
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qrain(i,j,k) = qrain_store(i,j,k)
              End do
            End do
          End do
        End if

        If (L_mcr_qgraup) then  ! Graupel
          Do k = 1, wet_levels
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                qgraup(i,j,k) = qgraup_store(i,j,k)
              End do
            End do
          End do
        End if

        ! At this point all the conversions are complete.
        ! d1(q) again contains the specific humidities before the
        ! atmos_physics1 call. q_star contains specific humidity
        ! increments from atmos_physics1. However, the actual
        ! calculations that were done within atmos_physics1 were done
        ! by using mixing ration variables.
        !
        !
        ! Finally, Deallocate the temporary variables

        deallocate(q_store)
        deallocate(qcl_store)
        deallocate(qcf_store)
        deallocate(qcf2_store)
        deallocate(qrain_store)
        deallocate(qgraup_store)
        deallocate(mix_v_star)
        deallocate(mix_cl_star)
        deallocate(mix_cf_star)
        deallocate(mix_cf2_star)
        deallocate(mix_rain_star)
        deallocate(mix_graup_star)

      ! -----------------------------------------------------------
      ! End of conversion to back to specific humidities
      ! -----------------------------------------------------------
      End if  ! l_mr_physics1

        ! Clear workspace for radiation
        DEALLOCATE ( ozone3D )
! Deallocate space for tropopause diagnostics
        DEALLOCATE (O3_trop_level)
        DEALLOCATE (O3_trop_height)
        DEALLOCATE (T_trop_level)
        DEALLOCATE (T_trop_height)

! Diagnostics STASHed for each section in Atmos_Physics1:
      If (L_radiation) then
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,1,STASHwork1,                              &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',4)
        DEALLOCATE (STASHwork1)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,2,STASHwork2,                              &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',4)
        DEALLOCATE (STASHwork2)
      Endif

      If (L_rain) then
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,4,STASHwork4,                              &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',4)
        DEALLOCATE (STASHwork4)
      Endif


      IF(L_gwd .or. l_use_ussp)then
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,6,STASHwork6,                              &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',4)
        DEALLOCATE (STASHwork6)
      endif

        IF(L_EMCORR) THEN
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,14,STASHwork14,                          &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork14)
        ENDIF   ! l_emcorr

      Else  ! L_physics =.false.
! initialise arrays that hold physics increments to zero, only needs
! doing at non-halo points
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              theta_star(i,j,k) = 0.0
            End Do
          End Do
        End Do
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              q_star(i,j,k) = 0.0
              qcl_star(i,j,k) = 0.0
              qcf_star(i,j,k) = 0.0
              cf_star(i,j,k) = 0.0
              cfl_star(i,j,k) = 0.0
              cff_star(i,j,k) = 0.0
            End Do
          End Do
        End Do
        ! Initialise additional microphysics variables if in use
        If (L_mcr_qcf2)   qcf2_star(1:row_length, 1:rows, :) = 0.0
        If (L_mcr_qrain)  qrain_star(1:row_length, 1:rows, :) = 0.0
        If (L_mcr_qgraup) qgraup_star(1:row_length, 1:rows, :) = 0.0
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_u(i,j,k) = 0.0
            End Do
          End Do
        End Do
        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = 0.0
            End Do
          End Do
        End Do

      End If ! L_Physics

      If (L_tracer .and. ErrorStatus == 0) then

! store physics changes for use in Sl_tracer2

! DEPENDS ON: tr_set_phys
        call TR_Set_Phys(                                               &
                         super_array_size, super_tracer_phys1,          &
                         L_CO2_interactive, CO2,                        &
                         L_Murk_advect, murk,                           &
                         L_Soot, soot_new,                              &
                                 soot_agd,                              &
                                 soot_cld,                              &
                         L_SULPC_SO2, SO2,                              &
                                      SO4_aitken,                       &
                                      so4_accu,                         &
                                      so4_diss,                         &
                         L_sulpc_nh3, nh3,                              &
                         L_sulpc_dms, dms,                              &
                         L_dust, DUST_DIV1,                             &
                                 DUST_DIV2,                             &
                                 DUST_DIV3,                             &
                                 DUST_DIV4,                             &
                                 DUST_DIV5,                             &
                                 DUST_DIV6,                             &
                         L_biomass, bmass_new,                          &
                                    bmass_agd,                          &
                                    bmass_cld,                          &
                         L_ocff, ocff_new,                              &
                                 ocff_agd,                              &
                                 ocff_cld,                              &
                         L_USE_CARIOLLE, OZONE_TRACER,                  &
                         tracer_phys1, tracer, tracer_ukca,             &
                         row_length, rows,                              &
                         model_levels, tr_levels, tr_vars, tr_ukca,     &
                         offx, offy, model_domain,                      &
                         .true., halo_i, halo_j                         &
                                               )

      end if  ! L_tracer and ErrorStatus == 0
     
      IF (L_mix_ratio) Then
        allocate ( mix_v (1-halo_i:row_length+halo_i,                   &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        allocate ( mix_cl(1-halo_i:row_length+halo_i,                   &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        allocate ( mix_cf(1-halo_i:row_length+halo_i,                   &
     &                   1-halo_j:rows+halo_j, wet_levels) )
        allocate ( mix_v_star (1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy,wet_levels) )
        allocate ( mix_cl_star (1-offx:row_length+offx,                 &
     &                         1-offy:rows+offy,wet_levels) )
        allocate ( mix_cf_star (1-offx:row_length+offx,                 &
     &                         1-offy:rows+offy,wet_levels) )
        IF(L_mcr_qcf2)then
          allocate ( mix_cf2      (1-halo_i:row_length+halo_i,          &
     &                             1-halo_j:rows+halo_j, wet_levels) )
          allocate ( mix_cf2_star (1-offx:row_length+offx,              &
     &                             1-offy:rows+offy,wet_levels) )
        else
          allocate ( mix_cf2      (1,1,1) )
          allocate ( mix_cf2_star (1,1,1) )
        endif
        IF(L_mcr_qrain)then
          allocate ( mix_rain     (1-halo_i:row_length+halo_i,          &
     &                             1-halo_j:rows+halo_j, wet_levels) )
          allocate ( mix_rain_star(1-offx:row_length+offx,              &
     &                             1-offy:rows+offy,wet_levels) )
        else
          allocate ( mix_rain      (1,1,1) )
          allocate ( mix_rain_star (1,1,1) )
        endif
        IF(L_mcr_qgraup)then
          allocate ( mix_graup     (1-halo_i:row_length+halo_i,         &
     &                             1-halo_j:rows+halo_j, wet_levels) )
          allocate ( mix_graup_star(1-offx:row_length+offx,             &
     &                             1-offy:rows+offy,wet_levels) )
        else
          allocate ( mix_graup      (1,1,1) )
          allocate ( mix_graup_star (1,1,1) )
        endif

!  convert q, qcl,qcf to mix_v, mix_cl,mix_cf
! DEPENDS ON: q_to_mix_halo
        call q_to_mix_halo (row_length, rows, wet_levels,               &
     &               halo_i, halo_j,                                    &
     &               Q, QCL, QCF,                                       &
     &               QCF2, QRAIN, QGRAUP,                               &
     &               L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,             &
     &               mix_v, mix_cl, mix_cf,                             &
     &               mix_cf2, mix_rain, mix_graup )

!  convert q, qcl,qcf _star to mix_v, mix_cl,mix_cf _star
!  q_star holds delta_q

! DEPENDS ON: calc_mix_star
        call calc_mix_star                                              &
     &                  (row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, offx, offy,                    &
     &                   Q, QCL, QCF,                                   &
     &                   QCF2, QRAIN, QGRAUP,                           &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   qcf2_star, qrain_star, qgraup_star,            &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   mix_v_star, mix_cl_star, mix_cf_star           &
     &                  ,mix_cf2_star, mix_rain_star, mix_graup_star    &
     &                   )

      Else

        allocate ( mix_v_star (1,1,1) )
        allocate ( mix_cl_star (1,1,1) )
        allocate ( mix_cf_star (1,1,1) )
        allocate ( mix_cf2_star (1,1,1) )
        allocate ( mix_rain_star (1,1,1) )
        allocate ( mix_graup_star (1,1,1) )

      End If    ! L_mix_ratio


! store physics changes for use in Sl_moist_conserve
      If (L_moist_nonhydro_conserve .and.                               &
     &      ErrorStatus == 0 )then
! At this point q_star holds the increment to q
        IF(l_mix_ratio)then
          allocate ( mix_v_phys1(1-offx:row_length+offx,                &
     &                     1-offy:rows+offy, wet_levels) )
          allocate ( mix_cl_phys1(1-offx:row_length+offx,               &
     &                       1-offy:rows+offy, wet_levels) )
          allocate ( mix_cf_phys1(1-offx:row_length+offx,               &
     &                       1-offy:rows+offy, wet_levels) )
          if(L_mcr_qcf2)then
            allocate ( mix_cf2_phys1(1-offx:row_length+offx,            &
     &                               1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_cf2_phys1(1,1,1) )
          endif
          if(L_mcr_qrain)then
            allocate ( mix_rain_phys1(1-offx:row_length+offx,           &
     &                                1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_rain_phys1(1,1,1) )
          endif
          if(L_mcr_qgraup)then
            allocate ( mix_graup_phys1(1-offx:row_length+offx,          &
     &                                 1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_graup_phys1(1,1,1) )
          endif
! DEPENDS ON: calc_mix_star
          call calc_mix_star                                            &
     &                  (row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j, offx, offy,                    &
     &                   Q, QCL, QCF,                                   &
     &                   QCF2, QRAIN, QGRAUP,                           &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   qcf2_star, qrain_star, qgraup_star,            &
!    &                   .false.,.false.,.false.,
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   mix_v_phys1, mix_cl_phys1, mix_cf_phys1        &
     &                   ,mix_cf2_phys1, mix_rain_phys1, mix_graup_phys1&
     &                   )
          if(l_pc2)then
          allocate ( cf_phys1 (1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
          allocate ( cfl_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
          allocate ( cff_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                cf_phys1(i,j,k)  = cf_star(i,j,k)
                cfl_phys1(i,j,k) = cfl_star(i,j,k)
                cff_phys1(i,j,k) = cff_star(i,j,k)
              End Do
            End Do
          End Do
          else
            allocate ( cf_phys1 (1,1,1) )
            allocate ( cfl_phys1 (1,1,1) )
            allocate ( cff_phys1 (1,1,1) )
         endif
        else
          allocate ( q_phys1(1-offx:row_length+offx,                    &
     &                       1-offy:rows+offy, wet_levels) )
          allocate ( qcl_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
          allocate ( qcf_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
         if(l_pc2)then
          allocate ( cf_phys1(1-offx:row_length+offx,                   &
     &                        1-offy:rows+offy, wet_levels) )
          allocate ( cfl_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
          allocate ( cff_phys1(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
         else
           allocate ( cf_phys1 (1,1,1) )
           allocate ( cfl_phys1 (1,1,1) )
           allocate ( cff_phys1 (1,1,1) )
         endif
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_phys1(i,j,k) = q_star(i,j,k)
                qcl_phys1(i,j,k) = qcl_star(i,j,k)
                qcf_phys1(i,j,k) = qcf_star(i,j,k)
              End Do
            End Do
          End Do
          if(L_pc2)then
           Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  cf_phys1(i,j,k)  = cf_star(i,j,k)
                  cfl_phys1(i,j,k) = cfl_star(i,j,k)
                  cff_phys1(i,j,k) = cff_star(i,j,k)
                End Do
              End Do
            End Do
          endif
         if(L_mcr_qcf2)then
           allocate ( qcf2_phys1(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
           Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  qcf2_phys1(i,j,k)  = qcf2_star(i,j,k)
                End Do
              End Do
            End Do
         else
           allocate ( qcf2_phys1(1,1,1) )
         endif
         if(L_mcr_qrain)then
           allocate ( qrain_phys1(1-offx:row_length+offx,               &
     &                            1-offy:rows+offy, wet_levels) )
           Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  qrain_phys1(i,j,k)  = qrain_star(i,j,k)
                End Do
              End Do
            End Do
         else
           allocate ( qrain_phys1(1,1,1) )
         endif
         if(L_mcr_qgraup)then
           allocate ( qgraup_phys1(1-offx:row_length+offx,              &
     &                             1-offy:rows+offy, wet_levels) )
           Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  qgraup_phys1(i,j,k)  = qgraup_star(i,j,k)
                End Do
              End Do
            End Do
         else
           allocate ( qgraup_phys1(1,1,1) )
         endif
        endif       !L_mix_ratio
      Else If ( NumCycles > 1 .and. L_pc2                               &
     &                           .and. ErrorStatus == 0 ) Then
!
! When cycling and PC2 is used cf_star etc need to be reset (at the
! beginning of each new cycle) to the value they had when they
! exited Physics1(). The following arrays hold these values.
!
        allocate ( cf_phys1(1-offx:row_length+offx,                     &
     &                      1-offy:rows+offy, wet_levels) )
        allocate ( cfl_phys1(1-offx:row_length+offx,                    &
     &                         1-offy:rows+offy, wet_levels) )
        allocate ( cff_phys1(1-offx:row_length+offx,                    &
     &                         1-offy:rows+offy, wet_levels) )
        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              cf_phys1(i,j,k)  = cf_star(i,j,k)
              cfl_phys1(i,j,k) = cfl_star(i,j,k)
              cff_phys1(i,j,k) = cff_star(i,j,k)
            End Do
          End Do
        End Do
      end if  !L_moist_nonhydro_conserve

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS1 Atmos_Phys1',6)
! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after Atmos_Phys1 ***'
        ENDIF ! L_print_pe .or. mype ==0
        R_w =0.0
#include "princnorm.h"
      endif !  L_print_L2norms
! ---------------------------------------------------------------

!
! When cycling the following variables need to be reset (at the
! beginning of each new cycle) to the value they had when they
! exited Physics1().
!

      If ( NumCycles > 1 ) Then

        Allocate( R_u_phys1(row_length,rows,model_levels) )
        Allocate( R_v_phys1(row_length,n_rows,model_levels) )
        Allocate( thetastar_phys1(row_length,rows, model_levels) )
        Allocate( qstar_phys1(row_length,rows,wet_levels) )
        Allocate( qclstar_phys1(row_length,rows,wet_levels) )
        Allocate( qcfstar_phys1(row_length,rows,wet_levels) )
        If (L_mcr_qcf2)                                                 &
     &    Allocate( qcf2_star_phys1(row_length,rows,wet_levels) )
        If (L_mcr_qrain)                                                &
     &    Allocate( qrain_star_phys1(row_length,rows,wet_levels) )
        If (L_mcr_qgraup)                                               &
     &    Allocate( qgraup_star_phys1(row_length,rows,wet_levels) )
        If ( L_mix_ratio ) Then
          Allocate ( mix_v_star_phys1(row_length,rows,wet_levels) )
          Allocate ( mix_cl_star_phys1(row_length,rows,wet_levels) )
          Allocate ( mix_cf_star_phys1(row_length,rows,wet_levels) )
          If ( L_mcr_qcf2 )                                             &
     &      Allocate( mix_cf2_star_phys1(row_length,rows,wet_levels) )
          If ( L_mcr_qrain )                                            &
     &      Allocate( mix_rain_star_phys1(row_length,rows,wet_levels) )
          If ( L_mcr_qgraup )                                           &
     &      Allocate( mix_graup_star_phys1(row_length,rows,wet_levels) )
        End If
        Allocate( area_cld_frac_phys1(row_length,rows,wet_levels) )
        Allocate( bulk_cld_frac_phys1(row_length,rows,wet_levels) )
        Allocate( bulk_cld_liq_phys1(row_length,rows, wet_levels) )
        Allocate( bulk_cld_fr_phys1(row_length,rows, wet_levels) )
        Allocate( ti_phys1(row_length, rows, nice) )
        Allocate( zh_phys1(row_length, rows) )
        Allocate( z0msea_phys1(row_length, rows) )
        Allocate( cca_phys1 (row_length, rows, n_cca_lev) )
        Allocate( ccb_phys1 (row_length, rows ) )
        Allocate( cct_phys1 (row_length, rows ) )
        If ( L_ctile ) Then
          Allocate( T_LAND_CTILE_PHYS1(ROW_LENGTH,ROWS) )
          Allocate( T_SICE_CTILE_PHYS1(ROW_LENGTH,ROWS) )
        End If
        Allocate( t_surf_phys1(row_length, rows) )
        Allocate( t_sf_tile_phys1(land_field,ntiles) )
        Allocate( snow_tile_phys1(land_field,ntiles) )
        Allocate( dolr_phys1(row_length,rows) )

        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_u_phys1(i,j,k)       = R_u(i,j,k)
              thetastar_phys1(i,j,k) = theta_star(i,j,k)
              qstar_phys1(i,j,k)     = q_star(i,j,k)
              qclstar_phys1(i,j,k)   = qcl_star(i,j,k)
              qcfstar_phys1(i,j,k)   = qcf_star(i,j,k)
            End Do
          End Do
        End Do

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v_phys1(i,j,k) = R_v(i,j,k)
            End Do
          End Do
        End Do

        Do k = wet_levels, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_u_phys1(i,j,k)       = R_u(i,j,k)
              thetastar_phys1(i,j,k) = theta_star(i,j,k)
            End Do
          End Do
        Enddo

        If ( L_mcr_qcf2 ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qcf2_star_phys1(i,j,k) = qcf2_star(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mcr_qrain ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qrain_star_phys1(i,j,k) = qrain_star(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mcr_qgraup ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qgraup_star_phys1(i,j,k) = qgraup_star(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mix_ratio ) Then

          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                mix_v_star_phys1(i,j,k)  = mix_v_star(i,j,k)
                mix_cl_star_phys1(i,j,k) = mix_cl_star(i,j,k)
                mix_cf_star_phys1(i,j,k) = mix_cf_star(i,j,k)
              End do
            End do
          End do

          If ( L_mcr_qcf2 ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_cf2_star_phys1(i,j,k) = mix_cf2_star(i,j,k)
                End Do
              End Do
            End do
          End If

          If ( L_mcr_qrain ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_rain_star_phys1(i,j,k) = mix_rain_star(i,j,k)
                End Do
              End Do
            End do
          End If

          If ( L_mcr_qgraup ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_graup_star_phys1(i,j,k) = mix_graup_star(i,j,k)
                End Do
              End Do
            End do
          End If

        End If

        Do k = 1, wet_levels
          Do j = j_start, j_stop
!CDIR NODEP
            Do i = i_start, i_stop
              bulk_cld_frac_phys1(i,j,k) = cf_bulk(i,j,k)
              bulk_cld_liq_phys1(i,j,k)  = cf_liquid(i,j,k)
              bulk_cld_fr_phys1(i,j,k)   = cf_frozen(i,j,k)
              area_cld_frac_phys1(i,j,k) = cf_area(i,j,k)
            Enddo
          Enddo
        End do

        Do k=1, n_cca_lev
          Do j=j_start, j_stop
!CDIR NODEP
             Do i=i_start, i_stop
               cca_phys1(i,j,k) = cca(i,j,k)
             Enddo
          Enddo
        Enddo

        Do j=j_start, j_stop
!CDIR NODEP
          Do i=i_start, i_stop
            ccb_phys1(i,j)=ccb(i,j)
            cct_phys1(i,j)=cct(i,j)
          Enddo
        Enddo


! The following piece of code is intended to keep surface
! vars at the same timelevel at each cycle.
        If (nice == 1) Then
! Set sea ice catagory scheme D1 pointers as catogories are not pres
! if nice = 1.
          p_ti        => ti
        else
          p_ti        => ti_cat
        End If

        Do k=1, nice
          Do j = 1, rows
!CDIR NODEP
            Do i = 1, row_length
               ti_phys1(i,j,k) = p_ti(i,j,k)
            Enddo
          Enddo
        Enddo

        Do j = 1, rows
!CDIR NODEP
          Do i = 1, row_length
             z0msea_phys1(i,j) = z0(i,j)
             zh_phys1(i,j) = zh(i,j)
          Enddo
        Enddo

        Do j = 1, rows
!CDIR NODEP
          Do i = 1, row_length
            t_surf_phys1(i,j) = TSTAR(i,j)
            dolr_phys1(i,j)   = dolr(i,j)
          Enddo
        Enddo

        Do i=1, land_field
           GS1(i) = GS(i)
        Enddo

        If ( L_ctile ) Then
          Do j = 1, rows
!CDIR NODEP
            Do i = 1, row_length
              t_land_ctile_phys1(i,j) = TSTAR_LAND(i,j)
              t_sice_ctile_phys1(i,j) = TSTAR_SICE(i,j)
            Enddo
          Enddo
        End If

        Do j = 1, ntiles
!CDIR NODEP
          Do i = 1, land_field
            t_sf_tile_phys1(i,j) = TSTAR_TILE(i+(j-1)*land_field)
            snow_tile_phys1(i,j) = SNODEP_TILE(i+(j-1)*land_field)
          Enddo
        Enddo

      End If  ! Num_Cycles > 1

      If ( L_new_tdisc ) Then

        Allocate( theta_np1(1-offx:row_length+offx,                     &
     &               1-offy:rows+offy, model_levels))
        Allocate( rho_np1 (1-offx:row_length+offx,                      &
     &       1-offy:rows+offy, model_levels) )
        Allocate( u_np1(1-offx:row_length+offx,                         &
     &       1-offy:rows+offy, model_levels ) )
        Allocate( v_np1(1-offx:row_length+offx,                         &
     &       1-offy:n_rows+offy, model_levels) )
        Allocate( w_np1(1-offx:row_length+offx,                         &
     &                  1-offy:rows+offy, 0:model_levels) )

        If ( L_mix_ratio ) Then

          Allocate ( mix_v_np1(1-offx:row_length+offx,                  &
     &                          1-offy:rows+offy, wet_levels) )
          Allocate ( mix_cl_np1(1-offx:row_length+offx,                 &
     &                          1-offy:rows+offy, wet_levels) )
          Allocate ( mix_cf_np1(1-offx:row_length+offx,                 &
     &                          1-offy:rows+offy, wet_levels) )

          If ( L_mcr_qcf2 ) Then
            Allocate( mix_cf2_np1(1-offx:row_length+offx,               &
     &                          1-offy:rows+offy, wet_levels) )
          Else
            Allocate( mix_cf2_np1(1,1,1) )
          End If

          If ( L_mcr_qrain ) Then
            Allocate( mix_rain_np1(1-offx:row_length+offx,              &
     &                          1-offy:rows+offy, wet_levels) )
          Else
            Allocate( mix_rain_np1(1,1,1) )
          End If

          If ( L_mcr_qgraup ) Then
            Allocate( mix_graup_np1(1-offx:row_length+offx,             &
     &                          1-offy:rows+offy, wet_levels) )
          Else
            Allocate( mix_graup_np1(1,1,1) )
          End If

          Allocate( q_np1(1,1,1) )
          Allocate( qcl_np1(1,1,1) )
          Allocate( qcf_np1(1,1,1) )
          Allocate ( qcf2_np1(1,1,1) )
          Allocate ( qrain_np1(1,1,1) )
          Allocate ( qgraup_np1(1,1,1) )

        Else

          Allocate( q_np1(1-offx:row_length+offx,                   &
     &                    1-offy:rows+offy, wet_levels) )
          Allocate( qcl_np1(1-offx:row_length+offx,                 &
     &                      1-offy:rows+offy, wet_levels) )
          Allocate( qcf_np1(1-offx:row_length+offx,                 &
     &                      1-offy:rows+offy, wet_levels) )

          If ( L_mcr_qcf2 ) Then
            Allocate ( qcf2_np1(1-offx:row_length+offx,                 &
     &                          1-offy:rows+offy, wet_levels) )
          Else
            Allocate ( qcf2_np1(1,1,1) )
          End If

          If ( L_mcr_qrain ) Then
            Allocate ( qrain_np1(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
          Else
            Allocate ( qrain_np1(1,1,1) )
          End If

          If ( L_mcr_qgraup ) Then
            Allocate ( qgraup_np1(1-offx:row_length+offx,               &
     &                            1-offy:rows+offy, wet_levels) )
          Else
            Allocate ( qgraup_np1(1,1,1) )
          End If

          Allocate ( mix_v_np1(1,1,1) )
          Allocate ( mix_cl_np1(1,1,1) )
          Allocate ( mix_cf_np1(1,1,1) )
          Allocate( mix_cf2_np1(1,1,1) )
          Allocate( mix_rain_np1(1,1,1) )
          Allocate( mix_graup_np1(1,1,1) )

        End If

      Else

        Allocate( theta_np1(1,1,1) )
        Allocate( u_np1(1,1,1) )
        Allocate( v_np1(1,1,1) )
        Allocate( w_np1(1,1,1) )
        Allocate( rho_np1 (1,1,1) )
        Allocate( q_np1(1,1,1) )
        Allocate( qcl_np1(1,1,1) )
        Allocate( qcf_np1(1,1,1) )
        Allocate( qcf2_np1(1,1,1) )
        Allocate( qrain_np1(1,1,1) )
        Allocate( qgraup_np1(1,1,1) )
        Allocate( mix_v_np1(1,1,1) )
        Allocate( mix_cl_np1(1,1,1) )
        Allocate( mix_cf_np1(1,1,1) )
        Allocate( mix_cf2_np1(1,1,1) )
        Allocate( mix_rain_np1(1,1,1) )
        Allocate( mix_graup_np1(1,1,1) )

      End If  ! L_new_tdisc

! ----------------------------------------------------------------------
! Section 2. Calculate new time-level estimates to Theta, q, qcl, qcf
!            by calling semi-Lagrangian advection routine.
!            Estimated values are returned in the _star variables.
! ----------------------------------------------------------------------
        ALLOCATE (exner_prime( 1-offx:row_length+offx, 1-offy:rows+offy,&
     &                         model_levels) )

        exner_prime(:,:,:) = 0.0

!
! Iterative SISL: Cycle through dynamics-physics to enable trajectory
! calculations from interpolated winds and utilize improved time
! discretization if requested.
!
      Do CycleNo = 1, NumCycles

! Restore phys1 variables to be used as predictors.
      If ( CycleNo > 1 ) Then

! reset weights after the first cycle
        alpha1 = alpha_1_2
        alpha2 = alpha_2_2
        alpha3 = alpha_3_2
        alpha4 = alpha_4_2

        Do k = 1, wet_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_u(i,j,k)        = R_u_phys1(i,j,k)
              theta_star(i,j,k) = thetastar_phys1(i,j,k)
              q_star(i,j,k)     = qstar_phys1(i,j,k)
              qcl_star(i,j,k)   = qclstar_phys1(i,j,k)
              qcf_star(i,j,k)   = qcfstar_phys1(i,j,k)
            End Do
          End Do
        Enddo

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = R_v_phys1(i,j,k)
            End Do
          End Do
        End Do

        Do k = wet_levels, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_u(i,j,k)        = R_u_phys1(i,j,k)
              theta_star(i,j,k) = thetastar_phys1(i,j,k)
            End Do
          End Do
        End Do

        If ( L_mcr_qcf2 ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qcf2_star(i,j,k) = qcf2_star_phys1(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mcr_qrain ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qrain_star(i,j,k) = qrain_star_phys1(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mcr_qgraup ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                qgraup_star(i,j,k) = qgraup_star_phys1(i,j,k)
              End Do
            End Do
          End do
        End If

        If ( L_mix_ratio ) Then

          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                mix_v_star(i,j,k)  = mix_v_star_phys1(i,j,k)
                mix_cl_star(i,j,k) = mix_cl_star_phys1(i,j,k)
                mix_cf_star(i,j,k) = mix_cf_star_phys1(i,j,k)
              End do
            End do
          End do

          If ( L_mcr_qcf2 ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_cf2_star(i,j,k) = mix_cf2_star_phys1(i,j,k)
                End Do
              End Do
            End do
          End If

          If ( L_mcr_qrain ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_rain_star(i,j,k) = mix_rain_star_phys1(i,j,k)
                End Do
              End Do
            End do
          End If

          If ( L_mcr_qgraup ) Then
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_graup_star(i,j,k) = mix_graup_star_phys1(i,j,k)
                End Do
              End Do
            End do
          End If

        End If ! L_mix_ratio

        Do k = 1, wet_levels
          Do j = j_start, j_stop
!CDIR NODEP
            Do i = i_start, i_stop
              cf_bulk(i,j,k)   = bulk_cld_frac_phys1(i,j,k)
              cf_liquid(i,j,k) = bulk_cld_liq_phys1(i,j,k)
              cf_frozen(i,j,k) = bulk_cld_fr_phys1(i,j,k)
              cf_area(i,j,k)   = area_cld_frac_phys1(i,j,k)
            Enddo
          Enddo
        End do

        If ( L_pc2 ) Then
          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                cf_star(i,j,k)  = cf_phys1(i,j,k)
                cfl_star(i,j,k) = cfl_phys1(i,j,k)
                cff_star(i,j,k) = cff_phys1(i,j,k)
              End Do
            End Do
          End Do
        End If


        Do k=1, n_cca_lev
          Do j=j_start, j_stop
!CDIR NODEP
             Do i=i_start, i_stop
               cca(i,j,k) = cca_phys1(i,j,k)
             Enddo
          Enddo
        Enddo

        Do j=j_start, j_stop
!CDIR NODEP
          Do i=i_start, i_stop
            ccb(i,j) = ccb_phys1(i,j)
            cct(i,j) = cct_phys1(i,j)
          Enddo
        Enddo

        Do k=1, nice
          Do j = 1, rows
!CDIR NODEP
            Do i = 1, row_length
               p_ti(i,j,k) = ti_phys1(i,j,k)
            Enddo
          Enddo
        Enddo

        Do j = 1, rows
!CDIR NODEP
          Do i = 1, row_length
            z0(i,j)  = z0msea_phys1(i,j)
            zh(i,j)  = zh_phys1(i,j)
          Enddo
        Enddo

        Do i=1, land_field
           GS(i) = GS1(i)
        Enddo

        Do j = 1, rows
!CDIR NODEP
          Do i = 1, row_length
            TSTAR(i,j) = t_surf_phys1(i,j)
            dolr(i,j) = dolr_phys1(i,j)
          Enddo
        Enddo

        If ( L_ctile ) Then
          Do j = 1, rows
!CDIR NODEP
            Do i = 1, row_length
              TSTAR_LAND(i,j) = t_land_ctile_phys1(i,j)
              TSTAR_SICE(i,j) = t_sice_ctile_phys1(i,j)
            Enddo
          Enddo
        End If

        Do j = 1, ntiles
!CDIR NODEP
          Do i = 1, land_field
            TSTAR_TILE(i+(j-1)*land_field)  = t_sf_tile_phys1(i,j)
            SNODEP_TILE(i+(j-1)*land_field) = snow_tile_phys1(i,j)
          Enddo
        Enddo

      Else

! umui user defined alphas at first cycle
        alpha1 = alpha_1
        alpha2 = alpha_2
        alpha3 = alpha_3
        alpha4 = alpha_4

      End If

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS2 S-L Advection',5)

! Note: The departure points which are passed into sl_full_wind are
! returned in depart_lambda for lambda, depart_phi for phi and
! depart_r_w for r. depart_r_theta is passed into sl_tracer1.

      If (ErrorStatus  ==  0 ) Then

! Apply diagnostics at last cycle only.
       If ( CycleNo == NumCycles ) Then

! Save current values to form diagnostics increments over advection
        IF( sf(185,12) ) THEN
          ALLOCATE ( u_incr_diagnostic(row_length,rows,model_levels) )
! Hold u increment (no halos needed)
          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                u_incr_diagnostic(i,j,k) =R_u(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                   ! on STASHflag

        IF( sf(186,12) ) THEN
          ALLOCATE ( v_incr_diagnostic(row_length,n_rows,model_levels) )
! Hold v increment (no halos needed)
          DO k=1,model_levels
            DO j=1,n_rows
              DO i=1,row_length
                v_incr_diagnostic(i,j,k) =R_v(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(181,12) ) THEN
          ALLOCATE ( T_incr_diagnostic(row_length,rows,model_levels) )
! note at this point theta_star holds all increments to theta
          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                T_incr_diagnostic(i,j,k) =theta_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(182,12) ) THEN
          ALLOCATE ( q_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point q_star holds all increments to q
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                q_incr_diagnostic(i,j,k) =q_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(183,12) ) THEN
          ALLOCATE ( qcl_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point qcl_star holds all increments to qcl
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                qcl_incr_diagnostic(i,j,k) =qcl_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(184,12) ) THEN
          ALLOCATE ( qcf_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point qcf_star holds all increments to qcf
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                qcf_incr_diagnostic(i,j,k) =qcf_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(192,12) ) THEN
          ALLOCATE ( cf_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point cf_star holds all increments to cf
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                cf_incr_diagnostic(i,j,k) =cf_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(193,12) ) THEN
          ALLOCATE ( cfl_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point cfl_star holds all increments to cfl
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                cfl_incr_diagnostic(i,j,k) =cfl_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(194,12) ) THEN
          ALLOCATE ( cff_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point cff_star holds all increments to cff
          DO k=1,wet_levels
            DO j=1,rows
              DO i=1,row_length
                cff_incr_diagnostic(i,j,k) =cff_star(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

       End If ! CycleNo == NumCycles

!      IF( sf(187,12) ) THEN    ! w increments
! Note no need to hold R_w before semi-lagrangian advection as unset.

        If ( CycleNo == 1 ) Then
        ALLOCATE ( depart_lambda(row_length,rows,model_levels) )
        ALLOCATE ( depart_phi(row_length,rows,model_levels) )
        ALLOCATE ( depart_r_w(row_length,rows,model_levels) )
        ALLOCATE ( depart_r_theta(row_length,rows,model_levels) )
        End If

! DEPENDS ON: timer
        If (Ltimer) Call timer ('SL_Thermo',3)

! DEPENDS ON: ni_sl_thermo
        Call NI_SL_Thermo(                                              &
     &                    moisture_array_size,                          &
     &                    theta, q, qcl, qcf, qcf2, qrain, qgraup,      &
     &                    mix_v, mix_cl, mix_cf,                        &
     &                    mix_cf2, mix_rain, mix_graup,                 &
     &                    cf_bulk, cf_liquid, cf_frozen,                &
     &                    q_star, qcl_star, qcf_star,                   &
     &                    qcf2_star, qrain_star, qgraup_star,           &
     &                    mix_v_star, mix_cf_star, mix_cl_star,         &
     &                    mix_cf2_star, mix_rain_star, mix_graup_star,  &
     &                    cf_star, cfl_star, cff_star,                  &
     &                    exner_star, theta_star, theta_np1,            &
     &                    w, w_adv, u_adv, v_adv,                       &
     &                    exner_theta_levels,                           &
     &                    pstar, p, p_theta_levels, rho,                &
     &             eta_rho_levels, eta_theta_levels, r_rho_levels,      &
     &             r_theta_levels, row_length, rows, n_rows,            &
     &                    model_levels, wet_levels,                     &
     &                    alpha_2, check_bottom_levels,                 &
     &                    interp_vertical_search_tol,                   &
     &                    first_constant_r_rho_level,                   &
     &                    delta_lambda, delta_phi,                      &
     &             glambda_p, phi_p, glambda_u, phi_v,                  &
     &             gdlambda_p, dphi_p, gdlambda_u, dphi_v, grecip_dlamp,&
     &             recip_dphip, grecip_dlamu, recip_dphiv,              &
     &             wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v,        &
     &             lambda_p_rm, lambda_p_rp, lambda_u_rm, lambda_u_rp,  &
     &             phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,              &
     &           recip_lambda_p_m, recip_lambda_p_0, recip_lambda_p_p,  &
     &           recip_lambda_p_p2, recip_lambda_u_m, recip_lambda_u_0, &
     &           recip_lambda_u_p, recip_lambda_u_p2, recip_phi_p_m,    &
     &           recip_phi_p_0, recip_phi_p_p, recip_phi_p_p2,          &
     &           recip_phi_v_m, recip_phi_v_0, recip_phi_v_p,           &
     &           recip_phi_v_p2, Base_lambda, base_phi, lambda_p_end,   &
     &           phi_p_end, dlambda_p_end, dphi_p_end, dphi_v_end,      &
     &           recip_dlam, recip_dphi, max_look,                      &
     &           look_lam, look_phi, halo_lam, halo_phi,                &
     &                 timestep, FV_cos_theta_latitude,                 &
     &                 cos_theta_latitude, sec_theta_latitude,          &
     &                 sin_theta_latitude, tan_theta_latitude,          &
     &                 cos_v_latitude, sec_v_latitude,                  &
     &                 sin_v_latitude, tan_v_latitude,                  &
     &                 sin_theta_longitude, cos_theta_longitude,        &
     &                 LAM_max_cfl, THETA_LBC, n_rims_to_do,            &
     &                 RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,           &
     &                 LENRIMA(fld_type_p,halo_type_extended,           &
     &                         rima_type_norm),                         &
     &                 LBC_SIZEA(1,fld_type_p,halo_type_extended,       &
     &                           rima_type_norm),                       &
     &                 LBC_STARTA(1,fld_type_p,halo_type_extended,      &
     &                            rima_type_norm),                      &
     &                 r_at_u, r_at_v,                                  &
     &                 mype, nproc, nproc_x, nproc_y,                   &
     &                 offx, offy, halo_i, halo_j, datastart,           &
     &                 g_i_pe, at_extremity,                            &
     &                 global_row_length, global_rows,                  &
     &                 gc_proc_row_group, gc_proc_col_group,            &
     &                 Pi, Depart_scheme, Depart_order,                 &
     &                 high_order_scheme(Theta_SL),                     &
     &                 monotone_scheme(Theta_SL),                       &
     &                 high_order_scheme(moist_SL),                     &
     &                 monotone_scheme(moist_SL),                       &
     &                 L_Ritchie_high, L_Ritchie_mono,                  &
     &                 Ritchie_high_order_scheme,                       &
     &                 Ritchie_monotone_scheme,                         &
     &                 model_domain, L_high(Theta_SL),                  &
     &                 L_mono(Theta_SL), thmono_levels,                 &
     &                 L_high(moist_SL), L_mono(moist_SL),              &
     &                 L_conserv(moist_SL), L_pc2,                      &
     &                 L_2d_sl_geometry, L_mix_ratio,                   &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 L_regular, L_sl_halo_reprod,                     &
     &                 L_fint_theta, L_lbc_new,                         &
     &                 L_new_tdisc, CycleNo,                            &
     &                 depart_r_theta,                                  &
     &                 depart_lambda, depart_phi, depart_r_w,           &
     &                 ErrorStatus )

        If (.not. L_moist_nonhydro_conserve .and.                       &
     &      .not. L_tracer .and. CycleNo == NumCycles ) Then
          DEALLOCATE (depart_r_theta)
        endif

! DEPENDS ON: timer
        If (Ltimer) Call timer ('SL_Thermo',4)

      End If        !     ErrorStatus  ==  0

! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after sl_thermo  ***'
        ENDIF ! L_print_pe .or. mype ==0
#include "princnorm.h"
      endif !  L_print_L2norms

! ----------------------------------------------------------------------
! Section 2.1 Calculate advective Momentum increments.
! ----------------------------------------------------------------------
      If (ErrorStatus  ==  0 ) Then

! DEPENDS ON: timer
        If (Ltimer) Call timer ('SL_Full_wind',3)

! DEPENDS ON: ni_sl_full_wind
        Call NI_SL_Full_wind(                                           &
     &                      u, u_np1, v, v_np1,                         &
     &                      w, w_np1,                                   &
     &                      u_adv, v_adv, w_adv,                        &
     &                      theta, theta_np1,                           &
     &                      exner_rho_levels,                           &
     &                      q, qcl, qcf,                                &
     &                      qcf2, qrain, qgraup,                        &
     &                      q_np1, qcl_np1, qcf_np1, qcf2_np1,          &
     &                      qrain_np1, qgraup_np1,                      &
     &                      mix_v, mix_cl, mix_cf,                      &
     &                      mix_v_np1, mix_cl_np1, mix_cf_np1,          &
     &                      mix_cf2, mix_rain, mix_graup,               &
     &                      mix_cf2_np1, mix_rain_np1, mix_graup_np1,   &
     &                      L_mix_ratio,                                &
     &                      L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
     &                      depart_lambda, depart_phi, depart_r_w,      &
     &                      r_theta_levels, r_rho_levels,               &
     &                      eta_theta_levels, eta_rho_levels,           &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      sin_theta_latitude, cos_v_latitude,         &
     &                      sec_v_latitude, sin_v_latitude,             &
     &                      tan_theta_latitude, tan_v_latitude,         &
     &                      cos_theta_longitude,                        &
     &                      sin_theta_longitude,                        &
     &                      f1_at_v, f2_at_u, f3_at_u, f3_at_v,         &
     &                      delta_lambda, delta_phi, timestep,          &
     &                      glambda_p, phi_p, glambda_u, phi_v,         &
     &                      gdlambda_p, dphi_p, gdlambda_u, dphi_v,     &
     &            grecip_dlamp, recip_dphip, grecip_dlamu, recip_dphiv, &
     &            wt_lambda_p, wt_phi_p, wt_lambda_u, wt_phi_v,         &
     &            lambda_p_rm, lambda_p_rp, lambda_u_rm, lambda_u_rp,   &
     &            phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,               &
     &            recip_lambda_p_m, recip_lambda_p_0, recip_lambda_p_p, &
     &            recip_lambda_p_p2, recip_lambda_u_m, recip_lambda_u_0,&
     &            recip_lambda_u_p, recip_lambda_u_p2, recip_phi_p_m,   &
     &            recip_phi_p_0, recip_phi_p_p, recip_phi_p_p2,         &
     &            recip_phi_v_m, recip_phi_v_0, recip_phi_v_p,          &
     &            recip_phi_v_p2, base_lambda, base_phi, lambda_p_end,  &
     &            phi_p_end, dlambda_p_end, dphi_p_end, dphi_v_end,     &
     &            recip_dlam, recip_dphi, max_look,                     &
     &            look_lam, look_phi, halo_lam, halo_phi,               &
     &                      alpha_3, alpha_4, LAM_max_cfl,              &
     &                      n_Y_arrays, n_Yw_arrays,                    &
     &                      n_Yd_arrays, n_Ydw_arrays,                  &
     &                      U_LBC,V_LBC,W_LBC,                          &
     &                      LENRIMA(1,1,rima_type_norm),                &
     &                      LBC_SIZEA(1,1,1,rima_type_norm),            &
     &                      LBC_STARTA(1,1,1,rima_type_norm),           &
     &                      RIMWIDTHA(rima_type_norm), n_rims_to_do,    &
     &                      Pi, Cp, epsilon, g,                         &
     &                      model_domain, row_length, rows, n_rows,     &
     &                      model_levels, wet_levels,                   &
     &                      Depart_scheme, Depart_order,                &
     &                      high_order_scheme(Wind_SL),                 &
     &                      monotone_scheme(Wind_SL),                   &
     &                      L_trivial_trigs,                            &
     &                      L_high(Wind_SL), L_mono(Wind_SL),           &
     &                      L_conserv(Wind_SL),                         &
     &                      L_Ritchie_high,                             &
     &                      L_Ritchie_mono,                             &
     &                      Ritchie_high_order_scheme,                  &
     &                      Ritchie_monotone_scheme,                    &
     &                      first_constant_r_rho_level,                 &
     &                      check_bottom_levels,                        &
     &                      interp_vertical_search_tol,                 &
     &                      r_at_u, r_at_v,                             &
     &                      mype, nproc, nproc_x, nproc_y,              &
     &                      offx, offy, halo_i, halo_j,                 &
     &                      global_row_length, global_rows,             &
     &                      datastart, at_extremity, g_i_pe,            &
     &                      gc_proc_row_group, gc_proc_col_group,       &
     &                      L_2d_sl_geometry, L_sl_halo_reprod,         &
     &                      L_free_slip, L_regular,                     &
     &                      L_qwaterload, L_interp_depart,              &
     &                      L_new_tdisc, CycleNo,                       &
     &                      R_u, R_v, R_w, ErrorStatus )

! DEPENDS ON: timer
        If (Ltimer) Call timer ('SL_Full_wind',4)

! ---------------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after SL_Full_wind  ***'
        ENDIF ! L_print_pe .or. mype ==0
#include "princnorm.h"
      endif !  L_print_L2norms
! ---------------------------------------------------------------

        IF ( .Not.SF(0,12) .and. CycleNo == NumCycles ) Then
          DEALLOCATE (depart_r_w)
          If (.not. L_moist_nonhydro_conserve .and.                     &
     &        .not. L_tracer )then
            DEALLOCATE (depart_phi)
            DEALLOCATE (depart_lambda)
          endif
        endif
#if defined(T3E)
        call barrier( )
#endif
      End If         !  ErrorStatus  ==  0

      IF(L_mix_ratio)then
!  convert mix_v, mix_cl,mix_cf to q, qcl,qcf

! DEPENDS ON: mix_to_q
        call mix_to_q                                                   &
     &                  (row_length, rows, wet_levels,                  &
     &                   offx, offy,                                    &
     &                   mix_v_star, mix_cl_star, mix_cf_star,          &
     &                   mix_cf2_star, mix_rain_star, mix_graup_star,   &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   qcf2_star, qrain_star, qgraup_star             &
     &                   )

      endif         !     L_mix_ratio

      If ( L_tracer .and. CycleNo == NumCycles ) Then

! DEPENDS ON: timer
        If (Ltimer) CALL TIMER('SL_tracer1',3)

! DEPENDS ON: sl_tracer1
        Call SL_tracer1(                                                &
     &                 super_array_size, eta_theta_levels,              &
     &                 r_rho_levels, r_theta_levels,                    &
     &                 PSTAR, P, P_THETA_LEVELS,                        &
     &                 RHO, l_tracer1_non_hydro,                        &
     &                 row_length, rows, n_rows, model_levels,          &
     &                 delta_lambda, delta_phi,                         &
     &                 glambda_p, phi_p, gdlambda_u, dphi_v,            &
     &                 grecip_dlamp, recip_dphip,                       &
     &                 lambda_p_rm, lambda_p_rp, phi_p_rm, phi_p_rp,    &
     &                 recip_lambda_p_m, recip_lambda_p_0,              &
     &                 recip_lambda_p_p, recip_lambda_p_p2,             &
     &                 recip_phi_p_m, recip_phi_p_0,                    &
     &                 recip_phi_p_p, recip_phi_p_p2,                   &
     &                 base_lambda, base_phi,                           &
     &                 recip_dlam, recip_dphi, max_look,                &
     &                 look_lam, look_phi, halo_lam, halo_phi,          &
     &                 FV_cos_theta_latitude,                           &
     &                 mype, nproc, nproc_x, nproc_y,                   &
     &                 halo_i, halo_j, datastart,                       &
     &                 g_i_pe, at_extremity,                            &
     &                 global_row_length, global_rows,                  &
     &                 gc_proc_row_group,                               &
     &                 gc_proc_col_group, offx, offy,                   &
     &                 L_regular, L_sl_halo_reprod,                     &
     &                 high_order_scheme(moist_SL),                     &
     &                 monotone_scheme(moist_SL),                       &
     &                 model_domain, L_high(moist_SL),                  &
     &                 L_mono(moist_SL),                                &
     &                 L_conserve_tracers,                              &
     &                 check_bottom_levels,                             &
     &                 interp_vertical_search_tol,                      &
     &                 first_constant_r_rho_level,                      &
     &                 depart_lambda, depart_phi, depart_r_theta,       &

!  d1(jxxx) holds time level n value plus physics1 increment
     &                 CO2, L_CO2_interactive,                          &
     &                 MURK, L_Murk_advect,                             &
     &                 DUST_DIV1,DUST_DIV2,                             &
     &                 DUST_DIV3,DUST_DIV4,                             &
     &                 DUST_DIV5,DUST_DIV6, L_DUST,                     &
     &                 SOOT_NEW, SOOT_AGD,                              &
     &                 SOOT_CLD, L_soot,                                &
     &                 BMASS_NEW, BMASS_AGD,                            &
     &                 BMASS_CLD, L_biomass,                            &
     &                 OCFF_NEW, OCFF_AGD, OCFF_CLD, L_OCFF,            &
     &                 SO2, SO4_AITKEN,                                 &
     &                 SO4_ACCU,                                        &
     &                 SO4_DISS, NH3, DMS,                              &
     &                 L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms,           &
     &                 TRACER, tr_levels, tr_vars,                      &
     &                 tracer_ukca, tr_ukca,                            &
     &                 L_USE_CARIOLLE,    OZONE_TRACER,                 &
     &                 ErrorStatus)

! DEPENDS ON: timer
        If (Ltimer) CALL TIMER('SL_tracer1',4)

! ----------------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after SL_tracer1  ***'
        ENDIF ! L_print_pe .or. mype ==0
#include "princnorm.h"
      endif !  L_print_L2norms
! ---------------------------------------------------------------

      end if  ! L_tracer

! ----------------------------------------------------------------------
! Section 2.21 Calculation of coefficients in turbulence scheme
! ----------------------------------------------------------------------
      If (L_subfilter_horiz .or. L_subfilter_vert) Then

        If (bl_levels  /=  model_levels - 1) then
          write(6,*)'BL_LEVELS =',bl_levels                             &
     &,             'MODEL_LEVELS =',model_levels
          ErrorStatus=123
! DEPENDS ON: ereport
          Call Ereport("ATM_STEP", ErrorStatus,                         &
     &        "The number of boundary layer levels must equal " //      &
     &        "the number of model levels minus one when the " //       &
     &        "3D subgrid turbulence scheme is selected.")
        End If

        If (diff_factor  <   0.0 .or. diff_factor  >   1.0) then
          write(6,*)'DIFF_FACTOR =',diff_factor
          ErrorStatus=123
! DEPENDS ON: ereport
          Call Ereport("ATM_STEP", ErrorStatus,                         &
     &        "DIFF_FACTOR must have a value greater than " //          &
     &        "zero and less than or equal to 1.0 so that the "//       &
     &        "numerical stability is maintained")
        End If

        If (mix_factor  <=  0.0 .or. mix_factor  >   1.0) then
          write(6,*)'MIX_FACTOR =',mix_factor
          ErrorStatus=123
! DEPENDS ON: ereport
          Call Ereport("ATM_STEP", ErrorStatus,                         &
     &        "MIX_FACTOR should have a value greater or equal " //     &
     &        "to zero and less than or equal to 1.0 ")
        End If

        If (turb_startlev_horiz > turb_endlev_horiz .OR.                &
     &      turb_startlev_vert > turb_endlev_vert) Then
          ErrorStatus=123
! DEPENDS ON: ereport
          Call Ereport("ATM_STEP", ErrorStatus,                         &
     &        "The start level for the turbulence scheme is " //        &
     &        "greater than the end level! ")
        End If

! The levels over which the turbulence scheme acts must be
! between 2 and model_levels-1

        If (turb_startlev_vert < 2) turb_startlev_vert = 2
        If (turb_endlev_vert > model_levels - 1)                        &
     &                           turb_endlev_vert = model_levels - 1

        If (turb_startlev_horiz < 2) turb_startlev_horiz = 2
        If (turb_endlev_horiz > model_levels - 1)                       &
     &                           turb_endlev_horiz = model_levels - 1

        If ( CycleNo == 1 ) Then
          ALLOCATE (visc_m(1-halo_i:row_length+halo_i                   &
     &                           , 1-halo_j:rows+halo_j, model_levels) )
          ALLOCATE (visc_h(1-halo_i:row_length+halo_i                   &
     &                           , 1-halo_j:rows+halo_j, model_levels) )
          ALLOCATE (RNEUTML(1:row_length, 1:rows,model_levels))
          ALLOCATE (shear(1-offx:row_length+offx, 1-offy:rows+offy      &
     &                                            , model_levels))
        End If

        Do k = 1, model_levels
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              visc_m(i,j,k)= 0.0
              visc_h(i,j,k)= 0.0
            End Do
          End Do
        End Do

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              shear(i,j,k)= 0.0
              RNEUTML(i,j,k) = 0.0
            End Do
          End Do
        End Do

        If (L_subfilter_vert) Then

          ALLOCATE (visc_BL_m(1:row_length, 1:rows, bl_levels))

          Do k = 1, bl_levels
            Do j = 1, rows
              Do i = 1, row_length
                visc_BL_m(i,j,k) = 0.0
              End Do
            End Do
          End Do

        Else

          ALLOCATE (visc_BL_m(1,1,1) )

        End If    ! L_subfilter_vert

!
! Calculate  lambda^2*S in TURB_Smagorinsky
!
! DEPENDS ON: turb_smagorinsky
        CALL TURB_Smagorinsky(                                          &
     &                     rows, row_length, n_rows                     &
     &,                    model_levels                                 &
     &,                    r_theta_levels, r_rho_levels                 &
     &,                    U, V, W                                      &
     &,                    visc_m, shear                                &
     &,                    Z0, RNEUTML, timestep                        &
     &,                    diff_factor, mix_factor, max_diff            &
     &,                    cos_theta_latitude                           &
     &,                    cos_v_latitude                               &
     &,                    delta_lambda, delta_phi)


        If (L_subfilter_vert) Then

          Do k = 2, bl_levels
            Do j = 1, rows
              Do i = 1, row_length
                visc_BL_m(i,j,k)= visc_m(i,j,k-1)
              End Do
            End Do
          End Do

        End If ! L_subfilter_vert

      Else

        If ( CycleNo == 1 ) Then
          ALLOCATE (visc_m(1,1,1) )
          ALLOCATE (visc_h(1,1,1) )
          ALLOCATE (visc_BL_m(1,1,1) )
          ALLOCATE (RNEUTML(1,1,1) )
          ALLOCATE (shear(1,1,1) )
        End If

      End If  !L_subfilter_horiz or L_subfilter_vert

! ---------------------------------------------------------------
! Section 2.3 Diagnostics at end of advection
! ----------------------------------------------------------------------
! Apply diagnostics at final cycle only
      If ( CycleNo == NumCycles ) Then

! DEPENDS ON: timer
      If (Ltimer) Call timer ('Diag_adv',3)

! section 12: 'dynamics advection' based quantities
      IF(      SF(0,12)                                                 &
                                ! Diagnostics required for this section
     &   .AND. ErrorStatus == 0) THEN

! Allocate diagnostic space for STASH
        ALLOCATE (STASHwork12(STASH_maxlen(12,A_im)))

! DEPENDS ON: diagnostics_adv
        CALL Diagnostics_adv(                                           &
     &            row_length, rows, n_rows, model_levels, wet_levels,   &
! primary wind fields:
     &            U, V,                                                 &
     &            THETA, Q, QCL, QCF,                                   &
     &            CF_BULK, CF_LIQUID, CF_FROZEN,                        &
! wind field increments after advection:
     &            R_u, R_v, R_w,                                        &
! wind field increments before advection (on stashflag):
     &            u_incr_diagnostic, v_incr_diagnostic,                 &
     &            T_incr_diagnostic, q_incr_diagnostic,                 &
     &            qcl_incr_diagnostic, qcf_incr_diagnostic,             &
     &            cf_incr_diagnostic, cfl_incr_diagnostic,              &
     &            cff_incr_diagnostic,                                  &
     &            theta_star, q_star, qcl_star, qcf_star,               &
     &            cf_star, cfl_star, cff_star,                          &
     &            EXNER_THETA_LEVELS,                                   &
! Departure points for w
     &            depart_lambda, depart_phi, depart_r_w,                &
     &            r_theta_levels,                                       &
#include "argsts.h"
     &            STASHwork12)

        DEALLOCATE (depart_r_w)
        If (.not. L_moist_nonhydro_conserve .and.                       &
     &      .not. L_tracer )then
          DEALLOCATE (depart_phi)
          DEALLOCATE (depart_lambda)
        endif

! Tidy allocatable arrays
        IF( sf(185,12) ) THEN
          DEALLOCATE ( u_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(186,12) ) THEN
          DEALLOCATE ( v_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(181,12) ) THEN
          DEALLOCATE ( T_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(182,12) ) THEN
          DEALLOCATE ( q_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(183,12) ) THEN
          DEALLOCATE ( qcl_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(184,12) ) THEN
          DEALLOCATE ( qcf_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(192,12) ) THEN
          DEALLOCATE ( cf_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(193,12) ) THEN
          DEALLOCATE ( cfl_incr_diagnostic )
        ENDIF  ! on STASHflag

        IF( sf(194,12) ) THEN
          DEALLOCATE ( cff_incr_diagnostic )
        ENDIF  ! on STASHflag
        
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',3)

! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,12,STASHwork12,                            &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('STASH',4)

        DEALLOCATE (STASHwork12)

      ENDIF !   SF(0,12)

! DEPENDS ON: timer
      If (Ltimer) Call timer ('Diag_adv',4)

      End If ! CycleNo == NumCycles

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS2 S-L Advection',6)

! ----------------------------------------------------------------------
! Section 3.0  Call Atmospheric Physics2
! ----------------------------------------------------------------------
      If (.NOT. L_run_with_physics2) Then
        L_physics_store=l_physics
        L_physics=.false.
      end if

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS3 Atmos_Phys2',5)

      If (ErrorStatus  ==  0 ) Then

        If (L_idealised_data) Then

          If ( CycleNo == 1 ) Then

            ALLOCATE (q_inc_subs(row_length, rows, model_levels))
            ALLOCATE (th_inc_subs(row_length, rows, model_levels))
            ALLOCATE (q_inc_ls(row_length, rows, model_levels))
            ALLOCATE (th_inc_ls(row_length, rows, model_levels))
            ALLOCATE (u_inc_dmp(row_length, rows, model_levels))
            ALLOCATE (q_inc_dmp(row_length, rows, model_levels))
            ALLOCATE (th_inc_dmp(row_length, rows, model_levels))
            ALLOCATE (v_inc_dmp(row_length, n_rows, model_levels))

          End If
! include idealised forcing for theta.

! DEPENDS ON: idl_force
          CALL IDL_Force(                                               &
     &                   R, g, kappa, pi, p_zero                        &
     &,                  row_length, rows, model_levels, timestep       &
     &,                  delta_lambda, delta_phi, model_domain          &
     &,                  lambda_half_width, phi_half_width              &
     &,                  p_max, p_top, p_bottom, p_arbitrary            &
     &,                  lambda_heat_centre, phi_heat_centre            &
     &,                  max_heat_per_day, newtonian_timescale          &
! Dynamical core settings
     &,                  SuHe_newtonian_timescale_ka                    &
     &,                  SuHe_newtonian_timescale_ks                    &
     &,                  SuHe_pole_equ_deltaT, SuHe_static_stab         &
     &,                  SuHe_level_weight, SuHe_sigma_cutoff           &
     &,                  L_SH_Williamson, SuHE_relax                    &
     &,                  L_damp, L_geo_for                              &
     &,                  L_bomex                                        &
     &,                  DMPTIM, HDMP, ZDMP                             &
     &,                  u_geo, v_geo                                   &
     &,                  offx, offy, datastart, at_extremity            &
     &,                  EXNER_THETA_LEVELS                             &
     &,                  P, theta_star                                  &
     &,                  P_THETA_LEVELS, PSTAR                          &
     &,                  THETA                                          &
     &,                  cos_theta_latitude,sin_theta_latitude          &
     &,                  cool_rate, theta_surface                       &
     &,                  timestep_number                                &
     &,                  max_model_levels, max_num_force_times          &
     &,                  Q, q_star                                      &
     &,                  U, V, R_u, R_v                                 &
     &,                  u_ref, v_ref, theta_ref, n_rows                &
     &,                  q_ref                                          &
     &,                  eta_theta_levels, eta_rho_levels               &
     &,                  height_domain                                  &
     &,                  global_row_length, global_rows                 &
     &,                  gc_all_proc_group, nproc                       &
     &,                  tforce_option, qforce_option, uvforce_option   &
     &,                  num_tforce_times, num_qforce_times             &
     &,                  num_uvforce_times                              &
     &,                  tforce_time_interval, qforce_time_interval     &
     &,                  uvforce_time_interval                          &
     &,                  tforce_data_modlev, qforce_data_modlev         &
     &,                  uforce_data_modlev, vforce_data_modlev         &
     &,                  r_rho_levels, r_theta_levels                   &
     &,                  halo_i, halo_j                                 &
     &,                  q_inc_subs, th_inc_subs                        &
     &,                  q_inc_ls, th_inc_ls                            &
     &,                  u_inc_dmp, q_inc_dmp, th_inc_dmp               &
     &,                  v_inc_dmp                                      &
     &,                  f3_at_u, f3_at_v                               &
     &,                  L_physics, problem_number, L_force, LTimer)


          If ( CycleNo == NumCycles ) Then
            DEALLOCATE (q_inc_subs)
            DEALLOCATE (th_inc_subs)
            DEALLOCATE (q_inc_ls)
            DEALLOCATE (th_inc_ls)
            DEALLOCATE (u_inc_dmp)
            DEALLOCATE (q_inc_dmp)
            DEALLOCATE (th_inc_dmp)
            DEALLOCATE (v_inc_dmp)
          End If

        End If      !  L_idealised_data

        If (L_tracer) then

! protect from multiple mem allocations
! save input fields to obtain Atmos_Physics2 increments

! tracers only at final cycle
          If ( CycleNo == NumCycles ) Then
! DEPENDS ON: tr_set_phys
            call TR_Set_Phys(                                           &
                         super_array_size, super_tracer_phys2,          &
                         L_CO2_interactive, CO2,                        &
                         L_Murk_advect, murk,                           &
                         L_Soot, soot_new,                              &
                                 soot_agd,                              &
                                 soot_cld,                              &
                         L_SULPC_SO2, SO2,                              &
                                      SO4_aitken,                       &
                                      so4_accu,                         &
                                      so4_diss,                         &
                         L_sulpc_nh3, nh3,                              &
                         L_sulpc_dms, dms,                              &
                         L_dust, DUST_DIV1,                             &
                                 DUST_DIV2,                             &
                                 DUST_DIV3,                             &
                                 DUST_DIV4,                             &
                                 DUST_DIV5,                             &
                                 DUST_DIV6,                             &
                         L_biomass, bmass_new,                          &
                                    bmass_agd,                          &
                                    bmass_cld,                          &
                         L_ocff, ocff_new,                              &
                                 ocff_agd,                              &
                                 ocff_cld,                              &
                         L_USE_CARIOLLE, OZONE_TRACER,                  &
                         tracer_phys1, tracer, tracer_ukca,             &
                         row_length, rows,                              &
                         model_levels, tr_levels, tr_vars, tr_ukca,     &
                         offx, offy, model_domain,                      &
                         .false., 0, 0                                 &
                         )
          End If ! CycleNo == NumCycles

        end if  ! L_tracer

! save star fields to obtain increments after call to Atmos_Physics2
        If (L_moist_nonhydro_conserve)then
! after SL_thermo q_star holds q_dash (latest estimate to q_(n+1))
          IF(l_mix_ratio)then
! protect from multiple mem allocations
            If ( CycleNo == 1 ) Then
            allocate ( mix_v_phys2(1-offx:row_length+offx,              &
     &                             1-offy:rows+offy, wet_levels) )
            allocate ( mix_cl_phys2(1-offx:row_length+offx,             &
     &                              1-offy:rows+offy, wet_levels) )
            allocate ( mix_cf_phys2(1-offx:row_length+offx,             &
     &                              1-offy:rows+offy, wet_levels) )
            if(L_mcr_qcf2)then
              allocate ( mix_cf2_phys2(1-offx:row_length+offx,          &
     &                                 1-offy:rows+offy, wet_levels) )
            else
              allocate ( mix_cf2_phys2(1,1,1) )
            endif
            if(L_mcr_qrain)then
              allocate ( mix_rain_phys2(1-offx:row_length+offx,         &
     &                                  1-offy:rows+offy, wet_levels) )
            else
              allocate ( mix_rain_phys2(1,1,1) )
            endif
            if(L_mcr_qgraup)then
              allocate ( mix_graup_phys2(1-offx:row_length+offx,        &
     &                                   1-offy:rows+offy, wet_levels) )
            else
              allocate ( mix_graup_phys2(1,1,1) )
            endif
            End If ! CycleNo == 1
! DEPENDS ON: q_to_mix
            call q_to_mix (row_length, rows, wet_levels,                &
     &               offx,offy     ,                                    &
     &               q_star, qcl_star, qcf_star,                        &
     &               qcf2_star, qrain_star, qgraup_star,                &
!    &               .false. ,.false. ,.false.,
     &               L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,             &
     &               mix_v_phys2, mix_cl_phys2, mix_cf_phys2,           &
     &               mix_cf2_phys2, mix_rain_phys2, mix_graup_phys2     &
     &               )
            if(L_pc2)then
! protect from multiple mem allocations
            If ( CycleNo == 1 ) Then
            allocate ( cf_phys2 (1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            allocate ( cfl_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            allocate ( cff_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            End If ! CycleNo == 1
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  cf_phys2(i,j,k)  = cf_star(i,j,k)
                  cfl_phys2(i,j,k) = cfl_star(i,j,k)
                  cff_phys2(i,j,k) = cff_star(i,j,k)
                End Do
              End Do
            End Do
            else
              If ( CycleNo == 1 ) Then
                allocate ( cf_phys2(1,1,1) )
                allocate ( cfl_phys2(1,1,1) )
                allocate ( cff_phys2(1,1,1) )
              End If
            endif
          else

! protect from multiple mem allocations
            If ( CycleNo == 1 ) Then
            allocate ( q_phys2(1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy, wet_levels) )
            allocate ( qcl_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            allocate ( qcf_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            if(L_pc2)then
            allocate ( cf_phys2(1-offx:row_length+offx,                 &
     &                          1-offy:rows+offy, wet_levels) )
            allocate ( cfl_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            allocate ( cff_phys2(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
            else
              allocate ( cf_phys2(1,1,1) )
              allocate ( cfl_phys2(1,1,1) )
              allocate ( cff_phys2(1,1,1) )
            endif
            End If ! CycleNo == 1
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  q_phys2(i,j,k) = q_star(i,j,k)
                  qcl_phys2(i,j,k) = qcl_star(i,j,k)
                  qcf_phys2(i,j,k) = qcf_star(i,j,k)
                End Do
              End Do
            End Do
            if(L_pc2)then
              Do k = 1, wet_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    cf_phys2(i,j,k)  = cf_star(i,j,k)
                    cfl_phys2(i,j,k) = cfl_star(i,j,k)
                    cff_phys2(i,j,k) = cff_star(i,j,k)
                  End Do
                End Do
              End Do
            endif
           if(L_mcr_qcf2)then
             If ( CycleNo == 1 )                                        &
     &       allocate ( qcf2_phys2(1-offx:row_length+offx,              &
     &                             1-offy:rows+offy, wet_levels) )
             Do k = 1, wet_levels
               Do j = 1, rows
                 Do i = 1, row_length
                   qcf2_phys2(i,j,k)  = qcf2_star(i,j,k)
                 End Do
               End Do
             End Do
           else
             If ( CycleNo == 1 ) allocate ( qcf2_phys2(1,1,1) )
           endif
           if(L_mcr_qrain)then
             If ( CycleNo == 1 )                                        &
     &       allocate ( qrain_phys2(1-offx:row_length+offx,             &
     &                              1-offy:rows+offy, wet_levels) )
             Do k = 1, wet_levels
               Do j = 1, rows
                 Do i = 1, row_length
                   qrain_phys2(i,j,k)  = qrain_star(i,j,k)
                 End Do
               End Do
             End Do
           else
             If ( CycleNo == 1 ) allocate ( qrain_phys2(1,1,1) )
           endif
           if(L_mcr_qgraup)then
             If ( CycleNo == 1 )                                        &
     &       allocate ( qgraup_phys2(1-offx:row_length+offx,            &
     &                               1-offy:rows+offy, wet_levels) )
             Do k = 1, wet_levels
               Do j = 1, rows
                 Do i = 1, row_length
                   qgraup_phys2(i,j,k)  = qgraup_star(i,j,k)
                 End Do
               End Do
             End Do
           else
             If ( CycleNo == 1 ) allocate ( qgraup_phys2(1,1,1) )
           endif
          endif            !L_mix_ratio

        end if  !L_moist_nonhydro_conserve

        If (L_Physics) Then

          If (.NOT. (L_mix_ratio .and. L_moist_nonhydro_conserve)       &
                     .and. CycleNo == 1) Then
            Allocate ( mix_v_phys2 (1-offx:row_length+offx,             &
                         1-offy:rows+offy, wet_levels) )
            Allocate ( mix_cl_phys2(1-offx:row_length+offx,             &
                         1-offy:rows+offy, wet_levels) )
            Allocate ( mix_cf_phys2(1-offx:row_length+offx,             &
                         1-offy:rows+offy, wet_levels) )
            If(L_mcr_qcf2)Then
              Allocate ( mix_cf2_phys2(1-offx:row_length+offx,          &
                                   1-offy:rows+offy, wet_levels) )
            Else
              Allocate ( mix_cf2_phys2(1,1,1) )
            End If
            If(L_mcr_qrain)Then
              Allocate ( mix_rain_phys2(1-offx:row_length+offx,         &
                                   1-offy:rows+offy, wet_levels) )
            Else
              Allocate ( mix_rain_phys2(1,1,1) )
            End If
            If(L_mcr_qgraup)Then
              Allocate ( mix_graup_phys2(1-offx:row_length+offx,        &
                                   1-offy:rows+offy, wet_levels) )
            Else
              Allocate ( mix_graup_phys2(1,1,1) )
            End If

! convert q_star,qcl_star,qcf_star to mix_v_phys2,mix_cl_phys2,mix_cf_phys2
! DEPENDS ON: q_to_mix
            call q_to_mix (row_length, rows, wet_levels,                &
                     offx,offy,                                         &
                     q_star, qcl_star, qcf_star,                        &
                     qcf2_star, qrain_star, qgraup_star,                &
                     L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,             &
                     mix_v_phys2, mix_cl_phys2, mix_cf_phys2,           &
                     mix_cf2_phys2, mix_rain_phys2, mix_graup_phys2     &
                     )
          End If ! (.NOT. (L_mix_ratio .and. L_moist_nonhydro_conserve))

! NB: the star variables and R_u and R_v have not been set in the
!     halo region yet.

! Apply diagnostics only at last cycle.
         If ( CycleNo == NumCycles ) Then

! Allocate diagnostic space for STASH
      IF(L_bl)then
          ALLOCATE (STASHwork3(STASH_maxlen(3,A_im)))
      endif
          IF(convection_option == 3) THEN
            ALLOCATE (STASHwork5(STASH_maxlen(5,A_im)))
          ENDIF   ! on convection_option
      If (L_hydrology)then
          ALLOCATE (STASHwork8(STASH_maxlen(8,A_im)))
      endif
          ALLOCATE (STASHwork9(STASH_maxlen(9,A_im)))
          ALLOCATE (STASHwork19(STASH_maxlen(19,A_im)))
          ALLOCATE (STASHwork26(STASH_maxlen(26,A_im)))

         End If  ! CycleNo == NumCycles

! DEPENDS ON: timer
          If (Ltimer) Call timer ('Atmos_Physics2',3)

! NB if you are changing the argument list to atmos_physics2, please
! do an equivalent change to the two places in routine scm_main to keep
! the single column model consistent.

          IF (nice  ==  1)THEN
! Set sea ice catagory scheme D1 pointers as catogories are not present
! if nice = 1.
            p_ti        => ti
            p_ice_fract => ice_fraction
            p_ice_thick => ice_thickness
          ELSE
            p_ti        => ti_cat
            p_ice_fract => ice_fract_cat
            p_ice_thick => ice_thick_cat
          ENDIF


          If (L_RHCPT) Then
!     Dimension diagnostic 3D RHcrit array
            rhc_row_length = row_length
            rhc_rows = rows
          Else
!     RHcrit will be a 1D parametrized array input from user interface
            rhc_row_length = 1
            rhc_rows = 1
          End if    ! L_RHCPT

         If ( CycleNo == 1 .and. L_physics ) Then
           ALLOCATE ( RHCPT(rhc_row_length, rhc_rows, wet_levels) )
         End If
       If ( CycleNo == 1 ) Then
         If (L_subfilter_horiz .or. L_subfilter_vert) then
           ALLOCATE (FM_3D(row_length,rows,BL_LEVELS))
           ALLOCATE (FH_3D(row_length,rows,BL_LEVELS))
         Else
           ALLOCATE (FM_3D(1,1,1))
           ALLOCATE (FH_3D(1,1,1))
         Endif
       End If

!kdcorbin, 05/10 - initialize mass on first timestep
if (l_co2_mass .or. l_tracer_mass) then
   if (l_tracer_mass .and. i_tracermass_start .lt. 1) then
       i_tracermass_start=1
   endif
   if (l_tracer_mass .and. i_tracermass_start .gt. tr_vars) then
       l_tracer_mass=.false.
       write(6,*) '***********************************'
       write(6,*) 'Tracer Mass Error: '
       write(6,*) '  Incorrect assignment of starting tracer'
       write(6,*) 'Turning off tracer mass fixer'
       write(6,*) 'Specify I_TRACERMASS_START in CNTLATM'
   endif

 if (first_atmstep_call) then
   !DEPENDS ON: TRACER_MASSINIT
   Call TRACER_MASSINIT(     &
  &     row_length, rows, model_levels,   &
  &     tr_levels, tr_vars,               &
  &     halo_i, halo_j, offx, offy, mype, timestep_number,  &
  &     r_theta_levels, r_rho_levels, exner_theta_levels,   &
  &     FV_cos_theta_latitude, delta_lambda, delta_phi,     &
  &     RHO, q, qcl, qcf,                                   &
!CO2 and tracer mass flags
  &     L_CO2_MASS,L_TRACER_MASS,                           &
!CO2 and tracer 3-D fields
  &     CO2, tracer,                                        &
  &     tmass)
   endif
endif

! NB if you are changing the argument list to atmos_physics2, please
! do an equivalent change in routine scm_main to keep the single column
! model consistent.  Note also that there are two calls to
! atmos_physics2 in scm_main.

! DEPENDS ON: atmos_physics2
          Call Atmos_Physics2(halo_i,halo_j,offx,offy,global_row_length &
! Parallel variables
     &, global_rows, gc_proc_row_group, gc_proc_col_group, at_extremity &
     &, nproc, nproc_x, nproc_y, neighbour, g_rows, g_row_length        &
     &,     g_datastart, mype, NumCycles, CycleNo                       &

! field dimensions etc.
     &,     row_length, rows, n_rows, land_field, model_levels, nice    &
     &,     wet_levels, bl_levels, st_levels, sm_levels, cloud_levels   &
     &,     land_ice_points, soil_points, n_cca_lev, ntiles, tr_levels  &
     &,     first_constant_r_rho_level, DIM_CS1, DIM_CS2                &

! IN Substepping information and switches
     &,     Num_Substeps, L_phys2_substep, L_regular                    &
     &,     l_mr_physics2, model_domain, L_dry                          &
     &,     FORMDRAG, OROG_DRAG_PARAM, LCAL360, L_emcorr, Ltimer        &
     &,     L_us_blsol,BL_OPTIONS,L_CLD_AREA,L_ACF_Cusack,L_ACF_Brooks  &
     &,     L_RHCPT,L_hydrology,L_bl,L_Murk,L_murk_advect,L_murk_source &
     &,     L_BL_TRACER_MIX,L_DUST,L_CAM_DUST,L_SULPC_SO2,L_SULPC_NH3   &
     &,     L_SBLeq,L_SBLco,L_sulpc_dms,L_soot,L_biomass,L_ocff         &
     &,     L_co2_interactive,L_ctile, L_co2_emits, L_use_bl_diag_term  &
     &,     L_pc2,L_pc2_reset,L_eacf,L_LAMBDAM2,L_FULL_LAMBDAS,L_ukca   &
     &,     L_sice_heatflux, L_OHADGEM1, L_USE_CARIOLLE                 &
     &,     l_anthrop_heat_src                                          &
! Model Parameters
     &,     alpha_Cd, Puns, Pstb, RHcrit, CO2_MMR, tr_vars, tr_ukca     &
     &,     Muw_SBL, Mwt_SBL, cloud_fraction_method, overlap_ice_liquid &
     &,     ice_fraction_method, ctt_weight, t_weight, qsat_fixed       &
     &,     sub_cld, x1i, x1ic, x1r, x2r, x4r, l_psd, ai, bi, aic, bic  &
     &,     lsp_ei,lsp_fi,lsp_eic,lsp_fic                               &
     &,     dbsdtbs_turb_0, Charnock, SeaSalinityFactor                 &
! Physical constants
     &,Lc,Lf,Cp, two_Omega,p_zero, kappa, R, g, Lapse, earth_radius, Pi &
! Vertical coordinate levels.
     &,r_rho_levels, r_theta_levels, r_at_u, r_at_v, unscaled_dry_rho   &
     &,     eta_theta_levels, eta_rho_levels, delta_lambda, delta_phi   &
     &,     gdlambda_p(lambda_start), dphi_p                            &
     &,     wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v                &
     &,     lat_rot_NP, long_rot_NP, f3_at_u                            &
! Variables required by STPH_SCV
! add mix ratios 
     &,     mix_v_phys2, mix_cl_phys2, mix_cf_phys2                     &
! Time stepping information and trig arrays
     &, timestep, I_year, I_day_number, I_hour, I_minute, I_second      &
     &, timestep_number, sin_theta_longitude, cos_theta_longitude       &
!     &,     FV_cos_theta_latitude,sec_theta_latitude,                   &
     &, sin_theta_latitude,FV_cos_theta_latitude,sec_theta_latitude,    &
! River routing
     &AOCPL_ROW_LENGTH,AOCPL_P_ROWS,L_RIVERS,XPA,XUA,XVA,YPA,YUA,YVA,   &
     &G_P_FIELD,G_R_FIELD,A_INTHD(16),lasize(1,fld_type_r,              &
     &halo_type_no_halo),lasize(2,fld_type_r,halo_type_no_halo),        &
     &glsize(1,fld_type_r),glsize(2,fld_type_r),RIVER_STEP,RIVER_VEL,   &
     &RIVER_MCOEF,RIV_DIRECTION,RIV_SEQUENCE,RIV_STORAGE,               &
!  Add inland basin outflow to arguments
     &  RIV_INLANDATM,L_INLAND,                                         &


! Grid-to-grid river routing
     &  RIV_IAREA,RIV_SLOPE, RIV_FLOWOBS1,RIV_INEXT                     &
     &  ,RIV_JNEXT,RIV_LAND,RIV_SUBSTORE,                               &
     &  RIV_SURFSTORE,RIV_FLOWIN,RIV_BFLOWIN,                           &

!
! diagnostic info
#include "argsts.h"
     &      STASHwork3,STASHwork5,STASHwork8,STASHwork9,STASHwork19,    &
     &      STASHwork26,                                                &
!
! SCM Diagnostics (dummy values in full UM)
     &  nSCMDpkgs, L_SCMDiags,                                          &
!
! Data Fields.
     & THETA, Q, QCL, QCF,RHO                                           &
     &,U, V, W, W_ADV,P                                                 &

     &,PSTAR, EXNER_RHO_LEVELS,EXNER_THETA_LEVELS                       &
     &,LAND, P_THETA_LEVELS                                             &
! variables for subgrid turbulence scheme
     &,visc_BL_m, FM_3D, FH_3D, L_subfilter_vert, L_subfilter_horiz     &
     &,L_subfilter_blend,max_diff, turb_startlev_vert,turb_endlev_vert  &
     &,BL_COEF_KM,BL_COEF_KH                                            &

! ancillary fields and fields needed to be kept from timestep to
! timestep

     &,land_index, land_ice_index, soil_index, CANOPY_WATER             &
     &,SNODEP, THERM_COND, THERM_CAP, VOL_SMC_CRIT                      &
     &,VOL_SMC_WILT, VOL_SMC_SAT, STHF, STHU                            &
     &,OROG_SIL,OROG_HO2,ICE_THICKNESS,ICE_FRACTION                     &
     &,U_SEA,V_SEA,U_0_P,V_0_P,CCA,CCB                                  &
     &,CCT, CCLWP, CCW_RAD, LCBASE, DEEP_SOIL_TEMP, p_ti, TI            &
     &,TSTAR,Z0, p_ice_fract, p_ice_thick                               &
     &,SAT_SOIL_COND,SAT_SOILW_SUCTION,CLAPP_HORN                       &
     &,SMCL, T1_SD, Q1_SD, ZH                                           &
     &,CF_AREA, CF_BULK, CF_LIQUID                                      &
     &,CF_FROZEN, ls_rain, ls_snow, micro_tends                         &
     &,photosynth_act_rad, rad_hr, surf_radflux, SOIL_CLAY              &
     &,SOIL_SILT,SOIL_SAND, DUST_MREL1,DUST_MREL2                       &
     &,DUST_MREL3,DUST_MREL4,DUST_MREL5,DUST_MREL6                      &
     &,SO2_HILEM, SO2_EM, NH3_EM, DMS_EM                                &
     &,SOOT_HILEM, SOOT_EM, OCFF_HILEM, OCFF_EM, CO2_EMITS, CO2FLUX     &

! tracer fluxes - kdcorbin, 05/10
     &, TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4          &
     &, TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8          &
     &, TRACER_FLUX9, TRACER_FLUX10,TRACER_FLUX11,TRACER_FLUX12         &
     &, TRACER_FLUX13,TRACER_FLUX14,TRACER_FLUX15,TRACER_FLUX16         &
     &, TRACER_FLUX17,TRACER_FLUX18,TRACER_FLUX19,TRACER_FLUX20         &
! CO2 global emissions for mass balance - kdcorbin, 05/10               &
     &, CO2EMITMASS                                                     &
! rml 1/7/13 flag for co2 flux into passive tracer
     &, L_CO2_TRACER                                                    &

! IN/OUT
     &,theta_star,q_star,qcl_star,qcf_star,cf_star,cfl_star,cff_star    &
     &,R_u, R_v, R_w                                                    &
     &,net_flux,net_mflux,murk,tracer,tracer_ukca                       &
     &, DUST_DIV1, DUST_DIV2, DUST_DIV3                                 &
     &, DUST_DIV4, DUST_DIV5, DUST_DIV6                                 &
     &, so2, dms, so4_aitken, so4_accu                                  &
     &, so4_diss, nh3, soot_new                                         &
     &, soot_agd, soot_cld, bmass_new                                   &
     &, bmass_agd, bmass_cld, ocff_new, ocff_agd, ocff_cld, co2         &
! IN/OUT STPH_RP
     &,G0_RP,par_mezcla                                                 &

! IN/OUT River routing
     &, TOT_SURFROFF, TOT_SUBROFF                                       &
!
! OUT Fields
     &, rho_km, cH, NTML, CUMULUS, NBDSC, NTDSC                         &
     &, rhcpt, rhc_row_length, rhc_rows                                 &

! Additional variables for MOSES II
     &, FRAC_TYP, DISTURB_VEG, CANHT_PFT, LAI_PFT                       &
     &, CAN_WATER_TILE, CATCH_TILE, CATCH_SNOW                          &
     &, SNOW_GRND, SNODEP_TILE, Z0_TILE, TSTAR_TILE                     &
     &, INFIL_TILE, RGRAIN_TILE, CS, GS                                 &
     &, co2_dim_row, co2_dim_len, l_neg_tstar, l_snow_albedo            &
     &, l_phenol, l_triffid, l_trif_eq, l_q10, A_INTHD(23)              &
     &, STEPim(atmos_im),phenol_period,A_INTHD(22),CAN_MODEL            &
     &, G_LF_PFT_ACC, G_PHLF_PFT_ACC, NPP_PFT_ACC                       &
     &, RSP_W_PFT_ACC, RSA                                              &
     &, land_pts_trif, npft_trif, dolr, LW_down, SW_tile                &
     &,FRAC_LAND,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE                        &
     &,SOIL_ALB,cos_zenith_angle,can_rad_mod,ilayers                    &
     &, RADNET_TILE                                                     &
!glr cable start-----------------
     &, surf_down_sw,alb_tile,l_tile_pts                                &
     &, lat,cos_theta_longitude,day,time_sec,SW_DOWN                    &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L,SNOW_RHO3L        &
     &, SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE,TSOIL_TILE,T_SURF_TILE &
     &, HCONS,SOIL_TYPE,VEG_TYPE                                        &
     &, SNOW_FLG3L,total_nsteps                                         &
     &           ,FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB               &
     &           ,TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB           &
     &           ,USTAR_CAB,SURF_HTF_CAB                                &
! Lestevens - need first_atmstep_call for wblake fix
     &, l_cable,wblake_ratio,WB_LAKE                                    &
     &, TOT_WBLAKE,TOT_SUBRUN                                           &
     &           ,TOT_ALB                                               &
     &           ,U_S_CAB,CH_CAB,CD_CAB                                 &
     &           ,CD,CHX                                                &
     &           ,TILE_PTS,TILE_INDEX                                   &  
     &, SNOW_AGE,RTSOIL_TILE                                            &
     &, GFLUX_TILE,SGFLUX_TILE                                          &
     &, TRANSP_TILE                                                     &
     &, CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER                     &
     &, NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE                           &
!glr cable end -----------------
 
! Additional variables required for large-scale hydrology:
     &, L_TOP,L_PDM,FEXP,GAMMA_INT,TI_MEAN,TI_SIG                       &
     &, FSFC_SAT,F_WETLAND,WATER_TABLE,STHZW                            &
     &, A_FSAT,C_FSAT,A_FWET,C_FWET,L_SOIL_SAT_DOWN                     &
! Cariolle ozone 
     &, OZONE_TRACER                                                    &

! SCM and idealised UM surface forcing parameters
     &,L_flux_bc,flux_e,flux_h,L_spec_z0,z0m_scm,z0h_scm                &
!
! error information
     &,                      ErrorStatus,                               &
      !end step of experiment
     endstep, mype)

! DEPENDS ON: timer
          If (Ltimer) Call timer ('Atmos_Physics2',4)

          If(.NOT. (L_mix_ratio .and. L_moist_nonhydro_conserve)        &
     &            .and. CycleNo == NumCycles)then
            Deallocate(mix_graup_phys2)
            Deallocate(mix_rain_phys2)
            Deallocate(mix_cf2_phys2)
            Deallocate(mix_cf_phys2)
            Deallocate(mix_cl_phys2)
            Deallocate(mix_v_phys2)
          End If

!-----------------------------
!Tracer/CO2 Mass Fixer
!kdcorbin, 05/10
!-----------------------------

if (l_co2_mass .or. l_tracer_mass) then
! DEPENDS ON: tracer_massfix
  call tracer_massfix(row_length,rows,model_levels                    &
             ,tr_levels,tr_vars,i_tracermass_start                    &
             ,halo_i, halo_j, offx, offy, mype                        &
             ,timestep,timestep_number                                &
             ,r_theta_levels, r_rho_levels, exner_theta_levels        &
             ,FV_cos_theta_latitude, delta_lambda, delta_phi          &
             ,RHO, q, qcl, qcf                                        &
             ,L_CO2_MASS,L_TRACER_MASS,CO2, tracer               &
             ,co2emitmass,tracer_flux1,tracer_flux2,tracer_flux3        &
             ,tracer_flux4,tracer_flux5,tracer_flux6,tracer_flux7     &
             ,tracer_flux8,tracer_flux9,tracer_flux10,tracer_flux11   &
             ,tracer_flux12,tracer_flux13,tracer_flux14,tracer_flux15 &
             ,tracer_flux16,tracer_flux17,tracer_flux18,tracer_flux19 &
             ,tracer_flux20,tmass)

endif 

!-----------------------
!Methane Loss via OH
!kdcorbin, 05/10
!-----------------------

if (tr_vars > 0) then
!if (l_methane_loss .or. l_mcf_loss) then
   allocate(lossrate(row_length,rows,tr_levels,tr_vars))
   lossrate=0.
endif

if (l_methane_loss) then

  if (mype .eq. 0) then
     write(6,*) 'Calculating Methane Loss via OH, O1D and CL'
  endif

  if ((i_methane_tracers .lt. 1) .or.  &
     (i_methane_tracers .gt. tr_vars)) then
      write(6,*) '**********************************'
      write(6,*) 'Methane Specification Error: '
      write(6,*) '  Incorrect assignment of methane tracers'
      write(6,*) 'Specify I_METHANE_TRACERS in CNTLATM'   
      write(6,*) 'Stopping model'
      stop
   endif

! DEPENDS ON: tracer_methaneloss
  call tracer_methaneloss(row_length,rows,model_levels,i_methane_tracers  &
             ,tr_levels,tr_vars,offx,offy,timestep       &
             ,theta,p,rho,tracer,oh,ho2,lossrate)
endif  !methane loss

!-----------------------
!Radon Decay
!kdcorbin, 05/10
!-----------------------
if (l_radon_decay) then
 
   if (mype .eq. 0) then
      write(6,*) 'Calculating Radon Decay'
   endif

  if ((i_radon_tracernumber .lt. 1) .or.  &
      (i_radon_tracernumber .gt. tr_vars)) then
      write(6,*) '**********************************'
      write(6,*) 'Radon Specification Error: '
      write(6,*) '  Incorrect assignment of radon tracer'
      write(6,*) 'Specify I_RADON_TRACERNUMBER in CNTLATM'   
      write(6,*) 'Stopping model'
      stop
   endif

! DEPENDS ON: tracer_radondecay
  call tracer_radondecay(row_length,rows,tr_levels,tr_vars  &
          ,i_radon_tracernumber,offx,offy,timestep,tracer)

endif

!-----------------------
!MCF Loss
!kdcorbin, 05/10
!-----------------------
if (l_mcf_loss) then

   if (mype .eq. 0) then
      write(6,*) 'Calculating MCF Loss'
   endif

  if ((i_mcf_tracernumber .lt. 1) .or.  &
      (i_mcf_tracernumber .gt. tr_vars)) then
      write(6,*) '**********************************'
      write(6,*) 'MCF Specification Error: '
      write(6,*) '  Incorrect assignment of MCF tracer'
      write(6,*) 'Specify I_MCF_TRACERNUMBER in CNTLATM'   
      write(6,*) 'Stopping model'
      stop
   endif

! DEPENDS ON: tracer_mcfloss
  call tracer_mcfloss(row_length,rows,model_levels   &
         ,i_mcf_tracernumber   &
         ,tr_levels,tr_vars,offx,offy,timestep       &
         ,halo_i,halo_j,r_theta_levels,r_rho_levels  &
         ,theta,p,rho,tracer,oh,h2o2_limit,o3_chem,lossrate)

endif

!------------------------
!Tracer Information Printout
!kdcorbin, 05/10
!-----------------------
if (tr_vars > 0) then
!if (l_methane_loss .or. l_mcf_loss) then
   !DEPENDS ON: tracer_massprint
   call tracer_massprint(row_length,rows,model_levels                    &
             ,tr_levels,tr_vars                                       &
             ,halo_i, halo_j, offx, offy, mype                        &
             ,timestep, timestep_number                               &
             ,r_theta_levels, r_rho_levels, exner_theta_levels        &
             ,FV_cos_theta_latitude, delta_lambda, delta_phi          &
             ,RHO, q, qcl, qcf, tracer                                &
             ,lossrate,tracer_flux1,tracer_flux2,tracer_flux3        &
             ,tracer_flux4,tracer_flux5,tracer_flux6,tracer_flux7     &
             ,tracer_flux8,tracer_flux9,tracer_flux10,tracer_flux11   &
             ,tracer_flux12,tracer_flux13,tracer_flux14,tracer_flux15 &
             ,tracer_flux16,tracer_flux17,tracer_flux18,tracer_flux19 &
             ,tracer_flux20)

   deallocate(lossrate)
endif


! ----------------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after  Atmos_Physics2 ***'
        ENDIF ! L_print_pe .or. mype ==0
#include "princnorm.h"
      endif !  L_print_L2norms
! ---------------------------------------------------------------

! Apply diagnostics only at last cycle.
         If ( CycleNo == NumCycles ) Then

! Diagnostics STASHed for each section in Atmos_Physics2:

! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
      IF(L_bl)then
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,3,STASHwork3,                            &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork3)
      endif

          IF(convection_option == 3) THEN
! DEPENDS ON: timer
            IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
            CALL STASH(a_sm,a_im,5,STASHwork5,                          &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
            IF (Ltimer) CALL TIMER('STASH',4)
            DEALLOCATE (STASHwork5)
          ENDIF   ! on convection_option

      If (L_hydrology)then
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,8,STASHwork8,                            &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork8)
      endif
!
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,9,STASHwork9,                            &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork9)
!
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,19,STASHwork19,                          &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork19)
!
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,26,STASHwork26,                          &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)
          DEALLOCATE (STASHwork26)
!
          End If  ! CycleNo == NumCycles

          If (.NOT. L_run_with_physics2 ) Then
            L_physics=L_physics_store
          End If

!L --------- UM Section 18---- Data Assimilation --------------------
! DEPENDS ON: timer
          IF (LTimer) CALL TIMER('AS18 Assimilation',5)

          IF(L_AC .AND. LASSIMILATION .AND. ErrorStatus == 0)THEN
! Do AC assimilation

! copy non halo values of _star moisture variables into WORK arrays
! for passing to AC_CTL where they can be interpreted more
! conveniently as 2-d arrays (and halos not required)
            ALLOCATE ( STASHwork18(STASH_maxlen(18,A_im)) )
            ALLOCATE ( work_q(row_length,rows,wet_levels) )
            ALLOCATE ( work_qcl(row_length,rows,wet_levels) )
            ALLOCATE ( work_qcf(row_length,rows,wet_levels) )
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  work_q(i,j,k)   = q_star(i,j,k)
                  work_qcl(i,j,k) = qcl_star(i,j,k)
                  work_qcf(i,j,k) = qcf_star(i,j,k)
                End Do
              End Do
            End Do

! DEPENDS ON: timer
            IF(LTIMER)  CALL TIMER('AC_CTL',3)

! DEPENDS ON: ac_ctl
            CALL AC_CTL(STASH_MAXLEN(18,atmos_im),theta_field_size,     &
     &                  wet_levels, model_levels,                       &
     &                  theta_star, work_q, work_qcl, work_qcf,         &
     & OBS_FLAG,OBS,obs_flag_len,obs_len,                               &
     & p, p_theta_levels, exner_theta_levels,                           &
     & r_theta_levels, FV_cos_theta_latitude,                           &
     & cf_area, cf_bulk, cf_liquid, cf_frozen,                          &
     & pstar, ntml, cumulus,                                            &
     & STASHwork18,                                                     &
#include "argduma.h"
#include "argsts.h"
! #include "arglndm.h" 
#include "argppx.h"
     &                  l_mr_acctl,ErrorStatus,CMessage)
! DEPENDS ON: timer
            IF(LTIMER) CALL TIMER('AC_CTL',4)

! DEPENDS ON: timer
            IF(LTIMER) CALL TIMER('STASH',3)

! moved from out of AC_CTL
! DEPENDS ON: stash
        CALL STASH(a_sm, a_im, 18, STASHWORK18,                         &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &                               ErrorStatus,CMessage)

! DEPENDS ON: timer
            IF(LTIMER) CALL TIMER('STASH',4)

! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

! copy back new values of _star moisture variables
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  q_star  (i,j,k) = work_q(i,j,k)
                  qcl_star(i,j,k) = work_qcl(i,j,k)
                  qcf_star(i,j,k) = work_qcf(i,j,k)
                End Do
              End Do
            End Do
            DEALLOCATE (work_q)
            DEALLOCATE (work_qcl)
            DEALLOCATE (work_qcf)
            DEALLOCATE (STASHwork18)
          END IF      ! LASSIMILATION and L_AC

! DEPENDS ON: timer
          IF (LTimer) CALL TIMER('AS18 Assimilation',6)

          If (.not. L_use_bl_diag_term) Then
! zero ch term
            Do k = 1, model_levels-1
              Do j = 1-offy, rows+offy
                Do i = 1-offx, row_length+offx
                  cH(i,j,k) = 0.
                End Do
              End Do
            End Do
            Do k = 0, bl_levels-1
              Do j = 1-offy, rows+offy
                Do i = 1-offx, row_length+offx
                  rho_km(i,j,k) = 0.
                End Do
              End Do
            End Do
          End If    ! L_use_bl_diag_term

       Else    ! L_physics .false.

! zero ch term
          Do k = 1, model_levels-1
            Do j = 1-offy, rows+offy
              Do i = 1-offx, row_length+offx
                cH(i,j,k) = 0.
              End Do
            End Do
          End Do
! zero rho_km term
          Do k = 0, bl_levels-1
            Do j = 1-offy, rows+offy
              Do i = 1-offx, row_length+offx
                rho_km(i,j,k) = 0.
              End Do
            End Do
          End Do

! Include friction terms

! Run with a positive timestep if integrating backwards.
          IF (L_Backwards) timestep = pos_timestep
!  When no physics applied then may need a simple friction
!   if non-inviscid (idealised_problem) then DO NOT apply friction

          if(problem_number  /=  idealised_problem) Then
!  standard simple friction being used
            if(problem_number  ==  dynamical_core) Then
! DEPENDS ON: idl_friction_suarez_held
              Call IDL_Friction_Suarez_Held(                            &
     &                         row_length, rows, n_rows                 &
     &,                        model_levels, timestep                   &
     &,                        model_domain                             &
     &,                        offx, offy, at_extremity                 &
     &,                        friction_level                           &
     &,                        base_frictional_timescale                &
     &,                        SuHe_sigma_cutoff, SuHe_fric             &
     &,                        P, PSTAR                                 &
     &,                        U, V                                     &
     &,                        R_u, R_v )
            Else     ! problem_number  /=  dynamical_core
              Do k = 1, bl_levels
                Do j = j_begin, j_end
                  Do i = 1, row_length
                    R_u(i,j,k) = R_u(i,j,k) - timestep *                &
     &                     friction_level(k) * u(i,j,k)
                  End Do
                End Do
                Do j = 1, n_rows
                  Do i = 1, row_length
                    R_v(i,j,k) = R_v(i,j,k) - timestep *                &
     &                     friction_level(k) * v(i,j,k)
                  End Do
                End Do
              End Do

            EndIf    ! problem_number  ==  dynamical_core

          EndIf    ! problem_number  /=  idealised_problem

! Go back to negative timestep if integrating backwards.
          IF (L_Backwards) timestep = neg_timestep

        End If !   ! L_physics

! remove any instability
        If (L_adjust_wet .or.                                           &
     &      L_dry .or.                                                  &
     &      ( .not. L_physics .and.                                     &
     &        (.not. L_idealised_data .or. problem_number  ==  1        &
     &        .or. problem_number  ==  2) ) )Then

! DEPENDS ON: timer
          If (Ltimer) Call timer ('dry_static_adj',3)
! DEPENDS ON: dry_static_adj
          Call dry_static_adj(                                          &
     &                       theta_star, q_star, epsilon,               &
     &                       rows, row_length, model_levels,            &
     &                       wet_levels, offx, offy,                    &
     &                       Instability_diagnostics,                   &
     &                       .false.)
! DEPENDS ON: timer
          If (Ltimer) Call timer ('dry_static_adj',4)

        End If   ! L_adjust_wet .or. L_dry ....

      End If   ! ErrorStatus = 0

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS3 Atmos_Phys2',6)

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS3 Diffusion',5)

      If (ErrorStatus  ==  0 ) Then

        If ( L_diff_ctl ) then

! Apply diagnostics only at last cycle.
         If ( CycleNo == NumCycles ) Then

! Save current values to form diagnostics increments over diffusion
          IF( sf(185,13)) THEN
            ALLOCATE ( u_incr_diagnostic(row_length,rows,model_levels) )
! Hold u increment (no halos needed)
            DO k=1,model_levels
              DO j=1,rows
                DO i=1,row_length
                  u_incr_diagnostic(i,j,k) =R_u(i,j,k)
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDIF                   ! on STASHflag

          IF( sf(186,13) ) THEN
            ALLOCATE ( v_incr_diagnostic(row_length,n_rows,             &
     &                                              model_levels) )
! Hold v increment (no halos needed)
            DO k=1,model_levels
              DO j=1,n_rows
                DO i=1,row_length
                  v_incr_diagnostic(i,j,k) =R_v(i,j,k)
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDIF                    ! on STASHflag

          IF( sf(181,13) ) THEN
            ALLOCATE ( T_incr_diagnostic(row_length,rows,model_levels) )
! note at this point theta_star holds all increments to theta
            DO k=1,model_levels
              DO j=1,rows
                DO i=1,row_length
                  T_incr_diagnostic(i,j,k) =theta_star(i,j,k)
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDIF                    ! on STASHflag

          IF( sf(182,13) ) THEN
            ALLOCATE ( q_incr_diagnostic(row_length,rows,wet_levels) )
! note at this point q_star holds all increments to q
            DO k=1,wet_levels
              DO j=1,rows
                DO i=1,row_length
                  q_incr_diagnostic(i,j,k) =q_star(i,j,k)
                ENDDO ! i
              ENDDO ! j
            ENDDO ! k
          ENDIF                    ! on STASHflag

         End If  ! CycleNo == NumCycles
          IF ( CycleNo == 1 ) Then
          IF( sf(201,13) ) THEN
            ALLOCATE ( w_local_mask(row_length,rows) )
          ENDIF                    ! on STASHflag

          Endif  ! CycleNo == 1

          If (L_subfilter_horiz) then
!
! visc_m is currently S * lengthscale^2. Now multiply this by
! stability fn. Both FH and FM are on w points at this point.
!
            Do k = 1, model_levels - 2

              If (k >= turb_startlev_horiz .AND.                        &
     &               k <= turb_endlev_horiz) Then

                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = visc_m(i,j,k)*FH_3D(i,j,k+1)
                    visc_m(i,j,k) = visc_m(i,j,k)*FM_3D(i,j,k+1)
                    visc_h(i,j,k) = min(visc_h(i,j,k), max_diff)
                    visc_m(i,j,k) = min(visc_m(i,j,k), max_diff)
                  End Do
                End Do

              Else

                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = 0.0
                    visc_m(i,j,k) = 0.0
                  End Do
                End Do

              End If

            End Do
!
! Set visc_m and visc_h at the top two levels
!
            Do k = model_levels - 1, model_levels

              If (turb_startlev_horiz <= k .AND.                        &
     &            turb_endlev_horiz >= k) Then
     
                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = visc_h(i,j,model_levels-2)
                    visc_m(i,j,k) = visc_m(i,j,model_levels-2)
                  End Do
                End Do

              Else

                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = 0.0
                    visc_m(i,j,k) = 0.0
                  End Do
                End Do

              End If

            End Do

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(visc_m,                                    &
     &                       row_length, rows, model_levels,            &
     &                       halo_i, halo_j, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(visc_h,                                    &
     &                       row_length, rows, model_levels,            &
     &                       halo_i, halo_j, fld_type_p, .false.)

            If (model_domain  ==  mt_lam .or. model_domain  ==          &
     &                                          mt_cyclic_lam ) Then

! DEPENDS ON: fill_external_halos
              Call Fill_external_halos (visc_m,                         &
     &                                  row_length, rows, model_levels, &
     &                                  halo_i, halo_j)

! DEPENDS ON: fill_external_halos
              Call Fill_external_halos (visc_h,                         &
     &                                  row_length, rows, model_levels, &
     &                                  halo_i, halo_j)

            End If

          Else If (L_subfilter_vert) then

! visc_m is currently S * lengthscale^2. Now multiply this by
! stability fn. Both FH and FM are on w points at this point.
!
            Do k = 1, model_levels - 2

              If (k >= turb_startlev_vert .AND.                        &
     &               k <= turb_endlev_vert) Then
! use stability functions from BL scheme

                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = visc_m(i,j,k)*FH_3D(i,j,k+1)
                    visc_m(i,j,k) = visc_m(i,j,k)*FM_3D(i,j,k+1)
                  End Do
                End Do

              Else

                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = 0.0
                    visc_m(i,j,k) = 0.0
                  End Do
                End Do

              End If

            End Do
!
! Set visc_m and visc_h at the top two levels
!
            Do k = model_levels - 1, model_levels
              If (turb_startlev_vert <= k .AND.                         &
     &            turb_endlev_vert >= k) Then 
                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = visc_h(i,j,model_levels-2)
                    visc_m(i,j,k) = visc_m(i,j,model_levels-2)
                  End Do
                End Do
              Else
                Do j = 1, rows
                  Do i = 1, row_length
                    visc_h(i,j,k) = 0.0
                    visc_m(i,j,k) = 0.0
                  End Do
                End Do
              End If
            End Do


! DEPENDS ON: swap_bounds
            Call Swap_Bounds(visc_m,                                    &
     &                       row_length, rows, model_levels,            &
     &                       halo_i, halo_j, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(visc_h,                                    &
     &                       row_length, rows, model_levels,            &
     &                       halo_i, halo_j, fld_type_p, .false.)

            If (model_domain == mt_lam .or. model_domain ==             &
     &                                        mt_cyclic_lam ) Then     

! DEPENDS ON: fill_external_halos
              Call Fill_external_halos (visc_m,                         &
     &                                  row_length, rows, model_levels, &
     &                                  halo_i, halo_j)

! DEPENDS ON: fill_external_halos
              Call Fill_external_halos (visc_h,                         &
     &                                  row_length, rows, model_levels, &
     &                                  halo_i, halo_j)

            End If


          End If  ! L_subfilter_horiz or L_subfilter_vert

! DEPENDS ON: ni_diff_ctl
          Call NI_Diff_Ctl(                                             &
     &                     L_diffusion, L_cdiffusion, L_subfilter_horiz,&
     &                     L_vertical_diffusion, L_divdamp,             &
     &                     L_ramp, ramp_lat_radians,                    &
     &                     L_Backwards, Ltimer,                         &
     &                     timestep, pos_timestep, neg_timestep,        &
     &                     THETA, W, Q,                                 &
     &                     QCL, QCF,                                    &
     &                     QCF2, QRAIN, QGRAUP,                         &
     &                     mix_v, mix_cl, mix_cf,                       &
     &                     mix_cf2, mix_rain, mix_graup,                &
     &                     L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,       &
     &                     U, V, RHO,                                   &
     &                     EXNER_RHO_LEVELS, W_ADV,                     &
     &                     r_theta_levels, r_rho_levels, r_at_u, r_at_v,&
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     sec_theta_latitude,                          &
     &                     cos_theta_latitude, sec_v_latitude,          &
     &                     FV_sec_theta_latitude, cos_v_latitude,       &
     &                     sin_theta_latitude, sin_v_latitude, pi,      &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity, gc_proc_row_group,             &
     &                     mype, nproc, nproc_x, nproc_y, neighbour,    &
     &                     delta_lambda, delta_phi, L_regular,          &
     &                     glambda_p(lambda_start), phi_p,              &
     &                     glambda_u(lambda_start), phi_v,              &
     &                     rows, row_length, n_rows,                    &
     &                     model_levels, wet_levels, model_domain,      &
     &                     global_row_length, global_rows,              &
     &                     diffusion_coefficient_thermo,                &
     &                     diffusion_coefficient_w,                     &
     &                     diffusion_coefficient_q,                     &
     &                     diffusion_coefficient_wind,                  &
     &                     diffusion_order_thermo,                      &
     &                     diffusion_order_w,                           &
     &                     diffusion_order_q,                           &
     &                     diffusion_order_wind,                        &
     &                     visc_m, visc_h,                              &
     &                     horizontal_level, tar_horizontal,            &
     &                     level_start_wind, level_stop_wind,           &
     &                     level_start_q, level_stop_q,                 &
     &                     level_start_theta, level_stop_theta,         &
     &                     L_tardiff_q, w_conv_limit, tardiffq_factor,  &
     &                     tardiffq_test, tardiffq_start, tardiffq_end, &
     &                     sf(201,13), w_local_mask,                    &
     &                     L_adjust_theta,                              &
     &                     adjust_theta_start, adjust_theta_end,        &
     &                     L_vdiff_uv, vdiffuv_test, epsilon,           &
     &                     vdiffuv_factor, vdiffuv_start, vdiffuv_end,  &
     &                     vert_diffusion_coeff_wind,                   &
     &                     vert_diffusion_coeff_q,                      &
     &                     vert_diffusion_coeff_theta,                  &
     &                     div_damp_coefficient,                        &
     &                     theta_star, R_w, q_star,                     &
     &                     qcl_star, qcf_star,                          &
     &                     qcf2_star, qrain_star, qgraup_star,          &
     &                     R_u, R_v,                                    &
     &                     L_mix_ratio)

! ----------------------------------------------------------------------
! Section 13.1 Diagnostics at from diffusion and divergence damping
! ----------------------------------------------------------------------
! Apply diagnostics only at last cycle.
        If ( CycleNo == NumCycles ) Then

! DEPENDS ON: timer
          If (Ltimer) Call timer ('Diag_dif',3)

! section 13:
          IF( SF(0,13)                                                  &
                               ! Diagnostics required for this section
     &        .AND. ErrorStatus == 0) THEN

! Allocate diagnostic space for STASH
            IF(.NOT. L_Filter)then
            ALLOCATE (STASHwork13(STASH_maxlen(13,A_im)))
            ENDIF

! DEPENDS ON: diagnostics_dif
            CALL Diagnostics_dif(                                       &
     &         row_length, rows, n_rows, model_levels, wet_levels,      &
! primary  fields:
     &         THETA, Q,                                                &
! wind field increments after  dif :
     &         R_u, R_v,                                                &
! wind field increments before diffusion (on stashflag):
     &         u_incr_diagnostic, v_incr_diagnostic,                    &
     & T_incr_diagnostic,q_incr_diagnostic,                             &
! variables for subgrid turbulence scheme
     &         L_subfilter_horiz, L_subfilter_vert, visc_m, visc_h,     &
     &         shear, RNEUTML, FM_3D, FH_3D,                            &
     &         BL_LEVELS,BL_COEF_KM,BL_COEF_KH,                         &
     &         w_local_mask,                                            &
     &         theta_star, q_star,                                      &
     &         EXNER_THETA_LEVELS,                                      &
#include "argsts.h"
     &         STASHwork13)

! Tidy allocatable arrays
            IF( sf(185,13) ) THEN
              DEALLOCATE ( u_incr_diagnostic )
            ENDIF  ! on STASHflag

            IF( sf(201,13) ) THEN
              DEALLOCATE ( w_local_mask )
            ENDIF  ! on STASHflag

            IF( sf(186,13) ) THEN
              DEALLOCATE ( v_incr_diagnostic )
            ENDIF  ! on STASHflag

            IF( sf(181,13) ) THEN
              DEALLOCATE ( T_incr_diagnostic )
            ENDIF  ! on STASHflag

            IF( sf(182,13) ) THEN
              DEALLOCATE ( q_incr_diagnostic )
            ENDIF  ! on STASHflag

            IF( .NOT. ( (Model_domain==mt_global ) .AND.               &
     &                (L_polar_filter_incs .or. L_filter_incs) ) )THEN

! DEPENDS ON: timer
            IF (Ltimer) CALL TIMER('STASH',3)

! DEPENDS ON: stash
            CALL STASH(a_sm,a_im,13,STASHwork13,                        &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &      ErrorStatus,Cmessage)

! DEPENDS ON: timer
            IF (Ltimer) CALL TIMER('STASH',4)

            DEALLOCATE (STASHwork13)
            
          ENDIF ! .NOT. ( (Model_domain==mt_global ) etc

          ENDIF !   SF(0,13)

! DEPENDS ON: timer
          If (Ltimer) Call timer ('Diag_dif',4)

          End If ! CycleNo == NumCycles
        Endif      ! L_diff_ctl

        If ( CycleNo == NumCycles ) Then
          DEALLOCATE (BL_COEF_KM)
          DEALLOCATE (BL_COEF_KH)
          DEALLOCATE (visc_m)
          DEALLOCATE (visc_h)
          DEALLOCATE (visc_BL_m)
          DEALLOCATE (RNEUTML)
          DEALLOCATE (shear)
          If (L_Physics)then
            DEALLOCATE (FM_3D)
            DEALLOCATE (FH_3D)
          End If     ! L_physics
        End If    !  CycleNo == NumCycles

      End If       ! ErrorStatus  ==  0

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AS3 Diffusion',6)
!----------------------------------------------------------------------

      If ( L_tracer .and. CycleNo == NumCycles ) Then

! Obtain increments from Atmos_Physics2
! DEPENDS ON: tr_reset
        call TR_Reset(                                                  &
                      super_array_size, super_tracer_phys2,             &
                      L_CO2_interactive, CO2,                           &
                      L_Murk_advect, murk,                              &
                      L_Soot, soot_new,                                 &
                              soot_agd,                                 &
                              soot_cld,                                 &
                      L_SULPC_SO2, SO2,                                 &
                                   SO4_aitken,                          &
                                   so4_accu,                            &
                                   so4_diss,                            &
                      L_sulpc_nh3, nh3,                                 &
                      L_sulpc_dms, dms,                                 &
                      L_dust, DUST_DIV1,                                &
                              DUST_DIV2,                                &
                              DUST_DIV3,                                &
                              DUST_DIV4,                                &
                              DUST_DIV5,                                &
                              DUST_DIV6,                                &
                      L_biomass, bmass_new,                             &
                                 bmass_agd,                             &
                                 bmass_cld,                             &
                      L_ocff, ocff_new,                                 &
                              ocff_agd,                                 &
                              ocff_cld,                                 &
                      L_USE_CARIOLLE, OZONE_TRACER,                     &
                      tracer_phys2, tracer, tracer_ukca,                &
                      row_length, rows,                                 &
                      model_levels, tr_levels, tr_vars, tr_ukca,        &
                      offx, offy                                        &
                                            )

      end If ! L_tracer .and. CycleNo == NumCycles

! Obtain increments from Atmos_Physics2 and diffusion

      If (L_moist_nonhydro_conserve)then
! after atmos_phys2 q_star holds q_dash (latest estimate to q_(n+1))
        IF(l_mix_ratio)then
          allocate ( mix_v_inter(1-offx:row_length+offx,                &
     &                           1-offy:rows+offy, wet_levels) )
          allocate ( mix_cl_inter(1-offx:row_length+offx,               &
     &                            1-offy:rows+offy, wet_levels) )
          allocate ( mix_cf_inter(1-offx:row_length+offx,               &
     &                            1-offy:rows+offy, wet_levels) )
          if(L_mcr_qcf2)then
            allocate ( mix_cf2_inter(1-offx:row_length+offx,            &
     &                               1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_cf2_inter(1,1,1) )
          endif
          if(L_mcr_qrain)then
            allocate ( mix_rain_inter(1-offx:row_length+offx,           &
     &                                1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_rain_inter(1,1,1) )
          endif
          if(L_mcr_qgraup)then
            allocate ( mix_graup_inter(1-offx:row_length+offx,          &
     &                                 1-offy:rows+offy, wet_levels) )
          else
            allocate ( mix_graup_inter(1,1,1) )
          endif
! DEPENDS ON: q_to_mix
          call q_to_mix (row_length, rows, wet_levels,                  &
     &                   offx,offy     ,                                &
     &                   q_star, qcl_star, qcf_star,                    &
     &               qcf2_star, qrain_star, qgraup_star,                &
!    &               .false.,.false.,.false.,
     &                L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,            &
     &               mix_v_inter, mix_cl_inter, mix_cf_inter,           &
     &               mix_cf2_inter, mix_rain_inter, mix_graup_inter     &
     &                   )
          do k=1, wet_levels
            do j=1,rows
              do i=1,row_length
                mix_v_phys2(i,j,k) = mix_v_inter(i,j,k)                 &
     &                             - mix_v_phys2(i,j,k)
                mix_cl_phys2(i,j,k) = mix_cl_inter(i,j,k)               &
     &                              - mix_cl_phys2(i,j,k)
                mix_cf_phys2(i,j,k) = mix_cf_inter(i,j,k)               &
     &                              - mix_cf_phys2(i,j,k)
              end do
            end do
          end do
          if(L_mcr_qcf2)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  mix_cf2_phys2(i,j,k) = mix_cf2_inter(i,j,k)           &
     &                                 - mix_cf2_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_mcr_qrain)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  mix_rain_phys2(i,j,k) = mix_rain_inter(i,j,k)         &
     &                                  - mix_rain_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_mcr_qgraup)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  mix_graup_phys2(i,j,k) = mix_graup_inter(i,j,k)       &
     &                                   - mix_graup_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_pc2)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  cf_phys2(i,j,k) = cf_star(i,j,k)                      &
     &                            - cf_phys2(i,j,k)
                  cfl_phys2(i,j,k) = cfl_star(i,j,k)                    &
     &                             - cfl_phys2(i,j,k)
                  cff_phys2(i,j,k) = cff_star(i,j,k)                    &
     &                             - cff_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          deallocate (mix_v_inter)
          deallocate (mix_cl_inter)
          deallocate (mix_cf_inter)
          deallocate (mix_cf2_inter)
          deallocate (mix_rain_inter)
          deallocate (mix_graup_inter)

        else

          Do k = 1, wet_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_phys2(i,j,k) = q_star(i,j,k)                          &
     &                         - q_phys2(i,j,k)
                qcl_phys2(i,j,k) = qcl_star(i,j,k)                      &
     &                           - qcl_phys2(i,j,k)
                qcf_phys2(i,j,k) = qcf_star(i,j,k)                      &
     &                           - qcf_phys2(i,j,k)
              End Do
            End Do
          End Do
          if(L_pc2)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  cf_phys2(i,j,k) = cf_star(i,j,k)                      &
     &                            - cf_phys2(i,j,k)
                  cfl_phys2(i,j,k) = cfl_star(i,j,k)                    &
     &                             - cfl_phys2(i,j,k)
                  cff_phys2(i,j,k) = cff_star(i,j,k)                    &
     &                             - cff_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_mcr_qcf2)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  qcf2_phys2(i,j,k) = qcf2_star(i,j,k)                  &
     &                              - qcf2_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_mcr_qrain)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  qrain_phys2(i,j,k) = qrain_star(i,j,k)                &
     &                               - qrain_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
          if(L_mcr_qgraup)then
            do k=1, wet_levels
              do j=1,rows
                do i=1,row_length
                  qgraup_phys2(i,j,k) = qgraup_star(i,j,k)              &
     &                               -  qgraup_phys2(i,j,k)
                End Do
              End Do
            End Do
          endif
        endif !L_mix_ratio

      end if  !L_moist_nonhydro_conserve

! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after diffusion  ***'
        ENDIF ! L_print_pe .or. mype ==0
#include "princnorm.h"
      endif !  L_print_L2norms
! ----------------------------------------------------------------------
! Section 4.  Form and solve Helmholtz equation and update variables.
! ----------------------------------------------------------------------
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AS4 Solver',5)

      If (ErrorStatus  ==  0) then

! Apply diagnostics only at last cycle.
        If ( CycleNo == NumCycles ) Then

! Save current values to form diagnostics increments over Solver
        IF( sf(185,10)) THEN
          ALLOCATE ( u_incr_diagnostic(row_length,rows,model_levels) )
! Hold u increment (no halos needed)
          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                u_incr_diagnostic(i,j,k) =R_u(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                   ! on STASHflag

        IF( sf(186,10) ) THEN
          ALLOCATE ( v_incr_diagnostic(row_length,n_rows,model_levels) )
! Hold v increment (no halos needed)
          DO k=1,model_levels
            DO j=1,n_rows
              DO i=1,row_length
                v_incr_diagnostic(i,j,k) =R_v(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        IF( sf(187,10) ) THEN
          ALLOCATE ( w_incr_diagnostic(row_length,rows,model_levels) )
! Hold w increment (no halos needed)
          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                w_incr_diagnostic(i,j,k) =R_w(i,j,k)
              ENDDO ! i
            ENDDO ! j
          ENDDO ! k
        ENDIF                    ! on STASHflag

        End If ! CycleNo == NumCycles
! set halos for q_star, qcl_star, qcf_star and theta_star

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   q_star, row_length, rows, wet_levels,          &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(q_star,row_length,rows,                &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   qcl_star, row_length, rows, wet_levels,        &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(qcl_star,row_length,rows,              &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   qcf_star, row_length, rows, wet_levels,        &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(qcf_star,row_length,rows,              &
     &                           wet_levels,Offx,Offy)
        If(L_mcr_qcf2)then
! DEPENDS ON: swap_bounds
          CALL Swap_Bounds(                                             &
     &                   qcf2_star, row_length, rows, wet_levels,       &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
          CALL FILL_EXTERNAL_HALOS(qcf2_star,row_length,rows,           &
     &                           wet_levels,Offx,Offy)
        endif
        If(L_mcr_qrain)then
! DEPENDS ON: swap_bounds
          CALL Swap_Bounds(                                             &
     &                   qrain_star, row_length, rows, wet_levels,      &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
          CALL FILL_EXTERNAL_HALOS(qrain_star,row_length,rows,          &
     &                           wet_levels,Offx,Offy)
        endif
        If(L_mcr_qgraup)then
! DEPENDS ON: swap_bounds
          CALL Swap_Bounds(                                             &
     &                   qgraup_star, row_length, rows, wet_levels,     &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
          CALL FILL_EXTERNAL_HALOS(qgraup_star,row_length,rows,         &
     &                           wet_levels,Offx,Offy)
        endif


! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   theta_star, row_length, rows, model_levels,    &
     &                   offx, offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(theta_star,row_length,rows,            &
     &                           model_levels,Offx,Offy)

        If (L_mix_ratio) Then

! DEPENDS ON: q_to_mix_halo
        call q_to_mix_halo (row_length, rows, wet_levels,               &
     &                 offx,offy     ,                                  &
     &                 q_star, qcl_star, qcf_star,                      &
     &                 qcf2_star,qrain_star, qgraup_star,               &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 mix_v_star, mix_cl_star, mix_cf_star,            &
     &                 mix_cf2_star, mix_rain_star, mix_graup_star      &
     &                 )

        End If       !  L_mix_ratio

! Calculate updated values of pressure gradient terms.
! DEPENDS ON: timer
        If (Ltimer) Call timer ('pg_update',3)

! DEPENDS ON: ni_pg_update
        Call NI_pg_update(                                              &
     &                 Q, QCL, QCF,                                     &
     &                 QCF2, QRAIN, QGRAUP,                             &
     &                 q_star, qcl_star, qcf_star,                      &
     &                 qcf2_star, qrain_star, qgraup_star,              &
     &                 mix_v, mix_cl, mix_cf,                           &
     &                 mix_cf2, mix_rain, mix_graup,                    &
     &                 mix_v_star, mix_cl_star, mix_cf_star,            &
     &                 mix_cf2_star, mix_rain_star, mix_graup_star,     &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 theta_star, THETA, EXNER_RHO_LEVELS,             &
     &                 r_theta_levels, r_rho_levels,                    &
     &                 sec_theta_latitude,                              &
     &                 delta_lambda, delta_phi, timestep,               &
     &                 grecip_dlamp(lambda_start), recip_dphip,         &
     &                 wt_lambda_p, wt_phi_p,                           &
     &                 alpha_3, alpha_4,                                &
     &                 Cp, epsilon,                                     &
     &                 row_length, rows, n_rows, model_levels,          &
     &                 wet_levels, model_domain,                        &
     &                 first_constant_r_rho_level,                      &
     &                 offx, offy, halo_i, halo_j, at_extremity,        &
     &                 CycleNo, L_new_tdisc, L_qwaterload,              &
     &                 R_u, R_v, R_w, L_mix_ratio, L_regular )

! DEPENDS ON: timer
        If (Ltimer) Call timer ('pg_update',4)

! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after pg_update  ***'
        ENDIF ! L_print_pe .or. mype ==0
#include "princnorm.h"
      endif !  L_print_L2norms
! ---------------------------------------------------------------

        If (model_domain  ==  mt_global ) Then
          If ( L_polar_filter_incs .or. L_filter_incs) Then

! DEPENDS ON: ni_filter_incs_ctl
            Call NI_filter_incs_Ctl(                                    &
     &                      THETA, theta_star, R_u, R_v, R_w,           &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      r_theta_levels, r_rho_levels,               &
     &                      r_at_u, r_at_v,                             &
     &                      delta_lambda, delta_phi,                    &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      cos_v_latitude, sec_v_latitude,             &
     &                      polar_filter_north_lat_limit,               &
     &                      polar_filter_south_lat_limit,               &
     &                      polar_filter_coefficient,                   &
     &                      polar_filter_n_sweeps,                      &
     &                      polar_filter_step_per_sweep,                &
     &                      polar_filter_lat_limit,                     &
     &                      max_121_rows, u_sweeps, v_sweeps,           &
     &                      global_u_filter, global_v_filter,           &
     &                      u_begin, u_end, v_begin, v_end,             &
     &                      diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
     &                      diff_coeff_thermo, diff_coeff_wind,         &
     &                      diff_order_thermo, diff_order_wind,         &
     &                      horizontal_level, mype,                     &
     &                      offx, offy, halo_i, halo_j,                 &
     &                      nproc, nproc_x, nproc_y, datastart,         &
     &                      neighbour, at_extremity, model_domain,      &
     &                      gc_proc_row_group, gc_proc_col_group,       &
     &                      L_polar_filter_incs, L_filter_incs,         &
     &                      L_pfcomb, L_pftheta, L_pfuv, L_pfw,         &
     &                      L_pofil_new, L_diff_incs, Ltimer,           &
     &                      exner_theta_levels,                         &
#include "argsts.h"
     &                      STASHwork13)

            IF( SF(0,13) ) THEN  ! Diagnostics required for this section

! DEPENDS ON: timer
              IF (Ltimer) CALL TIMER('STASH',3)

! DEPENDS ON: stash
              CALL STASH(a_sm,a_im,13,STASHwork13,                      &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &        ErrorStatus,Cmessage)

! DEPENDS ON: timer
              IF (Ltimer) CALL TIMER('STASH',4)

              DEALLOCATE (STASHwork13)

            ENDIF !   SF(0,13)

          End If   ! L_polar_filter_incs .or. L_filter_incs

! calculate R_u at the poles.

! DEPENDS ON: polar_vector_wind_n
          Call Polar_vector_wind_n(                                     &
     &                       R_v,                                       &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, model_levels, mag_vector_np,       &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       offx, offy, global_row_length,             &
     &                       gc_proc_row_group, at_extremity)

          If (at_extremity(PSouth) ) Then
            Do k = 1,model_levels
              Do i = 1, row_length
                R_u(i,1,k) = - mag_vector_sp(k) *                       &
     &                         sin ( lambda_a(i) - dir_vector_sp(k) )
              End Do
            End Do
          End If    ! at_extremity(PSouth)
          If (at_extremity(PNorth) ) Then
            Do k = 1,model_levels
              Do i = 1, row_length
                R_u(i,rows,k) = mag_vector_np(k) *                      &
     &                         sin ( lambda_a(i) - dir_vector_np(k) )
              End Do
            End Do
          End If    ! at_extremity(PNorth)

        End If    !model_domain

! Do halo swops for * variables

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   R_u, row_length, rows, model_levels,           &
     &                   offx, offy, fld_type_u, .true.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(R_u,row_length,rows,                   &
     &                           wet_levels,Offx,Offy)

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   R_v, row_length, n_rows, model_levels,         &
     &                   offx, offy, fld_type_v, .true.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(R_v,row_length,n_rows,                 &
     &                           wet_levels,Offx,Offy)

      End If        !   ErrorStatus  ==  0

! ----------------------------------------------------------------------
! Section 4.1 Form and solve Helmholtz equation, returning all
!             corrections that are required.
! ----------------------------------------------------------------------
      If (ErrorStatus  ==  0) then

! in the idealised versions the boundary layer coefficients Km are
! replaced by the frictional timescale.

        If (model_domain ==  mt_lam) Then

          L_do_halos=.FALSE.
          L_do_boundaries=.TRUE.

          If (RIM_STEPSA  ==  0) Then
            increment_factor=0.0
          Else
            increment_factor=1.0/                                       &
     &      (RIM_STEPSA-MOD(Timestep_Number-1,RIM_STEPSA))
          End If

          lbc_size=LENRIMA(fld_type_p,halo_type_extended,               &
     &             rima_type_norm)

          If (L_lbc_new) Then
! Obtain Exner tendency to pass to the solver to apply lbc

! Apply additional balance to EXNER_LBC_TEND, RHO_LBC_TEND 
! and W_LBC_TEND only on timestep after up_bound has been called.

          If ( MOD(Timestep_Number-1,RIM_STEPSA) == 0) Then 
            If (L_LBC_balance) Then

! DEPENDS ON: balance_lbc_values
              CALL BALANCE_LBC_VALUES(                                  &
     &        EXNER_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND,             &
     &        Q_LBC_TEND, W_LBC_TEND, W_ADV_LBC_TEND,                   &
     &        U_LBC_TEND, V_LBC_TEND,                                   &
     &        R_RHO_LEVELS, R_THETA_LEVELS,                             &
     &        ROW_LENGTH, ROWS, model_levels, wet_levels, HALO_I,HALO_J,&
     &        LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),    &
     &        LENRIMA(fld_type_u,halo_type_extended,rima_type_norm),    &
     &        LENRIMA(fld_type_v,halo_type_extended,rima_type_norm),    &
     &        LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm),&
     &        LBC_STARTA(1,fld_type_u,halo_type_extended,rima_type_norm),&
     &        LBC_STARTA(1,fld_type_v,halo_type_extended,rima_type_norm),&
     &        RIMWIDTHA, N_RIMS_TO_DO, RIMWEIGHTSA, AT_EXTREMITY,       &
     &        DELTA_PHI, DELTA_LAMBDA,                                  &
     &        BASE_PHI, BASE_LAMBDA,                                    &
     &        TWO_OMEGA, DATASTART,                                     &
     &        LAT_ROT_NP,                                               &
     &        GLOBAL_ROW_LENGTH, GLOBAL_ROWS, L_int_uvw_lbc             &
     &        )

            End If  !  L_LBC_balance
          End If  !  MOD(Timestep_Number-1,RIM_STEPSA) == 0

          DO k = 1, MODEL_LEVELS
            DO i=1,lbc_size
              EXNER_LBC_REAL_TEND(i,k) = increment_factor *             &
     &                         ( EXNER_LBC_TEND(i,k) - EXNER_LBC(i,k) )
            END DO
          END DO

          Else  ! Original lbc lgorithm requires winds at boundaries

! If old lbcs set the outer boundaries of R_u, R_v and R_w to the 
! difference between the boundary at time level n+1 and time level n

          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              W_LBC_REAL_TEND(i,k) = increment_factor *                 &
     &             (W_LBC_TEND(i,k) - W_LBC(i,k))
            END DO
          END DO

          lbc_size=LENRIMA(fld_type_u,halo_type_extended,               &
     &             rima_type_norm)
          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              U_LBC_REAL_TEND(i,k) = increment_factor *                 &
     &             (U_LBC_TEND(i,k) - U_LBC(i,k))
            END DO
          END DO
          lbc_size=LENRIMA(fld_type_v,halo_type_extended,               &
     &             rima_type_norm)

          DO k=1,MODEL_LEVELS
            DO i=1,lbc_size
              V_LBC_REAL_TEND(i,k) = increment_factor *                 &
     &             (V_LBC_TEND(i,k) - V_LBC(i,k))
            END DO
          END DO
! U
! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH,ROWS,Offx,Offy,                                    &
     &    MODEL_LEVELS,fld_type_u,R_u,                                  &
     &    LENRIMA(fld_type_u,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_u,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_u,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    U_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm),n_rims_to_do,                       &
     &    RIMWEIGHTSA,AT_EXTREMITY,                                     &
     &    L_do_boundaries,L_do_halos)

! V
! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH,N_ROWS,Offx,Offy,                                  &
     &    MODEL_LEVELS,fld_type_v,R_v,                                  &
     &    LENRIMA(fld_type_v,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_v,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_v,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    V_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm),n_rims_to_do,                       &
     &    RIMWEIGHTSA,AT_EXTREMITY,                                     &
     &    L_do_boundaries,L_do_halos)

! W
! DEPENDS ON: set_lateral_boundaries
          CALL SET_LATERAL_BOUNDARIES(                                  &
     &    ROW_LENGTH,ROWS,0,0,                                          &
     &    MODEL_LEVELS,fld_type_p,R_w,                                  &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    W_LBC_REAL_TEND,                                              &
     &    RIMWIDTHA(rima_type_norm),n_rims_to_do,                       &
     &    RIMWEIGHTSA,AT_EXTREMITY,                                     &
     &    L_do_boundaries,L_do_halos)

          End If !  L_lbc_new

        Else If (model_domain  ==  mt_cyclic_lam) Then
! store old time level u and v in R_u and R_v at boundaries
          If (at_extremity(PSouth) ) Then
            Do k = 1, model_levels
              Do i = 1-offx, row_length+offx
                R_u(i,1,k)=U(i,1,k)
                R_v(i,1,k)=V(i,1,k)
              End Do
            End Do
            Do k = 1, model_levels
              Do i = 1, row_length
                R_w(i,1,k)=W(i,1,k)
              End Do
            End Do
          End If    ! at_extremity(PSouth)
          If (at_extremity(PNorth) ) Then
            Do k = 1, model_levels
              Do i = 1-offx, row_length+offx
                R_u(i,rows,k)= U(i,rows,k)
                R_v(i,n_rows,k)=V(i,n_rows,k)
              End Do
            End Do
            Do k = 1, model_levels
              Do i = 1, row_length
                R_w(i,rows,k)=W(i,rows,k)
              End Do
            End Do
          End If    ! at_extremity(PNorth)

! set increments to u, v and w at boundaries.
          If (at_extremity(PSouth) ) Then

            Do k = 1, model_levels
              Do i = 1-offx, row_length+offx
                R_u(i,1,k)=U(i,1,k)-R_u(i,1,k)
                R_v(i,1,k)=V(i,1,k)-R_v(i,1,k)
              End Do
            End Do
            Do k = 1, model_levels
              Do i = 1, row_length
                R_w(i,1,k)=W(i,1,k)-R_w(i,1,k)
              End Do
            End Do
          End If    ! at_extremity(PSouth)
          If (at_extremity(PNorth) ) Then
            Do k = 1, model_levels
              Do i = 1-offx, row_length+offx
                R_u(i,rows,k) = U(i,rows,k) - R_u(i,rows,k)
                R_v(i,n_rows,k) = V(i,n_rows,k) - R_v(i,n_rows,k)
              End Do
            End Do
            Do k = 1, model_levels
              Do i = 1, row_length
                R_w(i,rows,k) = W(i,rows,k) - R_w(i,rows,k)
              End Do
            End Do
          End If    ! at_extremity(PNorth)

        Endif  ! model_domain

! Note:
! R_u, R_v, R_w on output contain u_prime, v_prime, w_prime

! DEPENDS ON: timer
        If (Ltimer) Call timer ('PE_Helmholtz',3)

        ALLOCATE ( dtheta_dr_term(row_length,rows,model_levels) )


        If ( CycleNo == 1 ) Then
          GCR_zero_guess_it = GCR_zero_init_guess
        Else
          GCR_zero_guess_it = ( .not. L_GCR_cycle_opt )                 &
     &                          .and. GCR_zero_init_guess
          If ( .NOT. L_GCR_cycle_opt ) Then
              exner_prime(:,:,:) = 0.0
          End If
        End If

        If ( CycleNo == 1 ) Then
          If ( GCR_use_tol_abs ) Then
            GCR_run_tol_abs = GCR_tol_abs
            GCR_run_tol_res = 0.0
          Else
            GCR_run_tol_res = GCR_tol_res
            GCR_run_tol_abs = 0.0
          End If
        Else
          If ( GCR_use_tol_abs ) Then
            GCR_run_tol_abs = GCR_tol_abs2
            GCR_run_tol_res = 0.0
          Else
            GCR_run_tol_res = GCR_tol_res2
            GCR_run_tol_abs = 0.0
          End If
        End If

! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms before solver  ***'
        ENDIF ! L_print_pe .or. mype ==0
#include "princnorm.h"
      endif !  L_print_L2norms
! ---------------------------------------------------------------

! DEPENDS ON: ni_pe_helmholtz
        Call NI_PE_Helmholtz(                                           &
     &                      U, V, W,                                    &
     &                      r_theta_levels,                             &
     &                      r_rho_levels, P, RHO,                       &
     &                      rho_np1, THETA,                             &
     &                      theta_star,theta_np1,                       &
     &                      Q, QCL, QCF,                                &
     &                      QCF2, QRAIN, QGRAUP,                        &
     &                      mix_v, mix_cl, mix_cf,                      &
     &                      mix_cf2, mix_rain, mix_graup,               &
     &                      mix_v_star, mix_cl_star, mix_cf_star,       &
     &                      mix_cf2_star, mix_rain_star, mix_graup_star,&
     &                      mix_v_np1, mix_cl_np1, mix_cf_np1,          &
     &                      mix_cf2_np1, mix_rain_np1, mix_graup_np1,   &
     &                      q_star, q_np1, qcl_star, qcf_star,          &
     &                      qcf2_star, qrain_star, qgraup_star,         &
     &                      qcl_np1, qcf_np1,                           &
     &                      qcf2_np1 , qrain_np1, qgraup_np1,           &
     &                      rho_Km, cH, G_term_tol,                     &
     &                      EXNER_RHO_LEVELS,                           &
     &                      EXNER_THETA_LEVELS,                         &
     &                      frictional_timescale,                       &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      cos_v_latitude, sec_v_latitude,             &
     &                      FV_cos_theta_latitude,                      &
     &                      FV_sec_theta_latitude,                      &
     &                      f3_at_u, f3_at_v,                           &
     &                      timestep, timestep_number,                  &
     &                      row_length, rows, n_rows,                   &
     &                      model_levels, wet_levels,                   &
     &                      bl_levels, L_print_L2helm, L_print_pe,      &
     &                      L_flush6, norm_lev_start, norm_lev_end,     &
! ---------------------------------------------------------------
     &                      delta_lambda, delta_phi,                    &
     &                      glambda_p(lambda_start), phi_p,             &
     &                      glambda_u(lambda_start), phi_v,             &
     &                      gdlambda_p(lambda_start), dphi_p,           &
     &                      gdlambda_u(lambda_start), dphi_v,           &
     &                      grecip_dlamp(lambda_start), recip_dphip,    &
     &                      grecip_dlamu(lambda_start), recip_dphiv,    &
     &                      wt_lambda_p, wt_phi_p,                      &
     &                      wt_lambda_u, wt_phi_v,                      &
     &                      GCR_max_iterations, GCR_diagnostics,        &
     &                      GCR_its_switch(CycleNo), GCR_its_avg_step,  &
     &                      GCR_max_its(CycleNo), GCR_min_its(CycleNo), &
     &                      GCR_sum_its(CycleNo), GCR_max_time(CycleNo),&
     &                      GCR_min_time(CycleNo),                      &
     &  GCR_run_tol_res, GCR_run_tol_abs, GCR_use_tol_abs,              &
     &                      GCR_zero_guess_it,                          &
     &                      GCR_use_residual_Tol,                       &
     &                      GCR_adi_add_full_soln, L_gcr_fast_x,        &
     &                      GCR_precon_option, GCR_ADI_Pseudo_timestep, &
     &                      GCR_n_ADI_pseudo_timesteps,                 &
     &                      eta_theta_levels, eta_rho_levels,           &
     &                      alpha1, alpha2, alpha3, alpha4,             &
     &                      alpha_Cd, kappa, Cp, R, Pi,                 &
     &                      epsilon, model_domain, L_physics,           &
     &                      GCR_Restart_value,                          &
     &                      first_constant_r_rho_level,                 &
     &                      first_constant_r_rho_level_m1,              &
     &                      R_u, R_v, R_w, exner_prime, dtheta_dr_term, &
     &                      EXNER_LBC_REAL_TEND,                        &
     &                      LENRIMA(1,1,rima_type_norm),                &
     &                      LBC_SIZEA(1,1,1,rima_type_norm),            &
     &                      LBC_STARTA(1,1,1,rima_type_norm),           &
     &                      RIMWIDTHA(rima_type_norm), RIMWEIGHTSA,     &
!!! parallel variables
     &                      mype, nproc, nproc_x, nproc_y,              &
     &                      halo_i, halo_j, datastart,                  &
     &                      L_regular, at_extremity,                    &
     &                      n_rims_to_do, offx, offy,                   &
     &                      gc_proc_row_group, gc_proc_col_group,       &
     &                      global_row_length, global_rows,             &
     &                      g_rows, g_row_length, g_datastart, CycleNo, &
     &                      L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,      &
     &                      L_mix_ratio, L_fint_theta, L_new_tdisc,     &
     &                      L_qwaterload, L_lbc_new                     &
     &                      )

! DEPENDS ON: timer
        If (Ltimer) Call timer ('PE_Helmholtz',4)

! If increments are swapped now then this avoids needing swap
! in update rho

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   R_u, row_length, rows, model_levels,           &
     &                   offx, offy, fld_type_u, .true.)

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
     &                   R_v, row_length, n_rows, model_levels,         &
     &                   offx, offy, fld_type_v, .true.)


!----------------------------------------------------------------------
! Section 4.1 Diagnostics from Solver (may need to move to after
!  rho update if wanted to output rho increments in future).
!----------------------------------------------------------------------
! Apply diagnostics only at lat cycle.
        If ( CycleNo == NumCycles ) Then

! DEPENDS ON: timer
        if(Ltimer) Call timer ('Diag_sol',3)
!
        IF (SF(0,10)                                                    &
                         ! diagnostics required for this section
     &      .and. ErrorStatus == 0) THEN

! Allocate diagnostic space for STASH
          ALLOCATE (STASHwork10(STASH_maxlen(10,A_im)))

! DEPENDS ON: diagnostics_solver
          CALL Diagnostics_solver(                                      &
     &       row_length,rows,n_rows,model_levels,                       &
! wind field increments after  solver :
     &       R_u,R_v,R_w,                                               &
! wind field increments before solver (on stashflag):
     &       u_incr_diagnostic,v_incr_diagnostic,w_incr_diagnostic,     &
#include "argsts.h"
     &       STASHwork10)

! Tidy allocatable arrays
          IF( sf(185,10) ) THEN
            DEALLOCATE ( u_incr_diagnostic )
          ENDIF  ! on STASHflag
          IF( sf(186,10) ) THEN
            DEALLOCATE ( v_incr_diagnostic )
          ENDIF  ! on STASHflag
          IF( sf(187,10) ) THEN
            DEALLOCATE ( w_incr_diagnostic )
          ENDIF  ! on STASHflag

! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',3)

! DEPENDS ON: stash
          CALL STASH(a_sm,a_im,10,STASHwork10,                          &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)

          DEALLOCATE (STASHwork10)

! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('STASH',4)

        ENDIF !   SF(0,10)

! DEPENDS ON: timer
        if(Ltimer) Call timer ('Diag_sol',4)

        End If ! CycleNo == NumCycles

      End If  !      ErrorStatus  ==  0

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AS4 Solver',6)

! ---------------------------------------------------------------
!    diagnostic printing of l2norms
      if( L_print_L2norms ) then
        IF( L_print_pe .or. mype ==0 ) then
          write(6,*)' ***   L2 norms after solver  ***'
        ENDIF ! L_print_pe .or. mype ==0
#include "princnorm.h"
      endif !  L_print_L2norms
! ---------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 5.0 Calculate rho at new time level, using flux form of
!             equation.
! ----------------------------------------------------------------------
! Call timer for updates code
! DEPENDS ON: timer
      If (Ltimer .and. CycleNo == 1 ) CALL TIMER('AS5-8 Updates',5)

      If (ErrorStatus  ==  0) then

        If ( CycleNo == 1 ) Then
          If ( sf(188,30) ) Then
            ALLOCATE ( inc_rho( 1-offx:row_length+offx,                 &
     &                        1-offy:rows+offy, model_levels) )
          Else
            ALLOCATE ( inc_rho( 1,1,1) )
          End If
        End If

! stash diag 30188 only needed at last cycle.
        L_do_inc_rho = (CycleNo == NumCycles) .and. sf(188,30)

! Set rho_n at
        If ( CycleNo == NumCycles .and. ( L_tracer .or.                 &
     &       L_moist_nonhydro_conserve .or. L_do_inc_rho ) ) Then

! store rho_n before updating by flux_rho
! DEPENDS ON: copy_field
      CALL COPY_FIELD(RHO,rho_n                                         &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels, model_levels, 1, model_levels       &
     &,               Offx, Offy, Offx, Offy                            &
     &,               fld_type_p, .true., .false., .false.)

! DEPENDS ON: fill_external_halos
          CALL FILL_EXTERNAL_HALOS(rho_n,row_length, rows,              &
     &                     model_levels,offx,offy)

        End If ! ! CycleNo == NumCycles .and. ...


        If (  CycleNo == NumCycles .OR. L_new_tdisc ) Then

! DEPENDS ON: ni_update_rho
          Call NI_Update_Rho(                                           &
     &                 rho, rho_n, rho_np1, inc_rho,                    &
     &                 u, v, w, R_u, R_v, R_w,                          &
     &                 q, qcl,qcf,                                      &
     &                 qcf2, qrain, qgraup,                             &
     &                 q_star, qcl_star, qcf_star,                      &
     &                 qcf2_star, qrain_star, qgraup_star,              &
     &                 q_np1, qcl_np1, qcf_np1,                         &
     &                 qcf2_np1, qrain_np1, qgraup_np1,                 &
     &                 mix_v, mix_cl, mix_cf,                           &
     &                 mix_cf2, mix_rain, mix_graup,                    &
     &                 mix_v_star, mix_cl_star, mix_cf_star,            &
     &                 mix_cf2_star, mix_rain_star, mix_graup_star,     &
     &                 mix_v_np1, mix_cl_np1, mix_cf_np1,               &
     &                 mix_cf2_np1, mix_rain_np1, mix_graup_np1,        &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 timestep, CycleNo, NumCycles,                    &
     &                 rows, n_rows, row_length,                        &
     &                 model_levels, wet_levels, model_domain,          &
     &                 first_constant_r_rho_level,                      &
     &                 alpha_1, alpha_2, n_rims_to_do,                  &
     &                 nproc, gc_proc_row_group,                        &
     &                 L_regular, at_extremity, global_row_length,      &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 cos_v_latitude, delta_lambda, delta_phi,         &
     &                 gdlambda_p(lambda_start), dphi_p,                &
     &                 grecip_dlamu(lambda_start), recip_dphiv,         &
     &                 wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v,    &
     &                 r_theta_levels, r_rho_levels, eta_theta_levels,  &
     &                 eta_rho_levels, FV_sec_theta_latitude,           &
     &                 wet_to_dry_n, wet_to_dry_np1,                    &
     &                 L_do_inc_rho, L_new_tdisc,                       &
     &                 L_mix_ratio, L_dry )

        End If !  CycleNo == NumCycles .OR. L_new_tdisc

        If ( CycleNo == NumCycles ) Then

          if(L_mix_ratio .and. .not. L_moist_nonhydro_conserve)then
            deallocate(mix_v)
            deallocate(mix_cl)
            deallocate(mix_cf)
            deallocate(mix_cf2)
            deallocate(mix_rain)
            deallocate(mix_graup)
          endif

          If (NumCycles > 1 .and. L_mix_ratio                           &
     &                      .and. .not. L_moist_nonhydro_conserve) Then
            deallocate( mix_v_star_phys1)
            deallocate( mix_cl_star_phys1)
            deallocate( mix_cf_star_phys1)
            If ( L_mcr_qcf2 ) deallocate( mix_cf2_star_phys1)
            If ( L_mcr_qrain ) deallocate( mix_rain_star_phys1)
            If ( L_mcr_qgraup ) deallocate( mix_graup_star_phys1)
          End If

        If (L_tracer) then
! DEPENDS ON: timer
          If (Ltimer) CALL TIMER('SL_tracer2',3)

! DEPENDS ON: sl_tracer2
          Call SL_tracer2(                                              &
     &                           super_array_size,                      &
     &                  super_tracer_phys1,super_tracer_phys2,          &
     &                 eta_theta_levels,                                &
     &                 r_rho_levels,  r_theta_levels,                   &
     &                 rho_n, RHO,                                      &
     &                 row_length, rows, model_levels,                  &
     &                 delta_lambda, delta_phi,                         &
     &                 glambda_p, phi_p, grecip_dlamp, recip_dphip,     & 
     &                 lambda_p_rm, lambda_p_rp, phi_p_rm, phi_p_rp,    &
     &                 recip_lambda_p_m, recip_lambda_p_0,              &
     &                 recip_lambda_p_p, recip_lambda_p_p2,             &
     &                 recip_phi_p_m, recip_phi_p_0,                    &
     &                 recip_phi_p_p, recip_phi_p_p2,                   &
     &                 base_lambda, base_phi,                           &
     &                 recip_dlam, recip_dphi, max_look,                &
     &                 look_lam, look_phi, halo_lam, halo_phi,          &
     &                 FV_cos_theta_latitude,                           &
     &                 wet_to_dry_n, wet_to_dry_np1,                    &
     &                 mype, nproc, nproc_x, nproc_y,                   &
     &                 halo_i, halo_j, datastart,                       &
     &                 g_i_pe, at_extremity,                            &
     &                 global_row_length,                               &
     &                 global_rows,                                     &
     &                 gc_all_proc_group,                               &
     &                 gc_proc_row_group,                               &
     &                 gc_proc_col_group, offx, offy,                   &
     &                 L_regular, L_sl_halo_reprod,                     &
     &                 high_order_scheme(moist_SL),                     &
     &                 monotone_scheme(moist_SL),                       &
     &                 model_domain, L_high(moist_SL),                  &
     &                 L_mono(moist_SL),                                &
     &                 L_conserve_tracers,                              &
     &                 check_bottom_levels,                             &
     &                 interp_vertical_search_tol,                      &
     &                 first_constant_r_rho_level,                      &
     &                 depart_lambda, depart_phi, depart_r_theta,       &
     &                 co2, L_CO2_interactive,                          &
     &                 Murk, L_Murk_advect,                             &
     &                 soot_new, soot_agd,                              &
     &                 soot_cld, L_soot,                                &
     &                 bmass_new, bmass_agd,                            &
     &                 bmass_cld, L_biomass,                            &
     &                 ocff_new, ocff_agd, ocff_cld, L_ocff,            &
     &                 DUST_DIV1,DUST_DIV2,                             &
     &                 DUST_DIV3,DUST_DIV4,                             &
     &                 DUST_DIV5,DUST_DIV6, L_DUST,                     &
     &                 so2, so4_aitken,                                 &
     &                 so4_accu,                                        &
     &                 so4_diss, nh3, dms,                              &
     &                 L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms,           &
     &                 tracer, tr_levels, tr_vars,                      &
     &                 tracer_ukca, tr_ukca,                            &
     &                 tracer_phys1,tracer_phys2,                       &
     &                 i_start, i_stop, j_start, j_stop,                &
     &                 L_USE_CARIOLLE, OZONE_TRACER,                    &
     &                 ErrorStatus)

!      Cariolle scheme is called to calculate the tracer ozone. All the tracers
!      are calculated at the end of the timestep. The ozone tracer 
!      calculated here will be used in the radiation scheme on the next timestep.

!      Insert if statement here to check that tracer that goes into this call 
!      ozone. Ozone may not be the only tracer.
            If (PrintStatus  >=  PrStatus_Normal .AND.                  &
                first_atmstep_call) then

               WRITE(6,*) 'Atm_Step: L_USE_CARIOLLE = ',L_USE_CARIOLLE
            End If

            If (L_USE_CARIOLLE) then
               
               If (PrintStatus  >=  PrStatus_Normal .AND.                  &
                first_atmstep_call) then

                    WRITE(6,*) 'Atm_Step: Calling Cariolle_o3_psc'
               End If

! DEPENDS ON: cariolle_o3_psc
               CALL cariolle_o3_psc (OZONE_TRACER,                     &
                      O3_PROD_LOSS,   O3_P_L_VMR,                      &
                      O3_VMR,         O3_P_L_TEMP,                     &
                      O3_TEMP,        O3_P_L_COLO3,                    &
                      O3_COLO3,                                        &
                      THETA,                                           &
                      P_THETA_LEVELS,                                  &
                      offx,offy,theta_off_size,                        &
                      rows,row_length,                                 &
                      timestep,                                        &
                      EXNER_THETA_LEVELS,                              &
                      model_levels)      
!    
! Halos updated
! DEPENDS ON: swap_bounds
               call Swap_Bounds(OZONE_TRACER,                           &
                        row_length, rows, model_levels,                 &
                        offx, offy, fld_type_p, .false.)
!DEPENDS ON: fill_external_halos
               CALL FILL_EXTERNAL_HALOS(OZONE_TRACER,row_length, rows,  &
                         model_levels,offx,offy)
            
            else
               If (PrintStatus  >=  PrStatus_Normal .AND.               &
                   first_atmstep_call) then

                   WRITE(6,*) 'Atm_Step: Cariolle scheme not called'
               End If
            End if

! DEPENDS ON: timer
          If (Ltimer) CALL TIMER('SL_tracer2',4)

          If (.not. L_moist_nonhydro_conserve)then
            DEALLOCATE (depart_r_theta)
            DEALLOCATE (depart_lambda)
            DEALLOCATE (depart_phi)
          endif
        end if  ! L_tracer

! ----------------------------------------------------------------------
! Section 5.2 Recalculate sl_thermo to allow rho at n+1
!             to be used in the conservation step
!            Estimated values are returned in the _phys1 variables.
! ----------------------------------------------------------------------
!  At this point, q_star, qcl_star and qcf_star in the boundary layer
!  have been changed by the implicit bl calculation but are not used
!  from now on.
        If (L_moist_nonhydro_conserve) then

! DEPENDS ON: ni_sl_moist
          Call NI_SL_moist(                                             &
     &                 moisture_array_size,                             &
     &                 q, qcl, qcf,                                     &
     &                 cf_bulk,cf_frozen,                               &
     &                 cf_liquid,                                       &
     &                 qcf2, qrain, qgraup,                             &
     &                 q_phys1, qcl_phys1, qcf_phys1,                   &
     &                 cf_phys1, cff_phys1, cfl_phys1,                  &
     &                 qcf2_phys1, qrain_phys1, qgraup_phys1,           &
     &                 mix_v, mix_cl, mix_cf,                           &
     &                 mix_cf2, mix_rain, mix_graup,                    &
     &                 mix_v_phys1, mix_cl_phys1, mix_cf_phys1,         &
     &                 mix_cf2_phys1, mix_rain_phys1, mix_graup_phys1,  &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 L_pc2, L_mix_ratio,                              &
     &                 eta_theta_levels,                                &
     &                 r_rho_levels,  r_theta_levels,                   &
     &                 exner_theta_levels,                              &
     &                 rho_n, rho,                                      &
     &                 depart_lambda, depart_phi, depart_r_theta,       &
     &                 row_length, rows, n_rows,                        &
     &                 model_levels, wet_levels,                        &
     &                 delta_lambda, delta_phi,                         &
     &                 base_lambda, base_phi,                           &
     &                 glambda_p, phi_p, grecip_dlamp, recip_dphip,     &
     &                 lambda_p_rm, lambda_p_rp, phi_p_rm, phi_p_rp,    &
     &                 recip_lambda_p_m, recip_lambda_p_0,              &
     &                 recip_lambda_p_p, recip_lambda_p_p2,             &
     &                 recip_phi_p_m, recip_phi_p_0,                    &
     &                 recip_phi_p_p, recip_phi_p_p2,                   &
     &                 recip_dlam, recip_dphi, max_look,                &
     &                 look_lam, look_phi, halo_lam, halo_phi,          &
     &                 FV_cos_theta_latitude,                           &
     &                 wet_to_dry_n, wet_to_dry_np1,                    &
     &                 L_regular, n_rims_to_do,                         &
     &                 mype, nproc, nproc_x, nproc_y,                   &
     &                 halo_i, halo_j, datastart,                       &
     &                 g_i_pe, at_extremity,                            &
     &                 global_row_length, global_rows,                  &
     &                 gc_proc_row_group,                               &
     &                 gc_proc_col_group, offx, offy,                   &
     &                 L_sl_halo_reprod,                                &
     &                 high_order_scheme(moist_SL),                     &
     &                 monotone_scheme(moist_SL),                       &
     &                 model_domain, L_high(moist_SL),                  &
     &                 L_mono(moist_SL), L_conserv(moist_SL),           &
     &                 check_bottom_levels, interp_vertical_search_tol, &
     &                 first_constant_r_rho_level,                      &
     &                 ErrorStatus )

          DEALLOCATE (depart_lambda)
          DEALLOCATE (depart_phi)
          DEALLOCATE (depart_r_theta)
          if(L_mix_ratio)then
            deallocate(mix_v)
            deallocate(mix_cl)
            deallocate(mix_cf)
            deallocate(mix_cf2)
            deallocate(mix_rain)
            deallocate(mix_graup)
          endif

          If ( NumCycles > 1 .and. L_mix_ratio ) Then
            deallocate( mix_v_star_phys1)
            deallocate( mix_cl_star_phys1)
            deallocate( mix_cf_star_phys1)
            If ( L_mcr_qcf2 ) deallocate( mix_cf2_star_phys1)
            If ( L_mcr_qrain ) deallocate( mix_rain_star_phys1)
            If ( L_mcr_qgraup ) deallocate( mix_graup_star_phys1)
          End If
        end if          ! L_moist_nonhydro_conserve

      End If ! NoCycles == NumCycles

      End If       ! L_Primitive .and. ErrorStatus  ==  0

! ----------------------------------------------------------------------
! Section 6.0 Calculate u, v, w, theta, q, qcl, qcf at new time level.
!             Extrapolate u, v, w, to time level n+1.5 for use in
!             advection step on next timestep.
! ----------------------------------------------------------------------
      If ( ErrorStatus  ==  0) Then

        If ( CycleNo == 1 ) Then

        if(sf(0,30)) then
          if (sf(182,30) .or. sf(183,30) .or. sf(184,30)                &
     &       .or. sf(189,30) .or. sf(190,30) .or. sf(191,30)) then
            ALLOCATE ( inc_q(1-halo_i:row_length+halo_i,                &
     &                       1-halo_j:rows+halo_j,wet_levels) )
            ALLOCATE ( inc_qcl(1-halo_i:row_length+halo_i,              &
     &                         1-halo_j:rows+halo_j,wet_levels) )
            ALLOCATE ( inc_qcf(1-halo_i:row_length+halo_i,              &
     &                         1-halo_j:rows+halo_j,wet_levels) )
            ALLOCATE ( inc_cf(1-halo_i:row_length+halo_i,               &
     &                       1-halo_j:rows+halo_j,wet_levels) )
            ALLOCATE ( inc_cfl(1-halo_i:row_length+halo_i,              &
     &                         1-halo_j:rows+halo_j,wet_levels) )
            ALLOCATE ( inc_cff(1-halo_i:row_length+halo_i,              &
     &                         1-halo_j:rows+halo_j,wet_levels) )

!Store values for increment calculations for moisture variables
            Do k = 1, wet_levels
              Do j = j_start, j_stop
                Do i = i_start, i_stop
                  inc_q(i,j,k)   = q(i,j,k)
                  inc_qcl(i,j,k) = qcl(i,j,k)
                  inc_qcf(i,j,k) = qcf(i,j,k)
                  inc_cf(i,j,k)  = cf_bulk(i,j,k)
                  inc_cfl(i,j,k) = cf_liquid(i,j,k)
                  inc_cff(i,j,k) = cf_frozen(i,j,k)
                End Do
              End Do
            End Do
          else
            ALLOCATE ( inc_q(1,1,1) )
            ALLOCATE ( inc_qcl(1,1,1) )
            ALLOCATE ( inc_qcf(1,1,1) )
            ALLOCATE ( inc_cf(1,1,1) )
            ALLOCATE ( inc_cfl(1,1,1) )
            ALLOCATE ( inc_cff(1,1,1) )
          endif
        endif       ! sf(0,30)

        If ( sf(181,30) ) then
          ALLOCATE ( inc_t( 1-offx:row_length+offx,                     &
     &                      1-offy:rows+offy, model_levels) )
        else
          ALLOCATE ( inc_t( 1,1,1) )
        endif

        L_do_inc_vels = sf(185,30) .or. sf(186,30) .or. sf(187,30)

        if ( L_do_inc_vels ) then
          ALLOCATE ( inc_u( 1-offx:row_length+offx,                     &
     &                      1-offy:rows+offy, model_levels) )
          ALLOCATE ( inc_v( 1-offx:row_length+offx,                     &
     &                      1-offy:n_rows+offy, model_levels) )
          ALLOCATE ( inc_w(row_length, rows, model_levels) )
        else
          ALLOCATE ( inc_u( 1,1,1) )
          ALLOCATE ( inc_v( 1,1,1) )
          ALLOCATE ( inc_w( 1,1,1) )
        endif

        End If ! CycleNo == 1

! DEPENDS ON: update_fields
        call update_fields(                                             &
     &                   NumCycles, CycleNo, L_new_tdisc,               &
     &                   L_mix_ratio, extrp_weight,                     &
     &                   EXNER_RHO_LEVELS,                              &
     &                   EXNER_THETA_LEVELS,                            &
     &                   U, u_np1, V, v_np1,                            &
     &                   W, w_np1,                                      &
     &                   U_ADV, V_ADV, W_ADV,                           &
     &                   THETA,Q,QCL,QCF,                               &
     &                   QCF2, QRAIN, QGRAUP,                           &
     &                   CF_BULK,CF_LIQUID,                             &
     &                   CF_FROZEN,                                     &
     &                   exner_prime, R_u, R_v, R_w,                    &
     &                   theta_star, theta_np1, dtheta_dr_term,         &
     &                   q_star, q_np1, qcl_star, qcl_np1, qcf_star,    &
     &                   qcf_np1, qcf2_star, qcf2_np1, qrain_star,      &
     &                   qrain_np1, qgraup_star, qgraup_np1,            &
     &                   cf_star, cfl_star, cff_star,                   &
     &                   row_length, rows, n_rows, model_levels,        &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   sin_theta_longitude, cos_theta_longitude,      &
     &                   mag_vector_np, dir_vector_np,                  &
     &                   mag_vector_sp, dir_vector_sp,                  &
     &                   global_row_length, gc_proc_row_group,          &
     &                   at_extremity, datastart,model_domain,          &
     &                   alpha_2 , timestep, wet_levels,                &
     &                   i_start, i_stop, j_start, j_stop,              &
     &                   j_begin, j_end,                                &
     &                   delta_lambda, glambda_p(lambda_start),         &
     &                   glambda_u(lambda_start),                       &
     &                   inc_t, inc_u, inc_v, inc_w,                    &
     &                   sf(181,30), L_do_inc_vels, L_pc2,              &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   L_regular )

        If ( NumCycles > 1 .and. CycleNo < NumCycles ) Then
!
! For consistency in the formulation, the _adv and _np1 vars need
! to be updated in the lateral boundary region.
!
! DEPENDS ON: swap_bounds
          call Swap_Bounds(U_ADV,                                       &
     &                   row_length, rows, model_levels,                &
     &                   halo_i, halo_j, fld_type_u, .true.)

! DEPENDS ON: swap_bounds
          call Swap_Bounds(V_ADV,                                       &
     &                   row_length, n_rows, model_levels,              &
     &                   halo_i, halo_j, fld_type_v, .true.)

! DEPENDS ON: swap_bounds
          call Swap_Bounds(W_ADV,                                       &
     &                   row_length, rows, model_levels+1,              &
     &                   halo_i, halo_j, fld_type_p, .false.)

          If ( L_new_tdisc ) Then

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(u_np1,                                     &
     &                   row_length, rows, model_levels,                &
     &                   offx, offy, fld_type_u, .true.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(v_np1,                                     &
     &                   row_length, n_rows, model_levels,              &
     &                   offx, offy, fld_type_v, .true.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(w_np1,                                     &
     &                   row_length, rows, model_levels+1,              &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(                                           &
     &                   theta_np1, row_length, rows, model_levels,     &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
            Call Swap_Bounds(                                           &
     &                   rho_np1, row_length, rows, model_levels,       &
     &                   offx, offy, fld_type_p, .false.)


!
! NOTE: mix_v_star etc will temporarily hold q_star etc. These will be
!       converted to mix_v_np1 at the end of the current iteration.
!
            If ( L_mix_ratio ) Then

              mix_v_star = q_star
              mix_cl_star = qcl_star
              mix_cf_star = qcf_star
              If ( L_mcr_qcf2 ) mix_cf2_star = qcf2_star
              If ( L_mcr_qrain ) mix_rain_star = qrain_star
              If ( L_mcr_qgraup ) mix_graup_star = qgraup_star

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   mix_v_star, row_length, rows, wet_levels,      &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   mix_cl_star, row_length, rows, wet_levels,     &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   mix_cf_star, row_length, rows, wet_levels,     &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qcf2 )                                         &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   mix_cf2_star, row_length, rows, wet_levels,    &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qrain )                                        &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   mix_rain_star, row_length, rows, wet_levels,   &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qgraup )                                       &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   mix_graup_star, row_length, rows, wet_levels,  &
     &                   offx, offy, fld_type_p, .false.)

            Else

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   q_np1, row_length, rows, wet_levels,           &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   qcl_np1, row_length, rows, wet_levels,         &
     &                   offx, offy, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
              Call Swap_Bounds(                                         &
     &                   qcf_np1, row_length, rows, wet_levels,         &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qcf2 )                                         &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   qcf2_np1, row_length, rows, wet_levels,        &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qrain )                                        &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   qrain_np1, row_length, rows, wet_levels,       &
     &                   offx, offy, fld_type_p, .false.)

              If ( L_mcr_qgraup )                                       &
! DEPENDS ON: swap_bounds
     &          Call Swap_Bounds(                                       &
     &                   qgraup_np1, row_length, rows, wet_levels,      &
     &                   offx, offy, fld_type_p, .false.)

            End If ! .NOT. L_mix_ratio

          End If ! If L_new_tdisc

! Now need to fix Uadv and temp n+1 dynamical vars at LB region.
          If ( model_domain == mt_lam                                   &
     &         .or.  model_domain == mt_cyclic_lam ) Then

            If ( L_new_tdisc ) Then
!
! Fill external haloes of theta_np1, q_np1,... following approach
! for theta_star, q_star.
!
               If ( .NOT. L_mix_ratio ) Then

! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(q_np1,row_length,rows,        &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(qcl_np1,row_length,rows,      &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(qcf_np1,row_length,rows,      &
     &                           wet_levels,Offx,Offy)

                 If ( L_mcr_qcf2 ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(qcf2_np1,row_length,rows,   &
     &                           wet_levels,Offx,Offy)
                 End If

                 If ( L_mcr_qrain ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(qrain_np1,row_length,rows,  &
     &                           wet_levels,Offx,Offy)
                 End If

                 If ( L_mcr_qgraup ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(qgraup_np1,row_length,rows, &
     &                           wet_levels,Offx,Offy)
                 End If

               Else
! mix_v_star etc temporarily hold q_star

! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(mix_v_star,row_length,rows,   &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(mix_cl_star,row_length,rows,  &
     &                           wet_levels,Offx,Offy)
! DEPENDS ON: FILL_EXTERNAL_HALOS
                 CALL FILL_EXTERNAL_HALOS(mix_cf_star,row_length,rows,  &
     &                           wet_levels,Offx,Offy)
                 If ( L_mcr_qcf2 ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(mix_cf2_star,row_length,    &
     &                           rows,wet_levels,Offx,Offy)
                 End If
                 If ( L_mcr_qrain ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(mix_rain_star,row_length,   &
     &                           rows,wet_levels,Offx,Offy)
                 End If
                 If ( L_mcr_qgraup ) Then
! DEPENDS ON: FILL_EXTERNAL_HALOS
                   CALL FILL_EXTERNAL_HALOS(mix_graup_star,row_length,  &
     &                           rows,wet_levels,Offx,Offy)
                 End If
               End If ! .NOT. L_mix_ratio

! DEPENDS ON: FILL_EXTERNAL_HALOS
               CALL FILL_EXTERNAL_HALOS(theta_np1,row_length,rows,      &
     &                         model_levels,Offx,Offy)

               If (RIM_STEPSA == 0) Then
                 increment_factor=0.0
               Else
                 increment_factor=1.0/                                  &
     &           (RIM_STEPSA-MOD(Timestep_Number-1,RIM_STEPSA))
               End If

               lbc_size=LENRIMA(fld_type_p,halo_type_extended,          &
     &                         rima_type_norm)
               DO k=1,MODEL_LEVELS
                 DO i=1,lbc_size
                   W_LBC_REAL_TEND(i,k)=W_LBC(i,k)                      &
     &          +  increment_factor*( W_LBC_TEND(i,k) - W_LBC(i,k) )
                 END DO
               END DO

               lbc_size=LENRIMA(fld_type_u,halo_type_extended,          &
     &                          rima_type_norm)
               DO k=1,MODEL_LEVELS
                 DO i=1,lbc_size
                   U_LBC_REAL_TEND(i,k)=U_LBC(i,k)                      &
     &        + increment_factor*( U_LBC_TEND(i,k) - U_LBC(i,k) )      
                 END DO
               END DO

               lbc_size=LENRIMA(fld_type_v,halo_type_extended,          &
     &                         rima_type_norm)
               DO k=1,MODEL_LEVELS
                 DO i=1,lbc_size
                   V_LBC_REAL_TEND(i,k)=V_LBC(i,k)                      &
     &        + increment_factor*( V_LBC_TEND(i,k) - V_LBC(i,k) )
                 END DO
               END DO

               If ( CycleNo == 1 ) Then
                 Allocate (                                             &
     &                    RHO_LBC_REAL_TEND(LENRIMA(fld_type_p,         &
     &                    halo_type_extended,rima_type_norm),           &
     &                    MODEL_LEVELS) )
               End If

               lbc_size=LENRIMA(fld_type_p,halo_type_extended,          &
     &                         rima_type_norm)
               DO k=1,MODEL_LEVELS
                 DO i=1,lbc_size
                 RHO_LBC_REAL_TEND(i,k)=RHO_LBC(i,k)                    &
     &        + increment_factor*(RHO_LBC_TEND(i,k) - RHO_LBC(i,k) )
                 END DO
               END DO

            Else

              If ( CycleNo == 1 ) Allocate ( RHO_LBC_REAL_TEND(1,1) )

            End If ! If L_new_tdisc

! DEPENDS ON: TIMER
            IF (LTIMER) CALL TIMER('UPDATE_LBC_ITERSL',3)

            If ( .NOT. L_mix_ratio ) Then
! DEPENDS ON: UPDATE_LBC_ITERSL
              CALL UPDATE_LBC_ITERSL(                                   &
     &           ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,        &
     &           OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY, L_new_tdisc,     &
     &           L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                       &
     &           L_mcr_qgraup_lbc, L_pc2_lbc,                           &
     &           RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                 &
     &           LENRIMA(1,1,rima_type_norm),                           &
     &           LBC_SIZEA(1,1,1,rima_type_norm),                       &
     &           LBC_STARTA(1,1,1,rima_type_norm),                      &
     &           THETA_LBC,Q_LBC,QCL_LBC,                               &
     &           QCF_LBC,QCF2_LBC,QRAIN_LBC,                            &
     &           QGRAUP_LBC, CF_BULK_LBC,CF_LIQUID_LBC,                 &
     &           CF_FROZEN_LBC, RHO_LBC_REAL_TEND, EXNER_LBC,           &
     &           U_LBC_REAL_TEND,V_LBC_REAL_TEND,W_LBC_REAL_TEND,       &
     &           U_ADV_LBC,V_ADV_LBC,W_ADV_LBC,                         &
     &           theta_np1, q_np1, qcl_np1, qcf_np1,                    &
     &           qcf2_np1, qrain_np1, qgraup_np1,                       &
     &           CF_BULK,CF_LIQUID,CF_FROZEN,                           &
     &           rho_np1,EXNER_RHO_LEVELS,                              &
     &           u_np1,v_np1,w_np1,                                     &
     &           U_ADV,V_ADV,W_ADV)          
            Else
! mix_v_star etc hold q_star etc temporarily
! DEPENDS ON: UPDATE_LBC_ITERSL
              CALL UPDATE_LBC_ITERSL(                                   &
     &           ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,        &
     &           OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY, L_new_tdisc,     &
     &           L_mcr_qcf2_lbc, L_mcr_qrain_lbc,                       &
     &           L_mcr_qgraup_lbc, L_pc2_lbc,                           &
     &           RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                 &
     &           LENRIMA(1,1,rima_type_norm),                           &
     &           LBC_SIZEA(1,1,1,rima_type_norm),                       &
     &           LBC_STARTA(1,1,1,rima_type_norm),                      &
     &           THETA_LBC,Q_LBC,QCL_LBC,                               &
     &           QCF_LBC,QCF2_LBC,QRAIN_LBC,                            &
     &           QGRAUP_LBC, CF_BULK_LBC,CF_LIQUID_LBC,                 &
     &           CF_FROZEN_LBC, RHO_LBC_REAL_TEND, EXNER_LBC,           &
     &           U_LBC_REAL_TEND,V_LBC_REAL_TEND,W_LBC_REAL_TEND,       &
     &           U_ADV_LBC,V_ADV_LBC,W_ADV_LBC,                         &
     &           theta_np1, mix_v_star, mix_cl_star, mix_cf_star,       &
     &           mix_cf2_star, mix_rain_star, mix_graup_star,           &
     &           CF_BULK,CF_LIQUID,CF_FROZEN,                           &
     &           rho_np1,EXNER_RHO_LEVELS,                              &
     &           u_np1,v_np1,w_np1,                                     &
     &           U_ADV,V_ADV,W_ADV)
            End If ! .NOT. L_mix_ratio

! DEPENDS ON: TIMER
            IF (LTIMER) CALL TIMER('UPDATE_LBC_ITERSL',4)

            If ( CycleNo == NumCycles-1 ) Deallocate(RHO_LBC_REAL_TEND)

          End If ! If model_domain == mt_lam ...

! Now convert qX(=mix_v_star,...) to mX and store at mix_v_np1 etc
          If ( L_mix_ratio .and. L_new_tdisc ) Then
            call q_to_mix_halo (row_length, rows, wet_levels,           &
     &                 offx, offy, mix_v_star, mix_cl_star, mix_cf_star,&
     &                 mix_cf2_star, mix_rain_star, mix_graup_star,     &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 mix_v_np1, mix_cl_np1, mix_cf_np1,               &
     &                 mix_cf2_np1, mix_rain_np1, mix_graup_np1         &
     &                 )
          End If

        End If ! NumCycles > 1 .and. CycleNo < NumCycles
!  LBC updating now done at end of time step     
! ---------------------------------------------------------------------
!   Section 6.1  Update lbcs for LAMs
! ---------------------------------------------------------------------

      If ( model_domain == mt_lam ) Then
      
        IF (L_Fixed_lbcs) Then
          L_update_lbcs = .false.
          L_apply_lbcs = .true.
        Else IF ( RIM_STEPSA == 0 ) Then
          L_update_lbcs = .false.
        Else IF ( L_lbc_new ) Then
          L_update_lbcs = .true.
        Else
          L_update_lbcs = .false.
        END IF !  L_Fixed_lbcs

        If ( L_update_lbcs ) Then
        
! DEPENDS ON: BOUNDVAL
              Call BOUNDVAL(LENRIMA(1,1,rima_type_norm),                &
     &                      L_mcr_qcf2_lbc, L_mcr_qrain_lbc,            &
     &                      L_mcr_qgraup_lbc, L_pc2_lbc,                &
     &                      L_murk_lbc, L_int_uvw_lbc,                  &
     &                      U_LBC, U_LBC_TEND,                          &
     &                      V_LBC, V_LBC_TEND,                          &
     &                      W_LBC, W_LBC_TEND,                          &
     &                      RHO_LBC, RHO_LBC_TEND,                      &
     &                      THETA_LBC, THETA_LBC_TEND,                  &
     &                      Q_LBC, Q_LBC_TEND,                          &
     &                      QCL_LBC, QCL_LBC_TEND,                      &
     &                      QCF_LBC, QCF_LBC_TEND,                      &
     &                      QCF2_LBC, QCF2_LBC_TEND,                    &
     &                      QRAIN_LBC, QRAIN_LBC_TEND,                  &
     &                      QGRAUP_LBC, QGRAUP_LBC_TEND,                &
     &                      CF_BULK_LBC, CF_BULK_LBC_TEND,              &
     &                      CF_LIQUID_LBC, CF_LIQUID_LBC_TEND,          &
     &                      CF_FROZEN_LBC, CF_FROZEN_LBC_TEND,          &
     &                      EXNER_LBC, EXNER_LBC_TEND,                  &
     &                      U_ADV_LBC, U_ADV_LBC_TEND,                  &
     &                      V_ADV_LBC, V_ADV_LBC_TEND,                  &
     &                      W_ADV_LBC, W_ADV_LBC_TEND,                  &
     &                      MURK_LBC, MURK_LBC_TEND,                    &
     &                      TRACER_LBC, TRACER_LBC_TEND,                &
     &                      1, 0, ErrorStatus, CMESSAGE)

         End If ! L_update_lbcs

        IF ( L_lbc_new ) Then

        !--------------------------------------------------------------
        !           Update primary fields with LAM LBC data
        !--------------------------------------------------------------

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UPDATE_LAM_LBCS',3)
! DEPENDS ON: update_lam_lbcs
        CALL UPDATE_LAM_LBCS(                                                &
     &    r_rho_levels, r_theta_levels,                                 &
     &    ROW_LENGTH,ROWS,N_ROWS,MODEL_LEVELS,WET_LEVELS,               &
     &    OFFX,OFFY,HALO_I,HALO_J,AT_EXTREMITY,                         &
     &    L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_pc2,                 &
     &    L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc, L_pc2_lbc, &
     &    L_murk, L_murk_lbc,                                           &
     &    L_LBC_balance, L_int_uvw_lbc,                                 &
     &    RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                        &
     &    LENRIMA(1,1,rima_type_norm),                                  &
     &    LBC_SIZEA(1,1,1,rima_type_norm),                              &
     &    LBC_STARTA(1,1,1,rima_type_norm),                             &
     &    THETA_LBC, Q_LBC, QCL_LBC,                                    &
     &    QCF_LBC, QCF2_LBC, QRAIN_LBC,                                 &
     &    QGRAUP_LBC, CF_BULK_LBC, CF_LIQUID_LBC,                       &
     &    CF_FROZEN_LBC, RHO_LBC,EXNER_LBC,                             &
     &    U_LBC, V_LBC, W_LBC,                                          &
     &    U_ADV_LBC, V_ADV_LBC, W_ADV_LBC,                              &
     &    MURK_LBC,                                                     &
     &    THETA, Q,QCL, QCF,                                            &
     &    QCF2, QRAIN, QGRAUP,                                          &
     &    CF_BULK, CF_LIQUID, CF_FROZEN,                                &
     &    RHO, EXNER_RHO_LEVELS,                                        &
     &    U, V, W, U_ADV, V_ADV, W_ADV, MURK,                           &
     &    DELTA_PHI, DELTA_LAMBDA,                                      &
     &    BASE_PHI, BASE_LAMBDA,                                        &
     &    TWO_OMEGA, DATASTART,                                         &
     &    lat_rot_NP,                                                   &
     &    GLOBAL_ROW_LENGTH, GLOBAL_ROWS                                &
     &     )

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UPDATE_LAM_LBCS',4)

        END IF ! L_lbc_new

      ENDIF     !   model_domain  ==  mt_lam

! Clear workspace
        DEALLOCATE (dtheta_dr_term)

      End If !  ErrorStatus == 0

      Enddo ! end iterations for trajectory calc

      If ( ErrorStatus == 0 ) Then

        Deallocate( u_np1 )
        Deallocate( v_np1 )
        Deallocate( w_np1 )
        Deallocate( theta_np1 )
        Deallocate( rho_np1 )
        Deallocate( q_np1 )
        Deallocate( qcl_np1 )
        Deallocate( qcf_np1 )
        Deallocate( qcf2_np1 )
        Deallocate( qrain_np1 )
        Deallocate( qgraup_np1 )
        Deallocate( mix_v_np1 )
        Deallocate( mix_cl_np1 )
        Deallocate( mix_cf_np1 )
        Deallocate( mix_cf2_np1 )
        Deallocate( mix_rain_np1 )
        Deallocate( mix_graup_np1 )
        Deallocate( mix_v_star)
        Deallocate( mix_cl_star)
        Deallocate( mix_cf_star)
        Deallocate( mix_cf2_star)
        Deallocate( mix_rain_star)
        Deallocate( mix_graup_star)
        If ( NumCycles > 1 ) Then
          Deallocate( R_u_phys1 )
          Deallocate( R_v_phys1 )
          Deallocate( thetastar_phys1 )
          Deallocate( qstar_phys1 )
          Deallocate( qclstar_phys1 )
          Deallocate( qcfstar_phys1 )
          If (L_mcr_qcf2)                                               &
     &      Deallocate( qcf2_star_phys1 )
          If (L_mcr_qrain)                                              &
     &      Deallocate( qrain_star_phys1 )
          If (L_mcr_qgraup)                                             &
     &      Deallocate( qgraup_star_phys1 )
          Deallocate( bulk_cld_frac_phys1 )
          Deallocate( bulk_cld_liq_phys1 )
          Deallocate( bulk_cld_fr_phys1 )
          Deallocate( area_cld_frac_phys1 )
          Deallocate( ti_phys1 )
          Deallocate( zh_phys1 )
          Deallocate( z0msea_phys1 )
          Deallocate( cca_phys1 )
          Deallocate( ccb_phys1 )
          Deallocate( cct_phys1 )
          If ( L_ctile ) Then
            Deallocate( T_LAND_CTILE_PHYS1 )
            Deallocate( T_SICE_CTILE_PHYS1 )
          End If
          Deallocate( t_surf_phys1 )
          Deallocate( t_sf_tile_phys1 )
          Deallocate( snow_tile_phys1 )
          Deallocate( dolr_phys1 )

        End If

!     Deallocate extra microphysics variables
        DEALLOCATE (qcf2_star)
        DEALLOCATE (qrain_star)
        DEALLOCATE (qgraup_star)
        DEALLOCATE (exner_prime)

        If (L_moist_nonhydro_conserve)then
          IF(l_mix_ratio)then
            allocate ( mix_v (1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_levels) )
            allocate ( mix_cl(1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_levels) )
            allocate ( mix_cf(1-halo_i:row_length+halo_i,               &
     &                        1-halo_j:rows+halo_j, wet_levels) )
            if(L_mcr_qcf2)then
              allocate ( mix_cf2(1-halo_i:row_length+halo_i,            &
     &                           1-halo_j:rows+halo_j, wet_levels) )
            else
              allocate ( mix_cf2(1,1,1) )
            endif
            if(L_mcr_qrain)then
              allocate ( mix_rain(1-halo_i:row_length+halo_i,           &
     &                            1-halo_j:rows+halo_j, wet_levels) )
            else
              allocate ( mix_rain(1,1,1) )
            endif
            if(L_mcr_qgraup)then
              allocate ( mix_graup(1-halo_i:row_length+halo_i,          &
     &                             1-halo_j:rows+halo_j, wet_levels) )
            else
              allocate ( mix_graup(1,1,1) )
            endif
            Do k = 1, wet_levels
              Do j = 1, rows
                Do i = 1, row_length
                  mix_v(i,j,k) = mix_v_phys1(i,j,k) + mix_v_phys2(i,j,k)
               mix_cl(i,j,k) = mix_cl_phys1(i,j,k) + mix_cl_phys2(i,j,k)
               mix_cf(i,j,k) = mix_cf_phys1(i,j,k) + mix_cf_phys2(i,j,k)
                End Do
              End Do
            End Do
            if(L_mcr_qcf2)then
              Do k = 1, wet_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    mix_cf2(i,j,k) = mix_cf2_phys1(i,j,k) +             &
     &                               mix_cf2_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            if(L_mcr_qrain)then
              Do k = 1, wet_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    mix_rain(i,j,k) = mix_rain_phys1(i,j,k) +           &
     &                                mix_rain_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            if(L_mcr_qgraup)then
              Do k = 1, wet_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    mix_graup(i,j,k) = mix_graup_phys1(i,j,k) +         &
     &                                 mix_graup_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
! DEPENDS ON: mix_to_q
            call mix_to_q(                                              &
     &                   row_length, rows, wet_levels,                  &
     &                   halo_i, halo_j,                                &
     &                   mix_v, mix_cl, mix_cf,                         &
     &                   mix_cf2, mix_rain, mix_graup,                  &
!    &                   .false.,.false.,.false.,
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   Q, QCL, QCF,                                   &
     &                   QCF2, QRAIN, QGRAUP  )

            if(L_pc2)then
              Do k = 1, wet_levels
                Do j = j_start, j_stop
                  Do i = i_start, i_stop
                    cf_bulk(i,j,k)   = cf_phys1(i,j,k) +cf_phys2(i,j,k)
                    cf_liquid(i,j,k) = cfl_phys1(i,j,k)+cfl_phys2(i,j,k)
                    cf_frozen(i,j,k) = cff_phys1(i,j,k)+cff_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif

            DEALLOCATE (mix_v_phys1)
            DEALLOCATE (mix_cl_phys1)
            DEALLOCATE (mix_cf_phys1)
            DEALLOCATE (mix_v_phys2)
            DEALLOCATE (mix_cl_phys2)
            DEALLOCATE (mix_cf_phys2)
            DEALLOCATE (mix_v)
            DEALLOCATE (mix_cl)
            DEALLOCATE (mix_cf)
            DEALLOCATE (cf_phys1)
            DEALLOCATE (cfl_phys1)
            DEALLOCATE (cff_phys1)
            DEALLOCATE (cf_phys2)
            DEALLOCATE (cfl_phys2)
            DEALLOCATE (cff_phys2)
            DEALLOCATE (mix_cf2_phys1)
            DEALLOCATE (mix_rain_phys1)
            DEALLOCATE (mix_graup_phys1)
            DEALLOCATE (mix_cf2_phys2)
            DEALLOCATE (mix_rain_phys2)
            DEALLOCATE (mix_graup_phys2)
            DEALLOCATE (mix_cf2)
            DEALLOCATE (mix_rain)
            DEALLOCATE (mix_graup)
          else
            Do k = 1, wet_levels
              Do j = j_start, j_stop
                Do i = i_start, i_stop
                  q(i,j,k) = q_phys1(i,j,k) + q_phys2(i,j,k)
                End Do
              End Do
            End Do

            Do k = 1, wet_levels
              Do j = j_start, j_stop
!CDIR NODEP
                Do i = i_start, i_stop
                 qcl(i,j,k) = qcl_phys1(i,j,k) + qcl_phys2(i,j,k)
                 qcf(i,j,k) = qcf_phys1(i,j,k) + qcf_phys2(i,j,k)
                End Do
              End Do
            End Do
            if(L_mcr_qcf2)then
              Do k = 1, wet_levels
                Do j = j_start, j_stop
                  Do i = i_start, i_stop
                    qcf2(i,j,k) = qcf2_phys1(i,j,k) + qcf2_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            if(L_mcr_qrain)then
              Do k = 1, wet_levels
                Do j = j_start, j_stop
                  Do i = i_start, i_stop
                  qrain(i,j,k) = qrain_phys1(i,j,k) + qrain_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            if(L_mcr_qgraup)then
              Do k = 1, wet_levels
                Do j = j_start, j_stop
                  Do i = i_start, i_stop
                qgraup(i,j,k) = qgraup_phys1(i,j,k) + qgraup_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            if(L_pc2)then
              Do k = 1, wet_levels
                Do j = j_start, j_stop
                  Do i = i_start, i_stop
                    cf_bulk(i,j,k)   = cf_phys1(i,j,k) +cf_phys2(i,j,k)
                    cf_liquid(i,j,k) = cfl_phys1(i,j,k)+cfl_phys2(i,j,k)
                    cf_frozen(i,j,k) = cff_phys1(i,j,k)+cff_phys2(i,j,k)
                  End Do
                End Do
              End Do
            endif
            DEALLOCATE (cf_phys1)
            DEALLOCATE (cfl_phys1)
            DEALLOCATE (cff_phys1)
            DEALLOCATE (cf_phys2)
            DEALLOCATE (cfl_phys2)
            DEALLOCATE (cff_phys2)
            DEALLOCATE (q_phys1)
            DEALLOCATE (qcl_phys1)
            DEALLOCATE (qcf_phys1)
            DEALLOCATE (q_phys2)
            DEALLOCATE (qcl_phys2)
            DEALLOCATE (qcf_phys2)
            DEALLOCATE (qcf2_phys1)
            DEALLOCATE (qcf2_phys2)
            DEALLOCATE (qrain_phys1)
            DEALLOCATE (qrain_phys2)
            DEALLOCATE (qgraup_phys1)
            DEALLOCATE (qgraup_phys2)
          endif    !L_mix_ratio

! Use q_pos to remove any negative qcf values
! qlimit being 0.0 try to do locally

! DEPENDS ON: q_pos_ctl
          call Q_Pos_Ctl(                                               &
     &                   QCF, row_length, rows, wet_levels,             &
     &                   global_row_length, global_rows,                &
     &                   mype, nproc, halo_i, halo_j,                   &
     &                   gc_all_proc_group,                             &
     &                   model_domain,                                  &
     &                   halo_type_extended, l_q_pos_local, 0.0         &
     &                   )

        Else If ( NumCycles >1 .and. L_pc2 ) Then
! If not L_moist_nonhydro_conserve but cycling and L_pc2
! is used deallocate cf_phys1, cfl_phys1, cff_phys1
          DEALLOCATE (cf_phys1)
          DEALLOCATE (cfl_phys1)
          DEALLOCATE (cff_phys1)
        End If ! L_moist_nonhydro_conserve



      End If        !  ErrorStatus  ==  0

! ----------------------------------------------------------------------
! section 6.2  Check for q below qlimit and reset
! ----------------------------------------------------------------------

      If(L_qpos)then

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('Q_Pos_Ctl',3)
! DEPENDS ON: q_pos_ctl
        call Q_Pos_Ctl(                                                 &
     &                   Q, row_length, rows, wet_levels,               &
     &                   global_row_length, global_rows,                &
     &                   mype, nproc, halo_i, halo_j,                   &
     &                   gc_all_proc_group,                             &
     &                   model_domain,                                  &
     &                   halo_type_extended, l_q_pos_local, qlimit      &
     &                   )
!
        If (L_pc2) then
! For the PC2 cloud scheme we also need to limit the liquid

! DEPENDS ON: q_pos_ctl
          call Q_Pos_Ctl(                                               &
     &                   QCL, row_length, rows, wet_levels,             &
     &                   global_row_length, global_rows,                &
     &                   mype, nproc, halo_i, halo_j,                   &
     &                   gc_all_proc_group,                             &
     &                   model_domain,                                  &
     &                   halo_type_extended, l_q_pos_local, 0.0         &
     &                   )
!
       End If       ! L_pc2
!
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('Q_Pos_Ctl',4)

      End If       ! L_qpos

      if(sf(0,30)) then
        if(sf(182,30).or.sf(183,30).or.sf(184,30).or.sf(189,30).or.     &
     &     sf(190,30).or.sf(191,30) ) then
! increment calculations for moisture variables
          Do k = 1, wet_levels
            Do j = j_start, j_stop
              Do i = i_start, i_stop
               inc_q(i,j,k)   = q(i,j,k)   - inc_q(i,j,k)
               inc_qcl(i,j,k) = qcl(i,j,k) - inc_qcl(i,j,k)
               inc_qcf(i,j,k) = qcf(i,j,k) - inc_qcf(i,j,k)
               inc_cf(i,j,k)  = cf_bulk(i,j,k)   - inc_cf(i,j,k)
               inc_cfl(i,j,k) = cf_liquid(i,j,k) - inc_cfl(i,j,k)
               inc_cff(i,j,k) = cf_frozen(i,j,k) - inc_cff(i,j,k)
              End Do
            End Do
          End Do
        endif
      endif      !    sf(0,30)

! ----------------------------------------------------------------------
! Section 7.0 Mean all polar variables on a level to remove deviation
!             of values due to rounding error.
!             Only called every polar_reset_timesteps
! ----------------------------------------------------------------------
      If (ErrorStatus  ==  0 .and. model_domain  ==  mt_global          &
     &    .and. L_polar_reset) Then

! need if test on correct timestep
        If (Timestep_Number  ==  1 .or. (                               &
     &      (Timestep_number-1)/polar_reset_timesteps*                  &
     &      polar_reset_timesteps  ==  Timestep_number-1 ) )            &
     &    Then

! DEPENDS ON: polar_reset_mean
          Call Polar_Reset_Mean(                                        &
     &                      EXNER_RHO_LEVELS,RHO,                       &
     &                      THETA,W,                                    &
     &                      Q, QCL, QCF,                                &
     &                      CF_BULK, CF_LIQUID,                         &
     &                      CF_FROZEN,                                  &
     &                      row_length, rows, model_levels,             &
     &                      wet_levels, global_row_length,              &
     &                      offx, offy, halo_i, halo_j,                 &
     &                      nproc, nproc_y, gc_proc_row_group,          &
     &                      at_extremity)
        End If
      End If

! ------------------------------------------------------------------
! Section 17  Aerosol Modelling - includes Sulphur Cycle and Soot
!
! ------------------------------------------------------------------
!
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AEROSOL MODELLING',5)
!
      IF (L_SULPC_SO2.OR.L_SOOT.OR.L_BIOMASS.OR.L_OCFF) THEN
!
! Allocate diagnostic space for STASH
        ALLOCATE (STASHwork17(STASH_maxlen(17,A_im)))
!
! Don't call Swap_bounds for fields used in Aero_Ctl
!
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('AEROSOL MODELLING',3)
!
        ALLOCATE(O3_MMR(row_length,rows,model_levels))
        ALLOCATE(H2O2_MMR(row_length,rows,model_levels))
        ALLOCATE(OH_conc(row_length,rows,model_levels))
        ALLOCATE(HO2_conc(row_length,rows,model_levels)) 
!

! DEPENDS ON: get_sulpc_oxidants
        CALL GET_SULPC_OXIDANTS(                                        &
     & L_SULPC_ONLINE_OXIDANTS, L_UKCA, L_UKCA_TROP, L_UKCA_TROPISOP,   &
     & L_UKCA_STRATTROP, first_atmstep_call,                            &
#include "arg_atm_fields.h"
#include "argd1.h"
     & O3_MMR, H2O2_MMR, OH_conc, HO2_conc)
!
! DEPENDS ON: aero_ctl
        CALL AERO_CTL(                                                  &
! Parallel variables
     &  halo_i, halo_j, offx, offy, global_row_length, global_rows      &
     &, gc_proc_row_group, gc_proc_col_group                            &
     &, at_extremity, nproc, nproc_x, nproc_y                           &
     &, neighbour, g_rows, g_row_length, g_datastart, mype              &
! model dimensions
     &, row_length, rows, n_rows, land_field                            &
     &, model_levels, wet_levels, bl_levels, n_cca_lev                  &
     &, theta_field_size                                                &
     &, salt_dim1, salt_dim2, salt_dim3                                 &
     &, aero_dim1, aero_dim2, aero_dim3                                 &
! Model switches
     &, model_domain, LCAL360, L_SEC_VAR, L_EqT, Ltimer                 &
! Model parameters
     &, Ntot_land, Ntot_sea                                             &
! Physical constants
     &, Lc, Lf, Cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse, earth_radius, Pi                                   &
! Co-ordinate information
     &, r_rho_levels, r_theta_levels                                    &
     &, eta_theta_levels, eta_rho_levels                                &
     &, delta_lambda, delta_phi                                         &
     &, lat_rot_NP, long_rot_NP                                         &
! Time stepping information
     &, timestep                                                        &
     &, I_year, I_day_number, I_hour, I_minute                          &
     &, I_second, timestep_number                                       &
     &, PREVIOUS_TIME                                                   &
     &, CALL_CHEM_FREQ                                                  &
! Trig arrays
     &, sin_theta_longitude, cos_theta_longitude                        &
     &, FV_cos_theta_latitude                                           &
! Grid-dependent arrays
     &,     f3_at_u, true_longitude, true_latitude                      &
!
! Data fields IN
     &, U, V, TSTAR, TSTAR_SEA                                          &
     &, THETA, Q, QCL, QCF                                              &
     &, RHO, LAND, FRAC_LAND, P_THETA_LEVELS                            &
     &, EXNER_RHO_LEVELS, EXNER_THETA_LEVELS                            &
     &, ICE_FRACTION, SNODEP                                            &
     &, CF_BULK                                                         &
     &, OH_conc, H2O2_MMR, HO2_conc, O3_MMR                             &
     &, SO2_EM, SO2_HILEM, SO2_NATEM                                    &
     &, DMS_EM, DMS_CONC, NH3_EM                                        &
     &, DMS_OFLUX                                                       &
     &, SOOT_EM, SOOT_HILEM, BMASS_EM, BMASS_HILEM, OCFF_EM, OCFF_HILEM &
     &, SO2_HIGH_LEVEL, SOOT_HIGH_LEVEL, BMASS_HIGH_LEVEL_1             &
     &, BMASS_HIGH_LEVEL_2, OCFF_HIGH_LEVEL, land_index                 &
! Logicals IN
     &, L_SULPC_SO2, L_SULPC_DMS, L_SULPC_OZONE                         &
     &, L_SULPC_SO2_O3_NONBUFFERED, L_SULPC_NH3                         &
     &, L_use_sulphate_sulpc, L_use_seasalt_sulpc, L_SOOT               &
     &, l_use_soot_sulpc, l_biomass, l_use_bmass_sulpc                  &
     &, l_ocff, l_use_ocff_sulpc                                        &
     &, L_SO2_SURFEM, L_SO2_HILEM, L_SO2_NATEM, L_DMS_EM                &
     &, L_DMS_em_inter, L_DMS_Liss_Merlivat                             &
     &, L_DMS_Wanninkhof, L_DMS_Nightingale                             &
     &, L_DMS_Ointer                                                    &
     &, L_NH3_EM, L_CTILE                                               &
     &, L_SOOT_SUREM, L_SOOT_HILEM, L_BMASS_SUREM, L_BMASS_HILEM        &
     &, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_BIOGENIC                      &
     &, L_USE_SEASALT_DIRECT, L_USE_SEASALT_INDIRECT                    &
     &, L_USE_SEASALT_AUTOCONV, L_DUST                                  &
!
! Data fields IN/OUT
     &, SO2, DMS                                                        &
     &, SO4_AITKEN, SO4_ACCU, SO4_DISS                                  &
     &, H2O2, NH3                                                       &
     &, SOOT_NEW, SOOT_AGD, SOOT_CLD                                    &
     &, BMASS_NEW, BMASS_AGD, BMASS_CLD                                 &
     &, OCFF_NEW, OCFF_AGD, OCFF_CLD                                    &
     &, biogenic                                                        &
!
! Data fields IN
     &, DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5           &
!
! Data fields OUT
! Diagnostic info
     &,                                                                 &
#include "argsts.h"
     &  STASHwork17                                                     &
! Error info
     &, ErrorStatus                                                     &
     & )
!
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('AEROSOL MODELLING',4)
!
! Don't call Swap_bounds for updated fields
!
! Diagnostics STASHed for Aerosol section 17
!
! DEPENDS ON: timer
        If (Ltimer) CALL TIMER('STASH',3)
!
! DEPENDS ON: stash
        CALL STASH(a_sm,a_im,17,STASHwork17,                            &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &    ErrorStatus,Cmessage)
!
! DEPENDS ON: timer
        If (Ltimer) CALL TIMER('STASH',4)
!
        DEALLOCATE (STASHwork17)
!
        IF(ALLOCATED(O3_MMR)) DEALLOCATE(O3_MMR)
        IF(ALLOCATED(H2O2_MMR)) DEALLOCATE(H2O2_MMR)
        IF(ALLOCATED(OH_conc)) DEALLOCATE(OH_conc)
        IF(ALLOCATED(HO2_conc)) DEALLOCATE(HO2_conc)

      END IF         ! END L_SULPC_SO2.OR.L_SOOT TEST
!
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('AEROSOL MODELLING',6)
!
! ----------------------------------------------------------------------
! Section 8.0 Update pressure to value at new-time level.
!             Calculate surface pressure, and exner and pressure at
!             theta levels.
!             If old lbcs then need to update exner on boundaries first.
! ----------------------------------------------------------------------
      If (ErrorStatus  ==  0 .and. model_domain  ==  mt_lam) Then

       If ( .NOT. L_lbc_new ) Then

! Calculate EXNER_LBCs at next time level

! Change halo size for exner_lbc
        lbc_size=LENRIMA(fld_type_p,halo_type_extended,rima_type_norm)

        If (RIM_STEPSA  ==  0) Then

          ! No LBC updating
          Do k=1,MODEL_LEVELS+1
            Do i=1,lbc_size
              EXNER_LBC_REAL_TEND(i,k)=EXNER_LBC(i,k)  
            End Do !i
          End Do !k

        Else If (MOD(Timestep_Number,RIM_STEPSA)  ==  0) Then

          ! End of current LBC period
          DO k=1,MODEL_LEVELS+1
            DO i=1,lbc_size
              EXNER_LBC_REAL_TEND(i,k)= EXNER_LBC_TEND(i,k)
            ENDDO !i
          ENDDO !k

        ELSE ! Just a normal timestep during a LBC period

          increment_factor=1.0/                                         &
     &      (RIM_STEPSA-MOD(Timestep_Number-1,RIM_STEPSA))

          DO k=1,MODEL_LEVELS+1
            DO i=1,lbc_size
              EXNER_LBC_REAL_TEND(i,k) = EXNER_LBC(i,k) +               &
     &                                              increment_factor *  &
     &                            (EXNER_LBC_TEND(i,k) - EXNER_LBC(i,k))
            ENDDO !i
          ENDDO !k

        ENDIF ! End of current LBC period?

!     n_rims_to_do set at begining of subroutine. Only wts=1 rims
!     updated here since derived p fields are only being used for
!     diagnostics. They get properly weighted at start of next timestep
        L_do_halos=.FALSE.
        L_do_boundaries=.TRUE.

! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
     &    ROW_LENGTH,ROWS,Offx,Offy,                                    &
     &    MODEL_LEVELS+1,fld_type_p,EXNER_RHO_LEVELS,                   &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm),   &
     &    halo_i, halo_j,                                               &
     &    EXNER_LBC_REAL_TEND,                                          &
     &    RIMWIDTHA(rima_type_norm),n_rims_to_do,                       &
     &    RIMWEIGHTSA,AT_EXTREMITY,                                     &
     &    L_do_boundaries,L_do_halos)

       End If ! .NOT. L_lbc_new

      Endif      ! ErrorStatus  ==  0 .and. model_domain  ==  mt_lam

      If (ErrorStatus  ==  0) Then

! Check for negative pressure if requested to
        If (Instability_diagnostics  >   0) Then
          Do k = 1, model_levels+1
            Do j = 1, rows
              Do i = 1, row_length
                ij = i+offx + (j+offy-1) * (row_length+2*offx)
                If (exner_rho_levels(((k-1)*theta_off_size)+ij) <  0.) Then
                  ErrorStatus = 123
                End If
              End Do
            End Do
          End Do
! ErrorStatus 123 message now generated outside of the loops.
          if ( ErrorStatus == 123 ) THEN
! DEPENDS ON: ereport
            Call Ereport("ATM_STEP", ErrorStatus,                       &
     &           "Negative pressure value about to be created" )
          endif
        End If       !    Instability_diagnostics  >   0

      End If      ! ErrorStatus  ==  0

      If ( ErrorStatus  ==  0) Then
! DEPENDS ON: consistent_pressure
        call     Consistent_Pressure (                                  &
     &           exner_rho_levels,                                      &
     &           offx,offy,halo_i,halo_J,                               &
     &           row_length,rows,model_levels,                          &
     &           kappa, g, r_theta_levels, r_rho_levels, rho, p_zero,   &
     &           p, pstar, p_theta_levels,exner_theta_levels)
!
        IF (L_physics) then
!
! Are we using the PC2 cloud scheme?
!
        If (L_pc2) then
!
! ----------------------------------------------------------------------
! PC2: Calculate condensation due to changes in temperature resulting
!    from adiabatic changes in pressure (mainly from vertical advection)
!    Also call checking routine from within this subroutine.
! ----------------------------------------------------------------------
!
          allocate(t_inc_pres  (row_length,rows,wet_levels))
          allocate(q_inc_pres  (row_length,rows,wet_levels))
          allocate(qcl_inc_pres(row_length,rows,wet_levels))
          allocate(qcf_inc_pres(row_length,rows,wet_levels))
          allocate(cf_inc_pres (row_length,rows,wet_levels))
          allocate(cfl_inc_pres(row_length,rows,wet_levels))
          allocate(cff_inc_pres(row_length,rows,wet_levels))
          allocate(t_dini  (row_length,rows,wet_levels))
          allocate(q_dini  (row_length,rows,wet_levels))
          allocate(qcl_dini(row_length,rows,wet_levels))
          allocate(qcf_dini(row_length,rows,wet_levels))
          allocate(cf_dini (row_length,rows,wet_levels))
          allocate(cfl_dini(row_length,rows,wet_levels))
          allocate(cff_dini(row_length,rows,wet_levels))

! DEPENDS ON: pc2_pressure_forcing
          Call pc2_pressure_forcing(                                    &
     &              halo_i, halo_j, offx, offy,                         &
     &              P,PSTAR,                                            &
     &              P_THETA_LEVELS, wet_levels,                         &
     &              row_length, rows, rhc_row_length, rhc_rows,         &
     &              timestep, rhcpt, THETA, CF_BULK,                    &
     &              CF_LIQUID, CF_FROZEN,                               &
     &              Q, QCL,QCF,                                         &
     &              exner_star, EXNER_THETA_LEVELS,                     &
     &              CCB, CUMULUS, rhts, tlts, qtts, ptts, CF_AREA,      &
     &              t_inc_pres, q_inc_pres, qcl_inc_pres, qcf_inc_pres, &
     &              cf_inc_pres, cfl_inc_pres,                          &
     &   cff_inc_pres, t_dini, q_dini, qcl_dini, qcf_dini,              &
     &   cf_dini, cfl_dini, cff_dini, l_mr_pc2, L_ACF_Cusack )

          deallocate(rhts)
          deallocate(tlts)
          deallocate(qtts)
          deallocate(ptts)

        End if  ! L_pc2

          If ( (.not. L_pc2) .or. l_pc2_reset ) Then
! ----------------------------------------------------------------------
! Call cloud scheme to make cloud consistent with moisture fields
! ----------------------------------------------------------------------

! DEPENDS ON: qt_bal_cld
            call qt_bal_cld(                                            &
     &        PSTAR,P_THETA_LEVELS,P,                                   &
     &        THETA,EXNER_THETA_LEVELS,                                 &
     &        Q,QCL,QCF,QCF2,                                           &
     &        rhcpt, rhc_row_length, rhc_rows, bl_levels,               &
     &        cloud_fraction_method,overlap_ice_liquid,                 &
     &        ice_fraction_method,ctt_weight,t_weight,                  &
     &        qsat_fixed,sub_cld,                                       &
     &        row_length,rows,model_levels,wet_levels,                  &
     &        offx,offy,halo_i,halo_j,                                  &
     &        delta_lambda, delta_phi,                                  &
     &        r_theta_levels, FV_cos_theta_latitude,                    &
     &        lc, cp, L_cld_area, L_ACF_Cusack, L_ACF_Brooks,           &
     &        L_eacf, L_mcr_qcf2,                                       &
     &        l_mr_qtbalcld,                                            &
     &        ntml, cumulus, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN,    &
     &        mype)

          Endif   ! L_pc2 and L_pc2_reset
          If (L_run_with_physics2) DEALLOCATE ( RHCPT )
        Endif   ! L_physics

        if (sf(181,30) ) then
! Store increment to T
          Do k = 1, model_levels
            Do j = j_start, j_stop
              Do i = i_start, i_stop
                inc_t(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k) &
     &            - inc_t(i,j,k)
              End Do
            End Do
          End Do
        endif     ! stashflag

      End If      !  ErrorStatus  ==  0

! DEPENDS ON: timer
      If (Ltimer) CALL TIMER('AS5-8 Updates',6)

! ----------------------------------------------------------------------
! Section 9.0a Calculate Mass and energy of atmosphere if required
!              using end of timestep values
!             Only done every energy correction period.
! ----------------------------------------------------------------------
      IF (L_EMCORR) THEN

! Set energy correction for use in section 30
        energy_corr_now=A_REALHD(rh_energy_corr)

        IF (LENERGY) THEN
! DEPENDS ON: timer
          IF(LTIMER) CALL TIMER('AS9 Energy mass   ',5)
! zero fields to be calculated

          tot_energy_final = 0.0
          tot_dry_mass_final = 0.0
          tot_moist_final = 0.0

! DEPENDS ON: eng_mass_diag
          Call eng_mass_diag (                                          &
! Parallel variables
     &                      halo_i, halo_j, offx, offy                  &
     &,                     global_row_length, gc_proc_row_group        &
     &,                     gc_proc_col_group                           &
     &,                     at_extremity, nproc, nproc_x, nproc_y       &
     &,                     neighbour                                   &
     &,                     mype                                        &
! model info
     &,                     row_length, rows, n_rows                    &
     &,                     model_domain                                &
     &,                     model_levels,wet_levels                     &
     &,                     r_theta_levels,r_rho_levels                 &
     &,                     delta_lambda,delta_phi                      &
! trig
     &,                     FV_cos_theta_latitude,cos_v_latitude        &
     &,                     cos_theta_longitude,sin_theta_longitude     &
! data arrays in only
     &,                     THETA , U, V                                &
     &,                     W, RHO , Q                                  &
     &,                     QCL, QCF                                    &
     &,                     EXNER_THETA_LEVELS                          &
!     &,                     d1(jexner_rho_levels(1))
!     &,                     d1(jp(1)), d1(jp_theta_levels(1))
!     &,                     d1(jpstar)
! sum of moist fluxes
     &,                     NET_MFLUX                                   &
     &,                     A_REALHD(rh_tot_mass_init)                  &
     &,                     A_REALHD(rh_tot_m_init)                     &
! logical to indicate mass and moist correction required
     &,                     Lmass_corr,Lqt_corr,Lemq_print              &
! energy correction timestep info
     &,                     a_energysteps,timestep                      &
! IN/OUT  results from calculations
     &,                     tot_energy_final,tot_dry_mass_final         &
     &,                     tot_moist_final)

! DEPENDS ON: timer
          If (Ltimer) Call timer ('cal_eng_mass_corr',3)

! DEPENDS ON: cal_eng_mass_corr
          Call cal_eng_mass_corr (                                      &
     &                      row_length, rows, mype                      &
     &,                     a_energysteps,timestep                      &
     &,                     delta_lambda,delta_phi                      &
     &,                     halo_i, halo_j                              &
     &,                     NET_FLUX                                    &
     &,                     A_REALHD(rh_tot_mass_init)                  &
     &,                     A_REALHD(rh_tot_energy_init)                &
     &,                     A_REALHD(rh_energy_corr)                    &
     &,                     tot_energy_final )

! DEPENDS ON: timer
          If (Ltimer) Call timer ('cal_eng_mass_corr',4)

! Swap initial energy and final energy.

          A_REALHD(rh_tot_energy_init) = tot_energy_final

! Swap initial moisture and final moisture.

          A_REALHD(rh_tot_m_init) = tot_moist_final

! DEPENDS ON: timer
          IF(LTIMER) CALL TIMER('AS9 Energy mass   ',6)

        ENDIF   ! LENERGY

      ELSE
! Set energy correction for use in section 30
        energy_corr_now=0.0
      ENDIF     ! L_EMCORR

! ----------------------------------------------------------------------
! Section ?: If required, call temporal filtering control routine.
! ----------------------------------------------------------------------
      Sec = STEPIM(A_IM) * INT(SECS_PER_STEPIM(A_IM))

      L_CallTFilt = ( L_IAU .AND. (Sec >= IAU_StartMin * 60)            &
     &                      .AND. (Sec <= IAU_EndMin   * 60) )          &
     &              .OR.                                                &
     &              ( L_TDF .AND. (Sec >= TDF_StartMin * 60)            &
     &                      .AND. (Sec <= TDF_EndMin   * 60) )

      IF (L_CallTFilt) THEN

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('TFILTCTL',3)

! DEPENDS ON: tfilt_cntl
        CALL TFilt_cntl (                                               &
      U, V, W, U_adv, V_adv, W_adv,                                     &
      Theta, Exner_rho_levels, Rho,                                     &
      Q, Qcl, Qcf, Murk, ozone_tracer,                                  &
      Deep_soil_temp,                                                   &
      P, P_theta_levels, Exner_theta_levels,                            &
      snodep,                                                           &
      cf_area, cf_bulk, cf_liquid, cf_frozen,                           &
      Pstar, Tstar, Tstar_tile,                                         &
#include "argduma.h"
#include "argcona.h"
#include "arglndm.h"
#include "argppx.h"
     &                    l_mr_tfiltctl,                                &
                                         ! In. Use mixing ratios
     &                    IAU_lookup,                                   &
                                      ! inout
     &                    D1_IAU_k4,                                    &
                                      ! inout
     &                    D1_TFilt )  ! inout

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('TFILTCTL',4)

      END IF      !  L_CallTFilt
!
! Are we using the PC2 cloud scheme to determine area cloud fraction?
!
        If (L_pc2 .and. .not. L_pc2_reset) then
          If (.not. L_cld_area) then
!
! ----------------------------------------------------------------------
! PC2: Set area cloud fraction to the bulk cloud fraction. Use the
!    D1 arrays directly
! ----------------------------------------------------------------------
!
          Do k = 1, wet_levels
            Do j = j_start, j_stop
              Do i = i_start, i_stop
                CF_AREA(i,j,k) = CF_BULK(i,j,k)
              End do
            End do
          End do

          Else If (L_cld_area) then

            If (L_ACF_Brooks) then
              Allocate ( cf_bulk_nohalo(row_length,rows,wet_levels) )
              Allocate ( cf_liquid_nohalo(row_length,rows,wet_levels) )
              Allocate ( cf_frozen_nohalo(row_length,rows,wet_levels) )

! Place bulk, liquid and frozen cloud fractions in halo-free arrays
! Use indexing over the full row and row_length (including any LAM
! boundary rim) since the call to ls_acf_brooks uses this indexing.
              Do k = 1, wet_levels
                Do j = 1, rows
                  Do i = 1, row_length
                    cf_bulk_nohalo(i,j,k)   = cf_bulk(i,j,k)
                    cf_liquid_nohalo(i,j,k) = cf_liquid(i,j,k)
                    cf_frozen_nohalo(i,j,k) = cf_frozen(i,j,k)
                  End Do
                End Do
              End Do

! DEPENDS ON: ls_acf_brooks
              Call LS_ACF_Brooks (                                      &
     &             halo_i, halo_j, Offx, Offy                           &
     &,            row_length, rows, model_levels, wet_levels           &
     &,            r_theta_levels, delta_lambda, delta_phi              &
     &,            FV_cos_theta_latitude                                &
     &,            cf_bulk_nohalo, cf_liquid_nohalo                     &
     &,            cf_frozen_nohalo, cumulus                            &
     &,            CF_AREA )

              Deallocate ( cf_bulk_nohalo )
              Deallocate ( cf_liquid_nohalo )
              Deallocate ( cf_frozen_nohalo )

            End If ! L_ACF_Brooks

          End If ! L_cld_area
!
        End if  ! L_pc2 .and. .not. L_pc2_reset

! ----------------------------------------------------------------------
! Section 9.0 Diagnostics at end of timestep
! ----------------------------------------------------------------------

! DEPENDS ON: timer
      IF(LTIMER) CALL TIMER('AS9 End TStep Diags',5)

! section 15: 'dynamics' based quantities
      IF(      SF(0,15)                                                 &
                                ! Diagnostics required for this section
     &   .AND. ErrorStatus == 0) THEN

! DEPENDS ON: st_diag1
        CALL St_diag1(STASH_maxlen(15,A_im),                            &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "arg_atm_fields.h"
#include "argcona.h"
#include "arglndm.h"
#include "argppx.h"
     & ErrorStatus,CMessage)

      ENDIF ! Diagnostics required for this section

! section 16: 'physics' based quantities
      IF(      SF(0,16)                                                 &
                                ! Diagnostics required for this section
     &   .AND. ErrorStatus == 0) THEN

! DEPENDS ON: st_diag2
        CALL St_diag2(STASH_maxlen(16,A_im),                            &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "arg_atm_fields.h"
#include "argcona.h"
#include "arglndm.h"
#include "argppx.h"
     & ErrorStatus,CMessage)

      ENDIF ! Diagnostics required for this section

      IF (L_pc2 .and. ErrorStatus == 0 .and. l_physics) then

        deallocate(t_inc_pres)
        deallocate(q_inc_pres)
        deallocate(qcl_inc_pres)
        deallocate(qcf_inc_pres)
        deallocate(cf_inc_pres)
        deallocate(cfl_inc_pres)
        deallocate(cff_inc_pres)
        deallocate(t_dini)
        deallocate(q_dini)
        deallocate(qcl_dini)
        deallocate(qcf_dini)
        deallocate(cf_dini)
        deallocate(cfl_dini)
        deallocate(cff_dini)

      End If  ! L_pc2 and ErrorStatus == 0

! section 30: climate diagnostics
      IF(      SF(0,30)                                                 &
                                ! Diagnostics required for this section
     &   .AND. ErrorStatus == 0) THEN
! size of diagnostic space
        ALLOCATE (STASHwork30(STASH_maxlen(30,A_im)))
! DEPENDS ON: st_diag3
        CALL St_diag3(STASHwork30,STASH_maxlen(30,A_im),                &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "arg_atm_fields.h"
#include "argcona.h"
#include "argppx.h"
     &    energy_corr_now,                                              &
     &    inc_u, inc_v, inc_w, inc_t,                                   &
     &    inc_q, inc_qcl, inc_qcf,                                      &
!PC2     &    inc_cf, inc_cfl, inc_cff,
     &    inc_rho,                                                      &
     &    ErrorStatus,CMessage)

        DEALLOCATE (STASHwork30) ! Clear space
        DEALLOCATE (inc_q)
        DEALLOCATE (inc_qcl)
        DEALLOCATE (inc_qcf)
        DEALLOCATE (inc_cf)
        DEALLOCATE (inc_cfl)
        DEALLOCATE (inc_cff)
      ENDIF ! Diagnostics required for this section

      If( ErrorStatus == 0) then
        DEALLOCATE (inc_rho)
        DEALLOCATE (inc_t)
        DEALLOCATE (inc_u)
        DEALLOCATE (inc_v)
        DEALLOCATE (inc_w)
      Endif

! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport 
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

! DEPENDS ON: timer
      IF(LTIMER) CALL TIMER('AS9 End TStep Diags',6)

! section 0: extraction of primary variables

! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,0,STASHWORK0_dummy,                          &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &           ErrorStatus,CMessage)
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',4)
! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,33,STASHWORK0_dummy,                         &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &           ErrorStatus,CMessage)
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',4)
! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF
! Section 34: extraction of UKCA tracer variables
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',3)
! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,34,STASHWORK0_dummy,                         &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &           ErrorStatus,CMessage)
! DEPENDS ON: timer
      IF (Ltimer) CALL TIMER('STASH',4)

! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

! ---------------------------------------------------------------
! Section 9.1  Optional diagnostic printing of w, divergence
!               lapse rate and bottom level theta
! ---------------------------------------------------------------
       if( L_diag_print .and. L_diag_print_ops ) then
! DEPENDS ON: print_ops_diag
         Call Print_ops_diag(                                           &
     &                  W, THETA, row_length, rows,                     &
     &                  model_levels, model_domain,                     &
     &                  global_row_length, global_rows,                 &
     &                  offx, offy, halo_i, halo_j,                     &
     &                  mype, nproc, gc_proc_row_group,                 &
     &                  timestep_number, print_step, diag_interval,     &
     &                  rpemax, rpemin, ipesum,                         &
     &                  L_print_pe, L_print_wmax, L_print_theta1,       &
     &                  max_w_run, min_theta1_run,                      &
     &                  time_w_max, time_theta1_min )
! DEPENDS ON: um_fort_flush
         if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
       elseif( L_diag_print ) then
! DEPENDS ON: print_diag
         Call Print_diag(                                               &
     &                  U, V, THETA,                                    &
     &                  RHO, W,                                         &
     &                  Q, QCL, QCF,                                    &
     &                  rows, n_rows, row_length,                       &
     &                  model_levels, wet_levels, model_domain,         &
     &                  global_row_length, global_rows,                 &
     &                  r_theta_levels, r_rho_levels,                   &
     &                  r_at_u, r_at_v , Pi,                            &
     &                  FV_sec_theta_latitude, FV_cos_theta_latitude,   &
     &                  cos_v_latitude,                                 &
     &                  cos_theta_longitude, sin_theta_longitude,       &
     &                  offx, offy, halo_i, halo_j,                     &
     &                  mype, nproc, at_extremity, datastart,           &
     &                  gc_proc_row_group, delta_lambda, delta_phi,     &
     &                  timestep_number, print_step, diag_interval,     &
     &                  rpemax, rpemin, ipesum, rpesum, w_print_limit,  &
     &                  L_print_pe, L_print_w,                          &
     &                  L_print_wmax, L_print_max_wind,                 &
     &                  L_print_div, L_print_lapse, L_print_theta1,     &
     &                  L_print_shear, L_diag_wind, L_diag_noise,       &
     &                  max_w_run, max_wind_run, min_theta1_run,        &
     &                  dtheta1_run, max_div_run, min_div_run,          &
     &                  min_lapse_run, max_shear_run, time_max_shear,   &
     &                  time_div_max, time_div_min, time_lapse_min,     &
     &                  time_w_max, time_max_wind, time_theta1_min,     &
     &                  max_KE_run, min_KE_run, max_noise_run,          &
     &                  time_KE_max, time_KE_min, time_noise_max )
! DEPENDS ON: um_fort_flush
         if ( L_flush6 ) call UM_FORT_FLUSH(6,info)
       endif     !  L_diag_print

! Logical first_atmstep_call is true on the first call to ATM_STEP
! and is set to false at the end of ATM_STEP (uses the SAVE command)

        first_atmstep_call = .false.

      RETURN
      END SUBROUTINE Atm_Step

#endif
