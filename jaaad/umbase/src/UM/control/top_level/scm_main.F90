

#if defined(SCMA)
! *********************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *********************************COPYRIGHT*********************************
!
!+ Routine to run the SCM.
!
! Subroutine Interface:
      Subroutine SCM_main(vert_lev, nfor, l_ts_log, ntrop, sec_day            &
               , land_points, nsprog, ntab, co2_dim_len, co2_dim_row          &
               , cloud_levels, model_levels, wet_model_levels, tr_levels      &
               , tr_vars, tr_ukca, st_levels, sm_levels                       &
               , boundary_layer_levels, ozone_levels, ntiles                  &
#include "swcarg3a.h"
#include "lwcarg3a.h"
               , atmos_sr, nsectp                                             &
               )

      Use global_SCMop
      Use scm_utils
      Use rad_switches_mod
      Use cv_cntl_mod, Only:                                                  &
          lcv_phase_lim, lcv_3d_cca, lcv_3d_ccw, lcv_ccrad, lcv_pc2_diag_sh
!dhb599 20110812: get SC from swsc.h (moved from 'downstairs')
#include "swsc.h"
      Implicit None

!
! Description:
!   Subroutine that runs the single column model.  This replaces
!         the previous version of scm_main that was the top level of
!         the model.
!
! Method:
!
! Owner: Ricky Wong
!
! Code description:
!   FORTRAN 77 + common extensions also in fortran 90.
!   This code is written to UM programming standards version 7.4

! Arguments with Intent IN. ie: Input variables.

      Integer, Intent(IN)      :: nfor
                 ! Number terms for observational forcing
      Integer, Intent(IN)      :: ntrop
                 ! Max number of levels in the troposphere
      Integer, Intent(IN)      :: CO2_DIM_LEN
      Integer, Intent(IN)      :: CO2_DIM_ROW
      Integer, Intent(IN)      :: CLOUD_LEVELS
      Integer, Intent(IN)      :: MODEL_LEVELS
      Integer, INtent(IN)      :: WET_MODEL_LEVELS
      Integer, Intent(IN)      :: TR_LEVELS
      Integer, Intent(IN)      :: TR_VARS
      Integer, Intent(IN)      :: TR_UKCA
      Integer, Intent(IN)      :: ST_LEVELS
      Integer, Intent(IN)      :: SM_LEVELS
      Integer, Intent(IN)      :: Boundary_layer_LEVELS
      Integer, Intent(IN)      :: OZONE_LEVELS
      Integer, Intent(IN)      :: NTILES
      Integer, Intent(IN)      :: NSECTP

      Integer, Intent(IN)      :: sec_day
      Integer, Intent(IN)      :: nsprog
                 ! no. of single level prognostics
      Integer, Intent(IN)      :: ntab
                 ! Dimension of array used in random
                 !  generator (Do not change this
                 !  value as it is hard coded into
                 !  the S_RANDOM deck)

      Character*200, Intent(IN) :: vert_lev
      Character*200             :: SW_spec_file
      Character*200             :: LW_spec_file

      Character(LEN=2), Intent(IN) :: ATMOS_SR(0:NSECTP)
      Logical         , Intent(IN) :: l_ts_log
                 ! Option to output timestep information

! Local variables

      Integer   :: ISTATUS
      Integer   :: height_gen_method
      Parameter(height_gen_method = 2)
      Integer  :: ICODE
      Integer  :: n_cca_levels
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      Integer  :: ErrorStatus
#endif


!=====================================================================


#include "nstypes.h"
#include "chsunits.h"
#include "cntlatm.h"
#include "cntlall.h"
#include "cmaxsize.h"
#include "s_vcoord.h"
#include "cruntimc.h"
#include "swopt3a.h"
#include "lwopt3a.h"
#include "ctlnl3a.h"
#include "c_a.h"
#include "c_pi.h"
#include "c_mdi.h"

! degrees to radians & vice versa
#include "c_g.h"
#include "c_lheat.h"
#include "c_r_cp.h"
#include "cconsts.h"
!dhb599 20110812: moved 'upstairs'
!#include "swsc.h"

! Comdecks for coriolis parameter
#include "c_omega.h"
#include "s_soilpr.h"
#include "c_epslon.h"
#include "c_lapse.h"

! Required for tropopause min/max tropopause level calculation
#include "trophgt1.h"

! for N_INTERNAL_MODEL in typsts.h
#include "csubmodl.h"
#include "typsts.h"
#include "s_maina.h"
#include "c_soilh.h"

!--------------------------------------------------------------------

      Integer nconvars
      Parameter(nconvars = 20*max_model_levels+12*max_wet_levels+       &
     &       max_tr_levels*max_tr_vars+3)

      Integer  nprimvars        ! minimum no. of  variables required
                                !  to restart from a dump and is
                                !  equal to :
      Parameter(nprimvars = 4*max_model_levels+ 4*max_wet_levels+       &
     &                     max_soil_temp_levs+max_soil_moist_levs+      &
     &                     max_nsprog)


!========================================================================
! Dummy/Fixed variables
!========================================================================
! These either have not impact on the SCM, or must be a set value
! for use with the SCM.

      Logical :: at_extremity(4) = .false.
                                      ! Indicates if this processor is at
                                      ! north, south, east or west of the
                                      ! processor grid

      Logical :: l_ohadgem1      = .false.
                                      ! HadGEM1 ocean science. Here for
                                      ! coupled runs.  Probably will not
                                      ! work with SCM.

      Integer :: neighbour(4)    = -1 ! Array with the IDs of four
                                      ! neighbours in the horizontal plane
      Integer :: g_datastart(3)  = 1
      Integer :: CycleNo         = 1  ! Define CycleNo default


      Integer, Parameter ::      &
        halo_i            = 0    &! Halo in i.
      , halo_j            = 0    &! Halo in j.
      , off_x             = 0    &! Small halo in i.
      , off_y             = 0    &! Small halo in j.
      , proc_row_group    = 1    &! Group id for processors on same row
      , proc_col_group    = 1    &! Group id for processors on same column
      , n_proc            = 1    &! Total number of processors
      , n_procx           = 1    &! Number of processors in longitude
      , n_procy           = 1    &! Number of processors in latitude
      , nice              = 1    &! No. of sea ice catagories
      , me                = 0     ! My processor number

      Real, Parameter ::         &
        lat_rot_NP        = 0.0  &! Diagnostic variable
      , long_rot_NP       = 0.0  &! Diagnostic variable
      , sd_orog_land      = 0.0  &! Ancillary fields and fields needed
      , orog_grad_xx_land = 0.0  &! to be kept from timestep to
      , orog_grad_xy_land = 0.0  &! timestep- in atmos_physics2 call
      , orog_grad_yy_land = 0.0

      Real, Parameter ::         &! STASH space for
        STASHwork3        = 0.0  &! Section 3  (b layer)
      , STASHwork5        = 0.0  &! Section 5  (convection)
      , STASHwork8        = 0.0  &! Section 8  (hydrology)
      , STASHwork9        = 0.0  &! Section 9  (LS Cloud)
      , STASHwork19       = 0.0  &! Section 19 (Veg)
      , STASHwork26       = 0.0   ! Section 26 (River Routing)


      ! Duplication of variables needed to pass UKCA gases to radiation
      ! scheme. This option will not work in the SCM.
      Integer, Parameter         :: ngrgas     = 8
      Integer, Dimension(ngrgas) :: grgas_addr = -1


      Integer ::            &
        land_pts_trif  = 1  &! For MOSES II to dimension land fields
      , npft_trif      = 1  &! For MOSES II to dimension pft  fields
      , co2_dim_lev    = 1   ! 3-D CO2 field passed to NI_rad_ctl


      ! Declare pointers for dummy arrays

      Real,    Pointer ::            &
        max_diff                     &! max diffusion coeff for run
      , lambda_p                     &! VarRes hor. co-ordinate information
      , lambda_u                     &! VarRes hor. co-ordinate information
      , dlambda_p                    &! VarRes hor. co-ordinate information
      , wt_lambda_p                  &! VarRes hor. co-ordinate information
      , wt_lambda_u                   ! VarRes hor. co-ordinate information


      Integer, Pointer ::            &
        g_p_field                    &! Size of global atmos field
      , g_r_field                    &! Size of global river field
      , a_steps_since_riv            &! Number of physics timsteps since last
                                      ! call to river routing
      , river_row_length             &! Local  river row length
      , river_rows                   &! Local  river rows
      , global_river_row_length      &! Global river row length
      , global_river_rows             ! Global river rows


      Real, Pointer, Dimension(:) :: &
        inlandout_atm                &! Inland basin flow (kg/m2/s)
      , tot_surf_runoff              &! Accumulated runoff over river
      , tot_sub_runoff               &! routing timestep (Kg/m2/s)
      , xpa                          &! Atmosphere TP long coordinates
      , xua                          &! Atmosphere U  long coordinates
      , xva                          &! Atmosphere V  long coordinates
      , ypa                          &! Atmosphere TP lat  coordinates
      , yua                          &! Atmosphere U  lat  coordinates
      , yva                           ! Atmosphere V  lat  coordinates


      Real, Pointer, Dimension(:,:) :: &
        r_area        &! Accumulated areas file
      , slope         &
      , flowobs1      &! Initialisation for flows
      , r_inext       &! x-coordinate of downstream grid point
      , r_jnext       &! y-coordinate of downstream grid point
      , r_land        &! Land/river/sea
      , substore      &! Routing sub_surface store (mm)
      , surfstore     &! Routing surface store (mm)
      , flowin        &! Rurface lateral inflow (mm)
      , bflowin       &! Rub-surface lateral inflow (mm)
      , trivdir       &! River direction file
      , trivseq       &! River sequence file
      , twatstor      &! Water storage
      , soil_clay     &! Fields for mineral dust
      , soil_silt     &! source flux calculations
      , soil_sand     &
      , dust_mrel1    &
      , dust_mrel2    &
      , dust_mrel3    &
      , dust_mrel4    &
      , dust_mrel5    &
      , dust_mrel6    &
      , phi_p         &! VarRes horizontal co-ordinate information
      , phi_v         &! VarRes horizontal co-ordinate information
      , dphi_p        &! VarRes horizontal co-ordinate information
      , wt_phi_p      &! VarRes horizontal co-ordinate information
      , wt_phi_v       ! VarRes horizontal co-ordinate information


      Real, Pointer, Dimension(:,:,:) :: &
        visc_BL_m     &! visc_m only on BL levels
      , fm_3d         &! Stability func. for momentum transport.
      , fh_3d         &! Stability func. for heat and moisture.
      , bl_coef_km    &! RHOKM from BL scheme
      , bl_coef_kh    &! RHOKH from BL scheme
      , dust_div1     &! Tracer variables
      , dust_div2     &
      , dust_div3     &
      , dust_div4     &
      , dust_div5     &
      , dust_div6

      ! Declare dummy target arrays

      Real, Target ::                                                         &
        rdum0                                                                 &
      , rdum1(row_length)                                                     &
      , rdum2(row_length+1)                                                   &
      , rdum3(row_length, rows)                                               &
      , rdum4(row_length, rows, model_levels)                                 &
      , rdum5(row_length, rows, boundary_layer_levels)                        &
      , rdum6(row_length, rows, boundary_layer_levels-1)                      &
      , dum_land(land_points)

      Integer, Target :: idum0

      Real ::                                                                 &
        frac_dummy(1,1)                                                       &
      , ES_SPACE_INTERP_dummy(4,row_length, rows)                             &
      , O3_trop_level           (row_length,rows)                             &
      , O3_trop_height          (row_length,rows)                             &
      , T_trop_level            (row_length,rows)                             &
      , T_trop_height           (row_length,rows)

!---------------------------------------------------------------------
!     &INOBSFOR information
!---------------------------------------------------------------------
!  we need to keep the variables that are read in from inobsfor and
!  reset them to their correct values at the start of each timestep so
!  use dummy equivalent variables.

      Real                                                           :: &
     &  tls         (row_length,rows,nfor,    model_levels)             &
     &, uls         (row_length,rows,nfor,    model_levels)             &
     &, vls         (row_length,rows,nfor,    model_levels)             &
     &, wls         (row_length,rows,nfor,    model_levels)             &
     &, qcl_inc     (row_length,rows,nfor,    model_levels)             &
     &, qcf_inc     (row_length,rows,nfor,    model_levels)             &
     &, q_star_keep (row_length,rows,nfor,wet_model_levels)

!--------------------------------------------------------------------
!     Time information
!--------------------------------------------------------------------
      Integer                                                           &
     &   ihour                                                          &
                                ! current hour
     &  ,imin                   ! current min

      Real                                                              &
     &   radiation_timestep                                             &
                                ! Timestep for radiation
     &  ,radiation_tstep_diag                                           &
                                  ! Diagnostic timestep
     &  ,radiation_tstep_prog     ! prognostic timestep

!---------------------------------------------------------------------
!     Initial day of the year and initial time (in seconds) in that
!     day
!---------------------------------------------------------------------

      Integer                                                           &
     &  dayno_init                                                      &
                                ! Initial day
     & ,time_initi              ! Initial time in seconds (integer)
      Real                                                              &
     &  time_init               ! Initial time in seconds (real)

!---------------------------------------------------------------------
!     Primary Model Variables plus P_star (UMDP No1)
!---------------------------------------------------------------------

      Real                                                              &
     &  r_rho_levels (row_length,rows, model_levels)                    &
     &, r_theta_levels (row_length,rows, 0:model_levels)

      Integer                                                           &
     &  iccb(row_length, rows)                                          &
                                ! Convective cloud base and top
     &  ,icct(row_length, rows) !  at levels 1 to model_levels

      Real                                                              &
     &   canopy_gb(land_points)                                         &
                                ! Canopy water content(kg/m2)
     &  ,cca(row_length, rows,model_levels)                             &
                                ! Convective cloud amount
     &  ,q(row_length, rows,wet_model_levels)                           &
                                ! Specific humidity (Kg Kg^-1)
     &  ,qcf(row_length, rows,wet_model_levels)                         &
                                ! Cloud ice content (Kg Kg^-1)
     &  ,qcl(row_length, rows,wet_model_levels)                         &
                                ! Cloud water content(Kg Kg^-1)
     &  ,inc_qcf(row_length, rows,wet_model_levels)                     &
                                ! Cloud ice content (Kg Kg^-1)
     &  ,inc_qcl(row_length, rows,wet_model_levels)                     &
                                ! Cloud water content(Kg Kg^-1)
     &  ,qcf2(row_length, rows,wet_model_levels)                        &
                                ! 2nd ice water content (Kg Kg^-1)
     &  ,qrain(row_length, rows,wet_model_levels)                       &
                                ! rain water content (Kg Kg^-1)
     &  ,qgraup(row_length, rows,wet_model_levels)                      &
                                ! graupel water content(Kg Kg^-1)
     &  ,mix_v(row_length, rows,wet_model_levels)                       &
                                ! Vapour mixing ratio (Kg Kg^-1)
     &  ,mix_cf(row_length, rows,wet_model_levels)                      &
                                ! Cloud ice content (Kg Kg^-1)
     &  ,mix_cl(row_length, rows,wet_model_levels)                      &
                                ! Cloud water content(Kg Kg^-1)
     &  ,mix_cf2(row_length, rows,wet_model_levels)                     &
                                ! 2nd ice water content (Kg Kg^-1)
     &  ,mix_rain(row_length, rows,wet_model_levels)                    &
                                ! rain water content (Kg Kg^-1)
     &  ,mix_graup(row_length, rows,wet_model_levels)                   &
                                ! graupel water content(Kg Kg^-1)
     &  ,exner_rho_levels( row_length, rows, model_levels+1)            &
     &  ,exner_theta_levels( row_length, rows, model_levels)            &
     &  ,exner_prime( row_length, rows, model_levels)                   &
                                ! increment to exner
     &  ,dtheta_dr_term( row_length, rows, model_levels)                &
     &  ,rho(row_length, rows, model_levels)                            &
     &  ,rho_n(row_length, rows, model_levels)                          &
     &  ,p(row_length, rows, model_levels+1)                            &
                                ! Pressure on rho levels
     &  ,p_theta_levels(row_length, rows, model_levels)                 &
                                !Pressure on theta levels
     &  ,rp(row_length, rows, model_levels+1)                           &
                                ! reciprocol pressure on rho levels
     &  ,rp_theta(row_length, rows, model_levels)                       &
                                !reciprocol pressure on theta levels
     &  ,p_star(row_length, rows)                                       &
                                ! Pressure at earth's surface
     &  ,smc(land_points)                                               &
                                ! Soil moisture content(Kg m^-2)
     &  ,smcl(land_points, sm_levels)                                   &
                                ! Soil moisture content in layers
                                !  (Kg m^-2)
     &  ,sthf(land_points, sm_levels)                                   &
                                ! Frozen soil moisture content
                                !  of each layer as a fraction of
                                !  saturation.
     &  ,sthu(land_points, sm_levels)                                   &
                                ! Unfrozen soil moisture content
                                !  of each layer as a fraction of
                                !  saturation.
     &  ,snodep(row_length, rows)                                       &
                                  ! Snow depth (Kg m^-2)
     &  ,t(row_length, rows,model_levels)                               &
                                           ! Temperature(K)
     &  ,t_deep_soil(land_points, st_levels)                            &
                                ! Deep soil temperatures (K)
                                !    top level not included, =surface
     &  ,tsi(row_length, rows)                                          &
                                ! Temperature of sea-ice
     &  ,tstar(row_length, rows)                                        &
                                ! Surface temperature (K)
     &  ,u(row_length, rows,model_levels)                               &
                                          ! Zonal wind (m s^-1)
     &  ,v(row_length, rows,model_levels)                               &
                                          ! Meridional wind (m s^-1)
     &  ,w(row_length, rows,0:model_levels)                             &
                                ! Vertical velocity (m s^-1)
     &  ,inc_u(row_length, rows,model_levels)                           &
                                              ! Zonal wind (m s^-1)
     &  ,inc_v(row_length, rows,model_levels)                           &
                                ! Meridional wind (m s^-1)
     &  ,inc_w(row_length, rows,0:model_levels)                         &
                                ! Vertical velocity (m s^-1)
     &  ,inc_rho(row_length, rows,model_levels)                         &
                                ! increment to rho over timestep
     &  ,inc_q(row_length, rows,wet_model_levels)                       &
                                ! increment to q over timestep
     &  ,inc_t(row_length, rows,model_levels)                           &
                                ! increment to t over timestep
     &  ,z0msea(row_length, rows)                                       &
                                  ! Sea surface roughness length
     &  ,w_adv(row_length, rows, 0:model_levels)                        &
                                  ! Advective w component of wind
     &  ,ozone_tracer(row_length, rows, model_levels)
                                  ! Ozone tracer


      Real ::                                                                 &
        theta_star     (row_length, rows,     model_levels)                   &
      , qcl_star       (row_length, rows, wet_model_levels)                   &
      , qcf_star       (row_length, rows, wet_model_levels)                   &
      , qcf2_star      (row_length, rows, wet_model_levels)                   &
      , qrain_star     (row_length, rows, wet_model_levels)                   &
      , qgraup_star    (row_length, rows, wet_model_levels)                   &
      , mix_v_star     (row_length, rows, wet_model_levels)                   &
      , mix_cl_star    (row_length, rows, wet_model_levels)                   &
      , mix_cf_star    (row_length, rows, wet_model_levels)                   &
      , mix_cf2_star   (row_length, rows, wet_model_levels)                   &
      , mix_rain_star  (row_length, rows, wet_model_levels)                   &
      , mix_graup_star (row_length, rows, wet_model_levels)                   &
      , cf_star        (row_length, rows, wet_model_levels)                   &
      , cfl_star       (row_length, rows, wet_model_levels)                   &
      , cff_star       (row_length, rows, wet_model_levels)                   &
      , uinc_geo       (row_length, rows,     model_levels)                   &
      , vinc_geo       (row_length, rows,     model_levels)

!---------------------------------------------------------------------
!     CCRad Prognostics
!---------------------------------------------------------------------
      Integer :: lcbase (row_length, rows)
                        ! Model level of lowest Convective cloud base
                        ! in profile(theta level) passed to radiation

      Real    :: ccw    (row_length, rows, wet_model_levels)
                        ! Convective Cloud Water profile passed to
                        ! radiation scheme (theta levels). (kg/kg)
                        ! Note: May be different to ccw produced by
                        !       Convection Scheme it Radiative
                        !       convective cloud decay is used.

!---------------------------------------------------------------------
!     PC2 cloud and condensation scheme
!---------------------------------------------------------------------
      Real                                                              &
     &  q_forcing   (row_length, rows, wet_model_levels)                &
     &, qcl_forcing (row_length, rows, wet_model_levels)                &
     &, t_forcing   (row_length, rows, wet_model_levels)                &
     &, p_forcing   (row_length, rows, wet_model_levels)                &
     &, q_earliest  (row_length, rows, wet_model_levels)                &
     &, qcl_earliest(row_length, rows, wet_model_levels)                &
     &, t_earliest  (row_length, rows, wet_model_levels)
!
!---------------------------------------------------------------------
!     Water
!---------------------------------------------------------------------

      Real                                                              &
     &   ls_rain(row_length, rows)                                      &
                                  ! Large scale rainfall rate (Kg*m^-2)
     &  ,ls_snow(row_length, rows)! Large scale snowfall rate
                                !  (Kg m^-2 s^-1)

!---------------------------------------------------------------------
!     Large scale statistical forcing
!---------------------------------------------------------------------

!     Random generator variables

      Integer                                                           &
     &  iv(ntab),iy,idum                                                &
                                ! Contains info on generator
     &  ,iseed                  ! Seed for random number generator

      Integer                                                           &
     &  dayno_wint              ! Day number relative to winter
                                !  solstice
      Real                                                              &
     &  ad(row_length, rows,wet_model_levels-1)                         &
                                ! Term a of equation 2.22
                                !  for dewpoint depression
     &  ,alfada(row_length, rows)                                       &
                                 ! Amplitude and mean of seasonal
     &  ,alfadb(row_length, rows)                                       &
                                 !  variation of tuning factor
     &  ,at(row_length, rows,model_levels-1)                            &
                                            ! Term a of equation 2.22
                                !  for dewpoint depression
     &  ,atime,btime                                                    &
                                ! Constants for calculating
                                !  annual cycle
     &  ,avn(row_length, rows,model_levels-1)                           &
                                             ! Term a of equation 2.22
     &  ,aw(row_length, rows,ntrop-1)                                   &
                                          !  for horiz. and vert. vel.
     &  ,cdbar(row_length, rows,wet_model_levels)                       &
                                ! Mean and SD of random variable
     &  ,cdsd(row_length, rows,wet_model_levels)                        &
                                !  for dew pt. depression
     &  ,ctbar(row_length, rows,model_levels)                           &
                                ! Mean and SD of random variable
     &  ,ctsd(row_length, rows,model_levels)                            &
                                             !  for temp.
     &  ,cvnbar(row_length, rows,model_levels)                          &
                                ! Mean and SD of random variable
     &  ,cvnsd(row_length, rows,model_levels)                           &
                                              !  for velocity VN
     &  ,cwbar(row_length, rows,ntrop)                                  &
                                       ! Mean and SD of random variable
     &  ,cwsd(row_length, rows,ntrop)                                   &
                                       !  for vertical velocity
     &  ,dbar(row_length, rows,wet_model_levels)                        &
                                ! Mean dewpoint depression at daycount
                                !  days from winter solstice (K)
     &  ,dbara(row_length, rows,wet_model_levels)                       &
                                ! Amplitude and mean of seasonal
     &  ,dbarb(row_length, rows,wet_model_levels)                       &
                                !  variation of mean dew pt.
                                !  depression (K)
     &  ,ddash(row_length, rows,wet_model_levels)                       &
                                ! Dew pt depression correction (K)
     &  ,deltan(row_length, rows)                                       &
                                   ! Radius of area (m)
     &  ,dgrada(row_length, rows,wet_model_levels)                      &
                                ! Amplitude and mean of seasonal
     &  ,dgradb(row_length, rows,wet_model_levels)                      &
                                !  variation of dew pt. depression
                                !  gradient (K km^-1)
     &  ,dsd(row_length, rows,wet_model_levels)                         &
                                                ! SD dewpoint depression
                                !  at daycount days
                                !  from winter solstice (K)
     &  ,pa(row_length, rows, model_levels+1)                           &
                                ! Amplitude and mean of seasonal
     &  ,pb(row_length, rows, model_levels+1)                           &
                                !  variation of pressure (Pa)
     &  ,px(row_length, rows,ntrop)                                     &
                                     ! Reciprocal log functions for
     &  ,py(row_length, rows,ntrop-1)                                   &
                                     !  calc. of vert. advection
                                !  used in eqns 2.12 and 2.13
     &  ,qr(row_length, rows,wet_model_levels,2)                        &
                                ! Randomly sampled specific
                                !  humidity (Kg Kg^-1)
     &  ,tbar(row_length, rows,model_levels)                            &
                                                 ! Mean temperature at
                                !  daycount days from
                                !  winter solstice (K)
     &  ,tbara(row_length, rows,model_levels)                           &
                                ! Amplitude and mean of seasonal
     &  ,tbarb(row_length, rows,model_levels)                           &
                                !  variation of temp. (K)
     &  ,tdash(row_length, rows,model_levels)                           &
                                              ! Temp correction (K)
     &  ,tgrada(row_length, rows,model_levels)                          &
                                ! Amplitude and mean of seasonal
     &  ,tgradb(row_length, rows,model_levels)                          &
                                !  variation of temp. gradient
                                !  (k Km^-1)
     &  ,tr(row_length, rows,model_levels,2)                            &
                                ! INOUT Randomly sampled temp. (K)
     &  ,tsd(row_length, rows,model_levels)                             &
                                             ! SD temp. at daycount days
                                ! from winter solstice (K)
     &  ,tsda(row_length, rows,model_levels)                            &
                                ! Amplitude and mean of seasonal
     &  ,tsdb(row_length, rows,model_levels)                            &
                                !  variation of SD of temp. (K)
     &  ,vnbar(row_length, rows,model_levels)                           &
                                ! Mean velocity VN at daycount days
                                !  from winter solstice
     &  ,vnbara(row_length, rows,model_levels)                          &
                                ! Amplitude and mean of seasonal
     &  ,vnbarb(row_length, rows,model_levels)                          &
                                !  variation of velocity VN (m s^-1)
     &  ,vnr(row_length, rows,model_levels,2)                           &
                                ! INOUT Randomly sampled horizontal
                                !    velocity (m s^-1)
     &  ,vnsd(row_length, rows,model_levels)                            &
                                ! SD velocity VN at daycount
                                !  days from winter solstice (m s^-1)
     &  ,vnsda(row_length, rows,model_levels)                           &
                                ! Amplitude and mean of seasonal
     &  ,vnsdb(row_length, rows,model_levels)                           &
                                !  variation of SD of velocity VN
                                !  (m s^-1)
     &  ,vpbar(row_length, rows,model_levels)                           &
                                ! Mean  velocity VP at daycount days
                                !   from winter solstice
     &  ,vpbara(row_length, rows,model_levels)                          &
                                ! Amplitude and mean of seasonal
     &  ,vpbarb(row_length, rows,model_levels)                          &
                                !  variation of velocity VP (m s^-1)
     &  ,vpr(row_length, rows,model_levels,2)                           &
                                ! INOUT Randomly sampled horizontal
                                !    velocity (m s^-1)
     &  ,wbar(row_length, rows,ntrop)                                   &
                                ! Mean vertical velocity at daycount
                                !  days from winter solstice (m s^-1)
     &  ,wbara(row_length, rows,ntrop)                                  &
                                ! Amplitude and mean of seasonal
     &  ,wbarb(row_length, rows,ntrop)                                  &
                                !  variation of SD of vert. vel.
                                !  (mb or HPa s^-1)
     &  ,wr(row_length, rows,ntrop,2)                                   &
                                ! INOUT Randomly sampled vertical
                                !    velocity (mb s^-1)
     &  ,wsd(row_length, rows,ntrop)                                    &
                                      ! SD vertical velocity
                                !  at daycount days from
                                !  winter solstice (m s^-1)
     &  ,wsda(row_length, rows,ntrop)                                   &
                                ! Amplitude and mean of seasonal
     &  ,wsdb(row_length, rows,ntrop)
                                !  variation of SD of vert. vel.
                                !  (mb s^-1)
                                !  roughness length (m)

!---------------------------------------------------------------------
!     Large scale observational forcing
!---------------------------------------------------------------------

!     Variables for diagnostic output for observational forcing

      Real                                                              &
     &  dap1(row_length, rows,36,model_levels)                          &
                                ! Instantaneous profiles
     &  ,dap2(row_length, rows,36,model_levels)                         &
                                ! Mean profiles
     &  ,dap3(row_length, rows,36,nfor-1,model_levels)                  &
                                ! Mean profiles - timeseries
     &  ,dab1(row_length, rows,44)                                      &
                                   ! Instantaneous budgets
     &  ,dab2(row_length, rows,44)                                      &
                                   ! Mean budgets
     &  ,dab3(row_length, rows,44,nfor-1)                               &
                                ! Mean budgets - timeseries
     &  ,deltap(row_length, rows, model_levels)                         &
                                ! Layer thickness (Pa)
     &  ,factor_rhokh(row_length, rows)
                                ! used to specify surface flux
                                !  from observation


      ! Dump array of restart variables
      Real resdump(row_length, rows,nprimvars)

!---------------------------------------------------------------------
!     Variables enabling diagnostic output
!---------------------------------------------------------------------

      integer main_diag_switch ! If zero -> no diagnostic output
      integer netcdf_chunksize ! The ChunkSize input to NF90_Create()
      Integer, Parameter :: nSCMDpkgs = 12 ! No of diags packages
      ! Note that if nSCMDpkgs is changed it should also be changed
      ! in PC2_PRES and ATMSTEP2

      Logical L_SCMDiags(nSCMDpkgs) ! Logicals for diagnostics packages

! Parameters used in calls to SCMoutput
#include "s_scmop.h"
! n_vis_thresh in c_visbty.h is required for the call to setup_diags
#include "c_visbty.h"

      integer site(row_length*rows)           ! SSFM site WMO number

! Arrays to store initial fields for calculating total increments and
! total increment fields
      Real                                                              &
     &  t_start(row_length,rows,model_levels)                           &
     &, u_start(row_length,rows,model_levels)                           &
     &, v_start(row_length,rows,model_levels)                           &
     &, w_start(row_length,rows,model_levels)                           &
     &, q_start(row_length,rows,wet_model_levels)                       &
     &, qcl_start(row_length,rows,wet_model_levels)                     &
     &, qcf_start(row_length,rows,wet_model_levels)                     &
     &, t_totalinc(row_length,rows,model_levels)                        &
     &, u_totalinc(row_length,rows,model_levels)                        &
     &, v_totalinc(row_length,rows,model_levels)                        &
     &, w_totalinc(row_length,rows,model_levels)                        &
     &, q_totalinc(row_length,rows,wet_model_levels)                    &
     &, qcl_totalinc(row_length,rows,wet_model_levels)                  &
     &, qcf_totalinc(row_length,rows,wet_model_levels)

!---------------------------------------------------------------------
!     Extra diagnostics for MOSES boundary layer code
!---------------------------------------------------------------------

      Real                                                              &
     &  resp_w_ft(row_length, rows,npft)
                                ! OUT Wood maintenance respiration
                                !     (Kg C m^-2 s^-1).

!---------------------------------------------------------------------
!     Clouds
!---------------------------------------------------------------------

      Real                                                              &
     &  layer_cloud(row_length, rows,wet_model_levels)                  &
                                ! layer cloud amount
                                !  (decimal fraction)
     &  ,ccwpin(row_length, rows)                                       &
                                 ! Condensed water path Kg m^-2
     &  ,RHTS(row_length,rows,wet_model_levels)                         &
                                 ! Initial relative humidity wrt TL
     &  ,qtts(row_length,rows,wet_model_levels)                         &
                                ! initial q_T
     &  ,tlts(row_length,rows,wet_model_levels)                         &
                                ! initial TL
     &  ,ptts(row_length,rows,wet_model_levels)                         &
                                ! initial pressure on theta levels
     &  ,TL(row_length,rows)                                            &
                                 ! Liquid water temperature (K)
     &  ,QSL_TL(row_length,rows) ! Qsat wrt liquid water at temp TL

!---------------------------------------------------------------------
!     Radiation
!---------------------------------------------------------------------

      Integer                                                           &
     &  daynumber                                                       &
                                ! Day in the year (default=1)
     &  ,year                   ! Year

      Real                                                              &
     &   lwout(row_length, rows,model_levels+1)
                                ! Longwave atmospheric heating
                                !  rates in levels 2,model_levels+1
                                !  (K/timestep). NET LW flux
                                !  in level 1
                                !  If sea point LWOUT(1) contains
                                !  net longwave flux over land
                                !  portion (land or land ice) and
                                !  LWSEA that flux over sea portion
!---------------------------------------------------------------------
!     Used in calculation of min/max tropopause
!---------------------------------------------------------------------
      Real r_ref_rho(model_levels)  ! height of model rho levels
                                    ! above surface.

!--------------------------------------------------------------------
!     Loop Counters and limits
!--------------------------------------------------------------------

      Integer                                                           &
     &  daycount                                                        &
                                ! Counts through days
     &  ,daysteps                                                       &
                                ! Number of timestep in a day
     &  ,full_daysteps                                                  &
                                ! Number of timesteps in a full day
     &  ,nstepsin                                                       &
                                ! Number of steps in final day
     &  ,total_nsteps                                                   &
                                ! Total number of steps in run
     &  ,i,j,j2,k                                                       &
                                ! General loop counters
                                !  array dumpmean
     &  ,itype                                                          &
     &  ,m1                                                             &
                                ! No. of dumps
     &  ,stepcount                                                      &
                                ! Counts through timesteps
!EAK 
     &  ,istep_cur              ! current time step  

        save  istep_cur

!--------------------------------------------------------------------
!     Miscellaneous
!---------------------------------------------------------------------

      CHARACTER*(80)                                                    &
     &       CMESSAGE              ! Error message if ICODE >0
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='SCM_MAIN')

      Character*8                                                       &
     &  time_string             ! String containing actual time

      Integer ::                                                        &
        Error_code                                                      &
      , day                                                             &
      , yearno

      Real                                                              &
     &  time_sec                                                        &
                                ! actual time of day in secs.
     &  ,modug(row_length, rows)                                        &
                                 ! Magnitude of Geostrophic wind
                                !  (m s^-1)
     &  ,f_coriolis(row_length, rows)                                   &
                                      ! 2*omega*sin(latitude)
     &  ,f3_at_u(row_length, rows)                                      &
     &  ,rccb(row_length, rows)                                         &
                                ! Real val. cloud base geoint only
     &  ,rcct(row_length, rows)                                         &
                                ! Real cloud top geoint only
     &  ,maxinc                                                         &
                                ! Maximum wind increment from geoinit
     &, SW_incs(row_length, rows, 0:model_levels+1)                     &
     &, LW_incs(row_length, rows, 0:model_levels)                       &
     &, dirpar_incs(row_length, rows)                                   &
                                       ! PAR flux variable
     &, t1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, q1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, SW_tile_rts(row_length*rows, ntiles)                            &
                                ! Surface net SW on land tiles

     &, long_rad(row_length, rows)                                      &
                                  ! longitude in radians used for
                                   ! diagnostic outputting
     &, lat_rad(row_length, rows)                                       &
                                 ! longitude in radians used for
                                   ! diagnostic outputting
     &  ,rhokh(row_length, rows,boundary_layer_levels)
                                ! Exchange coeffs for moisture.
                                ! surface:out of SF_EXCH contains
                                ! contains only RHOKH.
                                ! above surface:out of KMKH cont-
                                ! ains GAMMA(1)*RHOKH(,1)*RDZ(,1)

      Real :: gridbox_area_m(row_length, rows)
                                ! Gridbox area in m^2

      Integer                                                           &
     &  land_points                                                     &
                     ! no of land points being processed - can be 0
     &  ,land_ice_points                                                &
                         ! no of ice land points
     &  ,soil_points                                                    &
                      ! no of soil points
     &  ,land_index(row_length*rows) !Shd be defined by land_points
                                     !but land_points has not yet been
                                     !given a value
      Integer :: land_ice_index(row_length*rows)
                                     !Shd be defined by land_ice_points
                                     !but land_points has not yet been
                                     !given a value
      Integer :: soil_index(row_length*rows)
                                     !Shd be defined by soil_points
                                     !but soil_points has not yet been
                                     !given a value

      ! Land point specific surface soil parameters
      Real :: b_exp  (land_points) ! Clapp-Hornberger exponent
      Real :: hcap   (land_points) ! Soil heat capacity
      Real :: hcon   (land_points) ! Soil thermal conductivity
      Real :: satcon (land_points) ! Saturated hydrological conductivity
      Real :: sathh  (land_points) ! Saturated soil water suction
      Real :: v_sat  (land_points) ! Volumetric SMC at saturation
      Real :: v_wilt (land_points) ! Volumetric SMC at wilting point
      Real :: v_crit (land_points) ! Volumetric SMC at critical point


      Real                                                              &
     &  area_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, bulk_cloud_fraction(row_length, rows, wet_model_levels)         &
     &, cloud_fraction_liquid(row_length, rows, wet_model_levels)       &
     &, cloud_fraction_frozen(row_length, rows, wet_model_levels)

      Real                                                              &
     &  energy_correction                                               &
     &, free_tracers( row_length, rows, tr_levels, tr_vars )            &
     &, ukca_tracers( row_length, rows, tr_levels, tr_ukca )

      Real                                                           :: &
     & biogenic(row_length, rows, model_levels)

! Local arrays for using the aerosol climatology for NWP

#include "arcl_dim.h"
#include "arcl_ids.h"

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

      Real :: unscl_dry_rho(row_length, rows, model_levels)
              ! unscaled dry density

      Real                                                              &
     &  rhokm(row_length, rows, 0:boundary_layer_levels-1)              &
     &, cH_term(row_length, rows, model_levels-1)


!     EAK
      logical, dimension(row_length*rows,NTILES) :: L_TILE_PTS

      Real                                                              &
     &  SNOW_DEPTH3L(row_length*rows,NTILES,3)                          &
     &, SNOW_MASS3L(row_length*rows,NTILES,3)                           &
     &, SNOW_COND(row_length*rows,NTILES,3)                             &
     &, SNOW_TMP3L(row_length*rows,NTILES,3)                            &
     &, SNOW_RHO3L(row_length*rows,NTILES,3)                            &
     &, SNOW_RHO1L(row_length*rows,NTILES)                              &
     &, SNAGE_TILE(row_length*rows,NTILES)                              &
     &, SMCL_TILE(row_length*rows,NTILES,sm_levels)                     &
     &, STHU_TILE(row_length*rows,NTILES,sm_levels)                     &
     &, STHF_TILE(row_length*rows,NTILES,sm_levels)                     &
     &, TSOIL_TILE(row_length*rows,NTILES,sm_levels)                    &
     &, T_SURF_TILE(row_length*rows,NTILES)                             &
     &, RTSOIL_TILE(row_length*rows,NTILES)                             &
     &, GFLUX_TILE(row_length*rows,NTILES)                              &
     &, SGFLUX_TILE(row_length*rows,NTILES)                             &
     &, HCONS(row_length*rows)

!     EAK
      Integer ISNOW_FLG3L(row_length*rows,NTILES)                      
      ! MRD - does this need to be initialised?
      Integer veg_type(row_length*rows)                      

      Real                                                              &
     &  photosynth_act_rad(row_length, rows)                            &
                                             ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
!     EAK
!     *** note in ni-rad_ctl & glue_rad alb_tile(land_field,ntiles) dim ***

     &, alb_tile(row_length*rows,ntiles,4)                              &
     &, surf_down_sw(row_length,rows,4)                                 &
     &,FTL_TILE_CAB(row_length*rows,NTILES)                             &
     &,LE_TILE_CAB(row_length*rows,NTILES)                              &
     &,LE_CAB(row_length*rows)                                          &
     &,TSTAR_TILE_CAB(row_length*rows,NTILES)                           &
     &,TSTAR_CAB(row_length*rows)                                       &
     &,SMCL_CAB(row_length*rows,SM_LEVELS)                             &
     &,TSOIL_CAB(row_length*rows,SM_LEVELS)                            &
     &,USTAR_CAB(row_length*rows)                                       &
     &,SURF_HTF_CAB(row_length*rows)                                    &
!sxy
     &, FTL_CAB(row_length*rows)                                        &
     &, TOT_ALB(row_length*rows,ntiles)                                 &
     &, CH_CAB(row_length*rows)                                         &
     &, CD_CAB(row_length*rows)                                         &
     &, U_S_CAB(row_length*rows)                                        &
     &, CD(row_length*rows)                                             &
     &, CH(row_length*rows)                                             &
     &, LAND_ALBEDO(row_length,rows,4)                                  & 
     &,SW_DOWN(row_length*rows)                                         &
                                   ! Surface downward SW radiation (W/m2).
     &,radnet_tile(row_length*rows,ntiles)                              &
!     &, rad_hr(row_length, rows, boundary_layer_levels, 2)
!                                  ! BL radiative heating rates
! dimesions of rad_hr changed in 5.4
     &, rad_hr(row_length, rows, boundary_layer_levels, 2)              &
                                  ! BL radiative heating rates
     &, micro_tends(row_length, rows, boundary_layer_levels, 2)         &
                                  ! BL Microphys tendencies
     &, surf_radflux(row_length, rows)                                  &
     &, dOLR(row_length, rows)                                          &
                                    ! TOA - surface upward LW
     &, SW_tile(row_length*rows, ntiles)                                &
                                    ! Surface net SW on land tiles
     &, cos_zenith_angle(row_length, rows)

      !=============================================
      ! Rate of changes in forcing increments (1/s)
      !=============================================
      Real                                                           :: &
     &  ch_t_inc (row_length,rows,nfor-1,    model_levels)              &
                                        ! Temperature
     &, ch_u_inc (row_length,rows,nfor-1,    model_levels)              &
                                        ! Zonal wind
     &, ch_v_inc (row_length,rows,nfor-1,    model_levels)              &
                                        ! Meridional wind
     &, ch_w_inc (row_length,rows,nfor-1,    model_levels)              &
                                        ! Vertical wind
     &, ch_q_star(row_length,rows,nfor-1,wet_model_levels)              &
                                        ! Specific Humidity
     &, ch_flux_e(row_length,rows,nfor-1)                               &
                                        ! Latent heat flux
     &, ch_flux_h(row_length,rows,nfor-1)
                                        ! Sensible heat flux

      Integer                                                           &
     & TILE_INDEX(row_length*rows,NTYPE) &! OUT Index of tile points         
     &,TILE_PTS(NTYPE)                    ! OUT Number of tile points

      Logical ::         &
        L_Rad_Step       &! True if radiation  timestep
      , L_Rad_Step_Diag  &! True if diagnostic timestep
      , L_Rad_Step_Prog  &! True if prognostic timestep
      , l_mr_pc2          ! True if using mixing ratios in PC2 section


      ! Land based variables

      Real ::                    &
        fexp       (land_points) &! Decay factor in Sat. conductivity in water
                                  ! table layer.
      , gamtot     (land_points) &! Integrated complete Gamma function.
      , ti_mean    (land_points) &! Mean topographic index.
      , ti_sig     (land_points) &! St dev. of topographic index. in water
                                  ! table layer
      , fsat       (land_points) &! Surface saturation fraction.
      , fwetl      (land_points) &! Wetland fraction.
      , zw         (land_points) &! Water table depth (m).
      , sthzw      (land_points) &! Soil moist fract. in deep-zw layer.
      , a_fsat     (land_points) &! Fitting parameter for Fsat in LSH model
      , c_fsat     (land_points) &! Fitting parameter for Fsat in LSH model
      , a_fwet     (land_points) &! Fitting parameter for Fwet in LSH model
      , c_fwet     (land_points) &! Fitting parameter for Fwet in LSH model
      , catch_snow (land_points) &! Snow capacity for NLT tile (kg/m2)
      , snow_grnd  (land_points)  ! Snow below canopy (kg/m2).


      Character(200) :: sdum0, sdum1

      Integer :: asteps_since_triffid
      Integer :: previous_time(7)      ! Time information for current timestep

      ! Trig arrays
      Real ::                                                                 &
        cos_latitude    (row_length, rows)                                    &
      , sec_latitude    (row_length, rows)                                    &
      , cos_longitude   (row_length, rows)                                    &
      , sin_latitude    (row_length, rows)                                    &
      , sin_longitude   (row_length, rows)                                    &
      , FV_cos_latitude (row_length, rows)


      Real :: delta_lambda, delta_phi


      !  Ancillary fields and fields needed to be kept from timestep to
      ! timestep

      Real ::                                &
        ice_fract_n(row_length, rows, nice)  &! Ice catergories
      , ice_thick_n(row_length, rows, nice)  &
      , ti_n       (row_length, rows)         ! Ice temperature as grib box
                                              ! mean

      ! Add new variables from extra args in calls to physics routines
      Real :: rhcpt(row_length, rows, wet_model_levels)

#include "s_mainn.h"
#include "s_maind.h"

      ! Assign Pointers to dummy target arrays
      g_p_field               => idum0
      g_r_field               => idum0
      a_steps_since_riv       => idum0
      river_row_length        => idum0
      river_rows              => idum0
      global_river_row_length => idum0
      global_river_rows       => idum0
      max_diff                => rdum0
      lambda_p                => rdum0
      phi_p                   => rdum3
      lambda_u                => rdum0
      phi_v                   => rdum3
      dlambda_p               => rdum0
      dphi_p                  => rdum3
      wt_lambda_p             => rdum0
      wt_lambda_u             => rdum0
      wt_phi_p                => rdum3
      wt_phi_v                => rdum3
      dust_div1               => rdum4
      dust_div2               => rdum4
      dust_div3               => rdum4
      dust_div4               => rdum4
      dust_div5               => rdum4
      dust_div6               => rdum4
      r_area                  => rdum3
      slope                   => rdum3
      flowobs1                => rdum3
      r_inext                 => rdum3
      r_jnext                 => rdum3
      r_land                  => rdum3
      substore                => rdum3
      surfstore               => rdum3
      flowin                  => rdum3
      bflowin                 => rdum3
      trivdir                 => rdum3
      trivseq                 => rdum3
      twatstor                => rdum3
      soil_clay               => rdum3
      soil_silt               => rdum3
      soil_sand               => rdum3
      dust_mrel1              => rdum3
      dust_mrel2              => rdum3
      dust_mrel3              => rdum3
      dust_mrel4              => rdum3
      dust_mrel5              => rdum3
      dust_mrel6              => rdum3
      visc_bl_m               => rdum5
      fm_3d                   => rdum5
      fh_3d                   => rdum5
      bl_coef_km              => rdum6
      bl_coef_kh              => rdum6
      inlandout_atm           => dum_land
      tot_surf_runoff         => dum_land
      tot_sub_runoff          => dum_land
      xpa                     => rdum2
      xua                     => rdum2
      xva                     => rdum2
      ypa                     => rdum1
      yua                     => rdum1
      yva                     => rdum2

      ! Initialise Dummy/Hardwired values

      ice_fract_n  = 0.0    ! Was uninitialised before
      ice_thick_n  = 0.0    ! Was uninitialised before
      ti_n         = 273.0  ! Ice temperature as grib box mean
      model_domain = 5      ! Specifies model type
      n_rims_to_do = 1      ! Size of LAM rim zone
      l_regular    = .true. ! True if NOT variable resolution
      NumCycles    = 1      ! cycling dynamics-physics to be compatible
                            ! with UM
      co2emitmass = 0.0

      turb_startlev_vert = 1
      turb_endlev_vert   = 1
      river_step         = 1.0
      river_vel          = 1.0
      river_mcoef        = 1.0
      fexp               = 1.0
      gamtot             = 1.0
      ti_mean            = 1.0
      ti_sig             = 1.0
      fsat               = 1.0
      fwetl              = 1.0
      zw                 = 1.0
      sthzw              = 1.0
      a_fsat             = 1.0
      c_fsat             = 1.0
      a_fwet             = 1.0
      c_fwet             = 1.0
      catch_snow         = 1.0
      snow_grnd          = 1.0


      ! Hardwired model switches for SCM
      l_mr_pc2     = .false. ! Use mixing ratios in PC2 section
      l_gwd        = .false. ! Atmos_physics1 call
      l_murk_bdry  = .false. ! Atmos_physics1 call
      lflux_reset  = .false. ! Atmos_physics1 call
      l_rad_deg    = .false. ! Atmos_physics1 call
      l_RHCPT      = .false. ! Atmos_physics2 call
      l_inland     = .false. ! Re-route inland basin flow to soil moisture

      ! Set dimensions of aerosol climatologies - aerosols are not used
      ! in the SCM so this are hardwired to be unused.
      !
      ! If aerosols are enabled in the furute then these should be correctly
      ! dimensioned using set_arcl_dimensions and set_arcl_clim as in
      ! atm_step
      n_arcl_species=0
      n_arcl_compnts=1   ! Set to one to prevent
                         ! out of bounds error in radiation

      L_USE_ARCL(IP_ARCL_BIOM) = .false.
      L_USE_ARCL(IP_ARCL_BLCK) = .false.
      L_USE_ARCL(IP_ARCL_SSLT) = .false.
      L_USE_ARCL(IP_ARCL_SULP) = .false.
      L_USE_ARCL(IP_ARCL_DUST) = .false.
      L_USE_ARCL(IP_ARCL_OCFF) = .false.
      L_USE_ARCL(IP_ARCL_DLTA) = .false.
      Allocate ( arcl(1,1,1,1) )


      ! Mid Latitude values, (ref J-C-Thelen)
      O3_trop_level  = 22
      O3_trop_height = 11000
      T_trop_level   = 22
      T_trop_height  = 11000


      ! Initialise dummy arguments
      rdum0    = 1.0
      rdum1    = 1.0
      rdum2    = 1.0
      rdum3    = 1.0

      rdum4    = 1.0
      rdum5    = 1.0
      rdum6    = 1.0

      idum0    = 1
      dum_land = 1.0

      Write(sdum0,*) nfor
      Write(sdum1,*) max_nfor

      ! Set convection switches in cv_cntl_mod
      lcv_phase_lim    = l_phase_lim
      lcv_3d_cca       = l_3d_cca
      lcv_3d_ccw       = l_3d_ccw
      lcv_ccrad        = l_ccrad
      lcv_pc2_diag_sh  = l_pc2_diag_sh

!      write(6,*)'XXX SCM_main: SC = ', SC

!=====================================================================
!  Set up vertical level information
!=====================================================================

      Open(10, File=vert_lev, IOstat=istatus, Status='old')

      If (istatus /= 0) Then

        Icode=500
        Write(Cmessage,*)  " Error opening " //TRIM(ADJUSTL(vert_lev))//     &
                           " file on unit 10"

! DEPENDS ON: ereport
        Call Ereport(RoutineName,Icode,Cmessage)

      End If

      READ(10,VERTLEVS)

      CLOSE(10)

!     Initialise error code

      Error_code = 0

!-------------------------------------------------------------------
!     Set namelist variables to default values.
!-------------------------------------------------------------------

! &INDATA
!=========

      gridbox_area = 100000
      lat          = 0
      long         = 0
      soil_type    = IMDI

! &INOBSFOR
!===========

      flux_h_scm        = 0.0
      flux_e_scm        = 0.0
      tstar_forcing_scm = 0.0

      q_star_scm = 0.0
      t_inc_scm  = 0.0
      u_inc_scm  = 0.0
      v_inc_scm  = 0.0
      w_inc_scm  = 0.0
      qcl_inc    = 0.0
      qcf_inc    = 0.0

! &INGEOFOR
!===========

      ug = 5.0
      vg = 5.0

! &INMOSES
!=========

      smi_opt    = IMDI
      frac_typ   = RMDI
      smcli      = RMDI
      sth        = RMDI

      gs         = 1.0
      fsmc       = RMDI
      canopy     = 0.0
      rgrain     = 50.0

      canht      = RMDI
      lai        = RMDI
      catch      = RMDI

      snow_tile  = 0.0
      z0_tile    = RMDI
      infil_tile = RMDI
      tstar_tile = RMDI


      ! 2A Triffid variables
      !=====================
      frac_disturb     = RMDI
      npp_ft_acc       = RMDI
      resp_w_ft_acc    = RMDI
      resp_s_acc       = RMDI
      cs               = RMDI
      g_leaf_phen_acc  = RMDI
      g_leaf_acc       = RMDI


! &INPROF
!=========

      theta_scm = 0.0
      p_in_scm  = 0.0
      ti_scm    = 0.0
      ui_scm    = 0.0
      vi_scm    = 0.0
      wi_scm    = 0.0
      qi_scm    = 0.0

      t_deep_soili = RMDI
      canopy_gbi   = RMDI
      smci         = RMDI

      di         = 0       ! These model variables are set up as
      ice_fract  = 0.0     ! constants for a land point.
      snodepi    = 0.0
      tstari     = 0.0
      z0mseai    = 0.0001
      z0m_scm    = 0.0
      z0h_scm    = 0.0
      ccai       = 0.0
      iccbi      = 0
      iccti      = 0
      u_0        = 0.0
      v_0        = 0.0
      tracer     = 0.0



! &RUNDATA
!==========

      ntml    = boundary_layer_levels
      zh      = 500.0
      ls_rain = 0.0
      ls_snow = 0.0

      so2        = 0.0
      so4_aitken = 0.0
      so4_accu   = 0.0
      so4_diss   = 0.0
      aerosol    = 0.0

      If (land_points > 0) Then
        albsoil = RMDI
      End If


! &RADCLOUD
!==========


      layer_cloud_rad = 0.0
      qcl_rad         = 0.0
      qcf_rad         = 0.0

      cca_rad    = 0.3132
      ccwpin_rad = 0.0
      ccwpin     = 1.0
      iccb_rad   = 2
      icct_rad   = 8


! Set radiation switches in rad_switches_mod
      lrad_ctile_fix     = l_ctile_fix
      lrad_cldtop_t_fix  = l_cldtop_t_fix
      lrad_ccrad         = l_ccrad
      lrad_3d_ccw        = l_3d_ccw
      lrad_ovrlap        = l_ovrlap
      lrad_ccw_scav      = l_ccw_scav
      lrad_emis_land_gen = l_emis_land_gen
      rad_emis_land_gen  = emis_land_gen

!---------------------------------------------------------------------
!     Initialise the array giving the unit nos. for output of the
!     initial data and anything else that needs setting
!---------------------------------------------------------------------

      asteps_since_triffid = 0
      omega                = 7.292116E-5
      two_omega            = 2 * 7.292116E-5

! Initialise

      rad_hr       = 0.0
      micro_tends  = 0.0
      surf_radflux = 0.0

      ! ilscnt hardwired initialisation to zero until further
      ! development to correct functionality. Also removed from
      ! &INOBSFOR namelist to prevent incorrect usage.

      ilscnt            = 0


!---------------------------------------------------------------------
!     Read in NAMELISTS
!---------------------------------------------------------------------
        Open(10, File='namelist.scm', IOSTAT=ISTATUS)
        If(Istatus  /=  0) Then
          Icode=500
          Write(6,*) ' ERROR OPENING FILE ON UNIT 10'
          Write(6,*) ' FILENAME = namelist.scm'
          Write(6,*) ' IOSTAT =',ISTATUS
          GOTO 999
        End If


      Read(10,INDATA)

      Read(10,RUNDATA)
      If (int(sec_day/timestep) /= float(sec_day)/timestep) Then
         ! Print a warning if we're going to run for more
         ! than a day.
         If (ndayin+float(nminin)/60./24.+                              &
     &        float(nsecin)/3600./24. >  1.) Then
           print*,'****************************************************'
           print*,'* Warning: your timestep is such that you have     *'
           print*,'* requested a non-integer number of steps per day. *'
           print*,'****************************************************'
         End If
      End If
      full_daysteps = int(sec_day/timestep)
      nstepsin = int((nminin*60 + nsecin)/timestep)
      total_nsteps = ndayin*full_daysteps+nstepsin
      scm_timestep = timestep


      Do i = 1, row_length
        Do j = 1, rows
          gridbox_area_m(i,j) = gridbox_area(i,j) * 1e+6
          deltan(i,j)         = sqrt(gridbox_area_m(i,j)/pi)
        End Do
      End Do

      Read(10,LOGIC)

      If    ((stats    .AND. obs)     .OR.(stats .AND. noforce)         &
     &  .OR. (stats    .AND. geoforce).OR.(obs   .AND. geoforce)        &
     &  .OR. (geoforce .AND. noforce) .OR.(obs   .AND. noforce)) Then

        write (*,'(4(TR1,A52/),4(TR1,A12,TR1,L1,T53,A1/),4(TR1,A52/))') &
     &        '===================================================='    &
     &,       '| Warning: More than one forcing logical is set to |'    &
     &,       '| true. This may produce unexpected results.       |'    &
     &,       '|                                                  |'    &
     &,       '|  noforce: ', noforce,                           '|'    &
     &,       '| geoforce: ', geoforce,                          '|'    &
     &,       '|    stats: ', stats,                             '|'    &
     &,       '|      obs: ', obs,                               '|'    &
     &,       '|                                                  |'    &
     &,       '===================================================='

      End If ! Test for multiple forcing options


      If     ((.NOT. stats)   .AND. (.NOT. obs)                         &
     &  .AND. (.NOT. noforce) .AND. (.NOT. geoforce)) Then

        write (*,'(4(TR1,A52/))')                                       &
     &        '===================================================='    &
     &,       '| No forcing logical set.                          |'    &
     &,       '| Setting noforce to true as default               |'    &
     &,       '===================================================='
        noforce = .TRUE.

      End If ! Test for at least one forcing option.

      If (L_use_dust .OR. L_dust) Then
        write (*,'(8(TR1,A52/))')                                       &
     &        '===================================================='    &
     &,       '| Mineral dust scheme not yet implemented in SCM.  |'    &
     &,       '| Setting following logicals to false:             |'    &
     &,       '|                                                  |'    &
     &,       '|  L_use_dust                                      |'    &
     &,       '|  L_dust                                          |'    &
     &,       '|                                                  |'    &
     &,       '===================================================='
        L_use_dust = .FALSE.
        L_dust     = .FALSE.

      End If ! Test for dust scheme options


      If    (L_use_soot_direct   .OR. L_use_soot_indirect               &
     &  .OR. L_use_soot_autoconv .OR. L_soot) Then
        write (*,'(8(TR1,A52/))')                                       &
     &        '===================================================='    &
     &,       '| Soot chemistry not yet implemented in the SCM.   |'    &
     &,       '| Setting following logicals to false:             |'    &
     &,       '|                                                  |'    &
     &,       '|  L_use_soot                                      |'    &
     &,       '|  L_soot                                          |'    &
     &,       '|                                                  |'    &
     &,       '===================================================='
        L_soot              = .FALSE.
        L_use_soot_direct   = .FALSE.
        L_use_soot_indirect = .FALSE.
        L_use_soot_autoconv = .FALSE.

      End If ! Test for soot chemistry options


      If    (l_biomass            .OR. l_use_bmass_direct               &
     &  .OR. l_use_bmass_indirect .OR. l_use_bmass_autoconv) Then

        write (*,'(10(TR1,A52/))')                                      &
     &        '===================================================='    &
     &,       '| Biomass scheme not available in the SCM.         |'    &
     &,       '| Setting following logicals to false:             |'    &
     &,       '|                                                  |'    &
     &,       '|  L_biomass                                       |'    &
     &,       '|  L_use_bmass_direct                              |'    &
     &,       '|  L_use_bmass_indirect                            |'    &
     &,       '|  L_use_bmass_autoconv                            |'    &
     &,       '|                                                  |'    &
     &,       '===================================================='
        L_biomass             = .FALSE.
        L_use_bmass_direct    = .FALSE.
        L_use_bmass_indirect  = .FALSE.
        L_use_bmass_autoconv  = .FALSE.

      End If

      If    (l_ocff               .OR. l_use_ocff_direct                &
     &  .OR. l_use_ocff_indirect  .OR. l_use_ocff_autoconv) Then

        write (*,'(10(TR1,A52/))')                                      &
     &        '===================================================='    &
     &,       '| Fossil-fuel OC scheme not available in the SCM.  |'    &
     &,       '| Setting following logicals to false:             |'    &
     &,       '|                                                  |'    &
     &,       '|  L_ocff                                          |'    &
     &,       '|  L_use_ocff_direct                               |'    &
     &,       '|  L_use_ocff_indirect                             |'    &
     &,       '|  L_use_ocff_autoconv                             |'    &
     &,       '|                                                  |'    &
     &,       '===================================================='
        L_ocff                = .FALSE.
        L_use_ocff_direct     = .FALSE.
        L_use_ocff_indirect   = .FALSE.
        L_use_ocff_autoconv   = .FALSE.

      End If

      ! No of levels for Convective Cloud Amount.
      If (lcv_3d_cca) Then
        n_cca_levels = wet_model_levels
      Else
        n_cca_levels = 1
      End If

      If (obs) Then

        If (nfor == -999) Then
          Icode=503
          Write(6,*) " "
          Write(6,*) "============================================="
          Write(6,*) " Number of observational forcings (nfor) in  "
          Write(6,*) " SCM namelist (&CNTLSCM) has not been set.   "
          Write(6,*) "============================================="
          Write(6,*) " "

          CLOSE(10)

          Write(Cmessage,*) ' Variable NFOR has not been set'

! DEPENDS ON: ereport
          Call Ereport(RoutineName,Icode,Cmessage)

        Else If (nfor > max_nfor) Then
          Icode=504
          Write(6,*) " "
          Write(6,*) "============================================="
          Write(6,*) " Specified nfor(" // TRIM(ADJUSTL(sdum0)) //              &
                     ") > max_nfor("    // TRIM(ADJUSTL(sdum1)) // ")."
          Write(6,*) " This WILL produce incorrect forcings.       "
          Write(6,*) " Please check your namelist.scm "
          Write(6,*) "============================================="
          Write(6,*) " "

          CLOSE(10)

          Write(Cmessage,*)                                                     &
                     " Specified nfor(" // TRIM(ADJUSTL(sdum0)) //              &
                     ") > max_nfor("    // TRIM(ADJUSTL(sdum1)) // ")."

! DEPENDS ON: ereport
          Call Ereport(RoutineName,Icode,Cmessage)

        End If
      End If !  obs

      !-------------------------------------------------------------------------
      ! Write out user notification of ancillary files
      !-------------------------------------------------------------------------

      ! Get the SW spectral file location from environment variable
      ! "UNIT57"
      CALL FORT_GET_ENV('UNIT57',6,SW_spec_file,200,istatus)
      If (istatus /= 0) Then
        icode=506
        Cmessage='Environment variable UNIT57 (SW spectral file '// &
     &     'location) not set.'
        Call Ereport(RoutineName,Icode,Cmessage)
      End If

      ! Get the LW spectral file location from environment variable
      ! "UNIT80"
      CALL FORT_GET_ENV('UNIT80',6,LW_spec_file,200,istatus)
      If (istatus /= 0) Then
        icode=507
        Cmessage='Environment variable UNIT80 (LW spectral file '// &
     &     'location) not set.'
        Call Ereport(RoutineName,Icode,Cmessage)
      End If

      Write(6,*) ' '
      Write(6,*) ' Using Ancillary files:'
      Write(6,'(A100)')                                                       &
                 ' ---------------------------------------------------'//     &
                 '----------------------------------------------------'
      Write(6,*) ' Vertical levels:'//TRIM(ADJUSTL(vert_lev))
      Write(6,*) ' SW spectral    :'//TRIM(ADJUSTL(sw_spec_file))
      Write(6,*) ' LW spectral    :'//TRIM(ADJUSTL(lw_spec_file))
      Write(6,'(A100)')                                                       &
                 ' ---------------------------------------------------'//     &
                 '----------------------------------------------------'
      Write(6,*) ' '


! When area cloud fraction is requested, either parametrisation may
! be on depending on the cloud scheme

      If (L_CLD_AREA) Then

! only one of the ACF options should be set but either may be chosen

        If ((.not. L_ACF_Cusack .and. .not. L_ACF_Brooks) .or.                &
            (      L_ACF_Cusack .and.       L_ACF_Brooks)) Then

! If they are both set the same (both true or both false)
! then default to using the (old) Cusack scheme
          WRITE(6,*) ' WARNING: Area cloud frac settings inconsistent'
          WRITE(6,*) ' One and only one of the L_ACF logicals should'
          WRITE(6,*) ' be set to true when L_CLD_AREA is true'

        End If

      End If


      ! Set no of land points (either 0 or 1)
      !======================================


      land_ice_points = 0
      soil_points     = 0
      k = 0

      Do j = 1, rows
        Do i = 1, row_length
          If (land_sea_mask(i,j)) Then
            k = k + 1
            land_index(k) = (j-1)*row_length + i

            If (land_ice_mask(i,j)) Then
              land_ice_points = land_ice_points + 1
              land_ice_index(land_ice_points) = (j-1)*row_length+i
            End If

            If (soil_mask(i,j)) Then
              soil_points = soil_points + 1
              soil_index(soil_points) = (j-1)*row_length+i
            End If

          End If

          If (.NOT. land_sea_mask(i,j)) Then
            If (land_ice_mask(i,j) .OR. soil_mask(i,j)) Then

              Write (*,'(9(TR1,A52/))')                                       &
                '===================================================='        &
              , '| LAND_ICE_MASK/SOIL_MASK in inconsistent with sea |'        &
              , '| point (land_sea_mask=.FALSE.)                    |'        &
              , '| Setting following logicals to false:             |'        &
              , '|                                                  |'        &
              , '|  land_ice_mask                                   |'        &
              , '|  soil_mask                                       |'        &
              , '|                                                  |'        &
              , '===================================================='

              land_ice_mask(i,j)    = .FALSE.
              soil_mask(i,j)        = .FALSE.

            End If
          End If

        End Do
      End Do

      If (land_points /= k) Then
        Icode = 507
        Write (6,'(4(TR1,A52/))')                                             &
          '===================================================='              &
        , '| Specified total number of land points and        |'              &
        , '| land_sea_mask are inconsisent                    |'              &
        , '===================================================='
        Write(6,*) " "
        Write(6,*) "Run ABORTED"
        Write(6,*) " "

        Write(Cmessage,*) ' Check that land_points and land_sea_mask in'//    &
                          ' namelist are consistent.'

! DEPENDS ON: ereport
        Call Ereport(RoutineName,Icode,Cmessage)

      End If

      Read(10,INMOSES)

      !-----------------------------------------------------------------------
      ! Read in initial profiles and reshape to correct array sizes as
      ! specified by &SIZES namelist. Note: model_levels_nml in &CNTLSCM in
      ! namelist.scm MUST match model_levels or incorrect forcing may occur.
      !-----------------------------------------------------------------------

      Read(10,INPROF)

      theta_scm  = RESHAPE(theta   ,(/row_length,rows,    model_levels  /))
      ui_scm     = RESHAPE(ui      ,(/row_length,rows,    model_levels  /))
      vi_scm     = RESHAPE(vi      ,(/row_length,rows,    model_levels  /))
      p_in_scm   = RESHAPE(p_in    ,(/row_length,rows,    model_levels+1/))
      wi_scm     = RESHAPE(wi      ,(/row_length,rows,    model_levels+1/))
      w_advi_scm = RESHAPE(w_advi  ,(/row_length,rows,    model_levels+1/))
      qi_scm     = RESHAPE(qi      ,(/row_length,rows,wet_model_levels  /))

      If (land_points > 0) Then

        canopy_gbi_scm   = RESHAPE(canopy_gbi   ,(/land_points           /))
        smci_scm         = RESHAPE(smci         ,(/land_points           /))
        soil_type_scm    = RESHAPE(soil_type    ,(/land_points           /))
        fsmc_scm         = RESHAPE(fsmc         ,(/land_points           /))
        frac_typ_scm     = RESHAPE(frac_typ     ,(/land_points, ntype    /))
        smcli_scm        = RESHAPE(smcli        ,(/land_points, sm_levels/))
        sth_scm          = RESHAPE(sth          ,(/land_points, sm_levels/))
        t_deep_soili_scm = RESHAPE(t_deep_soili ,(/land_points, st_levels/))

        Do j=1,ntype
          Do i=1, land_points
            If (frac_typ_scm(i,j) == RMDI) Then
              Write (*,'(2(TR1,A52/), (TR1,A13,I2,A30,T53,A1/),'//            &
                       ' 2(TR1,A52/))')                                       &
                '===================================================='        &
              , '| Land tile fractions (frac_typ) has not been      |'        &
              , '| fully set. ',ntype,' fractions must be specified','|'      &
              , '| for each land point in namelist.                 |'        &
              , '===================================================='
              Stop
            End If
          End Do
        End Do

        Do i=1, land_points
          If (SUM(frac_typ_scm(i,:)) /= 1.0) Then
            Write (*,'(4(TR1,A52/))')                                         &
                '===================================================='        &
              , '| Total tile_fractions (frac_typ) must sum to 1.0  |'        &
              , '| for each land point.                             |'        &
              , '===================================================='
            Stop
          End If

          If (soil_type_scm(i) == IMDI) Then
            Write (*,'(4(TR1,A52/))')                                         &
                '===================================================='        &
              , '| Soil types (soil_type) must be specified by user |'        &
              , '| for each land point.                             |'        &
              , '===================================================='
            Stop
          End If
        End Do
      End If

      If (geoforce) Then
        READ(10,INGEOFOR)
      End If

      If (obs) Then

        !---------------------------------------------------------------------
        ! Read obs forcing data and Re-shape arrays to ensure correct forcing
        ! of the SCM.  Nfor is specified in the &CNTLSCM namelist. This MUST
        ! match the number of observational forcings in the namelist.scm file
        ! in order to avoid incorrect forcing data.
        !---------------------------------------------------------------------

        READ(10,INOBSFOR)

        t_inc_scm  = RESHAPE(t_inc ,(/row_length,rows,nfor,    model_levels/))
        u_inc_scm  = RESHAPE(u_inc ,(/row_length,rows,nfor,    model_levels/))
        v_inc_scm  = RESHAPE(v_inc ,(/row_length,rows,nfor,    model_levels/))
        w_inc_scm  = RESHAPE(w_inc ,(/row_length,rows,nfor,    model_levels/))
        q_star_scm = RESHAPE(q_star,(/row_length,rows,nfor,wet_model_levels/))

        tstar_forcing_scm =                                                   &
                     RESHAPE(tstar_forcing,(/row_length,rows,nfor/))
        flux_e_scm = RESHAPE(flux_e       ,(/row_length,rows,nfor/))
        flux_h_scm = RESHAPE(flux_h       ,(/row_length,rows,nfor/))


      End If

      CLOSE(10)


      !-----------------------------------------------------------------------
      ! Code taken from SETCONA
      ! Calculation of min_trop_level and max_trop_level
      ! NOTE: min and max_trop_level are used if climatological aerosols
      ! are chosen.
      !-----------------------------------------------------------------------
      !   The tropopause diagnosed for radiative purposes divides theta-levels
      !   considered to lie in the stratosphere from those considered to lie
      !   in the troposphere: the tropopause is therefore taken as a
      !   rho-level. This level is constrained to lie between heights of
      !   z_min_trop and z_max_trop. The level is used in the setting of the
      !   background aerosol climatology and of ozone profiles, subject to the
      !   choice of appropriate options; additionally it is used in the
      !   calculation of diagnostics defined at the tropopause.
      !
      !   Start at the second rho-level because the first is omitted from the
      !   grid seen in the physics.
      !-----------------------------------------------------------------------
      !   This code is only used if min_trop_level and max_trop_level in
      !   &RUNDATA namelist are set both set to 0.  Otherwise values in scm
      !   namelist are used.
      !-----------------------------------------------------------------------

      If ((min_trop_level == 0 .AND. max_trop_level == 0) .OR.                &
           min_trop_level <  0  .OR. max_trop_level <  0) Then

        min_trop_level=2

        Do k=1, model_levels
          r_ref_rho(k) = eta_rho(k)*z_top_of_model
        End Do

        Do ; If ((r_ref_rho(min_trop_level) >= z_min_trop) .OR.               &
                (min_trop_level == model_levels)) EXIT
             min_trop_level = min_trop_level+1
        End Do

        max_trop_level=min_trop_level
        Do ; If ( (r_ref_rho(max_trop_level) > z_max_trop) .OR.               &
                (max_trop_level == model_levels) ) EXIT
             max_trop_level = max_trop_level+1
        End Do

        max_trop_level = max_trop_level-1
      End If
!
!----------------------------------------------------------------------------

!     Derive the initial daynumber in the year and the initial time
!     from the UM data structure supplied

! DEPENDS ON: inittime
      Call INITTIME(year_init, month_init, day_init, hour_init,min_init       &
         , sec_init, dayno_init, time_initi, lcal360)

      time_init = time_initi    ! type real expected elsewhere.

!     Initial tape daynumber in the year is the same as the initial
!     daynumber

      tapeday_init = dayno_init

      PREVIOUS_TIME(1) = year_init
      PREVIOUS_TIME(2) = month_init
      PREVIOUS_TIME(3) = day_init
      PREVIOUS_TIME(4) = hour_init
      PREVIOUS_TIME(5) = min_init
      PREVIOUS_TIME(6) = sec_init
      PREVIOUS_TIME(7) = dayno_init


!---------------------------------------------------------------------
!     Set longitude to zero so that diagnostics apply to local time
!     rather than GMT and convert to radians
!---------------------------------------------------------------------

      Do j = 1, rows
        Do i = 1, row_length
          lat_rad(i,j) = pi_over_180 * lat(i,j)
          If (local_time) Then
            long_rad(i,j) = 0.0
          Else
            long_rad(i,j) = pi_over_180 * long(i,j)
          End If
        End Do
      End Do


      ! Calculate trig arrays

      Do j = 1, rows
        Do i = 1, row_length
          cos_latitude    (i,j) = cos(lat_rad(i,j))
          sec_latitude    (i,j) = 1.0/cos_latitude(i,j)
          cos_longitude   (i,j) = cos(long_rad(i,j))
          sin_latitude    (i,j) = sin(lat_rad(i,j))
          sin_longitude   (i,j) = sin(long_rad(i,j))
          FV_cos_latitude (i,j) = cos_latitude(i,j)
        End Do
      End Do

      ! Initialise various arrays

      SW_incs     = 0.0
      LW_incs     = 0.0
      dirpar_incs = 0.0
      t1_sd       = 0.0
      q1_sd       = 0.0

      area_cloud_fraction   = 0.0
      bulk_cloud_fraction   = 0.0
      cloud_fraction_liquid = 0.0
      cloud_fraction_frozen = 0.0


!
! Set h_sect - for use in choosing radiation scheme.
! Code copied from set_h_sect (readlsa2.dk)
!
      Do i = 0,nsectp
        H_SECT(i)(1:1) = '0'
        H_SECT(i)(2:3) = ATMOS_SR(i)
      End Do

! calculate radiation timestep

#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)

      Write(6,*) "WARNING: Variable ntrad is no longer used in"
      Write(6,*) "versions 3C or 3Z. The timestepping is now"
      Write(6,*) "controlled through the namelist NLSTCATM."

      If ( (ntrad/=a_sw_radstep_diag).or.                               &
     &   (ntrad/=a_sw_radstep_prog) ) Then

        ErrorStatus=715

! DEPENDS ON: ereport
        Call Ereport("S_MAIN", ErrorStatus,                             &
     &         "ERROR: ntrad is not compatible with variables"  //      &
     &         "a_sw_radstep_diag or a_sw_radstep_prog defined" //      &
     &         "in namelist NLSTCATM. ntrad is obsolete in"    //       &
     &         "versions 3C and 3Z of radiation code."
      End If

      radiation_tstep_diag = timestep * a_sw_radstep_diag
      radiation_tstep_prog = timestep * a_sw_radstep_prog
#else
      radiation_timestep = ntrad * timestep
#endif
! In SCM A_LW_SEGMENTS should be set to 1 and segment sizes for
! LW and SW should be unset (i.e. -99)
        A_LW_Segments = 1
        A_LW_Seg_size = -99
        A_SW_Seg_size = -99

! print num of substeps
        If ( L_phys2_substep )                                          &
     &    print *, 'Num of Substeps:',Num_Substeps


!---------------------------------------------------------------------
!     Set L_Triffid to false because it is not currently available
!     for use in the SCM (UM vn6.2). See Chris D Jones of the
!     Terrestrial Carbon Cycle Group, Hadley Centre if you wish to
!     use dynamical vegetation.
!---------------------------------------------------------------------
      If (L_Triffid) Then
        Write(6,*) '*************************************************'
        Write(6,*) '* Warning: you have requested 19_2A interactive *'
        Write(6,*) '* vegetation distribution, which is not         *'
        Write(6,*) '* currently available in the SCM.               *'
        Write(6,*) '* Setting L_Triffid to false                    *'
        Write(6,*) '*************************************************'
        L_Triffid = .FALSE.
      End If
!---------------------------------------------------------------------
!     Calculate values of r_theta_levels and r_rho_levels
!---------------------------------------------------------------------

! DEPENDS ON: calc_levels
      Call calc_levels                                                  &
! Input data
     &    (eta_theta, eta_rho, z_top_of_model, orog,                    &
     &     height_gen_method, boundary_layer_levels, model_levels,      &
     &     first_constant_r_rho_level, rows, row_length,                &
! Output data
     &     r_theta_levels, r_rho_levels)

!---------------------------------------------------------------------
! Calculate values of delta_lambda and delta_phi
!---------------------------------------------------------------------
! The SCM sets delta_lambda and delta_phi to 0 as a default.
! Calculate these, instead, from gridbox_area_m

       delta_lambda = 0
       delta_phi    = 0

#if !defined(SSFM)

       Do j = 1, rows
         Do i = 1, row_length
           delta_lambda = SQRT ( gridbox_area_m(i,j) /                  &
                        ( r_theta_levels(i,j,0) * r_theta_levels(i,j,0) &
                         * FV_cos_latitude(i,j) ) )
           delta_phi = delta_lambda
         End Do
       End Do

#endif
!---------------------------------------------------------------------
!     Calculate values of exner
!---------------------------------------------------------------------

! DEPENDS ON: calc_press
      Call calc_press                                                   &
! Input data
     &    (model_levels, wet_model_levels, rows, row_length, p_in_scm,  &
     &     theta_scm, qi_scm, r_theta_levels, r_rho_levels, .true.,     &
     &     .true.,                                                      &
! In/Out
     &     rho,                                                         &
! Output data
     &     exner_theta_levels, exner_rho_levels, p_theta_levels,        &
     &     rp, rp_theta, p_star)


!----------------------------------------------------------------------
! convert thetai to ti
!----------------------------------------------------------------------

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            ti_scm(i,j,k) = theta_scm(i,j,k)*exner_theta_levels(i,j,k)
          End Do
        End Do
      End Do

!---------------------------------------------------------------------
!     Set up the unit nos. for output.
!---------------------------------------------------------------------


      nout(1)=6
      If (test)          nout(2) = 22
      If (prindump_step) nout(3) = 30
      If (prindump_day)  nout(4) = 31

      If (prindump_days) Then
        If (dump_days(1) > 1) nout(5) = 32
        If (dump_days(2) > 1) nout(6) = 33
        If (dump_days(3) > 1) nout(7) = 34
        If (dump_days(4) > 1) nout(8) = 35
      End If

      If (grafdump_step) nout(9)  = 37
      If (grafdump_day)  nout(10) = 38

      If (grafdump_days) Then
        If (dump_days(1) > 1) nout(11) = 39
        If (dump_days(2) > 1) nout(12) = 40
        If (dump_days(3) > 1) nout(13) = 41
        If (dump_days(4) > 1) nout(14) = 42
      End If

      If (obs .AND. prindump_obs) Then
        Do i = 1, 5
          nout(i+14) = 42 + i
        End Do
      End If

!---------------------------------------------------------------------
!     Write out initial data for run to standard output and to all
!     the units to which diagnostics will be written.
!---------------------------------------------------------------------

! DEPENDS ON: print_initdata
      Call PRINT_INITDATA (                                                   &
        row_length, rows, land_points, model_levels, wet_model_levels,        &
        ozone_levels, nfor, boundary_layer_levels, st_levels, sm_levels,      &
        ntrop, year_init, dayno_init, ndayin, nminin, nsecin,                 &
        timestep, ntrad, a_sw_radstep_prog, a_sw_radstep_diag, lat,           &
        long, ancyc, local_time, change_clim, exname_in, runno_in,            &
        exname_out, runno_out, soil_type_scm, land_sea_mask, obs,             &
        geoforce, geoinit, stats, noforce, tapein, tapeout, smi_opt,          &
        smcli_scm, fsmc_scm, sth_scm, ug, vg, conv_mode, altdat,              &
        t_inc_scm, tstar_forcing_scm, q_star_scm, u_inc_scm, v_inc_scm,       &
        w_inc_scm, ichgf, ilscnt, flux_h_scm, flux_e_scm, ui_scm,             &
        vi_scm, wi_scm, ti_scm, qi_scm, ccai, iccbi, iccti, p_in_scm,         &
        canopy_gbi_scm, smci_scm, snodepi, t_deep_soili_scm, tstari,          &
        z0mseai, ntrad1, radcloud_fixed, cca_rad, iccb_rad, icct_rad,         &
        ccwpin_rad, layer_cloud_rad(1,1:rows,1:wet_model_levels),             &
        qcl_rad(1,1:rows,1:wet_model_levels), time_init,                      &
        ozone(1,1:rows,1:ozone_levels), nout, 19)


!---------------------------------------------------------------------
!     large scale cloud
!---------------------------------------------------------------------


      ! qcl/qcf values initialised in RUN_INIT, except for when
      ! geostropic forcing is used.

      layer_cloud = 0.0
      qcl         = 1.0e-2
      qcf         = 1.0e-2

      ! Initialise extra microphysics variables to zero as they are not
      ! set in RUN_INIT. Forcing options are currently not available for
      ! these variables.

      qcf2        = 0.0
      qrain       = 0.0
      qgraup      = 0.0
      qcf2_star   = 0.0
      qrain_star  = 0.0
      qgraup_star = 0.0


      daynumber = dayno_init
      year = 1

!=====================================================================
!     Options to set initial profiles
!=====================================================================
!     (i)   Observational large scale forcing (OBS=TRUE of
!           Namelist LOGIC)
!           Initial data is then from namelist INPROF
!     (ii)  Statistical large scale forcing (STATS=TRUE of
!           Namelist LOGIC)
!           Initial data can either be derived from climate datasets
!           using subroutine INITSTAT or set from namelist
!           INPROF (set ALTDAT=TRUE in namelist LOGIC)
!     (iii) No large-scale forcing initial data is set from namelist
!           INPROF
!     (iv)  Continuation from previous run stored on tape
!     (Set TAPEIN=TRUE in namelist LOGIC).  All other initial data
!     is overwritten
!=====================================================================

! DEPENDS ON: run_init
      Call RUN_INIT(                                                          &
        row_length, rows, model_levels, wet_model_levels, land_points,        &
        nfor, boundary_layer_levels, st_levels, sm_levels, ntiles,            &
        ntrop, nprimvars, stats, obs, prindump_obs, noforce, altdat,          &
        land_sea_mask, tapein, tapeout, l_climat_aerosol, l_use_dust,         &
        l_use_sulpc_direct, ltimer, l_ch4_lw, l_n2o_lw, l_cfc11_lw,           &
        l_cfc12_lw, l_cfc113_lw, l_hcfc22_lw, l_hfc125_lw, l_hfc134a_lw,      &
        l_o2_sw ,l_use_soot_direct, l_use_bmass_direct, l_use_ocff_direct,    &
        l_use_seasalt_direct, l_murk_rad, l_use_aod, l_use_biogenic,          &
        l_use_arclbiom, l_use_arclblck, l_use_arclsslt, l_use_arclsulp,       &
        l_use_arcldust, l_use_arclocff, l_use_arcldlta,                       &
        smi_opt, smcli_scm, fsmc_scm, sth_scm, geoforce, geoinit, ug, vg,     &
        year_init, dayno_init, lcal360, ichgf, timestep, ndayin,              &
        resdump_days, soil_type_scm, p_in_scm, smci_scm, canopy_gbi_scm,      &
        snodepi, tstari, t_deep_soili_scm, z0mseai, ui_scm ,vi_scm, wi_scm,   &
        ti_scm, qi_scm, w_advi_scm, ccai, iccbi, iccti, time_init,            &
        tapeday_init, exname_in, exname_out, runno_in, runno_out,             &
        theta_scm, theta_star, u, v, w, t, q, w_adv,                          &
        flux_h_scm, flux_e_scm, u_inc_scm, v_inc_scm, w_inc_scm,              &
        t_inc_scm, tstar_forcing_scm, q_star_scm, exner_rho_levels,           &
        exner_theta_levels, p, p_theta_levels, r_theta_levels,                &
        r_rho_levels, rho, eta_theta(1:model_levels+1),                       &
        eta_rho(1:model_levels), orog, height_gen_method,                     &
        ch_flux_h, ch_flux_e, ch_u_inc, ch_v_inc, ch_w_inc,                   &
        ch_t_inc, ch_q_star, dap1, dap2, dap3, dab1, dab2, dab3,              &
        deltap, p_star, smc, smcl, canopy_gb, snodep, tstar, tsi,             &
        t_deep_soil, sthu, sthf, frac_typ, canht, lai, catch, infil_tile,     &
        z0_tile, catch_snow,can_model, z0msea, zh, cca, n_cca_levels,         &
        iccb, icct, layer_cloud, qcf, qcl, dayno_wint, alfada, alfadb,        &
        atime, btime, lat, long, dbara, dbarb, dgrada, dgradb, pa, pb, rp,    &
        rp_theta, tbara,   tbarb, tgrada, tgradb, tsda, tsdb,                 &
        vnbara, vnbarb, vnsda, vnsdb,                                         &
        vpbara, vpbarb, wbara, wbarb, wsda, wsdb,                             &
        iv, ntab, iy, idum, iseed, resdump, rhcrit,                           &
        b_exp, hcap, hcon, satcon, sathh, v_sat, v_wilt, v_crit)




        !==========================================================
        ! Layer_cloud is initialised in RUN_INIT,
        ! area_cloud_fraction and bulk_cloud_fraction initialised
        ! to the same value as layer_cloud
        !==========================================================
        Do  i = 1, row_length
          Do j = 1, rows
            Do  k = 1, wet_model_levels
              area_cloud_fraction(i,j,k) = layer_cloud(i,j,k)
              bulk_cloud_fraction(i,j,k) = layer_cloud(i,j,k)
            End Do  ! k
          End Do  ! j
        End Do  ! i


!  Initialise all surf temperature variables to TSTARI from
!  namelist &INPROF
        Do  i = 1, row_length
          Do j = 1, rows

            If (tstar_sea(i,j)  == -999.0) Then
                tstar_sea(i,j)  = tstar(i,j)
            End If

            If (tstar_sice(i,j) == -999.0) Then
                tstar_sice(i,j) = tstar(i,j)
            End If

            If (tstar_land(i,j) == -999.0) Then
                tstar_land(i,j) = tstar(i,j)
            End If

            Do itype=1,ntype
              If (tstar_tile(i,itype) == RMDI) Then
                tstar_tile(i,itype) = tstar(i,j)
              End If
            End Do
          End Do
        End Do


      If (Obs) Then

        ! Reset tls values here as  they are changed in FORCING
        ! These _tls values store forcing values if OBS as the _incs
        ! get overwritten by the values from Atmos_Physics1
        ! They store the Atmos_physics1 increments if STATS

        Do k = 1, model_levels
          Do j2 = 1, nfor
            Do j = 1, rows
              Do i = 1, row_length
                tls(i,j,j2,k)= t_inc_scm(i,j,j2,k)
                uls(i,j,j2,k)= u_inc_scm(i,j,j2,k)
                vls(i,j,j2,k)= v_inc_scm(i,j,j2,k)
                wls(i,j,j2,k)= w_inc_scm(i,j,j2,k)
              End Do
            End Do
          End Do
        End Do

        Do  k = 1, wet_model_levels
          Do j2 = 1, nfor
            Do j = 1, rows
              Do  i = 1, row_length
                q_star_keep(i,j,j2,k) = q_star_scm(i,j,j2,k)
              End Do
            End Do
          End Do
        End Do

      End If  ! Obs

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     For geostrophic forcing
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      Do i = 1, row_length
        Do j = 1, rows
          f_coriolis(i,j) = 2.0 * omega * sin(lat(i,j) * pi_over_180)
          f3_at_u(i,j) = two_omega * sin_latitude(i,j)
        End Do
      End Do

      If (geoinit .and. geoforce) Then
        Do i = 1, row_length
          Do j = 1, rows
            modug(i,j) = sqrt(ug(i,j)*ug(i,j) + vg(i,j)*vg(i,j))
            Do k = 1, model_levels
              u(i,j,k) = ug(i,j)
              v(i,j,k) = vg(i,j)
            End Do
          End Do
        End Do

!---------------------------------------------------------------------
!       Form restart dump
!---------------------------------------------------------------------

        Do i = 1, row_length
          Do j = 1, rows
            Do k = 1, nprimvars
              resdump(i,j,k) = 0.0
            End Do
          End Do
        End Do

! DEPENDS ON: restart_dump
        Call RESTART_DUMP(                                                    &
          row_length, rows, model_levels, wet_model_levels, nprimvars,        &
          land_points, boundary_layer_levels, st_levels, sm_levels,           &
          n_cca_levels, land_sea_mask, resdump, u, v, w, t, theta_scm,        &
          q, qcl, qcf, layer_cloud, p, rho, t_deep_soil, smc, canopy_gb,      &
          snodep, tstar, zh, z0msea, cca, iccb, icct, smcl)

        Do i = 1, row_length
          Do j = 1, rows
            rccb(i,j) = Real(iccb(i,j))
            rcct(i,j) = Real(icct(i,j))
          End Do
        End Do

        stepcount = 0
        maxinc = 9999.0

        Do While (maxinc  >   0.1                                       &
     &    .AND.  stepcount  <   full_daysteps)

          stepcount = stepcount + 1
          daycount  = 1

! DEPENDS ON: timecalc
          Call TIMECALC(year_init, dayno_init, time_init, timestep      &
     &      ,time_string, lcal360, yearno, day, time_sec,               &
     &      previous_time, ihour, imin)


!         Call the pre-physics routine to set up variables

! DEPENDS ON: pre_physics
          CALL PRE_PHYSICS(                                             &
! ARGUMENTS IN
     &    row_length, rows, model_levels, wet_model_levels, qcl, qcf,   &
     &    ug, vg, f_coriolis, timestep, obs, geoforce, lcal360,         &
     &    co2start, co2rate, co2end, daycount, stepcount, ntrad1,       &
     &    ntrad, a_sw_radstep_diag, a_sw_radstep_prog,                  &
     &    l_triffid, npft,                                              &
! Arguments in/out
     &    u, v, npft_trif,                                              &
! Arguments with intent out
     &    co2_mmr, .false., .false., .false.)

!         Now call the physics routines with flags set to
!         only enable the boundary layer call.

#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
           ntrad=a_sw_radstep_prog
#endif
! DEPENDS ON: atmos_physics2
          CALL Atmos_Physics2(                                          &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, row_length, rows                  &
     &, proc_row_group, proc_col_group, at_extremity                    &
     &, n_proc, n_procx, n_procy,  neighbour, rows                      &
     &, row_length, g_datastart, me, NumCycles, CycleNo                 &

! model dimensions.
     &, row_length, rows, rows, land_points, model_levels               &
     &, nice, wet_model_levels, boundary_layer_levels                   &
     &, st_levels, sm_levels, cloud_levels, land_ice_points             &
     &, soil_points, n_cca_levels, ntiles, tr_levels                    &
     &, first_constant_r_rho_level, dim_cs1, dim_cs2                    &

!IN Substepping Information
     &, Num_Substeps, L_phys2_substep                                   &

! Model switches
     &, l_regular, l_mr_physics2, model_domain, L_dry,                  &
     &  FORMDRAG, OROG_DRAG_PARAM, LCAL360, .false., Ltimer             &
     &, L_us_blsol, BL_OPTIONS, L_cld_area, L_ACF_Cusack, L_ACF_Brooks  &
     &, L_RHCPT, L_hydrology, L_bl, L_MURK, L_MURK_ADVECT, L_MURK_SOURCE&
     &, L_BL_TRACER_MIX, L_DUST, L_CAM_DUST, L_sulpc_so2 , L_sulpc_nh3  &
     &, L_SBLeq, L_SBLco, L_sulpc_dms,  L_soot, L_biomass, L_ocff       &
     &, L_co2_interactive, l_ctile, l_co2_emits, L_use_bl_diag_term     &
     &, L_pc2, L_pc2_reset, L_eacf, L_LAMBDAM2, L_FULL_LAMBDAS, L_ukca  &
     &, L_sice_heatflux, L_OHADGEM1, L_USE_CARIOLLE                     &
     &, l_anthrop_heat_src                                              &

! model Parameters
     &, alpha_cd, Puns, Pstb, rhcrit, CO2_MMR                           &
     &, tr_vars, tr_ukca, Muw_SBL, Mwt_SBL, cloud_fraction_method       &
     &, overlap_ice_liquid, ice_fraction_method, ctt_weight             &
     &, t_weight, qsat_fixed, sub_cld, x1i, x1ic, x1r, x2r, x4r, l_psd  &
     &, ai, bi, aic, bic, lsp_ei, lsp_fi, lsp_eic, lsp_fic              &
     &, dbsdtbs_turb_0, Charnock, SeaSalinityFactor                     &

! Physical constants
     &, lc,lf,cp,two_Omega,p_zero,kappa,R,g,lapse,earth_radius,Pi       &

! in coordinate information
     &, r_rho_levels, r_theta_levels, r_rho_levels, r_rho_levels        &
     &, unscl_dry_rho                                                   &
     &, eta_theta(1:model_levels+1), eta_rho(1:model_levels)            &
     &, delta_lambda, delta_phi, dlambda_p, dphi_p                      &
     &, wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v                    &
     &, lat_rot_NP, long_rot_NP, f3_at_u                                &

! Variables required by STPH_SCV
     &, mix_v_star, mix_cl_star, mix_cf_star                            &

! in time stepping information.
     &, timestep, yearno, day, ihour, imin, time_sec, stepcount         &

! trig arrays
     &, sin_longitude, cos_longitude, sin_latitude, FV_cos_latitude, sec_latitude     &

! River routing
     &, row_length, row_length, l_rivers, xpa, xua, xva, ypa, yua, yva  &
     &, g_p_field, g_r_field, a_steps_since_riv, river_row_length       &
     &, river_rows, global_river_row_length, global_river_rows          &
     &, river_step, river_vel, river_mcoef                              &
     &, trivdir, trivseq, twatstor, inlandout_atm, l_inland,            &

! Grid-to-grid river routing
     &  r_area, slope, flowobs1, r_inext, r_jnext, r_land,              &
     &  substore, surfstore, flowin, bflowin,                           &

! diagnostic info
#include "argsts.h"
     &  STASHwork3,STASHwork5,STASHwork8,STASHwork9,STASHwork19         &
     &, STASHwork26                                                     &
!
! SCM diagnostics
     &, nSCMDpkgs, L_SCMDiags                                           &

! in data fields.
     &, theta_scm, q, qcl, qcf, rho, u, v, w, w_adv, p, p_star          &
     &, exner_rho_levels, exner_theta_levels, land_sea_mask             &
     &, p_theta_levels                                                  &

! Variables for subgrid turbulence scheme                               &
     &, visc_BL_m, fm_3d, fh_3d, .false., .false., .false., max_diff    &
     &, turb_startlev_vert, turb_endlev_vert, bl_coef_km, bl_coef_kh    &

! ancillary fields and fields needed to be kept from timestep to
! timestep

     &, land_index, land_ice_index, soil_index, canopy_gb, snodep, hcon &
     &, hcap, v_crit, v_wilt, v_sat, sthf, sthu, sil_orog_land          &
     &, ho2r2_orog, di, ice_fract, u_0, v_0, u_0, v_0, cca, iccb, icct  &
     &, cclwp, ccw, lcbase, t_deep_soil, tsi, ti_n                      &
     &, tstar, z0msea, ice_fract_n, ice_thick_n, satcon, sathh, b_exp   &
     &, smcl, t1_sd, q1_sd, zh, area_cloud_fraction, bulk_cloud_fraction&
     &, cloud_fraction_liquid, cloud_fraction_frozen, ls_rain, ls_snow  &
     &, micro_tends,photosynth_act_rad,rad_hr,surf_radflux              &
     &, soil_clay,soil_silt,soil_sand,dust_mrel1,dust_mrel2,dust_mrel3  &
     &, dust_mrel4,dust_mrel5,dust_mrel6                                &
     &, so2_high_level, so2_em, nh3_em, dms_em, soot_hilem              &
     &, soot_em, ocff_hilem, ocff_em, co2_emits, co2flux                &

! tracer fluxes - kdcorbin, 05/10
     &, tracer_flux1, tracer_flux2, tracer_flux3, tracer_flux4          &
     &, tracer_flux5, tracer_flux6, tracer_flux7, tracer_flux8          &
     &, tracer_flux9, tracer_flux10,tracer_flux11,tracer_flux12         &
     &, tracer_flux13,tracer_flux14,tracer_flux15,tracer_flux16         &
     &, tracer_flux17,tracer_flux18,tracer_flux19,tracer_flux20         &
     &, CO2EMITMASS                                                     &
! rml 1/7/13 flag for co2 flux into passive tracer
     &, L_CO2_TRACER                                                    &

! in/out
     &, theta_star                                                      &
     &, q_star_scm(1:row_length,1:rows,nfor:nfor,1:wet_model_levels)    &
     &, qcl_star, qcf_star, cf_star, cfl_star, cff_star                 &
     &, u_inc_scm(1:row_length,1:rows,nfor:nfor,1:model_levels)         &
     &, v_inc_scm(1:row_length,1:rows,nfor:nfor,1:model_levels)         &
     &, rdum4, sum_eng_fluxes, sum_moist_flux                           &

! In/Out tracer fields
     &, aerosol(1,1:rows,1:wet_model_levels), free_tracers              &
     &, ukca_tracers                                                    &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, so2, dms,  so4_aitken, so4_accu, so4_diss,  nh3, soot_new       &
     &, soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld           &
     &, ocff_new, ocff_aged, ocff_cld, co2                              &

! IN/OUT STPH_RP
     &, g0_rp, par_mezcla                                               &

! IN/OUT RIVERS
     &, tot_surf_runoff, tot_sub_runoff                                 &

! out fields
     &, rhokm, cH_term, ntml, cumulus, nbdsc, ntdsc                     &
     &, rhcpt,row_length, rows                                          &

! Additional variables for MOSES II
     &, frac_typ, frac_disturb, canht, lai                              &
     &, canopy(1:land_points,1:ntiles)                                  &
     &, catch(1:land_points,1:ntiles), catch_snow, snow_grnd, snow_tile &
     &, z0_tile, tstar_tile, infil_tile(1:land_points,1:ntiles)         &
     &, rgrain(1:land_points,1:ntiles)                                  &
     &, cs, gs, co2_dim_row, co2_dim_len, l_neg_tstar, l_snow_albedo    &
     &, .false., .false., l_trif_eq, L_Q10, asteps_since_triffid        &
     &, stepcount, phenol_period, triffid_period, can_model, g_leaf_acc &
     &, g_leaf_phen_acc, npp_ft_acc, resp_w_ft, resp_s_acc              &
     &, land_pts_trif, npft_trif, dolr, lw_down, sw_tile, fland_ctile   &
     &, tstar_land, tstar_sea, tstar_sice                               &
     &, albsoil, cos_zenith_angle, can_rad_mod, ilayers                 &

     &,radnet_tile                                                      &

!    EAK
!    IN
     &, surf_down_sw,alb_tile,l_tile_pts                                &
!     &, surf_down_sw,alb_tile,cos_zenith_angle                          &
     &, lat,long,day,time_sec,SW_DOWN                                   &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L                   &
     &, SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE             &
     &, TSOIL_TILE,T_SURF_TILE,HCONS                                    &
     &, SOIL_TYPE,VEG_TYPE                                              &
     &, ISNOW_FLG3L,total_nsteps                                        &
     &           ,FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB               &
     &           ,TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB           &
     &           ,USTAR_CAB,SURF_HTF_CAB                                &
     &, l_cable                                                         &
!sxy 
     &           ,TOT_ALB                                               &
     &           ,U_S_CAB,CH_CAB,CD_CAB                                 &
     &           ,CD,CH                                                 &
     &           ,TILE_PTS,TILE_INDEX                                   & 
     &           ,SNAGE_TILE,RTSOIL_TILE                                &
     &           ,GFLUX_TILE,SGFLUX_TILE                                &  

! Additional variables required for large-scale hydrology:
     &, L_TOP,L_PDM,FEXP,GAMTOT,TI_MEAN,TI_SIG, FSAT, FWETL, ZW         &
     &, STHZW,A_FSAT,C_FSAT,A_FWET,C_FWET,L_SOIL_SAT_DOWN               &
! Cariolle ozone
     &, OZONE_TRACER                                                    &

! Prescribed surface forcing and roughness lengths
     &, L_flux_bc                                                       &
     &, flux_e_scm(1,1:rows,nfor:nfor),flux_h_scm(1,1:rows,nfor:nfor)   &
     &, L_spec_z0,z0m_scm,z0h_scm                                       &

! Additional variables for SCM diagnostics
     &, nfor, ntrop, time_string, daycount, layer_cloud                 &
     &, tstar, nconvars, ntrad, ntrad1, conv_mode                       &

! error information
     &, Error_code  )


!         Calculate max increment and copy winds for safe-keeping
          maxinc = 0.0
          Do i = 1, row_length
            Do j = 1, rows
              Do k = 1, model_levels
                maxinc = max(maxinc,(u_inc_scm(i,j,nfor,k)**2 +         &
     &                               v_inc_scm(i,j,nfor,k)**2))
                ui_scm(i,j,k) = u(i,j,k)
                vi_scm(i,j,k) = v(i,j,k)
              End Do
              maxinc = sqrt(maxinc)/(f_coriolis(i,j)*                   &
     &                                  timestep*modug(i,j))
            End Do
          End Do


!-----------------------------------------------------------------------------
!         Copy initial data back from DUMP.
!-----------------------------------------------------------------------------

! DEPENDS ON: dumpinit
          Call DUMPINIT(                                                      &
            row_length, rows, nprimvars, land_points, model_levels,           &
            wet_model_levels, boundary_layer_levels, st_levels, sm_levels,    &
            ntrop, n_cca_levels, land_sea_mask, resdump, u, v, w, t,          &
            theta_scm, q, qcl, qcf, layer_cloud, p, rho, t_deep_soil, smc,    &
            canopy_gb, snodep, tstar, zh, z0msea, cca, rccb, rcct, smcl)

          Do i = 1, row_length
            Do j = 1, rows
              iccb(i,j) = int(rccb(i,j))
              icct(i,j) = int(rcct(i,j))

!             Copy saved U,V back

              Do k = 1, model_levels
                u(i,j,k) = ui_scm(i,j,k)
                v(i,j,k) = vi_scm(i,j,k)
              End Do
            End Do
          End Do
        End Do                   ! maxinc and stepcount < daysteps

        Do i = 1, row_length
          Do j = 1, rows
            Do k = 1, nprimvars
              resdump(i,j,k) = 0.0
            End Do
          End Do
        End Do

        Write (6,*) "Geostrophic wind initialised."
        Write (6,*) "max relative wind change at end=",maxinc,          &
     &    " after ",stepcount," steps"
        daynumber = dayno_init
        year = 1

      End If                     ! Geostrophic forcing initialising.

!---------------------------------------------------------------------
!     Initialise the output diagnostic system
!---------------------------------------------------------------------
! DEPENDS ON: setup_diags
      call setup_diags(row_length,rows,model_levels,            & ! IN
     &     wet_model_levels,boundary_layer_levels,sm_levels,    & ! IN
     &     st_levels,land_points,ntiles,n_vis_thresh,           & ! IN
     &     cloud_levels,                                        & ! IN
     &     total_nsteps,timestep,full_daysteps,ntrad,           & ! IN
     &     a_sw_radstep_prog,a_sw_radstep_diag,ntrad1,          & ! IN
     &     daycount,stepcount,Num_Substeps,main_diag_switch,    & ! IN
     &     netcdf_chunksize,nSCMDpkgs,L_SCMDiags,               & ! IN
     &     SCMop)                                                 ! OUT

! If PC2 is off then L_SCMDiags(SCMDiag_pc2) must be false
      If (.NOT. L_PC2) Then
        If (L_SCMDiags(SCMDiag_pc2)) Then
          write(6,*) '***********************************'
          write(6,*) '* Warning: you have requested PC2 *'
          write(6,*) '* diagnostics but PC2 is off,     *'
          write(6,*) '* resetting L_SCMDiag_PC2 logical *'
          write(6,*) '***********************************'
        End If
        L_SCMDiags(SCMDiag_pc2) = .FALSE.
      End If

! The availability of Surface based diagnostics packages is determined
! by the surface type.
!
! Surface Package     (SCMDiags_surf) - always available
! Land Points Package (SCMDiags_land) - Only if land_sea_mask TRUE
! Sea Points Package  (SCMDiags_sea)  - Only if land_sea_mask FALSE
!
! This only works in the SCM because rows, row_length are size 1 here

      Do j = 1, rows
        Do i = 1, row_length
          If (land_sea_mask(i,j)) Then

!           ! land point
            If (L_SCMDiags(SCMDiag_sea)) Then
              write(6,*) '********************************************'
              write(6,*) '* Warning: you have requested sea          *'
              write(6,*) '* diagnostics but this is a land point,    *'
              write(6,*) '* resetting L_SCMDiag_sea logical to false *'
              write(6,*) '********************************************'
            End If

            L_SCMDiags(SCMDiag_sea) = .FALSE.

          Else

!           ! sea point
            If (L_SCMDiags(SCMDiag_land)) Then
              write(6,*)'*********************************************'
              write(6,*)'* Warning: you have requested land          *'
              write(6,*)'* diagnostics but this is a sea point,      *'
              write(6,*)'* resetting L_SCMDiag_land logical to false *'
              write(6,*)'*********************************************'
            End If

            L_SCMDiags(SCMDiag_land) = .FALSE.

          End If ! land_sea_mask
        End Do
      End Do

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     LOOP OVER DAYS
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      ! Timestepping proper is about to begin, switch the
      ! diagnostic system on (as long as at least one stream is open)
      if (main_diag_switch /= 0.and.                                    &
     &        any(SCMop%strm%switch /= 0)) Then
            ! any() is a Fortran90 intrinsic funtion
         SCMop%on=.true.
      End If


! EAK
!      print *,'2scm_mainPREVIOUS_TIME',PREVIOUS_TIME,timestep,daysteps
!      print *,'2scm_main time0',year_init, dayno_init, time_init
!      print *,'2scm_main time',stepcount,total_nsteps,nstepsin
!      print *,'2scm_main time1',daycount,daysteps ,full_daysteps
!      print *,'2scm_main time2',stepcount,ihour,imin,day,time_sec

      SNOW_DEPTH3L = 0.0
      SNOW_MASS3L = 0.0
      SNOW_COND = SNOW_HCON
      SNOW_TMP3L = 273.16
      SNOW_RHO3L = 120.0
      SNOW_RHO1L = 120.0
      STHU_TILE = 0.0
      T_SURF_TILE = 0.0
      HCONS = 0.0
      ISNOW_FLG3L = 0
      ! CABLE tiled soil temperature and moisture initialised the 
      ! same on all tiles
      Do itype=1,ntype
         tsoil_tile(:,itype,:) = t_deep_soil(:,:)
         smcl_tile(:,itype,:) = smcl(:,:)
         sthf_tile(:,itype,:) = sthf(:,:)
      End Do

!     Cable diagnostics for output
      FTL_TILE_CAB = 0.
      FTL_CAB = 0.
      LE_TILE_CAB = 0.
      LE_CAB = 0.
      TSTAR_TILE_CAB = 273.16
      TSTAR_CAB = 273.16
      TSOIL_CAB = 273.16
      SMCL_CAB = 100.0
      USTAR_CAB = 100.0
      SURF_HTF_CAB = 0.0

      istep_cur = 0

      Do daycount = 1, ndayin+1

        istep_cur = istep_cur + 1

!---------------------------------------------------------------------
!       Reset sinusoidal distribution every 10 days if climate
!         stats required
!---------------------------------------------------------------------

        If (stats .and. (daycount  ==  1                                &
     &    .or. (ancyc .and. mod(daycount, change_clim)  ==  0))) Then
! DEPENDS ON: statday
          Call STATDAY(                                                 &
!         IN leading dimensions of arrays
     &      row_length,rows, model_levels, wet_model_levels,ntrop,      &
!
     &      atime,btime,dayno_wint,deltan,daycount,                     &
     &      tbara,tbarb,tsda,tsdb,dbara,dbarb,vnbara,vnbarb,            &
     &      vnsda,vnsdb,vpbara,vpbarb,wbara,wbarb,wsda,wsdb,            &
     &      alfada,alfadb,pa,pb,p,tgrada,                               &
     &      tgradb,dgrada,dgradb,cort,cord,corvn,corw,                  &
     &      tdash,ddash,ctbar,ctsd,at,cdbar,cdsd,ad,                    &
     &      cvnbar,cvnsd,avn,cwbar,cwsd,aw,                             &
     &      tbar,tsd,dbar,dsd,                                          &
     &      vnbar,vnsd,vpbar,wbar,wsd)

!---------------------------------------------------------------------
!     Calculate values of exner
!---------------------------------------------------------------------

! DEPENDS ON: calc_press
        Call Calc_press                                                 &
! Input data
     &    (model_levels, wet_model_levels, rows, row_length, p,         &
     &     theta_scm, q, r_theta_levels, r_rho_levels, .true., .false., &
! In/Out
     &     rho,                                                         &
! Output data
     &     exner_theta_levels, exner_rho_levels, p_theta_levels,        &
     &     rp, rp_theta, p_star)

!---------------------------------------------------------------------
!         Initialise PX and PY arrays for calculation of
!           vertical fluxes later
!---------------------------------------------------------------------

          Do i = 1, row_length
            Do j = 1, rows
              Do k = 1, ntrop
                px(i,j,k) = 1. / alog(p(i,j,k+1)/ p(i,j,k))
              End Do

              Do k = 1, ntrop-1
                py(i,j,k) = 1. / alog(p(i,j,k+2)/ p(i,j,k))
              End Do
            End Do
          End Do
        End If                   ! stats

!=====================================================================
!       Options for forcing
!---------------------------------------------------------------------
!
!       Observational forcing : use observed values of T,Q,U,V
!       to calculate their increments due to large scale effects.
! OR
!       Statistical forcing : take random sample from Normal
!       (Gaussian) distribution with mean and sd climlogical
!       average to calculate increments to T and Q due to large scale
!       effects.
!=====================================================================
!
!       Loop over timesteps
!
!
!       If it is the last day in the run and a full day is not
!       required, loop over the number of timesteps required
!       otherwise Do the full number of timesteps in a day.

        If (daycount  ==  ndayin+1                                      &
     &    .and. nstepsin  /=  full_daysteps) Then
          daysteps = nstepsin
        Else
          daysteps = full_daysteps
        End If

        Do stepcount = 1, daysteps

          ! Update local_timestep_count
          scm_timestep_count = stepcount + ((daycount-1)*daysteps)

!---------------------------------------------------------------------
!         VARIABLE VALUES: q         is vapour         at timelevel n
!                          qcl       is liquid         at timelevel n
!                          t         is temperature    at timelevel n
!                          theta_scm is potential temp at timelevel n
!---------------------------------------------------------------------


          If (main_diag_switch /= 0) Then
!
!-----------------------------------------------------------------------
!           SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------
            If (L_SCMDiags(SCMDiag_PC2)) Then

! DEPENDS ON: scmoutput
              Call SCMoutput(q,                                         &
                   'q_timen','Vapour at start of timestep','kg/kg',     &
                   t_inst,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(qcl,                                       &
                   'qcl_timen','Liquid at start of timestep','kg/kg',   &
                   t_inst,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(t,                                         &
                   't_timen','Temperature at start of timestep','K',    &
                   t_inst,d_all,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(theta_scm,                                 &
                   'th_timen',                                          &
                   'Potential temperature at start of timestep','K',    &
                   t_inst,d_all,default_streams,'',RoutineName)

            End If ! L_SCMDiags(SCMDiag_PC2)

          End If ! main_diag_switch /= 0

          !
          ! Reset the wind increments to zero -
          ! The other increments are reset in physics1
          !

          Do k = 1, model_levels
            Do j2 = 1, nfor
              Do j = 1, rows
                Do i = 1, row_length
                  u_inc_scm(i,j,j2,k)= 0.0
                  v_inc_scm(i,j,j2,k)= 0.0
                  w_inc_scm(i,j,j2,k)= 0.0
                End Do
              End Do
            End Do
          End Do

!---------------------------------------------------------------------
!         Calculate the year (in run) and actual time and day number
!         for labelling  of the diagnostics only.
!---------------------------------------------------------------------

! DEPENDS ON: timecalc
          Call TIMECALC(year_init, dayno_init, time_init, timestep      &
     &      ,time_string, lcal360, yearno, day, time_sec,               &
     &       previous_time, ihour,imin)

!         If there is no annual cycle, the year and day numbers
!         will be the init ones.
          If (.not. ancyc) Then
            year = 1
            day = dayno_init
          End If

          !-------------------------------------------------------------
          ! Convert temperature to potential temperature
          !         t_inc_scm   to theta_star
          !-------------------------------------------------------------

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                theta_scm(i,j,k)  = t(i,j,k)/exner_theta_levels(i,j,k)
                theta_star(i,j,k) = t_inc_scm(i,j,nfor,k)/              &
     &                              exner_theta_levels(i,j,k)
              End Do
            End Do
          End Do
!
          If (L_pc2) Then
!
!--------------------------------------------------------------------
!  Calculate initial relative humidity w.r.t. Liquid temperature TL
!  Is used by PC2 routines at end of timestep. Jeremy Price Feb 2005.
!--------------------------------------------------------------------
            Do k=1,wet_model_levels
              Do j=1,rows
                Do i=1,row_length
                  TL(i,j)=T(i,j,k)-(LC/CP)*QCL(i,j,k)
                End Do ! i
              End Do ! j
! DEPENDS ON: qsat_wat
              CALL QSAT_WAT(qsl_tl,tl,p_theta_levels(1,1,k),            &
     &                      row_length*rows)
              Do j=1,rows
                Do i=1,row_length
                  RHTS(i,j,k)=(Q(i,j,k)+QCL(i,j,k))/QSL_TL(i,j)
                  tlts(i,j,k)=TL(i,j)
                  ptts(i,j,k)=p_theta_levels(i,j,k)
                  qtts(i,j,k)=Q(i,j,k)+QCL(i,j,k)
                End Do ! i
              End Do ! j
            End Do ! k
!
          End If  ! L_pc2
!
!
!---------------------------------------------------------------------
! Save initial fields for later to calculate total increments
!---------------------------------------------------------------------
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                t_start(i,j,k)   = t(i,j,k)
                u_start(i,j,k)   = u(i,j,k)
                v_start(i,j,k)   = v(i,j,k)
                w_start(i,j,k)   = w(i,j,k)
              End Do
            End Do
          End Do
          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                q_start(i,j,k)   = q(i,j,k)
                qcl_start(i,j,k) = qcl(i,j,k)
                qcf_start(i,j,k) = qcf(i,j,k)
              End Do
            End Do
          End Do

!---------------------------------------------------------------------
!       Call physics one before advection step
!---------------------------------------------------------------------

! DEPENDS ON: pre_physics
          Call PRE_PHYSICS(                                             &
! arguments in
     &    row_length, rows, model_levels, wet_model_levels, qcl, qcf,   &
     &    ug, vg, f_coriolis, timestep, obs, geoforce, lcal360,         &
     &    co2start, co2rate, co2end, daycount, stepcount, ntrad1,       &
     &    ntrad, a_sw_radstep_diag, a_sw_radstep_prog,                  &
     &    l_triffid, npft,                                              &
! arguments in/out
     &    u, v, npft_trif,                                              &
! arguments with intent out
     &    co2_mmr, L_rad_step, L_rad_step_prog, L_rad_step_diag)

#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
!
! Set flags for forcing to false
!
          L_forcing=.false.
#endif

!------------------------------------------------------------------|
! UM5.x timestepping stores all wind increments in u_inc_scm and   |
! v_inc_scm. These are then added to u,v in atmos_physics2.        |
!------------------------------------------------------------------|

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length

               !---------------------------------------------------|
               ! Calculate increments from geostrophic forcing     |
               ! These are then added to u_inc_scm and v_inc_scm   |
               ! before atmos_physics2 as atmos_physics1 sets      |
               ! u_inc_scm and v_inc_scm to zero.                  |
               !---------------------------------------------------|

                uinc_geo(i,j,k) = u(i,j,k)-u_start(i,j,k)
                vinc_geo(i,j,k) = v(i,j,k)-v_start(i,j,k)

               !--------------------------------|
               ! Reset u,v back to time-level n |
               !--------------------------------|

                u(i,j,k) = u_start(i,j,k)
                v(i,j,k) = v_start(i,j,k)

              End Do  ! i
            End Do  ! j
          End Do  ! k

          ! ----------------------------------
          ! Convert to mixing ratios if needed
          ! ----------------------------------
          If (l_mr_physics1) Then

            Write(6,*) 'Convert to mixing ratios'

            ! Convert to mixing ratios
! DEPENDS ON: q_to_mix
            Call q_to_mix (row_length, rows, wet_model_levels,          &
     &           halo_i, halo_j,                                        &
     &           q, qcl, qcf,                                           &
     &           qcf2, qrain, qgraup,                                   &
     &           L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                 &
     &           mix_v, mix_cl, mix_cf,                                 &
     &           mix_cf2, mix_rain, mix_graup                           &
     &          )

            ! Now place mixing ratio values back into d1

            Do k = 1, wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length

                  q(i,j,k)     = mix_v(i,j,k)     ! Vapour
                  qcl(i,j,k)   = mix_cl(i,j,k)    ! Liquid
                  qcf(i,j,k)   = mix_cf(i,j,k)    ! Ice
                  qcf2(i,j,k)  = mix_cf2(i,j,k)   ! Ice 2
                  qrain(i,j,k) = mix_rain(i,j,k)  ! Rain
                  qgraup(i,j,k)= mix_graup(i,j,k) ! Graupel

                End Do
              End Do
            End Do

          End If  ! l_mr_physics1

!=======================================================================

! DEPENDS ON: atmos_physics1
          Call Atmos_Physics1 (                                         &

! Parallel variables
     &  halo_i,halo_j,off_x,off_y,row_length, rows,  proc_row_group     &
     &, proc_col_group, at_extremity, n_proc, n_procx, n_procy          &
     &, neighbour, rows, row_length, g_datastart, me                    &

! model dimensions.
     &, row_length, rows, rows, land_points, model_levels               &
     &, wet_model_levels, boundary_layer_levels, st_levels, sm_levels   &
     &, Ozone_levels, cloud_levels, land_ice_points, soil_points        &
     &, n_cca_levels, ntiles, salt_dim1, salt_dim2, salt_dim3           &
     &, tr_levels, tr_vars, co2_dim_len, co2_dim_row, co2_dim_lev       &
     &, n_arcl_species, n_arcl_compnts, i_arcl_compnts                  &

! Model switches
     &, model_domain, l_regular, L_SEC_VAR, L_EqT                       &
     &, L_rad_step, L_Rad_Step_diag,L_Rad_Step_prog, L_Forcing          &
     &, L_Timestep, L_Radiance, L_Wenyi, LCAL360, L_microphy            &
     &, L_emcorr, L_climat_aerosol, Ltimer, L_gwd, l_use_ussp           &
     &, l_taus_scale, l_fix_gwsatn, l_gwd_40km, l_ussp_opaque           &
     &, sat_scheme, l_use_clearrh, l_ssice_albedo, l_rhcpt, l_murk      &
     &, L_MURK_SOURCE                                                   &
     &, L_MURK_BDRY, L_MURK_RAD, L_DUST, l_sulpc_so2                    &
     &, l_sulpc_nh3, l_soot,  L_biomass, l_ocff, l_co2_interactive      &
     &, l_co2_radiation                                                 &

! Setting lflux_reset to false but, if using energy correction, will
! need to be worked out every timestep
     &, LFLUX_RESET, L_CLIM_AERO_HGT, L_HadGEM1_Clim_Aero, L_USE_DUST   &
     &, L_use_sulphate_autoconv, L_auto_debias, L_use_seasalt_autoconv  &
     &, L_use_seasalt_indirect, L_use_seasalt_direct                    &
     &, L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a              &
     &, L_snow_albedo, l_ctile, L_radiation, L_rain                     &
     &, L_INHOM_CLOUD, l_use_biogenic, L_use_sulpc_direct               &
     &, L_use_soot_direct, L_use_soot_indirect, L_use_soot_autoconv     &
     &, L_use_bmass_direct, L_use_bmass_indirect, L_use_bmass_autoconv  &
     &, L_use_ocff_direct, L_use_ocff_indirect, L_use_ocff_autoconv     &
     &, L_use_sulpc_indirect_SW, L_use_sulpc_indirect_LW, L_pc2         &
     &, L_eacf, L_mr_physics1, l_cry_agg_dep, l_droplet_settle          &
     &, L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_use_methox             &
     &, L_rad_deg, .false., l_use_stochem_ch4, L_it_melting, L_ukca     &
     &, L_USE_ARCL, L_USE_SPEC_SEA, L_MOD_BARKER_ALBEDO, L_MOD_K_FLUX   &

! model Parameters
     &, rhcrit, cw_sea, cw_land                                         &
     &, A_SW_segments, A_SW_seg_size, A_LW_segments, A_LW_seg_size, 0, 0&
     &, 0, 0, 0, 0, aero_bl_levels                                      &
     &, INHOM_CLOUD_SW, INHOM_CLOUD_LW, DP_CORR_STRAT, DP_CORR_CONV     &
     &, CO2_MMR, alpham, alphac, alphab, dtice                          &
     &, dt_bare,dalb_bare_wet,pen_rad_frac,SW_beta                      &
     &, min_trop_level, max_trop_level, 0.0, 0.0, 0.0                   &
     &, O2MMR, N2Ommr, CH4mmr, C11mmr, C12mmr                           &
     &, C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR                       &
! Pass down grgas_addr to atmos_physics1
     &, ngrgas, grgas_addr                                              &
     &, cloud_fraction_method,overlap_ice_liquid                        &
     &, ice_fraction_method,ctt_weight,t_weight,qsat_fixed,sub_cld      &
     &, dbsdtbs_turb_0                                                  &
     &, L_seq_mcr,L_autoc_3b,L_autolim_3b,L_autoconv_murk               &
     &, ec_auto,N_drop_land                                             &
     &, N_drop_sea,N_drop_land_cr,N_drop_sea_cr,Ntot_land, Ntot_sea     &
     &, x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic                        &
     &, lsp_ei,lsp_fi,lsp_eic,lsp_fic                                   &
! parameter for stochastic physics random parameters2
     &, m_ci                                                            &

! Physical constants
     &, lc,lf,cp,two_Omega,p_zero,kappa,R,g,lapse,earth_radius,Pi       &

! in coordinate information
     &, r_rho_levels, r_theta_levels, r_rho_levels, r_rho_levels        &
     &, eta_theta(1:model_levels+1), eta_rho(1:model_levels)            &
! atmos_physics1 call - dummy argument to replace VarRes coordinate info
!  order is delta_lambda, delta_phi, lat_rot_NP, long_rot_NP
     &, 0, 0, 0, 0                                                      &

! in time stepping information.
     &, timestep,radiation_timestep,radiation_tstep_diag                &
     &, radiation_tstep_prog, yearno, day, ihour                        &
     &, imin, idum0, stepcount, PREVIOUS_TIME                           &
     &, istep_cur                                                       &

! trig arrays
     &, sin_longitude, cos_longitude, FV_cos_latitude, sin_latitude     &


! grid-dependent arrays
     &, f_coriolis, long_rad, lat_rad,                                  &

! diagnostic info
#include "argsts.h"
     &  0, 0, 0, 0, 0                                                   &
!
! SCM diagnostics
     &, nSCMDpkgs, L_SCMDiags                                           &

! in data fields.
     &, theta_scm, q, qcl, qcf, qcf2, qrain, qgraup, rho,  u, v, p      &
     &, p_star, exner_rho_levels, exner_theta_levels, land_sea_mask     &
     &, p_theta_levels, fland_ctile, frac_dummy                         &

! ancillary fields and fields needed to be kept from timestep to
! timestep

     &, land_index, rgrain(1:land_points,1:ntiles), soot, ntml, cumulus &
     &, ice_fract, cca, iccb                                            &
     &, icct, cclwp, ccw, lcbase, tstar, tstar_land                     &
     &, tstar_sea, tstar_sice, sice_alb, land_alb, snodep, snodep       &
     &, ozone(1:row_length,1:rows,1:ozone_levels), SW_incs, LW_incs     &
     &, dirpar_incs, O3_trop_level, O3_trop_height, T_trop_level        &
     &, T_trop_height, zh, 0.0, 0.0, 0.0, 0.0, area_cloud_fraction      &
     &, bulk_cloud_fraction, cloud_fraction_liquid                      &
     &, cloud_fraction_frozen, aerosol_em, arcl, albsoil, lai           &
     &, snow_tile, frac_typ, tstar_tile, z0_tile                        &
     &, dOLR_rts, LW_down, SW_tile_rts, ES_SPACE_INTERP_dummy, .true.   &
     &, rdum4                                                           &
                ! this is for ch4_stochem
     &, cos_zenith_angle,can_rad_mod                                    &

! in/out
     &, theta_star                                                      &
     &, q_star_scm(1:row_length,1:rows,nfor:nfor,1:wet_model_levels)    &
     &, qcl_star, qcf_star, qcf2_star, qrain_star, qgraup_star          &
     &, cf_star, cfl_star, cff_star                                     &
     &, u_inc_scm(1:row_length,1:rows,nfor:nfor,1:model_levels)         &
     &, v_inc_scm(1:row_length,1:rows,nfor:nfor,1:model_levels)         &
     &, energy_correction, sum_eng_fluxes, sum_moist_flux               &
     &, aerosol(1:row_length,1:rows,1:wet_model_levels)                 &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, so2(1:row_length,1:rows,1:wet_model_levels)                     &
     &, so4_aitken(1:row_length,1:rows,1:wet_model_levels)              &
     &, so4_accu(1:row_length,1:rows,1:wet_model_levels)                &
     &, so4_diss(1:row_length,1:rows,1:wet_model_levels)                &
     &, nh3(1:row_length,1:rows,1:model_levels)                         &
     &, soot_new, soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld &
     &, ocff_new, ocff_aged, ocff_cld                                   &
     &, co2, free_tracers, biogenic, ASTEPS_SINCE_TRIFFID               &

!    EAK
!    IN
     &, surf_down_sw,alb_tile                                           &
     !sxy 
!     &, day,L_TILE_PTS,TILE_PTS,SM_LEVELS,TILE_INDEX                    &
     &, day,TILE_PTS,SM_LEVELS,TILE_INDEX                               &
     &, SNOW_TMP3L,SNOW_RHO1L,TSOIL_TILE,ISNOW_FLG3L,LAND_ALBEDO        &
     &, l_cable                                                         & 

! out fields
     &, ls_rain, ls_snow, micro_tends, unscl_dry_rho                    &
     &, photosynth_act_rad, rad_hr, surf_radflux, dOLR, SW_tile         &

! Section information
     &,     maxsects,h_sect                                             &

! error information
     &, Error_code  )

          !-------------------------------
          ! Convert to specific humidities
          !-------------------------------

          If (l_mr_physics1) Then
            write(6,*) 'In l_mr_physics1 branch of reconversion'

            Do k = 1, wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length

                  ! Copy q and qstar variables to mix variables
                  mix_v_star(i,j,k)  = q_star_scm(i,j,nfor,k)
                                                          ! Vapour
                  mix_cl_star(i,j,k) = qcl_star(i,j,k)    ! Liquid
                  mix_cf_star(i,j,k) = qcf_star(i,j,k)    ! Ice
                  mix_cf2_star(i,j,k) = qcf2_star(i,j,k)  ! Ice2
                  mix_rain_star(i,j,k) = qrain_star(i,j,k)! Rain
                  mix_graup_star(i,j,k) = qgraup_star(i,j,k)
                                                          ! Graupel

                End Do
              End Do
            End Do

            ! Convert mixing ratios back to specific humidities
! DEPENDS ON: mix_to_q
            Call mix_to_q (row_length, rows, wet_model_levels,          &
     &                   halo_i, halo_j,                                &
     &                   mix_v, mix_cl, mix_cf,                         &
     &                   mix_cf2, mix_rain, mix_graup,                  &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   q, qcl, qcf,                                   &
     &                   qcf2, qrain, qgraup                            &
     &                   )

            ! Convert mixing ratio increments (mix_star) back
            ! to specific humidity increments (q_star_scm)
! DEPENDS ON: calc_q_star
            Call calc_q_star (row_length, rows, wet_model_levels,       &
     &                   halo_i, halo_j, off_x, off_y,                  &
     &                   mix_v, mix_cl, mix_cf,                         &
     &                   mix_cf2, mix_rain, mix_graup,                  &
     &                   mix_v_star, mix_cl_star, mix_cf_star,          &
     &                   mix_cf2_star, mix_rain_star, mix_graup_star,   &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   q, qcl, qcf,                                   &
     &                   qcf2, qrain, qgraup,                           &
     &  q_star_scm(1:row_length,1:rows,nfor:nfor,1:wet_model_levels),   &
     &                   qcl_star, qcf_star,                            &
     &                   qcf2_star, qrain_star, qgraup_star             &
     &                   )

          End If  ! l_mr_physics1


!------------------------------------------------------------------|
! Add on increments from geostrophic wind forcing.                 |
!------------------------------------------------------------------|
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                u_inc_scm(i,j,nfor,k) = u_inc_scm(i,j,nfor,k)           &
     &                                + uinc_geo(i,j,k)
                v_inc_scm(i,j,nfor,k) = v_inc_scm(i,j,nfor,k)           &
     &                                + vinc_geo(i,j,k)
              End Do  ! i
            End Do  ! j
          End Do  ! k

!---------------------------------------------------------------------
!       Convert theta_star to t_inc_scm for call to forcing
!---------------------------------------------------------------------

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                t_inc_scm(i,j,nfor,k) = theta_star(i,j,k)               &
     &                                * exner_theta_levels(i,j,k)
              End Do
            End Do
          End Do

! set qcl_inc and qcf_inc to increments

          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                qcl_inc(i,j,nfor,k) = qcl_star(i,j,k)
                qcf_inc(i,j,nfor,k) = qcf_star(i,j,k)
              End Do
            End Do
          End Do
!
!---------------------------------------------------------------------
!         VARIABLE VALUES:
!         qcl_inc and qcl_star = liq. water incs      (atmos_physics1)
!         q_star_scm           = vapour incs          (atmos_physics1)
!         t_inc_scm            = temp. incs           (atmos_physics1)
!         theta_star           = pot. temp.incs       (atmos_physics1)
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!         PC2 section. Store increments from atmos_physics1. We can
!            then work out the forcing by subtraction of these
!            _earliest values from the values returned after forcing.
!---------------------------------------------------------------------
!
          If (L_pc2) Then
            Do k = 1, wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  q_earliest  (i,j,k) = q_star_scm(i,j,nfor,k)
                  qcl_earliest(i,j,k) = qcl_inc   (i,j,nfor,k)
                  t_earliest  (i,j,k) = t_inc_scm (i,j,nfor,k)
                End Do ! i
              End Do ! j
            End Do ! k
          End If  ! L_pc2
!
          If (main_diag_switch /= 0) Then
!
!-----------------------------------------------------------------------
!           SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------
            If (L_SCMDiags(SCMDiag_PC2)) Then

! DEPENDS ON: scmoutput
              Call SCMoutput(q_earliest,                                &
                  'dq_earliest','q_star vapour incs from atmos_phys1',  &
                   'kg/kg',t_inst,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(qcl_earliest,                              &
                   'dqcl_earliest',                                     &
                   'qcl_inc liq water incs atmos_phys1',                &
                   'kg/kg',t_inst,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(t_earliest,                                &
                   'dt_earliest','t_inc temp incs from atmos_phys1',    &
                   'K',t_inst,d_all,default_streams,'',RoutineName)

            End If ! L_SCMDiags(SCMDiag_PC2)

          End If ! main_diag_switch /= 0

!---------------------------------------------------------------------
!         Include any forcing required
!---------------------------------------------------------------------

          If (stats .or. obs) Then
! DEPENDS ON: forcing
            Call FORCING(                                               &
     &        row_length, rows, model_levels, wet_model_levels,         &
     &        nfor, boundary_layer_levels, st_levels, sm_levels, ntrop, &
     &        sec_day,stats, obs, prindump_obs, prinstat,               &
     &        dayno_wint, daycount, daysteps, stepcount,                &
     &        timestep, ichgf, ad, at, avn, aw, cdbar, cdsd, ctbar,     &
     &        ctsd, cvnbar, cvnsd, cwbar, cwsd, dbar, dsd, ddash,       &
     &        deltan, p, rp, px, py, tbar, tdash, tsd, vnbar, vnsd,     &
     &        vpbar, wbar, wsd, t, q, u, v, w,                          &
     &        nSCMDpkgs, L_SCMDiags,relaxT_eurocs,relaxq_eurocs,        &
     &        relaxuv_eurocs,tau_T_eurocs, tau_q_eurocs, tau_uv_eurocs, &
     &        qr, tr, vnr, vpr, wr, qcl, qcf,                           &
     &        flux_h_scm, flux_e_scm, t_inc_scm, tstar_forcing_scm,     &
     &        tstar, q_star_scm, u_inc_scm, v_inc_scm, w_inc_scm,       &
     &        r_theta_levels, exner_theta_levels,                       &

! Add variables for qcf and qcl forcing
     &        qcl_inc, qcf_inc,                                         &
     &        ch_flux_h, ch_flux_e, ch_t_inc, ch_q_star, ch_u_inc,      &
     &        ch_v_inc, ch_w_inc, dap1, dab1,                           &
     &        ui_scm, vi_scm, ti_scm, qi_scm, ilscnt,                   &
     &        rhokh, factor_rhokh, iv, ntab, iy, idum, q_star_keep,     &
     &        uls, vls, wls, tls, l_windrlx, l_vertadv, tau_rlx         &
     &        )
          End If
!
!---------------------------------------------------------------------
!         VARIABLE VALUES:
!         qcl_inc    = liquid water increments        (atmos_physics1)
!                    + forcing increments.
!         t_inc_scm  = temperature increments         (atmos_physics1)
!                    + forcing increments
!         q_star_scm = vapour increments              (atmos_physics1)
!                    + forcing increments
!         wls(nfor)  = the current vertical velocity
!---------------------------------------------------------------------
!         Update model vertical velocities to be consistent with
!         subsidence forcing
!---------------------------------------------------------------------

! Want latest w whether vertical advection on or off.
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                w(i,j,k)     = wls(i,j,nfor,k)
                w_adv(i,j,k) = w(i,j,k)
              End Do ! i
            End Do ! j
          End Do ! k

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!         This is a PC2 section. Calculate the forcing of vapour,
!         liquid, temperature and pressure.
!---------------------------------------------------------------------
!
          If (L_pc2) Then
            Do k = 1, wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
! Calculate the forcing of the variables
                  q_forcing  (i,j,k) = q_star_scm(i,j,nfor,k)           &
     &                               - q_earliest(i,j,k)
                  qcl_forcing(i,j,k) = qcl_inc(i,j,nfor,k)              &
     &                               - qcl_earliest(i,j,k)
                  t_forcing  (i,j,k) = t_inc_scm(i,j,nfor,k)            &
     &                               - t_earliest(i,j,k)
                  p_forcing  (i,j,k) = 0.0
                End Do ! i
              End Do ! j
            End Do ! k
          End If  ! L_pc2
!
!---------------------------------------------------------------------
!         Increments need to be converted to increment+value for t,q
!         qcl and qcf
!---------------------------------------------------------------------

!  Before going into physics2, t_inc_scm, q_inc, qcl_inc and qcf_inc are
!  converted to increment plus value whilst the winds stay as increments

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            t_inc_scm(i,j,nfor,k) = t_inc_scm(i,j,nfor,k) + t(i,j,k)
          End Do
        End Do
      End Do

      Do k = 1, wet_model_levels
        Do j = 1, rows
          Do i = 1, row_length
            q_star_scm(i,j,nfor,k) = q_star_scm(i,j,nfor,k) + q(i,j,k)
            qcl_star(i,j,k)        = qcl_inc(i,j,nfor,k)    + qcl(i,j,k)
            qcf_star(i,j,k)        = qcf_inc(i,j,nfor,k)    + qcf(i,j,k)

            ! At present qcf2_inc etc are not used as there is
            ! no forcing option.

            qcf2_star(i,j,k)   = qcf2_star(i,j,k)   + qcf2(i,j,k)
            qrain_star(i,j,k)  = qrain_star(i,j,k)  + qrain(i,j,k)
            qgraup_star(i,j,k) = qgraup_star(i,j,k) + qgraup(i,j,k)

            cf_star (i,j,k) = cf_star (i,j,k)                           &
     &                      + bulk_cloud_fraction(i,j,k)
            cfl_star(i,j,k) = cfl_star(i,j,k)                           &
     &                      + cloud_fraction_liquid(i,j,k)
            cff_star(i,j,k) = cff_star(i,j,k)                           &
     &                      + cloud_fraction_frozen(i,j,k)
          End Do
        End Do
      End Do

!---------------------------------------------------------------------
!       Convert t_inc_scm back to theta_star for call to Atmos_Physics2
!       theta_star, q_star_scm, qcl_star and qcf_star now all store
!       Increment + Whole Values
!---------------------------------------------------------------------

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                theta_star(i,j,k) = t_inc_scm(i,j,nfor,k)               &
     &                            / exner_theta_levels(i,j,k)
              End Do
            End Do
          End Do

            ! Convert X_star to mixing ratios
! DEPENDS ON: q_to_mix
            call q_to_mix (row_length, rows, wet_model_levels,          &
     &           halo_i, halo_j,                                        &
     &           q_star, qcl_star, qcf_star,                            &
     &           qcf2_star, qrain_star, qgraup_star,                    &
     &           L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                 &
     &           mix_v_star, mix_cl_star, mix_cf_star,                  &
     &           mix_cf2_star, mix_rain_star, mix_graup_star            &
     &          )
!
!---------------------------------------------------------------------
!         VARIABLE VALUES:
!         qcl_star   = liq. water       (at timelevel n)
!                    + liq. water incs  (atmos_physics1 & forcing)
!         t_inc_scm  = temp.            (at timelevel n)
!                    + temp. incs       (atmos_physics1 & forcing)
!         q_star_scm = vapour           (at timelevel n)
!                    + vapour incs      (atmos_physics1 & forcing)
!         theta_star = pot. temp.       (at timelevel n)
!                    + pot. temp. incs  (atmos_physics1 & forcing)
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!       Call physics2
!---------------------------------------------------------------------

#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
            ntrad=a_sw_radstep_prog
#endif
! DEPENDS ON: atmos_physics2
          CALL Atmos_Physics2(                                          &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, row_length, rows                  &
     &, proc_row_group, proc_col_group, at_extremity                    &
     &, n_proc, n_procx, n_procy,  neighbour, rows                      &
     &, row_length, g_datastart, me, NumCycles, CycleNo                 &

! model dimensions.
     &, row_length, rows, rows, land_points, model_levels               &
     &, nice, wet_model_levels, boundary_layer_levels                   &
     &, st_levels, sm_levels, cloud_levels, land_ice_points             &
     &, soil_points, n_cca_levels, ntiles, tr_levels                    &
     &, first_constant_r_rho_level, dim_cs1, dim_cs2                    &

!IN Substepping Information
     &, Num_Substeps, L_phys2_substep                                   &

! Model switches
     &, l_regular, l_mr_physics2, model_domain, L_dry,                  &
     &  FORMDRAG, OROG_DRAG_PARAM, LCAL360, L_emcorr, Ltimer            &
     &, L_us_blsol, BL_OPTIONS, L_cld_area, L_ACF_Cusack, L_ACF_Brooks  &
     &, L_RHCPT, L_hydrology, L_bl, L_MURK, L_MURK_ADVECT, L_MURK_SOURCE&
     &, L_BL_TRACER_MIX, L_DUST, L_CAM_DUST, L_sulpc_so2 , L_sulpc_nh3  &
     &, L_SBLeq, L_SBLco, L_sulpc_dms, L_soot, L_biomass, L_ocff        &
     &, L_co2_interactive, l_ctile, l_co2_emits, L_use_bl_diag_term     &
     &, L_pc2, L_pc2_reset, L_eacf, L_LAMBDAM2, L_FULL_LAMBDAS, L_ukca  &
     &, L_sice_heatflux, L_OHADGEM1, L_USE_CARIOLLE                     &
     &, l_anthrop_heat_src                                              &

! model Parameters
     &, alpha_cd, Puns, Pstb, rhcrit, CO2_MMR                           &
     &, tr_vars, tr_ukca, Muw_SBL, Mwt_SBL, cloud_fraction_method       &
     &, overlap_ice_liquid, ice_fraction_method, ctt_weight             &
     &, t_weight, qsat_fixed, sub_cld, x1i, x1ic, x1r, x2r, x4r, l_psd  &
     &, ai, bi, aic, bic, lsp_ei, lsp_fi, lsp_eic, lsp_fic              &
     &, dbsdtbs_turb_0, Charnock, SeaSalinityFactor                     &

! Physical constants
     &, lc,lf,cp,two_Omega,p_zero,kappa,R,g,lapse,earth_radius,Pi       &

! in coordinate information
     &, r_rho_levels, r_theta_levels, r_rho_levels, r_rho_levels        &
     &, unscl_dry_rho                                                   &
     &, eta_theta(1:model_levels+1), eta_rho(1:model_levels)            &
     &, delta_lambda, delta_phi, dlambda_p, dphi_p                      &
     &, wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v                    &
     &, lat_rot_NP, long_rot_NP, f3_at_u                                &

! Variables required by STPH_SCV
     &, mix_v_star, mix_cl_star, mix_cf_star                            &

! in time stepping information.
     &, timestep, yearno, day, ihour, imin, idum0, stepcount            &


! trig arrays
     &, sin_longitude, cos_longitude, sin_latitude, FV_cos_latitude, sec_latitude     &

! River routing
     &, row_length, row_length, l_rivers, xpa,xua,xva,ypa,yua,yva       &
     &, g_p_field, g_r_field, a_steps_since_riv, river_row_length       &
     &, river_rows, global_river_row_length, global_river_rows          &
     &, river_step, river_vel, river_mcoef                              &
     &, trivdir, trivseq, twatstor, inlandout_atm, l_inland,            &

! Grid-to-grid river routing
     &  r_area, slope, flowobs1, r_inext, r_jnext, r_land,              &
     &  substore, surfstore, flowin, bflowin,                           &

! diagnostics info
#include "argsts.h"
     &  STASHwork3,STASHwork5,STASHwork8,STASHwork9,STASHwork19         &
     &, STASHwork26                                                     &
!
! SCM diagnostics
     &, nSCMDpkgs, L_SCMDiags                                           &

! in data fields.
     &, theta_scm, q, qcl, qcf, rho, u, v, w, w_adv, p, p_star          &
     &, exner_rho_levels, exner_theta_levels                            &
     &, land_sea_mask, p_theta_levels                                   &

! Variables for subgrid turbulence scheme
     &, visc_BL_m, fm_3d, fh_3d, .false., .false., .false., max_diff    &
     &, turb_startlev_vert, turb_endlev_vert, bl_coef_km, bl_coef_kh    &

! ancillary fields and fields needed to be kept from timestep to
! timestep

     &, land_index, land_ice_index, soil_index, canopy_gb               &
     &, snodep, hcon, hcap, v_crit, v_wilt, v_sat , sthf                &
     &, sthu, sil_orog_land                                             &
     &, ho2r2_orog, di, ice_fract, u_0, v_0, u_0, v_0, cca, iccb, icct  &
     &, cclwp, ccw, lcbase, t_deep_soil, tsi, ti_n                      &
     &, tstar, z0msea, ice_fract_n, ice_thick_n, satcon, sathh, b_exp   &
     &, smcl, t1_sd, q1_sd, zh, area_cloud_fraction, bulk_cloud_fraction&
     &, cloud_fraction_liquid, cloud_fraction_frozen, ls_rain, ls_snow  &
     &, micro_tends, photosynth_act_rad, rad_hr, surf_radflux           &
     &, soil_clay,soil_silt,soil_sand,dust_mrel1,dust_mrel2,dust_mrel3  &
     &, dust_mrel4,dust_mrel5,dust_mrel6                                &
     &, so2_high_level, so2_em, nh3_em, dms_em, soot_hilem              &
     &, soot_em, ocff_hilem, ocff_em, co2_emits, co2flux                &

! tracer fluxes - kdcorbin, 05/10
     &, tracer_flux1, tracer_flux2, tracer_flux3, tracer_flux4          &
     &, tracer_flux5, tracer_flux6, tracer_flux7, tracer_flux8          &
     &, tracer_flux9, tracer_flux10,tracer_flux11,tracer_flux12         &
     &, tracer_flux13,tracer_flux14,tracer_flux15,tracer_flux16         &
     &, tracer_flux17,tracer_flux18,tracer_flux19,tracer_flux20         &
     &, CO2EMITMASS                                                     &
! rml 1/7/13 flag for co2 flux into passive tracer
     &, L_CO2_TRACER                                                    &

! in/out
     &, theta_star                                                      &
     &, q_star_scm(1:row_length,1:rows,nfor:nfor,1:wet_model_levels)    &
     &, qcl_star, qcf_star, cf_star, cfl_star, cff_star                 &
     &, u_inc_scm(1:row_length,1:rows,nfor:nfor,1:model_levels)         &
     &, v_inc_scm(1:row_length,1:rows,nfor:nfor,1:model_levels)         &
     &, rdum4, sum_eng_fluxes, sum_moist_flux                           &

! In/Out tracer fields
     &, aerosol(1:row_length,1:rows,1:wet_model_levels), FREE_TRACERS   &
     &, ukca_tracers                                                    &
     &, dust_div1, dust_div2, dust_div3, dust_div4, dust_div5,dust_div6 &
     &, SO2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new         &
     &, soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld           &
     &, ocff_new, ocff_aged, ocff_cld, co2                              &

! IN/OUT STPH_RP AND RIVERS
     &, g0_rp, par_mezcla                                               &

! IN/OUT RIVERS
     &, tot_surf_runoff, tot_sub_runoff                                 &

! out fields
     &, rhokm, cH_term, ntml, cumulus, nbdsc, ntdsc                     &
     &, rhcpt,row_length, rows                                          &

! Additional variables for MOSES II
     &, frac_typ, frac_disturb, canht, lai                              &
     &, canopy(1:land_points,1:ntiles)                                  &
     &, catch(1:land_points,1:ntiles), catch_snow, snow_grnd, snow_tile &
     &, z0_tile, tstar_tile, infil_tile(1:land_points,1:ntiles)         &
     &, rgrain(1:land_points,1:ntiles)                                  &
     &, cs, gs, co2_dim_row, co2_dim_len, l_neg_tstar, l_snow_albedo    &
     &, .false., .false., l_trif_eq, L_Q10, asteps_since_triffid        &
     &, stepcount, phenol_period, triffid_period, can_model, g_leaf_acc &
     &, g_leaf_phen_acc, npp_ft_acc, resp_w_ft, resp_s_acc              &
     &, land_pts_trif, npft_trif, dolr, lw_down, sw_tile, fland_ctile   &
     &, tstar_land, tstar_sea, tstar_sice                               &
     &, albsoil, cos_zenith_angle, can_rad_mod, ilayers                 &

     &,radnet_tile                                                      & 

!    EAK
!    IN
     &, surf_down_sw,alb_tile,l_tile_pts                                &
!     &, surf_down_sw,alb_tile,cos_zenith_angle                          &
     &, lat,long,day,time_sec,SW_DOWN                                   &
     &, SNOW_DEPTH3L,SNOW_MASS3L,SNOW_COND,SNOW_TMP3L                   &
     &, SNOW_RHO3L,SNOW_RHO1L,SMCL_TILE,STHU_TILE,STHF_TILE             &
     &, TSOIL_TILE,T_SURF_TILE,HCONS                                    &
     &, SOIL_TYPE,VEG_TYPE                                              &
     &, ISNOW_FLG3L,total_nsteps                                        &
     &           ,FTL_TILE_CAB,FTL_CAB,LE_TILE_CAB,LE_CAB               &
     &           ,TSTAR_TILE_CAB,TSTAR_CAB,SMCL_CAB,TSOIL_CAB           &
     &           ,USTAR_CAB,SURF_HTF_CAB                                &
     &, l_cable                                                         &
!sxy
     &           ,TOT_ALB                                               &
     &           ,U_S_CAB,CH_CAB,CD_CAB                                 &
     &           ,CD,CH                                                 &
     &           ,TILE_PTS,TILE_INDEX                                   &
     &           ,SNAGE_TILE,RTSOIL_TILE                                &
     &           ,GFLUX_TILE,SGFLUX_TILE                                &

! Additional variables required for large-scale hydrology:
     &, L_TOP, L_PDM, FEXP, GAMTOT, TI_MEAN, TI_SIG, FSAT, FWETL, ZW    &
     &, STHZW,A_FSAT,C_FSAT,A_FWET,C_FWET,L_SOIL_SAT_DOWN               &
! Cariolle ozone
     &, OZONE_TRACER                                                    &

! Prescribed surface forcing and roughness lengths
     &, L_flux_bc                                                       &
     &, flux_e_scm(1,1:rows,nfor:nfor),flux_h_scm(1,1:rows,nfor:nfor)   &
     &, L_spec_z0,z0m_scm,z0h_scm                                       &

! Additional variables for SCM diagnostics
     &, nfor, ntrop, time_string, daycount, layer_cloud                 &
     &, tstar, nconvars, ntrad, ntrad1, conv_mode                       &

! error information
     &, Error_code  )
!
!---------------------------------------------------------------------
!         VARIABLE VALUES
!         qcl_star   = liq. water          (at timelevel n)
!                    + liq. water incs     (atmos_physics1, forcing &
!                                           atmos_physics2)
!         q_star_scm = vapour              (at timelevel n)
!                    + vapour incs         (atmos_physics1, forcing &
!                                           atmos_physics2)
!         theta_star = pot. temp.          (at timelevel n)
!                    + pot. temp. incs     (atmos_physics1, forcing &
!                                           atmos_physics2)
!---------------------------------------------------------------------
!

          If (stats .OR. obs .OR. noforce) Then

! rho and pressure out of sync with T, but can't update with IdealGL
! settings i.e. P=rho.R.T as this will cause a crash in the physics
! routines. I assume that Pressure and density need to be updated
! using exner_prime


          ! Update theta_scm and q and winds
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length

                theta_scm(i,j,k) = theta_star(i,j,k)

                !  Convert theta_scm back to t
                t(i,j,k) = theta_scm(i,j,k) * exner_theta_levels(i,j,k)

                ! add wind increments to winds
                u(i,j,k) = u(i,j,k) + u_inc_scm(i,j,nfor,k)
                v(i,j,k) = v(i,j,k) + v_inc_scm(i,j,nfor,k)

                ! Vertical wind either prescribed for vertical advection
                ! forcing and updated after s_forcing or left at initial
                ! value.

              End Do
            End Do
          End Do

          Do k = 1, wet_model_levels
            Do j = 1, rows
              Do i = 1, row_length
                q      (i,j,k) = q_star_scm  (i,j,nfor,k)
                qcl    (i,j,k) = qcl_star    (i,j,k)
                qcf    (i,j,k) = qcf_star    (i,j,k)
                qcf2   (i,j,k) = qcf2_star   (i,j,k)
                qrain  (i,j,k) = qrain_star  (i,j,k)
                qgraup (i,j,k) = qgraup_star (i,j,k)

                If (L_pc2) Then
                  bulk_cloud_fraction   (i,j,k) = cf_star (i,j,k)
                  cloud_fraction_liquid (i,j,k) = cfl_star(i,j,k)
                  cloud_fraction_frozen (i,j,k) = cff_star(i,j,k)
                End If  ! L_pc2

              End Do
            End Do
          End Do

          End If   ! stats.or.obs.or.noforce
!
!---------------------------------------------------------------------
!         VARIABLE VALUES FOR THE NON-PC2 SITUATION.
!         (FOR PC2 SEE NEXT SECTION OF CODE FOR ADDITIONAL TERMS)
!
!         qcl       = liq. water  (at timelevel n+1)
!         q         = vapour      (at timelevel n+1)
!         theta_scm = pot. temp.  (at timelevel n+1)
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
!         PC2 section
!---------------------------------------------------------------------
!
          If (L_pc2) Then
!
! The PC2 section needs to
! 1) increment condensation and cloud fractions due to the forcing
! 2) initiate cloud.
! 3) set area_cloud_fraction to bulk_cloud_fraction and correct
!    theta_scm.
!
! 1. Calculate condensation increments which result from the forcing
! by using the homogeneous forcing approach
!
! DEPENDS ON: pc2_homog_plus_turb
            CALL PC2_HOMOG_PLUS_TURB(                                   &
     &        p_theta_levels(1,1,1:wet_model_levels),                   &
     &        wet_model_levels,                                         &
     &        row_length, rows, timestep, t, bulk_cloud_fraction,       &
     &        cloud_fraction_liquid, cloud_fraction_frozen,             &
     &        q, qcl, t_forcing,                                        &
     &        q_forcing, qcl_forcing, p_forcing,                        &
     &        0.0, 0.0, l_mr_pc2)
!
!
! 1a. We have already applied the forcing to q, qcl and t in the
!     forcing section. PC2_homog_plus_turb has now done this again
!     so we need to subtract off these increments. Just the
!     condensation will therefore remain.
!
            Do k = 1, wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  q(i,j,k)   = q(i,j,k)   - q_forcing(i,j,k)
                  qcl(i,j,k) = qcl(i,j,k) - qcl_forcing(i,j,k)
                  t(i,j,k)   = t(i,j,k)   - t_forcing(i,j,k)
                End Do ! i
              End Do ! j
            End Do ! k
!
!
! 2. Initiate cloud if necessary and check for sensible values.
!    Use a dummy argument to receive the increment information.
!
! DEPENDS ON: pc2_initiation_ctl
            Call pc2_initiation_ctl(                                    &
     &        0, 0, 0, 0                                                &
     &,       row_length, rows, row_length, rows                        &
     &,       model_levels, wet_model_levels, .false.                   &
     &,       l_mr_pc2, L_ACF_CUSACK, timestep                          &
     &,       nSCMDpkgs,L_SCMDiags                                      &
     &,       T, q, qcl, qcf                                            &
     &,       bulk_cloud_fraction,cloud_fraction_liquid                 &
     &,       cloud_fraction_frozen, rhts, tlts,qtts,ptts               &
     &,       area_cloud_fraction,p,p_star,p_theta_levels(1,1,1)        &
     &,       iccb, cumulus                                             &
     &,       rhcpt,rdum4,rdum4,rdum4,rdum4,rdum4,rdum4,rdum4           &
     &        )
!
!
!
           If (.not. L_CLD_AREA) Then

! 3. Set area_cloud_fraction to bulk_cloud_fraction and update
!    theta_scm from new temperature
!
            Do k = 1, wet_model_levels
              Do j = 1, rows
                Do i = 1, row_length
                  area_cloud_fraction(i,j,k) = bulk_cloud_fraction(i,j,k)
                  theta_scm(i,j,k) = t(i,j,k)/exner_theta_levels(i,j,k)
                End Do ! i
              End Do ! j
            End Do ! k

            Else If (L_CLD_AREA) Then

              If (L_ACF_Brooks) Then

! DEPENDS ON: ls_acf_brooks
                 Call LS_ACF_Brooks (                                   &
     &             halo_i, halo_j, off_x, off_y                         &
     &,            row_length, rows, model_levels, wet_model_levels     &
     &,            r_theta_levels, delta_lambda, delta_phi              &
     &,            FV_cos_latitude                                      &
     &,            bulk_cloud_fraction, cloud_fraction_liquid           &
     &,            cloud_fraction_frozen                                &
     &,            cumulus                                              &
     &,            area_cloud_fraction )

           End If ! L_ACF_Brooks

         End If ! L_CLD_AREA

!
          End If  ! L_pc2
!
!---------------------------------------------------------------------
!         End of PC2 section.
!
!         Variables q, qcl, t, theta_scm and cloud fractions are now
!         all at timelevel n+1
!---------------------------------------------------------------------
!
          If (main_diag_switch /= 0) Then
!
!-----------------------------------------------------------------------
!           SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------
            If (L_SCMDiags(SCMDiag_PC2)) Then

! DEPENDS ON: scmoutput
              Call SCMoutput(qcl,                                       &
                   'qcl_n1_notpc2','qcl at end of timestep before pc2', &
                   'kg/kg',t_inst,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(q,                                         &
                   'q_n1_notpc2','q at end of timestep before pc2',     &
                    'kg/kg',t_inst,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(theta_scm,                                 &
                  'th_n1_notpc2','theta at end of timestep before pc2', &
                  'kg/kg',t_inst,d_all,default_streams,'',RoutineName)

            End If ! L_SCMDiags(SCMDiag_PC2)

          End If ! main_diag_switch /= 0
!---------------------------------------------------------------------
!         Calculate some final diagnostics, and write the lot
!         out if needs be.
!---------------------------------------------------------------------

          If (main_diag_switch /= 0) Then
!
!-----------------------------------------------------------------------
!           SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------
            If (L_SCMDiags(SCMDiag_PC2)) Then

! DEPENDS ON: scmoutput
              Call SCMoutput(qcl,                                       &
                   'qcl_n1_afterpc2',                                   &
                   'qcl at end of timestep after pc2','kg/kg',          &
                   t_inst,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(q,                                         &
                   'q_n1_afterpc2',                                     &
                   'q at end of timestep after pc2','kg/kg',            &
                   t_inst,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(theta_scm,                                 &
                   'th_n1_afterpc2',                                    &
                   'theta at end of timestep after pc2','K',            &
                   t_inst,d_all,default_streams,'',RoutineName)

            End If ! L_SCMDiags(SCMDiag_PC2)
!
!-----------------------------------------------------------------------
!           SCM Increments Diagnostics Package
!-----------------------------------------------------------------------
            If (L_SCMDiags(SCMDiag_incs)) Then

              Do k=1,model_levels
                Do j=1,rows
                  Do i=1,row_length
                    t_totalinc(i,j,k)=(t(i,j,k) - t_start(i,j,k))
                    u_totalinc(i,j,k)=(u(i,j,k) - u_start(i,j,k))
                    v_totalinc(i,j,k)=(v(i,j,k) - v_start(i,j,k))
                    w_totalinc(i,j,k)=(v(i,j,k) - w_start(i,j,k))
                  End Do
                End Do
              End Do
              Do k=1,wet_model_levels
                Do j=1,rows
                  Do i=1,row_length
                    q_totalinc(i,j,k)=(q(i,j,k)     - q_start(i,j,k))
                    qcl_totalinc(i,j,k)=(qcl(i,j,k) - qcl_start(i,j,k))
                    qcf_totalinc(i,j,k)=(qcf(i,j,k) - qcf_start(i,j,k))
                  End Do
                End Do
              End Do

! DEPENDS ON: scmoutput
              Call SCMoutput(t_totalinc,                                &
                   'dt_total',                                          &
                   'Total increment to T','K',                          &
                   t_avg,d_all,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(u_totalinc,                                &
                   'du_total',                                          &
                   'Total increment to u','m/s',                        &
                   t_avg,d_all,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(v_totalinc,                                &
                   'dv_total',                                          &
                   'Total increment to v','m/s',                        &
                   t_avg,d_all,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(w_totalinc,                                &
                   'dw_total',                                          &
                   'Total increment to w','m/s',                        &
                   t_avg,d_all,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(q_totalinc,                                &
                   'dq_total',                                          &
                   'Total increment to q','kg/kg',                      &
                   t_avg,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(qcl_totalinc,                              &
                   'dqcl_total',                                        &
                   'Total increment to qcl','kg/kg',                    &
                   t_avg,d_wet,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
              Call SCMoutput(qcf_totalinc,                              &
                   'dqcf_total',                                        &
                   'Total increment to qcf','kg/kg',                    &
                   t_avg,d_wet,default_streams,'',RoutineName)

              If (geoforce) Then
! DEPENDS ON: scmoutput
                Call SCMoutput(uinc_geo,                                &
                   'du_geo',                                            &
                   'Geostrophic forcing increment to u','m/s',          &
                   t_avg,d_all,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
                Call SCMoutput(vinc_geo,                                &
                   'dv_geo',                                            &
                   'Geostrophic forcing increment to v','m/s',          &
                   t_avg,d_all,default_streams,'',RoutineName)
              End If

            End If ! L_SCMDiags(SCMDiag_incs)


             ! Store some diagnostic
! DEPENDS ON: dgnstcs_scm_main
             Call dgnstcs_scm_main(                                           &
               row_length, rows, land_points, model_levels,                   &
               wet_model_levels, sm_levels, st_levels, ntype,                 &
               r_theta_levels, r_rho_levels, rho, timestep, u, v, T,          &
               theta_scm, q, qcl, qcf, layer_cloud, cca, ccw,                 &
               t_deep_soil, p_star, tstar, smc, canopy_gb, snodep, zh,        &
               z0msea, smcl, sthu, sthf, gs, lw_incs, surf_radflux,           &
               photosynth_act_rad, tstar_tile, aerosol, p_theta_levels,       &
               p, iccb, icct, w, w_adv, area_cloud_fraction,                  &
               bulk_cloud_fraction, cloud_fraction_liquid,                    &
               cloud_fraction_frozen, cclwp,                                  &
               ftl_tile_cab, ftl_cab, le_tile_cab, le_cab,                    &
               tstar_tile_cab, tstar_cab, smcl_cab, tsoil_cab,                &
               ustar_cab, surf_htf_cab, snow_rho1l, tot_alb,                  &
               u_s_cab, ch_cab, cd_cab,                                       &
               cd,ch, sw_down, radnet_tile, snow_tile, snow_tmp3l,            &
               snow_rho3l, snow_depth3l, lw_down,                             &
               nSCMDpkgs, L_SCMDiags)

             ! Initialise the diagnostic output files if this is the
             ! end of the first timestep for which the system was on
             If (SCMop%first_pass) Then
                SCMop%first_pass=.false.

                ! The list of diagnostics should be finalised now, so
                ! dump_streams_init can be called.
! DEPENDS ON: dump_streams_init
                Call DUMP_STREAMS_INIT(SCMop,                     & ! INOUT
     &             row_length,rows,model_levels,wet_model_levels, & ! IN
     &             boundary_layer_levels,cloud_levels,            & ! IN
     &             ozone_levels,st_levels,sm_levels,ntiles,       & ! IN
     &             year_init,month_init,day_init,                 & ! IN
     &             hour_init,min_init,sec_init,                   & ! IN
     &             timestep,ndayin,nminin,nsecin,sec_day,         & ! IN
     &             ndayin*full_daysteps+nstepsin,                 & ! IN
     &             a_sw_radstep_prog,a_sw_radstep_diag,           & ! IN
     &             ntrad,                                         & ! IN
     &             z_top_of_model,first_constant_r_rho_level,     & ! IN
     &             eta_theta,eta_rho,orog,r_theta_levels,         & ! IN
     &             r_rho_levels,netcdf_chunksize)                   ! IN

                !=======================================================
                ! ADD REMINDER TO CHECK NFOR IS CORRECT
                ! Added here it doesn't get push off screen
                ! and missed by user
                !=======================================================

                Write(6,*) "========================================"
                Write(6,*) " Namelist NFOR = " // TRIM(ADJUSTL(sdum0))
                Write(6,*) " NOTE: INCORRECT SPECIFICATION OF NFOR  "
                Write(6,*) "       WILL PRODUCE UNINTENDED RESULTS. "
                Write(6,*) "========================================"
                Write(6,*) " "


             End If

             ! Write output diagnostics to file(s)
! DEPENDS ON: dump_streams
             Call DUMP_STREAMS(SCMop,day,time_sec,                      &
     &            row_length,rows,model_levels,dayno_init,              &
     &            int(stepcount*timestep),site,lat,long,                &
     &            time_initi,year,lcal360)

          End If ! main_diag_switch /= 0

          If (l_ts_log) Then
            Call scm_message('Complete')
          End If

        End Do   ! stepcount

      End Do   ! daycount

      ! Close the output files
! DEPENDS ON: dump_streams_end
      Call DUMP_STREAMS_END(SCMop)


      ! Deallocate Arrays
      Deallocate(arcl)    ! Aerosol climatology array for NWP


      Write(6,*) '------------------------'
      Write(6,*) 'Run completed'
      Write(6,*) '------------------------'


999   Return
      END SUBROUTINE SCM_main

#endif
