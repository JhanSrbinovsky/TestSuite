
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE CONV_DIAG-----------------------------------------------
!!!
!!!  PURPOSE : TO DIAGNOSE CONVECTION OCCURRENCE AND TYPE
!!!
!!!  SUITABLE FOR SINGLE COLUMN MODEL USE
!!!
!!!  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
!!!
!!!  DOCUMENTATION :
!!!
!!!
!!!--------------------------------------------------------------------

!    ARGUMENTS---------------------------------------------------------
!
      SUBROUTINE CONV_DIAG(                                             &

! IN values defining field dimensions and subset to be processed :
     & row_length, rows, n_rows, halo_i, halo_j, off_x, off_y           &
     &, global_row_length, proc_row_group, at_extremity                 &

! IN values defining vertical grid of model atmosphere :
     &, model_domain, bl_levels, model_levels, wet_model_levels         &
     &, p, P_theta_lev, exner_rho, r_rho_levels,r_theta_levels          &
     &, rho_only, rho_theta, z_full, z_half                             &
     &, sin_theta_longitude,cos_theta_longitude                         &

! IN Model switches
     &, l_mixing_ratio                                                  &

! IN Cloud data :
     &, qcf, qcl, cloud_fraction                                        &

! IN everything not covered so far :
     &, pstar, q, theta, exner_theta_levels, u, v, u_0_p, v_0_p         &
     &, L_flux_bc, flux_e, flux_h, L_spec_z0, z0h_scm                   &
     &, tstar, land_mask, timestep                                      &
     &, w_copy                                                          &

! SCM Diagnostics (dummy values in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! OUT data required elsewhere in UM system :
     &, zh,zhpar,z_lcl,z_lcl_uv,delthvu,ql_ad,ntml,ntpar,NLCL,cumulus   &
     &, L_shallow,l_congestus, l_congestus2, cin, CAPE                  &

     &, ERROR                                                           &
     & )

      Use cv_run_mod, Only:                                             &
          icvdiag

      IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
! (a) Defining horizontal grid and subset thereof to be processed.

      integer, intent(in) ::                                            &
     &  row_length                                                      &
                           ! Local number of points on a row
     &, rows                                                            &
                           ! Local number of rows in a theta field
     &, n_rows                                                          &
                           ! Local number of rows in a v field
     &, halo_i                                                          &
                           ! Size of halo in i direction.
     &, halo_j                                                          &
                           ! Size of halo in j direction.
     &, off_x                                                           &
                           ! Size of small halo in i
     &, off_y                                                           &
                           ! Size of small halo in j.
     &, global_row_length                                               &
                           ! number of points on a row
     &, proc_row_group     ! Group id for processors on the same row

      Logical,  intent(in) ::                                           &
     &  at_extremity(4)    ! Indicates If this processor is at north,
                           !  south,east or west of the processor grid

! (b) Defining vertical grid of model atmosphere.


      integer, intent(in) ::                                            &
     &  model_domain                                                    &
     &, bl_levels                                                       &
                           ! Max. no. of "boundary" levels allowed.
     &, model_levels                                                    &
                           ! number of model levels
     &, wet_model_levels
                           ! number of wet model levels

      real, intent(in) ::                                               &
     &  P(1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels)    &
                   ! pressure on rho levels  (Pa)
     &, P_theta_lev(row_length, rows, model_levels)                     &
                   ! pressure on theta levels  (Pa)
     &, exner_rho(1-off_x:row_length+off_x,1-off_y:rows+off_y, model_levels) &
                   ! Exner pressure on rho levels 
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j,model_levels)                 &
                   ! radius of rho levels (m)
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
                   ! radius of theta levels (m)
     &, rho_only(row_length,rows,model_levels)                          &
                                               ! density (kg/m3)
     &, rho_theta(row_length,rows,model_levels-1)                       &
                                               ! rho th lev (kg/m3)
     &, z_full(row_length,rows,model_levels)                            &
                                               ! height th lev (m)
     &, z_half(row_length,rows,model_levels)                            &
                                               ! height rho lev (m)
     &, sin_theta_longitude (row_length,rows)                           &
                                                 ! sin(long)
     &, cos_theta_longitude (row_length,rows)    ! cos(long)

      Logical,  intent(in) ::                                           &
     &  l_mixing_ratio                                                  &
                           ! true moisture input as mixing ratios
                           ! false moisture input as specific humidity
     &, L_flux_bc                                                       &
                           ! true if SCM using specified surface fluxes
     &, L_spec_z0          ! true if roughness length has been specified

! (c) Cloud data.

      real, intent(in) ::                                               &
     & QCF(row_length,rows,wet_model_levels)                            &
!                                 ! Cloud ice (kg per kg air)
     &,QCL(row_length,rows,wet_model_levels)                            &
!                                 ! Cloud liquid water (kg/kg air).
     &,cloud_fraction(row_length, rows, wet_model_levels)
!                                 ! Cloud fraction

! (d) Atmospheric + any other data not covered so far, incl control.

      real, intent(in) ::                                               &
     &  pstar(row_length, rows)                                         &
                                       ! Surface pressure (Pascals).
     &, Q(row_length,rows,wet_model_levels)                             &
                                             ! water vapour (kg/kg)
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
                                       ! Theta (Kelvin)
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
                                       ! exner pressure theta lev (Pa)
     &, U(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               & 
                                       !  W'ly wind component (m/s).
     &, V(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               & 
                                       !  S'ly wind component (m/s).
     &, U_0_P(row_length,rows)                                          & 
                                  ! W'ly component of surface current
                                  !    (metres per second) on P-grid.
     &, V_0_P(row_length,rows)                                          & 
                                  ! S'ly component of surface current
                                  !    (metres per second) on P-grid.
     &,flux_e(row_length,rows)                                          &
                                 ! Specified surface
                                 !    latent heat flux (W/m^2)
     &,flux_h(row_length,rows)                                          &
                                 ! Specified surface
                                 !    sensible heat fluxes (in W/m2)
     &,z0h_scm(row_length,rows)  ! Namelist input z0h (if >0)

      real, intent(inout) ::                                            &
     &  tstar(row_length,rows)    ! Surface temperature (= top soil
                                  !  layer temperature) (K).
                                  ! NOTE only inout for SCM
      logical, intent(in) ::                                            &
     & land_mask(row_length, rows)  ! T If land, F Elsewhere.

      real, intent(in) ::                                               &
     &  timestep                                                        &
                                     ! timestep (seconds).
     &, w_copy(row_length,rows,                                         &
     &          0:model_levels)      ! vertical velocity (m/s)

! Additional variables for SCM diagnostics which are dummy in full UM
      integer                                                           &
     &  nSCMDpkgs          ! No of diagnostics packages
!
      logical                                                           &
     &  L_SCMDiags(nSCMDpkgs) ! Logicals for diagnostics packages
!
      real, intent(inout) ::                                            &
     & zh(row_length,rows)    ! Height above surface of top
                              !  of boundary layer (metres).
      real, intent(out) ::                                              &
     & zhpar(row_length,rows)                                           &
                              ! Height of max parcel ascent (m)
     &,z_lcl(row_length,rows)                                           &
                              ! Height of lIfting condensation level (m)
     &,z_lcl_uv(row_length,rows)                                        &
                                 ! Height of lIfting condensation
                                 !     level on uv grid (m)
     &,delthvu(row_length,rows)                                         &
                                 ! Integral of undilute parcel buoyancy
                                 ! over convective cloud layer
                                 ! (for convection scheme)
     &,ql_ad(row_length,rows)                                           & 
                                 ! adiabatic liquid water content at
                                 ! inversion or cloud top (kg/kg)
                                 ! NOT USED for 4A code
     &, cape(row_length, rows)                                          &
                                 ! CAPE from parcel ascent (m2/s2)
     &, cin(row_length, rows)    ! CIN from parcel ascent (m2/s2)

      integer, intent(out)  ::                                          &
     & ntml(row_length,rows)                                            &
                               ! Number of model levels in the
                               ! turbulently mixed layer.
     &,ntpar(row_length,rows)                                           &
                               ! Max levels for parcel ascent
     &,nlcl(row_length,rows)   ! No. of model layers below the
                               ! lifting condensation level.

      logical, intent(out) ::                                           &
     & cumulus(row_length,rows)                                         &
                                  ! Logical indicator for convection
     &,L_shallow(row_length,rows)                                       &
                                  ! Logical indicator for shallow Cu
     &,L_congestus(row_length,rows)                                     &
                                     ! Logical indicator for congestus Cu
     &,L_congestus2(row_length,rows) ! Logical ind 2 for congestus Cu

      integer, intent(out)  ::                                          &
     & ERROR          ! 0 - no error in this routine

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

      Character*(*), Parameter ::  RoutineName = 'conv_diag'

      integer ::  &
     & I,J        & ! LOCAL Loop counter (horizontal field index).
     &, ii        & ! Local compressed array counter.
     &, K         & ! LOCAL Loop counter (vertical level index).
     &,MBL        & ! Maximum number of model layers allowed in the
!                   ! mixing layer; set to bl_levels-1.
     &,nunstable    ! total number of unstable points

      real ::                                                           &
     &  mag_vector_np(1)                                                &
     &, dir_vector_np(1)                                                &
     &, mag_vector_sp(1)                                                &
     &, dir_vector_sp(1)

! uncompressed arrays - all points

      real ::                                                           &
     & qs_star(row_length, rows)                                        &
                                 ! Saturated sp humidity at sea surface
!                                ! 1 virtual temperature (K).
     &,fb_surf(row_length, rows)                                        &
                                 ! Change in theta_v from surface
!                                ! to layer 1 (note diff from BL)
     &,U_P(row_length, rows)                                            &
                                 ! U(1) on P-grid.
     &,V_P(row_length, rows)                                            &
                                 ! V(1) on P-grid.
     &,dqsdt(row_length,rows)                                           &
                                 ! d(QSAT)/dT 
     &,Z0H(row_length,rows)      ! roughness length (m) (SCM only)

! Used in calculation to decide on unstable points 

      real ::          &
     & theta1          &  ! Potential temperature in layer 1
     &,Ushear          &  ! U wind shear from level 1 to surface
     &,Vshear          &  ! V wind shear from level 1 to surface
     &,wshr1           &  ! magnitude of surface wind shear
     &,wshr2           &  ! (wshr1)**2
     &,rhostar         &  ! surface air density
     &,theta_star      &  ! theta at surface
     &,wthvbar         &  ! surface buoyancy flux
     &,cd              &  ! bulk transfer coefficient for momentum
!                         ! (Approximation. Also assume CH(heat)=cd .)
     &,ustar           &  ! surface friction velocity
     &,ustar_eff       &  ! effective surface friction velocity
     &,w_m             &  ! 
     &,w_s_cubed       &  !
     &,wstar              ! convective velocity scale

! compressed arrays store only values for unstable columns

      integer ::                   &
     &  index_i(row_length*rows)   & ! column number of unstable points
     &, index_j(row_length*rows)     ! row number of unstable points

! compressed copys for unstable points - names as original array plus _c

      integer ::                                                        &
     & ntml_c(row_length*rows)                                          &
     &,nlcl_c(row_length*rows) 

      real ::                                                           &
     &  q_c(row_length*rows, model_levels)                              &
     &, qcl_c(row_length*rows, model_levels)                            &
     &, qcf_c(row_length*rows, model_levels)                            &
     &, z_full_c(row_length*rows, model_levels)                         &
     &, z_half_c(row_length*rows, model_levels)                         &
     &, exner_theta_levels_c(row_length*rows, model_levels)             &
     &, P_theta_lev_c(row_length*rows, model_levels)                    &
     &, r_theta_levels_c( row_length *rows, model_levels)               &
     &, z_lcl_c(row_length*rows)                                        &
     &, zh_c(row_length*rows)                                           &
     &, delthvu_c(row_length*rows)                                      &
     &, cape_c(row_length*rows)                                         &
     &, cin_c(row_length*rows) 

! Arrays only used for unstable calculations

      logical ::                                                        &
     & topbl(row_length*rows)                                           &
                                 ! Flag set when top of boundary layer
!                                ! is reached.
     &,topprof(row_length*rows)                                         &
                                 ! Flag set when top of ascent
!                                ! is reached.
     &,above_lcl(row_length*rows)                                       &
                                 ! Flag set when parcel above LCL.
     &,shmin(row_length*rows)                                           
                                 ! Flag for finding min in parcel
!                                ! buoyancy below 3km (for shallow Cu)

      real ::                                                          &
     & T(row_length*rows, model_levels)                                &
                                    ! temperature (from theta)
     &,TL(row_length*rows, model_levels)                               &
                                    ! Ice/liquid water temperature,
                                    ! but replaced by T in LS_CLD.
     &,QW(row_length*rows, wet_model_levels)                           &
                                    ! Total water content
     &,SVL(row_length*rows, wet_model_levels)                          &
                                    ! Liquid/frozen water virtual
                                    ! static energy over CP.
     &,env_svl(row_length*rows,wet_model_levels)                       &
                                    ! Density (virtual) static energy
                                    ! over CP for layer.
     &,par_svl(row_length*rows,wet_model_levels)                       &
                                    ! Density (virtual) static energy
                                    ! over CP of parcel for level.
     &,QS(row_length*rows, wet_model_levels)                           
!                                   ! Saturated sp humidity at pressure and
!                                   ! temperature of sucessive levels.
      real ::                           &
     & T_LCL(row_length*rows)           & ! Temperature at lifting 
                                          ! condensation level.
     &,P_LCL(row_length*rows)           & ! Pressure at lifting condensation
                                          ! level.
     &,sl_plume(row_length*rows)        & ! Liquid/frozen water static energy
                                          ! over CP for a plume rising without
                                          ! dilution from level 1.
     &,QW_plume(row_length*rows)        & ! QW for a plume rising without
                                          ! dilution from level 1.
     &,th_ref(row_length*rows)          & ! theta - reference value
     &,t_ref(row_length*rows)           & ! t - reference value
     &,qsat_lev(row_length*rows)        & ! qsat for reference temperature
     &,Dt_dens_parc_T(row_length*rows)  & ! t_dens_parc-t_dens_env at ntpar
     &,Dt_dens_parc_Tmin(row_length*rows) & ! t_dens_parc-t_dens_env at kshmin
     &,thv_pert(row_length*rows)        & ! threshold thv of parcel
     &,dtv_min(row_length*rows)         & ! min TV of parcel in cld layer
     &,tv1_sd(row_length*rows)          & ! Approx to standard dev of level
                                          ! 1 virtual temperature (K).
     &,ENV_SVL_KM1(row_length*rows)     & ! Density (virtual) static energy
                                          ! over CP for last layer considered.
     &,PAR_SVL_KM1(row_length*rows)     & ! Density (virtual) static energy
                                          ! over CP of parcel at last level
                                          ! considered.
! Added for improved parcel top - mainly used for finding an ascent
! capped by an inversion.
     &,Dt_dens_parc_T2(row_length*rows) & ! 2nd copy of Dt_dens_parc_T
     &,dtv_min2(row_length*rows)        & ! 2nd copy min TV of parcel
     &,delthvu2(row_length*rows)        & !  2nd copy
     &,zh2(row_length*rows)             & !  2nd copy
     &,max_buoy(row_length*rows)          ! max parcel buoyancy 

! arrays holding various key model levels

      integer ::                                                       &
     & kshmin(row_length*rows)                                         &
                                   ! Position of buoyancy minimum above
!                                  ! topbl (used for shallow Cu diag)
     &,kcucheck(row_length*rows)                                       &
                                   ! Position of level just below 2.5km
!                                  ! (used for gradient check to
!                                  !  diagnose convection)
     &,k_plume(row_length*rows)                                        &
                                   ! start level for surface-driven plume
     &,k_max(row_length*rows)                                          &
                                   ! level of max parcel buoyancy
     &,k_neutral(row_length*rows)                                      &
                                   ! level of neutral parcel buoyancy
     &,k_inv(row_length*rows)                                          &
                                   ! level from inversion testing
     &,freeze_lev(row_length*rows) ! freezing level

!
! parcel calculation
!
      real ::                                                          &
     & T_parc(row_length*rows, wet_model_levels)                       &
                         ! Temperature of parcel.
     &,t_dens_parc(row_length*rows, wet_model_levels)                  &
                         ! Density potential temperature of parcel.
     &,t_dens_env(row_length*rows, wet_model_levels)                   &
                         ! Density potential temperature of environment.
     &,denv_bydz(row_length*rows, wet_model_levels)                    &
                         ! Gradient of density potential
                         ! temperature in the environment.
     &,dpar_bydz(row_length*rows, wet_model_levels)                    &
                         ! Gradient of density potential
                         ! temperature of the parcel.
     &,buoyancy(row_length*rows, wet_model_levels)                     &
                         ! undilute parcel buoyancy (K)
     &,buoyancy_dil(row_length*rows, wet_model_levels)                 &
                         ! dilute parcel buoyancy (K)
     &,th_par_km1(row_length*rows)
                         ! parcel theta at level below

      real ::           &
     &  z_surf          &  ! approx height of top of surface layer
     &, vap_press       &  ! Vapour pressure.
     &, grad_cld        &  ! SVL gradient in layer above LCL.
     &, grad_sub        &  ! SVL gradient in layer below LCL.
     &, q_liq_env       &  ! Condensed water content of environment.
     &, dq_sat_env      &  ! DQSAT/DT for environment
     &, Lrcp_const      &  ! lc or lc+lf over cp
     &, lrcp_const_parc &  ! lc or lc+lf over cp
     &, L_const         &  ! lc or lc+lf
     &, dz              &  ! layer depth
     &, z_pr            &  ! used in estimating th_ref at next level
     &, th_par          &  ! theta value for parcel
     &, inc             &  ! CIN/CAPE increment for layer
     &, dq_sat_par      &  ! dqsat/dT for parcel
     &, temp_parc       &  ! average temperature of parcel after entrainment
     &, q_vap_parc      &  ! Vapour content of undilute parcel
     &, q_liq_parc      &  ! Liquid water content of undilute parcel
     &, qcl_parc        &  ! parcel qcl - dilute parcel cal
     &, qcf_parc        &  ! parcel qcf - dilute parcel cal 
     &, factor          &  ! multiplying factor 
     &, denv_bydz2      &  ! Gradient of density potential
!                          ! temperature in the environment.
     &, dpar_bydz2         ! Gradient of density potential
!                          ! temperature of the parcel.
     
! required for average w calculation

      real ::                                                          &
     &  dmass_theta(row_length*rows,model_levels)                      &
                     ! r**2rho*dr on theta levels
     &, w_avg(row_length*rows)                                         &
                                      ! mean w over layer (m/s)
     &, w_avg2(row_length*rows)                                        &
                                      ! mean w over layer (m/s)
     &, mass(row_length*rows)         ! mass for column

! Arrays added for dilute parcel calculation
      
      real ::                                                          &
     &  entrain_fraction(row_length*rows,model_levels)                 &
                     ! fraction of environmental air to mix with parcel
     &, qw_parc(row_length*rows,model_levels)                          &
                     ! parcel total water
     &, sl_parc(row_length*rows,model_levels)                          &
                     ! parcel SL
     &, ql_parc(row_length*rows,model_levels)                          &
                     ! parcel water
     &, t_parc_dil(row_length*rows,model_levels)                       &
                     ! dilute parcel temeperature
     &, th_ref_dil(row_length*rows)                                    &
                     ! dilute parcel reference potential temperature
     &, th_par_km_dil(row_length*rows)                                 &
                     ! dilute parcel reference potential temperature 2nd
     &, t_ref_dil(row_length*rows)                                     &
                     ! dilute parcel reference temperature 
     &, qsat_lev_dil(row_length*rows)                                  &
                     ! qsat for dilute parcel
     &, t_dens_parc_dil(row_length*rows,model_levels)                  &
                     ! dilute parcel t_dens
     &, ql_parc_dil(row_length*rows,model_levels)                      
                     ! dilute parcel liquid water


     logical ::     &
     & l_keep_water & ! if true keeps water loading in plume
                      ! false removed if water exceeds 1g/kg
     &,l_wtest      & ! do w test
     &,l_tv1sd_new    ! Use new adiabatic parcel "perturbation" 
!
! Model constants:
!

! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_A----------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Replace variable A by more meaningful name for
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Convert to Fixed/Free format. P. Selwood

      ! Mean radius of Earth in m.
      Real, Parameter  :: Earth_Radius = 6371000.

!*----------------------------------------------------------------------
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
!!----------------------------------------------------------------------
!!!-----------COMDECK C_ROUGH FOR SUBROUTINE SF_EXCH----------
! Sea ice parameters
! Z0FSEA = roughness length for free convective heat and moisture
!          transport over the sea (m).
!          DUMMY VARIABLE - Only used in 7A boundary layer scheme
! Z0HSEA = roughness length for heat and moisture transport
!          over the sea (m).
! Z0MIZ  = roughness length for heat, moisture and momentum over
!          the Marginal Ice Zone (m).
! Z0SICE = roughness length for heat, moisture and momentum over
!          sea-ice (m).
      REAL :: Z0HSEA
      REAL :: Z0MIZ
      REAL :: Z0SICE

      COMMON  /RUN_BLICE/Z0HSEA,Z0MIZ,Z0SICE

!!----------------------------------------------------------------------
! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end
! Description:
!  Parameters for the calculation of T at the lifting condensation level.
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
! See Bolton 1980 Mon Wea Rev 108, P1046-1053 for details of empirical 
! relationship

      REAL,PARAMETER::                                                  &
     &  a_bolton = 55.0                                                 &
     &, b_bolton = 2840.0                                               &
     &, c_bolton = 3.5                                                  &
     &, d_bolton = 4.805
      
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end
! ENTCNST start

!  6.2   28/02/05  Modifications for STPH_RP
      ! coefficients used in calculation of entrainment rate
      REAL,PARAMETER:: AE1     = 1.0
      REAL,PARAMETER:: AE2     = 1.5
      REAL:: ENTCOEF
      COMMON/STPH_RP_1/entcoef
      REAL,PARAMETER:: SH_FAC  = 1.0

! ENTCNST end
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

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

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end

      real, parameter ::                                                &
     & GRCP=G/CP                                                        &
                                ! Adiabatic lapse rate.
     &,LCRCP=LC/CP                                                      &
                                ! Latent heat of condensation / CP.
     &,LFRCP=LF/CP                                                      &
                                ! Latent heat of fusion / CP.
     &,LS=LC+LF                                                         &
                                ! Latent heat of sublimation.
     &,LSRCP=LS/CP                                                      &
                                ! Latent heat of sublimation / CP.
     &,ra2=1./(Earth_Radius*Earth_Radius)  ! 1/(a*a)

      real, parameter ::      &
       a_plume = 0.2          & ! Minimum initial parcel dtheta
      ,b_plume = 3.26         & ! Used in initial parcel perturbation cal.
      ,max_t_grad = 1.0E-3    & ! Used in initial parcel perturbation cal.
      ,SC_CFTOL   = 0.1         ! Cloud fraction required for a
                                ! cloud layer to be diagnosed.

!-----------------------------------------------------------------------
! Mixing ratio, r,  versus specific humidity, q
!
! In most cases the expression to first order are the same
!
!  Tl = T - (lc/cp)qcl - [(lc+lf)/cp]qcf
!  Tl = T - (lc/cp)rcl - [(lc+lf)/cp]rcf  - equally correct definition
!
! thetav = theta(1+cvq)         accurate
!        = theta(1+r/epsilon)/(1+r) ~ theta(1+cvr) approximate
! 
! svl = (Tl+gz/cp)*(1+(1/epsilon-1)qt)   
!     ~ (Tl+gz/cp)*(1+(1/epsilon-1)rt)
!
! dqsat/dT = epsilon*Lc*qsat/(R*T*T) 
! drsat/dT = epsilon*Lc*rsat/(R*T*T)  equally approximate
!
! Only altering the expression for vapour pressure
!
!  e = qp/epsilon       - approximation 
!  e = rp/(epsilon+r)   - accurate    
!-----------------------------------------------------------------------
! 1.0 Initialisation
!-----------------------------------------------------------------------

      ERROR = 0

      l_tv1sd_new   = .false.


!-----------------------------------------------------------------------
! 1.1 Verify grid/subset definitions.
!-----------------------------------------------------------------------

      If ( bl_levels <  1 .OR. ROWS <  1 .OR. MODEL_levels <  1 .OR.    &
     &      wet_model_levels <  1 ) then
        ERROR = 1
        goto 9999

      End If

!-----------------------------------------------------------------------
! 1.2 initialisation of output arrays
!-----------------------------------------------------------------------
      Do j=1,rows
        Do i=1,row_length

          cumulus(i,j)     = .false.
          L_shallow(i,j)   = .false.
          L_congestus(i,j) = .false.
          L_congestus2(i,j) = .false.
          ntml(i,j) = 1
          nlcl(i,j) = 1
          delthvu(i,j)  = 0.0
          ql_ad(i,j)    = 0.0
          CAPE(i,j)     = 0.0
          CIN(i,j)      = 0.0
          z_lcl(i,j) = z_half(i,j,nlcl(i,j)+1)

! Set LCL for UV grid to level below z_lcl (note that (at vn5.1) UV
! levels are half levels (below) relative to P-grid. Consistent with
! BL scheme's treatment of staggered vertical grid.)

          z_lcl_uv(i,j)= z_full(i,j,nlcl(i,j))

        End Do
      End Do




!-----------------------------------------------------------------------
! set variables
!-----------------------------------------------------------------------
!  Set MBL, "maximum number of boundary levels" for the purposes of
!  boundary layer height calculation.

      MBL = bl_levels - 1

!-----------------------------------------------------------------------
! 1.3 Put level 1 wind components on P grid
!-----------------------------------------------------------------------

! DEPENDS ON: u_to_p
        Call U_TO_P(U,row_length,rows,1,                                &
     &            off_x, off_y, model_domain,                           &
     &            at_extremity, u_p)

! DEPENDS ON: v_to_p
        Call V_TO_P(V,row_length,rows,n_rows,1,                         &
     &            off_x, off_y, model_domain,                           &
     &            at_extremity, v_p)

      If(model_domain  ==  mt_global) then

! Overwrite values of U_P, V_P at the poles with the magnitude of
! the vector wind.

! DEPENDS ON: polar_vector_wind_n
        Call Polar_vector_wind_n(                                       &
     &                       v,                                         &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, 1, mag_vector_np,                  &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       off_x, off_y, global_row_length,           &
     &                       proc_row_group, at_extremity)
        If (at_extremity(PSouth) ) Then
            DO I=1,row_length
              V_P(i,1) = mag_vector_sp(1)
              U_P(i,1) = 0.0
            End Do
        End If
        If (at_extremity(PNorth) ) Then
            DO I=1,row_length
              V_P(i,rows) = mag_vector_np(1)
              U_P(i,rows) = 0.0
            End Do
        End If

      End If

!-----------------------------------------------------------------------
! 1.4 Set appropriate roughness length
!-----------------------------------------------------------------------
      ! Default z0h for standard tstar boundary condition
      DO J=1, rows
      DO I=1, row_length
        If (land_mask(I,j)) then
!         ! Approximate z0h for land as 0.1.
          z0h(i,j) = 0.1
        Else
          z0h(i,j) = z0hsea
        End If
      End Do ! I
      End Do ! J

      If ( L_spec_z0 ) then
        ! Code to use z0h_scm if Namelist specifies it
        DO J=1, rows
        DO I=1, row_length
          If ( z0h_SCM(i,j)  >   0.0 ) then
            z0h(i,j) = z0h_SCM(i,j)
          End If ! z0h_scm
        End Do ! I
        End Do ! J
      End If

!-----------------------------------------------------------------------
! 1.5 Surface buoyancy flux and calculation of unstable points
!-----------------------------------------------------------------------

      IF ( .NOT. L_flux_bc) THEN        ! used by most UM runs

! DEPENDS ON: qsat_mix
       Call qsat_mix(qs_star,tstar,pstar,row_length*rows,l_mixing_ratio)

       ii = 0
       Do j=1,rows
       Do i=1,row_length

!-----------------------------------------------------------------------
! Calculate the surface buoyancy flux
! Uses approximation for unstable CD as 1.5*neutral value (defined
! Garratt p54) and that CH=CD. Approximate Z0H for land as 0.1.
!-----------------------------------------------------------------------

         theta1 = theta(i,j,1)
         rhostar = pstar(i,j) / ( R*tstar(i,j) )

         Ushear = U_P(i,j) - U_0_P(i,j)
         Vshear = V_P(i,j) - V_0_P(i,j)
         wshr2 = MAX (1.0E-6 , Ushear*Ushear + Vshear*Vshear)
         wshr1 = SQRT(wshr2)
         cd  = 1.5 * ( vkman/LOG(z_full(i,j,1)/z0h(i,j)) )**2

         If (land_mask(i,j)) then         ! land
            wthvbar = wshr1 * CD                                          &
     &       * ( TSTAR(I,j)*((100000.0/PSTAR(I,j))**kappa) - theta1 )
         Else
            wthvbar = wshr1 * CD                                          &
     &       * ( TSTAR(I,j)*((100000.0/PSTAR(I,j))**kappa) - theta1     &
     &         + 0.61*theta1*(QS_STAR(I,j)-Q(I,j,1)) )
         End If

         ustar = SQRT(cd  * wshr2)
         fb_surf(i,j) = G * wthvbar /                                   &
     &                       ( rhostar * theta1*(1.0+0.61*Q(i,j,1)) )

         If (fb_surf(i,j)  >   0.0) then
          ii= ii+1
          tv1_sd(ii) = 1.93 * wthvbar / ( rhostar * ( 75.0 *           &
     &            fb_surf(i,j) + ustar*ustar*ustar)**(1.0/3.0) )
         End If

       End Do ! I
       End Do ! J

      Else ! if L_flux_bc       (used by some SCM runs)

!       !------------------------------------------
!       ! If specified surface fluxes are required
!       !------------------------------------------
        Do J=1, rows
        Do I=1, row_length
          ! For taylor expansion about T0=SL(K=1)
          tstar(I,J) = theta(i,j,1) * exner_theta_levels(i,j,1)         &
     &                                        +GRCP*z_full(I,J,1)
        End Do
        End Do

! DEPENDS ON: qsat_mix
        Call qsat_mix(qs_star,tstar,pstar,row_length*rows,l_mixing_ratio)

! dqsat/dT - approximation same expression whether specific humidity or 
!            mixing ratio.

        DO J=1, rows
        DO I=1, row_length
          dqsdt(I,J) = (epsilon * LC * qs_star(I,J))                    &
     &                 / ( R * tstar(I,J) * tstar(I,J) )
        End Do
        End Do

       If (icvdiag < 2 ) then

        ii = 0
        Do j=1,rows
        Do i=1,row_length

         Ushear = U_P(I,j) - U_0_P(I,j)
         Vshear = V_P(I,j) - V_0_P(I,j)
!        ! Need to have a higher minimum wind speed limit with
!        ! specified fluxes in order not to generate huge TSTAR
         wshr2 = MAX (0.1, Ushear*Ushear + Vshear*Vshear)
         wshr1 = SQRT(wshr2)

         ! Calculate WTHV from namelist flux_h and flux_e (in W/m2)
         wthvbar = ((flux_h(I,j)/CP)+0.61*(flux_e(I,j)/LC))             &
     &           * ((100000.0/PSTAR(I,j))**kappa)

         cd  = 1.5 * ( vkman/LOG(z_full(I,j,1)/z0h(i,j)) )**2

         theta1 = theta(I,j,1)

         ! Taylor expansion for qsat(T*) about SL(k=1)
         TSTAR(I,j) = ( theta1 + (wthvbar/(wshr1*CD))                   &
     &              -   0.61*theta1                                     &
     &              *   (QS_STAR(I,j)-Q(I,j,1)-DQSDT(I,j)*TSTAR(I,j)))  &
     &              /   ( (100000.0/PSTAR(I,j))**kappa +                &
     &                       0.61*theta1*DQSDT(I,j) )

         RHOSTAR = PSTAR(I,j) / ( R*TSTAR(I,j) )

         USTAR = SQRT(CD * wshr2)
         FB_SURF(I,j) = G * wthvbar /                                   &
     &                       ( RHOSTAR * theta1*(1.0+0.61*Q(I,j,1)) )


         IF (FB_SURF(I,j)  >   0.0) THEN
           ii = ii+1
           tv1_sd(ii) = 1.93 * wthvbar / ( rhostar * ( 75.0 *          &
     &                   fb_surf(I,j) + ustar*ustar*ustar)**(1.0/3.0) )
         End If

        End Do ! I
        End Do ! J
       
       Else    ! (icvdiag >=2 )

        ii=0
        Do j=1,rows
        Do i=1,row_length

         Ushear = U_P(I,j) - U_0_P(I,j)
         Vshear = V_P(I,j) - V_0_P(I,j)
!        ! Need to have a higher minimum wind speed limit with
!        ! specified fluxes in order not to generate huge TSTAR
         wshr2 = MAX (0.1, Ushear*Ushear + Vshear*Vshear)
         wshr1 = SQRT(wshr2)

         ! Calculate WTHV from namelist flux_h and flux_e (in W/m2)
         wthvbar = ((flux_h(I,j)/CP)+0.61*(flux_e(I,j)/LC))             &
     &           * ((100000.0/PSTAR(I,j))**kappa)

         CD = 1.5 * ( VKMAN/log(z_full(I,j,1)/Z0H(i,j)) )**2

         theta1 = theta(I,j,1)
! new 
!         WSTAR = ( ZH(I,j) * G * 
!     &             WTHVBAR / ( THETA1*(1.0+0.61*Q(I,j,1)) ) )**(1./3.)

         !  Use "effective" ustar, allowing for convective eddies 
!         USTAR_EFF = SQRT( CD*WSHR2 + BETA*BETA*WSTAR*WSTAR )

         ! Taylor expansion for qsat(T*) about SL(k=1)
!new         TSTAR(I,j) = ( theta1 + (wthvbar/(ustar_eff1*CD))          &
         TSTAR(I,j) = ( theta1 + (wthvbar/(wshr1*CD))                   &
     &              -   0.61*theta1                                     &
     &              *   (QS_STAR(I,j)-Q(I,j,1)-DQSDT(I,j)*TSTAR(I,j)))  &
     &              /   ( (100000.0/PSTAR(I,j))**kappa +                &
     &                       0.61*theta1*DQSDT(I,j) )

         RHOSTAR = PSTAR(I,j) / ( R*TSTAR(I,j) )

         USTAR = SQRT(CD * wshr2)
         FB_SURF(I,j) = G * wthvbar /                                   &
     &                       ( RHOSTAR * theta1*(1.0+0.61*Q(I,j,1)) )

         IF (FB_SURF(I,j)  >   0.0) THEN
          ii=ii+1
           If (l_tv1sd_new) then 
             w_s_cubed = 0.25 * zh(i,j) * fb_surf(i,j)
! Std: W_S_CUBED = 75.0*FB_SURF(I,j), ie. ZH=300.
             w_m       = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
             tv1_sd(ii) = 1.93 * wthvbar/( rhostar * w_m )           
 
           Else
! old method assumes W_S_CUBED = 75.0*FB_SURF(I,j), ie. ZH=300.
             TV1_SD(ii) = 1.93 * wthvbar / ( RHOSTAR * ( 75.0 *        &
     &                   FB_SURF(I,j) + USTAR*USTAR*USTAR)**(1.0/3.0) )
           End If
         End If

        End Do ! I
        End Do ! J
  
       End If  ! test on icvdiag

      End If !  L_flux_bc

!-----------------------------------------------------------------------
! 2.0 Decide on unstable points  ( fb_surf > 0.0 )
!     Only work on these points for the rest of the calculations.
!-----------------------------------------------------------------------

      nunstable = 0           ! total number of unstable points
      Do j=1,rows
      Do i=1,row_length
        If ( fb_surf(i,j)  >   0.0 ) then
          nunstable = nunstable + 1
          index_i(nunstable) = i
          index_j(nunstable) = j
        End If 
      End Do
      End Do

!-----------------------------------------------------------------------
! 2.1 initialise just unstable points
!-----------------------------------------------------------------------
      Do ii=1,nunstable

          shmin(ii)   = .false.
          topbl(ii)   = .false.
          topprof(ii) = .false.
          kshmin(ii)     = 1
          kcucheck(ii)   = 1
          k_max(ii)      = 1
          k_neutral(ii)  = 1
          k_inv(ii)      = 1
          freeze_lev(ii) = 1
          dtv_min(ii)   = 0.0
          dtv_min2(ii)  = 0.0
          delthvu2(ii)  = 0.0
          max_buoy(ii)  = 0.0
          ntml_c(ii)    = 1
          nlcl_c(ii)    = 1
          delthvu_c(ii) = 0.0
          cape_c (ii)   = 0.0
          cin_c  (ii)   = 0.0

      End Do

! compress boundary layer depth

      Do ii=1,nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        zh_c(ii)  = zh(i,j)
        z_lcl_c(ii) = z_half(i,j,nlcl_c(ii)+1)
      End Do



!-----------------------------------------------------------------------
! 2.1 Calculate mass of layers 
!-----------------------------------------------------------------------
      Do K=1, model_levels-1
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          dmass_theta(ii,k) = rho_theta(i,j,k)*                         &
     &                            (r_theta_levels(i,j,k)/earth_radius)  &
     &                           *(r_theta_levels(i,j,k)/earth_radius)  &
     &              *(r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))
        End Do
      End Do


!-----------------------------------------------------------------------
! 3.0 Calculate various quantities required by the parcel calculations
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 3.1 Section BL.1 Calculate T at old time level.
!-----------------------------------------------------------------------

      Do k = 1, model_levels
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          T(ii,k) = theta(i,j,k) * exner_theta_levels(i,j,k)

            ! initialise t_parc at all points
            ! added for safety of qsat_mix calls later

          t_parc(ii,k) = t(ii,k)
          q_c(ii,k)   = q(i,J,K)
          qcl_c(ii,k) = qcl(i,J,K)
          qcf_c(ii,k) = qcf(i,J,K)
          z_full_c(ii,k) = z_full(i,j,k)  
          z_half_c(ii,k) = z_half(i,j,k)  
          exner_theta_levels_c(ii,k) = exner_theta_levels(i,j,k)
          P_theta_lev_c(ii,k) = P_theta_lev(i,j,k)
          r_theta_levels_c(ii,k) =r_theta_levels (i,j,k)
        End Do
      End Do

!-----------------------------------------------------------------------
! 3.2 Calculate total water content, QW and Liquid water temperature, TL
!     Definitions for Qw and Tl the same whether q, qcl, qcf  are 
!     specific humidities  or mixing ratio.
!     
!-----------------------------------------------------------------------

      DO K=1,wet_model_levels
        Do ii=1, nunstable

          QW(ii,K) = q_c(ii,K) + qcl_c(ii,k) + qcf_c(ii,k)
                                ! BL doc P243.10
          TL(ii,K) = T(ii,K) - LCRCP*qcl_c(ii,k) - LSRCP*qcf_c(ii,k)
                                ! BL doc P243.9

!
! Calculate SVL: conserved variable  a form of moist static energy /cp 
!       svl = (Tl+gz/cp)*(1+(1/epsilon-1)qt)  - specific humidity
!

          SVL(ii,K) = ( TL(ii,K) + GRCP * z_full_c(ii,K) )             &
     &                               * ( 1.0 + c_virtual*QW(ii,K) )

        End Do
      End Do



!-----------------------------------------------------------------------
! 4.0 Parcel calculation - various options available 
!-----------------------------------------------------------------------
!  icvdiag option         explaination
!     0              original undilute parcel calculation    
!     1              undilute parcel calculation (improved version).
!     2 and above    dilute parcel calculation
!        2           1/z entrainment rates ?
!        3            0.55/z
!=======================================================================
! Choice of diagnosis method
!=======================================================================
      IF (ICVDIAG == 0) THEN   ! old code Still required at present
!=======================================================================
!-----------------------------------------------------------------------
!! 0.1 Set TOPBL to .FALSE. and calculate boundary layer top using
!!     a non-local, plume method.
! ----------------------------------------------------------------------
!
      Do ii=1, nunstable

! Initialise CAPE and CIN
        CAPE_c(ii) = 0.0
        CIN_c(ii)  = 0.0
        
        k_plume(ii) = 1

!-----------------------------------------------------------------------
! Only perform parcel ascent if unstable
! Start plume ascent from grid-level above top of surface layer, taken
! to be at a height, z_surf, given by 0.1*ZH
!-----------------------------------------------------------------------
        z_surf = 0.1 * zh_c(ii)

        Do while( z_full_c(ii,k_plume(ii))  <   z_surf .and.           &
!                   ! not reached z_surf
     &              SVL(ii,k_plume(ii)+1)  <   SVL(ii,k_plume(ii)) )
!                   ! not reached inversion

            k_plume(ii) = k_plume(ii) + 1

        End Do
      End Do       ! ii loop 

 

      Do ii=1, nunstable
        sl_plume(ii) = TL(ii,k_plume(ii))                              &
     &                        + GRCP * z_full_c(ii,k_plume(ii))
        thv_pert(ii) = max( a_plume,                                   &
     &                  min( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) )
        qw_plume(ii) = QW(ii,k_plume(ii))
      End Do      ! ii loop 

!-----------------------------------------------------------------------
!! 0.2 Calculate temperature and pressure of lifting condensation level
!!     using approximations from Bolton (1980)
!-----------------------------------------------------------------------
!
!   vapour pressure e ~ qp/epsilon   q specific humidity
!   vapour pressure e ~ qp/(epsilon+q)   q mixing ratio

      If (l_mixing_ratio) then
        Do ii=1, nunstable
! expression for mixing ratio
          vap_press = 0.01*Q_c(ii,k_plume(ii)) *                       &
     &                                 P_theta_lev_c(ii,k_plume(ii))   &
     &                      / (EPSILON+Q_c(ii,k_plume(ii)) )
          If (vap_press  >   0.0) then
           T_LCL(ii) = a_bolton + b_bolton/                            &
     &                      (c_bolton*LOG(T(ii,k_plume(ii)))           &
     &                                    - LOG(vap_press) - d_bolton )

           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
     &                ( T_LCL(ii) / T(ii,k_plume(ii)) )**(1.0/kappa)
          Else
           i = index_i(ii)   
           j = index_j(ii)   
           P_LCL(ii) = pstar(i,j)
          End If
 
        End Do

      Else       ! not l_mixing_ratio
! expression for specific humidity
        Do ii=1, nunstable
          vap_press = Q_c(ii,k_plume(ii)) *                            &
     &         P_theta_lev_c(ii,k_plume(ii)) / ( 100.0*EPSILON )
          If (vap_press  >   0.0) then
           T_LCL(ii) = a_bolton + b_bolton/                            &
     &                      (c_bolton*LOG(T(ii,k_plume(ii)))           &
     &                                    - LOG(vap_press) - d_bolton )

           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
     &                ( T_LCL(ii) / T(ii,k_plume(ii)) )**(1.0/kappa)
          Else
            i = index_i(ii)   
            j = index_j(ii)   
            P_LCL(ii) = pstar(i,j)
          End If
 
        End Do

      End If ! test on l_mixing_ratio  

!-----------------------------------------------------------------------
! calculate parcel water by linearising qsat about the environmental
! temperature
! Note: the following calculation is for parcel in level 1:
!
        k = 1 
! DEPENDS ON: qsat_mix
        Call QSAT_MIX(QS(1,K),T(1,K),P_theta_lev_c(1,K),               &
     &                                     nunstable,l_mixing_ratio)

        Do ii=1, nunstable

         IF(T(ii,K) >  TM) THEN
           DQ_SAT_ENV = epsilon*LC*QS(ii,K)/(R*T(ii,K)**2.0)

           Q_LIQ_PARC = MAX( 0.0, ( QW_plume(ii) - QS(ii,K) -          &
     &       DQ_SAT_ENV*( SL_plume(ii)-GRCP*z_full_c(ii,K)-T(ii,K) )   &
     &                             ) / (1.0+LCRCP*DQ_SAT_ENV) )

           Q_LIQ_ENV = MAX( 0.0, ( QW(ii,K) - QS(ii,K) -DQ_SAT_ENV*    &
     &        ( TL(ii,K) - T(ii,K) ) ) / (1.0+LCRCP*DQ_SAT_ENV) )
!
! add on the difference in the environment's ql as calculated by the
! UM cloud scheme (using some RH_CRIT value) and what it
! would be if RH_CRIT=1. This then imitates partial condensation
! in the parcel.
!
           Q_LIQ_PARC = Q_LIQ_PARC + QCL_c(ii,K) + QCF_c(ii,K)         &
     &                     - Q_LIQ_ENV
           T_PARC(ii,k)=                                               &
     &            SL_plume(ii)-GRCP*z_full_c(ii,K)+LCRCP*Q_LIQ_PARC
         Else
           DQ_SAT_ENV=epsilon*LS*QS(ii,K)/(R*T(ii,K)**2.0)

           Q_LIQ_PARC = MAX( 0.0, ( QW_plume(ii) - QS(ii,K) -          &
     &       DQ_SAT_ENV*( SL_plume(ii)-GRCP*z_full_c(ii,K)-T(ii,K) )   &
     &                             ) / (1.0+LSRCP*DQ_SAT_ENV) )

           Q_LIQ_ENV = MAX( 0.0, ( QW(ii,K) - QS(ii,K) -DQ_SAT_ENV*    &
     &         ( TL(ii,K) - T(ii,K) ) ) / (1.0+LSRCP*DQ_SAT_ENV) )

! add on difference in environment's ql between RH_CRIT and RH_CRIT=1

           Q_LIQ_PARC = Q_LIQ_PARC + QCL_c(ii,K) + QCF_c(ii,K)         &
     &                     - Q_LIQ_ENV

           T_PARC(ii,k) =                                              &
     &            SL_plume(ii)-GRCP*z_full_c(ii,K)+LSRCP*Q_LIQ_PARC
         End If     ! test on TM

         Q_VAP_PARC=QW_plume(ii)-Q_LIQ_PARC

         T_DENS_PARC(ii,k) =                                           &
     &         T_PARC(ii,k) *(1.0+C_VIRTUAL*Q_VAP_PARC-Q_LIQ_PARC)
         T_DENS_ENV(ii,k)=T(ii,K)*                                     &
     &                (1.0+C_VIRTUAL*Q_c(ii,K)-QCL_c(ii,K)-QCF_c(ii,K))

         ENV_SVL_KM1(ii) = T_DENS_ENV(ii,k)+GRCP*(z_full_c(ii,K))
         PAR_SVL_KM1(ii) = T_DENS_PARC(ii,k)+GRCP*(z_full_c(ii,K))
       End Do         ! ii loop


!
! Reset zh  (at this point in the code ntml is initialised as =1)
!
      do j=1,rows
        do i=1,row_length
          zh(i,j) = z_half(i,j,ntml(i,j)+1)
        End Do
      End Do
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        zh_c(ii) = zh(i,j)   
      End Do
!
!-----------------------------------------------------------------------
! Find NLCL
!-----------------------------------------------------------------------
      DO  K = 2,wet_model_levels

        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
           IF ( P_LCL(ii)  <   P(I,j,K) ) THEN
             NLCL(I,j) = K-1
             NLCL_c(ii) = K-1
             Z_LCL(I,j) = Z_HALF_c(ii,NLCL(I,j)+1)
             Z_LCL_UV(I,j) = z_full_c(ii,NLCL(I,j))
           End If
        End Do
      End Do

!-----------------------------------------------------------------------
!! 0.3 Now compare plume s_VL with each model layer s_VL in turn to
!!     find the first time that plume has negative buoyancy.
!-----------------------------------------------------------------------
!
      DO  K = 2,wet_model_levels

! DEPENDS ON: qsat_mix
        Call QSAT_MIX(QS(1,K),T(1,K),P_theta_lev_c(1,K),               &
     &                                       nunstable,l_mixing_ratio)

!CDIR nodep
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   

!-----------------------------------------------------------------------
! Only perform parcel ascent if unstable
!-----------------------------------------------------------------------

!  Find level just below 2.5km (for use in Cu diagnosis)

          If ( z_full_c(ii,K)  >   2500.0                              &
     &         .and. kcucheck(ii)  ==  1 ) kcucheck(ii) = K-1
!
!-----------------------------------------------------------------------
! Set flag to true when level BELOW is above the lcl (in case P_LCL is
! in the upper half of layer NLCL+1). Implies ABOVE_LCL at NLCL+3.
!-----------------------------------------------------------------------
          IF ( k > 2 .and. P_LCL(ii)  <   P(I,j,K-1) ) THEN
            ABOVE_LCL(ii)=.FALSE.
          Else IF ( k == 2 .and. P_LCL(ii)  <   Pstar(I,j) ) THEN
            ABOVE_LCL(ii)=.FALSE.
          Else
            ABOVE_LCL(ii)=.TRUE.
          End If
!
!         !-----------------------------------------------------------
!         ! calculate parcel potential temperature by linearising
!         ! q_sat about the environmental temperature.
!         !-----------------------------------------------------------
!
          IF(T(ii,K) >  TM) THEN
            DQ_SAT_ENV=epsilon*LC*QS(ii,K)/(R*T(ii,K)**2.0)
            Q_LIQ_PARC = MAX( 0.0, ( QW_plume(ii) - QS(ii,K) -         &
     &        DQ_SAT_ENV*( SL_plume(ii)-GRCP*z_full_c(ii,K)-T(ii,K) )  &
     &                             ) / (1.0+LCRCP*DQ_SAT_ENV) )
            Q_LIQ_ENV = MAX( 0.0, ( QW(ii,K) - QS(ii,K) -DQ_SAT_ENV*   &
     &         ( TL(ii,K) - T(ii,K) ) ) / (1.0+LCRCP*DQ_SAT_ENV) )
! add on the difference in the environment's ql as calculated by the
! partial condensation scheme (using some RH_CRIT value) and what it
! would be if RH_CRIT=1. This then imitates partial condensation
! in the parcel.
            Q_LIQ_PARC = Q_LIQ_PARC + QCL_c(ii,K) + QCF_c(ii,K)        &
     &                     - Q_LIQ_ENV
            T_PARC(ii,k) =                                             &
     &             SL_plume(ii)-GRCP*z_full_c(ii,K)+LCRCP*Q_LIQ_PARC
          Else       
            DQ_SAT_ENV=epsilon*LS*QS(ii,K)/(R*T(ii,K)**2.0)
            Q_LIQ_PARC = MAX( 0.0, ( QW_plume(ii) - QS(ii,K) -         &
     &        DQ_SAT_ENV*( SL_plume(ii)-GRCP*z_full_c(ii,K)-T(ii,K) )  &
     &                             ) / (1.0+LSRCP*DQ_SAT_ENV) )
            Q_LIQ_ENV = MAX( 0.0, ( QW(ii,K) - QS(ii,K) -DQ_SAT_ENV*   &
     &         ( TL(ii,K) - T(ii,K) ) ) / (1.0+LSRCP*DQ_SAT_ENV) )
! add on difference in environment's ql between RH_CRIT and RH_CRIT=1
            Q_LIQ_PARC = Q_LIQ_PARC + QCL_c(ii,K) + QCF_c(ii,K)        &
     &                     - Q_LIQ_ENV
            T_PARC(ii,k) =                                             &
     &             SL_plume(ii)-GRCP*z_full_c(ii,K)+LSRCP*Q_LIQ_PARC
          End If       ! test on TM
          Q_VAP_PARC=QW_plume(ii)-Q_LIQ_PARC
!
          T_DENS_PARC(ii,k)=T_PARC(ii,k)*                              &
     &                       (1.0+C_VIRTUAL*Q_VAP_PARC-Q_LIQ_PARC)
          T_DENS_ENV(ii,k)=T(ii,K)*                                    &
     &               (1.0+C_VIRTUAL*Q_c(ii,K)-QCL_c(ii,K)-QCF_C(ii,K))

!         !-------------------------------------------------------------
!         ! Find vertical gradients in parcel and environment SVL
!         ! (using values from level below (i.e. K-1)).
!         !-------------------------------------------------------------
          DPAR_BYDZ2 = (T_DENS_PARC(ii,k) + GRCP*z_full_c(ii,K) -      &
     &                 PAR_SVL_KM1(ii)) /                              &
     &              (z_full_c(ii,K) - z_full_c(ii,K-1))
          DENV_BYDZ2 = (T_DENS_ENV(ii,k) + GRCP*z_full_c(ii,K) -       &
     &                 ENV_SVL_KM1(ii)) /                              &
     &              (z_full_c(ii,K) - z_full_c(ii,K-1))

          IF ( TOPBL(ii) .AND. ( DENV_BYDZ2  <   DPAR_BYDZ2 .OR.       &
     &               PAR_SVL_KM1(ii) - ENV_SVL_KM1(ii)  <=  0.0 )      &
     &                              .AND. .NOT. SHMIN(ii) ) THEN
              SHMIN(ii) = .TRUE.
              DT_DENS_PARC_TMIN(ii) = PAR_SVL_KM1(ii) -                &
     &                                        ENV_SVL_KM1(ii)
              KSHMIN(ii) = K-1
          End If
          IF ( .NOT.TOPBL(ii) .AND. K  >   K_plume(ii) .AND.           &
     &  (  ( T_DENS_PARC(ii,k)-T_DENS_ENV(ii,k)  <=  - THV_PERT(ii) )  &
!
!                      plume non buoyant
!
     &    .OR. (ABOVE_LCL(ii) .AND. (DENV_BYDZ2  >   1.25*DPAR_BYDZ2)) &
!
!                      or environmental virtual temperature gradient
!                      significantly larger than parcel gradient
!                      above lifting condensation level
!
     &           .OR. (K  >   wet_model_levels-1)                      &
!                      or reached top of model
     &         )                                                       &
     &         ) THEN
!
            TOPBL(ii) = .TRUE.
            ZH(I,j) = Z_HALF_c(ii,K)
            NTML(i,j) = K-1
            DT_DENS_PARC_T(ii) = PAR_SVL_KM1(ii) - ENV_SVL_KM1(ii)
            IF ( DELTHVU(I,j)  >   0.0) THEN
! compensate for any negative buoyancy of parcel in cloud layer
              DELTHVU(I,j) = DELTHVU(I,j) - DTV_MIN(ii) *              &
     &                                    ( ZH(I,j) - Z_LCL(I,j) )
            End If

! Added to estimate undilute CAPE and CIN. Note for this option the CAPE
! may not always be that of an ascent reaching the level of neutral buoyancy.

            inc = g * (t_dens_parc(ii,k) - t_dens_env(ii,k))           &
                  *(z_half_c(ii,k+1) - z_half_c(ii,k))/t_dens_env(ii,k)

            If (inc < 0.0) Then
              CIN_c(ii) = CIN_c(ii) + inc
            End If          

            CAPE_c(ii) = CAPE_c(ii) + inc
           
          End If           ! test on topbl

          ENV_SVL_KM1(ii) = T_DENS_ENV(ii,k) + GRCP*z_full_c(ii,K)
          PAR_SVL_KM1(ii) = T_DENS_PARC(ii,k) + GRCP*z_full_C(ii,K)
!
          IF (K  >   NLCL(I,j) .AND. .NOT. TOPBL(ii)) THEN
            DTV_MIN(ii) = MIN( DTV_MIN(ii),                            &
     &                     (T_DENS_PARC(ii,k)-T_DENS_ENV(ii,k)) *      &
     &                     ( (100000.0/P_theta_lev_c(ii,K))**kappa ) )
            DELTHVU(I,j) = DELTHVU(I,j) + (T_DENS_PARC(ii,k)           &
     &                             -T_DENS_ENV(ii,k)) *                &
     &                     ( (100000.0/P_theta_lev_c(ii,K))**kappa ) * &
     &                     ( Z_HALF_c(ii,K+1) - Z_HALF_c(ii,K) )
          End If
!

        End Do   ! ii loop
      End Do     ! k loop
!-----------------------------------------------------------------------
!! 0.4 Save parcel ascent top: this will be used to allow mixing and
!!     entrainment into decoupled Sc of single layer thickness when it
!!     occurs above Cu.
!-----------------------------------------------------------------------
      do j=1,rows
      do i=1,row_length
        ZHPAR(I,j) = ZH(I,j)
        NTPAR(I,j) = NTML(I,j)
      End Do
      End Do

!-----------------------------------------------------------------------
!     Test height derived above against lifting condensation level
!-----------------------------------------------------------------------
      Do ii=1,nunstable
          i = index_i(ii)   
          j = index_j(ii)   
!-----------------------------------------------------------------------
!     Check lifting condensation levels against height of parcel ascent,
!     if lifting condensation level lower than parcel ascent, and is
!     within BL_LEVELS, then decide
!     on type of cloudy layer. If lifting condensation level at or below
!     low grid point, assume fog layer and turbulent mixing. For
!     gradient tests assume any if LCL and top of parcel ascent is less
!     than two levels then stratocumulus.
!-----------------------------------------------------------------------
        IF ( NTML(I,j)-NLCL(I,j)  >=  2                                &
     &                        .AND. NLCL(I,j)  >   K_plume(ii)         &
     &                              .AND. NLCL(I,j)  <   MBL-1 ) THEN
!-----------------------------------------------------------------------
!     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
!     For stratocumulus top of mixed layer = ZH
!     For cumulus top of mixed layer = ZLCL
!     NTML >= MBL indicates convection.
!     Diagnosis is done by comparing gradients
!-----------------------------------------------------------------------
          IF (NTML(I,j)  >=  MBL) THEN
            CUMULUS(I,j) = .TRUE.
          Else

! TEMPORARY: SHOULD REALLY BE DONE WITH SCALE HEIGHT
            IF (NTML(I,j)  >   KCUCHECK(ii)                            &
     &             .AND. NLCL(I,j)  <=  KCUCHECK(ii)-2) THEN

              GRAD_CLD =  ABS( QW(ii,KCUCHECK(ii)) -                   &
     &                                  QW(ii,NLCL(I,j)) ) /           &
     &          ( z_full_C(ii,KCUCHECK(ii)) - z_full_c(ii,NLCL(I,j)) )
            Else
              GRAD_CLD =  ABS( QW(ii,NTML(I,j)) -                      &
     &                                  QW(ii,NLCL(I,j)) ) /           &
     &              ( z_full_c(ii,NTML(I,j)) - z_full_C(ii,NLCL(I,j)) )
            End If

            GRAD_SUB =  ABS( QW(ii,NLCL(I,j)) -                        &
     &                                QW(ii,K_plume(ii)) ) /           &
     &           ( z_full_c(ii,NLCL(I,j)) - z_full_c(ii,K_plume(ii)) )

            IF (GRAD_CLD  >   1.10*GRAD_SUB) THEN
!-----------------------------------------------------------------------
!     Not well mixed, however it is possible that the depth of a well
!     mixed boundary layer has increased but not yet been mixed yet so
!     test gradient from next level down.
!-----------------------------------------------------------------------

! TEMPORARY: SHOULD REALLY BE DONE WITH SCALE HEIGHT
              IF (NTML(I,j)  <=  KCUCHECK(ii)) THEN
              GRAD_CLD =  ABS( QW(ii,NTML(I,j)-1) - QW(ii,NLCL(I,j)) ) &
     &           /( z_full_c(ii,NTML(I,j)-1) - z_full_c(ii,NLCL(I,j)) )
              End If

              IF ( GRAD_CLD  >   1.10*GRAD_SUB) THEN
!-----------------------------------------------------------------------
!      Diagnose a cumulus layer
!-----------------------------------------------------------------------
                CUMULUS(I,j) = .TRUE.
              End If
            Else    

! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
! and not yet been mixed (so could have been erroneously identified as
! well-mixed)

! First check using level below (recalculate GRAD_SUB)

              IF (NLCL(I,j) - K_plume(ii)  >=  2) THEN
                 GRAD_SUB =  ABS( QW(ii,NLCL(I,j)-1) -                 &
     &                                QW(ii,K_plume(ii)) ) /           &
     &           ( z_full_c(ii,NLCL(I,j)-1) - z_full_c(ii,K_plume(ii)) )

                 IF ( GRAD_CLD  >   1.10*GRAD_SUB) THEN
                   CUMULUS(I,j) =.TRUE.
                 End If

              End If

! If still diagnosing well-mixed, check using level above
! (recalculate GRAD_CLD)

              IF (.NOT. CUMULUS(I,j) ) THEN

               IF (NTML(I,j)  >   KCUCHECK(ii)                         &
     &             .AND. NLCL(I,j)  <=  KCUCHECK(ii)-2) THEN

                GRAD_CLD =  ABS( QW(ii,KCUCHECK(ii)) -                 &
     &                                  QW(ii,NLCL(I,j)+1) ) /         &
     &         ( z_full_c(ii,KCUCHECK(ii)) - z_full_c(ii,NLCL(I,j)+1) )
               Else
                GRAD_CLD =  ABS( QW(ii,NTML(I,j)) -                    &
     &                                  QW(ii,NLCL(I,j)+1) ) /         &
     &            ( z_full_c(ii,NTML(I,j)) - z_full_c(ii,NLCL(I,j)+1) )
               End If

               IF ( GRAD_CLD  >   1.10*GRAD_SUB) THEN
                 CUMULUS(I,j) =.TRUE.
               End If

              End If    ! not cumulus
            End If      ! cloud gradient test
          End If        ! ntml > mbl test
        End If        !  level test

!-----------------------------------------------------------------------
!      Check that a cumulus layer has not been erroneously diagnosed in
!      a deep cloudy region
!-----------------------------------------------------------------------
        K=NLCL(I,j)
        IF ( LAND_MASK(I,j) .AND. CUMULUS(I,j) .AND.                   &
     &                                 NTPAR(I,j)  <   MBL ) THEN
          DO WHILE ( K  <=  NTPAR(I,j) .AND. CLOUD_FRACTION(I,j,K)     &
     &                                            >=  SC_CFTOL )
            K = K + 1
          End Do
          IF (K  ==  NTPAR(I,j)+1) CUMULUS(I,j) = .FALSE.
        End If


        IF ( CUMULUS(I,j) ) THEN

!-----------------------------------------------------------------------
!       If cumulus has been diagnosed, determine whether it is shallow
!       or deep convection
!-----------------------------------------------------------------------
          IF ( SHMIN(ii) ) THEN

            IF ( w_copy(I,j,BL_LEVELS)  <   0.0 .AND.                  &
     &           (z_full_c(ii,NTPAR(I,j))  <=  2500.0 .OR.             &
     &            T(ii,NTPAR(I,j))  >=  TM)                            &
     &       .AND. (z_full_c(ii,KSHMIN(ii)) - z_full_c(ii,NTPAR(I,j))) &
     &          <=  1.25*(ZHPAR(I,j) - Z_LCL(I,j)) .AND.               &
     &         DT_DENS_PARC_TMIN(ii)  <=  0.55*DT_DENS_PARC_T(ii) )    &
     &      THEN

              L_SHALLOW(I,j) = .TRUE.
            End If

          End If

!-----------------------------------------------------------------------
!      Set mixed layer depth to Z_LCL
!-----------------------------------------------------------------------
          IF (P_LCL(ii)  <   (P_theta_lev(I,j,NLCL(I,j)+1))) THEN
!-----------------------------------------------------------------------
!      If LCL is diagnosed in the upper half of the layer set Z_LCL to
!      the height of the upper layer interface
!      (in code above LCL is always set to the lower interface).
!-----------------------------------------------------------------------
             NLCL(I,j) = NLCL(I,j)+1
             Z_LCL(I,j) = Z_HALF(I,j,NLCL(I,j)+1)
             Z_LCL_UV(I,j)= z_full(I,j,NLCL(I,j))
          End If
          ZH(I,j) = Z_LCL(I,j)
          NTML(I,j) = NLCL(I,j)

!      If CUMULUS has been diagnosed but DELTHVU is negative, reset
!      CUMULUS and L_SHALLOW to FALSE but leave ZH and NTML at LCL

            IF (DELTHVU(I,j)  <=  0.0) THEN

              CUMULUS(I,j)   = .FALSE.
              L_SHALLOW(I,j) = .FALSE.

            End If

        Else

!-----------------------------------------------------------------------
!      If not cumulus, reset parameters to within BL_LEVELS
!-----------------------------------------------------------------------
          IF (NTML(I,j)  >   MBL) THEN
            NTML(I,j) = MBL
            NTPAR(I,j) = MBL
            ZH(I,j) = Z_HALF(I,j,MBL+1)
            ZHPAR(I,j) = ZH(I,j)
          End If
          IF (NLCL(I,j)  >   MBL) THEN
              NLCL(I,j) = MBL
              Z_LCL(I,j) = ZH(I,j)
              Z_LCL_UV(I,j)=z_full(I,j,MBL-1)
          End If

        End If       ! cumulus test
      End Do

! Expand CAPE and CIN values back to full arrays
      Do ii=1,nunstable
        i = index_i(ii)
        j = index_j(ii)
        CAPE(i,j) = CAPE_c(ii)
        CIN(i,j)  = CIN_c(ii)
      End Do



!=======================================================================
! Use  new code  to calculate both parcel and ntpar differently
! This code is taken from 5A scheme
!=======================================================================
      Else if (icvdiag == 1) then

      Do ii=1, nunstable

        k_plume(ii) = 1

!-----------------------------------------------------------------------
! Only perform parcel ascent If unstable
! Start plume ascent from grid-level above top of surface layer, taken
! to be at a height, z_surf, given by 0.1*zh
!-----------------------------------------------------------------------
         z_surf = 0.1 * zh_c(ii)

         Do while( z_full_c(ii,k_plume(ii))  <   z_surf .and.           &
!                   ! not reached z_surf
     &              SVL(ii,k_plume(ii)+1)  <   SVL(ii,k_plume(ii)) )
!                   ! not reached inversion

            k_plume(ii) = k_plume(ii) + 1

         End Do
      End Do

      Do ii=1, nunstable

         sl_plume(ii) = TL(ii,k_plume(ii))                              &
     &                        + GRCP * z_full_c(ii,k_plume(ii))
         thv_pert(ii) = max( a_plume,                                   &
     &                  min( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) )
         qw_plume(ii) = QW(ii,k_plume(ii))

! Added for more acturate parcel cal later

         th_ref(ii) = tl(ii,k_plume(ii))                                &
     &                        /exner_theta_levels_c(ii,k_plume(ii))
         th_par_km1(ii) = th_ref(ii)

      End Do
!
!-----------------------------------------------------------------------
!! 0.2 Calculate temperature and pressure of lifting condensation level
!!     using approximations from Bolton (1980)
!-----------------------------------------------------------------------
!
!   vapour pressure e ~ qp/epsilon   q specific humidity
!   vapour pressure e ~ qp/(epsilon+q)   q mixing ratio

      If (l_mixing_ratio) then
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
! expression for mixing ratio
          vap_press = 0.01*Q_c(ii,k_plume(ii)) *                       &
     &                                 P_theta_lev_c(ii,k_plume(ii))   &
     &                      / (epsilon+Q_c(ii,k_plume(ii)) )
          If (vap_press  >   0.0) then
           T_LCL(ii) = a_bolton + b_bolton/                            &
     &                      (c_bolton*LOG(T(ii,k_plume(ii)))           &
     &                                    - LOG(vap_press) - d_bolton )
           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
     &                ( T_LCL(ii) / T(ii,k_plume(ii)) )**(1.0/kappa)
          Else
           P_LCL(ii) = pstar(i,j)
          End If
 
        End Do

      Else
! expression for specific humidity
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          vap_press = Q_c(ii,k_plume(ii)) *                            &
     &         P_theta_lev_c(ii,k_plume(ii)) / ( 100.0*epsilon )
          If (vap_press  >   0.0) then
           T_LCL(ii) = a_bolton + b_bolton/                            &
     &                      (c_bolton*LOG(T(ii,k_plume(ii)))           &
     &                                    - LOG(vap_press) - d_bolton )
           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
     &                ( T_LCL(ii) / T(ii,k_plume(ii)) )**(1.0/kappa)
          Else
            P_LCL(ii) = pstar(i,j)
          End If
 
        End Do

      End If ! test on l_mixing_ratio  
!
! Reset zh  (at this point in the code ntml is initialised as =1)
!
      do j=1,rows
        do i=1,row_length
          zh(i,j) = z_half(i,j,ntml(i,j)+1)
        End Do
      End Do
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        zh_c(ii) = zh(i,j)   
        zh2(ii)  = zh_c(ii)
      End Do
!
!-----------------------------------------------------------------------
! Find NLCL
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)    either
!     + + + + + + + +   lcl, Plcl, not a model level        lower part
!    ---------------   p      nlcl , p_theta(nlcl+1)         of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!     + + + + + + + +   lcl, Plcl, not a model level
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)     or
!                                                          upper part
!    ---------------   p      nlcl , p_theta_lev(nlcl+1)        of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
      DO  K = 2,wet_model_levels
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If ( P_LCL(ii)  <   P(i,j,K) ) then
! compressed copies
            nlcl_c(ii) = K-1
            z_lcl_c(ii)    = z_half_c(ii,nlcl_c(ii)+1)

! expand to full arrays
            nlcl(i,j) = K-1
            z_lcl(i,j)    = z_half_c(ii,nlcl_c(ii)+1)
            z_lcl_uv(i,j) = z_full_c(ii,nlcl_c(ii))
          End If     ! test on p_lcl
        End Do       ! ii loop
      End Do         ! k loop

!-----------------------------------------------------------------------
! 4.0 Parcel ascent - only perform parcel ascent If unstable
!-----------------------------------------------------------------------
! Note initial parcel conditions different from those used in the G-R
! mass flux convection scheme.
!
!-----------------------------------------------------------------------
! Calculate parcel water by linearising qsat about the Parcel's
! temperature extrapolated up to the next grid_level.
!----------------------------------------------------------------------

       DO  K = 1,wet_model_levels

! Require t_ref on all point for qsat call
         Do ii=1, nunstable
           t_ref(ii) = th_ref(ii)*exner_theta_levels_c(ii,k)
         End Do

! DEPENDS ON: qsat_mix
         Call qsat_mix(qsat_lev,t_ref,p_theta_lev_c(1,k)               &
     &                                    ,nunstable,l_mixing_ratio)

         Do ii=1, nunstable
           If(T_ref(ii) >  TM) then
             lrcp_const = lcrcp
             l_const    = lc
           Else
             lrcp_const = lsrcp
             l_const    = ls
           End If

! dqsat/dT - same whether q specific humidity or mixing ratio

           dq_sat_env = epsilon*l_const*qsat_lev(ii)/(R*T_ref(ii)**2.0)

           q_liq_parc = max( 0.0, ( qw_plume(ii) - QSat_lev(ii)        &
     &       -dq_sat_env*( sl_plume(ii)-grcp*z_full_c(ii,K)-T_ref(ii) )&
     &                             ) / (1.0+lrcp_const*dq_sat_env) )
           q_liq_env  = max( 0.0, ( qw(ii,K) - QSat_lev(ii)            &
     &        -dq_sat_env*( TL(ii,K)               - T_ref(ii) )       &
     &                             ) / (1.0+Lrcp_const*dq_sat_env) )
!
! add on the difference in the environment's ql as calculated by the
! UM cloud scheme (using some RH_CRIT value) and what it
! would be If RH_CRIT=1. This then imitates partial condensation
! in the parcel.
!
           q_liq_parc = q_liq_parc + qcl_c(ii,k)                       &
     &                            + qcf_c(ii,k)- q_liq_env
           T_parc(ii,k)=sl_plume(ii)-GRCP*z_full_c(ii,K)               &
     &                                  +lrcp_const*q_liq_parc

! May need to recalculate if T_parc is > Tm and T_ref < Tm

           If (T_ref(ii) <= TM.and.T_parc(ii,k) >  TM) then

! recalculate using corrected latent heats
             lrcp_const_parc = lcrcp
             q_liq_parc = max( 0.0, ( qw_plume(ii) - Qsat_lev(ii)      &
     &      -dq_sat_env*( sl_plume(ii)-grcp*z_full_c(ii,K)-T_ref(ii) ) &
     &                     ) / (1.0+lrcp_const_parc*dq_sat_env) )
             q_liq_parc = q_liq_parc + qcl_c(ii,k)                     &
     &                               + qcf_c(ii,k)- q_liq_env

! revised at parcel calculation

             T_parc(ii,k)=sl_plume(ii)-GRCP*z_full_c(ii,K)             &
     &                                  +lrcp_const_parc*q_liq_parc

           End If

           q_vap_parc=qw_plume(ii)-q_liq_parc
!
           t_dens_parc(ii,k)=T_parc(ii,k)*                             &
     &                    (1.0+c_virtual*q_vap_parc-q_liq_parc)


           t_dens_env(ii,k)=T(ii,K)*                                   &
     &                 (1.0+c_virtual*Q_c(ii,K)-qcl_c(ii,k)-qcf_c(ii,k))

           buoyancy(ii,k) = t_dens_parc(ii,k) - t_dens_env(ii,k)

           env_svl(ii,k) = t_dens_env(ii,k)  + GRCP*z_full_c(ii,K)
           par_svl(ii,k) = t_dens_parc(ii,k) + GRCP*z_full_c(ii,K)

           If (k >= 2) then

!         !-------------------------------------------------------------
!         ! Find vertical gradients in parcel and environment SVL
!         ! (using values from level below (i.e. K-1)).
!         !-------------------------------------------------------------

             dz = z_full_c(ii,K) - z_full_c(ii,K-1)

             dpar_bydz(ii,k) = (par_svl(ii,k) - par_svl(ii,k-1))/dz
             denv_bydz(ii,k) = (env_svl(ii,k) - env_svl(ii,k-1))/dz

           End If   ! test on k

! calculate t_ref for next level
          If (k >  1 .and. k <   wet_model_levels-1) then
             z_pr = (z_full_c(ii,k+1)-z_full_c(ii,k))                  &
     &                         /(z_full_c(ii,k)-z_full_c(ii,K-1))
             th_par = T_parc(ii,k)/exner_theta_levels_c(ii,k)
             th_ref(ii) = th_par*(1.+z_pr) - th_par_km1(ii)*z_pr

! Check sensible value otherwise set to previous reference value
! Problems can occur near top of model where calculation are nolonger 
! important.
             If (th_ref(ii) < 0.0) then
               th_ref(ii) = th_par_km1(ii)
             End If
             If (th_par > 0.0) then   
               th_par_km1(ii) = th_par
             End If
  
          End If

        End Do    ! ii loop
      End Do      ! level loop

!-----------------------------------------------------------------------
! tests on parcel ascent
!-----------------------------------------------------------------------
!   Now compare plume s_VL with each model layer s_VL in turn to
!     find the first time that plume has negative buoyancy.
!-----------------------------------------------------------------------

      Do  K = 2,wet_model_levels

!-----------------------------------------------------------------------
! Only perform tests if parcel ascent If unstable
!-----------------------------------------------------------------------
        Do ii=1,nunstable

!  Find level just below 2.5km (for use in Cu diagnosis)

          If ( z_full_c(ii,K)  >   2500.0                              &
     &         .and. kcucheck(ii)  ==  1 ) kcucheck(ii) = K-1

! freezing level

          If (t(ii,k) <  TM.and.t(ii,k-1) >= TM) then  
            If (freeze_lev(ii) == 1) then
              freeze_lev(ii) = k    
            End If
          End If

!-----------------------------------------------------------------------
! Set flag to true when level below is at least one level above the lcl
! and above the lcl transition zone
! Code implies ABOVE_LCL at NLCL+3 or greater.

          If (k-1 >  nlcl_c(ii)+1                                      &
     &                 .and. z_full_c(ii,k-1) >  1.1*z_lcl_c(ii)) then
            above_lcl(ii)=.true.
          Else
            above_lcl(ii)=.false.
          End If

!-----------------------------------------------------------------------
! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
!-----------------------------------------------------------------------
! Not reached LNB continue testing

          If ( .not.topprof(ii).and.k >  k_plume(ii) )then

            If (buoyancy(ii,k) >  max_buoy(ii)) then
              max_buoy(ii) = buoyancy(ii,k)
              k_max(ii)    = k  
            End If 

! Is parcel still buoyant ?

            If ( (buoyancy(ii,k)  <=  - thv_pert(ii))                  &
!                      or reached top of model
     &           .OR. (K  >   wet_model_levels-1)  ) then

              k_neutral(ii) = k-1
              topprof(ii) = .true.
              zh_c(ii) = z_half_c(ii,K)

! Buoyancy at last buoyant level

              Dt_dens_parc_T(ii) = buoyancy(ii,k-1)

              If ( delthvu_c(ii)  >   0.0) then
! compensate for any negative buoyancy of parcel in cloud layer
                delthvu_c(ii) = delthvu_c(ii) - dtv_min(ii) *          &
     &                                ( z_half_c(ii,K) - z_lcl_c(ii) )
              End If                                                     
            End If
          End If

!-----------------------------------------------------------------------
! Tests applied once found top of parcel ascent.
! Aim - to establish if the ascent has an inversion above the top
!       i.e. the ascent may indicate shallow /congestus convection.
! Sets indicator shmin = .true. if conditions met and stops testing.
!
! Conditions are ;
! either  denv/dz(k) < dpar/dz(k)
!   or    par_svl(k-1) -env_svl(k-1) <= 0.0
!
!-----------------------------------------------------------------------

          If ( topbl(ii) .and. ( denv_bydz(ii,k)  <   dpar_bydz(ii,k)  &
     &                .OR. buoyancy(ii,k-1)  <=  0.0 )                 &
     &                              .and. .not. shmin(ii) ) then
            shmin(ii) = .TRUE.
            Dt_dens_parc_TMIN(ii) = buoyancy(ii,k-1)
            kshmin(ii) = K-1
          End If

!-----------------------------------------------------------------------
! Tests applied to find parcel top
!
!-----------------------------------------------------------------------

          If ( .not.topbl(ii) .and. K  >   k_plume(ii) .and.           &
     & (  ( buoyancy(ii,k) <=  - thv_pert(ii)).OR.                     &
!
!                      plume non buoyant
!
     & (above_lcl(ii).and.(denv_bydz(ii,k) >  1.25*dpar_bydz(ii,k)))   &
!
!                      or environmental virtual temperature gradient
!                      signIficantly larger than parcel gradient
!                      above lIfting condensation level
!
     &           .OR. (K  >   wet_model_levels-1)                      &
!                      or reached top of model
     &         )                                                       &
     &         ) then
!
            topbl(ii) = .TRUE.
            zh2(ii) = z_half_c(ii,K)
            k_inv(ii) = K-1

            Dt_dens_parc_T2(ii) = buoyancy(ii,k-1)
            If ( delthvu2(ii)  >   0.0) then
! compensate for any negative buoyancy of parcel in cloud layer
              delthvu2(ii) = delthvu2(ii) - dtv_min2(ii) *             &
     &                              ( z_half_c(ii,k) - z_lcl_c(ii) )
            End If
          End If          ! test on .not.topbl

!-----------------------------------------------------------------------
! While doing parcel ascent
! (a) find minimum buoyancy
! (b) integrate CAPE over the ascent
!-----------------------------------------------------------------------

          If (K > nlcl_c(ii) .and. k < wet_model_levels ) then

            inc = g  * buoyancy(ii,k)                             &
     &          * (z_half_c(ii,K+1) - z_half_c(ii,K))/t_dens_env(ii,k)

            !---------------------------------------------------------- 
            ! If not reached an inversion or level of neutral buoyancy
            !---------------------------------------------------------- 

            If (.not. topbl(ii)) then

              dtv_min2(ii) = MIN( dtv_min2(ii),                        &
     &                       buoyancy(ii,k)/exner_theta_levels_c(ii,k) )

              delthvu2(ii) = delthvu2(ii) +                            &
     &              buoyancy(ii,k)*(z_half_c(ii,K+1) - z_half_c(ii,K)) & 
     &                        /exner_theta_levels_c(ii,k)
            End If

            !---------------------------------------------------------- 
            ! If not reached level of neutral buoyancy (top of ascent)
            !---------------------------------------------------------- 

            If (.not. topprof(ii)) then
                
              ! Note only calculating CIN and CAPE from ascents reaching
              ! level of neutral buoyancy. This may not always correspond 
              ! to the diagnosed top for the convection scheme.

              If (inc <  0.0) then
                 CIN_c(ii)  = CIN_c(ii) + inc
              End If
              CAPE_c(ii) = CAPE_c(ii) + inc

              dtv_min(ii) = MIN( dtv_min(ii),                         &
     &                     buoyancy(ii,k)/exner_theta_levels_c(ii,k)  )

              delthvu_c(ii) = delthvu_c(ii) +                         &
     &            buoyancy(ii,k)*(z_half_c(ii,K+1) - z_half_c(ii,K))  & 
     &                        /exner_theta_levels_c(ii,k)
 
            End If    ! test on topprof

          End If

!-----------------------------------------------------------------------
        End Do   ! ii loop

      End Do     ! level loop

!-----------------------------------------------------------------------
! Average vertical velocity over a layer  - required for shallow
!   convection test.
!-----------------------------------------------------------------------
! Layer from top in cloud w to value 1500km above cloud top?
!-----------------------------------------------------------------------
      Do ii=1,nunstable
        w_avg(ii) = 0.0
        mass(ii)  = 0.0
      End Do
      Do k=1,model_levels-1
        Do ii=1,nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If (k >= k_neutral(ii).and.                                     &
     &       z_full_c(ii,k) <= (z_half_c(ii,k_neutral(ii)+1)+1500.)) then

            mass(ii)  = mass(ii) + dmass_theta(ii,k)
            w_avg(ii) = w_avg(ii) + w_copy(i,j,k)*dmass_theta(ii,k)

          End If
        End Do
      End Do
      Do ii=1,nunstable
        if (mass(ii)  >  0.0 ) then
          w_avg(ii) = w_avg(ii)/mass(ii)
        endif
      End Do

      Do ii=1,nunstable
        w_avg2(ii) = 0.0
        mass(ii)   = 0.0
      End Do
      Do k=1,model_levels-1
        Do ii=1,nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If (k >= k_inv(ii) .and.                                     &
     &       z_full_c(ii,k) <= (z_half_c(ii,k_inv(ii)+1)+1500.)) then

            mass(ii)   = mass(ii) + dmass_theta(ii,k)
            w_avg2(ii) = w_avg2(ii) + w_copy(i,j,k)*dmass_theta(ii,k)
          End If
        End Do
      End Do
      Do ii=1,nunstable
        If (mass(ii)  >  0.0 ) then
          w_avg2(ii) = w_avg2(ii)/mass(ii)
        End If
      End Do

!-----------------------------------------------------------------------
! Default parcel top properties are assumed to be those when the
! ascent reaches the level of neutral buoyancy. These may not be those
! required in the case of shallow convection.
! Shallow convection requires the possible identifcation of an inversion
! at the top of the ascent. This may not be detected by the LNB test.
! The gradient tests are designed to detect the shallow top.
!-----------------------------------------------------------------------
! Modify top based on topbl test if ascent is likely to be shallow

      Do ii=1,nunstable

        If (shmin(ii) ) then    ! found an inversion    
! points where k_inv not the same as k_neutral and level below freezing
! may be shallow or congestus or deep

          If (k_inv(ii) == k_neutral(ii)) then
!  Both methods give same answer for top level leave shmin set
            ntml_c(ii) = k_neutral(ii)

! Inversion top lower than level of neutral buoyancy.
! Check also, either below freezing level or less than 2500m for shallow 
! convection.
          Else if ((k_inv(ii) <  freeze_lev(ii) .or.                   &
     &                     z_full_C(ii,k_inv(ii)+1)  <=  2500.0 )      &
     &             .and. k_inv(ii) <  k_neutral(ii) )then     


            If ( (z_full_c(ii,kshmin(ii)) - z_full_c(ii,k_inv(ii)))    &
     &         <=  1.25*(z_half_c(ii,k_inv(ii)+1) - z_lcl_c(ii)).and.  &
     &         (dt_dens_parc_tmin(ii)  <=  0.55*dt_dens_parc_t2(ii))   &
     &         .and.     (w_avg2(ii)  <   0.0)  ) then

! May be shallow or congestus
! set values to those found from inversion testing
               ntml_c(ii)  = k_inv(ii)
               delthvu_c(ii) = delthvu2(ii)
               zh_c(ii)    = zh2(ii)
               w_avg(ii)   = w_avg2(ii)
               dt_dens_parc_t(ii) = dt_dens_parc_t2(ii)

            Else   ! Assume not shallow or congestus
               ntml_c(ii) = k_neutral(ii)
               shmin(ii) = .false.  ! inversion top found not good 
                                   ! don't do shallow tests
            Endif
          Else   ! Assume deep  and therefore top LNB
              ntml_c(ii) = k_neutral(ii)
              shmin(ii) = .false.  ! inversion top found not good 
                                   ! don't do shallow tests
          Endif

        Else    !  No inversion found  i.e. shmin=false
          ntml_c(ii) = k_neutral(ii)
        Endif   ! shmin test

      End Do    ! ii loop

! Is there any need to recalculate w_avg  for any points ?
! Assuming not but may need to check this.
!-----------------------------------------------------------------------
!!     Save parcel ascent top: this will be used to allow mixing and
!!     entrainment into decoupled Sc of single layer thickness when it
!!     occurs above Cu.
!-----------------------------------------------------------------------
! Expand back up to full arrays 

      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        ntpar(i,j) = ntml_c(ii)
        zh(i,j)    = zh_c(ii)
        ntml(i,j)  = ntml_c(ii)
        nlcl(i,j)  = nlcl_c(ii)
        delthvu(i,j) = delthvu_c(ii)
        CAPE(i,j)    = CAPE_c(ii)
        CIN(i,j)     = CIN_c(ii)
      End Do

      Do j = 1,rows
        Do i = 1, row_length
          zhpar(i,j) = zh(i,j)
        End Do
      End Do
!-----------------------------------------------------------------------
!     Test height derived above against lifting condensation level
!-----------------------------------------------------------------------

! loop over only unstable 

!CDIR NODEP
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
!-----------------------------------------------------------------------
!     Check lifting condensation levels against height of parcel ascent,
!     If lifting condensation level lower than parcel ascent, and is
!     within bl_levels, then decide
!     on type of cloudy layer. If lifting condensation level at or below
!     low grid point, assume fog layer and turbulent mixing. For
!     gradient tests assume any If LCL and top of parcel ascent is less
!     than two levels then stratocumulus.
!-----------------------------------------------------------------------
! 4A code  test ntml-nlcl >= 2 (5A code requires 3 cloud levels)

        If ( ntml(i,j)-nlcl(i,j)  >=  2                                 &
     &                        .and. nlcl(i,j)  >   k_plume(ii)          &
     &                              .and. nlcl(i,j)  <   MBL-1 ) then
!-----------------------------------------------------------------------
!     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
!     For stratocumulus top of mixed layer = zh
!     For cumulus top of mixed layer = ZLCL
!     New test zhpar >= 3000m replaces (NTML >= MBL) i.e. parcel top 
!     is above boundary layer indicates convection.
!     Diagnosis is done by comparing gradients
!-----------------------------------------------------------------------

          If (zhpar(I,j) >= 3000.0) then
            cumulus(i,j) = .TRUE.
          Else


! Current test is against a height of ~<2.5km
! This could be replaced by a scale height if a suitable method
! for determining a sensible height was possible from profile/cumulus
! depth information available in this routine
! Assume moisture gradient tests stay the same whether specific or 
! mixing ratio though imply slightly different moisture gradients.

            If (ntml(i,j)  >   kcucheck(ii)                            &
     &             .and. nlcl(i,j)  <=  kcucheck(ii)-2) then

              grad_cld =  ABS( QW(ii,kcucheck(ii)) -                    &
     &                                  QW(ii,nlcl_c(ii)) ) /           &
     &          ( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)) )
            Else
              grad_cld =  ABS( QW(ii,ntml_c(ii)) -                      &
     &                                  QW(ii,nlcl_c(ii)) ) /           &
     &              ( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)) )
            End If

            grad_sub =  ABS( QW(ii,nlcl_c(ii)) -                        &
     &                                QW(ii,k_plume(ii)) ) /            &
     &           ( z_full_c(ii,nlcl_c(ii)) - z_full_c(ii,k_plume(ii)) )

            If (grad_cld  >   1.10*grad_sub) then
!-----------------------------------------------------------------------
!     Not well mixed, however it is possible that the depth of a well
!     mixed boundary layer has increased but not yet been mixed yet so
!     test gradient from next level down.
!     Note typical cumulus profiles are expected to have a fairly
!     uniform q profile from the surface to the cloud base and then a
!     decreasing profile of q above this in the cloud. Typical the
!     decreasing gradient from the cloud base to 2.5km will be the
!     order of > 1.10 the below cloud value.
!-----------------------------------------------------------------------

! test against a height ~ 2.5km

              If (ntml_c(ii)  <=  kcucheck(ii)) then
              grad_cld =  ABS( QW(ii,ntml_c(ii)-1) -                    &
     &                                  QW(ii,nlcl_c(ii)) ) /           &
     &           ( z_full_c(ii,ntml_c(ii)-1) - z_full_C(ii,nlcl_c(ii)) )
              End If

              If ( grad_cld  >   1.10*grad_sub) then
!-----------------------------------------------------------------------
!      Diagnose a cumulus layer
!-----------------------------------------------------------------------
                cumulus(i,j) = .TRUE.
              End If
            Else

! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
! and not yet been mixed (so could have been erroneously identIfied as
! well-mixed)

! First check using level below (recalculate grad_sub)

              If (nlcl_c(ii) - k_plume(ii)  >=  2) then
                 grad_sub =  ABS( QW(ii,nlcl(i,j)-1) -                 &
     &                                QW(ii,k_plume(ii)) ) /           &
     &          ( z_full_c(ii,nlcl_c(ii)-1) - z_full_c(ii,k_plume(ii)) )

                 If ( grad_cld  >   1.10*grad_sub) then
                   cumulus(i,j) =.TRUE.
                 End If

              End If

! If still diagnosing well-mixed, check using level above
! (recalculate grad_cld)

              If (.not. cumulus(i,j) ) then

               If (ntml_c(ii)  >   kcucheck(ii)                        &
     &             .and. nlcl_c(ii)  <=  kcucheck(ii)-2) then

                grad_cld =  ABS( QW(ii,kcucheck(ii)) -                 &
     &                                  QW(ii,nlcl_c(ii)+1) ) /        &
     &         ( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               Else
                grad_cld =  ABS( QW(ii,ntml_c(ii)) -                   &
     &                                  QW(ii,nlcl_c(ii)+1) ) /        &
     &            ( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               End If

               If ( grad_cld  >   1.10*grad_sub) then
                 cumulus(i,j) =.TRUE.
               End If

              End If
            End If
          End If
        End If
      End Do       ! ii loop 


!-----------------------------------------------------------------------
!      Check that a cumulus layer has not been erroneously diagnosed in
!      a deep cloudy region
!      As the above checks are done on the total water rather than q it
!      is possible the conditions can be met in areas where the level of
!      prognostic qcl or qcf is high. The type of mistake is only
!      thought to occur over land.
!-----------------------------------------------------------------------
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        K=nlcl(i,j)

        If ( land_mask(i,j) .and. cumulus(i,j) .and.                    &
     &                                 ntpar(i,j)  <   MBL ) then
          Do while( K  <=  ntpar(i,j) .and. cloud_fraction(i,j,K)       &
     &                                            >=  SC_CFTOL )
            K = K + 1
          End Do
          If (K  ==  ntpar(i,j)+1) cumulus(i,j) = .false.
        End If

      End Do       ! ii loop 

!-----------------------------------------------------------------------
! Original shallow diagnosis no congestus diagnosis
!-----------------------------------------------------------------------
!CDIR NODEP
        Do ii=1,nunstable
          i = index_i(ii) 
          j = index_j(ii) 
          If ( cumulus(i,j) ) then

!-----------------------------------------------------------------------
!       If cumulus has been diagnosed, determine whether it is shallow
!       or deep convection
!-----------------------------------------------------------------------
! Conditions for shallow convection
!    wadv < 0.0            (descending air)
!   top of parcel ascent < 2500. or T (top of parcel) > TM
!   height of min buoyancy (above Bl) - height of parcel top T level
!      <1.25(height parcel top - z lifting condensation level)
!   t_dens_parc -t_dens_env at kshmin <0.55t_dens_parc -t_dens_env
!   at ntpar
!
!   The last 2 conditions are looking for a strong inversion at the top
!   of the shallow cumulus.
!-----------------------------------------------------------------------
           If ( shmin(ii) ) then

!changeed by hongyan            If ( w_avg(ii)  <   0.0 .and.                              &
!            If (                               &
! sza: discussion in couple meeting recommands not use above correction now
            If ( w_avg(ii)  <   0.0 .and.                              &
     &        (z_full_c(ii,ntpar(i,j))  <=  2500.0 .OR.                &
     &            T(ii,ntpar(i,j))  >=  TM)                            &
     &       .and. (z_full_c(ii,kshmin(ii)) - z_full_c(ii,ntpar(i,j))) &
     &          <=  1.25*(zhpar(i,j) - z_lcl_c(ii)) .and.              &
     &         Dt_dens_parc_Tmin(ii)  <=  0.55*Dt_dens_parc_T(ii) )    &
     &      then

              L_shallow(i,j) = .TRUE.
! may be problem with ntpar diagnosis for deep if wadv test sets
! L_shallow  false

            End If

           End If       ! test on shmin

!-----------------------------------------------------------------------
!      Set mixed layer depth to z_lcl
!-----------------------------------------------------------------------
           If (P_LCL(ii)  <   (P_theta_lev_c(ii,nlcl(i,j)+1))) then
!-----------------------------------------------------------------------
!      If LCL is diagnosed in the upper half of the layer set z_lcl to
!      the height of the upper layer interface
!      (in code above LCL is always set to the lower interface).
!-----------------------------------------------------------------------
             nlcl(i,j) = nlcl(i,j)+1
             z_lcl(i,j) = z_half_c(ii,nlcl(i,j)+1)
             z_lcl_uv(i,j)= z_full_c(ii,nlcl(i,j))
           End If
           zh(i,j) = z_lcl(i,j)
           ntml(i,j) = nlcl(i,j)

!      If cumulus has been diagnosed but delthvu is negative, reset
!      cumulus and L_shallow to FALSE but leave zh and ntml at LCL

           If (delthvu(i,j)  <=  0.0) then

             cumulus(i,j)   = .false.
             L_shallow(i,j) = .false.

           End If

        Else      ! not cumulus

!-----------------------------------------------------------------------
!      If not cumulus, reset parameters to within bl_levels
!-----------------------------------------------------------------------
          If (ntml(i,j)  >   MBL) then
            ntml(i,j)  = MBL
            ntpar(i,j) = MBL
            zh(i,j)    = z_half(i,j,MBL+1)
            zhpar(i,j) = zh(i,j)
          End If
          If (nlcl(i,j)  >   MBL) then
            nlcl(i,j)    = MBL
            z_lcl(i,j)   = zh(i,j)
            z_lcl_uv(i,j)= z_full_c(ii,MBL-1)
          End If

        End If        ! test on cumulus

      End Do          ! ii loop  

!=======================================================================
! Option 2. dilute parcel - deep entrainment rates ? 
!=======================================================================
! Calculation of LCL unchanged
! Parcel ascent now dilute
! 
!=======================================================================

      Else if (icvdiag >= 2 .and. icvdiag < 4) then

        l_keep_water  = .false. ! water loading not kept, reduced if >1g/kg
        l_wtest  = .true.       ! w test
! Monsoon test
!        l_wtest  = .false.       ! w test
     
!-----------------------------------------------------------------------
! Set entrainment fractions
!-----------------------------------------------------------------------
      If (icvdiag == 2) then

!  0.55/z entrainment rates eg Jakob & Siebesma   e*dz

        Do k= 1,wet_model_levels-1
          Do ii=1, nunstable
            i = index_i(ii)   
            j = index_j(ii)   
            entrain_fraction(ii,k) = 0.55*                             &
     &                 (r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))     &
     &                 /z_full_c(ii,k)
          End do
        End do 

      Else If (icvdiag ==3) then

! 1/z entrainment rates - similar to many CRM results for deep and 
!                         shallow convection.

        Do k= 1,wet_model_levels-1
          Do ii=1, nunstable
            i = index_i(ii)   
            j = index_j(ii)   
            entrain_fraction(ii,k) = 1.0 *                             &
     &                 (r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))     &
     &                 /z_full_c(ii,k)

          End do
        End do 

      Else
! Deep entrainment rates NOT being used at present.
!      Do k= 1,wet_model_levels-1
!        Do ii=1, nunstable
!          i = index_i(ii)   
!          j = index_j(ii)   
!            entrain_fraction(ii,k) = -1.*entcoef*ae2*p_theta_lev_c(ii,k)* &
!     &                  (p(i,j,k+1)-p(i,j,k))/( pstar(i,j)*pstar(i,j))
!        End do
!      End do 
        write(6,*) ' not defined'  
      Endif

      k=wet_model_levels
        Do ii=1, nunstable
            entrain_fraction(ii,k) = 0.0
        End do

!-----------------------------------------------------------------------
! Work out initial parcel properties and LCL
!-----------------------------------------------------------------------
      Do ii=1, nunstable

        k_plume(ii) = 1
!-----------------------------------------------------------------------
! Only perform parcel ascent If unstable
! Start plume ascent from grid-level above top of surface layer, taken
! to be at a height, z_surf, given by 0.1*zh
!-----------------------------------------------------------------------
        z_surf = 0.1 * zh_c(ii)

        Do while( z_full_c(ii,k_plume(ii))  <   z_surf .and.           &
!                   ! not reached z_surf
     &              SVL(ii,k_plume(ii)+1)  <   SVL(ii,k_plume(ii)) )
!                   ! not reached inversion

          k_plume(ii) = k_plume(ii) + 1

        End Do
      End Do

      Do ii=1, nunstable
         sl_plume(ii) = TL(ii,k_plume(ii))                             &
     &                        + GRCP * z_full_c(ii,k_plume(ii))
         thv_pert(ii) = max( a_plume,                                  &
     &                 min( max_t_grad*zh_c(ii), b_plume*tv1_sd(ii) ) )
         qw_plume(ii) = QW(ii,k_plume(ii))

! Added for more acturate parcel cal later

         th_ref(ii) = tl(ii,k_plume(ii))                               &
     &                        /exner_theta_levels_c(ii,k_plume(ii))
         th_par_km1(ii) = th_ref(ii)

! 2nd set for undilute calculation
         th_ref_dil(ii)    = th_ref(ii)
         th_par_km_dil(ii) = th_ref_dil(ii)

      End Do

!-----------------------------------------------------------------------
! Calculate temperature and pressure of lifting condensation level
!     using approximations from Bolton (1980)
!-----------------------------------------------------------------------
!
!   vapour pressure e ~ qp/epsilon   q specific humidity
!   vapour pressure e ~ qp/(epsilon+q)   q mixing ratio

      If (l_mixing_ratio) then      ! expression for mixing ratio

        Do ii=1, nunstable
          vap_press = 0.01*q_c(ii,k_plume(ii)) *                       &
     &                                 P_theta_lev_c(ii,k_plume(ii))   &
     &                      / (epsilon+q_c(ii,k_plume(ii)) )
          If (vap_press  >   0.0) then
           T_LCL(ii) = a_bolton + b_bolton/                            &
     &                      (c_bolton*LOG(T(ii,k_plume(ii)))           &
     &                                    - LOG(vap_press) - d_bolton )

           P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                 &
     &                ( T_LCL(ii) / T(ii,k_plume(ii)) )**(1.0/kappa)
          Else
            i = index_i(ii)   
            j = index_j(ii)   
            P_LCL(ii) = pstar(i,j)
          End If
        End Do

      Else                           ! expression for specific humidity

        Do ii=1, nunstable
          vap_press = q_c(ii,k_plume(ii)) *                            &
     &                P_theta_lev_c(ii,k_plume(ii)) / ( 100.0*epsilon )

          If (vap_press  >   0.0) then
            T_LCL(ii) = a_bolton + b_bolton/                           &
     &                      (c_bolton*LOG(T(ii,k_plume(ii)))           &
     &                                    - LOG(vap_press) - d_bolton )

            P_LCL(ii) = P_theta_lev_c(ii,k_plume(ii)) *                &
     &                ( T_LCL(ii) / T(ii,k_plume(ii)) )**(1.0/kappa)
          Else
            i = index_i(ii)   
            j = index_j(ii)   
            P_LCL(ii) = pstar(i,j)
          End If
        End Do

      End If ! test on l_mixing_ratio  

!
! Reset zh  (at this point in the code ntml is initialised as =1)
!
      Do j=1,rows
        Do i=1,row_length
          zh(i,j) = z_half(i,j,ntml(i,j)+1)
        End Do
      End Do
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        zh_c(ii) = zh(i,j)   
      End Do


!
!-----------------------------------------------------------------------
! Find NLCL
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)    either
!     + + + + + + + +   lcl, Plcl, not a model level        lower part
!    ---------------   p      nlcl , p_theta(nlcl+1)         of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
!
!    ---------------   p      nlcl+1  , p_theta(nlcl+2)
!     + + + + + + + +   lcl, Plcl, not a model level
!
!    - - - - - - - -   uv     nlcl+1,  z_lcl , p(nlcl+1)     or
!                                                          upper part
!    ---------------   p      nlcl , p_theta_lev(nlcl+1)        of layer
!
!    - - - - - - - -   uv     nlcl   p(nlcl)
!
!-----------------------------------------------------------------------
      DO  K = 2,wet_model_levels
        Do ii=1, nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If ( P_LCL(ii)  <   P(i,j,K) ) then
! compressed copies
            nlcl_c(ii) = K-1
            z_lcl_c(ii)    = z_half_c(ii,nlcl_c(ii)+1)

! expand to full arrays
            nlcl(i,j) = K-1
            z_lcl(i,j)    = z_half_c(ii,nlcl_c(ii)+1)
            z_lcl_uv(i,j) = z_full_c(ii,nlcl_c(ii))
          End If
        End Do       !  ii loop 
      End Do         !  k  loop

!-----------------------------------------------------------------------
! Parcel ascent - only perform parcel ascent If unstable
!-----------------------------------------------------------------------
! Parcel starts from level k_plume and is lifted up
! Dilute parcel ascent - mix in environmental air above lifting 
! condensation level
! 
!  Initial testing want tparc_undilute and t_parc (dilute)
!-----------------------------------------------------------------------
! Calculate parcel water by linearising qsat about the Parcel's
! temperature extrapolated up to the next grid_level.
!----------------------------------------------------------------------

      DO  K = 1,wet_model_levels

! Require t_ref on all point for qsat call
        Do ii=1, nunstable
          t_ref(ii)     = th_ref(ii)*exner_theta_levels_c(ii,k)
          t_ref_dil(ii) = th_ref_dil(ii)*exner_theta_levels_c(ii,k)
        End Do

! DEPENDS ON: qsat_mix
        call qsat_mix(qsat_lev,t_ref,p_theta_lev_c(1,k)                &
     &                                       ,nunstable,l_mixing_ratio)

! Second set
! DEPENDS ON: qsat_mix
        call qsat_mix(qsat_lev_dil,t_ref_dil,p_theta_lev_c(1,k)        &
     &                                       ,nunstable,l_mixing_ratio)

        Do ii=1, nunstable

! Undilute parcel calculation required as a reference

          If(T_ref(ii) >  TM) then
            lrcp_const = lcrcp
            l_const    = lc
          Else
            lrcp_const = lsrcp
            l_const    = ls
          End If

          dq_sat_par = epsilon*l_const*qsat_lev(ii)/(R*T_ref(ii)**2.0)

          q_liq_parc = MAX( 0.0, ( qw_plume(ii) - QSat_lev(ii)         &
     &       -dq_sat_par*( sl_plume(ii)-grcp*z_full_c(ii,K)-T_ref(ii) )&
     &                             ) / (1.0+lrcp_const*dq_sat_par) )

          q_liq_env  = MAX( 0.0, ( qw(ii,K) - QSat_lev(ii)             &
     &        -dq_sat_par*( TL(ii,K)               - T_ref(ii) )       &
     &                             ) / (1.0+Lrcp_const*dq_sat_par) )
!
! add on the difference in the environment's ql as calculated by the
! UM cloud scheme (using some RH_CRIT value) and what it
! would be If RH_CRIT=1. This then imitates partial condensation
! in the parcel.
!
          ql_parc(ii,k) = q_liq_parc + qcl_c(ii,k)                     &
     &                                   + qcf_c(ii,k)- q_liq_env
          T_PARC(ii,k)=sl_plume(ii)-GRCP*z_full_c(ii,K)                &
     &                                  +lrcp_const*q_liq_parc

! May need to recalculate if T_parc is > Tm and T_ref < Tm

          If (T_ref(ii) <= TM.and.T_parc(ii,k) >  TM) then

! recalculate using corrected latent heats
            lrcp_const_parc = lcrcp

            q_liq_parc = MAX( 0.0, ( qw_plume(ii) - Qsat_lev(ii)       &
     &      -dq_sat_env*( sl_plume(ii)-grcp*z_full_c(ii,K)-T_ref(ii) ) &
     &                     ) / (1.0+lrcp_const_parc*dq_sat_env) )

            ql_parc(ii,k) = q_liq_parc + qcl_c(ii,k)                   &
     &                                   + qcf_c(ii,k)- q_liq_env
! revised at parcel calculation

            T_PARC(ii,k)=sl_plume(ii)-GRCP*z_full_c(ii,K)              &
     &                                  +lrcp_const_parc*ql_parc(ii,k)

          End If

          q_vap_parc=qw_plume(ii)-ql_parc(ii,k)
!          
          t_dens_parc(ii,k)=T_PARC(ii,k)*                              &
     &                    (1.0+c_virtual*q_vap_parc-ql_parc(ii,k))


! calculate t_ref for next level
          If (k >  1 .and. k <   wet_model_levels-1) then
            z_pr = (z_full_c(ii,k+1)-z_full_c(ii,k))                   &
     &                         /(z_full_c(ii,k)-z_full_c(ii,K-1))
            th_par = t_parc(ii,k)/exner_theta_levels_c(ii,k)
            th_ref(ii) = th_par*(1.+z_pr) - th_par_km1(ii)*z_pr

! Check sensible value otherwise set to previous reference value
! Problems can occur near top of model where calculation are nolonger 
! important.
            If (th_ref(ii) < 0.0) then
              th_ref(ii) = th_par_km1(ii)
            End If
            If (th_par > 0.0) then   
              th_par_km1(ii) = th_par
            End If
          End If

! dilute parcel same as undilute parcel

          If (k <= nlcl_c(ii)) then

            sl_parc(ii,k) = sl_plume(ii)
            qw_parc(ii,k) = qw_plume(ii)
            t_parc_dil(ii,k)  = t_parc(ii,k)
            ql_parc_dil(ii,k) = ql_parc(ii,k)         
            t_dens_parc_dil(ii,k)= t_dens_parc(ii,k)
            th_ref_dil(ii) = th_ref(ii)
            th_par_km_dil(ii) = th_par_km1(ii)

          Else

! Dilute  parcel ascents now required

            If(t_ref_dil(ii) >  TM) then
              lrcp_const = lcrcp
              l_const    = lc
            Else
              lrcp_const = lsrcp
              l_const    = ls
            End If

!-----------------------------------------------------------------------
! Dilute parcel
!-----------------------------------------------------------------------
! Mix in entrain_fraction from environmental air from level below and 
! raise this to current level.
! Assume mix in fraction of mass from environment.
! Estimate parcel properties after mixing air from environment with 
! parcel. Temperature given approximately by average

            temp_parc = (t_parc_dil(ii,k-1)                            &
     &                          + entrain_fraction(ii,k)*t(ii,k-1))    &
     &                   /(1.+entrain_fraction(ii,k))

     
            qw_parc(ii,k) = (qw_parc(ii,k-1) +                 &
     &                         entrain_fraction(ii,k)*qw(ii,k-1))      & 
     &                    /(1.+entrain_fraction(ii,k)) 

            qcl_parc = (ql_parc_dil(ii,k-1)   +                        &    
     &                         entrain_fraction(ii,k)*qcl_c(ii,k-1))   & 
     &                    /(1.+entrain_fraction(ii,k)) 

            qcf_parc = (0.0     +                                      &    
     &                         entrain_fraction(ii,k)*qcf_c(ii,k-1))   & 
     &                    /(1.+entrain_fraction(ii,k)) 

! All condensed water either ice or liquid based on t_ref ?
            sl_parc(ii,k) = temp_parc - lrcp_const*(qcl_parc+qcf_parc) &
     &                           +grcp*z_full_c(ii,k-1) 

            dq_sat_par = epsilon*l_const*qsat_lev(ii)/(R*t_ref_dil(ii)**2.0)

            q_liq_parc = MAX( 0.0, ( qw_parc(ii,k) - qsat_lev_dil(ii)    &
     &     -dq_sat_par*( sl_parc(ii,k)-grcp*z_full_c(ii,K)-t_ref_dil(ii))&
     &                             ) / (1.0+lrcp_const*dq_sat_par) )

            q_liq_env  = MAX( 0.0, ( qw(ii,K) - qsat_lev_dil(ii)       &
     &        -dq_sat_par*( TL(ii,K)               - t_ref_dil(ii) )   &
     &                             ) / (1.0+Lrcp_const*dq_sat_par) )
!
! add on the dIfference in the environment's ql as calculated by the
! UM cloud scheme (using some RH_CRIT value) and what it
! would be If RH_CRIT=1. This then imitates partial condensation
! in the parcel.
!
            ql_parc_dil(ii,k) = q_liq_parc + qcl_c(ii,k)               &
     &                                     + qcf_c(ii,k)- q_liq_env

            t_parc_dil(ii,k)=sl_parc(ii,k)-GRCP*z_full_c(ii,K)         &
     &                                  +lrcp_const*ql_parc_dil(ii,k)

! May need to recalculate if T_parc is > Tm and T_ref < Tm

            If (t_ref_dil(ii) <= TM.and.t_parc_dil(ii,k) >  TM) then

! recalculate using corrected latent heats
               lrcp_const_parc = lcrcp

               q_liq_parc = MAX( 0.0, ( qw_parc(ii,k) - qsat_lev_dil(ii)     &
     &         -dq_sat_env*(sl_parc(ii,k)-grcp*z_full_c(ii,K)-t_ref_dil(ii)) &
     &                     ) / (1.0+lrcp_const_parc*dq_sat_env) )

               ql_parc_dil(ii,k) = q_liq_parc + qcl_c(ii,k)                  &
     &                               + qcf_c(ii,k)- q_liq_env

! revised at parcel calculation

              t_parc_dil(ii,k)=sl_parc(ii,k)-GRCP*z_full_c(ii,K)       &
     &                               +lrcp_const_parc*ql_parc_dil(ii,k)

            End if   ! test on t_ref

            q_vap_parc=qw_parc(ii,k)-ql_parc_dil(ii,k)

            If (.not.l_keep_water) then
!  water removed from parcel after condesation
              If (ql_parc_dil(ii,k).gt.0.001) then
                ql_parc_dil(ii,k) = 0.001
                qw_parc(ii,k) = q_vap_parc + ql_parc_dil(ii,k)
              End If
            End If
            t_dens_parc_dil(ii,k)=t_parc_dil(ii,k)*                    &
     &                    (1.0+c_virtual*q_vap_parc-ql_parc_dil(ii,k))

! calculate dilute t_ref for next level
            If (k >  1 .and. k <   wet_model_levels-1) then
              z_pr = (z_full_c(ii,k+1)-z_full_c(ii,k))                 &
     &                         /(z_full_c(ii,k)-z_full_c(ii,K-1))

              th_par = t_parc_dil(ii,k)/exner_theta_levels_c(ii,k)
              th_ref_dil(ii) = th_par*(1.+z_pr) - th_par_km_dil(ii)*z_pr
! Check new reference sensible
              If (th_ref_dil(ii) < 0.0) then
                th_ref_dil(ii) = th_par_km_dil(ii)
              End If
              If (th_par > 0.0) then   
                th_par_km_dil(ii) = th_par
              End If
            End If      ! k level test

          End If   ! test on LCL
 
          t_dens_env(ii,k)=T(ii,K)*                                    &
     &                  (1.0+c_virtual*Q_c(ii,K)-qcl_c(ii,k)-qcf_c(ii,k))

          buoyancy(ii,k)     = t_dens_parc(ii,k)     - t_dens_env(ii,k)
          buoyancy_dil(ii,k) = t_dens_parc_dil(ii,k) - t_dens_env(ii,k)

          env_svl(ii,k) = t_dens_env(ii,k)      + GRCP*z_full_c(ii,K)
          par_svl(ii,k) = t_dens_parc_dil(ii,k) + GRCP*z_full_c(ii,K)

          If (k >= 2) then

!         !-------------------------------------------------------------
!         ! Find vertical gradients in parcel and environment SVL
!         ! (using values from level below (i.e. K-1)).
!         !-------------------------------------------------------------

            dz = z_full_c(ii,K) - z_full_c(ii,K-1)

            dpar_bydz(ii,k) = (par_svl(ii,k) - par_svl(ii,k-1))/dz
            denv_bydz(ii,k) = (env_svl(ii,k) - env_svl(ii,k-1))/dz

          End If   ! test on k

        End Do    ! ii loop
      End Do      ! level loop

!-----------------------------------------------------------------------
! tests on parcel buoyancy
!-----------------------------------------------------------------------
!   Now compare plume s_VL with each model layer s_VL in turn to
!     find the first time that plume has negative buoyancy.
!-----------------------------------------------------------------------

      DO  K = 2,wet_model_levels

        Do ii=1,nunstable
!-----------------------------------------------------------------------
! Only perform tests if parcel ascent If unstable
!-----------------------------------------------------------------------

!  Find level just below 2.5km (for use in Cu diagnosis)

          If ( z_full_c(ii,K)  >   2500.0                             &
     &         .and. kcucheck(ii)  ==  1 ) kcucheck(ii) = K-1

! freezing level

          If (t(ii,k) <  TM.and.t(ii,k-1) >= TM) then  
            If (freeze_lev(ii) == 1) then
              freeze_lev(ii) = k    
            End If
          End If

!-----------------------------------------------------------------------
! No flag for above_lcl required. Reduce thv_pert by a factor dependent
! on height relative to LCL.

          If (k-1 >  nlcl_c(ii)+1                                        &
     &                   .and. z_full_c(ii,k-1) >  1.1*z_lcl_c(ii)) then

! decrease thv_pert by exp(-(z-zlcl)/2000.)

            factor = exp( (z_lcl_c(ii)-z_full_c(ii,k))*1.e-3)

! set to zero if z-zlcl >1000?

            If ((z_full_c(ii,k)-z_lcl_c(ii)).gt. 1000.) then
              factor =0.0
            End If 

          Else
            factor= 1.0  
          End If

!-----------------------------------------------------------------------
! Level of neutral buoyancy (LNB) & maximum buoyancy level below this
!-----------------------------------------------------------------------
! Not reached LNB continue testing

          If ( .not.topprof(ii).and.k >  k_plume(ii) )then
            If (buoyancy_dil(ii,k) >  max_buoy(ii)) then
              max_buoy(ii) = buoyancy_dil(ii,k)
              k_max(ii)    = k  
            End If 

! Is parcel still buoyant ?

            If ( (buoyancy_dil(ii,k)  <=  - thv_pert(ii)*factor)       &
!                      or reached top of model
     &           .OR. (K  >   wet_model_levels-1)  ) then

              k_neutral(ii) = k-1
              topprof(ii) = .true.
              zh_c(ii) = z_half_c(ii,K)

! Buoyancy at last buoyant level

              Dt_dens_parc_T(ii) = buoyancy_dil(ii,k-1)

              If ( delthvu_c(ii)  >   0.0) then
! compensate for any negative buoyancy of parcel in cloud layer
                delthvu_c(ii) = delthvu_c(ii) - dtv_min(ii) *          &
     &                                ( z_half_c(ii,K) - z_lcl_C(ii) )
              End If                                                     
            End If
          End If

!-----------------------------------------------------------------------
! While doing parcel ascent
! (a) find minimum buoyancy
! (b) integrate CAPE over the ascent
!-----------------------------------------------------------------------

          If (K > nlcl_c(ii) .and. k < wet_model_levels ) then

! Not reached top of ascent
            If (.not. topprof(ii)) then
              dtv_min(ii) = MIN( dtv_min(ii),                          &
     &                       buoyancy(ii,k)/exner_theta_levels_c(ii,k) )

! undilute value
              delthvu_c(ii) = delthvu_c(ii) + buoyancy(ii,k)*          &
     &                     ( z_half_c(ii,K+1) - z_half_c(ii,K) )       &
     &                        /exner_theta_levels_c(ii,k)

! calculation of CIN and CAPE from profiles - undilute values 
              inc =g * buoyancy(ii,k)                                  &
     &          * (z_half_c(ii,K+1) - z_half_c(ii,K))/t_dens_env(ii,k)

              If (inc <  0.0) then
                 CIN_c(ii)  = CIN_c(ii) + inc
              End if

              CAPE_c(ii) = CAPE_c(ii) + inc

            End If    ! test on topprof

          End If      ! test on level k

!-----------------------------------------------------------------------

        End Do   ! ii loop
      End Do

!-----------------------------------------------------------------------
! Default parcel top properties are assumed to be those when the
! ascent reaches the level of neutral buoyancy. These may not be those
! required in the case of shallow convection.
! Shallow convection requires the possible identifcation of an inversion
! at the top of the ascent. This may not be detected by the LNB test.
! The gradient tests are designed to detect the shallow top.
!-----------------------------------------------------------------------

      Do ii=1,nunstable
        ntml_c(ii) = k_neutral(ii)
      End Do
!-----------------------------------------------------------------------
!!     Save parcel ascent top: this will be used to allow mixing and
!!     entrainment into decoupled Sc of single layer thickness when it
!!     occurs above Cu.
!-----------------------------------------------------------------------
! Expand back up to full arrays 

      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        ntpar(i,j) = ntml_c(ii)
        zh(i,j)    = zh_c(ii)
        ntml(i,j)  = ntml_c(ii)
        nlcl(i,j)  = nlcl_c(ii)
        delthvu(i,j) = delthvu_c(ii)
        CAPE(i,j)    = CAPE_c(ii)
        CIN(i,j)     = CIN_c(ii)
      End Do
! set zhpar(i,j)
      Do j = 1,rows
        Do i = 1, row_length
          zhpar(i,j) = zh(i,j)
        End Do
      End Do

!-----------------------------------------------------------------------
! Average vertical velocity over a layer  - required for shallow
!   convection test.
!-----------------------------------------------------------------------
! Layer from top in cloud w to value 1500km above cloud top?
!-----------------------------------------------------------------------
      If (l_wtest) then
      Do ii=1,nunstable
        w_avg(ii) = 0.0
        mass(ii)  = 0.0
      End Do
      Do k=1,model_levels-1
        Do ii=1,nunstable
          i = index_i(ii)   
          j = index_j(ii)   
          If (k >= ntml_c(ii).and.                                     &
     &       z_full_c(ii,k) <= (z_half_c(ii,ntml_c(ii)+1)+1500.)) then

            mass(ii)  = mass(ii) + dmass_theta(ii,k)
            w_avg(ii) = w_avg(ii) + w_copy(i,j,k)*dmass_theta(ii,k)
          End If
        End Do
      End Do
      Do ii=1,nunstable
        If (mass(ii)  >  0.0 ) then
          w_avg(ii) = w_avg(ii)/mass(ii)
        End If
      End Do

      End If     ! test on w
!-----------------------------------------------------------------------
!     Test height derived above against lifting condensation level
!-----------------------------------------------------------------------
!CDIR NODEP
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
!-----------------------------------------------------------------------
!     Check lifting condensation levels against height of parcel ascent,
!     If lifting condensation level lower than parcel ascent, and is
!     within bl_levels, then decide
!     on type of cloudy layer. If lifting condensation level at or below
!     low grid point, assume fog layer and turbulent mixing. For
!     gradient tests assume any If LCL and top of parcel ascent is less
!     than two levels then stratocumulus.
!-----------------------------------------------------------------------
! 4A scheme  ntpar-ntml >=2  (must be 2 cloud levels, 5A code requires 3)

        If ( ntml(i,j)-nlcl(i,j)  >=  2                                 &
     &                        .and. nlcl(i,j)  >   k_plume(ii)          &
     &                              .and. nlcl(i,j)  <   MBL-1 ) then
!-----------------------------------------------------------------------
!     Cloudy boundary layer, diagnose whether stratocumulus or cumulus.
!     For stratocumulus top of mixed layer = zh
!     For cumulus top of mixed layer = ZLCL
!     New test zhpar >= 3000m replaces (NTML >= MBL) i.e. parcel top 
!     is above boundary layer indicates convection.
!     Diagnosis is done by comparing gradients
!-----------------------------------------------------------------------

          If (zhpar(I,j) >= 3000.0) then
            cumulus(i,j) = .TRUE.
          Else

! Current test is against a height of ~<2.5km
! This could be replaced by a scale height if a suitable method
! for determining a sensible height was possible from profile/cumulus
! depth information available in this routine

            If (ntml(i,j)  >   kcucheck(ii)                             &
     &             .and. nlcl(i,j)  <=  kcucheck(ii)-2) then

              grad_cld = ABS( QW(ii,kcucheck(ii)) - QW(ii,nlcl_c(ii)) ) &
     &          /( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)) )
            Else
              grad_cld = ABS( QW(ii,ntml_c(ii)) - QW(ii,nlcl_c(ii)) )   &
     &          /( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)) )
            End If

            grad_sub   =  ABS( QW(ii,nlcl_c(ii)) - QW(ii,k_plume(ii)) ) &
     &           /( z_full_c(ii,nlcl_c(ii)) - z_full_c(ii,k_plume(ii)) )

            If (grad_cld  >   1.10*grad_sub) then
!-----------------------------------------------------------------------
!     Not well mixed, however it is possible that the depth of a well
!     mixed boundary layer has increased but not yet been mixed yet so
!     test gradient from next level down.
!     Note typical cumulus profiles are expected to have a fairly
!     uniform q profile from the surface to the cloud base and then a
!     decreasing profile of q above this in the cloud. Typical the
!     decreasing gradient from the cloud base to 2.5km will be the
!     order of > 1.10 the below cloud value.
!-----------------------------------------------------------------------

! test against a height ~ 2.5km

              If (ntml_c(ii)  <=  kcucheck(ii)) then
              grad_cld = ABS( QW(ii,ntml_c(ii)-1) - QW(ii,nlcl_c(ii)) ) &
     &           /( z_full_c(ii,ntml_c(ii)-1) - z_full_C(ii,nlcl_c(ii)) )
              End If

              If ( grad_cld  >   1.10*grad_sub) then
!-----------------------------------------------------------------------
!      Diagnose a cumulus layer
!-----------------------------------------------------------------------
                cumulus(i,j) = .TRUE.
              End If

            Else

! Diagnosed well-mixed, but now check that LCL hasn't risen or fallen
! and not yet been mixed (so could have been erroneously identIfied as
! well-mixed)

! First check using level below (recalculate grad_sub)

              If (nlcl_c(ii) - k_plume(ii)  >=  2) then

                 grad_sub = ABS( QW(ii,nlcl(i,j)-1) - QW(ii,k_plume(ii)) ) &
     &           /( z_full_c(ii,nlcl_c(ii)-1) - z_full_c(ii,k_plume(ii)) )

                 If ( grad_cld  >   1.10*grad_sub) then
                   cumulus(i,j) =.TRUE.
                 End If

              End If

! If still diagnosing well-mixed, check using level above
! (recalculate grad_cld)

              If (.not. cumulus(i,j) ) then

               If (ntml_c(ii)  >   kcucheck(ii)                             &
     &             .and. nlcl_c(ii)  <=  kcucheck(ii)-2) then

                grad_cld = ABS( QW(ii,kcucheck(ii)) - QW(ii,nlcl_c(ii)+1) ) &
     &           /( z_full_c(ii,kcucheck(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               Else
                grad_cld = ABS( QW(ii,ntml_c(ii)) - QW(ii,nlcl_c(ii)+1) )   &
     &            /( z_full_c(ii,ntml_c(ii)) - z_full_c(ii,nlcl_c(ii)+1) )
               End If

               If ( grad_cld  >   1.10*grad_sub) then
                 cumulus(i,j) =.TRUE.
               End If

              End If   ! not cumulus
            End If     ! test on cloud gradient
          End If       ! test on cloud top height
        End If         ! tests on nlcl
      End Do           ! ii loop 

!-----------------------------------------------------------------------
!      Check that a cumulus layer has not been erroneously diagnosed in
!      a deep cloudy region
!      As the above checks are done on the total water rather than q it
!      is possible the conditions can be met in areas where the level of
!      prognostic qcl or qcf is high. The type of mistake is only
!      thought to occur over land.
!-----------------------------------------------------------------------
      Do ii=1, nunstable
        i = index_i(ii)   
        j = index_j(ii)   
        K=nlcl(i,j)

        If ( Land_MASK(i,j) .and. cumulus(i,j) .and.                    &
     &                                 ntpar(i,j)  <   MBL ) then
          Do while( K  <=  ntpar(i,j) .and. cloud_fraction(i,j,K)       &
     &                                            >=  SC_CFTOL )
            K = K + 1
          End Do
          If (K  ==  ntpar(i,j)+1) cumulus(i,j) = .false.
        End If
      End Do       ! ii loop 

!-----------------------------------------------------------------------
! Original shallow diagnosis no congestus diagnosis
!-----------------------------------------------------------------------
!CDIR NODEP
        Do ii=1,nunstable
          i = index_i(ii) 
          j = index_j(ii) 
          If ( cumulus(i,j) ) then

!-----------------------------------------------------------------------
!       If cumulus has been diagnosed, determine whether it is shallow
!       or deep convection
!-----------------------------------------------------------------------
! Conditions for shallow convection 
! 
!   top of parcel ascent < 2500. or T (top of parcel) > TM
!-----------------------------------------------------------------------

            If ( z_full_c(ii,ntpar(i,j))  <=  2500.0 .OR.               &
      &            T(ii,ntpar(i,j))  >=  TM ) then

              If (l_wtest) then
! Only shallow if descending air
                If (w_avg(ii).lt.0.0) then
                  L_shallow(i,j) = .TRUE.
                End If
                  
              Else
                L_shallow(i,j) = .TRUE.
              End If
            End If

!-----------------------------------------------------------------------
!      Set mixed layer depth to z_lcl
!-----------------------------------------------------------------------
          If (P_LCL(ii)  <   (P_theta_lev_c(ii,nlcl(i,j)+1))) then
!-----------------------------------------------------------------------
!      If LCL is diagnosed in the upper half of the layer set z_lcl to
!      the height of the upper layer interface
!      (in code above LCL is always set to the lower interface).
!-----------------------------------------------------------------------
            nlcl(i,j)    = nlcl(i,j)+1
            z_lcl(i,j)   = z_half_c(ii,nlcl(i,j)+1)
            z_lcl_uv(i,j)= z_full_c(ii,nlcl(i,j))
          End If
          zh(i,j)   = z_lcl(i,j)
          ntml(i,j) = nlcl(i,j)


!      If cumulus has been diagnosed but delthvu is negative, reset
!      cumulus and L_shallow to FALSE but leave zh and NTML at LCL
!      Need undilute CAPE to be > 0

          If (delthvu(i,j)<=  0.0 ) then

            cumulus(i,j)   = .false.
            L_shallow(i,j) = .false.

          End If

        Else      ! not cumulus

!-----------------------------------------------------------------------
!      If not cumulus, reset parameters to within bl_levels
!-----------------------------------------------------------------------
          If (ntml(i,j)  >   MBL) then
            ntml(i,j)  = MBL
            ntpar(i,j) = MBL
            zh(i,j)    = z_half(i,j,MBL+1)
            zhpar(i,j) = zh(i,j)
          End If
          If (nlcl(i,j)  >   MBL) then
            nlcl(i,j)    = MBL
            z_lcl(i,j)   = zh(i,j)
            z_lcl_uv(i,j)= z_full_c(ii,MBL-1)
          End If

        End If        ! test on cumulus

      End Do          ! ii loop  

!=======================================================================
! Option ?. New diagnosis - yet to be written  
!=======================================================================

      Else 
        write(6,*) 'Choice of convective diagnosis not allowed ',icvdiag
        error = 2       ! should force model to stop
!=======================================================================
!    End of choice of diagnosis code
!=======================================================================
      End If

!-----------------------------------------------------------------------
 9999  CONTINUE  ! Branch for error exit.

      RETURN
      END SUBROUTINE CONV_DIAG
