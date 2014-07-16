#if defined(A05_4A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Mid-level convection scheme
!

      SUBROUTINE  MID_CONV(nbl,nlev,ntra,n_cca_lev,npnts,trlev,         &
     &                       bland,W_MAX,exner_layer_centres,           &
     &                       exner_layer_boundaries,flg_up_flx,         &
     &                       flg_up_flx_half,                           &
     &                       flg_entr_up, flg_detr_up, flg_dwn_flx,     &
     &                       flg_entr_dwn,flg_detr_dwn,                 &
     &                       L_calc_dxek, L_q_interact,                 &
     &                       L_tracer,midtrig,ntml,ntpar,               &
     &                       pstar,p_layer_centres,                     &
     &                       p_layer_boundaries,r_theta,r_rho,          &
     &                       z_theta, z_rho,rho,rho_theta,              &
     &                        q,q1_sd,t1_sd,th,                         &
     &                       timestep, u,v,recip_pstar,qse,ad_on,       &
     &                       mdet_on, ent_on, sdet_on, new_termc,       &
     &                       bulk_cf,cf_frozen,cf_liquid,qcf,           &
     &                       qcl,tracer,cape_out,cclwp,ccw,             &
     &                       dbcfbydt,dcffbydt,dcflbydt,dqbydt,dqcfbydt,&
     &                       dqclbydt,dthbydt,                          &
     &                       dubydt,dvbydt,dtrabydt,                    &
     &                       detrain_up,detrain_dwn,                    &
     &                       entrain_up,entrain_dwn,                    &
     &                       iccb,icct,lcca,lcclwp,                     &
     &                       lcbase,lctop,rain,snow,rain_3d,snow_3d,    &
     &                       up_flux,up_flux_half,                      &
     &                       dwn_flux,tcw,l_mid_all,cca_2d,             &
     &                       rbuoy_p_out,the_out,thp_out,qe_out,qp_out, &
     &                       uw_mid,vw_mid                              &
     &                       )

!
! Purpose:
!   Mid level convection scheme - works on all points.
!
!   Called by GLUE_CONV.
!
! Current owners of code: R A Stratton
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
!
      Use cv_cntl_mod, Only:                                            &
          lcv_ccrad

      Use cv_run_mod, Only:                                             &
          l_mom, l_eman_dd, cape_opt, cape_ts_w, cape_min,              &
          w_cape_limit, cape_timescale, mid_cmt_opt,                    &
          mid_cnv_pmin

      IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Auguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!

      integer, intent(in) :: nbl      ! No. of boundary layer levels

      integer, intent(in) :: nlev     ! No. of model layers

      integer, intent(in) :: ntra     ! No. of tracer fields

      integer, intent(in) :: n_cca_lev! No. of convective cloud
                                      ! amount levels (1 for 2D,
                                      ! nlevs for 3D)

      integer, intent(in) :: npnts    ! No. of deep convection points

      integer, intent(in) :: trlev    ! No. of model levels on which
                                      ! tracers are included

      logical, intent(in) :: bland(npnts) ! Land/sea mask

      real, intent(in)    :: exner_layer_centres(npnts,0:nlev) !Exner

      real, intent(in)    :: exner_layer_boundaries(npnts,0:nlev)
                                      ! Exner at half level above
                                      ! exner_layer_centres

      logical, intent(in) :: flg_up_flx ! STASH flag for updraught
                                      ! mass flux

      logical, intent(in) :: flg_up_flx_half ! STASH flag for updraught
                                      ! mass flux on rho levels

      logical, intent(in) :: flg_entr_up ! STASH flag for updraught
                                      ! entrainment

      logical, intent(in) :: flg_detr_up ! STASH flag for updraught
                                      ! detrainment

      logical, intent(in) :: flg_dwn_flx ! STASH flag for downdraught
                                      ! mass flux

      logical, intent(in) :: flg_entr_dwn ! STASH flag for downdraught
                                      ! entrainment

      logical, intent(in) :: flg_detr_dwn ! STASH flag for downdraught
                                      ! detrainment

      logical, intent(in) :: L_calc_dxek ! Switch for calculation of
                                      ! condensate increment

      logical, intent(in) :: L_q_interact ! Switch allows overwriting
                                      ! parcel variables when
                                      ! calculating condensate incr.

      logical, intent(in) :: L_tracer ! Switch for inclusion of tracers

      integer, intent(in) :: midtrig(npnts) ! Lowest trigger level
                                      ! for convection

      integer, intent(in) :: ntml(npnts) ! Top level of surface mixed
                                      ! layer defined relative to
                                      ! theta,q grid

      integer, intent(in) :: ntpar(npnts) ! Top level of initial parcel
                                      ! ascent in BL scheme defined
                                      ! relative to theta,q grid

      real, intent(in)    :: pstar(npnts) ! Surface pressure (Pa)

      real, intent(in)    :: p_layer_centres(npnts,0:nlev) ! Pressure
                                      ! (Pa)


      real, intent(in)    :: p_layer_boundaries(npnts,0:nlev) ! Pressure
                                      ! at half level above
                                      ! p_layer_centres (Pa)

      real, intent(in)    :: r_theta(npnts,0:nlev) ! radius of theta
                                                 ! levels (m)

      real, intent(in)    :: r_rho(npnts,nlev) ! radius of rho levels (m)

      real, intent(in)    :: z_theta(npnts,nlev) ! height of theta
                                                 ! levels (m)

      real, intent(in)    :: z_rho(npnts,nlev) ! height of rho levels
                                               ! (m)

      real, intent(in)    :: rho(npnts,nlev) ! density on rho levels
                                             ! (kg/m3)

      real, intent(in)    :: rho_theta(npnts,nlev) ! density on theta levels
                                                   ! (kg/m3)

      real, intent(in)    :: q(npnts,nlev) ! Model mixing ratio (kg/kg)

      real, intent(in)    :: q1_sd(npnts) ! Standard deviation of
                                      ! turbulent flucts. of layer 1 q
                                      ! (kg/kg)

      real, intent(in)    :: t1_sd(npnts) ! Standard deviation of
                                      ! turbulent flucts. of layer 1
                                      ! temp. (K)

      real, intent(in)    :: th(npnts,nlev) !Model potential
                                      ! temperature (K)

      real, intent(in)    :: timestep ! Model timestep (s)

      real, intent(in)    :: u(npnts,nlev) !Model U field (m/s)

      real, intent(in)    :: v(npnts,nlev) !Model V field (m/s)

      real, intent(in) :: W_MAX(npnts)    !  max w in column
                                  ! for use in scale dependent
                                  ! cape timescale

      real, intent(in) :: recip_pstar(npnts) ! Reciprocal of pstar array

      real,intent(in) :: qse(npnts,nlev) ! Saturation mixing ratio of
                                      ! cloud environment (kg/kg)
!
      integer, intent(in) :: ad_on      !flag for adaptive detrainment
                                        !0 = off, 1 = on

      integer, intent(in) :: mdet_on    !flag for adaptive detrainment
                                        !0 = off, 1 = on

      integer, intent(in) :: ent_on     !flag for adaptive entrainment
                                        !0 = off, 1 = on

      integer, intent(in) :: sdet_on    !flag for smoothed forced detrainment
                                        !0 = off, 1 = on

      integer, intent(in) :: new_termc  !flag for simplified
                                        !termination of convection
                                        !0 = off, 1 = on



!
! Arguments with intent INOUT:
!

      real, intent(inout) :: bulk_cf(npnts,nlev) ! Bulk total cloud
                                      ! volume ( )

      real, intent(inout) :: cf_frozen(npnts,nlev) ! Frozen water cloud
                                      ! volume ( )

      real, intent(inout) :: cf_liquid(npnts,nlev) ! Liq water cloud
                                      ! volume ( )

      real, intent(inout) :: qcf(npnts,nlev) ! Ice condensate mix ratio
                                      ! (kg/kg)

      real, intent(inout) :: qcl(npnts,nlev) ! Liq condensate mix ratio
                                      ! (kg/kg)

      real, intent(inout) :: tracer(npnts,trlev,ntra) !Model tracer
                                      ! fields (kg/kg)

      real, intent(inout) :: tcw(npnts) ! Total condensed water(kg/m2/s)


!
! Arguments with intent OUT:
!

      real, intent(out) :: cape_out(npnts) ! Saved convective available
                                      ! potential energy for diagnostic
                                      ! output (J/kg)

      real, intent(out) :: cclwp(npnts) ! Condensed water path (k/m2)

      real, intent(out) :: ccw(npnts,nlev) ! Convective cloud liquid

                                      ! water on model levels (g/kg)
      real, intent(out) :: dbcfbydt(npnts,nlev) ! Increments to
                                      ! total cld volume due to
                                      ! convection(/s)

      real, intent(out) :: dcffbydt(npnts,nlev) ! Increments to ice
                                      ! cloud volume due to convection
                                      ! (/s)

      real, intent(out) :: dcflbydt(npnts,nlev) ! Increments to liq
                                      ! cloud volume due to convection
                                      ! (/s)

      real, intent(out) :: dqbydt(npnts,nlev) ! Increments to q due to
                                      ! convection (kg/kg/s)

      real, intent(out) :: dqcfbydt(npnts,nlev) ! Increments to ice
                                      ! condensate due to convection
                                      ! (kg/kg/s)

      real, intent(out) :: dqclbydt(npnts,nlev) ! Increments to liq
                                      ! condensate due to convection
                                      ! (kg/kg/s)

      real, intent(out) :: dthbydt(npnts,nlev) ! Increments to potential
                                      ! temp. due to convection (K/s)

      real, intent(out) :: dubydt(npnts,nlev+1) ! Increments to U due
                                      ! to CMT (m/s2)

      real, intent(out) :: dvbydt(npnts,nlev+1) ! Increments to V due
                                      ! to CMT (m/s2)

      real, intent(out) :: dtrabydt(npnts,nlev,ntra) !Increment to
                                      ! tracer convection (kg/kg/s)

      real, intent(out) :: detrain_up(npnts,nlev) ! Fractional
                                      ! detrainment rate into updraughts
                                      ! (Pa/s)

      real, intent(out) :: detrain_dwn(npnts,nlev) ! Fractional
                                      ! detrainment rate into
                                      ! downdraughts (Pa/s)

      real, intent(out) :: entrain_up(npnts,nlev) ! Fractional
                                      ! entrainment rate into updraughts
                                      ! (Pa/s)

      real, intent(out) :: entrain_dwn(npnts,nlev) ! Fractional
                                      ! entrainment rate into
                                      ! downdraughts (Pa/s)

      integer, intent(out) :: iccb(npnts) ! Convective cloud base
                                      ! level (m)

      integer, intent(out) :: icct(npnts) ! Convective cloud top
                                      ! level (m)

      real, intent(out) :: lcca(npnts) ! Lowest conv. cloud amt. (%)

      real, intent(out) :: lcclwp(npnts) ! Condensed water path for
                                      ! lowest conv. cld. level (kg/m2)

      integer, intent(out) :: lcbase(npnts) ! Lowest conv. cloud base
                                      ! level (m)

      integer, intent(out) :: lctop(npnts) ! Lowest conv. cloud top
                                      ! level (m)

      real, intent(out) :: rain(npnts) ! Surface convective rainfall
                                      ! (kg/m2/s)

      real, intent(out) :: snow(npnts) ! Surface convective snowfall
                                      ! (kg/m2/s)

      real, intent(out) :: rain_3d(npnts,nlev) ! Convective rainfall flux
                                      ! (kg/m2/s)

      real, intent(out) :: snow_3d(npnts,nlev) ! Convective snowfall flux
                                      ! (kg/m2/s)

      real, intent(out) :: up_flux(npnts,nlev) ! Updraught mass flux
                                      ! (Pa/s)

      real, intent(out) :: up_flux_half(npnts,nlev)
                                      ! Updraught mass flux
                                      ! on rho levels (Pa/s)

      real, intent(out) :: dwn_flux(npnts,nlev) ! Downdraught mass
                                      ! flux (Pa/s)

      logical, intent(out) :: l_mid_all(npnts)
                                    ! Points where mid level convection
                                      ! occrus at some level
      real,intent(inout) :: cca_2d(npnts) !2D convective cloud amount(%)

!
! Adaptive detrainment output variables
!
      real, intent(out) :: rbuoy_p_out(npnts,nlev)  !buoyancy excess

      real, intent(out) :: the_out(npnts,nlev)    !th_E in parcel
                                                  !routine

      real, intent(out) :: thp_out(npnts,nlev)    !th_P in parcel
                                                  !routine

      real, intent(out) :: qe_out(npnts,nlev)     !q_E in parcel
                                                  !routine

      real, intent(out) :: qp_out(npnts,nlev)     !q_P in parcel
                                                  !routine

! CMT diagnostics
      real, intent(out) :: &
      uw_mid(npnts,nlev)   & ! U component of stress from mid-level convection
                             ! (kg/m/s2)
     ,vw_mid(npnts,nlev)     ! V component of stress from mid-level convection
                             ! (kg/m/s2)


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

!
! Adaptive detrainment output variables
!
      real :: rbuoy_p_here(npnts)       !buoyancy excess

      real :: the_here(npnts)         !th_E in parcel routine

      real :: thp_here(npnts)         !th_P in parcel routine

      real :: qe_here(npnts)          !q_E in parcel routine

      real :: qp_here(npnts)          !q_P in parcel routine

      real :: rbuoy_p_old(npnts)      !buoyancy excess on previous
                                      !level

      real :: zk(npnts)               !heights for use in calc

      real :: zkp12(npnts)            !of moist static energy

      real :: zkp1(npnts)

      integer :: index1(npnts),index2(npnts)

      integer :: ncposs               ! No. of points which may convect

      integer :: nconv                ! No. of convecting points

      real :: amdetk(npnts)           ! Mixing detrainment coefficient
                                      ! at level k multiplied by
                                      ! appropriate layer thickness

      real :: cape(npnts)             ! Convective available potential
                                      ! energy (J/kg)

      real :: dcpbydt(npnts)          ! Rate of change of cape (J/kg/s)

      real :: depth(npnts)            ! Depth of convective cloud (m)

      real :: delexkp1(npnts)         ! Difference in exner ratio
                                      ! across layer k+1

      real :: dqsthk(npnts)           ! Gradient of saturation mixing
                                      ! ratio of cloud environment with
                                      ! theta in layer k (kg/kg/K)

      real :: dqsthkp1(npnts)         ! Gradient of saturation mixing
                                      ! ratio of cloud environment with
                                      ! theta in layer k+1 (kg/kg/K)

      real :: eminds(npnts)           ! Minimum buoyancy for convection
                                      ! to initiate from level k
                                      ! (Kelvin)

      real :: ekp14(npnts)            ! Entrainment coefficients at
                                      ! level k+1/4 multiplied by
                                      ! appropriate layer thickness
                                      !(Dimensionless)

      real :: ekp34(npnts)            ! Entrainment coefficients at
                                      ! level k+3/4 multiplied by
                                      ! appropriate layer thickness
                                      !(dimensionless)

      real :: ekm14(npnts)            ! Entrainment coefficients at
                                      ! level k-1+1/4 multiplied by
                                      ! appropriate layer thickness
                                      !(dimensionless)

      real :: exk(npnts)              ! Exner ratio at layer k

      real :: exkp1(npnts)            ! Exner ratio at layer k+1

      real :: flxmax(npnts)           ! Maximum initial convective
                                      ! mass flux (Pa/s)

      real :: flx_init(npnts)         ! Initial mass flux at cloud base
                                      ! (Pa/s)

      real :: flx_init_new(npnts)     ! flx_init scaled to destroy cape
                                      ! over timescale cape_timescale (Pa/s)

      real :: flxmax_init(npnts)      ! Maximum possible initial mass
                                      ! flux (limited to the mass in
                                      ! the initial convecting layer
                                      ! in Pa/s)

      real :: max_cfl(npnts)          ! Max cfl ratio over a convecting
                                      ! layer

!      real :: p_lcl(npnts)            ! Pressure at LCL (Pa)

      real :: precip(npnts,nlev)      ! Amount of precip from each layer
                                      ! from each layer (kg/m/s)

!      real :: pt(npnts)               ! Temporary store for P in calc.
                                      ! of sat. mixing ratio. (Pa)

      real :: pk(npnts)               ! Pressure at midpoint of layer
                                      ! k (Pa)

      real :: pkp1(npnts)             ! Pressure at midpoint of layer
                                      ! k+1 (Pa)

      real :: delpk(npnts)            ! Pressure difference over layer
                                      ! k (Pa)

      real :: delpkp1(npnts)          ! Pressure difference over layer
                                      ! k+1 (Pa)

      real :: delpkp12(npnts)         ! Pressure difference between
                                      ! layers k and k+1 (Pa)

      real :: delp_uv_k(npnts)        ! Pressure difference across uv
                                      ! layer k (Pa)

      real :: delp_uv_kp1(npnts)      ! Pressure difference across uv
                                      ! layer k+1 (Pa)

      real :: rhum(npnts)             ! Dummy relative humidity
                                      ! (only used on shallow points)

!      real :: tt(npnts)               ! Temporary store for T in calc.
                                      ! of saturation mixing ratio. (K)

!      real :: ttkm1(npnts)            ! Temporary store for T in layer
                                      ! k-1 for use in freezing level
                                      ! calc. for anvil. (K)

      logical :: L_shallow(npnts)     ! Dummy variable (=.F.)

      logical :: L_mid(npnts)         ! Dummy variable (=.T.)

      logical :: cumulus(npnts)       ! Dummy variable (=.F.)

      logical :: bgmk(npnts)          ! Mask for points where parcel in
                                      ! layer k is saturated

      logical :: bwater(npnts,2:nlev) ! Mask for points at which
                                      ! condensate is liquid

      logical :: bwk(npnts)            !mask for liquid condensate on k
      logical :: bwkp1(npnts)          !mask for liquid condensate on
                                       !k+1

      logical :: bterm(npnts)         ! Mask for points which have
                                      ! stopped convecting

      logical :: bconv(npnts)         ! Mask for points at which
                                      ! convection is occurring

      logical :: bcposs(npnts)        ! Mask for points passing
                                      ! initial stability test

      logical :: binit(npnts)         ! Mask for points initiating

!
! Parcel variables
!

      real :: qpi(npnts)              ! Initial parcel mixing ratio
                                      !(kg/kg)

      real :: qp(npnts,nlev)          ! Parcel mixing ratio (kg/kg)

      real :: thpi(npnts)             ! Initial parcel potential temp.
                                      !(K)

      real :: thp(npnts,nlev)         ! Parcel potential temp (K)

      real :: up(npnts,nlev)          ! Parcel U (m/s)

      real :: vp(npnts,nlev)          ! Parcel V (m/s)

      real :: trap(npnts,nlev,ntra)   ! Tracer content of parcel
                                      ! (kg/kg)
      real :: expi(npnts)
                                      ! Initial parcel exner pressure

      real :: xpk(npnts,nlev)         ! Parcel cloud water (kg/kg)

      real :: flx(npnts,nlev)         ! Parcel massflux (Pa/s)

      real :: xsbmin_v(npnts,nlev)    ! Minmum parcel buoyancy excess

      real :: thpixs_v(npnts,nlev)    ! Theta parcel excess (K)

      real :: qpixs_v(npnts,nlev)     ! Q parcel excess(kg/kg)

      real :: dpmin(npnts)            ! work array for parcel excess cal

!
! PC2
!
      real :: qclp(npnts,nlev)         ! Parcel liquid condensate mixing
                                      ! ratio in layer k (kg/kg)

      real :: qcfp(npnts,nlev)         ! Parcel frozen condensate mixing
                                      ! ratio in layer k (kg/kg)
!
! Parameters
!

      real, parameter :: CFL_LIMIT = 1.0 ! Max CFL ratio allowed

!
! CMT variables
!

      integer :: kterm(npnts)         ! Level index for termination of

      real :: eflux_u_ud(npnts)       ! Vertical eddy flux of momentum
                                      ! due to UD at top of layer
                                      ! (Pa m/s2)

      real :: eflux_v_ud(npnts)       ! Vertical eddy flux of momentum
                                      ! due to UD at bottom of layer
                                      ! (Pa m/s2)

      real :: flxkp12(nlev,npnts+1)   ! Mass flux on half level (Pa/s)

      Logical :: l_mom_gk             ! true if Gregory-Kershaw CMT


!
! Cape scaling/closure variables
!

      integer :: start_lev(npnts)     !Level at which convection
                                      !initiates

      integer :: start_lev3c(npnts)   ! PC2 Compressed convection
                                      ! initiation level

      integer :: det_lev(npnts)       ! Level at which split final
                                      ! detrainment last occurred

      integer :: nterm                ! No. of points where conv.
                                      ! has terminated

      integer :: index_nterm(npnts)   ! Index for points where conv.
                                      ! has terminated

      real :: tempnum                 ! Temporary variable for storage

      real :: scale_f(npnts)          ! scaling factor

      real :: dthef(npnts)            ! Theta increment from convection
                                      ! in model level at which split
                                      ! final detrainment last occurred
                                      ! (K/s)

      real :: dqf(npnts)              ! Specific humidity increment
                                      ! from convection in model level
                                      ! at which split final detrainment
                                      ! last occurred (kg/kg/s)

      real :: dqclf(npnts)            ! As dqf but for qcl (kg/kg/s)
      real :: dqcff(npnts)            ! As dqf but for qcf (kg/kg/s)
      real :: dcflf(npnts)            ! As dqf but for cfl (/s)
      real :: dcfff(npnts)            ! As dqf but for cff (/s)
      real :: dbcff(npnts)            ! As dqf but for bcf (/s)

      real :: duef(npnts)             ! As for dthef but for U
                                      !increments (m/s2)

      real :: dvef(npnts)             ! As for dthef but for V
                                      !increments (m/s2)

      real :: dtraef(npnts,ntra)      ! As for dthef but for tracer
                                      ! increments (kg/kg/s)

      real :: cape_ts_new             ! Used as variable in RH-based
                                      ! closure

      real :: relh(npnts)             ! RH integral (average when
                                      ! convection terminates)

      real :: dptot(npnts)            ! Delta P integral


!
! Downdraught scheme variables
!

      integer :: npossdd              ! Max. no. of downdraughts
                                      ! possible

      integer :: nnodd                ! No. of downdraughts not possible

      integer :: index_possdd(npnts)  ! Index of downdraughts possible

      integer :: index_nodd(npnts)    ! Index of downdraughts not
                                      ! possible

      real :: deltap_cld              ! Pressure thickness of convective
                                      ! cloud (Pa)

!
! Local compressed arrays
!

      logical :: bconv_c2(npnts)

      logical :: bgmkp1_c(npnts), bgmkp1_c2(npnts) ! Mask for points
                                      ! where parcel in layer k+1
                                      ! is saturated

      logical :: bwk_c(npnts), bwk_c2(npnts) ! bwater mask in layer k

      logical :: bwkp1_c(npnts), bwkp1_c2(npnts) ! bwater mask in layer
                                      ! k+1

      real :: deltak_c2(npnts)        ! Parcel forced detrainment rate
                                      ! in layer k multiplied by
                                      ! appropriate layer thickness

      real :: dqek_c2(npnts)          ! Increment to q due to
                                      ! convection in layer k (kg/kg)

      real :: dqekp1_c2(npnts)        ! Increment to q due to
                                      ! convection in layer k+1 (kg/kg)

      real :: dthek_c2(npnts)         ! Increment to potential temp.
                                      ! due to convection in layer k

      real :: dthekp1_c2(npnts)       ! Increment to potential temp.
                                      ! due to convection in layer k+1

      real :: dtraek_c2(npnts,ntra)   ! Increment to model tracer due
                                      ! to conv. at level k (kg/kg/s)

      real :: dtraekp1_c2(npnts,ntra) ! Increment to model tracer due
                                      ! to conv. at level k+1 (kg/kg/s)

      real :: duek_c2(npnts)          ! Increment to model U in layer k
                                      ! due to CMT (m/s2)

      real :: duekp1_c2(npnts)        ! Increment to model U in layer
                                      ! k+1 due to CMT (m/s2)

      real :: dvek_c2(npnts)          ! Increment to model V in layer k
                                      ! due to CMT (m/s2)

      real :: dvekp1_c2(npnts)        ! Increment to model V in layer
                                      ! k+1 due to CMT (m/s2)

      real :: flxk_c(npnts), flxk_c2(npnts) !Parcel mass flux in layer k
                                      ! (Pa/s)

      real :: flxkp12_c2(npnts)       ! Half level mass flux (Pa/s)

      real :: prekp1_c2(npnts)        ! Precip. from parcel as it rises
                                      ! from layer k to k+1 (kg/m2/s)

      real :: qpk_c(npnts), qpk_c2(npnts) ! Parcel mixing ratio in
                                      ! layer k(kg/kg)

      real :: qpk(npnts)              !ad. entrain.

      real :: qpkp1_c(npnts), qpkp1_c2(npnts) ! Parcel mixing ratio
                                      ! in layer k+1 (kg/kg)

      real :: qek_c(npnts), qek_c2(npnts) ! Env. mixing ratio in
                                      ! layer k (kg/kg)

      real :: qek(npnts)             !for ad entrain.

      real :: qekp1_c(npnts), qekp1_c2(npnts) ! Env. mixing ratio in
                                      ! layer k+1 (kgkg-1)

      real :: qekp1(npnts)             !for ad entrain.

      real :: qsek_c2(npnts)          ! Saturation mixing ratio of
                                      ! cld. env. in layer k (kg/kg)

      real :: qsek(npnts)             !for ad entrain.

      real :: qsekp1_c(npnts), qsekp1_c2(npnts) ! Saturation mixing
                                      ! ratio of cld. env. in layer k+1
                                      ! (kg/kg)

      real :: qsekp1(npnts)             !for ad entrain.

      real :: thek_c(npnts), thek_c2(npnts) ! Env. potential temp
                                      ! in layer k (K)

      real :: thek(npnts)             !for ad entrain.

      real :: thekp1_c(npnts), thekp1_c2(npnts) ! Env. potential temp i
                                      ! in layer k (K)

      real :: thekp1(npnts)             !for ad entrain.

      real :: thpk_c(npnts), thpk_c2(npnts) ! Parcel potential temp
                                      ! in layer k (K)

      real :: thpk(npnts)             !for ad entrain.

      real :: thpkp1_c(npnts), thpkp1_c2(npnts)! Parcel potential temp
                                      ! in layer k (K)

      real :: traek_c(npnts,ntra), traek_c2(npnts,ntra) ! Tracer content
                                      ! cld. env. in layer k (kgkg-1)

      real :: traekp1_c(npnts,ntra), traekp1_c2(npnts,ntra) ! Tracer
                                      ! content of cloud env.
                                      ! in layer k+1 (kg/kg)

      real :: trapk_c(npnts,ntra), trapk_c2(npnts,ntra) ! Tracer cont.
                                      ! of parcel in layer k (kg/kg)

      real :: trapkp1_c(npnts,ntra), trapkp1_c2(npnts,ntra) ! Tracer
                                  ! cont.of parcel in layer k+1 (kg/kg)

      real :: rbuoy_c(npnts), rbuoy_c2(npnts)! Buoyancy of parcel
                                      ! at k+1 (Kelvin)

      real :: uek_c(npnts), uek_c2(npnts) ! Model U field on layer k
                                      ! (m/s)

      real :: uekp1_c(npnts), uekp1_c2(npnts)! Model U field on layer
                                      ! k+1 (m/s)

      real :: vek_c(npnts), vek_c2(npnts) ! Model V field on layer k
                                      ! (m/s)

      real :: vekp1_c(npnts), vekp1_c2(npnts) ! Model V field on layer
                                      ! k+1 (m/s)

      real :: upk_c(npnts), upk_c2(npnts) ! Parcel U in layer k
                                      ! after entrainment (m/s)


      real :: upkp1_c(npnts), upkp1_c2(npnts) ! Parcel U in layer k+1
                                      ! after entrainment (m/s)

      real :: vpk_c(npnts), vpk_c2(npnts) ! Parcel V in layer k
                                      ! after entrainment (m/s)

      real :: vpkp1_c(npnts), vpkp1_c2(npnts) ! Parcel V in layer k+1
                                      ! after entrainment (m/s)

      real :: xsqkp1_c(npnts), xsqkp1_c2(npnts) ! Excess water vapour
                                      ! in parcel at k+1 (kg/kg)

      real :: qclek_c(npnts), qclek_c2(npnts) ! Environment liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qclekp1_c(npnts), qclekp1_c2(npnts) ! Environment liquid
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qcfek_c(npnts), qcfek_c2(npnts) ! Environment frozen
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qcfekp1_c(npnts), qcfekp1_c2(npnts) ! Environment frozen
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qclpk_c(npnts), qclpk_c2(npnts) ! Parcel liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qclpkp1_c(npnts), qclpkp1_c2(npnts) ! Parcel liquid
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qcfpk_c(npnts), qcfpk_c2(npnts) ! Parcel frozen
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qcfpkp1_c(npnts), qcfpkp1_c2(npnts) ! Parcel frozen
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: cflek_c2(npnts),cflekp1_c2(npnts)
                                      ! Environment liquid water
                                      ! cloud volume ( )

      real :: cffek_c2(npnts),cffekp1_c2(npnts)
                                      ! Environment frozen water
                                      ! cloud volume ( )

      real :: bcfek_c2(npnts),bcfekp1_c2(npnts)
                                      ! Environment bulk total
                                      ! cloud volume ( )

      real :: dqclek_c2(npnts),dqclekp1_c2(npnts)
                                      ! Environment increments
                                      ! to liquid condensate mixing
                                      ! ratio to convection (kg/kg/s)

      real :: dqcfek_c2(npnts),dqcfekp1_c2(npnts)
                                      ! Environment increments
                                      ! to frozen condensate mixing
                                      ! ratio to convection (kg/kg/s)

      real :: dcflek_c2(npnts),dcflekp1_c2(npnts)
                                      ! Environment increments
                                      ! to liquid water cloud volume due
                                      ! to convection (/s)

      real :: dcffek_c2(npnts),dcffekp1_c2(npnts)
                                      ! Environment increments
                                      ! to frozen water cloud volume due
                                      ! to convection (/s)

      real :: dbcfek_c2(npnts),dbcfekp1_c2(npnts)
                                      ! Environment increments
                                      ! to bulk total cloud volume due
                                      ! to convection (/s)

      real :: amdetk_c2(npnts)
      logical :: bgmk_c2(npnts)
      logical :: bland_c2(npnts)
      logical :: blowst_c2(npnts)
      logical :: bterm_c2(npnts)
      real :: cape_c2(npnts)
      real :: cca_2d_c2(npnts)
      real :: cclwp_c2(npnts)
      real :: ccw_c2(npnts)
      logical :: cumulus_c(npnts), cumulus_c2(npnts)
      real :: dcpbydt_c2(npnts)
      real :: delexkp1_c2(npnts)
      real :: delpk_c2(npnts)
      real :: delpkp1_c2(npnts)
      real :: delp_uv_k_c2(npnts)
      real :: delp_uv_kp1_c2(npnts)
      real :: depth_c2(npnts)
      real :: dptot_c2(npnts)
      real :: dqsthkp1_c2(npnts)
      real :: dqsthk_c2(npnts)
      real :: eflux_u_ud_c2(npnts)
      real :: eflux_v_ud_c2(npnts)
      real :: ekp14_c(npnts),ekp14_c2(npnts)
      real :: ekp34_c(npnts),ekp34_c2(npnts)
      real :: eminds_c(npnts)
      real :: exk_c2(npnts)
      real :: exkp1_c(npnts),exkp1_c2(npnts)
      real :: expi_c2(npnts)
      integer :: icct_c2(npnts)
      integer :: iccb_c2(npnts)
      real :: lcclwp_c2(npnts)
      integer :: lctop_c2(npnts)
      integer :: lcbase_c2(npnts)
      real :: lcca_c2(npnts)
      logical :: L_shallow_c2(npnts)
      logical :: L_mid_c2(npnts)
      real :: max_cfl_c2(npnts)
      real :: pk_c(npnts),pk_c2(npnts)
      real :: pkp1_c(npnts),pkp1_c2(npnts)
      real :: pstar_c2(npnts)
      real :: q1_sd_c2(npnts)
      real :: qpi_c2(npnts)
      real :: qpixs_v_c2(npnts)
      real :: rbuoy_p_here_c2(npnts)
      real :: the_here_c2(npnts)
      real :: thp_here_c2(npnts)
      real :: qe_here_c2(npnts)
      real :: qp_here_c2(npnts)
      real :: rbuoy_p_old_c2(npnts)
      real :: relh_c2(npnts)
      real :: tcw_c2(npnts)
      real :: thpi_c2(npnts)
      real :: thpixs_v_c2(npnts)
      real :: t1_sd_c2(npnts)
      real :: xpk_c(npnts),xpk_c2(npnts)
      real :: xsbmin_v_c2(npnts)

      integer ::        &
       nmid             & ! total number of points where mid-level
                          ! convection occurs
      ,nmax_layers        ! Maximum number of allow mid-level layers

! Arrays required by Emanuel downdraughts
!
      Integer ::                                                        &
     &  kterm_mid(npnts)                                                &
                             ! termination level of highest mid level
     &, kterm_max


!
! Loop counters
!

      integer :: i,j,k,ktra,kt,kk

!
! Model constants:
!

#include "parxs.h"
#include "c_epslon.h"
#include "c_r_cp.h"
#include "xsbmin.h"
#include "mparb.h"
#include "c_lheat.h"
#include "massfc.h"
#include "c_0_dg_c.h"
#include "c_g.h"
#include "c_mass.h"
#include "delthst.h"
#include "wcape_spd_up.h"


!-----------------------------------------------------------------------
! Decide whether Gregory-Kershaw CMT scheme is to be used if l_mom is true.
! The Gregory-Kershaw scheme requires calculations in the main plume ascent
! loop whereas the alternative diffusive scheme is called after the plume
! calculation.

      If (mid_cmt_opt == 1 ) then
         l_mom_GK = .false.       ! Use diffusive scheme
      Else
         l_mom_GK = l_mom         ! Use Gregory-Kershaw scheme
      End If

!
!initialise SCM diagnostics
!
      Do k = 1,nlev
        Do i = 1,npnts
          rbuoy_p_out(i,k)=0.0
          the_out(i,k)=th(i,k)
          thp_out(i,k)=th(i,k)
          qe_out(i,k)=q(i,k)
          qp_out(i,k)=q(i,k)
        End Do
      End Do
!
!
! Initialise logicals
!

      Do i = 1,npnts
!        blowst(i)    = .true.
        bterm(i)     = .false.
        bconv(i)     = .false.
        bcposs(i)    = .false.
        cumulus(i)   = .false.
        L_shallow(i) = .false.
        L_mid(i)     = .true.
        L_mid_all(i) = .false.

! Zero bottom level mass flux array as never set by levels loop.

        flx(i,1) = 0.0
      End Do

!-----------------------------------------------------------------------
! 1.0  Create saturation mixing ratio arrays
!-----------------------------------------------------------------------

!
! Re-calculate XSBMIN and THPIXS constants based on layer
! thickness (Pa).
!
      Do k = 1,nlev-1
          Do i = 1,npnts
            dpmin(i) = min( ((p_layer_centres(i,k) -                    &
     &                         p_layer_centres(i,k+1))/5000.),1.0)

            xsbmin_v(i,k) = dpmin(i) *0.2

            If (midtrig(i)  ==  ntml(i)                                 &
     &                     .or. midtrig(i)  ==  ntml(i)+1) then

              thpixs_v(i,k) = dpmin(i) * THPIXS_DEEP
              qpixs_v(i,k)  = QPIXS_DEEP
            else
              thpixs_v(i,k) = dpmin(i) * THPIXS_MID
              qpixs_v(i,k)  = QPIXS_MID
            End If
          End Do
       End Do  ! nlev

!
! Set bwater=.true. on points where water will condense rather than
! ice.
! SUBROUTINE FLAG_WET
! UM Documentation paper 27, section (2B)
!

! DEPENDS ON: flag_wet
      Call FLAG_WET(bwater,th,exner_layer_centres,npnts,npnts,nlev)


!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 2.1  Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
!      increment arrays and cca
!-----------------------------------------------------------------------

      Do k = 1,nlev
        Do i = 1,npnts
          precip(i,k) = 0.0
          ccw(i,k)    = 0.0
         dthbydt(i,k) = 0.0
         dqbydt(i,k)  = 0.0
        End Do
      End Do

      If (L_mom) then

        l_mom_gk = .true.        ! use Gregory-Kershaw CMT

        Do k = 1,nlev+1
          Do i = 1,npnts
            dubydt(i,k) = 0.0
            dvbydt(i,k) = 0.0
          End Do
        End Do
        Do k = 1,nlev
          Do i = 1,npnts
            uw_mid(i,k) = 0.0
            vw_mid(i,k) = 0.0
          End Do
        End Do
      End If  ! L_mom

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,nlev
            Do i = 1,npnts
              dtrabydt(i,k,ktra) = 0.0
            End Do
          End Do
          Do i = 1,npnts
            dtraef(i,ktra) = 0.0
          End Do
        End Do
      End If  ! L_tracer

!-----------------------------------------------------------------------
! 2.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------

      If (flg_up_flx) then
        Do k = 1,nlev
          Do i = 1,npnts
            up_flux(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_up_flx_half) then
        Do k = 1,nlev
          Do i = 1,npnts
            up_flux_half(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_dwn_flx) then
        Do k = 1,nlev
          Do i = 1,npnts
            dwn_flux(i,k) = 0.0
          End Do
        End Do
      End If
! Now required in all cases
      Do k = 1,nlev
        Do i = 1,npnts
          entrain_up(i,k) = 0.0
        End Do
      End Do

      If (flg_detr_up) then
        Do k = 1,nlev
          Do i = 1,npnts
            detrain_up(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_entr_dwn) then
        Do k = 1,nlev
          Do i = 1,npnts
            entrain_dwn(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_detr_dwn) then
        Do k = 1,nlev
          Do i = 1,npnts
            detrain_dwn(i,k) = 0.0
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------


      If (lcv_ccrad) Then
        ! Zeroes values to remove mid-level dependence on shallow and deep
        ! cloud
        Do i = 1,npnts
          iccb   (i) = 0
          icct   (i) = 0
          tcw    (i) = 0.0
          cca_2d (i) = 0.0
        End Do
      End If      ! lcv_ccrad

      Do i = 1,npnts
        cclwp(i)  = 0.0
        lcca(i)   = 0.0
        lctop(i)  = 0
        lcbase(i) = 0
        lcclwp(i) = 0.0
      End Do

!-----------------------------------------------------------------------
! 2.4  Initialise gridbox mean diagnostics (now in glue)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2.5  Initialise diagnostics for scaling and closure calculations
!-----------------------------------------------------------------------

      Do i = 1,npnts
        flx_init(i)     = 0.0
        flx_init_new(i) = 0.0
        cape(i)      = 0.0
        cape_out(i)  = 0.0
        dcpbydt(i)   = 0.0
        max_cfl(i)   = 0.0
        det_lev(i)   = 0
        start_lev(i) = 0
        dthef(i)      = 0.0
        dqf(i)        = 0.0
        dqclf(i)        = 0.0
        dqcff(i)        = 0.0
        dcflf(i)        = 0.0
        dcfff(i)        = 0.0
        dbcff(i)        = 0.0
        duef(i)       = 0.0
        dvef(i)       = 0.0
        relh(i)       = 0.0
        dptot(i)      = 0.0

!-----------------------------------------------------------------------
! 2.6  Initialise eddy flux arrays for updraught
!-----------------------------------------------------------------------

        eflux_u_ud(i) = 0.0
        eflux_v_ud(i) = 0.0

!-----------------------------------------------------------------------
! 2.7  Initialise surface precipitation arrays
!-----------------------------------------------------------------------

        rain(i) = 0.0
        snow(i) = 0.0

      End Do

!-----------------------------------------------------------------------
! 2.8 Initialise PC2 arrays
!-----------------------------------------------------------------------

      Do k = 1,nlev
        Do i = 1,npnts
          dqclbydt(i,k)  = 0.0
          dqcfbydt(i,k)  = 0.0
          dbcfbydt(i,k)  = 0.0
          dcflbydt(i,k)  = 0.0
          dcffbydt(i,k)  = 0.0
        End Do
      End Do

!
! Set rhum to an arbitrary value since it is not used in LAYER_CN
! for mid-level convection. Therefore value not level dependent.
!

        Do i = 1,npnts
          rhum(i) = 0.0
        End Do

!
! set SCM adaptive diagnostics for level k = 1
!
      Do i = 1,npnts
        rbuoy_p_out(i,1) = 0.0
        the_out(i,1) = th(i,1)
        thp_out(i,1) = th(i,1)
        qe_out(i,1) = q(i,1)
        qp_out(i,1) = q(i,1)
        rbuoy_p_old(i) = 0.0
      End Do
!initialise ekm14
      Do i =1, npnts
         ekm14(i) =0.0
      End do

!Initialise adaptive entrainment variables
!Note that unlike p_layer_boundaries, where k indexing is offset
!by one compared to the dynamics numbering, z retains the numbering
!convention for dynamics variables i.e. for theta levels, k->k
!and for rho levels k+1/2 -> k+1

      Do i = 1, npnts
        thek(i)=th(i,1)
        qek(i)=q(i,1)
        qsek(i)=qse(i,1)
        thpk(i)=thp(i,1)
        qpk(i)=qp(i,1)
        zk(i) = z_theta(i,1)
      End Do
      Do i = 1, npnts
        thekp1(i)=th(i,2)
        qekp1(i)=q(i,2)
        qsekp1(i)=qse(i,2)
        bwk(i)=bwater(i,2)
        bwkp1(i)=bwater(i,2)
        zkp12(i)=z_rho(i, 2)
        zkp1(i)=z_theta(i, 2)
      End Do

!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------

      Do k = 2,nlev-1

!
!  Initialise adaptive diagnostics for this level
!
      Do i = 1,npnts
        rbuoy_p_here(i) =0.0
        the_here(i) = th(i,k)
        thp_here(i) = th(i,k)
        qe_here(i) = q(i,k)
        qp_here(i) = q(i,k)
        rbuoy_p_here_c2(i) =0.0
        the_here_c2(i) = 0.0
        thp_here_c2(i) = 0.0
        qe_here_c2(i) = 0.0
        qp_here_c2(i) = 0.0
      End Do
!Initialise adaptive entrainment variables
      Do i = 1, npnts
        thek(i)=th(i,k)
        qek(i)=q(i,k)
        qsek(i)=qse(i,k)
        thekp1(i)=th(i,k+1)
        qekp1(i)=q(i,k+1)
        qsekp1(i)=qse(i,k+1)
        bwk(i)=bwater(i,k)
        bwkp1(i)=bwater(i,k+1)
 !Note that unlike p_layer_boundaries, where k indexing is offset
!by one compared to the dynamics numbering, z retains the numbering
!convention for dynamics variables i.e. for theta levels, k->k
!and for rho levels k+1/2 -> k+1
!check this with Rachel
        zk(i) = z_theta(i,k)
        zkp12(i)=z_rho(i, k+1)
        zkp1(i)=z_theta(i, k+1)
      End Do
      If(k  ==  2) then      !set to environmental values
        Do i = 1,npnts
          thpk(i)=th(i,2)
          qpk(i)=q(i,2)
        End Do
      Else
        Do i = 1,npnts
          thpk(i)=thp(i,k)
          qpk(i)=qp(i,k)
        End Do
      End if
!
! Initialise binit i.e. set as no convection initialised in layer
! at start of this levels calculations
!

        Do i = 1,npnts
          binit(i) = .false.
        End Do

!-----------------------------------------------------------------------
! 3.1  Calculate layer dependent constants (pressure,
!      layer thickness, entrainment coefficients, detrainment
!      coefficients)
!-----------------------------------------------------------------------

! DEPENDS ON: layer_cn
        Call LAYER_CN(k,npnts,npnts,nlev,                               &
     &                exner_layer_boundaries,                           &
     &                exner_layer_centres,                              &
     &                p_layer_boundaries,                               &
     &                p_layer_centres,                                  &
     &                pstar,pk,pkp1,delpk,delpkp1,                      &
     &                delpkp12,ekp14,ekp34,amdetk,exk,exkp1,            &
     &                delexkp1,delp_uv_k,delp_uv_kp1,recip_pstar,       &
     &                rhum,L_shallow,ntml,ntpar                         &
     &                , cumulus, mdet_on,ent_on                         &
     &                , bconv                                           &
     &                ,thek,qek, qsek, thekp1,qekp1,qsekp1              &
     &                ,thpk, qpk                                        &
     &                ,bwk,bwkp1,ekm14                                  &
     &                ,zk, zkp12, zkp1                                  &
     &                )

! Set ekm14 for next pass through loop
         Do i = 1, npnts
           ekm14(i) = ekp14(i)
         End Do
!
! Calculate dqs/dth for layers k and k+1 (subroutine DQS_DTH)
!

        If (k == 2) then
! DEPENDS ON: dqs_dth
          Call DQS_DTH(dqsthk,k,th(1,k),qse(1,k),exk,npnts)
        else
          Do i = 1,npnts
            dqsthk(i) = dqsthkp1(i)
          End Do
        End If

! DEPENDS ON: dqs_dth
        Call DQS_DTH(dqsthkp1,k+1,th(1,k+1),qse(1,k+1),exkp1,npnts)

!
! Set other grid dependent constants
!

        Do i = 1,npnts

!
! Maximum initial convective mass flux
!

          flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)

!
! Minimum buoyancy for convection to start from level k
!

          eminds(i) = MPARB * delpkp12(i) * recip_pstar(i)

        End Do
!
! Set initial parcel properties (theta,q,tracer,momentum) if convection
! is not occurring at level k
! qclp and qcfp zero at this point but will add an initial
! value at cloud base
!
        Do i = 1,npnts
          If ( .not. bconv(i)) then
            expi(i)  = exk(i)
            xpk(i,k)  = 0.0
            qclp(i,k) = 0.0
            qcfp(i,k) = 0.0
            flx(i,k) = 0.0
            bgmk(i)  = .false.
            depth(i) = 0.0
            thpi(i)  = th(i,k) + thpixs_v(i,k)
            thp(i,k)  = thpi(i)
            qpi(i)   = q(i,k) + qpixs_v(i,k)
            qp(i,k)   = qpi(i)
            If (L_mom_gk) then
              up(i,k)   = u(i,k)
              vp(i,k)   = v(i,k)
            End If
          End If  !not bconv
        End Do

        If (L_tracer) then
          Do ktra=1,ntra
            Do i = 1,npnts
              If ( .not. bconv(i)) then
                 trap(i,k,ktra)  = tracer(i,k,ktra)
              End If  !not bconv
            End Do    ! npnts
          End Do
        End If

!
! Carry out initial test to see if convection is possible from layer
! k to k+1. Set bcposs = .T. if
! 1. the point was convecting (bconv = .T.) and did not terminate
! in the previous layer  OR
! 2. k > ntml+1 and the stability is low enough at k+1
!  (At levels above nbl this is the same as midtrig(i)=k)

        Do i = 1,npnts
          bcposs(i) = bconv(i) .or.                                     &
                      (p_layer_centres(i,k) > mid_cnv_pmin .and.        &
                      k > ntml(i) + 1 .and.                             &
                      (( th(i,k) * (1.0 + C_VIRTUAL * q(i,k)) -         &
                      th(i,k+1) * (1.0 + C_VIRTUAL * q(i,k+1)) +        &
                      DELTHST + MAX(0.0,(q(i,k)-qse(i,k+1))) *          &
                      (LC/(CP * exkp1(i)))) > 0.)) .or.                 &
                      (k == ntml(i) .and. ntml(i) == nbl-1)

        End Do  ! npnts

!
! Calculate number of points which may convect (ncposs) and
! set compression indices (index1)
!

        ncposs = 0
        Do i = 1,npnts
          If (bcposs(i)) then
            ncposs = ncposs + 1
            index1(ncposs) = i
          End If
        End Do

!
! Compress points where convection may occur
!

        If (ncposs  >   0) then
          Do i = 1,ncposs
            thek_c(i)   = th(index1(i),k)
            thekp1_c(i) = th(index1(i),k+1)
            qek_c(i)    = q(index1(i),k)
            qekp1_c(i)  = q(index1(i),k+1)
            qsekp1_c(i) = qse(index1(i),k+1)
            thpk_c(i)   = thp(index1(i),k)
            qpk_c(i)    = qp(index1(i),k)
            xpk_c(i)    = xpk(index1(i),k)
            bwk_c(i)    = bwater(index1(i),k)
            bwkp1_c(i)  = bwater(index1(i),k+1)
            pk_c(i)     = pk(index1(i))
            pkp1_c(i)   = pkp1(index1(i))
            ekp14_c(i)  = ekp14(index1(i))
            ekp34_c(i)  = ekp34(index1(i))
            exkp1_c(i)  = exkp1(index1(i))
            bgmkp1_c(i) = bgmk(index1(i)) ! bgmk into lift_par,
            cumulus_c(i) = .false.        ! bgmkp1_c out
            eminds_c(i) = eminds(index1(i))
          End Do
          If (L_q_interact) then
           Do i = 1,ncposs
            qclek_c(i)    = qcl(index1(i),k)
            qclekp1_c(i)  = qcl(index1(i),k+1)
            qcfek_c(i)    = qcf(index1(i),k)
            qcfekp1_c(i)  = qcf(index1(i),k+1)
            qclpk_c(i)    = qclp(index1(i),k)
            qcfpk_c(i)    = qcfp(index1(i),k)
           End Do
          End If
          If (L_mom_gk) then
           Do i = 1,ncposs
            uek_c(i)    = u(index1(i),k)
            uekp1_c(i)  = u(index1(i),k+1)
            vek_c(i)    = v(index1(i),k)
            vekp1_c(i)  = v(index1(i),k+1)
            upk_c(i)    = up(index1(i),k)
            vpk_c(i)    = vp(index1(i),k)
           End Do
          Endif
          If (L_tracer) then
          Do ktra = 1,ntra
            Do i = 1,ncposs
              traek_c(i,ktra) = tracer(index1(i),k,ktra)
              traekp1_c(i,ktra) = tracer(index1(i),k+1,ktra)
              trapk_c(i,ktra) = trap(index1(i),k,ktra)
            End Do
          End Do
          End If
        End If  ! ncposs>0


!-----------------------------------------------------------------------
! 3.2  Lift parcel from layer k to layer k+1
!-----------------------------------------------------------------------

! DEPENDS ON: lift_par
        Call LIFT_PAR(ncposs,npnts,thpkp1_c,qpkp1_c,xsqkp1_c,           &
                 bgmkp1_c,bwkp1_c,bwk_c,thpk_c,qpk_c,xpk_c,thekp1_c,    &
                 qekp1_c,thek_c,qek_c,qsekp1_c,                         &
                 qclpkp1_c,qclpk_c,qclekp1_c,qclek_c,l_q_interact,      &
                 qcfpkp1_c,qcfpk_c,qcfekp1_c,qcfek_c,                   &
                 pk_c,pkp1_c,exkp1_c,ekp14_c,ekp34_c,L_mom_gk,upkp1_c,  &
                 vpkp1_c,upk_c,vpk_c,uek_c,uekp1_c,vek_c,vekp1_c,       &
                 L_tracer,ntra,trapkp1_c,trapk_c,traekp1_c,             &
                 traek_c,cumulus_c)

!
! Loop over points which may convect (ncposs)
!
! NEC compiler directive
!CDIR NODEP

        Do i = 1,ncposs

!
! Calculate buoyancy (virt. potential temp.) of parcel in layer k+1
!

          rbuoy_c(i) = thpkp1_c(i) * (1.0 + C_VIRTUAL *qpkp1_c(i))      &
     &                  - thekp1_c(i) * (1.0 + C_VIRTUAL *qekp1_c(i))

!
! Allow parcel to convect at midtrig or above if
!1. it has not convected in the past (bconv=.F.,bterm=.F.) AND
!2. k >= midtrig AND
!3. the buoyancy of parcel at k+1 > min buoyancy required for
!   convection + XSBMIN
! Flag convecting points with logical array bconv
!

          If ( .not. bconv(index1(i)) .and. .not. bterm(index1(i)) .and.&
     &   ( k  >   midtrig(index1(i)).or.                                &
     &     (k == ntml(index1(i)).and.ntml(index1(i)) == nbl-1)) .and.   &
     &           rbuoy_c(i)  >   (eminds_c(i) + XSBMIN))then
            bconv(index1(i)) = .true.
            start_lev(index1(i)) = k
            binit(index1(i)) = .true.
            l_mid_all(index1(i)) = .true. ! mid level convection

          End If


!
! Set parcel mass flux (UM Documentation paper 27, section 1.5)
! If mass flux out of the initial layer is greater than the mass flux
! of the layer over the timestep then limit mass flux to mass of layer.
!

          If (bconv(index1(i)).and.k == start_lev(index1(i))) then
            flxk_c(i) = 1.0e-3*pstar(index1(i)) * (D_MID + C_MID *      &
     &                    pstar(index1(i)) * ((rbuoy_c(i) - XSBMIN)     &
     &                     / delpkp12(index1(i))))
            If (flxk_c(i)  >   flxmax(index1(i))) then
              flxk_c(i) = flxmax(index1(i))
            End If
            If (flg_up_flx) then
              up_flux(index1(i),k)=flxk_c(i)
            End If

!
! Write compressed mass flux back to full array
!

            flx(index1(i),k) = flxk_c(i)

!
! At ntml set mixing detrainment rate equal to zero
! Store diagnostics linked to initial convective mass flux for
! calculation of final closure.
!
            amdetk(index1(i))      = 0.0
            flx_init(index1(i))    = flxk_c(i)
            flxmax_init(index1(i)) = flxmax(index1(i))

            If (L_q_interact) then
!
!           Initialize QCLP(*,K) and QCFP(*,K) at start level and
!           perform the Parcel Lift again for such points: needed
!           for PC2.   Duplicates code from
!           Convection Parcel Lifting Scheme, LIFPAR.
!           Assume "if k eq ntml" equiv "if bwork(i,3) in convec4a
!           Use non-compressed arrays for variables not updated
!           in or since LIFT_PAR

              qclpk_c(i)  = qcl(index1(i),k) +                          &
     &                     (1. + ekp14(index1(i))) *                    &
     &                     (qcl(index1(i),k+1) - qcl(index1(i),k))
              qclpkp1_c(i) = ( qclpk_c(i) +                             &
     &               (ekp14(index1(i)) * qcl(index1(i),k)) +            &
     &               (ekp34(index1(i)) *  (1. + ekp14(index1(i)))       &
     &                 * qcl(index1(i),k+1)) )/                         &
     &          ( (1. + ekp14(index1(i))) * (1. + ekp34(index1(i))) )
!
              qcfpk_c(i) = qcf(index1(i),k) +                           &
     &                    (1. + ekp14(index1(i))) *                     &
     &                    (qcf(index1(i),k+1) - qcf(index1(i),k))
              qcfpkp1_c(i) = ( qcfpk_c(i) +                             &
     &               (ekp14(index1(i)) * qcf(index1(i),k)) +            &
     &               (ekp34(index1(i)) * (1. + ekp14(index1(i)))        &
     &                 * qcf(index1(i),k+1)) )/                         &
     &          ( (1. + ekp14(index1(i))) * (1. + ekp34(index1(i))) )

            End If ! l_q_interact

          End If

        End Do  !ncposs loop

! L_q_interact_if0:
      If (L_Q_INTERACT) then
! ----------------------------------------------------------------------
!     Follow-on calculation from QCLP(*,K) and QCFP(*,K) initialization.
!       Duplicates code from Convection Parcel Lifting Scheme, LIFPAR.
! ----------------------------------------------------------------------
! If test for convection different from dp/sh conv.
! Use bconv and k eq start_lev as above

        Do I=1, ncposs
          If (bconv(index1(i)).and.k == start_lev(index1(i))) then
! ----------------------------------------------------------------------
!       CURRENTLY MIXED PHASE PARCEL IS FORBIDDEN. MELT OR FREEZE THE
!       ENTRAINED LAYER CLOUD AND ADJUST PARCEL TEMPERATURE ACCORDINGLY.
! ----------------------------------------------------------------------

            If (bwater(index1(i),k+1) .and. qcfpkp1_c(i)  >   0.0) then
              qclpkp1_c(i) = qclpkp1_c(i) + qcfpkp1_c(i)
              thpkp1_c(i)  = ( thpkp1_c(i) - ( qcfpkp1_c(i) * LF /      &
     &                                   (CP * exkp1(index1(i))) ) )
              qcfpkp1_c(i) = 0.0
            Else If (.not. bwater(index1(i),k+1) .and.                  &
     &               qclpkp1_c(i)  >   0.0) then
              qcfpkp1_c(i) = qclpkp1_c(i) + qcfpkp1_c(i)
              thpkp1_c(i)  = ( thpkp1_c(i) + ( qclpkp1_c(i) * LF /      &
     &                                   (CP * exkp1(index1(i))) ) )
              qclpkp1_c(i) = 0.0
            End If
       End If
      End Do
      End If  ! L_q_interact_if0

!
! Calculate number of points which are convecting  (nconv)
! set compression indices (index2).
!

        nconv = 0
        Do i = 1,ncposs
          If (bconv(index1(i))) then
            nconv = nconv + 1
            index2(nconv) = i
          End If
        End Do

!
! Second compression to form arrays of length nconv to be passed
! into CONVEC2
!

!
! Input variables to CONVEC2
!

        If (nconv  >   0) then
          Do i = 1,nconv
            thek_c2(i)   = thek_c(index2(i))
            thekp1_c2(i) = thekp1_c(index2(i))
            qek_c2(i)    = qek_c(index2(i))
            qekp1_c2(i)  = qekp1_c(index2(i))
            dqsthkp1_c2(i) = dqsthkp1(index1(index2(i)))
            qsekp1_c2(i)   = qsekp1_c(index2(i))
            pstar_c2(i)    = pstar(index1(index2(i)))
            thpkp1_c2(i)   = thpkp1_c(index2(i))
            qpkp1_c2(i)    = qpkp1_c(index2(i))
            xsqkp1_c2(i) = xsqkp1_c(index2(i))
            rbuoy_c2(i)  = rbuoy_c(index2(i))
            qsek_c2(i)   = qse(index1(index2(i)),k)
            dqsthk_c2(i) = dqsthk(index1(index2(i)))
            thpi_c2(i)   = thpi(index1(index2(i)))
            qpi_c2(i)    = qpi(index1(index2(i)))
            expi_c2(i)   = expi(index1(index2(i)))
            bconv_c2(i)  = bconv(index1(index2(i)))
            bwk_c2(i)    = bwk_c(index2(i))
            bwkp1_c2(i)  = bwkp1_c(index2(i))
            bgmkp1_c2(i) = bgmkp1_c(index2(i))
            bland_c2(i)  = bland(index1(index2(i)))
            blowst_c2(i) = binit(index1(index2(i)))
            L_shallow_c2(i) = .false.
            L_mid_c2(i)     = .true.
            cumulus_c2(i)   = .false.
            ekp14_c2(i)  = ekp14_c(index2(i))
            ekp34_c2(i)  = ekp34_c(index2(i))
            amdetk_c2(i)    = amdetk(index1(index2(i)))
            pk_c2(i)     = pk_c(index2(i))
            pkp1_c2(i)   = pkp1_c(index2(i))
            exk_c2(i)    = exk(index1(index2(i)))
            exkp1_c2(i)  = exkp1_c(index2(i))
            delexkp1_c2(i)    = delexkp1(index1(index2(i)))
            delpk_c2(i)       = delpk(index1(index2(i)))
            delpkp1_c2(i)     = delpkp1(index1(index2(i)))
            delp_uv_k_c2(i)   = delp_uv_k(index1(index2(i)))
            delp_uv_kp1_c2(i) = delp_uv_kp1(index1(index2(i)))
            t1_sd_c2(i)     = t1_sd(index1(index2(i)))
            q1_sd_c2(i)     = q1_sd(index1(index2(i)))
            xsbmin_v_c2(i)  = xsbmin_v(index1(index2(i)),k)
          End Do
          If (L_q_interact) then
           Do i=1,nconv
!           PC2 variables
            qclek_c2(i)   = qclek_c(index2(i))
            qclekp1_c2(i) = qclekp1_c(index2(i))
            qcfek_c2(i)   = qcfek_c(index2(i))
            qcfekp1_c2(i) = qcfekp1_c(index2(i))
            qclpk_c2(i)   = qclpk_c(index2(i))
            qcfpk_c2(i)   = qcfpk_c(index2(i))
!           (PC2) Compress input cloud fields
            cflek_c2(i)   = cf_liquid(index1(index2(i)),k)
            cflekp1_c2(i) = cf_liquid(index1(index2(i)),k+1)
            cffek_c2(i)   = cf_frozen(index1(index2(i)),k)
            cffekp1_c2(i) = cf_frozen(index1(index2(i)),k+1)
            bcfek_c2(i)   = bulk_cf(index1(index2(i)),k)
            bcfekp1_c2(i) = bulk_cf(index1(index2(i)),k+1)
!           Compress convective base indicator
            start_lev3c(I) = start_lev(index1(index2(i)))
           End Do
          End If
          If (L_mom_gk) then
           Do i=1,nconv
            uek_c2(i)    = uek_c(index2(i))
            uekp1_c2(i)  = uekp1_c(index2(i))
            vek_c2(i)    = vek_c(index2(i))
            vekp1_c2(i)  = vekp1_c(index2(i))
            upkp1_c2(i)    = upkp1_c(index2(i))
            vpkp1_c2(i)    = vpkp1_c(index2(i))
           End Do
          End If
          If (L_tracer) then
          Do ktra = 1,ntra
            Do i = 1,nconv
              traek_c2(i,ktra)   = traek_c(index2(i),ktra)
              traekp1_c2(i,ktra) = traekp1_c(index2(i),ktra)
              trapkp1_c2(i,ktra) = trapkp1_c(index2(i),ktra)
            End Do
          End Do
          End If

!
! Input/output variables to/from CONVEC2
!

          Do i = 1,nconv
            thpk_c2(i)   = thpk_c(index2(i))
            qpk_c2(i)    = qpk_c(index2(i))
            xpk_c2(i)    = xpk_c(index2(i))
            flxk_c2(i)   = flx(index1(index2(i)),k)
            bgmk_c2(i)   = bgmk(index1(index2(i)))
            bterm_c2(i)  = .false.
            dthek_c2(i)  = dthbydt(index1(index2(i)),k)
            dqek_c2(i)   = dqbydt(index1(index2(i)),k)
            dthekp1_c2(i)  = dthbydt(index1(index2(i)),k+1)
            dqekp1_c2(i)   = dqbydt(index1(index2(i)),k+1)
            tcw_c2(i)     = tcw(index1(index2(i)))
            depth_c2(i)   = depth(index1(index2(i)))
            cclwp_c2(i)   = cclwp(index1(index2(i)))
            cape_c2(i)    = cape(index1(index2(i)))
            dcpbydt_c2(i) = dcpbydt(index1(index2(i)))
            relh_c2(i)    = relh(index1(index2(i)))
            rbuoy_p_here_c2(i)=rbuoy_p_here(index1(index2(i)))
            the_here_c2(i)=the_here(index1(index2(i)))
            thp_here_c2(i)=thp_here(index1(index2(i)))
            qe_here_c2(i)=qe_here(index1(index2(i)))
            qp_here_c2(i)=qp_here(index1(index2(i)))
            rbuoy_p_old_c2(i)=rbuoy_p_old(index1(index2(i)))
            dptot_c2(i)   = dptot(index1(index2(i)))
            thpixs_v_c2(i) = thpixs_v(index1(index2(i)),k)
            qpixs_v_c2(i)  = qpixs_v(index1(index2(i)),k)
            iccb_c2(i)     = iccb(index1(index2(i)))
            icct_c2(i)     = icct(index1(index2(i)))
            cca_2d_c2(i)   = cca_2d(index1(index2(i)))
            lcca_c2(i)     = lcca(index1(index2(i)))
            lcbase_c2(i)   = lcbase(index1(index2(i)))
            lctop_c2(i)    = lctop(index1(index2(i)))
            lcclwp_c2(i)   = lcclwp(index1(index2(i)))
          End Do
          If (L_q_interact) then
           Do i = 1,nconv
!           PC2 variables
            qclpkp1_c2(i) = qclpkp1_c(index2(i))
            qcfpkp1_c2(i) = qcfpkp1_c(index2(i))
!           Compress increment fields
            dqclek_c2(i) = dqclbydt(index1(index2(i)),k)
            dqclekp1_c2(i) = dqclbydt(index1(index2(i)),k+1)
            dqcfek_c2(i) = dqcfbydt(index1(index2(i)),k)
            dqcfekp1_c2(i) = dqcfbydt(index1(index2(i)),k+1)
            dcflek_c2(i) = dcflbydt(index1(index2(i)),k)
            dcflekp1_c2(i) = dcflbydt(index1(index2(i)),k+1)
            dcffek_c2(i) = dcffbydt(index1(index2(i)),k)
            dcffekp1_c2(i) = dcffbydt(index1(index2(i)),k+1)
            dbcfek_c2(i) = dbcfbydt(index1(index2(i)),k)
            dbcfekp1_c2(i) = dbcfbydt(index1(index2(i)),k+1)
           End Do
          End If
          If (L_mom_gk) then
            Do i = 1,nconv
              upk_c2(i)    = upk_c(index2(i))
              vpk_c2(i)    = vpk_c(index2(i))
              duek_c2(i)   = dubydt(index1(index2(i)),k)
              dvek_c2(i)   = dvbydt(index1(index2(i)),k)
              duekp1_c2(i)   = dubydt(index1(index2(i)),k+1)
              dvekp1_c2(i)   = dvbydt(index1(index2(i)),k+1)
              eflux_u_ud_c2(i) = eflux_u_ud(index1(index2(i)))
              eflux_v_ud_c2(i) = eflux_v_ud(index1(index2(i)))
            End Do
          End If
          If (L_tracer) then
          Do ktra = 1,ntra
            Do i = 1,nconv
             trapk_c2(i,ktra) = trapk_c(index2(i),ktra)
            dtraek_c2(i,ktra) = dtrabydt(index1(index2(i)),k,ktra)
            dtraekp1_c2(i,ktra) = dtrabydt(index1(index2(i)),k+1,ktra)
            End Do
          End Do
          End If


!-----------------------------------------------------------------------
! 3.3  Calculate the rest of the parcel ascent  and the effect of
!      convection on the large-scale atmosphere.
!
!      Subroutine CONVEC2
!
!      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
!-----------------------------------------------------------------------

! DEPENDS ON: convec2
        Call CONVEC2(nconv,npnts,nlev,k,thek_c2,thekp1_c2,qek_c2,       &
     &               qekp1_c2,qclek_c2,qclekp1_c2,qcfek_c2,qcfekp1_c2,  &
     &               cflek_c2,cflekp1_c2,cffek_c2,cffekp1_c2,bcfek_c2,  &
     &               bcfekp1_c2,qsekp1_c2,dqsthkp1_c2,pstar_c2,         &
     &               thpk_c2,qpk_c2,qclpk_c2,qcfpk_c2,thpkp1_c2,        &
     &               qpkp1_c2,qclpkp1_c2,qcfpkp1_c2,start_lev3c,        &
     &               xsqkp1_c2,rbuoy_c2,qsek_c2,dqsthk_c2,thpi_c2,      &
     &               qpi_c2,expi_c2,xpk_c2,flxk_c2,bwk_c2,bwkp1_c2,     &
     &               bgmkp1_c2,bgmk_c2,blowst_c2,                       &
     &               bland_c2,bterm_c2,depth_c2,prekp1_c2,              &
     &               dthek_c2,dqek_c2,dqclek_c2,dqcfek_c2,dcflek_c2,    &
     &               dcffek_c2,dbcfek_c2,dthekp1_c2,dqekp1_c2,          &
     &               dqclekp1_c2,dqcfekp1_c2,dcflekp1_c2,dcffekp1_c2,   &
     &               dbcfekp1_c2,bconv_c2,cca_2d_c2,iccb_c2,icct_c2,    &
     &               tcw_c2,ekp14_c2,ekp34_c2,amdetk_c2,pk_c2,pkp1_c2,  &
     &               exk_c2, exkp1_c2,delexkp1_c2,delpk_c2,delpkp1_c2,  &
     &               cclwp_c2,ccw_c2,lcca_c2,lcbase_c2,lctop_c2,        &
     &               lcclwp_c2,t1_sd_c2,                                &
     &               q1_sd_c2,L_mom_gk,uek_c2,uekp1_c2,vek_c2,vekp1_c2, &
     &               upk_c2,vpk_c2,                                     &
     &               upkp1_c2,vpkp1_c2,duek_c2,duekp1_c2,               &
     &               dvek_c2,dvekp1_c2,                                 &
     &               eflux_u_ud_c2,eflux_v_ud_c2,                       &
     &               delp_uv_k_c2, delp_uv_kp1_c2,                      &
     &               thpixs_v_c2,qpixs_v_c2,                            &
     &               xsbmin_v_c2,L_shallow_c2,L_mid_c2,                 &
     &               L_tracer,ntra,traek_c2,traekp1_c2,                 &
     &               trapk_c2,trapkp1_c2,dtraek_c2,dtraekp1_c2,cape_c2, &
     &               dcpbydt_c2,                                        &
     &               max_cfl_c2,timestep,rbuoy_p_here_c2,               &
     &               the_here_c2,thp_here_c2,qe_here_c2,qp_here_c2,     &
     &               rbuoy_p_old_c2, ad_on, sdet_on, new_termc,         &
     &               deltak_c2,l_calc_dxek,l_q_interact,                &
     &               flxkp12_c2,cumulus_c2,relh_c2,dptot_c2)

!
! Calculate fractional entrainment rate for level k.
! If convection has terminated (bterm=.T.) then set
! fractional entrainment rate for k+1 to zero.
!


         Do i = 1,nconv
           entrain_up(index1(index2(i)),k) = (1.0 - deltak_c2(i)) *     &
                      (1.0 - amdetk_c2(i)) * (ekp14_c2(i) + ekp34_c2(i) &
                      * (1.0 + ekp14_c2(i))) * flx(index1(index2(i)),k)
           If (bterm_c2(i)) then
             entrain_up(index1(index2(i)),k+1) = 0.0
           End If
         End Do


!
! Calculate fractional detrainment rate for level k
! (and k+1 if bterm=.T.)
!

        If (flg_detr_up) then
          Do i = 1,nconv
            detrain_up(index1(index2(i)),k) = -(amdetk_c2(i)            &
     &                    + deltak_c2(i) * (1.0 - amdetk_c2(i)))        &
     &                  * flx(index1(index2(i)),k)
            If (bterm_c2(i)) then
              detrain_up(index1(index2(i)),k+1) =                       &
     &                 -(1.0 - deltak_c2(i)) * flx(index1(index2(i)),k)
            End If
          End Do
        End If
        End If !if nconv>0

!
! Write CONVEC2 compressed output arrays back to full fields
!


        Do i = 1,npnts
          thp(i,k+1)    = 0.0
          qp(i,k+1)     = 0.0
          xpk(i,k+1)    = 0.0
          flx(i,k+1)    = 0.0
          depth(i)      = 0.0
          precip(i,k+1) = 0.0
          qclp(i,k+1)   = 0.0
          qcfp(i,k+1)   = 0.0
          bgmk(i)       = .false.
          bterm(i)      = .false.
        End Do

        If (L_tracer) then
          Do ktra = 1,ntra
            Do i = 1,npnts
              trap(i,k+1,ktra) = 0.0
            End Do
          End Do
        End If

        If (L_mom_gk) then
          Do i = 1,npnts
            up(i,k+1) = 0.0
            vp(i,k+1) = 0.0
          End Do
         End If

        If (nconv  >   0) then
        Do i = 1,nconv
          thp(index1(index2(i)),k+1) = thpkp1_c2(i)
          qp(index1(index2(i)),k+1)  = qpkp1_c2(i)
          xpk(index1(index2(i)),k+1) = xpk_c2(i)
          flx(index1(index2(i)),k+1) = flxk_c2(i)
          depth(index1(index2(i)))   = depth_c2(i)
          precip(index1(index2(i)),k+1) = prekp1_c2(i)
          bgmk(index1(index2(i)))    = bgmk_c2(i)
          bterm(index1(index2(i)))   = bterm_c2(i)
          dthbydt(index1(index2(i)),k) = dthek_c2(i)
          dqbydt(index1(index2(i)),k)  = dqek_c2(i)
          dthbydt(index1(index2(i)),k+1) = dthekp1_c2(i)
          dqbydt(index1(index2(i)),k+1)  = dqekp1_c2(i)
          cca_2d(index1(index2(i)))  = cca_2d_c2(i)
          tcw(index1(index2(i)))     = tcw_c2(i)
          iccb(index1(index2(i)))    = iccb_c2(i)
          icct(index1(index2(i)))    = icct_c2(i)
          cclwp(index1(index2(i)))   = cclwp_c2(i)
          lcca(index1(index2(i)))    = lcca_c2(i)
          lcbase(index1(index2(i)))  = lcbase_c2(i)
          lctop(index1(index2(i)))   = lctop_c2(i)
          lcclwp(index1(index2(i)))  = lcclwp_c2(i)
          ccw(index1(index2(i)),k+1) = ccw_c2(i)
          cape(index1(index2(i)))    = cape_c2(i)
          dcpbydt(index1(index2(i))) = dcpbydt_c2(i)
          rbuoy_p_here(index1(index2(i))) = rbuoy_p_here_c2(i)
          the_here(index1(index2(i))) = the_here_c2(i)
          thp_here(index1(index2(i))) = thp_here_c2(i)
          qe_here(index1(index2(i)))  = qe_here_c2(i)
          qp_here(index1(index2(i)))  = qp_here_c2(i)
        End Do
        Do i = 1,nconv
          max_cfl(index1(index2(i))) =                                  &
     &                   MAX(max_cfl(index1(index2(i))),max_cfl_c2(i))
          relh(index1(index2(i)))  = relh_c2(i)
          dptot(index1(index2(i))) = dptot_c2(i)
        End Do
        If (L_q_interact) then
         Do i = 1,nconv
          qclp(index1(index2(i)),k+1)   = qclpkp1_c2(i)
          qcfp(index1(index2(i)),k+1)   = qcfpkp1_c2(i)
          dqclbydt(index1(index2(i)),k)   = dqclek_c2(i)
          dqcfbydt(index1(index2(i)),k)   = dqcfek_c2(i)
          dcflbydt(index1(index2(i)),k)   = dcflek_c2(i)
          dcffbydt(index1(index2(i)),k)   = dcffek_c2(i)
          dbcfbydt(index1(index2(i)),k)   = dbcfek_c2(i)
          dqclbydt(index1(index2(i)),k+1) = dqclekp1_c2(i)
          dqcfbydt(index1(index2(i)),k+1) = dqcfekp1_c2(i)
          dcflbydt(index1(index2(i)),k+1) = dcflekp1_c2(i)
          dcffbydt(index1(index2(i)),k+1) = dcffekp1_c2(i)
          dbcfbydt(index1(index2(i)),k+1) = dbcfekp1_c2(i)
         End Do
        End If
        If (L_mom_gk) then
          Do i = 1,nconv
            up(index1(index2(i)),k+1) = upk_c2(i)
            vp(index1(index2(i)),k+1) = vpk_c2(i)
            dubydt(index1(index2(i)),k) = duek_c2(i)
            dvbydt(index1(index2(i)),k) = dvek_c2(i)
            dubydt(index1(index2(i)),k+1) = duekp1_c2(i)
            dvbydt(index1(index2(i)),k+1) = dvekp1_c2(i)
            eflux_u_ud(index1(index2(i))) = eflux_u_ud_c2(i)
            eflux_v_ud(index1(index2(i))) = eflux_v_ud_c2(i)
            flxkp12(k,index1(index2(i)))  = flxkp12_c2(i)
          End Do
        End If
        If  (L_tracer) then
          Do i = 1,nconv
            Do ktra = 1,ntra
              trap(index1(index2(i)),k+1,ktra)     = trapk_c2(i,ktra)
              dtrabydt(index1(index2(i)),k,ktra)   = dtraek_c2(i,ktra)
              dtrabydt(index1(index2(i)),k+1,ktra) = dtraekp1_c2(i,ktra)
            End Do
          End Do
        End If
        If (flg_up_flx) then
          Do i = 1,nconv
            up_flux(index1(index2(i)),k+1) = flxk_c2(i)
          End Do
        End If
        If (flg_up_flx_half) then
          Do i = 1,nconv
            up_flux_half(index1(index2(i)),k+1) =  flxkp12_c2(i)
          End Do
        End If
        End If     ! nconv >0

!
!   Write adaptive diagnostics for this level to full array for output
!
        Do i = 1,npnts
          rbuoy_p_out(i,k) = rbuoy_p_here(i)
          the_out(i,k) = the_here(i)
          thp_out(i,k) = thp_here(i)
          qe_out(i,k) = qe_here(i)
          qp_out(i,k) = qp_here(i)
        End Do

!
!Write rbuoy for previous level for next pass through loop
        Do i = 1,npnts
          rbuoy_p_old(i) = rbuoy_p_here(i)
        End Do


!-----------------------------------------------------------------------
! 3.4  Cape and CFL scaling - adjust initial mass flux so that cape is
!      removed by convection over timescale cape_timescale.
!-----------------------------------------------------------------------
!
! Set up integer nterm which is the total number of points where
! convection has terminated.
! Index to full array (npnts) with index_nterm
!

        nterm = 0
        Do i = 1,npnts
          If (bterm(i)) then
            nterm = nterm + 1
            index_nterm(nterm) = i
          End If
        End Do

        If(cape_opt  ==  0) then
!default 4a convection scheme - RH-based CAPE closure
          If (nterm >  0) then

! NEC compiler directive
!CDIR  NODEP
            Do j = 1,nterm
              i = index_nterm(j)
              If (dcpbydt(i)  >   0.0) then
! RH CAPE for mid level convection
                cape_ts_new =                                           &
     &          MIN(MAX(900.0*(1.0 - relh(i)/dptot(i))/0.1,60.0)        &
     &          ,cape_timescale)

                flx_init_new(i) = flx_init(i)*cape(i)                   &
     &                            /(cape_ts_new*dcpbydt(i))

                If (flx_init_new(i)  >   flxmax_init(i)) then
                  flx_init_new(i) = flxmax_init(i)
                End If

!
! Scale max_cfl with cape scale
!

                max_cfl(i) = max_cfl(i) * flx_init_new(i)/flx_init(i)
              else
                flx_init_new(i) = flx_init(i)
              End If  ! dcpbydt > 0
            End Do !nterm
          End If !nterm  >   0

        Else If(cape_opt  ==  1) then
! modified 4a convection scheme - RH-based CAPE closure, if thickness
! greater than 150hPa, timescale limited to timestep
! (operational version)

          If (nterm >  0) then

! NEC compiler directive
!CDIR  NODEP
            Do j = 1,nterm
              i = index_nterm(j)
              If (dcpbydt(i)  >   0.0) then
                If (p_layer_centres(i,start_lev(i))-                    &
     &            p_layer_centres(i,k)  >   15000.0)  then
!
! Only use RH cape if thickness of convective layer is
! greater than 150hPa
!
                  cape_ts_new =                                         &
     &              MIN(MAX(cape_timescale*(1.0-relh(i)/dptot(i))/0.4          &
     &              ,timestep)                                          &
     &              ,cape_timescale)
                Else
                  cape_ts_new = cape_timescale
                End If

                flx_init_new(i) = flx_init(i)*cape(i)                   &
     &                            /(cape_ts_new*dcpbydt(i))

                If (flx_init_new(i)  >   flxmax_init(i)) then
                  flx_init_new(i) = flxmax_init(i)
                End If

!
! Scale max_cfl with cape scale
!

                max_cfl(i) = max_cfl(i) * flx_init_new(i)/flx_init(i)
              else
                flx_init_new(i) = flx_init(i)
              End If  ! dcpbydt > 0
            End Do !nterm
          End If !nterm  >   0



        Else If(cape_opt  ==  2) then
! Switch off RH-based CAPE closure, no test on w_max
          If (nterm >  0) then

! NEC compiler directive
!CDIR  NODEP
            Do j = 1,nterm
              i = index_nterm(j)
              If (dcpbydt(i)  >   0.0) then
! RH CAPE for mid level convection
                cape_ts_new = cape_timescale

                flx_init_new(i) = flx_init(i)*cape(i)                   &
     &                            /(cape_ts_new*dcpbydt(i))

                If (flx_init_new(i)  >   flxmax_init(i)) then
                  flx_init_new(i) = flxmax_init(i)
                End If

!
! Scale max_cfl with cape scale
!

                max_cfl(i) = max_cfl(i) * flx_init_new(i)/flx_init(i)
              else
                flx_init_new(i) = flx_init(i)
              End If  ! dcpbydt > 0
            End Do !nterm
          End If !nterm  >   0
        End if  !end if cape_opt = 0,1 or 2


        If ( nterm > 0 .and. w_cape_limit < 1000.0 ) then
!  This section includes test on w_max

          IF (cape_opt  ==  3 ) THEN

!CDIR  NODEP
            Do j = 1, nterm
              i = index_nterm(j)
              If ( dcpbydt(i) > 0.0 ) then
! new denominator introduced at vn6.6
                if ( W_MAX(I) > w_cape_limit ) then
                    cape_ts_new =   CAPE_TS_w * w_cape_limit/           &
     &            (w_cape_limit + (W_MAX(I)-w_cape_limit)*wcape_fac)
                else
                  cape_ts_new = cape_timescale
                endif !  W_MAX(I) > w_cape_limit

                flx_init_new(i) = flx_init(i) * cape(i) /               &
     &                                    (cape_ts_new * dcpbydt(i))
                If (flx_init_new(i) > flxmax_init(i)) then
                  flx_init_new(i) = flxmax_init(i)
                End If ! flx_init_new(i) > flxmax_init(i)
! Scale max_cfl with cape scale
                max_cfl(i) = max_cfl(i) * flx_init_new(i)               &
     &          / flx_init(i)
              else
                flx_init_new(i) = flx_init(i)
              End If  ! dcpbydt(i) > 0.0
            End Do ! j = 1, nterm

          ELSE IF (cape_opt  ==  4 ) THEN

!CDIR  NODEP
            Do j = 1, nterm
              i = index_nterm(j)
              If ( dcpbydt(i) > 0.0 ) then

                if ( W_MAX(I) > w_cape_limit ) then
                  cape_ts_new =   CAPE_TS_w * w_cape_limit/W_MAX(I)
                else
                  cape_ts_new = CAPE_TIMESCALE * CAPE(I) / CAPE_min +          &
     &                        CAPE_TIMESCALE * EXP( - CAPE(I) / CAPE_min)
                endif !  W_MAX(I) > w_cape_limit

                flx_init_new(i) = flx_init(i) * cape(i) /               &
     &                                    (cape_ts_new * dcpbydt(i))
                If (flx_init_new(i) > flxmax_init(i)) then
                  flx_init_new(i) = flxmax_init(i)
                End If ! flx_init_new(i) > flxmax_init(i)
! Scale max_cfl with cape scale
                max_cfl(i) = max_cfl(i) * flx_init_new(i)               &
     &                   / flx_init(i)
              else
                flx_init_new(i) = flx_init(i)
              End If  ! dcpbydt(i) > 0.0
            End Do ! j = 1, nterm

          ELSE IF (cape_opt  ==  5 ) THEN

!CDIR  NODEP
            Do j = 1, nterm
              i = index_nterm(j)
              If ( dcpbydt(i) > 0.0 ) then

                if ( W_MAX(I) > w_cape_limit ) then
                  cape_ts_new =   CAPE_TS_w * w_cape_limit/W_MAX(I)
                else
                  if ( relh(i) / dptot(i) >= 0.75 ) then
                    CAPE_TS_NEW = CAPE_TIMESCALE *                             &
     &                        ( 0.2373 / (relh(i) / dptot(i))**5)
                  else
                    CAPE_TS_NEW = CAPE_TIMESCALE
                  endif ! relh(i) / dptot(i) >= 0.75
                endif !  W_MAX(I) > w_cape_limit

                flx_init_new(i) = flx_init(i) * cape(i) /               &
     &                                    (cape_ts_new * dcpbydt(i))
                If (flx_init_new(i) > flxmax_init(i)) then
                  flx_init_new(i) = flxmax_init(i)
                End If ! flx_init_new(i) > flxmax_init(i)
! Scale max_cfl with cape scale
                max_cfl(i) = max_cfl(i) * flx_init_new(i)               &
     &                      / flx_init(i)
              else
                flx_init_new(i) = flx_init(i)
              End If  ! dcpbydt(i) > 0.0
            End Do ! j = 1, nterm

          End If ! cape_opt = 3, 4 or 5 (and test on w_max)

        Else If ( nterm > 0 ) then
!  This section as above except no test on w_max

          IF ( cape_opt  ==  3 ) THEN

!CDIR  NODEP
            Do j = 1, nterm
              i = index_nterm(j)
              If ( dcpbydt(i) > 0.0 ) then

                flx_init_new(i) = flx_init(i) * cape(i) /               &
     &                                    (cape_timescale * dcpbydt(i))

                If (flx_init_new(i) > flxmax_init(i)) then
                  flx_init_new(i) = flxmax_init(i)
                End If ! flx_init_new(i) > flxmax_init(i)
! Scale max_cfl with cape scale
                max_cfl(i) = max_cfl(i) * flx_init_new(i)               &
     &               / flx_init(i)
              else
                flx_init_new(i) = flx_init(i)
              End If  ! dcpbydt(i) > 0.0
            End Do ! j = 1, nterm

          ELSEIF ( cape_opt  ==  4 ) THEN

!CDIR  NODEP
            Do j = 1, nterm
              i = index_nterm(j)
              If ( dcpbydt(i) > 0.0 ) then

                cape_ts_new = CAPE_TIMESCALE * CAPE(I) / CAPE_min +            &
     &                       CAPE_TIMESCALE * EXP( - CAPE(I) /CAPE_min)

                flx_init_new(i) = flx_init(i) * cape(i) /               &
     &                                    ( cape_ts_new * dcpbydt(i))
                If (flx_init_new(i) > flxmax_init(i)) then
                  flx_init_new(i) = flxmax_init(i)
                End If ! flx_init_new(i) > flxmax_init(i)
! Scale max_cfl with cape scale
                max_cfl(i) = max_cfl(i) * flx_init_new(i)               &
     &                   / flx_init(i)
              else
                flx_init_new(i) = flx_init(i)
              End If  ! dcpbydt(i) > 0.0
            End Do ! j = 1, nterm

          ELSEIF (cape_opt  ==  5 ) THEN

!CDIR  NODEP
            Do j = 1, nterm
              i = index_nterm(j)
              If ( dcpbydt(i) > 0.0 ) then

                if ( relh(i) / dptot(i) >= 0.75 ) then
                  CAPE_TS_NEW = CAPE_TIMESCALE *                               &
     &                        ( 0.2373 / (relh(i) / dptot(i))**5)
                else
                  CAPE_TS_NEW = CAPE_TIMESCALE
                endif ! relh(i) / dptot(i) >= 0.75

                flx_init_new(i) = flx_init(i) * cape(i) /               &
     &                                    (cape_ts_new * dcpbydt(i))
                If (flx_init_new(i) > flxmax_init(i)) then
                  flx_init_new(i) = flxmax_init(i)
                End If ! flx_init_new(i) > flxmax_init(i)
! Scale max_cfl with cape scale
                max_cfl(i) = max_cfl(i) * flx_init_new(i)               &
     &                     / flx_init(i)
              else
                flx_init_new(i) = flx_init(i)
              End If  ! dcpbydt(i) > 0.0
            End Do ! j = 1, nterm

          End If ! cape_opt = 3,4 or 5 (no test on w_max)

        End If ! nterm > 0 .and. w_cape_limit < 1000.0

        If ( nterm > 0 ) then

!
! Work out scaled mass flux needed to keep cfl ratio below limit.
! L_CAPE assumed to be true.

          Do j = 1,nterm
            i = index_nterm(j)
            max_cfl(i) = max_cfl(i) * timestep

            If (max_cfl(i)  >   CFL_LIMIT) then
                flx_init_new(i) = flx_init_new(i) * CFL_LIMIT           &
     &                              / max_cfl(i)
            else
                flx_init_new(i) = flx_init_new(i)
            End If

            If (flx_init_new(i)  >   flxmax_init(i)) then
              flx_init_new(i) = flxmax_init(i)
            End If
            max_cfl(i) = 0.0
!
! Scale cloud fraction
!
            If (flx_init_new(i)  >   0.0) then
              scale_f(i) = flx_init_new(i) / flx_init(i)
              cca_2d(i) = cca_2d(i) + 0.06 * LOG(scale_f(i))

!
! Check scaled cloud fraction not smaller than minimum value
! (2.0E-5) or greater than unity.
!

              cca_2d(i) = max(2.0E-5,cca_2d(i))
              If (cca_2d(i)  >   1.0) then
                cca_2d(i) = 1.0
              End If
            End If

          End Do  ! nterm

!
! Carry out cape and cfl scaling
!

          Do kt = 2, k+1
! NEC compiler directive
!CDIR  NODEP
            Do j = 1,nterm
            i = index_nterm(j)
              If (kt  >=  start_lev(i) .and. flx_init_new(i)            &
     &                                            >   0.0) then

                If (kt  ==  det_lev(i)) then
                  dthbydt(i,kt)  = ((dthbydt(i,kt) - dthef(i))          &
     &                                       * scale_f(i)) + dthef(i)
                  dqbydt(i,kt)   = ((dqbydt(i,kt) - dqf(i))             &
     &                                       * scale_f(i)) + dqf(i)

                  If (L_q_interact) then  ! PC2
                    dqclbydt(i,kt)   = ((dqclbydt(i,kt) - dqclf(i))     &
     &                                       * scale_f(i)) + dqclf(i)
                    dqcfbydt(i,kt)   = ((dqcfbydt(i,kt) - dqcff(i))     &
     &                                       * scale_f(i)) + dqcff(i)
                    dcflbydt(i,kt)   = ((dcflbydt(i,kt) - dcflf(i))     &
     &                                       * scale_f(i)) + dcflf(i)
                    dcffbydt(i,kt)   = ((dcffbydt(i,kt) - dcfff(i))     &
     &                                       * scale_f(i)) + dcfff(i)
                    dbcfbydt(i,kt)   = ((dbcfbydt(i,kt) - dbcff(i))     &
     &                                       * scale_f(i)) + dbcff(i)
                  End If

                  If (l_mom_GK) then
                    dubydt(i,kt) = ((dubydt(i,kt) - duef(i))            &
     &                                       * scale_f(i)) + duef(i)
                    dvbydt(i,kt) = ((dvbydt(i,kt) - dvef(i))            &
     &                                       * scale_f(i)) + dvef(i)

                  End If

                  If (L_tracer) then
                    Do ktra = 1,ntra
                      dtrabydt(i,kt,ktra) = ((dtrabydt(i,kt,ktra) -     &
     &                                   dtraef(i,ktra)) * scale_f(i))  &
     &                                        + dtraef(i,ktra)
                    End Do
                  End If
                else
                  dthbydt(i,kt) = dthbydt(i,kt) * scale_f(i)
                  dqbydt(i,kt)  = dqbydt(i,kt) * scale_f(i)

                  If (L_q_interact) then  ! PC2
                    dqclbydt(i,kt)  = dqclbydt(i,kt) * scale_f(i)
                    dqcfbydt(i,kt)  = dqcfbydt(i,kt) * scale_f(i)
                    dcflbydt(i,kt)  = dcflbydt(i,kt) * scale_f(i)
                    dcffbydt(i,kt)  = dcffbydt(i,kt) * scale_f(i)
                    dbcfbydt(i,kt)  = dbcfbydt(i,kt) * scale_f(i)
                  End If

                  If (l_mom_GK) then
                    dubydt(i,kt)  = dubydt(i,kt) * scale_f(i)
                    dvbydt(i,kt)  = dvbydt(i,kt) * scale_f(i)
                    If (kt <  k+1) then
                      flxkp12(kt,i) = flxkp12(kt,i) * scale_f(i)
                    Endif
                  End If
                  If (L_tracer) then
                    Do ktra = 1,ntra
                    dtrabydt(i,kt,ktra) =dtrabydt(i,kt,ktra)*scale_f(i)
                    End Do
                  End If

                End If !kt=det_lev
                flx(i,kt)    = flx(i,kt) * scale_f(i)
                precip(i,kt) = precip(i,kt) * scale_f(i)

                If (flg_up_flx) then
                  up_flux(i,kt) = flx(i,kt)
                End If
                If (kt <  k+1) then
                  If (flg_up_flx_half) then
                    up_flux_half(i,kt+1) = flxkp12(kt,i)
                  End If
                End If

                entrain_up(i,kt) = entrain_up(i,kt) * scale_f(i)

                If (flg_detr_up) then
                  detrain_up(i,kt) = detrain_up(i,kt) * scale_f(i)
                End If

              End If ! kt>start_lev and flx_init_new >0
            End Do  ! j loop
          End Do  ! kt loop

!
! Set cape scaling parameters for next level
!

! NEC compiler directive
!CDIR  NODEP
          Do j = 1,nterm
            i = index_nterm(j)
            dthef(i) = dthbydt(i,k+1)
            dqf(i)   = dqbydt(i,k+1)
            If (L_q_interact) then  ! PC2
              dqclf(i)   = dqclbydt(i,k+1)
              dqcff(i)   = dqcfbydt(i,k+1)
              dcflf(i)   = dcflbydt(i,k+1)
              dcfff(i)   = dcffbydt(i,k+1)
              dbcff(i)   = dbcfbydt(i,k+1)
            End If
            det_lev(i)= k+1
            If (l_mom_GK) then
              duef(i)   = dubydt(i,k+1)
              dvef(i)   = dvbydt(i,k+1)
            End If
          End Do  ! nterm loop
          If (L_tracer) then
            Do ktra = 1,ntra
              Do j = 1,nterm
                i = index_nterm(j)
                dtraef(i,ktra) = dtrabydt(i,k+1,ktra)
              End Do  ! nterm loop
            End Do
          End If

!-----------------------------------------------------------------------
! 3.5  Downdraft calculation - on all points where convection is
!      terminating.
!
!      Subroutine DD_CALL
!
!      UM Documentation Paper 27, part 2
!
!-----------------------------------------------------------------------

        If (L_eman_dd) then

! save termination level

           Do j=1,nterm
             i = index_nterm(j)
             kterm_mid(i) = k
           End do

        Else        ! original down draughts

        npossdd = 0
        nnodd = 0
        Do j = 1,nterm
          i = index_nterm(j)
          tempnum = 0.0
          If (iccb(i)  >   0) then
            deltap_cld = p_layer_centres(i,iccb(i))                     &
     &                               - p_layer_centres(i,k)
            Do kt = iccb(i), k+1
              tempnum = tempnum + precip(i,kt)
            End Do
          else
            deltap_cld = 0.0
          End If

!
! Downdraughts possible if pressure thickness of convective
! cloud (deltap_cld) is greater than 15000m, the point is saturated
! and the precip. in the layer is greater than a threashold
! value (1E-12).
!

          If (deltap_cld  >   15000.0 .and. bgmk(i) .and.               &
     &                                     tempnum  >   1E-12) then
            npossdd = npossdd + 1
            index_possdd(npossdd) = i
          else
            nnodd = nnodd + 1
            index_nodd(nnodd) = i
          End If
        End Do  ! nterm loop

!
! If some downdraughts are possible (npossdd > 0) then call
! downdraught code
!

        If (npossdd  >   0) then
! DEPENDS ON: dd_call
          Call DD_CALL (npnts,npnts,k,thp(1,1),qp(1,1),th(1,1),         &
     &                  q(1,1),dthbydt(1,1),dqbydt(1,1),flx(1,1),       &
     &                  pstar,                                          &
     &                  p_layer_boundaries,                             &
     &                  p_layer_centres,                                &
     &                  exner_layer_boundaries,                         &
     &                  exner_layer_centres,                            &
     &                  precip(1,1), rain,snow, rain_3d, snow_3d,       &
     &                  iccb, icct, bwater(1,2),                        &
     &                  bterm,bgmk,timestep,cca_2d,npossdd,             &
     &                  L_tracer,ntra,trap,tracer,dtrabydt,nlev,trlev,  &
     &                  recip_pstar,dwn_flux,flg_dwn_flx,entrain_dwn,   &
     &                  flg_entr_dwn,detrain_dwn,flg_detr_dwn,          &
     &                  L_shallow,index_possdd)


        End If
!
! Surface precipitation calculation for points where downdraught not
! possible.
!

        If (nnodd  >   0) then
! DEPENDS ON: evap_bcb_nodd
          Call EVAP_BCB_NODD(npnts,npnts,nnodd,k,th,q,dthbydt,          &
     &                       dqbydt,exner_layer_centres,                &
     &                       exner_layer_boundaries,rain, snow,         &
     &                       rain_3d, snow_3d, precip,                  &
     &                       pstar,p_layer_centres,p_layer_boundaries,  &
     &                       bwater(1,2),cca_2d,iccb,timestep,          &
     &                       index_nodd)

        End If

        End if  ! Emanuel test
!
! Adjust cloud base, top and amount to prevent errors occurring in
! radiation scheme when iccb = icct (this happens when convection
! saturates upon forced detrainment).
!


        ! Prevent model from wiping out single layer cloud. Problems with the
        ! radiation scheme may well have occured if iccb == icct, however the
        ! radiation scheme expects the cloud base and top boundaries not
        ! layers as specified by iccb and icct in the convection scheme so
        ! that even with single layer clouds iccb and icct that are passed to
        ! the radiation scheme should not be equal to each other.

        If (.NOT. lcv_ccrad) Then
          Do j=1, nterm
            i = index_nterm(j)

            If (iccb(i) == icct(i)) Then
              iccb   (i) = 0
              icct   (i) = 0
              cca_2d (i) = 0.0
              tcw    (i) = 0.0
              cclwp  (i) = 0.0
            End If

            If (lcbase(i) == lctop(i)) Then
              lcbase (i) = 0
              lctop  (i) = 0
              lcca   (i) = 0.0
              lcclwp (i) = 0.0
            End If

          End Do
        End If      ! lcv_ccrad

!
! If convection has terminated write cape to diagnostic output
! variable (cape_out). Note if more than one lot of mid-level
! convection or if mid above deep or shallow then stores highest
! convection results only.
! Set kterm array which holds the level index for termination
! of convection.
! Zero integrals as convection terminated. Note may convect again at
! same locations at higher levels in the atmosphere.
!

! NEC compiler directive
!CDIR  NODEP
          Do j = 1,nterm
            i=index_nterm(j)
            cape_out(i) = cape(i)
            dcpbydt(i) = 0.0
            cape(i) = 0.0
            kterm(i) = k
            kterm_mid(i) = k
            bconv(i) = .false.
            bterm(i) = .false.     ! reset bterm to false
            start_lev(i) = 0
            relh(i)     = 0.0      ! needed for RH CAPE
            dptot(i)    = 0.0      !
          End Do

        End If        ! nterm >0

!-----------------------------------------------------------------------
! 3.6  End of main loop over levels
!-----------------------------------------------------------------------

      End Do

      If (lcv_ccrad) Then

        ! (Was copied out of scaling if test to ensure these limits at all
        ! times, not just when cca_2d is scaled)

        ! Check scaled cloud fraction not smaller than minimum value
        ! (2.0E-5) or greater than unity.

        Do i=1, npnts
          If (cca_2d(i) < 0.0) Then
            cca_2d(i) = MAX(2.0E-5, cca_2d(i))
            If (cca_2d(i) > 1.0) Then
              cca_2d(i) = 1.0
            End If
          End If
        End Do

      End If      ! lcv_ccrad

!-----------------------------------------------------------------------
! Emanuel down draughts - one call for all points with mid-level conv
!-----------------------------------------------------------------------

      If (L_eman_dd) then

! how many points have mid level convection?
! What is the maximum termination level ?

         kterm_max = 2
         nmid = 0
         Do i=1,npnts
           If (l_mid_all(i)) then
              nmid = nmid +1
              if (kterm_mid(i) >  kterm_max) then
                 kterm_max = kterm_mid(i)
              endif

           Endif
         Enddo

! Call routine to compress to just those points with convection
! and call Emanuel downdraughts then expand back to full grid.

! DEPENDS ON: eman_cex
         call eman_cex(npnts, nmid, kterm_max, nlev, trlev ,ntra        &
     &,                  kterm_mid, l_mid_all, l_tracer                 &
     &,                  exner_layer_centres,exner_layer_boundaries     &
     &,                  p_layer_centres, p_layer_boundaries            &
     &,                  timestep                                       &
     &,                  th, q, qse, tracer ,precip                     &
     &,                  dthbydt, dqbydt, dtrabydt                      &
     &,                  rain, snow, dwn_flux)

      End if

!-----------------------------------------------------------------------
! 4.0  Diffusive Convective momentum tranport option
!-----------------------------------------------------------------------

       If (l_mom .and. mid_cmt_opt == 1) then

! Maximum number of mid-level layers, required for holding cloud bases
! and top of each layer. Note each layer must be at least 2 levels and have at
! least one level between. The bottom model level is not used and don't expect
! convection in any stratospheric levels so trying estimating by dividing by 6.
! This should give a number bigger than required but not excessive.

         nmax_layers = nlev/6

         nmid = 0
         Do i=1,npnts
           If (l_mid_all(i)) then
             nmid = nmid +1
             index1(nmid) = i
           Endif
         Enddo
! Only call if there are points with mid_level convection.
         If (nmid > 0) then

! DEPENDS ON: mid_conv_dif_cmt
           call mid_conv_dif_cmt(npnts, nmid, nlev, nmax_layers,        &
                             index1, l_mid_all,                         &
                             timestep,                                  &
                             u, v, r_theta, r_rho,                      &
                             z_theta, z_rho, rho, rho_theta,            &
                             p_layer_boundaries,                        &
                             flx, entrain_up,                           &
                             dubydt, dvbydt, uw_mid, vw_mid )

         End If  ! test on mid points
       End If  ! test on mid_cmt_opt 1

!-----------------------------------------------------------------------
! 5.0  End Subroutine
!-----------------------------------------------------------------------

      Return
      END SUBROUTINE MID_CONV
#endif
