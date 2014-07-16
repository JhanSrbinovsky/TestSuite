#if defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Deep convection scheme
!


      SUBROUTINE DEEP_CONV(nbl,nlev,ntra,n_cca_lev,n_dp,trlev,          &
                             bland, delthvu, exner_layer_centres,       &
                             exner_layer_boundaries,flg_up_flx,         &
                             flg_entr_up, flg_detr_up, flg_dwn_flx,     &
                             flg_entr_dwn,flg_detr_dwn,flg_uw_deep,     &
                             flg_vw_deep, L_calc_dxek,                  &
                             L_q_interact,L_tracer,ntml,ntpar,          &
                             pstar,p_layer_centres,                     &
                             p_layer_boundaries,z_theta,z_rho,          &
                             r_theta,r_rho,rho_theta,rho,               &
                             r2rho_th,r2rho,dr_across_th,dr_across_rh,  &
                             q,q1_sd,t1_sd,th,                          &
                             timestep,u,v,uw0,vw0,w_max,wstar,zlcl_uv,  &
                             freeze_lev,recip_pstar,qse,                &
                             ad_on,mdet_on,ent_on,sdet_on,new_termc,    &
                             bulk_cf,cf_frozen,cf_liquid,qcf,           &
                             qcl,tracer,cape_out,cclwp,ccw,cca,         &
                             dbcfbydt,dcffbydt,dcflbydt,dqbydt,dqcfbydt,&
                             dqclbydt,dthbydt,                          &
                             dubydt,dvbydt,dtrabydt,                    &
                             detrain_up,detrain_dwn,                    &
                             entrain_up,entrain_dwn,                    &
                             iccb,icct,lcca,lcclwp,                     &
                             lcbase,lctop,rain,snow,                    &
                             rain_3d, snow_3d, up_flux,                 &
                             dwn_flux,uw_deep,vw_deep,kterm,tcw,cca_2d, &
                             rbuoy_p_out,the_out,thp_out,qe_out,qp_out  &
                             )

!
! Purpose:
!   Deep convection scheme - works on points diagnosed as deep in
!   subroutine CONV_DIAG.
!
!   Called by GLUE_CONV.
!
! Current owners of code: R A Stratton
!
! History:
! Version     Date     Comment
! -------     ----     -------
!
!   6.2      07/11/05  New code 5A scheme compatible with Emanuel
!                       downdraughts. R A Stratton
!            02/12/05  Added A Maidens adaptive changes.
!            08/12/05  Added D Wilson's chnages(adw2602) to 4A.
!   6.2      21/11/05  Pass through tice and qstice. Damian Wilson
!   6.2      25/12/05  Make CAPE_TIMESCALE function of length scale
!                      and/or w_max        Nigel Roberts/Anna Maidens
!   6.4      18/12/06  Add R_det as intent in. A. Maidens
!   6.4      04/12/06  Alter CMT to use kterm rather than ntpar. R A Stratton
!   6.4      05/02/07  Improve comments on if test on cape_opt 
!                      A. Maidens
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
      Use cv_cntl_mod, Only:                                            &
          lcv_3d_cca, lcv_ccrad

      Use cv_run_mod, Only:                                             &
          l_mom, l_eman_dd, cape_opt, cape_ts_w, cape_min,              &
          w_cape_limit, cape_timescale, deep_cmt_opt, bl_cnv_mix,       &
          cca2d_dp_opt, cca_dp_knob, ccw_dp_knob

      IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
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

      integer, intent(in) :: n_dp     ! No. of deep convection points

      integer, intent(in) :: trlev    ! No. of model levels on which
                                      ! tracers are included

      logical, intent(in) :: bland(n_dp) ! Land/sea mask

      real, intent(in)    :: delthvu(n_dp) ! a measure of CAPE used to cal wcld

      real, intent(in)    :: exner_layer_centres(n_dp,0:nlev) !Exner

      real, intent(in)    :: exner_layer_boundaries(n_dp,0:nlev)
                                      ! Exner at half level above
                                      ! exner_layer_centres

      logical, intent(in) :: flg_up_flx ! STASH flag for updraught
                                      ! mass flux

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

      logical, intent(in) :: flg_uw_deep ! STASH flag for x-comp deep
                                      ! stress

      logical, intent(in) :: flg_vw_deep ! STASH flag for y-comp deep

      logical, intent(in) :: L_calc_dxek ! Switch for calculation of
                                      ! condensate increment

      logical, intent(in) :: L_q_interact ! Switch allows overwriting
                                      ! parcel variables when
                                      ! calculating condensate incr.

      logical, intent(in) :: L_tracer ! Switch for inclusion of tracers

      integer, intent(in) :: ntml(n_dp) ! Top level of surface mixed
                                      ! layer defined relative to
                                      ! theta,q grid

      integer, intent(in) :: ntpar(n_dp) ! Top level of initial parcel
                                      ! ascent in BL scheme defined
                                      ! relative to theta,q grid

      real, intent(in)    :: pstar(n_dp) ! Surface pressure (Pa)

      real, intent(in)    :: p_layer_centres(n_dp,0:nlev) ! Pressure
                                      ! (Pa)


      real, intent(in)    :: p_layer_boundaries(n_dp,0:nlev) ! Pressure
                                      ! at half level above
                                      ! p_layer_centres (Pa)

      real, intent(in)    ::    &
        z_theta(n_dp,nlev)      & ! height of theta levels (m)
      , z_rho(n_dp,nlev)        & ! height of rho levels (m)
      , r_theta(n_dp,0:nlev)    & ! radius of theta levels (m)
      , r_rho(n_dp,nlev)        & ! radius of rho levels (m)
      , rho_theta(n_dp,nlev)    & ! density for theta lev (kg/m3)
      , rho(n_dp,nlev)          & ! density for rho lev (kg/m3)
      , r2rho_th(n_dp,nlev)     & ! radius**2 density for theta lev (kg/m)
      , r2rho(n_dp,nlev)        & ! radius**2 density for rho lev (kg/m)
      , dr_across_th(n_dp,nlev) & ! thickness of theta levels (m)
      , dr_across_rh(n_dp,nlev)   ! thickness of rho levels (m)

      real, intent(in)    :: q(n_dp,nlev) ! Model mixing ratio (kg/kg)

      real, intent(in)    :: q1_sd(n_dp) ! Standard deviation of
                                      ! turbulent flucts. of layer 1 q
                                      ! (kg/kg)

      real, intent(in)    :: t1_sd(n_dp) ! Standard deviation of
                                      ! turbulent flucts. of layer 1
                                      ! temp. (K)

      real, intent(in)    :: th(n_dp,nlev) !Model potential
                                      ! temperature (K)

      real, intent(in)    :: timestep ! Model timestep (s)

      real, intent(in)    :: u(n_dp,nlev) !Model U field (m/s)

      real, intent(in)    :: v(n_dp,nlev) !Model V field (m/s)

      real, intent(in)    :: uw0(n_dp) ! U-comp of surface stress
                                      ! (N/m2)

      real, intent(in)    :: vw0(n_dp) ! V-comp of surface stress
                                      ! (N/m2)

      real, intent(in)    :: W_MAX(n_dp) ! max w in column

      real, intent(in)    :: wstar(n_dp) ! Convective velocity scale
                                         ! (m/s)

      real, intent(in)    :: zlcl_uv(n_dp) !Lifting condensation level
                                      ! defined for the uv grid (m)

      integer, intent(in) :: freeze_lev(n_dp) ! Level index for
                                               ! freezing level

      real, intent(in) :: recip_pstar(n_dp) ! Reciprocal of pstar array

      real, intent(in) :: qse(n_dp,nlev) ! Saturation mixing ratio of
                                         ! cloud environment (kg/kg)

      integer, intent(in) :: ad_on      !flag for adaptive detrainment
                                        !0 = off, 1 = on

      integer, intent(in) :: mdet_on    !flag for adaptive mixing
                                        ! detrainment
                                        !0 = off, 1 = on

      integer, intent(in) :: ent_on     !flag for adaptive entraiment
                                        !0 = off, 1 = on

      integer, intent(in) :: sdet_on    !flag for smoothed forced detrainment
                                        !0 = off, 1 = on

      integer, intent(in) :: new_termc  !flag for simplified termination
                                        !0 = off, 1 = on

!
! Arguments with intent INOUT:
!

      real, intent(inout) :: bulk_cf(n_dp,nlev) ! Bulk total cloud
                                      ! volume ( )

      real, intent(inout) :: cf_frozen(n_dp,nlev) ! Frozen water cloud
                                      ! volume ( )

      real, intent(inout) :: cf_liquid(n_dp,nlev) ! Liq water cloud
                                      ! volume ( )

      real, intent(inout) :: qcf(n_dp,nlev) ! Ice condensate mix ratio
                                      ! (kg/kg)

      real, intent(inout) :: qcl(n_dp,nlev) ! Liq condensate mix ratio
                                      ! (kg/kg)

      real, intent(inout) :: tracer(n_dp,trlev,ntra) !Model tracer
                                      ! fields (kg/kg)

!
! Arguments with intent OUT:
!

      real, intent(out) :: cape_out(n_dp) ! Saved convective available
                                      ! potential energy for diagnostic
                                      ! output (J/kg)

      real, intent(out) :: cclwp(n_dp) ! Condensed water path (k/m2)

      real, intent(out) :: ccw(n_dp,nlev) ! Convective cloud liquid
                                       ! water on model levels (g/kg)

      real, intent(out) :: cca(n_dp,nlev) ! Convective cloud amount
                                       !  on model levels (fraction)

      real, intent(out) :: dbcfbydt(n_dp,nlev) ! Increments to
                                      ! total cld volume due to
                                      ! convection(/s)

      real, intent(out) :: dcffbydt(n_dp,nlev) ! Increments to ice
                                      ! cloud volume due to convection
                                      ! (/s)

      real, intent(out) :: dcflbydt(n_dp,nlev) ! Increments to liq
                                      ! cloud volume due to convection
                                      ! (/s)

      real, intent(out) :: dqbydt(n_dp,nlev) ! Increments to q due to
                                      ! convection (kg/kg/s)

      real, intent(out) :: dqcfbydt(n_dp,nlev) ! Increments to ice
                                      ! condensate due to convection
                                      ! (kg/kg/s)

      real, intent(out) :: dqclbydt(n_dp,nlev) ! Increments to liq
                                      ! condensate due to convection
                                      ! (kg/kg/s)

      real, intent(out) :: dthbydt(n_dp,nlev) ! Increments to potential
                                      ! temp. due to convection (K/s)

      real, intent(out) :: dubydt(n_dp,nlev+1) ! Increments to U due
                                      ! to CMT (m/s2)

      real, intent(out) :: dvbydt(n_dp,nlev+1) ! Increments to V due
                                      ! to CMT (m/s2)

      real, intent(out) :: dtrabydt(n_dp,nlev,ntra) !Increment to tracer
                                      ! due to convection (kg/kg/s)

      real, intent(out) :: detrain_up(n_dp,nlev) ! Fractional
                                      ! detrainment rate into updraughts
                                      ! (Pa/s)

      real, intent(out) :: detrain_dwn(n_dp,nlev) ! Fractional
                                      ! detrainment rate into
                                      ! downdraughts (Pa/s)


      real, intent(out) :: entrain_up(n_dp,nlev) ! Fractional
                                      ! entrainment rate into updraughts
                                      ! (Pa/s)

      real, intent(out) :: entrain_dwn(n_dp,nlev) ! Fractional
                                      ! entrainment rate into
                                      ! downdraughts (Pa/s)

      integer, intent(out) :: iccb(n_dp) ! Convective cloud base
                                      ! level (m)

      integer, intent(out) :: icct(n_dp) ! Convective cloud top
                                      ! level (m)

      real, intent(out) :: lcca(n_dp) ! Lowest conv. cloud amt. (%)

      real, intent(out) :: lcclwp(n_dp) ! Condensed water path for
                                      ! lowest conv. cld. level (kg/m2)

      integer, intent(out) :: lcbase(n_dp) ! Lowest conv. cloud base
                                      ! level (m)

      integer, intent(out) :: lctop(n_dp) ! Lowest conv. cloud top
                                      ! level (m)

      real, intent(out) :: rain(n_dp) ! Surface convective rainfall
                                      ! (kg/m2/s)

      real, intent(out) :: snow(n_dp) ! Surface convective snowfall
                                      ! (kg/m2/s)

      real, intent(out) :: rain_3d(n_dp,nlev) ! Convective rainfall flux
                                              ! (kg/m2/s)

      real, intent(out) :: snow_3d(n_dp,nlev) ! Convective snowfall flux
                                              ! (kg/m2/s)

      real, intent(out) :: up_flux(n_dp,nlev) ! Updraught mass flux
                                      ! (Pa/s)

      real, intent(out) :: dwn_flux(n_dp,nlev) ! Downdraught mass
                                      ! flux (Pa/s)

      real, intent(out) :: uw_deep(n_dp,nlev) ! X-comp. of stress
                                      ! from deep convection
                                      !(kg/m/s2)

      real, intent(out) :: vw_deep(n_dp,nlev) ! Y-comp. of stress
                                      ! from deep convection
                                      !(kg/m/s2)

      integer, intent(out) :: kterm(n_dp) ! Level at which deep
                                          ! convection terminates,
                                          ! required by mid level scheme
      real, intent(out) :: tcw(n_dp)   ! Total condensed water(kg/m2/s)
                                       ! required by mid-level CCA cal.

      real, intent(out) :: cca_2d(n_dp) ! 2D convective cloud amount (%)

!
! Adaptive detrainment output variables
!
      real, intent(out) :: rbuoy_p_out(n_dp,nlev)  !buoyancy excess

      real, intent(out) :: the_out(n_dp,nlev)    !th_E in parcel routine

      real, intent(out) :: thp_out(n_dp,nlev)    !th_P in parcel routine

      real, intent(out) :: qe_out(n_dp,nlev)     !q_E in parcel routine

      real, intent(out) :: qp_out(n_dp,nlev)     !q_P in parcel routine


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
!
! Adaptive detrainment output variables
!
      real :: rbuoy_p_here(n_dp)       !buoyancy excess

      real :: the_here(n_dp)         !th_E in parcel routine

      real :: thp_here(n_dp)         !th_P in parcel routine

      real :: qe_here(n_dp)          !q_E in parcel routine

      real :: qp_here(n_dp)          !q_P in parcel routine

      real :: rbuoy_p_old(n_dp)       !buoyancy excess on previous level

      real :: zk(n_dp)               !heights for use in calc
      real :: zkp12(n_dp)            !of moist static energy
      real :: zkp1(n_dp)


      integer :: index1(n_dp),index2(n_dp)

      integer :: ncposs               ! No. of points which may convect

      integer :: nconv                ! No. of convecting points

      real :: amdetk(n_dp)            ! Mixing detrainment coefficient
                                      ! at level k multiplied by
                                      ! appropriate layer thickness

      real :: b_calc                  ! Coefficient in thpert calc.

      real :: c_calc                  ! Coefficient in thpert calc.

      real :: cape(n_dp)              ! Convective available potential
                                      ! energy (J/kg)

      real :: dq_sat_env              ! dqsat/dT  (kg/kg/K)

      real :: dcpbydt(n_dp)           ! Rate of change of cape (J/kg/s)

      real :: depth(n_dp)             ! Depth of convective cloud (m)

      real :: delexkp1(n_dp)          ! Difference in exner ratio
                                      ! across layer k+1

      real :: dqsthk(n_dp)            ! Gradient of saturation mixing
                                      ! ratio of cloud environment with
                                      ! theta in layer k (kg/kg/K)

      real :: dqsthkp1(n_dp)          ! Gradient of saturation mixing
                                      ! ratio of cloud environment with
                                      ! theta in layer k+1 (kg/kg/K)

      real :: ekp14(n_dp)             ! Entrainment coefficients at
                                      ! level k+1/4 multiplied by
                                      ! appropriate layer thickness
                                      !(dimensionless)

      real :: ekp34(n_dp)             ! Entrainment coefficients at
                                      ! level k+3/4 multiplied by
                                      ! appropriate layer thickness
                                      !(dimensionless)

      real :: ekm14(n_dp)             ! Entrainment coefficients at
                                      ! level k-1+1/4 multiplied by
                                      ! appropriate layer thickness
                                      !(dimensionless)
      real :: exk(n_dp)               ! Exner ratio at layer k

      real :: exkp1(n_dp)             ! Exner ratio at layer k+1

      real :: flxmax(n_dp)            ! Maximum initial convective
                                      ! mass flux (Pa/s)

      real :: flx_init(n_dp)          ! Initial mass flux at cloud base
                                      ! (Pa/s)

      real :: flx_init_new(n_dp)      ! flx_init scaled to destroy cape
                                      ! over timescale cape_timescale (Pa/s)

      real :: flxmax_init(n_dp)       ! Maximum possible initial mass
                                      ! flux (limited to the mass in
                                      ! the initial convecting layer
                                      ! in Pa/s)

      real :: max_cfl(n_dp)           ! Max cfl ratio over a convecting
                                      ! layer

      real :: p_lcl(n_dp)             ! Pressure at LCL (Pa)

      real :: precip(n_dp,nlev)       ! Amount of precip from each layer
                                      ! from each layer (kg/m/s)

!      real :: pt(n_dp)                ! Temporary store for P in calc.
                                      ! of sat. mixing ratio. (Pa)

      real :: pk(n_dp)                ! Pressure at midpoint of layer
                                      ! k (Pa)

      real :: pkp1(n_dp)              ! Pressure at midpoint of layer
                                      ! k+1 (Pa)

      real :: delpk(n_dp)             ! Pressure difference over layer
                                      ! k (Pa)

      real :: delpkp1(n_dp)           ! Pressure difference over layer
                                      ! k+1 (Pa)

      real :: delpkp12(n_dp)          ! Pressure difference between
                                      ! layers k and k+1 (Pa)

      real :: delp_uv_k(n_dp)         ! Pressure difference across uv
                                      ! layer k (Pa)

      real :: delp_uv_kp1(n_dp)       ! Pressure difference across uv
                                      ! layer k+1 (Pa)

      real :: q_lcl(n_dp)             ! Mixing ratio at LCL (kg/kg)

      real :: qse_lcl(n_dp)           ! Saturated q at LCL (kg/kg)

      real :: rhum(n_dp)              ! Dummy relative humidity
                                      ! (only used on shallow points)

      real :: t_lcl(n_dp)             ! Temperature at LCL (K)

      real :: th_lcl(n_dp)            ! Theta at LCL (K)

      real :: thv_pert                ! Theta_v parcel pertubation (K)

      real :: thpert(n_dp)            ! Theta parcel pertubation (K)
      
      real :: qpert(n_dp)             ! q parcel pertubation (kg/kg)

!      real :: tt(n_dp)                ! Temporary store for T in calc.
!                                      ! of saturation mixing ratio. (K)

!      real :: ttkm1(n_dp)             ! Temporary store for T in layer
                                      ! k-1 for use in freezing level
                                      ! calc. for anvil. (K)

      integer :: start_lev3c(n_dp)    ! PC2 Compressed convection
                                      ! initiation level

      logical :: L_shallow(n_dp)      ! Dummy variable (=.F.)

      logical :: L_mid(n_dp)          ! Dummy variable (=.F.)

      logical :: cumulus(n_dp)        ! Dummy variable (=.T.)

      logical :: bgmk(n_dp)           ! Mask for points where parcel in
                                      ! layer k is saturated

      logical :: bwater(n_dp,2:nlev)  ! Mask for points at which
                                      ! condensate is liquid

      logical :: bwk(n_dp)            !mask for liquid condensate on k

      logical :: bwkp1(n_dp)          !mask for liquid condensate on k+1

      logical :: blowst(n_dp)         ! Dummy variable indicating low
                                      ! enough stability for convection
                                      ! to occur

      logical :: bterm(n_dp)          ! Mask for points which have
                                      ! stopped convecting

      logical :: bconv(n_dp)          ! Mask for points at which
                                      ! convection is occurring

      logical :: bcposs(n_dp)         ! Mask for points passing
                                      ! initial stability test

!
! Parcel variables
!

      real :: qpi(n_dp)               ! Initial parcel mixing ratio
                                      !(kg/kg)

      real :: qp(n_dp,nlev)           ! Parcel mixing ratio (kg/kg)

      real :: thpi(n_dp)              ! Initial parcel potential temp.
                                      !(K)

      real :: thp(n_dp,nlev)          ! Parcel potential temp (K)

      real :: trap(n_dp,nlev,ntra)    ! Tracer content of parcel
                                      ! (kg/kg)
      real :: expi(n_dp)
                                      ! Initial parcel exner pressure

      real :: xpk(n_dp,nlev)          ! Parcel cloud water (kg/kg)

      real :: flx(n_dp,nlev)          ! Parcel massflux (Pa/s)

      real :: xsbmin_v(n_dp,nlev)     ! Minmum parcel buoyancy excess

      real :: thpixs_v(n_dp,nlev)     ! Theta parcel excess (K)

      real :: qpixs_v(n_dp,nlev)      ! Q parcel excess(kg/kg)

!
! PC2
!

      real :: qclp(n_dp,nlev)         ! Parcel liquid condensated mixing
                                      ! ratio in layer k (kg/kg)

      real :: qcfp(n_dp,nlev)         ! Parcel frozen condensated mixing
                                      ! ratio in layer k (kg/kg)

!
! Parameters
!

      real, parameter :: CFL_LIMIT = 1.0 ! Max CFL ratio allowed

!
! CMT variables  - those used depend on scheme 
!
! Required by Gregory-Kershaw scheme operating in plume calculation

      Integer ::         &
       nstart(n_dp)        ! Level for start of plume

      Real ::            &
       eflux_u_ud(n_dp)  & ! Vertical eddy flux of momentum due to UD at
                           !  top of layer (Pa m/s)
      ,eflux_v_ud(n_dp)  & ! Vertical eddy flux of momentum due to UD at
                           ! bottom of layer (Pa m/s)
      ,up(n_dp,nlev)     & ! Parcel U (m/s)
      ,vp(n_dp,nlev)     & ! Parcel V (m/s)
      ,zsurf(n_dp)         ! Height of start of plume = 0.1*zlcl

! Required by Turbulence base scheme called after plume calculation

      Integer ::             &
       nlcl_uv(n_dp+1)       & ! Level index for LCL
      ,ntop_uv(n_dp+1)       & ! Level index for top of layer
      ,n_0degc(n_dp+1)       & ! Level index for zero degrees
      ,cu_term(n_dp)         & ! Indicies for CMT subs
      ,cu_tend(n_dp)           ! Indicies for CMT subs

      Real ::                &  
       mass_dwn(nlev,n_dp+1) & ! Downdraught mass flux (Pa/s)
      ,p_uv(nlev,n_dp+1)     & ! Pressure of model level (Pa)
      ,phalf_uv(nlev,n_dp+1) & ! Pressure of half level (Pa)
      ,plcl_uv(n_dp+1)       & ! Pressure at LCL (Pa)
      ,ptop_uv(n_dp+1)       & ! Pressure at top of cloud layer (Pa)
      ,p_0degc_uv(n_dp+1)    & ! Pressure of zero degree level (Pa)
      ,rho_uv(nlev,n_dp+1)   & ! Density on uv level (kg/m3)
      ,visc(nlev,n_dp+1)     & ! CMT eddy viscosity (m2/s)
      ,uw(nlev,n_dp+1)       & ! U- comp stress profile (N/m2)
                               ! (units vary through calls)
      ,vw(nlev,n_dp+1)       & ! V-comp stress profile (N/m2)
      ,uw_base(nlev,n_dp+1)  & ! Cloud base U stress (N/m2)
      ,vw_base(nlev,n_dp+1)  & ! Cloud base V stress (N/m2)
      ,ue_p(nlev,n_dp+1)     & ! Environment U profile (m/s)
      ,ve_p(nlev,n_dp+1)       ! Environment V profile (m/s)

      real :: exk_temp                ! Temporary exner

! Required by all version of CMT

      real :: flxkp12(nlev,n_dp+1)    ! Mass flux on half level (Pa/s)

      real :: mb(n_dp+1)              ! Cloud base mass flux (Pa/s)

      Logical :: l_mom_gk             ! true if Gregory-Kershaw CMT required

!
! Cape scaling/closure variables
!

      integer :: det_lev(n_dp)        ! Level at which split final
                                      ! detrainment last occurred

      integer :: nterm                ! No. of points where conv.
                                      ! has terminated

      integer :: index_nterm(n_dp)    ! Index for points where conv.
                                      ! has terminated

      real :: tempnum                 ! Temporary variable for storage

      real :: scale_f(n_dp)           ! scale factor

      real :: cape_ts_new             ! Used as variable in RH-based
                                      ! closure

      real :: relh(n_dp)              ! RH integral (average when
                                      ! convection terminates)

      real :: dptot(n_dp)             ! Delta P integral

!
! Downdraught scheme variables
!

      integer :: npossdd              ! Max. no. of downdraughts
                                      ! possible

      integer :: nnodd                ! No. of downdraughts not possible

      integer :: index_possdd(n_dp)   ! Index of downdraughts possible

      integer :: index_nodd(n_dp)     ! Index of downdraughts not
                                      ! possible
      integer :: kmax_term            ! maximum termination level + 1

      real :: deltap_cld              ! Pressure thickness of convective
                                      ! cloud (Pa)

!
! Local compressed arrays
!

      logical :: bconv_c2(n_dp)

      logical :: bgmkp1_c(n_dp), bgmkp1_c2(n_dp) ! Mask for points
                                      ! where parcel in layer k+1
                                      ! is saturated

      logical :: bwk_c(n_dp), bwk_c2(n_dp) ! bwater mask in layer k

      logical :: bwkp1_c(n_dp), bwkp1_c2(n_dp) ! bwater mask in layer
                                      ! k+1

      real :: deltak_c2(n_dp)         ! Parcel forced detrainment rate
                                      ! in layer k multiplied by
                                      ! appropriate layer thickness

      real :: dqek_c2(n_dp)           ! Increment to q due to
                                      ! convection in layer k (kg/kg)

      real :: dqekp1_c2(n_dp)         ! Increment to q due to
                                      ! convection in layer k+1 (kg/kg)

      real :: dthek_c2(n_dp)          ! Increment to potential temp.
                                      ! due to convection in layer k

      real :: dthekp1_c2(n_dp)        ! Increment to potential temp.
                                      ! due to convection in layer k+1

      real :: dtraek_c2(n_dp,ntra)    ! Increment to model tracer due
                                      ! to conv. at level k (kg/kg/s)

      real :: dtraekp1_c2(n_dp,ntra)  ! Increment to model tracer due
                                      ! to conv. at level k+1 (kg/kg/s)

      real :: duek_c2(n_dp)           ! Increment to model U in layer k
                                      ! due to CMT (m/s2)

      real :: duekp1_c2(n_dp)         ! Increment to model U in layer
                                      ! k+1 due to CMT (m/s2)

      real :: dvek_c2(n_dp)           ! Increment to model V in layer k

      real :: dvekp1_c2(n_dp)         ! Increment to model V in layer
                                      ! k+1 due to CMT (m/s2)

      real :: flxk_c(n_dp), flxk_c2(n_dp) !Parcel mass flux in layer k
                                      ! (Pa/s)

      real :: flxkp12_c2(n_dp)        ! Half level mass flux (Pa/s)

      real :: prekp1_c2(n_dp)         ! Precip. from parcel as it rises
                                      ! from layer k to k+1 (kg/m2/s)

      real :: qpk_c(n_dp), qpk_c2(n_dp) ! Parcel mixing ratio in
                                      ! layer k(kg/kg)

      real :: qpk(n_dp)
      real :: qpkp1_c(n_dp), qpkp1_c2(n_dp) ! Parcel mixing ratio
                                      ! in layer k+1 (kg/kg)

      real :: qek_c(n_dp), qek_c2(n_dp) ! Env. mixing ratio in
                                      ! layer k (kg/kg)

      real :: qek(n_dp)
      real :: qekp1_c(n_dp), qekp1_c2(n_dp) ! Env. mixing ratio in
                                      ! layer k+1 (kgkg-1)

      real :: qekp1(n_dp)
      real :: qsek_c2(n_dp)           ! Saturation mixing ratio of
                                      ! cld. env. in layer k (kg/kg)

      real :: qsek(n_dp)
      real :: qsekp1_c(n_dp), qsekp1_c2(n_dp) ! Saturation mixing ratio
                                      ! of cld. env. in layer k+1
                                      ! (kg/kg)

      real :: qsekp1(n_dp)
      real :: thek_c(n_dp), thek_c2(n_dp) ! Env. potential temp
                                      ! in layer k (K)

      real :: thek(n_dp)
      real :: thekp1_c(n_dp), thekp1_c2(n_dp) ! Env. potential temp i
                                      ! in layer k (K)

      real :: thekp1(n_dp)
      real :: thpk_c(n_dp), thpk_c2(n_dp) ! Parcel potential temp
                                      ! in layer k (K)

      real :: thpk(n_dp)
      real :: thpkp1_c(n_dp), thpkp1_c2(n_dp)! Parcel potential temp
                                      ! in layer k (K)

      real :: thpkp1(n_dp)
      real :: traek_c(n_dp,ntra), traek_c2(n_dp,ntra) ! Tracer content
                                      ! cld. env. in layer k (kgkg-1)

      real :: traekp1_c(n_dp,ntra), traekp1_c2(n_dp,ntra) ! Tracer
                                      ! content of cloud env.
                                      ! in layer k+1 (kg/kg)

      real :: trapk_c(n_dp,ntra), trapk_c2(n_dp,ntra) ! Tracer cont.
                                      ! of parcel in layer k (kg/kg)

      real :: trapkp1_c(n_dp,ntra), trapkp1_c2(n_dp,ntra) ! Tracer cont.
                                      ! of parcel in layer k+1 (kg/kg)

      real :: rbuoy_c(n_dp), rbuoy_c2(n_dp) ! Buoyancy of parcel at k+1
                                      ! (Kelvin)

      real :: uek_c(n_dp), uek_c2(n_dp) ! Model U field on layer k
                                      ! (m/s)

      real :: uekp1_c(n_dp), uekp1_c2(n_dp)! Model U field on layer
                                      ! k+1 (m/s)

      real :: vek_c(n_dp), vek_c2(n_dp) ! Model V field on layer k
                                      ! (m/s)

      real :: vekp1_c(n_dp), vekp1_c2(n_dp) ! Model V field on layer
                                      ! k+1 (m/s)

      real :: upk_c(n_dp), upk_c2(n_dp) ! Parcel U in layer k
                                      ! after entrainment (m/s)


      real :: upkp1_c(n_dp), upkp1_c2(n_dp) ! Parcel U in layer k+1
                                      ! after entrainment (m/s)

      real :: vpk_c(n_dp), vpk_c2(n_dp) ! Parcel V in layer k
                                      ! after entrainment (m/s)

      real :: vpkp1_c(n_dp), vpkp1_c2(n_dp) ! Parcel V in layer k+1
                                      ! after entrainment (m/s)

      real :: xsqkp1_c(n_dp), xsqkp1_c2(n_dp) ! Excess water vapour
                                      ! in parcel at k+1 (kg/kg)

! PC2 compression arrays

      real :: qclek_c(n_dp), qclek_c2(n_dp) ! Environment liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qclekp1_c(n_dp), qclekp1_c2(n_dp) ! Environment liquid
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qcfek_c(n_dp), qcfek_c2(n_dp) ! Environment frozen
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qcfekp1_c(n_dp), qcfekp1_c2(n_dp) ! Environment frozen
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qclpk_c(n_dp), qclpk_c2(n_dp) ! Parcel liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qclpkp1_c(n_dp), qclpkp1_c2(n_dp) ! Parcel liquid
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qcfpk_c(n_dp), qcfpk_c2(n_dp) ! Parcel frozen
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qcfpkp1_c(n_dp), qcfpkp1_c2(n_dp) ! Parcel frozen
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: cflek_c2(n_dp),cflekp1_c2(n_dp) ! Environment liquid water
                                      ! cloud volume ( )

      real :: cffek_c2(n_dp),cffekp1_c2(n_dp) ! Environment frozen water
                                      ! cloud volume ( )

      real :: bcfek_c2(n_dp),bcfekp1_c2(n_dp) ! Environment bulk total
                                      ! cloud volume ( )

      real :: dqclek_c2(n_dp),dqclekp1_c2(n_dp) ! Environment increments
                                      ! to liquid condensate mixing
                                      ! ratio to convection (kg/kg/s)

      real :: dqcfek_c2(n_dp),dqcfekp1_c2(n_dp) ! Environment increments
                                      ! to frozen condensate mixing
                                      ! ratio to convection (kg/kg/s)

      real :: dcflek_c2(n_dp),dcflekp1_c2(n_dp) ! Environment increments
                                      ! to liquid water cloud volume due
                                      ! to convection (/s)

      real :: dcffek_c2(n_dp),dcffekp1_c2(n_dp) ! Environment increments
                                      ! to frozen water cloud volume due
                                      ! to convection (/s)

      real :: dbcfek_c2(n_dp),dbcfekp1_c2(n_dp) ! Environment increments
                                      ! to bulk total cloud volume due
                                      ! to convection (/s)

      real :: amdetk_c2(n_dp)
      logical :: bgmk_c2(n_dp)
      logical :: bland_c2(n_dp)
      logical :: blowst_c2(n_dp)
      logical :: bterm_c2(n_dp)
      real :: cape_c2(n_dp)
      real :: cca_2d_c2(n_dp)
      real :: cclwp_c2(n_dp)
      real :: ccw_c2(n_dp)
      logical :: cumulus_c(n_dp), cumulus_c2(n_dp)
      real :: dcpbydt_c2(n_dp)
      real :: delexkp1_c2(n_dp)
      real :: delpk_c2(n_dp)
      real :: delpkp1_c2(n_dp)
      real :: delp_uv_k_c2(n_dp)
      real :: delp_uv_kp1_c2(n_dp)
      real :: depth_c2(n_dp)
      real :: dptot_c2(n_dp)
      real :: dqsthkp1_c2(n_dp)
      real :: dqsthk_c2(n_dp)
      real :: eflux_u_ud_c2(n_dp)
      real :: eflux_v_ud_c2(n_dp)
      real :: ekp14_c(n_dp),ekp14_c2(n_dp)
      real :: ekp34_c(n_dp),ekp34_c2(n_dp)
      real :: exk_c2(n_dp)
      real :: exkp1_c(n_dp),exkp1_c2(n_dp)
      real :: expi_c2(n_dp)
      integer :: icct_c2(n_dp)
      integer :: iccb_c2(n_dp)
      real :: lcclwp_c2(n_dp)
      integer :: lctop_c2(n_dp)
      integer :: lcbase_c2(n_dp)
      real :: lcca_c2(n_dp)
      logical :: L_shallow_c2(n_dp)
      logical :: L_mid_c2(n_dp)
      real :: max_cfl_c2(n_dp)
      real :: pk_c(n_dp),pk_c2(n_dp)
      real :: pkp1_c(n_dp),pkp1_c2(n_dp)
      real :: pstar_c2(n_dp)
      real :: q1_sd_c2(n_dp)
      real :: qpi_c2(n_dp)
      real :: qpixs_v_c2(n_dp)
      real :: relh_c2(n_dp)
      real :: rbuoy_p_here_c2(n_dp)
      real :: the_here_c2(n_dp)
      real :: thp_here_c2(n_dp)
      real :: qe_here_c2(n_dp)
      real :: qp_here_c2(n_dp)
      real :: rbuoy_p_old_c2(n_dp)
      real :: tcw_c2(n_dp)
      real :: thpi_c2(n_dp)
      real :: thpixs_v_c2(n_dp)
      real :: t1_sd_c2(n_dp)
      real :: xpk_c(n_dp),xpk_c2(n_dp)
      real :: xsbmin_v_c2(n_dp)
      logical :: b_nodd(n_dp)   ! points with no downdraught
      logical :: b_dd(n_dp)     ! points with downdraught on termination
!
! required by water conservation check
!
      real ::                                                           &
     & qMinInColumn(n_dp)                                               &
                              ! Minimum value for q in column(kg/kg)
     &,temp1(n_dp)            ! work array

! required by CMT

      real ::   &
        wcld(n_dp)    &  ! Convective veloicty scale
       ,zlcl(n_dp)       ! lifting condensation level

      real, parameter :: QMIN = 1.0E-8 ! Global minimum allowed Q

      !===============================================================
      ! CCRad Variables local variables
      !===============================================================

      real             :: a_land = 0.3   ! parameters used to relate
      real             :: a_sea  = 0.3   ! cca_2d of cld to its
      real             :: b_land = 0.025 ! precipitation rate
      real             :: b_sea  = 0.025 !

      INTEGER, Parameter   :: cca2d_total_condensed_water = 0
      INTEGER, Parameter   :: cca2d_grant_lock            = 1 ! Shallow cnv
      INTEGER, Parameter   :: cca2d_srf_precip_rate       = 1 ! mid/deep cnv


      !===============================================================
      ! End CCRad Variables local variables
      !===============================================================
!
! Loop counters
!

      integer :: i,j,k,ktra,kt

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
#include "wcape_spd_up.h"
!
!initialise SCM diagnostics
!
      Do k = 1,nlev
        Do i = 1,n_dp
          rbuoy_p_out(i,k)=0.0
          the_out(i,k)=th(i,k)
          thp_out(i,k)=th(i,k)
          qe_out(i,k)=q(i,k)
          qp_out(i,k)=q(i,k)
        End Do
      End Do
!
! Initialise logicals
!

      Do i = 1,n_dp
        blowst(i)    = .true.
        bterm(i)     = .false.
        bconv(i)     = .false.
        bcposs(i)    = .false.
        cumulus(i)   = .true.
        L_shallow(i) = .false.
        L_mid(i)     = .false.
        b_nodd(i)    = .false.
        b_dd(i)      = .false.
      End Do

!-----------------------------------------------------------------------
! 1.0  Create saturation mixing ratio arrays (now in glue)
!-----------------------------------------------------------------------
!
! Re-calculate XSBMIN and THPIXS constants based on layer thickness (Pa)
!
      Do k = 1,nlev-1
          Do i = 1,n_dp
            xsbmin_v(i,k) = min( ((p_layer_centres(i,k) -               &
     &              p_layer_centres(i,k+1))/5000.),1.0) *0.2

            thpixs_v(i,k) = min( ((p_layer_centres(i,k) -               &
     &              p_layer_centres(i,k+1))/5000.),1.0) * THPIXS_DEEP

            qpixs_v(i,k)  = QPIXS_DEEP
          End Do
      End Do  ! nlev


!
! Calculate cloud base mass flux
!

        Do i = 1,n_dp
          mb(i) = C_MASS * wstar(i)
        End Do
!
! Define the LCL at the half level above ntml. Find environmental
! T at p_lcl by approximating theta there with
! th(i,k) + constant*(th(i,k+1)-th(i,k))  where constant is tunable.
! Similarly for q.
!

        Do i = 1,n_dp
          k=ntml(i)
          p_lcl(i)  = p_layer_boundaries(i,k)
          th_lcl(i) = th(i,k) + 0.1 * (th(i,k+1) - th(i,k))
          t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / PREF)**KAPPA)
          q_lcl(i)  = q(i,k) + 0.1 * (q(i,k+1) - q(i,k))
        End Do
!
! Calculate saturation mixing ratio at LCL
!

! DEPENDS ON: qsat_mix
      Call QSAT_mix(qse_lcl,t_lcl,p_lcl,n_dp,.false.)


!-----------------------------------------------------------------------
! Initialize arrays required for Convective Momentum Transport(CMT)
!-----------------------------------------------------------------------
      If (l_mom) Then
        Select Case (deep_cmt_opt)
        
        Case(2)         ! Gregory-Kershaw CMT

          ! need level near surface for initial parcel U & V values
          ! zsurf = 0.1*z_lcl

          Do i = 1,n_dp
            zsurf(i)  = 0.1*z_rho(i,ntml(i))
          End Do
          Do k=nlev-1,1,-1
            Do i = 1,n_dp
             If (zsurf(i) <= z_theta(i,k)) Then
               nstart(i) = k
             End If
            End Do
          End Do
          l_mom_gk = .true.
        
        Case Default    ! (0/1) Alan Grant's eddy viscosity based CMT 
!
! Note: In terms of array indices p and phalf follow the convention
!       used in the boundary layer scheme. phalf(k,*) refers to the
!       lower boundary of uv layer k. This follows the convention for
!       um UM4.5 and before
!
!       Also note that p_layer_boundaries(0) and p_layer_centres(0)
!       = pstar, so p_uv(k,1) and phalf_uv(k,1) will be equal.
!
!       Because of the definition of nlcl, the pressure of the top of
!       the mixed layer is phalf_uv(nlcl,*)

 
          k=1
          Do i = 1,n_dp
            p_uv(k,i)     = p_layer_boundaries(i,k-1)
            phalf_uv(k,i) = p_layer_centres(i,k-1)
            ue_p(k,i)     = u(i,k)
            ve_p(k,i)     = v(i,k)
            nlcl_uv(i)    = ntml(i) + 1
            n_0degc(i)    = freeze_lev(i)
            kterm(i)=0
          End Do

          Do i = 1,n_dp
            Do k = 2,nlev
              p_uv(k,i)     = p_layer_boundaries(i,k-1)
              phalf_uv(k,i) = p_layer_centres(i,k-1)
              ue_p(k,i)     = u(i,k)
              ve_p(k,i)     = v(i,k)
              exk_temp      = (p_uv(k,i)/pref)**kappa
              rho_uv(k,i)   = 2.0 * p_uv(k,i) / (R * exk_temp *        &
                              (th(i,k-1) + th(i,k)))
            End Do
            plcl_uv(i)      = phalf_uv(nlcl_uv(i),i)
            p_0degc_uv(i)   = phalf_uv(n_0degc(i),i)
            rho_uv(1,i)     = rho_uv(2,i)
          End Do

          l_mom_gk = .false.
 
        End Select      ! test on deep_cmt_opt

      Else

! Initialise variable
        l_mom_gk = .false.

      End If     !L_mom

!-----------------------------------------------------------------------
! Calculate theta and q pertubation (pertubation is based on
! environment buoyancy gradient)
! Re-set th and q xs's at ntml
!

      Do i = 1,n_dp
        If (t_lcl(i) >  TM) then
          dq_sat_env = EPSILON * LC * qse_lcl(i)                        &
     &                  / (R * t_lcl(i) * t_lcl(i))
        else
          dq_sat_env = EPSILON * (LC+LF) * qse_lcl(i)                   &
     &                  / (R * t_lcl(i) * t_lcl(i))
        End If

        b_calc   = t_lcl(i) * C_VIRTUAL * dq_sat_env + 1.0              &
     &             + C_VIRTUAL * qse_lcl(i)

        thv_pert = -0.5 * (th(i,ntml(i)+1)                              &
     &                   * (1.0+C_VIRTUAL * q(i,ntml(i)+1)) -           &
     &                 th(i,ntml(i)) * (1.0 + C_VIRTUAL                 &
     &                   * q(i,ntml(i)))) + 0.5

        c_calc   = th_lcl(i) * C_VIRTUAL * (qse_lcl(i)                  &
     &                   - q_lcl(i)) - thv_pert

        thpert(i) = -c_calc / b_calc   ! ignore term in thpert**2

        thpixs_v(i,ntml(i)) = thpert(i)

        qpert(i)  = qse_lcl(i) + ((p_lcl(i) / PREF)                     &
     &                         **KAPPA) * thpert(i) * dq_sat_env        &
     &                           - q_lcl(i)

        qpixs_v(i,ntml(i))  = qpert(i)

      End Do ! n_dp


!
! Set bwater=.true. on points where water will condense rather than
! ice.
! SUBROUTINE FLAG_WET
! UM Documentation paper 27, section (2B)
!

! DEPENDS ON: flag_wet
      Call FLAG_WET(bwater,th,exner_layer_centres,n_dp,n_dp,nlev)

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 2.1  Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
!      increment arrays and cca
!      and dqcl/dt, dqcf/dt, dcfl/dt, dcff/dt, dbcf/dt
!-----------------------------------------------------------------------

      Do k = 1,nlev
        Do i = 1,n_dp
          precip(i,k) = 0.0
          ccw(i,k)    = 0.0
          dthbydt(i,k) = 0.0
          dqbydt(i,k)  = 0.0
          dqclbydt(i,k) = 0.0
          dqcfbydt(i,k) = 0.0
          dbcfbydt(i,k) = 0.0
          dcflbydt(i,k) = 0.0
          dcffbydt(i,k) = 0.0
        End Do
      End Do

      If (L_mom) then
        Do k = 1,nlev+1
          Do i = 1,n_dp
            dubydt(i,k) = 0.0
            dvbydt(i,k) = 0.0
          End Do
        End Do
      End If  ! L_mom

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,nlev
            Do i = 1,n_dp
              dtrabydt(i,k,ktra) = 0.0
            End Do
          End Do
        End Do
      End If  ! L_tracer

!-----------------------------------------------------------------------
! 2.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------

      If (flg_up_flx) then
        Do k = 1,nlev
          Do i = 1,n_dp
            up_flux(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_dwn_flx) then
        Do k = 1,nlev
          Do i = 1,n_dp
            dwn_flux(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_entr_up) then
        Do k = 1,nlev
          Do i = 1,n_dp
            entrain_up(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_detr_up) then
        Do k = 1,nlev
          Do i = 1,n_dp
            detrain_up(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_entr_dwn) then
        Do k = 1,nlev
          Do i = 1,n_dp
            entrain_dwn(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_detr_dwn) then
        Do k = 1,nlev
          Do i = 1,n_dp
            detrain_dwn(i,k) = 0.0
          End Do
        End Do
      End If

      If (L_mom) then
        If (flg_uw_deep) then
          Do k = 1,nlev
            Do i = 1,n_dp
              uw_deep(i,k) = 0.0
            End Do
          End Do
        End If
        If (flg_vw_deep) then
          Do k = 1,nlev
            Do i = 1,n_dp
              vw_deep(i,k) = 0.0
            End Do
          End Do
        End If
      End If  ! L_mom

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------

      Do i = 1,n_dp
        cca_2d(i) = 0.0
        iccb(i)   = 0
        icct(i)   = 0
        tcw(i)    = 0.0
        cclwp(i)  = 0.0
        lcca(i)   = 0.0
        lctop(i)  = 0
        lcbase(i) = 0
        lcclwp(i) = 0.0
      End Do

!-----------------------------------------------------------------------
! 2.4  Initialise gridbox mean diagnostics
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2.5  Initialise diagnostics for scaling and closure calculations
!-----------------------------------------------------------------------

      Do i = 1,n_dp
        flx_init(i)     = 0.0
        flx_init_new(i) = 0.0
        cape(i)      = 0.0
        cape_out(i)  = 0.0
        dcpbydt(i)   = 0.0
        max_cfl(i)   = 0.0
        det_lev(i)   = 0
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

!
! Set rhum to an arbitrary value since it is not used in LAYER_CN
! for deep convection. Moved outside level loop as values unchanged.
!

!        rhum(i) = 0.0
        rhum(i) = q(i,1)/qse(i,1)
      End Do

!
! set SCM adaptive diagnostics for level k = 1
!
      Do i = 1,n_dp
        rbuoy_p_out(i,1) = 0.0
        the_out(i,1) = th(i,1)
        thp_out(i,1) = th(i,1)
        qe_out(i,1) = q(i,1)
        qp_out(i,1) = q(i,1)
      End Do
!Also, initialise rbuoy_p_old, i.e. for previous level,
! to 0.0 for level1
      Do i = 1,n_dp
        rbuoy_p_old(i) = 0.0
      End Do
!initialise ekm14
      Do i =1, n_dp
         ekm14(i) =0.0
      End do
!Initialise adaptive entrainment variables
!intitilaise to level 2 'cos that's where parcel lift starts from
      Do i = 1, n_dp
        thek(i)=th(i,2)
        qek(i)=q(i,2)
        qsek(i)=qse(i,2)
        thekp1(i)=th(i,2)
        qekp1(i)=q(i,2)
        qsekp1(i)=qse(i,2)
        thpk(i)=th(i,2)
        qpk(i)=q(i,2)
        bwk(i)=bwater(i,2)
        bwkp1(i)=bwater(i,2)
!Note that unlike p_layer_boundaries, where k indexing is offset
!by one compared to the dynamics numbering, z retains the numbering
!convention for dynamics variables i.e. for theta levels, k->k
!and for rho levels k+1/2 -> k+1
!check this with Rachel
        zk(i) = z_theta(i,2)
        zkp12(i)=z_rho(i, 3)
        zkp1(i)=z_theta(i, 3)
      End Do

!intialise parcel values over all levels to make sure we have
!non-garbage values at points which don't convect (do not intialise
!inside loop, 'cos some values are set at the end of each pass at k+1
!and must not be overwritten at k at the start of the next pass).
       Do k = 2, nlev -1
         Do i = 1, n_dp
           qp(i,k) = q(i,k)
           thp(i,k) = th(i,k)
         End Do
       End Do


!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------

      Do k = 2,nlev-1

        Do i = 1, n_dp
          rhum(i) = q(i,k)/qse(i,k)
        End Do

!
!  Initialise adaptive diagnostics for this level
!
        Do i = 1,n_dp
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
        Do i = 1, n_dp
          thek(i)=th(i,k)
          qek(i)=q(i,k)
          qsek(i)=qse(i,k)
          thekp1(i)=th(i,k+1)
          qekp1(i)=q(i,k+1)
          qsekp1(i)=qse(i,k+1)
!note for all levels above k=2, thp(i, k) for this current pass
!will have been set as thp(i, k+1) at end of previous pass through
!loop, and similarly for qp (where point has not convected before
!will be set to environmental value)
          If(k  ==  2) then      !set to environmental values
            thpk(i)=th(i,2)
            qpk(i)=q(i,2)
          Else
            thpk(i)=thp(i,k)
            qpk(i)=qp(i,k)
          End if
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
!
! Set initial parcel properties for thp, qp if convection
! is not occurring at level k
!
        Do i = 1,n_dp
! not convecting and not convected in column before
!NB Calc for thp, qp have been moved to before call to LAYER_CN for
!adaptive mod,
!and settings for thpk, qpk added, but should not change other calcs
!which are carried out as before after call to layer_cn.

          If ( .not. bconv(i).and.det_lev(i) == 0) then
            thpi(i)  = th(i,k) + thpixs_v(i,k)
            thp(i,k) = thpi(i)
            thpk(i)=thp(i,k)         !set thpk for first convecting
                                     !level
            qpi(i)   = q(i,k) + qpixs_v(i,k)
            qp(i,k)  = qpi(i)
            qpk(i) = qp(i,k)         !set qpk for first convecting level
          End If
        End Do  ! n_dp


!-----------------------------------------------------------------------
! 3.1  Calculate layer dependent constants (pressure,
!      layer thickness, entrainment coefficients, detrainment
!      coefficients)
!-----------------------------------------------------------------------

! DEPENDS ON: layer_cn
        Call LAYER_CN(k,n_dp,nlev                                       &
      ,                   mdet_on, ent_on                               &
      ,                   ntml,ntpar                                    &
      ,                   .false.,.false.,.true.                        &
      ,                   bconv,bwk,bwkp1                               &
      ,                   exner_layer_boundaries                        &
      ,                   exner_layer_centres                           &
      ,                   p_layer_boundaries,p_layer_centres            &
      ,                   recip_pstar,rhum, zk, zkp12, zkp1             &
      ,                   thek, qek,qsek, thekp1,qekp1,qsekp1           &
      ,                   thpk,qpk ,ekm14                               &
      ,                   pkp1,delpkp1,exkp1                            &
      ,                   pk,delpk,delpkp12,exk,delexkp1                &
      ,                   delp_uv_k, delp_uv_kp1                        &
      ,                   ekp14,ekp34,amdetk)

! Set ekm14 for next pass through loop
         Do i = 1, n_dp
           ekm14(i) = ekp14(i)
         End Do

!
! Calculate dqs/dth for layers k and k+1 (subroutine DQS_DTH)
!

        If (k == 2) then
! DEPENDS ON: dqs_dth
          Call DQS_DTH(dqsthk,k,th(1,k),qse(1,k),exk,n_dp)
        else
          Do i = 1,n_dp
            dqsthk(i) = dqsthkp1(i)
          End Do
        End If

! DEPENDS ON: dqs_dth
        Call DQS_DTH(dqsthkp1,k+1,th(1,k+1),qse(1,k+1),exkp1,n_dp)

!
! Set other grid dependent constants
!

        Do i = 1,n_dp

!
! Maximum initial convective mass flux
!

          flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)
        End Do  ! n_dp

!
! Set initial parcel properties (theta,q,tracer,momentum) if convection
! is not occurring at level k
!
        Do i = 1,n_dp
! not convecting and not convected in column before
! PC2 qclp and qcfp zero at this point but will add an initial
! value at cloud base

          If ( .not. bconv(i).and.det_lev(i) == 0) then
            expi(i)  = exk(i)
            xpk(i,k)  = 0.0
            qclp(i,k) = 0.0
            qcfp(i,k) = 0.0
            flx(i,k) = 0.0
            bgmk(i)  = .false.
            depth(i) = 0.0

            If (L_mom_gk) then  ! Gregory Kershaw CMT

! Set initial parcel values at cloud base to values of near surface winds
                up(i,k) = u(i,nstart(i))
                vp(i,k) = v(i,nstart(i))

            End If
          End If
        End Do  ! n_dp

        If (L_tracer) then
          Do ktra=1,ntra
            Do i = 1,n_dp
              If ( .not. bconv(i)) then
                 trap(i,k,ktra)  = tracer(i,k,ktra)
              End If  !not bconv
            enddo
          End Do
        End If

!
! Carry out initial test to see if convection is possible from layer
! k to k+1. Set bcposs = .T. if
! 1. the point was convecting (bconv = .T.) and did not terminate
! in the previous layer  OR
! 2. k = ntml
!
        Do i = 1,n_dp

          bcposs(i) = bconv(i) .or. k  ==  ntml(i)

        End Do  ! n_dp

!
! Calculate number of points which may convect (ncposs) and
! set compression indices (index1)
!

        ncposs = 0
        Do i = 1,n_dp
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
            cumulus_c(i)= .true.          ! bgmkp1_c out
            uek_c(i)    = u(index1(i),k)
            uekp1_c(i)  = u(index1(i),k+1)
            vek_c(i)    = v(index1(i),k)
            vekp1_c(i)  = v(index1(i),k+1)
            upk_c(i)    = up(index1(i),k)
            vpk_c(i)    = vp(index1(i),k)
          End Do
          If (L_q_interact) then
          Do i = 1,ncposs
!           PC2 variables
            qclek_c(i)    = qcl(index1(i),k)
            qclekp1_c(i)  = qcl(index1(i),k+1)
            qcfek_c(i)    = qcf(index1(i),k)
            qcfekp1_c(i)  = qcf(index1(i),k+1)
            qclpk_c(i)    = qclp(index1(i),k)
            qcfpk_c(i)    = qcfp(index1(i),k)
          End Do
          End If
          If (L_tracer) then
          Do ktra = 1,ntra
            Do i = 1,ncposs
              trapk_c(i,ktra) = trap(index1(i),k,ktra)
              traek_c(i,ktra) = tracer(index1(i),k,ktra)
              traekp1_c(i,ktra) = tracer(index1(i),k+1,ktra)
            End Do
          End Do
          End If
        End If  ! ncposs>0


!-----------------------------------------------------------------------
! 3.2  Lift parcel from layer k to layer k+1
!-----------------------------------------------------------------------

! DEPENDS ON: lift_par
        Call LIFT_PAR(ncposs,n_dp,thpkp1_c,qpkp1_c,xsqkp1_c,            &
                 bgmkp1_c,bwkp1_c,bwk_c,thpk_c,qpk_c,xpk_c,thekp1_c,    &
                 qekp1_c,thek_c,qek_c,qsekp1_c,                         &
                 qclpkp1_c,qclpk_c,qclekp1_c,qclek_c,l_q_interact,      &
                 qcfpkp1_c,qcfpk_c,qcfekp1_c,qcfek_c,                   &
                 pk_c,pkp1_c,exkp1_c,ekp14_c,ekp34_c,L_mom_gk,upkp1_c,  &
                 vpkp1_c,upk_c,vpk_c,uek_c,uekp1_c,vek_c,vekp1_c,       &
                 L_tracer,ntra,trapkp1_c,trapk_c,traekp1_c,             &
                 traek_c)

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
! Allow parcel to convect from ntml.
! Flag convecting points with logical array bconv
!

          If (k  ==  ntml(index1(i))) then
            bconv(index1(i)) = .true.  ! convection active
            blowst(index1(i)) = .true. ! convection initialised in layer

!
! Set parcel mass flux (UM Documentation paper 27, section 1.5)
! If mass flux out of the initial layer is greater than the mass flux
! of the layer over the timestep then limit mass flux to mass of layer.
!

            flxk_c(i) = mb(index1(i)) * G *                             &
     &                      p_layer_centres(index1(i),k) / (R *         &
     &                        thpk_c(i) * (p_layer_centres(index1(i),k) &
     &                           / PREF)**KAPPA)

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

          Else     ! k=ntml test
            blowst(index1(i))=.false. ! not initialised in this layer
          End If   ! k=ntml test

!
! Reset threshold for forced detrainment
! to the initial (positive or negative) buoyancy (limit positive buoy.
! threshold to XSBMIN fn(delta P)), ONLY for first 5 levels of lift
!

          If (k  >=  ntml(index1(i)) .and.                              &
     &        k  <=  ntml(index1(i)) + 4) then

            xsbmin_v(index1(i),k) = min ( xsbmin_v(index1(i),k),        &
     &      - 0.5 * ( th(index1(i),ntml(index1(i))+1) *                 &
     &      (1.0 + C_VIRTUAL * q(index1(i),ntml(index1(i)) + 1))        &
     &      - th(index1(i),ntml(index1(i))) * (1.0 +                    &
     &      C_VIRTUAL * q(index1(i),ntml(index1(i))))  ) + 0.5 )

          End If

        End Do  !ncposs

! L_q_interact_if0:
      If (L_Q_INTERACT) then
! ----------------------------------------------------------------------
!     Follow-on calculation from QCLP(*,K) and QCFP(*,K) initialization.
!       Duplicates code from Convection Parcel Lifting Scheme, LIFPAR.
! ----------------------------------------------------------------------
! ntml is level from which convection starts [dp/sh conv only]

        Do I=1, ncposs
          If (k  ==  ntml(index1(i))) then
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
            uek_c2(i)    = uek_c(index2(i))
            uekp1_c2(i)  = uekp1_c(index2(i))
            vek_c2(i)    = vek_c(index2(i))
            vekp1_c2(i)  = vekp1_c(index2(i))
            dqsthkp1_c2(i) = dqsthkp1(index1(index2(i)))
            qsekp1_c2(i)   = qsekp1_c(index2(i))
            pstar_c2(i)    = pstar(index1(index2(i)))
            thpkp1_c2(i)   = thpkp1_c(index2(i))
            qpkp1_c2(i)    = qpkp1_c(index2(i))
            upkp1_c2(i)    = upkp1_c(index2(i))
            vpkp1_c2(i)    = vpkp1_c(index2(i))
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
            blowst_c2(i) = blowst(index1(index2(i)))
            L_shallow_c2(i) = .false.
            L_mid_c2(i)     = .false.
            cumulus_c2(i)   = .true.
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
          End Do
          If (L_q_interact) then
          Do i = 1,nconv
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
            start_lev3c(I) = ntml(index1(index2(i)))
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
          Do i = 1,nconv
            tcw_c2(i)    = tcw(index1(index2(i)))
            depth_c2(i)  = depth(index1(index2(i)))
            cclwp_c2(i)  = cclwp(index1(index2(i)))
            cape_c2(i)    = cape(index1(index2(i)))
            dcpbydt_c2(i) = dcpbydt(index1(index2(i)))
            relh_c2(i)    = relh(index1(index2(i)))
            dptot_c2(i)   = dptot(index1(index2(i)))
            rbuoy_p_here_c2(i)=rbuoy_p_here(index1(index2(i)))
            the_here_c2(i)=the_here(index1(index2(i)))
            thp_here_c2(i)=thp_here(index1(index2(i)))
            qe_here_c2(i)=qe_here(index1(index2(i)))
            qp_here_c2(i)=qp_here(index1(index2(i)))
            rbuoy_p_old_c2(i)=rbuoy_p_old(index1(index2(i)))
            thpixs_v_c2(i) = thpixs_v(index1(index2(i)),k)
            qpixs_v_c2(i)  = qpixs_v(index1(index2(i)),k)
            xsbmin_v_c2(i) = xsbmin_v(index1(index2(i)),k)
            iccb_c2(i)     = iccb(index1(index2(i)))
            icct_c2(i)     = icct(index1(index2(i)))
            cca_2d_c2(i)   = cca_2d(index1(index2(i)))
            lcca_c2(i)=lcca(index1(index2(i)))
            lcbase_c2(i)=lcbase(index1(index2(i)))
            lctop_c2(i)=lctop(index1(index2(i)))
            lcclwp_c2(i)=lcclwp(index1(index2(i)))
          End Do
          If (L_tracer) then
          Do ktra = 1,ntra
            Do i = 1,nconv
              trapk_c2(i,ktra) = trapk_c(index2(i),ktra)
            dtraek_c2(i,ktra) = dtrabydt(index1(index2(i)),k,ktra)
            dtraekp1_c2(i,ktra) = dtrabydt(index1(index2(i)),k+1,ktra)
            End Do
          End Do
          End If
        End If  ! nconv>0

!-----------------------------------------------------------------------
! 3.3  Calculate the rest of the parcel ascent  and the effect of
!      convection on the large-scale atmosphere.
!
!      Subroutine CONVEC2
!
!      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
!-----------------------------------------------------------------------

! DEPENDS ON: convec2
        Call CONVEC2(nconv,n_dp,nlev,k,thek_c2,thekp1_c2,qek_c2,        &
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
     &               dcpbydt_c2,max_cfl_c2,timestep,                    &
     &               rbuoy_p_here_c2,                                   &
     &               the_here_c2,thp_here_c2,qe_here_c2,qp_here_c2,     &
     &               rbuoy_p_old_c2, ad_on, sdet_on, new_termc,         &
     &               deltak_c2,l_calc_dxek,l_q_interact,                &
     &               flxkp12_c2,cumulus_c2,relh_c2,dptot_c2             &
     &               )

!
! Calculate fractional entrainment rate for level k.
! If convection has terminated (bterm=.T.) then set
! fractional entrainment rate for k+1 to zero.
!

        If (nconv  >   0) then
        If (flg_entr_up) then
          Do i = 1,nconv
            entrain_up(index1(index2(i)),k) = (1.0 - deltak_c2(i)) *    &
     &           (1.0 - amdetk_c2(i)) * (ekp14_c2(i) + ekp34_c2(i)      &
     &           * (1.0 + ekp14_c2(i))) * flx(index1(index2(i)),k)
            If (bterm_c2(i)) then
              entrain_up(index1(index2(i)),k+1) = 0.0
            End If
          End Do
        End If

!
! Calculate fractional detrainment rate for level k
! (and k+1 if bterm=.T.)
!

        If (flg_detr_up) then
          Do i = 1,nconv
            detrain_up(index1(index2(i)),k) = -(amdetk_c2(i)            &
     &               + deltak_c2(i) * (1.0 - amdetk_c2(i)))             &
     &                      * flx(index1(index2(i)),k)
            If (bterm_c2(i)) then
              detrain_up(index1(index2(i)),k+1) =                       &
     &           -(1.0 - deltak_c2(i)) * flx(index1(index2(i)),k)
            End If
          End Do
        End If
        End If

!
! Write CONVEC2 compressed output arrays back to full fields
!


        Do i = 1,n_dp
          xpk(i,k+1)    = 0.0
          flx(i,k+1)    = 0.0
          depth(i)      = 0.0
          precip(i,k+1) = 0.0
          qclp(i,k+1)   = 0.0
          qcfp(i,k+1)   = 0.0
        End Do
        Do i = 1,n_dp
          bgmk(i)       = .false.
          bterm(i)      = .false.
        End Do

        If (L_tracer) then
          Do ktra = 1,ntra
            Do i = 1,n_dp
              trap(i,k+1,ktra) = 0.0
            End Do
          End Do
        End If

        If (L_mom_gk) then
          Do i = 1,n_dp
            up(i,k+1) = 0.0
            vp(i,k+1) = 0.0
          End Do
        End If

        If (nconv  >   0)then
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
          qe_here(index1(index2(i))) = qe_here_c2(i)
          qp_here(index1(index2(i))) = qp_here_c2(i)
         End Do
         Do i = 1,nconv
               max_cfl(index1(index2(i))) =                             &
     &                   MAX(max_cfl(index1(index2(i))),max_cfl_c2(i))
               relh(index1(index2(i)))  = relh_c2(i)
               dptot(index1(index2(i))) = dptot_c2(i)
         End Do

        If (L_q_interact) then
         Do i = 1,nconv
!          PC2 variables
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
        If (L_mom_gk) then ! Gregory Kershaw CMT

          Do i = 1,nconv
            up(index1(index2(i)),k+1) = upk_c2(i)
            vp(index1(index2(i)),k+1) = vpk_c2(i)
          End Do

          Do i = 1,nconv
            dubydt(index1(index2(i)),k) = duek_c2(i)
            dvbydt(index1(index2(i)),k) = dvek_c2(i)
            dubydt(index1(index2(i)),k+1) = duekp1_c2(i)
            dvbydt(index1(index2(i)),k+1) = dvekp1_c2(i)
            eflux_u_ud(index1(index2(i))) = eflux_u_ud_c2(i)
            eflux_v_ud(index1(index2(i))) = eflux_v_ud_c2(i)
          End Do
        End If  ! L_mom_gk

        If (L_mom) then    ! needed for all versions
          Do i = 1,nconv
           flxkp12(k,index1(index2(i)))  = flxkp12_c2(i)
          End Do
        End If  ! L_mom

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
          Do i = 1 ,nconv
            up_flux(index1(index2(i)),k+1) = flxk_c2(i)
          End Do
        End If
        End If
!
!   Write adaptive diagnostics for this level to full array for output
!
        Do i = 1,n_dp
          rbuoy_p_out(i,k) = rbuoy_p_here(i)
          the_out(i,k) = the_here(i)
          thp_out(i,k) = thp_here(i)
          qe_out(i,k) = qe_here(i)
          qp_out(i,k) = qp_here(i)
        End Do


!   Write rbuoy for this level to rbuoy_p_old for previous level
!   for use in parcel
!
        Do i = 1,n_dp
          rbuoy_p_old(i) = rbuoy_p_here(i)
        End Do


!-----------------------------------------------------------------------
! 3.4  Cape and CFL scaling - adjust initial mass flux so that cape is
!      removed by convection over timescale cape_timescale.
!-----------------------------------------------------------------------
!
! Set up integer nterm which is the total number of points where
! convection has terminated.
! Index to full array (n_dp) with index_nterm
!

        nterm = 0
        Do i = 1,n_dp
          If (bterm(i)) then
            nterm = nterm + 1
            index_nterm(nterm) = i
          End If
        End Do

! Note nterm is often small so vectorisation of the loop over nterm
! will not always result in faster code.

        If (nterm >  0) then

         If (cape_opt == 0) then
!default 4a convection scheme - RH-based CAPE closure

          Do j = 1,nterm
            i = index_nterm(j)
            If (dcpbydt(i)  >   0.0) then
              cape_ts_new =                                             &
     &        MIN(MAX(900.0*(1.0 - relh(i)/dptot(i))/0.1,60.0)          &
     &        ,cape_timescale)

              flx_init_new(i) = flx_init(i)*cape(i)                     &
     &                            /(cape_ts_new*dcpbydt(i))

              If (flx_init_new(i)  >   flxmax_init(i)) then
                flx_init_new(i) = flxmax_init(i)
              End If

!
! Scale max_cfl with cape scale
!

              max_cfl(i) = max_cfl(i) * flx_init_new(i) / flx_init(i)
            else
              flx_init_new(i) = flx_init(i)
            End If  ! dcpbydt > 0

! If convection has terminated write cape to diagnostic output
! variable (cape_out).

            cape_out(i) = cape(i)
            cape(i)     = 0.0
            dcpbydt(i)  = 0.0
          End Do   ! nterm

         Else If(cape_opt  ==  1) then

! modified 4a convection scheme - RH-based CAPE closure, timescale
! limited to timestep

            Do j = 1,nterm
              i = index_nterm(j)
              If (dcpbydt(i)  >   0.0) then
                  cape_ts_new =                                         &
     &            MIN(MAX(cape_timescale*(1.0-relh(i)/dptot(i))/0.4            &
     &            ,timestep)                                            &
     &            ,cape_timescale)

                flx_init_new(i) = flx_init(i)*cape(i)                   &
     &                              /(cape_ts_new*dcpbydt(i))

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
              cape_out(i) = cape(i)
              cape(i)     = 0.0
              dcpbydt(i)  = 0.0
            End Do   ! nterm

         Else If(cape_opt  ==  2) then
! Switch off RH-based CAPE closure, no test on w_max

            Do j = 1,nterm
              i = index_nterm(j)
              If (dcpbydt(i)  >   0.0) then
                cape_ts_new =cape_timescale
                flx_init_new(i) = flx_init(i)*cape(i)                   &
     &                              /(cape_ts_new*dcpbydt(i))

                If (flx_init_new(i)  >   flxmax_init(i)) then
                  flx_init_new(i) = flxmax_init(i)
                End If
! Scale max_cfl with cape scale
!

                max_cfl(i) = max_cfl(i) * flx_init_new(i)/flx_init(i)
              else
                flx_init_new(i) = flx_init(i)
              End If  ! dcpbydt > 0

! If convection has terminated write cape to diagnostic output
! variable (cape_out).

              cape_out(i) = cape(i)
              cape(i)     = 0.0
              dcpbydt(i)  = 0.0
            End Do   ! nterm

        End if    ! test on cape_opt = 0,1,2


        If((cape_opt  == 3) .or. (cape_opt  ==  4) .or.                 &
     &     (cape_opt  ==  5)) then

          If ( w_cape_limit < 1000.0 ) then
!  This section includes test on w_max
            IF (cape_opt  == 3) THEN
              Do j = 1, nterm
                i = index_nterm(j)
                If ( dcpbydt(i) > 0.0 ) then
! new denominator introduced at vn6.6
                  if ( W_MAX(I) > w_cape_limit ) then
                    cape_ts_new =   CAPE_TS_w * w_cape_limit/           &
     &              (w_cape_limit+ (W_MAX(I)-w_cape_limit)*wcape_fac)
                  else
                    cape_ts_new = cape_timescale
                  endif !  W_MAX(I) > w_cape_limit


                  flx_init_new(i) = flx_init(i) * cape(i) /             &
     &                               ( cape_ts_new * dcpbydt(i) )

                  If (flx_init_new(i) > flxmax_init(i)) then
                    flx_init_new(i) = flxmax_init(i)
                  End If  ! flx_init_new(i) > flxmax_init(i)

! Scale max_cfl with cape scale

                  max_cfl(i) = max_cfl(i) * flx_init_new(i)             &
     &                   / flx_init(i)
                else    !dcpbydt(i) <= 0.0
                  flx_init_new(i) = flx_init(i)
                End If  ! dcpbydt(i) > 0.0

! If convection has terminated write cape to diagnostic output
! variable (cape_out).

                cape_out(i) = cape(i)
                cape(i)     = 0.0
                dcpbydt(i)  = 0.0

              End Do  ! j = 1, nterm

            ELSE IF ( cape_opt == 4) THEN
              Do j = 1, nterm
                i = index_nterm(j)
                If ( dcpbydt(i) > 0.0 ) then

                  if ( W_MAX(I) > w_cape_limit ) then
                    cape_ts_new =   CAPE_TS_w * w_cape_limit/W_MAX(I)
                  else ! do not test on w_max
                    cape_ts_new = CAPE_TIMESCALE * CAPE(I) / CAPE_min +        &
     &                            CAPE_TIMESCALE * EXP( - CAPE(I) /CAPE_min)
                  endif !  W_MAX(I) > w_cape_limit

                  flx_init_new(i) = flx_init(i) * cape(i) /             &
     &                               ( cape_ts_new * dcpbydt(i) )

                  If (flx_init_new(i) > flxmax_init(i)) then
                    flx_init_new(i) = flxmax_init(i)
                  End If  ! flx_init_new(i) > flxmax_init(i)

! Scale max_cfl with cape scale

                  max_cfl(i) = max_cfl(i) * flx_init_new(i)             &
     &                        / flx_init(i)
                else    !dcpbydt(i) <= 0.0
                  flx_init_new(i) = flx_init(i)
                End If  ! dcpbydt(i) > 0.0

! If convection has terminated write cape to diagnostic output
! variable (cape_out).

                cape_out(i) = cape(i)
                cape(i)     = 0.0
                dcpbydt(i)  = 0.0

              End Do  ! j = 1, nterm

            ELSE IF ( cape_opt == 5 ) THEN
              Do j = 1, nterm
                i = index_nterm(j)
                If ( dcpbydt(i) > 0.0 ) then

                  if ( W_MAX(I) > w_cape_limit ) then
                    cape_ts_new =   CAPE_TS_w * w_cape_limit /W_MAX(I)
                  else
                    if ( relh(i) / dptot(i) >= 0.75 ) then
                      CAPE_TS_NEW = CAPE_TIMESCALE *                           &
     &                            ( 0.2373 / (relh(i) / dptot(i))**5)
                    else
                      CAPE_TS_NEW = CAPE_TIMESCALE
                    endif ! relh(i) / dptot(i) >= 0.75
                  endif !  W_MAX(I) > w_cape_limit

                  flx_init_new(i) = flx_init(i) * cape(i) /             &
     &                               ( cape_ts_new * dcpbydt(i) )

                  If (flx_init_new(i) > flxmax_init(i)) then
                    flx_init_new(i) = flxmax_init(i)
                  End If  ! flx_init_new(i) > flxmax_init(i)

! Scale max_cfl with cape scale

                  max_cfl(i) = max_cfl(i) * flx_init_new(i)             &
     &                   / flx_init(i)
                else    !dcpbydt(i) <= 0.0
                  flx_init_new(i) = flx_init(i)
                End If  ! dcpbydt(i) > 0.0

! If convection has terminated write cape to diagnostic output
! variable (cape_out).

                cape_out(i) = cape(i)
                cape(i)     = 0.0
                dcpbydt(i)  = 0.0

              End Do  ! j = 1, nterm


            End If ! CAPE_opt 3,4,5

          Else
!  This section as above except no test on w_max

            IF ( cape_opt == 3 ) THEN
              Do j = 1, nterm
                i = index_nterm(j)
                If ( dcpbydt(i) > 0.0 ) then

                  flx_init_new(i) = flx_init(i) * cape(i) /             &
     &                               ( cape_timescale * dcpbydt(i) )

                  If (flx_init_new(i) > flxmax_init(i)) then
                    flx_init_new(i) = flxmax_init(i)
                  End If  ! flx_init_new(i) > flxmax_init(i)

! Scale max_cfl with cape scale

                  max_cfl(i) = max_cfl(i) * flx_init_new(i)/flx_init(i)
                else    !dcpbydt(i) <= 0.0
                  flx_init_new(i) = flx_init(i)
                End If  ! dcpbydt(i) > 0.0

! If convection has terminated write cape to diagnostic output
! variable (cape_out).

                cape_out(i) = cape(i)
                cape(i)     = 0.0
                dcpbydt(i)  = 0.0

              End Do  ! j = 1, nterm

            ELSE IF ( cape_opt == 4 ) THEN
              Do j = 1, nterm
                i = index_nterm(j)
                If ( dcpbydt(i) > 0.0 ) then

                  cape_ts_new = CAPE_TIMESCALE * CAPE(I) / CAPE_min +          &
     &                            CAPE_TIMESCALE * EXP( - CAPE(I) / CAPE_min )

                  flx_init_new(i) = flx_init(i) * cape(i) /             &
     &                               ( cape_ts_new * dcpbydt(i) )

                  If (flx_init_new(i) > flxmax_init(i)) then
                    flx_init_new(i) = flxmax_init(i)
                  End If  ! flx_init_new(i) > flxmax_init(i)

! Scale max_cfl with cape scale

                  max_cfl(i) = max_cfl(i) * flx_init_new(i)/flx_init(i)
                else    !dcpbydt(i) <= 0.0
                  flx_init_new(i) = flx_init(i)
                End If  ! dcpbydt(i) > 0.0

! If convection has terminated write cape to diagnostic output
! variable (cape_out).

                cape_out(i) = cape(i)
                cape(i)     = 0.0
                dcpbydt(i)  = 0.0

              End Do  ! j = 1, nterm

            ELSE IF ( cape_opt == 5 ) THEN
              Do j = 1, nterm
                i = index_nterm(j)
                If ( dcpbydt(i) > 0.0 ) then

                  if ( relh(i) / dptot(i) >= 0.75 ) then
                    CAPE_TS_NEW = CAPE_TIMESCALE *                             &
     &                            ( 0.2373 / (relh(i) / dptot(i))**5)
                  else
                    CAPE_TS_NEW = CAPE_TIMESCALE
                  endif ! relh(i) / dptot(i) >= 0.75

                  flx_init_new(i) = flx_init(i) * cape(i) /             &
     &                               ( cape_ts_new * dcpbydt(i) )

                  If (flx_init_new(i) > flxmax_init(i)) then
                    flx_init_new(i) = flxmax_init(i)
                  End If  ! flx_init_new(i) > flxmax_init(i)

! Scale max_cfl with cape scale

                  max_cfl(i) = max_cfl(i) * flx_init_new(i)/flx_init(i)
                else    !dcpbydt(i) <= 0.0
                  flx_init_new(i) = flx_init(i)
                End If  ! dcpbydt(i) > 0.0

! If convection has terminated write cape to diagnostic output
! variable (cape_out).

                cape_out(i) = cape(i)
                cape(i)     = 0.0
                dcpbydt(i)  = 0.0

              End Do  ! j = 1, nterm


            End If ! cape_opt test 3,4,5

          End If !  w_cape_limit < 1000.0
        End If ! cape_opt = 3 or 4 or 5

!
! Work out scaled mass flux needed to keep cfl ratio below limit.
! Note CAPE closure assumed

          Do j = 1,nterm
            i = index_nterm(j)
            max_cfl(i) = max_cfl(i) * timestep

            If (max_cfl(i)  >   CFL_LIMIT) then
              flx_init_new(i) = flx_init_new(i) * CFL_LIMIT             &
     &                                       / max_cfl(i)
            else
              flx_init_new(i) = flx_init_new(i)
            End If

            If (flx_init_new(i)  >   flxmax_init(i)) then
              flx_init_new(i) = flxmax_init(i)
            End If
            max_cfl(i) = 0.0
          End Do      ! j (nterm)

          !
          ! Scale cloud fraction
          !
          If (lcv_ccrad) Then

            Do j = 1,nterm
              i = index_nterm(j)

              If (flx_init_new(i) > 0.0) Then
                scale_f(i) = flx_init_new(i) / flx_init(i)
                cca_2d(i)  = cca_2d(i) + 0.06 * LOG(scale_f(i))

                ! set flx_init to the new value to provide the real initial mass
                ! flux in all conditions
                flx_init(i) = flx_init_new(i)
              End If

              ! Check scaled cloud fraction not smaller than minimum value
              ! (2.0E-5) or greater than unity.
              !
              ! (Was moved out of scaling if test to ensure these limits
              ! at all times, not just when cca_2d is scaled)

              cca_2d(i) = MAX(2.0E-5, cca_2d(i))
              cca_2d(i) = MIN(1.0E+0, cca_2d(i))


            End Do      ! j (nterm)

          Else      ! original

            Do j = 1,nterm
              i = index_nterm(j)
              If (flx_init_new(i) > 0.0) Then
                scale_f(i) = flx_init_new(i) / flx_init(i)
                cca_2d(i)  = cca_2d(i) + 0.06 * LOG(scale_f(i))

                !
                ! Check scaled cloud fraction not smaller than minimum value
                ! (2.0E-5) or greater than unity.
                !
                cca_2d(i) = MAX(2.0E-5, cca_2d(i))

                If (cca_2d(i) > 1.0) Then
                  cca_2d(i) = 1.0
                End If  

                ! set flx_init to the new value to provide the real initial mass
                ! flux in all conditions
                flx_init(i) = flx_init_new(i)

              End If
            End Do      ! j (nterm)
          End If      ! lcv_ccrad


!
! Carry out cape and cfl scaling
!

          Do kt = 2, k+1
            Do j = 1,nterm
              i = index_nterm(j)
              If (kt  >=  ntml(i) .and. flx_init_new(i) >   0.0) then

                dthbydt(i,kt) = dthbydt(i,kt) * scale_f(i)
                dqbydt(i,kt)  = dqbydt(i,kt) * scale_f(i)
                If (L_q_interact) then  ! PC2
                  dqclbydt(i,kt)  = dqclbydt(i,kt) * scale_f(i)
                  dqcfbydt(i,kt)  = dqcfbydt(i,kt) * scale_f(i)
                  dcflbydt(i,kt)  = dcflbydt(i,kt) * scale_f(i)
                  dcffbydt(i,kt)  = dcffbydt(i,kt) * scale_f(i)
                  dbcfbydt(i,kt)  = dbcfbydt(i,kt) * scale_f(i)
                End If
                If (L_mom_gk) then
                  dubydt(i,kt)  = dubydt(i,kt) * scale_f(i)
                  dvbydt(i,kt)  = dvbydt(i,kt) * scale_f(i)
                End If
                If (L_mom) then     ! required for all versions
                  If (kt <  k+1) then
                    flxkp12(kt,i) = flxkp12(kt,i) * scale_f(i)
                  Endif
                End If
                If (L_tracer) then
                  Do ktra = 1,ntra
                   dtrabydt(i,kt,ktra) = dtrabydt(i,kt,ktra)*scale_f(i)
                  End Do
                End If

                flx(i,kt)    = flx(i,kt) * scale_f(i)
                precip(i,kt) = precip(i,kt) * scale_f(i)

                If (flg_up_flx) then
                  up_flux(i,kt) = flx(i,kt)
                End If
                If (flg_entr_up) then
                  entrain_up(i,kt) = entrain_up(i,kt) * scale_f(i)
                End If
                If (flg_detr_up) then
                  detrain_up(i,kt) = detrain_up(i,kt) * scale_f(i)
                End If


             End If !kt>ntml and flx_init_new>0
            End Do  ! j loop
          End Do  ! kt loop

!
! Set cape scaling parameters for next level
!

          Do j = 1,nterm
            i = index_nterm(j)
            det_lev(i)= k+1
          End Do  ! nterm loop

!-----------------------------------------------------------------------
! 3.5  Downdraft calculation - on all points where convection is
!      terminating.
!
!      Subroutine DD_CALL
!
!      UM Documentation Paper 27, part 2
!
!-----------------------------------------------------------------------

        npossdd = 0
        nnodd = 0
        Do j = 1,nterm
          i = index_nterm(j)
          tempnum = 0.0
          If (iccb(i)  >   0) then
            deltap_cld = p_layer_centres(i,iccb(i))                     &
     &                           - p_layer_centres(i,k)
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
! Set logical for use later

          If (deltap_cld  >   15000.0 .and. bgmk(i) .and.               &
     &                                tempnum  >   1E-12) then
            b_dd(i) = .true.
          else
            b_nodd(i) = .true.
          End If
        End Do  ! nterm

!
! If convection has terminated write cape to diagnostic output
! variable (cape_out).
! Set kterm array which holds the level index for termination
! of convection.
!

          Do j = 1,nterm
            i=index_nterm(j)
            if (k >= ntml(i)) then
              kterm(i) = k
            endif
            bconv(i) = .false.
          End Do
        End If

!-----------------------------------------------------------------------
! 3.6  End of main loop over levels
!-----------------------------------------------------------------------

      End Do

!-----------------------------------------------------------------------
! 4.0 Down draughts  - now 2 options
!                      original code
!                      Emanuel down draught code
! Note the level at which deep convection terminates has been stored
! in the above updraught loop as Kterm.
!-----------------------------------------------------------------------

      If (l_eman_dd) then


! Work out maximum termination level
          kmax_term = 2
          Do i = 1,n_dp
            If(kterm(i) >  kmax_term) then
              kmax_term = kterm(i)
            End if
          End do

! DEPENDS ON: eman_dd
        call eman_dd (n_dp,kmax_term,nlev,trlev,ntra                    &
     &,                      kterm,l_tracer                             &
     &,                      exner_layer_centres,exner_layer_boundaries &
     &,                      p_layer_centres, p_layer_boundaries        &
     &,                      timestep, th, q, qse, tracer, precip       &
     &,                      dthbydt, dqbydt, dtrabydt                  &
     &,                      rain, snow ,dwn_flux                       &
     &                     )

      Else         ! original down draught code
!-----------------------------------------------------------------------
! 4.1  Downdraft calculation - on all points where convection is
!      terminating.
!
!      Subroutine DD_ALL_CALL
!
!      UM Documentation Paper 27, part 2
!
!-----------------------------------------------------------------------

        npossdd = 0
        Do i = 1,n_dp
          if (b_dd(i)) then
            npossdd = npossdd +1
            index_possdd(npossdd) = i
          End If
        End do

        If (npossdd  >   0) then

! Work out maximum termination level
          kmax_term = 2
          Do i = 1,npossdd
            If(kterm(index_possdd(i)) >  kmax_term) then
              kmax_term = kterm(index_possdd(i))
            End if
          End do

! DEPENDS ON: dd_all_call
         call DD_ALL_CALL (n_dp,npossdd,kmax_term,nlev,trlev,ntra       &
     &,                      kterm, iccb, icct, index_possdd, l_tracer  &
     &,                      flg_dwn_flx, flg_entr_dwn, flg_detr_dwn    &
     &,                      bwater(1,2)                                &
     &,                      exner_layer_centres,exner_layer_boundaries &
     &,                      p_layer_centres, p_layer_boundaries,pstar  &
     &,                      recip_pstar,timestep , cca_2d              &
     &,                      thp, qp, th, q, trap,tracer, flx,precip    &
     &,                      dthbydt, dqbydt, dtrabydt                  &
     &,                      rain, snow , rain_3d, snow_3d, dwn_flux    &
     &,                      entrain_dwn, detrain_dwn)


         End if

!-----------------------------------------------------------------------
! 4.2 Surface precipitation calculation for terminating points with
!     no downdraught (moved outside level loop) ie do this calculation
!     on all points at the end.
!-----------------------------------------------------------------------
! Points where no downdraught possible
        nnodd = 0
        Do i = 1,n_dp

          if (b_nodd(i)) then
            nnodd = nnodd +1
            index_nodd(nnodd) = i
          End If
        End do

        If (nnodd  >   0) then

! Work out maximum termination level
          kmax_term = 2
          Do i = 1,nnodd
            If(kterm(index_nodd(i)) >  kmax_term) then
              kmax_term = kterm(index_nodd(i))
            End if
          End do
! Only add 1 if kmax_term is less than model levels
! (which should be true).
          if (kmax_term  <  nlev ) then
            kmax_term = kmax_term + 1
          End if

!
! Surface precipitation calculation
!

! DEPENDS ON: evap_bcb_nodd_all
          Call EVAP_BCB_NODD_ALL(n_dp,nnodd,kmax_term,kterm             &
     &,                      iccb, index_nodd, bwater(1,2)              &
     &,                      exner_layer_centres,exner_layer_boundaries &
     &,                      p_layer_centres, p_layer_boundaries,pstar  &
     &,                      timestep , cca_2d, th, q, precip           &
     &,                      dthbydt, dqbydt, rain, snow                &
     &,                      rain_3d, snow_3d )
        End If

      Endif        ! test on down draught type

!-----------------------------------------------------------------------
! 4.3 Adjust cloud base, top and amount to prevent errors occurring in
! radiation scheme when iccb = icct (this happens when convection
! saturates upon forced detrainment).
!-----------------------------------------------------------------------
!
      If (.NOT. lcv_ccrad) Then
        Do i=1, n_dp

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


!-----------------------------------------------------------------------
! 5.0  Convective Momentum Transport (if L_mom = .true.)
!-----------------------------------------------------------------------

      If (L_mom) then
          
        Select Case (deep_cmt_opt)
        Case (2)       ! Gregory-Kershaw deep CMT

          ! Do nothing here as calculated in parcel ascent earlier  

        Case Default   ! (0/1) Alan Grant's Eddy viscosity CMT
    
        ! altered to use kterm instead of ntpar
          Do i = 1,n_dp
            If (kterm(i) >=  nlcl_uv(i)) then
              ntop_uv(i)    = kterm(i) + 1
            Else     ! case where deep convection fails
        ! I think in this case the cloud base mass flux will be zero so
        ! there will be no CMT. (The value will not matter)
              ntop_uv(i)    = ntpar(i) + 1
            End If

            ptop_uv(i)    = phalf_uv(ntop_uv(i),i)
          End Do

          nterm = 0
!
! Set cloud base mass flux equal to mass flux at half level below
! the LCL.
!
          Do i = 1,n_dp
            If (kterm(i)  >=  nlcl_uv(i)) then
              nterm = nterm + 1
              cu_term(nterm) = i
              cu_tend(nterm) = i
              mb(i) = flxkp12(nlcl_uv(i),i)
              Do j = 1,nlev
                flxkp12(j,i) = 0.0
              End Do
            End If

! initialise output arrays as lower level subroutines don't set all
! values

            Do j = 1,nlev
              uw(j,i)=0.0
              vw(j,i)=0.0
              visc(j,i)=0.0
            End Do

          End Do

          If (nterm  >   0) then
! DEPENDS ON: cmt_mass
            Call CMT_MASS(n_dp,n_dp,nlev,nterm, cu_term,mb,             &
                          p_0degc_uv,kterm,plcl_uv,ptop_uv,n_0degc,     &
                          nlcl_uv,ntop_uv,phalf_uv,p_uv,cu_tend,        &
                          !OUTPUT
                          flxkp12,mass_dwn,visc)

! DEPENDS ON: deep_grad_stress
            Call DEEP_GRAD_STRESS(n_dp,n_dp,n_dp,nlev,                    &
                                  nlcl_uv,nterm,cu_term,                  &
                                  cu_tend,ue_p,ve_p,visc,phalf_uv,p_uv,   &
                                  rho_uv,timestep,ntop_uv,                &
                                  !OUTPUT
                                  uw,vw)


! DEPENDS ON: deep_ngrad_stress
            Call DEEP_NGRAD_STRESS(n_dp,n_dp,n_dp,nterm,nlev,             &
                             nlcl_uv,cu_term,cu_tend,cu_tend,pstar,uw0,   &
                                   vw0,zlcl_uv,ue_p,ve_p,visc,flxkp12,    &
                                   p_uv,phalf_uv,rho_uv,timestep,ntop_uv, &
                                   flg_uw_deep,flg_vw_deep,               &
                                   !INPUT/OUTPUT
                                   uw,vw,                                 &
                                   uw_base,vw_base,uw_deep,vw_deep)

! DEPENDS ON: deep_cmt_incr
            Call DEEP_CMT_INCR(n_dp,n_dp,n_dp,nlev,nlcl_uv,               &
                               nterm,cu_term,cu_tend,zlcl_uv,phalf_uv,    &
                               p_uv,rho_uv,uw_base,vw_base,uw,vw,ntop_uv, &
                               ! OUTPUT
                               dubydt,dvbydt)

          End If  ! nterm > 0

        Case (3,4)        ! New Turbulence scheme using heights 

          nterm = 0   ! count of number of deep points which actually convected
          Do i = 1,n_dp
            If (kterm(i)  >=  nlcl_uv(i)) then
              nterm = nterm + 1
              cu_term(nterm) = i

! Use CAPE scaled mass flux as initial mass flux rather than CRM derived value
! to be consistent with thermodynamic part of convection.
              mb(i) = flxkp12(nlcl_uv(i),i)
              zlcl(i) = z_rho(i,ntml(i))
 
! Cloud velocity scale - derived from CRM simulations
!             wcld = (C_mass*wstar*CAPE)**(1/3)

              wcld(i) = (delthvu(i) * C_MASS * wstar(i) * G / (th(i,ntml(i)) &
                    * (1.0 + C_VIRTUAL * q(i,ntml(i)))))**0.3333
            End If
          End Do 

          If (nterm > 0) Then
! DEPENDS ON: deep_turb_cmt
            call deep_turb_cmt (n_dp, nterm, nlev, deep_cmt_opt,         &
                          ntml, kterm,cu_term,freeze_lev,                &
                          timestep,                                      &
                          uw0, vw0, mb, wcld, wstar ,zlcl,               &
                          flx,                                           &
                          r_rho, r_theta, z_rho, z_theta,rho,rho_theta,  &
                          r2rho, r2rho_th, dr_across_th, dr_across_rh,   &
                          u, v,                                          &
                          dubydt, dvbydt, uw_deep, vw_deep)
  
          End If           ! nterm > 0

        End Select           ! deep_cmt_opt

      End If  ! L_mom

!-----------------------------------------------------------------------
! 6.0  Energy correction calculation - removed as old code not correct
!     for new dynmaics grid ( attempts to correct this give problems).
!     UM documentation paper 27 - section 12.
!-----------------------------------------------------------------------
      Do i = 1,n_dp
          index1(i) = i
      End Do

!      call Cor_engy(n_dp,n_dp,n_dp,nlev,dthbydt,dqbydt,snow
!     &               ,exner_layer_centres,p_layer_boundaries,index1)

!      call Cor_engy(n_dp,n_dp,n_dp,nlev,index1,.false.,r2rho_th
!     &               ,dr_across_th,dqbydt,rain,snow
!     &               ,exner_layer_centres,dthbydt)

!-----------------------------------------------------------------------
! 7.0  Total water conservation  - also works on whole column
!-----------------------------------------------------------------------

! only check columns where convection has occurred.


        Do i = 1,n_dp
          qMinInColumn(i) = q(i,nlev)
        End Do
        Do k = 1,nlev-1
          Do i = 1,n_dp
            If (q(i,k)  <   qMinInColumn(i)) then
              qMinInColumn(i) = q(i,k)
            End If
          End Do
        End Do

!
! Ensure Q does not go below global allowed minimum (QMIN)
!
        Do i = 1,n_dp
          qMinInColumn(i)=MAX(QMIN,qMinInColumn(i))
        End Do

!
! Apply an artificial upwards flux from k-1 level to ensure Q
! remians above minimum value in the column.
!

        Do k = nlev,2,-1
          Do i = 1,n_dp
            If (dqbydt(i,k) /= 0.0) then
              temp1(i)=q(i,k) + dqbydt(i,k) * timestep
              If (temp1(i)  <   qMinInColumn(i)) then

                dqbydt(i,k-1) = dqbydt(i,k-1) -                         &
     &              ((qMinInColumn(i) - q(i,k)) / timestep-dqbydt(i,k)) &
     &               * (r2rho_th(i,k)*dr_across_th(i,k))                &
     &               / (r2rho_th(i,k-1)*dr_across_th(i,k-1))

                dqbydt(i,k) = (qMinInColumn(i) - q(i,k)) / timestep
              End If
            End If
          End Do ! n_dp loop
        End Do  ! nlev
!
! check negative q
!
        k=1
          Do i = 1,n_dp
            temp1(i)=q(i,k) + dqbydt(i,k) * timestep
            If (temp1(i)  <   qMinInColumn(i)) then
              write(6,*) ' negative q deep',i,temp1(i),dqbydt(i,k)
            End If
          End Do ! n_dp loop

!-----------------------------------------------------------------------
! 8.0  Mixing of the convective increments in the boundary
!      layer.
!-----------------------------------------------------------------------

      If (bl_cnv_mix == 1) then
!      Mixes the increments from the initial parcel perturbation throughout
!      the subcloud layer if this option is selected.

! DEPENDS ON: mix_ipert
        Call MIX_IPERT(n_dp, nlev, nbl, ntml, p_layer_boundaries,       &
                     exner_layer_centres, dthbydt, dqbydt, flx_init,    &
                     thpert, qpert)
      ELSE

!      Mixes convective increments in the boundary
!      layer (essentially distributes incr. at ntml over layers 1 to
!      ntml e.g. incr(1) = incr(2) = incr(ntml)/ntml)
!      Works on boundary layer - columns integrals involved.

! DEPENDS ON: mix_inc
      Call MIX_INC(n_dp,n_dp,n_dp,nlev,nbl,ntml,dthbydt,dqbydt,         &
                   dubydt,dvbydt,L_tracer,ntra,dtrabydt,                &
                   p_layer_boundaries,p_layer_centres,index1)

      End If
!-----------------------------------------------------------------------
! 9.0  3D - Convective cloud amount assumed 3d required ie L_3d_cca
!      is true in old code
!-----------------------------------------------------------------------
! Initialise output array

      Do k = 1,nlev
        Do i = 1,n_dp
          cca(i,k) = 0.0
        End Do
      End Do


!-----------------------------------------------------------------------------
! 9.1 CCRad - Calculate CCA for Deep levels only
!-----------------------------------------------------------------------------

      If (Lcv_ccrad) Then

        !-----------------------------------------------------------------
        ! 9.11 Calculate CCA_2D of Deep Cloud
        !-----------------------------------------------------------------
        Select Case (cca2d_dp_opt)
          Case(cca2d_srf_precip_rate)
            Do i=1, n_dp
              If (iccb(i) /= 0) Then ! Deep convection was successful
        
                ! Use determination of CCA_2D based on surface precip. rate
                ! from deep cloud.
 
                ! NOTE: at present the a_ and b_ parameters for LAND and SEA
                !       are equal, so there will be no difference between
                !       land and sea points.

                If (bland(i)) Then ! LAND POINT
                  If ((rain(i) + snow(i)) > 0.0) Then
                    cca_2d(i)  = a_land                                    &
                                  + b_land                                 &
                                  * ALOG(24.0*3600.0 * (rain(i)+snow(i)))
                  End If
                Else                  ! SEA POINT
                  If ((rain(i) + snow(i)) > 0.0) Then
                    cca_2d(i)  = a_sea                                     &
                                  + b_sea                                  &
                                  * ALOG(24.0*3600.0 * (rain(i)+snow(i)))
                  End If
                End If      ! bland
              End If      ! iccb
            End Do      ! i (n_dp)
          Case(cca2d_total_condensed_water)
            ! cca_2d_dp left unchanged from code, which is based on
            ! TCW (Total Condensed Water) (This is a rate)

        End Select



        Do i=1, n_dp
          !-------------------------------------          
          ! Apply Deep CCA Tuning factor
          !-------------------------------------
          cca_2d(i) = cca_2d(i)*cca_dp_knob

          !-------------------------------------
          ! Make sure cca_2d_dp within limits
          !-------------------------------------

          cca_2d(i) = MAX(0.0,    cca_2d(i))
          cca_2d(i) = MIN(1.0E+0, cca_2d(i))

        End Do      ! i (n_dp)


        !---------------------------------------------------------------------
        ! 1.42 Apply CCA_2D to 3d cloud profile
        !---------------------------------------------------------------------

        If (Lcv_3d_cca) Then

          ! Apply anvil scheme to deep cloud
! DEPENDS ON: CALC_3D_CCA 
          Call calc_3d_cca(n_dp, n_dp, nlev, nbl, iccb, icct            & 
             , p_layer_boundaries, freeze_lev, cca_2d, cca              &
             , z_theta, z_rho, l_q_interact, .false., l_shallow)

          ! NOTE: iccb, icct are layer centres (theta levels) at this
          !        point.

        Else    ! Do we really want this option ?

          ! Apply cca_2d to all levels from deep base to deep top
          Do i=1, n_dp
            Do k=iccb(i), icct(i)
              cca(i,k) = cca_2d(i)
            End Do
          End Do
        End If      ! Lcv_3d_cca

        !---------------------------------------------------------------------
        ! 1.43 Apply ccw_dp_knob
        !---------------------------------------------------------------------
        Do k=1, nlev
          Do i=1, n_dp
            ccw(i,k) = ccw(i,k)*ccw_dp_knob
          End Do
        End Do

      Else        ! Not CCRAD 

      ! Assume 3D anvil cloud required. 

! DEPENDS ON: CALC_3D_CCA 
        Call CALC_3D_CCA(n_dp,n_dp,nlev,nbl,iccb,icct,                  &
                       p_layer_boundaries,freeze_lev,                   &
                       cca_2d,cca,z_theta,z_rho,                        &
                       l_q_interact, .false., l_shallow)

      End If      ! lcv_ccrad
!-----------------------------------------------------------------------
! 10.0  PC2 homogenous forcing ?
!-----------------------------------------------------------------------
! 10.0  End Subroutine
!-----------------------------------------------------------------------

      Return
      END SUBROUTINE DEEP_CONV
#endif
