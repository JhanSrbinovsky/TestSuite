#if defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Congestus convection scheme
!

      SUBROUTINE CONGEST_CONV(nbl,nlev,ntra,n_cca_lev,n_cg,trlev,       &
     &                       bland,delthvu,exner_layer_centres,         &
     &                       exner_layer_boundaries,flg_up_flx,         &
     &                       flg_entr_up, flg_detr_up, flg_dwn_flx,     &
     &                       flg_entr_dwn,flg_detr_dwn,flg_uw_shall,    &
     &                       flg_vw_shall,L_calc_dxek,L_q_interact,     &
     &                       L_tracer,ntml,ntpar,                       &
     &                       pstar,p_layer_centres,                     &
     &                       p_layer_boundaries,z_theta,z_rho,          &
     &                       r2rho_th, dr_across_th,                    &
     &                       q,q1_sd,t1_sd,th,timestep,u,v,uw0,vw0,     &
     &                       wstar,wthvs,zlcl_uv,ztop_uv,freeze_lev,    &
     &                       recip_pstar,qse,                           &
     &                       ad_on,mdet_on,ent_on,sdet_on,new_termc,    &
     &                       bulk_cf,cf_frozen,cf_liquid,qcf,           &
     &                       qcl,tracer,cape_out,cclwp,ccw,cca,         &
     &                       dbcfbydt,dcffbydt,dcflbydt,dqbydt,dqcfbydt,&
     &                       dqclbydt,dthbydt,                          &
     &                       dubydt,dvbydt,dtrabydt,                    &
     &                       detrain_up,detrain_dwn,                    &
     &                       entrain_up,entrain_dwn,                    &
     &                       iccb,icct,lcca,lcclwp,                     &
     &                       lcbase,lctop,rain,snow,up_flux,            &
     &                       dwn_flux,uw_shall,vw_shall,kterm,          &
     &                       tcw,cca_2d,                                &
     &                       rbuoy_p_out,the_out,thp_out,qe_out,qp_out)


!
! Purpose:
! Congestus convection scheme - works on points diagnosed as congestus
! in  subroutine CONV_DIAG. (At present as shallow routine.)
!
!   Called by GLUE_CONV.
!
! Current owners of code: R A Stratton
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!

      Use cv_run_mod, Only:                                             &
          l_mom

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

      integer, intent(in) :: n_cg     ! No. of shallow convection points

      integer, intent(in) :: trlev    ! No. of model levels on which
                                      ! tracers are included

      logical, intent(in) :: bland(n_cg) ! Land/sea mask

      real, intent(in)    :: delthvu(n_cg) !Integral of undilute parcel
                                      ! buoyancy over convective cloud
                                      ! layer (Kelvin m)

      real, intent(in)    :: exner_layer_centres(n_cg,0:nlev) !Exner

      real, intent(in)    :: exner_layer_boundaries(n_cg,0:nlev)
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

      logical, intent(in) :: flg_uw_shall ! STASH flag for x-comp
                                      ! shallow stress

      logical, intent(in) :: flg_vw_shall ! STASH flag for y-comp
                                      ! shallow stress

      logical, intent(in) :: L_calc_dxek ! Switch for calculation of
                                      ! condensate increment

      logical, intent(in) :: L_q_interact ! Switch allows overwriting
                                      ! parcel variables when
                                      ! calculating condensate incr.

      logical, intent(in) :: L_tracer ! Switch for inclusion of tracers

      integer, intent(in) :: ntml(n_cg) ! Top level of surface mixed
                                      ! layer defined relative to
                                      ! theta,q grid

      integer, intent(in) :: ntpar(n_cg) ! Top level of initial parcel
                                      ! ascent in BL scheme defined
                                      ! relative to theta,q grid

      real, intent(in)    :: pstar(n_cg) ! Surface pressure (Pa)

      real, intent(in)    :: p_layer_centres(n_cg,0:nlev) ! Pressure
                                      ! (Pa)


      real, intent(in)    :: p_layer_boundaries(n_cg,0:nlev) ! Pressure
                                      ! at half level above
                                      ! p_layer_centres (Pa)

! Note heights passed in but not currently used - will be
! required by new turbulence based scheme therefore been added to
! arguement list ready for future developments.

      real, intent(in)    ::                                            &
     &  z_theta(n_cg,nlev)                                              &
                                      ! height of theta levels (m)
     &, z_rho(n_cg,nlev)                                                &
                                      ! height of rho levels (m)
     &, r2rho_th(n_cg,nlev)                                             &
                                      ! r2*rho theta levels (kg/m)
     &, dr_across_th(n_cg,nlev)       ! dr across theta levels (m)

      real, intent(in)    :: q(n_cg,nlev) ! Model mixing ratio (kg/kg)

      real, intent(in)    :: q1_sd(n_cg) ! Standard deviation of
                                      ! turbulent flucts. of layer 1 q
                                      ! (kg/kg)

      real, intent(in)    :: t1_sd(n_cg) ! Standard deviation of
                                      ! turbulent flucts. of layer 1
                                      ! temp. (K)

      real, intent(in)    :: th(n_cg,nlev) !Model potential
                                      ! temperature (K)

      real, intent(in)    :: timestep    ! Model timestep (s)

      real, intent(in)    :: u(n_cg,nlev) !Model U field (m/s)

      real, intent(in)    :: v(n_cg,nlev) !Model V field (m/s)

      real, intent(in)    :: uw0(n_cg) ! U-comp of surface stress
                                      ! (N/m2)

      real, intent(in)    :: vw0(n_cg) ! V-comp of surface stress
                                      ! (N/m2)

      real, intent(in)    :: wstar(n_cg) ! Convective velocity scale
                                         ! (m/s)

      real, intent(in)    :: wthvs(n_cg) ! Surface flux of THV (Pa m/s2)


      real, intent(in)    :: zlcl_uv(n_cg) !Lifting condensation level
                                      ! defined for the uv grid (m)

      real, intent(in)    :: ztop_uv(n_cg) ! Top of cloud layer
                                      ! defined for the uv
                                      ! grid (m)

      integer, intent(in) :: freeze_lev(n_cg) ! Level index for freezing
                                               ! level

      real, intent(in) :: recip_pstar(n_cg)  ! Reciprocal of pstar array

      real, intent(in) :: qse(n_cg,nlev) ! Saturation mixing ratio of
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

      real, intent(inout) :: bulk_cf(n_cg,nlev) ! Bulk total cloud
                                      ! volume ( )

      real, intent(inout) :: cf_frozen(n_cg,nlev) ! Frozen water cloud
                                      ! volume ( )

      real, intent(inout) :: cf_liquid(n_cg,nlev) ! Liq water cloud
                                      ! volume ( )

      real, intent(inout) :: qcf(n_cg,nlev) ! Ice condensate mix ratio
                                      ! (kg/kg)

      real, intent(inout) :: qcl(n_cg,nlev) ! Liq condensate mix ratio
                                      ! (kg/kg)

      real, intent(inout)    :: tracer(n_cg,trlev,ntra) !Model tracer
                                      ! fields (kg/kg)

!
! Arguments with intent OUT:
!

      real, intent(out) :: cape_out(n_cg) ! Saved convective available
                                      ! potential energy for diagnostic
                                      ! output (J/kg)
      real, intent(out) :: cclwp(n_cg) ! Condensed water path (k/m2)

      real, intent(out) :: ccw(n_cg,nlev) ! Convective cloud liquid
                                      ! water on model levels (g/kg)
      real, intent(out) :: cca(n_cg,nlev) ! Convective cloud amount
                                      ! on model levels (fraction)
      real, intent(out) :: dbcfbydt(n_cg,nlev) ! Increments to
                                      ! total cld volume due to
                                      ! convection(/s)

      real, intent(out) :: dcffbydt(n_cg,nlev) ! Increments to ice
                                      ! cloud volume due to convection
                                      ! (/s)

      real, intent(out) :: dcflbydt(n_cg,nlev) ! Increments to liq
                                      ! cloud volume due to convection
                                      ! (/s)

      real, intent(out) :: dqbydt(n_cg,nlev) ! Increments to q due to
                                      ! convection (kg/kg/s)

      real, intent(out) :: dqcfbydt(n_cg,nlev) ! Increments to ice
                                      ! condensate due to convection
                                      ! (kg/kg/s)

      real, intent(out) :: dqclbydt(n_cg,nlev) ! Increments to liq
                                      ! condensate due to convection
                                      ! (kg/kg/s)

      real, intent(out) :: dthbydt(n_cg,nlev) ! Increments to potential
                                      ! temp. due to convection (K/s)

      real, intent(out) :: dubydt(n_cg,nlev+1) ! Increments to U due
                                      ! to CMT (m/s2)

      real, intent(out) :: dvbydt(n_cg,nlev+1) ! Increments to V due
                                      ! to CMT (m/s2)

      real, intent(out) :: dtrabydt(n_cg,nlev,ntra) !Increment to tracer
                                      ! due to convection (kg/kg/s)

      real, intent(out) :: detrain_up(n_cg,nlev) ! Fractional
                                      ! detrainment rate into updraughts
                                      ! (Pa/s)

      real, intent(out) :: detrain_dwn(n_cg,nlev) ! Fractional
                                      ! detrainment rate into
                                      ! downdraughts (Pa/s)


      real, intent(out) :: entrain_up(n_cg,nlev) ! Fractional
                                      ! entrainment rate into updraughts
                                      ! (Pa/s)

      real, intent(out) :: entrain_dwn(n_cg,nlev) ! Fractional
                                      ! entrainment rate into
                                      ! downdraughts (Pa/s)

      integer, intent(out) :: iccb(n_cg) ! Convective cloud base
                                      ! level (m)

      integer, intent(out) :: icct(n_cg) ! Convective cloud top
                                      !level (m)

      real, intent(out) :: lcca(n_cg) ! Lowest conv. cloud amt. (%)

      real, intent(out) :: lcclwp(n_cg) ! Condensed water path for
                                      ! lowest conv. cld. level (kg/m2)

      integer, intent(out) :: lcbase(n_cg) ! Lowest conv. cloud base
                                      ! level (m)

      integer, intent(out) :: lctop(n_cg) ! Lowest conv. cloud top
                                      ! level (m)

      real, intent(out) :: rain(n_cg) ! Surface convective rainfall
                                      ! (kg/m2/s)

      real, intent(out) :: snow(n_cg) ! Surface convective snowfall
                                      ! (kg/m2/s)

      real, intent(out) :: up_flux(n_cg,nlev) ! Updraught mass flux
                                      ! (Pa/s)

      real, intent(out) :: dwn_flux(n_cg,nlev) ! Downdraught mass
                                      ! flux (Pa/s)

      real, intent(out) :: uw_shall(n_cg,nlev) ! X-comp. of stress
                                      ! from shallow convection
                                      !(kg/m/s2)

      real, intent(out) :: vw_shall(n_cg,nlev) ! Y-comp. of stress
                                      ! from shallow convection
                                      !(kg/m/s2)
      integer, intent(out) :: kterm(n_cg) ! termination level

      real, intent(out) :: tcw(n_cg)  ! Total condensed water(kg/m2/s)

      real, intent(out) :: cca_2d(n_cg) ! 2D convective cloud amount (%)

!
! Adaptive detrainment output variables
!
      real, intent(out) :: rbuoy_p_out(n_cg,nlev)  !buoyancy excess

      real, intent(out) :: the_out(n_cg,nlev)    !th_E in parcel routine

      real, intent(out) :: thp_out(n_cg,nlev)    !th_P in parcel routine

      real, intent(out) :: qe_out(n_cg,nlev)     !q_E in parcel routine

      real, intent(out) :: qp_out(n_cg,nlev)     !q_P in parcel routine


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

!
! Adaptive detrainment output variables
!
      real :: rbuoy_p_here(n_cg)       !buoyancy excess

      real :: the_here(n_cg)         !th_E in parcel routine

      real :: thp_here(n_cg)         !th_P in parcel routine

      real :: qe_here(n_cg)          !q_E in parcel routine

      real :: qp_here(n_cg)          !q_P in parcel routine

      real :: rbuoy_p_old(n_cg)       !buoyancy excess on previous level

      real :: zk(n_cg)               !heights for use in calc
      real :: zkp12(n_cg)            !of moist static energy
      real :: zkp1(n_cg)


      integer :: index1(n_cg),index2(n_cg)

      integer :: ncposs               ! No. of points which may convect

      integer :: nconv                ! No. of convecting points

      real :: amdetk(n_cg)            ! Mixing detrainment coefficient
                                      ! at level k multiplied by
                                      ! appropriate layer thickness

      real :: b_calc                  ! Coefficient in thpert calc.

      real :: c_calc                  ! Coefficient in thpert calc.

      real :: cape(n_cg)              ! Convective available potential
                                      ! energy (J/kg)

      real :: dq_sat_env              ! dqsat/dT  (kg/kg/K)

      real :: dcpbydt(n_cg)           ! Rate of change of cape (J/kg/s)

      real :: depth(n_cg)             ! Depth of convective cloud (m)

      real :: delexkp1(n_cg)          ! Difference in exner ratio
                                      ! across layer k+1

      real :: dqsthk(n_cg)            ! Gradient of saturation mixing
                                      ! ratio of cloud environment with
                                      ! theta in layer k (kg/kg/K)

      real :: dqsthkp1(n_cg)          ! Gradient of saturation mixing
                                      ! ratio of cloud environment with
                                      ! theta in layer k+1 (kg/kg/K)


      real :: ekp14(n_cg)             ! Entrainment coefficients at
                                      ! level k+1/4 multiplied by
                                      ! appropriate layer thickness
                                      !(dimensionless)

      real :: ekp34(n_cg)             ! Entrainment coefficients at
                                      ! level k+3/4 multiplied by
                                      ! appropriate layer thickness
                                      !(dimensionless)

      real :: ekm14(n_cg)             ! Entrainment coefficients at
                                      ! level k-1+1/4 multiplied by
                                      ! appropriate layer thickness
                                      !(dimensionless)

      real :: exk(n_cg)               ! Exner ratio at layer k

      real :: exkp1(n_cg)             ! Exner ratio at layer k+1

      real :: flxmax(n_cg)            ! Maximum initial convective
                                      ! mass flux (Pa/s)

      real :: flx_init(n_cg)          ! Initial mass flux at cloud base
                                      ! (Pa/s)

      real :: flx_init_new(n_cg)      ! flx_init scaled (Pa/s)

      real :: flxmax_init(n_cg)       ! Maximum possible initial mass
                                      ! flux (limited to the mass in
                                      ! the initial convecting layer
                                      ! in Pa/s)

      real :: max_cfl(n_cg)           ! Max cfl ratio over a convecting
                                      ! layer

      real :: p_lcl(n_cg)             ! Pressure at LCL (Pa)

      real :: precip(n_cg,nlev)       ! Amount of precip from each layer
                                      ! from each layer (kg/m/s)

!      real :: pt(n_cg)                ! Temporary store for P in calc.
                                      ! of sat. mixing ratio. (Pa)

      real :: pk(n_cg)                ! Pressure at midpoint of layer
                                      ! k (Pa)

      real :: pkp1(n_cg)              ! Pressure at midpoint of layer
                                      ! k+1 (Pa)

      real :: delpk(n_cg)             ! Pressure difference over layer
                                      ! k (Pa)

      real :: delpkp1(n_cg)           ! Pressure difference over layer
                                      ! k+1 (Pa)

      real :: delpkp12(n_cg)          ! Pressure difference between
                                      ! layers k and k+1 (Pa)

      real :: delp_uv_k(n_cg)         ! Pressure difference across uv
                                      ! layer k (Pa)

      real :: delp_uv_kp1(n_cg)       ! Pressure difference across uv
                                      ! layer k+1 (Pa)

      real :: q_lcl(n_cg)             ! Mixing ratio at LCL (kg/kg)

      real :: qse_lcl(n_cg)           ! Saturated q at LCL (kg/kg)

      real :: rhum(n_cg)              ! Dummy relative humidity
                                      ! (only used on shallow points)

      real :: t_lcl(n_cg)             ! Temperature at LCL (K)

      real :: th_lcl(n_cg)            ! Theta at LCL (K)

      real :: thv_pert                ! Theta_v parcel pertubation (K)

      real :: thpert                  ! Theta parcel pertubation (K)

!      real :: tt(n_cg)                ! Temporary store for T in calc.
                                      ! of saturation mixing ratio. (K)

!      real :: ttkm1(n_cg)             ! Temporary store for T in layer
                                      ! k-1 for use in freezing level
                                      ! calc. for anvil. (K)

      integer :: start_lev3c(n_cg)    ! Compressed convection
                                      ! initiation level

      real :: wsc(n_cg)               ! Convective velocity scale (m/s)

      real :: wsc_o_mb(n_cg)          ! Convective velocity scale /mb

      logical :: L_shallow(n_cg)      ! Dummy variable (=.T.)

      logical :: L_mid(n_cg)          ! Dummy variable (=.F.)

      logical :: cumulus(n_cg)        ! Dummy variable (=.T.)

      logical :: bgmk(n_cg)           ! Mask for points where parcel in
                                      ! layer k is saturated

      logical :: bwater(n_cg,2:nlev)  ! Mask for points at which
                                      ! condensate is liquid

      logical :: bwk(n_cg)            !mask for liquid condensate on k

      logical :: bwkp1(n_cg)          !mask for liquid condensate on k+1

      logical :: blowst(n_cg)         ! Dummy variable indicating low
                                      ! enough stability for convection
                                      ! to occur

      logical :: bterm(n_cg)          ! Mask for points which have
                                      ! stopped convecting

      logical :: bconv(n_cg)          ! Mask for points at which
                                      ! convection is occurring

      logical :: bcposs(n_cg)         ! Mask for points passing
                                      ! initial stability test

!
! Parcel variables
!

      real :: qpi(n_cg)               ! Initial parcel mixing ratio
                                      !(kg/kg)

      real :: qp(n_cg,nlev)           ! Parcel mixing ratio (kg/kg)

      real :: thpi(n_cg)              ! Initial parcel potential temp.
                                      !(K)

      real :: thp(n_cg,nlev)          ! Parcel potential temp (K)

      real :: up(n_cg,nlev)           ! Parcel U (m/s)

      real :: vp(n_cg,nlev)           ! Parcel V  (m/s)

      real :: trap(n_cg,nlev,ntra)    ! Tracer content of parcel
                                      ! (kg/kg)

      real :: expi(n_cg)              ! Initial parcel exner pressure

      real :: xpk(n_cg,nlev)          ! Parcel cloud water (kg/kg)

      real :: flx(n_cg,nlev)          ! Parcel massflux (Pa/s)

      real :: xsbmin_v(n_cg,nlev)     ! Minmum parcel buoyancy excess

      real :: thpixs_v(n_cg,nlev)     ! Theta parcel excess (K)

      real :: qpixs_v(n_cg,nlev)      ! Q parcel excess(kg/kg)

!
! PC2
!
      real :: qclp(n_cg,nlev)         ! Parcel liquid condensated mixing
                                      ! ratio in layer k (kg/kg)

      real :: qcfp(n_cg,nlev)         ! Parcel frozen condensated mixing
                                      ! ratio in layer k (kg/kg)

!
! Parameters
!
      real, parameter :: CFL_LIMIT = 1.0 ! Max CFL ratio allowed

!
! CMT variables
!

!      integer :: kterm(n_cg)          ! Level index for termination of
                                      ! convection

      integer :: nlcl_uv(n_cg+1)      ! Level index for LCL

      integer :: ntop_uv(n_cg+1)      ! Level index for top of layer

      integer :: n_0degc(n_cg+1)      ! Level index for zero degrees

      integer :: cu_term(n_cg),cu_tend(n_cg) !Indicies for CMT subs

      real :: exk_temp                ! Temporary exner

      real :: eflux_u_ud(n_cg)        ! Vertical eddy flux of momentum
                                      ! due to UD at top of layer
                                      ! (Pa m/s2)

      real :: eflux_v_ud(n_cg)        ! Vertical eddy flux of momentum
                                      ! due to UD at bottom of layer
                                      ! (Pa m/s2)

      real :: flxkp12(nlev,n_cg+1)    ! Mass flux on half level (Pa/s)

      real :: mb(n_cg+1)              ! Cloud base mass flux (Pa/s)

      real :: p_uv(nlev,n_cg+1)       ! Pressure of model level (Pa)

      real :: phalf_uv(nlev,n_cg+1)   ! Pressure of half level (Pa)

      real :: plcl_uv(n_cg+1)         ! Pressure at LCL (Pa)

      real :: ptop_uv(n_cg+1)         ! Pressure at top of cloud layer
                                      ! (Pa)

      real :: p_0degc_uv(n_cg+1)      ! Pressure of zero degree level
                                      ! (Pa)

      real :: rho_uv(nlev,n_cg+1)     ! Density on uv level (kg/m3)

      real :: uw(nlev,n_cg+1)         ! U- comp stress profile (N/m2)
                                      ! (units change through calls)

      real :: ue_p(nlev,n_cg+1)       ! Environment U profile (m/s)

      real :: vw(nlev,n_cg+1)         ! V-comp stress profile (N/m2)

      real :: ve_p(nlev,n_cg+1)       ! Environment V profile (m/s)

      real :: zcld(n_cg)              ! Depth of cloud layer (m)

      logical :: l_mom_GK     ! Set to true if using Gregory-Kershaw CMT scheme

!
! CFL scaling variables
!

      integer :: det_lev(n_cg)        ! Level at which split final
                                      ! detrainment last occurred

      integer :: nterm                ! No. of points where conv.
                                      ! has terminated

      integer :: index_nterm(n_cg)    ! Index for points where conv.
                                      ! has terminated

      real :: tempnum                 ! Temporary variable for storage

      real :: scale_f(n_cg)           ! store scaling factor

!
! Downdraught scheme variables
!

      integer :: nnodd                ! No. of downdraughts not possible

      integer :: index_nodd(n_cg)     ! Index of downdraughts not
                                      ! possible
      integer :: npossdd              ! No. downdraughts possible

      integer :: index_possdd(n_cg)   ! Index of downdraughts

      integer :: kmax_term            ! maximum termination level + 1

      real :: deltap_cld              ! pressure thickness of convective
                                      ! cloud (Pa)

!
! Limit nlev loop to those levels actually required using ntpar
! diagnosed in conv_diag
!
      real :: ntpar_max          ! max ntpar value

!
! parameters ect for qmin checks
!
      real, parameter :: QMIN = 1.0E-8 ! Global minimum allowed Q

      real :: qMinInColumn(n_cg)     ! Minimum value for q in column
                                      ! (kg/kg)
      real :: temp1(n_cg)            ! work array

!
! Local compressed arrays
!

      logical :: bconv_c2(n_cg)

      logical :: bgmkp1_c(n_cg), bgmkp1_c2(n_cg) ! Mask for points
                                      ! where parcel in layer k+1
                                      ! is saturated

      logical :: bwk_c(n_cg), bwk_c2(n_cg) ! bwater mask in layer k

      logical :: bwkp1_c(n_cg), bwkp1_c2(n_cg) ! bwater mask in layer
                                      ! k+1

      real :: deltak_c2(n_cg)         ! Parcel forced detrainment rate
                                      ! in layer k multiplied by
                                      ! appropriate layer thickness

      real :: dqek_c2(n_cg)           ! Increment to q due to
                                      ! convection in layer k (kg/kg)

      real :: dqekp1_c2(n_cg)         ! Increment to q due to
                                      ! convection in layer k+1 (kg/kg)

      real :: dthek_c2(n_cg)          ! Increment to potential temp.
                                      ! due to convection in layer k

      real :: dthekp1_c2(n_cg)        ! Increment to potential temp.
                                      ! due to convection in layer k+1

      real :: dtraek_c2(n_cg,ntra)    ! Increment to model tracer due
                                      ! to conv. at level k (kg/kg/s)

      real :: dtraekp1_c2(n_cg,ntra)  ! Increment to model tracer due
                                      ! to conv. at level k+1 (kg/kg/s)

      real :: duek_c2(n_cg)           ! Increment to model U in layer k
                                      ! due to CMT (m/s2)

      real :: duekp1_c2(n_cg)         ! Increment to model U in layer
                                      ! k+1 due to CMT (m/s2)

      real :: dvek_c2(n_cg)           ! Increment to model V in layer k

      real :: dvekp1_c2(n_cg)         ! Increment to model V in layer
                                      ! k+1 due to CMT (m/s2)

      real :: flxk_c(n_cg), flxk_c2(n_cg) !Parcel mass flux in layer k
                                      ! (Pa/s)

      real :: flxkp12_c2(n_cg)        ! Half level mass flux (Pa/s)

      real :: prekp1_c2(n_cg)         ! Precip. from parcel as it rises
                                      ! from layer k to k+1 (kg/m2/s)

      real :: qpk_c(n_cg), qpk_c2(n_cg) ! Parcel mixing ratio in
                                      ! layer k(kg/kg)

      real :: qpk(n_cg)
      real :: qpkp1_c(n_cg), qpkp1_c2(n_cg) ! Parcel mixing ratio
                                      ! in layer k+1 (kg/kg)

      real :: qek_c(n_cg), qek_c2(n_cg) ! Env. mixing ratio in
                                      ! layer k (kg/kg)

      real :: qek(n_cg)
      real :: qekp1_c(n_cg), qekp1_c2(n_cg) ! Env. mixing ratio in
                                      ! layer k+1 (kgkg-1)

      real :: qekp1(n_cg)
      real :: qsek_c2(n_cg)           ! Saturation mixing ratio of
                                      ! cld. env. in layer k (kg/kg)

      real :: qsek(n_cg)
      real :: qsekp1_c(n_cg), qsekp1_c2(n_cg) ! Saturation mixing ratio
                                      ! of cld. env. in layer k+1
                                      ! (kg/kg)

      real :: qsekp1(n_cg)
      real :: thek_c(n_cg), thek_c2(n_cg) ! Env. potential temp
                                      ! in layer k (K)

      real :: thek(n_cg)
      real :: thekp1_c(n_cg), thekp1_c2(n_cg) ! Env. potential temp i
                                      ! in layer k (K)

      real :: thekp1(n_cg)
      real :: thpk_c(n_cg), thpk_c2(n_cg) ! Parcel potential temp
                                      ! in layer k (K)

      real :: thpk(n_cg)
      real :: thpkp1_c(n_cg), thpkp1_c2(n_cg)! Parcel potential temp
                                      ! in layer k (K)

      real :: thpkp1(n_cg)
      real :: traek_c(n_cg,ntra), traek_c2(n_cg,ntra) ! Tracer content
                                      ! cld. env. in layer k (kgkg-1)

      real :: traekp1_c(n_cg,ntra), traekp1_c2(n_cg,ntra) ! Tracer
                                      ! content of cloud env.
                                      ! in layer k+1 (kg/kg)

      real :: trapk_c(n_cg,ntra), trapk_c2(n_cg,ntra) ! Tracer cont.
                                      ! of parcel in layer k (kg/kg)

      real :: trapkp1_c(n_cg,ntra), trapkp1_c2(n_cg,ntra) ! Tracer cont.
                                      ! of parcel in layer k+1 (kg/kg)

      real :: rbuoy_c(n_cg), rbuoy_c2(n_cg) ! Buoyancy of parcel at k+1
                                      ! (Kelvin)

      real :: uek_c(n_cg), uek_c2(n_cg) ! Model U field on layer k
                                      ! (m/s)

      real :: uekp1_c(n_cg), uekp1_c2(n_cg)! Model U field on layer
                                      ! k+1 (m/s)

      real :: vek_c(n_cg), vek_c2(n_cg) ! Model V field on layer k
                                      ! (m/s)

      real :: vekp1_c(n_cg), vekp1_c2(n_cg) ! Model V field on layer
                                      ! k+1 (m/s)

      real :: upk_c(n_cg), upk_c2(n_cg) ! Parcel U in layer k
                                      ! after entrainment (m/s)


      real :: upkp1_c(n_cg), upkp1_c2(n_cg) ! Parcel U in layer k+1
                                      ! after entrainment (m/s)

      real :: vpk_c(n_cg), vpk_c2(n_cg) ! Parcel V in layer k
                                      ! after entrainment (m/s)

      real :: vpkp1_c(n_cg), vpkp1_c2(n_cg) ! Parcel V in layer k+1
                                      ! after entrainment (m/s)

      real :: xsqkp1_c(n_cg), xsqkp1_c2(n_cg) ! Excess water vapour
                                      ! in parcel at k+1 (kg/kg)

! PC2 compression arrays

      real :: qclek_c(n_cg), qclek_c2(n_cg) ! Environment liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qclekp1_c(n_cg), qclekp1_c2(n_cg) ! Environment liquid
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qcfek_c(n_cg), qcfek_c2(n_cg) ! Environment frozen
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qcfekp1_c(n_cg), qcfekp1_c2(n_cg) ! Environment frozen
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qclpk_c(n_cg), qclpk_c2(n_cg) ! Parcel liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qclpkp1_c(n_cg), qclpkp1_c2(n_cg) ! Parcel liquid
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qcfpk_c(n_cg), qcfpk_c2(n_cg) ! Parcel frozen
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qcfpkp1_c(n_cg), qcfpkp1_c2(n_cg) ! Parcel frozen
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: cflek_c2(n_cg),cflekp1_c2(n_cg) ! Environment liquid water
                                      ! cloud volume ( )

      real :: cffek_c2(n_cg),cffekp1_c2(n_cg) ! Environment frozen water
                                      ! cloud volume ( )

      real :: bcfek_c2(n_cg),bcfekp1_c2(n_cg) ! Environment bulk total
                                      ! cloud volume ( )

      real :: dqclek_c2(n_cg),dqclekp1_c2(n_cg) ! Environment increments
                                      ! to liquid condensate mixing
                                      ! ratio to convection (kg/kg/s)

      real :: dqcfek_c2(n_cg),dqcfekp1_c2(n_cg) ! Environment increments
                                      ! to frozen condensate mixing
                                      ! ratio to convection (kg/kg/s)

      real :: dcflek_c2(n_cg),dcflekp1_c2(n_cg) ! Environment increments
                                      ! to liquid water cloud volume due
                                      ! to convection (/s)

      real :: dcffek_c2(n_cg),dcffekp1_c2(n_cg) ! Environment increments
                                      ! to frozen water cloud volume due
                                      ! to convection (/s)

      real :: dbcfek_c2(n_cg),dbcfekp1_c2(n_cg) ! Environment increments
                                      ! to bulk total cloud volume due
                                      ! to convection (/s)

      real :: amdetk_c2(n_cg)
      logical :: bgmk_c2(n_cg)
      logical :: bland_c2(n_cg)
      logical :: blowst_c2(n_cg)
      logical :: bterm_c2(n_cg)
      real :: cape_c2(n_cg)
      real :: cca_2d_c2(n_cg)
      real :: cclwp_c2(n_cg)
      real :: ccw_c2(n_cg)
      logical :: cumulus_c(n_cg), cumulus_c2(n_cg)
      real :: dcpbydt_c2(n_cg)
      real :: delexkp1_c2(n_cg)
      real :: delpk_c2(n_cg)
      real :: delpkp1_c2(n_cg)
      real :: delp_uv_k_c2(n_cg)
      real :: delp_uv_kp1_c2(n_cg)
      real :: depth_c2(n_cg)
      real :: dptot_c2(n_cg)
      real :: dqsthkp1_c2(n_cg)
      real :: dqsthk_c2(n_cg)
      real :: eflux_u_ud_c2(n_cg)
      real :: eflux_v_ud_c2(n_cg)
      real :: ekp14_c(n_cg),ekp14_c2(n_cg)
      real :: ekp34_c(n_cg),ekp34_c2(n_cg)
      real :: exk_c2(n_cg)
      real :: exkp1_c(n_cg),exkp1_c2(n_cg)
      real :: expi_c2(n_cg)
      integer :: icct_c2(n_cg)
      integer :: iccb_c2(n_cg)
      real :: lcclwp_c2(n_cg)
      integer :: lctop_c2(n_cg)
      integer :: lcbase_c2(n_cg)
      real :: lcca_c2(n_cg)
      logical :: L_shallow_c2(n_cg)
      logical :: L_mid_c2(n_cg)
      real :: max_cfl_c2(n_cg)
      real :: pk_c(n_cg),pk_c2(n_cg)
      real :: pkp1_c(n_cg),pkp1_c2(n_cg)
      real :: pstar_c2(n_cg)
      real :: q1_sd_c2(n_cg)
      real :: qpi_c2(n_cg)
      real :: qpixs_v_c2(n_cg)
      real :: relh_c2(n_cg)
      real :: rbuoy_p_here_c2(n_cg)
      real :: the_here_c2(n_cg)
      real :: thp_here_c2(n_cg)
      real :: qe_here_c2(n_cg)
      real :: qp_here_c2(n_cg)
      real :: rbuoy_p_old_c2(n_cg)
      real :: tcw_c2(n_cg)
      real :: thpi_c2(n_cg)
      real :: thpixs_v_c2(n_cg)
      real :: t1_sd_c2(n_cg)
      real :: xpk_c(n_cg),xpk_c2(n_cg)
      real :: xsbmin_v_c2(n_cg)
      logical :: b_nodd(n_cg)   ! points with no downdraught
      logical :: b_dd(n_cg)     ! points with downdraught on termination

!
! Loop counters
!

      integer :: i,i2,j,k,ktra,kt

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

!
!initialise SCM diagnostics
!
      Do k = 1,nlev
        Do i = 1,n_cg
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

      Do i = 1,n_cg
        blowst(i)    = .true.
        bterm(i)     = .false.
        bconv(i)     = .false.
        bcposs(i)    = .false.
        cumulus(i)   = .true.
        L_shallow(i) = .true.
        L_mid(i)     = .false.
        b_nodd(i)    = .false.
        b_dd(i)      = .false.
      End Do

      l_mom_GK = .false.    ! not using Gregory-Kershaw scheme

!-----------------------------------------------------------------------
! 1.0  Create saturation mixing ratio arrays
!-----------------------------------------------------------------------
! Re-calculate XSBMIN and THPIXS constants based on layer thickness (Pa)
!

      Do k = 1,nlev-1
          Do i = 1,n_cg
            xsbmin_v(i,k) = min( ((p_layer_centres(i,k) -               &
     &              p_layer_centres(i,k+1))/5000.),1.0) *0.2

            thpixs_v(i,k) = min( ((p_layer_centres(i,k) -               &
     &              p_layer_centres(i,k+1))/5000.),1.0)                 &
     &              * THPIXS_SHALLOW

            qpixs_v(i,k)  = QPIXS_SHALLOW
          End Do
      End Do  ! nlev


!
! Calculate convective velocity scale and cloud base mass flux
!

        Do i = 1,n_cg
          wsc(i) = (delthvu(i) * C_MASS * wstar(i) * G / (th(i,ntml(i)) &
     &               * (1.0 + C_VIRTUAL * q(i,ntml(i)))))**0.3333
          mb(i)  = C_MASS * wstar(i)
          zcld(i) = ztop_uv(i) - zlcl_uv(i)
          wsc_o_mb(i) = wsc(i)/mb(i)
        End Do
!
! Define the LCL at the half level above ntml. Find environmental
! T at p_lcl by approximating theta there with
! th(i,k) + constant*(th(i,k+1)-th(i,k))  where constant is tunable.
! Similarly for q.
!

      Do i = 1,n_cg
        k =ntml(i)
        p_lcl(i)  = p_layer_boundaries(i,k)
        th_lcl(i) = th(i,k) + 0.1 * (th(i,k+1) - th(i,k))
        t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / PREF)**KAPPA)
        q_lcl(i)  = q(i,k) + 0.1 * (q(i,k+1) - q(i,k))
      End Do
!
! Calculate saturation mixing ratio at LCL
!

! DEPENDS ON: qsat_mix
      Call QSAT_mix(qse_lcl,t_lcl,p_lcl,n_cg,.false.)

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
!



!
! Calculate theta and q pertubations (pertubation is based on
! environment buoyancy gradient)
! Reset th and q xs's at ntml
!

      Do i = 1,n_cg
          If (t_lcl(i) >  TM) then
            dq_sat_env = EPSILON * LC * qse_lcl(i)                      &
     &                       / (R * t_lcl(i) * t_lcl(i))
          else
            dq_sat_env = EPSILON * (LC+LF) * qse_lcl(i)                 &
     &                       / (R * t_lcl(i) * t_lcl(i))
          End If

          b_calc   = t_lcl(i) * C_VIRTUAL * dq_sat_env + 1.0            &
     &               + C_VIRTUAL * qse_lcl(i)

          thv_pert = -0.17 * wthvs(i) / mb(i)                           &
     &               + (th(i,ntml(i)+1) * (1.0 + C_VIRTUAL              &
     &                 * q(i,ntml(i)+1)) - th(i,ntml(i))                &
     &                   * (1.0 + C_VIRTUAL * q(i,ntml(i))))

          c_calc   = th_lcl(i) * C_VIRTUAL * (qse_lcl(i) - q_lcl(i))    &
     &               - thv_pert

          thpert   = -c_calc / b_calc  !ignore term in THPERT**2

          thpixs_v(i,ntml(i)) = thpert

          qpixs_v(i,ntml(i))  = qse_lcl(i) + ((p_lcl(i) / PREF)         &
     &                         **KAPPA) * thpert * dq_sat_env           &
     &                           - q_lcl(i)

      End Do !n_cg

!
! Set bwater=.true. on points where water will condense rather than
! ice.
! SUBROUTINE FLAG_WET
! UM Documentation paper 27, section (2B)
!

! DEPENDS ON: flag_wet
      Call FLAG_WET(bwater,th,exner_layer_centres,n_cg,n_cg,nlev)

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 2.1  Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
!      increment arrays and cca
!      and dqcl/dt, dqcf/dt, dcfl/dt, dcff/dt, dbcf/dt
!-----------------------------------------------------------------------

      Do k = 1,nlev
        Do i = 1,n_cg
          precip(i,k) = 0.0
          ccw(i,k)    = 0.0
          xpk(i,k)    = 0.0
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
          Do i = 1,n_cg
            dubydt(i,k) = 0.0
            dvbydt(i,k) = 0.0
          End Do
        End Do
      End If  ! L_mom

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,nlev
            Do i = 1,n_cg
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
          Do i = 1,n_cg
            up_flux(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_dwn_flx) then
        Do k = 1,nlev
          Do i = 1,n_cg
            dwn_flux(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_entr_up) then
        Do k = 1,nlev
          Do i = 1,n_cg
            entrain_up(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_detr_up) then
        Do k = 1,nlev
          Do i = 1,n_cg
            detrain_up(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_entr_dwn) then
        Do k = 1,nlev
          Do i = 1,n_cg
           entrain_dwn(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_detr_dwn) then
        Do k = 1,nlev
          Do i = 1,n_cg
            detrain_dwn(i,k) = 0.0
          End Do
        End Do
      End If

      If (L_mom) then
        If (flg_uw_shall) then
          Do k = 1,nlev
            Do i = 1,n_cg
              uw_shall(i,k) = 0.0
            End Do
          End Do
        End If
        If (flg_vw_shall) then
          Do k = 1,nlev
            Do i = 1,n_cg
              vw_shall(i,k) = 0.0
            End Do
          End Do
        End If
      End If  ! L_mom

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------
      Do i = 1,n_cg
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
! 2.4  Initialise gridbox mean diagnostics - done in glue routine
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 2.5  Initialise diagnostics for scaling calculations
!-----------------------------------------------------------------------

      Do i = 1,n_cg
        flx_init(i)     = 0.0
        flx_init_new(i) = 0.0
        cape(i)      = 0.0
        cape_out(i)  = 0.0
        dcpbydt(i)   = 0.0
        max_cfl(i)   = 0.0
        det_lev(i)   = 0

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

!
! set SCM adaptive diagnostics for level k = 1
!
      Do i = 1,n_cg
        rbuoy_p_out(i,1) = 0.0
        the_out(i,1) = th(i,1)
        thp_out(i,1) = th(i,1)
        qe_out(i,1) = q(i,1)
        qp_out(i,1) = q(i,1)
      End Do
! Also, initialise rbuoy_p_old, i.e. for previous level,
! to 0.0 for level1
      Do i = 1,n_cg
        rbuoy_p_old(i) = 0.0
      End Do
!initialise ekm14
      Do i =1, n_cg
         ekm14(i) =0.0
      End do
!Initialise adaptive entrainment variables
!intitilaise to level 2 'cos that's where parcel lift starts from
      Do i = 1, n_cg
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
         Do i = 1, n_cg
           qp(i,k) = q(i,k)
           thp(i,k) = th(i,k)
         End Do
       End Do


!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------
! To reduce cost limit level loop to NTPAR_MAX maximum NTPAR value.
! NTPAR is the top of the parcel ascent for shallow convection.

      ntpar_max=0
      Do i = 1,n_cg
        If (ntpar(i) >  ntpar_max) then
          ntpar_max=ntpar(i)
        End If
      End Do
! Add a level (needed to ensure tracers give same answers on different
! CPU configurations). Not sure why.
      ntpar_max=ntpar_max+1

! This test should not really be required as don't expect any
! shallow convection to reach model top unless a very funny set of
! model levels (very shallow atmosphere).

      If (ntpar_max == nlev) then
        ntpar_max=nlev-1
      End If

!      Do k = 2,nlev-1       ! original level loop

      Do k = 2,ntpar_max

!
! Set relative humidity in layer k (rhum)
!

        Do i = 1,n_cg
          rhum(i) = q(i,k) / qse(i,k)
        End Do

!
!  Initialise adaptive diagnostics for this level
!
      Do i = 1,n_cg
        rbuoy_p_here(i) =0.0
        the_here(i) = th(i,k)
        qe_here(i) = q(i,k)
        rbuoy_p_here_c2(i) =0.0
        the_here_c2(i) = 0.0
        thp_here_c2(i) = 0.0
        qe_here_c2(i) = 0.0
        qp_here_c2(i) = 0.0
      End Do

!Initialise adaptive entrainment variables
      Do i = 1, n_cg
        thek(i)=th(i,k)
        qek(i)=q(i,k)
        qsek(i)=qse(i,k)
        thekp1(i)=th(i,k+1)
        qekp1(i)=q(i,k+1)
        qsekp1(i)=qse(i,k+1)
        If(k  ==  2) then      !set to environmental values
          thpk(i)=th(i,2)
          qpk(i)=q(i,2)
        !only if bconv is true do we have non-zero value for thp(i,k)
        !in this iteration of the main loop - so set to environment
        !value
        Else if (thp(i,k)  <   1) then
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

        zk(i) = z_theta(i,k)
        zkp12(i)=z_rho(i, k+1)
        zkp1(i)=z_theta(i, k+1)
      End Do
!
! Set initial parcel properties for thp, qp if convection
! is not occurring at level k
!
        Do i = 1,n_cg
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
        End Do  ! n_cg


!-----------------------------------------------------------------------
! 3.1  Calculate layer dependent constants (pressure,
!      layer thickness, entrainment coefficients, detrainment
!      coefficients)
!-----------------------------------------------------------------------

! DEPENDS ON: layer_cn
        Call LAYER_CN(k,n_cg,nlev                                       &
      ,               mdet_on, ent_on                                   &
      ,               ntml,ntpar                                        &
      ,               .false.,.true.,.false.                            &
      ,               bconv,bwk,bwkp1                                   &
      ,               exner_layer_boundaries                            &
      ,               exner_layer_centres                               &
      ,               p_layer_boundaries,p_layer_centres                &
      ,               recip_pstar,rhum, zk, zkp12, zkp1                 &
      ,               thek, qek,qsek, thekp1,qekp1,qsekp1               &
      ,               thpk,qpk ,ekm14                                   &
      ,               pkp1,delpkp1,exkp1                                &
      ,               pk,delpk,delpkp12,exk,delexkp1                    &
      ,               delp_uv_k, delp_uv_kp1                            &
      ,               ekp14,ekp34,amdetk)
                         

! Set ekm14 for next pass through loop
         Do i = 1, n_cg
           ekm14(i) = ekp14(i)
         End Do

!
! Calculate dqs/dth for layers k and k+1 (subroutine DQS_DTH)
!

        If (k == 2) then
! DEPENDS ON: dqs_dth
          Call DQS_DTH(dqsthk,k,th(1,k),qse(1,k),exk,n_cg)
        else
          Do i = 1,n_cg
            dqsthk(i) = dqsthkp1(i)
          End Do
        End If

! DEPENDS ON: dqs_dth
        Call DQS_DTH(dqsthkp1,k+1,th(1,k+1),qse(1,k+1),exkp1,n_cg)

!
! Set other grid dependent constants
!

        Do i = 1,n_cg

!
! Maximum initial convective mass flux
!

          flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)

        End Do
!
! Set initial parcel properties (theta,q,tracer,momentum) if convection
! is not occurring at level k
!

        Do i = 1,n_cg
! not convecting and not convected in column before
! PC2 qclp and qcfp zero at this point but will add an initial
! value at cloud base
          If ( .not. bconv(i).and.det_lev(i) == 0) then
            expi(i)  = exk(i)
            xpk(i,k) = 0.0
            qclp(i,k) = 0.0
            qcfp(i,k) = 0.0
            flx(i,k) = 0.0
            bgmk(i)  = .false.
            depth(i) = 0.0
            thpi(i)  = th(i,k) + thpixs_v(i,k)
            thp(i,k) = thpi(i)
            qpi(i)   = q(i,k) + qpixs_v(i,k)
            qp(i,k)  = qpi(i)
            If (L_mom_GK) then
              up(i,k) = u(i,k)
              vp(i,k) = v(i,k)
            End If
          End If
        End Do  ! n_cg
        If (L_tracer) then
          Do ktra=1,ntra
            Do i = 1,n_cg
              If ( .not. bconv(i)) then
                 trap(i,k,ktra)  = tracer(i,k,ktra)
              End If  !not bconv
            enddo
          End Do
        End If

!
! Scale entrainment coefficients with cloud base mass flux
! and convective velocity scale
!
        Do i = 1,n_cg

          If (k  >=  ntml(i)) then
            ekp14(i)  = ekp14(i) * wsc_o_mb(i)
            ekp34(i)  = ekp34(i) * wsc_o_mb(i)
            amdetk(i) = amdetk(i) * wsc_o_mb(i)
          End If

!
! Carry out initial test to see if convection is possible from layer
! k to k+1. Set bcposs = .T. if
! 1. the point was convecting (bconv = .T.) and did not terminate
! in the previous layer  OR
! 2. k = ntml
!

          bcposs(i) = bconv(i) .or. k  ==  ntml(i)

        End Do  ! n_cg

!
! Calculate number of points which may convect (ncposs) and
! set compression indices (index1)
!

        ncposs = 0
        Do i = 1,n_cg
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
            cumulus_c(i) = .true.         ! bgmkp1_c out
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
        Call LIFT_PAR(ncposs,n_cg,thpkp1_c,qpkp1_c,xsqkp1_c,            &
     &           bgmkp1_c,bwkp1_c,bwk_c,thpk_c,qpk_c,xpk_c,thekp1_c,    &
     &           qekp1_c,thek_c,qek_c,qsekp1_c,                         &
     &           qclpkp1_c,qclpk_c,qclekp1_c,qclek_c,l_q_interact,      &
     &           qcfpkp1_c,qcfpk_c,qcfekp1_c,qcfek_c,                   &
     &           pk_c,pkp1_c,exkp1_c,ekp14_c,ekp34_c,L_mom_GK,upkp1_c,  &
     &           vpkp1_c,upk_c,vpk_c,uek_c,uekp1_c,vek_c,vekp1_c,       &
     &           L_tracer,ntra,trapkp1_c,trapk_c,traekp1_c,             &
     &           traek_c,cumulus_c)

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
            bconv(index1(i)) = .true.
!          End If

!
! Set parcel mass flux
! UM Documentation paper 27, section 1.5
!

!          If (bconv(index1(i)).and.k == ntml(index1(i))) then
! as k=ntml
            flxk_c(i) = mb(index1(i)) * G *                             &
     &                      p_layer_centres(index1(i),k) / (R *         &
     &                        thpk_c(i) * (p_layer_centres(index1(i),k) &
     &                          / PREF)**KAPPA)
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

            blowst(index1(i)) = .true.

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

          else
            blowst(index1(i)) = .false.   ! not initial layer

          End If


!
! Reset threashold for forced detrainment to the initial (negative)
! buoyancy
!

          xsbmin_v(index1(i),k) =                                       &
     &                      -0.17 * wthvs(index1(i)) / mb(index1(i))    &
     &               + (th(index1(i),ntml(index1(i))+1) *               &
     &                (1.0 + C_VIRTUAL * q(index1(i),ntml(index1(i))+1))&
     &               - th(index1(i),ntml(index1(i))) * (1.0 + C_VIRTUAL &
     &                    * q(index1(i),ntml(index1(i)))))

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
            blowst_c2(i) = blowst(index1(index2(i)))
            L_shallow_c2(i) = .true.
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
            uek_c2(i)    = uek_c(index2(i))
            uekp1_c2(i)  = uekp1_c(index2(i))
            vek_c2(i)    = vek_c(index2(i))
            vekp1_c2(i)  = vekp1_c(index2(i))
            upkp1_c2(i)    = upkp1_c(index2(i))
            vpkp1_c2(i)    = vpkp1_c(index2(i))
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
            upk_c2(i)    = upk_c(index2(i))
            vpk_c2(i)    = vpk_c(index2(i))
          End Do
          Do i = 1,nconv
            xpk_c2(i)    = xpk_c(index2(i))
            flxk_c2(i)   = flx(index1(index2(i)),k)
            bgmk_c2(i)   = bgmk(index1(index2(i)))
            bterm_c2(i)  = .false.
            dthek_c2(i)  = dthbydt(index1(index2(i)),k)
            dqek_c2(i)   = dqbydt(index1(index2(i)),k)
            dthekp1_c2(i)  = dthbydt(index1(index2(i)),k+1)
            dqekp1_c2(i)   = dqbydt(index1(index2(i)),k+1)
            tcw_c2(i)    = tcw(index1(index2(i)))
            depth_c2(i)  = depth(index1(index2(i)))
            cclwp_c2(i)  = cclwp(index1(index2(i)))
            cape_c2(i)    = cape(index1(index2(i)))
            dcpbydt_c2(i) = dcpbydt(index1(index2(i)))
            relh_c2(i)    = 0.0 ! dummy variable
            dptot_c2(i)   = 0.0 ! dummy variable
            rbuoy_p_here_c2(i)=rbuoy_p_here(index1(index2(i)))
            the_here_c2(i)=the_here(index1(index2(i)))
            thp_here_c2(i)=thp_here(index1(index2(i)))
            qe_here_c2(i)=qe_here(index1(index2(i)))
            qp_here_c2(i)=qp_here(index1(index2(i)))
            rbuoy_p_old_c2(i)=rbuoy_p_old(index1(index2(i)))
            eflux_u_ud_c2(i) = eflux_u_ud(index1(index2(i)))
            eflux_v_ud_c2(i) = eflux_v_ud(index1(index2(i)))
            thpixs_v_c2(i) = thpixs_v(index1(index2(i)),k)
            qpixs_v_c2(i)  = qpixs_v(index1(index2(i)),k)
            xsbmin_v_c2(i) = xsbmin_v(index1(index2(i)),k)
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
          If (L_mom_GK) then
            Do i = 1,nconv
              duek_c2(i)   = dubydt(index1(index2(i)),k)
              dvek_c2(i)   = dvbydt(index1(index2(i)),k)
              duekp1_c2(i)   = dubydt(index1(index2(i)),k+1)
              dvekp1_c2(i)   = dvbydt(index1(index2(i)),k+1)
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

! Part of BTERM comes from values in array before convec2.
! Original code set bterm_c2 according to this test before call to
! convec2.
          Do i = 1,nconv
            If (k  ==  ntpar(index1(index2(i))))then
              bterm_c2(i) = .true.
            End If
          End Do


!-----------------------------------------------------------------------
! 3.3  Calculate the rest of the parcel ascent  and the effect of
!      convection on the large-scale atmosphere.
!
!      Subroutine CONVEC2
!
!      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
!-----------------------------------------------------------------------

! DEPENDS ON: convec2
        Call CONVEC2(nconv,n_cg,nlev,k,thek_c2,thekp1_c2,qek_c2,        &
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
     &               q1_sd_c2,L_mom_GK,uek_c2,uekp1_c2,vek_c2,vekp1_c2, &
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
     &               max_cfl_c2,timestep,                               &
     &               rbuoy_p_here_c2,                                   &
     &               the_here_c2,thp_here_c2,qe_here_c2,qp_here_c2,     &
     &               rbuoy_p_old_c2, ad_on, sdet_on,                    &
     &               new_termc, deltak_c2,                              &
     &               l_calc_dxek,l_q_interact,flxkp12_c2,cumulus_c2,    &
     &               relh_c2,dptot_c2)

!
! Calculate fractional entrainment rate for level k.
! If convection has terminated (bterm=.T.) then set
! fractional entrainment rate for k+1 to zero.
!

      If (flg_entr_up) then
        Do i = 1,nconv
          entrain_up(index1(index2(i)),k) = (1.0 - deltak_c2(i)) *      &
     &                (1.0 - amdetk_c2(i)) * (ekp14_c2(i) + ekp34_c2(i) &
     &                * (1.0 + ekp14_c2(i))) * flx(index1(index2(i)),k)
          If (bterm_c2(i)) then
            entrain_up(index1(index2(i)),k+1) = 0.0
          End If
        End Do
      End If

!
! Calculate fractional detrainment rate for level k
!(and k+1 if bterm=.T.)
!

      If (flg_detr_up) then
        Do i = 1,nconv
          detrain_up(index1(index2(i)),k) = -(amdetk_c2(i)              &
     &                    + deltak_c2(i) * (1.0 - amdetk_c2(i)))        &
     &                  * flx(index1(index2(i)),k)
           If (bterm_c2(i)) then
             detrain_up(index1(index2(i)),k+1) =                        &
     &                 -(1.0 - deltak_c2(i)) * flx(index1(index2(i)),k)
           End If
         End Do
       End If

       End If   ! nconv > 0

!
! Write CONVEC2 compressed output arrays back to full fields
!

        Do i = 1,n_cg
          thp(i,k+1)    = 0.0
          qp(i,k+1)     = 0.0
          xpk(i,k+1)    = 0.0
          flx(i,k+1)    = 0.0
          depth(i)      = 0.0
          precip(i,k+1) = 0.0
          qclp(i,k+1)   = 0.0
          qcfp(i,k+1)   = 0.0
        End Do
        Do i = 1,n_cg
          bgmk(i)       = .false.
          bterm(i)      = .false.
        End Do

        If (L_tracer) then
          Do ktra = 1,ntra
            Do i = 1,n_cg
              trap(i,k+1,ktra) = 0.0
            End Do
          End Do
        End If

        If (L_mom_GK) then
          Do i = 1,n_cg
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
        End Do
        Do i = 1,nconv
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
          max_cfl(index1(index2(i))) =                                  &
     &                   MAX(max_cfl(index1(index2(i))),max_cfl_c2(i))
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
        If (L_mom_GK) then
          Do i = 1,nconv
            dubydt(index1(index2(i)),k) = duek_c2(i)
            dvbydt(index1(index2(i)),k) = dvek_c2(i)
            dubydt(index1(index2(i)),k+1) = duekp1_c2(i)
            dvbydt(index1(index2(i)),k+1) = dvekp1_c2(i)
            eflux_u_ud(index1(index2(i))) = eflux_u_ud_c2(i)
            eflux_v_ud(index1(index2(i))) = eflux_v_ud_c2(i)
          End Do
        End If
        If (L_mom) then      ! require whatever scheme used
          Do i = 1,nconv
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
        End if      ! nconv > 0
!
!   Write adaptive diagnostics for this level to full array for output
!
        Do i = 1,n_cg
          rbuoy_p_out(i,k) = rbuoy_p_here(i)
          the_out(i,k) = the_here(i)
          thp_out(i,k) = thp_here(i)
          qe_out(i,k) = qe_here(i)
          qp_out(i,k) = qp_here(i)
        End Do
!   Write rbuoy for this level to rbuoy_p_old for previous level
!   for use in parcel
!
        Do i = 1,n_cg
          rbuoy_p_old(i) = rbuoy_p_here(i)
        End Do


!-----------------------------------------------------------------------
! 3.4  CFL scaling
!-----------------------------------------------------------------------
!
! Set up integer nterm which is the total number of points where
! convection has terminated.
! Index to full array (n_cg) with index_nterm
!

        nterm = 0
        Do i = 1,n_cg
          If (bterm(i)) then
            nterm = nterm + 1
            index_nterm(nterm) = i
          End If
        End Do

        If (nterm >  0) then

!
! Work out scaled mass flux needed to keep cfl ratio below limit.
! Note L_CAPE not applied to shallow convection

        Do j = 1,nterm
          i = index_nterm(j)

          max_cfl(i) = max_cfl(i) * timestep
          If (max_cfl(i)  >   CFL_LIMIT) then
              flx_init_new(i) = flx_init(i) * CFL_LIMIT                 &
     &                              / max_cfl(i)
          else
              flx_init_new(i) = flx_init(i)
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

        End Do


        Do kt = 2, k+1
          Do j = 1,nterm
            i = index_nterm(j)
            If (kt  >=  ntml(i) .and. flx_init_new(i) >   0.0) then

              dthbydt(i,kt) = dthbydt(i,kt) * scale_f(i)
              dqbydt(i,kt)  = dqbydt(i,kt) * scale_f(i)
              If (L_q_interact) then  ! PC2
                dqclbydt(i,kt)  = dqclbydt(i,kt) * scale_f(i)
                dqcfbydt(i,kt)  = dqcfbydt(i,kt) * scale_f(i)
              End If
              If (L_mom_GK) then
                dubydt(i,kt)  = dubydt(i,kt) * scale_f(i)
                dvbydt(i,kt)  = dvbydt(i,kt) * scale_f(i)
              End If
              If (L_tracer) then
              Do ktra = 1,ntra
               dtrabydt(i,kt,ktra) = dtrabydt(i,kt,ktra)*scale_f(i)
              End Do
              End If

              flx(i,kt)    = flx(i,kt) *  scale_f(i)
              precip(i,kt) = precip(i,kt) *  scale_f(i)

              If (flg_up_flx) then
                up_flux(i,kt) = flx(i,kt)
              End If
              If (flg_entr_up) then
                entrain_up(i,kt) = entrain_up(i,kt) *  scale_f(i)
              End If
              If (flg_detr_up) then
                detrain_up(i,kt) = detrain_up(i,kt) *  scale_f(i)
              End If

            End If !kt >ntml and flx_init_new >0
          End Do  ! j loop
        End Do  ! kt loop

!
! Set final detrainment level (but not used).
!

          Do j = 1,nterm
            i = index_nterm(j)
            det_lev(i)= k+1
          End Do  ! nterm loop

!-----------------------------------------------------------------------
! 3.5  Downdraught calculation - on all points where convection is
!      terminating. Downdraughts are possible for some deeper shallow
!      convection.
!
!      Subroutine DD_CALL
!
!      UM Documentation Paper 27, part 2
!
!-----------------------------------------------------------------------

        npossdd=0
        nnodd = 0

        Do i = 1,nterm
          i2=index_nterm(i)
          tempnum=0.0
          If(iccb(i2) >  0) then
            deltap_cld=p_layer_centres(i2,iccb(i2))                     &
     &                                   -p_layer_centres(i2,k)
            Do kt=iccb(i2),k+1
              tempnum=tempnum+precip(i2,kt)
            End do
          Else
            deltap_cld = 0.0
          End if

! Set logicals for use later

          If (deltap_cld >  15000.0.and.bgmk(i2)                        &
     &                            .and.tempnum >  1e-12) then
             b_dd(i2) = .true.
          Else
             b_nodd(i2) = .true.
          End if
        End Do  ! nterm loop

!
! If convection has terminated write cape to diagnostic output
! variable (cape_out).
! Set kterm array which holds the level index for termination
! of convection.
!
          Do j = 1,nterm
            i=index_nterm(j)
            cape_out(i) = cape(i)
            dcpbydt(i) = 0.0
            cape(i) = 0.0
            kterm(i) = k
            bconv(i) = .false.
          End Do

        End If  ! nterm > 0

!-----------------------------------------------------------------------
! 3.6  End of main loop over levels
!-----------------------------------------------------------------------

      End Do

!-----------------------------------------------------------------------
! 4.0 All shallow convection will terminate at some level. This level
!     has been stored in the main level loop.
!     The convection will either have a down draught or none will be
!     possible.
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
        Do i = 1,n_cg
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
         call DD_ALL_CALL (n_cg,npossdd,kmax_term,nlev,trlev,ntra       &
     &,                      kterm, iccb, icct, index_possdd, l_tracer  &
     &,                      flg_dwn_flx, flg_entr_dwn, flg_detr_dwn    &
     &,                      bwater(1,2)                                &
     &,                      exner_layer_centres,exner_layer_boundaries &
     &,                      p_layer_centres, p_layer_boundaries,pstar  &
     &,                      recip_pstar,timestep , cca_2d              &
     &,                      thp, qp, th, q, trap,tracer, flx,precip    &
     &,                      dthbydt, dqbydt, dtrabydt, rain, snow      &
     &,                      dwn_flux, entrain_dwn, detrain_dwn)

         End if

!-----------------------------------------------------------------------
! 4.2 Surface precipitation calculation for terminating points with
!     no downdraught (moved outside level loop) ie do this calculation
!     on all points at the end.
!-----------------------------------------------------------------------
! Points where no downdraught possible
        nnodd = 0
        Do i = 1,n_cg

          If (b_nodd(i)) then
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
! Only add 1 if kmax_term is less than model levels (which should be
! true).
          If (kmax_term  <  nlev ) then
             kmax_term = kmax_term + 1
          End if

!
! Surface precipitation calculation
!
! DEPENDS ON: evap_bcb_nodd_all
          Call EVAP_BCB_NODD_ALL(n_cg,nnodd,kmax_term,kterm,iccb        &
      ,                      index_nodd,bwater(1,2),exner_layer_centres &
      ,                      exner_layer_boundaries,p_layer_centres     &
      ,                      p_layer_boundaries,pstar,timestep,cca_2d   &
      ,                      th,q,precip,dthbydt,dqbydt,rain,snow)


        End If

!
! Adjust cloud base, top and amount to prevent errors occurring in
! radiation scheme when iccb = icct (this happens when convection
! saturates upon forced detrainment).
!

         Do i = 1,n_cg
          If (iccb(i)  ==  icct(i)) then
            iccb(i)   = 0
            icct(i)   = 0
            cca_2d(i) = 0.0
            tcw(i)    = 0.0
            cclwp(i)  = 0.0
          End If
          If (lcbase(i)  ==  lctop(i)) then
            lcbase(i) = 0
            lctop(i)  = 0
            lcca(i)   = 0.0
            lcclwp(i) = 0.0
          End If
         End Do

!-----------------------------------------------------------------------
! 5.0  Convective Momentum Transport (if L_mom = .T.)
!-----------------------------------------------------------------------

      If (L_mom) then

!
! Initialize arrays required for Convective Momentum Transport(CMT)
!

        k=1
        Do i = 1,n_cg
          p_uv(k,i)     = p_layer_boundaries(i,k-1)
          phalf_uv(k,i) = p_layer_centres(i,k-1)
          ue_p(k,i)     = u(i,k)
          ve_p(k,i)     = v(i,k)
        End Do

        Do i = 1,n_cg
          nlcl_uv(i)    = ntml(i) + 1
          n_0degc(i)    = freeze_lev(i)
        End Do

        Do i = 1,n_cg
          Do k = 2,nlev
            p_uv(k,i)     = p_layer_boundaries(i,k-1)
            phalf_uv(k,i) = p_layer_centres(i,k-1)
            ue_p(k,i)     = u(i,k)
            ve_p(k,i)     = v(i,k)
            exk_temp      = (p_uv(k,i)/PREF)**KAPPA
            rho_uv(k,i)   = 2.0 * p_uv(k,i) / (R * exk_temp *           &
     &                      (th(i,k-1) + th(i,k)))
          End Do
          plcl_uv(i)      = phalf_uv(nlcl_uv(i),i)
          ptop_uv(i)      = phalf_uv(ntop_uv(i),i)
          p_0degc_uv(i)   = phalf_uv(n_0degc(i),i)
          rho_uv(1,i)     = rho_uv(2,i)
        End Do

! altered to use kterm instead of ntpar
        Do i=1,n_cg
          If (kterm(i) >= nlcl_uv(i)) then
            ntop_uv(i) = kterm(i) +1
          Else     ! case where congestus convection fails
! I think in this case the mass flux will be zero so no CMT
            ntop_uv(i) = ntpar(i) + 1
          End If
          ptop_uv(i) = phalf_uv(ntop_uv(i),i)
        End Do             

! Calculate CMT for required points
 
        nterm = 0

        Do i = 1, n_cg
          nterm = nterm + 1
          cu_term(nterm) = i
          cu_tend(nterm) = i
        End Do

! Note using shallow CMT assumptions but may be using top from kterm
! May not be a sensible choice.

        If (nterm  >   0) then

! DEPENDS ON: shallow_grad_stress
          Call SHALLOW_GRAD_STRESS(n_cg,n_cg,nterm,nlev,cu_term,        &
     &                             nlcl_uv,ntop_uv,mb,wsc,wstar,zcld,   &
     &                             plcl_uv,ptop_uv,p_uv,phalf_uv,       &
     &                             rho_uv,ue_p,ve_p,timestep,           &
                                   ! IN
     &                             uw,vw)

! DEPENDS ON: shallow_base_stress
          Call SHALLOW_BASE_STRESS(n_cg,n_cg,n_cg,nlev,nterm,cu_term,   &
     &                             cu_tend,nlcl_uv,ntop_uv,mb,wsc,      &
     &                             zlcl_uv,zcld,uw0,vw0,plcl_uv,        &
     &                             ptop_uv,ue_p,ve_p,phalf_uv,p_uv,     &
     &                             rho_uv,timestep,flg_uw_shall,        &
     &                             flg_vw_shall,                        &
                                   ! INOUT
     &                             uw,vw,                               &
                                   ! OUT
     &                             uw_shall,vw_shall)


! DEPENDS ON: shallow_cmt_incr
          Call SHALLOW_CMT_INCR(n_cg,n_cg,n_cg,nlev,nterm,cu_term,      &
     &                          cu_tend,nlcl_uv,ntop_uv,uw,vw,phalf_uv, &
     &                          rho_uv,zlcl_uv,                         &
                                !OUT
     &                          dubydt,dvbydt)

        End If  ! nterm>0
      End If ! L_mom

!-----------------------------------------------------------------------
! 6.0  Energy correction calculation - removed as old code not correct
!     for new dynamics grid ( attempts to correct this give problems).
!     UM documentation paper 27 - section 12.
!-----------------------------------------------------------------------
      Do i = 1,n_cg
          index1(i) = i
      End Do

!-----------------------------------------------------------------------
! 7.0  Total water conservation  - also works on whole column
!-----------------------------------------------------------------------

! only check columns where convection has occurred.


        Do i = 1,n_cg
          qMinInColumn(i) = q(i,nlev)
        End Do
        Do k = 1,nlev-1
          Do i = 1,n_cg
            If (q(i,k)  <   qMinInColumn(i)) then
              qMinInColumn(i) = q(i,k)
            End If
          End Do
        End Do

!
! Ensure Q does not go below global allowed minimum (QMIN)
!
        Do i = 1,n_cg
          qMinInColumn(i)=MAX(QMIN,qMinInColumn(i))
        End Do

!
! Apply an artificial upwards flux from k-1 level to ensure Q
! remians above minimum value in the column.
!

        Do k = nlev,2,-1
          Do i = 1,n_cg
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
          End Do ! n_s loop
        End Do  ! nlev
!
! check negative q
!
        k=1
          Do i = 1,n_cg
            temp1(i)=q(i,k) + dqbydt(i,k) * timestep
            If (temp1(i)  <   qMinInColumn(i)) then
              write(6,*) ' negative q deep',i,temp1(i),dqbydt(i,k)
            End If
          End Do ! n_cg loop

!-----------------------------------------------------------------------
! 8.0  Subroutine MIX_INC mixes convective increments in the boundary
!      layer (essentially distributes incr. at ntml over layers 1 to
!      ntml e.g. incr(1) = incr(2) = incr(ntml)/ntml)
!      Works on boundary layer - columns integrals involved.
!-----------------------------------------------------------------------

! DEPENDS ON: mix_inc
      Call MIX_INC(n_cg,n_cg,n_cg,nlev,nbl,ntml,                        &
                   dthbydt,dqbydt,dubydt,dvbydt,L_tracer,ntra,dtrabydt, &
                   p_layer_boundaries,p_layer_centres,index1)


!-----------------------------------------------------------------------
! 9.0  Calculate convective cloud amount on model levels - no anvils 
!-----------------------------------------------------------------------
! Initialise output array

      Do k = 1,n_cca_lev
        Do i = 1,n_cg
          cca(i,k)        = 0.0
        End Do
      End Do


! Note if (l_q_interact) PC2 then remains set to zero

      If (.not.l_q_interact) Then
        Do k = 1,n_cca_lev
          Do i = 1,n_cg
            If (k >= iccb(i) .and. k <= icct(i)) then  
              cca(i,k) = cca_2d(i)
            End If
          End Do    
        End Do  
      End If


!-----------------------------------------------------------------------
! 10.0  End Subroutine
!-----------------------------------------------------------------------

      Return
      END SUBROUTINE CONGEST_CONV
#endif
