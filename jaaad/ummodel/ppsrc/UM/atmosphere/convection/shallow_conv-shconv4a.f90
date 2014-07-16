
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Shallow convection scheme
!

      SUBROUTINE SHALLOW_CONV(nbl,nlev,ntra,n_cca_lev,n_sh,trlev,       &
     &                       bland,delthvu,exner_layer_centres,         &
     &                       exner_layer_boundaries,flg_up_flx,         &
     &                       flg_up_flx_half,                           &
     &                       flg_entr_up, flg_detr_up, flg_dwn_flx,     &
     &                       flg_entr_dwn,flg_detr_dwn,flg_uw_shall,    &
     &                       flg_vw_shall, L_calc_dxek,                 &
     &                       L_q_interact, L_tracer, ntml, ntpar,       &
     &                       pstar,p_layer_centres,                     &
     &                       p_layer_boundaries,z_theta,z_rho,          &
     &                       r2rho_th, dr_across_th,                    &
     &                       q,q1_sd,t1_sd,th,timestep,u,v,uw0,vw0,     &
     &                       wstar,wthvs,zlcl_uv,                       &
     &                       ztop_uv,freeze_lev,recip_pstar,qse,ad_on,  &
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
     &                       up_flux, up_flux_half,                     &
     &                       dwn_flux,uw_shall,vw_shall,tcw,cca_2d,     &
     &                       rbuoy_p_out,the_out,thp_out,qe_out,qp_out  &
     &                       )


!
! Purpose:
!   Shallow convection scheme - works on points diagnosed as shallow in
!   subroutine CONV_DIAG.
!
!   Called by GLUE_CONV.
!
! Current owners of code: R A Stratton
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!

      Use cv_cntl_mod, Only:                                            &
          lcv_ccrad

      Use cv_run_mod, Only:                                             &
          l_mom, sh_pert_opt, bl_cnv_mix
              

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

      integer, intent(in) :: n_sh     ! No. of shallow convection points

      integer, intent(in) :: trlev    ! No. of model levels on which
                                      ! tracers are included

      logical, intent(in) :: bland(n_sh) ! Land/sea mask

      real, intent(in)    :: delthvu(n_sh) !Integral of undilute parcel
                                      ! buoyancy over convective cloud
                                      ! layer (Kelvin m)

      real, intent(in)    :: exner_layer_centres(n_sh,0:nlev) !Exner

      real, intent(in)    :: exner_layer_boundaries(n_sh,0:nlev)
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

      integer, intent(in) :: ntml(n_sh) ! Top level of surface mixed
                                      ! layer defined relative to
                                      ! theta,q grid

      integer, intent(in) :: ntpar(n_sh) ! Top level of initial parcel
                                      ! ascent in BL scheme defined
                                      ! relative to theta,q grid

      real, intent(in)    :: pstar(n_sh) ! Surface pressure (Pa)

      real, intent(in)    :: p_layer_centres(n_sh,0:nlev) ! Pressure
                                      ! (Pa)


      real, intent(in)    :: p_layer_boundaries(n_sh,0:nlev) ! Pressure
                                      ! at half level above
                                      ! p_layer_centres (Pa)

! Note heights passed in but not currently used - will be
! required by new turbulence based scheme therefore been added to
! arguement list ready for future developments.

      real, intent(in)    :: z_theta(n_sh,nlev) ! height of theta levels
                                       ! (m)
      real, intent(in)    :: z_rho(n_sh,nlev)                           &
                                              ! height of rho levels (m)
     &, r2rho_th(n_sh,nlev)                                             &
                                 ! r**2*rho on theta levels (kg/m)
     &, dr_across_th(n_sh,nlev)  ! thickness of theta levels (m)

      real, intent(in)    :: q(n_sh,nlev) ! Model mixing ratio (kg/kg)

      real, intent(in)    :: q1_sd(n_sh) ! Standard deviation of
                                      ! turbulent flucts. of layer 1 q
                                      ! (kg/kg)

      real, intent(in)    :: t1_sd(n_sh) ! Standard deviation of
                                      ! turbulent flucts. of layer 1
                                      ! temp. (K)

      real, intent(in)    :: th(n_sh,nlev) !Model potential
                                      ! temperature (K)

      real, intent(in)    :: timestep    ! Model timestep (s)

      real, intent(in)    :: u(n_sh,nlev) !Model U field (m/s)

      real, intent(in)    :: v(n_sh,nlev) !Model V field (m/s)

      real, intent(in)    :: uw0(n_sh) ! U-comp of surface stress
                                      ! (N/m2)

      real, intent(in)    :: vw0(n_sh) ! V-comp of surface stress
                                      ! (N/m2)

      real, intent(in)    :: wstar(n_sh) ! Convective velocity scale
                                         ! (m/s)

      real, intent(in)    :: wthvs(n_sh) ! Surface flux of THV (Pa m/s2)


      real, intent(in)    :: zlcl_uv(n_sh) !Lifting condensation level
                                      ! defined for the uv grid (m)

      real, intent(in)    :: ztop_uv(n_sh) ! Top of cloud layer
                                      ! defined for the uv
                                      ! grid (m)

      integer, intent(in) :: freeze_lev(n_sh) ! Level index for freezing
                                               ! level

      real, intent(in) :: recip_pstar(n_sh)  ! Reciprocal of pstar array

      real, intent(in) :: qse(n_sh,nlev) ! Saturation mixing ratio of
                                      ! cloud environment (kg/kg)

      integer, intent(in) :: ad_on   !flag for adaptive detrainment

      integer, intent(in) :: mdet_on !flag for adaptive detrainment

      integer, intent(in) :: ent_on  !Flag for adaptive entrainment

      integer, intent(in) :: sdet_on    !flag for smoothed forced detrainment
                                        !0 = off, 1 = on

      integer, intent(in) :: new_termc !flag for simplified termination
                                       !of convection

!
! Arguments with intent INOUT:
!

      real, intent(inout) :: bulk_cf(n_sh,nlev) ! Bulk total cloud
                                      ! volume ( )

      real, intent(inout) :: cf_frozen(n_sh,nlev) ! Frozen water cloud
                                      ! volume ( )

      real, intent(inout) :: cf_liquid(n_sh,nlev) ! Liq water cloud
                                      ! volume ( )

      real, intent(inout) :: qcf(n_sh,nlev) ! Ice condensate mix ratio
                                      ! (kg/kg)

      real, intent(inout) :: qcl(n_sh,nlev) ! Liq condensate mix ratio
                                      ! (kg/kg)

      real, intent(inout)    :: tracer(n_sh,trlev,ntra) !Model tracer
                                      ! fields (kg/kg)

!
! Arguments with intent OUT:
!

      real, intent(out) :: cape_out(n_sh) ! Saved convective available
                                      ! potential energy for diagnostic
                                      ! output (J/kg)
      real, intent(out) :: cclwp(n_sh) ! Condensed water path (k/m2)

      real, intent(out) :: ccw(n_sh,nlev) ! Convective cloud liquid

                                      ! water on model levels (g/kg)
      real, intent(out) :: dbcfbydt(n_sh,nlev) ! Increments to
                                      ! total cld volume due to
                                      ! convection(/s)

      real, intent(out) :: dcffbydt(n_sh,nlev) ! Increments to ice
                                      ! cloud volume due to convection
                                      ! (/s)

      real, intent(out) :: dcflbydt(n_sh,nlev) ! Increments to liq
                                      ! cloud volume due to convection
                                      ! (/s)

      real, intent(out) :: dqbydt(n_sh,nlev) ! Increments to q due to
                                      ! convection (kg/kg/s)

      real, intent(out) :: dqcfbydt(n_sh,nlev) ! Increments to ice
                                      ! condensate due to convection
                                      ! (kg/kg/s)

      real, intent(out) :: dqclbydt(n_sh,nlev) ! Increments to liq
                                      ! condensate due to convection
                                      ! (kg/kg/s)

      real, intent(out) :: dthbydt(n_sh,nlev) ! Increments to potential
                                      ! temp. due to convection (K/s)

      real, intent(out) :: dubydt(n_sh,nlev+1) ! Increments to U due
                                      ! to CMT (m/s2)

      real, intent(out) :: dvbydt(n_sh,nlev+1) ! Increments to V due
                                      ! to CMT (m/s2)

      real, intent(out) :: dtrabydt(n_sh,nlev,ntra) !Increment to tracer
                                      ! due to convection (kg/kg/s)

      real, intent(out) :: detrain_up(n_sh,nlev) ! Fractional
                                      ! detrainment rate into updraughts
                                      ! (Pa/s)

      real, intent(out) :: detrain_dwn(n_sh,nlev) ! Fractional
                                      ! detrainment rate into
                                      ! downdraughts (Pa/s)


      real, intent(out) :: entrain_up(n_sh,nlev) ! Fractional
                                      ! entrainment rate into updraughts
                                      ! (Pa/s)

      real, intent(out) :: entrain_dwn(n_sh,nlev) ! Fractional
                                      ! entrainment rate into
                                      ! downdraughts (Pa/s)

      integer, intent(out) :: iccb(n_sh) ! Convective cloud base
                                      ! level (m)

      integer, intent(out) :: icct(n_sh) ! Convective cloud top
                                      !level (m)

      real, intent(out) :: lcca(n_sh) ! Lowest conv. cloud amt. (%)

      real, intent(out) :: lcclwp(n_sh) ! Condensed water path for
                                      ! lowest conv. cld. level (kg/m2)

      integer, intent(out) :: lcbase(n_sh) ! Lowest conv. cloud base
                                      ! level (m)

      integer, intent(out) :: lctop(n_sh) ! Lowest conv. cloud top
                                      ! level (m)

      real, intent(out) :: rain(n_sh) ! Surface convective rainfall
                                      ! (kg/m2/s)

      real, intent(out) :: snow(n_sh) ! Surface convective snowfall
                                      ! (kg/m2/s)
      real, intent(out) :: rain_3d(n_sh,nlev) ! Convective rainfall flux
                                      ! (kg/m2/s)

      real, intent(out) :: snow_3d(n_sh,nlev) ! Convective snowfall flux
                                      ! (kg/m2/s)



      real, intent(out) :: up_flux(n_sh,nlev) ! Updraught mass flux
                                      ! (Pa/s)

      real, intent(out) :: up_flux_half(n_sh,nlev)
                                      ! Updraught mass flux
                                      ! (Pa/s) on rho levels

      real, intent(out) :: dwn_flux(n_sh,nlev) ! Downdraught mass
                                      ! flux (Pa/s)

      real, intent(out) :: uw_shall(n_sh,nlev) ! X-comp. of stress
                                      ! from shallow convection
                                      !(kg/m/s2)

      real, intent(out) :: vw_shall(n_sh,nlev) ! Y-comp. of stress
                                      ! from shallow convection
                                      !(kg/m/s2)
      real, intent(out) :: tcw(n_sh)  ! Total condensed water(kg/m2/s)

      real, intent(out) :: cca_2d(n_sh) ! 2D convective cloud amount (%)


!
! Adaptive detrainment output variables
!
      real, intent(out) :: rbuoy_p_out(n_sh,nlev)  !buoyancy excess

      real, intent(out) :: the_out(n_sh,nlev)    !th_E in parcel
                                                 !routine

      real, intent(out) :: thp_out(n_sh,nlev)    !th_P in parcel
                                                 !routine

      real, intent(out) :: qe_out(n_sh,nlev)     !q_E in parcel routine

      real, intent(out) :: qp_out(n_sh,nlev)     !q_P in parcel routine


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
!
! Adaptive detrainment output variables
!
      real :: rbuoy_p_here(n_sh)       !buoyancy excess

      real :: the_here(n_sh)         !th_E in parcel routine

      real :: thp_here(n_sh)         !th_P in parcel routine

      real :: qe_here(n_sh)          !q_E in parcel routine

      real :: qp_here(n_sh)          !q_P in parcel routine

      real :: rbuoy_p_old(n_sh)       !buoyancy excess from previous k

      real :: zk(n_sh)               !heights for use in calc

      real :: zkp12(n_sh)            !of moist static energy

      real :: zkp1(n_sh)

      integer :: index1(n_sh),index2(n_sh)

      integer :: ncposs               ! No. of points which may convect

      integer :: nconv                ! No. of convecting points

      real :: amdetk(n_sh)            ! Mixing detrainment coefficient
                                      ! at level k multiplied by
                                      ! appropriate layer thickness

      real :: b_calc                  ! Coefficient in thpert calc.

      real :: c_calc                  ! Coefficient in thpert calc.

      real :: cape(n_sh)              ! Convective available potential
                                      ! energy (J/kg)

      real :: dq_sat_env              ! dqsat/dT  (kg/kg/K)

      real :: dcpbydt(n_sh)           ! Rate of change of cape (J/kg/s)

      real :: depth(n_sh)             ! Depth of convective cloud (m)

      real :: delexkp1(n_sh)          ! Difference in exner ratio
                                      ! across layer k+1

      real :: dqsthk(n_sh)            ! Gradient of saturation mixing
                                      ! ratio of cloud environment with
                                      ! theta in layer k (kg/kg/K)

      real :: dqsthkp1(n_sh)          ! Gradient of saturation mixing
                                      ! ratio of cloud environment with
                                      ! theta in layer k+1 (kg/kg/K)


      real :: ekp14(n_sh)             ! Entrainment coefficients at
                                      ! level k+1/4 multiplied by
                                      ! appropriate layer thickness
                                      !(dimensionless)

      real :: ekp34(n_sh)             ! Entrainment coefficients at
                                      ! level k+3/4 multiplied by
                                      ! appropriate layer thickness
                                      !(dimensionless)

      real :: ekm14(n_sh)             ! Entrainment coefficients at
                                      ! level k-1+1/4 multiplied by
                                      ! appropriate layer thickness
                                      !(dimensionless)

      real :: exk(n_sh)               ! Exner ratio at layer k

      real :: exkp1(n_sh)             ! Exner ratio at layer k+1

      real :: flxmax(n_sh)            ! Maximum initial convective
                                      ! mass flux (Pa/s)

      real :: flx_init(n_sh)          ! Initial mass flux at cloud base
                                      ! (Pa/s)

      real :: flx_init_new(n_sh)      ! flx_init scaled (Pa/s)

      real :: flxmax_init(n_sh)       ! Maximum possible initial mass
                                      ! flux (limited to the mass in
                                      ! the initial convecting layer
                                      ! in Pa/s)

      real :: max_cfl(n_sh)           ! Max cfl ratio over a convecting
                                      ! layer

      real :: p_lcl(n_sh)             ! Pressure at LCL (Pa)

      real :: precip(n_sh,nlev)       ! Amount of precip from each layer
                                      ! from each layer (kg/m/s)

!      real :: pt(n_sh)                ! Temporary store for P in calc.
                                      ! of sat. mixing ratio. (Pa)

      real :: pk(n_sh)                ! Pressure at midpoint of layer
                                      ! k (Pa)

      real :: pkp1(n_sh)              ! Pressure at midpoint of layer
                                      ! k+1 (Pa)

      real :: delpk(n_sh)             ! Pressure difference over layer
                                      ! k (Pa)

      real :: delpkp1(n_sh)           ! Pressure difference over layer
                                      ! k+1 (Pa)

      real :: delpkp12(n_sh)          ! Pressure difference between
                                      ! layers k and k+1 (Pa)

      real :: delp_uv_k(n_sh)         ! Pressure difference across uv
                                      ! layer k (Pa)

      real :: delp_uv_kp1(n_sh)       ! Pressure difference across uv
                                      ! layer k+1 (Pa)

      real :: q_lcl(n_sh)             ! Mixing ratio at LCL (kg/kg)

      real :: qse_lcl(n_sh)           ! Saturated q at LCL (kg/kg)

      real :: rhum(n_sh)              ! Dummy relative humidity
                                      ! (only used on shallow points)

      real :: t_lcl(n_sh)             ! Temperature at LCL (K)

      real :: th_lcl(n_sh)            ! Theta at LCL (K)

      real :: dthv_ma                 ! Moist adiabtic change in thv
                                      ! from ntml to ntml+1 (K)

      real :: thv_pert(n_sh)          ! Theta_v parcel pertubation (K)

      real :: thpert(n_sh)            ! Theta parcel pertubation (K)
      
      real :: qpert(n_sh)             ! q parcel pertubation (kg/kg)

      real :: rho_k                   ! density on level k

!      real :: tt(n_sh)                ! Temporary store for T in calc.
                                      ! of saturation mixing ratio. (K)

!      real :: ttkm1(n_sh)             ! Temporary store for T in layer
                                      ! k-1 for use in freezing level
                                      ! calc. for anvil. (K)

      integer :: start_lev3c(n_sh)    ! Compressed convection
                                      ! initiation level

      real :: wsc(n_sh)               ! Convective velocity scale (m/s)

      real :: wsc_o_mb(n_sh)          ! Convective velocity scale /mb

      logical :: L_shallow(n_sh)      ! Dummy variable (=.T.)

      logical :: L_mid(n_sh)          ! Dummy variable (=.F.)

      logical :: cumulus(n_sh)        ! Dummy variable (=.T.)

      logical :: bgmk(n_sh)           ! Mask for points where parcel in
                                      ! layer k is saturated

      logical :: bwater(n_sh,2:nlev)  ! Mask for points at which
                                      ! condensate is liquid
      logical :: bwk(n_sh)            !mask for liquid condensate on k

      logical :: bwkp1(n_sh)          !mask for liquid condensate on
                                      !k+1

      logical :: blowst(n_sh)         ! Dummy variable indicating low
                                      ! enough stability for convection
                                      ! to occur

      logical :: bterm(n_sh)          ! Mask for points which have
                                      ! stopped convecting

      logical :: bconv(n_sh)          ! Mask for points at which
                                      ! convection is occurring

      logical :: bcposs(n_sh)         ! Mask for points passing
                                      ! initial stability test

!
! Parcel variables
!

      real :: qpi(n_sh)               ! Initial parcel mixing ratio
                                      !(kg/kg)

      real :: qp(n_sh,nlev)           ! Parcel mixing ratio (kg/kg)

      real :: thpi(n_sh)              ! Initial parcel potential temp.
                                      !(K)

      real :: thp(n_sh,nlev)          ! Parcel potential temp (K)

      real :: up(n_sh,nlev)           ! Parcel U (m/s)

      real :: vp(n_sh,nlev)           ! Parcel V  (m/s)

      real :: trap(n_sh,nlev,ntra)    ! Tracer content of parcel
                                      ! (kg/kg)

      real :: expi(n_sh)              ! Initial parcel exner pressure

      real :: xpk(n_sh,nlev)          ! Parcel cloud water (kg/kg)

      real :: flx(n_sh,nlev)          ! Parcel massflux (Pa/s)

      real :: xsbmin_v(n_sh,nlev)     ! Minmum parcel buoyancy excess

      real :: thpixs_v(n_sh,nlev)     ! Theta parcel excess (K)

      real :: qpixs_v(n_sh,nlev)      ! Q parcel excess(kg/kg)

!
! PC2
!
      real :: qclp(n_sh,nlev)         ! Parcel liquid condensated mixing
                                      ! ratio in layer k (kg/kg)

      real :: qcfp(n_sh,nlev)         ! Parcel frozen condensated mixing
                                      ! ratio in layer k (kg/kg)

!
! Parameters
!
      real, parameter :: CFL_LIMIT = 1.0 ! Max CFL ratio allowed

!
! CMT variables
!

      integer :: kterm(n_sh)          ! Level index for termination of
                                      ! convection

      integer :: nlcl_uv(n_sh+1)      ! Level index for LCL

      integer :: ntop_uv(n_sh+1)      ! Level index for top of layer

      integer :: n_0degc(n_sh+1)      ! Level index for zero degrees

      integer :: cu_term(n_sh),cu_tend(n_sh) !Indicies for CMT subs

      real :: exk_temp                ! Temporary exner

      real :: eflux_u_ud(n_sh)        ! Vertical eddy flux of momentum
                                      ! due to UD at top of layer
                                      ! (Pa m/s2)

      real :: eflux_v_ud(n_sh)        ! Vertical eddy flux of momentum
                                      ! due to UD at bottom of layer
                                      ! (Pa m/s2)

      real :: flxkp12(nlev,n_sh+1)    ! Mass flux on half level (Pa/s)

      real :: mb(n_sh+1)              ! Cloud base mass flux (Pa/s)

      real :: p_uv(nlev,n_sh+1)       ! Pressure of model level (Pa)

      real :: phalf_uv(nlev,n_sh+1)   ! Pressure of half level (Pa)

      real :: plcl_uv(n_sh+1)         ! Pressure at LCL (Pa)

      real :: ptop_uv(n_sh+1)         ! Pressure at top of cloud layer
                                      ! (Pa)

      real :: p_0degc_uv(n_sh+1)      ! Pressure of zero degree level
                                      ! (Pa)

      real :: rho_uv(nlev,n_sh+1)     ! Density on uv level (kg/m3)

      real :: uw(nlev,n_sh+1)         ! U- comp stress profile (N/m2)
                                      ! (units change through calls)

      real :: ue_p(nlev,n_sh+1)       ! Environment U profile (m/s)

      real :: vw(nlev,n_sh+1)         ! V-comp stress profile (N/m2)

      real :: ve_p(nlev,n_sh+1)       ! Environment V profile (m/s)

      real :: zcld(n_sh)              ! Depth of cloud layer (m)

      logical :: l_mom_gk             ! true if Gregory-Kershaw CMT

!
! CFL scaling variables
!

      integer :: det_lev(n_sh)        ! Level at which split final
                                      ! detrainment last occurred

      integer :: nterm                ! No. of points where conv.
                                      ! has terminated

      integer :: index_nterm(n_sh)    ! Index for points where conv.
                                      ! has terminated

      real :: tempnum                 ! Temporary variable for storage

      real :: scale_f(n_sh)           ! store scaling factor

!
! Downdraught scheme variables
!

      integer :: nnodd                ! No. of downdraughts not possible

      integer :: index_nodd(n_sh)     ! Index of downdraughts not
                                      ! possible
      integer :: npossdd              ! No. downdraughts possible

      integer :: index_possdd(n_sh)   ! Index of downdraughts

      integer :: kmax_term            ! maximum termination level + 1

      real :: deltap_cld              ! pressure thickness of convective
                                      ! cloud (Pa)

!
! Limit nlev loop to those levels actually required using ntpar
! diagnosed in conv_diag
!
      real :: ntpar_max          ! max ntpar value

!
! Local compressed arrays
!

      logical :: bconv_c2(n_sh)

      logical :: bgmkp1_c(n_sh), bgmkp1_c2(n_sh) ! Mask for points
                                      ! where parcel in layer k+1
                                      ! is saturated

      logical :: bwk_c(n_sh), bwk_c2(n_sh) ! bwater mask in layer k

      logical :: bwkp1_c(n_sh), bwkp1_c2(n_sh) ! bwater mask in layer
                                      ! k+1

      real :: deltak_c2(n_sh)         ! Parcel forced detrainment rate
                                      ! in layer k multiplied by
                                      ! appropriate layer thickness
                                      
      real :: dqek_c2(n_sh)           ! Increment to q due to
                                      ! convection in layer k (kg/kg)

      real :: dqekp1_c2(n_sh)         ! Increment to q due to
                                      ! convection in layer k+1 (kg/kg)

      real :: dthek_c2(n_sh)          ! Increment to potential temp.
                                      ! due to convection in layer k

      real :: dthekp1_c2(n_sh)        ! Increment to potential temp.
                                      ! due to convection in layer k+1

      real :: dtraek_c2(n_sh,ntra)    ! Increment to model tracer due
                                      ! to conv. at level k (kg/kg/s)

      real :: dtraekp1_c2(n_sh,ntra)  ! Increment to model tracer due
                                      ! to conv. at level k+1 (kg/kg/s)

      real :: duek_c2(n_sh)           ! Increment to model U in layer k
                                      ! due to CMT (m/s2)

      real :: duekp1_c2(n_sh)         ! Increment to model U in layer
                                      ! k+1 due to CMT (m/s2)

      real :: dvek_c2(n_sh)           ! Increment to model V in layer k

      real :: dvekp1_c2(n_sh)         ! Increment to model V in layer
                                      ! k+1 due to CMT (m/s2)

      real :: flxk_c(n_sh), flxk_c2(n_sh) !Parcel mass flux in layer k
                                      ! (Pa/s)

      real :: flxkp12_c2(n_sh)        ! Half level mass flux (Pa/s)

      real :: prekp1_c2(n_sh)         ! Precip. from parcel as it rises
                                      ! from layer k to k+1 (kg/m2/s)

      real :: qpk_c(n_sh), qpk_c2(n_sh) ! Parcel mixing ratio in
                                      ! layer k(kg/kg)
      real :: qpk(n_sh)               !ad. entrain.

      real :: qpkp1_c(n_sh), qpkp1_c2(n_sh) ! Parcel mixing ratio
                                      ! in layer k+1 (kg/kg)

      real :: qek_c(n_sh), qek_c2(n_sh) ! Env. mixing ratio in
                                      ! layer k (kg/kg)

      real :: qek(n_sh)               !ad. entrain.

      real :: qekp1_c(n_sh), qekp1_c2(n_sh) ! Env. mixing ratio in
                                      ! layer k+1 (kgkg-1)

      real :: qekp1(n_sh)               !ad. entrain.

      real :: qsek_c2(n_sh)           ! Saturation mixing ratio of
                                      ! cld. env. in layer k (kg/kg)

      real :: qsek(n_sh)              !ad. entrain.

      real :: qsekp1_c(n_sh), qsekp1_c2(n_sh) ! Saturation mixing ratio
                                      ! of cld. env. in layer k+1
                                      ! (kg/kg)

      real :: qsekp1(n_sh)            !ad. entrain.

      real :: thek_c(n_sh), thek_c2(n_sh) ! Env. potential temp
                                      ! in layer k (K)

      real :: thek(n_sh)              !ad. entrain.

      real :: thekp1_c(n_sh), thekp1_c2(n_sh) ! Env. potential temp i
                                      ! in layer k (K)

      real :: thekp1(n_sh)            !ad. entrain.

      real :: thpk_c(n_sh), thpk_c2(n_sh) ! Parcel potential temp
                                      ! in layer k (K)

      real :: thpk(n_sh)              !ad. entrain.

      real :: thpkp1_c(n_sh), thpkp1_c2(n_sh)! Parcel potential temp
                                      ! in layer k (K)

      real :: traek_c(n_sh,ntra), traek_c2(n_sh,ntra) ! Tracer content
                                      ! cld. env. in layer k (kgkg-1)

      real :: traekp1_c(n_sh,ntra), traekp1_c2(n_sh,ntra) ! Tracer
                                      ! content of cloud env.
                                      ! in layer k+1 (kg/kg)

      real :: trapk_c(n_sh,ntra), trapk_c2(n_sh,ntra) ! Tracer cont.
                                      ! of parcel in layer k (kg/kg)

      real :: trapkp1_c(n_sh,ntra), trapkp1_c2(n_sh,ntra) ! Tracer cont.
                                      ! of parcel in layer k+1 (kg/kg)

      real :: rbuoy_c(n_sh), rbuoy_c2(n_sh) ! Buoyancy of parcel at k+1
                                      ! (Kelvin)

      real :: uek_c(n_sh), uek_c2(n_sh) ! Model U field on layer k
                                      ! (m/s)

      real :: uekp1_c(n_sh), uekp1_c2(n_sh)! Model U field on layer
                                      ! k+1 (m/s)

      real :: vek_c(n_sh), vek_c2(n_sh) ! Model V field on layer k
                                      ! (m/s)

      real :: vekp1_c(n_sh), vekp1_c2(n_sh) ! Model V field on layer
                                      ! k+1 (m/s)

      real :: upk_c(n_sh), upk_c2(n_sh) ! Parcel U in layer k
                                      ! after entrainment (m/s)


      real :: upkp1_c(n_sh), upkp1_c2(n_sh) ! Parcel U in layer k+1
                                      ! after entrainment (m/s)

      real :: vpk_c(n_sh), vpk_c2(n_sh) ! Parcel V in layer k
                                      ! after entrainment (m/s)

      real :: vpkp1_c(n_sh), vpkp1_c2(n_sh) ! Parcel V in layer k+1
                                      ! after entrainment (m/s)

      real :: xsqkp1_c(n_sh), xsqkp1_c2(n_sh) ! Excess water vapour
                                      ! in parcel at k+1 (kg/kg)

! PC2 compression arrays

      real :: qclek_c(n_sh), qclek_c2(n_sh) ! Environment liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qclekp1_c(n_sh), qclekp1_c2(n_sh) ! Environment liquid
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qcfek_c(n_sh), qcfek_c2(n_sh) ! Environment frozen
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qcfekp1_c(n_sh), qcfekp1_c2(n_sh) ! Environment frozen
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qclpk_c(n_sh), qclpk_c2(n_sh) ! Parcel liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qclpkp1_c(n_sh), qclpkp1_c2(n_sh) ! Parcel liquid
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: qcfpk_c(n_sh), qcfpk_c2(n_sh) ! Parcel frozen
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

      real :: qcfpkp1_c(n_sh), qcfpkp1_c2(n_sh) ! Parcel frozen
                                      ! condensate mixing ratio in
                                      ! layer k+1 (kg/kg)

      real :: cflek_c2(n_sh),cflekp1_c2(n_sh) ! Environment liquid water
                                      ! cloud volume ( )

      real :: cffek_c2(n_sh),cffekp1_c2(n_sh) ! Environment frozen water
                                      ! cloud volume ( )

      real :: bcfek_c2(n_sh),bcfekp1_c2(n_sh) ! Environment bulk total
                                      ! cloud volume ( )

      real :: dqclek_c2(n_sh),dqclekp1_c2(n_sh) ! Environment increments
                                      ! to liquid condensate mixing
                                      ! ratio to convection (kg/kg/s)

      real :: dqcfek_c2(n_sh),dqcfekp1_c2(n_sh) ! Environment increments
                                      ! to frozen condensate mixing
                                      ! ratio to convection (kg/kg/s)

      real :: dcflek_c2(n_sh),dcflekp1_c2(n_sh) ! Environment increments
                                      ! to liquid water cloud volume due
                                      ! to convection (/s)

      real :: dcffek_c2(n_sh),dcffekp1_c2(n_sh) ! Environment increments
                                      ! to frozen water cloud volume due
                                      ! to convection (/s)

      real :: dbcfek_c2(n_sh),dbcfekp1_c2(n_sh) ! Environment increments
                                      ! to bulk total cloud volume due
                                      ! to convection (/s)

      real :: amdetk_c2(n_sh)
      logical :: bgmk_c2(n_sh)
      logical :: bland_c2(n_sh)
      logical :: blowst_c2(n_sh)
      logical :: bterm_c2(n_sh)
      real :: cape_c2(n_sh)
      real :: cca_2d_c2(n_sh)
      real :: cclwp_c2(n_sh)
      real :: ccw_c2(n_sh)
      logical :: cumulus_c(n_sh), cumulus_c2(n_sh)
      real :: dcpbydt_c2(n_sh)
      real :: delexkp1_c2(n_sh)
      real :: delpk_c2(n_sh)
      real :: delpkp1_c2(n_sh)
      real :: delp_uv_k_c2(n_sh)
      real :: delp_uv_kp1_c2(n_sh)
      real :: depth_c2(n_sh)
      real :: dptot_c2(n_sh)
      real :: dqsthkp1_c2(n_sh)
      real :: dqsthk_c2(n_sh)
      real :: eflux_u_ud_c2(n_sh)
      real :: eflux_v_ud_c2(n_sh)
      real :: ekp14_c(n_sh),ekp14_c2(n_sh)
      real :: ekp34_c(n_sh),ekp34_c2(n_sh)
      real :: exk_c2(n_sh)
      real :: exkp1_c(n_sh),exkp1_c2(n_sh)
      real :: expi_c2(n_sh)
      integer :: icct_c2(n_sh)
      integer :: iccb_c2(n_sh)
      real :: lcclwp_c2(n_sh)
      integer :: lctop_c2(n_sh)
      integer :: lcbase_c2(n_sh)
      real :: lcca_c2(n_sh)
      logical :: L_shallow_c2(n_sh)
      logical :: L_mid_c2(n_sh)
      real :: max_cfl_c2(n_sh)
      real :: pk_c(n_sh),pk_c2(n_sh)
      real :: pkp1_c(n_sh),pkp1_c2(n_sh)
      real :: pstar_c2(n_sh)
      real :: q1_sd_c2(n_sh)
      real :: qpi_c2(n_sh)
      real :: qpixs_v_c2(n_sh)
      real :: relh_c2(n_sh)
      real :: rbuoy_p_here_c2(n_sh)
      real :: the_here_c2(n_sh)
      real :: thp_here_c2(n_sh)
      real :: qe_here_c2(n_sh)
      real :: qp_here_c2(n_sh)
      real :: rbuoy_p_old_c2(n_sh)
      real :: tcw_c2(n_sh)
      real :: thpi_c2(n_sh)
      real :: thpixs_v_c2(n_sh)
      real :: t1_sd_c2(n_sh)
      real :: xpk_c(n_sh),xpk_c2(n_sh)
      real :: xsbmin_v_c2(n_sh)
      logical :: b_nodd(n_sh)   ! points with no downdraught
      logical :: b_dd(n_sh)     ! points with downdraught on termination

!
! Loop counters
!

      integer :: i,i2,j,k,ktra,kt

!
! Model constants:
!

! PARXS start
      ! initial excess potential temperature (k) and mixing ratio
      ! (kg/kg) for deep convection
      REAL, PARAMETER :: THPIXS_DEEP= 0.2
      REAL, PARAMETER :: QPIXS_DEEP =0.0

      ! initial excess potential temperature (k) and mixing ratio
      ! (kg/kg) for shallow convection
      REAL, PARAMETER :: THPIXS_SHALLOW = 0.2
      REAL, PARAMETER :: QPIXS_SHALLOW  = 0.0

      ! initial excess potential temperature (k) and mixing ratio
      ! (kg/kg) for mid-level convection
      REAL, PARAMETER :: THPIXS_MID= 0.2
      REAL, PARAMETER :: QPIXS_MID =0.0
! PARXS end
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
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
! XSBMIN start
      ! minimum excess buoyancy to continue parcel ascent (K)
      REAL,PARAMETER:: XSBMIN = 0.2
! XSBMIN end
      REAL MPARB  !  MINIMUM (PARCEL BUOYANCY/LAYER THICKNESS) (K/PA)
      PARAMETER (MPARB = 1.0)
!
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
      REAL C_DEEP,                                                      &
                      ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     &     D_DEEP     ! MASS FLUX FROM PARCEL BUOYANCY FOR DEEP
                      ! CONVECTION
!
      REAL C_SHALLOW,                                                   &
                      ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     &     D_SHALLOW  ! MASS FLUX FROM PARCEL BUOYANCY FOR SHALLOW
                      ! CONVECTION
!
      REAL C_MID,                                                       &
                      ! CONSTANTS USED TO DETERMINE INITIAL CONVECTIVE
     &     D_MID      ! MASS FLUX FROM PARCEL BUOYANCY FOR MID-LEVEL
                      ! CONVECTION
!
      PARAMETER (C_DEEP = 5.17E-4, D_DEEP = 0.0)
      PARAMETER (C_SHALLOW = 5.17E-4, D_SHALLOW = 0.0)
      PARAMETER (C_MID = 5.17E-4, D_MID = 0.0)
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!C_MASS start

! coefficient relating sub-cloud convective velocity scale to cumulus
! mass flux for shallow convection

      REAL,PARAMETER:: C_MASS=0.03

! C_MASS end

!
!initialise SCM diagnostics
!
      Do k = 1,nlev
        Do i = 1,n_sh
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

      Do i = 1,n_sh
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

!-----------------------------------------------------------------------
! 1.0  Create saturation mixing ratio arrays
!-----------------------------------------------------------------------
! Re-calculate XSBMIN and THPIXS constants based on layer thickness (Pa)
!

      Do k = 1,nlev-1
          Do i = 1,n_sh
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

        Do i = 1,n_sh
          wsc(i) = (delthvu(i) * C_MASS * wstar(i) * G / (th(i,ntml(i)) &
     &               * (1.0 + C_VIRTUAL * q(i,ntml(i)))))**0.3333
          mb(i)  = C_MASS * wstar(i)
          zcld(i) = ztop_uv(i) - zlcl_uv(i)
          wsc_o_mb(i) = wsc(i)/mb(i)
        End Do
!
! Define the LCL
!
        If ( Sh_pert_opt == 0) Then
!
!         Define the LCL at the half level above ntml. Find 
!         environmental T at p_lcl by approximating theta there with
!         th(i,k) + constant*(th(i,k+1)-th(i,k))  where constant is 
!         tunable.  Similarly for q.
!
          Do i = 1,n_sh
            k =ntml(i)
            p_lcl(i)  = p_layer_boundaries(i,k)
            th_lcl(i) = th(i,k) + 0.1 * (th(i,k+1) - th(i,k))
            t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / PREF)**KAPPA)
            q_lcl(i)  = q(i,k) + 0.1 * (q(i,k+1) - q(i,k))
          End Do

        Else  ! Sh_pert_opt = 1
!
!         Define the LCL at the half level above ntml. Find 
!         environmental T at p_lcl by approximating theta there with
!         th(i,k) Similarly for q.
!
          Do i = 1,n_sh
            k =ntml(i)
            p_lcl(i)  = p_layer_boundaries(i,k)
            th_lcl(i) = th(i,k)
            t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / PREF)**KAPPA)
            q_lcl(i)  = q(i,k)
          End Do

        End if
!
! Calculate saturation mixing ratio at LCL
!

! DEPENDS ON: qsat
      Call QSAT(qse_lcl,t_lcl,p_lcl,n_sh)

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

      If (L_mom) then

        l_mom_gk = .false.      ! not using Gregory-Kershaw CMT
                                ! Shallow convection uses turbulence based CMT
!
! Initialize arrays required for Convective Momentum Transport(CMT)
!

        k=1
        Do i = 1,n_sh
          p_uv(k,i)     = p_layer_boundaries(i,k-1)
          phalf_uv(k,i) = p_layer_centres(i,k-1)
          ue_p(k,i)     = u(i,k)
          ve_p(k,i)     = v(i,k)
        End Do

        Do i = 1,n_sh
          nlcl_uv(i)    = ntml(i) + 1
          ntop_uv(i)    = ntpar(i) + 1
          n_0degc(i)    = freeze_lev(i)
        End Do

        Do i = 1,n_sh
          Do k = 2,nlev
            p_uv(k,i)     = p_layer_boundaries(i,k-1)
            phalf_uv(k,i) = p_layer_centres(i,k-1)
            ue_p(k,i)     = u(i,k)
            ve_p(k,i)     = v(i,k)
            exk_temp      = (p_uv(k,i)/PREF)**0.288
            rho_uv(k,i)   = 2.0 * p_uv(k,i) / (R * exk_temp *           &
     &                      (th(i,k-1) + th(i,k)))
          End Do
          plcl_uv(i)      = phalf_uv(nlcl_uv(i),i)
          ptop_uv(i)      = phalf_uv(ntop_uv(i),i)
          p_0degc_uv(i)   = phalf_uv(n_0degc(i),i)
          rho_uv(1,i)     = rho_uv(2,i)
        End Do
      End If     !L_mom


!
! Calculate theta and q pertubations (pertubation is based on
! environment buoyancy gradient)
! Reset th and q xs's at ntml
!

      Do i = 1,n_sh

         k = ntml(i)

         If (t_lcl(i) >  TM) then
            dq_sat_env = EPSILON * LC * qse_lcl(i)                      &
     &                       / (R * t_lcl(i) * t_lcl(i))
!           estimate of moist adiabatic lapse rate
            dthv_ma    = ( (LC/CP) - (1.+C_VIRTUAL)*th(i,k) )*          &
     &                    dq_sat_env*(G/CP)/(1.+(LC/CP)*dq_sat_env)
          else
            dq_sat_env = EPSILON * (LC+LF) * qse_lcl(i)                 &
     &                       / (R * t_lcl(i) * t_lcl(i))
!           estimate of moist adiabatic lapse rate (in K/m)
            dthv_ma    = ( ((LC+LF)/CP) - (1.+C_VIRTUAL)*th(i,k) )*     &
     &                    dq_sat_env*(G/CP)/(1.+((LC+LF)/CP)*dq_sat_env)
          End If

          b_calc   = t_lcl(i) * C_VIRTUAL * dq_sat_env + 1.0            &
     &               + C_VIRTUAL * qse_lcl(i)

!
!         Calculate theta_v perturbation:
!
          If ( Sh_pert_opt == 0) Then

            thv_pert(i) = -0.17 * wthvs(i) / mb(i)                      &
     &               + (th(i,k+1) * (1.0 + C_VIRTUAL                    &
     &                   * q(i,k+1)) - th(i,k)                          &
     &                   * (1.0 + C_VIRTUAL * q(i,k)))

          Else  ! Sh_pert_opt = 1
!
!           First convert moist adiabatic lapse rate to thv difference 
!           between levels k and k+1
!
            rho_k = p_layer_centres(i,k) /                              & 
     &          (R * th(i,k) * (p_layer_centres(i,k)/ PREF)**KAPPA)
            dthv_ma = -dthv_ma*                                         &
     &         (p_layer_centres(i,k+1)-p_layer_centres(i,k)) / (rho_k*G)
!
!           Make perturbation relative to a target lapse rate (namely 
!           0.6*dthv_ma, which is approximately what is seen in LES)
!
            thv_pert(i) = -0.17 * wthvs(i) / mb(i) +  0.6*dthv_ma       &
     &                    - (  th(i,k+1)*(1.0 + C_VIRTUAL*q(i,k+1))     &
     &                       - th(i,k)  *(1.0 + C_VIRTUAL*q(i,k))  )

          End if

          c_calc   = th_lcl(i) * C_VIRTUAL * (qse_lcl(i) - q_lcl(i))    &
     &               - thv_pert(i)

          thpert(i)   = -c_calc / b_calc  !ignore term in THPERT**2

          thpixs_v(i,k) = thpert(i)

          qpert(i)  = qse_lcl(i) + ((p_lcl(i) / PREF)                   &
     &                         **KAPPA) * thpert(i) * dq_sat_env        &
     &                           - q_lcl(i)

          qpixs_v(i,ntml(i))  = qpert(i)

      End Do !n_sh

!
! Set bwater=.true. on points where water will condense rather than
! ice.
! SUBROUTINE FLAG_WET
! UM Documentation paper 27, section (2B)
!

! DEPENDS ON: flag_wet
      Call FLAG_WET(bwater,th,exner_layer_centres,n_sh,n_sh,nlev)

!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! 2.1  Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
!      increment arrays and cca
!      and dqcl/dt, dqcf/dt, dcfl/dt, dcff/dt, dbcf/dt
!-----------------------------------------------------------------------

      Do k = 1,nlev
        Do i = 1,n_sh
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
          Do i = 1,n_sh
            dubydt(i,k) = 0.0
            dvbydt(i,k) = 0.0
          End Do
        End Do
      End If  ! L_mom

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,nlev
            Do i = 1,n_sh
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
          Do i = 1,n_sh
            up_flux(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_up_flx_half) then
        Do k = 1,nlev
          Do i = 1,n_sh
            up_flux_half(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_dwn_flx) then
        Do k = 1,nlev
          Do i = 1,n_sh
            dwn_flux(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_entr_up) then
        Do k = 1,nlev
          Do i = 1,n_sh
            entrain_up(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_detr_up) then
        Do k = 1,nlev
          Do i = 1,n_sh
            detrain_up(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_entr_dwn) then
        Do k = 1,nlev
          Do i = 1,n_sh
           entrain_dwn(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_detr_dwn) then
        Do k = 1,nlev
          Do i = 1,n_sh
            detrain_dwn(i,k) = 0.0
          End Do
        End Do
      End If

      If (L_mom) then
        If (flg_uw_shall) then
          Do k = 1,nlev
            Do i = 1,n_sh
              uw_shall(i,k) = 0.0
            End Do
          End Do
        End If
        If (flg_vw_shall) then
          Do k = 1,nlev
            Do i = 1,n_sh
              vw_shall(i,k) = 0.0
            End Do
          End Do
        End If
      End If  ! L_mom

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------
      Do i = 1,n_sh
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

      Do i = 1,n_sh
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
      Do i = 1,n_sh
        rbuoy_p_here(i) = 0.0
        the_here(i) = th(i,1)
        thp_here(i) = th(i,1)
        qe_here(i) = q(i,1)
        qp_here(i) = q(i,1)
        rbuoy_p_old(i) = 0.0
      End Do
!initialise ekm14
      Do i =1, n_sh
         ekm14(i) =0.0
      End do

!Initialise adaptive entrainment variables
      Do i = 1, n_sh
        thek(i)=th(i,1)
        qek(i)=q(i,1)
        qsek(i)=qse(i,1)
        thekp1(i)=th(i,2)
        qekp1(i)=q(i,2)
        qsekp1(i)=qse(i,2)
        thpk(i)=thp(i,1)
        qpk(i)=qp(i,1)
        bwk(i)=bwater(i,2)
        bwkp1(i)=bwater(i,2)
!Note that unlike p_layer_boundaries, where k indexing is offset
!by one compared to the dynamics numbering, z retains the numbering
!convention for dynamics variables i.e. for theta levels, k->k
!and for rho levels k+1/2 -> k+1
!check this with Rachel
        zk(i) = z_theta(i,1)
        zkp12(i)=z_rho(i, 2)
        zkp1(i)=z_theta(i, 2)
      End Do

!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------
! To reduce cost limit level loop to NTPAR_MAX maximum NTPAR value.
! NTPAR is the top of the parcel ascent for shallow convection.

      If (ad_on == 1) Then   ! Adaptive forced detrainment on
                             ! No limit on convection top
        ntpar_max = nlev-3   ! What is a sensible value to have here?

      Else                   ! Top limited
                         
        ntpar_max=0
        Do i = 1,n_sh
          If (ntpar(i) >  ntpar_max) then
            ntpar_max=ntpar(i)
          End If
        End Do
      End If

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

! initialize SCM diagnostics for this pass through the loop
!NB do not re-initialise rbuoy_p_old 'cos value from k-1 needed
!instead, re-initialize at end of loop after it's been used
!but before it is re-set
      Do i = 1,n_sh
        rbuoy_p_here(i) = 0.0
        the_here(i) = th(i,k)
        thp_here(i) = th(i,k)
        qe_here(i) = q(i,k)
        qp_here(i) = q(i,k)
        rbuoy_p_here_c2(i) = 0.0
        the_here_c2(i) = 0.0
        thp_here_c2(i) = 0.0
        qe_here_c2(i) = 0.0
        qp_here_c2(i) = 0.0
      End Do
!Initialise adaptive entrainment variables
      Do i = 1, n_sh
        thek(i)=th(i,k)
        qek(i)=q(i,k)
        qsek(i)=qse(i,k)
        thekp1(i)=th(i,k+1)
        qekp1(i)=q(i,k+1)
        qsekp1(i)=qse(i,k+1)
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
! Set relative humidity in layer k (rhum)
!

        Do i = 1,n_sh
          rhum(i) = q(i,k) / qse(i,k)
        End Do

!-----------------------------------------------------------------------
! 3.1  Calculate layer dependent constants (pressure,
!      layer thickness, entrainment coefficients, detrainment
!      coefficients)
!-----------------------------------------------------------------------

! DEPENDS ON: layer_cn
        Call LAYER_CN(k,n_sh,n_sh,nlev,                                 &
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
     &                ,thek, qek, qsek, thekp1,qekp1,qsekp1             &
     &                ,thpk, qpk,                                       &
     &                bwk,bwkp1,ekm14                                   &
     &                ,zk, zkp12, zkp1                                  &
     &                )

! Set ekm14 for next pass through loop
         Do i = 1, n_sh
           ekm14(i) = ekp14(i)
         End Do

! Calculate dqs/dth for layers k and k+1 (subroutine DQS_DTH)
!

        If (k == 2) then
! DEPENDS ON: dqs_dth
          Call DQS_DTH(dqsthk,k,th(1,k),qse(1,k),exk,n_sh)
        else
          Do i = 1,n_sh
            dqsthk(i) = dqsthkp1(i)
          End Do
        End If

! DEPENDS ON: dqs_dth
        Call DQS_DTH(dqsthkp1,k+1,th(1,k+1),qse(1,k+1),exkp1,n_sh)

!
! Set other grid dependent constants
!

        Do i = 1,n_sh

!
! Maximum initial convective mass flux
!

          flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)

        End Do
!
! Set initial parcel properties (theta,q,tracer,momentum) if convection
! is not occurring at level k
!

        Do i = 1,n_sh
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
            If (L_mom_gk) then
              up(i,k) = u(i,k)
              vp(i,k) = v(i,k)
            End If
          End If
        End Do  ! n_sh
        If (L_tracer) then
          Do ktra=1,ntra
            Do i = 1,n_sh
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
        Do i = 1,n_sh

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

        End Do  ! n_sh

!
! Calculate number of points which may convect (ncposs) and
! set compression indices (index1)
!
        ncposs = 0
        Do i = 1,n_sh
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
        Call LIFT_PAR(ncposs,n_sh,thpkp1_c,qpkp1_c,xsqkp1_c,            &
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
! Reset threashold for forced detrainment to the initial 
! (potentially negative) buoyancy
!
          xsbmin_v(index1(i),k) = thv_pert(index1(i))

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
            rbuoy_p_here_c2(i)=rbuoy_p_here(index1(index2(i)))
            the_here_c2(i)=the_here(index1(index2(i)))
            thp_here_c2(i)=thp_here(index1(index2(i)))
            qe_here_c2(i)=qe_here(index1(index2(i)))
            qp_here_c2(i)=qp_here(index1(index2(i)))
            rbuoy_p_old_c2(i)=rbuoy_p_old(index1(index2(i)))
            dptot_c2(i)   = 0.0 ! dummy variable
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
          If (L_mom_gk) then
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

! Force shallow convection to stop at ntpar unless using adaptive forced 
! detrainment.

          If (ad_on == 0) Then
           
! Part of BTERM comes from values in array before convec2.
! Original code set bterm_c2 according to this test before call to
! convec2.
            Do i = 1,nconv
              If (k  ==  ntpar(index1(index2(i))))then
                bterm_c2(i) = .true.
              End If
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
        Call CONVEC2(nconv,n_sh,nlev,k,thek_c2,thekp1_c2,qek_c2,        &
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
     &               upk_c2,vpk_c2,upkp1_c2,vpkp1_c2,duek_c2,duekp1_c2, &
     &               dvek_c2,dvekp1_c2,eflux_u_ud_c2,eflux_v_ud_c2,     &
     &               delp_uv_k_c2, delp_uv_kp1_c2,                      &
     &               thpixs_v_c2,qpixs_v_c2,                            &
     &               xsbmin_v_c2,L_shallow_c2,L_mid_c2,                 &
     &               L_tracer,ntra,traek_c2,traekp1_c2,                 &
     &               trapk_c2,trapkp1_c2,dtraek_c2,dtraekp1_c2,cape_c2, &
     &               dcpbydt_c2,max_cfl_c2,timestep,rbuoy_p_here_c2,    &
     &               the_here_c2,thp_here_c2,qe_here_c2,qp_here_c2,     &
     &               rbuoy_p_old_c2, ad_on, sdet_on, new_termc,         &
     &               deltak_c2, l_calc_dxek,l_q_interact,               &
     &               flxkp12_c2,cumulus_c2,relh_c2,dptot_c2)

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

        Do i = 1,n_sh
          thp(i,k+1)    = 0.0
          qp(i,k+1)     = 0.0
          xpk(i,k+1)    = 0.0
          flx(i,k+1)    = 0.0
          depth(i)      = 0.0
          precip(i,k+1) = 0.0
          qclp(i,k+1)   = 0.0
          qcfp(i,k+1)   = 0.0
        End Do
        Do i = 1,n_sh
          bgmk(i)       = .false.
          bterm(i)      = .false.
        End Do

        If (L_tracer) then
          Do ktra = 1,ntra
            Do i = 1,n_sh
              trap(i,k+1,ktra) = 0.0
            End Do
          End Do
        End If

        If (L_mom_gk) then
          Do i = 1,n_sh
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
          qe_here(index1(index2(i)))  = qe_here_c2(i)
          qp_here(index1(index2(i)))  = qp_here_c2(i)
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
        If (L_mom_gk) then
          Do i = 1,nconv
            dubydt(index1(index2(i)),k) = duek_c2(i)
            dvbydt(index1(index2(i)),k) = dvek_c2(i)
            dubydt(index1(index2(i)),k+1) = duekp1_c2(i)
            dvbydt(index1(index2(i)),k+1) = dvekp1_c2(i)
            eflux_u_ud(index1(index2(i))) = eflux_u_ud_c2(i)
            eflux_v_ud(index1(index2(i))) = eflux_v_ud_c2(i)
          End Do
        End If
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
          Do i = 1,nconv
            up_flux(index1(index2(i)),k+1) = flxk_c2(i)
          End Do
        End If
!NB up_flux_half is k+1 because flxkp12 is level above
! p_layer_centres on k which is z_rho on k+1
        If (flg_up_flx_half) then
          Do i = 1,nconv
            up_flux_half(index1(index2(i)),k+1) = flxkp12_c2(i)
          End Do
        End If
        End if      ! nconv > 0

!
!   Write adaptive diagnostics for this level to full array for output
!
        Do i = 1,n_sh
          rbuoy_p_out(i,k) = rbuoy_p_here(i)
          the_out(i,k) = the_here(i)
          thp_out(i,k) = thp_here(i)
          qe_out(i,k) = qe_here(i)
          qp_out(i,k) = qp_here(i)
        End Do

!re-initialize r_buoy_p_old ready for next level
      Do i = 1,n_sh
        rbuoy_p_old(i) = 0.0
      End do


!  write rbuoy_here to rbuoy_p_old for next pass through loop
        Do i = 1,n_sh
           rbuoy_p_old(i) = rbuoy_p_here(i)
        End Do
!-----------------------------------------------------------------------
! 3.4  CFL scaling
!-----------------------------------------------------------------------
!
! Set up integer nterm which is the total number of points where
! convection has terminated.
! Index to full array (n_sh) with index_nterm
!

        nterm = 0
        Do i = 1,n_sh
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
              
              ! set the flx_init to the new value to provide the real initial mass
              ! flux under all conditions              
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

        Else ! original code

          Do j = 1,nterm
            i = index_nterm(j)
            If (flx_init_new(i)  >   0.0) Then
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
             
              ! set the flx_init to the new value to provide the real initial mass
              ! flux under all conditions              
              flx_init(i) = flx_init_new(i)
      
            End If
          End Do      ! j (nterm)

        End If      ! lcv_ccrad


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
              If (L_mom) then
                dubydt(i,kt)  = dubydt(i,kt) * scale_f(i)
                dvbydt(i,kt)  = dvbydt(i,kt) * scale_f(i)
                  If (kt <  k+1) then
                    flxkp12(kt,i) = flxkp12(kt,i) * scale_f(i)
                  End if
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
              If (kt <  k+1) then
                If (flg_up_flx_half) then
                  up_flux_half(i,kt+1) = flxkp12(kt,i)
                End If
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
        Do i = 1,n_sh
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
         call DD_ALL_CALL(n_sh,npossdd,kmax_term,nlev,trlev,ntra        &
                        , kterm,iccb,icct,index_possdd,l_tracer         &
                        , flg_dwn_flx,flg_entr_dwn,flg_detr_dwn         &
                        , bwater(1,2),exner_layer_centres               &
                        , exner_layer_boundaries,p_layer_centres        &
                        , p_layer_boundaries,pstar,recip_pstar,timestep &
                        , cca_2d,thp,qp,th,q,trap,tracer,flx,precip     &
                        , dthbydt,dqbydt,dtrabydt,rain,snow,rain_3d     &
                        , snow_3d,dwn_flux,entrain_dwn,detrain_dwn)

         End if

!-----------------------------------------------------------------------
! 4.2 Surface precipitation calculation for terminating points with
!     no downdraught (moved outside level loop) ie do this calculation
!     on all points at the end.
!-----------------------------------------------------------------------
! Points where no downdraught possible
        nnodd = 0
        Do i = 1,n_sh

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
            End If
          End do
! Only add 1 if kmax_term is less than model levels (which should be
! true).
          If (kmax_term  <  nlev ) then
             kmax_term = kmax_term + 1
          End If

!
! Surface precipitation calculation
!
! DEPENDS ON: evap_bcb_nodd_all
          Call EVAP_BCB_NODD_ALL(n_sh,nnodd,kmax_term,kterm,iccb        &
     &,                      index_nodd,bwater(1,2)                     &
     &,                      exner_layer_centres,exner_layer_boundaries &
     &,                      p_layer_centres, p_layer_boundaries,pstar  &
     &,                      timestep , cca_2d, th, q, precip           &
     &,                      dthbydt, dqbydt                            &
     &,                      rain, snow, rain_3d, snow_3d)


        End If

!
! Adjust cloud base, top and amount to prevent errors occurring in
! radiation scheme when iccb = icct (this happens when convection
! saturates upon forced detrainment).
!
         
        If (.NOT. lcv_ccrad) Then
          Do i=1, n_sh
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
! 5.0  Convective Momentum Transport (if L_mom = .T.)
!-----------------------------------------------------------------------

      If (L_mom) then
        nterm = 0

        Do i = 1, n_sh
          nterm = nterm + 1
          cu_term(nterm) = i
          cu_tend(nterm) = i
        End Do

        If (nterm  >   0) then

! DEPENDS ON: shallow_grad_stress
          Call SHALLOW_GRAD_STRESS(n_sh,n_sh,nterm,nlev,cu_term,        &
     &                             nlcl_uv,ntop_uv,mb,wsc,wstar,zcld,   &
     &                             plcl_uv,ptop_uv,p_uv,phalf_uv,       &
     &                             rho_uv,ue_p,ve_p,timestep,           &
                                   ! IN
     &                             uw,vw)

! DEPENDS ON: shallow_base_stress
          Call SHALLOW_BASE_STRESS(n_sh,n_sh,n_sh,nlev,nterm,cu_term,   &
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
          Call SHALLOW_CMT_INCR(n_sh,n_sh,n_sh,nlev,nterm,cu_term,      &
     &                          cu_tend,nlcl_uv,ntop_uv,uw,vw,phalf_uv, &
     &                          rho_uv,zlcl_uv,                         &
                                !OUT
     &                          dubydt,dvbydt)

        End If  ! nterm>0
      End If ! L_mom


!-----------------------------------------------------------------------
! 6.0  Mix the increments from the initial parcel perturbation throughout
!      the subcloud layer if this option is selected.
!-----------------------------------------------------------------------
      If (bl_cnv_mix == 1) then
      
! DEPENDS ON: mix_ipert
        Call MIX_IPERT(n_sh, nlev, nbl, ntml, p_layer_boundaries, &
                     exner_layer_centres, dthbydt, dqbydt, flx_init, &
                     thpert, qpert)
      End If

!-----------------------------------------------------------------------
! 7.0  End Subroutine
!-----------------------------------------------------------------------

      Return
      END SUBROUTINE SHALLOW_CONV
