
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+ Gather-scatter routine for deep and shallow convection points
!

      SUBROUTINE GLUE_CONV(NP_FIELD,NPNTS,NLEV,NBL,TH,Q,QCL,QCF         &
     &,                 CF_LIQUID,CF_FROZEN,BULK_CF,PSTAR               &
     &,                 BLAND,U,V                                       &
     &,                 TRACER, DTHBYDT, DQBYDT, DQCLBYDT, DQCFBYDT     &
     &,                 DCFLBYDT, DCFFBYDT, DBCFBYDT, DUBYDT, DVBYDT    &
     &,                 RAIN, SNOW, RAIN_3D, SNOW_3D                    &
     &,                 CCA,ICCB,ICCT,CCLWP,CCW                         &
     &,                 LCBASE,LCTOP,LCCA,freeze_lev                    &
     &,                 l_mid_all,kterm_deep                            &
     &,                 precip_deep,precip_shall,precip_mid,precip_cong &
     &,                 wstar_dn,wstar_up,mb1,mb2,kterm_congest         &
     &,                 N_CUMULUS,UW0,VW0,W_MAX,zlcl,ZLCL_UV,ZTOP_UV    &
     &,                 LCCLWP,CAPE_OUT                                 &
     &,                 N_DP,n_cg,N_SH                                  &
     &,                 r_rho,r_theta,rho,rho_theta                     &
     &,                 EXNER_LAYER_BOUNDARIES                          &
     &,                 EXNER_LAYER_CENTRES                             &
     &,                 P_LAYER_BOUNDARIES                              &
     &,                 P_LAYER_CENTRES                                 &
     &,                 z_theta, z_rho                                  &
     &,                 TIMESTEP,T1_SD,Q1_SD                            &
     &,                 NTML,NTPAR,NBDSC,NTDSC,L_SHALLOW_BL             &
     &,                 l_pc2_diag_sh_pts, l_congestus                  &
     &,                 CUMULUS_BL,WSTAR,WTHVS,DELTHVU_BL,ql_ad,ftl,fqt &
     &,                 L_TRACER,NTRA,TRLEV,N_CCA_LEV                   &
     &,                 l_mixing_ratio,L_CALC_DXEK,L_Q_INTERACT         &
     &,                 up_flux_half,flg_up_flx_half                    &
     &,                 UP_FLUX,FLG_UP_FLX,DWN_FLUX,FLG_DWN_FLX         &
     &,                 ENTRAIN_UP,FLG_ENTR_UP,DETRAIN_UP               &
     &,                 FLG_DETR_UP,ENTRAIN_DWN,FLG_ENTR_DWN            &
     &,                 DETRAIN_DWN,FLG_DETR_DWN,UW_DEEP,FLG_UW_DEEP    &
     &,                 VW_DEEP,FLG_VW_DEEP,UW_SHALL,FLG_UW_SHALL       &
     &,                 VW_SHALL,FLG_VW_SHALL                           &
     &,                 flg_wqt_flux_sh,flg_wthetal_flux_sh             &
     &,                 flg_wthetav_flux_sh,flg_wql_flux_sh             &
     &,           flg_mf_deep,flg_mf_congest,flg_mf_shall,flg_mf_midlev &
     &,           flg_dt_deep,flg_dt_congest,flg_dt_shall,flg_dt_midlev &
     &,           flg_dq_deep,flg_dq_congest,flg_dq_shall,flg_dq_midlev &
     &,           flg_du_deep,flg_du_congest,flg_du_shall,flg_du_midlev &
     &,           flg_dv_deep,flg_dv_congest,flg_dv_shall,flg_dv_midlev &
     &,         wqt_flux_sh,wthetal_flux_sh,wthetav_flux_sh,wql_flux_sh &
     &,           mf_deep,mf_congest,mf_shall,mf_midlev                 &
     &,           dt_deep,dt_congest,dt_shall,dt_midlev                 &
     &,           dq_deep,dq_congest,dq_shall,dq_midlev                 &
     &,           du_deep,du_congest,du_shall,du_midlev                 &
     &,           dv_deep,dv_congest,dv_shall,dv_midlev                 &
     &                  )

! Purpose:
!   Gather-scatter routine for deep and shallow convection points.
!
!   Called by Load_bal_conv
!
! Current owners of code: R A Stratton
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
!
      Use cv_cntl_mod, Only:                                            &
          lcv_3d_cca, lcv_ccrad
      Use cv_run_mod,  Only:                                            &
          l_mom,        l_convcld_hadgem1, adapt,       termconv,       &
          cca2d_sh_opt, cca2d_md_opt,      cca2d_dp_opt,                &
          cca_sh_knob,  cca_md_knob,       cca_dp_knob,                 &
          ccw_sh_knob,  ccw_md_knob,       ccw_dp_knob,                 &
          bl_cnv_mix

      IMPLICIT NONE

! Model Constants required

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


!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!

      integer, intent(in) :: np_field ! No of points in a full field
                                     !(note all multi-dimensional fields
                                     !  passed to routine MUST be
                                     !   dimensioned with np_field)

      integer, intent(in) :: npnts    ! No. of points in segment

      integer, intent(in) :: nlev     ! No. of model layers

      integer, intent(in) :: ntml(npnts) ! Top level of surface mixed
                                      ! layer defined relative to
                                      ! theta,q grid

      integer, intent(in) :: ntpar(npnts) ! Top level of initial parcel
                                      ! ascent in BL scheme defined
                                      ! relative to theta,q grid

      integer, intent(in) :: n_cca_lev! No. of convective cloud
                                      ! amount levels (1 for 2D,
                                      ! nlevs for 3D)

      integer, intent(in) :: nbl      ! No. of boundary layer levels

      integer, intent(in) :: ntra     ! No. of tracer fields

      integer, intent(in) :: trlev    ! No. of model levels on which
                                      ! tracers are included

      integer, intent(in) :: n_dp     ! Number of deep points

      integer, intent(in) :: n_cg     ! Number of congestus points

      integer, intent(in) :: n_sh     ! Number of shallow points

      real, intent(in) :: delthvu_bl(npnts) !Integral of undilute parcel
                                      ! buoyancy over convective cloud
                                      ! layer (Kelvin m)

      real, intent(in) ::     &
        ql_ad(npnts)          & ! adiabatic liquid water content
      , ftl(npnts)            & ! Surface sensible heat flux from BL
                                ! (W/m2) i.e. cp*rho*w'tl'
      , fqt(npnts)              ! Total water flux from surface
                                ! (kg/m2/s) i.e. rho*w'qT'

      real, intent(in) :: exner_layer_centres(np_field,0:nlev) !Exner

      real, intent(in) :: exner_layer_boundaries(np_field,0:nlev) !Exner
                                      ! at half level above
                                      ! exner_layer_centres

      real, intent(in) :: pstar(npnts) ! Surface pressure (Pa)

      real, intent(in) :: p_layer_centres(np_field,0:nlev) !Pressure(Pa)


      real, intent(in) :: p_layer_boundaries(np_field,0:nlev) ! Pressure
                                      ! at half level above
                                      ! p_layer_centres (Pa)

      real, intent(in) :: z_theta(np_field,nlev) ! height of theta
                                      ! levels above surface (m)

      real, intent(in) :: z_rho(np_field,nlev)    ! height of rho levels
                                      ! above surface (m)


      real, intent(in) :: t1_sd(np_field) ! Standard deviation of
                                      ! turbulent fluctuations of
                                      ! layer 1 temp. (K)

      real, intent(in) :: q1_sd(np_field) ! Standard deviation of
                                      ! turbulent fluctuations of
                                      ! layer 1 q (kg/kg)

      real, intent(in) :: th(np_field,nlev) !Model potential
                                      ! temperature (K)

      real, intent(in) :: q(np_field,nlev) ! Model mixing ratio (kg/kg)

      real, intent(in) :: timestep    ! Model timestep (s)

      real, intent(in) :: uw0(np_field) ! U-comp of surface stress(N/m2)

      real, intent(in) :: vw0(np_field) ! V-comp of surface stress(N/m2)

      real, intent(in) :: W_MAX(np_field) ! max w in column

      real, intent(in) :: u(np_field,nlev) !Model U field (m/s)

      real, intent(in) :: v(np_field,nlev) !Model V field (m/s)

      real, intent(in) :: wstar(npnts) ! Convective velocity scale
                                       ! (m/s)

      real, intent(in) :: wthvs(npnts) ! Surface flux of THV  (Pa m/s2)

      real, intent(in) :: zlcl(npnts) ! Lifting condensation level accurate
                                      ! height (m) (NOT USED)

      real, intent(in) :: zlcl_uv(npnts) !Lifting condensation level
                                      ! defined for the uv grid (m)

      real, intent(in) :: ztop_uv(npnts) ! Top of cloud layer
                                      ! defined for the uv
                                      ! grid (m)

      real, intent(in) :: r_rho(np_field,nlev)      ! radius rho lev (m)
      real, intent(in) :: r_theta(np_field,0:nlev)  ! radius th lev (m)
      real, intent(in) :: rho(np_field,nlev)        ! density kg/m3
      real, intent(in) :: rho_theta(np_field,nlev)  ! th lev kg/m3

      logical, intent(in) :: L_tracer ! Switch for inclusion of tracers

      logical, intent(in) :: l_mixing_ratio 
                        ! true for moisture inputs as mixing ratio
                        ! false for moisture inputs as specific humidity

      logical, intent(in) :: L_shallow_bl(npnts) ! Shallow cumulus
                                      ! indicator

      logical, intent(in) :: L_pc2_diag_sh_pts(npnts) ! Carry
                                      ! diagnostic shallow convective
                                      ! information for PC2

      logical, intent(in) :: L_congestus(npnts) ! congestus cumulus

      logical, intent(in) :: cumulus_bl(npnts) ! Cumulus indicator

      logical, intent(in) :: L_calc_dxek ! Switch for calculation of
                                      ! condensate increment

      logical, intent(in) :: L_q_interact ! Switch allows overwriting
                                      ! parcel variables when
                                      ! calculating condensate incr.

      logical, intent(in) :: flg_up_flx ! STASH flag for updraught
                                      ! mass flux

      logical, intent(in) :: flg_up_flx_half ! STASH flag for updraught
                                      ! mass flux on half levels

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
                                      ! stress

      logical, intent(in) :: flg_uw_shall ! STASH flag for x-comp
                                      ! shallow stress

      logical, intent(in) :: flg_vw_shall ! STASH flag for y-comp
                                      ! shallow stress

      logical, intent(in) :: bland(npnts) ! Land/sea mask

!
! Flags for flux diagnostics - only possible from turbulence based
! schemes
!
      logical, intent(in) ::                                            &
     &  flg_wqt_flux_sh                                                 &
                               ! STASH flag for wqt_flux_sh
     &, flg_wql_flux_sh                                                 &
                               ! STASH flag for wql_flux_sh
     &, flg_wthetal_flux_sh                                             &
                               ! STASH flag for wthetal_flux_sh
     &, flg_wthetav_flux_sh    ! STASH flag for wthetav_flux_sh

!
! Arguments with intent INOUT:
!

      real, intent(inout) :: qcl(np_field,nlev) ! Liq condensate mix
                                      ! ratio (kg/kg)

      real, intent(inout) :: qcf(np_field,nlev) ! Ice condensate mix
                                      ! ratio (kg/kg)

      real, intent(inout) :: cf_liquid(np_field,nlev)
                                ! Liq water cloud volume (fraction)

      real, intent(inout) :: cf_frozen(np_field,nlev)
                                ! Frozen water cloud volume (fraction?)

      real, intent(inout) :: bulk_cf(np_field,nlev) ! Bulk total cloud
                                       ! volume ( )

      real, intent(inout) :: tracer(np_field,trlev,ntra) !Model tracer
                                      ! fields (kg/kg)

!
! Arguments with intent OUT:
!

      real, intent(out) :: dqclbydt(np_field,nlev) ! Increments to liq
                                      ! condensate due to convection
                                      ! (kg/kg/s)

      real, intent(out) :: dqcfbydt(np_field,nlev) ! Increments to ice
                                      ! condensate due to convection
                                      ! (kg/kg/s)

      real, intent(out) :: dcflbydt(np_field,nlev) ! Increments to liq
                                      ! cloud volume due to convection
                                      ! (/s)

      real, intent(out) :: dcffbydt(np_field,nlev) ! Increments to ice
                                      ! cloud volume due to convection
                                      ! (/s)

      real, intent(out) :: dbcfbydt(np_field,nlev) ! Increments to
                                      ! total cld volume due to
                                      ! convection(/s)

      real, intent(out) :: dthbydt(np_field,nlev) ! Increments to
                               ! potential temp. due to convection (K/s)

      real, intent(out) :: dqbydt(np_field,nlev) ! Increments to q due
                                      ! to convection (kg/kg/s)

      real, intent(out) :: dubydt(np_field,nlev+1) ! Increments to U due
                                      ! to CMT (m/s2)

      real, intent(out) :: dvbydt(np_field,nlev+1) ! Increments to V due
                                      ! to CMT (m/s2)

      real, intent(out) :: rain(npnts) ! Surface convective rainfall
                                      ! (kg/m2/s)

      real, intent(out) :: snow(npnts) ! Surface convective snowfall
                                      ! (kg/m2/s)

      real, intent(out) :: rain_3d(npnts,nlev) ! convective rainfall flux
                                      ! (kg/m2/s)

      real, intent(out) :: snow_3d(npnts,nlev) ! convective snowfall flux
                                      ! (kg/m2/s)

      real, intent(out) :: cca(np_field,n_cca_lev) ! Convective cloud
                                      ! amount (%)

      integer, intent(out) :: iccb(npnts) ! Convective cloud base
                                      ! level (m)

      integer, intent(out) :: icct(npnts) ! Convective cloud top
                                      ! level (m)

      real, intent(out) :: cclwp(npnts) ! Condensed water path (k/m2)

      real, intent(out) :: ccw(np_field,nlev) ! Convective cloud liquid
                                      ! water on model levels (g/kg)

      real, intent(out) :: lcca(npnts) ! Lowest conv. cloud amt. (%)

      integer, intent(out) :: lcbase(npnts) ! Lowest conv. cloud base
                                      ! level (m)

      integer, intent(out) :: lctop(npnts) ! Lowest conv. cloud top
                                      ! level (m)

      integer, intent(out) :: freeze_lev(npnts) !index for freezing lev

      integer, intent(out) :: kterm_deep(npnts) ! index deep conv

      real, intent(out) :: precip_deep(npnts) ! deep precip (kg/m2/s)

      real, intent(out) :: precip_shall(npnts) ! shallow precip(kg/m2/s)

      real, intent(out) :: precip_mid(npnts) ! mid precip (kg/m2/s)
      real, intent(out) :: precip_cong(npnts) !congestus precip(kg/m2/s)

      real, intent(out) :: lcclwp(npnts) ! Condensed water path for
                                      ! lowest conv. cld. level (kg/m2)

      real, intent(out) :: up_flux(np_field,nlev) ! Updraught mass flux
                                      ! (Pa/s)

      real, intent(out) :: up_flux_half(np_field,nlev)
                                      ! Updraught mass flux
                                      ! on half levels (Pa/s)

      real, intent(out) :: dwn_flux(np_field,nlev) ! Downdraught mass
                                      ! flux (Pa/s)

      real, intent(out) :: entrain_up(np_field,nlev) ! Fractional
                                      ! entrainment rate into updraughts
                                      ! (Pa/s)

      real, intent(out) :: detrain_up(np_field,nlev) ! Fractional
                                      ! detrainment rate into updraughts
                                      ! (Pa/s)

      real, intent(out) :: entrain_dwn(np_field,nlev) ! Fractional
                                      ! entrainment rate into
                                      ! downdraughts (Pa/s)

      real, intent(out) :: detrain_dwn(np_field,nlev) ! Fractional
                                      ! detrainment rate into
                                      ! downdraughts (Pa/s)

      real, intent(out) :: uw_deep(np_field,nlev) ! X-comp. of stress
                                      ! from deep convection
                                      !(kg/m/s2)

      real, intent(out) :: vw_deep(np_field,nlev) ! Y-comp. of stress
                                      ! from deep convection
                                      !(kg/m/s2)

      real, intent(out) :: uw_shall(np_field,nlev) ! X-comp. of stress
                                      ! from shallow convection
                                      !(kg/m/s2)

      real, intent(out) :: vw_shall(np_field,nlev) ! Y-comp. of stress
                                      ! from shallow convection
                                      !(kg/m/s2)

!      real, intent(out) ::  &  ! Intend to make a diagnostic in future
      real ::               &
      uw_mid(np_field,nlev) & ! U comp of stress from mid convection (kg/m/s2)
     ,vw_mid(np_field,nlev)   ! V comp of stress from mid convection (kg/m/s2)

      real, intent(out) :: cape_out(npnts) ! Saved convective available
                                      ! potential energy for diagnostic
                                      ! output (Jkg-1)
!
! Flux diagnostics - only avialable from turbulence based schemes
!
      real, intent(out) ::                                              &
     &  wqt_flux_sh(np_field,nlev)                                      &
                                   ! w'qt' flux
     &, wql_flux_sh(np_field,nlev)                                      &
                                   ! w'ql' flux
     &, wthetal_flux_sh(np_field,nlev)                                  &
                                       ! w'thetal' flux
     &, wthetav_flux_sh(np_field,nlev)                                  &
                                       ! w'thetav' flux
     &, wstar_dn(npnts)                                                 &
                                !
     &, wstar_up(npnts)                                                 &
                                !
     &, mb1(npnts)                                                      &
                           !
     &, mb2(npnts)                                                      &
                           !
     &, mf_deep(np_field,nlev),mf_congest(np_field,nlev)                &
     &, mf_shall(np_field,nlev),mf_midlev(np_field,nlev)                &
     &, dt_deep(np_field,nlev),dt_congest(np_field,nlev)                &
     &, dt_shall(np_field,nlev),dt_midlev(np_field,nlev)                &
     &, dq_deep(np_field,nlev),dq_congest(np_field,nlev)                &
     &, dq_shall(np_field,nlev),dq_midlev(np_field,nlev)                &
     &, du_deep(np_field,nlev),du_congest(np_field,nlev)                &
     &, du_shall(np_field,nlev),du_midlev(np_field,nlev)                &
     &, dv_deep(np_field,nlev),dv_congest(np_field,nlev)                &
     &, dv_shall(np_field,nlev),dv_midlev(np_field,nlev)

      integer, intent(out) ::                                           &
     &  kterm_congest(npnts)   ! termination level for congestus

! flag not required for 4A scheme used by 5A
      logical, intent(in) ::                                            &
     &            flg_mf_deep,flg_mf_congest,flg_mf_shall,flg_mf_midlev &
     &,           flg_dt_deep,flg_dt_congest,flg_dt_shall,flg_dt_midlev &
     &,           flg_dq_deep,flg_dq_congest,flg_dq_shall,flg_dq_midlev &
     &,           flg_du_deep,flg_du_congest,flg_du_shall,flg_du_midlev &
     &,           flg_dv_deep,flg_dv_congest,flg_dv_shall,flg_dv_midlev
!-----------------------------------------------------------------------
! Redundant arguments
!-----------------------------------------------------------------------

      integer, intent(in) :: n_cumulus
      real, intent(in)    :: ntdsc
      real, intent(in)    :: nbdsc

! local required in move mixing to glue

      real :: dtrabydt(npnts,nlev,ntra) ! Increment to tracer due to
                                      ! convection (kg/kg/s)

!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to DEEP convection scheme.
! Arrays are identified by underscore DP (_dp) and are of length
! n_dp where n_dp is the number of points diagnosed as deep in the
! boundary layer diagnosis routine. For full desriptions of variables
! see above.
!-----------------------------------------------------------------------

      integer :: dp_new_termc !flag for simplified termination of
                              !convection
                              !0 = off, 1 = on

      integer :: dp_on        !flag for adaptive applied to deep conv
                              !0 = off, 1 = on

      integer :: mdet_dp_on   !flag for adaptive mixing detrainment
                              !applied to deep conv
                              !0 = off, 1 = on

      integer :: dp_ent_on    !flag for adaptive mod for entrainment
                              !0 = off, 1 = on
                              
      integer :: dp_sdet_on   !flag for smoothed forced detrainment
                              !0 = off, 1 = on                        

      integer :: dpi(n_dp)        ! index for deep points in full grid

      integer :: ntml_dp(n_dp)

      integer :: ntpar_dp(n_dp)

      logical :: bland_dp(n_dp)

      real ::                                                           &
       pstar_dp(n_dp)                                                   &
     , recip_pstar_dp(n_dp)                                             &
     , t1_sd_dp(n_dp)                                                   &
     , q1_sd_dp(n_dp)                                                   &
     , uw0_dp(n_dp)                                                     &
     , vw0_dp(n_dp)                                                     &
     , zlcl_uv_dp(n_dp)                                                 &
     , wstar_dp(n_dp)                                                   &
     , delthvu_dp(n_dp)

      real ::                                                           &
       p_layer_centres_dp(n_dp,0:nlev)                                  &
     , p_layer_boundaries_dp(n_dp,0:nlev)                               &
     , exner_layer_centres_dp(n_dp,0:nlev)                              &
     , exner_layer_boundaries_dp(n_dp,0:nlev)                           &
     , r2rho_dp(n_dp,nlev)                                              &
     , r2rho_th_dp(n_dp,nlev)                                           &
     , rho_dp(n_dp,nlev)                                                &
     , rho_theta_dp(n_dp,nlev)                                          &
     , dr_across_th_dp(n_dp,nlev)                                       &
     , dr_across_rh_dp(n_dp,nlev)                                       &
     , z_theta_dp(n_dp,nlev)                                            &
     , z_rho_dp(n_dp,nlev)                                              &
     , r_theta_dp(n_dp,0:nlev)                                          &
     , r_rho_dp(n_dp,nlev)

      real :: W_MAX_dp(n_dp)

      real :: u_dp(n_dp,nlev)

      real :: v_dp(n_dp,nlev)

      real :: th_dp(n_dp,nlev)

      real :: q_dp(n_dp,nlev)

      real :: tracer_dp(n_dp,trlev,ntra)

      real :: dthbydt_dp(n_dp,nlev)

      real :: dqbydt_dp(n_dp,nlev)

      real :: dubydt_dp(n_dp,nlev+1)

      real :: dvbydt_dp(n_dp,nlev+1)

      real :: rain_dp(n_dp)

      real :: snow_dp(n_dp)

      real :: rain_3d_dp(n_dp,nlev)

      real :: snow_3d_dp(n_dp,nlev)

      real :: tcw_dp(n_dp)

      integer :: iccb_dp(n_dp)

      integer :: icct_dp(n_dp)

      real :: cclwp_dp(n_dp)

      real :: ccw_dp(n_dp,nlev)

      real :: lcca_dp(n_dp)

      integer :: lcbase_dp(n_dp)

      integer :: lctop_dp(n_dp)

      integer :: freeze_lev_dp(n_dp)

      real :: lcclwp_dp(n_dp)

      real :: up_flux_dp(n_dp,nlev)

      real :: up_flux_half_dp(n_dp,nlev)

      real :: dwn_flux_dp(n_dp,nlev)

      real :: entrain_up_dp(n_dp,nlev)

      real :: detrain_up_dp(n_dp,nlev)

      real :: entrain_dwn_dp(n_dp,nlev)

      real :: detrain_dwn_dp(n_dp,nlev)

      real :: uw_deep_dp(n_dp,nlev)

      real :: vw_deep_dp(n_dp,nlev)

      real :: cape_out_dp(n_dp)

      real :: qse_dp(n_dp,nlev)

      real :: qcl_dp(n_dp,nlev)

      real :: qcf_dp(n_dp,nlev)

      real :: cf_liquid_dp(n_dp,nlev)

      real :: cf_frozen_dp(n_dp,nlev)

      real :: bulk_cf_dp(n_dp,nlev)

      real :: dqclbydt_dp(n_dp,nlev)

      real :: dqcfbydt_dp(n_dp,nlev)

      real :: dcflbydt_dp(n_dp,nlev)

      real :: dcffbydt_dp(n_dp,nlev)

      real :: dbcfbydt_dp(n_dp,nlev)

      real :: dtrabydt_dp(n_dp,nlev,ntra) ! Increment to tracer due to
                                          ! convection (kg/kg/s)
      integer :: kterm_dp(n_dp)    ! required by mid level scheme

      real :: cca_2d_dp(n_dp)      ! required by mid level scheme

!  compressed SCM diagnostics for adaptive modset

      real :: rbuoy_p_out_dp(n_dp,nlev)

      real :: the_out_dp(n_dp,nlev)

      real :: thp_out_dp(n_dp,nlev)

      real :: qe_out_dp(n_dp,nlev)

      real :: qp_out_dp(n_dp,nlev)
!temporary local diagnostics
      integer :: thp_out_copy_count
      integer :: thp_out_md_count
      integer :: th_count

!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to SHALLOW convection scheme.
! Arrays are identified by underscore SH (_sh) and are of length
! n_sh where n_sh is the number of points diagnosed as shallow in the
! boundary layer diagnosis routine.For full desriptions of variables
! see above.
!-----------------------------------------------------------------------

      integer :: sh_new_termc !flag for simplified termination of
                              !convection

      integer :: sh_on        !flag for adaptive applied to shallow
                              !conv

      integer :: mdet_sh_on   !flag for adaptive mixing detrainment
                              !applied to shallow conv

      integer :: sh_ent_on    !flag for adaptive applied to entrainment
      
      integer :: sh_sdet_on   !flag for smoothed forced detrainment
                              !0 = off, 1 = on                        
      
      integer :: shi(n_sh)      ! index for shallow points in full grid

      integer :: ntml_sh(n_sh)

      integer :: ntpar_sh(n_sh)

      real :: pstar_sh(n_sh)

      real :: recip_pstar_sh(n_sh)

      real :: p_layer_centres_sh(n_sh,0:nlev)

      real :: p_layer_boundaries_sh(n_sh,0:nlev)

      logical :: bland_sh(n_sh)

      real :: exner_layer_centres_sh(n_sh,0:nlev)

      real :: exner_layer_boundaries_sh(n_sh,0:nlev)

      real :: z_theta_sh(n_sh,nlev), z_rho_sh(n_sh,nlev)                &
     &, r2rho_th_sh(n_sh,nlev)                                          &
     &, dr_across_th_sh(n_sh,nlev)

      real :: delthvu_sh(n_sh)

      real :: t1_sd_sh(n_sh)

      real :: q1_sd_sh(n_sh)

      real :: uw0_sh(n_sh)

      real :: vw0_sh(n_sh)

      real :: wstar_sh(n_sh)

      real :: wthvs_sh(n_sh)

      real :: zlcl_uv_sh(n_sh)

      real :: ztop_uv_sh(n_sh)

      real :: u_sh(n_sh,nlev)

      real :: v_sh(n_sh,nlev)

      real :: th_sh(n_sh,nlev)

      real :: q_sh(n_sh,nlev)

      real :: tracer_sh(n_sh,trlev,ntra)

      real :: dthbydt_sh(n_sh,nlev)

      real :: dqbydt_sh(n_sh,nlev)

      real :: dubydt_sh(n_sh,nlev+1)

      real :: dvbydt_sh(n_sh,nlev+1)

      real :: rain_sh(n_sh)

      real :: snow_sh(n_sh)

      real :: rain_3d_sh(n_sh,nlev)

      real :: snow_3d_sh(n_sh,nlev)

      real :: tcw_sh(n_sh)

      integer :: iccb_sh(n_sh)

      integer :: icct_sh(n_sh)

      real :: cclwp_sh(n_sh)

      real :: ccw_sh(n_sh,nlev)

      real :: lcca_sh(n_sh)

      integer :: lcbase_sh(n_sh)

      integer :: lctop_sh(n_sh)

      integer :: freeze_lev_sh(n_sh)

      real :: lcclwp_sh(n_sh)

      real :: up_flux_sh(n_sh,nlev)

      real :: up_flux_half_sh(n_sh,nlev)

      real :: dwn_flux_sh(n_sh,nlev)

      real :: entrain_up_sh(n_sh,nlev)

      real :: detrain_up_sh(n_sh,nlev)

      real :: entrain_dwn_sh(n_sh,nlev)

      real :: detrain_dwn_sh(n_sh,nlev)

      real :: uw_shall_sh(n_sh,nlev)

      real :: vw_shall_sh(n_sh,nlev)

      real :: cape_out_sh(n_sh)

      real :: qse_sh(n_sh,nlev)

      real :: qcl_sh(n_sh,nlev)

      real :: qcf_sh(n_sh,nlev)

      real :: cf_liquid_sh(n_sh,nlev)

      real :: cf_frozen_sh(n_sh,nlev)

      real :: bulk_cf_sh(n_sh,nlev)

      real :: dqclbydt_sh(n_sh,nlev)

      real :: dqcfbydt_sh(n_sh,nlev)

      real :: dcflbydt_sh(n_sh,nlev)

      real :: dcffbydt_sh(n_sh,nlev)

      real :: dbcfbydt_sh(n_sh,nlev)
      real :: dtrabydt_sh(n_sh,nlev,ntra) ! Increment to tracer due to
                                     ! convection (kg/kg/s)
      real :: cca_2d_sh(n_sh)        ! required by mid level scheme

!  compressed SCM diagnostics for adaptive modset

      real :: rbuoy_p_out_sh(n_sh,nlev)

      real :: the_out_sh(n_sh,nlev)

      real :: thp_out_sh(n_sh,nlev)

      real :: qe_out_sh(n_sh,nlev)

      real :: qp_out_sh(n_sh,nlev)

!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to MID-LEVEL convection scheme.
! Arrays are identified by underscore MD (_md) and are of length
! npnts where npnts is the total number of points in the grid (since
! mid-level convection may occur on any point previously diagnosed as
! shallow or deep).
!-----------------------------------------------------------------------

!  compressed SCM diagnostics for adaptive modset

      real :: rbuoy_p_out_md(npnts,nlev)

      real :: the_out_md(npnts,nlev)

      real :: thp_out_md(npnts,nlev)

      real :: qe_out_md(npnts,nlev)

      real :: qp_out_md(npnts,nlev)

      integer :: md_new_termc !flag for simplified termination of
                              !convection
                              !0 = off, 1 = on

      integer :: md_on        !flag for adaptive applied to mid conv
                              !0 = off, 1 = on

      integer :: mdet_md_on   !flag for adaptive mixing detrainment
                               !0 = off, 1 = on
                             !applied to mid conv

      integer :: md_ent_on    !Flag for adaptive mod for entrainment

      integer :: md_sdet_on   !flag for smoothed forced detrainment
                              !0 = off, 1 = on                        

      integer :: midtrig(npnts)   ! Level at which mid level convection
                                  ! may start

      real :: p_layer_centres_md(npnts,0:nlev)

      real :: p_layer_boundaries_md(npnts,0:nlev)

      real :: exner_layer_centres_md(npnts,0:nlev)

      real :: exner_layer_boundaries_md(npnts,0:nlev)

      real :: z_theta_md(npnts,nlev), z_rho_md(npnts,nlev)
      real :: r_theta_md(npnts,0:nlev), r_rho_md(npnts,nlev)

      real :: u_md(npnts,nlev)

      real :: v_md(npnts,nlev)

      real :: th_md(npnts,nlev)

      real :: q_md(npnts,nlev)

      real :: tracer_md(npnts,trlev,ntra)

      real :: dthbydt_md(npnts,nlev)

      real :: dqbydt_md(npnts,nlev)

      real :: dubydt_md(npnts,nlev+1)

      real :: dvbydt_md(npnts,nlev+1)

      real :: rain_md(npnts)

      real :: snow_md(npnts)

      real :: rain_3d_md(npnts,nlev)

      real :: snow_3d_md(npnts,nlev)

      real :: tcw_md(npnts)

      integer :: iccb_md(npnts)

      integer :: icct_md(npnts)

      real :: cclwp_md(npnts)

      real :: ccw_md(npnts,nlev)

      real :: lcca_md(npnts)

      integer :: lcbase_md(npnts)

      integer :: lctop_md(npnts)

      integer :: freeze_lev_md(npnts)

      real :: lcclwp_md(npnts)

      real :: up_flux_md(npnts,nlev)

      real :: up_flux_half_md(npnts,nlev)

      real :: dwn_flux_md(npnts,nlev)

      real :: entrain_up_md(npnts,nlev)

      real :: detrain_up_md(npnts,nlev)

      real :: entrain_dwn_md(npnts,nlev)

      real :: detrain_dwn_md(npnts,nlev)

      real :: cape_out_md(npnts)

      real :: qcl_md(npnts,nlev)

      real :: qcf_md(npnts,nlev)

      real :: cf_liquid_md(npnts,nlev)

      real :: cf_frozen_md(npnts,nlev)

      real :: bulk_cf_md(npnts,nlev)

      real :: dqclbydt_md(npnts,nlev)

      real :: dqcfbydt_md(npnts,nlev)

      real :: dcflbydt_md(npnts,nlev)

      real :: dcffbydt_md(npnts,nlev)

      real :: dbcfbydt_md(npnts,nlev)
      real :: dtrabydt_md(npnts,nlev,ntra) ! Increment to tracer due to
                                      ! convection (kg/kg/s)
      real :: qse_md(npnts,nlev)      ! Saturation mixing ratio of
                                      ! cloud environment (kg/kg)
      real :: cca_2d_md(npnts)        ! required by mid level scheme

      logical :: l_mid_all(npnts)     ! true if mid level convection


!-----------------------------------------------------------------------
! Mixing ratio / specific humidity conversions
!-----------------------------------------------------------------------
!  The full model holds q, qcl and qcf, all of which should be used in the
!  conversions and evaluation of qt. (High resolution versions can also
!  hold other moist variables but these will be ignored here.) The current
!  mass flux scheme without PC2 assume qt=q i.e. it ignores the presence
!  of qcl and qcf. This assumption is also true of the new shallow
!  turbulence based scheme. This makes it slightly more difficult to decide
!  how to do the conversions.
!
!  Let m indicate mixing ratio and q specific then:
!
!  If PC2                        else
!      qt = qv + qcl + qcf            qt = qv
!      mt = mv + mcl + mcf            mt = mv
!
!  Not PC2
!      dqcl = 0.0,      dqcf = 0.0
!      dmcl = 0.0,      dmcf = 0.0
!
!  Before scheme
!  mv(t) = qv(t)/(1-qt(t))          or qv(t) = mv(t)/(1+mt(t))
!
!  after scheme
!  mv(t+dt) = qv(t+dt)/(1-qt(t+dt)) or qv(t+dt) = mv(t+dt)/(1+mt(t+dt))
!
!  where
!  dmv= mv(t+dt) -  mv(t)           or dqv  =  qv(t+dt)- qv(t)
!
! How to convert increments
!
!  dmv = (dqv + qv*dqt - qt*dqv)/[(1-qt)(1-qt-dqt)]
!
! if qt = qv  then dmv = dqv /[(1-qt)(1-qt-dqt)]
!
!   or
!
!  dqv = (dmv + mt*dmv - mv*dmt)/[(1+mt)(1+mt+dmt)]
!
! if mt = mv  then dqv = dmv /[(1+mt)(1+mt+dmt)]
!
! Conversion of qsat
! ------------------
! rsat - mixing ratio saturation value
!
! qsat=rsat/(1+rsat)      and rsat=qsat(1-qsat)
!
!-----------------------------------------------------------------------
!  Arrays required for mixing ratio to specific conversions

      real :: qt(npnts,nlev)    ! total water variable either mixing
                                ! ratio or specific humidity (kg/kg)
                                ! qt = q + qcl + qcf
      real ::                &
        dqt                  &  ! total moisture increment per timestep
      , dqtt                 &  ! total moisture increment per s
      , denom                   ! 1./denominator
!-----------------------------------------------------------------------
!  Arrays required by 3d CCA calculation

      real :: cca_2d(npnts)        ! 2d convective cloud

! Saturation mixxing ratio calulation and 1/pstar

      real :: recip_pstar(npnts)      ! Reciprocal of pstar array

      real :: qse(npnts,nlev)         ! Saturation mixing ratio of
                                      ! cloud environment (kg/kg)
      real :: pt(npnts)               ! Temporary store for P in calc.
                                      ! of sat. mixing ratio. (Pa)
      real :: tt(npnts)               ! Temporary store for T in calc.
                                      ! of saturation mixing ratio. (K)
      real :: ttkm1(npnts)            ! Temporary store for T in layer
                                      ! k-1 for use in freezing level
                                      ! calc. for anvil. (K)


!  compressed SCM diagnostics for adaptive modset

      real :: rbuoy_p_out(npnts,nlev)

      real :: the_out(npnts,nlev)

      real :: thp_out(npnts,nlev)

      real :: qe_out(npnts,nlev)

      real :: qp_out(npnts,nlev)

! copies so we can see if mid convection modifies the profile

      real :: rbuoy_p_out_copy(npnts,nlev)

      real :: the_out_copy(npnts,nlev)

      real :: thp_out_copy(npnts,nlev)

      real :: qe_out_copy(npnts,nlev)

      real :: qp_out_copy(npnts,nlev)

!  arrays required by grid calculations

      real ::                    &
       r2rho_th(npnts,nlev)      & ! radius**2 density theta lev (kg/m)
      ,r2rho(npnts,nlev)         & ! radius**2 density rho lev (kg/m)
      ,dr_across_th(npnts,nlev)  & ! thick of theta levels (m)
      ,dr_across_rh(npnts,nlev)    ! thick of rho levels (m)


      !===============================================================
      ! CCRad Variables
      !===============================================================
      real             :: overlap_fac(n_sh)
                                      ! Factor designed to improve
                                      ! shallow Cu cover by allowing
                                      ! for non-vertical clouds.
      real             :: zcld(n_sh)  ! Cloud depth, m (for shallow).

      real             :: mb          ! Used in calculation of CCA_2d
      real             :: wsc         ! of shallow cloud via climate
      real             :: zpr         ! method (BL Fluxes) only if
                                      !   lcv_ccrad       = T .AND.
                                      !   cca2d_sh_opt  = 1
                                      
      integer          :: dum_iccb    ! Holding variables for
      integer          :: dum_icct    ! identifying cld base and tops
                                      ! on gridpoints with more than 1
                                      ! (possibly 2) mid-level
                                      ! convective events
                                      !
      integer          :: n_mdcld     ! Number of gridpoints with
                                      ! single mid-level cloud banks
                                      !
      integer          :: n_mdcld_mult   ! Number of gridpoints with
                                         ! multiple mid-level cloud
                                         ! banks
                                         !
      real             :: a_land = 0.3   ! parameters used to relate
      real             :: a_sea  = 0.3   ! cca_2d of cld to its
      real             :: b_land = 0.025 ! precipitation rate
      real             :: b_sea  = 0.025 !


      !===============================================================
      ! Allocatable arrays, because we do not know how many gridpoints
      ! will contain mid-level until we test for them.
      !===============================================================

      integer                            :: dum1(npnts)
      integer                            :: dum2(npnts)

      integer, dimension(:), allocatable :: mdcldi
                                      ! INDICES in full array of
                                      ! gridpoints single mid-level
                                      ! cloud banks
      integer, dimension(:), allocatable :: mdcldi_mult
                                      ! index in full array of
                                      ! gridpoints multiple mid-level
                                      ! cloud banks


      !===============================================================
      ! Compressed arrays for gridpoints which have one bank of
      ! mid-level cloud and multiple banks of mid-level cloud.
      ! Requires the use of compressed allocatable arrays because the
      ! location and number of gridpoints with mid-level has not been
      ! diagnosed yet.
      !===============================================================

      !===============================================================
      ! For multiple mid-level cloud
      !===============================================================
      integer, dimension(:)  , allocatable ::  iccb_md_c
      integer, dimension(:)  , allocatable ::  icct_md_c
      integer, dimension(:)  , allocatable ::  freeze_lev_md_c
      real   , dimension(:)  , allocatable ::  cca_2d_md_c
      real   , dimension(:,:), allocatable ::  cca_md_c
      real   , dimension(:,:), allocatable ::  ccw_md_c
      real   , dimension(:,:), allocatable ::  z_theta_md_c
      real   , dimension(:,:), allocatable ::  z_rho_md_c
      real   , dimension(:,:), allocatable ::  p_lyr_bnds_md_c

      !===============================================================
      ! For multiple mid-level cloud
      !===============================================================
      integer, dimension(:)  , allocatable ::  iccb_md_mult_c
      integer, dimension(:)  , allocatable ::  icct_md_mult_c
      integer, dimension(:)  , allocatable ::  freeze_lev_md_mult_c
      real   , dimension(:)  , allocatable ::  cca_2d_md_mult_c
      real   , dimension(:,:), allocatable ::  cca_md_mult_c
      real   , dimension(:,:), allocatable ::  ccw_md_mult_c
      real   , dimension(:,:), allocatable ::  z_theta_md_mult_c
      real   , dimension(:,:), allocatable ::  z_rho_md_mult_c
      real   , dimension(:,:), allocatable ::  p_lyr_bnds_md_mult_c

      real :: cca_dp(n_dp , nlev) ! Arrays to hold cloud fractions  
      real :: cca_sh(n_sh , nlev) ! of shallow/mid/deep convection
      real :: cca_md(npnts, nlev) !

      INTEGER, Parameter   :: cca2d_total_condensed_water = 0
      INTEGER, Parameter   :: cca2d_grant_lock            = 1 ! Shallow cnv
      INTEGER, Parameter   :: cca2d_srf_precip_rate       = 1 ! mid/deep cnv

      !===============================================================
      ! End CCRad Variables local variables
      !===============================================================

!  arrays required by energy correction

      integer :: index1(npnts)
      integer :: nconv_all         ! Total number of points convecting

!   required by check on -ve q

      real :: qMinInColumn(npnts)     ! Minimum value for q in column
                                      ! (kg/kg)
      real :: temp1(npnts)            ! work array

! array required by tracers at end

      real :: limited_step(npnts)     ! Reduced step size for tracer
                                      ! mixing

      real :: step_test(npnts)        ! Work array used in reducing step

!      real :: reduction_factor(npnts,ntra) ! Timestep reduction factor
                                      !for tracers

!
! Parameters
!

      real, parameter :: QMIN = 1.0E-8 ! Global minimum allowed Q

      real, parameter :: SAFETY_MARGIN = 1.0E-100 ! Small number used in
                                      ! tracer step reduction



!
! Loop counters
!

      integer :: i,j,k,ktra,idp


!
! External routines called:
!
      EXTERNAL                                                          &
     &  DEEP_CONV,SHALLOW_CONV,MID_CONV,COR_ENGY,MIX_INC,QSAT_MIX

!-----------------------------------------------------------------------
! Initialise output variables
!-----------------------------------------------------------------------


! initialise 3d rain varibles
!
      Do k=1, nlev
        Do i=1, n_sh
          rain_3d_sh(i,k) = 0.0
          snow_3d_sh(i,k) = 0.0
        End Do
      End Do

      Do k=1, nlev
        Do i=1, n_dp
          rain_3d_dp(i,k) = 0.0
          snow_3d_dp(i,k) = 0.0
        End Do
      End Do

      Do k=1, nlev
        Do i=1, npnts
          rain_3d_md(i,k) = 0.0
          snow_3d_md(i,k) = 0.0
        End Do
      End Do

      Do i = 1,npnts
        cape_out(i)                = 0.0
        cclwp(i)                   = 0.0
        rain(i)                    = 0.0
        snow(i)                    = 0.0
        precip_deep(i)             = 0.0
        precip_shall(i)            = 0.0
        precip_mid(i)              = 0.0
        precip_cong(i)             = 0.0
      End Do

 ! initialize adaptive en/detrainment variables
 ! full arrays initialized to environmental profile
      Do k = 1, nlev
        Do i = 1, npnts
          rbuoy_p_out(i,k)=0.0
          the_out(i,k)=th(i,k)
          thp_out(i,k)=th(i,k)
          qe_out(i,k)=q(i,k)
          qp_out(i,k)=q(i,k)
          rbuoy_p_out_copy(i,k)=0.0
          the_out_copy(i,k)=th(i,k)
          thp_out_copy(i,k)=th(i,k)
          qe_out_copy(i,k)=q(i,k)
          qp_out_copy(i,k)=q(i,k)
          rbuoy_p_out_md(i,k)=0.0
          the_out_md(i,k)=th(i,k)
          thp_out_md(i,k)=th(i,k)
          qe_out_md(i,k)=q(i,k)
          qp_out_md(i,k)=q(i,k)
        End Do
      End Do

 ! compressed arrays initialised to zero (will be
 ! initialised to compressed versions of environmental
 ! profile within DPCONV, SHCONV
      Do k = 1, nlev
        Do i = 1, n_dp
          rbuoy_p_out_dp(i,k)=0.0
          the_out_dp(i,k)=0.0
          thp_out_dp(i,k)=0.0
          qe_out_dp(i,k)=0.0
          qp_out_dp(i,k)=0.0
        End Do
      End Do
      Do k = 1, nlev
        Do i = 1, n_sh
          rbuoy_p_out_sh(i,k)=0.0
          the_out_sh(i,k)=0.0
          thp_out_sh(i,k)=0.0
          qe_out_sh(i,k)=0.0
          qp_out_sh(i,k)=0.0
        End Do
      End Do

      Do k = 1,nlev
        Do i = 1,npnts
          dqbydt(i,k)              = 0.0
          dthbydt(i,k)             = 0.0
          dqclbydt(i,k)            = 0.0
          dqcfbydt(i,k)            = 0.0
          dcflbydt(i,k)            = 0.0
          dcffbydt(i,k)            = 0.0
          dbcfbydt(i,k)            = 0.0
          ccw(i,k)                 = 0.0
          rbuoy_p_out(i,k) = 0.0
          the_out(i,k) = th(i,k)
          thp_out(i,k) = th(i,k)
          qe_out(i,k) = q(i,k)
          qp_out(i,k) = q(i,k)
        End Do
      End Do

! Initialisation of qt used for mixing/specific conversions
! Definition dependent on PC2.

      If (L_q_interact) Then
        Do k=1,nlev
          Do i=1,npnts
            qt(i,k) = q(i,k) + qcl(i,k) + qcf(i,k)
          End Do
        End Do
      Else
        Do k=1,nlev
          Do i=1,npnts
            qt(i,k) = q(i,k)
          End Do
        End Do
      End If

      If (L_mom) then
        Do k=1,nlev+1
          Do i = 1,npnts
            dubydt(i,k)              = 0.0
            dvbydt(i,k)              = 0.0
          End Do
        End Do
      End If
      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,nlev
            Do i = 1,npnts
              dtrabydt(i,k,ktra)  = 0.0   ! tracer increments
            End Do
          End Do
        End Do
      End If

! Initialise diagnostics if requested

      If (flg_up_flx) then
        Do k = 1,nlev
          Do i = 1,npnts
            up_flux(i,k)             = 0.0
          End Do
        End Do
      End If
      If (flg_up_flx_half) then
        Do k = 1,nlev
          Do i = 1,npnts
            up_flux_half(i,k)             = 0.0
          End Do
        End Do
      End If
      If (flg_dwn_flx) then
        Do k = 1,nlev
          Do i = 1,npnts
            dwn_flux(i,k)            = 0.0
          End Do
        End Do
      End If
      If (flg_entr_up) then
        Do k = 1,nlev
          Do i = 1,npnts
            entrain_up(i,k)          = 0.0
          End Do
        End Do
      End If
      If (flg_detr_up) then
        Do k = 1,nlev
          Do i = 1,npnts
            detrain_up(i,k)          = 0.0
          End Do
        End Do
      End If
      If (flg_entr_dwn) then
        Do k = 1,nlev
          Do i = 1,npnts
            entrain_dwn(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_detr_dwn) then
        Do k = 1,nlev
          Do i = 1,npnts
            detrain_dwn(i,k)         = 0.0
          End Do
        End Do
      End If
!
! Extra variable initialisation for safety
!

      Do i = 1,npnts
        iccb(i)     = 0
        icct(i)     = 0
        lcca(i)     = 0.0
        lcbase(i)   = 0
        lctop(i)    = 0
        lcclwp(i)   = 0.0
      End Do

      If (flg_uw_deep) then
        Do k = 1,nlev
          Do i = 1,npnts
            uw_deep(i,k)  = 0.0
          End Do
        End Do
      End If
      If (flg_vw_deep) then
        Do k = 1,nlev
          Do i = 1,npnts
            vw_deep(i,k)  = 0.0
          End Do
        End Do
      End If
      If (flg_uw_shall) then
        Do k = 1,nlev
          Do i = 1,npnts
            uw_shall(i,k) = 0.0
          End Do
        End Do
      End If
      If (flg_vw_shall) then
        Do k = 1,nlev
          Do i = 1,npnts
            vw_shall(i,k) = 0.0
          End Do
        End Do
      End If

!
! Required to get same convective cloud as old scheme
! Need to pass values from deep and shallow to mid level.
! Currently conversion of 2d to 3d makes use of total 2d from
! shallow/deep and mid in a column. If in future the conversion
! does not work on a column total this could be removed.
!
      Do i=1,npnts
        tcw_md(i) = 0.0
        iccb_md(i) = 0
        icct_md(i) = 0
        cca_2d_md(i) = 0.0
        kterm_deep(i) = 0
      End Do

      If (lcv_ccrad) Then

        ! Initialise mid-level arrays
        Do i=1, npnts
          lcbase_md(i) = 0
          lctop_md(i)  = 0
          l_mid_all(i) = .FALSE.
          cca_2d_md(i) = 0.0
        End Do

        Do k=1, nlev
          Do i=1, npnts
            ccw_md(i,k) = 0.0
            cca_md(i,k) = 0.0
          End Do
        End Do


        ! Initialise shallow arrays
        Do i=1, n_sh
          lcbase_sh(i) = 0
          iccb_sh(i)   = 0
          icct_sh(i)   = 0
          cca_2d_sh(i) = 0.0
        End Do

        Do k=1, nlev
          Do i=1, n_sh
            ccw_sh(i,k) = 0.0
            cca_sh(i,k) = 0.0
          End Do
        End Do


        ! Initialise deep arrays
        Do i=1, n_dp
          lcbase_dp(i) = 0
          iccb_dp(i)   = 0
          icct_dp(i)   = 0
          cca_2d_dp(i) = 0.0
        End Do

        Do k=1, nlev
          Do i=1, n_dp
            ccw_dp(i,k) = 0.0
            cca_dp(i,k) = 0.0
          End Do
        End Do

      End If      ! lcv_ccrad

! Initialize flags to control whether adaptive scheme is on or off
! Default initialised to off, reset just before call to routine if req.
      sh_on = 0
      md_on = 0
      dp_on = 0
      mdet_dp_on = 0
      mdet_md_on = 0
      mdet_sh_on = 0
      sh_ent_on=0
      md_ent_on=0
      dp_ent_on=0
      sh_new_termc=0
      md_new_termc=0
      dp_new_termc=0
      sh_sdet_on=0
      md_sdet_on=0
      dp_sdet_on=0

!Set flags for adaptive scheme depending on values passed in from umui
!flags all initialized to zero, so no changes needed if adapt  ==  0
      if(adapt  ==  1) then !HadGEM1a
        md_on = 1
        dp_on = 1
        mdet_dp_on = 1
        mdet_md_on = 1
      else if(adapt  ==  2) then  !convection test code
        md_on = 1
        dp_on = 1
        mdet_dp_on = 1
        mdet_md_on = 1
        md_ent_on=1
        dp_ent_on=1
      else if(adapt  ==  3) then  !operational
        mdet_dp_on = 1
        dp_on = 1
      else if(adapt  ==  4) then  ! Possible HadGEM3 as option 1 plus shallow
        md_on = 1
        dp_on = 1
        mdet_dp_on = 1
        mdet_md_on = 1
        sh_on = 1
!        mdet_sh_on = 1   Mixing detrainment not set on as has no impact
      else if(adapt  ==  5) then !adapt det + smoothed forced det for deep + mid
        md_on = 1
        dp_on = 1
        mdet_dp_on = 1
        mdet_md_on = 1
        md_sdet_on=1
        dp_sdet_on=1
      else if(adapt  ==  6) then ! as 5 + shallow
        sh_on = 1
        md_on = 1
        dp_on = 1
!       mdet_sh_on = 1  Mixing detrainment not set on as has no impact 
        mdet_md_on = 1
        mdet_dp_on = 1
        sh_sdet_on=1
        md_sdet_on=1
        dp_sdet_on=1
      end if

      if(termconv  ==  0) then
        md_new_termc=0
        dp_new_termc=0
      else if(termconv  ==  1) then
        md_new_termc=1
        dp_new_termc=1
        If (adapt == 4 .or. adapt == 6) Then
          sh_new_termc=1          ! use new termination condition
        End If
      end if
!-----------------------------------------------------------------------
!  grid information
!-----------------------------------------------------------------------
! Calculate quantities involving density etc for future use

      Do k=1,nlev
        Do i=1,npnts
          dr_across_rh(i,k) = r_theta(i,k) - r_theta(i,k-1)
          r2rho(i,k)        = r_rho(i,k)*r_rho(i,k)*rho(i,k)
        End Do
      End Do
!
! rho_theta only set for nlev-1
!
      k=1     ! bottom theta level thicker
        Do i=1,npnts
          dr_across_th(i,k) = r_rho(i,k+1) - r_theta(i,0)
          r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k)
        End Do

      Do k=2,nlev-1
        Do i=1,npnts
          dr_across_th(i,k) = r_rho(i,k+1) - r_rho(i,k)
          r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k)
        End Do
      End Do

      k=nlev     ! top layer  (hope not used ?
                 !             assume density as layer below)
        Do i=1,npnts
          dr_across_th(i,k) = r_theta(i,nlev) - r_rho(i,k)
          r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k-1)
        End Do


!-----------------------------------------------------------------------
! 1.0 Section to calculate fields used by all convection types.
! Note cheaper to do here than in individual routines.
! Create saturation mixing ratio arrays  & calculate freeze level
!-----------------------------------------------------------------------

!
! Calculate 1/pstar and initialize freeze_lev array.
!

      Do i = 1,npnts
        recip_pstar(i)=1.0 / pstar(i)
        freeze_lev(i) = 1
      End Do

!
! Loop over levels
!

      Do k = 1,nlev

!
! Find freezing level
!

        If (k  ==  1) then
          Do i = 1,npnts
            tt(i) = th(i,k) * exner_layer_centres(i,k)
            pt(i) = p_layer_centres(i,k)
!
! Commented out as initialisation sets freeze_lev to 1.
! Code left incase altered in future.
!
!            If (tt(i)  <   TM) then
!              freeze_lev(i) = k
!            End If

          End Do
        else
          Do i = 1,npnts
            ttkm1(i) = tt(i)
            tt(i) = th(i,k) * exner_layer_centres(i,k)
            pt(i) = p_layer_centres(i,k)
            If (tt(i)  <   TM .and. ttkm1(i)  >=  TM) then
              If (freeze_lev(i) == 1) then
                 freeze_lev(i) = k
              End if
            End If
          End Do
        End If

!
! Calculate saturation specific humidity/mixing ratio
!

! DEPENDS ON: qsat_mix
        Call QSAT_mix(qse(1,k),tt,pt,npnts,l_mixing_ratio)

      End Do  ! nlev

!-----------------------------------------------------------------------
! 1.0 DEEP Convection
! 1.1 Compress input variable arrays for deep convection scheme to
!     length n_dp (deep points only)
!-----------------------------------------------------------------------
      If (n_dp  >   0) then
        j = 0
        Do i = 1,npnts
          If (cumulus_bl(i).and..not.l_shallow_bl(i)) then
            j                        = j+1
            dpi(j)                   = i
          End If
        End Do
!
! In only variables
!
        Do j=1,n_dp
          bland_dp(j)           = bland(dpi(j))
          ntml_dp(j)            = ntml(dpi(j))
          ntpar_dp(j)           = ntpar(dpi(j))
          pstar_dp(j)           = pstar(dpi(j))
          recip_pstar_dp(j)     = recip_pstar(dpi(j))
          q1_sd_dp(j)           = q1_sd(dpi(j))
          t1_sd_dp(j)           = t1_sd(dpi(j))
          uw0_dp(j)             = uw0(dpi(j))
          vw0_dp(j)             = vw0(dpi(j))
          W_MAX_dp(j)             = W_MAX(dpi(j))
          wstar_dp(j)           = wstar(dpi(j))
          zlcl_uv_dp(j)         = zlcl_uv(dpi(j))
          freeze_lev_dp(j)      = freeze_lev(dpi(j))
          delthvu_dp(j)         = delthvu_bl(dpi(j))
        End Do

        Do k = 0,nlev
          Do j=1,n_dp
            p_layer_centres_dp(j,k)    = p_layer_centres(dpi(j),k)
            p_layer_boundaries_dp(j,k) = p_layer_boundaries(dpi(j),k)

            exner_layer_centres_dp(j,k)    = exner_layer_centres(dpi(j),k)
            exner_layer_boundaries_dp(j,k) = exner_layer_boundaries(dpi(j),k)
            r_theta_dp(j,k)  = r_theta(dpi(j),k)
          End Do
        End Do
        Do k = 1,nlev
          Do j=1,n_dp
            u_dp(j,k)           = u(dpi(j),k)
            v_dp(j,k)           = v(dpi(j),k)
            th_dp(j,k)          = th(dpi(j),k)
            z_theta_dp(j,k)     = z_theta(dpi(j),k)
            z_rho_dp(j,k)       = z_rho(dpi(j),k)
            r_rho_dp(j,k)       = r_rho(dpi(j),k)
            r2rho_th_dp(j,k)     = r2rho_th(dpi(j),k)
            r2rho_dp(j,k)        = r2rho(dpi(j),k)
            rho_theta_dp(j,k)    = rho_theta(dpi(j),k)
            rho_dp(j,k)          = rho(dpi(j),k)
            dr_across_th_dp(j,k) = dr_across_th(dpi(j),k)
            dr_across_rh_dp(j,k) = dr_across_rh(dpi(j),k)
          End Do
        End Do
!
! G-R mass flux scheme requires input of specific humidity
!
        If (l_mixing_ratio) Then

          Do k = 1,nlev
            Do j = 1,n_dp
              q_dp(j,k)   = q(dpi(j),k)  /(1.0+qt(dpi(j),k))
              qse_dp(j,k) = qse(dpi(j),k)/(1.0+qse(dpi(j),k))
              qcl_dp(j,k) = qcl(dpi(j),k)/(1.0+qt(dpi(j),k))
              qcf_dp(j,k) = qcf(dpi(j),k)/(1.0+qt(dpi(j),k))
            End Do
          End Do

        Else   ! Input is specific humidity therefore no problems

          Do k = 1,nlev
            Do j = 1,n_dp
              q_dp(j,k)   = q(dpi(j),k)
              qse_dp(j,k) = qse(dpi(j),k)
              qcl_dp(j,k) = qcl(dpi(j),k)
              qcf_dp(j,k) = qcf(dpi(j),k)
            End Do
          End Do

        End If ! (l_mixing_ratio)

!
! In/out variables
!

        Do k = 1,nlev
          Do j=1,n_dp
            cf_liquid_dp(j,k)   = cf_liquid(dpi(j),k)
            cf_frozen_dp(j,k)   = cf_frozen(dpi(j),k)
            bulk_cf_dp(j,k)     = bulk_cf(dpi(j),k)
          End Do
        End Do

        If (L_tracer) then
          Do ktra = 1,ntra
            Do k = 1,trlev
               Do j=1,n_dp
                tracer_dp(j,k,ktra)  = tracer(dpi(j),k,ktra)
               End Do
            End Do
            Do k = 1,nlev
               Do j=1,n_dp
                dtrabydt_dp(j,k,ktra)  = 0.0
               End Do
            End Do
          End Do
        End If


!-----------------------------------------------------------------------
! 1.2 Call deep convection code
!-----------------------------------------------------------------------
! DEPENDS ON: deep_conv
      Call DEEP_CONV(                                                   &
                       !IN
                       nbl,nlev,ntra,n_cca_lev,n_dp,trlev,              &
                       bland_dp, delthvu_dp, exner_layer_centres_dp,    &
                       exner_layer_boundaries_dp,flg_up_flx,            &
                       flg_up_flx_half,                                 &
                       flg_entr_up, flg_detr_up, flg_dwn_flx,           &
                       flg_entr_dwn,flg_detr_dwn,flg_uw_deep,           &
                       flg_vw_deep, L_calc_dxek, L_q_interact,          &
                       L_tracer, ntml_dp, ntpar_dp,                     &
                       pstar_dp,p_layer_centres_dp,                     &
                       p_layer_boundaries_dp,                           &
                       z_theta_dp, z_rho_dp, r_theta_dp, r_rho_dp,      &
                       rho_theta_dp,rho_dp,r2rho_th_dp,r2rho_dp,        &
                       dr_across_th_dp,dr_across_rh_dp,                 &
                       q_dp,q1_sd_dp,                                   &
                       t1_sd_dp,th_dp,timestep,u_dp,v_dp,               &
                       uw0_dp, vw0_dp, W_MAX_dp, wstar_dp,              &
                       zlcl_uv_dp ,freeze_lev_dp,                       &
                       recip_pstar_dp,qse_dp, dp_on,mdet_dp_on,         &
                       dp_ent_on, dp_sdet_on, dp_new_termc,             &
                       !INOUT
                       bulk_cf_dp,cf_frozen_dp,cf_liquid_dp,            &
                       qcf_dp,qcl_dp,tracer_dp,                         &
                       !OUT
                       cape_out_dp,cclwp_dp,ccw_dp,                     &
                       dbcfbydt_dp,dcffbydt_dp,dcflbydt_dp,             &
                       dqbydt_dp,dqcfbydt_dp,dqclbydt_dp,               &
                       dthbydt_dp,dubydt_dp,dvbydt_dp,dtrabydt_dp,      &
                       detrain_up_dp,detrain_dwn_dp,entrain_up_dp,      &
                       entrain_dwn_dp,                                  &
                       iccb_dp,icct_dp,                                 &
                       lcca_dp,lcclwp_dp,lcbase_dp,lctop_dp,            &
                       rain_dp, snow_dp, rain_3d_dp, snow_3d_dp,        &
                       up_flux_dp, up_flux_half_dp,                     &
                       dwn_flux_dp,uw_deep_dp,vw_deep_dp,kterm_dp,      &
                       tcw_dp,cca_2d_dp,                                &
                       rbuoy_p_out_dp,the_out_dp,thp_out_dp,            &
                       qe_out_dp,qp_out_dp                              &
                       )

!-----------------------------------------------------------------------
! 1.3 Write data from deep convection points to full arrays
!-----------------------------------------------------------------------

      Do i = 1,n_dp
        cape_out(dpi(i))           = cape_out_dp(i)
        cclwp(dpi(i))              = cclwp_dp(i)
        iccb(dpi(i))               = iccb_dp(i)
        icct(dpi(i))               = icct_dp(i)
        lcca(dpi(i))               = lcca_dp(i)
        lcbase(dpi(i))             = lcbase_dp(i)
        lctop(dpi(i))              = lctop_dp(i)
        lcclwp(dpi(i))             = lcclwp_dp(i)
        rain(dpi(i))               = rain_dp(i)
        snow(dpi(i))               = snow_dp(i)
        precip_deep(dpi(i))        = rain_dp(i) + snow_dp(i)
        kterm_deep(dpi(i))         = kterm_dp(i)
      End Do

      !=====================================================================
      ! CCRad
      !=====================================================================
      ! tcw_md & cca_2d_md are passed to MID_CONV, any extra cloud from
      ! MID_CONV is added to tcw_md & cca_2d_md.  If this is used, (in the
      ! case of mid-level cld appearing above deep cld) condensed water from
      ! deep_cloud may affect higher mid-level cloud because it is included
      ! in tcw_md & cca_2d_md.
      !
      ! So zero mid-level and do not overwrite mid-level cloud variables
      ! before entering MID-CONV to remove this interdependence

      If (.NOT. lcv_ccrad) Then
        Do i=1, n_dp
          tcw_md    (dpi(i)) = tcw_dp   (i)
          iccb_md   (dpi(i)) = iccb_dp  (i)
          icct_md   (dpi(i)) = icct_dp  (i)
          cca_2d_md (dpi(i)) = cca_2d_dp(i)
        End Do
      End If ! lcv_ccrad


      If (L_mom) then
        Do k=1,nlev+1
          Do i = 1,n_dp
            dubydt(dpi(i),k)         = dubydt_dp(i,k)
            dvbydt(dpi(i),k)         = dvbydt_dp(i,k)
          End Do
        End Do
      End If

      Do k = 1,nlev
        Do i = 1,n_dp
          dthbydt(dpi(i),k)        = dthbydt_dp(i,k)
          dcflbydt(dpi(i),k)       = dcflbydt_dp(i,k)
          dcffbydt(dpi(i),k)       = dcffbydt_dp(i,k)
          dbcfbydt(dpi(i),k)       = dbcfbydt_dp(i,k)
          ccw(dpi(i),k)            = ccw_dp(i,k)
          rbuoy_p_out(dpi(i),k)      = rbuoy_p_out_dp(i,k)
          the_out(dpi(i),k)        = the_out_dp(i,k)
          thp_out(dpi(i),k)        = thp_out_dp(i,k)
! output remains in mixing/specific as scheme not converted
          qe_out(dpi(i),k)         = qe_out_dp(i,k)
          qp_out(dpi(i),k)         = qp_out_dp(i,k)
        End Do
      End Do

!
! G-R outputs specific humidity
!
      If (l_mixing_ratio) Then  ! Requires conversion

        If (l_q_interact) Then  ! PC2

          Do k = 1,nlev
            Do i = 1,n_dp
              dqtt = dqbydt_dp(i,k)+dqclbydt_dp(i,k)+dqcfbydt_dp(i,k)
              dqt  = dqtt*timestep
              denom = 1.0/((1.0-qt(dpi(i),k))*(1.0-qt(dpi(i),k)-dqt))
              dqbydt(dpi(i),k)   =  denom *                               &
                     ( dqbydt_dp(i,k)*(1.0-qt(dpi(i),k))+q_dp(i,k)*dqtt )
              dqclbydt(dpi(i),k) =  denom *                               &
                   ( dqclbydt_dp(i,k)*(1.0-qt(dpi(i),k))+qcl_dp(i,k)*dqtt )
              dqcfbydt(dpi(i),k) =  denom *                               &
                   ( dqcfbydt_dp(i,k)*(1.0-qt(dpi(i),k))+qcf_dp(i,k)*dqtt )
            End Do
          End Do

        Else                    ! Not PC2
                                ! No qcl and qcf increments anyway
          Do k = 1,nlev
            Do i = 1,n_dp
              dqt   = dqbydt_dp(i,k)*timestep
              denom = 1.0/((1.0-qt(dpi(i),k))*(1.0-qt(dpi(i),k)-dqt))
              dqbydt(dpi(i),k)    =  dqbydt_dp(i,k) *denom
              dqclbydt(dpi(i),k)  = 0.0
              dqcfbydt(dpi(i),k)  = 0.0
            End Do
          End Do

        End If                  ! test on PC2

      Else        ! output is specific humidity therefore no problems

        Do k = 1,nlev
          Do i = 1,n_dp
            dqbydt(dpi(i),k)   = dqbydt_dp(i,k)
            dqclbydt(dpi(i),k) = dqclbydt_dp(i,k)
            dqcfbydt(dpi(i),k) = dqcfbydt_dp(i,k)
          End Do
        End Do

      End If      ! Test on l_mixing_ratio


      If (flg_up_flx) then
        Do k = 1,nlev
          Do i = 1,n_dp
            up_flux(dpi(i),k)        = up_flux_dp(i,k)
          End Do
        End Do
      End If
      If (flg_up_flx_half) then
        Do k = 1,nlev
          Do i = 1,n_dp
            up_flux_half(dpi(i),k)        = up_flux_half_dp(i,k)
          End Do
        End Do
      End If
      If (flg_dwn_flx) then
        Do k = 1,nlev
          Do i = 1,n_dp
            dwn_flux(dpi(i),k)       = dwn_flux_dp(i,k)
          End Do
        End Do
      End If
      If (flg_entr_up) then
        Do k = 1,nlev
          Do i = 1,n_dp
            entrain_up(dpi(i),k)     = entrain_up_dp(i,k)
          End Do
        End Do
      End If
      If (flg_detr_up) then
        Do k = 1,nlev
          Do i = 1,n_dp
            detrain_up(dpi(i),k)     = detrain_up_dp(i,k)
          End Do
        End Do
      End If
      If (flg_entr_dwn) then
        Do k = 1,nlev
          Do i = 1,n_dp
            entrain_dwn(dpi(i),k)    = entrain_dwn_dp(i,k)
          End Do
        End Do
      End If
      If (flg_detr_dwn) then
        Do k = 1,nlev
          Do i = 1,n_dp
            detrain_dwn(dpi(i),k)    = detrain_dwn_dp(i,k)
          End Do
        End Do
      End If
      If (flg_uw_deep) then
        Do k = 1,nlev
          Do i = 1,n_dp
            uw_deep(dpi(i),k)        = uw_deep_dp(i,k)
          End Do
        End Do
      End If
      If (flg_vw_deep) then
        Do k = 1,nlev
          Do i = 1,n_dp
            vw_deep(dpi(i),k)        = vw_deep_dp(i,k)
          End Do
        End Do
      End If

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,nlev
            Do i = 1, n_dp
              dtrabydt(dpi(i),k,ktra)  = dtrabydt_dp(i,k,ktra)
            End Do
          End Do
        End Do
      End If

      ! Merge dp 3d rain & snow profiles

      Do k=1, nlev
        Do i=1, n_dp
          rain_3d(dpi(i),k) = rain_3d_dp(i,k)
          snow_3d(dpi(i),k) = snow_3d_dp(i,k)
        End Do
      End Do         

!-----------------------------------------------------------------------------
! 1.4 CCRad - Calculate CCA fow Deep levels only
!-----------------------------------------------------------------------------

      If (Lcv_ccrad) Then

        !-----------------------------------------------------------------
        ! 1.41 Calculate CCA_2D of Deep Cloud
        !-----------------------------------------------------------------

        Select Case (cca2d_dp_opt)
          Case(cca2d_srf_precip_rate)
            Do i=1, n_dp
              If (iccb_dp(i) /= 0) Then ! Deep convection was successful

                ! Use determination of CCA_2D_dp based on surface precip. rate
                ! from deep cloud.
 
                ! NOTE: at present the a_ and b_ parameters for LAND and SEA
                !       are equal, so there will be no difference between
                !       land and sea points.
 

                If (bland_dp(i)) Then ! LAND POINT
                  If ((rain_dp(i) + snow_dp(i)) > 0.0) Then
                    cca_2d_dp(i)  = a_land                                    &
                                  + b_land                                    &
                                  * ALOG(86400.0                              &
                                  * (rain_dp(i)+snow_dp(i)))
                  End If
                Else                  ! SEA POINT
                  If ((rain_dp(i) + snow_dp(i)) > 0.0) Then
                    cca_2d_dp(i)  = a_sea                                     &
                                  + b_sea                                     &
                                  * ALOG(86400.0                              &
                                  * (rain_dp(i)+snow_dp(i)))
                  End If
                End If      ! bland_dp
              End If      ! iccb_dp
            End Do      ! i (n_dp)

          Case (cca2d_total_condensed_water)
            ! cca_2d_dp left unchanged from code, which is based on
            ! TCW (Total Condensed Water) (This is a rate)

        End Select


        Do i=1, n_dp
          !-------------------------------------          
          ! Apply Deep CCA Tuning factor
          !-------------------------------------
          cca_2d_dp(i) = cca_2d_dp(i)*cca_dp_knob

          !-------------------------------------
          ! Make sure cca_2d_dp within limits
          !-------------------------------------
          cca_2d_dp(i) = MAX(0.0,    cca_2d_dp(i))
          cca_2d_dp(i) = MIN(1.0E+0, cca_2d_dp(i))

        End Do      ! i (n_dp)


        !---------------------------------------------------------------------
        ! 1.42 Apply CCA_2D to 3d cloud profile
        !---------------------------------------------------------------------

        If (Lcv_3d_cca) Then

          ! Apply anvil scheme to deep cloud
          Call calc_3d_cca(n_dp, n_dp, nlev, nbl, iccb_dp, icct_dp            &
             , p_layer_boundaries_dp, freeze_lev_dp, cca_2d_dp, cca_dp        &
             , z_theta_dp, z_rho_dp, l_q_interact, .false., l_pc2_diag_sh_pts)

          ! NOTE: iccb_dp, icct_dp are layer centres (theta levels) at this
          !        point.

        Else

          ! Apply cca_2d_dp to all levels from deep base to deep top
          Do i=1, n_dp
            Do k=iccb_dp(i), icct_dp(i)
              cca_dp(i,k) = cca_2d_dp(i)
            End Do
          End Do

        End If      ! Lcv_3d_cca

        !---------------------------------------------------------------------
        ! 1.43 Apply ccw_dp_knob
        !---------------------------------------------------------------------
        Do k=1, nlev
          Do i=1, n_dp
            ccw_dp(i,k) = ccw_dp(i,k)*ccw_dp_knob
          End Do
        End Do


        !---------------------------------------------------------------------
        ! 1.44 Apply CCA_dp/CCW_dp to full field ..... where ccw > 0.0
        !---------------------------------------------------------------------

        Do k=1, nlev
          Do i=1, n_dp
            If (iccb_dp(i) /= 0) Then ! Deep conv. was successful

              cca(dpi(i),k) = cca_dp(i,k)
              ccw(dpi(i),k) = ccw_dp(i,k)

            End If      ! iccb_dp
          End Do      ! i
        End Do      ! k

      End If      ! lcv_ccrad

      End If !n_dp > 0


!-----------------------------------------------------------------------
! 2.0 SHALLOW convection
! 2.1 Compress input variable arrays for shallow convection scheme to
!     length n_sh (shallow points only)
!-----------------------------------------------------------------------
      If (n_sh  >   0) then
        j = 0
        Do i = 1,npnts
          If (cumulus_bl(i).and.l_shallow_bl(i)) then
            j                        = j+1
            shi(j)                   = i
          End If
        End Do
!
! In only variables
!
        Do j = 1,n_sh

          bland_sh(j)           = bland(shi(j))
          delthvu_sh(j)         = delthvu_bl(shi(j))
          ntml_sh(j)            = ntml(shi(j))
          ntpar_sh(j)           = ntpar(shi(j))
          pstar_sh(j)           = pstar(shi(j))
          recip_pstar_sh(j)     = recip_pstar(shi(j))
          q1_sd_sh(j)           = q1_sd(shi(j))
          t1_sd_sh(j)           = t1_sd(shi(j))
          uw0_sh(j)             = uw0(shi(j))
          vw0_sh(j)             = vw0(shi(j))
          wstar_sh(j)           = wstar(shi(j))
          wthvs_sh(j)           = wthvs(shi(j))
          zlcl_uv_sh(j)         = zlcl_uv(shi(j))
          ztop_uv_sh(j)         = ztop_uv(shi(j))
          freeze_lev_sh(j)      = freeze_lev(shi(j))
        End Do

        Do k = 0,nlev
          Do j = 1,n_sh
            p_layer_centres_sh(j,k)                                     &
     &                            = p_layer_centres(shi(j),k)
            p_layer_boundaries_sh(j,k)                                  &
     &                            = p_layer_boundaries(shi(j),k)
            exner_layer_centres_sh(j,k)                                 &
     &                            = exner_layer_centres(shi(j),k)
            exner_layer_boundaries_sh(j,k)                              &
     &                            = exner_layer_boundaries(shi(j),k)
          End Do
        End Do

        Do k = 1,nlev
          Do j = 1,n_sh
            u_sh(j,k)           = u(shi(j),k)
            v_sh(j,k)           = v(shi(j),k)
            th_sh(j,k)          = th(shi(j),k)
            z_theta_sh(j,k)     = z_theta(shi(j),k)
            z_rho_sh(j,k)       = z_rho(shi(j),k)
          End Do
        End Do

! Moisture - scheme requires input of specific humidity

        If (l_mixing_ratio) Then
!  Conversion required
          Do k = 1,nlev
            Do j = 1,n_sh
              q_sh(j,k)   = q(shi(j),k)  /(1.0+qt(shi(j),k))
              qse_sh(j,k) = qse(shi(j),k)/(1.0+qse(shi(j),k))
              qcl_sh(j,k) = qcl(shi(j),k)/(1.0+qt(shi(j),k))
              qcf_sh(j,k) = qcf(shi(j),k)/(1.0+qt(shi(j),k))
            End Do
          End Do

        Else     ! Input is specific humidity therefore no problems

          Do k = 1,nlev
            Do j = 1,n_sh
              q_sh(j,k)   = q(shi(j),k)
              qse_sh(j,k) = qse(shi(j),k)
              qcl_sh(j,k) = qcl(shi(j),k)
              qcf_sh(j,k) = qcf(shi(j),k)
            End Do
          End Do

        End If   ! Test on l_mixing_ratio
!
! In/out variables
!

        Do k = 1,nlev
          Do j = 1,n_sh
            cf_liquid_sh(j,k)   = cf_liquid(shi(j),k)
            cf_frozen_sh(j,k)   = cf_frozen(shi(j),k)
            bulk_cf_sh(j,k)     = bulk_cf(shi(j),k)
          End Do
        End Do

        If (L_tracer) then
          Do ktra = 1,ntra
            Do k = 1,trlev
              Do j = 1,n_sh
                tracer_sh(j,k,ktra)  = tracer(shi(j),k,ktra)
              End Do
            End Do
            Do k = 1,nlev
              Do j = 1,n_sh
                dtrabydt_sh(j,k,ktra)  = 0.0
              End Do
            End Do
          End Do
        End If


!-----------------------------------------------------------------------
! 2.2 Call shallow convection code
!-----------------------------------------------------------------------
! DEPENDS ON: shallow_conv
       Call SHALLOW_CONV(                                               &
                           !IN
     &                     nbl,nlev,ntra,n_cca_lev,n_sh,trlev,          &
     &                     bland_sh,delthvu_sh,exner_layer_centres_sh,  &
     &                     exner_layer_boundaries_sh,                   &
     &                     flg_up_flx,flg_up_flx_half,                  &
     &                     flg_entr_up, flg_detr_up, flg_dwn_flx,       &
     &                     flg_entr_dwn,flg_detr_dwn,flg_uw_shall,      &
     &                     flg_vw_shall, L_calc_dxek,                   &
     &                     L_q_interact, L_tracer, ntml_sh, ntpar_sh,   &
     &                     pstar_sh,p_layer_centres_sh,                 &
     &                     p_layer_boundaries_sh,z_theta_sh,z_rho_sh,   &
     &                     r2rho_th_sh, dr_across_th_sh,                &
     &                     q_sh,q1_sd_sh,t1_sd_sh,th_sh,timestep,       &
     &                     u_sh,v_sh,uw0_sh,vw0_sh,wstar_sh,wthvs_sh,   &
     &                     zlcl_uv_sh,ztop_uv_sh,freeze_lev_sh,         &
     &                     recip_pstar_sh,qse_sh, sh_on,mdet_sh_on,     &
     &                     sh_ent_on, sh_sdet_on, sh_new_termc,         &
                           !INOUT
     &                     bulk_cf_sh,cf_frozen_sh,cf_liquid_sh,        &
     &                     qcf_sh,qcl_sh,tracer_sh,                     &
                           !OUT
     &                     cape_out_sh,cclwp_sh,ccw_sh,                 &
     &                     dbcfbydt_sh,dcffbydt_sh,dcflbydt_sh,         &
     &                     dqbydt_sh,dqcfbydt_sh,dqclbydt_sh,           &
     &                     dthbydt_sh,dubydt_sh,dvbydt_sh,dtrabydt_sh,  &
     &                     detrain_up_sh,detrain_dwn_sh,entrain_up_sh,  &
     &                     entrain_dwn_sh,                              &
     &                     iccb_sh,icct_sh,                             &
     &                     lcca_sh,lcclwp_sh,lcbase_sh,lctop_sh,        &
     &                     rain_sh, snow_sh, rain_3d_sh, snow_3d_sh,    &
     &                     up_flux_sh, up_flux_half_sh,                 &
     &                     dwn_flux_sh,uw_shall_sh,vw_shall_sh,         &
     &                     tcw_sh,cca_2d_sh                             &
     &                     ,rbuoy_p_out_sh,the_out_sh,thp_out_sh        &
     &                     ,qe_out_sh,qp_out_sh                         &
     &                     )

!-----------------------------------------------------------------------
! 2.3 Write data from shallow convection points to full arrays
!-----------------------------------------------------------------------

      Do i = 1,n_sh
        cape_out(shi(i))           = cape_out_sh(i)
        cclwp(shi(i))              = cclwp_sh(i)
        iccb(shi(i))               = iccb_sh(i)
        icct(shi(i))               = icct_sh(i)
        lcca(shi(i))               = lcca_sh(i)
        lcbase(shi(i))             = lcbase_sh(i)
        lctop(shi(i))              = lctop_sh(i)
        lcclwp(shi(i))             = lcclwp_sh(i)
        rain(shi(i))               = rain_sh(i)
        snow(shi(i))               = snow_sh(i)
        precip_shall(shi(i))       = rain_sh(i) + snow_sh(i)
      End Do

      !=====================================================================
      ! CCRAD
      !=====================================================================
      ! tcw_md & cca_2d_md are passed to MID_CONV, any extra cloud from
      ! MID_CONV is added to tcw_md & cca_2d_md.  If this is used, in the
      ! case of mid-level cld appearing above shallow cld, condensed water
      ! from shallow_cloud may affect higher mid-level cloud because it is
      ! included in tcw_md & cca_2d_md.

      ! So zero mid-level and do not overwrite mid-level cloud variables
      ! before entering MID-CONV to remove this interdependence

      If (.NOT. Lcv_ccrad) Then
        Do i=1, n_sh          
          tcw_md    (shi(i)) = tcw_sh    (i)
          iccb_md   (shi(i)) = iccb_sh   (i)
          icct_md   (shi(i)) = icct_sh   (i)
          cca_2d_md (shi(i)) = cca_2d_sh (i)
        End Do
      End If


      If (L_mom) then
        Do k=1,nlev+1
          Do i = 1,n_sh
            dubydt(shi(i),k)         = dubydt_sh(i,k)
            dvbydt(shi(i),k)         = dvbydt_sh(i,k)
          End Do
        End Do
      End If

      Do k = 1,nlev
        Do i = 1,n_sh
          dthbydt(shi(i),k)        = dthbydt_sh(i,k)
          dcflbydt(shi(i),k)       = dcflbydt_sh(i,k)
          dcffbydt(shi(i),k)       = dcffbydt_sh(i,k)
          dbcfbydt(shi(i),k)       = dbcfbydt_sh(i,k)
          ccw(shi(i),k)            = ccw_sh(i,k)
          rbuoy_p_out(shi(i),k)      = rbuoy_p_out_sh(i,k)
          the_out(shi(i),k)        = the_out_sh(i,k)
          thp_out(shi(i),k)        = thp_out_sh(i,k)
          qe_out(shi(i),k)         = qe_out_sh(i,k)
          qp_out(shi(i),k)         = qp_out_sh(i,k)
        End Do
      End Do

!
! G-R requires input of specific humidity
!
      If (l_mixing_ratio) Then ! requires conversion

       If (l_q_interact) Then  ! PC2

        Do k = 1,nlev
          Do i = 1,n_sh
            dqtt  = dqbydt_sh(i,k)+dqclbydt_sh(i,k)+dqcfbydt_sh(i,k)
            dqt   = dqtt*timestep
            denom = 1.0/((1.0-qt(shi(i),k))*(1.0-qt(shi(i),k)-dqt))
            dqbydt(shi(i),k)   = denom *                                   &
                     ( dqbydt_sh(i,k)*(1.0-qt(shi(i),k))+q_sh(i,k)*dqtt )
            dqclbydt(shi(i),k) = denom *                                   &
                   ( dqclbydt_sh(i,k)*(1.0-qt(shi(i),k))+qcl_sh(i,k)*dqtt )
            dqcfbydt(shi(i),k) = denom *                                   &
                   ( dqcfbydt_sh(i,k)*(1.0-qt(shi(i),k))+qcf_sh(i,k)*dqtt )
          End Do
        End Do

       Else                    ! Not PC2
                               ! No qcl and qcf increments anyway
        Do k = 1,nlev
          Do i = 1,n_sh
            dqt   = dqbydt_sh(i,k)*timestep
            denom = 1.0/((1.0-qt(shi(i),k))*(1.0-qt(shi(i),k)-dqt))
            dqbydt(shi(i),k)    = dqbydt_sh(i,k) *denom
            dqclbydt(shi(i),k)  = 0.0
            dqcfbydt(shi(i),k)  = 0.0
          End Do
        End do

       End If                  ! end test on PC2

      Else   ! output is specific humidity therefore no problems

        Do k = 1,nlev
          Do i = 1,n_sh
            dqbydt(shi(i),k)   = dqbydt_sh(i,k)
            dqclbydt(shi(i),k) = dqclbydt_sh(i,k)
            dqcfbydt(shi(i),k) = dqcfbydt_sh(i,k)
          End Do
        End Do

      End If      ! Test on l_mixing_ratio

      If (flg_up_flx) then
        Do k = 1,nlev
          Do i = 1,n_sh
            up_flux(shi(i),k)        = up_flux_sh(i,k)
          End Do
        End Do
      End If
      If (flg_up_flx_half) then
        Do k = 1,nlev
          Do i = 1,n_sh
            up_flux_half(shi(i),k)   = up_flux_half_sh(i,k)
          End Do
        End Do
      End If
      If (flg_dwn_flx) then
        Do k = 1,nlev
          Do i = 1,n_sh
            dwn_flux(shi(i),k)       = dwn_flux_sh(i,k)
          End Do
        End Do
      End If
      If (flg_entr_up) then
         Do k = 1,nlev
          Do i = 1,n_sh
           entrain_up(shi(i),k)     = entrain_up_sh(i,k)
          End Do
        End Do
      End If
      If (flg_detr_up) then
        Do k = 1,nlev
          Do i = 1,n_sh
            detrain_up(shi(i),k)     = detrain_up_sh(i,k)
          End Do
        End Do
      End If
      If (flg_entr_dwn) then
        Do k = 1,nlev
          Do i = 1,n_sh
            entrain_dwn(shi(i),k)    = entrain_dwn_sh(i,k)
          End Do
        End Do
      End If
      If (flg_detr_dwn) then
        Do k = 1,nlev
          Do i = 1,n_sh
            detrain_dwn(shi(i),k)    = detrain_dwn_sh(i,k)
          End Do
        End Do
      End If
      If (flg_uw_shall) then
        Do k = 1,nlev
          Do i = 1,n_sh
            uw_shall(shi(i),k)       = uw_shall_sh(i,k)
          End Do
        End Do
      End If
      If (flg_vw_shall) then
        Do k = 1,nlev
          Do i = 1,n_sh
            vw_shall(shi(i),k)       = vw_shall_sh(i,k)
          End Do
        End Do
      End If

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,nlev
            Do i = 1,n_sh
              dtrabydt(shi(i),k,ktra)  = dtrabydt_sh(i,k,ktra)
            End Do
          End Do
        End Do
      End If

      ! Merge sh 3d rain & snow profiles

      Do k=1,nlev
        Do i=1,n_sh
          rain_3d(shi(i),k) = rain_3d_sh(i,k)
          snow_3d(shi(i),k) = snow_3d_sh(i,k)
        End Do
      End Do

!-----------------------------------------------------------------------------
! 2.4 CCRad - Calculate CCA fow shallow levels only
!-----------------------------------------------------------------------------

        If (lcv_ccrad) Then

          Do i=1, n_sh

            overlap_fac(i) = 0.0
            zcld(i)        = 0.0

            If (iccb_sh(i) /= 0) Then ! Shallow convection occured

              !---------------------------------------------------------------
              ! Grant and Lock (2004) LES show mb/wsc nicely scales the cloud
              ! fraction profiles but not the TCA.  Also the UM overlap
              ! assumption in radiation is maximal.  This implies significant
              ! underestimate of TCA. So include a further parametrization of
              ! Cu "overlap", based again on the LES of Grant and Lock (2004).
              ! This increases cca_2d proportional to the ratio of the cloud
              ! to sub-cloud layer depths.  In order to preserve the grid-box
              ! cloud water, ccw will be divided by the same factor.
              !---------------------------------------------------------------

              zcld(i) = z_rho_sh(i,ntpar_sh(i)+1) -                           &
                        z_rho_sh(i, ntml_sh(i)+1)

              overlap_fac(i) = 2.0*zcld(i)                                    &
                             / z_rho_sh(i,ntml_sh(i)+1)
            End If     ! iccb_sh
          End Do     ! n_sh

          !---------------------------------------------------------------
          ! 2.41 Calculate CCA
          !---------------------------------------------------------------
          Select Case (cca2d_sh_opt)
             
            Case(cca2d_grant_lock)
              Do i=1, n_sh
                If (iccb_sh(i) /= 0) Then ! Shallow convection occured
                  mb  = 0.03*wstar_sh(i)

                  wsc = (  delthvu_sh(i)*mb*g                                 &
                      / (th_sh(i,ntml_sh(i))*(1.0+0.61*q_sh(i,ntml_sh(i)))))  &
                      ** 0.3333

                  cca_2d_sh(i) = 2.0*mb/wsc

                  mb  = 0.0
                  wsc = 0.0
                End If     ! iccb_sh
              End Do     ! n_sh

            Case(cca2d_total_condensed_water)
              ! cca_2d_sh is left unchanged from that calculated in the
              ! code, which is based on TCW (Total Condensed Water)
              ! (TCW is a rate) 

          End Select


          Do i=1, n_sh
            !---------------------------------------
            ! Apply Shallow CCA Tuning factor
            !---------------------------------------
            cca_2d_sh(i) = cca_sh_knob*cca_2d_sh(i)

            !---------------------------------------
            ! Make sure cca_2d_sh remains in limits
            !---------------------------------------
            cca_2d_sh(i) = MAX(0.0,    cca_2d_sh(i))
            cca_2d_sh(i) = MIN(1.0E+0, cca_2d_sh(i))

            overlap_fac(i) = MAX( 0.5, overlap_fac(i) )
            overlap_fac(i) = MIN( 5.0, overlap_fac(i) )

            If (overlap_fac(i)*cca_2d_sh(i) > 0.99) Then
              overlap_fac(i) = 0.99/cca_2d_sh(i)
            End If
          End Do      ! i (n_sh)


          !-------------------------------------------------------------------
          ! 2.42 Fill cca_sh with cca_2d_sh where non-zero ccw_sh
          !-------------------------------------------------------------------
          Do k=1, nlev
            Do i=1, n_sh
              If (iccb_sh(i) /= 0) Then ! Shallow convection occured

                ! Apply ccw_sh_knob
                ccw_sh(i,k) = ccw_sh_knob*ccw_sh(i,k)

                If (ccw_sh(i,k) > 0.0) Then
                  zpr = (z_rho_sh(i,k) -                                      &
                         z_rho_sh(i,ntml_sh(i)+1))                            &
                      /  zcld(i)

                  ! Apply Shape-function
                  !
                  ! Apply overlap_fac to cca, also preserving grid-box water
                  ! by dividing ccw by overlap_fac, at least at cloud-base

                  ccw_sh(i,k)  = ccw_sh(i,k)/overlap_fac(i)
                  zpr          = MIN(1.0,zpr)
                  cca_sh(i,k)  = overlap_fac(i)*cca_2d_sh(i)                  &
                               * 0.25*( 1.0 + 3.0*exp(-5.0*zpr) )

                End If       ! ccw_sh
              End If       ! iccb_sh
            End Do       ! i (n_sh)
          End Do       ! k (nlev)

          !-------------------------------------------------------------------
          ! 2.43 Apply CCA/CCW to full field .....(on shallow points)
          !-------------------------------------------------------------------

          Do k=1, nlev
            Do i=1, n_sh
              If (iccb_sh(i) /= 0) Then ! Shallow convection occured
                cca(shi(i),k) = cca_sh(i,k)
                ccw(shi(i),k) = ccw_sh(i,k)
              End If
            End Do      ! i
          End Do      ! k
        End If      ! lcv_ccrad
      End If     ! n_sh > 0

! make copies of rbuoy etc. for checking purposes
!Used after mid convection to see if mid has updated values
      Do k = 1, nlev
        Do i = 1, npnts
          rbuoy_p_out_copy(i,k) = rbuoy_p_out(i,k)
          the_out_copy(i,k) = the_out(i,k)
          thp_out_copy(i,k) = thp_out(i,k)
          qe_out_copy(i,k) = qe_out(i,k)
          qp_out_copy(i,k) = qp_out(i,k)
        End Do
      End Do
!-----------------------------------------------------------------------
! 3.0 MID-LEVEL Convection
! 3.1 Set lowest level that mid level convection can trigger
!-----------------------------------------------------------------------

      idp=0
      Do i = 1,npnts
        If (.not. cumulus_bl(i)) then
          midtrig(i) = ntml(i) + 1
          If (ntml(i)  ==  nbl) then
            midtrig(i) = ntml(i)
          End If
        else
!
!  Cumulus points ie deep or shallow convection has occurred
!
          midtrig(i) = ntpar(i) + 2

!
! NTPAR has a maximum value which can be less than the top level for
! deep. Deep convection may terminate above this. Need to prevent mid
! level convection occurring at the same levels as deep.
!
          if(.not.l_shallow_bl(i)) then   ! deep points
            idp=idp+1
            If (kterm_dp(idp) >  ntpar(i)+2) then
               midtrig(i) = kterm_dp(idp)
            End If
          End If  ! deep points
        End If
      End Do

!-----------------------------------------------------------------------
! 3.2 Copy all input arrays to arrays ending in _md for passing to
!     mid-level scheme
!-----------------------------------------------------------------------

      Do k = 0,nlev
        Do i = 1,npnts
          p_layer_centres_md(i,k)                                       &
     &                                 = p_layer_centres(i,k)
          p_layer_boundaries_md(i,k)                                    &
     &                                 = p_layer_boundaries(i,k)
          exner_layer_centres_md(i,k)                                   &
     &                                 = exner_layer_centres(i,k)
          exner_layer_boundaries_md(i,k)                                &
     &                                 = exner_layer_boundaries(i,k)
          r_theta_md(i,k)     = r_theta(i,k)
        End Do
      End Do
      Do k = 1,nlev
        Do i = 1,npnts
          u_md(i,k)           = u(i,k)
          v_md(i,k)           = v(i,k)
          th_md(i,k)          = th(i,k)
          z_theta_md(i,k)     = z_theta(i,k)
          z_rho_md(i,k)       = z_rho(i,k)
          r_rho_md(i,k)       = r_rho(i,k)
          cf_liquid_md(i,k)   = cf_liquid(i,k)
          cf_frozen_md(i,k)   = cf_frozen(i,k)
          bulk_cf_md(i,k)     = bulk_cf(i,k)
        End Do
      End Do

! Moisture - scheme requires input of specific humidity

      If (l_mixing_ratio) Then
!  Conversion required
        Do k = 1,nlev
          Do i = 1,npnts
            q_md(i,k)   = q(i,k)  /(1.0+qt(i,k))
            qse_md(j,k) = qse(j,k)/(1.0+qse(j,k))
            qcl_md(i,k) = qcl(i,k)/(1.0+qt(i,k))
            qcf_md(i,k) = qcf(i,k)/(1.0+qt(i,k))
          End Do
        End Do

      Else       ! input is specific humidity therefore no problems

        Do k = 1,nlev
          Do j = 1,npnts
            q_md(j,k)   = q(j,k)
            qse_md(j,k) = qse(j,k)
            qcl_md(j,k) = qcl(j,k)
            qcf_md(j,k) = qcf(j,k)
          End Do
        End Do

      End If     ! Test on l_mixing_ratio

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,trlev
            Do i = 1,npnts
              tracer_md(i,k,ktra)  = tracer(i,k,ktra)
            End Do
          End Do
          Do k = 1,nlev
            Do i = 1,npnts
              dtrabydt_md(i,k,ktra)  = 0.0
            End Do
          End Do
        End Do
      End If

!-----------------------------------------------------------------------
! 3.3 Call mid-level convection code
!-----------------------------------------------------------------------

! DEPENDS ON: mid_conv
      Call MID_CONV(                                                    &
                      !IN
     &                nbl,nlev,ntra,n_cca_lev,npnts,trlev,              &
     &                bland,W_MAX,exner_layer_centres_md,               &
     &                exner_layer_boundaries_md,flg_up_flx,             &
     &                flg_up_flx_half,                                  &
     &                flg_entr_up, flg_detr_up, flg_dwn_flx,            &
     &                flg_entr_dwn,flg_detr_dwn,                        &
     &                L_calc_dxek, L_q_interact,                        &
     &                L_tracer, midtrig, ntml,                          &
     &                ntpar,pstar,p_layer_centres_md,                   &
     &                p_layer_boundaries_md,                            &
     &                r_theta_md, r_rho_md,                             &
     &                z_theta_md, z_rho_md,rho,rho_theta,q_md,q1_sd,    &
     &                t1_sd,th_md,timestep,u_md,v_md,recip_pstar,qse,   &
     &                md_on, mdet_md_on, md_ent_on, md_sdet_on,         &
     &                md_new_termc,                                     &
                      !INOUT
     &                bulk_cf_md,cf_frozen_md,cf_liquid_md,             &
     &                qcf_md,qcl_md,tracer_md,                          &
                      !OUT
     &                cape_out_md,cclwp_md,ccw_md,                      &
     &                dbcfbydt_md,dcffbydt_md,dcflbydt_md,              &
     &                dqbydt_md,dqcfbydt_md,dqclbydt_md,                &
     &                dthbydt_md,dubydt_md,dvbydt_md,dtrabydt_md,       &
     &                detrain_up_md,detrain_dwn_md,entrain_up_md,       &
     &                entrain_dwn_md,                                   &
     &                iccb_md,icct_md,                                  &
     &                lcca_md,lcclwp_md,lcbase_md,lctop_md,             &
     &                rain_md, snow_md, rain_3d_md, snow_3d_md,         &
     &                up_flux_md, up_flux_half_md,                      &
     &                dwn_flux_md,tcw_md,l_mid_all,cca_2d_md            &
     &                ,rbuoy_p_out_md,the_out_md,thp_out_md             &
     &                ,qe_out_md,qp_out_md                              &
     &                ,uw_mid, vw_mid                                   &
     &                )



!-----------------------------------------------------------------------
! 3.4 Write data from mid-level convection to full arrays
!-----------------------------------------------------------------------

      Do i = 1,npnts
!
! Cloud variables - only overwrite deep or shallow values if
! iccb_md  and icct_md > 0
!
        If (iccb_md(i)  >   0 .and. icct_md(i)  >   0) then
          iccb(i)                    = iccb_md(i)
          icct(i)                    = icct_md(i)
        End If
!
! Overwrite lowest cloud values only if cumulus = .F.
!
       If (.not. cumulus_bl(i)) then
          lcca(i)                    = lcca_md(i)
          lcbase(i)                  = lcbase_md(i)
          lctop(i)                   = lctop_md(i)
          lcclwp(i)                  = lcclwp_md(i)
       End If
!
! Write remaining data to full arrays
!
        cape_out(i)                = cape_out(i) + cape_out_md(i)
        cclwp(i)                   = cclwp(i)    + cclwp_md(i)
        rain(i)                    = rain(i)     + rain_md(i)
        snow(i)                    = snow(i)     + snow_md(i)
        precip_mid(i)              = rain_md(i)  + snow_md(i)
        cca_2d(i)                  = cca_2d_md(i)

        !=====================================================================
        ! NOTE: At this point if Lcv_ccrad = T, then cca_2d_md is ONLY equal to
        !       that from mid-level cloud on a given grid point.  The original
        !       code I.E. Lcv_ccrad = F means that cca_2d_md will include that
        !       from sh/dp aswell.
        !=====================================================================

      End Do

      ! Merge md 3d rain & snow profiles
      Do k=1,nlev
        Do i=1,npnts
          rain_3d(i,k) = rain_3d(i,k) + rain_3d_md(i,k)
          snow_3d(i,k) = snow_3d(i,k) + snow_3d_md(i,k)
        End Do
      End Do


      If (L_mom) then
        Do k=1,nlev+1
          Do i = 1,npnts
            dubydt(i,k)            = dubydt(i,k) + dubydt_md(i,k)
            dvbydt(i,k)            = dvbydt(i,k) + dvbydt_md(i,k)
          End Do
        End Do
      End If

      If (.NOT. lcv_ccrad) Then
        Do k = 1,nlev
          Do i = 1,npnts
            ccw(i,k)               = ccw(i,k)      + ccw_md(i,k)
          End Do
        End Do
      End If

      Do k = 1,nlev
        Do i = 1,npnts
          dthbydt(i,k)             = dthbydt(i,k)  + dthbydt_md(i,k)
          dcflbydt(i,k)            = dcflbydt(i,k) + dcflbydt_md(i,k)
          dcffbydt(i,k)            = dcffbydt(i,k) + dcffbydt_md(i,k)
          dbcfbydt(i,k)            = dbcfbydt(i,k) + dbcfbydt_md(i,k)




! with mid level arrays from parcel, check that values
! have changed from the values they were initialised to
! so that we only overwrite at points where mid level convection
! has taken place and don't inadvertently over-write shallow and
! deep values with the original model profiles.
!If test should be on parcel values, 'cos environment values
!shouldn't change
!looking for points where shallow and deep didn't affect original
!profile, but mid did
          If ((rbuoy_p_out_copy(i,k)  ==  0.0)                          &
     &        .and. (rbuoy_p_out_md(i,k)  /=  0.0)) then
            rbuoy_p_out(i,k) = rbuoy_p_out_md(i,k)
          End If

          If ((thp_out_copy(i,k)  ==  th(i,k)) .AND.                    &
     &        (thp_out_md(i,k)  /=  th(i,k))) then
            the_out(i,k) = the_out_md(i,k)
            thp_out(i,k) = thp_out_md(i,k)
          End If
          If ((qp_out_copy(i,k)  ==  q(i,k)) .and.                      &
     &        (qp_out_md(i,k)  /=  q(i,k))) then
            qe_out(i,k) = qe_out_md(i,k)
            qp_out(i,k) = qp_out_md(i,k)
          End If

        End Do
      End Do

!
! G-R outputs specific humidity
!
      If (l_mixing_ratio) Then  ! Requires conversion

        If (l_q_interact) Then  ! PC2

          Do k = 1,nlev
            Do i = 1,npnts
              dqtt= dqbydt_md(i,k)+dqclbydt_md(i,k)+dqcfbydt_md(i,k)
              dqt = dqtt*timestep
              denom = 1.0/((1.0-qt(i,k))*(1.0-qt(i,k)-dqt))
              dqbydt(i,k)   =  dqbydt(i,k) + denom *                     &
                         ( dqbydt_md(i,k)*(1.0-qt(i,k))+q(i,k)*dqtt )
              dqclbydt(i,k) = dqclbydt(i,k) + denom *                    &
                         ( dqclbydt_md(i,k)*(1.0-qt(i,k))+qcl(i,k)*dqtt )
              dqcfbydt(i,k) = dqcfbydt(i,k) +  denom *                   &
                         ( dqcfbydt_md(i,k)*(1.0-qt(i,k))+qcf(i,k)*dqtt )
            End Do
          End Do

        Else                    ! Not PC2
                                ! No qcl and qcf increments anyway
          Do k = 1,nlev
            Do i = 1,npnts
              dqt   = dqbydt_md(i,k)*timestep
              denom = 1.0/((1.0-qt(i,k))*(1.0-qt(i,k)-dqt))
              dqbydt(i,k)    = dqbydt(i,k) + dqbydt_md(i,k)*denom
              dqclbydt(i,k)  = 0.0
              dqcfbydt(i,k)  = 0.0
            End Do
          End do

        End If                  ! test on PC2

      Else        ! output is specific humidity therefore no problems

        Do k = 1,nlev
          Do i = 1,npnts
            dqbydt(i,k)   = dqbydt(i,k)   + dqbydt_md(i,k)
            dqclbydt(i,k) = dqclbydt(i,k) + dqclbydt_md(i,k)
            dqcfbydt(i,k) = dqcfbydt(i,k) + dqcfbydt_md(i,k)
          End Do
        End Do

      End If      ! Test on l_mixing_ratio


      If (flg_up_flx) then
        Do k = 1,nlev
          Do i = 1,npnts
             up_flux(i,k)           = up_flux(i,k) + up_flux_md(i,k)
          End Do
        End Do
      End If
      If (flg_up_flx_half) then
        Do k = 1,nlev
          Do i = 1,npnts
             up_flux_half(i,k)= up_flux_half(i,k)                       &
     &             + up_flux_half_md(i,k)
          End Do
        End Do
      End If
      If (flg_dwn_flx) then
        Do k = 1,nlev
          Do i = 1,npnts
            dwn_flux(i,k)          = dwn_flux(i,k) + dwn_flux_md(i,k)
          End Do
        End Do
      End If
      If (flg_entr_up) then
        Do k = 1,nlev
          Do i = 1,npnts
            entrain_up(i,k)        = entrain_up(i,k)                    &
     &                                      + entrain_up_md(i,k)
           End Do
        End Do
      End If
      If (flg_detr_up) then
        Do k = 1,nlev
          Do i = 1,npnts
            detrain_up(i,k)        = detrain_up(i,k)                    &
     &                                      + detrain_up_md(i,k)
          End Do
        End Do
      End If
      If (flg_entr_dwn) then
        Do k = 1,nlev
          Do i = 1,npnts
            entrain_dwn(i,k)       = entrain_dwn(i,k)                   &
     &                                      + entrain_dwn_md(i,k)
          End Do
        End Do
      End If
      If (Flg_detr_dwn) then
        Do k = 1,nlev
          Do i = 1,npnts
            detrain_dwn(i,k)       = detrain_dwn(i,k)                   &
     &                                      + detrain_dwn_md(i,k)
          End Do
        End Do
      End If
        If (L_tracer) then
          Do ktra = 1,ntra
            Do k = 1,nlev
              Do i = 1,npnts
                dtrabydt(i,k,ktra)  = dtrabydt(i,k,ktra)                &
     &                                     + dtrabydt_md(i,k,ktra)
              End Do
            End Do
          End Do
        End If


!-----------------------------------------------------------------------------
! 3.5 Calculate CCA fow MID-CONVECTION levels only
!-----------------------------------------------------------------------------

      If (lcv_ccrad) Then


        n_mdcld      = 0    ! Initialise no of points with mid-cld
        n_mdcld_mult = 0    ! Initialise no of points with multiple mid-cld

        Do i=1, npnts
          If ((l_mid_all(i)) .AND. (lcbase_md(i) /= 0)) Then

            ! Mid-level cloud present
            ! Create index of gridpoints with single/multiple mid-level cloud

            If (lcbase_md(i) /= iccb_md(i)) Then

              n_mdcld_mult       = n_mdcld_mult + 1
              dum2(n_mdcld_mult) = i

            Else If (lcbase_md(i) == iccb_md(i)) Then

              n_mdcld       = n_mdcld + 1
              dum1(n_mdcld) = i
            End If

          End If
        End Do


        !-----------------------------------------------------------------
        ! Calculate CCA_2D of Mid Cloud
        !-----------------------------------------------------------------
        Select Case (cca2d_md_opt)
          Case(cca2d_srf_precip_rate)
            Do i=1, npnts
              If ((l_mid_all(i)) .AND. (lcbase_md(i) /= 0)) Then
                !-------------------------------------------------------------
                ! CCA_2D_md based on surface precipitation rate from mid cloud
                !-------------------------------------------------------------
 
                ! NOTE: at present the a_ and b_ parameters for LAND and SEA
                !       are equal, so there will be no difference between
                !       land and sea points.
  
                If (bland(i)) Then  ! LAND POINT

                  If ((rain_md(i) + snow_md(i)) > 0.0) Then
                    cca_2d_md(i)  = a_land                                    &
                                  + b_land                                    &
                                  * ALOG(24.0*3600.0                          &
                                  *(rain_md(i) + snow_md(i)))
                  End If

                Else                ! SEA-POINT

                  If ((rain_md(i) + snow_md(i)) > 0.0) Then
                    cca_2d_md(i)  = a_sea                                     &
                                  + b_sea                                     &
                                  * ALOG(24.0*3600.0                          &
                                  * (rain_md(i) + snow_md(i)))
                  End If

                End If      ! bland
              End If      ! mid-level cloud present
            End Do      ! npnts

          Case(cca2d_total_condensed_water)
            ! cca_2d_md left unchanged from code, which is based on
            ! TCW (Total Condensed Water) (This is a rate)

        End Select      ! cca2d_md_opt
  

        Do i=1, npnts
          ! Apply Mid-Level CCA Tuning factor
          !----------------------------------------
          cca_2d_md(i) = cca_md_knob*cca_2d_md(i)

          ! Make sure cca_2d_md is within limits
          !----------------------------------------
          cca_2d_md(i) = MAX(0.0,    cca_2d_md(i))
          cca_2d_md(i) = MIN(1.0E+0, cca_2d_md(i))

        End Do      ! i (npnts)


        !---------------------------------------------------------------------
        ! 3.53 Apply CCA_3D
        !---------------------------------------------------------------------
        ! A value of cca_2d has been calculated, though we need to know which
        ! levels to apply it to. So we identify which layers require CCA based
        ! on layers (above shallow/deep, if they have occurred) with non-zero
        ! ccw.
        !---------------------------------------------------------------------
 
        !=====================================================================
        ! Single mid-level cloud bank gridpoints
        !=====================================================================
 
        If (n_mdcld > 0) Then

          ! Resize compressed arrays for single Mid-level events

          ALLOCATE (mdcldi(n_mdcld))

          ! Now we know the number and indexes of gridpoints with
          ! single/multiple mid-level cloud, we can now assign the compressed
          ! mid-level arrays to pass to the anvil scheme.

          ALLOCATE (iccb_md_c       (n_mdcld       ))
          ALLOCATE (icct_md_c       (n_mdcld       ))
          ALLOCATE (freeze_lev_md_c (n_mdcld       ))
          ALLOCATE (cca_2d_md_c     (n_mdcld       ))
          ALLOCATE (cca_md_c        (n_mdcld,  nlev))
          ALLOCATE (ccw_md_c        (n_mdcld,  nlev))
          ALLOCATE (z_theta_md_c    (n_mdcld,  nlev))
          ALLOCATE (z_rho_md_c      (n_mdcld,  nlev))
          ALLOCATE (p_lyr_bnds_md_c (n_mdcld,0:nlev))

          Do i=1, n_mdcld
            mdcldi          (i) = dum1(i)
            iccb_md_c       (i) = iccb_md    (mdcldi(i))
            icct_md_c       (i) = icct_md    (mdcldi(i))
            freeze_lev_md_c (i) = freeze_lev (mdcldi(i))
            cca_2d_md_c     (i) = cca_2d_md  (mdcldi(i))
          End Do

          Do k=1, nlev
            Do i=1, n_mdcld
              cca_md_c        (i,k) = 0.0
              ccw_md_c        (i,k) = ccw_md            (mdcldi(i),k)
              z_theta_md_c    (i,k) = z_theta           (mdcldi(i),k)
              z_rho_md_c      (i,k) = z_rho             (mdcldi(i),k)
              p_lyr_bnds_md_c (i,k) = p_layer_boundaries(mdcldi(i),k)
            End Do
          End Do

          !-------------------------------------------------------------------
          ! End Resizing compressed arrays for single Mid-level events. There
          ! is only one mid-level cloud bank if a mid-event has occurred, so
          ! does not require checking to see if there are multiple mid-level
          ! cloud banks.
          !-------------------------------------------------------------------

          If (Lcv_3d_cca) Then
! DEPENDS ON: calc_3d_cca
            Call calc_3d_cca(n_mdcld, n_mdcld, nlev, nbl, iccb_md_c           &
               , icct_md_c, p_lyr_bnds_md_c, freeze_lev_md_c, cca_2d_md_c     &
               , cca_md_c, z_theta_md_c, z_rho_md_c, l_q_interact, .false.    &
               , l_pc2_diag_sh_pts)
          Else

            !-----------------------------------------------------------------
            ! Do not use anvil scheme and apply cca_2d_md_c to all cloud
            ! levels between mid-level base and top
            !-----------------------------------------------------------------
            Do k=1, nlev
              Do i=1, n_mdcld
                If (ccw_md_c(i,k) > 0) Then
                  cca_md_c(i,k) = cca_2d_md_c(i)
                Else
                  cca_md_c(i,k) = 0.0
                End If
              End Do      ! i (n_mdcld)
            End Do      ! k (nlev)
          End If     ! Lcv_3d_cca

          !-------------------------------------------------------------------
          ! Merge cca_md/ccw_md to full cca array and scale ccw_md_c by
          ! ccw_md_knob
          !-------------------------------------------------------------------

          Do k=1,nlev
            Do i=1, n_mdcld
              cca(mdcldi(i),k) = cca(mdcldi(i),k)                             &
                               + cca_md_c(i,k)
              ccw(mdcldi(i),k) = ccw(mdcldi(i),k)                             &
                               + ccw_md_c(i,k)*ccw_md_knob
            End do      ! i (n_mdcld)
          End do      ! k (nlev)

          ! Deallocate compressed single mid-level cloud arrays

          DEALLOCATE (mdcldi)
          DEALLOCATE (iccb_md_c)
          DEALLOCATE (icct_md_c)
          DEALLOCATE (freeze_lev_md_c)
          DEALLOCATE (cca_2d_md_c)
          DEALLOCATE (cca_md_c)
          DEALLOCATE (ccw_md_c)
          DEALLOCATE (z_theta_md_c)
          DEALLOCATE (z_rho_md_c)
          DEALLOCATE (p_lyr_bnds_md_c)

        End If ! n_mdcld > 0

        !=====================================================================
        ! Multiple mid-level cloud bank gridpoints
        !=====================================================================

        If (n_mdcld_mult > 0) Then
          ! Resize index array for gridpoints with multiple Mid-level cloud
          ! banks and apply indices

          ALLOCATE(mdcldi_mult(n_mdcld_mult))

          !-------------------------------------------------------------------
          ! Allocate arrays with multiple-mid level cloud gridpoints
          !-------------------------------------------------------------------
          ALLOCATE (iccb_md_mult_c       (n_mdcld_mult       ))
          ALLOCATE (icct_md_mult_c       (n_mdcld_mult       ))
          ALLOCATE (freeze_lev_md_mult_c (n_mdcld_mult       ))
          ALLOCATE (cca_2d_md_mult_c     (n_mdcld_mult       ))
          ALLOCATE (cca_md_mult_c        (n_mdcld_mult,  nlev))
          ALLOCATE (ccw_md_mult_c        (n_mdcld_mult,  nlev))
          ALLOCATE (z_theta_md_mult_c    (n_mdcld_mult,  nlev))
          ALLOCATE (z_rho_md_mult_c      (n_mdcld_mult,  nlev))
          ALLOCATE (p_lyr_bnds_md_mult_c (n_mdcld_mult,0:nlev))

          Do i=1, n_mdcld_mult
            mdcldi_mult          (i) = dum2(i)
            iccb_md_mult_c       (i) = 0
            icct_md_mult_c       (i) = 0
            freeze_lev_md_mult_c (i) = freeze_lev (mdcldi_mult(i))
            cca_2d_md_mult_c     (i) = cca_2d_md  (mdcldi_mult(i))
          End do

          Do k=1, nlev
            Do i=1, n_mdcld_mult
              cca_md_mult_c       (i,k) = 0.0
              ccw_md_mult_c       (i,k) = ccw_md            (mdcldi_mult(i),k)
              z_theta_md_mult_c   (i,k) = z_theta           (mdcldi_mult(i),k)
              z_rho_md_mult_c     (i,k) = z_rho             (mdcldi_mult(i),k)
              p_lyr_bnds_md_mult_c(i,k) = p_layer_boundaries(mdcldi_mult(i),k)
            End Do
          End Do

          ! For multiple mid-level cloud banks, increment up model levels
          ! through compressed ccw_md_c array. Enter CAL3DCCA on locating
          ! cloud/bases of multiple clouds.

          Do k=2, nlev
            Do i=1, n_mdcld_mult
 
              !---------------------------------------------------------------
              ! Check for cloud base
              !---------------------------------------------------------------
              If ((ccw_md_mult_c(i,k) >  0.0) .AND.                           &
                  (iccb_md_mult_c(i)  == 0)   .AND.                           &
                  (icct_md_mult_c(i)  == 0)) Then
                iccb_md_mult_c(i) = k
              End if

              !---------------------------------------------------------------
              ! Check for cloud top
              !---------------------------------------------------------------
              If ((ccw_md_mult_c(i,k)   <= 0.0) .AND.                         &
                  (ccw_md_mult_c(i,k-1) >  0.0) .AND.                         &
                  (iccb_md_mult_c(i)    /= 0)   .AND.                         &
                  (icct_md_mult_c(i)    == 0)) Then
                icct_md_mult_c(i) = k-1
              End if

              !---------------------------------------------------------------
              ! Check for for anvil if both a cloud base and top found
              !---------------------------------------------------------------
              If (iccb_md_mult_c(i) /= 0 .AND.                                &
                 icct_md_mult_c(i) /= 0) Then

                ! Apply CCA to vertically continuous cloud
                If (Lcv_3d_cca) Then
! DEPENDS ON: calc_3d_cca
                  Call calc_3d_cca(1, 1, nlev, nbl, iccb_md_mult_c(i)         &
                     , icct_md_mult_c(i), p_lyr_bnds_md_mult_c(i,:)           &
                     , freeze_lev_md_mult_c(i), cca_2d_md_mult_c(i)           &
                     , cca_md_mult_c(i,:), z_theta_md_mult_c(i,:)             &
                     , z_rho_md_mult_c(i,:), l_q_interact, .false.            &
                     , l_pc2_diag_sh_pts)
                Else

                  ! Copy cca_2d_md_c to all levels from cloud base to cloud
                  ! top of located cloud

                  Do j = iccb_md_mult_c(i), icct_md_mult_c(i)
                    cca_md_mult_c(i,j) = cca_2d_md_mult_c(i)
                  End Do      ! j
                End If      ! Lcv_3d_cca

                iccb_md_mult_c(i) = 0
                icct_md_mult_c(i) = 0

              End If      ! Test on iccb_md_c(i) and icct_md_c(i)
            End Do      ! i (n_mdcld_mult)
          End Do      ! k (nlev)


          !-------------------------------------------------------------------
          ! Merge cca_md/ccw_md to cca full array and scale ccw_md_c by
          ! ccw_md_knob
          !-------------------------------------------------------------------

          Do k = 1, nlev
            Do i=1, n_mdcld_mult
              cca(mdcldi_mult(i),k) = cca(mdcldi_mult(i),k)                   &
                                    + cca_md_mult_c(i,k)
              ccw(mdcldi_mult(i),k) = ccw(mdcldi_mult(i),k)                   &
                                    + ccw_md_mult_c(i,k)*ccw_md_knob
            End Do      ! i (nlev)
          End Do      ! k (n_mdcld_mult)


          !-------------------------------------------------------------------
          ! Deallocate arrays
          !-------------------------------------------------------------------
          DEALLOCATE (mdcldi_mult)
          DEALLOCATE (iccb_md_mult_c)
          DEALLOCATE (icct_md_mult_c)
          DEALLOCATE (freeze_lev_md_mult_c)
          DEALLOCATE (cca_2d_md_mult_c)
          DEALLOCATE (cca_md_mult_c)
          DEALLOCATE (ccw_md_mult_c)
          DEALLOCATE (z_theta_md_mult_c)
          DEALLOCATE (z_rho_md_mult_c)
          DEALLOCATE (p_lyr_bnds_md_mult_c)

        End If      ! n_mdcld_mult
      End If      ! lcv_ccrad

!
! ---------------------------------------------------------------------
!write adaptive scm diagnostics
!commented out because current version causes problems with
!substepping - code still available for future development work.
! ---------------------------------------------------------------------
!#if defined(SCMA)
!      sname='rbuoy_p_out'
!      lname='buoyancy excess as used in Parcel'
!      units='K '
!      call SCMoutput(rbuoy_p_out,
!     &               sname,lname,units,t_inst,d_all,stream(3),'')
!!
!      sname='the_out'
!      lname='th_E as used in Parcel'
!      units='K '
!      call SCMoutput(the_out,
!     &               sname,lname,units,t_inst,d_all,stream(3),'')
!!
!      sname='thp_out'
!      lname='th_P as used in Parcel'
!      units='K '
!      call SCMoutput(the_out,
!     &               sname,lname,units,t_inst,d_all,stream(3),'')
!!
!      sname='qe_out'
!      lname='q_E as used in Parcel'
!      units='kg/kg'
!      call SCMoutput(qe_out,
!     &               sname,lname,units,t_inst,d_all,stream(3),'')
!!
!      sname='qp_out'
!      lname='q_P as used in Parcel'
!      units='kg/kg '
!      call SCMoutput(qp_out,
!     &               sname,lname,units,t_inst,d_all,stream(3),'')
!
!#endif
!
! ---------------------------------------------------------------------
! 4.0 Energy correction calculation - involves integral over
!     whole column therefore this should not be placed in the
!     separate calls to deep, shallow and mid. Some locations will
!     have both shallow and mid level or deep and mid level convection
!
!     UM documentation paper 27 - section 12.
! ---------------------------------------------------------------------
! First work out which points convection has occurred at .

      nconv_all=0
      Do i = 1,npnts
        If (cumulus_bl(i).or.l_mid_all(i)) then
          nconv_all = nconv_all + 1
          index1(nconv_all) = i
        End If
      End Do

      If (nconv_all >  0) then

! DEPENDS ON: cor_engy
        call Cor_engy(np_field,npnts,nconv_all,nlev,dthbydt,dqbydt,snow &
     &               ,exner_layer_centres,p_layer_boundaries,index1)

!-----------------------------------------------------------------------
! 5.0  Total water conservation  - also works on whole column
!-----------------------------------------------------------------------

! only check columns where convection has occurred.


        Do j = 1,nconv_all
          i=index1(j)
          qMinInColumn(j) = q(i,nlev)
        End Do
        Do k = 1,nlev-1
          Do j = 1,nconv_all
            i=index1(j)
            If (q(i,k)  <   qMinInColumn(j)) then
              qMinInColumn(j) = q(i,k)
            End If
          End Do
        End Do

!
! Ensure Q does not go below global allowed minimum (QMIN)
!
        Do j = 1,nconv_all
          qMinInColumn(j)=MAX(QMIN,qMinInColumn(j))
        End Do

!
! Apply an artificial upwards flux from k-1 level to ensure Q
! remians above minimum value in the column.
!

        Do k = nlev,2,-1
! NEC SX6 compiler directive so that next loop vectorized
!CDIR NODEP
          Do j = 1,nconv_all
            i=index1(j)
            temp1(j)=q(i,k) + dqbydt(i,k) * timestep

            If (temp1(j)  <   qMinInColumn(j)) then

              dqbydt(i,k-1) = dqbydt(i,k-1) -                           &
     &              ((qMinInColumn(j) - q(i,k)) / timestep-dqbydt(i,k)) &
     &               * ( - (p_layer_boundaries(i,k) -                   &
     &                      p_layer_boundaries(i,k-1)))                 &
     &               / ( - (p_layer_boundaries(i,k-1) -                 &
     &                      p_layer_boundaries(i,k-2)))

              dqbydt(i,k) = (qMinInColumn(j) - q(i,k)) / timestep
            End If
          End Do ! nconv_all loop
        End Do  ! nlev

!-----------------------------------------------------------------------
! 6.0  Subroutine MIX_INC mixes convective increments in the boundary
!      layer (essentially distributes incr. at ntml over layers 1 to
!      ntml e.g. incr(1) = incr(2) = incr(ntml)/ntml)
!      Works on boundary layer - columns integrals involved.
!-----------------------------------------------------------------------
        IF (bl_cnv_mix == 0) then
! DEPENDS ON: mix_inc
        Call MIX_INC(np_field,npnts,nconv_all,nlev,nbl,ntml,            &
                     dthbydt,dqbydt,dubydt,dvbydt,L_tracer,ntra,        &
                     dtrabydt,p_layer_boundaries,                       &
                     p_layer_centres,index1)


        END If
      End If    ! test on nconv_all
!-----------------------------------------------------------------------
! 7.0  Calculate convective cloud amount on model levels if
!      Lcv_3d_cca = .true.
!-----------------------------------------------------------------------
! Initialise output array

      If (.NOT. lcv_ccrad) Then
        Do k=1, n_cca_lev
          Do i=1, npnts
            cca(i,k) = 0.0
          End Do
        End Do
      End If      ! lcv_ccrad

!
! HadGEM1 3-d cloud is done in control therefore don't call here
! This is messy but unfortunately has to be done to get the same results
! as the original modset.
!
      If (.NOT. l_convcld_hadgem1) Then
        If (.NOT. lcv_ccrad) Then
          If (lcv_3d_cca) Then
            ! Resetting of tower cloud to zero in PC2 is performed within
            ! calc_3d_cca.

! DEPENDS ON: calc_3d_cca
        Call CALC_3D_CCA(np_field,npnts,nlev,nbl,iccb,icct,           &
                         p_layer_boundaries,freeze_lev,               &
                         cca_2d,cca,z_theta,z_rho,l_q_interact,       &
                         .true., l_pc2_diag_sh_pts)

          Else

            Do i = 1,npnts
              If (l_q_interact) Then
                ! Need to reset 2D convective cloud to zero for PC2
                cca(i,1) = 0.0
              Else
                ! Not PC2. Copy the 2D conv cloud amount to the CCA field.
                cca(i,1) = cca_2d(i)
              End If
            End Do      ! i (npnts)

          End If      ! lcv_3d_cca
        End If      ! lcv_ccrad
      End If      ! l_convcld_hadgem1

!-----------------------------------------------------------------------
! 8.0  Update tracer field
!      More efficient to do here rather than in subroutines.
!      This seems to be an expensive bit of code on NEC.
!      Changed to operate on just convective points.
!-----------------------------------------------------------------------

      If (L_tracer.and.(nconv_all >  0)) then

!
! Adjust timestep to prevent any negative values invading the tracer
! fields (adjusted timestep is a function of geograhical  location and
! tracer type.
!
        Do ktra=1,ntra
! initialise step to timestep
          Do i = 1,nconv_all
            limited_step(i) =timestep
          End Do

          Do k=1,nlev
            Do j = 1,nconv_all
              i=index1(j)
! negative increments  may have a problem
              If (dtrabydt(i,k,ktra)  <   0.0 ) then

                step_test(j) = (0.9999 * abs(tracer(i,k,ktra))) /       &
     &                    (abs(dtrabydt(i,k,ktra)) + SAFETY_MARGIN)

                If (step_test(j)   <   limited_step(j) ) then
! then increment is bigger than tracer and timestep needs to be reduced
                  limited_step (j) = step_test(j)
                End If

             End If
            End Do
          End Do
!
! Update tracer field using limited_step.
!
          Do k = 1,nlev
! NEC compiler directive
!CDIR NODEP
            Do j = 1,nconv_all
             i=index1(j)
            tracer(i,k,ktra) = tracer(i,k,ktra) +                       &
     &                         dtrabydt(i,k,ktra) * limited_step(j)
            End Do
          End Do

        End Do  ! ktra loop
      End If   ! L_tracer

!-----------------------------------------------------------------------

      Return
      END SUBROUTINE GLUE_CONV
