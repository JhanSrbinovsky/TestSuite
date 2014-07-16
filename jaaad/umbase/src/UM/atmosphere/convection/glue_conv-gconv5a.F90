#if defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
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
     &,                 N_CUMULUS,UW0,VW0,w_max,zlcl,ZLCL_UV,ZTOP_UV    &
     &,                 LCCLWP,CAPE_OUT                                 &
     &,                 N_DP,n_cg, N_SH                                 &
     &,                 r_rho,r_theta,rho, rho_theta                    &
     &,                 EXNER_LAYER_BOUNDARIES                          &
     &,                 EXNER_LAYER_CENTRES                             &
     &,                 P_LAYER_BOUNDARIES                              &
     &,                 P_LAYER_CENTRES                                 &
     &,                 z_theta, z_rho                                  &
     &,                 TIMESTEP,T1_SD,Q1_SD                            &
     &,                 NTML,NTPAR,NBDSC,NTDSC,L_SHALLOW_BL             &
     &,                 l_pc2_diag_sh_pts,l_congestus                   &
     &,                 CUMULUS_BL,WSTAR,WTHVS,DELTHVU_BL,ql_ad,ftl,fqt &
     &,                 L_TRACER, NTRA, TRLEV, N_CCA_LEV                &
     &,                 l_mixing_ratio, L_CALC_DXEK, L_Q_INTERACT       &
     &,                 UP_FLUX_HALF,FLG_UP_FLX_HALF                    &
     &,                 UP_FLUX,FLG_UP_FLX,DWN_FLUX,FLG_DWN_FLX         &
     &,                 ENTRAIN_UP,FLG_ENTR_UP,DETRAIN_UP               &
     &,                 FLG_DETR_UP,ENTRAIN_DWN,FLG_ENTR_DWN            &
     &,                 DETRAIN_DWN,FLG_DETR_DWN,UW_DEEP,FLG_UW_DEEP    &
     &,                 VW_DEEP,FLG_VW_DEEP,UW_SHALL,FLG_UW_SHALL       &
     &,                 VW_SHALL,FLG_VW_SHALL                           &
     &,          flg_wqt_flux_sh,flg_wthetal_flux_sh                    &
     &,          flg_wthetav_flux_sh,flg_wql_flux_sh                    &
     &,          flg_mf_deep,flg_mf_congest,flg_mf_shall,flg_mf_midlev  &
     &,          flg_dt_deep,flg_dt_congest,flg_dt_shall,flg_dt_midlev  &
     &,          flg_dq_deep,flg_dq_congest,flg_dq_shall,flg_dq_midlev  &
     &,          flg_du_deep,flg_du_congest,flg_du_shall,flg_du_midlev  &
     &,          flg_dv_deep,flg_dv_congest,flg_dv_shall,flg_dv_midlev  &
     &,         wqt_flux_sh,wthetal_flux_sh,wthetav_flux_sh,wql_flux_sh &
     &,          mf_deep,mf_congest,mf_shall,mf_midlev                  &
     &,          dt_deep,dt_congest,dt_shall,dt_midlev                  &
     &,          dq_deep,dq_congest,dq_shall,dq_midlev                  &
     &,          du_deep,du_congest,du_shall,du_midlev                  &
     &,          dv_deep,dv_congest,dv_shall,dv_midlev                  &
     &                  )

! Purpose:
!  Gather-scatter routine for deep and shallow convection points.
!  Interface to deep, shallow and mid level convection
!
!   Called by Load_bal_conv
!
! Current owners of code: Convection code owner
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
      Use cv_cntl_mod, Only:                                            &
          lcv_3d_cca, lcv_ccrad

      Use cv_run_mod,  Only:                                            &
          l_mom, adapt, iconv_shallow, iconv_congestus, iconv_mid,      &
          iconv_deep,   termconv,                                       &
          cca2d_sh_opt, cca2d_md_opt,  cca2d_dp_opt,                    &
          cca_sh_knob,  cca_md_knob,   cca_dp_knob,                     &
          ccw_sh_knob,  ccw_md_knob,   ccw_dp_knob 

      IMPLICIT NONE

! Model Constants required

#include "c_0_dg_c.h"
#include "c_g.h"

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

      real, intent(in) :: &
        ql_ad(npnts)      & ! adiabatic liquid water content at inversion (kg/kg)
      , ftl(npnts)        & ! Surface sensible heat flux from BL
                            ! (W/m2) i.e. cp*rho*w'tl'
      , fqt(npnts)          ! Total water flux from surface (kg/m2/s)
                            ! i.e. rho*w'qT'
   

      real, intent(in) :: r_rho(np_field,nlev)    ! radius rho levels(m)
      real, intent(in) :: r_theta(np_field,0:nlev)  ! theta levels (m)
      real, intent(in) :: rho(np_field,nlev)     ! density kg/m3
      real, intent(in) :: rho_theta(np_field,nlev) ! th lev kg/m3

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
                                      ! height (m) not a model level.

      real, intent(in) :: zlcl_uv(npnts) ! Lifting condensation level
                                      ! defined for the uv grid (m)

      real, intent(in) :: ztop_uv(npnts) ! Top of cloud layer
                                      ! defined for the uv
                                      ! grid (m)

      logical, intent(in) :: L_mixing_ratio
                               ! true - input moisture mixing ratio
                               ! false - input moisture specific humidity

      logical, intent(in) :: L_tracer ! Switch for inclusion of tracers

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
                                    ! mass flux on rho (dummy variable)

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
! flags to control model level diagnostics of mass flux, and increments
! from the individual convection schemes.
!
      logical, intent(in) ::                                            &
     &           flg_mf_deep,flg_mf_congest,flg_mf_shall,flg_mf_midlev  &
     &,          flg_dt_deep,flg_dt_congest,flg_dt_shall,flg_dt_midlev  &
     &,          flg_dq_deep,flg_dq_congest,flg_dq_shall,flg_dq_midlev  &
     &,          flg_du_deep,flg_du_congest,flg_du_shall,flg_du_midlev  &
     &,          flg_dv_deep,flg_dv_congest,flg_dv_shall,flg_dv_midlev

!
! Flags for turbulence diagnostics
!
      logical, intent(in) ::                                            &
     &  flg_wqt_flux_sh                                                 &
                             ! STASH flag for wqt flux
     &, flg_wthetal_flux_sh                                             &
                             ! STASH flag for wthetal flux
     &, flg_wthetav_flux_sh                                             &
                             ! STASH flag for wthetav flux
     &, flg_wql_flux_sh      ! STASH flag for wql flux

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
      real, intent(out) :: precip_cong(npnts) ! congest precip (kg/m2/s)

      real, intent(out) ::                                              &
     &  wstar_dn(npnts)                                                 &
                         ! subcloud layer convective velocity scale(m/s)
     &, wstar_up(npnts)                                                 &
                         ! cumulus layer convective velocity scale (m/s)
     &, mb1(npnts)                                                      &
                         ! cloud base mass flux from wstar_dn (m/s)
     &, mb2(npnts)       ! cloud base mass flux for cloud layer (m/s)

      integer, intent(out) :: kterm_congest(npnts) ! termination level
!                                                      for congestus

      real, intent(out) :: lcclwp(npnts) ! Condensed water path for
                                      ! lowest conv. cld. level (kg/m2)

!
! Meaning of this diagnostic depends on scheme
!  Plume model - updraught mass flux  (Pa/s)
!  turbulence model - mass flux (not exactly updraught) (m/s)

      real, intent(out) :: up_flux(np_field,nlev) ! mass flux

      real, intent(out) :: up_flux_half(np_field,nlev) !mass flux on rho
                           !dummy variable not used in turbulence

!
! Diagnostics with no meaning for turbulence based schemes
!
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
!
! Diagnostics relating to momentum fluxes
!
      real, intent(out) ::                                              &
     &  uw_deep(np_field,nlev)                                          &
                               ! X-comp. of stress from deep convection
                               !(kg/m/s2)
     &, vw_deep(np_field,nlev)                                          &
                               ! Y-comp. of stress from deep convection
                               !(kg/m/s2)
     &, uw_shall(np_field,nlev)                                         &
                                ! X-comp. of stress from shallow
                                ! convection (kg/m/s2)
     &, vw_shall(np_field,nlev) ! Y-comp. of stress from shallow
                                ! convection (kg/m/s2)

!      real, intent(out) ::  &  ! plan to make output diagnostics in future
      real ::               &
      uw_mid(np_field,nlev) & ! U comp of stress from mid convection (kg/m/s2)
     ,vw_mid(np_field,nlev)   ! V comp of stress from mid convection (kg/m/s2)

      real, intent(out) :: cape_out(npnts) ! Saved convective available
                                      ! potential energy for diagnostic
                                      ! output (Jkg-1)

! Fluxes from turbulence based convection schemes

      real, intent(out) ::                                              &
     &  wqt_flux_sh(np_field,nlev)                                      &
                                        ! w'qt' flux (m/s kg/kg)
     &, wthetal_flux_sh(np_field,nlev)                                  &
                                        ! w'thetal' flux  (m/s K)
     &, wthetav_flux_sh(np_field,nlev)                                  &
                                        ! w'thetav' flux  (m/s K)
     &, wql_flux_sh(np_field,nlev)                                      &
                                        ! w'ql' flux  (m/s kg/kg)
     &, mf_deep(np_field,nlev)                                          &
                                     ! mass flux deep
     &, mf_congest(np_field,nlev)                                       &
                                     ! mass flux congestus
     &, mf_shall(np_field,nlev)                                         &
                                     ! mass flux shallow
     &, mf_midlev(np_field,nlev)                                        &
                                     ! mass flux mid-lev
     &, dt_deep(np_field,nlev)                                          &
                                     ! dt increment deep   (K/s)
     &, dt_congest(np_field,nlev)                                       &
                                     ! dt increment congestus (K/s)
     &, dt_shall(np_field,nlev)                                         &
                                     ! dt increment shallow (K/s)
     &, dt_midlev(np_field,nlev)                                        &
                                     ! dt increment mid-level (K/s)
     &, dq_deep(np_field,nlev)                                          &
                                     ! dq increment deep (kg/kg/s)
     &, dq_congest(np_field,nlev)                                       &
                                     ! dq increment congestus (kg/kg/s)
     &, dq_shall(np_field,nlev)                                         &
                                     ! dq increment shallow (kg/kg/s)
     &, dq_midlev(np_field,nlev)                                        &
                                     ! dq increment mid-level (kg/kg/s)
     &, du_deep(np_field,nlev+1)                                        &
                                       ! du increment deep (m/s)
     &, du_congest(np_field,nlev+1)                                     &
                                       ! du increment congestus (m/s)
     &, du_shall(np_field,nlev+1)                                       &
                                       ! du increment shallow (m/s)
     &, du_midlev(np_field,nlev+1)                                      &
                                       ! du increment mid-level (m/s)
     &, dv_deep(np_field,nlev+1)                                        &
                                       ! dv increment deep (m/s)
     &, dv_congest(np_field,nlev+1)                                     &
                                       ! dv increment congestus (m/s)
     &, dv_shall(np_field,nlev+1)                                       &
                                       ! dv increment shallow (m/s)
     &, dv_midlev(np_field,nlev+1)     ! dv increment mid-level (m/s)


!-----------------------------------------------------------------------
! Redundant arguments
!-----------------------------------------------------------------------

      integer, intent(in) :: n_cumulus
      real, intent(in)    :: ntdsc
      real, intent(in)    :: nbdsc


! local required in move mixing to glue

      real :: dtrabydt(npnts,nlev,ntra) ! Increment to tracer due to
                                      ! convection (kg/kg/s)

      Character*(*), Parameter ::  RoutineName = 'glue_conv'

!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to DEEP convection scheme.
! Arrays are identified by underscore DP (_dp) and are of length
! n_dp where n_dp is the number of points diagnosed as deep in the
! boundary layer diagnosis routine. For full desriptions of variables
! see above.
!-----------------------------------------------------------------------

      integer :: dpi(n_dp)        ! index for deep points in full grid

      integer ::                                                        &
     &  ntml_dp(n_dp)                                                   &
     &, ntpar_dp(n_dp)

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
     , r_rho_dp(n_dp,nlev)                                              &
     , r_theta_dp(n_dp,0:nlev)     

      logical :: bland_dp(n_dp)

      real ::                                                           &
     &  u_dp(n_dp,nlev)                                                 &
     &, v_dp(n_dp,nlev)                                                 &
     &, th_dp(n_dp,nlev)                                                &
     &, q_dp(n_dp,nlev)                                                 &
     &, qse_dp(n_dp,nlev)

      real :: tracer_dp(n_dp,trlev,ntra)
! increments
      real ::                                                           &
     &  dthbydt_dp(n_dp,nlev)                                           &
     &, dqbydt_dp(n_dp,nlev)                                            &
     &, dubydt_dp(n_dp,nlev+1)                                          &
     &, dvbydt_dp(n_dp,nlev+1)

! output variables
      real ::                                                           &
     &  rain_dp(n_dp)                                                   &
     &, snow_dp(n_dp)                                                   &
     &, rain_3d_dp(n_dp,nlev)                                           &
     &, snow_3d_dp(n_dp,nlev)                                           &
     &, tcw_dp(n_dp)                                                    &
     &, cclwp_dp(n_dp)                                                  &
     &, lcca_dp(n_dp)                                                   &
     &, lcclwp_dp(n_dp)                                                 &
     &, cape_out_dp(n_dp)

      integer ::                                                        &
     &  iccb_dp(n_dp)                                                   &
     &, icct_dp(n_dp)                                                   &
     &, lcbase_dp(n_dp)                                                 &
     &, lctop_dp(n_dp)                                                  &
     &, freeze_lev_dp(n_dp)

      real ::                                                           &
     &  ccw_dp(n_dp,nlev)                                               &
     &, up_flux_dp(n_dp,nlev)                                           &
     &, dwn_flux_dp(n_dp,nlev)                                          &
     &, entrain_up_dp(n_dp,nlev)                                        &
     &, detrain_up_dp(n_dp,nlev)                                        &
     &, entrain_dwn_dp(n_dp,nlev)                                       &
     &, detrain_dwn_dp(n_dp,nlev)                                       &
     &, uw_deep_dp(n_dp,nlev)                                           &
     &, vw_deep_dp(n_dp,nlev)                                           &
     &, W_MAX_dp(n_dp)

!PC2
      real ::                                                           &
     &  qcl_dp(n_dp,nlev)                                               &
     &, qcf_dp(n_dp,nlev)                                               &
     &, cf_liquid_dp(n_dp,nlev)                                         &
     &, cf_frozen_dp(n_dp,nlev)                                         &
     &, bulk_cf_dp(n_dp,nlev)                                           &
     &, dqclbydt_dp(n_dp,nlev)                                          &
     &, dqcfbydt_dp(n_dp,nlev)                                          &
     &, dcflbydt_dp(n_dp,nlev)                                          &
     &, dcffbydt_dp(n_dp,nlev)                                          &
     &, dbcfbydt_dp(n_dp,nlev)

      real :: dtrabydt_dp(n_dp,nlev,ntra) ! Increment to tracer due to
                                          ! convection (kg/kg/s)
      integer :: kterm_dp(n_dp)    ! required by mid level scheme

      real ::                 &
          cca_2d_dp(n_dp)     & ! required by mid level scheme
      ,   cca_dp(n_dp,nlev)     ! 3d CCA     

      integer :: dp_new_termc !flag for simplified termination of conv
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

      integer :: ishall_precip  ! flag for precip in shallow turbulence 
                                ! scheme 
                                !  0 - no precip
                                !  1 - precip

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
! LOCAL compressed arrays to be passed to congestus convection scheme.
! Arrays are identified by underscore cg (_cg) and are of length
! n_cg where n_cg is the number of points diagnosed as congestus in the
! boundary layer diagnosis routine.For full desriptions of variables
! see above.
!-----------------------------------------------------------------------

      integer :: cgi(n_cg)     ! index for congestus points in full grid

      integer ::                                                        &
     &  ntml_cg(n_cg)                                                   &
     &, ntpar_cg(n_cg)

      real ::                                                           &
     &  pstar_cg(n_cg)                                                  &
     &, recip_pstar_cg(n_cg)                                            &
     &, delthvu_cg(n_cg)                                                &
     &, t1_sd_cg(n_cg)                                                  &
     &, q1_sd_cg(n_cg)                                                  &
     &, uw0_cg(n_cg)                                                    &
     &, vw0_cg(n_cg)                                                    &
     &, wstar_cg(n_cg)                                                  &
     &, wthvs_cg(n_cg)                                                  &
     &, zlcl_uv_cg(n_cg)                                                &
     &, ztop_uv_cg(n_cg)

      real ::                                                           &
     &  p_layer_centres_cg(n_cg,0:nlev)                                 &
     &, p_layer_boundaries_cg(n_cg,0:nlev)                              &
     &, exner_layer_centres_cg(n_cg,0:nlev)                             &
     &, exner_layer_boundaries_cg(n_cg,0:nlev)                          &
     &, z_theta_cg(n_cg,nlev)                                           &
     &, z_rho_cg(n_cg,nlev)                                             &
     &, u_cg(n_cg,nlev)                                                 &
     &, v_cg(n_cg,nlev)                                                 &
     &, th_cg(n_cg,nlev)                                                &
     &, q_cg(n_cg,nlev)                                                 &
     &, ccw_cg(n_cg,nlev)                                               &
     &, qse_cg(n_cg,nlev)                                               &
     &, rho_cg(n_cg,nlev)                                               &
     &, r_rho_cg(n_cg,nlev),r_theta_cg(n_cg,0:nlev)                     &
     &, r2rho_th_cg(n_cg,nlev)                                          &
                                    ! radius**2 density theta lev (kg/m)
     &, r2rho_cg(n_cg,nlev)                                             &
                                    ! radius**2 density rho lev (kg/m)
     &, dr_across_th_cg(n_cg,nlev)                                      &
                                    ! thickness of theta levels (m)
     &, dr_across_rh_cg(n_cg,nlev)                                      &
                                    ! thickness of rho levels (m)
     &, rho_th_cg(n_cg,nlev)        ! rho on theta levels (kg/m3)

      logical :: bland_cg(n_cg)

      integer :: cg_new_termc !flag for simplified termination of conv
                              !0 = off, 1 = on

      integer :: cg_on        !flag for adaptive applied to deep conv
                              !0 = off, 1 = on

      integer :: mdet_cg_on   !flag for adaptive mixing detrainment
                              !applied to deep conv
                              !0 = off, 1 = on

      integer :: cg_ent_on    !flag for adaptive mod for entrainment
                              !0 = off, 1 = on

      integer :: cg_sdet_on   !flag for smoothed forced detrainment
                              !0 = off, 1 = on                        
                              
! tracers
      real :: tracer_cg(n_cg,trlev,ntra)                                &
     &, dtrabydt_cg(n_cg,nlev,ntra)

      real ::                                                           &
     &  dthbydt_cg(n_cg,nlev)                                           &
     &, dqbydt_cg(n_cg,nlev)                                            &
     &, dubydt_cg(n_cg,nlev+1)                                          &
     &, dvbydt_cg(n_cg,nlev+1)

      real ::                                                           &
     &  rain_cg(n_cg)                                                   &
     &, snow_cg(n_cg)                                                   &
     &, tcw_cg(n_cg)                                                    &
     &, cclwp_cg(n_cg)                                                  &
     &, lcca_cg(n_cg)                                                   &
     &, lcclwp_cg(n_cg)                                                 &
     &, cape_out_cg(n_cg)
      real ::             &
       cca_2d_cg(n_cg)    & ! required by mid level scheme
      ,cca_cg(n_cg,nlev)    ! 3d CCA     

      integer ::                                                        &
     &  iccb_cg(n_cg)                                                   &
     &, icct_cg(n_cg)                                                   &
     &, lcbase_cg(n_cg)                                                 &
     &, lctop_cg(n_cg)                                                  &
     &, freeze_lev_cg(n_cg)                                             &
     &, kterm_cg(n_cg)    ! required by mid level scheme

! diagnostics
      real ::                                                           &
     &  up_flux_cg(n_cg,nlev)                                           &
     &, dwn_flux_cg(n_cg,nlev)                                          &
     &, entrain_up_cg(n_cg,nlev)                                        &
     &, detrain_up_cg(n_cg,nlev)                                        &
     &, entrain_dwn_cg(n_cg,nlev)                                       &
     &, detrain_dwn_cg(n_cg,nlev)                                       &
     &, uw_shall_cg(n_cg,nlev)                                          &
     &, vw_shall_cg(n_cg,nlev)

!PC2
      real ::                                                           &
     &  qcl_cg(n_cg,nlev)                                               &
     &, qcf_cg(n_cg,nlev)                                               &
     &, cf_liquid_cg(n_cg,nlev)                                         &
     &, cf_frozen_cg(n_cg,nlev)                                         &
     &, bulk_cf_cg(n_cg,nlev)                                           &
     &, dqclbydt_cg(n_cg,nlev)                                          &
     &, dqcfbydt_cg(n_cg,nlev)                                          &
     &, dcflbydt_cg(n_cg,nlev)                                          &
     &, dcffbydt_cg(n_cg,nlev)                                          &
     &, dbcfbydt_cg(n_cg,nlev)

! may be required for turbulence version
      real ::                                                           &
     &  wstar_dn_cg(n_cg)                                               &
     &, wstar_up_cg(n_cg)                                               &
     &, mb1_cg(n_cg)                                                    &
     &, mb2_cg(n_cg)                                                    &
     &, wthetal_cg(n_cg,nlev)                                           &
     &, wthetav_cg(n_cg,nlev)                                           &
     &, wql_cg(n_cg,nlev)                                               &
     &, wqt_cg(n_cg,nlev)

!  compressed SCM diagnostics for adaptive modset

      real :: rbuoy_p_out_cg(n_cg,nlev)

      real :: the_out_cg(n_cg,nlev)

      real :: thp_out_cg(n_cg,nlev)

      real :: qe_out_cg(n_cg,nlev)

      real :: qp_out_cg(n_cg,nlev)


!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to SHALLOW convection scheme.
! Arrays are identified by underscore SH (_sh) and are of length
! n_sh where n_sh is the number of points diagnosed as shallow in the
! boundary layer diagnosis routine.For full desriptions of variables
! see above.
!-----------------------------------------------------------------------

      integer :: shi(n_sh)      ! index for shallow points in full grid

! compressed inputs

      integer :: ntml_sh(n_sh)                                          &
     &, ntpar_sh(n_sh)

      real ::                                                           &
     &  pstar_sh(n_sh)                                                  &
     &, recip_pstar_sh(n_sh)                                            &
     &, delthvu_sh(n_sh)                                                &
     &, ql_ad_sh(n_sh)                                                  &
     &, t1_sd_sh(n_sh)                                                  &
     &, q1_sd_sh(n_sh)                                                  &
     &, uw0_sh(n_sh)                                                    &
     &, vw0_sh(n_sh)                                                    &
     &, wstar_sh(n_sh)                                                  &
     &, wthvs_sh(n_sh)                                                  &
     &, zlcl_uv_sh(n_sh)                                                &
     &, zlcl_sh(n_sh)                                                   &
     &, ztop_uv_sh(n_sh)

      real ::                                                           &
     &  p_layer_centres_sh(n_sh,0:nlev)                                 &
     &, p_layer_boundaries_sh(n_sh,0:nlev)                              &
     &, exner_layer_centres_sh(n_sh,0:nlev)                             &
     &, exner_layer_boundaries_sh(n_sh,0:nlev)                          &
     &, z_theta_sh(n_sh,nlev), z_rho_sh(n_sh,nlev)                      &
     &, rho_sh(n_sh,nlev)                                               &
     &, r_rho_sh(n_sh,nlev),r_theta_sh(n_sh,0:nlev)                     &
     &, r2rho_th_sh(n_sh,nlev)                                          &
                                    ! radius**2 density theta lev (kg/m)
     &, r2rho_sh(n_sh,nlev)                                             &
                                    ! radius**2 density rho lev (kg/m)
     &, dr_across_th_sh(n_sh,nlev)                                      &
                                    ! thickness of theta levels (m)
     &, dr_across_rh_sh(n_sh,nlev)                                      &
                                    ! thickness of rho levels (m)
     &, rho_th_sh(n_sh,nlev)        ! rho on theta levels (kg/m3)

      logical :: bland_sh(n_sh)

      real ::                                                           &
     &  u_sh(n_sh,nlev)                                                 &
     &, v_sh(n_sh,nlev)                                                 &
     &, th_sh(n_sh,nlev)                                                &
     &, q_sh(n_sh,nlev)                                                 &
     &, qse_sh(n_sh,nlev)

! output
      real ::                                                           &
     &  dthbydt_sh(n_sh,nlev)                                           &
     &, dqbydt_sh(n_sh,nlev)                                            &
     &, dubydt_sh(n_sh,nlev+1)                                          &
     &, dvbydt_sh(n_sh,nlev+1)

      real ::                 &
        rain_sh(n_sh)         &
      , snow_sh(n_sh)         &
      , rain_3d_sh(n_sh,nlev) &
      , snow_3d_sh(n_sh,nlev) &
      , tcw_sh(n_sh)          &
      , cclwp_sh(n_sh)        &
      , lcca_sh(n_sh)         &
      , lcclwp_sh(n_sh)       &
      , cape_out_sh(n_sh)     &     
      , cca_2d_sh(n_sh)       &   ! required by mid level scheme
      , cca_sh(n_sh,nlev)         ! 3d CCA     

      real ::                                                           &
     &  wqt_sh(n_sh,nlev)                                               &
     &, wthetal_sh(n_sh,nlev)                                           &
     &, wthetav_sh(n_sh,nlev)                                           &
     &, wql_sh(n_sh,nlev)                                               &
     &, wstar_dn_sh(n_sh)                                               &
     &, wstar_up_sh(n_sh)                                               &
     &, mb1_sh(n_sh)                                                    &
     &, mb2_sh(n_sh)

      integer :: iccb_sh(n_sh)                                          &
     &, icct_sh(n_sh)                                                   &
     &, lcbase_sh(n_sh)                                                 &
     &, lctop_sh(n_sh)                                                  &
     &, freeze_lev_sh(n_sh)

      real ::                                                           &
     &  ccw_sh(n_sh,nlev)                                               &
     &, up_flux_sh(n_sh,nlev)                                           &
     &, dwn_flux_sh(n_sh,nlev)                                          &
     &, entrain_up_sh(n_sh,nlev)                                        &
     &, detrain_up_sh(n_sh,nlev)                                        &
     &, entrain_dwn_sh(n_sh,nlev)                                       &
     &, detrain_dwn_sh(n_sh,nlev)                                       &
     &, uw_shall_sh(n_sh,nlev)                                          &
     &, vw_shall_sh(n_sh,nlev)

! PC2  input
      real ::                                                           &
     &  qcl_sh(n_sh,nlev)                                               &
     &, qcf_sh(n_sh,nlev)                                               &
     &, cf_liquid_sh(n_sh,nlev)                                         &
     &, cf_frozen_sh(n_sh,nlev)                                         &
     &, bulk_cf_sh(n_sh,nlev)
! PC2 output
      real ::                                                           &
     &  dqclbydt_sh(n_sh,nlev)                                          &
     &, dqcfbydt_sh(n_sh,nlev)                                          &
     &, dcflbydt_sh(n_sh,nlev)                                          &
     &, dcffbydt_sh(n_sh,nlev)                                          &
     &, dbcfbydt_sh(n_sh,nlev)

! Tracers in and out
      real :: tracer_sh(n_sh,trlev,ntra)
      real :: dtrabydt_sh(n_sh,nlev,ntra) ! Increment to tracer due to
                                     ! convection (kg/kg/s)

      integer :: sh_new_termc !flag for simplified termination of
                              ! convection

      integer :: sh_on        !flag for adaptive applied to shallow conv

      integer :: mdet_sh_on   !flag for adaptive mixing detrainment
                              !applied to shallow conv

      integer :: sh_ent_on    !flag for adaptive applied to entrainment

      integer :: sh_sdet_on   !flag for smoothed forced detrainment
                              !0 = off, 1 = on                        
      
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


      integer :: midtrig(npnts)   ! Level at which mid level convection
                                  ! may start

      real ::               &
       cca_2d_md(npnts)     & ! required by mid level scheme
      ,cca_md(npnts,nlev)     ! 3d CCA     

      real ::                                                           &
     &  p_layer_centres_md(npnts,0:nlev)                                &
     &, p_layer_boundaries_md(npnts,0:nlev)                             &
     &, exner_layer_centres_md(npnts,0:nlev)                            &
     &, exner_layer_boundaries_md(npnts,0:nlev)                         &
     &, z_theta_md(npnts,nlev), z_rho_md(npnts,nlev)                    &
     &, r_theta_md(npnts,0:nlev), r_rho_md(npnts,nlev)                  &
     &, u_md(npnts,nlev)                                                &
     &, v_md(npnts,nlev)                                                &
     &, th_md(npnts,nlev)                                               &
     &, q_md(npnts,nlev)

! output variables
      real ::                                                           &
     &  dthbydt_md(npnts,nlev)                                          &
     &, dqbydt_md(npnts,nlev)                                           &
     &, dubydt_md(npnts,nlev+1)                                         &
     &, dvbydt_md(npnts,nlev+1)

      logical :: l_mid_all(npnts)     ! true if mid level convection

      real :: rain_md(npnts)                                            &
     &, snow_md(npnts)                                                  &
     &, rain_3d_md(npnts,nlev)                                          &
     &, snow_3d_md(npnts,nlev)                                          &
     &, tcw_md(npnts)                                                   &
     &, cclwp_md(npnts)                                                 &
     &, ccw_md(npnts,nlev)                                              &
     &, lcca_md(npnts)                                                  &
     &, lcclwp_md(npnts)                                                &
     &, cape_out_md(npnts)

      integer ::                                                        &
     &  iccb_md(npnts)                                                  &
     &, icct_md(npnts)                                                  &
     &, lcbase_md(npnts)                                                &
     &, lctop_md(npnts)                                                 &
     &, freeze_lev_md(npnts)

      real ::                                                           &
     &  up_flux_md(npnts,nlev)                                          &
     &, dwn_flux_md(npnts,nlev)                                         &
     &, entrain_up_md(npnts,nlev)                                       &
     &, detrain_up_md(npnts,nlev)                                       &
     &, entrain_dwn_md(npnts,nlev)                                      &
     &, detrain_dwn_md(npnts,nlev)

! PC2
      real ::                                                           &
     &  qcl_md(npnts,nlev)                                              &
     &, qcf_md(npnts,nlev)                                              &
     &, cf_liquid_md(npnts,nlev)                                        &
     &, cf_frozen_md(npnts,nlev)                                        &
     &, bulk_cf_md(npnts,nlev)                                          &
     &, qse_md(npnts,nlev)

      real ::                                                           &
     &  dqclbydt_md(npnts,nlev)                                         &
     &, dqcfbydt_md(npnts,nlev)                                         &
     &, dcflbydt_md(npnts,nlev)                                         &
     &, dcffbydt_md(npnts,nlev)                                         &
     &, dbcfbydt_md(npnts,nlev)

! tracers
      real :: tracer_md(npnts,trlev,ntra)
      real :: dtrabydt_md(npnts,nlev,ntra) ! Increment to tracer due to
                                           ! convection (kg/kg/s)
      real :: rbuoy_p_out_md(npnts,nlev)

      real :: the_out_md(npnts,nlev)

      real :: thp_out_md(npnts,nlev)

      real :: qe_out_md(npnts,nlev)

      real :: qp_out_md(npnts,nlev)

      integer :: md_new_termc !flag for simplified termination of
                              !convection 0 = off, 1 = on

      integer :: md_on        !flag for adaptive applied to mid conv
                              !0 = off, 1 = on

      integer :: mdet_md_on   !flag for adaptive mixing detrainment
                              !0 = off, 1 = on
                              !applied to mid conv

      integer :: md_ent_on    !Flag for adaptive mod for entrainment
      
      integer :: md_sdet_on   !flag for smoothed forced detrainment
                              !0 = off, 1 = on                        

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
!  Let m indicate mixing ratio and q specific then
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
!   or
!
!  dqv = (dmv + mt*dmv - mv*dmt)/[(1+mt)(1+mt+dmt)]
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

! Saturation mixing ratio calulation and 1/pstar

      real ::                                                           &
     &  recip_pstar(npnts)                                              &
                                ! Reciprocal of pstar array
     &, qse(npnts,nlev)                                                 &
                                ! Saturation mixing ratio of cloud
                                ! cloud environment (kg/kg)
     &, pt(npnts)                                                       &
                                ! Temporary store for P in calc. of sat.
                                ! mixing ratio. (Pa)
     &, tt(npnts)                                                       &
                                ! Temporary store for T in calc.
                                ! of saturation mixing ratio. (K)
     &, ttkm1(npnts)            ! Temporary store for T in layer
                                ! k-1 for use in freezing level
                                ! calc. for anvil. (K)

!  arrays required by grid calculations

      real ::                    &
       r2rho_th(npnts,nlev)      & ! radius**2 density theta lev (kg/m)
      ,r2rho(npnts,nlev)         & ! radius**2 density rho lev (kg/m)
      ,dr_across_th(npnts,nlev)  & ! thickness of theta levels (m)
      ,dr_across_rh(npnts,nlev)    ! thickness of rho levels (m)

!  arrays required by for identifying convective points

      integer :: index1(npnts)
      integer :: nconv_all         ! Total number of points convecting

! array required by tracers at end

      real :: limited_step(npnts)     ! Reduced step size for tracer
                                      ! mixing

      real :: step_test(npnts)        ! Work array used in reducing step

!
! Parameters
!

      real, parameter :: QMIN = 1.0E-8 ! Global minimum allowed Q

      real, parameter :: SAFETY_MARGIN = 1.0E-100 ! Small number used in
                                      ! tracer step reduction

!---------------------------------------------------------------------
!
! Loop counters
!

      integer :: i,j,k,ktra,idp,idp2

#if defined(SCMA)
      real :: rout(npnts)
      real :: r2out(npnts,nlev+1)
      character (len=15) :: sname
      character (len=100) :: lname
      character (len=11) :: units
! INOUT SCMop is declared in here  - for SCM diagnostics
#include "s_scmop.h"
#endif

!-----------------------------------------------------------------------
! Initialise output variables
!-----------------------------------------------------------------------
#if defined(SCMA)
! extra diagnostics
       write(6,*) 'GCONV5A : deep ',n_dp,' shall ',n_sh,' cong ',n_cg

#endif

!
! Initialise 3d rain variables
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
        wstar_dn(i)                = 0.0
        wstar_up(i)                = 0.0
        mb1(i)                     = 0.0
        mb2(i)                     = 0.0
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

      Do k = 1, nlev
        Do i = 1, n_cg
          rbuoy_p_out_cg(i,k)=0.0
          the_out_cg(i,k)=0.0
          thp_out_cg(i,k)=0.0
          qe_out_cg(i,k)=0.0
          qp_out_cg(i,k)=0.0
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
          cca(i,k)                 = 0.0
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
          Do i = 1,npnts
            qt(i,k) = q(i,k) + qcl(i,k) + qcf(i,k)
          End Do
        End Do
      Else
        Do k=1,nlev
          Do i = 1,npnts
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
            up_flux_half(i,k)        = 0.0
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
      If (flg_mf_deep) then
        Do k = 1,nlev
          Do i = 1,npnts
            mf_deep(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_mf_congest) then
        Do k = 1,nlev
          Do i = 1,npnts
            mf_congest(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_mf_shall) then
        Do k = 1,nlev
          Do i = 1,npnts
            mf_shall(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_mf_midlev) then
        Do k = 1,nlev
          Do i = 1,npnts
            mf_midlev(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dt_deep) then
        Do k = 1,nlev
          Do i = 1,npnts
            dt_deep(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dt_congest) then
        Do k = 1,nlev
          Do i = 1,npnts
            dt_congest(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dt_shall) then
        Do k = 1,nlev
          Do i = 1,npnts
            dt_shall(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dt_midlev) then
        Do k = 1,nlev
          Do i = 1,npnts
            dt_midlev(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dq_deep) then
        Do k = 1,nlev
          Do i = 1,npnts
            dq_deep(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dq_congest) then
        Do k = 1,nlev
          Do i = 1,npnts
            dq_congest(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dq_shall) then
        Do k = 1,nlev
          Do i = 1,npnts
            dq_shall(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dq_midlev) then
        Do k = 1,nlev
          Do i = 1,npnts
            dq_midlev(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_du_deep) then
        Do k = 1,nlev+1
          Do i = 1,npnts
            du_deep(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_du_congest) then
        Do k = 1,nlev+1
          Do i = 1,npnts
            du_congest(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_du_shall) then
        Do k = 1,nlev+1
          Do i = 1,npnts
            du_shall(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_du_midlev) then
        Do k = 1,nlev+1
          Do i = 1,npnts
            dv_midlev(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dv_deep) then
        Do k = 1,nlev+1
          Do i = 1,npnts
            dv_deep(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dv_congest) then
        Do k = 1,nlev+1
          Do i = 1,npnts
            dv_congest(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dv_shall) then
        Do k = 1,nlev+1
          Do i = 1,npnts
            dv_shall(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_dv_midlev) then
        Do k = 1,nlev+1
          Do i = 1,npnts
            dv_midlev(i,k)         = 0.0
          End Do
        End Do
      End If

      If (flg_wqt_flux_sh) then
        Do k = 1,nlev
          Do i = 1,npnts
            wqt_flux_sh(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_wql_flux_sh) then
        Do k = 1,nlev
          Do i = 1,npnts
            wql_flux_sh(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_wthetal_flux_sh) then
        Do k = 1,nlev
          Do i = 1,npnts
            wthetal_flux_sh(i,k)         = 0.0
          End Do
        End Do
      End If
      If (flg_wthetav_flux_sh) then
        Do k = 1,nlev
          Do i = 1,npnts
            wthetav_flux_sh(i,k)         = 0.0
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
        kterm_congest(i) = 0
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
      cg_on = 0
      md_on = 0
      dp_on = 0
      mdet_dp_on = 0
      mdet_md_on = 0
      mdet_sh_on = 0
      mdet_cg_on = 0
      sh_ent_on=0
      cg_ent_on=0
      md_ent_on=0
      dp_ent_on=0
      sh_new_termc=0
      cg_new_termc=0
      md_new_termc=0
      dp_new_termc=0
      sh_sdet_on=0
      cg_sdet_on=0
      md_sdet_on=0
      dp_sdet_on=0

      if(adapt  ==  0) then !default Gregory-Rowntree
        md_on = 0
        dp_on = 0
        mdet_dp_on = 0
        mdet_md_on = 0
        md_ent_on=0
        dp_ent_on=0
      else if(adapt  ==  1) then !HadGEM1a
        md_on = 1
        dp_on = 1
        mdet_dp_on = 1
        mdet_md_on = 1
        md_ent_on=0
        dp_ent_on=0
      else if(adapt  ==  2) then !test setting
        md_on = 1
        dp_on = 1
        mdet_dp_on = 1
        mdet_md_on = 1
        md_ent_on=1
        dp_ent_on=1
      else if(adapt  ==  3) then !operational
        md_on = 0
        dp_on = 1
        mdet_dp_on = 1
        mdet_md_on = 0
        md_ent_on=0
        dp_ent_on=0
      else if(adapt  ==  4) then  ! Possible HadGEM3 as option 1 plus shallow
        md_on = 1
        dp_on = 1
        mdet_dp_on = 1
        mdet_md_on = 1
        sh_on = 1
!        mdet_sh_on = 1    mixing detrainment not set on as has no impact.
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
        If (adapt == 4) Then
          sh_new_termc=1       ! use new termination condition for shallow
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
! Calculate saturation specific humidity/mixing ratio  lq_mix=.false.
!

! DEPENDS ON: qsat_mix
        Call QSAT_mix(qse(1,k),tt,pt,npnts,.false.)

      End Do  ! nlev

!-----------------------------------------------------------------------
! 1.0 DEEP Convection
! 1.1 Compress input variable arrays for deep convection scheme to
!     length n_dp (deep points only)
!-----------------------------------------------------------------------

      If (n_dp  >   0 .and. iconv_deep >  0 ) then
        j = 0
        Do i = 1,npnts
          If (cumulus_bl(i).and..not.l_shallow_bl(i).and.               &
     &                                       .not.l_congestus(i)) then
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
          W_MAX_dp(j)           = W_MAX(dpi(j))
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


        If (iconv_deep == 1) Then  ! G-R mass flux scheme

          If (l_mixing_ratio) Then  !  Conversion required

            Do k = 1,nlev
              Do j = 1,n_dp
                q_dp(j,k)   = q(dpi(j),k)  /(1.0+qt(dpi(j),k))
                qse_dp(j,k) = qse(dpi(j),k)/(1.0+qse(dpi(j),k))
                qcl_dp(j,k) = qcl(dpi(j),k)/(1.0+qt(dpi(j),k))
                qcf_dp(j,k) = qcf(dpi(j),k)/(1.0+qt(dpi(j),k))
              End Do
            End do

          Else         ! input is specific humidity therefore no problems

            Do k = 1,nlev
              Do j = 1,n_dp
                q_dp(j,k)   = q(dpi(j),k)
                qse_dp(j,k) = qse(dpi(j),k)
                qcl_dp(j,k) = qcl(dpi(j),k)
                qcf_dp(j,k) = qcf(dpi(j),k)
              End Do
            End do

          End If      ! Test on l_mixing_ratio

        Else        ! Future schemes to be code in mixing ratio

          If (l_mixing_ratio) Then   ! Input as required

            Do k = 1,nlev
              Do j = 1,n_dp
                q_dp(j,k)   = q(dpi(j),k)
                qse_dp(j,k) = qse(dpi(j),k)
                qcl_dp(j,k) = qcl(dpi(j),k)
                qcf_dp(j,k) = qcf(dpi(j),k)
              End Do
            End do

          Else      !  Conversion required

            Do k = 1,nlev
              Do j = 1,n_dp
                q_dp(j,k)   = q(dpi(j),k)  /(1.0-qt(dpi(j),k))
                qse_dp(j,k) = qse(dpi(j),k)/(1.0-qse(dpi(j),k))
                qcl_dp(j,k) = qcl(dpi(j),k)/(1.0-qt(dpi(j),k))
                qcf_dp(j,k) = qcf(dpi(j),k)/(1.0-qt(dpi(j),k))
              End Do
            End do

          End If    ! Test on l_mixing_ratio

        End If   ! test on deep scheme
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

      If (iconv_deep == 1) then  ! 4A like Gregory Rowntree deep conv

! DEPENDS ON: deep_conv
        Call DEEP_CONV(                                                 &
                       !IN
                       nbl,nlev,ntra,n_cca_lev,n_dp,trlev,              &
                       bland_dp, delthvu_dp, exner_layer_centres_dp,    &
                       exner_layer_boundaries_dp,flg_up_flx,            &
                       flg_entr_up, flg_detr_up, flg_dwn_flx,           &
                       flg_entr_dwn,flg_detr_dwn,flg_uw_deep,           &
                       flg_vw_deep, L_calc_dxek, L_q_interact,          &
                       L_tracer, ntml_dp, ntpar_dp,                     &
                       pstar_dp,p_layer_centres_dp,                     &
                       p_layer_boundaries_dp,                           &
                       z_theta_dp,z_rho_dp, r_theta_dp, r_rho_dp,       &
                       rho_theta_dp,rho_dp,r2rho_th_dp,r2rho_dp,        &
                       dr_across_th_dp,dr_across_rh_dp,                 &
                       q_dp,q1_sd_dp,                                   &
                       t1_sd_dp,th_dp,timestep,u_dp,v_dp,               &
                       uw0_dp,vw0_dp,w_max_dp,wstar_dp,                 &
                       zlcl_uv_dp,freeze_lev_dp,                        &
                       recip_pstar_dp,qse_dp,                           &
                       dp_on, mdet_dp_on, dp_ent_on,                    &
                       dp_sdet_on, dp_new_termc,                        &
                       !INOUT
                       bulk_cf_dp,cf_frozen_dp,cf_liquid_dp,            &
                       qcf_dp,qcl_dp,tracer_dp,                         &
                       !OUT
     &                 cape_out_dp,cclwp_dp,ccw_dp,cca_dp,              &
                       dbcfbydt_dp,dcffbydt_dp,dcflbydt_dp,             &
                       dqbydt_dp,dqcfbydt_dp,dqclbydt_dp,               &
                       dthbydt_dp,dubydt_dp,dvbydt_dp,dtrabydt_dp,      &
                       detrain_up_dp,detrain_dwn_dp,entrain_up_dp,      &
                       entrain_dwn_dp,                                  &
                       iccb_dp,icct_dp,                                 &
                       lcca_dp,lcclwp_dp,lcbase_dp,lctop_dp,            &
                       rain_dp,snow_dp, rain_3d_dp, snow_3d_dp,         &
                       up_flux_dp,                                      &
                       dwn_flux_dp,uw_deep_dp,vw_deep_dp,kterm_dp,      &
                       tcw_dp,cca_2d_dp,                                &
                       rbuoy_p_out_dp,the_out_dp,thp_out_dp,            &
                       qe_out_dp,qp_out_dp                              &
                       )

      Else       ! new turbulence scheme ?

!        Call DEEP_TURB_CONV(
!                       !IN
!     &                 nbl,nlev,ntra,n_cca_lev,n_dp,trlev,
!     &                 bland_dp,                
!     &                 exner_layer_centres_dp,
!     &                 )


      End If

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
        kterm_deep(dpi(i))        = kterm_dp(i)
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
          dqbydt(dpi(i),k)         = dqbydt_dp(i,k)
          dthbydt(dpi(i),k)        = dthbydt_dp(i,k)
          dqclbydt(dpi(i),k)       = dqclbydt_dp(i,k)
          dqcfbydt(dpi(i),k)       = dqcfbydt_dp(i,k)
          dcflbydt(dpi(i),k)       = dcflbydt_dp(i,k)
          dcffbydt(dpi(i),k)       = dcffbydt_dp(i,k)
          dbcfbydt(dpi(i),k)       = dbcfbydt_dp(i,k)
          ccw(dpi(i),k)            = ccw_dp(i,k)
          cca(dpi(i),k)            = cca_dp(i,k)
          rbuoy_p_out(dpi(i),k)    = rbuoy_p_out_dp(i,k)
          the_out(dpi(i),k)        = the_out_dp(i,k)
          thp_out(dpi(i),k)        = thp_out_dp(i,k)
          qe_out(dpi(i),k)         = qe_out_dp(i,k)
          qp_out(dpi(i),k)         = qp_out_dp(i,k)
        End Do
      End Do


      If (iconv_deep == 1) Then  !G-R mass flux scheme
        If (l_mixing_ratio) Then ! requires conversion

          If (l_q_interact) Then   ! PC2

            Do k = 1,nlev
              Do i = 1,n_dp
                dqtt  = dqbydt_dp(i,k)+dqclbydt_dp(i,k)+dqcfbydt_dp(i,k)
                dqt   = dqtt*timestep
                denom = 1.0/((1.0-qt(dpi(i),k))*(1.0-qt(dpi(i),k)-dqt))
                dqbydt(dpi(i),k)   = denom *                               &
                      ( dqbydt_dp(i,k)*(1.0-qt(dpi(i),k))+q_dp(i,k)*dqtt )
                dqclbydt(dpi(i),k) = denom *                               &
                    ( dqclbydt_dp(i,k)*(1.0-qt(dpi(i),k))+qcl_dp(i,k)*dqtt )
                dqcfbydt(dpi(i),k) = denom *                               &
                    ( dqcfbydt_dp(i,k)*(1.0-qt(dpi(i),k))+qcf_dp(i,k)*dqtt )
              End Do
            End do

          Else                     ! Not PC2
                                   ! No qcl and qcf increments anyway
            Do k = 1,nlev
              Do i = 1,n_dp
                dqt   = dqbydt_dp(i,k)*timestep
                denom = 1./((1.-qt(dpi(i),k))*(1.-qt(dpi(i),k)-dqt))
                dqbydt(dpi(i),k)    = dqbydt_dp(i,k)*denom
                dqclbydt(dpi(i),k)  = 0.0
                dqcfbydt(dpi(i),k)  = 0.0
              End Do
            End Do

          End If                   ! end test on PC2

        Else   ! output is specific humidity therefore no problems

          Do k = 1,nlev
            Do i = 1,n_dp
              dqbydt(dpi(i),k)   = dqbydt_dp(i,k)
              dqclbydt(dpi(i),k) = dqclbydt_dp(i,k)
              dqcfbydt(dpi(i),k) = dqcfbydt_dp(i,k)
            End Do
          End do

        End If      ! Test on l_mixing_ratio

      Else         ! Future schemes output as  mixing ratio

        If (l_mixing_ratio) Then ! ok
          Do k = 1,nlev
            Do i = 1,n_dp
              dqbydt(dpi(i),k)   = dqbydt_dp(i,k)
              dqclbydt(dpi(i),k) = dqclbydt_dp(i,k)
              dqcfbydt(dpi(i),k) = dqcfbydt_dp(i,k)
            End Do
          End Do

        Else   ! output is specific humidity

          If (l_q_interact) Then      ! PC2

            Do k = 1,nlev
              Do i = 1,n_dp
                dqtt= dqbydt_dp(i,k)+dqclbydt_dp(i,k)+dqcfbydt_dp(i,k)
                dqt = dqtt*timestep
                denom = 1.0/((1.0+qt(dpi(i),k))*(1.0+qt(dpi(i),k)+dqt))
                dqbydt(dpi(i),k)   = denom *                             &
                     (dqbydt_dp(i,k)*(1.0+qt(dpi(i),k))-q_dp(i,k)*dqtt)
                dqclbydt(dpi(i),k) = denom *                             &
                   (dqclbydt_dp(i,k)*(1.0+qt(dpi(i),k))-qcl_dp(i,k)*dqtt)
                dqcfbydt(dpi(i),k) = denom *                             &
                   (dqcfbydt_dp(i,k)*(1.0+qt(dpi(i),k))-qcf_dp(i,k)*dqtt)
              End Do
            End Do

          Else                         ! Not PC2
                                       ! No qcl and qcf increments
            Do k = 1,nlev
              Do i = 1,n_dp
                dqt   = dqbydt_dp(i,k)*timestep
                denom = 1.0/((1.0+qt(dpi(i),k))*(1.0+qt(dpi(i),k)+dqt))
                dqbydt(dpi(i),k)   =  dqbydt_dp(i,k)*denom
                dqclbydt(dpi(i),k) = 0.0
                dqcfbydt(dpi(i),k) = 0.0
              End Do
            End Do

          End If                       ! end test on PC2

        End If      ! Test on l_mixing_ratio
      End If        ! test on deep scheme

      If (flg_up_flx) then
        Do k = 1,nlev
          Do i = 1,n_dp
            up_flux(dpi(i),k)        = up_flux_dp(i,k)
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
      If (flg_mf_deep) then
        Do k = 1,nlev
          Do i = 1,n_dp
            mf_deep(dpi(i),k)        = up_flux_dp(i,k)
          End Do
        End Do
      End If
      If (flg_dt_deep) then
        Do k = 1,nlev
          Do i = 1,n_dp
            dt_deep(dpi(i),k)        = dthbydt_dp(i,k)
          End Do
        End Do
      End If
      If (flg_dq_deep) then
        Do k = 1,nlev
          Do i = 1,n_dp
            dq_deep(dpi(i),k)        = dqbydt_dp(i,k)
          End Do
        End Do
      End If
      If (L_mom) then
       If (flg_du_deep) then
        Do k = 1,nlev+1
          Do i = 1,n_dp
            du_deep(dpi(i),k)        = dubydt_dp(i,k)
          End Do
        End Do
       End If
       If (flg_dv_deep) then
        Do k = 1,nlev+1
          Do i = 1,n_dp
            dv_deep(dpi(i),k)        = dvbydt_dp(i,k)
          End Do
        End Do
       End If
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
! 1.4 CCRad - Apply CCA_dp/CCW_dp to full field ..... where ccw > 0.0
!-----------------------------------------------------------------------------

      If (Lcv_ccrad) Then

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

      If (n_sh  >   0 .and. iconv_shallow >  0) then
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
          ql_ad_sh(j)           = ql_ad(shi(j))
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
            p_layer_centres_sh(j,k)     = p_layer_centres(shi(j),k)
            p_layer_boundaries_sh(j,k)  = p_layer_boundaries(shi(j),k)
            exner_layer_centres_sh(j,k) = exner_layer_centres(shi(j),k)
            exner_layer_boundaries_sh(j,k)                              &
     &                               = exner_layer_boundaries(shi(j),k)
            r_theta_sh(j,k)             = r_theta(shi(j),k)
          End Do
        End Do

        Do k = 1,nlev
          Do j = 1,n_sh
            u_sh(j,k)           = u(shi(j),k)
            v_sh(j,k)           = v(shi(j),k)
            th_sh(j,k)          = th(shi(j),k)
            z_theta_sh(j,k)     = z_theta(shi(j),k)
            z_rho_sh(j,k)       = z_rho(shi(j),k)
            rho_sh(j,k)         = rho(shi(j),k)
            rho_th_sh(j,k)      = rho_theta(shi(j),k)
            r2rho_sh(j,k)       = r2rho(shi(j),k)
            r2rho_th_sh(j,k)    = r2rho_th(shi(j),k)
            r_rho_sh(j,k)       = r_rho(shi(j),k)
            dr_across_rh_sh(j,k)       = dr_across_rh(shi(j),k)
            dr_across_th_sh(j,k)       = dr_across_th(shi(j),k)
          End Do
        End Do
!
! moisture  - input depends on scheme
!
        If (iconv_shallow == 1) Then     ! G-R scheme
!
! G-R requires input of specific humidity
!
           If (l_mixing_ratio) Then

             Do k = 1,nlev
               Do j = 1,n_sh
                 q_sh(j,k)   = q(shi(j),k)  /(1.0+qt(shi(j),k))
                 qse_sh(j,k) = qse(shi(j),k)/(1.0+qse(shi(j),k))
                 qcl_sh(j,k) = qcl(shi(j),k)/(1.0+qt(shi(j),k))
                 qcf_sh(j,k) = qcf(shi(j),k)/(1.0+qt(shi(j),k))
               End Do
             End Do

           Else       ! Input is specific humidity therefore no problems

             Do k = 1,nlev
               Do j = 1,n_sh
                 q_sh(j,k)   = q(shi(j),k)
                 qse_sh(j,k) = qse(shi(j),k)
                 qcl_sh(j,k) = qcl(shi(j),k)
                 qcf_sh(j,k) = qcf(shi(j),k)
               End Do
             End Do

           End If     ! Test on l_mixing_ratio

         Else         ! turbulence scheme requires mixing ratio
!
! turbulence scheme requires mixing ratios as input
!
           If (l_mixing_ratio) Then

             Do k = 1,nlev
               Do j = 1,n_sh
                 q_sh(j,k)   = q(shi(j),k)
                 qse_sh(j,k) = qse(shi(j),k)
                 qcl_sh(j,k) = qcl(shi(j),k)
                 qcf_sh(j,k) = qcf(shi(j),k)
               End Do
             End Do

           Else       ! Require conversion

             Do k = 1,nlev
               Do j = 1,n_sh
                 q_sh(j,k)     = q(shi(j),k)/(1.0 -qt(shi(j),k))

! check for negative values
                 If (q_sh(j,k) <  0.0) Then
                   write(6,*) 'Problem q_mix -ve ',q_sh(j,k),' j,k',j,k
                   q_sh(j,k) = 0.0
                 End If
                 If (qse(shi(j),k) <  1.0) Then
                   qse_sh(j,k) = qse(shi(j),k)/(1.-qse(shi(j),k))
                 Else
! Appears to be activated more frequently at 70 levels near the top
! of the model
! Only print if q actually a larger value ?
                   If (q_sh(j,k) > 5.e-6) then
                     write(6,*) 'Problem qsat 1 ',qse(shi(j),k),' j,k',j,k  &
                         ,p_layer_centres(shi(j),k),th(shi(j),k),q_sh(j,k)
                   End If
                   qse_sh(j,k) = 10000.    ! very big value
                 End If     ! Test on qse

! Currently turbulence shallow scheme does not alter or use qcl or qcf
! so no checks are being applied to the values.
                 qcl_sh(j,k) = qcl(shi(j),k)/(1.0-qt(shi(j),k))
                 qcf_sh(j,k) = qcf(shi(j),k)/(1.0-qt(shi(j),k))
               End Do  ! j
             End Do    ! k

           End If      ! Test on l_mixing_ratio

         End If        ! test on scheme

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

       If (iconv_shallow == 1) then !  G-R mass flux type scheme

! DEPENDS ON: shallow_conv
         Call SHALLOW_CONV(                                             &
                           !IN
     &                     nbl,nlev,ntra,n_cca_lev,n_sh,trlev,          &
     &                     bland_sh,                                    &
     &                     delthvu_sh,exner_layer_centres_sh,           &
     &                     exner_layer_boundaries_sh,                   &
     &                     flg_up_flx,                                  &
     &                     flg_entr_up, flg_detr_up, flg_dwn_flx,       &
     &                     flg_entr_dwn,flg_detr_dwn,flg_uw_shall,      &
     &                     flg_vw_shall, L_calc_dxek,                   &
     &                     L_q_interact,L_tracer,ntml_sh,ntpar_sh,      &
     &                     pstar_sh,p_layer_centres_sh,                 &
     &                     p_layer_boundaries_sh,z_theta_sh,z_rho_sh,   &
     &                     r2rho_th_sh,dr_across_th_sh,                 &
     &                     q_sh,q1_sd_sh,t1_sd_sh,th_sh,timestep,       &
     &                     u_sh,v_sh,uw0_sh,vw0_sh,wstar_sh,wthvs_sh,   &
     &                     zlcl_uv_sh,ztop_uv_sh,freeze_lev_sh,         &
     &                     recip_pstar_sh,qse_sh,                       &
     &                     sh_on, mdet_sh_on, sh_ent_on,                &
     &                     sh_sdet_on, sh_new_termc,                    &
                           !INOUT
     &                     bulk_cf_sh,cf_frozen_sh,cf_liquid_sh,        &
     &                     qcf_sh,qcl_sh,tracer_sh,                     &
                           !OUT
     &                     cape_out_sh,cclwp_sh,ccw_sh,cca_sh,          &
     &                     dbcfbydt_sh,dcffbydt_sh,dcflbydt_sh,         &
     &                     dqbydt_sh,dqcfbydt_sh,dqclbydt_sh,           &
     &                     dthbydt_sh,dubydt_sh,dvbydt_sh,dtrabydt_sh,  &
     &                     detrain_up_sh,detrain_dwn_sh,entrain_up_sh,  &
     &                     entrain_dwn_sh,                              &
     &                     iccb_sh,icct_sh,                             &
     &                     lcca_sh,lcclwp_sh,lcbase_sh,lctop_sh,        &
     &                     rain_sh,snow_sh, rain_3d_sh, snow_3d_sh,     &
     &                     up_flux_sh,                                  &
     &                     dwn_flux_sh,uw_shall_sh,vw_shall_sh,         &
     &                     tcw_sh,cca_2d_sh,                            &
     &                     rbuoy_p_out_sh,the_out_sh,thp_out_sh,        &
     &                     qe_out_sh,qp_out_sh                          &
     &                     )


       Else if (iconv_shallow == 2 .or. iconv_shallow == 3) then

! turbulence based shallow convection scheme
!              2 - version with no shallow precipitation
!              3 - version with shallow precipitation
         If (iconv_shallow == 2) then
           ishall_precip = 0
         Else
           ishall_precip = 1
         End If

! Information only required by turbulence version

         Do j = 1,n_sh
           zlcl_sh(j)            = zlcl(shi(j))
         End Do


! DEPENDS ON: shallow_turb_conv
         Call SHALLOW_TURB_CONV(                                        &
                   !IN
                            nlev, n_sh, ntra,trlev                      &
      ,                     ntml_sh, ntpar_sh                           &
      ,                     ishall_precip, L_calc_dxek, L_q_interact    &
      ,                     L_tracer, flg_uw_shall, flg_vw_shall        &
      ,                     bland_sh, timestep                          &
      ,                     delthvu_sh,ql_ad_sh,uw0_sh,vw0_sh           &
      ,                     wstar_sh,wthvs_sh,zlcl_sh,zlcl_uv_sh        &
      ,                     ztop_uv_sh                                  &
      ,                     pstar_sh,p_layer_centres_sh                 &
      ,                     p_layer_boundaries_sh                       &
      ,                     exner_layer_centres_sh                      &
      ,                     exner_layer_boundaries_sh                   &
      ,                     z_theta_sh,z_rho_sh,r_rho_sh,r_theta_sh     &
      ,                     rho_sh,rho_th_sh,r2rho_sh,r2rho_th_sh       &
      ,                     dr_across_rh_sh,dr_across_th_sh             &
      ,                     q_sh,th_sh,u_sh,v_sh                        &
      ,                     qse_sh                                      &
                   !INOUT
      ,                     bulk_cf_sh,cf_frozen_sh,cf_liquid_sh        &
      ,                     qcf_sh,qcl_sh,tracer_sh                     &
                   !OUT
      ,                     cape_out_sh,cclwp_sh,ccw_sh,cca_sh          &
      ,                     dbcfbydt_sh,dcffbydt_sh,dcflbydt_sh         &
      ,                     dqbydt_sh,dqcfbydt_sh,dqclbydt_sh           &
      ,                     dthbydt_sh,dubydt_sh,dvbydt_sh,dtrabydt_sh  &
      ,                     iccb_sh,icct_sh,lcca_sh,lcbase_sh,lctop_sh  &
      ,                     rain_sh,snow_sh,up_flux_sh                  &
      ,                     uw_shall_sh,vw_shall_sh                     &
      ,                     wqt_sh,wthetal_sh,wthetav_sh,wql_sh         &
      ,                     wstar_up_sh,mb1_sh,mb2_sh                   &
      ,                     tcw_sh,cca_2d_sh                            &
                           )

        End If        ! test on shallow convection scheme

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
        wstar_dn(shi(i))           = wstar_sh(i)  ! wstar_dn
        wstar_up(shi(i))           = wstar_up_sh(i)
        mb1(shi(i))                = mb1_sh(i)
        mb2(shi(i))                = mb2_sh(i)
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
          cca(shi(i),k)            = cca_sh(i,k)

          rbuoy_p_out(shi(i),k)    = rbuoy_p_out_sh(i,k)
          the_out(shi(i),k)        = the_out_sh(i,k)
          thp_out(shi(i),k)        = thp_out_sh(i,k)
          qe_out(shi(i),k)         = qe_out_sh(i,k)
          qp_out(shi(i),k)         = qp_out_sh(i,k)

        End Do
      End Do
!
! moisture  - output depends on scheme
!
        If (iconv_shallow == 1) then     ! G-R scheme
!
! G-R requires input of specific humidity
!
          If (l_mixing_ratio) Then ! requires conversion

            If (l_q_interact) Then    ! PC2

              Do k = 1,nlev
                Do i = 1,n_sh
                  dqtt  = dqbydt_sh(i,k)+dqclbydt_sh(i,k)+dqcfbydt_sh(i,k)
                  dqt   = dqtt*timestep
                  denom = 1.0/((1.0-qt(shi(i),k))*(1.0-qt(shi(i),k)-dqt))
                  dqbydt(shi(i),k)    = denom *                            &
                      ( dqbydt_sh(i,k)*(1.0-qt(shi(i),k))+q_sh(i,k)*dqtt )
                  dqclbydt(shi(i),k)  = denom *                            &
                    ( dqclbydt_sh(i,k)*(1.0-qt(shi(i),k))+qcl_sh(i,k)*dqtt )
                  dqcfbydt(shi(i),k)  = denom *                            &
                    ( dqcfbydt_sh(i,k)*(1.0-qt(shi(i),k))+qcf_sh(i,k)*dqtt )
                End Do
              End Do

            Else                      ! Not PC2
                                      ! No qcl and qcf increments anyway
              Do k = 1,nlev
                Do i = 1,n_sh
                  dqt = dqbydt_sh(i,k)*timestep
                  denom = 1.0/((1.0-qt(shi(i),k))*(1.0-qt(shi(i),k)-dqt))
                  dqbydt(shi(i),k)    = dqbydt_sh(i,k)* denom
                  dqclbydt(shi(i),k)  = 0.0
                  dqcfbydt(shi(i),k)  = 0.0
                End Do
              End Do
            End If                    ! End test on PC2

           Else   ! output is specific humidity therefore no problems

             Do k = 1,nlev
               Do i = 1,n_sh
                 dqbydt(shi(i),k)   = dqbydt_sh(i,k)
                 dqclbydt(shi(i),k) = dqclbydt_sh(i,k)
                 dqcfbydt(shi(i),k) = dqcfbydt_sh(i,k)
               End Do
             End Do

           End If     ! Test on l_mixing_ratio

         Else         ! turbulence scheme requires mixing ratio
!
! turbulence scheme outputs mixing ratios
!
           If (l_mixing_ratio) Then

             Do k = 1,nlev
               Do i = 1,n_sh
                 dqbydt(shi(i),k)   = dqbydt_sh(i,k)
                 dqclbydt(shi(i),k) = dqclbydt_sh(i,k)
                 dqcfbydt(shi(i),k) = dqcfbydt_sh(i,k)
               End Do
             End Do

           Else    ! Requires conversion increment to mixing ratio not q
             If (l_q_interact) then     ! PC2
! At present turbulence scheme has no dqcl and dqcf so dqt=dq
               Do k = 1,nlev
                 Do i = 1,n_sh
                   dqtt = dqbydt_sh(i,k)
                   dqt = dqtt*timestep
                   denom = 1./((1.+qt(shi(i),k))*(1.+qt(shi(i),k)+dqt))
                   dqbydt(shi(i),k)   = denom *                           &
                       (dqbydt_sh(i,k)*(1.+qt(shi(i),k))-q_sh(i,k)*dqtt)
                   dqclbydt(shi(i),k) = denom *                           &
                     (dqclbydt_sh(i,k)*(1.+qt(shi(i),k))-qcl_sh(i,k)*dqtt)
                   dqcfbydt(shi(i),k) = denom *                           &
                     (dqcfbydt_sh(i,k)*(1.+qt(shi(i),k))-qcf_sh(i,k)*dqtt)
                 End do
               End do

             Else                       ! Not PC2
                                        !  No qcl and qcf increments
               Do k = 1,nlev
                 Do i = 1,n_sh
                   dqt   = dqbydt_sh(i,k)*timestep
                   denom = 1.0/((1.0+qt(shi(i),k))*(1.0+qt(shi(i),k)+dqt))
                   dqbydt(shi(i),k)   =  dqbydt_sh(i,k)*denom
                   dqclbydt(shi(i),k) = 0.0
                   dqcfbydt(shi(i),k) = 0.0
                 End Do
               End do
             End If                     ! Test on PC2

           End If       ! Test on l_mixing_ratio

         End If        ! test on scheme
!
! diagnostics output
!
      If (flg_up_flx) then
        Do k = 1,nlev
          Do i = 1,n_sh
            up_flux(shi(i),k)        = up_flux_sh(i,k)
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

      If (flg_wqt_flux_sh) then
        Do k = 1,nlev
          Do i = 1,n_sh
            wqt_flux_sh(shi(i),k)       = wqt_sh(i,k)
          End Do
        End Do
      End If
      If (flg_wql_flux_sh) then
        Do k = 1,nlev
          Do i = 1,n_sh
            wql_flux_sh(shi(i),k)       = wql_sh(i,k)
          End Do
        End Do
      End If
      If (flg_wthetal_flux_sh) then
        Do k = 1,nlev
          Do i = 1,n_sh
            wthetal_flux_sh(shi(i),k)    = wthetal_sh(i,k)
          End Do
        End Do
      End If
      If (flg_wthetav_flux_sh) then
        Do k = 1,nlev
          Do i = 1,n_sh
            wthetav_flux_sh(shi(i),k)    = wthetav_sh(i,k)
          End Do
        End Do
      End If

      If (flg_mf_shall) then
        Do k = 1,nlev
          Do i = 1,n_sh
            mf_shall(shi(i),k)        = up_flux_sh(i,k)
          End Do
        End Do
      End If
      If (flg_dt_shall) then
        Do k = 1,nlev
          Do i = 1,n_sh
            dt_shall(shi(i),k)        = dthbydt_sh(i,k)
          End Do
        End Do
      End If
      If (flg_dq_shall) then
        Do k = 1,nlev
          Do i = 1,n_sh
            dq_shall(shi(i),k)        = dqbydt_sh(i,k)
          End Do
        End Do
      End If
      If (L_mom) then
       If (flg_du_shall) then
        Do k = 1,nlev+1
          Do i = 1,n_sh
            du_shall(shi(i),k)        = dubydt_sh(i,k)
          End Do
        End Do
       End If
       If (flg_dv_shall) then
        Do k = 1,nlev+1
          Do i = 1,n_sh
            dv_shall(shi(i),k)        = dvbydt_sh(i,k)
          End Do
        End Do
       End If
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
! 2.4 CCRad - Apply CCA/CCW to full field .....(on shallow points)
!-----------------------------------------------------------------------------

        If (lcv_ccrad) Then

          Do k=1, nlev
            Do i=1, n_sh
              If (iccb_sh(i) /= 0) Then ! Shallow convection occured
                cca(shi(i),k) = cca_sh(i,k)
                ccw(shi(i),k) = ccw_sh(i,k)
              End If
            End Do      ! i
          End Do      ! k

        End If      ! lcv_ccrad

#if defined(SCMA)
! extra 5A diagnostics

      Else               ! SCM diagnostics if no shallow

! output zero
      Do i=1,npnts
        rout(i)=0.0
      End Do
      Do k=1,nlev+1
      Do i=1,npnts
        r2out(i,k)=0.0
      End Do
      End Do

        If (iconv_shallow == 2 .or. iconv_shallow == 3) then

!  cloud base fluxes etc
      sname='mb_cb'
      lname='cloud base mass flux'
      units='m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

      sname='mb_new_cb'
      lname='revised cloud base mass flux'
      units='m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

      sname='wstar_up'
      lname='wstar_up'
      units='m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

 
      sname='wthetal_cb'
      lname='wthetal at cloud base'
      units='Km/s'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)
      sname='wqt_cb'
      lname='wqt at cloud base'
      units='Kg/kg m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

      sname='wql_cb'
      lname='wql at cloud base'
      units='Kg/kg m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

      sname='wql_inv'
      lname='wql at inversion'
      units='Kg/kg m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

      sname='wqr_inv'
      lname='wqr at inversion'
      units='Kg/kg m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

      sname='dthetav_cb'
      lname='dthetav across cloud base'
      units='K'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

      sname='dtheta_cb'
      lname='dtheta across cloud base'
      units='K'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)
      sname='drv_cb'
      lname='drv across cloud base'
      units='kg/kg'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

      sname='wthetav_m_cb'
      lname='wthetav on lower side of cloud base'
      units='Km/s'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

      sname='wqt_inv'
      lname='wqt at inversion'
      units='kg/kg m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

      sname='wthetal_inv'
      lname='wthetal at inversion'
      units='kg/kg m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(rout,sname,lname,units,t_inst,d_point              &
         , default_streams,'',RoutineName)

      sname='dwqldz'
      lname='dql/dz  on model levels'
      units='Kg/kg s'
! DEPENDS ON: scmoutput
      call SCMoutput(r2out,                                             &
                     sname,lname,units,t_inst,d_wet,default_streams,'', &
                     RoutineName)

      sname='qlup'
      lname='qlup  on model levels'
      units='Kg/kg '
! DEPENDS ON: scmoutput
      call SCMoutput(r2out,                                             &
                     sname,lname,units,t_inst,d_wet,default_streams,'', &
                     RoutineName)
      sname='wup'
      lname='wup  on model levels'
      units='m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(r2out,                                             &
                     sname,lname,units,t_inst,d_wet,default_streams,'', &
                     RoutineName)

      sname='wthetav'
      lname='wthetav  on model levels'
      units='K m/s '
! DEPENDS ON: scmoutput
      call SCMoutput(r2out,                                             &
                     sname,lname,units,t_inst,d_wet,default_streams,'', &
                     RoutineName)

      sname='wql'
      lname='wql flux in convective cloud'
      units=' '
! DEPENDS ON: scmoutput
      call SCMoutput(r2out,                                             &
                     sname,lname,units,t_inst,d_wet,default_streams,'', &
                     RoutineName)

      sname='dqsat_dt'
      lname='dqsat_dt'
      units='K m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(r2out,                                             &
                     sname,lname,units,t_inst,d_wet,default_streams,'', &
                     RoutineName)
      sname='wthetal'
      lname='wthetal'
      units='K m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(r2out,                                             &
                     sname,lname,units,t_inst,d_wet,default_streams,'', &
                     RoutineName)
      sname='wqt'
      lname='wqt'
      units='Kg/kg m/s'
! DEPENDS ON: scmoutput
      call SCMoutput(r2out,                                             &
                     sname,lname,units,t_inst,d_wet,default_streams,'', &
                     RoutineName)
      sname='wqr'
      lname='wqr on uv levels'
      units='Kg/kgm/s'
! DEPENDS ON: scmoutput
      call SCMoutput(r2out,                                             &
                     sname,lname,units,t_inst,d_wet,default_streams,'', &
                     RoutineName)
      sname='precip_prod'
      lname='precip production on theta levels'
      units='/s'
! DEPENDS ON: scmoutput
      call SCMoutput(r2out,                                             &
                     sname,lname,units,t_inst,d_wet,default_streams,'', &
                     RoutineName)

        End If ! type of shallow convection
#endif


      End If ! n_sh > 0

!-----------------------------------------------------------------------
! 3.0 CONGESTUS Convection   - call depending on switches
!-----------------------------------------------------------------------

      If (iconv_congestus >  0) then

        If (n_cg  >   0) then   ! there must be congestus points

!-----------------------------------------------------------------------
! 3.1 Compress input variable arrays for congestus convection scheme to
!     length n_cg (congestus points only)
!-----------------------------------------------------------------------

          j = 0
          Do i = 1,npnts
            If (cumulus_bl(i).and.l_congestus(i)) then
              j                        = j+1
              cgi(j)                   = i
            End If
          End Do
!
! In only variables
!
          Do j = 1,n_cg
            bland_cg(j)           = bland(cgi(j))
            delthvu_cg(j)         = delthvu_bl(cgi(j))
            ntml_cg(j)            = ntml(cgi(j))
            ntpar_cg(j)           = ntpar(cgi(j))
            pstar_cg(j)           = pstar(cgi(j))
            recip_pstar_cg(j)     = recip_pstar(cgi(j))
            q1_sd_cg(j)           = q1_sd(cgi(j))
            t1_sd_cg(j)           = t1_sd(cgi(j))
            uw0_cg(j)             = uw0(cgi(j))
            vw0_cg(j)             = vw0(cgi(j))
            wstar_cg(j)           = wstar(cgi(j))
            wthvs_cg(j)           = wthvs(cgi(j))
            zlcl_uv_cg(j)         = zlcl_uv(cgi(j))
            ztop_uv_cg(j)         = ztop_uv(cgi(j))
            freeze_lev_cg(j)      = freeze_lev(cgi(j))
            kterm_cg(j) = 0     ! initialise array Will only get set
                                ! if call deep scheme
          End Do

          Do k = 0,nlev
            Do j = 1,n_cg
              p_layer_centres_cg(j,k)    = p_layer_centres(cgi(j),k)
              p_layer_boundaries_cg(j,k) = p_layer_boundaries(cgi(j),k)
              exner_layer_centres_cg(j,k)= exner_layer_centres(cgi(j),k)
              exner_layer_boundaries_cg(j,k)                            &
     &                                = exner_layer_boundaries(cgi(j),k)
              r_theta_cg(j,k)            = r_theta(cgi(j),k)
            End Do
          End Do

          Do k = 1,nlev
            Do j = 1,n_cg
              u_cg(j,k)           = u(cgi(j),k)
              v_cg(j,k)           = v(cgi(j),k)
              th_cg(j,k)          = th(cgi(j),k)
              z_theta_cg(j,k)     = z_theta(cgi(j),k)
              z_rho_cg(j,k)       = z_rho(cgi(j),k)
              rho_cg(j,k)         = rho(cgi(j),k)
              rho_th_cg(j,k)      = rho_theta(cgi(j),k)
              r2rho_cg(j,k)       = r2rho(cgi(j),k)
              r2rho_th_cg(j,k)    = r2rho_th(cgi(j),k)
              r_rho_cg(j,k)       = r_rho(cgi(j),k)
              dr_across_rh_cg(j,k)       = dr_across_rh(cgi(j),k)
              dr_across_th_cg(j,k)       = dr_across_th(cgi(j),k)
            End Do
          End Do
!
! moisture  - input depends on scheme
!
        If (iconv_congestus == 1) Then     ! G-R scheme
!
! G-R requires input of specific humidity
!
           If (l_mixing_ratio) Then

             Do k = 1,nlev
               Do j = 1,n_cg
                 q_cg(j,k)   = q(cgi(j),k)  /(1.0+qt(cgi(j),k))
                 qse_cg(j,k) = qse(cgi(j),k)/(1.0+qse(cgi(j),k))
                 qcl_cg(j,k) = qcl(cgi(j),k)/(1.0+qt(cgi(j),k))
                 qcf_cg(j,k) = qcf(cgi(j),k)/(1.0+qt(cgi(j),k))
               End do
             End do

           Else   ! input is specific humidity therefore no problems

             Do k = 1,nlev
               Do j = 1,n_cg
                 q_cg(j,k)   = q(cgi(j),k)
                 qse_cg(j,k) = qse(cgi(j),k)
                 qcl_cg(j,k) = qcl(cgi(j),k)
                 qcf_cg(j,k) = qcf(cgi(j),k)
               End Do
             End Do

           End If     ! l_mixing_ratio

         Else         ! turbulence scheme requires mixing ratio
!
! turbulence scheme requires mixing ratios as input
!
           If (l_mixing_ratio) Then

             Do k = 1,nlev
               Do j = 1,n_cg
                 q_cg(j,k)   = q(cgi(j),k)
                 qse_cg(j,k) = qse(cgi(j),k)
                 qcl_cg(j,k) = qcl(cgi(j),k)
                 qcf_cg(j,k) = qcf(cgi(j),k)
               End Do
             End Do

           Else       ! Require conversion

             Do k = 1,nlev
               Do j = 1,n_cg

                 q_cg(j,k)     = q(cgi(j),k)/(1.0 -qt(cgi(j),k))

! check for negative values
                 If (q_cg(j,k) <  0.0) Then
                   write(6,*) 'Problem q_mix -ve ',q_cg(j,k),' j,k',j,k
                   q_cg(j,k) = 0.0
                 End If
                 If (qse(cgi(j),k) <  1.0) Then
                   qse_cg(j,k) = qse(cgi(j),k)/(1.-qse(cgi(j),k))
                 Else
                   write(6,*) 'Problem qsat 1 ',qse(cgi(j),k),' j,k',j,k &
                         ,p_layer_centres(cgi(j),k),th(cgi(j),k),q_cg(j,k)
                   qse_cg(j,k) = 10000.    ! very big value
                 End If

! No checks on qcl and qcf as not expected to be used or altered.

                 qcl_cg(j,k) = qcl(cgi(j),k)/(1.0-qt(cgi(j),k))
                 qcf_cg(j,k) = qcf(cgi(j),k)/(1.0-qt(cgi(j),k))
               End Do
             End Do

           End If      ! Test on l_mixing_ratio

         End If        ! Test on scheme (iconv_congestus)

!
! In/out variables
!
          Do k = 1,nlev
            Do j = 1,n_cg
              cf_liquid_cg(j,k)   = cf_liquid(cgi(j),k)
              cf_frozen_cg(j,k)   = cf_frozen(cgi(j),k)
              bulk_cf_cg(j,k)     = bulk_cf(cgi(j),k)
            End Do
          End Do

          If (L_tracer) then
            Do ktra = 1,ntra
              Do k = 1,trlev
                Do j = 1,n_cg
                  tracer_cg(j,k,ktra)  = tracer(cgi(j),k,ktra)
                End Do
              End Do
              Do k = 1,nlev
                Do j = 1,n_cg
                  dtrabydt_cg(j,k,ktra)  = 0.0
                End Do
              End Do
            End Do
          End If

!-----------------------------------------------------------------------
! 3.2 Call congestus convection code
!-----------------------------------------------------------------------

          If (iconv_congestus == 1) then   ! G-R type scheme

! new mass flux routine (at present a copy of shallow scheme)

! DEPENDS ON: congest_conv
            Call congest_conv(nbl,nlev,ntra,n_cca_lev,n_cg,trlev        &
     &,                       bland_cg                                  &
     &,                       delthvu_cg,exner_layer_centres_cg         &
     &,                       exner_layer_boundaries_cg,flg_up_flx      &
     &,                       flg_entr_up, flg_detr_up, flg_dwn_flx     &
     &,                       flg_entr_dwn,flg_detr_dwn,flg_uw_shall    &
     &,                       flg_vw_shall, L_calc_dxek, L_q_interact   &
     &,                       L_tracer, ntml_cg, ntpar_cg               &
     &,                       pstar_cg,p_layer_centres_cg               &
     &,                       p_layer_boundaries_cg,z_theta_cg,z_rho_cg &
     &,                       r2rho_th_cg, dr_across_th_cg              &
     &,                       q_cg,q1_sd_cg,t1_sd_cg,th_cg,timestep     &
     &,                       u_cg,v_cg,uw0_cg,vw0_cg,wstar_cg,wthvs_cg &
     &,                       zlcl_uv_cg,ztop_uv_cg,freeze_lev_cg       &
     &,                       recip_pstar_cg,qse_cg                     &
     &,                       cg_on, mdet_cg_on, cg_ent_on              &
     &,                       cg_sdet_on, cg_new_termc                  &
     &,                     bulk_cf_cg,cf_frozen_cg,cf_liquid_cg,qcf_cg &
     &,                     qcl_cg,tracer_cg,cape_out_cg,cclwp_cg,ccw_cg&
     &,                     cca_cg                                      &
     &,                     dbcfbydt_cg,dcffbydt_cg,dcflbydt_cg         &
     &,                     dqbydt_cg,dqcfbydt_cg,dqclbydt_cg,dthbydt_cg&
     &,                     dubydt_cg,dvbydt_cg,dtrabydt_cg             &
     &,                     detrain_up_cg,detrain_dwn_cg                &
     &,                     entrain_up_cg,entrain_dwn_cg                &
     &,                     iccb_cg,icct_cg,lcca_cg,lcclwp_cg,lcbase_cg &
     &,                     lctop_cg,rain_cg,snow_cg,up_flux_cg         &
     &,                     dwn_flux_cg,uw_shall_cg,vw_shall_cg,kterm_cg&
     &,                     tcw_cg,cca_2d_cg                            &
     &,                     rbuoy_p_out_cg,the_out_cg,thp_out_cg        &
     &,                     qe_out_cg,qp_out_cg )


!         Else If (iconv_congestus == 2) then  ! turbulence scheme
!
!            Call TURB_CONGESTUS_CONV(
!
!        Not available yet - still to be written
!
          End If   ! test of type of congestus scheme

!-----------------------------------------------------------------------
! 3.3 Write data from congestus convection points to full arrays
!-----------------------------------------------------------------------
!
          Do i = 1,n_cg
            cape_out(cgi(i))           = cape_out_cg(i)
            cclwp(cgi(i))              = cclwp_cg(i)
            iccb(cgi(i))               = iccb_cg(i)
            icct(cgi(i))               = icct_cg(i)
            lcca(cgi(i))               = lcca_cg(i)
            lcbase(cgi(i))             = lcbase_cg(i)
            lctop(cgi(i))              = lctop_cg(i)
            lcclwp(cgi(i))             = lcclwp_cg(i)
            rain(cgi(i))               = rain_cg(i)
            snow(cgi(i))               = snow_cg(i)
            precip_cong(cgi(i))        = rain_cg(i) + snow_cg(i)
            wstar_dn(cgi(i))           = wstar_cg(i)  ! wstar_dn
            wstar_up(cgi(i))           = wstar_up_cg(i)
            mb1(cgi(i))                = mb1_cg(i)
            mb2(cgi(i))                = mb2_cg(i)
! may be required if congestus not forced to stop at ntpar
            kterm_congest(cgi(i))      = kterm_cg(i)

! required at present to replicate old code
            tcw_md(cgi(i))            = tcw_cg(i)
            iccb_md(cgi(i))           = iccb_cg(i)
            icct_md(cgi(i))           = icct_cg(i)
            cca_2d_md(cgi(i))         = cca_2d_cg(i)

          End Do

          If (L_mom) then
            Do k=1,nlev+1
              Do i = 1,n_cg
                dubydt(cgi(i),k)         = dubydt_cg(i,k)
                dvbydt(cgi(i),k)         = dvbydt_cg(i,k)
              End Do
            End Do
          End If

          Do k = 1,nlev
            Do i = 1,n_cg
              dthbydt(cgi(i),k)        = dthbydt_cg(i,k)
              dcflbydt(cgi(i),k)       = dcflbydt_cg(i,k)
              dcffbydt(cgi(i),k)       = dcffbydt_cg(i,k)
              dbcfbydt(cgi(i),k)       = dbcfbydt_cg(i,k)
              ccw(cgi(i),k)            = ccw_cg(i,k)
              cca(cgi(i),k)            = cca_cg(i,k)
              rbuoy_p_out(cgi(i),k)    = rbuoy_p_out_cg(i,k)
              the_out(cgi(i),k)        = the_out_cg(i,k)
              thp_out(cgi(i),k)        = thp_out_cg(i,k)
              qe_out(cgi(i),k)         = qe_out_cg(i,k)
              qp_out(cgi(i),k)         = qp_out_cg(i,k)
            End Do
          End Do
!
! moisture  - output depends on scheme
!
        If (iconv_congestus == 1) then     ! G-R scheme
!
! Outputs specific humidity increments
!
          If (l_mixing_ratio) Then ! requires conversion

            If (l_q_interact) Then    ! PC2

              Do k = 1,nlev
                Do i = 1,n_cg
                  dqtt  = dqbydt_cg(i,k)+dqclbydt_cg(i,k)+dqcfbydt_cg(i,k)
                  dqt   = dqtt*timestep
                  denom = 1.0/((1.0-qt(cgi(i),k))*(1.0-qt(cgi(i),k)-dqt))
                  dqbydt(cgi(i),k)   = denom *                              &
                     ( dqbydt_cg(i,k)*(1.0-qt(cgi(i),k))+q_cg(i,k)*dqtt )
                  dqclbydt(cgi(i),k) = denom *                              &
                   ( dqclbydt_cg(i,k)*(1.0-qt(cgi(i),k))+qcl_cg(i,k)*dqtt )
                  dqcfbydt(cgi(i),k) = denom *                              &
                   ( dqcfbydt_cg(i,k)*(1.0-qt(cgi(i),k))+qcf_cg(i,k)*dqtt )
                End Do
              End Do

            Else                      ! Not PC2
                                      ! No qcl and qcf increments anyway
              Do k = 1,nlev
                Do i = 1,n_cg
                  dqt = dqbydt_cg(i,k)*timestep
                  denom = 1.0/((1.0-qt(cgi(i),k))*(1.0-qt(cgi(i),k)-dqt))
                  dqbydt(cgi(i),k)    =  dqbydt_cg(i,k)*denom
                  dqclbydt(cgi(i),k)  = 0.0
                  dqcfbydt(cgi(i),k)  = 0.0
                End Do
              End Do
            End If                     ! End test on PC2

          Else   ! output is specific humidity therefore no problems

             Do k = 1,nlev
               Do i = 1,n_cg
                 dqbydt(cgi(i),k)   = dqbydt_cg(i,k)
                 dqclbydt(cgi(i),k) = dqclbydt_cg(i,k)
                 dqcfbydt(cgi(i),k) = dqcfbydt_cg(i,k)
               End Do
             End Do

           End If     ! l_mixing_ratio

         Else         ! turbulence scheme
!
! turbulence scheme outputs mixing ratios
!
           If (l_mixing_ratio) Then

             Do k = 1,nlev
               Do i = 1,n_cg
                 dqbydt(cgi(i),k)   = dqbydt_cg(i,k)
                 dqclbydt(cgi(i),k) = dqclbydt_cg(i,k)
                 dqcfbydt(cgi(i),k) = dqcfbydt_cg(i,k)
               End Do
             End Do

           Else    ! Requires conversion increment to mixing ratio not q

             If (l_q_interact) Then       ! PC2

               Do k = 1,nlev
                 Do i = 1,n_cg
                   dqtt  = dqbydt_cg(i,k)+dqclbydt_cg(i,k)+dqcfbydt_cg(i,k)
                   dqt   = dqtt*timestep
                   denom = 1.0/((1.0+qt(cgi(i),k))*(1.0+qt(cgi(i),k)+dqt))
                   dqbydt(cgi(i),k)   = denom *                           &
                        (dqbydt_cg(i,k)*(1.0+qt(cgi(i),k))-q_cg(i,k)*dqtt)
                   dqclbydt(cgi(i),k) = denom *                           &
                      (dqclbydt_cg(i,k)*(1.0+qt(cgi(i),k))-qcl_cg(i,k)*dqtt)
                   dqcfbydt(cgi(i),k) = denom *                           &
                      (dqcfbydt_cg(i,k)*(1.0+qt(cgi(i),k))-qcf_cg(i,k)*dqtt)
                 End Do
               End Do

             Else                         ! Not PC2
                                          ! No qcl and qcf increments
               Do k = 1,nlev
                 Do i = 1,n_cg
                   dqt   = dqbydt_cg(i,k)*timestep
                   denom = 1.0/((1.0+qt(cgi(i),k))*(1.0+qt(cgi(i),k)+dqt))
                   dqbydt(cgi(i),k)   = dqbydt_cg(i,k)*denom
                   dqclbydt(cgi(i),k) = 0.0
                   dqcfbydt(cgi(i),k) = 0.0
                 End Do
               End Do
             End If                       ! End test on PC2

           End If     ! Test on l_mixing_ratio

         End If        ! Test on scheme (iconv_congestus)

          If (flg_up_flx) then
            Do k = 1,nlev
              Do i = 1,n_cg
                up_flux(cgi(i),k)        = up_flux_cg(i,k)
              End Do
            End Do
          End If

      If (flg_mf_congest) then
        Do k = 1,nlev
          Do i = 1,n_cg
            mf_congest(cgi(i),k)        = up_flux_cg(i,k)
          End Do
        End Do
      End If
      If (flg_dt_congest) then
        Do k = 1,nlev
          Do i = 1,n_cg
            dt_congest(cgi(i),k)        = dthbydt_cg(i,k)
          End Do
        End Do
      End If
      If (flg_dq_congest) then
        Do k = 1,nlev
          Do i = 1,n_cg
            dq_congest(cgi(i),k)        = dqbydt_cg(i,k)
          End Do
        End Do
      End If
      If (L_mom) then
       If (flg_du_congest) then
        Do k = 1,nlev+1
          Do i = 1,n_cg
            du_congest(cgi(i),k)        = dubydt_cg(i,k)
          End Do
        End Do
       End If
       If (flg_dv_congest) then
        Do k = 1,nlev+1
          Do i = 1,n_cg
            dv_congest(cgi(i),k)        = dvbydt_cg(i,k)
          End Do
        End Do
       End If
      End If
          If (L_tracer) then
            Do ktra = 1,ntra
              Do k = 1,nlev
                Do i = 1,n_cg
                  dtrabydt(cgi(i),k,ktra)  = dtrabydt_cg(i,k,ktra)
                End Do
              End Do
            End Do
          End If

         End If    ! number of congestus points > 0

       End If      ! iconv_congestus > 0

!-----------------------------------------------------------------------
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
! 4.0 MID-LEVEL Convection
! 4.1 Set lowest level that mid level convection can trigger
!-----------------------------------------------------------------------

      If (iconv_mid >  0) then
        idp=0
        idp2=0
        Do i = 1,npnts
          If (.not. cumulus_bl(i)) then
            midtrig(i) = ntml(i) + 1
            If (ntml(i)  ==  nbl) then
              midtrig(i) = ntml(i)
            End If
          Else
!
!  Cumulus points ie deep or shallow convection has occurred
!
            midtrig(i) = ntpar(i) + 2

!
! NTPAR has a maximum value which can be less than the top level for
! deep. Deep convection may terminate above this. Need to prevent mid
! level convection occurring at the same levels as deep.
!
            If(.not.l_shallow_bl(i).and..not.L_congestus(i)) then
              idp=idp+1
              If (kterm_dp(idp) >  ntpar(i)+2) then
                 midtrig(i) = kterm_dp(idp)
              End If
            End If  ! deep points

! if congestus calls deep scheme can terminate above ntpar
            If(L_congestus(i)) then
              idp2=idp2+1
              If (kterm_cg(idp2) >  ntpar(i)+2) then
                 midtrig(i) = kterm_cg(idp2)
              End If
            End If

          End If
        End Do

!-----------------------------------------------------------------------
! 4.2 Copy all input arrays to arrays ending in _md for passing to
!     mid-level scheme
!-----------------------------------------------------------------------

        Do k = 0,nlev
          Do i = 1,npnts
            p_layer_centres_md(i,k)        = p_layer_centres(i,k)
            p_layer_boundaries_md(i,k)     = p_layer_boundaries(i,k)
            exner_layer_centres_md(i,k)    = exner_layer_centres(i,k)
            exner_layer_boundaries_md(i,k) = exner_layer_boundaries(i,k)
            r_theta_md(i,k)                = r_theta(i,k)
          End Do
        End Do
        Do k = 1,nlev
          Do i = 1,npnts
            u_md(i,k)           = u(i,k)
            v_md(i,k)           = v(i,k)
            th_md(i,k)          = th(i,k)
            cf_liquid_md(i,k)   = cf_liquid(i,k)
            cf_frozen_md(i,k)   = cf_frozen(i,k)
            bulk_cf_md(i,k)     = bulk_cf(i,k)
            z_theta_md(i,k)     = z_theta(i,k)
            z_rho_md(i,k)       = z_rho(i,k)
            r_rho_md(i,k)       = r_rho(i,k)
          End Do
        End Do
!
! moisture  - input depends on scheme
!

        If (iconv_mid == 1) Then     ! G-R scheme
!
! G-R requires input of specific humidity
!
          If (l_mixing_ratio) Then

            Do k = 1,nlev
              Do j = 1,npnts
                q_md(j,k)   = q(j,k)  /(1.0+qt(j,k))
                qse_md(j,k) = qse(j,k)/(1.0+qse(j,k))
                qcl_md(j,k) = qcl(j,k)/(1.0+qt(j,k))
                qcf_md(j,k) = qcf(j,k)/(1.0+qt(j,k))
              End Do
            End do

          Else   ! input is specific humidity therefore no problems

            Do k = 1,nlev
              Do j = 1,npnts
                q_md(j,k)   = q(j,k)
                qse_md(j,k) = qse(j,k)
                qcl_md(j,k) = qcl(j,k)
                qcf_md(j,k) = qcf(j,k)
              End Do
            End do

          End If      ! Test on l_mixing_ratio

        Else         ! future turbulence scheme requires mixing ratio
!
! turbulence scheme requires mixing ratios as input
!
          If (l_mixing_ratio) Then

            Do k = 1,nlev
              Do j = 1,npnts
                q_md(j,k)   = q(j,k)
                qse_md(j,k) = qse(j,k)
                qcl_md(j,k) = qcl(j,k)
                qcf_md(j,k) = qcf(j,k)
              End Do
            End do

          Else       ! Require conversion

            Do k = 1,nlev
              Do j = 1,npnts

                q_md(j,k)     = q(j,k)/(1.0 -qt(j,k))

! check for negative values
                If (q_md(j,k) <  0.0) then
                  write(6,*) 'Problem q_mix -ve ',q_md(j,k),' j,k',j,k
                  q_md(j,k) = 0.0
                End If
                If (qse(j,k) <  1.0) then
                  qse_md(j,k) = qse(j,k)/(1.-qse(j,k))
                Else
                  write(6,*) 'Problem qsat 1 ',qse(j,k),' j,k',j,k   &
                             ,p_layer_centres(j,k),th(j,k),q(j,k)
                  qse_md(j,k) = 10000.    ! very big value
                End If

! Not sure whether qcl and qcf will be used or altered and don't
! know whether additional tests will be required if converting.

                qcl_md(j,k) = qcl(j,k)/(1.0-qt(j,k))
                qcf_md(j,k) = qcf(j,k)/(1.0-qt(j,k))
              End Do
            End do

          End If     ! Test on l_mixing_ratio

        End If        ! test on scheme (iconv_mid)

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
! 4.3 Call mid-level convection code
!-----------------------------------------------------------------------

        If (iconv_mid == 1) then ! Gregory-Rowntree mid level convection

! DEPENDS ON: mid_conv
          Call MID_CONV(                                                &
                      !IN
     &                nbl,nlev,ntra,n_cca_lev,npnts,trlev,              &
     &                bland,w_max,exner_layer_centres_md,               &
     &                exner_layer_boundaries_md,flg_up_flx,             &
     &                flg_entr_up, flg_detr_up, flg_dwn_flx,            &
     &                flg_entr_dwn, flg_detr_dwn, L_calc_dxek,          &
     &                L_q_interact, L_tracer, midtrig, ntml, ntpar,     &
     &                freeze_lev,                                       &
     &                pstar,p_layer_centres_md,                         &
     &                p_layer_boundaries_md,r_theta_md,r_rho_md,        &
     &                z_theta,z_rho,rho,rho_theta,                      &
     &                r2rho_th,dr_across_th,                            &
     &                q_md,q1_sd,t1_sd,th_md,timestep,                  &
     &                u_md,v_md,recip_pstar,qse,                        &
     &                md_on, mdet_md_on, md_ent_on,                     &
     &                md_sdet_on, md_new_termc,                         &
                      !INOUT
     &                bulk_cf_md,cf_frozen_md,cf_liquid_md,             &
     &                qcf_md,qcl_md,tracer_md,                          &
                      !OUT
     &                cape_out_md,cclwp_md,ccw_md,cca_md,               &
     &                dbcfbydt_md,dcffbydt_md,dcflbydt_md,              &
     &                dqbydt_md,dqcfbydt_md,dqclbydt_md,                &
     &                dthbydt_md,dubydt_md,dvbydt_md,dtrabydt_md,       &
     &                detrain_up_md,detrain_dwn_md,entrain_up_md,       &
     &                entrain_dwn_md,                                   &
     &                iccb_md,icct_md,                                  &
     &                lcca_md,lcclwp_md,lcbase_md,lctop_md,             &
     &                rain_md,snow_md,rain_3d_md,snow_3d_md,            &
     &                up_flux_md,                                       &
     &                dwn_flux_md,tcw_md,l_mid_all,cca_2d_md,           &
     &                rbuoy_p_out_md,the_out_md,thp_out_md,             &
     &                qe_out_md,qp_out_md,                              &
     &                uw_mid,vw_mid                                    &
     &                )

         Else        ! turbulence based alternative
           write(6,*) ' New mid-level scheme not avilable yet '
!          call mid_turb_conv( ) ?
         End If

!-----------------------------------------------------------------------
! 4.4 Write data from mid-level convection to full arrays
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

          !===================================================================
          ! NOTE: At this point if Lcv_ccrad = T, then cca_2d_md is ONLY equal
          !       to that from mid-level cloud on a given grid point.  The
          !       original code I.E. Lcv_ccrad = F means that cca_2d_md will
          !       include that from sh/dp aswell.
          !===================================================================

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
              ccw(i,k)            = ccw(i,k)      + ccw_md(i,k)
            End Do
          End Do
        End If

        Do k = 1,nlev
          Do i = 1,npnts
            ccw(i,k)              = ccw(i,k)     + ccw_md(i,k)
            cca(i,k)              = cca(i,k)     + cca_md(i,k)
            dthbydt(i,k)          = dthbydt(i,k) + dthbydt_md(i,k)
            dcflbydt(i,k)         = dcflbydt(i,k) + dcflbydt_md(i,k)
            dcffbydt(i,k)         = dcffbydt(i,k) + dcffbydt_md(i,k)
            dbcfbydt(i,k)         = dbcfbydt(i,k) + dbcfbydt_md(i,k)



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
! moisture  - output depends on scheme
!
        If (iconv_mid == 1) Then     ! G-R scheme
!
! Output in specific humidity
!
          If (l_mixing_ratio) Then ! requires conversion

            If (l_q_interact) Then    ! PC2

              Do k = 1,nlev
                Do i = 1,npnts
                 dqtt  = dqbydt_md(i,k)+dqclbydt_md(i,k)+dqcfbydt_md(i,k)
                 dqt   = dqtt*timestep
                 denom = 1.0/((1.0-qt(i,k))*(1.0-qt(i,k)-dqt))
                 dqbydt(i,k)   = dqbydt(i,k)   + denom *                    &
                          ( dqbydt_md(i,k)*(1.-qt(i,k))+q(i,k)*dqtt )
                 dqclbydt(i,k) = dqclbydt(i,k) + denom *                    &
                        ( dqclbydt_md(i,k)*(1.-qt(i,k))+qcl(i,k)*dqtt )
                 dqcfbydt(i,k) = dqcfbydt(i,k) + denom *                    &
                        ( dqcfbydt_md(i,k)*(1.-qt(i,k))+qcf(i,k)*dqtt )
                End Do
              End do

            Else                       ! Not PC2
                                       ! No qcl and qcf increments anyway
              Do k = 1,nlev
                Do i = 1,npnts
                  dqt = dqbydt_md(i,k)*timestep
                  denom = 1.0/((1.0-qt(i,k))*(1.0-qt(i,k)-dqt))
                  dqbydt(i,k)    =  dqbydt(i,k) + dqbydt_md(i,k)*denom
                  dqclbydt(i,k)  = 0.0
                  dqcfbydt(i,k)  = 0.0
                End Do
              End do

            End If                     ! End test on PC2

          Else   ! output is specific humidity therefore no problems

             Do k = 1,nlev
               Do i = 1,npnts
                 dqbydt(i,k)  = dqbydt(i,k) + dqbydt_md(i,k)
                 dqclbydt(i,k) = dqclbydt(i,k) + dqclbydt_md(i,k)
                 dqcfbydt(i,k) = dqcfbydt(i,k) + dqcfbydt_md(i,k)
               End Do
             End do

          End If     ! Test on scheme (iconv_mid)

        Else         ! turbulence scheme requires mixing ratio
!
! turbulence scheme outputs mixing ratios
!
          If (l_mixing_ratio) Then

            Do k = 1,nlev
              Do i = 1,npnts
                dqbydt(i,k)   = dqbydt(i,k) + dqbydt_md(i,k)
                dqclbydt(i,k) = dqclbydt(i,k) + dqclbydt_md(i,k)
                dqcfbydt(i,k) = dqcfbydt(i,k) + dqcfbydt_md(i,k)
              End Do
            End do

          Else    ! Requires conversion increment to mixing ratio not q

            If (l_q_interact) Then    ! PC2

              Do k = 1,nlev
                Do i = 1,npnts
                  dqtt  = dqbydt_md(i,k)+dqclbydt_md(i,k)+dqcfbydt_md(i,k)
                  dqt   = dqtt*timestep
                  denom = 1.0/((1.0+qt(i,k))*(1.0+qt(i,k)+dqt))
                  dqbydt(i,k)   = dqbydt(i,k)   + denom *                &
                         (dqbydt_md(i,k)*(1.0+qt(i,k))-q_md(i,k)*dqtt)
                  dqclbydt(i,k) = dqclbydt(i,k) + denom *                &
                       (dqclbydt_md(i,k)*(1.0+qt(i,k))-qcl(i,k)*dqtt)
                  dqcfbydt(i,k) = dqcfbydt(i,k) + denom *                &
                       (dqcfbydt_md(i,k)*(1.0+qt(i,k))-qcf(i,k)*dqtt)
                End Do
              End Do

            Else                      ! Not PC2
                                       ! No qcl and qcf increments
              Do k = 1,nlev
                Do i = 1,npnts
                  dqt   = dqbydt_md(i,k)*timestep
                  denom = 1.0/((1.0+qt(i,k))*(1.0+qt(i,k)+dqt))
                  dqbydt(i,k)   = dqbydt(i,k) + dqbydt_md(i,k)*denom
                  dqclbydt(i,k) = 0.0
                  dqcfbydt(i,k) = 0.0
                End Do
              End Do

            End If                    ! End test on PC2

          End If  ! Test on l_mixing_ratio

        End If        ! test on scheme (iconv_mid)

        If (flg_up_flx) then
          Do k = 1,nlev
            Do i = 1,npnts
               up_flux(i,k)    = up_flux(i,k) + up_flux_md(i,k)
            End Do
          End Do
        End If
        If (flg_dwn_flx) then
          Do k = 1,nlev
            Do i = 1,npnts
              dwn_flux(i,k)    = dwn_flux(i,k) + dwn_flux_md(i,k)
            End Do
          End Do
        End If
        If (flg_entr_up) then
          Do k = 1,nlev
            Do i = 1,npnts
              entrain_up(i,k)  = entrain_up(i,k)+ entrain_up_md(i,k)
            End Do
          End Do
        End If
        If (flg_detr_up) then
          Do k = 1,nlev
            Do i = 1,npnts
              detrain_up(i,k)  = detrain_up(i,k) + detrain_up_md(i,k)
            End Do
          End Do
        End If
        If (flg_entr_dwn) then
          Do k = 1,nlev
            Do i = 1,npnts
              entrain_dwn(i,k) = entrain_dwn(i,k) + entrain_dwn_md(i,k)
            End Do
          End Do
        End If
        If (Flg_detr_dwn) then
          Do k = 1,nlev
            Do i = 1,npnts
              detrain_dwn(i,k) = detrain_dwn(i,k) + detrain_dwn_md(i,k)
            End Do
          End Do
        End If
        If (flg_mf_midlev) then
          Do k = 1,nlev
            Do i = 1,npnts
              mf_midlev(i,k)        = up_flux_md(i,k)
            End Do
          End Do
        End If
        If (flg_dt_midlev) then
          Do k = 1,nlev
            Do i = 1,npnts
              dt_midlev(i,k)        = dthbydt_md(i,k)
            End Do
          End Do
        End If
        If (flg_dq_midlev) then
          Do k = 1,nlev
            Do i = 1,npnts
              dq_midlev(i,k)        = dqbydt_md(i,k)
            End Do
          End Do
        End If
        If (L_mom) then
         If (flg_du_midlev) then
           Do k = 1,nlev+1
             Do i = 1,npnts
               du_midlev(i,k)        = dubydt_md(i,k)
             End Do
           End Do
         End If
         If (flg_dv_midlev) then
           Do k = 1,nlev+1
             Do i = 1,npnts
               dv_midlev(i,k)        = dvbydt_md(i,k)
             End Do
           End Do
         End If
        End If
        If (L_tracer) then
          Do ktra = 1,ntra
            Do k = 1,nlev
              Do i = 1,npnts
                dtrabydt(i,k,ktra) = dtrabydt(i,k,ktra)                 &
     &                                          + dtrabydt_md(i,k,ktra)
              End Do
            End Do
          End Do
        End If

!-----------------------------------------------------------------------------
! 4.5 CCRad - Apply CCA_md/CCW_md to full field .....
!-----------------------------------------------------------------------------

        If (lcv_ccrad) Then

            Do k=1,nlev
              Do i=1, npnts
                cca(i,k) = cca(i,k) + cca_md(i,k)  
                ccw(i,k) = ccw(i,k) + ccw_md(i,k)
              End do      ! i (npnts)
            End do      ! k (nlev)

        End If      ! lcv_ccrad

      End If     ! test on iconv_mid

!-----------------------------------------------------------------------
! 5.0 Work out which points convection has occurred at.
!-----------------------------------------------------------------------

      nconv_all=0
      Do i = 1,npnts
        If (cumulus_bl(i).or.l_mid_all(i)) then
          nconv_all = nconv_all + 1
          index1(nconv_all) = i
        End If
      End Do

!-----------------------------------------------------------------------
! 6.0  Update tracer field
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
#endif
