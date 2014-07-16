#if defined(A05_0A)
!
! SUBROUTINE GLUE_CONV--------------------------------------------------------
!
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!
!+ Gather-scatter routine for deep and shallow convection points
!

SUBROUTINE glue_conv(np_field, npnts, nlev, nbl, th, q, qcl, qcf, cf_liquid   &
         , cf_frozen, bulk_cf, pstar, bland, u, v, tracer, dthbydt, dqbydt    &
         , dqclbydt, dqcfbydt, dcflbydt, dcffbydt, dbcfbydt, dubydt, dvbydt   &
         , rain, snow, rain_3d, snow_3d, cca, iccb, icct, cclwp, ccw          &
         , lcbase, lctop, lcca, freeze_lev, l_mid_all, kterm_deep             &
         , precip_deep, precip_shall, precip_mid, precip_cong                 &
         , wstar_dn, wstar_up, mb1, mb2, kterm_congest, n_cumulus, uw0, vw0   &
         , w_max, zlcl, zlcl_uv, ztop_uv, lcclwp, cape_out, n_dp, n_cg, n_sh  &
         , r_rho, r_theta, rho, rho_theta                                     &
         , exner_layer_boundaries, exner_layer_centres, p_layer_boundaries    &
         , p_layer_centres, z_theta, z_rho, timestep, t1_sd, q1_sd, ntml      &
         , ntpar, nbdsc, ntdsc, l_shallow_bl, l_pc2_diag_sh_pts               &
         , l_congestus, cumulus_bl, wstar                                     &
         , wthvs, delthvu_bl, ql_ad, ftl, fqt, l_tracer, ntra, trlev          &
         , n_cca_lev, l_mixing_ratio, l_calc_dxek, l_q_interact               &
         , up_flux_half, flg_up_flx_half, up_flux, flg_up_flx                 &
         , dwn_flux,     flg_dwn_flx                                          &
         , entrain_up,  flg_entr_up,  detrain_up,  flg_detr_up                &
         , entrain_dwn, flg_entr_dwn, detrain_dwn, flg_detr_dwn               &
         , uw_deep,     flg_uw_deep,  vw_deep,     flg_vw_deep                &
         , uw_shall,    flg_uw_shall, vw_shall,    flg_vw_shall               &
         , flg_wqt_flux_sh, flg_wthetal_flux_sh, flg_wthetav_flux_sh          &
         , flg_wql_flux_sh                                                    &
         , flg_mf_deep, flg_mf_congest, flg_mf_shall, flg_mf_midlev           &
         , flg_dt_deep, flg_dt_congest, flg_dt_shall, flg_dt_midlev           &
         , flg_dq_deep, flg_dq_congest, flg_dq_shall, flg_dq_midlev           &
         , flg_du_deep, flg_du_congest, flg_du_shall, flg_du_midlev           &
         , flg_dv_deep, flg_dv_congest, flg_dv_shall, flg_dv_midlev           &
         , wqt_flux_sh, wthetal_flux_sh, wthetav_flux_sh, wql_flux_sh         &
         , mf_deep, mf_congest, mf_shall, mf_midlev, dt_deep, dt_congest      &
         , dt_shall, dt_midlev, dq_deep, dq_congest, dq_shall, dq_midlev      &
         , du_deep, du_congest, du_shall, du_midlev, dv_deep, dv_congest      &
         , dv_shall, dv_midlev                                                &
         )
!
! Description:
!   Dummy stub version Glue_Conv to be called when Convection Scheme routines
!   are disabled.
!
!   Called by Load_bal_conv/ni_conv_ctl
!
! Current code owner: R A Stratton
!
! Code Description:
!  Language: FORTRAN 95
!  This code is written to UMDP3 v8 programming standards
!
! ----------------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------------
!
! Arguments with intent IN:
!
  INTEGER, INTENT(IN) :: np_field
                                      ! No. of points in a full field (NOTE:
                                      ! all multi-dimensional fields passed in
                                      ! MUST be dimensioned with np_field).
  INTEGER, INTENT(IN) :: npnts        ! No. of points in segment.
  INTEGER, INTENT(IN) :: nlev         ! No. of model layers.
  INTEGER, INTENT(IN) :: n_cca_lev    ! No. of convective cloud amount levels:
                                      !                     1 for 2D
                                      !                 nlevs for 3D)
  INTEGER, INTENT(IN) :: nbl          ! No. of boundary layer levels.
  INTEGER, INTENT(IN) :: ntra         ! No. of tracer fields.
  INTEGER, INTENT(IN) :: trlev        ! No. of model levels on which tracers
                                      ! are included.
  INTEGER, INTENT(IN) :: n_dp         ! Number of deep points.
  INTEGER, INTENT(IN) :: n_cg         ! Number of congestus points.
  INTEGER, INTENT(IN) :: n_sh         ! Number of shallow points.
  INTEGER, INTENT(IN) :: ntml(npnts)  ! Top level of surface mixed layer
                                      ! defined relative to theta, q grid.
  INTEGER, INTENT(IN) :: ntpar(npnts) ! Top level of initial parcel ascent in
                                      ! BL scheme defined relative to theta  
                                      ! ,q grid.

  REAL, INTENT(IN) :: timestep           ! Model timestep (s)

  REAL, INTENT(IN) :: pstar      (npnts) ! Surface pressure (Pa)
  REAL, INTENT(IN) :: wstar      (npnts) ! Convective velocity scale (m/s)
  REAL, INTENT(IN) :: wthvs      (npnts) ! Surface flux of THV (Pa m/s2)
  REAL, INTENT(IN) :: ql_ad      (npnts) ! Adiabatic liquid water content.
  REAL, INTENT(IN) :: ftl        (npnts) ! Surface sensible heat flux from BL
                                         ! (W/m2)
  REAL, INTENT(IN) :: fqt        (npnts) ! Total water flux from surface(kg/m2/s)
  REAL, INTENT(IN) :: zlcl       (npnts) ! Lifting condensation level (not a
                                         ! model level) (m)
  REAL, INTENT(IN) :: zlcl_uv    (npnts) ! Lifting condensation level defined
                                         ! for the uv grid (m)
  REAL, INTENT(IN) :: ztop_uv    (npnts) ! Top of cloud layer defined for the
                                         ! uv grid (m)
  REAL, INTENT(IN) :: delthvu_bl (npnts) ! Integral of undilute parcel
                                         ! buoyancy over convective cloud
                                         ! layer (Kelvin m)

  REAL, INTENT(IN) :: uw0     (np_field) ! U-comp of surface stress(N/m2)
  REAL, INTENT(IN) :: vw0     (np_field) ! V-comp of surface stress(N/m2)
  REAL, INTENT(IN) :: w_max   (np_field) ! Max. w in column

  REAL, INTENT(IN) :: t1_sd   (np_field) ! Standard deviation of turbulent
                                         ! fluctuations of layer 1 temp. (K)
  REAL, INTENT(IN) :: q1_sd   (np_field) ! Standard deviation of turbulent
                                         ! fluctuations of layer 1 q (kg/kg)

  REAL, INTENT(IN) :: th       (np_field,   nlev) ! Model potential temp. (K)
  REAL, INTENT(IN) :: q        (np_field,   nlev) ! Model mixing ratio (kg/kg)
  REAL, INTENT(IN) :: u        (np_field,   nlev) ! Model U field (m/s)
  REAL, INTENT(IN) :: v        (np_field,   nlev) ! Model V field (m/s)
  REAL, INTENT(IN) :: r_rho    (np_field,   nlev) ! Radius rho lev (m)
  REAL, INTENT(IN) :: r_theta  (np_field, 0:nlev) ! Radius th lev (m)
  REAL, INTENT(IN) :: rho      (np_field,   nlev) ! Density, rho lvl.(kg/m3)
  REAL, INTENT(IN) :: rho_theta(np_field,   nlev) ! Density, theta lvl.(kg/m3)
  REAL, INTENT(IN) :: z_theta  (np_field,   nlev) ! Height of theta levels
                                                  ! above surface (m)
  REAL, INTENT(IN) :: z_rho    (np_field,   nlev) ! Height of rho   levels
                                                  ! above surface (m)

  REAL, INTENT(IN) :: exner_layer_centres    (np_field, 0:nlev)
                             ! Exner
  REAL, INTENT(IN) :: p_layer_centres        (np_field, 0:nlev)
                             ! Pressure(Pa)
  REAL, INTENT(IN) :: exner_layer_boundaries (np_field, 0:nlev)
                             ! Exner at half level above exner_layer_centres
  REAL, INTENT(IN) :: p_layer_boundaries     (np_field, 0:nlev)
                             ! Pressure at half level above
                             ! p_layer_centres (Pa)

  LOGICAL, INTENT(IN) :: l_tracer     ! Switch for inclusion of tracers
  LOGICAL, INTENT(IN) :: flg_entr_up  ! STASH flag for updraught entrainment
  LOGICAL, INTENT(IN) :: flg_detr_up  ! STASH flag for updraught detrainment
  LOGICAL, INTENT(IN) :: flg_dwn_flx  ! STASH flag for downdraught mass flux
  LOGICAL, INTENT(IN) :: flg_entr_dwn ! STASH flag for downdraught entrainment
  LOGICAL, INTENT(IN) :: flg_detr_dwn ! STASH flag for downdraught detrainment
  LOGICAL, INTENT(IN) :: flg_uw_deep  ! STASH flag for x-comp deep stress
  LOGICAL, INTENT(IN) :: flg_vw_deep  ! STASH flag for y-comp deep stress
  LOGICAL, INTENT(IN) :: flg_uw_shall ! STASH flag for x-comp shallow stress
  LOGICAL, INTENT(IN) :: flg_vw_shall ! STASH flag for y-comp shallow stress
  LOGICAL, INTENT(IN) :: flg_up_flx   ! STASH flag for updraught mass flux
  LOGICAL, INTENT(IN) :: flg_up_flx_half
                                      ! STASH flag for updraught mass flux on
                                      ! half levels

  LOGICAL, INTENT(IN) :: l_mixing_ratio 
                                      ! .TRUE.  for moisture inputs as
                                      !         mixing ratio
                                      ! .FALSE. for moisture inputs as
                                      !         specific humidity
  LOGICAL, INTENT(IN) :: l_calc_dxek  ! Switch for calculation of condensate
                                      ! increment.
  LOGICAL, INTENT(IN) :: l_q_interact ! Switch allows overwriting parcel
                                      ! variables when calculating condensate
                                      ! increment.


  LOGICAL, INTENT(IN) :: bland        (npnts) ! Land/sea mask
  LOGICAL, INTENT(IN) :: l_shallow_bl (npnts) ! Shallow cumulus indicator
  LOGICAL, INTENT(IN) :: l_pc2_diag_sh_pts(npnts) ! Shallow cumulus PC2
  LOGICAL, INTENT(IN) :: l_congestus  (npnts) ! Congestus cumulus
  LOGICAL, INTENT(IN) :: cumulus_bl   (npnts) ! Cumulus indicator

!
! Flags for flux diagnostics - only possible from turbulence based schemes
!

  LOGICAL, INTENT(IN) :: flg_wqt_flux_sh      ! STASH flag for wqt_flux_sh
  LOGICAL, INTENT(IN) :: flg_wql_flux_sh      ! STASH flag for wql_flux_sh
  LOGICAL, INTENT(IN) :: flg_wthetal_flux_sh  ! STASH flag for wthetal_flux_sh
  LOGICAL, INTENT(IN) :: flg_wthetav_flux_sh  ! STASH flag for wthetav_flux_sh

!
! Flags for convection scheme (5A scheme)
!
  INTEGER, INTENT(IN) :: iccb       (npnts)   ! Conv. cld base level (m)
  INTEGER, INTENT(IN) :: icct       (npnts)   ! Conv. cld top level (m)
  INTEGER, INTENT(IN) :: lcbase     (npnts)   ! Lowest conv. cld base level(m)
  INTEGER, INTENT(IN) :: lctop      (npnts)   ! Lowest conv. cld  top level(m)
  INTEGER, INTENT(IN) :: freeze_lev (npnts)   ! Index for freezing level
  INTEGER, INTENT(IN) :: kterm_deep (npnts)   ! Index for deep conv.
                                              ! termination


  REAL, INTENT(IN) :: tracer    (np_field, trlev, ntra)
                                                 ! Model tracer fields(kg/kg)
  REAL, INTENT(IN) :: qcl       (np_field, nlev) ! Liq condensate mix ratio
                                                 ! (kg/kg)
  REAL, INTENT(IN) :: qcf       (np_field, nlev) ! Ice condensate mix ratio
                                                 ! (kg/kg)
  REAL, INTENT(IN) :: cf_liquid (np_field, nlev) ! Liq water cld vol.(frac.)
  REAL, INTENT(IN) :: cf_frozen (np_field, nlev) ! Ice water cld vol.(frac?)
  REAL, INTENT(IN) :: bulk_cf   (np_field, nlev) ! Bulk total cloud volume
  REAL, INTENT(IN) :: dqclbydt  (np_field, nlev) ! Incs to liq condensate
                                                 ! due to convection (kg/kg/s)
  REAL, INTENT(IN) :: dqcfbydt  (np_field, nlev) ! Incs to ice condensate
                                                 ! due to convection (kg/kg/s)
  REAL, INTENT(IN) :: dcflbydt  (np_field, nlev) ! Incs to liq cloud volume
                                                 ! due to convection (/s)
  REAL, INTENT(IN) :: dcffbydt  (np_field, nlev) ! Incs to ice cloud volume
                                                 ! due to convection (/s)
  REAL, INTENT(IN) :: dbcfbydt  (np_field, nlev) ! Incs to total cloud volume
                                                 ! due to convection(/s)
  REAL, INTENT(IN) :: dthbydt   (np_field, nlev) ! Incs to potential temp.
                                                 ! due to convection (K/s)
  REAL, INTENT(IN) :: dqbydt    (np_field, nlev) ! Incs to q due to
                                                 ! convection (kg/kg/s)
  REAL, INTENT(IN) :: dubydt    (np_field, nlev+1)
                                                 ! Incs to U due to CMT (m/s2)
  REAL, INTENT(IN) :: dvbydt    (np_field, nlev+1)
                                                 ! Increments to V due to CMT
                                                 ! (m/s2)

  REAL, INTENT(IN) :: rain         (npnts) ! Srf. conv. rainfall (kg/m2/s)
  REAL, INTENT(IN) :: snow         (npnts) ! Srf. conv. snowfall (kg/m2/s)
  REAL, INTENT(IN) :: rain_3d      (npnts) ! 3D   conv. rainfall (kg/m2/s)
  REAL, INTENT(IN) :: snow_3d      (npnts) ! 3D   conv. snowfall (kg/m2/s)
  REAL, INTENT(IN) :: cclwp        (npnts) ! Condensed water path (k/m2)
  REAL, INTENT(IN) :: lcca         (npnts) ! Lowest conv. cloud amt. (%)
  REAL, INTENT(IN) :: precip_deep  (npnts) ! Deep    precipitation (kg/m2/s)
  REAL, INTENT(IN) :: precip_shall (npnts) ! Shallow precipitation (kg/m2/s)
  REAL, INTENT(IN) :: precip_mid   (npnts) ! Mid-lvl precipitation (kg/m2/s)
  REAL, INTENT(IN) :: precip_cong  (npnts) ! Congestus precip.     (kg/m2/s)
  REAL, INTENT(IN) :: lcclwp       (npnts) ! Condensed water path for lowest
                                           ! conv. cld. level (kg/m2)
  REAL, INTENT(IN) :: cape_out     (npnts) ! Saved convective available
                                           ! potential energy for diagnostic
                                           ! output (Jkg-1)


  REAL, INTENT(IN) :: cca         (np_field, n_cca_lev)
                                                   ! Conv. cloud amount (%)
  REAL, INTENT(IN) :: ccw         (np_field, nlev) ! Conv. cloud liquid water
                                                   ! on model levels (g/kg)
  REAL, INTENT(IN) :: up_flux     (np_field, nlev) ! Updraught mass flux
                                                   ! (Pa/s)
  REAL, INTENT(IN) :: up_flux_half(np_field, nlev) ! Updraught mass flux on
                                                   ! half levels
                                                   ! (Pa/s)
  REAL, INTENT(IN) :: dwn_flux    (np_field, nlev) ! Downdraught mass flux
                                                   ! (Pa/s)
  REAL, INTENT(IN) :: entrain_up  (np_field, nlev) ! Fractional entrainment
                                                   ! rate into updraughts
                                                   ! (Pa/s)
  REAL, INTENT(IN) :: detrain_up  (np_field, nlev) ! Fractional detrainment
                                                   ! rate into updraughts
                                                   ! (Pa/s)
  REAL, INTENT(IN) :: entrain_dwn (np_field, nlev) ! Fractional entrainment
                                                   ! rate into Downdraughts
                                                   ! (Pa/s)
  REAL, INTENT(IN) :: detrain_dwn (np_field, nlev) ! Fractional detrainment
                                                   ! rate into downdraughts
                                                   ! (Pa/s)
  REAL, INTENT(IN) :: uw_deep     (np_field, nlev) ! X-comp. of stress from
                                                   ! deep convection(kg/m/s2)
  REAL, INTENT(IN) :: vw_deep     (np_field, nlev) ! Y-comp. of stress from
                                                   ! deep convection(kg/m/s2)
  REAL, INTENT(IN) :: uw_shall    (np_field, nlev) ! X-comp. of stress from
                                                   ! shallow convection
                                                   !(kg/m/s2)
  REAL, INTENT(IN) :: vw_shall    (np_field, nlev) ! Y-comp. of stress from
                                                   ! shallow convection
                                                   !(kg/m/s2)

!
! Flux diagnostics - only avialable from turbulence based schemes
!

  REAL, INTENT(IN) :: mb1             (npnts)
  REAL, INTENT(IN) :: mb2             (npnts)
  REAL, INTENT(IN) :: wstar_dn        (npnts)
  REAL, INTENT(IN) :: wstar_up        (npnts)
  REAL, INTENT(IN) :: wqt_flux_sh     (np_field, nlev) ! w'qt' flux
  REAL, INTENT(IN) :: wql_flux_sh     (np_field, nlev) ! w'ql' flux
  REAL, INTENT(IN) :: wthetal_flux_sh (np_field, nlev) ! w'thetal' flux
  REAL, INTENT(IN) :: wthetav_flux_sh (np_field, nlev) ! w'thetav' flux
  REAL, INTENT(IN) :: mf_shall        (np_field, nlev)
  REAL, INTENT(IN) :: dt_shall        (np_field, nlev)
  REAL, INTENT(IN) :: dq_shall        (np_field, nlev)
  REAL, INTENT(IN) :: du_shall        (np_field, nlev)
  REAL, INTENT(IN) :: dv_shall        (np_field, nlev)
  REAL, INTENT(IN) :: dt_congest      (np_field, nlev)
  REAL, INTENT(IN) :: dq_congest      (np_field, nlev)
  REAL, INTENT(IN) :: du_congest      (np_field, nlev)
  REAL, INTENT(IN) :: dv_congest      (np_field, nlev)
  REAL, INTENT(IN) :: mf_congest      (np_field, nlev)
  REAL, INTENT(IN) :: mf_midlev       (np_field, nlev)
  REAL, INTENT(IN) :: dt_midlev       (np_field, nlev)
  REAL, INTENT(IN) :: dq_midlev       (np_field, nlev)
  REAL, INTENT(IN) :: du_midlev       (np_field, nlev)
  REAL, INTENT(IN) :: dv_midlev       (np_field, nlev)
  REAL, INTENT(IN) :: mf_deep         (np_field, nlev)
  REAL, INTENT(IN) :: dt_deep         (np_field, nlev)
  REAL, INTENT(IN) :: dq_deep         (np_field, nlev)
  REAL, INTENT(IN) :: du_deep         (np_field, nlev)
  REAL, INTENT(IN) :: dv_deep         (np_field, nlev)

  INTEGER, INTENT(IN) :: kterm_congest(npnts)  ! termination level for
                                               ! congestus
  LOGICAL, INTENT(IN) :: l_mid_all(npnts)      ! .TRUE. if mid level
                                               ! convection

!
! flag not required for 4A scheme used by 5A
!

  LOGICAL, INTENT(IN) :: flg_mf_deep
  LOGICAL, INTENT(IN) :: flg_mf_congest
  LOGICAL, INTENT(IN) :: flg_mf_shall
  LOGICAL, INTENT(IN) :: flg_mf_midlev
  LOGICAL, INTENT(IN) :: flg_dt_deep
  LOGICAL, INTENT(IN) :: flg_dt_congest
  LOGICAL, INTENT(IN) :: flg_dt_shall
  LOGICAL, INTENT(IN) :: flg_dt_midlev
  LOGICAL, INTENT(IN) :: flg_dq_deep
  LOGICAL, INTENT(IN) :: flg_dq_congest
  LOGICAL, INTENT(IN) :: flg_dq_shall
  LOGICAL, INTENT(IN) :: flg_dq_midlev
  LOGICAL, INTENT(IN) :: flg_du_deep
  LOGICAL, INTENT(IN) :: flg_du_congest
  LOGICAL, INTENT(IN) :: flg_du_shall
  LOGICAL, INTENT(IN) :: flg_du_midlev
  LOGICAL, INTENT(IN) :: flg_dv_deep
  LOGICAL, INTENT(IN) :: flg_dv_congest
  LOGICAL, INTENT(IN) :: flg_dv_shall
  LOGICAL, INTENT(IN) :: flg_dv_midlev

!
! Redundant arguments
! -------------------

  INTEGER, INTENT(IN) :: n_cumulus
  REAL,    INTENT(IN) :: ntdsc
  REAL,    INTENT(IN) :: nbdsc





! Local Variables
!-----------------

  CHARACTER(Len=9 )   :: Message
  CHARACTER(Len=52)   :: RoutineName
  INTEGER             :: ErrorStat        ! Return code:
                                          !   0 = Normal exit
                                          ! +ve = Fatal Error
                                          ! -ve = Warning






!-----------------------------------------------------------------------------
! Code Statements
!-----------------------------------------------------------------------------
        
  ErrorStat   = 1
  RoutineName = 'GLUE_CONV'
  Message     = 'Convection Scheme Routines unavailable - see output.'

  WRITE (6,*) '**ERROR**: GLUE_CONV called but is unavailable.  '
  WRITE (6,*) '  Sections 5: Convection Scheme is required '


! DEPENDS ON: ereport
  CALL ereport(RoutineName, ErrorStat, Message)
        

!-----------------------------------------------------------------------------


  RETURN

!-----------------------------------------------------------------------------
END SUBROUTINE GLUE_CONV
#endif
