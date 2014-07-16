!======================== COMDECK UMSCALAR ========================
!LL
!LL  Model            Modification history:
!LL  version  date
!LL   4.0   01/06/95 Insert ice variables which are to be read in via
!LL                  UI rather than set locally.  C.Cooper.
!LL
!LL   4.1   03/04/96 Insert variables required for Gent and McWilliams
!LL                  scheme and for tapering of isopycnals  C.Roberts
!LL   4.4   16/09/97  ! Include variables for depth in quadratic Large
!LL   4.4   15/06/97 Add Free Surface solution timestep and
!LL                  weighting factor used with the delplus
!LL                  operator in TROPIC            R.Lenton
!LL   4.5   05/08/97 Removed old boundary logicals and added constants
!LL                  for Gill ocean boundary routine and fieldcodes.
!LL                  C.G. Jones
!LL   4.5   3.11.98   Set variables for Med/Hud outflow   M. Roberts
!LL   4.5    10/11/98   Add variables for biharmonic momentum and GM
!LL                     and Large scheme into common block
!LL   4.5   3.9.98   Add new variables for HADCM4 sea-ice scheme
!LL                  and remove cavitating fluid variables.
!LL                  Doug Cresswell and Jonathan Gregory
!LL 5.2  31/07/00   Add logicals to control print out of land/sea
!LL                 mask and streamfn indices by set_constants_ocean
!LL                 M J Bell & R S Hill
!LL   5.2   03/08/00 Add parameters for use in RELAX. Also moved
!LL                  MXSCAN and MSCAN here from COMOCFLW.  D.Storkey
!LL   5.3   04/10/01 Slab model now uses this include file for defining
!LL                  variables but does not use the SCALAR common block
!LL                                                          K.Williams
!LL   5.3   17/09/01 Add parameter for tracer flux diagnostics. S.Spall
!LL   5.3   20/09/01 Remove C14to12_ATM_0. S.Spall
!   5.4    Apr 2002  Use general PI rather than yet another locally
!                    defined version. Define constants as parameters.
!                    R. Hill
!LL   5.5   28/02/03 Added multiple ice category parameters: HS_MIN,
!LL                  HFRAZILMIN and ICEPUNY, and ice ridging parameters:
!LL                  GSTAR, HSTAR, CS_RDG, CF_RDG and FSNOW_RDG.
!   5.5    Feb 2003  Add parameter for initial atmospheric carbon 14
!                    concentration, PC14_ATM_0. S.Liddicoat
!   6.0    Sep 2003  Add MCPHEE_COEFF and MCPHEE_MINUSTAR for
!                     ocean-ice heat flux parameterization.  A.McLaren
!   6.1    27/07/04  Add HANEY_SALT. A. Hines
!   6.1    May 2004  Add ISTRHSM_LAT      A.McLaren
!!   6.2    Nov 2005 Add BH. M Roberts, R Hill
!   6.2    Dec 2005  Add ICE_VEL_FACTOR.  A.McLaren
!   6.2    22/05/06  Minor FCM continuation fix. P.Selwood
!   6.2    29/11/05  Add TMIN1_MAX,ATHKDFT_MIN and ATHKDFT_MAX
!                    (no longer hardwired). C. Harris
!LL

#include "c_pi.h"
      LOGICAL ROWPRT,                                                   &
                        ! Switch for printouts in ROWCALC
     &        TSUVPRT,                                                  &
                        !   "     "  printout of T,S,U,V in UMWRITOT
     &        ANCILPRT  !   "     "     "     "  ancil. fields " "
      LOGICAL L_OPRT_LSM ! print out land/sea mask
      LOGICAL L_OPRT_SF_INDX ! print out streamfn indices
      INTEGER ISTRT,ISTOP,                                              &
     &        MSCAN,                                                    &
                        ! Counter for number of scans of SOR in RELAX
     &        MXSCAN,                                                   &
                        ! Max number of scans in SOR of RELAX
     &        MIN_SCANS ! Min number of scans in SOR of RELAX
!
!         DECLARE REAL NUMBER VARIABLES
!
      REAL KAPPA0_SI,                                                   &
                             ! Vert. diffusivity at surface
     &     DKAPPA_DZ_SI,                                                &
                             ! Vert. diffusivity rate of increase
     &     KAPPA_MIN_FRAC,                                              &
                             ! Factor to reduce surf. vert. diff. to
     &     KAPPA_LOW_DEPTH   !Depth over which to reduce vert. diff.

!
      REAL AH,                                                          &
                             ! Coeff of horizontal mixing of T
     & BH,                                                              &
            ! Coeff mixing for biharmonic T
     &     DTBT,                                                        &
                             ! Free Surface barotropic timestep
     &     DTSF,                                                        &
                             ! Time step on STREAM FUNCTION
     &     DTTS,                                                        &
                             ! Time step on TRACERS
     &     DTUV,                                                        &
                             ! Time step on internal mode of U,V
     &     C2DTBT,                                                      &
                             ! 2 * time step for free surface solution
     &     C2DTSF,                                                      &
                             ! 2 * time step on STREAM FUNCTION
     &     C2DTTS,                                                      &
                             ! 2 * time step for TRACERS
     &     C2DTUV,                                                      &
                             ! 2 * time step on internal mode of U,V
     &     OMEGA,                                                       &
                             ! Rate of rotation of coordinate system.
     &     ACOR,                                                        &
                             ! Used in coriolis term weighting
     &     CRIT,                                                        &
                             ! Criterion for relaxation convergence
     &     SOR,                                                         &
                             ! Coeff of over-relaxation (1.5-1.8)
     &     FIRST_GUESS_FRAC,                                            &
                             ! factor multiplying first guess in RELAX
     &     GRAV_SI,                                                     &
                             ! Accel due to gravity SI units
     &     GRAV,                                                        &
                             ! Accel due to gravity cgs units
     &     RHO_WATER_SI,                                                &
                             ! Density 1026 kg/m**3
     &     SPECIFIC_HEAT_SI,                                            &
                             ! S.H. SI units
     &     CONV,                                                        &
                             ! SI to cgs in wind stress calcs.
     &     RADIAN,                                                      &
                             ! 57.29578 degrees
     &     RADIUS                                                       &
                             ! Radius of Earth
     &    ,WGHT_DELPLUS      ! Weighting for delplus in TROPIC
!
      REAL AM                                                           &
                             ! Coeff of horizontal mixing of U,V
     &    ,AM0                                                          &
                             ! } Coeffs for lateral viscosity
     &    ,AM1               ! } variation with latitude
!
       REAL CD_BOT  ! bottom friction (dimensionless) drag coeff
      REAL FKPH                                                         &
                             ! Coeff of vertical mixing of T
     &    ,FKPM                                                         &
                             ! Coeff of vertical mixing of U,V
     &    ,FNUB_SI                                                      &
                             ! }Constants used in calculating
     &    ,FNU0_SI                                                      &
                             ! }Richardson number dependent viscosity
     &    ,FKAPB_SI                                                     &
                             ! }
     &    ,STABLM_SI                                                    &
                             ! }
     &    ,GNUMINC_SI                                                   &
                             ! Min values of vertical
     &    ,GNUMINT_SI        ! diffusivity (level 1)
      REAL                                                              &
     & max_qLarge_depth                                                 &
                         ! max depth allowed for quadratic Large
     &,crit_Ri                                                          &
                ! critical Richardson no used for depth of
                ! quadratic Large scheme
     &,MAX_LARGE_DEPTH,CRIT_RI_FL  ! Equivalent numbers for Full Large
!
      REAL AHI1_SI,                                                     &
                             ! } Coeffs for isopycnal diffusion coeff
     &     AHI2_SI,                                                     &
                             ! } allowing it to vary with depth
     &     AHI3_SI,                                                     &
                             ! }
     &     SLOPE_MAX         ! Max slope allowed (isopyc diffusncalcs)
      REAL TMIN1_MAX,ATHKDFT_MIN,ATHKDFT_MAX
!
      REAL GAMMA_BBL,                                                   &
     & AH_Sigma,                                                        &
     & AH_Sigma_min
      REAL ATHKDF1_SI,                                                  &
                             ! } Coeffs for Gent and McWilliams scheme
     &     ATHKDF2_SI,                                                  &
                             ! } thickness diffusion allowing it to
     &     ATHKDF3_SI,                                                  &
                             ! } vary with depth
     &     dslope,                                                      &
                             !  Values required for tapering of
     &     slopec            !  isopycnals
      REAL athkdf_bi                                                    &
                            ! Coeff for biharmonic GM mixing
     &    ,BM               ! Coeff mixing for biharmonic U,V
      REAL HANEY_SI                                                     &
                             ! Haney temperature coeff (W/m2/K)
     &    ,HANEY_SALT        ! Haney salt coefficient
!
      REAL ETA1_SI,                                                     &
                             ! } Depth decay of solar radiation bands
     &     ETA2_SI,                                                     &
                             ! }
     &     RSOL              ! ratio of long/short wave radiation
      Integer KFIX           ! max penetration depth of solar radiation
!
      REAL LAMDA,                                                       &
                             ! Frctn wind mixing energy availbl to ocean
     &     DELTA_SI,                                                    &
                             ! Decay scale of wind mixing energy
     &     EPSILON,                                                     &
                             ! Proportion of convectivelyreleased
                             ! energy available for mixing
     &     DELPSF            ! Non-solar density change (cm3/s2)
!
      REAL                                                              &
     &     OBDY_GILL_MU,                                                &
                              ! weighting value for Gill bdy routine
     &     OBDY_GILL_LAMDA    ! constant for Gill ocean bdy routine

      REAL TFREEZE,                                                     &
                             ! Freezing pt of sea water
     &     EDDYDIFFN,                                                   &
                             ! Ocean/ice diffusion constant N hemisphr
     &     EDDYDIFFS,                                                   &
                             ! Ocean/ice diffusion constant S hemisphr
     &     SALICE,                                                      &
                             ! Salinity of sea ice
     &     QFUSION,                                                     &
                             ! Latent heat of fusion sea ice (J/kg)
     &     H0,                                                          &
                             ! Min local depth newly formed ice (m)
     &     AMXSOUTH,                                                    &
                             ! Max ice fraction S hemisphere
     &     AMXNORTH,                                                    &
                             ! Max ice fraction N hemisphere
     &     AICEMIN,                                                     &
                             ! Min ice fraction for icy points
     &     HICEMIN,                                                     &
                             ! Min ice depth for icy points
     &     RHOICE,                                                      &
                             ! Density of ice
     &     RHOSNOW,                                                     &
                             ! Density of snow
     &     RHOWATER,                                                    &
                             ! Density of water
     &     AH_ICE,                                                      &
                             ! Ice diffusion coefficient
     &     quad_ice_drag,                                               &
                           ! coefft of quadratic ice-ocean drag
     &     hicestop,                                                    &
                           ! max hice at which convergence is allowed
     &     hiceslow,                                                    &
                           ! min hice at which convergence is impeded
     &     aicemizfry,                                                  &
                            ! aice below which GBM O2I is constant
     &     MCPHEE_COEFF,                                                &
                            ! Ocean-ice heat transfer coefficient
     &     MCPHEE_MINUSTAR  ! Min friction vel for McPhee scheme
      REAL :: dte
      REAL :: cw
      REAL :: sw
      REAL :: ecc2
      REAL :: eyc
      REAL :: cstar
      REAL :: pstar
      REAL :: ISTRHSM_LAT ! Latitude above which the ice strength is
                          !   smoothed
      REAL :: HS_MIN      ! Min allowable snow thickness
      REAL :: HFRAZILMIN  ! Min thickness of new frazil ice
      REAL :: ICEPUNY     ! Tiny value used to check if ice present
      REAL :: GSTAR       ! Used to compute participation function
      REAL :: HSTAR       ! Determines mean thickness of ridged ice (m)
      REAL :: CF_RDG      ! Ratio of ridging work to PE change in ridging
      REAL :: CS_RDG      ! Fraction of shear energy contrbtng to ridging
      REAL :: FSNOW_RDG   ! Snow fraction that survives in ridging
      REAL :: ice_vel_factor ! scaling of ice vel from jmtm2 to jmtm1
!
      REAL TFREEZ,                                                      &
                             ! } Freezing pt of sea water and
     &     TCHECK            ! } value for .LE. test

      INTEGER FLUX_TRACER       ! Tracer to use for diagnostics
!
      REAL PCO2_ATM_0        ! Initial atmospheric pCO2
      REAL PC14_ATM_0        ! Initial atmospheric carbon 14 conc.

      REAL CSR_JFU0                                                     &
                     ! CSR(JFU0) needed for mpp filtering code
     &,    CSTR_JFT0 ! CSTR(JFT0) needed for mpp filtering code
!
      REAL tendfrc           ! Fraction of box mixed by Med outflow
       real SALREF,SALLOW,SALUP
!
      REAL med_in,med_out                                               &
     &    ,hud_in,hud_out   ! Hudson Bay inflow goes in at depth
      REAL LAMBDA_LARGE

      ! parameters for the implicit free surface code
      REAL :: ALPHA_FS
      REAL :: GAMMA_FS
      REAL :: THETA_FS
      REAL :: GAM_FS
      REAL :: GCOR
      REAL :: APGR
      REAL :: AREA2
#if !defined(SLAB)
      COMMON /SCALAR/ DTTS,DTUV,DTBT,DTSF,C2DTTS,C2DTUV,C2DTSF,C2DTBT   &
     &,AH,AM                                                            &
     &,AM0,AM1                                                          &
     &,FIRST_GUESS_FRAC,MSCAN,MXSCAN,MIN_SCANS                          &
     &,ACOR,CRIT,SOR,OMEGA,RADIUS,WGHT_DELPLUS,GRAV,RADIAN              &
     &,GRAV_SI                                                          &
     &,RHO_WATER_SI,SPECIFIC_HEAT_SI                                    &
     &,ROWPRT,TSUVPRT,ANCILPRT,L_OPRT_LSM,L_OPRT_SF_INDX                &
     &,ISTRT,ISTOP                                                      &
     &,FKPH,FKPM,CD_BOT                                                 &
     &,AHI1_SI,AHI2_SI,AHI3_SI,SLOPE_MAX,TMIN1_MAX                      &
     &,ATHKDFT_MIN,ATHKDFT_MAX                                          &
     &,ATHKDF1_SI,ATHKDF2_SI,ATHKDF3_SI,dslope,slopec                   &
     &,athkdf_bi,BM,BH                                                  &
     &,FNUB_SI,FNU0_SI,KAPPA0_SI,DKAPPA_DZ_SI,STABLM_SI                 &
     &,KAPPA_MIN_FRAC,KAPPA_LOW_DEPTH,GNUMINC_SI,GNUMINT_SI             &
                              ! Min values of vert diffusivity (level 1)
     &,max_qLarge_depth,crit_Ri,MAX_LARGE_DEPTH,CRIT_RI_FL              &
     &,HANEY_SI, HANEY_SALT                                             &
                            ! Haney coefficients
     &,ETA1_SI,ETA2_SI,RSOL,KFIX                                        &
     &,LAMDA,DELTA_SI,EPSILON                                           &
                                      ! Mixed layer params
     &,DELPSF                                                           &
               ! 'Non-solar density change' (mixed layer) (units cm3/s2)
     &,OBDY_GILL_MU,OBDY_GILL_LAMDA                                     &
     &,TFREEZE,EDDYDIFFN,EDDYDIFFS,SALICE,QFUSION,H0,AMXSOUTH,AMXNORTH  &
     &,AICEMIN,HICEMIN,RHOICE,RHOSNOW,RHOWATER                          &
     &,ah_ice,quad_ice_drag,hicestop,hiceslow,aicemizfry                &
     &,MCPHEE_COEFF,MCPHEE_MINUSTAR                                     &
     &,dte,cw,sw,ecc2,eyc,cstar,pstar,ISTRHSM_LAT                       &
     &,HS_MIN,HFRAZILMIN,ICEPUNY,GSTAR,HSTAR,CF_RDG,CS_RDG,FSNOW_RDG    &
     &,ice_vel_factor                                                   &
     &,TFREEZ,TCHECK                                                    &
                        ! Freezing pt and value for .LE. test
     &,pco2_atm_0                                                       &
                        ! Initial atmospheric pCO2
     &,FLUX_TRACER                                                      &
                        ! Tracer to use for diagnostics
     &,CSR_JFU0,CSTR_JFT0                                               &
     &,tendfrc,med_in,med_out,hud_in,hud_out                            &
     &,ALPHA_FS,GAMMA_FS,THETA_FS,GAM_FS,GCOR,APGR,AREA2                &

     &,SALREF,SALLOW,SALUP                                              &
     &,LAMBDA_LARGE                                                     &
     &,GAMMA_BBL                                                        &
     &,AH_Sigma                                                         &
     &,AH_Sigma_min
#endif


!===================================================================
