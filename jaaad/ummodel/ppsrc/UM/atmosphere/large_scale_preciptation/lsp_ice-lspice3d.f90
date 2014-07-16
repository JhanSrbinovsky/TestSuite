
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Precipitation microphysics calculations.
! Subroutine Interface:
      SUBROUTINE LSP_ICE(                                               &
     &  P,RHODZ,deltaz,rhodz_dry,rhodz_moist,                           &
     &  TIMESTEPFIXED,n_iterations,POINTS,L_MURK,RHCPT,                 &
     &  L_USE_SULPHATE_AUTOCONV, L_AUTO_DEBIAS,                         &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_psd, L_it_melting,     &
     &  l_non_hydrostatic, l_mixing_ratio, l_cry_agg_dep,               &
     &  SO4_ACC,SO4_DIS,SO4_AIT,                                        &
     &  L_BIOMASS_CCN, BMASS_AGD, BMASS_CLD,                            &
     &  L_OCFF_CCN, OCFF_AGD, OCFF_CLD,                                 &
     &  L_SEASALT_CCN, SEA_SALT_FILM, SEA_SALT_JET, SALT_DIM_ICE,       &
     &  L_biogenic_CCN, biogenic, biog_dim_ice,                         &
     &  SNOW_DEPTH,LAND_FRACT,AEROSOL,                                  &
     &  L_seq_mcr,L_autoc_3b,l_autolim_3b,L_autoconv_murk,              &
     &  l_droplet_settle, ec_auto,N_drop_land,                          &
     &  N_drop_sea,N_drop_land_cr,N_drop_sea_cr,Ntot_land, Ntot_sea,    &
     &  ai, bi, aic, bic,                                               &
     &  QCF,QCL,Q,QCF2,QRAIN,QGRAUP,                                    &
     &  RAINRATE, VF_RAIN, SNOW_AGG, VF_AGG,                            &
     &  SNOW_CRY, VF_CRY, GRAUPRATE, VF_GRAUP, droplet_flux,            &
     &  FRAC_ICE_ABOVE,CTTEMP,RAINFRAC,                                 &
     &  T,CFKEEP,CFLIQKEEP,CFICEKEEP,BLAND,CX,CONSTP                    &
     &, PSDEP,PSAUT,PSACW,PSACR,PSACI,PSMLT,PSMLTEVP                    &
     &, PRAUT,PRACW,PREVP                                               &
     &, PGAUT,PGACW,PGACS,PGMLT                                         &
     &, PIFRW,PIPRM,PIDEP,PIACW,PIACR,PIMLT,PIMLTEVP                    &
     &, PIFALL,PSFALL,PRFALL,PGFALL,PLSET,PLEVPSET                      &
     & )

      IMPLICIT NONE
!
! Description:
!   Updates ice, liquid and vapour contents, temperature and
!   cloud fractions as a result of microphysical processes.
!
! Method:
!   Calculates transfers of water between vapour, ice/snow,
!   cloud liquid water, rain and graupel.
!
!   Processes included are:
!   - Fall of hydrometeor into and out of the layer (sedimentation)
!   - Homogenous and heterogenous nucleation of ice;
!   - Deposition and sublimation of ice/snow;
!   - Autoconversion of ice->snow, snow->graupel, liquid->rain
!   - Collection processes
!   - Melting and evaporation
!
!   This is described in Unified Model Documentation Paper 26.
!
!   There are a number of different options for prognostic
!   hydrometeor variables. Each is independent of the others.
!   - Second prognostic cloud ice variables
!      Active if L_mcr_qcf2=.True. (in CNTLATM namelist)
!      The code supports the use of a second cloud ice prognostic
!      variable so that both cloud ice aggregates (QCF/QCF_AGG)
!      and cloud ice pristine crystals (QCF2/QCF_CRY) can be
!      represented and advected separately.
!      If False, then there is a diagnostic split within each level
!      at each timestep into ice crystals and snow aggregates.
!   - Prognostic rain
!      Active if L_mcr_qrain=.True. (in CNTLATM namelist)
!      If False, then rain is treated diagnostically.
!   - Prognostic graupel
!      Active if L_mcr_qgraup=.True. (in CNTLATM namelist)
!      If False, then graupel is not represented.
!
! Current Code Owner: Jonathan Wilkinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables:
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
! --------------------------COMDECK C_LSPMIC----------------------------
! SPECIFIES MICROPHYSICAL PARAMETERS FOR AUTOCONVERSION, HALLETT MOSSOP
! PROCESS, ICE NUCLEATION. ALSO SPECIFIES NUMBER OF ITERATIONS OF
! THE MICROPHYSICS AND ICE CLOUD FRACTION METHOD
! ----------------------------------------------------------------------
!
! History:
!
! Version    Date     Comment
! -------    ----     -------
!   5.4    16/08/02   Correct comment line, add PC2 parameters and
!                     move THOMO to c_micro     Damian Wilson
!   6.0    11/08/03   Correct value of wind_shear_factor for PC2
!                                                          Damian Wilson
!   6.2    17/11/05   Remove variables that are now in UMUI. D. Wilson
!   6.2    03/02/06   Include droplet settling logical. Damian Wilson
!
! ----------------------------------------------------------------------
!      AUTOCONVERSION TERMS
! ----------------------------------------------------------------------
!
!     LOGICAL, PARAMETER :: L_AUTOCONV_MURK is set in UMUI
! Set to .TRUE. to calculate droplet concentration from MURK aerosol,
! which will override L_USE_SULPHATE_AUTOCONV (second indirect effect
! of sulphate aerosol). If both are .FALSE., droplet concentrations
! from comdeck C_MICRO are used to be consistent with the values
      ! used in the radiation scheme.

      ! This next set of parameters is to allow the 3B scheme to
      ! be replicated at 3C/3D
        ! Inhomogeneity factor for autoconversion rate
        REAL,PARAMETER:: INHOMOG_RATE=1.0

        ! Inhomogeneity factor for autoconversion limit
        REAL,PARAMETER:: INHOMOG_LIM=1.0

        ! Threshold droplet radius for autoconversion
        REAL,PARAMETER:: R_THRESH=7.0E-6
      ! End of 3B repeated code

      !Do not alter R_AUTO and N_AUTO since these values are effectively
      ! hard wired into a numerical approximation in the autoconversion
      ! code. EC_AUTO will be multiplied by CONSTS_AUTO

      ! Threshold radius for autoconversion
      REAL, PARAMETER :: R_AUTO=20.0E-6

      ! Critical droplet number for autoconversion
      REAL, PARAMETER :: N_AUTO=1000.0

      ! Collision coalesence efficiency for autoconversion
!      REAL, PARAMETER :: EC_AUTO is set in UMUI

      ! The autoconversion powers define the variation of the rate with
      ! liquid water content and droplet concentration. The following are
      ! from Tripoli and Cotton

      !  Dependency of autoconversion rate on droplet concentration
      REAL, PARAMETER :: POWER_DROPLET_AUTO=-0.33333

      ! Dependency of autoconversion rate on water content
      REAL, PARAMETER :: POWER_QCL_AUTO=2.33333

      ! Dependency of autoconversion rate on air density
      REAL, PARAMETER :: power_rho_auto=1.33333

      ! CONSTS_AUTO = (4 pi)/( 18 (4 pi/3)^(4/3)) g /  mu (rho_w)^(1/3)
      ! See UM documentation paper 26, equation P26.132

      ! Combination of physical constants
      REAL, PARAMETER :: CONSTS_AUTO=5907.24

      ! Quantites for calculation of drop number by aerosols.
      ! Need only set if L_AUTOCONV_MURK=.TRUE.  See file C_VISBTY

      ! Scaling concentration (m-3) in droplet number concentration
      REAL, PARAMETER :: N0_MURK=500.0E6

      ! Scaling mass (kg/kg) in droplet number calculation from aerosols
      REAL, PARAMETER :: M0_MURK=1.458E-8

      ! Power in droplet number calculation from aerosols
      REAL, PARAMETER :: POWER_MURK=0.5

      ! Ice water content threshold for graupel autoconversion (kg/m^3)
      REAL, PARAMETER :: AUTO_GRAUP_QCF_THRESH = 3.E-4

      ! Temperature threshold for graupel autoconversion (degC)
      REAL, PARAMETER :: AUTO_GRAUP_T_THRESH = -4.0

      ! Temperature threshold for graupel autoconversion
      REAL, PARAMETER :: AUTO_GRAUP_COEFF = 0.5

      !-----------------------------------------------------------------
      ! Iterations of microphysics
      !-----------------------------------------------------------------

      ! Number of iterations in microphysics.
      INTEGER,PARAMETER :: LSITER=1
      ! Advise 1 iteration for every 10 minutes or less of timestep.

      !-----------------------------------------------------------------
      ! Nucleation of ice
      !-----------------------------------------------------------------

      ! Note that the assimilation scheme uses temperature thresholds
      ! in its calculation of qsat.

      ! Nucleation mass
      REAL, PARAMETER :: M0=1.0E-12

      ! Maximum Temp for ice nuclei nucleation (deg C)
      REAL, PARAMETER :: TNUC=-10.0

      ! Maximum temperature for homogenous nucleation is now in c_micro
      ! so that it is available to code outside of section A04.

      !  1.0/Scaling quantity for ice in crystals
      REAL, PARAMETER :: QCF0=1.0E4       ! This is an inverse quantity

      ! Minimum allowed QCF after microphysics
      REAL,PARAMETER:: QCFMIN=1.0E-8

      ! 1/scaling temperature in aggregate fraction calculation
      REAL, PARAMETER :: T_SCALING=0.0384

      !  Minimum temperature limit in calculation  of N0 for ice (deg C)
      REAL, PARAMETER :: T_AGG_MIN=-45.0

      !-----------------------------------------------------------------
      ! Hallett Mossop process
      !-----------------------------------------------------------------

      ! Switch off Hallett Mossop in this version but allow
      ! functionality

      ! Min temp for production of Hallett Mossop splinters (deg C)
      REAL, PARAMETER :: HM_T_MIN=-8.0

      ! Max temp for production of Hallett Mossop splinters (deg C)
      REAL, PARAMETER :: HM_T_MAX=-273.0
      ! REAL, PARAMETER :: HM_T_MAX=-3.0

      !  Residence distance for Hallett Mossop splinters (1/deg C)
      REAL, PARAMETER :: HM_DECAY=1.0/7.0

      ! Reciprocal of scaling liquid water content for HM process
      REAL, PARAMETER :: HM_RQCL=1.0/0.1E-3

      !-----------------------------------------------------------------
      ! PC2 Cloud Scheme Terms
      !-----------------------------------------------------------------

      ! Specifies the ice content (in terms of a fraction of qsat_liq)
      ! that corresponds to a factor of two reduction in the width of
      ! the vapour distribution in the liquid-free part of the gridbox.
      REAL, PARAMETER :: ICE_WIDTH=0.04

      ! Parameter that governs the rate of spread of ice cloud fraction
      ! due to windshear
      REAL, PARAMETER :: WIND_SHEAR_FACTOR = 1.5E-4

!
! Description:
!
!  Contains various cloud droplet parameters, defined for
!  land and sea areas.
!
!  NTOT_* is the total number concentration (m-3) of cloud droplets;
!  KPARAM_* is the ratio of the cubes of the volume-mean radius
!                                           and the effective radius;
!  DCONRE_* is the effective radius (m) for deep convective clouds;
!  DEEP_CONVECTION_LIMIT_* is the threshold depth (m) bewteen shallow
!                                          and deep convective cloud.
!
! Current Code Owner: Andy Jones
!
! History:
!
! Version   Date     Comment
! -------   ----     -------
!    1     040894   Original code.    Andy Jones
!  5.2     111000   Updated in line with Bower et al. 1994 (J. Atmos.
!                   Sci., 51, 2722-2732) and subsequent pers. comms.
!                   Droplet concentrations now as used in HadAM4.
!                                     Andy Jones
!  5.4     02/09/02 Moved THOMO here from C_LSPMIC.      Damian Wilson
!  6.2     17/11/05 Remove variables that are now in UMUI. D. Wilson
!
!     REAL,PARAMETER:: NTOT_LAND is set in UMUI
!     REAL,PARAMETER:: NTOT_SEA is set in UMUI
      REAL,PARAMETER:: KPARAM_LAND = 0.67
      REAL,PARAMETER:: KPARAM_SEA = 0.80
      REAL,PARAMETER:: DCONRE_LAND = 9.5E-06
      REAL,PARAMETER:: DCONRE_SEA = 16.0E-06
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_LAND = 500.0
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_SEA = 1500.0
!
! Maximum Temp for homogenous nucleation (deg C)
      REAL,PARAMETER:: THOMO = -40.0
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
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
! --------------------------COMDECK C_LSPDIF---------------------------
! input variables
!
! Values of reference variables
      REAL,PARAMETER:: AIR_DENSITY0=1.0             ! kg m-3
      REAL,PARAMETER:: AIR_VISCOSITY0=1.717E-5      ! kg m-1 s-1
      REAL,PARAMETER:: AIR_CONDUCTIVITY0=2.40E-2    ! J m-1 s-1 K-1
      REAL,PARAMETER:: AIR_DIFFUSIVITY0=2.21E-5     ! m2 s-1
      REAL,PARAMETER:: AIR_PRESSURE0=1.0E5          ! Pa

! Values of diffusional growth parameters
      ! Terms in deposition and sublimation
      REAL,PARAMETER:: APB1=(LC+LF)**2 * EPSILON /(R*AIR_CONDUCTIVITY0)
      REAL,PARAMETER:: APB2=(LC+LF) / AIR_CONDUCTIVITY0
      REAL,PARAMETER:: APB3=R/(EPSILON*AIR_PRESSURE0*AIR_DIFFUSIVITY0)
      ! Terms in evap of melting snow and rain
      REAL,PARAMETER:: APB4=LC**2*EPSILON/(R*AIR_CONDUCTIVITY0)
      REAL,PARAMETER:: APB5=LC /AIR_CONDUCTIVITY0
      REAL,PARAMETER:: APB6=R/(EPSILON*AIR_PRESSURE0*AIR_DIFFUSIVITY0)

! Values of numerical approximation to wet bulb temperature
      ! Numerical fit to wet bulb temperature
      REAL,PARAMETER:: TW1=1329.31
      REAL,PARAMETER:: TW2=0.0074615
      REAL,PARAMETER:: TW3=0.85E5
      ! Numerical fit to wet bulb temperature
      REAL,PARAMETER:: TW4=40.637
      REAL,PARAMETER:: TW5=275.0

! Ventilation parameters
      REAL,PARAMETER:: SC=0.6
! f(v)  =  VENT_ICE1 + VENT_ICE2  Sc**(1/3) * Re**(1/2)
      REAL,PARAMETER:: VENT_ICE1=0.65
      REAL,PARAMETER:: VENT_ICE2=0.44
! f(v)  =  VENT_RAIN1 + VENT_RAIN2  Sc**(1/3) * Re**(1/2)
      REAL,PARAMETER:: VENT_RAIN1=0.78
      REAL,PARAMETER:: VENT_RAIN2=0.31
! c_lspdif will call c_r_cp, c_epslon and c_lheat
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
!
! Subroutine arguments
!
! Obtain the size for CX and CONSTP
! Start C_LSPSIZ
! Description: Include file containing idealised forcing options
! Author:      R. Forbes
!
! History:
! Version  Date      Comment
! -------  ----      -------
!   6.1    01/08/04  Increase dimension for rain/graupel.  R.Forbes
!   6.2    22/08/05  Include the step size between ice categories.
!                                                   Damian Wilson

! Sets up the size of arrays for CX and CONSTP
      REAL CX(100),CONSTP(100)
      INTEGER,PARAMETER:: ice_type_offset=20

! End C_LSPSIZ
!
      LOGICAL                                                           &
                      !, INTENT(IN)
     &  L_MURK                                                          &
!         Murk aerosol is a valid quantity
     & ,L_USE_SULPHATE_AUTOCONV                                         &
!         Switch to use sulphate aerosol in the calculation of cloud
!         droplet number concentration in the autoconversion section,
!         i.e. activate the second indirect ("lifetime") effect.
     & ,L_SEASALT_CCN                                                   &
!         Switch to supplement sulphate aerosol with sea-salt if the
!         second indirect effect is active.
     & ,L_BIOMASS_CCN                                                   &
!         Switch to supplement sulphate aerosol with biomass smoke
!         aerosol if the second indirect effect is active.
     & ,L_OCFF_CCN                                                      &
!         Switch to supplement sulphate aerosol with fossil-fuel organic
!         carbon aerosol if the second indirect effect is active.
     & ,L_BIOGENIC_CCN                                                  &
!         Switch to supplement sulphate aerosol with biogenic aerosol
!         if the second indirect effect is active.
     & ,L_AUTO_DEBIAS                                                   &
!         Switch to apply de-biasing correction to autoconversion rate.
     &, L_mcr_qcf2                                                      &
                       ! true if using 2nd cloud ice prognostic
     &, L_mcr_qrain                                                     &
                       ! true if using prognostic rain
     &, L_mcr_qgraup                                                    &
                       ! true if using prognostic graupel
     &, L_psd                                                           &
                       ! true if using generic ice particle size distn
     &, L_it_melting                                                    &
                       ! true if using iterative melting
     &, L_non_hydrostatic                                               &
                           ! Use non-hydrostatic
                           ! calculation of model layer depths
     &, L_mixing_ratio     ! q is a mixing ratio
!
      INTEGER                                                           &
                      !, INTENT(IN)
     &  POINTS                                                          &
!         Number of points to be processed.
     &, n_iterations                                                    &
                           ! Number of iterations for iterative melting
     & ,SALT_DIM_ICE                                                    &
!        Number of points for sea-salt arrays (either POINTS or 1)
     &, BIOG_DIM_ICE
!        Number of points for biogenic array (either POINTS or 1)
!

      LOGICAL                                                           &
                      !, INTENT(IN)
     &  BLAND(POINTS)
!         Land/sea mask
!
      REAL                                                              &
                      !, INTENT(IN)
     &  TIMESTEPFIXED                                                   &
!         Timestep of physics in model (s).
     &, CFLIQKEEP(POINTS)                                               &
!         Liquid cloud fraction in this layer (no units).
     &, CFICEKEEP(POINTS)                                               &
!         Frozen cloud fraction in this layer (no units).
     &, CFKEEP(POINTS)                                                  &
!         Total cloud fraction in this layer (no units).
     &, P(POINTS)                                                       &
!         Air pressure at this level (Pa).
     &, RHODZ(POINTS)                                                   &
!         Air mass p.u.a. in this layer (kg per sq m).
     &, deltaz(points)                                                  &
                           ! Thickness of layer (m)
     &, rhodz_dry(points)                                               &
                           ! Dry air density * layer thickness (kg m-2)
     &, rhodz_moist(points)                                             &
                           ! Moist air density * layer thick.  (kg m-2)
     &, SO4_ACC(POINTS)                                                 &
!         Sulphur cycle variable
     &, SO4_DIS(POINTS)                                                 &
!         Sulphur cycle variable
     &, SO4_AIT(POINTS)                                                 &
!         Sulphur cycle variable
     &, BMASS_AGD(POINTS)                                               &
!         Aged biomass smoke mass mixing ratio
     &, BMASS_CLD(POINTS)                                               &
!         In-cloud biomass smoke mass mixing ratio
     &, OCFF_AGD(POINTS)                                                &
!         Aged fossil-fuel organic carbon mass mixing ratio
     &, OCFF_CLD(POINTS)                                                &
!         In-cloud fossil-fuel organic carbon mass mixing ratio
     &, SNOW_DEPTH(POINTS)                                              &
!         Snow depth for aerosol amount (m)
     &, LAND_FRACT(POINTS)                                              &
!         Land fraction
     &, AEROSOL(POINTS)                                                 &
!         Aerosol mass (ug/kg)
     &, RHCPT(POINTS)                                                   &
!         Critical relative humidity of all points for cloud formation.
     &, SEA_SALT_FILM(SALT_DIM_ICE)                                     &
!         Film-mode sea-salt aerosol number concentration (m-3)
     &, SEA_SALT_JET(SALT_DIM_ICE)                                      &
!         Jet-mode sea-salt aerosol number concentration (m-3)
     &, BIOGENIC(BIOG_DIM_ICE)
!         biogenic aerosl mass mixing ratio

      REAL                                                              &
                      !, INTENT(INOUT)
     &  Q(POINTS)                                                       &
!         Specific humidity at this level (kg water per kg air).
     &, QCF(POINTS)                                                     &
!         Cloud ice (kg water per kg air).
     &, QCL(POINTS)                                                     &
!         Cloud liquid water (kg water per kg air).
     &, QCF2(POINTS)                                                    &
!         Second cloud ice (kg water per kg air).
     &, QRAIN(POINTS)                                                   &
!         Rain (kg water per kg air).
     &, QGRAUP(POINTS)                                                  &
!         Graupel (kg water per kg air).
     &, T(POINTS)
!         Temperature at this level (K).

!     Hydrometeor flux between layers (kg m-2 s-1)
!     On input: Hydrometeor flux entering this layer from above.
!     On output: Hydrometeor flux leaving this layer.
!     Note: If only one ice prognostic is active (L_mcr_qcf2=.F.)
!     then SNOW_AGG contains all the ice/snow and SNOW_CRY is zero.
!
      REAL                                                              &
     &  SNOW_AGG(POINTS)                                                &
                              ! snow aggregates
     &, SNOW_CRY(POINTS)                                                &
                              ! ice crystals
     &, RAINRATE(POINTS)                                                &
                              ! rain
     &, GRAUPRATE(POINTS)                                               &
                              ! graupel
     &, droplet_flux(points)  ! droplets

!     Hydrometeor mass-weighted fall velocities (m s-1)
!     On input: Fall velocity of hydrometeor entering layer.
!     On Output: Fall velocity of hydrometeor leaving layer.
!     Note: If only one ice prognostic is active, then only
!     VF_AGG is used.
!
      REAL                                                              &
     &  VF_AGG(POINTS)                                                  &
                           ! snow aggregates
     &, VF_CRY(POINTS)                                                  &
                           ! ice crystals
     &, VF_RAIN(POINTS)                                                 &
                           ! rain
     &, VF_GRAUP(POINTS)   ! graupel

!     Cloud/precipitation sub-grid fraction variables
!
      REAL                                                              &
     &  CTTEMP(POINTS)                                                  &
!         Ice cloud top temperature (K)
     &, RAINFRAC(POINTS)                                                &
!         Rain fraction (no units)
     &, FRAC_ICE_ABOVE(POINTS)
!         Fraction of ice in layer above (no units)
!
      Logical                                                           &
                      !, Intent(IN)
     &  L_seq_mcr                                                       &
                               ! Use sequential updating
     &, L_autoc_3b                                                      &
                               ! Use 3B autoconversion method
     &, l_autolim_3b                                                    &
                               ! Use fixed 3B values for the
                               ! autoconversion limit 
     &, L_autoconv_murk                                                 &
                               ! Use murk aerosol to calc. drop number
     &, L_cry_agg_dep                                                   &
                               ! Limit the supersaturation that can be
                               ! removed by deposition depending on
                               ! amount of ice in each category
     &, l_droplet_settle       ! Allow cloud droplets to settle

      Real                                                              &
                      !, Intent(IN)
     &  ec_auto                                                         &
                               ! Collision coalescence efficiency
     &, N_drop_land                                                     &
                               ! Dummy in 3D version
     &, N_drop_sea                                                      &
                               ! Dummy in 3D version
     &, N_drop_land_cr                                                  &
                               ! Dummy in 3D version
     &, N_drop_sea_cr                                                   &
                               ! Dummy in 3C version
     &, Ntot_land                                                       &
                               ! Number of droplets over land / m-3
     &, Ntot_sea                                                        &
                               ! Number of droplets over sea / m-3
     &, ai, bi, aic, bic
                               ! Ice particle mass-size relationships
!
! Microphysical process rate diagnostics
      REAL, Intent(InOut) ::                                            &
     &  PSDEP(POINTS)                                                   &
                       ! Deposition of vapour to snow aggregates
     &, PSAUT(POINTS)                                                   &
                       ! Autoconversion of aggregates from crystals
     &, PSACW(POINTS)                                                   &
                       ! Accretion of liq. water by snow aggregates
     &, PSACR(POINTS)                                                   &
                       ! Collection of rain by snow aggregates
     &, PSACI(POINTS)                                                   &
                       ! Collection of ice crystals by aggregates
     &, PSMLT(POINTS)                                                   &
                       ! Melting of snow aggregates
     &, PSMLTEVP(POINTS)  ! Evaporation of melting aggregates
      REAL, Intent(InOut) ::                                            &
     &  PRAUT(POINTS)                                                   &
                       ! Autoconversion of cloud drops to rain
     &, PRACW(POINTS)                                                   &
                       ! Accretion of liq. water by rain
     &, PREVP(POINTS)  ! Evaporation of rain
      REAL, Intent(InOut) ::                                            &
     &  PGAUT(POINTS)                                                   &
                       ! Autoconversion of graupel from aggregates
     &, PGACW(POINTS)                                                   &
                       ! Accretion of liq. water by graupel
     &, PGACS(POINTS)                                                   &
                       ! Collection of snow aggregates by graupel
     &, PGMLT(POINTS)  ! Melting of graupel
      REAL, Intent(InOut) ::                                            &
     &  PIFRW(POINTS)                                                   &
                       ! Homogeneous freezing nucleation
     &, PIPRM(POINTS)                                                   &
                       ! Heterogeneous (primary) nucleation
     &, PIDEP(POINTS)                                                   &
                       ! Deposition of vapour to ice crystals
     &, PIACW(POINTS)                                                   &
                       ! Accretion of liq. water by ice crystals
     &, PIACR(POINTS)                                                   &
                       ! Collection of rain by ice crystals
     &, PIMLT(POINTS)                                                   &
                       ! Melting of ice crystals
     &, PIMLTEVP(POINTS)  ! Evaporation of melting ice crystals
       REAL, Intent(InOut) ::                                           &
     &  PIFALL(POINTS)                                                  &
                       ! Sedimentation of ice crystals
     &, PSFALL(POINTS)                                                  &
                       ! Sedimentation of aggregates
     &, PRFALL(POINTS)                                                  &
                       ! Sedimentation of rain
     &, PGFALL(POINTS) ! Sedimentation of graupel
       REAL, Intent(InOut) ::                                           &
     &  PLSET(POINTS)                                                   &
                       ! Droplet settling of liquid water
     &, PLEVPSET(POINTS) ! Evaporated settled droplets
!
!  Local parameters
      REAL LCRCP,LFRCP,LSRCP,RHO1,ONE_OVER_EPSILON,ONE_OVER_ZERODEGC
      PARAMETER(                                                        &
     &  LCRCP=LC/CP                                                     &
!         Latent heat of condensation / cp (K).
     &, LFRCP=LF/CP                                                     &
!         Latent heat of fusion / cp (K).
     &, LSRCP=LCRCP+LFRCP                                               &
!         Sum of the above (S for Sublimation).
     &, RHO1=1.0                                                        &
!         Reference density of air (kg m-3)
     &, ONE_OVER_EPSILON=1.0/EPSILON                                    &
!        Inverse of epsilon to speed up calculations
     &, ONE_OVER_ZERODEGC=1.0/ZERODEGC                                  &
!        Inverse of zero degrees Celsius to speed up code (K-1)
     &  )
!
!  Local scalars and dynamic arrays
!
      INTEGER                                                           &
     &  I                                                               &
!         Loop counter (horizontal field index).
     &, J                                                               &
!         Counter for the iterations
     &, KK                                                              &
!         Variable for condensed points compression
     &, KK2                                                             &
!         Another condensed points compression variable
     &, SEA_SALT_PTR                                                    &
!         Pointer for sea-salt arrays
     &, BIOMASS_PTR                                                     &
!         Pointer for biomass smoke arrays
     &, OCFF_PTR                                                        &
!         Pointer for fossil-fuel organic carbon arrays
     &, BIOGENIC_PTR
!         Pointer for biogenic aerosol arrays
!
      REAL                                                              &
     &  QS(POINTS)                                                      &
!         Saturated sp humidity for (T,p) in layer (kg kg-1)
     &, QSL(POINTS)                                                     &
!         Saturated sp humidity for (T,p) in layer
!         wrt water at all temps (kg kg-1)
     &, QS_CTT(POINTS)                                                  &
!         Saturated sp humidity for (CTTEMP,p) in layer (kg kg-1)
     &, QSL_CTT(POINTS)
!         Saturated sp humidity for (CTTEMP,p) in layer (kg kg-1)
!         wrt water at all temps

!     Cumulative fall out of hydrometeor within iterations (kg m-2 s-1)
      REAL                                                              &
     &  SNOWT_AGG(POINTS)                                               &
     &, SNOWT_CRY(POINTS)                                               &
     &, RAINRATET(POINTS)                                               &
     &, GRAUPRATET(POINTS)

      REAL                                                              &
     &  RHO(POINTS)                                                     &
!         Density of air in the layer (kg m-3).
     &, RHOR(POINTS)                                                    &
!         1.0/RHO to speed up calculations (kg-1 m3).
     &, VTEMP                                                           &
!         Virtual temperature as at start of loop (K).
     &, ESI(POINTS)                                                     &
!         saturation vapour pressure (wrt ice below zero Celsius)(Pa)
     &, ESW(POINTS)                                                     &
!         saturation vapour pressure (wrt water at all temperatures)(Pa)
     &, ESI_CTT(POINTS)                                                 &
!         saturation vapour pressure at cloud top temp (wrt ice) (Pa)
     &, ESW_CTT(POINTS)                                                 &
!         saturation vapour pressure at cloud top temp (wrt water) (Pa)
     &, DQI(POINTS)                                                     &
!         increment to/from ice/snow (kg kg-1)
     &, DQIL(POINTS)                                                    &
!         increment to/from cloud water (kg kg-1)
     &, DPR(POINTS)
!         increment to/from rain (kg m-2 s-1)

      REAL                                                              &
     &  CFICE(POINTS)                                                   &
!         fraction of ice inferred for the microphysics (no units).
     &, CFICEI(POINTS)                                                  &
!         inverse of CFICE (no units)
     &, CFLIQ(POINTS)                                                   &
!         liquid cloud fraction for the microphysics
     &, CF(POINTS)                                                      &
!         total cloud fraction for the microphysics (no units)
     &, DHI(POINTS)                                                     &
!         CFL limit (s m-1)
     &, DHIR(POINTS)                                                    &
!         1.0/DHI (m s-1)
     &, DHILSITERR(POINTS)
!         1.0/(DHI*LSITER) (m s-1)

!     Fallspeeds (m s-1)
      REAL                                                              &
     &  FQI_AGG(POINTS)                                                 &
                                ! snow aggregates
     &, FQI_CRY(POINTS)                                                 &
                                ! ice crystals
     &, FQI_RAIN(POINTS)                                                &
                                ! rain
     &, FQI_GRAUP(POINTS)       ! graupel

!     Saved hydrometeor fluxes out of layer (kg m-2 s-1)
      REAL                                                              &
     &  FQIRQI_AGG                                                      &
                                ! snow aggregates
     &, FQIRQI_CRY                                                      &
                                ! ice crystals
     &, FQIRQI_RAIN                                                     &
                                ! rain
     &, FQIRQI_GRAUP            ! graupel

!     Saved fluxes out of layer from layer above (kg m-2 s-1)
      REAL                                                              &
     &  FQIRQI2_AGG(POINTS)                                             &
                                ! snow aggregates
     &, FQIRQI2_CRY(POINTS)                                             &
                                ! ice crystals
     &, FQIRQI2_RAIN(POINTS)                                            &
                                ! rain
     &, FQIRQI2_GRAUP(POINTS)   ! graupel

      REAL                                                              &
     &  QCLNEW(POINTS)                                                  &
!         updated liquid cloud in implicit calculations (kg kg-1)
     &, TEMP7(POINTS)                                                   &
!         temporary variable
     &, TEMPW(POINTS)                                                   &
!         temporary for vapour calculations
     &, PR02(POINTS)                                                    &
!         term in evaporation of rain
     &, PR04(POINTS)                                                    &
!         square of pr02
     &, QC(POINTS)                                                      &
!         term in autoconversion of cloud to rain (kg kg-1)
     &, APLUSB(POINTS)                                                  &
!         denominator in deposition or evaporation of ice
     &, CORR(POINTS)                                                    &
!         density correction for fall speed (no units)
     &, ROCOR(POINTS)                                                   &
!         density correction for fall speed (no units)
     &, VS1(POINTS)                                                     &
!         Mean fall speed of snow (m s-1)
     &, VI1(POINTS)                                                     &
!         Mean fall speed of ice crystals (m s-1)
     &, VG1(POINTS)                                                     &
!         Mean fall speed of graupel (m s-1)
     &, FV1(POINTS)
!         Mean velocity difference between rain and snow (m s-1)

      REAL                                                              &
     &  LAMFAC1(POINTS)                                                 &
!         Expression containing calculations with lambda
     &, LAMS1(POINTS)                                                   &
!         Inverse lambda in snow exponential distribution (m)
     &, LAMI1(POINTS)                                                   &
!         Inverse lambda in ice crystal exponential distribution (m)
     &, LAMG1(POINTS)
!         Inverse lambda in graupel exponential distribution (m)

      REAL                                                              &
     &  TIMESTEP                                                        &
!         Timestep of each iteration (s)
     &, CORR2(POINTS)                                                   &
!         Temperature correction of viscosity etc. (no units)
     &, RHNUC                                                           &
!         Relative humidity required for nucleation (no units)
     &, TCG(POINTS)                                                     &
!         Temperature Factor for aggregate size distribution (no units)
     &, TCGI(POINTS)                                                    &
!         Inverse of TCG (no units)
     &, TCGC(POINTS)                                                    &
!         Temperature Factor for crystal size distribution (no units)
     &, TCGCI(POINTS)                                                   &
!         Inverse of TCGC (no units)
     &, TCGG(POINTS)                                                    &
!         Temperature Factor for graupel size distribution (no units)
     &, TCGGI(POINTS)                                                   &
!         Inverse of TCGC (no units)
     &, RATEQS(POINTS)                                                  &
!         Sub grid model variable (no units)
     &, HM_NORMALIZE                                                    &
!         Normalization for Hallett Mossop process (no units)
     &, HM_RATE
!         Increase in deposition due to Hallett Mossop process(no units)

      REAL                                                              &
     &  AREA_LIQ(POINTS)                                                &
!         Liquid only area of gridbox (no units)
     &, AREA_MIX(POINTS)                                                &
!         Mixed phase area of gridbox (no units)
     &, AREA_ICE(POINTS)                                                &
!         Ice only area of gridbox (no units)
     &, AREA_CLEAR(POINTS)                                              &
!         Cloud free area of gridbox (no units)
     &, RAIN_LIQ(POINTS)                                                &
!         Overlap fraction of gridbox between rain and liquid cloud
     &, RAIN_MIX(POINTS)                                                &
!         Overlap fraction of gridbox between rain and mixed phase cloud
     &, RAIN_ICE(POINTS)                                                &
!         Overlap fraction of gridbox between rain and ice cloud
     &, RAIN_CLEAR(POINTS)                                              &
!         Overlap fraction of gridbox between rain and no cloud
     &, Q_ICE(POINTS)                                                   &
!         Vapour content in the ice only part of the grid box (kg kg-1)
     &, Q_CLEAR(POINTS)                                                 &
!         Vapour content in the cloud free part of the grid box(kg kg-1)
     &, QCF_AGG(POINTS)                                                 &
!         QCF in the form of aggregates (kg kg-1)
     &, QCF_CRY(POINTS)                                                 &
!         QCF in the form of crystals (kg kg-1)
     &, FRAC_AGG(POINTS)
!         Fraction of aggregates (no units)

      REAL                                                              &
     &  TEMP1(POINTS)                                                   &
     &, TEMP2(POINTS)                                                   &
     &, TEMP3(POINTS)                                                   &
!         Temporary arrays for T3E vector functions
     &, POWER
!         Power for T3E vector functions

      REAL                                                              &
     &  QCFAUTORATE(POINTS)                                             &
!         Rate constant for ice to snow autoconversion
     &, QCFAUTOLIM(POINTS)                                              &
!         Limit for ice to snow autoconversion
     &, QCF_TOT(POINTS)                                                 &
!         Total amount of ice (crystals+aggregates)
     &, Frac_cry_dep(points)                                            &
!         Fraction of supersaturation that can be removed by crystals
     &, Frac_agg_dep(points)
!         Fraction of supersaturation that can be removed by aggregates

      REAL                                                              &
     &  LHEAT_CORREC_LIQ(POINTS)                                        &
!         Reduction factor in evaporation limits because of latent heat
     &, LHEAT_CORREC_ICE(POINTS)                                        &
!         Reduction factor in evaporation limits because of latent heat
     &, OVERHANG                                                        &
!         Amount of overhang of ice cloud in the layer above to that
!         below
     &, DQI_DEP(POINTS)                                                 &
     &, DQI_SUB(POINTS)                                                 &
     &, TEMPW_DEP(POINTS)                                               &
     &, TEMPW_SUB(POINTS)                                               &
     &, WIDTH(POINTS)                                                   &
     &, Q_ICE_1(POINTS), Q_ICE_2(POINTS)                                &
     &, AREA_ICE_1(POINTS), AREA_ICE_2(POINTS)                          &
!         Subgrid splitting for deposition term
     &, TEMPR                                                           &
!         Temporary for optimization
     &, qcft(points)   ! Holds total ice content qcf_cry + qcf_agg
!

! Cloud and rain fraction transfer rate diagnostics
      REAL                                                              &
     & cf_transfer_diag(points)                                         &
                                  ! Dummy to receive cf increment
     &,cfl_transfer_diag(points)                                        &
                                  ! Dummy to receive cfl increments
     &,cff_transfer_diag(points)                                        &
                                  ! Dummy to receive cff increments
     &,rf_transfer_diag(points)   ! Dummy to receive rainfrac increment

      LOGICAL,PARAMETER:: l_update_cf=.true. ! Update cloud fractions
      LOGICAL,PARAMETER:: l_seq=.true.       ! Sequential updating
      Logical l_crystals                     ! A crystal category exists

! Variables for calls to lsp_collection
      Logical                                                           &
     & l_use_area                                                       &
!       Use the ice partition to calculate transfer rates rather than
!       assuming a uniform distribution through the gridbox
     &,l_no_t_check
!       Do not check that the temperature is below zero degrees Celsius

      Integer                                                           &
     & ice_type1, ice_type2
!       Category of ice involved in collision process:
!       0 - crystals; 1 - aggregates; 3 - graupel.
!
!  Function calls
      REAL, EXTERNAL :: NUMBER_DROPLET
!
!- End of header
!
! ======================================================================
!       Initialize variables and calculate microphysical quantities
! ======================================================================
! DEPENDS ON: lsp_init
       Call lsp_init(points, timestepfixed, timestep                    &
     &,              l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup              &
     &,              l_psd, l_crystals                                  &
     &,              l_non_hydrostatic, l_mixing_ratio, l_cry_agg_dep   &
     &,              cx, constp, T, p, cttemp                           &
     &,              deltaz, rhodz,rhodz_dry,rhodz_moist, rho,rhor      &
     &,              dhi, dhir, dhilsiterr                              &
     &,              q, qcl, qcf, qcf2, qrain, qgraup, rainrate         &
     &,              qcf_agg, qcf_cry                                   &
     &,              qcf_tot, frac_agg, frac_cry_dep, frac_agg_dep      &
     &,              qs, qsl, esi, esw                                  &
     &,              psdep, psaut, psacw, psacr, psaci, psmlt           &
     &,              psmltevp, psfall, pifrw, piprm, pidep              &
     &,              piacw, piacr, pimlt, pimltevp, pifall              &
     &,              praut, pracw, prevp, prfall, plset, plevpset       &
     &,              pgaut, pgacw, pgacs, pgmlt, pgfall                 &
     &,              cf_transfer_diag, cfl_transfer_diag                &
     &,              cff_transfer_diag, rf_transfer_diag                &
     &,              snow_cry, snow_agg, snowt_cry, snowt_agg           &
     &,              rainratet, graupratet                              &
     &,              lheat_correc_liq, lheat_correc_ice                 &
     &,              corr, corr2, rocor                                 &
     &,              tcg, tcgi, tcgc, tcgci, tcgg, tcggi                &
     &              )

! ======================================================================
!       Start iterations
! ======================================================================
       DO J=1,LSITER

       Do i=1,points

         ! -------------------------------------------------------------
         ! Check that ice cloud fraction is sensible.
         ! -------------------------------------------------------------

         If (l_update_cf) Then
           CFICE(I)  = MAX(CFICEKEEP(I),0.001)
         Else
           CFICE(I)=MAX(MAX(CFICEKEEP(I),0.001),FRAC_ICE_ABOVE(I))
         End if
         CFICEI(I) = 1.0/CFICE(I)
         CF(I)     = CFKEEP(I)
         CFLIQ(I)  = CFLIQKEEP(I)
         CF(I)     = MIN( MAX(CF(I),CFICE(I)) ,(CFICE(I)+CFLIQ(I)) )

         ! -------------------------------------------------------------
         ! Calculate overlaps of liquid, ice and rain fractions
         ! -------------------------------------------------------------

         AREA_ICE(I)   = MAX(CF(I)-CFLIQ(I),0.0)
         AREA_CLEAR(I) = MAX(1.0-CF(I),0.0)
!
         ! -------------------------------------------------------------
         ! Update ice cloud top temperature if no ice falling in
         ! -------------------------------------------------------------

         IF (SNOW_CRY(I)+SNOW_AGG(I)  <=  0.0) THEN
           CTTEMP(I) = T(I)
         ENDIF

       ENDDO
!
! ======================================================================
!        Droplet settling
! ======================================================================
         If (l_droplet_settle) Then

! DEPENDS ON: lsp_settle
           Call lsp_settle(points, timestep, lsiter                     &
     &,                    q, qcl, T, droplet_flux                      &
     &,                    l_update_cf, l_seq, bland                    &
     &,                    cfliq, rho, rhor, corr2, lcrcp               &
     &,                    dhi, dhir                                    &
     &,                    Ntot_land, Ntot_sea                          &
     &,                    plset, plevpset                              &
     &                    )

         End if  ! l_droplet_settle

! ======================================================================
!        Sedimentation of ice, rain and graupel
! ======================================================================
         If (l_it_melting) Then
           ! Iterate the fall-out and the melting terms

! ======================================================================
!        TO TRY AND FIX 3D BUG: LSP_SUBGRID MOVED HERE - FRJW
! ======================================================================
! ======================================================================
!        Subgrid-scale set-up calculations and tidy-ups
! ======================================================================
! DEPENDS ON: lsp_subgrid
         Call lsp_subgrid(points                                        &
     &,                  q, qcf_cry, qcf_agg, qcf_tot, T                &
     &,                  qsl, qs                                        &
     &,                  q_ice, q_clear, q_ice_1, q_ice_2               &
     &,                  area_liq,area_mix,area_ice,area_clear          &
     &,                  area_ice_1, area_ice_2                         &
     &,                  rain_liq,rain_mix,rain_ice,rain_clear          &
     &,                  cf, cfliq, cfice, cficei                       &
     &,                  frac_ice_above                                 &
     &,                  cfkeep, cficekeep, rainfrac                    &
     &,                  lsrcp, rhcpt, l_update_cf                      &
     &                  )


! DEPENDS ON: lsp_it_fall_melt
           Call lsp_it_fall_melt(points, timestep                       &
     &,                  T, p, q, q_ice, qsl                            &
     &,                  qcf_cry, qcf_agg, frac_agg, qrain, qgraup      &
     &,                  snow_agg, snow_cry, rainrate, grauprate        &
     &,                  snowt_agg, snowt_cry, rainratet, graupratet    &
     &,                  vf_agg, vf_cry, vf_rain, vf_graup              &
     &,                  area_liq, area_mix, area_clear, area_ice       &
     &,                  cfice, cficei, frac_ice_above                  &
     &,                  cfkeep, cfliqkeep, cficekeep                   &
     &,                  rainfrac,rain_liq,rain_mix,rain_ice,rain_clear &
     &,                  rho, rhor, tcg, tcgi, tcgc, tcgci, tcgg, tcggi &
     &,                  corr,corr2,rocor, dhi, dhir, cx, constp        &
     &,                  ai, bi, aic, bic, lfrcp                        &
     &,                  l_update_cf, l_seq                             &
     &,                  l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup          &
     &,                  l_psd, l_crystals                              &
     &,                  pifall,psfall,prfall,pgfall,pimlt,psmlt,pgmlt  &
     &,                  lsiter, n_iterations                           &
     &,                  cf_transfer_diag, cff_transfer_diag            &
     &,                  rf_transfer_diag                               &
     &                          )

         Else  ! l_it_melting

! DEPENDS ON: lsp_fall
           Call lsp_fall(points, timestep                               &
     &,                  qcf_cry, qcf_agg, frac_agg, qrain, qgraup, T   &
     &,                  snow_agg, snow_cry, rainrate, grauprate        &
     &,                  snowt_agg, snowt_cry, rainratet, graupratet    &
     &,                  vf_agg, vf_cry, vf_rain, vf_graup              &
     &,                  area_clear, area_ice, cfice, cficei            &
     &,                  frac_ice_above, cfkeep, cfliqkeep, cficekeep   &
     &,                  rho, rhor, tcgi, tcgci                         &
     &,                  corr, dhi, dhir, cx, constp, ai, bi, aic, bic  &
     &,                  l_update_cf, l_seq                             &
     &,                  l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup          &
     &,                  l_psd, l_crystals                              &
     &,                  pifall, psfall, prfall, pgfall                 &
     &,                  lsiter                                         &
     &,                  cf_transfer_diag, cff_transfer_diag            &
     &                  )

! ======================================================================
!        TO TRY AND FIX 3D BUG: LSP_SUBGRID MOVED HERE - FRJW
! ======================================================================
! ======================================================================
!        Subgrid-scale set-up calculations and tidy-ups
! ======================================================================
! DEPENDS ON: lsp_subgrid
         Call lsp_subgrid(points                                        &
     &,                  q, qcf_cry, qcf_agg, qcf_tot, T                &
     &,                  qsl, qs                                        &
     &,                  q_ice, q_clear, q_ice_1, q_ice_2               &
     &,                  area_liq,area_mix,area_ice,area_clear          &
     &,                  area_ice_1, area_ice_2                         &
     &,                  rain_liq,rain_mix,rain_ice,rain_clear          &
     &,                  cf, cfliq, cfice, cficei                       &
     &,                  frac_ice_above                                 &
     &,                  cfkeep, cficekeep, rainfrac                    &
     &,                  lsrcp, rhcpt, l_update_cf                      &
     &                  )

         End If  ! l_it_melting

! ======================================================================
!        HOMOGENEOUS (PIFRW) AND HETEROGENEOUS NUCLEATION (PIPRM)
! ======================================================================
         If (l_crystals) then
           ! Call nucleation with the crystals ice category (qcf_cry)
! DEPENDS ON: lsp_nucleation
           Call lsp_nucleation(points, timestep                         &
     &,                  q, qcl, qcf_cry, T                             &
     &,                  qs, qsl                                        &
     &,                  cfliq, cfice                                   &
     &,                  area_liq, area_mix                             &
     &,                  cfkeep, cfliqkeep, cficekeep                   &
     &,                  rain_liq, rain_mix                             &
     &,                  rhor, lheat_correc_ice                         &
     &,                  lfrcp, lsrcp                                   &
     &,                  l_update_cf, l_seq                             &
     &,                  piprm, pifrw, lsiter                           &
     &,                  cf_transfer_diag, cfl_transfer_diag            &
     &,                  cff_transfer_diag                              &
     &                  )
         Else
           ! Call nucleation with the only ice category (qcf_agg)
! DEPENDS ON: lsp_nucleation
           Call lsp_nucleation(points, timestep                         &
     &,                  q, qcl, qcf_agg, T                             &
     &,                  qs, qsl                                        &
     &,                  cfliq, cfice                                   &
     &,                  area_liq, area_mix                             &
     &,                  cfkeep, cfliqkeep, cficekeep                   &
     &,                  rain_liq, rain_mix                             &
     &,                  rhor, lheat_correc_ice                         &
     &,                  lfrcp, lsrcp                                   &
     &,                  l_update_cf, l_seq                             &
     &,                  piprm, pifrw, lsiter                           &
     &,                  cf_transfer_diag, cfl_transfer_diag            &
     &,                  cff_transfer_diag                              &
     &                  )    
         End if  ! l_crystals 

! ======================================================================
!             DEPOSITION/SUBLIMATION OF ICE CRYSTALS (PIDEP)
! ======================================================================
         If (l_crystals) then

         Do i=1,points
           ! Calculate total ice content
           qcft(i) = qcf_cry(i) + qcf_agg(i)
         End do

! DEPENDS ON: lsp_deposition
         Call lsp_deposition(points, timestep                           &
     &,                  q, qcl, qcf_cry, qcft, T, p, frac_cry_dep      &
     &,                  q_ice_1, q_ice_2                               &
     &,                  area_ice_1, area_ice_2                         &
     &,                  esi, qs, qsl                                   &
     &,                  area_mix, cfliq, cfice, cficei                 &
     &,                  cfkeep, cfliqkeep, cficekeep                   &
     &,                  rho, tcgc, tcgci                               &
     &,                  corr2, rocor, lheat_correc_ice                 &
     &,                  cx, constp, ai, bi, aic, bic, lfrcp, lsrcp, 0  &
     &,                  l_update_cf, l_seq, .false.                    &
     &,                  pidep, lsiter                                  &
     &,                  cf_transfer_diag, cfl_transfer_diag            &
     &,                  cff_transfer_diag                              &
     &                  )

         End if  ! l_crystals
! ======================================================================
!           DEPOSITION/SUBLIMATION OF SNOW AGGREGATES (PSDEP)
! ======================================================================
         Do i=1,points
           ! Calculate total ice content
           qcft(i) = qcf_cry(i) + qcf_agg(i)
         End do

! DEPENDS ON: lsp_deposition
         Call lsp_deposition(points, timestep                           &
     &,                  q, qcl, qcf_agg, qcft, T, p, frac_agg_dep      &
     &,                  q_ice_1, q_ice_2                               &
     &,                  area_ice_1, area_ice_2                         &
     &,                  esi, qs, qsl                                   &
     &,                  area_mix, cfliq, cfice, cficei                 &
     &,                  cfkeep, cfliqkeep, cficekeep                   &
     &,                  rho, tcg, tcgi                                 &
     &,                  corr2, rocor, lheat_correc_ice                 &
     &,                  cx, constp, ai, bi, aic, bic, lfrcp, lsrcp, 1  &
     &,                  l_update_cf, l_seq, l_psd                      &
     &,                  psdep, lsiter                                  &
     &,                  cf_transfer_diag, cfl_transfer_diag            &
     &,                  cff_transfer_diag                              &
     &                  )

         If (l_mcr_qcf2) Then  ! Note that l_mcr_qcf2==true implies
                               ! that l_crystals==true. It is *intended*
                               ! in parametrization development that
                               ! l_mcr_qcf2==true implies l_psd==false.
! ======================================================================
!          AUTOCONVERSION OF ICE CRYSTALS TO AGGREGATES (PSAUT)
! ======================================================================
! DEPENDS ON: lsp_snow_autoc
           Call lsp_snow_autoc(points, timestep,                        &
                         qcf_cry, qcf_agg, T, CTTemp,                   &
                         m0, T_scaling, qcf0,                           &
                         l_seq,                                         &
                         psaut, lsiter                                  &
                        )

! ======================================================================
!           COLLECTION OF CRYSTALS BY AGGREGATES (PSACI)
! ======================================================================
             l_use_area=.false.
             l_no_t_check=.true.
             ice_type1=1  ! aggregates
             ice_type2=0  ! crystals
! DEPENDS ON: lsp_collection
             Call lsp_collection(points, timestep,                      &
                         qcf_agg, qcf_cry, T,                           &
                         area_mix, area_ice, cficei,                    &
!                        cf, cff,
                         rho, rhor, m0, tcg, tcgi, tcgc, tcgci,         &
                         corr, cx, constp, ai, bi, aic, bic,            &
                         ice_type1,ice_type2,                           &
                         l_update_cf, l_seq, .false.,                   &
                         l_use_area, l_no_t_check,                      &
                         psaci, lsiter                                  &
!                        cf_transfer_diag, cff_transfer_diag
                        )

           End if  ! l_mcr_qcf2

! ======================================================================
!              RIMING OF ICE CRYSTALS BY CLOUD WATER (PIACW)
! ======================================================================
         If (l_crystals) then
! DEPENDS ON: lsp_riming
         Call lsp_riming(points, timestep                               &
     &,                  qcl, qcf_cry, T                                &
     &,                  area_liq, area_mix, cfliq, cficei              &
!    &,                  cfkeep, cfliqkeep, cficekeep
     &,                  rho, m0, tcgc, tcgci, corr                     &
     &,                  cx, constp, ai, bi, aic, bic, lfrcp, 0         &
     &,                  l_update_cf, l_seq, .false.                    &
     &,                  piacw, lsiter                                  &
!    &,                  cf_transfer_diag, cfl_transfer_diag
!    &,                  cff_transfer_diag
     &                  )
         End if  ! l_crystals

! ======================================================================
!             RIMING OF SNOW AGGREGATES BY CLOUD WATER (PSACW)
! ======================================================================
! DEPENDS ON: lsp_riming
         Call lsp_riming(points, timestep                               &
     &,                  qcl, qcf_agg, T                                &
     &,                  area_liq, area_mix, cfliq, cficei              &
!    &,                  cfkeep, cfliqkeep, cficekeep
     &,                  rho, m0, tcg, tcgi, corr                       &
     &,                  cx, constp, ai, bi, aic, bic, lfrcp, 1         &
     &,                  l_update_cf, l_seq, l_psd                      &
     &,                  psacw, lsiter                                  &
!    &,                  cf_transfer_diag, cfl_transfer_diag
!    &,                  cff_transfer_diag
     &                  )

         If (l_mcr_qgraup) Then
! ======================================================================
!             AUTOCONVERSION OF SNOW TO GRAUPEL  (PGAUT)   
! ======================================================================

! DEPENDS ON: lsp_graup_autoc
           Call lsp_graup_autoc(points, timestep,                       &
                         qcf_agg, qgraup, T, rho,                       &
                         auto_graup_qcf_thresh, auto_graup_T_thresh,    &
                         auto_graup_coeff,                              &
                         l_seq,                                         &
                         psacw, psdep, pgaut, lsiter                    &
                        )

! ======================================================================
!              RIMING OF GRAUPEL BY CLOUD WATER (PGACW)
! ======================================================================
           ! Graupel is not included in ice cloud fraction
! DEPENDS ON: lsp_riming
           Call lsp_riming(points, timestep                             &
     &,                    qcl, qgraup, T                               &
     &,                    area_liq, area_mix, cfliq, cficei            &
!    &,                    cfkeep, cfliqkeep, cficekeep
     &,                    rho, m0, tcgg, tcggi, corr                   &
     &,                    cx, constp, ai, bi, aic, bic, lfrcp, 3       &
     &,                    .false., l_seq, .false.                      &
     &,                    pgacw, lsiter                                &
!    &,                    cf_transfer_diag, cfl_transfer_diag
!    &,                    cff_transfer_diag
     &                  )

! ======================================================================
!           COLLECTION OF SNOW BY GRAUPEL (PGACS)
! ======================================================================
           l_use_area=.false.
           l_no_t_check=.false.
           ice_type1=3  ! graupel
           ice_type2=1  ! aggregates
! DEPENDS ON: lsp_collection
           Call lsp_collection(points, timestep,                        &
                           qgraup, qcf_agg, T,                          &
                           area_mix, area_ice, cficei,                  &
!                          cf, cff,
                           rho, rhor, m0, tcgg, tcggi, tcg, tcgi,       &
                           corr, cx, constp, ai, bi, aic, bic,          &
                           ice_type1, ice_type2,                        &
                           l_update_cf, l_seq, l_psd,                   &
                           l_use_area, l_no_t_check,                    &
                           pgacs, lsiter                                &
!                          cf_transfer_diag, cff_transfer_diag
                        )

         End if  ! l_mcr_graup

! ======================================================================
!               COLLECTION OF RAIN BY ICE CRYSTALS (PIACR)
! ======================================================================
         If (l_crystals) then
! DEPENDS ON: lsp_capture
         Call lsp_capture(points, timestep                              &
     &,                   qcf_cry, qrain, T                             &
     &,                   area_liq, area_mix, area_ice, cficei          &
     &,                   rainfrac,rain_liq,rain_mix,rain_ice,rain_clear&
!    &,                   cfkeep, cficekeep
     &,                   rho, rhor, m0, tcgc, tcgci, corr, dhilsiterr  &
     &,                   cx, constp, ai, bi, aic, bic, lfrcp, 0        &
     &,                   l_update_cf, l_seq, .false.                   &
     &,                   l_mcr_qrain, piacr, lsiter                    &
!    &,                   cf_transfer_diag, cff_transfer_diag
     &,                   rf_transfer_diag                              &
     &                   )
         End if

! ======================================================================
!               COLLECTION OF RAIN BY SNOW AGGREGATES (PSACR)
! ======================================================================
! DEPENDS ON: lsp_capture
         Call lsp_capture(points, timestep                              &
     &,                   qcf_agg, qrain, T                             &
     &,                   area_liq, area_mix, area_ice, cficei          &
     &,                   rainfrac,rain_liq,rain_mix,rain_ice,rain_clear&
!    &,                   cfkeep, cficekeep
     &,                   rho, rhor, m0, tcg, tcgi, corr, dhilsiterr    &
     &,                   cx, constp, ai, bi, aic, bic, lfrcp, 1        &
     &,                   l_update_cf, l_seq, l_psd                     &
     &,                   l_mcr_qrain, psacr, lsiter                    &
!    &,                   cf_transfer_diag, cff_transfer_diag
     &,                   rf_transfer_diag                              &
     &                   )

! ======================================================================
!                 EVAPORATION OF MELTING ICE CRYSTALS (PIMLTEVP)
! ======================================================================
         If (l_crystals) then

         Do i=1,points
           ! Calculate total ice content
           qcft(i) = qcf_cry(i) + qcf_agg(i)
         End do

! DEPENDS ON: lsp_evap_snow
         Call lsp_evap_snow(points, timestep                            &
     &,                     q, q_ice, qcf_cry, qcft, T, p, esw, qsl     &
     &,                     area_ice, cficei, cfkeep, cficekeep         &
     &,                     rho, tcgc, tcgci                            &
     &,                     corr2, rocor, lheat_correc_liq              &
     &,                     cx, constp, ai, bi, aic, bic, lsrcp, 0      &
     &,                     l_update_cf, l_seq, .false., pimltevp,lsiter&
     &,                     cf_transfer_diag,cff_transfer_diag          &
     &                     )
         End if  ! l_crystals

! ======================================================================
!                 EVAPORATION OF MELTING SNOW AGGREGATES (PSMLTEVP)
! ======================================================================
         Do i=1,points
           ! Calculate total ice content
           qcft(i) = qcf_cry(i) + qcf_agg(i)
         End do

! DEPENDS ON: lsp_evap_snow
         Call lsp_evap_snow(points, timestep                            &
     &,                     q, q_ice, qcf_agg, qcft, T, p, esw, qsl     &
     &,                     area_ice, cficei, cfkeep, cficekeep         &
     &,                     rho, tcg, tcgi                              &
     &,                     corr2, rocor, lheat_correc_liq              &
     &,                     cx, constp,  ai, bi, aic, bic, lsrcp, 1     &
     &,                     l_update_cf, l_seq, l_psd, psmltevp, lsiter &
     &,                     cf_transfer_diag,cff_transfer_diag          &
     &                     )

         If (.not. l_it_melting) Then
         ! Perform the melting steps here

! ======================================================================
!                    MELTING OF ICE CRYSTALS (PIMLT)
! ======================================================================
         If (l_crystals) then

         Do i=1,points
           ! Calculate total ice content
           qcft(i) = qcf_cry(i) + qcf_agg(i)
         End do

! DEPENDS ON: lsp_melting
         Call lsp_melting(points, timestep                              &
     &,                   q, q_ice, qcf_cry, qcft, qrain, qsl, T, p     &
     &,                   area_liq, area_mix, area_ice                  &
     &,                   cfice, cficei, frac_ice_above                 &
     &,                   rainfrac, rain_liq, rain_mix                  &
     &,                   rain_ice, rain_clear, cfkeep, cficekeep       &
     &,                   rho, rhor, m0, tcg, tcgi, corr2, rocor        &
     &,                   cx, constp, ai, bi, aic, bic, lfrcp, 0        &
     &,                   l_update_cf, l_seq, .false.                   &
     &,                   pimlt, lsiter                                 &
     &,                   cf_transfer_diag, cff_transfer_diag           &
     &,                   rf_transfer_diag                              &
     &                   )
         End if  ! l_crystals

! ======================================================================
!                    MELTING OF SNOW AGGREGATES (PSMLT)
! ======================================================================
         Do i=1,points
           ! Calculate total ice content
           qcft(i) = qcf_cry(i) + qcf_agg(i)
         End do

! DEPENDS ON: lsp_melting
         Call lsp_melting(points, timestep                              &
     &,                   q, q_ice, qcf_agg, qcft, qrain, qsl, T, p     &
     &,                   area_liq, area_mix, area_ice                  &
     &,                   cfice, cficei, frac_ice_above                 &
     &,                   rainfrac, rain_liq, rain_mix                  &
     &,                   rain_ice, rain_clear, cfkeep, cficekeep       &
     &,                   rho, rhor, m0, tcgc, tcgci, corr2, rocor      &
     &,                   cx, constp, ai, bi, aic, bic, lfrcp, 1        &
     &,                   l_update_cf, l_seq, l_psd                     &
     &,                   psmlt, lsiter                                 &
     &,                   cf_transfer_diag, cff_transfer_diag           &
     &,                   rf_transfer_diag                              &
     &                   )

! ======================================================================
!                    MELTING OF GRAUPEL (PGMLT)
! ======================================================================
         If (L_mcr_qgraup) Then
         ! Graupel does not update cloud fractions so there is no need
         ! to update qcft (it is not used)
! DEPENDS ON: lsp_melting
           Call lsp_melting(points, timestep                            &
     &,                     q, q_ice, qgraup, qcft, qrain, qsl, T, p    &
     &,                     area_liq, area_mix, area_ice                &
     &,                     cfice, cficei, frac_ice_above               &
     &,                     rainfrac, rain_liq, rain_mix                &
     &,                     rain_ice, rain_clear, cfkeep, cficekeep     &
     &,                     rho, rhor, m0, tcgg, tcggi, corr2, rocor    &
     &,                     cx, constp, ai, bi, aic, bic, lfrcp, 3      &
     &,                     .false., l_seq, .false.                     &
     &,                     pgmlt, lsiter                               &
     &,                     cf_transfer_diag, cff_transfer_diag         &
     &,                     rf_transfer_diag                            &
     &                     )
         End If  ! L_mcr_qgraup

         End If  ! .not. l_it_melting

! ======================================================================
!                   EVAPORATION OF RAINDROPS (PREVP)
! ======================================================================
! DEPENDS ON: lsp_evap
         Call lsp_evap(points, timestep, p, q, qrain, T                 &
     &,                q_ice, q_clear                                   &
     &,                area_liq, area_mix, area_ice, area_clear         &
     &,                rainfrac, rain_liq, rain_mix                     &
     &,                rain_ice, rain_clear                             &
     &,                rho, corr, corr2, rocor                          &
     &,                dhilsiterr, cx, constp, lcrcp, lheat_correc_liq  &
     &,                qsl, esw, l_seq, l_mcr_qrain                     &
     &,                prevp, rf_transfer_diag, lsiter                  &
     &                )

! ======================================================================
!             ACCRETION OF CLOUD DROPLETS ON RAINDROPS (PRACW)
! ======================================================================
! DEPENDS ON: lsp_accretion
         Call lsp_accretion(points, timestep, qcl, qrain                &
!    &,                     area_liq, area_mix, area_ice
     &,                     cfliq, rainfrac, rain_liq, rain_mix         &
!    &,                     rain_ice, rain_clear, cfkeep, cfliqkeep
     &,                     rho, corr, dhilsiterr, cx, constp           &
     &,                     l_update_cf, l_seq, l_mcr_qrain             &
     &,                     pracw, lsiter                               &
!    &,                     cf_transfer_diag, cfl_transfer_diag
!    &,                     rf_transfer_diag
     &                     )

! ======================================================================
!              AUTOCONVERSION OF CLOUD LIQUID TO RAIN (PRAUT)
! ======================================================================
! DEPENDS ON: lsp_autoc
           Call lsp_autoc(points, timestep, qcl, qrain, T, p            &
     &,                  cfliq, rhcpt                                   &
!    &,                  cfkeep, cfliqkeep
     &,                  area_liq, area_mix, area_ice, rainfrac         &
     &,                  rain_liq, rain_mix, rain_ice, rain_clear       &
     &,                  rho, rhor, corr2, ec_auto, lcrcp               &
     &,                  l_update_cf, l_seq, l_auto_debias, l_murk      &
     &,                  l_use_sulphate_autoconv, l_seasalt_ccn         &
     &,                  l_biomass_ccn, l_ocff_ccn, l_biogenic_ccn      &
     &,                  l_autoconv_murk                                &
     &,                  l_autoc_3b, l_autolim_3b, l_mixing_ratio       &
     &,                  lsiter, praut, rf_transfer_diag                &
!    &,                  cf_transfer_diag, cfl_transfer_diag
     &,                  aerosol, SO4_ait, SO4_acc, SO4_dis             &
     &,                  sea_salt_film, sea_salt_jet, bmass_agd         &
     &,                  bmass_cld, ocff_agd, ocff_cld, biogenic        &
     &,                  snow_depth, land_fract                         &
     &,                  Ntot_land, Ntot_sea                            &
     &                  )

! ----------------------------------------------------------------------
!  End loop over iterations.
! ----------------------------------------------------------------------

       END DO ! Iters_do1
!
! ======================================================================
!
!
!                  COPY ICE/SNOW VARIABLES AND FLUXES
!                        TO OUTPUT VARIABLES
!
!
! ======================================================================
! Copy contents of ice/snow variables and fluxes to output variables
! to fall into next layer down
! ----------------------------------------------------------------------

       DO I=1,POINTS

        IF (L_mcr_qcf2) THEN ! two ice prognostics
          QCF(I)  = QCF_AGG(I)
          QCF2(I) = QCF_CRY(I)
          SNOW_CRY(I) = SNOWT_CRY(I)
          SNOW_AGG(I) = SNOWT_AGG(I)
        ELSE ! only one ice prognostic, put all snow in to snow_agg
          QCF(I) = QCF_CRY(I) + QCF_AGG(I)
          SNOW_CRY(I) = 0.0  ! Redundant variable
          SNOW_AGG(I) = SNOWT_CRY(I) + SNOWT_AGG(I)
        ENDIF ! on L_mcr_qcf2

        IF (L_mcr_qrain)  RAINRATE(I) = RAINRATET(I)
        IF (L_mcr_qgraup) GRAUPRATE(I) = GRAUPRATET(I)

      END DO ! Points

! ======================================================================
!              NUMERICAL TIDYING UP OF SMALL VALUES
! ======================================================================
! DEPENDS ON: lsp_tidy
      Call lsp_tidy(points, timestep, lsiter, q,qcl,qcf,qcf2, qrain, T  &
     &,            area_liq, area_mix, area_ice                         &
     &,            cfice, cficei, cfkeep, cfliqkeep, cficekeep          &
     &,            rainfrac, rain_liq, rain_mix, rain_ice, rain_clear   &
     &,            q_ice, qs, qsl, snow_agg, snow_cry                   &
     &,            rho, rhor, p                                         &
     &,            cttemp,dhi,dhilsiterr,frac_ice_above                 &
     &,            lcrcp, lfrcp, lsrcp                                  &
     &,            l_update_cf, l_seq, l_mcr_qcf2, l_it_melting         &
     &,            psdep, pidep, psmlt, pimlt, prevp                    &
     &,            cf_transfer_diag,cfl_transfer_diag                   &
     &,            cff_transfer_diag, rf_transfer_diag                  &
     &             )

      DO I=1,POINTS
        !------------------------------------------------
        ! Now update fraction of ice in layer above
        ! for next layer down
        !------------------------------------------------
        FRAC_ICE_ABOVE(I)=CFICEKEEP(I)

      END DO ! Points


! ======================================================================
!       IF DIAGNOSTIC RAIN, CONVERT MASS (kg/kg) TO FLUX (kg/m2/s)
! ======================================================================
      IF (.NOT. L_mcr_qrain) THEN

        DO I=1,POINTS
          If (l_non_hydrostatic) Then
            If (l_mixing_ratio) Then
              ! Use mixing ratio formulation
              rainrate(i) = qrain(i)*rhodz_dry(i)/timestepfixed
            Else
              ! Use specific humidity formulation
              rainrate(i) = qrain(i)*rhodz_moist(i)/timestepfixed
            End if  ! l_mixing_ratio
          Else
            ! Use hydrostatic formulation
            rainrate(i) = qrain(i)*rhodz(i)/timestepfixed
          End if  ! l_non_hydrostatic
        END DO

      ENDIF ! on prognostic rain mixing ratio

! ----------------------------------------------------------------------
!   End of the LSP_ICE subroutine
! ----------------------------------------------------------------------
      RETURN
      END SUBROUTINE LSP_ICE
