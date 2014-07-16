#if defined(A04_3B)
! --------------------------COMDECK C_LSPMIC----------------------------
! SPECIFIES MICROPHYSICAL PARAMETERS FOR AUTOCONVERSION, HALLETT MOSSOP
! PROCESS, ICE NUCLEATION. ALSO SPECIFIES NUMBER OF ITERATIONS OF
! THE MICROPHYSICS.
! ----------------------------------------------------------------------
!
! History:
!
! Version    Date     Comment
! -------    ----     -------
!   5.3    24/10/01   Removed hardwired switch for second indirect
!                     effect: now set via user interface.
!                                                        Andy Jones
!   6.1    01/08/04   Include graupel autoconversion params. R.M.Forbes
!   6.2    17/11/05   Remove variables that are now in UMUI. D. Wilson
!   6.2    03/02/06   Include droplet settling logical. Damian Wilson

! ----------------------------------------------------------------------
!      AUTOCONVERSION TERMS
! ----------------------------------------------------------------------

! WARNING. BE AWARE THAT DROPLET CONCENTRATION IS ALSO DEFINED IN THE
! RADIATION SCHEME. ARE YOU HAPPY THAT THE VALUES ARE CONSISTENT?

      ! Inhomogeneity factor for autoconversion rate
      REAL,PARAMETER:: INHOMOG_RATE=1.0

      ! Inhomogeneity factor for autoconversion limit
      REAL,PARAMETER:: INHOMOG_LIM=1.0

      ! Collision collection coefficient
!     REAL,PARAMETER:: EC_AUTO is set in UMUI

! The compilation is being picky about using non integer powers in
! parameter statements. The best I can do at the moment is to
! directly define 1/cubed roots of droplet concentrations.
! Just be careful you remember to change all of these.
      ! Droplet concentration over land
!     REAL,PARAMETER:: N_DROP_LAND is set in UMUI

      ! (N_DROP_LAND)^(-1/3)
!     REAL,PARAMETER:: N_DROP_LAND_CR is set in UMUI

      ! Droplet concentration over sea
!     REAL,PARAMETER:: N_DROP_SEA is set in UMUI

      ! (N_DROP_SEA)^(-1/3)
!     REAL,PARAMETER:: N_DROP_SEA_CR is set in UMUI

      ! Threshold droplet radius for autoconversion
      REAL,PARAMETER:: R_THRESH=7.0E-6

! ----------------------------------------------------------------------
!     ITERATIONS OF MICROPHYSICS
! ----------------------------------------------------------------------

      INTEGER,PARAMETER:: LSITER=1 !   Number of iterations in microphysi
      ! Advise 1 iteration for every 10 minutes or less of timestep.

      INTEGER,PARAMETER:: ADV_TYPE=2 ! Vertical advection method.
      ! ADV_TYPE=1: Original formulation
      ! ADV_TYPE=2: Revised formulation has better fall through layers

! ----------------------------------------------------------------------
!     NUCLEATION OF ICE
! ----------------------------------------------------------------------

! Note that the assimilation scheme uses temperature thresholds
! in its calculation of qsat.

      ! Nucleation mass
      REAL,PARAMETER:: M0=1.0E-12
      ! Maximum Temp for ice nuclei nucleation (deg C)
      REAL,PARAMETER:: TNUC=-10.0
      ! Maximum Temp for homogenous nucleation (deg C)
      REAL,PARAMETER:: THOMO=-40.0

! ----------------------------------------------------------------------
!     HALLETT MOSSOP PROCESS
! ----------------------------------------------------------------------

      ! Min temp for production of Hallett Mossop splinters (deg C)
      REAL,PARAMETER:: HM_T_MIN=-8.0


      ! Max temp for production of Hallett Mossop splinters (deg C)
      ! Switch off Hallett Mossop in this version but allow functionality
      REAL,PARAMETER:: HM_T_MAX=-273.0
!      REAL,PARAMETER:: HM_T_MAX=-3.0

      ! Residence distance for Hallett Mossop splinters (1/deg C)
      REAL,PARAMETER:: HM_DECAY=1.0/7.0

      ! Reciprocal of scaling liquid water content for HM process
      REAL,PARAMETER:: HM_RQCL=1.0/0.1E-3

      ! Minimum allowed QCF after microphysics
      REAL,PARAMETER:: QCFMIN=1.0E-8

! C_LSPMIC end
#endif
#if defined(A04_3C) || defined(A04_3D)
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

#endif
