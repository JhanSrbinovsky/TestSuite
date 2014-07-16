
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Evaporation of melting snow
! Subroutine Interface:
      SUBROUTINE LSP_EVAP_SNOW(                                         &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, q, q_ice, qcf, qcft, T, p                                       &
                                          ! Water contents, temp, pres
     &, esw, qsl                                                        &
                                          ! Saturated quantities
     &, area_ice, cficei                                                &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
     &, cf, cff                                                         &
                                          ! Current cloud fractions for
                                          ! updating
     &, rho, tcg, tcgi                                                  &
                                          ! Parametrization information
     &, corr2, rocor, lheat_correc_liq                                  &
     &, cx, constp, ai, bi, aic, bic, lsrcp, ice_type                   &
                                          ! Microphysical information
     &, l_update_cf, l_seq, l_psd                                       &
                                          ! Code options
     &, ptransfer, iterations                                           &
                                          ! Mass transfer diagnostic
     &, cftransfer,cfftransfer                                          &
                                          ! Cloud transfer diagnostics
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of sublimation of melting
!   snow
!
! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles in a distribution of vapour
!
! Current Owner of Code: Jonathan Wilkinson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!
! Source for vapour. Sink for ice.
!
! Subroutine Arguments
!
! Obtain the size for cx and constp and define ice_type_offset
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
!
      Integer, Intent(In) ::                                            &
     &  points                                                          &
                          ! Number of points to calculate
     &, iterations                                                      &
                          ! Number of microphysics iterations
                          ! (for rate diagnostic)
     &, ice_type          ! Type of ice (0 - crystals, 1 - aggregates)
!
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, p(points)                                                       &
                          ! Air pressure / N m-2
     &, q_ice(points)                                                   &
                          ! Vapour content in ice partition / kg kg-1
     &, esw(points)                                                     &
                          ! Saturated vapour pres. over liquid / N m-2
     &, qsl(points)                                                     &
                          ! Saturated humidity wrt liquid / kg kg-1
     &, area_ice(points)                                                &
                          ! Fraction of gridbox with ice but no liquid
     &, cficei(points)                                                  &
                          ! 1/Fraction of gridbox with ice cloud
     &, rho(points)                                                     &
                          ! Air density / kg m-3
!    &, m0                ! Seed ice water content / kg kg-1 (c_lspmic)
     &, tcg(points)                                                     & 
                          ! T dependent function in ice size dist'n
     &, tcgi(points)                                                    &
                          ! 1/tcg (no units)
     &, corr2(points)                                                   &
                          ! Temperature correction factor (no units)
     &, rocor(points)                                                   &
                          ! Combined fall and corr2 correction factor
     &, lheat_correc_liq(points)                                        &
                                 ! Liquid latent heat correction factor
     &, ai, bi, aic, bic                                                &
                          ! Ice mass size relationship. m(D) = AI D^BI
     &, lsrcp             ! Latent heat of sublimation
                          ! /heat capacity of air (cP) / K
!
      Real, Intent(InOut) ::                                            &
     &  q(points)                                                       &
                          ! Vapour content / kg kg-1
     &, qcf(points)                                                     &
                          ! Ice water content in ice category to be
!                           updated    / kg kg-1
     &, qcft(points)                                                    &
                          ! Ice water in all ice categories
!                           (for cloud fraction calculations)
     &, T(points)                                                       &
                          ! Temperature / K
     &, cf(points)                                                      &
                          ! Current cloud fraction
     &, cff(points)                                                     &
                          ! Current ice cloud fraction
     &, ptransfer(points)  ! Mass deposited in this timestep / kg kg-1
!
      Real, Intent(InOut) ::                                            &
     &  cftransfer(points)                                              &
                           ! Cloud fraction increment this tstep
     &, cfftransfer(points)! Ice cloud fraction inc this tstep
!
      Logical, Intent(In) ::                                            &
     &  l_update_cf                                                     &
                          ! Update cloud fractions
     &, l_seq                                                           &
                          ! Carry out sequential updating
     &, l_psd
                          ! Use generic ice particle size distribution
!
! Local Variables
!
      Integer                                                           &
     &  i                                                               &
                          ! Loop counter for points
     &, cry_offset        ! Index offset for ice crystals

      Real                                                              &
     &  pr02(points)                                                    &
                          ! Temporary in calculation of PSD slopes
     &, pr04(points)                                                    &
                          ! Temporary in calculation of PSD slopes
     &, dpr(points)                                                     &
                          ! Temporary in calculating ice transfers
     &, deltacf(points)                                                 &
                          ! Change in cf across timestep
     &, tempw(points)                                                   &
                          ! Available subsaturation for evaporation
     &, m_1(points)                                                     &
                          ! 1st moment of the generic ice size distribtn
     &, m_0p5_dp3(points)
                          ! 1+(di+1)/2 moment of the generic ice PSD

      !-----------------------------------------------
      ! Select appropriate ice parametrization (see c_lspsiz)
      !-----------------------------------------------
      cry_offset=ice_type*ice_type_offset

      ! Use the generic ice particle size distribution
      ! Calculate the 1st (cx(84)) moment and the 1+0.5(di+1) (cx(85))
      ! moment of the ice particle size distribution
      If (l_psd) then

! DEPENDS ON: lsp_moments
            Call lsp_moments(points,rho,T,qcf,cficei,                   &
     &                       ai,bi,cx(84),m_1)
! DEPENDS ON: lsp_moments
            Call lsp_moments(points,rho,T,qcf,cficei,                   &
     &                       ai,bi,cx(85),m_0p5_dp3)

      End if

      Do i=1, points

        If (qcf(i) >  m0 .and. T(i) >  zerodegc) Then

          !-----------------------------------------------
          ! Diffusional growth parameters
          !-----------------------------------------------
          pr04(i) = ((apb4-apb5*T(i))*esw(i)+apb6*p(i)*T(i)**3)

          !-----------------------------------------------
          ! Calculate transfer rates
          !-----------------------------------------------
          If (l_psd) then
            ! Use generic particle size distribution
            ! constp(83) = 2 pi axial_ratio_correction
            ! constp(84) = ventilation coefficient 1
            ! constp(85) = ventilation coefficient 2
            !        * Sc^(1/3)*ci^0.5/viscosity0^0.5
            dpr(i) = constp(83) * timestep * T(i)**2 * esw(i)           &
    &                           * (constp(84)*m_1(i)*corr2(i)           &
    &                           +  constp(85)*rocor(i)*m_0p5_dp3(i))    &
    &                           / (qsl(i) * rho(i) * pr04(i))

          Else
            ! Use particle size distribution based on intercepts
            pr02(i)=rho(i)*qcf(i)*cficei(i)*constp(5+cry_offset)*tcgi(i)
            dpr(i)=tcg(i)*constp(6+cry_offset)*T(i)**2*esw(i)*timestep* &
     &      (constp(7+cry_offset)*corr2(i)*pr02(i)**cx(4+cry_offset)    &
     &       +constp(8+cry_offset)*rocor(i)*pr02(i)**cx(5+cry_offset))  &
     &       /(qsl(i)*rho(i)*pr04(i))

          End if  ! l_psd

          !-----------------------------------------------
          ! Limit transfers
          !-----------------------------------------------
          ! tempw is the subsaturation that is available
          tempw(i) = area_ice(i) * (qsl(i) - q_ice(i))
          dpr(i) = dpr(i) * tempw(i)
          dpr(i) = max(min(dpr(i),tempw(i)*lheat_correc_liq(i)),0.0)

          ! Limit on the amount of ice available
          dpr(i) = min(dpr(i),qcf(i))

          !-----------------------------------------------
          ! Store process rate
          !-----------------------------------------------
          ptransfer(i) = ptransfer(i) + dpr(i)/(timestep*iterations)

          !-----------------------------------------------
          ! Store and update ice cloud fraction and total fractions
          !-----------------------------------------------
          If (l_update_cf) Then
            deltacf(i) = (-cff(i)*dpr(i)/qcft(i))
            cftransfer(i) =cftransfer(i) +deltacf(i)                    &
     &                       /(timestep*iterations)
            cfftransfer(i)=cfftransfer(i)+deltacf(i)                    &
     &                       /(timestep*iterations)

            If (l_seq) Then
              cf(i)  = cf(i)  + deltacf(i)
              cff(i) = cff(i) + deltacf(i)
            End if  ! l_seq

          End if  ! l_update_cf

          !-----------------------------------------------
          ! Update values of ice and vapour
          !-----------------------------------------------
          If (l_seq) Then
            qcf(i) = qcf(i) - dpr(i)
            q(i)   = q(i)   + dpr(i)
            T(i)   = T(i)   - dpr(i)*lsrcp
          End if  ! l_seq

        End If  ! qcf gt m0 etc.

      End do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_EVAP_SNOW
