
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Evaporation of rain
! Subroutine Interface:
      SUBROUTINE LSP_EVAP(                                              &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, p, q, qrain, T, q_ice, q_clear                                  &
                                          ! Water contents and temp
     &, area_liq, area_mix                                              &
                                          ! Cloud fraction partitions
     &, area_ice, area_clear                                            &
     &, rainfrac, rain_liq, rain_mix                                    &
                                          ! Rain fractions for updating
     &, rain_ice, rain_clear                                            &
     &, rho, corr, corr2,rocor,dhilsiterr                               &
                                          ! Parametrization information
     &, cx, constp                                                      &
                                          ! Microphysical information
     &, lcrcp, lheat_correc_liq, qsl, esw                               &
     &, l_seq, l_mcr_qrain                                              &
                                          ! Code options
     &, ptransfer, rftransfer, iterations                               &
                                          ! Mass transfer diagnostic
     &  )
!
      Implicit None
!
! Purpose:
!   Update rain and water vapour due to evaporation
!
! Method:
!   Integrate evaporation rates of a single raindrop over the
!   raindrop size distribution
!
! Current Owner of Code: Jonathan Wilkinson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!
!
! Subroutine Arguments
!
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
      Integer, Intent(In) ::                                            &
     &  points                                                          &
                          ! Number of points to calculate
     &, iterations        ! Number of microphysics iterations
                          ! (for rate diagnostic)

      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, p(points)                                                       &
                          ! Air pressure / N m-2
     &, q_ice(points)                                                   &
                          ! Local vapour in ice partition / kg kg-1
     &, q_clear(points)                                                 &
                          ! Local vapour in clear partition / kg kg-1
     &, area_liq(points)                                                &
                          ! Fraction of gridbox with liquid but no ice
     &, area_mix(points)                                                &
                          ! Fraction of gridbox with liquid and ice
     &, area_ice(points)                                                &
                          ! Fraction of gridbox with liquid and ice
     &, area_clear(points)                                              &
                          ! Fraction of gridbox with no condensate
     &, lcrcp                                                           &
                          ! Latent heat of condensation/cP / K
     &, dhilsiterr(points)                                              &
                          ! layer thickness/timestep / m s-1
     &, corr(points)                                                    &
                          ! Fall speed correction factor due to
                          !   air density changes (no units)
     &, corr2(points)                                                   &
                          ! Air diffusivity correction (no units)
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, lheat_correc_liq(points)                                        &
                                  ! Correction of subsaturation due
                                  ! to latent heat
     &, qsl(points)                                                     &
                          ! Saturated spec humidity wrt liquid / kg kg-1
     &, esw(points)                                                     &
                          ! Saturated vapour pressure wrt liquid / N m-2
     &, rocor(points)     ! Air density and viscosity correction factor
                          !   (no units)

      Real, Intent(InOut) ::                                            &
     &  qrain(points)                                                   &
                          ! Rain water content / kg kg-1
     &, q(points)                                                       &
                          ! Vapour content / kg kg-1
     &, T(points)                                                       &
                          ! Temperature / K
     &, ptransfer(points)                                               &
                          ! Evaporation rate / kg kg-1 s-1
     &, rftransfer(points)                                              &
                          ! Rate of change of rain fraction / s-1
     &, rainfrac(points)                                                &
                          ! Rain fraction (no units)
     &, rain_liq(points)                                                &
                          ! Overlap of rain with liquid (no units)
     &, rain_mix(points)                                                &
                          ! Overlap of rain with mixed phase region
     &, rain_ice(points)                                                &
                          ! Overlap of rain with ice
     &, rain_clear(points)! Overlap of rain with clear sky

!     Real, Intent(Out) ::

      Logical, Intent(In) ::                                            &
     &  l_seq                                                           &
                          ! Sequential updating
     &, l_mcr_qrain       ! Rain is a prognostic variable
!
! Local Variables
!
      Integer                                                           &
     &  i                 ! Loop counter

      Real                                                              &
     &  dpr(points)                                                     &
                          ! Amount of mass evaporated / kg kg-1
     &, pr04(points)                                                    &
                          ! Temporary in evaporation rate calc.
     &, lamr1(points)                                                   &
                          ! Powers of drop size distribution slope
     &, lamr2(points)                                                   &
                          !   lambda.
     &, temp7(points)                                                   &
                          ! Subsaturation in gridbox / kg kg-1
     &, rainfracnew(points) ! Updated rain fraction

      Do i=1,points
        !-----------------------------------------------
        ! If there is only a small amount of rain present, evaporate
        ! it completely to avoid numerical problems.
        !-----------------------------------------------
        If (qrain(i)  <   qcfmin) Then

          ! Evaporate all this rain
          dpr(i) = qrain(i)
          If (l_seq) Then
            T(i)   = T(i) - lcrcp * dpr(i)
            q(i)   = q(i) + dpr(i)
            qrain(i) = 0.0
          End if

          ! Store evaporation rate
          ptransfer(i) = ptransfer(i) + dpr(i)/(timestep*iterations)

          ! Update rain fractions
          rftransfer(i)  = rftransfer(i)                                &
     &                   -rainfrac(i)/(timestep*iterations)
          If (l_seq) Then
            rainfrac(i)  = 0.0
            rain_liq(i)  = 0.0
            rain_mix(i)  = 0.0
            rain_ice(i)  = 0.0
            rain_clear(i)= 0.0
          End if  ! l_seq

        End if  ! qrain lt qcfmin

        If (qrain(i)  >   0.0 .and. rainfrac(i)  >   0.0) Then
          !-----------------------------------------------
          ! Calculate evaporation parameters
          !-----------------------------------------------
          pr04(i) = ((apb4-apb5*T(i))*esw(i)+apb6*p(i)*T(i)**3)

          !-----------------------------------------------
          ! Calculate powers of drop size distribution slopes
          !-----------------------------------------------
          If (L_mcr_qrain) Then
            ! Rain is a mixing ratio - use fall speeds as parametrized
            lamr1(i) = qrain(i) * constp(50) * rho(i) / (rainfrac(i))
            lamr2(i) = lamr1(i) ** (cx(47)*cx(52))
            lamr1(i) = lamr1(i) ** (cx(49)*cx(52))
          Else
            ! Rain is a diagnostic quantity
            ! - use fall speeds to allow sedimentation into next layer
            lamr1(i) = qrain(i) * rho(i) * dhilsiterr(i) /              &
     &                 (constp(42) * corr(i) * rainfrac(i))
            lamr2(i) = lamr1(i) ** (cx(47)*cx(48))
            lamr1(i) = lamr1(i) ** (cx(49)*cx(48))
          End If  ! l_mcr_qrain

          !-----------------------------------------------
          ! Calculate evaporation rate
          !-----------------------------------------------
          dpr(i) = constp(46) * T(i)**2 * esw(i) * timestep
          dpr(i) = dpr(i) * ( (constp(47) * corr2(i) * lamr2(i))        &
     &               + (constp(48) * rocor(i) * lamr1(i)) )

          !-----------------------------------------------
          ! Calculate transfers
          !-----------------------------------------------
          ! Calculate gridbox mean supersaturation
          temp7(i) = (q_ice(i) - qsl(i)) * rain_ice(i)                  &
     &           +(q_clear(i) - qsl(i)) * rain_clear(i)

          ! Limit on the gridbox mean supersaturation
          dpr(i) = dpr(i) * max(-temp7(i)*lheat_correc_liq(i),0.0)      &
     &           /(qsl(i) * rho(i) * pr04(i) + dpr(i))

          ! Limit on the amount of rain available
          dpr(i) = min( dpr(i), qrain(i) *                              &
     &            (rain_ice(i) + rain_clear(i)) / rainfrac(i))

          !-----------------------------------------------
          ! Store process rate (kg/kg/s)
          !-----------------------------------------------
          ptransfer(i) = ptransfer(i) + dpr(i)/(timestep*iterations)

          !-----------------------------------------------
          ! Update values of rain, vapour and temperature
          !-----------------------------------------------
          If (l_seq) Then
            qrain(i) = qrain(i) - dpr(i)
            q(i)     = q(i)     + dpr(i)
            T(i)     = T(i)     - dpr(i) * LCRCP
          End If

          !-----------------------------------------------
          ! Calculate new rain fraction
          !-----------------------------------------------
!         These are commented out to ensure that rain fraction converges
!         to a 0 or 1 behaviour when RHcrit tends to 1.
!          rainfracnew(i) = rainfrac(i)*qrain(i)/(qrain(i)+dpr(i))
!          rftransfer(i)  = rftransfer(i)
!     &       + (rainfracnew(i) - rainfrac(i)) / (timestep*iterations)
!
          !-----------------------------------------------
          ! Update rain fractions
          !-----------------------------------------------
!         These are commented out to ensure that rain fraction converges
!         to a 0 or 1 behaviour when RHcrit tends to 1.
!          If (l_seq) Then
!            rainfrac(i)  = rainfracnew(i)
!            rain_liq(i)  = min(area_liq(i) , rainfrac(i))
!            rain_mix(i)  = min(area_mix(i) , rainfrac(i)-rain_liq(i))
!            rain_ice(i)  =
!     &          min(area_ice(i) , rainfrac(i)-rain_liq(i)-rain_mix(i))
!            rain_clear(i)= rainfrac(i)-rain_liq(i)
!     &                    -rain_mix(i)-rain_ice(i)
!          End If  ! l_seq

        End If  !  qrain gt 0 etc.

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_EVAP
