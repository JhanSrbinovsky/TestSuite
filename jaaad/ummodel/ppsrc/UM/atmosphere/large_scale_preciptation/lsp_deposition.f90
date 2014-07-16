
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Deposition of ice particles
! Subroutine Interface:
      SUBROUTINE LSP_DEPOSITION(                                        &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, q, qcl, qcf, qcft, T, p, frac_in_category                       &
                                          ! Water contents, temp, pres
     &, q_ice_1, q_ice_2                                                & 
                                          ! Subgrid-scale water c'tents
     &, area_ice_1, area_ice_2                                          &
                                          ! Subgrid-scale areas
     &, esi, qs, qsl                                                    &
                                          ! Saturated quantities
     &, area_mix, cfliq, cfice, cficei                                  &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
     &, cf, cfl, cff                                                    &
                                          ! Current cloud fractions for
                                          ! updating
     &, rho, tcg, tcgi                                                  &
                                          ! Parametrization information
     &, corr2, rocor, lheat_correc_ice                                  &
     &, cx, constp, ai, bi, aic, bic, lfrcp, lsrcp, ice_type            &
                                          ! Microphysical information
     &, l_update_cf, l_seq, l_psd                                       &
                                          ! Code options
     &, ptransfer, iterations                                           &
                                          ! Mass transfer diagnostic
     &, cftransfer,cfltransfer,cfftransfer                              &
                                          ! Cloud transfer diagnostics
     &  )


!dhb599 20110601: new cloud fraction increment scaling factor:
      use auscom_cpl_data_mod, only : xfactor

!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of ice deposition and
!   sublimation
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
! Source for ice. Sink for liquid water and vapour.
! Hallett Mossop provess enhances growth when turned on -
! Check value of HM_T_MAX in c_lspmic.
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
     &, esi(points)                                                     &
                          ! Saturated vapour pressure over ice / N m-2
     &, qs(points)                                                      &
                          ! Saturated humidity wrt ice / kg kg-1
     &, qsl(points)                                                     &
                          ! Saturated humidity wrt liquid / kg kg-1
     &, q_ice_1(points)                                                 &
                          ! Mean vapour in ice only region / kg kg-1
     &, q_ice_2(points)                                                 &
                          ! Mean vapour in clear region / kg kg-1
     &, area_ice_1(points)                                              &
                          ! Ice only area that is growing by deposition
     &, area_ice_2(points)                                              &
                          ! Ice only area that is subliming
     &, area_mix(points)                                                &
                          ! Fraction of gridbox with mixed phase cloud
     &, cfliq(points)                                                   &
                          ! Fraction of gridbox with liquid cloud
     &, cfice(points)                                                   &
                          ! Fraction of gridbox with ice cloud
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
     &, lheat_correc_ice(points)                                        &
                          ! Ice latent heat correction factor
     &, ai, bi, aic, bic                                                &
                          ! Ice mass size relationship. m(D) = AI D^BI
     &, frac_in_category(points)                                        &
                          ! Fraction of the ice that is in this category
     &, lfrcp                                                           &
                          ! Latent heat of fusion
                          ! /heat capacity of air (cP) / K
     &, lsrcp             ! Latent heat of sublimation/cP / K

!
      Real, Intent(InOut) ::                                            &
     &  q(points)                                                       &
                          ! Vapour content / kg kg-1
     &, qcl(points)                                                     &
                          ! Liquid water content / kg kg-1
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
     &, cfl(points)                                                     &
                          ! Current liquid cloud fraction
     &, cff(points)                                                     &
                          ! Current ice cloud fraction
     &, ptransfer(points)  ! Mass deposited in this timestep / kg kg-1
!
      Real, Intent(InOut) ::                                            &
     &  cftransfer(points)                                              &
                           ! Cloud fraction increment this tstep
     &, cfltransfer(points)                                             &
                           ! Liquid cloud fraction inc this tstep
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
     &  aplusb(points)                                                  &
                          ! A+B terms in diffusional growth equation
     &, tempw_dep(points)                                               &
                          ! Temporary for available moisture
     &, tempw_sub(points)                                               &
                          ! Temporary for saturation deficit
     &, pr02(points)                                                    &
                          ! Temporary in calculation of PSD slopes
     &, lamr1(points)                                                   &
                          ! Power of PSD slope
     &, lamr2(points)                                                   &
                          ! Power of PSD slope
     &, dqi(points)                                                     &
                          ! Temporary in calculating ice transfers
     &, dqil(points)                                                    &
                          ! Temporary in calculating liquid transfers
     &, dqi_dep(points)                                                 &
                          ! Deposition amount / kg kg-1
     &, dqi_sub(points)                                                 &
                          ! Sublimation amount / kg kg-1
     &, cfltemp(points)                                                 &
                          ! Temporary in calc of cfl and qcl change
     &, deltacf(points)                                                 &
                          ! Change in cf across timestep
     &, deltacfl(points)  ! Change in cfl across timestep

      Real                                                              &
     &  hm_rate                                                         &
                          ! Hallett-Mossop rate multiplier
     &, hm_normalize                                                    &
                          ! Normalization for T func. in HM process
     &, m_1(points)                                                     &
                          ! 1st moment of generic particle size dist. 
     &, m_0p5_dp3(points)
                          ! 1+(di+1)/2 moment of generic PSD 

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

      !-----------------------------------------------
      ! Hallett Mossop process normalisation
      !-----------------------------------------------
      hm_normalize=1.0/(1.0-exp((hm_t_min-hm_t_max)*hm_decay))

      Do i=1, points

        If (qcf(i) >  m0 .and. T(I) <  zerodegc) then

          !-----------------------------------------------
          ! Diffusional growth parameters
          !-----------------------------------------------
          aplusb(i) = (apb1-apb2*t(i)) * esi(i)
          aplusb(i) = aplusb(i) + (T(i)**3) * p(i) * apb3

          !-----------------------------------------------
          ! Moisture available from subgrid scale calculation
          !-----------------------------------------------
          tempw_dep(i) = qsl(i) * area_mix(i)                           &
     &                 + min(q_ice_1(i),qsl(i)) * area_ice_1(i)         &
     &                 - qs(i) * (area_mix(i) + area_ice_1(i))
          tempw_sub(i) = (q_ice_2(i) - qs(i)) * area_ice_2(i)

          ! Only allow the sub saturation to be reduced by the
          ! fraction of ice in this category
          tempw_sub(i) = tempw_sub(i)*frac_in_category(i)

          !-----------------------------------------------
          ! Calculate transfer rates
          !-----------------------------------------------
          If (l_psd) then
            ! Use generic particle size distribution
            ! constp(83) = 2 pi axial_ratio_correction
            ! constp(84) = ventilation coefficient 1
            ! constp(85) = ventilation coefficient 2
            !        * Sc^(1/3)*ci^0.5/viscosity0^0.5 
            dqi(i) = constp(83) * T(I)**2 * esi(I)                      &
     &                          * (constp(84)*m_1(i)*corr2(i)           &
     &                          +  constp(85)*rocor(i)*m_0p5_dp3(i))    &
     &                          / (qs(i) * aplusb(i) * rho(i))

          Else
            ! Use particle size distribution based on intercepts
            pr02(i)=rho(i)*qcf(i)*cficei(i)*constp(5+cry_offset)*tcgi(i)
            lamr1(I) = pr02(i)**cx(4+cry_offset)
            lamr2(I) = pr02(i)**CX(5+cry_offset)
            dqi(i) = tcg(i) * constp(6+cry_offset) * T(I)**2 *          &
     &               esi(I) * (constp(7+cry_offset) * corr2(i) *        &
     &               lamr1(i) + constp(8+cry_offset) * rocor(i) *       &
     &               lamr2(i)) / (qs(i) * aplusb(i) * rho(i))

          End if  ! l_psd

          dqi_dep(i) = dqi(i) * tempw_dep(i)
          dqi_sub(i) = dqi(I) * tempw_sub(i)


          If (dqi_dep(i) >  0.0) then  ! Limits depend on whether
                                       ! deposition or sublimation

            !-----------------------------------------------
            ! Deposition is occuring
            !-----------------------------------------------
            If (ice_type  ==  0) then  ! Only for crystals

              !-----------------------------------------------
              ! Hallett Mossop Enhancement
              !-----------------------------------------------
              If ( (T(i)-zerodegc)  >=  hm_t_max) then  ! no enhancement
                hm_rate=0.0
              Elseif ((T(i)-zerodegc)  <   hm_t_max                     &
     &          .and. (T(I)-zerodegc)  >   hm_t_min) then
                ! Some enhancement occurs between temperature limits
                hm_rate = (1.0-exp( (T(I)-zerodegc-hm_t_max)*hm_decay) )&
     &                    * hm_normalize
              Else  ! Some enhancement at lower temperatures
                hm_rate = exp( (T(i)-zerodegc - hm_t_min) * hm_decay )
              End if

              ! Calculate enhancement factor for HM process
              hm_rate = 1.0 + hm_rate * qcl(i) * hm_rqcl
              dqi_dep(i) = dqi_dep(i) * hm_rate

            End if  ! ice_type  ==  0

            !-----------------------------------------------
            ! Molecular diffusion to/from a surface is more efficient
            ! when a particle is at a molecular step. This is more
            ! likely for sublimation. For growth, reduce rate by 10%
            ! ----------------------------------------------
             dqi_dep(i)=0.9*dqi_dep(i)

            !-----------------------------------------------
            ! Latent heat correction (equivalent to aL in LS cloud)
            !-----------------------------------------------
            tempw_dep(i)=tempw_dep(i)*lheat_correc_ice(i)

            ! Only allow the supersaturation to be reduced by
            ! the fraction of ice in this category
            tempw_dep(i) = tempw_dep(i)*frac_in_category(i)

            !-----------------------------------------------
            ! Calculate available moisture and transfer amount.
            !-----------------------------------------------
            If (cfliq(i) >  0.0) then
            ! Include liquid water contribution. Ignore latent
            ! heat correction from freezing liquid.
              tempw_dep(i)=tempw_dep(i)+qcl(i)*area_mix(i)/cfliq(i)     &
     &                    +max((q_ice_1(i)-qsl(i)),0.0)*area_ice_1(i)
            End if
            dqi_dep(I)=min(dqi_dep(I)*timestep,tempw_dep(i))

!!fra298
!            !-----------------------------------------------
!            ! Update cloud fractions due to deposition in ice only region
!            !-----------------------------------------------
            If (cfliq(i) ==  0.0) then
            If (area_ice_1(i) > 0.0) then
            If (l_update_cf) Then
              deltacf(i) = area_ice_1(i) * ( sqrt( max(                 &
     &                 1.0+dqi_dep(i)*cfice(i)/(qcft(i)*area_ice_1(i))  &
     &                 , 0.0) ) - 1.0 )
!dhb599 20110601--scaling (down) this new cloud fraction increment as per fra298
!(in the purose of cooling down the system by reducing high cloud coverage)
!              write(6,*)'XXX lsp_deposition: xfactor = ', xfactor
              deltacf(i) = deltacf(i) * xfactor
              If (l_seq) then
                cff(i) = cff(i) + deltacf(i)
                cf (i) = cf (i) + deltacf(i)
                if (cff(i) > 1.0) then
                 deltacf(i)=1.0-cff(i)
                 cff(i) = 1.0
                endif 
                if (cf(i) > 1.0) then
                 cf(i) = 1.0
                endif
              End if
              cftransfer(i) =cftransfer(i) +deltacf(i)                  &
     &                        /(timestep*iterations)
              cfftransfer(i)=cfftransfer(i)+deltacf(i)                  &
     &                        /(timestep*iterations)
            End if  ! l_update_cf
            End if  ! areaice1 gt 0 
            End if  ! cfliq eq 0.
!!!fra298 end

          End if  ! dqi_dep gt 0.0

          If (dqi_sub(i)  <   0.0) then
            !-----------------------------------------------
            ! Sublimation is occuring
            !-----------------------------------------------
            ! Limits are spare moisture capacity and QCF
            ! outside liquid cloud
            dqi_sub(i) = max( max( dqi_sub(i) * timestep                &
     &                           , tempw_sub(i) * lheat_correc_ice(i) ) &
     &                 , -( qcf(i) * area_ice_2(i) * cficei(i) ))

            !-----------------------------------------------
            ! Update cloud fractions
            !-----------------------------------------------
            If (l_update_cf) Then
              deltacf(i) = area_ice_2(i) * ( sqrt( max(                 &
     &                 1.0+dqi_sub(i)*cfice(i)/(qcft(i)*area_ice_2(i))  &
     &                 , 0.0) ) - 1.0 )
              If (l_seq) then
                cff(i) = cff(i) + deltacf(i)
                cf (i) = cf (i) + deltacf(i)
              End if
              cftransfer(i) =cftransfer(i) +deltacf(i)                  &
     &                        /(timestep*iterations)
              cfftransfer(i)=cfftransfer(i)+deltacf(i)                  &
     &                        /(timestep*iterations)
            End if  ! l_update_cf
          End if  ! dqi_sub lt 0.0

          !-----------------------------------------------
          ! Adjust ice content
          !-----------------------------------------------
          If (l_seq) then
            qcf(i) = qcf(i) + dqi_dep(i) + dqi_sub(i)
          End if

          !-----------------------------------------------
          ! Calculate liquid water change
          !-----------------------------------------------
          If (cfliq(i) >  0.0 .and. area_mix(i) >  0.0                  &
     &            .and. qcl(i) >  0.0) then
            ! Deposition removes some liquid water content
            ! First estimate of the liquid water removed is explicit
            dqil(i) = max (min ( dqi_dep(i)*area_mix(i)                 &
     &                /(area_mix(i)+area_ice_1(i)),                     &
     &                qcl(i)*area_mix(i)/cfliq(i)) ,0.0)

          !-----------------------------------------------
          ! Update liquid cloud fraction (and new liquid water est.)
          !-----------------------------------------------
            If (l_update_cf) Then
              ! First estimate of the liquid cloud fraction is based
              ! on this explicit estimate of liquid lost by deposition
              cfltemp(i)=cfl(i)*sqrt(max(1.0-dqil(i)/qcl(i),0.0))
              ! Now form a half timestep estimate of the proportion
              ! of the depositing volume which contains liquid cloud
              cfltemp(i)=0.5*max(area_mix(i)-(cfl(i)-cfltemp(i)),0.0)   &
     &                         /(area_mix(i)+area_ice_1(i))             &
     &                 + 0.5*area_mix(I)/ (area_mix(i)+area_ice_1(i))
              ! Recalculate an estimate of the liquid water removed
              dqil(i)=max( min( dqi_dep(i)*cfltemp(i),                  &
     &                          qcl(i)*area_mix(I)/cfliq(i)),0.0)

              ! Update liquid cloud fraction and transfer rate
              deltacfl(i) = cfl(i)
              If (l_seq) then
                cfl(i) = cfl(i) * sqrt(max(1.0-dqil(i)/qcl(i),0.0))
              End if
              cfltransfer(i) = cfltransfer(i) + (cfl(i) - deltacfl(i))  &
     &                         / (timestep*iterations)
            End If  ! l_update_cf

          Else
            ! Deposition does not remove any liquid water content
            dqil(i)=0.0
          End if ! cfliq gt 0.0 etc.

          !-----------------------------------------------
          ! Adjust liquid and vapour contents (liquid adjusts first)
          !-----------------------------------------------
          If (l_seq) Then
            qcl(i) = qcl(i) - dqil(i)  ! Bergeron Findeisen acts first
            T(i) = T(i) + lfrcp * dqil(i)
            dqi(i) = dqi_dep(i) + dqi_sub(i)- dqil(i)

            q(i) = q(i) - dqi(i)
            T(i) = T(i) + lsrcp * dqi(i)
          End if  ! l_seq

          !-----------------------------------------------
          ! Store depostion/sublimation rate
          !-----------------------------------------------
          ptransfer(i) = ptransfer(i) + (dqi_dep(i) + dqi_sub(i))       &
     &                   / (timestep*iterations)

        End if  ! qcf  >   m0 etc.

      End do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_DEPOSITION
