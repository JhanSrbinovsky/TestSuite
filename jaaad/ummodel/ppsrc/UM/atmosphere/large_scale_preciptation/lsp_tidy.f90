
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Riming of ice particles
! Subroutine Interface:
      SUBROUTINE LSP_TIDY(                                              &
     &  points, timestep, iterations                                    &
                                          ! Number of points and tstep
     &, q, qcl, qcf, qcf2, qrain, T                                     &
                                          ! Water contents and temp
     &, area_liq, area_mix, area_ice                                    &
                                          ! Cloud fraction information
     &, cfice, cficei                                                   &
                                          ! at start of microphysics ts
     &, cf, cfl, cff                                                    &
                                          ! Current cloud fractions for
                                          ! updating
     &, rainfrac, rain_liq, rain_mix                                    &
                                          ! Rain fractions
     &, rain_ice, rain_clear                                            &
     &, q_ice, qs, qsl, snow_agg, snow_cry                              &
                                          ! Other water contents
     &, rho, rhor, p                                                    &
                                          ! Other model prognostics
     &, cttemp,dhi,dhilsiterr,frac_ice_above                            &
                                             ! Other information
     &, lcrcp, lfrcp, lsrcp                                             &
                                          ! Microphysical information
     &, l_update_cf, l_seq, l_mcr_qcf2                                  &
                                          ! Code options
     &, l_it_melting                                                    &
     &, psdep, pidep, psmlt, pimlt, prevp                               &
                                          ! Mass transfer diagnostics
     &, cftransfer,cfltransfer,cfftransfer                              &
                                          ! Cloud transfer diagnostics
     &, rftransfer                                                      &
                                          ! Rain transfer diagnostics
     &  )
!
      Implicit None
!
! Purpose:
!   Tidy up small numerical values after large-scale precipitation.
!   Ideally, this code would be unnecessary, but in reality there
!   are always going to be small values left over numerically that
!   should be reset.
!
! Method:
!   1. Evaporate rain amounts
!   2. Evaporate small ice amounts
!   3. Reset cloud fractions
!   4. Melt small snow amounts
!
! Current Owner of Code: Jonathan Wilkinson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!
! Subroutine Arguments
!
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
                          ! Number of points
     &, iterations        ! Number of iterations

      Real, Intent(In) ::                                               &
     &  timestep                                                        &
                          ! Timestep of one iteration / s
     &, qcl(points)                                                     &
                          ! Liquid water content / kg kg-1
     &, area_liq(points)                                                &
                          ! Fraction of gridbox with liquid-only cloud
     &, area_mix(points)                                                &
                          ! Fraction of gridbox with mixed phase cloud
     &, area_ice(points)                                                &
                          ! Fraction of gridbox with ice-only cloud
     &, cfice(points)                                                   &
                          ! Fraction of gridbox with ice cloud
     &, cficei(points)                                                  &
                          ! 1/Fraction of gridbox with ice cloud
     &, q_ice(points)                                                   &
                          ! Vapour content in ice cloud / kg kg-1
     &, qs(points)                                                      &
                          ! Saturated water content wrt ice / kg kg-1
     &, qsl(points)                                                     &
                          ! Saturated water content wrt liquid / kg kg-1
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, rhor(points)                                                    &
                          ! 1 / air density / m3 kg-1
     &, p(points)                                                       &
                          ! Air pressure / N m-2
     &, dhi(points)                                                     &
                          ! Timestep/layer thickness / s m-1
     &, dhilsiterr(points)                                              &
                          ! 1/(dhi*iterations) / m s-1
     &, frac_ice_above(points)                                          &
                               ! Ice cloud fraction in layer above
     &, lcrcp                                                           &
                          ! Latent heat of condensation/cP / K
     &, lfrcp                                                           &
                          ! Latent heat of fusion/cP / K
     &, lsrcp             ! Latent heat of sublimation/cP / K

      Real, Intent(InOut) ::                                            &
     &  q(points)                                                       &
                          ! Vapour mixing ratio / kg kg-1
     &, qcf(points)                                                     &
                          ! Aggregates mixing ratio / kg kg-1
     &, qcf2(points)                                                    &
                          ! Crystals mixing ratio / kg kg-1
     &, qrain(points)                                                   &
                          ! Rain mixing ratio / kg kg-1
     &, T(points)                                                       &
                          ! Temperature / K
     &, cttemp(points)                                                  &
                          ! Ice-cloud top temperature / K
     &, cf(points)                                                      &
                          ! Current cloud fraction
     &, cfl(points)                                                     &
                          ! Current liquid cloud fraction
     &, cff(points)                                                     &
                          ! Current ice cloud fraction
     &, rainfrac(points)                                                &
                          ! Rain fraction (no units)
     &, rain_liq(points)                                                &
                          ! Overlap of rain with liquid (no units)
     &, rain_mix(points)                                                &
                          ! Overlap of rain with mixed phase region
     &, rain_ice(points)                                                &
                          ! Overlap of rain with ice
     &, rain_clear(points)                                              &
                          ! Overlap of rain with clear sky
     &, snow_agg(points)                                                &
                          ! Aggregate snowfall rate / kg m-2 s-1
     &, snow_cry(points)                                                &
                          ! Crystal snowfall rate / kg m-2 s-1
     &, psdep(points)                                                   &
                          ! Deposition of aggregates diag. / kg kg s-1
     &, pidep(points)                                                   &
                          ! Deposition of crystals diag. / kg kg s-1
     &, psmlt(points)                                                   &
                          ! Melting of aggregates diagnostic / kg kg s-1
     &, pimlt(points)                                                   &
                          ! Melting of crystals diagnostic / kg kg s-1
     &, prevp(points)                                                   &
                          ! Evaporation of rain diagnostic / kg kg s-1
     &, cftransfer(points)                                              &
                           ! Rate of change of bulk cloud frac / s-1
     &, cfltransfer(points)                                             &
                           ! Rate of change of liquid cloud frac / s-1
     &, cfftransfer(points)                                             &
                           ! Rate of change of ice cloud frac / s-1
     &, rftransfer(points)! Rate of change of rain fraction / s-1

!      Real, Intent(Out) ::

      Logical, Intent(In) ::                                            &
     &  l_update_cf                                                     &
                          ! Update cloud fractions
     &, l_seq                                                           &
                          ! Carry out sequential updating
     &, l_mcr_qcf2                                                      &
                          ! Use crystals as a prognostic
     &, l_it_melting      ! Iterative melting being used

!
! Local Variables
!
      Integer                                                           &
     &  i                 ! Loop counter for points

      Real                                                              &
     &  dpr(points)                                                     &
                          ! Mass transfer / kg kg-1
     &, temp7(points)                                                   &
                          ! Wet bulb temperature / deg C
     &, tempw(points)                                                   &
                          ! Saturation excess / kg kg-1
     &, cfnew(points)                                                   &
                          ! New cloud fraction
     &, cflnew(points)                                                  &
                          ! New liquid cloud fraction
     &, cffnew(points)                                                  &
                          ! New ice cloud fraction
     &, t_melt            ! Emergency melting temperature

      Do i=1,points

        If ((qrain(i)  <=  qcfmin .and. qcl(i) <= 0.0)                  &
     &       .or. qrain(i)  <   0.0) Then
          !-----------------------------------------------
          ! 1. If there is a tiny rain amount and no liquid then
          ! evaporate the rain
          !-----------------------------------------------
          dpr(i) = qrain(i)

          If (l_seq) Then
            ! Update prognostics
            q(i)   = q(i) + dpr(i)
            T(i)   = T(i) - dpr(i) * lcrcp
            qrain(i) = 0.0

            ! Update rain fractions
            rainfrac(i) = 0.0
            rain_liq(i) = 0.0
            rain_mix(i) = 0.0
            rain_ice(i) = 0.0
            rain_clear(i) = 0.0

          End If  ! l_seq

          ! Update evaporation rate
          rftransfer(i) = rftransfer(i)                                 &
     &                  - rainfrac(i)/(timestep*iterations)
          prevp(i) = prevp(i) + dpr(i)/(timestep*iterations)

        End If  ! rain lt qcfmin etc.

        !-----------------------------------------------
        ! 2. Evaporate small ice amounts
        !-----------------------------------------------
        If (L_mcr_qcf2) Then
          ! Evaporate small aggregate and crystal amounts separately

          !-----------------------------------------------
          ! 2a. Aggregates
          !-----------------------------------------------
          If ((qcf(i)  <   qcfmin) .and.                                &
     &      (T(i)  >   zerodegc .or.                                    &
     &      (q_ice(i)  <=  qs(i) .and. area_mix(i)  <=  0.0)            &
     &      .or. qcf(i) <  0.0)) then
            ! Ice is very small and T>0 and ice is not growing by
            ! deposition, so evaporate it.
            dpr(i) = qcf(i)

            ! Update prognostics
            If (l_seq) Then
              q(i)   = q(i) + dpr(i)
              T(i)   = T(i) - lsrcp * dpr(i)
              qcf(i) = 0.0
            End If  ! l_seq

            ! Update deposition rate
            psdep(i) = psdep(i) - dpr(i)/(timestep*iterations)

          End if  ! qcf < qcfmin etc

          !-----------------------------------------------
          ! 2b. Crystals
          !-----------------------------------------------
          If ((qcf2(i)  <   qcfmin) .and.                               &
     &      (T(i)  >   zerodegc .or.                                    &
     &      (q_ice(i)  <=  qs(i) .and. area_mix(i)  <=  0.0)            &
     &      .or. qcf2(i)  <   0.0)) Then
            ! Ice is very small and T>0 and ice is not growing by
            ! deposition, so evaporate it.
            dpr(i) = qcf2(i)

            ! Update prognostics
            If (l_seq) Then
              q(i)   = q(i) + dpr(i)
              T(i)   = T(i) - lsrcp * dpr(i)
              qcf2(i)= 0.0
            End If  ! l_seq

            ! Update deposition rate
            pidep(i) = pidep(i) - dpr(i)/(timestep*iterations)

          End If  ! qcf < 0 etc.

          If (qcf(i)  <=  0.0 .and. qcf2(i)  <=  0.0) Then
            !----------------------------------------------
            ! 2c. Update ice cloud amounts
            !----------------------------------------------
            If (l_update_cf) Then
              cfftransfer(i) = -cff(i)/(timestep*iterations)
              cftransfer(i)  = (cfl(i)-cff(i))/(timestep*iterations)
              If (l_seq) Then
                cff(i) = 0.0
                cf(i)  = cfl(i)
              End If  ! l_seq
            End If  ! l_update_cf

            !----------------------------------------------
            ! 2d. Update cloud top temperature
            !----------------------------------------------
            cttemp(i) = T(i)
          End If  ! qcf eq 0 etc

        Else  ! l_mcr_qcf2

          !----------------------------------------------
          ! 2e. One prognostic ice content is active
          !----------------------------------------------
          If ((qcf(i)  <   qcfmin) .and.                                &
     &      (T(i)  >   zerodegc .or.                                    &
     &      (q_ice(i)  <=  qs(i) .and. area_mix(i)  <=  0.0)            &
     &      .or. qcf(i) <= 0.0) )  Then
            ! Ice is very small and T>0 and ice is not growing by
            ! deposition, so evaporate it.
            dpr(i) = qcf(i)

            ! Update prognostics
            If (l_seq) Then
              q(i)   = q(i) + dpr(i)
              T(i)   = T(i) - lsrcp * dpr(i)
              qcf(i) = 0.0
            End If  ! l_seq

            ! Update deposition rate
            psdep(i) = psdep(i) - dpr(i)/(timestep*iterations)

            ! Update ice cloud amounts
            If (l_update_cf) Then
              cfftransfer(i) = -cff(i)/(timestep*iterations)
              cftransfer(i)  = (cfl(i)-cff(i))/(timestep*iterations)
              If (l_seq) Then
                cff(i) = 0.0
                cf(i)  = cfl(i)
              End If  ! l_seq
            End If  ! l_update_cf

            ! Update cloud top temperature
            cttemp(i) = T(i)

          End If  ! qcf lt qcfmin etc

        End If  ! L_mcr_qcf2

        !------------------------------------------------
        ! 3. Limit cloud fractions to physically reasonable values
        !------------------------------------------------
        ! It isn't clear that this ought to be a parallel calculation
        ! but it is allowed to be for the present time.
        If (l_update_cf) Then

          ! Calculate new cloud fraction values
          cffnew(i) = max(min(cff(i),1.0),0.0)
          cflnew(i) = max(min(cfl(i),1.0),0.0)
          cfnew (i) = max( max(cflnew(i),cffnew(i)) , cf(i) )
          cfnew (i) = min( min(cflnew(i)+cffnew(i),1.0), cfnew(i) )

          ! Calculate transfer rates
          cfftransfer(i) = cfftransfer(i)                               &
     &                   + (cffnew(i) - cff(i)) / (timestep*iterations)
          cfltransfer(i) = cfltransfer(i)                               &
     &                   + (cflnew(i) - cfl(i)) / (timestep*iterations)
          cftransfer(i)  = cftransfer(i)                                &
     &                   + (cfnew(i)  - cf (i)) / (timestep*iterations)

          ! Update cloud fractions
          If (l_seq) then
            cff(i) = cffnew(i)
            cfl(i) = cflnew(i)
            cf(i)  = cfnew(i)
          End if  ! l_seq

        End If  ! l_update_cf

        !------------------------------------------------
        ! 4. Emergency melting of snow to avoid excess snow at
        !    warm temperatures
        !------------------------------------------------
        If (l_it_melting) then
          ! Use a warmer temperature if iterative melting is active
          t_melt = zerodegc + 2.0
        Else
          t_melt = zerodegc
        End if  ! l_it_melting

        !------------------------------------------------
        ! Melt snow_agg first
        !------------------------------------------------
        If (snow_agg(i)  >   0.0 .and. T(i)  >   t_melt) Then

          ! Numerical approximation of wet bulb temperature excess
          tempw(i) = area_ice(i) * max(qsl(i)-q_ice(i),0.0)*cficei(i)
          temp7(i) = T(i) - zerodegc                                    &
     &               -tempw(i)*(tw1+tw2*(P(i)-tw3)-tw4*(T(i)-tw5))
          temp7(i) = max(temp7(i),0.0)

          ! Calculate transfer rate
          dpr(i) = temp7(i) / lfrcp ! Rate based on Tw excess
          ! Limit to the amount of snow available
          dpr(i) = min(dpr(i) , snow_agg(i)                             &
     &                        * dhi(i)*iterations*rhor(i) )

          ! Add to melting rate
          psmlt(i) = psmlt(i) + dpr(i)/(timestep*iterations)

          ! Update values of snow and rain
          If (l_seq) Then
            snow_agg(i) = snow_agg(i) - dpr(i)*rho(i)*dhilsiterr(i)
            qrain(i)    = qrain(i)    + dpr(i)
            T(i)        = T(i)        - dpr(i) * lfrcp
          End If  ! l_seq

        End If  !  snow_agg gt 0 etc

        !------------------------------------------------
        ! Melt snow_cry next
        !------------------------------------------------
        If (l_mcr_qcf2 .and.                                            &
     &      snow_cry(i)  >   0.0 .and. T(i)  >   t_melt) Then

          ! Numerical approximation of wet bulb temperature excess
          tempw(i) = area_ice(i) * max(qsl(i)-q_ice(i),0.0)*cficei(i)
          temp7(i) = T(i) - zerodegc                                    &
     &               - tempw(i)*(TW1+TW2*(P(i)-TW3)-TW4*(T(i)-TW5))
          temp7(i) = max(temp7(i),0.0)

          ! Calculate transfer rate
          dpr(i) = temp7(i) / lfrcp ! Rate based on Tw excess
          ! Limit to the amount of snow available
          dpr(i) = min(dpr(i) , snow_cry(i)                             &
     &                        * dhi(i)*iterations*rhor(i) )

          ! Add to melting rate
          pimlt(i) = pimlt(i) + dpr(i)/(timestep*iterations)

          ! Rain fraction will take on the value of the ice fraction
          rftransfer(i) = rftransfer(i) +                               &
     &      max( max(frac_ice_above(i),cff(i))-rainfrac(i) , 0.0)       &
     &      / (timestep*iterations)

          ! Update values of snow and rain
          If (l_seq) Then
            snow_cry(I) = snow_cry(i) - dpr(i)*rho(i)*dhilsiterr(i)
            qrain(i)    = qrain(i)    + dpr(i)
            T(i)        = T(i)        - dpr(i) * lfrcp

            ! Update rain fractions
            rainfrac(i) = max(rainfrac(i),cfice(i))
            rain_liq(i) = min(area_liq(i),rainfrac(i))
            rain_mix(i) = min(area_mix(i),rainfrac(i)-rain_liq(i))
            rain_ice(i) =                                               &
     &           min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
            rain_clear(i) =                                             &
     &           rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i)
          End If  ! l_seq

        End If  ! snow_cry lt 0 etc

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_TIDY
