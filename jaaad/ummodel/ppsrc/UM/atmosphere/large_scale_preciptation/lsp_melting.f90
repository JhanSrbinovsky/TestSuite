
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Melting of ice particles
! Subroutine Interface:
      SUBROUTINE LSP_MELTING(                                           &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, q, q_ice, qcf, qcft, qrain, qsl                                 &
                                          ! Water contents
     &, T, p                                                            &
                                          ! Temperature and pressure
     &, area_liq, area_mix, area_ice                                    &
                                          ! Partition information
     &, cfice, cficei, frac_ice_above                                   &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
     &, rainfrac, rain_liq, rain_mix                                    &
                                          ! Rain fraction information
     &, rain_ice, rain_clear                                            &
     &, cf, cff                                                         &
                                          ! Current cloud fractions for
!                                         ! updating
     &, rho, rhor, m0, tcg, tcgi                                        &
                                          ! Parametrization information
     &, corr2, rocor                                                    &
     &, cx, constp, ai, bi, aic, bic, lfrcp, ice_type                   &
                                          ! Microphysical information
     &, l_update_cf, l_seq, l_psd                                       &
                                          ! Code options
     &, ptransfer, iterations                                           &
                                          ! Mass transfer diagnostic
     &, cftransfer,cfftransfer                                          &
                                          ! Cloud transfer diagnostics
     &, rftransfer                                                      &
                                          ! Rain fraction transfer
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of melting of ice particles
!
! Method:
!   Solve implicitly the microphysical transfer equation for a
!   specified distribution of ice particles melting at an
!   approximated wet-bulb temperature
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
     &, ice_type          ! Type of ice (0 - crystals, 1 - aggregates
                          !              3 - graupel)
!
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, q(points)                                                       &
                          ! Gridbox mean vapour content / kg kg-1
     &, q_ice(points)                                                   &
                          ! Vapour content in ice partition / kg kg-1
     &, qcft(points)                                                    &
                          ! Total ice water content / kg kg-1
     &, qsl(points)                                                     &
                          ! Saturated mixing ratio wrt liquid / kg kg-1
     &, p(points)                                                       &
                          ! Pressure / N m-2
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
     &, frac_ice_above(points)                                          &
                               ! Ice fraction in layer above this layer
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, rhor(points)                                                    &
                          ! 1/Air density / m3 kg-1
     &, corr2(points)                                                   &
                          ! Temperature dependency of diffusion params
     &, rocor(points)                                                   &
                          ! Combination of corr and corr2
     &, m0                                                              &
                          ! Seed ice water content / kg kg-1
     &, tcg(points)                                                     & 
                          ! T dependent function in ice size dist'n
     &, tcgi(points)                                                    &
                          ! 1/tcg (no units)
     &, ai, bi, aic, bic                                                &
                          ! Ice mass size relationship. m(D) = AI D^BI
     &, lfrcp             ! Latent heat of fusion
                          ! / heat capacity of air / K
!
      Real, Intent(InOut) ::                                            &
     &  qcf(points)                                                     &
                          ! Ice water content of this category / kg kg-1
     &, qrain(points)                                                   &
                          ! Rain mixing ratio / kg kg-1
     &, T(points)                                                       &
                          ! Temperature / K
     &, cf(points)                                                      &
                          ! Current cloud fraction
     &, cff(points)                                                     &
                          ! Current ice cloud fraction
     &, rainfrac(points)                                                &
                          ! Rain fraction
     &, rain_liq(points)                                                &
                          ! Overlap fraction of rain and liquid
     &, rain_mix(points)                                                &
                          ! Overlap fraction of rain and mixed phase
     &, rain_ice(points)                                                &
                          ! Overlap fraction of rain and ice
     &, rain_clear(points)                                              &
                          ! Overlap fraction of rain and clear sky
     &, ptransfer(points)                                               &
                          ! Mass rimed in this timestep / kg kg-1
     &, rftransfer(points)                                              &
                          ! Transfer of rain fraction this timestep
     &, cftransfer(points)                                              &
                           ! Cumulative cloud fraction increment
     &, cfftransfer(points)! Cumulative liquid cloud fraction increment
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
!
      Real                                                              &
     &  tempw(points)                                                   &
                          ! Average supersaturation in gridbox
     &, temp7(points)                                                   &
                          ! Wet bulb temperature minus 0 C / deg C
     &, pr02(points)                                                    &
                          ! Temporary in melting mass calculation
     &, dpr(points)                                                     &
                          ! Mass melted this timestep / kg kg-1
     &, delta_cf(points)                                                &
                          ! Cloud fraction transfered (no units)
     &, delta_cff(points)                                               &
                          ! Ice cloud fraction transfered (no units)
     &, rf_final(points)                                                &
                          ! Final rain fraction (no units)
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

      Do i = 1, points

        If (qcf(i) >  m0 .and. T(i) >  zerodegc) Then

          !-----------------------------------------------
          ! Calculate the average wet-bulb temperature (tempw)
          !-----------------------------------------------
          ! First calculate the average supersaturation
          If (ice_type  ==  3) then
            ! Graupel does not use ice cloud fraction
            tempw(i) = max(qsl(i)-q(i),0.0)
          Else
            ! Use the full ice cloud fraction information
            tempw(i)=area_ice(i)*max(qsl(i)-q_ice(i),0.0)*cficei(i)
          End if

          ! Now calculate the average amount (temp7) by which the
          ! wet-bulb temperature exceeds 0 degrees C.
          temp7(i) = T(i) - zerodegc - tempw(i)                         &
     &             * (tw1 + tw2 * (p(i)-tw3) - tw4*(T(i)-tw5) )
          temp7(i) = max(temp7(i),0.0)

          !-----------------------------------------------
          ! Calculate melting rate
          !-----------------------------------------------
          If (l_psd) then
            ! Use generic particle size distribution
            ! constp(84) = ventilation coefficient 1
            ! constp(85) = ventilation coefficient 2
            ! constp(90) = 2 pi axial_ratio_correction
            !              * air_conductivity / Lf
            dpr(i) = rhor(i) * constp(90) * timestep *                  &
     &               ( constp(84)*m_1(i)*corr2(i) +                     &
     &                 constp(85)*rocor(i)*m_0p5_dp3(i) )
          Else
            ! Use particle size distribution based on intercepts
            If (ice_type  ==  3) then
              ! Graupel does not use ice cloud fraction
              pr02(i) = rho(i)*qcf(i)*constp(5+cry_offset)*tcgi(i)
            Else
              ! Use the full ice cloud information
              pr02(i) = rho(i)*qcf(i)*cficei(i) *                       &
     &                  constp(5+cry_offset)*tcgi(i)
            End if  ! ice_type==3

            dpr(i)  = tcg(i) * constp(14+cry_offset) * timestep *       &
     &                (constp(7+cry_offset) * corr2(i) *                &
     &                pr02(i)**cx(4+cry_offset)                         &
     &              + constp(8+cry_offset) * rocor(i) *                 &
     &                pr02(i)**cx(5+cry_offset)) * rhor(i)

          End if  ! l_psd

          !-----------------------------------------------
          ! Calculate transfer
          !-----------------------------------------------
          ! Solve implicitly in terms of temperature
          dpr(i) = temp7(i) * (1.0-1.0/(1.0+dpr(i)*lfrcp))/lfrcp

          ! Limit mass transfer to the mass available
          dpr(i) = min(dpr(i),qcf(i))

          ptransfer(i) = ptransfer(i) + dpr(i) / (timestep*iterations)

          !------------------------------------------------
          ! Adjust ice and rain contents
          !------------------------------------------------
          If (l_seq) then
            qcf(i)   = qcf(i)   - dpr(i)
            qrain(i) = qrain(i) + dpr(i)
            T(i)     = T(i)     - dpr(i) * lfrcp
          End if

          !------------------------------------------------
          ! Update rain fractions
          !------------------------------------------------
          If (dpr(i)  >   0.0) Then
            rf_final(i) = max(max(rainfrac(i),cff(i)),frac_ice_above(i))

            rftransfer(i) = rftransfer(i) + (rf_final(i) - rainfrac(i))

            If (L_seq) then
              rainfrac(i)= rf_final(i)
              rain_liq(i)= min(area_liq(i),rainfrac(i))
              rain_mix(i)= min(area_mix(i),rainfrac(i)-rain_liq(i))
              rain_ice(i)=                                              &
     &           min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
              rain_clear(i)=                                            &
     &             rainfrac(i) - rain_liq(i)-rain_mix(i)-rain_ice(i)
            End if  ! L_seq

          End If  ! dpr  >   0

          !------------------------------------------------
          ! Update cloud fractions
          !------------------------------------------------
          If (l_update_cf) then

            ! Assume that ice cloud fraction is reduced in
            ! proportion to the mass removed
            delta_cff(i) = -cff(i) * dpr(i) / qcft(i)

            ! Assume a random overlap between the removed ice
            ! cloud and the liquid cloud
            delta_cf(i)  = delta_cff(i) * area_ice(i) * cficei(i)

            ! Calculate transfer rate diagnostics
            cftransfer(i)  = cftransfer(i)                              &
     &                       + delta_cf(i)  / (timestep*iterations)
            cfftransfer(i) = cfftransfer(i)                             &
     &                       + delta_cff(i) / (timestep*iterations)

            If (l_seq) then
              cf(i)  = cf (i) + delta_cf (i)
              cff(i) = cff(i) + delta_cff(i)
            End if  ! l_seq

          End if  ! l_update_cf

        End If  ! qcf  >   m0 etc.

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_MELTING
