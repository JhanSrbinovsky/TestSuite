
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Capture of raindrops by ice
! Subroutine Interface:
      SUBROUTINE LSP_CAPTURE(                                           &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, qcf, qrain, T                                                   &
                                          ! Water contents and temp
     &, area_liq, area_mix, area_ice                                    &
     &, cficei                                                          &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
     &, rainfrac, rain_liq, rain_mix                                    &
                                          ! Rain fraction information
     &, rain_ice, rain_clear                                            &
!    &, cf, cff                           ! Current cloud fractions for
!                                         ! updating
     &, rho, rhor, m0, tcg, tcgi                                        &
                                          ! Parametrization information
     &, corr, dhilsiterr                                                &
     &, cx, constp, ai, bi, aic, bic, lfrcp , ice_type                  &
                                          ! Microphysical information
     &, l_update_cf, l_seq, l_psd, l_mcr_qrain                          &
                                          ! Code options
     &, ptransfer, iterations                                           &
                                          ! Mass transfer diagnostic
!    &, cftransfer,cfftransfer            ! Cloud transfer diagnostics
     &, rftransfer                                                      &
                                          ! Rain fraction transfer
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of capture of raindrops
!   by ice particles
!
! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles sweeping out a specified
!   distribution of raindrops.
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
     &, area_liq(points)                                                &
                          ! Fraction of gridbox with liquid-only cloud
     &, area_mix(points)                                                &
                          ! Fraction of gridbox with mixed phase cloud
     &, area_ice(points)                                                &
                          ! Fraction of gridbox with ice-only cloud
     &, cficei(points)                                                  &
                          ! 1/Fraction of gridbox with ice cloud
!    &, cf(points)        ! Current cloud fraction
!    &, cff(points)       ! Current ice cloud fraction
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, rhor(points)                                                    &
                          ! 1/Air density / m3 kg-1
     &, m0                                                              &
                          ! Seed ice water content / kg kg-1
     &, tcg(points)                                                     & 
                          ! T dependent function in ice size dist'n
     &, tcgi(points)                                                    &
                          ! 1/tcg (no units)
     &, corr(points)                                                    &
                          ! Fall velocity correction factor (no units)
     &, dhilsiterr(points)                                              &
                          ! Depth of layer / timestep  / m s-1
     &, ai, bi, aic, bic                                                &
                          ! Ice mass size relationship m(D)=ai D^bi
     &, lfrcp             ! Latent heat of fusion
                          ! / heat capacity of air / K
!
      Real, Intent(InOut) ::                                            &
     &  qcf(points)                                                     &
                          ! Ice water content    / kg kg-1
     &, qrain(points)                                                   &
                          ! Rain mixing ratio / kg kg-1
     &, T(points)                                                       &
                          ! Temperature / K
     &, rainfrac(points)                                                &
                          ! Rain fraction
     &, rain_liq(points)                                                &
                          ! Overlap fraction of rain and liquid
     &, rain_mix(points)                                                &
                          ! Overlap frac of rain with mixed phase cloud
     &, rain_ice(points)                                                &
                          ! Overlap fraction of rain and ice
     &, rain_clear(points)                                              &
                          ! Overlap fraction of rain and clear-sky
     &, ptransfer(points)                                               &
                          ! Mass rimed in this timestep / kg kg-1
     &, rftransfer(points)! Transfer of rain fraction this timestep
!    &, cftransfer(points) ! Cumulative cloud fraction increment
!    &, cfftransfer(points)! Cumulative ice cloud fraction increment
!
      Logical, Intent(In) ::                                            &
     &  l_update_cf                                                     &
                          ! Update cloud fractions
     &, l_seq                                                           &
                          ! Carry out sequential updating
     &, l_psd                                                           &
                          ! Use generic ice particle size distribution
     &, l_mcr_qrain       ! Use prognostic rain mixing ratio
!
! Local Variables
!
      Integer                                                           &
     &  i                                                               &
                          ! Loop counter for points
     &, cry_offset        ! Index offset for ice crystals
!
      Real                                                              &
     &  dpr(points)                                                     &
                          ! Transfer of mixing ratio  / kg kg-1
     &, vr1(points)                                                     &
                          ! Average fall speed of raindrops  / m s-1
     &, vi1(points)                                                     &
                          ! Average fall speed of ice particles  / m s-1
     &, fv1(points)                                                     &
                          ! Average fall speed difference between
                          ! ice particles and raindrops  / m s-1
     &, lamr1(points)                                                   &
                          ! Reciprocal of slope parameter in raindrop
                          ! size distribution  / m
     &, lami1(points)                                                   &
                          ! Reciprocal of slope parameter in ice
                          ! particle size distribution  / m
     &, lamfac1(points)                                                 &
                          ! Combination of lamr1 and lamr2
     &, rf_final(points)                                                &
                          ! Rain fraction at end of the timestep
     &, m_0(points), m_1(points), m_2(points)                           &
                          ! Moments of the ice particle size distributn.
     &, m_bi_di(points)
                          ! bi+di moment of the generic ice size distn

      !-----------------------------------------------
      ! Select appropriate ice parametrization (see c_lspsiz)
      !-----------------------------------------------
      cry_offset=ice_type*ice_type_offset

      !-----------------------------------------------
      ! Calculate moments of size distribution if appropriate 
      !----------------------------------------------- 
      If (l_psd) then
        ! Calculate the 0th, 1st, 2nd and bi+di (cx(82)) moments of the
        ! ice particle size distribution 
! DEPENDS ON: lsp_moments 
        Call lsp_moments(points,rho,T,qcf,cficei,                       & 
     &                   ai,bi,0.0,m_0)
! DEPENDS ON: lsp_moments 
        Call lsp_moments(points,rho,T,qcf,cficei,                       & 
     &                   ai,bi,1.0,m_1)
! DEPENDS ON: lsp_moments 
        Call lsp_moments(points,rho,T,qcf,cficei,                       & 
     &                   ai,bi,2.0,m_2)
! DEPENDS ON: lsp_moments 
        Call lsp_moments(points,rho,T,qcf,cficei,                       & 
     &                   ai,bi,cx(82),m_bi_di)
      End if

      Do i = 1, points

        If (qcf(i)  >   m0 .and. qrain(i)  >   0.0                      &
     &      .and. (rain_mix(i)+rain_ice(i))  >   0.0                    &
     &      .and. T(i)  <   zerodegc) Then

          !-----------------------------------------------
          ! Calculate rain mass-weighted fallspeed
          !-----------------------------------------------
          If (L_mcr_qrain) then
            ! rain is a mixing ratio (kg kg-1)
            vr1(i) = (constp(41)/6.0) * corr(i) *                       &
     &             ( rho(i)*constp(50)*qrain(i)/rainfrac(i) )**cx(51)

          Else
            ! rain is a flux (kg m-2 s-1)
            vr1(i) = corr(i) * constp(41)/6.0 *                         &
     &               (qrain(i) * rho(i) * dhilsiterr(i)                 &
     &               /( rainfrac(i)*constp(42)*corr(i)) )**cx(41)

          End if  ! L_mcr_qrain

          !-----------------------------------------------
          ! Calculate ice mass-weighted fallspeed
          !-----------------------------------------------
          If (l_psd) then
            ! Use the generic PSD
            ! constp(82) = ci*ai
            vi1(i) = constp(82) * corr(i) * m_bi_di(i)                  &
     &                 / (rho(i) * qcf(i) * cficei(i))

          Else
            vi1(i) = constp(4+cry_offset) * corr(i) *                   &
     &             (rho(i) * qcf(i) * cficei(i)                         &
     &             * constp(5+cry_offset) * tcgi(i))**cx(3+cry_offset)
          End if  ! l_psd

          !-----------------------------------------------
          ! Estimate the mean absolute differences in velocities
          !-----------------------------------------------
          fv1(i) = max(abs(vr1(i)-vi1(i)),(vr1(i)+vi1(i))/8.0)

          !-----------------------------------------------
          ! Calculate reciprocal of lambda for rain
          !-----------------------------------------------
          If (l_mcr_qrain) then
            ! rain is a mixing ratio (kg kg-1)
            lamr1(i) = (rho(i) * constp(50) * qrain(i) / rainfrac(i))   &
     &                 **(cx(52))

          Else
            ! rain is a flux (kg m-2 s-1)
            lamr1(I) = (qrain(i) * rho(i) * dhilsiterr(i)               &
     &                 /( rainfrac(i)*constp(42)*corr(i)) )**(cx(42))

          End If  ! l_mcr_qrain

          !-----------------------------------------------
          ! Calculate reciprocal of lambda for ice crystals
          !-----------------------------------------------
          If (.not. l_psd) then
            lami1(i) = (rho(i) * qcf(i) * cficei(i)                     &
     &             *constp(5+cry_offset)*tcgi(I))**(-cx(7+cry_offset))
          End if

          !------------------------------------------------
          ! Calculate transfer
          !------------------------------------------------
          If (l_psd) then
            ! Use the generic ice particle size distribution
            ! constp(86)=pi**2/24 x1r rho_water
            ! constp(87)=gamma(6+x4r)
            ! constp(88)=2 gamma(5+x4r)
            ! constp(89)=gamma(4+x4r)
            dpr(i) = constp(86) * fv1(i) *  timestep * rhor(i) *        &
     &              (rain_mix(i)+rain_ice(i)) * lamr1(i)**(-cx(46)) *   &
     &              (lamr1(i)**cx(43) * constp(87) * m_0(i) +           &
     &               lamr1(i)**cx(44) * constp(88) * m_1(i) +           &
     &               lamr1(i)**cx(45) * constp(89) * m_2(i) )
          Else
            lamfac1(i) = constp(10+cry_offset) * constp(43) *           &
     &               (lamr1(i)**cx(43) * lami1(i)**cx(8+cry_offset)) +  &
     &               constp(11+cry_offset) * constp(44) *               &
     &               (lamr1(i)**cx(44) * lami1(i)**cx(9+cry_offset)) +  &
     &               constp(12+cry_offset) * constp(45) *               &
     &               (lamr1(i)**cx(45) * lami1(i)**cx(10+cry_offset))

            dpr(i) = tcg(i) * constp(13+cry_offset) *                   &
     &             lami1(i)**(-cx(11+cry_offset)) *                     &
     &             lamr1(i)**(-cx(46)) * fv1(i) * lamfac1(i) *          &
     &             timestep * rhor(i) * (rain_mix(i)+rain_ice(i))
          End if  ! l_psd

          ! Limit transfer to the mass of rain that is available
          dpr(i) = min(dpr(i),qrain(i) *                                &
     &                      (rain_mix(i)+rain_ice(i))/rainfrac(i))

          ptransfer(i) = ptransfer(i) + dpr(i) / (timestep*iterations)

          !------------------------------------------------
          ! Adjust ice and rain contents
          !------------------------------------------------
          If (l_seq) then
            qcf(i)   = qcf(i)   + dpr(i)
            qrain(i) = qrain(i) - dpr(i)
            T(i)     = T(i)     + dpr(i) * lfrcp
          End if

          !------------------------------------------------
          ! Update cloud fractions
          !------------------------------------------------
!         These are commented out since there is currently no
!         cloud fraction update associated with the capture term.
!          cf_transfer_rate(i)  = 0.0 / (timestep*iterations)
!          cff_transfer_rate(i) = 0.0 / (timestep*iterations)
!
!          If (l_update_cf) then
!            cftransfer(i)  = cftransfer(i)  + cf_transfer_rate(i)
!            cfftransfer(i) = cfftransfer(i) + cff_transfer_rate(i)
!
!            If (l_seq) then
!              cf(i)  = cf(i)  +cf_transfer_rate(i) *timestep*iterations
!              cff(i) = cff(i) +cff_transfer_rate(i)*timestep*iterations
!            End if  ! l_seq
!
!          End if  ! l_update_cf

           !-----------------------------------------------
           ! Update rain fractions
           !-----------------------------------------------
!          These are commented out to ensure that rainfraction tends
!          to a 0 or 1 behaviour as RHcrit tends to 1.
!           rf_final(i) = rainfrac(i) * qrain(i) / (qrain(i)+dpr(i))
!           rftransfer(i) = rftransfer(i) + (rf_final(i) - rainfrac(i))
!
!           If (L_seq) then
!             rainfrac(i)= rf_final(i)
!             rain_liq(i)= min(area_liq(i),rainfrac(i))
!             rain_mix(i)= min(area_mix(i),rainfrac(i)-rain_liq(i))
!             rain_ice(i)=
!     &          min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
!             rain_clear(i)=
!     &            rainfrac(i) - rain_liq(i)-rain_mix(i)-rain_ice(i)
!           End if  ! L_seq

         End If ! qcf(i) >  m0 etc.

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_CAPTURE
