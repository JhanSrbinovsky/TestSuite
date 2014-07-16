#if defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Initialisation of variables
! Subroutine Interface:
      SUBROUTINE LSP_INIT(                                              &
     &  points                                                          &
                                     ! Number of points
     &, timestepfixed, timestep                                         &
                                     ! Timesteps
     &, l_mcr_qcf2, l_mcr_qrain                                         &
                                     ! Microphysics logicals
     &, l_mcr_qgraup, l_psd, l_crystals                                 &
     &, l_non_hydrostatic, l_mixing_ratio, l_cry_agg_dep                &
     &, cx, constp                                                      &
                                     ! Microphysical array sizes
     &, T, p, cttemp                                                    &
                                     ! Temperature and pressure
     &, deltaz,rhodz,rhodz_dry,rhodz_moist                              &
                                           ! Air density information
     &, rho, rhor                                                       &
     &, dhi, dhir, dhilsiterr                                           &
                                     ! Tstep and layer thick. ratios
     &, q, qcl, qcf, qcf2, qrain, qgraup, rainrate                      &
                                                   ! Water contents
     &, qcf_agg, qcf_cry, qcf_tot, frac_agg, frac_cry_dep, frac_agg_dep &
     &, qs, qsl, esi, esw                                               &
                                     ! Saturation water quantities
     &, psdep, psaut, psacw, psacr                                      &
                                     ! Aggregate transfer rates
     &, psaci, psmlt, psmltevp,psfall                                   &
     &, pifrw, piprm, pidep, piacw                                      &
                                     ! Ice transfer rates
     &, piacr, pimlt, pimltevp,pifall                                   &
     &, praut, pracw, prevp, prfall                                     &
                                     ! Rain transfer rates
     &, plset, plevpset                                                 &
                                     ! Droplet settling transfers
     &, pgaut, pgacw, pgacs, pgmlt                                      &
                                     ! Graupel transfer rates
     &, pgfall                                                          &
     &, cf_transfer_diag, cfl_transfer_diag                             &
     &, cff_transfer_diag, rf_transfer_diag                             &
     &, snow_cry, snow_agg, snowt_cry                                   &
                                     ! Precipitation rates
     &, snowt_agg, rainratet, graupratet                                &
     &, lheat_correc_liq, lheat_correc_ice                              &
                                           ! Latent heat corrections
     &, corr, corr2, rocor                                              &
                                     ! Fall speed and diffusivity
                                     ! corrections
     &, tcg,tcgi,tcgc,tcgci,tcgg,tcggi                                  &
                                           ! Ice PSD intercepts
     &  )
!
      Implicit None
!
! Purpose:
!   Perform initialization of variables required for the
!   microphysics scheme.
!
! Method:
!   Set variables to their initial values.
!
! Current Owner of Code: Jonathan Wilkinson
!
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
!
      Integer, Intent(In) ::                                            &
     &  points         ! Number of points
!
      Real, Intent(In) ::                                               &
     &  timestepfixed                                                   &
                       ! Physics timestep / s
     &, T(points)                                                       &
                       ! Temperature / K
     &, p(points)                                                       &
                       ! Pressure / N m-2
     &, q(points)                                                       &
                       ! Vapour content / kg kg-1
     &, qcl(points)                                                     &
                       ! Liquid water content / kg kg-1
     &, qcf(points)                                                     &
                       ! Ice aggregate content / kg kg-1
     &, qcf2(points)                                                    &
                       ! Ice crystal content / kg kg-1
     &, cttemp(points)                                                  &
                       ! Cloud top temperature / K
     &, rainrate(points)                                                &
                         ! Rainrate into layer / kg m-2 s-1
     &, deltaz(points)                                                  &
                       ! Layer thickness / m
     &, rhodz(points)                                                   &
                       ! Air density * deltaz based on hydrostatic
                       ! assumption / kg m-2
     &, rhodz_dry(points)                                               &
                           ! Air density*deltaz for dry air based on
                           ! non-hydrostatic assumption / kg m-2
     &, rhodz_moist(points)! Air density*deltaz for moist air based on
!

      Real, Intent(InOut) ::                                            &
     &  qgraup(points)                                                  &
                       ! Graupel content / kg kg-1
     &, snow_cry(points)                                                &
                          ! Ice crystal precip into layer / kg m-2 s-1
     &, snow_agg(points)  ! Ice agg. precip into layer / kg m-2 s-1

      Real, Intent(InOut) ::                                              &
        ! Microphysical process rate diagnostics / kg kg-1 s-1
     &  psdep(points)                                                   &
                       ! Deposition of vapour to snow aggregates
     &, psaut(points)                                                   &
                       ! Autoconversion of aggregates from crystals
     &, psacw(points)                                                   &
                       ! Accretion of liq. water by snow aggregates
     &, psacr(points)                                                   &
                       ! Collection of rain by snow aggregates
     &, psaci(points)                                                   &
                       ! Collection of ice crystals by aggregates
     &, psmlt(points)                                                   &
                       ! Melting of snow aggregates
     &, psmltevp(points)                                                &
                          ! Evaporation of melting aggregates
     &, psfall(points)                                                  &
                       ! Fall-out of snow aggregates

     &, praut(points)                                                   &
                       ! Autoconversion of cloud drops to rain
     &, pracw(points)                                                   &
                       ! Accretion of liq. water by rain
     &, prevp(points)                                                   &
                       ! Evaporation of rain
     &, prfall(points)                                                  &
                       ! Fall-out of rain

     &, plset(points)                                                   &
                       ! Droplet settling of liquid water
     &, plevpset(points)                                                &
                          ! Evaporated settled droplets

     &, pgaut(points)                                                   &
                       ! Autoconversion of graupel from aggregates
     &, pgacw(points)                                                   &
                       ! Accretion of liq. water by graupel
     &, pgacs(points)                                                   &
                       ! Collection of snow aggregates by graupel
     &, pgmlt(points)                                                   &
                       ! Melting of graupel
     &, pgfall(points)                                                  &
                       ! Fall-out of graupel

     &, pifrw(points)                                                   &
                       ! Homogeneous freezing nucleation
     &, piprm(points)                                                   &
                       ! Heterogeneous (primary) nucleation
     &, pidep(points)                                                   &
                       ! Deposition of vapour to ice crystals
     &, piacw(points)                                                   &
                       ! Accretion of liq. water by ice crystals
     &, piacr(points)                                                   &
                       ! Collection of rain by ice crystals
     &, pimlt(points)                                                   &
                       ! Melting of ice crystals
     &, pimltevp(points)                                                &
                          ! Evaporation of melting ice crystals
     &, pifall(points)
                       ! Fall-out of ice crystals
!
      Real, Intent(Out) ::                                              &
     &  timestep                                                        &
                          ! Timestep of each iteration / s
     &, snowt_cry(points)                                               &
                          ! Ice crystal precip out of layer / kg m-2 s-1
     &, snowt_agg(points)                                               &
                          ! Ice agg. precip out of layer / kg m-2 s-1
     &, rainratet(points)                                               &
                          ! Rain precip rate out of layer / kg m-2 s-1
     &, graupratet(points)                                              &
                          ! Graupel precip rate out of layer/ kg m-2 s-1
     &, qcf_agg(points)                                                 &
                          ! Ice aggregate mixing ratio / kg kg-1
     &, qcf_cry(points)                                                 &
                          ! Ice crystal mixing ratio / kg kg-1
     &, qcf_tot(points)                                                 &
                          ! Total ice content / kg kg-1
     &, frac_agg(points)                                                &
                          ! Fraction of ice that is aggregates
     &, Frac_cry_dep(points)                                            &
                          ! Fraction of supersaturation that can be
                          ! removed by crystals
     &, Frac_agg_dep(points)                                            &
                          ! Fraction of supersaturation that can be
                          !removed by aggregates
     &, qrain(points)                                                   &
                          ! Rain water content / kg kg-1
     &, qs(points)                                                      &
                          ! Saturated humidity wrt ice / kg kg-1
     &, qsl(points)                                                     &
                          ! Saturated humidity wrt liquid / kg kg-1
     &, rho(points)                                                     &
                       ! Air density / kg m-3
     &, rhor(points)                                                    &
                       ! 1 / air density / m kg-1
     &, esi(points)                                                     &
                       ! Vapour pressure wrt ice / N m-2
     &, esw(points)                                                     &
                       ! Vapour pressure wrt liquid / N m-2
     &, lheat_correc_liq(points)                                        &
                                 ! Latent heat correction for
                                 ! liquid (no units)
     &, lheat_correc_ice(points)                                        &
                                 ! Latent heat correction for
                                 ! ice (no units)
     &, dhi(points)                                                     &
                          ! Timestep / deltaz / s m-1
     &, dhir(points)                                                    &
                          ! deltaz / timestep / m s-1
     &, dhilsiterr(points)                                              &
                          ! deltaz / (timestep*lsiter) / m s-1
     &, corr(points)                                                    &
                          ! Fall speed correction factor (no units)
     &, corr2(points)                                                   &
                          ! Diffusivity correction factor (no units)
     &, rocor(points)                                                   &
                          ! sqrt(rho*corr*corr2) (no units)
     &, tcg(points)                                                     &
                          ! Temperature dependent aggregate PSD
                          ! intercept factor (no units)
     &, tcgi(points)                                                    &
                          ! 1 / tcg (no units)
     &, tcgc(points)                                                    &
                          ! Temperature dependent crystal PSD
                          ! intercept factor (no units)
     &, tcgci(points)                                                   &
                          ! 1 / tcgc (no units)
     &, tcgg(points)                                                    &
                          ! Temperature dependent graupel PSD
                          ! intercept factor (no units)
     &, tcggi(points)                                                   &
                          ! 1 / tcgg (no units)
     &, cf_transfer_diag(points)                                        &
                          ! Diagnostic for change in cloud frac
     &, cfl_transfer_diag(points)                                       &
                          ! Diagnostic for change in liq cloud frac
     &, cff_transfer_diag(points)                                       &
                          ! Diagnostic for change in ice cloud frac
     &, rf_transfer_diag(points)
                          ! Diagnostic for change in rain fraction

      Logical, Intent(In) ::                                            &
     &  l_mcr_qcf2                                                      &
                       ! Use second ice prognostic
     &, l_mcr_qrain                                                     &
                       ! Use prognostic rain
     &, l_mcr_qgraup                                                    &
                       ! Use prognostic graupel
     &, l_psd                                                           &
                       ! Use generic ice particle size distribution
     &, l_non_hydrostatic                                               &
                           ! Use non-hydrostatic method of calculating
                           ! level separation
     &, l_mixing_ratio                                                  &
                       ! Use mixing ratio formulation
     &, L_cry_agg_dep
                       ! Limit the supersaturation that can be removed
                       ! by deposition depending on amount of ice in
                       ! each category

      Logical, Intent(Out) ::                                           &
     &  l_crystals
                       ! Two ice categories exist

#include "c_0_dg_c.h"
#include "c_lspmic.h"
#include "c_r_cp.h"
#include "c_epslon.h"
#include "c_lheat.h"
! Array sizes for cx and constp
#include "c_lspsiz.h"
!
! Local Variables
!
      Integer                                                           &
     &  i              ! Loop counter

      Real                                                              &
     &  qc_tot(points) ! Total condensate / kg kg-1

      Real, Parameter ::                                                &
     &  one_over_zerodegc = 1.0/zerodegc                                &
     &, rho1 = 1.0                                                      &
                                     ! Ref. air density / kg m-3
     &, one_over_epsilon = 1.0/epsilon

      !-----------------------------------------------
      ! Set up short timestep
      !-----------------------------------------------
      timestep = timestepfixed/lsiter

      Do i=1,points

      !-----------------------------------------------
      ! Set fluxes to zero
      !-----------------------------------------------
        snowt_cry(i)  = 0.0
        snowt_agg(i)  = 0.0
        rainratet(i)  = 0.0
        graupratet(i) = 0.0

      !-----------------------------------------------
      ! Initialize transfer diagnostics to zero
      !-----------------------------------------------
        pifrw(i) = 0.0
        piprm(i) = 0.0
        pidep(i) = 0.0
        piacw(i) = 0.0
        piacr(i) = 0.0
        pimlt(i) = 0.0
        pimltevp(i) = 0.0
        pifall(i)= 0.0

        psdep(i) = 0.0
        psaut(i) = 0.0
        psacw(i) = 0.0
        psacr(i) = 0.0
        psaci(i) = 0.0
        psmlt(i) = 0.0
        psmltevp(i) = 0.0
        psfall(i)= 0.0

        praut(i) = 0.0
        pracw(i) = 0.0
        prevp(i) = 0.0
        prfall(i)= 0.0

        plset(i) = 0.0
        plevpset(i) = 0.0

        pgaut(i) = 0.0
        pgacw(i) = 0.0
        pgacs(i) = 0.0
        pgmlt(i) = 0.0
        pgfall(i)= 0.0

        cf_transfer_diag(i) = 0.0
        cfl_transfer_diag(i)= 0.0
        cff_transfer_diag(i)= 0.0
        rf_transfer_diag(i) = 0.0

      End Do

      !-----------------------------------------------
      ! Set mixing ratios of graupel to zero if not used.
      ! If not done then (with high optimisation) this can be any
      ! old rubbish from memory and you get d_qgraup_dt being
      ! calculated as nan-nan in later routines, which causes a crash
      !-----------------------------------------------
      If (.not. l_mcr_qgraup) Then
        Do i=1,points
          qgraup(i)   = 0.0
        End Do
      End If  ! .not. l_mcr_qgraup

      !-----------------------------------------------
      ! Split ice crystal and snow aggregate categories
      !-----------------------------------------------
      l_crystals = l_mcr_qcf2 .or. (.not. l_psd)

      If (L_mcr_qcf2) then
        ! If L_mcr_qcf2 is true then there are two ice/snow
        ! prognostics (qcf and qcf2), so copy them to
        ! qcf_cry and qcf_agg

        Do i=1,points
          qcf_cry(i) = qcf2(i)
          qcf_agg(i) = qcf(i)
        End Do

      Else if (.not. l_crystals) then
        ! All ice is placed in a single ice category (aggregates).

        Do i=1,points
          frac_agg(i) = 1.0
          qcf_cry(i)  = 0.0
          qcf_agg(i)  = qcf(i)
          snow_cry(i) = 0.0
        End Do

      Else
        ! Split the one ice/snow
        ! prognostic (qcf) diagnostically into ice crystals
        ! (qcf_cry) and snow aggregates (qcf_agg)

        Do i=1,points

#if defined(VECTLIB)
          frac_agg(i) = -T_scaling*max((T(i)-cttemp(i)),0.0)            &
     &                  *max(qcf(i)*qcf0,0.0)
        End Do
! DEPENDS ON: exp_v
        call exp_v(points,frac_agg,frac_agg)
        Do i=1,points
          frac_agg(i) = max(1.0-frac_agg(i) , 0.0)
#else
          frac_agg(i) = max(1.0-exp(-T_scaling                          &
     &                              *max((T(i)-cttemp(i)),0.0)          &
     &                * max(qcf(i),0.0)*qcf0) , 0.0)
#endif

          ! Allocate ice content to crystals and aggregates
          qcf_cry(i) = qcf(i) * (1.0-frac_agg(i))
          qcf_agg(i) = qcf(i) * frac_agg(i)

          ! Assume falling snow is partitioned into crystals and
          ! aggregates. Snow_agg contains total snow on input
          snow_cry(i) = snow_agg(i) * (1.0-frac_agg(i))
          snow_agg(i) = snow_agg(i) * frac_agg(i)

        End Do

      End If ! L_mcr_qcf2

      !-----------------------------------------------
      ! Calculate total ice content
      !-----------------------------------------------
      Do i=1,points
        qcf_tot(i) = qcf_cry(i) + qcf_agg(i)
      End Do

      !-----------------------------------------------
      ! Calculate the fraction of ice in the crystals and
      ! aggregates category that is allowed to remove
      ! supersaturation
      !-----------------------------------------------
      If (l_cry_agg_dep) then
        If (l_mcr_qcf2) then
          ! Use the partition given by the prognostic ice categories.
          ! Here we add a small amount of ice to the crystal category
          ! to ensure no problems when the ice contents are zero.
          Do i=1,points
            frac_cry_dep(i) = (qcf_cry(i)+m0)/(Max(qcf_tot(i),0.0)+m0)
            frac_agg_dep(i) = qcf_agg(i)/(Max(qcf_tot(i),0.0)+m0)
          End do
        Else  ! l_mcr_qcf2
          ! Use the diagnostic partition function
          Do i=1,points
            frac_cry_dep(i)=1.0-frac_agg(i)
            frac_agg_dep(i)=frac_agg(i)
          End do
        End if
      Else
        ! Set the fractions to 1 to maintain bit reproducibility
        Do i=1,points
          frac_cry_dep(i)=1.0
          frac_agg_dep(i)=1.0
        End do
      End if

      !-----------------------------------------------
      ! If rain is a diagnostic, convert flux (kg m-2 s-1)
      ! to mass (kg kg-1)
      !-----------------------------------------------
      If (.not. l_mcr_qrain) Then

        ! Rain is diagnostic
        Do i=1,points
          If (rainrate(i)  >   0.0) then
            If (L_non_hydrostatic) Then
              If (L_mixing_ratio) Then
                ! Mixing ratio formulation
                qrain(i) = rainrate(i) * timestepfixed/rhodz_dry(i)
              Else
                ! Specific humidity formulation
                qrain(i) = rainrate(i) * timestepfixed/rhodz_moist(i)
              End if  ! l_mixing_ratio
            Else
              ! Use hydrostatic formulation
              qrain(i) = rainrate(i) * timestepfixed/rhodz(i)
            End If  ! l_non_hydrostatic
          Else
            qrain(i) = 0.0
          End If  ! rainrate ne 0
        End Do

      End If  ! .not. l_mcr_qrain

      !-----------------------------------------------
      ! Calculate saturation specific humidities
      !-----------------------------------------------
      ! Qsat with respect to ice
! DEPENDS ON: qsat_mix
      Call qsat_mix(qs,T,p,points,l_mixing_ratio)

      ! Qsat with respect to liquid water
! DEPENDS ON: qsat_wat_mix
      Call qsat_wat_mix(qsl,T,p,points,l_mixing_ratio)

      Do i=1,points
      !-----------------------------------------------
      ! Calculate saturation vapour pressures
      !-----------------------------------------------
        esi(i) = qs (i)*p(i)*one_over_epsilon
        esw(i) = qsl(i)*p(i)*one_over_epsilon

        !-----------------------------------------------
        ! Calculate density of air
        !-----------------------------------------------
        If (l_non_hydrostatic) then
          If (l_mixing_ratio) then
            ! rho is the dry density
            rho(i) = rhodz_dry(i) / deltaz(i)
          Else
            ! rho is the moist density
            rho(i) = rhodz_moist(i) / deltaz(i)
          End if  ! l_mixing_ratio
        Else
          ! Use the pressure to retrieve the moist density.
          ! An exact expression for the air density is
          ! rho = p / (1+c_virtual q - qcl - qcf) / RT but
          ! the approximation below is more numerically
          ! stable.
          !
          ! First calculate total condensate mass
          qc_tot(i) = qcl(i)+qcf_tot(i)
          If (l_mcr_qrain) then
            qc_tot(i) = qc_tot(i) + qrain(i)
          End If
          If (l_mcr_qgraup) then
            qc_tot(i) = qc_tot(i) + qgraup(i)
          End If

          rho(i) = p(i)*(1.-c_virtual*q(i)+qc_tot(i)) / (R*T(i))
        End if  ! l_non_hydrostatic

        ! Calculate the inverse of the air density
        rhor(i) = 1.0/rho(i)

        !-----------------------------------------------
        ! Estimate latent heat correction to rate of evaporation
        !-----------------------------------------------
        lheat_correc_liq(i) = 1.0/(1.0+epsilon*lc**2*qsl(i)             &
     &                           /(cp*R*T(i)**2))
        lheat_correc_ice(i) = 1.0/(1.0+epsilon*(lc+lf)**2*qs(i)         &
     &                           /(cp*R*T(i)**2))

        !-----------------------------------------------
        ! Calculate CFL timestep divided by level separation
        !-----------------------------------------------
        If (l_non_hydrostatic) Then
          ! Use the formulation based on the heights of the levels
          dhi(i) = timestep/deltaz(i)
        Else
          ! Use the formulation based on the pressure difference
          ! across a layer.
          dhi(i)        = timestep * rho(i)/rhodz(i)
        End if

        ! Calculate the inverse
        dhir(i)       = 1.0/dhi(i)
        ! Calculate the inverse for the long timestep
        dhilsiterr(i) = 1.0/(dhi(i)*lsiter)

      End Do

      !-----------------------------------------------
      ! Correction factors due to air density and temperature
      !-----------------------------------------------
#if defined(VECTLIB)
      Do i=1, points

        ! Correction of fall speeds
        If (l_mixing_ratio .and. l_non_hydrostatic) then
          corr(i) = rho1*deltaz(i) / rhodz_moist(i)
        Else
          corr(i) = rho1*rhor(i)
        End if

        ! Correction factor in viscosity etc. due to temperature
        corr2(i) = (T(i)*one_over_zerodegc)

      End Do

      Call powr_v(points,corr,0.4,corr)
      Call powr_v(points,corr2,1.5,corr2)

      Do i=1,points
        corr2(i) = corr2(i)*(393.0/(T(i)+120.0))

        ! Combined correction factor
        If (l_mixing_ratio .and. l_non_hydrostatic) then
          rocor(i) = rhodz_moist(i) / deltaz(i) * corr(i) * corr2(i)
        Else
          rocor(i) = rho(i)*corr(i)*corr2(i)
        End if

      End Do

      Call powr_v(points,rocor,0.5,rocor)
#else
      Do i=1,points

        ! Correction of fall speeds
        If (l_mixing_ratio .and. l_non_hydrostatic) then
          corr(i) = (rho1*deltaz(i) / rhodz_moist(i))**0.4
        Else
          corr(i) = (rho1*rhor(i))**0.4
        End if

        ! Correction factor in viscosity etc. due to temperature
        corr2(i)=(T(i)*one_over_zerodegc)**1.5 * (393.0/(T(i)+120.0))

        ! Combined correction factor
        If (l_mixing_ratio .and. l_non_hydrostatic) then
          rocor(i) = sqrt(rhodz_moist(i)/deltaz(i)*corr(i)*corr2(i))
        Else
          rocor(i) = sqrt(rho(i)*corr(i)*corr2(i))
        End if

      End Do
#endif

      !-----------------------------------------------
      ! Calculate ice particle size distributions
      !-----------------------------------------------
#if defined(VECTLIB)
      Do i=1,points

        ! Calculate a temperature factor for N0crystals
        tcgc(i) = -cx(12)*max(T(i)-zerodegc,T_agg_min)
        ! Calculate a temperature factor for N0aggregates
        tcg(i)  = -cx(32)*max(T(i)-zerodegc,T_agg_min)

        ! Define inverse of TCG values
        tcgci(i) = -tcgc(i)
        tcgi(i)  = -tcg(i)

      End Do

! DEPENDS ON: exp_v
      Call exp_v(points,tcg,tcg)
! DEPENDS ON: exp_v
      Call exp_v(points,tcgc,tcgc)
! DEPENDS ON: exp_v
      Call exp_v(points,tcgi,tcgi)
! DEPENDS ON: exp_v
      Call exp_v(points,tcgci,tcgci)

#else
      Do i=1,points

        ! Calculate a temperature factor for N0crystals
        tcgc(i) = exp(-cx(12)*max(T(i)-zerodegc,T_agg_min))
        ! Calculate a temperature factor for N0aggregates
        tcg(i)  = exp(-cx(32)*max(T(i)-zerodegc,T_agg_min))

        ! Define inverse of TCG values
        tcgci(i) = 1.0/tcgc(i)
        tcgi(i)  = 1.0/tcg(i)

      End Do
#endif
      !-----------------------------------------------
      ! Calculate graupel size distributions
      !-----------------------------------------------
      Do i=1,points
        tcgg(i)  = 1.0
        tcggi(i) = 1.0/tcgg(i)
      End Do

      Return  ! End of the subroutine
      END SUBROUTINE LSP_INIT
#endif
