#if defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Iterative melting interface
! Subroutine Interface:
      SUBROUTINE LSP_IT_FALL_MELT(                                      &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, T, p, q, q_ice, qsl                                             &
                                          ! Temperature and vapour
     &, qcf_cry, qcf_agg, frac_agg                                      &
                                          ! Ice contents
     &, qrain, qgraup                                                   &
                                          ! Rain and graupel contents
     &, snow_agg,snow_cry,rainrate,grauprate                            &
                                             ! Sedimentation into layer
     &, snowt_agg,snowt_cry,rainratet,graupratet                        &
                                                 ! Sedim. out of layer
     &, vf_agg, vf_cry, vf_rain, vf_graup                               &
                                          ! Fall speeds of hydrometeors
     &, area_liq, area_mix, area_clear                                  &
                                          ! Cloud fraction information
     &, area_ice,cfice,cficei                                           &
                                          ! at start of microphy. ts
     &, frac_ice_above                                                  &
     &, cf, cfl, cff                                                    &
                                          ! Current cloud fractions for
                                          ! updating
     &, rainfrac, rain_liq, rain_mix                                    &
                                          ! Rain fraction information
     &, rain_ice, rain_clear                                            &
     &, rho, rhor                                                       &
     &, tcg,tcgi,tcgc,tcgci,tcgg,tcggi                                  &
                                          ! Parametrization information
     &, corr, corr2, rocor, dhi, dhir                                   &
                                          ! Parametrization information
     &, cx, constp, ai, bi, aic, bic, lfrcp                             &
                                          ! Microphysical information
     &, l_update_cf, l_seq                                              &
                                          ! Code options
     &, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_psd, l_crystals        &
                                              ! Microphysics options
     &, pifall, psfall, prfall, pgfall                                  &
                                          ! Mass transfer diagnostics
     &, pimlt, psmlt, pgmlt                                             &
                                          ! Mass transfer diagnostics
     &, iterations, n_melt_iterations                                   &
                                          ! Iterations of microphysics
     &, cftransfer,cfftransfer,rftransfer                               &
                                          ! Cloud transfer diagnostics
     &  )
!
      Implicit None
!
! Purpose:
!   Call fall out and melting terms in an iterative way
!
! Method:
!   Iterate over the fall-out and melting terms for a prescribed
!   number of iterations.
!
! Current Owner of Code: Jonathan Wilkinson
!
! History:
! Version   Date     Comment
!  6.2    19-12-05   Original Code (Damian Wilson)
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
#include "c_lspsiz.h"
#include "c_lspmic.h"
!
      Integer, Intent(In) ::                                            &
     &  points                                                          &
                          ! Number of points to calculate
     &, iterations                                                      &
                          ! Number of microphysics iterations
     &, n_melt_iterations ! Number of iterative melting iterations
!
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, q(points)                                                       &
                          ! Vapour content / kg kg-1
     &, q_ice(points)                                                   &
                          ! Vapour content in ice only region / kg kg-1
     &, qsl(points)                                                     &
                          ! Sat. vapour content wrt liquid / kg kg-1
     &, p(points)                                                       &
                          ! Pressure / N m-2
     &, frac_agg(points)                                                &
                          ! Fraction of ice mass that is aggregates
     &, snow_cry(points)                                                &
                          ! Crystal flux into layer / kg m-2 s-1
     &, snow_agg(points)                                                &
                          ! Aggregate flux into layer / kg m-2 s-1
     &, rainrate(points)                                                &
                          ! Rain flux into layer / kg m-2 s-1
     &, grauprate(points)                                               &
                          ! Graupel flux into layer / kg m-2 s-1
     &, area_liq(points)                                                &
                          ! Frac of gridbox with liq cloud but not ice
     &, area_mix(points)                                                &
                          ! Fraction of gridbox with ice and liq cloud
     &, area_clear(points)                                              &
                          ! Fraction of gridbox with clear sky
     &, area_ice(points)                                                &
                          ! Frac of gridbox with ice cloud but not liq
     &, cfice(points)                                                   &
                          ! Fraction of gridbox with ice cloud
     &, cficei(points)                                                  &
                          ! 1/Fraction of gridbox with ice cloud
     &, frac_ice_above(points)                                          &
                               ! Ice cloud fraction in layer above
     &, rainfrac(points)                                                &
                          ! Fraction of gridbox with rain
     &, rain_liq(points)                                                &
                          ! Fraction of gridbox with rain and liq cloud
     &, rain_mix(points)                                                &
                          ! Frac of gbox with rain and mixed phase cloud
     &, rain_ice(points)                                                &
                          ! Fraction of gridbox with rain and ice cloud
     &, rain_clear(points)                                              &
                          ! Fraction of gridbox with rain but no cloud
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, rhor(points)                                                    &
                          ! 1 / Air density / m3 kg-1
     &, dhi(points)                                                     &
                          ! Timestep / thickness of model layer / s m-1
     &, dhir(points)                                                    &
                          ! 1/dhi / m s-1
!    &, m0                ! Seed ice water / kg kg-1 (in c_lspmic)
!    &, wind_shear_factor ! Degree of overhang of ice fraction between
                          ! two model layers (in c_lspmic)
     &, tcg(points)                                                     & 
                          ! T dependent func in ice agg. size dist'n
     &, tcgi(points)                                                    &
                          ! 1/tcg (no units)
     &, tcgc(points)                                                    & 
                          ! T dependent func in ice crystal size dist'n
     &, tcgci(points)                                                   &
                          ! 1/tcgc (no units)
     &, tcgg(points)                                                    & 
                          ! T dependent func in graupel size dist'n
     &, tcggi(points)                                                   &
                          ! 1/tcgg (no units)
     &, corr(points)                                                    &
                          ! Air density fall speed correction (no units)
     &, corr2(points)                                                   &
                          ! Temperature correcn to diffusivity (no units
     &, rocor(points)                                                   &
                          ! sqrt(rho*corr*corr2)
     &, lfrcp                                                           &
                          ! Latent heat of fusion
                          !   / heat capacity of air / K
     &, ai, bi, aic, bic
                          ! Ice mass size relationships. m(D)=ai D^bi
!
      Real, Intent(InOut) ::                                            &
     &  T(points)                                                       &
                          ! Temperature / K
     &, qcf_cry(points)                                                 &
                          ! Ice crystal mixing ratio / kg kg-1
     &, qcf_agg(points)                                                 &
                          ! Ice aggregate mixing ratio / kg kg-1
     &, qrain(points)                                                   &
                          ! Rain mixing ratio / kg kg-1
     &, qgraup(points)                                                  &
                          ! Graupel mixing ratio / kg kg-1
     &, snowt_cry(points)                                               &
                          ! Snowfall rate out of this layer / kg m-2 s-1
     &, snowt_agg(points)                                               &
                            ! for crystals and aggregates
     &, rainratet(points)                                               &
                          ! Rain rate out of this layer / kg m-2 s-1
     &, graupratet(points)                                              &
                          ! Graupel rate out of this layer / kg m-2 s-1
     &, vf_cry(points)                                                  &
                          ! On input: Fall speed of hydrometeors
     &, vf_agg(points)                                                  &
                                    ! entering the current layer / m s-1
     &, vf_rain(points)                                                 &
                          ! On output: Fall speed of hydrometeors
     &, vf_graup(points)                                                &
                                    ! leaving the current layer / m s-1
     &, cf(points)                                                      &
                          ! Current cloud fraction
     &, cfl(points)                                                     &
                          ! Current liquid cloud fraction
     &, cff(points)       ! Current ice cloud fraction

      Real, Intent(InOut) ::                                            &
     &  pifall(points)                                                  &
                          ! Rate of change of crystal, aggregate,
     &, psfall(points)                                                  &
                          ! rain and graupel mixing ratios due
     &, prfall(points)                                                  &
                          ! to sedimentation / kg kg-1 s-1
     &, pgfall(points)                                                  &
     &, pimlt(points)                                                   &
                          ! Melting rate of crystals / kg kg-1 s-1
     &, psmlt(points)                                                   &
                          ! Melting rate of snow / kg kg-1 s-1
     &, pgmlt(points)                                                   &
                          ! Melting rate of graupel / kg kg-1 s-1
     &, cftransfer(points)                                              &
                          ! Cloud fraction increment this tstep
     &, cfftransfer(points)                                             &
                           ! Ice cloud fraction inc this tstep
     &, rftransfer(points)! Rain fraction increment this tstep
!
      Logical, Intent(In) ::                                            &
     &  l_update_cf                                                     &
                          ! Update cloud fractions
     &, l_seq                                                           &
                          ! Carry out sequential updating
     &, l_mcr_qcf2                                                      &
                          ! Two ice prognostics are active
     &, l_mcr_qrain                                                     &
                          ! Prognostic rain is active
     &, l_mcr_qgraup                                                    &
                          ! Graupel is active
     &, l_psd                                                           &
                          ! Use generic ice particle size distribution
     &, l_crystals
                          ! Use the ice crystal category
!
! Local Variables
!
      Integer                                                           &
     &  i                                                               &
                          ! Loop counter for points
     &, k                                                               &
                          ! Loop counter for melting iterations
     &, total_iterations  ! Net number of iterations
!
      Real                                                              &
     &  qcft(points)                                                    &
                          ! Total ice content / kg kg-1
     &, dhi_it(points)                                                  &
                          ! Iterated tstep / layer thick. / s m-1
     &, dhir_it(points)                                                 &
                          ! Layer thick. / iterated tstep / m s-1
     &, timestep_it       ! Iterated timestep / s
!
      !-----------------------------------------------
      ! Set up short timestep
      !-----------------------------------------------
      timestep_it = timestep / n_melt_iterations
      total_iterations = iterations * n_melt_iterations

      Do i = 1, points
        ! Update dhi values for short timestep
        dhi_it(i) = dhi(i) / n_melt_iterations
        dhir_it(i) = dhir(i) * n_melt_iterations
      End Do

      !-----------------------------------------------
      ! Call fall and melting for the iterated points
      !-----------------------------------------------
      Do k = 1, n_melt_iterations

! DEPENDS ON: lsp_fall
        Call lsp_fall(points, timestep_it                               &
     &,                  qcf_cry, qcf_agg, frac_agg, qrain, qgraup, T   &
     &,                  snow_agg, snow_cry, rainrate, grauprate        &
     &,                  snowt_agg, snowt_cry, rainratet, graupratet    &
     &,                  vf_agg, vf_cry, vf_rain, vf_graup              &
     &,                  area_clear, area_ice, cfice, cficei            &
     &,                  frac_ice_above, cf, cfl, cff                   &
     &,                  rho, rhor, tcgi, tcgci                         &
     &,                  corr, dhi_it, dhir_it, cx, constp              &
     &,                  ai, bi, aic, bic                               &
     &,                  l_update_cf, l_seq                             &
     &,                  l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup          &
     &,                  l_psd, l_crystals                              &
     &,                  pifall, psfall, prfall, pgfall                 &
     &,                  total_iterations                               &
     &,                  cftransfer, cfftransfer                        &
     &                  )

        Do i=1,points
          ! Calculate total ice content
          qcft(i) = qcf_cry(i) + qcf_agg(i)
        End do

        If (l_crystals) then
! DEPENDS ON: lsp_melting
          Call lsp_melting(points, timestep_it                          &
     &,                   q, q_ice, qcf_cry, qcft, qrain, qsl, T, p     &
     &,                   area_liq, area_mix, area_ice                  &
     &,                   cfice, cficei, frac_ice_above                 &
     &,                   rainfrac, rain_liq, rain_mix                  &
     &,                   rain_ice, rain_clear, cf, cff                 &
     &,                   rho, rhor, m0, tcg, tcgi, corr2, rocor        &
     &,                   cx, constp, ai, bi, aic, bic, lfrcp, 0        &
     &,                   l_update_cf, l_seq, .false.                   &
     &,                   pimlt, total_iterations                       &
     &,                   cftransfer, cfftransfer, rftransfer           &
     &                   )
        End if  ! l_crystals

        Do i=1,points
          ! Calculate total ice content
          qcft(i) = qcf_cry(i) + qcf_agg(i)
        End do

! DEPENDS ON: lsp_melting
        Call lsp_melting(points, timestep_it                            &
     &,                   q, q_ice, qcf_agg, qcft, qrain, qsl, T, p     &
     &,                   area_liq, area_mix, area_ice                  &
     &,                   cfice, cficei, frac_ice_above                 &
     &,                   rainfrac, rain_liq, rain_mix                  &
     &,                   rain_ice, rain_clear, cf, cff                 &
     &,                   rho, rhor, m0, tcgc, tcgci, corr2, rocor      &
     &,                   cx, constp, ai, bi, aic, bic, lfrcp, 1        &
     &,                   l_update_cf, l_seq, l_psd                     &
     &,                   psmlt, total_iterations                       &
     &,                   cftransfer, cfftransfer, rftransfer           &
     &                   )

        If (L_mcr_qgraup) Then
          ! Graupel does not update cloud fractions so there is no need
          ! to update qcft (it is not used)
! DEPENDS ON: lsp_melting
          Call lsp_melting(points, timestep_it                          &
     &,                     q, q_ice, qgraup, qcft, qrain, qsl, T, p    &
     &,                     area_liq, area_mix, area_ice                &
     &,                     cfice, cficei, frac_ice_above               &
     &,                     rainfrac, rain_liq, rain_mix                &
     &,                     rain_ice, rain_clear, cf, cff               &
     &,                     rho, rhor, m0, tcgg, tcggi, corr2, rocor    &
     &,                     cx, constp, ai, bi, aic, bic, lfrcp, 3      &
     &,                     .false., l_seq, .false.                     &
     &,                     pgmlt, total_iterations                     &
     &,                     cftransfer, cfftransfer, rftransfer         &
     &                     )
        End If  ! L_mcr_qgraup

      End Do  ! n_melt_iterations

      Return  ! End of the subroutine
      END SUBROUTINE LSP_IT_FALL_MELT
#endif
