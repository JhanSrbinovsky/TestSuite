#if defined(A04_3D)
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
#include "c_lspsiz.h"
#include "c_0_dg_c.h"
#include "c_lspmic.h"
#include "c_lspdif.h"
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
#endif
