#if defined(A04_3D)
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

#if defined(ACCESS)
!dhb599 20110601: new cloud fraction increment scaling factor:
      use auscom_cpl_data_mod, only : xfactor
#endif
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
#if defined(ACCESS)
!dhb599 20110601--scaling (down) this new cloud fraction increment as per fra298
!(in the purose of cooling down the system by reducing high cloud coverage)
!              write(6,*)'XXX lsp_deposition: xfactor = ', xfactor
              deltacf(i) = deltacf(i) * xfactor
#endif
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
#endif
