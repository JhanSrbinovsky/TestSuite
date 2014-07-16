#if defined(A04_3D)
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
#include "c_0_dg_c.h"
#include "c_lspmic.h"
#include "c_lspdif.h"
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
#endif
