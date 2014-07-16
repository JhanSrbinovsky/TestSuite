#if defined(A04_3D)
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
#include "c_lspsiz.h"
#include "c_lspdif.h"
#include "c_lspmic.h"
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
#endif
