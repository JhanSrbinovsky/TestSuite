#if defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Accretion of droplets by raindrops
! Subroutine Interface:
      SUBROUTINE LSP_ACCRETION(                                         &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, qcl, qrain                                                      &
                                          ! Water contents
!    &, area_liq, area_mix, area_ice      ! Partition information
     &, cfliq                                                           &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
     &, rainfrac, rain_liq, rain_mix                                    &
                                          ! Rain fraction information
!    &, rain_ice, rain_clear
!    &, cf, cfl                           ! Current cloud fractions for
!                                         ! updating
     &, rho, corr, dhilsiterr                                           &
                                          ! Parametrization information
     &, cx, constp                                                      &
                                          ! Microphysical information
     &, l_update_cf, l_seq, l_mcr_qrain                                 &
                                          ! Code options
     &, ptransfer, iterations                                           &
                                          ! Mass transfer diagnostic
!    &, cftransfer,cfltransfer            ! Cloud transfer diagnostics
!    &, rftransfer                        ! Rain fraction transfer
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of accretion of cloud
!   droplets by rain drops
!
! Method:
!   Solve implicitly the microphysical transfer equation for a
!   specified distribution of raindrops sweeping out a volume
!   of air uniformally inhabited by cloud water droplets.
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
#include "c_lspsiz.h"
!
      Integer, Intent(In) ::                                            &
     &  points                                                          &
                          ! Number of points to calculate
     &, iterations        ! Number of microphysics iterations
                          ! (for rate diagnostic)
!
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, rain_liq(points)                                                &
                          ! Overlap fraction of rain and liquid
     &, rain_mix(points)                                                &
                          ! Overlap fraction of rain and mixed phase
     &, rainfrac(points)                                                &
                          ! Rain fraction
     &, cfliq(points)                                                   &
                          ! Liquid cloud fraction at start of timestep
!    &, cf(points)        ! Current cloud fraction
!    &, cfl(points)       ! Current liquid cloud fraction
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, corr(points)                                                    &
                          ! Fall velocity correction factor (no units)
     &, dhilsiterr(points)! Depth of layer / timestep  / m s-1
!
      Real, Intent(InOut) ::                                            &
     &  qcl(points)                                                     &
                          ! Liquid water content / kg kg-1
     &, qrain(points)                                                   &
                          ! Rain mixing ratio / kg kg-1
!    &, rain_ice(points)  ! Overlap fraction of rain and ice
!    &, rain_clear(points)! Overlap fraction of rain and clear sky
     &, ptransfer(points) ! Mass rimed in this timestep / kg kg-1
!    &, rftransfer(points)! Transfer of rain fraction this timestep
!    &, cftransfer(points) ! Cumulative cloud fraction increment
!    &, cfltransfer(points)! Cumulative liquid cloud fraction increment
!
      Logical, Intent(In) ::                                            &
     &  l_update_cf                                                     &
                          ! Update cloud fractions
     &, l_seq                                                           &
                          ! Carry out sequential updating
     &, l_mcr_qrain       ! Use prognostic rain mixing ratio
!
! Local Variables
!
      Integer                                                           &
     &  i                 ! Loop counter for points
!
      Real                                                              &
     &  dpr(points)                                                     &
                          ! Transfer of mixing ratio  / kg kg-1
     &, temp1(points)                                                   &
                          ! Temporary in rate calculation
     &, temp7(points)                                                   &
                          ! Rain free liquid cloud fraction (no units)
     &, qclnew(points)    ! Updated value of liquid water / kg kg-1

      Do i=1,points

        If (rainfrac(i) >  0.0 .and. qrain(i) >  0.0                    &
     &      .and. cfliq(i) >  0.0) then

          If (L_mcr_qrain) Then
            ! Use the prognostic rain formulation
            temp1(i) = qrain(i)*rho(i)*constp(50)/(rainfrac(i))
          Else
            ! Use the diagnostic rain formulation
            temp1(i) = qrain(i)*rho(i)*dhilsiterr(i)                    &
     &                 /(rainfrac(i)*constp(42)*corr(i))
          End If  ! L_mcr_qrain

          !-----------------------------------------------
          ! Calculate the new local value of qcl
          !-----------------------------------------------
          If (L_mcr_qrain) Then
            qclnew(i)=qcl(i)/((cfliq(i)+cfliq(i)*constp(49)*corr(i)     &
     &                               *timestep*temp1(i)**cx(53)))
          Else
            qclnew(i)=qcl(i)/((cfliq(i)+cfliq(i)*constp(49)*corr(i)     &
     &                               *timestep*temp1(i)**cx(50)))
          End if

          !-----------------------------------------------
          ! Convert qclnew to a gridbox mean
          !-----------------------------------------------
          ! temp7 is the rain free region of liquid cloud
          temp7(i) = max(cfliq(i)-rain_liq(i)-rain_mix(i),0.0)
          qclnew(i) = qclnew(i)*(rain_liq(i)+rain_mix(i))               &
     &               +qcl(i)/cfliq(i)*temp7(i)

          !-----------------------------------------------
          ! Calculate transfer
          !-----------------------------------------------
          dpr(i) = qcl(i) - qclnew(i)
          ptransfer(i) = ptransfer(i) + dpr(i) / (timestep*iterations)

          !------------------------------------------------
          ! Adjust liquid and rain contents
          !------------------------------------------------
          If (L_seq) Then
            qrain(i) = qrain(i) + dpr(i)
            qcl(i)   = qcl(i)   - dpr(i)
          End If

          !------------------------------------------------
          ! Update cloud fractions
          !------------------------------------------------
!         These are commented out since there is currently no
!         cloud fraction update associated with the accretion term.
!          cf_transfer_rate(i)  = 0.0 / (timestep*iterations)
!          cfl_transfer_rate(i) = 0.0 / (timestep*iterations)
!
!          If (l_update_cf) then
!            cftransfer(i)  = cftransfer(i)  + cf_transfer_rate(i)
!            cfltransfer(i) = cfltransfer(i) + cfl_transfer_rate(i)
!
!            If (l_seq) then
!              cf(i)  = cf(i)  +cf_transfer_rate(i) *timestep*iterations
!              cfl(i) = cfl(i) +cfl_transfer_rate(i)*timestep*iterations
!            End if  ! l_seq
!
!          End if  ! l_update_cf

          !-----------------------------------------------
          ! Update rain fractions
          !-----------------------------------------------
!         These are commented out since there is currently no
!         rain fraction update associated with the accretion term.
!          rf_final(i) = rainfrac(i)
!          rftransfer(i) = rftransfer(i) + (rf_final(i) - rainfrac(i))
!
!          If (L_seq) then
!            rainfrac(i)= rf_final(i)
!            rain_liq(i)= min(area_liq(i),rainfrac(i))
!            rain_mix(i)= min(area_mix(i),rainfrac(i)-rain_liq(i))
!            rain_ice(i)=
!     &         min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
!            rain_clear(i)=
!     &           rainfrac(i) - rain_liq(i)-rain_mix(i)-rain_ice(i)
!          End if  ! L_seq

        End If  ! rainfrac gt 0 etc.

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_ACCRETION
#endif
