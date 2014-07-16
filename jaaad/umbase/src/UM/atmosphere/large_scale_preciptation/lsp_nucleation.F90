#if defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Nucleation of ice particles
! Subroutine Interface:
      SUBROUTINE LSP_NUCLEATION(                                        &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, q, qcl, qcf, T                                                  &
                                          ! Water contents, temperature
     &, qs, qsl                                                         &
                                          ! Saturated quantities
     &, cfliq, cfice                                                    &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
     &, area_liq, area_mix, cf, cfl, cff                                &
                                          ! Current cloud fractions for
                                          ! updating
     &, rain_liq, rain_mix                                              &
                                          ! Overlaps of cloud with rain
     &, rhor, lheat_correc_ice                                          &
                                          ! Parametrization information
     &, lfrcp, lsrcp                                                    &
                                          ! Microphysical information
     &, l_update_cf, l_seq                                              &
     &, hettransfer,homtransfer,iterations                              &
                                          ! Mass transfer diagnostics
     &, cftransfer,cfltransfer,cfftransfer                              &
                                          ! Cloud transfer diagnostics
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of homogeneous and
!   heterogeneous ice nucleation
!
! Method:
!   Homogeneous nucleation converts all liquid to ice with a temperature
!   threshold. Heterogeneous nucleation converts a seed amount of liquid
!   or vapour to ice.
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
#include "c_micro.h"
!
      Integer, Intent(In) ::                                            &
     &  points                                                          &
                          ! Number of points to calculate
     &, iterations        ! Number of microphysics iterations
                          ! (for rate diagnostic)
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, qs(points)                                                      &
                          ! Saturated humidity wrt ice / kg kg-1
     &, qsl(points)                                                     &
                          ! Saturated humidity wrt liquid / kg kg-1
     &, cfliq(points)                                                   &
                          ! Fraction of gridbox with liquid cloud
     &, cfice(points)                                                   &
                          ! Fraction of gridbox with ice cloud
     &, area_liq(points)                                                &
                          ! Fraction of gridbox with only liquid cloud
     &, area_mix(points)                                                &
                          ! Fraction of gridbox with mixed phase cloud
     &, rain_liq(points)                                                &
                          ! Frac. of gridbox with rain and liquid cld.
     &, rain_mix(points)                                                &
                          ! Frac. of grbx. with rain and mixed ph. cld
     &, rhor(points)                                                    &
                          ! 1/Air density / kg m-3
     &, lheat_correc_ice(points)                                        &
                                 ! Ice latent heat correction factor
!    &, m0                ! Seed ice water content / kg kg-1 c_lspmic
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
     &, T(points)                                                       &
                          ! Temperature / K
     &, cf(points)                                                      &
                          ! Current cloud fraction
     &, cfl(points)                                                     &
                          ! Current liquid cloud fraction
     &, cff(points)                                                     &
                          ! Current ice cloud fraction
     &, homtransfer(points)                                             &
                           ! Mass homog. nucleated this ts / kg kg-1
     &, hettransfer(points)! Mass het. nucleated this ts / kg kg-1
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
     &, l_seq             ! Carry out sequential updating
!
! Local Variables
!
      Integer                                                           &
     &  i                 ! Loop counter for points

      Real                                                              &
     &  dqi(points)                                                     &
                          ! Temporary in calculating ice transfers
     &, dqil(points)                                                    &
                          ! Temporary in calculating liquid transfers
     &, rhnuc             ! Nucleation relative humidity threshold

      Do i=1, points

        !-----------------------------------------------
        ! Homogeneous nucleation - freeze all liquid if T < threshold
        !-----------------------------------------------
        If (T(i) <  (zerodegc+thomo) .and. qcl(i) >  0.0) Then

          !-----------------------------------------------
          ! Update cloud fractions
          !-----------------------------------------------
          If (l_update_cf) Then
            cfftransfer(i)=cfftransfer(i)                               &
     &                    +(cf(i)-cff(i))/(timestep*iterations)
            cfltransfer(i)=cfltransfer(i)                               &
     &                    -cfl(i)/(timestep*iterations)
            If (l_seq) Then
              cff(i)=cf(i)
              cfl(i)=0.0
            End if  ! l_seq
          End If  ! l_update_cf

          !-----------------------------------------------
          ! Update water contents
          !-----------------------------------------------
          homtransfer(i) = homtransfer(i) + qcl(i)/(timestep*iterations)
          If (l_seq) Then
            qcf(i) = qcf(i) + qcl(i)
            T(i)   = T(i)   +lfrcp*qcl(i)
            qcl(i) = 0.0
          End if  ! l_seq

        End if ! T lt 0+thomo etc.

        !-----------------------------------------------
        ! Heterogeneous nucleation for T < tnuc
        !-----------------------------------------------
        If (T(i) <  (zerodegc+tnuc) .and. area_liq(i) >  0.0            &
     &      .and. cfliq(i)  >   0.0 ) then

          !-----------------------------------------------
          ! Calculate number and mixing ratio of nucleated crystals
          !-----------------------------------------------
          dqi(i)=min(0.01*exp(-0.6*(T(i)-zerodegc)),1.0E5)
          dqi(i)=m0 * dqi(i) * rhor(i)

          !-----------------------------------------------
          ! How much moisture is available for ice formation
          !-----------------------------------------------
          ! Firstly, calculate the threshold relative humidity
          rhnuc=(188.92+2.81*(T(i)-zerodegc)                            &
     &                +0.013336*(T(i)-zerodegc)**2)*0.01
          rhnuc=min(rhnuc,1.0)-0.1
          rhnuc=max(qsl(i)*rhnuc,qs(i))

          ! Next calculate the available water
          dqil(i)=(qcl(i)/cfliq(i)+qsl(i)-rhnuc)
          dqi(i) = max(area_liq(i)*                                     &
     &                 min(dqi(i),dqil(i)*lheat_correc_ice(i)),0.0)
          qcf(i) = qcf(i)+dqi(i)

          !-----------------------------------------------
          ! Store nucleation rate
          !-----------------------------------------------
          hettransfer(i) = hettransfer(i)+dqi(i)/(timestep*iterations)

          !-----------------------------------------------
          ! Calculate mass transfers
          !-----------------------------------------------
          If (l_seq) Then
            ! Firstly take mass from liquid water
            dqil(i) = min(dqi(i),qcl(i)*area_liq(i)/cfliq(i))
            qcl(i)  = qcl(i)-dqil(i)
            T(i)    = T(i)+lfrcp*dqil(i)
            ! If more ice is needed then the mass comes from vapour
            dqi(i)  = dqi(i)-dqil(i)
            T(i)    = T(i)+lsrcp*dqi(i)
            q(i)    = q(i)-dqi(i)
          End if  ! l_seq

          !-----------------------------------------------
          ! Udate cloud fractions.
          !-----------------------------------------------
          If (l_update_cf) then
            cfftransfer(i) = cfftransfer(i)                             &
     &                     + (cf(i) - cff(i))/(timestep*iterations)
            If (l_seq) then
              cff(i)      = cf(i)
            End if  ! l_seq
          End if  ! l_update cf

        End if  ! On temperature threshold

      End do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_NUCLEATION
#endif
