#if defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Subgrid-scale set ups and checks
! Subroutine Interface:
      SUBROUTINE LSP_SUBGRID(                                           &
     &  points                                                          &
                                          ! Number of points
     &, q, qcf_cry, qcf_agg, qcftot, T                                  &
                                          ! Water contents and temp
     &, qsl, qs                                                         &
                                          ! Saturated water contents
     &, q_ice, q_clear, q_ice_1, q_ice_2                                &
                                          ! Local vapour contents
     &, area_liq,area_mix,area_ice,area_clear                           &
                                              ! Cloud frac partitions
     &, area_ice_1, area_ice_2                                          &
                                          ! Subdivision of area_ice
     &, rain_liq,rain_mix,rain_ice,rain_clear                           &
                                              ! Rain overlap partitions
     &, cftot, cfliq, cfice, cficei                                     &
                                          ! Cloud fractions for
     &, frac_ice_above                                                  &
                                          ! partition calculations
     &, cf, cff, rainfrac                                               &
                                          ! Cloud and rain fractions
                                          ! for updating
     &, lsrcp                                                           &
                                          ! Latent heat of sublim./cp
     &, rhcpt                                                           &
                                          ! RH crit values
     &, l_update_cf                                                     &
                                          ! Code options
     &  )
!
      Implicit None
!
! Purpose:
!   Perform the subgrid-scale setting up calculations
!
! Method:
!   Parametrizes the width of the vapour distribution in the part
!   of the gridbox which does not have liquid water present.
!   Calculates the overlaps within each gridbox between  the cloud
!   fraction prognostics and rainfraction diagnostic.
!
! Current Owner of Code: Jonathan Wilkinson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!
! The subgrid calculations are a necessary step to calculating the
! subsequent deposition and sublimation transfers and for setting up
! the partition information that is used by the subsequent transfers.
!
! Subroutine Arguments
!
#include "c_lspmic.h"
#include "c_0_dg_c.h"
!
      Integer, Intent(In) ::                                            &
     &  points            ! Number of points to calculate
!
      Real, Intent(In) ::                                               &
     &  qs(points)                                                      &
                          ! Saturated mixing ratio wrt ice
     &, qsl(points)                                                     &
                          ! Saturated mixing ratio wrt liquid
     &, cfliq(points)                                                   &
                          ! Fraction of gridbox with liquid cloud
     &, rainfrac(points)                                                &
                          ! Fraction of gridbox containing rain
     &, frac_ice_above(points)                                          &
                               ! Ice cloud in level above this one
     &, lsrcp                                                           &
                          ! Latent heat of sublimation
                          ! / heat capacity of air / K
     &, rhcpt(points)     ! RH crit values
!
      Real, Intent(InOut) ::                                            &
     &  q(points)                                                       &
                          ! Vapour content / kg kg-1
     &, qcf_cry(points)                                                 &
                          ! Ice crystal content / kg kg-1
     &, qcf_agg(points)                                                 &
                          ! Ice aggregate content / kg kg-1
     &, qcftot(points)                                                  &
                          ! Total ice content before advection / kg kg-1
     &, T(points)                                                       &
                          ! Temperature / K
     &, cf(points)                                                      &
                          ! Current cloud fraction
     &, cff(points)       ! Current ice cloud fraction
!
      Real, Intent(Out) ::                                              &
     &  q_clear(points)                                                 &
                          ! Local vapour in clear-sky region / kg kg-1
     &, q_ice(points)                                                   &
                          ! Local vapour in ice-only region  / kg kg-1
     &, q_ice_1(points)                                                 &
                          ! Local vapour in ice-only regions that are:
     &, q_ice_2(points)                                                 &
                          !   1, depositing; and 2, subliming.
     &, cftot(points)                                                   &
                          ! Modified cloud fraction for partition calc.
     &, cfice(points)                                                   &
                          ! Modified ice cloud frac. for partition calc.
     &, cficei(points)                                                  &
                          ! 1/cfice
     &, area_liq(points)                                                &
                          ! Frac of gridbox with liquid cloud but no ice
     &, area_mix(points)                                                &
                          ! Frac of gridbox with liquid and ice cloud
     &, area_ice(points)                                                &
                          ! Frac of gridbox with ice cloud but no liquid
     &, area_clear(points)                                              &
                          ! Frac of gridbox with no cloud
     &, area_ice_1(points)                                              &
                          ! Frac of gridbox where ice-only cloud is:
     &, area_ice_2(points)                                              &
                          !  1, depositing; and 2, subliming.
     &, rain_liq(points)                                                &
                          ! Frac of gbox with rain and liquid but no ice
     &, rain_mix(points)                                                &
                          ! Frac of gbox with rain and liquid and ice
     &, rain_ice(points)                                                &
                          ! Frac of gbox with rain and ice but no liquid
     &, rain_clear(points)! Frac of gbox with rain but no condensate

      Logical, Intent(In) ::                                            &
     &  l_update_cf       ! Update cloud fractions
!
! Local Variables
!
      Integer                                                           &
     &  i                 ! Loop counter

      Real                                                              &
     &  tempw(points)                                                   &
                          ! Vapour content in ice and clear partitions
     &, temp7(points)                                                   &
                          ! Temporary in width of PDF calculation
     &, width(points)     ! Full width of vapour distribution in ice and
                          ! clear sky.

      Do i = 1, points

        !-----------------------------------------------
        ! Check that ice cloud fraction is sensible.
        !-----------------------------------------------
        ! Difference between the way PC2 and non-PC2 code operates
        ! is kept here in order to be tracable across model versions.
        ! However, perhaps the code ought to be the same.
        If (l_update_cf) then
          ! 0.001 is to avoid divide by zero problems
          cfice(i) = max(cff(i),0.001)
          cficei(i) = 1.0/cfice(i)
          cftot(i)=cf(i)
          cftot(i)=min( max(cftot(i),cfice(i)) ,(cfice(i)+cfliq(i)) )
        Else
          ! 0.01 is to avoid divide by zero problems
          cfice(i) = max(max(cff(i),0.01),frac_ice_above(i))
          cficei(i) = 1.0/cfice(i)
          cftot(i)=cf(i)
        End If ! l_ update_cf

        ! -----------------------------------------------
        ! Calculate overlaps of liquid, ice and rain fractions
        ! -----------------------------------------------
        area_liq(i) = max(cftot(i)-cfice(i),0.0)
        area_mix(i) = max(cfice(i)+cfliq(i)-cftot(i),0.0)
        area_ice(i) = max(cftot(i)-cfliq(i),0.0)
        area_clear(i) = max(1.0-cftot(i),0.0)
        rain_liq(i) = max(min(area_liq(i),rainfrac(i)),0.0)
        rain_mix(i) = max(min(area_mix(i),rainfrac(i)-rain_liq(i)),0.0)
        rain_ice(i) =                                                   &
     &    max(min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i)),0.0)
        rain_clear(i) =                                                 &
     &    max(rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i),0.0)

        If (cfliq(i)  <   1.0) Then

          ! -----------------------------------------------
          ! Calculate width of vapour dist. in ice and clear region
          ! -----------------------------------------------
          ! tempw is the mean vapour content in the ice only and clear
          ! sky partitions
          tempw(i) = (q(i) - cfliq(i)*qsl(i)) / (1.0 - cfliq(i))
          temp7(i) = ice_width * qsl(i)
          ! 0.001 is to avoid divide by zero problems
          width(i) = 2.0 *(1.0-rhcpt(i))*qsl(i)                         &
     &                   *max(  (1.0-0.5*qcftot(i)/temp7(i)),0.001  )
          ! The full width cannot be greater than 2q because otherwise
          ! part of the gridbox would have negative q. Also ensure that
          ! the full width is not zero (possible if rhcpt is 1).
          width(i) = min(width(i) , max(2.0*q(i),0.001*qs(i) )    )

          ! -----------------------------------------------
          ! Calculate vapour contents in ice only and clear regions
          ! -----------------------------------------------
          If (area_ice(i)  >   0.0) Then
             q_clear(i) = tempw(i) - 0.5*width(i) * area_ice(i)
             q_ice(i) = (q(i)-cfliq(i)*qsl(i)-area_clear(i)*q_clear(i)) &
     &                / area_ice(i)
          Else
             q_clear(i) = tempw(i)
             q_ice(i) = 0.0               ! q_ice is a dummy value here
          End If  ! area_ice gt 0

        Else ! cf_liq lt 1

          ! -----------------------------------------------
          ! Specify dummy values for q_clear and q_ice
          ! -----------------------------------------------
          width(i)   = 1.0
          q_clear(i) = 0.0
          q_ice(i)   = 0.0

        End If ! cf_liq lt 1

        ! -------------------------------------------------
        ! Remove any small amount of ice to be tidy.
        ! -------------------------------------------------
        ! If QCF is less than QCFMIN and isn't growing by deposition
        ! (assumed to be given by RHCPT) then evaporate it.
        If ((qcf_cry(i)+qcf_agg(i)) <  qcfmin) Then
          If (T(i) >  zerodegc .or.                                     &
     &       (q_ice(i)  <=  qs(i) .and. area_mix(i)  <=  0.0)           &
     &       .or. (qcf_cry(i)+qcf_agg(i)) <  0.0)  Then
            q(i) = q(i) +qcf_cry(i)+qcf_agg(i)
            T(i) = T(i) - lsrcp * (qcf_cry(i)+qcf_agg(i))
            qcf_cry(i)=0.0
            qcf_agg(i)=0.0
          End If ! T gt 0 etc.
        End If ! qcf_cry+qcf_agg lt qcfmin

        ! -------------------------------------------------
        ! First estimate of partition sizes for ice sublimation
        ! and deposition and vapour contents within these partitions
        ! -------------------------------------------------
        If (q_ice(i)  >   qs(i)) Then
          ! First estimate is to use a deposition process
          area_ice_1(i) = area_ice(i)
          area_ice_2(i) = 0.0
          q_ice_1(i)    = q_ice(i)
          q_ice_2(i)    = qs(i)       ! Dummy value

        Else ! q_ice gt qs
          ! First estimate is to use a sublimation process
          area_ice_1(i) = 0.0
          area_ice_2(i) = area_ice(I)
          q_ice_1(i)    = qs(i)       ! Dummy value
          q_ice_2(i)    = q_ice(i)

        End If ! q_ice gt qs

        ! -------------------------------------------------
        ! Detailed estimate of partition sizes for ice sublimation
        ! and deposition and vapour contents within these partitions
        ! -------------------------------------------------
        If (area_ice(i)  >   0.0) Then
        ! Temp7 is the estimate of the proportion of the gridbox
        ! which contains ice and has local q > than qs (wrt ice)
        temp7(i) = 0.5*area_ice(i) + (q_ice(i)-qs(i)) / width(i)

          If (temp7(i) >  0.0 .and. temp7(i) <  area_ice(i)) Then
            ! Calculate sizes of regions and q in each region
            ! These overwrite previous estimates
            area_ice_1(i) = temp7(i)
            area_ice_2(I) = area_ice(i) - area_ice_1(i)
            q_ice_1(i) = qs(i) + 0.5 * area_ice_1(i) * width(i)
            q_ice_2(i) = qs(i) - 0.5 * area_ice_2(i) * width(i)
          End If ! temp7 gt 0 etc.

        End If ! area_ice gt 0

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_SUBGRID
#endif
