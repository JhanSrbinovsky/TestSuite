#if defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Autoconversion of liquid to rain
! Subroutine Interface:
      SUBROUTINE LSP_AUTOC(                                             &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, qcl, qrain, T, p                                                &
                                          ! Water contents, temp and p
     &, cfliq, rhcpt                                                    &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
!    &, cf, cfl                           ! Current cloud fractions for
!                                         ! updating
     &, area_liq, area_mix, area_ice                                    &
                                          ! Cloud fraction overlaps
     &, rainfrac, rain_liq, rain_mix                                    &
                                          ! Rain fractions for updating
     &, rain_ice, rain_clear                                            &
     &, rho, rhor, corr2, ec_auto                                       &
                                          ! Parametrization information
     &, lcrcp                                                           &
                                          ! Microphysical information
     &, l_update_cf, l_seq, l_auto_debias                               &
                                          ! Code options
     &, l_murk, l_use_sulphate_autoconv                                 &
     &, l_seasalt_ccn, l_biomass_ccn, l_ocff_ccn, l_biogenic_ccn        &
     &, l_autoconv_murk, l_autoc_3b, l_autolim_3b, l_mixing_ratio       &
     &, iterations, ptransfer, rftransfer                               &
                                          ! Mass and rainfrac transfer
!    &, cftransfer,cfltransfer            ! Cloud transfer diagnostics
     &, aerosol, SO4_ait, SO4_acc, SO4_dis                              &
                                          ! Murk and sulphate aerosol
     &, sea_salt_film, sea_salt_jet                                     &
                                          ! Sea salt aerosol
     &, bmass_agd, bmass_cld                                            &
                                          ! Biomass aerosol
     &, ocff_agd, ocff_cld                                              &
                                          ! Fossil-fuel org carb aerosol
     &, biogenic                                                        &
                                          ! Biogenic aerosol
     &, snow_depth, land_fract                                          &
                                          ! Land surface information
     &, Ntot_land, Ntot_sea                                             &
                                          ! Fixed value droplet concs.
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud and rain due to autoconversion of liquid to rain
!
! Method:
!   Use a rate based on a power law of liquid content and droplet
!   number with a minimum liquid content for autoconversion.
!
! Current Owner of Code: Jonathan Wilkinson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!
!   There are many routes through this code depending upon the
!   selected options. Vector options have been removed since
!   they are very difficult to trace through.
!   1) Calculate number of droplets based upon either a
!      prescribed number or the amount of aerosol (murk or
!      sulphate/sea-salt/biomass/fossil-fuel organic carbon) in 
!      the model.
!   2) Calculate the autoconversion rate based upon a
!      Tripoli and Cotton power law formulation.
!   3) Calculate the autoconversion limit based either on
!      the concentration of droplets above 20 um diameter or
!      on a Tripoli and Cotton type argument
!
! Subroutine Arguments
!
! Autoconversion parameters are contained in c_lspmic.h
#include "c_lspmic.h"
#include "c_pi.h"
#include "c_epslon.h"
#include "c_r_cp.h"
#include "c_lheat.h"
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
     &, p(points)                                                       &
                          ! Air pressure / N m-2
     &, T(points)                                                       &
                          ! Temperature / K
     &, aerosol(points)                                                 &
                          ! Mass of aerosol / ug kg-1
     &, cfliq(points)                                                   &
                          ! Fraction of gridbox with liquid cloud
!    &, cf(points)        ! Bulk cloud fraction for updating
!    &, cfl(points)       ! Liquid cloud fraction for updating
     &, area_liq(points)                                                &
                          ! Fraction of gridbox with liquid but no ice
     &, area_mix(points)                                                &
                          ! Fraction of gridbox with liquid and ice
     &, area_ice(points)                                                &
                          ! Fraction of gridbox with ice but no liquid
     &, rhcpt(points)                                                   &
                          ! Rhcrit on each point (no units)
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, rhor(points)                                                    &
                          ! 1 / air density / m3 kg-1
     &, corr2(points)                                                   &
                          ! Temperature correction factor (no units)
     &, ec_auto                                                         & 
                          ! Collision - coalescence eff'cy (no units)
     &, lcrcp                                                           &
                          ! Latent heat of condensation/cP / K
     &, SO4_ait(points)                                                 &
                          ! Aitken mode sulphate concentration
     &, SO4_acc(points)                                                 &
                          ! Accumulation mode sulphate conc
     &, SO4_dis(points)                                                 &
                          !
     &, sea_salt_film(points)                                           &
                              ! Film mode sea salt conc
     &, sea_salt_jet(points)                                            &
                              ! Jet mode sea salt conc
     &, bmass_agd(points)                                               &
                          ! Aged biomass smoke mass mixing ratio
     &, bmass_cld(points)                                               &
                          ! In-cloud smoke mass mixing ratio
     &, ocff_agd(points)                                                &
                          ! Aged fossil-fuel org carb mass mixing ratio
     &, ocff_cld(points)                                                &
                          ! In-cloud fossil-fuel organic carbon mmr
     &, biogenic(points)                                                &
                          ! biogenic mass mixing ratio
     &, snow_depth(points)                                              &
                          ! Snow depth / m
     &, land_fract(points)                                              &
                          ! Fraction of gridbox with land
     &, Ntot_land                                                       &
                          ! Droplet number conc. over land / m-3
     &, Ntot_sea
                          ! Droplet number conc. over sea / m-3
!
      Real, Intent(InOut) ::                                            &
     &  qcl(points)                                                     &
                          ! Liquid water content / kg kg-1
     &, qrain(points)                                                   &
                          ! Rain water content / kg kg-1
     &, ptransfer(points)                                               &
                          ! Autoconversion rate / kg kg-1 s-1
!    &, cftransfer(points) ! Rate of change of bulk cloud frac / s-1
!    &, cfltransfer(points)! Rate of change of liquid cloud frac / s-1
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
     &  l_murk                                                          &
                          ! Murk aerosol is active in the model
     &, l_autoconv_murk                                                 &
                          ! Use murk aerosol for droplet number
     &, l_use_sulphate_autoconv                                         &
                                ! Use sulphate aerosol to calc number
     &, l_seasalt_ccn                                                   &
                          ! Sea salt contributes to droplet number
     &, l_biomass_ccn                                                   &
                          ! Biomass contributes to droplet number
     &, l_ocff_ccn                                                      &
                          ! Fossil-fuel organic carbon contributes to
                          ! droplet number
     &, l_biogenic_ccn                                                  &
                          ! Biogenic contributes to droplet number
     &, l_auto_debias                                                   &
                          ! Debias autoconversion for subgrid
                          !   distribution of moisture
     &, l_seq                                                           &
                          ! Sequential updating
     &, l_update_cf                                                     &
                          ! Update cloud fractions
     &, l_autoc_3b                                                      &
                          ! Use 3B autoconversion method
     &, l_autolim_3b                                                    &
                          ! Use fixed 3B values for the autoconv limit 
     &, l_mixing_ratio
                          ! Use mixing ratio
                          ! (otherwise specific humidities)
!
! Local Variables
!
      Integer                                                           &
     &  kk                                                              &
                          ! Index for condensed points
     &, i                                                               &
                          ! Loop counter
     &, sea_salt_ptr                                                    &
                          ! Pointer for sea salt array
     &, biomass_ptr                                                     &
                          ! Pointer for biomass array
     &, ocff_ptr                                                        &
                          ! Pointer for fossil-fuel org carbon array
     &, biogenic_ptr
                          ! Pointer for biogenic aerosol arrays
!
      Real                                                              &
     &  n_drop(points)                                                  &
                          ! Number of droplets / m-3
     &, dpr(points)                                                     &
                          ! Amount of mass autoconverted / kg kg-1
     &, qsl_tl(points)                                                  &
                          ! Saturated humidity wrt liquid at
                          !  temperature TL / kg kg-1
     &, T_L(points)                                                     &
                          ! Liquid temperature / K
     &, r_mean(points)                                                  &
                          ! Mean radius of droplets / m
     &, r_mean0                                                         &
                          ! Factor in calculation of r_mean
     &, a_factor(points)                                                &
                          ! 3 r_mean / m
     &, b_factor(points)                                                &
                          ! 0.5 n_drop / m-3
     &, n_gt_20(points)                                                 &
                          ! Concentration of droplets with
                          !  diameter gt 20 um / m-3
     &, autorate(points)                                                &
                          ! Rate constants for autoc.
     &, autolim(points)                                                 &
                          ! Autoconversion threshold / kg kg-1
     &, qc(points)                                                      &
                          ! Maximum amout of liquid to remove / kg kg-1
!    &, delta_cf(points)  ! Cloud fraction transferred (no units)
!    &, delata_cfl(points)! Liquid cloud fraction transferred (no units)
     &, rainfracnew(points)                                             &
                          ! Updated rain fraction
     &, rho_1(points)
                          ! rho^(power_rho_auto - power_qcl_auto + 1)

      Real                                                              &
                          ! For debiasing code
     &  alpha_l                                                         &
                          ! dqsat/dT at T_L / kg kg-1 K-1
     &, a_L                                                             &
                          ! 1 / (1 + L/cp alpha)
     &, sigma_s                                                         &
                          ! Width of moisture PDF
     &, g_L                                                             &
                          ! Factor in debiasing calculation
     &, gacb                                                            &
                          ! Factor in debiasing calculation
     &, ac_factor
                          ! Multiplying factor for debiasing
!
!  Function calls
      Real, External :: Number_Droplet

      !-----------------------------------------------
      ! Part 1. Calculate droplet number concentration
      !-----------------------------------------------
      kk=0 ! Reset index counter for condensed points

      If (l_autoconv_murk .and. l_murk) then

        !-----------------------------------------------
        ! Calculate droplet concentration using murk aerosol
        !-----------------------------------------------
        Do i=1,points

          If (qcl(i) >  0.0.and.cfliq(I) >  0.0) then
            ! Autoconversion is active
            kk=kk+1
            !-----------------------------------------------
            ! Convert aerosol mass to aerosol number.
            ! This follows subroutine VISBTY
            !-----------------------------------------------
            n_drop(kk) = max(aerosol(i)/m0_murk*1.0E-9, 0.0001)
            ! 1.0E-9 converts from ug/kg to kg/kg

            !-----------------------------------------------
            ! Calculation of the droplet number
            !-----------------------------------------------
            n_drop(kk) = n0_murk * n_drop(kk)**power_murk
          End If  ! qcl >  0
        End Do  ! points

      Else if (l_use_sulphate_autoconv) then

        !-----------------------------------------------
        ! Full droplet number calculation
        ! Set up pointers for droplet number calculation
        !-----------------------------------------------
        Do i=1, points

          If (qcl(i)  >   0.0 .and. cfliq(i)  >   0.0) Then
            ! Autoconversion is active

            ! Do we include sea salt aerosol?
            If (l_seasalt_ccn) Then
              sea_salt_ptr = i
            Else
              sea_salt_ptr = 1
            End If

            ! Do we include biomass aerosol?
            If (l_biomass_ccn) Then
              biomass_ptr = i
            Else
              biomass_ptr = 1
            End if
            
            ! Do we include fossil-fuel organic carbon aerosol?
            If (l_ocff_ccn) Then
              ocff_ptr = i
            Else
              ocff_ptr = 1
            End if
            

            ! Do we include biogenic aerosol?
            IF (L_biogenic_CCN) THEN
              biogenic_ptr=i
            ELSE
              biogenic_ptr=1
            ENDIF

            ! Increment condensed points counter
            kk = kk + 1

            !-----------------------------------------------
            ! Calculate droplet number
            !-----------------------------------------------
! DEPENDS ON: number_droplet
            n_drop(kk) = number_droplet(l_use_sulphate_autoconv,        &
     &                     .false., SO4_ait(i), SO4_acc(i), SO4_dis(i), &
     &                     l_seasalt_ccn, sea_salt_film(sea_salt_ptr),  &
     &                     sea_salt_jet(sea_salt_ptr),                  &
     &                     L_biogenic_CCN, biogenic(biogenic_ptr),      &
     &                     l_biomass_ccn,                               &
     &                     bmass_agd(biomass_ptr),                      &
     &                     bmass_cld(biomass_ptr),                      &
     &                     l_ocff_ccn,                                  &
     &                     ocff_agd(ocff_ptr), ocff_cld(ocff_ptr),      &
     &                     rho(i),                                      &
     &                     snow_depth(i), land_fract(i),                &
     &                     Ntot_land, Ntot_sea)
          End If  ! qcl gt 0

        End Do  ! points

      Else

        !-----------------------------------------------
        ! Simple land / sea n_drop calculation reproduced from n_drop
        !-----------------------------------------------
        Do i=1,points
          If (qcl(i) >  0.0.and.cfliq(i) >  0.0) then
            kk = kk + 1

            If (land_fract(i) >=  0.5) then
              n_drop(kk)=Ntot_land
            Else
              n_drop(kk)=Ntot_sea
            End if

          End if
        End do

      End If  ! l_autoconv_murk

      !-----------------------------------------------
      ! Part 2: Calculating the autoconversion rate and limit
      !-----------------------------------------------
      If (l_autoc_3b) Then

        !-----------------------------------------------
        ! Use the original 3B method for autoconversion
        !-----------------------------------------------
        kk=0

        Do i=1, points

          ! The numbers 5907.24 and 4188.79 represent combinations
          ! of physical constants. Do NOT change them.
          ! Please refer to UMDP26
          If (qcl(i) >  0.0.and.cfliq(i) >  0.0) Then

            kk=kk+1
            If (n_drop(kk) > 0.0) then
              Autorate(i) = 5907.24*ec_auto*inhomog_rate                &
     &                           *n_drop(kk)**(-1.0/3.0)
              Autolim(i)  = 4188.79*r_thresh**3*n_drop(kk)*inhomog_lim
            Else
              ! No droplets so no autoconversion
              Autorate(i) = 0.0
              Autolim(i)  = 0.0
            End if  ! n_drop >0

          End if  ! qcl > 0

        End do  ! points

      Else  ! l_auto_3b
        !-----------------------------------------------
        ! Use the later 3C/3D method for autoconversion
        !-----------------------------------------------

        kk=0
        ! Calculate multiplying factor for droplet size
        r_mean0=(27.0/(80.0*pi*1000.0))**(-1.0/3.0)

        Do i=1,points
          If (qcl(i) >  0.0.and.cfliq(i) >  0.0) Then

            ! Autoconversion is active
            kk=kk+1

            ! Calculate inverse of mean droplet size
            r_mean(kk) = r_mean0 * (rho(i) * qcl(i)                     &
     &                   / (cfliq(i) * n_drop(kk))) ** (-1.0/3.0)

            ! Calculate numerical factors
            b_factor(kk) = 3.0 * r_mean(kk)
            a_factor(kk) = n_drop(kk) * 0.5

            !-----------------------------------------------
            ! Calculate droplet number concentration greater
            ! than threshold radius
            !-----------------------------------------------
            n_gt_20(i) = (b_factor(kk)**2 * a_factor(kk)) * r_auto**2   &
     &                                  * exp(-b_factor(kk) * r_auto)   &
     &            + (2.0 * b_factor(kk) * a_factor(kk)) * r_auto        &
     &                                  * exp(-b_factor(kk) * r_auto)   &
     &            + (2.0 * a_factor(kk))                                &
     &                                  * exp(-b_factor(kk) * r_auto)

            !-----------------------------------------------
            ! Test to see if there is a sufficient concentration of
            ! droplets with diameters > threshold for autoconversion
            ! to proceed.
            !-----------------------------------------------
            If (n_gt_20(i)  >=  n_auto) Then
              ! Calculate autoconversion rate
              autorate(i) = consts_auto * ec_auto                       &
     &                    * n_drop(kk) ** power_droplet_auto
            Else
              ! No autoconversion
              autorate(i)=0.0
            End If ! n_gt_20 ge n_auto

            !-----------------------------------------------
            ! Calculate the autoconversion limit
            !-----------------------------------------------
            ! Calculate value of local qcl at which the droplet
            ! concentration with radii greater than 20um will
            ! fall below a threshold (1000 m-3). This is a
            ! hardwired numerical approximation
            autolim(i)=(6.20E-31*n_drop(kk)**3)-(5.53E-22*n_drop(kk)**2)&
     &               +(4.54E-13*n_drop(kk))+(3.71E-6)-(7.59/n_drop(kk))

          End If  ! qcl(i) >  0.0
        End Do  ! points

      End if ! l_autoc_3b

      If (l_autolim_3b) then

        !-----------------------------------------------
        ! Overwrite the autoconversion limit with 3B values
        !-----------------------------------------------      
        Do i=1,points
          If (land_fract(i) .ge. 0.5) Then  
            autolim(i)=8.621e-4
          Else                           
            autolim(i)=2.155e-4
          End if  ! land_fract ge 0.5  
        End Do

      End if ! l_autolim3b

      !-----------------------------------------------
      ! Part 3: Optionally debias the autoconversion rate
      !-----------------------------------------------
      If (l_auto_debias) Then

        !-----------------------------------------------
        ! Calculate qsat(TL) in order to allow subgrid
        ! debiasing of the autoconversion rate.
        !-----------------------------------------------
        Do i=1,points
          T_L(i) = T(i) - (lcrcp * qcl(i) )
        End Do
! DEPENDS ON: qsat_wat_mix
        Call qsat_wat_mix(qsl_tl,T_L,p,points,l_mixing_ratio)

        kk=0
        Do i=1, points
          If (qcl(i) >  0.0.and.cfliq(i) >  0.0) Then

            ! Autoconversion is active
            kk=kk+1          

            !-----------------------------------------------
            ! De-bias the autoconversion rate based on a
            ! triangular Smith PDF, following Wood et al
            ! (Atmos. Res., 65, 109-128, 2002).
            !-----------------------------------------------

            If (qcl(i)  <   1.0E-15) Then
              ! Water contents are small so set factor to 1
              ac_factor=1.0

            Else  ! qcl lt 1e-15
              alpha_l = epsilon * lc * qsl_tl(i) / ( R * T_L(i)**2 )
              a_L = 1.0 / (1.0+(lcrcp*alpha_l))
              sigma_s = (1.0 - rhcpt(i)) * a_L * qsl_tl(i) / sqrt(6.0)
              g_l = 1.15 * (power_qcl_auto-1.0) * sigma_s
              gacb = exp(-1.0 * qcl(i) / g_l)

              ! Calculate autoconversion rate multiplication factor
              ac_factor = max(                                          &
     &                  (cfliq(i)**(power_qcl_auto-1.0))*1.0/(1.0-gacb) &
     &                  ,1.0)

            End If  ! qcl lt 1e-15

            !-----------------------------------------------
            ! Apply the debiasing factor
            !-----------------------------------------------
            autorate(i) = ac_factor * autorate(i)

          End if ! qcl > 0

        End do  ! points

      End If  ! l_auto_debias

      !-----------------------------------------------
      ! Part 4. Calculate the autoconversion
      !-----------------------------------------------
      Do i=1, points

        !-----------------------------------------------
        ! Set the dependence of autoconversion on air density
        !-----------------------------------------------
        ! power_rho_auto and power_qcl_auto are set in c_lspmic.
        rho_1(i) = rho(i) ** (power_rho_auto-power_qcl_auto+1.0)

        If (qcl(i) >  0.0.and.cfliq(i) >  0.0) Then

          ! Calculate maximum amount of liquid that can be removed
          ! from the grid box
          qc(i) = min( autolim(i) * cfliq(i) * rhor(i) , qcl(i) )

          !-----------------------------------------------
          ! Calculate the autoconversion amount
          !-----------------------------------------------
          dpr(i) = min(autorate(i)                                      &
     &              *(rho(i)*qcl(i)/cfliq(i))**(power_qcl_auto-1.0)     &
     &              *timestep*qcl(i)*rho_1(i)/corr2(i),qcl(i)-qc(i))

          !-----------------------------------------------
          ! Update liquid water content and rain
          !-----------------------------------------------
          If (l_seq) Then
            qcl(i)   = qcl(i)   - dpr(i)
            qrain(i) = qrain(i) + dpr(i)
          End if

          !-----------------------------------------------
          ! Store process rate
          !-----------------------------------------------
          ptransfer(i) = ptransfer(i) + dpr(i) / (timestep*iterations)

          If (dpr(i) >  0.0) Then
            !-----------------------------------------------
            ! Calculate change in rain fraction
            !-----------------------------------------------
            rainfracnew(i) = max(rainfrac(i),cfliq(i))
            rftransfer(i) = rftransfer(i)                               &
     &        + (rainfracnew(i) - rainfrac(i)) / (timestep*iterations)

            !-----------------------------------------------
            ! Update rain fractions
            !-----------------------------------------------
            If (l_seq) Then
              rainfrac(i)   = rainfracnew(i)
              rain_liq(i)   = min(area_liq(i),rainfrac(i))
              rain_mix(i)   = min(area_mix(i),rainfrac(i)-rain_liq(i))
              rain_ice(i)   =                                           &
     &            min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
              rain_clear(i) =                                           &
     &            rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i)
            End If  ! l_seq

          End If  ! dpr gt 0

          !-----------------------------------------------
          ! Update cloud fractions
          !-----------------------------------------------
!          If (l_update_cf) Then
!            These are commented out since there is currently no
!            cloud fraction update associated with the autoconversion.
!            If (l_seq) Then
!              cf(i)  = cf(i) + delta_cf(i)
!              cfl(i) = cfl(i)+ delta_cfl(i)
!            End If  ! l_seq
!            cftransfer(i)  = cftransfer(i)  + delta_cf(i)
!     &                                      / (timestep*iterations)
!            cfltransfer(i) = cfltransfer(i) + delta_cfl(i)
!     &                                      / (timestep*iterations)
!          End If  ! l_update_cf

        End If  ! qcl gt 0
      End Do  ! points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_AUTOC
#endif
