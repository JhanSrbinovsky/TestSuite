
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
! --------------------------COMDECK C_LSPMIC----------------------------
! SPECIFIES MICROPHYSICAL PARAMETERS FOR AUTOCONVERSION, HALLETT MOSSOP
! PROCESS, ICE NUCLEATION. ALSO SPECIFIES NUMBER OF ITERATIONS OF
! THE MICROPHYSICS AND ICE CLOUD FRACTION METHOD
! ----------------------------------------------------------------------
!
! History:
!
! Version    Date     Comment
! -------    ----     -------
!   5.4    16/08/02   Correct comment line, add PC2 parameters and
!                     move THOMO to c_micro     Damian Wilson
!   6.0    11/08/03   Correct value of wind_shear_factor for PC2
!                                                          Damian Wilson
!   6.2    17/11/05   Remove variables that are now in UMUI. D. Wilson
!   6.2    03/02/06   Include droplet settling logical. Damian Wilson
!
! ----------------------------------------------------------------------
!      AUTOCONVERSION TERMS
! ----------------------------------------------------------------------
!
!     LOGICAL, PARAMETER :: L_AUTOCONV_MURK is set in UMUI
! Set to .TRUE. to calculate droplet concentration from MURK aerosol,
! which will override L_USE_SULPHATE_AUTOCONV (second indirect effect
! of sulphate aerosol). If both are .FALSE., droplet concentrations
! from comdeck C_MICRO are used to be consistent with the values
      ! used in the radiation scheme.

      ! This next set of parameters is to allow the 3B scheme to
      ! be replicated at 3C/3D
        ! Inhomogeneity factor for autoconversion rate
        REAL,PARAMETER:: INHOMOG_RATE=1.0

        ! Inhomogeneity factor for autoconversion limit
        REAL,PARAMETER:: INHOMOG_LIM=1.0

        ! Threshold droplet radius for autoconversion
        REAL,PARAMETER:: R_THRESH=7.0E-6
      ! End of 3B repeated code

      !Do not alter R_AUTO and N_AUTO since these values are effectively
      ! hard wired into a numerical approximation in the autoconversion
      ! code. EC_AUTO will be multiplied by CONSTS_AUTO

      ! Threshold radius for autoconversion
      REAL, PARAMETER :: R_AUTO=20.0E-6

      ! Critical droplet number for autoconversion
      REAL, PARAMETER :: N_AUTO=1000.0

      ! Collision coalesence efficiency for autoconversion
!      REAL, PARAMETER :: EC_AUTO is set in UMUI

      ! The autoconversion powers define the variation of the rate with
      ! liquid water content and droplet concentration. The following are
      ! from Tripoli and Cotton

      !  Dependency of autoconversion rate on droplet concentration
      REAL, PARAMETER :: POWER_DROPLET_AUTO=-0.33333

      ! Dependency of autoconversion rate on water content
      REAL, PARAMETER :: POWER_QCL_AUTO=2.33333

      ! Dependency of autoconversion rate on air density
      REAL, PARAMETER :: power_rho_auto=1.33333

      ! CONSTS_AUTO = (4 pi)/( 18 (4 pi/3)^(4/3)) g /  mu (rho_w)^(1/3)
      ! See UM documentation paper 26, equation P26.132

      ! Combination of physical constants
      REAL, PARAMETER :: CONSTS_AUTO=5907.24

      ! Quantites for calculation of drop number by aerosols.
      ! Need only set if L_AUTOCONV_MURK=.TRUE.  See file C_VISBTY

      ! Scaling concentration (m-3) in droplet number concentration
      REAL, PARAMETER :: N0_MURK=500.0E6

      ! Scaling mass (kg/kg) in droplet number calculation from aerosols
      REAL, PARAMETER :: M0_MURK=1.458E-8

      ! Power in droplet number calculation from aerosols
      REAL, PARAMETER :: POWER_MURK=0.5

      ! Ice water content threshold for graupel autoconversion (kg/m^3)
      REAL, PARAMETER :: AUTO_GRAUP_QCF_THRESH = 3.E-4

      ! Temperature threshold for graupel autoconversion (degC)
      REAL, PARAMETER :: AUTO_GRAUP_T_THRESH = -4.0

      ! Temperature threshold for graupel autoconversion
      REAL, PARAMETER :: AUTO_GRAUP_COEFF = 0.5

      !-----------------------------------------------------------------
      ! Iterations of microphysics
      !-----------------------------------------------------------------

      ! Number of iterations in microphysics.
      INTEGER,PARAMETER :: LSITER=1
      ! Advise 1 iteration for every 10 minutes or less of timestep.

      !-----------------------------------------------------------------
      ! Nucleation of ice
      !-----------------------------------------------------------------

      ! Note that the assimilation scheme uses temperature thresholds
      ! in its calculation of qsat.

      ! Nucleation mass
      REAL, PARAMETER :: M0=1.0E-12

      ! Maximum Temp for ice nuclei nucleation (deg C)
      REAL, PARAMETER :: TNUC=-10.0

      ! Maximum temperature for homogenous nucleation is now in c_micro
      ! so that it is available to code outside of section A04.

      !  1.0/Scaling quantity for ice in crystals
      REAL, PARAMETER :: QCF0=1.0E4       ! This is an inverse quantity

      ! Minimum allowed QCF after microphysics
      REAL,PARAMETER:: QCFMIN=1.0E-8

      ! 1/scaling temperature in aggregate fraction calculation
      REAL, PARAMETER :: T_SCALING=0.0384

      !  Minimum temperature limit in calculation  of N0 for ice (deg C)
      REAL, PARAMETER :: T_AGG_MIN=-45.0

      !-----------------------------------------------------------------
      ! Hallett Mossop process
      !-----------------------------------------------------------------

      ! Switch off Hallett Mossop in this version but allow
      ! functionality

      ! Min temp for production of Hallett Mossop splinters (deg C)
      REAL, PARAMETER :: HM_T_MIN=-8.0

      ! Max temp for production of Hallett Mossop splinters (deg C)
      REAL, PARAMETER :: HM_T_MAX=-273.0
      ! REAL, PARAMETER :: HM_T_MAX=-3.0

      !  Residence distance for Hallett Mossop splinters (1/deg C)
      REAL, PARAMETER :: HM_DECAY=1.0/7.0

      ! Reciprocal of scaling liquid water content for HM process
      REAL, PARAMETER :: HM_RQCL=1.0/0.1E-3

      !-----------------------------------------------------------------
      ! PC2 Cloud Scheme Terms
      !-----------------------------------------------------------------

      ! Specifies the ice content (in terms of a fraction of qsat_liq)
      ! that corresponds to a factor of two reduction in the width of
      ! the vapour distribution in the liquid-free part of the gridbox.
      REAL, PARAMETER :: ICE_WIDTH=0.04

      ! Parameter that governs the rate of spread of ice cloud fraction
      ! due to windshear
      REAL, PARAMETER :: WIND_SHEAR_FACTOR = 1.5E-4

!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
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
