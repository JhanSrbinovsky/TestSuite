
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
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
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

!
! Description:
!
!  Contains various cloud droplet parameters, defined for
!  land and sea areas.
!
!  NTOT_* is the total number concentration (m-3) of cloud droplets;
!  KPARAM_* is the ratio of the cubes of the volume-mean radius
!                                           and the effective radius;
!  DCONRE_* is the effective radius (m) for deep convective clouds;
!  DEEP_CONVECTION_LIMIT_* is the threshold depth (m) bewteen shallow
!                                          and deep convective cloud.
!
! Current Code Owner: Andy Jones
!
! History:
!
! Version   Date     Comment
! -------   ----     -------
!    1     040894   Original code.    Andy Jones
!  5.2     111000   Updated in line with Bower et al. 1994 (J. Atmos.
!                   Sci., 51, 2722-2732) and subsequent pers. comms.
!                   Droplet concentrations now as used in HadAM4.
!                                     Andy Jones
!  5.4     02/09/02 Moved THOMO here from C_LSPMIC.      Damian Wilson
!  6.2     17/11/05 Remove variables that are now in UMUI. D. Wilson
!
!     REAL,PARAMETER:: NTOT_LAND is set in UMUI
!     REAL,PARAMETER:: NTOT_SEA is set in UMUI
      REAL,PARAMETER:: KPARAM_LAND = 0.67
      REAL,PARAMETER:: KPARAM_SEA = 0.80
      REAL,PARAMETER:: DCONRE_LAND = 9.5E-06
      REAL,PARAMETER:: DCONRE_SEA = 16.0E-06
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_LAND = 500.0
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_SEA = 1500.0
!
! Maximum Temp for homogenous nucleation (deg C)
      REAL,PARAMETER:: THOMO = -40.0
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
