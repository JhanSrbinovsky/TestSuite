
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Iterative melting interface
! Subroutine Interface:
      SUBROUTINE LSP_IT_FALL_MELT(                                      &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, T, p, q, q_ice, qsl                                             &
                                          ! Temperature and vapour
     &, qcf_cry, qcf_agg, frac_agg                                      &
                                          ! Ice contents
     &, qrain, qgraup                                                   &
                                          ! Rain and graupel contents
     &, snow_agg,snow_cry,rainrate,grauprate                            &
                                             ! Sedimentation into layer
     &, snowt_agg,snowt_cry,rainratet,graupratet                        &
                                                 ! Sedim. out of layer
     &, vf_agg, vf_cry, vf_rain, vf_graup                               &
                                          ! Fall speeds of hydrometeors
     &, area_liq, area_mix, area_clear                                  &
                                          ! Cloud fraction information
     &, area_ice,cfice,cficei                                           &
                                          ! at start of microphy. ts
     &, frac_ice_above                                                  &
     &, cf, cfl, cff                                                    &
                                          ! Current cloud fractions for
                                          ! updating
     &, rainfrac, rain_liq, rain_mix                                    &
                                          ! Rain fraction information
     &, rain_ice, rain_clear                                            &
     &, rho, rhor                                                       &
     &, tcg,tcgi,tcgc,tcgci,tcgg,tcggi                                  &
                                          ! Parametrization information
     &, corr, corr2, rocor, dhi, dhir                                   &
                                          ! Parametrization information
     &, cx, constp, ai, bi, aic, bic, lfrcp                             &
                                          ! Microphysical information
     &, l_update_cf, l_seq                                              &
                                          ! Code options
     &, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_psd, l_crystals        &
                                              ! Microphysics options
     &, pifall, psfall, prfall, pgfall                                  &
                                          ! Mass transfer diagnostics
     &, pimlt, psmlt, pgmlt                                             &
                                          ! Mass transfer diagnostics
     &, iterations, n_melt_iterations                                   &
                                          ! Iterations of microphysics
     &, cftransfer,cfftransfer,rftransfer                               &
                                          ! Cloud transfer diagnostics
     &  )
!
      Implicit None
!
! Purpose:
!   Call fall out and melting terms in an iterative way
!
! Method:
!   Iterate over the fall-out and melting terms for a prescribed
!   number of iterations.
!
! Current Owner of Code: Jonathan Wilkinson
!
! History:
! Version   Date     Comment
!  6.2    19-12-05   Original Code (Damian Wilson)
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
! Start C_LSPSIZ
! Description: Include file containing idealised forcing options
! Author:      R. Forbes
!
! History:
! Version  Date      Comment
! -------  ----      -------
!   6.1    01/08/04  Increase dimension for rain/graupel.  R.Forbes
!   6.2    22/08/05  Include the step size between ice categories.
!                                                   Damian Wilson

! Sets up the size of arrays for CX and CONSTP
      REAL CX(100),CONSTP(100)
      INTEGER,PARAMETER:: ice_type_offset=20

! End C_LSPSIZ






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
      Integer, Intent(In) ::                                            &
     &  points                                                          &
                          ! Number of points to calculate
     &, iterations                                                      &
                          ! Number of microphysics iterations
     &, n_melt_iterations ! Number of iterative melting iterations
!
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
     &, q(points)                                                       &
                          ! Vapour content / kg kg-1
     &, q_ice(points)                                                   &
                          ! Vapour content in ice only region / kg kg-1
     &, qsl(points)                                                     &
                          ! Sat. vapour content wrt liquid / kg kg-1
     &, p(points)                                                       &
                          ! Pressure / N m-2
     &, frac_agg(points)                                                &
                          ! Fraction of ice mass that is aggregates
     &, snow_cry(points)                                                &
                          ! Crystal flux into layer / kg m-2 s-1
     &, snow_agg(points)                                                &
                          ! Aggregate flux into layer / kg m-2 s-1
     &, rainrate(points)                                                &
                          ! Rain flux into layer / kg m-2 s-1
     &, grauprate(points)                                               &
                          ! Graupel flux into layer / kg m-2 s-1
     &, area_liq(points)                                                &
                          ! Frac of gridbox with liq cloud but not ice
     &, area_mix(points)                                                &
                          ! Fraction of gridbox with ice and liq cloud
     &, area_clear(points)                                              &
                          ! Fraction of gridbox with clear sky
     &, area_ice(points)                                                &
                          ! Frac of gridbox with ice cloud but not liq
     &, cfice(points)                                                   &
                          ! Fraction of gridbox with ice cloud
     &, cficei(points)                                                  &
                          ! 1/Fraction of gridbox with ice cloud
     &, frac_ice_above(points)                                          &
                               ! Ice cloud fraction in layer above
     &, rainfrac(points)                                                &
                          ! Fraction of gridbox with rain
     &, rain_liq(points)                                                &
                          ! Fraction of gridbox with rain and liq cloud
     &, rain_mix(points)                                                &
                          ! Frac of gbox with rain and mixed phase cloud
     &, rain_ice(points)                                                &
                          ! Fraction of gridbox with rain and ice cloud
     &, rain_clear(points)                                              &
                          ! Fraction of gridbox with rain but no cloud
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, rhor(points)                                                    &
                          ! 1 / Air density / m3 kg-1
     &, dhi(points)                                                     &
                          ! Timestep / thickness of model layer / s m-1
     &, dhir(points)                                                    &
                          ! 1/dhi / m s-1
!    &, m0                ! Seed ice water / kg kg-1 (in c_lspmic)
!    &, wind_shear_factor ! Degree of overhang of ice fraction between
                          ! two model layers (in c_lspmic)
     &, tcg(points)                                                     & 
                          ! T dependent func in ice agg. size dist'n
     &, tcgi(points)                                                    &
                          ! 1/tcg (no units)
     &, tcgc(points)                                                    & 
                          ! T dependent func in ice crystal size dist'n
     &, tcgci(points)                                                   &
                          ! 1/tcgc (no units)
     &, tcgg(points)                                                    & 
                          ! T dependent func in graupel size dist'n
     &, tcggi(points)                                                   &
                          ! 1/tcgg (no units)
     &, corr(points)                                                    &
                          ! Air density fall speed correction (no units)
     &, corr2(points)                                                   &
                          ! Temperature correcn to diffusivity (no units
     &, rocor(points)                                                   &
                          ! sqrt(rho*corr*corr2)
     &, lfrcp                                                           &
                          ! Latent heat of fusion
                          !   / heat capacity of air / K
     &, ai, bi, aic, bic
                          ! Ice mass size relationships. m(D)=ai D^bi
!
      Real, Intent(InOut) ::                                            &
     &  T(points)                                                       &
                          ! Temperature / K
     &, qcf_cry(points)                                                 &
                          ! Ice crystal mixing ratio / kg kg-1
     &, qcf_agg(points)                                                 &
                          ! Ice aggregate mixing ratio / kg kg-1
     &, qrain(points)                                                   &
                          ! Rain mixing ratio / kg kg-1
     &, qgraup(points)                                                  &
                          ! Graupel mixing ratio / kg kg-1
     &, snowt_cry(points)                                               &
                          ! Snowfall rate out of this layer / kg m-2 s-1
     &, snowt_agg(points)                                               &
                            ! for crystals and aggregates
     &, rainratet(points)                                               &
                          ! Rain rate out of this layer / kg m-2 s-1
     &, graupratet(points)                                              &
                          ! Graupel rate out of this layer / kg m-2 s-1
     &, vf_cry(points)                                                  &
                          ! On input: Fall speed of hydrometeors
     &, vf_agg(points)                                                  &
                                    ! entering the current layer / m s-1
     &, vf_rain(points)                                                 &
                          ! On output: Fall speed of hydrometeors
     &, vf_graup(points)                                                &
                                    ! leaving the current layer / m s-1
     &, cf(points)                                                      &
                          ! Current cloud fraction
     &, cfl(points)                                                     &
                          ! Current liquid cloud fraction
     &, cff(points)       ! Current ice cloud fraction

      Real, Intent(InOut) ::                                            &
     &  pifall(points)                                                  &
                          ! Rate of change of crystal, aggregate,
     &, psfall(points)                                                  &
                          ! rain and graupel mixing ratios due
     &, prfall(points)                                                  &
                          ! to sedimentation / kg kg-1 s-1
     &, pgfall(points)                                                  &
     &, pimlt(points)                                                   &
                          ! Melting rate of crystals / kg kg-1 s-1
     &, psmlt(points)                                                   &
                          ! Melting rate of snow / kg kg-1 s-1
     &, pgmlt(points)                                                   &
                          ! Melting rate of graupel / kg kg-1 s-1
     &, cftransfer(points)                                              &
                          ! Cloud fraction increment this tstep
     &, cfftransfer(points)                                             &
                           ! Ice cloud fraction inc this tstep
     &, rftransfer(points)! Rain fraction increment this tstep
!
      Logical, Intent(In) ::                                            &
     &  l_update_cf                                                     &
                          ! Update cloud fractions
     &, l_seq                                                           &
                          ! Carry out sequential updating
     &, l_mcr_qcf2                                                      &
                          ! Two ice prognostics are active
     &, l_mcr_qrain                                                     &
                          ! Prognostic rain is active
     &, l_mcr_qgraup                                                    &
                          ! Graupel is active
     &, l_psd                                                           &
                          ! Use generic ice particle size distribution
     &, l_crystals
                          ! Use the ice crystal category
!
! Local Variables
!
      Integer                                                           &
     &  i                                                               &
                          ! Loop counter for points
     &, k                                                               &
                          ! Loop counter for melting iterations
     &, total_iterations  ! Net number of iterations
!
      Real                                                              &
     &  qcft(points)                                                    &
                          ! Total ice content / kg kg-1
     &, dhi_it(points)                                                  &
                          ! Iterated tstep / layer thick. / s m-1
     &, dhir_it(points)                                                 &
                          ! Layer thick. / iterated tstep / m s-1
     &, timestep_it       ! Iterated timestep / s
!
      !-----------------------------------------------
      ! Set up short timestep
      !-----------------------------------------------
      timestep_it = timestep / n_melt_iterations
      total_iterations = iterations * n_melt_iterations

      Do i = 1, points
        ! Update dhi values for short timestep
        dhi_it(i) = dhi(i) / n_melt_iterations
        dhir_it(i) = dhir(i) * n_melt_iterations
      End Do

      !-----------------------------------------------
      ! Call fall and melting for the iterated points
      !-----------------------------------------------
      Do k = 1, n_melt_iterations

! DEPENDS ON: lsp_fall
        Call lsp_fall(points, timestep_it                               &
     &,                  qcf_cry, qcf_agg, frac_agg, qrain, qgraup, T   &
     &,                  snow_agg, snow_cry, rainrate, grauprate        &
     &,                  snowt_agg, snowt_cry, rainratet, graupratet    &
     &,                  vf_agg, vf_cry, vf_rain, vf_graup              &
     &,                  area_clear, area_ice, cfice, cficei            &
     &,                  frac_ice_above, cf, cfl, cff                   &
     &,                  rho, rhor, tcgi, tcgci                         &
     &,                  corr, dhi_it, dhir_it, cx, constp              &
     &,                  ai, bi, aic, bic                               &
     &,                  l_update_cf, l_seq                             &
     &,                  l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup          &
     &,                  l_psd, l_crystals                              &
     &,                  pifall, psfall, prfall, pgfall                 &
     &,                  total_iterations                               &
     &,                  cftransfer, cfftransfer                        &
     &                  )

        Do i=1,points
          ! Calculate total ice content
          qcft(i) = qcf_cry(i) + qcf_agg(i)
        End do

        If (l_crystals) then
! DEPENDS ON: lsp_melting
          Call lsp_melting(points, timestep_it                          &
     &,                   q, q_ice, qcf_cry, qcft, qrain, qsl, T, p     &
     &,                   area_liq, area_mix, area_ice                  &
     &,                   cfice, cficei, frac_ice_above                 &
     &,                   rainfrac, rain_liq, rain_mix                  &
     &,                   rain_ice, rain_clear, cf, cff                 &
     &,                   rho, rhor, m0, tcg, tcgi, corr2, rocor        &
     &,                   cx, constp, ai, bi, aic, bic, lfrcp, 0        &
     &,                   l_update_cf, l_seq, .false.                   &
     &,                   pimlt, total_iterations                       &
     &,                   cftransfer, cfftransfer, rftransfer           &
     &                   )
        End if  ! l_crystals

        Do i=1,points
          ! Calculate total ice content
          qcft(i) = qcf_cry(i) + qcf_agg(i)
        End do

! DEPENDS ON: lsp_melting
        Call lsp_melting(points, timestep_it                            &
     &,                   q, q_ice, qcf_agg, qcft, qrain, qsl, T, p     &
     &,                   area_liq, area_mix, area_ice                  &
     &,                   cfice, cficei, frac_ice_above                 &
     &,                   rainfrac, rain_liq, rain_mix                  &
     &,                   rain_ice, rain_clear, cf, cff                 &
     &,                   rho, rhor, m0, tcgc, tcgci, corr2, rocor      &
     &,                   cx, constp, ai, bi, aic, bic, lfrcp, 1        &
     &,                   l_update_cf, l_seq, l_psd                     &
     &,                   psmlt, total_iterations                       &
     &,                   cftransfer, cfftransfer, rftransfer           &
     &                   )

        If (L_mcr_qgraup) Then
          ! Graupel does not update cloud fractions so there is no need
          ! to update qcft (it is not used)
! DEPENDS ON: lsp_melting
          Call lsp_melting(points, timestep_it                          &
     &,                     q, q_ice, qgraup, qcft, qrain, qsl, T, p    &
     &,                     area_liq, area_mix, area_ice                &
     &,                     cfice, cficei, frac_ice_above               &
     &,                     rainfrac, rain_liq, rain_mix                &
     &,                     rain_ice, rain_clear, cf, cff               &
     &,                     rho, rhor, m0, tcgg, tcggi, corr2, rocor    &
     &,                     cx, constp, ai, bi, aic, bic, lfrcp, 3      &
     &,                     .false., l_seq, .false.                     &
     &,                     pgmlt, total_iterations                     &
     &,                     cftransfer, cfftransfer, rftransfer         &
     &                     )
        End If  ! L_mcr_qgraup

      End Do  ! n_melt_iterations

      Return  ! End of the subroutine
      END SUBROUTINE LSP_IT_FALL_MELT
