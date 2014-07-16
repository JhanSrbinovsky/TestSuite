
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Large-scale precipitation scheme. Advection through falling of ice
! and rain
! Subroutine Interface:
      SUBROUTINE LSP_FALL(                                              &
     &  points, timestep                                                &
                                          ! Number of points and tstep
     &, qcf_cry_0, qcf_agg_0, frac_agg                                  &
                                          ! Ice contents
     &, qrain_0,qgraup_0, T                                             &
                                          ! Rain and graupel contents
     &, snow_agg,snow_cry,rainrate,grauprate                            &
                                             ! Sedimentation into layer
     &, snowt_agg,snowt_cry,rainratet,graupratet                        &
                                                 ! Sedim. out of layer
     &, vf_agg, vf_cry, vf_rain, vf_graup                               &
                                          ! Fall speeds of hydrometeors
     &, area_clear,area_ice,cfice,cficei                                &
                                          ! Cloud fraction information
     &, frac_ice_above                                                  &
                                            ! at start of microphy. ts
     &, cf, cfl, cff                                                    &
                                          ! Current cloud fractions for
                                          ! updating
     &, rho, rhor, tcgi, tcgci                                          &
                                          ! Parametrization information
     &, corr, dhi, dhir                                                 &
                                          ! Parametrization information
     &, cx, constp, ai, bi, aic, bic                                    &
                                          ! Microphysical information
     &, l_update_cf, l_seq                                              &
                                          ! Code options
     &, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_psd, l_crystals        &
                                              ! Microphysics options
     &, d_qcf_cry_dt,d_qcf_agg_dt                                       &
                                          ! Mass transfer diagnostics
     &, d_qrain_dt,d_qgraup_dt                                          &
                                          ! Mass transfer diagnostics
     &, iterations                                                      &
                                          ! Iterations of microphysics
     &, cftransfer,cfftransfer                                          &
                                          ! Cloud transfer diagnostics
     &  )
!
      Implicit None
!
! Purpose:
!   Update cloud prognostics as a result of ice particle and
!   raindrop fall
!
! Method:
!   Calculate particle fall speeds and solve the advection equation
!   for mixing ratios following Rotstayns method.
!
! Current Owner of Code: Jonathan Wilkinson
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: UMDP 26.
!
! Both small and large ice, and raindrops, will fall, hence require
! advection downwards. We include both code for the two-ice prognostics
! and the single ice prognostic that is diagnostically split. Although
! the advection methods are very similar, they are different enough
! (in their calculation of fall-speed of the ice from above) to
! use different branches of code.
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
     &, iterations        ! Number of microphysics iterations
                          ! (for rate diagnostic)
!
      Real, Intent(In) ::                                               &
     &  Timestep                                                        &
                          ! Timestep / s
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
     &, area_clear(points)                                              &
                          ! Fraction of gridbox with clear sky
     &, area_ice(points)                                                &
                          ! Frac of gridbox with ice but not liquid
     &, cfice(points)                                                   &
                          ! Fraction of gridbox with ice cloud
     &, cficei(points)                                                  &
                          ! 1/Fraction of gridbox with ice cloud
     &, frac_ice_above(points)                                          &
                               ! Ice cloud fraction in layer above
     &, rho(points)                                                     &
                          ! Air density / kg m-3
     &, rhor(points)                                                    &
                          ! 1 / Air density / m3 kg-1
     &, dhi(points)                                                     &
                          ! Timestep / thickness of model layer / s m-1
     &, dhir(points)                                                    &
                          ! 1/dhi / m s-1
!    &, m0                ! Seed ice water content / kg kg-1 (c_lspmic)
!    &, wind_shear_factor ! Degree of overhang of ice fraction between
                          ! two model layers (in c_lspmic)
!    &, tcg(points)       ! T dependent func. in ice agg. size dist'n
     &, tcgi(points)                                                    &
                          ! 1/tcg (no units)
!    &, tcgc(points)      ! T dependent func. in ice crystal size dist'n
     &, tcgci(points)                                                   &
                          ! 1/tcgc (no units)
     &, corr(points)                                                    &
                          ! Air density fall speed correction (no units)
     &, ai, bi, aic, bic
                          ! Ice mass size relationship m(D)=ai D^bi
!
      Real, Intent(InOut) ::                                            &
     &  qcf_cry_0(points)                                               &
                          ! Ice crystal mixing ratio / kg kg-1
     &, qcf_agg_0(points)                                               &
                          ! Ice aggregate mixing ratio / kg kg-1
     &, qrain_0(points)                                                 &
                          ! Rain mixing ratio / kg kg-1
     &, qgraup_0(points)                                                &
                          ! Graupel mixing ratio / kg kg-1
     &, T(points)                                                       &
                          ! Temperature / K
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
!
      Real, Intent(InOut) ::                                            &
     &  d_qcf_cry_dt(points)                                            &
                             ! Rate of change of crystal, aggregate,
     &, d_qcf_agg_dt(points)                                            &
                               ! rain and graupel mixing ratios due
     &, d_qrain_dt(points)                                              &
                               ! to sedimentation / kg kg-1 s-1
     &, d_qgraup_dt(points)                                             &
     &, cftransfer(points)                                              &
                           ! Cloud fraction increment this tstep
     &, cfftransfer(points)! Ice cloud fraction inc this tstep
!
      Logical, Intent(In) ::                                            &
     &  l_update_cf                                                     &
                          ! Update cloud fractions
     &, l_seq                                                           &
                          ! Carry out sequential updating
     &, l_psd                                                           &
                          ! Use generic ice particle size distribution
     &, l_crystals                                                      &
                          ! Use crystals ice category
     &, l_mcr_qcf2                                                      &
                          ! Two ice prognostics are active
     &, l_mcr_qrain                                                     &
                          ! Prognostic rain is active
     &, l_mcr_qgraup      ! Graupel is active
!
! Local Variables
!
      Integer                                                           &
     &  i                 ! Loop counter for points

      Real                                                              &
     &  qcf_cry(points)                                                 &
                          ! Working copy of qcf_cry_0 / kg kg-1
     &, qcf_agg(points)                                                 &
                          ! Working copy of qcf_agg_0 / kg kg-1
     &, qrain(points)                                                   &
                          ! Working copy of q_rain / kg kg-1
     &, qgraup(points)                                                  &
                          ! Working copy of q_graup / kg kg-1
     &, fqi_agg(points)                                                 &
                          ! Fall speed of aggregates out of the
                          ! current layer / m s-1
     &, fqi_cry(points)                                                 &
                          ! Fall speed of crystals out of the
                          ! current layer / m s-1
     &, fqirqi2_agg(points)                                             &
                            ! Fraction of aggregate and crystal mass
     &, fqirqi2_cry(points)                                             &
                              ! that remains in the current layer
                              ! after sedimentation
     &, fqirqi_agg                                                      &
                          ! Flux of aggregates out of layer / m s-1
     &, fqirqi_cry                                                      &
                          ! Flux of crystals out of layer / m s-1
     &, fqi_rain(points)                                                &
                          ! Bulk fall speed of rain and graupel
     &, fqi_graup(points)                                               &
                            ! out of the current layer / m s-1
     &, MixRatio_FromAbove(points)                                      &
                                   ! Estimate of the mixing ratio
                          ! of ice in the layer above / kg kg-1
     &, temp3(points)                                                   &
                          ! Fraction of layer ice falls through
     &, overhang                                                        &
                          ! Ice cloud fraction that overhangs
                          ! the current layer from above
     &, deltacf(points)                                                 &
                          ! Change in cloud fraction across timestep
     &, deltacff(points)                                                &
                          ! Change in ice cloud fraction across tstep
     &, m_bi_di(points)
                          ! bi+di'th moment of the ice particle dist.
!!fra298
      real save_vf_agg(points), save_vf_cry(points)

      Logical                                                           &
     &  l_separate_advection

      !-----------------------------------------------
      ! Set advection code method
      !-----------------------------------------------
      l_separate_advection = l_mcr_qcf2

      !-----------------------------------------------
      ! Calculate moment of size distribution if appropriate 
      !----------------------------------------------- 
      If (l_psd) then
        ! Calculate the bi+di (cx(82)) moment of the
        ! ice particle size distribution 
! DEPENDS ON: lsp_moments 
        Call lsp_moments(points,rho,T,qcf_agg_0,cficei,                 & 
     &                   ai,bi,cx(82),m_bi_di) 
      End if

      Do i=1,points
        !-----------------------------------------------
        ! Copy input water contents (allows sequential updates)
        !-----------------------------------------------
        qcf_cry(i) = qcf_cry_0(i)
        qcf_agg(i) = qcf_agg_0(i)
        qrain(i)   = qrain_0(i)
        qgraup(i)  = qgraup_0(i)
!!fra298
        ! Save the fall speeds coming in from layer above before 
        ! they get overwritten with fall speed passed to layer below.
        ! Will be used in updating falling ice cloud fraction
        save_vf_agg(i)=vf_agg(i)
        save_vf_cry(i)=vf_cry(i)

        If (qcf_agg(i) >  m0) Then
          !-----------------------------------------------
          ! Estimate fall speed out of this layer
          !-----------------------------------------------
          If (l_psd) then
            ! Use the generic PSD
            ! constp(82) = ci*ai
            fqi_agg(i) = constp(82) * corr(i) * m_bi_di(i)              &
     &                 / (rho(i) * qcf_agg(i) * cficei(i))
          Else
            fqi_agg(i) = constp(24) * corr(i) *                         &
     &      (rho(i)*qcf_agg(i)*constp(25)*tcgi(i)*cficei(i))**cx(23)
          End if  ! l_psd

        Else
          !-----------------------------------------------
          ! Ice content is small so it is numerically best to
          ! assume there is no ice here, so the fall speed is zero.
          !-----------------------------------------------
          fqi_agg(i)=0.0
        End if  ! qcf_agg gt 0

        If (qcf_cry(i) >  m0 .and. l_crystals) Then
         !-----------------------------------------------
         ! Estimate fall speed out of this layer
         !-----------------------------------------------
          fqi_cry(i) = constp(4) * corr(i)*                             &
     &    (rho(i)*qcf_cry(i)*constp(5)*tcgci(i)*cficei(i))**cx(3)
        Else
          !-----------------------------------------------
          ! Ice content is small so it is numerically best to
          ! assume there is no ice here, so the fall speed is zero.
          !-----------------------------------------------
          fqi_cry(i)=0.0
        End if  ! qcf_cry gt 0

      End Do  ! Points

      If (.not. l_separate_advection) Then

        Do i=1,points
          !-----------------------------------------------
          ! Calculate fall speed from above as an average of the
          ! crystal and aggregate fall speeds.
          !-----------------------------------------------
          MixRatio_FromAbove(i) = (snow_agg(i) + snow_cry(i))           &
     &                             * dhi(i)*rhor(i)
          If ((qcf_cry(i)+qcf_agg(i)+MixRatio_FromAbove(i)) >  m0) Then

            !-----------------------------------------------
            ! Make a linear combination of fall speeds between this
            ! layer and the layer above to aid numerical solution
            !-----------------------------------------------
            fqi_cry(i)= ( fqi_cry(i)*(qcf_cry(i)+qcf_agg(i))            &
     &                    + vf_agg(i)*MixRatio_FromAbove(i) )           &
     &                  / (qcf_cry(i)+qcf_agg(i)+MixRatio_FromAbove(i))
            fqi_agg(i)= ( fqi_agg(i)*(qcf_cry(i)+qcf_agg(i))            &
     &                    + vf_agg(i)*MixRatio_FromAbove(i) )           &
     &                  / (qcf_cry(i)+qcf_agg(i)+MixRatio_FromAbove(i))

            !-----------------------------------------------
            ! Fall speed of ice to pass to layer below
            !-----------------------------------------------
            vf_agg(i) = fqi_cry(i) * (1.0-frac_agg(i))                  &
     &                + fqi_agg(i) *      frac_agg(i)

          Else  ! qcf gt m0
            ! Little ice so set fall speed to zero
            vf_agg(i)=0.0

          End if  ! qcf gt m0

          !-----------------------------------------------
          ! Solve for fraction of ice that remains in the same layer
          !-----------------------------------------------
          fqirqi2_agg(i) = exp(-fqi_agg(i)*dhi(i))
          fqirqi2_cry(i) = exp(-fqi_cry(i)*dhi(i))

          If (fqi_agg(i)  >   0.0) Then
            !-----------------------------------------------
            ! Advect aggregates
            !-----------------------------------------------
            fqirqi_agg = snow_agg(i) + dhir(i) *                        &
     &                  (rho(i)*qcf_agg(i)-snow_agg(i)/fqi_agg(i))      &
     &                  * (1.0-fqirqi2_agg(i))
            qcf_agg(i) = snow_agg(i) * rhor(i) / fqi_agg(i)             &
     &                  * (1.0-fqirqi2_agg(i))                          &
     &                  + qcf_agg(i)*fqirqi2_agg(i)
          Else
            !-----------------------------------------------
            ! No fall of ice out of the layer
            !-----------------------------------------------
            fqirqi_agg = 0.0
            qcf_agg(i) = snow_agg(i)*rhor(i)*dhi(i)
          End if  ! fqi_agg gt 0

          If (fqi_cry(i)  >   0.0 .and. l_crystals) Then
            !-----------------------------------------------
            ! Advect crystals
            !-----------------------------------------------
            fqirqi_cry = snow_cry(i) + dhir(i) *                        &
     &                  (rho(i)*qcf_cry(i)-snow_cry(i)/fqi_cry(i))      &
     &                  * (1.0-fqirqi2_cry(i))
            qcf_cry(i) = snow_cry(i) * rhor(i) / fqi_cry(i)             &
     &                  * (1.0-fqirqi2_cry(i))                          &
     &                  + qcf_cry(i)*fqirqi2_cry(i)
          Else
            !-----------------------------------------------
            ! No fall of ice out of the layer
            !-----------------------------------------------
            fqirqi_cry = 0.0
            qcf_cry(i) = snow_cry(i)*rhor(i)*dhi(i)
          End if  ! fqi_cry gt 0

          ! --------------------------------------------------
          ! Snow is used to save flux out of layer
          ! --------------------------------------------------
          snowt_cry(i) = snowt_cry(i) + fqirqi_cry/iterations
          snowt_agg(i) = snowt_agg(i) + fqirqi_agg/iterations

        End Do  ! Points

      Else  ! l_separate_advection
        !--------------------------------------------------
        ! Prognostic Ice Crystal and Snow Aggregate Sedimentation
        !--------------------------------------------------
        ! Ice Crystals
        If (l_crystals) then
! DEPENDS ON: lsp_sedim_eulexp
          Call lsp_sedim_eulexp(                                        &
     &      iterations, points, m0, dhi, dhir, rho, rhor                &
     &,     snow_cry, fqi_cry, qcf_cry, vf_cry                          &
     &,     snowt_cry)
        End if  ! l_crystals

        ! Aggregates
! DEPENDS ON: lsp_sedim_eulexp
        Call lsp_sedim_eulexp(                                          &
     &    iterations, points, m0, dhi, dhir, rho, rhor                  &
     &,   snow_agg, fqi_agg, qcf_agg, vf_agg                            &
     &,   snowt_agg)

      End if ! L_separate_advection

      !--------------------------------------------------
      ! Rain Sedimentation
      !--------------------------------------------------
      If (L_mcr_qrain) Then

        Do i=1,points

          If (qrain(i)  >   m0) Then
            !--------------------------------------------------
            ! Estimate fall speed out of this layer (FQI_RAIN)
            !--------------------------------------------------
            fqi_rain(i) = constp(41) * corr(i) / 6.0 *                  &
     &        (rho(i) * qrain(i) * constp(50)) ** CX(51)
              ! The rain flux should be divided by the rain fraction,
              ! but this causes problems with the way rain fraction
              ! is currently formulated. Commented out for now.
!    &        (rho(i) * qrain(i) * constp(50) / rainfrac(i)) ** CX(51)
          Else
            ! Fall speed is set to zero
            fqi_rain(i) = 0.0
          End If  ! qrain gt m0

        End Do  ! Points

! DEPENDS ON: lsp_sedim_eulexp
        Call lsp_sedim_eulexp(                                          &
     &    iterations, points, m0, dhi, dhir, rho, rhor                  &
     &,   rainrate, fqi_rain, qrain, vf_rain                            &
     &,   rainratet)

      End if  ! L_mcr_qrain

      !--------------------------------------------------
      ! Graupel Sedimentation
      !--------------------------------------------------
      If (L_mcr_qgraup) Then

        Do i=1,points

          If (qgraup(i)  >   m0) then
            !--------------------------------------------------
            ! Estimate fall speed out of this layer (FQI_GRAUP)
            !--------------------------------------------------
            fqi_graup(i) = constp(64) * corr(i) *                       &
     &        (rho(i) * qgraup(i) * constp(65)) ** cx(63)
          Else
            ! Fall speed is set to zero
            fqi_graup(i) = 0.0
          End If  ! qgraup gt m0

        End Do  ! Points

! DEPENDS ON: lsp_sedim_eulexp
        Call lsp_sedim_eulexp(                                          &
     &     iterations, points, m0, dhi, dhir, rho, rhor                 &
     &,    grauprate, fqi_graup, qgraup, vf_graup                       &
     &,    graupratet)

      End if  !  L_mcr_qgraup
      
      Do i=1,points

        If (l_update_cf) Then
          !-----------------------------------------------
          ! Update cloud fractions
          !-----------------------------------------------

          If ((qcf_cry(i)+qcf_agg(i))  >   0.0) Then
            !-----------------------------------------------
            ! Calculate fraction of a layer the ice has fallen
            !-----------------------------------------------
            If (L_mcr_qcf2) then
!!fra298              temp3(i) = dhi(i) * (vf_agg(i)*qcf_agg(i)                 &
!!     &            + vf_cry(i)*qcf_cry(i))/(qcf_cry(i)+qcf_agg(i))
              temp3(i) = dhi(i) * (save_vf_agg(i)*qcf_agg(i)          &
                  + save_vf_cry(i)*qcf_cry(i))/(qcf_cry(i)+qcf_agg(i)) 
            Else
!!              temp3(i) = dhi(i) * vf_agg(i)
              temp3(i) = dhi(i) * save_vf_agg(i)
            End if
            !-----------------------------------------------
            ! Calculate the amount of cloud overhang between levels
            !-----------------------------------------------
            overhang = max(frac_ice_above(i)-cff(i),0.0)
!!fra298            If (temp3(i)  >   0.0) Then
!!              overhang = overhang + wind_shear_factor / temp3(i)        &
!!     &                          * timestep
!!            End if

            !-----------------------------------------------
            ! Physially limit the amount of overhang. Note there is
            ! no limit on the fraction of ice in the layer above
            !-----------------------------------------------
            overhang = min(overhang,1.0-cff(i))
            ! Now limit the fall out quantity
            temp3(i) = min(max(temp3(i),0.0),1.0)

            !-----------------------------------------------
            ! Calculate change in ice cloud fraction
            !-----------------------------------------------
            
            temp3(i) = temp3(i) * overhang
            deltacff(i) = temp3(i)

            If (temp3(i)  <=  0.0) Then
              !-----------------------------------------------
              ! Total cloud fraction will be reduced
              !-----------------------------------------------
              deltacf(i) = temp3(i)*area_ice(i)*cficei(i) !Random o'lap
            Else If (cfice(i)  <   1.0) Then
              !-----------------------------------------------
              ! Total cloud fraction will be increased
              !-----------------------------------------------
!             deltacf(i) = temp3(i)*area_clear(i) /(1.0-cfice(i))
!                                                        !Random o'lap
              deltacf(i) = (min(temp3(i),area_clear(i))) !Minimum o'lap
            End if

          Else  ! qcf gt 0
            ! Set ice cloud fraction to zero and total cloud
            ! fraction to the liquid cloud fraction
            deltacff(i) = -cff(i)
            deltacf(i)  = (cfl(i) - cf(i))
          End if  ! qcf gt 0

        End if  ! l_update_cf

        !-----------------------------------------------
        ! Form transfer diagnostics
        !-----------------------------------------------
        
        d_qcf_cry_dt(i) = d_qcf_cry_dt(i) +                             &
     &    (qcf_cry(i) - qcf_cry_0(i)) / (timestep*iterations)
        d_qcf_agg_dt(i) = d_qcf_agg_dt(i) +                             &
     &    (qcf_agg(i) - qcf_agg_0(i)) / (timestep*iterations)
        d_qrain_dt(i) = d_qrain_dt(i) +                                 &
     &    (qrain(i) - qrain_0(i)) / (timestep*iterations)
        d_qgraup_dt(i) = d_qgraup_dt(i) +                               &
     &    (qgraup(i) - qgraup_0(i)) / (timestep*iterations)

        !-----------------------------------------------
        ! Update water contents
        !-----------------------------------------------
        If (l_seq) then
          qcf_cry_0(i) = qcf_cry(i)
          qcf_agg_0(i) = qcf_agg(i)
          qrain_0(i)   = qrain(i)
          qgraup_0(i)  = qgraup(i)
        End if

        !-----------------------------------------------
        ! Update cloud fractions
        !-----------------------------------------------
        If (l_update_cf) Then
          If (l_seq) then
            cf(i) = cf(i) + deltacf(i)
            cff(i) = cff(i) + deltacff(i)
          End if
          cftransfer(i) = cftransfer(i) + deltacf(i)                    &
     &                    /(timestep*iterations)
          cfftransfer(i)= cfftransfer(i)+ deltacff(i)                   &
     &                    /(timestep*iterations)
        End if

      End Do  ! Points

      Return  ! End of the subroutine
      END SUBROUTINE LSP_FALL
