#if defined(A04_3B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE LSP_ICE------------------------------------------------
!
!  Purpose: Form or augment ice at the expense of cloud water or
!           vapour in one model layer.
!           Also perform flux divergence of falling ice and rain,
!           Evaporation and melting of snow,
!           Evaporation of rain, formation of rain.
!           This is the principal subroutine of the 3B large scale
!           precipitation scheme.
!
! S Ballard   <- programmer
! D Wilson    <- programmer
!
! Current Owner of Code: Jonathan Wilkinson
!
!  Documentation: Unified Model Documentation Paper No 26.
!
!  Arguments:-----------------------------------------------------------
      SUBROUTINE LSP_ICE(                                               &
     &  P,RHODZ, deltaz, rhodz_dry, rhodz_moist,                        &
     &  TIMESTEPFIXED,n_iterations,POINTS,L_MURK,                       &
     &  RHCRIT, L_USE_SULPHATE_AUTOCONV, L_AUTO_DEBIAS,                 &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_psd, L_it_melting,     &
     &  l_non_hydrostatic, l_mixing_ratio, l_cry_agg_dep,               &
     &  SO4_ACC,SO4_DIS,SO4_AIT,                                        &
     &  L_BIOMASS_CCN, BMASS_AGD, BMASS_CLD,                            &
     &  L_OCFF_CCN, OCFF_AGD, OCFF_CLD,                                 &
     &  L_SEASALT_CCN, SEA_SALT_FILM, SEA_SALT_JET, SALT_DIM_ICE,       &
     &  L_biogenic_CCN, biogenic, biog_dim_ice,                         &
     &  SNOW_DEPTH,LAND_FRACT,AEROSOL,                                  &
     &  L_seq_mcr,L_autoc_3b,L_autolim_3b,L_autoconv_murk,              &
     &  l_droplet_settle, ec_auto,N_drop_land,                          &
     &  N_drop_sea,N_drop_land_cr,N_drop_sea_cr,Ntot_land, Ntot_sea,    &
     &  ai, bi, aic, bic,                                               &
     &  QCF,QCL,Q,QCF2,QRAIN,QGRAUP,                                    &
     &  RAIN, VF_RAIN, SNOW, VF,                                        &
     &  SNOW2, VF2, GRAUPRATE, VF_GRAUP, droplet_flux,                  &
     &  FRAC_ICE_ABOVE,CTTEMP,RAINFRAC,                                 &
     &  T,CFKEEP,CFLIQ,CFICE,BLAND,CX,CONSTP                            &
     &, PSDEP,PSAUT,PSACW,PSACR,PSACI,PSMLT,PSMLTEVP                    &
     &, PRAUT,PRACW,PREVP                                               &
     &, PGAUT,PGACW,PGACS,PGMLT                                         &
     &, PIFRW,PIPRM,PIDEP,PIACW,PIACR,PIMLT,PIMLTEVP                    &
     &, PIFALL,PSFALL,PRFALL,PGFALL,PLSET,PLEVPSET                      &
     & )
      IMPLICIT NONE
!
#include "c_0_dg_c.h"
! c_lspdif will call c_r_cp, c_epslon and c_lheat
#include "c_lspmic.h"
#include "c_lspdif.h"

!
      INTEGER                                                           &
                      !, INTENT(IN)
     & POINTS                                                           &
!        Number of points to be processed.
     &, n_iterations                                                    &
                               ! Dummy in 3B version of the code
     &,SALT_DIM_ICE                                                     &
!        Number of points for sea-salt arrays (either POINTS or 1)
     &,biog_dim_ice
!        Number of points for sea-salt arrays (either POINTS or 1)

      REAL                                                              &
                      !, INTENT(IN)
     & TIMESTEPFIXED
!        Timestep of physics in model (s).
!
      REAL                                                              &
                      !, INTENT(IN)
     &  CFLIQ(POINTS)                                                   &
!         Liquid cloud fraction in this layer.
     &, CFICE(POINTS)                                                   &
!         Frozen cloud fraction in this layer.
     &, P(POINTS)                                                       &
!         Air pressure at this level (Pa).
     &, RHODZ(POINTS)                                                   &
!         Air mass p.u.a. in this layer (kg per sq m).
     &, deltaz(points)                                                  &
                           ! Depth of layer / m
     &, rhodz_dry(points)                                               &
                           ! Density of dry air / kg m-3
     &, rhodz_moist(points)! Density of moist air / kg m-3
!
      REAL                                                              &
                      !, INTENT(INOUT)
     & Q(POINTS)                                                        &
!        Specific humidity at this level (kg wat per kg air).
     &,QCF(POINTS)                                                      &
!        Cloud ice (kg water per kg air).
     &,QCL(POINTS)                                                      &
!        Cloud liquid water (kg water per kg air).
     &, QCF2(POINTS)                                                    &
                       ! Dummy in 3B version of the code
     &, QRAIN(POINTS)                                                   &
                       ! Dummy in 3B version of the code
     &, QGRAUP(POINTS)                                                  &
                       ! Dummy in 3B version of the code
     &,T(POINTS)                                                        &
!        Temperature at this level (K).
     &,RAIN(POINTS)                                                     &
!        On input: Rate of rainfall entering this layer from above.
!        On output: Rate of rainfall leaving this layer.
!                   (kg per sq m per s).
     &,SNOW(POINTS)                                                     &
     &,SNOW2(POINTS)                                                    &
                       ! Dummy in 3B version of the code
!        On input: Rate of snowfall entering this layer from above.
!        On Output: Rate of snowfall leaving this layer.
!                    (kg per sq m per s).
     &,VF(POINTS)                                                       &
     &,VF2(POINTS)                                                      &
                    ! Dummy in 3B version of the code
!        On input: Fall speed of ice into layer from above.
!        On output: Fall speed of ice into layer below.
!                   (m per s).
     &,VF_RAIN(POINTS)                                                  &
                          ! Dummy in 3B version of the code
     &,GRAUPRATE(POINTS)                                                &
                          ! Dummy in 3B version of the code
     &,VF_GRAUP(POINTS)                                                 &
                          ! Dummy in 3B version of the code
     &,droplet_flux(points)                                             &
                               ! Dummy in 3B version of the code
     &,RHCRIT(POINTS)                                                   &
!        Critical humidity for cloud formation.
     &,SO4_ACC(POINTS)                                                  &
                               ! Dummy in 3B version of the code
     &,SO4_DIS(POINTS)                                                  &
                               ! Dummy in 3B version of the code
     &,SO4_AIT(POINTS)                                                  &
                               ! Dummy in 3B version of the code
     &,BMASS_AGD(POINTS)                                                &
                               ! Dummy in 3B version of the code
     &,BMASS_CLD(POINTS)                                                &
                               ! Dummy in 3B version of the code
     &,OCFF_AGD(POINTS)                                                 &
                               ! Dummy in 3B version of the code
     &,OCFF_CLD(POINTS)                                                 &
                               ! Dummy in 3B version of the code
     &,SNOW_DEPTH(POINTS)                                               &
                               ! Dummy in 3B version of the code
     &,LAND_FRACT(POINTS)                                               &
                               ! Dummy in 3B version of the code
     &,AEROSOL(POINTS)                                                  &
                               ! Dummy in 3B version of the code
     &,FRAC_ICE_ABOVE(POINTS)                                           &
                               ! Dummy in 3B version of the code
     &,CTTEMP(POINTS)                                                   &
                               ! Dummy in 3B version of the code
     &,RAINFRAC(POINTS)                                                 &
                               ! 1 if raining at end of routine
!                              ! 0 if not raining
     &,CFKEEP(POINTS)                                                   &
                               ! Dummy in 3B version of the code
     &,SEA_SALT_FILM(SALT_DIM_ICE)                                      &
                                   ! Dummy in 3B version
     &,SEA_SALT_JET(SALT_DIM_ICE)                                        &
                                   ! Dummy in 3B version
     &,biogenic(biog_dim_ice)
                                   ! Dummy in 3B version
!
! Process rate diagnostics (Dummy arrays in 3B version of code)
      REAL, Intent(InOut) ::                                            &
     &  PSDEP(POINTS)                                                   &
                       ! Deposition of vapour to snow aggregates
     &, PSAUT(POINTS)                                                   &
                       ! Autoconversion of aggregates from crystals
     &, PSACW(POINTS)                                                   &
                       ! Accretion of liq. water by snow aggregates
     &, PSACR(POINTS)                                                   &
                       ! Collection of rain by snow aggregates
     &, PSACI(POINTS)                                                   &
                       ! Collection of ice crystals by aggregates
     &, PSMLT(POINTS)                                                   &
                       ! Melting of snow aggregates
     &, PSMLTEVP(POINTS)  ! Evaporation of melting aggregates
      REAL, Intent(InOut) ::                                            &
     &  PRAUT(POINTS)                                                   &
                       ! Autoconversion of cloud drops to rain
     &, PRACW(POINTS)                                                   &
                       ! Accretion of liq. water by rain
     &, PREVP(POINTS)  ! Evaporation of rain
      REAL, Intent(InOut) ::                                            &
     &  PGAUT(POINTS)                                                   &
                       ! Autoconversion of graupel from aggregates
     &, PGACW(POINTS)                                                   &
                       ! Accretion of liq. water by graupel
     &, PGACS(POINTS)                                                   &
                       ! Collection of snow aggregates by graupel
     &, PGMLT(POINTS)  ! Melting of graupel
      REAL, Intent(InOut) ::                                            &
     &  PIFRW(POINTS)                                                   &
                       ! Homogeneous freezing nucleation
     &, PIPRM(POINTS)                                                   &
                       ! Heterogeneous (primary) nucleation
     &, PIDEP(POINTS)                                                   &
                       ! Deposition of vapour to ice crystals
     &, PIACW(POINTS)                                                   &
                       ! Accretion of liq. water by ice crystals
     &, PIACR(POINTS)                                                   &
                       ! Collection of rain by ice crystals
     &, PIMLT(POINTS)                                                   &
                       ! Melting of ice crystals
     &, PIMLTEVP(POINTS)  ! Evaporation of melting ice crystals
      REAL, Intent(InOut) ::                                            &
     &  PIFALL(POINTS)                                                  &
                       ! Sedimentation of ice crystals
     &, PSFALL(POINTS)                                                  &
                       ! Sedimentation of aggregates
     &, PRFALL(POINTS)                                                  &
                       ! Sedimentation of rain
     &, PGFALL(POINTS) ! Sedimentation of graupel
      REAL, Intent(InOut) ::                                            &
     &  PLSET(POINTS)                                                   &
                         ! Droplet settling of liquid water
     &, PLEVPSET(POINTS) ! Evaporated settled droplets
!
      LOGICAL                                                           &
                      !, INTENT(IN)
     &  BLAND(POINTS)                                                   &
!         Land/sea mask
     &, L_MURK                                                          &
                               ! Dummy in 3B version of the code
     &, L_USE_SULPHATE_AUTOCONV                                         &
                               ! Dummy in 3B version of the code
     &, L_BIOMASS_CCN                                                   &
                               ! Dummy in 3B version of the code
     &, L_OCFF_CCN                                                      &
                               ! Dummy in 3B version of the code
     &, L_SEASALT_CCN                                                   &
                               ! Dummy in 3B version of the code
     &, L_biogenic_CCN                                                  &
                               ! Dummy in 3B version of the code
     &, L_AUTO_DEBIAS                                                   &
                               ! Dummy in 3B version of the code
     &, L_mcr_qcf2                                                      &
                               ! Dummy in 3B version of the code
     &, L_mcr_qrain                                                     &
                               ! Dummy in 3B version of the code
     &, L_mcr_qgraup                                                    &
                               ! Dummy in 3B version of the code
     &, L_psd                                                           &
                               ! Dummy in 3B version of the code
     &, L_it_melting                                                    &
                               ! Use iterative melting
     &, L_non_hydrostatic                                               &
                               ! Use non hydrostatic layer masses
     &, L_mixing_ratio                                                  &
                               ! Use mixing ratio formulation
     &, l_cry_agg_dep                                                   &
                               ! Dummy in 3B version of the code
     &, l_droplet_settle       
                               ! Dummy in 3B version of the code

      Logical                                                           &
                      !, Intent(IN)
     &  L_seq_mcr                                                       &
                               ! Dummy in 3B version
     &, L_autoc_3b                                                      &
                               ! Dummy in 3B version
     &, l_autolim_3b                                                    &
                               ! Dummy in 3B version
     &, L_autoconv_murk        ! Dummy in 3B version

      Real                                                              &
                      !, Intent(IN)
     &  ec_auto                                                         &
                               ! Collision coalescence efficiency
     &, N_drop_land                                                     &
                               ! Number of droplets over land / m-3
     &, N_drop_sea                                                      &
                               ! Number of droplets over sea / m-3
     &, N_drop_land_cr                                                  &
                               ! N_drop_land ^ (-1/3) / m
     &, N_drop_sea_cr                                                   &
                               ! N_drop_sea ^ (-1/3) / m
     &, Ntot_land                                                       &
                               ! Dummy in 3B version
     &, Ntot_sea                                                        &
                               ! Dummy in 3B version
     &, ai, bi, aic, bic
                               ! Dummy in 3B version

!  Workspace usage: 3 real arrays---------------------------------------
      REAL                                                              &
     &  QS(POINTS)                                                      &
!         Saturated sp humidity for (T,p) in layer
     &, QSL(POINTS)                                                     &
!         Saturated sp humidity for (T,p) in layer
!         wrt water at all temps
     &, SNOWT(POINTS)
!         Cumulative fall out of snow within iterations.
!
! external subprograms are called --------------------------------------
!
!  Local (derived) physical constants ----------------------------------
      REAL LCRCP,LFRCP,LSRCP,CONW
      PARAMETER(                                                        &
     &  LCRCP=LC/CP                                                     &
!         Latent heat of condensation / Cp (K).
     &, LFRCP=LF/CP                                                     &
!         Latent heat of fusion / Cp (K).
     &, LSRCP=LCRCP+LFRCP                                               &
!         Sum of the above (S for Sublimation).
     &, CONW=R/(EPSILON*LC)                                             &
!         Constant in wet bulb temperature calculation.
     &  )
!
! ----------------------------------------------------------------------
!  1   Define local scalars.
! ----------------------------------------------------------------------
      INTEGER                                                           &
     &  I                                                               &
! Loop counter (horizontal field index).
     &, J                                                               &
! Counter for the iterations
     &, K                                                               &
! Counter for the LSITER2 loop
     &, LSITER2
! Number of times the advection and melting sections are iterated.
!
!       Reals effectively expanded to workspace by the Cray (using
!       vector registers).
!
      REAL                                                              &
!       Real workspace.  At end of DO loop, contains :-
     &  RHO(POINTS)                                                     &
!         Density of air in the layer.
     &, RHOR(POINTS)                                                    &
!         1.0/RHO to speed up calculations.
     &, VTEMP                                                           &
!         Virtual temperature as at start of loop.
     &, TEMPC                                                           &
!         temperature degree C as at start of loop.
     &, t_rapid_melt                                                    &
                      ! Rapid melting temperature
     &, ESI                                                             &
!         saturation vapour pressure (wrt ice below zero)
     &, ESW                                                             &
!         saturation vapour pressure (wrt water at all T)
     &, DQI                                                             &
!         increment to/from ice/snow
     &, DQIL                                                            &
!         increment to/from cloud water
     &, DPR                                                             &
!         increment to/from rain
     &, CFICETEMP                                                       &
!         fraction of ice inferred for fall speed calculations.
     &, FQI                                                             &
!         fallspeed for ice
     &, DHI(POINTS)                                                     &
!         CFL limit
     &, DHIR(POINTS)                                                    &
!         1.0/DHI
     &, DHILSITERR                                                      &
!         1.0/(DHI*LSITER)
     &, FQIRQI                                                          &
!         saved flux of ice out of layer
     &, QUP                                                             &
!         updated ice for long timestep
     &, QCLNEW                                                          &
!         updated liquid cloud in implicit calculations
     &, TEMP7                                                           &
!         term in melting
     &, PR02                                                            &
!         term in evaporation of rain
     &, PR04                                                            &
!         square of pr02
     &, QC                                                              &
!         term in autoconversion of cloud to rain
     &, APLUSB                                                          &
!         denominator in deposition or evaporation of ice
     &, CORR                                                            &
!         density correction for fall speed
     &, ROCOR                                                           &
!         density correction for fall speed
     &, RHO1                                                            &
!         surface air density
     &, VR1                                                             &
!         Mean fall speed of rain
     &, VS1                                                             &
!         Mean fall speed of snow
     &, LAMR1                                                           &
!         Inverse lambda in rain exponential distribution
     &, LAMFAC1                                                         &
!         Expression containing calculations with lambda
     &, LAMS1                                                           &
!         lambda in snow exponential distribution
     &, FV1                                                             &
!         Mean velocity difference between rain and snow
     &, TIMESTEP                                                        &
!         Timestep of each iteration
     &, CORR2                                                           &
!         Temperature correction of viscosity etc.
     &, RHNUC                                                           &
!         Relative humidity required for nucleation
     &, TCG                                                             &
!         Temperature Factor for X1I in Cox Golding calculation
     &, TCGI                                                            &
!         Inverse of TCG
     &, RATEQ                                                           &
!         Constant effecting rate of deposition/sublimation of ice
     &, RATEQCF                                                         &
!         Constant representing effect of sub grid distribution of ice
     &, RATEQS(POINTS)                                                  &
!         Critical humidity for ice deposition
     &, RATEQSL                                                         &
!         Critical humidity for rain evaporation
     &, HM_NORMALIZE                                                    &
!         Normalization for Hallett Mossop process
     &, HM_RATE
!         Increase in deposition due to Hallett Mossop process
! Obtain the size for CX and CONSTP
#include "c_lspsiz.h"

      Real                                                              &
     &  autorate_land                                                   &
                               ! Autoconversion rate consts over land
     &, autorate_sea                                                    &
                               ! Autoconversion rate constants over sea
     &, autolim_land                                                    &
                               ! Autoconversn limit constants over land
     &, autolim_sea            ! Autoconversn limit constants over sea

! ----------------------------------------------------------------------
!  Calculate autoconversion parameters
! ----------------------------------------------------------------------
! The numbers 5907.24 and 4188.79 represent combinations of
! physical constants. Do NOT change them.
! Please refer to UMDP26

      Autorate_land = 5907.24*ec_auto*inhomog_rate*n_drop_land_cr
      Autorate_sea  = 5907.24*ec_auto*inhomog_rate*n_drop_sea_cr
      Autolim_land  = 4188.79*r_thresh**3*n_drop_land*inhomog_lim
      Autolim_sea   = 4188.79*r_thresh**3*n_drop_sea*inhomog_lim

! ----------------------------------------------------------------------
!  2.1 Start the microphysics calculations
! ----------------------------------------------------------------------
! Set up the iterations
       TIMESTEP=TIMESTEPFIXED/LSITER
! Set up sub grid scale constants. For no sub grid scale variability
! use the set up RATEQ=1.0, RATEQCF=0.0, RATEQS=1.0.
! Ideally represent these in a comdeck but require RHCRIT
       RATEQ=1.0
       RATEQCF=0.0
! Set up SNOWT(I) to be zero for all I.
! Points_do1:
       DO I=1,POINTS
         SNOWT(I)=0.0
         RATEQS(I)=RHCRIT(I)
       END DO ! Points_do1
!
! ----------------------------------------------------------------------
!  2.2 Start iterating.
! ----------------------------------------------------------------------
! Iters_do1:
       DO J=1,LSITER
! ----------------------------------------------------------------------
!  2.3  Calculate sat humidity mixing ratio(respect to ice below zero)
! ----------------------------------------------------------------------
! DEPENDS ON: qsat_mix
       call qsat_mix(qs,T,p,points,l_mixing_ratio)
! ----------------------------------------------------------------------
!  2.4 Calculate sat humidity mixing ratio(respect to water at all temps
! ----------------------------------------------------------------------
! DEPENDS ON: qsat_wat_mix
       call qsat_wat_mix(qsl,T,p,points,l_mixing_ratio)
! ----------------------------------------------------------------------
!  2.5 Start loop over points.
! ----------------------------------------------------------------------
! Points_do2:
       DO I=1,POINTS
!
! ----------------------------------------------------------------------
!  3.1 Calculate density of air, RHO, via virtual temperature.
! ----------------------------------------------------------------------
         If (l_non_hydrostatic) then
           If (l_mixing_ratio) then
             ! rho is the dry density
             rho(i) = rhodz_dry(i) / deltaz(i)
           Else
             ! rho is the moist density
             rho(i) = rhodz_moist(i) / deltaz(i)
           End if  ! l_mixing_ratio
         Else
           ! Use the pressure to retrieve the moist density
           ! through the equation of state
           VTEMP=T(I)*(1.+C_VIRTUAL*Q(I)-QCL(I)-QCF(I)) ! Virtual temp
           RHO(I)=P(I)/(R*VTEMP)
         End if

         RHOR(I)=1.0/RHO(I)
         RHO1=1.0
! Correction factor of fall speeds etc. due to density.

        ! Correction of fall speeds
        If (l_mixing_ratio .and. l_non_hydrostatic) then
          corr = (rho1*deltaz(i) / rhodz_moist(i))**0.4
        Else
          corr = (rho1*rhor(i))**0.4
        End if

! Correction factor in viscosity etc. due to temperature.
         CORR2=(T(I)/273.0)**1.5 * (393.0/(T(I)+120.0))

        ! Combined correction factor
        If (l_mixing_ratio .and. l_non_hydrostatic) then
          rocor = sqrt(rhodz_moist(i) / deltaz(i) * corr * corr2)
        Else
          rocor = sqrt(rho(i)*corr*corr2)
        End if

! Calculate RATEQS as a function of RHCRIT and CFICE(I)
         RATEQS(I)=RHCRIT(I)+CFICE(I)*(1.0-RHCRIT(I))
         RATEQSL=RHCRIT(I)+CFLIQ(I)*(1.0-RHCRIT(I))
!
! ----------------------------------------------------------------------
!  3.2 set T deg C sat vapour pressures in N/m2
! ----------------------------------------------------------------------
         TEMPC=T(I)-ZERODEGC
         ESI=QS(I)*P(I)/EPSILON
         ESW=QSL(I)*P(I)/EPSILON
! Calculate a temperature factor for N0snow. CX(13)=1.0 if there is a
! temperature dependence, and 0.0 if there is not.
         TCG=EXP(-CX(13)*TEMPC/8.18)
         TCGI=1.0/TCG
!
! ----------------------------------------------------------------------
!  4   Check that ice cloud fraction is sensible.
! ----------------------------------------------------------------------
         CFICETEMP=CFICE(I)
! The possibility exists in multiple iterations that ice cloud fraction
! is equal to zero but nucleation and deposition from the last iteration
! has produced a finite ice content. Hence this section produces a fix
! which will stop the scheme crashing. Only need to use for more than 1
! iteration
         IF (LSITER >  1) THEN
           IF (QCF(I) >  0.0 .AND. CFICE(I) <= 0.1) THEN
             CFICETEMP=MAX(CFLIQ(I),0.1)
           END IF
         END IF
!
! ----------------------------------------------------------------------
!  5.1 Calculate fallspeeds for ice and water and limit
!      to CFL criteria for time to fall through layer.
! ----------------------------------------------------------------------
! Estimate fall speed out of this layer. We want to avoid advecting
! very small amounts of snow between layers, as this can cause numerical
! problems in other routines, so if QCF is smaller than a single
! nucleation mass per metre cubed don't advect it.
         IF (QCF(I) >  M0) THEN
!
! Estimate the mean fall speed across the entire gridbox.
! Use a top hat distrubution within the gridbox
           FQI=CONSTP(3)*CORR*                                          &
     &         (RHO(I)*QCF(I)*CONSTP(1)*TCGI/CFICETEMP)**CX(3)
         ELSE
! QCF is smaller than zero so set fall speed to zero
           FQI=0.0
! Endif for calculation of fall speed
         END IF

         ! Calculate CFL quantity of timestep over level separation
         If (l_non_hydrostatic) then
           ! Use the formulation based on the heights of the levels
           dhi(i) = timestep/deltaz(i)
         Else
           ! Use the formulation based on the pressure difference
           ! across a layer.
           dhi(i)        = timestep * rho(i)/rhodz(i)
         End if

! Define DHIR and DHILSITERR to speed up calculations.
         DHIR(I)=1.0/DHI(I)
         DHILSITERR=1.0/(DHI(I)*LSITER)
!
! Calculate the additional iterations required by the advection
! and melting schemes.
!
         IF (QCF(I) >  M0.AND.T(I) >  ZERODEGC.AND.                     &
     &       T(I) <  (ZERODEGC+2.0) .and. l_it_melting) THEN
           LSITER2 = VF(I)*DHI(I) + 1
           ! Limit the number of iterations
           IF (LSITER2  >   5) THEN
             LSITER2 = 5
           ENDIF
         ELSE
           ! Do not iterate
           LSITER2 = 1
         ENDIF
!
! Now start the additional iterations loop
         DO K=1,LSITER2
!
!
! See if fall speed is small enough that not all the ice falls out of
! the layer.
!
           IF(VF(I) <= (DHIR(I)*LSITER2))THEN
! short timestep solution
           VF(I)=FQI
             IF(VF(I) <= (DHIR(I)*LSITER2))THEN
! flux out is just represented by the fall speed estimated above
             FQIRQI=FQI*RHO(I)*QCF(I)
           ELSE
! cannot allow more ice to leave than was already there
               FQIRQI=RHO(I)*QCF(I)*DHIR(I)*LSITER2
           ENDIF
! calculate new ice content in this layer by flux divergence
             QCF(I)=QCF(I)+(SNOW(I)-FQIRQI)*DHI(I)/LSITER2*RHOR(I)
         ELSE
! long timestep case
           QUP = SNOW(I)*RHOR(I)/VF(I)
             FQIRQI = SNOW(I)-(RHODZ(I)*(QUP-QCF(I))/TIMESTEP*LSITER2)
           QCF(I) = QUP
! VF must be as least as great as the fall velocity of the current layer
           IF(VF(I) <  FQI) VF(I)=FQI
!
!
         END IF
! ----------------------------------------------------------------------
!  5.2 Snow is used to save fall out of layer
!      for calculation of fall into next layer
! ----------------------------------------------------------------------
         SNOWT(I)=SNOWT(I)+FQIRQI/(LSITER*LSITER2)
!
         If (L_it_melting) Then
           ! Move melting section to after the fall out section
! ----------------------------------------------------------------------
!  11  Melting of snow - explicit
!      USE WET BULB TEMP (DEG.C) IN SNOW MELT CALC.
!      Use a numerical approximation.
! ----------------------------------------------------------------------

           IF(QCF(I) >  M0.AND.T(I) >  ZERODEGC)THEN
             TEMPC=T(I)-ZERODEGC
! An approximate calculation of wet bulb temperature
             TEMP7=TEMPC-RATEQ*                                         &
     &            MAX((RATEQS(I)*QSL(I)+RATEQCF*QCF(I)-Q(I)),0.0)       &
     &           *(TW1+TW2*(P(I)-TW3) - TW4*(T(I)-TW5) )
             TEMP7=MAX(TEMP7,0.0)
! End of wet bulb temp formulations.
             PR02=RHO(I)*QCF(I)*CONSTP(1)*TCGI
             DPR=TCG*CONSTP(7)*TIMESTEP/LSITER2*                        &
     &              (0.65*CONSTP(13)*CORR2*PR02**CX(1)                  &
     &           + CONSTP(6)*ROCOR*PR02**CX(2))*RHOR(I)
! Solve implicitly in terms of temperature
             DPR=TEMP7*(1.0-1.0/(1.0+DPR*LFRCP))/LFRCP
             DPR=MIN(DPR,QCF(I))
! Update values of ice and Rain
            QCF(I)=QCF(I)-DPR
             RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR
             T(I)=T(I)-LFRCP*DPR
! ENDIF for melting snow
           END IF
!
         End If ! Iterative melting
!
!
!
! End of loop over additional iterations
         ENDDO
! ----------------------------------------------------------------------
!        Transfer processes only active at T less than 0 deg C
! ----------------------------------------------------------------------
         IF(T(I) <  ZERODEGC) THEN
! ----------------------------------------------------------------------
!  6   Nucleation of snow - explicit
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
!  6.1 Homogenous nucleation takes place at temperatures less than THOMO
! ----------------------------------------------------------------------
            IF (T(I) <  (ZERODEGC+THOMO)) THEN
! Turn all liquid to ice
              QCF(I)=QCF(I)+QCL(I)
              T(I)=T(I)+LFRCP*QCL(I)
              QCL(I)=0.0
            END IF
! ----------------------------------------------------------------------
!  6.2 Heteorgenous nucleation occurs for temps less than TNUC deg C
! ----------------------------------------------------------------------
           IF (T(I) <  (ZERODEGC+TNUC)) THEN
! Calculate number of active ice nucleii
             DQI=MIN(0.01*EXP(-0.6*TEMPC),1.0E5)
! Each nucleus can grow to arbitary mass of M0 kg
             DQI=M0*DQI*RHOR(I)
! RHNUC represents how much moisture is available for ice formation.
             RHNUC=(188.92+2.81*(T(I)-ZERODEGC)                         &
     &       +0.013336*(T(I)-ZERODEGC)**2-10.0)*0.01
             RHNUC=MIN(RHNUC,1.0)
! Predict transfer of mass to ice.
             DQI=MAX(MIN(DQI,Q(I)+QCL(I)                                &
     &           -RATEQS(I)*MAX(QSL(I)*RHNUC,QS(I))),0.0)
             QCF(I)=QCF(I)+DQI
! This comes initially from liquid water
             DQIL=MIN(DQI,QCL(I))
             QCL(I)=QCL(I)-DQIL
             T(I)=T(I)+LFRCP*DQIL
! If more moisture is required then nucleation removes from vapour.
             DQI=DQI-DQIL
             T(I)=T(I)+LSRCP*DQI
             Q(I)=Q(I)-DQI
! END IF for nucleation
           END IF
! ----------------------------------------------------------------------
!  7   Deposition/Sublimation of snow - explicit
! ----------------------------------------------------------------------
           IF(QCF(I) >  M0) THEN
! Calculate transfer rate as a function of QCF and T
             PR02=RHO(I)*QCF(I)*CONSTP(1)*TCGI
             APLUSB=(APB1-APB2*T(I))*ESI
             APLUSB=APLUSB+(T(I)**3)*P(I)*APB3
             DQI=TCG*CONSTP(5)*T(I)**2*ESI*RATEQ*                       &
     &         (MIN((Q(I)+QCL(I)),QSL(I))-RATEQCF*QCF(I)-               &
     &          RATEQS(I)*QS(I))*                                       &
     &       (0.65*CONSTP(13)*CORR2*PR02**CX(1)+CONSTP(6)*ROCOR*        &
     &       PR02**CX(2))/(QS(I)*APLUSB*RHO(I))
! Limits depend on whether deposition or sublimation occurs
             IF (DQI >  0.0) THEN
! Deposition is occuring. Limit is available moisture.
               DQI=MIN(DQI*TIMESTEP,                                    &
     &         (Q(I)+QCL(I)-RATEQCF*QCF(I)-RATEQS(I)*QS(I))/            &
     &         (1.0+RATEQCF))
             ELSE
! Sublimation is occuring. Limits are spare moisture capacity and QCF
               DQI=MAX(MAX(DQI*TIMESTEP,                                &
     &           (Q(I)+QCL(I)-RATEQCF*QCF(I)-RATEQS(I)*QS(I))           &
     &                            /(1.0+RATEQCF)),-QCF(I))
             END IF
! Adjust ice content
             QCF(I)=QCF(I)+DQI
             DQIL=MAX(MIN(DQI,QCL(I)),0.0)
! Adjust liquid content (deposits before vapour by Bergeron Findeison
!  process).
             QCL(I)=QCL(I)-DQIL
             T(I)=T(I)+LFRCP*DQIL
             DQI=DQI-DQIL
! Adjust vapour content
             Q(I)=Q(I)-DQI
             T(I)=T(I)+LSRCP*DQI
! END IF for QCF >  M0.
           END IF
! ----------------------------------------------------------------------
!  8   Riming of snow by cloud water -implicit in QCL
! ----------------------------------------------------------------------
           IF (QCF(I) >  M0.AND.QCL(I) >  0.0) THEN
               QCLNEW=QCL(I)/(1.0+CONSTP(4)*TCG*CORR*TIMESTEP*          &
     &         (RHO(I)*QCF(I)*CONSTP(1)*TCGI)**CX(4))
! Recalculate water contents
               QCF(I)=QCF(I)+(QCL(I)-QCLNEW)
               T(I)=T(I)+LFRCP*(QCL(I)-QCLNEW)
               QCL(I)=QCLNEW
! END IF for QCF >  M0.AND.QCL(I) >  0.0
           END IF
! ----------------------------------------------------------------------
!  9   Capture of rain by snow - implicit in rain
!      This is a candidate to be removed to speed things up.
! ----------------------------------------------------------------------
           IF (RAIN(I) >  0.0.AND.QCF(I) >  M0) THEN
! Calculate velocities
             VR1=CORR*CONSTP(11)/6.0*(RAIN(I)/(CONSTP(8)*CORR))**CX(5)
             VS1=CONSTP(3)*CORR*(RHO(I)*QCF(I)*CONSTP(1)*TCGI)**CX(3)
! Estimate the mean absolute differences in velocities.
             FV1=MAX(ABS(VR1-VS1),(VR1+VS1)/8.0)
! Calculate functions of slope parameter lambda
! MOD fdm2f404 changes next 3 active lines
! LAMR1 and LAMS1 are INVERSE lambda values
             LAMR1=(RAIN(I)/(CONSTP(8)*CORR))**(CX(10))
             LAMS1=(RHO(I)*QCF(I)*CONSTP(1)*TCGI)**(-CX(6))
             LAMFAC1=CONSTP(16)*(LAMR1**6.0*LAMS1**CX(16)) +            &
     &               CONSTP(15)*(LAMR1**5.0*LAMS1**CX(15)) +            &
     &               CONSTP(14)*(LAMR1**4.0*LAMS1**CX(14))
! Calculate transfer
!fdm2f404 change      DPR=TCG*CONSTP(9)*LAMS1**CX(8)*LAMR1**CX(9)*FV1*
             DPR=TCG*CONSTP(9)*LAMS1**(-CX(8))*LAMR1**(-CX(9))*FV1*     &
     &       LAMFAC1*TIMESTEP*RHOR(I)
             DPR=MIN(DPR,RAIN(I)*(DHI(I)*LSITER)*RHOR(I))
! Adjust ice and rain contents
             QCF(I)=QCF(I)+DPR
             RAIN(I)=RAIN(I)-DPR*RHO(I)*DHILSITERR
             T(I)=T(I)+LFRCP*DPR
!      Endif for RAIN >  0.0 in capture term
           END IF
! ----------------------------------------------------------------------
!      End of transfer processes only active at T less than 0 deg C
! ----------------------------------------------------------------------
         END IF
! ----------------------------------------------------------------------
!  10  Evaporate melting snow - implicit in subsaturation
! ----------------------------------------------------------------------
         IF(QCF(I) >  M0.AND.T(I) >  ZERODEGC)THEN
! Calculate transfer as a function of QCF, T and specific humidity
           PR02=RHO(I)*QCF(I)*CONSTP(1)*TCGI
           PR04=((APB4-APB5*T(I))*ESW+APB6*P(I)*T(I)**3)
           DPR=TCG*RATEQ*CONSTP(5)*T(I)**2*ESW*TIMESTEP*                &
     &     (0.65*CONSTP(13)*CORR2*PR02**CX(1)                           &
     &      +CONSTP(6)*ROCOR*PR02**CX(2))/(QSL(I)*RHO(I)*PR04)
           DPR=DPR*MAX((RATEQS(I)*QSL(I)+                               &
     &                  RATEQCF*QCF(I)-Q(I)-QCL(I)),0.0)                &
     &                /(1.0+DPR*(1.0+RATEQCF))
! Extra check to see we don't get a negative QCF
           DPR=MIN(DPR,QCF(I))
! Update values of ice and vapour
           QCF(I)=QCF(I)-DPR
           Q(I)=Q(I)+DPR
           T(I)=T(I)-DPR*LSRCP
         END IF
!
! If iterative melting is active then we have already done this term.
!
         If (.not. L_it_melting) Then
! ----------------------------------------------------------------------
!  11  Melting of snow - explicit
!      USE WET BULB TEMP (DEG.C) IN SNOW MELT CALC.
!      Use a numerical approximation.
! ----------------------------------------------------------------------
         IF(QCF(I) >  M0.AND.T(I) >  ZERODEGC)THEN
           TEMPC=T(I)-ZERODEGC
! An approximate calculation of wet bulb temperature
           TEMP7=TEMPC-RATEQ*MAX((RATEQS(I)*QSL(I)+                     &
     &             RATEQCF*QCF(I)-Q(I)),0.0)                            &
     &           *(TW1+TW2*(P(I)-TW3) - TW4*(T(I)-TW5) )
           TEMP7=MAX(TEMP7,0.0)
! End of wet bulb temp formulations.
           PR02=RHO(I)*QCF(I)*CONSTP(1)*TCGI
           DPR=TCG*CONSTP(7)*TIMESTEP*(0.65*CONSTP(13)*CORR2*PR02**CX(1)&
     &         + CONSTP(6)*ROCOR*PR02**CX(2))*RHOR(I)
! Solve implicitly in terms of temperature
           DPR=TEMP7*(1.0-1.0/(1.0+DPR*LFRCP))/LFRCP
           DPR=MIN(DPR,QCF(I))
! Update values of ice and Rain
           QCF(I)=QCF(I)-DPR
           RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR
           T(I)=T(I)-LFRCP*DPR
! ENDIF for melting snow
         END IF
!
         End If  ! .not. L_it_melting
!
! ----------------------------------------------------------------------
!  12  Evaporation of rain - implicit in subsaturation
! ----------------------------------------------------------------------
         IF(RAIN(I) >  0.0)THEN
           PR04=((APB4-APB5*T(I))*ESW+APB6*P(I)*T(I)**3)
! New, consistent evaporation method, with rain fall speed relationship.
!fdm2f404 mod    LAMR1=(CONSTP(8)*CORR/RAIN(I))**(CX(10))
! LAMR1 is now an INVERSE lambda value
           LAMR1=(RAIN(I)/(CONSTP(8)*CORR))**(CX(10))

           DPR=CONSTP(2)*T(I)**2*ESW*TIMESTEP
!fdm2f404 mod           DPR=DPR*( (0.78*CORR2/LAMR1**CX(12))
!     &               + (CONSTP(12)*ROCOR/(LAMR1**CX(11))) )
           DPR=DPR*( (0.78*CORR2*LAMR1**CX(12))                         &
     &               + (CONSTP(12)*ROCOR*(LAMR1**CX(11))) )

           DPR=DPR*RATEQ/(QSL(I)*RHO(I)*PR04)
! Calculate transfers.
           DPR=DPR*MAX((RATEQSL*QSL(I)-Q(I)-QCL(I)),0.0)/(1.0+DPR)
           DPR=DPR*RHO(I)*DHILSITERR
           DPR=MIN(DPR,RAIN(I))
! Update values of rain and vapour
           RAIN(I)=RAIN(I)-DPR
           Q(I)=Q(I)+DPR*DHI(I)*LSITER*RHOR(I)
           T(I)=T(I)-DPR*LCRCP*DHI(I)*LSITER*RHOR(I)
! END IF for evaporation of rain.
         END IF
! ----------------------------------------------------------------------
!  13  Accretion of cloud on rain - implicit in liquid water content
! ----------------------------------------------------------------------
           IF(RAIN(I) >  0.0.AND.QCL(I) >  0.0)THEN
! New accretion formulation.
             PR02=RAIN(I)/(CONSTP(8)*CORR)
             QCLNEW=QCL(I)/((1.0+CONSTP(10)*CORR*TIMESTEP*PR02**CX(7)))
! Now calculate increments to rain.
             RAIN(I)=RAIN(I)+(QCL(I)-QCLNEW)*RHO(I)*DHILSITERR
             QCL(I)=QCL(I)-(QCL(I)-QCLNEW)
! END IF for accretion of cloud on rain.
           END IF
! ----------------------------------------------------------------------
!  14  Autoconversion of cloud to rain - explicit
! ----------------------------------------------------------------------
           IF (QCL(I) >  0.0.AND.CFLIQ(I) >  0.0) THEN
! Use a liquid cloud fraction here as this term is very non-linear
! The section below is a simple way of proceeding.
             IF (BLAND(I)) THEN
! Land point
               QC=MIN(AUTOLIM_LAND*CFLIQ(I)*RHOR(I),QCL(I))
               DPR=MIN(AUTORATE_LAND*(RHO(I)*QCL(I)/CFLIQ(I))**1.333    &
     &                 *TIMESTEP*QCL(I)/CORR2,QCL(I)-QC)
             ELSE
! Sea point
               QC=MIN(AUTOLIM_SEA*CFLIQ(I)*RHOR(I),QCL(I))
               DPR=MIN(AUTORATE_SEA*(RHO(I)*QCL(I)/CFLIQ(I))**1.333     &
     &                 *TIMESTEP*QCL(I)/CORR2,QCL(I)-QC)
             END IF

! End of calculation of autoconversion amount DPR
             QCL(I)=QCL(I)-DPR
             RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR
! ENDIF for autoconversion.
            END IF
! Update rain fraction for diagnostic
         IF (RAIN(I) >  0.0) THEN
           RAINFRAC(I)=1.0
         ENDIF
! ----------------------------------------------------------------------
!  15  Now continue the loops over points and iterations.
! ----------------------------------------------------------------------
! Continue DO loop over points
      END DO ! Points_do2
! Continue DO loop over iterations
      END DO ! Iters_do1
!
! Copy contents of SNOWT to SNOW, to fall into next layer down
! Points_do3
      DO I=1,POINTS
        SNOW(I)=SNOWT(I)
! ----------------------------------------------------------------------
!  16 Melt any SNOW which has reached here, as long as T is large enough
! ----------------------------------------------------------------------
! Only use if long timestep case. In which case melt the excess snow
! which falls straight through a layer. Use DQI variable to save space.
! DQI APPROXIMATELY represents the excess SNOW.
           DQI=MIN(VF(I)*DHI(I)-1.0,1.0)
           IF (DQI >  0.0) THEN
! Long timestep case
             TEMPC=T(I)-ZERODEGC
             ! Define rapid melting temperature
             If (l_it_melting) then
               t_rapid_melt = zerodegc + 2.0
             Else
               t_rapid_melt = zerodegc
             End If
             If (snow(i)  >   0.0 .and. T(i) >   t_rapid_melt) then
! Numerical approximation of wet bulb temperature.
               TEMP7=TEMPC-RATEQ*MAX((RATEQS(I)*QSL(I)                  &
     &           +RATEQCF*QCF(I)-Q(I)),0.0)*(TW1+TW2*                   &
     &           (P(I)-TW3) - TW4*(T(I)-TW5) )
               TEMP7=MAX(TEMP7,0.0)
! End of wet bulb calculation
               DPR=TEMP7/(LFRCP*LSITER)
               DPR=MIN(DPR,SNOW(I)*DHI(I)*RHOR(I)*DQI)
! Update values of snow and rain
               SNOW(I)=SNOW(I)-DPR*RHO(I)*DHIR(I)
               RAIN(I)=RAIN(I)+DPR*RHO(I)*DHIR(I)
               T(I)=T(I)-LFRCP*DPR*LSITER
! END IF for long timestep
             END IF
! END IF for melting of excess snow.
           END IF
! ----------------------------------------------------------------------
!  17  Remove any small amount of QCF which is left over to be tidy.
!      If QCF is less than 1e-10 kg/kg (hardwired) and isn't growing
!      by deposition (assumed to be given by RHCRIT) then remove it.
! ----------------------------------------------------------------------
!           DQI=M0*RHOR(I)*MIN( 0.01*EXP(-0.6*TEMPC),1.0E5 )
!           IF (QCF(I) <  MIN( MAX(M0*RHOR(I),DQI),1.0E-5*QS(I) ) ) THEN
           IF (QCF(I) <  1.0E-10.AND.                                   &
     &     (T(I) >  ZERODEGC .OR. (Q(I)+QCL(I)  <=  RHCRIT(I)*QS(I))    &
     &     .OR. QCF(I) <  0.0)  )  THEN
             Q(I)=Q(I)+QCF(I)
             T(I)=T(I)-LSRCP*QCF(I)
             QCF(I)=0.0
           END IF
! END DO for melting of excess snow loop over points.
      END DO ! Points_do3
! ----------------------------------------------------------------------
!  18  End of the LSP_ICE subroutine
! ----------------------------------------------------------------------
      RETURN
      END SUBROUTINE LSP_ICE
#endif
