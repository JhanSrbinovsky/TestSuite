#if defined(A04_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Precipitation microphysics calculations.
! Subroutine Interface:
      SUBROUTINE LSP_ICE(                                               &
     &  P,RHODZ, deltaz, rhodz_dry, rhodz_moist,                        &
     &  TIMESTEPFIXED,n_iterations,POINTS,L_MURK,                       &
     &  RHCPT,                                                          &
     &  L_USE_SULPHATE_AUTOCONV, L_AUTO_DEBIAS,                         &
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
     &  RAIN, VF_RAIN, SNOW_AGG, VF,                                    &
     &  SNOW_CRY, VF_CRY, GRAUPRATE, VF_GRAUP, droplet_flux,            &
     &  FRAC_ICE_ABOVE,CTTEMP,RAINFRAC,                                 &
     &  T,CFKEEP,CFLIQ,CFICEKEEP,BLAND,CX,CONSTP                        &
     &, PSDEP,PSAUT,PSACW,PSACR,PSACI,PSMLT,PSMLTEVP                    &
     &, PRAUT,PRACW,PREVP                                               &
     &, PGAUT,PGACW,PGACS,PGMLT                                         &
     &, PIFRW,PIPRM,PIDEP,PIACW,PIACR,PIMLT,PIMLTEVP                    &
     &, PIFALL,PSFALL,PRFALL,PGFALL,PLSET,PLEVPSET                      &
     & )

      IMPLICIT NONE
!
! Description:
!   Updates ice, liquid and vapour contents and temperature as a
!   result of microphysical processes.
!
! Method:
!   Calculates transfers of water between vapour, ice, cloud liquid
!   and rain. Advects ice downwards. Processes included are:
!   Fall of ice into and out of the layer;
!   Homogenous and heterogenous nucleation of ice;
!   Deposition and sublimation of ice;
!   Riming; riming by supercooled raindrops;
!   Melting of ice; Evaporation of rain; accretion;
!   Autoconversion of liquid water to rain.
!   This is described in Unified Model Documentation Paper 26.
!
!   Microphysics options:
!   - Second prognostic cloud ice variables
!      Active if L_mcr_qcf2=.True.
!      The code supports the use of a second cloud ice prognostic
!      variable so that both cloud ice aggregates (QCF and QCF_AGG)
!      and cloud ice pristine crystals (QCF2 and QCF_CRY) can be
!      represented and advected.
!      At UM5.5 this code is still experimental.
!   - Prognostic rain (controlled by L_mcr_qrain)
!     At UM5.5 no microphysical process/sedimentation code is present.
!   - Prognostic graupel (controlled by L_mcr_qrain)
!     At UM5.5 no microphysical process/sedimentation code is present.
!
! Current Code Owner: Jonathan Wilkinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables:
#include "c_0_dg_c.h"
#include "c_lspmic.h"
#include "c_micro.h"
#include "c_lspdif.h"
! c_lspdif will call c_r_cp, c_epslon and c_lheat
#include "c_pi.h"
!
! Subroutine arguments
!
! Obtain the size for CX and CONSTP
#include "c_lspsiz.h"
!
      LOGICAL                                                           &
                      !, INTENT(IN)
     &  L_MURK                                                          &
!         Murk aerosol is a valid quantity
     & ,L_USE_SULPHATE_AUTOCONV                                         &
!         Switch to use sulphate aerosol in the calculation of cloud
!         droplet number concentration in the autoconversion section,
!         i.e. activate the second indirect ("lifetime") effect.
     & ,L_SEASALT_CCN                                                   &
!         Switch to supplement sulphate aerosol with sea-salt if the
!         second indirect effect is active.
     & ,L_BIOMASS_CCN                                                   &
!         Switch to supplement sulphate aerosol with biomass smoke
!         aerosol if the second indirect effect is active.
     & ,L_OCFF_CCN                                                      &
!         Switch to supplement sulphate aerosol with fossil-fuel organic
!         carbon aerosol if the second indirect effect is active.
     & ,L_biogenic_CCN                                                  &
!         Switch to supplement sulphate aerosol with biogenic aerosol
!         if the second indirect effect is active.
     & ,L_AUTO_DEBIAS                                                   &
!         Switch to apply de-biasing correction to autoconversion rate.
     &, L_mcr_qcf2                                                      &
                       ! true if using 2nd cloud ice prognostic
     &, L_mcr_qrain                                                     &
                       ! true if using prognostic rain
     &, L_mcr_qgraup                                                    &
                       ! true if using prognostic graupel
     &, L_psd                                                           &
                       ! Dummy in 3C version of large-scale precip
     &, L_it_melting                                                    &
                       ! true if using iterative melting
     &, L_non_hydrostatic                                               &
                               ! Use non hydrostatic layer masses
     &, L_mixing_ratio                                                  &
                       ! Use mixing ratio formulation
     &, l_cry_agg_dep                                                   &
                       ! Dummy in 3C version of the code
     &, l_droplet_settle       
                       ! Dummy in 3C version of the code

      INTEGER                                                           &
                      !, INTENT(IN)
     &  POINTS                                                          &
!         Number of points to be processed.
     &, n_iterations                                                    &
                           ! Dummy in 3C version of the code
     & ,SALT_DIM_ICE                                                    &
!        Number of points for sea-salt arrays (either POINTS or 1)
     & ,biog_dim_ice
!        Number of points for biogenic array (either POINTS or 1)
!
      REAL                                                              &
                      !, INTENT(IN)
     &  TIMESTEPFIXED                                                   &
!         Timestep of physics in model (s).
     &, CFLIQ(POINTS)                                                   &
!         Liquid cloud fraction in this layer (no units).
     &, CFICEKEEP(POINTS)                                               &
!         Frozen cloud fraction in this layer (no units).
     &, CFKEEP(POINTS)                                                  &
!         Total cloud fraction in this layer (no units).
     &, P(POINTS)                                                       &
!         Air pressure at this level (Pa).
     &, RHODZ(POINTS)                                                   &
!         Air mass p.u.a. in this layer (kg per sq m).
     &, deltaz(points)                                                  &
                           ! Depth of layer / m
     &, rhodz_dry(points)                                               &
                           ! Density of dry air / kg m-3
     &, rhodz_moist(points)                                             &
                           ! Density of moist air / kg m-3
     &, SO4_ACC(POINTS)                                                 &
!         Sulphur cycle variable
     &, SO4_DIS(POINTS)                                                 &
!         Sulphur cycle variable
     &, SO4_AIT(POINTS)                                                 &
!         Sulphur cycle variable
     &, BMASS_AGD(POINTS)                                               &
!         Aged biomass smoke mass mixing ratio
     &, BMASS_CLD(POINTS)                                               &
!         In-cloud biomass smoke mass mixing ratio
     &, OCFF_AGD(POINTS)                                                &
!         Aged fossil-fuel organic carbon mass mixing ratio
     &, OCFF_CLD(POINTS)                                                &
!         In-cloud fossil-fuel organic carbon mass mixing ratio
     &, SNOW_DEPTH(POINTS)                                              &
!         Snow depth for aerosol amount (m)
     &, LAND_FRACT(POINTS)                                              &
!         Land fraction
     &, AEROSOL(POINTS)                                                 &
!         Aerosol mass (ug/kg)
     &, RHCPT(POINTS)                                                   &
!         Critical relative humidity of all points for cloud formation.
     &, SEA_SALT_FILM(SALT_DIM_ICE)                                     &
!         Film-mode sea-salt aerosol number concentration (m-3)
     &, SEA_SALT_JET(SALT_DIM_ICE)                                      &
!         Jet-mode sea-salt aerosol number concentration (m-3)
     &, biogenic(biog_dim_ice)
!         Biogenic aerosol m.m.r.
!
      REAL                                                              &
                      !, INTENT(INOUT)
     &  Q(POINTS)                                                       &
!         Specific humidity at this level (kg water per kg air).
     &, QCF(POINTS)                                                     &
!         Cloud ice (kg water per kg air).
     &, QCL(POINTS)                                                     &
!         Cloud liquid water (kg water per kg air).
     &, QCF2(POINTS)                                                    &
!         Second cloud ice (kg water per kg air).
     &, QRAIN(POINTS)                                                   &
!         Rain (kg water per kg air).
     &, QGRAUP(POINTS)                                                  &
!         Graupel (kg water per kg air).
     &, T(POINTS)                                                       &
!         Temperature at this level (K).
     &, RAIN(POINTS)                                                    &
!         On input: Rate of rainfall entering this layer from above.
!         On output: Rate of rainfall leaving this layer.
!                   (kg m-2 s-1).
     &, SNOW_AGG(POINTS)                                                &
!         On input: Rate of snow aggregates entering layer from above.
!         On output: Rate of snow aggregates leaving this layer.
!                   (kg m-2 s-1). If only one ice prognostic is active
!                   then this variable contains all the snow.
     &, SNOW_CRY(POINTS)                                                &
!         On input: Rate of snow crystals entering layer from above.
!         On Output: Rate of snow crystals leaving this layer.
!                   (kg m-2 s-1). Only non-zero if two ice
!                   prognostics in use.
     &, VF(POINTS)                                                      &
!         On input: Fall velocity of snow aggregates entering layer.
!         On Output: Fall velocity of snow aggregates leaving layer.
!                   (m s-1). If only one ice prognostic is active
!                   then this is the velocity of all falling snow.
     &, VF_CRY(POINTS)                                                  &
!         On input: Fall velocity of snow crystals entering this layer.
!         On Output: Fall velocity of snow crystals leaving this layer.
!                   (m s-1). Only used if two ice prognostics in use.
     &, VF_RAIN(POINTS)                                                 &
                           ! Dummy in 3C version of the code
     &, GRAUPRATE(POINTS)                                               &
                           ! Dummy in 3C version of the code
     &, VF_GRAUP(POINTS)                                                &
                           ! Dummy in 3C version of the code
     &, droplet_flux(points)                                            &
                             ! Dummy in 3C version of the code
     &, CTTEMP(POINTS)                                                  &
!         Ice cloud top temperature (K)
     &, RAINFRAC(POINTS)                                                &
!         Rain fraction (no units)
     &, FRAC_ICE_ABOVE(POINTS)
!         Fraction of ice in layer above (no units)
!
! Process rate diagnostics (Dummy arrays in 3C version of code)
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
     &  BLAND(POINTS)
!         Land/sea mask
!
      Logical                                                           &
                      !, Intent(IN)
     &  L_seq_mcr                                                       &
                               ! Dummy in 3C version
     &, L_autoc_3b                                                      &
                               ! Dummy in 3C version
     &, l_autolim_3b                                                    &
                               ! Dummy in 3C version
     &, L_autoconv_murk        ! Use murk aerosol to calc. drop number

      Real                                                              &
                      !, Intent(IN)
     &  ec_auto                                                         &
                               ! Collision coalescence efficiency
     &, N_drop_land                                                     &
                               ! Dummy in 3C version
     &, N_drop_sea                                                      &
                               ! Dummy in 3C version
     &, N_drop_land_cr                                                  &
                               ! Dummy in 3C version
     &, N_drop_sea_cr                                                   &
                               ! Dummy in 3C version
     &, Ntot_land                                                       &
                               ! Number of droplets over land / m-3
     &, Ntot_sea                                                        &
                               ! Number of droplets over sea / m-3
     &, ai, bi, aic, bic
                               ! Dummy in 3C version

!  Local parameters
      REAL LCRCP,LFRCP,LSRCP,RHO1,ONE_OVER_EPSILON,ONE_OVER_ZERODEGC
      PARAMETER(                                                        &
     &  LCRCP=LC/CP                                                     &
!         Latent heat of condensation / cp (K).
     &, LFRCP=LF/CP                                                     &
!         Latent heat of fusion / cp (K).
     &, LSRCP=LCRCP+LFRCP                                               &
!         Sum of the above (S for Sublimation).
     &, RHO1=1.0                                                        &
!         Reference density of air (kg m-3)
     &, ONE_OVER_EPSILON=1.0/EPSILON                                    &
!        Inverse of epsilon to speed up calculations
     &, ONE_OVER_ZERODEGC=1.0/ZERODEGC                                  &
!        Inverse of zero degrees Celsius to speed up code (K-1)
     &  )
!
!  Local scalars and dynamic arrays
!
      INTEGER                                                           &
     &  I                                                               &
!         Loop counter (horizontal field index).
     &, J                                                               &
!         Counter for the iterations
     &, K                                                               &
!         Counter for the LSITER2 loop
     &, LSITER2                                                         &
!         Number of times advection and melting sections are iterated
     &, KK                                                              &
!         Variable for condensed points compression
     &, KK2                                                             &
!         Another condensed points compression variable
     &, SEA_SALT_PTR                                                    &
!         Pointer for sea-salt arrays
     &, BIOMASS_PTR                                                     &
!         Pointer for biomass smoke arrays
     &, OCFF_PTR                                                        &
!         Pointer for fossil-fuel organic carbon arrays
     &, biogenic_ptr
!         Pointer for biogenic aerosol arrays
!
      REAL                                                              &
     &  QS(POINTS)                                                      &
!         Saturated sp humidity for (T,p) in layer (kg kg-1)
     &, QSL(POINTS)                                                     &
!         Saturated sp humidity for (T,p) in layer
!         wrt water at all temps (kg kg-1)
     &, QS_CTT(POINTS)                                                  &
!         Saturated sp humidity for (CTTEMP,p) in layer (kg kg-1)
     &, QSL_CTT(POINTS)                                                 &
!         Saturated sp humidity for (CTTEMP,p) in layer (kg kg-1)
!         wrt water at all temps
     &, SNOWT_AGG(POINTS)                                               &
!         Cumulative fall out of snow aggregates within iterations.
!         (kg m-2 s-1). If only one ice prognostic is active then
!         this variable contains all the snow.
     &, SNOWT_CRY(POINTS)                                               &
!         Cumulative fall out of snow crystals within iterations.
!         (kg m-2 s-1). Only non-zero if two ice prognostics in use.
     &, NUMBER_DROPLET                                                  &
!         Droplet concentration calculated using a function (m-3)
     &,  RHO(POINTS)                                                    &
!         Density of air in the layer (kg m-3).
     &, RHOR(POINTS)                                                    &
!         1.0/RHO to speed up calculations (kg-1 m3).
     &, VTEMP                                                           &
!         Virtual temperature as at start of loop (K).
     &, ESI(POINTS)                                                     &
!         saturation vapour pressure (wrt ice below zero Celsius)(Pa)
     &, ESW(POINTS)                                                     &
!         saturation vapour pressure (wrt water at all temperatures)(Pa)
     &, ESI_CTT(POINTS)                                                 &
!         saturation vapour pressure at cloud top temp (wrt ice) (Pa)
     &, ESW_CTT(POINTS)                                                 &
!         saturation vapour pressure at cloud top temp (wrt water) (Pa)
     &, DQI                                                             &
!         increment to/from ice/snow (kg kg-1)
     &, DQIL                                                            &
!         increment to/from cloud water (kg kg-1)
     &, DPR                                                             &
!         increment to/from rain (kg m-2 s-1)
     &, CFICE(POINTS)                                                   &
!         fraction of ice inferred for the microphysics (no units).
     &, CFICEI(POINTS)                                                  &
!         inverse of CFICE (no units)
     &, CF(POINTS)                                                      &
!         total cloud fraction for the microphysics (no units)
     &, FQI_AGG(POINTS)                                                 &
!         fallspeed for aggregates (m s-1)
     &, FQI_CRY(POINTS)                                                 &
!         fallspeed for aggregates (m s-1)
     &, DHI(POINTS)                                                     &
!         CFL limit (s m-1)
     &, DHIR(POINTS)                                                    &
!         1.0/DHI (m s-1)
     &, DHILSITERR(POINTS)                                              &
!         1.0/(DHI*LSITER) (m s-1)
     &, FQIRQI_AGG                                                      &
!         saved flux of ice out of layer (kg m-2 s-1)
     &, FQIRQI2_AGG(POINTS)                                             &
!         saved flux of ice out of layer from layer above (kg m-2 s-1)
     &, FQIRQI_CRY                                                      &
!         saved flux of ice out of layer (kg m-2 s-1)
     &, FQIRQI2_CRY(POINTS)                                             &
!         saved flux of ice out of layer from layer above (kg m-2 s-1)
     &, QCLNEW                                                          &
!         updated liquid cloud in implicit calculations (kg kg-1)
     &, TEMP7                                                           &
!         temporary variable
     &, TEMPW                                                           &
!         temporary for vapour calculations
     &, t_rapid_melt                                                    &
                      ! Rapid melting temperature
     &, PR02                                                            &
!         term in evaporation of rain
     &, PR04                                                            &
!         square of pr02
     &, QC                                                              &
!         term in autoconversion of cloud to rain (kg kg-1)
     &, APLUSB                                                          &
!         denominator in deposition or evaporation of ice
     &, CORR(POINTS)                                                    &
!         density correction for fall speed (no units)
     &, ROCOR(POINTS)                                                   &
!         density correction for fall speed (no units)
     &, VR1                                                             &
!         Mean fall speed of rain (m s-1)
     &, VS1                                                             &
!         Mean fall speed of snow (m s-1)
     &, LAMR1                                                           &
!         Inverse lambda in rain exponential distribution (m)
     &, LAMR2                                                           &
!         Inverse lambda in rain exponential distribution (m)
     &, LAMFAC1                                                         &
!         Expression containing calculations with lambda
     &, LAMS1                                                           &
!         Inverse lambda in snow exponential distribution (m)
     &, FV1                                                             &
!         Mean velocity difference between rain and snow (m s-1)
     &, TIMESTEP                                                        &
!         Timestep of each iteration (s)
     &, CORR2(POINTS)                                                   &
!         Temperature correction of viscosity etc. (no units)
     &, RHNUC                                                           &
!         Relative humidity required for nucleation (no units)
     &, TCG(POINTS)                                                     &
!         Temperature Factor for aggregate size distribution (no units)
     &, TCGI(POINTS)                                                    &
!         Inverse of TCG (no units)
     &, TCGC(POINTS)                                                    &
!         Temperature Factor for crystal size distribution (no units)
     &, TCGCI(POINTS)                                                   &
!         Inverse of TCGC (no units)
     &, RATEQS(POINTS)                                                  &
!         Sub grid model variable (no units)
     &, HM_NORMALIZE                                                    &
!         Normalization for Hallett Mossop process (no units)
     &, HM_RATE                                                         &
!         Increase in deposition due to Hallett Mossop process(no units)
     &, AREA_LIQ(POINTS)                                                &
!         Liquid only area of gridbox (no units)
     &, AREA_MIX(POINTS)                                                &
!         Mixed phase area of gridbox (no units)
     &, AREA_ICE(POINTS)                                                &
!         Ice only area of gridbox (no units)
     &, AREA_CLEAR(POINTS)                                              &
!         Cloud free area of gridbox (no units)
     &, RAIN_LIQ(POINTS)                                                &
!         Overlap fraction of gridbox between rain and liquid cloud
     &, RAIN_MIX(POINTS)                                                &
!         Overlap fraction of gridbox between rain and mixed phase cloud
     &, RAIN_ICE(POINTS)                                                &
!         Overlap fraction of gridbox between rain and ice cloud
     &, RAIN_CLEAR(POINTS)                                              &
!         Overlap fraction of gridbox between rain and no cloud
     &, Q_ICE(POINTS)                                                   &
!         Vapour content in the ice only part of the grid box (kg kg-1)
     &, Q_CLEAR(POINTS)                                                 &
!         Vapour content in the cloud free part of the grid box(kg kg-1)
     &, QCF_AGG(POINTS)                                                 &
!         QCF in the form of aggregates (kg kg-1)
     &, QCF_CRY(POINTS)                                                 &
!         QCF in the form of crystals (kg kg-1)
     &, FRAC_AGG(POINTS)                                                &
!         Fraction of aggregates (no units)
     &, TEMP1(POINTS)                                                   &
     &, TEMP2(POINTS)                                                   &
     &, TEMP3(POINTS)                                                   &
     &, TEMP4(POINTS)                                                   &
!         Temporary arrays for T3E vector functions
     &, POWER                                                           &
!         Power for T3E vector functions
     &, N_CCN(POINTS)                                                   &
!         Cloud condensation nuclei concentration (m-3)
     &, N_DROP(POINTS)                                                  &
!         Droplet concentration (m-3)
     &, A_FACTOR(POINTS)                                                &
     &, B_FACTOR(POINTS)                                                &
!         Numerical factors in autoconversion calculation
     &, R_MEAN(POINTS)                                                  &
!         Mean droplet radius (m)
     &, R_MEAN0                                                         &
!         Constant in the calculation of R_MEAN
     &, N_GT_20                                                         &
!         Concentration of particles greater than a threshold
!         radius in size (m-3)
     &, AUTOLIM                                                         &
!         Minimum water content for autoconversion (kg m-3)
     &, AUTORATE                                                        &
!         Rate constant for autoconversion
     &, AC_FACTOR                                                       &
!         Autoconversion de-biasing factor
     &, T_L(POINTS), QSL_TL(POINTS)                                     &
     &, ALPHA_L, A_L, SIGMA_S, G_L, GACB                                &
!         Subsidiary terms used in the autoconversion de-biasing
     &, QCFAUTOLIM                                                      &
!         Minimum ice content for autoconversion (kg/kg)
     &, QCFAUTORATE                                                     &
!         Rate constant for ice autoconversion
     &, QCF_TOT(POINTS)                                                 &
!         Total amount of ice (crystals+aggregates)
     &, LHEAT_CORREC_LIQ(POINTS)                                        &
!         Reduction factor in evaporation limits because of latent heat
     &, LHEAT_CORREC_ICE(POINTS)                                        &
!         Reduction factor in evaporation limits because of latent heat
     &, TEMPR                                                           &
!         Temporary for optimization
     &, DQI_DEP, DQI_SUB, TEMPW_DEP, TEMPW_SUB, WIDTH(POINTS)           &
     &, Q_ICE_1(POINTS), Q_ICE_2(POINTS)                                &
     &, AREA_ICE_1(POINTS), AREA_ICE_2(POINTS)
!         Subgrid splitting for deposition term
!
!  Function and Subroutine calls
      EXTERNAL NUMBER_DROPLET
!
!- End of header
!
! ----------------------------------------------------------------------
!  2.1 Set up some variables
! ----------------------------------------------------------------------
! Set up the iterations
       TIMESTEP=TIMESTEPFIXED/LSITER

       DO I=1,POINTS
         ! Set up SNOWT to be zero for all I
         SNOWT_CRY(I) = 0.0
         SNOWT_AGG(I) = 0.0
       END DO

! Set qcf for crystals and aggregates depending on
! whether there are one or two ice prognostics

      IF (L_mcr_qcf2) THEN ! two ice prognostics

        DO I=1,POINTS
          QCF_CRY(I) = QCF2(I)
          QCF_AGG(I) = QCF(I)
        END DO

      ELSE ! only one ice prognostic, split diagnostically

       DO I=1,POINTS ! Points 0
! Work out fraction of ice in aggregates
#if defined(VECTLIB)
         FRAC_AGG(I)=-T_SCALING*MAX((T(I)-CTTEMP(I)),0.0)               &
     &               *MAX(QCF(I)*QCF0,0.0)
#else
         FRAC_AGG(I)=MAX(1.0-EXP(-T_SCALING*MAX((T(I)-CTTEMP(I)),0.0)   &
     &               *MAX(QCF(I),0.0)*QCF0) , 0.0)
#endif
       END DO ! Points 0
#if defined(VECTLIB)
! DEPENDS ON: exp_v
       CALL EXP_V(POINTS,FRAC_AGG,FRAC_AGG)
#endif
! Points_do1:
       DO I=1,POINTS
#if defined(VECTLIB)
         FRAC_AGG(I)=MAX(1.0-FRAC_AGG(I) , 0.0)
#endif
! Allocate ice content to crystals and aggregates
         QCF_CRY(I)=QCF(I)*(1.0-FRAC_AGG(I))
         QCF_AGG(I)=QCF(I)*FRAC_AGG(I)
         ! Assume falling snow is partitioned into crystals and
         ! aggregates, SNOW_AGG contains total snow on input
         SNOW_CRY(I)=SNOW_AGG(I)*(1.0-FRAC_AGG(I))
         SNOW_AGG(I)=SNOW_AGG(I)*FRAC_AGG(I)

! The compiler should unroll the above loop so break it here
       END DO
      ENDIF ! on L_mcr_qcf2 (number of ice prognostics)

      DO I=1,POINTS
        QCF_TOT(I) = QCF_CRY(I) + QCF_AGG(I)
      END DO

! Set up Hallett Mossop calculation
       HM_NORMALIZE=1.0/(1.0-EXP((HM_T_MIN-HM_T_MAX)*HM_DECAY))
! Set up autoconversion factor
       R_MEAN0=(27.0/(80.0*PI*1000.0))**(-1.0/3.0)
! ----------------------------------------------------------------------
!  2.2 Start iterating.
! ----------------------------------------------------------------------
! Iters_do1:
       DO J=1,LSITER
! ----------------------------------------------------------------------
!  2.3  Calculate saturation specific humidities
! ----------------------------------------------------------------------
! Qsat with respect to ice
! DEPENDS ON: qsat_mix
       call qsat_mix(qs,T,p,points,l_mixing_ratio)
!
! Qsat with respect to liquid water
! DEPENDS ON: qsat_wat_mix
       call qsat_wat_mix(qsl,T,p,points,l_mixing_ratio)
! ----------------------------------------------------------------------
!  2.4 Start loop over points.
! ----------------------------------------------------------------------
! Points_do2:
       DO I=1,POINTS
!
! ----------------------------------------------------------------------
!  3.1 Calculate density of air and density dependent correction factors
! ----------------------------------------------------------------------
         ESI(I)=QS(I)*P(I)*ONE_OVER_EPSILON
         ESW(I)=QSL(I)*P(I)*ONE_OVER_EPSILON

        !-----------------------------------------------
        ! Calculate density of air
        !-----------------------------------------------
        If (l_non_hydrostatic) then
          If (l_mixing_ratio) then
            ! rho is the dry density
            rho(i) = rhodz_dry(i) / deltaz(i)
          Else
            ! rho is the moist density
            rho(i) = rhodz_moist(i) / deltaz(i)
          End if  ! l_mixing_ratio
        Else
          ! Use the pressure to retrieve the moist density.
          ! An exact expression for the air density is
          ! rho = p / (1+c_virtual q - qcl - qcf) / RT but
          ! the approximation below is more numerically
          ! stable.
          rho(i) = p(i)*(1.-c_virtual*q(i)+qcl(i)+qcf_tot(i))           &
     &           /(R*T(i))
        End if  ! l_non_hydrostatic

         RHOR(I)=1.0/RHO(I)
! Estimate latent heat correction to rate of evaporation etc.
         LHEAT_CORREC_LIQ(I)=1.0/(1.0+EPSILON*LC**2*QSL(I)              &
     &                            /(CP*R*T(I)**2))
         LHEAT_CORREC_ICE(I)=1.0/(1.0+EPSILON*(LC+LF)**2*QS(I)          &
     &                            /(CP*R*T(I)**2))
       END DO
       DO I=1,POINTS
#if defined(VECTLIB)
! Correction factor of fall speeds etc. due to density.
         If (l_mixing_ratio .and. l_non_hydrostatic) then
           corr(i) = (rho1*deltaz(i) / rhodz_moist(i))
         Else
           corr(i) = (rho1*rhor(i))
         End if
! Correction factor in viscosity etc. due to temperature.
         CORR2(I)=(T(I)*ONE_OVER_ZERODEGC)
#else
! Correction factor of fall speeds etc. due to density.
         If (l_mixing_ratio .and. l_non_hydrostatic) then
           corr(i) = (rho1*deltaz(i) / rhodz_moist(i))**0.4
         Else
           corr(i) = (rho1*rhor(i))**0.4
         End if
! Correction factor in viscosity etc. due to temperature.
         CORR2(I)=(T(I)*ONE_OVER_ZERODEGC)**1.5 * (393.0/(T(I)+120.0))
#endif
       ENDDO
#if defined(VECTLIB)
! Use T3E vector functions to speed up code
       CALL POWR_V(POINTS,CORR,0.4,CORR)
       CALL POWR_V(POINTS,CORR2,1.5,CORR2)
#endif
       DO I=1,POINTS
! ----------------------------------------------------------------------
!  3.2 Calculate particle size distributions. Use vector functions
!       if available because they are much faster
! ----------------------------------------------------------------------
#if defined(VECTLIB)
! Complete calculation of CORR2
         CORR2(I)=CORR2(I)*(393.0/(T(I)+120.0))
! Calculate a temperature factor for N0crystals.
         TCGC(I)=-CX(12)*MAX(T(I)-ZERODEGC,T_AGG_MIN)
! Calculate a temperature factor for N0aggregates.
         TCG(I)=-CX(32)*MAX(T(I)-ZERODEGC,T_AGG_MIN)
! Define inverse of TCG values to speed up calculations
         TCGCI(I)=-TCGC(I)
         TCGI(I)=-TCG(I)
! Crystal temperature dependence is the same
! Break loop as above section can be pipelined
       END DO
       DO I=1,POINTS
! Combined correction factor
         If (l_mixing_ratio .and. l_non_hydrostatic) then
           rocor(i) = rhodz_moist(i) / deltaz(i) * corr(i) * corr2(i)
         Else
           rocor(i) = rho(i)*corr(i)*corr2(i)
         End if
#else
! Combined correction factor
         If (l_mixing_ratio .and. l_non_hydrostatic) then
           rocor(i) = sqrt(rhodz_moist(i)/deltaz(i)*corr(i)*corr2(i))
         Else
           rocor(i) = sqrt(rho(i)*corr(i)*corr2(i))
         End if
! Calculate a temperature factor for N0aggregates. CX(32)>0.0
! if there is a temperature dependence, and 0.0 if there is not.
         TCG(I)=EXP(-CX(32)*MAX(T(I)-ZERODEGC,T_AGG_MIN))
! Temperature dependence for N0crystals
         TCGC(I)=EXP(-CX(12)*MAX(T(I)-ZERODEGC,T_AGG_MIN))
! Define inverse of TCG values to speed up calculations
         TCGI(I)=1.0/TCG(I)
         TCGCI(I)=1.0/TCGC(I)
#endif
       ENDDO
#if defined(VECTLIB)
! Use T3E vector functions to speed up code
! DEPENDS ON: exp_v
         CALL EXP_V(POINTS,TCG,TCG)
! DEPENDS ON: exp_v
         CALL EXP_V(POINTS,TCGC,TCGC)
! DEPENDS ON: exp_v
         CALL EXP_V(POINTS,TCGI,TCGI)
! DEPENDS ON: exp_v
         CALL EXP_V(POINTS,TCGCI,TCGCI)
         CALL POWR_V(POINTS,ROCOR,0.5,ROCOR)
#endif
       DO I=1,POINTS
!
! ----------------------------------------------------------------------
!  4.1 Check that ice cloud fraction is sensible.
! ----------------------------------------------------------------------
         CFICE(I)=MAX(MAX(CFICEKEEP(I),0.01),FRAC_ICE_ABOVE(I))
         CFICEI(I)=1.0/CFICE(I)
         CF(I)=CFKEEP(I)
         CF(I)=MAX(CF(I),CFICE(I))
! Break loop to aid efficient pipelining of calculations
       END DO
       DO I=1,POINTS
         RATEQS(I)=RHCPT(I)   ! Sub grid parameter
         FRAC_ICE_ABOVE(I)=CFICEKEEP(I)
!
! ----------------------------------------------------------------------
!  4.2 Calculate overlaps of liquid, ice and rain fractions
! ----------------------------------------------------------------------
         AREA_LIQ(I)=CF(I)-CFICE(I)
         AREA_MIX(I)=CFICE(I)+CFLIQ(I)-CF(I)
         AREA_ICE(I)=CF(I)-CFLIQ(I)
         AREA_CLEAR(I)=1.0-CF(I)
! Break loop again since the above can be pipelined and the next loop
! the compiler should unroll
       END DO
       DO I=1,POINTS
         RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
         RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
         RAIN_ICE(I)=                                                   &
     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
         RAIN_CLEAR(I)=RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
       END DO
       DO I=1,POINTS
! ----------------------------------------------------------------------
!  4.3 Calculate vapour contents in ice only and clear regions
! ----------------------------------------------------------------------
         IF (CFLIQ(I)  <   1.0) THEN
! TEMPW is the mean vapour content in the ice only and clear regions
           TEMPW=(Q(I) - CFLIQ(I)*QSL(I)) / (1.0 - CFLIQ(I))
           TEMP7 = ICE_WIDTH * QSL(I)
! 0.001 is to avoid divide by zero problems
           WIDTH(I) = 2.0 *(1.0-RATEQS(I))*QSL(I)                       &
     &                    *MAX(  (1.0-0.5*QCF(I)/TEMP7),0.001  )
! The full width cannot be greater than 2q because otherwise part of
! gridbox would have negative q.
! However, we must also ensure that a q of zero cannot produce an error.
           WIDTH(I) = MIN(WIDTH(I) ,    MAX(2.0*Q(I),0.001*QS(I) )    )
           IF (AREA_ICE(I)  >   0.0) THEN
             Q_CLEAR(I)= TEMPW - 0.5*WIDTH(I)*AREA_ICE(I)
              Q_ICE(I)=(Q(I)-CFLIQ(I)*QSL(I)-AREA_CLEAR(I)*Q_CLEAR(I))  &
     &                / AREA_ICE(I)
           ELSE
             Q_CLEAR(I) = TEMPW
             Q_ICE(I)=0.0
           END IF
         ELSE
! Q_CLEAR and Q_ICE are undefined. Set them to zero
           WIDTH(I)=1.0
           Q_CLEAR(I)=0.0
           Q_ICE(I)=0.0
         ENDIF
! ----------------------------------------------------------------------
!  4.4 Update ice cloud top temperature if no ice falling in
! ----------------------------------------------------------------------
         IF (SNOW_CRY(I)+SNOW_AGG(I)  <=  0.0) THEN
           CTTEMP(I)=T(I)
         ENDIF
! ----------------------------------------------------------------------
!  4.5  Remove any small amount of ice to be tidy.
!      If QCF is less than QCFMIN and isn't growing
!      by deposition (assumed to be given by RHCPT) then remove it.
! ----------------------------------------------------------------------
         IF (QCF_TOT(I)  <   QCFMIN) THEN
           IF (T(I) >  ZERODEGC .OR.                                    &
     &       (Q_ICE(I)  <=  QS(I) .AND. AREA_MIX(I)  <=  0.0)           &
     &       .OR. QCF_TOT(I)  <   0.0)  THEN
             Q(I)=Q(I)+QCF_TOT(I)
             T(I)=T(I)-LSRCP*QCF_TOT(I)
             QCF(I)     = 0.0
             QCF_TOT(I) = 0.0
             QCF_AGG(I) = 0.0
             QCF_CRY(I) = 0.0
           END IF
         END IF
       END DO
! ----------------------------------------------------------------------
!  5.1   Falling ice is advected downwards
! ----------------------------------------------------------------------
! Estimate fall speed out of this layer. We want to avoid advecting
! very small amounts of snow between layers, as this can cause numerical
! problems in other routines, so if QCF is smaller than a single
! nucleation mass per metre cubed don't advect it.
#if defined(VECTLIB)
       KK=0
       KK2=0
       DO I=1,POINTS
         IF (QCF_AGG(I) >  M0) THEN
         KK=KK+1
!
! Estimate the mean fall speed across the entire gridbox.
! Use a top hat distrubution within the gridbox
           TEMP1(KK)=                                                   &
     &     (RHO(I)*QCF_AGG(I)*CONSTP(25)*TCGI(I)*CFICEI(I))
         ENDIF
         IF (QCF_CRY(I) >  M0) THEN
         KK2=KK2+1
!
! Estimate the mean fall speed across the entire gridbox.
! Use a top hat distrubution within the gridbox
           TEMP2(KK2)=                                                  &
     &     (RHO(I)*QCF_CRY(I)*CONSTP(5)*TCGCI(I)*CFICEI(I))
         ENDIF
       ENDDO
       CALL POWR_V(KK,TEMP1,CX(23),TEMP1)
       CALL POWR_V(KK2,TEMP2,CX(3),TEMP2)
! Set condensed points back to zero
       KK=0
       KK2=0
#endif
       DO I=1,POINTS
         IF (QCF_AGG(I) >  M0) THEN
#if defined(VECTLIB)
           KK=KK+1
           FQI_AGG(I)=CONSTP(24)*CORR(I)*TEMP1(KK)
#else
           FQI_AGG(I)=CONSTP(24)*CORR(I)*                               &
     &     (RHO(I)*QCF_AGG(I)*CONSTP(25)*TCGI(I)*CFICEI(I))**CX(23)
#endif
         ELSE
! QCF is smaller than zero so set fall speed to zero
           FQI_AGG(I)=0.0
! Endif for calculation of fall speed
         END IF
         IF (QCF_CRY(I) >  M0) THEN
#if defined(VECTLIB)
           KK2=KK2+1
           FQI_CRY(I)=CONSTP(4)*CORR(I)*TEMP2(KK2)
#else
           FQI_CRY(I)=CONSTP(4)*CORR(I)*                                &
     &       (RHO(I)*QCF_CRY(I)*CONSTP(5)*TCGCI(I)*CFICEI(I))**CX(3)
#endif
         ELSE
! QCF is smaller than zero so set fall speed to zero
           FQI_CRY(I)=0.0
! Endif for calculation of fall speed
         END IF
!
! ----------------------------------------------------------------------
!  5.2 Calculate CFL quantity of timestep over level separation.
! ----------------------------------------------------------------------
         If (l_non_hydrostatic) then
           ! Use the formulation based on the heights of the levels
           dhi(i) = timestep/deltaz(i)
         Else
           ! Use the formulation based on the pressure difference
           ! across a layer.
           dhi(i)        = timestep * rho(i)/rhodz(i)
         End if

! Define DHIR and DHILSITERR(I) to speed up calculations.
         DHIR(I)=1.0/DHI(I)
         DHILSITERR(I)=1.0/(DHI(I)*LSITER)
!
! ----------------------------------------------------------------------
!  5.2a Adjust fall speeds to make a linear combination of the fall
!       speed from the layer above. This ensures that the fall speed of
!       ice in this layer will not be calculated as zero even though
!       there is ice falling into it from above.
! ----------------------------------------------------------------------
        ! If using two cloud ice prognostics, calculate fallspeeds
        !  separately for crystals and aggregates, in which case
        !  VF represents aggregates, VF_CRY represents crystals
        ! Otherwise calculate an average fallspeed represented by VF

        IF (L_mcr_qcf2) THEN

         ! Crystals
         TEMP7 = SNOW_CRY(I)*DHI(I)*RHOR(I)
         IF (QCF_CRY(I)+TEMP7  >   M0) THEN
           TEMPR = 1.0/(QCF_CRY(I)+TEMP7)
           FQI_CRY(I)=(FQI_CRY(I)*QCF_CRY(I)+VF_CRY(I)*TEMP7)*TEMPR
           VF_CRY(I) = FQI_CRY(I)
         ELSE
           VF_CRY(I)=0.0
         ENDIF

         ! Aggregates
         TEMP7 = SNOW_AGG(I)*DHI(I)*RHOR(I)
         IF (QCF_AGG(I)+TEMP7  >   M0) THEN
           TEMPR = 1.0/(QCF_AGG(I)+TEMP7)
           FQI_AGG(I)=(FQI_AGG(I)*QCF_AGG(I)+VF(I)*TEMP7)*TEMPR
           VF(I) = FQI_AGG(I)
         ELSE
           VF(I)=0.0
         ENDIF

        ELSE ! L_mcr_qcf2=.F.-> use totals to estimate fallspeed

         TEMP7 = (SNOW_AGG(I) + SNOW_CRY(I)) *DHI(I)*RHOR(I)
         IF((QCF_CRY(I)+QCF_AGG(I)+TEMP7)  >   M0) THEN
           TEMPR=1.0 / (QCF_CRY(I)+QCF_AGG(I)+TEMP7)
           FQI_CRY(I)=(FQI_CRY(I)*(QCF_CRY(I)+QCF_AGG(I))+VF(I)*TEMP7)  &
     &                * TEMPR
           FQI_AGG(I)=(FQI_AGG(I)*(QCF_CRY(I)+QCF_AGG(I))+VF(I)*TEMP7)  &
     &                * TEMPR
           VF(I)=FQI_CRY(I)*(1.0-FRAC_AGG(I))+FQI_AGG(I)*FRAC_AGG(I)
         ELSE
           VF(I)=0.0
         ENDIF

        ENDIF ! on L_mcr_qcf2
!
! ----------------------------------------------------------------------
!  5.2b Calculate the additional iterations required by the advection
! and melting schemes.
! ----------------------------------------------------------------------
!
         IF ((QCF_CRY(I)+QCF_AGG(I)) >  M0.AND.T(I) >  ZERODEGC.AND.    &
     &       T(I) <  (ZERODEGC+2.0) .and. l_it_melting) THEN
           LSITER2 = MAX(FQI_CRY(I),FQI_AGG(I))*DHI(I) + 1
           IF (LSITER2  >   5) THEN
             LSITER2 = 5
           ENDIF
         ELSE
           LSITER2 = 1
         ENDIF
!
! Now start the additional iterations loop
         DO K=1,LSITER2
!
! ----------------------------------------------------------------------
!  5.3 Analytical Solution
! ----------------------------------------------------------------------
! FQIRQI2 is used as a temporary for the exponential function
! Calculate values point by point
         FQIRQI2_AGG(I)=EXP(-FQI_AGG(I)*DHI(I)/LSITER2)
         FQIRQI2_CRY(I)=EXP(-FQI_CRY(I)*DHI(I)/LSITER2)
! Advect aggregates
         IF (FQI_AGG(I)  >   0.0) THEN
           ! Assume falling snow is partitioned into crystals and
           ! aggregates
           FQIRQI_AGG = SNOW_AGG(I)+DHIR(I)*LSITER2*                    &
     &                  (RHO(I)*QCF_AGG(I)-SNOW_AGG(I)/FQI_AGG(I))      &
     &                  * (1.0-FQIRQI2_AGG(I))
           QCF_AGG(I) = SNOW_AGG(I)*RHOR(I)/FQI_AGG(I)                  &
     &                  * (1.0-FQIRQI2_AGG(I))                          &
     &                  + QCF_AGG(I)*FQIRQI2_AGG(I)
         ELSE  ! No fall of QCF out of the layer
           FQIRQI_AGG = 0.0
           QCF_AGG(I) = SNOW_AGG(I)*RHOR(I)*DHI(I)/LSITER2
         END IF

! Advect crystals
         IF (FQI_CRY(I)  >   0.0) THEN
           FQIRQI_CRY = SNOW_CRY(I)+DHIR(I)*LSITER2*                    &
     &                  (RHO(I)*QCF_CRY(I)-SNOW_CRY(I)/FQI_CRY(I))      &
     &                  * (1.0-FQIRQI2_CRY(I))
           QCF_CRY(I) = SNOW_CRY(I)*RHOR(I)/FQI_CRY(I)                  &
     &                  * (1.0-FQIRQI2_CRY(I))                          &
     &                  + QCF_CRY(I)*FQIRQI2_CRY(I)
         ELSE  ! No fall of QCF out of the layer
           FQIRQI_CRY=0.0
           QCF_CRY(I)=SNOW_CRY(I)*RHOR(I)*DHI(I)/LSITER2
         END IF
! No need to compute fall speed out of the layer in this method.
! ----------------------------------------------------------------------
!  5.4 Snow is used to save fall out of layer
!      for calculation of fall into next layer
! ----------------------------------------------------------------------
           SNOWT_CRY(I) = SNOWT_CRY(I) + (FQIRQI_CRY)/(LSITER*LSITER2)
           SNOWT_AGG(I) = SNOWT_AGG(I) + (FQIRQI_AGG)/(LSITER*LSITER2)
!
! ----------------------------------------------------------------------
!  11  Melting of snow. Uses a numerical approximation to the wet bulb
!      temperature.
! ----------------------------------------------------------------------
         If (L_it_melting) Then
         ! Do iterative melting

         IF(QCF_CRY(I) >  M0.AND.T(I) >  ZERODEGC)THEN
! ----------------------------------------------------------------------
!  11.1 Crystals first
! ----------------------------------------------------------------------
! An approximate calculation of wet bulb temperature
! TEMPW represents AVERAGE supersaturation in the ice
! Strictly speaking, we need to do two melting calculations,
! for the two different wet bulb temperatures.
           TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0)*CFICEI(I)
           TEMP7=T(I)-ZERODEGC-TEMPW                                    &
     &           *(TW1+TW2*(P(I)-TW3) - TW4*(T(I)-TW5) )
           TEMP7=MAX(TEMP7,0.0)
! End of wet bulb temp formulations.
           PR02=RHO(I)*MAX(QCF_CRY(I),0.0)*CFICEI(I)*CONSTP(5)*TCGCI(I)
           DPR=TCGC(I)*CONSTP(14)*TIMESTEP/LSITER2*                     &
     &            (CONSTP(7)*CORR2(I)*PR02**CX(4)                       &
     &         + CONSTP(8)*ROCOR(I)*PR02**CX(5))*RHOR(I)
! Solve implicitly in terms of temperature
           DPR=TEMP7*(1.0-1.0/(1.0+DPR*LFRCP))/LFRCP
           DPR=MIN(DPR,QCF_CRY(I))
! Update values of ice and Rain
           QCF_CRY(I)=QCF_CRY(I)-DPR
           RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
           T(I)=T(I)-LFRCP*DPR
           IF (DPR >  0.0) THEN
             RAINFRAC(I)=MAX(RAINFRAC(I),CFICE(I))
! Update rain fractions
             RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
             RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
             RAIN_ICE(I)=                                               &
     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
             RAIN_CLEAR(I)=                                             &
     &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
           ENDIF
! ENDIF for melting snow
         END IF
!
! ----------------------------------------------------------------------
!  11.1 Aggregates next
! ----------------------------------------------------------------------
         IF(QCF_AGG(I) >  M0.AND.T(I) >  ZERODEGC)THEN
! An approximate calculation of wet bulb temperature
! TEMPW represents AVERAGE supersaturation in the ice
! Strictly speaking, we need to do two melting calculations,
! for the two different wet bulb temperatures.
           TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0)*CFICEI(I)
           TEMP7=T(I)-ZERODEGC-TEMPW                                    &
     &           *(TW1+TW2*(P(I)-TW3) - TW4*(T(I)-TW5) )
           TEMP7=MAX(TEMP7,0.0)
! End of wet bulb temp formulations.
           PR02=RHO(I)*MAX(QCF_AGG(I),0.0)*CFICEI(I)*CONSTP(25)*TCGI(I)
           DPR=TCG(I)*CONSTP(34)*TIMESTEP/LSITER2*                      &
     &            (CONSTP(27)*CORR2(I)*PR02**CX(24)                     &
     &         + CONSTP(28)*ROCOR(I)*PR02**CX(25))*RHOR(I)
! Solve implicitly in terms of temperature
           DPR=TEMP7*(1.0-1.0/(1.0+DPR*LFRCP))/LFRCP
           DPR=MIN(DPR,QCF_AGG(I))
! Update values of ice and Rain
           QCF_AGG(I)=QCF_AGG(I)-DPR
           RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
           T(I)=T(I)-LFRCP*DPR
           IF (DPR >  0.0) THEN
             RAINFRAC(I)=MAX(RAINFRAC(I),CFICE(I))
! Update rain fractions
             RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
             RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
             RAIN_ICE(I)=                                               &
     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
             RAIN_CLEAR(I)=                                             &
     &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
           ENDIF
! ENDIF for melting snow
         END IF

         End If  ! l_it_melting
!
! End of loop over additional iterations
         ENDDO
!
! ----------------------------------------------------------------------
!  6.1 Homogenous nucleation takes place at temperatures less than THOMO
! ----------------------------------------------------------------------
            IF (T(I) <  (ZERODEGC+THOMO)) THEN
! Turn all liquid to ice
              QCF_CRY(I)=QCF_CRY(I)+QCL(I)
              T(I)=T(I)+LFRCP*QCL(I)
              QCL(I)=0.0
            END IF
! ----------------------------------------------------------------------
!  6.2 Heteorgenous nucleation occurs for temps less than TNUC deg C
! ----------------------------------------------------------------------
           IF (T(I) <  (ZERODEGC+TNUC).AND.                             &
     &       QCF_TOT(I) <  (1.0E5*M0*RHOR(I))) THEN
! Calculate number of active ice nucleii.
             DQI=MIN(0.01*EXP(-0.6*(T(I)-ZERODEGC)),1.0E5)
! Each nucleus can grow to arbitary mass of M0 kg
             DQI=M0*DQI*RHOR(I)
! DQI is amount of nucleation
             DQI=MAX(DQI-QCF_TOT(I),0.0)
             IF (DQI >  0.0) THEN
! RHNUC represents how much moisture is available for ice formation.
             RHNUC=(188.92+2.81*(T(I)-ZERODEGC)                         &
     &       +0.013336*(T(I)-ZERODEGC)**2)*0.01
             RHNUC=MIN(RHNUC,1.0)-0.1
! Predict transfer of mass to ice.
             RHNUC=MAX(QSL(I)*RHNUC,QS(I))
! DQIL is amount of moisture available
             DQIL=(QCL(I)+QSL(I)-RHNUC)*CFLIQ(I)                        &
     &            +MAX(Q_ICE(I)-RHNUC,0.0)*AREA_ICE(I)                  &
     &            +MAX(Q_CLEAR(I)-RHNUC,0.0)*AREA_CLEAR(I)
             DQI=MAX(MIN(DQI,DQIL*LHEAT_CORREC_ICE(I)),0.0)
             QCF_CRY(I)=QCF_CRY(I)+DQI
! This comes initially from liquid water
             DQIL=MIN(DQI,QCL(I))
             QCL(I)=QCL(I)-DQIL
             T(I)=T(I)+LFRCP*DQIL
! If more moisture is required then nucleation removes from vapour.
             DQI=DQI-DQIL
             T(I)=T(I)+LSRCP*DQI
             Q(I)=Q(I)-DQI
! END IFs for nucleation
             END IF
           END IF
!
! ----------------------------------------------------------------------
!  7   Deposition/Sublimation of snow.
!      Hallett Mossop process enhances growth.
! ----------------------------------------------------------------------
! Assume we can't calculate a meaningful width
           IF (Q_ICE(I)  >   QS(I)) THEN
! Deposition
             AREA_ICE_1(I) = AREA_ICE(I)
             AREA_ICE_2(I) = 0.0
             Q_ICE_1(I)    = Q_ICE(I)
             Q_ICE_2(I)    = QS(I)       ! Dummy value
           ELSE
! Sublimation
             AREA_ICE_1(I) = 0.0
             AREA_ICE_2(I) = AREA_ICE(I)
             Q_ICE_1(I)    = QS(I)       ! Dummy value
             Q_ICE_2(I)    = Q_ICE(I)
           END IF
           IF (AREA_ICE(I)  >   0.0) THEN
             TEMP7 = 0.5*AREA_ICE(I) + (Q_ICE(I)-QS(I)) / WIDTH(I)
! Temp7 is now the estimate of the proportion of the gridbox which
! contains ice and has local q greater than saturation wrt ice
             IF (TEMP7  >   0.0 .AND. TEMP7  <   AREA_ICE(I)) THEN
! Calculate values of q in each region.
! These overwrite previous estimates.
               AREA_ICE_1(I) = TEMP7
               AREA_ICE_2(I) = AREA_ICE(I) - AREA_ICE_1(I)
               Q_ICE_1(I)=QS(I) + 0.5 * AREA_ICE_1(I) * WIDTH(I)
               Q_ICE_2(I)=QS(I) - 0.5 * AREA_ICE_2(I) * WIDTH(I)
             END IF
           END IF
         ENDDO
!
#if defined(VECTLIB)
         KK=0
         DO I=1,POINTS
           IF (QCF_CRY(I) >  M0.AND.T(I) <  ZERODEGC) THEN
             KK=KK+1
! ----------------------------------------------------------------------
!  7.1 Crystals first
! ----------------------------------------------------------------------
! Calculate transfer rate as a function of QCF and T
             TEMP1(KK)=RHO(I)*QCF_CRY(I)*CFICEI(I)*CONSTP(5)*TCGCI(I)
           END IF
         ENDDO
! Take powers of function of ice content
         CALL POWR_V(KK,TEMP1,CX(5),TEMP2)
         CALL POWR_V(KK,TEMP1,CX(4),TEMP1)
         KK=0
#endif
         DO I=1,POINTS
           IF (QCF_CRY(I) >  M0.AND.T(I) <  ZERODEGC) THEN
             TEMP3(I)=QCF_CRY(I)+QCF_AGG(I)
! Diffusional parameters
             APLUSB=(APB1-APB2*T(I))*ESI(I)
             APLUSB=APLUSB+(T(I)**3)*P(I)*APB3
! Moisture available from subgrid scale calculation
             TEMPW_DEP = QSL(I)*AREA_MIX(I)                             &
     &                 + MIN(Q_ICE_1(I),QSL(I)) * AREA_ICE_1(I)         &
     &                 - QS(I) * (AREA_MIX(I) + AREA_ICE_1(I))
             TEMPW_SUB = (Q_ICE_2(I) - QS(I)) * AREA_ICE_2(I)
#if defined(VECTLIB)
             KK=KK+1
             LAMR1=TEMP1(KK)
             LAMR2=TEMP2(KK)
#else
             PR02=RHO(I)*MAX(QCF_CRY(I),0.0)                            &
     &                  *CFICEI(I)*CONSTP(5)*TCGCI(I)
             LAMR1=PR02**CX(4)
             LAMR2=PR02**CX(5)
#endif
! Transfer rate
             DQI=TCGC(I)*CONSTP(6)*T(I)**2*ESI(I)*                      &
     &       (CONSTP(7)*CORR2(I)*LAMR1+CONSTP(8)*ROCOR(I)*              &
     &       LAMR2)/(QS(I)*APLUSB*RHO(I))
             DQI_DEP = DQI * TEMPW_DEP
             DQI_SUB = DQI * TEMPW_SUB
! Limits depend on whether deposition or sublimation occurs
             IF (DQI_DEP >  0.0) THEN
! Deposition is occuring.
! Hallett Mossop Enhancement
               IF ( (T(I)-ZERODEGC)  >=  HM_T_MAX) THEN
! Temperature is greater than maximum threshold for HM.
                 HM_RATE=0.0
               ELSEIF ((T(I)-ZERODEGC)  <   HM_T_MAX                    &
! Temperature is between HM thresholds
     &           .AND. (T(I)-ZERODEGC)  >   HM_T_MIN) THEN
                 HM_RATE=(1.0-EXP( (T(I)-ZERODEGC-HM_T_MAX)*HM_DECAY) ) &
     &             *HM_NORMALIZE
               ELSE
! Temperature is less than minimum threshold for HM.
                 HM_RATE=EXP( (T(I)-ZERODEGC-HM_T_MIN)*HM_DECAY)
               ENDIF
! Calculate enhancement factor for HM process.
               HM_RATE=1.0+HM_RATE*QCL(I)*HM_RQCL
! The molecular diffusion to or from the surface is more efficient
! when a particle is at a molecular step. This is more likely when
! a particle is subliming. For growth, reduce the rate by 10 percent.
               HM_RATE=0.9*HM_RATE
               TEMPW_DEP=TEMPW_DEP*LHEAT_CORREC_ICE(I)
! Calculate Transfer. Limit is available moisture.
               IF (CFLIQ(I) >  0.0) THEN
! Add on liquid water contribution to available moisture
! The latent heat correction at this stage becomes very tedious to
! calculate. Freezing the liquid should raise qsat and hence allow
! for some evaporation of liquid itself. It is best to assume that
! no latent heat effect from the freezing should be employed.
                 TEMPW_DEP=TEMPW_DEP+QCL(I)*AREA_MIX(I)/CFLIQ(I)        &
     &                    +MAX((Q_ICE_1(I)-QSL(I)),0.0)*AREA_ICE_1(I)
               ENDIF
               DQI_DEP=MIN(DQI_DEP*TIMESTEP*HM_RATE,TEMPW_DEP)
             END IF
             IF (DQI_SUB <  0.0) THEN
! Sublimation is occuring. Limits are spare moisture capacity and QCF
! outside the liquid cloud
               DQI_SUB=MAX(MAX                                          &
     &                 (DQI_SUB*TIMESTEP,TEMPW_SUB*LHEAT_CORREC_ICE(I)) &
     &                 ,-(QCF_CRY(I) * AREA_ICE_2(I) * CFICEI(I) ))
             END IF
! Adjust ice content
             QCF_CRY(I)=QCF_CRY(I)+DQI_DEP+DQI_SUB
!
             IF (CFLIQ(I) >  0.0 .AND. AREA_MIX(I) >  0.0               &
     &            .AND. QCL(I) >  0.0) THEN
! Deposition removes some liquid water content.
!
! First estimate of the liquid water removed is explicit
               DQIL=MAX(MIN( DQI_DEP*AREA_MIX(I)                        &
     &                /(AREA_MIX(I)+AREA_ICE_1(I)),                     &
     &                QCL(I)*AREA_MIX(I)/CFLIQ(I)),0.0)
             ELSE
! Deposition does not remove any liquid water content
               DQIL=0.0
             ENDIF
! Adjust liquid content (deposits before vapour by Bergeron Findeison
!  process).
             QCL(I)=QCL(I)-DQIL
             T(I)=T(I)+LFRCP*DQIL
             DQI=DQI_DEP+DQI_SUB-DQIL
! Adjust vapour content
             Q(I)=Q(I)-DQI
             T(I)=T(I)+LSRCP*DQI
! END IF for QCF >  M0.
           END IF
!
         ENDDO
! ----------------------------------------------------------------------
!  7.2 Now aggregates
! ----------------------------------------------------------------------
#if defined(VECTLIB)
         KK=0
         DO I=1,POINTS
           IF (QCF_AGG(I) >  M0.AND.T(I) <  ZERODEGC) THEN
             KK=KK+1
! Calculate transfer rate as a function of QCF and T
             TEMP1(KK)=RHO(I)*QCF_AGG(I)*CFICEI(I)*CONSTP(25)*TCGI(I)
           END IF
         ENDDO
! Take powers of function of ice content
         CALL POWR_V(KK,TEMP1,CX(25),TEMP2)
         CALL POWR_V(KK,TEMP1,CX(24),TEMP1)
         KK=0
#endif
         DO I=1,POINTS
           IF (QCF_AGG(I) >  M0.AND.T(I) <  ZERODEGC) THEN
             TEMP3(I)=QCF_CRY(I)+QCF_AGG(I)
! Diffusional parameters
             APLUSB=(APB1-APB2*T(I))*ESI(I)
             APLUSB=APLUSB+(T(I)**3)*P(I)*APB3
! Moisture available from subgrid scale calculation
             TEMPW_DEP = QSL(I)*AREA_MIX(I)                             &
     &                 + MIN(Q_ICE_1(I),QSL(I)) * AREA_ICE_1(I)         &
     &                 - QS(I) * (AREA_MIX(I) + AREA_ICE_1(I))
             TEMPW_SUB = (Q_ICE_2(I) - QS(I)) * AREA_ICE_2(I)
#if defined(VECTLIB)
             KK=KK+1
             LAMR1=TEMP1(KK)
             LAMR2=TEMP2(KK)
#else
             PR02=RHO(I)*MAX(QCF_AGG(I),0.0)                            &
     &                  *CFICEI(I)*CONSTP(25)*TCGI(I)
             LAMR1=PR02**CX(24)
             LAMR2=PR02**CX(25)
#endif
! Transfer rate
             DQI=TCGC(I)*CONSTP(26)*T(I)**2*ESI(I)*                     &
     &       (CONSTP(27)*CORR2(I)*LAMR1+CONSTP(28)*ROCOR(I)*            &
     &       LAMR2)/(QS(I)*APLUSB*RHO(I))
             DQI_DEP = DQI * TEMPW_DEP
             DQI_SUB = DQI * TEMPW_SUB
! Limits depend on whether deposition or sublimation occurs
             IF (DQI_DEP >  0.0) THEN
! Deposition is occuring.
!
! The molecular diffusion to or from the surface is more efficient
! when a particle is at a molecular step. This is more likely when
! a particle is subliming. For growth, reduce the rate by 10 percent.
               DQI_DEP=0.9*DQI_DEP
               TEMPW_DEP=TEMPW_DEP*LHEAT_CORREC_ICE(I)
! Calculate Transfer. Limit is available moisture.
               IF (CFLIQ(I) >  0.0) THEN
! Add on liquid water contribution to available moisture
                 TEMPW_DEP=TEMPW_DEP+QCL(I)*AREA_MIX(I)/CFLIQ(I)        &
     &                    +MAX((Q_ICE_1(I)-QSL(I)),0.0)*AREA_ICE_1(I)
               ENDIF
               DQI_DEP=MIN(DQI_DEP*TIMESTEP,TEMPW_DEP)
             END IF
             IF (DQI_SUB  <   0.0) THEN
! Sublimation is occuring. Limits are spare moisture capacity and QCF
! outside liquid cloud
               DQI_SUB=MAX(MAX(                                         &
     &                 DQI_SUB*TIMESTEP,TEMPW_SUB*LHEAT_CORREC_ICE(I))  &
     &                    ,-(QCF_AGG(I) * AREA_ICE_2(I) * CFICEI(I) ))
             END IF
! Adjust ice content
             QCF_AGG(I)=QCF_AGG(I)+DQI_SUB+DQI_DEP
!
             IF (CFLIQ(I) >  0.0 .AND. AREA_MIX(I) >  0.0               &
     &            .AND. QCL(I) >  0.0) THEN
! Deposition removes some liquid water content.
!
! First estimate of the liquid water removed is explicit
               DQIL=MAX(MIN( DQI_DEP*AREA_MIX(I)                        &
     &                /(AREA_MIX(I)+AREA_ICE_1(I)),                     &
     &                QCL(I)*AREA_MIX(I)/CFLIQ(I)),0.0)
             ELSE
! Deposition does not remove any liquid water content
               DQIL=0.0
             ENDIF
! Adjust liquid content (deposits before vapour by Bergeron Findeison
!  process).
             QCL(I)=QCL(I)-DQIL
             T(I)=T(I)+LFRCP*DQIL
             DQI=DQI_DEP+DQI_SUB-DQIL
! Adjust vapour content
             Q(I)=Q(I)-DQI
             T(I)=T(I)+LSRCP*DQI
! END IF for QCF >  M0.
           END IF

! ----------------------------------------------------------------------
!  7.5  Autoconversion of ice crystals to aggregates
! ----------------------------------------------------------------------
           IF (L_mcr_qcf2 .AND. QCF_CRY(I) >  M0) THEN

             ! Simple explicit Kessler type param. of autoconversion
             ! Autoconversion rate from Lin et al. (1983)
             ! QCFAUTORATE = 0.005*EXP(0.025*(T(I)-273.16))

             ! Set autoconversion limit to emulate split-ice scheme
             QCFAUTOLIM = (QCF_AGG(I)+QCF_CRY(I))                       &
     &               *MAX(EXP(-T_SCALING*MAX((T(I)-CTTEMP(I)),0.0)      &
     &               *MAX(QCF_AGG(I)+QCF_CRY(I),0.0)*QCF0) , 0.0)

             ! Set rate to emulate spilt-ice scheme, i.e. infinite
             QCFAUTORATE = 1./TIMESTEP

             QC  = MIN(QCFAUTOLIM,QCF_CRY(I))
             DPR = MIN(QCFAUTORATE*TIMESTEP*(QCF_CRY(I)-QC)             &
     &                 ,QCF_CRY(I)-QC)

             ! End of calculation of autoconversion amount DPR
             QCF_CRY(I) = QCF_CRY(I)-DPR
             QCF_AGG(I) = QCF_AGG(I)+DPR

           END IF  ! on autoconversion of crystals to agg


! ----------------------------------------------------------------------
!      Transfer processes only active at T less than 0 deg C
! ----------------------------------------------------------------------
         IF(T(I) <  ZERODEGC) THEN
!
! ----------------------------------------------------------------------
!  8   Riming of snow by cloud water -implicit in QCL
! ----------------------------------------------------------------------
           IF (QCF_CRY(I) >  M0.AND.QCL(I) >  0.0                       &
     &         .AND.AREA_MIX(I) >  0.0.AND.CFLIQ(I) >  0.0) THEN
! ----------------------------------------------------------------------
!  8.1 Crystals first
! ----------------------------------------------------------------------
! Calculate water content of mixed phase region
             QCLNEW=QCL(I)/(CFLIQ(I)+CFLIQ(I)*CONSTP(9)*TCGC(I)*CORR(I) &
     &                *TIMESTEP*(RHO(I)*QCF_CRY(I)*CFICEI(I)            &
     &                *CONSTP(5)*TCGCI(I))**CX(6))
! Convert to new grid box total water content
             QCLNEW=QCL(I)*AREA_LIQ(I)/CFLIQ(I)+QCLNEW*AREA_MIX(I)
! Recalculate water contents
               QCF_CRY(I)=QCF_CRY(I)+(QCL(I)-QCLNEW)
               T(I)=T(I)+LFRCP*(QCL(I)-QCLNEW)
               QCL(I)=QCLNEW
! END IF for QCF >  M0.AND.QCL(I) >  0.0
           END IF
!
           IF (QCF_AGG(I) >  M0.AND.QCL(I) >  0.0                       &
     &         .AND.AREA_MIX(I) >  0.0.AND.CFLIQ(I) >  0.0) THEN
! ----------------------------------------------------------------------
!  8.2 Aggregates next
! ----------------------------------------------------------------------
! Calculate water content of mixed phase region
              QCLNEW=QCL(I)/(CFLIQ(I)+CFLIQ(I)*CONSTP(29)*TCG(I)*CORR(I)&
     &                *TIMESTEP*(RHO(I)*QCF_AGG(I)*CFICEI(I)            &
     &                *CONSTP(25)*TCGI(I))**CX(26))
! Convert to new grid box total water content
             QCLNEW=QCL(I)*AREA_LIQ(I)/CFLIQ(I)+QCLNEW*AREA_MIX(I)
! Recalculate water contents
               QCF_AGG(I)=QCF_AGG(I)+(QCL(I)-QCLNEW)
               T(I)=T(I)+LFRCP*(QCL(I)-QCLNEW)
               QCL(I)=QCLNEW
! END IF for QCF >  M0.AND.QCL(I) >  0.0
           END IF
!
! ----------------------------------------------------------------------
!  9   Capture of rain by snow
! ----------------------------------------------------------------------
           IF (QCF_CRY(I) >  M0 .AND. RAIN(I)  >   0.0                  &
     &          .AND. (RAIN_MIX(I)+RAIN_ICE(I)) >  0.0) THEN
! ----------------------------------------------------------------------
!  9.1 Crystals first
! ----------------------------------------------------------------------
! Calculate velocities
             VR1=CORR(I)*CONSTP(41)/6.0*                                &
     &              (RAIN(I)/(RAINFRAC(I)*CONSTP(42)*CORR(I)))**CX(41)
             VS1=CONSTP(4)*CORR(I)*(RHO(I)*QCF_CRY(I)*CFICEI(I)         &
     &           *CONSTP(5)*TCGCI(I))**CX(3)
! Estimate the mean absolute differences in velocities.
             FV1=MAX(ABS(VR1-VS1),(VR1+VS1)/8.0)
! Calculate functions of slope parameter lambda
             LAMR1=(RAIN(I)/(RAINFRAC(I)*CONSTP(42)*CORR(I)))**(CX(42))
             LAMS1=(RHO(I)*QCF_CRY(I)*CFICEI(I)                         &
     &               *CONSTP(5)*TCGCI(I))**(-CX(7))
             LAMFAC1=CONSTP(10)*CONSTP(43)*                             &
     &                 (LAMR1**CX(43)*LAMS1**CX(8)) +                   &
     &               CONSTP(11)*CONSTP(44)*                             &
     &                 (LAMR1**CX(44)*LAMS1**CX(9)) +                   &
     &               CONSTP(12)*CONSTP(45)*                             &
     &                 (LAMR1**CX(45)*LAMS1**CX(10))
! Calculate transfer
          DPR=TCGC(I)*CONSTP(13)*LAMS1**(-CX(11))*LAMR1**(-CX(46))*FV1* &
     &       LAMFAC1*TIMESTEP*RHOR(I)*(RAIN_MIX(I)+RAIN_ICE(I))
             DPR=MIN(DPR,RAIN(I)*(DHI(I)*LSITER)*RHOR(I)                &
     &                   *(RAIN_MIX(I)+RAIN_ICE(I))/RAINFRAC(I))
! Adjust ice and rain contents
             QCF_CRY(I)=QCF_CRY(I)+DPR
             RAIN(I)=RAIN(I)-DPR*RHO(I)*DHILSITERR(I)
             T(I)=T(I)+LFRCP*DPR
             RAINFRAC(I)=RAINFRAC(I)*RAIN(I)/                           &
     &                   (RAIN(I)+DPR*RHO(I)*DHILSITERR(I))
! Update rain fractions
             RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
             RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
             RAIN_ICE(I)=                                               &
     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
             RAIN_CLEAR(I)=                                             &
     &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
! Endif for RAIN >  0.0 in capture term
           END IF
!
! ----------------------------------------------------------------------
!  9.2  Aggregates next
! ----------------------------------------------------------------------
           IF (QCF_AGG(I) >  M0 .AND. RAIN(I)  >   0.0                  &
     &          .AND. (RAIN_MIX(I)+RAIN_ICE(I)) >  0.0) THEN
! Calculate velocities
             VR1=CORR(I)*CONSTP(41)/6.0*                                &
     &              (RAIN(I)/(RAINFRAC(I)*CONSTP(42)*CORR(I)))**CX(41)
             VS1=CONSTP(24)*CORR(I)*(RHO(I)*QCF_AGG(I)*CFICEI(I)        &
     &           *CONSTP(25)*TCGI(I))**CX(23)
! Estimate the mean absolute differences in velocities.
             FV1=MAX(ABS(VR1-VS1),(VR1+VS1)/8.0)
! Calculate functions of slope parameter lambda
             LAMR1=(RAIN(I)/(RAINFRAC(I)*CONSTP(42)*CORR(I)))**(CX(42))
             LAMS1=(RHO(I)*QCF_AGG(I)*CFICEI(I)                         &
     &             *CONSTP(25)*TCGI(I))**(-CX(27))
             LAMFAC1=CONSTP(30)*CONSTP(43)*                             &
     &                 (LAMR1**CX(43)*LAMS1**CX(28)) +                  &
     &               CONSTP(31)*CONSTP(44)*                             &
     &                 (LAMR1**CX(44)*LAMS1**CX(29)) +                  &
     &               CONSTP(32)*CONSTP(45)*                             &
     &                 (LAMR1**CX(45)*LAMS1**CX(30))
! Calculate transfer
          DPR=TCG(I)*CONSTP(33)*LAMS1**(-CX(31))*LAMR1**(-CX(46))*FV1*  &
     &       LAMFAC1*TIMESTEP*RHOR(I)*(RAIN_MIX(I)+RAIN_ICE(I))
             DPR=MIN(DPR,RAIN(I)*(DHI(I)*LSITER)*RHOR(I)                &
     &                   *(RAIN_MIX(I)+RAIN_ICE(I))/RAINFRAC(I))
! Adjust ice and rain contents
             QCF_AGG(I)=QCF_AGG(I)+DPR
             RAIN(I)=RAIN(I)-DPR*RHO(I)*DHILSITERR(I)
             T(I)=T(I)+LFRCP*DPR
             RAINFRAC(I)=RAINFRAC(I)*RAIN(I)/                           &
     &                   (RAIN(I)+DPR*RHO(I)*DHILSITERR(I))
! Update rain fractions
             RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
             RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
             RAIN_ICE(I)=                                               &
     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
             RAIN_CLEAR(I)=                                             &
     &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
!      Endif for RAIN >  0.0 in capture term
           END IF
! ----------------------------------------------------------------------
!      End of transfer processes only active at T less than 0 deg C
! ----------------------------------------------------------------------
         END IF
! ----------------------------------------------------------------------
!  10  Evaporate melting snow
! ----------------------------------------------------------------------
         IF(QCF_CRY(I) >  M0.AND.T(I) >  ZERODEGC)THEN
! ----------------------------------------------------------------------
!  10.1 Crystals first
! ----------------------------------------------------------------------
! Calculate transfer as a function of QCF, T and specific humidity
           PR02=RHO(I)*MAX(QCF_CRY(I),0.0)*CFICEI(I)*CONSTP(5)*TCGCI(I)
           PR04=((APB4-APB5*T(I))*ESW(I)+APB6*P(I)*T(I)**3)
           DPR=TCGC(I)*CONSTP(6)*T(I)**2*ESW(I)*TIMESTEP*               &
     &     (CONSTP(7)*CORR2(I)*PR02**CX(4)                              &
     &      +CONSTP(8)*ROCOR(I)*PR02**CX(5))/(QSL(I)*RHO(I)*PR04)
! TEMPW is the subsaturation in the ice region
           TEMPW=AREA_ICE(I)*(QSL(I)-Q_ICE(I))
           DPR=DPR*TEMPW
           DPR=MAX(MIN(DPR,TEMPW*LHEAT_CORREC_LIQ(I)),0.0)
! Extra check to see we don't get a negative QCF
           DPR=MIN(DPR,QCF_CRY(I))
! Update values of ice and vapour
           QCF_CRY(I)=QCF_CRY(I)-DPR
           Q(I)=Q(I)+DPR
           T(I)=T(I)-DPR*LSRCP
         END IF
!
         IF(QCF_AGG(I) >  M0.AND.T(I) >  ZERODEGC)THEN
! ----------------------------------------------------------------------
!  10.2 Aggregates next
! ----------------------------------------------------------------------
! Calculate transfer as a function of QCF, T and specific humidity
           PR02=RHO(I)*MAX(QCF_AGG(I),0.0)*CFICEI(I)*CONSTP(25)*TCGI(I)
           PR04=((APB4-APB5*T(I))*ESW(I)+APB6*P(I)*T(I)**3)
           DPR=TCG(I)*CONSTP(26)*T(I)**2*ESW(I)*TIMESTEP*               &
     &     (CONSTP(27)*CORR2(I)*PR02**CX(24)                            &
     &      +CONSTP(28)*ROCOR(I)*PR02**CX(25))/(QSL(I)*RHO(I)*PR04)
! TEMPW is the subsaturation in the ice region
           TEMPW=AREA_ICE(I)*(QSL(I)-Q_ICE(I))
           DPR=DPR*TEMPW
           DPR=MAX(MIN(DPR,TEMPW*LHEAT_CORREC_LIQ(I)),0.0)
! Extra check to see we don't get a negative QCF
           DPR=MIN(DPR,QCF_AGG(I))
! Update values of ice and vapour
           QCF_AGG(I)=QCF_AGG(I)-DPR
           Q(I)=Q(I)+DPR
           T(I)=T(I)-DPR*LSRCP
         END IF
!
!
! If iterative melting is active then we have already done this term.
!
         If (.not. l_it_melting) Then
! ----------------------------------------------------------------------
!  11  Melting of snow. Uses a numerical approximation to the wet bulb
!      temperature.
! ----------------------------------------------------------------------
         IF(QCF_CRY(I) >  M0.AND.T(I) >  ZERODEGC)THEN
! ----------------------------------------------------------------------
!  11.1 Crystals first
! ----------------------------------------------------------------------
! An approximate calculation of wet bulb temperature
! TEMPW represents AVERAGE supersaturation in the ice
! Strictly speaking, we need to do two melting calculations,
! for the two different wet bulb temperatures.
           TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0)*CFICEI(I)
           TEMP7=T(I)-ZERODEGC-TEMPW                                    &
     &           *(TW1+TW2*(P(I)-TW3) - TW4*(T(I)-TW5) )
           TEMP7=MAX(TEMP7,0.0)
! End of wet bulb temp formulations.
           PR02=RHO(I)*QCF_CRY(I)*CFICEI(I)*CONSTP(5)*TCGCI(I)
           DPR=TCGC(I)*CONSTP(14)*TIMESTEP*                             &
     &            (CONSTP(7)*CORR2(I)*PR02**CX(4)                       &
     &         + CONSTP(8)*ROCOR(I)*PR02**CX(5))*RHOR(I)
! Solve implicitly in terms of temperature
           DPR=TEMP7*(1.0-1.0/(1.0+DPR*LFRCP))/LFRCP
           DPR=MIN(DPR,QCF_CRY(I))
! Update values of ice and Rain
           QCF_CRY(I)=QCF_CRY(I)-DPR
           RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
           T(I)=T(I)-LFRCP*DPR
           IF (DPR >  0.0) THEN
             RAINFRAC(I)=MAX(RAINFRAC(I),CFICE(I))
! Update rain fractions
             RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
             RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
             RAIN_ICE(I)=                                               &
     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
             RAIN_CLEAR(I)=                                             &
     &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
           ENDIF
! ENDIF for melting snow
         END IF
!
! ----------------------------------------------------------------------
!  11.1 Aggregates next
! ----------------------------------------------------------------------
         IF(QCF_AGG(I) >  M0.AND.T(I) >  ZERODEGC)THEN
! An approximate calculation of wet bulb temperature
! TEMPW represents AVERAGE supersaturation in the ice
! Strictly speaking, we need to do two melting calculations,
! for the two different wet bulb temperatures.
           TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0)*CFICEI(I)
           TEMP7=T(I)-ZERODEGC-TEMPW                                    &
     &           *(TW1+TW2*(P(I)-TW3) - TW4*(T(I)-TW5) )
           TEMP7=MAX(TEMP7,0.0)
! End of wet bulb temp formulations.
           PR02=RHO(I)*QCF_AGG(I)*CFICEI(I)*CONSTP(25)*TCGI(I)
           DPR=TCG(I)*CONSTP(34)*TIMESTEP*                              &
     &            (CONSTP(27)*CORR2(I)*PR02**CX(24)                     &
     &         + CONSTP(28)*ROCOR(I)*PR02**CX(25))*RHOR(I)
! Solve implicitly in terms of temperature
           DPR=TEMP7*(1.0-1.0/(1.0+DPR*LFRCP))/LFRCP
           DPR=MIN(DPR,QCF_AGG(I))
! Update values of ice and Rain
           QCF_AGG(I)=QCF_AGG(I)-DPR
           RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
           T(I)=T(I)-LFRCP*DPR
           IF (DPR >  0.0) THEN
             RAINFRAC(I)=MAX(RAINFRAC(I),CFICE(I))
! Update rain fractions
             RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
             RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
             RAIN_ICE(I)=                                               &
     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
             RAIN_CLEAR(I)=                                             &
     &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
           ENDIF
! ENDIF for melting snow
         END IF
!
         End If  ! l_it_melting
!
! ----------------------------------------------------------------------
!  12.-1 Evaporate all small rain contents to avoid numerical problems
! ----------------------------------------------------------------------
         IF (RAIN(I)  <   ( QCFMIN*DHILSITERR(I)*RHO(I) ) ) THEN
! Evaporate all this rain
           DPR=RAIN(I)*RHOR(I)*DHI(I)*LSITER
           T(I)=T(I)-LCRCP*DPR
           Q(I)=Q(I)+DPR
           RAIN(I)=0.0
! Update rain fractions
           RAINFRAC(I)=0.0
           RAIN_LIQ(I)=0.0
           RAIN_MIX(I)=0.0
           RAIN_ICE(I)=0.0
           RAIN_CLEAR(I)=0.0
         ENDIF
! ----------------------------------------------------------------------
!      Break loop at this point to use vector multiplication
! ----------------------------------------------------------------------
       ENDDO
!
! ----------------------------------------------------------------------
!  12  Evaporation of rain - implicit in subsaturation
! ----------------------------------------------------------------------
#if defined(VECTLIB)
! Use a condensed points calculation to speed up the code
       KK=0
       DO I=1,POINTS
          IF(RAIN(I) >  0.0.AND.RAINFRAC(I) >  0.0)THEN
            KK=KK+1
            TEMP1(KK)=RAIN(I)/(RAINFRAC(I)*CONSTP(42)*CORR(I))
          ENDIF
       ENDDO
! Define a joint power as real to real operations
! are expensive
       POWER=(CX(47)*CX(48))
       CALL POWR_V(KK,TEMP1,POWER,TEMP2)
       POWER=(CX(49)*CX(48))
       CALL POWR_V(KK,TEMP1,POWER,TEMP1)
! Set condensed points counter back to zero
       KK=0
#endif
       DO I=1,POINTS
         IF(RAIN(I) >  0.0.AND.RAINFRAC(I) >  0.0)THEN
           PR04=((APB4-APB5*T(I))*ESW(I)+APB6*P(I)*T(I)**3)
#if defined(VECTLIB)
! Copy temp values to LAMR2 and LAMR1
           KK=KK+1
           LAMR2=TEMP2(KK)
           LAMR1=TEMP1(KK)
#else
! Define LAMR1 and LAMR2
           LAMR1=RAIN(I)/(CONSTP(42)*CORR(I)*RAINFRAC(I))
           LAMR2=LAMR1**(CX(47)*CX(48))
           LAMR1=LAMR1**(CX(49)*CX(48))
#endif
! New, consistent evaporation method, with rain fall speed relationship.
           DPR=CONSTP(46)*T(I)**2*ESW(I)*TIMESTEP
           DPR=DPR*( (CONSTP(47)*CORR2(I)*LAMR2)                        &
     &               + (CONSTP(48)*ROCOR(I)*LAMR1) )
! Calculate transfers.
! TEMP7 is grid box mean supersaturation. This provides a limit.
           TEMP7=(Q_ICE(I)-QSL(I))*RAIN_ICE(I)                          &
     &           +(Q_CLEAR(I)-QSL(I))*RAIN_CLEAR(I)
           DPR=DPR*MAX(-TEMP7*LHEAT_CORREC_LIQ(I),0.0)                  &
     &           /(QSL(I)*RHO(I)*PR04+DPR)
           DPR=DPR*RHO(I)*DHILSITERR(I)
! Another limit is on the amount of rain available
           DPR=MIN(DPR,RAIN(I)*                                         &
     &            (RAIN_ICE(I)+RAIN_CLEAR(I))/RAINFRAC(I))
! Update values of rain
           RAIN(I)=RAIN(I)-DPR
! Update vapour and temperature
           Q(I)=Q(I)+DPR*DHI(I)*LSITER*RHOR(I)
           T(I)=T(I)-DPR*LCRCP*DHI(I)*LSITER*RHOR(I)
           RAINFRAC(I)=RAINFRAC(I)*RAIN(I)/(RAIN(I)+DPR)
! Update rain fractions
           RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
           RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
           RAIN_ICE(I)=                                                 &
     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
           RAIN_CLEAR(I)=RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
! END IF for evaporation of rain.
         END IF
       END DO
!
! ----------------------------------------------------------------------
!  13  Accretion of cloud on rain - implicit in liquid water content
! ----------------------------------------------------------------------
       KK=0
       DO I=1,POINTS
         IF(RAINFRAC(I) >  0.0.AND.RAIN(I) >  0.0                       &
     &      .AND.CFLIQ(I) >  0.0)THEN
           KK=KK+1
! New accretion formulation.
             TEMP1(KK)=RAIN(I)/(RAINFRAC(I)*CONSTP(42)*CORR(I))
! QCLNEW is the in cloud water content left
         END IF
       END DO
#if defined(VECTLIB)
           CALL POWR_V(KK,TEMP1,CX(50),TEMP1)
#endif
       KK=0
       DO I=1,POINTS
         IF(RAINFRAC(I) >  0.0.AND.RAIN(I) >  0.0                       &
     &      .AND.CFLIQ(I) >  0.0)THEN
           KK=KK+1
#if defined(VECTLIB)
           QCLNEW=QCL(I)/((CFLIQ(I)+CFLIQ(I)*CONSTP(49)*CORR(I)         &
     &                             *TIMESTEP*TEMP1(KK)))
#else
           QCLNEW=QCL(I)/((CFLIQ(I)+CFLIQ(I)*CONSTP(49)*CORR(I)         &
     &                             *TIMESTEP*TEMP1(KK)**CX(50)))
#endif
! Convert QCLNEW to a gridbox mean
! TEMP7 is overlap of liquid cloud but no rain
           TEMP7=MAX(CFLIQ(I)-RAIN_LIQ(I)-RAIN_MIX(I),0.0)
           QCLNEW=QCLNEW*(RAIN_LIQ(I)+RAIN_MIX(I))+QCL(I)/CFLIQ(I)*TEMP7
! Now calculate increments to rain.
           RAIN(I)=RAIN(I)+(QCL(I)-QCLNEW)*RHO(I)*DHILSITERR(I)
           QCL(I)=QCL(I)-(QCL(I)-QCLNEW)
! END IF for accretion of cloud on rain.
         END IF
       END DO
!
! ----------------------------------------------------------------------
!  14  Autoconversion of cloud to rain
! ----------------------------------------------------------------------
! This looks complicated because there are several routes through
! depending on whether we are using a second indirect effect,
! have a land or sea point or are using the T3E optimizations. Look
! at UMDP 26 to see what is going on.
!
       KK=0 ! Index counter
! ----------------------------------------------------------------------
!  14.1 Calculate the droplet concentration.
! ----------------------------------------------------------------------
! There are several possible methods.
!
       IF (L_AUTOCONV_MURK .AND. L_MURK) THEN ! Use murk aerosol
         DO I=1,POINTS ! Loop b2
           IF (QCL(I) >  0.0.AND.CFLIQ(I) >  0.0) THEN ! Proceed
             KK=KK+1
! Convert aerosol mass to droplet number. See subroutine VISBTY
             N_DROP(KK) = MAX(AEROSOL(I)/M0_MURK*1.0E-9, 0.0001)
!              1.0E-9 converts from ug/kg to kg/kg
#if defined(VECTLIB)
           END IF ! Proceed
         END DO ! Loop b2
         CALL POWR_V(KK,N_DROP,POWER_MURK,N_DROP)
         DO I=1,KK
           N_DROP(I)=N0_MURK*N_DROP(I)
         END DO
#else
             N_DROP(KK)=N0_MURK*N_DROP(KK)**POWER_MURK
           END IF ! Proceed
         END DO ! Loop b2
#endif
!
       ELSE ! Calculate number of droplets either from the sulphate
!             (and optionally other) aerosol(s) or from a fixed value
!
         KK=0
         DO I=1, POINTS ! Loop x
           IF (QCL(I)  >   0.0 .AND. CFLIQ(I)  >   0.0) THEN ! Proceed
             IF (L_SEASALT_CCN) THEN
               SEA_SALT_PTR=I
             ELSE
               SEA_SALT_PTR=1
             ENDIF
             IF (L_BIOMASS_CCN) THEN
               BIOMASS_PTR=I
             ELSE
               BIOMASS_PTR=1
             ENDIF
             IF (L_OCFF_CCN) THEN
               OCFF_PTR=I
             ELSE
               OCFF_PTR=1
             ENDIF
             IF (L_biogenic_CCN) THEN
               biogenic_ptr=I
             ELSE
               biogenic_ptr=1
             ENDIF
             KK=KK+1
! DEPENDS ON: number_droplet
               N_DROP(KK)=NUMBER_DROPLET(L_USE_SULPHATE_AUTOCONV        &
     &                         , .FALSE., SO4_AIT(I), SO4_ACC(I)        &
     &                         , SO4_DIS(I), L_SEASALT_CCN              &
     &                         , SEA_SALT_FILM(SEA_SALT_PTR)            &
     &                         , SEA_SALT_JET(SEA_SALT_PTR)             &
     &                         , L_biogenic_CCN, biogenic(biogenic_ptr) &
     &                         , L_BIOMASS_CCN                          &
     &                         , BMASS_AGD(BIOMASS_PTR)                 &
     &                         , BMASS_CLD(BIOMASS_PTR)                 &
     &                         , L_OCFF_CCN                             &
     &                         , OCFF_AGD(OCFF_PTR)                     &
     &                         , OCFF_CLD(OCFF_PTR)                     &
     &                         , RHO(I), SNOW_DEPTH(I), LAND_FRACT(I)   &
     &                         , Ntot_land, Ntot_sea )
           ENDIF ! Proceed
         END DO ! Loop x
!
      ENDIF
! ----------------------------------------------------------------------
!  14.2 Calculate whether a sufficient concentration of large droplets
!      is present to allow autoconversion. Calculate its rate.
! ----------------------------------------------------------------------
#if defined(VECTLIB)
       KK=0
       DO I=1,POINTS ! Loop d
            IF (QCL(I) >  0.0.AND.CFLIQ(I) >  0.0) THEN ! Proceed
              KK=KK+1
! Calculate intermediate temporaries for vector operations
              TEMP1(KK)=(RHO(I)*QCL(I)/CFLIQ(I))
              TEMP2(KK)=TEMP1(KK)/N_DROP(KK)
              TEMP4(KK)=N_DROP(KK)
            END IF ! Proceed
       END DO ! Loop d
       TEMP7=POWER_QCL_AUTO-1.0
! Vector real to real operations
       CALL POWR_V(KK,TEMP1,TEMP7,TEMP1)
       CALL POWR_V(KK,TEMP2,-0.33333,TEMP2)
       CALL POWR_V(KK,TEMP4,POWER_DROPLET_AUTO,TEMP4)
#endif
       IF (L_AUTO_DEBIAS) THEN
          DO I=1,POINTS
             T_L(I)=T(I)-(LCRCP*QCL(I))
          END DO
! DEPENDS ON: qsat_wat
          CALL QSAT_WAT(QSL_TL,T_L,P,POINTS)
       ENDIF
       KK=0
       DO I=1,POINTS ! Loop e
         IF (QCL(I) >  0.0.AND.CFLIQ(I) >  0.0) THEN ! Proceed
           KK=KK+1
#if defined(VECTLIB)
! Calculate inverse of droplet mean size
           R_MEAN(KK)=R_MEAN0*TEMP2(KK)
! Calculate numerical factors
           B_FACTOR(KK)=3.0*R_MEAN(KK)
           A_FACTOR(KK)=N_DROP(KK)*0.5
           TEMP3(KK)=-B_FACTOR(KK)*R_AUTO
         END IF ! Proceed
       END DO ! Loop e
! Vector exponential function
! DEPENDS ON: exp_v
       CALL EXP_V(KK,TEMP3,TEMP3)
       KK=0
       DO I=1,POINTS ! Loop e2
         IF (QCL(I) >  0.0.AND.CFLIQ(I) >  0.0) THEN ! Proceed
           KK=KK+1
! Calculate droplet number concentration greater than threshold radius
           N_GT_20=(A_FACTOR(KK)*B_FACTOR(KK)**2)*R_AUTO**2*TEMP3(KK)   &
     &         +(2.*A_FACTOR(KK)*B_FACTOR(KK))*R_AUTO*TEMP3(KK)         &
     &         +(2.*A_FACTOR(KK))*TEMP3(KK)
! Test to see if there is a sufficient concentration of
! droplets >20um for autoconversion to proceed (i.e. >1 N_AUTO):
           IF (N_GT_20  >=  N_AUTO) THEN ! Number
! Calculate an autoconversion rate constant
              AUTORATE=CONSTS_AUTO*EC_AUTO*TEMP4(KK)
           ELSE ! Number
              AUTORATE=0.0
           END IF ! Number
!
#else
! Non T3E code. This is rather easier to follow.
!
! Calculate inverse of mean droplet size
           R_MEAN(KK)=R_MEAN0*(RHO(I)*QCL(I)                            &
     &                /(CFLIQ(I)*N_DROP(KK)))**(-1.0/3.0)
! Calculate numerical factors
           B_FACTOR(KK)=3.0*R_MEAN(KK)
           A_FACTOR(KK)=N_DROP(KK)*0.5
! Calculate droplet number concentration greater than threshold radius
           N_GT_20=(B_FACTOR(KK)**2*A_FACTOR(KK))*R_AUTO**2             &
     &                                   *EXP(-B_FACTOR(KK)*R_AUTO)     &
     &            +(2.*B_FACTOR(KK)*A_FACTOR(KK))*R_AUTO                &
     &                                   *EXP(-B_FACTOR(KK)*R_AUTO)     &
     &            +(2.*A_FACTOR(KK))                                    &
     &                                   *EXP(-B_FACTOR(KK)*R_AUTO)
! Test to see if there is a sufficient concentration of
! droplets >20um for autoconversion to proceed (i.e. > N_AUTO):
           IF (N_GT_20  >=  N_AUTO) THEN ! Number
              AUTORATE=CONSTS_AUTO*EC_AUTO                              &
     &                 *N_DROP(KK)**POWER_DROPLET_AUTO
           ELSE ! Number
              AUTORATE=0.0
           END IF ! Number
!
#endif
! ----------------------------------------------------------------------
!  14.2.1 Optionally de-bias the autoconversion rate.
! ----------------------------------------------------------------------
! If this option is active, the autoconversion rate is corrected ("de-
! biased") based on Wood et al. (Atmos. Res., 65, 109-128, 2002).
! This correction should only be used if the second indirect (or "life-
! time") effect has been selected (i.e. cloud droplet number is calcul-
! ated interactively) and if the diagnostic RHcrit scheme is on (i.e. an
! interactive measure of cloud inhomogeneity is available).
!
           IF (L_AUTO_DEBIAS) THEN

             IF (QCL(I)  <   1.0E-15) THEN  ! Don't bother for very
                                            ! small LWCs - avoids
               AC_FACTOR=1.0                ! AC_FACTOR --> infinity.

             ELSE

               ALPHA_L=EPSILON*LC*QSL_TL(I)/(R*T_L(I)**2)
               A_L=1.0/(1.0+(LCRCP*ALPHA_L))
               SIGMA_S=(1.0-RHCPT(I))*A_L*QSL_TL(I)/SQRT(6.0)
               G_L=1.15*(POWER_QCL_AUTO-1.0)*SIGMA_S
               GACB=EXP(-1.0*QCL(I)/G_L)

               AC_FACTOR=MAX(                                           &
     &                  (CFLIQ(I)**(POWER_QCL_AUTO-1.0))*1.0/(1.0-GACB) &
     &                  ,1.0)

             ENDIF

             AUTORATE=AC_FACTOR*AUTORATE

           ENDIF
!
! ----------------------------------------------------------------------
!  14.3 How much water content can autoconversion remove?
! ----------------------------------------------------------------------
! Calculate value of local liquid water content at which the droplet
! concentration with radii greater than 20um will
! fall below a threshold (1000 m-3), and so determine the minimum
! liquid water content AUTOLIM: This is a numerical approximation
! which ideally could be expressed in terms of N_AUTO and R_AUTO
! but isn't.
           AUTOLIM=(6.20E-31*N_DROP(KK)**3)-(5.53E-22*N_DROP(KK)**2)    &
     &               +(4.54E-13*N_DROP(KK))+(3.71E-6)-(7.59/N_DROP(KK))
! Calculate maximum amount of liquid that can be removed from the
! grid box
           QC=MIN(AUTOLIM*CFLIQ(I)*RHOR(I),QCL(I))
! Calculate autoconversion amount (finally!)
#if defined(VECTLIB)
           DPR=MIN(AUTORATE*TEMP1(KK)                                   &
     &              *TIMESTEP*QCL(I)/CORR2(I),QCL(I)-QC)
#else
           DPR=MIN(AUTORATE                                             &
     &              *(RHO(I)*QCL(I)/CFLIQ(I))**(POWER_QCL_AUTO-1.0)     &
     &              *TIMESTEP*QCL(I)/CORR2(I),QCL(I)-QC)
#endif
! Update liquid water content and rain
           QCL(I)=QCL(I)-DPR
           RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
! Update rain fractions
           IF (DPR >  0.0) THEN ! Rain fraction
             RAINFRAC(I)=MAX(RAINFRAC(I),CFLIQ(I))
             RAIN_LIQ(I)=MIN(AREA_LIQ(I),RAINFRAC(I))
             RAIN_MIX(I)=MIN(AREA_MIX(I),RAINFRAC(I)-RAIN_LIQ(I))
             RAIN_ICE(I)=                                               &
     &            MIN(AREA_ICE(I),RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I))
             RAIN_CLEAR(I)=                                             &
     &            RAINFRAC(I)-RAIN_LIQ(I)-RAIN_MIX(I)-RAIN_ICE(I)
           ENDIF ! Rain fraction
         END IF ! Proceed
! ----------------------------------------------------------------------
!  14.4 Remove small rain amounts.
! ----------------------------------------------------------------------
         IF (RAIN(I) <= (QCFMIN*DHILSITERR(I)*RHO(I))                   &
                                                      ! Tidy up rain
     &       .AND.QCL(I) <= 0.0) THEN
           RAINFRAC(I)=0.0
           Q(I)=Q(I)+RAIN(I)*DHI(I)*RHOR(I)*LSITER
           T(I)=T(I)-RAIN(I)*LCRCP*DHI(I)*LSITER*RHOR(I)
           RAIN(I)=0.0
! Update rain fractions
           RAINFRAC(I)=0.0
           RAIN_LIQ(I)=0.0
           RAIN_MIX(I)=0.0
           RAIN_ICE(I)=0.0
           RAIN_CLEAR(I)=0.0
         ENDIF ! Tidy up rain
! ----------------------------------------------------------------------
!  15  Now continue the loops over points and iterations.
! ----------------------------------------------------------------------
! Continue DO loop over points
       END DO ! Loop e or e2
! Continue DO loop over iterations
       END DO ! Iters_do1
!
! Copy contents of SNOWT to SNOW, to fall into next layer down
! Points_do3
      DO I=1,POINTS
        IF (L_mcr_qcf2) THEN ! two ice prognostics
          QCF(I)  = QCF_AGG(I)
          QCF2(I) = QCF_CRY(I)
          SNOW_CRY(I) = SNOWT_CRY(I)
          SNOW_AGG(I) = SNOWT_AGG(I)
        ELSE ! only one ice prognostic, put all snow in to snow_cry
          QCF(I) = QCF_CRY(I) + QCF_AGG(I)
          SNOW_CRY(I) = 0.0  ! Redundant variable
          SNOW_AGG(I) = SNOWT_CRY(I) + SNOWT_AGG(I)
        ENDIF ! on L_mcr_qcf2
! ----------------------------------------------------------------------
!  16.1 Remove any small amount of QCF which is left over to be tidy.
!       If QCF is less than QCFMIN and isn't growing
!       by deposition (assumed to be given by RHCPT) then remove it.
! ----------------------------------------------------------------------
        IF (L_mcr_qcf2) THEN

          ! Aggregates
          IF (QCF(I) <  QCFMIN) THEN
            IF (T(I) >  ZERODEGC .OR.                                   &
     &        (Q_ICE(I)  <=  QS(I) .AND. AREA_MIX(I)  <=  0.0)          &
     &        .OR. QCF(I) <  0.0) THEN
              Q(I)=Q(I)+QCF(I)
              T(I)=T(I)-LSRCP*QCF(I)
              QCF(I)=0.0
            END IF
          END IF
          ! Crystals
          IF (QCF2(I) <  QCFMIN) THEN
            IF (T(I) >  ZERODEGC .OR.                                   &
     &        (Q_ICE(I)  <=  QS(I) .AND. AREA_MIX(I)  <=  0.0)          &
     &        .OR. QCF2(I) <  0.0) THEN
              Q(I)=Q(I)+QCF2(I)
              T(I)=T(I)-LSRCP*QCF2(I)
              QCF2(I)=0.0
              ! Update ice cloud top temperature
              IF (QCF(I) == 0.0) CTTEMP(I)=T(I)
            END IF
          END IF

        ELSE  ! only one prognostic is active

          IF (QCF(I) <  QCFMIN) THEN
            IF (T(I) >  ZERODEGC .OR.                                   &
     &        (Q_ICE(I)  <=  QS(I) .AND. AREA_MIX(I)  <=  0.0)          &
     &        .OR. QCF(I) <  0.0) THEN
              Q(I)=Q(I)+QCF(I)
              T(I)=T(I)-LSRCP*QCF(I)
              QCF(I)=0.0
              ! Update ice cloud top temperature
              CTTEMP(I)=T(I)
            END IF
          END IF

        END IF  ! on L_mcr_qcf2

! ----------------------------------------------------------------------
!  16.2 Emergency melting of any excess snow to avoid surface snowfall
!       at high temperatures.
! ----------------------------------------------------------------------
        ! Define rapid melting temperature
        If (l_it_melting) then
          t_rapid_melt = zerodegc + 2.0
        Else
          t_rapid_melt = zerodegc
        End If
!
        ! Melt SNOW_AGG first
        IF (SNOW_AGG(I)  >   0.0 .AND. T(I)  >   (t_rapid_melt)) THEN
! Numerical approximation of wet bulb temperature.
! Similar to the melting calculation in section 11.
          TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0)*CFICEI(I)
          TEMP7=T(I)-ZERODEGC-TEMPW*(TW1+TW2*(P(I)-TW3)-TW4*(T(I)-TW5))
          TEMP7=MAX(TEMP7,0.0)
! End of wet bulb calculation
! Remember that DHI uses the shortened timestep which is why LSITER
! appears in the following statements.
          DPR=TEMP7/(LFRCP*LSITER)
          DPR = MIN(DPR,SNOW_AGG(I)*DHI(I)*RHOR(I))
! Update values of snow and rain
          SNOW_AGG(I) = SNOW_AGG(I) - DPR*RHO(I)*DHIR(I)
          RAIN(I)=RAIN(I)+DPR*RHO(I)*DHIR(I)
          T(I)=T(I)-LFRCP*DPR*LSITER
! END IF for emergency melting of snow
        END IF

        ! If ice crystals prognostic is active, also melt snow_cry
        IF (L_mcr_qcf2 .AND.                                            &
     &      SNOW_CRY(I)  >   0.0 .AND. T(I)  >   (t_rapid_melt)) THEN
          ! Numerical approximation of wet bulb temperature.
          ! Similar to the melting calculation in section 11.
          TEMPW=AREA_ICE(I)*MAX(QSL(I)-Q_ICE(I),0.0)*CFICEI(I)
          TEMP7=T(I)-ZERODEGC-TEMPW*(TW1+TW2*(P(I)-TW3)-TW4*(T(I)-TW5))
          TEMP7=MAX(TEMP7,0.0)
          ! End of wet bulb calculation
          ! Remember that DHI uses the shortened timestep which is
          ! why LSITER appears in the following statements.
          DPR = TEMP7/(LFRCP*LSITER)
          DPR = MIN(DPR, SNOW_CRY(I)*DHI(I)*RHOR(I))
          ! Update values of snow and rain
          SNOW_CRY(I) = SNOW_CRY(I) - DPR*RHO(I)*DHIR(I)
          RAIN(I)     = RAIN(I)     + DPR*RHO(I)*DHIR(I)
          T(I)        = T(I)        - LFRCP*DPR*LSITER
        END IF  ! on emergency melting of SNOW_CRY
!
! END DO for tidying up small amounts of ice and snow.
      END DO ! Points_do3
! ----------------------------------------------------------------------
!  17  End of the LSP_ICE subroutine
! ----------------------------------------------------------------------
      RETURN
      END SUBROUTINE LSP_ICE
#endif
