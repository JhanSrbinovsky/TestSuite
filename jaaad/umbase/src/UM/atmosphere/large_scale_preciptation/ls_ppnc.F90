#if defined(A04_3B) || defined(A04_3C) || defined(A04_3D)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  SUBROUTINES LS_PPN and LS_PPNC------------------------------------
!LL  Purpose:
!LL          LS_PPN and LS_PPNC:
!LL           Calculate large-scale (dynamical) precipitation.
!LL           LS_PPNC is the gather/scatter routine which then
!LL           calls LSP_ICE.
!LL  Note: in all cases, level counters (incl subscripts) run from 1
!LL        (lowest model layer) to Q_LEVELS (topmost "wet" model
!LL        layer) - it is assumed that the bottom Q_LEVELS layers are
!LL        the "wet" layers.
!LL
!LL  Put through fpp on Cray.  Activate *IF definition CRAY if running
!LL  on the Cray.
!LL
!LL  Modification History from Version 4.4
!LL     Version    Date
!LL       4.5      March 98          New Deck        Damian Wilson
!
!         5.1      16-12-99   Change to allow 3D RHcrit diagnostic
!                             in place of 1D parameter. AC Bushell
!         5.2      25/10/00   Reintroduction of tracers. P.Selwood
!         5.2      21-11-00   Allow interfacting with 3C scheme.
!                             Return diagnostics. Damian Wilson
!LL
!         5.2      28-11-00   Pass down sea-salt arrays.
!                                                       A. Jones
!         5.3      29-08-01   Introduce logical to control effect of
!                             sulphate aerosol on autoconversion.
!                                                       A. Jones
!         5.3      24-09-01   Tidy up code (mostly redundant) relating
!                             to the sulphur cycle.     A. Jones
!         5.4      30-08-02   Remove lscav_agedsoot, which is now
!                             calculated in microphys_ctl.  P. Davison
!         5.4      27-07-02   Change cloud fractions to prognostics
!                             for PC2                   D. Wilson
!
!         5.4      25-04-02   Pass land fraction down to LSP_ICE.
!                                                       A. Jones
!         5.5      03-02-03   Pass extra microphysics variables
!                             down to LSP_ICE.           R.M.Forbes
!         5.5      03-02-03   Include extra microphysics variables
!                             (qcf2,qrain,qgraup)          R.M.Forbes
!         6.1      01-08-04   Include variables for prognostic rain
!                             and graupel.             R.M.Forbes
!         6.1      07-04-04   Add biomass smoke to call to LSP_ICE.
!                                                          A. Jones
!         6.1      07-04-04   Pass switch for autoconversion de-biasing
!                             down to LSP_ICE.              A. Jones
!         6.2      22-08-05   Remove commented out code. P.Selwood.
!         6.2      23-11-05   Pass through precip variables from
!                             UMUI. Damian Wilson
!         6.2      11-01-06   Include non-hydrostatic calculation
!                             of air density and model thickness.
!                                                       D. Wilson
!         6.2      31-01-06   Pass process rate diags through. R.Forbes
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL
!LL  Documentation: UM Documentation Paper 26.
!LL
!*L  Arguments:---------------------------------------------------------
!*LL  SUBROUTINE LS_PPNC------------------------------------------------
!*L  Arguments:---------------------------------------------------------
      SUBROUTINE LS_PPNC(                                               &
     & IX,N,TIMESTEP,n_iterations,LSRAIN,LSSNOW,LSSNOW2,LSGRAUP         &
     &,droplet_flux                                                     &
     &,CF,CFL,CFF,QCF,QCL,T                                             &
     &,QCF2,QRAIN,QGRAUP                                                &
      ! Process rate diagnostics
     &,PSDEP,PSAUT,PSACW,PSACR,PSACI,PSMLT,PSMLTEVP                     &
     &,PRAUT,PRACW,PREVP                                                &
     &,PGAUT,PGACW,PGACS,PGMLT                                          &
     &,PIFRW,PIPRM,PIDEP,PIACW,PIACR,PIMLT,PIMLTEVP                     &
     &,PIFALL,PSFALL,PRFALL,PGFALL,PLSET,PLEVPSET                       &
      ! Process rate diagnostic switches
     &,L_PSDEP_diag,L_PSAUT_diag,L_PSACW_diag,L_PSACR_diag              &
     &,L_PSACI_diag,L_PSMLT_diag,L_PSMLTEVP_diag                        &
     &,L_PRAUT_diag,L_PRACW_diag,L_PREVP_diag                           &
     &,L_PGAUT_diag,L_PGACW_diag,L_PGACS_diag,L_PGMLT_diag              &
     &,L_PIFRW_diag,L_PIPRM_diag,L_PIDEP_diag,L_PIACW_diag              &
     &,L_PIACR_diag,L_PIMLT_diag,L_PIMLTEVP_diag                        &
     &,L_PIFALL_diag,L_PSFALL_diag,L_PRFALL_diag,L_PGFALL_diag          &
     &,L_PLSET_diag,L_PLEVPSET_diag                                     &
     &,L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_psd, L_it_melting       &
     &,l_non_hydrostatic, l_mixing_ratio,l_cry_agg_dep,l_droplet_settle &
     &,SO2, L_SULPC_SO2                                                 &
     &,NH3, L_SULPC_NH3                                                 &
     &,SO4_AIT, SO4_ACC, SO4_DIS                                        &
     &,AGED_SOOT, L_SOOT                                                &
     &,BMASS_AGD, BMASS_CLD, L_BIOMASS_CCN                              &
     &,OCFF_AGD, OCFF_CLD, L_OCFF_CCN                                   &
     &,AEROSOL,L_MURK, l_pc2                                            &
     &,LSCAV_SO2, LSCAV_NH3,LSCAV_SO4AIT,LSCAV_SO4ACC,LSCAV_SO4DIS      &
     &,SNOW_DEPTH, LAND_FRACT                                           &
     &,L_USE_SULPHATE_AUTOCONV                                          &
     &,L_AUTO_DEBIAS                                                    &
     &,SEA_SALT_FILM, SEA_SALT_JET                                      &
     &,L_SEASALT_CCN, salt_dim1, salt_dim2, salt_dim_ice                &
     &,L_biogenic_CCN, biogenic, biog_dim_ice                           &
     &,Q,p_theta_levels,layer_thickness,deltaz,rhodz_dry,rhodz_moist    &
     &,row_length,rows,rhc_row_length,rhc_rows                          &
     &,BLAND,CW_SEA,CW_LAND                                             &
     &,RHCRIT                                                           &
     &,VFALL,VFALL2,VFALL_RAIN,VFALL_GRAUP                              &
     &,FRAC_ICE_ABOVE,CTTEMP,RAINFRAC,CX,CONSTP                         &
     &,L_seq_mcr,L_autoc_3b,L_autolim_3b,L_autoconv_murk                &
     &,ec_auto,N_drop_land                                              &
     &,N_drop_sea,N_drop_land_cr,N_drop_sea_cr,Ntot_land, Ntot_sea      &
     &,ai,bi,aic,bic                                                    &
     &)

      IMPLICIT NONE

      INTEGER                                                           &
     & N                                                                &
                ! IN Number of points where pptn non-zero from above
!                    or where CF>CFMIN
     &,n_iterations                                                     &
                     ! Number of iterations for iterative melting
     &,row_length,rows                                                  &
     &,rhc_row_length,rhc_rows                                          &
     &,IX(row_length*rows,2)                                            &
                                ! IN gather/scatter index
     &,salt_dim1                                                        &
                     ! Dimensions of input sea-salt arrays - either
     &,salt_dim2                                                        &
                     ! (row_length,rows) or (1,1), depending on
!                    ! L_SEASALT_CCN being T or F respectively.
     &,salt_dim_ice                                                     &
                     ! Dimension to use for call to LSP_ICE.
     &,biog_dim_ice  !     "     "   "   "   "   "     "   .
      REAL                                                              &
     &CF(row_length,rows)                                               &
                              ! IN Cloud fraction.
     &,CFL(row_length,rows)                                             &
                              ! IN Cloud liquid fraction.
     &,CFF(row_length,rows)                                             &
                              ! IN Cloud ice fraction.
! CF, CFL and CFF are IN/OUT if the PC2 cloud scheme is in use.
!
     &, p_theta_levels(row_length,rows)                                 &
     &, layer_thickness(row_length,rows)                                &
                                          ! IN thickness of layer (Pa)
     &, deltaz(row_length,rows)                                         &
                                          ! IN thickness of layer (m)
     &, rhodz_dry(row_length,rows)                                      &
                                          ! Dry air density
                                          ! * layer thickness (kg m-2)
     &, rhodz_moist(row_length,rows)                                    &
                                          ! Moist air density
                                          ! * layer thickness (kg m-2)

     &,RHCRIT(rhc_row_length,rhc_rows)                                  &
!                       IN Critical humidity for cloud formation.
     &,TIMESTEP                                                         &
                      ! IN Timestep (sec).
     &,CW_SEA                                                           &
                      ! IN threshold cloud liquid water content over sea
!                          for conversion to ppn (kg water per m**3).
     &,CW_LAND        ! IN threshold cloud liq. water content over land
!                          for conversion to ppn (kg water per m**3).

      LOGICAL, INTENT(IN) ::                                            &
     &     BLAND(row_length,rows)                                       &
                          !IN Land/sea mask
     &,    L_MURK                                                       &
                          !IN Aerosol needs scavenging.
     &,    L_SULPC_SO2                                                  &
                          !IN Sulphur Cycle on, tracers scavngd if T
     &,    L_SULPC_NH3                                                  &
                          !IN Sulphur cycle HN3 on, tracers scaved if T
     &,    L_SOOT                                                       &
                          !IN Soot cycle on, tracer scavenged if T
     &,    L_USE_SULPHATE_AUTOCONV                                      &
                                    !IN Sulphate aerosols used in
!                                   !   autoconversion if T (i.e.
!                                   !   second indirect effect)
     &,    L_SEASALT_CCN                                                &
                                    !IN Sea-salt aerosols used for
!                                   !   second indirect effect if T
     &,    L_BIOMASS_CCN                                                &
                                    !IN Biomass smoke aerosols used for
!                                   !   second indirect effect if T
     &,    L_OCFF_CCN                                                   &
                                    !IN Fossil-fuel organic carbon
                                    !   aerosol used for second indirect
                                    !   effect if T
     &,    L_biogenic_CCN                                               &
                                    !IN Biogenic aerosol used for
!                                   !   second indirect effect if T
     &,    L_AUTO_DEBIAS                                                &
!                                   !IN Use autoconversion de-biasing
     &,    L_pc2                                                        &
                                    !IN Use PC2 cloud and condensation
     &,    L_mcr_qcf2                                                   &
                                    !IN Use second prognostic ice if T
     &,    L_mcr_qrain                                                  &
                                    !IN Use prognostic rain if T
     &,    L_mcr_qgraup                                                 &
                                    !IN Use prognostic graupel if T
     &,    L_psd                                                        &
                          !IN Use generic ice particle size distribution
     &,    L_it_melting                                                 &
                                    !IN Use iterative melting
     &,    L_non_hydrostatic                                            &
                                    !IN Use non-hydrostatic
                                    !calculation of model layer depths
     &,    L_mixing_ratio                                               &
                                    !IN q is a mixing ratio
     &,    L_cry_agg_dep                                                &
                                    !IN Limit the supersaturation 
                                    !that can be removed by deposition 
                                    !depending on amount of ice in 
                                    !each category.
     &,    l_droplet_settle         !Allow cloud droplets to settle

      REAL, INTENT(INOUT) ::                                            &
     &     Q(row_length,rows)                                           &
                            ! INOUT Specific humidity (kg water/kg air).
     &     ,QCF(row_length,rows)                                        &
                            ! INOUT Cloud ice (kg per kg air).
     &     ,QCL(row_length,rows)                                        &
                            ! INOUT Cloud liquid water (kg per kg air).
     &     ,QCF2(row_length,rows)                                       &
                            ! INOUT Cloud ice2 (kg per kg air).
     &     ,QRAIN(row_length,rows)                                      &
                            ! INOUT Rain water (kg per kg air).
     &     ,QGRAUP(row_length,rows)                                     &
                            ! INOUT Graupel water (kg per kg air).
     &     ,T(row_length,rows)                                          &
                            ! INOUT Temperature (K).
     &     ,AEROSOL(row_length,rows)                                    &
                            ! INOUT Aerosol (K).
     &     ,LSRAIN(row_length,rows)                                     &
                            !INOUT Surface rainfall rate (kg m^-2 s^-1).
     &     ,LSSNOW(row_length,rows)                                     &
                            !INOUT Surface snowfall rate (kg m^-2 s^-1).
     &     ,LSSNOW2(row_length,rows)                                    &
                            !INOUT layer snowfall rate (kg m^-2 s^-1).
     &     ,LSGRAUP(row_length,rows)                                    &
                            !INOUT layer graupelfall rate (kg m^-2 s^-1)
     &     ,droplet_flux(row_length,rows)                               &
                            !INOUT water droplet flux / kg m^-2 s^-1
     &     ,CTTEMP(row_length,rows)                                     &
                            ! INOUT Ice cloud top temperature (K)
     &     ,RAINFRAC(row_length,rows)                                   &
                            ! INOUT Rain fraction.
     &     ,FRAC_ICE_ABOVE(row_length,rows)                             &
                            ! INOUT Ice fraction from layer above
!                                   water.
     &,VFALL(row_length,rows)                                           &
                                     ! INOUT fall velocity of ice (m per
     &,VFALL2(row_length,rows)                                          &
                                      ! INOUT fall vel. of rain (m/s)
     &,VFALL_RAIN(row_length,rows)                                      &
                                      ! INOUT fall vel. of rain (m/s)
     &,VFALL_GRAUP(row_length,rows)   ! INOUT fall vel. of graupel (m/s)

      REAL,INTENT(INOUT) ::                                             &
                            !INOUT S Cycle tracers & scavngd amoun
     &    SO2(row_length, rows)                                         &
     &   ,NH3(row_length, rows)                                         &
     &   ,SO4_AIT(row_length, rows)                                     &
     &   ,SO4_ACC(row_length, rows)                                     &
     &   ,SO4_DIS(row_length, rows)                                     &
     &   ,BMASS_AGD(row_length, rows)                                   &
     &   ,BMASS_CLD(row_length, rows)                                   &
     &   ,OCFF_AGD(row_length, rows)                                    &
     &   ,OCFF_CLD(row_length, rows)                                    &
     &   ,LSCAV_SO2(row_length, rows)                                   &
     &   ,LSCAV_NH3(row_length, rows)                                   &
     &   ,LSCAV_SO4AIT(row_length, rows)                                &
     &   ,LSCAV_SO4ACC(row_length, rows)                                &
     &   ,LSCAV_SO4DIS(row_length, rows)

      REAL                                                              &
                              !INOUT Soot cycles tracers and scav. amnt
     &    AGED_SOOT(row_length, rows)

!
      REAL                                                              &
     &     SNOW_DEPTH(row_length, rows)                                 &
     &   , LAND_FRACT(row_length, rows)
!
      REAL                                                              &
     &     SEA_SALT_FILM(salt_dim1, salt_dim2)                          &
     &   , SEA_SALT_JET(salt_dim1, salt_dim2)
!
      REAL                                                              &
     &     biogenic(row_length, rows)
!
      Logical                                                           &
                      !, Intent(IN)
     &      L_seq_mcr                                                   &
                               ! Use sequential updating of mphys
     &,     L_autoc_3b                                                  &
                               ! Use 3B autoconversion method
     &,     l_autolim_3b                                                &
                               ! Use fixed 3B values for the
                               ! autoconversion limit
     &,     L_autoconv_murk    ! Use murk aerosol to calc. drop number

      Real                                                              &
                      !, Intent(IN)
     &      ec_auto                                                     &
                               ! Collision coalescence efficiency
     &     ,N_drop_land                                                 &
                               ! Number of droplets over land / m-3
     &     ,N_drop_sea                                                  &
                               ! Number of droplets over sea / m-3
     &     ,N_drop_land_cr                                              &
                               ! N_drop_land ^ (-1/3) / m
     &     ,N_drop_sea_cr                                               &
                               ! N_drop_sea ^ (-1/3) / m
     &     ,Ntot_land                                                   &
                               ! Number of droplets over land / m-3
     &     ,Ntot_sea                                                    &
                               ! Number of droplets over sea / m-3
     &     ,ai, bi, aic, bic                                            
                               ! Ice particle mass size parameters

! Microphysical process rate diagnostics (2D arrays on one level)
! Note: These arrays will only increase memory usage and are
!       only referenced if the particular diagnostic is active

      REAL, Intent(InOut) ::                                            &
     &  PSDEP(row_length,rows)                                          &
                                 ! Deposition of vapour to snow agg.
     &, PSAUT(row_length,rows)                                          &
                                 ! Autoconversion of aggregates from cry
     &, PSACW(row_length,rows)                                          &
                                 ! Accretion of liq. water by snow agg.
     &, PSACR(row_length,rows)                                          &
                                 ! Collection of rain by snow aggregates
     &, PSACI(row_length,rows)                                          &
                                 ! Collection of ice crystals by agg.
     &, PSMLT(row_length,rows)                                          &
                                 ! Melting of snow aggregates
     &, PSMLTEVP(row_length,rows)! Evaporation of melting aggregates
      REAL, Intent(InOut) ::                                            &
     &  PRAUT(row_length,rows)                                          &
                                 ! Autoconversion of cloud drops to rain
     &, PRACW(row_length,rows)                                          &
                                 ! Accretion of liq. water by rain
     &, PREVP(row_length,rows)   ! Evaporation of rain
      REAL, Intent(InOut) ::                                            &
     &  PGAUT(row_length,rows)                                          &
                                 ! Autoconversion of graupel from agg.
     &, PGACW(row_length,rows)                                          &
                                 ! Accretion of liq. water by graupel
     &, PGACS(row_length,rows)                                          &
                                 ! Collection of snow agg. by graupel
     &, PGMLT(row_length,rows)   ! Melting of graupel
      REAL, Intent(InOut) ::                                            &
     &  PIFRW(row_length,rows)                                          &
                                 ! Homogeneous freezing nucleation
     &, PIPRM(row_length,rows)                                          &
                                 ! Heterogeneous (primary) nucleation
     &, PIDEP(row_length,rows)                                          &
                                 ! Deposition of vapour to ice crystals
     &, PIACW(row_length,rows)                                          &
                                 ! Accretion of liq. water by ice cry.
     &, PIACR(row_length,rows)                                          &
                                 ! Collection of rain by ice crystals
     &, PIMLT(row_length,rows)                                          &
                                 ! Melting of ice crystals
     &, PIMLTEVP(row_length,rows)! Evaporation of melting ice crystals
      REAL, Intent(InOut) ::                                            &
     &  PIFALL(row_length,rows)                                         &
                                 ! Sedimentation of ice crystals
     &, PSFALL(row_length,rows)                                         &
                                 ! Sedimentation of aggregates
     &, PRFALL(row_length,rows)                                         &
                                 ! Sedimentation of rain
     &, PGFALL(row_length,rows)  ! Sedimentation of graupel
      REAL, Intent(InOut) ::                                            &
     &  PLSET(row_length,rows)                                          &
                                 ! Droplet settling of liquid water
     &, PLEVPSET(row_length,rows)! Evaporated settled droplets

! Microphysical process rate diagnostic logical switches
      LOGICAL, Intent(In) ::                                            &
     &  L_PSDEP_diag                                                    &
                       ! Deposition of vapour to snow agg.
     &, L_PSAUT_diag                                                    &
                       ! Autoconversion of aggregates from cry
     &, L_PSACW_diag                                                    &
                       ! Accretion of liq. water by snow agg.
     &, L_PSACR_diag                                                    &
                       ! Collection of rain by snow aggregates
     &, L_PSACI_diag                                                    &
                       ! Collection of ice crystals by agg.
     &, L_PSMLT_diag                                                    &
                       ! Melting of snow aggregates
     &, L_PSMLTEVP_diag! Evaporation of melting aggregates
      LOGICAL, Intent(In) ::                                            &
     &  L_PRAUT_diag                                                    &
                       ! Autoconversion of cloud drops to rain
     &, L_PRACW_diag                                                    &
                       ! Accretion of liq. water by rain
     &, L_PREVP_diag   ! Evaporation of rain
      LOGICAL, Intent(In) ::                                            &
     &  L_PGAUT_diag                                                    &
                       ! Autoconversion of graupel from agg.
     &, L_PGACW_diag                                                    &
                       ! Accretion of liq. water by graupel
     &, L_PGACS_diag                                                    &
                       ! Collection of snow agg. by graupel
     &, L_PGMLT_diag   ! Melting of graupel
      LOGICAL, Intent(In) ::                                            &
     &  L_PIFRW_diag                                                    &
                       ! Homogeneous freezing nucleation
     &, L_PIPRM_diag                                                    &
                       ! Heterogeneous (primary) nucleation
     &, L_PIDEP_diag                                                    &
                       ! Deposition of vapour to ice crystals
     &, L_PIACW_diag                                                    &
                       ! Accretion of liq. water by ice cry.
     &, L_PIACR_diag                                                    &
                       ! Collection of rain by ice crystals
     &, L_PIMLT_diag                                                    &
                       ! Melting of ice crystals
     &, L_PIMLTEVP_diag! Evaporation of melting ice crystals
      LOGICAL, Intent(In) ::                                            &
     &  L_PIFALL_diag                                                   &
                       ! Sedimentation of ice crystals
     &, L_PSFALL_diag                                                   &
                       ! Sedimentation of aggregates
     &, L_PRFALL_diag                                                   &
                       ! Sedimentation of rain
     &, L_PGFALL_diag  ! Sedimentation of graupel
      LOGICAL, Intent(In) ::                                            &
     &  L_PLSET_diag                                                    &
                       ! Droplet settling of liquid water
     &, L_PLEVPSET_diag! Evaporated settled droplets

!*L  Workspace usage ---------------------------------------------------
!
      REAL                                                              &
!ajm     & PSTAR_C(N)        ! gathered Surface pressure (Pa).
     &CF_C(N)                                                           &
                        ! gathered Cloud fraction.
     &,Q_C(N)                                                           &
                         ! gathered Specific humidity (kg water/kg air).
     &,QCF_C(N)                                                         &
                         ! gathered Cloud ice (kg per kg air).
     &,QCL_C(N)                                                         &
                         ! gathered Cloud liquid water (kg per kg air).
     &,QCF2_C(N)                                                        &
                          ! gathered cloud ice2 (kg per kg air).
     &,QRAIN_C(N)                                                       &
                           ! gathered rain (kg per kg air).
     &,QGRAUP_C(N)                                                      &
                            ! gathered graupel (kg per kg air).
     &,T_C(N)                                                           &
                         ! gathered Temperature (K).
     &,AERO_C(N)                                                        &
                         ! gathered Aerosol.
     &,LSRAIN_C(N)                                                      &
                   !gathered Surface rainfall rate (kg per sq m per s).
     &,LSSNOW_C(N)                                                      &
                   !gathered Surface snowfall rate (kg per sq m per s).
     &,LSSNOW2_C(N)                                                     &
                    !gathered layer snowfall rate (kg per sq m per s).
     &,LSGRAUP_C(N)                                                     &
                    !gathered layer graupel fall rate (kg/sq m/s)
     &,droplet_flux_c(N)                                                &
                         ! gathered water droplet flux / kg m-2 s-1
     &,CTTEMP_C(N)                                                      &
                              !gathered ice cloud top temperature.
     &,RAINFRAC_C(N)                                                    &
                              !gathered rain fraction.
     &,FRAC_ICE_ABOVE_C(N)                                              &
                              !gathered fraction of ice in layer above
     &,CFL_C(N)                                                         &
                         ! gathered Cloud liquid fraction.
     &,CFF_C(N)                                                         &
                         ! gathered Cloud ice fraction.
     &,VFALL_C(N)                                                       &
                         ! gathered fall velocity (m per s).
     &,VFALL2_C(N)                                                      &
                          ! gathered fall velocity for qcf2 (m per s).
     &,VFALL_RAIN_C(N)                                                  &
                          ! gathered fall velocity for qcf2 (m per s).
     &,VFALL_GRAUP_C(N)                                                 &
                          ! gathered fall velocity for qcf2 (m per s).
     &,SEA_SALT_FILM_C(salt_dim_ice)                                    &
                                    ! gathered film-mode sea-salt (m-3)
     &,SEA_SALT_JET_C(salt_dim_ice)                                     &
                                    ! gathered jet-mode sea-salt (m-3)
     &,biogenic_C(biog_dim_ice)                                         &
                                    ! gathered biogenic aerosol (m.m.r.)
     &,RHC_C(N)          ! gathered RH_crit value at points.

      REAL                                                              &
                         ! gathered sulphate aerosol arrays
     &    SO4_AIT_C(N)                                                  &
     &   ,SO4_ACC_C(N)                                                  &
     &   ,SO4_DIS_C(N)
!
      REAL                                                              &
                         ! gathered biomass aerosol arrays
     &    BMASS_AGD_C(N)                                                &
     &   ,BMASS_CLD_C(N)
!
      REAL                                                              &
                         ! gathered fossil-fuel org carb aerosol arrays
     &    OCFF_AGD_C(N)                                                 &
     &   ,OCFF_CLD_C(N)

!
      REAL                                                              &
     &   SNOW_DEPTH_C(N)                                                &
     &  ,LAND_FRACT_C(N)

!
      REAL                                                              &
     & RHODZ(N)                                                         &
                      ! WORK Used for air mass p.u.a. in successive
!                            layers.
     &,deltaz_c(n)                                                      &
                       ! Thickness of layer (m)
     &,rhodz_dry_c(n)                                                   &
                       ! Dry air density * layer thickness (kg m-2)
     &,rhodz_moist_c(n)                                                 &
                       ! Moist air density * layer thickness (kg m-2)

     &,P(N)           ! WORK Used for pressure at successive levels.
      LOGICAL BLAND_C(N)          ! gathered land/sea mask
!
! Microphysical process rate diagnostics (compressed arrays)
      REAL                                                              &
     &  PSDEP_C(N)                                                      &
                    ! Deposition of vapour to snow aggregates
     &, PSAUT_C(N)                                                      &
                    ! Autoconversion of aggregates from crystals
     &, PSACW_C(N)                                                      &
                    ! Accretion of liq. water by snow aggregates
     &, PSACR_C(N)                                                      &
                    ! Collection of rain by snow aggregates
     &, PSACI_C(N)                                                      &
                    ! Collection of ice crystals by aggregates
     &, PSMLT_C(N)                                                      &
                    ! Melting of snow aggregates
     &, PSMLTEVP_C(N)  ! Evaporation of melting aggregates
      REAL                                                              &
     &  PRAUT_C(N)                                                      &
                    ! Autoconversion of cloud drops to rain
     &, PRACW_C(N)                                                      &
                    ! Accretion of liq. water by rain
     &, PREVP_C(N)  ! Evaporation of rain
      REAL                                                              &
     &  PGAUT_C(N)                                                      &
                    ! Autoconversion of graupel from aggregates
     &, PGACW_C(N)                                                      &
                    ! Accretion of liq. water by graupel
     &, PGACS_C(N)                                                      &
                    ! Collection of snow aggregates by graupel
     &, PGMLT_C(N)  ! Melting of graupel
      REAL                                                              &
     &  PIFRW_C(N)                                                      &
                    ! Homogeneous freezing nucleation
     &, PIPRM_C(N)                                                      &
                    ! Heterogeneous (primary) nucleation
     &, PIDEP_C(N)                                                      &
                    ! Deposition of vapour to ice crystals
     &, PIACW_C(N)                                                      &
                    ! Accretion of liq. water by ice crystals
     &, PIACR_C(N)                                                      &
                    ! Collection of rain by ice crystals
     &, PIMLT_C(N)                                                      &
                    ! Melting of ice crystals
     &, PIMLTEVP_C(N)  ! Evaporation of melting ice crystals
      REAL                                                              &
     &  PIFALL_C(N)                                                     &
                    ! Sedimentation of ice crystals
     &, PSFALL_C(N)                                                     &
                    ! Sedimentation of aggregates
     &, PRFALL_C(N)                                                     &
                    ! Sedimentation of rain
     &, PGFALL_C(N) ! Sedimentation of graupel
      REAL                                                              &
     &  PLSET_C(N)                                                      &
                    ! Droplet settling of liquid water
     &, PLEVPSET_C(N) ! Evaporated settled droplets

! Call size of CX and CONSTP
#include "c_lspsiz.h"
! Call comdeck containing ls ppn scavenging coeffs for Sulphur Cycle
#include "c_sullsp.h"
!
!  External subroutines called -----------------------------------------
      EXTERNAL LSP_ICE,LSP_SCAV
!     &        ,SLSPSCV
!*----------------------------------------------------------------------
!  Physical constants -------------------------------------------------
#include "c_g.h"
      REAL P1UPONG
      PARAMETER (                                                       &
     & P1UPONG=1./G                                                     &
                           ! One upon g (sq seconds per m).
     &)
!  Define local variables ----------------------------------------------
      INTEGER                                                           &
     & MULTRHC         ! Zero if (rhc_row_length*rhc_rows) le 1, else 1
!
      INTEGER I,II,IJ                                                   &
                            ! Loop counters: I - horizontal field index;
     &         ,IRHI,IRHJ   !           : IRHI,IRHJ-indices for RHcrit.
!
!-----------------------------------------------------------------------
!L Internal structure.
!L 1. gather variables using index
!-----------------------------------------------------------------------
      IF ( (rhc_row_length * rhc_rows)  >   1) THEN
        MULTRHC = 1
      ELSE
        MULTRHC = 0
      END IF
!
      DO I=1,N
         II=IX(I,1)
         IJ=IX(I,2)
        IRHI = (MULTRHC * (II - 1)) + 1
        IRHJ = (MULTRHC * (IJ - 1)) + 1
        LSRAIN_C(I)=LSRAIN(II,IJ)
        LSSNOW_C(I)=LSSNOW(II,IJ)
        LSSNOW2_C(I) = LSSNOW2(II,IJ)
        LSGRAUP_C(I) = LSGRAUP(II,IJ)
        droplet_flux_c(I) = droplet_flux(II,IJ)
        CTTEMP_C(I)=CTTEMP(II,IJ)
        RAINFRAC_C(I)=RAINFRAC(II,IJ)
        FRAC_ICE_ABOVE_C(I)=FRAC_ICE_ABOVE(II,IJ)
!ajm        PSTAR_C(I) =PSTAR(II,IJ)
        P(I)=p_theta_levels(II,IJ)
        RHODZ(I)=-P1UPONG*layer_thickness(II,IJ)
        deltaz_c(i)    = deltaz(ii,ij)
        rhodz_dry_c(i)   = rhodz_dry(ii,ij)
        rhodz_moist_c(i) = rhodz_moist(ii,ij)
        BLAND_C(I) =BLAND(II,IJ)
        CF_C(I)=CF(II,IJ)
        CFL_C(I)=CFL(II,IJ)
        CFF_C(I)=CFF(II,IJ)
        QCF_C(I)=QCF(II,IJ)
        QCL_C(I)=QCL(II,IJ)
        Q_C(I)=Q(II,IJ)
        T_C(I)=T(II,IJ)
        IF (L_mcr_qcf2) QCF2_C(I) = QCF2(II,IJ)
        IF (L_mcr_qrain) QRAIN_C(I) = QRAIN(II,IJ)
        IF (L_mcr_qgraup) QGRAUP_C(I) = QGRAUP(II,IJ)
        IF (L_MURK) AERO_C(I)=AEROSOL(II,IJ)
        VFALL_C(I)=VFALL(II,IJ)
        VFALL2_C(I) = VFALL2(II,IJ)
        VFALL_RAIN_C(I) = VFALL_RAIN(II,IJ)
        VFALL_GRAUP_C(I) = VFALL_GRAUP(II,IJ)
        RHC_C(I)=RHCRIT(IRHI,IRHJ)
        IF (L_USE_SULPHATE_AUTOCONV) THEN
          SO4_AIT_C(I)=SO4_AIT(II,IJ)
          SO4_ACC_C(I)=SO4_ACC(II,IJ)
          SO4_DIS_C(I)=SO4_DIS(II,IJ)
        END IF
        IF (L_SEASALT_CCN) THEN
          SEA_SALT_FILM_C(I)=SEA_SALT_FILM(II,IJ)
          SEA_SALT_JET_C(I)=SEA_SALT_JET(II,IJ)
        ELSE
          SEA_SALT_FILM_C(1)=SEA_SALT_FILM(1,1)
          SEA_SALT_JET_C(1)=SEA_SALT_JET(1,1)
        ENDIF
        IF (L_BIOMASS_CCN) THEN
          BMASS_AGD_C(I)=BMASS_AGD(II,IJ)
          BMASS_CLD_C(I)=BMASS_CLD(II,IJ)
        ELSE
          BMASS_AGD_C(I)=BMASS_AGD(1,1)
          BMASS_CLD_C(I)=BMASS_CLD(1,1)
        ENDIF
        IF (L_biogenic_CCN) THEN
          biogenic_C(I)=biogenic(II,IJ)
        ELSE
          biogenic_C(1)=biogenic(1,1)
        ENDIF
        IF (L_OCFF_CCN) THEN
          OCFF_AGD_C(I)=OCFF_AGD(II,IJ)
          OCFF_CLD_C(I)=OCFF_CLD(II,IJ)
        ELSE
          OCFF_AGD_C(I)=OCFF_AGD(1,1)
          OCFF_CLD_C(I)=OCFF_CLD(1,1)
        ENDIF
        SNOW_DEPTH_C(I)=SNOW_DEPTH(II,IJ)
        LAND_FRACT_C(I)=LAND_FRACT(II,IJ)

        ! Process diagnostic arrays are initialised to zero in LSP_ICE
      END DO ! Loop over points
!
!-----------------------------------------------------------------------
! ICE FORMATION/EVAPORATION/MELTING
! WATER CLOUD AND RAIN FORMATION/EVAPORATION
!-----------------------------------------------------------------------
! The call to LSP_ICE replaces the calls to LSP_EVAP, LSPFRMT
! and LSP_FORM.
! CFL_C contains cloud fraction for ice
! CFF_C contains cloud fraction for water
! DEPENDS ON: lsp_ice
          CALL LSP_ICE(P,RHODZ,deltaz_c,rhodz_dry_c,rhodz_moist_c,      &
     &      TIMESTEP,n_iterations,N,L_MURK,RHC_C,                       &
     &      L_USE_SULPHATE_AUTOCONV, L_AUTO_DEBIAS,                     &
     &      L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_psd, L_it_melting, &
     &      l_non_hydrostatic, l_mixing_ratio, l_cry_agg_dep,           &
     &      SO4_ACC_C,SO4_DIS_C,SO4_AIT_C,                              &
     &      L_BIOMASS_CCN, BMASS_AGD_C, BMASS_CLD_C,                    &
     &      L_OCFF_CCN, OCFF_AGD_C, OCFF_CLD_C,                         &
     &      L_SEASALT_CCN, SEA_SALT_FILM_C,                             &
     &      SEA_SALT_JET_C, salt_dim_ice,                               &
     &      L_biogenic_CCN, biogenic_C, biog_dim_ice,                   &
     &      SNOW_DEPTH_C,LAND_FRACT_C,AERO_C,                           &
     &      L_seq_mcr,L_autoc_3b,L_autolim_3b,L_autoconv_murk,          &
     &      l_droplet_settle, ec_auto,                                  &
     &      N_drop_land, N_drop_sea,N_drop_land_cr,N_drop_sea_cr,       &
     &      Ntot_land, Ntot_sea, ai, bi, aic, bic,                      &
     &      QCF_C,QCL_C,Q_C,QCF2_C,QRAIN_C,QGRAUP_C,                    &
     &      LSRAIN_C,VFALL_RAIN_C, LSSNOW_C,VFALL_C,                    &
     &      LSSNOW2_C,VFALL2_C, LSGRAUP_C,VFALL_GRAUP_C, droplet_flux_c,&
     &      FRAC_ICE_ABOVE_C,CTTEMP_C,RAINFRAC_C,                       &
     &      T_C,CF_C,CFL_C,CFF_C,BLAND_C,CX,CONSTP                      &
     &,     PSDEP_C,PSAUT_C,PSACW_C,PSACR_C,PSACI_C,PSMLT_C,PSMLTEVP_C  &
     &,     PRAUT_C,PRACW_C,PREVP_C                                     &
     &,     PGAUT_C,PGACW_C,PGACS_C,PGMLT_C                             &
     &,     PIFRW_C,PIPRM_C,PIDEP_C,PIACW_C,PIACR_C,PIMLT_C,PIMLTEVP_C  &
     &,     PIFALL_C,PSFALL_C,PRFALL_C,PGFALL_C,PLSET_C,PLEVPSET_C      &
     & )
!-----------------------------------------------------------------------
!L 3.4 Lose aerosol by scavenging: call LSP_SCAV
!-----------------------------------------------------------------------
!
      IF (L_MURK)  THEN
! DEPENDS ON: lsp_scav
        CALL LSP_SCAV(TIMESTEP,N,LSRAIN_C,LSSNOW_C,AERO_C)
      ENDIF
!
!-----------------------------------------------------------------------
!L 4  Scatter back arrays which will have been changed.
!L
!-----------------------------------------------------------------------
!

      DO I=1,N
         II=IX(I,1)
         IJ=IX(I,2)
        T(II,IJ)=T_C(I)
        Q(II,IJ)=Q_C(I)
        QCF(II,IJ)=QCF_C(I)
        QCL(II,IJ)=QCL_C(I)
        IF (L_mcr_qcf2) QCF2(II,IJ) = QCF2_C(I)
        IF (L_mcr_qrain) QRAIN(II,IJ) = QRAIN_C(I)
        IF (L_mcr_qgraup) QGRAUP(II,IJ) = QGRAUP_C(I)
        IF (L_MURK) AEROSOL(II,IJ)=AERO_C(I)
        LSRAIN(II,IJ)=LSRAIN_C(I)
        LSSNOW(II,IJ)=LSSNOW_C(I)
        LSSNOW2(II,IJ)=LSSNOW2_C(I)
        LSGRAUP(II,IJ)=LSGRAUP_C(I)
        droplet_flux(ii,ij)=droplet_flux_c(i)
        CTTEMP(II,IJ)=CTTEMP_C(I)
        RAINFRAC(II,IJ)=RAINFRAC_C(I)
        FRAC_ICE_ABOVE(II,IJ)=FRAC_ICE_ABOVE_C(I)
        If (L_pc2) then
          CFF(II,IJ)=CFF_C(I)
          CFL(II,IJ)=CFL_C(I)
          CF(II,IJ)=CF_C(I)
        End If  ! L_pc2

        VFALL(II,IJ)=VFALL_C(I)
        VFALL2(II,IJ)=VFALL2_C(I)
        VFALL_RAIN(II,IJ)=VFALL_RAIN_C(I)
        VFALL_GRAUP(II,IJ)=VFALL_GRAUP_C(I)
        ! Only store process rates in array for diagnostic
        ! if a particular diagnostic is requested,
        ! otherwise overwriting will occur
        ! (space for the 3D array in MCR_CTL is only allocated
        ! if the diagnostic is active, to save memory)
        If (L_PSDEP_diag) PSDEP(II,IJ) = PSDEP_C(I)
        If (L_PSAUT_diag) PSAUT(II,IJ) = PSAUT_C(I)
        If (L_PSACW_diag) PSACW(II,IJ) = PSACW_C(I)
        If (L_PSACR_diag) PSACR(II,IJ) = PSACR_C(I)
        If (L_PSACI_diag) PSACI(II,IJ) = PSACI_C(I)
        If (L_PSMLT_diag) PSMLT(II,IJ) = PSMLT_C(I)
        If (L_PSMLTEVP_diag) PSMLTEVP(II,IJ) = PSMLTEVP_C(I)
        If (L_PRAUT_diag) PRAUT(II,IJ) = PRAUT_C(I)
        If (L_PRACW_diag) PRACW(II,IJ) = PRACW_C(I)
        If (L_PREVP_diag) PREVP(II,IJ) = PREVP_C(I)
        If (L_PGAUT_diag) PGAUT(II,IJ) = PGAUT_C(I)
        If (L_PGACW_diag) PGACW(II,IJ) = PGACW_C(I)
        If (L_PGACS_diag) PGACS(II,IJ) = PGACS_C(I)
        If (L_PGMLT_diag) PGMLT(II,IJ) = PGMLT_C(I)
        If (L_PIFRW_diag) PIFRW(II,IJ) = PIFRW_C(I)
        If (L_PIPRM_diag) PIPRM(II,IJ) = PIPRM_C(I)
        If (L_PIDEP_diag) PIDEP(II,IJ) = PIDEP_C(I)
        If (L_PIACW_diag) PIACW(II,IJ) = PIACW_C(I)
        If (L_PIACR_diag) PIACR(II,IJ) = PIACR_C(I)
        If (L_PIMLT_diag) PIMLT(II,IJ) = PIMLT_C(I)
        If (L_PIMLTEVP_diag) PIMLTEVP(II,IJ) = PIMLTEVP_C(I)
        If (L_PIFALL_diag) PIFALL(II,IJ) = PIFALL_C(I)
        If (L_PSFALL_diag) PSFALL(II,IJ) = PSFALL_C(I)
        If (L_PRFALL_diag) PRFALL(II,IJ) = PRFALL_C(I)
        If (L_PGFALL_diag) PGFALL(II,IJ) = PGFALL_C(I)
        If (L_PLSET_diag) PLSET(II,IJ) = PLSET_C(I)
        If (L_PLEVPSET_diag) PLEVPSET(II,IJ) = PLEVPSET_C(I)
      END DO ! Loop over points
!
!
      RETURN
      END SUBROUTINE LS_PPNC
#endif
