! TYPPTRO start
! needs TYPSIZE included beforehand
!LL History:
!LL Version  Date  Comment
!LL  3.4   18/5/94 Remove sea ice flux correction. J F Thomson.
!LL  4.4   15/06/97 introduce pointers for the free surface solution.
!LL                                                         R.Lenton
!LL  4.5  04/08/97 Add pointers for ocean boundary data in D1 array
!LL                                                C.G. Jones
!LL  4.5    1/07/98 Add new pointer to atmospheric CO2 field. C.D.Jones
!LL  5.3   14/09/01 Add a new pointer for Levitus nutreint. S.Spall
!LL  5.3   15/10/01 Move ocean boundary data (to TYPBNDO) M J Bell
!LL  5.5   28/02/03 Add new pointers for multiple ice category fields
!LL                                                       A.McLaren
!    5.5   27/02/03 Upgrade to FOAM data assimilation. A. Hines
!LL  5.5   28/02/03 Add pointers to Helium source fields. S.Liddicoat
!LL  6.1   20/08/04 Add pointers to local flux fields.  A. Sellar
!    6.1   02/09/04   Implement Pressure Correction scheme. M. Martin.
!    6.2   24/02/06 Add pointer for 10 windspeed A2O     J.Gunson
!    6.2   23/11/03 Modifications for vertical velocity fix. C. Harris
#if defined(OCEAN)
! Pointers for OCEAN model variables. Configuration dependent.

      ! Ocean primary array variables (depends on resolution)
      INTEGER :: joc_tracer(NT,2)        ! Start of each tracer

      ! Scalar variables
      INTEGER :: joc_u(2)                ! Baroclinic u
      INTEGER :: joc_v(2)                ! Baroclinic v
      INTEGER :: joc_stream(2)           ! Stream function

      INTEGER :: joc_FWFLX       ! Fresh water flux (TAU)
      INTEGER :: joc_FWFLXb      ! Fresh water flux (TAU-1)
      INTEGER :: joc_tend(2)             ! Stream func tendency
      INTEGER :: joc_pmepassb  ! pmepass from previous timestep
      INTEGER :: joc_dusp      ! unsmoothed fsh tendency
      INTEGER :: joc_eta     ! Surface elevation current T step
      INTEGER :: joc_etab    ! Surface elevation previous T step

      ! depth integrated x-comp of b'tropic velocity at current T step
      INTEGER :: joc_ubt

      ! depth integrated x-comp of b'tropic velocity - previous b'tropic
      !  T step
      INTEGER :: joc_ubtbbt

      ! depth integrated y-comp of b'tropic velocity at current T step
      INTEGER :: joc_vbt

      ! depth integrated y-comp of b'tropic velocity - previous b'tropic
      ! T step
      INTEGER :: joc_vbtbbt

      ! depth integrated x-comp of b'tropic velocity - previous b'clinic
      ! T step
      INTEGER :: joc_ubtbbc

      ! depth integrated y-comp of b'tropic velocity - previous b'clinic
      ! T step
      INTEGER :: joc_vbtbbc

      INTEGER :: joc_mld                 ! Mixed layer depth
      INTEGER :: joc_athkdft             ! thickness diff coeff
      INTEGER :: joc_varchl1             ! variable C:Chl (phytopl.)
      INTEGER :: joc_varchl2             ! variable C:Chl (diatoms)

      ! Sea ice primary
      INTEGER :: joc_snow                ! Snow depth over ice
      INTEGER :: joc_mischt              ! Misc heat for ice
      INTEGER :: joc_htotoi              ! Ocean to ice heat
      INTEGER :: joc_salinc              ! Ice inc to salinity
      INTEGER :: joc_icy                 ! Logical for ice
      INTEGER :: joc_isx                 ! Ice/ocean stress
      INTEGER :: joc_isy                 ! Ice/ocean stress
      INTEGER :: joc_icecon              ! Ice concentration
      INTEGER :: joc_icedep              ! Mean ice depth
      INTEGER :: joc_iceu                ! Dyn ice u
      INTEGER :: joc_icev                ! Dyn ice v
      INTEGER :: joc_isig11ne        !Internal ice stress cpt
      INTEGER :: joc_isig11se        !Internal ice stress cpt
      INTEGER :: joc_isig11sw        !Internal ice stress cpt
      INTEGER :: joc_isig11nw        !Internal ice stress cpt
      INTEGER :: joc_isig12ne        !Internal ice stress cpt
      INTEGER :: joc_isig12se        !Internal ice stress cpt
      INTEGER :: joc_isig12sw        !Internal ice stress cpt
      INTEGER :: joc_isig12nw        !Internal ice stress cpt
      INTEGER :: joc_isig22ne        !Internal ice stress cpt
      INTEGER :: joc_isig22se        !Internal ice stress cpt
      INTEGER :: joc_isig22sw        !Internal ice stress cpt
      INTEGER :: joc_isig22nw        !Internal ice stress cpt
      INTEGER :: joc_snown           ! Snow depth over ice for each cat
      INTEGER :: joc_iceconn         ! Ice concentration for each cat
      INTEGER :: joc_icedepn         ! Mean ice depth for each category

      ! Ocean ancillary
      INTEGER :: joc_taux                ! Zonal windstress
      INTEGER :: joc_tauy                ! Merid windstress
      INTEGER :: joc_wme                 ! Wind mixing energy
      INTEGER :: joc_w10                 ! 10m windspeed
      INTEGER :: joc_surfp               ! Surface pressure
      INTEGER :: joc_solar               ! Solar heating
      INTEGER :: joc_heat                ! Non-penetrative heat
      INTEGER :: joc_ple                 ! Precip less evap
      INTEGER :: joc_river               ! River outflow
      INTEGER :: joc_watop               ! Water optical prop
      INTEGER :: joc_ddepd1              ! Dust depn div 1 : ocean
      INTEGER :: joc_ddepd2              ! Dust depn div 2 : ocean
      INTEGER :: joc_ddepd3              ! Dust depn div 3 : ocean
      INTEGER :: joc_ddepd4              ! Dust depn div 4 : ocean
      INTEGER :: joc_ddepd5              ! Dust depn div 5 : ocean
      INTEGER :: joc_ddepd6              ! Dust depn div 6 : ocean
      INTEGER :: joc_atmco2              ! atmospheric co2 conc
      INTEGER :: joc_hep                 ! Depth of mid-ocean ridge srce
      INTEGER :: joc_he3                 ! Helium-3 concentration
      INTEGER :: joc_he4                 ! Helium-4 concentration

      ! Sea ice ancillary
      INTEGER :: joc_solice              ! Solar radn over ice
      INTEGER :: joc_snowrate            ! Snowfall
      INTEGER :: joc_sublim              ! Sublimation
      INTEGER :: joc_topmelt             ! Top melt from ice
      INTEGER :: joc_botmelt             ! Bottom melt from ice
      INTEGER :: joc_topmeltn            ! Top melt from ice for each cat
      INTEGER :: joc_botmeltn            ! Bottom melt from ice for each
                                         ! category

      ! Fluxes local over ice/leads fraction
      INTEGER :: joc_solar_l             ! Solar heating (local)
      INTEGER :: joc_heat_l              ! Non-penetrative heat (local)
      INTEGER :: joc_sublim_l            ! Sublimation (local)
      INTEGER :: joc_topmelt_l           ! Top melt from ice (local)
      INTEGER :: joc_botmelt_l           ! Bottom melt from ice (local)
      INTEGER :: joc_ple_l               ! Precip less evap (local)
      INTEGER :: joc_wme_l               ! Wind mixing energy (local)

      ! Ocean flux correction (ancillary)
      INTEGER :: joc_climsst             ! Reference surf temp
      INTEGER :: joc_climsal             ! Ref surf sal'ty
      INTEGER :: joc_climair             ! Reference air temp
      INTEGER :: joc_climicedep          ! Reference ice depth
      INTEGER :: joc_anom_heat           ! Heat flux correction
      INTEGER :: joc_anom_salt           ! Salinity flux corrn

      ! User ancillaries
      INTEGER :: jousr_anc1              ! User ancillary 1
      INTEGER :: jousr_anc2              ! User ancillary 2
      INTEGER :: jousr_anc3              ! User ancillary 3
      INTEGER :: jousr_anc4              ! User ancillary 4
      INTEGER :: jousr_anc5              ! User ancillary 5
      INTEGER :: jousr_anc6              ! User ancillary 6
      INTEGER :: jousr_anc7              ! User ancillary 7
      INTEGER :: jousr_anc8              ! User ancillary 8
      INTEGER :: jousr_anc9              ! User ancillary 9
      INTEGER :: jousr_anc10             ! User ancillary 10
      INTEGER :: jousr_mult1             ! multi-lev user ancil 1
      INTEGER :: jousr_mult2             ! multi-lev user ancil 2
      INTEGER :: jousr_mult3             ! multi-lev user ancil 3
      INTEGER :: jousr_mult4             ! multi-lev user ancil 4

      ! Ocean housekeeping
      INTEGER :: joc_index_comp          ! Compress array index
      INTEGER :: joc_index_exp           ! Expanded array index
      INTEGER :: joc_index_start         ! Rows and levels
      INTEGER :: joc_no_seapts           ! Number of comp pts
      INTEGER :: joc_no_segs             ! No of segs in comp

      ! Ocean analysis increments
      INTEGER :: joc_temp_incs           ! T analysis increments
      INTEGER :: joc_sal_incs            ! S analysis increments
      ! Ocean Temperature and Salinity bias fields
      INTEGER :: joc_t_bias
      INTEGER :: joc_s_bias

      ! Rigid lid pressure
      INTEGER :: joc_rlp                 ! Rigid lid pressure

      COMMON/CARGPT_OCEAN/                                              &
     &  joc_u, joc_v,joc_stream,                                        &
     &  joc_FWFLX,joc_FWFLXb,joc_tend,joc_pmepassb,joc_dusp,            &
     &  joc_eta,joc_etab,joc_ubt,joc_ubtbbt,                            &
     &  joc_vbt,joc_vbtbbt,joc_ubtbbc,joc_vbtbbc,                       &
     &  joc_mld,joc_athkdft,joc_snow,joc_mischt,joc_varchl1,joc_varchl2,&
     &  joc_htotoi, joc_salinc, joc_icy, joc_icecon, joc_icedep,        &
     &  joc_iceu,joc_icev,                                              &
     &  joc_isig11ne,                                                   &
     &  joc_isig11se,                                                   &
     &  joc_isig11sw,                                                   &
     &  joc_isig11nw,                                                   &
     &  joc_isig12ne,                                                   &
     &  joc_isig12se,                                                   &
     &  joc_isig12sw,                                                   &
     &  joc_isig12nw,                                                   &
     &  joc_isig22ne,                                                   &
     &  joc_isig22se,                                                   &
     &  joc_isig22sw,                                                   &
     &  joc_isig22nw,                                                   &
     &  joc_snown,joc_iceconn,joc_icedepn,joc_topmeltn,joc_botmeltn,    &
     &  joc_solar_l,joc_heat_l,joc_sublim_l,joc_topmelt_l,              &
     &  joc_botmelt_l,joc_ple_l,joc_wme_l,                              &
     &  joc_isx,joc_isy, joc_taux, joc_tauy,                            &
     &  joc_wme, joc_w10, joc_surfp, joc_solar, joc_heat,               &
     &  joc_ple, joc_river, joc_watop, joc_solice, joc_snowrate,        &
     &  joc_ddepd1, joc_ddepd2, joc_ddepd3, joc_ddepd4, joc_ddepd5,     &
     &  joc_ddepd6,                                                     &
     &  joc_hep, joc_he3, joc_he4,                                      &
     &  joc_sublim,                                                     &
     &  joc_atmco2,                                                     &
     &  joc_climsst, joc_climsal, joc_climair,                          &
     &  joc_climicedep, joc_anom_heat, joc_anom_salt,                   &
     &  joc_topmelt, joc_botmelt, joc_index_comp, joc_index_exp,        &
     &  joc_index_start, joc_no_seapts, joc_no_segs,                    &
     &  jousr_anc1, jousr_anc2, jousr_anc3, jousr_anc4, jousr_anc5,     &
     &  jousr_anc6, jousr_anc7, jousr_anc8, jousr_anc9, jousr_anc10,    &
     &  jousr_mult1, jousr_mult2, jousr_mult3, jousr_mult4,             &
     &  joc_temp_incs, joc_sal_incs, joc_rlp, joc_t_bias, joc_s_bias

#endif
! TYPPTRO end
