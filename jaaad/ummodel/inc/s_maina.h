!include file: s_maina.h
! Description:
!   Defines additional namelist data for the SCM.
!
! Current Code Owner: R Wong
!

#if defined(ATMOS)

      Integer :: row_length
      Integer :: rows
      Integer :: max_land_points

      Parameter (row_length = 1)
      Parameter (rows = 1)
      Parameter (max_land_points = row_length*rows)      

! Used for declarations
      Integer  :: max_no_ntiles
      Parameter (max_no_ntiles = 20)
      Integer  :: max_soil_temp_levs
      Parameter (max_soil_temp_levs = 10)
      Integer  :: max_soil_moist_levs
      Parameter (max_soil_moist_levs = 10)
      Integer  :: max_wet_levels
      Parameter (max_wet_levels = 100)
      Integer  :: max_tr_levels
      Parameter (max_tr_levels = 100)
      Integer  :: max_tr_vars
      Parameter (max_tr_vars = 30)
      Integer  :: max_ozone_levels
      Parameter (max_ozone_levels = 100)
      Integer  :: max_nfor
      Parameter (max_nfor = 500)
      Integer  :: max_nsprog
      Parameter (max_nsprog = 15)
! Soil carbon dimensions hardwired to 1 because Triffid disabled
      Integer  :: dim_cs1
      Parameter (dim_cs1 = 1)
      Integer  :: dim_cs2
      Parameter (dim_cs2 = 1)

!---------------------------------------------------------------------
!     &INDATA information
!---------------------------------------------------------------------

      Integer ::                                                        &
         soil_type(max_land_points)                                     &
      ,  soil_type_scm(land_points)                                     &
                                ! Soil type code 1 to 4
                                !  1 Ice
                                !  2 Fine
                                !  3 Medium
                                !  4 Coarse
      ,  salt_dim1, salt_dim2, salt_dim3
                                ! Dimensions of sea-salt aerosol arrays

      Real                                                              &
     &   gridbox_area(row_length, rows)                                 &
                                ! (km^2)
                                ! Global dimensions = the price of
                                ! fish (huge to minimise advective
     &  ,lat(row_length, rows)                                          &
                                          ! Lat. of gridpoint chosen
                                !  Read automatically from climate
                                !  dataset if STATS forcing chosen
     &  ,long(row_length, rows)
                                          ! Long. of gridpoint chosen
                                !  Read automatically from climate
                                !  dataset if STATS forcing chosen

      Logical                                                           &
     &  gather  ! true if running on the Cray

      Integer                                                           &
     &  year_init                                                       &
                                ! Initial year
     &  ,month_init                                                     &
                                ! Initial month
     &  ,day_init                                                       &
                                ! Initial day
     &  ,hour_init                                                      &
                                ! Initial hour
     &  ,min_init                                                       &
                                ! Initial minute
     &  ,sec_init                                                       &
                                ! Initial second
     &  ,tapeyear_init                                                  &
                                ! Initial year for tape input
     &  ,tapemonth_init                                                 &
                                ! Initial month for tape input
     &  ,tapeday_init                                                   &
                                ! Initial day for tape input
     &  ,tapehour_init                                                  &
                                ! Initial hour for tape input
     &  ,tapemin_init                                                   &
                                ! Initial minute for tape input
     &  ,tapesec_init           ! Initial second for tape input

!---------------------------------------------------------------------
!     &RUNDATA information
!---------------------------------------------------------------------

      CHARACTER*2 modelID       ! Forcing model identifier
      Character*8                                                       &
     &  exname_in                                                       &
                                ! Name of expt. to be read
                                !  from previous run stored
                                !  on tape up to 6 Characters
     &  ,exname_out             ! Name of expt. to be written
                                !  to tape up to 6 Characters
      Integer                                                           &
     &  change_clim                                                     &
                                ! No. of days between
                                !  changes of climatological data
                                !  (default=10)
     &, dump_step                                                       &
                                ! No. of timesteps for which mean dump
                                ! is required (default of 0 will cause
                                ! dump_step to be set such that time
                                ! between dumps is close to 1200s)
     &, dump_days(4)                                                    &
                                ! No. of days for which mean
                                !  dump is required (default=1)
     &, ndayin                                                          &
                                ! No. of days requested in run
     &, nminin                                                          &
                                ! No. of minutes requested in run
     &, nsecin                                                          &
                                ! No. of seconds requested in run
     &, resdump_days                                                    &
                                ! frequency of dumps for restart
     &, runno_in                                                        &
                                ! Number of run to be read
                                !  from previous run stored on tape
     &, runno_out                                                       &
                                ! Number of run to be written to tape
     &, min_trop_level                                                  &
                        ! first reference model level about 700hpa
     &, max_trop_level                                                  &
                        ! first reference model level about 50hpa
     &, ntml(row_length, rows)                                          &
                                ! top level of surface mixed layer
     &, nbdsc(row_length, rows)                                         &
                                ! Bottom level of decoupled SC layer
     &, ntdsc(row_length, rows)                                         &
                                ! Top level of decoupled SC layer
     &, ntrad                                                           &
                                ! No. of timesteps between
                                ! calls to radiation
     &, ntrad1
                                ! Timestep of 1st Radiation call

     Real :: albsoil(row_length*rows)
                                ! Soil albedo, User Specified
                                ! In future, should be land_points

     Real :: & ! Moved to here from scm_main because they are in
               ! namelists
     &  sum_eng_fluxes(row_length,rows)                                 &
                                         ! sum atmosphere fluxes
     &, sum_moist_flux(row_length,rows)                                
                                         ! sum moist fluxes
      Real                                                              &
     &  timestep                ! Model timestep for all physics
                                !  subroutines except radiation
      Real                                                              &
     &  co2start                                                        &
     &, co2end                                                          &
     &, co2rate                                                         &
     &, so2_em    ( row_length, rows )                                  &
                                         ! so2 surface emissions
     &, nh3_em    ( row_length, rows )                                  &
                                         ! nh3 surface emissions
     &, dms_em    ( row_length, rows )                                  &
                                         ! dms surface emissions
     &, soot_em   ( row_length, rows )                                  &
                                         ! surface soot emissions
     &, soot_hilem ( row_length, rows )                                 &
     &, bmass_em   ( row_length, rows )                                 &
                                         ! surface biomass emissions
     &, bmass_hilem ( row_length, rows )                                &
                                         ! high level biomass emiss
     &, ocff_em    ( row_length, rows )                                 &
                                         ! surface fossil-fuel OC emiss
     &, ocff_hilem ( row_length, rows )                                 &
                                         ! high level fossil-fuel OC em
     &, soot(row_length, rows)                                          &
                                         ! Snow soot
     &, cort(row_length, rows)                                          &
                                         ! Vertical correlation coeff.
                                         !  for temp.
     &, cord(row_length, rows)                                          &
                                ! Vertical correlation coeff for
                                ! dew pt. depression
     &, corvn(row_length, rows)                                         &
                                ! Vertical correlation coeff.
                                ! for velocity VN
     &, corw(row_length, rows)                                          &
                                ! Vertical correlation coeff.
                                ! for vertical velocity
     &  ,ozone(row_length, rows, max_ozone_levels)                      &
     &  ,cclwp(row_length, rows)                                        &
                                ! Condensed water path Kg m^-2
     &  ,orog(row_length, rows)                                         &
                                ! Orogrophy
     &  ,so2(row_length, rows,max_wet_levels)                           &
                                  ! Sulphur Cycle tracers for wet
     &  ,so4_aitken(row_length, rows,max_wet_levels)                    &
                                                     !scavenging.
     &  ,so4_accu(row_length, rows,max_wet_levels)                      &
                                                     !
     &  ,so4_diss(row_length, rows,max_wet_levels)                      &
                                                     !
     &  ,aerosol(row_length, rows,max_model_levels)                     &
                                ! Aerosol values ; only used if
                                !    l_murk=.true. ; default .false.
     &, aerosol_em(row_length,rows,max_model_levels)                    &
     &, nh3( row_length, rows, max_model_levels )                       &
     &, dms( row_length, rows, max_model_levels )                       &
     &, soot_new( row_length, rows, max_model_levels )                  &
     &, soot_aged( row_length, rows, max_model_levels )                 &
     &, soot_cld( row_length, rows, max_model_levels )                  &
     &, bmass_new( row_length, rows, max_model_levels )                 &
                                                         ! biomass
     &, bmass_aged( row_length, rows, max_model_levels )                &
                                                         ! aerosol
     &, bmass_cld( row_length, rows, max_model_levels )                 &
                                                         ! arrays
     &, ocff_new( row_length, rows, max_model_levels )                  &
                                                         ! fossil-fuel
     &, ocff_aged( row_length, rows, max_model_levels )                 &
                                                         ! OC aerosol
     &, ocff_cld( row_length, rows, max_model_levels )                  &
                                                         ! arrays
     &, co2( row_length, rows, max_model_levels )                       &
     &, zh(row_length, rows)                                            &
                                ! Height above surface of top
                                !  of boundary layer (m)

     &, dOLR_rts(row_length, rows)                                      &
                                         ! TOA - surface upward LW
     &, fland_ctile(row_length, rows)                                   &
                                      !Land fraction on land points.

     &, tstar_land(row_length, rows)                                    &
!                               ! Land mean sfc temperature (K)
     &, tstar_sea(row_length, rows)                                     &
!                               ! Open sea sfc temperature (K).
     &, tstar_sice(row_length, rows)                                    &
!                               ! Sea-ice sfc temperature (K).
     &, land_alb(row_length, rows)                                      &
!                               ! Mean land albedo.
     &, sice_alb(row_length, rows)                                      &
!                               ! Sea-ice albedo.
     &, co2_emits(row_length, rows )                                    &
     &, co2flux(row_length, rows )

!---------------------------------------------------------------------
!      Tracer Flux Information - kdcorbin, 04/10
!---------------------------------------------------------------------

     Real :: tracer_flux1(row_length,rows),  tracer_flux2(row_length,rows), &
	&  tracer_flux3(row_length,rows),  tracer_flux4(row_length,rows), &
	&  tracer_flux5(row_length,rows),  tracer_flux6(row_length,rows), &
	&  tracer_flux7(row_length,rows),  tracer_flux8(row_length,rows), &
	&  tracer_flux9(row_length,rows),  tracer_flux10(row_length,rows),&
	&  tracer_flux11(row_length,rows), tracer_flux12(row_length,rows),&
	&  tracer_flux13(row_length,rows), tracer_flux14(row_length,rows),&
	&  tracer_flux15(row_length,rows), tracer_flux16(row_length,rows),&
	&  tracer_flux17(row_length,rows), tracer_flux18(row_length,rows),&
	&  tracer_flux19(row_length,rows), tracer_flux20(row_length,rows)

!---------------------------------------------------------------------
!     &LOGIC information
!---------------------------------------------------------------------

      Logical                                                           &
     &  test                                                            &
                                ! T if detailed sub-timestep
                                !   diagnostics required
     &, altdat                                                          &
                                ! T if alternative initial profiles
                                !   of T,Q,U and V are to be input
     &, ancyc                                                           & 
                                ! T if annual cycle req'd
                                !   (ie. radiation input then
                                !   varies throughout year)
     &, geoforce                                                        &
                                ! T if geostrophic forcing.
     &, geoinit                                                         &
                                ! T if initialising dump to
                                !   geostrophic.
     &, grafdump_day                                                    &
                                ! T if graphical dump of mean dai
                                !   values required
     &, grafdump_days                                                   &
                                ! T if graphical dump of mean val
                                !   over DUMP_DAYS required
     &, grafdump_step                                                   &
                                ! T if graphical dump of mean val
                                !   required each DUMP_STEP
     &, noforce                                                         &
                                ! T if no large-scale forcing
                                !   is required
     &, obs                                                             &
                                ! T if observational
                                !   large-scale forcing used
     &, prindump_day                                                    &
                                ! T if printout of mean daily
                                !   dump required
     &, prindump_days                                                   &
                                ! T if printout of mean dump
                                !   over DUMP_DAYS required
     &, prindump_obs                                                    &
                                ! T if printout of observational
                                !   diagnostics required every OBS_
                                !   timesteps
     &, prindump_step                                                   &
                                ! T if printout of mean dump
                                !   required each DUMP_STEP
     &, prinstat                                                        &
                                ! T if printout of stats forcing
                                !   required every timestep
     &, radcloud_fixed                                                  &
                                ! T if cloud required fixed for
                                !   radiation
     &, stats                                                           &
                                ! T if statistical
                                !   large-scale forcing used
     &, land_sea_mask(row_length, rows)                                 &
                                        ! T for a land point
     &, land_ice_mask(row_length, rows)                                 &
                                        ! T for an ice land point
     &, soil_mask(row_length, rows)                                     &
                                     ! T for a soil point
     &, tapein                                                          &
                                ! T if initial data is to be read
                                !   from previous run stored on
                                !   tape
     &, tapeout                                                         &
                                ! T if restart information
                                !   plus diagnostic output to be
                                !   stored on tape routines
     &, cumulus (row_length, rows)                                      &
                                   ! bl convection flag
     &, local_time                                                      &
                                ! T if diagnostics required
                                !   for local time rather than GMT
     &, l_do_t                                                          &
                                ! keep t increments
     &, l_do_inc_vels                                                   &
                                ! keep velocity increments
     &, l_do_inc_q                                                      &
                                ! keep moisture increments
     &, L_flux_bc                                                       &
                                ! T if prescribed surface fluxes to be
                                ! used, otherwise T-star forcing used.
     &, relaxT_eurocs                                                   &
                                !eurocs forcing relax to T profile 
                                !switched on
     &, relaxq_eurocs                                                   &
                                !eurocs forcing relax to q profile 
                                !switched on
     &, relaxuv_eurocs                                                  
                                !eurocs forcing relax to u,v profiles 
                                !switched on

!-----------------------------------------------------------------------------
!     &INMOSES information
!-----------------------------------------------------------------------------



      Integer :: smi_opt  ! Option to define method of initialisation of soil
                          ! moisture content:    0 : Use smcli
                          !                      1 : Use fsmc
                          !                      2 : Use sth

      ! In Future, Surface variables should be set-up to be
      ! land-points to be consistent with full model routines.
      Real ::                          &
        gs   (row_length,rows)          ! Stomatal conductance

      Real :: & ! Soil Moisture Stress
        fsmc (max_land_points)                          &
      , fsmc_scm (land_points)    

      Real :: & ! Fractions of surface types
        frac_typ (max_land_points*ntype)                &
      , frac_typ_scm (land_points,ntype) 

      Real :: & ! Initial SMCL profile(Kg/m^2)
        smcli     (max_land_points*max_soil_moist_levs) &
      , smcli_scm (land_points, sm_levels)

      Real :: & ! Total soil moisture in layers as fraction of saturation
        sth       (max_land_points*max_soil_moist_levs) & 
      , sth_scm   (land_points, sm_levels)



      Real ::                                  &
        canht    (max_land_points, npft)       & ! Canopy height(m)
      , lai      (max_land_points, npft)         ! Leaf area index


      Real ::                                 &
        z0_tile    (max_land_points, max_no_ntiles)  &! Tile roughness lengths (m).
      , rgrain     (max_land_points, max_no_ntiles)  &! Snow grain size (microns)
      , infil_tile (max_land_points, max_no_ntiles)  &! Max. surface infiltration
      , snow_tile  (max_land_points, max_no_ntiles)  &! Lying snow on tiles (MOSES 2)
      , tstar_tile (max_land_points, max_no_ntiles)  &! Surface tile temperature
      , catch      (max_land_points, max_no_ntiles)  &! Surf/canopy water capacity
                                               ! (snow-free land tiles)
                                               ! (Kg/m^2)
      , canopy     (max_land_points, max_no_ntiles)   ! Surf/canopy water
                                               ! (snow-free land tiles)
                                               ! (Kg/m^2)


      ! 2A Triffid variables
      !=====================
      Real ::                                  &
        frac_disturb    (max_land_points)      &! Fraction of gridbox in which
                                                ! vegetation is disturbed.
      , npp_ft_acc      (max_land_points)          &! Acc. NPP_FT
      , resp_w_ft_acc   (max_land_points, npft)    &! Acc. RESP_W_FT
      , resp_s_acc      (max_land_points, dim_cs1) &! Acc. RESP_S
      , cs              (max_land_points, dim_cs1) &! Soil carbon (kg C/m2)
      , g_leaf_phen_acc (max_land_points, npft)    &! Acc. g_leaf
      , g_leaf_acc      (max_land_points, npft)     ! Acc. leaf turnover rate
                                                    ! with phenology

      Real :: CO2_3D(row_length, rows)   ! IN 3D CO2 field if required.
      Real :: LW_down(row_length, rows)  ! Surface downward LW


!---------------------------------------------------------------------
!     &INGEOFOR information
!---------------------------------------------------------------------

      Real                                                              &
     &  ug(row_length, rows)                                            &
                                ! Geostrophic U velocity (m s^-1)
     &, vg(row_length, rows)    ! Geostrophic V velocity (m s^-1)


!---------------------------------------------------------------------
!     &INOBSFOR information
!---------------------------------------------------------------------
      Logical                                                           &
     &  L_windrlx, L_vertadv
      Real                                                              &
     &  tau_rlx

      Real                                                              &
     & tau_T_eurocs
                                !relaxation timescale for T for EUROCS 
                                !forcing
      Real                                                              &
     & tau_q_eurocs
                                !relaxation timescale for q for EUROCS 
                                !forcing
      Real                                                              &
     & tau_uv_eurocs
                                !relaxation timescale for u,v for EUROCS 
                                !forcing


      Integer                                                           &
     &  ichgf                                                           &
                                ! No. of timesteps between change in
                                !  observational forcing
     &  ,ilscnt                 ! Counts for observational forcing


     !----------------------------------------------------------------
     ! Surface forcings
     !----------------------------------------------------------------
      Real                                                           :: &
     &  flux_e_scm       (row_length,rows,nfor)                         & 
                                           ! Evaporation fluxes.
     &, flux_h_scm       (row_length,rows,nfor)                         &
                                           ! Heat fluxes.
     &, tstar_forcing_scm(row_length,rows,nfor)     
                                           ! Sea Surface temperature


     !----------------------------------------------------------------
     ! Increments due to large-scale horizontal and vertical advection
     !----------------------------------------------------------------

      Real                                                           :: &
     &  t_inc_scm (row_length,rows,nfor,    model_levels)               &
                                           ! Temperature increment
                                           ! (K s^-1 day^-1)
     &, u_inc_scm (row_length,rows,nfor,    model_levels)               &
                                           ! Absolute Zonal wind
                                           ! (m/s)
     &, v_inc_scm (row_length,rows,nfor,    model_levels)               & 
                                           ! Absolute Meridional wind
                                           ! (m/s)
     &, w_inc_scm (row_length,rows,nfor,    model_levels)               &
                                           ! Absolute Vertical wind
                                           ! (m/s)
     &, q_star_scm(row_length,rows,nfor,wet_model_levels) 
                                           ! Specific humidity increment
                                           ! (Kg Kg^-1 s^-1 day^-1)




!---------------------------------------------------------------------
!     &INPROF information
!---------------------------------------------------------------------

      Integer                                                           &
     &  iccbi(row_length, rows)                                         &
                                ! Convective cloud base  and top
     &  ,iccti(row_length, rows)!  (model levels)


      Real :: & ! Initial deep soil temps. (K)
        t_deep_soili     (max_land_points*max_soil_temp_levs)           &
      , t_deep_soili_scm (land_points, st_levels)

      Real :: & ! Initial canopy water content (Kg/m2)
        canopy_gbi     (max_land_points)                                &
      , canopy_gbi_scm (land_points)

      Real :: & ! Initial soil moisture content (Kg/m2)
        smci     (max_land_points)                                      &
      , smci_scm (land_points)

      Real                                                              &
     &   ccai(row_length, rows)                                         &
                                  ! Convective cloud amnt.
                                !  (decimal fraction)
     &  ,snodepi(row_length, rows)                                      &
                                   ! Initial snow depth (Kg m^-2)
     &  ,tstari(row_length, rows)                                       &
                                 ! Initial surface temp. (K)
     &  ,z0mseai(row_length, rows)                                      &
                                ! Initial sea surface roughness
                                ! length (m)
     &  ,z0m_scm(row_length, rows)                                      &
                                ! Fixed value sea surface roughness
                                ! length for momentum (m)
     &  ,z0h_scm(row_length, rows)                                      &
                                ! Fixed value sea surface roughness
                                ! length for heat (m)
     &  ,sil_orog_land(row_length, rows)                                &
                                ! Silhouette area of unresolved
                                !  orography per unit horizontal
                                !  area on land points only.
     &  ,ho2r2_orog(row_length, rows)                                   &
                                ! Standard Deviation of orography
                                !  equivalent to peak to trough
                                !  height of unresolved orography
                                !  divided by 2SQRT(2) on land
                                !  points only (m)
     &  ,ice_fract(row_length, rows)                                    &
                                ! Fraction of grid box covered by
                                !  sea ice(decimal fraction)
     &  ,di(row_length, rows)                                           &
                                ! Equivalent thickness of sea-ice (m)
     &  ,u_0(row_length, rows)                                          &
                                ! Westerly & easterly component of
     &  ,v_0(row_length, rows)                                          &
                                !  surface current (metres per second)

     &  ,tracer(row_length, rows,max_tr_levels,max_tr_vars)             
                                ! Model tracer fields (Kg Kg^-1)
                                ! Pressure on rho levels
      !-----------------------------------------------------------------
      ! Initial Profiles
      !-----------------------------------------------------------------
      Real                                                           :: &
     &  qi_scm     (row_length,rows,wet_model_levels)                   &
                                ! Specific humidity, (Kg/Kg)
     &, theta_scm  (row_length,rows,    model_levels)                   &
                                ! Potential temperature, (K)
     &, ti_scm     (row_length,rows,    model_levels)                   & 
                                ! Temperature, (K)
     &, ui_scm     (row_length,rows,    model_levels)                   &
                                ! Absolute zonal wind, (m/s)
     &, vi_scm     (row_length,rows,    model_levels)                   &
                                ! Absolute meridional wind, (m/s)
     &, wi_scm     (row_length,rows,    model_levels+1)                 &
                                ! Vertical velocity (m/s)
     &, p_in_scm   (row_length,rows,    model_levels+1)                 & 
                                ! Pressure (Pa)
     &, w_advi_scm (row_length,rows,  0:model_levels)    

      !----------------------------------------------------------------
      ! Dummy arrays to ensure above forcing data is read in correctly.
      ! variable names have been kept the same so user namelists are 
      ! unaffected.
      !----------------------------------------------------------------

      Real                                                           :: &
     &  qi    (row_length*rows* max_wet_levels     ) = 0.0              & 
     &, theta (row_length*rows* max_model_levels   ) = 0.0              &
     &, ui    (row_length*rows* max_model_levels   ) = 0.0              &
     &, vi    (row_length*rows* max_model_levels   ) = 0.0              &
     &, wi    (row_length*rows*(max_model_levels+1)) = 0.0              &
     &, p_in  (row_length*rows*(max_model_levels+1)) = 0.0              &
     &, w_advi(row_length*rows*(max_model_levels+1)) = 0.0          

      Real                                                           :: &     
     &  q_star(row_length*rows*max_nfor* max_wet_levels   ) = 0.0       &
     &, t_inc (row_length*rows*max_nfor* max_model_levels ) = 0.0       &
     &, u_inc (row_length*rows*max_nfor* max_model_levels ) = 0.0       &
     &, v_inc (row_length*rows*max_nfor* max_model_levels ) = 0.0       &
     &, w_inc (row_length*rows*max_nfor*(max_model_levels)) = 0.0

      Real                                                           :: &
     &  flux_e       (row_length*rows*max_nfor) = 0.0                   &
     &, flux_h       (row_length*rows*max_nfor) = 0.0                   &
     &, tstar_forcing(row_length*rows*max_nfor) = 0.0

!---------------------------------------------------------------------
!     Variables for diagnostic output
!---------------------------------------------------------------------

      Integer                                                           &
     &  nout(19)                ! Output units to receive printout of
                                ! initial data by PRINT_INITDATA

!---------------------------------------------------------------------
!     Cloud variable values if clouds to be fixed for radiation
!---------------------------------------------------------------------

      Integer                                                           &
     &  iccb_rad(row_length, rows)                                      &
                                ! Model level of base of convective
                                !  cloud
     &  ,icct_rad(row_length, rows)
                                ! Model level of top of convective
                                !  cloud
      Real                                                              &
     &  cca_rad(row_length, rows)                                       &
                                ! Convective cloud amount (fraction)
     &  ,layer_cloud_rad(row_length, rows,max_wet_levels)               &
                                ! Layer cloud amount (fraction)
     &  ,qcl_rad(row_length, rows,max_wet_levels)                       &
                                ! Total cloud water and ice content
                                !  over cloud(Kg Kg^-1)
     &  ,qcf_rad(row_length, rows,max_wet_levels)                       &
                                ! Tet to zero as user will usually
                                !  input combined cloud water and
                                !  ice content over cloud(Kg Kg^-1)
                                !  in QCL_RAD.
     &  ,ccwpin_rad(row_length, rows)
                                ! Convective water path over cloud
                                !  only  (Kg m^-2).
!---------------------------------------------------------------------
!     conv_mode determines actions for convection scheme
!     0 = Run normally
!     1 = Run for diagnostics every
!     radiation timestep but save dump state
!     (except CCA,ICCB and ICCT)
!     2 = Don't run
      integer conv_mode

!---------------------------------------------------------------------
!     &DIAGS is defined in routine SETUP_DIAGS
!---------------------------------------------------------------------

#endif
