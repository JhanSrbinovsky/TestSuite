#if defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!     SUBROUTINE RUN_INIT---------------------------------------------
!
!     Purpose: Called by SCMMAIN (Single Column Model main routine) to
!     Do the initialisations (previously done in SCMMAIN).
!
!     Code Description:
!     Language - FORTRAN 77
!
!     Author: C. Bunton
!
!     Documentation: Single Column Model Guide - J. Lean
!=====================================================================
!     OPTIONS TO SET INITIAL PROFILES
!=====================================================================
! (i)   Observational large scale forcing (OBS=TRUE of namelist LOGIC)
!         Initial data is then from namelist INPROF
! (ii)  Statistical large scale forcing (STATS=TRUE of namelist LOGIC)
!         Initial data can either be derived from climate datasets
!         using subroutine INITSTAT or set from namelist
!         INPROF (set ALTDAT=TRUE in namelist LOGIC)
! (iii) No large-scale forcing initial data is set fron namelist
!         INPROF
! (iv)  Continuation from previous run stored on tape
!         (Set TAPEIN=TRUE in namelist LOGIC).  All other initial data
!         is overwritten
!=====================================================================
!
      Subroutine RUN_INIT(                                              &
!     IN leading dimensions of arrays
     &  row_length, rows, model_levels, nwet, land_pts                  &
     &  ,nfor, nbl_levs, nsoilt_levs, nsoilm_levs, ntiles, ntrop        &
!     IN dimension of dump array.
     &  ,nprimvars, stats, obs, prindump_obs, noforce, altdat           &
     &  ,land_mask, tapein, tapeout                                     &
     &  ,l_climat_aerosol, l_use_dust, l_use_sulpc_direct, ltimer       &
     &  ,l_ch4, l_n2o, l_cfc11, l_cfc12                                 &
     &  ,l_cfc113, l_hcfc22, l_hfc125, l_hfc134a                        &
     &  ,l_o2                                                           &
     &  ,l_use_soot_direct, l_use_bmass_direct, l_use_ocff_direct       &
     &  ,l_use_seasalt_direct, l_murk_rad, l_use_aod, l_use_biogenic    &
     &  ,l_use_arclbiom, l_use_arclblck, l_use_arclsslt, l_use_arclsulp &
     &  ,l_use_arcldust, l_use_arclocff, l_use_arcldlta                 &
     &  ,smi_opt, smcli, fsmc, sth                                      &
     &  ,geoforce, geoinit, ug, vg                                      &
     &  ,year_init, dayno_init, lcal360, ichgf, timestep, ndayin        &
     &  ,resdump_days, soil_type                                        &
     &  ,p_in, smci, canopy_gbi, snodepi, tstari, t_deep_soili          &
     &  ,z0mseai, ui, vi, wi, ti, qi, w_advi, ccai, iccbi, iccti        &
     &  ,time_init, tapeday_init                                        &
     &  ,exname_in, exname_out, runno_in, runno_out, theta, theta_star  &
     &  ,u, v, w, t, q, w_adv, flux_h, flux_e, uls, vls, wls, tls       &
     &  ,tstar_forcing, qls                                             &
     &  ,exner_rho_levels, exner_theta_levels, p                        &
     &  ,p_theta_levels, r_theta_levels, r_rho_levels, rho              &
     &  ,eta_theta_levels, eta_rho_levels, orog, height_gen_method      &
     &  ,ch_flux_h, ch_flux_e, ch_uls, ch_vls, ch_wls, ch_tls           &
     &  ,ch_qls, dap1, dap2, dap3, dab1, dab2, dab3, deltap ,pstar      &
     &  ,smc, smcl, canopy_gb, snodep, tstar, tsi, t_deep_soil          &
     &  ,sthu, sthf, frac_typ, canht, lai, catch, infil_tile            &
     &  ,z0_tile, catch_snow, can_model                                 &
     &  ,z0msea, zh, cca, n_cca_levels, iccb, icct, layer_cloud, qcf    &
     &  ,qcl, dayno_wint, alfada, alfadb, atime, btime, lat, long       &
     &  ,dbara, dbarb, dgrada, dgradb, pa, pb, rp, rp_theta             &
     &  ,tbara, tbarb, tgrada, tgradb, tsda, tsdb                       &
     &  ,vnbara, vnbarb, vnsda, vnsdb                                   &
     &  ,vpbara, vpbarb, wbara, wbarb, wsda, wsdb                       &
     &  ,iv, ntab, iy, idum, iseed, resdump                             &
     &  ,rhcrit                                                         &
     &  , b_exp, hcap, hcon, satcon, sathh, v_sat, v_wilt, v_crit)

      Implicit none

#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
!
! If new version of radiation code is required then we
! need to use the modules for improved time-stepping
!
      Use CONTROL_STRUC
      Use SW_CONTROL_STRUCT
      Use LW_CONTROL_STRUCT
      Use SPEC_SW_LW
#endif
      Integer                                                           &
     &  Row_length                                                      &
                                ! IN leading x dimension of SCM arrays.
     &  ,Rows                                                           &
                                ! IN leading y dimension of SCM arrays.
     &  ,model_levels                                                   &
                                       ! IN no of levels.
     &  ,nwet                                                           &
                                ! IN no of model levels in which Q is
                                !  set.
     &  ,land_pts                                                       &
                                ! IN Number of land points to be
                                !    processed.
     &  ,nfor                                                           &
                                ! IN Number terms for observational
                                !  forcing
     &  ,nbl_levs                                                       &
                                ! IN Number of Boundary layer levels
     &  ,nsoilt_levs                                                    &
                                ! IN Number of soil temperature
                                !  levels
     &  ,nsoilm_levs                                                    &
                                ! IN Number of soil moisture levels

     &  ,ntiles                                                         &
                                ! IN Number of surface tiles
     &  ,ntrop                                                          &
                                ! IN Max number of levels in the
                                !  troposphere
     &  ,nprimvars                                                      &
                                ! IN minimum no. of variables
                                !  required to restart from a dump.
     &  ,n_cca_levels                                                   &
                                ! IN No of levels for cca
     &  ,height_gen_method

      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='S_RUNINI')

      Real                                                              &
     &  eta_rho_levels(model_levels)                                    &
     &, eta_theta_levels(0:model_levels)                                &
     &, r_rho_levels (row_length,rows, model_levels)                    &
     &, r_theta_levels (row_length,rows, 0:model_levels)                &
     &, z_top_of_model                                                  &
     &, orog(row_length, rows) !Height above mean sea level in metres
                               !used to derive r_theta_levels and
                               !r_rho_level
!
!     Comdecks
!
!     SCM specific :

! Surface scheme.
! Soil Parameters for land surface
#include "s_soilpr.h"
!     Others :
#include "c_densty.h"
#include "cmaxsize.h"
! Radiation LW nad SW spectral tables
#include "cconsts.h"
! Actual soil layer thicknesses(MOSES)
#include "soil_thick.h"
! MOSES II Size Parameters
#include "nstypes.h"
! Various soil parameters
#include "c_soilh.h"
! Comdecks for coriolis parameter
#include "c_omega.h"
#include "c_pi.h"

#include "c_mdi.h"



!
!---------------------------------------------------------------------
!      Primary Model Variables plus T (UMDP No1)
!---------------------------------------------------------------------
!
      Integer                                                           &
     &  iccb(row_length, rows)                                          &
                                    ! Convective cloud base and top
     &  ,icct(row_length, rows)     !  at levels 1 to model_levels


      Real                                                              &
     &   canopy_gb(land_pts)                                            &
                                    ! Canopy water content (kg/m2)

     &  ,cca(row_length, rows, n_cca_levels)                            &
                                    ! Convective cloud amount
     &  ,pstar(row_length, rows)                                        & 
                                    ! Pressure at earth's surface
                                !  (Pa not hPa)
     &  ,q(row_length, rows,nwet)                                       &
                                    ! Specific humidity (Kg Kg^-1)
     &  ,qcf(row_length, rows,nwet)                                     &
                                    ! Cloud ice content (Kg Kg^-1)
     &  ,qcl(row_length, rows,nwet)                                     &
                                    ! Cloud water content(Kg Kg^-1)
     &  ,rccb(row_length, rows)                                         &
                                    ! Convective cloud base and top
     &  ,rcct(row_length, rows)                                         &
                                    !  at levels 1 to model_levels
                                !  real values for DUMP purposes
     &  ,smc(land_pts)                                                  &
                                    ! Soil moisture content(Kg m^-2)
     &  ,smcl(land_pts, nsoilm_levs)                                    &
                                ! Soil moisture content in layers
                                !  (Kg m^-2)
     &  ,sthf(land_pts, nsoilm_levs)                                    &
                                ! INOUT Frozen soil moisture
                                !  content of each layer as a
                                !  fraction of saturation.
     &  ,sthu(land_pts, nsoilm_levs)                                    &
                                ! INOUT Unfrozen soil moisture
                                !  content of each layer as a fraction
                                !  of saturation. (Kg m^-2)
     &  ,snodep(row_length, rows)                                       &
                                    ! Snow depth (Kg m^-2)
     &  ,t(row_length, rows,model_levels)                               &
                                           ! Temperature(K)
     &  ,t_deep_soil(land_pts, nsoilt_levs)                             &
                                ! Deep soil temperatures (K)
                                !  top level not included,=surface
     &  ,theta(row_length, rows,model_levels)                           &
                                    ! Potential temperature (K)
     &  ,theta_star(row_length, rows,model_levels)                      &
                                    ! Potential temp. increments(K)
     &  ,tsi(row_length, rows)                                          &
                                    ! Temperature of sea-ice
     &  ,tstar(row_length, rows)                                        &
                                    ! Surface temperature (K)
     &  ,u(row_length, rows,model_levels)                               &
     &  ,v(row_length, rows,model_levels)                               &
                                    ! Zonal,Meridional wind (m s^-1)
     &  ,w(row_length, rows,0:model_levels)                             &
     &  ,w_adv(row_length, rows,0:model_levels)                         &
     &  ,z0msea(row_length, rows)                                       &
                                    ! Sea surface roughness length
     &  ,zh(row_length, rows)                                           &
                                    ! Height above surface of top
                                ! of boundary layer (m)
     &  ,exner_rho_levels( row_length, rows, model_levels+1)            &
     &  ,exner_theta_levels( row_length, rows, model_levels)            &
     &  ,deltap(row_length, rows, model_levels)                         &
                                    ! Layer Thickness
     &  ,p(row_length, rows, model_levels+1)                            &
                                             ! Pressure on rho levels
     &  ,p_in(row_length, rows, model_levels+1)                         &
                                               ! Pressure on rho levels
     &  ,p_theta_levels(row_length, rows, model_levels)                 &
                                    !Pressure on theta
                                    ! levels
     &  ,rp(row_length, rows, model_levels+1)                           &
                                    ! Reciprocol pressure on rho levels
     &  ,rp_theta(row_length, rows, model_levels)                       &
                                    ! Reciprocol pressure on theta
                                    ! levels
     &  ,rho(row_length, rows, model_levels)                            &
     &  ,layer_cloud(row_length, rows,nwet)
                                ! layer cloud amount (decimal
                                !  fraction)

      Integer                                                           &
     &  year_init                                                       &
                                ! IN Initial year
     &  ,dayno_init                                                     &
                                ! IN Initial day in year
     &  ,tapeday_init           ! IN Initial day for tape input
                                !    ie last day on tape + 1

      Logical                                                           &
     &  lcal360                 ! IN ? 360 days year ?


!
!---------------------------------------------------------------------
!     Tape information
!---------------------------------------------------------------------
!
      Integer                                                           &
     &  runno                                                           &
                                ! Run no. of expt on tape
     &  ,tapeday                                                        &
                                ! Tape year day
     &  ,tapedump_no            ! No. of dumps on tape
      Character*8                                                       &
     &  exname                  ! Name of expt. on tape
      Real                                                              &
     &  resdump(row_length, rows,nprimvars)
                                ! DUMP array of restart variables
!
      Real                                                              &
     &   lat(row_length, rows)                                          &
                                ! Lat. of gridpoint chosen
                                !  Read automatically from climate
                                !  dataset if STATS forcing chosen
     &  ,long(row_length, rows)                                         &
                                    ! Long. of gridpoint chosen
                                ! Read automatically from climate
                                !  dataset if STATS forcing chosen
     &  ,time_init                                                      &
                                ! Initial time in seconds
     &  ,andayy                                                         &
                                ! LOC No. of days in 1 year.
                                !  (for one year effects)
     &  ,flux_e(row_length, rows,nfor)                                  &
     &  ,flux_h(row_length, rows,nfor)                                  &
                                          !
     &  ,ch_flux_e(row_length, rows,nfor-1)                             &
                                    ! Change per sec in FLUX_E,FLUX_H
     &  ,ch_flux_h(row_length, rows,nfor-1)                             &
     &  ,tls(row_length, rows,nfor,model_levels)                        &
     &  ,tstar_forcing(row_length, rows,nfor)                           &
                                ! Temp increment due to large-scale
                                !  horizontal and vertical advection
                                !  (K s^-1 day^-1)
     &  ,ch_tls(row_length, rows,nfor-1,model_levels)                   &
                                ! Change per sec in Temp increment
     &  ,qls(row_length, rows,nfor,nwet)                                &
                                    ! Specific humidity increment
                                !  due to large-scale horizontal
                                !  and vertical advection
                                !  (Kg Kg^-1 s^-1 day^-1)
     &  ,ch_qls(row_length, rows,nfor-1,nwet)                           &
                                  ! Change per sec in Specific humidity
     &  ,uls(row_length, rows,nfor,model_levels)                        &
                                    ! Zonal and meridional wind
     &  ,vls(row_length, rows,nfor,model_levels)                        &
                                    !  increment due to large-scale
                                !  horizontal and vertical
                                !  advection (m s^-1 day^-1)
     &  ,wls(row_length, rows,nfor,model_levels)                        &
     &  ,ch_uls(row_length, rows,nfor-1,model_levels)                   &
                                ! Change per sec in Zonal and merid
     &  ,ch_vls(row_length, rows,nfor-1,model_levels)                   &
                                    !  wind increm.
     &  ,ch_wls(row_length, rows,nfor-1,model_levels)
!
!
!     &INDATA
!
      Integer                                                           &
     &  soil_type(land_pts)
                                ! Soil type code:
                                !  1 Ice
                                !  2 Fine
                                !  3 Medium
                                !  4 Coarse
!
!     &INGEOFOR
!
      Real                                                              &
     &  UG(row_length, rows)                                            &
                                    ! Geostrophic U velocity (m s^-1)
     &  ,VG(row_length, rows)       ! Geostrophic V velocity (m s^-1)
!     &   ,QG(nwet)             ! Initial Moisture profile for
                                ! Geo. force
!
!     &INPROF
!
      Integer                                                           &
     &  iccbi(row_length, rows)                                         &
                                    ! Convective cloud base
     &  ,iccti(row_length, rows)    !  and top (model levels)

      Real                                                              &
     &   canopy_gbi(land_pts)                                           &
                                    ! Initial canopy water content
                                    ! (kg/m^2)
     &  ,ccai(row_length, rows)                                         &
                                    ! Convective cloud amnt.
     &  ,qi(row_length, rows,nwet)                                      &
                                    ! Initial specific humidity
                                !  (Kg Kg^-1)
     &  ,smci(land_pts)                                                 &
                                    ! Initial soil moisture content
                                !  (Kg m^-2)
     &  ,snodepi(row_length, rows)                                      &
                                    ! Initial snow depth (Kg m^-2)
     &  ,t_deep_soili(land_pts, nsoilt_levs)                            &
                                ! Initial deep soil temps. (K)
     &  ,ti(row_length, rows,model_levels)                              &
                                           ! Initial temp. profile  (K)
     &  ,tstari(row_length, rows)                                       &
                                    ! Initial surface temp. (K)
     &  ,ui(row_length, rows,model_levels)                              &
                                           ! Initial zonal and merid.
     &  ,vi(row_length, rows,model_levels)                              &
                                           !  wind comps. (m s^-1)
     &  ,wi(row_length, rows,0:model_levels)                            &
                                             !
     &  ,w_advi(row_length, rows,0:model_levels)                        &
                                                 !
     &  ,z0mseai(row_length, rows)  ! Initial sea surface roughness

!
!      &INMOSES
!
      Real                                                              &
        fsmc(land_pts)                                                  &
                                ! Soil moisture stress to initialise
                                ! SMCL
      , smcli(land_pts, nsoilm_levs)                                    &
                                ! Initial values for SMCL (Kg m^-2)
      , sth(land_pts, nsoilm_levs)
                                ! Total soil moisture in layers as a
                                !  fraction of saturation
      Integer :: smi_opt

      Integer, Parameter :: smc_spec = 0 ! Initial smcl specifed
                                         ! via namelist 
      Integer, Parameter :: smc_fsmc = 1 ! Initial smcl calculated
                                         ! via fsmc
      Integer, Parameter :: smc_sth  = 2 ! Initial smcl calculated
                                         ! via sth

      INTEGER                                                           &
     & CAN_MODEL

      REAL ::                                &
        frac_typ   (land_pts, ntype)         & ! IN Fractions of surface
                                               !    types.
      , canht      (row_length*rows, npft)   & ! IN Canopy Height (m)
      , lai        (row_length*rows, npft)   & ! IN Leaf Area Index
      , catch      (row_length*rows, ntiles) & ! OUT Surface/canopy water
                                               !     capacity of snow-free
                                               !     land tiles (kg/m2).
      , infil_tile (row_length*rows, ntiles) & ! OUT Maximum surface
                                               !     infiltration rate for
                                               !     each tile (kg/m2/s).
      , z0_tile    (row_length*rows, ntiles) & ! OUT Roughness length for
                                               !     each tile (m).
      , catch_snow (row_length*rows)           ! OUT Snow capacity for NLT
                                               !     tile (kg/m2).
!
!     &LOGIC
!
      Logical                                                           &
     &  altdat                                                          &
                                  ! T if alternative initial profiles
                                !  T,Q,U and V are to be input
     &  ,geoforce                                                       &
                                ! T if geostrophic forcing.
     &  ,geoinit                                                        &
                                ! T if initialising dump to
                                !  geostrophic.
     &  ,land_mask(row_length, rows)                                    &
                                    ! T for a land point
     &  ,noforce                                                        &
                                ! T if no large-scale forcing
                                !  is required
     &  ,obs                                                            &
                                ! T if observational
                                !  large-scale forcing used
     &  ,prindump_obs                                                   &
                                ! T if printout of observational
                                !  diagnostics required every
                                !  OBS_PRINT timesteps
     &  ,stats                                                          &
                                ! T if statistical large-scale
                                !  forcing used
     &  ,tapein                                                         &
                                ! T if initial data is to be read
                                !  from previous run stored on tape
     &  ,tapeout                                                        &
                                ! T if restart information plus
                                !  diagnostic output to be stored on
                                !  tape
     &  ,l_climat_aerosol                                               &
     &  ,l_use_sulpc_direct                                             &
     &  ,ltimer                                                         &
     &  ,l_ch4, l_n2o, l_cfc11, l_cfc12                                 &
     &  ,l_cfc113, l_hcfc22, l_hfc125, l_hfc134a                        &
     &  ,l_o2                                                           &
     &  ,l_use_soot_direct                                              &
                                ! Flag to use sulphur cycle for
                                !  direct effect
     &  ,l_use_seasalt_direct                                           &
     &  ,l_murk_rad                                                     &
     &  ,l_use_aod                                                      &
                                ! T if at least one of the aerosol
                                ! optical depth diags is requested
     &  ,l_use_biogenic                                                 &
     &  ,l_use_arclbiom, l_use_arclblck, l_use_arclsslt, l_use_arclsulp &
     &  ,l_use_arcldust, l_use_arclocff, l_use_arcldlta                 &
                                ! logicals for biogenic and other
                                ! aerosol climatologies
     &  ,l_use_dust                                                     &
     &  ,l_use_bmass_direct                                             &
     &  ,l_use_ocff_direct
!
!     &RUNDATA
!
      Character*8                                                       &
     &  exname_in                                                       &
                                ! Name of expt. to be read from
                                 !  previous run stored on tape up to
                                !  characters
     &  ,exname_out              ! Name of expt. to be written to tap
                                !  up to 6 characters
      Integer                                                           &
     &  ndayin                                                          &
                                ! No. of days in integration
     &  ,resdump_days                                                   &
                                ! frequency of dumps for restart
     &  ,runno_in                                                       &
                                ! Number of run to be read from
                                !  previous run stored on tape
     &  ,runno_out               ! Number of run to be written to tap
      Real                                                              &
     &  timestep                ! Model timestep for all physics
                                !  subroutines except radiation
!
!---------------------------------------------------------------------
!     Large scale observational forcing
!---------------------------------------------------------------------
!
!     Variables for diagnostic output for observational forcing
!
      Real                                                              &
     &  dap1(row_length, rows,36,model_levels)                          &
     &  ,dap2(row_length, rows,36,model_levels)                         &
     &  ,dap3(row_length, rows,36,nfor-1,model_levels)                  &
     &  ,dab1(row_length, rows,44)                                      &
     &  ,dab2(row_length, rows,44)                                      &
     &  ,dab3(row_length, rows,44,nfor-1)
!
!---------------------------------------------------------------------
!     Large scale statistical forcing
!---------------------------------------------------------------------
!
!     Random generator variables
!
      Integer                                                           &
     &  ntab                                                            &
                                  ! IN Dimension of array used in rand
                                !  generator.
     &  ,iv(ntab),iy,idum                                               &
                                ! On exit contains info on generator
     &  ,iseed                  ! Seed for random number generator

      Integer                                                           &
     &  dayno_wint                                                      &
                                ! Day number relative to winter
                                !  solstice
     &  ,ichgf                  ! No. of timesteps between change in
                                !  observational forcing
      Real                                                              &
     &  alfada(row_length, rows)                                        &
                                    ! Amplitude and mean of seasonal
     &  ,alfadb(row_length, rows)                                       &
                                    !  variation of tuning factor
     &  ,atime,btime                                                    &
                                ! Constants for calculating annual
                                !  cycle
     &  ,dbara(row_length, rows,nwet)                                   &
                                    ! Amplitude and mean of seasonal
     &  ,dbarb(row_length, rows,nwet)                                   &
                                    !  variation of mean dew pt.
                                !  depression (K)
     &  ,dgrada(row_length, rows,nwet)                                  &
                                    ! Amplitude and mean of seasonal
     &  ,dgradb(row_length, rows,nwet)                                  &
                                    !  variation of dew pt. depression
                                !  gradient (K km^-1)
     &  ,pa(row_length, rows, model_levels+1)                           &
                                    ! Amplitude and mean of seasonal
     &  ,pb(row_length, rows, model_levels+1)                           &
                                    !  variation of pressure
     &  ,tbara(row_length, rows,model_levels)                           &
                                    ! Amplitude and mean of seasonal
     &  ,tbarb(row_length, rows,model_levels)                           &
                                    !  variation of temp. (K)
     &  ,tgrada(row_length, rows,model_levels)                          &
                                    ! Amplitude and mean of seasonal
     &  ,tgradb(row_length, rows,model_levels)                          &
                                    !  variation of temp. gradient
                                !  (K Km^-1)
     &  ,tsda(row_length, rows,model_levels)                            &
                                    ! Amplitude and mean of seasonal
     &  ,tsdb(row_length, rows,model_levels)                            &
                                    !  variation of SD of temp. (K)
     &  ,vnbara(row_length, rows,model_levels)                          &
                                    ! Amplitude and mean of seasonal
     &  ,vnbarb(row_length, rows,model_levels)                          &
                                  !  variation of velocity VN (m s^-1)
     &  ,vnsda(row_length, rows,model_levels)                           &
                                    ! Amplitude and mean of seasonal
     &  ,vnsdb(row_length, rows,model_levels)                           &
                                    !  variation of SD of velocity VN
                                !  (m s^-1)
     &  ,vpbara(row_length, rows,model_levels)                          &
                                    ! Amplitude and mean of seasonal
     &  ,vpbarb(row_length, rows,model_levels)                          &
                                  !  variation of velocity VP (m s^-1)
     &  ,wbara(row_length, rows,ntrop)                                  &
                                    ! Amplitude and mean of seasonal
     &  ,wbarb(row_length, rows,ntrop)                                  &
                                    !  variation of SD of vert. vel.
                                !  (mb or HPa s^-1)
     &  ,wsda(row_length, rows,ntrop)                                   &
                                    ! Amplitude and mean of seasonal
     &  ,wsdb(row_length, rows,ntrop)
                                    !  variation of SD of vert. vel.
                                !  (mb s^-1)
                                !  roughness length (m)

!
      Real rhcrit(nwet)         ! Critical humidity for cloud
                                !  formation.

!---------------------------------------------------------------------
!     Local
!---------------------------------------------------------------------
!
      Integer                                                           &
     &  i, j, k, l, m, jlev                                             &
                                    ! Loop counters
     &  ,nresdump                                                       &
                                ! no. of restart dumps
     &  ,soil_index(row_length, rows)                                   &
                                    ! Index on land points
     &  ,temp_basis_days                                                &
                                ! for computation of andayy
     &  ,temp_basis_secs                                                &
     &  ,andayy1, andayy2                                               &
                                ! To calculate andday using time2sec
     &  ,dummy
      Character*80 cmessage     ! error message

      Integer :: ErrorStatus
      
      Real  ::                                                          &
        fsmc_limit                                                      &
                                ! limit for soil moisture init. FSMC
      , fsmc_min(nsoilp)                                                &
      , fsmc_max(nsoilp)                                                &
      , dzsoil_n(row_length, rows, nsoilm_levs)
                                ! Thickness of soil layers (m).

      Real, Intent(out) ::                                              &
        b_exp  (land_pts)                                               &
                                ! Land-point Clapp-Hornberger exponent
      , hcap   (land_pts)                                               &
                                ! Land-point soil heat capacity
      , hcon   (land_pts)                                               &
                                ! Land-point soil thermal conductivity
      , satcon (land_pts)                                               &
                                ! Land-point saturated hydrological
                                ! conductivity
      , sathh  (land_pts)                                               &
                                ! Land-point saturated soil water suction
      , v_sat  (land_pts)                                               &
                                ! Land-point volumetric soil moisture
                                ! content at saturation
      , v_wilt (land_pts)                                               &
                                ! Land-point volumetric soil moisture
                                ! content at wilting point
      , v_crit (land_pts)
                                ! Land-point volumetric soil moisture
                                ! content at the critical point


      INTEGER                                                           &
     & TILE_PTS(NTYPE)                                                  &
                                    ! LOCAL Number of land points which
!                                   !       include the nth surface type
     &,TILE_INDEX(LAND_PTS,NTYPE)! LOCAL Indices of land points which
!                                   !       include the nth surface type

!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!     Define site specific soil parameters and Initialise SWNOCZ
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!

        Do i=1, land_pts

          b_exp  (i) = b_exp_typ  (soil_type(i))
          hcap   (i) = hcap_typ   (soil_type(i))
          hcon   (i) = hcon_typ   (soil_type(i))
          satcon (i) = satcon_typ (soil_type(i))
          sathh  (i) = sathh_typ  (soil_type(i))

          ! Need these to calculate soil_layer_moisture (smcl)
          ! which comes from the dump at 5.2 for other model
          ! configurations.

          v_sat  (i) = v_sat_typ  (soil_type(i))
          v_wilt (i) = v_wilt_typ (soil_type(i))
          v_crit (i) = v_crit_typ (soil_type(i))

        End Do

        ! Need this to calculate soil_layer_moisture (smcl) which
        ! comes from the dump at 5.2 for other model configurations.
        Do k=1, nsoilt_levs
          Do j=1, rows
            Do i=1, row_length
              dzsoil_n(i,j,k) = dzsoil(k)
            End Do
          End Do
        End Do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!---------------------------------------------------------------------
!     Calculate number of days in year
!---------------------------------------------------------------------
!
! DEPENDS ON: time2sec
      Call time2sec(year_init,  1,1,0,0,0,0,0                           &
     &  ,andayy1, dummy, lcal360)
! DEPENDS ON: time2sec
      Call time2sec(year_init+1,1,1,0,0,0,0,0                           &
     &  ,andayy2, dummy, lcal360)
      andayy = andayy2 - andayy1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If (stats) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
!       Calculate day number relative to winter solstice (ie day 351)
!---------------------------------------------------------------------
!
        If (dayno_init  <   351) then
          dayno_wint = dayno_init + 9
        else
          dayno_wint = dayno_init - 351
        Endif
!
!---------------------------------------------------------------------
!       Derive initial data from climate datasets
!---------------------------------------------------------------------
       Do k = 1,model_levels+1
         Do j = 1, rows
           Do i = 1, row_length
             p(i,j,k) = p_in(i,j,k)
           Enddo
         Enddo
       Enddo
!
! DEPENDS ON: initstat
        Call INITSTAT(                                                  &
     &    row_length, rows, model_levels, nwet, ntrop,                  &
     &    andayy, dayno_wint, q, t, lat, long,                          &
     &    p_in, pa, pb, alfada, alfadb, tbara, tbarb,                   &
     &    tsda, tsdb, tgrada, tgradb, dbara, dbarb,                     &
     &    dgrada, dgradb, vnbara, vnbarb, vnsda, vnsdb,                 &
     &    vpbara, vpbarb, wbara, wbarb,                                 &
     &    wsda, wsdb, atime, btime, p_theta_levels)
!
       Do k = 1,model_levels
         Do j = 1, rows
           Do i = 1, row_length
             p_in(i,j,k) = p(i,j,k)
           Enddo
         Enddo
       Enddo
!---------------------------------------------------------------------
!       Initialise random generator
!---------------------------------------------------------------------
!
! DEPENDS ON: g05cbe
        Call G05CBE(iseed)
      Endif                     ! stats
!
!---------------------------------------------------------------------
!     Set initial data from &INPROF
!---------------------------------------------------------------------
!
! Routine is not be entered for geoforce !
      Do j = 1, rows
        Do i = 1, row_length
          If (stats .or. obs .or. noforce .or. geoforce) then
!  Do for stats as well for now as we have no start data otherwise.
            Do k = 1, model_levels
              u(i,j,k) = ui(i,j,k)
              v(i,j,k) = vi(i,j,k)
              w(i,j,k) = wi(i,j,k)
              w_adv(i,j,k) = w(i,j,k)
              t(i,j,k) = ti(i,j,k)
            Enddo
            w(i,j,0) = wi(i,j,0)
            w_adv(i,j,0) = w(i,j,0)
              Do k = 1, nwet
                q(i,j,k) = qi(i,j,k)
              Enddo
          Endif                   ! (obs .or. noforce)
        If (stats .and. altdat) then
            Do k = 1, model_levels
              t(i,j,k) = ti(i,j,k)
          Enddo
          Do k = 1, nwet
              q(i,j,k) = qi(i,j,k)
          Enddo
        Endif                   ! (stats .and. altdat)
        If (geoforce .and. geoinit) then
            Do k = 1, model_levels
              u(i,j,k) = ug(i,j)
              v(i,j,k) = vg(i,j)
          Enddo
        Endif                   ! geoforce and geoinit
        Enddo
      Enddo                     ! i
!
!---------------------------------------------------------------------
!     Calculate rates of change for large scale observational forcing
!---------------------------------------------------------------------
!
      If (obs) then
        Do j = 1, rows
          Do i = 1, row_length
            tstar(i,j) = tstar_forcing(i,j,1)
            Do l = 1, (nfor-1)
              tstar_forcing(i,j,l) = (tstar_forcing(i,j,l+1)            &
     &              -tstar_forcing(i,j,l))/(ichgf * timestep)
              ch_flux_h(i,j,l) = (flux_h(i,j,l+1) - flux_h(i,j,l))      &
     &          /              (ichgf * timestep)
              ch_flux_e(i,j,l) = (flux_e(i,j,l+1) - flux_e(i,j,l))      &
     &        /              (ichgf * timestep)
              Do  k = 1, model_levels
                ch_tls(i,j,l,k) = (tls(i,j,l+1,k) - tls(i,j,l,k))       &
     &        /              (ichgf * timestep)
                ch_uls(i,j,l,k) = (uls(i,j,l+1,k) - uls(i,j,l,k))       &
     &          /             (ichgf * timestep)
                ch_vls(i,j,l,k) = (vls(i,j,l+1,k) - vls(i,j,l,k))       &
     &          /             (ichgf * timestep)
                ch_wls(i,j,l,k) = (wls(i,j,l+1,k) - wls(i,j,l,k))       &
     &          /             (ichgf * timestep)
            Enddo               ! k
            Enddo                 ! l
            Do l = 1, (nfor-1)
            Do k = 1, nwet
                ch_qls(i,j,l,k) = (qls(i,j,l+1,k) - qls(i,j,l,k))       &
     &          /             (ichgf * timestep)
            Enddo               ! k
            Enddo                 ! l
            Do k = 1, model_levels
              tls(i,j,nfor,k)=tls(i,j,1,k)
              uls(i,j,nfor,k)=uls(i,j,1,k)
              vls(i,j,nfor,k)=vls(i,j,1,k)
              wls(i,j,nfor,k)=wls(i,j,1,k)
            end do
            Do k = 1, nwet
              qls(i,j,nfor,k)=qls(i,j,1,k)
            end do
            tstar_forcing(i,j,nfor) = tstar(i,j)
            flux_h(i,j,nfor)        = flux_h(i,j,1)
            flux_e(i,j,nfor)        = flux_e(i,j,1)

        Enddo                   ! i
        Enddo                    ! j
      Endif
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      If (tapein) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!---------------------------------------------------------------------
!       Read tape data if required to carry on from previous run
!---------------------------------------------------------------------
!
        Read (50) exname,runno,tapedump_no
        Write (6,205)exname,runno,tapedump_no
 205    Format                                                          &
     &    (' from tape header'/,' expt. name is                '        &
     &    ,a6,/,' run no. is                   ',i4,                    &
     &    /,' no. of dumps on tape are     ',i4)
!
!---------------------------------------------------------------------
!       Check for correct data set
!---------------------------------------------------------------------
!
        If (exname  ==  exname_in .and. runno  ==  runno_in) then
!
!---------------------------------------------------------------------
!         Look for correct day - tapeday_init input in namelist
!         INDATA.
!---------------------------------------------------------------------
!
          Do i = 1, tapedump_no
            If (stats) then
              Read (50) tapeday, resdump, iv, iy, idum
            elseif (obs) then
              Read (50) tapeday, resdump
            Endif
            Write (6,200) tapeday, tapeday_init
 200        Format (' tape year day= ',i4,'      start year day= ',i4)
            If (tapeday  ==  (tapeday_init-1) ) goto 999
!
!           If the end of the tape is reached and the specified day
!           tapeday_init not found o/p error message and stop run.
!
            If (i  ==  tapedump_no) then
              Write (6,201) exname, runno
 201          Format                                                    &
     &          (' initial day not found on data set m20.',a6,          &
     &          '.run',i3/)
              Stop
            Endif
          Enddo
 999      close (50)
!
!---------------------------------------------------------------------
!         Read initial data from tape in DUMP format.
!---------------------------------------------------------------------
!
! DEPENDS ON: dumpinit
          Call dumpinit(                                                      &
            row_length, rows, nprimvars, land_pts, model_levels, nwet,        &
            nbl_levs, nsoilt_levs, nsoilm_levs, ntrop, n_cca_levels,          &
            resdump, u, v, w, t, theta, q, qcl, qcf, layer_cloud,             &
            p, rho, t_deep_soil, smc, canopy_gb, snodep,                      &
            tstar, zh, z0msea,                                                &
            cca, rccb, rcct, smcl)
!
!         Sometimes (when initial wind in dump is arbitrary) we want
!         to reset the wind to geostrophic. To do this set geoinit to
!         true in logic namelist
!
          Do j = 1, rows
            Do i = 1, row_length
            If (geoinit .and. geoforce) then
                Do k = 1, model_levels
                  u(i,j,k) = ug(i,j)
                  v(i,j,k) = vg(i,j)
              Enddo
            Endif
              iccb(i,j) = int(rccb(i,j))
              icct(i,j) = int(rcct(i,j))
              tsi(i,j) = tstar(i,j)
            Enddo
          Enddo
!
!---------------------------------------------------------------------
!         G05CGE restores the state of the basic generator
!         routine G05DDE following the call to G05CFE at the end of
!         each day.
!---------------------------------------------------------------------
!
! DEPENDS ON: g05cge
          Call G05CGE(idum,iv,iy)
        else
          Write (6,203) exname_in, runno_in
 203      Format (' initial data set m20',a6,'.run',i3,' not found'/)
          Stop
        Endif                   ! (exname  ==  exname_in
                                ! .and. runno  ==  runno_in)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else                      ! not tapein
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Set initial values if no tape data to be used
!---------------------------------------------------------------------
!

        If (land_pts > 0) Then 
!-----------------------------------------------------------------------
! Call TILEPTS to initialise TILE_PTS and TILE_INDEX
!-----------------------------------------------------------------------
! DEPENDS ON: tilepts
          CALL TILEPTS(LAND_PTS,frac_typ,TILE_PTS,TILE_INDEX)
!-----------------------------------------------------------------------
! Initialise tiled and gridbox mean vegetation parameters
!-----------------------------------------------------------------------
          WRITE(6,*) 'INITVEG: CALLING SPARM'
! DEPENDS ON: sparm
          CALL SPARM (LAND_PTS,NTYPE,CAN_MODEL,TILE_PTS,TILE_INDEX,     &
                      frac_typ, canht, lai, satcon,                     &
                      catch_snow, catch, infil_tile, z0_tile )

        End If

        Do k=1, model_levels+1
          Do j=1, rows
            Do i=1, row_length
              p(i,j,k) = p_in(i,j,k)
            End Do
          End Do
        End Do

        Do k=1, n_cca_levels
          Do j=1, rows
            Do i=1, row_length
              cca(i,j,k) = ccai(i,j)
            End Do
          End Do
        End Do

        Do j=1, rows
          Do i=1, row_length
            tstar  (i,j) = tstari  (i,j)
            tsi    (i,j) = tstar   (i,j)
            z0msea (i,j) = z0mseai (i,j)
            iccb   (i,j) = iccbi   (i,j)
            icct   (i,j) = iccti   (i,j)
            snodep (i,j) = snodepi (i,j) ! Includes snow depth on land/sea
          End Do
        End Do


        If (land_pts > 0) Then

          !--------------------------------------------------------------
          ! Initialise the canopy water
          !--------------------------------------------------------------
          Do i=1, land_pts
            If (canopy_gbi(i) == RMDI) Then
              write (*,'(5(TR1,A52/))')                                       &
                '===================================================='        &
              , '| Initial Gridbox Mean Canopy Water (canopy_gbi)   |'        &
              , '| must be specified for each land point in         |'        &
              , '| namelist.                                        |'        &
              , '===================================================='
              Stop                    
            Else
              canopy_gb(i) = canopy_gbi(i)
            End If
          End Do    
       
          !--------------------------------------------------------------
          ! Initialise the soil moisture content (SMC)
          !--------------------------------------------------------------
          Do i=1, land_pts
            If (smci(i) == RMDI) Then
              write (*,'(4(TR1,A52/))')                                       &
                '===================================================='        &
              , '| Initial Soil Moisture Content (smci) must be     |'        &
              , '| specified for each land point in namelist.       |'        &
              , '===================================================='
              Stop                    
            Else
              smc(i) = smci(i)
            End If
          End Do

          !--------------------------------------------------------------
          ! Initialise the deep soil temperature
          !--------------------------------------------------------------
          Do k=1, nsoilm_levs
            Do i=1, land_pts
              If (t_deep_soili(i,k) == RMDI) Then
                write (*,'(2(TR1,A52/),(TR1,A15,I2,A18,T53,A1/),'//           &
                         ' 2(TR1,A52/))')                                     &
                  '===================================================='      &
                , '| Initial Deep Soil Temperature (t_deep_soili)     |'      &
                , '| must specify ',nsoilt_levs,' soil temperatures','|'      &
                , '| for each land point in namelist.                 |'      &
                , '===================================================='
                Stop                    
              Else
                t_deep_soil(i,k) = t_deep_soili(i,k)
              End If
            End Do    
          End Do               ! nsoilm_levs

          !--------------------------------------------------------------
          ! Initialise the soil moisture in root zone layers (SMCL)
          !--------------------------------------------------------------
 
          Select Case (smi_opt)
            Case(smc_spec)
              ! Soil moisture (smcli) to be specified by user namelist. 

              ! Check for missing data
              Do k=1, nsoilm_levs
                Do i=1, land_pts
                  If (smcli(i,k) == RMDI) Then
                    write (*,'(2(TR1,A52/),'//                                &
                             '  (TR1,A25,I2,A21,T53,A1/),'//                  &
                             ' 2(TR1,A52/))')                                 &
                      '===================================================='  &
                    , '| Initial Soil Moisture Content in layers (smcli)  |'  &
                    , '| must specify values on ',nsoilm_levs                 &
                    ,                            ' soil moisture levels','|'  &
                    , '| for each land point in namelist.                 |'  &
                    , '===================================================='
                    Stop                    
                  Else
                    smcl(i,k) = smcli(i,k)
                  End If
                End Do
              End Do

            Case(smc_fsmc)
              ! Calculate max & min limits of FSMC
              Do i=1, nsoilp
                fsmc_min(i) = -  v_wilt_typ(i)                                &
                              / (v_crit_typ(i) - v_wilt_typ(i))
                fsmc_max(i) =   (v_sat_typ(i)  - v_wilt_typ(i))               &
                              / (v_crit_typ(i) - v_wilt_typ(i))
              End Do

              ! Soil moisture to be initialised by soil stress factor
              ! (FSMC) from namelist.


              Do i=1, land_pts

                ! Check for missing data
                If (fsmc(i) == RMDI) Then
                  write (*,'(4(TR1,A52/))')                                   &
                    '===================================================='    &
                  , '| Soil moisture stress factor (fsmc) must be       |'    &
                  , '| specified for each land point in namelist.       |'    &
                  , '===================================================='
                  Stop
                End If

                ! Check fsmc is within range for given soil type
                If ((fsmc(i) < fsmc_min(soil_type(i))) .or.                   &
                    (fsmc(i) > fsmc_max(soil_type(i)))) Then

                  write (*,'(5(TR1,A52/),'//                                  &
                           ' 4(TR1,A19,F6.3,A10,F6.3,T53,A1/),'//             &
                           ' 2(TR1,A52/))')                                   &
                    '===================================================='    &
                  , '| Soil moisture stress factor (fsmc) is outside    |'    &
                  , '| acceptable range. For the specified soil types,  |'    &
                  , '| fsmc should be in ranges:                        |'    &
                  , '|                                                  |'    &
                  , '|  1) Ice        : ',fsmc_min(1),' < fsmc < '            &
                                         ,fsmc_max(1),                 '|'    &
                  , '|  2) Clay       : ',fsmc_min(2),' < fsmc < '            &
                                         ,fsmc_max(2),                 '|'    & 
                  , '|  3) Loam       : ',fsmc_min(3),' < fsmc < '            &
                                         ,fsmc_max(3),                 '|'    &
                  , '|  4) Loamy Sand : ',fsmc_min(4),' < fsmc < '            &
                                         ,fsmc_max(4),                 '|'    &
                  , '|                                                  |'    &
                  , '===================================================='

                  Stop
                End If
              End Do

              Do k=1, nsoilm_levs
                Do i=1, land_pts
                  smcl(i,k) = rho_water*dzsoil(k)                             &
                            * ( fsmc(i)*v_crit(i) + (1-fsmc(i))*v_wilt(i) )
                End Do
              End Do

            Case(smc_sth)
              ! Soil moisture is to be initialised by total soil
              ! moisture as a fraction of saturation (STH).

              Do k=1, nsoilm_levs
                Do i=1, land_pts

                  ! Check for missing data
                  If (sth(i,k) == RMDI) Then
                    write (*,'(2(TR1,A52/),'//                                &
                             '  (TR1,A37,I2,A5,T53,A1/),'//                   &
                             ' 2(TR1,A52/))')                                 &
                      '===================================================='  &
                    , '| Total Soil Moisture (sth) as a fraction of       |'  &
                    , '| saturation, must specify values on ',nsoilm_levs     &
                    ,                                            ' soil','|'  &
                    , '| moisture levels for each land point in namelist. |'  &
                    , '===================================================='
                    Stop
                  End If

                  ! Check sth range
                  If ((sth(i,k) < 0.0) .or.                                   &
                      (sth(i,k) > 1.0)) Then
                    write (*,'(5(TR1,A52/))')                                 &
                      '===================================================='  &
                    , '| Specified values for Total Soil Moisture (sth)   |'  &
                    , '| as a fraction of saturation, must be within the  |'  &
                    , '| range: 0.0 ~ 1.0                                 |'  &
                    , '===================================================='
                    Stop
                  End If

                  smcl(i,k) = sth(i,k) * rho_water * dzsoil(k) * v_sat(i)

                End Do
              End Do

            Case(IMDI)

              write (*,'(5(TR1,A52/))')                                       &
                '===================================================='        &
              , '| An initialization option for soil moisture       |'        &
              , '| content in layers (smi_opt) must be specified    |'        &
              , '| when running over land points                    |'        &
              , '===================================================='
              Stop

          End Select
        End If                    ! land_pts

      Endif                     ! tapein

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if defined(A08_7A)
!---------------------------------------------------------------------
!     Initialise the stomatol conductance for MOSES.
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!     Initialise the frozen and unfrozen soil moisture for MOSES
!---------------------------------------------------------------------
! DEPENDS ON: freeze_soil
      Call FREEZE_SOIL (land_pts, nsoilm_levs, b_exp, dzsoil, sathh,    &
                        smcl, t_deep_soil, v_sat, sthu, sthf)
#endif

!
!---------------------------------------------------------------------
!     Calculate pressure, exner_theta_levels and pstar
!---------------------------------------------------------------------
!
! DEPENDS ON: calc_press
      Call Calc_press                                                   &
! Input data
     &    (model_levels, nwet, rows, row_length, p, theta,              &
     &      q, r_theta_levels, r_rho_levels, .true., .true.,            &
! In/Out
     &     rho,                                                         &
! Output data
     &     exner_theta_levels, exner_rho_levels, p_theta_levels,        &
     &     rp, rp_theta, pstar)

!
!---------------------------------------------------------------------
!     Calculate DELTAP
!---------------------------------------------------------------------
!
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            deltap(i,j,k) = p(i,j,k+1)-p(i,j,k)
          Enddo
        End do
      Enddo
!---------------------------------------------------------------------
!     Initialise cloud water (QCL,QCF)
!---------------------------------------------------------------------
!
      If (.not. tapein) then
        Do jlev = 1, nwet
! DEPENDS ON: initqlcf
          Call INITQLCF                                                 &
     &      (p_theta_levels(1,1,jlev),rhcrit,q(1,1,jlev),               &
     &      t(1,1,jlev),model_levels,                                   &
     &      row_length, rows,layer_cloud(1,1,jlev),                     &
     &      qcf(1,1,jlev),qcl(1,1,jlev),nbl_levs,jlev)
        Enddo
      End If

!
!---------------------------------------------------------------------
!       Convert temperature to potential temperature and tls to
!       theta_star
!---------------------------------------------------------------------
!
! Set T to Theta and TLS increment to theta_star !!!!
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta(i,j,k) = t(i,j,k)/exner_theta_levels(i,j,k)
            theta_star(i,j,k) = tls(i,j,nfor,k) /                       &
     &                            exner_theta_levels(i,j,k)
          Enddo
        Enddo
      Enddo
!
!---------------------------------------------------------------------
!     Zero diagnostics for observational forcing
!---------------------------------------------------------------------
!
      If (obs .and. prindump_obs) then
        Do i = 1, row_length
          Do m = 1, rows
            Do k = 1, model_levels
            Do j = 1, 36
                dap1(i,m,j,k) = 0.0
                dap2(i,m,j,k) = 0.0
              Do  l = 1, nfor-1
                  dap3(i,m,j,k,l) = 0.0
              Enddo
            Enddo
          Enddo
          Do j = 1, 44
              dab1(i,m,j) = 0.0
              dab2(i,m,j) = 0.0
            Do k = 1, nfor-1
                dab3(i,m,j,k) = 0.0
          Enddo
          Enddo
        Enddo                   ! i
        Enddo                    ! m
      Endif                     ! obs

!
!---------------------------------------------------------------------
!     Radiation - LWLKIN in UM code deck LWTRAN, SWLKIN in UM code
!     SWTRAN.
!     Initialise the LW and SW spectral band tables - standard
!     Radiation code.
!     For Edwards-Slingo radiation code 3A use R2_LW_SPECIN and
!     R2_SW_SPECIN .
!     (in deck SPIN3A) to pick up the spectral files from units 57
!     (SW) and 80 (LW).
!---------------------------------------------------------------------
!     Longwave
      ErrorStatus = 0
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      Write(6,*) "WARNING:  Names of the spectral files are now"
      Write(6,*) "read in through namelists R2SWCLNL and R2LWCLNL"

         Do j=1,n_lwcall
! DEPENDS ON: r2_lw_specin
           Call r2_lw_specin(icode, cmessage                            &
     &        , lw_control(j)%spectral_file                             &
     &        , lw_control(j)%l_ch4, lw_control(j)%l_n2o                &
     &        , lw_control(j)%l_cfc11, lw_control(j)%l_cfc12            &
     &        , lw_control(j)%l_cfc113                                  &
     &        , lw_control(j)%l_hcfc22, lw_control(j)%l_hfc125          &
     &        , lw_control(j)%l_hfc134A                                 &
     &        , l_climat_aerosol, l_use_dust, l_use_arcldust            &
     &        , l_use_sulpc_direct, l_use_arclsulp                      &
     &        , l_use_soot_direct, l_use_arclblck                       &
     &        , l_use_bmass_direct, l_use_arclbiom                      &
     &        , l_use_seasalt_direct, l_use_arclsslt                    &
     &        , l_use_ocff_direct, l_use_arclocff                       &
     &        , l_use_biogenic, l_use_arcldlta                          &
     &        , l_murk_rad                                              &
     &        , l_use_aod, lw_control(j)%l_gas                          &
     &        , lw_control(j)%l_continuum, lw_control(j)%l_drop         &
     &        , lw_control(j)%l_aerosol, lw_control(j)%l_ice            &
     &        , lw_spectrum(j)                                          &
     &        )
          Enddo
#else
! DEPENDS ON: r2_lw_specin
      Call r2_lw_specin(ErrorStatus, cmessage                           &
     &  , l_ch4, l_n2o, l_cfc11, l_cfc12                                &
     &  , l_cfc113, l_hcfc22, l_hfc125, l_hfc134a                       &
     &  , l_climat_aerosol, l_use_dust, l_use_arcldust                  &
     &  , l_use_sulpc_direct, l_use_arclsulp                            &
     &  , l_use_soot_direct, l_use_arclblck                             &
     &  , l_use_bmass_direct, l_use_arclbiom                            &
     &  , l_use_seasalt_direct, l_use_arclsslt                          &
     &  , l_use_ocff_direct, l_use_arclocff                             &
     &  , l_use_biogenic, l_use_arcldlta                                &
     &  , l_murk_rad, l_use_aod)
#endif
      If (ErrorStatus  /=  0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,ErrorStatus,Cmessage)
      Endif
!     Shortwave
      ErrorStatus = 0
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      Write(6,*) "WARNING:  Names of the spectral files are now  "
      Write(6,*) "read in through namelists R2SWCLNL and R2LWCLNL"

         Do j=1, n_swcall
! DEPENDS ON: r2_sw_specin
           Call r2_sw_specin(icode, cmessage                            &
     &        , sw_control(j)%spectral_file                             &
     &        , sw_control(j)%l_o2                                      &
     &        , l_climat_aerosol, l_use_dust, l_use_arcldust            &
     &        , l_use_sulpc_direct, l_use_arclsulp                      &
     &        , l_use_soot_direct, l_use_arclblck                       &
     &        , l_use_bmass_direct, l_use_arclbiom                      &
     &        , l_use_seasalt_direct, l_use_arclsslt                    &
     &        , l_use_ocff_direct, l_use_arclocff                       &
     &        , l_use_biogenic, l_use_arcldlta                          &
     &        , l_murk_rad                                              &
     &        , sw_control(j)%l_rayleigh, sw_control(j)%l_gas           &
     &        , sw_control(j)%l_continuum, sw_control(j)%l_drop         &
     &        , sw_control(j)%l_aerosol, sw_control(j)%l_ice            &
     &        , sw_spectrum(j)                                          &
     &        )
         Enddo
#else
! DEPENDS ON: r2_sw_specin
      Call r2_sw_specin(ErrorStatus, cmessage                           &
     &  , l_o2                                                          &
     &  , l_climat_aerosol, l_use_dust, l_use_arcldust                  &
     &  , l_use_sulpc_direct, l_use_arclsulp                            &
     &  , l_use_soot_direct, l_use_arclblck                             &
     &  , l_use_bmass_direct, l_use_arclbiom                            &
     &  , l_use_seasalt_direct, l_use_arclsslt                          &
     &  , l_use_ocff_direct, l_use_arclocff                             &
     &  , l_use_biogenic, l_use_arcldlta                                &
     &  , l_murk_rad)
#endif



      If (ErrorStatus  /=  0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,ErrorStatus,Cmessage)
      Endif
!
!---------------------------------------------------------------------
!     Write restart dump information to tape
!---------------------------------------------------------------------
!
      If (tapeout) then
        nresdump = int( ndayin / resdump_days)
        If (mod(ndayin, resdump_days)  /=  0) nresdump = nresdump + 1
        Write (55) exname_out,runno_out,nresdump
      Endif

      Return
      END SUBROUTINE RUN_INIT
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#endif
