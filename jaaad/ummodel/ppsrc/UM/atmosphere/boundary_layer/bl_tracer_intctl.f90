
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      Subroutine BL_TRACER_INTCTL(                                      &
     &        row_length, rows, bl_levels, off_x, off_y, at_extremity   &
     &,       dtrdz_charney_grid,r_rho_levels,r_theta_levels            &
     &,       halo_i,halo_j                                             &
     &,       tr_levels, tr_vars, model_levels                          &
     &,       DIM_CS2                                                   &
     &,       alpha_cd, rhokh_mix                                       &
     &,       p_star, p, timestep                                       &
! Control logicals
     &,       L_MURK,L_MURK_ADVECT,L_BL_TRACER_MIX,L_DUST,L_CAM_DUST    &
     &,       L_SULPC_SO2, l_sulpc_nh3, l_sulpc_dms, l_soot, l_biomass  &
     &,       l_ocff, l_co2_interactive                                 &
     &,       l_co2_emits                                               &
! rml 2/7/13
     &,       L_CO2_TRACER                                              &
     &,       L_ukca, L_USE_CARIOLLE                                    &
! Fields to mix
     &,       murk, free_tracers, OZONE_TRACER                          &
! Mineral Dust
     &,        DUST_DIV1,DUST_DIV2,DUST_DIV3                            &
     &,        DUST_DIV4,DUST_DIV5,DUST_DIV6                            &
! Sulphur cycle
     &,       so2, dms, so4_aitken, so4_accu, so4_diss, nh3             &
! Soot cycle
     &,       soot_new, soot_aged, soot_cld                             &
! Biomass aerosol
     &,       bmass_new, bmass_agd, bmass_cld                           &
! Fossil-fuel organic carbon aerosol
     &,       ocff_new, ocff_agd, ocff_cld                              &
! Carbon cycle
     &,       co2, co2_emits, co2flux, npp, resp_s                      &
     &,       co2_flux_tot, land_co2                                    &
! Tracer fluxes - kdcorbin, 05/10
     &,       tracer_flux1, tracer_flux2, tracer_flux3, tracer_flux4    &
     &,       tracer_flux5, tracer_flux6, tracer_flux7, tracer_flux8    &
     &,       tracer_flux9, tracer_flux10,tracer_flux11,tracer_flux12   &
     &,       tracer_flux13,tracer_flux14,tracer_flux15,tracer_flux16   &
     &,       tracer_flux17,tracer_flux18,tracer_flux19,tracer_flux20   &     
! Emissions fields
     &, DUST_FLUX                                                       &
     &,       so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em    &
     &,       ocff_hilem, ocff_em, drydep_str                           &
     &,       kent, we_lim, t_frac, zrzi                                &
     &,       kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                &
     &,       zh, zhsc, z_half                                          &
! IN for dry deposition of tracers
     &,       rho_aresist, aresist, resist_b                            &
     &,       R_B_DUST,TSTAR                                            &
     &,       land_points, land_index, ice_frac                         &
! IN variables for MOSES II
     &,       ntype, ntiles, tile_pts, tile_index, tile_frac            &
     &,       canopy, catch, snow_tile, gc                              &
     &,       aresist_tile, resist_b_tile, FLANDG                       &
     &,                                                                 &
! STASH related variables
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &        stashwork                                                 &
     &        )
!
! ---------------------------------------------------------------------
! Purpose: Intermediate control routine to call BL_TRMIX_DD for tracer
!          variables that require boundary layer mixing and/or dry
!          deposition with MOSES II surface scheme.
!
!          Called by IMP_CTL2.
!
! Current code owner: M Woodage
!
! History:
! Version   Date   Comment
! -------   ----   -------
!  5.3   20/10/01  New deck.                               M Woodage
!  5.4   29/08/02  Extra variables passed through for changes
!                  to tr_mix.                   P. Davison.
!    5.5  26/02/03  Pass 3D CO2 to tracer mixing routine. C Jones
!  5.5   05/02/03  Pass through extra arguments for biomass aerosol
!                  scheme.                           P Davison
!  5.5   19/02/03  Remove redundent RML code. Adrian Lock
!    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
!   6.2   08/07/05  Added logical and call to copy BL variables to
!                   ukca module.          F. O'Connor
!  6.2  01/03/06  pass CO2 fluxes back up for use in diagnostics
!                 routines.                             C.D. Jones
!  6.1   11/01/06  Remove MOSES I code                A.P.Lock
!
!
!
! Code description:
!   Language: FORTRAN 77  + common extensions
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      Implicit None

!C_DUST_NDIV.............................................................
! Description: Contains parameters for mineral dust code
! Current Code Owner: Stephanie Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  5.5      12/02/03  Original Code.   Stephanie Woodward
!
! Declarations:
!
      INTEGER NDIV        ! number of particle size divisions
      PARAMETER (NDIV = 6)
!.....................................................................
!
! Parallel setup variables
      Integer, Intent(In) ::                                            &
     &  off_x                                                           &
                              ! Size of small halo in i
     &, off_y                                                           &
                              ! Size of small halo in j.
     &, halo_i                                                          &
                              ! Size of halo in i direction
     &, halo_j                ! Size of halo in j direction

       Logical, Intent(In) ::                                           &
     &  at_extremity(4)   ! Indicates if this processor is at north,
                          ! south, east or west of the processor grid

! Model dimensions
      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, bl_levels                                                       &
     &, tr_levels                                                       &
                          ! number of free tracer levels
     &, tr_vars                                                         &
                          ! number of free tracer variables
     &, DIM_CS2         ! soil carbon dimension


! Model switches
      Logical, Intent(In) ::                                            &
     &  L_Murk                                                          &
                          ! Switch for (visibility) aerosol
     &, L_murk_advect                                                   &
                          ! Switch for advecting aerosol
     &, L_bl_tracer_mix                                                 &
                          ! Switch for BL mixing of free tracers
     &, L_DUST                                                          &
                          ! Switch for Mineral Dust
     &, L_CAM_DUST                                                      &
                          !Old version of dust_uplift scheme used in CAM NWP models
     &, L_sulpc_so2                                                     &
                          ! Switch for Sulphur Cycle
     &, L_sulpc_nh3                                                     &
                          ! NH3 included in Sulphur Cycle
     &, L_sulpc_dms                                                     &
                          ! DMS included in Sulphur Cycle
     &, L_soot                                                          &
                          ! Switch for Soot Cycle
     &, L_biomass                                                       &
                          ! Switch for Biomass aerosol
     &, L_ocff                                                          &
                          ! Switch for Fossil-fuel OC aerosol
     &, L_co2_interactive                                               &
                          ! Switch for interactive CO2
     &, L_co2_emits                                                     &
                          ! Switch for anthro CO2 emissions
! rml 2/7/13
     &, L_CO2_tracer                                                    &
                          ! Switch to put co2 flux into passive tracer
     &, L_ukca                                                          &
                          ! switch for UKCA
     &, L_USE_CARIOLLE    
                          ! Switch for the cariolle ozone tracer scheme
!
! Model parameters
      Real, Intent(In) ::                                               &
     &  timestep                                                        &
     &, alpha_cd(bl_levels)                                             &
     &, rhokh_mix (row_length, rows, bl_levels)                         &
     &, drydep_str(row_length, rows)                                    &
     &, we_lim(row_length,rows,3)                                       &
                                      !rho*entrainment rate implied by
                                      !placing of subsidence
     &, zrzi(row_length,rows,3)                                         &
                                      !(z-z_base)/(z_i-z_base)
     &, t_frac(row_length,rows,3)                                       &
                                      !a fraction of the timestep
     &, we_lim_dsc(row_length,rows,3)                                   &
                                      !rho*entrainment rate implied by
                                      !  placing of subsidence
     &, zrzi_dsc(row_length,rows,3)                                     &
                                      !(z-z_base)/(z_i-z_base)
     &, t_frac_dsc(row_length,rows,3)                                   &
                                      !a fraction of the timestep
     &, z_half(row_length,rows,BL_LEVELS)                               &
                                      !Z_HALF(*,K) is height of half
                                      !    level k-1/2.
     &, zhsc(row_length,rows)                                           &
                                      !Top of decoupled layer
     &, zh  (row_length,rows)                                           &
                                      !Top of surface mixed layer
     &, dtrdz_charney_grid(row_length,rows,BL_LEVELS)
                                      ! For tracer mixing

      Integer, Intent(In) ::                                            &
     &  kent(row_length,rows)                                           &
                                      !grid-level of SML inversion
     &, kent_dsc(row_length,rows)                                       &
                                      !grid-level of DSC inversion
     &, land_points                                                     &
                                      !No.of land points being processed
     &, land_index(land_points)                                         &
                                      !set from land_sea_mask
! For MOSES II
     &, ntype                                                           &
                                      !No. of tile types
     &, ntiles                                                          &
                                      !No.of land-surface tiles
     &, tile_pts(ntype)                                                 &
                                      !No.of tile points
     &, tile_index(land_points,ntype) !Index of tile points
!
! Required fields (input)
      Real, Intent(In) ::                                               &
     &  p_star(row_length, rows)                                        &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:bl_levels)              &
                                      ! Heights of theta levels
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, bl_levels)
                                      ! Heights of rho levels


! Mineral Dust Emissions fields (input)
      Real, Intent(In) ::                                               &
     & DUST_FLUX(ROW_LENGTH, ROWS, NDIV)! dust emission flux

! Emissions fields (input)
      Real, Intent(In) ::                                               &
     &  so2_hilem ( row_length, rows )                                  &
     &, so2_em    ( row_length, rows )                                  &
     &, nh3_em    ( row_length, rows )                                  &
     &, dms_em    ( row_length, rows )                                  &
     &, soot_hilem( row_length, rows )                                  &
     &, soot_em   ( row_length, rows )                                  &
     &, ocff_hilem( row_length, rows )                                  &
     &, ocff_em   ( row_length, rows )                                  &
     &, co2_emits ( row_length, rows )                                  &
                                        ! anthro CO2 emissions
     &, co2flux   ( row_length, rows )  ! ocean  CO2 emissions

! Tracer fluxes for atmospheric tracers - kdcorbin, 05/10
      Real, Intent(In) ::                                               &
     &   tracer_flux1(row_length,rows),  tracer_flux2(row_length,rows)   &
     &,  tracer_flux3(row_length,rows),  tracer_flux4(row_length,rows)   &
     &,  tracer_flux5(row_length,rows),  tracer_flux6(row_length,rows)   &
     &,  tracer_flux7(row_length,rows),  tracer_flux8(row_length,rows)   &
     &,  tracer_flux9(row_length,rows),  tracer_flux10(row_length,rows)  &
     &,  tracer_flux11(row_length,rows), tracer_flux12(row_length,rows)  &
     &,  tracer_flux13(row_length,rows), tracer_flux14(row_length,rows)  &
     &,  tracer_flux15(row_length,rows), tracer_flux16(row_length,rows)  &
     &,  tracer_flux17(row_length,rows), tracer_flux18(row_length,rows)  &
     &,  tracer_flux19(row_length,rows), tracer_flux20(row_length,rows)

!
! For dry deposition of tracers (input)
      Real, Intent(In) ::                                               &
     &  rho_aresist(row_length,rows)                                    &
                                       !RHOSTAR*CD_STD*VSHR
     &, aresist(row_length,rows)                                        &
                                       !1/(CD_STD*VSHR)
     &, resist_b(row_length,rows)                                       &
                                       !(1/CH-1/CD_STD)/VSHR
     &, R_B_DUST(ROW_LENGTH,ROWS,NDIV)                                  &
                                       ! surf layer res for dust
     &, TSTAR(ROW_LENGTH,ROWS)                                          &
                               ! temperature
     &, ice_frac(row_length,rows)                                       &
                                       !sea_ice fractn (ice/qrclim.ice)
!
! For MOSES II
     &, tile_frac(land_points,ntiles)                                   &
                                       !fractional coverage for each
!                                       surface tile
     &, canopy(land_points,ntiles)                                      &
                                       !Surface/canopy water (kg/m2)
     &, catch(land_points,ntiles)                                       &
                                       !Surface/canopy water capacity
!                                       of snow-free land tiles (kg/m2)
     &, snow_tile(land_points,ntiles)                                   &
                                       !Snow on tiles (kg/m2).
     &, gc(land_points,ntiles)                                          &
                                       !Stomatal conductance to evapn
!                                       for land tiles (m/s)
     &, aresist_tile(land_points,ntiles)                                &
!                                      ! 1/(CD_STD*VSHR) on land tiles
     &, resist_b_tile(land_points,ntiles)                               &
!                                    !(1/CH-1/CD_STD)/VSHR on land tiles
     &, FLANDG(row_length,rows)                                         &
                                     !land fraction of grid box
!                                       (0 if all sea, 1 if all land)
     &, npp(land_points)                                                &
                                     ! net primary productivity
     &, resp_s(DIM_CS2)                                                 &
                                       ! soil respiration
     &, co2_flux_tot(row_length,rows)                                   &
                                       ! OUT total CO2 flux
     &, land_co2(land_points)          ! OUT terrestrial CO2 flux
! Tracer fields for mixing
      Real, Intent(InOut) ::                                            &
     &  murk        ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,free_tracers( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y,                         &
     &                tr_levels, tr_vars)

      Real, Intent(InOut) ::                                            &
     &  DUST_DIV1   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV2   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV3   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV4   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV5   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV6   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )

      Real, Intent(InOut) ::                                            &
     &  so2         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,dms         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_aitken  ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_accu    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_diss    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,nh3         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )
!
      Real, Intent(InOut) ::                                            &
     &  soot_new    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,soot_aged   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,soot_cld    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_new   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_agd   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_cld   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_new    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_agd    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_cld    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,co2         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,OZONE_TRACER( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          

      Real, Intent(InOut) ::                                            &
     &  stashwork(*)

! Argument comdeck declaration
! CSUBMODL start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.5    07/04/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!  1. Internal model and submodel dump partition identifiers - fixed
!     for all experiments.
! CSMID start
!
! Description:
!    Hold parameters defining internal model identifiers and submodel
!    data partition (ie main D1 data array and consequent dump), both
!    short and long form.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.3    26/10/93   M. Carter. Part of an extensive mod that:
!                    1.Removes the limit on primary STASH item numbers.
!                    2.Removes the assumption that (section,item)
!                      defines the sub-model.
!                    3.Thus allows for user-prognostics.
!                    Add index to submodel home dump.
! 3.5    13/03/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
! 6.0    02/07/03   Add X_IM and X_SM for small exec.      E.Leung
!
! Declarations:
!
!   Hold parameters defining internal model identifiers and submodel
!   data partition (ie main D1 data array and consequent dump), both
!   short and long form
      ! Internal models
      INTEGER,PARAMETER:: A_IM      = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: ATMOS_IM  = 1 ! Atmosphere internal model
      INTEGER,PARAMETER:: O_IM      = 2 ! Ocean internal model
      INTEGER,PARAMETER:: OCEAN_IM  = 2 ! Ocean internalmodel
      INTEGER,PARAMETER:: S_IM      = 3 ! Slab internal model
      INTEGER,PARAMETER:: SLAB_IM   = 3 ! Slab internal model
      INTEGER,PARAMETER:: W_IM      = 4 ! Wave internal model
      INTEGER,PARAMETER:: WAVE_IM   = 4 ! Wave internal model
      INTEGER,PARAMETER:: I_IM      = 5 ! Sea=ice internal model
      INTEGER,PARAMETER:: SEAICE_IM = 5 ! Sea=ice internal model
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_IM      = 6 ! ND internal model
      INTEGER,PARAMETER:: NATMOS_IM = 6 ! ND internal model
      ! Small Executables
      INTEGER,PARAMETER:: X_IM      = 7 ! SX indicator

      ! Submodels
      INTEGER,PARAMETER:: A_SM      = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: ATMOS_SM  = 1 ! Atmosphere submodel
      INTEGER,PARAMETER:: O_SM      = 2 ! Ocean submodel
      INTEGER,PARAMETER:: OCEAN_SM  = 2 ! Ocean submodel
      INTEGER,PARAMETER:: W_SM      = 4 ! Wave submodel
      INTEGER,PARAMETER:: WAVE_SM   = 4 ! Wave submodel
      ! New dynamics (Charney-Phillips grid)
      INTEGER,PARAMETER:: N_SM      = 6 ! ND submodel
      INTEGER,PARAMETER:: NATMOS_SM = 6 ! ND submodel
      ! Small Executables
      INTEGER,PARAMETER:: X_SM      = 7 ! SX indicator

! CSMID end

!
!  2. Maximum internal model/submodel array sizes for this version.
!
! CSUBMAX start
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.5    13/07/95   Original code. D.M. Goddard
! 4.0     3/11/95   Reduce max internal model, submodel from 10 to 4
!                   to save space in model. At 4.0 the max no of
!                   supported models is 3, 1 slot is reserved for
!                   expansion. Rick Rawlins.
!  4.1  21/02/96  Wave model introduced as 4th sub-model.  RTHBarnes
!
! Declarations:
!
!
!  1. Maximum internal model/submodel array sizes for this version.
!
      ! Max no. of internal models
      INTEGER,PARAMETER:: N_INTERNAL_MODEL_MAX=4

      ! Max no. of submodel dump partitions
      INTEGER,PARAMETER:: N_SUBMODEL_PARTITION_MAX=4

      ! Max value of internal model id
      INTEGER,PARAMETER:: INTERNAL_ID_MAX=N_INTERNAL_MODEL_MAX

      ! Max value of submodel dump id
      INTEGER,PARAMETER:: SUBMODEL_ID_MAX=N_SUBMODEL_PARTITION_MAX

! CSUBMAX end
!
!  3. Lists of internal models and their submodel dump partitions -
!     initialised by the user interface - experiment specific.
      INTEGER :: N_INTERNAL_MODEL          ! No. of internal models
      INTEGER :: N_SUBMODEL_PARTITION      ! No. of submodel partitions

      ! Internal models
      INTEGER :: INTERNAL_MODEL_LIST(N_INTERNAL_MODEL_MAX)

      ! Submodel identifier for each internal model in list
      INTEGER :: SUBMODEL_FOR_IM    (N_INTERNAL_MODEL_MAX)

      ! Submodel number for each submodel id
      INTEGER :: SUBMODEL_FOR_SM(N_INTERNAL_MODEL_MAX)

      ! Namelist for information in 3.
      NAMELIST/NSUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,          &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM

      ! 4. Lists calculated in model from user interface supplied arrays
      ! experiment specific.

      ! No of internal models in each submodel partition indexed by sm
      !  identifier
      INTEGER :: N_INTERNAL_FOR_SM(SUBMODEL_ID_MAX)

      ! List of  submodel partition identifiers
      INTEGER :: SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION_MAX)

      ! Submodel partition identifier indexed by internal model identifie
      INTEGER :: SUBMODEL_PARTITION_INDEX(INTERNAL_ID_MAX)

      ! Sequence number of internal model indexed by internal model
      ! identifier: required to map from id to STASH internal model
      ! sequence
      INTEGER :: INTERNAL_MODEL_INDEX(INTERNAL_ID_MAX)


      ! Last internal model within a submodel partition if .TRUE.,
      ! indexed by internal model id.
      LOGICAL :: LAST_IM_IN_SM(INTERNAL_ID_MAX)

      ! Common block for information in 3. and 4.
      COMMON/SUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,             &
     &  INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM,SUBMODEL_FOR_SM,            &
     &  N_INTERNAL_FOR_SM,SUBMODEL_PARTITION_LIST,                      &
     &  SUBMODEL_PARTITION_INDEX,                                       &
     &  INTERNAL_MODEL_INDEX,                                           &
     &  LAST_IM_IN_SM

!
!  5. Time information specifying coupling frequencies between internal
!     models and submodels, and multipliers, indexed by sequence of
!     internal models and submodels (ie left to right along node tree).
!     {Not required at this release}.
!
! Namelists for information in 5. {Not required at this release}
!
!
!  6. Lists of coupling nodes defining coupling frequencies between
!     internal models and between submodel partitions. (Not defined
!     yet at this release).
!CALL CNODE
!
!  7. Variables dealing with general coupling switches at the control
!     level. {These will require revision at the next release when
!     coupling between internal models is dealt with more generally.
!     Logicals below are set in routine SETGRCTL.}

      ! new internal model next group of timesteps if .true.
      LOGICAL :: new_im

      ! new submodel dump  next group of timesteps if .true.
      LOGICAL :: new_sm

      COMMON/CSUBMGRP/new_im,new_sm

      INTEGER SUBMODEL_IDENT
      COMMON/SUBMODID/SUBMODEL_IDENT
! CSUBMODL end
! TYPSTS starts
! CSUBMODL must be included before this file
!Applicable to all configurations (except MOS variables)
!STASH related variables for describing output requests and space
!management.
!LL
!LL   AUTHOR            Rick Rawlins
!LL
!LL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
!LL VERSION  DATE
!LL   3.2             Code creation for Dynamic allocation
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL   3.5  Apr. 95   Sub-Models project.
!LL                  Dimensioning of various STASH arrays altered in
!LL                  accordance with internal model separation scheme.
!LL                  Arrays PPXREF, INDEX_PPXREF deleted as they are no
!LL                  longer required.
!LL                  S.J.Swarbrick
!LL
!
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS

      INTEGER :: MOS_MASK_LEN         ! Size of bit mask for MOS

      COMMON/DSIZE_AO/  MOS_MASK_LEN

! TYPSTSZ end
!LL  Comdeck: CPPXREF --------------------------------------------------
!LL
!LL  Purpose: Holds PARAMETER definitions to describe the structure of
!LL           each STASHmaster file record plus some valid entries.
!LL
!LL  Author    Dr T Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                  Add a PPXREF record for model number.
!LL  4.0   26/07/95  T.Johns.  Add codes for real/int/log data types.
!LL  3.5   10/3/94   Sub-Models project:
!LL                 List of PPXREF addressing codes augmented, in order
!LL                 to include all of the pre_STASH master information
!LL                 in the new PPXREF file.
!LL                 PPXREF_CODELEN increased to 38.
!LL                 PPXREF_IDLEN deleted - no longer relevant.
!LL                   S.J.Swarbrick
!LL  4.1   June 96  Wave model parameters included.
!LL                 ppx_ address parameters adjusted to allow for
!LL                  reading option code as 4x5 digit groups.
!LL                   S.J.Swarbrick
!LL  5.0   29/06/99  Add halo type parameter for new dynamics.
!LL                  New grid codes for LAM boundary conditions
!LL                  D.M. Goddard
!LL  5.1   07/03/00  Fixed/Free format conversion
!LL  5.2   19/09/00  Added ppx_atm_lbc_orog descriptor   P.Burton
!LL  5.3   21/08/01  Added ocean lbc descriptors.   M. J. Bell
!LL  5.3   23/07/01  Add valid pp_lbvc codes referenced in UM. R Rawlins
!LL  5.5   30/01/03  Option code increase from 20 to 30 digits thus
!LL                  requiring option code address range increase by
!LL                  2 so all subsequent addressing codes need to be
!LL                  increased by 2 to make a gap.
!LL                  W Roseblade
!LL
!LL  Logical components covered: C40
!LL
!-----------------------------------------------------------------------
! Primary file record definition
      ! length of ID in a record
      Integer, Parameter :: PPXREF_IDLEN      = 2

      ! total length of characters *WARNING* must be multiple of 4
      ! to avoid overwriting
      Integer, Parameter :: PPXREF_CHARLEN    = 36

      ! number of packing profiles
      Integer, Parameter :: PPXREF_PACK_PROFS = 10

      ! total length of codes = no. of codes (excluding profs)
      ! + pack_profs
      Integer, Parameter :: PPXREF_CODELEN    = 33 + PPXREF_PACK_PROFS

! Derived file record sizes
      ! Assume that an integer is at least 4 bytes long. Wastes some
      ! space on an 8 byte machine.
      ! ppx_charword = 9.
      Integer, Parameter :: PPX_CHARWORD      = ((PPXREF_CHARLEN+3)/4)

      ! read buffer record length
      Integer, Parameter :: PPX_RECORDLEN = PPX_CHARWORD+PPXREF_CODELEN
!
!-----------------------------------------------------------------------
! Addressing codes within PPXREF
      Integer, Parameter ::  ppx_model_number   = 1  ! Model number
                                                     ! address
      Integer, Parameter ::  ppx_section_number = 2  ! Section number
                                                     ! address
      Integer, Parameter ::  ppx_item_number    = 3  ! Item number
                                                     ! address
      Integer, Parameter ::  ppx_version_mask   = 4  ! Version mask
                                                     ! address
      Integer, Parameter ::  ppx_space_code     = 5  ! Space code
                                                     ! address
      Integer, Parameter ::  ppx_timavail_code  = 6  ! Time availability
                                                     !  code  address
      Integer, Parameter ::  ppx_grid_type      = 7  ! Grid type code
                                                     ! address
      Integer, Parameter ::  ppx_lv_code        = 8  ! Level type code
                                                     ! address
      Integer, Parameter ::  ppx_lb_code        = 9  ! First level code
                                                     !  address
      Integer, Parameter ::  ppx_lt_code        =10  ! Last level code
                                                     ! address
      Integer, Parameter ::  ppx_lev_flag       =11  ! Level compression
                                                     !  flag  address
      Integer, Parameter ::  ppx_opt_code       =12  ! Sectional option
                                                     ! code  address
      Integer, Parameter ::  ppx_pt_code        =18  ! Pseudo dimension
                                                     ! type  address
      Integer, Parameter ::  ppx_pf_code        =19  ! First pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_pl_code        =20  ! Last pseudo dim
                                                     ! code  address
      Integer, Parameter ::  ppx_ptr_code       =21  ! Section 0 point-
                                                     ! back code address
      Integer, Parameter ::  ppx_dump_packing   =22  ! Dump packing code
                                                     ! address
      Integer, Parameter ::  ppx_lbvc_code      =23  ! PP LBVC code
                                                     ! address
      Integer, Parameter ::  ppx_rotate_code    =24  ! Rotation code
                                                     ! address
      Integer, Parameter ::  ppx_field_code     =25  ! PP field code
                                                     ! address
      Integer, Parameter ::  ppx_user_code      =26  ! User code address
      Integer, Parameter ::  ppx_meto8_levelcode=27  ! CF level code
                                                     ! address
      Integer, Parameter ::  ppx_meto8_fieldcode=28  ! CF field code
                                                     ! address
      Integer, Parameter ::  ppx_cf_levelcode   =27
      Integer, Parameter ::  ppx_cf_fieldcode   =28
      Integer, Parameter ::  ppx_base_level     =29  ! Base level code
                                                     ! address
      Integer, Parameter ::  ppx_top_level      =30  ! Top level code
                                                     ! address
      Integer, Parameter ::  ppx_ref_lbvc_code  =31  ! Ref level LBVC
                                                     ! code address
      Integer, Parameter ::  ppx_data_type      =32  ! Data type code
                                                     ! address
      Integer, Parameter ::  ppx_halo_type      =33
      Integer, Parameter ::  ppx_packing_acc    =34  ! Packing accuracy
                                                     ! code  address
      Integer, Parameter ::  ppx_pack_acc       =34  ! Must be last:


                                                 ! multiple pack_acc to
                                                 ! fill up remaining
                                                 ! array elements


!-------------------------------------------------------------------
! Valid grid type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_atm_nonstd=0      ! Non-standard atmos
                                                  ! grid
      Integer, Parameter :: ppx_atm_tall=1        ! All T points (atmos)
      Integer, Parameter :: ppx_atm_tland=2       ! Land-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tsea=3        ! Sea-only T points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_tzonal=4      ! Zonal field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_tmerid=5      ! Merid field at T
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_uall=11       ! All u points (atmos)
      Integer, Parameter :: ppx_atm_uland=12      ! Land-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_usea=13       ! Sea-only u points
                                                  ! (atmos)
      Integer, Parameter :: ppx_atm_uzonal=14     ! Zonal field at u
                                                  ! points  (atmos)
      Integer, Parameter :: ppx_atm_umerid=15     ! Merid field at u
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_scalar=17     ! Scalar (atmos)
      Integer, Parameter :: ppx_atm_cuall=18      ! All C-grid (u)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_cvall=19      ! All C-grid (v)
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_compressed=21 ! Compressed land
                                                  ! points (atmos)
      Integer, Parameter :: ppx_atm_ozone=22      ! Field on ozone
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_river=23      ! River routing
                                                  ! grid (atmos)
      Integer, Parameter :: ppx_atm_rim=25        ! Rim type field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_theta=26  ! All T points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_u=27      ! All u points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_v=28      ! All v points
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_atm_lbc_orog=29   ! Orography field
                                                  ! (LAM BCs atmos)
      Integer, Parameter :: ppx_ocn_nonstd=30     ! Non-standard ocean
                                                  ! grid
      Integer, Parameter :: ppx_ocn_tcomp=31      ! Compressed T points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_ucomp=32      ! Compressed u points
                                                  !  (ocean)
      Integer, Parameter :: ppx_ocn_tall=36       ! All T points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_uall=37       ! All u points incl.
                                                  ! cyclic  (ocean)
      Integer, Parameter :: ppx_ocn_cuall=38      ! All C-grid (u)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_cvall=39      ! All C-grid (v)
                                                  ! points (ocean)
      Integer, Parameter :: ppx_ocn_tfield=41     ! All non-cyclic T
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_ufield=42     ! All non-cyclic u
                                                  ! points  (ocean)
      Integer, Parameter :: ppx_ocn_tzonal=43     ! Zonal n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_uzonal=44     ! Zonal n-c field at
                                                  ! u points (ocean)
      Integer, Parameter :: ppx_ocn_tmerid=45     ! Merid n-c field at
                                                  ! T points  (ocean)
      Integer, Parameter :: ppx_ocn_umerid=46     ! Merid n-c field at
                                                  ! u points  (ocean)
      Integer, Parameter :: ppx_ocn_scalar=47     ! Scalar (ocean)
      Integer, Parameter :: ppx_ocn_rim=51        ! Rim type field
                                                  ! (LAM BCs ocean)
      Integer, Parameter :: ppx_ocn_lbc_theta=52  ! Ocean rim fields
      Integer, Parameter :: ppx_ocn_lbc_u=53      ! on T & U grids
      Integer, Parameter :: ppx_wam_all=60        ! All points (wave
                                                  ! model)
      Integer, Parameter :: ppx_wam_sea=62        ! Sea points only
                                                  ! (wave model)
      Integer, Parameter :: ppx_wam_rim=65        ! Rim type field
                                                  ! (LAM BCs wave)

!--------------------------------------------------------------------
! Valid rotation type codes
!--------------------------------------------------------------------
      Integer, Parameter :: ppx_unrotated=0       ! Unrotated output
                                                  ! field
      Integer, Parameter :: ppx_elf_rotated=1     ! Rotated ELF field

!-------------------------------------------------------------------
! Valid level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_full_level=1      ! Model full level
      Integer, Parameter :: ppx_half_level=2      ! Model half level
      Integer, Parameter :: ppx_rho_level=1       ! Model rho level
      Integer, Parameter :: ppx_theta_level=2     ! Model theta level
      Integer, Parameter :: ppx_single_level=5    ! Model single level
      Integer, Parameter :: ppx_soil_level=6      ! Deep Soil level

!-------------------------------------------------------------------
! Valid data type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_type_real=1       ! Real data type
      Integer, Parameter :: ppx_type_int=2        ! Integer data type
      Integer, Parameter :: ppx_type_log=3        ! Logical data type

!-------------------------------------------------------------------
! Valid meto8 level type codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_meto8_surf=9999   ! MetO8 surface type
                                                  ! code

!-------------------------------------------------------------------
! Valid dump packing codes
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_pack_off=0        ! Field not packed
                                                  ! (ie. 64 bit)
      Integer, Parameter :: ppx_pack_32=-1        ! Field packed to
                                                  ! 32 bit in  dump
      Integer, Parameter :: ppx_pack_wgdos=1      ! Field packed by
                                                  ! WGDOS method
      Integer, Parameter :: ppx_pack_cfi1=11      ! Field packed using
                                                  ! CFI1  (ocean)

!-------------------------------------------------------------------
! Add valid lbvc codes referenced in model (pp header output labels)
!-------------------------------------------------------------------
      Integer, Parameter :: ppx_lbvc_height  =  1 ! height
      Integer, Parameter :: ppx_lbvc_depth   =  2 ! depth (ocean)
      Integer, Parameter :: ppx_lbvc_pressure=  8 ! pressure
      Integer, Parameter :: ppx_lbvc_theta   = 19 ! potential T
      Integer, Parameter :: ppx_lbvc_hybrid  = 65 ! hybrid height(atmos)
      Integer, Parameter :: ppx_lbvc_PV      = 82 ! potential vorticity
      Integer, Parameter :: ppx_lbvc_surface =129 ! surface
! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      INTEGER :: MOS_OUTPUT_LENGTH
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP,MOS_OUTPUT_LENGTH

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
      INTEGER:: MOS_MASK(MOS_MASK_LEN)
! TYPSTS end
!
      External BL_TRMIX_DD

! DEPENDS ON: bl_trmix_dd
      Call BL_TRMIX_DD(                                                 &
     &        row_length, rows, bl_levels, off_x, off_y, at_extremity   &
     &,       dtrdz_charney_grid,r_rho_levels,r_theta_levels            &
     &,       halo_i,halo_j                                             &
     &,       tr_levels, tr_vars, model_levels                          &
     &,       DIM_CS2                                                   &
     &,       alpha_cd, rhokh_mix                                       &
     &,       p_star, p, timestep                                       &
! Control logicals
     &,       L_MURK,L_MURK_ADVECT,L_BL_TRACER_MIX,L_DUST,L_CAM_DUST    &
     &,       L_SULPC_SO2,l_sulpc_nh3, l_sulpc_dms, l_soot, l_biomass   &
     &,       l_ocff, l_co2_interactive                                 &
! rml 2/7/13 added flag for putting co2 flux into passive tracer
     &,       l_co2_emits, L_CO2_TRACER, L_USE_CARIOLLE                 &
! Fields to mix
     &,       murk, free_tracers, OZONE_TRACER                          &
! Mineral Dust
     &,       DUST_DIV1,DUST_DIV2,DUST_DIV3                             &
     &,       DUST_DIV4,DUST_DIV5,DUST_DIV6                             &
! Sulphur cycle
     &,       so2, dms, so4_aitken, so4_accu, so4_diss, nh3             &
! Soot cycle
     &,       soot_new, soot_aged, soot_cld                             &
! Biomass aerosol
     &,       bmass_new, bmass_agd, bmass_cld                           &
! Fossil-fuel organic carbon aerosol
     &,       ocff_new, ocff_agd, ocff_cld                              &
! Carbon cycle
     &,       co2, co2_emits, co2flux, npp, resp_s                      &
     &,       co2_flux_tot, land_co2                                    &
! Tracer fluxes - kdcorbin, 05/10
     &,       tracer_flux1, tracer_flux2, tracer_flux3, tracer_flux4    &
     &,       tracer_flux5, tracer_flux6, tracer_flux7, tracer_flux8    &
     &,       tracer_flux9, tracer_flux10,tracer_flux11,tracer_flux12   &
     &,       tracer_flux13,tracer_flux14,tracer_flux15,tracer_flux16   &
     &,       tracer_flux17,tracer_flux18,tracer_flux19,tracer_flux20   &
! Emissions fields
     &,       DUST_FLUX                                                 &
     &,       so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em    &
     &,       ocff_hilem, ocff_em, drydep_str                           &
     &,       kent, we_lim, t_frac, zrzi                                &
     &,       kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                &
     &,       zh, zhsc, z_half                                          &
! IN for dry deposition of tracers
     &,       rho_aresist, aresist, resist_b                            &
     &,       R_B_DUST,TSTAR                                            &
     &,       land_points, land_index, ice_frac                         &
! IN variables for MOSES II
     &,       ntype, ntiles, tile_pts, tile_index, tile_frac            &
     &,       canopy, catch, snow_tile, gc                              &
     &,       aresist_tile, resist_b_tile, FLANDG                       &
     &,                                                                 &
! STASH related variables
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &        stashwork                                                 &
     &        )
!
      RETURN
      END SUBROUTINE BL_TRACER_INTCTL
!
