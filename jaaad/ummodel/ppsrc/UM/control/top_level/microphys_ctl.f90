

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


! Subroutine microphys_ctl

      Subroutine microphys_ctl (                                        &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &

! model dimensions.
     &, row_length, rows, rhc_row_length, rhc_rows, n_rows, land_points &
     &, model_levels, wet_model_levels, bl_levels                       &
     &, lspice_dim1,lspice_dim2,lspice_dim3                             &
     &, dst_levels, dsm_levels, Ozone_levels, cloud_levels              &
     &, land_ice_points, soil_points                                    &
     &, n_cca_levels                                                    &
     &, salt_dim1, salt_dim2, salt_dim3                                 &

! Model switches
     &, model_domain, L_CAL360, Ltimer, L_RHCPT, L_Murk, L_DUST         &
     &, l_sulpc_so2, l_sulpc_nh3, l_soot, l_biomass, l_ocff             &
     &, L_use_sulphate_autoconv,  L_auto_debias                         &
     &, L_seasalt_CCN, L_soot_CCN, L_bmass_CCN, L_ocff_CCN              &
     &, L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup, L_it_melting             &
     &, l_eacf, l_mixing_ratio, l_cry_agg_dep, l_droplet_settle         &
     &, L_pc2, L_use_biogenic                                           &

! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &

! in coordinate information
     &, r_rho_levels, r_theta_levels, r_at_u, r_at_v                    &
     &, eta_theta_levels, eta_rho_levels                                &
     &, delta_lambda, delta_phi                                         &
     &, lat_rot_NP, long_rot_NP                                         &

! in time stepping information.
     &, timestep, radiation_timestep                                    &
     &, radiation_tstep_diag                                            &
     &, radiation_tstep_prog                                            &
     &, val_year, val_day_number, val_hour, val_minute                  &
     &, val_second, timestep_number                                     &

! Model parameters
     &, rhcrit, cw_sea, cw_land                                         &
     &, cloud_fraction_method,overlap_ice_liquid                        &
     &, ice_fraction_method,ctt_weight,t_weight                         &
     &, qsat_fixed,sub_cld                                              &
     &, L_seq_mcr,L_autoc_3b,L_autolim_3b,L_autoconv_murk               &
     &, ec_auto,N_drop_land                                             &
     &, N_drop_sea,N_drop_land_cr,N_drop_sea_cr,Ntot_land, Ntot_sea     &
     &, x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic                        &
     &, lsp_ei,lsp_fi,lsp_eic,lsp_fic                                   &

! Primary fields passed in
     &, T_n, q_n, qcl_n, qcf_n, qcf2_n, qrain_n, qgraup_n               &
     &, cf_n, cfl_n, cff_n                                              &
     &, snow_depth                                                      &
     &, land_sea_mask, ice_fract                                        &
     &, p_layer_centres, p_layer_boundaries, exner_theta_levels         &
     &, rho_r2                                                          &
     &, AEROSOL                                                         &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, SO2, NH3, SO4_AITKEN, SO4_ACCU, SO4_DISS                        &
     &, aged_soot, cloud_soot, fresh_bmass, aged_bmass, cloud_bmass     &
     &, aged_ocff, cloud_ocff, biogenic                                 &

! Other fields passed in
     &, energy_correction, cca, ccb, cct, p, p_star, p_theta_levels     &
     &, exner_rho_levels                                                &
     &, ntml, cumulus                                                   &
     &, fland, land_index                                               &
     &, u_1, v_1                                                        &

! Variables for stochastic physics random parameters
     &, M_CI                                                            &

! diagnostic info
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &  STASHwork4                                                      &
!
! SCM diagnostics switches (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! Increment fields passed in/out
     &, T_inc, q_inc, qcl_inc, qcf_inc                                  &
     &, qcf2_inc, qrain_inc, qgraup_inc                                 &
     &, cf_inc, cfl_inc, cff_inc                                        &

! Fields required elsewhere
     &, ls_rain, ls_snow, micro_tends                                   &

! Section information
     &, maxsects, h_sect                                                &

! error information
     &, Error_code  )

! purpose: Interface to Atmospheric Physics parametrizations.
!          This version interfaces to cloud scheme, boundary layer,
!          hydrology, large scale rain, radiation, convection.
!
! but Not:
!          gravity wave drag, and optionally energy
!          correction. Diagnostics also missing.
!
!          IN/OUT etc intents to be added later.
!
! code description:
!   language: fortran 77 + cray extensions
!   this code is written to umdp3 programming standards.
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
!C_DUSTSCAV.............................................................
! Description: Contains  mineral dust scavenging coefficients
! Current Code Owner: Stephanie Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
!  1      12/02/03  Original Code.   Stephanie Woodward
!
! Mineral dust scavenging coefficients
!
      REAL :: KRAIN_DUST(NDIV)  ! scav. coeff. for rain
      REAL :: KSNOW_DUST(NDIV)  ! scav. coeff. for snow
!
      DATA KRAIN_DUST/2.e-5,2.e-5,3.e-5,6.e-5,4.e-4,4.e-4/
      DATA KSNOW_DUST/2.e-5,2.e-5,3.e-5,6.e-5,4.e-4,4.e-4/
!.......................................................................
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------


! arguments with intent in. ie: input variables.

! Parallel setup variables
      Integer                                                           &
     &  halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, global_row_length                                               &
                           ! number of points on a row
     &, global_rows                                                     &
                           ! NUMBER OF global rows
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbours
                             ! in the horizontal plane
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)                                      &
     &, me         ! My processor number

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
!                          south, east or west of the processor grid
!                          (array index PNorth etc. from parparm.h)

      Logical                                                           &
     &  L_murk                                                          &
                    ! true for aerosol calculations
     &, L_DUST                                                          &
               ! true for mineral dust calcs
     &, l_sulpc_so2                                                     &
                    ! true for sulphur cycle calcs
     &, l_sulpc_nh3                                                     &
                    ! true for sulphur/nh3 calcs
     &, l_soot                                                          &
                    ! true for soot cycles calcs
     &, l_biomass                                                       &
                    ! true for biomass aerosol calcs
     &, l_ocff                                                          &
                    ! true for fossil-fuel organic carbon calcs
     &, l_eacf                                                          &
                    ! true for empirically adjusted cloud frac.
     &, l_mixing_ratio                                                  &
                    ! true for mixing ratio formulation
     &, l_cry_agg_dep                                                   &
                    ! Limit the supersaturation that can be
                    ! removed by deposition depending on
                    ! amount of ice in each category
     &, l_droplet_settle 
                    ! Allow cloud droplets to settle

! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, rhc_row_length                                                  &
                       ! = row_length if L_RHCPT true, 1 otherwise.
     &, rhc_rows                                                        &
                       ! = rows       if L_RHCPT true, 1 otherwise.
     &, n_rows                                                          &
     &, land_points                                                     &
                    ! IN No.of land points being processed, can be 0.
     &, model_levels                                                    &
     &, wet_model_levels                                                &
     &, bl_levels                                                       &
     &, dst_levels                                                      &
                    ! number of deep soil temperature levels
     &, dsm_levels                                                      &
                    ! number of deep soil moisture levels
     &, Ozone_levels                                                    &
     &, cloud_levels                                                    &
     &, land_ice_points                                                 &
                        ! number of land ice points
     &, soil_points                                                     &
                        ! number of soil points
     &, n_cca_levels                                                    &
                      ! Number of levels for conv cloud
                      ! amount: 1 for 2D, nlevs for 3D.
     &, lspice_dim1                                                     &
                      ! Dimensions for 3D diagnostic arrays.
     &, lspice_dim2                                                     &
                      ! These are set to 1 in order to save
     &, lspice_dim3                                                     &
                      ! memory if the diagnostics are not used.
     &, salt_dim1                                                       &
                      !
     &, salt_dim2                                                       &
                      ! Dimensions for sea-salt aerosol arrays
     &, salt_dim3     !

! Model switches
      Integer                                                           &
     &  model_domain

      Logical                                                           &
     &  Ltimer                                                          &
                 ! true then output some timing information
     &, L_CAL360                                                        &
                        ! true if using 360 day calender
     &, L_RHCPT                                                         &
                 ! Switch for 3D diagnosed RHcrit not 1D parameter
     &, L_use_sulphate_autoconv                                         &
                                ! Switch for sulphate aerosol for CCN
     &, L_seasalt_CCN                                                   &
                      ! Switch for sea-salt aerosol param. for CCN
     &, L_soot_CCN                                                      &
                    ! Switch for soot aerosol for CCN
     &, L_bmass_CCN                                                     &
                    ! Switch for biomass aerosol for CCN
     &, L_ocff_CCN                                                      &
                 ! Switch for fossil-fuel organic carbon aerosol for CCN
     &, L_use_biogenic                                                  &
                        ! Switch biogenic aerosol for CCN
     &, L_auto_debias                                                   &
                      ! Switch for autoconversion de-biasing scheme.
     &, L_pc2                                                           &
                 ! Use PC2 cloud scheme
     &, L_mcr_qcf2                                                      &
                     ! Switch for second ice variable
     &, L_mcr_qrain                                                     &
                     ! Switch for prognostic rain
     &, L_mcr_qgraup                                                    &
                     ! Switch for prognostic graupel
     &, L_psd                                                           &
                     ! Use generic ice particle size distribution
     &, L_it_melting
                     ! Use iterative melting

!  Additional arguments for SCM diagnostics which are dummy in full UM
      Integer                                                           &
     &  nSCMDpkgs             ! No of SCM diagnostics packages

      Logical                                                           &
     &  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

! physical constants
      Real                                                              &
     &  lc, lf, cp                                                      &
     &, two_Omega                                                       & 
                        ! twice Earth's rotation rate
     &, p_zero                                                          &
     &, kappa                                                           &
     &, R, g, Lapse_Rate, earth_radius, Pi

! Model parameters
      Real                                                              &
     &  cw_sea                                                          &
                        ! IN threshold cloud liquid water content
                        !    over sea for conversion to ppn
                        !   (kg water per m**3)
     &, cw_land                                                         &
                        ! IN threshold cloud liquid water content
                        !    over land for conversion to ppn
                        !    (kg water per m**3)
     &, rhcrit(wet_model_levels)  ! IN Critical relative humidity.
                                  ! the values need to be tuned
                                  ! for the given set of levels.
! Variable used in stochastic physics random parameters
      Real                                                              &
     & M_CI            ! Used to modify ice fall speed 

! time information for current timestep
      Real                                                              &
     &  timestep                                                        &
     &, radiation_timestep                                              &
     &, radiation_tstep_diag                                            &
     &, radiation_tstep_prog

      Integer                                                           &
     &  val_year                                                        &
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second                                                      &
     &, timestep_number

! Primary fields passed in
      Real, Intent (InOut) ::                                           &
     &  T_n(row_length, rows, model_levels)                             &
     &, q_n(row_length, rows, wet_model_levels)                         &
     &, qcl_n(row_length, rows, wet_model_levels)                       &
     &, qcf_n(row_length, rows, wet_model_levels)                       &
     &, qcf2_n(row_length, rows, wet_model_levels)                      &
     &, qrain_n(row_length, rows, wet_model_levels)                     &
     &, qgraup_n(row_length, rows, wet_model_levels)                    &
     &, cf_n(row_length, rows, wet_model_levels)                        &
     &, cfl_n(row_length, rows, wet_model_levels)                       &
     &, cff_n(row_length, rows, wet_model_levels)                       &
     &, snow_depth(row_length, rows)                                    &
     &, ice_fract(row_length, rows)                                     &
     &, fland(land_points)                                              &
     &, u_1(salt_dim1, salt_dim2)                                       &
                                   ! Lowest layer u-wind component
     &, v_1(salt_dim1, salt_dim2)  ! Lowest layer v-wind component

      Integer                                                           &
     &  land_index(land_points)

      logical                                                           &
     &  land_sea_mask(row_length, rows)

      Real                                                              &
     &  energy_correction

      Real                                                              &
     &  p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, p_layer_centres(row_length, rows, 0:model_levels)               &
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = pstar, and at top = 0.
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y, model_levels)            &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                   1-off_y:rows+off_y, model_levels)              &
     &, rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &          model_levels)                                           &
                               ! Density*radius^2
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, p_theta_levels(1-off_x:row_length+off_x,                        &
     &                   1-off_y:rows+off_y, model_levels)              &
     &, p_star(row_length, rows)

      Integer                                                           &
     &  ntml(row_length, rows)   ! Height of diagnosed BL top

      Logical                                                           &
     &  cumulus(row_length, rows)  ! Logical indicator for convection

      Real, Intent(InOut) ::                                            &
                                     ! tracer variables
     &  aerosol   ( 1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, DUST_DIV1(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div1
     &, DUST_DIV2(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div2
     &, DUST_DIV3(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div3
     &, DUST_DIV4(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div4
     &, DUST_DIV5(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div5
     &, DUST_DIV6(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &       MODEL_LEVELS)                                              &
                           !dust mmr in div6
     &, so2       ( 1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, nh3       ( 1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, so4_aitken( 1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, so4_accu  ( 1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, so4_diss  ( 1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, aged_soot ( 1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, cloud_soot( 1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, fresh_bmass(1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, aged_bmass( 1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, cloud_bmass(1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, aged_ocff(  1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )            &
     &, cloud_ocff( 1 - off_x : row_length + off_x,                     &
     &              1 - off_y : rows + off_y, model_levels )

! Convection
      Real                                                              &
     &  cca (row_length, rows, n_cca_levels)

      Integer                                                           &
     &  ccb (row_length, rows)                                          &
     &, cct (row_length, rows)

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
! Diagnostics info
      REAL                                                              &
     & STASHwork4(*)     ! STASH workspace

!
! Co-ordinate arrays
      Real                                                              &
           ! local vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
     &  r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)

      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP

      Logical                                                           &
                        !, Intent(IN)
     &    L_seq_mcr                                                     &
                               ! Use sequential updating of mphys
     &   ,L_autoc_3b                                                    &
                               ! Use 3B autoconversion method
     &   ,L_autolim_3b                                                  &
                               ! Use fixed 3B values for the
                               ! autoconversion limit
     &   ,L_autoconv_murk      ! Use murk aerosol to calc. drop number

      Integer                                                           &
                        !, Intent(IN)
     &    cloud_fraction_method                                         &
                                 ! Method for calculating
                               ! total cloud fraction
     &   ,ice_fraction_method  ! Method for calculating ice cloud frac.

      Real                                                              &
                        !, Intent(IN)
     &    overlap_ice_liquid                                            &
                               ! Overlap between ice and liquid phases
     &   ,ctt_weight                                                    &
                               ! Weighting of cloud top temperature
     &   ,t_weight                                                      &
                               ! Weighting of local temperature
     &   ,qsat_fixed                                                    &
                               ! Fixed value of saturation humidity
     &   ,sub_cld              ! Scaling parameter

      Real                                                              &
                        !, Intent(IN)
     &    ec_auto                                                       &
                               ! Collision coalescence efficiency
     &   ,N_drop_land                                                   &
                               ! Number of droplets over land / m-3
     &   ,N_drop_sea                                                    &
                               ! Number of droplets over sea / m-3
     &   ,N_drop_land_cr                                                &
                               ! N_drop_land ^ (-1/3) / m
     &   ,N_drop_sea_cr                                                 &
                               ! N_drop_sea ^ (-1/3) / m
     &   ,Ntot_land                                                     &
                               ! Number of droplets over land / m-3
     &   ,Ntot_sea             ! Number of droplets over sea / m-3

      Real                                                              &
                        !, Intent(IN)
     &    x1i                                                           &
                    ! Intercept of aggregate size distribution
     &   ,x1ic                                                          &
                    ! Intercept of crystal size distribution
     &   ,x1r                                                           &
                    ! Intercept of raindrop size distribution
     &   ,x2r                                                           &
                    ! Scaling parameter of raindrop size distribn
     &   ,x4r                                                           &
                    ! Shape parameter of raindrop size distribution
     &   ,ai,  bi                                                       &
                    ! Ice aggregate mass-size relationship m(D)=AI D^BI
     &   ,aic, bic                                                      &
                    ! Ice crystal mass-size relationship m(D)=AIC D^BIC
     &   ,lsp_ei,  lsp_fi                                               &
                    ! Ice aggregate Best number and Reynolds number 
                    ! relationship: Re(D) =LSP_EI Be^LSP_FI
     &   ,lsp_eic, lsp_fic                                                      
                    ! Ice crystal Best number and Reynolds number 
                    ! relationship: Re(D) =LSP_EIC Be^LSP_FIC

! arguments with intent in/out. ie: input variables changed on output.
      Real, Intent (InOut) ::                                           &
     &  T_inc(row_length, rows, model_levels)                           &
     &, q_inc(row_length, rows, wet_model_levels)                       &
     &, qcl_inc(row_length, rows, wet_model_levels)                     &
     &, qcf_inc(row_length, rows, wet_model_levels)                     &
     &, qcf2_inc(row_length, rows, wet_model_levels)                    &
     &, qrain_inc(row_length, rows, wet_model_levels)                   &
     &, qgraup_inc(row_length, rows, wet_model_levels)                  &
     &, cf_inc(row_length, rows, wet_model_levels)                      &
     &, cfl_inc(row_length, rows, wet_model_levels)                     &
     &, cff_inc(row_length, rows, wet_model_levels)

      Integer                                                           &
     &  Error_code

! arguments with intent out. ie: output variables.
      Real                                                              &
     &  ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, micro_tends(row_length, rows, bl_levels, 2)
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)

      Integer                                                           &
     &   maxsects

      Character*3                                                       &
     &   h_sect(0:maxsects)

! local variables.
! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

      Character*(*), Parameter ::  RoutineName = 'microphys_ctl'
!
! loop counters
      Integer                                                           &
     &  i, j, k                                                         &
     &, IDIV !loop counter for dust divisions

! Diagnostic switches

! Local variables
      Logical                                                           &
     &  ext_LAND(0:rhc_row_length+1,0:rhc_rows+1)
!
      Real                                                              &
     & land_fract(row_length, rows)
      Real                                                              &
     &  RHCPT(rhc_row_length, rhc_rows, wet_model_levels)               &
     &, LS_RAIN3D(LSPICE_DIM1,LSPICE_DIM2,LSPICE_DIM3)                  &
! rainfall rate out of each model level for diagonstic
     &, LS_SNOW3D(LSPICE_DIM1,LSPICE_DIM2,LSPICE_DIM3)                  &
! snowfall rate out of each model level for diagonstic
     &, RAINFRAC3D(LSPICE_DIM1,LSPICE_DIM2,LSPICE_DIM3)                 &
! Rain fraction for diagnostic
     &, RNOUT_TRACER(LSPICE_DIM1,LSPICE_DIM2)                           &
! Total tracer amount scavenged by rainout (kg m-2)
     &, LSCAV_TR(LSPICE_DIM1,LSPICE_DIM2)                               &
! Total tracer amount scavenged by washout (kg m-2)
     &, rnout_soot(LSPICE_DIM1,LSPICE_DIM2)                             &
! Total soot amount scavenged by rainout (kg m-2)
     &, lscav_soot(LSPICE_DIM1,LSPICE_DIM2)                             &
! Total soot amount scavenged by washout (kg m-2)
     &, rnout_bmass(LSPICE_DIM1,LSPICE_DIM2)                            &
! Total biomass amount scavenged by rainout (kg m-2)
     &, lscav_bmass(LSPICE_DIM1,LSPICE_DIM2)                            &
! Total biomass amount scavenged by washout (kg m-2)
     &, rnout_ocff(LSPICE_DIM1, LSPICE_DIM2)                            &
! Total fossil-fuel organic carbon amount scavenged by rainout (kg m-2)
     &, lscav_ocff(LSPICE_DIM1, LSPICE_DIM2)
! Total fossil-fuel organic carbon amount scavenged by washout (kg m-2)

      Real                                                              &
     &  T_work(row_length, rows, model_levels)                          &
     &, q_work(row_length, rows, wet_model_levels)                      &
     &, qcl_work(row_length, rows, wet_model_levels)                    &
     &, qcf_work(row_length, rows, wet_model_levels)                    &
     &, cf_work(row_length, rows, wet_model_levels)                     &
     &, cfl_work(row_length, rows, wet_model_levels)                    &
     &, cff_work(row_length, rows, wet_model_levels)                    &
     &, ext_p_layer_centres(0:rhc_row_length+1,0:rhc_rows+1,            &
     &                                         0:wet_model_levels)      &
     &, ext_TL(0:rhc_row_length+1, 0:rhc_rows+1, wet_model_levels)      &
     &, ext_QL(0:rhc_row_length+1, 0:rhc_rows+1, wet_model_levels)      &
     &, ext_QCF(0:rhc_row_length+1,0:rhc_rows+1, wet_model_levels)      &
     &, ext_ICE_FRAC(0:rhc_row_length+1,0:rhc_rows+1)                   &
     &, ext_land_fract(0:rhc_row_length+1,0:rhc_rows+1)                 &
     &, sea_salt_film(salt_dim1, salt_dim2, salt_dim3)                  &
!          Film-mode sea-salt aerosol number concentration
     &, sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)                   &
!          Jet-mode sea-salt aerosol number concentration
     &, height(salt_dim1, salt_dim2, salt_dim3)                         &
!          Layer-centre height above surface
     &, biogenic(row_length, rows, model_levels)                        &
!          m.m.r. of biogenic aerosols for CCN
     &, DM(ROW_LENGTH, ROWS, MODEL_LEVELS) ! mass air p.u. area in lev

      ! Local work arrays for additional microphysics fields if in use
      Real, Dimension (:,:,:), Allocatable ::                           &
     &  qcf2_work, qrain_work, qgraup_work
      Real                                                              &
     &  Tinc_np1(row_length, rows)                                      &
     &, qinc_np1(row_length, rows)                                      &
     &, qclinc_np1(row_length, rows)                                    &
     &, qcfinc_np1(row_length, rows)
!
! Local work array for increment diagnostics
      Real                                                              &
     &  work_3d(row_length, rows, model_levels)

! Scavenged tracers (in column) for diagnostics
      Real                                                              &
     &  lscav_so2( row_length, rows )                                   &
     &, lscav_nh3( row_length, rows )                                   &
     &, lscav_so4ait( row_length, rows )                                &
     &, lscav_so4acc( row_length, rows )                                &
     &, lscav_so4dis( row_length, rows )                                &
     &, LSCAV_DUST_ALL(ROW_LENGTH,ROWS,NDIV)                            &
     &, DUST_ALL(ROW_LENGTH,ROWS,MODEL_LEVELS,NDIV)


       REAL                                                             &
     & RAINRATE                                                         &
                 !rate corrected for possible -ive precip
     &,SNOWRATE                                                         &
                 !rate corrected for possible -ive precip
     &,RATE                                                             &
                 !
     &,DELTA_DUST !

! Microphysical process rate diagnostics
! Note: These arrays will only increase memory usage and are
!       only referenced if the particular diagnostic is active
!       (i.e. if the associated logical is .true., set from SF
!       STASH flag.

      ! Local work arrays for additional microphysics fields if in use
      Real, Dimension (:,:,:), Allocatable ::                           &
     &  PSDEP,PSAUT,PSACW,PSACR,PSACI,PSMLT,PSMLTEVP                    &
     &, PRAUT,PRACW,PREVP                                               &
     &, PGAUT,PGACW,PGACS,PGMLT                                         &
     &, PIFRW,PIPRM,PIDEP,PIACW,PIACR,PIMLT,PIMLTEVP                    &
     &, PIFALL,PSFALL,PRFALL,PGFALL,PLSET,PLEVPSET

! Microphysical process rate diagnostic logical switches
      LOGICAL                                                           &
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
      LOGICAL                                                           &
     &  L_PRAUT_diag                                                    &
                       ! Autoconversion of cloud drops to rain
     &, L_PRACW_diag                                                    &
                       ! Accretion of liq. water by rain
     &, L_PREVP_diag   ! Evaporation of rain
      LOGICAL                                                           &
     &  L_PGAUT_diag                                                    &
                       ! Autoconversion of graupel from agg.
     &, L_PGACW_diag                                                    &
                       ! Accretion of liq. water by graupel
     &, L_PGACS_diag                                                    &
                       ! Collection of snow agg. by graupel
     &, L_PGMLT_diag   ! Melting of graupel
      LOGICAL                                                           &
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
      LOGICAL                                                           &
     &  L_PIFALL_diag                                                   &
                       ! Sedimentation of ice crystals
     &, L_PSFALL_diag                                                   &
                       ! Sedimentation of aggregates
     &, L_PRFALL_diag                                                   &
                       ! Sedimentation of rain
     &, L_PGFALL_diag  ! Sedimentation of graupel
      LOGICAL                                                           &
     &  L_PLSET_diag                                                    &
                       ! Droplet settling of liquid water
     &, L_PLEVPSET_diag! Evaporated settled droplets

! External Routines:
      External ls_calc_rhcrit
      External ls_cld
      External ls_ppn
      External timer
      External set_seasalt
      External rainout_intctl
      External sl3dwash
      External nh3dwash
      EXTERNAL MASS_CALC
      External diagnostics_lsrain


! ----------------------------------------------------------------------
! Section Microphysics. Call microphys_ctl routine
! ----------------------------------------------------------------------
! Call timer for cloud code
! DEPENDS ON: timer
      If (Ltimer) Call timer ('LS Cloud',3)

! Store values in work arrays
      Do k = 1, wet_model_levels
        Do j = 1, rows
          Do i = 1, row_length
            If (L_pc2) Then
! No diagnostic cloud scheme is required. Store temperature and
! vapour content in T_work and q_work
              T_work(i,j,k)   = t_n(i,j,k)
              q_work(i,j,k)   = q_n(i,j,k)
            Else  ! L_pc2
! inputs to diagnostic cloud scheme are Tl_star, qT_star
! outputs from diag. cloud scheme are cloud_fraction, T, q, qcl, qcf
! output T is held in theta_star
! Calculate Tl, store in T_work
! Calculate qT, store in q_work
              T_work(i,j,k) = t_n(i,j,k) -                              &
     &                        (lc * qcl_n(i,j,k) ) / Cp
              q_work(i,j,k) = q_n(i,j,k) + qcl_n(i,j,k)
            End If  ! L_pc2
            qcl_work(i,j,k) = qcl_n(i,j,k)
            qcf_work(i,j,k) = qcf_n(i,j,k)
            cf_work(i,j,k)  = cf_n(i,j,k)
            cfl_work(i,j,k) = cfl_n(i,j,k)
            cff_work(i,j,k) = cff_n(i,j,k)
          End Do
        End Do
      End Do

      If (L_mcr_qcf2) Then  ! Second cloud ice variable in use
        Allocate ( qcf2_work(row_length, rows, wet_model_levels) )
        ! Add qcf and qcf2 to get total ice for cloud scheme
        ! qcf_work and qcf2_work are reset after call to ls_cld
        qcf_work(:,:,:) = qcf_n(:,:,:) + qcf2_n(:,:,:)
      Else
        Allocate ( qcf2_work(1,1,1) )
      End If

      If (L_mcr_qrain) Then  ! Prognostic rain in use
        Allocate ( qrain_work(row_length, rows, wet_model_levels) )
        qrain_work(:,:,:) = qrain_n(:,:,:)
      Else
        Allocate ( qrain_work(1,1,1) )
      End If

      If (L_mcr_qgraup) Then  ! Prognostic graupel in use
        Allocate ( qgraup_work(row_length, rows, wet_model_levels) )
        qgraup_work(:,:,:) = qgraup_n(:,:,:)
      Else
        Allocate ( qgraup_work(1,1,1) )
      End If

! Dry level T values (only required for diagnostics):
      Do k = wet_model_levels+1,model_levels
        Do j = 1, rows
          Do i = 1, row_length
            T_work(i,j,k) = t_n(i,j,k)
          End Do
        End Do
      End Do

      L_PIFRW_diag    = SF(240,4)
      L_PIPRM_diag    = SF(241,4)
      L_PIDEP_diag    = SF(243,4)
      L_PSDEP_diag    = SF(245,4)
      L_PIACW_diag    = SF(247,4)
      L_PSACW_diag    = SF(248,4)
      L_PIACR_diag    = SF(249,4)
      L_PSACR_diag    = SF(250,4)
      L_PIMLTEVP_diag = SF(251,4)
      L_PSMLTEVP_diag = SF(252,4)
      L_PIMLT_diag    = SF(253,4)
      L_PSMLT_diag    = SF(254,4)
      L_PSAUT_diag    = SF(255,4)
      L_PSACI_diag    = SF(256,4)
      L_PRAUT_diag    = SF(257,4)
      L_PRACW_diag    = SF(258,4)
      L_PREVP_diag    = SF(259,4)
      L_PGAUT_diag    = SF(260,4)
      L_PGACW_diag    = SF(261,4)
      L_PGACS_diag    = SF(262,4)
      L_PGMLT_diag    = SF(263,4)
      L_PIFALL_diag   = SF(265,4)
      L_PSFALL_diag   = SF(266,4)
      L_PRFALL_diag   = SF(267,4)
      L_PGFALL_diag   = SF(268,4)
      L_PLSET_diag    = SF(269,4)
      L_PLEVPSET_diag = SF(270,4)

      ! Allocate arrays for required microphysics diagnostics
      ! Homogeneous freezing nucl.
      If (L_PIFRW_diag) Then
        Allocate ( PIFRW(row_length,rows,wet_model_levels) )
        PIFRW(:,:,:) = 0.0
      Else
        Allocate ( PIFRW(1,1,1) )
      End If

      ! Heterogeneous nucl.
      If (L_PIPRM_diag) Then
        Allocate ( PIPRM(row_length,rows,wet_model_levels) )
        PIPRM(:,:,:) = 0.0
      Else
        Allocate ( PIPRM(1,1,1) )
      End If

      ! Deposition of vapour to ice
      If (L_PIDEP_diag) Then
        Allocate ( PIDEP(row_length,rows,wet_model_levels) )
        PIDEP(:,:,:) = 0.0
      Else
        Allocate ( PIDEP(1,1,1) )
      End If

      ! Deposition of vapour to snow
      If (L_PSDEP_diag) Then
        Allocate ( PSDEP(row_length,rows,wet_model_levels) )
        PSDEP(:,:,:) = 0.0
      Else
        Allocate ( PSDEP(1,1,1) )
      End If

      ! Accretion of liq. by ice
      If (L_PIACW_diag) Then
        Allocate ( PIACW(row_length,rows,wet_model_levels) )
        PIACW(:,:,:) = 0.0
      Else
        Allocate ( PIACW(1,1,1) )
      End If

      ! Accretion of liq. by snow
      If (L_PSACW_diag) Then
        Allocate ( PSACW(row_length,rows,wet_model_levels) )
        PSACW(:,:,:) = 0.0
      Else
        Allocate ( PSACW(1,1,1) )
      End If

      ! Collection of rain by ice
      If (L_PIACR_diag) Then
        Allocate ( PIACR(row_length,rows,wet_model_levels) )
        PIACR(:,:,:) = 0.0
      Else
        Allocate ( PIACR(1,1,1) )
      End If

      ! Collection of rain by snow
      If (L_PSACR_diag) Then
        Allocate ( PSACR(row_length,rows,wet_model_levels) )
        PSACR(:,:,:) = 0.0
      Else
        Allocate ( PSACR(1,1,1) )
      End If

      ! Evaporation of melting ice
      If (L_PIMLTEVP_diag) Then
        Allocate ( PIMLTEVP(row_length,rows,wet_model_levels) )
        PIMLTEVP(:,:,:) = 0.0
      Else
        Allocate ( PIMLTEVP(1,1,1) )
      End If

      ! Evap. of melting aggregates
      If (L_PSMLTEVP_diag) Then
        Allocate ( PSMLTEVP(row_length,rows,wet_model_levels) )
        PSMLTEVP(:,:,:) = 0.0
      Else
        Allocate ( PSMLTEVP(1,1,1) )
      End If

      ! Melting of ice crystals
      If (L_PIMLT_diag) Then
        Allocate ( PIMLT(row_length,rows,wet_model_levels) )
        PIMLT(:,:,:) = 0.0
      Else
        Allocate ( PIMLT(1,1,1) )
      End If

      ! Melting of snow aggregates
      If (L_PSMLT_diag) Then
        Allocate ( PSMLT(row_length,rows,wet_model_levels) )
        PSMLT(:,:,:) = 0.0
      Else
        Allocate ( PSMLT(1,1,1) )
      End If

      ! Autoconversion of snow
      If (L_PSAUT_diag) Then
        Allocate ( PSAUT(row_length,rows,wet_model_levels) )
        PSAUT(:,:,:) = 0.0
      Else
        Allocate ( PSAUT(1,1,1) )
      End If

      ! Collection of ice crystals
      If (L_PSACI_diag) Then
        Allocate ( PSACI(row_length,rows,wet_model_levels) )
        PSACI(:,:,:) = 0.0
      Else
        Allocate ( PSACI(1,1,1) )
      End If

      ! Autoconversion of cloud
      If (L_PRAUT_diag) Then
        Allocate ( PRAUT(row_length,rows,wet_model_levels) )
        PRAUT(:,:,:) = 0.0
      Else
        Allocate ( PRAUT(1,1,1) )
      End If

      ! Accretion of liq. by rain
      If (L_PRACW_diag) Then
        Allocate ( PRACW(row_length,rows,wet_model_levels) )
        PRACW(:,:,:) = 0.0
      Else
        Allocate ( PRACW(1,1,1) )
      End If

      ! Evaporation of rain
      If (L_PREVP_diag) Then
        Allocate ( PREVP(row_length,rows,wet_model_levels) )
        PREVP(:,:,:) = 0.0
      Else
        Allocate ( PREVP(1,1,1) )
      End If

      ! Autoconversion of graupel
      If (L_PGAUT_diag) Then
        Allocate ( PGAUT(row_length,rows,wet_model_levels) )
        PGAUT(:,:,:) = 0.0
      Else
        Allocate ( PGAUT(1,1,1) )
      End If

      ! Accretion of liq. by graup
      If (L_PGACW_diag) Then
        Allocate ( PGACW(row_length,rows,wet_model_levels) )
        PGACW(:,:,:) = 0.0
      Else
        Allocate ( PGACW(1,1,1) )
      End If

      ! Collection of snow by graup
      If (L_PGACS_diag) Then
        Allocate ( PGACS(row_length,rows,wet_model_levels) )
        PGACS(:,:,:) = 0.0
      Else
        Allocate ( PGACS(1,1,1) )
      End If

      ! Melting of graupel
      If (L_PGMLT_diag) Then
        Allocate ( PGMLT(row_length,rows,wet_model_levels) )
        PGMLT(:,:,:) = 0.0
      Else
        Allocate ( PGMLT(1,1,1) )
      End If

      ! Sedimentation of ice crystals
      If (L_PIFALL_diag) Then
        Allocate ( PIFALL(row_length,rows,wet_model_levels) )
        PIFALL(:,:,:) = 0.0
      Else
        Allocate ( PIFALL(1,1,1) )
      End If

      ! Sedimentation of aggregates
      If (L_PSFALL_diag) Then
        Allocate ( PSFALL(row_length,rows,wet_model_levels) )
        PSFALL(:,:,:) = 0.0
      Else
        Allocate ( PSFALL(1,1,1) )
      End If

      ! Sedimentation of rain
      If (L_PRFALL_diag) Then
        Allocate ( PRFALL(row_length,rows,wet_model_levels) )
        PRFALL(:,:,:) = 0.0
      Else
        Allocate ( PRFALL(1,1,1) )
      End If

      ! Sedimentation of graupel
      If (L_PGFALL_diag) Then
        Allocate ( PGFALL(row_length,rows,wet_model_levels) )
        PGFALL(:,:,:) = 0.0
      Else
        Allocate ( PGFALL(1,1,1) )
      End If

      ! Droplet settling of liquid water
      If (L_PLSET_diag) Then
        Allocate ( PLSET(row_length,rows,wet_model_levels) )
        PLSET(:,:,:) = 0.0
      Else
        Allocate ( PLSET(1,1,1) )
      End If

      ! Evaporated settled droplets
      If (L_PLEVPSET_diag) Then
        Allocate ( PLEVPSET(row_length,rows,wet_model_levels) )
        PLEVPSET(:,:,:) = 0.0
      Else
        Allocate ( PLEVPSET(1,1,1) )
      End If

!
! ----------------------------------------------------------------------
! Section CLD.1.a Calculate diagnostic RHcrit or read as namelist param.
! ----------------------------------------------------------------------
!
! Lrhcpt_if1:
      If (L_RHCPT) Then
!       RHCRIT is 3D diagnosed variable
! Wet_mlev_do1:
        Do k = 1, wet_model_levels
! Rhc_rows_do1:
          Do j = 1, rhc_rows
! Rhc_rowlen_do1:
            Do i = 1, rhc_row_length
              ext_p_layer_centres(i,j,k) = p_layer_centres(i,j,k)
              If (L_pc2) Then
                ext_TL(i,j,k) = T_work(i,j,k)                           &
     &                        - (lc * qcl_work(i,j,k) ) / Cp
                ext_QL(i,j,k) = q_work(i,j,k) + qcl_work(i,j,k)
              Else  ! L_pc2
              ext_TL(i,j,k) = T_work(i,j,k)
              ext_QL(i,j,k) = q_work(i,j,k)
              End If  ! L_pc2
                ext_QCF(i,j,k) = qcf_work(i,j,k)
            End Do ! Rhc_rowlen_do1
          End Do ! Rhc_rows_do1
        End Do ! Wet_mlev_do1
!
! Rhc_rows_do2:
        Do j = 1, rhc_rows
! Rhc_rowlen_do2:
          Do i = 1, rhc_row_length
            ext_p_layer_centres(i,j,0) = p_layer_centres(i,j,0)
            ext_LAND(i,j) = land_sea_mask(i,j)
            ext_ICE_FRAC(i,j) = ice_fract(i,j)
            ext_land_fract(i,j) = 0.0
          End Do ! Rhc_rowlen_do2
        End Do ! Rhc_rows_do2
        Do k = 1, land_points
          j = (land_index(k)-1)/row_length + 1
          i = land_index(k) - (j-1)*row_length
          ext_land_fract(i,j) = fland(k)
        End Do

!
! Synchronize haloes.
!
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &     ext_p_layer_centres, rhc_row_length, rhc_rows,               &
     &         wet_model_levels+1, 1, 1, fld_type_p,  .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_p_layer_centres,rhc_row_length,    &
     &                     rhc_rows,                                    &
     &                     wet_model_levels+1,1,1)
!
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   ext_TL, rhc_row_length, rhc_rows,              &
     &                   wet_model_levels, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_TL,rhc_row_length,                 &
     &                     rhc_rows,                                    &
     &                     wet_model_levels,1,1)
!
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   ext_QL, rhc_row_length, rhc_rows,              &
     &                   wet_model_levels, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_QL,rhc_row_length,                 &
     &                     rhc_rows,                                    &
     &                     wet_model_levels,1,1)
!
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   ext_QCF, rhc_row_length, rhc_rows,             &
     &                   wet_model_levels, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_QCF,rhc_row_length,                &
     &                     rhc_rows,                                    &
     &                     wet_model_levels,1,1)
!
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   ext_LAND, rhc_row_length, rhc_rows,            &
     &                   1, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_land,rhc_row_length,               &
     &                     rhc_rows,                                    &
     &                     1,1,1)

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   ext_land_fract, rhc_row_length, rhc_rows,      &
     &                   1, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_land_fract, rhc_row_length,        &
     &                     rhc_rows,                                    &
     &                     1,1,1)
!
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   ext_ICE_FRAC, rhc_row_length, rhc_rows,        &
     &                   1, 1, 1, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        CALL FILL_EXTERNAL_HALOS(ext_ICE_FRAC,rhc_row_length,           &
     &                     rhc_rows,                                    &
     &                     1,1,1)


!
! DEPENDS ON: ls_calc_rhcrit
        Call ls_calc_rhcrit( ext_p_layer_centres                        &
!            Array dimensions
     &   , wet_model_levels, rhc_row_length, rhc_rows, bl_levels        &
     &   , global_row_length                                            &
!            Prognostic Fields
     &   , ext_TL, ext_QL, ext_QCF,ext_LAND,ext_land_fract,ext_ICE_FRAC &
!            Logical controls
     &   , l_mixing_ratio                                               &
!            Output
     &   , RHCPT)
!
      Else
!       RHCRIT is 1D Parameter read in from namelist
        Do k = 1, wet_model_levels
          RHCPT(1,1,k) = rhcrit(k)
        End Do
      End if  ! Lrhcpt_if1
!
      If (.not. L_pc2) Then

! ----------------------------------------------------------------------
! Section CLD.1.b Calculate (2A) large-scale cloud fraction/ condensate.
! ----------------------------------------------------------------------
!
! DEPENDS ON: ls_cld
      Call ls_cld( p_layer_centres(1,1,1), RHCPT,                       &
     &             wet_model_levels, bl_levels, row_length, rows,       &
     &             rhc_row_length, rhc_rows,                            &
     &             cloud_fraction_method,overlap_ice_liquid,            &
     &             ice_fraction_method,ctt_weight,t_weight,             &
     &             qsat_fixed,sub_cld,                                  &
     &             ntml, cumulus, L_eacf, l_mixing_ratio, T_work,       &
     &             cf_work, q_work, qcf_work, qcl_work,                 &
     &             cfl_work, cff_work, error_code )

      End If  ! .not. L_pc2
!
! Call timer for cloud code
! DEPENDS ON: timer
      If (Ltimer) Call timer ('LS Cloud',4)
!
      If (error_code  ==  0 ) Then

! ----------------------------------------------------------------------
! Section LSP.1 Call Large Scale precipitation scheme.
! ----------------------------------------------------------------------
! Call timer for LS rain code
! DEPENDS ON: timer
        If (Ltimer) Call timer ('LS Rain ',3)

! initialise large scale rain and snow to zero.
        Do j = 1, rows
          Do i = 1, row_length
            ls_rain(i,j) = 0.
            ls_snow(i,j) = 0.
          End Do
        End Do

! Set global land fraction
        Do j = 1, rows
          Do i = 1, row_length
            land_fract(i,j) = 0.0
          End Do
        End Do
        Do k = 1, land_points
          j = (land_index(k)-1)/row_length + 1
          i = land_index(k) - (j-1)*row_length
          land_fract(i,j) = fland(k)
        End Do

! Calculate sea-salt concentration if required

        If (L_seasalt_CCN) then

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                height(i,j,k) = r_theta_levels(i,j,k)                   &
     &                                          -r_theta_levels(i,j,0)
              End Do
            End Do
          End Do

! DEPENDS ON: set_seasalt
          Call set_seasalt(u_1, v_1, height, land_fract, ice_fract      &
     &                    , row_length, rows, model_levels              &
     &                    , bl_levels, sea_salt_film, sea_salt_jet)

        Else
          sea_salt_film(1, 1, 1)=0.0
          sea_salt_jet(1, 1, 1)=0.0
        End if

      If (L_mcr_qcf2) Then  ! Second cloud ice variable in use
        ! Store ice variables in qcf_work and qcf2_work for LS_PPN
        qcf_work(:,:,:)  = qcf_n(:,:,:)
        qcf2_work(:,:,:) = qcf2_n(:,:,:)
      End If

! theta_star holds current estimate of temperature at new time level.
! 3A
! DEPENDS ON: ls_ppn
        Call ls_ppn(                                                    &
     &              halo_i, halo_j, off_x, off_y,                       &
     &              p_layer_boundaries, p_layer_centres(1,1,1),         &
     &              timestep, land_sea_mask, cw_sea, cw_land,           &
     &              L_seq_mcr,L_autoc_3b,L_autolim_3b,                  &
     &              L_autoconv_murk,ec_auto,                            &
     &              N_drop_land,N_drop_sea,N_drop_land_cr,              &
     &              N_drop_sea_cr,Ntot_land, Ntot_sea,                  &
     &              x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic,           &
     &              lsp_ei,lsp_fi,lsp_eic,lsp_fic,                      & 
     &              cf_work, cfl_work, cff_work,                        &
     &              RHCPT, wet_model_levels, model_levels, bl_levels,   &
     &             lspice_dim1,lspice_dim2,lspice_dim3,                 &
     &              rho_r2, r_rho_levels, r_theta_levels,               &
     &              q_work, qcf_work, qcl_work, T_work,                 &
     &              qcf2_work, qrain_work, qgraup_work,                 &
     &              L_mcr_qcf2,L_mcr_qrain,L_mcr_qgraup,L_it_melting,   &
     &              L_use_sulphate_autoconv, L_auto_debias,             &
     &              l_mixing_ratio, l_cry_agg_dep, l_droplet_settle,    &
     &              sea_salt_film, sea_salt_jet,                        &
     &              L_seasalt_CCN, salt_dim1, salt_dim2, salt_dim3,     &
     &              L_use_biogenic, biogenic,                           &
     &              snow_depth, land_fract,                             &
     &              so2( 1:row_length, 1:rows, 1:wet_model_levels),     &
     &              l_sulpc_so2,                                        &
     &              nh3( 1:row_length, 1:rows, 1:wet_model_levels ),    &
     &              l_sulpc_nh3,                                        &
     &              so4_aitken(1:row_length, 1:rows,1:wet_model_levels),&
     &              so4_accu( 1:row_length, 1:rows, 1:wet_model_levels),&
     &              so4_diss( 1:row_length, 1:rows, 1:wet_model_levels),&
     &              aged_soot( 1:row_length, 1:rows,1:wet_model_levels),&
     &              l_soot,                                             &
     &              aged_bmass(1:row_length,1:rows,1:wet_model_levels), &
     &              cloud_bmass(1:row_length,1:rows,1:wet_model_levels),&
     &              L_bmass_CCN,                                        &
     &              aged_ocff(1:row_length,1:rows,1:wet_model_levels),  &
     &              cloud_ocff(1:row_length,1:rows,1:wet_model_levels), &
     &              L_ocff_CCN,                                         &
     &              aerosol( 1:row_length, 1:rows, 1:wet_model_levels), &
     &              l_murk, l_pc2,                                      &
     &              ls_rain, ls_snow,                                   &
     &              lscav_so2, lscav_nh3,                               &
     &              lscav_so4ait, lscav_so4acc,                         &
     &              lscav_so4dis,                                       &
     &              LS_RAIN3D,LS_SNOW3D,RAINFRAC3D,                     &
     &              row_length, rows,                                   &
     &              rhc_row_length, rhc_rows,                           &
      ! Process rate diagnostics
     &  PSDEP,PSAUT,PSACW,PSACR,PSACI,PSMLT,PSMLTEVP                    &
     &, PRAUT,PRACW,PREVP                                               &
     &, PGAUT,PGACW,PGACS,PGMLT                                         &
     &, PIFRW,PIPRM,PIDEP,PIACW,PIACR,PIMLT,PIMLTEVP                    &
     &, PIFALL,PSFALL,PRFALL,PGFALL,PLSET,PLEVPSET                      &
      ! Process rate diagnostic switches
     &, L_PSDEP_diag,L_PSAUT_diag,L_PSACW_diag,L_PSACR_diag             &
     &, L_PSACI_diag,L_PSMLT_diag,L_PSMLTEVP_diag                       &
     &, L_PRAUT_diag,L_PRACW_diag,L_PREVP_diag                          &
     &, L_PGAUT_diag,L_PGACW_diag,L_PGACS_diag,L_PGMLT_diag             &
     &, L_PIFRW_diag,L_PIPRM_diag,L_PIDEP_diag,L_PIACW_diag             &
     &, L_PIACR_diag,L_PIMLT_diag,L_PIMLTEVP_diag                       &
     &, L_PIFALL_diag,L_PSFALL_diag,L_PRFALL_diag,L_PGFALL_diag         &
     &, L_PLSET_diag,L_PLEVPSET_diag                                    &
      ! Variables for stochastic physics random parameters
     &, M_CI,                                                           &
     &              maxsects, h_sect, Error_code )

! Calculate increment fields
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              T_inc(i,j,k) = T_work(i,j,k) - T_n(i,j,k) + T_inc(i,j,k)
              q_inc(i,j,k) = q_work(i,j,k) - q_n(i,j,k) + q_inc(i,j,k)
              qcl_inc(i,j,k) = qcl_work(i,j,k) - qcl_n(i,j,k)           &
     &                       + qcl_inc(i,j,k)
              qcf_inc(i,j,k) = qcf_work(i,j,k) - qcf_n(i,j,k)           &
     &                       + qcf_inc(i,j,k)
              cf_inc(i,j,k) = cf_work(i,j,k)   - cf_n(i,j,k)            &
     &                       + cf_inc(i,j,k)
              cfl_inc(i,j,k) = cfl_work(i,j,k)   - cfl_n(i,j,k)         &
     &                       + cfl_inc(i,j,k)
              cff_inc(i,j,k) = cff_work(i,j,k)   - cff_n(i,j,k)         &
     &                       + cff_inc(i,j,k)
            End Do
          End Do
        End Do
        Do k = 1, bl_levels
          Do j = 1, rows
            Do i = 1, row_length
!             ! Calculate and save microphys increments to TL and QW,
!             ! required in boundary layer.
              micro_tends(i,j,k,1) = T_inc(i,j,k)                       &
     &                              - (      lc*qcl_inc(i,j,k) ) / Cp   &
     &                              - ( (lc+lf)*qcf_inc(i,j,k) ) / Cp
              micro_tends(i,j,k,2) = q_inc(i,j,k) + qcl_inc(i,j,k)      &
     &                                            + qcf_inc(i,j,k)
              micro_tends(i,j,k,1) = micro_tends(i,j,k,1)/timestep
              micro_tends(i,j,k,2) = micro_tends(i,j,k,2)/timestep
            enddo
          enddo
        enddo

      If (L_mcr_qcf2)                                                   &
                       ! Second cloud ice variable in use
     &  qcf2_inc(:,:,:) = qcf2_work(:,:,:) - qcf2_n(:,:,:)              &
     &                  + qcf2_inc(:,:,:)

      If (L_mcr_qrain)                                                  &
                        ! Prognostic rain in use
     &  qrain_inc(:,:,:) = qrain_work(:,:,:) - qrain_n(:,:,:)           &
     &                   + qrain_inc(:,:,:)

      If (L_mcr_qgraup)                                                 &
                         ! Prognostic graupel in use
     &  qgraup_inc(:,:,:) = qgraup_work(:,:,:) - qgraup_n(:,:,:)        &
     &                    + qgraup_inc(:,:,:)

! Call timer for LS rain code
! DEPENDS ON: timer
        If (Ltimer) Call timer ('LS Rain ',4)
      End If ! on error code zero

! ----------------------------------------------------------------------
! Section LSP.1.1 Tracer scavenging
! ----------------------------------------------------------------------
!
!     Call timer for LS scavenging code
! DEPENDS ON: timer
        If (Ltimer) Call timer ('LS Scavenging',3)

! below cloud scavenging of dust
!
! Scavenging of mineral dust
! (at present this is simplistic scheme based on 4.5 code)
!
      IF (L_DUST) THEN



! DEPENDS ON: mass_calc
        CALL MASS_CALC(                                                 &
     & ROW_LENGTH, ROWS, MODEL_LEVELS, WET_MODEL_LEVELS,                &
     & R_RHO_LEVELS(1:row_length,1:rows,1:model_levels),                &
     & R_THETA_LEVELS(1:row_length,1:rows,0:model_levels),              &
     & TIMESTEP, RHO_R2(1:row_length,1:rows,1:model_levels),            &
     & Q_N, QCL_N, QCF_N,                                               &
     & DM )


!       put dust arrays together for simplicity
        DO IDIV = 1,NDIV
          DO J = 1,ROWS
            DO I = 1,ROW_LENGTH
              LSCAV_DUST_ALL(I,J,IDIV) = 0.
            ENDDO
          ENDDO
        ENDDO

        DO K = 1,MODEL_LEVELS
          DO J = 1,ROWS
            DO I = 1,ROW_LENGTH
              DUST_ALL(I,J,K,1)=DUST_DIV1(I,J,K)
              DUST_ALL(I,J,K,2)=DUST_DIV2(I,J,K)
              DUST_ALL(I,J,K,3)=DUST_DIV3(I,J,K)
              DUST_ALL(I,J,K,4)=DUST_DIV4(I,J,K)
              DUST_ALL(I,J,K,5)=DUST_DIV5(I,J,K)
              DUST_ALL(I,J,K,6)=DUST_DIV6(I,J,K)
             ENDDO
           ENDDO
         ENDDO


         DO K = WET_MODEL_LEVELS,1,-1
           DO J = 1, ROWS
             DO I = 1,ROW_LENGTH
!            deal with possible -ive precip
               RAINRATE = MAX(LS_RAIN3D(I,J,K),0.)
               SNOWRATE = MAX(LS_SNOW3D(I,J,K),0.)
!              calc proportion of dust mixing ratio scavenged
!CDIR UNROLL=NDIV
               DO IDIV = 1, NDIV
                 RATE = ( RAINRATE * KRAIN_DUST(IDIV) +                 &
     &                 SNOWRATE * KSNOW_DUST(IDIV) ) *                  &
     &                 3600.0 * TIMESTEP
                 DELTA_DUST=DUST_ALL(I,J,K,IDIV)*(1.0-1.0/(1.0+ RATE ))
!                calc mass of dust removed
                 LSCAV_DUST_ALL(I,J,IDIV)=LSCAV_DUST_ALL(I,J,IDIV) +    &
     &           DELTA_DUST * DM (I,J,K)
!                decrement mixing ratio
                 DUST_ALL(I,J,K,IDIV)=DUST_ALL(I,J,K,IDIV) - DELTA_DUST
               ENDDO !IDIV
             ENDDO !ROW_LENGTH
           ENDDO !ROWS
         ENDDO ! level 1

!       put newly calculated values into main arrays

        DO K = 1,MODEL_LEVELS
          DO J = 1,ROWS
            DO I = 1,ROW_LENGTH
              DUST_DIV1(I,J,K)=DUST_ALL(I,J,K,1)
              DUST_DIV2(I,J,K)=DUST_ALL(I,J,K,2)
              DUST_DIV3(I,J,K)=DUST_ALL(I,J,K,3)
              DUST_DIV4(I,J,K)=DUST_ALL(I,J,K,4)
              DUST_DIV5(I,J,K)=DUST_ALL(I,J,K,5)
              DUST_DIV6(I,J,K)=DUST_ALL(I,J,K,6)
            ENDDO
          ENDDO
        ENDDO

      ENDIF !L_DUST
!
!
!
!
      If (l_sulpc_so2) Then
!
! DEPENDS ON: rainout_intctl
        Call rainout_intctl(row_length, rows,                           &
     &               off_x, off_y,                                      &
     &               halo_i, halo_j,                                    &
     &               model_levels, wet_model_levels,                    &
     &               r_rho_levels, r_theta_levels,                      &
     &               rho_r2, q_n,                                       &
     &               qcf_work, qcl_work,                                &
     &               qcf_n, qcl_n,                                      &
     &               ls_rain3d, ls_snow3d,                              &
     &               timestep,                                          &
     &               so4_diss, so4_accu,                                &
     &               RNOUT_TRACER)
!
! DEPENDS ON: sl3dwash
        Call sl3dwash(row_length, rows,                                 &
     &                off_x, off_y,                                     &
     &                halo_i, halo_j,                                   &
     &                model_levels, wet_model_levels,                   &
     &                r_rho_levels, r_theta_levels,                     &
     &                timestep,                                         &
     &                rho_r2,                                           &
     &                q_n, qcl_n, qcf_n, so2,                           &
     &                ls_rain3d,                                        &
     &                LSCAV_TR)
!
          If (l_sulpc_nh3) Then
!
! DEPENDS ON: nh3dwash
          Call nh3dwash(row_length, rows,                               &
     &                  off_x, off_y,                                   &
     &                  halo_i, halo_j,                                 &
     &                  model_levels, wet_model_levels,                 &
     &                  r_rho_levels, r_theta_levels,                   &
     &                  timestep,                                       &
     &                  rho_r2,                                         &
     &                  q_n, qcl_n, qcf_n, nh3,                         &
     &                  ls_rain3d,                                      &
     &                  LSCAV_NH3)
!
          End If ! On test for l_sulpc_nh3
!
      End If ! On test for sulphur cycle
!
      If (l_soot) then  ! If soot modelling is included
!
! DEPENDS ON: rainout_intctl
        Call rainout_intctl(row_length, rows,                           &
     &               off_x, off_y,                                      &
     &               halo_i, halo_j,                                    &
     &               model_levels, wet_model_levels,                    &
     &               r_rho_levels, r_theta_levels,                      &
     &               rho_r2, q_n,                                       &
     &               qcf_work, qcl_work,                                &
     &               qcf_n, qcl_n,                                      &
     &               ls_rain3d, ls_snow3d,                              &
     &               timestep,                                          &
     &               cloud_soot, aged_soot,                             &
     &               rnout_soot)

!
!       LS washout of soot was neglected at 4.5, although code was
!       allegedly written to perform the calculation. We can't use
!       the same treatment for soot as we use for the S cycle.
!       For now, we will again neglect this process, but note that
!       we may wish to include it in the future.
!
        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_soot(i,j)=0.0
          End Do
        End Do
!
      End If  ! l_soot
!
      If (l_biomass) then  ! If biomass modelling is included
!
! DEPENDS ON: rainout_intctl
        Call rainout_intctl(row_length, rows,                           &
     &               off_x, off_y,                                      &
     &               halo_i, halo_j,                                    &
     &               model_levels, wet_model_levels,                    &
     &               r_rho_levels, r_theta_levels,                      &
     &               rho_r2, q_n,                                       &
     &               qcf_work, qcl_work,                                &
     &               qcf_n, qcl_n,                                      &
     &               ls_rain3d, ls_snow3d,                              &
     &               timestep,                                          &
     &               cloud_bmass, aged_bmass,                           &
     &               rnout_bmass)

!
!       As with soot, LS washout of biomass smoke aerosol is
!       currently neglected. Again, though, we may wish to include
!       it in the future.
!
        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_bmass(i,j)=0.0
          End Do
        End Do
!
      End If  ! l_biomass
!
      If (l_ocff) then  ! If fossil-fuel org carb modelling is included
!
! DEPENDS ON: rainout_intctl
        Call rainout_intctl(row_length, rows,                           &
     &               off_x, off_y,                                      &
     &               halo_i, halo_j,                                    &
     &               model_levels, wet_model_levels,                    &
     &               r_rho_levels, r_theta_levels,                      &
     &               rho_r2, q_n,                                       &
     &               qcf_work, qcl_work,                                &
     &               qcf_n, qcl_n,                                      &
     &               ls_rain3d, ls_snow3d,                              &
     &               timestep,                                          &
     &               cloud_ocff, aged_ocff,                             &
     &               rnout_ocff)

!
!       As with soot and biomass, LS washout of fossil-fuel organic
!       carbon aerosol is currently neglected. Again, though, we may 
!       wish to include it in the future.
!
        Do i=1,lspice_dim1
          Do j=1,lspice_dim2
            lscav_ocff(i,j)=0.0
          End Do
        End Do
!
      End If  ! l_ocff
!
!     Call timer for LS scavenging code
! DEPENDS ON: timer
        If (Ltimer) Call timer ('LS Scavenging',4)
!

! ----------------------------------------------------------------------
! Section LSP.2 Output Diagnostics
! ----------------------------------------------------------------------

! Check that microphysics diagnostics requested this timestep
      If (error_code  ==  0 .and. sf(0,4)) Then
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',3)
! DEPENDS ON: diagnostics_lsrain
        Call diagnostics_lsrain(                                        &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels                           &
     &,                     lspice_dim1,lspice_dim2,lspice_dim3         &
     &,                      timestep                                   &
     &,                      at_extremity                               &
     &,                      L_DUST                                     &
     &,                      p_layer_centres                            &
     &,                      T_work, q_work, qcl_work, qcf_work         &
     &,                      cf_work, cfl_work, cff_work                &
     &,                      T_n, q_n, qcl_n, qcf_n                     &
     &,                      cf_n, cfl_n, cff_n                         &
     &,                      ls_rain, ls_snow                           &
     &,                      ls_rain3d,ls_snow3d,rainfrac3d             &
     &,                      RNOUT_TRACER,LSCAV_DUST_ALL,LSCAV_TR       &
     &,                      LSCAV_NH3                                  &
     &,                      rnout_soot, lscav_soot                     &
     &,                      rnout_bmass, lscav_bmass                   &
     &,                      rnout_ocff, lscav_ocff                     &
     &,                      l_mixing_ratio                             &
     &,                      PSDEP,PSAUT,PSACW,PSACR                    &
     &,                      PSACI,PSMLT,PSMLTEVP                       &
     &,                      PRAUT,PRACW,PREVP                          &
     &,                      PGAUT,PGACW,PGACS,PGMLT                    &
     &,                      PIFRW,PIPRM,PIDEP,PIACW                    &
     &,                      PIACR,PIMLT,PIMLTEVP                       &
     &,                      PIFALL,PSFALL,PRFALL,PGFALL                &
     &,                      PLSET,PLEVPSET                             &
     &       ,                                                          &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork4                                                       &
     &       )

! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',4)
      End If   ! on error_code and sf(0,4)

      ! Deallocate additional microphysics work arrays
      Deallocate ( qcf2_work )
      Deallocate ( qrain_work )
      Deallocate ( qgraup_work )
      ! Deallocate additional microphysics work arrays
      Deallocate ( PSDEP,PSAUT,PSACW,PSACR,PSACI,PSMLT,PSMLTEVP         &
     &, PRAUT,PRACW,PREVP                                               &
     &, PGAUT,PGACW,PGACS,PGMLT                                         &
     &, PIFRW,PIPRM,PIDEP,PIACW,PIACR,PIMLT,PIMLTEVP                    &
     &, PIFALL,PSFALL,PRFALL,PGFALL,PLSET,PLEVPSET )

! end of routine microphys_ctl
      Return
      END SUBROUTINE microphys_ctl
