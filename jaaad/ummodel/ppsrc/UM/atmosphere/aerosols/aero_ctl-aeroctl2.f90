
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE AERO_CTL(                                              &
!
! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, proc_col_group, at_extremity, n_proc, n_procx   &
     &, n_procy, neighbour, g_rows, g_row_length, g_datastart, me       &
! model dimensions
     &, row_length, rows, n_rows, land_points                           &
     &, model_levels, wet_model_levels, bl_levels, n_cca_levels         &
     &, theta_field_size                                                &
     &, salt_dim1, salt_dim2, salt_dim3                                 &
     &, aero_dim1, aero_dim2, aero_dim3                                 &
! Model switches
     &, model_domain, L_CAL360, L_SEC_VAR, L_EqT, Ltimer                &
! Model parameters
     &, Ntot_land, Ntot_sea                                             &
! Physical constants
     &, lc, lf, cp, two_Omega, p_zero, kappa                            &
     &, R, g, Lapse_Rate, earth_radius, Pi                              &
! Co-ordinate information
     &, r_rho_levels, r_theta_levels                                    &
     &, eta_theta_levels, eta_rho_levels                                &
     &, delta_lambda, delta_phi                                         &
     &, lat_rot_NP, long_rot_NP                                         &
! Time stepping information
     &, timestep                                                        &
     &, val_year, val_day_number, val_hour, val_minute                  &
     &, val_second, timestep_number                                     &
     &, PREVIOUS_TIME                                                   &
     &, CALL_CHEM_FREQ                                                  &
! Trig arrays
     &, sin_theta_longitude, cos_theta_longitude                        &
     &, FV_cos_theta_latitude                                           &
! Grid-dependent arrays
     &, f3_at_u, true_longitude, true_latitude                          &
!
! Data fields IN
     &, u, v, Tstar, Tstar_sea                                          &
     &, theta, q, qcl, qcf                                              &
     &, rho, land_mask, fland_ctile                                     &
     &, p_theta_levels, exner_rho_levels, exner_theta_levels            &
     &, ice_fract, snow_depth                                           &
     &, CLOUDF_halos                                                    &
     &, OH_CONC, H2O2_LMT, HO2_CONC, O3                                 &
     &, SO2_SURFEM, SO2_HILEM, SO2_NATEM                                &
     &, DMS_em_ancil, DMS_conc, NH3_EM                                  &
     &, DMS_Ointer                                                      &
     &, SOOT_SUREM, SOOT_HILEM, BMASS_SUREM, BMASS_HILEM, OCFF_SUREM    &
     &, OCFF_HILEM, SO2_HIGH_LEVEL, SOOT_HIGH_LEVEL, BMASS_HIGH_LEVEL_1 &
     &, BMASS_HIGH_LEVEL_2, OCFF_HIGH_LEVEL, land_index                 &
! Logicals IN
     &, L_SULPC_SO2, L_SULPC_DMS, L_SULPC_OZONE                         &
     &, L_SULPC_SO2_O3_NONBUFFERED, L_SULPC_NH3                         &
     &, L_sulphate_CCN, L_seasalt_CCN, L_SOOT, L_soot_CCN               &
     &, L_biomass, L_biomass_CCN, L_ocff, L_ocff_CCN                    &
     &, L_SO2_SURFEM, L_SO2_HILEM, L_SO2_NATEM, L_DMS_EM                &
     &, L_DMS_em_inter, L_DMS_Liss_Merlivat                             &
     &, L_DMS_Wanninkhof, L_DMS_Nightingale                             &
     &, L_DMS_Ointer                                                    &
     &, L_NH3_EM, L_ctile                                               &
     &, L_SOOT_SUREM, L_SOOT_HILEM, L_BMASS_SUREM, L_BMASS_HILEM        &
     &, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_BIOGENIC                      &
     &, L_USE_SEASALT_DIRECT, L_USE_SEASALT_INDIRECT                    &
     &, L_USE_SEASALT_AUTOCONV, L_DUST                                  &
!
! Data fields IN/OUT
     &, SO2, DMS                                                        &
     &, SO4_AIT, SO4_ACC, SO4_DIS                                       &
     &, H2O2_MXR, NH3                                                   &
     &, SOOT_NEW, SOOT_AGD, SOOT_CLD                                    &
     &, BMASS_NEW, BMASS_AGD, BMASS_CLD                                 &
     &, OCFF_NEW, OCFF_AGD, OCFF_CLD                                    &
     &, BIOGENIC                                                        &
!
! Data fields IN
     &, DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5           &
!
! Data fields OUT
! Diagnostic info
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     &  STASHwork17                                                     &
     &, Error_code                                                      &
     & )
!
!---------------------------------------------------------------------
! Purpose: Interface for Aerosol Modelling, to include Sulphur Cycle
!          and SOOT modelling
!
! Level 2 control routine
!
! Current owners of code:                C E Johnson
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
! System components covered:
!
! System task:
!
! Documentation: Not yet available
!
!-----------------------------------------------------------------
!
      IMPLICIT NONE
!
! Arguments with intent IN:
!
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
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same col
     &, global_rows                                                     &
                           ! NUMBER OF global rows
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbour
                             ! in the horizontal plane
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)                                      &
     &, me         ! My processor number
!
      Logical                                                           &
     &  at_extremity(4) ! Indicates if this processor is at north, sout
                        ! east or west of the processor grid
!
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

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

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end
!
! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, land_points                                                     &
     &, model_levels                                                    &
     &, wet_model_levels                                                &
     &, bl_levels                                                       &
     &, n_cca_levels                                                    &
                        ! No. conv cloud levels (1 if 2D, nlevs if 3D)
     &, theta_field_size                                                &
     &,salt_dim1                                                        &
                                           !dimensns of seasalt array
     &,salt_dim2                                                        &
     &,salt_dim3                                                        &
     &,aero_dim1                                                        &
                                           !row_length or 1
     &,aero_dim2                                                        &
                                           !rows or 1
     &,aero_dim3                           !model_levels or 1
! Model switches
      Integer                                                           &
     &  model_domain
      Logical                                                           &
     &  L_CAL360                                                        &
                        ! T if using 360 day calendar
     &, L_SEC_VAR                                                       & 
                        ! if T include secular varn of earth's orbit,
     &, L_EqT                                                           &
                        ! if T include equn of time         in SOLPOS
     &, Ltimer                                                           
                        ! if T then output some timing information

! Physical constants
      Real                                                              &
     &  lc, lf, cp                                                      &
     &, two_Omega                                                       & 
                        ! twice Earth's rotation rate
     &, p_zero                                                          &
     &, kappa                                                           &
     &, R, g, Lapse_Rate, earth_radius, Pi
! Co-ordinate arrays
      Real                                                              &
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, delta_lambda                                                    &
     &, delta_phi                                                       &
     &, lat_rot_NP                                                      &
     &, long_rot_NP

! Trig arrays
      Real                                                              &
     &  cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)
! Grid-dependent arrays
      Real                                                              &
     &  f3_at_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y)          &
     &, true_longitude(row_length, rows)                                &
     &, true_latitude(row_length, rows)                                 &

     &, FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)
! Time stepping information
      Real                                                              &
     &  timestep                       !atmosphere model timetsep
      Integer                                                           &
     &  val_year                                                        &
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second                                                      &
     &, timestep_number                                                 &
     &, PREVIOUS_TIME(7)                                                &
     &, CALL_CHEM_FREQ                  !frequency of calling chemistry
!                                         per atmos phys timestep
!
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
      Real                                                              &
     &  STASHwork17(*)  ! STASH workspace for section 17 (Aero_Ctl)
!
      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &                                                 model_levels)    &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &                                                 model_levels)    &
     &, theta(1-off_x:row_length+off_x,1-off_y:rows+off_y,              &
     &                                                 model_levels)    &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &                                             wet_model_levels)    &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &                                             wet_model_levels)    &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &                                             wet_model_levels)    &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &                                                 model_levels)    &
     &, p_theta_levels(1-off_x:row_length+off_x,                        &
     &                 1-off_y:rows+off_y, model_levels)                &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                           1-off_y:rows+off_y, model_levels+1)    &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                             1-off_y:rows+off_y, model_levels)    &
     &, ice_fract(row_length, rows)                                     &
     &, snow_depth(row_length, rows)                                    &
     &, land_fract(row_length, rows)                                    &
     &, Tstar(row_length, rows)                                         &
     &, Tstar_sea(row_length, rows)                                     &
     &, fland_ctile(land_points)                                        &
     &, CLOUDF_halos(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,  &
     &                                        wet_model_levels)         &
!
     &,OH_CONC(row_length,rows,model_levels)                            &
     &,HO2_CONC(row_length,rows,model_levels)                           &
     &,H2O2_LMT(row_length,rows,model_levels)                           &
     &,O3(row_length,rows,model_levels)                                 &
     &,SO2_SURFEM(row_length,rows)                                      &
                                          !SO2 emiss at surface
     &,SO2_HILEM(row_length,rows)                                       &
                                          !SO2 emiss at chimney level
     &,SO2_NATEM(row_length,rows,model_levels)                          &
                                                !Volcanic SO2 emiss
     &,DMS_em_ancil(row_length,rows)                                    &
                                          !Ancillary DMS emiss (surf)
     &,DMS_conc(row_length,rows)                                        &
                                          !Seawater DMS conc (n mol l-1)
     &,NH3_EM(row_length,rows)                                          &
                                         !NH3 emiss (surface)
     &,SOOT_SUREM(row_length,rows)                                      &
                                        !SOOT emiss at surface
     &,SOOT_HILEM(row_length,rows)                                      &
                                        !SOOT emiss at chimney lev
     &,BMASS_SUREM(row_length,rows)                                     &
                                        !Biomass surface emissions
     &,BMASS_HILEM(row_length,rows)                                     &
                                        !Biomass high level emissions
     &,OCFF_SUREM(row_length,rows)                                      &
                                        !OCFF emiss at surface
     &,OCFF_HILEM(row_length,rows)
                                        !OCFF emiss at chimney lev

!
      Real                                                              &
     & Ntot_land                                                        &
                           ! Number of droplets over land / m-3
     &,Ntot_sea            ! Number of droplets over sea / m-3

      Integer                                                           &
     & SO2_HIGH_LEVEL                                                   &
                                         !Model level for chimney emiss
     &,SOOT_HIGH_LEVEL                                                  &
                                         !Model level for chimney emiss
     &,BMASS_HIGH_LEVEL_1                                               &
                                         !Lowest and highest model
     &,BMASS_HIGH_LEVEL_2                                               &
                                         !levels with biomass emiss
     &,OCFF_HIGH_LEVEL                                                  &
                                         !Model level for chimney emiss
     &,land_index(land_points)
!
      LOGICAL                                                           &
     & L_SULPC_SO2                                                      &
                                         !T if Sulphur Cycle required
     &,L_SULPC_DMS                                                      &
                                         !T if DMS chemistry required
     &,L_SULPC_OZONE                                                    &
                                         !T if O3 oxidn required
     &,L_SULPC_SO2_O3_NONBUFFERED                                       &
                                         !T if SO2+O3 reaction is NOT
                                         !  to be buffered by NH3.
     &,L_SULPC_NH3                                                      &
                                         !T if NH3 buffering required
                                         !(always T if L_SULPC_OZONE
                                         ! is T)
     &,L_sulphate_CCN                                                   &
                                         !T if sulphate used for CCN
     &,L_seasalt_CCN                                                    &
                                         !T if sea-salt used for CCN
     &,L_soot_CCN                                                       &
                                         !T if soot used for CCN
     &,L_biomass_CCN                                                    &
                                         !T if biomass used for CCN
     &,L_ocff_CCN                                                       &
                                         !T if OCFF used for CCN
     &,L_ctile                                                          &
                                         !T if coastal tiling is on
     &,L_SOOT                                                           &
                                         !T if SOOT modelling required
     &,L_BIOMASS                                                        &
                                         !T if biomass modelling reqd
     &,L_OCFF                                                           &
                                         !T if OCFF modelling required
     &,L_SO2_SURFEM                                                     &
                                         !T if surface SO2 ems present
     &,L_SO2_HILEM                                                      &
                                         !T if high lev SO2 em present
     &,L_SO2_NATEM                                                      &
                                         !T if volcanic SO2 ems present
     &,L_DMS_EM                                                         &
                                         !T if DMS emiss present (surf)
     &,L_DMS_em_inter                                                   &
                                         !T if interactive DMS emiss
     &,L_DMS_Ointer                                                     &
                                         !T if ocean DMS emiss model
     &,L_DMS_Liss_Merlivat                                              &
                                         !Switches to determine which
     &,L_DMS_Wanninkhof                                                 &
                                         !  scheme to use for inter-
     &,L_DMS_Nightingale                                                &
                                         !  active DMS emissions
     &,L_NH3_EM                                                         &
                                         !T if NH3 emiss present (surf)
     &,L_SOOT_SUREM                                                     &
                                         !T if surface SOOT ems present
     &,L_SOOT_HILEM                                                     &
                                         !T if high lev SOOT ems presnt
     &,L_BMASS_SUREM                                                    &
                                         !T if sfc biomass ems present
     &,L_BMASS_HILEM                                                    &
                                         !T if hi lev bmass ems present
     &,L_OCFF_SUREM                                                     &
                                         !T if surface OCFF ems present
     &,L_OCFF_HILEM                                                     &
                                         !T if high lev OCFF ems presnt
     &,L_USE_BIOGENIC                                                   &
                                         !T if using biogenics for CCN
     &,LAND_MASK(row_length,rows)                                       &    
                                         !T IF LAND, F IF SEA
     &,L_DUST                                                           &
                                         !T if mineral dust used
     &,L_USE_SEASALT_DIRECT                                             &
                                         !T if SS dir. rad. effect.
     &,L_USE_SEASALT_INDIRECT                                           &
                                         !T if SS 1st indir. effect
     &,L_USE_SEASALT_AUTOCONV
                                         !T if SS 2nd indir. effect
!
! Arguments with intent IN/OUT:
      REAL                                                              &
     & SO2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &                                model_levels)                     &
                                                         !mmr S in SO2
     &,DMS(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &                                model_levels)                     &
                                                         !mmr S in DMS
     &,SO4_AIT(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                model_levels)                     &
                                                         !mmr S in AIT
     &,SO4_ACC(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                model_levels)                     &
                                                         !mmr S in ACC
     &,SO4_DIS(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                model_levels)                     &
                                                         !mmr S in DIS
     &,NH3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &                                model_levels)                     &
                                                         !mmr N in NH3
     &,H2O2_MXR(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                         !mmr H2O2
     &,SOOT_NEW(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr fresh soot
     &,SOOT_AGD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr aged soot
     &,SOOT_CLD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr sootincloud
     &,BMASS_NEW(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                model_levels)                     &
                                                      !mmr fresh smoke
     &,BMASS_AGD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                model_levels)                     &
                                                      !mmr aged smoke
     &,BMASS_CLD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                model_levels)                     &
                                                      !mmr cloud smoke
     &,OCFF_NEW(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr fresh OCFF
     &,OCFF_AGD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr aged OCFF
     &,OCFF_CLD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &                                model_levels)                     &
                                                      !mmr OCFF incloud
     &,BIOGENIC(row_length, rows, model_levels)       !mmr biogenics
!
! Arguments with intent IN:
      REAL, INTENT(IN) ::                                               &
     & DUST_DIV1(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                model_levels)                     &   
                                                      !mmr Dust div 1
     &,DUST_DIV2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &                                model_levels)                     &      
                                                      !mmr Dust div 2
     &,DUST_DIV3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &         
     &                                model_levels)                     &   
                                                      !mmr Dust div 3
     &,DUST_DIV4(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &         
     &                                model_levels)                     &
                                                      !mmr Dust div 4
     &,DUST_DIV5(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &         
     &                                model_levels)                        
                                                      !mmr Dust div 5
! arguments with intent OUT:
      Integer                                                           &
     &  Error_code
!
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL                                                              &
     &     EVAPTAU,                                                     &
                          ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,                                                      &
                          ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,                                                 &
                          ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,                                                  &
                                  ! high pressure reaction rate limit
     &     K_DMS_OH,                                                    &
                                  ! reaction rate for DMS+OH  cc/mcl/s
     &      K4_CH3SO2_O3,                                               &
                             ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
     &      K5_CH3SO3_HO2,                                              &
                             ! Rate coeff for CH3SO3+HO2 -> MSA+O2
     &      RMM_O3,                                                     &
                             ! relative molecular mass O3
     &     BRAT_SO2,                                                    &
                                  ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,                                                    &
                                  ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,                                                    &
                                 ! no. of molecules in 1 mole
     &     RMM_H2O2,                                                    &
                                 ! relative molecular mass H2O2 kg/mole
     &     RMM_AIR,                                                     &
                                 ! relative molecular mass dry air
     &     RMM_W,                                                       &
                                 ! relative molecular mass water
     &     RELM_S_H2O2,                                                 &
                                 ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,                                                   &
                              ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
     &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
     &     T1,                                                          &
                                !
     &     FC,                                                          &
                                ! parameters for interpolation between
     &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
     &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,                                                    &
                                !
     &     CLOUDTAU,                                                    &
                                  ! air parcel lifetime in cloud
     &     CHEMTAU,                                                     &
                                  ! chem lifetime in cloud before oxidn
     &     O3_MIN,                                                      &
                              ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (                                                       &
     &           EVAPTAU = 300.0,                                       &
                                              ! secs  (=5 mins)
     &             NUCTAU = 30.0,                                       &
                                          ! secs
     &       DIFFUSE_AIT = 1.7134E-9,                                   &
                                             ! sq m/s
     &        K_SO2OH_HI = 2.0E-12,                                     &
                                       ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,                                    &
                                          ! cc/mcl/s
     &       K4_CH3SO2_O3 = 1.0E-14,                                    &
                                        ! cc/mcl/s
     &      K5_CH3SO3_HO2 = 4.0E-11,                                    &
     &             RMM_O3 = 4.8E-2,                                     &
                                        ! kg/mole
     &          BRAT_SO2 = 0.9,                                         &
     &           BRAT_MSA = 1.0-BRAT_SO2,                               &
     &           AVOGADRO = 6.022E23,                                   &
                                          ! per mole
     &           RMM_H2O2 = 3.40E-2,                                    &
                                          ! kg/mole
     &            RMM_AIR = 2.896E-2,                                   &
                                          ! kg/mole
     &              RMM_W = 1.8E-2,                                     &
                                          ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,                                 &
     &           RELM_S_2N = 3.206/2.80,                                &
     &               PARH = 3.3,                                        &
     &                K1 = 4.0E-31,                                     &
                                       ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,                                      &
                                          ! K
     &                FC = 0.45,                                        &
                                        ! from STOCHEM model
     &              FAC1 = 1.1904,                                      &
                                    ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.2E-13,                                    &
                                          ! cc/mcl/s
     &                 K3 = 1.9E-33,                                    &
                                          ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,                                    &
                                          ! cc/mcl
     &                 T2 = 600.0,                                      &
                                          ! K
     &                 T3 = 890.0,                                      &
                                          ! K
     &                 T4 = 2200.0,                                     &
                                          ! K
     &           CLOUDTAU = 1.08E4,                                     &
                                          ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,                                      &
                                          ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,                                    &
                                        !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8                                      &
                                          ! kg/kg
     &          )
!
      REAL RAD_AIT,                                                     &
                            ! median radius of Aitken mode particles
     &     DIAM_AIT,                                                    &
                            !   "    diameter    "
     &     RAD_ACC,                                                     &
                            ! median radius of acccumulation mode
     &     DIAM_ACC,                                                    &
                            !   "    diameter    "
     &     CHI,                                                         &
                            ! mole fraction of S in particle
     &     RHO_SO4,                                                     &
                            ! density of  SO4 particle
     &     SIGMA,                                                       &
                            ! standard devn of particle size distn
!                                 for accumulation mode
     &     E_PARM,                                                      &
                            ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (                                                       &
     &           RAD_AIT = 6.5E-9,                                      &
                                             ! m
     &          DIAM_AIT = 2.0*RAD_AIT,                                 &
     &           RAD_ACC = 95.0E-9,                                     &
                                             ! m
     &          DIAM_ACC = 2.0*RAD_ACC,                                 &
     &               CHI = 32.0/132.0,                                  &
     &           RHO_SO4 = 1769.0,                                      &
                                              ! kg/m3
     &             SIGMA = 1.4,                                         &
     &            E_PARM = 0.9398,                                      &
     &          NUM_STAR = 1.0E6                                        &
                                             ! m-3
     &          )
!
      REAL BOLTZMANN       !Boltzmanns constant.
      REAL MFP_REF         !Reference value of mean free path.
      REAL TREF_MFP        !Reference temperature for mean free path.
      REAL PREF_MFP        !Reference pressure for mean free path.
      REAL SIGMA_AIT       !Geometric standard deviation of the Aitken
!                             mode distribution.
!
      PARAMETER (BOLTZMANN = 1.3804E-23)  ! J K-1
      PARAMETER (MFP_REF = 6.6E-8                                       &
                                          ! m
     &        ,  TREF_MFP = 293.15                                      &
                                          ! K
     &        ,  PREF_MFP = 1.01325E5)    ! Pa
      PARAMETER (SIGMA_AIT = 1.30)
!
!*---------------------------------------------------------------------
!
! Variables with intent OUT (diagnostics):
      REAL                                                              &
     & MSA(row_length,rows,model_levels)                                &
                                                  !mmr S in MSA
     &,F_DMS_TO_SO2(row_length,rows,model_levels)                       &
                                                  !frac oxid DMS to SO2
     &,F_DMS_TO_SO4(row_length,rows,model_levels)                       &
                                                  !frac oxid DMS to SO4
     &,F_DMS_TO_MSA(row_length,rows,model_levels)                       &
                                                  !frac oxid DMS to MSA
     &,NH3_DEP(row_length,rows,model_levels)                            &
                                                  !NH3 depleted
     &,DELTAS_DRY(row_length,rows,model_levels)                         &
                                                   !SO2 dry ox per ts
     &,DELTAS_WET(row_length,rows,model_levels)                         &
                                                   !SO2 wet ox by H2O2
     &,DELTAS_WET_O3(row_length,rows,model_levels)                      &
                                                   !SO2 wet ox by O3
     &,DELTAS_TOT(row_length,rows,model_levels)                         &
                                                   !total SO2 ox per ts
     &,DELTAS_DMS(row_length,rows,model_levels)                         &
                                                   !DMS dry ox per ts
     &,DELTAS_EVAP(row_length,rows,model_levels)                        &
                                                   !SO4_DIS released by
!                             evapn of cloud droplets to SO4_ACC per ts
     &,DELTAS_NUCL(row_length,rows,model_levels)                        &
                                                   !SO4_ACC transfd by
!                                          nucleation to SO4_DIS per ts
     &,DELTAS_DIFFUSE(row_length,rows,model_levels)                     &
                                                    !SO4_AIT transfd to
!                                           SO4_DIS by diffusion per ts
     &,DELTAS_MERGE(row_length,rows,model_levels)                       &
     &,DELTAS_COAG(row_length,rows,model_levels)                        &
     &,PSI(row_length,rows,model_levels)
!
! Local variables
!
      Integer i,j,k,n                                                   &
     &, i_start                                                         &
                        ! Row start point for polar row tidying
     &, istat           ! Status (error code) indicator
!
      REAL CHEMSTEP                         ! chemistry timestep
! Diagnostic increments generated by each call of SULPHR
      REAL                                                              &
     & MSA_inc(row_length,rows,model_levels)                            &
     &,NH3_DEP_inc(row_length,rows,model_levels)
!
!
      REAL T(row_length,rows,model_levels)                              &
                                            ! Temp (K) on theta levels
     &,    T_surf(row_length,rows)          ! Surface temperature (K)
!
! For calculation of cos zenith angle
      REAL COSZA2D(row_length,rows)      !cos zenith angle
      Real                                                              &
     &  Sindec                                                          &
                                         ! sin(solar declination)
     &, hour_angle(row_length,rows)                                     &
     &, SCS                                                             &
                                         ! solar constant scaling factor
     &, seconds_since_midnight                                          &
     &, Eq_Time                                                         &
                                         ! The Equation of time
     &, day_fraction(row_length,rows)                                   &
     &, sin_true_latitude(row_length, rows)
!
! For sea-salt aerosol and DMS emissions calculations
      Real                                                              &
     &  u_1(aero_dim1,aero_dim2)                                        &
                                              ! surface wind u
     &, v_1(aero_dim1,aero_dim2)                                        &
                                              ! surface wind v
     &, u_1_mean                                                        &
                                              ! mean of u_1 for N Pole
     &, v_1_mean                                                        &
                                              ! mean of v_1 for N Pole
     &, sea_salt_film(salt_dim1,salt_dim2,salt_dim3)                    &
                                                     ! Sea-salt
     &, sea_salt_jet(salt_dim1,salt_dim2,salt_dim3)                     &
                                                     !    aerosol.
     &, height(aero_dim1,aero_dim2,aero_dim3)                           &
                                              ! Layer centre heights
     &, DMS_em_inter(aero_dim1,aero_dim2)                               &
                                              ! Interactive DMS emissns
     &, DMS_Ointer(aero_dim1,aero_dim2)                                 &
                                            ! Interactive DMS emissns
     &, DMS_Ointer_mean                                                 &
                                            ! N Pole mean DMS_Ointer
     &, DMS_em_combined(row_length,rows)      ! Combined DMS emissns
!
! Increments from soot ageing and diffusional scavenging
      Real                                                              &
     &  delta_agesoot(row_length,rows,model_levels)                     &
                                   ! Increment from soot ageing
     & ,delta_sootdiffscav(row_length,rows,model_levels)
                                   ! Increment from soot diff. scav.
!
! Increments from biomass smoke ageing and nucleation scavenging
      Real                                                              &
     &  delta_agebmass(row_length,rows,model_levels)                    &
                                    ! Increment from smoke ageing
     & ,delta_bmassnuclscav(row_length,rows,model_levels)
                                    ! Increment from smoke nucl. scav.
! For emissions of biomass smoke
      Integer biomass_level         ! Model level of biomass emissions
      Real biomass_per_level(row_length,rows)
                                    ! Quantity of biomass smoke emitted
                                    ! per model level
!
! Increments from fossil-fuel organic carbon ageing and nucl. scavenging
      Real                                                              &
     &  delta_ageocff(row_length, rows, model_levels)                   &
                                    ! Increment from OCFF ageing
     & ,delta_ocffnuclscav(row_length, rows, model_levels)
                                    ! Increment from OCFF nucl. scav.

!
! For cloud fraction without halos
      Real                                                              &
     &  cloudf(row_length,rows,wet_model_levels)
!
      LOGICAL L_SULPC_NEWDMS        ! TRUE if new DMS scheme required
!
      PARAMETER (L_SULPC_NEWDMS = .TRUE.)
!
! For droplet number calculation
      Real n_droplet(row_length, rows, wet_model_levels)                &
     &,    rho_air                                                      &
     &,    number_droplet                                               &
     &,    ss_film, ss_jet                                              &
     &,    sulp_ait, sulp_acc, sulp_dis                                 &
     &,    agd_bmass, cld_bmass                                         &
     &,    agd_ocff, cld_ocff                                           &
     &,    biogenic_ccn
!
! For PM10 & PM2.5 calculation
      Real                                                              &
     &  PM10(row_length, rows, model_levels)                            &
     &, PM2p5(row_length, rows, model_levels)                           &
     &, PM10_SO4 (row_length, rows, model_levels)                       &
     &, PM2p5_SO4(row_length, rows, model_levels)                       &
     &, PM10_BC (row_length, rows, model_levels)                        &
     &, PM2p5_BC(row_length, rows, model_levels)                        &
     &, PM10_BB (row_length, rows, model_levels)                        &
     &, PM2p5_BB(row_length, rows, model_levels)                        &
     &, PM10_OCFF (row_length, rows, model_levels)                      &
     &, PM2p5_OCFF(row_length, rows, model_levels)                      &
     &, PM10_SOA (row_length, rows, model_levels)                       &
     &, PM2p5_SOA(row_length, rows, model_levels)                       &
     &, PM10_SS (row_length, rows, model_levels)                        &
     &, PM2p5_SS(row_length, rows, model_levels)                        &
     &, PM10_DUST (row_length, rows, model_levels)                      &
     &, PM2p5_DUST(row_length, rows, model_levels)
!
! External functions and subroutines
      External TRSRCE, SOLPOS, SOLANG, SET_SEASALT, SULPHR              &
     &, AGESOOT, SOOTDIFFSCAV                                           &
     &, DMS_FLUX                                                        &
     &, AGEBMASS, BMASSNUCLSCAV                                         &
     &, AGEOCFF, OCFFNUCLSCAV                                           &
     &, gcg_rvecsumr, gcg_rvecsumf                                      &
     &, diagnostics_aero                                                &
     &, number_droplet                                                  & 
     &, pm10_pm2p5
!
! Set up global land fraction field
!
      Do j = 1, rows
        Do i = 1, row_length
          land_fract(i,j) = 0.0
        End Do
      End Do

      If (L_ctile) Then
        Do k = 1,land_points
          j = (land_index(k)-1)/row_length + 1
          i = land_index(k) - (j-1)*row_length
          land_fract(i,j) = fland_ctile(k)
        End Do
      Else
        Do j = 1, rows
          Do i = 1, row_length
            If (LAND_MASK(i,j)) Then
              land_fract(i,j) = 1.0
            Else
              land_fract(i,j) = 0.0
            End If
          End Do
        End Do
      End If

!
! Copy cloud fraction into an array without halos and set to zero if
! smaller than machine precision
          Do k = 1,wet_model_levels
            Do j = 1,rows
              Do i = 1,row_length
                cloudf(i,j,k)=cloudf_halos(i,j,k)
                if (cloudf(i,j,k) <  epsilon(cloudf(i,j,k)))            &
     &            cloudf(i,j,k)=0.0
                if (cloudf(i,j,k) >  1.0) then
                  cloudf(i,j,k) = 1.0
                end if
              End Do
            End Do
          End Do
!
! If sea-salt aerosol or interactive DMS emissions are to be
! used then calculate surface windspeed and height information
!
      If (L_seasalt_CCN .OR. L_DMS_em_inter) Then
!
        Do j = 1,rows
          Do i = 1,row_length
            u_1(i,j) = u(i,j,1)
            v_1(i,j) = v(i,j,1)
          End Do
        End Do
!
! Tidy up at North Pole; not done at South Pole as it's land.
! Polar values are mean of 1st interior row:
        If (at_extremity(PNorth)) Then
! Start point of first interior (i.e. non-polar) row:
          i_start  = (rows - 2) * row_length + 1
! Sum over points on PEs in order along first interior row:
          Call gcg_rvecsumr(row_length*rows, row_length, i_start, 1,    &
     &                      u_1, proc_row_group, istat, u_1_mean)
          Call gcg_rvecsumr(row_length*rows, row_length, i_start, 1,    &
     &                      v_1, proc_row_group, istat, v_1_mean)
          u_1_mean=u_1_mean/global_row_length
          v_1_mean=v_1_mean/global_row_length
          Do i = 1, row_length
            u_1(i, rows) = u_1_mean
            v_1(i, rows) = v_1_mean
          End Do
        Endif
!
          Do k = 1, aero_dim3
            Do j = 1, aero_dim2
              Do i = 1, aero_dim1
                height(i,j,k) = r_theta_levels(i,j,k)                   &
     &                                          -r_theta_levels(i,j,0)
              End Do
            End Do
          End Do
!
      End If
      If (L_DMS_Ointer) Then
!
! Tidy up ocean DMS flux at North Pole, as done above for u_1,v_1.
! Polar values are mean of 1st interior row:
        If (at_extremity(PNorth)) Then
! Start point of first interior (i.e. non-polar) row:
          i_start  = (rows - 2) * row_length + 1
! Sum over points on PEs in order along first interior row:
          Call gcg_rvecsumr(row_length*rows, row_length, i_start, 1,    &
     &         dms_Ointer, proc_row_group, istat, dms_Ointer_mean)
          dms_Ointer_mean=dms_Ointer_mean/global_row_length
          Do i = 1, row_length
            dms_Ointer(i, rows) = dms_Ointer_mean
          End Do
        Endif
      Endif ! L_DMS_Ointer
!

!
      If (L_seasalt_CCN) Then
! DEPENDS ON: set_seasalt
          CALL SET_SEASALT(u_1, v_1, height, land_fract, ice_fract      &
     &                    , row_length, rows, model_levels              &
     &                    , bl_levels, sea_salt_film, sea_salt_jet)
      Else
        sea_salt_film(1,1,1) = 0.0
        sea_salt_jet(1,1,1) = 0.0
      Endif
!
! Calculate T from theta if soot, S cycle, biomass or OCFF switched on
! and droplet numbers for diffusional scavenging.
!
      IF (L_SULPC_SO2 .OR. L_SOOT .OR. L_BIOMASS .OR. L_OCFF) THEN
        Do k=1,model_levels
          Do j=1,rows
            Do i=1,row_length
              T(i,j,k)=theta(i,j,k)*exner_theta_levels(i,j,k)
            End Do
          End Do
        End Do
        Do k=1,wet_model_levels
          Do j=1,rows
            Do i=1,row_length
              rho_air=p_theta_levels(i,j,k)/                            &
     &                (R * T(i,j,k))
!
              If (L_seasalt_CCN) Then
                ss_film=sea_salt_film(i,j,k)
                ss_jet=sea_salt_jet(i,j,k)
              Else
                ss_film=0.0
                ss_jet=0.0
              End If
!
              If (L_sulphate_CCN) Then
                sulp_ait=so4_ait(i,j,k)
                sulp_acc=so4_acc(i,j,k)
                sulp_dis=so4_dis(i,j,k)
              Else
                sulp_ait=0.0
                sulp_acc=0.0
                sulp_dis=0.0
              End If
!
              If (L_biomass_CCN) Then
                agd_bmass=bmass_agd(i,j,k)
                cld_bmass=bmass_cld(i,j,k)
              Else
                agd_bmass=0.0
                cld_bmass=0.0
              End If
!
              If (L_ocff_CCN) Then
                agd_ocff=ocff_agd(i,j,k)
                cld_ocff=ocff_cld(i,j,k)
              Else
                agd_ocff=0.0
                cld_ocff=0.0
              End If
!
              If (L_USE_BIOGENIC) Then
                biogenic_ccn=BIOGENIC(i,j,k)
              Else
                biogenic_ccn=0.0
              End If
!
! DEPENDS ON: number_droplet
              n_droplet(i,j,k)=NUMBER_DROPLET(                          &
     &                 L_sulphate_CCN,.FALSE.                           &
     &                ,sulp_ait,sulp_acc,sulp_dis                       &
     &                ,L_seasalt_CCN,ss_film,ss_jet                     &
     &                ,L_USE_BIOGENIC,biogenic_ccn                      &
     &                ,L_biomass_CCN,agd_bmass,cld_bmass                &
     &                ,L_ocff_CCN,agd_ocff,cld_ocff                     &
     &                ,rho_air,snow_depth(i,j),land_fract(i,j)          &
     &                ,Ntot_land,Ntot_sea                               &
     &                 )
!
            End Do
          End Do
        End Do
      END IF
!
! Soot cycle code.
! Calculate emissions of soot from surface and chimney level, ageing
! of soot and diffusional scavenging of aged soot to soot-in-cloud.
!
      IF (L_SOOT) THEN  ! If soot modelling is included

        If (L_SOOT_SUREM) then  ! If surface soot emissions are included
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &    rows, row_length, off_x, off_y, halo_i, halo_j                &
     &,   model_levels, wet_model_levels                                &
     &,   halo_i, halo_j                                                &
     &,   r_rho_levels, r_theta_levels                                  &
     &,   theta, q , qcl , qcf , exner_rho_levels, rho                  &
     &,   SOOT_NEW(:,:,1), SOOT_SUREM, 1                                &
     &,   TIMESTEP, 1, 1, 0.0                                           &
     &    )
        End If  ! L_SOOT_SUREM
!
        If (L_SOOT_HILEM) then  ! If chimney level soot emissions
                                ! are included
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &    rows, row_length, off_x, off_y, halo_i, halo_j                &
     &,   model_levels, wet_model_levels                                &
     &,   halo_i, halo_j                                                &
     &,   r_rho_levels, r_theta_levels                                  &
     &,   theta, q , qcl , qcf , exner_rho_levels, rho                  &
     &,   SOOT_NEW(:,:,SOOT_HIGH_LEVEL), SOOT_HILEM                     &
     &,   SOOT_HIGH_LEVEL, TIMESTEP, 1, 1, 0.0                          &
     &    )
        End If  ! L_SOOT_HILEM
!
! Calculate quantity of fresh soot converted to aged soot
!
! DEPENDS ON: agesoot
        CALL AGESOOT(                                                   &
     &  row_length, rows, off_x, off_y,                                 &
     &  model_levels, timestep,                                         &
     &  soot_new,                                                       &
     &  delta_agesoot                                                   &
     &  )
!
! Calculate quantity of aged soot scavenged to cloud soot
!
! DEPENDS ON: sootdiffscav
        CALL SOOTDIFFSCAV(                                              &
     &  rows, row_length, off_x, off_y, halo_i, halo_j,                 &
     &  model_levels, wet_model_levels, timestep,                       &
     &  cloudf, qcl, qcf, p_theta_levels, t,                            &
     &  n_droplet,                                                      &
     &  soot_agd, soot_cld,                                             &
     &  delta_sootdiffscav                                              &
     &  )
!
! Update soot arrays with increments from ageing and
! diffusional scavenging
!
        Do k=1,model_levels
          Do j=1,rows
            Do i=1,row_length
              soot_new(i,j,k)=soot_new(i,j,k) - delta_agesoot(i,j,k)
              soot_agd(i,j,k)=soot_agd(i,j,k) + delta_agesoot(i,j,k)    &
     &                               - delta_sootdiffscav(i,j,k)
              soot_cld(i,j,k)=soot_cld(i,j,k)                           &
     &                               + delta_sootdiffscav(i,j,k)
            End Do
          End Do
        End Do
!
      End If  ! L_SOOT
!
! End of soot cycle code
!
! Biomass aerosol code.
! Calculate emissions of smoke from surface and high level, ageing
! of smoke and nucleation scavenging of aged smoke to smoke-in-cloud.
!
      IF (L_BIOMASS) THEN  ! If biomass smoke modelling is included

        If (L_BMASS_SUREM) then  ! If surface smoke emissions included
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &    rows, row_length, off_x, off_y, halo_i, halo_j                &
     &,   model_levels, wet_model_levels                                &
     &,   halo_i, halo_j                                                &
     &,   r_rho_levels, r_theta_levels                                  &
     &,   theta, q , qcl , qcf , exner_rho_levels, rho                  &
     &,   BMASS_NEW(:,:,1), BMASS_SUREM, 1                              &
     &,   TIMESTEP, 1, 1, 0.0                                           &
     &    )
        End If  ! L_BMASS_SUREM
!
        If (L_BMASS_HILEM) then  ! If high level smoke emissions
                                 ! are included
!
! Check that the range of emission levels is correctly specified
! and, if so, emit equal quantities of smoke on all model levels
! between bmass_high_level_1 and bmass_high_level_2.
!
          If (bmass_high_level_1  >   bmass_high_level_2 .or.           &
     &        bmass_high_level_1  >   model_levels .or.                 &
     &        bmass_high_level_2  >   model_levels) then
            Write(6,*) 'Aero_Ctl: Invalid range of biomass emission '
            Write(6,*) 'levels specified.'
            Write(6,*) 'Lowest level: ',bmass_high_level_1
            Write(6,*) 'Highest level: ',bmass_high_level_2
          Else
!
            Do j=1,rows
              Do i=1,row_length
                biomass_per_level(i,j) = bmass_hilem(i,j)/              &
     &          ((bmass_high_level_2 - bmass_high_level_1) + 1)
              End Do
            End Do
!
            Do biomass_level = bmass_high_level_1, bmass_high_level_2
!
! DEPENDS ON: trsrce
              CALL TRSRCE(                                              &
     &        rows, row_length, off_x, off_y, halo_i, halo_j            &
     &,       model_levels, wet_model_levels                            &
     &,       halo_i, halo_j                                            &
     &,       r_rho_levels, r_theta_levels                              &
     &,       theta, q , qcl , qcf , exner_rho_levels, rho              &
     &,       bmass_new(:,:,biomass_level), biomass_per_level           &
     &,       biomass_level, timestep, 1, 1, 0.0                        &
     &        )
!
            End Do  ! biomass_level
!
          End If  ! Test on emissions levels
!
        End If  ! L_BMASS_HILEM
!
! Calculate quantity of fresh smoke converted to aged smoke
!
! DEPENDS ON: agebmass
        CALL AGEBMASS(                                                  &
     &  row_length, rows, off_x, off_y,                                 &
     &  model_levels, timestep,                                         &
     &  bmass_new,                                                      &
     &  delta_agebmass                                                  &
     &  )
!
! Calculate quantity of aged smoke scavenged to cloud smoke
!
! DEPENDS ON: bmassnuclscav
        CALL BMASSNUCLSCAV(                                             &
     &  rows, row_length, off_x, off_y, halo_i, halo_j,                 &
     &  model_levels, wet_model_levels, timestep,                       &
     &  cloudf, qcl, qcf,                                               &
     &  bmass_agd, bmass_cld,                                           &
     &  delta_bmassnuclscav                                             &
     &  )
!
! Update smoke arrays with increments from ageing and
! nucleation scavenging
!
        Do k=1,model_levels
          Do j=1,rows
            Do i=1,row_length
              bmass_new(i,j,k)=bmass_new(i,j,k) - delta_agebmass(i,j,k)
!    Simulate the condensation of VOCs onto aged biomass by
!    increasing the mass transferred upon ageing by (8.75/5.4),
!    the ratio of the fraction of BC in fresh (8.75%) and aged
!    (5.4%) biomass aerosol:
              bmass_agd(i,j,k)=bmass_agd(i,j,k)                         &
     &                        + ((8.75/5.4)*delta_agebmass(i,j,k))      &
     &                               - delta_bmassnuclscav(i,j,k)
              bmass_cld(i,j,k)=bmass_cld(i,j,k)                         &
     &                               + delta_bmassnuclscav(i,j,k)
            End Do
          End Do
        End Do
!
      End If  ! L_BIOMASS
!
! End of biomass aerosol code
!
!
      IF ( L_SULPC_SO2 ) THEN
!
! Calculate cos zenith angle for SULPH2 by calling SOLPOS and SOLANG
! Calculate number of seconds since midnight to the beginning of
! the timetsep.
!
      seconds_since_midnight = Real( PREVIOUS_TIME(4) * 3600            &
     &         + PREVIOUS_TIME(5) * 60  + PREVIOUS_TIME(6))

! Calculate true latitude at each grid-point.
        Do j = 1, rows
          Do i = 1, row_length
            sin_true_latitude(i,j) = f3_at_u(i,j) / two_Omega
          End Do
        End Do
!
! DEPENDS ON: solpos
        CALL SOLPOS (PREVIOUS_TIME(7), PREVIOUS_TIME(1),                &
     &               L_CAL360, L_SEC_VAR, L_EqT, Eq_Time, Sindec, SCS)
!
! DEPENDS ON: solang
        CALL SOLANG(                                                    &
! input constants
     &                Sindec, seconds_since_midnight,                   &
     &                timestep, Eq_Time,                                &
! row and column dependent constants
     &                sin_true_latitude,                                &
     &                true_longitude,                                   &
! size variables
     &                row_length*rows,                                  &
! output fields
     &                day_fraction, COSZA2D, hour_angle )
!
!
! Calculate interactive DMS emissions if required
!
       If (L_DMS_em_inter) Then

         Do j = 1, rows
           Do i = 1, row_length
             If (L_ctile) Then
               T_surf(i, j) = Tstar_sea(i, j)
             Else
               T_surf(i, j) = Tstar(i, j)
             End If
           End Do
         End Do

! DEPENDS ON: dms_flux
         CALL DMS_FLUX (                                                &
! size variables
     &                row_length, rows                                  &
! input fields
     &,               u_1, v_1                                          &
     &,               height(:,:,1)                                     &
     &,               T_surf                                            &
     &,               land_fract                                        &
     &,               DMS_conc                                          &
! logical switches
     &,               L_DMS_Liss_Merlivat                               &
     &,               L_DMS_Wanninkhof                                  &
     &,               L_DMS_Nightingale                                 &
! output field
     &,               DMS_em_inter )

       End If
!
!
! Zero output diagnostics for full model timestep
      Do k=1,model_levels
        Do j=1,rows
          Do i=1,row_length
            MSA(i,j,k) = 0.0
            NH3_DEP(i,j,k) = 0.0
          End Do
        End Do
      End Do
!
! Calculate length of chemistry timestep and use it to control input
! of emissions and S Chemistry
!
      CHEMSTEP=timestep/CALL_CHEM_FREQ
!
        Do n=1,CALL_CHEM_FREQ
!
! Call TRSRCE to insert emissions
!
      If (L_SO2_SURFEM) Then         ! Insert surface SO2 emiss
! DEPENDS ON: trsrce
        CALL TRSRCE(                                                    &
     &  rows, row_length, off_x, off_y, halo_i, halo_j                  &
     &, model_levels, wet_model_levels                                  &
     &, halo_i, halo_j                                                  &
     &, r_rho_levels, r_theta_levels                                    &
     &, theta, q , qcl , qcf , exner_rho_levels, rho                    &
     &, SO2(:,:,1), SO2_SURFEM, 1                                       &
     &, CHEMSTEP, 1, 1, 0.0                                             &
     & )
      End If
!
      If (L_SO2_HILEM) Then          ! Insert chimney SO2 emiss
! DEPENDS ON: trsrce
        CALL TRSRCE(                                                    &
     &  rows, row_length, off_x, off_y, halo_i, halo_j                  &
     &, model_levels, wet_model_levels                                  &
     &, halo_i, halo_j                                                  &
     &, r_rho_levels, r_theta_levels                                    &
     &, theta, q , qcl , qcf , exner_rho_levels, rho                    &
     &, SO2(:,:,SO2_HIGH_LEVEL), SO2_HILEM, SO2_HIGH_LEVEL              &
     &, CHEMSTEP, 1, 1, 0.0                                             &
     & )
      End If
!
      If (L_SO2_NATEM) Then          ! Insert volcanic SO2 emiss
        Do k= 1,model_levels
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &  rows, row_length, off_x, off_y, halo_i, halo_j                  &
     &, model_levels, wet_model_levels                                  &
     &, halo_i, halo_j                                                  &
     &, r_rho_levels, r_theta_levels                                    &
     &, theta, q , qcl , qcf , exner_rho_levels, rho                    &
     &, SO2(:,:,k), SO2_NATEM(:,:,k), k                                 &
     &, CHEMSTEP, 1, 1, 0.0                                             &
     & )
        End Do
      End If
!
! Merge ocean DMS emissions with land DMS emissions
      If (L_DMS_Ointer) Then
        DMS_em_inter(:,:) = DMS_Ointer(:,:)
      End If

      If (L_DMS_em_inter .OR. L_DMS_Ointer) Then
        Do j = 1, rows
          Do i = 1, row_length
            DMS_em_combined(i, j)                                       &
     &          = (land_fract(i, j) * DMS_em_ancil(i, j))               &
     &          + ( ((1.0 - land_fract(i, j)) * DMS_em_inter(i, j))     &
     &             * (1.0 - ice_fract(i, j)) )
          End Do
        End Do
      Else
        If (L_DMS_em) Then     ! Just copy over the standard ancil
          Do j = 1, rows
            Do i = 1, row_length
              DMS_em_combined(i, j) = DMS_em_ancil(i, j)
            End Do
          End Do
        End If
      End If
!
      If (L_DMS_EM) Then             ! Insert DMS emiss (surface)
! DEPENDS ON: trsrce
        CALL TRSRCE(                                                    &
     &  rows, row_length, off_x, off_y, halo_i, halo_j                  &
     &, model_levels, wet_model_levels                                  &
     &, halo_i, halo_j                                                  &
     &, r_rho_levels, r_theta_levels                                    &
     &, theta, q , qcl , qcf , exner_rho_levels, rho                    &
     &, DMS(:,:,1), DMS_em_combined, 1                                  &
     &, CHEMSTEP, 1, 1, 0.0                                             &
     & )
      End If
!
      If (L_NH3_EM) Then             ! Insert NH3 emiss (surface)
! DEPENDS ON: trsrce
        CALL TRSRCE(                                                    &
     &  rows, row_length, off_x, off_y, halo_i, halo_j                  &
     &, model_levels, wet_model_levels                                  &
     &, halo_i, halo_j                                                  &
     &, r_rho_levels, r_theta_levels                                    &
     &, theta, q , qcl , qcf , exner_rho_levels, rho                    &
     &, NH3(:,:,1), NH3_EM, 1                                           &
     &, CHEMSTEP, 1, 1, 0.0                                             &


     & )
      End If
!
!
! DEPENDS ON: sulphr
      CALL SULPHR(                                                      &
! Arguments IN
     &             halo_i, halo_j, off_x, off_y                         &
     &,            row_length, rows                                     &
     &,            model_levels, wet_model_levels                       &
     &, theta_field_size                                                &
     &,            CHEMSTEP                                             &
     &,            CLOUDF, COSZA2D                                      &
     &,            p_theta_levels, T, q, qcl, qcf                       &
     &,            OH_CONC, H2O2_LMT, HO2_CONC, O3                      &
     &,            n_droplet                                            &
     &,            L_SULPC_DMS, L_SULPC_NEWDMS                          &
     &,            L_SULPC_OZONE                                        &
     &,            L_SULPC_SO2_O3_NONBUFFERED                           &
     &,            L_SULPC_NH3                                          &
! Arguments IN/OUT
     &,            SO2, DMS, SO4_AIT, SO4_ACC, SO4_DIS                  &
     &,            NH3, H2O2_MXR                                        &
! Arguments OUT (diagnostics)
     &,            MSA_inc                                              &
     &,            NH3_DEP_inc                                          &
     &,            F_DMS_TO_SO2, F_DMS_TO_SO4, F_DMS_TO_MSA             &
     &,            DELTAS_DRY, DELTAS_WET, DELTAS_WET_O3                &
     &,            DELTAS_TOT, DELTAS_DMS                               &
     &,            DELTAS_EVAP, DELTAS_NUCL, DELTAS_DIFFUSE             &
     &,            DELTAS_MERGE                                         &
     &,            DELTAS_COAG, PSI                                     &
     &            )
!
! Add diagnostic increments from SULPHR to total for model timestep
!
      Do k=1,model_levels
        Do j=1,rows
          Do i=1,row_length
            MSA(i,j,k)=MSA(i,j,k)+MSA_inc(i,j,k)
          NH3_DEP(i,j,k)=NH3_DEP(i,j,k)+NH3_DEP_inc(i,j,k)
          End Do
        End Do
      End Do
!
!
        End Do            ! End CALL_CHEM_FREQ loop
!
      ENDIF               ! End L_SULPC_SO2 test

!
! OCFF cycle code.
! Calculate emissions of ocff from surface and chimney level, ageing
! of ocff and nucleation scavenging of aged ocff to ocff-in-cloud.
!
      IF (L_OCFF) THEN  ! If ocff modelling is included

        If (L_OCFF_SUREM) then  ! If surface ocff emissions are included
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &    rows, row_length, off_x, off_y, halo_i, halo_j                &
     &,   model_levels, wet_model_levels                                &
     &,   halo_i, halo_j                                                &
     &,   r_rho_levels, r_theta_levels                                  &
     &,   theta, q , qcl , qcf , exner_rho_levels, rho                  &
     &,   OCFF_NEW(:,:,1), OCFF_SUREM, 1                                &
     &,   TIMESTEP, 1, 1, 0.0                                           &
     &    )
        End If  ! L_OCFF_SUREM
!
        If (L_OCFF_HILEM) then  ! If chimney level ocff emissions
                                ! are included
! DEPENDS ON: trsrce
          CALL TRSRCE(                                                  &
     &    rows, row_length, off_x, off_y, halo_i, halo_j                &
     &,   model_levels, wet_model_levels                                &
     &,   halo_i, halo_j                                                &
     &,   r_rho_levels, r_theta_levels                                  &
     &,   theta, q , qcl , qcf , exner_rho_levels, rho                  &
     &,   OCFF_NEW(:,:,OCFF_HIGH_LEVEL), OCFF_HILEM                     &
     &,   OCFF_HIGH_LEVEL, TIMESTEP, 1, 1, 0.0                          &
     &    )
        End If  ! L_OCFF_HILEM
!
! Calculate quantity of fresh ocff converted to aged ocff
!
! DEPENDS ON: ageocff
        CALL AGEOCFF(                                                   &
     &  row_length, rows, off_x, off_y,                                 &
     &  model_levels, timestep,                                         &
     &  ocff_new,                                                       &
     &  delta_ageocff                                                   &
     &  )
!
! Calculate quantity of aged ocff scavenged to cloud ocff
!
! DEPENDS ON: ocffnuclscav
        CALL OCFFNUCLSCAV(                                              &
     &  rows, row_length, off_x, off_y, halo_i, halo_j,                 &
     &  model_levels, wet_model_levels, timestep,                       &
     &  cloudf, qcl, qcf,                                               &
     &  ocff_agd, ocff_cld,                                             &
     &  delta_ocffnuclscav                                              &
     &  )
!
! Update ocff arrays with increments from ageing and
! nucleation scavenging
!
        Do k=1,model_levels
          Do j=1,rows
            Do i=1,row_length
              ocff_new(i,j,k)=ocff_new(i,j,k) - delta_ageocff(i,j,k)
              ocff_agd(i,j,k)=ocff_agd(i,j,k) + delta_ageocff(i,j,k)    &
     &                               - delta_ocffnuclscav(i,j,k)
              ocff_cld(i,j,k)=ocff_cld(i,j,k)                           &
     &                               + delta_ocffnuclscav(i,j,k)
            End Do
          End Do
        End Do
!
      End If  ! L_OCFF
!
! End of ocff cycle code
!

!
! PM10 / PM2.5 code: If any of the diagnostics on PM10 and PM2.5
! mass concentrations (item numbers 220-235, Sect. 17) is requested
! then call the subroutine that considers all aerosol species and modes, 
! to calculate their contributions to PM concentrations

      If (sf(220,17) .OR. sf(221,17) .OR. sf(222,17) .OR.               & 
          sf(223,17) .OR. sf(224,17) .OR. sf(225,17) .OR.               &  
          sf(226,17) .OR. sf(227,17) .OR. sf(228,17) .OR.               & 
          sf(229,17) .OR. sf(230,17) .OR. sf(231,17) .OR.               & 
          sf(232,17) .OR. sf(233,17) .OR. sf(234,17) .OR.               & 
          sf(235,17)) then

! DEPENDS ON: pm10_pm2p5
        Call pm10_pm2p5 (                                               &
          off_x, off_y,                                                 &
          row_length, rows,                                             &
          model_levels,                                                 &
          salt_dim1, salt_dim2, salt_dim3,                              &
          L_SULPC_SO2, L_SOOT, L_BIOMASS,                               & 
          L_OCFF, L_USE_BIOGENIC, L_DUST,                               &   
          L_seasalt_CCN, L_USE_SEASALT_DIRECT,                          & 
          L_USE_SEASALT_INDIRECT, L_USE_SEASALT_AUTOCONV,               &   
          p_theta_levels, T,                                            &
          SO4_AIT, SO4_ACC, SOOT_NEW, SOOT_AGD, BMASS_NEW, BMASS_AGD,   &
          OCFF_NEW, OCFF_AGD, BIOGENIC, sea_salt_film, sea_salt_jet,    &
          DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5,        &
          PM10, PM2p5,                                                  &
          PM10_SO4, PM2p5_SO4,                                          &
          PM10_BC, PM2p5_BC,                                            &
          PM10_BB, PM2p5_BB,                                            &
          PM10_OCFF, PM2p5_OCFF,                                        &
          PM10_SOA, PM2p5_SOA,                                          &
          PM10_SS, PM2p5_SS,                                            &
          PM10_DUST, PM2p5_DUST)
      
      End if

!
! Call diagnostic routine
!
! DEPENDS ON: diagnostics_aero
      Call diagnostics_aero(                                            &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels                           &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      timestep                                   &
     &,                      at_extremity                               &
!
     &,                      L_SULPC_SO2, L_DMS_em                      &
     &,                      L_SULPC_DMS, L_SULPC_NEWDMS                &
     &,                      L_SULPC_OZONE, L_SULPC_NH3                 &
     &,                      L_SOOT                                     &
     &,                      MSA, NH3_DEP                               &
     &,                      DMS_em_combined                            &
     &,                      DELTAS_DMS                                 &
     &,                      F_DMS_TO_SO2                               &
     &,                      F_DMS_TO_SO4                               &
     &,                      F_DMS_TO_MSA                               &
     &,                      DELTAS_DRY                                 &
     &,                      DELTAS_WET                                 &
     &,                      DELTAS_WET_O3                              &
     &,                      DELTAS_EVAP                                &
     &,                      DELTAS_NUCL                                &
     &,                      DELTAS_DIFFUSE                             &
     &,                      DELTAS_COAG                                &
     &,                      DELTAS_MERGE                               &
     &,                      PSI                                        &
     &,                      PM10, PM2p5                                &
     &,                      PM10_SO4, PM2p5_SO4                        &
     &,                      PM10_BC, PM2p5_BC                          &
     &,                      PM10_BB, PM2p5_BB                          &
     &,                      PM10_OCFF, PM2p5_OCFF                      &
     &,                      PM10_SOA, PM2p5_SOA                        &
     &,                      PM10_SS, PM2p5_SS                          &
     &,                      PM10_DUST, PM2p5_DUST                      &
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork17                                                      &
     &  )
      RETURN
      END SUBROUTINE AERO_CTL
!
