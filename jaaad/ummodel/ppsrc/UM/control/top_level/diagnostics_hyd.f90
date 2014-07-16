
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine diagnostics_hyd

      Subroutine diagnostics_hyd(                                       &
     &                       row_length, rows, model_levels             &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      land_points, dsm_levels                    &
!  Add inland basin outflow to arguments
     &,                      land_index,inlandout_atm                   &
     &,                      smc, surf_roff, sub_surf_roff              &
     &,                      snow_depth_land, snow_melt, canopy, t_soil &
     &,                      soil_layer_moisture                        &
     &,                      ntiles, snomlt_surf_htf, sthu, sthf        &
     &,                      tot_tfall, snow_tile, melt_tile            &
     &,                      rgrain, land_sea_mask                      &
     &,                      dun_roff, drain, qbase, qbase_zw           &
     &,                      fch4_wetl, fexp, gamtot, ti_mean, ti_sig   &
     &,                      fsat, fwetl, zw, sthzw                     &
     &,                      timestep                                   &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork                                                        &
     &     )

! Description:
!   Calculates hydrology-related diagnostics (held in STASH section 8).
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Intercepted arrays and diagnostic arrays are input from the
!   hydrology routines called previously. Each diagnostic is simply
!   copied into the STASHwork array to be passed on to STASH for
!   output processing.
!
!   Diagnostics currently available (in order calculated):
!   Item  Description
!    208  Soil moisture content
!     23  Snow depth
!    201  Snow melt
!    234  Surface run-off
!    235  Sub-surface run-off
!    225  Soil temperature
!    223  Soil layer moisture
!    204  Surface run-off (accumulated)
!    205  Sub-surface run-off (accumulated)
!    209  Canopy water content
!    202  Land snow melt heat flux
!    231  Land snow melt rate
!    233  Canopy throughfall rate
!    236  Snow amount on tiles
!    237  Snow melt rate on tiles
!    238  Snow grain size on tiles
!    229  Unfrozen soil moisture fraction
!    230  Frozen soil moisture fraction
!    239  Baseflow
!    240  Dunne Runoff
!    241  Baseflow from water table layer
!    242  Wetland methane flux
!    252  Drainage out of "nshyd"th model layer
!    245  Inland basin outflow on atmos grid
!    252  Drainage out of bottom "nshyd"th soil layer (currently layer 4).

! History:
! Version   Date     Comment
! ----     -------     -------
! 5.0  30/11/99 Original version. J-C Thil.
! 5.1  09/12/99 Add error trapping. Replace rmdi_pp by rmdi.
!               Rick Rawlins
! 5.1  27/03/00 Add new diagnostics. J-C Thil.
!  5.1  28/04/00   Access land points via proper 2D indexes instead
!                  of the much objectionable out of bounds first
!                  index & second index set to 1. JC Thil
! 5.2  12/10/00 Add diagnostics 201,209. D Matthews
! 5.3  26/09/01 Add additional hydrology diagnostics.    M. Best
! 5.3  27/07/01 Correct extraction to be over soil not model levels
!               for (8,223),(8,225). R Rawlins
! 5.5  04/02/03 Add additional hydrology diagnostics.    Nic Gedney
!   6.2   21/2/06  Re-route outflow from inland basins to soil moisture
!                  P. Falloon

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, number_format                                                   &
                         ! switch controlling number format diagnostics
                         ! are written out in. See PP_WRITE for details.
     &, model_domain     ! indicator as to model type, ie global, lam

      Integer                                                           &
     &  global_row_length                                               &
                            !IN. NUMBER OF points on a global row
     &, global_rows                                                     &
                            !IN. NUMBER OF global rows
     &, me                                                              &
                            !IN. Processor number
     &, halo_i                                                          &
                            !IN. size of large halo in x direction
     &, halo_j                                                          &
                            !IN. size of large halo in y direction
     &, off_x                                                           &
                            !IN. size of small halo in x direction
     &, off_y                                                           &
                            !IN. size of small halo in y direction
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        &
     &, g_datastart(3,0:n_proc-1)                                       &
     &, land_points                                                     &
                    ! No.of land points being processed, can be 0.
     &, dsm_levels                                                      &
     &, ntiles      ! No. of land-surface tiles ( MOSES II )

      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP


      Real                                                              &
     &  timestep

! Primary Arrays used in all models
      Integer                                                           &
     &  land_index(land_points)      ! set from land_sea_mask

      Real                                                              &
     &  snow_depth_land(land_points)                                    &
                                     !
     &, snow_melt(land_points)                                          &
                                ! snowmelt (kg/m2/s)
     &, canopy (land_points)                                            &
                              ! canopy water content (kg/m2)
     &, smc(land_points)                                                &
                            ! available soil moisture in the
!                                 rootzone (kg/m2).
     &, surf_roff(land_points)                                          &
                                ! surface runoff (kg/m2/s).
     &, sub_surf_roff(land_points)                                      &
                                    ! sub-surface runoff
! Declare inland basin outflow variable
     &, inlandout_atm(land_points)                                      &
                                    !inland basin outflow

     &, t_soil(land_points,dsm_levels)                                  &
     &, soil_layer_moisture(land_points,dsm_levels)                     &
     &, snomlt_surf_htf(row_length, rows)                               &
     &, sthu(land_points,dsm_levels)                                    &
                                       ! Unfrozen soil moisture
!                content of each layer as a fraction of saturation
     &, sthf(land_points,dsm_levels)                                    &
                                       ! Frozen soil moisture content
!                of each layer as a fraction of saturation
     &, tot_tfall(land_points)                                          &
                                       ! total throughfall (kg/m2/s)
     &, snow_tile(land_points,ntiles)                                   &
                                       ! Lying snow on tiles (kg/m2)
     &, melt_tile(land_points,ntiles)                                   &
                                       ! Snowmelt on tiles (kg/m2/s)
     &, rgrain(land_points,ntiles)                                      &
                                       ! Snow grain size (microns)

! Additional variables required for large-scale hydrology:
     &, qbase(land_points)                                              &
                                    ! Baseflow (kg/m2/s).
     &, dun_roff(land_points)                                           &
                                    ! Dunne runoff (kg/m2/s).
     &, qbase_zw(land_points)                                           &
                                    ! Baseflow from Zw layer (kg/m2/s).
     &, drain(land_points)                                              &
                                    ! Drainage out of "nshyd" later
!                                   ! (kg/m2/s).
     &, fch4_wetl(land_points)                                          &
                                    ! Wetland methane flux (kg C/m2/s).
     &, TI_MEAN(LAND_POINTS)                                            &
                                    ! Mean topographic index.
     &, TI_SIG(LAND_POINTS)                                             &
                                    ! Std. dev. of topographic index.
     &, FEXP(LAND_POINTS)                                               &
                                    ! Decay factor in Sat. Conductivity
!                                   !   in water table layer.
     &, GAMTOT(LAND_POINTS)                                             &
                                    ! Integrated complete Gamma
!                                   !   function.
     &, FSAT(LAND_POINTS)                                               &
                                    ! Surface saturation fraction.
     &, FWETL(LAND_POINTS)                                              &
                                    ! Wetland fraction.
     &, ZW(LAND_POINTS)                                                 &
                                    ! Water table depth (m).
     &, STHZW(LAND_POINTS)          ! Saturation fraction in water
!                                   !   table layer.


      Logical                                                           &
     &  land_sea_mask(row_length, rows)

      Integer                                                           &
     & PSLEVEL                                                          &
                     !  loop counter for pseudolevels
     &,PSLEVEL_OUT                                                      &
                     !  index for pseudolevels sent to STASH
     &,LEVEL                                                            &
                     !  loop counter for levels
     &,LEVEL_OUT     !  index for levels sent to STASH

      Logical                                                           &
     & PLLTILE(NTILES)                                                  &
                          ! pseudolevel list for surface types
     &,LIST(DSM_LEVELS)   ! level list for soil levels


! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Local variables

      Integer                                                           &
     &  i, j, k, l                                                      &
     &,    icode                ! Return code  =0 Normal exit  >1 Error

      Character*80  cmessage
      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_hyd')

      Integer                                                           &
     &  im_index        ! internal model index

      Real                                                              &
     &  interp_data(row_length,rows)                                    &
     &, interp_data_3(row_length,rows,dsm_levels)

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

! Diagnostic variables
       Real                                                             &
     &  STASHwork(*)    ! STASH workspace

      External                                                          &
     &  copydiag, copydiag_3d                                           &
     &  ,Ereport

! ----------------------------------------------------------------------
! Section 1.  Initialisation.
! ----------------------------------------------------------------------


      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)

! ----------------------------------------------------------------------
! Soil moisture content
! ----------------------------------------------------------------------
! Item 8 208  smc

      If (sf(208,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = smc(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag (STASHwork(si(208,8,im_index)),interp_data,      &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,208,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 208)"
            goto 9999
         End if

      End if


! ----------------------------------------------------------------------
! Snow depth
! ----------------------------------------------------------------------
! Item 8 023 snow_depth_land

      If (sf(023,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = snow_depth_land(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(023,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,023,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 023)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Snow melt
! ----------------------------------------------------------------------
! Item 201

      If (sf(201,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = snow_melt(l) * timestep
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(201,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,201,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 201)"
            goto 9999
         End if

      End if


! ----------------------------------------------------------------------
! Surface run-off.
! ----------------------------------------------------------------------
! Item 234  surf_roff

      If (sf(234,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = surf_roff(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(234,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,234,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 234)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Sub-surface run-off.
! ----------------------------------------------------------------------
! Item 235  sub_surf_roff

      If (sf(235,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = sub_surf_roff(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(235,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,235,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 235)"
            goto 9999
         End if

      End if


! ----------------------------------------------------------------------
!  Soil temperature
! ----------------------------------------------------------------------
! Item 8 225 t_soil

      If (sf(225,8)) Then

         Do k = 1, dsm_levels
            Do j= 1, rows
               Do i = 1, row_length
                  interp_data_3(i,j,k) = rmdi
               End Do
            End Do

            Do l = 1, land_points
               j=(land_index(l)-1)/row_length + 1
               i=land_index(l) - (j-1)*row_length
               interp_data_3(i,j,k) = t_soil(l,k)
            End Do

         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(225,8,im_index)),interp_data_3,  &
     &        row_length,rows,dsm_levels,0,0,0,0,                       &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,225,8,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,8,225,                                           &
     &        icode,cmessage)



         If (icode  >   0) then
            cmessage="Error in copydiag_3d( item 225)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Soil layer moisture
! ----------------------------------------------------------------------
! Item 8 223 soil_layer_moisture

      If (sf(223,8)) Then

         Do k = 1, dsm_levels
            Do j= 1, rows
               Do i = 1, row_length
                  interp_data_3(i,j,k) = rmdi
               End Do
            End Do
             Do l = 1, land_points
               j=(land_index(l)-1)/row_length + 1
               i=land_index(l) - (j-1)*row_length
               interp_data_3(i,j,k) =                                   &
     &              soil_layer_moisture(l,k)
             End Do
         End Do

! DEPENDS ON: copydiag_3d
         Call copydiag_3d(STASHwork(si(223,8,im_index)),interp_data_3,  &
     &        row_length,rows,dsm_levels,0,0,0,0,                       &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,223,8,im_index)),len_stlist,           &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,8,223,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 223)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Surface run-off (accumulated)
! ----------------------------------------------------------------------
! Item 204  surf_roff

      If (sf(204,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do i = 1, land_points
            interp_data(land_index(i),1) = surf_roff(i) * timestep
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(204,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,204,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 204 )"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Sub-surface run-off (accumulated)
! ----------------------------------------------------------------------
! Item 205  sub_surf_roff

      If (sf(205,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do i = 1, land_points
            interp_data(land_index(i),1) = sub_surf_roff(i)             &
     &     *  timestep
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(205,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,205,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 205)"
            goto 9999
         endif
      endif

! ----------------------------------------------------------------------
! Canopy water content
! ----------------------------------------------------------------------
! Item 209

      If (sf(209,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do i = 1, land_points
            interp_data(land_index(i),1) = canopy(i)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(209,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,209,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 209)"
            goto 9999
         endif
      endif


! ----------------------------------------------------------------------
! Land snow melt heat flux (W/m2)
! ----------------------------------------------------------------------
! Item 202 SNOMLT_SURF_HTF

      IF (SF(202,8)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(202,8,im_index)),                              &
     &       snomlt_surf_htf,                                           &
     &       land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Land snow melt heat rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 231 snow_melt

      IF (SF(231,8)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(231,8,im_index)),                              &
     &       snow_melt,                                                 &
     &       land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Canopy throughfall rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 233 tot_tfall

      IF (SF(233,8)) THEN
! DEPENDS ON: from_land_points
        CALL FROM_LAND_POINTS (                                         &
     &      STASHWORK(SI(233,8,im_index)),                              &
     &       tot_tfall,                                                 &
     &       land_sea_mask,row_length*rows,land_points)
      END IF


! ----------------------------------------------------------------------
! Snow amount on tiles (Kg/m2)
! ----------------------------------------------------------------------
! Item 236 snow_tile

      IF (SF(236,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,236,8,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_pseudo_list(item 236 = snow_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(236,8,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),snow_tile(1,PSLEVEL_OUT),            &
     &           land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Snow melt rate on tiles (Kg/m2)
! ----------------------------------------------------------------------
! Item 237 melt_tile

      IF (SF(237,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,237,8,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_pseudo_list(item 237 = melt_tile)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(237,8,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),melt_tile(1,PSLEVEL_OUT),            &
     &           land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Snow grain size on tiles
! ----------------------------------------------------------------------
! Item 238 rgrain

      IF (SF(238,8)) THEN
! DEPENDS ON: set_pseudo_list
        CALL SET_PSEUDO_LIST(NTILES,LEN_STLIST,                         &
     &       STLIST(1,STINDEX(1,238,8,im_index)),                       &
     &       PLLTILE,STASH_PSEUDO_LEVELS,NUM_STASH_PSEUDO,              &
     &       ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_pseudo_list(item 238 = rgrain)"
            goto 9999
        END IF
        PSLEVEL_OUT=0
        DO PSLEVEL=1,NTILES
          IF (PLLTILE(PSLEVEL)) THEN
            PSLEVEL_OUT=PSLEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(238,8,im_index)+(PSLEVEL_OUT-1)            &
     &           *row_length*rows),rgrain(1,PSLEVEL_OUT),               &
     &           land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF



! ----------------------------------------------------------------------
! Unfrozen soil moisture fraction
! ----------------------------------------------------------------------
! Item 229 sthu

      IF (SF(229,8)) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(DSM_LEVELS,LEN_STLIST,                     &
     &       STLIST(1,STINDEX(1,229,8,im_index)),                       &
     &       LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_levels_list(item 229 = sthu)"
            goto 9999
        END IF
        LEVEL_OUT=0
        DO LEVEL=1,DSM_LEVELS
          IF(LIST(LEVEL)) THEN
            LEVEL_OUT=LEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(229,8,im_index)+(LEVEL_OUT-1)              &
     &           *row_length*rows),sthu(1,LEVEL_OUT),                   &
     &           land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF


! ----------------------------------------------------------------------
! Frozen soil moisture fraction
! ----------------------------------------------------------------------
! Item 230 sthf

      IF (SF(230,8)) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(DSM_LEVELS,LEN_STLIST,                     &
     &       STLIST(1,STINDEX(1,230,8,im_index)),                       &
     &       LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
            cmessage=                                                   &
     &    "daghyd  : error in set_levels_list(item 230 = sthf)"
            goto 9999
        END IF
        LEVEL_OUT=0
        DO LEVEL=1,DSM_LEVELS
          IF(LIST(LEVEL)) THEN
            LEVEL_OUT=LEVEL_OUT+1
! DEPENDS ON: from_land_points
            CALL FROM_LAND_POINTS (                                     &
     &          STASHWORK(SI(230,8,im_index)+(LEVEL_OUT-1)              &
     &           *row_length*rows),sthf(1,LEVEL_OUT),                   &
     &           land_sea_mask,row_length*rows,land_points)
          END IF
        END DO
      END IF


! ----------------------------------------------------------------------
! Baseflow
! ----------------------------------------------------------------------
! Item 239  baseflow

      If (sf(239,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = qbase(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(239,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,239,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 239)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Dunne Runoff
! ----------------------------------------------------------------------
! Item 240 Dunne Runoff

      If (sf(240,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = dun_roff(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(240,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,240,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 240)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Baseflow from water table layer
! ----------------------------------------------------------------------
! Item 241  baseflow from water table layer

      If (sf(241,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = qbase_zw(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(241,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,241,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 241)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Wetland methane flux
! ----------------------------------------------------------------------
! Item 242  Wetland methane flux

      If (sf(242,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fch4_wetl(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(242,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,242,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 242)"
            goto 9999
         End if

      End if
      If (sf(243,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = ti_mean(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(243,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,243,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 243)"
            goto 9999
         End if
      End if
      If (sf(244,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = ti_sig(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(244,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,244,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 244)"
            goto 9999
         End if
      End if
      If (sf(251,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fexp(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(251,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,251,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 251)"
            goto 9999
         End if
      End if
      If (sf(246,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = gamtot(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(246,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,246,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 246)"
            goto 9999
         End if
      End if
      If (sf(247,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fsat(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(247,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,247,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 247)"
            goto 9999
         End if
      End if
      If (sf(248,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = fwetl(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(248,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,248,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 248)"
            goto 9999
         End if
      End if
      If (sf(249,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = zw(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(249,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,249,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 249)"
            goto 9999
         End if
      End if
      If (sf(250,8)) Then
         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do
         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = sthzw(l)
         End Do
! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(250,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,250,                                           &
     &        icode,cmessage)
         If (icode  >   0) then
            cmessage="Error in copydiag( item 250)"
            goto 9999
         End if
      End if
! ----------------------------------------------------------------------
! Drainage out of  "nshyd"th model layer
! ----------------------------------------------------------------------
! Item 252  drainage

      If (sf(252,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = drain(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(252,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,252,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 252)"
            goto 9999
         End if

      End if

! ----------------------------------------------------------------------
! Output inland basin outflow on atmosphere grid

! --------------------------------------------------------------------
! Inland basin outflow (atmos grid)
! --------------------------------------------------------------------
! Item 245  Inland basin outflow

      If (sf(245,8)) Then

         Do j = 1, rows
            Do i = 1, row_length
               interp_data(i,j) = rmdi
            End Do
         End Do

         Do l = 1, land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            interp_data(i,j) = inlandout_atm(l)
         End Do

! DEPENDS ON: copydiag
         Call copydiag(STASHwork(si(245,8,im_index)),interp_data,       &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,8,245,                                           &
     &        icode,cmessage)

         If (icode  >   0) then
            cmessage="Error in copydiag( item 245)"
            goto 9999
         End if

      End if
! -------------------------------------------
 9999 continue
      If(icode /= 0) Then
! DEPENDS ON: ereport
        Call Ereport(RoutineName,icode,Cmessage)
      Endif

      Return
      END SUBROUTINE diagnostics_hyd

