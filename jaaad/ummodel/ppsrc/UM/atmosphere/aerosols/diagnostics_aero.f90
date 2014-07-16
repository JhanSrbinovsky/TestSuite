
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      Subroutine diagnostics_aero(                                      &
     &                       row_length, rows, model_levels             &
     &,                      wet_model_levels                           &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      timestep                                   &
     &,                      at_extremity                               &
     &,                      L_SULPC_SO2, L_DMS_em                      &
     &,                      L_SULPC_DMS, L_SULPC_NEWDMS                &
     &,                      L_SULPC_OZONE, L_SULPC_NH3                 &
     &,                      L_SOOT                                     &
     &,                      MSA, NH3_DEP                               &
     &,                      DMS_emiss                                  &
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
     &,                      PM10,      PM2p5                           &
     &,                      PM10_SO4,  PM2p5_SO4                       &
     &,                      PM10_BC,   PM2p5_BC                        &
     &,                      PM10_BB,   PM2p5_BB                        &
     &,                      PM10_OCFF, PM2p5_OCFF                      &
     &,                      PM10_SOA,  PM2p5_SOA                       &
     &,                      PM10_SS,   PM2p5_SS                        &
     &,                      PM10_DUST, PM2p5_DUST                      &
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork                                                        &
     &  )
!
!----------------------------------------------------------------------
! Purpose:  Calculates diagnostics from section 17 AERO_CTL2 routine
!           and outputs them.
!
! Current owner of code:      C E Johnson
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
! System Components covered:
!
! System task:
!
! Documentation:  UMDP 20
!
!-----------------------------------------------------------------------
!
      Implicit None
!
! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Parameters
!
      Integer                                                           &
     &  row_length                                                      &
                            ! number of points on a row
     &, rows                                                            &
                            ! number of rows in a theta field
     &, n_rows                                                          &
                            ! number of rows in a v field
     &, model_levels                                                    &
                            ! number of model levels
     &, wet_model_levels    ! number of model levels where moisture
!
      Integer                                                           &
     &  global_row_length                                               &
                            ! NUMBER OF points on a global row
     &, global_rows                                                     &
                            ! NUMBER OF global rows
     &, me                                                              &
                            ! Processor number
     &, halo_i                                                          &
                            ! size of large halo in x direction
     &, halo_j                                                          &
                            ! size of large halo in y direction
     &, off_x                                                           &
                            ! size of small halo in x direction
     &, off_y                                                           &
                            ! size of small halo in y direction
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, g_datastart (3,0:n_proc-1)

      Real                                                              &
     &  timestep
!
      Logical                                                           &
     &  L_SULPC_SO2                                                     &
                          ! T if S Cycle on
     &, L_SULPC_DMS                                                     &
                          ! T if DMS included
     &, L_DMS_em                                                        &
                          ! T if DMS emissions used
     &, L_SULPC_NEWDMS                                                  &
                          ! T if new DMS scheme used (requires OZONE)
     &, L_SULPC_OZONE                                                   &
                          ! T if OZONE field present
     &, L_SULPC_NH3                                                     &
                          ! T if NH3 field present
     &, L_SOOT            ! T if SOOT modelling on
!
! Arguments with intent IN/OUT (diagnostics):
      REAL                                                              &
     & MSA(row_length,rows,model_levels)                                &
                                                  !mmr S in MSA
     &,NH3_DEP(row_length,rows,model_levels)                            &
                                                  !NH3 depleted
     &,DMS_emiss(row_length,rows)                 !DMS emiss (kgSm-2s-1)
      REAL                                                              &
     & DELTAS_DMS(row_length,rows,model_levels)                         &
     &,F_DMS_TO_SO2(row_length,rows,model_levels)                       &
     &,F_DMS_TO_SO4(row_length,rows,model_levels)                       &
     &,F_DMS_TO_MSA(row_length,rows,model_levels)                       &
     &,DELTAS_DRY(row_length,rows,model_levels)                         &
     &,DELTAS_WET(row_length,rows,model_levels)                         &
     &,DELTAS_WET_O3(row_length,rows,model_levels)                      &
     &,DELTAS_EVAP(row_length,rows,model_levels)                        &
     &,DELTAS_NUCL(row_length,rows,model_levels)                        &
     &,DELTAS_DIFFUSE(row_length,rows,model_levels)                     &
     &,DELTAS_COAG(row_length,rows,model_levels)                        &
     &,DELTAS_MERGE(row_length,rows,model_levels)                       &
     &,PSI(row_length,rows,model_levels)                                &
     &,PM10(row_length,rows,model_levels)                               &
                                                  !PM10 (ug m-3)
     &,PM2p5(row_length,rows,model_levels)                              &  
                                                  !PM2.5 (ug m-3)
     &,PM10_SO4 (row_length, rows, model_levels)                        &
     &,PM2p5_SO4(row_length, rows, model_levels)                        &
                                    !Sulphate contributions to PM concs.
     &,PM10_BC (row_length, rows, model_levels)                         &
     &,PM2p5_BC(row_length, rows, model_levels)                         &
                                    !Black carbon contrib. to PM concs.
     &,PM10_BB (row_length, rows, model_levels)                         &
     &,PM2p5_BB(row_length, rows, model_levels)                         &
                                    !Biomass aerosol contrib to PM concs.
     &,PM10_OCFF (row_length, rows, model_levels)                       &
     &,PM2p5_OCFF(row_length, rows, model_levels)                       &
                                    !OCFF contributions to PM concs.
     &,PM10_SOA (row_length, rows, model_levels)                        &
     &,PM2p5_SOA(row_length, rows, model_levels)                        &
                                    !SOA contributions to PM concs.
     &,PM10_SS (row_length, rows, model_levels)                         &
     &,PM2p5_SS(row_length, rows, model_levels)                         &
                                    !Sea-salt contributions to PM concs.
     &,PM10_DUST (row_length, rows, model_levels)                       &
     &,PM2p5_DUST(row_length, rows, model_levels)
                                    !Dust contributions to PM concs.
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
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------

! Diagnostic variables
      Real                                                              &
     & STASHwork(*)     ! STASH workspace for section 17

! Local variables
      Integer                                                           &
     & i, j, level, k                                                   &
                             ! loop counters
     &,    icode             ! Return code  =0 Normal exit  >1 Error

      Integer sect,item    ! STASH section, item no.s
      Parameter (sect = 17) !  for aero_ctl (S Cycle or soot)

      Real work_3d(row_length,rows,model_levels) ! work space

      Character*80  cmessage

      Character(*) RoutineName
      Parameter ( RoutineName='diagnostics_aero')

      Integer                                                           &
     &  im_index        ! internal model index

! External routines
      External                                                          &
     &  copydiag_3d                                                     &
     &, copydiag                                                        &
     &, ereport

      icode = 0 ! Initialise error status
      im_index = internal_model_index(atmos_im)
!
! Copy diagnostic information to STASHwork for STASH processing
!
! Write MSA to STASH if DMS included
!
      If (L_SULPC_DMS) Then
!
        item = 203                          !MSA
        IF (icode <= 0 .and. sf(item,sect)) THEN
!
! Convert to flux per sec
          Do level=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                work_3d(i,j,level) = MSA(i,j,level)/timestep
              End Do
            End Do
          End Do
!
! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(item,sect,im_index)),          &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

          If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 203)"//cmessage
          End if
!
        END IF
!
      End If
!
!
! Write NH3_DEP to STASH if NH3 included
!
      If (L_SULPC_NH3) Then
!
        item = 204                          !NH3_DEP
        IF (icode <= 0 .and. sf(item,sect)) THEN
!
! Convert to flux per sec
          Do level=1,model_levels
            Do j=1,rows
              Do i=1,row_length
                work_3d(i,j,level) = NH3_DEP(i,j,level)/timestep
              End Do
            End Do
          End Do
!
! DEPENDS ON: copydiag_3d
          Call copydiag_3d (stashwork(si(item,sect,im_index)),          &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

          If (icode >  0) Then
            cmessage=": error in copydiag_3d(item 204)"//cmessage
          End if
!
        END IF
!
      End If         ! End L_SULPC_NH3 condn
!
!
!
! Diagnose DMS emissions if requested
!
      If (L_DMS_em) Then
!
        item = 205                          !DMS emissions
        IF (icode <= 0 .and. sf(item,sect)) THEN
!
! DEPENDS ON: copydiag
          Call copydiag (stashwork(si(item,sect,im_index)),             &
     &        DMS_emiss,                                                &
     &        row_length,rows,0,0,0,0, at_extremity,                    &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)

          If (icode >  0) Then
            cmessage=": Error in copydiag (item 205)"//cmessage
          End if
!
        END IF
!
      End If         ! End L_DMS_em condn

      IF (L_SULPC_SO2 .AND. L_SULPC_DMS) THEN

        item = 206
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DMS(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 206)"//cmessage
          ENDIF

        ENDIF

        item = 207
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! get the fraction and convert to flux per second
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DMS(i,j,level) *            &
     &              F_DMS_TO_SO2(i,j,level) / timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 207)"//cmessage
          ENDIF

        ENDIF

        item = 208
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! get the fraction and convert to flux per second
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DMS(i,j,level) *            &
     &             F_DMS_TO_SO4(i,j,level) *                            &
     &             (1.0E00 - PSI(i,j,level))/ timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 208)"//cmessage
          ENDIF

        ENDIF

        item = 209
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! get the fraction and convert to flux per second
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DMS(i,j,level) *            &
     &             F_DMS_TO_SO4(i,j,level) *                            &
     &             PSI(i,j,level) / timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 209)"//cmessage
          ENDIF

        ENDIF

        item = 210
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DRY(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 210)"//cmessage
          ENDIF

        ENDIF

        item = 211
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_WET(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 211)"//cmessage
          ENDIF

        ENDIF

        item = 212
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_WET_O3(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 212)"//cmessage
          ENDIF

        ENDIF

        item = 213
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_EVAP(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 213)"//cmessage
          ENDIF

        ENDIF

        item = 214
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_NUCL(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 214)"//cmessage
          ENDIF

        ENDIF

        item = 215
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DIFFUSE(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 215)"//cmessage
          ENDIF

        ENDIF

        item = 216
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_COAG(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 216)"//cmessage
          ENDIF

        ENDIF

        item = 217
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! get the fraction and convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DRY(i,j,level)*             &
     &             (1.0E00-PSI(i,j,level))/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 217)"//cmessage
          ENDIF

        ENDIF

        item = 218
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! get the fraction and convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_DRY(i,j,level)*             &
     &             PSI(i,j,level)/timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 218)"//cmessage
          ENDIF

        ENDIF

        item = 219
        IF (icode <= 0 .and. sf(item,sect)) THEN

          ! convert to flux per sec
          DO level=1, model_levels
            DO j=1, rows
              DO i=1, row_length
                work_3d(i,j,level) = DELTAS_MERGE(i,j,level)            &
     &             /timestep
              ENDDO
            ENDDO
          ENDDO

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(stashwork(si(item,sect,im_index)),           &
     &        work_3d,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        icode,cmessage)
          IF(icode  >   0) THEN
            cmessage=": error in copydiag_3d(item 219)"//cmessage
          ENDIF

        ENDIF


      ENDIF
      
! ---------------------------------------------------------------------------
!
!     Diagnose PM10, PM2.5 and the contributions of the different
!     aerosols species to them if requested. Note that PM10 & PM2.5
!     can be calculated as long as any of the aerosol species is used. 
     
! Diagnose PM10 if requested:
!
      item = 220
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10,                                                     &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 220)"//cmessage
        End if
      END IF
   
! Diagnose PM2.5 if requested:
!
      item = 221
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5,                                                    &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 221)"//cmessage
        End if
      END IF

! Diagnose PM10 & PM2.5 concs. due to different aerosol species
! if requested:
!
      item = 222
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_SO4,                                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 222)"//cmessage
        End if
      END IF
!
      item = 223
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_SO4,                                                &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 223)"//cmessage
        End if
      END IF
!
      item = 224
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_BC,                                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 224)"//cmessage
        End if
      END IF
!
      item = 225
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_BC,                                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 225)"//cmessage
        End if
      END IF
!
      item = 226
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_BB,                                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 226)"//cmessage
        End if
      END IF
!
      item = 227
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_BB,                                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 227)"//cmessage
        End if
      END IF
!
      item = 228
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_OCFF,                                                &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 228)"//cmessage
        End if
      END IF
!
      item = 229
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_OCFF,                                               &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 229)"//cmessage
        End if
      END IF
!
      item = 230
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_SOA,                                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 230)"//cmessage
        End if
      END IF
!
      item = 231
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_SOA,                                                &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 231)"//cmessage
        End if
      END IF
!
      item = 232
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_SS,                                                  &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 232)"//cmessage
        End if
      END IF
!
      item = 233
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_SS,                                                 &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 233)"//cmessage
        End if
      END IF
!
      item = 234
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM10_DUST,                                                &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 234)"//cmessage
        End if
      END IF
!
      item = 235
      IF (icode <= 0 .and. sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
              PM2p5_DUST,                                               &
              row_length,rows,model_levels,0,0,0,0, at_extremity,       &
              stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
              stash_levels,num_stash_levels+1,                          &
              atmos_im,sect,item,                                       &
              icode,cmessage)
        If (icode > 0) then
          cmessage=": error in copydiag_3d (item 235)"//cmessage
        End if
      END IF

! ---------------------------------------------------------------------------

      IF (icode > 0) THEN
! DEPENDS ON: ereport
        CALL ereport ('DIAGNOSTICS_AERO', icode, cmessage)
      END IF
     
      RETURN
      END SUBROUTINE diagnostics_aero
!
