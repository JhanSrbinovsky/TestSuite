

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine NI_gwd_ctl

      Subroutine NI_gwd_ctl (                                           &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y, global_row_length, global_rows    &
     &, proc_row_group, at_extremity, n_proc, n_procx, n_procy          &
     &, neighbour, g_rows, g_row_length, g_datastart, me                &

! model dimensions.
     &, row_length, rows, n_rows, land_points                           &
     &, model_levels                                                    &

! Model switches
     &, model_domain, gw_kay, gwd_frc                                   &
     &, Ltimer, l_gwd, L_use_ussp                                       &
     &, l_taus_scale, l_fix_gwsatn, l_gwd_40km, l_ussp_opaque           &
     &, sat_scheme, gwd_fsat                                            &

! trig arrays
     &, sin_theta_longitude, sin_theta_latitude                         &

! in coordinate information
     &, r_rho_levels, r_theta_levels                                    &
     &, eta_theta_levels, eta_rho_levels                                &

! in time stepping information.
     &, timestep, timestep_number                                       &

! diagnostic info
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork                                                        &

! in data fields.
     &, u, v                                                            &
     &, land_sea_mask, p_layer_boundaries                               &
     &, rho, theta_latest, sd_orog_land, orog_grad_xx_land              &
     &, orog_grad_xy_land, orog_grad_yy_land, land_index                &

! in/out
     &, R_u, R_v                                                        &

! error information
     &, Error_code  )

! purpose: Interface to Atmospheric Physics GWD Schemes.
!         Scheme 1: Orographic flow blocking and gravity wave scheme.
!         Scheme 2: Non-orographic ultra-simple spectral gw scheme.
!
!
! current code owner: S. Webster
!
! history:
! Version   Date     Comment
! ----   -------   -------
! 5.0  30/11/99 Original version. J-C Thil.
! 5.1  09/12/99 Call diagnostics routine only when diagnostics
!               requested + comment changes. Rick Rawlins
! 5.1 15/03/00 Correct diagnostic dimensioning and remove now
!              superfluous lists. R Rawlins
! 5.1  28/02/00 Provide explicit model increments as STASH output
!               diagnostics. R Rawlins
! 5.1  10/04/00 Remove setting GWD gw_p_ref in namelist. S. Webster
!
! 5.2  15/11/00 Replace call to gwave with call to gwd_glue to allow
!               3B or 4A schemes to be called. Update call to daggwd
!               to include all new diagnostics. Allocate memory
!               for all diagnostics only when required. Pass 4A
!               switches down from cruntimc.
!                                                         S. Webster
! 5.3  17/10/01 Changes required for Single Column Model
!                                             Z. Gardner
! 5.3  16/10/01 Pass l_use_ussp thru' routine rather than l_gwd_wake
!                                                         S. Webster
! 5.3  11/10/01 Add arguments to diagnostic_gwd to calculate
!               orographic standard deviation. D.M. Goddard
! 5.4  28/08/02 Add arrays for numerical limiter diagnostics. S.Webster
!
!  5.4   05/09/02   Add diagnostics for spectral (non-orographic)
!               gravity wave forcing scheme (GW_USSP).      Adam Scaife
!
! 5.5  25/02/03 Remove 3B GWD code and hence call 4A code from this
!               deck rather than via glue routines. Also remove other
!               redundant arrays.                            S. Webster
! 6.2  21/02/06 Pass l_taus_scale and l_fix_gwsatn through to GWAVE4A
!               and l_ussp_opaque thru' to GW_USSP         S. Webster
!
! 6.2  22/06/06 Fix dependency on diagnostics_gwd for SCM.
!               P.Selwood
! 6.2  21/02/06 Add orographic surface stress and
!               sigma_xx, xy and yy diagnostics.             S. Webster
! 6.3  11/10/06 Fix unmatched endif in SCMA extraction.  R Barnes
!
! 6.4  19/05/06 Correction of USSP diagnostic flux output.
!                                                   A.C. Bushell
!
! code description:
!   language: fortran 77 + cray extensions
!   this code is written to umdp3 programming standards.

      Implicit None

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
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, global_rows                                                     &
                           ! NUMBER OF global rows
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
     &  at_extremity(4)  ! Indicates if this processor is at north, sout
                         ! east or west of the processor grid


! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, model_levels                                                    &
     &, land_points

! Model switches
      Integer                                                           &
     &  model_domain                                                    &
     &, sat_scheme

      Logical                                                           &
     &  Ltimer                                                          &
                 ! true then output some timing information
     &, L_gwd                                                           &
                    ! switch for the orographic GWD scheme
     &, L_use_ussp                                                      &
                    ! switch for the non-orographic USSP scheme
     &, L_taus_scale                                                    &
                      ! switch to scale orog surf stress by Froude No.
     &, L_fix_gwsatn                                                    &
                      ! switch to invoke minor bug fixes in gwsatn4a
     &, L_gwd_40km                                                      &
                      ! switch to turn off orographic GWD above 40km     
     &, L_ussp_opaque ! switch to change lid condition in gw_ussp

! model parameters
      Real                                                              &
     &  timestep                                                        &
     &, GW_kay                                                          &
     &, gwd_frc                                                         &
     &, gwd_fsat

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

! Diagnostics info
       Real                                                             &
     & STASHwork(*) ! STASH workspace for section 6 (GW Drag)


! Data arrays

      Integer                                                           &
     &  land_index (land_points)      ! set from land_sea_mask

      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &      model_levels)    ! density*r*r

      Real                                                              &
     &  theta_latest(row_length, rows, model_levels)

      Real                                                              &
     &  p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, sd_orog_land (land_points)                                      &
                                   ! orog/qrparm.orog.stdev
     &, orog_grad_xx_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmaxx
     &, orog_grad_xy_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmaxy
     &, orog_grad_yy_land(land_points) ! orog/qrparm.orog.sigmayy

      logical                                                           &
     &  land_sea_mask(row_length, rows)

! Co-ordinate arrays
      Real                                                              &
           ! local vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

      Real                                                              &
     &  sin_theta_longitude (row_length, rows)                          &
     &, sin_theta_latitude  (row_length, rows)

! time information for current timestep
      Integer                                                           &
     &  timestep_number


! arguments with intent in/out. ie: input variables changed on output.

! arguments with intent out. ie: output variables.

      Real                                                              &
     &  R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        model_levels)

      Integer                                                           &
     &  Error_code

! local variables
      Logical                                                           &
     &  stress_ud_on                                                    &
     &, stress_vd_on                                                    &
     &, stress_ud_satn_on                                               &
     &, stress_vd_satn_on                                               &
     &, stress_ud_wake_on                                               &
     &, stress_vd_wake_on                                               &
     &, du_dt_satn_on                                                   &
     &, dv_dt_satn_on                                                   &
     &, du_dt_wake_on                                                   &
     &, dv_dt_wake_on                                                   &
     &, GWSPEC_EFLUX_ON                                                 &
     &, GWSPEC_SFLUX_ON                                                 &
     &, GWSPEC_WFLUX_ON                                                 &
     &, GWSPEC_NFLUX_ON                                                 &
     &, GWSPEC_EWACC_ON                                                 &
     &, GWSPEC_NSACC_ON                                                 &
     &, u_s_d_on                                                        &
     &, v_s_d_on                                                        &
     &, nsq_s_d_on                                                      &
     &, fr_d_on                                                         &
     &, bld_d_on                                                        &
     &, bldt_d_on                                                       &
     &, num_lim_d_on                                                    &
     &, num_fac_d_on                                                    &
     &, tausx_d_on                                                      &
     &, tausy_d_on                                                      &
     &, taus_scale_d_on                                                 &
     &, L_u_incr_gwd                                                    &
     &, L_v_incr_gwd

      Integer                                                           &
     & i,j,k       ! loop counters


! Local data arrays

! Allocatable arrays for diagnostic variables - required to save memory
! when diagnostic not requested
      Real,Dimension(:,:,:),Allocatable::                               &
     &  u_incr_diagnostic                                               &
                              ! u wind increment for STASH
     &, v_incr_diagnostic                                               &
                              ! v wind increment for STASH
     &, stress_ud                                                       &
     &, stress_vd                                                       &
     &, stress_ud_satn                                                  &
     &, stress_vd_satn                                                  &
     &, stress_ud_wake                                                  &
     &, stress_vd_wake                                                  &
     &, du_dt_satn                                                      &
     &, dv_dt_satn                                                      &
     &, du_dt_wake                                                      &
     &, dv_dt_wake                                                      &
     &, GWSPEC_EFLUX                                                    &
     &, GWSPEC_SFLUX                                                    &
     &, GWSPEC_WFLUX                                                    &
     &, GWSPEC_NFLUX                                                    &
     &, GWSPEC_EWACC                                                    &
     &, GWSPEC_NSACC

      Real,Dimension(:,:),Allocatable::                                 &
     &  u_s_d                                                           &
     &, v_s_d                                                           &
     &, nsq_s_d                                                         &
     &, fr_d                                                            &
     &, bld_d                                                           &
     &, bldt_d                                                          &
     &, num_lim_d                                                       &
     &, num_fac_d                                                       &
     &, tausx_d                                                         &
     &, tausy_d                                                         &
     &, taus_scale_d


! Diagnostic land_point array sizes
      Integer                                                           &
     & points_stress_ud                                                 &
     &,points_stress_vd                                                 &
     &,points_stress_ud_satn                                            &
     &,points_stress_vd_satn                                            &
     &,points_stress_ud_wake                                            &
     &,points_stress_vd_wake                                            &
     &,points_du_dt_satn                                                &
     &,points_dv_dt_satn                                                &
     &,points_du_dt_wake                                                &
     &,points_dv_dt_wake                                                &
     &,points_u_s_d                                                     &
     &,points_v_s_d                                                     &
     &,points_nsq_s_d                                                   &
     &,points_fr_d                                                      &
     &,points_bld_d                                                     &
     &,points_bldt_d                                                    &
     &,points_num_lim_d                                                 &
     &,points_num_fac_d                                                 &
     &,points_tausx_d                                                   &
     &,points_tausy_d                                                   &
     &,points_taus_scale_d

! Local arrays holding information to be passed between physics
! routines.

! Diagnostics controlled by Diagnostic switches


! External Routines:
      External timer, g_wave
      External diagnostics_gwd

! ----------------------------------------------------------------------
! Section GWD.1 Set stash diagnostic switches
! ----------------------------------------------------------------------

! General case of the atmosphere model, ie : with stash.
      GWSPEC_EFLUX_ON = sf(101,6)
      GWSPEC_SFLUX_ON = sf(102,6)
      GWSPEC_WFLUX_ON = sf(103,6)
      GWSPEC_NFLUX_ON = sf(104,6)
      GWSPEC_EWACC_ON = sf(105,6)
      GWSPEC_NSACC_ON = sf(106,6)
      L_u_incr_gwd      = sf(185,6)
      L_v_incr_gwd      = sf(186,6)
      stress_ud_on      = sf(201,6)
      stress_vd_on      = sf(202,6)
      du_dt_satn_on     = sf(207,6)
      dv_dt_satn_on     = sf(208,6)
      u_s_d_on          = sf(214,6)
      v_s_d_on          = sf(215,6)
      nsq_s_d_on        = sf(216,6)
      fr_d_on           = sf(217,6)
      bld_d_on          = sf(218,6)
      bldt_d_on         = sf(222,6)
      stress_ud_satn_on = sf(223,6)
      stress_vd_satn_on = sf(224,6)
      stress_ud_wake_on = sf(227,6)
      stress_vd_wake_on = sf(228,6)
      du_dt_wake_on     = sf(231,6)
      dv_dt_wake_on     = sf(232,6)
      num_lim_d_on      = sf(233,6)
      num_fac_d_on      = sf(234,6)
      tausx_d_on        = sf(235,6)
      tausy_d_on        = sf(236,6)
      taus_scale_d_on   = sf(237,6)


! ----------------------------------------------------------------------
! Section GWD.1
! ----------------------------------------------------------------------

      If (error_code  ==  0 ) Then
! DEPENDS ON: timer
        If (Ltimer) Call timer ('GW Drag ',3)
! Save R_u    before updating
      If ( L_u_incr_gwd) Then  ! STASHflag set
        Allocate ( u_incr_diagnostic(row_length,rows,model_levels) )

        Do k=1,model_levels
         Do j=1,rows
          Do i=1,row_length
            u_incr_diagnostic(i,j,k) = R_u(i,j,k)
          Enddo ! i
         Enddo ! j
        Enddo ! k

      Endif                    ! on STASHflag

! Save R_v    before updating
      If ( L_v_incr_gwd) Then  ! STASHflag set
        Allocate ( v_incr_diagnostic(row_length,n_rows,model_levels) )

        Do k=1,model_levels
         Do j=1,n_rows
          Do i=1,row_length
            v_incr_diagnostic(i,j,k) = R_v(i,j,k)
          Enddo ! i
         Enddo ! j
        Enddo ! k

      Endif                    ! on STASHflag

      If ( stress_ud_on ) Then  ! STASHflag set
        Allocate ( stress_ud(row_length,rows,0:model_levels) )
        points_stress_ud = land_points
      Else
        points_stress_ud = 1
      Endif

      If ( stress_vd_on ) Then  ! STASHflag set
        Allocate ( stress_vd(row_length,n_rows,0:model_levels) )
        points_stress_vd = land_points
      Else
        points_stress_vd = 1
      Endif

      If ( du_dt_satn_on ) Then  ! STASHflag set
        Allocate ( du_dt_satn(row_length,rows,model_levels) )
        points_du_dt_satn = land_points
      Else
        points_du_dt_satn = 1
      Endif

      If ( dv_dt_satn_on ) Then  ! STASHflag set
        Allocate ( dv_dt_satn(row_length,n_rows,model_levels) )
        points_dv_dt_satn = land_points
      Else
        points_dv_dt_satn = 1
      Endif

      If ( u_s_d_on ) Then  ! STASHflag set
        Allocate ( u_s_d(row_length,rows) )
        points_u_s_d = land_points
      Else
        points_u_s_d = 1
      Endif

      If ( v_s_d_on ) Then  ! STASHflag set
        Allocate ( v_s_d(row_length,rows) )
        points_v_s_d = land_points
      Else
        points_v_s_d = 1
      Endif

      If ( nsq_s_d_on ) Then  ! STASHflag set
        Allocate ( nsq_s_d(row_length,rows) )
        points_nsq_s_d = land_points
      Else
        points_nsq_s_d = 1
      Endif

      If ( fr_d_on ) Then  ! STASHflag set
        Allocate ( fr_d(row_length,rows) )
        points_fr_d = land_points
      Else
        points_fr_d = 1
      Endif

      If ( bld_d_on ) Then  ! STASHflag set
        Allocate ( bld_d(row_length,rows) )
        points_bld_d = land_points
      Else
        points_bld_d = 1
      Endif

      If ( bldt_d_on ) Then  ! STASHflag set
        Allocate ( bldt_d(row_length,rows) )
        points_bldt_d = land_points
      Else
        points_bldt_d = 1
      Endif

      If ( num_lim_d_on ) Then  ! STASHflag set
        Allocate ( num_lim_d(row_length,rows) )
        points_num_lim_d = land_points
      Else
        points_num_lim_d = 1
      Endif

      If ( num_fac_d_on ) Then  ! STASHflag set
        Allocate ( num_fac_d(row_length,rows) )
        points_num_fac_d = land_points
      Else
        points_num_fac_d = 1
      Endif

      If ( stress_ud_satn_on ) Then  ! STASHflag set
        Allocate ( stress_ud_satn(row_length,rows,0:model_levels) )
        points_stress_ud_satn = land_points
      Else
        points_stress_ud_satn = 1
      Endif

      If ( stress_vd_satn_on ) Then  ! STASHflag set
        Allocate ( stress_vd_satn(row_length,n_rows,0:model_levels) )
        points_stress_vd_satn = land_points
      Else
        points_stress_vd_satn = 1
      Endif

      If ( stress_ud_wake_on ) Then  ! STASHflag set
        Allocate ( stress_ud_wake(row_length,rows,0:model_levels) )
        points_stress_ud_wake = land_points
      Else
        points_stress_ud_wake = 1
      Endif

      If ( stress_vd_wake_on ) Then  ! STASHflag set
        Allocate ( stress_vd_wake(row_length,n_rows,0:model_levels) )
        points_stress_vd_wake = land_points
      Else
        points_stress_vd_wake = 1
      Endif

      If ( du_dt_wake_on ) Then  ! STASHflag set
        Allocate ( du_dt_wake(row_length,rows,model_levels) )
        points_du_dt_wake = land_points
      Else
        points_du_dt_wake = 1
      Endif

      If ( dv_dt_wake_on ) Then  ! STASHflag set
        Allocate ( dv_dt_wake(row_length,n_rows,model_levels) )
        points_dv_dt_wake = land_points
      Else
        points_dv_dt_wake = 1
      Endif

      If ( tausx_d_on ) Then  ! STASHflag set
        Allocate ( tausx_d(row_length,rows) )
        points_tausx_d = land_points
      Else
        points_tausx_d = 1
      Endif

      If ( tausy_d_on ) Then  ! STASHflag set
        Allocate ( tausy_d(row_length,n_rows) )
        points_tausy_d = land_points
      Else
        points_tausy_d = 1
      Endif

      If ( taus_scale_d_on ) Then  ! STASHflag set
        Allocate ( taus_scale_d(row_length,rows) )
        points_taus_scale_d = land_points
      Else
        points_taus_scale_d = 1
      Endif

      IF ( GWSPEC_EFLUX_ON ) Then  ! STASHflag set
        Allocate (GWSPEC_EFLUX(row_length,rows,model_levels))
      ENDIF
      IF ( GWSPEC_SFLUX_ON ) Then  ! STASHflag set
        Allocate (GWSPEC_SFLUX(row_length,n_rows,model_levels))
      ENDIF
      IF ( GWSPEC_WFLUX_ON ) Then  ! STASHflag set
        Allocate (GWSPEC_WFLUX(row_length,rows,model_levels))
      ENDIF
      IF ( GWSPEC_NFLUX_ON ) Then  ! STASHflag set
        Allocate (GWSPEC_NFLUX(row_length,n_rows,model_levels))
      ENDIF
      IF ( GWSPEC_EWACC_ON ) Then  ! STASHflag set
       Allocate (GWSPEC_EWACC(row_length,rows,model_levels))
      ENDIF
      IF ( GWSPEC_NSACC_ON ) Then  ! STASHflag set
       Allocate (GWSPEC_NSACC(row_length,n_rows,model_levels))
      ENDIF


!-------------------------------------------------------------------
!L Section GWD.2  Call Gravity Wave Drag Scheme version 4A
!-------------------------------------------------------------------

       If ( L_gwd ) Then

! DEPENDS ON: g_wave
        CALL G_WAVE(                                                    &
     &   theta_latest, u, v, row_length, rows, n_rows,                  &
     &   off_x, off_y, halo_i, halo_j,                                  &
     &   model_domain, at_extremity, model_levels,                      &
     &   rho, r_rho_levels, r_theta_levels, sd_orog_land,               &
     &   orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land,       &
     &   land_index, land_points, timestep, gw_kay, gwd_frc,            &
     &   r_u, r_v, l_taus_scale, l_fix_gwsatn, l_gwd_40km,              &
     &   sat_scheme, gwd_fsat,                                          &
     &   stress_ud      , stress_ud_on     , points_stress_ud     ,     &
     &   stress_vd      , stress_vd_on     , points_stress_vd     ,     &
     &   stress_ud_satn , stress_ud_satn_on, points_stress_ud_satn,     &
     &   stress_vd_satn , stress_vd_satn_on, points_stress_vd_satn,     &
     &   stress_ud_wake , stress_ud_wake_on, points_stress_ud_wake,     &
     &   stress_vd_wake , stress_vd_wake_on, points_stress_vd_wake,     &
     &   du_dt_satn     , du_dt_satn_on    , points_du_dt_satn    ,     &
     &   dv_dt_satn     , dv_dt_satn_on    , points_dv_dt_satn    ,     &
     &   du_dt_wake     , du_dt_wake_on    , points_du_dt_wake    ,     &
     &   dv_dt_wake     , dv_dt_wake_on    , points_dv_dt_wake    ,     &
     &   u_s_d          , u_s_d_on         , points_u_s_d         ,     &
     &   v_s_d          , v_s_d_on         , points_v_s_d         ,     &
     &   nsq_s_d        , nsq_s_d_on       , points_nsq_s_d       ,     &
     &   fr_d           , fr_d_on          , points_fr_d          ,     &
     &   bld_d          , bld_d_on         , points_bld_d         ,     &
     &   bldt_d         , bldt_d_on        , points_bldt_d        ,     &
     &   num_lim_d      , num_lim_d_on     , points_num_lim_d     ,     &
     &   num_fac_d      , num_fac_d_on     , points_num_fac_d     ,     &
     &   tausx_d        , tausx_d_on       , points_tausx_d       ,     &
     &   tausy_d        , tausy_d_on       , points_tausy_d       ,     &
     &   taus_scale_d   , taus_scale_d_on  , points_taus_scale_d  ,     &
     &   error_code)

      End If

      If ( L_use_ussp ) Then

! DEPENDS ON: gw_ussp
      CALL GW_USSP(MODEL_LEVELS, MODEL_DOMAIN, ROWS, N_ROWS,            &
     &    OFF_X, OFF_Y, HALO_I, HALO_J, ROW_LENGTH,                     &
     &    R_RHO_LEVELS, R_THETA_LEVELS, P_LAYER_BOUNDARIES,             &
     &    R_U,R_V, L_USSP_OPAQUE,                                       &
     &    SIN_THETA_LONGITUDE, SIN_THETA_LATITUDE,                      &
     &    THETA_LATEST, RHO, TIMESTEP, U, V, AT_EXTREMITY,              &
     &    GWSPEC_EFLUX,GWSPEC_SFLUX,GWSPEC_WFLUX,GWSPEC_NFLUX,          &
     &    GWSPEC_EWACC,GWSPEC_NSACC,ETA_THETA_LEVELS,                   &
     &    GWSPEC_EFLUX_ON,GWSPEC_SFLUX_ON,GWSPEC_WFLUX_ON,              &
     &    GWSPEC_NFLUX_ON,GWSPEC_EWACC_ON,GWSPEC_NSACC_ON)

      End If


! DEPENDS ON: timer
        If (Ltimer) Call timer ('GW Drag ',4)

! ----------------------------------------------------------------------
! Section GWD.3 Call GWD diagnostics
! ----------------------------------------------------------------------

        If(sf(0,6)) Then ! diagnostics requested this timestep
! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',3)

! DEPENDS ON: diagnostics_gwd
          Call diagnostics_gwd(                                         &
     &                       row_length, rows, model_levels             &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, off_x, off_y, me           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length, g_datastart          &
     &,                      at_extremity                               &
     &,                      u, v, R_u, R_v                             &
     &,                      u_incr_diagnostic,v_incr_diagnostic        &
     &,                      stress_ud     ,  stress_vd                 &
     &,                      stress_ud_satn,  stress_vd_satn            &
     &,                      stress_ud_wake,  stress_vd_wake            &
     &,                      du_dt_satn    ,  dv_dt_satn                &
     &,                      du_dt_wake   ,  dv_dt_wake                 &
     &,                      u_s_d, v_s_d, nsq_s_d                      &
     &,                      num_lim_d, num_fac_d                       &
     &,                      fr_d, bld_d, bldt_d                        &
     &, tausx_d, tausy_d, taus_scale_d                                  &
     &, sd_orog_land, land_sea_mask, land_points                        &
     &, orog_grad_xx_land, orog_grad_xy_land, orog_grad_yy_land         &
     &, GWSPEC_EFLUX, GWSPEC_SFLUX, GWSPEC_WFLUX                        &
     &, GWSPEC_NFLUX, GWSPEC_EWACC, GWSPEC_NSACC,                       &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
! 6.1 26/10/04  Reduce continuation lines. R Hill
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,   &
     & STASH_PSEUDO_LEVELS,STASH_SERIES, STASH_SERIES_INDEX, MOS_MASK,  &
! ARGSTS end
     & STASHwork                                                        &
     & )

      If ( L_u_incr_gwd      ) Deallocate ( u_incr_diagnostic )
      If ( L_v_incr_gwd      ) Deallocate ( v_incr_diagnostic )
      If ( stress_ud_on      ) Deallocate ( stress_ud         )
      If ( stress_vd_on      ) Deallocate ( stress_vd         )
      If ( du_dt_satn_on     ) Deallocate ( du_dt_satn        )
      If ( dv_dt_satn_on     ) Deallocate ( dv_dt_satn        )
      If ( u_s_d_on          ) Deallocate ( u_s_d             )
      If ( v_s_d_on          ) Deallocate ( v_s_d             )
      If ( nsq_s_d_on        ) Deallocate ( nsq_s_d           )
      If ( fr_d_on           ) Deallocate ( fr_d              )
      If ( bld_d_on          ) Deallocate ( bld_d             )
      If ( bldt_d_on         ) Deallocate ( bldt_d            )
      If ( num_lim_d_on      ) Deallocate ( num_lim_d         )
      If ( num_fac_d_on      ) Deallocate ( num_fac_d         )
      If ( stress_ud_satn_on ) Deallocate ( stress_ud_satn    )
      If ( stress_vd_satn_on ) Deallocate ( stress_vd_satn    )
      If ( stress_ud_wake_on ) Deallocate ( stress_ud_wake    )
      If ( stress_vd_wake_on ) Deallocate ( stress_vd_wake    )
      If ( du_dt_wake_on     ) Deallocate ( du_dt_wake        )
      If ( dv_dt_wake_on     ) Deallocate ( dv_dt_wake        )
      If ( tausx_d_on        ) Deallocate ( tausx_d           )
      If ( tausy_d_on        ) Deallocate ( tausy_d           )
      If ( taus_scale_d_on   ) Deallocate ( taus_scale_d      )
      IF ( GWSPEC_EFLUX_ON   ) Deallocate ( GWSPEC_EFLUX      )
      IF ( GWSPEC_SFLUX_ON   ) Deallocate ( GWSPEC_SFLUX      )
      IF ( GWSPEC_WFLUX_ON   ) Deallocate ( GWSPEC_WFLUX      )
      IF ( GWSPEC_NFLUX_ON   ) Deallocate ( GWSPEC_NFLUX      )
      IF ( GWSPEC_EWACC_ON   ) Deallocate ( GWSPEC_EWACC      )
      IF ( GWSPEC_NSACC_ON   ) Deallocate ( GWSPEC_NSACC      )

! DEPENDS ON: timer
        If (Ltimer) Call timer ('Diags   ',4)
        Endif            ! on sf(0,6)

      End If ! on error code equal to zero

! end of routine NI_gwd_ctl
      Return
      END SUBROUTINE NI_gwd_ctl
