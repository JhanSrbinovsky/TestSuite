
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine INANCILA
!LL
!LL Purpose : Takes as input,the code defining the frequency of update
!LL           of ancillary fields as set by the user interface.
!LL           Converts them into a list of numbers of timesteps after
!LL           which each field must be updated, and calculates the
!LL           frequency with which this list must be interrogated.
!LL           Where the update interval is in months or years,
!LL           the check will be carried out each day. The physical
!LL           files required are also determined by input code,
!LL           and the headers and lookup tables are read into
!LL           the arguments FIXHD,INTHD,LOOKUP which are in
!LL           COMMON/ANCILHDA/ of calling routine INANCCTL.
!LL           Indexes for each possible ancillary field are set up in
!LL           COMMON/IXANCILA/
!LL
!LL Level 2 Control routine for CRAY YMP
!LL
!LL CW, DR      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  22/02/93  Changes to add 2 SLAB fields (STASH items 178,179)
!LL                  - to be updated from existing atmosphere files.
!LL   3.3  22/11/93  Add aerosol ancillary fields.  R T H Barnes.
!LL   3.3  21/12/93  Fix put in to prevent array 'out of bounds'
!LL                  problem in section 1.6. Problem to be investigated
!LL                  for 3.4 D. Robinson.
!LL   3.4  16/06/94  DEF CAL360 replaced by LOGICAL LCAL360
!LL                                                   S.J.Swarbrick
!LL  3.4  05/09/94  Add murk and user ancillary fields.  RTHBarnes.
!LL   3.4  22/06/94  Array 'out of bounds' problem solved. D. Robinson
!LL  3.4   11/10/94   Part of modset which sorts out some handling
!LL                   of unset data by recon_dump.
!LL                   Necessary to port model to a T3D.
!LL                   Author D.M. Goddard
!LL   3.5  24/03/95    Changed OPEN to FILE_OPEN  P.Burton
!LL  3.5  24/07/95  Check fields for updating have valid address. RTHB
!    4.0  01/09/95  Add diagnostic information to output about
!                   ozone ancillary fields and test correct ozone
!                   data provided.  D. Goddard & D. Robinson
!LL  4.0  10/10/95  Set LOOKUP(45) in ancillary files. D. Robinson.
!LL
!LL  4.0  29/09/95  Need extra rewind of namelist file. RTHBarnes.
!LL  4.0  05/08/95  Temporary solution to get round problem of
!LL                 no. of soil moisture levels being hard-wired
!LL                 to no. of deep soil temperature levels
!LL                 This causes a problem with introduction of
!LL                 Penman-Monteith BL code at 4.0 - use if test
!LL                 on number of deep soil temperature
!LL                 levels which is set to 4 for Penman-Monteith code
!LL                 (set to 3 for all other BL options)
!LL                 Permanent solution suggested for 4.1
!LL                 search on C**** for comments
!LL                 J Smith
!LL  4.0  06/01/96  SI array received for two internal models (atmos
!LL                 and slab) in argument list. Hardwire processing of
!LL                 slab ancillary field (item code 177) to use
!LL                 SI_SLAB. D. Robinson
!    4.1  03/05/96  Use READHEAD to read in ancillary file headers.
!                   D. Robinson
!    4.1  18/06/96  Changes to cope with changes in STASH addressing
!                   Author D.M. Goddard.
!LL  4.1  22/05/96  Call new CANC* comdecks. Use new arrays in
!LL                 CANCFLDA. Cater for new sulphur ancillary files.
!LL                 Remove hardwired fix for slab ancillary fields
!LL                 introduced at 4.0 D. Robinson.
!LL  4.4  28/07/97  Add LAMIPII to namelist for special updating of
!LL                 ice in AMIP II runs. R A Stratton
!LL  4.4  16/09/97  Set number of headers for multi-pseudo-level
!LL                 ancillary fields for surface and vegetation types.
!LL                                              Richard Betts
!LL  4.4  09/09/97  New namelist UPANCA for updating information.
!LL                 D. Robinson.
!LL  4.4  10/09/97  Check calendar indicator in Anc File. D Robinson.
!    4.5  22/10/98  Set LEVELS array for new user hulti-layer
!                   ancillary fields
!                   Author D.M Goddard
!LL  4.5  19/01/98  Remove SOIL_VARS and VEG_VARS. D. Robinson.
!LL  4.5  05/05/98  Improve error message for missing files. R. Rawlins
!LL  5.1  20/07/00  Remove potential OOB reference from STASHANCIL
!LL                                                     D.Robinson
!LL  5.2  09/03/01  Initialise FIELDCODE to zero. D.Robinson
!    5.3  23/01/02  Update vertical levels checking for ozone files.
!                   Dave Robinson
!LL  5.3  19/06/01  Added code for tropopause-based ozone. D.Tan
!    5.4  27/08/02  Update vertical levels checking for murk and
!                   mulit_level ancillary files. D Robinson.
!    6.1  18/08/04  Fix possible out of bounds error. S.Wilson
!    6.1  27/07/04  Correct dimensioning of LEVDEPC. D Robinson
!    6.2  22/08/05  Remove RECON def. P.Selwood
!LL
!LL System components covered : C710
!LL
!LL System task : C7
!LL
!LL Documentation : Unified Model Documentation Paper No C7
!LL                 Version No 4  dated 15/06/90
!LLEND
      SUBROUTINE INANCILA(LEN_FIXHD,LEN_INTHD,LEN_REALHD,               &
                                                           !Intent (In)
     &                    LEN1_LEVDEPC,LEN2_LEVDEPC,                    &
     &                    FIXHD,INTHD,REALHD,LOOKUP,                    &
     &                    A_FIXHD,A_REALHD,A_LEVDEPC,                   &
     &                    NDATASETS,NLOOKUPS,FTNANCIL,                  &
     &                    LOOKUP_START,LEN1_LOOKUP,ROW_LENGTH,          &
     &                    P_ROWS,U_ROWS,P_LEVELS,                       &
     &                    TR_LEVELS,ST_LEVELS,SM_LEVELS,                &
     &                    OZONE_LEVELS,tpps_ozone_levels,TITLE,         &
!kdcorbin, 08/10 - added variables to call
     &                    SI,NITEM,NSECT,N_MODEL,                       &
     &                    SI_ATMOS,SI_SLAB,SILEN,                       &
     &                    ANCILLARY_STEPS,STEPS_PER_HR,                 &
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &                    ICODE,CMESSAGE,LCAL360)         ! Intent (Out)

      USE CSENARIO_MOD
      IMPLICIT NONE

      LOGICAL LCAL360  ! Logical switch for 360-day calendar

      INTEGER                                                           &
     &        LEN_FIXHD,                                                &
                               ! Length of header blocks in ancillary
!                              ! data sets
     &        LEN_INTHD,                                                &
                               !
     &        LEN_REALHD,                                               &
                               !
     &        LEN1_LEVDEPC,                                             &
                               ! Dimension of LEVDEPC in model
     &      LEN2_LEVDEPC                                                &
     &     ,ANCILLARY_STEPS,                                            &
     &        STEPS_PER_HR


      INTEGER                                                           &
     &        NDATASETS,                                                &
                               ! No of physical files
     &        NLOOKUPS,                                                 &
                               ! No of lookups required(set by User I.)
     &                IOUNIT,                                           &
     &        FTNANCIL(NDATASETS),                                      &
                                   ! Fortran nos of physical files
     &        LOOKUP_START(NDATASETS),                                  &
                                      !start of each individual lookup
!                                     !in overall LOOKUP array
     &        LEN1_LOOKUP,                                              &
                               ! Length of PP header
     &        ROW_LENGTH,                                               &
                               ! Atmosphere model dimensions
     &        P_ROWS,                                                   &
                               ! No. of rows for pressure-type variables
     &        U_ROWS,                                                   &
                               ! No. of rows for wind-type variables
     &        P_LEVELS,                                                 &
                               ! No. of pressure levels
     &        TR_LEVELS,                                                &
                               ! No. of tracer levels
     &        FILE_LEVELS,                                              &
                               ! Number of levels of data in files
!                              ! contining multi-level data.
     &        ST_LEVELS,                                                &
                               ! No. of soil temperature levels
     &        SM_LEVELS,                                                &
                               ! No. of soil moisture levels
     &        OZONE_LEVELS,                                             &
                                 ! No. of ozone levels
     &        tpps_ozone_levels                                         &
                                 ! No of ozone levs in TppsOzon dataset

!      For atmos only runs SI_SLAB is a copy of SI_ATMOS
!      SI_SLAB is only used in SLAB runs.

     &       ,SILEN                                                     &
                                ! Length for SI_ATMOS/SLAB arrays
     &       ,SI_ATMOS(SILEN)                                           &
                                ! ) STASHin addresses of atmos and
     &       ,SI_SLAB(SILEN)    ! ) slab ancillary fields.

     !Added SI array - kdcorbin, 05/10
      INTEGER :: NITEM,NSECT,N_MODEL
      INTEGER :: SI(NITEM,0:NSECT,N_MODEL)

      !Added filename for ancillary data - kdcorbin, 05/10
      INTEGER, PARAMETER :: Max_Filename_Len=120
      CHARACTER*Max_Filename_Len :: AncFileName
      INTEGER :: len_anc_filename

      CHARACTER*80 TITLE(NDATASETS) ! Titles of each dataset

      Logical :: l_vert_mismatch    ! T : Vertical levels mismatch

      INTEGER                                                           &
     &        FIXHD(LEN_FIXHD,NDATASETS),                               &
                                         ! Overall Fixed header array
     &        A_FIXHD(LEN_FIXHD),                                       &
                                  ! Fixed header for Dump
     &        INTHD(LEN_INTHD,NDATASETS),                               &
                                         ! Overall Integer header array
     &        LOOKUP(LEN1_LOOKUP,NLOOKUPS),                             &
                                           !Overall Lookup array
     &        ICODE            ! Return code =0 Normal Exit
!                              !             >0 Error

      REAL                                                              &
     &      REALHD(LEN_REALHD,NDATASETS),                               &
                                         !
     &      A_REALHD(LEN_REALHD),                                       &
                                 !
     &      A_LEVDEPC(LEN1_LEVDEPC,LEN2_LEVDEPC),                       &
     &      LEVDEPC( (P_LEVELS+1)*4 ) ! Space to hold level dependent
                                      ! constants from Anc File

      CHARACTER*100                                                     &
     &        CMESSAGE         ! Out error message if I>0
      Character (Len=*), Parameter :: RoutineName = 'INANCILA'

      INTEGER :: isec, iind  !kdcorbin, 05/10

! Comdecks:----------------------------------------------------------
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
! COMDECK PPXLOOK
! Description:
!
!   Declares ppxref look-up arrays used by the UM and associated
!    arrays and parameters.
!   Comdecks CSUBMODL,CPPXREF must be *CALLed before this
!    comdeck
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       May. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.0       Dec. 95   Replace dynamic dim ppxRecs with
!                     NUM_DIAG_MAX in PPXC   N. Farnon
! 4.1       July 96   *CALL VERSION introduced - NUM_DIAG_MAX made
!                      equal to NDIAGP.
!                     NUM_USR_DIAG_MAX increased from 200 to 300
!                      (just in case).
! 4.4       03/11/97  Removed MKPPXRF *DEF references. K Rogers
! 4.4       04/11/97  Changed -RECON def line to allow for other small
!                     execs which had used the RECON def. K Rogers
!
! Declarations:

! Global parameters:
! VERSION STASH parameter definitions
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Mar. 95   Original code.  S.J.Swarbrick
! 4.0                                 S.J.Swarbrick
! 4.1       Apr. 96   Rationalise MDI  S.J.Swarbrick
!  4.1  29/02/96  Increase OUTFILE_E.  RTHBarnes.
!  4.2  27/11/96  mpp code : Increase NELEMP   P.Burton
!  4.3  04/06/97  Increase NELEMP for D1 addressing S.D.Mullerworth
!  4.4  04/06/97  Increase NELEMP for sampling offset. S.D.Mullerworth
!  4.5  28/01/98  Increade NELEMP for mpp code.   P.Burton
!  4.5  18/09/98  Modify name of VERSION common block to stop potential
!                 clashes with Fortran variable names          P.Burton
!  4.5  30/09/98  Increase NRECDP from 600 to 800. D. Robinson.
!  5.2  29/01/01  OUTFILE_E changed. Adam Clayton
!  5.5  20/02/03  Increased size of STASH_SET.  P.Dando
!  6.1  03/08/04  Increase size of NPSLEVP and NPSLISTP
!                 (Pseudo Levels)  Anthony A. Dickinson
!  6.1  04/08/04  Increase size of NDIAGP W Roseblade.
!  6.2  31/03/06  Increase size of NDIAGP again.
!                 R Sempers (frpz)
!  6.2  06/04/06  Increased size of OUTFILE_E   T. Edwards
!  6.2  03/02/06  Increase NRECDP to 1500. T Johns
!
      ! Max. no. of STASH sections  per internal model (44 in practice)
      INTEGER,PARAMETER :: NSECTP=99
      ! Max. no. of STASH items per section
      INTEGER,PARAMETER :: NITEMP=999
      ! Max. no. of STASH list records (prognostic + diagnostic)
      INTEGER,PARAMETER :: NRECDP=1500
      ! Max. no. of output times tables in STASHC
      INTEGER,PARAMETER :: NTIMEP=100
      ! Max. no. of time profiles in STASHC
      INTEGER,PARAMETER :: NPROFTP=100
      ! Max. no. of domain profiles/levels lists in STASHC (used for
      ! both)
      INTEGER,PARAMETER :: NPROFDP=100
      ! Max. total no. of time series in STASHC
      INTEGER,PARAMETER :: NTimSerP=1500
      ! Max. no. time series per domain profile
      INTEGER,PARAMETER :: tsdp=250
      ! Max. no. of useage profiles in STASHC
      INTEGER,PARAMETER :: NPROFUP=40
      ! Max. no. of levels in a levels list
      INTEGER,PARAMETER :: NLEVP=50
      ! Max. no. of pseudo levels in a  pseudo levels list
      INTEGER,PARAMETER :: NPSLEVP=100
      ! Max. no. of pseudo levels lists in STASHC
      INTEGER,PARAMETER :: NPSLISTP=100
      ! Max. no. non-blank records in PPXREF file
      INTEGER,PARAMETER :: NDIAGP=2600
      INTEGER,PARAMETER :: NDIAGPM=NRECDP  ! Same as NRECDP
      INTEGER,PARAMETER :: NELEMP=33
      INTEGER,PARAMETER :: NLEVP_S=NLEVP*6+1
      INTEGER,PARAMETER :: NLEVLSTSP=NPROFDP
      INTEGER,PARAMETER :: NMEANP=4  ! No. of meaning periods
      ! OUTFILE_S, OUTFILE_L and OUTFILE_E must be consistent with
      ! NUNITS and NUNITS_LEN in file CHSUNITS.
      ! Ranges of output file numbers
      INTEGER,PARAMETER :: OUTFILE_S=20
      INTEGER,PARAMETER :: OUTFILE_E=161
      INTEGER,PARAMETER :: OUTFILE_L=OUTFILE_E-OUTFILE_S+1
!Global scalar:
      CHARACTER(LEN=80) :: STASH_SET     !Names of stasets files
!Common block:
      COMMON/common_VERSION/ STASH_SET
! VERSION end
! No. of STASH items per section
      INTEGER      PPXREF_ITEMS
        PARAMETER (PPXREF_ITEMS    =NITEMP)
! No. of STASH sections per internal model
      INTEGER      PPXREF_SECTIONS
        PARAMETER (PPXREF_SECTIONS =NSECTP-55)
! Max. number of non-null records in ppxref file (>1200)
      INTEGER      NUM_DIAG_MAX
        PARAMETER (NUM_DIAG_MAX    =NDIAGP)
! Max. number of user-defined ppxref records allowed
      INTEGER      NUM_USR_DIAG_MAX
        PARAMETER (NUM_USR_DIAG_MAX=450)

! No. of ppxref records read into PPXI,PPXC (for dyn. allocation)
      INTEGER      ppxRecs

! Global arrays:
! ppxref look-up arrays
      INTEGER   PPXI(ppxRecs,PPXREF_CODELEN)
      CHARACTER PPXC(NUM_DIAG_MAX,PPXREF_CHARLEN)
! Arrays for temporary storage of user-ppxref records -
!   used to transfer these records from STASH_PROC into U_MODEL
      INTEGER   PPXI_U(NUM_USR_DIAG_MAX,PPXREF_CODELEN)
      CHARACTER PPXC_U(NUM_USR_DIAG_MAX,PPXREF_CHARLEN)
! Array of flags to indicate origin of ppxref record
! 'P' for ppxref file; 'U' for user-stash master file
      CHARACTER OriginFlag(NUM_DIAG_MAX)
! Array of indices to identify which ppxref record corresponds to
!   any given row of PPXI, PPXC
      INTEGER   RowIndex(NUM_DIAG_MAX)
! Pointer array for PPXI, PPXC arrays
      INTEGER PPXPTR                                                    &
     & (N_INTERNAL_MODEL    ,0:PPXREF_SECTIONS ,PPXREF_ITEMS)

! Common block:
      COMMON/PPX_INT/ RowIndex,PPXI_U
      COMMON/PPX_CHA/ OriginFlag,PPXC_U
! - End --------------------------------------------------------------
! MODEL Defines model-dependent quantities used by data addressing and
! STASH
!
! Files CSUBMODL and VERSION must be included before this one
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Mar. 95   Original code.  S.J.Swarbrick
! 4.1       Apr. 96   Generalisation and incorporation of
!                      wave model     S.J.Swarbrick
! 4.2       28/11/96  mpp code : Added variables
!                                global_LPRIM and global_LDUMP to
!                                provide relevant information for
!                                 global (dump)data.  P.Burton
! 4.5       29/07/98  Remove redundant code. Processing for Boundary
!                     files moved to INTF_CTL. D. Robinson.
! 5.3       19/06/01  Add ZonAvTppsOzone for processing of
!                     tropopause-based ozone.     Dave Tan
! 6.1       30/03/04  Put free tracers into section 33.  R Barnes
! 6.1       29/09/04  Increase size of OAFLD and OASFLDID
!                     to correspond to UMUI changes
!                       Anthony A. Dickinson
! 6.2       12/08/05  Remove RECON defs and fix continuations.
!                     P.Selwood
! 6.2       10/11/05  Set up section 34 for UKCA tracers.  R Barnes
!
      INTEGER, PARAMETER :: AASSETS    = 9
      INTEGER, PARAMETER :: MEAD_TYPES = 4
      INTEGER, PARAMETER :: A_MAX_TRVARS=150 !Max.no.of tracers allowed
      INTEGER, PARAMETER :: A_MAX_UKCAVARS=150 ! Max.no.of UKCA allowed
      INTEGER, PARAMETER :: MAX_AOBS=100

      REAL :: H_A_EWSPACE
      REAL :: H_A_NSSPACE
      REAL :: H_A_FIRSTLAT
      REAL :: H_A_FIRSTLONG
      REAL :: H_A_POLELAT
      REAL :: H_A_POLELONG

      INTEGER :: H_A_GROUP
      INTEGER :: H_OROG_ROUGH
      INTEGER :: A_ASSMGRPS
      INTEGER :: NUM_PVPR

      LOGICAL :: A_RECON
      LOGICAL :: H_OROG_GRAD
      LOGICAL :: ATMODS
      LOGICAL :: CMODS
      LOGICAL :: LMESO

      LOGICAL :: TRACER_A (0:A_MAX_TRVARS)
      LOGICAL :: TR_UKCA_A (0:A_MAX_UKCAVARS)
      LOGICAL :: AASSET   (AASSETS)
      INTEGER :: AOBINC   (MAX_AOBS)
      INTEGER :: AOBGRP   (MAX_AOBS)
      INTEGER :: AASPF    (AASSETS)
      INTEGER :: AASPL    (AASSETS)
      INTEGER :: RUN_TARGET_END( 6)

      COMMON/MODELA/ H_A_EWSPACE,H_A_NSSPACE,H_A_FIRSTLAT,H_A_FIRSTLONG,&
     &  H_A_POLELAT,H_A_POLELONG,A_ASSMGRPS,NUM_PVPR ,A_RECON,H_A_GROUP,&
     &  H_OROG_GRAD,ATMODS,CMODS,LMESO,TRACER_A,TR_UKCA_A,              &
     &  AASSET,AASPF,AASPL

!Total data length for primary fields for each submodel data partition
      INTEGER      LPRIM(N_SUBMODEL_PARTITION_MAX)
! Global (ie. dump on disk) version of LPRIM
      INTEGER      global_LPRIM(N_SUBMODEL_PARTITION_MAX)
!Total data length for primary fields for each internal model
      INTEGER      LPrimIM(N_INTERNAL_MODEL_MAX)
!Total data length for diagnostic flds for each submodel data partition
! Global (ie. dump on disk) version of LPrimIM
      INTEGER      global_LPrimIM(N_INTERNAL_MODEL_MAX)
      INTEGER      LDUMP(N_SUBMODEL_PARTITION_MAX)
! Global (ie. dump on disk) version of LDUMP
      INTEGER      global_LDUMP(N_SUBMODEL_PARTITION_MAX)
!Total data length for diagnostic flds for each internal model
      INTEGER      LDumpIM(N_INTERNAL_MODEL_MAX)
! Global (ie. dump on disk) version of LDumpIM
      INTEGER      global_LDumpIM(N_INTERNAL_MODEL_MAX)
!Total data length for secondary flds for each submodel data partition
      INTEGER      LSECD(N_SUBMODEL_PARTITION_MAX)
!Total data length for secondary flds for each internal model
      INTEGER      LSecdIM(N_INTERNAL_MODEL_MAX)
!Total workspace length for each submodel data partition
      INTEGER      LWORK(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each submodel data partition
      INTEGER      NHeadSub(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each internal model
      INTEGER      NHEAD(N_INTERNAL_MODEL_MAX)
!Total length of extra space for each submod. data part.
      INTEGER      LEXTRA(N_SUBMODEL_PARTITION_MAX)
!Data length for dual-time level ocean fields
      INTEGER      LPRIM_O2
      INTEGER      ITEM_MAX_REQ
      INTEGER      ITEM_MAX_ALL

      INTEGER      NRECS_S
      INTEGER      NTIMES_S
      INTEGER      NSERBLK_S
      INTEGER      NSERREC_S
      INTEGER      NLEVL_S
      INTEGER      NMAXLEV_S
      INTEGER      NPSLISTS_S
      INTEGER      NMAXPSL_S
      INTEGER      NHEAD_FILE(OUTFILE_S:OUTFILE_E)
      LOGICAL      LSTUSER

      COMMON/STRET/                                                     &
     &  LPRIM,LDUMP,LSECD,LWORK,NHEAD,LEXTRA,LPRIM_O2,LPrimIM,LDumpIM,  &
     &  LSecdIM,NHeadSub,ITEM_MAX_REQ,ITEM_MAX_ALL,NSERBLK_S,NSERREC_S, &
     &  NLEVL_S,NMAXLEV_S,NPSLISTS_S,NMAXPSL_S,LSTUSER,NRECS_S,NTIMES_S,&
     &  NHEAD_FILE,                                                     &
     &  global_LPRIM,global_LPrimIM,global_LDUMP,global_LDumpIM
      CHARACTER*1  H_ATMOS
      CHARACTER*1  H_OCEAN
      CHARACTER*1  H_SLAB
      CHARACTER*1  H_WAVE
      CHARACTER*1  H_FLOOR
      CHARACTER*1  H_STRAT
      CHARACTER*1  H_SLAB_CAL
      CHARACTER*1  H_TOTEM
      CHARACTER*1  H_GLOBAL(N_INTERNAL_MODEL_MAX         )
      INTEGER      H_VERS  (N_INTERNAL_MODEL_MAX,0:NSECTP)

      COMMON/CHOICE/ H_ATMOS,H_OCEAN,H_SLAB,H_WAVE,H_GLOBAL,H_SLAB_CAL, &
     &  H_TOTEM,H_FLOOR,H_STRAT

      COMMON/HVERS/ H_VERS

      REAL H_O_EWSPACE ,H_O_NSSPACE
      REAL H_O_FIRSTLAT,H_O_FIRSTLONG
      REAL H_O_POLELAT ,H_O_POLELONG

      INTEGER H_O_PTSPROW
      INTEGER N_COMP_O
      INTEGER H_NSIDEIMTO       ,H_NSIDEJMTO
      INTEGER SEAICE_TYPE       ,OCEAN_BASINS

      LOGICAL COX_Z,COX_Y,COX_P,COX_L,COX_PMSL
      LOGICAL COX_O,COX_X
      LOGICAL COX_1234
      LOGICAL COX_LCASE_I
      LOGICAL COX_LCASE_C,COX_OCARB
      LOGICAL TRACER_O(0:18)

      CHARACTER*1 O_ASSM_FIELDS(6)

      COMMON/MODELO/                                                    &
     &  H_O_EWSPACE,H_O_NSSPACE,H_O_FIRSTLAT,H_O_FIRSTLONG,H_O_POLELAT, &
     &  H_O_POLELONG,H_O_PTSPROW,N_COMP_O,H_NSIDEIMTO,H_NSIDEJMTO,      &
     &  SEAICE_TYPE,OCEAN_BASINS,COX_Z,COX_Y,COX_P,COX_L,COX_1234,      &
     &  COX_PMSL,COX_O,COX_X,COX_LCASE_I,COX_LCASE_C,COX_OCARB,TRACER_O,&
     &  O_ASSM_FIELDS

! These are set in SETMODL:
      INTEGER MEAN_NUMBER(N_INTERNAL_MODEL_MAX)
      COMMON/MODLMEAN/ MEAN_NUMBER

      REAL    H_W_EWSPACE ,H_W_NSSPACE
      REAL    H_W_FIRSTLAT,H_W_FIRSTLONG

      COMMON/MODELW/ H_W_EWSPACE ,H_W_NSSPACE,H_W_FIRSTLAT,H_W_FIRSTLONG

! Variables read in by namelist and used in SETMODL
      INTEGER      OCAAA   ,OCAAO   ,OCAAW
      INTEGER      NROWSO  ,NCOLSO  ,NLEVSO
      INTEGER      NROWSW  ,NCOLSW
      INTEGER      NWTRAIN
      REAL         EWSPACEA,NSSPACEA
      REAL         EWSPACEO,NSSPACEO
      REAL         EWSPACEW,NSSPACEW
      REAL         FRSTLATA,FRSTLONA
      REAL         FRSTLATO,FRSTLONO
      REAL         FRSTLATW,FRSTLONW

      LOGICAL      ZonAvOzone
      LOGICAL      ZonAvTppsOzone !! for tropopause-based ozone
      INTEGER      IVDF
      REAL         LATS
      REAL         LONS
      INTEGER      LWBND
      INTEGER      LWINC
      INTEGER      NECF(50)
      INTEGER      OASFLDID(7)

      INTEGER      OASLEV(6) ! dimensioned by max no of O-Assm groups
      INTEGER      OBAS
      INTEGER      OBS
      INTEGER      OCALB
      INTEGER      OCBOHaney
      INTEGER      OICE
      INTEGER      OIDYN
      INTEGER      OMP(4)
      REAL         POLELATA
      REAL         POLELONA
      REAL         POLELATO
      REAL         POLELONO
      INTEGER      PSA
      INTEGER      StLevGWdrag
      INTEGER      SWBND
      INTEGER      SWINC
      INTEGER      TCA(A_MAX_TRVARS)
      INTEGER      TC_UKCA(A_MAX_UKCAVARS)
      INTEGER      TCO(29)
      INTEGER      BotVDiffLev
      INTEGER      TopVDiffLev


      COMMON/STSHCOMM/                                                  &
     &  RUN_TARGET_END,                                                 &
     &  OCAAA,EWSPACEA,POLELATA,FRSTLATA,LATS,                          &
     &  NSSPACEA,POLELONA,FRSTLONA,LONS,                                &
     &  OCAAO,EWSPACEO,POLELATO,FRSTLATO,NCOLSO,NLEVSO,                 &
     &  NSSPACEO,POLELONO,FRSTLONO,NROWSO,                              &
     &  OCAAW,EWSPACEW,FRSTLATW,NCOLSW,                                 &
     &  NSSPACEW,FRSTLONW,NROWSW,NWTRAIN,                               &
     &  SWBND,LWBND,SWINC,LWINC,                                        &
     &  ZonAvOzone ,AOBINC   , ZonavTppsOzone,                          &
     &  StLevGWdrag,AOBGRP,                                             &
     &  BotVDiffLev,TopVDiffLev,OCALB,TCA,TC_UKCA,OIDYN,OBAS,OCBOHaney, &
     &  OBS,OICE,IVDF,PSA,NECF,OASLEV,TCO,OMP,OASFLDID

      CHARACTER(LEN=2) :: ATMOS_SR(0:NSECTP)
      CHARACTER(LEN=2) :: OCEAN_SR(0:NSECTP)
      CHARACTER(LEN=2) :: SLAB_SR (0:NSECTP)
      CHARACTER(LEN=2) :: WAVE_SR (0:NSECTP)
      CHARACTER(LEN=2) :: INDEP_SR(0:NSECTP)

      CHARACTER(LEN=1) :: BSPMSL
      CHARACTER(LEN=1) :: CCEW
      CHARACTER(LEN=1) :: FLOOR
      CHARACTER(LEN=1) :: IDO
      CHARACTER(LEN=1) :: LOSSM
      CHARACTER(LEN=1) :: MLMO
      CHARACTER(LEN=1) :: OAFLD(7)
      CHARACTER(LEN=1) :: OCARB
      CHARACTER(LEN=1) :: OROGR
      CHARACTER(LEN=1) :: OSFC
      CHARACTER(LEN=1) :: SCAL
      CHARACTER(LEN=1) :: SSTAnom
      CHARACTER(LEN=1) :: SWMCR
      CHARACTER(LEN=1) :: TOTAE
      CHARACTER(LEN=1) :: TOTEM
      CHARACTER(LEN=1) :: UPD175
      CHARACTER(LEN=1) :: MESO

      COMMON/STSHCHAR/                                                  &
     &  BSPMSL, CCEW, INDEP_SR, FLOOR, IDO, LOSSM, ATMOS_SR, MLMO,      &
     &  OAFLD, OCARB, OROGR, OSFC, OCEAN_SR, SCAL, SSTAnom, SWMCR,      &
     &  TOTAE, TOTEM, SLAB_SR, UPD175, MESO, WAVE_SR

      NAMELIST/STSHCOMP/                                                &
     & RUN_TARGET_END,                                                  &
     &  INDEP_SR    ,ATMOS_SR    ,OCEAN_SR ,SLAB_SR ,WAVE_SR,           &
     &  OCAAA       ,EWSPACEA    ,POLELATA ,FRSTLATA,LATS   ,           &
     &               NSSPACEA    ,POLELONA ,FRSTLONA,LONS   ,           &
     &  OCAAO       ,EWSPACEO    ,POLELATO ,FRSTLATO,NCOLSO ,NLEVSO,    &
     &               NSSPACEO    ,POLELONO ,FRSTLONO,NROWSO ,           &
     &  OCAAW       ,EWSPACEW    ,          FRSTLATW,NCOLSW ,           &
     &               NSSPACEW    ,          FRSTLONW,NROWSW ,           &
     &  SWBND       ,LWBND       ,SWINC    ,LWINC   ,OROGR  ,           &
     &  ZonAvOzone  , ZonAvTppsOzone,SWMCR       ,MESO     ,            &
     &  StLevGWdrag ,BotVDiffLev, TopVDiffLev,                          &
     &  OCALB       ,FLOOR       ,AOBINC   ,TOTAE   ,TOTEM  ,TCA,       &
     &  TC_UKCA     ,SSTAnom     ,SCAL     ,AOBGRP   ,                  &
     &  NECF        ,BSPMSL      ,CCEW              ,UPD175 ,           &
     &  OIDYN       ,OBAS        ,OCBOHaney,OBS     ,OICE   ,IVDF, IDO, &
     &  OCARB       ,MLMO        ,PSA      ,OSFC    ,                   &
     &  LOSSM       ,OASLEV      ,OAFLD    ,TCO     ,                   &
     &  OMP         ,OASFLDID
! MODEL end
!*L------------------ COMDECK LOOKADD ----------------------------------
!LL
!LL Purpose : Contains information about the format
!LL           of the PP header
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   4.0  12/09/95   Change NPERIODS to LBUSER3, BRSVD1 to BULEV,
!LL                   BRSVD2 to BHULEV and definitions for BRLEV and
!LL                   BHRLEV. Corresponding changes made to STWORK1A
!LL                   and PPHEAD1A. (Andrew Brady)
!LL  4.0  12/10/95  Change item 45 from lbuser7 to model_code. RTHBarnes
!LL  5.1  17/04/00    Fixed/Free format. P.Selwood.
!LL  5.2  25/09/00    Add LBCC_xxxx variables for the compressed
!LL                   LBC LOOKUP array                  P.Burton
!LL
!LL Programming standard :
!LL
!LL Logical components covered : F092
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LLEND -----------------------------------------------------------------
!

! Validity time
      Integer, Parameter :: LBYR   =1   ! Year
      Integer, Parameter :: LBMON  =2   ! Month
      Integer, Parameter :: LBDAT  =3   ! Day of month
      Integer, Parameter :: LBHR   =4   ! Hour
      Integer, Parameter :: LBMIN  =5   ! Minute
      Integer, Parameter :: LBDAY  =6   ! Day number

! Data time
      Integer, Parameter :: LBYRD  =7   ! Year
      Integer, Parameter :: LBMOND =8   ! Month
      Integer, Parameter :: LBDATD =9   ! Day of month
      Integer, Parameter :: LBHRD  =10  ! Hour
      Integer, Parameter :: LBMIND =11  ! Minute
      Integer, Parameter :: LBDAYD =12  ! Day number

      Integer, Parameter :: LBTIM  =13  ! Time indicator
      Integer, Parameter :: LBFT   =14  ! Forcast period (hours)
      Integer, Parameter :: LBLREC =15  ! Length of data record
      Integer, Parameter :: LBCODE =16  ! Grid type code
      Integer, Parameter :: LBHEM  =17  ! Hemisphere indicator
      Integer, Parameter :: LBROW  =18  ! Number of rows in grid
      Integer, Parameter :: LBNPT  =19  ! Number of points per row
      Integer, Parameter :: LBEXT  =20  ! Length of extra data
      Integer, Parameter :: LBPACK =21  ! Packing method indicator
      Integer, Parameter :: LBREL  =22  ! Header release number
      Integer, Parameter :: LBFC   =23  ! Field code
      Integer, Parameter :: LBCFC  =24  ! Second field code
      Integer, Parameter :: LBPROC =25  ! Processing code
      Integer, Parameter :: LBVC   =26  ! Vertical coordinate type
      Integer, Parameter :: LBRVC  =27  ! Coordinate type for reference
                                        ! level

      Integer, Parameter :: LBEXP  =28  ! Experiment number
      Integer, Parameter :: LBEGIN =29  ! Start record
      Integer, Parameter :: LBNREC =30  ! No of records-Direct access
                                        ! only
      Integer, Parameter :: LBPROJ =31  ! Met-O-8 projection number
      Integer, Parameter :: LBTYP  =32  ! Met-O-8 field type
      Integer, Parameter :: LBLEV  =33  ! Met-O-8 level code
      Integer, Parameter :: LBRSVD1=34  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD2=35  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD3=36  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD4=37  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBSRCE =38  ! =1111 to indicate following
                                        ! apply to UM
      Integer, Parameter :: DATA_TYPE =39  ! Indicator for real/int or
                                           ! timeseries
      Integer, Parameter :: NADDR  =40  ! Start address in DATA_REAL or
                                        ! DATA_INT
      Integer, Parameter :: LBUSER3=41  ! Free for user-defined function
      Integer, Parameter :: ITEM_CODE =42  !Stash code
      Integer, Parameter :: LBPLEV =43  ! Pseudo-level indicator (if
                                        ! defined)
      Integer, Parameter :: LBUSER6=44  ! Free for user-defined function
      Integer, Parameter :: MODEL_CODE =45 ! internal model identifier

      Integer, Parameter :: BULEV  =46  ! Upper level boundary
      Integer, Parameter :: BHULEV =47  ! Upper level boundary
      Integer, Parameter :: BRSVD3 =48  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BRSVD4 =49  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BDATUM =50  ! Datum value
      Integer, Parameter :: BACC   =51  ! (Packed fields) Packing
                                        ! accuracy
      Integer, Parameter :: BLEV   =52  ! Level
      Integer, Parameter :: BRLEV  =53  ! Lower level boundary
      Integer, Parameter :: BHLEV  =54  ! (Hybrid levels) A-level of
                                        ! value
      Integer, Parameter :: BHRLEV =55  ! Lower level boundary
      Integer, Parameter :: BPLAT  =56  ! Real latitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BPLON  =57  ! Real longitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BGOR   =58  ! Grid orientation
      Integer, Parameter :: BZY    =59  ! Zeroth latitude
      Integer, Parameter :: BDY    =60  ! Latitude interval
      Integer, Parameter :: BZX    =61  ! Zeroth longitude
      Integer, Parameter :: BDX    =62  ! Longitude interval
      Integer, Parameter :: BMDI   =63  ! Missing data indicator
      Integer, Parameter :: BMKS   =64  ! M,K,S scaling factor

      Integer, Parameter :: LBCC_LBYR   = 1  ! Year
      Integer, Parameter :: LBCC_LBMON  = 2  ! Month
      Integer, Parameter :: LBCC_LBDAT  = 3  ! Day of the month
      Integer, Parameter :: LBCC_LBHR   = 4  ! Hour
      Integer, Parameter :: LBCC_LBMIN  = 5  ! Minute
      Integer, Parameter :: LBCC_LBDAY  = 6  ! Day number
      Integer, Parameter :: LBCC_LBEGIN = 7  ! Start record
      Integer, Parameter :: LBCC_NADDR  = 8  ! Start address of DATA
! Mapping of MPP_LOOKUP; analogous to mapping in PP header

      Integer, Parameter :: P_NADDR=1    ! Address on local PE
      Integer, Parameter :: P_LBLREC=2   ! Local length of record

!*----------------------------------------------------------------------
! NADDR IS LOCATION IN PP-HEADER (LOOKUP) FOR START POSN OF VARIABLE
! ITEM_CODE is the location in PP header for a code defined as
!           (section number)*1000+item number
! DATA_TYPE is the location in the PP header defining data as REAL or
!           INTEGER.
! LBNPT is the location defining the number of points per row
!
!*L--------------------COMDECK  CANCILA ---------------------------
!
! Purpose : Contains index blocks for control of update of
!           ancillary fields.
!
! System component F0171
!
! History :
! Version   Date   Comment
! -------   ----   -------
!  3.4   23/06/94  Update comments. D. Robinson
!  3.4   04/10/94  Increase NANCIL_FIELDS from 43 to 71. RTHBarnes
!  4.1   22/05/96  Move NANCIL_FIELDS to comdeck CANCMAXA. D. Robinson.
!  4.4   28/07/97  Add LAMIPII to common block. R A Stratton
!  6.2   22/08/05  Remove un-needed ampersand. P.Selwood
!
! -------------------------------------------------------------------
!
! CANCMAXA Store maximum total no of atmosphere/slab ancillary fields.
!
! History :
! Version   Date   Comment
! -------   ----   -------
!  4.1   20/05/96  New comdeck. Increase from 71 to 81. D. Robinson.
!  4.3   18/3/97  Increase from 81 to 82.  William Ingram.
!  4.4   16/9/97  Increase from 82 to 89.  Richard Betts.
!  4.5   05/05/98 Increase from 89 to 109. D.M. Goddard
!  5.3   19/06/01 Increase from 109 to 110 Dave Tan
!  5.5   30/12/02 Increase from 110 to 122 Dave Robinson
!  6.1   07/04/04 Increase from 122 to 123 Andy Jones
!  6.1   08/11/04 Increase from 123 to 126 R.Sharp
!  7.3   31/05/10 Increase from 187 to 207 kdcorbin
!
! -------------------------------------------------------------------
! Type Declarations

      INTEGER NANCIL_FIELDS  ! No of Atmosphere & Slab Ancillary fields
      PARAMETER (NANCIL_FIELDS = 207)
! CANCMAXA end
! Type Declarations

      INTEGER                                                           &
     &  FILEANCIL,                                                      &
                         ! File number associated with ancillary fields
     &  NLOOKUP,                                                        &
                         ! Position of ancillary field in lookup tables.
     &  LOOKUP_STEP,                                                    &
                         ! Interval between PP Headers refering to
!                        ! to the same ancillary fields at diferent time
     &  LEVELS,                                                         &
                         ! Number of levels of data in each ancillary
!                        ! field.
     &  STASHANCIL,                                                     &
                         ! Stash codes for ancillary files
     &  D1_ANCILADD      ! Address of ancillary field in main data block


      COMMON/IXANCILA/ FILEANCIL(NANCIL_FIELDS),                        &
     &           NLOOKUP(NANCIL_FIELDS),                                &
     &           LOOKUP_STEP(NANCIL_FIELDS),                            &
     &           LEVELS(NANCIL_FIELDS),                                 &
     &           STASHANCIL(NANCIL_FIELDS),                             &
     &           D1_ANCILADD(NANCIL_FIELDS)

!*L---------- Control data calculated from NAMELIST-------------------
      LOGICAL                                                           &
     &         UPDATE                                                   &
     &,      L_SSTANOM                                                  &
                                ! Indicator if SST anom to be formed
                                ! (RECON=T) or used (-DEF,RECON)
     & ,     LAMIPII            ! True if special AMIP II updating

      INTEGER  FIELDCODE,                                               &
     &         STEPS

    !kdcorbin, 05/10
    INTEGER :: ANC_FILE_FINPUT
    CHARACTER*120 :: ANC_FILE_FNAME    
    
!*----------------------------------------------------------------------
      COMMON/CTANCILA/                                                  &
     &                L_SSTANOM,LAMIPII,                                &
     &         FIELDCODE(2,NANCIL_FIELDS),                              &
     &         STEPS(NANCIL_FIELDS),UPDATE(NANCIL_FIELDS),              &
               ANC_FILE_FINPUT(NANCIL_FIELDS),ANC_FILE_FNAME(NANCIL_FIELDS)
!kdcorbin, 05/10 - added anc_file_finput and anc_file_fname
!------------------------ nstypes.h ----------------------------------
!jhan:further renovation of ths file may be necessary params are dependent on dataset
!jhan: ALSO nstypes_cable.h should be unecessary nsoil/soil is only difference
      !--- Number of non-vegetation surface types
      Integer, Parameter :: NNVG  = 4

      !--- Number of plant functional types.
      Integer, Parameter :: NPFT  = 13
      
      !--- Number of surface types.
      Integer, Parameter :: NTYPE =17 
      
      !--- Index of the surface type 'Soil'
      !Integer, Parameter :: SOIL  = 16 
      !dhb599, 20110615: change made as per Peter Vohralik, item 1:
      Integer, Parameter :: SOIL  = 14

!--- Land surface types :
!--- original veg. tiles 
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!--- for testing these tiles are set = 1:5 
!     6 - Broadleaf Tree
!     7 - Needleleaf Tree
!     8 - C3 Grass
!     9 - C4 Grass
!    10 - Shrub
!--- for testing these tiles are set = 0
!    11 - 0 
!    11 - 0
!    11 - 0
!--- original non-veg tiles moved to these indices
!     14 - Urban
!     15 - Water
!     16 - Soil
!     17 - Ice


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
!*L --------------------- Comdeck: CENVIR   ----------------------------
!
!   Purpose: COMDECK defining Character enviroment variables used
!            by portable IO to open and close files
!
!   Author : R A Stratton      Date : 22/10/92
!
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL 3.2     28/05/93  Add file BAS_IND at unit number 58. M.Carter.
!LL
!LL 3.1     15/01/93  Increase no. of unit nos. from 1-99  to 1-199
!LL                   Dummy names have been set up temporarily for
!LL                   files 104-119. R.Rawlins
!LL
!LL 3.3     09/03/94  Separate data statements into COMDECK
!LL                   CENVIRDT. Also includes mods originally
!LL                   in RB221193 : Add source terms at unit no.110
!LL                   P.Burton and R.T.H Barnes
!LL

!   Vn3.0  12/02/93 - Environment variables PERTURB and TRANSP put in
!                     positions 37 and 97 respectively in character
!                     array FT_ENVIRON, and the appropriate character
!                     lengths put in LEN_FT_ENVIR. C. S. Douglas
!
! Type declarations
!
      CHARACTER*8 FT_ENVIRON(199)  ! Array holding enviroment variables
!                                  for filenames
      INTEGER     LEN_FT_ENVIR(199) ! character length of each variable
!


!
!Common Blocks for character and integer arrays
!
      COMMON/CENVIR/FT_ENVIRON
      COMMON/CLENVIR/LEN_FT_ENVIR
!
! ------------------------ Comdeck PARVARS -------------------------
! Parameters and common blocks required by the mpp-UM
! Added new comdeck AMAXSIZE required for new arrays in PARCOMM
! Add non-mpp option
!                                                      P.Burton
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
!====================== COMDECK AMAXSIZE ========================
! Description
!   This comdeck provides parameters giving the maximum likely sizes
!   of key UM resolution variables, useful for sizing static arrays.
!
!   History:
!   Model    Date     Modification history
!  version
!   4.2      18/11/96 New comdeck created.  P.Burton
!   4.3      24/01/97 Define MaxFieldSize to be a quarter of the
!                     SHMEM common block size.         P.Burton
!   4.4      3/7/97   Add MaxFieldSizeMes. Deborah Salmond
!   4.5     12/01/98  Added new variables, and changed sizes to
!                     correspond to global hi-res forecast - current
!                     largest configuration.                P.Burton
!                     Changed MAX_SHMEM_COMMON_SIZE to 3000000
!                     required for operational data assimilation.
!                                                           P.Burton
!   5.0     29/04/99  Changed variable names:
!                       P_ROWS_MAX -> ROWS_MAX
!                       P_LEVELS_MAX -> MODEL_LEVELS_MAX
!                       Q_LEVELS_MAX -> WET_LEVELS_MAX
!                       MaxHaloSize -> MaxHaloArea
!                     Removed variable:
!                       HALO_MAX (use PARPARM Max_Halo_Size instead)
!    5.0   29/04/99  Remove mpp #define
!    5.3   05/12/01  Remove MaxFieldSize, MaxFieldSizeMes and
!                    Max3DFieldSize.  S.Cusack
!    5.5   22/01/03  Increase ROW_LENGTH_MAX and HORIZ_DIM_MAX
!                    from 432 to 548. D Robinson.
!    6.1   31/08/04  Allow up to 100 levels.  R.Barnes
!    6.2   13/02/06  Increase max values of row_length and
!                    rows to cope with FOAM high res, as well
!                    as Global N320 and NAE.  M Martin.
!    6.2   24/11/05  Use max function for horiz_dim_max. R Barnes
!    6.2     11/01/06 Remove max_shmem_common_size here and
!                     in rdobs2.   Camilla Mathison/R Barnes
!

! Maximum sector size for I/O
      INTEGER,PARAMETER:: IO_SECTOR_SIZE_MAX=4096
      INTEGER,PARAMETER:: ROW_LENGTH_MAX   = 840 ! Maximum row length
      INTEGER,PARAMETER:: ROWS_MAX         = 600 ! Max no of rows

      ! MAX(ROW_LENGTH_MAX,ROWS_MAX)
      INTEGER,PARAMETER:: HORIZ_DIM_MAX=MAX(ROW_LENGTH_MAX,ROWS_MAX)

      INTEGER,PARAMETER:: MODEL_LEVELS_MAX = 100 ! Max no of total levels
      INTEGER,PARAMETER:: WET_LEVELS_MAX   = 100 ! Max no of wet levels
      INTEGER, PARAMETER :: Max2DFieldSize = ROW_LENGTH_MAX*ROWS_MAX +  &
     &  IO_SECTOR_SIZE_MAX
      INTEGER, PARAMETER :: MaxHaloArea    = HORIZ_DIM_MAX*Max_Halo_Size
!========================== COMDECK PARCOMM ====================
!
! *** NOTE : This comdeck requires comdeck PARPARM to be *CALLed
!            first.
!
!   Description:
!
!   This COMDECK contains COMMON blocks for the mpp-UM
!
!
!   Two COMMON blocks are defined:
!     i)  UM_PARVAR holds information required by the
!         Parallel Unified Model itself
!     ii) MP_PARVAR holds information required by the interface to
!         the Message Passing Software used by the PUM
!
!   Key concepts used in the inline documentation are:
!     o global data - the entire data domain processed by the UM
!     o LOCAL data - the fragment of the global data which is
!       stored by this particular process
!     o PERSONAL data - the fragment of the LOCAL data which is
!       updated by this particular process
!     o HALO data - a halo around the PERSONAL data which forms
!       the LOCAL data
!
!     Acronyms used:
!     LPG - Logical Process Grid, this is the grid of logical
!           processors; each logical processor handles one of the
!           decomposed parts of the global data. It does not
!           necessarily represent a physical grid of processors.
!
!   History:
!
!   4.1      27/1/96  New comdeck based on second section of
!                     old PARVARS.   P.Burton
!   4.2     19/08/96  Removed some unused variables, and added
!                     current_decomp_type variable to allow use
!                     of flexible decompositions.
!                     Added nproc_max to indicate the max. number
!                     of processors used for mpp-UM
!                                                      P.Burton
!   5.0     12/04/99  - Added halosize array to common block
!                     - Added halo_i and halo_j to common block
!                     - Added fld_type dimension to glsize
!                     - Added halo type dimension to lasize
!                     - Added fld_type dimension to lasize
!                     - Replaced blsizep/blsizeu by blsize with
!                       extra fld_type dimension
!                     - Replace attop etc. with at_extremity
!                     - Added g_pe_index to common block
!                                                      P.Burton
!   5.1     22/05/00  Removed DATA statement and put in BLKDATA
!                                                      P.Burton
!   5.1     26/01/00  - Renamed g_pe_index -> g_pe_index_EW
!                     - Added g_pe_index_NS
!                                                     P.Burton
!   5.2     02/08/00  Added g_at_extremity        P.Burton
!   5.3     14/09/01  Added sb_model_domain   P.Burton
!   5.5     06/08/00  Modification for parallelisation of WAM.
!                     Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!   5.5     30/01/03  Generalised datastart   P.Selwood.
!   5.5     07/02/03  SX now uses PARCOMM instead of SXCOMM    E.Leung
!   6.0     18/09/03  F90-fy continuation lines.               P.Dando
!   6.2     23/11/05  Removed all references to the wavemodel.
!                     T.Edwards
!
! -------------------- COMMON BLOCKS --------------------
!
! =======================================================
! Common block for the Parallel Unified Model
! =======================================================

      INTEGER :: first_comp_pe       ! top left pe in LPG
      INTEGER :: last_comp_pe        ! bottom right pe in LPG
      INTEGER :: current_decomp_type ! current decomposition type
      INTEGER :: Offx                ! standard halo size in East-West
      INTEGER :: Offy                ! standard halo size in North-South
      INTEGER :: halo_i              ! extended halo size in East-West
      INTEGER :: halo_j              ! extended halo size in North-South
      INTEGER :: halosize(Ndim_max,NHalo_max) ! available halo sizes
      INTEGER :: glsize(Ndim_max,Nfld_max) ! global data size
      INTEGER :: lasize(Ndim_max,Nfld_max,NHalo_max) ! local data size
      INTEGER :: blsize(Ndim_max,Nfld_max) ! personal data size

      ! position of personal data in global data (in terms of standard
      ! Fortran array notation
      INTEGER :: datastart(Ndim_max)

      ! Generalised version of datastart for *all* fieldtypes
      INTEGER :: datastart_f(Ndim_max,Nfld_max)

      INTEGER :: gridsize(Ndim_max)  ! size of the LPG in each dimension

      ! position of this process in the LPG 0,1,2,...,nproc_x-1 etc.
      INTEGER :: gridpos(Ndim_max)

      INTEGER :: sb_model_domain


      ! logicals indicating if a processor is at the edge of the LPG
      LOGICAL :: at_extremity(4)

      COMMON /UM_PARVAR/                                                &
     &  first_comp_pe, last_comp_pe, current_decomp_type, Offx, Offy,   &
     &  halo_i, halo_j, halosize, glsize, lasize, blsize, datastart,    &
     &  datastart_f, gridsize, gridpos                                  &
     &,                 at_extremity,sb_model_domain

      ! Common block for the Message Passing Software

      ! type of boundary (cyclic or static) in each direction
      INTEGER :: bound(Ndim_max)

      ! global copy of local data size
      INTEGER :: g_lasize(Ndim_max,Nfld_max,NHalo_max,0:maxproc)

      ! global copy of personal data area
      INTEGER :: g_blsize(Ndim_max,Nfld_max,0:maxproc)

      ! global copy of datastart
      INTEGER :: g_datastart(Ndim_max,0:maxproc)

      ! global copy of datastart_f
      INTEGER :: g_datastart_f(Ndim_max,Nfld_max,0:maxproc)

      INTEGER :: g_gridpos(Ndim_max,0:maxproc) ! global copy of gridpos

      ! Which processor column a given point is in: 0 -> nproc_x-1
      INTEGER :: g_pe_index_EW(1-Max_Halo_Size:                         &
     &  ROW_LENGTH_MAX+Max_Halo_Size)

      ! Which processor row a given point is in: 0 -> nproc_y-1
      INTEGER :: g_pe_index_NS(1-Max_Halo_Size:ROWS_MAX+Max_Halo_Size)

      INTEGER :: nproc      ! number of processors in current decomp
      INTEGER :: mype      ! number of this processor (starting from 0)
      INTEGER :: nproc_max  ! maximum number of processors
      INTEGER :: nproc_x    ! number of processors in x-direction
      INTEGER :: nproc_y    ! number of processors in y-direction

      ! array with the tids of the four neighbours in the horizontal
      ! plane
      INTEGER :: neighbour(4)

      INTEGER :: gc_proc_row_group  ! GID for procs along a proc row
      INTEGER :: gc_proc_col_group  ! GID for procs along a proc col
      INTEGER :: gc_all_proc_group  ! GID for all procs

      ! at_extremity for each processor
      LOGICAL :: g_at_extremity(4,0:maxproc)

      COMMON /MP_PARVAR/                                                &
     &  bound,                                                          &
     &  g_lasize,g_blsize,g_datastart, g_datastart_f, g_gridpos,        &
     &  g_pe_index_EW,g_pe_index_NS,                                    &
     &  nproc_max,nproc_x,nproc_y,                                      &
     &  neighbour,gc_proc_row_group,                                    &
     &  gc_proc_col_group, gc_all_proc_group                            &
     &  ,nproc,mype                                                     &
     &, g_at_extremity

! PARCOMM end
! CPRINTST defines print status for standard output messages

      ! Minimum output, only essential messages
      INTEGER,PARAMETER :: PrStatus_Min    = 1

      ! Normal informative messages + warnings
      INTEGER,PARAMETER :: PrStatus_Normal = 2

      ! Operational status, all informative messages
      INTEGER,PARAMETER :: PrStatus_Oper   = 3

      ! All informative + extra diagnostic messages
      INTEGER,PARAMETER :: PrStatus_Diag   = 4

      INTEGER PrintStatus ! Control volume of standard output messages
      COMMON/PrintSt/PrintStatus

! CPRINTST end


! Comdecks for ancillary files/fields.
! CANCFLDA List of Ancillary Fields - Atmosphere Stash Codes,
! Model Codes and Logical file Numbers
!
! -------------------------------------------------------------------
      INTEGER :: ITEM_CODES_ANCIL(NANCIL_FIELDS)  ! Stash Codes
      INTEGER :: MODEL_CODES_ANCIL(NANCIL_FIELDS) ! Model Codes
      INTEGER :: ANCIL_FILE_NO(NANCIL_FIELDS)     ! Logical file numbers

! -----------------------------------------
! Note 127:154 not used in UM ; set to zero
! -----------------------------------------

      DATA ITEM_CODES_ANCIL(1:100)/                                     &
     &  30,  33,  34,  35,  36,  37,  60,   0,  23,  20,                &
     &  40,  41,   0,  43,  44,   0,  46,  47,  50,  51,                &
     &  52,  53,  54,  55,  56,  26,  31,  24,  32,  28,                &
     &  29,  93,   0,   0,  48,   9,   0,   0,  58,  59,                &
     &  88,  87,  85,  57,  90,  17,  18, 301, 302, 303,                &
     & 304, 305, 306, 307, 308, 309, 310, 311, 312, 313,                &
     & 314, 315, 316, 317, 318, 319, 320, 127, 128, 129,                &
     &   0, 121, 122, 123, 124, 125, 126, 251, 207, 208,                &
     & 209, 160, 216, 217, 218, 213, 219, 220, 223, 321,                &
     & 322, 323, 324, 325, 326, 327, 328, 329, 330, 331/

    !kdcorbin, 05/10 - added 20 tracer flux variables
      DATA ITEM_CODES_ANCIL(101:NANCIL_FIELDS)/                         &
     & 332, 333, 334, 335, 336, 337, 338, 339, 340, 341,                &
     & 505, 418, 419, 420, 421, 422, 423, 424, 425, 426,                &
     & 130, 131, 132, 153, 151, 152,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   5,   6, 351, 352, 353, 354,                &
     & 355, 356, 357, 358, 359, 360, 361, 362, 363, 364,                &
     & 365, 366, 367, 368, 369, 370, 371, 480, 481, 482,                &
     & 483, 484, 485, 486, 487, 134, 135,3100,3101,3102,                &
     &3103,3104,3105,3106,3107,3108,3109,3110,3111,3112,                &
     &3113,3114,3115,3116,3117,3118,3119/

      DATA MODEL_CODES_ANCIL(1:100) /                                   &
     &   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,                &
     &   1,   1,   1,   0,   1,   0,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   0,   0,   1,   1,   0,   0,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1/

      DATA MODEL_CODES_ANCIL(101:NANCIL_FIELDS) /                       &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,                &
     &   1,   1,   1,   1,   1,   1,   1/
 
      DATA ANCIL_FILE_NO(1:100) /                                       &
     &   9,  10,  10,  10,  10,  10,   1,   0,   2,   3,                &
     &   4,   4,   0,   4,   4,   0,   4,   4,   5,   5,                &
     &   5,   5,   5,   5,   5,   5,   7,   6,   7,   8,                &
     &   8,   9,   0,   0,   4,   2,   0,   0,  12,  12,                &
     &  13,  13,  13,  14,  14,  10,  10,  15,  15,  15,                &
     &  15,  15,  15,  15,  15,  15,  15,  15,  15,  15,                &
     &  15,  15,  15,  15,  15,  15,  15,  12,  23,  23,                &
     &   0,  17,  18,  18,  18,  18,  12,  24,   4,   5,                &
     &   5,  19,  20,  21,  21,  21,  22,   4,   4,  16,                &
     &  16,  16,  16,  16,  16,  16,  16,  16,  16,  16/

    !kdcorbin, 05/10 - added 20 tracer flux variables
      DATA ANCIL_FILE_NO(101:NANCIL_FIELDS) /                           &
     &  16,  16,  16,  16,  16,  16,  16,  16,  16,  25,                &
     &  26,  27,  27,  27,  27,  27,  27,  27,  27,  27,                &
     &  28,  28,  29,  30,  31,  31,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,                &
     &   0,   0,   0,   0,  10,  10,  38,  39,  39,  39,                &
     &  40,  40,  41,  41,  42,  42,  42,  43,  43,  43,                &
     &  43,  43,  43,  44,  44,  44,  45,  46,  46,  46,                &
     &  46,  46,  46,  46,  46,  47,  47,  48,  49,  50,                &
     &  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,                &
     &  61,  62,  63,  64,  65,  66,  67/
! CANCFLDA end

!L External subroutines called:

      EXTERNAL                                                          &
     &        FILE_OPEN,                                                &
     &        READ_FLH, READHEAD, SETPOS

!L Namelist input

      NAMELIST/ANCILCTA/L_SSTANOM,LAMIPII

!     UPANCA Namelist
      INTEGER                                                           &
     &   ANC_REF_NO                                                     &
                          ! Ancil Ref. No : See comdeck CANCFLDA
     &  ,PERIOD                                                         &
                          ! Period of Updating Interval (Y/M/D/H)
     &  ,INTERVAL         ! Updating Interval

     !Added finput flag and filename to namelist info - kdcorbin, 04/10
      INTEGER :: FINPUT
      CHARACTER*Max_Filename_Len :: FNAME


      NAMELIST /UPANCA/ ANC_REF_NO,PERIOD,INTERVAL,FINPUT,FNAME

! Local Variables

      INTEGER                                                           &
     &        I,                                                        &
                               !
     &        ITEM,                                                     &
                               !
     &        J,                                                        &
                               !
     &        J1,                                                       &
                               !
     &        K,                                                        &
                               !
     &        LEN_IO,                                                   &
                               !
     &        LOOKUPS,                                                  &
                               !
     &        NFTIN,                                                    &
                               ! Current FTN number for ancillary data
     &        START_BLOCK                                               &
                               !
     &       ,STASH_CODE                                                &
                               ! Stash item code
     &       ,NREC_A,NREC_S                                             &
                               ! No of atmos & slab records
     &       ,STASH_ADDR                                                &
                               ! Stash address
     &       ,DUMMY                                                     &
                               !
     &       ,N_ANC_UPD        ! No of ancillaries to be updated
      DATA DUMMY /1/

      CHARACTER*8 CPERIOD      ! PERIOD in characters.
      LOGICAL                                                           &
     &        LFILE            !

!     SoilDepths : position of soil depths in level dependent constants
      Integer, Parameter :: SoilDepths = 4

      REAL P1,P2
      LOGICAL LNER
      LNER(P1,P2) = ((ABS(P1-P2))  >   (1.E-6*ABS(P1+P2)))

!L Internal Structure

      ICODE=0
      CMESSAGE=' '
      IOUNIT=0

!
!L  1.  Initialisation for atmosphere model

      DO I=1,NANCIL_FIELDS
        FILEANCIL(I) =ANCIL_FILE_NO(I)
        STASHANCIL(I)=ITEM_CODES_ANCIL(I)
      ENDDO
      FIELDCODE(:,:) = 0


! Set default values

      L_SSTANOM=.FALSE.
      LAMIPII=.FALSE.

!L  Read in control information from namelist

        REWIND 5
      READ(5,ANCILCTA)

!     Initialise FIELDCODE from Namelist UPANCA
      N_ANC_UPD = 0
      FINPUT = 0
      FNAME = ''
      DO I=1,NANCIL_FIELDS
        READ (5,UPANCA,ERR=101,END=101)
        FIELDCODE(1,ANC_REF_NO) = PERIOD
        FIELDCODE(2,ANC_REF_NO) = INTERVAL
        ANC_FILE_FINPUT(ANC_REF_NO) = FINPUT  !kdcorbin, 04/10
        ANC_FILE_FNAME(ANC_REF_NO) = FNAME  !kdcorbin, 04/10
        N_ANC_UPD = N_ANC_UPD+1
      ENDDO

 101  CONTINUE
      WRITE (6,*) ' '
      WRITE (6,*) N_ANC_UPD,' Atmos & Slab Ancillaries to be updated.'
      DO I=1,NANCIL_FIELDS
        IF (FIELDCODE(1,I) >  0) THEN
        IF (FIELDCODE(1,I) == 1) CPERIOD=' Years'
        IF (FIELDCODE(1,I) == 2) CPERIOD=' Months'
        IF (FIELDCODE(1,I) == 3) CPERIOD=' Days'
        IF (FIELDCODE(1,I) == 4) CPERIOD=' Hours'
        WRITE (6,*) 'Anc Ref No ',I,' Stash code ',ITEM_CODES_ANCIL(I), &
     &  ' Interval ',FIELDCODE(2,I),CPERIOD
        ENDIF
      ENDDO
      WRITE (6,*) ' '

! Check that ancillary field has valid address (>1) before proceding
!  to try and update it.  If not, switch off updating via FIELDCODE.
      DO I=1,NANCIL_FIELDS
      IF (STASHANCIL(I)  >   0) THEN
        if (model_codes_ancil(i) == slab_im) then
          stash_addr = si_slab(stashancil(i))
         else

          ! Including alternate sections from 0 to add 
          !   tracer fluxes - kdcorbin, 05/10
          !stash_addr = si_atmos(stashancil(i))
          isec=stashancil(i)/1000
          iind=stashancil(i)-isec*1000
          stash_addr = SI(iind,isec,1)

        endif
      ELSE
        stash_addr=0
      ENDIF
        IF (stash_addr  <=  1) THEN
          IF (FIELDCODE(1,I) >  0) THEN
           WRITE(6,*)' INANCILA: update requested for item ',i,         &
     &     ' STASHcode ',stashancil(i),' but prognostic address not set'
            WRITE(6,*)' FIELDCODE values reset to zeroes'
            FIELDCODE(1,I) = 0
            FIELDCODE(2,I) = 0
          END IF
        END IF
      END DO

!L  1.1 Set number of steps after which each ancillary field is updated
!       Zero is used for fields not to be updated

      DO I=1,NANCIL_FIELDS
        STEPS(I)=0
        IF (FIELDCODE(1,I) == 4)THEN
          STEPS(I)=FIELDCODE(2,I)*STEPS_PER_HR
        END IF
        IF (FIELDCODE(1,I) == 3) THEN
          STEPS(I)=FIELDCODE(2,I)*24*STEPS_PER_HR
        END IF

      IF (LCAL360) THEN
        IF (FIELDCODE(1,I) == 2) THEN
          STEPS(I)=FIELDCODE(2,I)*30*24*STEPS_PER_HR
        END IF
        IF (FIELDCODE(1,I) == 1) THEN
          STEPS(I)=FIELDCODE(2,I)*360*24*STEPS_PER_HR
        END IF
      ELSE
! Gregorian calender:
! If update interval is months or years, test each day. Further testing
! done in REPLANCA.

        IF (FIELDCODE(1,I) == 1.OR.FIELDCODE(1,I) == 2)THEN
         STEPS(I)=24*STEPS_PER_HR
        END IF
      END IF

      END DO

!L  1.2 Set master number of steps ANCILLARY_STEPS at which
!L      individual switches are tested.

!   Find first active field

      DO I=1,NANCIL_FIELDS
        IF (STEPS(I) >  0) THEN
          ANCILLARY_STEPS=STEPS(I)
          GOTO 121
        END IF
      END DO

! No above fields found

      ANCILLARY_STEPS=0

      GOTO 900
121   ITEM=I

!L      Set ANCILLARY_STEPS to lowest common denominater of
!L      frequencies for active fields

      DO I=ITEM+1,NANCIL_FIELDS
        IF (STEPS(I) <  ANCILLARY_STEPS                                 &
     &      .AND. STEPS(I) >  0) THEN
          IF (MOD(ANCILLARY_STEPS,STEPS(I)) == 0) THEN
            ANCILLARY_STEPS=STEPS(I)
          ELSE
            J1=STEPS(I)-1
            DO J=J1,1,-1
              IF ((MOD(ANCILLARY_STEPS,J) == 0).AND.                    &
     &           (MOD(STEPS(I),J) == 0)) THEN
                 GOTO 124
              ENDIF
            END DO
124         ANCILLARY_STEPS = J
          END IF
        END IF
      END DO

!L 1.2.4 Sea surface temperature must be updated when sea ice is update

      IF (STEPS(27) >  0.AND.STEPS(28) <= 0) THEN
         STEPS(28)=1
      END IF


!L 1.3 Set number of headers for each ancillary field

      DO I=1,NANCIL_FIELDS
        LEVELS(I)=1
!   Multilayer hydrology
        IF(I == 36)LEVELS(I)=SM_LEVELS
!   Multilayer aerosols
        IF(I >= 41.AND.I <= 43) LEVELS(I)=TR_LEVELS
!   Multilayer murk concentration and source
        IF(I >= 44.AND.I <= 45) LEVELS(I)=P_LEVELS
!   Multilayer user ancillaries
        IF(I >= 90.AND.I <= 109) LEVELS(I)=P_LEVELS
!   Multi-level ancillaries for sulphur cycle
        IF (I == 72) LEVELS(I) = P_LEVELS
        IF (I == 73) LEVELS(I) = P_LEVELS
        IF (I == 74) LEVELS(I) = P_LEVELS
        IF (I == 75) LEVELS(I) = P_LEVELS
        IF (I == 76) LEVELS(I) = P_LEVELS
        IF (I == 82) LEVELS(I) = NSULPAT
        IF (I == 83) LEVELS(I) = NTYPE
        IF (I == 84) LEVELS(I) = NPFT
        IF (I == 85) LEVELS(I) = NPFT
!   Multi-level ancillaries aerosol climatology
        IF(I >= 157.AND.I <= 177) LEVELS(I)=P_LEVELS
        IF(I >= 178.AND.I <= 185) LEVELS(I)=P_LEVELS
      END DO

      LEVELS(7)=OZONE_LEVELS
!! consider do a check for l_use_tpps_ozone and if set then
!! set levels(110)=tpps_ozone_levels
      LEVELS(10)=ST_LEVELS


!L 1.4 Read headers

      LOOKUPS=0

      DO I=1,NDATASETS

!  Initialise LOOKUP_START (=0 implies file I not required)
        LOOKUP_START(I)=0

!L Check whether each physical file is needed

        LFILE=.FALSE.
        DO 141 J=1,NANCIL_FIELDS


          IF (FILEANCIL(J) == I.AND.STEPS(J) >  0) THEN

            LFILE=.TRUE.

            !kdcorbin, 05/10 - added check for ancillary file name
            FINPUT=ANC_FILE_FINPUT(J)
            FNAME=ANC_FILE_FNAME(J)

          END IF
141     CONTINUE

        IF(LFILE) THEN

      ! Open the File
      ! Changed from opening with environmental variable to filename
      ! kdcorbin, 05/10
 
      NFTIN=FTNANCIL(I)

      if (ANC_FILE_FINPUT(NFTIN) == 1) Then
         AncFileName=ANC_FILE_FNAME(NFTIN)
      else
         CALL Fort_Get_Env(FT_ENVIRON(NFTIN),LEN_FT_ENVIR(NFTIN),  &
             AncFileName,Max_Filename_Len,icode)
      endif

        Write(6,*) ''
        Write(6,*) 'Field: ',NFTIN,FT_ENVIRON(NFTIN)
        Write(6,*) 'Opening Anc File: ',trim(AncFileName)

        len_anc_filename=len_trim( AncFileName)
! DEPENDS ON: file_open
        CALL FILE_OPEN(NFTIN,AncFileName,                      &
     &                 len_anc_filename,0,1,ICODE)

      !Original File Open:
      !! DEPENDS ON: file_open
      !  CALL FILE_OPEN(NFTIN,FT_ENVIRON(NFTIN),                         &
      !               LEN_FT_ENVIR(NFTIN),0,0,ICODE)

        IF(ICODE /= 0)THEN
          CMESSAGE='INANCLA: Error opening file'
          write(6,*) 'INANCILA: Error opening file on unit ',NFTIN,     &
     &               ' accessed from env.var.: ',FT_ENVIRON(NFTIN)
          RETURN
        ENDIF
! DEPENDS ON: setpos
        CALL SETPOS(NFTIN,0,ICODE)

!       Read in fixed header to get array dimensions
! DEPENDS ON: read_flh
        CALL READ_FLH(NFTIN,FIXHD(1,I),LEN_FIXHD,ICODE,CMESSAGE)
        IF (ICODE >  0) THEN
          WRITE (6,*) ' Error in reading fixed header for file ',I
          GO TO 9999   !  Return
        ENDIF

!       Check for negative dimensions
        IF (FIXHD(101,I) <= 0) FIXHD(101,I)=1
        IF (FIXHD(106,I) <= 0) FIXHD(106,I)=1
        IF (FIXHD(111,I) <= 0) FIXHD(111,I)=1
        IF (FIXHD(112,I) <= 0) FIXHD(112,I)=1
        IF (FIXHD(151,I) <= 0) FIXHD(151,I)=1
        IF (FIXHD(152,I) <= 0) FIXHD(152,I)=1
        IF (FIXHD(161,I) <= 0) FIXHD(161,I)=1

! Set start position of boundary fields for file
        LOOKUP_START(I)=LOOKUPS+1

        IF (LOOKUPS+FIXHD(152,I) >  NLOOKUPS) THEN
          WRITE (6,*) 'No room in LOOKUP table for Ancillary File ',I
          CMESSAGE='INANCILA: Insufficient space for LOOKUP headers'
          ICODE=14
          GO TO 9999   !  Return
        END IF

! DEPENDS ON: setpos
        CALL SETPOS(NFTIN,0,ICODE)
        IF (ICODE >  0) THEN
          WRITE (6,*) ' ERROR in SETPOS called from INANCA1A'
          WRITE (6,*) ' SETPOS attempted with Unit No ',NFTIN
          CMESSAGE = 'INANCA1A : ERROR in SETPOS'
          GO TO 9999    !   Return
        ENDIF

! DEPENDS ON: readhead
        CALL READHEAD(NFTIN,                                            &
     &                FIXHD(1,I),LEN_FIXHD,                             &
     &                INTHD(1,I),FIXHD(101,I),                          &
     &                REALHD(1,I),FIXHD(106,I),                         &
     &                LEVDEPC,FIXHD(111,I),FIXHD(112,I),                &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,DUMMY,                                &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                DUMMY,DUMMY,                                      &
     &                LOOKUP(1,LOOKUPS+1),FIXHD(151,I),FIXHD(152,I),    &
     &                FIXHD(161,I),                                     &
! COMDECK ARGPPX
! Description:
!   Argument comdeck for passing arrays which hold the
!   ppxref information.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.4  03/11/97  Removed MKPPXRF *DEF reference. K Rogers
! 4.4  04/11/97  Changed -RECON def line to allow for other small
!                execs which had used the RECON def. K Rogers
! 6.2 30/06/06 Rearrange statements to avoid too many continuation lines.
!                        J. M. Rodriguez
!
     & PPXI,PPXC,ppxRecs, PPXPTR,                                       &
! End of comdeck
     &                START_BLOCK,ICODE,CMESSAGE)

        IF (ICODE >  0) THEN
           WRITE(6,*) 'ERROR in READHEAD for Ancillary File ',I
           WRITE(6,*) 'Unit Number ',NFTIN
           GO TO 9999   !   Return
        ENDIF

!     Check calendar indicator
        IF ((     LCAL360 .and. FIXHD(8,I) /= 2) .or.                   &
     &      (.not.LCAL360 .and. FIXHD(8,I) /= 1) ) THEN
          ICODE=100+I
          CMESSAGE='INANCILA : Wrong calendar set in Ancillary File'
          WRITE (6,*) ' ******** Error in INANCILA ********'
          WRITE (6,*) ' Wrong calendar setting in Ancillary File ',I
          IF (LCAL360) THEN
            WRITE (6,*) ' Model run is set up for 360 day calendar.'
            WRITE (6,*) ' Ancillary File is for 365 day calendar.'
          ELSE
            WRITE (6,*) ' Model run is set up for 365 day calendar.'
            WRITE (6,*) ' Ancillary File is for 360 day calendar.'
          ENDIF
          WRITE (6,*) ' Rerun with correct ancillary file.'
          GO TO 9999   !  Return
        ENDIF

        FILE_LEVELS=1

        IF(I == 1) THEN
          FILE_LEVELS=OZONE_LEVELS
        ELSE IF(I == 2) THEN
          FILE_LEVELS=SM_LEVELS
! This is the maximum value that might be present on the ancillary
! file if it includes soil moisture in layers; otherwise only single
! level data is present and PR_FIXHD will not check value since
! FIXHD(110) will be zero
        ELSE IF(I == 3) THEN
            FILE_LEVELS=ST_LEVELS
        ELSE IF(I == 13) THEN   ! for multilevel aerosols
            FILE_LEVELS=TR_LEVELS
        ELSE IF(I == 14.or.I == 16) THEN   ! for murk and user ancil.
            FILE_LEVELS=P_LEVELS
        ELSE IF(I == 17.or.I == 18) THEN
!           multi-level sulphur cycle ancillary files.
            FILE_LEVELS=P_LEVELS
        else if (i == 25) then
! tropopause-based ozone file with tpps_ozone_levels
           file_levels=tpps_ozone_levels
        END IF


!L 1.4.2 Buffer in integer constants

           IF(FIXHD(100,I) >  0) THEN

! Check for error in file pointers

! Check validity of integer data and print out information
! All files except ozone should contain full fields

            IF(INTHD(6,I) /= ROW_LENGTH) THEN
! Ozone may contain zonal mean data
! also applies to tropopause-based ozone -- no 25.
              IF(.not.((I == 1) .or. (i  ==  25))                       &
     &        .OR. INTHD(6,I) /= 1) THEN
                ICODE=4
                CMESSAGE='INANCILA:integer header error'
                WRITE(6,*) ' INTHD(6) : ',INTHD(6,I),' ?'
                RETURN
              END IF
            END IF

            IF(INTHD(7,I) /= P_ROWS.AND.(I == 9.AND.INTHD               &
     &        (7,I) /= U_ROWS)) THEN
              ICODE=5
              CMESSAGE='INANCILA:integer header error'
              WRITE(6,*) ' INTHD(7) : ',INTHD(7,I),' ?'
              RETURN
            END IF

            IF (I == 1 .or. i  ==  25) THEN  ! Ozone or tpps-ozone file
              WRITE (6,*) ' '
              IF (INTHD(6,I) == 1)THEN
                WRITE (6,*) ' OZONE file contains zonal mean data for ',&
     &          INTHD(6,I),' points x ',INTHD(7,I),' rows'
              ELSEIF (INTHD(6,I) == ROW_LENGTH)THEN
                WRITE (6,*) ' OZONE file contains full fields for ',    &
     &          INTHD(6,I),' points x ',INTHD(7,I),' rows'
              ENDIF
! Check that correct ozone file has been provided.
              IF ((ZonAvOzone .and. i  ==  1) .or.                      &
     &          (ZonAvTppsOzone .and. i  ==  25)) THEN
!! Where is ZonAvOzone (and ZonAvTppsOzone) defined.
                IF (INTHD(6,I) /= 1) THEN
                  WRITE (6,*) ' Zonal Ozone Data is expected',          &
     &            ' for 1 point x ',P_ROWS,' rows'
                  ICODE = 51
                  CMESSAGE = 'INANCA1A : Wrong Ozone data provided.'
                  GO TO 9999   !  Return
                ENDIF
              ELSE
                IF (INTHD(6,I) /= ROW_LENGTH) THEN
                  WRITE (6,*) ' Ozone Data is expected for ',           &
     &            ROW_LENGTH,' points x ',P_ROWS,' rows.'
                  ICODE = 52
                  CMESSAGE = 'INANCA1A : Wrong Ozone data provided.'
                  GO TO 9999   !  Return
                ENDIF
              ENDIF
            ENDIF

          END IF

!L 1.4.3 Buffer in real constants

          IF(FIXHD(105,I) >  0) THEN

! Check validity of real header and print out information

           DO J=1,6
             IF(REALHD(J,I) >  (A_REALHD(J)+0.1).OR.                    &
     &         REALHD(J,I) <  (A_REALHD(J)-0.1))THEN
             IF(I /= 1.OR.(J /= 1.AND.J /= 4))THEN
               WRITE(6,*)(REALHD(K,I),K=1,6),(A_REALHD(K),K=1,6)
               ICODE=8
               CMESSAGE='INANCILA: REAL header Error.'
               RETURN
             END IF
             END IF
           END DO

         END IF

!L 1.4.4 Buffer in level dependent constants if required
!        Not retained in model after initial check

         IF(FIXHD(110,I) >  0) THEN

!L Only files 1 (Ozone), and 3 (Soil temperature)should contain multi
!L level data. File 2 (Soil moisture,snow depth,fractional snow time
!L and soil moisture in layers) may possibly also have multi level data.
!L FILES 13,14,16 (aerosols, murkiness, user ancil.) may also have
!L  multi level data.
!! Files 25 (TppsOzon) may also contain multi-level data
!! File 46 has multilevel data for cariolle ozone scheme

           If (I == 1   .or.                                            &
                                !  Ozone File
     &         I == 14  .or.                                            &
                                !  Murkiness File
     &         I == 16  .or.                                            &
                                !  User Ancillary File                  
     &         I == 46) Then    !  cariolle ozone File  

! Check that ancillary file is set up for correct vertical levels

            If (fixhd(111,I)-1 /= p_levels) Then
              icode=110
              write (CMESSAGE,*) ' Ancillary File set up for wrong',    &
     &        ' no of model levels. Anc ',fixhd(111,I)-1,               &
     &        ' Model ',p_levels
! DEPENDS ON: ereport
              Call Ereport ( RoutineName, icode, Cmessage )
            End If

            l_vert_mismatch = .false.

! Check eta_theta and eta_rho

            Do j=1,p_levels+1
              If (LNER( LEVDEPC(J), A_LEVDEPC(J,1) )) Then
                l_vert_mismatch = .true.
                exit
              End If
            End Do

            Do j=1,p_levels
              If (LNER( LEVDEPC(FIXHD(111,I)+J), A_LEVDEPC(J,2) )) Then
                l_vert_mismatch = .true.
                exit
              End If
            End Do

! Abort if there is a mis-match

            If (l_vert_mismatch) then
              write (6,*) 'Mismatch in vertical levels between model ', &
     &                    'and Ancillary File.'
              write (6,*) 'Anc File : ',title(i)
              write (6,*) 'Eta_Theta - Model'
              write (6,'(5F10.7)') (a_levdepc(k,1),k=1,p_levels+1)
              write (6,*) 'Eta_Theta - Anc File'
              write (6,'(5F10.7)') (levdepc(k),k=1,p_levels+1)
              write (6,*) 'Eta_Rho   - Model'
              write (6,'(5F10.7)') (a_levdepc(k,2),k=1,p_levels)
              write (6,*) 'Eta_Rho   - Anc File'
              write (6,'(5F10.7)') (levdepc(p_levels+1+k),k=1,p_levels)
                   ICODE=11
              Write (CMESSAGE,*) 'Mismatch in LEVDEPC ',                &
     &        'between model and Ancillary File.'
! DEPENDS ON: ereport
              Call Ereport ( RoutineName, icode, Cmessage )
            End If

           else if (i  ==  25) then !! tropopause-based ozone
             !! no checks to run....

           Else If (I == 2) Then  !  Soil Moisture File

             if (PrintStatus >= PrStatus_Diag .and. mype == 0 )then
               write (6,*)
               write (6,*) 'SoilDepths = ',SoilDepths
               write (6,*) 'SM_Levels  = ',sm_levels
               do j=1,sm_levels
                 write (6,*) 'model ',A_LEVDEPC(J,SoilDepths),          &
     &                       ' anc ',LEVDEPC(fixhd(111,I)*3+J)
               enddo
             endif

! Check Soil moisture levels

             Do J=1,SM_LEVELS
               If (LNER(LEVDEPC(fixhd(111,I)*3+J),                      &
     &                  A_LEVDEPC(J,SoilDepths))) Then
                 ICODE=12
                 CMESSAGE='INANCILA: error in LEVDEPC.'
                 RETURN
               End If
             End Do

           Else If (I == 3) Then  !  Deep Soil Temperature File

             If (PrintStatus >= PrStatus_Diag .and. mype == 0) Then
               write (6,*)
               write (6,*) 'SoilDepths = ',SoilDepths
               write (6,*) 'st_levels  = ',st_levels
               do j=1,st_levels
                 write (6,*) 'model ',A_LEVDEPC(J,SoilDepths),          &
     &                       ' anc ',LEVDEPC(fixhd(111,I)*3+J)
               End Do
             End If

! Check Deep Soil levels

             Do J=1,ST_LEVELS
               If (LNER(LEVDEPC(fixhd(111,I)*3+J),                      &
     &                  A_LEVDEPC(J,SoilDepths))) Then
                 ICODE=13
                 CMESSAGE='INANCILA: error in LEVDEPC.'
                 RETURN
               End If
             End Do

!L If aerosol file, check against model levels

           ELSE IF (I == 13) THEN

             DO J=1,TR_LEVELS
               DO J1=1,4
                 IF(LNER(LEVDEPC(J+(J1-1)*FIXHD(111,I)),A_LEVDEPC       &
     &                   (J,J1))) THEN
      WRITE(6,*)'Error in level dependent constants:Level=',J
                   WRITE(6,*)'Position=',J1
                   WRITE(6,*)'Value in model =',A_LEVDEPC(J,J1)
                   WRITE(6,*)'Value in ancillary data =',LEVDEPC(J+     &
     &                             (J1-1)*FIXHD(111,I))
                   ICODE=16
               CMESSAGE='INANCILA: error in LEVDEPC.'
                   RETURN
                 END IF
               END DO
             END DO

           END IF  !  If I

         END IF  !  If Fixhd(110,I) > 0

!L 1.4.5 Buffer in lookup table
! Set start position of boundary fields for file

         IF(FIXHD(150,I) >  0) THEN


           NREC_A = 0
           NREC_S = 0
           DO J = 1,FIXHD(152,I)
             IF (LOOKUP(MODEL_CODE,LOOKUPS+J)  ==  0 .or.               &
     &           LOOKUP(MODEL_CODE,LOOKUPS+J)  ==  imdi) THEN
               STASH_CODE = LOOKUP(ITEM_CODE,LOOKUPS+J)
               IF ((STASH_CODE >= 177 .and. STASH_CODE <= 179) .or.     &
     &             (STASH_CODE >= 210 .and. STASH_CODE <= 212)) THEN
                 LOOKUP(MODEL_CODE,LOOKUPS+J) = slab_im
                 NREC_S = NREC_S+1
               ELSE
                 LOOKUP(MODEL_CODE,LOOKUPS+J) = atmos_im
                 NREC_A = NREC_A+1
               END IF
             END IF
           END DO
           IF (NREC_A >  0) THEN
             WRITE (6,*) ' '
             WRITE (6,*) ' INANCA1A : submodel_id in ',NREC_A,          &
     &       ' records set to atmos_im in ancillary file ',I
           ENDIF
           IF (NREC_S >  0) THEN
             WRITE (6,*) ' '
             WRITE (6,*) ' INANCA1A : submodel_id in ',NREC_S,          &
     &       ' records set to slab_im in ancillary file ',I
           ENDIF

         END IF

         LOOKUPS=LOOKUPS+FIXHD(152,I)

       ELSE

!L  If file not required, zero fixed length header
         DO J=1,LEN_FIXHD
      FIXHD(J,I)=0
         END DO

         LOOKUP_START(I)=LOOKUPS+1
       END IF

      END DO

!L 1.5 Set positions in main data blocks

      DO I=1,NANCIL_FIELDS
        IF (STASHANCIL(I)  >   0) THEN
        IF (MODEL_CODES_ANCIL(I) == SLAB_IM) THEN
          D1_ANCILADD(I)=SI_SLAB(STASHANCIL(I))
        ELSE
          !Using si array directly to allow all sections - kdcorbin, 05/10
          !D1_ANCILADD(I)=SI_ATMOS(STASHANCIL(I))
          isec=stashancil(i)/1000
          iind=stashancil(i)-isec*1000
          D1_ANCILADD(I) = SI(iind,isec,1)
        ENDIF
        ELSE
          D1_ANCILADD(I)=0
        ENDIF
      ENDDO

!L 1.51 If a request is made to update a field, ensure that space for
!L     that field has been allocted in D1.

      DO I=1,NANCIL_FIELDS
        IF((FIELDCODE(1,I) >  0).AND.(D1_ANCILADD(I) <= 1)) THEN
          WRITE(6,*)' An address in D1 has not been set for ancillary   &
     & field number ',I
          ICODE=30
          CMESSAGE='INANCILA: updating for ancillary field is requested &
     & but no space has been allocated in D1'
          RETURN
        ENDIF
      END DO

!L 1.6 Set positions of data

      DO I=1,NANCIL_FIELDS
      NLOOKUP(I) =0
      LOOKUP_STEP(I)=0

! If LOOKUP_START=0 for file FILEANCIL(I), no fields required.
        IF   (STASHANCIL(I)  >   0) THEN
        IF (LOOKUP_START(FILEANCIL(I)) >  0) THEN

        DO J=LOOKUP_START(FILEANCIL(I)),LOOKUPS

          IF (LOOKUP(ITEM_CODE,J) == STASHANCIL(I)) THEN
            NLOOKUP(I)=J-LOOKUP_START(FILEANCIL(I))+1
            GOTO 161
          END IF

        END DO

! Find second occurence of data to set LOOKUP_STEP

161     LOOKUP_STEP(I)=0


        IF(J <  LOOKUPS) THEN

          DO J1=J+LEVELS(I),LOOKUPS
            IF (LOOKUP(ITEM_CODE,J1) == STASHANCIL(I)) THEN
              LOOKUP_STEP(I)=J1-NLOOKUP(I)-LOOKUP_START(FILEANCIL(I))+1
              GOTO 164
            END IF
          END DO
164      CONTINUE
        END IF

        END IF
        ENDIF

      END DO

!L SET LEVELS=2 FOR ICE FRACTION AND SNOW DEPTH, TO INDICATE PRESCENCE
!L fractional time fields

      LEVELS(27)=2
      LEVELS(9)=2

 900  CONTINUE
 9999 CONTINUE
      RETURN
      END SUBROUTINE INANCILA
