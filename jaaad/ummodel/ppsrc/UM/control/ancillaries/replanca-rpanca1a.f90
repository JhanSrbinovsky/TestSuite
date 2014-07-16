
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine REPLANCA ---------------------------------------------
!LL
!LL Purpose:  Updates ancillary fields as requested in FIELDCODE array.
!LL   Tests whether update is required for each field, allowing for
!LL   dependencies between fields. Uses LOOKUP array to find data for
!LL   appropriate time, reads a record and checks for current data
!LL   type. Reads second record if time interpolation required. Updates
!LL   the field. Under DEF RECON, the interface to the routine is
!LL   modified for use in the reconfiguration rather than the model.
!LL   Under DEF CAL360 the 360 day rather than the Gregorian calender
!LL   is used.
!LL
!LL Level 2 control routine for CRAY YMP
!LL
!LL C.Wilson    <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  22/02/93  Changes to allow updating of SLAB ref SST and ice
!LL                 ancillary fields (items 178,179) from SST/ice files.
!LL                 Correct bug if SST updated but not ice fraction.
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.2  15/04/93  Remove misleading warning messages if no time
!LL                  interpolation of SST
!LL   3.3  08/02/94  Modify calls to TIME2SEC/SEC2TIME to output/input
!LL                  elapsed times in days & secs, for portability. TCJ
!LL  3.3  22/11/93  Source term and aerosol ancillary fields added.
!LL   3.3  17/11/93 Initializing Integer UPDATE_MONTHS (N.Farnon)
!LL   3.3  08/12/93  Extra argument for READFLDS. D. Robinson
!LL   3.4  17/06/94  DEF CAL360 replaced by LOGICAL LCAL360
!LL                  Argument LCAL360 passed to SEC2TIME, TIME2SEC
!LL                                                   S.J.Swarbrick
!LL  3.4  20/07/94  Improve time interpolation by using reference time
!LL                  for ancillary updating.   R.T.H.Barnes.
!LL  3.4  05/09/94  Add murk and user ancillary fields.  R.T.H.Barnes.
!LL   3.4  18/05/94  Allow prognostic slabtemp under sea ice.J Thomson
!LL   3.4  31/8/94   More error trapping.        (William Ingram)
!LL  4.0  06/09/95  Only print time interpolation diagnostics when it
!LL                 is really done.  RTHBarnes.
!LL   4.0  08/09/95   Allow time interpolation of ozone fields and
!LL                   cater for zonal/full fields. D. Robinson
!LL   4.0  29/11/95   Set land points to zero for surface currents
!LL                   and Heat convergence fields. D. Robinson
!     4.1  18/06/96   Changes to cope with changes in STASH addressing
!                     Author D.M. Goddard.
!LL   4.1  22/05/96   Replace list of ancillary fields with call to
!LL                   comdeck CANCLSTA. D. Robinson.
!LL   4.2  08/11/96   Initialise fields to ensure haloes contain data
!LL                   for time interpolation for mpp runs. D. Robinson
!     4.1  16/12/96   Check ancillary files for non-constant polar rows,
!                     in reconfiguration only. Correct if LPOLARCHK=T
!                     Author D.M. Goddard
!LL   4.4  07/07/97   Alter SST and ice updating for AMIPII runs
!LL                   R A Stratton
!     4.4  13/11/97   Ancilary fields 72 - 89 no longer used for
!                     for user defined ancillaries. Code altered to
!                     ensure correct treatment of these fields.
!                     Author D.M. Goddard
!LL   4.4  25/07/97   (Reconfiguration only). Prevent failure when
!LL                   non-constant polar values for ancillary files are
!LL                   corrected.                         R. Rawlins
!LL   4.5  22/04/98   Add control of new NH3, soot aerosol emission
!LL                   ancillary fields. Plus minor message changes.
!LL                   R.Rawlins
!     4.5  22/10/98   Increase number of user ancillary fields by
!                     deleting existing four fields 68 - 72 and
!                     adding twenty to end 90 - 109
!                     Author D.M. Goddard
!LL   4.5  22/01/98   Correct level of second field read in for time
!LL                   interpolation.  D. Robinson
!LL   5.2  18/10/00   Account for U_FIELD and V_FIELD being different
!LL                   sizes when reading currents into D1 which are
!LL                   used, for example, by the slab model sea ice
!LL                   code.                               K.D.Williams
!LL   5.2  27/09/00   Correct setting of Artic ice depth for SN grid.
!LL   5.3  04/10/01   Removed land masking for ocn currents K.Williams
!     5.3     06/01   Add coastal tiling. Nic Gedney
!     5.4  27/08/02   Removed code which used slab sea ice fractional
!                     time ancil when slab SSTs were updated K.Williams
!LL   5.4  04/09/02   Allow time interpolation of effective vegetation
!LL                   parameter ancillaries.  M. Best.
!     5.5  05/02/03   Add control of biomass smoke emissions and mineral
!                     dust parent soil properties ancillary fields.
!                                                          P Davison
!     6.1  07/04/04   Add control for seawater DMS concentration.
!                                                            A. Jones
!     6.1  08/11/04   Add check for River Routing fields. R.Sharp
!     6.2  22/08/05   Fix operators and remove RECON. P.Selwood.
!     6.2  10/03/06   Allow for a single updating of soil moisture.
!                     Send level info back up. Clive Jones
!LL
!LL Programing standard : UMDP no 3, version no 2, dated 07/09/90
!LL
!LL Logical component covered : C71
!LL
!LL System task : C7
!LL
!LL   External Documentation: UMDP no C7
!LL
!LLEND-------------------------------------------------------------

       SUBROUTINE REPLANCA(I_YEAR,I_MONTH,I_DAY,I_HOUR,                 &
     &                     I_MINUTE,I_SECOND,I_DAY_NUMBER,              &
     &                     ANCIL_REFTIME,OFFSET_STEPS,                  &
     &                     P_FIELD,P_ROWS,U_FIELD,V_FIELD,D1,LAND,      &
     &                     A_STEP,LAND_FIELD,STEPS_PER_HR,              &
     &                     FLAND_CTILE,                                 &
     &                     TSTAR_LAND_CTILE,TSTAR_SEA_CTILE,            &
     &                     TSTAR_SICE_CTILE,                            &
     &                     ICE_FRACTION,TSTAR,TSTAR_ANOM,               &
     &                     SM_LEVELS,DZ_SOIL,SMC_UPDATED,               &
     &                     NS_SPACE,FIRST_LAT,                          &
     &                     LEN1_LOOKUP,LEN_FIXHD,LEN_INTHD,             &
     &                     LEN_REALHD,LEN_D1,FIXHD,INTHD,REALHD,        &
     &                     LOOKUP,RLOOKUP,FTNANCIL,LOOKUP_START,        &
     &                     NDATASETS,NLOOKUPS,                          &
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
     &                     ICODE,CMESSAGE,LCAL360)        ! Intent Out

      IMPLICIT NONE

      LOGICAL LCAL360

      INTEGER                                                           &
     &       I_YEAR,                                                    &
                                ! Curent Model Time
     &       I_MONTH,                                                   &
                                !   "      "     "
     &       I_DAY,                                                     &
                                !   "      "     "
     &       I_HOUR,                                                    &
                                !   "      "     "
     &       I_MINUTE,                                                  &
                                !   "      "     "
     &       I_SECOND,                                                  &
                                !   "      "     "
     &       I_DAY_NUMBER,                                              &
     &       ANCIL_REFTIME(6),                                          &
                                ! Reference time for ancillary updating
     &       OFFSET_STEPS,                                              &
                                ! Offset in timesteps of ref. from basis


     &       A_STEP,LAND_FIELD,STEPS_PER_HR,                            &


     &       P_FIELD,                                                   &
                                ! Size of horizontal fields
     &       P_ROWS,                                                    &
                                !
     &       U_FIELD,                                                   &
                                !   "  "      "         "
     &       V_FIELD,                                                   &
                                !   "  "      "         "
     &       NDATASETS,                                                 &
                                ! Number of ancillary datasets
     &       NLOOKUPS,                                                  &
                                ! Number of lookup tables
     &       LEN_D1             ! Size of primary data array

      INTEGER                                                           &
     &       LEN1_LOOKUP,                                               &
                                ! First dimension of lookup table
     &       LEN_FIXHD,                                                 &
                                ! Length of headers in data sets
     &       LEN_INTHD,                                                 &
                                !
     &       LEN_REALHD,                                                &
                         !
     &       FIXHD(LEN_FIXHD,NDATASETS),                                &
                                          ! Data set headers
     &       INTHD(LEN_INTHD,NDATASETS),                                &
                                          !
     &       LOOKUP(LEN1_LOOKUP,NLOOKUPS),                              &
                                          ! Data set lookup tables
     &       FTNANCIL(NDATASETS),                                       &
                                          ! FTN numbers of data sets
     &       LOOKUP_START(NDATASETS),                                   &
                                          ! Start of lookup tables
!                                         ! referring to each data set.
     &       SM_LEVELS          ! number of soil levels

      REAL                                                              &
     &       D1(LEN_D1),                                                &
                                !INOUT  Primary data array used to hold
!                               !       all fields except TSTAR and
!                               !       ICE_FRACTION
     &       ICE_FRACTION(P_FIELD),                                     &
                                    !INOUT  Ice frac of sea part of grid
!                                   !       box, updated if requested
     &       FLAND_CTILE(LAND_FIELD),                                   &
!                                   !IN  Fractional land on land pts.
     &       FLAND_G(P_FIELD),                                          &
                                 !WORK Frac land over all points.
     &       TSTAR(P_FIELD),                                            &
                                 !INOUT  TSTAR:updated if requested
     &       TSTAR_LAND_CTILE(P_FIELD),                                 &
!                                !INOUT  as above, but for land.
     &       TSTAR_SEA_CTILE(P_FIELD),                                  &
!                                !INOUT  as above, but for open sea.
     &       TSTAR_SICE_CTILE(P_FIELD),                                 &
!                                !INOUT  as above, but for sea-ice.
     &       TSTAR_ANOM(P_FIELD),                                       &
                                 !INOUT  SST anomaly,formed in recon;
                                 !       added if requested in model run
     &       REALHD(LEN_REALHD,NDATASETS),                              &
     &       RLOOKUP(LEN1_LOOKUP,NLOOKUPS)                              &
     &       ,NS_SPACE                                                  &
                               ! NS latitude spacing
     &       ,FIRST_LAT                                                 &
                               ! latitude of first gridpoint
     &       ,DZ_SOIL(SM_LEVELS) !OUT soil thicknesses

      LOGICAL                                                           &
     &       LAND(P_FIELD),                                             &
                                 ! WORK LAND mask
     &       SEA(P_FIELD),                                              &
                                 ! WORK SEA mask
     &       LTSTAR_SICE,                                               &
                                 ! IN TRUE if TSTAR_SICE has been read i
!                                ! from input dump.
!                                ! If FALSE set to TSTAR_SEA.
     &       SMC_UPDATED         ! OUT T if smc updated

      INTEGER                                                           &
     &       ICODE                                                      &
                       ! Return code
     &      ,IOUNIT       !OUT I/O unit passed out in RECON mode

      CHARACTER*(80)                                                    &
     &       CMESSAGE  ! Error message
!*
! Comdecks:------------------------------------------------------------
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
! ----------------------- header file: CNTLATM  -----------------------
! Description: Control variables for the Atmosphere model (read only).
!              Contains logical switches for science options needed
!              for addressing calculations and intermediate control.
!              Predominantly used for holding logical flags set up by
!              the UMUI - read in by a namelist, but can also hold
!              derived control flags set by the model.
!              [Note that CRUNTIMC holds accompanying run-time
!              constants.]
!
! Author : R. Rawlins
!
! History:
! Version  Date      Comment.
!  5.0 20/04/99  Extensive revision to earlier comdeck for conversion
!                to C-P 'C' dynamics grid. R. Rawlins
!  5.1  4/02/00  Restore energy correction switches removed at UM5.0
!                Also added additional switches for mass and moisture
!                R A Stratton.
!  5.1 13/04/00  IAU control moved to CTFilt. Adam Clayton
!  5.2 22/08/00  Reinstate Murk and Tracer switches. P.Selwood.
!  5.2 15/11/00  Reintroduce logicals for MOSES 2 and triffid. M.Best.
!  5.3   12/10/01   Remove hard-wired L_trivial_trigs.
!                   put problem_number in CNTLATM.    Terry Davies
!  5.3 27/07/01  Add logical switch L_MURK_RAD   S. Cusack
!  5.3    06/01  Introduce and hardwire logical for setting of
!                leads temperature.                         Nic Gedney
!  5.3 15/10/01  Added L_USE_METHOX switch. David Jackson
!  5.3 29/08/01  Sulphate and sea-salt logicals split into separate
!                versions for different processes.  A. Jones
!  5.3 09/04/01  Add logical switch for spatial degradation of
!                E-S radiation calculations.             S. Cusack
!
!  5.3 19/06/01   Stuff to handle tropopause-based ozone added
!                 -- see Radiation block              Dave Tan
!  5.4 15/08/02   Reconfiguration use changes. P.Selwood.
!  5.4  2/08/02  Add logical switch for PC2 cloud and condensation
!                scheme.                              Damian Wilson
!  5.4 24/10/02  Moved L_GWD and L_USE_USSP from CRUNTIMC. S. Webster
!  5.4 11/03/02  Remove comment lines from after #include
!                                                 S. Carroll
!  5.5 06/02/03  River routing support. P.Selwood.
!  5.5 05/02/03  Add logicals for biomass aerosol scheme     P Davison

!  5.5 17/02/03  Add two large-scale hydrology logicals.
!                L_TOP and L_PDM.                  Nic Gedney
!  5.5 08/01/03  Remove L_3DVAR, L_4DVAR, L_3DVAR_BG, LSINGLE_HYDROL,
!                L_H2_SULPH and L_LSPICE_BDY. D Robinson
!  5.5 21/01/03  Move L_USE_TPPS_OZONE to NLSTCATM namelist. D Robinson
!  5.5 13/03/03  Move I_OZONE_INT here from CRUNTIMC.
!                                                  Jean-Claude Thelen
!  5.5 03/02/03  Include L_mcr logicals.    R.M.Forbes
!  5.5 19/02/03  Remove redundant L_BL_LSPICE and L_RMBL   A.Lock
!  6.0 30/07/03  Include l_pc2_lbc for cloud frac. lbcs. Damian Wilson
!  6.0 19/08/03  Add runtime controls for 4 physics sections
!  6.1  02/08/04  Add logicals for stochem coupling. R Barnes
!  6.1 13/05/04  Add super_array_size variable                 A.Malcolm
!  6.1 07/04/04  Add logical for autoconversion de-biasing.   A. Jones
!  6.1 07/04/04  Add logicals for interactive DMS emissions.  A. Jones
!  6.2 25/01/06  Add iteration count logical test variable    T.Edwards
!  6.2 15/11/05  Add logical for aerosol optical depth     N. Bellouin
!  6.2 23/11/05  Add logical for RH and hygroscopic aerosols N.Bellouin
!  6.2 01/03/06  Add L_UPPER_STRATQ switch - David Jackson
!  6.2 06/01/06  Add logicals for seaice albedo options. J.Ridley
!  6.2 23/02/06  Add logicals for UKCA sub-model.  F.O'Connor
!  6.1 07/04/04  Add logical for RothC temperature function.  C.D. Jones
!  6.2 21/07/05  Add moisture_array_size variable              A.Malcolm
!  6.2 01/10/05  Include L_murk_lbc logical.    R.M.Forbes
!  6.2 07/11/05  Add L_MOD_BARKER_ALBEDO and L_USE_SPEC_SEA
!                to NLSTCATM.        James Manners
!  6.2 25/01/06  Add L_INHOM_CLOUD to NLSTCATM.    James Manners
!  6.2 09/03/06  Add logicals for inland basin rerouting. P.Falloon
!  6.2 24/02/06  Add logical for 10m windspeed pass A2O  J.Gunson

!  6.2 24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!  6.2 07/11/05 Add L_use_orog_corr to NLSTCATM.     James Manners
!  6.2 24/02/06  Add logical to allow DMS flux from ocean model J.Gunson
!  6.4 19/01/07 Removed a comment relating to A05_3c scheme. R A Stratton
!  6.4 08/01/06 Include Brooks cloud fraction logicals. Damian Wilson
!  6.4 10/01/07 Include mixing ratio control logical. Damian Wilson
!  6.6.2 10/06/09 Logicals for reading solar/volcanic forcing. Chris Jones
!------------------   General:  -------------------------------------
      INTEGER Model_domain        ! Domain of atmosphere model:
!                                   global,LAM,cyclic LAM,single column
! Model_domain meaningful names
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end

      INTEGER :: problem_number      ! type of problem to be solved

      INTEGER MAXSECTS            ! Max. no. of code sections
      PARAMETER (MAXSECTS=99)
      CHARACTER*3 H_SECT(0:MAXSECTS) ! Array of code section-versions

      ! Physics:   ------------------------------------
      LOGICAL :: l_ssice_albedo     ! Sea-ice albedo affected by snow

      LOGICAL :: l_sice_heatflux    ! Semi-impl update of sea-ice temp

      LOGICAL :: l_sice_meltponds ! Sea-ice albedo affected by meltponds

      LOGICAL :: l_sice_scattering  ! Sea-ice albedo affected scattering

      LOGICAL :: l_sice_hadgem1a ! HadGEM1 sea-ice albedo bug corrected

      LOGICAL :: L_NEG_TSTAR        ! Test for -ve surface temperature.
!
      ! Use sulphate aerosol in microphysics
      LOGICAL :: L_USE_SULPHATE_AUTOCONV

      ! Use sea-salt aerosol in microphysics
      LOGICAL :: L_USE_SEASALT_AUTOCONV

      ! Use soot aerosol in microphysics
      LOGICAL :: L_USE_SOOT_AUTOCONV

      ! Use biomass aerosol in microphysics
      LOGICAL :: L_USE_BMASS_AUTOCONV
      
      ! Use fossil-fuel organic carbon in microphysics
      LOGICAL :: L_USE_OCFF_AUTOCONV
      
      ! Use autoconversion de-biasing scheme in microphysics
      LOGICAL :: L_AUTO_DEBIAS
      ! Use sulphate aerosol no. in S-cycle
      LOGICAL :: L_USE_SULPHATE_SULPC

      ! Use sea-salt aerosol no. in S-cycle
      LOGICAL :: L_USE_SEASALT_SULPC

      ! Use soot aerosol no. in S-cycle
      LOGICAL :: L_USE_SOOT_SULPC

      ! Use biomass aerosol no. in S-cycle
      LOGICAL :: L_USE_BMASS_SULPC
      
      ! Use fossil-organic carbon aerosol no. in S-cycle
      LOGICAL :: L_USE_OCFF_SULPC
      
      ! Energy correction:
      LOGICAL :: L_emcorr    ! T: turns on energy correction code
      LOGICAL :: LMASS_corr  ! T: Apply mass correction
      LOGICAL :: LQT_corr    ! T: Apply total moisture correction
      LOGICAL :: LEMQ_print  ! T: Print additional info from em code
      LOGICAL :: LENERGY     ! T: if timestep to cal energy correction
      LOGICAL :: LFLUX_RESET ! T: if timestep to reset flux array in D1

      ! number of model timesteps per energy correction period.
      INTEGER :: A_ENERGYSTEPS

      ! Radiation:

      LOGICAL :: L_radiation     !  F: Turns off radiation code
      LOGICAL :: L_MICROPHY           !  Microphysics in sw rad scheme

!     Use mineral dust in radiation calculations
      LOGICAL :: L_USE_DUST

!     Use biogenic aerosol in radiation code
      LOGICAL :: L_USE_BIOGENIC

      ! Use SO4 aerosol from sulphur cycle for direct/indirect effect
      ! in radiation, the latter for both SW and LW.
      LOGICAL :: L_USE_SULPC_DIRECT
      LOGICAL :: L_USE_SULPC_INDIRECT_SW
      LOGICAL :: L_USE_SULPC_INDIRECT_LW
      ! Indirect radiative effect of sea-salt
      LOGICAL :: L_USE_SEASALT_INDIRECT

      ! Direct radiative effect of sea-salt
      LOGICAL :: L_USE_SEASALT_DIRECT
      LOGICAL :: L_USE_SOOT_DIRECT  ! direct radiative effects of soot
      LOGICAL :: L_USE_SOOT_INDIRECT  ! indirect effects of soot
      ! Use biomass aerosol for direct/indirect effect in radiation.
      LOGICAL :: L_USE_BMASS_DIRECT
      LOGICAL :: L_USE_BMASS_INDIRECT
      
      ! Use fossil-fuel organic carbon aerosol for direct/indirect
      ! effect in radiation
      LOGICAL :: L_USE_OCFF_DIRECT
      LOGICAL :: L_USE_OCFF_INDIRECT
      
!     Use aerosol climatologies in radiation instead of prognostic variables
!     Set on a species by species basis
      LOGICAL :: L_USE_ARCLBIOM   ! biomass burning aerosol
      LOGICAL :: L_USE_ARCLBLCK   ! black carbon
      LOGICAL :: L_USE_ARCLSSLT   ! sea salt
      LOGICAL :: L_USE_ARCLSULP   ! sulpahtes
      LOGICAL :: L_USE_ARCLDUST   ! mineral dust
      LOGICAL :: L_USE_ARCLOCFF   ! organic carbon (fossil fuel)
      LOGICAL :: L_USE_ARCLDLTA   ! delta aerosol

      ! Aerosol optical depth diagnostic was requested
      LOGICAL :: L_USE_AOD

      ! Clear-sky relative humidity is to be used instead of
      ! grid-box mean RH, for hygroscopic aerosols
      LOGICAL :: L_USE_CLEARRH

      ! controls the use of spatial degradation of radiation calc.
      LOGICAL :: L_rad_deg
      LOGICAL :: L_CTILE       ! Switch for coastal tiling.
!                              ! If TRUE then land and sea can
!                              ! coexist in the same gridbox.
!                              ! If FALSE the land fraction
!                              ! must be equal to 0 to 1
      LOGICAL :: L_TOP         ! If TRUE then TOPMODEL-based
!                              ! hydrology scheme.
      LOGICAL :: L_PDM         ! If TRUE then PDM-based
!                              ! hydrology scheme.
      LOGICAL :: L_SOIL_SAT_DOWN ! If TRUE then super-saturated soil
!                                ! moisture moves downward, else upward

      INTEGER :: H_SWBANDS    ! Number of shortwave radiation bands
      INTEGER :: H_LWBANDS    ! Number of longwave radiation bands
      INTEGER :: A_SW_RADSTEP ! Number of advection steps per SW step
      INTEGER :: A_LW_RADSTEP ! Number of advection steps per LW step
! Number of advection steps per prognostic/diagnostic SW and LW step.
!'Prognostic' and 'Diagnostic' refer to the frequency of the calls
! to radiation code in the Unified Model.
! In the case of time stepping prognostic and diagnostic refer to the
! slow and fast radiative timestep respectively. In the case of radiative
! forcing they refer to the prognostic and diagnostic calls to radiation.
      INTEGER :: A_SW_RADSTEP_DIAG
! Number of advection steps per 'fast' SW step (3C)
      INTEGER :: A_LW_RADSTEP_DIAG
! Number of advection steps per 'fast' LW step (3C)
      INTEGER :: A_SW_RADSTEP_PROG
! Number of advection steps per 'slow' LW step (3C)
      INTEGER :: A_LW_RADSTEP_PROG
! Number of advection steps per 'slow' LW step (3C)

      INTEGER :: i_ozone_int  ! Option for interpolation of ozone

      ! Cloud:

      LOGICAL:: L_CLD_AREA           ! controls cloud area parametriz.
      LOGICAL:: L_ACF_CUSACK         ! ... whether to have Cusack
      LOGICAL:: L_ACF_BROOKS         ! ... or Brooks
      LOGICAL:: L_PC2                ! controls PC2 cloud scheme
      LOGICAL:: L_PC2_RESET          ! run PC2 scheme diagnostically
      LOGICAL:: L_PC2_LBC            ! LBC's contain cloud fractions
      LOGICAL:: L_PC2_DIAG_SH        ! Use diagnostic convective shallow cloud
                                     ! in PC2.

      ! Assimilation:

      ! Switches for assm mode
      LOGICAL:: L_AC

      ! T: Use V_INT_TP to output Temp on model levels.
      LOGICAL:: L_VINT_TP

      ! UM6.5 - MODEL_ANALYSIS_HRS replaced by MODEL_ANALYSIS_MINS - 
      !         change A_ASSIM_START_HR and A_ASSIM_END_HR to 
      !         A_ASSIM_START_MIN, A_ASSIM_END_MIN
      !         so that all three variables  have the same flexibility
      ! Time at which data assimilation starts (Hours after Basis Time)
      INTEGER :: A_ASSIM_START_MIN
      ! Time at which data assimilation  ends
      INTEGER :: A_ASSIM_END_MIN

      ! Switch for assimilation mode
      CHARACTER(LEN=5) :: A_ASSIM_MODE
      
      !---  PMSL diagnostic  ---
      ! Orographic height threshold for new pmsl calculation
      ! from ATMOS_STASH_Misc in UMUI for Calc_NPMSL routine
      REAL :: NPMSL_HEIGHT

      ! Switch for interpolated winds in lbcs
      LOGICAL :: L_int_uvw_lbc  ! .true. for advecting winds
                                ! interpolated in boundary zone
      
      !---  Tracers ---
      ! Aerosol

      LOGICAL:: L_MURK          !      :Total aerosol field
      LOGICAL:: L_MURK_ADVECT   !      :Aerosol advection
      LOGICAL:: L_MURK_SOURCE   !Bndry :Aerosol source & sink terms
      LOGICAL:: L_MURK_BDRY     !Layer :UK Mes bndry model
      LOGICAL:: L_BL_TRACER_MIX !model :Bndry layer tracer mixing
      LOGICAL :: L_MURK_RAD
      LOGICAL :: L_murk_lbc    !  Murk aerosol lbcs active

      ! For Aero_Ctl (Sulphur cycle or Soot)

      INTEGER CALL_CHEM_FREQ     !No. times chem called per atm tstep

!     Mineral dust aerosol
      LOGICAL :: L_DUST
!     Use old version of dust_uplift scheme used in CAM NWP models
      LOGICAL :: L_CAM_DUST

      ! Sulphur cycle

      LOGICAL :: L_SULPC_SO2   ! S Cycle: SO2 MMR included
      LOGICAL :: L_SULPC_DMS   ! S Cycle: DMS MMR included
      LOGICAL :: L_SULPC_OZONE ! S Cycle: Ozone included for oxidation 
                               !          of DMS and SO2
      LOGICAL :: L_SULPC_SO2_O3_NONBUFFERED ! S Cycle: SO2+O3 reaction NOT
                                            ! buffered by NH3.
      LOGICAL :: L_SULPC_ONLINE_OXIDANTS ! Sulphur Cycle : Use online
                                         ! oxidants from UKCA
      LOGICAL :: L_SO2_SURFEM  ! SO2 Surface Emissions
      LOGICAL :: L_SO2_HILEM   ! SO2 High Level Emissions
      LOGICAL :: L_SO2_NATEM   ! SO2 Natural Emissions
      LOGICAL :: L_DMS_EM      ! DMS Emissions
      LOGICAL :: L_DMS_EM_INTER      ! Interactive DMS Emissions
      LOGICAL :: L_DMS_Ointer        ! DMS emissions from ocean model
      LOGICAL :: L_DMS_Liss_Merlivat ! Switches to determine which
      LOGICAL :: L_DMS_Wanninkhof    !    scheme to use for interactive
      LOGICAL :: L_DMS_Nightingale   !    sea-air exchange of DMS
      LOGICAL :: L_SULPC_NH3   ! S Cycle: NH3 tracer included
      LOGICAL :: L_NH3_EM      ! S Cycle: NH3 emiss included

      ! Soot cycle

      LOGICAL :: L_SOOT                ! Soot included
      LOGICAL :: L_SOOT_SUREM          ! surface Soot emiss included
      LOGICAL :: L_SOOT_HILEM          ! elevated Soot emiss included

      ! Biomass aerosol

      LOGICAL :: L_BIOMASS             ! Biomass aerosol included
      LOGICAL :: L_BMASS_SUREM         ! Sfc biomass emiss included
      LOGICAL :: L_BMASS_HILEM         ! Elevated bmass emiss included
      
      ! Fossil-fuel organic carbon aerosol
      
      LOGICAL :: L_OCFF                ! OCFF aerosol included
      LOGICAL :: L_OCFF_SUREM          ! Surface OCFF emissions included
      LOGICAL :: L_OCFF_HILEM          ! Elevated OCFF emiss included
      
      ! Carbon cycle

      ! Interactive 3D CO2 field for use with carbon cycle model
      LOGICAL :: L_CO2_INTERACTIVE
      ! Switch for Radiation Interaction with CO2 - kdcorbin, 06/10
      LOGICAL :: L_CO2_RADIATION   
      ! Switch for CABLE - kdcorbin, 03/10
      LOGICAL :: l_cable
      
      ! Switch for calculating CO2/tracer atmospheric mass - kdcorbin, 05/10
      LOGICAL :: L_TRACER_MASS, L_CO2_MASS
      INTEGER :: I_TRACERMASS_START
      ! Switch for running passive tracers using CO2 fluxes - rml, 1/7/13
      LOGICAL :: L_CO2_TRACER

      ! Switch for calculating methane atmospheric loss - kdcorbin, 05/10
      LOGICAL :: L_METHANE_LOSS
      INTEGER :: I_METHANE_TRACERS

      ! Switch for calculating mcf atmospheric/ocean loss - kdcorbin, 05/10
      LOGICAL :: L_MCF_LOSS
      INTEGER :: I_MCF_TRACERNUMBER

      ! Switch for calculating radon decay - kdcorbin, 05/10
      LOGICAL :: L_RADON_DECAY
      INTEGER :: I_RADON_TRACERNUMBER

      ! CO2 Mass - kdcorbin, 05/10
      REAL :: CO2EMITMASS

     ! Tracer Mass - kdcorbin, 05/10
      REAL :: TMASS(21)

      LOGICAL :: L_CO2_EMITS          ! Include surface emissions

      ! 10m windspeed for air/sea gas flux calculations
      LOGICAL :: L_PSSWND          ! pass 10m windspeed to the ocean

      ! Dust deposition for ocean biology
      LOGICAL :: L_DUST2OCN        ! pass dust dep fields to the ocean

      LOGICAL :: L_Q10                  ! control T fn for soil resp

      ! Switch for turning off boundary layer code

      LOGICAL :: L_bl            !  F: Turns off boundary layer code
      ! MOSES II and Triffid logicals--------------------

      LOGICAL :: L_VEG_FRACS          ! Switch for vegetation fractions
      LOGICAL :: L_SNOW_ALBEDO        ! Prognostic snow albedo
      LOGICAL :: L_TRIFFID            ! Switch for interactive veg model
      LOGICAL :: L_PHENOL             ! Switch for leaf phenology

      ! Switch for running TRIFFID in equilibrium mode
      LOGICAL :: L_TRIF_EQ

      ! Switch for starting NRUN mid-way through a TRIFFID calling
      ! period
      LOGICAL :: L_NRUN_MID_TRIF
      LOGICAL :: L_DISTURB      ! Switch for anthropogenic disturbance

      INTEGER :: CAN_MODEL ! Switch for thermal vegetation canopy

      ! Vegetation:

      ! Update frequency for leaf phenology (days)
      INTEGER :: PHENOL_PERIOD

      INTEGER :: TRIFFID_PERIOD ! Update frequency for TRIFFID (days)

      ! Switch for anthropogenic heat source 
      LOGICAL :: l_anthrop_heat_src 

      ! Hardwire until re-assessment of whether these need to be
      ! re-introduced as UMUI-set switches.
! RR old switches needed for addressing but should be defunct?
! RR - can be set with parameters if needed in the interim. ->

      ! Large scale precipitation:

      LOGICAL :: L_rain          !  F: Turns off precipitation code

      ! 'New' cloud/precip microphysics, Defunct, only mixed phase
      ! phys supported
      LOGICAL, PARAMETER :: L_LSPICE    =.true.

      ! Microphysics complexity
      LOGICAL :: L_mcr_qcf2    !  Include second ice variable
      LOGICAL :: L_mcr_qrain   !  Include prognostic rain
      LOGICAL :: L_mcr_qgraup  !  Include prognosic graupel
      LOGICAL :: L_mcr_qcf2_lbc    !  Second ice variable lbcs active
      LOGICAL :: L_mcr_qrain_lbc   !  Prognostic rain lbcs active
      LOGICAL :: L_mcr_qgraup_lbc  !  Prognosic graupel lbcs active

      ! Controls the use of new RHcrit parametrization, option in Sec 9
      ! vn 2A.
      LOGICAL :: L_RHCPT

! Logicals for different radiation packages
      LOGICAL :: L_FORCING
! Calculate radiative forcings (3C)
      LOGICAL :: L_RADIANCE
! Calculate radiances          (3C)
      LOGICAL :: L_TIMESTEP
! Use new timestepping scheme  (3C)
      LOGICAL :: L_WENYI       
! Include Wenyi's pressure & temperature scaling (3A/3C)

      ! Convection:

      LOGICAL :: L_3D_CCA             ! Use 3D conv cloud amount
      LOGICAL :: L_PHASE_LIM          ! Limits phase change of precip
                                      ! in convective downdraught
      LOGICAL :: L_CCRAD              ! Main logical, will remove code
                                      ! connected with CCRad
                                      ! (including bugfixes)
      LOGICAL :: L_3D_CCW             ! Requires l_ccrad=.TRUE.
                                      ! .TRUE. : Radiation to use 3d ccw
                                      ! profile passed to it from
                                      ! convection.
                                      ! .FALSE.: Radiation constructs
                                      ! mean CCW profile from cclwp,ccb
                                      ! and cct as in original.

      ! Timestep frequency for calling convection
      ! Hardwired to calling every timestep
      INTEGER,PARAMETER :: A_CONV_STEP = 1

      ! GWD scheme:
      LOGICAL :: L_GWD        ! Use SSO drag scheme
      LOGICAL :: L_USE_USSP   ! Use spectral GWD scheme

      ! Radiation:

      ! Changes to open sea albedo for HadGEM1
      LOGICAL :: L_MOD_BARKER_ALBEDO ! Modified Barker albedo
      LOGICAL :: L_USE_SPEC_SEA      ! Spectr. dep. sea albedos

      ! Use modulus of fluxes to remove negative effective extinctions
      LOGICAL :: L_MOD_K_FLUX

      ! Fix the selection of fractional sea points in LW radiation
      LOGICAL :: L_CTILE_FIX

      ! Fix instability in quadratic correction to LW source term
      LOGICAL :: L_QUAD_SRC_FIX

      ! Scale the condensed water content to simulate
      ! inhomogeneous clouds
      LOGICAL :: L_INHOM_CLOUD

! Orography correction to SW radiation
      LOGICAL :: L_use_orog_corr    !  Find gradients from mean orog
      LOGICAL :: L_use_grad_corr    !  Use ancillary X & Y gradients

      ! Tropopause-based Ozone Scheme
      LOGICAL :: L_use_tpps_ozone   !  Use TPPS ozone scheme

! Methane oxidation
      REAL    :: Z_TOP
      LOGICAL :: L_USE_METHOX

! STOCHEM coupling to radiation
      LOGICAL :: L_USE_STOCHEM_CH4   ! for methane
      LOGICAL :: L_USE_STOCHEM_O3    ! for ozone
      ! River Routing
      LOGICAL :: L_RIVERS
      LOGICAL :: L_INLAND   ! control rerouting of inland basin water
      REAL    :: RIVER_STEP

      ! Hydrology:

      LOGICAL :: L_hydrology     !  F: Turns off hydrology code

! Max humidity in STRATQ
      LOGICAL :: L_UPPER_STRATQ

      LOGICAL, PARAMETER :: LMOSES        =.TRUE.  ! MOSES hydrology
      LOGICAL :: L_ICOUNT       !  T: Output iteration counts

      ! Mixing ratios:

      Logical :: l_mr_physics1            ! Use mixing ratio in
                                          ! atmos_physics1
      Logical :: l_mr_physics2            ! Use mixing ratio in
                                          ! atmos_physics2
! Stochastic Physics Random Parameters      
      LOGICAL :: L_RPSEED_READ  !  T: Read in previously specified seed
      LOGICAL :: L_RPSEED_WRITE !  T: WRITE out seed


! Ozone tracer as input to radiation scheme      
      LOGICAL :: L_USE_CARIOLLE
      LOGICAL :: L_USE_OZONEINRAD

! RR old switches---------------------------------------- <-

! OASIS coupling
      LOGICAL :: L_OASIS   ! OASIS coupling switch
      LOGICAL :: L_COUPLE_MASTER    ! Couple through master PE
      INTEGER :: OASIS_COUPLE_FREQ  ! Coupling frequency in
                                    ! number of timesteps. 

!     Logicals for UK Chemistry and Aerosols (UKCA) Model

      LOGICAL :: L_ukca           ! True when UKCA is switched on

! Natural climate forcing
      LOGICAL :: L_SCVARY            ! time varying solar forcing
      LOGICAL :: L_VOLCTS            ! time varying volcanic forcing

      COMMON/ CNTLCATM/                                                 &
     &  Model_domain,L_emcorr,                                          &
     &  L_OASIS, OASIS_COUPLE_FREQ, L_COUPLE_MASTER,                    &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                          &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,              &
     &  problem_number,                                                 &
     &  L_SNOW_ALBEDO,l_ssice_albedo,L_MICROPHY,H_SWBANDS,H_LWBANDS,    &
     &  A_SW_RADSTEP,A_LW_RADSTEP,                                      &
     &  A_SW_RADSTEP_DIAG,A_LW_RADSTEP_DIAG,                            &
     &  A_SW_RADSTEP_PROG,A_LW_RADSTEP_PROG,                            &
     &  L_CLD_AREA,L_ACF_CUSACK,L_ACF_BROOKS,L_PC2,L_PC2_RESET,         &
     &  L_PC2_LBC,L_PC2_diag_sh,L_rad_deg, L_AC ,L_VINT_TP,             &
     &  A_ASSIM_START_MIN, A_ASSIM_END_MIN,                             &
     &  NPMSL_HEIGHT,                                                   &
     &  LENERGY, LFLUX_RESET, LMASS_CORR , LQT_CORR, LEMQ_PRINT,        &
     &  A_ENERGYSTEPS,L_GWD,L_USE_USSP,                                 &
     &  L_Murk, L_MURK_ADVECT, L_MURK_SOURCE, L_MURK_BDRY,              &
     &  L_MURK_RAD, L_murk_lbc, L_BL_TRACER_MIX, L_int_uvw_lbc,         &
     &  L_DUST,L_CAM_DUST,L_SULPC_SO2,L_SULPC_DMS,L_SULPC_OZONE,        &
     &  L_SULPC_SO2_O3_NONBUFFERED,L_SO2_SURFEM, L_SO2_HILEM,           &
     &  L_SO2_NATEM,L_DMS_EM, L_DMS_EM_INTER,                           &
     &  L_SULPC_ONLINE_OXIDANTS,                                        &
     &  L_DMS_Ointer,                                                   &
     &  L_DMS_Liss_Merlivat, L_DMS_Wanninkhof, L_DMS_Nightingale,       &
     &  L_SULPC_NH3, L_NH3_EM, L_SOOT, L_SOOT_SUREM, L_SOOT_HILEM,      &
     &  L_BIOMASS, L_BMASS_SUREM, L_BMASS_HILEM,                        &
     &  L_PSSWND, L_DUST2OCN,                                           &
     &  L_RHCPT, L_CCRAD, L_3D_CCA, L_3D_CCW, L_PHASE_LIM,              &
     &  L_CO2_INTERACTIVE,L_CO2_EMITS,l_cable,                          &
! rml 1/7/13
     &  L_CO2_TRACER,                                                   &
!kdcorbin, 08/10
     &  L_CO2_RADIATION,L_TRACER_MASS,L_CO2_MASS,I_TRACERMASS_START,    &
     &  L_METHANE_LOSS,I_METHANE_TRACERS,L_MCF_LOSS,I_MCF_TRACERNUMBER, &
     &  L_RADON_DECAY,I_RADON_TRACERNUMBER,                             &
     &  L_Q10, L_NEG_TSTAR, L_VEG_FRACS, L_TRIFFID, L_PHENOL,           &
     &  L_TRIF_EQ, L_NRUN_MID_TRIF, L_DISTURB,                          &
     &  CAN_MODEL, PHENOL_PERIOD, TRIFFID_PERIOD,                       &
     &  L_USE_SEASALT_INDIRECT, L_USE_BIOGENIC,                         &
     &  L_USE_SEASALT_DIRECT, L_USE_DUST, L_USE_SULPC_INDIRECT_SW,      &
     &  L_USE_SULPC_INDIRECT_LW, L_USE_SULPC_DIRECT,                    &
     &  L_USE_SULPHATE_AUTOCONV, L_USE_SEASALT_AUTOCONV, L_AUTO_DEBIAS, &
     &  L_USE_SULPHATE_SULPC, L_USE_SEASALT_SULPC,                      &
     &  L_OCFF, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_OCFF_AUTOCONV,        &
     &  L_USE_OCFF_SULPC, L_USE_OCFF_DIRECT, L_USE_OCFF_INDIRECT,       &
     &  L_USE_STOCHEM_CH4, L_USE_STOCHEM_O3,                            &
     &  L_MOD_BARKER_ALBEDO, L_USE_SPEC_SEA, L_MOD_K_FLUX, L_CTILE_FIX, &
     &  L_QUAD_SRC_FIX, L_USE_TPPS_OZONE, I_OZONE_INT,                  &
     &  L_ICOUNT,                                                       &
     &  L_use_orog_corr, L_use_grad_corr,                               &
     &  CALL_CHEM_FREQ, L_USE_SOOT_DIRECT, L_USE_SOOT_INDIRECT,         &
     &  L_USE_SOOT_AUTOCONV, L_USE_SOOT_SULPC, L_USE_BMASS_DIRECT,      &
     &  L_USE_BMASS_INDIRECT, L_USE_BMASS_AUTOCONV, L_USE_BMASS_SULPC,  &
     &  L_USE_ARCLBIOM, L_USE_ARCLBLCK,  L_USE_ARCLSSLT,                &
     &  L_USE_ARCLSULP, L_USE_ARCLDUST,  L_USE_ARCLOCFF, L_USE_ARCLDLTA,&
     &  L_ukca,                                                         &
     &  L_CTILE, L_RIVERS,L_INLAND, RIVER_STEP,                         &
     &  l_sice_meltponds, l_sice_scattering, l_sice_hadgem1a,           &
     &  L_USE_METHOX,Z_TOP,l_mr_physics1,l_mr_physics2,                 &
     &  L_INHOM_CLOUD,                                                  &
     &  L_radiation,L_FORCING,L_TIMESTEP,                               &
     &  L_RADIANCE, L_WENYI, L_bl, L_rain, L_hydrology,                 &
     &  L_TOP,L_PDM,L_USE_AOD,L_USE_CLEARRH,L_UPPER_STRATQ,             &
     &  L_sice_heatflux,L_SOIL_SAT_DOWN,                                &
     &  L_SCVARY,L_VOLCTS,                                              &
     &  l_anthrop_heat_src,                                             &
     &  L_USE_CARIOLLE,L_USE_OZONEINRAD,                                & 
     &  L_RPSEED_READ, L_RPSEED_WRITE,                                  &
     ! Character variables need to be at the end.
     &  H_SECT, A_ASSIM_MODE

      NAMELIST/NLSTCATM/                                                &
     &  Model_domain,L_emcorr,                                          &
     &  L_OASIS, OASIS_COUPLE_FREQ, L_COUPLE_MASTER,                    &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                          &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,              &
     &  problem_number,                                                 &
     &  L_SNOW_ALBEDO,l_ssice_albedo,L_MICROPHY,H_SWBANDS,H_LWBANDS,    &
     &  A_SW_RADSTEP,A_LW_RADSTEP,                                      &
     &  A_SW_RADSTEP_DIAG,A_LW_RADSTEP_DIAG,                            &
     &  A_SW_RADSTEP_PROG,A_LW_RADSTEP_PROG,                            &
     &  L_CLD_AREA,L_ACF_CUSACK,L_ACF_BROOKS,L_PC2,L_PC2_RESET,         &
     &  L_PC2_LBC,L_PC2_diag_sh, L_rad_deg, L_AC, L_VINT_TP,            &
     &  A_ASSIM_START_MIN, A_ASSIM_END_MIN,                             &
     &  NPMSL_HEIGHT,                                                   &
     &  LMASS_CORR , LQT_CORR, LEMQ_PRINT ,A_ENERGYSTEPS,               &
     &  L_GWD,L_USE_USSP,                                               &
     &  L_Murk, L_MURK_ADVECT, L_MURK_SOURCE, L_MURK_BDRY,              &
     &  L_MURK_RAD, L_murk_lbc, L_BL_TRACER_MIX, L_int_uvw_lbc,         &
     &  L_DUST,L_CAM_DUST,L_SULPC_SO2,L_SULPC_DMS,L_SULPC_OZONE,        &
     &  L_SULPC_SO2_O3_NONBUFFERED, L_SO2_SURFEM, L_SO2_HILEM,          &
     &  L_SO2_NATEM,                                                    &
     &  L_SULPC_ONLINE_OXIDANTS,                                        &
     &  L_DMS_EM, L_DMS_EM_INTER,                                       &
     &  L_DMS_Ointer,                                                   &
     &  L_DMS_Liss_Merlivat, L_DMS_Wanninkhof, L_DMS_Nightingale,       &
     &  L_SULPC_NH3, L_NH3_EM, L_SOOT, L_SOOT_SUREM, L_SOOT_HILEM,      &
     &  L_BIOMASS, L_BMASS_SUREM, L_BMASS_HILEM,                        &
     &  L_PSSWND, L_DUST2OCN,                                           &
     &  L_RHCPT, L_CCRAD, L_3D_CCA, L_3D_CCW, L_PHASE_LIM,              &
     &  L_CO2_INTERACTIVE,L_CO2_EMITS,l_cable,                          &
! rml 1/7/13
     &  L_CO2_TRACER,                                                   &
!kdcorbin, 08/10
     &  L_CO2_RADIATION,L_TRACER_MASS,L_CO2_MASS,I_TRACERMASS_START,    &
     &  L_METHANE_LOSS,I_METHANE_TRACERS,L_MCF_LOSS,I_MCF_TRACERNUMBER, &
     &  L_RADON_DECAY,I_RADON_TRACERNUMBER,                             &  
     &  L_Q10, L_NEG_TSTAR, L_VEG_FRACS, L_TRIFFID, L_PHENOL,           &
     &  L_TRIF_EQ, L_NRUN_MID_TRIF, L_DISTURB,                          &
     &  CAN_MODEL, PHENOL_PERIOD, TRIFFID_PERIOD,                       &
     &  L_USE_SEASALT_INDIRECT, L_USE_BIOGENIC,                         &
     &  L_USE_SEASALT_DIRECT, L_USE_DUST, L_USE_SULPC_INDIRECT_SW,      &
     &  L_USE_SULPC_INDIRECT_LW, L_USE_SULPC_DIRECT,                    &
     &  L_USE_SULPHATE_AUTOCONV, L_USE_SEASALT_AUTOCONV, L_AUTO_DEBIAS, &
     &  L_USE_SULPHATE_SULPC, L_USE_SEASALT_SULPC,                      &
     &  L_OCFF, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_OCFF_AUTOCONV,        &
     &  L_USE_OCFF_SULPC, L_USE_OCFF_DIRECT, L_USE_OCFF_INDIRECT,       &
     &  L_USE_STOCHEM_CH4, L_USE_STOCHEM_O3,                            &
     &  L_MOD_BARKER_ALBEDO, L_USE_SPEC_SEA, L_MOD_K_FLUX, L_CTILE_FIX, &
     &  L_QUAD_SRC_FIX, L_USE_TPPS_OZONE, I_OZONE_INT,                  &
     &  L_ICOUNT,                                                       &
     &  L_use_orog_corr, L_use_grad_corr,                               &
     &  CALL_CHEM_FREQ, L_USE_SOOT_DIRECT, L_USE_SOOT_INDIRECT,         &
     &  L_USE_SOOT_AUTOCONV, L_USE_SOOT_SULPC, L_USE_BMASS_DIRECT,      &
     &  L_USE_BMASS_INDIRECT, L_USE_BMASS_AUTOCONV, L_USE_BMASS_SULPC,  &
     &  L_USE_ARCLBIOM, L_USE_ARCLBLCK,  L_USE_ARCLSSLT,                &
     &  L_USE_ARCLSULP, L_USE_ARCLDUST,  L_USE_ARCLOCFF, L_USE_ARCLDLTA,&
     &  L_ukca,                                                         &
     &  L_CTILE, L_RIVERS,L_INLAND, RIVER_STEP,                         &
     &  l_sice_meltponds, l_sice_scattering, l_sice_hadgem1a,           &
     &  L_USE_METHOX,Z_TOP,l_mr_physics1,l_mr_physics2,                 &
     &  L_INHOM_CLOUD,                                                  &
     &  L_radiation,L_FORCING,L_TIMESTEP,                               &
     &  L_RADIANCE, L_WENYI, L_bl, L_rain, L_hydrology,                 &
     &  L_TOP,L_PDM,L_USE_AOD,L_USE_CLEARRH,L_UPPER_STRATQ,             &
     &  L_SCVARY,L_VOLCTS,                                              &
     &  L_sice_heatflux, L_SOIL_SAT_DOWN,                               &
     &  l_anthrop_heat_src,                                             &
     &  L_USE_CARIOLLE,L_USE_OZONEINRAD,                                &
     &  L_RPSEED_READ, L_RPSEED_WRITE,                                  &
     &  H_SECT, A_ASSIM_MODE

      ! Control switches derived and set within the model, hence not
      ! passed in from the UMUI via namelist:

      ! Radiation
      LOGICAL :: L_SW_RADIATE  ! Activate SW radiation this timestep
      LOGICAL :: L_LW_RADIATE  ! Activate LW radiation this timestep
      LOGICAL :: L_SW_RADIATE_DIAG
! Activate fast SW radiation this timestep (3C)
      LOGICAL :: L_LW_RADIATE_DIAG
! Activate fast LW radiation this timestep (3C)
      LOGICAL :: L_SW_RADIATE_PROG
! Activate slow SW radiation this timestep (3C)
      LOGICAL :: L_LW_RADIATE_PROG
! Activate slow LW radiation this timestep (3C)
      LOGICAL :: Lexpand_ozone ! Convert zonal mean ozone to field

      ! convert zonal mean tpps ozone to field
      LOGICAL :: Lexpand_tpps_ozone

      INTEGER :: I_tpps_ozone_opts ! options for tropopause-based ozone

      ! size of super array holding all tracers
      Integer ::   super_array_size
      Integer ::   moisture_array_size

      COMMON/CNTLCATM2/                                                 &
       !  super tracer array size
     &  super_array_size, moisture_array_size,                          &
     &  L_SW_RADIATE,L_LW_RADIATE,                                      &
     &  L_SW_RADIATE_DIAG,L_LW_RADIATE_DIAG,                            &
     &  L_SW_RADIATE_PROG,L_LW_RADIATE_PROG,                            &
     &  Lexpand_ozone,                                                  &
     & Lexpand_tpps_ozone, I_tpps_ozone_opts
! options for tropopause-based ozone



      LOGICAL :: LTLEADS ! Switch for Leads temperature.
!                      ! If FALSE, they are assumed to be TFS
!                      ! Else they are prognostic.
! Default setting: leads temperatures are set to TFS
! HARDWIRE LTLEADS
      PARAMETER(LTLEADS =.FALSE.)
!*L -----------------COMDECK PHYSCONS----------------------------------
!
!  Purpose : contains physical constants required by the whole of the
!            model. It is made up of individual COMDECKS for sets of
!            of related constants, each routine can access one or
!            several of these COMDECKS seperately
!  System Component : F07
!  System task : Z
! END
!*----------------------------------------------------------------------
!*L------------------COMDECK C_A----------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Replace variable A by more meaningful name for
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Convert to Fixed/Free format. P. Selwood

      ! Mean radius of Earth in m.
      Real, Parameter  :: Earth_Radius = 6371000.

!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
!*L------------------COMDECK C_OMEGA------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable two_omega for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.3   12/10/01  two_omega initialised in SETCON     Terry Davies
!OMEGA IS MAGNITUDE OF EARTH'S ANGULAR VELOCITY
      REAL OMEGA                                                        &
     &,two_omega

       Common/Omega/Omega, two_omega
!  Angular speed of Earth's rotation Omega to be initialised in SETCON
!*----------------------------------------------------------------------
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
! C_KT_FT start

      REAL,PARAMETER:: KT2MS=1852.0/3600.0 ! Knots to m/s conversion
      REAL,PARAMETER:: FT2M =0.3048        ! Feet to meters conversion

! C_KT_FT end
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
!*L------------------COMDECK C_LAPSE ----------------------------------
      Real, Parameter :: lapse      = 0.0065  ! Near surface lapse rate
      Real, Parameter :: lapse_trop = 0.002   ! Tropopause lapse rate
!*---------------------------------------------------------------------
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

!*L   Subroutines called;

      EXTERNAL                                                          &
     &       TIME2SEC,                                                  &
     &       READFLDS,                                                  &
     &       T_INT,                                                     &

     &       TO_LAND_POINTS,                                            &
     &       SEC2TIME,TIME_DF,                                          &

     &       T_INT_C

!*
!*L   Local real arrays

      REAL                                                              &
     &       ANCIL1(P_FIELD),                                           &
                                ! Buffers to hold values of ancillary
!                               ! data for time interpolation.
     &       ANCIL2(P_FIELD),                                           &
                                !
     &       ANCIL_DATA(P_FIELD),                                       &
                                 ! Field of ancillary data held prior
!                               ! to selective updating.
     &       SNOW_CHANGE(P_FIELD),                                      &
                                  ! Fractional time of change of
!                               ! snow cover
     &       ICE_EXTENT(P_FIELD,2),                                     &
                                   ! Fractional time of change
!                               ! of ice cover
     &       PRES_VALUE(P_FIELD)                                        &
                                 ! Prescribed value of data when
!                               ! controlling field is zero.
     &,      NO_ICE_EXTENT(P_FIELD)                                     &
                                   ! Indicator for no sea ice
!                               ! =0 if ice cover
     &,       TSTAR_LAND(P_FIELD)                                       &
                                 !Temporary store for land surface temp.
     &,       TSTAR_SEA(P_FIELD)                                        &
                                 !as above, but for open sea.
     &,       TSTAR_SICE(P_FIELD)                                       &
                                 !as above, but for sea-ice.
     &,       TSTAR_SSI(P_FIELD) !as above, but for sea mean.
!*
!     Local variables

      INTEGER                                                           &
     &       I,                                                         &
                                !
     &       I1,                                                        &
                                !
     &       I2,                                                        &
                                !
     &       I3,                                                        &
     &       ID,                                                        &
                                !
     &       IM,                                                        &
                                !
     &       IY,                                                        &
                                !
     &       K,                                                         &
     &       L,                                                         &
                                    ! Land index
     &       FIELD,                                                     &
                                ! Current field number.
     &       FILE               !

      INTEGER                                                           &
     &       INTERVAL,                                                  &
                                ! Interval between data times
     &       STEP,                                                      &
                                ! Number of data times skipped.
     &       MONTHS,                                                    &
                                ! Used in calculation of position
!                               ! of data required.
     &       HOURS,                                                     &
                                !
     &       PERIOD,                                                    &
                                ! Period of periodic data
     &       START_MONTH,                                               &
                                !
     &       LEVEL,                                                     &
                                !
     &       NFTIN,                                                     &
                                ! Current FTN number for ancillary field
     &       ANCIL_REF_DAYS,                                            &
                                ! Ancil.reference time in whole days
     &       ANCIL_REF_SECS,                                            &
                                ! Ancil.reference time in extra seconds
     &       DAY,SEC,                                                   &
                                ! Times relative to reference time
     &       DAY1,SEC1,                                                 &
                                ! Times relative to reference time
     &       INCR_SEC,                                                  &
                                ! Increment in sec
     &       LEN                                                        &
     &      ,IEND                                                       &
     &      ,II,ROW_LENGTH,J
      INTEGER                                                           &
     &       I_YEAR1,                                                   &
                                 ! Copy of Curent Model Time year
     &       I_MONTH1,                                                  &
                                 !   "      "     "          month
     &       I_DAY1,                                                    &
                                 !   "      "     "          day
     &       I_HOUR1             !   "      "     "          hour

      INTEGER                                                           &
     &       UPDATE_MONTHS      ! update frequency (months) if Gregorian
      LOGICAL                                                           &
     &       LGREG_MONTHLY      ! True for Gregorian monthly updating
!
! *IF -DEF,CAL360
!
      INTEGER                                                           &
     &       I_YEAR_BASIS,                                              &
                                      ! Basis Model Time
     &       I_MONTH_BASIS,                                             &
                                      !   "     "     "
     &       I_DAY_BASIS,                                               &
                                      !   "     "     "
     &       I_HOUR_BASIS,                                              &
                                      !   "     "     "
     &       I_MINUTE_BASIS,                                            &
                                      !   "     "     "
     &       I_SECOND_BASIS,                                            &
                                      !   "     "     "
     &       I_DAY_NUMBER_BASIS
!
! *ENDIF
!

      INTEGER                                                           &
     &       I_YEAR_REF,                                                &
                                      ! Reference Time
     &       I_MONTH_REF,                                               &
                                      !    "       "
     &       I_DAY_REF,                                                 &
                                      !    "       "
     &       I_HOUR_REF,                                                &
                                      !    "       "
     &       I_MINUTE_REF,                                              &
                                      !    "       "
     &       I_SECOND_REF             !    "       "


      LOGICAL                                                           &
     &       LINTERPOLATE,                                              &
                                ! Indicates whether time
!                               ! interpolation needed.
     &       LT_INT_C,                                                  &
                                ! Indicates use of controlled time
!                               ! interpolation
     &       LMISMATCH,                                                 &
                                ! Used in header checks
     &       LICE_FRACTION,                                             &
                                !
     &       LSNOW_DEPTH,                                               &
                                !
     &       SINGLE_TIME,                                               &
                                ! Indicates that only one time is
!                               ! available in data set
     &       PERIODIC,                                                  &
                                ! Data set is periodic
     &       REGULAR                                                    &
                                ! Interval between data times in
!                               ! dataset is regular in model timesteps.
     &       ,LICE_DEPTH

      REAL                                                              &
     &       ZERO,                                                      &
                              !
     &       TIME1,                                                     &
                              ! Times if data used in time interpolation
     &       TIME2,                                                     &
                              !
     &       TIME                                                       &
                        !Target time for time interpolation
     &       ,LAT_P          ! latitude of point
     LOGICAL LAMIP2X ! Local flag

!L    Internal structure

!L  List of Atmosphere & Slab Ancillary fields.
! CANCLSTA start
!
! Purpose : Cross-Reference List of Ancillary Fields
!           Atmosphere. Ancillary Reference numbers,
!           Stash Codes, Model Codes and Logical file Numbers
!
! User changes:
!   05/10 kdcorbin - added tracer fluxes
!
! -------------------------------------------------------------------
!   Column A  : Ancillary Reference Number
!          B  : Internal Model Number
!          C  : Stash Item code
!          D  : Logical file number
!          E  : Field Description
!
!   A  B    C   D  E
!   ----------------------------------------------------------------
!   1  1   30   9  Land Sea Mask
!   2  1   33  10  Orography
!   3  1   34  10  Orographic Variance
!   4  1   35  10  Orographic gradient XX
!   5  1   36  10  Orographic gradient XY
!   6  1   37  10  Orographic gradient YY
!   7  1   60   1  Ozone
!   8              Not used
!   9  1   23   2  Snow Depth
!  10  1   20   3  Deep Soil Temperature
!  11  1   40   4  Vol SMC at Wilting
!  12  1   41   4  Vol SMC at Critical Point
!  13              Not used
!  14  1   43   4  Vol SMC at Saturation
!  15  1   44   4  Saturated Soil Conductivity
!  16              Not used
!  17  1   46   4  Thermal Capacity
!  18  1   47   4  Thermal Conductivity
!  19  1   50   5  Vegetation Fraction
!  20  1   51   5  Root Depth
!  21  1   52   5  Snow Free Surface Albedo
!  22  1   53   5  Deep Snow Surface Albedo
!  23  1   54   5  Surface Resistance to Evaporation
!  24  1   55   5  Surface Capacity
!  25  1   56   5  Infiltration Factor
!  26  1   26   5  Surface Roughness (vegetation)
!  27  1   31   7  Sea-ice Fraction
!  28  1   24   6  Sea-surface Temperature
!  29  1   32   7  Sea-ice Thickness
!  30  1   28   8  Surface Currents : u-component
!  31  1   29   8  Surface Currents : v-component
!  32  1   93   9  Runoff coastal outflow point
!  33              Not used (slab model)
!  34              Not used
!  35  1   48   4  Saturated soil water suction
!  36  1    9   2  Soil moisture in layers
!  37              Not used (SLAB) 
!  38              Not used (SLAB) 
!  39  1   58  12  Sulphur dioxide emission
!  40  1   59  12  Dimethyl sulphur emission
!  41  1   88  13  Sulphate aerosol mass mixing ratio
!  42  1   87  13  Sulphuric acid aerosol mass mixing ratio
!  43  1   85  13  Soot aerosol mass mixing ratio
!  44  1   57  14  Multi-level murk source term emission
!  45  1   90  14  Multi-level murk concentration
!  46  1   17  10  Silhouette area of orography (orog. roughness)
!  47  1   18  10  Peak to trough height (for orog. roughness scheme)
!  48  1  301  15  User ancillary field 1
!  49  1  302  15  User ancillary field 2
!  50  1  303  15  User ancillary field 3
!  51  1  304  15  User ancillary field 4
!  52  1  305  15  User ancillary field 5
!  53  1  306  15  User ancillary field 6
!  54  1  307  15  User ancillary field 7
!  55  1  308  15  User ancillary field 8
!  56  1  309  15  User ancillary field 9
!  57  1  310  15  User ancillary field 10
!  58  1  311  15  User ancillary field 11
!  59  1  312  15  User ancillary field 12
!  60  1  313  15  User ancillary field 13
!  61  1  314  15  User ancillary field 14
!  62  1  315  15  User ancillary field 15
!  63  1  316  15  User ancillary field 16
!  64  1  317  15  User ancillary field 17
!  65  1  318  15  User ancillary field 18
!  66  1  319  15  User ancillary field 19
!  67  1  320  15  User ancillary field 20
!  68  1  127  12  NH3 (ammonia) aerosol emission
!  69  1  128  23  Surface fresh soot aerosol emission
!  70  1  129  23  Elevated fresh soot aerosol emission
!  71              Not used
!  72  1  121  17  Natural Sulphur dioxide emissions
!  73  1  122  18  OH concentrations
!  74  1  123  18  HO2 concentrations
!  75  1  124  18  H2O2 concentrations
!  76  1  125  18  Ozone (CHEM) concentrations
!  77  1  126  12  Sulphur dioxide high level emission
!  78  1  251  24  Surface CO2 emissions
!  79  1  207   4  Clapp-Hornberger parameter
!  80  1  208   5  Leaf area index of vegetated fraction
!  81  1  209   5  Canopy height of vegetated fraction
!  82  1  160  19  Aerosol data for radiative forcing of climate change
!  83  1  216  20  Initial fractions of surface types
!  84  1  217  21  Initial leaf area index of plant functional types
!  85  1  218  21  Initial canopy height of plant functional types
!  86  1  213  21  Initial gridbox mean canopy conductance
!  87  1  219  22  Fraction of vegetation subject to disturbance
!  88  1  220   4  Snow free albedo of bare soil
!  89  1  223   4  Soil carbon content
!  90  1  321  16  User ancillary multi 1
!  91  1  322  16  User ancillary multi 2
!  92  1  323  16  User ancillary multi 3
!  93  1  324  16  User ancillary multi 4
!  94  1  325  16  User ancillary multi 5
!  95  1  326  16  User ancillary multi 6
!  96  1  327  16  User ancillary multi 7
!  97  1  328  16  User ancillary multi 8
!  98  1  329  16  User ancillary multi 9
!  99  1  330  16  User ancillary multi 10
! 100  1  331  16  User ancillary multi 11
! 101  1  332  16  User ancillary multi 12
! 102  1  333  16  User ancillary multi 13
! 103  1  334  16  User ancillary multi 14
! 104  1  335  16  User ancillary multi 15
! 105  1  336  16  User ancillary multi 16
! 106  1  337  16  User ancillary multi 17
! 107  1  338  16  User ancillary multi 18
! 108  1  339  16  User ancillary multi 19
! 109  1  340  16  User ancillary multi 20
! 110  1  341  25  Tropopause-based Ozone
! 111  1  505  26  Land Fraction
! 112  1  418  27  Dust parent soil clay fraction
! 113  1  419  27  Dust parent soil silt fraction
! 114  1  420  27  Dust parent soil sand fraction
! 115  1  421  27  Dust soil mass fraction div 1
! 116  1  422  27  Dust soil mass fraction div 2
! 117  1  423  27  Dust soil mass fraction div 3
! 118  1  424  27  Dust soil mass fraction div 4
! 119  1  425  27  Dust soil mass fraction div 5
! 120  1  426  27  Dust soil mass fraction div 6
! 121  1  130  28  Biomass surface emissions
! 122  1  131  28  Biomass elevated emissions
! 123  1  132  29  DMS concentration in seawater
! 124  1  153  30  River Water Storage
! 125  1  151  31  Riv Channel Sequence
! 126  1  152  31  Riv Channel Direction
!
! 127-154 not used in UM
!
! 155  1    5  10  Orographic X gradient
! 156  1    6  10  Orographic Y gradient 
! 157  1  351  38  Climatalogical biogenic aerosol mmr
! 158  1  352  39  Clim Biomass-burning (fresh) mmr
! 159  1  353  39  Clim Biomass-burning (aged) mmr
! 160  1  354  39  Clim Biomass-burning (in-cloud) mmr
! 161  1  355  40  Clim Black Carbon (fresh) mmr
! 162  1  356  40  Clim Black Carbon (aged) mmr
! 163  1  357  41  Clim Sea-salt (film mode) npm3
! 164  1  358  41  Clim Sea-salt (jet mode) npm3
! 165  1  359  42  Clim Sulphate (accumulation mode) mmr
! 166  1  360  42  Clim Sulphate (Aitken mode) mmr
! 167  1  361  42  Clim Sulphate (dissolved) mmr
! 168  1  362  43  Clim Dust size division 1 mmr
! 169  1  363  43  Clim Dust size division 2 mmr
! 170  1  364  43  Clim Dust size division 3 mmr
! 171  1  365  43  Clim Dust size division 4 mmr
! 172  1  366  43  Clim Dust size division 5 mmr
! 173  1  367  43  Clim Dust size division 6 mmr
! 174  1  368  44  Clim Organic Carbon from Fossil Fuels (fresh) mmr
! 175  1  369  44  Clim Organic Carbon from Fossil Fuels (aged) mmr
! 176  1  370  44  Clim Organic Carbon from Fossil Fuels (in-cloud) mmr
! 177  1  371  45  Clim Delta Aerosol mmr
! 178  1  480  46  Prognostic Ozone Tracer Cariolle Scheme
! 179  1  481  46  Cariolle Ozone Production - Loss (P-L)
! 180  1  482  46  Cariolle Ozone P-L wrt Ozone Mix Ratio
! 181  1  483  46  Cariolle Ozone Volume Mixing Ratio
! 182  1  484  46  Cariolle Ozone P-L wrt Temperature
! 183  1  485  46  Cariolle Ozone Clim Temp
! 184  1  486  46  Cariolle Ozone P-L wrt Ozone above a point
! 185  1  487  46  Cariolle Ozone Column above a point
! 186  1  134  47  Surface fresh fossil-fuel OC aerosol emissions
! 187  1  135  47  Elevated fresh fossil-fuel OC aerosol emissions
! 188  1  100  48  Flux of Tracer 1
! 189  1  100  49  Flux of Tracer 2
! 190  1  100  50  Flux of Tracer 3
! 191  1  100  51  Flux of Tracer 4
! 192  1  100  52  Flux of Tracer 5
! 193  1  100  53  Flux of Tracer 6
! 194  1  100  54  Flux of Tracer 7
! 195  1  100  55  Flux of Tracer 8
! 196  1  100  56  Flux of Tracer 9
! 197  1  100  57  Flux of Tracer 10
! 198  1  100  58  Flux of Tracer 11
! 199  1  100  59  Flux of Tracer 12
! 200  1  100  60  Flux of Tracer 13
! 201  1  100  61  Flux of Tracer 14
! 202  1  100  62  Flux of Tracer 15
! 203  1  100  63  Flux of Tracer 16
! 204  1  100  64  Flux of Tracer 17
! 205  1  100  65  Flux of Tracer 18
! 206  1  100  66  Flux of Tracer 19
! 207  1  100  67  Flux of Tracer 20
!  ------------------------------------------------------------------
! CANCLSTA end

!L  1.  Initialisation for atmosphere
      ICODE=0
      IOUNIT=0
      SMC_UPDATED=.FALSE.
      UPDATE_MONTHS=0
      INCR_SEC = 0

!     Set up surface temperatures:

       IF(L_CTILE)THEN
         DO I=1,P_FIELD
            TSTAR_LAND(I)=TSTAR_LAND_CTILE(I)
            TSTAR_SEA(I)=TSTAR_SEA_CTILE(I)
            TSTAR_SICE(I)=TSTAR_SICE_CTILE(I)
            IF(ICE_FRACTION(I) <= 0.0)THEN
              TSTAR_SSI(I)=TSTAR_SEA(I)
            ELSE
              TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                &
     &          +(1.0-ICE_FRACTION(I))*TSTAR_SEA(I)
            ENDIF
         ENDDO
       ELSE
         DO I=1,P_FIELD
            TSTAR_LAND(I)=TSTAR(I)
            TSTAR_SSI(I)=TSTAR(I)
         ENDDO
       ENDIF


!     Initialise ANCIL1/2. Includes Halos for mpp runs.
      DO I=1,P_FIELD
        ANCIL1(I)=0.0
        ANCIL2(I)=0.0
      ENDDO
!L  1.1 Set logical UPDATE for each ancillary field independently

      DO FIELD=1,NANCIL_FIELDS


        UPDATE(FIELD)=.FALSE.
        IF(STEPS(FIELD) /= 0) THEN
!         UPDATE(FIELD)=MOD(A_STEP,STEPS(FIELD)) == 0
          UPDATE(FIELD)=(MOD(A_STEP+OFFSET_STEPS,STEPS(FIELD)) == 0     &
     &                   .OR.A_STEP == 0)                               &
     &                    .AND.FIELDCODE(1,FIELD) >  0                  &
     &                     .AND.D1_ANCILADD(FIELD) >  1
        END IF

!L  1.05 Copy ancillary updating reference time to local variables
      I_YEAR_REF   = ANCIL_REFTIME(1)
      I_MONTH_REF  = ANCIL_REFTIME(2)
      I_DAY_REF    = ANCIL_REFTIME(3)
      I_HOUR_REF   = ANCIL_REFTIME(4)
      I_MINUTE_REF = ANCIL_REFTIME(5)
      I_SECOND_REF = ANCIL_REFTIME(6)
!L       and convert to reference days & secs
! DEPENDS ON: time2sec
            CALL TIME2SEC(I_YEAR_REF,I_MONTH_REF,I_DAY_REF,             &
     &                    I_HOUR_REF,I_MINUTE_REF,I_SECOND_REF,         &
     &                    0,0,ANCIL_REF_DAYS,ANCIL_REF_SECS,LCAL360)

!
      IF (.NOT. LCAL360) THEN

!L  1.11 Set logical UPDATE for Gregorian calender updates at monthly
!L       or yearly intervals. NB STEPS value set to 1 day in INANCILA
        IF(FIELDCODE(1,FIELD) == 1.OR.FIELDCODE(1,FIELD) == 2) THEN
          MONTHS=I_MONTH+I_YEAR*12-(I_MONTH_REF+I_YEAR_REF*12)
          UPDATE_MONTHS= FIELDCODE(2,FIELD)*                            &
     &     ((3-FIELDCODE(1,FIELD))/2 *12+ 1-(3-FIELDCODE(1,FIELD))/2)
          UPDATE(FIELD)=MOD(MONTHS,UPDATE_MONTHS) == 0.AND.I_DAY == 1
        END IF
      END IF !  (.NOT.LCAL360)
!


      END DO

!L 1.2 Allow for dependencies between fields
! Sea surface temperature must be updated when sea ice is updated

      UPDATE(28)=UPDATE(27).OR.UPDATE(28)

! Both surface current components must be updated together

      UPDATE(30)=UPDATE(30).OR.UPDATE(31)
      UPDATE(31)=UPDATE(30)

!L Select method of time interpolation for SST. The interpolation
!L allows for sea ice if ice data is available at the same times
!L as the temperature data. Otherwise linear interpolation is used.

      LT_INT_C=.TRUE.

      IF(UPDATE(28)) THEN
      IF(FIXHD(10,FILEANCIL(27)) == 0) LT_INT_C=.FALSE.
        IF(LT_INT_C) THEN
        DO I=21,41
          IF(FIXHD(I,FILEANCIL(27)) /= FIXHD(I,                         &
     &      FILEANCIL(28))) THEN
            LT_INT_C=.FALSE.
            WRITE(6,*)' WARNING:controlled time interpolation for SST', &
     &      ' not available: Mismatch in SST and SEA-ICE ancillary data'&
     &     ,' times in FIXED HEADER'
            WRITE(6,*)' position=',I,' SEA-ICE=',FIXHD(I,FILEANCIL(27))
            WRITE(6,*)' position=',I,' SST=',FIXHD(I,FILEANCIL(28))
          END IF
        END DO
        ENDIF
      END IF


! Read in fractional land field
! Set up global fractional land field
         IF(L_CTILE)THEN
           L=0
           DO I=1,P_FIELD
             FLAND_G(I)=0.0
             IF(LAND(I))THEN
               L=L+1
               FLAND_G(I)=FLAND_CTILE(L)
            ENDIF
           ENDDO
         ELSE
           DO I=1,P_FIELD
             IF(LAND(I))THEN
               FLAND_G(I)=1.0
             ELSE
               FLAND_G(I)=0.0
            ENDIF
           ENDDO
         ENDIF
!

        DO I=1,P_FIELD
          SEA(I)=.FALSE.
          IF(FLAND_G(I) <  1.0)SEA(I)=.TRUE.
        ENDDO

!L Loop over ancillary fields(atmosphere)

      DO FIELD=1,NANCIL_FIELDS
       ! Turn this off for all but SST and sea-ice
        LAMIP2X = LAMIPII .and. (FIELD == 27 .or. FIELD == 28 .or. FIELD == 29)

       ! Turn this off for all but SST and sea-ice
        LAMIP2X = LAMIPII .and. (FIELD == 27 .or. FIELD == 28 .or. FIELD == 29)

        LICE_DEPTH=field == 29  ! required for LAMIPII

      IF (UPDATE(FIELD)) THEN  ! (1st level IF)
        FILE=FILEANCIL(FIELD)
        NFTIN=FTNANCIL(FILE)

       IF(LICE_DEPTH.AND.LAMIP2X) THEN

! Uses ice fraction set earlier in field loop.
! WARNING this will fail if the order of ancillary fields is ever
! changed so that ice-depth preceeds ice fraction
! Note : For complete sea ice cover
!        Arctic ice depth    = 2m
!        Antarctic ice depth = 1m
! For ice concentrations less than 1. ice depth is 1 or 2 times conc.
! This results in similar values to those from runs using ancillary
! files containing ice depths set to 1 or 2m.

          ROW_LENGTH=P_FIELD/P_ROWS
          DO I=1,P_ROWS
! work out latitude in radians
            LAT_P=FIRST_LAT+NS_SPACE*(I+datastart(2)-Offy-1)
            DO J=1,ROW_LENGTH
              II=J+(I-1)*ROW_LENGTH
              ANCIL_DATA(II)=0.0
              IF (ICE_FRACTION(II) >  0.0) THEN
                IF (LAT_P >  0.0) THEN   ! Arctic ice depth
                  ANCIL_DATA(II)=2.*ICE_FRACTION(II)
                ELSE                     ! Antarctic ice depth
                  ANCIL_DATA(II)=1.*ICE_FRACTION(II)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
!L     Sea ice thickness
!L       Update over all sea points (all sea ice points are the only
!L       ones strictly required, but this cannot be determined easily)

          DO I=1,P_FIELD
            IF(SEA(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            END IF
          END DO
       ELSE
!     Update required for field

        WRITE(6,*)'REPLANCA: UPDATE REQUIRED FOR FIELD',FIELD

          IF ( FIXHD(10,FILE)  <   0 .OR. FIXHD(10,FILE)  >   2 ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Error in fixed header(10) of ancillary&
     & file                           '
            RETURN
          ENDIF

!L    Check whether more than one data time available in data set

        SINGLE_TIME=FIXHD(10,FILE) == 0

!L    Set default values for time interpolation

        LINTERPOLATE=.TRUE.
        IF(SINGLE_TIME) THEN
          LINTERPOLATE=.FALSE.
        END IF

        IF (FIELD >  9 .AND. FIELD <  19) THEN
          LINTERPOLATE=.FALSE.
        END IF

!L 2.1 Find position of input record

!L    Default settings of search parameters if only one time present

        IF(SINGLE_TIME) THEN
          STEP=0
        ELSE


          LGREG_MONTHLY=.FALSE.
!
      IF (.NOT. LCAL360) THEN
          IF(FIELDCODE(1,FIELD) == 1.OR.FIELDCODE(1,FIELD) == 2) THEN
            LGREG_MONTHLY=.TRUE.
            UPDATE_MONTHS= FIELDCODE(2,FIELD)*                          &
     &      ((3-FIELDCODE(1,FIELD))/2 *12+ 1-(3-FIELDCODE(1,FIELD))/2)
          END IF
      END IF
!


          PERIODIC=FIXHD(10,FILE) == 2
          REGULAR=.TRUE.

!
      IF (.NOT. LCAL360) THEN
          REGULAR=FIXHD(35,FILE) == 0.AND.FIXHD(36,FILE) == 0
! i.e. data at intervals of days/hours & non-periodic
          IF(PERIODIC) REGULAR=REGULAR.AND.FIXHD(37,FILE) == 0
! i.e. data at intervals of hours & periodic
      END IF
!

!         Error checking on time information.

          IF ( FIXHD(35,FILE)  <   0 .OR.                               &
     &         FIXHD(36,FILE)  <   0 .OR. FIXHD(36,FILE)  >   12 .OR.   &
     & REGULAR .AND. ( FIXHD(37,FILE)  <   0 .OR. FIXHD(37,FILE)  >   31&
     &  .OR. FIXHD(38,FILE)  <   0 .OR. FIXHD(38,FILE)  >   24 ) ) THEN
!           FIXHD(39-40) are not used by REPLANCA.
!           FIXHD(35-37) have already been used if not CAL360.
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Error in validity time interval given &
     &in ancillary file (FIXHD(35-38))'
            RETURN
          ENDIF

          IF ( FIXHD(21,FILE)  <   0 .AND. .NOT. PERIODIC               &
     &  .OR. .NOT. ( REGULAR .AND. PERIODIC ) .AND.                     &
!    !  If it is REGULAR & PERIODIC more detailed check is applied below
     &     ( FIXHD(22,FILE)  <   0 .OR. FIXHD(22,FILE)  >   12 .OR.     &
     &       FIXHD(23,FILE)  <   0 .OR. FIXHD(23,FILE)  >   31 .OR.     &
     &       FIXHD(24,FILE)  <   0 .OR. FIXHD(24,FILE)  >   24 .OR.     &
     &       FIXHD(25,FILE)  <   0 .OR. FIXHD(25,FILE)  >   60 .OR.     &
     &       FIXHD(26,FILE)  <   0 .OR. FIXHD(26,FILE)  >   60 ) ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Error in first validity time given in &
     & ancillary file (FIXHD(21-26))  '
            RETURN
          ENDIF

          IF(.NOT.PERIODIC) THEN

!L            If data taken from full time series of input data.

! DEPENDS ON: time2sec
            CALL TIME2SEC(I_YEAR,I_MONTH,I_DAY,I_HOUR                   &
     &                    ,I_MINUTE,I_SECOND                            &
     &                    ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC        &
     &                    ,LCAL360)


!L Adjust time to middle of updating interval

            IF(.NOT.LGREG_MONTHLY) THEN
              SEC=SEC+STEPS(FIELD)*1800/STEPS_PER_HR

!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
              IF (A_STEP == 0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
! DEPENDS ON: time_df
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF

            ELSE
              IM=MOD(I_MONTH+UPDATE_MONTHS-1,12) + 1
              IY=I_YEAR+(I_MONTH+UPDATE_MONTHS-1)/12
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,I_DAY,I_HOUR                          &
     &                    ,I_MINUTE,I_SECOND                            &
     &                    ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY1,SEC1      &
     &                    ,LCAL360)
              IF (MOD(DAY+DAY1,2) == 0) THEN
                DAY=(DAY+DAY1)/2
                SEC=(SEC+SEC1)/2
              ELSE
                DAY=(DAY+DAY1-1)/2
                SEC=(SEC+SEC1+86400)/2
              ENDIF
!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
              IF (A_STEP == 0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
! DEPENDS ON: time_df
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF
            ENDIF


            IF(REGULAR) THEN
!L 2.1.1  Standard cases:360 day calender;
!L 2.1.1  or Gregorian calendar with
!L        interval between data times in days or hours
!L        updating interval may be regular in model timesteps,
!L        or (LGREG_MONTHLY=T) irregular in model timesteps,

              HOURS=SEC/3600+DAY*24
!L FInd time(in hours) of first ancillary data on file
! DEPENDS ON: time2sec
              CALL TIME2SEC(FIXHD(21,FILE),FIXHD(22,FILE),              &
     &                   FIXHD(23,FILE),FIXHD(24,FILE),                 &
     &                   FIXHD(25,FILE),FIXHD(26,FILE),                 &
     &                   ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,         &
     &                   LCAL360)
              HOURS=HOURS-SEC/3600-DAY*24

              IF(HOURS <  0) THEN
                ICODE=400+FIELD
           CMESSAGE='REPLANCA: Current time precedes start time of data'
                RETURN
              END IF

!L FInd interval(in hours) between ancillary data on file
              INTERVAL=FIXHD(35,FILE)*8640+FIXHD(36,FILE)*720+          &
     &                FIXHD(37,FILE)*24+FIXHD(38,FILE)

! Do not interpolate in time if data time exactly matches model time

              IF(MOD(HOURS,INTERVAL) == 0) THEN
                LINTERPOLATE=.FALSE.
              END IF

              STEP=HOURS/INTERVAL
              TIME=REAL(HOURS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            ELSE

!L 2.1.2 Gregorian calender;ancillary data interval is in months or
!L       years,which is irregular in model timesteps.
!L original code is inaccurate for this section - corrected code under
!L LAMIPII makes use of dates in lookup headers
!L For a real calendar year the mid-point of each month is different
!L in terms of its hour and day. The old inaccurate method assumes
!L the hour and day are taken from the fixhd values. These are only
!L usually correct for the first month on the ancillary file.


!L Adjust YMD time to middle of updating interval

              I_YEAR1=I_YEAR
              I_MONTH1=I_MONTH
              I_DAY1=I_DAY
              I_HOUR1=I_HOUR
! DEPENDS ON: sec2time
              CALL SEC2TIME(DAY,SEC,ANCIL_REF_DAYS,ANCIL_REF_SECS,      &
     &                     I_YEAR,I_MONTH,I_DAY,                        &
     &                     I_HOUR,I_MINUTE,I_SECOND,I_DAY_NUMBER,       &
     &                     LCAL360)


!L FInd interval(in months) between ancillary data on file
              INTERVAL=FIXHD(35,FILE)*12+FIXHD(36,FILE)
              MONTHS=I_YEAR*12+I_MONTH
              START_MONTH=FIXHD(21,FILE)*12+FIXHD(22,FILE)
              MONTHS=MONTHS-START_MONTH
!  Check for time within month
           IF (LAMIP2X) THEN   ! corrected code uses pp header
              STEP=MONTHS/INTERVAL
              I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
              I1=I2+LOOKUP_START(FILE)-1
! Check against day and hour of actual lookup header not first field
              IF((I_DAY*24+I_HOUR) <                                    &
     &           (LOOKUP(3,I1)*24+LOOKUP(4,I1))) THEN
                MONTHS=MONTHS-1
              END IF
           ELSE              ! old less accurate code uses FIXHD
              IF((I_DAY*24+I_HOUR) <                                    &
     &           (FIXHD(23,FILE)*24+FIXHD(24,FILE))) THEN
                MONTHS=MONTHS-1
              END IF
           ENDIF ! LAMIP2X

              IF(MONTHS <  0) THEN
                ICODE=400+FIELD
           CMESSAGE='REPLANCA: Current time precedes start time of data'
                RETURN
              END IF


!L Adjust YMD time back to start of updating interval

              I_YEAR=I_YEAR1
              I_MONTH=I_MONTH1
              I_DAY=I_DAY1
              I_HOUR=I_HOUR1



              STEP=MONTHS/INTERVAL

           IF (LAMIP2X) THEN       ! corrected code
              TIME=REAL(SEC)/3600+REAL(DAY*24)
! correct calculation of dates uses lookup table dates not fixhd date
              I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
              I1=I2+LOOKUP_START(FILE)-1
              I_YEAR1=lookup(1,i1)
              I_MONTH1=lookup(2,i1)
              I_DAY1=lookup(3,i1)
              I_HOUR1=lookup(4,i1)
! DEPENDS ON: time2sec
              CALL TIME2SEC(I_YEAR1,I_MONTH1,I_DAY1,I_HOUR1,            &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME1=REAL(SEC)/3600+REAL(DAY*24)
! I1+1 correct pointer to next field as only one field in ancil file
              I_YEAR1=lookup(1,i1+1)
              I_MONTH1=lookup(2,i1+1)
              I_DAY1=lookup(3,i1+1)
              I_HOUR1=lookup(4,i1+1)
! DEPENDS ON: time2sec
              CALL TIME2SEC(I_YEAR1,I_MONTH1,I_DAY1,I_HOUR1,            &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME2=REAL(SEC)/3600+REAL(DAY*24)

           ELSE   ! LAMIP2X test - old inaccurate code using FIXHD
! NB INTERVAL may be > 1 month
              MONTHS=STEP*INTERVAL
! Calculate data times for time interpolation
              TIME=REAL(SEC)/3600+REAL(DAY*24)
              IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
              IY=FIXHD(21,FILE)+(MONTHS+FIXHD(22,FILE)-1)/12
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),        &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME1=REAL(SEC)/3600+REAL(DAY*24)
              IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
              IY=FIXHD(21,FILE)+(MONTHS+INTERVAL+FIXHD(22,FILE)-1)/12
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),        &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME2=REAL(SEC)/3600+REAL(DAY*24)
           ENDIF     ! end LAMIP2X test

! Do not interpolate in time if data time exactly matches model time

              IF(TIME == TIME1) THEN
                LINTERPOLATE=.FALSE.
              END IF

            ENDIF ! End of REGULAR/not REGULAR

          ELSE  ! PERIODIC data

!L 2.2   If data is taken from ancillary periodic data.

! DEPENDS ON: time2sec
            CALL TIME2SEC(I_YEAR,I_MONTH,I_DAY,I_HOUR,                  &
     &                     I_MINUTE,I_SECOND,                           &
     &                     ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,       &
     &                     LCAL360)


!L Adjust time to middle of updating interval

            IF(.NOT.LGREG_MONTHLY) THEN
              SEC=SEC+STEPS(FIELD)*1800/STEPS_PER_HR

!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
              IF (A_STEP == 0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
! DEPENDS ON: time_df
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF

            ELSE
              IM=MOD(I_MONTH+UPDATE_MONTHS-1,12) + 1
              IY=I_YEAR+(I_MONTH+UPDATE_MONTHS-1)/12
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,I_DAY,I_HOUR                          &
     &                    ,I_MINUTE,I_SECOND                            &
     &                    ,ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY1,SEC1      &
     &                    ,LCAL360)
              IF (MOD(DAY+DAY1,2) == 0) THEN
                DAY=(DAY+DAY1)/2
                SEC=(SEC+SEC1)/2
              ELSE
                DAY=(DAY+DAY1-1)/2
                SEC=(SEC+SEC1+86400)/2
              ENDIF
!  If start-up, adjust for offset of reference time from initial time,
!  & update with values for half a period before first standard update.
              IF (A_STEP == 0) THEN
                DAY1 = DAY
                SEC1 = SEC
             INCR_SEC=-3600*MOD(OFFSET_STEPS,STEPS(FIELD))/STEPS_PER_HR
! DEPENDS ON: time_df
                CALL TIME_DF(DAY1,SEC1,0,INCR_SEC,DAY,SEC)
              END IF
            ENDIF


!L Adjust YMD time to middle of updating interval

            I_YEAR1=I_YEAR
            I_MONTH1=I_MONTH
            I_DAY1=I_DAY
            I_HOUR1=I_HOUR
! DEPENDS ON: sec2time
            CALL SEC2TIME(DAY,SEC,ANCIL_REF_DAYS,ANCIL_REF_SECS,        &
     &                     I_YEAR,I_MONTH,I_DAY,                        &
     &                     I_HOUR,I_MINUTE,I_SECOND,I_DAY_NUMBER,       &
     &                     LCAL360)



            IF (REGULAR) THEN
!L 2.2.1 Standard cases:1) 360 day calender, with allowed periods of
!L       1 day, 1 month or 1 year;
!L
!L       2) Gregorian calender with update in hours,and period of
!L       data 1 day.
!L
!L       For both updating interval and number of
!L       data times to be skipped in data set calculated in hours.

              HOURS=SEC/3600+DAY*24
              INTERVAL=FIXHD(35,FILE)*8640+FIXHD(36,FILE)*720+          &
     &                FIXHD(37,FILE)*24+FIXHD(38,FILE)

              PERIOD=INTHD(3,FILE)*INTERVAL

!L   Do not allow non-standard periods
      IF (LCAL360) THEN
              IF(PERIOD /= 8640.AND.PERIOD /= 720.AND.PERIOD /= 24)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              ENDIF
      ELSE
              IF(PERIOD /= 24)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              ENDIF
      ENDIF
              IF(PERIOD == 24)THEN
! Ancillary data interval in hour(s), period is 1 day

                IY=I_YEAR
                IM=I_MONTH
                ID=I_DAY
                IF(I_HOUR <  FIXHD(24,FILE)) HOURS=HOURS+24

              ELSE IF(PERIOD == 720)THEN
! Ancillary data interval in day(s) or hours , period is 1 month

                IY=I_YEAR
                IM=I_MONTH
                ID=FIXHD(23,FILE)
                IF((I_DAY*24+I_HOUR) <                                  &
     &             (FIXHD(23,FILE)*24+FIXHD(24,FILE)))                  &
     &           HOURS=HOURS+720

              ELSE IF(PERIOD == 8640)THEN
! Ancillary data interval in month(s)or days or hours, period is 1 year

                IY=I_YEAR
                IM=FIXHD(22,FILE)
                ID=FIXHD(23,FILE)
                IF((I_MONTH*720+I_DAY*24+I_HOUR) <                      &
     &          (FIXHD(22,FILE)*720+FIXHD(23,FILE)*24+FIXHD(24,FILE)))  &
     &           HOURS=HOURS+8640

              END IF

! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,ID,FIXHD(24,FILE),                    &
     &                     FIXHD(25,FILE),FIXHD(26,FILE),               &
     &                     ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,       &
     &                     LCAL360)
              HOURS=HOURS-SEC/3600-DAY*24

! Do not interpolate in time if data time exactly matches model time

              IF(MOD(HOURS,INTERVAL) == 0) THEN
                LINTERPOLATE=.FALSE.
              END IF
              STEP=HOURS/INTERVAL
              TIME=REAL(HOURS)
              TIME1=STEP*INTERVAL
              TIME2=(STEP+1)*INTERVAL

            ELSE  ! non regular case

!L 2.2.2 Gregorian calender,and data interval is in months,
!L       period is 1 year
!L       Updating interval and number of data times to be skipped
!L       calculated in months.

              TIME=REAL(SEC)/3600+REAL(DAY*24)
              INTERVAL=FIXHD(36,FILE)+FIXHD(35,FILE)*12
              PERIOD=INTHD(3,FILE)*INTERVAL
              IF(PERIOD /= 12)THEN
                ICODE=600+FIELD
           CMESSAGE='REPLANCA: Non-standard period for periodic data'
                RETURN
              ENDIF
!  Difference between date now (month) & first date ancil file (month)
              MONTHS=I_MONTH-FIXHD(22,FILE)
              ! required to avoid tripping -check uninit when
              ! writing out hours below
              HOURS=INT(TIME)

           IF (LAMIP2X) THEN ! correct code to use lookup header dates
! Correctly use day and hour from lookup header not fixhd which
! contains values for first field on ancillary file only.
             step=months/INTERVAL
             I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*step
             I1=I2+LOOKUP_START(FILE)-1
!  Check for time within month - using ppheader information
             IF((I_DAY*24+I_HOUR) <  (lookup(3,i1)*24+lookup(4,i1))) THEN
                 MONTHS=MONTHS-1
             END IF
             IF(MONTHS <  0) THEN
                MONTHS=MONTHS+12
             END IF
! recalculate STEP
             STEP=MONTHS/INTERVAL
! NB INTERVAL may be > 1 month
             MONTHS=STEP*INTERVAL
             IY=I_YEAR
             IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
             IF(IM >  I_MONTH) IY=IY-1
             I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP
             I1=I2+LOOKUP_START(FILE)-1
! DEPENDS ON: time2sec
             CALL TIME2SEC(IY,IM,lookup(3,i1),lookup(4,i1),             &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)
             TIME1=REAL(SEC)/3600+REAL(DAY*24)
!  Calculate  TIME2 for second ancillary data time
!  set IY correctly for time interpolation calculations
             IY=I_YEAR
             IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
             IF(IM <  I_MONTH) IY=IY+1
             I1=(IM-1)/INTERVAL
             I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*I1
             I1=I2+LOOKUP_START(FILE)-1
! DEPENDS ON: time2sec
             CALL TIME2SEC(IY,IM,lookup(3,i1),lookup(4,i1),             &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,LCAL360)
             TIME2=REAL(SEC)/3600+REAL(DAY*24)

           ELSE   ! original code inaccurate use of FIXHD dates
!  Check for time within month
              IF((I_DAY*24+I_HOUR) <                                    &
     &           (FIXHD(23,FILE)*24+FIXHD(24,FILE))) THEN
                MONTHS=MONTHS-1
              END IF
              IF(MONTHS <  0) THEN
                MONTHS=MONTHS+12
              END IF

              STEP=MONTHS/INTERVAL
! NB INTERVAL may be > 1 month
              MONTHS=STEP*INTERVAL
!  Calculate TIME1 for first ancillary data time
!  set IY correctly for time interpolation calculations
              IY=I_YEAR
              IM=MOD(FIXHD(22,FILE)+MONTHS-1,12)+1
              IF(IM >  I_MONTH) IY=IY-1
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),        &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME1=REAL(SEC)/3600+REAL(DAY*24)
!  Calculate  TIME2 for second ancillary data time
!  set IY correctly for time interpolation calculations
              IY=I_YEAR
              IM=MOD(FIXHD(22,FILE)+MONTHS+INTERVAL-1,12)+1
              IF(IM <  I_MONTH) IY=IY+1
! DEPENDS ON: time2sec
              CALL TIME2SEC(IY,IM,FIXHD(23,FILE),FIXHD(24,FILE),        &
     &              FIXHD(25,FILE),FIXHD(26,FILE),                      &
     &              ANCIL_REF_DAYS,ANCIL_REF_SECS,DAY,SEC,              &
     &              LCAL360)
              TIME2=REAL(SEC)/3600+REAL(DAY*24)
           ENDIF  ! end LAMIP2X test

! Do not interpolate in time if data time exactly matches model time

              IF(TIME == TIME1) THEN
                LINTERPOLATE=.FALSE.
              END IF

            ENDIF  ! regular/non-regular


!L Adjust YMD time back to start of updating interval

            I_YEAR=I_YEAR1
            I_MONTH=I_MONTH1
            I_DAY=I_DAY1
            I_HOUR=I_HOUR1


          ENDIF  ! non-periodic/periodic

        IF (LINTERPOLATE) THEN
        WRITE(6,*)' REPLANCA - time interpolation for field ',field
        WRITE(6,*)' time,time1,time2 ',time,time1,time2
        WRITE(6,*)' hours,int,period ',hours,interval,period
        END IF

        END IF ! singletime/non-singletime

!L 2.3   Check STASH Code

        I2=NLOOKUP(FIELD)+LOOKUP_STEP(FIELD)*STEP

        I1=LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)-1)

        LMISMATCH=.FALSE.
        WRITE(6,*)' Information used in checking ancillary data set:',  &
     &  ' position of lookup table in dataset:',I2
        WRITE(6,*)' Position of first lookup table referring to ',      &
     &  'data type ',NLOOKUP(FIELD)
        WRITE(6,*)' Interval between lookup tables referring to data ', &
     &  'type ', LOOKUP_STEP(FIELD),' Number of steps', STEP
        WRITE(6,*)' STASH code in dataset ',I1,                         &
     &  '  STASH code requested ',STASHANCIL(FIELD)
        WRITE(6,*)'''Start'' position of lookup tables for dataset ',   &
     &  'in overall lookup array ' ,LOOKUP_START(FILE)

        IF(I1 /= STASHANCIL(FIELD)) THEN
        WRITE(6,*)I1,STASHANCIL(FIELD),FIELD
          LMISMATCH=.TRUE.
        END IF

!L Error exit if checks fail

        IF(LMISMATCH) THEN
          ICODE=200+FIELD
         CMESSAGE='REPLANCA: PP HEADERS ON ANCILLARY FILE DO NOT MATCH'
         RETURN
        END IF

        IF(LINTERPOLATE.AND..NOT.SINGLE_TIME) THEN
!L Check time interpolation factors
          IF(TIME <  TIME1.OR.TIME >  TIME2) THEN
           WRITE(6,*)' Information used in interpolation/replacement:'
           WRITE(6,*)' Time of first data=', TIME1
           WRITE(6,*)' Validity Time for update=', TIME
           WRITE(6,*)' Time of second data=', TIME2

           ICODE=500+FIELD
           CMESSAGE='REPLANCA: TIME INTERPOLATION ERROR'
           RETURN
          END IF
        END IF

!L 3   Loop over levels of ancillary data for field I
!L Reset pointer for dataset


!L Includes loop over X and Y components of surface currents

         LICE_FRACTION=FIELD == 27
         LSNOW_DEPTH=FIELD == 9
         LICE_DEPTH=FIELD == 29

        DO 30 LEVEL=1,LEVELS(FIELD)

!L Do not go through loop for ice edge or snow edge

        IF((LICE_FRACTION.OR.LSNOW_DEPTH).AND.LEVEL == 2) THEN
          GOTO 30
        END IF

!L 3.1 Read data for single level of ancillary field.

        IF(.NOT.LICE_FRACTION) THEN
! AMIPII case ice depth field not read from ancillary file
         IF(.NOT.(LICE_DEPTH.and.LAMIP2X)) THEN
! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,1,I2,LOOKUP(1,LOOKUP_START(FILE)),        &
     &                  LEN1_LOOKUP,ANCIL1,P_FIELD,FIXHD(1,FILE),       &
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
     &                  ICODE,CMESSAGE)

         ENDIF
          IF(ICODE /= 0)THEN
            ICODE=FIELD+100
            IOUNIT=NFTIN
           CMESSAGE='REPLANCA :I/O ERROR '
            RETURN
          END IF

        ELSE

!L If ice-fraction,read fractional time field as well
!L       UNLESS IT IS A SINGLE TIME FIELD
!L If snow-depth,read fractional time field as well only if time
!L interpolation required.

      IF(.NOT.SINGLE_TIME.and..NOT.LAMIP2X) THEN
         IF(LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)) == 38) THEN
! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,2,I2,LOOKUP(1,LOOKUP_START(FILE)),        &
     &                  LEN1_LOOKUP,ICE_EXTENT,P_FIELD,FIXHD(1,FILE),   &
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
     &                  ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
            ICODE=FIELD+100
            IOUNIT=NFTIN
           CMESSAGE='REPLANCA :I/O ERROR '
            RETURN
          END IF

         ELSE
           ICODE=FIELD+100
           IOUNIT=NFTIN
            CMESSAGE='REPLANCA :ICE CHANGE DATA MISSING'
            RETURN
         END IF
        ELSE    ! single time or LAMIP2X - ie no time change field
! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,1,I2,LOOKUP(1,LOOKUP_START(FILE)),        &
     &                  LEN1_LOOKUP,ICE_EXTENT,P_FIELD,FIXHD(1,FILE),   &
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
     &                  ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
            ICODE=FIELD+100
            IOUNIT=NFTIN
            CMESSAGE='REPLANCA: I/O ERROR'
            RETURN
          ENDIF
        END IF
      ENDIF

        IF(LSNOW_DEPTH.AND.LINTERPOLATE) THEN
      IF(LOOKUP(ITEM_CODE,I2+LOOKUP_START(FILE)) == 27) THEN

! DEPENDS ON: readflds
           CALL READFLDS(NFTIN,1,I2+1,LOOKUP(1,LOOKUP_START(FILE)),     &
     &                   LEN1_LOOKUP,SNOW_CHANGE,P_FIELD,FIXHD(1,FILE), &
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
     &                   ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
             ICODE=FIELD+100
             IOUNIT=NFTIN
            CMESSAGE='REPLANCA :I/O ERROR '
             RETURN
           END IF

         ELSE
           ICODE=FIELD+100
           IOUNIT=NFTIN
           CMESSAGE='REPLANCA :SNOW CHANGE DATA MISSING'
           RETURN
         END IF
        END IF

!L If sea surface temperature or other ice fields, read ice fraction
!L and fractional time field if not already pressent and if required
!L by time interpolation.  Similar if SLAB ref SST or ice depth needed.

        IF(FIELD == 29.OR.(FIELD == 28.AND.LT_INT_C).OR.                &
     &     FIELD == 38)                                                 &
     &    THEN

         IF(.NOT.UPDATE(27)) THEN
          I3 = NLOOKUP(27) + LOOKUP_STEP(27)*STEP + LOOKUP_START(       &
     &       FILEANCIL(27))
          IF ( LOOKUP(ITEM_CODE,I3)  ==  38 ) THEN

! DEPENDS ON: readflds
            CALL READFLDS(FTNANCIL(FILEANCIL(27)),2,                    &
     &                    NLOOKUP(27)+LOOKUP_STEP(27)*STEP,             &
     &                    LOOKUP(1,LOOKUP_START(FILEANCIL(27))),        &
     &                    LEN1_LOOKUP,ICE_EXTENT,                       &
     &                    P_FIELD,FIXHD(1,FILEANCIL(27)),               &
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
     &                    ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
              ICODE=FIELD+100
              IOUNIT=NFTIN
             CMESSAGE='REPLANCA :I/O ERROR '
              RETURN
            END IF
          IF ( RLOOKUP(BMDI,I3-1)  /=  RMDI ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti&
     &mes of sea-ice chge not standard'
            RETURN
          ENDIF


          ELSE
            ICODE=FIELD+100
            IOUNIT=NFTIN
            CMESSAGE='REPLANCA :ICE FIELD DATA MISSING'
            RETURN
          END IF
         END IF
        END IF

!L 3.3 If time interpolation required, read second record

        IF(LINTERPOLATE) THEN

          I1=I2+ LOOKUP_STEP(FIELD)
          IF(I1 <= FIXHD(152,FILE)) THEN

! AMIP II and ice depth don't read in ice depth field
          IF (.NOT.(LAMIP2X.and.LICE_DEPTH)) THEN

! DEPENDS ON: readflds
            CALL READFLDS(NFTIN,1,I1,LOOKUP(1,LOOKUP_START(FILE)),      &
     &                    LEN1_LOOKUP,ANCIL2,P_FIELD,FIXHD(1,FILE),     &
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
     &                    ICODE,CMESSAGE)
          ENDIF
          IF(ICODE /= 0)THEN
              ICODE=FIELD+300
              IOUNIT=NFTIN
              CMESSAGE='REPLANCA :I/O ERROR '
              RETURN
            END IF

          ELSE !end of data on file

!L  If end of data has been reached go back to the start.If data is
!L  periodic.
!L  Otherwise cancel time interpolation

            IF(PERIODIC) THEN

              I1 = NLOOKUP(FIELD) + LEVEL - 1

! DEPENDS ON: readflds
              CALL READFLDS(NFTIN,1,I1,LOOKUP(1,LOOKUP_START(FILE)),    &
     &                      LEN1_LOOKUP,ANCIL2,P_FIELD,FIXHD(1,FILE),   &
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
     &                      ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
                ICODE=FIELD+300
                IOUNIT=NFTIN
               CMESSAGE='REPLANCA :I/O ERROR '
                RETURN
              END IF
            ELSE
              LINTERPOLATE=.FALSE.
            END IF
          END IF! End of position on file test

          ICODE=0
        END IF ! End LINTERPOLATE

!L 3.4 Perform time interpolation

        IF(LINTERPOLATE) THEN

          ZERO=0.0

!L Select appropriate time interpolation for each field
!  Snowdepth: set equal to zero if no snow cover

          IF(LSNOW_DEPTH) THEN
            DO I=1,P_FIELD
              PRES_VALUE(I)=ZERO
            END DO

! For the call to T_INT_C, need to know BMDI is OK for SNOW_CHANGE
!  which was read in from position I2+1.
          IF ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I2)  /=  RMDI ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti&
     &mes of snow change non-standard '
            RETURN
          ENDIF

! DEPENDS ON: t_int_c
            CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,         &
     &           TIME,P_FIELD,SNOW_CHANGE,ANCIL1,PRES_VALUE)

! Ice fraction: ice depth set equal to zero if no ice

          ELSE IF(FIELD == 27.OR.FIELD == 29.OR.FIELD == 38) THEN
            IF(FIELD == 27) THEN
! For the call to T_INT_C, need to know BMDI is OK for ICE_EXTENT(1,2)
!  which was read in from position I1+1
          IF(.NOT.LAMIP2X) THEN
          IF ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I1)  /=  RMDI ) THEN
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: RMDI in lookup of ancillary file of ti&
     &mes of sea-ice chge non-standard'
            RETURN
          ENDIF
          ENDIF

             IF (LAMIP2X) THEN
! linear uncontrolled time interpolation
! DEPENDS ON: t_int
              CALL T_INT (ICE_EXTENT,TIME1,ANCIL2,TIME2,ANCIL_DATA,     &
     &             TIME,P_FIELD)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.

              DO I=1,P_FIELD
                IF (ANCIL_DATA(I) <  0.3) ANCIL_DATA(I)=0.0
                IF (ANCIL_DATA(I) >  1.0) ANCIL_DATA(I)=1.0
              ENDDO

             ELSE       ! non AMIPII option
              DO I=1,P_FIELD
                PRES_VALUE(I)=0
              END DO

! DEPENDS ON: t_int_c
              CALL T_INT_C (ICE_EXTENT,TIME1,ANCIL2,TIME2,ANCIL_DATA,   &
     &             TIME,P_FIELD,ICE_EXTENT(1,2),ICE_EXTENT,PRES_VALUE)

             ENDIF     ! end AMIPII test

            ELSE IF (FIELD == 29.OR.FIELD == 38) THEN

              DO I=1,P_FIELD
                PRES_VALUE(I)=0
              END DO

! DEPENDS ON: t_int_c
              CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,       &
     &             TIME,P_FIELD,ICE_EXTENT(1,2),ICE_EXTENT,PRES_VALUE)


            END IF


! Sea surface temperature, set equal to TFS if ice present

          ELSE IF (FIELD == 28.AND.LT_INT_C) THEN
           IF (LAMIP2X) THEN

! DEPENDS ON: t_int
            CALL T_INT (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,           &
     &              TIME,P_FIELD)
! remove any T below TFS
            DO I=1,P_FIELD
              IF (ANCIL_DATA(i) <  TFS)  ANCIL_DATA(I)=TFS
            ENDDO

           ELSE     ! non AMIPII option

            IF(.NOT.LTLEADS)THEN
            DO I=1,P_FIELD
                PRES_VALUE(I)=TFS

! Set no_ice_extent indicator for controlled SST interpolation
                IF(ICE_EXTENT(I,1) == 0) THEN
                  NO_ICE_EXTENT(I)=1.0
                ELSE
                  NO_ICE_EXTENT(I)=0.0
                ENDIF
            END DO

! DEPENDS ON: t_int_c
            CALL T_INT_C (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,         &
     &           TIME,P_FIELD,ICE_EXTENT(1,2),NO_ICE_EXTENT,PRES_VALUE)
            ELSE
! DEPENDS ON: t_int
            CALL T_INT (ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,           &
     &              TIME,P_FIELD)
            ENDIF

           ENDIF   ! end AMIPII test
! Otherwise linear interpolation in time, unless missing data indicator
! present at either time.

          ELSE

! Time interpolation checks the data against the standard missing data
!   indicator - check that the field is labelled as using the same one.
!  (It is to have the right I1 here that I3 is used above.)
          IF ( RLOOKUP(BMDI,LOOKUP_START(FILE)+I1-1)  /=  RMDI .OR.     &
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I2-1)  /=  RMDI ) THEN
            WRITE (6, *) 'LOOKUPS:',                                    &
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I1-1),                   &
     &         RLOOKUP(BMDI,LOOKUP_START(FILE)+I2-1)
            ICODE = 700 + FIELD
            CMESSAGE = 'REPLANCA: Missing data indicator in lookup of an&
     &cillary file is non-standard    '
            RETURN
          ENDIF

          LEN=P_FIELD
!L  Ozone, test for zonal mean or full field
          IF(FIELD == 7) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
              LEN=P_ROWS
            END IF
!   Cariolle ozone, test for zonal mean or full field. 
!   Currently same test as for conventional ozone.
          ELSE IF (FIELD >= 178.AND.FIELD <= 185) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
              LEN=P_ROWS
            END IF
          END IF

! DEPENDS ON: t_int
            CALL T_INT(ANCIL1,TIME1,ANCIL2,TIME2,ANCIL_DATA,            &
     &                 TIME,LEN)

          END IF ! End Lsnow_depth

! If no interpolation, copy data into final array

        ELSE ! no interpolation
         IF(LICE_FRACTION) THEN
          IF (LAMIP2X) THEN
          DO I=1,P_FIELD

          ANCIL_DATA(I)=ICE_EXTENT(I,1)

! For AMIP II strictly ice concentrations should range between
! 0.0 and 1.0 but because of assumptions on T* made by the boundary
! layer and radiation schemes ice concentrations are restricted to
! 0.3 to 1.0. This will allow SSTs in areas of less than 30% ice to
! be used rather than TFS=-1.8C.

             IF (ANCIL_DATA(I) <  0.3) ANCIL_DATA(I)=0.0
             IF (ANCIL_DATA(I) >  1.0) ANCIL_DATA(I)=1.0

          ENDDO
          ELSE           ! non AMIP II option
            DO I=1,P_FIELD
             ANCIL_DATA(I)=ICE_EXTENT(I,1)
            ENDDO
          ENDIF           ! end of AMIPII test
         ELSE IF (LAMIP2X.AND.FIELD == 28) THEN
          DO I=1,P_FIELD
            ANCIL_DATA(I)=ANCIL1(I)
            IF (ANCIL_DATA(I) <  TFS) ANCIL_DATA(I)=TFS
          ENDDO
         ELSE
          DO I=1,P_FIELD
            ANCIL_DATA(I)=ANCIL1(I)
          END DO
         ENDIF
        END IF !End interpolate/no interpolate

!L 3.5 Updating action for each field at each level
!L     Fields replaced except that Sea Surface Temperature may be
!L     incremented. Take apropriate action for each field.

        IF(FIELD <= 2.OR.FIELD == 7.OR.FIELD == 39.OR.FIELD == 40       &
     &  .OR.FIELD == 41.OR.FIELD == 42.OR.FIELD == 43                   &
     &  .OR.FIELD == 44.OR.FIELD == 45                                  &
                                          ! multi-level murk
     &  .OR.(FIELD >= 48 .AND. FIELD <= 67 .AND. L_UKCA)                &
                                          ! single-level user ancillaries
     &  .OR.(FIELD >= 68.AND.FIELD <= 70)                               &
                                          !NH3,soot aerosol emissions
     &  .OR.(FIELD >= 72.AND.FIELD <= 77)                               &
                                          !Sulphur cycle
     &  .OR.FIELD == 78                                                 &
                                          !CO2 EMISSIONS
     &  .OR.FIELD == 82                                                 &
                                          !HADCM2 sulphate aerosol
     &  .OR.(FIELD >= 90.AND.FIELD <= 109)                              &
                                           !multi-level user ancillaries
     &  .OR.(FIELD >= 112.AND.FIELD <= 120)                             &
                                            !mineral dust fields
     &  .OR.(FIELD >= 121.AND.FIELD <= 122)                             &
                                            !Biomass emissions
     &  .OR.FIELD == 123                                                &
                                           !Seawater DMS concentration
     &  .OR.(FIELD >= 157.AND.FIELD <= 177)                             &
                                           !Aerosol climatologies
     &  .OR.(FIELD >=178 .AND.FIELD<=185)                               &
                                           !Cariolle ozone ancillaries
     &  .OR.(FIELD >= 186.AND.FIELD <= 187)                             &
                                           !OCFF emissions
     &  )THEN

!L 3.5.0 Updates at all points

          LEN=P_FIELD
!L  Ozone, test for zonal mean or full field
          IF(FIELD == 7) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
              LEN=P_ROWS
            END IF

!   Cariolle ozone, test for zonal mean or full field. 
!   Currently same test as for conventional ozone.
          ELSE IF (FIELD >= 178.AND.FIELD <= 185) THEN
            IF(LOOKUP(LBNPT,LOOKUP_START(FILE)+I2-1) == 1) THEN
              LEN=P_ROWS
            END IF
          END IF

          DO I=1,LEN
            D1(D1_ANCILADD(FIELD)+I-1+(LEVEL-1)*LEN)=ANCIL_DATA(I)
          END DO

!L 3.5.1 Updates over all land points

        ELSEIF((FIELD > 2.AND.FIELD < 7)                                &
     &   .OR.(FIELD > 7.AND.FIELD < 27)                                 &
     &   .OR.(FIELD == 32).OR.(FIELD >= 34.AND.FIELD <= 36)             &
     &   .OR.(FIELD >= 48 .AND. FIELD <= 67 .AND. .NOT. L_UKCA)         &
                                      ! single level user ancillaries
     &   .OR.(FIELD >= 46.AND.FIELD <= 47)                              &
                                                 !Orographic roughness
     &   .OR.(FIELD >= 155.AND.FIELD <= 156)                            &
                                      ! Orographic X & Y gradients
     &   .OR.(FIELD >= 79.AND.FIELD <= 81)                              &
                                                 !MOSES-I
     &   .OR.(FIELD >= 83.AND.FIELD <= 89)) THEN !MOSES-II



!L If not reconfiguration, set snowdepth values at all land points
!L Reset TSTAR to TM if snow cover present

          IF(LSNOW_DEPTH) THEN
            DO I=1,P_FIELD
              IF(LAND(I)) THEN
                D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
                IF(TSTAR_LAND(I) >  TM.AND.ANCIL_DATA(I) >  0.0) THEN
                  TSTAR_LAND(I)=TM
                END IF
              END IF
            END DO

!L Set all other fields , which are stored at land points only

          ELSE
!L If field is a single time soil moisture, only update if
!L model time matches data time and then deactivate to prevent
!L any further updating

            IF(FIELD == 36.AND.FIXHD(10,FILE) == 0) THEN
              IF(LOOKUP(LBYR,LOOKUP_START(FILE)) == I_YEAR.AND.         &
     &           LOOKUP(LBMON,LOOKUP_START(FILE)) == I_MONTH.AND.       &
     &           LOOKUP(LBDAT,LOOKUP_START(FILE)) == I_DAY.AND.         &
     &           LOOKUP(LBHR,LOOKUP_START(FILE)) == I_HOUR) THEN
                WRITE(6,*) 'Updating soil moisture at ',                &
     &                       I_YEAR,I_MONTH,I_DAY,I_HOUR,               &
     &          ' for level ',RLOOKUP(BLEV,LOOKUP_START(FILE)+LEVEL-1)
                DZ_SOIL(LEVEL)=RLOOKUP(BLEV,LOOKUP_START(FILE)+LEVEL-1)
! DEPENDS ON: to_land_points
                CALL TO_LAND_POINTS(ANCIL_DATA,D1(D1_ANCILADD(FIELD)+   &
     &                             (LEVEL-1)*LAND_FIELD),LAND,P_FIELD,I)
! Switch off to prevent further attempts to update
                FIELDCODE(1,FIELD)=0
                STEPS(FIELD)=0
! Set flag to indicate that soil moisture has been updated
                SMC_UPDATED=.TRUE.
              ELSE
                WRITE(6,*) 'Update of soil moisture skipped'
              END IF

            ELSE
! other fields
! DEPENDS ON: to_land_points
              CALL TO_LAND_POINTS(ANCIL_DATA,D1(D1_ANCILADD(FIELD)+     &
     &                           (LEVEL-1)*LAND_FIELD),LAND,P_FIELD,I)
            END IF

          END IF



!L 3.5.2 Ice fraction


        ELSE IF(FIELD == 27) THEN
          DO I=1,P_FIELD
            ICE_FRACTION(I)=0.
            IF (SEA(I)) THEN
              ICE_FRACTION(I)=ANCIL_DATA(I)
            END IF
          END DO

!L Reduce TSTAR to TFS where ice fraction greater than zero
! Required at present because radiation and boundary layer codes
! assume T* is TFS and ignore any value set in TSTAR.

          IF(.NOT.LTLEADS)THEN
            DO I=1,P_FIELD
              IF(ICE_FRACTION(I) >  0.0) THEN
                TSTAR_SSI(I)=AMIN1(TSTAR_SSI(I),TFS)
              ENDIF
            END DO
          ENDIF

!L 3.5.3 Sea surface temperatures for atmosphere, allow fields to be
!L       incremented rather than replaced

        ELSE IF (FIELD == 28) THEN


          DO I=1,P_FIELD
            IF(L_CTILE.OR.ICE_FRACTION(I) == 0.0)THEN
              IF (SEA(I)) THEN
                IF (L_SSTANOM) THEN
                TSTAR_SEA(I)=ANCIL_DATA(I)+TSTAR_ANOM(I)
                ELSE
                  TSTAR_SEA(I)=ANCIL_DATA(I)
                END IF
                IF(ICE_FRACTION(I) == 0.0)TSTAR_SSI(I)=TSTAR_SEA(I)
              END IF
            END IF
          END DO

!L 3.5.3.1 Reference SSTs for SLAB model

        ELSE IF (FIELD == 37) THEN

          DO I=1,P_FIELD
            IF (SEA(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            ELSE
              D1(D1_ANCILADD(FIELD)+I-1)=RMDI

            ENDIF
          END DO

!L 3.5.4 Sea ice thickness/Reference seaice thickness for SLAB
!L       Update over all sea points (all sea ice points are the only
!L       ones strictly required, but this cannot be determined easily)

        ELSE IF (FIELD == 29.OR.FIELD == 38) THEN

          DO I=1,P_FIELD
            IF(SEA(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            END IF
          END DO

!L 3.5.5 Surface currents

        ELSE IF (FIELD == 30) THEN
          DO I=1,U_FIELD
            D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
          END DO

        ELSE IF (FIELD == 31) THEN
          DO I=1,V_FIELD
            D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
          END DO
!L 3.5.6 Heat convergence (slab model)
!L       Update over all non-land points

        ELSE IF (FIELD == 33) THEN

          DO I=1,P_FIELD
            IF(SEA(I)) THEN
              D1(D1_ANCILADD(FIELD)+I-1)=ANCIL_DATA(I)
            ELSE
              D1(D1_ANCILADD(FIELD)+I-1) = 0.0
            END IF
          END DO

        ELSE IF (FIELD >= 124.AND.FIELD <= 126) THEN   ! Riv Routing
          ICODE=750
          CMESSAGE='REPLANCA: ERROR Trying to use Riv Route Ancils'
!         There is no code yet to support in model updateing of the
!         River Routing fields.
          RETURN

! Tracer Fluxes - kdcorbin, 05/10
        ELSE IF (FIELD .ge. 188 .and. FIELD .lt. 208) Then
            DO I=1,P_FIELD
              D1(D1_ANCILADD(FIELD)+I-1) = ANCIL_DATA(I)
            ENDDO

        ELSE

        WRITE(6,*)' REPLANCA: ERROR - FIELD ',FIELD,                    &
     &  ' omitted from update block'

        END IF !End tests on FIELD numbers

!L End loop over levels

      I2=I2+1

 30   CONTINUE

!L End loop over ancillary fields (atmosphere)
       ENDIF ! LAMIP2X and ice depth test

      END IF    ! End UPDATE(field) test     level 1 IF


      END DO


      IF(L_CTILE)THEN
        DO I=1,P_FIELD
          IF(SEA(I).AND.ICE_FRACTION(I) >  0.0)THEN
            IF(LTLEADS.OR.LAMIP2X)THEN

              TSTAR_SSI(I)=ICE_FRACTION(I)*TSTAR_SICE(I)                &
     &          +(1.-ICE_FRACTION(I))*TSTAR_SEA(I)

            ELSE

              TSTAR_SEA(I)=TFS
              TSTAR_SICE(I)=(TSTAR_SSI(I)                               &
     &          -(1.-ICE_FRACTION(I))*TSTAR_SEA(I))/ICE_FRACTION(I)

            ENDIF
          ENDIF
!
          TSTAR(I)=FLAND_G(I)*TSTAR_LAND(I)                             &
     &      +(1.-FLAND_G(I))*TSTAR_SSI(I)
        ENDDO
      ELSE
        DO I=1,P_FIELD
          IF(LAND(I))THEN
            TSTAR(I)=TSTAR_LAND(I)
          ELSE
            TSTAR(I)=TSTAR_SSI(I)
          ENDIF
        ENDDO
      ENDIF

!     Set up surface temperatures:
      IF(L_CTILE)THEN
        DO I=1,P_FIELD
          TSTAR_LAND_CTILE(I)=TSTAR_LAND(I)
          TSTAR_SEA_CTILE(I)=TSTAR_SEA(I)
          ! The use of TSTAR_SICE appears to cause problems in 
          ! some configurations (e.g. seasonal). Possibly because
          ! of the use of inconsistent ancillary fields.
          ! Hence this is commented out but retained for reference.
          ! [See also equivalent change in replanca-rcf_replanca.F90]
          ! TSTAR_SICE_CTILE(I)=TSTAR_SICE(I)
        ENDDO
      ENDIF

!
900   RETURN
      END SUBROUTINE REPLANCA

