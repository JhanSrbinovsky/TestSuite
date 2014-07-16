
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Construct preliminary STASH list of user requests
!
! Subroutine Interface:

      SUBROUTINE PRELIM(NRECS,                                          &
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
     &                 NTIMES,NLEVELS,ErrorStatus,CMESSAGE)
      IMPLICIT NONE

!  Description:
!  Constructs a preliminary STASH list of user requests. Uses interim
!  pointer system, by means of the "extra entry" NELEMP+1 in the LIST_S
!  array. At this stage, the input levels encompass all possible levels.
!  Called by STPROC.
!
!  Method:
!
!  Current code owner:  UM System Team
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!   4.1     Apr. 96    Numerous improvements associated with wave model,
!                       correction of output-times table processing,
!                       comprehensive soft-abort system, etc.
!                                             S.J.Swarbrick
!   4.4     Sep. 97    Allow offset for sampling frequency
!                      S.D. Mullerworth
!LL 4.4    21/11/96   Allow daily mean timeseries. R.A.Stratton
!   4.4     Oct. 97    Added checking of error returns from TOTIMP.
!                         Shaun de Witt
!   4.5    18/11/98    Allow new sampling frequencies for vegetation.
!                      Richard Betts
!   5.2  15/11/00      Correct st_end_time_code when using time
!                      processing                    P.Selwood
!   5.0    22/11/99    Correct STASH input code for fields in
!                      secondary space - first step to enabling
!                      STASH output of these fields. R. Rawlins
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!   5.1    15/05/00    S-N ordering consistency correction. R Rawlins
!   5.3    23/10/01    Minor enhancement to error reporting. R Rawlins
!   5.5    17/02/03    Upgrade Wave model from 4.1 to 5.5. D.Holmes-Bell
!   6.0    18/08/03    Check that requested diagnostic exists in
!                      STASHmaster.  E. Leung
!   6.0  15/12/03  Update IOPN inline with
!                  stashmaster 30 digit option codes.  M.Hughes
!   6.2  23/11/05  Removed all references to the wavemodel.
!                  T.Edwards
!   6.2   16/11/05     Fix offsets for radiation/convection diagnostics.
!                      Improve warning reporting. P.Selwood
!
! 6.2  25/11/05 Functionality for improved time stepping, radiative
!               forcing and radiance code added for versions 3C
!               and 3Z of radiation code             (J.-C. Thelen)
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!
!  Global variables:

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
! TYPSIZE start
!   Description:
!     This file contains sizes needed for dynamic allocation of
!   main data arrays within the model. Sizes read in from the user
!   interface via NAMELISTs are passed by /COMMON/. Other control
!   sizes that are fundamental in the definition of data structures
!   are assigned by PARAMETER statements.
!
      ! All sizes
      ! Not dependent on sub-model
      ! DATA IN NAMLST#x MEMBER OF THE JOB LIBRARY
      ! atmos START
      ! Main sizes of fields for each submodel
      ! Grid-related sizes for ATMOSPHERE submodel.
      INTEGER:: ROW_LENGTH           ! IN: No of points per local row
      INTEGER:: global_ROW_LENGTH    ! IN: Points per global row
      INTEGER:: ROWS                 ! IN: No of local (theta) rows
      INTEGER:: global_ROWS          ! IN: No of global (theta) rows
      INTEGER:: MODEL_LEVELS         ! IN: No of model levels
      INTEGER:: LAND_FIELD           ! IN: No of land points in field
      INTEGER:: NTILES               ! IN: No of land surface tiles

      ! Physics-related sizes for ATMOSPHERE submodel
      INTEGER:: WET_LEVELS          ! IN: No of moist-levels
      INTEGER:: CLOUD_LEVELS        ! IN: No of cloud-levels
      INTEGER:: ST_LEVELS           ! IN: No of soil temperature levels
      INTEGER:: SM_LEVELS           ! IN: No of soil moisture levels
      INTEGER:: BL_LEVELS           ! IN: No of boundary-layer-levels
      INTEGER :: OZONE_LEVELS       ! IN: No of ozone-levels
      INTEGER :: TPPS_OZONE_LEVELS  ! IN: No of tropopause-ozone-levels
      INTEGER :: RIVER_ROWS         ! IN: No of rows for river routing
      INTEGER :: RIVER_ROW_LENGTH   ! IN: Row length for river routing
      ! Dynamics-related sizes for ATMOSPHERE submodel

      INTEGER:: TR_LEVELS            ! IN: No of tracer-levels
      INTEGER:: TR_VARS              ! IN: No of passive tracers
      INTEGER:: TR_UKCA              ! IN: No of UKCA tracers

      ! Dynamics output diagnostic-related sizes for ATMOSPHERE submodel
      INTEGER:: THETA_PV_P_LEVS   ! IN: No of levels requested for pvort

      ! For Small executables
      INTEGER:: TOT_LEVELS
      ! Assimilation-related sizes for ATMOSPHERE submodel
      INTEGER :: N_AOBS           ! IN: No. of atmos observation types

      ! Grid related sizes for data structure
      ! Data structure sizes for ATMOSPHERE submodel
      INTEGER:: A_PROG_LOOKUP     ! IN: No of prognostic fields
      INTEGER:: A_PROG_LEN        ! IN: Total length of prog fields
      INTEGER:: A_LEN_INTHD       ! IN: Length of INTEGER header
      INTEGER:: A_LEN_REALHD      ! IN: Length of REAL header
      INTEGER:: A_LEN2_LEVDEPC    ! IN: No of LEVEL-dependent arrays
      INTEGER:: A_LEN2_ROWDEPC    ! IN: No of ROW-dependent arrays
      INTEGER:: A_LEN2_COLDEPC    ! IN: No of COLUMN-dependent arrays
      INTEGER:: A_LEN2_FLDDEPC    ! IN: No of FIELD arrays
      INTEGER:: A_LEN_EXTCNST     ! IN: No of EXTRA scalar constants
      INTEGER:: A_LEN_CFI1        ! IN: Length of compressed fld index 1
      INTEGER:: A_LEN_CFI2        ! IN: Length of compressed fld index 2
      INTEGER:: A_LEN_CFI3        ! IN: Length of compressed fld index 3
      ! atmos end

      ! OCEAN start
! TYPOCPAR start
!  History:
!  Version   Date     Comment
!  -------   ----     -------
!    4.4   15.06.97   Add free surface scalar R.Lenton
!     5.1   07.01.00   Invert_ocean to false with New Dynamics. JC Thil.
!     5.4   29.08.02   Add N_STRAIT and N_STRAIT_CLM. D. Storkey
!     5.5   03.01.03   Remove typocbas. R. Hill
!
      ! Grid related sizes for OCEAN model
      INTEGER::LSEG!IN:Max no of sets of start/end indices for vorticity
      INTEGER:: NISLE            ! IN: No of islands
      INTEGER:: ISEGM            ! IN: Max no of island segments per box
      INTEGER :: N_STRAIT       ! IN: No of pairs of Strait exchange pts
      INTEGER :: N_STRAIT_CLM   ! IN: No of Strait pts set by climate
                                !     values
      INTEGER:: O_LEN_COMPRESSED ! IN: No of ocean points in 3D field
      INTEGER:: LSEGC            ! IN: No of island basins for mead calc
      ! No of start/end indicies for the free surface solution
      INTEGER:: LSEGFS ! IN

      ! Fourier filtering for OCEAN submodel
      INTEGER :: LSEGF    ! IN: max. no of sets of indices for filtering
      INTEGER :: JFRST    ! IN: first J row of T to be filtered

      ! filtering is done on T with a low pass cut off to make the zonal
      ! dimension of the box filtered effectively the same as that of
      ! the boxes on row JFT0
      INTEGER :: JFT0     ! IN:

      INTEGER :: JFT1     ! IN: last J row of T in SH to be filtered
      INTEGER :: JFT2     ! IN: first J row of T in NH to be filtered
      INTEGER :: JFU0     ! IN: same function as JFT0 but for U,V
      INTEGER :: JFU1     ! IN: last J row of U,V in SH to be filtered
      INTEGER :: JFU2     ! IN: first J row of U,V in NH to be filtered

      ! Variables derived from those above
      INTEGER :: IMU     ! IN: total number of U,V grid boxes zonally
      INTEGER :: IMTP1   ! IN: IMT+1
      INTEGER :: IMTM1   ! IN: IMT-1
      INTEGER :: IMTM2   ! IN: IMT-2
      INTEGER :: IMUM1   ! IN: IMU-1
      INTEGER :: IMUM2   ! IN: IMU-2
      INTEGER :: JMTP1   ! IN: JMT+1
      INTEGER :: JMTM1   ! IN: JMT-1
      INTEGER :: JMTM2   ! IN: JMT-2
      INTEGER :: JSCAN   ! IN: JMTM2+1
      INTEGER :: KMP1    ! IN: KM+1
      INTEGER :: KMP2    ! IN: KM+2
      INTEGER :: KMM1    ! IN: KM-1
      INTEGER :: NSLAB   ! IN: no of words in one slab
      INTEGER :: JSKPT   ! IN: no of rows of T and U,V not filtered in
      INTEGER :: JSKPU   ! IN: low and mid latitudes + 1
      INTEGER :: NJTBFT  ! IN: no of J rows to be filtered on T
      INTEGER :: NJTBFU  ! IN: no of J rows to be filtered on U,V
      INTEGER :: IMTKM   ! IN: IMT*KM
      INTEGER :: NTMIN2  ! IN: maximum of NT or 2
      INTEGER :: IMTD2   ! IN: IMT/2
      INTEGER :: LQMSUM  ! IN: IMTD2*(IMT-IMTD2)
      INTEGER :: LHSUM   ! IN: IMT*IMTP1/2
      INTEGER :: IMTX8   ! IN: IMT*8
      INTEGER :: IMTIMT  ! IN: IMT*IMT
      INTEGER :: IMROT   ! X dimension for Coriolis array
      INTEGER :: JMROT   ! Y dimension for Coriolis array
      INTEGER :: IMBC    ! No of columns in boundary field array
      INTEGER :: JMBC    ! No of rows in boundary field array
      INTEGER :: KMBC    ! No of levels in boundary field array
      INTEGER :: NTBC    ! No of tracers in boundary field array
      INTEGER :: JMMD    ! No of rows for mead diagnostic basin indices
      INTEGER :: LDIV    ! No of divisions mead basin indices

      ! Grid-related switches for OCEAN submodel
      LOGICAL :: CYCLIC_OCEAN        ! IN: TRUE if CYCLIC E-W boundary
      LOGICAL :: GLOBAL_OCEAN        ! IN: TRUE if global domain

      ! TRUE if ocean grid NS-inverted cf atmos
      LOGICAL, PARAMETER :: INVERT_OCEAN=.FALSE.
      ! User interface limit for tracers
      ! Max no. tracers in STASHMASTER
      INTEGER, PARAMETER :: O_MAX_TRACERS=20
! COMOCPAR start
      COMMON /COMOCPAR/ GLOBAL_OCEAN, CYCLIC_OCEAN,                     &
     &  LSEG,NISLE,ISEGM,N_STRAIT,N_STRAIT_CLM,O_LEN_COMPRESSED,LSEGC,  &
     &  LSEGFS,LSEGF,JFRST,JFT0,JFT1,JFT2,JFU0,JFU1,JFU2,IMU,IMTP1,     &
     &  IMTM1,IMTM2,IMUM1,IMUM2,JMTP1,JMTM1,JMTM2,JSCAN,KMP1,KMP2,KMM1, &
     &  NSLAB,JSKPT,JSKPU,NJTBFT,NJTBFU,IMTKM,NTMIN2,                   &
     &  IMTD2,LQMSUM,LHSUM,IMTX8,IMTIMT,                                &
     &  IMROT,JMROT,IMBC,JMBC,KMBC,NTBC,JMMD,LDIV
! COMOCPAR end
! TYPOCPAR end
! TYPOCBAS Physics-related sizes for OCEAN submodel
      INTEGER ::  NT ! IN: No of ocean tracers (inc T,S)
      ! Grid related sizes for OCEAN model
      INTEGER :: IMT  ! IN: No of points per row (incl wrap)
      INTEGER :: JMT  ! IN: No of tracer rows
      INTEGER :: KM   ! IN: No of tracer levels
      INTEGER :: NT_UI     ! Copy of NT
      INTEGER :: IMT_UI    ! Copy of IMT
      INTEGER :: JMT_UI    ! Copy of JMT
      INTEGER :: KM_UI     ! Copy of KM
      INTEGER :: NICE      ! IN: No. of sea ice thickness categories
! COMOCBAS start
      COMMON /COMOCBAS/ NT_UI, IMT_UI, JMT_UI, KM_UI                    &
     &        ,NT, IMT, JMT, KM                                         &
     &                 ,NICE
! COMOCBAS end
! TYPOCBAS end
! TYPOASZ sizes for dynamic allocation of ocean assim.
! 5.2 11/08/00  JO_NMAX_OBS_ICE introduced. JO_MAX_OBS_VAL and
!               JO_NMAX_OBS_ICE put into COMOCASZ. M J Bell
      INTEGER :: JO_MAX_OBS_VAL !max number of values in OBS array

      ! max no of flds reqd at once in sea ice assim
      INTEGER :: JO_NMAX_OBS_ICE

      ! length of climate/covariances array
      INTEGER, PARAMETER:: JO_LEN_COV = 1

      ! max number of columns in climate grid
      INTEGER,PARAMETER:: JO_MAX_COLS_C = 1

      ! max number of rows in climate grid
      INTEGER,PARAMETER:: JO_MAX_ROWS_C = 1

      ! max number of levels in climate grid
      INTEGER,PARAMETER:: JO_MAX_LEVS_C = 1

! COMOCASZ start
      COMMON /COMOCASZ/ JO_MAX_OBS_VAL, JO_NMAX_OBS_ICE
! COMOCASZ end
! TYPOASZ end
      ! OCEAN end

      !  WAVE sub-model start
      ! WAVE end

      ! Grid related sizes for COUPLING between atmos and OCEAN
      ! submodels [For mpp, sizes are 'global' values over all
      ! PEs.Also needed for river routing]
      INTEGER:: AOCPL_IMT                ! Ocean rowlength
      INTEGER:: AOCPL_JMT                ! Ocean no. of rows
      INTEGER:: AOCPL_ROW_LENGTH         ! Atmos rowlength
      INTEGER:: AOCPL_P_ROWS             ! Atmos no. of p rows

      COMMON/SIZE_AOCPL/                                                &
     &  AOCPL_IMT, AOCPL_JMT, AOCPL_ROW_LENGTH, AOCPL_P_ROWS

      ! Data structure sizes for ATMOSPHERE ANCILLARY file control
      ! routines
      INTEGER :: NANCIL_LOOKUPSA  ! IN: Max no of fields to be read

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER::N_INTF_A          ! No of atmosphere interface areas
      INTEGER::MAX_INTF_MODEL_LEVELS ! Max no of model levs in all areas
      INTEGER::MAX_LBCROW_LENGTH ! Max no of lbc row length in all areas
      INTEGER::MAX_LBCROWS ! Max no of lbc rows in all areas
      INTEGER::TOT_LEN_INTFA_P   ! Total length of interface p grids.
      INTEGER::TOT_LEN_INTFA_U    ! Total length of interface u grids.
      INTEGER::U_FIELD_INTFA      ! Length of Model U field (= U_FIELD)

      !  Data structure sizes for ATMOSPHERE BOUNDARY file control
      ! routines
! PARAMETERs defining the RIM_TYPE characteristics

!  History:
!  Date      Vn     Modification
!  31/10/01  5.3    Reset Nrima_max to 1. Change rima_type_orog
!                   from 2 to 1. D. Robinson

      INTEGER, PARAMETER :: Nrima_max = 1

      INTEGER, PARAMETER :: rima_type_norm=1  ! Normal field
      INTEGER, PARAMETER :: rima_type_orog=1  ! Orography field

! At 5.3 rima_type_orog=2 => rima_type_orog=1. This means that
! Orography LBCs has the same rim type as all prognostic LBCs. The
! value was changed rather than remove all occurences from the code to
! retain the functionality if it ever needs to be restored and to
! simplify the number of changes required.
      INTEGER :: RIMWIDTHA(Nrima_max)
      INTEGER :: NRIM_TIMESA      ! IN: Max no of timelevels in rim flds

      ! Data structure sizes for atmos & OCEAN BOUNDARY file control
      !routines

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)

      INTEGER :: PP_LEN_INTHD   ! IN: Length of PP file integer header
      INTEGER :: PP_LEN_REALHD  ! IN: Length of PP file real    header

      ! Other sizes passed from namelist into common blocks
      COMMON/NLSIZES/                                                   &
     &  ROW_LENGTH,global_ROW_LENGTH,ROWS,global_ROWS,                  &
     &  LAND_FIELD,MODEL_LEVELS,WET_LEVELS,                             &
     &  NTILES,                                                         &
     &  CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS,           &
     &  OZONE_LEVELS,TPPS_OZONE_LEVELS,TR_VARS,TR_UKCA,                 &
     &  RIVER_ROWS, RIVER_ROW_LENGTH,                                   &

     &  THETA_PV_P_LEVS, N_AOBS,                                        &

     &  A_PROG_LOOKUP,A_PROG_LEN,                                       &
     &  A_LEN_INTHD,A_LEN_REALHD,                                       &
     &  A_LEN2_LEVDEPC,A_LEN2_ROWDEPC,A_LEN2_COLDEPC,                   &
     &  A_LEN2_FLDDEPC,A_LEN_EXTCNST,                                   &
     &  A_LEN_CFI1,A_LEN_CFI2,A_LEN_CFI3,                               &

     &  NANCIL_LOOKUPSA,                                                &

     &  N_INTF_A, MAX_INTF_MODEL_LEVELS, MAX_LBCROW_LENGTH,             &
     &  MAX_LBCROWS, TOT_LEN_INTFA_P,                                   &
     &  TOT_LEN_INTFA_U, U_FIELD_INTFA,                                 &

     &  RIMWIDTHA, NRIM_TIMESA,                                         &

     &  PP_LEN_INTHD,PP_LEN_REALHD

      !-----------------------------------------------------------------
      ! data in STASHC#x member of the job library

      ! Data structure sizes for ATMOSPHERE submodel (config dependent)
      INTEGER:: A_LEN2_LOOKUP   ! IN: Total no of fields (incl diags)
      INTEGER:: A_LEN_DATA      ! IN: Total no of words of data
      INTEGER:: A_LEN_D1        ! IN: Total no of words in atmos D1
      ! Data structure sizes for SLAB submodel (config dependent)
      INTEGER:: S_LEN2_LOOKUP   !IN: Tot no of fields (incl diags)
      INTEGER:: S_LEN_DATA      !IN: Tot no of words of data
      ! Data structure sizes for OCEAN submodel (config dependent)
      INTEGER:: O_LEN2_LOOKUP     ! IN: Total no of fields (incl diags)
      INTEGER:: O_LEN_DATA        ! IN: Total no of words of data
      INTEGER:: O_LEN_DUALDATA    ! IN: Words of data at 2 time levels
      INTEGER:: O_LEN_D1          ! IN: Total no of words in ocean D1
      ! Data structure sizes for WAVE submodel (config dependent)
      INTEGER:: W_LEN2_LOOKUP     ! IN: Total no of fields (incl diags)
      INTEGER:: W_LEN_DATA        ! IN: Total no of words of data
      INTEGER:: W_LEN_D1          ! IN: Total no of words in atmos D1

      ! Size of main data array for this configuration

      INTEGER:: LEN_TOT             ! IN: Length of D1 array
      INTEGER:: N_OBJ_D1_MAX         ! IN: No of objects in D1 array

      COMMON/STSIZES/                                                   &
     &  S_LEN2_LOOKUP,S_LEN_DATA,                                       &
     &  A_LEN2_LOOKUP,A_LEN_DATA,A_LEN_D1,                              &
     &  O_LEN2_LOOKUP,O_LEN_DATA,O_LEN_DUALDATA,O_LEN_D1,               &
     &  W_LEN2_LOOKUP,W_LEN_DATA,W_LEN_D1,                              &
     &  LEN_TOT,N_OBJ_D1_MAX
      ! global (ie. dump version) of *_LEN_DATA
      INTEGER:: global_A_LEN_DATA
      INTEGER:: global_O_LEN_DATA
      INTEGER :: global_W_LEN_DATA ! global (ie. dump version) of
                                   !                        W_LEN_DATA

      COMMON /MPP_STSIZES_extra/ global_A_LEN_DATA,global_O_LEN_DATA    &
     &      ,global_W_LEN_DATA
      ! Sizes of Stash Auxillary Arrays and associated index arrays
      ! Initialised in UMINDEX and UMINDEX_A/O/W
      INTEGER:: LEN_A_IXSTS
      INTEGER:: LEN_A_SPSTS
      INTEGER:: LEN_O_IXSTS
      INTEGER:: LEN_O_SPSTS
      INTEGER:: LEN_W_IXSTS
      INTEGER:: LEN_W_SPSTS

      COMMON /DSIZE_STS/                                                &
     &  LEN_A_IXSTS, LEN_A_SPSTS,LEN_O_IXSTS, LEN_O_SPSTS,              &
     &  LEN_W_IXSTS, LEN_W_SPSTS
!     From 4.5, the number of land points is computed for each
!     PE before the addressing section. All prognostics on land
!     points in the D1 space are now dimensioned by the local
!     no of land points rather than the global no of land points.

      INTEGER:: global_land_field    !  Global no of land points
      INTEGER:: local_land_field     !  Local no of land points
      COMMON /mpp_landpts/ global_land_field,local_land_field
      ! ----------------------------------------------------------------
      ! extra variables not passed through user interface

      ! fundamental data sizes :
      ! Fundamental parameter  sizes of data structure
      ! Sizes applicable to all configurations (HISTORY FILE)

      ! Length of history file in dump
      INTEGER, PARAMETER :: LEN_DUMPHIST = 0

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)
      ! Length of dump fixed header
      INTEGER, PARAMETER :: LEN_FIXHD = 256

      ! Size of a single LOOKUP header
      INTEGER, PARAMETER :: LEN1_LOOKUP  = 64
      INTEGER, PARAMETER :: MPP_LEN1_LOOKUP= 2

      ! Size of compressed LBC LOOKUP (only used internally and
      ! contains just the items which change between each set of LBCs
      INTEGER, PARAMETER :: LEN1_LBC_COMP_LOOKUP = 8

      ! Sizes applicable to all configurations (STASH)
      ! Moved to <typstsz/typstsz.h>

      INTEGER:: INTF_LEN2_LEVDEPC !1st dim of interface out lev dep cons
      INTEGER:: INTF_LEN2_ROWDEPC !2nd dim of interface out Row dep cons
      INTEGER:: INTF_LEN2_COLDEPC !2nd dim of interface out Col dep cons
      
      COMMON /DSIZE/                                                    &
     &  INTF_LEN2_LEVDEPC,INTF_LEN2_ROWDEPC,INTF_LEN2_COLDEPC
      ! sub-model atmosphere   :
      ! Data structure sizes derived from grid size
      INTEGER:: A_LEN1_LEVDEPC ! IN: 1st dim of level  dep const
      INTEGER:: A_LEN1_ROWDEPC ! IN: 1st dim of row    dep const
      INTEGER:: A_LEN1_COLDEPC ! IN: 1st dim of column dep const
      INTEGER:: A_LEN1_FLDDEPC ! IN: 1st dim of field  dep const

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER :: INTF_LOOKUPSA        ! No of interface lookups.
      INTEGER :: OLD_INTF_LOOKUPSA    ! No of interface lookups
                                      ! for old (4.5) LBCs.
      COMMON /DSIZE_A/                                                  &
     &  A_LEN1_LEVDEPC,A_LEN1_FLDDEPC,A_LEN1_ROWDEPC,A_LEN1_COLDEPC,    &
     &  INTF_LOOKUPSA,OLD_INTF_LOOKUPSA

      ! sub-model atmosphere   : derived sizes
      ! derived from model grid/levels. Arakawa B-grid

                                  ! Size of fields on THETA grid:
      INTEGER :: THETA_FIELD_SIZE     ! IN: with no halos
      INTEGER :: THETA_OFF_SIZE       ! IN: with simple halos
      INTEGER :: THETA_HALO_SIZE      ! IN: with extended halos

                                  ! Size of fields on U grid:
      INTEGER :: U_FIELD_SIZE         ! IN: with no halos
      INTEGER :: U_OFF_SIZE           ! IN: with simple halos
      INTEGER :: U_HALO_SIZE          ! IN: with extended halos

                                  ! Size of fields on V grid
      INTEGER :: V_FIELD_SIZE         ! IN: with no halos
      INTEGER :: V_OFF_SIZE           ! IN: with simple halos
      INTEGER :: V_HALO_SIZE          ! IN: with extended halos

      INTEGER :: N_ROWS               ! IN: No of V-rows
      INTEGER :: N_CCA_LEV            ! IN: No of CCA Levels
      COMMON/DRSIZE_A/                                                  &
     &  N_ROWS,N_CCA_LEV,THETA_FIELD_SIZE,U_FIELD_SIZE,V_FIELD_SIZE,    &
     &  THETA_OFF_SIZE,U_OFF_SIZE,V_OFF_SIZE,                           &
     &  THETA_HALO_SIZE,U_HALO_SIZE,V_HALO_SIZE
      ! boundary updating      : derived values
      ! Variables describing the Atmosphere Lateral Boundary Conditions
      ! Local (per processor) information

      ! Size of atmos LBC for given field type, halo type and rimwidth
      ! type
      INTEGER:: LENRIMA(Nfld_max,NHalo_max,Nrima_max)

      ! Size of given side (PNorth,PEast,PSouth and PWest), field type,
                        ! halo type and rimwidth type
      INTEGER:: LBC_SIZEA(4,Nfld_max,NHalo_max,Nrima_max)

                        ! Start of a given side within the LBC
      INTEGER:: LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max)

      ! Start of a given side within the LBC on a given processor
      INTEGER:: g_LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max,0:Maxproc-1)

      ! and global (within the file) information

      ! Size of atmos LBC on disk for given field type, halo type and
      ! rimwidth type
      INTEGER:: global_LENRIMA(Nfld_max,NHalo_max,Nrima_max)

                        ! Size of given side, field type and halo type
      INTEGER:: global_LBC_SIZEA(4,Nfld_max,NHalo_max,Nrima_max)

                        ! Start of a given side within the LBC
      INTEGER:: global_LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max)

      ! Variables describing the Ocean Lateral Boundary Conditions
      INTEGER:: LENRIMO                ! Size of ocean LBC (theta)
      INTEGER:: LENRIMO_U              ! Size of ocean LBC (velocity)

      ! Variables that may be needed for vn5.2 but have not yet been
      ! dealt with at vn5.1
      INTEGER:: RIMFLDSA
      INTEGER:: RIMFLDSO
      INTEGER:: global_LENRIMDATA_A
      INTEGER:: LENRIMDATA_A
      INTEGER:: LENRIMDATA_O
      INTEGER:: BOUNDFLDS
      INTEGER:: RIM_LOOKUPSA
      INTEGER:: RIM_LOOKUPSO
      INTEGER:: BOUND_LOOKUPSA
      INTEGER:: BOUND_LOOKUPSO
      INTEGER:: BOUND_LOOKUPSW
      INTEGER :: RIM_LOOKUPSW
      INTEGER :: LENRIMDATA_W
      INTEGER :: global_LENRIMDATA_W
      COMMON/DRSIZ_BO/                                                  &
      ! Atmosphere variables
     &  LENRIMA, LBC_SIZEA, LBC_STARTA, g_LBC_STARTA,                   &
     &  global_LENRIMA,global_LBC_SIZEA,global_LBC_STARTA,              &
      ! Wave model variables
     & RIM_LOOKUPSW, LENRIMDATA_W, global_LENRIMDATA_W,                 &
      ! Ocean variables
     &  LENRIMO, LENRIMO_U,                                             &
      ! Variables still to be dealt with
     &  RIMFLDSA,RIMFLDSO,BOUNDFLDS,RIM_LOOKUPSA,RIM_LOOKUPSO,          &
     &  BOUND_LOOKUPSA,BOUND_LOOKUPSO,BOUND_LOOKUPSW,                   &
     &  global_LENRIMDATA_A,                                            &
     &  LENRIMDATA_A,LENRIMDATA_O
! TYPSIZE end
! COMDECK CSTASH
! Description:
!  Contains variables and arrays involved in STASH request specification
!  and STASH processing in the UM. Includes namelist STASH.
!   NOTE: comdeck VERSION should be *CALLed before this comdeck.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Mar. 95   Original code.  S.J.Swarbrick
! 4.0       Sept.95   Original code.  S.J.Swarbrick
! 4.1  06/02/96  Comdeck renamed from STASH to CSTASH to avoid clashes
!                 with deck STASH1 in html searches.  RTHBarnes.
! 4.1       May 96    Add array MODL_T - for correct processing
!                      of output times tables  S.J.Swarbrick
! 4.4       Sep 97    Add IOFF_T to allow offset for sampling
!                     S.D.Mullerworth
! 5.0       23/06/99  Added halo_type information from ppxref file
!                                                         P.Burton
! 5.5       28/01/03  Change IOPN(4) to IOPN(6) to cater
!                     for 30 digit option codes.
!                     W Roseblade
!
! Declarations:
! Imported global variables:
!    None, but see note above.

! Global parameters:

! Global scalars:
      INTEGER      NDIAG   ! No. of diagnostics
      INTEGER      NTPROF  ! No. of time profiles
      INTEGER      NSERIES ! No. of stash time series
      INTEGER      NDPROF  ! No. of domain profiles
      INTEGER      NUPROF  ! No. of useage profiles

! Global dynamic arrays:

!   STASH specification table (JSTASH file):
!   NDIAGPM set in comdeck VERSION
      INTEGER      MODL_B(NDIAGPM)  ! Internal model no.
      INTEGER      ISEC_B(NDIAGPM)  ! Section
      INTEGER      ITEM_B(NDIAGPM)  ! Item
      INTEGER      ITIM_B(NDIAGPM)  ! Time profile number
      INTEGER      IDOM_B(NDIAGPM)  ! Domain profile number
      INTEGER      IUSE_B(NDIAGPM)  ! Useage profile number

!   Time profile information:

      CHARACTER*8  TIMPRO(NPROFTP)         ! Name of profile
      INTEGER      ITYP_T(NPROFTP)         ! Type of profile
      INTEGER      INTV_T(NPROFTP)         ! Time Interval
      CHARACTER*2  UNT1_T(NPROFTP)         ! Units for time interval
      INTEGER      ISAM_T(NPROFTP)         ! Sampling period
      CHARACTER*2  UNT2_T(NPROFTP)         ! Units for sampling period
      INTEGER      IOPT_T(NPROFTP)         ! Output option
      INTEGER      ISTR_T(NPROFTP)         ! Output Start time
      INTEGER      IEND_T(NPROFTP)         ! Output End time
      INTEGER      IFRE_T(NPROFTP)         ! Output frequency
      INTEGER      IOFF_T(NPROFTP)         ! Offset for sampling
      CHARACTER*2  UNT3_T(NPROFTP)         ! Units for output times
      INTEGER      ITIM_T(NPROFTP)         ! No. of times in times table
      INTEGER      ISER_T(NTIMEP ,NPROFTP) ! Times table (with units)
      INTEGER      MODL_T(NPROFTP)         ! Indicates internal model
                                           !  for each times table

!   Domain profile information:

      CHARACTER*8 DOMPRO  (NPROFDP)           ! Name of domain profile
      INTEGER     IOPL_D  (NPROFDP)           ! Levels option
      INTEGER     LEVB_D  (NPROFDP)           ! Bottom level
      INTEGER     LEVT_D  (NPROFDP)           ! Top level
      INTEGER     IOPA_D  (NPROFDP)           ! Area option
      INTEGER     INTH_D  (NPROFDP)           ! North boundary
      INTEGER     ISTH_D  (NPROFDP)           ! South boundary
      INTEGER     IEST_D  (NPROFDP)           ! East boundary
      INTEGER     IWST_D  (NPROFDP)           ! West boundary
      INTEGER     IMSK_D  (NPROFDP)           ! Mask type
      INTEGER     IMN_D   (NPROFDP)           ! Meaning option
      INTEGER     IWT_D   (NPROFDP)           ! Weighting option
      CHARACTER*1 TS_D    (NPROFDP)           ! Time series profile
      INTEGER     IG_TS
      INTEGER     I1_TS
      INTEGER     I51_TS
      INTEGER     BLIM_TS (NTimSerP)
      INTEGER     TLIM_TS (NTimSerP)
      REAL        BLIMR_TS(NTimSerP)
      REAL        TLIMR_TS(NTimSerP)
      INTEGER     NLIM_TS (NTimSerP)
      INTEGER     SLIM_TS (NTimSerP)
      INTEGER     ELIM_TS (NTimSerP)
      INTEGER     WLIM_TS (NTimSerP)
      INTEGER     ILEV_D  (NPROFDP)           ! Output levels code
      INTEGER     LEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      REAL       RLEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      INTEGER     PLT_D   (NPROFDP)
      INTEGER     PLLEN_D (NPROFDP)
      INTEGER     PLPOS_D (NPROFDP)
      INTEGER     PSLIST_D(NPSLEVP ,NPSLISTP)
      INTEGER     NPSLISTS
      EQUIVALENCE        (RLEVLST_D,LEVLST_D)

! Useage information:

      CHARACTER*8 USEPRO(NPROFUP)   ! Name of useage profile
      INTEGER     LOCN_U(NPROFUP)   ! Storage location of profile
      INTEGER     IUNT_U(NPROFUP)   ! Unit no.

! Information from ppxref file:

      INTEGER      MODEL_ST       ! Internal model number
      INTEGER      ISPACE         ! Space code
      INTEGER      ITIMA          ! Time availability code
      INTEGER      IGP            ! Grid of data code
      INTEGER      ILEV           ! Level type code
      INTEGER      IBOT           ! First level code
      INTEGER      ITOP           ! Last level code
      INTEGER      IFLAG          ! Level compression flag
      INTEGER      IOPN(6)        ! Sectional option code
      INTEGER      VMSK           ! Integer equiv of bin vers mask
      INTEGER      IPSEUDO        ! Pseudo dimension type
      INTEGER      IPFIRST        ! First pseudo dim code
      INTEGER      IPLAST         ! Last pseudo dim code
      INTEGER      PTR_PROG       ! Section zero point back
      INTEGER      HALO_TYPE      ! Type of halo the field has

! PP output file units
      INTEGER      PPlen2LkUp(OUTFILE_S:OUTFILE_E)
      CHARACTER*1  FTOutUnit (OUTFILE_S:OUTFILE_E)

! User ppxref files
      INTEGER      N_USTASH        ! Number of user ppxref files
      INTEGER      NRECS_USTASH    ! Total no. of user stash records
      CHARACTER*8  USTSFILS(20)    ! Names of user ppxref files
      NAMELIST/USTSNUM /N_USTASH,NRECS_USTASH,USTSFILS

! COMMON blocks:
      COMMON/STCHA/ TIMPRO,UNT1_T,UNT2_T,UNT3_T,DOMPRO,TS_D,            &
     &  USEPRO,USTSFILS,FTOutUnit

      COMMON/STSH/                                                      &
     &  NDIAG   ,MODL_B  ,ISEC_B ,ITEM_B  ,ITIM_B  ,IDOM_B  ,IUSE_B,    &
     &  NTPROF  ,ITYP_T  ,INTV_T ,ISAM_T  ,ITIM_T  ,                    &
     &  IOPT_T  ,ISTR_T  ,IEND_T ,IFRE_T  ,IOFF_T, ISER_T  ,MODL_T  ,   &
     &  NDPROF  ,IOPL_D  ,LEVB_D ,                                      &
     &  IOPA_D  ,INTH_D  ,ISTH_D ,IEST_D  ,IWST_D  ,                    &
     &  IMSK_D  ,IMN_D   ,IWT_D  ,                                      &
     &  LEVT_D  ,LEVLST_D,                                              &
     &  PLT_D   ,PLLEN_D ,PLPOS_D,PSLIST_D,NPSLISTS,                    &
     &  BLIM_TS ,TLIM_TS ,BLIMR_TS,TLIMR_TS,IG_TS   ,I1_TS   ,          &
     &  NLIM_TS ,SLIM_TS ,ELIM_TS ,WLIM_TS ,I51_TS  ,NSERIES ,          &
     &  NUPROF  ,LOCN_U  ,IUNT_U ,                                      &
     &  MODEL_ST,ISPACE  ,ITIMA  ,IGP     ,                             &
     &  ILEV    ,IBOT    ,ITOP   ,IFLAG   ,IOPN    ,VMSK    ,           &
     &  IPSEUDO ,IPFIRST ,IPLAST ,PTR_PROG, HALO_TYPE,                  &
     & N_USTASH,NRECS_USTASH,                                           &
     & PPlen2LkUp

! CSTASH end
! Description:
!   Contains variables and arrays involved in STASH
!   processing in the UM.
!   NOTE: comdecks CSUBMODEL and VERSION must be
!        *CALLed before this comdeck.
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Mar. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 5.2       25/08/00  Add another level of info (section number) to
!                     the D1_PADDR array
!                                                          P.Burton
! 6.2       03/02/06  Increase Max_D1_Len to 1500. T Johns
!
! Declarations:
! Imported global variables:
!    None, but see note above.

! Global arrays:
!   Output levels lists
!     List type (real/int)
      CHARACTER*1 LLISTTY  (NPROFDP*6            )
!     Real levels
      REAL        RLEVLST_S(NLEVP_S  ,  NLEVLSTSP)
!     Integer (i.e. model) levels
      INTEGER      LEVLST_S(NLEVP_S  ,  NLEVLSTSP)
!   STASH lengths and addresses
      INTEGER IN_S    (2,N_INTERNAL_MODEL_MAX,0:NSECTP,NITEMP)
!   STASH list index
      INTEGER INDX_S  (2,N_INTERNAL_MODEL_MAX,0:NSECTP,NITEMP)

!   STASH list array (extra row only for internal indexing in
!                   processing routines)
      INTEGER LIST_S  (NELEMP+1             , NRECDP   )
!   Output times tables
      INTEGER ITIM_S  (NTIMEP              ,2*NPROFTP+2)
!   Start addresses for pp headers
      INTEGER PPIND_S (N_INTERNAL_MODEL_MAX,  NITEMP   )
!   Time series block information
!     No. of records in a block
      INTEGER NRECS_TS(NPROFDP                         )
!     Start position of block
      INTEGER NPOS_TS (NPROFDP                         )
!   lengths of pseudo-levels lists
      INTEGER LENPLST (NPSLISTP                        )

      EQUIVALENCE(LEVLST_S,RLEVLST_S)

!     Set up preliminary array for addressing D1:
!     Number of items of info needed for each object and likely maximum
!     number of objects in D1 - this can be increased if necessary

      Integer, Parameter :: D1_Items_Prel = 5
      Integer, Parameter :: Max_D1_Len    = 1500

      ! Names of items

      Integer, Parameter :: d1_type       = 1 ! Prognostic, diagnostic
                                              ! or other
      Integer, Parameter :: d1_im         = 2 ! Internal model id
      Integer, Parameter :: d1_extra_info = 3 ! Progs and other :-
                                              ! PPXREF item no
                                              ! Diags :-
                                              ! Stash list item no
      Integer, Parameter :: d1_levs       = 4 ! No of levels
      Integer, Parameter :: d1_sect       = 5 ! Section No

      ! Types of items for d1_type

      Integer, Parameter :: Prog     = 0
      Integer, Parameter :: Diag     = 1
      Integer, Parameter :: Seco     = 2
      Integer, Parameter :: Extra_d1 = 3

      ! Stores number of objects in D1
      INTEGER      N_OBJ_D1(N_SUBMODEL_PARTITION_MAX)

!     Preliminary array for addressing D1. Holds the minimum amount of
!     info required for order of objects of D1; this amount of info is
!     enough to obtain any other required info from stashlist or ppxref

      INTEGER :: D1_PADDR(D1_ITEMS_PREL,MAX_D1_LEN,                     &
     &  N_SUBMODEL_PARTITION_MAX)

       COMMON/CHARLIST/  LLISTTY
      COMMON/STEXTEND/ LIST_S,INDX_S,ITIM_S,IN_S,PPIND_S,LEVLST_S,      &
     &  NRECS_TS,NPOS_TS,LENPLST
       COMMON/D1_PRELIM/ D1_PADDR, N_OBJ_D1

! STEXTEND end
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
! STPARAM
!
!  Purpose: Meaningful PARAMETER names for STASH processing routines.
!           Both a long name and short name have been declared, to
!           reduce the existence of "magic" numbers in STASH.
!           Format is that first the address of the item is declare in
!           both long and short form. example is;
!             integer st_item_code,s_item  !Item number (declaration)
!             parameter(st_item_code=3,s_item=3)
!
!  Author:   S.Tett             Date:           22 January 1991
!
!  Model            Modification history from model version 3.0:
! version  Date
!   3.5    Mar. 95  Sub-models project.
!                   st_model_code=28 added to STLIST addresses
!                                   S.J.Swarbrick
!   4.2    27/11/96 mpp code: Added new stlist "magic numbers" :
!                   st_dump_output_length, st_dump_output_addr
!                                                       P.Burton
!   4.4    23/09/97 Add st_offset_code to the STASH list
!                   S.D. Mullerworth
!    4.4  02/12/96 Time mean timeseries added R A Stratton.
!    4.5  23/01/98 Added new stlist magic number
!                  st_dump_level_output_length
!    4.5  23/01/98 A
!    5.5  28/02/03 Original Modifications for WAM. M.Holt
!         06/08/00 Modification for parallelisation of WAM.
!                          Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!    6.0  08/09/03 Add st_riv_grid 23. C.Bunton
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!  Logical components covered: D70
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!--------------------------------------------------------------

      ! Internal model number address
      INTEGER,PARAMETER:: st_model_code = 28
      INTEGER,PARAMETER:: s_modl        = 28

      ! Section Number address
      INTEGER,PARAMETER:: st_sect_no_code = 2
      INTEGER,PARAMETER:: s_sect          = 2
      INTEGER,PARAMETER:: st_sect_code    = 2

      INTEGER,PARAMETER:: st_item_code=1,s_item=1 ! Item number address

      ! Processing Code address
      INTEGER,PARAMETER:: st_proc_no_code=3,s_proc=3

      ! subsidiary codes for st_proc_no_code now
      INTEGER,PARAMETER:: st_replace_code=1
      INTEGER,PARAMETER:: st_accum_code=2
      INTEGER,PARAMETER:: st_time_mean_code=3
      INTEGER,PARAMETER:: st_time_series_code=4
      INTEGER,PARAMETER:: st_max_code=5
      INTEGER,PARAMETER:: st_min_code=6
      INTEGER,PARAMETER:: st_append_traj_code=7
      INTEGER,PARAMETER:: st_time_series_mean=8
      INTEGER,PARAMETER:: st_variance_code=9

      ! Frequency (Input & output) addres
      INTEGER,PARAMETER:: st_freq_code=4,s_freq=4

      ! Offset for sampling
      INTEGER,PARAMETER:: st_offset_code=30,s_offs=30

      ! start timestep address
      INTEGER,PARAMETER:: st_start_time_code=5,s_times=5

      ! end timestep address
      INTEGER,PARAMETER:: st_end_time_code=6,s_timee=6

      ! period in timesteps address
      INTEGER,PARAMETER:: st_period_code=7,s_period=7

      ! infinite end/period value
      INTEGER,PARAMETER:: st_infinite_time=-1

      INTEGER,PARAMETER:: st_end_of_list=-1 !end-of-list marker in times

      ! grid point stuff
      ! gridpoint info address
      INTEGER,PARAMETER:: st_gridpoint_code=8,s_grid=8

      ! now subsid grid point stuff
      ! no masking done
      INTEGER,PARAMETER:: stash_null_mask_code=1,s_nomask=1

      ! land mask conds
      INTEGER,PARAMETER:: stash_land_mask_code=2,s_lndms=2

      ! sea mask code
      INTEGER,PARAMETER:: stash_sea_mask_code=3,s_seams =3

      ! processing options

      ! size of block for gridpoint code
      INTEGER,PARAMETER:: block_size=10

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: extract_base=block_size*0

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: extract_top=block_size*1

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_base=block_size*1

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_top=block_size*2

      ! max code for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_base=block_size*2

      ! base codes for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_top=block_size*3

      ! max code for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_base=block_size*3

      ! base codes for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_top=block_size*4

      ! max code for field mean subroutine
      INTEGER,PARAMETER:: field_mean_base=block_size*4

      ! base codes for field mean subroutine
      INTEGER,PARAMETER:: field_mean_top=block_size*5

      ! max code for global mean subroutine
      INTEGER,PARAMETER:: global_mean_base=block_size*5

      ! base codes for global mean subroutine
      INTEGER,PARAMETER:: global_mean_top=block_size*6

      ! Weighting

      ! weighting info address
      INTEGER,PARAMETER:: st_weight_code=9,s_weight=9

      INTEGER,PARAMETER:: stash_weight_null_code  =0,s_noweight  =0
      INTEGER,PARAMETER:: stash_weight_area_code  =1,s_areaweight=1
      INTEGER,PARAMETER:: stash_weight_volume_code=2,s_volweight =2
      INTEGER,PARAMETER:: stash_weight_mass_code  =3,s_massweight=3

      ! Domain definition

      ! row addresses
      INTEGER,PARAMETER:: st_north_code=12,s_north=12
      INTEGER,PARAMETER:: st_south_code=13,s_south=13
      INTEGER,PARAMETER:: st_west_code =14,s_west =14
      INTEGER,PARAMETER:: st_east_code =15,s_east =15

      ! Levels

      ! input bottom level address
      INTEGER,PARAMETER:: st_input_bottom=10,s_bottom =10

      ! special code
      INTEGER,PARAMETER:: st_special_code=100,s_special=100

      ! input top level address
      INTEGER,PARAMETER:: st_input_top=11,s_top=11

      ! output bottom level address
      INTEGER,PARAMETER:: st_output_bottom=21,s_outbot=21

      ! output top level address
      INTEGER,PARAMETER:: st_output_top=22,s_outtop=22

      INTEGER,PARAMETER:: st_model_level_code=1,s_model=1

      ! code for pressure leve
      INTEGER,PARAMETER:: st_pressure_level_code=2,s_press=2

      ! code for height levels
      INTEGER,PARAMETER:: st_height_level_code=3,s_height=3

      ! input code addres
      INTEGER,PARAMETER:: st_input_code=16,s_input=16

      ! input length of diagnostic address
      INTEGER,PARAMETER:: st_input_length=17,s_length=17

      ! output code address
      INTEGER,PARAMETER:: st_output_code=18,s_output=18

      ! Pointer to D1 addressing information
      ! Pos of item in D1 for relevant submodel
      INTEGER,PARAMETER:: st_position_in_d1=29,st_d1pos=29

      ! Output destination options

      INTEGER,PARAMETER:: st_dump=1
      INTEGER,PARAMETER:: st_secondary=2

      ! output length of diagnostic address
      INTEGER,PARAMETER:: st_output_length=19,s_outlen=19
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

      ! ptr to dump lookup header address
      INTEGER,PARAMETER:: st_lookup_ptr=23

      ! ptr into stash_series where control data address
      INTEGER,PARAMETER:: st_series_ptr=24

      ! subsid stuff for time series
      INTEGER,PARAMETER:: series_grid_type=1
      INTEGER,PARAMETER:: series_grid_code=0
      INTEGER,PARAMETER:: series_long_code=1
      INTEGER,PARAMETER:: series_size=2
      INTEGER,PARAMETER:: series_proc_code=3
      INTEGER,PARAMETER:: series_north=4
      INTEGER,PARAMETER:: series_south=5
      INTEGER,PARAMETER:: series_west=6
      INTEGER,PARAMETER:: series_east=7
      INTEGER,PARAMETER:: series_list_start=8
      INTEGER,PARAMETER:: series_list_end=9
      INTEGER,PARAMETER:: record_size=9

      ! Miscellaneous parameters

      ! system/user tag field in stlist address
      INTEGER,PARAMETER:: st_macrotag=25

      ! Pseudo-level list pointers

      ! pseudo-levels input list address
      INTEGER,PARAMETER:: st_pseudo_in=26

      ! pseudo-levels output list address
      INTEGER,PARAMETER:: st_pseudo_out=27

      ! Internal horizontal gridtype codes common to all diagnostics

      INTEGER,PARAMETER:: st_tp_grid =1 ! T-p grid
      INTEGER,PARAMETER:: st_uv_grid =2 ! u-v grid
      INTEGER,PARAMETER:: st_cu_grid =3 ! C-grid u point
      INTEGER,PARAMETER:: st_cv_grid =4 ! C-grid v point
      INTEGER,PARAMETER:: st_zt_grid =5 ! Zonal T-grid
      INTEGER,PARAMETER:: st_zu_grid =6 ! Zonal u-grid
      INTEGER,PARAMETER:: st_mt_grid =7 ! Meridional T-grid
      INTEGER,PARAMETER:: st_mu_grid =8 ! Meridional u-grid
      INTEGER,PARAMETER:: st_riv_grid= 23    ! river_routing grid
      INTEGER,PARAMETER:: st_scalar  =9 ! Scalar (ie. single value)
      INTEGER,PARAMETER:: st_wam_all= 60    ! Wam Field on Full Grid
      INTEGER,PARAMETER:: st_wam_sea= 62    ! Wam Field on Sea Points

! STPARAM end

! Subroutine arguments


!   Scalar arguments with intent(out):

      INTEGER NRECS
      INTEGER NTIMES
      INTEGER NLEVELS ! Total no. of sets of levs for diags (inpt+outp)
      CHARACTER*80 CMESSAGE

! ErrorStatus:
      INTEGER ErrorStatus

! Local scalars:
      LOGICAL      MODEL_LEV
      LOGICAL      LMASK
      LOGICAL      LMEAN
      LOGICAL      LOFFSET
      INTEGER      TOTIMP
      INTEGER      I
      INTEGER      IBOT1
      INTEGER      IDIAG
      INTEGER      IDOMLEV
      INTEGER      IDOM_L
      LOGICAL      LDUM
      INTEGER      IFIRST
      INTEGER      IFIRST1
      INTEGER      ILAST
      INTEGER      ILAST1
      INTEGER      IM
      INTEGER      IMD
      INTEGER      IPLOF
      INTEGER      MODL_L
      INTEGER      ISEC_L
      INTEGER      ITEM_L
      INTEGER      ITIM_L
      INTEGER      ITIM
      INTEGER      ITOP1
      INTEGER      IUSE_L
      INTEGER      IX1
      INTEGER      IX2
      INTEGER      IY1
      INTEGER      IY2
      INTEGER      JLEV
      INTEGER      LEV_OFFSET
      INTEGER      LBVC
      INTEGER      IMAX          ! to find max of times-table
      INTEGER      ITIMLST       ! column of times-table
      INTEGER      item_chk

      CHARACTER (LEN=256)          :: CMESSAGE2
      CHARACTER (LEN=*), PARAMETER :: RoutineName='PRELIM'

! Function and subroutine calls:
      LOGICAL  DISCT_LEV
      INTEGER  EXPPXI
      EXTERNAL EXPPXI,LEVSRT,TSTMSK,LLTORC,LEVCOD,PSLCOM,PSLIMS

!- End of Header ------------------------------------------------------

! 0.1  Store output-times tables in array ITIM_S

      IF(NTIMES == 0) THEN
      DO I=1,NPROFTP
        IF (IOPT_T(I) == 2.AND.MODL_T(I) >  0) THEN
! Profile has output times list
!  MODL_T(I) labels internal model for times list
          DO ITIM=1,ITIM_T(I)
! DEPENDS ON: totimp
             ITIM_S(ITIM,I)=TOTIMP(ISER_T(ITIM,I),UNT3_T(I),MODL_T(I))
             if (ITIM_S(ITIM,I)  ==  -999) then
                 ErrorStatus = 100
                 write (cmessage,'(a,a,i3)')                            &
     &           'PRELIM:TOTIMP:Error in time period conversion',       &
     &           ' output times table no.=',i
                 write(6,*) cmessage
                 GOTO 9999
              endif
          END DO
          ITIM_S(ITIM_T(I)+1,I)=-1
        ELSE
          ITIM_S(1,I)=-1
        END IF
      END DO
      NTIMES=NPROFTP
      END IF

! 0.2  Store output levels lists in array LEVLST_S

      LEV_OFFSET=NLEVELS ! Initialised to 0 before entering this routine

! Loop over domain profiles in STASH basis file
      DO I=1,NDPROF
        IF (LEVB_D(I) == -1) THEN
! There is a levels list for this dom prof
          IF (IOPL_D(I) == 1.OR.IOPL_D(I) == 2.OR.                      &
     &                          IOPL_D(I) == 6    ) THEN
! Levs list contains model levs - list type is integer
             LLISTTY(I+LEV_OFFSET)='I'
          ELSE
! Not model levs - list type real
             LLISTTY(I+LEV_OFFSET)='R'
          END IF
! LEVT_D(I) = no. of levs in list 'I'
          LEVLST_S(1,I+LEV_OFFSET)=LEVT_D(I)

! Levels list 'I' was read into (R)LEVLST_D(J,I), J=1,LEVT_D(I),
!  by RDBASIS.
!  Transfer this levels list to (R)LEVLST_S(J,I+LEV_OFFSET),
!  J=2,LEVT_D(I)+1.

          DO JLEV=1,LEVT_D(I)
            IF (IOPL_D(I) == 1.OR.IOPL_D(I) == 2.OR.                    &
     &                            IOPL_D(I) == 6    ) THEN
!         Model levels
               LEVLST_S(JLEV+1,I+LEV_OFFSET)= LEVLST_D(JLEV,I)
            ELSE IF (IOPL_D(I) /= 5) THEN
!         Real levels
              RLEVLST_S(JLEV+1,I+LEV_OFFSET)=RLEVLST_D(JLEV,I)
            END IF
          END DO

          IPLOF=I+LEV_OFFSET

!   Sort this levels list into correct order (if not already in order)
! DEPENDS ON: levsrt
          CALL LEVSRT( LLISTTY(  IPLOF), LEVLST_S(1,IPLOF),             &
     &                LEVLST_S(2,IPLOF),RLEVLST_S(2,IPLOF))
        ELSE
! No levels list, i.e., the output from this diag. is on a
!    contiguous range of model levels
          LEVLST_S(1,I+LEV_OFFSET)=0
        END IF
      END DO  !  Domain profiles

      NLEVELS=NDPROF+LEV_OFFSET  ! NDPROF = no. of sets of input levels

      IF(NLEVELS >  NLEVLSTSP) THEN
        WRITE(6,*)                                                      &
     &  'PRELIM: TOO MANY LEVELS LISTS, ARRAYS OVERWRITTEN'
        CMESSAGE=                                                       &
     & 'PRELIM: TOO MANY LEVELS LISTS, ARRAYS OVERWRITTEN'
        GO TO 9999
      END IF

! Section 1. MAIN LOOP - loop over diag requests in STASH basis file

      IF(NDIAG >  0) THEN

      DO IDIAG=1,NDIAG

      MODL_L=MODL_B(IDIAG)
      ISEC_L=ISEC_B(IDIAG)
      ITEM_L=ITEM_B(IDIAG)
      IDOM_L=IDOM_B(IDIAG)
      IUSE_L=IUSE_B(IDIAG)
      ITIM_L=ITIM_B(IDIAG)

      item_chk=0
! DEPENDS ON: exppxi
      item_chk=EXPPXI(MODL_L,ISEC_L,ITEM_L,ppx_item_number,             &
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
     &                ErrorStatus,CMESSAGE)
      IF(ITEM_CHK /= ITEM_L)THEN
        WRITE(CMESSAGE2,*)'Diagnostic discarded ',                      &
     &       ' model ',modl_l,' section ',isec_l,' item',item_l,        &
     &       ' No stashmaster record'
        ERRORSTATUS = -10
! DEPENDS ON: ereport
        CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)

        ! Make this item null
        ITEM_B(IDIAG)=0
        GOTO 999
      ENDIF
      IF(ITIM_L /= 0) THEN       ! If the diag is not a null request

! Section 1.0  Extract data required for STASH processing from PPXI

        IF(NRECS == NRECDP) THEN
          WRITE(CMESSAGE2,*)                                            &
     &   'TOO MANY STASH LIST ENTRIES, REQUEST DENIED',                 &
     &   ' (M,S,I)',MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS = -20
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          GOTO 999
        END IF

! DEPENDS ON: exppxi
        VMSK    = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_version_mask ,      &
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
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ISPACE  = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_space_code   ,      &
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
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ITIMA   = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_timavail_code,      &
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
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IGP     = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_grid_type    ,      &
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
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ILEV    = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_lv_code      ,      &
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
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IBOT    = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_lb_code      ,      &
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
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        ITOP    = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_lt_code      ,      &
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
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IFLAG   = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_lev_flag     ,      &
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
     &                                          ErrorStatus,CMESSAGE)
        DO I=1,6


! DEPENDS ON: exppxi
        IOPN(I) = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_opt_code+I-1 ,      &
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
     &                                          ErrorStatus,CMESSAGE)
        END DO
! DEPENDS ON: exppxi
        IPSEUDO = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_pt_code      ,      &
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
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPFIRST = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_pf_code      ,      &
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
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        IPLAST  = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_pl_code      ,      &
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
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        PTR_PROG= EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_ptr_code     ,      &
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
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
        LBVC    = EXPPXI(MODL_L ,ISEC_L ,ITEM_L,ppx_lbvc_code    ,      &
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
     &                                          ErrorStatus,CMESSAGE)

! Check availability of diagnostic
! DEPENDS ON: tstmsk
        CALL TSTMSK(MODL_L,ISEC_L,LMASK,LDUM,ErrorStatus,CMESSAGE)
        IF(.NOT.LMASK) THEN
          WRITE(CMESSAGE2,*)                                            &
     &   'DIAGNOSTIC NOT AVAILABLE TO THIS VERSION ',                   &
     &   'REQUEST DENIED ',                                             &
     &   '(M,S,I)',                                                     &
     &                MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS = -30
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          GOTO 999
        END IF

        NRECS=NRECS+1

        LIST_S(st_model_code  ,NRECS)= MODL_L
        LIST_S(st_sect_no_code,NRECS)= ISEC_L
        LIST_S(st_item_code   ,NRECS)= ITEM_L
! Prelim pointer for 'child' records
        LIST_S(NELEMP+1       ,NRECS)= NRECS
        LIST_S(st_lookup_ptr  ,NRECS)=-1

! Set input code for STASH requests:
!  =0 Use primary or secondary field:       D1(SI(item,section,model))
!  =1 Use field in diagnostic space: STASHwork(SI(item,section,model))
!  =-j Use diagnostic at D1(LIST_S(st_output_addr,j))
        IF( (ISPACE == 2).OR.(ISPACE == 4)                              &
     &  .OR.(ISPACE == 7).OR.(ISPACE == 8) .OR.(ISPACE == 9)) THEN
          LIST_S(st_input_code,NRECS)=0
        ELSE
          LIST_S(st_input_code,NRECS)=1
        END IF

        IF((ITIMA >= 5).AND.(ITIMA <= 12)) THEN
          LMEAN=.TRUE.
        ELSE
          LMEAN=.FALSE.
        END IF


! 1.1   Expand the domain profile ---------------------------

!   Averaging and Weighting
        IM=IMSK_D(IDOM_L)
        IF ((IGP ==  2).OR.(IGP ==  3)    .OR.                          &
     &      (IGP == 12).OR.(IGP == 13))   THEN
! Diags only available over land/sea
          IF((IMSK_D(IDOM_L)  ==  1)      .AND.                         &
     &       (IGP == 3.OR.IGP == 13))     THEN
! Diag requested over land+sea, only available over sea
            IM=3
          ELSE IF((IMSK_D(IDOM_L)  ==  1) .AND.                         &
     &            (IGP == 2.OR.IGP == 12))THEN
! Diag requested over land+sea, only available over land
            IM=2
          ELSE IF((IMSK_D(IDOM_L)  ==  2) .AND.                         &
     &            (IGP == 3.OR.IGP == 13))THEN
! Diag requested over land, only available over sea
            WRITE(6,*)'PRELIM: CHANGED TO SEA DIAG'
            WRITE(6,*) 'MODEL,SECTION,ITEM ',                           &
     &                  MODL_L,ISEC_L,ITEM_L
            IM=3
          ELSE IF((IMSK_D(IDOM_L)  ==  3) .AND.                         &
     &            (IGP == 2.OR.IGP == 12))THEN
! Diag requested over sea, only available over land
            WRITE(6,*)'PRELIM: CHANGED TO LAND DIAG'
            WRITE(6,*) 'MODEL,SECTION,ITEM ',                           &
     &                  MODL_L,ISEC_L,ITEM_L
            IM=2
          END IF
        END IF

        LIST_S(st_gridpoint_code,NRECS)=IM+10*IMN_D(IDOM_L)
        LIST_S(st_weight_code   ,NRECS)=      IWT_D(IDOM_L)

!   Horizontal area
!    - convert lat/long spec to row/column numbers if appropriate;
!    - convert lat/long spec to equatorial lat/long if appropriate.
        IF(IOPA_D(IDOM_L) == 1) THEN
! Full domain
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,90,-90,0,360,                                 &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 2 ) THEN
! N Hemis
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,90,0,0,360,                                   &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 3 ) THEN
! S Hemis
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,0,-90,0,360,                                  &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 4 ) THEN
! 90N-30N
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,90,30,0,360,                                  &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 5 ) THEN
! 30S-90S
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,-30,-90,0,360,                                &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 6 ) THEN
! 30N-00N
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,30,00,0,360,                                  &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 7 ) THEN
! 00S-30S
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,00,-30,0,360,                                 &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 8 ) THEN
! 30N-30S
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,30,-30,0,360,                                 &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 9 ) THEN
! Other lat/long spec
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,INTH_D(IDOM_L),ISTH_D(IDOM_L),                &
     &                    IWST_D(IDOM_L),IEST_D(IDOM_L),                &
     &         LIST_S(st_south_code,NRECS),LIST_S(st_north_code,NRECS), &
     &         LIST_S(st_west_code,NRECS),LIST_S(st_east_code,NRECS))
        ELSE IF(IOPA_D(IDOM_L) == 10) THEN
! Grid point spec
! DEPENDS ON: lltorc
          CALL LLTORC(IGP,90,-90,0,360,IY1,IY2,IX1,IX2)
          LIST_S(st_north_code,NRECS)=MIN(INTH_D(IDOM_L),IY2)
          LIST_S(st_south_code,NRECS)=MIN(ISTH_D(IDOM_L),IY2)
          LIST_S(st_west_code ,NRECS)=MIN(IWST_D(IDOM_L),IX2)
          LIST_S(st_east_code ,NRECS)=MIN(IEST_D(IDOM_L),IX2)
        ELSE
          WRITE(CMESSAGE2,*) 'INVALID DOMAIN AREA OPTION=',             &
     &                      IOPA_D(IDOM_L),                             &
     &                      '(M,S,I)',                                  &
     &                      MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS = -35
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          NRECS=NRECS-1
          GOTO 999
        END IF

! Input level setting
! DEPENDS ON: disct_lev
        MODEL_LEV=DISCT_LEV(ILEV,ErrorStatus,CMESSAGE)
        IF (MODEL_LEV) THEN
! Model levels
! Set bottom level
! DEPENDS ON: levcod
          CALL LEVCOD(IBOT,IBOT1,ErrorStatus,CMESSAGE)
! Set top level
! DEPENDS ON: levcod
          CALL LEVCOD(ITOP,ITOP1,ErrorStatus,CMESSAGE)
! Contig. range of model levels
          IF(IFLAG == 0) THEN
            LIST_S(st_input_bottom,NRECS)=IBOT1
            LIST_S(st_input_top   ,NRECS)=ITOP1
! Non-contig. levels list
          ELSE IF(IFLAG == 1) THEN
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 1
          END IF
        ELSE
! Non-model levels
          IF(ILEV == 3) THEN
!  Pressure levels
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 2
          ELSE IF(ILEV == 4) THEN
!  Height levels
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 3
          ELSE IF(ILEV == 5) THEN
!  Special levels
            LIST_S(st_input_bottom,NRECS)=100
            LIST_S(st_input_top   ,NRECS)=LBVC
          ELSE IF(ILEV == 7) THEN
!  Theta levels
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 4
          ELSE IF(ILEV == 8) THEN
!  PV levels
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 5
          ELSE IF(ILEV == 9) THEN
!  Cloud threshold levels
            LIST_S(st_input_bottom,NRECS)=-1
            LIST_S(st_input_top   ,NRECS)= 6
          END IF
        END IF

! Output level specification
! DEPENDS ON: disct_lev
        MODEL_LEV=DISCT_LEV(ILEV,ErrorStatus,CMESSAGE)
        IF (MODEL_LEV) THEN
! Model levels
          IF (LEVB_D(IDOM_L) >= 0) THEN
! Contiguous range of model levels
            LIST_S(st_output_bottom,NRECS)=MAX(LEVB_D(IDOM_L),IBOT1)
            LIST_S(st_output_top   ,NRECS)=MIN(LEVT_D(IDOM_L),ITOP1)
            IF ((LEVB_D(IDOM_L) <  IBOT1).OR.                           &
     &          (LEVT_D(IDOM_L) >  ITOP1)) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'DIAGNOSTIC HAS LEVEL RANGE OUT OF BOUNDS; CORRECTED ',      &
     &     '(M,S,I) ',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-40
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            END IF
            IF ( (  TS_D(IDOM_L) ==   'Y').AND.                         &
     &          ((LEVB_D(IDOM_L) <  IBOT1).OR.                          &
     &           (LEVT_D(IDOM_L) >  ITOP1))    ) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'TIME SERIES DOMAIN',                                      &
     &       'HAS INCONSISTENT LEVELS; DIAGNOSTIC IGNORED',             &
     &       ' (M,S,I) ',                                               &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-50
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
            IF ((LEVT_D(IDOM_L) <  IBOT1).OR.                           &
     &          (LEVB_D(IDOM_L) >  ITOP1)) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'DIAGNOSTIC HAS TOP/BOT LEVELS INCONSISTENT; DIAG IGNORED',&
     &       ' (M,S,I) ',                                               &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-60
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          ELSE
! Non-contig. list of model levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=1
          END IF
        ELSE
! Non-model levels
          IF(ILEV == 5) THEN
! Special level
            LIST_S(st_output_bottom,NRECS)=100
            LIST_S(st_output_top   ,NRECS)=LBVC
          ELSE IF(ILEV == 3) THEN
! Pressure levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=2
          ELSE IF(ILEV == 4) THEN
! Height levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=3
          ELSE IF(ILEV == 7 ) THEN
! Theta levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=4
          ELSE IF(ILEV == 8 ) THEN
! PV levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=5
          ELSE IF(ILEV == 9 ) THEN
! Cloud threshold levels
            LIST_S(st_output_bottom,NRECS)=-(IDOM_L+LEV_OFFSET)
            LIST_S(st_output_top   ,NRECS)=6
          ELSE
            WRITE(CMESSAGE2,*) 'DOMAIN LEVEL OPTION=',IOPL_D(IDOM_L),   &
     &                        'DIAG IGNORED. (M,S,I) ',                 &
     &                        MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-70
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999
          END IF
        END IF

! Output pseudo-levels level setting
        IF(IPSEUDO /= PLT_D(IDOM_L)) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'DIAGNOSTIC HAS ',                                           &
     &     'INVALID PSEUDO LEVEL TYPE; IGNORED.',                       &
     &     ' (M,S,I)',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-80
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999
        END IF
        LIST_S(st_pseudo_in,NRECS)=0  !(This is set in INPUTL)
        IF(IPSEUDO >  0) THEN
! Pseudo levels list for this diagnostic
            LIST_S(st_pseudo_out,NRECS)=PLPOS_D(IDOM_L)
            LENPLST(PLPOS_D(IDOM_L))   =PLLEN_D(IDOM_L)
            IFIRST=PSLIST_D(1,PLPOS_D(IDOM_L))
            ILAST =PSLIST_D(PLLEN_D(IDOM_L),PLPOS_D(IDOM_L))
! Check pseudo level limits
! DEPENDS ON: pslims
            CALL PSLIMS(IPFIRST,IPLAST,IFIRST1,ILAST1)
            IF(IFIRST <  IFIRST1) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'DIAGNOSTIC HAS ',                                         &
     &       'FIRST PSEUDO LEVEL TOO LOW; IGNORED.',                    &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-90
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
            IF(ILAST >  ILAST1) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'DIAGNOSTIC HAS ',                                         &
     &       'LAST PSEUDO LEVEL TOO HIGH; IGNORED',                     &
     &       ' (M,S,I)',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-95
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
        ELSE
            LIST_S(st_pseudo_out,NRECS)=0
        END IF

! Time-series domain profiles
        IF(TS_D(IDOM_L) == 'Y') THEN
! Pointer for location of time series
            LIST_S(st_series_ptr,NRECS)=NPOS_TS(IDOM_L)
        ELSE
          LIST_S(st_series_ptr,NRECS)=0
        END IF

! 1.2   Expand the useage profile --------------------------

        IF (LOCN_U(IUSE_L) == 5) THEN                    ! PP file

          IF(LMEAN) THEN
            LIST_S(st_output_code,NRECS)=-27
            LIST_S(st_macrotag,NRECS)=0
          ELSE
            WRITE(6,*)                                                  &
     &     'MESSAGE FROM ROUTINE PRELIM: DIAGNOSTIC REQUEST HAS ',      &
     &     'OUTPUT DESTINATION CODE 5 (CLIMATE MEAN PP FILE) ',         &
     &     'BUT DIAGNOSTIC IS NOT A CLIMATE MEAN; REQUEST IGNORED'

            WRITE(CMESSAGE2,*)                                          &
     &     'DIAGNOSTIC IS NOT CLIMATE MEAN, BUT SENT TO MEAN FILE:',    &
     &     'IGNORED. (M,S,I) ',                                         &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-100
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999
          END IF

        ELSE IF (LMEAN) THEN

            WRITE(6,*)                                                  &
     &     'MESSAGE FROM ROUTINE PRELIM: DIAGNOSTIC REQUEST IS A ',     &
     &     'CLIMATE MEAN - SHOULD HAVE OUTPUT DESTINATION CODE 5 ',     &
     &     '(CLIMATE MEAN PP FILE); REQUEST IGNORED'

            WRITE (CMESSAGE2,*)                                         &
     &      'DIAGNOSTIC IS CLIMATE MEAN, INCORRECT DESTINATION:',       &
     &      'IGNORED. (M,S,I) ',                                        &
     &       MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-110
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999

        ELSE IF (LOCN_U(IUSE_L) == 3) THEN                ! PP file

          LIST_S(st_output_code,NRECS)=-IUNT_U(IUSE_L)
          IF (UNT1_T(ITIM_L)=="DU" .AND. UNT1_T(ITIM_L)=="DU") THEN
             ! Special tag for dump frequency means
             LIST_S(st_macrotag,NRECS)=-999
          ELSE
             LIST_S(st_macrotag,NRECS)=0
          ENDIF

        ELSE IF (LOCN_U(IUSE_L) == 1) THEN ! Dump store: set user tag

          LIST_S(st_output_code,NRECS)=1
          LIST_S(st_macrotag,NRECS)=IUNT_U(IUSE_L)

        ELSE IF (LOCN_U(IUSE_L) == 6) THEN ! Secondary dump store:
                                           !             set user tag
          LIST_S(st_output_code,NRECS)=2
          LIST_S(st_macrotag,NRECS)=IUNT_U(IUSE_L)

        ELSE IF (LOCN_U(IUSE_L) == 2) THEN ! Climate mean: tag set
                                           !   1000*(time mean tag)
          LIST_S(st_output_code,NRECS)=1
          LIST_S(st_macrotag,NRECS)=IUNT_U(IUSE_L)*1000

        ELSE IF (LOCN_U(IUSE_L) == 4)THEN  ! Printed output

          LIST_S(st_output_code,NRECS)=7
          LIST_S(st_macrotag,NRECS)=0

        ELSE

          WRITE(CMESSAGE2,*)                                            &
     &    'INVALID USEAGE OPTION=', LOCN_U(IUSE_L),                     &
     &    ': DIAGNOSTIC IGNORED. (M,S,I) ',                             &
     &     MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS=-120
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          NRECS=NRECS-1
          GOTO 999

        END IF

! 1.3   Expand the time profile ------------------------------

! Initialise as single time field

!   Set time processing record

        IF (LMEAN) THEN
          IF (ITYP_T(ITIM_L) /= 1) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'CLIMATE MEANS MUST NOT BE TIME PROCESSED.',                 &
     &     '(M,S,I) ',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-130
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          END IF
          LIST_S(st_proc_no_code,NRECS)=1
        ELSE
          LIST_S(st_proc_no_code,NRECS)=ITYP_T(ITIM_L)
        END IF

! Initialise offset to 0
            LIST_S(st_offset_code,NRECS)=0
!   Set period record

        IF (ITYP_T(ITIM_L) == 1.OR.LMEAN) THEN        ! No period
          LIST_S(st_period_code,NRECS)=0
        ELSE IF((INTV_T(ITIM_L) == -1).AND.                             &
     &          (ITYP_T(ITIM_L) == 2)) THEN
          LIST_S(st_period_code,NRECS)=-1
        ELSE
          LIST_S(st_period_code,NRECS)=                                 &
! DEPENDS ON: totimp
     &           TOTIMP(INTV_T(ITIM_L),UNT1_T(ITIM_L),MODL_L)
          if (LIST_S(st_period_code,NRECS)  ==  -999) then
              ErrorStatus = 101
              write (cmessage,'(a,a,i2,a,i3,a,i3)')                     &
     &        'PRELIM:TOTIMP:Error in time period conversion',          &
     &        ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
              write(6,*) cmessage
              GOTO 9999
           endif
        END IF

        IF (LMEAN.AND.(IOPT_T(ITIM_L) /= 1)) THEN
          WRITE(CMESSAGE2,*)                                            &
     &   'CLIMATE MEANS MUST USE STANDARD FREQUENCY. DIAG IGNORED.',    &
     &   '(M,S,I) ',                                                    &
     &     MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS=-140
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          NRECS=NRECS-1
          GOTO 999
        END IF

        IF(IOPT_T(ITIM_L) == 1) THEN
!Regular output times
          LIST_S(st_freq_code,NRECS)=                                   &
! DEPENDS ON: totimp
     &           TOTIMP(IFRE_T(ITIM_L),UNT3_T(ITIM_L),MODL_L)
          if (LIST_S(st_freq_code,NRECS)  ==  -999) then
              ErrorStatus = 102
              write (cmessage,'(a,a,i2,a,i3,a,i3)')                     &
     &        'PRELIM:TOTIMP:Error in time period conversion',          &
     &        ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
              write(6,*) cmessage
              GOTO 9999
           endif
          LIST_S(st_start_time_code,NRECS)=                             &
! DEPENDS ON: totimp
     &           TOTIMP(ISTR_T(ITIM_L),UNT3_T(ITIM_L),MODL_L)
          if (LIST_S(st_start_time_code,NRECS)  ==  -999) then
              ErrorStatus = 103
              write (cmessage,'(a,a,i2,a,i3,a,i3)')                     &
     &        'PRELIM:TOTIMP:Error in time period conversion',          &
     &        ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
              write(6,*) cmessage
              GOTO 9999
           endif
          IF(IEND_T(ITIM_L) == -1) THEN
             LIST_S(st_end_time_code,NRECS)=-1
          ELSE
             LIST_S(st_end_time_code,NRECS)=                            &
! DEPENDS ON: totimp
     &           TOTIMP(IEND_T(ITIM_L),UNT3_T(ITIM_L),MODL_L)
          ENDIF
          if (LIST_S(st_end_time_code,NRECS)  ==  -999) then
              ErrorStatus = 104
              write (cmessage,'(a,a,i2,a,i3,a,i3)')                     &
     &        'PRELIM:TOTIMP:Error in time period conversion',          &
     &        ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
              write(6,*) cmessage
              GOTO 9999
           endif

!   Set end time to -1 if output requested to end of run

!   Correct start time for radiation, periodic convection, leaf
!   phenology and vegetation competition
          IF((ITIMA == 2).AND.(A_LW_RADSTEP /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_LW_RADSTEP)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
          ELSE IF((ITIMA == 3).AND.(A_SW_RADSTEP /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_SW_RADSTEP)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
          ELSE IF((ITIMA == 13).AND.(A_CONV_STEP /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),A_CONV_STEP)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
          ELSE IF((ITIMA == 14).AND.(PHENOL_PERIOD /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),PHENOL_PERIOD)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
          ELSE IF((ITIMA == 15).AND.(TRIFFID_PERIOD /= 1)) THEN
            IMD=MOD(LIST_S(st_start_time_code,NRECS),TRIFFID_PERIOD)
            LIST_S(st_start_time_code,NRECS)=                           &
     &      LIST_S(st_start_time_code,NRECS)+1-IMD
            LOFFSET=.TRUE.
          ELSE
            LOFFSET=.FALSE.
          END IF
        ELSE IF(IOPT_T(ITIM_L) == 2) THEN
!List of specified output times
            LIST_S(st_freq_code,NRECS)=-ITIM_L
        ELSE
          WRITE(CMESSAGE2,*)                                            &
     &    'INVALID OUTPUT TIMES CODE. DIAG IGNORED.',                   &
     &    '(M,S,I) ',                                                   &
     &     MODL_L,ISEC_L,ITEM_L

          ERRORSTATUS=-150
! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
          NRECS=NRECS-1
          GOTO 999
        END IF

        IF (LMEAN) LIST_S(st_freq_code,NRECS)=1

        IF ((LIST_S(st_proc_no_code,NRECS) >  1).AND.                   &
     &      (LIST_S(st_proc_no_code,NRECS) <= 6)) THEN
! Other than single time field
          IF(NRECS >= NRECDP) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'TOO MANY S_LIST REQUESTS. REQUEST IGNORED',                 &
     &     '(M,S,I) ',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-160
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999
          END IF

          DO I=1,NELEMP+1          ! Copy stash list forward
            LIST_S(I,NRECS+1)=LIST_S(I,NRECS)
          END DO

          IF(LOFFSET) THEN         ! Rad or conv timesteps,
                                   !       1 alresdy added
            LIST_S(st_start_time_code,NRECS+1)=                         &
     &      LIST_S(st_start_time_code,NRECS+1)-1
            IF (LIST_S(st_period_code,NRECS) /= -1) THEN
! Offsets are added to start time
             LIST_S(st_offset_code,NRECS)=                              &
! DEPENDS ON: totimp
     &           TOTIMP(IOFF_T(ITIM_L),UNT2_T(ITIM_L),MODL_L)
             if (LIST_S(st_offset_code,NRECS)  ==  -999) then
               ErrorStatus = 1
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               write(6,*) cmessage
               GOTO 9999
              endif
              LIST_S(st_start_time_code,NRECS)=                         &
     &        LIST_S(st_start_time_code,NRECS)-                         &
     &        LIST_S(st_period_code,NRECS) +                            &
     &        LIST_S(st_offset_code,NRECS)
            ELSE
              LIST_S(st_start_time_code,NRECS)=1
            END IF

          ELSE

            IF (LIST_S(st_period_code,NRECS) /= -1) THEN
! Offsets are added to start time
             LIST_S(st_offset_code,NRECS)=                              &
! DEPENDS ON: totimp
     &           TOTIMP(IOFF_T(ITIM_L),UNT2_T(ITIM_L),MODL_L)
             if (LIST_S(st_offset_code,NRECS)  ==  -999) then
               ErrorStatus = 1
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               write(6,*) cmessage
               GOTO 9999
              endif
              LIST_S(st_start_time_code,NRECS)=                         &
     &        LIST_S(st_start_time_code,NRECS)-                         &
     &        LIST_S(st_period_code,NRECS)+1+                           &
     &        LIST_S(st_offset_code,NRECS)
            ELSE
              LIST_S(st_start_time_code,NRECS)=1
            END IF

          END IF

          IF(LIST_S(st_start_time_code,NRECS) <  1) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'DIAGNOSTIC START TIME BEFORE PERIOD, SETTING TO 1.',        &
     &     '(M,S,I) ',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-170
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            LIST_S(st_start_time_code,NRECS)=1
          END IF

! Check if offset corresponds to a valid timestep.
! If not, reject the diagnostic.
          SELECT CASE (ITIMA)
            CASE (2)   ! Long-Wave Radiation
              IF (MOD(LIST_S(st_offset_code,NRECS),A_LW_RADSTEP)        &
     &            /= 0)                                                 &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH LW_RAD STEP.',          &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-180
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF

            CASE (3)   ! Short-Wave Radiation
              IF (MOD(LIST_S(st_offset_code,NRECS),A_SW_RADSTEP)        &
     &            /= 0)                                                 &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH SW_RAD STEP.',          &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-190
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF
            CASE (13)   ! Convection
              IF (MOD(LIST_S(st_offset_code,NRECS),A_CONV_STEP)/= 0)    &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH CONVECT STEP.',         &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-200
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF

            CASE (14)   ! Leaf Phenology
              IF (MOD(LIST_S(st_offset_code,NRECS),PHENOL_PERIOD)/= 0)  &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH PHENOL PERIOD.',        &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-210
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF

            CASE (15)   ! Triffid
              IF (MOD(LIST_S(st_offset_code,NRECS),TRIFFID_PERIOD)/= 0) &
     &        THEN
                WRITE(CMESSAGE2,*)                                      &
     &              'OFFSET DOES NOT AGREE WITH TRIFFID PERIOD.',       &
     &              ' DIAG IGNORED. (M,S,I)',                           &
     &               MODL_L,ISEC_L,ITEM_L

                ERRORSTATUS=-220
! DEPENDS ON: ereport
                CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
                NRECS=NRECS-1
                GOTO 999
              END IF
          END SELECT

          LIST_S(st_proc_no_code ,NRECS+1)=1

          LIST_S(st_input_bottom ,NRECS+1)=                             &
     &    LIST_S(st_output_bottom,NRECS  )

          LIST_S(st_input_top    ,NRECS+1)=                             &
     &    LIST_S(st_output_top   ,NRECS  )

          LIST_S(st_input_code   ,NRECS+1)=-NRECS
          LIST_S(st_output_code  ,NRECS  )=1
          LIST_S(st_series_ptr   ,NRECS+1)=0
          LIST_S(NELEMP+1        ,NRECS+1)=NRECS+1

          LIST_S(st_freq_code,NRECS)=                                   &
                                                    ! Frequency
! DEPENDS ON: totimp
     &    TOTIMP(ISAM_T(ITIM_L),UNT2_T(ITIM_L),MODL_L)
          if (LIST_S(st_freq_code,NRECS)  ==  -999) then
             ErrorStatus = 105
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               write(6,*) cmessage
             GOTO 9999
          endif

!   Correct frequency for radiation, periodic convection, leaf
!   phenology and vegetation competition

          IF (ITIMA == 2) THEN
            IF (LIST_S(st_freq_code,NRECS) == 1) THEN
              LIST_S(st_freq_code,NRECS)=A_LW_RADSTEP
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),A_LW_RADSTEP) /= 0) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'INCORRECT SAMPLING FOR LW_RADSTEP. FREQ=',                &
     &        LIST_S(st_freq_code,NRECS), ':IGNORED.',                  &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-225
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          ELSE IF(ITIMA == 3) THEN
            IF (LIST_S(st_freq_code,NRECS) == 1) THEN
              LIST_S(st_freq_code,NRECS)=A_SW_RADSTEP
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),A_SW_RADSTEP) /= 0) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'INCORRECT SAMPLING FOR SW_RADSTEP. FREQ=',                &
     &        LIST_S(st_freq_code,NRECS), ':IGNORED.',                  &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-230
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          ELSE IF(ITIMA == 13) THEN
            IF (LIST_S(st_freq_code,NRECS) == 1) THEN
              LIST_S(st_freq_code,NRECS)=A_CONV_STEP
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),A_CONV_STEP) /= 0) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'INCORRECT SAMPLING FOR CONV_STEP. FREQ=',                 &
     &        LIST_S(st_freq_code,NRECS), ':IGNORED.',                  &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-240
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          ELSE IF(ITIMA == 14) THEN
            IF (LIST_S(st_freq_code,NRECS) == 1) THEN
              LIST_S(st_freq_code,NRECS)=PHENOL_PERIOD
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),PHENOL_PERIOD) /= 0) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'INCORRECT SAMPLING FOR PHENOL_PERIOD. FREQ=',             &
     &        LIST_S(st_freq_code,NRECS), ':IGNORED.',                  &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-250
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          ELSE IF(ITIMA == 15) THEN
            IF (LIST_S(st_freq_code,NRECS) == 1) THEN
              LIST_S(st_freq_code,NRECS)=TRIFFID_PERIOD
            ELSE IF                                                     &
     &      (MOD(LIST_S(st_freq_code,NRECS),TRIFFID_PERIOD) /= 0) THEN
              WRITE(CMESSAGE2,*)                                        &
     &       'INCORRECT SAMPLING FOR TRIFFID_PERIOD. FREQ=',            &
     &        LIST_S(st_freq_code,NRECS), ':IGNORED.',                  &
     &       '(M,S,I) ',                                                &
     &        MODL_L,ISEC_L,ITEM_L

              ERRORSTATUS=-260
! DEPENDS ON: ereport
              CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
              NRECS=NRECS-1
              GOTO 999
            END IF
          END IF

    ! For the NRECS item an end_time_code needs to be set if we
    ! are dealing with a times table rather than  regular diagn.
    ! This should be the maximum timestep in the time list. The list
    ! should be ready sorted (and thus maximum is last member) but
    ! will run through and find maximum to be on the safe side.

          IMAX = 0
          ITIMLST = -LIST_S(st_freq_code,NRECS+1)

          IF (ITIMLST  >   0 .and. ITIMLST /= 9999 ) THEN      ! List *not* regular
            DO I = 1, ITIM_T(ITIMLST)
              IF (IMAX  <   ITIM_S(I, ITIMLST)) THEN
                IMAX = ITIM_S(I, ITIMLST)
              END IF
            END DO

            LIST_S(st_end_time_code,NRECS) = IMAX

          END IF

!   Period

          IF ((INTV_T(ITIM_L) == -1).AND.(ITYP_T(ITIM_L) == 2)) THEN
            LIST_S(st_period_code,NRECS)=-1
          ELSE
            LIST_S(st_period_code,NRECS)=                               &
! DEPENDS ON: totimp
     &      TOTIMP(INTV_T(ITIM_L),UNT1_T(ITIM_L),MODL_L)
            if (LIST_S(st_period_code,NRECS)  ==  -999) then
               ErrorStatus = 106
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               write(6,*) cmessage
               GOTO 9999
            endif
          END IF

!   Add the record - unless the output destination is the dump,
!                      and output at the accumulating period
          IF (    LOCN_U(IUSE_L) >  2                                   &
     &      .OR.                                                        &
     &       ( (LIST_S(st_freq_code  ,NRECS+1) /=                       &
     &          LIST_S(st_period_code,NRECS  ))                         &
     &                                        .AND.                     &
     &         (LIST_S(st_start_time_code,NRECS+1) /=                   &
     &          LIST_S(st_end_time_code  ,NRECS+1))   )                 &
     &        )THEN
! No tag for parent, except for dump frequency means
             IF (UNT1_T(ITIM_L)=="DU" .AND. UNT1_T(ITIM_L)=="DU") THEN
                ! Special tag for dump frequency means
                LIST_S(st_macrotag,NRECS)=-999
             ELSE
                LIST_S(st_macrotag,NRECS)=0
             ENDIF
            NRECS=NRECS+1
          END IF

        ELSE IF (LIST_S(st_proc_no_code,NRECS) == 8) THEN
! Option of "daily" mean timeseries

          IF(NRECS >= NRECDP) THEN
            WRITE(CMESSAGE2,*)                                          &
     &     'TOO MANY S_LIST REQUESTS. REQUEST IGNORED',                 &
     &     '(M,S,I) ',                                                  &
     &      MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-270
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            NRECS=NRECS-1
            GOTO 999
          END IF

! Special case where 2 extra records required
!  Record 1 - time mean only no spatial processing
!  Record 2 - timeseries formed extracting from record 1
!  Record 3 - extract timeseries from dump ie record 2

          DO I=1,NELEMP+1          ! Copy stash list forward
            LIST_S(I,NRECS+1)=LIST_S(I,NRECS)
            LIST_S(I,NRECS+2)=LIST_S(I,NRECS)
          END DO

          IF(LOFFSET) THEN         ! Rad or conv timesteps,
                                   !       1 already added
            LIST_S(st_start_time_code,NRECS+2)=                         &
     &      LIST_S(st_start_time_code,NRECS+2)-1
            IF (LIST_S(st_period_code,NRECS) /= -1) THEN
              LIST_S(st_start_time_code,NRECS)=                         &
     &        LIST_S(st_start_time_code,NRECS)-                         &
     &        LIST_S(st_period_code,NRECS)
            ELSE
              LIST_S(st_start_time_code,NRECS)=1
            END IF

          ELSE

            IF (LIST_S(st_period_code,NRECS) /= -1) THEN
              LIST_S(st_start_time_code,NRECS)=                         &
     &        LIST_S(st_start_time_code,NRECS)-                         &
     &        LIST_S(st_period_code,NRECS)+1
            ELSE
              LIST_S(st_start_time_code,NRECS)=1
            END IF

          END IF

          IF(LIST_S(st_start_time_code,NRECS) <  1) THEN
            WRITE(CMESSAGE2,*)                                          &
     &       'START TIME BEFORE PERIOD, SETTING TO 1',                  &
     &       '(M,S,I) ',MODL_L,ISEC_L,ITEM_L

            ERRORSTATUS=-280
! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ERRORSTATUS, CMESSAGE2)
            LIST_S(st_start_time_code,NRECS)=1
          END IF

          LIST_S(st_proc_no_code ,NRECS)=3    ! time mean
          LIST_S(st_proc_no_code ,NRECS+1)=8  ! timseries special case
          LIST_S(st_proc_no_code ,NRECS+2)=1  !  extract

! Reset first record to no area weight or spatial processing
! ie first record just controls time meaning

          LIST_S(st_gridpoint_code,NRECS)=1
          LIST_S(st_weight_code,NRECS)=0

          LIST_S(st_input_bottom ,NRECS+1)=                             &
     &    LIST_S(st_output_bottom,NRECS  )
          LIST_S(st_input_bottom ,NRECS+2)=                             &
     &    LIST_S(st_output_bottom,NRECS+1)

          LIST_S(st_input_top    ,NRECS+1)=                             &
     &    LIST_S(st_output_top   ,NRECS  )
          LIST_S(st_input_top    ,NRECS+2)=                             &
     &    LIST_S(st_output_top   ,NRECS+1)

          LIST_S(st_input_code   ,NRECS+1)=-NRECS
          LIST_S(st_input_code   ,NRECS+2)=-NRECS-1
          LIST_S(st_output_code  ,NRECS  )=1
          LIST_S(st_output_code  ,NRECS+1)=1  ! dump
          LIST_S(st_series_ptr   ,NRECS+2)=0
          LIST_S(st_series_ptr   ,NRECS)=0
          LIST_S(NELEMP+1        ,NRECS+1)=NRECS+1
          LIST_S(NELEMP+1        ,NRECS+2)=NRECS+2
!  definition 8 implies frequency of time mean over every timestep
          LIST_S(st_freq_code,NRECS)=1
          LIST_S(st_freq_code,NRECS+1)=                                 &
                                                      ! Frequency
! DEPENDS ON: totimp
     &    TOTIMP(ISAM_T(ITIM_L),UNT2_T(ITIM_L),MODL_L)
          IF (LIST_S(st_freq_code,NRECS)  ==  -999) THEN
             ErrorStatus = 107
             write (cmessage,'(a,a,i2,a,i3,a,i3)')                      &
     &       'PRELIM:TOTIMP:Error in time period conversion',           &
     &       ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
             GOTO 9999
          ENDIF


!   Correct frequency for radiation, periodic convection, leaf
!   phenology and vegetation competition

          IF (ITIMA == 2) THEN
              LIST_S(st_freq_code,NRECS)=A_LW_RADSTEP
          ELSE IF(ITIMA == 3) THEN
              LIST_S(st_freq_code,NRECS)=A_SW_RADSTEP
          ELSE IF(ITIMA == 13) THEN
              LIST_S(st_freq_code,NRECS)=A_CONV_STEP
          ELSE IF(ITIMA == 14) THEN
              LIST_S(st_freq_code,NRECS)=PHENOL_PERIOD
          ELSE IF(ITIMA == 15) THEN
              LIST_S(st_freq_code,NRECS)=TRIFFID_PERIOD
          END IF

!   Period
! time mean over sampling period
            LIST_S(st_period_code,NRECS)=                               &
! DEPENDS ON: totimp
     &      TOTIMP(ISAM_T(ITIM_L),UNT2_T(ITIM_L),MODL_L)
            IF (LIST_S(st_period_code,NRECS)  ==  -999) THEN
               ErrorStatus = 108
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               GOTO 9999
            ENDIF

! period for timeseries recycle period
            LIST_S(st_period_code,NRECS+1)=                             &
! DEPENDS ON: totimp
     &      TOTIMP(INTV_T(ITIM_L),UNT1_T(ITIM_L),MODL_L)
            IF (LIST_S(st_period_code,NRECS+1)  ==  -999) THEN
               ErrorStatus = 109
               write (cmessage,'(a,a,i2,a,i3,a,i3)')                    &
     &         'PRELIM:TOTIMP:Error in time period conversion',         &
     &         ' model=',MODL_L,' section=',ISEC_L,' item=',ITEM_L
               GOTO 9999
            ENDIF


! st_start_time for 2 record should be period for first record
! unless offset from start of run. Note value independent of logical
!  OFFSET
!
            IF (LIST_S(st_period_code,NRECS) /= -1) THEN
              IF (LOFFSET) THEN
               LIST_S(st_start_time_code,NRECS+1)=                      &
     &          LIST_S(st_start_time_code,NRECS+1) -                    &
     &          LIST_S(st_period_code,NRECS+1) +                        &
     &          LIST_S(st_freq_code,NRECS+1) - 1
              ELSE
               LIST_S(st_start_time_code,NRECS+1)=                      &
     &          LIST_S(st_start_time_code,NRECS+1) -                    &
     &          LIST_S(st_period_code,NRECS+1) +                        &
     &          LIST_S(st_freq_code,NRECS+1)
              ENDIF
            ELSE
              LIST_S(st_start_time_code,NRECS+1)=1
            END IF


!   Add both record
            LIST_S(st_macrotag,NRECS)=0
            NRECS=NRECS+2

        END IF       ! Other than single time field

      END IF         ! Diag request not null - ITIM_L /= 0
 999  CONTINUE
      END DO         ! Loop over diagnostic requests

      END IF         ! NDIAG >  0

! DEPENDS ON: pslcom
      CALL PSLCOM(NRECS)    ! Compress out unused pseudo levels lists

 9999 RETURN
      END SUBROUTINE PRELIM

!- End of Subroutine code -------------------------------------------


!+ Compress out unused pseudo levels lists

!- End of Subroutine code ------------------------------------------
