
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Read the basis file information
! Subroutine Interface:

      SUBROUTINE RDBASIS                                                &
     & (IU,CMESSAGE,ErrorStatus)
      IMPLICIT NONE

! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!   4.0     Sept.95    Original code.  S.J.Swarbrick
!  4.0  18/10/95  Add ErrorStatus to GET_FILE call. RTHBarnes
!   4.1     June 96    DISCT_LEV function used to check for model
!                         levels       S.J.Swarbrick
!   4.4     Sep. 97    Add IOFF to namelist for offset to sampling
!                      frequency. S.D. Mullerworth
!   4.5     July 98    Remove call to INTRFACE (A Van der Wal)
!   4.5    30/10/97    Read stash data on PE 0 for the T3E
!                      and distribute it.
!                        Author: Bob Carruthers, Cray Research
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!    6.0  11/09/03  Add file='empty' for IBM (provided by Zoe Chaplin)
!                                                           P.Dando
!    6.2  23/11/05  Removed DIAG190 define and diagnostic. T. Edwards
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:

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
! COMDECK LENFIL
! Description:
! Defines character string length used for file names in STASH request
! processing routines
!
! Current code owner: S.J.Swarbrick
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 3.5       Apr. 95   Original code.  S.J.Swarbrick
! 4.0       Oct. 95                   S.J.Swarbrick
! 4.1       Apr. 96   Tidy            S.J.Swarbrick
!
      CHARACTER*80 FILE
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
!========================== COMDECK COCNINDX ==========================
!
! Description:
!
!       This comdeck contains all the indices and row-wise loop
!       control variables required by the ocean mpp code.
!
      ! Note: All variables prefixed "J_" contain values which
      ! take account of halo sizes. Eg: for the 3 row domain defined
      ! by JST = 10 and JFIN = 12, with a halo of 2 rows, J_1
      ! will be 3, J_JMT will be 5.

      INTEGER :: J_1     ! Local value of loop control for J = 1, n
      INTEGER :: J_2     !   "     "    "   "     "     "  J = 2, n
      INTEGER :: J_3     !   "     "    "   "     "     "  J = 3, n
      INTEGER :: J_JMT   !   "     "    "   "     "     "  J = n, JMT
      INTEGER :: J_JMTM1 !   "     "    "   "     "     "  J = n, JMTM1
      INTEGER :: J_JMTM2 !   "     "    "   "     "     "  J = n, JMTM2
      INTEGER :: J_JMTP1 !   "     "    "   "     "     "  J = n, JMTP1
      INTEGER :: JMT_NOHALO  ! New variable needed for SWAP_BOUNDS!
      INTEGER :: JMTM1_NOHALO
      INTEGER :: JST     ! First row this process considers (no halo)
      INTEGER :: JFIN    ! Last   "    "     "     "        "
      INTEGER :: J_FROM_LOC       ! Local value of start index
      INTEGER :: J_TO_LOC         ! Local value of end index
      INTEGER :: JMT_GLOBAL       ! Global value of JMT
      INTEGER :: JMTM1_GLOBAL     ! Global value of JMT - 1
      INTEGER :: JMTM2_GLOBAL     ! Global value of JMT - 2
      INTEGER :: JMTP1_GLOBAL     ! Global value of JMT + 1
      INTEGER :: J_OFFSET         ! Start row - 1
      INTEGER :: O_MYPE           ! MYPE value in arg lists for ocean
      INTEGER :: O_EW_HALO        ! EW_HALO for ocean arg lists
      INTEGER :: O_NS_HALO        ! NS_HALO for ocean arg lists
      INTEGER :: J_PE_JSTM1
      INTEGER :: J_PE_JSTM2
      INTEGER :: J_PE_JFINP1
      INTEGER :: J_PE_JFINP2
      INTEGER :: O_NPROC          ! NPROC for ocean
      INTEGER :: imout(4),jmout(4)! i,j indices for pts in Med outflow
      INTEGER :: J_PE_IND_MED(4)  ! no for each PE in Med outflow
      INTEGER :: NMEDLEV          ! no of levels for Med outflow
      INTEGER :: lev_med  ! level at which deep flow from Med occurs
      INTEGER :: lev_hud  ! level at which deep flow into Hudson Bay
      INTEGER :: imout_hud(4),jmout_hud(4)  ! Hudson Bay i,j
      INTEGER :: J_PE_IND_HUD(4)  ! PE's involved in Hudson Bay outflow

      ! last level for which there is inflow to Mediterranean
      INTEGER :: med_topflow
      INTEGER :: JMT_MAX ! Highest value of JMT on any PE
      REAL :: LAT_CHECK  ! Save NSSPACEO for consistency checks.

      COMMON /COCNINDX/                                                 &
     &  J_1, J_2, J_3,                                                  &
     &  J_JMT, J_JMTM1, J_JMTM2, J_JMTP1,                               &
     &  JST, JFIN, JMT_NOHALO,JMTM1_NOHALO, J_FROM_LOC, J_TO_LOC,       &
     &  JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL,                         &
     &  JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO,           &
     &  J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2,               &
     &  O_NPROC,JMT_MAX,                                                &
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,                               &
     &  lev_med,lev_hud,imout_hud,jmout_hud,J_PE_IND_HUD,med_topflow    &
     &  ,LAT_CHECK

! COCNINDX end

! Subroutine arguments

!   Scalar arguments with intent(in):
      INTEGER IU       ! Unit no. of stash basis file
      INTEGER IE
!   Scalar arguments with intent(out):
      CHARACTER*80 CMESSAGE      ! Error return message

!   Error status:
      INTEGER        ErrorStatus ! Error return code

! Local variables:
      LOGICAL MODEL_LEV  !TRUE for model levels
      INTEGER I,J,L
      INTEGER IDUM
      INTEGER IOSTAT
      INTEGER NtsRecs    !Counter for time series records

!   Namelist STASHNUM
      INTEGER            NUM_REQ,NUM_DOM,NUM_TIM,NUM_USE
      NAMELIST/STASHNUM/ NUM_REQ,NUM_DOM,NUM_TIM,NUM_USE

!   Namelist STREQ: STASH requests
      INTEGER         IMOD,ISEC,ITEM,IDOM,ITIM,IUSE
      NAMELIST/STREQ/ IMOD,ISEC,ITEM,IDOM,ITIM,IUSE

!   Namelist TIME: Time profiles
      CHARACTER*8 NAME
      CHARACTER*2 UNT1,UNT2,UNT3
      INTEGER     ITYP,ISAM,INTV,IOPT
      INTEGER     ISTR,IEND,IFRE,IOFF,ITIMES,ISER(NTIMEP)
      NAMELIST/TIME/ITYP,ISAM,INTV,UNT1  ,UNT2,UNT3,IOPT                &
     &             ,ISTR,IEND,IFRE,IOFF,ITIMES,ISER,NAME

!   Namelist DOMAIN: Domain profiles
      INTEGER     IOPL           !Level type code
      INTEGER     ILEVS          !Flag for range/selected model levels
      INTEGER     LEVB,LEVT      !Bottom/top levels for range
      INTEGER     PLT            !Pseudo level type code
      INTEGER     IOPA                !Horizontal domain type code
      INTEGER     INTH,ISTH,IEST,IWST !Horiz domain limits (IOPA=9,10)
      INTEGER     IMSK                !Grid type code
      INTEGER     IMN                 !Spatial meaning code
      INTEGER     IWT                 !Weighting code
      INTEGER     LEVLST (NLEVP)      !Levels lists array: integer
      REAL        RLEVLST(NLEVP)      !                    real
      INTEGER     PSLIST (NPSLEVP)    !                    pseudo
      CHARACTER*1 TS                    !Flag for time series domain
      INTEGER     TSNUM                 !No. of time ser doms in prof
      INTEGER     TBLIM (tsdp),TTLIM (tsdp) !TS dom limits (top/bottom)
      REAL        TBLIMR(tsdp),TTLIMR(tsdp) !ditto for real levels
      INTEGER     TNLIM (tsdp),TSLIM (tsdp) !TS dom limits (N/S)
      INTEGER     TWLIM (tsdp),TELIM (tsdp) !TS dom limits (E/W)

      NAMELIST/DOMAIN/IOPL ,ILEVS ,LEVB  ,LEVT  ,PLT    ,IOPA ,IMSK ,   &
     &                IMN  ,IWT   ,TS    ,LEVLST,RLEVLST,NAME ,         &
     &                INTH ,ISTH  ,IEST  ,IWST  ,PSLIST ,               &
     &                TSNUM,TBLIM ,TTLIM ,TNLIM ,TSLIM  ,TELIM,TWLIM,   &
     &                      TBLIMR,TTLIMR

!   Namelist USE: Useage profiles
      INTEGER      LOCN,IUNT
      NAMELIST/USE/LOCN,IUNT,NAME

! Function and subroutine calls:
      LOGICAL  DISCT_LEV
      EXTERNAL GET_FILE

!- End of Header ------------------------------------------------------

! Initialisation
      NDIAG   =0
      NTPROF  =0
      NDPROF  =0
      NUPROF  =0
      NUM_REQ =0
      NUM_TIM =0
      NUM_DOM =0
      NUM_USE =0
      NtsRecs =0

      DO I = 1,NDIAGPM
        MODL_B(I)=0
        ISEC_B(I)=0
        ITEM_B(I)=0
        ITIM_B(I)=0
        IDOM_B(I)=0
        IUSE_B(I)=0
      END DO

        ITIMES   =0
      DO I=1,NPROFTP
        TIMPRO(I)='        '
        ITYP_T(I)=0
        INTV_T(I)=0
        UNT1_T(I)='  '
        ISAM_T(I)=0
        UNT2_T(I)='  '
        IOPT_T(I)=0
        ISTR_T(I)=0
        IEND_T(I)=0
        IFRE_T(I)=0
        IOFF_T(I)=0
        UNT3_T(I)='  '
        ITIM_T(I)=0
        MODL_T(I)=0
        DO J = 1,NTIMEP
          ISER_T(J,I)=0
        END DO
      END DO

      DO I=1,NPROFDP
        DOMPRO  (I)='        '
        IOPL_D  (I)=0
        LEVB_D  (I)=0
        LEVT_D  (I)=0
        PLT_D   (I)=0
        IOPA_D  (I)=0
        INTH_D  (I)=0
        ISTH_D  (I)=0
        IEST_D  (I)=0
        IWST_D  (I)=0
        IMSK_D  (I)=0
        IMN_D   (I)=0
        IWT_D   (I)=0
        TS_D    (I)=' '
        PLLEN_D (I)=0
        PLPOS_D (I)=0
        ILEV_D  (I)=0
      END DO
      DO I = 1,NLEVP
      DO J = 1,NPROFDP
        LEVLST_D (I,J)=0
        RLEVLST_D(I,J)=0
      END DO
      END DO
      DO I = 1,NPSLEVP
      DO J = 1,NPSLISTP
        PSLIST_D(I,J)=0
      END DO
      END DO

      DO I = 1,NPROFUP
        USEPRO(I)='        '
        LOCN_U(I)=0
        IUNT_U(I)=0
      END DO

      CALL GET_FILE(IU,FILE,80,ErrorStatus)   ! Get name for stash file

! Rewind stash control file
      REWIND(IU)

!STASH control file header namelist
  99  READ(IU,STASHNUM)
      IF (NUM_REQ == -1) GOTO 999
      NDIAG =NDIAG +NUM_REQ
      NTPROF=NTPROF+NUM_TIM
      NDPROF=NDPROF+NUM_DOM
      NUPROF=NUPROF+NUM_USE
      IF (NDIAG  >  NRECDP ) THEN
        WRITE(6,*) 'NUMBER OF DIAGNOSTIC REQUESTS EXCEEDS LIMIT OF ',   &
     &  NRECDP ,' SOME HAVE BEEN IGNORED'
        GO TO 999
      END IF
      IF (NDPROF >  NPROFDP) THEN
        WRITE(6,*) 'ERROR IN STASHC:'
        WRITE(6,*) 'NUMBER OF DOMAIN PROFILES EXCEEDS LIMIT OF ',       &
     &  NPROFDP
        WRITE(6,*) 'ARRAYS WILL BE OVERWRITTEN'
        ErrorStatus=1
        GO TO 9999
      END IF
      IF (NUPROF >  NPROFUP) THEN
        WRITE(6,*) 'ERROR IN STASHC:'
        WRITE(6,*) 'NUMBER OF USEAGE PROFILES EXCEEDS LIMIT OF ',       &
     &  NPROFUP
        WRITE(6,*) 'ARRAYS WILL BE OVERWRITTEN'
        ErrorStatus=1
        GO TO 9999
      END IF
      IF (NTPROF >  NPROFTP) THEN
        WRITE(6,*) 'ERROR IN STASHC:'
        WRITE(6,*) 'NUMBER OF TIME PROFILES EXCEEDS LIMIT OF ',         &
     &  NPROFTP
        WRITE(6,*) 'ARRAYS WILL BE OVERWRITTEN'
        ErrorStatus=1
        GO TO 9999
      END IF

!STASH request namelists
      IF (NUM_REQ >  0) THEN
      DO I = NDIAG-NUM_REQ+1,NDIAG
        IMOD=0
        ISEC=0
        ITEM=0
        ITIM=0
        IDOM=0
        IUSE=0
        READ(IU,STREQ)
        MODL_B(I)=IMOD
        ISEC_B(I)=ISEC
        ITEM_B(I)=ITEM
        IF (ITIM /= 0) THEN
          ITIM     =ITIM+NTPROF-NUM_TIM
        END IF
        IDOM     =IDOM+NDPROF-NUM_DOM
        IUSE     =IUSE+NUPROF-NUM_USE
        ITIM_B(I)=ITIM
        IDOM_B(I)=IDOM
        IUSE_B(I)=IUSE
      END DO
      END IF

!Time profile namelists
      IF (NUM_TIM >  0) THEN
      DO I = NTPROF-NUM_TIM+1,NTPROF
        NAME='        '
        ISAM=0
        INTV=0
        IOPT=0
        ISTR=0
        IEND=0
        IFRE=0
        IOFF=0
        DO J = 1,NTIMEP
          ISER(J)=0
        END DO
!  Read namelist
        READ(IU,TIME)
        TIMPRO(I)=NAME
        ITYP_T(I)=ITYP
        IF (ITYP /= 1) THEN
!  Diagnostic output is time-processed
          ISAM_T(I)=ISAM  !Sampling frequency
          UNT2_T(I)=UNT2
          INTV_T(I)=INTV  !Processing interval
          UNT1_T(I)=UNT1
        END IF
!  Diag. output time option
        IOPT_T(I)=IOPT
        UNT3_T(I)=UNT3
        IF      (IOPT == 1) THEN
!    Regular output time interval
          ISTR_T(I)=ISTR
          IEND_T(I)=IEND
          IFRE_T(I)=IFRE
          IOFF_T(I)=IOFF
        ELSE IF (IOPT == 2) THEN
!    Specified list of output times
!     Length of times table
          ITIM_T(I)=ITIMES
!     Internal model label for times table
          IF (NDIAG >  0) THEN
            MODL_T(I)=MODL_B(NDIAG)
          END IF
!     Times table
          DO J = 1,ITIMES
            ISER_T(J,I)=ISER(J)
          END DO
        END IF
      END DO
      END IF

!Domain profile namelists
      IF (NUM_DOM >  0) THEN
      DO I = NDPROF-NUM_DOM+1,NDPROF
!Initialise
        NAME ='        '
        IOPL =0
        LEVB =0
        LEVT =0
        ILEVS=0
          LEVLST (1)= IMDI
         RLEVLST (1)= RMDI
        DO J = 2,NLEVP
          LEVLST (J)= 0
          RLEVLST(J)= 0.0
        END DO
        DO J = 1,NPSLEVP
          PSLIST(J)=0
        END DO
        DO J = 1,tsdp
          TBLIM (J)=0
          TTLIM (J)=0
          TBLIMR(J)=0.
          TTLIMR(J)=0.
          TNLIM (J)=0
          TSLIM (J)=0
          TELIM (J)=0
          TWLIM (J)=0
        END DO
!Read namelist
        READ(IU,DOMAIN)
!Check for errors in levels lists
! DEPENDS ON: disct_lev
        MODEL_LEV=DISCT_LEV(IOPL,ErrorStatus,CMESSAGE)
        IF (MODEL_LEV) THEN
! Model levels
          IF (ILEVS == 2) THEN
            IF ( LEVLST(1) == IMDI ) THEN
              WRITE(6,*)                                                &
     &       'ERROR,RDBASIS: LEVELS LIST IN DOMAIN PROFILE '            &
     &       ,I,' HAS NO ENTRIES'
              CMESSAGE='ERROR,RDBASIS: LEVELS LIST HAS NO ENTRIES'
              ErrorStatus=1
              GO TO 9999
            END IF
          END IF
        ELSE IF (IOPL /= 5) THEN
          IF (RLEVLST(1) == RMDI) THEN
            WRITE(6,*)                                                  &
     &     'ERROR,RDBASIS: LEVELS LIST IN DOMAIN PROFILE '              &
     &     ,I,' HAS NO ENTRIES'
            CMESSAGE='ERROR,RDBASIS: LEVELS LIST HAS NO ENTRIES'
            ErrorStatus=1
            GO TO 9999
          END IF
        END IF
!Profile name
        DOMPRO(I)=NAME
!Store level type code in IOPL_D
        IOPL_D(I)=IOPL
! DEPENDS ON: disct_lev
        MODEL_LEV=DISCT_LEV(IOPL,ErrorStatus,CMESSAGE)
        IF (MODEL_LEV) THEN
!Integer levels
          ILEV_D(I)=ILEVS
          IF (ILEVS == 1) THEN
!  Range of model levels
            LEVB_D(I)=LEVB
            LEVT_D(I)=LEVT
          END IF
          IF (ILEVS == 2) THEN
!  List of selected model levels
            LEVB_D(I)=-1
            DO J=1,NLEVP
              LEVLST_D(J,I) = LEVLST(J)
              IF (LEVLST(J) >  0) THEN
!  Store no. of levels in LEVT_D(I)
                LEVT_D(I)=LEVT_D(I)+1
              END IF
            END DO
          END IF
        ELSE IF (IOPL /= 5) THEN
!Real levels
          LEVB_D(I)=-1
          DO J=1,NLEVP
            RLEVLST_D(J,I) = RLEVLST(J)
              IF (RLEVLST(J) >  0.0) THEN
!  Store no. of levels in LEVT_D(I)
                LEVT_D(I)=LEVT_D(I)+1
              END IF
          END DO
        END IF
!Store pseudo level type code in PLT_D
        PLT_D (I)=PLT
        IF (PLT >  0) THEN
!Domain profile 'I' has pseudo levels list
!    Count total no. of pseudo levels lists
        NPSLISTS = NPSLISTS + 1
!    Store list in column 'NPSLISTS' of PSLIST_D
          DO J=1,NPSLEVP
            PSLIST_D(J,NPSLISTS) = PSLIST (J)
!    PPLEN_D(I) stores length of ps.lev.list for domain 'I'
            IF (PSLIST(J) >  0) THEN
              PLLEN_D (I)        = PLLEN_D(I) + 1
            END IF
          END DO
!    PLPOS(I) stores the column no. in PSLIST_D for dom. prof. 'I'
          PLPOS_D(I) = NPSLISTS
        END IF
!Store horizontal domain type in IOPA_D
        IOPA_D(I)=IOPA
        IF (IOPA == 9.OR.IOPA == 10) THEN
!    Specified area
          INTH_D(I)=INTH
          ISTH_D(I)=ISTH
          IEST_D(I)=IEST
          IWST_D(I)=IWST
        END IF
        IMSK_D(I)=IMSK ! Gridpoint option
        IMN_D (I)=IMN  ! Meaning option
        IWT_D (I)=IWT  ! Weighting option
        TS_D  (I)=TS   ! Time series domain
        IF (TS_D(I)  ==  'Y') THEN
!This domain profile has a block of time series domains
!  Store time series data for current profile in _TS arrays
          NSERIES    = NSERIES+1        ! Time series block number:
          NPOS_TS(I) = NSERIES          !      used as a pointer
          NRECS_TS(NSERIES) = TSNUM     ! No. of records in ts block
          NSERREC_S  = NSERREC_S+TSNUM  ! Cumulative total ts records
          IF (NSERREC_S <= NTimSerP) THEN
            DO J = 1,TSNUM
              IF (J <= tsdp) THEN
                NtsRecs = NtsRecs+1
! DEPENDS ON: disct_lev
                MODEL_LEV=DISCT_LEV(IOPL,ErrorStatus,CMESSAGE)
                IF (MODEL_LEV) THEN
                  BLIM_TS (NtsRecs)= TBLIM (J)
                  TLIM_TS (NtsRecs)= TTLIM (J)
                ELSE IF (IOPL /= 5) THEN
                  BLIMR_TS(NtsRecs)= TBLIMR(J)
                  TLIMR_TS(NtsRecs)= TTLIMR(J)
                END IF
                NLIM_TS(NtsRecs) = TNLIM(J)
                SLIM_TS(NtsRecs) = TSLIM(J)
                ELIM_TS(NtsRecs) = TELIM(J)
                WLIM_TS(NtsRecs) = TWLIM(J)
              ELSE
                WRITE(6,*)                                              &
     &         'MESSAGE FROM ROUTINE RDBASIS: ',                        &
     &         'no. of time series in domain profile ',I,               &
     &         ' exceeds allowed limit of ',tsdp,' some will be',       &
     &         ' ignored'
              END IF
            END DO
          ELSE
            WRITE(6,*)                                                  &
     &     'TIMSER: total no. of time series requested exceeds ',       &
     &     'allowed limit of ',NTimSerP,'; some will be ignored.'
          END IF
        ELSE
          NPOS_TS  (I)=-1
        END IF
      END DO
      END IF
      NSERBLK_S = NSERIES

!Useage profile namelists
      IF (NUM_USE >  0) THEN
        NAME='        '
        LOCN=0
        IUNT=0
      DO I = NUPROF-NUM_USE+1,NUPROF
        READ(IU,USE)
        USEPRO(I)=NAME
        LOCN_U(I)=LOCN
        IUNT_U(I)=IUNT
      END DO
      END IF
      GO TO 99

 999  CONTINUE

!Initialise model config. arrays before reading STSHCOMP
      DO I = 0,NSECTP
        ATMOS_SR(I)='  '
        OCEAN_SR(I)='  '
         SLAB_SR(I)='  '
         WAVE_SR(I)='  '
        INDEP_SR(I)='  '
      END DO
      READ(5,STSHCOMP)

      CLOSE(UNIT=IU,STATUS='KEEP',IOSTAT=IOSTAT)

 9999 RETURN
      END SUBROUTINE RDBASIS

!- End of subroutine code ------------------------------------------
