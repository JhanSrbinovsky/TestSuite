
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Convert lat/long specification to row/column specification
! Subroutine Interface:

      SUBROUTINE LLTORC(GRID_CODE,                                      &
     &                  NORTH_LAT_IN,SOUTH_LAT_IN,                      &
     &                  WEST_LONG_IN,EAST_LONG_IN,                      &
     &                  START_ROW,END_ROW,START_COL,END_COL             &



     &                 )

      IMPLICIT NONE

! Description:
!   Uses the gridpoint code, the lat/long spec of the required area,
!   and the model area sizes to calculate the model row/column
!   numbers for the area.
!   Fix for the ocean model, since N is at the bottom.
!   Called by PRELIM, INPUTL, ADDRES.
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:

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
! DECOMPTP comdeck
!
! Description
!
! Magic numbers indicating decomposition types.
! These numbers are used to index the arrays defined in the
! DECOMPDB comdeck, and are required as an argument to
! the CHANGE_DECOMPOSITION subroutine.
!
! Current code owner : P.Burton
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 4.2       19/08/96  Original code.   P.Burton
! 4.3       17/02/97  Added new ocean decomposition decomp_nowrap_ocean
!                     which does not contain extra wrap points at
!                     start and end of row.                  P.Burton
! 5.5       04/08/00  Modification for parallelisation of WAM
!                   Author:  Bob Carruthers, Cray UK Inc(D.Holmes-Bell)

! Magic Numbers indicating decomposition types

      INTEGER                                                           &
     &  max_decomps                                                     &
                               ! maximum number of decompositions
     &, decomp_unset                                                    &
                               ! no decomposition selected
     &, decomp_standard_atmos                                           &
                               ! standard 2D atmosphere
!                              ! decomposition
     &, decomp_standard_ocean                                           &
                               ! standard 1D ocean decomposition
     &, decomp_nowrap_ocean                                             &
                               ! 1D ocean without extra wrap-around
!                              ! points at ends of each row
     &, decomp_smexe                                                    &
     &, decomp_standard_wave   ! standard 1D WAM Wave Model
                               ! decomposition

      PARAMETER (                                                       &
     &  max_decomps=5                                                   &
     &, decomp_unset=-1                                                 &
     &, decomp_standard_atmos=1                                         &
     &, decomp_standard_ocean=2                                         &
     &, decomp_nowrap_ocean=3                                           &
     &, decomp_smexe=4                                                  &
     &, decomp_standard_wave=5)

! End of DECOMPTP comdeck
! DECOMPDB comdeck
!
! Description:
!
! DECOMPDB comdeck (Decomposition Database) contains information
! describing the various decompositions used by the mpp-UM
! The CHANGE_DECOMPOSITION subroutine can be used to select
! a particular decomposition (which copies the appropriate
! decomposition information into the PARVARS common block).
!
! Requires comdeck PARVARS to be *CALLed before it.
!
! Current code owner : P.Burton
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 4.2       19/08/96  Original code.   P.Burton
! 5.0       12/04/99  - added dimension Nfld_max to decomp_db_glsize
!                     - added dimension NHalo_max to decomp_db_halosize
!                     - added dimension NHalo_max to decomp_db_g_lasize
!                     - added dimension Nfld_max to
!                       decomp_db_g_lasize
!                     - replace blsizep/u variables by blsize with new
!                       Nfld_max dimension
!                                                            P.Burton
! 5.1       27/01/00  - changed g_pe_index to g_pe_index_EW
!                     - added g_pe_index_NS
!                                                  P.Burton
! 5.3       27/01/00  - Moved data statement into blkdata. A Van der Wal
! 5.3       14/09/01  Added model_domain variable.  P.Burton
! 5.5       30/01/03  Generalised datastart. P.Selwood.

! Common blocks containing information about each decomposition
! (For description of variables see the PARVARS comdeck)

      INTEGER :: decomp_db_bound(Ndim_max,max_decomps)                  &
     &, decomp_db_sb_model_domain(max_decomps)
      INTEGER :: decomp_db_glsize(Ndim_max,Nfld_max,max_decomps)
      INTEGER :: decomp_db_gridsize(Ndim_max,max_decomps)
      INTEGER :: decomp_db_g_lasize(Ndim_max,Nfld_max,                  &
     &                     NHalo_max,0:maxproc,max_decomps)
      INTEGER :: decomp_db_g_blsize(Ndim_max,Nfld_max,0:maxproc,        &
     &                          max_decomps)
      INTEGER :: decomp_db_g_datastart(Ndim_max,0:maxproc,max_decomps)
      INTEGER :: decomp_db_g_datastart_f(Ndim_max,Nfld_max,0:maxproc,   &
     &                                   max_decomps)
      INTEGER :: decomp_db_g_gridpos(Ndim_max,0:maxproc,max_decomps)
      INTEGER :: decomp_db_g_pe_index_EW(1-Max_Halo_Size:               &
     &  ROW_LENGTH_MAX+Max_Halo_Size, max_decomps)
      INTEGER :: decomp_db_g_pe_index_NS(1-Max_Halo_Size:               &
     &  ROWS_MAX+Max_Halo_Size, max_decomps)
      INTEGER :: decomp_db_halosize(Ndim_max,NHalo_max,max_decomps)
      INTEGER :: decomp_db_neighbour(4,max_decomps)
      INTEGER :: decomp_db_first_comp_pe(max_decomps)
      INTEGER :: decomp_db_last_comp_pe(max_decomps)
      INTEGER :: decomp_db_nproc(max_decomps)
      INTEGER :: decomp_db_gc_proc_row_group(max_decomps)
      INTEGER :: decomp_db_gc_proc_col_group(max_decomps)
      INTEGER :: decomp_db_gc_all_proc_group(max_decomps)

      ! indicates if a decomposition has been initialised

      LOGICAL :: decomp_db_set(max_decomps)

      COMMON /DECOMP_DATABASE/                                          &
     &  decomp_db_bound , decomp_db_sb_model_domain , decomp_db_glsize  &
     & ,decomp_db_g_lasize , decomp_db_gridsize,                        &
     &  decomp_db_g_blsize,                                             &
     &  decomp_db_g_datastart , decomp_db_g_datastart_f,                &
     &  decomp_db_g_gridpos,                                            &
     &  decomp_db_g_pe_index_EW,decomp_db_g_pe_index_NS,                &
     &  decomp_db_halosize , decomp_db_neighbour,                       &
     &  decomp_db_first_comp_pe , decomp_db_last_comp_pe,               &
     &  decomp_db_nproc,                                                &
     &  decomp_db_gc_proc_row_group , decomp_db_gc_proc_col_group,      &
     &  decomp_db_gc_all_proc_group,                                    &
     &  decomp_db_set

! End of DECOMPDB comdeck
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
! PARAMETERs defining the GRID_TYPE characteristics

! General
      INTEGER, PARAMETER :: gt_unset = -1
      ! Any value which is unset

! MODEL_TYPE
      INTEGER, PARAMETER :: gt_atmos = 1
      ! Atmosphere field

      INTEGER, PARAMETER :: gt_ocean = 2
      ! Ocean field

      INTEGER, PARAMETER :: gt_wave = 3
      ! Wave field

! CONTENT
      INTEGER, PARAMETER :: gt_thetamass = 1
      ! Contains theta or mass points

      INTEGER, PARAMETER :: gt_velocity = 2
      ! Contains velocity (B grid U,V) points

      INTEGER, PARAMETER :: gt_U_C = 3
      ! Contains U points on C grid

      INTEGER, PARAMETER :: gt_V_C = 4
      ! Contains V points on C grid

      INTEGER, PARAMETER :: gt_hybrid = 5
      ! Points on none of the above
      INTEGER, PARAMETER :: gt_river = 6
      ! River routing grid

! COVERAGE
      INTEGER, PARAMETER :: gt_allpts = 1
      ! All points

      INTEGER, PARAMETER :: gt_land = 2
      ! Land points

      INTEGER, PARAMETER :: gt_sea = 3
      ! Sea points

! DOMAIN
      INTEGER, PARAMETER :: gt_full = 1
      ! Full field

      INTEGER, PARAMETER :: gt_zonal = 2
      ! Zonal field

      INTEGER, PARAMETER :: gt_meridional = 3
      ! Meridional field

      INTEGER, PARAMETER :: gt_ozone = 4
      ! Ozone field

      INTEGER, PARAMETER :: gt_scalar = 5
      ! Single point

      INTEGER, PARAMETER :: gt_compressed = 6
      ! Compressed points

      INTEGER, PARAMETER :: gt_LBC = 7
      ! Lateral Boundary Condition Field

! CYCLIC

      INTEGER, PARAMETER :: gt_nocyclic = 1
      ! No cyclic columns

      INTEGER, PARAMETER :: gt_optcyclic = 2
      ! Optional cyclic columns

      INTEGER, PARAMETER :: gt_cyclic = 3
      ! Includes cyclic columns

! Subroutine arguments:


      INTEGER                                                           &
     &  GRID_CODE                                                       &
                    ! IN : Grid type code from STASHmaster

     &, NORTH_LAT_IN                                                    &
                       ! IN : Latitude of Northern boundary
     &, SOUTH_LAT_IN                                                    &
                       ! IN : Latitude of Southern boundary
     &, WEST_LONG_IN                                                    &
                       ! IN : Longitude of Western boundary
     &, EAST_LONG_IN                                                    &
                       ! IN : Longitude of Eastern boundary

     &, START_ROW                                                       &
                    ! OUT : Row number of start of area
     &, END_ROW                                                         &
                    ! OUT : Row number of end of area
     &, START_COL                                                       &
                    ! OUT : Column number of start of area
     &, END_COL     ! OUT : Column number of end of area


! Local variables
      INTEGER                                                           &
     &  MODEL_TYPE                                                      &
                    ! model type of grid
     &, CONTENT                                                         &
                    ! content type of grid
     &, COVERAGE                                                        &
                    ! coverage of grid
     &, DOMAIN                                                          &
                    ! domain of grid
     &, CYCLIC                                                          &
                    ! does grid contain cyclic wrap columns
     &, fld_type    ! P, U or V field type

      INTEGER                                                           &
     &  NORTH_LAT                                                       &
                   ! Modifiable copies
     &, SOUTH_LAT                                                       &
                   ! of the input arguments
     &, WEST_LONG                                                       &
                   ! so that the format can
     &, EAST_LONG                                                       &
                   ! be changed if necessary

     &, NROWS                                                           &
                   ! number of rows
     &, NCOLS      ! number of columns

      LOGICAL                                                           &
     &  LAM_MODEL  ! True if this is a LAM configuration

      REAL                                                              &
     &  POLE_LAT                                                        &
                   ! Latitude of rotated pole (LAM only)
     &, POLE_LONG  ! Longitude of rotated pole (LAM only)

      INTEGER                                                           &
     &  ROW_ORDERING                                                    &
                      ! ordering North->South or South->North
     &, North_to_South                                                  &
                       ! indicates North->South ordering
     &, South_to_North ! indicates South->North ordering

      PARAMETER                                                         &
     &  (North_to_South=1,South_to_North=2)

      REAL                                                              &
     &  START_LAT                                                       &
                   ! Starting latitude
     &, START_LONG                                                      &
                   ! Starting longitude
     &, DEL_LAT                                                         &
                   ! Latitude grid spacing
     &, DEL_LONG                                                        &
                   ! Longitude grid spacing
     &, R_START_ROW                                                     &
                    ! first row number
     &, R_END_ROW                                                       &
                    ! last row number
     &, R_START_COL                                                     &
                    ! first column number
     &, R_END_COL   ! last column number

! External and subroutine calls

      EXTERNAL LLTOLL,LLTOEQ
      EXTERNAL GET_FLD_TYPE,EREPORT
      INTEGER GET_FLD_TYPE

!-----------------------------------------------------------------------

! Copy input arguments
      NORTH_LAT=NORTH_LAT_IN
      SOUTH_LAT=SOUTH_LAT_IN
      WEST_LONG=WEST_LONG_IN
      EAST_LONG=EAST_LONG_IN

      IF (WEST_LONG  <   0) WEST_LONG=WEST_LONG+360
      IF (EAST_LONG  <=  0) EAST_LONG=EAST_LONG+360

! Get information about grid type

! DEPENDS ON: gt_decode
      CALL GT_DECODE(GRID_CODE,                                         &
     &               MODEL_TYPE,CONTENT,COVERAGE,DOMAIN,CYCLIC)

! Get information of field type for PARVARS variables

! DEPENDS ON: get_fld_type
      fld_type=GET_FLD_TYPE(GRID_CODE)

! Find start latitude and longitude of full grid, grid spacing and
! size of grid

      IF (MODEL_TYPE  ==  gt_atmos) THEN

        START_LAT=H_A_FIRSTLAT
        START_LONG=H_A_FIRSTLONG
        DEL_LAT=H_A_NSSPACE
        DEL_LONG=H_A_EWSPACE
        NROWS=decomp_db_glsize(2,fld_type,decomp_standard_atmos)
        NCOLS=decomp_db_glsize(1,fld_type,decomp_standard_atmos)
        ROW_ORDERING=South_to_North

        IF (H_GLOBAL(A_IM) == 'N') THEN   ! LAM Configuration
          LAM_MODEL=.TRUE.
          POLE_LAT=H_A_POLELAT
          POLE_LONG=H_A_POLELONG
        ELSE
          LAM_MODEL=.FALSE.
        ENDIF

      ENDIF

! Ensure that DEL_LAT is always positive
! (ie. doesn't imply a direction)
      IF (DEL_LAT  <   0) DEL_LAT=-DEL_LAT

! This assumes the mass grid. The start latitude and longitude
! may be offset for velocity/U/V fields

      IF (CONTENT  ==  gt_velocity) THEN  ! B grid U/C
        START_LAT=START_LAT+DEL_LAT/2
        START_LONG=START_LONG+DEL_LONG/2
      ELSEIF (CONTENT  ==  gt_U_C) THEN   ! C grid U
        START_LONG=START_LONG+DEL_LONG/2
      ELSEIF( CONTENT  ==  gt_V_C) THEN   ! C grid V
        START_LAT=START_LAT-DEL_LAT/2
      ENDIF

! This assumes full domain. Now take account of zonal,meridional
! and scalar fields

      IF (DOMAIN  ==  gt_zonal) THEN
        NCOLS=1
      ELSEIF (DOMAIN  ==  gt_meridional) THEN
        NROWS=1
      ELSEIF (DOMAIN  ==  gt_scalar) THEN
        NCOLS=1
        NROWS=1
      ELSEIF (DOMAIN  ==  gt_ozone) THEN
        IF (ZonAvOzone) NCOLS=1
      ENDIF

! This is the global sizes - this may be all we need

      IF ((NORTH_LAT  ==  90) .AND. (SOUTH_LAT  ==  -90) .AND.          &
     &    (WEST_LONG  ==  0) .AND. (EAST_LONG  ==  360)) THEN

        START_ROW=1
        END_ROW=NROWS
        START_COL=1
        END_COL=NCOLS

      ELSE ! Not a simple case

! If this is a LAM configuration we need to transform the latitudes
! and longitudes relative to the rotated pole

        IF (LAM_MODEL) THEN
! DEPENDS ON: lltoll
          CALL LLTOLL(                                                  &
     &        NORTH_LAT,SOUTH_LAT,EAST_LONG,WEST_LONG,                  &
     &        POLE_LAT,POLE_LONG)
        ENDIF

! Make sure that DEL_LAT has a sign consistent with the
! grid ordering. Further up we have ensured that DEL_LAT
! is positive. Now we change it negative if necessary.
!               ( North->South => -ive
!                 South->North => +ive )

        IF (ROW_ORDERING  ==  North_to_South) DEL_LAT=-DEL_LAT

! Calculate the start and end rows

        IF (ROW_ORDERING  ==  North_to_South) THEN

          R_START_ROW=1.0   + (NORTH_LAT-START_LAT)/DEL_LAT
          R_END_ROW  =1.999 + (START_LAT-SOUTH_LAT)/DEL_LAT

        ELSE ! South->North

          R_START_ROW=1.0   + (SOUTH_LAT-START_LAT)/DEL_LAT
          R_END_ROW  =1.999 + (NORTH_LAT-START_LAT)/DEL_LAT

        ENDIF

        START_ROW=MAX(1.0,R_START_ROW)
        END_ROW=MIN(REAL(NROWS),R_END_ROW)

        IF (START_LONG  <   0) START_LONG=START_LONG+360.0

        IF ((REAL(START_LONG) + DEL_LONG*(NCOLS-1))  <=  360.0)         &
     &  THEN  ! If the total model domain doesn't cross
              ! the zero meridian

          R_START_COL=1.0 +   (WEST_LONG-START_LONG)/DEL_LONG
          R_END_COL=  1.999 + (EAST_LONG-START_LONG)/DEL_LONG

        ELSE ! model domain crosses the zero meridian

          IF (REAL(WEST_LONG)  <   START_LONG) THEN
            ! If the start is before the zero meridian

            R_START_COL=1.0  + (WEST_LONG + 360.0 - START_LONG)/        &
     &                         DEL_LONG

          ELSE ! the start lies after the meridian

            R_START_COL=1.0 + (REAL(WEST_LONG)-START_LONG)/             &
     &                        DEL_LONG

          ENDIF

          IF (REAL(EAST_LONG)  <   START_LONG) THEN
            ! If the end is after the zero meridian

            R_END_COL=1.0 + (EAST_LONG + 360.0 - START_LONG)/           &
     &                      DEL_LONG

          ELSE ! the end lies before the zero meridian

            R_END_COL=1.0 + (REAL(EAST_LONG) - START_LONG)/             &
     &                      DEL_LONG

          ENDIF

        ENDIF

        START_COL=MIN(REAL(NCOLS),MAX(1.0,R_START_COL))
        END_COL=MIN(REAL(NCOLS),MAX(1.0,R_END_COL))

      ENDIF

      RETURN
      END SUBROUTINE LLTORC

!- End of subroutine code ----------------------------------------


!+Use subroutine LLTOEQ to convert to equatorial lat/lon for LAM
! Subroutine Interface:


!- End of subroutine code ---------------------------------------------