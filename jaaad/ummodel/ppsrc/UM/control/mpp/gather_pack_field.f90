

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Gathers a field from many processors to one processor and packs
!
! Subroutine Interface:
        SUBROUTINE GATHER_PACK_FIELD(                                   &
                                LOCAL_FIELD,GLOBAL_FIELD,               &
                                LOCAL_ROW_LEN,LOCAL_ROWS,               &
                                GLOBAL_ROW_LEN,GLOBAL_ROWS,             &
                                GRID_TYPE,HALO_TYPE,                    &
                                GLOBAL_GATHER_PE,PROC_GROUP,            &
                                PACKING, IM_IDENT, LRLE, PACKING_TYPE,  &
                                NUM_OUT,                                &
                                COMP_ACCRCY, RMDI)

      IMPLICIT NONE

!
! Description:
!  Takes a model field that has been decomposed over a group of
!  processors, and gathers the data together so that one processor
!  contains the entire global field. Optionally WGDOS packs the
!  gathered data as it is sent to get parallelism of the packing
!  process.
!
! Method:
!  In the most simple situation,
!  a send and receive map is constructed which instructs the GCOM
!  permute operation to do a gather from all processors in the
!  group to the GATHER_PE.
!
!  In the more ususal version, there is a 2 stage gather. The
!  first gathers a "block_factor" number of rows together. These
!  are packed (optionally) and then a final gather gets the full
!  field into place. This will require some adjustment to ensure
!  the packed data is correct.
!
!
! Current Code Owner: Paul Selwood
!
! Subroutine Arguments:

        INTEGER, INTENT(IN) ::                                          &
        LOCAL_ROW_LEN                                                   &
                         ! IN length of rows in local part of field
      , LOCAL_ROWS                                                      &
                           ! IN number of rows in local part of field
      , GLOBAL_ROW_LEN                                                  &
                         ! IN length of rows in global field
      , GLOBAL_ROWS                                                     &
                           ! IN number of rows in global field
      , GRID_TYPE                                                       &
                           ! IN type (P,U or V) of grid
      , HALO_TYPE                                                       &
                           ! IN halo type (hence width) of grid
      , GLOBAL_GATHER_PE                                                &
                         ! IN processor to gather global field to
      , PROC_GROUP         ! IN group ID of processors involved here

!
! Optional Arguments to handle the COEX packing if necessary
!
      LOGICAL, INTENT(IN) ::                                            &
        PACKING                                                         &
                           ! IN: Set .true. if packing of the input
                           !     field is to be packed
      , LRLE               ! IN: True if Run Length Encoding is required

      INTEGER, INTENT(IN) ::                                            &
        IM_IDENT           ! IN: Internal model identifier

      INTEGER, INTENT(INOUT) ::                                         &
        PACKING_TYPE       ! IN/OUT: This flag is zero on input,
                           !         then stash packing is selected,
                           !         and the routine returns the
                           !         packing flag.
                           !
                           !         If the variable is set to 2 on
                           !         input then 32-bit packing for
                           !         dumpfiles is selected

      INTEGER, INTENT(OUT) ::                                           &
        NUM_OUT            ! OUT: Number of 32-bit IBM words in the
                           !      Packed field for WDGOS packing

      INTEGER, INTENT(IN) ::                                            &
        COMP_ACCRCY        ! IN: Packing Accuracy in Power of 2

      REAL, INTENT(IN) ::                                               &
        RMDI               ! IN: Missing data indicator
!
! Remaining Non-Optional Arguments
!
      REAL, INTENT(IN) ::                                               &
        LOCAL_FIELD(LOCAL_ROW_LEN*LOCAL_ROWS)
                           ! IN local part of field

      REAL, INTENT(OUT) ::                                              &
        GLOBAL_FIELD(GLOBAL_ROW_LEN*GLOBAL_ROWS)
                           ! OUT (on PE GATHER_PE) global field


! Parameters and Common blocks

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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! GC - General Communication primitives package. For use on
! multiprocessor shared memory and message passing systems.
!
!
! LICENSING TERMS
!
!  GC is provided free of charge. Unless otherwise agreed with SINTEF,
!  use and redistribution in source and binary forms are permitted
!  provided that
!
!      (1) source distributions retain all comments appearing within
!          this file header, and
!
!      (2) distributions including binaries display the following
!          acknowledgement:
!
!              "This product includes software developed by SINTEF.",
!
!          in the documentation or other materials provided with the
!          distribution and in all advertising materials mentioning
!          features or use of this software.
!
!  The name of SINTEF may not be used to endorse or promote products
!  derived from this software without specific prior written
!  permission.  SINTEF disclaims any warranty that this software will
!  be fit for any specific purposes. In no event shall SINTEF be liable
!  for any loss of performance or for indirect or consequential damage
!  or direct or indirect injury of any kind. In no case shall SINTEF
!  be liable for any representation or warranty make to any third party
!  by the users of this software.
!
!
! Fortran header file. PLEASE use the parameter variables in user
! routines calling GC and NOT the numeric values. The latter are
! subject to change without further notice.
!
!---------------------------------------------- ------------------------
! $Id: gps0h501,v 1.6 2000/04/17 10:05:47 t11ps Exp $
! (C) Jorn Amundsen, Roar Skaalin, SINTEF Industrial Mathematics.

!    4.4   30/09/97  Added code to permit the SHMEM/NAM timeout
!                    value to be set from a shell variable.
!                      Author: Bob Carruthers  Cray Research.
!    5.1   17/04/00  Fixed/Free format. P.Selwood.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!     GC general options
      INTEGER, PARAMETER :: GC_OK         =     0
      INTEGER, PARAMETER :: GC_FAIL       =    -1
      INTEGER, PARAMETER :: GC_NONE       =     0
      INTEGER, PARAMETER :: GC_ANY        =    -1
      INTEGER, PARAMETER :: GC_DONTCARE   =    -1
      INTEGER, PARAMETER :: GC_SHM_DIR    =     1
      INTEGER, PARAMETER :: GC_SHM_SAFE   =     2
      INTEGER, PARAMETER :: GC_NAM_TIMEOUT=     4
      INTEGER, PARAMETER :: GC_SHM_GET    = -9999
      INTEGER, PARAMETER :: GC_SHM_PUT    = -9998
      INTEGER, PARAMETER :: GC_USE_GET    = -9999
      INTEGER, PARAMETER :: GC_USE_PUT    = -9998

!     GC functions
      INTEGER GC_COMLEN, GC_ISIZE, GC_RSIZE, GC_ME, GC_NPROC

!     GC groups (GCG) support
      INTEGER, PARAMETER :: GC_ALLGROUP = 0
      INTEGER, PARAMETER :: GCG_ALL = GC_ALLGROUP

!     GC groups (GCG) functions
      INTEGER GCG_ME

!     GC reserved message tags
      INTEGER, PARAMETER :: GC_MTAG_LOW   = 999999901
      INTEGER, PARAMETER :: GC_MTAG_HIGH  = 999999999

!     GCG_RALLETOALLE index parameters
      INTEGER, PARAMETER :: S_DESTINATION_PE = 1
      INTEGER, PARAMETER :: S_BASE_ADDRESS_IN_SEND_ARRAY = 2
      INTEGER, PARAMETER :: S_NUMBER_OF_ELEMENTS_IN_ITEM = 3
      INTEGER, PARAMETER :: S_STRIDE_IN_SEND_ARRAY = 4
      INTEGER, PARAMETER :: S_ELEMENT_LENGTH = 5
      INTEGER, PARAMETER :: S_BASE_ADDRESS_IN_RECV_ARRAY = 6
      INTEGER, PARAMETER :: S_STRIDE_IN_RECV_ARRAY = 7

      INTEGER, PARAMETER :: R_SOURCE_PE = 1
      INTEGER, PARAMETER :: R_BASE_ADDRESS_IN_RECV_ARRAY = 2
      INTEGER, PARAMETER :: R_NUMBER_OF_ELEMENTS_IN_ITEM = 3
      INTEGER, PARAMETER :: R_STRIDE_IN_RECV_ARRAY = 4
      INTEGER, PARAMETER :: R_ELEMENT_LENGTH = 5
      INTEGER, PARAMETER :: R_BASE_ADDRESS_IN_SEND_ARRAY = 6
      INTEGER, PARAMETER :: R_STRIDE_IN_SEND_ARRAY = 7

! Local variables

      INTEGER                                                           &
         send_map(7,1)                                                  &
      ,  receive_map(7,MAXPROC)                                         &
      ,  n_mess_to_recv                                                 &
      ,  n_mess_to_send                                                 &
      ,  send_map_2(7,1)                                                &
      ,  receive_map_2(7,MAXPROC)                                       &
      ,  n_mess_to_recv_2                                               &
      ,  n_mess_to_send_2

      INTEGER                                                           &
        old_GLOBAL_ROW_LEN                                              &
                              ! value on last call
      , old_GLOBAL_ROWS                                                 &
                              ! value on last call
      , old_PROC_GROUP                                                  &
                              ! value on last call
      , old_GATHER_PE                                                   &
                              ! value on last call
      , old_DECOMP                                                      &
                              ! value on last call
      , old_GRID_TYPE                                                   &
                              ! value on last call
      , old_HALO_TYPE         ! value on last call

        SAVE send_map,n_mess_to_send,receive_map,n_mess_to_recv,        &
           old_GLOBAL_ROW_LEN,old_GLOBAL_ROWS,old_PROC_GROUP,           &
           old_GATHER_PE,old_DECOMP,                                    &
           old_GRID_TYPE,old_HALO_TYPE
      DATA old_GLOBAL_ROW_LEN,old_GLOBAL_ROWS,old_PROC_GROUP,           &
           old_GATHER_PE,old_DECOMP,                                    &
           old_GRID_TYPE,old_HALO_TYPE                                  &
         / -1234, -1234, -1234, -1234, -1234, -1234, -1234/

      INTEGER                                                           &
        iproc                                                           &
      , jproc                                                           &
      , kproc                                                           &
      , lproc                                                           &
      , info                                                            &
      , flag                                                            &
      , gather_pe                                                       &
                             ! Local gather PE for a group/block of
                             ! rows
      , row_start_pe                                                    &
                             ! First PE in a block - may or may not
                             ! the gather PE
      , data_address                                                    &
      , data_size                                                       &
      , block_pe                                                        &
                             ! The PE holding the current block of
                             ! rows on the GLOBAL_GATHER_PE
      , pes_per_block                                                   & 
                             ! Number of PE's in a block (nproc_x*
                             ! block_factor)
      , n_local_rows                                                    &
                             ! Number of rows per block
      , packed_buffer_size                                              &
                             ! Size of the buffer to hold packed data
      , length_fullwrd                                                  &
                             ! Size of full word on this machine 64-bit
                             ! for Cray
      , local_packing_type   ! Local copy of the packing_type, which
                             ! maybe not present

      INTEGER :: icode       ! error code
 
      INTEGER, ALLOCATABLE :: icomp(:)
                             ! Local array to hold compressed data

! The block_factor is the number of rows that are gathered together
! for packing purposes before the final gather. This is best set to
! 1 where maximum parallelism can be exploited (eg on the NEC SX-6).
! Massively parallel systems may benefit from a higher value.
      INTEGER, PARAMETER   :: block_factor = 1
 
      REAL  :: my_global_buffer(global_row_len*global_rows+2)
 
      CHARACTER (LEN=*), PARAMETER :: ROUTINENAME='GATHER_PACK_FIELD'
      CHARACTER (LEN=80)           :: CMESSAGE

!-------------------------------------------------------

! Compute the local gather PE - head of row PE
      pes_per_block=block_factor*nproc_x
      gather_pe=(mype/pes_per_block)*pes_per_block
      row_start_pe=gather_pe
 
! Check if the global_gather_pe is in the block - if so
! make sure it is the gather_pe
      IF(global_gather_pe >= gather_pe .AND.                            &
         global_gather_pe <= MIN(nproc, gather_pe+pes_per_block)-1) THEN
        gather_pe=global_gather_pe
      END IF
 
! 0.0 Can we use the same send/receive map that we calculated
!     last time round?

      IF ((GLOBAL_ROW_LEN  /=  old_GLOBAL_ROW_LEN) .OR.                 &
          (GLOBAL_ROWS     /=  old_GLOBAL_ROWS   ) .OR.                 &
          (PROC_GROUP      /=  old_PROC_GROUP    ) .OR.                 &
          (GATHER_PE       /=  old_GATHER_PE     ) .OR.                 &
          (GRID_TYPE       /=  old_GRID_TYPE     ) .OR.                 &
          (HALO_TYPE       /=  old_HALO_TYPE     ) .OR.                 &
          (current_decomp_type  /=  old_DECOMP  )) THEN

!       Different arguments from the last call so we need
!       to calculate a new send/receive map
        IF (GRID_TYPE  ==  fld_type_unknown) THEN
          WRITE(6,*) 'GATHER_PACK_FIELD : Bad field type'
          WRITE(6,*) 'Field will not be scattered.'
          CMESSAGE='GATHER_PACK_FIELD : Bad field type'
          ICODE=1

! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME,ICODE,CMESSAGE)
        END IF


! 2.0 Set up send map

        send_map(S_DESTINATION_PE,1) = GATHER_PE

        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,1) =                      &
          halosize(2,halo_type)*LOCAL_ROW_LEN+                          &
          1+halosize(1,halo_type)

        send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)=blsize(2,grid_type)

        send_map(S_STRIDE_IN_SEND_ARRAY,1) = LOCAL_ROW_LEN

        send_map(S_ELEMENT_LENGTH,1) =                                  &
          LOCAL_ROW_LEN-2*halosize(1,halo_type)

        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) =                      &
          datastart_f(1,grid_type) +                                    &
         (datastart_f(2,grid_type)-1)*GLOBAL_ROW_LEN

        send_map(S_STRIDE_IN_RECV_ARRAY,1) = GLOBAL_ROW_LEN

        n_mess_to_send=1

! 3.0 Set up the receive map (for PE GATHER_PE only)

        n_mess_to_recv=0

        IF (mype  ==  GATHER_PE) THEN
 
! Loop over PE's in this block
          DO jproc=row_start_pe,                                        &
           min(nproc, row_start_pe+pes_per_block)-1
            iproc=jproc-row_start_pe

            receive_map(R_SOURCE_PE,iproc+1) = jproc

            receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) =         &
                g_datastart_f(1,grid_type,jproc)+                       &
                (g_datastart_f(2,grid_type,jproc)-1)                    &
                 *GLOBAL_ROW_LEN

            receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) =         &
                g_blsize(2,grid_type,jproc)

            receive_map(R_STRIDE_IN_RECV_ARRAY,iproc+1) =               &
              GLOBAL_ROW_LEN

            receive_map(R_ELEMENT_LENGTH,iproc+1) =                     &
                g_blsize(1,grid_type,jproc)

            receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) =         &
              halosize(2,halo_type)*                                    &
                g_lasize(1,grid_type,halo_type,jproc)+                  &
              halosize(1,halo_type)+1

            receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) =               &
                g_lasize(1,grid_type,halo_type,jproc)

            n_mess_to_recv=n_mess_to_recv+1

          END DO
        END IF

        old_GLOBAL_ROW_LEN=GLOBAL_ROW_LEN
        old_GLOBAL_ROWS=GLOBAL_ROWS
        old_PROC_GROUP=PROC_GROUP
        old_GATHER_PE=GATHER_PE
        old_DECOMP=current_decomp_type
        old_GRID_TYPE=GRID_TYPE
        old_HALO_TYPE=HALO_TYPE

      END IF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

      flag=0  ! This is currently ignored at GCG v1.1
      info=GC_NONE

      local_packing_type = 0
 
! Only the gather PE's need to do anything now
!
!      Run length encoding applies only to unpacked ocean fields, but
!      most ocean fields remain unpacked even when packing profiles
!      5 or 6 are set. Hence selecting, for example, both packing
!      profile 5 and run length encoding makes sense.
 
      IF(PACKING .AND. GLOBAL_ROWS  >=  2) THEN
        ! 1. Climate wgdos packing has been selected for current
        !    file stream via UMUI
        IF(COMP_ACCRCY  <=  -99 .AND. LRLE .AND.                        &
          im_ident  ==  ocean_im )THEN
           ! 2. STASH packing profile for the field is -99
           ! 3. Run Length Encoding has been selected
           ! 4. Submodel is Ocean
           PACKING_TYPE = 4
        ELSE IF(COMP_ACCRCY  >   -99) THEN
           ! 2. STASH packing profile for the field is set.
           PACKING_TYPE = 1
        END IF
      ELSE
         ! 1. Packing may or may not have been selected for current
         !    file stream. This section of code is not executed when
         !    GRIB packing selected.
         IF (LRLE .AND. im_ident  ==  ocean_im )THEN
           ! 2. Run Length Encoding has been selected
           ! 3. Submodel is Ocean
           PACKING_TYPE = 4
         END IF
      END IF

      num_out=0
      local_packing_type = packing_type

! If there is no packing, then use the output buffer
! if I am the gather PE, or if it has been reduced to a single
! stage gather.  Otherwise, use a temporary buffer.
      IF(local_packing_type /= 1 .AND.                                  &
       ((gather_pe == global_gather_pe) .OR.                            &
        (block_factor >= nproc_y))) THEN

        CALL GCG_RALLTOALLE(LOCAL_FIELD,send_map,n_mess_to_send,        &
                            LOCAL_ROW_LEN*LOCAL_ROWS,                   &
                            GLOBAL_FIELD,receive_map,n_mess_to_recv,    &
                            GLOBAL_ROW_LEN*GLOBAL_ROWS,                 &
                            PROC_GROUP,flag,icode)
      ELSE

        CALL GCG_RALLTOALLE(LOCAL_FIELD,send_map,                       &
                            n_mess_to_send,                             &
                            LOCAL_ROW_LEN*LOCAL_ROWS,                   &
                            my_global_buffer,receive_map,               &
                            n_mess_to_recv,                             &
                            GLOBAL_ROW_LEN*GLOBAL_ROWS,                 &
                            PROC_GROUP,flag,icode)
      END IF
 
      IF (ICODE /= 0) THEN
        CMESSAGE = ' Error in GCG_RALLTOALLE'

! DEPENDS ON: ereport
        CALL EREPORT(ROUTINENAME,ICODE,CMESSAGE)
      END IF

 
! Now we must check if the packing_type is 1
      IF(local_packing_type == 1) THEN
 
! Only gather_pe's need to do anything
        IF(mype == gather_pe) THEN
 
! Work out how much local data we have
          n_local_rows=0
          DO kproc=1, block_factor
            IF(row_start_pe+(kproc-1)*nproc_x <  nproc) THEN
              n_local_rows=n_local_rows+                                &
               g_blsize(2,grid_type,row_start_pe+(kproc-1)*nproc_x)
            END IF
          END DO
 
! Setup a buffer for the packed data
          packed_buffer_size=n_local_rows*global_row_len+2
          ALLOCATE (icomp(packed_buffer_size))
 
! Pack the data
          icode=0
          length_fullwrd=64
 
! Check if the2 stage gather has been eliminated - if so
! put the data straight into the output array
          IF(block_factor >= nproc_y) THEN
! DEPENDS ON: coex
            CALL coex(                                                  &
               my_global_buffer(                                        &
                         g_datastart_f(1,grid_type,row_start_pe)+       &
                         (g_datastart_f(2,grid_type,row_start_pe)-1)*   &
                                global_row_len),                        &
                                global_row_len*global_rows,             &
                                global_field, packed_buffer_size,       &
                                global_row_len, n_local_rows,           &
                                num_out,                                &
                                comp_accrcy, .true., rmdi,              &
                                length_fullwrd,                         &
                                icode, cmessage)
 
! Still doing 2 stage gather
          ELSE
 
! DEPENDS ON: coex
            CALL coex(                                                  &
               my_global_buffer(g_datastart_f(1,grid_type,row_start_pe)+&
                         (g_datastart_f(2,grid_type,row_start_pe)-1)*   &
                                global_row_len),                        &
                                global_row_len*global_rows,             &
                                icomp, packed_buffer_size,              &
                                global_row_len, n_local_rows,           &
                                num_out,                                &
                                comp_accrcy, .true., rmdi,              &
                                length_fullwrd,                         &
                                icode, cmessage)
          END IF

          IF (ICODE /= 0) THEN
            cmessage='Error in COEX'

! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ICODE, CMESSAGE)
          END IF


        END IF ! mype == gather_pe
 
      END IF !packing_type is 1
 
      IF(local_packing_type == 1) THEN

        IF(mype == gather_pe) THEN

          IF(block_factor <  nproc_y) THEN

            CALL gc_rsend(1001+mype, (num_out+1)/2, global_gather_pe,   &
             info,                                                      &
             my_global_buffer(                                          &
             g_datastart_f(1,grid_type,row_start_pe)+                   &
            (g_datastart_f(2,grid_type,row_start_pe)-1)*global_row_len),&
             icomp)

          END IF ! are we doing 2 stage gather?

        END IF ! mype == gather_pe

      END IF !  packing_type is 1
 
      IF(local_packing_type == 1) THEN

! Now pack the into the output buffer on the global_gather_pe
        IF(mype == global_gather_pe) THEN

! Check if we are doing 2 stage gather
          IF(block_factor <  nproc_y) THEN
 
! Loop over each processors Contribution
            DO jproc=0, nproc_y-1, block_factor
              block_pe=jproc*nproc_x
 
! Preserve the row_start_pe for this block of data
              iproc=block_pe

! If the global_gather_pe is in the current block, then the block_pe
! needs to be set to global_gather_pe, not the row/block leader
!
              IF(row_start_pe == block_pe) THEN
                block_pe=global_gather_pe
              END IF
 
              n_local_rows=0
              DO kproc=1, block_factor
                IF(iproc+(kproc-1)*nproc_x <  nproc) THEN
                  n_local_rows=n_local_rows+                            &
                   g_blsize(2,grid_type,iproc+(kproc-1)*nproc_x)
                END IF
              END DO
 
              CALL gc_rrecv(1001+block_pe,                              &
               n_local_rows*global_row_len+2,                           &
               block_pe, info,                                          &
               my_global_buffer(                                        &
               g_datastart_f(1,grid_type,iproc)+                        &
               (g_datastart_f(2,grid_type,iproc)-1)*global_row_len),    &
               icomp)

! DEPENDS ON: unite_coex_files
              CALL unite_coex_files(                                    &
                my_global_buffer(g_datastart_f(1,grid_type,iproc)+      &
               (g_datastart_f(2,grid_type,iproc)-1)*global_row_len),    &
                global_field, num_out, iproc)
            END DO

          END IF ! are we doing 2 stage gather?

        END IF ! mype == global_gather_pe

      ELSE

! Normal Gather Operation, without Packing
! 5.0 Set up the second send map

        IF(block_factor <  nproc_y) THEN

          IF(mype == gather_pe .AND. mype /= global_gather_pe) THEN

            send_map_2(S_DESTINATION_PE,1) = global_gather_pe

            send_map_2(S_BASE_ADDRESS_IN_SEND_ARRAY,1) =                &
                  g_datastart_f(1,grid_type,mype)+                      &
                  (g_datastart_f(2,grid_type,mype)-1)*GLOBAL_ROW_LEN

            send_map_2(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)=0
            DO kproc=1, block_factor
              IF(row_start_pe+(kproc-1)*nproc_x <  nproc) THEN
                send_map_2(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)=             &
                  send_map_2(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)+           &
                  g_blsize(2,grid_type,mype+(kproc-1)*nproc_x)
              END IF
            END DO

            send_map_2(S_STRIDE_IN_SEND_ARRAY,1) =                      &
                  GLOBAL_ROW_LEN

            send_map_2(S_ELEMENT_LENGTH,1) =                            &
                  GLOBAL_ROW_LEN

            send_map_2(S_BASE_ADDRESS_IN_RECV_ARRAY,1) =                &
                  g_datastart_f(1,grid_type,mype)+                      &
                  (g_datastart_f(2,grid_type,mype)-1)*GLOBAL_ROW_LEN

            send_map_2(S_STRIDE_IN_RECV_ARRAY,1) =                      &
                  GLOBAL_ROW_LEN

            n_mess_to_send_2=1

          ELSE

            n_mess_to_send_2=0

          END IF

! 6.0 Set up the second receive map (for PE GLOBAL_GATHER_PE only)

          n_mess_to_recv_2=0

          IF(mype == global_gather_pe) THEN

            iproc=0
            DO jproc=0, nproc_y-1, block_factor
 
! Compute the block PE for this group of rows, and check it is not
! not the global_gather_pe
              block_pe=jproc*nproc_x
              lproc=block_pe
 
! Check if the global_gather_pe is in the block - if so
! make sure it is the block_pe
              IF(global_gather_pe >= block_pe .AND.                     &
                global_gather_pe <=                                     &
                MIN(nproc, block_pe+pes_per_block)-1) THEN
                 block_pe=global_gather_pe
              END IF
 
              IF (block_pe  /=  global_gather_pe) THEN

                receive_map_2(R_SOURCE_PE,iproc+1) = block_pe

                receive_map_2(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) =   &
                  g_datastart_f(1,grid_type,block_pe)+                  &
                  (g_datastart_f(2,grid_type,block_pe)-1)*GLOBAL_ROW_LEN

                receive_map_2(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = 0
                DO kproc=1, block_factor
                  IF(lproc+(kproc-1)*nproc_x <  nproc) THEN
                    receive_map_2(R_NUMBER_OF_ELEMENTS_IN_ITEM,         &
                                  iproc+1) =                            &
                     receive_map_2(R_NUMBER_OF_ELEMENTS_IN_ITEM,        &
                                  iproc+1) +                            &
                     g_blsize(2,grid_type,lproc+(kproc-1)*              &
                              nproc_x)
                  END IF
                END DO

                receive_map_2(R_STRIDE_IN_RECV_ARRAY,iproc+1) =         &
                  GLOBAL_ROW_LEN

                receive_map_2(R_ELEMENT_LENGTH,iproc+1) =               &
                  GLOBAL_ROW_LEN

                receive_map_2(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) =   &
                  g_datastart_f(1,grid_type,block_pe)+                  &
                  (g_datastart_f(2,grid_type,block_pe)-1)*GLOBAL_ROW_LEN

                receive_map_2(R_STRIDE_IN_SEND_ARRAY,iproc+1) =         &
                  GLOBAL_ROW_LEN

                iproc=iproc+1
                n_mess_to_recv_2=n_mess_to_recv_2+1

              END IF ! block_pe  /=  global_gather_pe

            END DO ! gather_pe's

          END IF ! mype == global_gather_pe

! 7.0 Do the exchange of data

          flag=0  ! This is currently ignored at GCG v1.1
          info=GC_NONE

          CALL GCG_RALLTOALLE(my_global_buffer,send_map_2,              &
                              n_mess_to_send_2,                         &
                              GLOBAL_ROW_LEN*GLOBAL_ROWS,               &
                              GLOBAL_FIELD,receive_map_2,               &
                              n_mess_to_recv_2,                         &
                              GLOBAL_ROW_LEN*GLOBAL_ROWS,               &
                              PROC_GROUP,flag,icode)

        END IF ! packing_type is 1

      END IF ! are we doing 2 stage gather
 
! Deallocate the temporary buffer for coex processed data
      IF(ALLOCATED(icomp)) DEALLOCATE(icomp)
 

      RETURN
      END SUBROUTINE GATHER_PACK_FIELD
