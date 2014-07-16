

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE WRITHEAD---------------------------------------
!LL
!LL  Purpose: Writes out model dump header records on unit NFTOUT &
!LL           checks model and dump dimensions for consistency.
!LL           32-bit IEEE output option supported
!LL           64-bit IEEE output option supported
!LL
!LL  Written by A. Dickinson 31/01/90
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.5    28/03/95 mpp code: New code for parallel I/O
!LL                                              P.Burton
!     4.1    18/06/96 Changes to cope with changes in STASH addressing
!                     Author D.M. Goddard.
!    4.1  21/05/96  Correct conversion of LOOKUP(65-128) in obs
!                   files for IEEE/64 bits. D Robinson
!     4.4    17/07/97 Introduce conversion from ieee to Cray PVP
!                     numbers and reintroduce functionality for
!                     PVP machines
!                     Author: D.M. Goddard
!     4.5    25/08/98 Correct conversion of LOOKUP(65-128) for
!                     cx files and OPS obstores.
!                     Author D.M. Goddard
!     5.1    10/04/00 New reconfiguration and stdout control using
!                     PrintStatus. P.Selwood.
!     5.3    22/11/01 Enable mpp as the only option for
!                     small executables         E.Leung
!     6.0    11/09/03 Replaced ABORT with call to EREPORT.    P.Dando
!     6.0    20/08/03 Porting to NEC  E.Leung
!     6.2    28/10/05 Fixes for -64e option and remove obsolete code.
!                     P.Selwood
!     6.2    06/12/05 Removed the DIAG81 CPP DEF. T.Edwards
!     6.2    15/08/05 Minor reworkings for FCM. P.Selwood
!LL
!LL  Programming standard:
!LL           Unified Model Documentation Paper No 3
!LL           Version No 1 15/1/90
!LL
!LL  System component: W30
!LL
!LL  System task: F3
!LL
!LL  Documentation:
!LL           Unified Model Documentation Paper No F3
!LL           Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE WRITHEAD(NFTOUT,FIXHD,LEN_FIXHD,                       &
                                                        ! Intent (In)
     &                    INTHD,LEN_INTHD,                              &
     &                    REALHD,LEN_REALHD,                            &
     &                    LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,            &
     &                    ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,            &
     &                    COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,            &
     &                    FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,            &
     &                    EXTCNST,LEN_EXTCNST,                          &
     &                    DUMPHIST,LEN_DUMPHIST,                        &
     &                    CFI1,LEN_CFI1,                                &
     &                    CFI2,LEN_CFI2,                                &
     &                    CFI3,LEN_CFI3,                                &
     &                    LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,      &
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
     &                    START_BLOCK,ICODE,CMESSAGE)   ! Intent (Out)


      IMPLICIT NONE

      INTEGER                                                           &
     & NFTOUT                                                           &
                     !IN Unit no of dump
     &,LEN_FIXHD                                                        &
                     !IN Length of fixed length header
     &,LEN_INTHD                                                        &
                     !IN Length of integer header
     &,LEN_REALHD                                                       &
                     !IN Length of real header
     &,LEN1_LEVDEPC                                                     &
                     !IN 1st dim of level dep consts
     &,LEN2_LEVDEPC                                                     &
                     !IN 2ndt dim of level dep consts
     &,LEN1_ROWDEPC                                                     &
                     !IN 1st dim of row dep consts
     &,LEN2_ROWDEPC                                                     &
                     !IN 2nd dim of row dep consts
     &,LEN1_COLDEPC                                                     &
                     !IN 1st dim of column dep consts
     &,LEN2_COLDEPC                                                     &
                     !IN 2nd dim of column dep consts
     &,LEN1_FLDDEPC                                                     &
                     !IN 1st dim of field dep consts
     &,LEN2_FLDDEPC                                                     &
                     !IN 2nd dim of field dep consts
     &,LEN_EXTCNST                                                      &
                     !IN Length of extra constants
     &,LEN_DUMPHIST                                                     &
                     !IN Length of history block
     &,LEN_CFI1                                                         &
                     !IN Length of comp field index 1
     &,LEN_CFI2                                                         &
                     !IN Length of comp field index 2
     &,LEN_CFI3                                                         &
                     !IN Length of comp field index 3
     &,LEN1_LOOKUP                                                      &
                     !IN 1st dim of lookup
     &,LEN2_LOOKUP   !IN 2nd dim of lookup

      INTEGER                                                           &
     & LEN_DATA                                                         &
                      !IN Length of model data
     &,START_BLOCK    !OUT Pointer to position of each block.
                      !Should point to start of model data block on exit

      INTEGER ::  ICODE      !OUT Return code; successful=0
                             !                 error > 0







      CHARACTER*(80)                                                    &
     & CMESSAGE       !OUT Error message if ICODE > 0

      CHARACTER (Len=*), Parameter  ::  RoutineName = 'WRITHEAD'

      INTEGER                                                           &
     & FIXHD(LEN_FIXHD)                                                 &
                        !IN Fixed length header
     &,INTHD(LEN_INTHD)                                                 &
                        !IN Integer header
     &,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                                  &
                                       !IN PP lookup tables
     &,CFI1(LEN_CFI1+1)                                                 &
                        !IN Compressed field index no 1
     &,CFI2(LEN_CFI2+1)                                                 &
                        !IN Compressed field index no 2
     &,CFI3(LEN_CFI3+1) !IN Compressed field index no 3

      REAL                                                              &
     & REALHD(LEN_REALHD)                                               &
                          !IN Real header
     &,LEVDEPC(1+LEN1_LEVDEPC*LEN2_LEVDEPC)                             &
                                            !IN Lev dep consts
     &,ROWDEPC(1+LEN1_ROWDEPC*LEN2_ROWDEPC)                             &
                                            !IN Row dep consts
     &,COLDEPC(1+LEN1_COLDEPC*LEN2_COLDEPC)                             &
                                            !IN Col dep consts
     &,FLDDEPC(1+LEN1_FLDDEPC*LEN2_FLDDEPC)                             &
                                            !IN Field dep consts
     &,EXTCNST(LEN_EXTCNST+1)                                           &
                                !IN Extra constants
     &,DUMPHIST(LEN_DUMPHIST+1) !IN History block


! Local arrays:------------------------------------------------
! None
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
! -------------------------------------------------------------
! External subroutines called:---------------------------------
      EXTERNAL IOERROR,POSERROR,PR_FIXHD,BUFFOUT
!*-------------------------------------------------------------
! Local variables:---------------------------------------------
      INTEGER LEN_IO
      REAL A
! -------------------------------------------------------------
! Comdecks:----------------------------------------------------------
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
! -------------------------------------------------------------
      ICODE=0
      CMESSAGE=' '

!L 1. Buffer out fixed length header record

! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,FIXHD(1),LEN_FIXHD,LEN_IO,A)


! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= LEN_FIXHD)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of fixed length header',A,LEN_IO       &
     &               ,LEN_FIXHD)
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=1
        RETURN
      ENDIF

      START_BLOCK=LEN_FIXHD+1

! Check validity of data and print out fixed header information
      IF(PrintStatus >= PrStatus_Oper) THEN

! DEPENDS ON: pr_fixhd
      CALL PR_FIXHD(FIXHD,LEN_FIXHD,LEN_INTHD,LEN_REALHD,LEN1_LEVDEPC   &
     &,LEN2_LEVDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC &
     &,LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,LEN_DUMPHIST,LEN_CFI1      &
     &,LEN_CFI2,LEN_CFI3,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA               &
     &,ICODE,CMESSAGE)

       IF(ICODE >  0) THEN
        CMESSAGE="WRITHEAD: Error returned from PR_FIXHD"
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       END IF
      END IF

!L 2. Buffer out integer constants

      IF(FIXHD(100) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(100) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('integer constants',START_BLOCK,100,FIXHD(100))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=2
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,INTHD(1),FIXHD(101),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(101))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of integer constants',A,LEN_IO,        &
     &               FIXHD(101))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=3
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(101)

      ENDIF

!L 3. Buffer out real constants

      IF(FIXHD(105) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(105) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('real constants',START_BLOCK,105,FIXHD(105))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=4
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

! Check for I/O errors
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,REALHD(1),FIXHD(106),LEN_IO,A)

       IF(A /= -1.0.OR.LEN_IO /= FIXHD(106))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of real constants',A,LEN_IO,           &
     &                FIXHD(106))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=5
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(106)

      ENDIF

!L 4. Buffer out level dependent constants

      IF(FIXHD(110) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(110) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('level dependent constants',                      &
     &  START_BLOCK,110,FIXHD(110))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=6
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,LEVDEPC(1),FIXHD(111)*FIXHD(112),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(111)*FIXHD(112))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of level dependent constants',A,LEN_IO,&
     &               FIXHD(111)*FIXHD(112))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=7
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(111)*FIXHD(112)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' LEVEL DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(111)*FIXHD(112)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 5. Buffer out row dependent constants

      IF(FIXHD(115) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(115) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('row dependent constants',                        &
     &  START_BLOCK,115,FIXHD(115))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=8
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,ROWDEPC(1),FIXHD(116)*FIXHD(117),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(116)*FIXHD(117))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of row dependent constants',A,LEN_IO,  &
     &                FIXHD(116)*FIXHD(117))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=9
        RETURN
      ENDIF


       START_BLOCK=START_BLOCK+FIXHD(116)*FIXHD(117)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' ROW DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(116)*FIXHD(117)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 6. Buffer out column dependent constants

      IF(FIXHD(120) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(120) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('column dependent constants',                     &
     &  START_BLOCK,120,FIXHD(120))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=10
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,COLDEPC(1),FIXHD(121)*FIXHD(122),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(121)*FIXHD(122))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of column dependent constants',A,LEN_IO&
     &              ,FIXHD(121)*FIXHD(122))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=11
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(121)*FIXHD(122)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COLUMN DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(121)*FIXHD(122)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 7. Buffer out constants stored as fields

      IF(FIXHD(125) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(125) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('fields of constants',                            &
     &  START_BLOCK,125,FIXHD(125))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=12
        RETURN
       ENDIF

! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,FLDDEPC(1),FIXHD(126)*FIXHD(127),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(126)*FIXHD(127))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of field dependent constants',A,LEN_IO,&
     &               FIXHD(126)*FIXHD(127))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=13
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(126)*FIXHD(127)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' FIELD DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(126)*FIXHD(127)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 8. Buffer out extra constants

      IF(FIXHD(130) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(130) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('extra constants',                                &
     &  START_BLOCK,130,FIXHD(130))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=14
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,EXTCNST(1),FIXHD(131),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(131))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out extra constants',A,LEN_IO,             &
     &               FIXHD(131))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=15
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(131)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' EXTRA CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(131)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 9. Buffer out temporary history block

      IF(FIXHD(135) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(135) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('history',                                        &
     &  START_BLOCK,136,FIXHD(136))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=16
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,DUMPHIST(1),FIXHD(136),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(136))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of history file',A,LEN_IO,             &
     &               FIXHD(136))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=17
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(136)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' TEMPORARY HISTORY BLOCK'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(136)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 10. Buffer out compressed field index1

      IF(FIXHD(140) >  0)THEN

! Check for error in file pointers

       IF(FIXHD(140) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('compressed field index1',                        &
     &  START_BLOCK,140,FIXHD(140))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=18
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF


! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,CFI1(1),FIXHD(141),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(141))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of compressed index1',A,LEN_IO,        &
     &               FIXHD(141))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=19
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(141)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 1'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(141)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 11. Buffer out compressed field index2

      IF(FIXHD(142) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(142) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('compressed field index2',                        &
     &  START_BLOCK,142,FIXHD(142))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=20
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,CFI2(1),FIXHD(143),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(143))THEN
! DEPENDS ON: ioerror
       CALL IOERROR('buffer out of compressed index2',A,LEN_IO,         &
     &               FIXHD(143))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=21
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(143)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 2'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(143)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 12. Buffer out compressed field index3

      IF(FIXHD(144) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(144) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('compressed field index3',                        &
     &  START_BLOCK,144,FIXHD(144))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=22
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,CFI3(1),FIXHD(145),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(145))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of compressed index3',A,LEN_IO,        &
     &               FIXHD(145))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=23
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(145)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 3'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(145)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 13. Buffer out lookup table

      IF(FIXHD(150) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(150) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('lookup table',                                   &
     &  START_BLOCK,150,FIXHD(150))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=24
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,LOOKUP(1,1),FIXHD(151)*FIXHD(152),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(151)*FIXHD(152))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of lookup table',A,LEN_IO,             &
     &               FIXHD(151)*FIXHD(152))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=25
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(151)*FIXHD(152)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' LOOKUP TABLE'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(151)*FIXHD(152)
      ENDIF ! if mype  ==  0
      End If

! No consistency checks for parallel code. The LOOKUP headers don't
! match the data layout in memory within a PE.

      ENDIF

      RETURN
      END SUBROUTINE WRITHEAD
