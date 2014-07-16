! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2000, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************

!+ Parallel UM: Perform 2D data decomposition

! Subroutine interface:
      SUBROUTINE DECOMPOSE_SMEXE(global_row_len, global_n_rows,         &
     &                           extended_halo_EW,extended_halo_NS,     &
     &                           TOT_LEVELS)
      IMPLICIT NONE

! Description:
! This routine performs a 2D decomposition - taking the global X
! (global_row_len) and Y (global_n_rows) data sizes and decomposing
! across 1 processor in both the X & Y direction for small executables.
! Local data size should be the same as Global field, as no halo
! is included.

! Method:
! Various data is stored in the COMMON block DECOMPDB for later use.


! Subroutine Arguments:

      INTEGER                                                           &

     &  global_row_len,                                                 &
                          ! IN  :number of E-W points of entire model
     &  global_n_rows,                                                  &
                          ! IN  :number of P rows of entire model
     &  extended_halo_EW,                                               &
                          ! IN  :EW halo size
     &  extended_halo_NS,                                               &
                          ! IN  :NS halo size
     &  tot_levels        ! IN  :total number of levels

! Parameters and Common blocks

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




! Local variables
      INTEGER                                                           &
     &  nproc_EW                                                        &
     &, nproc_NS                                                        &
     &, iproc                                                           &
     &, iproc_x                                                         &
     &, iproc_y                                                         &
     &, ifld                                                            &
     &, ihalo                                                           &
     &, idim                                                            &
     &, ipt                                                             &
     &, irest                                                           &
     &, start                                                           &
     &, size                                                            &
     &, prow_N                                                          &
     &, prow_S                                                          &
     &, info                                                            &
     &, in_atm_decomp

      PARAMETER (  nproc_EW=1                                           &
     &,            nproc_NS=1)

      LOGICAL                                                           &
     &  at_north

      INTEGER                                                           &
     &  size_x(0:nproc_EW-1)                                            &
     &, size_y(0:nproc_NS-1)                                            &
     &, start_x(0:nproc_EW-1)                                           &
     &, start_y(0:nproc_NS-1)

! Error reporting
      INTEGER       ErrorStatus ! =0 normal exit; >0 error exit
      CHARACTER*256 Cmessage    ! Error message
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='DECOMPOSE_SMEXE')

! ------------------------------------------------------------------
! 0.0 Check for valid decomposition
! ------------------------------------------------------------------

! No need to check for valid decomposition as nproc_EW, nproc_NS
! extended_halo_EW & extended_halo_NS are set

! ------------------------------------------------------------------
      decomp_db_halosize(1,halo_type_single,decomp_smexe) = 1
      decomp_db_halosize(2,halo_type_single,decomp_smexe) = 1
      decomp_db_halosize(3,halo_type_single,decomp_smexe) = 0

      decomp_db_halosize(1,halo_type_extended,decomp_smexe) =           &
     &  extended_halo_EW
      decomp_db_halosize(2,halo_type_extended,decomp_smexe) =           &
     &  extended_halo_NS
      decomp_db_halosize(3,halo_type_extended,decomp_smexe) =           &
     &  0

      decomp_db_halosize(1,halo_type_no_halo,decomp_smexe) = 0
      decomp_db_halosize(2,halo_type_no_halo,decomp_smexe) = 0
      decomp_db_halosize(3,halo_type_no_halo,decomp_smexe) = 0


! ------------------------------------------------------------------
! 1.0 Set up global data size
! ------------------------------------------------------------------

      decomp_db_glsize(1,fld_type_p,decomp_smexe) =                     &
     &  global_row_len
      decomp_db_glsize(2,fld_type_p,decomp_smexe) =                     &
     &  global_n_rows
      decomp_db_glsize(3,fld_type_p,decomp_smexe) =                     &
     &  tot_levels

      decomp_db_glsize(1,fld_type_u,decomp_smexe) =                     &
     &  global_row_len
      decomp_db_glsize(2,fld_type_u,decomp_smexe) =                     &
     &  global_n_rows
      decomp_db_glsize(3,fld_type_u,decomp_smexe) =                     &
     &  tot_levels

      decomp_db_glsize(1,fld_type_v,decomp_smexe) =                     &
     &  global_row_len
      decomp_db_glsize(2,fld_type_v,decomp_smexe) =                     &
     &  GLOBAL_N_ROWS
      decomp_db_glsize(3,fld_type_v,decomp_smexe) =                     &
     &  tot_levels

! ------------------------------------------------------------------
! 2.0 Calculate decomposition
! ------------------------------------------------------------------


! select processors to use for the data decomposition
      decomp_db_nproc(decomp_smexe)=nproc_EW*nproc_NS
      decomp_db_first_comp_pe(decomp_smexe) = 0
      decomp_db_last_comp_pe(decomp_smexe) =                            &
     &  decomp_db_nproc(decomp_smexe)-1

!     Set the grid size

      decomp_db_gridsize(1,decomp_smexe) = nproc_EW
      decomp_db_gridsize(2,decomp_smexe) = nproc_NS
      decomp_db_gridsize(3,decomp_smexe) = 1


! Work out the decomposition in the East-West direction. As far as
! possible each processor has the same local row length. However, if
! this is not possible, the extra points are distributed symetrically
! such that each processor has the same number of points as the
! processor on the opposite side of the globe.

      start=1
      size=global_row_len/nproc_EW  ! local data size on each processor
                                    ! assuming nproc_EW divides exactly
                                    ! into global_row_len.
      irest=global_row_len-(size*nproc_EW)
                                    ! If it doesn't divide exactly then
                                    ! irest contains the number of left
                                    ! over points that need to be
                                    ! allocated to processors

      DO iproc=1,nproc_EW
        start_x(iproc-1)=start

        IF (iproc  <=  nproc_EW/2) THEN
          ! processor in the first half of the row
          IF (iproc  <=  (irest/2)) THEN
            size_x(iproc-1)=size+1 ! gets one of the extra points
          ELSE
            size_x(iproc-1)=size   ! gets the standard row length
          ENDIF
        ELSE
          ! processor in the second half of the row
          IF (iproc-(nproc_EW/2)  <=  (irest/2)) THEN
            size_x(iproc-1)=size+1 ! gets one of the extra points
          ELSE
            size_x(iproc-1)=size   ! gets the standard row length
          ENDIF
        ENDIF

        start=start+size_x(iproc-1)

      ENDDO

! Work out the decomposition in the North-South direction. As far as
! possible each processor has the same number of rows. However, if this
! is not possible, the extra rows are distributed thus:
! - an extra row is given to the Northern most processor
! - the remaining extra rows are distributed symetrically around the
!   equator, starting at the processor(s) closest to the equator.

      start=1
      size=global_n_rows/nproc_NS  ! local data size on each processor
                                   ! assuming nproc_NS divides exactly
                                   ! into global_n_rows
      irest=global_n_rows-(size*nproc_NS)
                                   ! If it doesn't divide exactly then
                                   ! irest contains the number of left
                                   ! over points that need to be
                                   ! allocated to processors

      DO iproc=1,nproc_NS
        size_y(iproc-1)=size
      ENDDO

      IF (irest  >=  1) THEN
        ! give Northern most processors an extra row
        size_y(nproc_NS-1)=size+1
        irest=irest-1
      ENDIF

! Set up pointers to processor rows to which we will add extra rows
! to. These start around the equator, and will work out towards
! the poles.

      IF (MOD(nproc_NS,2)  ==  0) THEN  ! Even number of NS processors
        prow_S=nproc_NS/2
        prow_N=prow_S+1
      ELSE  ! Odd number of NS processors
        prow_S=(nproc_NS/2)+1
        prow_N=prow_S
      ENDIF

      DO WHILE (irest  >=  1)

        IF (prow_N  ==  prow_S) THEN
          size_y(prow_N-1)=size+1
          irest=irest-1
        ELSE
          size_y(prow_S-1)=size+1
          irest=irest-1
          IF (irest  >=  1) THEN
            size_y(prow_N-1)=size+1
            irest=irest-1
          ENDIF
        ENDIF

        prow_S=MAX(1,prow_S-1)
        prow_N=MIN(nproc_NS,prow_N+1)

      ENDDO

      DO iproc=1,nproc_NS
        start_y(iproc-1)=start
        start=start+size_y(iproc-1)
      ENDDO

! Set the local data shape and offsets of each processor

      DO iproc_y=0,nproc_NS-1

        IF (iproc_y  ==  (nproc_NS-1)) THEN
          at_north=.TRUE.  ! Doing the nothernmost processors
        ELSE
         at_north=.FALSE. ! Not doing the northernmost processors
        ENDIF

        DO iproc_x=0,nproc_EW-1

          iproc=decomp_db_first_comp_pe(decomp_smexe)+                  &
     &          iproc_x+(iproc_y*nproc_EW)

! Set the position in the Logical Processor Grid

          decomp_db_g_gridpos(1,iproc,decomp_smexe)=iproc_x
          decomp_db_g_gridpos(2,iproc,decomp_smexe)=iproc_y
          decomp_db_g_gridpos(3,iproc,decomp_smexe)=0

! Set the number of local datapoints (blsize) on the processor

!  Fields on P grid:
          decomp_db_g_blsize(1,fld_type_p,iproc,                        &
     &                       decomp_smexe)=                             &
     &      size_x(iproc_x)
          decomp_db_g_blsize(2,fld_type_p,iproc,                        &
     &                       decomp_smexe)=                             &
     &      size_y(iproc_y)
          decomp_db_g_blsize(3,fld_type_p,iproc,                        &
     &                       decomp_smexe)=                             &
     &      tot_levels

! Fields on U grid:
          decomp_db_g_blsize(1,fld_type_u,iproc,                        &
     &                       decomp_smexe)=                             &
     &      size_x(iproc_x)
          decomp_db_g_blsize(2,fld_type_u,iproc,                        &
     &                       decomp_smexe)=                             &
     &      size_y(iproc_y)
          decomp_db_g_blsize(3,fld_type_u,iproc,                        &
     &                       decomp_smexe)=                             &
     &      tot_levels

! Fields on V grid:
          decomp_db_g_blsize(1,fld_type_v,iproc,                        &
     &                       decomp_smexe)=                             &
     &      size_x(iproc_x)
          IF (at_north) THEN
            decomp_db_g_blsize(2,fld_type_v,iproc,                      &
     &                         decomp_smexe)=                           &
     &      size_y(iproc_y)
          ELSE
            decomp_db_g_blsize(2,fld_type_v,iproc,                      &
     &                         decomp_smexe)=                           &
     &      size_y(iproc_y)
          ENDIF
          decomp_db_g_blsize(3,fld_type_v,iproc,                        &
     &                       decomp_smexe)=                             &
     &      tot_levels

! Set the number of points including the halos on the processor

          DO ihalo=1,NHalo_max
            DO ifld=1,Nfld_max
              DO idim=1,Ndim_max
                decomp_db_g_lasize(idim,ifld,ihalo,iproc,               &
     &                             decomp_smexe)=                       &
     &            decomp_db_g_blsize(idim,ifld,iproc,                   &
     &                               decomp_smexe)+                     &
     &            2*decomp_db_halosize(idim,ihalo,                      &
     &                               decomp_smexe)
              ENDDO  ! idim
            ENDDO  ! ifld
          ENDDO  ! ihalo

! Set the starting point in the global domain

          decomp_db_g_datastart(1,iproc,decomp_smexe)=                  &
     &      start_x(iproc_x)
          decomp_db_g_datastart(2,iproc,decomp_smexe)=                  &
     &      start_y(iproc_y)
          decomp_db_g_datastart(3,iproc,decomp_smexe)=1

          decomp_db_g_datastart_f(1,fld_type_p,iproc,                   &
     &                decomp_smexe)=start_x(iproc_x)
          decomp_db_g_datastart_f(2,fld_type_p,iproc,                   &
     &                decomp_smexe)=start_y(iproc_y)
          decomp_db_g_datastart_f(3,fld_type_p,iproc,                   &
     &                decomp_smexe)=1

          decomp_db_g_datastart_f(1,fld_type_u,iproc,                   &
     &                decomp_smexe)=start_x(iproc_x)
          decomp_db_g_datastart_f(2,fld_type_u,iproc,                   &
     &                decomp_smexe)=start_y(iproc_y)
          decomp_db_g_datastart_f(3,fld_type_u,iproc,                   &
     &                decomp_smexe)=1

          decomp_db_g_datastart_f(1,fld_type_v,iproc,                   &
     &                decomp_smexe)=start_x(iproc_x)
          decomp_db_g_datastart_f(2,fld_type_v,iproc,                   &
     &                decomp_smexe)=start_y(iproc_y)
          decomp_db_g_datastart_f(3,fld_type_v,iproc,                   &
     &                decomp_smexe)=1

          decomp_db_g_datastart_f(1,fld_type_r,iproc,                   &
     &                decomp_smexe)=0
          decomp_db_g_datastart_f(2,fld_type_r,iproc,                   &
     &                decomp_smexe)=0
          decomp_db_g_datastart_f(3,fld_type_r,iproc,                   &
     &                decomp_smexe)=1
        ENDDO ! iproc_x
      ENDDO ! iproc_y

! Set up the pe_index_EW array - for each point along a global row
! it indicates the PE index (along the processor row) which
! contains that point

      DO iproc_x=0,nproc_EW-1
        DO ipt=decomp_db_g_datastart(1,iproc_x,decomp_smexe),           &
     &         decomp_db_g_datastart(1,iproc_x,decomp_smexe)+           &
     &         size_x(iproc_x)
          decomp_db_g_pe_index_EW(ipt,decomp_smexe)=iproc_x
        ENDDO
      ENDDO

! And fill in the halos at either end

      DO ipt=1-extended_halo_EW,1
        decomp_db_g_pe_index_EW(ipt,decomp_smexe)=0
      ENDDO

      DO ipt=global_row_len+1,global_row_len+1+extended_halo_EW
        decomp_db_g_pe_index_EW(ipt,decomp_smexe)=nproc_x-1
      ENDDO

! Now set up the pe_index_NS_array - for each point along a global
! North-South column it indicates the PE index (along the processor
! column) which contains that point

      DO iproc_y=0,nproc_NS-1
        DO ipt=decomp_db_g_datastart(2,iproc_y*nproc_EW,                &
     &                               decomp_smexe),                     &
     &         decomp_db_g_datastart(2,iproc_y*nproc_EW,                &
     &                               decomp_smexe)+                     &
     &         size_y(iproc_y)
          decomp_db_g_pe_index_NS(ipt,decomp_smexe)=iproc_y
        ENDDO
      ENDDO

! And fill in the halos at either end

      DO ipt=1-extended_halo_NS,1
        decomp_db_g_pe_index_NS(ipt,decomp_smexe)=0
      ENDDO

      DO ipt=global_n_rows+1,global_n_rows+1+extended_halo_NS
        decomp_db_g_pe_index_NS(ipt,decomp_smexe)=nproc_y-1
      ENDDO

! ------------------------------------------------------------------
! 3.0 Return the new data sizes and exit subroutine
! ------------------------------------------------------------------

! Set up the GCOM groups:

! 1) Group of all processors on my row

      IF ( decomp_db_gridsize(2,decomp_smexe)  ==  1)                   &
     & THEN
       decomp_db_gc_proc_row_group(decomp_smexe)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,                                  &
     &    decomp_db_g_gridpos(2,mype,decomp_smexe),                     &
     &    info,                                                         &
     &    decomp_db_gc_proc_row_group(decomp_smexe))
      ENDIF

! 2) Group of all processors on my column

      IF ( decomp_db_gridsize(1,decomp_smexe)  ==  1)                   &
     & THEN
        decomp_db_gc_proc_col_group(decomp_smexe)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,                                  &
     &    decomp_db_g_gridpos(1,mype,decomp_smexe),                     &
     &    info,                                                         &
     &    decomp_db_gc_proc_col_group(decomp_smexe))
      ENDIF

! 3) Group of all processors in the atmosphere model
      IF (decomp_db_nproc(decomp_smexe)  ==  nproc_max)                 &
     & THEN
        decomp_db_gc_all_proc_group(decomp_smexe)=GCG_ALL
      ELSE
        IF ((mype  >=  decomp_db_first_comp_pe(decomp_smexe))           &
     &    .AND.                                                         &
     &     (mype  <=  decomp_db_last_comp_pe(decomp_smexe)) )           &
     &   THEN
          in_atm_decomp=1
        ELSE
          in_atm_decomp=0
        ENDIF

        CALL GCG_SPLIT(mype,nproc_max,in_atm_decomp,info,               &
     &    decomp_db_gc_all_proc_group(decomp_smexe))
      ENDIF

! Set logical indicating this decomposition has been initialised
! and is now ready for use

      decomp_db_set(decomp_smexe)=.TRUE.

      RETURN

      END SUBROUTINE DECOMPOSE_SMEXE
