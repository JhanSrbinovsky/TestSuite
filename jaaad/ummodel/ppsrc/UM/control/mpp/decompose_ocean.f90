

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel UM: Perform data decomposition for ocean model
!
! Subroutine Interface:
      SUBROUTINE DECOMPOSE_OCEAN(global_row_len,global_n_rows,          &
     &                           tot_levels,model_type,                 &
     &                           nproc_EW, nproc_NS,                    &
     &                           extended_halo_EW,                      &
     &                           extended_halo_NS,                      &
     &                           local_row_len,local_n_rows,            &
     &                           l_ocyclic)
      IMPLICIT NONE
!
! Desciption:
! This routine currently performs a 1D North-South decomposition on the
! ocean model. nproc_EW is currently ignored.
! The decomposition strategy is much the same as the atmosphere's -
! First try and divide the rows equally between processors, and then
! distribute any left over rows to the processors, starting from the
! top.
!
! Method:
! The local data sizes are calculated and sotred in the COMMON block
! DECOMPDB. The boundary conditions are set (cyclic in East/West
! direction if *DEF,global
!
! Current Code Owner : Paul Burton
!
! History:
!  Model    Date     Modification history from model version 4.2
!  version
!  4.2      21/8/96  New deck created for mpp ocean model.  P.Burton
!  4.3      17/02/97 Added initialisation of new ocean decompositon
!                    decomp_nowrap_ocean - which does not include
!                    the wrap around points on the ends of rows.
!                    This requires passing in the l_ocyclic variable
!                    to indicate if these points are present.
!                                                         P.Burton
!  5.0      09/08/99 Reflect changes made to Atmosphere version of
!                    this routine in order to cater for mpp control
!                    code changes introduced in connection with the
!                    new atmosphere dynamics developments.
!                                            R. Hill
!  5.1      27/01/99 Added g_pe_index_EW/NS arrays    P.Burton
!  5.1      19/04/00  Correct lasize for velocity grid.
!                     Include set up of lasize values for
!                     halo_type_no_halo conditions. R. Hill
!  5.3      14/09/01 Added model_domain variable    P.Burton
!  5.3               Remove restriction on odd no of
!                    points E-W. That's not an ocean
!                    limitation.   R. Hill
!  5.4    Apr 2002      Fix bug in "extended halo type". Although
!                       not specifically used by the Ocean, these
!                       values are required in certain generic STASH
!                       routines when processing sub-areas.
!                       Also add g_pe_index_EW/NS arrays for non-wrap
!                       grids to re-enable sub-domain STASH extraction
!                       overlooked by 5.1 mpp developments.    R. Hill
!  5.5      25/03/03 Setup of datastart_f (field based datastart).
!                    P.Selwood.
!  6.0  17/09/03  Add def for new NEC opt section c96_1c. R Barnes
!
! Subroutine arguments:

      INTEGER                                                           &

     &  global_row_len                                                  &
                        ! IN :  number of E-W points of entire model
     &, global_n_rows                                                   &
                        ! IN :  number of N-S points of entire model
     &, tot_levels                                                      &
                        ! IN :  total number of levels
     &, model_type                                                      &
                        ! IN  : type (Global,LAM etc) of model
     &, nproc_EW                                                        &
                        ! IN :  number of processors to decompose E-W
     &, nproc_NS                                                        &
                        ! IN :  number of processors to decompose N-S
     &, extended_halo_EW                                                &
                         ! IN  : size of extended EW halo
     &, extended_halo_NS                                                &
                         ! IN  : size of extended NS halo
     &, local_row_len                                                   &
                        ! OUT : local number of E-W points
     &, local_n_rows    ! OUT : local number of N-S points
!                       ! local_row_len and local_n_rows include
!                       ! any halos

      LOGICAL                                                           &

     &  l_ocyclic       ! IN : true if extra wrap points are present
!                       !      at the start/ends of rows

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

! Local variables
      INTEGER                                                           &
     &  iproc                                                           &
     &, iproc_x                                                         &
     &, iproc_y                                                         &
     &, ifld                                                            &
     &, ihalo                                                           &
     &, idim                                                            &
     &, ipt                                                             &
     &, irest                                                           &
     &, jrest                                                           &
     &, start                                                           &
     &, size                                                            &
     &, prow_N                                                          &
     &, prow_S                                                          &
     &, info                                                            &
     &, in_ocn_decomp

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
      PARAMETER (   RoutineName='DECOMPOSE_OCEAN')


! ------------------------------------------------------------------
! 0.0 Check for valid decomposition
! ------------------------------------------------------------------

      IF (nproc_EW*nproc_NS  >   Maxproc) THEN
        ErrorStatus=1
        WRITE(Cmessage,                                                 &
     &    '("Cannot run with decomposition ",I3," x ",I3,               &
     &      " (",I3,") processors. ",                                   &
     &      "Maxproc is ",I3," processors.")')                          &
     &       nproc_EW,nproc_NS,nproc_EW*nproc_NS,Maxproc
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      IF ((nproc_EW  /=  1) .AND. (MOD(nproc_EW,2)  /=  0)) THEN
        ErrorStatus=2
        WRITE(Cmessage,                                                 &
     &    '("Cannot run with an odd (",I3,") number of processors ",    &
     &      "in the East-West direction.")') nproc_EW
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

!      No restriction on odd numbers of E-W points

      IF (extended_halo_EW  >   Max_Halo_Size) THEN
        ErrorStatus=4
        WRITE(Cmessage,                                                 &
     &    '("East-West extended halo size (",I2,") is too large.",      &
     &      "The maximum permitted size is Max_Halo_Size=",I2)')        &
     &    extended_halo_EW,Max_Halo_Size
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      IF (extended_halo_NS  >   Max_Halo_Size) THEN
        ErrorStatus=4
        WRITE(Cmessage,                                                 &
     &    '("North-South extended halo size (",I2,") is too large.",    &
     &      "The maximum permitted size is Max_Halo_Size=",I2)')        &
     &    extended_halo_NS,Max_Halo_Size
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF
! ------------------------------------------------------------------
      decomp_db_sb_model_domain(decomp_standard_ocean)=model_type

! Halo Sizes
      ! Reminder: The Ocean only decomposes in the S-N direction so
      ! no halos are required E-W.
      decomp_db_halosize(1,halo_type_single,decomp_standard_ocean) = 0
      decomp_db_halosize(2,halo_type_single,decomp_standard_ocean) = 1
      decomp_db_halosize(3,halo_type_single,decomp_standard_ocean) = 0

      ! As far as the ocean is concerned, we dont worry about
      ! extended halo sizes. Just use the standard settings.
      decomp_db_halosize(1,halo_type_extended,decomp_standard_ocean) =  &
     &  0
      decomp_db_halosize(2,halo_type_extended,decomp_standard_ocean) =  &
     &  1
      decomp_db_halosize(3,halo_type_extended,decomp_standard_ocean) =  &
     &  0

      ! As far as the ocean is concerned, no_halo means exactly that
      ! ... set accordingly.
      decomp_db_halosize(1,halo_type_no_halo,decomp_standard_ocean) = 0
      decomp_db_halosize(2,halo_type_no_halo,decomp_standard_ocean) = 0
      decomp_db_halosize(3,halo_type_no_halo,decomp_standard_ocean) = 0

! Size of global data

      decomp_db_glsize(1,fld_type_p,decomp_standard_ocean) =            &
     &  global_row_len
      decomp_db_glsize(2,fld_type_p,decomp_standard_ocean) =            &
     &  global_n_rows
      decomp_db_glsize(3,fld_type_p,decomp_standard_ocean) =            &
     &  tot_levels

      decomp_db_glsize(1,fld_type_u,decomp_standard_ocean) =            &
     &  global_row_len
      decomp_db_glsize(2,fld_type_u,decomp_standard_ocean) =            &
     &  global_n_rows-1
      decomp_db_glsize(3,fld_type_u,decomp_standard_ocean) =            &
     &  tot_levels

      decomp_db_glsize(1,fld_type_v,decomp_standard_ocean) =            &
     &  global_row_len
      decomp_db_glsize(2,fld_type_v,decomp_standard_ocean) =            &
     &  global_n_rows-1
      decomp_db_glsize(3,fld_type_v,decomp_standard_ocean) =            &
     &  tot_levels

! Make sure there's actually enough work for all the processors to do

      IF (nproc_NS  >   global_n_rows) THEN
        IF (mype  ==  0) THEN
          WRITE(6,*) 'Warning : Ocean model has more processors than ', &
     &               'rows. Reducing nproc_y to ',global_n_rows
        ENDIF
        nproc_NS=global_n_rows
      ENDIF

      decomp_db_nproc(decomp_standard_ocean)=nproc_NS
      decomp_db_first_comp_pe(decomp_standard_ocean) = 0
      decomp_db_last_comp_pe(decomp_standard_ocean) =                   &
     &  decomp_db_nproc(decomp_standard_ocean)-1

! Set the size of the Logical Processor Grid (LPG)

      decomp_db_gridsize(1,decomp_standard_ocean) = nproc_EW  ! =1
      decomp_db_gridsize(2,decomp_standard_ocean) = nproc_NS
      decomp_db_gridsize(3,decomp_standard_ocean) = 1

! Calculate processor specific information.

      DO iproc=decomp_db_first_comp_pe(decomp_standard_ocean),          &
     &         decomp_db_last_comp_pe(decomp_standard_ocean)
!       ! Loop over all processors in this decomposition

! NB : Although the decomposition is currently only N-S, all
! the code is included to allow an E-W decomposition too.
! All that is required is to supply nproc_NS > 1.

! Calculate the position in the LPG:
        decomp_db_g_gridpos(3,iproc,decomp_standard_ocean) = 0
        decomp_db_g_gridpos(2,iproc,decomp_standard_ocean) =            &
     &    iproc / decomp_db_gridsize(1,decomp_standard_ocean)
        decomp_db_g_gridpos(1,iproc,decomp_standard_ocean) =            &
     &    iproc - decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)*   &
     &            decomp_db_gridsize(1,decomp_standard_ocean)

! Calculate the local data sizes for processor iproc

! East-West decomposition

        decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)=   &
     &    decomp_db_glsize(1,fld_type_p,decomp_standard_ocean) /        &
     &    decomp_db_gridsize(1,decomp_standard_ocean)

        decomp_db_g_blsize(1,fld_type_u,iproc,decomp_standard_ocean)=   &
     &    decomp_db_glsize(1,fld_type_u,decomp_standard_ocean) /        &
     &    decomp_db_gridsize(1,decomp_standard_ocean)

        decomp_db_g_blsize(1,fld_type_v,iproc,decomp_standard_ocean)=   &
     &    decomp_db_glsize(1,fld_type_v,decomp_standard_ocean) /        &
     &    decomp_db_gridsize(1,decomp_standard_ocean)


        irest = decomp_db_glsize(1,fld_type_p,decomp_standard_ocean)-   &
     &  decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)*   &
     &          decomp_db_gridsize(1,decomp_standard_ocean)


        decomp_db_g_datastart(1,iproc,decomp_standard_ocean) =          &
     &    decomp_db_g_gridpos(1,iproc,decomp_standard_ocean)*           &
     & decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)+1

        IF (decomp_db_g_gridpos(1,iproc,decomp_standard_ocean)  <       &
     &      irest) THEN

          decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)= &
     &   decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)+1

          decomp_db_g_datastart(1,iproc,decomp_standard_ocean) =        &
     &      decomp_db_g_datastart(1,iproc,decomp_standard_ocean) +      &
     &      decomp_db_g_gridpos(1,iproc,decomp_standard_ocean)

        ELSE
          decomp_db_g_datastart(1,iproc,decomp_standard_ocean) =        &
     &      decomp_db_g_datastart(1,iproc,decomp_standard_ocean) +      &
     &      irest
        ENDIF

        decomp_db_g_datastart_f(1,fld_type_p,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(1,iproc,decomp_standard_ocean)
        decomp_db_g_datastart_f(1,fld_type_u,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(1,iproc,decomp_standard_ocean)
        decomp_db_g_datastart_f(1,fld_type_v,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(1,iproc,decomp_standard_ocean)

        decomp_db_g_lasize(1,fld_type_p,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean) +  &
     &  2*decomp_db_halosize(1,halo_type_single,decomp_standard_ocean)
        ! U grid, imt values
        decomp_db_g_lasize(1,fld_type_u,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &   decomp_db_g_lasize(1,fld_type_p,halo_type_single,              &
     &                     iproc,decomp_standard_ocean)

        ! V grid, imt values
        decomp_db_g_lasize(1,fld_type_v,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_lasize(1,fld_type_p,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)


!-----------------------------------------
! Set up column-wise lasize for no_halo conditions
!-----------------------------------------
        decomp_db_g_lasize(1,fld_type_p,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean) +  &
     & 2*decomp_db_halosize(1,halo_type_no_halo,decomp_standard_ocean)

        ! U grid, imt values
        decomp_db_g_lasize(1,fld_type_u,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &   decomp_db_g_lasize(1,fld_type_p,halo_type_no_halo,             &
     &                     iproc,decomp_standard_ocean)

        ! V grid, imt values
        decomp_db_g_lasize(1,fld_type_v,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_lasize(1,fld_type_p,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)



! North-South decomposition

        decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean) =  &
     &    decomp_db_glsize(2,fld_type_p,decomp_standard_ocean) /        &
     &    decomp_db_gridsize(2,decomp_standard_ocean)

        jrest = decomp_db_glsize(2,fld_type_p,decomp_standard_ocean)-   &
     &   decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)*  &
     &          decomp_db_gridsize(2,decomp_standard_ocean)

        decomp_db_g_datastart(2,iproc,decomp_standard_ocean) =          &
     &    decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)*           &
     &  decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)+1

        IF (decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)  <       &
     &      jrest) THEN

          decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)= &
     &  decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)+1

          decomp_db_g_datastart(2,iproc,decomp_standard_ocean) =        &
     &      decomp_db_g_datastart(2,iproc,decomp_standard_ocean) +      &
     &      decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)

        ELSE
          decomp_db_g_datastart(2,iproc,decomp_standard_ocean) =        &
     &      decomp_db_g_datastart(2,iproc,decomp_standard_ocean) +      &
     &      jrest
        ENDIF

        decomp_db_g_datastart_f(2,fld_type_p,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(2,iproc,decomp_standard_ocean)
        decomp_db_g_datastart_f(2,fld_type_u,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(2,iproc,decomp_standard_ocean)
        decomp_db_g_datastart_f(2,fld_type_v,iproc,                     &
     &                          decomp_standard_ocean) =                &
     &    decomp_db_g_datastart(2,iproc,decomp_standard_ocean)

        decomp_db_g_lasize(2,fld_type_p,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean) +  &
     &    2*decomp_db_halosize(2,halo_type_single,decomp_standard_ocean)


        ! U velocity grid in the S-N direction.
        decomp_db_g_lasize(2,fld_type_u,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &     decomp_db_g_lasize(2,fld_type_p,halo_type_single,            &
     &                     iproc,decomp_standard_ocean)

        ! V velocity grid in the S-N direction.
        decomp_db_g_lasize(2,fld_type_v,halo_type_single,               &
     &                     iproc,decomp_standard_ocean)=                &
     &     decomp_db_g_lasize(2,fld_type_p,halo_type_single,            &
     &                     iproc,decomp_standard_ocean)
!-----------------------------------------
! Set up row-wise lasize for no_halo conditions
!-----------------------------------------

        decomp_db_g_lasize(2,fld_type_p,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &  decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean) +  &
     & 2*decomp_db_halosize(2,halo_type_no_halo,decomp_standard_ocean)


        ! U velocity grid in the S-N direction.
        decomp_db_g_lasize(2,fld_type_u,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &     decomp_db_g_lasize(2,fld_type_p,halo_type_no_halo,           &
     &                     iproc,decomp_standard_ocean)

        ! V velocity grid in the S-N direction.
        decomp_db_g_lasize(2,fld_type_v,halo_type_no_halo,              &
     &                     iproc,decomp_standard_ocean)=                &
     &     decomp_db_g_lasize(2,fld_type_p,halo_type_no_halo,           &
     &                     iproc,decomp_standard_ocean)
! No decomposition in the vertical

        decomp_db_g_datastart(3,iproc,decomp_standard_ocean) = 1
        decomp_db_g_datastart_f(3,fld_type_p,iproc,                     &
     &                          decomp_standard_ocean) = 1
        decomp_db_g_datastart_f(3,fld_type_u,iproc,                     &
     &                          decomp_standard_ocean) = 1
        decomp_db_g_datastart_f(3,fld_type_v,iproc,                     &
     &                          decomp_standard_ocean) = 1
        decomp_db_g_blsize(3,fld_type_p,iproc,decomp_standard_ocean)=   &
     &    tot_levels
        DO ihalo=1,NHalo_max
           DO ifld=1,Nfld_max
              decomp_db_g_lasize(3,ifld,ihalo,                          &
     &                iproc,decomp_standard_ocean) = tot_levels
           ENDDO
        ENDDO



! One less U/V row at North. lasize does not reflect this - velocity
! grid sizes are the same as tracer grid sizes, row-wise for mpp
! purposes (mainly in STASH). Effectively, lasize on the northern-most
! PE is over-dimensioned by 1 row.

        decomp_db_g_blsize(1,fld_type_u,iproc,decomp_standard_ocean) =  &
     &    decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)

        decomp_db_g_blsize(1,fld_type_v,iproc,decomp_standard_ocean) =  &
     &    decomp_db_g_blsize(1,fld_type_p,iproc,decomp_standard_ocean)

        IF (  decomp_db_g_gridpos(2,iproc,decomp_standard_ocean)        &
     &   ==  (decomp_db_gridsize(2,decomp_standard_ocean)-1)) THEN

          decomp_db_g_blsize(2,fld_type_u,iproc,decomp_standard_ocean)= &
     &    decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)-1

          decomp_db_g_blsize(2,fld_type_v,iproc,decomp_standard_ocean)= &
     &    decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)-1


        ELSE

          decomp_db_g_blsize(2,fld_type_u,iproc,decomp_standard_ocean)= &
     &    decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)

          decomp_db_g_blsize(2,fld_type_v,iproc,decomp_standard_ocean)= &
     &    decomp_db_g_blsize(2,fld_type_p,iproc,decomp_standard_ocean)


        ENDIF

        decomp_db_g_blsize(3,fld_type_u,iproc,decomp_standard_ocean) =  &
     &    decomp_db_g_blsize(3,fld_type_p,iproc,decomp_standard_ocean)

        decomp_db_g_blsize(3,fld_type_v,iproc,decomp_standard_ocean) =  &
     &    decomp_db_g_blsize(3,fld_type_p,iproc,decomp_standard_ocean)

      ENDDO  ! loop over processors

! Set up the pe_index_EW_array - for each point along a global
! row it indicates the PE index (along the processor
! row) which contains that point

      DO iproc_x=0,nproc_EW-1
        DO ipt=decomp_db_g_datastart(1,iproc_x,decomp_standard_ocean),  &
     &         decomp_db_g_datastart(1,iproc_x,decomp_standard_ocean)+  &
     &         decomp_db_g_blsize(1,fld_type_p,iproc_x,                 &
     &                            decomp_standard_ocean)
          decomp_db_g_pe_index_EW(ipt,decomp_standard_ocean)=iproc_x
          decomp_db_g_pe_index_EW(ipt,decomp_nowrap_ocean)=iproc_x
        ENDDO
      ENDDO

! And fill in the halos at either end

      DO ipt=1-extended_halo_EW,1
        decomp_db_g_pe_index_EW(ipt,decomp_standard_ocean)=0
        decomp_db_g_pe_index_EW(ipt,decomp_nowrap_ocean)=0
      ENDDO

      DO ipt=global_row_len+1,global_row_len+1+extended_halo_EW
        decomp_db_g_pe_index_EW(ipt,decomp_standard_ocean)=nproc_x-1
        decomp_db_g_pe_index_EW(ipt,decomp_nowrap_ocean)=nproc_x-1
      ENDDO

! Now set up the pe_index_NS_array - for each point along a global
! North-South column it indicates the PE index (along the processor
! column) which contains that point

      DO iproc_y=0,nproc_NS-1
        DO ipt=decomp_db_g_datastart(2,iproc_y*nproc_EW,                &
     &                               decomp_standard_ocean),            &
     &         decomp_db_g_datastart(2,iproc_y*nproc_EW,                &
     &                               decomp_standard_ocean)+            &
     &         decomp_db_g_blsize(2,fld_type_p,iproc_y*nproc_EW,        &
     &                            decomp_standard_ocean)
          decomp_db_g_pe_index_NS(ipt,decomp_standard_ocean)=iproc_y
          decomp_db_g_pe_index_NS(ipt,decomp_nowrap_ocean)=iproc_y
        ENDDO
      ENDDO

! And fill in the halos at either end

      DO ipt=1-extended_halo_NS,1
        decomp_db_g_pe_index_NS(ipt,decomp_standard_ocean)=0
        decomp_db_g_pe_index_NS(ipt,decomp_nowrap_ocean)=0
      ENDDO

      DO ipt=global_n_rows+1,global_n_rows+1+extended_halo_NS
        decomp_db_g_pe_index_NS(ipt,decomp_standard_ocean)=nproc_y-1
        decomp_db_g_pe_index_NS(ipt,decomp_nowrap_ocean)=nproc_y-1
      ENDDO
! Set up the boundary types

      decomp_db_bound(1,decomp_standard_ocean) = BC_CYCLIC
!       ! Cyclic East-West boundaries
      decomp_db_bound(2,decomp_standard_ocean) = BC_STATIC
!       ! No North-South wrap around
      decomp_db_bound(3,decomp_standard_ocean) = BC_STATIC
!       ! No vertical wrap around

! And set up the neighbour array

! DEPENDS ON: set_neighbour
      CALL SET_NEIGHBOUR(                                               &
     &  decomp_standard_ocean)

! Set up the GCOM groups

! 1) Group of all processors on my row

      IF ( decomp_db_gridsize(2,decomp_standard_ocean)  ==  1)          &
     & THEN
       decomp_db_gc_proc_row_group(decomp_standard_ocean)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,                                  &
     &    decomp_db_g_gridpos(2,mype,decomp_standard_ocean),            &
     &    info,                                                         &
     &    decomp_db_gc_proc_row_group(decomp_standard_ocean))
      ENDIF

! 2) Group of all processors on my column

      IF ( decomp_db_gridsize(1,decomp_standard_ocean)  ==  1)          &
     & THEN
        decomp_db_gc_proc_col_group(decomp_standard_ocean)=GCG_ALL
      ELSE
        CALL GCG_SPLIT(mype,nproc_max,                                  &
     &    decomp_db_g_gridpos(1,mype,decomp_standard_ocean),            &
     &    info,                                                         &
     &    decomp_db_gc_proc_col_group(decomp_standard_ocean))
      ENDIF

! 3) Group of all processors in the atmosphere model
      IF (decomp_db_nproc(decomp_standard_ocean)  ==  nproc_max)        &
     & THEN
        decomp_db_gc_all_proc_group(decomp_standard_ocean)=GCG_ALL
      ELSE
        IF ((mype  >=  decomp_db_first_comp_pe(decomp_standard_ocean))  &
     &    .AND.                                                         &
     &     (mype  <=  decomp_db_last_comp_pe(decomp_standard_ocean) ))  &
     &  THEN
          in_ocn_decomp=1
        ELSE
          in_ocn_decomp=0
        ENDIF

        CALL GCG_SPLIT(mype,nproc_max,in_ocn_decomp,info,               &
     &    decomp_db_gc_all_proc_group(decomp_standard_ocean))
      ENDIF

! Set logical indicating this decomposition has been initialised
! and is now ready for use

      decomp_db_set(decomp_standard_ocean)=.TRUE.

! Initialise decomp_nowrap_ocean which doesn't contain extra wrap
! points at start and end of each row
! Mostly it is a straight copy of the original ocean decomposition

      decomp_db_sb_model_domain(decomp_nowrap_ocean)=                   &
     &  decomp_db_sb_model_domain(decomp_standard_ocean)
      DO idim=1,Ndim_max
        decomp_db_bound(idim,decomp_nowrap_ocean)=                      &
     &    decomp_db_bound(idim,decomp_standard_ocean)
        decomp_db_glsize(idim,fld_type_p,decomp_nowrap_ocean)=          &
     &    decomp_db_glsize(idim,fld_type_p,decomp_standard_ocean)
        decomp_db_glsize(idim,fld_type_u,decomp_nowrap_ocean)=          &
     &    decomp_db_glsize(idim,fld_type_u,decomp_standard_ocean)
        decomp_db_glsize(idim,fld_type_v,decomp_nowrap_ocean)=          &
     &    decomp_db_glsize(idim,fld_type_v,decomp_standard_ocean)


        decomp_db_gridsize(idim,decomp_nowrap_ocean)=                   &
     &    decomp_db_gridsize(idim,decomp_standard_ocean)
        decomp_db_halosize(idim,halo_type_single,decomp_nowrap_ocean)=  &
     &  decomp_db_halosize(idim,halo_type_single,decomp_standard_ocean)
      ENDDO

      DO iproc=decomp_db_first_comp_pe(decomp_standard_ocean),          &
     &         decomp_db_last_comp_pe(decomp_standard_ocean)
         DO ihalo=1,NHalo_max
            DO ifld=1,Nfld_max
               DO idim=1,Ndim_max
                  decomp_db_g_lasize(idim,ifld,ihalo                    &
     &                       ,iproc,decomp_nowrap_ocean)=               &
     &                    decomp_db_g_lasize(idim,ifld,ihalo            &
     &                         ,iproc,decomp_standard_ocean)
               ENDDO
            ENDDO
         ENDDO

         DO ifld=1,Nfld_max
            DO idim=1,Ndim_max
             decomp_db_g_blsize(idim,ifld,iproc,decomp_nowrap_ocean)=   &
     &       decomp_db_g_blsize(idim,ifld,iproc,decomp_standard_ocean)
            ENDDO
         ENDDO


         DO idim=1,Ndim_max
            decomp_db_g_datastart(idim,iproc,decomp_nowrap_ocean)=      &
     &      decomp_db_g_datastart(idim,iproc,decomp_standard_ocean)

           decomp_db_g_datastart_f(idim,fld_type_p,iproc,               &
     &                             decomp_nowrap_ocean)=                &
     &      decomp_db_g_datastart_f(idim,fld_type_p,iproc,              &
     &                             decomp_standard_ocean)
           decomp_db_g_datastart_f(idim,fld_type_u,iproc,               &
     &                             decomp_nowrap_ocean)=                &
     &      decomp_db_g_datastart_f(idim,fld_type_u,iproc,              &
     &                             decomp_standard_ocean)
           decomp_db_g_datastart_f(idim,fld_type_v,iproc,               &
     &                             decomp_nowrap_ocean)=                &
     &      decomp_db_g_datastart_f(idim,fld_type_v,iproc,              &
     &                             decomp_standard_ocean)

            decomp_db_g_gridpos(idim,iproc,decomp_nowrap_ocean)=        &
     &      decomp_db_g_gridpos(idim,iproc,decomp_standard_ocean)
         ENDDO
      ENDDO

      DO idim=1,4
        decomp_db_neighbour(idim,decomp_nowrap_ocean)=                  &
     &    decomp_db_neighbour(idim,decomp_standard_ocean)
      ENDDO

      decomp_db_first_comp_pe(decomp_nowrap_ocean)=                     &
     &  decomp_db_first_comp_pe(decomp_standard_ocean)
      decomp_db_last_comp_pe(decomp_nowrap_ocean)=                      &
     &  decomp_db_last_comp_pe(decomp_standard_ocean)
      decomp_db_nproc(decomp_nowrap_ocean)=                             &
     &  decomp_db_nproc(decomp_standard_ocean)
      decomp_db_gc_proc_row_group(decomp_nowrap_ocean)=                 &
     &  decomp_db_gc_proc_row_group(decomp_standard_ocean)
      decomp_db_gc_proc_col_group(decomp_nowrap_ocean)=                 &
     &  decomp_db_gc_proc_col_group(decomp_standard_ocean)
      decomp_db_gc_all_proc_group(decomp_nowrap_ocean)=                 &
     &  decomp_db_gc_all_proc_group(decomp_standard_ocean)

      IF (l_ocyclic) THEN
! Make modifications to the decompositions to remove the point at
! the beginning and end of each row
        decomp_db_glsize(1,fld_type_p,decomp_nowrap_ocean)=             &
     &    decomp_db_glsize(1,fld_type_p,decomp_nowrap_ocean)-2

        decomp_db_glsize(1,fld_type_u,decomp_nowrap_ocean)=             &
     &    decomp_db_glsize(1,fld_type_u,decomp_nowrap_ocean)-2

        decomp_db_glsize(1,fld_type_v,decomp_nowrap_ocean)=             &
     &    decomp_db_glsize(1,fld_type_v,decomp_nowrap_ocean)-2

      DO iproc=decomp_db_first_comp_pe(decomp_standard_ocean),          &
     &         decomp_db_last_comp_pe(decomp_standard_ocean)

          IF (decomp_db_g_gridpos(1,iproc,decomp_nowrap_ocean)          &
     &         ==  0) THEN  ! this processor at left of LPG

            DO ihalo=1,Nhalo_max
               decomp_db_g_lasize(1,fld_type_p,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)=         &
     &         decomp_db_g_lasize(1,fld_type_p,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)-1
               decomp_db_g_lasize(1,fld_type_u,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)=         &
     &         decomp_db_g_lasize(1,fld_type_u,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)-1
               decomp_db_g_lasize(1,fld_type_v,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)=         &
     &         decomp_db_g_lasize(1,fld_type_v,ihalo,                   &
     &                              iproc,decomp_nowrap_ocean)-1
            ENDDO


            decomp_db_g_blsize(1,fld_type_p,iproc,decomp_nowrap_ocean)= &
     &      decomp_db_g_blsize(1,fld_type_p,iproc,decomp_nowrap_ocean)-1

            decomp_db_g_blsize(1,fld_type_u,iproc,decomp_nowrap_ocean)= &
     &      decomp_db_g_blsize(1,fld_type_u,iproc,decomp_nowrap_ocean)-1

            decomp_db_g_blsize(1,fld_type_v,iproc,decomp_nowrap_ocean)= &
     &      decomp_db_g_blsize(1,fld_type_v,iproc,decomp_nowrap_ocean)-1


          ELSE  ! processor not at left of LPG

            decomp_db_g_datastart(1,iproc,decomp_nowrap_ocean)=         &
     &        decomp_db_g_datastart(1,iproc,decomp_nowrap_ocean)-1
            decomp_db_g_datastart_f(1,fld_type_p,iproc,                 &
     &                              decomp_nowrap_ocean)=               &
     &       decomp_db_g_datastart_f(1,fld_type_p,iproc,                &
     &                               decomp_nowrap_ocean) - 1
            decomp_db_g_datastart_f(1,fld_type_u,iproc,                 &
     &                              decomp_nowrap_ocean)=               &
     &       decomp_db_g_datastart_f(1,fld_type_u,iproc,                &
     &                               decomp_nowrap_ocean) - 1
            decomp_db_g_datastart_f(1,fld_type_v,iproc,                 &
     &                              decomp_nowrap_ocean)=               &
     &       decomp_db_g_datastart_f(1,fld_type_v,iproc,                &
     &                               decomp_nowrap_ocean) - 1

          ENDIF

          IF (decomp_db_g_gridpos(1,iproc,decomp_nowrap_ocean)          &
     &        ==  (decomp_db_gridsize(1,decomp_nowrap_ocean)-1) )       &
     &    THEN  ! this processor at right of LPG

            DO ihalo=1,Nhalo_max
               decomp_db_g_lasize(1,fld_type_p,ihalo,                   &
     &                           iproc,decomp_nowrap_ocean)=            &
     &         decomp_db_g_lasize(1,fld_type_p,ihalo,                   &
     &                            iproc,decomp_nowrap_ocean)-1

               decomp_db_g_lasize(1,fld_type_u,ihalo,                   &
     &                           iproc,decomp_nowrap_ocean)=            &
     &         decomp_db_g_lasize(1,fld_type_u,ihalo,                   &
     &                            iproc,decomp_nowrap_ocean)-1

               decomp_db_g_lasize(1,fld_type_v,ihalo,                   &
     &                           iproc,decomp_nowrap_ocean)=            &
     &         decomp_db_g_lasize(1,fld_type_v,ihalo,                   &
     &                            iproc,decomp_nowrap_ocean)-1
            ENDDO

            decomp_db_g_blsize(1,fld_type_p,iproc                       &
     &                              ,decomp_nowrap_ocean)=              &
     &                decomp_db_g_blsize(1,fld_type_p,iproc             &
     &                                       ,decomp_nowrap_ocean)-1

            decomp_db_g_blsize(1,fld_type_u,iproc                       &
     &                              ,decomp_nowrap_ocean)=              &
     &                decomp_db_g_blsize(1,fld_type_u,iproc             &
     &                              ,decomp_nowrap_ocean)-1

            decomp_db_g_blsize(1,fld_type_v,iproc                       &
     &                              ,decomp_nowrap_ocean)=              &
     &                decomp_db_g_blsize(1,fld_type_v,iproc             &
     &                              ,decomp_nowrap_ocean)-1

          ENDIF

        ENDDO

      ENDIF

! Finally, indicate this decomposition has been initialised

      decomp_db_set(decomp_nowrap_ocean)=.TRUE.
! And return the new horizontal dimensions


       local_row_len=decomp_db_g_lasize(1,fld_type_p,halo_type_single,  &
     &                     mype,decomp_standard_ocean)

       local_n_rows=decomp_db_g_lasize(2,fld_type_p,halo_type_single,   &
     &                     mype,decomp_standard_ocean)

      RETURN
      END SUBROUTINE DECOMPOSE_OCEAN

