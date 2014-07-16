
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Block Data Subprogram : BLKDATA
!  
!   Purpose : Holds DATA for any variables that are in common blocks,
!             so that they are initialised in only one place.
!  
!   Written by P.Burton
!  

      BLOCK DATA BLKDATA

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
!====================== COMDECK CNTL_IO ========================
! Description:
!
!     Defines the sector size for well-formed transfers on Cray
!     Research systems.  Disk addresses must start on a sector
!     boundary, and transfers must be a number of sectors.  Disk
!     word addresses start at 0.
!
!     On the T3E, well-formed transfers must also start on a
!     cache-line boundary in memory.
!
!   4.3    30/04/97  New deck       B. Carruthers, Cray Research
!   4.4    27/10/97  Remove DATA statement. C.P. Jones
!   5.1    07/06/00  Upon VAR request, provide alternative to the
!                    common statement for um_sector_size.
!                    JC Thil
!
!
      INTEGER UM_SECTOR_SIZE    ! Sector size on disk for I/O
!
      COMMON / CNTL_IO / UM_SECTOR_SIZE
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
!include file: decompdt.h
!========================== COMDECK DECOMPDT ====================
! Description:
!
!     This COMDECK contains data initialisation for variables in
!     DECOMPDB.
!
!   History:
!
!   5.3    05/06/01  New deck.  A. Van der Wal
!

      DATA decomp_db_set / max_decomps * .FALSE. /

! ---------------------- End of comdeck DECOMPDT ----------------------
! ACPARM start
! parameters used for dimensioning small permanent arrays
! values are maximum likely to avoid recompile on resolution change
! the actual dimensions are passed as arguments
! 5.0 17/06/99 Replace amaxsize by amaxsacp for expediency - to retain
!              old, invalid sizes temporarily. R. Rawlins
! ----------------------- include file: AMAXSACM -----------------------
! Description: Quick fix replacement of AMAXSIZE reference in ACPARM.
!              Contains sizes superseded at 5.0, but still embedded
!              throughout AC assimilation code. This allows AC
!              routines to compile without having to re-analyse and
!              correct every routine at this stage. This file should
!              only be of transient use as an interim solution.
!
! Current Code Owner: R. Rawlins
!         5.2   30/11/00 remove ROW_LENGTH_MAX and HORIZ_DIM_MAX
!                        already in amaxsize (now called in ac_ctl)
!                                          B Macpherson
!         6.1   31/08/04 Allow up to 100 levels.  R.Barnes
!         6.2   25/11/05 Set p_rows_max for N320.  R.Barnes

      INTEGER,PARAMETER:: P_ROWS_MAX = 481 ! Max number of rows
      INTEGER,PARAMETER:: P_LEVELS_MAX = 100 ! Max no. of total levels
      INTEGER,PARAMETER:: Q_LEVELS_MAX = 100 ! Max no. of wet levels

! AMAXSACP end
      INTEGER,PARAMETER:: NOBTYPMX = 182

! ACPARM end
!-----------------------------------------------------------------------
!LCOMDECK COMOBS
! Description
!   This comdeck provides parameters relating to observations in
!   atmospheric AC assimilation (see DOCOBS for detail)
!
!   History:
!   Model    Date     Modification history
!  version
!   4.4      3/7/97   Add MAX_NUM_ACOB_FILES,PER_FILE_TNDVMAX S. Bell
!   4.5      5/8/98   Increase USED_FILES size S. Bell
!   5.3    07/06/01   Move data statements to blkdata.  A van der Wal
!L-------------------------------------------------------------------
      INTEGER NDATAVMX, NOBLEVMX
      INTEGER MAX_NUM_ACOB_FILES
      INTEGER NUM_OB_FILE_TYPES
      PARAMETER (NDATAVMX = 6+3*P_LEVELS_MAX)
      PARAMETER (NOBLEVMX = P_LEVELS_MAX+1)
      PARAMETER (MAX_NUM_ACOB_FILES=100)
      PARAMETER (NUM_OB_FILE_TYPES = 10)
      INTEGER NOBTYP,NDVHDR,MAXNLEV1,                                   &
     & OBSTYP(NOBTYPMX),NOBLEV(NOBTYPMX),NDATAV(NOBTYPMX),              &
     & NERLEV1(NOBTYPMX),NOBS(NOBTYPMX),OBLEVTYP(NOBTYPMX),             &
     & MDISPOBT(NOBTYPMX),OBS_NO_ST(NOBTYPMX),                          &
     & OBS_REF_YY, OBS_REF_MM, OBS_REF_DD, OBS_REF_HH, OBS_REF_MIN

      REAL           MISSD,                                             &
     &               OBLEVELS(NOBLEVMX,NOBTYPMX),                       &
     &               OBLAYERB(NOBLEVMX+1,NOBTYPMX),                     &
     &               TIMEINT,TIMENEXT,                                  &
     &               OBS_LAT_N, OBS_LAT_S, OBS_LONG_E, OBS_LONG_W

      INTEGER PER_FILE_TNDVMAX(MAX_NUM_ACOB_FILES)

      CHARACTER*30 OB_FILE_TYPE(NUM_OB_FILE_TYPES)

      CHARACTER*256 USED_FILES(MAX_NUM_ACOB_FILES)
      INTEGER FILENAME_LEN(MAX_NUM_ACOB_FILES)
      INTEGER NUM_USED_FILES

      COMMON /COMOBS/ NOBTYP,NDVHDR,MAXNLEV1,                           &
     & OBSTYP,NOBLEV,NDATAV,NERLEV1,NOBS,OBLEVTYP,                      &
     & MDISPOBT,OBS_NO_ST,MISSD,OBLEVELS,OBLAYERB,                      &
     & OBS_REF_YY, OBS_REF_MM, OBS_REF_DD, OBS_REF_HH, OBS_REF_MIN,     &
     & TIMEINT, TIMENEXT,                                               &
     & OBS_LAT_N, OBS_LAT_S, OBS_LONG_E, OBS_LONG_W, PER_FILE_TNDVMAX,  &
     & FILENAME_LEN,NUM_USED_FILES,                                     &
     ! Character variables at the end
     & OB_FILE_TYPE,USED_FILES
!-----------------------------------------------------------------------
!
!  Description: This comdeck defines the constants for the 3B and 4A
!               versions of the Gravity Wave Drag Code. These are
!               tuneable parameters but are unlikely to be changed.
!
!  History:
!  Version    Date     Comment
!  -------    ----     -------
!    3.4     18/10/94  Original Version    J.R. Mitchell
!    4.3      7/03/97  Remove KAY_LEE (now set in RUNCNST) S.Webster
!    4.5     03/08/98  Add GAMMA_SATN (Used in 06_3B). D. Robinson
!    5.0     14/07/99  Remove redundant switches/variables: only keep
!                      s-version 06_3B parameters. R Rawlins
!    5.2     15/11/00  Set parameters for the 4A scheme.
!                                                  Stuart Webster.
!    5.3     16/10/01  Partition 3B and 4A scheme parameters.
!                                                  Stuart Webster.
!     5.3     21/09/01  Set parameters for the spectral (middle
!                       atmosphere) gravity wave forcing scheme.
!                       Warner and McIntyre, J. Atm. Sci., 2001,
!                       Scaife et al., Geophys. Res. Lett., 2000
!                       Scaife et al, J. Atm. Sci., 2002 give
!                       further details.
!                                                    Adam Scaife
!    5.4     Alters c_gwave.h include file to change from launch level
!            to launch eta to make less model level dependent.
!                                                    Adam Scaife
!    5.5     25/02/03  Remove 3B GWD parameter settings. Stuart Webster
!
!    6.2     16/06/05  Move CCL parameter to gw_ussp. Andrew Bushell
!
!
      ! Number of standard deviations above the mean orography of top
      !  of sub-grid mountains
      REAL,PARAMETER :: NSIGMA = 2.5

      ! Switch to determine which wave saturation test is used
      INTEGER,PARAMETER :: Amplitude_saturation = 1
      INTEGER,PARAMETER :: Stress_saturation    = 0

! SPECTRAL GRAVITY WAVE SCHEME PARAMETERS:

! LMINL = 1/max wavelength at launch
! LSTAR = 1/characteristic wavelength
! ETALAUNCH = eta for model launch

      REAL,PARAMETER:: LMINL = 1./20000
      REAL,PARAMETER:: LSTAR = 1./4300
      REAL,PARAMETER:: ETALAUNCH = 0.045
! AERCMP3A start
!     ------------------------------------------------------------------
!     MODULE TO SET INDICES OF AEROSOL COMPONENTS.
!   4.5   Aug 1998     Set indices for two soot aerosol species.
!                                                 Luke Robinson
!
!
!   5.2   Dec 2000     Set indices for two sea-salt aerosol modes.
!                                                 Andy Jones
!
!   5.4   May 2002     Correct value of NPD_AEROSOL_COMPONENT
!                      and add warning comment.
!                                                 Andy Jones
!   5.5   Feb 2003     Set indices for two modes of biomass
!                      aerosol.                   P Davison
!   5.5   Feb 2003     Set indices for mineral dust aerosol.
!                                        S Woodward
!
      ! maximum number of aerosol components
      ! N.B: this must be at least as large as the
      !      largest value in the list below
      INTEGER,PARAMETER:: NPD_AEROSOL_COMPONENT=28

      INTEGER,PARAMETER:: IP_WATER_SOLUBLE=1
      INTEGER,PARAMETER:: IP_DUST_LIKE=2
      INTEGER,PARAMETER:: IP_OCEANIC=3
      INTEGER,PARAMETER:: IP_SOOT=4
      INTEGER,PARAMETER:: IP_ASH=5
      INTEGER,PARAMETER:: IP_SULPHURIC=6
      INTEGER,PARAMETER:: IP_ACCUM_SULPHATE=10
      INTEGER,PARAMETER:: IP_AITKEN_SULPHATE=11
      INTEGER,PARAMETER:: IP_FRESH_SOOT=12
      INTEGER,PARAMETER:: IP_AGED_SOOT=13
      INTEGER,PARAMETER:: IP_SEASALT_FILM=15
      INTEGER,PARAMETER:: IP_SEASALT_JET=16
      INTEGER,PARAMETER:: IP_DUST_1=17
      INTEGER,PARAMETER:: IP_DUST_2=18
      INTEGER,PARAMETER:: IP_DUST_3=19
      INTEGER,PARAMETER:: IP_DUST_4=20
      INTEGER,PARAMETER:: IP_DUST_5=21
      INTEGER,PARAMETER:: IP_DUST_6=22
      INTEGER,PARAMETER:: IP_BIOMASS_1=23
      INTEGER,PARAMETER:: IP_BIOMASS_2=24
      INTEGER,PARAMETER:: IP_BIOGENIC=25
      INTEGER,PARAMETER:: IP_OCFF_FRESH=26
      INTEGER,PARAMETER:: IP_OCFF_AGED=27
      INTEGER,PARAMETER:: IP_DELTA=28
! AERCMP3A end
! ENTCNST start

!  6.2   28/02/05  Modifications for STPH_RP
      ! coefficients used in calculation of entrainment rate
      REAL,PARAMETER:: AE1     = 1.0
      REAL,PARAMETER:: AE2     = 1.5
      REAL:: ENTCOEF
      COMMON/STPH_RP_1/entcoef
      REAL,PARAMETER:: SH_FAC  = 1.0

! ENTCNST end
! ------------------------------------------------------------------
! Description
!   This comdeck includes the initialization of entcoef for STPH_RP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.2   18/07/06   New include for the
!                    Stochastic Physics Random Parameters scheme
!                    (STPH_RP). A. Arribas
!
! ------------------------------------------------------------------
!

      DATA ENTCOEF /3.0/
! C_LSPDRP start

      ! Microphysics parameters

      ! Drop size distribution for rain: N(D) =  N0 D^m exp(-lambda D)
      ! where N0 = X1R lambda^X2R  and m=X4R
!     REAL, PARAMETER :: X1R is set in the UMUI
!     REAL, PARAMETER :: X2R is set in the UMUI
!     REAL, PARAMETER :: X4R is set in the UMUI

      ! Drop size distribution for graupel: N(D) =  N0 D^m exp(-lambda D)
      ! where N0 = X1G lambda^X2G  and m=X4G
      REAL, PARAMETER :: X1G=5.E25
      REAL, PARAMETER :: X2G=-4.0
      REAL, PARAMETER :: X4G=2.5

      ! Particle size distribution for ice: N(D) = N0 D^m exp(-lambda D)
      ! where N0 = X1I TCG lambda^X2I, m=X4I and TCG=exp(- X3I T[deg C])
!     REAL, PARAMETER :: X1I is set in the UMUI
      REAL, PARAMETER :: X2I=0.0
      REAL, PARAMETER :: X3I=0.1222
      REAL, PARAMETER :: X4I=0.0
!     REAL, PARAMETER :: X1IC is set in the UMUI
      REAL, PARAMETER :: X2IC=0.0
      REAL, PARAMETER :: X3IC=0.1222
      REAL, PARAMETER :: X4IC=0.0

      ! Mass diameter relationship for graupel:  m(D) = AG D^BG
      REAL, PARAMETER :: AG=261.8
      REAL, PARAMETER :: BG=3.0

      ! Mass diameter relationship for ice:  m(D) = AI D^BI
      ! These are set in the UMUI (at 6.6). Values at 6.5 are listed below.
      ! Recommended values for the generic particle size distribution are
      ! from Brown and Francis and are
      ! AI=AIC=1.85E-2, BI=BIC=1.9. If l_calcfall is changed from .true.
      ! then the generic psd values should be set below to ci=cic=8.203
      ! and di=dic=0.2888
      ! REAL, PARAMETER :: AI=0.0444
      ! REAL, PARAMETER :: BI=2.1
      ! REAL, PARAMETER :: AIC=0.587
      ! REAL, PARAMETER :: BIC=2.45

      ! The area diameter relationships are only used if
      ! L_CALCFALL=.TRUE.
      ! Area diameter relationship for ice:  Area(D) = RI D^SI
      REAL, PARAMETER :: RI=0.131
      REAL, PARAMETER :: SI=1.88
      REAL, PARAMETER :: RIC=0.131
      REAL, PARAMETER :: SIC=1.88

      ! The Best/Reynolds relationships are only used if
      ! L_CALCFALL=.TRUE.
      ! Relationship between Best number and Reynolds number:
! Re(D) =LSP_EI(C) Be^LSP_FI(C)
      ! These values are set in the UMUI, but the default values are
      ! listed below. N.B. these have been renamed for VN7.3 with 
      ! 'LSP_' added from the previous versions to avoid conflicts 
      ! later in the code. 
      ! REAL, PARAMETER :: LSP_EI=0.2072  Set in the UMUI
      ! REAL, PARAMETER :: LSP_FI=0.638   Set in the UMUI
      ! REAL, PARAMETER :: LSP_EIC=0.2072 Set in the UMUI
      ! REAL, PARAMETER :: LSP_FIC=0.638  Set in the UMUI


      ! The fall speeds of ice particles are only used if
      ! L_CALCFALL=.FALSE.
      ! Fall speed diameter relationships for ice:
      ! vt(D) = CI D^DI
      REAL, PARAMETER :: CI0=14.3
      REAL, PARAMETER :: DI0=0.416
      REAL, PARAMETER :: CIC0=74.5
      REAL, PARAMETER :: DIC0=0.640

      ! Axial ratio (c-axis divided by a-axis) ESTIMATES. These are not
      ! consistent with those from the area diameter relationships.
      REAL, PARAMETER :: AR=1.0
      REAL, PARAMETER :: ARC=1.0

      ! Fall speed diameter relationship for rain: vt(D) = CR D^DR
      REAL, PARAMETER :: CR=386.8
      REAL, PARAMETER :: DR=0.67

      ! Fall speed diameter relationship for graupel: vt(D) = CG D^DG
      REAL, PARAMETER :: CG=253.0
      REAL, PARAMETER :: DG=0.734

      ! Do we wish to calculate the ice fall velocities?
      ! TRUE if calculate speeds, FALSE if specify speeds
      LOGICAL, PARAMETER :: L_CALCFALL=.TRUE.

!*L --------------------- Comdeck: CENVIRDT ---------------------------
!
!   Purpose: Data statements for character enviroment variables used
!            by portable IO to open/close files (links with CENVIR)
!
!   Author : R A Stratton      Date : 22/10/92
!
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  3.4   1/8/94     Revised Obs file spec + CX files: Stuart Bell
!    3.4   18/5/94    Correct misspelling of CURNTIN and change length
!                     to match. J F Thomson
!    4.0  22/9/95     Units for spectral data added. J. M. Edwards
!LL  4.1  11/03/96  Introduce Wave sub-model.  RTHBarnes.
!    4.1 26/02/96     New env. variables for sulphur ancillary files.
!                     Rename SOURCES to SULPEMIS. D. Robinson.
!      4.3  18/3/97   Add aerosol fcgs of climate change. William Ingram
!    4.4  4/7/97    Add ANLINCR at 108. Chris Jones/Stuart Bell
!LL  4.4 12/9/97      New ancillary file environment variables for
!LL                   initial surface type fracs, initial vegetation
!LL                   state and vegetation disturbance.
!LL                                                  R.A.Betts
!    4.4 28/08/97     Move CACHED from logical unit no 3 to 138 in
!                     order to release a Fortran unit [for OASIS].
!                     R.Rawlins
!    4.5 22/04/98     Add new ancillary file for CO2 emissions:
!                     CO2EMITS - in I/O unit 118. Chris Jones.
!    4.5 22/04/98     Add new ancillary file for soot emissions:
!                     SOOTEMIS - in I/O unit 139. R.Rawlins
!    4.5 29/07/98     Move ALABCOU1/2/3/4 from 101-103 to 140-143.
!                     Add ALABCOU5/6/7/8 to 144-147. D. Robinson.
!    4.5 17/08/98     Remove OLABCOUT from Unit 90. Add OLABCOU1/2/3/4
!                     to 101-103. D. Robinson.
!    5.1 13/04/00     TDF_dump replaces FILE107.
!                     IAU_inc  replaces ANLINCR. Adam Clayton
!    5.2 21/08/00     Add an extra op macro for VAR plus 1 user pp outpu
!                     stream. R Rawlins
!    5.2 18/01/01     Add VERT_LEV. D. Robinson
!    5.3 26/10/01     Add LANDFRAC. Remove OBS06-OBS10, CX01-CX10.
!                     D Robinson
!    5.3 24/10/01     Add IDEALISE.                 Andy Malcolm
!    5.3 19/06/2001   Added file (unit no 119) for TPPSOZON. Dave Tan
!    5.4 29/08/02     Add ALABCIN1/2 to unit no 125/126. Remove ALABCIN
!                     from 95. D Robinson
!    5.4 29/08/02     A
!    5.5 17/02/03     Add Wave model interface & boundary files.
!                                                    D.Holmes-Bell
!    5.5 30/12/02     Add DUSTSOIL/BIOMASS to unit no 75/76.
!                     and RIVSTOR/RIVSEQ/RIVDIR to 77-79. D Robinson
!    6.1 07/04/04     Add DMSCONC to unit no 95.   A. Jones
!    6.1 08/11/04     Alter names of River Routing files. R.Sharp
!    6.2 26/01/06     Add ICFILE to unit no 152. T. Edwards
!    6.2 17/03/06     Add SURFEMIS,AIRCREMS,STRATEMS,EXTRAEMS,RADONEMS
!                     to unit nos 130-134. D.Robinson
!    6.2 24/00/05     Change STRATOUT to PPSMC and MESOUT to PPSCREEN,
!                     files for soil moisture nudging scheme. Clive Jones
!LL
!LL DATA statements for COMDECK CENVIR

      DATA FT_ENVIRON/                                                  &
     &  'PPXREF  ','PPXREFU ','        ','STASHCTL','        ',         &
                                                                !  1- 5
     &  '        ','OUTPUT2 ','        ','        ','XHIST   ',         &
     &  'IHIST   ','THIST   ','        ','ERRFLAG ','CACHE1  ',         &
     &  'CACHE2  ','AOTRANS ','ASWAP   ','OSWAP   ','AINITIAL',         &
     &  'ASTART  ','        ','APSUM1  ','APSTMP1 ','        ',         &
     &  '        ','AOMEAN  ','ATMANL  ','        ','OZONE   ',         &
     &  'SMCSNOWD','DSOILTMP','SOILTYPE','VEGTYPE ','SSTIN   ',         &
     &  'SICEIN  ','PERTURB ','CURNTIN ','        ','OINITIAL',         &
     &  'OSTART  ','        ','OPSUM1  ','OPSTMP1 ','        ',         &
     &  '        ','OCNANL  ','ATRACER ','OTRACER ','WFIN    ',         &
     &  'HFLUXIN ','PMEIN   ','ICEFIN  ','AIRTMP  ','SALINITY',         &
     &  'FLUXCORR','SWSPECTD','BAS_IND ','SLABHCON','PP0     ',         &
     &  'PP1     ','PP2     ','PP3     ','PP4     ','PP5     ',         &
     &  'PP6     ','PP7     ','PP8     ','PP9     ','OBS01   ',         &
                                                               !66-70
     &  'OBS02   ','OBS03   ','OBS04   ','OBS05   ','DUSTSOIL',         &
                                                               !71-75
     &  'BIOMASS ','RIVSTOR ','RIVCHAN ','RIVER2A ','LWSPECTD',         &
                                                               !76-80
     &  'WAVEOUT ','SURGEOUT','PPSCREEN','PPSMC   ','WFOUT   ',         &
     &  'HFLUXOUT','PMEOUT  ','ICFOUT  ','MOSOUT  ','VERT_LEV',         &
     &  'SSTOUT  ','SICEOUT ','CURNOUT ','        ','DMSCONC ',         &
                                                               !91-95
     &  'OROG    ','TRANSP  ','OLABCIN ','OCNDEPTH',                    &
     &  'OLABCOU1','OLABCOU2','OLABCOU3','OLABCOU4','FILE104 ',         &
     &  'FILE105 ','IDEALISE','TDF_dump','IAU_inc ','MURKFILE',         &
     &  'SULPEMIS','USRANCIL','USRMULTI','OUSRANCL','OUSRMULT',         &
     &  'SO2NATEM','CHEMOXID','AEROFCG ','CO2EMITS','TPPSOZON',         &
     &  'LANDFRAC','WLABCOU1','WLABCOU2','WLABCOU3','WLABCOU4',         &
                                                               !120-124
     &  'ALABCIN1','ALABCIN2','        ','OCFFEMIS','HORZGRID',         &
                                                               !125-129
     &  'SURFEMIS','AIRCREMS','STRATEMS','EXTRAEMS','RADONEMS',         &
                                                               !130-134
     &  'FRACINIT','VEGINIT ','DISTURB ','CACHED  ','SOOTEMIS',         &
                                                               !135-139
     &  'ALABCOU1','ALABCOU2','ALABCOU3','ALABCOU4','ALABCOU5',         &
                                                               !140-144
     &  'ALABCOU6','ALABCOU7','ALABCOU8','CARIOLO3','RPSEED  ',         & 
                                                               !145-149
     &  'PPVAR   ','PP10    ','ICFILE  ','VAR_GRID','ARCLBIOG',         &
                                                               !150-154
     &  'ARCLBIOM','ARCLBLCK','ARCLSSLT','ARCLSULP','ARCLDUST',         &
                                                               !155-159
     &  'ARCLOCFF','ARCLDLTA','        ','        ','        ',         &
                                                               !160-164
     &          35*'        '                                           &
     & /
!
      DATA LEN_FT_ENVIR/6,7,0,8,0, 0,7,0,0,5,                           &
                                                    !  1-10
     &                  5,5,0,7,6, 6,7,5,5,8,                           &
                                                    ! 11-20
     &                  6,0,6,7,0, 0,6,6,0,5,                           &
                                                    ! 21-30
     &                  8,8,8,7,5, 6,7,7,0,8,                           &
                                                    ! 31-40
     &                  6,0,6,7,0, 0,6,7,7,4,                           &
                                                    ! 41-50
     &                  7,5,6,6,8, 8,8,7,8,3,                           &
                                                    ! 51-60
     &                  3,3,3,3,3, 3,3,3,3,5,                           &
                                                    ! 61-70
     &                  5,5,5,5,8, 7,7,7,7,8,                           &
                                                    ! 71-80
     &                  7,8,8,5,5, 8,6,6,6,8,                           &
                                                    ! 81-90
     &                  6,7,7,0,7, 4,6,7,8,                             &
                                                    ! 91-99
     &                  8,8,8,8,7, 8,8,8,7,8,                           &
                                                    ! 100-109
     &                  8,8,8,8,8, 8,8,7,8,8,                           &
                                                    ! 110-119
     &                  8,8,8,8,8, 8,8,0,8,8,                           &
                                                    ! 120-129
     &                  8,8,8,8,8, 8,7,7,6,8,                           &
                                                    ! 130-139
     &                  8,8,8,8,8, 8,8,8,8,6,                           &
                                                    ! 140-149
     &                  5,4,6,8,8, 8,8,8,8,8,                           &
                                                    ! 150-159
     &                  8,8,0,0,0, 0,0,0,0,0,                           &
                                                    ! 160-169
     &                  30*0/                       ! 170-199

!End of COMDECK CENVIRDT

!
!====================== COMDECK CNTLIODT ========================
! Description:
!
!     Defines the sector size for well-formed transfers on Cray
!     Research systems.  Disk addresses must start on a sector
!     boundary, and transfers must be a number of sectors.  Disk
!     word addresses start at 0.
!
!     On the T3E, well-formed transfers must also start on a
!     cache-line boundary in memory.
!
!   4.4    24/10/97  New deck       C.P. Jones
!   4.5    02/10/98  Increase from 512 to 2048. D. Robinson.
!
!
      DATA UM_SECTOR_SIZE/2048/
!include file: comobsdt.h
!========================== COMDECK COMOBSDT ====================
! Description:
!
!     This COMDECK contains data initialisation for variables in
!     COMOBS.
!
!   History:
!
!   5.3    05/06/01  New deck.  A. Van der Wal
!

      DATA OB_FILE_TYPE/                                                &
     & 'Surface                       ',                                &
     & 'Sonde                         ',                                &
     & 'Aircraft                      ',                                &
     & 'Sat120                        ',                                &
     & 'Sat500                        ',                                &
     & 'GLOSS                         ',                                &
     & 'Satwind                       ',                                &
     & 'Scatwind                      ',                                &
     & 'MOPS                          ',                                &
     & 'Test                          '/
! ---------------------- End of comdeck COMOBSDT ----------------------
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


! CMAXSIZE maximum sizes for dimensioning arrays
! of model constants whose sizes are configuration dependent. This
! allows constants to be read in from a NAMELIST file and maintain
! the flexibility of dynamic allocation for primary variables. The
! maximum sizes should agree with the maximum sizes implicit in the
! front-end User Interface.

!
!  Model            Modification history:
! version  Date
! 3.2  26/03/93  New COMDECK. Author R.Rawlins
! 3.4  06/08/94: Parameter MAX_NO_OF_SEGS used to dimension addresses
!                in macro-tasked calls to SWRAD, LWRAD & CONVECT.
!                Authors: A.Dickinson, D.Salmond, Reviewer: R.Barnes
! 3.5  22/05/95  Add MAX_N_INTF. D. Robinson
! 4.5  29/07/98  Increase MAX_N_INTF/MAX_N_INTF_A to 8. D. Robinson.
! 5.0  20/04/99  Changes for conversion to C-P C dynamics grid.
!                R. Rawlins
!  6.1   04/08/04  Add diffusion variable max_power     Terry Davies
! 6.2  25/12/05  Add max_updiff_levels/max_sponge_width   Terry Davies

      INTEGER,PARAMETER::max_model_levels = 100 ! Maximum no. of levels

      ! Max levels in boundary layer
      INTEGER,PARAMETER:: max_bl_levels = max_model_levels

      ! Max size of alpha_Cd
      INTEGER,PARAMETER :: max_number_alpha_cds = max_bl_levels

      ! Max no. of levels for pvort output
      INTEGER,PARAMETER :: MAX_REQ_THPV_LEVS = max_model_levels

      ! Max no. 1-2-1 rows in polar filter
      INTEGER,PARAMETER ::  max_121_rows =  8
      ! 0 is used for horizontal diffusion pointer

      ! Max no. of levels (from top) to apply upper level diffusion
      INTEGER,PARAMETER ::  max_updiff_levels = 10

      ! Max size of any sponge zones
      INTEGER,PARAMETER ::  max_sponge_width = 10

      ! Max size of look-up tables for searches
      INTEGER,PARAMETER ::  max_look = 2048

      ! Max no. of atmos interface areas
      INTEGER,PARAMETER :: MAX_N_INTF_A =  8

      ! Max no. of points in LBC      
      INTEGER,PARAMETER :: max_intf_lbcrow_length = 1000
      INTEGER,PARAMETER :: max_intf_lbcrows = 1000
        
      ! Max no. of atmos interface levels
      INTEGER,PARAMETER :: MAX_INTF_LEVELS = max_model_levels

      ! Maximum number of physics segments
      INTEGER,PARAMETER :: MAX_NO_OF_SEGS = 200
      ! MAX_N_INTF/MAX_N_INTF_A to be sorted out in later version
      ! Max no. of interface areas
      INTEGER, PARAMETER :: MAX_N_INTF =  8
! CMAXSIZE end
! Data Statements for PARCOMM COMMON block

! History:
!   5.1      22/05/00  New include file created.          P.Burton
!   5.5      07/02/03  Initialise g_at_extremity for SX     E.Leung

      DATA current_decomp_type/-1/  ! set the initial decomposition
!                                   ! to an "unset" value

      END BLOCKDATA BLKDATA
