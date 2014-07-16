
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE PP_FILE -----------------------------------------
!LL
!LL  Purpose:- To output a field to a PP_FILE
!LL
!LL  Tested under compiler CFT77
!LL  Tested under OS version 5.1
!LL
!LL TJ, RR      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.2  19/04/93  Code for new real missing data indicator (TCJ).
!LL  3.4   04/08/94  No packing indicator change from -26 to -99  PJS
!LL  4.1   22/11/96  Modify I/O calls for mpp use  P.Burton
!LL  4.3   30/04/97  Added code to use UM_SECTOR_SIZE to make transfers
!LL                  well-formed.
!LL                  B. Carruthers  Cray Research.
!LL  4.4   16/06/97  Add processing after the write, so
!LL                  that all the processors know the answer
!LL                    Author: Bob Carruthers, Cray Rsearch.
!LL  4.5   28/05/98  Code for parallel processing in COEX Packing
!LL                    Author: Paul Burton & Bob Carruthers
!LL  5.2   11/02/00 Code for Run Length Encoding  I. Edmond
!LL
!LL   5.1 May. 00  Coex now returns an error code set greater than 0
!LL                if invalid data is being stored. (i.e if the number
!LL                of bits per value used in the wgdos packed data is
!LL                >31). PPFILE checks the return code from coex and
!LL                throws an exception by returning the error code from
!LL                coex (summed over all pes). Ian Edmond
!    5.4    11/04/02 Set up extra data vectors for variable
!                    horizontal grids (only ocean catered for
!                    at this stage). Also prevent RLE when
!                    extra data vectors are present.  R. Hill
!    5.5    14/01/03 Ensure WFIO output buffer is extended
!                    to acommodate extra data. R. Hill
!    5.5    20/02/03 Allowed for inclusion of PARVARS comdeck
!                    on non-T3E mpp platforms to provide variables
!                    needed by COMVGRID comdeck.      P.Dando
!    6.0    15/01/03  Replace calls to SETPOS and BUFFOUT with
!                     buffered versions for PP files.
!                     JC Rioual, P. Selwood
!    6.1    13/07/04 Add packing in 2 stage gather.
!                    P.Selwood/B.Carruthers (Cray)
!    6.1    28/08/04 Buffered IO in C. JC Rioual P Selwood
!   6.2    19/01/05  Allow RLE to be active even if EXTRA data
!                    is present in a record (but not other
!                    forms of packing which could still result
!                    in corruption of the extra data) R. Hill

!LL  Programming standard: U M DOC  Paper NO. 4,
!LL
!LL  Logical components covered C4
!LL
!LL  Project TASK: C4
!LL
!LL  External documentation  C4
!LL
!LLEND-------------------------------------------------------------

!
!*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE PP_FILE(PPFIELD,LENBUF,NUM_WORDS,RMDI,COMP_ACCRCY,     &
     &PPHORIZ_OUT,UNITPP,IWA,N_COLS_OUT,N_ROWS_OUT,PACKING,             &
     & im_ident,LRLE,                                                   &
     &     PACKING_TYPE,current_io_pe,LEN_EXTRA                         &
     &    ,SROW_OUT,WCOL_OUT,ICODE,CMESSAGE)
      IMPLICIT NONE


      INTEGER LEN_EXTRA ! IN size of expected extra data
      INTEGER SROW_OUT  ! 1st southern row to output
      INTEGER WCOL_OUT  ! 1st western column to output
      CHARACTER*(80) CMESSAGE !OUT OUT MESSAGE FROM ROUTINE
!
      LOGICAL                                                           &
     &  PACKING                                                         &
                           !IN OVERALL Packing switch (T if pckng reqd)
     &,  LRLE

      INTEGER                                                           &
     &  ICODE                                                           &
                           !IN    RETURN CODE FROM ROUTINE
     &, LENBUF                                                          &
                           !IN     LENGTH OFF PP BUFFER
     &, UNITPP                                                          &
                           !IN     OUTPUT PP UNIT NUMBER
     &, LEN_IO                                                          &
                           !NOT USED, BUT NEEDED FOR BUFFOUT CALL
     &, im_ident                                                        &
                           !IN    Internal model identifier
     &, current_io_pe                                                   &
                           !IN     PE which will do the I/O
     &, ocode                                                           &
                        !Copy of the input value of ICODE
     &, ncode                                                           &
                        ! seperate item for each processor
     &, istat           ! Isum return status, GC_OK on success.
!
      INTEGER                                                           &
     &  N_ROWS_OUT                                                      &
                      !IN   PPHORIZ_OUT=N_ROWS_OUT*N_COLS_OUT
     &, N_COLS_OUT                                                      &
                      !IN    PPHORIZ_OUT=N_COLS_OUT*N_ROWS_OUT
     &, NUM_OUT                                                         &
                      !IN    NUMBER OF COMPRESSED (32 BIT) WORDS
     &, COMP_ACCRCY                                                     &
                      !IN    PACKING ACCURACY IN POWER OF 2
     &, U_ROWS                                                          &
                      !IN    NO OF U,V, ROWS
     &, P_ROWS                                                          &
                      !IN    PRESS/TEMP ROWS
     &, PPHORIZ_OUT                                                     &
                      !IN    SIZE OF OUTPUT FIELD
     &, NUM_WORDS                                                       &
                      !IN    NUMBER OF 64 BIT WORDS WORTH OF DATA
     &, PACKING_TYPE  ! OUT set to 1 if WGDOS packing else set to zero.
!
      REAL                                                              &
     &  PPFIELD(PPHORIZ_OUT)                                            &
                               !INOUT ARRAY TO STORE PPDATA
     &, BUFOUT(LENBUF+LEN_EXTRA)                                        &
                                  !OUTPUT PP BUFFER (ROUNDED UP)
     &, RMDI                   !IN     Missing data indicator
!
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

! Description: This include file contains information needed when
!              generating variable horizontal grid data in the
!              STASH extra data vector. Introduced UM 5.4 - R. Hill
!===================================================================
      LOGICAL :: X_VAR_GRID ! Whether variable grid in E-W direction
      LOGICAL :: Y_VAR_GRID ! and/or in S-N direction

      INTEGER :: VAR_GRID_TYPE ! 0 = none
                               ! 1 = T grid
                               ! 2 = U/V grid

      ! Grid boundaries for T and U,V
      REAL :: X_BOUNDARY(ROW_LENGTH_MAX+1,2)
      REAL :: Y_BOUNDARY(ROWS_MAX+1,2)

      ! Grid Points for T and U,V
      REAL :: X_GRID(ROW_LENGTH_MAX,2)
      REAL :: Y_GRID(ROWS_MAX,2)

      COMMON /OVARGRID/ X_VAR_GRID,Y_VAR_GRID                           &
     & ,X_BOUNDARY,Y_BOUNDARY,X_GRID,Y_GRID,VAR_GRID_TYPE

      ! The following parameters correspond to the extra data
      ! vector descriptors expected, for e.g., in PV-WAVE
      ! plotting routines (e.g. decode_extra.pro). There are
      ! numerous other areas of code where these integer
      ! descriptors must be handled (e.g. FIELDCOS, PPI2H, FTT)
      ! So it is not a trivial matter to introduce new code
      ! descriptors. Furthermore, ieee -32 will destroy these
      ! integers so PP data must always be processed via the
      ! long winded route: QXFIELDCOS -> PUTIBM -> FTT/PPI2H.
      ! (This thoroughly unsatisfactory state of affairs may
      ! be correctable with developments to ieee and convpp).
      INTEGER,PARAMETER :: x_coord_vector=1
                                     ! Indicates that an extra
                                     ! data vector gives LBNPT
                                     ! x-coordinate values
      INTEGER,PARAMETER :: y_coord_vector=2
                                     ! Indicates that an extra
                                     ! data vector gives LBROW
                                     ! Y-coordinate values
      INTEGER,PARAMETER :: x_lbnd_vector=12
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! x-boundary values
      INTEGER,PARAMETER :: x_ubnd_vector=13
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! x-boundary values
      INTEGER,PARAMETER :: y_lbnd_vector=14
                                     ! Indicates that an extra
                                     ! data vector gives lower
                                     ! y-boundary values
      INTEGER,PARAMETER :: y_ubnd_vector=15
                                     ! Indicates that an extra
                                     ! data vector gives upper
                                     ! y-boundary values

!
!dir$ cache_align bufout
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
! CHSUNITS define the number of i/o units
!
!  Author : R A Stratton
!
!  Model            Modification history:
! version  date
!   3.1  03/02/93   Introduced at version 3.1
!   4.1  21/02/96   Increase no.of i/o units to accommodate wave
!                   sub-model.  RTHBarnes.
!   5.2  21/08/00   Add an extra op macro for VAR plus 1 user pp
!                   output stream. R Rawlins
!   6.2  19/01/06   Increased NUNITS to 152 to accomodate extra
!                   diagnostic
!
! Project task:
!
!  Documentation:  Unified Model Documentation Paper
!                  H- History Bricks
!
! ---------------------------------------------------------------

      ! These values must be consistent with OUTFILE_S, OUTFILE_L
      ! and OUTFILE_E in file VERSION.
      INTEGER,PARAMETER::NUNITS=161   ! No. of I/O units
      ! length of most unit no arrays
      INTEGER,PARAMETER::NUNITS_LEN=NUNITS-19

      ! The above parameter statements must not be altered without
      ! considering the effect on the following HISTORY files CHISTO,
      ! CLFHIST and IHISTO.
      ! This file must always preceed the above history file
      ! New file environment variable names may need to be added to
      ! CLFHIST and/or CENVIRDT (usually both) depending on manner of
      ! I/O.
! CHSUNITS end
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
!*---------------------------------------------------------------------

!*L  WORKSPACE USAGE:-------------------------------------------------
!   DEFINE LOCAL WORKSPACE ARRAYS: 1 REAL ARRAY
!   AT FULL FIELD LENGTH
!
!*---------------------------------------------------------------------
!
!*L EXTERNAL SUBROUTINES CALLED---------------------------------------
      EXTERNAL SETPOS,COEX,BUFFOUT
!*------------------------------------------------------------------
!L  MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS-1) * ROWLENGTH
!L---------------------------------------------------------------------
!----------------------------------------------------------------------
!    DEFINE LOCAL VARIABLES
      INTEGER                                                           &
     &  ML                                                              &
                      !     LONGITUDE COUNTER
     &, JL                                                              &
                      !     LATITUDE COUNTER
     &, IWA                                                             &
                      !     RECORD NUMBER
     &, II                                                              &
                      !     COUNTER
     &, LENGTH_FULLWRD                                                  &
                      !     LENGTH IN BITS OF FULLWORD VAR
     &, LEN_BUF_WORDS !     NUM_WORDS ROUNDED BY 512 AND ACTUALLY

      INTEGER                                                           &
     &  JJ            !     Local counter

      REAL                                                              &
     &  IX            !     RETURN VALUE FROM UNIT COMMAND

      LOGICAL                                                           &
     &  UV                 !

!    Check if an error has already been encountered, and get out
!    if it has.
      ocode = 0
      IF (icode  >   0) then
         goto 999
      ELSE IF (icode  <   0)then
         ocode = icode
         icode = 0
      END IF

!
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS IPPLOOK GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
!======================================================================
      LENGTH_FULLWRD=64   !   LENGTH IN BITS OF FULLWORD VAR
!L   At this point packing,if required,will be done using the WGDOS
!L   method of packing.
      PACKING_TYPE=0
! Note the value of -26 corresponds to -15 (F) in ppxref.
! The packing acuracy is scaled to allow greater accuracy.
! Packing will only be attempted if there are at least 2 points per row
! in the PPfield.
!
!      Run length encoding applies only to unpacked ocean fields, but
!      most ocean fields remain unpacked even when packing profiles
!      5 or 6 are set. Hence selecting, for example, both packing
!      profile 5 and run length encoding makes sense.
      IF(PACKING .AND. N_COLS_OUT  >=  2) THEN
        ! 1. Climate wgdos packing has been selected for current file
        !    stream via UMUI
        IF(COMP_ACCRCY  <=  -99 .AND. LRLE .AND.                        &
     &     im_ident  ==  ocean_im )THEN
           ! 2. STASH packing profile for the field is -99
           ! 3. Run Length Encoding has been selected
           ! 4. Submodel is Ocean
           PACKING_TYPE = 4
        ELSE IF(COMP_ACCRCY  >   -99) THEN
           ! 2. STASH packing profile for the field is set.
           PACKING_TYPE = 1
        ENDIF
      ELSE
         ! 1. Packing may or may not have been selected for current
         !    file stream. This section of code is not executed when
         !    GRIB packing selected.
         IF (LRLE .AND.                                                 &
     &       im_ident  ==  ocean_im) THEN
           ! 2. Run Length Encoding has been selected
           ! 3. Submodel is Ocean
           ! This should be safe at least as far as output
           ! generation goes - the issue will be post
           ! processing functionality - does that have
           ! all the necessary intelligence to unpack
           ! RLE fields AND process extra data
           PACKING_TYPE = 4
         ENDIF
      END IF
!
      IF(PACKING_TYPE == 1)THEN
! DEPENDS ON: coex
        CALL COEX(PPFIELD,PPHORIZ_OUT,BUFOUT,LENBUF,N_COLS_OUT,         &
     &  N_ROWS_OUT,NUM_OUT,COMP_ACCRCY,.TRUE.,RMDI,LENGTH_FULLWRD,      &
     &  icode,cmessage)

        NUM_WORDS=(NUM_OUT+1)/2 ! Round up to the nearest 64 Bit CRAY Wd
!  COEX returns the number of IBM words needed to hold the packed data
!                             ~~~
        LEN_BUF_WORDS=((NUM_WORDS+um_sector_size-1)/um_sector_size)*    &
     &   um_sector_size
      ELSE IF(PACKING_TYPE == 4)THEN
        if (mype  ==  current_io_pe) then
! DEPENDS ON: runlen_encode
          CALL RUNLEN_ENCODE(PPFIELD,PPHORIZ_OUT,BUFOUT,PPHORIZ_OUT,    &
     &                     NUM_OUT,RMDI,ICODE,CMESSAGE)
        ! Size of run length encoded data is greater than unpacked
        ! field therefore leave field unpacked.
          if (NUM_OUT  >=  PPHORIZ_OUT) then
            PACKING_TYPE = 0
            DO JJ=1,PPHORIZ_OUT
              BUFOUT(JJ) = PPFIELD(JJ)
            END DO
            NUM_WORDS=PPHORIZ_OUT
            LEN_BUF_WORDS=LENBUF
          else
            NUM_WORDS=NUM_OUT
            LEN_BUF_WORDS=((NUM_WORDS+um_sector_size-1)/um_sector_size)*&
     &      um_sector_size
          endif

        end if
      ELSE  ! No packing required.
        DO 1 JJ=1,PPHORIZ_OUT
        BUFOUT(JJ) = PPFIELD(JJ)
    1   CONTINUE
        NUM_WORDS=PPHORIZ_OUT
        LEN_BUF_WORDS=LENBUF
      ENDIF

      IF (VAR_GRID_TYPE >  0) THEN

         ! If we have a variable horizontal grid then we
         ! must add the grid data to the end of the PP record
         ! as "EXTRA DATA".
! DEPENDS ON: extra_variable_grid
         CALL EXTRA_VARIABLE_GRID(BUFOUT(NUM_WORDS+1),LEN_EXTRA         &
     &            ,SROW_OUT,N_ROWS_OUT                                  &
     &            ,WCOL_OUT,N_COLS_OUT)

         ! Adjust output buffer size for WFIO.
         NUM_WORDS = NUM_WORDS + LEN_EXTRA

         IF ((PACKING_TYPE == 1).OR.(PACKING_TYPE == 4)) THEN
           LEN_BUF_WORDS=((NUM_WORDS+um_sector_size-1)/um_sector_size)* &
     &                     um_sector_size ! No of words output
          ELSE
             LEN_BUF_WORDS=LENBUF+LEN_EXTRA
         ENDIF
      ENDIF ! If variable grid data needed adding

      if (mype  ==  current_io_pe) then
      DO JJ=NUM_WORDS+1,LEN_BUF_WORDS
        BUFOUT(JJ)= 0.0
      ENDDO
      CALL SETPOS_single(UNITPP,IWA,ICODE)
      CALL BUFFOUT_single(UNITPP,BUFOUT(1),LEN_BUF_WORDS,LEN_IO,IX)

!     WRITE(6,102) IWA,LEN_BUF_WORDS
  100 FORMAT(//,32X,'   ARRAY BUFOUT AT END OF PPOUT ',//,32(10F8.0/))
  102 FORMAT(' FROM PP_FILE    IWA  LEN_BUF_WORDS ',2I12)
!
      IF (IX /= -1.0.OR.LEN_IO /= LEN_BUF_WORDS) THEN
! DEPENDS ON: ioerror
        CALL IOERROR('Buffer out Data Field',IX,LEN_IO,                 &
     &                LEN_BUF_WORDS)
        CMESSAGE='PPFILE  : I/O error - PP Data Field Output'
        ICODE=7
        RETURN
      ENDIF
      endif ! (mype  ==  current_io_pe)
!     If we have found an error, leave it in icode.  If no error
!     occurred then check if the original input value of icode was
!     non-zero (a previous untrapped error/warning), and copy this
!     back into ICODE before eaving the routine.
      IF (icode  ==  0 .and. ocode  /=  0) then
         icode = ocode
      END IF
  999 CONTINUE
      RETURN
      END SUBROUTINE PP_FILE
