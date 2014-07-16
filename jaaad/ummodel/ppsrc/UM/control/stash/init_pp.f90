
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: INIT_PP  -------------------------------------------------
!LL
!LL  Purpose: Initialises direct access PP files at the start of
!LL           the run.  NB: Sequential PP files need no initialisation.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: D401
!LL
!LL  Project task:
!LL
!LL  External documentation: On-line UM document C61 - Zonal mean
!LL                          calculations.
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE INIT_PP ( FTN_UNIT,FILE_TYPE_LETTER,                   &
     &                     LEN1_LOOKUP,PP_LEN2_LOOKUP,FIXHD,            &
     &                     INTHD,REALHD,LEVDEPC,ROWDEPC,COLDEPC,        &
     &                     LEN_FIXHD,LEN_INTHD,                         &
     &                     LEN_REALHD,LEN1_LEVDEPC,LEN2_LEVDEPC,        &
     &                     LEN1_ROWDEPC,LEN2_ROWDEPC,                   &
     &                     LEN1_COLDEPC,LEN2_COLDEPC,                   &
     &                     PP_LEN_INTHD, PP_LEN_REALHD,                 &
     &                     ICODE,CMESSAGE)


      USE FIELD_BUFF_MOD, ONLY :                                        &
     &    INIT_FXH,                                                     &
     &    ATTACH_FXH,                                                   &
     &    INIT_IPPLOOK,                                                 &
     &    ATTACH_IPPLOOK


      IMPLICIT NONE
!
      CHARACTER*1                                                       &
     &    FILE_TYPE_LETTER    ! IN  - File type (p-PP, b-bndry)
      INTEGER                                                           &
     &    FTN_UNIT                                                      &
                              ! IN  - Fortran unit number
     &,   LEN1_LOOKUP                                                   &
                              ! IN  - Size of PP header
     &,   PP_LEN2_LOOKUP                                                &
                              ! IN  - Max allowable fields
     &,   LEN_FIXHD                                                     &
                              ! IN    LENGTH OF FIXED CONSTANTS
     &,   LEN_INTHD                                                     &
                              ! IN    LENGTH OF INTEGER CONSTANTS
     &,   LEN_REALHD                                                    &
                              ! IN    LENGTH OF REAL CONSTANTS
     &,   LEN1_LEVDEPC                                                  &
                              ! IN    LENGTH OF 1st Dim of lev depndt
     &,   LEN2_LEVDEPC                                                  &
                              ! IN    LENGTH OF 2nd Dim of lev depndt
     &,   LEN1_ROWDEPC                                                  &
                              ! IN    LENGTH OF 1st Dim of row depndt
     &,   LEN2_ROWDEPC                                                  &
                              ! IN    LENGTH OF 2nd Dim of row depndt 
     &,   LEN1_COLDEPC                                                  &
                              ! IN    LENGTH OF 1st Dim of col depndt
     &,   LEN2_COLDEPC                                                  &
                              ! IN    LENGTH OF 2nd Dim of col depndt
     &,   ICODE                                                         &
                              ! OUT - Error exit code
     &,   PP_LEN_INTHD                                                  &
                              ! IN - Length of PP FILE integer header
     &,   PP_LEN_REALHD       ! IN - Length of PP FILE real header
!
!
      INTEGER                                                           &
     &    FIXHD(LEN_FIXHD)                                              &
                                    ! IN    ARRAY OF FIXED CONSTANTS
     &,   INTHD(LEN_INTHD)                                              &
                                    ! IN    ARRAY OF integer CONSTANTS
     &,   LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC)                            &
                                              ! IN LEV DEP CONSTANTS
     &,   ROWDEPC(LEN1_ROWDEPC*LEN2_ROWDEPC)                            &
                                              ! IN ROW DEP CONSTANTS
     &,   COLDEPC(LEN1_COLDEPC*LEN2_COLDEPC)                            &
                                              ! IN COL DEP CONSTANT                                                                                  
     &,   PP_INTHD(PP_LEN_INTHD)                                        &
                                    ! OUT   ARRAY of integer constants
     &,   PP_LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC)                         &
                                                 ! OUT Level dep cts
     &,   PP_ROWDEPC(LEN1_ROWDEPC*LEN2_ROWDEPC)                         &
                                                 ! OUT Row dep cts
     &,   PP_COLDEPC(LEN1_COLDEPC*LEN2_COLDEPC)  ! OUT Col dep cts


      INTEGER, POINTER :: PP_FIXHD(:)





!
      REAL                                                              &
     &    REALHD(LEN_REALHD)                                            &
                                    ! IN    ARRAY OF REAL CONSTANTS
     &,   PP_REALHD(PP_LEN_REALHD)  ! OUT   ARRAY OF REAL CONSTANTS
!
      CHARACTER*80                                                      &
     &    CMESSAGE            ! OUT - Error message
!
!*----------------------------------------------------------------------
!
!  External subroutines
!
      EXTERNAL SETPOS,IOERROR,POSERROR,BUFFOUT,FLUSH_BUFFER
!
!  Local variables
!

      INTEGER, POINTER   :: IPPLOOK(:,:)
      INTEGER, PARAMETER :: current_io_pe=0

      INTEGER            :: DUMMY
      INTEGER            :: STEP



!
!dir$ cache_align pp_fixhd, pp_inthd, pp_realhd, pp_levdepc, 
!     pp_rowdepc, pp_coldepc, ipplook
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
      INTEGER                                                           &
     &       II,JJ,IWA,IX,LEN_IO,START_BLOCK  !
      REAL A_IO
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

!L----------------------------------------------------------------------
!L 1. Reserve space
!L
       NULLIFY(IPPLOOK)
       STEP = 1
       CALL INIT_IPPLOOK(IPPLOOK, FTN_UNIT, LEN1_LOOKUP,                &
     &                   PP_LEN2_LOOKUP, DUMMY, STEP)

!L----------------------------------------------------------------------
!L 1.1 Set up FIXED header record for the PP FILE
!L
! Attach fixed length header
      CALL INIT_FXH(FTN_UNIT,LEN_FIXHD)
      CALL ATTACH_FXH(PP_FIXHD,FTN_UNIT)
      DO 3 II=1,LEN_FIXHD
      PP_FIXHD(II)=FIXHD(II)
    3 CONTINUE
      IF (FILE_TYPE_LETTER == 'p' .OR.                                  &
     &    FILE_TYPE_LETTER == 'c') THEN
        PP_FIXHD(5)=3
      ELSEIF (FILE_TYPE_LETTER == 'b') THEN
        PP_FIXHD(5)=5
      ELSE
        ICODE=100
        CMESSAGE='INIT_PP  : Unknown output file type letter'
        RETURN
      ENDIF
      PP_FIXHD(101)=PP_LEN_INTHD
      PP_FIXHD(105)=PP_FIXHD(100)+PP_FIXHD(101)
      PP_FIXHD(106)=PP_LEN_REALHD
      PP_FIXHD(110)=PP_FIXHD(105)+PP_FIXHD(106)
      PP_FIXHD(111)=LEN1_LEVDEPC
      PP_FIXHD(112)=LEN2_LEVDEPC
      PP_FIXHD(115)=0
      PP_FIXHD(116)=rmdi
      PP_FIXHD(117)=rmdi
      PP_FIXHD(120)=0
      PP_FIXHD(121)=rmdi
      PP_FIXHD(122)=rmdi
      PP_FIXHD(125)=0
      PP_FIXHD(126)=rmdi
      PP_FIXHD(127)=rmdi
      PP_FIXHD(130)=0
      PP_FIXHD(131)=rmdi
      PP_FIXHD(135)=0
      PP_FIXHD(136)=rmdi
      PP_FIXHD(140)=0
      PP_FIXHD(141)=rmdi
      PP_FIXHD(142)=0
      PP_FIXHD(143)=rmdi
      PP_FIXHD(144)=0
      PP_FIXHD(145)=rmdi
      PP_FIXHD(150)=PP_FIXHD(110)+ PP_FIXHD(111)*PP_FIXHD(112)
      IF (LEN2_ROWDEPC > 0) THEN
        PP_FIXHD(115)=PP_FIXHD(110)+ PP_FIXHD(111)*PP_FIXHD(112)
        PP_FIXHD(116)=LEN1_ROWDEPC 
        PP_FIXHD(117)=LEN2_ROWDEPC
        PP_FIXHD(150)=PP_FIXHD(115)+ PP_FIXHD(116)*PP_FIXHD(117)  
      END IF 
      IF (LEN2_COLDEPC > 0) THEN 
        PP_FIXHD(120)=PP_FIXHD(115)+ PP_FIXHD(116)*PP_FIXHD(117)
        PP_FIXHD(121)=LEN1_COLDEPC
        PP_FIXHD(122)=LEN2_COLDEPC
        PP_FIXHD(150)=PP_FIXHD(120)+ PP_FIXHD(121)*PP_FIXHD(122) 
      END IF
      PP_FIXHD(151)=LEN1_LOOKUP
      PP_FIXHD(152)=PP_LEN2_LOOKUP
      pp_fixhd(160)=                                                    &
                         ! make sure the data starts on a sector bndry
     & ((pp_fixhd(150)+pp_len2_lookup*len1_lookup-1+um_sector_size-1)/  &
     & um_sector_size)*um_sector_size+1
!L----------------------------------------------------------------------
!L 1.2 Set up INTEGER constants record for the PP FILE
!L
      IF(PP_FIXHD(5) <= 2) THEN !  set all values initially to MDI
        DO II=1,PP_LEN_INTHD
          PP_INTHD(II)=INTHD(21)
        ENDDO
      ELSE
        DO II=1,PP_LEN_INTHD
          PP_INTHD(II)=IMDI
        ENDDO
      ENDIF

      PP_INTHD(6)=INTHD(6)
      PP_INTHD(7)=INTHD(7)
      PP_INTHD(8)=INTHD(8)
      PP_INTHD(9)=INTHD(9)
      PP_INTHD(10)=INTHD(10)
      PP_INTHD(13)=INTHD(13)
      PP_INTHD(17)=INTHD(17)
      PP_INTHD(24)=INTHD(24)
      PP_INTHD(25)=INTHD(25)
      PP_INTHD(28)=INTHD(28)
!L----------------------------------------------------------------------
!L 1.3 Set up REAL constants record for the PP FILE
!L
      DO II = 1, PP_LEN_REALHD
        PP_REALHD(II) = RMDI   ! Set all values to RMDI initially
      END DO
      PP_REALHD(1)=REALHD(1)
      PP_REALHD(2)=REALHD(2)
      PP_REALHD(3)=REALHD(3)
      PP_REALHD(4)=REALHD(4)
! Set to RMDI for VR      
      IF (LEN2_ROWDEPC > 0 .AND. LEN2_COLDEPC > 0) THEN
        PP_REALHD(1) = RMDI
        PP_REALHD(2) = RMDI     
        PP_REALHD(3) = RMDI    
        PP_REALHD(4) = RMDI       
      ENDIF      
      PP_REALHD(5)=REALHD(5)
      PP_REALHD(6)=REALHD(6)
      PP_REALHD(16)=REALHD(16)
      PP_REALHD(17)=REALHD(17)

!L----------------------------------------------------------------------
!L 1.4 Set up LEVEL/ROW/COL DEPENDANT constants record for the PP FILE
!L
      DO 5 II=1,LEN1_LEVDEPC*LEN2_LEVDEPC
      PP_LEVDEPC(II)=LEVDEPC(II)
    5 CONTINUE
      
      DO II=1,LEN1_ROWDEPC*LEN2_ROWDEPC
         PP_ROWDEPC(II)=ROWDEPC(II)
      END DO
      
      DO II=1,LEN1_COLDEPC*LEN2_COLDEPC
         PP_COLDEPC(II)=COLDEPC(II)
      END DO
       
!L----------------------------------------------------------------------
!L 2.1 BUFFER OUT Header Records starting with the FIXED LENGTH
!L
! DEPENDS ON: buffout
      CALL BUFFOUT(FTN_UNIT,PP_FIXHD(1),LEN_FIXHD,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
           CALL IOERROR('bufferout of fixed length header',A_IO,LEN_IO, &
     &                    LEN_FIXHD)
           CMESSAGE='INIT_PP:I/O error'
           ICODE=1
           RETURN
        ENDIF
      START_BLOCK=LEN_FIXHD+1
!L----------------------------------------------------------------------
!L 2.2 BUFFER OUT Integer Constants
!L

      IF(FIXHD(100) >  0) THEN  ! Any integer constants to output ?

! Check for error in file pointers

!        WRITE(6,*)  'START_BLOCK FIXHD(100)'
!        WRITE(6,*)   START_BLOCK
!        WRITE(6,*)   FIXHD(100)
!        WRITE(6,*)   FTN_UNIT
         IF(FIXHD(100) /= START_BLOCK) THEN  ! Check start address
! DEPENDS ON: poserror
            CALL POSERROR('integer constants',START_BLOCK,100,          &
     &      PP_FIXHD(100))
            CMESSAGE='INIT_PP:  Addressing conflict'
            ICODE=2
            RETURN
         END IF

! DEPENDS ON: buffout
         CALL BUFFOUT (FTN_UNIT,PP_INTHD(1),PP_FIXHD(101),LEN_IO,A_IO)

! Check for I/O errors

         IF(A_IO /= -1.0.OR.LEN_IO /= PP_FIXHD(101)) THEN
! DEPENDS ON: ioerror
            CALL IOERROR('buffer out of integer constants',A_IO,LEN_IO  &
     &     ,PP_FIXHD(101))
            CMESSAGE='INIT_PP: I/O Error'
            ICODE=3
            RETURN
         END IF

         START_BLOCK=START_BLOCK+PP_FIXHD(101)

      END IF

!L----------------------------------------------------------------------
!L 2.3 BUFFER OUT Real Constants
!L

      IF(PP_FIXHD(105) >  0) THEN   ! Any real constants to output ?

! Check for error in file pointers

        IF(PP_FIXHD(105) /= START_BLOCK) THEN
! DEPENDS ON: poserror
          CALL POSERROR('real constants',START_BLOCK,100,PP_FIXHD(105))
          CMESSAGE='INIT_PP: Addressing conflict'
          ICODE=4
          RETURN
        END IF

! DEPENDS ON: buffout
        CALL BUFFOUT(FTN_UNIT,PP_REALHD(1),PP_FIXHD(106),LEN_IO,A_IO)

! Check for I/O errors

        IF(A_IO /= -1.0.OR.LEN_IO /= PP_FIXHD(106)) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of real constants',A_IO,LEN_IO       &
     &                 ,PP_FIXHD(106))
          CMESSAGE='INIT_PP: I/O Error'
          ICODE=5
          RETURN
        END IF

        START_BLOCK=START_BLOCK+PP_FIXHD(106)

      END IF

!L----------------------------------------------------------------------
!L 2.4.1 BUFFER OUT Level Dependant Constants.
!L

      IF(PP_FIXHD(112) >  0) THEN ! Any level dependant constants ?

! Check for error in file pointers

         IF(PP_FIXHD(110) /= START_BLOCK) THEN
! DEPENDS ON: poserror
            CALL POSERROR('real constants',START_BLOCK,100,             &
     &                     PP_FIXHD(110))
            CMESSAGE='INIT_PP: Addressing conflict'
            ICODE=6
            RETURN
         END IF

! DEPENDS ON: buffout
         CALL BUFFOUT (FTN_UNIT,PP_LEVDEPC(1)                           &
     &              ,PP_FIXHD(111)*PP_FIXHD(112),LEN_IO,A_IO)

! Check for I/O errors

         IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(111)*PP_FIXHD(112)      &
     &        ))THEN
! DEPENDS ON: ioerror
           CALL IOERROR('buffer out of lev dep constants',A_IO,LEN_IO   &
     &            ,PP_FIXHD(111))
           CMESSAGE='INIT_PP: I/O Error'
           ICODE=7
           RETURN
         END IF

         START_BLOCK=START_BLOCK+ PP_FIXHD(111)*PP_FIXHD(112)

      END IF
!L----------------------------------------------------------------------
!L 2.4.2 BUFFER OUT Row Dependant Constants.
!L

      IF(PP_FIXHD(115) >  0) THEN ! Any row dependant constants ?

! Check for error in file pointers

        IF(PP_FIXHD(115) /= START_BLOCK) THEN
! DEPENDS ON: poserror
          CALL POSERROR('real constants',START_BLOCK,100,               &
     &                   PP_FIXHD(115))
          CMESSAGE='INIT_PP: Addressing conflict'
          ICODE=6
          RETURN
        END IF

! DEPENDS ON: buffout
        CALL BUFFOUT (FTN_UNIT,PP_ROWDEPC(1),                           &
     &                PP_FIXHD(116)*PP_FIXHD(117),LEN_IO,A_IO)

! Check for I/O errors

        IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(116)*PP_FIXHD(117)       &
     &       ))THEN

! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of row dep constants',A_IO,LEN_IO,   &
     &                  PP_FIXHD(116))
          CMESSAGE='INIT_PP: I/O Error'
          ICODE=7
          RETURN
        END IF

        START_BLOCK=START_BLOCK+ PP_FIXHD(116)*PP_FIXHD(117)

      END IF
!L----------------------------------------------------------------------
!L 2.4.3 BUFFER OUT Col Dependant Constants.
!L

      IF(PP_FIXHD(120) >  0) THEN ! Any col dependant constants ?

! Check for error in file pointers

        IF(PP_FIXHD(120) /= START_BLOCK) THEN
! DEPENDS ON: poserror
          CALL POSERROR('real constants',START_BLOCK,100,               &
     &                   PP_FIXHD(120))
          CMESSAGE='INIT_PP: Addressing conflict'
          ICODE=6
          RETURN
        END IF

! DEPENDS ON: buffout
        CALL BUFFOUT (FTN_UNIT,PP_COLDEPC(1),                           &
     &                PP_FIXHD(121)*PP_FIXHD(122),LEN_IO,A_IO)

! Check for I/O errors

        IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(121)*PP_FIXHD(122)       &
     &       ))THEN

! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of col dep constants',A_IO,LEN_IO ,  &
     &                  PP_FIXHD(121))
          CMESSAGE='INIT_PP: I/O Error'
          ICODE=7
          RETURN
        END IF

        START_BLOCK=START_BLOCK+ PP_FIXHD(121)*PP_FIXHD(122)

      END IF      
!L----------------------------------------------------------------------
!L 2.5 BUFFER OUT Lookup Table
!L
!     IWA= 0
!     CALL SETPOS(FTN_UNIT,3,IWA,ICODE)
           IF(PP_FIXHD(152) >  0) THEN

! Check for error in file pointers

             IF(PP_FIXHD(150) /= START_BLOCK) THEN
! DEPENDS ON: poserror
               CALL POSERROR('lookup table',START_BLOCK,100,            &
     &              PP_FIXHD(150))
               CMESSAGE='INIT_PP: Addressing conflict'
               ICODE=8
               RETURN
             END IF

! DEPENDS ON: buffout
      CALL BUFFOUT (FTN_UNIT,                                           &
     &              IPPLOOK,LEN1_LOOKUP*PP_LEN2_LOOKUP,LEN_IO,A_IO)

!
! Check for I/O errors

            IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(151)*PP_FIXHD(152))) &
     &          THEN
! DEPENDS ON: ioerror
              CALL IOERROR('buffer out of PP LOOKUP TABLE ',A_IO,LEN_IO &
     &            ,PP_FIXHD(152))
              CMESSAGE='INIT_PP: I/O Error'
              ICODE=9
              RETURN
            END IF
!
! Clear file buffer : force last buffer to be written to file
!  to avoid problems with continuation runs following hard failures.
!
      CALL FLUSH_BUFFER(FTN_UNIT,ICODE)
      IF(ICODE /= 0) THEN
         CMESSAGE='INIT_PP: Problem flushing buffer'
         ICODE=10
         RETURN
      ENDIF
!
            START_BLOCK=START_BLOCK+(PP_FIXHD(151)*PP_FIXHD(152))
!
! If we are the current I/O PE, we need to update our copy
! of the LOOKUP Table disk address
!
      STEP = 2
      CALL INIT_IPPLOOK(IPPLOOK, FTN_UNIT, DUMMY, DUMMY,                &
     &                  PP_FIXHD(150) - 1, STEP)
      NULLIFY(PP_FIXHD)

          END IF
      RETURN
      END SUBROUTINE INIT_PP
