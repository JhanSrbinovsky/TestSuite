
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE oasis3_puta2o(                                         &
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Start arg_atm_fields.h

! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "artptra.h" and "argptra.h".
!
! Current Code Owner: A. Treshansky
!


! End arg_atm_fields.h
  &  PUT_STEP,cmessage)

  IMPLICIT NONE

  !
  ! Description:
  ! Dummy version of oasis3_puta2o - should never be called.
  !
  ! Author: R. Hill
  !
  ! Current Code Owner : R. Hill
  !
  !=====================================================================

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
!     o GLOBAL data - the entire data domain processed by the UM
!     o LOCAL data - the fragment of the GLOBAL data which is
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
      ! ATMOS START
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

      ! Data structure sizes for ATMOS & OCEAN BOUNDARY file control
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
! TYPD1 Common block containing the ALT_N_SUBMODEL_PARTITION variables
! CALTSUBM
! TYPD1 needs access to N_SUBMODEL_PARTITION/_MAX in CSUBMODL. However,
! they are not always called in the same decks and in the right order.
! Therefore, copy the values to another file and include it from TYPD1

      INTEGER ALT_N_SUBMODEL_PARTITION

      INTEGER, PARAMETER :: ALT_N_SUBMODEL_PARTITION_MAX=4

      COMMON/CALTSUBM/ALT_N_SUBMODEL_PARTITION
! CALTSUBM end
! This file needs TYPSIZE included first

      REAL    ::  D1(LEN_TOT)       ! IN/OUT: Main data array
      LOGICAL :: LD1(LEN_TOT)       ! IN/OUT: Main data array (logical)
      INTEGER :: ID1(LEN_TOT)       ! I/OUT: Main data array (integer)

! D1_ADDR start
      ! Information for accessing D1 addressing array
      ! Number of items of info needed for each object and maximum
      ! number of objects in D1 -

      ! Number of items of information in D1 addressing array
      INTEGER,PARAMETER:: D1_LIST_LEN=17

! Names of items in D1 addressing array. Update D1_LIST_LEN above if
! items added

      ! Prognostic, Diagnostic, Secondary or other
      INTEGER,PARAMETER:: d1_object_type    = 1 ! Internal model id
      INTEGER,PARAMETER:: d1_imodl          = 2  ! Internal model id
      INTEGER,PARAMETER:: d1_section        = 3  ! Section
      INTEGER,PARAMETER:: d1_item           = 4  ! Item
      INTEGER,PARAMETER:: d1_address        = 5  ! Address in D1
      INTEGER,PARAMETER:: d1_length         = 6  ! Record length
      INTEGER,PARAMETER:: d1_grid_type      = 7  ! Grid type
      INTEGER,PARAMETER:: d1_no_levels      = 8  ! Number of levels

      ! Stash list number for diags. -1 for progs
      INTEGER,PARAMETER:: d1_stlist_no      = 9

      ! Pointer to dump header lookup table
      INTEGER,PARAMETER:: d1_lookup_ptr     = 10

      INTEGER,PARAMETER:: d1_north_code     = 11 ! Northern row
      INTEGER,PARAMETER:: d1_south_code     = 12 ! Southern row
      INTEGER,PARAMETER:: d1_east_code      = 13 ! Eastern row
      INTEGER,PARAMETER:: d1_west_code      = 14 ! Western row
      INTEGER,PARAMETER:: d1_gridpoint_code = 15 ! gridpoint info
      INTEGER,PARAMETER:: d1_proc_no_code   = 16 ! Processing Code
      INTEGER,PARAMETER:: d1_halo_type      = 17 ! Halo width type

      ! Types of items for d1_type

      INTEGER,PARAMETER:: prognostic = 0
      INTEGER,PARAMETER:: diagnostic = 1
      INTEGER,PARAMETER:: secondary  = 2
      INTEGER,PARAMETER:: other      = 3

! D1_ADDR end
      ! D1 addressing array and number of objects in each submodel
      INTEGER :: D1_ADDR(D1_LIST_LEN,N_OBJ_D1_MAX,                      &
     &  ALT_N_SUBMODEL_PARTITION)

      INTEGER :: NO_OBJ_D1(ALT_N_SUBMODEL_PARTITION_MAX)

      COMMON/common_D1_ADDRESS/ NO_OBJ_D1
! TYPD1 end
! CTRACERA start
!  Vn    Date    Modification History
! 6.1  23/06/04  Prognostic tracers now in section 33, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.
! 6.2  13/07/05  Also increase A_MAX_TRVARS to 150. R Barnes.
! 6.2  10/11/05  UKCA tracers put into section 34, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.

      ! First atmospheric tracer (STASH No)
      INTEGER,PARAMETER:: A_TRACER_FIRST = 1

      ! Last atmospheric tracer  (STASH No)
      INTEGER,PARAMETER:: A_TRACER_LAST = 150

      ! Maximum number of atmospheric tracers
      INTEGER,PARAMETER:: A_MAX_TRVARS  = 150
!     First UKCA tracer (STASH No)
      INTEGER,PARAMETER:: A_UKCA_FIRST = 1
!     Last UKCA tracer  (STASH No)
      INTEGER,PARAMETER:: A_UKCA_LAST = 150
!     Maximum number of UKCA tracers
      INTEGER,PARAMETER:: A_MAX_UKCAVARS  = 150

      ! Index to relative position.
      ! A_TR_INDEX(N) gives position in JTRACER for tracer number N.
      ! Set in SET_ATM_POINTERS.
      ! A_TR_INDEX(N) is the position, in the list of tracers
      ! actually present in D1, that tracer number N (in the list
      ! of all tracers selectable from the user interface) occupies,
      ! if it is present.
      ! If tracer number N is absent then A_TR_INDEX(N) is undefined.

      INTEGER :: A_TR_INDEX(A_MAX_TRVARS)
      INTEGER :: A_UKCA_INDEX(A_MAX_UKCAVARS)   ! Similar to A_TR_INDEX

      COMMON/ATRACER/A_TR_INDEX,A_UKCA_INDEX

! CTRACERA end
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!

! Start typ_atm_fields.h

! Description:
!  Contains set of atmospheric fields to be used as arguments to subroutines
!  without referring explicitly to D1 or "jpointers".
!  This file should replace "typptra.h", and requires that "ctracera.h" is
!  also included.  By 7.1 these should all be their natural shape instead
!  of all being 1D arrays.
!
! Current Code Owner: A. Treshansky
!

      ! 1.1: Data variables stored in primary space.
      real :: U(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: V(1-offx:row_length+offx,1-offy:n_rows+offy,model_levels)
      real :: W(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)
      real :: RHO(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: THETA(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      real :: Q(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: QCL(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: QCF(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: QCF2(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: QRAIN(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: QGRAUP(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)

! EXNER_RHO_LEVELS is still 1D; change for 7.1
      ! Exner pressure on rho levels
      real :: EXNER_RHO_LEVELS(((2*offx)+row_length)*((2*offy)+rows)*(model_levels+1))   

      real :: U_ADV(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: V_ADV(1-offx:row_length+offx,1-offy:n_rows+offy,model_levels)
      real :: W_ADV(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

! P is still 1D; change for 7.1
      ! 1.2: Data variables stored in secondary space.
      real :: P(((2*offx)+row_length)*((2*offy)+rows)*(model_levels+1))

      ! Pressure on theta levels
      real :: P_THETA_LEVELS(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      ! Exner pressure on theta levels
      real :: EXNER_THETA_LEVELS(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      ! 1.3: Cloud Fields
      real :: CCW_RAD(row_length,rows,wet_levels)
      ! CCA's size is dependant on L_3D_CCA
      ! N_CCA_LEV will be set to the correct value (either wet_levels or 1)
      real :: CCA(row_length,rows,n_cca_lev)
      real :: CF_AREA(row_length,rows,wet_levels)
      real :: CF_BULK(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: CF_LIQUID(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)
      real :: CF_FROZEN(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,wet_levels)

      ! 1.4: Soil Ancillary fields
      real :: DEEP_SOIL_TEMP(land_field,sm_levels)
      real :: SMCL(land_field,sm_levels)
      real :: STHU(land_field,sm_levels)
      real :: STHF(land_field,sm_levels)

      ! 1.5: Radiation Increments
      real :: SW_INCS(row_length,rows,0:model_levels+1)  ! SW radiation increments
      real :: LW_INCS(row_length,rows,0:model_levels)    ! LW radiation increments
! PAR radiation increment
      real :: DIRPAR(row_length,rows)

! OZONE fields are still 1D; change for 7.1
      ! 1.6: Ozone
      real :: O3(row_length*rows*ozone_levels)
!  tropopause-based ozone
      real :: TPPSOZONE(row_length*rows*ozone_levels)
!  Cariolle ozone tracer variables
      real :: OZONE_TRACER(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: O3_PROD_LOSS(1,rows,model_levels)
      real :: O3_P_L_VMR(1,rows,model_levels)
      real :: O3_VMR(1,rows,model_levels)
      real :: O3_P_L_TEMP(1,rows,model_levels)
      real :: O3_TEMP(1,rows,model_levels)
      real :: O3_P_L_COLO3(1,rows,model_levels)
      real :: O3_COLO3(1,rows,model_levels)
      
! TRACERS are still 1D; change for 7.1
      ! 1.7: Tracer and aerosol fields
      ! TRACERS are dealt w/ differently
      ! these are the maximum sizes:
      real :: TRACER(tr_levels*theta_off_size*(a_tracer_last-a_tracer_first))
      real :: TRACER_UKCA(tr_levels*theta_off_size*(a_ukca_last-a_ukca_first))
      real :: TRACER_LBC(LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_vars)
      real :: TRACER_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels*tr_vars)

      real :: MURK_SOURCE(row_length,rows,model_levels)
      real :: MURK(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV3(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV4(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV5(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DUST_DIV6(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SO2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: DMS(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SO4_AITKEN(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SO4_ACCU(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SO4_DISS(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: H2O2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: NH3(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SOOT_NEW(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SOOT_AGD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SOOT_CLD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: BMASS_NEW(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: BMASS_AGD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: BMASS_CLD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: OCFF_NEW(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: OCFF_AGD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: OCFF_CLD(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: SO2_NATEM(row_length,rows,model_levels)
      real :: OH(row_length,rows,model_levels)
      real :: HO2(row_length,rows,model_levels)
      real :: H2O2_LIMIT(row_length,rows,model_levels)
      real :: O3_CHEM(row_length,rows,model_levels)
      real :: CO2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: CH4_STOCH(row_length,rows,model_levels)
      real :: O3_STOCH(row_length,rows,model_levels)

! USER_MULT<N> are still 1D; change for 7.1
! 1.8: Multi-level user ancillary fields
      real :: USER_MULT1(row_length*rows*model_levels)
      real :: USER_MULT2(row_length*rows*model_levels)
      real :: USER_MULT3(row_length*rows*model_levels)
      real :: USER_MULT4(row_length*rows*model_levels)
      real :: USER_MULT5(row_length*rows*model_levels)
      real :: USER_MULT6(row_length*rows*model_levels)
      real :: USER_MULT7(row_length*rows*model_levels)
      real :: USER_MULT8(row_length*rows*model_levels)
      real :: USER_MULT9(row_length*rows*model_levels)
      real :: USER_MULT10(row_length*rows*model_levels)
      real :: USER_MULT11(row_length*rows*model_levels)
      real :: USER_MULT12(row_length*rows*model_levels)
      real :: USER_MULT13(row_length*rows*model_levels)
      real :: USER_MULT14(row_length*rows*model_levels)
      real :: USER_MULT15(row_length*rows*model_levels)
      real :: USER_MULT16(row_length*rows*model_levels)
      real :: USER_MULT17(row_length*rows*model_levels)
      real :: USER_MULT18(row_length*rows*model_levels)
      real :: USER_MULT19(row_length*rows*model_levels)
      real :: USER_MULT20(row_length*rows*model_levels)

!     CABLE
      real :: TSOIL_TILE(land_field,ntiles,st_levels)
      real :: SMCL_TILE(land_field,ntiles,sm_levels)
! Not used MRD
!      real :: STHU_TILE(land_field,ntiles,sm_levels)
      real :: STHF_TILE(land_field,ntiles,sm_levels)
      real :: SNOW_DEPTH3L(land_field,ntiles,3)
      real :: SNOW_MASS3L(land_field,ntiles,3)
      real :: SNOW_TMP3L(land_field,ntiles,3)
      real :: SNOW_RHO3L(land_field,ntiles,3)
      real :: SNOW_RHO1L(land_field,ntiles)
      real :: SNOW_AGE(land_field,ntiles)
      real :: SNOW_FLG3L(land_field,ntiles)

! Lestevens Sept 2012 - For CASA-CNP
      real :: CPOOL_TILE(land_field,ntiles,10)
      real :: NPOOL_TILE(land_field,ntiles,10)
      real :: PPOOL_TILE(land_field,ntiles,12)
      real :: SOIL_ORDER(land_field)
      real :: NIDEP(land_field)
      real :: NIFIX(land_field)
      real :: PWEA(land_field)
      real :: PDUST(land_field)
      real :: GLAI(land_field,ntiles)
      real :: PHENPHASE(land_field,ntiles)

      real :: cable_lai(land_field,ntiles)
      
      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
      real :: OROG_LBC(LENRIMA(fld_type_p,halo_type_extended,1))
      real :: U_LBC(LENRIMA(fld_type_u,halo_type_extended,1),model_levels)
      real :: V_LBC(LENRIMA(fld_type_v,halo_type_extended,1),model_levels)
      real :: W_LBC(LENRIMA(fld_type_p,halo_type_extended,1),0:model_levels)
      real :: RHO_LBC(LENRIMA(fld_type_p,halo_type_extended,1),model_levels)
      real :: THETA_LBC(LENRIMA(fld_type_p,halo_type_extended,1),model_levels)
      real :: Q_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCL_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCF_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCF2_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QRAIN_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QGRAUP_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_BULK_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_LIQUID_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_FROZEN_LBC(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: EXNER_LBC(LENRIMA(fld_type_p,halo_type_extended,1),model_levels+1)
      real :: U_ADV_LBC(LENRIMA(fld_type_u,halo_type_extended,1),model_levels)
      real :: V_ADV_LBC(LENRIMA(fld_type_v,halo_type_extended,1),model_levels)
      real :: W_ADV_LBC(LENRIMA(fld_type_p,halo_type_extended,1),0:model_levels)
      real :: MURK_LBC(LENRIMA(fld_type_p,halo_type_single,1),model_levels)
      
      ! Lateral Boundary Condition tendencies
      real :: U_LBC_TEND(LENRIMA(fld_type_u,halo_type_extended,1),model_levels)
      real :: V_LBC_TEND(LENRIMA(fld_type_v,halo_type_extended,1),model_levels)
      real :: W_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),0:model_levels)
      real :: RHO_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),model_levels)
      real :: THETA_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),model_levels)
      real :: Q_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCL_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCF_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QCF2_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QRAIN_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: QGRAUP_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_BULK_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_LIQUID_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: CF_FROZEN_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),wet_levels)
      real :: EXNER_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),model_levels+1)
      real :: U_ADV_LBC_TEND(LENRIMA(fld_type_u,halo_type_extended,1),model_levels)
      real :: V_ADV_LBC_TEND(LENRIMA(fld_type_v,halo_type_extended,1),model_levels)
      real :: W_ADV_LBC_TEND(LENRIMA(fld_type_p,halo_type_extended,1),0:model_levels)
      real :: MURK_LBC_TEND(LENRIMA(fld_type_p,halo_type_single,1),model_levels)

      ! 2: Scalar Variables

      ! 2.1: Data variables stored in primary space.
      real :: TSTAR(row_length,rows)
      real :: LAND(row_length,rows)
      real :: TSTAR_ANOM(row_length,rows)
      ! 2.15: Fields for coastal tiling
      real :: FRAC_LAND(land_field)
      real :: TSTAR_LAND(row_length,rows)
      real :: TSTAR_SEA(row_length,rows)
      real :: TSTAR_SICE(row_length,rows)
      ! Set pointers for sea-ice and land albedos
      real :: SICE_ALB(row_length,rows)
      real :: LAND_ALB(row_length,rows)

      ! 2.2: Data variables stored in secondary space.

      real :: PSTAR(row_length,rows)

      ! 2.3: Cloud fields
      real :: LCBASE(row_length,rows)
      real :: CCB(row_length,rows)
      real :: CCT(row_length,rows)
      real :: CCLWP(row_length,rows)

      ! 2.4: Boundary layer fields

      real :: ZH(row_length,rows)

      ! Standard deviation of turbulent fluctuations of layer 1 temperature
      real :: T1_SD(row_length,rows)

      ! Standard deviation of turbulent fluctuations of layer 1 humidity
      real :: Q1_SD(row_length,rows)

      ! Number of model levels in the  turbulently mixed layer
      real :: NTML(row_length,rows)

      ! Top level for turb mixing in any decoupled Sc layer
      real :: NTDSC(row_length,rows)

      ! Bottom level for turb mixing in any decoupled Sc layer
      real :: NBDSC(row_length,rows)

      real :: CUMULUS(row_length,rows)

      ! 2.4: Soil Ancillary fields

      real :: SAT_SOILW_SUCTION(land_field)
      real :: THERM_CAP(land_field)
      real :: THERM_COND(land_field)
      real :: VOL_SMC_CRIT(land_field)
      real :: VOL_SMC_WILT(land_field)
      real :: VOL_SMC_SAT(land_field)
      real :: SAT_SOIL_COND(land_field)
      real :: CLAPP_HORN(land_field)

      ! 2.5: Vegetation Ancillary fields

      real :: Z0(row_length,rows)
      real :: CANOPY_WATER(land_field)
      real :: SURF_CAP(row_length,rows)
      real :: SURF_RESIST(land_field)
      real :: ROOT_DEPTH(land_field)
      real :: INFILT(land_field)
      real :: VEG_FRAC(land_field)
      real :: LAI(land_field)
      real :: CANHT(land_field)
      real :: SFA(land_field)
      real :: MDSA(land_field)
      real :: GS(land_field)

      ! 2.6: Orographic Ancillary fields

      real :: OROGRAPHY(row_length,rows)
      real :: OROG_SD(land_field)
      real :: OROG_SIL(land_field)
      real :: OROG_HO2(land_field)
      real :: OROG_GRAD_X(land_field)
      real :: OROG_GRAD_Y(land_field)
      real :: OROG_GRAD_XX(land_field)
      real :: OROG_GRAD_XY(land_field)
      real :: OROG_GRAD_YY(land_field)

      ! 2.7: Sea/Sea Ice

      real :: U_SEA(row_length,rows)
      real :: V_SEA(row_length,n_rows)
      real :: U_0_P(row_length,rows)
      real :: V_0_P(row_length,rows)

! these fields are targets b/c there are pointers w/in ATM_STEP
! which point to one or another of them depending on the ice category selection
      real, target :: TI(row_length,rows,1)
      real, target :: ICE_FRACTION(row_length,rows,1)
      real, target :: ICE_THICKNESS(row_length,rows,1)
      real, target :: TI_CAT(row_length,rows,nice) 
      real, target :: ICE_FRACT_CAT(row_length,rows,nice)
      real, target :: ICE_THICK_CAT(row_length,rows,nice)

      ! 2.8: Snow

      real :: SNODEP(row_length,rows)
      real :: SNODEP_SEA(row_length,rows)
      real :: SNODEP_SEA_CAT(row_length,rows,nice)
      real :: CATCH_SNOW(land_field)
      real :: SNOW_GRND(land_field)
      ! SNSOOT may not be used as of vn6.6
      real :: SNSOOT(row_length,rows)  ! Snow soot content

      ! 2.9: aerosol emission fields, including mineral dust parent soil props

      real :: SOIL_CLAY(row_length,rows)
      real :: SOIL_SILT(row_length,rows)
      real :: SOIL_SAND(row_length,rows)
      real :: DUST_MREL1(row_length,rows)
      real :: DUST_MREL2(row_length,rows)
      real :: DUST_MREL3(row_length,rows)
      real :: DUST_MREL4(row_length,rows)
      real :: DUST_MREL5(row_length,rows)
      real :: DUST_MREL6(row_length,rows)

      real :: SO2_EM(row_length,rows)
      real :: DMS_EM(row_length,rows)
      real :: SO2_HILEM(row_length,rows)
      real :: NH3_EM(row_length,rows)
      real :: SOOT_EM(row_length,rows)
      real :: SOOT_HILEM(row_length,rows)
      real :: BMASS_EM(row_length,rows)
      real :: BMASS_HILEM(row_length,rows)
      real :: OCFF_EM(row_length,rows)
      real :: OCFF_HILEM(row_length,rows)
      real :: DMS_CONC(row_length,rows)
      real :: DMS_OFLUX(row_length,rows)

! USER_ANC<N> fields are still 1D; change for 7.1
      ! 2.10: User ancillary fields
      real :: USER_ANC1(row_length*rows*1)
      real :: USER_ANC2(row_length*rows*1)
      real :: USER_ANC3(row_length*rows*1)
      real :: USER_ANC4(row_length*rows*1)
      real :: USER_ANC5(row_length*rows*1)
      real :: USER_ANC6(row_length*rows*1)
      real :: USER_ANC7(row_length*rows*1)
      real :: USER_ANC8(row_length*rows*1)
      real :: USER_ANC9(row_length*rows*1)
      real :: USER_ANC10(row_length*rows*1)
      real :: USER_ANC11(row_length*rows*1)
      real :: USER_ANC12(row_length*rows*1)
      real :: USER_ANC13(row_length*rows*1)
      real :: USER_ANC14(row_length*rows*1)
      real :: USER_ANC15(row_length*rows*1)
      real :: USER_ANC16(row_length*rows*1)
      real :: USER_ANC17(row_length*rows*1)
      real :: USER_ANC18(row_length*rows*1)
      real :: USER_ANC19(row_length*rows*1)
      real :: USER_ANC20(row_length*rows*1)

      ! Tracer fluxes - kdcorbin, 05/10
      real :: TRACER_FLUX1(row_length,rows)
      real :: TRACER_FLUX2(row_length,rows)
      real :: TRACER_FLUX3(row_length,rows)
      real :: TRACER_FLUX4(row_length,rows)
      real :: TRACER_FLUX5(row_length,rows)
      real :: TRACER_FLUX6(row_length,rows)
      real :: TRACER_FLUX7(row_length,rows)
      real :: TRACER_FLUX8(row_length,rows)
      real :: TRACER_FLUX9(row_length,rows)
      real :: TRACER_FLUX10(row_length,rows)
      real :: TRACER_FLUX11(row_length,rows)
      real :: TRACER_FLUX12(row_length,rows)
      real :: TRACER_FLUX13(row_length,rows)
      real :: TRACER_FLUX14(row_length,rows)
      real :: TRACER_FLUX15(row_length,rows)
      real :: TRACER_FLUX16(row_length,rows)
      real :: TRACER_FLUX17(row_length,rows)
      real :: TRACER_FLUX18(row_length,rows)
      real :: TRACER_FLUX19(row_length,rows)
      real :: TRACER_FLUX20(row_length,rows)

      !   2.11: Store arrays for energy correction calculation
      real :: NET_FLUX(row_length,rows)
      real :: NET_MFLUX(row_length,rows)

      !   2.12: Tiled Vegetation and Triffid fields
      real :: FRAC_TYP(land_field)
      real :: FRAC_CON1(land_field)  ! fraction of broadleaf tree
      real :: FRAC_CON2(land_field)  ! fraction of needleleaf tree
      real :: FRAC_CON3(land_field)  ! fraction of C3 grass
      real :: FRAC_CON4(land_field)  ! fraction of C4 grass
      real :: FRAC_CON5(land_field)  ! fraction of shrub
      real :: FRAC_CON6(land_field)  ! fraction of urban
      real :: FRAC_CON7(land_field)  ! fraction of water
      real :: FRAC_CON8(land_field)  ! fraction of soil
      real :: FRAC_CON9(land_field)  ! fraction of ice
      real :: LAI_PFT(land_field)
      real :: CANHT_PFT(land_field)
      real :: DISTURB_VEG(land_field)
      real :: SOIL_ALB(land_field)
      real, target :: SOIL_CARB(land_field)
      real, target :: SOIL_CARB1(land_field)
      real, target :: SOIL_CARB2(land_field)
      real, target :: SOIL_CARB3(land_field)
      real, target :: SOIL_CARB4(land_field)
      real :: NPP_PFT_ACC(land_field)
      real :: G_LF_PFT_ACC(land_field)
      real :: G_PHLF_PFT_ACC(land_field)
      real :: RSP_W_PFT_ACC(land_field)
      real, target :: RSP_S_ACC(land_field)
      real, target :: RSP_S_ACC1(land_field)
      real, target :: RSP_S_ACC2(land_field)
      real, target :: RSP_S_ACC3(land_field)
      real, target :: RSP_S_ACC4(land_field)
      real :: CAN_WATER_TILE(land_field)
      real :: CATCH_TILE(land_field)
      real :: INFIL_TILE(land_field)
      real :: RGRAIN_TILE(land_field)
      real :: SNODEP_TILE(land_field)
      real :: TSTAR_TILE(land_field) 
      real :: Z0_TILE(land_field)  
      real :: DOLR_FIELD(row_length,rows) 
      real :: LW_DOWN(row_length,rows) 
      real :: SW_TILE_RTS(land_field)

!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!      real :: TSLAB(row_length*rows*1)
!      real :: TCLIM(row_length*rows*1)
!      real :: HCLIM(row_length*rows*1)
!      real :: CHEAT(row_length*rows*1)
!      real :: OIFLX(row_length*rows*1)
!      real :: UICE(row_length*rows*1)
!      real :: VICE(row_length*n_rows*1)
!      real :: SIG11NE(row_length*rows*1)
!      real :: SIG11SE(row_length*rows*1) 
!      real :: SIG11SW(row_length*rows*1)
!      real :: SIG11NW(row_length*rows*1)
!      real :: SIG12NE(row_length*rows*1)
!      real :: SIG12SE(row_length*rows*1)
!      real :: SIG12NW(row_length*rows*1)
!      real :: SIG22NE(row_length*rows*1)
!      real :: SIG22SE(row_length*rows*1)
!      real :: SIG22SW(row_length*rows*1)
!      real :: SIG22NW(row_length*rows*1)

!   2.14: Carbon cycle fields
      real :: CO2FLUX(row_length,rows)
      real :: CO2_EMITS(row_length,rows)

!   2.15: Fields carried forward from previous version.
!         May not be required
!      real, pointer :: SURF_RESIST_NIT(:)  ! Surface resistance on
!                                    ! non-ice tiles
!      real, pointer :: ROOT_DEPTH_NIT(:)   ! Root depth on non-ice tiles
!      real, pointer :: Z0V_TYP(:)          ! Vegetative Roughness length on
!                                    ! tiles
!      real, pointer :: ICE_EDGE(:)
!      real, pointer :: OROG_TENDENCY(:)    ! Orographic tendencies
!      real, pointer :: OROG_SD_TENDENCY(:) ! Orographic variable tendency

      ! Pointers for ATMOSPHERE model constants. Scalars only.
      ! Addresses in level dependent constants array.
!      real, pointer :: ETATHETA(:)
!      real, pointer :: ETARHO(:)
!      real, pointer :: RHCRIT(:)
!      real, pointer :: SOIL_THICKNESS(:)
! these fields are still 1D; change for 7.1
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
      real :: zseak_theta(0:model_levels)
      real :: Ck_theta(0:model_levels)
      real :: zseak_rho(model_levels)
      real :: Ck_rho(model_levels)   

      ! 2.16: Fields for large-scale hydrology scheme.
      real :: TI_MEAN(land_field)
      real :: TI_SIG(land_field)
      real :: FEXP(land_field)
      real :: GAMMA_INT(land_field)
      real :: WATER_TABLE(land_field)
      real :: FSFC_SAT(land_field)
      real :: F_WETLAND(land_field)

      real :: STHZW(land_field)
      real :: A_FSAT(land_field)
      real :: C_FSAT(land_field)
      real :: A_FWET(land_field)
      real :: C_FWET(land_field)

      ! 2.17: Fields for River routing.
      real :: RIV_SEQUENCE(river_row_length,river_rows)
      real :: RIV_DIRECTION(river_row_length,river_rows)
      real :: RIV_STORAGE(river_row_length,river_rows)
      real :: TOT_SURFROFF(river_row_length,river_rows)
      real :: TOT_SUBROFF(river_row_length,river_rows)
      real :: RIV_INLANDATM(land_field)
      ! Fields for grid-to-grid river routing (river routing 2A)
      real :: RIV_IAREA(row_length,rows)      ! Drainage area
      real :: RIV_SLOPE(row_length,rows)      ! Grid-cell slope
      real :: RIV_FLOWOBS1(row_length,rows)   ! Initial values of flow
      real :: RIV_INEXT(row_length,rows)      ! Flow direction (x)
      real :: RIV_JNEXT(row_length,rows)      ! Flow direction (y)
      real :: RIV_LAND(row_length,rows)       ! Land-type (land/river/sea)
      real :: RIV_SUBSTORE(row_length,rows)   ! Subsurface storage
      real :: RIV_SURFSTORE(row_length,rows)  ! Surface storage
      real :: RIV_FLOWIN(row_length,rows)     ! Surface inflow
      real :: RIV_BFLOWIN(row_length,rows)    ! Subsurface inflow

! Fields used when coupling using OASIS.
      real :: C_SOLAR(row_length,rows)         ! CPL solar radn
      real :: C_BLUE(row_length,rows)          ! CPL blue radn
      real :: C_DOWN(row_length,rows)          ! CPL downward radn
      real :: C_LONGWAVE(row_length,rows)      ! CPL lw radn
      real :: C_TAUX(row_length,rows)          ! CPL taux 
      real :: C_TAUY(row_length,rows)          ! CPL tauy 
      real :: C_WINDMIX(row_length,rows)       ! CPL WME   
      real :: C_SENSIBLE(row_length,rows)      ! CPL sensible ht flx
      real :: C_SUBLIM(row_length,rows)        ! CPL sublim rate
      real :: C_EVAP(row_length,rows)          ! CPL Evap rate
      real :: C_BOTMELTN(row_length,rows,nice) ! CPL Multi-cat bmlt
      real :: C_TOPMELTN(row_length,rows,nice) ! CPL Multi-cat tmlt
      real :: C_LSRAIN(row_length,rows)        ! CPL Lg scl rain rate
      real :: C_LSSNOW(row_length,rows)        ! CPL Lg scl snow rate
      real :: C_CVRAIN(row_length,rows)        ! CPL Cnvctv rain rate
      real :: C_CVSNOW(row_length,rows)        ! CPL Cnvctv snur rate
      real :: C_RIVEROUT(row_length,rows)      ! CPL Riv outflow->ocn                     

      real :: C_PRESS(row_length,rows)         ! CPL Surf pressure->ocn                     
      real :: C_U10(row_length,rows)         ! CPL 10m wind speed
      real :: C_V10(row_length,rows)         ! CPL 10m wind speed


! UKCA fields are still 1D; change for 7.1
      ! UKCA oxidant fields
      real :: OH_UKCA(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: HO2_UKCA(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: H2O2_UKCA(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      real :: O3_UKCA(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      ! Aerosol climatologies
      real :: ARCLBIOG_BG(row_length,rows,model_levels)
      real :: ARCLBIOM_FR(row_length,rows,model_levels)
      real :: ARCLBIOM_AG(row_length,rows,model_levels)
      real :: ARCLBIOM_IC(row_length,rows,model_levels)
      real :: ARCLBLCK_FR(row_length,rows,model_levels)
      real :: ARCLBLCK_AG(row_length,rows,model_levels)
      real :: ARCLSSLT_FI(row_length,rows,model_levels)
      real :: ARCLSSLT_JT(row_length,rows,model_levels)
      real :: ARCLSULP_AC(row_length,rows,model_levels)
      real :: ARCLSULP_AK(row_length,rows,model_levels)
      real :: ARCLSULP_DI(row_length,rows,model_levels)
      real :: ARCLDUST_B1(row_length,rows,model_levels)
      real :: ARCLDUST_B2(row_length,rows,model_levels)
      real :: ARCLDUST_B3(row_length,rows,model_levels)
      real :: ARCLDUST_B4(row_length,rows,model_levels)
      real :: ARCLDUST_B5(row_length,rows,model_levels)
      real :: ARCLDUST_B6(row_length,rows,model_levels)
      real :: ARCLOCFF_FR(row_length,rows,model_levels)
      real :: ARCLOCFF_AG(row_length,rows,model_levels)
      real :: ARCLOCFF_IC(row_length,rows,model_levels)
      real :: ARCLDLTA_DL(row_length,rows,model_levels)

! End typ_atm_fields.h

  !     Subroutine arguments
  LOGICAL :: PUT_STEP       ! Proper data is to be put
  CHARACTER*(*) :: cmessage ! OUT - Error return message

  CHARACTER(Len=52)   :: Message
  CHARACTER(Len=20)   :: RoutineName
  INTEGER             :: ErrorStat        ! Return code:
  !   0 = Normal exit
  ! +ve = Fatal Error
  ! -ve = Warning

  ErrorStat   = 1
  RoutineName = 'oasis3_puta2o_stub'
  Message     = 'OASIS3 Routines unavailable - see output.'

  WRITE (6,*) '**ERROR**: oasis3_puta2o called but is unavailable.'
  WRITE (6,*) 'Check OASIS3 cpp key is set'

  ! DEPENDS ON: ereport
  CALL ereport(RoutineName, ErrorStat, Message)

END SUBROUTINE oasis3_puta2o
