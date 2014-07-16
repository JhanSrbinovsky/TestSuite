
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: U_MODEL -----------------------------------------------
!LL
!LL  Purpose: High level control program for the Unified Model
!LL           (master routine).  Calls lower level control routines
!LL           according to top level switch settings. Called by
!LL           top level routine UMSHELL which provides dimension sizes
!LL           for dynamic allocation of data arrays.
!LL
!LL  Tested under compiler:   sxmpif90
!LL  Tested under OS version: 
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 8 (01/06/2007)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
      SUBROUTINE U_MODEL(                                               &
     &       NFT,NFTU,                                                  &
! ARGSZSP super arrays lengths not dependent on sub-models
     &  SPD1_LEN,SPSTS_LEN,SPBND_LEN,                                   &
! ARGSZSP end
! ARGSZSPA super arrays lengths (atmosphere)
     &  A_SPDUM_LEN,A_SPPTR_LEN,A_SPCON_LEN,A_SPINF_LEN,A_SPANC_LEN,    &
     &  A_SPBND_LEN,A_SPSTS_LEN,                                        &
! ARGSZSPA end
! 5.5 28/02/03 coupling arrays required for if defined RIVERS C.Bunton
! ARGSZSPC super arrays lengths (atmosphere-ocean coupled)
     &  AO_SPCPL_LEN,                                                   &
! ARGSZSPC end
     & ppxRecs, isSTASH) ! FLUME-STASH

! D1 replacement module
use atm_fields_mod

! OASIS Modules

      USE oasis3_atm_data_mod
      USE auscom_cpl_data_mod, ONLY : l_auscom, access_tfs, ocn_sss


! FLUME-STASH




Use flumerun       
      IMPLICIT NONE

!*L  Interface and arguments: ------------------------------------------
!L       Sizes of super arrays
! TYPSZSP super arrays lengths not dependent on sub-models
      INTEGER :: SPD1_LEN
      INTEGER :: SPSTS_LEN
      INTEGER :: SPBND_LEN
! TYPSZSP end
! TYPSZSPA super arrays lengths (atmosphere)
      INTEGER :: A_SPDUM_LEN
      INTEGER :: A_SPPTR_LEN
      INTEGER :: A_SPCON_LEN
      INTEGER :: A_SPINF_LEN
      INTEGER :: A_SPANC_LEN
      INTEGER :: A_SPBND_LEN
      INTEGER :: A_SPSTS_LEN
! TYPSZSPA end
! TYPSZSPC uper arrays lengths (atmosphere-ocean coupled)
      INTEGER :: AO_SPCPL_LEN
! TYPSZSPC end
!L
!L       Model sizes
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
!L
!L       Addresses of component arrays within super arrays
!---------------------Start of SPINDEX--------------------------------
!L
!L --------------- D1    array  -----------------------------------
!L (now includes extra copies of atmos/ocean D1 for mpp)
!L
      INTEGER     IXD1_LEN       ! No. of arrays
      PARAMETER(  IXD1_LEN =     4)
      INTEGER     IXD1           ! Addresses of arrays
      COMMON/    CIXD1/      IXD1(IXD1_LEN)
!L
!L --------------- Dump headers--------------------------
!L
!L atmos
      INTEGER   A_IXDUM_LEN       ! No. of arrays
      PARAMETER(A_IXDUM_LEN = 14)
      INTEGER   A_IXDUM           ! Addresses of arrays
      COMMON/  CA_IXDUM/ A_IXDUM(A_IXDUM_LEN)
      INTEGER   A_IXSTS_LEN
      PARAMETER(A_IXSTS_LEN = 12)
      INTEGER   A_IXSTS
      COMMON/  CA_IXSTS/ A_IXSTS(A_IXSTS_LEN)
!L
!L --------------- STASH arrays -----------------------------------
!L
      INTEGER     IXSTS_LEN       ! No. of arrays
      PARAMETER(  IXSTS_LEN = 14)
      INTEGER     IXSTS           ! Addresses of arrays
      COMMON/    CIXSTS/      IXSTS(IXSTS_LEN)
!L
!L --------------- Pointers in D1 array and row/level dependent ---
!L --------------- constants
!L atmos
      INTEGER   A_IXPTR_LEN       ! No. of arrays
      PARAMETER(A_IXPTR_LEN = 184)  !CABLE_LAI
      INTEGER   A_IXPTR           ! Addresses of arrays
      COMMON/  CA_IXPTR/ A_IXPTR(A_IXPTR_LEN)
!L
!L --------------- Pre-calculated arrays of constants -------------
!L
!L atmos
      INTEGER   A_IXCON_LEN       ! No. of arrays
! A_IXCON_LEN is so low b/c the bulk of constants that previously were
! stored in um_index.a are now allocated directly in SETCONA 
      PARAMETER(A_IXCON_LEN = 3)
      INTEGER   A_IXCON           ! Addresses of arrays
      COMMON/  CA_IXCON/ A_IXCON(A_IXCON_LEN)
!L
!L --------------- Headers for output interface datasets (boundary
!L                 conditions out)
!L atmos
      INTEGER   A_IXINF_LEN       ! No. of arrays
      PARAMETER(A_IXINF_LEN = 35)
      INTEGER   A_IXINF           ! Addresses of arrays
      COMMON/  CA_IXINF/ A_IXINF(A_IXINF_LEN)
!L
!L --------------- Headers for ancillary files -------------------
!L
!L atmos
      INTEGER   A_IXANC_LEN       ! No. of arrays
      PARAMETER(A_IXANC_LEN = 4)
      INTEGER   A_IXANC           ! Addresses of arrays
      COMMON/  CA_IXANC/ A_IXANC(A_IXANC_LEN)
!L
!L --------------- Headers from input boundary files -------------
!L
!L NOT SUB-MODEL DEPENDENT
      INTEGER     IXBND_LEN       ! No. of arrays
      PARAMETER(  IXBND_LEN = 1)
      INTEGER     IXBND           ! Addresses of arrays
      COMMON/    CIXBND/ IXBND(IXBND_LEN)
!L
!L atmos
      INTEGER   A_IXBND_LEN       ! No. of arrays
      PARAMETER(A_IXBND_LEN = 5)
      INTEGER   A_IXBND           ! Addresses of arrays
      COMMON/  CA_IXBND/ A_IXBND(A_IXBND_LEN)
!L
!L --------------- Constant arrays needed for atmosphere-ocean----
!L --------------- coupling
!L
      INTEGER   AO_IXCPL_LEN      ! No. of arrays
      PARAMETER(AO_IXCPL_LEN = 10)
      INTEGER   AO_IXCPL          ! Addresses of arrays
      COMMON/ CAO_IXCPL/ AO_IXCPL(AO_IXCPL_LEN)
!L
!-------------------End of SPINDEX------------------------------------
!L
      LOGICAL isSTASH   ! true for FLUME parallel-STASH process
!L
!*----------------------------------------------------------------------
!
!  Common blocks
!
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
!*L --------------------- Comdeck: CHISTORY ----------------------------
!LL
!LL  Purpose: COMMON block for history data needed by top level (C0)
!LL           routines, and passed from run to run.  Mostly set by
!LL           the User Interface.
!LL
!LL           Note that CHISTORY *CALLs ALL individual history comdecks
!LL
!LL  Author : A. Sangster
!LL
!LL  Model            Modification history
!LL version  Date
!LL  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!LL                 contents.  RTHBarnes.
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LLEND----------------------------------------------------------------
!*
!CC   *CALL CHSUNITS
! ----------------------- Comdeck: IHISTO   ----------------------------
! Description: COMDECK defining Integer History variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!

      ! Array containing model data time (Same as MODEL_BASIS_TIME/MODEL
      ! ANALYSIS_HRS depending whether before/after assimilation)
      INTEGER :: MODEL_DATA_TIME(6)

      ! Indicator for next mean period to be processed
      INTEGER :: RUN_MEANCTL_RESTART

      ! Indicator of operational run type
      INTEGER :: RUN_INDIC_OP

      ! Final target date for the run
      INTEGER :: RUN_RESUBMIT_TARGET(6)

      ! Last field written/read per FT unit
      INTEGER ::FT_LASTFIELD(20:NUNITS)

! History Common Block for overall model integers variables.

      COMMON /IHISTO/                                                   &
     &  MODEL_DATA_TIME,                                                &
     &  RUN_MEANCTL_RESTART, RUN_INDIC_OP,                              &
     &  RUN_RESUBMIT_TARGET, FT_LASTFIELD

      NAMELIST /NLIHISTO/                                               &
     &  MODEL_DATA_TIME,                                                &
     &  RUN_MEANCTL_RESTART, RUN_INDIC_OP,                              &
     &  RUN_RESUBMIT_TARGET, FT_LASTFIELD

! IHISTO end
! ----------------------- Comdeck: CHISTO   ----------------------------
! Description: COMDECK defining Character History variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.1  18/04/96  Add RUN_IN for qxhistreport.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
      CHARACTER(LEN=10) :: RUN_HIST_TYPE       ! Type of history file
      CHARACTER*8  RUN_TYPE            ! Type of run
      CHARACTER*14 RUN_COMPCODE        ! Run completion code
      CHARACTER*14 RUN_LAST_MEAN       ! Last mean dump created by run
      ! Appears Unused                 ! for pp fields
      CHARACTER*1  RUN_MEANS_TO_DO     ! Flag indicating the run stopped
                                       ! before creating next mean dump
      CHARACTER*1  RUN_OCEAN_FIRST     ! Flag set to true if ocean to be
                                       ! run first
      CHARACTER*8  RUN_JOB_NAME        ! Jobname this run
      CHARACTER*5  RUN_ID              ! Expt./Job id for this run
      CHARACTER*1  RUN_RESUBMIT        ! Flag controlling auto resubmit
      CHARACTER*12 RUN_RESUBMIT_Q      ! Job queue to which resubmit run
      CHARACTER*20 RUN_RESUBMIT_TIME   ! Time at which run resubmits
      CHARACTER*6  RUN_RESUBMIT_CPU    ! Time limit for resubmitted job
      CHARACTER*6  RUN_RESUBMIT_MEMORY ! Resubmitted job's memory limit
      CHARACTER*2  RUN_RESUBMIT_PRTY   ! Resubmitted job intra q prty
      CHARACTER*8  RUN_RESUBMIT_JOBNAME! Resubmitted jobname
      CHARACTER*1  FT_ACTIVE(20:NUNITS) ! "Y" if file partly written

      ! History Common Block for overall model character variables.

      COMMON /CHISTO/                                                   &
     &  RUN_HIST_TYPE, RUN_TYPE, RUN_COMPCODE, RUN_LAST_MEAN,           &
     &  RUN_MEANS_TO_DO, RUN_OCEAN_FIRST, RUN_JOB_NAME, RUN_ID,         &
     &  RUN_RESUBMIT, RUN_RESUBMIT_Q, RUN_RESUBMIT_TIME,                &
     &  RUN_RESUBMIT_CPU, RUN_RESUBMIT_MEMORY, RUN_RESUBMIT_PRTY,       &
     & RUN_RESUBMIT_JOBNAME, FT_ACTIVE

      NAMELIST /NLCHISTO/                                               &
     &  RUN_HIST_TYPE, RUN_TYPE, RUN_COMPCODE, RUN_LAST_MEAN,           &
     &  RUN_MEANS_TO_DO, RUN_OCEAN_FIRST, RUN_JOB_NAME, RUN_ID,         &
     &  RUN_RESUBMIT, RUN_RESUBMIT_Q, RUN_RESUBMIT_TIME,                &
     &  RUN_RESUBMIT_CPU, RUN_RESUBMIT_MEMORY, RUN_RESUBMIT_PRTY,       &
     & RUN_RESUBMIT_JOBNAME, FT_ACTIVE

! CHISTO end
! ----------------------- Comdeck: IHISTG   ----------------------------
! Description: COMDECK defining Integer History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      ! No. of tsteps completed this run
      INTEGER :: LENGTH(N_INTERNAL_MODEL_MAX)

      ! Model end time this run
      INTEGER :: ACTUAL_ENDT(6,N_INTERNAL_MODEL_MAX)

      ! These 2 appears to be purely diagnostic, and not really used.

      ! History block copy of A/O_STEP held in file CTIME
      INTEGER :: H_STEPim(N_INTERNAL_MODEL_MAX)

      ! No of steps in coupling period
      INTEGER :: H_GROUPim(N_INTERNAL_MODEL_MAX)

      ! No of means activated
      INTEGER :: MEAN_OFFSETim(N_INTERNAL_MODEL_MAX)

      ! Offset between MEAN_REFTIME and model basis time(in model dumps)
      INTEGER :: OFFSET_DUMPSim(N_INTERNAL_MODEL_MAX)

      ! No of mean periods chosen
      INTEGER :: MEAN_NUMBERim(N_INTERNAL_MODEL_MAX)

      ! Indicators used to correct logical units are used for
      ! atmos/ocean partial sum dump I/O
      INTEGER :: RUN_MEANCTL_INDICim(4,N_INTERNAL_MODEL_MAX)

      ! History Common Block for generic model integer variables.

      COMMON /IHISTG/                                                   &
     &  H_STEPim, H_GROUPim, MEAN_OFFSETim, OFFSET_DUMPSim,             &
     & MEAN_NUMBERim, RUN_MEANCTL_INDICim

      NAMELIST /NLIHISTG/                                               &
     &  H_STEPim, H_GROUPim, MEAN_OFFSETim, OFFSET_DUMPSim,             &
     & MEAN_NUMBERim, RUN_MEANCTL_INDICim

! IHISTG end
! ----------------------- Comdeck: CHISTG   ----------------------------
! Description: COMDECK defining Character History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.4  30/05/97  Added vars LASTATMim, CURRATMim, LASTDMPim.  K Rogers
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      CHARACTER*14 END_DUMPim(N_INTERNAL_MODEL_MAX)!most recent dumpname
      CHARACTER*80 RESTARTim(N_INTERNAL_MODEL_MAX) !current restart dump
      CHARACTER*14 SAFEDMPim(N_INTERNAL_MODEL_MAX)
! Name of old safe restart dump
      CHARACTER*14 NEWSAFEim(N_INTERNAL_MODEL_MAX)
! Name of new safe restart dump
      CHARACTER*14 LASTATMim(N_INTERNAL_MODEL_MAX) ! Keep name of last
!                                                  ! atmos restart dump
!                                                  ! until ocean dump
      CHARACTER*14 CURRATMim(N_INTERNAL_MODEL_MAX) ! Keep name of
!                                                  ! current atmos
!                                                  ! restart dump
      CHARACTER*14 LASTDMPim(N_INTERNAL_MODEL_MAX) ! Keep name of last
!                                                  ! atmos/ocean dumps
!                                                  ! until meaning done

!
!
! History Common Block for generic model characters variables.
!
      COMMON /CHISTG/                                                   &
     &  END_DUMPim, RESTARTim,                                          &
     &  SAFEDMPim, NEWSAFEim, LASTATMim, CURRATMim, LASTDMPim

      NAMELIST /NLCHISTG/                                               &
     &  END_DUMPim, RESTARTim,                                          &
     &  SAFEDMPim, NEWSAFEim, LASTATMim, CURRATMim, LASTDMPim

! CHISTG end
!*L --------------------- Comdeck: CLFHIST  ----------------------------
!LL
!LL  Purpose: COMDECK defining unit numbers relevant to history file
!LL           and variables used to hold the logical to physical
!LL           file associations made within the model
!LL
!LL  Author : A. Sangster
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LL                  Version 5  18/6/90
!LL
!LL  Model             Modification history from model version 3.0
!LL version  Date
!LL
!LL  3.4  30/09/94  Add files MURKFILE,OUSRANCL,OUSRMULT at 109,113,114
!LL  3.4  05/09/94  Add files USRANCIL,USRMULTI at unit nos. 111,112.
!LL
!LL  3.3  22/11/93  Add file SOURCES at unit number 110. R.T.H.Barnes.
!LL 3.2     28/05/93  Add file BAS_IND at unit number 58. M.Carter.
!LL  Vn3.0  12/02/93 - Variables PERTURB and TRANSP equivalenced to unit
!LL                    numbers 37, and 97 respectively. C.S. Douglas
!LL  3.4  1/8/94     Revised Obs file specification: Stuart Bell
!LL  3.5  01/05/95  Sub-models stage 1: History/control files. RTHBarnes
!    4.0  22/09/95  Added units for Spectral data for Radiation scheme.
!                                        (J. M. Edwards)
!LL  4.1  11/03/96  Introduce Wave sub-model.  RTHBarnes.
!    4.1  26/02/96  Associate new env. variables SO2NATEM and CHEMOXID
!                   with unit nos. 115 & 116. Rename SOURCES to
!                   SULPEMIS. D. Robinson.
!  4.3   18/3/97  Add aerosol forcings of climate change.  Will Ingram
!  4.4   4/7/97   Add ANLINCR  Chris Jones/Stuart Bell
!LL  4.4   12/9/97  Associate ancillary file EVs for initial surface
!LL                 type fracs, initial vegetation state and vegetation
!LL                 disturbance with unit no.s 135-137 R. Betts
!LL  4.4  17/10/97  Associate env var. CACHED with Unit 138. D Robinson
!LL  4.5  22/04/98  Add new ancillary file for soot emissions:
!LL                 SOOTEMIS - in I/O unit 139. R.Rawlins
!LL  4.5  29/07/98  Add new variables ALABCOU5/6/7/8. D. Robinson.
!LL  4.5  17/08/98  Add new variables OLABCOU1/2/3/4. Remove
!LL                 OLABCOUT. D. Robinson.
!LL  5.1  13/04/00  TDF_dump added. ANLINCR changed to IAU_inc.
!LL                 Adam Clayton
!LL  5.2  21/08/00  Add an extra op macro for VAR plus 1 user pp
!LL                 output stream. R Rawlins
!LL  5.2  18/01/01  Add VERT_LEV and attach to unit 90. D. Robinson
!LL  5.3  26/10/01  Add LANDFRAC and attach to unit 120.
!LL                 Free units 75-79 (OBS06-OBS10). D. Robinson
!    5.3  24/10/01  Add IDEALISE and attach to unit 106. A. Malcolm
!    5.3  14/11/00  Added TPPSOZON (tropopause-based ozone). Dave Tan
!    5.4  29/08/02  Add ALABCIN1/2 & attach to unit 125/126. D Robinson
!    5.5  17/02/03  Add Wave model boundary & interface files
!                                                 D.Holmes-Bell
!    5.5  30/12/02  Add DUSTSOIL/BIOMASS and attach to unit 75/76.
!                   RIVSTOR/RIVSEQ/RIVDIR to units 77-79.
!                   D Robinson
!    6.1  07/04/04  Add DMSCONC to unit 95.        A. Jones
!    6.1  08/11/04  Alter names of River Routing files. R.Sharp
!    6.2  26/01/06  Include iteration count file name string
!                                                    T. Edwards
!    6.2 24/00/05   Change STRATOUT to PPSMC and MESOUT to PPSCREEN,
!                   files for soil moisture nudging scheme. Clive Jones
!    6.2  17/03/06  Add SURFEMIS/AIRCREMS/STRATEMS/EXTRAEMS/RADONEMS
!                   filenames to unit nos 130-134 for UKCA emissions
!                   files. WINITIAL/WSTART/WRESTART/WAVANL/WAVANCIN
!                   removed. F. O'Connor
!LL
!LL  Type declarations
!LL
!LL
!LL  Logical Filenames used in the model
!LL
      CHARACTER*80 HKFILE,PPXREF,CONFIG,STASHCTL,NAMELIST,OUTPUT,       &
     &             OUTPUT2,MCTL,ICTL,PHIST,IHIST,THIST,FTXX,            &
     &             CACHE1,CACHE2,ASWAP,OSWAP,AOTRANS,                   &
     &             AINITIAL,ASTART,ARESTART,AOPSUM1,AOPSUM2,AOPSUM3,    &
     &             AOPSUM4,AOMEAN,SSU,                                  &
     &             OZONE,SMCSNOWD,DSOILTMP,SOILTYPE,VEGTYPE,SSTIN,      &
     &             SICEIN,PERTURB,MASK,                                 &
     &             OINITIAL,OSTART,ORESTART,AOPSTMP1,AOPSTMP2,AOPSTMP3, &
     &             AOPSTMP4,                                            &
     &             WFIN,HFLUXIN,PMEIN,ICEFIN,AIRTMP,                    &
     &             SWSPECTD,                                            &
     &             PP0,PP1,PP2,PP3,PP4,PP5,PP6,PP7,PP8,PP9,             &
     &             PPVAR,PP10,                                          &
     &             OBS01,OBS02,OBS03,OBS04,OBS05,                       &
     &             DUSTSOIL,BIOMASS,RIVSTOR,RIVCHAN,RIVER2A,            &
     &             SURFEMIS, AIRCREMS, STRATEMS, EXTRAEMS, RADONEMS,    &
     &             LWSPECTD,WAVEOUT,SURGEOUT,PPSCREEN,PPSMC,WFOUT,      &
     &          HFLUXOUT,FLXCROUT,PMEOUT,ICEFOUT,MOSOUT,SSTOUT,SICEOUT, &
     &             CURNTOUT,DMSCONC,OROG,OLABCIN,OCNDEPTH,CURNTIN,      &
     &             FLUXCORR,SLABHCON,ATMANL,OCNANL,BAS_IND,             &
     &             TRANSP,ATRACER,OTRACER,SULPEMIS,USRANCIL,USRMULTI,   &
     &             OUSRANCL,OUSRMULT,MURKFILE,                          &
     &             ALABCIN1,ALABCIN2,                                   &
     &             ALABCOU1,ALABCOU2,ALABCOU3,ALABCOU4,                 &
     &             ALABCOU5,ALABCOU6,ALABCOU7,ALABCOU8,CARIOLO3,        &
     &             OLABCOU1,OLABCOU2,OLABCOU3,OLABCOU4,                 &
     &             WLABCOU1,WLABCOU2,WLABCOU3,WLABCOU4,HORZGRID,        &
     &             TDF_dump,IAU_inc,                                    &
     &             LANDFRAC,                                            &
     &             SO2NATEM,CHEMOXID,AEROFCG,FRACINIT,VEGINIT,DISTURB,  &
     &             CACHED,SOOTEMIS,                                     &
     &             CO2EMITS,TPPSOZON,                                   &
     &             VERT_LEV,VAR_GRID,                                   &
     &             IDEALISE,ICFILE,                                     &
     &             ARCLBIOG,ARCLBIOM,ARCLBLCK,ARCLSSLT,ARCLSULP,        &
     &             ARCLDUST,ARCLOCFF,ARCLDLTA,RPSEED,OCFFEMIS

!
      CHARACTER*80 MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                               ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        PHIST_UNIT,                                               &
                                 ! Permanent history file unit
     &        IHIST_UNIT,                                               &
                                 ! Interim history file unit
     &        THIST_UNIT,                                               &
                                 ! Temporary history file unit
     &        FTXX_UNIT,                                                &
                                 ! Logical/physical file associations
     &        HKFILE_UNIT        ! Operational houskeeping file unit
!*
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(PHIST_UNIT =10)
      PARAMETER(IHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)
      PARAMETER(FTXX_UNIT  =13)
!
! Namelist of all permissible logical files.
!
      NAMELIST / NLCFILES /                                             &
     &             HKFILE,PPXREF,CONFIG,STASHCTL,NAMELIST,OUTPUT,       &
     &             OUTPUT2,MCTL,ICTL,PHIST,IHIST,THIST,FTXX,            &
     &             CACHE1,CACHE2,ASWAP,OSWAP,AOTRANS,                   &
     &             AINITIAL,ASTART,ARESTART,AOPSUM1,AOPSUM2,AOPSUM3,    &
     &             AOPSUM4,AOMEAN,SSU,                                  &
     &             OZONE,SMCSNOWD,DSOILTMP,SOILTYPE,VEGTYPE,SSTIN,      &
     &             SICEIN,PERTURB,MASK,                                 &
     &             OINITIAL,OSTART,ORESTART,AOPSTMP1,AOPSTMP2,AOPSTMP3, &
     &             AOPSTMP4,                                            &
     &             WFIN,HFLUXIN,PMEIN,ICEFIN,AIRTMP,                    &
     &             SWSPECTD,                                            &
     &             PP0,PP1,PP2,PP3,PP4,PP5,PP6,PP7,PP8,PP9,             &
     &             PPVAR,PP10,                                          &
     &             OBS01,OBS02,OBS03,OBS04,OBS05,                       &
     &             DUSTSOIL,BIOMASS,RIVSTOR,RIVCHAN,RIVER2A,            &
     &             SURFEMIS, AIRCREMS, STRATEMS, EXTRAEMS, RADONEMS,    &
     &             LWSPECTD,WAVEOUT,SURGEOUT,PPSCREEN,PPSMC,WFOUT,      &
     &          HFLUXOUT,FLXCROUT,PMEOUT,ICEFOUT,MOSOUT,SSTOUT,SICEOUT, &
     &             CURNTOUT,DMSCONC,OROG,OLABCIN,OCNDEPTH,CURNTIN,      &
     &             FLUXCORR,SLABHCON,ATMANL,OCNANL,BAS_IND,             &
     &             TRANSP,ATRACER,OTRACER,SULPEMIS,USRANCIL,USRMULTI,   &
     &             OUSRANCL,OUSRMULT,MURKFILE,                          &
     &             ALABCIN1,ALABCIN2,                                   &
     &             ALABCOU1,ALABCOU2,ALABCOU3,ALABCOU4,                 &
     &             ALABCOU5,ALABCOU6,ALABCOU7,ALABCOU8,CARIOLO3,        &
     &             OLABCOU1,OLABCOU2,OLABCOU3,OLABCOU4,                 &
     &             WLABCOU1,WLABCOU2,WLABCOU3,WLABCOU4,HORZGRID,        &
     &             TDF_dump,IAU_inc,                                    &
     &             LANDFRAC,                                            &
     &             SO2NATEM,CHEMOXID,AEROFCG,FRACINIT,VEGINIT,DISTURB,  &
     &             CACHED,SOOTEMIS,                                     &
     &             CO2EMITS,TPPSOZON,                                   &
     &             VERT_LEV,VAR_GRID,                                   &
     &             IDEALISE,ICFILE,                                     &
     &             ARCLBIOG,ARCLBIOM,ARCLBLCK,ARCLSSLT,ARCLSULP,        &
     &             ARCLDUST,ARCLOCFF,ARCLDLTA,RPSEED,OCFFEMIS

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(PHIST      ,MODEL_FT_UNIT(10) ), &
     &(IHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(FTXX      ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &(AOTRANS   ,MODEL_FT_UNIT(17) ),(ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(VEGTYPE   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
     &                                (LWSPECTD   ,MODEL_FT_UNIT(80) ), &
     &(WAVEOUT   ,MODEL_FT_UNIT(81) ),(SURGEOUT   ,MODEL_FT_UNIT(82) ), &
     &(PPSCREEN  ,MODEL_FT_UNIT(83) ),(PPSMC      ,MODEL_FT_UNIT(84) ), &
     &(WFOUT     ,MODEL_FT_UNIT(85) ),(HFLUXOUT   ,MODEL_FT_UNIT(86) ), &
     &(PMEOUT    ,MODEL_FT_UNIT(87) ),(ICEFOUT    ,MODEL_FT_UNIT(88) ), &
     &(MOSOUT    ,MODEL_FT_UNIT(89) ),(VERT_LEV   ,MODEL_FT_UNIT(90) ), &
     &(SSTOUT    ,MODEL_FT_UNIT(91) ),(SICEOUT    ,MODEL_FT_UNIT(92) ), &
     &(CURNTOUT  ,MODEL_FT_UNIT(93) ),(FLXCROUT   ,MODEL_FT_UNIT(94) ), &
     &(DMSCONC   ,MODEL_FT_UNIT(95) ),(OROG       ,MODEL_FT_UNIT(96) ), &
     &(TRANSP    ,MODEL_FT_UNIT(97) ),(OLABCIN    ,MODEL_FT_UNIT(98) ), &
     &(OCNDEPTH  ,MODEL_FT_UNIT(99) ),                                  &
     &(OLABCOU1  ,MODEL_FT_UNIT(100)),(OLABCOU2   ,MODEL_FT_UNIT(101)), &
     &(OLABCOU3  ,MODEL_FT_UNIT(102)),(OLABCOU4   ,MODEL_FT_UNIT(103)), &
     &(IDEALISE  ,MODEL_FT_UNIT(106)),(TDF_dump   ,MODEL_FT_UNIT(107)), &
     &(IAU_inc   ,MODEL_FT_UNIT(108)),(MURKFILE   ,MODEL_FT_UNIT(109)), &
     &(SULPEMIS  ,MODEL_FT_UNIT(110)),(USRANCIL   ,MODEL_FT_UNIT(111)), &
     &(USRMULTI  ,MODEL_FT_UNIT(112)),(OUSRANCL   ,MODEL_FT_UNIT(113)), &
     &(OUSRMULT  ,MODEL_FT_UNIT(114)),(SO2NATEM   ,MODEL_FT_UNIT(115)), &
     &(CHEMOXID  ,MODEL_FT_UNIT(116)),(AEROFCG    ,MODEL_FT_UNIT(117)), &
     &(CO2EMITS  ,MODEL_FT_UNIT(118)),(TPPSOZON   ,MODEL_FT_UNIT(119)), &
     &(LANDFRAC  ,MODEL_FT_UNIT(120)),(WLABCOU1   ,MODEL_FT_UNIT(121)), &
     &(WLABCOU2  ,MODEL_FT_UNIT(122)),(WLABCOU3   ,MODEL_FT_UNIT(123)), &
     &(WLABCOU4  ,MODEL_FT_UNIT(124)),(ALABCIN1   ,MODEL_FT_UNIT(125)), &
     &(ALABCIN2  ,MODEL_FT_UNIT(126)),                                  &
     &(OCFFEMIS  ,MODEL_FT_UNIT(128)),(HORZGRID   ,MODEL_FT_UNIT(129)), &
     &(SURFEMIS  ,MODEL_FT_UNIT(130)),(AIRCREMS   ,MODEL_FT_UNIT(131)), &
     &(STRATEMS  ,MODEL_FT_UNIT(132)),(EXTRAEMS   ,MODEL_FT_UNIT(133)), &
     &(RADONEMS  ,MODEL_FT_UNIT(134)),(FRACINIT   ,MODEL_FT_UNIT(135)), &
     &(VEGINIT   ,MODEL_FT_UNIT(136)),(DISTURB    ,MODEL_FT_UNIT(137)), &
     &(CACHED    ,MODEL_FT_UNIT(138)),(SOOTEMIS   ,MODEL_FT_UNIT(139)), &
     &(ALABCOU1  ,MODEL_FT_UNIT(140)),(ALABCOU2   ,MODEL_FT_UNIT(141)), &
     &(ALABCOU3  ,MODEL_FT_UNIT(142)),(ALABCOU4   ,MODEL_FT_UNIT(143)), &
     &(ALABCOU5  ,MODEL_FT_UNIT(144)),(ALABCOU6   ,MODEL_FT_UNIT(145)), &
     &(ALABCOU7  ,MODEL_FT_UNIT(146)),(ALABCOU8   ,MODEL_FT_UNIT(147)), &
     &(CARIOLO3  ,MODEL_FT_UNIT(148)),(RPSEED     ,MODEL_FT_UNIT(149)), &
     &(PPVAR     ,MODEL_FT_UNIT(150)),(PP10       ,MODEL_FT_UNIT(151)), &
     &(ICFILE    ,MODEL_FT_UNIT(152)),(VAR_GRID   ,MODEL_FT_UNIT(153)), &
     &(ARCLBIOG  ,MODEL_FT_UNIT(154)),(ARCLBIOM   ,MODEL_FT_UNIT(155)), &
     &(ARCLBLCK  ,MODEL_FT_UNIT(156)),(ARCLSSLT   ,MODEL_FT_UNIT(157)), &
     &(ARCLSULP  ,MODEL_FT_UNIT(158)),(ARCLDUST   ,MODEL_FT_UNIT(159)), &
     &(ARCLOCFF  ,MODEL_FT_UNIT(160)),(ARCLDLTA   ,MODEL_FT_UNIT(161))
!
!LL  Comdeck: CCONTROL -------------------------------------------------
!LL
!LL  Purpose: COMMON block for top level switches and 2nd level switches
!LL           needed by the top level (C0) and 2nd level routines, but
!LL           not held in the history COMMON block.
!LL
!LL Logical components covered :
!LL
!LL External documentation: Unified Model documentation paper No
!LL                         Version
!LL
!LLEND ---------------------------------------------------------------

! ----------------------- Comdeck: CNTLALL  ----------------------------
! Description: COMDECK defining Control variables for the
!              model overall.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  16/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.0  25/10/95  Add user switch CONTROL_RESUBMIT. RTHBarnes
!  4.4  28/07/97  Add user switch LCLIMREALYR. M Gallani
!  4.4  11/10/97  Add logical switch L_AO_D1_MEMORY. D. Robinson.
!  5.2  14/11/00  Enable Ocean Run Length Encoding. Ian Edmond
!  5.3  25/09/01  Add switch L_IO_Timer. P.Selwood.
!  5.3  18/09/01  Add FT_LASTSTEP. David Baker
!  5.4  17/09/02  Num_ALBCs and ALBC2_StartTime_mins added. Adam Clayton
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
      ! Array holding original data time (prior to assimilation)
      INTEGER:: MODEL_BASIS_TIME(6)

      ! Model analysis time in hours since Basis Time
      ! UM6.5 - Replace MODEL_ANALYSIS_HRS by MODEL_ANALYSIS_MINS 
      !         MODEL_ANALYSIS_HRS changed to REAL
      REAL :: MODEL_ANALYSIS_HRS
      INTEGER :: MODEL_ANALYSIS_MINS

      INTEGER:: MODEL_HRS_PER_GROUP! No. of hours in coupling period
      INTEGER:: NCPU               ! No of CPUs assigned to the program
      INTEGER:: ANCIL_REFTIME(6)   ! Ref. time for updating ancillaries
      INTEGER:: FT_PLOTSEL(60:69)  ! interval for plotting pp file
      INTEGER:: RUN_TARGET_END(6)  ! Target end time for this run

      INTEGER:: Num_ALBCs            ! Number of atmos boundary files
      INTEGER:: ALBC2_StartTime_mins ! VT of first block of data in 2nd
                                     ! atmos boundary file, in minutes
                                     ! from start of run

      ! Increment to be added on each resubmission of the job.
      INTEGER:: RUN_RESUBMIT_INC(6)

      ! Number of field headers reserved for non-mean PPfiles on each
      ! unit
      INTEGER:: PP_LEN2_LOOK(20:NUNITS)

      ! Internally defined PP packing code
      INTEGER:: PP_PACK_CODE(20:NUNITS)

      ! Frequency of initialisation of FTunit
      INTEGER:: FT_STEPS(20:NUNITS)
      INTEGER :: FT_FIRSTSTEP(20:NUNITS)  ! ... starting at step number .
      INTEGER :: FT_LASTSTEP(20:NUNITS)    ! ... ending at step number ..
      LOGICAL:: LATMOSNEXT,LOCEANNEXT ! Flags to select atmosphere/ocean
      LOGICAL:: LPP                   ! Activate PPCTL
      LOGICAL:: LPP_SELECT(20:NUNITS) ! Activate PP init on unit
      LOGICAL:: LDUMP                 ! Activate DUMPCTL
      LOGICAL:: LMEAN                 ! Activate MEANCTL
      LOGICAL:: LHISTORY              ! Update TEMP history file
      LOGICAL:: LPRINT                ! Activate PRINTCTL
      LOGICAL:: LINTERFACE            ! Activate GEN_INTF
      LOGICAL:: LEXIT                 ! Activate EXITCHEK
      LOGICAL:: LJOBRELEASE           ! Activate JOBCTL
      LOGICAL:: LMEANPR(4)            ! Select printed diags from means
      LOGICAL:: LANCILLARY            ! Activate UP_ANCIL
      LOGICAL:: LBOUNDARY             ! Activate UP_BOUND
      LOGICAL:: LASSIMILATION         ! Activate assimilation
      LOGICAL:: LCAL360               ! 360-day calendar
      LOGICAL:: LTIMER                ! Activate detailed TIMER routine
      LOGICAL:: L_AO_D1_MEMORY  ! T : D1 copied to memory for AO coupling
      LOGICAL:: LCLIMREALYR           ! Real-period climate means
      LOGICAL:: LRLE                  ! Indicates Run Length Encoding
      LOGICAL :: L_IO_TIMER              ! Activate IO Timer.

      CHARACTER*4  EXPT_ID          ! Unique alphanumeric serial number
!                                   ! associated with model
!                                   ! (Non-Operational expts)
!                                   !
!                                   ! Operational run name
!                                   ! (Operational expts)
      CHARACTER*8  EXPT_ALIAS       ! Non unique user defined expt name
      CHARACTER*1  JOB_ID           ! Unique alphanumeric job identifier
!                                   ! used for networking
      CHARACTER*4  EXPT_ID_IN       ! Experiment ID of driving model if
!                                   ! limited-area run
      CHARACTER(LEN=1) :: JOB_ID_IN        ! Job ID of driving model if
!                                   ! limited-area run
      CHARACTER*14 MODEL_STATUS     ! Operational or NonOperational
      CHARACTER*14 MODEL_ASSIM_MODE ! Atmosphere,Ocean,Coupled or None
      CHARACTER*17 TIME_CONVENTION  ! Relative, Timestep, Absolute_long,
!                                    Absolute_standard or Absolute_short
      CHARACTER*1  FT_WSSEND(60:69) ! "Y" if file to be sent to HP
!
      CHARACTER*1 TYPE_LETTER_1(20:NUNITS) ! File type letter #1
      CHARACTER*1 TYPE_LETTER_2(20:NUNITS) ! File type letter #2
      CHARACTER*1 TYPE_LETTER_3(20:NUNITS) ! File type letter #3
!
      CHARACTER*1  FT_INPUT (20:NUNITS) ! "Y" if input file on unit
      CHARACTER*1  FT_OUTPUT(20:NUNITS) ! "Y" if output file on unit
      CHARACTER*1  FT_SELECT(20:NUNITS) ! "Y" if file selected for post
!                                          processing request.
      CHARACTER*1  FT_ARCHSEL(20:NUNITS) ! "Y" if file to be archived.
!
      CHARACTER*10 RUN_ASSIM_MODE      ! cf MODEL_ASSIM_MODE (Oper use)
      CHARACTER*1  CONTROL_RESUBMIT    ! User flag for auto resubmit

      NAMELIST / NLSTCALL /                                             &
     &  MODEL_BASIS_TIME, MODEL_ANALYSIS_MINS,                          &
     &  MODEL_HRS_PER_GROUP,                                            &
     &  NCPU, ANCIL_REFTIME, FT_PLOTSEL, RUN_TARGET_END,                &
     &  Num_ALBCs, ALBC2_StartTime_mins,                                &
     &  RUN_RESUBMIT_INC, PP_LEN2_LOOK, PP_PACK_CODE,                   &
     &  FT_STEPS, FT_FIRSTSTEP, FT_LASTSTEP,                            &
     &  LATMOSNEXT, LOCEANNEXT, LPP, LPP_SELECT, LDUMP, LMEAN,          &
     &  LHISTORY, LPRINT, LINTERFACE, LEXIT, LJOBRELEASE,               &
     &  LMEANPR, LANCILLARY, LBOUNDARY, LASSIMILATION,                  &
     &  LCAL360, LTIMER, L_AO_D1_MEMORY, LRLE,                          &
     &  LCLIMREALYR, L_IO_TIMER,                                        &
     &  EXPT_ID, JOB_ID, EXPT_ID_IN, JOB_ID_IN,                         &
     &  EXPT_ALIAS, MODEL_STATUS, MODEL_ASSIM_MODE,                     &
     &  TIME_CONVENTION, FT_WSSEND,                                     &
     &  TYPE_LETTER_1, TYPE_LETTER_2, TYPE_LETTER_3,                    &
     &  FT_INPUT, FT_OUTPUT, FT_SELECT, FT_ARCHSEL,                     &
     &  RUN_ASSIM_MODE, CONTROL_RESUBMIT
      COMMON / CNTLCALL /                                               &
     &  MODEL_BASIS_TIME, MODEL_ANALYSIS_MINS,                          &
     &  MODEL_HRS_PER_GROUP,                                            &
     &  NCPU, ANCIL_REFTIME, FT_PLOTSEL, RUN_TARGET_END,                &
     &  Num_ALBCs, ALBC2_StartTime_mins,                                &
     &  RUN_RESUBMIT_INC, PP_LEN2_LOOK, PP_PACK_CODE,                   &
     &  FT_STEPS, FT_FIRSTSTEP, FT_LASTSTEP,                            &
     &  LATMOSNEXT, LOCEANNEXT, LPP, LPP_SELECT, LDUMP, LMEAN,          &
     &  LHISTORY, LPRINT, LINTERFACE, LEXIT, LJOBRELEASE,               &
     &  LMEANPR, LANCILLARY, LBOUNDARY, LASSIMILATION,                  &
     &  LCAL360, LTIMER, L_AO_D1_MEMORY, LRLE,                          &
     &  LCLIMREALYR, L_IO_TIMER,                                        &
! Character variables at the end of the common block
     &  EXPT_ID, JOB_ID, EXPT_ID_IN, JOB_ID_IN,                         &
     &  EXPT_ALIAS, MODEL_STATUS, MODEL_ASSIM_MODE,                     &
     &  TIME_CONVENTION, FT_WSSEND,                                     &
     &  TYPE_LETTER_1, TYPE_LETTER_2, TYPE_LETTER_3,                    &
     &  FT_INPUT, FT_OUTPUT, FT_SELECT, FT_ARCHSEL,                     &
     &  RUN_ASSIM_MODE, CONTROL_RESUBMIT
! ----------------------- Comdeck: CNTLGEN  ----------------------------
! Description: COMDECK defining Control variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! Author : R.T.H.Barnes
!
! History:
! Version  Date      Comment.
!  3.5  28/03/95  Sub-Models stage 1: revise History and Control file
!                 contents.  RTHBarnes.
!  4.0   3/11/95  Move character array MEANSim to the end of the
!                 common block to ensure that it starts correctly on a
!                 word boundary. [No problem is apparent on the Cray
!                 if N_INTERNAL_MODEL_MAX is an even no.]
!                 Rick Rawlins
!  4.1  03/04/96  Add new array DUMP_PACKim. D. Robinson
!  4.5  10/11/98  Increase number of dumps allowed at irregular
!                 timesteps from 10 to 40: Move lengths into
!                 CNTLGEN. R Rawlins
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!

      ! Max no. of irregular times for dumps
      INTEGER, PARAMETER :: DUMPTIMES_LEN1 = 40

      ! No. of areas of zonal mean prints
      INTEGER, PARAMETER :: PRINTFREQ_LEN1 = 5

      ! No. of time intervals for climate meaning
      INTEGER, PARAMETER :: MEANFREQ_LEN1 = 4

      ! Max no. of irregular times for job release
      INTEGER, PARAMETER :: JOBREL_LEN1 = 10

      INTEGER:: STEPS_PER_PERIODim(N_INTERNAL_MODEL_MAX)
      INTEGER:: SECS_PER_PERIODim(N_INTERNAL_MODEL_MAX)

      ! Number of advection timesteps between checks for model exit
      INTEGER:: EXITFREQim(N_INTERNAL_MODEL_MAX)

      ! Number of steps between atmosphere restart dumps
      INTEGER:: DUMPFREQim(N_INTERNAL_MODEL_MAX)

      ! Archiving frequency  for atmos dumps
      INTEGER:: ARCHDUMP_FREQim(N_INTERNAL_MODEL_MAX)

      ! Timesteps (from start of run) at which restart dumps are written
      INTEGER:: DUMPTIMESim(DUMPTIMES_LEN1,N_INTERNAL_MODEL_MAX)

      ! Indicators for mean dump frequency
      INTEGER:: MEANFREQim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Switches for mean dump arch.
      INTEGER:: MEANARCHim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! PP field selectors
      INTEGER:: PPSELECTim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Switches for pp field archive
      INTEGER:: ARCHPPSELim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Switches for chart plotting
      INTEGER:: PLOTSELim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Number of field headers to reserve for internal model mean
      ! PPfiles
      INTEGER:: PP_LEN2_MEANim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Reference time for production of means
      INTEGER:: MEAN_REFTIMEim(6,N_INTERNAL_MODEL_MAX)

      ! Indicators of zonal mean print frequency
      INTEGER:: PRINTFREQim(PRINTFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      ! Step numbers  at which to release user-specified scripts

      INTEGER:: JOBREL_STEPim(JOBREL_LEN1,N_INTERNAL_MODEL_MAX)

      ! Offset for dump archiving
      INTEGER:: ARCHDUMP_OFFSETim(N_INTERNAL_MODEL_MAX)

      ! Unit reserved for mean PPs
      INTEGER:: FT_MEANim(N_INTERNAL_MODEL_MAX)

      ! Packing indicator for dumps
      INTEGER:: DUMP_PACKim(N_INTERNAL_MODEL_MAX)

      ! "Y" if mean file to be sent to HP
      CHARACTER(LEN=1) :: MEANWSim(MEANFREQ_LEN1,N_INTERNAL_MODEL_MAX)

      LOGICAL:: LLBOUTim(N_INTERNAL_MODEL_MAX)  ! Lateral b.c.'s
      LOGICAL:: LANCILim(N_INTERNAL_MODEL_MAX)  ! Ancillary files

      NAMELIST / NLSTCGEN /                                             &
     &  STEPS_PER_PERIODim, SECS_PER_PERIODim,                          &
     &  EXITFREQim, DUMPFREQim,                                         &
     &  ARCHDUMP_FREQim, DUMPTIMESim, PPSELECTim, PLOTSELim,            &
     &  ARCHPPSELim, MEANARCHim, MEANFREQim, MEAN_REFTIMEim,            &
     &  PRINTFREQim,  JOBREL_STEPim, ARCHDUMP_OFFSETim, PP_LEN2_MEANim, &
     &  FT_MEANim,                                                      &
     &  DUMP_PACKim,                                                    &
     & MEANWSim, LLBOUTim, LANCILim

      COMMON / CNTLCGEN /                                               &
     &  STEPS_PER_PERIODim, SECS_PER_PERIODim,                          &
     &  EXITFREQim, DUMPFREQim,                                         &
     &  ARCHDUMP_FREQim, DUMPTIMESim, PPSELECTim, PLOTSELim,            &
     &  ARCHPPSELim, MEANARCHim, MEANFREQim, MEAN_REFTIMEim,            &
     &  PRINTFREQim,  JOBREL_STEPim, ARCHDUMP_OFFSETim, PP_LEN2_MEANim, &
     &  FT_MEANim,                                                      &
     &  DUMP_PACKim,                                                    &
     &  LLBOUTim, LANCILim,                                             &
     &  MEANWSim

! CNTLGEN end
! ----------------------- header file: CNTLATM  -----------------------
! Description: Control variables for the Atmosphere model (read only).
!              Contains logical switches for science options needed
!              for addressing calculations and intermediate control.
!              Predominantly used for holding logical flags set up by
!              the UMUI - read in by a namelist, but can also hold
!              derived control flags set by the model.
!              [Note that CRUNTIMC holds accompanying run-time
!              constants.]
!
! Author : R. Rawlins
!
! History:
! Version  Date      Comment.
!  5.0 20/04/99  Extensive revision to earlier comdeck for conversion
!                to C-P 'C' dynamics grid. R. Rawlins
!  5.1  4/02/00  Restore energy correction switches removed at UM5.0
!                Also added additional switches for mass and moisture
!                R A Stratton.
!  5.1 13/04/00  IAU control moved to CTFilt. Adam Clayton
!  5.2 22/08/00  Reinstate Murk and Tracer switches. P.Selwood.
!  5.2 15/11/00  Reintroduce logicals for MOSES 2 and triffid. M.Best.
!  5.3   12/10/01   Remove hard-wired L_trivial_trigs.
!                   put problem_number in CNTLATM.    Terry Davies
!  5.3 27/07/01  Add logical switch L_MURK_RAD   S. Cusack
!  5.3    06/01  Introduce and hardwire logical for setting of
!                leads temperature.                         Nic Gedney
!  5.3 15/10/01  Added L_USE_METHOX switch. David Jackson
!  5.3 29/08/01  Sulphate and sea-salt logicals split into separate
!                versions for different processes.  A. Jones
!  5.3 09/04/01  Add logical switch for spatial degradation of
!                E-S radiation calculations.             S. Cusack
!
!  5.3 19/06/01   Stuff to handle tropopause-based ozone added
!                 -- see Radiation block              Dave Tan
!  5.4 15/08/02   Reconfiguration use changes. P.Selwood.
!  5.4  2/08/02  Add logical switch for PC2 cloud and condensation
!                scheme.                              Damian Wilson
!  5.4 24/10/02  Moved L_GWD and L_USE_USSP from CRUNTIMC. S. Webster
!  5.4 11/03/02  Remove comment lines from after #include
!                                                 S. Carroll
!  5.5 06/02/03  River routing support. P.Selwood.
!  5.5 05/02/03  Add logicals for biomass aerosol scheme     P Davison

!  5.5 17/02/03  Add two large-scale hydrology logicals.
!                L_TOP and L_PDM.                  Nic Gedney
!  5.5 08/01/03  Remove L_3DVAR, L_4DVAR, L_3DVAR_BG, LSINGLE_HYDROL,
!                L_H2_SULPH and L_LSPICE_BDY. D Robinson
!  5.5 21/01/03  Move L_USE_TPPS_OZONE to NLSTCATM namelist. D Robinson
!  5.5 13/03/03  Move I_OZONE_INT here from CRUNTIMC.
!                                                  Jean-Claude Thelen
!  5.5 03/02/03  Include L_mcr logicals.    R.M.Forbes
!  5.5 19/02/03  Remove redundant L_BL_LSPICE and L_RMBL   A.Lock
!  6.0 30/07/03  Include l_pc2_lbc for cloud frac. lbcs. Damian Wilson
!  6.0 19/08/03  Add runtime controls for 4 physics sections
!  6.1  02/08/04  Add logicals for stochem coupling. R Barnes
!  6.1 13/05/04  Add super_array_size variable                 A.Malcolm
!  6.1 07/04/04  Add logical for autoconversion de-biasing.   A. Jones
!  6.1 07/04/04  Add logicals for interactive DMS emissions.  A. Jones
!  6.2 25/01/06  Add iteration count logical test variable    T.Edwards
!  6.2 15/11/05  Add logical for aerosol optical depth     N. Bellouin
!  6.2 23/11/05  Add logical for RH and hygroscopic aerosols N.Bellouin
!  6.2 01/03/06  Add L_UPPER_STRATQ switch - David Jackson
!  6.2 06/01/06  Add logicals for seaice albedo options. J.Ridley
!  6.2 23/02/06  Add logicals for UKCA sub-model.  F.O'Connor
!  6.1 07/04/04  Add logical for RothC temperature function.  C.D. Jones
!  6.2 21/07/05  Add moisture_array_size variable              A.Malcolm
!  6.2 01/10/05  Include L_murk_lbc logical.    R.M.Forbes
!  6.2 07/11/05  Add L_MOD_BARKER_ALBEDO and L_USE_SPEC_SEA
!                to NLSTCATM.        James Manners
!  6.2 25/01/06  Add L_INHOM_CLOUD to NLSTCATM.    James Manners
!  6.2 09/03/06  Add logicals for inland basin rerouting. P.Falloon
!  6.2 24/02/06  Add logical for 10m windspeed pass A2O  J.Gunson

!  6.2 24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!  6.2 07/11/05 Add L_use_orog_corr to NLSTCATM.     James Manners
!  6.2 24/02/06  Add logical to allow DMS flux from ocean model J.Gunson
!  6.4 19/01/07 Removed a comment relating to A05_3c scheme. R A Stratton
!  6.4 08/01/06 Include Brooks cloud fraction logicals. Damian Wilson
!  6.4 10/01/07 Include mixing ratio control logical. Damian Wilson
!  6.6.2 10/06/09 Logicals for reading solar/volcanic forcing. Chris Jones
!------------------   General:  -------------------------------------
      INTEGER Model_domain        ! Domain of atmosphere model:
!                                   global,LAM,cyclic LAM,single column
! Model_domain meaningful names
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

      INTEGER :: problem_number      ! type of problem to be solved

      INTEGER MAXSECTS            ! Max. no. of code sections
      PARAMETER (MAXSECTS=99)
      CHARACTER*3 H_SECT(0:MAXSECTS) ! Array of code section-versions

      ! Physics:   ------------------------------------
      LOGICAL :: l_ssice_albedo     ! Sea-ice albedo affected by snow

      LOGICAL :: l_sice_heatflux    ! Semi-impl update of sea-ice temp

      LOGICAL :: l_sice_meltponds ! Sea-ice albedo affected by meltponds

      LOGICAL :: l_sice_scattering  ! Sea-ice albedo affected scattering

      LOGICAL :: l_sice_hadgem1a ! HadGEM1 sea-ice albedo bug corrected

      LOGICAL :: L_NEG_TSTAR        ! Test for -ve surface temperature.
!
      ! Use sulphate aerosol in microphysics
      LOGICAL :: L_USE_SULPHATE_AUTOCONV

      ! Use sea-salt aerosol in microphysics
      LOGICAL :: L_USE_SEASALT_AUTOCONV

      ! Use soot aerosol in microphysics
      LOGICAL :: L_USE_SOOT_AUTOCONV

      ! Use biomass aerosol in microphysics
      LOGICAL :: L_USE_BMASS_AUTOCONV
      
      ! Use fossil-fuel organic carbon in microphysics
      LOGICAL :: L_USE_OCFF_AUTOCONV
      
      ! Use autoconversion de-biasing scheme in microphysics
      LOGICAL :: L_AUTO_DEBIAS
      ! Use sulphate aerosol no. in S-cycle
      LOGICAL :: L_USE_SULPHATE_SULPC

      ! Use sea-salt aerosol no. in S-cycle
      LOGICAL :: L_USE_SEASALT_SULPC

      ! Use soot aerosol no. in S-cycle
      LOGICAL :: L_USE_SOOT_SULPC

      ! Use biomass aerosol no. in S-cycle
      LOGICAL :: L_USE_BMASS_SULPC
      
      ! Use fossil-organic carbon aerosol no. in S-cycle
      LOGICAL :: L_USE_OCFF_SULPC
      
      ! Energy correction:
      LOGICAL :: L_emcorr    ! T: turns on energy correction code
      LOGICAL :: LMASS_corr  ! T: Apply mass correction
      LOGICAL :: LQT_corr    ! T: Apply total moisture correction
      LOGICAL :: LEMQ_print  ! T: Print additional info from em code
      LOGICAL :: LENERGY     ! T: if timestep to cal energy correction
      LOGICAL :: LFLUX_RESET ! T: if timestep to reset flux array in D1

      ! number of model timesteps per energy correction period.
      INTEGER :: A_ENERGYSTEPS

      ! Radiation:

      LOGICAL :: L_radiation     !  F: Turns off radiation code
      LOGICAL :: L_MICROPHY           !  Microphysics in sw rad scheme

!     Use mineral dust in radiation calculations
      LOGICAL :: L_USE_DUST

!     Use biogenic aerosol in radiation code
      LOGICAL :: L_USE_BIOGENIC

      ! Use SO4 aerosol from sulphur cycle for direct/indirect effect
      ! in radiation, the latter for both SW and LW.
      LOGICAL :: L_USE_SULPC_DIRECT
      LOGICAL :: L_USE_SULPC_INDIRECT_SW
      LOGICAL :: L_USE_SULPC_INDIRECT_LW
      ! Indirect radiative effect of sea-salt
      LOGICAL :: L_USE_SEASALT_INDIRECT

      ! Direct radiative effect of sea-salt
      LOGICAL :: L_USE_SEASALT_DIRECT
      LOGICAL :: L_USE_SOOT_DIRECT  ! direct radiative effects of soot
      LOGICAL :: L_USE_SOOT_INDIRECT  ! indirect effects of soot
      ! Use biomass aerosol for direct/indirect effect in radiation.
      LOGICAL :: L_USE_BMASS_DIRECT
      LOGICAL :: L_USE_BMASS_INDIRECT
      
      ! Use fossil-fuel organic carbon aerosol for direct/indirect
      ! effect in radiation
      LOGICAL :: L_USE_OCFF_DIRECT
      LOGICAL :: L_USE_OCFF_INDIRECT
      
!     Use aerosol climatologies in radiation instead of prognostic variables
!     Set on a species by species basis
      LOGICAL :: L_USE_ARCLBIOM   ! biomass burning aerosol
      LOGICAL :: L_USE_ARCLBLCK   ! black carbon
      LOGICAL :: L_USE_ARCLSSLT   ! sea salt
      LOGICAL :: L_USE_ARCLSULP   ! sulpahtes
      LOGICAL :: L_USE_ARCLDUST   ! mineral dust
      LOGICAL :: L_USE_ARCLOCFF   ! organic carbon (fossil fuel)
      LOGICAL :: L_USE_ARCLDLTA   ! delta aerosol

      ! Aerosol optical depth diagnostic was requested
      LOGICAL :: L_USE_AOD

      ! Clear-sky relative humidity is to be used instead of
      ! grid-box mean RH, for hygroscopic aerosols
      LOGICAL :: L_USE_CLEARRH

      ! controls the use of spatial degradation of radiation calc.
      LOGICAL :: L_rad_deg
      LOGICAL :: L_CTILE       ! Switch for coastal tiling.
!                              ! If TRUE then land and sea can
!                              ! coexist in the same gridbox.
!                              ! If FALSE the land fraction
!                              ! must be equal to 0 to 1
      LOGICAL :: L_TOP         ! If TRUE then TOPMODEL-based
!                              ! hydrology scheme.
      LOGICAL :: L_PDM         ! If TRUE then PDM-based
!                              ! hydrology scheme.
      LOGICAL :: L_SOIL_SAT_DOWN ! If TRUE then super-saturated soil
!                                ! moisture moves downward, else upward

      INTEGER :: H_SWBANDS    ! Number of shortwave radiation bands
      INTEGER :: H_LWBANDS    ! Number of longwave radiation bands
      INTEGER :: A_SW_RADSTEP ! Number of advection steps per SW step
      INTEGER :: A_LW_RADSTEP ! Number of advection steps per LW step
! Number of advection steps per prognostic/diagnostic SW and LW step.
!'Prognostic' and 'Diagnostic' refer to the frequency of the calls
! to radiation code in the Unified Model.
! In the case of time stepping prognostic and diagnostic refer to the
! slow and fast radiative timestep respectively. In the case of radiative
! forcing they refer to the prognostic and diagnostic calls to radiation.
      INTEGER :: A_SW_RADSTEP_DIAG
! Number of advection steps per 'fast' SW step (3C)
      INTEGER :: A_LW_RADSTEP_DIAG
! Number of advection steps per 'fast' LW step (3C)
      INTEGER :: A_SW_RADSTEP_PROG
! Number of advection steps per 'slow' LW step (3C)
      INTEGER :: A_LW_RADSTEP_PROG
! Number of advection steps per 'slow' LW step (3C)

      INTEGER :: i_ozone_int  ! Option for interpolation of ozone

      ! Cloud:

      LOGICAL:: L_CLD_AREA           ! controls cloud area parametriz.
      LOGICAL:: L_ACF_CUSACK         ! ... whether to have Cusack
      LOGICAL:: L_ACF_BROOKS         ! ... or Brooks
      LOGICAL:: L_PC2                ! controls PC2 cloud scheme
      LOGICAL:: L_PC2_RESET          ! run PC2 scheme diagnostically
      LOGICAL:: L_PC2_LBC            ! LBC's contain cloud fractions
      LOGICAL:: L_PC2_DIAG_SH        ! Use diagnostic convective shallow cloud
                                     ! in PC2.

      ! Assimilation:

      ! Switches for assm mode
      LOGICAL:: L_AC

      ! T: Use V_INT_TP to output Temp on model levels.
      LOGICAL:: L_VINT_TP

      ! UM6.5 - MODEL_ANALYSIS_HRS replaced by MODEL_ANALYSIS_MINS - 
      !         change A_ASSIM_START_HR and A_ASSIM_END_HR to 
      !         A_ASSIM_START_MIN, A_ASSIM_END_MIN
      !         so that all three variables  have the same flexibility
      ! Time at which data assimilation starts (Hours after Basis Time)
      INTEGER :: A_ASSIM_START_MIN
      ! Time at which data assimilation  ends
      INTEGER :: A_ASSIM_END_MIN

      ! Switch for assimilation mode
      CHARACTER(LEN=5) :: A_ASSIM_MODE
      
      !---  PMSL diagnostic  ---
      ! Orographic height threshold for new pmsl calculation
      ! from ATMOS_STASH_Misc in UMUI for Calc_NPMSL routine
      REAL :: NPMSL_HEIGHT

      ! Switch for interpolated winds in lbcs
      LOGICAL :: L_int_uvw_lbc  ! .true. for advecting winds
                                ! interpolated in boundary zone
      
      !---  Tracers ---
      ! Aerosol

      LOGICAL:: L_MURK          !      :Total aerosol field
      LOGICAL:: L_MURK_ADVECT   !      :Aerosol advection
      LOGICAL:: L_MURK_SOURCE   !Bndry :Aerosol source & sink terms
      LOGICAL:: L_MURK_BDRY     !Layer :UK Mes bndry model
      LOGICAL:: L_BL_TRACER_MIX !model :Bndry layer tracer mixing
      LOGICAL :: L_MURK_RAD
      LOGICAL :: L_murk_lbc    !  Murk aerosol lbcs active

      ! For Aero_Ctl (Sulphur cycle or Soot)

      INTEGER CALL_CHEM_FREQ     !No. times chem called per atm tstep

!     Mineral dust aerosol
      LOGICAL :: L_DUST
!     Use old version of dust_uplift scheme used in CAM NWP models
      LOGICAL :: L_CAM_DUST

      ! Sulphur cycle

      LOGICAL :: L_SULPC_SO2   ! S Cycle: SO2 MMR included
      LOGICAL :: L_SULPC_DMS   ! S Cycle: DMS MMR included
      LOGICAL :: L_SULPC_OZONE ! S Cycle: Ozone included for oxidation 
                               !          of DMS and SO2
      LOGICAL :: L_SULPC_SO2_O3_NONBUFFERED ! S Cycle: SO2+O3 reaction NOT
                                            ! buffered by NH3.
      LOGICAL :: L_SULPC_ONLINE_OXIDANTS ! Sulphur Cycle : Use online
                                         ! oxidants from UKCA
      LOGICAL :: L_SO2_SURFEM  ! SO2 Surface Emissions
      LOGICAL :: L_SO2_HILEM   ! SO2 High Level Emissions
      LOGICAL :: L_SO2_NATEM   ! SO2 Natural Emissions
      LOGICAL :: L_DMS_EM      ! DMS Emissions
      LOGICAL :: L_DMS_EM_INTER      ! Interactive DMS Emissions
      LOGICAL :: L_DMS_Ointer        ! DMS emissions from ocean model
      LOGICAL :: L_DMS_Liss_Merlivat ! Switches to determine which
      LOGICAL :: L_DMS_Wanninkhof    !    scheme to use for interactive
      LOGICAL :: L_DMS_Nightingale   !    sea-air exchange of DMS
      LOGICAL :: L_SULPC_NH3   ! S Cycle: NH3 tracer included
      LOGICAL :: L_NH3_EM      ! S Cycle: NH3 emiss included

      ! Soot cycle

      LOGICAL :: L_SOOT                ! Soot included
      LOGICAL :: L_SOOT_SUREM          ! surface Soot emiss included
      LOGICAL :: L_SOOT_HILEM          ! elevated Soot emiss included

      ! Biomass aerosol

      LOGICAL :: L_BIOMASS             ! Biomass aerosol included
      LOGICAL :: L_BMASS_SUREM         ! Sfc biomass emiss included
      LOGICAL :: L_BMASS_HILEM         ! Elevated bmass emiss included
      
      ! Fossil-fuel organic carbon aerosol
      
      LOGICAL :: L_OCFF                ! OCFF aerosol included
      LOGICAL :: L_OCFF_SUREM          ! Surface OCFF emissions included
      LOGICAL :: L_OCFF_HILEM          ! Elevated OCFF emiss included
      
      ! Carbon cycle

      ! Interactive 3D CO2 field for use with carbon cycle model
      LOGICAL :: L_CO2_INTERACTIVE
      ! Switch for Radiation Interaction with CO2 - kdcorbin, 06/10
      LOGICAL :: L_CO2_RADIATION   
      ! Switch for CABLE - kdcorbin, 03/10
      LOGICAL :: l_cable
      
      ! Switch for calculating CO2/tracer atmospheric mass - kdcorbin, 05/10
      LOGICAL :: L_TRACER_MASS, L_CO2_MASS
      INTEGER :: I_TRACERMASS_START
      ! Switch for running passive tracers using CO2 fluxes - rml, 1/7/13
      LOGICAL :: L_CO2_TRACER

      ! Switch for calculating methane atmospheric loss - kdcorbin, 05/10
      LOGICAL :: L_METHANE_LOSS
      INTEGER :: I_METHANE_TRACERS

      ! Switch for calculating mcf atmospheric/ocean loss - kdcorbin, 05/10
      LOGICAL :: L_MCF_LOSS
      INTEGER :: I_MCF_TRACERNUMBER

      ! Switch for calculating radon decay - kdcorbin, 05/10
      LOGICAL :: L_RADON_DECAY
      INTEGER :: I_RADON_TRACERNUMBER

      ! CO2 Mass - kdcorbin, 05/10
      REAL :: CO2EMITMASS

     ! Tracer Mass - kdcorbin, 05/10
      REAL :: TMASS(21)

      LOGICAL :: L_CO2_EMITS          ! Include surface emissions

      ! 10m windspeed for air/sea gas flux calculations
      LOGICAL :: L_PSSWND          ! pass 10m windspeed to the ocean

      ! Dust deposition for ocean biology
      LOGICAL :: L_DUST2OCN        ! pass dust dep fields to the ocean

      LOGICAL :: L_Q10                  ! control T fn for soil resp

      ! Switch for turning off boundary layer code

      LOGICAL :: L_bl            !  F: Turns off boundary layer code
      ! MOSES II and Triffid logicals--------------------

      LOGICAL :: L_VEG_FRACS          ! Switch for vegetation fractions
      LOGICAL :: L_SNOW_ALBEDO        ! Prognostic snow albedo
      LOGICAL :: L_TRIFFID            ! Switch for interactive veg model
      LOGICAL :: L_PHENOL             ! Switch for leaf phenology

      ! Switch for running TRIFFID in equilibrium mode
      LOGICAL :: L_TRIF_EQ

      ! Switch for starting NRUN mid-way through a TRIFFID calling
      ! period
      LOGICAL :: L_NRUN_MID_TRIF
      LOGICAL :: L_DISTURB      ! Switch for anthropogenic disturbance

      INTEGER :: CAN_MODEL ! Switch for thermal vegetation canopy

      ! Vegetation:

      ! Update frequency for leaf phenology (days)
      INTEGER :: PHENOL_PERIOD

      INTEGER :: TRIFFID_PERIOD ! Update frequency for TRIFFID (days)

      ! Switch for anthropogenic heat source 
      LOGICAL :: l_anthrop_heat_src 

      ! Hardwire until re-assessment of whether these need to be
      ! re-introduced as UMUI-set switches.
! RR old switches needed for addressing but should be defunct?
! RR - can be set with parameters if needed in the interim. ->

      ! Large scale precipitation:

      LOGICAL :: L_rain          !  F: Turns off precipitation code

      ! 'New' cloud/precip microphysics, Defunct, only mixed phase
      ! phys supported
      LOGICAL, PARAMETER :: L_LSPICE    =.true.

      ! Microphysics complexity
      LOGICAL :: L_mcr_qcf2    !  Include second ice variable
      LOGICAL :: L_mcr_qrain   !  Include prognostic rain
      LOGICAL :: L_mcr_qgraup  !  Include prognosic graupel
      LOGICAL :: L_mcr_qcf2_lbc    !  Second ice variable lbcs active
      LOGICAL :: L_mcr_qrain_lbc   !  Prognostic rain lbcs active
      LOGICAL :: L_mcr_qgraup_lbc  !  Prognosic graupel lbcs active

      ! Controls the use of new RHcrit parametrization, option in Sec 9
      ! vn 2A.
      LOGICAL :: L_RHCPT

! Logicals for different radiation packages
      LOGICAL :: L_FORCING
! Calculate radiative forcings (3C)
      LOGICAL :: L_RADIANCE
! Calculate radiances          (3C)
      LOGICAL :: L_TIMESTEP
! Use new timestepping scheme  (3C)
      LOGICAL :: L_WENYI       
! Include Wenyi's pressure & temperature scaling (3A/3C)

      ! Convection:

      LOGICAL :: L_3D_CCA             ! Use 3D conv cloud amount
      LOGICAL :: L_PHASE_LIM          ! Limits phase change of precip
                                      ! in convective downdraught
      LOGICAL :: L_CCRAD              ! Main logical, will remove code
                                      ! connected with CCRad
                                      ! (including bugfixes)
      LOGICAL :: L_3D_CCW             ! Requires l_ccrad=.TRUE.
                                      ! .TRUE. : Radiation to use 3d ccw
                                      ! profile passed to it from
                                      ! convection.
                                      ! .FALSE.: Radiation constructs
                                      ! mean CCW profile from cclwp,ccb
                                      ! and cct as in original.

      ! Timestep frequency for calling convection
      ! Hardwired to calling every timestep
      INTEGER,PARAMETER :: A_CONV_STEP = 1

      ! GWD scheme:
      LOGICAL :: L_GWD        ! Use SSO drag scheme
      LOGICAL :: L_USE_USSP   ! Use spectral GWD scheme

      ! Radiation:

      ! Changes to open sea albedo for HadGEM1
      LOGICAL :: L_MOD_BARKER_ALBEDO ! Modified Barker albedo
      LOGICAL :: L_USE_SPEC_SEA      ! Spectr. dep. sea albedos

      ! Use modulus of fluxes to remove negative effective extinctions
      LOGICAL :: L_MOD_K_FLUX

      ! Fix the selection of fractional sea points in LW radiation
      LOGICAL :: L_CTILE_FIX

      ! Fix instability in quadratic correction to LW source term
      LOGICAL :: L_QUAD_SRC_FIX

      ! Scale the condensed water content to simulate
      ! inhomogeneous clouds
      LOGICAL :: L_INHOM_CLOUD

! Orography correction to SW radiation
      LOGICAL :: L_use_orog_corr    !  Find gradients from mean orog
      LOGICAL :: L_use_grad_corr    !  Use ancillary X & Y gradients

      ! Tropopause-based Ozone Scheme
      LOGICAL :: L_use_tpps_ozone   !  Use TPPS ozone scheme

! Methane oxidation
      REAL    :: Z_TOP
      LOGICAL :: L_USE_METHOX

! STOCHEM coupling to radiation
      LOGICAL :: L_USE_STOCHEM_CH4   ! for methane
      LOGICAL :: L_USE_STOCHEM_O3    ! for ozone
      ! River Routing
      LOGICAL :: L_RIVERS
      LOGICAL :: L_INLAND   ! control rerouting of inland basin water
      REAL    :: RIVER_STEP

      ! Hydrology:

      LOGICAL :: L_hydrology     !  F: Turns off hydrology code

! Max humidity in STRATQ
      LOGICAL :: L_UPPER_STRATQ

      LOGICAL, PARAMETER :: LMOSES        =.TRUE.  ! MOSES hydrology
      LOGICAL :: L_ICOUNT       !  T: Output iteration counts

      ! Mixing ratios:

      Logical :: l_mr_physics1            ! Use mixing ratio in
                                          ! atmos_physics1
      Logical :: l_mr_physics2            ! Use mixing ratio in
                                          ! atmos_physics2
! Stochastic Physics Random Parameters      
      LOGICAL :: L_RPSEED_READ  !  T: Read in previously specified seed
      LOGICAL :: L_RPSEED_WRITE !  T: WRITE out seed


! Ozone tracer as input to radiation scheme      
      LOGICAL :: L_USE_CARIOLLE
      LOGICAL :: L_USE_OZONEINRAD

! RR old switches---------------------------------------- <-

! OASIS coupling
      LOGICAL :: L_OASIS   ! OASIS coupling switch
      LOGICAL :: L_COUPLE_MASTER    ! Couple through master PE
      INTEGER :: OASIS_COUPLE_FREQ  ! Coupling frequency in
                                    ! number of timesteps. 

!     Logicals for UK Chemistry and Aerosols (UKCA) Model

      LOGICAL :: L_ukca           ! True when UKCA is switched on

! Natural climate forcing
      LOGICAL :: L_SCVARY            ! time varying solar forcing
      LOGICAL :: L_VOLCTS            ! time varying volcanic forcing

      COMMON/ CNTLCATM/                                                 &
     &  Model_domain,L_emcorr,                                          &
     &  L_OASIS, OASIS_COUPLE_FREQ, L_COUPLE_MASTER,                    &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                          &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,              &
     &  problem_number,                                                 &
     &  L_SNOW_ALBEDO,l_ssice_albedo,L_MICROPHY,H_SWBANDS,H_LWBANDS,    &
     &  A_SW_RADSTEP,A_LW_RADSTEP,                                      &
     &  A_SW_RADSTEP_DIAG,A_LW_RADSTEP_DIAG,                            &
     &  A_SW_RADSTEP_PROG,A_LW_RADSTEP_PROG,                            &
     &  L_CLD_AREA,L_ACF_CUSACK,L_ACF_BROOKS,L_PC2,L_PC2_RESET,         &
     &  L_PC2_LBC,L_PC2_diag_sh,L_rad_deg, L_AC ,L_VINT_TP,             &
     &  A_ASSIM_START_MIN, A_ASSIM_END_MIN,                             &
     &  NPMSL_HEIGHT,                                                   &
     &  LENERGY, LFLUX_RESET, LMASS_CORR , LQT_CORR, LEMQ_PRINT,        &
     &  A_ENERGYSTEPS,L_GWD,L_USE_USSP,                                 &
     &  L_Murk, L_MURK_ADVECT, L_MURK_SOURCE, L_MURK_BDRY,              &
     &  L_MURK_RAD, L_murk_lbc, L_BL_TRACER_MIX, L_int_uvw_lbc,         &
     &  L_DUST,L_CAM_DUST,L_SULPC_SO2,L_SULPC_DMS,L_SULPC_OZONE,        &
     &  L_SULPC_SO2_O3_NONBUFFERED,L_SO2_SURFEM, L_SO2_HILEM,           &
     &  L_SO2_NATEM,L_DMS_EM, L_DMS_EM_INTER,                           &
     &  L_SULPC_ONLINE_OXIDANTS,                                        &
     &  L_DMS_Ointer,                                                   &
     &  L_DMS_Liss_Merlivat, L_DMS_Wanninkhof, L_DMS_Nightingale,       &
     &  L_SULPC_NH3, L_NH3_EM, L_SOOT, L_SOOT_SUREM, L_SOOT_HILEM,      &
     &  L_BIOMASS, L_BMASS_SUREM, L_BMASS_HILEM,                        &
     &  L_PSSWND, L_DUST2OCN,                                           &
     &  L_RHCPT, L_CCRAD, L_3D_CCA, L_3D_CCW, L_PHASE_LIM,              &
     &  L_CO2_INTERACTIVE,L_CO2_EMITS,l_cable,                          &
! rml 1/7/13
     &  L_CO2_TRACER,                                                   &
!kdcorbin, 08/10
     &  L_CO2_RADIATION,L_TRACER_MASS,L_CO2_MASS,I_TRACERMASS_START,    &
     &  L_METHANE_LOSS,I_METHANE_TRACERS,L_MCF_LOSS,I_MCF_TRACERNUMBER, &
     &  L_RADON_DECAY,I_RADON_TRACERNUMBER,                             &
     &  L_Q10, L_NEG_TSTAR, L_VEG_FRACS, L_TRIFFID, L_PHENOL,           &
     &  L_TRIF_EQ, L_NRUN_MID_TRIF, L_DISTURB,                          &
     &  CAN_MODEL, PHENOL_PERIOD, TRIFFID_PERIOD,                       &
     &  L_USE_SEASALT_INDIRECT, L_USE_BIOGENIC,                         &
     &  L_USE_SEASALT_DIRECT, L_USE_DUST, L_USE_SULPC_INDIRECT_SW,      &
     &  L_USE_SULPC_INDIRECT_LW, L_USE_SULPC_DIRECT,                    &
     &  L_USE_SULPHATE_AUTOCONV, L_USE_SEASALT_AUTOCONV, L_AUTO_DEBIAS, &
     &  L_USE_SULPHATE_SULPC, L_USE_SEASALT_SULPC,                      &
     &  L_OCFF, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_OCFF_AUTOCONV,        &
     &  L_USE_OCFF_SULPC, L_USE_OCFF_DIRECT, L_USE_OCFF_INDIRECT,       &
     &  L_USE_STOCHEM_CH4, L_USE_STOCHEM_O3,                            &
     &  L_MOD_BARKER_ALBEDO, L_USE_SPEC_SEA, L_MOD_K_FLUX, L_CTILE_FIX, &
     &  L_QUAD_SRC_FIX, L_USE_TPPS_OZONE, I_OZONE_INT,                  &
     &  L_ICOUNT,                                                       &
     &  L_use_orog_corr, L_use_grad_corr,                               &
     &  CALL_CHEM_FREQ, L_USE_SOOT_DIRECT, L_USE_SOOT_INDIRECT,         &
     &  L_USE_SOOT_AUTOCONV, L_USE_SOOT_SULPC, L_USE_BMASS_DIRECT,      &
     &  L_USE_BMASS_INDIRECT, L_USE_BMASS_AUTOCONV, L_USE_BMASS_SULPC,  &
     &  L_USE_ARCLBIOM, L_USE_ARCLBLCK,  L_USE_ARCLSSLT,                &
     &  L_USE_ARCLSULP, L_USE_ARCLDUST,  L_USE_ARCLOCFF, L_USE_ARCLDLTA,&
     &  L_ukca,                                                         &
     &  L_CTILE, L_RIVERS,L_INLAND, RIVER_STEP,                         &
     &  l_sice_meltponds, l_sice_scattering, l_sice_hadgem1a,           &
     &  L_USE_METHOX,Z_TOP,l_mr_physics1,l_mr_physics2,                 &
     &  L_INHOM_CLOUD,                                                  &
     &  L_radiation,L_FORCING,L_TIMESTEP,                               &
     &  L_RADIANCE, L_WENYI, L_bl, L_rain, L_hydrology,                 &
     &  L_TOP,L_PDM,L_USE_AOD,L_USE_CLEARRH,L_UPPER_STRATQ,             &
     &  L_sice_heatflux,L_SOIL_SAT_DOWN,                                &
     &  L_SCVARY,L_VOLCTS,                                              &
     &  l_anthrop_heat_src,                                             &
     &  L_USE_CARIOLLE,L_USE_OZONEINRAD,                                & 
     &  L_RPSEED_READ, L_RPSEED_WRITE,                                  &
     ! Character variables need to be at the end.
     &  H_SECT, A_ASSIM_MODE

      NAMELIST/NLSTCATM/                                                &
     &  Model_domain,L_emcorr,                                          &
     &  L_OASIS, OASIS_COUPLE_FREQ, L_COUPLE_MASTER,                    &
     &  L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,                          &
     &  L_mcr_qcf2_lbc, L_mcr_qrain_lbc, L_mcr_qgraup_lbc,              &
     &  problem_number,                                                 &
     &  L_SNOW_ALBEDO,l_ssice_albedo,L_MICROPHY,H_SWBANDS,H_LWBANDS,    &
     &  A_SW_RADSTEP,A_LW_RADSTEP,                                      &
     &  A_SW_RADSTEP_DIAG,A_LW_RADSTEP_DIAG,                            &
     &  A_SW_RADSTEP_PROG,A_LW_RADSTEP_PROG,                            &
     &  L_CLD_AREA,L_ACF_CUSACK,L_ACF_BROOKS,L_PC2,L_PC2_RESET,         &
     &  L_PC2_LBC,L_PC2_diag_sh, L_rad_deg, L_AC, L_VINT_TP,            &
     &  A_ASSIM_START_MIN, A_ASSIM_END_MIN,                             &
     &  NPMSL_HEIGHT,                                                   &
     &  LMASS_CORR , LQT_CORR, LEMQ_PRINT ,A_ENERGYSTEPS,               &
     &  L_GWD,L_USE_USSP,                                               &
     &  L_Murk, L_MURK_ADVECT, L_MURK_SOURCE, L_MURK_BDRY,              &
     &  L_MURK_RAD, L_murk_lbc, L_BL_TRACER_MIX, L_int_uvw_lbc,         &
     &  L_DUST,L_CAM_DUST,L_SULPC_SO2,L_SULPC_DMS,L_SULPC_OZONE,        &
     &  L_SULPC_SO2_O3_NONBUFFERED, L_SO2_SURFEM, L_SO2_HILEM,          &
     &  L_SO2_NATEM,                                                    &
     &  L_SULPC_ONLINE_OXIDANTS,                                        &
     &  L_DMS_EM, L_DMS_EM_INTER,                                       &
     &  L_DMS_Ointer,                                                   &
     &  L_DMS_Liss_Merlivat, L_DMS_Wanninkhof, L_DMS_Nightingale,       &
     &  L_SULPC_NH3, L_NH3_EM, L_SOOT, L_SOOT_SUREM, L_SOOT_HILEM,      &
     &  L_BIOMASS, L_BMASS_SUREM, L_BMASS_HILEM,                        &
     &  L_PSSWND, L_DUST2OCN,                                           &
     &  L_RHCPT, L_CCRAD, L_3D_CCA, L_3D_CCW, L_PHASE_LIM,              &
     &  L_CO2_INTERACTIVE,L_CO2_EMITS,l_cable,                          &
! rml 1/7/13
     &  L_CO2_TRACER,                                                   &
!kdcorbin, 08/10
     &  L_CO2_RADIATION,L_TRACER_MASS,L_CO2_MASS,I_TRACERMASS_START,    &
     &  L_METHANE_LOSS,I_METHANE_TRACERS,L_MCF_LOSS,I_MCF_TRACERNUMBER, &
     &  L_RADON_DECAY,I_RADON_TRACERNUMBER,                             &  
     &  L_Q10, L_NEG_TSTAR, L_VEG_FRACS, L_TRIFFID, L_PHENOL,           &
     &  L_TRIF_EQ, L_NRUN_MID_TRIF, L_DISTURB,                          &
     &  CAN_MODEL, PHENOL_PERIOD, TRIFFID_PERIOD,                       &
     &  L_USE_SEASALT_INDIRECT, L_USE_BIOGENIC,                         &
     &  L_USE_SEASALT_DIRECT, L_USE_DUST, L_USE_SULPC_INDIRECT_SW,      &
     &  L_USE_SULPC_INDIRECT_LW, L_USE_SULPC_DIRECT,                    &
     &  L_USE_SULPHATE_AUTOCONV, L_USE_SEASALT_AUTOCONV, L_AUTO_DEBIAS, &
     &  L_USE_SULPHATE_SULPC, L_USE_SEASALT_SULPC,                      &
     &  L_OCFF, L_OCFF_SUREM, L_OCFF_HILEM, L_USE_OCFF_AUTOCONV,        &
     &  L_USE_OCFF_SULPC, L_USE_OCFF_DIRECT, L_USE_OCFF_INDIRECT,       &
     &  L_USE_STOCHEM_CH4, L_USE_STOCHEM_O3,                            &
     &  L_MOD_BARKER_ALBEDO, L_USE_SPEC_SEA, L_MOD_K_FLUX, L_CTILE_FIX, &
     &  L_QUAD_SRC_FIX, L_USE_TPPS_OZONE, I_OZONE_INT,                  &
     &  L_ICOUNT,                                                       &
     &  L_use_orog_corr, L_use_grad_corr,                               &
     &  CALL_CHEM_FREQ, L_USE_SOOT_DIRECT, L_USE_SOOT_INDIRECT,         &
     &  L_USE_SOOT_AUTOCONV, L_USE_SOOT_SULPC, L_USE_BMASS_DIRECT,      &
     &  L_USE_BMASS_INDIRECT, L_USE_BMASS_AUTOCONV, L_USE_BMASS_SULPC,  &
     &  L_USE_ARCLBIOM, L_USE_ARCLBLCK,  L_USE_ARCLSSLT,                &
     &  L_USE_ARCLSULP, L_USE_ARCLDUST,  L_USE_ARCLOCFF, L_USE_ARCLDLTA,&
     &  L_ukca,                                                         &
     &  L_CTILE, L_RIVERS,L_INLAND, RIVER_STEP,                         &
     &  l_sice_meltponds, l_sice_scattering, l_sice_hadgem1a,           &
     &  L_USE_METHOX,Z_TOP,l_mr_physics1,l_mr_physics2,                 &
     &  L_INHOM_CLOUD,                                                  &
     &  L_radiation,L_FORCING,L_TIMESTEP,                               &
     &  L_RADIANCE, L_WENYI, L_bl, L_rain, L_hydrology,                 &
     &  L_TOP,L_PDM,L_USE_AOD,L_USE_CLEARRH,L_UPPER_STRATQ,             &
     &  L_SCVARY,L_VOLCTS,                                              &
     &  L_sice_heatflux, L_SOIL_SAT_DOWN,                               &
     &  l_anthrop_heat_src,                                             &
     &  L_USE_CARIOLLE,L_USE_OZONEINRAD,                                &
     &  L_RPSEED_READ, L_RPSEED_WRITE,                                  &
     &  H_SECT, A_ASSIM_MODE

      ! Control switches derived and set within the model, hence not
      ! passed in from the UMUI via namelist:

      ! Radiation
      LOGICAL :: L_SW_RADIATE  ! Activate SW radiation this timestep
      LOGICAL :: L_LW_RADIATE  ! Activate LW radiation this timestep
      LOGICAL :: L_SW_RADIATE_DIAG
! Activate fast SW radiation this timestep (3C)
      LOGICAL :: L_LW_RADIATE_DIAG
! Activate fast LW radiation this timestep (3C)
      LOGICAL :: L_SW_RADIATE_PROG
! Activate slow SW radiation this timestep (3C)
      LOGICAL :: L_LW_RADIATE_PROG
! Activate slow LW radiation this timestep (3C)
      LOGICAL :: Lexpand_ozone ! Convert zonal mean ozone to field

      ! convert zonal mean tpps ozone to field
      LOGICAL :: Lexpand_tpps_ozone

      INTEGER :: I_tpps_ozone_opts ! options for tropopause-based ozone

      ! size of super array holding all tracers
      Integer ::   super_array_size
      Integer ::   moisture_array_size

      COMMON/CNTLCATM2/                                                 &
       !  super tracer array size
     &  super_array_size, moisture_array_size,                          &
     &  L_SW_RADIATE,L_LW_RADIATE,                                      &
     &  L_SW_RADIATE_DIAG,L_LW_RADIATE_DIAG,                            &
     &  L_SW_RADIATE_PROG,L_LW_RADIATE_PROG,                            &
     &  Lexpand_ozone,                                                  &
     & Lexpand_tpps_ozone, I_tpps_ozone_opts
! options for tropopause-based ozone



      LOGICAL :: LTLEADS ! Switch for Leads temperature.
!                      ! If FALSE, they are assumed to be TFS
!                      ! Else they are prognostic.
! Default setting: leads temperatures are set to TFS
! HARDWIRE LTLEADS
      PARAMETER(LTLEADS =.FALSE.)
! ----------------------- Comdeck: CNTLOCN  ----------------------------
! Description: COMDECK defining Control variables for the Ocean
!              internal model.
!   This comdeck contains logical variables which are used on the
!   control of certain sections of Ocean model code
!   They replace the previous method of controlling code using *IF DEFs.
!
! Author : R.T.H.Barnes & R.Hill
!
!    Documentation:  Unified Model Documentation Paper
!                    H- History Bricks
!
!   Type declarations
!
      INTEGER :: O_CLM_START_HR    ! Time ocean climate increments start
      INTEGER :: O_CLM_END_HR      ! Time ocean climate increments end
      INTEGER :: O_INT_CLM_INC     ! # ocean steps  } climate incs.
      INTEGER :: O_INT_ANA_STP     ! # between      } analysis steps

      ! # ocean steps between fwd evolution  of bathys and tesacs
      INTEGER :: O_INT_EVO_BTS

      ! # ocean steps between re-calculation of future bathys and tesacs
      ! valid at this hour

      INTEGER :: O_INT_VRY_BTS

      INTEGER :: O_INT_WTS_ACC    ! # ocean steps betwn accumulating wts

      INTEGER :: O_INT_OBS_FRSH   ! # ocean  } reading new OBS files
      INTEGER :: O_INT_OBS_OUT    ! # steps  } outputting new OBS files
      INTEGER :: O_INT_OBS_STR    ! # between} caching OBS array
      INTEGER :: O_INT_FLD_STR    ! #        } caching model fields

      ! Time at which data assimilation starts (Hours after Basis Time)
      INTEGER :: O_ASSIM_START_HR

      ! Time at which data assimilation ends (Hours after Basis Time)
      INTEGER :: O_ASSIM_END_HR

      INTEGER :: O_ASSIM_ANAL_PER ! Period between analyses (Hours)
      INTEGER :: O_ASSIM_1ST_ANAL    ! First analysis time (Hours)
      LOGICAL :: L_FLUXCORR   ! Heat & water flux correction
      LOGICAL :: L_OGLOBAL    ! Global ocean
      LOGICAL :: L_ICEEVP    ! Use Elastic-Viscous-Plastic Ice dynamics
      LOGICAL :: L_ICESSTILT ! Include sea-surface tilt forcing on ice
      LOGICAL :: L_ICYNPOL   ! North pole fix for use with dynamic ice.
      LOGICAL :: L_ICEFREEDR  ! Free Drift Sea Ice model
      LOGICAL :: L_ICESIMPLE  ! Simple Advection Sea Ice model
      LOGICAL :: L_HADCM4O2I  ! HADCM4 version of ocean-to-ice heat flux
      LOGICAL :: L_IFILTER    ! Filter ice velocities
      LOGICAL :: L_IMCPHEE    ! McPhee ocean-to-ice heat flux
      LOGICAL :: L_IHANEY     ! Haney Forcing Ice
      LOGICAL :: L_ICEITD     ! Use ice thickness distribution
                              ! (i.e. multiple ice categories)
      LOGICAL :: L_ICONCHECK  ! Check ice conservation
                              ! (always FALSE if L_ICEITD=FALSE)
      LOGICAL :: L_ISTRHSM    ! Smooth ice strength
      LOGICAL :: L_ISTRH_PSTAR! If true, use Hibler 79 ice strength
                              ! formula (uses pstar), else use
                              ! Rothrock 75 ice strength formula.
      LOGICAL :: L_ICEFLUXSC  ! Scale A-O coupling of ice related fluxes
                              !    using ice concentrations
      LOGICAL :: L_OADGHR2    ! Ocean assimilation diagnostics
      LOGICAL :: L_OBDY_NORTH   ! Update northern lateral boundary
      LOGICAL :: L_OBDY_SOUTH   ! Update southern lateral boundary
      LOGICAL :: L_OBDY_EAST    ! Update eastern lateral boundary
      LOGICAL :: L_OBDY_WEST    ! Update western lateral boundary
      LOGICAL :: L_OGILL_LBCS   ! Use the Gill boundary scheme
      LOGICAL :: L_OFRS_LBCS    ! Use the FRS boundary scheme
      LOGICAL :: L_OSTVNS_LBCS  ! Use the Stevens boundary scheme
      LOGICAL :: L_OBDY_TRACER  ! Update the tracers
      LOGICAL :: L_OBDY_UV      ! Update the velocities
      LOGICAL :: L_OBDY_STREAM  ! Update the stream functions
      LOGICAL :: L_OBDY_ICE     ! Update ice fields (snow, aice, hice)
!  Start of switches for ocean biogeochemistry and air-sea gas flux
      LOGICAL :: L_OCARBON    ! Carbon cycle model
      ! interactive 3D CO2 field for use with carbon cycle model
      LOGICAL :: L_CO2O_INTERACTIVE
      LOGICAL :: L_OCARB14    ! Calculate atmospheric C12/C14 ratio
      LOGICAL :: L_OCSCEN     ! have a scenario for atmosphere CO2
      LOGICAL :: L_OANCACO2   ! read atmospheric CO2 from ancillary
      LOGICAL :: L_OEXTRAC    ! have two carbon tracers for scenario
      LOGICAL :: L_OCO2OCMIP  ! use carbo chem equm consts from OCMIP
      LOGICAL :: L_OALKSAL    ! base alkalinity on salin (if no bio)
      LOGICAL :: L_OVIRTUAL   ! include virtual surface fluxes
      LOGICAL :: L_OBIOLOGY   ! Effect of phytoplankton on carbon cycle
      LOGICAL :: L_ONUTRELAX  ! relaxation of nutrient to levitus
      LOGICAL :: L_OFEBIO     ! switch for iron-limitation scheme
      LOGICAL :: L_DUST2OCN_O ! get dust deposition from the atmosphere
      LOGICAL :: L_O2BLM4BIO  ! use 2-band light model for bio PrimProd
      LOGICAL :: L_OBIOSTD    ! use standard HadOCC biology
      LOGICAL :: L_OBIODETC   ! separate detrl C,N varbls in std HadOCC
      LOGICAL :: L_OBIODOM    ! use DOM biology
      LOGICAL :: L_OBIODTM    ! switch for Diatom model
      LOGICAL :: L_SWTCHGRZ   ! switch for Fasham switching grazer
      LOGICAL :: L_OFRATIO    ! do ammonium calculations
      LOGICAL :: L_OSRFLX     ! put detritus reaching bottom in surface
      LOGICAL :: L_OBIOTLIM   ! have temperature limitation of phyto.
      LOGICAL :: L_OSHLWCO3   ! even shallow water columns form CaCO3
      LOGICAL :: L_OCTOCHL    ! variable C:Chl for phytoplankton
      LOGICAL :: L_OADVCHL    ! advect chlorophyll as a tracer
      LOGICAL :: L_OCCHLPANC  ! Read carbon:chl ratio from ancillary
      LOGICAL :: L_OOXYGEN    ! Include Oxygen tracer
      LOGICAL :: L_ODMS       ! calculate DMS and flux to the atmosphere
      LOGICAL :: L_OBDIAGS2   ! use alternative bio-diagnostics (DTM)
      LOGICAL :: L_OEXTRACERS ! Include extra tracers
      LOGICAL :: L_OC14TRAC   ! Include Carbon 14 tracer
      LOGICAL :: L_OBOMC14    ! Include bomb Carbon 14 tracer
      LOGICAL :: L_OCFC1112   ! run CFCs as tracers
      LOGICAL :: L_OHELIUM    ! Include Helium-3 and Helium-4 tracers
      LOGICAL :: L_ANCWND_O   ! read 10m windspeed from an ancillary
      LOGICAL :: L_PSSWND_O   ! 10m windspeed passed from atmosphere
      LOGICAL :: L_ICECMSK    ! read an ice mask from an ancillary
      LOGICAL :: L_OCHLANC    ! Read surface chlorophyll from ancillary
      LOGICAL :: L_OLISS      ! Liss & Merlivat wind mixing of tracers
      LOGICAL :: L_OLISS660   ! in Liss/Mer normalise to 660 (old,wrong)
      LOGICAL :: L_OWKHOF     ! Use Wanninkhof 92 piston vel scheme
      LOGICAL :: L_ONGALE     ! Use Nightingale & al piston vel scheme
      LOGICAL :: L_OMNTHWND   ! Use monthly winds in piston vel calc
      LOGICAL :: L_ODLYWND    ! Use daily winds in piston vel calc
!  End of switches for ocean biogeochemistry and air-sea gas flux
      LOGICAL :: L_OCNASSM    ! Activate ocean assimilation
      LOGICAL :: L_OCYCLIC    ! Cyclic boundary conditions
      LOGICAL :: L_OFILTER    ! Fourier filtering for high latitudes
      LOGICAL :: L_OFILTHARD ! Extra stringency on F. filtering
      LOGICAL :: L_OFILTTROP  ! F. filtering on FS barotropic velocities
      LOGICAL :: L_OFILTLBAL  ! Control F. filter load balancing
      LOGICAL :: L_OFREESFC   ! Use free surface conditions
      LOGICAL :: L_OFSMARGINAL   ! Control IFS marginal seas height
      LOGICAL :: L_FLUXD
      LOGICAL :: L_OHANEY     ! Haney Forcing heat/fresh water fluxes
      LOGICAL :: L_OHMEAD     ! Mead tracer transport diagnostics
      LOGICAL :: L_OICECOUP   ! Coupled model with Sea Ice
      LOGICAL :: L_OFLXNRM    ! Flux inputs normalised over sea-ice
      LOGICAL :: L_OIMPDIF    ! CN vertical diffusion scheme
      LOGICAL :: L_OISLANDS   ! Include Island Routines
      LOGICAL :: L_OISOPYC    ! Isopycnal diffusion scheme
      LOGICAL :: L_OLATVISC   ! Latitude dependent viscosity
      LOGICAL :: L_OANIVISC   ! Anisotropic viscosity (as in GloSea)
      LOGICAL :: L_OMIXLAY    ! Wind mixing of tracers-mixed layer scheme
      LOGICAL :: L_ONOCLIN    ! Barotropic solution
      LOGICAL :: L_ONOPOLO    ! No sea ice at North Pole
      LOGICAL :: L_OPENBC     ! Read in lateral boundary fields
      LOGICAL :: L_ORICHARD   ! Evaluate & use Richardson No.
      LOGICAL :: L_OROTATE    ! Coriolis force calculation
      LOGICAL :: L_OSOLAR     ! Calc solar penetration for given water ty
      LOGICAL :: L_OSOLARAL   ! Calc sol. pen. - simplified layer structu
      LOGICAL :: L_OSYMM      ! Symmetric boundary conditions
      LOGICAL :: L_OVARYT     ! Varying time step with depth
      LOGICAL :: L_ORIVERS    ! River run-off routines
      LOGICAL :: L_SEAICE     ! Include Sea Ice model
      LOGICAL :: L_TRANGRID   ! Spatial interp. in coupled model
      LOGICAL :: L_OCONJ     ! Whether to use conjugate gradient solver
      LOGICAL :: L_UPWIND     ! Upwind differencing for tracer advection
      LOGICAL :: L_OPRINT     ! Whether to print incidental ocean info
      LOGICAL :: L_ODELPLUS   !
      LOGICAL :: L_OTROPIC    !
      LOGICAL :: L_OISOMOM
      LOGICAL :: L_OISOGMSKEW
      LOGICAL :: L_OISOGM
      LOGICAL :: L_OBIHARMGM
      LOGICAL :: L_OBIGMCUBE  ! Cubic cos(lat) term for biharm GM
      LOGICAL :: L_OVISHADGEM1
      ! Mediterranean outflow - 288*144 and 96*73 grids only - uses
      ! hardwired gridpoint nos
      LOGICAL :: L_OMEDOUT
      LOGICAL :: L_OCONVROUS  ! Roussenov convective adjustment
      LOGICAL :: L_OEXTRAP ! Extrapolation of vertical density gradients
      LOGICAL :: L_OISOPYCGM  ! Gent and McWilliams eddy parametrisation
      LOGICAL :: L_OISOTAPER  ! Tapering of isopycnal diffusion
      LOGICAL :: L_OVISBECK    ! Visbeck scheme
      LOGICAL :: L_OVISPLUS    ! Enhanced Visbeck for high lat damping
      LOGICAL :: L_OQLARGE     ! Quadratic Large scheme
      LOGICAL :: L_OFULARGE   ! FULL LARGE SCHEME
      LOGICAL :: L_OPANDP     ! RI-DEPENDENT VERT MIX SCHEMES
      LOGICAL :: L_OSTATEC    ! DENSITY CHOICE FOR RI-CALC
      LOGICAL :: L_OUSTARWME  ! WME OR WSTRESS TO FIND USTAR

      LOGICAL :: L_OZVRT      ! barotropic vorticity diagnostic switch
                           ! set by OCN_FOR_STEP (not in namelist)
      LOGICAL :: L_SLOPEMAX   ! Selects SLOPE_MAX isopycnal diffusion
      LOGICAL :: L_COXCNVC    ! Selects original Cox convection scheme
      LOGICAL :: L_OMEDADV
      LOGICAL :: L_OHUDOUT
      LOGICAL :: L_OSTRAIT    ! T => Strait exchange flows parametrised
      LOGICAL :: L_OSTR_CLM   ! T => some strait pts set by climate
                              !      values
      LOGICAL :: L_REFSAL
      LOGICAL :: L_SALFLUXFIX
      LOGICAL :: L_INLANSEA
      LOGICAL :: L_OBOTFRI
      LOGICAL :: L_OEOS25
      LOGICAL :: L_OBIMOM  ! biharmonic momentum diffusion
      LOGICAL :: L_OBISURF  ! biharmonic tracer diffusion in top layers
      LOGICAL :: L_OBMCUBE ! Cubic cos(lat) term for biharm mom diff
      LOGICAL :: L_OBULKRI
      LOGICAL :: L_OWINDMIX
      LOGICAL :: L_OBULKMAXMLD
      LOGICAL :: L_OBDBBL
      LOGICAL :: L_OBIAS
      LOGICAL :: L_ORLP       ! Select rigid lid pressure calculation
      ! Additions to CCONTROL for ocean assimilation

      LOGICAL :: LAS_CLM_INC   ! make increments to relax to climate
      LOGICAL :: LAS_ADD_INC   ! add or subtract analysis increments
      LOGICAL :: LAS_ANA_STP   ! calculate analysis increments
      LOGICAL :: LAS_EVO_BTS   ! evolve bathy and tesac obs 1 step
      LOGICAL :: LAS_VRY_BTS   ! estimate bathys and tesacs at this hour
      LOGICAL :: LAS_WTS_ACC   ! evolve accumulated weights
      LOGICAL :: LAS_OBS_FRSH  ! to refresh main OBS data set
      LOGICAL :: LAS_OBS_OUT   ! output ACOBS file for incremented obs
      LOGICAL :: LAS_FLD_STR   ! output model fields to cache store
      LOGICAL :: LAS_OBS_STR   ! output obs to cache store
      LOGICAL :: L_OMNRLP      ! =T if mean rlp to be read from dump
      LOGICAL :: L_OHADGEM1      ! controls HADGEM1 specific code


      NAMELIST / NLSTCOCN /                                             &
     &  O_CLM_START_HR, O_CLM_END_HR, O_INT_CLM_INC, O_INT_ANA_STP,     &
     &  O_INT_EVO_BTS, O_INT_VRY_BTS, O_INT_WTS_ACC, O_INT_OBS_FRSH,    &
     &  O_INT_OBS_OUT, O_INT_OBS_STR, O_INT_FLD_STR,                    &
     &  O_ASSIM_START_HR, O_ASSIM_END_HR, O_ASSIM_ANAL_PER,             &
     &  O_ASSIM_1ST_ANAL,L_FLUXCORR,L_OGLOBAL,L_ISTRHSM,L_ISTRH_PSTAR,  &
     &  L_IFILTER,L_ICEFLUXSC,                                          &
     &  L_ICEEVP,L_ICESSTILT,L_ICYNPOL,L_ICEFREEDR, L_ICESIMPLE,        &
     &  L_IHANEY, L_ICEITD, L_ICONCHECK,L_HADCM4O2I,L_IMCPHEE,          &
     &  L_OADGHR2,L_OBDY_NORTH,L_OBDY_SOUTH,L_OBDY_EAST,L_OBDY_WEST,    &
     &  L_OGILL_LBCS,L_OFRS_LBCS,L_OSTVNS_LBCS,                         &
     &  L_OBDY_TRACER,L_OBDY_UV,L_OBDY_STREAM,L_OBDY_ICE,               &
!  Start of switches for ocean biogeochemistry and air-sea gas flux
     &  L_OCARBON, L_CO2O_INTERACTIVE, L_OCARB14, L_OCSCEN, L_OANCACO2, &
     &  L_OEXTRAC, L_OCO2OCMIP, L_OALKSAL, L_OVIRTUAL,                  &
     &  L_OBIOLOGY, L_ONUTRELAX, L_OFEBIO, L_DUST2OCN_O, L_O2BLM4BIO,   &
     &  L_OBIOSTD, L_OBIODETC, L_OBIODOM, L_OBIODTM, L_SWTCHGRZ,        &
     &  L_OFRATIO, L_OSRFLX, L_OBIOTLIM, L_OSHLWCO3,                    &
     &  L_OCTOCHL, L_OADVCHL, L_OCCHLPANC,                              &
     &  L_OOXYGEN, L_ODMS, L_OBDIAGS2,                                  &
     &  L_OEXTRACERS, L_OC14TRAC, L_OBOMC14, L_OCFC1112, L_OHELIUM,     &
     &  L_ANCWND_O, L_PSSWND_O, L_ICECMSK, L_OCHLANC,                   &
     &  L_OLISS, L_OLISS660, L_OWKHOF, L_ONGALE, L_OMNTHWND, L_ODLYWND, &
!  End of switches for ocean biogeochemistry and air-sea gas flux
     &  L_OCNASSM, L_OCYCLIC, L_OFILTER, L_OFREESFC, L_OFSMARGINAL,     &
     &  L_FLUXD, L_OFILTLBAL, L_OFILTHARD, L_OFILTTROP,                 &
     &  L_OHANEY, L_OHMEAD, L_OICECOUP, L_OFLXNRM,                      &
     &  L_OIMPDIF, L_OISLANDS, L_OISOPYC, L_OLATVISC, L_OANIVISC,       &
     &  L_OMIXLAY, L_ONOCLIN, L_ONOPOLO, L_OPENBC,                      &
     &  L_ORICHARD, L_OROTATE, L_OSOLAR, L_OSOLARAL,                    &
     &  L_OSYMM, L_OVARYT, L_ORIVERS, L_SEAICE, L_OCONJ,                &
     &  L_TRANGRID, L_UPWIND, L_OPRINT,                                 &
     &  L_ODELPLUS, L_OTROPIC,                                          &
     &  L_OBIGMCUBE,                                                    &
     &  L_OISOMOM,L_OISOGMSKEW,L_OISOGM,L_OBIHARMGM,L_OVISHADGEM1,      &
     &  L_OMEDOUT,                                                      &
     &  L_OCONVROUS,                                                    &
     &  L_OEXTRAP,L_OISOPYCGM,L_OISOTAPER,                              &
     &  L_OVISBECK,                                                     &
     &  L_OVISPLUS,                                                     &
     & L_OBIMOM, L_OBISURF,                                             &
     &  L_OBMCUBE,                                                      &
     &  L_OQLARGE,                                                      &
     &  L_OMEDADV,L_OHUDOUT, L_OSTRAIT, L_OSTR_CLM,                     &
     &  L_REFSAL,L_SALFLUXFIX,L_INLANSEA,L_OBOTFRI,                     &
     &  L_OEOS25,                                                       &
     &  L_OFULARGE,L_OPANDP,L_OSTATEC,L_OUSTARWME,                      &
     &  L_SLOPEMAX,L_COXCNVC,                                           &
     &  L_OBDBBL,                                                       &
     &  L_OBIAS,                                                        &
     &  L_ORLP,                                                         &
      ! additions for control of ocean assimilation
     &  LAS_ADD_INC,LAS_CLM_INC,LAS_ANA_STP,                            &
     &  LAS_EVO_BTS,LAS_VRY_BTS,LAS_WTS_ACC,                            &
     &  LAS_OBS_FRSH,LAS_OBS_OUT,LAS_FLD_STR,LAS_OBS_STR,               &
     &  L_OMNRLP, L_OHADGEM1

      COMMON / CNTLCOCN /                                               &
     &  O_CLM_START_HR, O_CLM_END_HR, O_INT_CLM_INC, O_INT_ANA_STP,     &
     &  O_INT_EVO_BTS, O_INT_VRY_BTS, O_INT_WTS_ACC, O_INT_OBS_FRSH,    &
     &  O_INT_OBS_OUT, O_INT_OBS_STR, O_INT_FLD_STR,                    &
     &  O_ASSIM_START_HR, O_ASSIM_END_HR, O_ASSIM_ANAL_PER,             &
     &  O_ASSIM_1ST_ANAL,L_FLUXCORR,L_OGLOBAL,L_ISTRHSM,L_ISTRH_PSTAR,  &
     &  L_IFILTER,L_ICEFLUXSC,                                          &

     &  L_ICEEVP,L_ICESSTILT,L_ICYNPOL,L_ICEFREEDR, L_ICESIMPLE,        &
     &  L_IHANEY, L_ICEITD, L_ICONCHECK, L_HADCM4O2I,L_IMCPHEE,         &
     &  L_OADGHR2,L_OBDY_NORTH,L_OBDY_SOUTH,L_OBDY_EAST,L_OBDY_WEST,    &
     &  L_OGILL_LBCS,L_OFRS_LBCS,L_OSTVNS_LBCS,                         &
     &  L_OBDY_TRACER,L_OBDY_UV,L_OBDY_STREAM,L_OBDY_ICE,               &
!  Start of switches for ocean biogeochemistry and air-sea gas flux
     &  L_OCARBON, L_CO2O_INTERACTIVE, L_OCARB14, L_OCSCEN, L_OANCACO2, &
     &  L_OEXTRAC, L_OCO2OCMIP, L_OALKSAL, L_OVIRTUAL,                  &
     &  L_OBIOLOGY, L_ONUTRELAX, L_OFEBIO, L_DUST2OCN_O, L_O2BLM4BIO,   &
     &  L_OBIOSTD, L_OBIODETC, L_OBIODOM, L_OBIODTM, L_SWTCHGRZ,        &
     &  L_OFRATIO, L_OSRFLX, L_OBIOTLIM, L_OSHLWCO3,                    &
     &  L_OCTOCHL, L_OADVCHL, L_OCCHLPANC,                              &
     &  L_OOXYGEN, L_ODMS, L_OBDIAGS2,                                  &
     &  L_OEXTRACERS, L_OC14TRAC, L_OBOMC14, L_OCFC1112, L_OHELIUM,     &
     &  L_ANCWND_O, L_PSSWND_O, L_ICECMSK, L_OCHLANC,                   &
     &  L_OLISS, L_OLISS660, L_OWKHOF, L_ONGALE, L_OMNTHWND, L_ODLYWND, &
!  End of switches for ocean biogeochemistry and air-sea gas flux
     &  L_OCNASSM, L_OCYCLIC, L_OFILTER, L_OFREESFC, L_OFSMARGINAL,     &
     &  L_FLUXD, L_OFILTLBAL, L_OFILTHARD, L_OFILTTROP,                 &
     &  L_OHANEY, L_OHMEAD, L_OICECOUP, L_OFLXNRM,                      &
     &  L_OIMPDIF, L_OISLANDS, L_OISOPYC, L_OLATVISC, L_OANIVISC,       &
     &  L_OMIXLAY, L_ONOCLIN, L_ONOPOLO, L_OPENBC,                      &
     &  L_ORICHARD, L_OROTATE, L_OSOLAR, L_OSOLARAL,                    &
     &  L_OSYMM, L_OVARYT, L_ORIVERS, L_SEAICE, L_OCONJ,                &
     &  L_TRANGRID, L_UPWIND, L_OPRINT,                                 &
     &  L_ODELPLUS, L_OTROPIC, L_OZVRT,                                 &
     &  L_OMEDOUT,L_OISOMOM,L_OISOGMSKEW,L_OISOGM,                      &
     &  L_OBIGMCUBE,                                                    &
     &  L_OBIHARMGM,L_OVISHADGEM1,                                      &
     &  L_OCONVROUS,                                                    &
     &  L_OEXTRAP,L_OISOPYCGM,L_OISOTAPER,                              &
     &  L_OVISBECK,                                                     &
     &  L_OVISPLUS,                                                     &
     & L_OBIMOM, L_OBISURF,                                             &
     &  L_OBMCUBE,                                                      &
     &  L_OQLARGE,                                                      &
     &  L_OMEDADV,L_OHUDOUT, L_OSTRAIT, L_OSTR_CLM,                     &
     &  L_REFSAL,L_SALFLUXFIX,L_INLANSEA,L_OBOTFRI,                     &
     &  L_OEOS25,                                                       &
     &  L_OFULARGE,L_OPANDP,L_OSTATEC,L_OUSTARWME,                      &
     &  L_OBULKRI,L_OWINDMIX,L_OBULKMAXMLD,                             &
     &  L_SLOPEMAX,L_COXCNVC,                                           &
     &  L_OBDBBL,                                                       &
     &  L_OBIAS,                                                        &
     &  L_ORLP,                                                         &
      ! additions for control of ocean assimilation
     &  LAS_ADD_INC,LAS_CLM_INC,LAS_ANA_STP,                            &
     &  LAS_EVO_BTS,LAS_VRY_BTS,LAS_WTS_ACC,                            &
     &  LAS_OBS_FRSH,LAS_OBS_OUT,LAS_FLD_STR,LAS_OBS_STR,               &
     &  L_OMNRLP, L_OHADGEM1

! CTIME ----------------------------------------------------
!
!  Purpose: Derived model time/step information including start/end
!           step numbers and frequencies (in steps) of interface field
!           generation, boundary field updating, ancillary field
!           updating; and assimilation start/end times.
!           NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
!           Also contains current time/date information, current
!           step number (echoed in history file) and steps-per-group.
!
!  Model            Modification history from model version 3.0:
! version  Date
!
!   3.1   13/02/93  Dimension arrays A_INTERFACE_STEPS/FSTEP/LSTEP
!                   D. Robinson
!   3.3  01/02/94  Add BASIS_TIME_DAYS to BASIS_TIME_SECS for revised
!                  (32-bit portable) model clock calculations. TCJ
!  3.4  13/12/94  Change COMMOM name from CTIME to CTIMED to satisfy
!                 DEC alpha compiler for portability.  N.Farnon.
!  3.5  12/04/95  Stage 1 submodel changes: move to dimensioning
!                 arrays by internal model. R.Rawlins
!  4.4  06/10/97  Data time of IAU dump added. Adam Clayton.
!  4.5  21/08/98  Remove redundant code. D. Robinson.
!  5.1  13/04/00  Instead of saving full IAU data time, save step on
!                 which data time must be reset during an IAU run.
!                 Adam Clayton
!  5.5  17/02/03  Upgrade Wave model from 4.1 to 5.5 D.Holmes-Bell
!
! Programming standard :
!
!  Logical components covered: C0
!
! Project task :
!
! External documentation: Unified Model documentation paper No:
!                         Version:
!
!END -----------------------------------------------------------------

      INTEGER :: I_YEAR               ! Current model time (years)
      INTEGER :: I_MONTH              ! Current model time (months)
      INTEGER :: I_DAY                ! Current model time (days)
      INTEGER :: I_HOUR               ! Current model time (hours)
      INTEGER :: I_MINUTE             ! Current model time (minutes)
      INTEGER :: I_SECOND             ! Current model time (seconds)
      INTEGER :: I_DAY_NUMBER         ! Current model time (day no)
      INTEGER :: PREVIOUS_TIME(7)     ! Model time at previous step
      INTEGER :: IAU_DTResetStep      ! Data time reset step for IAU run

      INTEGER :: BASIS_TIME_DAYS  ! Integral no of days to basis time
      INTEGER :: BASIS_TIME_SECS  ! No of seconds-in-day at basis time
      INTEGER :: O_CLM_FIRSTSTEP  ! First } step for ocean climate
      INTEGER :: O_CLM_LASTSTEP   ! Last  } increments

      LOGICAL :: L_C360DY

! UM6.5MODEL_ANALYSIS_HRS changed to REAL - 
!   requires FORECAST_HRS and DATA_MINUS_BASIS_HRS to REAL also 
      REAL    :: FORECAST_HRS     ! Hours since Data Time (ie T+nn)
      REAL    :: DATA_MINUS_BASIS_HRS ! Data time - basis time (hours)

      COMMON /CTIMED/ I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,    &
     &  I_DAY_NUMBER,PREVIOUS_TIME,                                     &
     &  BASIS_TIME_DAYS,BASIS_TIME_SECS,                                &
     &  FORECAST_HRS,DATA_MINUS_BASIS_HRS,                              &
     &  IAU_DTResetStep,                                                &
     &  O_CLM_FIRSTSTEP,   O_CLM_LASTSTEP, L_C360DY

      INTEGER :: STEPim(INTERNAL_ID_MAX)  ! Step no since basis time
      INTEGER :: GROUPim(INTERNAL_ID_MAX) ! Number of steps per group

      ! Finish step number this run
      INTEGER :: TARGET_END_STEPim(INTERNAL_ID_MAX)

      REAL :: SECS_PER_STEPim(INTERNAL_ID_MAX) ! Timestep length in secs

      ! Frequency of interface field generation in steps
      INTEGER :: INTERFACE_STEPSim(MAX_N_INTF,INTERNAL_ID_MAX)

      ! Start steps for interface field generation
      INTEGER :: INTERFACE_FSTEPim(MAX_N_INTF,INTERNAL_ID_MAX)

      ! End steps for interface field generation
      INTEGER :: INTERFACE_LSTEPim(MAX_N_INTF,INTERNAL_ID_MAX)

      ! Frequency of  updating boundary fields in steps
      INTEGER :: BOUNDARY_STEPSim(INTERNAL_ID_MAX)

      ! No of steps from boundary data prior to basis time to model
      ! basis time
      INTEGER :: BNDARY_OFFSETim(INTERNAL_ID_MAX)

      ! Lowest frequency for updating of ancillary fields in steps
      INTEGER :: ANCILLARY_STEPSim(INTERNAL_ID_MAX)

      ! Start steps for assimilation
      INTEGER :: ASSIM_FIRSTSTEPim(INTERNAL_ID_MAX)

      ! Number of assimilation steps to analysis
      INTEGER :: ASSIM_STEPSim(INTERNAL_ID_MAX)

      ! Number of assimilation steps after analysis
      INTEGER :: ASSIM_EXTRASTEPSim(INTERNAL_ID_MAX)

      COMMON/CTIMEE/                                                    &
     &  STEPim,GROUPim,TARGET_END_STEPim,INTERFACE_STEPSim,             &
     &  INTERFACE_FSTEPim,INTERFACE_LSTEPim,BOUNDARY_STEPSim,           &
     &  BNDARY_OFFSETim,ANCILLARY_STEPSim,ASSIM_FIRSTSTEPim,            &
     &  ASSIM_STEPSim,ASSIM_EXTRASTEPSim,SECS_PER_STEPim

! CTIME end
!*L------------------ COMDECK CINTFA ----------------------------------
!L CMAXSIZE should be called first.
!
!   Contains Variables, Headers and Index blocks for control of
!   generation of boundary information for the limited area model.
!
!   Interfaces to all other models are handled by STASH, and there is
!   no explicit coding written for them in the model.
!
!Interface variables initialised through INTFCNSTA
!namelist read in the interface control routine INTF_CTL.
!L
!L 29/07/98  CINTF comdeck renamed to CINTFA. New arrays LBC_STREAM_A
!L           and LBC_UNIT_NO_A added. INTF_AK/BK/AKH/BKH removed - now
!L           in ARGINFA/TYPINFA. D. Robinson.
!L 10/11/00  5.2 Add Intf_ExtHalo_NS, Intf_ExtHalo_EW, Intf_RimW_Orog,
!L               LBC_ND, LBC_Z_TOP_MODEL, LBC_BL_LEVELS,
!L               LBC_FIRST_R_RHO and LBC_Q_MIN. D.Robinson
!  22/10/01  5.3 Remove Intf_RimW_Orog. D. Robinson
!L 18/09/01  5.3 Add A_INTF_FREQ_MN,A_INTF_FREQ_SC. Peter Clark
!L
      INTEGER                                                           &
     &  INTF_ROW_LENGTH                                                 &
                         ! Interface field row length
     & ,INTF_P_ROWS                                                     &
                         ! Interface field no of rows
     & ,INTF_P_LEVELS                                                   &
                         ! Interface field no of levels
     & ,INTF_Q_LEVELS                                                   &
                         ! Interface field no of wet levels
     & ,INTF_TR_LEVELS                                                  &
                         ! Interface field no of tracer levels
     & ,INTFWIDTHA                                                      &
                         ! Width of interface zone (atmosphere)
     & ,Intf_ExtHalo_NS                                                 &
                         ! Extended Halo in NS direction
     & ,Intf_ExtHalo_EW                                                 &
                         ! Extended Halo in EW direction
     & ,LBC_ND                                                          &
                         ! LBCs for old UM (0) or ND (1)
     & ,A_INTF_START_HR                                                 &
                         ! ) Start and End time in
     & ,A_INTF_FREQ_HR                                                  &
                         ! ) hours, Frequency in h,m,s for which
     & ,A_INTF_FREQ_MN                                                  &
                         ! ) atmosphere interface data
     & ,A_INTF_FREQ_SC                                                  &
                         ! ) is to be generated.
     & ,A_INTF_END_HR                                                   &
                         ! )
     & ,LEN_INTFA_P                                                     &
                         ! Length of interface p field
     & ,LEN_INTFA_U                                                     &
                         ! Length of interface u field
     & ,LEN_INTFA_DATA                                                  &
                         ! Length of interface data
     & ,INTF_PACK                                                       &
                         ! Packing Indicator for boundary data
     & ,LBC_STREAM_A                                                    &
                         ! Output streams in UMUI
     & ,LBC_UNIT_NO_A                                                   &
                         ! Unit Nos for Atmos Boundary Dataset
     & ,LBC_FIRST_R_RHO                                                 &
                         ! First rho level at which height is constant
     & ,LBC_BL_LEVELS                                                   &
                         ! No of Boundary Layer levels
!
! Following 3 variables not in common ; in namelist
     & ,INTF_METH_LEV_CALC(MAX_N_INTF_A)                                &
!                              !Method of calculating Eta level (ETAK)
!                              !from layers (ETAH)
     & ,INTF_MAX_SIG_HLEV(MAX_N_INTF_A)                                 &
!                              !level below which sigma coordinates used
     & ,INTF_MIN_PRS_HLEV(MAX_N_INTF_A)
!                              !level above which pressure coordinates

      REAL                                                              &
     &  INTF_EWSPACE                                                    &
                         ! E-W grid spacing (degrees)
     & ,INTF_NSSPACE                                                    &
                         ! N-S grid spacing (degrees)
     & ,INTF_FIRSTLAT                                                   &
                         ! Latitude of first row (degrees)
     & ,INTF_FIRSTLONG                                                  &
                         ! Longitude of first row (degrees)
     & ,INTF_POLELAT                                                    &
                         ! Real latitude of coordinate pole (degrees)
     & ,INTF_POLELONG                                                   &
                         ! Real longitude of coordinate pole (degrees)
     & ,LBC_Z_TOP_MODEL                                                 &
                         ! Height of top of model
     & ,LBC_Q_MIN                                                       &
                         ! Minimum value for q
!
! VarRes grid spacing
      , LAMBDA_INTF_P(MAX_INTF_LBCROW_LENGTH, MAX_N_INTF_A)             &
      , LAMBDA_INTF_U(MAX_INTF_LBCROW_LENGTH, MAX_N_INTF_A)             &    
      , PHI_INTF_P(MAX_INTF_LBCROWS, MAX_N_INTF_A)                      &
      , PHI_INTF_V(MAX_INTF_LBCROWS, MAX_N_INTF_A)                      &
      , LAMBDA_LBC_P(MAX_INTF_LBCROW_LENGTH)                            &
      , LAMBDA_LBC_U(MAX_INTF_LBCROW_LENGTH)                            &    
      , PHI_LBC_P(MAX_INTF_LBCROWS)                                     &
      , PHI_LBC_V(MAX_INTF_LBCROWS)

! Following variable not in common ; in namelist
      REAL INTF_ETAH(MAX_INTF_LEVELS+1,MAX_N_INTF_A)
!                          !Eta values at model layer boundaries ETAKH

      LOGICAL                                                           &
     &  INTF_VERT_INTERP                                                &
                         ! Switch to request vertical interpolation
     & ,LNEWBND          ! True for initialising new boundary data file

! Switch for variable resolution LBC output
      LOGICAL  INTF_L_VAR_LBC(MAX_N_INTF_A)

! Files for VERTLEVS namelist     
      Character(Len=80) :: INTF_VERTLEVS

! Files for HorzGrid namelist  
      Character(Len=80) :: INTF_HorzGrid(MAX_N_INTF_A)
!*----------------------------------------------------------------------
      COMMON /INTFCTL_ATMOS/                                            &
     &  INTF_EWSPACE(MAX_N_INTF_A)    ,INTF_NSSPACE(MAX_N_INTF_A)       &
     & ,INTF_FIRSTLAT(MAX_N_INTF_A)   ,INTF_FIRSTLONG(MAX_N_INTF_A)     &
     & ,INTF_POLELAT(MAX_N_INTF_A)    ,INTF_POLELONG(MAX_N_INTF_A)      &
     & ,INTF_ROW_LENGTH(MAX_N_INTF_A) ,INTF_P_ROWS(MAX_N_INTF_A)        &
     & ,INTF_P_LEVELS(MAX_N_INTF_A)   ,INTF_Q_LEVELS(MAX_N_INTF_A)      &
     & ,INTF_TR_LEVELS(MAX_N_INTF_A)  ,INTFWIDTHA(MAX_N_INTF_A)         &
     & ,Intf_ExtHalo_NS(Max_N_Intf_A) ,Intf_ExtHalo_EW(Max_N_Intf_A)    &
     & ,LBC_ND(Max_N_Intf_A)                                            &
     & ,A_INTF_START_HR(MAX_N_INTF_A) ,A_INTF_FREQ_HR(MAX_N_INTF_A)     &
     & ,A_INTF_FREQ_MN(MAX_N_INTF_A)  ,A_INTF_FREQ_SC(MAX_N_INTF_A)     &
     & ,A_INTF_END_HR(MAX_N_INTF_A)   ,LEN_INTFA_P(MAX_N_INTF_A)        &
     & ,LEN_INTFA_U(MAX_N_INTF_A)     ,LEN_INTFA_DATA(MAX_N_INTF_A)     &
     & ,LNEWBND(MAX_N_INTF_A)         ,INTF_VERT_INTERP(MAX_N_INTF_A)   &
     & ,INTF_PACK(MAX_N_INTF_A)       ,LBC_STREAM_A(MAX_N_INTF_A)       &
     & ,LBC_UNIT_NO_A(MAX_N_INTF_A)   ,LBC_FIRST_R_RHO(MAX_N_INTF_A)    &
     & ,LBC_BL_LEVELS(MAX_N_INTF_A)   ,LBC_Z_TOP_MODEL(MAX_N_INTF_A)    &
     & ,INTF_VERTLEVS(MAX_N_INTF_A)   ,LBC_Q_MIN                        &
     & ,INTF_L_VAR_LBC                ,INTF_HORZGRID                    &
     & ,LAMBDA_INTF_P                 ,LAMBDA_INTF_U                    &
     & ,PHI_INTF_P                    ,PHI_INTF_V
!---------------------------------------------------------------------
!
! Description:
!   This comdeck declares an integer variable 'ModelType' whose value
!   determines whether a model run is global, limited area or zonal.
!   The values of ModelType associated with each of the run types are
!   defined by integer parameters which are also declared below.
!    ModelType is set in subroutine SETLOGIC.
!
! Current Code Owner: S.J.Swarbrick
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   3.4    29/9/94  Original code. S.J.Swarbrick
!
      ! Value used to represent the global model configuration
      INTEGER,PARAMETER:: GlobalModel=1

      ! Value used to represent the limited area model configuration
      INTEGER,PARAMETER:: LimitedAreaModel=2

      ! Value used to represent the 'periodic in x' model config
      ! INTEGER,PARAMETER:: ZonalModel=2

! Global scalars:
      INTEGER     ModelType  ! Integer switch which is equated to one
                             ! of the above parameters in a model run,
                             ! and so determines the configuration

! COMMON blocks:
      COMMON /RunType/ ModelType

! C_GLOBAL end
!L
!  Dynamic allocation of ppxref look-up arrays and declaration of
!                                            ppxref pointer array.
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
! Data kind parameters
! Start C_KINDS

! Description:
!   Contains parameters defining kinds for 32 and 64 integers
!   and reals.
!
! Current Code Owner: Paul Selwood
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   6.1  10/03/04   Original code. Paul Selwood

! Parameters for 32 and 64 bit kinds

! Precision and range for 64 bit real
      Integer, Parameter :: prec64  = 15
      Integer, Parameter :: range64 = 307

! Precision and range for 32 bit real
      Integer, Parameter :: prec32  = 6
      Integer, Parameter :: range32 = 37

! Range for integers
      Integer, Parameter :: irange64=15
      Integer, Parameter :: irange32=9

! Kind for 64 bit real
      Integer, Parameter :: real64  = selected_real_kind(prec64,range64)
! Kind for 32 bit real
      Integer, Parameter :: real32  = selected_real_kind(prec32,range32)
! Kind for 64 bit integer
      Integer, Parameter :: integer64 = selected_int_kind(irange64)
! Kind for 32 bit integer
      Integer, Parameter :: integer32 = selected_int_kind(irange32)

! Kinds for 64 and 32 bit logicals. Note that there is no
! "selected_logical_kind", but using the equivalent integer kind is a
! workaround that works on every platform we've tested.
      Integer, Parameter :: logical64 = integer64
      Integer, Parameter :: logical32 = integer32

! End C_KINDS
!
!L
!L  DYNAMIC ALLOCATION OF SUPER ARRAYS:
!L
!L       Main D1 data array
! TYPSPD1 super array of D1 (NOTE: only single component array)
      REAL :: SPD1(SPD1_LEN)
! TYPSPD1 end
!L
!L       STASH related arrays
! TYPSPST super array of STASH arrays
      REAL :: SPSTS(SPSTS_LEN)
! TYPSPST end
!L
!L       Dump headers and lookups
! TYPSPDUA super array of dump arrays (atmosphere)
      REAL :: A_SPDUM(A_SPDUM_LEN)

      ! super array of auxilliary stash arrays (atmos)
      REAL :: A_SPSTS(A_SPSTS_LEN)
! TYPSPDUA end
!L
!L       Pointers (addresses) of model variables and constants
! TYPSPPTA super array of pointers (atmosphere)
      REAL :: A_SPPTR(A_SPPTR_LEN)
! TYPSPPTA end
!L Maximum sizes of fields limited by User Interface
!L CMAXSIZE now included earlier in routine
!L
!L       Model derived constants arrays
! TYPSPCOA super array of constants arrays (atmosphere)
      REAL :: A_SPCON(A_SPCON_LEN)
! TYPSPCOA end
!L
!L       Generation of output interface fields
! TYPSPINA super array of output interface arrays (atmosphere)
      REAL :: A_SPINF(A_SPINF_LEN)
! TYPSPINA end
!L
!L       Updating of model from ancillary files
! TYPSPANA super array of ancillary file arrays (atmosphere)
      REAL :: A_SPANC(A_SPANC_LEN)
! TYPSPANA end
!L
!L       Boundary updating for Limited Area Models
! TYPSPBO  super array of input boundary arrays (not sub-model specific)
      REAL :: SPBND(SPBND_LEN)
! TYPSPBO end
! TYPSPBOA super array of input boundary arrays (atmosphere)
      REAL A_SPBND(A_SPBND_LEN)
! TYPSPBOA end
!L
!L       Coupled model arrays (atmosphere-ocean)
! 5.5 28/02/03 coupling arrays required for if defined RIVERS C.Bunton
! TYPSPCPL super array of coupling arrays (atmosphere-ocean)
      REAL :: AO_SPCPL(AO_SPCPL_LEN)
! TYPSPCPL end
!L
!L       Control for temporal filtering.
! ----------------------- header file: CTFilt  -----------------------
!
! Description:
!
!   Parameters and variables for Incremental Analysis Update (IAU) and
!   Temporal Digital Filtering (TDF) schemes.
!
!
! Current Code Owner: Adam Clayton.
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      INTEGER :: TDF_unit,                                              &
                             ! Unit number for TDF output dump.
     &           IAU_unit    ! Unit number for IAU increment.

      PARAMETER (TDF_unit = 107)
      PARAMETER (IAU_unit = 108)

      INTEGER :: LeaveAlone,                                            &
                             ! \  Options for filtering
     &           DirectFilt,                                            &
                             !  > of advected wind
     &           CopyFiltUVW ! /  components.

      PARAMETER (LeaveAlone  = 1)
      PARAMETER (DirectFilt  = 2)
      PARAMETER (CopyFiltUVW = 3)

      REAL    :: q_CC_tol,                                              &
                             ! qCL/qCF tolerance for cloud clearing.
     &           Weight_LL   ! Lower limit for absolute value of
                             ! filter weight.

      PARAMETER (q_CC_tol  = 1.0E-12  )
      PARAMETER (Weight_LL = 0.0000001)

      REAL    :: q_min       ! Minimum value allowed for q after
                             ! addition of q increments.

      REAL    :: oz_min      ! Minimum value allowed for ozone after
                             ! addition of ozone increments

      INTEGER :: MaxNumWeights ! Maximum number of filter weights.

      PARAMETER (MaxNumWeights = 1000)

      REAL    :: IAU_Weights(MaxNumWeights),                            &
                                             ! IAU filter weights array.
     &           TDF_Weights(MaxNumWeights)  ! TDF filter weights array.

      INTEGER :: IAU_FixHd(LEN_FIXHD) ! Fixed-length header of
                                      ! IAU increment dump.

      INTEGER :: D1_TFilt_len,                                          &
                                 ! D1 length for unpacked fields
     &           D1_IAU_k4_len,                                         &
                                 ! D1 length for packed IAU fields
     &           IAU_Len1Lookup,                                        &
                                 ! \ IAU increment
     &           IAU_Len2Lookup  ! / lookup dimensions

      LOGICAL :: L_Pack_D1_IAU   ! Hold IAU fields in packed form?

      INTEGER :: IAU_NumFldCodes
      PARAMETER (IAU_NumFldCodes = 13)

      ! Codes of fields that may be read in from the IAU file:
      INTEGER :: IAU_FldCodes(IAU_NumFldCodes)
      PARAMETER (IAU_FldCodes = (/ 2,   3,   4,   10,  12, 90,          &
     &                             150, 253, 254, 255, 407, 18001,480 /))

      ! Local lengths of fields (0 if field not required):
      INTEGER :: IAU_LocFldLens(IAU_NumFldCodes)

      ! Field descriptions:
      CHARACTER(7) :: IAU_FldDescs(IAU_NumFldCodes)
      PARAMETER  (IAU_FldDescs = (/ 'u      ', 'v      ', 'theta  ',    &
     &                              'q      ', 'qCF    ', 'aerosol',    &
     &                              'w      ', 'rho    ', 'qCL    ',    &
     &                              'exner  ', 'p      ', 'qT     ',     &
     &                              'ozone  ' /))

  ! IAU namelist variables:
  ! -----------------------

      LOGICAL L_IAU               ! Activate IAU scheme?

      INTEGER IAU_StartMin,                                             &
                                  ! Start minute of filtering period.
     &        IAU_EndMin,                                               &
                                  ! End   minute of filtering period.
     &        IAU_ApexMin         ! Apex minute for triangular filter.

      REAL    IAU_Cutoff_period,                                        &
                                  ! Filter cut-off period in hours.
     &        IAU_SBE_period      ! Stop band edge period in hours.

      LOGICAL L_IAU_CalcExnerIncs ! Use p increments to calculate
                                  ! exner increments?

      LOGICAL L_IAU_CalcThetaIncs ! Calculate theta increments using
                                  ! exner and q increments?

      LOGICAL L_IAU_CalcRhoIncs   ! Calculate rho increments using
                                  ! exner, theta and (possibly) q
                                  ! increments?

      LOGICAL L_IAU_IncTStar      ! If set, add level-one temperature
                                  ! increments to surface temperature
                                  ! and top-level soil temperature.

      LOGICAL L_IAU_ResetPoles    ! If set, reset polar rows of
                                  ! relevant increment fields to
                                  ! their mean values.

      LOGICAL L_IAU_RemoveSS      ! Remove supersaturation wrt water?

      LOGICAL L_IAU_CallStratQ    ! Reset stratospheric humidities at
                                  ! end of IAU insertion period?

      LOGICAL L_IAU_Diags         ! If set, write out IAU diagnostics.

      LOGICAL L_IAU_LowMem        ! If set, activate low memory (but
                                  ! high IO) version of the IAU code.

      LOGICAL L_IAU_DumpTS0State  ! If set, write out model state
                                  ! immediately after timestep-zero
                                  ! call to TFilt_cntl.

      LOGICAL L_IAU_CalcCloudIncs ! If set, calculate q, qcl & Cl
                                  ! from q or qT increments.

      LOGICAL L_IAU_IncrementIce  ! If this and L_IAU_CalcCloudIncs,
                                  ! calculate qcf & Cf

      LOGICAL L_IAU_ScaleCloud    ! If set, scale qcl, qcf, Cl & Cf
                                  ! increments to be in physical bounds
                                  
      LOGICAL L_IAU_UPPER_THETA   ! If set, then constrain the upper theta 
                                  ! increments.

      LOGICAL L_IAU_SetOzoneMin  ! If set, reset ozone to oz_min        
                                 ! in IAU if ozone was negative

      CHARACTER*10                                                      &
     &        IAU_FilterType      ! Filter type.

      REAL    IAU_LL_strat_pv,                                          &
                                  ! Lower-limit for strat ABS(PV).
     &        IAU_UL_strat_p,                                           &
                                  ! Upper-limit for strat pressure.
     &        IAU_LL_trop_p       ! Lower-limit for trop  pressure.

  ! TDF namelist variables:
  ! -----------------------

      LOGICAL L_TDF               ! Activate TDF scheme?

      INTEGER TDF_StartMin,                                             &
                                  ! Start minute of filtering period.
     &        TDF_EndMin,                                               &
                                  ! End   minute of filtering period.
     &        TDF_ApexMin         ! Apex minute for triangular filter.

      REAL    TDF_Cutoff_period,                                        &
                                  ! Filter cut-off period in hours.
     &        TDF_SBE_period      ! Stop band edge period in hours.

      LOGICAL L_TDF_FilterQ,                                            &
                                  ! Filter q?
     &        L_TDF_FilterQCL,                                          &
                                  ! Filter qCL?
     &        L_TDF_FilterQCF     ! Filter qCF?

      LOGICAL L_TDF_ModifyCloud   ! If set, modify cloud variables in
                                  ! TDF dump so that it becomes suitable
                                  ! for starting forecasts including
                                  ! physics.

      INTEGER TDF_AdvWindOpt      ! Filtering option for advected winds.

      LOGICAL L_TDF_CallStratQ    ! Reset stratospheric humidities in
                                  ! TDF dump?

      CHARACTER*10                                                      &
     &        TDF_FilterType      ! Filter type.


      COMMON / CTFilt /                                                 &
     & q_min, oz_min, IAU_Weights, TDF_Weights, IAU_FixHd,              &
     & D1_TFilt_len, D1_IAU_k4_len, IAU_Len1Lookup, IAU_Len2Lookup,     &
     & L_Pack_D1_IAU, IAU_LocFldLens,                                   &
     & L_IAU, IAU_StartMin, IAU_EndMin,                                 &
     & IAU_ApexMin, IAU_Cutoff_period, IAU_SBE_period,                  &
     & L_IAU_CalcCloudIncs, L_IAU_IncrementIce, L_IAU_ScaleCloud,       &
     & L_IAU_CalcExnerIncs, L_IAU_CalcThetaIncs, L_IAU_CalcRhoIncs,     &
     & L_IAU_IncTStar, L_IAU_ResetPoles, L_IAU_RemoveSS,                &
     & L_IAU_CallStratQ, L_IAU_Diags, L_IAU_LowMem, L_IAU_DumpTS0State, &
     & L_IAU_UPPER_THETA,                                               &
     & IAU_LL_strat_pv, IAU_UL_strat_p, IAU_LL_trop_p,                  &
     & L_TDF, TDF_StartMin, TDF_EndMin,                                 &
     & TDF_ApexMin, TDF_Cutoff_period, TDF_SBE_period,                  &
     & L_TDF_FilterQ, L_TDF_FilterQCL, L_TDF_FilterQCF,                 &
     & L_TDF_ModifyCloud, TDF_AdvWindOpt,                               &
     & L_TDF_CallStratQ, L_IAU_SetOzoneMin,                             &
     ! Character variables at the end of common block
     & IAU_FilterType, TDF_FilterType

      NAMELIST / RUN_TFilt /                                            &
     & L_IAU, IAU_StartMin, IAU_EndMin,                                 &
     & IAU_ApexMin, IAU_Cutoff_period, IAU_SBE_period,                  &
     & L_IAU_CalcCloudIncs, L_IAU_IncrementIce, L_IAU_ScaleCloud,       &
     & L_IAU_CalcExnerIncs, L_IAU_CalcThetaIncs, L_IAU_CalcRhoIncs,     &
     & L_IAU_IncTStar, L_IAU_ResetPoles, L_IAU_RemoveSS,                &
     & L_IAU_CallStratQ, L_IAU_Diags, L_IAU_LowMem, L_IAU_DumpTS0State, &
     & L_IAU_UPPER_THETA,                                               &
     & IAU_LL_strat_pv, IAU_UL_strat_p, IAU_LL_trop_p,                  &
     & L_TDF, TDF_StartMin, TDF_EndMin,                                 &
     & TDF_ApexMin, TDF_Cutoff_period, TDF_SBE_period,                  &
     & L_TDF_FilterQ, L_TDF_FilterQCL, L_TDF_FilterQCF,                 &
     & L_TDF_ModifyCloud, TDF_AdvWindOpt,                               &
     & L_TDF_CallStratQ, L_IAU_SetOzoneMin,                             &
     & IAU_FilterType, TDF_FilterType
!  Sizes for allocation of AC scheme arrays
! CSIZEOBS start
! The variables involved in dimensioning observational data arrays
! in the data assimilation section are stored in this comdeck.

      ! For ATMOSPHERE assimilations the values are computed in the
      ! initialisation routine INITAC and then passed to the main
      ! routine AC.

      INTEGER :: A_MAX_NO_OBS    !  No of observations in AC Obs files.
      INTEGER :: A_MAX_OBS_SIZE  !  No of obs values in AC Obs files.

      COMMON /CSIZEOBS/ A_MAX_NO_OBS, A_MAX_OBS_SIZE
! CSIZEOBS end
      INTEGER :: obs_flag_len,obs_len
      INTEGER, ALLOCATABLE, DIMENSION(:) :: OBS_FLAG
      REAL,    ALLOCATABLE, DIMENSION(:) :: OBS
!L       Control for boundary updating.
! CBOUND start

!  History:
!  Date      Vn     Modification
!  31/10/01  5.3    Remove RIMWEIGHTS_OROG. D. Robinson
!  17/09/02  5.4    Variables for controlling 2nd bndy file. A Clayton
!  17/02/03  5.5    Allow Wave model to use boundary code. D.Holmes-Bell
!  20/01/06  6.2    Add Current_LBC_Step. Dave Robinson
!  01/03/06  6.2    Remove RIMWEIGHTSW. Dave Robinson

      ! These 3 arrays are set by namelist read in IN_BOUND hence
      ! cannot be in argument list nor in COMMON if array lengths are
      ! passed variables. Only way seems to be to set MAX allowed
      ! lengths consistent with User Interface so that can be in COMMON
      INTEGER, PARAMETER :: MAX_BND_FLDS=4
      INTEGER, PARAMETER :: MAX_RIMWIDTH=10
      INTEGER :: BOUND_FIELDCODE(MAX_BND_FLDS)  ! Set by NAMELIST
      REAL :: RIMWEIGHTSA(MAX_RIMWIDTH)      ! Set by NAMELIST
      REAL :: RIMWEIGHTSO(MAX_RIMWIDTH)      ! Set by NAMELIST

      ! Variable for controlling LBC updating
      ! - initialised in INBOUNDA and updated in BOUNDVAL

      Integer  :: Current_LBC_Step     ! Timestep at which LBCs were
                                       ! last updated.

      ! Variables for controlling 2nd atmos boundary file. All are
      ! calculated within the code from other data.
      INTEGER :: ALBC_num              ! Number of atmos boundary file
                                       ! currently in use
      INTEGER :: ALBC2_StartTime_steps ! VT of first block of data in
                                       ! 2nd atmos boundary file, in
                                       ! steps from start of run
      INTEGER :: ALBC_SwapStep         ! Step on which to swap to 2nd
                                       ! atmos boundary file

      COMMON /BOUND_CT/ BOUND_FIELDCODE,                                &
     &                  RIMWEIGHTSA, RIMWEIGHTSO,                       &
     &                  ALBC_num, ALBC2_StartTime_steps, ALBC_SwapStep, &
     &                  Current_LBC_Step

! CBOUND end
!L
!
!  Local variables
!
      INTEGER internal_model    ! Work - Internal model identifier
      INTEGER internal_model_prev!Work - Previous internal model ident
      INTEGER submodel          ! Work - Submodel id for dump partition
      INTEGER submodel_prev     ! Work - Previous submodel dump id
      INTEGER NGROUP            ! Work - Number of steps in "group"
      INTEGER MEANLEV           ! Work - Mean level indicator
      INTEGER IABORT            ! Work - Internal return code
      INTEGER I_STEP            ! Work - Loop counter over timesteps
      INTEGER G_THETA_FIELD_SIZE                                        &
                                   ! Sizes for mpp dynamic allocation
     &       ,G_IMTJMT          ! in A-O coupling routines
!
! River routing
      INTEGER G_RIVER_FIELD_SIZE   ! Sizes for mpp dynamic allocation
!
      LOGICAL LEXITNOW          ! Work - Immediate exit indicator
      CHARACTER*14 PPNAME       ! Work - Dummy PP filename
      INTEGER NFT           ! Unit no. for standard STASHmaster files
      INTEGER NFTU          ! Do. user STASH files (for GET_FILE)
      INTEGER RowNumber     ! Row no. counter for PPXI, PPXC arrays
      INTEGER I,J,K         ! Loop counters
      INTEGER CO2_DIMA,                                                 &
                                   ! CO2 array dimensions
     &        CO2_DIMO,                                                 &
     &        CO2_DIMO2
      INTEGER DMS_DIMA,                                                 &
                                   ! DMS array dimensions
     &        DMS_DIMO,                                                 &
     &        DMS_DIMO2
      Integer info   ! Return code from GCom routines

      integer len_runid       !  No of chars in RUNID

      character*4 runtype_char  !  Run Type (ie. NRUN, CRUN)

! 3-D fields of species to be passed down to radiation
      INTEGER, PARAMETER :: ngrgas = 8
      INTEGER, SAVE :: grgas_addr(ngrgas)

      ! Lookup tables of IAU increment file:
      INTEGER IAU_lookup(IAU_Len1Lookup, IAU_Len2Lookup)

      ! Arrays for storing IAU/TDF fields, one for packed fields, one
      ! for unpacked fields:
      REAL(KIND=real32) D1_IAU_k4(D1_IAU_k4_len)
      REAL         D1_TFilt (D1_TFilt_len )
      LOGICAL PUT_STEP ! True when we're going to do a put
                       ! for OASIS4 purposes, false if we use dummy data.
      LOGICAL GET_STEP ! True if we need to get data from OASIS4, false
                       ! if we perform a get but discard results.
      LOGICAL CPL_UPDATE_STEP ! True when we need to make sure
                              ! the D1 prognostic data is updated
                              ! in preparation for dump creation or
                              ! coupling actions.  
      INTEGER OASIS_COUPLE_TS ! OASIS_COUPLE_FREQ (hours) converted to
                              ! a number of timesteps.     
 
      INTEGER stepno
      INTEGER oasis_error  ! Just a flag
      LOGICAL L_OASIS3, L_OASIS4

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
! Error reporting
      INTEGER       ICODE       ! =0 normal exit; >0 error exit
      CHARACTER*256 Cmessage    ! Error message
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='U_MODEL')

! Needed for non-coupled access initialisation
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------


      ICODE=0
      CMESSAGE=''

!L----------------------------------------------------------------------
!L 0. Start Timer call for U_MODEL (NB: not conditional on LTIMER)
!L
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('U_MODEL ',3)
      END IF

! DEPENDS ON: lbc_coup_setup
      Call lbc_coup_setup(LTIMER)

      ICODE=0

!  Routine GETPPX_PART reads those ppxref records which correspond to
!  entries in the stash list into the ppx look-up arrays PPXI, PPXC.
!  It also sets the ppx pointer array PPXPTR. The lengths of PPXI, PPXC
!  have been dynamically allocated to the value of ppxRecs.

! Initialise row number in PPXI, PPXC arrays
      RowNumber = 1

! Initialise lookup and pointer array
      DO I=1,ppxRecs
        DO J=1,PPXREF_CODELEN
          PPXI(I,J)=0
        END DO
        DO J=1,PPXREF_CHARLEN
          PPXC(I,J) = ' '
        END DO
      END DO
      DO I = 1,N_INTERNAL_MODEL
        DO J   = 0,PPXREF_SECTIONS
          DO K = 1,PPXREF_ITEMS
            PPXPTR(I,J,K)=0
          END DO
        END DO
      END DO

! Read in STASHmaster records
      IF (INTERNAL_MODEL_INDEX(A_IM) >  0) THEN
! DEPENDS ON: getppx_part
      CALL GETPPX_PART(NFT,NFTU,'STASHmaster_A',A_IM,RowNumber,         &
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
     &                        ICODE,CMESSAGE)
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
      END IF

!L----------------------------------------------------------------------
!L 1. General initialisation of control and physical data blocks
!L

      ICODE=0

! DEPENDS ON: initial
      CALL INITIAL(                                                     &
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
! ARGSZSP super arrays lengths not dependent on sub-models
     &  SPD1_LEN,SPSTS_LEN,SPBND_LEN,                                   &
! ARGSZSP end
! ARGSZSPA super arrays lengths (atmosphere)
     &  A_SPDUM_LEN,A_SPPTR_LEN,A_SPCON_LEN,A_SPINF_LEN,A_SPANC_LEN,    &
     &  A_SPBND_LEN,A_SPSTS_LEN,                                        &
! ARGSZSPA end
! 5.5 28/02/03 coupling arrays required for if defined RIVERS C.Bunton
! ARGSZSPC super arrays lengths (atmosphere-ocean coupled)
     &  AO_SPCPL_LEN,                                                   &
! ARGSZSPC end
! ARGSP super arrays not dependent on sub-models
     &  SPD1,SPSTS,SPBND,                                               &
! ARGSP end
! ARGSPA super arrays  (atmosphere)
     &  A_SPDUM,A_SPPTR,A_SPCON,A_SPINF,A_SPANC,A_SPBND,A_SPSTS,        &
! ARGSPA end
! 5.5 28/02/03 coupling arrays required for if defined RIVERS C.Bunton
! ARGSPC super array    (atmosphere-ocean coupled)
     &  AO_SPCPL,                                                       &
! ARGSPC end
     &     IAU_lookup,                                                  &
     &     D1_IAU_k4,                                                   &
     &     D1_TFilt,                                                    &
     &     ngrgas,grgas_addr,                                           &
     &     internal_model,submodel,NGROUP,MEANLEV,                      &
     &     isSTASH)  ! FLUME-STASH
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)

! DEPENDS ON: icopenoutput
      CALL ICOPENOUTPUT(runtype_char)

!  Allocate AC scheme arrays using sizes from AC_INIT
      IF (L_AC) THEN
        obs_flag_len = A_MAX_NO_OBS
        ALLOCATE (OBS_FLAG(obs_flag_len))
! 2048 gives enough space in WORK in RDOBS2 when no or very few obs.
        obs_len = A_MAX_OBS_SIZE+2048
        ALLOCATE (OBS(obs_len))
      write(6,*)'U_MODEL - OBS arrays allocated with sizes ',           &
     & A_MAX_NO_OBS,A_MAX_OBS_SIZE+2048
      ELSE
        A_MAX_NO_OBS = 1
        A_MAX_OBS_SIZE = 1
        obs_flag_len = A_MAX_NO_OBS
        ALLOCATE (OBS_FLAG(obs_flag_len))
        obs_len = A_MAX_OBS_SIZE
        ALLOCATE (OBS(obs_len))
      write(6,*)'U_MODEL - OBS arrays allocated with length 1'
      END IF

!L----------------------------------------------------------------------
!L 2. Check for nothing-to-do
!L
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITCHEK',3)
! DEPENDS ON: exitchek
      CALL EXITCHEK( internal_model, LEXITNOW)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITCHEK',4)
      IF (LEXITNOW) GO TO 9999


!L----------------------------------------------------------------------
!L 3. Start group of timesteps

       ! This is the best place to initialise the oasis3 
       ! grid details and transient variables because by this
       ! stage in the code, we should have read the start data and we'll
       ! have the things we need such as land-sea masks available
       IF (L_OASIS .and. l_auscom) THEN

          PUT_STEP=.FALSE.
          GET_STEP=.FALSE.
          CPL_UPDATE_STEP=.FALSE. 

          ! Work out the frequency of coupling timesteps
          OASIS_COUPLE_TS = STEPS_PER_PERIODim(1) *                        &
     &                  (OASIS_COUPLE_FREQ*3600)/SECS_PER_PERIODim(1)

! gol124: auscom coupling
          WRITE(6,*)'U_MODEL: STEPS_PER_PERIODim=',STEPS_PER_PERIODim(1)
          WRITE(6,*)'U_MODEL: OASIS_COUPLE_FREQ=',OASIS_COUPLE_FREQ
          WRITE(6,*)'U_MODEL: SECS_PER_PERIODim=',SECS_PER_PERIODim(1)
          WRITE(6,*)'U_MODEL: OASIS_COUPLE_TS=',OASIS_COUPLE_TS

          ! Save our model start date for comparison with PRISM start 
          ! date. The easiest thing is to copy it to an array which we 
          ! use in a special module as follows:
          ! um_start_date(1:6) = MODEL_BASIS_TIME(1:6)
          ! But for now this is deactivated until we sort out similar
          ! checks in other components. 


          L_OASIS3=.TRUE.
          ! Set timestep for use in OASIS4 calls to be
          ! real version of the standard atmos timestep
          PRISM_Timestep = 1.0*SECS_PER_STEPim(atmos_sm)
          IF (l_auscom) THEN
             WRITE(6,*)'U_MODEL: Calling GRID64_OASIS3'

! DEPENDS ON: OASIS3_GRID64
             CALL OASIS3_GRID64(OASIS_comp_id,MODEL_BASIS_TIME,ICODE)

             WRITE(6,*) 'U_MODEL: PRISM_TIMESTEP=',PRISM_Timestep

             WRITE(6,*) 'U_MODEL: Completed GRID64_OASIS3',ICODE
          else
             WRITE(6,*)'U_MODEL: Skipping GRID64_OASIS3'
          endif


       ELSE
          ! Should come from c_0_dg_c.h 
          access_tfs = 271.35
          ocn_sss=.false.
          ! initilise um variables by "input_atm.nml" file
! DEPENDS ON: OASIS3_UMVARS_INIT
          CALL OASIS3_UMVARS_INIT()


       END IF

       ! Extra initialisation for running coupled model stand-alone
       if (.not. L_AUSCOM ) then
          ocn_sss = .false.
          access_tfs = tfs
       end if

!L----------------------------------------------------------------------
!L 3. Start group of timesteps
!L
   1  CONTINUE
!L----------------------------------------------------------------------
!L 3.1. Start main timestep loop
!L
!L 3.1.1 Increment model time ..
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INCRTIME',3)
! DEPENDS ON: incrtime
      CALL INCRTIME (                                                   &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
     &       internal_model,ICODE,CMESSAGE)
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('INCRTIME',4)

! Keep tabs on PRISM PUT/GET timesteps.
! At the moment we say a put and get timestep are one and the same.
! We're also hard coding to do operations every 48th timestep i.e.
! assuming a half hour atmos TS!. 
      IF (L_OASIS .and. l_auscom) THEN

         ! Is this a genuine exchange timestep.
! gol124: auscom coupling
         PUT_STEP = (MOD(STEPim(atmos_im),OASIS_COUPLE_TS).EQ.0)
         GET_STEP = (MOD(STEPim(atmos_im),OASIS_COUPLE_TS).EQ.1)
         CPL_UPDATE_STEP=(MOD(STEPim(atmos_im),OASIS_COUPLE_TS).EQ.0)
         stepno = STEPim(atmos_im)
         if (PUT_STEP ) write(6,*) "time_step to couple PUT:",stepno
         if (GET_STEP ) write(6,*) "time_step to couple GET:",stepno
         if (CPL_UPDATE_STEP ) write(6,*) "time_step to couple UPDATE:",stepno


         ! Perform coupling exchanges relating to TS-1
         ! avoid unnecessary calls to prism and scatter
         IF (L_OASIS3 .and. get_step) THEN
      
! DEPENDS ON: timer
            IF (LTIMER) CALL TIMER('oasis3get',3)
    
! DEPENDS ON: oasis3_geto2a
            CALL oasis3_geto2a(                                         &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
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


      ! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP, &
      ! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV, &
      ! 1.2: Data variables stored in secondary space.
       P, & 
      ! Pressure on theta levels
       P_THETA_LEVELS, &
      ! Exner pressure on theta levels
       EXNER_THETA_LEVELS, &
      ! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN, &
      ! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF, &
      ! 1.5: Radiation Increments
       SW_INCS, LW_INCS, &
! PAR radiation increment
       DIRPAR, &
      ! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP, &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3, &
!  tropopause-based ozone
       TPPSOZONE, &
      ! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK, &       
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6, &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3, &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD, &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, &
       CO2, CH4_STOCH, O3_STOCH, &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5, &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10, &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15, &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20, &
! CABLE
!       TSOIL_TILE, SMCL_TILE, STHU_TILE, STHF_TILE, SNOW_DEPTH3L,       &
       TSOIL_TILE, SMCL_TILE,  STHF_TILE, SNOW_DEPTH3L,                 &
       SNOW_MASS3L, SNOW_TMP3L, SNOW_RHO3L, SNOW_RHO1L, SNOW_AGE,       & 
       SNOW_FLG3L,                                                      &
! Lestevens Sept 2012 - adding progs for CASACNP
       CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
       NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
       cable_lai,                                                       &

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC, &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC, &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC, &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC, TRACER_LBC, &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND, &
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND, &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND, &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND, &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND, &
       MURK_LBC_TEND, TRACER_LBC_TEND, &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB, &
      ! 2.2: Data variables stored in secondary space.
       PSTAR, &
      ! 2.3: Cloud fields
       CCB, CCT, CCLWP, &
      ! 2.4: Boundary layer fields
       ZH, & 
      ! Standard deviation of turbulent fluctuations of layer 1
       T1_SD, &
      ! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD, &
      ! Number of model levels in the  turbulently mixed layer
       NTML, &
      ! Top level for turb mixing in any decoupled Sc layer
       NTDSC, &
      ! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS, & 
      ! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT, &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN, &
      ! 2.5: Vegetation Ancillary fields
       CANOPY_WATER, SURF_CAP, SURF_RESIST, ROOT_DEPTH, INFILT, &
       VEG_FRAC, LAI, CANHT, Z0, SFA, MDSA, GS, &
      ! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2, &
       OROG_GRAD_X, OROG_GRAD_Y, &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY, &
      ! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS, &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, &
      ! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT, &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND, &
       DUST_MREL1, DUST_MREL2, DUST_MREL3, &
       DUST_MREL4, DUST_MREL5, DUST_MREL6, &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM, &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX, &
      ! Tracer Fluxes - kdcorbin, 05/10
       TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4, &
       TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8, &
       TRACER_FLUX9, TRACER_FLUX10, TRACER_FLUX11, TRACER_FLUX12, &
       TRACER_FLUX13, TRACER_FLUX14, TRACER_FLUX15, TRACER_FLUX16, &
       TRACER_FLUX17, TRACER_FLUX18, TRACER_FLUX19, TRACER_FLUX20, &

      ! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5, &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10, &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15, &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20, &
      !   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX, &
      !   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5, &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9, &
       LAI_PFT, CANHT_PFT, DISTURB_VEG, &
       SOIL_ALB, SOIL_CARB, &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4, &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC, &
       RSP_W_PFT_ACC, RSP_S_ACC, &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4, &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE, &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, &
       DOLR_FIELD, &
       LW_DOWN, SW_TILE_RTS, &
!! MODIFIED BY AT
!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!       TSLAB, TCLIM, HCLIM, CHEAT, OIFLX, UICE, VICE, &
!       SIG11NE, SIG11SE, SIG11SW, SIG11NW, &
!       SIG12NE, SIG12SE, SIG12SW, SIG12NW, &
!       SIG22NE, SIG22SE, SIG22SW, SIG22NW, &
!! END MODIFIED BY AT
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS, &
!   2.15: Fields carried forward from previous version.
!         May not be required
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho, &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP, &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND, &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET, &
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE, &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM, &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT, &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN, &
       C_SOLAR,C_BLUE,C_DOWN,C_LONGWAVE,C_TAUX,C_TAUY,C_WINDMIX, &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_BOTMELTN,C_TOPMELTN,C_LSRAIN, &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_PRESS, C_U10, C_V10, &
! UKCA oxidant fields
        OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA, & 
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR, &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK, &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4, &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC, &
       ARCLDLTA_DL, &
! Convective Cloud Fields
       LCBASE, CCW_RAD, &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM, &

! End arg_atm_fields.h
     &               GET_STEP,cmessage)

! DEPENDS ON: timer
            IF (LTIMER) CALL TIMER('oasis3get',4)


         END IF  ! L_OASIS3=true
      END IF

!L 3.1.2 .. set timestep control switches
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('SETTSCTL',3)
! DEPENDS ON: settsctl
      CALL SETTSCTL (                                                   &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! 6.2 30/06/06 Reduce continuation lines so HadGEM cna compile. R Barnes
!L --------------- STASH arrays             -------------
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
!L
!L 14/6/94  DEF LBOUTA replaced by atmos-V3.4   S.J.Swarbrick
!L 29/7/98  Add 4 new arguments. D. Robinson.
!L 10/11/00 5.2 Add arguments 26 & 27. D.Robinson
!L 18/08/04 6.1 Add arguments 28-33. D.Robinson
!L
!L --------------- Interface arrays out(atmosphere)-------------
A_SPINF(A_IXINF( 1)),A_SPINF(A_IXINF( 2)),A_SPINF(A_IXINF( 3)),A_SPINF         &
(A_IXINF( 4)),A_SPINF(A_IXINF( 5)),A_SPINF(A_IXINF( 6)),A_SPINF(A_IXINF( 7)),  &
A_SPINF(A_IXINF( 8)),A_SPINF(A_IXINF( 9)),A_SPINF(A_IXINF(10)),A_SPINF         &
(A_IXINF(11)),A_SPINF(A_IXINF(12)),A_SPINF(A_IXINF(13)),A_SPINF(A_IXINF(14)),  &
A_SPINF(A_IXINF(15)),A_SPINF(A_IXINF(16)),A_SPINF(A_IXINF(17)),A_SPINF         &
(A_IXINF(18)),A_SPINF(A_IXINF(19)),A_SPINF(A_IXINF(20)),A_SPINF(A_IXINF(21)),  &
A_SPINF(A_IXINF(22)),A_SPINF(A_IXINF(23)),A_SPINF(A_IXINF(24)),A_SPINF         &
(A_IXINF(25)),A_SPINF(A_IXINF(26)),A_SPINF(A_IXINF(27)),A_SPINF(A_IXINF(28)),  &
A_SPINF(A_IXINF(29)),A_SPINF(A_IXINF(30)),A_SPINF(A_IXINF(31)),A_SPINF         &
(A_IXINF(32)),A_SPINF(A_IXINF(33)),A_SPINF(A_IXINF(34)),A_SPINF(A_IXINF(35)),  &
     &         internal_model,.FALSE.,MEANLEV,ICODE,CMESSAGE,&
     &         isSTASH)       ! FLUME-STASH
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('SETTSCTL',4)
!L 3.1.3 If PPfile initialisation time call PP control routine
!L          for instantaneous data (MEANLEV=0)
! FLUME-STASH - Flume_run is set to false in non-flume runs
      IF (LPP .and. .not.Flume_run) THEN
        write(6,*)  "PPCTL_REINIT timestep", STEPim(atmos_im)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('PPCTL_REINIT',3)
! DEPENDS ON: ppctl_reinit
        CALL PPCTL_REINIT(                                              &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
!L
!L 14/6/94  DEF LBOUTA replaced by atmos-V3.4   S.J.Swarbrick
!L 29/7/98  Add 4 new arguments. D. Robinson.
!L 10/11/00 5.2 Add arguments 26 & 27. D.Robinson
!L 18/08/04 6.1 Add arguments 28-33. D.Robinson
!L
!L --------------- Interface arrays out(atmosphere)-------------
A_SPINF(A_IXINF( 1)),A_SPINF(A_IXINF( 2)),A_SPINF(A_IXINF( 3)),A_SPINF         &
(A_IXINF( 4)),A_SPINF(A_IXINF( 5)),A_SPINF(A_IXINF( 6)),A_SPINF(A_IXINF( 7)),  &
A_SPINF(A_IXINF( 8)),A_SPINF(A_IXINF( 9)),A_SPINF(A_IXINF(10)),A_SPINF         &
(A_IXINF(11)),A_SPINF(A_IXINF(12)),A_SPINF(A_IXINF(13)),A_SPINF(A_IXINF(14)),  &
A_SPINF(A_IXINF(15)),A_SPINF(A_IXINF(16)),A_SPINF(A_IXINF(17)),A_SPINF         &
(A_IXINF(18)),A_SPINF(A_IXINF(19)),A_SPINF(A_IXINF(20)),A_SPINF(A_IXINF(21)),  &
A_SPINF(A_IXINF(22)),A_SPINF(A_IXINF(23)),A_SPINF(A_IXINF(24)),A_SPINF         &
(A_IXINF(25)),A_SPINF(A_IXINF(26)),A_SPINF(A_IXINF(27)),A_SPINF(A_IXINF(28)),  &
A_SPINF(A_IXINF(29)),A_SPINF(A_IXINF(30)),A_SPINF(A_IXINF(31)),A_SPINF         &
(A_IXINF(32)),A_SPINF(A_IXINF(33)),A_SPINF(A_IXINF(34)),A_SPINF(A_IXINF(35)),  &
     &         internal_model,PPNAME,ICODE,CMESSAGE)
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('PPCTL_REINIT',4)
      END IF  ! FLUME-STASH

!L       Integrate atmosphere or ocean by 1 timestep
          IF (internal_model == atmos_im) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('ATM_STEP',3)
! Synchronize before the timestep starts
      CALL GC_GSYNC(nproc,info)

! River routing
      IF (L_RIVERS) THEN
! Get 'global' atmos horizontal domain sizes from database
! in DECOMPDB to set dynamic allocation in ATM_STEP for River routing
! on PE 0                                                   .
        G_THETA_FIELD_SIZE=                                             &
     &  decomp_db_glsize(1,fld_type_p,decomp_standard_atmos) *          &
     &  decomp_db_glsize(2,fld_type_p,decomp_standard_atmos)
        G_RIVER_FIELD_SIZE=                                             &
     &  decomp_db_glsize(1,fld_type_r,decomp_standard_atmos) *          &
     &  decomp_db_glsize(2,fld_type_r,decomp_standard_atmos)
      ELSE
        G_THETA_FIELD_SIZE=1
        G_RIVER_FIELD_SIZE=1
      END IF
!

! DEPENDS ON: atm_step
         CALL ATM_STEP (                                                &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! History:
! Version  Date    Comment
!  5.1   13/12/99  New header file : last 3 items of artcona.h
!                  JC thil
!  6.1   09/1//04  Tidy up to avoid too many continuation lines
!                     Anthony A. Dickinson
!L --------------- Derived constants (atmosphere)--------
     & A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
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


      ! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP, &
      ! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV, &
      ! 1.2: Data variables stored in secondary space.
       P, & 
      ! Pressure on theta levels
       P_THETA_LEVELS, &
      ! Exner pressure on theta levels
       EXNER_THETA_LEVELS, &
      ! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN, &
      ! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF, &
      ! 1.5: Radiation Increments
       SW_INCS, LW_INCS, &
! PAR radiation increment
       DIRPAR, &
      ! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP, &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3, &
!  tropopause-based ozone
       TPPSOZONE, &
      ! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK, &       
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6, &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3, &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD, &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, &
       CO2, CH4_STOCH, O3_STOCH, &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5, &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10, &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15, &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20, &
! CABLE
!       TSOIL_TILE, SMCL_TILE, STHU_TILE, STHF_TILE, SNOW_DEPTH3L,       &
       TSOIL_TILE, SMCL_TILE,  STHF_TILE, SNOW_DEPTH3L,                 &
       SNOW_MASS3L, SNOW_TMP3L, SNOW_RHO3L, SNOW_RHO1L, SNOW_AGE,       & 
       SNOW_FLG3L,                                                      &
! Lestevens Sept 2012 - adding progs for CASACNP
       CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
       NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
       cable_lai,                                                       &

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC, &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC, &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC, &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC, TRACER_LBC, &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND, &
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND, &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND, &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND, &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND, &
       MURK_LBC_TEND, TRACER_LBC_TEND, &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB, &
      ! 2.2: Data variables stored in secondary space.
       PSTAR, &
      ! 2.3: Cloud fields
       CCB, CCT, CCLWP, &
      ! 2.4: Boundary layer fields
       ZH, & 
      ! Standard deviation of turbulent fluctuations of layer 1
       T1_SD, &
      ! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD, &
      ! Number of model levels in the  turbulently mixed layer
       NTML, &
      ! Top level for turb mixing in any decoupled Sc layer
       NTDSC, &
      ! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS, & 
      ! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT, &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN, &
      ! 2.5: Vegetation Ancillary fields
       CANOPY_WATER, SURF_CAP, SURF_RESIST, ROOT_DEPTH, INFILT, &
       VEG_FRAC, LAI, CANHT, Z0, SFA, MDSA, GS, &
      ! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2, &
       OROG_GRAD_X, OROG_GRAD_Y, &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY, &
      ! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS, &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, &
      ! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT, &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND, &
       DUST_MREL1, DUST_MREL2, DUST_MREL3, &
       DUST_MREL4, DUST_MREL5, DUST_MREL6, &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM, &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX, &
      ! Tracer Fluxes - kdcorbin, 05/10
       TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4, &
       TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8, &
       TRACER_FLUX9, TRACER_FLUX10, TRACER_FLUX11, TRACER_FLUX12, &
       TRACER_FLUX13, TRACER_FLUX14, TRACER_FLUX15, TRACER_FLUX16, &
       TRACER_FLUX17, TRACER_FLUX18, TRACER_FLUX19, TRACER_FLUX20, &

      ! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5, &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10, &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15, &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20, &
      !   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX, &
      !   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5, &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9, &
       LAI_PFT, CANHT_PFT, DISTURB_VEG, &
       SOIL_ALB, SOIL_CARB, &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4, &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC, &
       RSP_W_PFT_ACC, RSP_S_ACC, &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4, &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE, &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, &
       DOLR_FIELD, &
       LW_DOWN, SW_TILE_RTS, &
!! MODIFIED BY AT
!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!       TSLAB, TCLIM, HCLIM, CHEAT, OIFLX, UICE, VICE, &
!       SIG11NE, SIG11SE, SIG11SW, SIG11NW, &
!       SIG12NE, SIG12SE, SIG12SW, SIG12NW, &
!       SIG22NE, SIG22SE, SIG22SW, SIG22NW, &
!! END MODIFIED BY AT
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS, &
!   2.15: Fields carried forward from previous version.
!         May not be required
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho, &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP, &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND, &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET, &
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE, &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM, &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT, &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN, &
       C_SOLAR,C_BLUE,C_DOWN,C_LONGWAVE,C_TAUX,C_TAUY,C_WINDMIX, &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_BOTMELTN,C_TOPMELTN,C_LSRAIN, &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_PRESS, C_U10, C_V10, &
! UKCA oxidant fields
        OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA, & 
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR, &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK, &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4, &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC, &
       ARCLDLTA_DL, &
! Convective Cloud Fields
       LCBASE, CCW_RAD, &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM, &

! End arg_atm_fields.h
!L --------------- Input boundary arrays ------------------------
     &  SPBND(  IXBND( 1)),                                             &
! 6.2 30/06/06 Reduce continuation lines so HadGEM cna compile. R Barnes
!L --------------- STASH arrays             -------------
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
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
! River routing
! ARTATCPL start
! Description: Super Arrays containing gridline coordinates for
! interpolation and area-averaging between atmosphere and
! river-routing grids (Part of ARTAOCPL.hdk)
! Author: C.Bunton 28.02.03
!
! History
! Version  Date    Comment
!  5.5  28/02/03  Original code. C.Bunton
!L --------------- (Atmosphere-TRIP) coupling arrays  -----------
!L ---------------Lat., Long. values of Atmosphere --------------
     &AO_SPCPL(AO_IXCPL( 5)),AO_SPCPL(AO_IXCPL( 6)),                    &
     &AO_SPCPL(AO_IXCPL( 7)),AO_SPCPL(AO_IXCPL( 8)),                    &
     &AO_SPCPL(AO_IXCPL( 9)),AO_SPCPL(AO_IXCPL( 10)),                   &
! END ARTATCPL
     & G_THETA_FIELD_SIZE,                                              &
     & G_RIVER_FIELD_SIZE,                                              &
     & OBS_FLAG,OBS,obs_flag_len,obs_len,                               &
!
     &      ngrgas,grgas_addr,                                          &
     &                     IAU_lookup,                                  &
                                       ! inout
     &                     D1_IAU_k4,                                   &
                                       ! inout
     &                     D1_TFilt,   &! inout
          !jhan:pass thru end step to CABLE
     &    TARGET_END_STEPim(a_im) )

! DEPENDS ON: timer
      IF (LTIMER) CALL TIMER('ATM_STEP',4)

!L Generate Atmosphere lateral boundary values

      If (linterface) Then

! DEPENDS ON: timer
                            If (LTIMER) Call TIMER ('GEN_INTF',3)

! DEPENDS ON: gen_intf
        Call GEN_INTF (                                                 &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! 6.2 30/06/06 Reduce continuation lines so HadGEM cna compile. R Barnes
!L --------------- STASH arrays             -------------
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add 6 super-array pointers for murk & user ancil. RTHB
!  4.1  04/12/95  Add 2 extra array pointers for STHU and STHF J.Smith
!  4.1  26/04/96  Add 12 extra array pointers for Sulphur Cycle   MJW
!  4.3  26/03/97  And one for HadCM2 sulphate loading pattern.  WJI
!  4.4  05/08/97  And one for convective cloud amount on model levs. JMG
!  4.5  04/03/98   Add 1 pointer for NH3 in S Cycle      M. Woodage
!                  Add 3 pointers for Soot Vars          M. Woodage
!  4.5  08/05/98   Add 16 pointers for new User Anc.     D. Goddard
!  4.5  15/07/98   Add pointer for carbon cycle.         C.D. Jones.
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  4.5  17/08/98   Remove pointer for Soil & Veg fields  D. Robinson
!  5.1  26/04/00   Added pointers for LBC fields   P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2  13/09/00   Removed redundent tracer fields P.Selwood
!  5.3  19/06/01   Add pointers for tropopauase-based ozone Dave Tan
!  5.5  03/02/03   Add pointers for additional microphysics prognostics
!                  and LBCs.                                 R.M.Forbes
!  6.0  30/07/03   Add pointers for cloud fractions       Damian Wilson
!  6.1  07/09/04  Add pointers for STOCHEM fields.  C. Johnson
!  6.1  09/1//04   Tidy up to avoid too many continuation lines
!                        Anthony A. Dickinson
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add pointers for murk aerosol LBCs.    R.M.Forbes
!  6.2  30/06/06   Horrible coding, I know, but absolutely essential to
!                  reduce continuation lines so HadGEM can compile. RTHB
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(94)
      ! tpps_ozone: A_IXPTR(95)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(59)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),                       &
      ! Biomass aerosol
A_SPPTR(A_IXPTR(96)),A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(99)),A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)), &
A_SPPTR(A_IXPTR(102)),A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),&
A_SPPTR(A_IXPTR(105)),A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(108)),A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),&
A_SPPTR(A_IXPTR(111)),A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),&
      ! Cloud fractions
A_SPPTR(A_IXPTR(114)),A_SPPTR(A_IXPTR(115)),A_SPPTR(A_IXPTR(116)),&
A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)),A_SPPTR(A_IXPTR(119)),&
! Add pointer for direct PAR flux
A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)),A_SPPTR(A_IXPTR(122)),&
A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)),A_SPPTR(A_IXPTR(125)),&
      ! Pointers for UKCA oxidant fields (126-129)
A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),A_SPPTR(A_IXPTR(128)),&
A_SPPTR(A_IXPTR(129)),                                            &
      ! Convective Cloud Fields
A_SPPTR(A_IXPTR(130)),A_SPPTR(A_IXPTR(131)),                      &
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),A_SPPTR(A_IXPTR(134)),&
A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),A_SPPTR(A_IXPTR(137)),&
A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),                      &
      ! Pointers for aerosol climatologies (140-160)
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
!Fossil-fuel organic carbon aerosol
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
! CABLE and CASA-CNP + CABLE_LAI 
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),A_SPPTR(A_IXPTR(178)),&
A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),A_SPPTR(A_IXPTR(181)),&
A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),A_SPPTR(A_IXPTR(184)),&
!L
!L 14/6/94  DEF LBOUTA replaced by atmos-V3.4   S.J.Swarbrick
!L 29/7/98  Add 4 new arguments. D. Robinson.
!L 10/11/00 5.2 Add arguments 26 & 27. D.Robinson
!L 18/08/04 6.1 Add arguments 28-33. D.Robinson
!L
!L --------------- Interface arrays out(atmosphere)-------------
A_SPINF(A_IXINF( 1)),A_SPINF(A_IXINF( 2)),A_SPINF(A_IXINF( 3)),A_SPINF         &
(A_IXINF( 4)),A_SPINF(A_IXINF( 5)),A_SPINF(A_IXINF( 6)),A_SPINF(A_IXINF( 7)),  &
A_SPINF(A_IXINF( 8)),A_SPINF(A_IXINF( 9)),A_SPINF(A_IXINF(10)),A_SPINF         &
(A_IXINF(11)),A_SPINF(A_IXINF(12)),A_SPINF(A_IXINF(13)),A_SPINF(A_IXINF(14)),  &
A_SPINF(A_IXINF(15)),A_SPINF(A_IXINF(16)),A_SPINF(A_IXINF(17)),A_SPINF         &
(A_IXINF(18)),A_SPINF(A_IXINF(19)),A_SPINF(A_IXINF(20)),A_SPINF(A_IXINF(21)),  &
A_SPINF(A_IXINF(22)),A_SPINF(A_IXINF(23)),A_SPINF(A_IXINF(24)),A_SPINF         &
(A_IXINF(25)),A_SPINF(A_IXINF(26)),A_SPINF(A_IXINF(27)),A_SPINF(A_IXINF(28)),  &
A_SPINF(A_IXINF(29)),A_SPINF(A_IXINF(30)),A_SPINF(A_IXINF(31)),A_SPINF         &
(A_IXINF(32)),A_SPINF(A_IXINF(33)),A_SPINF(A_IXINF(34)),A_SPINF(A_IXINF(35)),  &
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
     &       submodel, ICode, CMessage)

! DEPENDS ON: timer
                            If (LTIMER) Call TIMER ('GEN_INTF',4)

        If (ICode /= 0) Then
! DEPENDS ON: ereport
          CALL Ereport(RoutineName, ICode ,Cmessage)
        End If

      End If  ! LINTERFACE

      IF (L_ukca) THEN
        WRITE(6,*) ' Before call to UKCA_MAIN1'
! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UKCA_MAIN1',5)
! DEPENDS ON: ukca_main1
        CALL UKCA_MAIN1(                                                &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! 6.2 30/06/06 Reduce continuation lines so HadGEM cna compile. R Barnes
!L --------------- STASH arrays             -------------
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
! History:
! Version  Date    Comment
!  5.1   13/12/99  New header file : last 3 items of artcona.h
!                  JC thil
!  6.1   09/1//04  Tidy up to avoid too many continuation lines
!                     Anthony A. Dickinson
!L --------------- Derived constants (atmosphere)--------
     & A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
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
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add 6 super-array pointers for murk & user ancil. RTHB
!  4.1  04/12/95  Add 2 extra array pointers for STHU and STHF J.Smith
!  4.1  26/04/96  Add 12 extra array pointers for Sulphur Cycle   MJW
!  4.3  26/03/97  And one for HadCM2 sulphate loading pattern.  WJI
!  4.4  05/08/97  And one for convective cloud amount on model levs. JMG
!  4.5  04/03/98   Add 1 pointer for NH3 in S Cycle      M. Woodage
!                  Add 3 pointers for Soot Vars          M. Woodage
!  4.5  08/05/98   Add 16 pointers for new User Anc.     D. Goddard
!  4.5  15/07/98   Add pointer for carbon cycle.         C.D. Jones.
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  4.5  17/08/98   Remove pointer for Soil & Veg fields  D. Robinson
!  5.1  26/04/00   Added pointers for LBC fields   P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2  13/09/00   Removed redundent tracer fields P.Selwood
!  5.3  19/06/01   Add pointers for tropopauase-based ozone Dave Tan
!  5.5  03/02/03   Add pointers for additional microphysics prognostics
!                  and LBCs.                                 R.M.Forbes
!  6.0  30/07/03   Add pointers for cloud fractions       Damian Wilson
!  6.1  07/09/04  Add pointers for STOCHEM fields.  C. Johnson
!  6.1  09/1//04   Tidy up to avoid too many continuation lines
!                        Anthony A. Dickinson
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add pointers for murk aerosol LBCs.    R.M.Forbes
!  6.2  30/06/06   Horrible coding, I know, but absolutely essential to
!                  reduce continuation lines so HadGEM can compile. RTHB
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(94)
      ! tpps_ozone: A_IXPTR(95)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(59)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),                       &
      ! Biomass aerosol
A_SPPTR(A_IXPTR(96)),A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(99)),A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)), &
A_SPPTR(A_IXPTR(102)),A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),&
A_SPPTR(A_IXPTR(105)),A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(108)),A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),&
A_SPPTR(A_IXPTR(111)),A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),&
      ! Cloud fractions
A_SPPTR(A_IXPTR(114)),A_SPPTR(A_IXPTR(115)),A_SPPTR(A_IXPTR(116)),&
A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)),A_SPPTR(A_IXPTR(119)),&
! Add pointer for direct PAR flux
A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)),A_SPPTR(A_IXPTR(122)),&
A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)),A_SPPTR(A_IXPTR(125)),&
      ! Pointers for UKCA oxidant fields (126-129)
A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),A_SPPTR(A_IXPTR(128)),&
A_SPPTR(A_IXPTR(129)),                                            &
      ! Convective Cloud Fields
A_SPPTR(A_IXPTR(130)),A_SPPTR(A_IXPTR(131)),                      &
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),A_SPPTR(A_IXPTR(134)),&
A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),A_SPPTR(A_IXPTR(137)),&
A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),                      &
      ! Pointers for aerosol climatologies (140-160)
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
!Fossil-fuel organic carbon aerosol
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
! CABLE and CASA-CNP + CABLE_LAI 
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),A_SPPTR(A_IXPTR(178)),&
A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),A_SPPTR(A_IXPTR(181)),&
A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),A_SPPTR(A_IXPTR(184)),&
     &   I)                  ! dummy to terminate call
! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UKCA_MAIN1',6)
      END IF

         END IF        ! internal_model = atmos_im

         IF (L_OASIS) THEN

            ! Applicable to oasis3 or OASIS4

            IF (CPL_UPDATE_STEP .and. l_auscom) THEN
!---------------------------------------------------------------------
! Ensure atmos coupling data in prognostic areas of D1 are up to date.
! The logic here is that CPL_UPDATE_STEP will be TRUE on the timestep
! BEFORE coupling is due to take place. 
!
! Newly generated coupling data is intercepted after ATM_STEP and 
! copied to  D1 prognostics.
!
! On the next timestep i.e. a coupling timestep, the prognostic
! contents will be sent to the other components in the coupling process.
!
! If there is no subsequent timestep (i.e. if this is the last model
! timestep then the D1 contents will be written to the dump
! ready for any future restart). There is no "end-of-model"
! coupling exchange. 
!---------------------------------------------------------------------

               If (ltimer) Call timer('updatecpl',3)

! DEPENDS ON: oasis_updatecpl
               CALL OASIS_UPDATECPL(                                           &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
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


      ! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP, &
      ! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV, &
      ! 1.2: Data variables stored in secondary space.
       P, & 
      ! Pressure on theta levels
       P_THETA_LEVELS, &
      ! Exner pressure on theta levels
       EXNER_THETA_LEVELS, &
      ! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN, &
      ! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF, &
      ! 1.5: Radiation Increments
       SW_INCS, LW_INCS, &
! PAR radiation increment
       DIRPAR, &
      ! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP, &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3, &
!  tropopause-based ozone
       TPPSOZONE, &
      ! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK, &       
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6, &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3, &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD, &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, &
       CO2, CH4_STOCH, O3_STOCH, &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5, &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10, &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15, &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20, &
! CABLE
!       TSOIL_TILE, SMCL_TILE, STHU_TILE, STHF_TILE, SNOW_DEPTH3L,       &
       TSOIL_TILE, SMCL_TILE,  STHF_TILE, SNOW_DEPTH3L,                 &
       SNOW_MASS3L, SNOW_TMP3L, SNOW_RHO3L, SNOW_RHO1L, SNOW_AGE,       & 
       SNOW_FLG3L,                                                      &
! Lestevens Sept 2012 - adding progs for CASACNP
       CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
       NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
       cable_lai,                                                       &

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC, &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC, &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC, &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC, TRACER_LBC, &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND, &
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND, &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND, &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND, &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND, &
       MURK_LBC_TEND, TRACER_LBC_TEND, &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB, &
      ! 2.2: Data variables stored in secondary space.
       PSTAR, &
      ! 2.3: Cloud fields
       CCB, CCT, CCLWP, &
      ! 2.4: Boundary layer fields
       ZH, & 
      ! Standard deviation of turbulent fluctuations of layer 1
       T1_SD, &
      ! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD, &
      ! Number of model levels in the  turbulently mixed layer
       NTML, &
      ! Top level for turb mixing in any decoupled Sc layer
       NTDSC, &
      ! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS, & 
      ! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT, &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN, &
      ! 2.5: Vegetation Ancillary fields
       CANOPY_WATER, SURF_CAP, SURF_RESIST, ROOT_DEPTH, INFILT, &
       VEG_FRAC, LAI, CANHT, Z0, SFA, MDSA, GS, &
      ! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2, &
       OROG_GRAD_X, OROG_GRAD_Y, &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY, &
      ! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS, &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, &
      ! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT, &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND, &
       DUST_MREL1, DUST_MREL2, DUST_MREL3, &
       DUST_MREL4, DUST_MREL5, DUST_MREL6, &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM, &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX, &
      ! Tracer Fluxes - kdcorbin, 05/10
       TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4, &
       TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8, &
       TRACER_FLUX9, TRACER_FLUX10, TRACER_FLUX11, TRACER_FLUX12, &
       TRACER_FLUX13, TRACER_FLUX14, TRACER_FLUX15, TRACER_FLUX16, &
       TRACER_FLUX17, TRACER_FLUX18, TRACER_FLUX19, TRACER_FLUX20, &

      ! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5, &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10, &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15, &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20, &
      !   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX, &
      !   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5, &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9, &
       LAI_PFT, CANHT_PFT, DISTURB_VEG, &
       SOIL_ALB, SOIL_CARB, &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4, &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC, &
       RSP_W_PFT_ACC, RSP_S_ACC, &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4, &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE, &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, &
       DOLR_FIELD, &
       LW_DOWN, SW_TILE_RTS, &
!! MODIFIED BY AT
!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!       TSLAB, TCLIM, HCLIM, CHEAT, OIFLX, UICE, VICE, &
!       SIG11NE, SIG11SE, SIG11SW, SIG11NW, &
!       SIG12NE, SIG12SE, SIG12SW, SIG12NW, &
!       SIG22NE, SIG22SE, SIG22SW, SIG22NW, &
!! END MODIFIED BY AT
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS, &
!   2.15: Fields carried forward from previous version.
!         May not be required
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho, &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP, &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND, &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET, &
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE, &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM, &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT, &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN, &
       C_SOLAR,C_BLUE,C_DOWN,C_LONGWAVE,C_TAUX,C_TAUY,C_WINDMIX, &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_BOTMELTN,C_TOPMELTN,C_LSRAIN, &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_PRESS, C_U10, C_V10, &
! UKCA oxidant fields
        OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA, & 
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR, &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK, &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4, &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC, &
       ARCLDLTA_DL, &
! Convective Cloud Fields
       LCBASE, CCW_RAD, &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM, &

! End arg_atm_fields.h
     &               cmessage)

               If (ltimer) Call timer('updatecpl',4)

            END IF

! gol124: auscom coupling
          IF (l_oasis3) THEN
! DEPENDS ON: timer
            IF (LTIMER) CALL TIMER('oasis3put',3)

! DEPENDS ON: oasis3_puta2o
            CALL oasis3_puta2o(                                         &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
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


      ! 1.1: Data variables stored in primary space.
       U, V, W, RHO, THETA, Q, QCL, QCF, QCF2, QRAIN, QGRAUP, &
      ! Exner pressure on rho levels
       EXNER_RHO_LEVELS, U_ADV, V_ADV, W_ADV, &
      ! 1.2: Data variables stored in secondary space.
       P, & 
      ! Pressure on theta levels
       P_THETA_LEVELS, &
      ! Exner pressure on theta levels
       EXNER_THETA_LEVELS, &
      ! 1.3: Cloud Fields
       CCA, CF_AREA, CF_BULK, CF_LIQUID, CF_FROZEN, &
      ! 1.4: Soil Ancillary fields
       DEEP_SOIL_TEMP, SMCL, STHU, STHF, &
      ! 1.5: Radiation Increments
       SW_INCS, LW_INCS, &
! PAR radiation increment
       DIRPAR, &
      ! 1.6: Ozone and cariolle ozone tracers
       O3, OZONE_TRACER,O3_PROD_LOSS,O3_P_L_VMR,O3_VMR,O3_P_L_TEMP, &
       O3_TEMP,O3_P_L_COLO3,O3_COLO3, &
!  tropopause-based ozone
       TPPSOZONE, &
      ! 1.7: Tracer and aerosol fields
       TRACER, TRACER_UKCA, MURK_SOURCE, MURK, &       
       DUST_DIV1, DUST_DIV2, DUST_DIV3, DUST_DIV4, DUST_DIV5, DUST_DIV6, &
       SO2, DMS, SO4_AITKEN, SO4_ACCU, SO4_DISS,  H2O2, NH3, &
       SOOT_NEW, SOOT_AGD, SOOT_CLD, BMASS_NEW, BMASS_AGD, BMASS_CLD, &
       SO2_NATEM, OH, HO2, H2O2_LIMIT, O3_CHEM, &
       CO2, CH4_STOCH, O3_STOCH, &
! 1.8: Multi-level user ancillary fields
       USER_MULT1, USER_MULT2, USER_MULT3, USER_MULT4, USER_MULT5, &
       USER_MULT6, USER_MULT7, USER_MULT8, USER_MULT9, USER_MULT10, &
       USER_MULT11, USER_MULT12, USER_MULT13, USER_MULT14, USER_MULT15, &
       USER_MULT16, USER_MULT17, USER_MULT18, USER_MULT19, USER_MULT20, &
! CABLE
!       TSOIL_TILE, SMCL_TILE, STHU_TILE, STHF_TILE, SNOW_DEPTH3L,       &
       TSOIL_TILE, SMCL_TILE,  STHF_TILE, SNOW_DEPTH3L,                 &
       SNOW_MASS3L, SNOW_TMP3L, SNOW_RHO3L, SNOW_RHO1L, SNOW_AGE,       & 
       SNOW_FLG3L,                                                      &
! Lestevens Sept 2012 - adding progs for CASACNP
       CPOOL_TILE,NPOOL_TILE,PPOOL_TILE,SOIL_ORDER,                     &
       NIDEP,NIFIX,PWEA,PDUST,GLAI,PHENPHASE,                           &
       cable_lai,                                                       &

      ! 1.9: Fields carried forward from previous version.
      ! Lateral Boundary Conditions
       OROG_LBC, U_LBC, V_LBC, W_LBC, RHO_LBC, THETA_LBC, &
       Q_LBC, QCL_LBC, QCF_LBC, QCF2_LBC, QRAIN_LBC, QGRAUP_LBC, &
       CF_BULK_LBC, CF_LIQUID_LBC, CF_FROZEN_LBC, EXNER_LBC, &
       U_ADV_LBC, V_ADV_LBC, W_ADV_LBC, MURK_LBC, TRACER_LBC, &       
      ! Lateral Boundary Condition tendencies
       U_LBC_TEND, V_LBC_TEND, W_LBC_TEND, RHO_LBC_TEND, THETA_LBC_TEND, &
       Q_LBC_TEND, QCL_LBC_TEND, QCF_LBC_TEND, QCF2_LBC_TEND, &
       QRAIN_LBC_TEND, QGRAUP_LBC_TEND, &
       CF_BULK_LBC_TEND, CF_LIQUID_LBC_TEND, CF_FROZEN_LBC_TEND, &
       EXNER_LBC_TEND, U_ADV_LBC_TEND, V_ADV_LBC_TEND, W_ADV_LBC_TEND, &
       MURK_LBC_TEND, TRACER_LBC_TEND, &
      ! 2: Scalar Variables
      ! 2.1: Data variables stored in primary space.
       TSTAR, LAND, TSTAR_ANOM, &
!   2.15: Fields for coastal tiling
       FRAC_LAND, TSTAR_LAND, TSTAR_SEA, TSTAR_SICE, &
! Set pointers for sea-ice and land albedos
       SICE_ALB, LAND_ALB, &
      ! 2.2: Data variables stored in secondary space.
       PSTAR, &
      ! 2.3: Cloud fields
       CCB, CCT, CCLWP, &
      ! 2.4: Boundary layer fields
       ZH, & 
      ! Standard deviation of turbulent fluctuations of layer 1
       T1_SD, &
      ! Standard deviation of turbulent fluctuations of layer 1 humidity
       Q1_SD, &
      ! Number of model levels in the  turbulently mixed layer
       NTML, &
      ! Top level for turb mixing in any decoupled Sc layer
       NTDSC, &
      ! Bottom level for turb mixing in any decoupled Sc layer
       NBDSC, CUMULUS, & 
      ! 2.4: Soil Ancillary fields
       SAT_SOILW_SUCTION, THERM_CAP, THERM_COND, VOL_SMC_CRIT, &
       VOL_SMC_WILT, VOL_SMC_SAT, SAT_SOIL_COND, CLAPP_HORN, &
      ! 2.5: Vegetation Ancillary fields
       CANOPY_WATER, SURF_CAP, SURF_RESIST, ROOT_DEPTH, INFILT, &
       VEG_FRAC, LAI, CANHT, Z0, SFA, MDSA, GS, &
      ! 2.6: Orographic Ancillary fields
       OROGRAPHY, OROG_SD, OROG_SIL, OROG_HO2, &
       OROG_GRAD_X, OROG_GRAD_Y, &
       OROG_GRAD_XX, OROG_GRAD_XY, OROG_GRAD_YY, &
      ! 2.7: Sea/Sea Ice
       U_SEA, V_SEA, U_0_P, V_0_P, ICE_FRACTION, ICE_THICKNESS, &
       TI, ICE_FRACT_CAT, ICE_THICK_CAT, TI_CAT, &
      ! 2.8: Snow
       SNODEP,SNODEP_SEA,SNODEP_SEA_CAT,CATCH_SNOW,SNOW_GRND,SNSOOT, &
! 2.9: aerosol emission fields, including mineral dust parent soil props
       SOIL_CLAY, SOIL_SILT, SOIL_SAND, &
       DUST_MREL1, DUST_MREL2, DUST_MREL3, &
       DUST_MREL4, DUST_MREL5, DUST_MREL6, &
       SO2_EM, DMS_EM, SO2_HILEM, NH3_EM, SOOT_EM, SOOT_HILEM, &
       BMASS_EM, BMASS_HILEM, DMS_CONC, DMS_OFLUX, &
      ! Tracer Fluxes - kdcorbin, 05/10
       TRACER_FLUX1, TRACER_FLUX2, TRACER_FLUX3, TRACER_FLUX4, &
       TRACER_FLUX5, TRACER_FLUX6, TRACER_FLUX7, TRACER_FLUX8, &
       TRACER_FLUX9, TRACER_FLUX10, TRACER_FLUX11, TRACER_FLUX12, &
       TRACER_FLUX13, TRACER_FLUX14, TRACER_FLUX15, TRACER_FLUX16, &
       TRACER_FLUX17, TRACER_FLUX18, TRACER_FLUX19, TRACER_FLUX20, &

      ! 2.10: User ancillary fields
       USER_ANC1, USER_ANC2, USER_ANC3, USER_ANC4, USER_ANC5, &
       USER_ANC6, USER_ANC7, USER_ANC8, USER_ANC9, USER_ANC10, &
       USER_ANC11, USER_ANC12, USER_ANC13, USER_ANC14, USER_ANC15, &
       USER_ANC16, USER_ANC17, USER_ANC18, USER_ANC19, USER_ANC20, &
      !   2.11: Store arrays for energy correction calculation
       NET_FLUX, NET_MFLUX, &
      !   2.12: Tiled Vegetation and Triffid fields
       FRAC_TYP, FRAC_CON1, FRAC_CON2, FRAC_CON3, FRAC_CON4, FRAC_CON5, &
       FRAC_CON6, FRAC_CON7, FRAC_CON8, FRAC_CON9, &
       LAI_PFT, CANHT_PFT, DISTURB_VEG, &
       SOIL_ALB, SOIL_CARB, &
       SOIL_CARB1, SOIL_CARB2, SOIL_CARB3, SOIL_CARB4, &
       NPP_PFT_ACC, G_LF_PFT_ACC, G_PHLF_PFT_ACC, &
       RSP_W_PFT_ACC, RSP_S_ACC, &
       RSP_S_ACC1, RSP_S_ACC2, RSP_S_ACC3, RSP_S_ACC4, &
       CAN_WATER_TILE, CATCH_TILE, INFIL_TILE, RGRAIN_TILE, &
       SNODEP_TILE, TSTAR_TILE, Z0_TILE, &
       DOLR_FIELD, &
       LW_DOWN, SW_TILE_RTS, &
!! MODIFIED BY AT
!! REMOVING SLAB AS PART OF VN7.0
!      !   2.13: Slab Model
!       TSLAB, TCLIM, HCLIM, CHEAT, OIFLX, UICE, VICE, &
!       SIG11NE, SIG11SE, SIG11SW, SIG11NW, &
!       SIG12NE, SIG12SE, SIG12SW, SIG12NW, &
!       SIG22NE, SIG22SE, SIG22SW, SIG22NW, &
!! END MODIFIED BY AT
!   2.14: Carbon cycle fields
       CO2FLUX, CO2_EMITS, &
!   2.15: Fields carried forward from previous version.
!         May not be required
      ! Definition of height(i,j,k) = zsea(k) + C(k)*zorog(i,j)
       zseak_theta, Ck_theta, zseak_rho, Ck_rho, &
!   2.16: Fields for large-scale hydrology scheme.
       TI_MEAN, TI_SIG, FEXP, &
       GAMMA_INT, WATER_TABLE, FSFC_SAT, F_WETLAND, &
       STHZW, A_FSAT, C_FSAT, A_FWET, C_FWET, &
!   2.17: Fields for River routing.
       RIV_SEQUENCE, RIV_DIRECTION, RIV_STORAGE, &
       TOT_SURFROFF, TOT_SUBROFF, RIV_INLANDATM, &
! Fields for grid-to-grid river routing (river routing 2A)
       RIV_IAREA, RIV_SLOPE, RIV_FLOWOBS1, RIV_INEXT, RIV_JNEXT, &
       RIV_LAND, RIV_SUBSTORE, RIV_SURFSTORE, RIV_FLOWIN, RIV_BFLOWIN, &
       C_SOLAR,C_BLUE,C_DOWN,C_LONGWAVE,C_TAUX,C_TAUY,C_WINDMIX, &
       C_SENSIBLE,C_SUBLIM,C_EVAP,C_BOTMELTN,C_TOPMELTN,C_LSRAIN, &
       C_LSSNOW,C_CVRAIN,C_CVSNOW,C_RIVEROUT,C_PRESS, C_U10, C_V10, &
! UKCA oxidant fields
        OH_UKCA, HO2_UKCA, H2O2_UKCA, O3_UKCA, & 
! Aerosol climatologies
       ARCLBIOG_BG, ARCLBIOM_FR, ARCLBIOM_AG, ARCLBIOM_IC, ARCLBLCK_FR, &
       ARCLBLCK_AG, ARCLSSLT_FI, ARCLSSLT_JT, ARCLSULP_AC, ARCLSULP_AK, &
       ARCLSULP_DI, ARCLDUST_B1, ARCLDUST_B2, ARCLDUST_B3, ARCLDUST_B4, &
       ARCLDUST_B5, ARCLDUST_B6, ARCLOCFF_FR, ARCLOCFF_AG, ARCLOCFF_IC, &
       ARCLDLTA_DL, &
! Convective Cloud Fields
       LCBASE, CCW_RAD, &
! Fossil-fuel organic carbon aerosol
       OCFF_NEW, OCFF_AGD, OCFF_CLD, OCFF_EM, OCFF_HILEM, &

! End arg_atm_fields.h
     &              PUT_STEP,cmessage)

! DEPENDS ON: timer
            IF (LTIMER) CALL TIMER('oasis3put',4)

         ! Advance date ready for next timestep (if there is one)
! DEPENDS ON: OASIS3_ADVANCE_DATE64
            CALL OASIS3_ADVANCE_DATE64(oasis_error)
           end if ! l_oasis3
         END IF  ! l_oasis


!L 3.1.4 If dump time, call dump control routine
          IF (LDUMP) THEN
            IF (LTIMER) THEN
! DEPENDS ON: timer
              CALL TIMER('DUMPCTL',5)
! DEPENDS ON: timer
              CALL TIMER('DUMPCTL ',3)
            END IF
! DEPENDS ON: dumpctl
            CALL DUMPCTL (                                              &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! History:
! Version  Date    Comment
!  5.1   13/12/99  New header file : last 3 items of artcona.h
!                  JC thil
!  6.1   09/1//04  Tidy up to avoid too many continuation lines
!                     Anthony A. Dickinson
!L --------------- Derived constants (atmosphere)--------
     & A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add 6 super-array pointers for murk & user ancil. RTHB
!  4.1  04/12/95  Add 2 extra array pointers for STHU and STHF J.Smith
!  4.1  26/04/96  Add 12 extra array pointers for Sulphur Cycle   MJW
!  4.3  26/03/97  And one for HadCM2 sulphate loading pattern.  WJI
!  4.4  05/08/97  And one for convective cloud amount on model levs. JMG
!  4.5  04/03/98   Add 1 pointer for NH3 in S Cycle      M. Woodage
!                  Add 3 pointers for Soot Vars          M. Woodage
!  4.5  08/05/98   Add 16 pointers for new User Anc.     D. Goddard
!  4.5  15/07/98   Add pointer for carbon cycle.         C.D. Jones.
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  4.5  17/08/98   Remove pointer for Soil & Veg fields  D. Robinson
!  5.1  26/04/00   Added pointers for LBC fields   P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2  13/09/00   Removed redundent tracer fields P.Selwood
!  5.3  19/06/01   Add pointers for tropopauase-based ozone Dave Tan
!  5.5  03/02/03   Add pointers for additional microphysics prognostics
!                  and LBCs.                                 R.M.Forbes
!  6.0  30/07/03   Add pointers for cloud fractions       Damian Wilson
!  6.1  07/09/04  Add pointers for STOCHEM fields.  C. Johnson
!  6.1  09/1//04   Tidy up to avoid too many continuation lines
!                        Anthony A. Dickinson
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add pointers for murk aerosol LBCs.    R.M.Forbes
!  6.2  30/06/06   Horrible coding, I know, but absolutely essential to
!                  reduce continuation lines so HadGEM can compile. RTHB
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(94)
      ! tpps_ozone: A_IXPTR(95)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(59)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),                       &
      ! Biomass aerosol
A_SPPTR(A_IXPTR(96)),A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(99)),A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)), &
A_SPPTR(A_IXPTR(102)),A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),&
A_SPPTR(A_IXPTR(105)),A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(108)),A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),&
A_SPPTR(A_IXPTR(111)),A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),&
      ! Cloud fractions
A_SPPTR(A_IXPTR(114)),A_SPPTR(A_IXPTR(115)),A_SPPTR(A_IXPTR(116)),&
A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)),A_SPPTR(A_IXPTR(119)),&
! Add pointer for direct PAR flux
A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)),A_SPPTR(A_IXPTR(122)),&
A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)),A_SPPTR(A_IXPTR(125)),&
      ! Pointers for UKCA oxidant fields (126-129)
A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),A_SPPTR(A_IXPTR(128)),&
A_SPPTR(A_IXPTR(129)),                                            &
      ! Convective Cloud Fields
A_SPPTR(A_IXPTR(130)),A_SPPTR(A_IXPTR(131)),                      &
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),A_SPPTR(A_IXPTR(134)),&
A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),A_SPPTR(A_IXPTR(137)),&
A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),                      &
      ! Pointers for aerosol climatologies (140-160)
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
!Fossil-fuel organic carbon aerosol
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
! CABLE and CASA-CNP + CABLE_LAI 
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),A_SPPTR(A_IXPTR(178)),&
A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),A_SPPTR(A_IXPTR(181)),&
A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),A_SPPTR(A_IXPTR(184)),&
! 6.2 30/06/06 Reduce continuation lines so HadGEM cna compile. R Barnes
!L --------------- STASH arrays             -------------
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
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
     &          submodel,MEANLEV,.false.,'           ',0,               &
     &          ICODE,CMESSAGE)
! DEPENDS ON: ereport
            IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)

            IF (LTIMER) THEN
! DEPENDS ON: timer
              CALL TIMER('DUMPCTL',4)
! DEPENDS ON: timer
              CALL TIMER('DUMPCTL ',6)
            END IF
!L 3.1.4.1 Update interim history file unless means are to follow
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',3)
            IF (.NOT.LMEAN) THEN
              IF (.NOT. (N_SUBMODEL_PARTITION  >  1 .AND. submodel      &
     &          /=  SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION) .AND. &
     &          steps_per_periodim(submodel)  /=                        &
     &          dumpfreqim(submodel) )) THEN
! DEPENDS ON: set_history_values
                CALL SET_HISTORY_VALUES
                CALL GC_SSYNC(nproc,info)
                IF (MYPE  ==  0) THEN
! DEPENDS ON: temphist
                  CALL TEMPHIST(PHIST_UNIT,ICODE,CMESSAGE)
! DEPENDS ON: ereport
            IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
                END IF
              END IF 
            END IF
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',4)
!L 3.1.4.2 If atmosphere timestep recalculate prognostic data and
!L         wrap-around fields using rounded off values
            IF (submodel == atmos_sm) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('RESETATM',3)
! DEPENDS ON: resetatm
       CALL RESETATM(                                                   &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add 6 super-array pointers for murk & user ancil. RTHB
!  4.1  04/12/95  Add 2 extra array pointers for STHU and STHF J.Smith
!  4.1  26/04/96  Add 12 extra array pointers for Sulphur Cycle   MJW
!  4.3  26/03/97  And one for HadCM2 sulphate loading pattern.  WJI
!  4.4  05/08/97  And one for convective cloud amount on model levs. JMG
!  4.5  04/03/98   Add 1 pointer for NH3 in S Cycle      M. Woodage
!                  Add 3 pointers for Soot Vars          M. Woodage
!  4.5  08/05/98   Add 16 pointers for new User Anc.     D. Goddard
!  4.5  15/07/98   Add pointer for carbon cycle.         C.D. Jones.
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  4.5  17/08/98   Remove pointer for Soil & Veg fields  D. Robinson
!  5.1  26/04/00   Added pointers for LBC fields   P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2  13/09/00   Removed redundent tracer fields P.Selwood
!  5.3  19/06/01   Add pointers for tropopauase-based ozone Dave Tan
!  5.5  03/02/03   Add pointers for additional microphysics prognostics
!                  and LBCs.                                 R.M.Forbes
!  6.0  30/07/03   Add pointers for cloud fractions       Damian Wilson
!  6.1  07/09/04  Add pointers for STOCHEM fields.  C. Johnson
!  6.1  09/1//04   Tidy up to avoid too many continuation lines
!                        Anthony A. Dickinson
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add pointers for murk aerosol LBCs.    R.M.Forbes
!  6.2  30/06/06   Horrible coding, I know, but absolutely essential to
!                  reduce continuation lines so HadGEM can compile. RTHB
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(94)
      ! tpps_ozone: A_IXPTR(95)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(59)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),                       &
      ! Biomass aerosol
A_SPPTR(A_IXPTR(96)),A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(99)),A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)), &
A_SPPTR(A_IXPTR(102)),A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),&
A_SPPTR(A_IXPTR(105)),A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(108)),A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),&
A_SPPTR(A_IXPTR(111)),A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),&
      ! Cloud fractions
A_SPPTR(A_IXPTR(114)),A_SPPTR(A_IXPTR(115)),A_SPPTR(A_IXPTR(116)),&
A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)),A_SPPTR(A_IXPTR(119)),&
! Add pointer for direct PAR flux
A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)),A_SPPTR(A_IXPTR(122)),&
A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)),A_SPPTR(A_IXPTR(125)),&
      ! Pointers for UKCA oxidant fields (126-129)
A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),A_SPPTR(A_IXPTR(128)),&
A_SPPTR(A_IXPTR(129)),                                            &
      ! Convective Cloud Fields
A_SPPTR(A_IXPTR(130)),A_SPPTR(A_IXPTR(131)),                      &
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),A_SPPTR(A_IXPTR(134)),&
A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),A_SPPTR(A_IXPTR(137)),&
A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),                      &
      ! Pointers for aerosol climatologies (140-160)
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
!Fossil-fuel organic carbon aerosol
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
! CABLE and CASA-CNP + CABLE_LAI 
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),A_SPPTR(A_IXPTR(178)),&
A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),A_SPPTR(A_IXPTR(181)),&
A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),A_SPPTR(A_IXPTR(184)),&
     &                                ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('RESETATM',4)
              IF (ICODE >  0) GO TO 9999
            END IF
          END IF ! L_DUMP

!L 3.1.5 If printed output time, call print control routine
          IF (LPRINT) THEN

          IF (PrintStatus >= PrStatus_Oper) THEN
         WRITE(6,*) RoutineName,':Warning, Printing of climate global ' &
     &             ,'and zonal diagnostics no longer supported'
          END IF  ! PrintStatus test

          END IF
!L 3.1.6 If interface generation time, generate interface fields
          IF ((internal_model == ocean_im) .or.                         &
     &               (internal_model == wave_im)) THEN

            IF (LINTERFACE) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('GEN_INTF',3)
! DEPENDS ON: gen_intf
              CALL GEN_INTF (                                           &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! 6.2 30/06/06 Reduce continuation lines so HadGEM cna compile. R Barnes
!L --------------- STASH arrays             -------------
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add 6 super-array pointers for murk & user ancil. RTHB
!  4.1  04/12/95  Add 2 extra array pointers for STHU and STHF J.Smith
!  4.1  26/04/96  Add 12 extra array pointers for Sulphur Cycle   MJW
!  4.3  26/03/97  And one for HadCM2 sulphate loading pattern.  WJI
!  4.4  05/08/97  And one for convective cloud amount on model levs. JMG
!  4.5  04/03/98   Add 1 pointer for NH3 in S Cycle      M. Woodage
!                  Add 3 pointers for Soot Vars          M. Woodage
!  4.5  08/05/98   Add 16 pointers for new User Anc.     D. Goddard
!  4.5  15/07/98   Add pointer for carbon cycle.         C.D. Jones.
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  4.5  17/08/98   Remove pointer for Soil & Veg fields  D. Robinson
!  5.1  26/04/00   Added pointers for LBC fields   P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2  13/09/00   Removed redundent tracer fields P.Selwood
!  5.3  19/06/01   Add pointers for tropopauase-based ozone Dave Tan
!  5.5  03/02/03   Add pointers for additional microphysics prognostics
!                  and LBCs.                                 R.M.Forbes
!  6.0  30/07/03   Add pointers for cloud fractions       Damian Wilson
!  6.1  07/09/04  Add pointers for STOCHEM fields.  C. Johnson
!  6.1  09/1//04   Tidy up to avoid too many continuation lines
!                        Anthony A. Dickinson
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add pointers for murk aerosol LBCs.    R.M.Forbes
!  6.2  30/06/06   Horrible coding, I know, but absolutely essential to
!                  reduce continuation lines so HadGEM can compile. RTHB
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(94)
      ! tpps_ozone: A_IXPTR(95)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(59)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),                       &
      ! Biomass aerosol
A_SPPTR(A_IXPTR(96)),A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(99)),A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)), &
A_SPPTR(A_IXPTR(102)),A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),&
A_SPPTR(A_IXPTR(105)),A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(108)),A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),&
A_SPPTR(A_IXPTR(111)),A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),&
      ! Cloud fractions
A_SPPTR(A_IXPTR(114)),A_SPPTR(A_IXPTR(115)),A_SPPTR(A_IXPTR(116)),&
A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)),A_SPPTR(A_IXPTR(119)),&
! Add pointer for direct PAR flux
A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)),A_SPPTR(A_IXPTR(122)),&
A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)),A_SPPTR(A_IXPTR(125)),&
      ! Pointers for UKCA oxidant fields (126-129)
A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),A_SPPTR(A_IXPTR(128)),&
A_SPPTR(A_IXPTR(129)),                                            &
      ! Convective Cloud Fields
A_SPPTR(A_IXPTR(130)),A_SPPTR(A_IXPTR(131)),                      &
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),A_SPPTR(A_IXPTR(134)),&
A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),A_SPPTR(A_IXPTR(137)),&
A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),                      &
      ! Pointers for aerosol climatologies (140-160)
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
!Fossil-fuel organic carbon aerosol
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
! CABLE and CASA-CNP + CABLE_LAI 
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),A_SPPTR(A_IXPTR(178)),&
A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),A_SPPTR(A_IXPTR(181)),&
A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),A_SPPTR(A_IXPTR(184)),&
!L
!L 14/6/94  DEF LBOUTA replaced by atmos-V3.4   S.J.Swarbrick
!L 29/7/98  Add 4 new arguments. D. Robinson.
!L 10/11/00 5.2 Add arguments 26 & 27. D.Robinson
!L 18/08/04 6.1 Add arguments 28-33. D.Robinson
!L
!L --------------- Interface arrays out(atmosphere)-------------
A_SPINF(A_IXINF( 1)),A_SPINF(A_IXINF( 2)),A_SPINF(A_IXINF( 3)),A_SPINF         &
(A_IXINF( 4)),A_SPINF(A_IXINF( 5)),A_SPINF(A_IXINF( 6)),A_SPINF(A_IXINF( 7)),  &
A_SPINF(A_IXINF( 8)),A_SPINF(A_IXINF( 9)),A_SPINF(A_IXINF(10)),A_SPINF         &
(A_IXINF(11)),A_SPINF(A_IXINF(12)),A_SPINF(A_IXINF(13)),A_SPINF(A_IXINF(14)),  &
A_SPINF(A_IXINF(15)),A_SPINF(A_IXINF(16)),A_SPINF(A_IXINF(17)),A_SPINF         &
(A_IXINF(18)),A_SPINF(A_IXINF(19)),A_SPINF(A_IXINF(20)),A_SPINF(A_IXINF(21)),  &
A_SPINF(A_IXINF(22)),A_SPINF(A_IXINF(23)),A_SPINF(A_IXINF(24)),A_SPINF         &
(A_IXINF(25)),A_SPINF(A_IXINF(26)),A_SPINF(A_IXINF(27)),A_SPINF(A_IXINF(28)),  &
A_SPINF(A_IXINF(29)),A_SPINF(A_IXINF(30)),A_SPINF(A_IXINF(31)),A_SPINF         &
(A_IXINF(32)),A_SPINF(A_IXINF(33)),A_SPINF(A_IXINF(34)),A_SPINF(A_IXINF(35)),  &
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
     &              submodel,ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('GEN_INTF',4)
! DEPENDS ON: ereport
            IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
            END IF

          END IF
!L 3.1.6.1 Release job to process output created so far, if selected
          IF (LJOBRELEASE) THEN
! DEPENDS ON: flush_all_pp
            CALL FLUSH_ALL_PP()
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('JOBCTL  ',3)
! DEPENDS ON: jobctl
            CALL JOBCTL(internal_model,ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('JOBCTL  ',4)
! DEPENDS ON: ereport
          IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
          END IF
!L 3.1.7 If partial sum/mean creation time, call means control routine
!L       (calls mean PPfield and diagnostic print routines internally)
          IF (LMEAN) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('MEANCTL ',3)
! DEPENDS ON: meanctl
            CALL MEANCTL (                                              &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add 6 super-array pointers for murk & user ancil. RTHB
!  4.1  04/12/95  Add 2 extra array pointers for STHU and STHF J.Smith
!  4.1  26/04/96  Add 12 extra array pointers for Sulphur Cycle   MJW
!  4.3  26/03/97  And one for HadCM2 sulphate loading pattern.  WJI
!  4.4  05/08/97  And one for convective cloud amount on model levs. JMG
!  4.5  04/03/98   Add 1 pointer for NH3 in S Cycle      M. Woodage
!                  Add 3 pointers for Soot Vars          M. Woodage
!  4.5  08/05/98   Add 16 pointers for new User Anc.     D. Goddard
!  4.5  15/07/98   Add pointer for carbon cycle.         C.D. Jones.
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  4.5  17/08/98   Remove pointer for Soil & Veg fields  D. Robinson
!  5.1  26/04/00   Added pointers for LBC fields   P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2  13/09/00   Removed redundent tracer fields P.Selwood
!  5.3  19/06/01   Add pointers for tropopauase-based ozone Dave Tan
!  5.5  03/02/03   Add pointers for additional microphysics prognostics
!                  and LBCs.                                 R.M.Forbes
!  6.0  30/07/03   Add pointers for cloud fractions       Damian Wilson
!  6.1  07/09/04  Add pointers for STOCHEM fields.  C. Johnson
!  6.1  09/1//04   Tidy up to avoid too many continuation lines
!                        Anthony A. Dickinson
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add pointers for murk aerosol LBCs.    R.M.Forbes
!  6.2  30/06/06   Horrible coding, I know, but absolutely essential to
!                  reduce continuation lines so HadGEM can compile. RTHB
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(94)
      ! tpps_ozone: A_IXPTR(95)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(59)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),                       &
      ! Biomass aerosol
A_SPPTR(A_IXPTR(96)),A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(99)),A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)), &
A_SPPTR(A_IXPTR(102)),A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),&
A_SPPTR(A_IXPTR(105)),A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(108)),A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),&
A_SPPTR(A_IXPTR(111)),A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),&
      ! Cloud fractions
A_SPPTR(A_IXPTR(114)),A_SPPTR(A_IXPTR(115)),A_SPPTR(A_IXPTR(116)),&
A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)),A_SPPTR(A_IXPTR(119)),&
! Add pointer for direct PAR flux
A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)),A_SPPTR(A_IXPTR(122)),&
A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)),A_SPPTR(A_IXPTR(125)),&
      ! Pointers for UKCA oxidant fields (126-129)
A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),A_SPPTR(A_IXPTR(128)),&
A_SPPTR(A_IXPTR(129)),                                            &
      ! Convective Cloud Fields
A_SPPTR(A_IXPTR(130)),A_SPPTR(A_IXPTR(131)),                      &
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),A_SPPTR(A_IXPTR(134)),&
A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),A_SPPTR(A_IXPTR(137)),&
A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),                      &
      ! Pointers for aerosol climatologies (140-160)
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
!Fossil-fuel organic carbon aerosol
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
! CABLE and CASA-CNP + CABLE_LAI 
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),A_SPPTR(A_IXPTR(178)),&
A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),A_SPPTR(A_IXPTR(181)),&
A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),A_SPPTR(A_IXPTR(184)),&
! 6.2 30/06/06 Reduce continuation lines so HadGEM cna compile. R Barnes
!L --------------- STASH arrays             -------------
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
! History:
! Version  Date    Comment
!  5.1   13/12/99  New header file : last 3 items of artcona.h
!                  JC thil
!  6.1   09/1//04  Tidy up to avoid too many continuation lines
!                     Anthony A. Dickinson
!L --------------- Derived constants (atmosphere)--------
     & A_SPCON(A_IXCON(1)),A_SPCON(A_IXCON(2)),A_SPCON(A_IXCON(3)),  &
!L
!L 14/6/94  DEF LBOUTA replaced by atmos-V3.4   S.J.Swarbrick
!L 29/7/98  Add 4 new arguments. D. Robinson.
!L 10/11/00 5.2 Add arguments 26 & 27. D.Robinson
!L 18/08/04 6.1 Add arguments 28-33. D.Robinson
!L
!L --------------- Interface arrays out(atmosphere)-------------
A_SPINF(A_IXINF( 1)),A_SPINF(A_IXINF( 2)),A_SPINF(A_IXINF( 3)),A_SPINF         &
(A_IXINF( 4)),A_SPINF(A_IXINF( 5)),A_SPINF(A_IXINF( 6)),A_SPINF(A_IXINF( 7)),  &
A_SPINF(A_IXINF( 8)),A_SPINF(A_IXINF( 9)),A_SPINF(A_IXINF(10)),A_SPINF         &
(A_IXINF(11)),A_SPINF(A_IXINF(12)),A_SPINF(A_IXINF(13)),A_SPINF(A_IXINF(14)),  &
A_SPINF(A_IXINF(15)),A_SPINF(A_IXINF(16)),A_SPINF(A_IXINF(17)),A_SPINF         &
(A_IXINF(18)),A_SPINF(A_IXINF(19)),A_SPINF(A_IXINF(20)),A_SPINF(A_IXINF(21)),  &
A_SPINF(A_IXINF(22)),A_SPINF(A_IXINF(23)),A_SPINF(A_IXINF(24)),A_SPINF         &
(A_IXINF(25)),A_SPINF(A_IXINF(26)),A_SPINF(A_IXINF(27)),A_SPINF(A_IXINF(28)),  &
A_SPINF(A_IXINF(29)),A_SPINF(A_IXINF(30)),A_SPINF(A_IXINF(31)),A_SPINF         &
(A_IXINF(32)),A_SPINF(A_IXINF(33)),A_SPINF(A_IXINF(34)),A_SPINF(A_IXINF(35)),  &
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
     &                  submodel,MEANLEV,ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('MEANCTL ',4)
            IF (ICODE >  0) THEN
! DEPENDS ON: del_hist
              CALL DEL_HIST(PHIST_UNIT)
      WRITE(6,*)'U_MODEL: interim history file deleted due to failu     &
     &re writing partial sum files'
! DEPENDS ON: ereport
              CALL Ereport(RoutineName,ICODE,Cmessage)
            END IF
!L 3.1.7.1 On successful completion, update interim history file
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',3)

            IF (.NOT. (N_SUBMODEL_PARTITION  >  1 .AND. submodel        &
     &          /=  SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION) .AND. &
     &          steps_per_periodim(submodel)  /=                        &
     &          dumpfreqim(submodel) )) THEN

! DEPENDS ON: set_history_values
              CALL SET_HISTORY_VALUES
              CALL GC_SSYNC(nproc,info)
              IF (MYPE  ==  0) THEN
! DEPENDS ON: temphist
                CALL TEMPHIST(PHIST_UNIT,ICODE,CMESSAGE)
              END IF
! DEPENDS ON: ereport
            IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
            END IF
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',4)
          END IF
!L 3.1.8 Update temporary history file if at a 'safe' restart point
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',3)
          IF (LHISTORY) THEN
!          In coupled model do not update history file until both
!          submodels have reached the safe restart point

            IF (.NOT. (N_SUBMODEL_PARTITION  >  1 .AND. submodel        &
     &          /=  SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION) .AND. &
     &          steps_per_periodim(submodel)  /=                        &
     &          dumpfreqim(submodel) )) THEN
! DEPENDS ON: set_history_values
              CALL SET_HISTORY_VALUES
            CALL GC_SSYNC(nproc,info)
            IF (MYPE  ==  0) THEN
! DEPENDS ON: temphist
            CALL TEMPHIST(THIST_UNIT,ICODE,CMESSAGE)
            END IF
          END IF
        END IF
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('TEMPHIST',4)
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
!L 3.1.9 If exit check time, check for immediate exit
          IF (LEXIT) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITCHEK',3)
! DEPENDS ON: exitchek
            CALL EXITCHEK(internal_model, LEXITNOW)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITCHEK',4)
            IF (LEXITNOW) THEN
              IF (.NOT.LDUMP) THEN

                IF (.NOT. (N_SUBMODEL_PARTITION  >  1 .AND. submodel    &
     &           /=  SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION).AND. &
     &          steps_per_periodim(submodel)  /=                        &
     &          dumpfreqim(submodel) )) THEN

! DEPENDS ON: set_history_values
                CALL SET_HISTORY_VALUES
                CALL GC_SSYNC(nproc,info)
                IF (MYPE  ==  0) THEN
! DEPENDS ON: temphist
                CALL TEMPHIST(PHIST_UNIT,ICODE,CMESSAGE)
                END IF
!               Exit model so no need to set ppflush
         END IF
              END IF
              GO TO 9999
      END IF
          END IF
!L 3.1.10 Update ancillary fields if necessary
          IF (LANCILLARY) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('UP_ANCIL',3)
! DEPENDS ON: up_ancil
         CALL UP_ANCIL (                                                &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add 6 super-array pointers for murk & user ancil. RTHB
!  4.1  04/12/95  Add 2 extra array pointers for STHU and STHF J.Smith
!  4.1  26/04/96  Add 12 extra array pointers for Sulphur Cycle   MJW
!  4.3  26/03/97  And one for HadCM2 sulphate loading pattern.  WJI
!  4.4  05/08/97  And one for convective cloud amount on model levs. JMG
!  4.5  04/03/98   Add 1 pointer for NH3 in S Cycle      M. Woodage
!                  Add 3 pointers for Soot Vars          M. Woodage
!  4.5  08/05/98   Add 16 pointers for new User Anc.     D. Goddard
!  4.5  15/07/98   Add pointer for carbon cycle.         C.D. Jones.
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  4.5  17/08/98   Remove pointer for Soil & Veg fields  D. Robinson
!  5.1  26/04/00   Added pointers for LBC fields   P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2  13/09/00   Removed redundent tracer fields P.Selwood
!  5.3  19/06/01   Add pointers for tropopauase-based ozone Dave Tan
!  5.5  03/02/03   Add pointers for additional microphysics prognostics
!                  and LBCs.                                 R.M.Forbes
!  6.0  30/07/03   Add pointers for cloud fractions       Damian Wilson
!  6.1  07/09/04  Add pointers for STOCHEM fields.  C. Johnson
!  6.1  09/1//04   Tidy up to avoid too many continuation lines
!                        Anthony A. Dickinson
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add pointers for murk aerosol LBCs.    R.M.Forbes
!  6.2  30/06/06   Horrible coding, I know, but absolutely essential to
!                  reduce continuation lines so HadGEM can compile. RTHB
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(94)
      ! tpps_ozone: A_IXPTR(95)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(59)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),                       &
      ! Biomass aerosol
A_SPPTR(A_IXPTR(96)),A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(99)),A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)), &
A_SPPTR(A_IXPTR(102)),A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),&
A_SPPTR(A_IXPTR(105)),A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(108)),A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),&
A_SPPTR(A_IXPTR(111)),A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),&
      ! Cloud fractions
A_SPPTR(A_IXPTR(114)),A_SPPTR(A_IXPTR(115)),A_SPPTR(A_IXPTR(116)),&
A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)),A_SPPTR(A_IXPTR(119)),&
! Add pointer for direct PAR flux
A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)),A_SPPTR(A_IXPTR(122)),&
A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)),A_SPPTR(A_IXPTR(125)),&
      ! Pointers for UKCA oxidant fields (126-129)
A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),A_SPPTR(A_IXPTR(128)),&
A_SPPTR(A_IXPTR(129)),                                            &
      ! Convective Cloud Fields
A_SPPTR(A_IXPTR(130)),A_SPPTR(A_IXPTR(131)),                      &
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),A_SPPTR(A_IXPTR(134)),&
A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),A_SPPTR(A_IXPTR(137)),&
A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),                      &
      ! Pointers for aerosol climatologies (140-160)
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
!Fossil-fuel organic carbon aerosol
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
! CABLE and CASA-CNP + CABLE_LAI 
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),A_SPPTR(A_IXPTR(178)),&
A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),A_SPPTR(A_IXPTR(181)),&
A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),A_SPPTR(A_IXPTR(184)),&
!L
!L --------------- Ancillary file arrays(atmosphere)-------------
     &A_SPANC(A_IXANC( 1)),A_SPANC(A_IXANC( 2)),A_SPANC(A_IXANC( 3)),   &
     &A_SPANC(A_IXANC( 4)),                                             &
     &                   submodel,                                      &
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
     &                   ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('UP_ANCIL',4)
! DEPENDS ON: ereport
      IF (ICODE  /=  0) CALL Ereport(RoutineName,ICODE,Cmessage)
          END IF
!L 3.1.11 Update boundary fields if necessary
          IF (LBOUNDARY) THEN
             ! DEPENDS ON: lbc_coup_update
             call lbc_coup_update( &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add 6 super-array pointers for murk & user ancil. RTHB
!  4.1  04/12/95  Add 2 extra array pointers for STHU and STHF J.Smith
!  4.1  26/04/96  Add 12 extra array pointers for Sulphur Cycle   MJW
!  4.3  26/03/97  And one for HadCM2 sulphate loading pattern.  WJI
!  4.4  05/08/97  And one for convective cloud amount on model levs. JMG
!  4.5  04/03/98   Add 1 pointer for NH3 in S Cycle      M. Woodage
!                  Add 3 pointers for Soot Vars          M. Woodage
!  4.5  08/05/98   Add 16 pointers for new User Anc.     D. Goddard
!  4.5  15/07/98   Add pointer for carbon cycle.         C.D. Jones.
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  4.5  17/08/98   Remove pointer for Soil & Veg fields  D. Robinson
!  5.1  26/04/00   Added pointers for LBC fields   P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2  13/09/00   Removed redundent tracer fields P.Selwood
!  5.3  19/06/01   Add pointers for tropopauase-based ozone Dave Tan
!  5.5  03/02/03   Add pointers for additional microphysics prognostics
!                  and LBCs.                                 R.M.Forbes
!  6.0  30/07/03   Add pointers for cloud fractions       Damian Wilson
!  6.1  07/09/04  Add pointers for STOCHEM fields.  C. Johnson
!  6.1  09/1//04   Tidy up to avoid too many continuation lines
!                        Anthony A. Dickinson
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add pointers for murk aerosol LBCs.    R.M.Forbes
!  6.2  30/06/06   Horrible coding, I know, but absolutely essential to
!                  reduce continuation lines so HadGEM can compile. RTHB
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(94)
      ! tpps_ozone: A_IXPTR(95)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(59)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),                       &
      ! Biomass aerosol
A_SPPTR(A_IXPTR(96)),A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(99)),A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)), &
A_SPPTR(A_IXPTR(102)),A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),&
A_SPPTR(A_IXPTR(105)),A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(108)),A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),&
A_SPPTR(A_IXPTR(111)),A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),&
      ! Cloud fractions
A_SPPTR(A_IXPTR(114)),A_SPPTR(A_IXPTR(115)),A_SPPTR(A_IXPTR(116)),&
A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)),A_SPPTR(A_IXPTR(119)),&
! Add pointer for direct PAR flux
A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)),A_SPPTR(A_IXPTR(122)),&
A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)),A_SPPTR(A_IXPTR(125)),&
      ! Pointers for UKCA oxidant fields (126-129)
A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),A_SPPTR(A_IXPTR(128)),&
A_SPPTR(A_IXPTR(129)),                                            &
      ! Convective Cloud Fields
A_SPPTR(A_IXPTR(130)),A_SPPTR(A_IXPTR(131)),                      &
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),A_SPPTR(A_IXPTR(134)),&
A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),A_SPPTR(A_IXPTR(137)),&
A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),                      &
      ! Pointers for aerosol climatologies (140-160)
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
!Fossil-fuel organic carbon aerosol
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
! CABLE and CASA-CNP + CABLE_LAI 
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),A_SPPTR(A_IXPTR(178)),&
A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),A_SPPTR(A_IXPTR(181)),&
A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),A_SPPTR(A_IXPTR(184)),&
!L --------------- Input boundary arrays ------------------------
     &  SPBND(  IXBND( 1)),                                             &
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
                  submodel,ICODE,CMESSAGE)
          END IF
!L
!L      End main timestep loop
!L----------------------------------------------------------------------


      GO TO 1
!
 9999 CONTINUE


!L----------------------------------------------------------------------
!L 4. Exit processing: Output error messages and perform tidy-up
!L
      IF (L_OASIS .and. l_auscom) THEN
! DEPENDS ON: oasis_tidy
         CALL oasis_tidy(                                        &
! History:
! Version  Date    Comment
!  4.2   11/10/96  Enable atmos-ocean coupling for mpp.
!                  (2): Swap D1 memory. Add copies of D1 for atmos and
!                  ocean. R.Rawlins
!L --------------- D1: main data array      -------------
!L ------ (with extra copy for logical values)-----------
     &  SPD1(IXD1( 1)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), SPD1(IXD1( 2)), &
! 6.2 30/06/06 Reduce continuation lines so HadGEM cna compile. R Barnes
!L --------------- STASH arrays             -------------
     &SPSTS(IXSTS(1)),SPSTS(IXSTS(2)),SPSTS(IXSTS(3)),SPSTS(IXSTS(4)),  &
     &SPSTS(IXSTS(5)),SPSTS(IXSTS(6)),SPSTS(IXSTS(7)),SPSTS(IXSTS(8)),SP&
     &STS(IXSTS(9)),SPSTS(IXSTS(10)),SPSTS(IXSTS(11)),SPSTS(IXSTS(12)), &
!L --------------- Dump headers (atmosphere)-------------
     &A_SPDUM(A_IXDUM( 1)),A_SPDUM(A_IXDUM( 2)),A_SPDUM(A_IXDUM( 3)),   &
     &A_SPDUM(A_IXDUM( 4)),A_SPDUM(A_IXDUM( 5)),A_SPDUM(A_IXDUM( 6)),   &
     &A_SPDUM(A_IXDUM( 7)),A_SPDUM(A_IXDUM( 8)),A_SPDUM(A_IXDUM( 9)),   &
     &A_SPDUM(A_IXDUM(10)),A_SPDUM(A_IXDUM(11)),A_SPDUM(A_IXDUM(12)),   &
     &A_SPDUM(A_IXDUM(13)),A_SPDUM(A_IXDUM(14)),a_ixsts, a_spsts,       &
! History:
! Version  Date     Comment
! -------  ----     -------
!  3.4  05/10/94  Add 6 super-array pointers for murk & user ancil. RTHB
!  4.1  04/12/95  Add 2 extra array pointers for STHU and STHF J.Smith
!  4.1  26/04/96  Add 12 extra array pointers for Sulphur Cycle   MJW
!  4.3  26/03/97  And one for HadCM2 sulphate loading pattern.  WJI
!  4.4  05/08/97  And one for convective cloud amount on model levs. JMG
!  4.5  04/03/98   Add 1 pointer for NH3 in S Cycle      M. Woodage
!                  Add 3 pointers for Soot Vars          M. Woodage
!  4.5  08/05/98   Add 16 pointers for new User Anc.     D. Goddard
!  4.5  15/07/98   Add pointer for carbon cycle.         C.D. Jones.
!  5.0    20/05/99 Extensive revision of pointers for new dynamics data
!                  variables. D.M. Goddard
!  4.5  17/08/98   Remove pointer for Soil & Veg fields  D. Robinson
!  5.1  26/04/00   Added pointers for LBC fields   P.Burton
!  5.2  25/09/00  Clear out RHcrit legacy code.    A.C.Bushell.
!
!  5.2  13/09/00   Removed redundent tracer fields P.Selwood
!  5.3  19/06/01   Add pointers for tropopauase-based ozone Dave Tan
!  5.5  03/02/03   Add pointers for additional microphysics prognostics
!                  and LBCs.                                 R.M.Forbes
!  6.0  30/07/03   Add pointers for cloud fractions       Damian Wilson
!  6.1  07/09/04  Add pointers for STOCHEM fields.  C. Johnson
!  6.1  09/1//04   Tidy up to avoid too many continuation lines
!                        Anthony A. Dickinson
!  6.2  01/03/06   Add pointer for direct component of total downward
!                  surface PAR flux.  M.G. Sanderson.
!  6.2  14/11/05   Add JTR_UKCA pointers for UKCA tracers. R Barnes
!  6.2  01/10/05   Add pointers for murk aerosol LBCs.    R.M.Forbes
!  6.2  30/06/06   Horrible coding, I know, but absolutely essential to
!                  reduce continuation lines so HadGEM can compile. RTHB
! --------------- Pointers in D1 & dep consts(atmosphere)-------
! Data variables stored in primary space
      ! Cloud Fields: A_IXPTR(16)- A_IXPTR(20)
      ! & Soil Fields:  A_IXPTR(21)- A_IXPTR(24)
      ! Radiation Increments and Ozone: A_IXPTR(25 - 27)
      ! Tracers and Aerosols: A_IXPTR(28 - 47)
A_SPPTR(A_IXPTR(1)),A_SPPTR(A_IXPTR(2)),A_SPPTR(A_IXPTR(3)),A_SPPTR(A_IXPTR(4))&
,A_SPPTR(A_IXPTR(5)),A_SPPTR(A_IXPTR(6)),A_SPPTR(A_IXPTR(7)),A_SPPTR           &
(A_IXPTR(8)),A_SPPTR(A_IXPTR(9)),A_SPPTR(A_IXPTR(10)),A_SPPTR(A_IXPTR(11)),    &
A_SPPTR(A_IXPTR(12)),A_SPPTR(A_IXPTR(13)),A_SPPTR(A_IXPTR(14)),A_SPPTR         &
(A_IXPTR(15)),A_SPPTR(A_IXPTR(16)),A_SPPTR(A_IXPTR(17)),A_SPPTR(A_IXPTR(18)),  &
A_SPPTR(A_IXPTR(19)),A_SPPTR(A_IXPTR(20)),A_SPPTR(A_IXPTR(21)),A_SPPTR         &
(A_IXPTR(22)),A_SPPTR(A_IXPTR(23)),A_SPPTR(A_IXPTR(24)),A_SPPTR(A_IXPTR(25)),  &
A_SPPTR(A_IXPTR(26)),A_SPPTR(A_IXPTR(27)),A_SPPTR(A_IXPTR(28)),A_SPPTR         &
(A_IXPTR(29)),A_SPPTR(A_IXPTR(30)),A_SPPTR(A_IXPTR(31)),A_SPPTR(A_IXPTR(32)),  &
A_SPPTR(A_IXPTR(33)),A_SPPTR(A_IXPTR(34)),A_SPPTR(A_IXPTR(35)),A_SPPTR         &
(A_IXPTR(36)),A_SPPTR(A_IXPTR(37)),A_SPPTR(A_IXPTR(38)),A_SPPTR(A_IXPTR(39)),  &
A_SPPTR(A_IXPTR(40)),A_SPPTR(A_IXPTR(41)),A_SPPTR(A_IXPTR(42)),A_SPPTR         &
(A_IXPTR(43)),A_SPPTR(A_IXPTR(44)),A_SPPTR(A_IXPTR(45)),A_SPPTR(A_IXPTR(46)),  &
A_SPPTR(A_IXPTR(47)),                                                          &
      ! User Ancillary fields: A_IXPTR(48)- A_IXPTR(67)
      ! Lateral Boundary Conditions and tendencies:
      ! A_IXPTR(68)- A_IXPTR(94)
      ! tpps_ozone: A_IXPTR(95)
A_SPPTR(A_IXPTR(48)),A_SPPTR(A_IXPTR(49)),A_SPPTR(A_IXPTR(50)),A_SPPTR         &
(A_IXPTR(51)),A_SPPTR(A_IXPTR(52)),A_SPPTR(A_IXPTR(53)),A_SPPTR(A_IXPTR(54)),  &
A_SPPTR(A_IXPTR(55)),A_SPPTR(A_IXPTR(56)),A_SPPTR(A_IXPTR(57)),A_SPPTR         &
(A_IXPTR(59)),A_SPPTR(A_IXPTR(59)),A_SPPTR(A_IXPTR(60)),A_SPPTR(A_IXPTR(61)),  &
A_SPPTR(A_IXPTR(62)),A_SPPTR(A_IXPTR(63)),A_SPPTR(A_IXPTR(64)),A_SPPTR         &
(A_IXPTR(65)),A_SPPTR(A_IXPTR(66)),A_SPPTR(A_IXPTR(67)),A_SPPTR(A_IXPTR(68)),  &
A_SPPTR(A_IXPTR(69)),A_SPPTR(A_IXPTR(70)),A_SPPTR(A_IXPTR(71)),A_SPPTR         &
(A_IXPTR(72)),A_SPPTR(A_IXPTR(73)),A_SPPTR(A_IXPTR(74)),A_SPPTR(A_IXPTR(75)),  &
A_SPPTR(A_IXPTR(76)),A_SPPTR(A_IXPTR(77)),A_SPPTR(A_IXPTR(78)),A_SPPTR         &
(A_IXPTR(79)),A_SPPTR(A_IXPTR(80)),A_SPPTR(A_IXPTR(81)),A_SPPTR(A_IXPTR(82)),  &
A_SPPTR(A_IXPTR(83)),A_SPPTR(A_IXPTR(84)),A_SPPTR(A_IXPTR(85)),A_SPPTR         &
(A_IXPTR(86)),A_SPPTR(A_IXPTR(87)),A_SPPTR(A_IXPTR(88)),A_SPPTR(A_IXPTR(89)),  &
A_SPPTR(A_IXPTR(90)),A_SPPTR(A_IXPTR(91)),A_SPPTR(A_IXPTR(92)),A_SPPTR         &
(A_IXPTR(93)),A_SPPTR(A_IXPTR(94)),A_SPPTR(A_IXPTR(95)),                       &
      ! Biomass aerosol
A_SPPTR(A_IXPTR(96)),A_SPPTR(A_IXPTR(97)),A_SPPTR(A_IXPTR(98)), &
      ! Additional microphysics fields and lbcs
A_SPPTR(A_IXPTR(99)),A_SPPTR(A_IXPTR(100)),A_SPPTR(A_IXPTR(101)), &
A_SPPTR(A_IXPTR(102)),A_SPPTR(A_IXPTR(103)),A_SPPTR(A_IXPTR(104)),&
A_SPPTR(A_IXPTR(105)),A_SPPTR(A_IXPTR(106)),A_SPPTR(A_IXPTR(107)),&
      ! Mineral dust
A_SPPTR(A_IXPTR(108)),A_SPPTR(A_IXPTR(109)),A_SPPTR(A_IXPTR(110)),&
A_SPPTR(A_IXPTR(111)),A_SPPTR(A_IXPTR(112)),A_SPPTR(A_IXPTR(113)),&
      ! Cloud fractions
A_SPPTR(A_IXPTR(114)),A_SPPTR(A_IXPTR(115)),A_SPPTR(A_IXPTR(116)),&
A_SPPTR(A_IXPTR(117)),A_SPPTR(A_IXPTR(118)),A_SPPTR(A_IXPTR(119)),&
! Add pointer for direct PAR flux
A_SPPTR(A_IXPTR(120)),A_SPPTR(A_IXPTR(121)),A_SPPTR(A_IXPTR(122)),&
A_SPPTR(A_IXPTR(123)),A_SPPTR(A_IXPTR(124)),A_SPPTR(A_IXPTR(125)),&
      ! Pointers for UKCA oxidant fields (126-129)
A_SPPTR(A_IXPTR(126)),A_SPPTR(A_IXPTR(127)),A_SPPTR(A_IXPTR(128)),&
A_SPPTR(A_IXPTR(129)),                                            &
      ! Convective Cloud Fields
A_SPPTR(A_IXPTR(130)),A_SPPTR(A_IXPTR(131)),                      &
      ! Ozone tracer and associated cariolle fields
A_SPPTR(A_IXPTR(132)),A_SPPTR(A_IXPTR(133)),A_SPPTR(A_IXPTR(134)),&
A_SPPTR(A_IXPTR(135)),A_SPPTR(A_IXPTR(136)),A_SPPTR(A_IXPTR(137)),&
A_SPPTR(A_IXPTR(138)),A_SPPTR(A_IXPTR(139)),                      &
      ! Pointers for aerosol climatologies (140-160)
A_SPPTR(A_IXPTR(140)),A_SPPTR(A_IXPTR(141)),A_SPPTR(A_IXPTR(142)),&
A_SPPTR(A_IXPTR(143)),A_SPPTR(A_IXPTR(144)),A_SPPTR(A_IXPTR(145)),&
A_SPPTR(A_IXPTR(146)),A_SPPTR(A_IXPTR(147)),A_SPPTR(A_IXPTR(148)),&
A_SPPTR(A_IXPTR(149)),A_SPPTR(A_IXPTR(150)),A_SPPTR(A_IXPTR(151)),&
A_SPPTR(A_IXPTR(152)),A_SPPTR(A_IXPTR(153)),A_SPPTR(A_IXPTR(154)),&
A_SPPTR(A_IXPTR(155)),A_SPPTR(A_IXPTR(156)),A_SPPTR(A_IXPTR(157)),&
A_SPPTR(A_IXPTR(158)),A_SPPTR(A_IXPTR(159)),A_SPPTR(A_IXPTR(160)),&
!Fossil-fuel organic carbon aerosol
A_SPPTR(A_IXPTR(161)),A_SPPTR(A_IXPTR(162)),A_SPPTR(A_IXPTR(163)),&
! CABLE and CASA-CNP + CABLE_LAI 
A_SPPTR(A_IXPTR(164)),A_SPPTR(A_IXPTR(165)),A_SPPTR(A_IXPTR(166)),&
A_SPPTR(A_IXPTR(167)),A_SPPTR(A_IXPTR(168)),A_SPPTR(A_IXPTR(169)),&
A_SPPTR(A_IXPTR(170)),A_SPPTR(A_IXPTR(171)),A_SPPTR(A_IXPTR(172)),&
A_SPPTR(A_IXPTR(173)),A_SPPTR(A_IXPTR(174)),A_SPPTR(A_IXPTR(175)),&
A_SPPTR(A_IXPTR(176)),A_SPPTR(A_IXPTR(177)),A_SPPTR(A_IXPTR(178)),&
A_SPPTR(A_IXPTR(179)),A_SPPTR(A_IXPTR(180)),A_SPPTR(A_IXPTR(181)),&
A_SPPTR(A_IXPTR(182)),A_SPPTR(A_IXPTR(183)),A_SPPTR(A_IXPTR(184)),&
     &        icode,cmessage)
      ENDIF


! DEPENDS ON: iccloseoutput
      CALL ICCLOSEOUTPUT()

!L 4.1 Exit processing: If abnormal completion, output error message
      IABORT = ICODE
      IF (ICODE /= 0) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EREPORT ',3)
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EREPORT ',4)
      END IF
!L 4.2 Exit processing: Perform tidy-up
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITPROC',3)
! DEPENDS ON: exitproc
      CALL EXITPROC(ICODE,CMESSAGE)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EXITPROC',4)
!L 4.3 Exit processing: If error in exit processing, output error mess
      IF (ICODE /= 0) THEN
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EREPORT ',3)
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)
! DEPENDS ON: timer
                                IF (LTIMER) CALL TIMER('EREPORT ',4)
      END IF


!L----------------------------------------------------------------------
!L 5. Complete Timer call and return
!L
      ICODE=IABORT
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('U_MODEL ',4)
      END IF
      RETURN
      END SUBROUTINE U_MODEL
