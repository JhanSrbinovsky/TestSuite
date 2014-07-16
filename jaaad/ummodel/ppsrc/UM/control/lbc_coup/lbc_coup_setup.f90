
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Setup for LBC coupling
Subroutine lbc_coup_setup(ltimer_lbc)

  Implicit None
  

! Description
!  Does the setup in preparation for LBC coupling.
!
! Method
!  The main function of this coupling code is to allow the driving model
! (generating the LBCs) and the LAM model (receiving the LBCs) to run in
! parallel. In a nutshell, this coupling ensures that the LBCs are
! available every time it reaches the stage when it needs to read in more
! LBCs. Communication between the global and LAM was done through a
! communication file - which the global wrote to and the LAM read from it.
!
! This functionality was used for a number of years in the operational
! suite to enable the global and the old UK mesoscale run alongside rather 
! than the UK Mes wait for the global to complete its run. It stopped being 
! used in the oper suite during 2006 when the NAE took over the old
! UK Mes, and the global/NAE main runs are not run in parallel. The
! global/NAE update runs are run in parallel but the NAE uses a LBC file
! from an earlier global run so no need for coupling.
!
! Current Code Owner:
!
! Code Description:
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.


  Logical,Intent(in) :: ltimer_lbc  ! Timer logical

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
  !L       Control for boundary updating.
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

  Integer           :: um_lbc_coup ! LBC Coupling Switch : 1/0 is on/off
  Character(len=80) :: filename    ! Filename of communication file.
  Character(len=8)  :: c_lbc_coup  ! Character variable to read env var
  Character(len=8)  :: ch_date2    ! Date returned from date_and_time
  Character(len=10) :: ch_time2    ! Time returned from date_and_time
  Character(len=5)  :: runid_char  ! RUNID for job
  Character(len=4)  :: runtype_char! Run Type (ie. NRUN, CRUN)
  Integer           :: lbc_ntimes  ! No of BCs in communication file.
  Integer           :: len_wait_tot! Total wait for availability of BCs
  Integer*8         :: isleep      ! Return value from SLEEP
  Integer*8         :: sleep       ! SLEEP function to make job wait
  Integer           :: j           ! Loop counter

  Logical           :: l_exist     !  T : Communication File exists
  Logical           :: l_active    !  T : Output stream active for LBCs.

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
  Integer           :: ICODE       ! =0 normal exit; >0 error exit
  Character(len=256):: Cmessage    ! Error message
  Character*(*)        RoutineName
  Parameter (RoutineName='LBC_SETUP_FILES')

! --------------------- COMDECK LBC_COUP --------------------------
!    Description:
!       This COMDECK stores the variables connected with the
!       Parallel Running between the Global and Mesoscale models.
!       The Mesoscale has to be held back until there are sufficient
!       Boundary Conditions (BCs) to proceed.
!
!   History:
!
!   Model    Date     Modification history
!  version
!   4.5    17/08/98 New comdeck. D. Robinson
!
      logical l_lbc_coup        ! T : Global/Meso Coupling on
      integer um_lbc_stream     ! Output Stream generating BCs.
      
      ! UM 6.5 -  MODEL_ANALYSIS_HRS changed to REAL - 
      !                requires LBC_FC_HRS to REAL also
      real    lbc_fc_hrs        ! Forecast time w.r.t analysis time
      integer um_lbc_wait       ! Wait time between re-tries if BCs
                                ! not available.
      integer um_lbc_wait_max   ! Maximum wait time.
      character*80 lbc_filename ! Name of file that communicates between
                                ! Global and Meso.

      COMMON /LBC_COUP/ l_lbc_coup, um_lbc_stream, lbc_fc_hrs,          &
     &                  um_lbc_wait, um_lbc_wait_max, lbc_filename
! ----------------------------------------------------------------------

  ! Find out if LBC Coupling has been switched on this run
  ! from the env. variable UM_LBC_COUP.

  Call fort_get_env('UM_LBC_COUP',11,c_lbc_coup,8,icode)
  If (icode /= 0) Then
     um_lbc_coup = 0 ! No coupling
     If (PrintStatus  >=  PrStatus_Normal) Then
        Write (6,*) ' Env Var UM_LBC_COUP not set.'
        Write (6,*) ' Setting UM_LBC_COUP to ',um_lbc_coup
     End If
  Else
     Read(c_lbc_coup,'(i8)') um_lbc_coup
     If (PrintStatus  >=  PrStatus_Normal) Then
        Write (6,*) ' UM_LBC_COUP is set to ',UM_LBC_COUP
     End If
  End If
  If (um_lbc_coup == 0 .Or. um_lbc_coup == 1) Then
     l_lbc_coup = um_lbc_coup == 1
  Else
     Write (6,*) ' Invalid value given to UM_LBC_COUP ',             &
          UM_LBC_COUP
     Write (6,*) ' Valid values are 0 or 1'
     Write (6,*) ' L_LBC_COUP set to F. No LBC Coupling ',           &
          'in this run.'
     cmessage = 'U_MODEL : Invalid value given to UM_LBC_COUP'
     icode = 100
! DEPENDS ON: ereport
     Call Ereport(RoutineName,ICODE,Cmessage)
  End If

  If (PrintStatus  >=  PrStatus_Normal) Then
     If (l_lbc_coup) Then
        Write (6,*) ' LBC COUPLING switched on in this run.'
     Else
        Write (6,*) ' LBC COUPLING switched off in this run.'
     End If
  End If

  If (l_lbc_coup) Then


     ! Find out which LBC output stream is providing the data
     ! from the env. variable UM_LBC_STREAM.

     Call fort_get_env('UM_LBC_STREAM',13,c_lbc_coup,8,icode)
     If (icode /= 0) Then
        um_lbc_stream = 0 ! No coupling
        If (PrintStatus  >=  PrStatus_Min) Then
           Write (6,*) ' gl : Env Var UM_LBC_STREAM not set.'
           Write (6,*) ' gl : Setting UM_LBC_STREAM to ',um_lbc_stream
        End If
     Else
        Read(c_lbc_coup,'(i8)') um_lbc_stream
        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' gl : UM_LBC_STREAM is set to ',UM_LBC_STREAM
        End If
     End If

     ! Check validity of UM_LBC_STREAM
     
     If (um_lbc_stream <  1.Or.um_lbc_stream >  max_n_intf_a) Then
        Write (6,*) ' gl : UM_LBC_STREAM = ',UM_LBC_STREAM,           &
             ' is an invalid value.'
        Write (6,*) ' gl : Valid values are 1-',MAX_N_INTF_A
        cmessage = 'U_MODEL : Invalid value given to UM_LBC_STREAM'
        icode = 101
! DEPENDS ON: ereport
        Call Ereport(RoutineName,ICODE,Cmessage)
     End If

     ! Check if this output stream is active.
     l_active = .False.
     Do j=1,n_intf_a
        l_active = l_active .Or. um_lbc_stream == lbc_stream_a(j)
     End Do
     If (.Not.l_active) Then
        Write (6,*) ' gl : Output LBC stream ',UM_LBC_STREAM,         &
             ' is inactive. Check UM_LBC_STREAM.'
        Write (6,*) ' gl : Active LBC streams are ',                  &
             (LBC_STREAM_A(j),j=1,n_intf_a)
        cmessage = 'U_MODEL : Output LBC stream is inactive.'
        icode = 101
! DEPENDS ON: ereport
        Call Ereport(RoutineName,ICODE,Cmessage)
     End If



  End If  !  if l_lbc_coup

  ICODE=0  ! (incase FORT_GET_ENV has left with non-zero value)

  If (l_lbc_coup) Then

     If (mype == 0) Then

        ! Get filename attached to Unit 190
        Call GET_FILE(190,filename,80,ICODE)

        If (icode /= 0) Then
           Write (6,*) ' gl : Problem with GET_FILE',                  &
                ' for Unit No 190.'
           Write (6,*) ' gl : Return code from GET_FILE ',icode
           Write (cmessage,*)                                          &
                'U_MODEL : Error in GET_FILE for Unit No 190.'
           icode = 102
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' gl : Filename for unit no 190 ',FILENAME
        End If

        ! Open the file with WRITE permission.
        Open(unit=190,FILE=filename,action="write",iostat=icode)
        
        If (icode /= 0) Then
           Write (6,*) ' gl : Problem with OPEN for Unit 190.'
           Write (6,*) ' gl : Return code from OPEN ',icode
           Write (cmessage,*)                                          &
                'U_MODEL : Problem with OPEN for Unit No 190.'
           icode = 103
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' gl : File(unit 190) has been opened.'
        End If

        ! Send info to meso that L_LBC_COUP=T in Global.
        
        lbc_ntimes = 1000 + um_lbc_stream
        Write (190,*) lbc_ntimes

! DEPENDS ON: um_fort_flush
        Call um_fort_flush (190,icode)

        If (icode /= 0) Then
           Write (6,*) 'Return Code from FLUSH ',icode
           icode = 104
           Write (cmessage,*) 'U_MODEL : Error flushing out '//        &
                'contents for Unit 190.'
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' gl : lbc_ntimes ', lbc_ntimes,                  &
                ' sent to LBC_FILE.'
        End If

        ! Unit 191 : File with information for operators to
        ! monitor progress with Boundary Conditions generated.
           
        ! Get filename attached to Unit 191
        Call GET_FILE(191,filename,80,ICODE)

        If (icode /= 0) Then
           Write (6,*) ' gl : Problem with GET_FILE for Unit 191.'
           Write (6,*) ' gl : Return code from GET_FILE ',icode
           Write (cmessage,*)                                          &
                'U_MODEL : Error in GET_FILE for Unit No 191.'
           icode = 105
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        If (PrintStatus  >=  PrStatus_Normal) Then
           Write (6,*) ' gl : Filename for unit no 191 ',FILENAME
        End If
           
        ! Open the file with WRITE permission only.
        Open(unit=191,FILE=filename,action="write",iostat=icode)

        If (icode /= 0) Then
           Write (6,*) ' gl : Problem with OPEN for Unit 191.'
           Write (6,*) ' gl : Return code from OPEN ',icode
           Write (cmessage,*)                                          &
                'U_MODEL : Problem with OPEN for Unit No 191.'
           icode = 106
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        Write (6,*) ' gl : File opened on unit 191.'

        ! Send RUNID and date to file on unit 191
        Call fort_get_env('RUNID',5,runid_char,5,icode)
        If (icode /= 0) Then
           Write (6,*) ' Problem with FORT_GET_ENV for RUNID.'
           Write (cmessage,*)                                          &
                'U_MODEL : Problem with FORT_GET_ENV for RUNID.'
           icode = 107
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        Call fort_get_env('TYPE',4,runtype_char,5,icode)
        If (icode /= 0) Then
           Write (6,*) ' Problem with FORT_GET_ENV for TYPE.'
           Write (cmessage,*)                                          &
                'U_MODEL : Problem with FORT_GET_ENV for TYPE.'
           icode = 108
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

        Call Date_and_time(ch_date2, ch_time2)
        Write (191,*) ' RUNID : ',Trim(runid_char),                   &
             '  RUN TYPE : ',Trim(runtype_char),                      &
             '  on ',ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4)
        
! DEPENDS ON: um_fort_flush
        Call um_fort_flush (191,icode)
        
        If (icode /= 0) Then
           Write (6,*) 'Return Code from FLUSH ',icode
           icode = 109
           Write (cmessage,*) 'U_MODEL : Error flushing out '//        &
                'contents for Unit 191.'
! DEPENDS ON: ereport
           Call Ereport(RoutineName,ICODE,Cmessage)
        End If

     End If  !  if mype == 0
        
  End If  !   if l_lbc_coup


  Return
End Subroutine lbc_coup_setup
