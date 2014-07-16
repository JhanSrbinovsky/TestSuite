
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine GRAVSETT ----------------------------------------------
!
! Purpose: To perform gravitational settlement of tracer particles
!          down to the lowest layer of the model.
!          This version allows tracers to fall through 1 or 2 layers.
!
! Current owners of code:                 S Woodward, M Woodage
!
! History:
! Version    Date     Comment
! -------    ----     -------
!   4.4    03/10/97   Original code        S Woodward, M Woodage
!
!   5.5    03/01/03   Modified for New Dynamics
!                     Removed if def 17-A, as required whenever dust
!                     code used             S Woodward
! 6.2      15/12/05  Correct dust diagnostics when substepping phys2.
!                                                        M. Diamantakis
!
! Code description:
!  Language: FORTRAN77 + extensions
!  Programming standard: UMDP 3 Vn 6
!
! System components covered:
!
! System task:
!
!Documentation: Not yet available
!
!-----------------------------------------------------------------------
!
      SUBROUTINE GRAVSETT(                                              &
     & ROW_LENGTH,ROWS,NLEVS,TRACFLD,DIAM,RHOP,P_LAYER_CENTRES,         &
     & P_LAYER_BOUNDARIES,T,TIMESTEP,NUM_SUBSTEPS,SUBSTEP_NUMBER,DRYDEP,&
     & L_CAM_DUST)
!
      IMPLICIT NONE
!
      INTEGER ROW_LENGTH           !IN row length
      INTEGER ROWS                 !IN number of rows
      INTEGER NLEVS                !IN number of model levels
      INTEGER NUM_SUBSTEPS         !IN number of phys2 substeps
      INTEGER SUBSTEP_NUMBER       !IN phys2 substep number
!
      REAL DIAM                    !IN tracer particle diameter
      REAL RHOP                    !IN tracer particle density
      REAL P_LAYER_CENTRES(ROW_LENGTH,ROWS,0:NLEVS)    !IN
      REAL P_LAYER_BOUNDARIES(ROW_LENGTH,ROWS,0:NLEVS) !IN
      REAL T(ROW_LENGTH,ROWS,NLEVS)!IN temperature
      REAL TIMESTEP                !IN timestep s
!
      REAL TRACFLD(ROW_LENGTH,ROWS,NLEVS) !IN/OUT tracer field
!
      REAL DRYDEP(ROW_LENGTH,ROWS) !IN/OUT dep flux from
                                   !       layer2(kg m-2 s-1)
      LOGICAL L_CAM_DUST           !IN Use old version of dust_uplift 
                                   !   scheme for use in CAM NWP models
                                   !   If this is the case, we allow
                                   !   tracers to settle directly from level 1
!

! Include COMDECKS
!

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
!
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL                                                              &
     &     EVAPTAU,                                                     &
                          ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,                                                      &
                          ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,                                                 &
                          ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,                                                  &
                                  ! high pressure reaction rate limit
     &     K_DMS_OH,                                                    &
                                  ! reaction rate for DMS+OH  cc/mcl/s
     &      K4_CH3SO2_O3,                                               &
                             ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
     &      K5_CH3SO3_HO2,                                              &
                             ! Rate coeff for CH3SO3+HO2 -> MSA+O2
     &      RMM_O3,                                                     &
                             ! relative molecular mass O3
     &     BRAT_SO2,                                                    &
                                  ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,                                                    &
                                  ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,                                                    &
                                 ! no. of molecules in 1 mole
     &     RMM_H2O2,                                                    &
                                 ! relative molecular mass H2O2 kg/mole
     &     RMM_AIR,                                                     &
                                 ! relative molecular mass dry air
     &     RMM_W,                                                       &
                                 ! relative molecular mass water
     &     RELM_S_H2O2,                                                 &
                                 ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,                                                   &
                              ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
     &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
     &     T1,                                                          &
                                !
     &     FC,                                                          &
                                ! parameters for interpolation between
     &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
     &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,                                                    &
                                !
     &     CLOUDTAU,                                                    &
                                  ! air parcel lifetime in cloud
     &     CHEMTAU,                                                     &
                                  ! chem lifetime in cloud before oxidn
     &     O3_MIN,                                                      &
                              ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (                                                       &
     &           EVAPTAU = 300.0,                                       &
                                              ! secs  (=5 mins)
     &             NUCTAU = 30.0,                                       &
                                          ! secs
     &       DIFFUSE_AIT = 1.7134E-9,                                   &
                                             ! sq m/s
     &        K_SO2OH_HI = 2.0E-12,                                     &
                                       ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,                                    &
                                          ! cc/mcl/s
     &       K4_CH3SO2_O3 = 1.0E-14,                                    &
                                        ! cc/mcl/s
     &      K5_CH3SO3_HO2 = 4.0E-11,                                    &
     &             RMM_O3 = 4.8E-2,                                     &
                                        ! kg/mole
     &          BRAT_SO2 = 0.9,                                         &
     &           BRAT_MSA = 1.0-BRAT_SO2,                               &
     &           AVOGADRO = 6.022E23,                                   &
                                          ! per mole
     &           RMM_H2O2 = 3.40E-2,                                    &
                                          ! kg/mole
     &            RMM_AIR = 2.896E-2,                                   &
                                          ! kg/mole
     &              RMM_W = 1.8E-2,                                     &
                                          ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,                                 &
     &           RELM_S_2N = 3.206/2.80,                                &
     &               PARH = 3.3,                                        &
     &                K1 = 4.0E-31,                                     &
                                       ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,                                      &
                                          ! K
     &                FC = 0.45,                                        &
                                        ! from STOCHEM model
     &              FAC1 = 1.1904,                                      &
                                    ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.2E-13,                                    &
                                          ! cc/mcl/s
     &                 K3 = 1.9E-33,                                    &
                                          ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,                                    &
                                          ! cc/mcl
     &                 T2 = 600.0,                                      &
                                          ! K
     &                 T3 = 890.0,                                      &
                                          ! K
     &                 T4 = 2200.0,                                     &
                                          ! K
     &           CLOUDTAU = 1.08E4,                                     &
                                          ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,                                      &
                                          ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,                                    &
                                        !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8                                      &
                                          ! kg/kg
     &          )
!
      REAL RAD_AIT,                                                     &
                            ! median radius of Aitken mode particles
     &     DIAM_AIT,                                                    &
                            !   "    diameter    "
     &     RAD_ACC,                                                     &
                            ! median radius of acccumulation mode
     &     DIAM_ACC,                                                    &
                            !   "    diameter    "
     &     CHI,                                                         &
                            ! mole fraction of S in particle
     &     RHO_SO4,                                                     &
                            ! density of  SO4 particle
     &     SIGMA,                                                       &
                            ! standard devn of particle size distn
!                                 for accumulation mode
     &     E_PARM,                                                      &
                            ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (                                                       &
     &           RAD_AIT = 6.5E-9,                                      &
                                             ! m
     &          DIAM_AIT = 2.0*RAD_AIT,                                 &
     &           RAD_ACC = 95.0E-9,                                     &
                                             ! m
     &          DIAM_ACC = 2.0*RAD_ACC,                                 &
     &               CHI = 32.0/132.0,                                  &
     &           RHO_SO4 = 1769.0,                                      &
                                              ! kg/m3
     &             SIGMA = 1.4,                                         &
     &            E_PARM = 0.9398,                                      &
     &          NUM_STAR = 1.0E6                                        &
                                             ! m-3
     &          )
!
      REAL BOLTZMANN       !Boltzmanns constant.
      REAL MFP_REF         !Reference value of mean free path.
      REAL TREF_MFP        !Reference temperature for mean free path.
      REAL PREF_MFP        !Reference pressure for mean free path.
      REAL SIGMA_AIT       !Geometric standard deviation of the Aitken
!                             mode distribution.
!
      PARAMETER (BOLTZMANN = 1.3804E-23)  ! J K-1
      PARAMETER (MFP_REF = 6.6E-8                                       &
                                          ! m
     &        ,  TREF_MFP = 293.15                                      &
                                          ! K
     &        ,  PREF_MFP = 1.01325E5)    ! Pa
      PARAMETER (SIGMA_AIT = 1.30)
!
!*---------------------------------------------------------------------
!
! External subroutines called
      EXTERNAL VGRAV
!
! Local variables
!
      INTEGER K                  !LOC loop counter for levels
      INTEGER J                  !LOC loop counter for points
      INTEGER I                  !LOC loop counter for points
!
      REAL VRHOCTIMESTEP(ROW_LENGTH,ROWS)  !  v*rho*tracer*deltat @lev
      REAL RHOK2(ROW_LENGTH,ROWS)  !  rho(lev+2)
      REAL RHOK1(ROW_LENGTH,ROWS)  !  rho(lev+1)
      REAL RHOK(ROW_LENGTH,ROWS)   ! rho(lev)
      REAL DZK(ROW_LENGTH,ROWS)    ! thickness of layer lev
      REAL DZK1(ROW_LENGTH,ROWS)   ! thickness of layer lev+1
      REAL DZK2(ROW_LENGTH,ROWS)   ! thickness of layer lev+2
      REAL VSTOKES(ROW_LENGTH,ROWS,NLEVS)!deposition velocity
!                                         (vstokes corrected)
      REAL MASSOUT2K2(ROW_LENGTH,ROWS) !flux falling 2 levs from lev k+2
      REAL MASSOUT1K2(ROW_LENGTH,ROWS) !flux falling 1 levs from lev k+2
      REAL MASSOUT2K1(ROW_LENGTH,ROWS) !flux falling 2 levs from lev k+1
      REAL MASSOUT1K1(ROW_LENGTH,ROWS) !flux falling 1 levs from lev k+1
      REAL MASSOUT2K(ROW_LENGTH,ROWS)  !flux falling 2 levs from lev k
      REAL MASSOUT1K(ROW_LENGTH,ROWS)  !flux falling 1 levs from lev k
      REAL DUMMY1(ROW_LENGTH,ROWS,NLEVS) !
      REAL DUMMY2(ROW_LENGTH,ROWS,NLEVS) !
!
!
! Calculate settlement velocity
!
!      CALL VGRAV(PFIELD,NLEVS,DIAM,RHOP,PSTAR,AK,BK,T,V,DUMMY1,DUMMY2,
!     &           FIRST_POINT,LAST_POINT)
! DEPENDS ON: vgrav
       CALL VGRAV(                                                      &
     &  ROW_LENGTH,ROWS,NLEVS,DIAM,RHOP,                                &
     &  P_LAYER_CENTRES(1:ROW_LENGTH,1:ROWS,1:NLEVS),T,                 &
     &  VSTOKES,DUMMY1,DUMMY2)

!
! Calculate new tracer mixing ratios
!
! Initialise deposition flux to zero
      IF ( SUBSTEP_NUMBER == 1 ) THEN
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            DRYDEP(I,J)=0.
          ENDDO
        ENDDO
      END IF
!
! Level 1 (K at start of loop)
!
      IF (L_CAM_DUST) THEN
!       If we are running with the CAM dust scheme, we allow
!       tracer to settle directly from level 1
        DO J = 1,ROWS
          DO I=1,ROW_LENGTH
            RHOK(I,J)=P_LAYER_CENTRES(I,J,1)/(R*T(I,J,1))
            DZK(I,J)=(P_LAYER_BOUNDARIES(I,J,0)-P_LAYER_BOUNDARIES(I,J,1))&
     &              /(RHOK(I,J)*G)
            MASSOUT2K(I,J)=0.
            MASSOUT1K(I,J)=RHOK(I,J)*TRACFLD(I,J,1)*                   &
     &         VSTOKES(I,J,1)*TIMESTEP 
          ENDDO               !I
        ENDDO                !J
      ELSE
!     If we are not running the CAM dust scheme, we do not allow
!     tracer to settle directly from level 1
        DO J = 1,ROWS
          DO I=1,ROW_LENGTH
            RHOK(I,J)=P_LAYER_CENTRES(I,J,1)/(R*T(I,J,1))
            DZK(I,J)=(P_LAYER_BOUNDARIES(I,J,0)-P_LAYER_BOUNDARIES(I,J,1))&
     &              /(RHOK(I,J)*G)
            MASSOUT2K(I,J)=0.
            MASSOUT1K(I,J)=0.
          ENDDO               !I
        ENDDO                !J
      ENDIF !L_CAM_DUST
!
! Level 2 (K+1 at start of loop)
!   NB  deposit tracer direct to ground from lev 2 if V high enough
!
      DO J = 1,ROWS
        DO I = 1,ROW_LENGTH
!
          RHOK1(I,J)=P_LAYER_CENTRES(I,J,2)/(R*T(I,J,2))
          DZK1(I,J)=(P_LAYER_BOUNDARIES(I,J,1)-                         &
     &     P_LAYER_BOUNDARIES(I,J,2))/(RHOK1(I,J)*G)
!
!   check for deposition :
          IF (VSTOKES(I,J,2)*TIMESTEP  >   DZK(I,J)) THEN
!       some tracer deposited onto ground
!
            IF (VSTOKES(I,J,2)*TIMESTEP  >   DZK1(I,J)+DZK(I,J)) THEN
!         all deposited to ground
              MASSOUT2K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*DZK1(I,J)
              MASSOUT1K1(I,J)=0.
            ELSE IF ( VSTOKES(I,J,2)*TIMESTEP  >   DZK1(I,J) ) THEN
!         some deposited to ground, some to layer 1

              MASSOUT2K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*                &
     &         (VSTOKES(I,J,2)*TIMESTEP-DZK(I,J))
              MASSOUT1K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*                &
     &         DZK1(I,J)-MASSOUT2K1(I,J)
            ELSE
!         some deposited to ground, some to layer1, some left in layer2
              MASSOUT2K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*                &
     &         (VSTOKES(I,J,2)*TIMESTEP-DZK(I,J))
              MASSOUT1K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*DZK(I,J)
            ENDIF
!
            DRYDEP(I,J)=DRYDEP(I,J)                                     &
     &                 +MASSOUT2K1(I,J)/(NUM_SUBSTEPS*TIMESTEP)
!
          ELSE
!         only falls into layer 1
            MASSOUT2K1(I,J)=0.
            IF ( VSTOKES(I,J,2)*TIMESTEP  >   DZK1(I,J)) THEN
!           all falls into layer 1
              MASSOUT1K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*DZK1(I,J)
            ELSE
!           some to layer 1 , some left in layer2
              MASSOUT1K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*VSTOKES(I,J,2)  &
     &         *TIMESTEP
            ENDIF
!
          ENDIF
!
        ENDDO                !END I LOOP
      ENDDO                 !END J LOOP
!
! Main loop through levels, from bottom up
!
      DO K = 1,NLEVS-2
!
        DO J = 1,ROWS
          DO I = 1,ROW_LENGTH
!
            RHOK2(I,J)=P_LAYER_CENTRES(I,J,K+2)/(R*T(I,J,K+2))
            DZK2(I,J)=(P_LAYER_BOUNDARIES(I,J,K+1)-                     &
     &       P_LAYER_BOUNDARIES(I,J,K+2))/(RHOK2(I,J)*G)
!
!       Calculate mass of tracer falling between levels
!
!        limit fall to 2 levs
           IF (VSTOKES(I,J,K+2)*TIMESTEP >  (DZK1(I,J)+DZK(I,J)))       &
     &      VSTOKES(I,J,K+2)=(DZK1(I,J)+DZK(I,J))/TIMESTEP
!
!          check how far tracer falls:
           IF ( VSTOKES(I,J,K+2)*TIMESTEP  >   DZK1(I,J) ) THEN
!          it falls through more than 1 layer
             IF ( VSTOKES(I,J,K+2)*TIMESTEP  >                          &
     &        (DZK2(I,J)+DZK1(I,J)) ) THEN
!             all into layer k
               MASSOUT2K2(I,J)=RHOK2(I,J)*TRACFLD(I,J,K+2)*DZK2(I,J)
               MASSOUT1K2(I,J)=0.
             ELSE IF ( VSTOKES(I,J,K+2)*TIMESTEP  >   DZK2(I,J) ) THEN
!            some into k+1, some into k
               MASSOUT2K2(I,J)=                                         &
     &          RHOK2(I,J)*TRACFLD(I,J,K+2)*(VSTOKES(I,J,K+2)*          &
     &          TIMESTEP-DZK1(I,J))
               MASSOUT1K2(I,J)=                                         &
     &          RHOK2(I,J)*TRACFLD(I,J,K+2)*DZK2(I,J)-MASSOUT2K2(I,J)
             ELSE
!            some left in k+2, some into k+1, some into k
               MASSOUT2K2(I,J)=                                         &
     &          RHOK2(I,J)*TRACFLD(I,J,K+2)*(VSTOKES(I,J,K+2)*          &
     &          TIMESTEP-DZK1(I,J))
               MASSOUT1K2(I,J)=RHOK2(I,J)*TRACFLD(I,J,K+2)*DZK1(I,J)
             ENDIF
!
           ELSE
!          falls no more than 1 layer
             MASSOUT2K2(I,J)=0.
             IF (VSTOKES(I,J,K+2)*TIMESTEP  >   DZK2(I,J)) THEN
!            all falls into layer k+1
               MASSOUT1K2(I,J)=RHOK2(I,J)*TRACFLD(I,J,K+2)*DZK2(I,J)
             ELSE
!            some falls into k+1, some left in k+2
               MASSOUT1K2(I,J)=RHOK2(I,J)*TRACFLD(I,J,K+2)*             &
     &          VSTOKES(I,J,K+2)*TIMESTEP
             ENDIF
!
          ENDIF
!
! Update tracer field
!
            TRACFLD(I,J,K)=TRACFLD(I,J,K)+(MASSOUT2K2(I,J)+             &
     &       MASSOUT1K1(I,J)-MASSOUT2K(I,J)-MASSOUT1K(I,J))/            &
     &       (RHOK(I,J)*DZK(I,J))
!
! Put k+2 vals in k+1's & k+1's in k's
            MASSOUT1K(I,J)=MASSOUT1K1(I,J)
            MASSOUT1K1(I,J)=MASSOUT1K2(I,J)
            MASSOUT2K(I,J)=MASSOUT2K1(I,J)
            MASSOUT2K1(I,J)=MASSOUT2K2(I,J)
            DZK(I,J)=DZK1(I,J)
            DZK1(I,J)=DZK2(I,J)
            RHOK(I,J)=RHOK1(I,J)
            RHOK1(I,J)=RHOK2(I,J)
!
          ENDDO           !END I LOOP
        ENDDO            !END J LOOP
!
      ENDDO              !END K LOOP
!
! Top 2 levels
!
      DO J=1,ROWS
        DO I=1,ROW_LENGTH
!
          TRACFLD(I,J,NLEVS-1)=TRACFLD(I,J,NLEVS-1)+                    &
     &     (MASSOUT1K1(I,J)-MASSOUT2K(I,J)-MASSOUT1K(I,J))/             &
     &     (RHOK(I,J)*DZK(I,J))
          TRACFLD(I,J,NLEVS)=TRACFLD(I,J,NLEVS)-                        &
     &     (MASSOUT2K1(I,J)+MASSOUT1K1(I,J))/                           &
     &                (RHOK1(I,J)*DZK1(I,J))
!
        ENDDO       !I
      ENDDO        !J
!
      RETURN
      END SUBROUTINE GRAVSETT
