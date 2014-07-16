
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine SETCONA ------------------------------------------------
!LL
!LL Purpose : to calculate additional constants derived from
!LL      information in the input dump and pass out into the model via
!LL      the argument lists in argcona and arglndm.
!LL
!LL Method:
!LL    0. Initialise vertical coordinates from eta values and surface
!LL       orography.
!LL    0.1 Initialise data arrays held in secondary dump space.
!LL    1. Trigonometric functions:
!LL        Convert from degrees to radians.
!LL        Calculate trig functions for this grid.
!LL        Update halos for trig functions.
!LL    2. Set up semi-lagrangian advection options.
!LL       Initialise land/soil points in index arrays.
!LL    3. Determine model level boundaries for L/M/H cloud diagnostics.
!LL    4. Add locally derived parameters.
!LL
!LL Language: FORTRAN 77 + common extensions also in Fortran 90.
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL version 7.2, dated 5/2/98
!LL
!LL Documentation : Unified Model Documentation Paper No P0
!LL
!LLEND----------------------------------------------------------------
!
!*L arguments

      SUBROUTINE Setcona                                                &
! Input data (including some logical fields)
     &     (eta_theta_in,eta_rho_in,Smvcst,Land_sea_mask,Orog,          &
     &      grad_x, grad_y,                                             &
     &      rho,exner_rho_levels,                                       &
     &      orog_lbc,exner_lbc,                                         &
! Input size and control variables
     &      LENRIMA,LBC_SIZEA,LBC_STARTA,                               &
     &      RIMWIDTHA, RIMWEIGHTSA,                                     &
     &      global_row_length, global_rows,                             &
     &      MODEL_LEVELS,ROWS,N_ROWS,ROW_LENGTH,                        &
     &      LAND_FIELD,WET_LEVELS,boundary_layer_levels,                &
     &      first_constant_r_rho_level,cloud_levels,z_top_of_model,     &
     &       tr_levels, tr_vars, tr_ukca,                               &
     &      height_gen_method,                                          &
! other constants
     &      Delta_lambda_in,Delta_phi_in,Base_phi_in,Base_lambda_in,    &
     &      lat_rot_NP_in,long_rot_NP_in,                               &
! VarRes grid info in degrees
     &      Lambda_p_in, Lambda_u_in, Phi_p_in, Phi_v_in,               &
     &      RowDepCStart,                                               &
! Initialise variables in secondary D1:
     &      exner_theta_levels,p_theta_levels,p,p_star,                 &
! Output data
! ARGCONA start
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add new field sin_u_latitude. J F Thomson.
!  5.0   02/06/99  Insert C-P C-grid constants. M L Gallani.
!  5.3   01/10/01  Add fields for chequerboard radiation. S Cusack
! 6.1  04/08/04  Add separate arrays diff_coeff_u, diff_coeff_v
!                                                     Terry Davies
! 6.2  05/01/06   Add true_latitude. Yongming Tang
! 6.2  14/12/06  Add separate arrays VarRes Array co-ordinates
!                                                       Terry Davies

! argcona.h contained constants for the atmosphere.
! As of vn6.6 these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, AD_MASK_TROP_MOD, ROT_COEFF_MOD

! ARGLNDM Constants for physics routines
     &  land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
     &     icode,Cmessage,isSTASH)  ! FLUME-STASH

      Use solinc_data, Only: L_orog, slope_angle
Use level_heights_Mod
Use trignometric_Mod
Use dyn_coriolis_Mod
Use dyn_var_res_Mod
Use diff_coeff_Mod
Use rad_mask_trop_Mod
Use rot_coeff_Mod
Use volcts_Mod
Use scvary_Mod

      IMPLICIT NONE

!     INCLUDED COMDECKS
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
! ----------------------- Header file CRUNTIMC  -----------------------
! Description: Run-time constants for the Atmosphere model (read only).
!              Contains variables that define parametrization values
!              chosen for atmosphere physics and dynamics schemes.
!              [Note that CNTLATM holds accompanying control switches
!              needed for addressing.]
!
! Current Code Owner: R. Rawlins
!
!------------------   Physics:   --------------------------------------
! Generalised physics switches:


! Boundary layer:

      LOGICAL :: l_use_bl_diag_term
      LOGICAL :: L_SBLeq
      LOGICAL :: L_SBLco
      LOGICAL :: L_LAMBDAM2    ! LambdaM=2*LambdaH (operational setting)
      LOGICAL :: L_FULL_LAMBDAS! Lambdas NOT reduced above NTML_LOCAL+1
      LOGICAL :: L_us_blsol    ! Switch for stable and non-oscillatory
                               !   BL vertical diffusion scheme
      REAL :: alpha_Cd (max_number_alpha_cds)
      REAL :: Muw_SBL
      REAL :: Mwt_SBL
      REAL :: Charnock
      REAL :: Z0HSEA           ! Roughness lengths for heat and moisture
      REAL :: Z0MIZ            !   over sea, marginal ice zone and
      REAL :: Z0SICE           !   sea ice (m)
      REAL :: ALBSNC_NVG(NNVG) ! Snow-covered albedo.
      REAL :: ALBSNF_NVG(NNVG) ! Snow-free albedo.
      REAL :: CATCH_NVG(NNVG)  ! Canopy capacity (kg/m2).
      REAL :: GS_NVG(NNVG)     ! Surface conductance (m/s).
      REAL :: INFIL_NVG(NNVG)  ! Infiltration enhancement factor.
      REAL :: ROOTD_NVG(NNVG)  ! Rootdepth (m).
      REAL :: Z0_NVG(NNVG)     ! Roughness length (m).
      REAL :: CH_NVG(NNVG)     ! "Canopy" heat capacity (J/K/m2)
      REAL :: VF_NVG(NNVG)     ! Fractional "canopy" coverage
      REAL :: Puns, Pstb ! parameters for uncond stable numerical solver
                         ! Puns : used in an unstable BL column
                         ! Pstb : used in an stable BL column
      REAL :: OROG_DRAG_PARAM  ! Drag coefficient for orographic form drag
      LOGICAL :: L_EMIS_LAND_GEN ! Aggregate land surface emissivity
      REAL :: EMIS_LAND_GEN      ! Aggregate land surface emissivity
      REAL :: WeightLouisToLong_in
!                        ! Weighting of Louis tails towards long
!                        ! tails: this should be removed in favour
!                        ! of placing RUN_BL in the module bl_option_mod 
!                        ! at a later release.

      INTEGER, DIMENSION(20) :: BL_OPTIONS
      INTEGER :: DECFIX
      INTEGER :: STOPWE_SBL
      INTEGER :: TRWEIGHTS1
      INTEGER :: FLUX_GRAD
      INTEGER :: NG_STRESS
      INTEGER :: SBL_OP
      INTEGER :: ISHEAR_BL
      INTEGER :: FORMDRAG, COR_UST, COR_MO_ITER
      INTEGER :: NON_LOCAL_BL
      INTEGER :: LOCAL_FA
      INTEGER :: PRANDTL
      INTEGER :: I_SCRN_T_DIAG
      INTEGER :: Buddy_sea
      INTEGER :: FD_stab_dep
      INTEGER :: Keep_Ri_FA
      INTEGER :: NL_BL_LEVELS

      INTEGER :: CAN_RAD_MOD              ! For radiation model (Pete Falloon)
      INTEGER :: ILAYERS                  ! For radiation model (Pete Falloon)

      COMMON  /RUN_BLVEG/ALBSNC_NVG,ALBSNF_NVG,CATCH_NVG,GS_NVG,        &
     &  INFIL_NVG,ROOTD_NVG,Z0_NVG,CH_NVG,VF_NVG,                       &
     &  CAN_RAD_MOD,ILAYERS
       NAMELIST/RUN_BLVEG/ALBSNC_NVG,ALBSNF_NVG,CATCH_NVG,GS_NVG,       &
     &  INFIL_NVG,ROOTD_NVG,Z0_NVG,CH_NVG,VF_NVG,                       &
     &  CAN_RAD_MOD,ILAYERS

      COMMON  /RUN_BLICE/Z0HSEA,Z0MIZ,Z0SICE
      NAMELIST /RUN_BLICE/Z0HSEA,Z0MIZ,Z0SICE
      REAL :: SeaSalinityFactor
!       Scaling of qsat allowing for salinity of sea water
      INTEGER :: ISeaZ0T
!       Option for specifying the thermal roughness length
!       at sea points
      INTEGER :: ISeaDynDiag
!       Option for dynamic diagnosis of boundary layer types
       ! Maximum and minimum values for the STPH_RP scheme
       ! Boundary Layer
       REAL         :: par_mezcla_max,par_mezcla,par_mezcla_min ! Max, 
                      ! mean and min value for the neutral mixing length
       REAL         :: G0_max,G0_RP,G0_min ! Max,mean and min values
                      ! for the flux profile parameter
       REAL         :: Charnock_max,Charnock_min ! Max and min values for
                      ! the charnock parameter

      COMMON  /RUN_BL/ FORMDRAG, OROG_DRAG_PARAM, l_use_bl_diag_term         &
     &  ,L_SBLeq, L_SBLco, alpha_Cd, Muw_SBL, Mwt_SBL, BL_OPTIONS            &
     &  ,par_mezcla_max, par_mezcla, par_mezcla_min                          &
     &  ,G0_max, G0_RP, G0_min, Charnock_max, Charnock, Charnock_min         &
     &  ,L_LAMBDAM2, L_FULL_LAMBDAS, SeaSalinityFactor, ISeaZ0T              &
     &  ,ISeaDynDiag ,L_us_blsol, Puns, Pstb, L_EMIS_LAND_GEN, EMIS_LAND_GEN &
     &  ,WeightLouisToLong_in

      NAMELIST/RUN_BL/ FORMDRAG, OROG_DRAG_PARAM                             &
     &  ,l_use_bl_diag_term, SBL_OP, COR_UST,  COR_MO_ITER                   &
     &  ,L_SBLeq, L_SBLco, alpha_Cd, Muw_SBL, Mwt_SBL, ISHEAR_BL, NG_STRESS  &
     &  ,DECFIX, STOPWE_SBL, TRWEIGHTS1, FLUX_GRAD                           &
     &  ,NON_LOCAL_BL, NL_BL_LEVELS                                          &
     &  ,par_mezcla_max, par_mezcla, par_mezcla_min                          &
     &  ,G0_max, G0_RP, G0_min, Charnock_max, Charnock, Charnock_min         &
     &  ,L_LAMBDAM2, L_FULL_LAMBDAS, LOCAL_FA, Prandtl, SeaSalinityFactor    &
     &  ,ISeaZ0T, ISeaDynDiag, Buddy_sea, L_us_blsol, Puns, Pstb             &
     &  ,FD_stab_dep, Keep_Ri_FA                                             &
     &  ,L_EMIS_LAND_GEN, EMIS_LAND_GEN, I_SCRN_T_DIAG                       &
     &  ,WeightLouisToLong_in

! Surface parameters for plant functional types

      REAL    :: ALBSNC_MAX(NPFT)
      REAL    :: ALBSNC_MIN(NPFT)
      REAL    :: ALBSNF_MAX(NPFT)
      REAL    :: DZ0V_DH(NPFT)
      REAL    :: CATCH0(NPFT)
      REAL    :: DCATCH_DLAI(NPFT)
      REAL    :: INFIL_F(NPFT)
      REAL    :: KEXT(NPFT)
      REAL    :: ROOTD_FT(NPFT)

      COMMON  /RUN_PFT/ALBSNC_MAX,ALBSNC_MIN,ALBSNF_MAX,DZ0V_DH,        &
     &  CATCH0,DCATCH_DLAI,INFIL_F,KEXT,ROOTD_FT
      NAMELIST/RUN_PFT/ALBSNC_MAX,ALBSNC_MIN,ALBSNF_MAX,DZ0V_DH,        &
     &  CATCH0,DCATCH_DLAI,INFIL_F,KEXT,ROOTD_FT


      ! Large scale precipitation:

      ! Threshold cloud liquid water content over land/sea for
      ! conversion to ppn (kg water per m**3)

      ! Defunct: not supported
      REAL,PARAMETER:: CW_LAND = 8.0E-4
      REAL,PARAMETER:: CW_SEA  = 2.0E-4


      LOGICAL :: L_cry_agg_dep ! Limit deposition to moisture available
      LOGICAL :: L_it_melting  ! Use Iterative melting
      LOGICAL :: L_seq_mcr     ! Use sequential microphysics updating
      LOGICAL :: L_autoc_3b    ! Use 3B autoconversion method
      LOGICAL :: L_autolim_3b  ! Use 3B autoconversion limit
      LOGICAL :: L_autoconv_murk  ! Use murk aerosol to calculate
                                  ! droplet number
      LOGICAL :: L_droplet_settle ! Allow cloud droplets to settle
      LOGICAL :: L_psd         ! Use generic ice particle size distn.

      REAL :: ec_auto          ! Collision / collection efficiency
      REAL :: N_drop_land      ! Droplet concentration over land
      REAL :: N_drop_sea       ! Droplet concentration over sea
      REAL :: Ntot_land        ! Droplet concentration over land
      REAL :: Ntot_sea         ! Droplet concentration over sea
      REAL :: N_drop_land_cr   ! N_drop_land ^ (-1/3)
      REAL :: N_drop_sea_cr    ! N_drop_sea ^ (-1/3)
      REAL :: X1R              ! Intercept parameter for raindrop
                               ! size distribution
      REAL :: X2R              ! Scaling parameter for rain DSD.
      REAL :: X4R              ! Shape factor for rain DSD.
      REAL :: X1I              ! Intercept parameter for aggregate
                               ! size distribution
      REAL :: X1IC             ! Intercept parameter for crystal
                               ! size distribution
      REAL :: AI,  BI          ! Aggregate mass-size params m(D)=AI D^BI
      REAL :: AIC, BIC         ! Crystal mass-size params m(D)=AIC D^BIC
      REAL :: LSP_EI, LSP_FI   ! Aggregate Best-Reynolds relationship
                               ! Re = LSP_EI BE^LSP_FI
      REAL :: LSP_EIC, LSP_FIC ! Crystal Best-Reynolds relationship
                               ! Re = LSP_EIC BE^LSP_FIC

      ! Critical RH for layer cloud formation
      REAL :: RHCRIT(max_model_levels)
       ! Maximum and minimum values for the STPH_RP scheme
       ! Large Scale Precipitation
      REAL         :: RHCRIT_max,RHCRIT_min ! Max and min values
                       ! for the critical relative humidity
      REAL         :: CI_max,CI_min ! Max and min values for the
                       ! CI parameter, which controls the ice fall speed
      REAL         :: M_CI_max,M_CI,M_CI_min ! Max, mean, min values for
                       ! the multiplication factor for CI and CIC in 3C
                       ! Microphysics
 
      COMMON  /RUN_Precip/RHCRIT                                        &
     & ,l_cry_agg_dep,L_it_melting,L_seq_mcr,L_autoc_3b,L_autolim_3b,   &
     &  l_psd,                                                          &
     &  l_autoconv_murk,ec_auto,N_drop_land,N_drop_sea,Ntot_land,       &
     &  Ntot_sea,N_drop_land_cr,N_drop_sea_cr,X1R,X2R,X4R,X1I,X1IC      &
     &  ,AI,BI,AIC,BIC,LSP_EI,LSP_FI,LSP_EIC,LSP_FIC                    &
     &  ,RHCRIT_max,RHCRIT_min                                          &
     &  ,CI_max,CI_min,M_CI_max,M_CI,M_CI_min,L_droplet_settle
      NAMELIST/RUN_Precip/RHCRIT                                        &
     & ,l_cry_agg_dep,L_it_melting,L_seq_mcr,L_autoc_3b,L_autolim_3b,   &
     &  l_psd,                                                          &
     &  l_autoconv_murk,ec_auto,N_drop_land,N_drop_sea,Ntot_land,       &
     &  Ntot_sea,N_drop_land_cr,N_drop_sea_cr,X1R,X2R,X4R,X1I,X1IC      &
     &  ,AI,BI,AIC,BIC,LSP_EI,LSP_FI,LSP_EIC,LSP_FIC                    &
     &  ,RHCRIT_max,RHCRIT_min                                          &
     &  ,CI_max,CI_min,M_CI_max,M_CI,M_CI_min,L_droplet_settle

      ! Large scale cloud:

      LOGICAL :: L_eacf                ! Use empirically adjusted
                                       ! cloud fraction

      INTEGER :: cloud_fraction_method ! Selects total cloud fraction
                                       ! calculation method
      INTEGER :: ice_fraction_method   ! Selects ice cloud fraction
                                       ! calculation method

      REAL    :: overlap_ice_liquid    ! Generic overlap parameter
                                       ! between ice and liquid phases
      REAL    :: ctt_weight            ! Cloud top temperature weight
      REAL    :: t_weight              ! Local temperature weight
      REAL    :: qsat_fixed            ! Prescribed qsat value
      REAL    :: sub_cld               ! Scaling factor
      REAL    :: dbsdtbs_turb_0        ! PC2 erosion rate / s-1

      COMMON  /RUN_Cloud/L_eacf,cloud_fraction_method,                  &
     &  overlap_ice_liquid,ice_fraction_method,ctt_weight,t_weight,     &
     &  qsat_fixed,sub_cld,dbsdtbs_turb_0

      NAMELIST/RUN_Cloud/L_eacf,cloud_fraction_method,                  &
     &  overlap_ice_liquid,ice_fraction_method,ctt_weight,t_weight,     &
     &  qsat_fixed,sub_cld,dbsdtbs_turb_0


      ! Radiation:

      ! True if climatological aerosol is included.
      LOGICAL :: L_climat_aerosol


      ! True to use real boundary layer heights to specify the boundary
      ! layer aerosol.
      LOGICAL :: L_clim_aero_hgt
      LOGICAL :: L_HadGEM1_Clim_Aero
!             Flag to use HadGEM1 setting for climatological aerosols

      LOGICAL :: L_SEC_VAR ! true if using time varying astronomy
      LOGICAL :: L_EqT     ! True if including the equation of time

      INTEGER :: A_SW_SEGMENTS ! No of batches used in shortwave code
      INTEGER :: A_SW_SEG_SIZE ! Size of sw batches. 
      INTEGER :: A_LW_SEGMENTS ! No of batches used in longwave code
      INTEGER :: A_LW_SEG_SIZE ! Size of lw batches.
      INTEGER :: AERO_BL_LEVELS !Common number of layers taken to be 
!                                occupied by the boundary-layer
!                                aerosol if the boundary layer
!                                depth is not used to determine the 
!                                number separately at each grid-point
!                                In previous versions of the code,
!                                this was taken to be BL_LEVELS 

      Logical :: l_ovrlap    ! Requires l_ccrad=.TRUE.
                             ! Allows Convective and LS Cloud to overlap 
                             ! for radiative impacts.
                             ! (Experimental, defaulted to FALSE,
                             ! requires a hand-edit to change)

      Logical :: l_ccw_scav  ! Requires l_ccrad=.TRUE. .AND. l_ovrlap=.TRUE.
                             ! Allows Convective Cloud Water (CCW) to
                             ! compensate for LS Cloud water in overlapping
                             ! LS/CCA fractions.
                             ! (Experimental, defaulted to FALSE, requires a
                             ! hand-edit to change)

      Logical :: l_cldtop_t_fix
                              ! .TRUE. Corrects Cloud Top Temperature in
                              !        radiation

      ! Valid options are specified in the include file o3intp.h.

      ! radiation block is already a separate module for some variables:
!
! Description:
!   Common Block for radiation
!   Contains:
!   Parameters used to increase the true surface albedo linearly from
!   a minimum value, alphaM at melting point, to a "cold-ice" value,
!   alphaC, over a temperature range, dTice
!   Mass Mixing Ratios of minor Gases N2O,CH4,CFC11,CFC12,O2
!   CFC113,HCFC22,HFC125, and HFC134A
!
! Current Code Owner: S Woodward
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.4      09/09/94 Original code. C Hewitt
! 3.4      02/10/94 Taken over, and MMRs added  S Woodward
! 4.4      18/09/97 Mass Mixing Ratio of Oxygen added
! 4.5      21/08/98 ALPHAB introduced and ALPHAM,ALPHAC reinterpreted
!          if l_ssice_albedo i.e. if snow on sea-ice affects albedo.
!          Jonathan Gregory
! 4.5      18/05/98 Mass Mixing Ratios of further halocarbons added.
!                                    J. M. Edwards
! 5.5      11/02/03 Added IS_NCOL    K.D.Williams
! 6.2      06/01/06 sea ice albedo parameters added  J.Ridley
!
! Declarations:
! Imported global variables (*CALLed COMDECKs etc...)

! Global parameters:
! Albedo of sea-ice at melting point (TM) if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at melting point (TM) if l_ssice_albedo
      REAL ALPHAM  ! "M" for "melting"
! Albedo of sea-ice at and below TM-DTICE if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at and below TM-DTICE if l_ssice_albedo
      REAL ALPHAC  ! "C" for "cold"
! Albedo of snow-free sea-ice if l_ssice_albedo
      REAL ALPHAB  ! "B" for "bare"
! Temperature range in which albedo of sea-ice, if .not.l_ssice_albedo,
! or of snow on sea-ice, if l_ssice_albedo, varies between its limits
      REAL DTICE
! Temperature range below TM over which meltponds form if
! l_sice_meltponds and l_ssice_albedo
      REAL DT_BARE
! Increment to albedo for each degree temperature rises above
! TM-DT_BARE. Only used if l_sice_meltponds and l_ssice_albedo
      REAL DALB_BARE_WET
! Fraction of SW radiation that penetrates seaice and scatters
! back causing an addition to the albedo. Only active if l_ssice_albedo
! and l_sice_scattering
      REAL PEN_RAD_FRAC
! attenutation parameter for SW in seaice which controls the
! additional albedo due to internal scattering
      REAL SW_BETA
      REAL N2OMMR  ! N2O mmr (Mass Mixing Ratio)
      REAL CH4MMR  ! CH4 mmr
      REAL C11MMR  ! CFC11 mmr
      REAL C12MMR  ! CFC12 mmr
      REAL O2MMR   ! O2 mmr
      REAL C113MMR      ! CFC113 mmr
      REAL HCFC22MMR    ! HCFC22 mmr
      REAL HFC125MMR    ! HFC125 mmr
      REAL HFC134AMMR   ! HFC134A mmr
! Number of columns used for internal sampling by the ISCCP simulator
      INTEGER IS_NCOL

! Global scalars:

! Global dynamic arrays:

! COMMON blocks:
      COMMON /RAD_COM/ ALPHAM, ALPHAC, ALPHAB, DTICE, N2OMMR, CH4MMR,   &
     & C11MMR, C12MMR, O2MMR, C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR,&
     & IS_NCOL, DT_BARE, DALB_BARE_WET, PEN_RAD_FRAC, SW_BETA

!- End of COMDECK declaration

      REAL :: CO2_MMR                ! CO2 concentration (if constant)

! DIMFIX3A defines internal dimensions tied to algorithms for
! two-stream radiation code, mostly for clouds

      ! number of components of clouds
      INTEGER,PARAMETER:: NPD_CLOUD_COMPONENT=4

      ! number of permitted types of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_TYPE=4

      ! number of permitted representations of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_REPRESENTATION=4

      ! number of overlap coefficients for clouds
      INTEGER,PARAMETER:: NPD_OVERLAP_COEFF=18

      ! number of coefficients for two-stream sources
      INTEGER,PARAMETER:: NPD_SOURCE_COEFF=2

      ! number of regions in a layer
      INTEGER,PARAMETER:: NPD_REGION=3

! DIMFIX3A end
      ! Scaling factors to simulate inhomogeneous cloud.
      REAL, DIMENSION(NPD_CLOUD_COMPONENT) :: INHOM_CLOUD_SW
      REAL, DIMENSION(NPD_CLOUD_COMPONENT) :: INHOM_CLOUD_LW


! sza: add fraction standard diviation of inhomogeneous cloud
!      parameter for triple_clouds
!
      Real FW_STD

      ! Decorrelation pressure scale for large scale cloud
      REAL :: DP_CORR_STRAT
      ! Decorrelation pressure scale for convective cloud
      REAL :: DP_CORR_CONV

      COMMON  /RUN_Radiation/L_climat_aerosol, L_clim_aero_hgt,         &
     &  L_HadGEM1_Clim_Aero, FW_STD,                                    &
     &  A_SW_SEGMENTS,A_SW_SEG_SIZE,A_LW_SEGMENTS,A_LW_SEG_SIZE,CO2_MMR,&
     &  L_SEC_VAR,L_EqT,INHOM_CLOUD_SW,INHOM_CLOUD_LW,DP_CORR_STRAT,    &
     &  DP_CORR_CONV,AERO_BL_LEVELS, l_ovrlap, l_ccw_scav, l_cldtop_t_fix

      NAMELIST/RUN_Radiation/L_climat_aerosol, L_clim_aero_hgt,         &
     &  L_HadGEM1_Clim_Aero, FW_STD,                                    &
     &  A_SW_SEGMENTS,A_SW_SEG_SIZE,A_LW_SEGMENTS,A_LW_SEG_SIZE,CO2_MMR,&
     &  L_SEC_VAR,L_EqT,INHOM_CLOUD_SW,INHOM_CLOUD_LW,DP_CORR_STRAT,    &
     &  DP_CORR_CONV, AERO_BL_LEVELS,l_ovrlap,l_ccw_scav,l_cldtop_t_fix,&
      ! list of variables declared in <rad_com/rad_com.h> to be added:
     &  ALPHAM, ALPHAC, ALPHAB, DTICE,                                  &
     &  DT_BARE,DALB_BARE_WET,PEN_RAD_FRAC,SW_BETA,                     &
     &  N2OMMR, CH4MMR, C11MMR, C12MMR,                                 &
     &  O2MMR, C113MMR, HCFC22MMR, HFC125MMR, HFC134AMMR, IS_NCOL

      ! Gravity wave drag:

      REAL :: KAY_GWAVE       ! Surface stress constant for GW drag

      REAL :: GWD_FRC         ! Critical Froude Number for 4A Scheme
       ! Maximum and minimum values for the STPH_RP scheme
       ! Gravity Wave drag
       REAL         :: GWD_FRC_max, GWD_FRC_min ! Max and min values
                       ! for the critical froud number
       REAL         :: KAY_GWAVE_max, KAY_GWAVE_min ! Max and min values
                       ! for the gravity wave constant

      REAL :: GWD_FSAT       ! Critical Froude number for wave breaking
                             ! used in the amplitude based saturation
                             ! test

      INTEGER :: SAT_SCHEME  ! Switch to determine whether to use
                             ! stress based    (sat_scheme=0) or
                             ! amplitude based (sat_scheme=1)
                             ! saturation test

      LOGICAL :: L_TAUS_SCALE  ! FROUDE NO. DEPENDENT SURFACE STRESS
      LOGICAL :: L_FIX_GWSATN  ! BUG FIXES TO GWSATN
      LOGICAL :: L_GWD_40km    ! Turn off orographic GWD above 40km
      LOGICAL :: L_USSP_OPAQUE ! OPAQUE LID FOR USSP SCHEME

      COMMON  /RUN_GWD/KAY_GWAVE, GWD_FRC, GWD_FSAT                     &
                      ,L_TAUS_SCALE, L_FIX_GWSATN, l_gwd_40km           &
                      ,L_USSP_OPAQUE, SAT_SCHEME                        &
                      ,GWD_FRC_max, GWD_FRC_min                         &
                      ,KAY_GWAVE_max, KAY_GWAVE_min

      NAMELIST/RUN_GWD/KAY_GWAVE, GWD_FRC, GWD_FSAT                     &
                      ,L_TAUS_SCALE, L_FIX_GWSATN, l_gwd_40km           &
                      ,L_USSP_OPAQUE, SAT_SCHEME                        &
                      ,GWD_FRC_max, GWD_FRC_min                         &
                      ,KAY_GWAVE_max, KAY_GWAVE_min

!------------------   End of Physics   ---------------------------------

      ! River Routing
      REAL :: RIVER_VEL       ! River velocity
      REAL :: RIVER_MCOEF     ! Meandering coefficient

      COMMON /RUN_RIVERS/ RIVER_VEL, RIVER_MCOEF
      NAMELIST /RUN_RIVERS/ RIVER_VEL, RIVER_MCOEF


      ! Aerosol Modelling

      INTEGER :: SO2_HIGH_LEVEL   ! Model level for chimney SO2 emiss
      INTEGER :: SOOT_HIGH_LEVEL  ! Model level for chimney soot emiss
      INTEGER :: BMASS_HIGH_LEVEL_1  ! Lowest and highest model levels
      INTEGER :: BMASS_HIGH_LEVEL_2  ! with biomass emissions.
      INTEGER :: OCFF_HIGH_LEVEL  ! Model level for chimney OCFF emiss

      ! improved nonhydrostatic  weights in tracer1 calculation
      LOGICAL :: L_tracer1_non_hydro

      COMMON  /RUN_Aerosol/                                             &
     & SO2_HIGH_LEVEL, SOOT_HIGH_LEVEL, BMASS_HIGH_LEVEL_1,             &
     & BMASS_HIGH_LEVEL_2, OCFF_HIGH_LEVEL, L_tracer1_non_hydro
      NAMELIST/RUN_Aerosol/                                             &
     & SO2_HIGH_LEVEL, SOOT_HIGH_LEVEL, BMASS_HIGH_LEVEL_1,             &
     & BMASS_HIGH_LEVEL_2, OCFF_HIGH_LEVEL, L_tracer1_non_hydro

!     Declarations for UKCA sub-model

      LOGICAL :: L_ukca_chem      ! True when UKCA chemistry is on
      LOGICAL :: L_ukca_family    ! True when using family chemistry
      LOGICAL :: L_ukca_advh2o    ! True when advecting H2O tracer
      LOGICAL :: L_ukca_phot2d    ! True when using 2D photolysis
      LOGICAL :: L_ukca_fastj     ! True when using Fastj photolysis
      LOGICAL :: L_ukca_trop      ! True for tropospheric chemistry
      LOGICAL :: L_ukca_tropisop  ! True for trop chemistry + isoprene
      LOGICAL :: L_ukca_strat     ! True for strat+reduced trop chemistry
      LOGICAL :: L_ukca_strattrop ! True for std strat+trop chemistry
      LOGICAL :: L_ukca_stratcfc  ! True for extended strat chemistry
      LOGICAL :: L_ukca_wachem    ! True for whole atm chemistry
      LOGICAL :: L_ukca_aerchem   ! True for trop+aerosol chemistry
      LOGICAL :: L_ukca_user      ! True for user-defined chemistry
      LOGICAL :: L_ukca_mode      ! True for UKCA-MODE aerosol scheme
      LOGICAL :: L_ukca_dust      ! True for UKCA dust scheme
      LOGICAL :: L_ukca_rnpb      ! True for Radon/Lead scheme
      LOGICAL :: L_ukca_o3budget  ! True for UKCA o3 budget
      LOGICAL :: L_ukca_qch4inter ! True for interact wetland CH4 ems
      LOGICAL :: L_ukca_isopinter ! True for interact isoprene ems
      LOGICAL :: L_ukca_terpinter ! True for interact terpene ems
      LOGICAL :: L_ukca_budget2   ! True for UKCA budget2 calculations
      LOGICAL :: L_ukca_qf11f12mbr  ! True for CFC emissions
      LOGICAL :: L_ukca_useumuivals ! True when using UMUI CFC values
      LOGICAL :: L_ukca_nat_sedi  ! True for NAT sedimentation
      LOGICAL :: L_ukca_het_psc   ! True for Het/PSC chemistry
      LOGICAL :: L_ukca_h2o_feedback ! True for H2O feedback from chem
      LOGICAL :: L_ukca_clbrcons  ! True for Cl/Br conservation
      LOGICAL :: L_ukca_usero3    ! True when using RO3 in photolysis
      LOGICAL :: L_ukca_useco3    ! True when using clim O3 in photolysis
      LOGICAL :: L_ukca_userelaxo3 ! True when using clim O3
!                                    above 0.3hPa in photolysis
      LOGICAL :: L_ukca_rado3   ! T when using UKCA O3 in radiation
      LOGICAL :: L_ukca_radch4  ! T when using UKCA CH4 in radiation
      LOGICAL :: L_ukca_radn2o  ! T when using UKCA N2O in radiation
      LOGICAL :: L_ukca_radf11  ! T when using UKCA CFC-11 in radn
      LOGICAL :: L_ukca_radf12  ! T when using UKCA CFC-12 in radn
      LOGICAL :: L_ukca_radf113 ! T when using UKCA CFC-113 in radn
      LOGICAL :: L_ukca_radf22  ! T when using UKCA HCFC-22 in radn
      LOGICAL :: L_ukca_radch2o ! T when using UKCA adv H2O in radn
      LOGICAL :: L_ukca_useoxid ! T when using UKCA oxidants in aer chem
      LOGICAL :: L_ukca_intdd   ! T when using interact dry deposition

! Tracers and chemistry integers:
      Integer :: jpctr          ! No. of chemical tracers
      Integer :: jpspec         ! No. of chemical species
      Integer :: jpbk           ! No. of bimolecular reactions
      Integer :: jptk           ! No. of termolecular reactions
      Integer :: jppj           ! No. of photolytic reactions
      Integer :: jphk           ! No. of heterogonous reactions
      Integer :: jpnr           ! jpbk + jptk + jppj + jphk
      Integer :: jpdd           ! No. of dry deposited species
      Integer :: jpdw           ! No. of wet deposited species
      Integer :: jpeq           ! No. of species in aqueous phase

! UKCA_MODE control features:
      LOGICAL :: L_ukca_primsu   ! T for primary sulphate aerosol emissions
      LOGICAL :: L_ukca_primss   ! T for primary sea-salt aerosol emissions
      LOGICAL :: L_ukca_primbcoc ! T for primary BC and OC aerosol emissions
      LOGICAL :: L_ukca_sedi     ! T for aerosol sedimentation
      LOGICAL :: L_ukca_nucl     ! T for aerosol sedimentation

      INTEGER :: i_mode_setup    ! Defines MODE aerosol scheme
      INTEGER :: i_mode_sizeprim ! Defines MODE size parameters
      INTEGER :: i_mode_nucscav  ! Defines nucleation scavenging method
      INTEGER :: i_mode_ddepaer  ! Defines dry deposition method
      INTEGER :: i_mode_nzts     ! No. of substeps for nucleation/sedimentation

      REAL :: mode_actdryr       ! Activation dry radius for MODE

      COMMON  /RUN_UKCA/ L_ukca_chem, L_ukca_family, L_ukca_advh2o,     &
     &         L_ukca_phot2d, L_ukca_fastj, L_ukca_trop,                &
     &         L_ukca_tropisop, L_ukca_strat, L_ukca_strattrop,         &
     &         L_ukca_stratcfc, L_ukca_wachem, L_ukca_aerchem,          &
     &         L_ukca_user, L_ukca_mode, L_ukca_dust,                   &
     &         L_ukca_rnpb, L_ukca_o3budget, L_ukca_qch4inter,          &
     &         L_ukca_isopinter, L_ukca_terpinter, L_ukca_budget2,      &
     &         L_ukca_qf11f12mbr, L_ukca_useumuivals, L_ukca_nat_sedi,  &
     &         L_ukca_het_psc, L_ukca_h2o_feedback, L_ukca_clbrcons,    &
     &         L_ukca_usero3, L_ukca_useco3, L_ukca_userelaxo3,         &
     &         L_ukca_rado3, L_ukca_radch4, L_ukca_radn2o,              &
     &         L_ukca_radf11, L_ukca_radf12, L_ukca_radf113,            &
     &         L_ukca_radf22, L_ukca_radch2o, L_ukca_useoxid,           &
     &         L_ukca_intdd, L_ukca_primsu, L_ukca_primss,              &
     &         L_ukca_primbcoc, L_ukca_sedi, L_ukca_nucl,               &
     &         jpctr, jpspec, jpbk, jptk, jppj, jphk, jpnr,             &
     &         jpdd, jpdw, jpeq,                                        &
     &         i_mode_setup, i_mode_sizeprim, i_mode_nucscav,           &
     &         i_mode_ddepaer, i_mode_nzts, mode_actdryr

      NAMELIST/RUN_UKCA/ L_ukca_chem, L_ukca_family, L_ukca_advh2o,     &
     &         L_ukca_phot2d, L_ukca_fastj, L_ukca_trop,                &
     &         L_ukca_tropisop, L_ukca_strat, L_ukca_strattrop,         &
     &         L_ukca_stratcfc, L_ukca_wachem, L_ukca_aerchem,          &
     &         L_ukca_user, L_ukca_mode, L_ukca_dust,                   &
     &         L_ukca_rnpb, L_ukca_o3budget, L_ukca_qch4inter,          &
     &         L_ukca_isopinter, L_ukca_terpinter, L_ukca_budget2,      &
     &         L_ukca_qf11f12mbr, L_ukca_useumuivals, L_ukca_nat_sedi,  &
     &         L_ukca_het_psc, L_ukca_h2o_feedback, L_ukca_clbrcons,    &
     &         L_ukca_usero3, L_ukca_useco3, L_ukca_userelaxo3,         &
     &         L_ukca_rado3, L_ukca_radch4, L_ukca_radn2o,              &
     &         L_ukca_radf11, L_ukca_radf12, L_ukca_radf113,            &
     &         L_ukca_radf22, L_ukca_radch2o, L_ukca_useoxid,           &
     &         L_ukca_intdd, L_ukca_primsu, L_ukca_primss,              &
     &         L_ukca_primbcoc, L_ukca_sedi, L_ukca_nucl,               &
     &         jpctr, jpspec, jpbk, jptk, jppj, jphk, jpnr,             &
     &         jpdd, jpdw, jpeq,                                        &
     &         i_mode_setup, i_mode_sizeprim, i_mode_nucscav,           &
     &         i_mode_ddepaer, i_mode_nzts, mode_actdryr

      INTEGER :: Instability_diagnostics  ! >0 if wanted, 0 otherwise
      INTEGER :: print_step    ! To control diagnostic printing interval
      INTEGER :: diag_interval ! diagnostic printing sampling frequency
      INTEGER :: norm_lev_start ! start level for norm diagnostics
      INTEGER :: norm_lev_end   ! end level for norm diagnostics
      INTEGER :: first_norm_print ! first timestep for norm printing
      INTEGER :: rpemax ! array size needed for diagnostic printing
      INTEGER :: rpemin ! array size needed for diagnostic printing
      INTEGER :: rpesum ! array size needed for diagnostic printing
      INTEGER :: ipesum ! array size needed for diagnostic printing
      INTEGER :: time_theta1_min      ! Timestep of min level 1 theta
      INTEGER :: time_w_max(max_model_levels) ! Timestep of max w
      INTEGER :: time_div_max(max_model_levels) ! Timestep of max div
      INTEGER :: time_div_min(max_model_levels) ! Timestep of min div
      INTEGER :: time_lapse_min(max_model_levels) ! Timestep of min
      INTEGER :: time_max_shear(max_model_levels) !Timestep max shear
      INTEGER :: time_max_wind(max_model_levels) ! Timestep of max wind
      INTEGER :: time_KE_max(max_model_levels) ! Timestep of max KE
      INTEGER :: time_KE_min(max_model_levels) ! Timestep of min KE
      INTEGER :: time_noise_max(max_model_levels) ! Timestep of max

      REAL:: frictional_timescale(max_model_levels) ! For idealised case
      REAL :: w_print_limit ! w Threshold for diagnostic printing
      REAL :: w_conv_limit  ! w Threshold for limiting convection
      REAL :: tropics_deg  ! define latitude for tropics
      REAL :: min_theta1_run                  ! Min theta level 1
      REAL :: dtheta1_run   ! Largest -ve delta theta at min theta1
      REAL :: max_w_run(max_model_levels)  ! Max w at a level
      REAL :: max_div_run(max_model_levels) ! Max divergence at a level
      REAL :: min_div_run(max_model_levels) ! Min divergence at a level
      REAL :: min_lapse_run(max_model_levels) ! Min dtheta/dz at a level
      REAL :: max_shear_run(max_model_levels) ! Max shear at a level
      REAL :: max_wind_run(max_model_levels) ! Max wind at a level
      REAL :: max_KE_run(max_model_levels)   ! Max KE at a level
      REAL :: min_KE_run(max_model_levels)   ! Min KE at a level
      REAL :: max_noise_run(max_model_levels) ! Max noise at a level

!     Problem_number not set here  ! Now controlled by namelist input
!     Instability_diagnostics      ! Now controlled by namelist input
!     frictional_timescale         ! Now intitialised in SETCONA

      LOGICAL :: L_idealised_data
      LOGICAL :: L_simple_friction
      LOGICAL :: L_print_w
      LOGICAL :: L_print_div
      LOGICAL :: L_diag_print_ops     ! diagnostic prints for ops
      LOGICAL :: L_print_pe     ! print diagnostics on all pe's if true
      LOGICAL :: L_print_shear        ! wind diagnostic prints
      LOGICAL :: L_print_max_wind     ! wind diagnostic prints
      LOGICAL :: L_diag_L2norms  ! l2norm diagnostic prints
      LOGICAL :: L_diag_L2helm   ! l2norm diagnostic prints from solver
      LOGICAL :: L_diag_wind     ! wind diagnostic prints
      LOGICAL :: L_diag_noise    ! w diagnostic prints
      LOGICAL :: L_flush6        ! if T then flush buffers on failure
      LOGICAL :: L_inviscid  ! Use prognostic inviscid theta at surface
      LOGICAL :: L_deep        ! Deep/shallow switch
      LOGICAL :: L_hydrostatic ! Hydrostatic option
      LOGICAL :: L_geostrophic ! Geostrophic balance at start
      LOGICAL :: L_solver      ! Solver active if true
      LOGICAL :: L_trap_errors ! Error trapping
      LOGICAL :: L_uv_zero     ! Set winds to zero
      LOGICAL :: L_diag_print  ! Print diagnostics
      LOGICAL :: L_print_lapse ! Print lapse_rate diagnostics
      LOGICAL :: L_print_wmax ! Print max w diagnostic
      LOGICAL :: L_print_theta1 ! Print level 1 theta diagnostic

!     L_idealised_data             ! Now controlled by namelist input

!------------------   Diagnostics:   --------------------------------

      COMMON  /RUN_Diagnostics/                                         &
     &  L_diag_print, L_diag_print_ops, L_flush6, print_step,           &
     &  diag_interval,w_print_limit, rpemax, rpemin, ipesum, rpesum,    &
     &  L_print_w, L_print_wmax,L_print_lapse, L_print_theta1,          &
     &  L_print_div, L_diag_wind, L_print_shear, L_print_max_wind,      &
     &  L_print_pe, L_diag_noise, L_diag_L2norms, L_diag_L2helm,        &
     &  norm_lev_start, norm_lev_end, first_norm_print,                 &
     &  max_w_run, min_theta1_run, dtheta1_run,                         &
     &  max_div_run, min_div_run, min_lapse_run,                        &
     &  max_shear_run, max_wind_run,                                    &
     &  max_noise_run, max_KE_run, min_KE_run,                          &
     &  time_KE_max, time_KE_min,                                       &
     &  time_w_max, time_div_max, time_div_min, time_lapse_min,         &
     &  time_max_shear, time_max_wind,                                  &
     &  time_theta1_min, time_noise_max

!------------------   Dynamics:   --------------------------------------

      ! Generalised dynamics:

      LOGICAL :: L_Physics   ! T: physics to be included
      LOGICAL :: L_Run_With_Physics2   ! T: physics2 to be included
      LOGICAL :: L_perturb_IC_theta    ! T: perturb theta on ts1   
      LOGICAL :: L_Backwards ! F: Integrate backwards without physics.
      LOGICAL :: L_Primitive ! T: Primitive Equation model
                   ! F: Semi-Geostrophic model (no dynamics)
      LOGICAL :: L_dry         ! T: run with no moisture
      LOGICAL :: L_adjust_wet  ! T: perform simple dry&moist adjustment
      LOGICAL :: L_mix_ratio

      LOGICAL :: L_LBC_balance  ! T: impose vertically balanced Exners
                                !    and rho in LBCs (set w=0)
                                ! F: leave Exners alone

      ! Switch for new (correct) lbcs developed in 2008
      LOGICAL :: L_lbc_new       ! .true. for new lbc switch

      LOGICAL :: L_free_slip    ! .true. for free-slip lower boundary

      ! T: reset polar values to mean value every polar_reset_timesteps
      LOGICAL, PARAMETER :: L_polar_reset = .false.

      ! interval for resetting polar to mean
      INTEGER, PARAMETER :: polar_reset_timesteps = 1

      INTEGER :: IntRand_seed

      ! Time weights for integration scheme:
      REAL :: alpha_1
      REAL :: alpha_2
      REAL :: alpha_3
      REAL :: alpha_4

!  variable resolution control
      LOGICAL :: L_regular  ! true if NOT variable resolution
      LOGICAL :: L_qwaterload  ! true if using adding waterloading terms

      INTEGER :: lam_var    ! number of variable res. lambda intervals
      INTEGER :: phi_var    ! number of variable res. phi intervals

      REAL :: var_ratio  ! grid-stretcing ratio for var grid
      REAL :: lam_ratio  ! scaling of original grid to high res grid
      REAL :: phi_ratio  ! scaling of original grid to high res grid
      REAL :: lam_frac  ! proportion of reg. points in West of domain
      REAL :: phi_frac  ! proportion of reg. points in South of domain
      REAL :: lambda_p_end
      REAL :: phi_p_end
      REAL :: lambda_u_end
      REAL :: phi_v_end
      REAL :: dlambda_p_end
      REAL :: dphi_p_end
      REAL :: dlambda_u_end
      REAL :: dphi_v_end

      ! GCR scheme= Generalized conjugate residual elliptic equation
      ! solver:
      LOGICAL :: GCR_use_tol_abs
      LOGICAL :: GCR_zero_init_guess
      LOGICAL :: GCR_use_residual_Tol

      ! T then use full equation on RHS on second and subsequent ADI
      ! timesteps
      LOGICAL :: GCR_adi_add_full_soln

      LOGICAL :: L_gcr_fast_x ! use faster non reproducible code
! Speed up solver convergence when dynamics-physics cycling
! is used.
       Logical :: L_GCR_cycle_opt

      LOGICAL :: L_interp_depart ! interpolate to find u,v departure pts

! ----------------------- include file: GCR_ITER_DIM -------------------
! Description: Introduce integer parameters to:
!              (1) declare the maximum allowed number of dynamics 
!                  cycles (timestep iterations)
!              (2) declare the number of GCR Iterations analysis steps
!
! Code Author: M. Diamantakis
! Current Code Owner: A. Malcom
!
      INTEGER, PARAMETER :: MAX_NUMCYCLES = 10 ! Max number of 
                                               ! dynamics cycles
      INTEGER, PARAMETER :: GCR_ANAL_STEPS = 3 ! Number of GCR Iterations
                                               ! analysis step

      INTEGER :: GCR_max_iterations
      INTEGER :: GCR_diagnostics    !RR replace by STASHflag?
      INTEGER :: GCR_its_avg_step(GCR_anal_steps)   
                                             ! GCR Iterations analysis 
                                             ! steps
      INTEGER, DIMENSION(max_numcycles) :: GCR_its_switch
                                             ! GCR Iterations analysis 
                                             ! switch
      INTEGER, DIMENSION(max_numcycles) :: GCR_max_its  
                                             ! GCR Max iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_min_its  
                                             ! GCR Min iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_max_time 
                                             ! Timestep number
                                             ! GCR Max iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_min_time 
                                             ! Timestep number
                                             ! GCR Min iterations
      INTEGER, DIMENSION(max_numcycles) :: GCR_sum_its  
                                             ! Sum iterations
                                             ! over test period
      INTEGER :: GCR_precon_option
      INTEGER :: GCR_n_ADI_pseudo_timesteps

      ! After how many iterations do we restart
      INTEGER :: GCR_Restart_value

      REAL :: GCR_tol_res
      REAL :: GCR_tol_abs
      REAL :: GCR_ADI_pseudo_timestep
      REAL :: G_term_tol
!      G_term_tol             ! Now controlled by namelist input
!
!  Parameters for dynamics-physics cycling
       INTEGER :: NumCycles
       LOGICAL :: L_new_tdisc
       REAL    :: extrp_weight
       REAL    :: GCR_tol_abs2
       REAL    :: GCR_tol_res2
       REAL    :: alpha_1_2, alpha_2_2, alpha_3_2, alpha_4_2
! Parameter for fully-interpolating theta option
       LOGICAL :: L_fint_theta
! Number of physics2 substeps (default is 1)
       INTEGER :: Num_Substeps
! Logical switch to enable phys2 substepping
       Logical :: L_phys2_substep
      COMMON  /RUN_Dyn/                                                 &
     &  L_idealised_data,                                               &
     &  L_Physics, L_phys2_substep, L_Backwards, L_Primitive,           &
     &  L_GCR_cycle_opt, L_Run_With_Physics2, L_perturb_IC_theta,       &
     &  L_mix_ratio, L_new_tdisc, L_fint_theta, IntRand_seed,           &
     &  L_free_slip,                                                    &
     &  L_simple_friction, L_inviscid, L_deep, L_hydrostatic,           &
     &  L_dry, L_adjust_wet, L_LBC_balance, L_lbc_new,                  &
     &  L_geostrophic, L_solver, L_trap_errors, L_uv_zero,              &
     &  L_regular, L_qwaterload, n_rims_to_do, lam_var, phi_var,        &
     &  var_ratio, lam_ratio, phi_ratio, lam_frac, phi_frac,            &
     &  lambda_p_end,  phi_p_end, lambda_u_end,  phi_v_end,             &
     &  dlambda_p_end,  dphi_p_end, dlambda_u_end,  dphi_v_end,         &
     &  alpha_1,alpha_2,alpha_3,alpha_4,                                &
     &  alpha_1_2, alpha_2_2, alpha_3_2, alpha_4_2,                     &
     &  Num_Substeps, NumCycles,                                        &
     &  extrp_weight,                                                   &
     &  GCR_diagnostics, GCR_max_its, GCR_min_its,                      &
     &  GCR_max_time, GCR_min_time,                                     &
     &  GCR_sum_its, GCR_its_switch, GCR_its_avg_step,                  &
     &  GCR_use_tol_abs,GCR_zero_init_guess,GCR_use_residual_Tol,       &
     &  GCR_adi_add_full_soln,L_gcr_fast_x,                             &
     &  L_interp_depart,                                                &
     &  GCR_max_iterations, GCR_precon_option,                          &
     &  GCR_n_ADI_pseudo_timesteps,GCR_Restart_value,GCR_tol_res,       &
     &  GCR_tol_res2,GCR_tol_abs,GCR_tol_abs2,GCR_ADI_pseudo_timestep,  &
     &  G_term_tol
      NAMELIST/RUN_Dyn/                                                 &
     &  L_Physics,L_phys2_substep,L_Backwards,L_Primitive,              &
     &  L_GCR_cycle_opt, L_Run_With_Physics2, L_perturb_IC_theta,       &
     &  L_mix_ratio, L_new_tdisc, L_fint_theta, IntRand_seed,           &
     &  L_free_slip, L_dry, L_adjust_wet,                               &
     &  L_LBC_balance, L_lbc_new,                                       &
     &  L_regular, L_qwaterload, n_rims_to_do, lam_var, phi_var,        &
     &  var_ratio, lam_ratio, phi_ratio, lam_frac, phi_frac,            &
     &  lambda_p_end,  phi_p_end, lambda_u_end,  phi_v_end,             &
     &  dlambda_p_end,  dphi_p_end, dlambda_u_end,  dphi_v_end,         &
     &  alpha_1,alpha_2,alpha_3,alpha_4,                                &
     &  alpha_1_2, alpha_2_2, alpha_3_2, alpha_4_2,                     &
     &  Num_Substeps, NumCycles,                                        &
     &  extrp_weight,                                                   &
     &  GCR_use_tol_abs,GCR_zero_init_guess,GCR_use_residual_Tol,       &
     &  GCR_adi_add_full_soln,L_gcr_fast_x,                             &
     &  L_interp_depart,                                                &
     &  GCR_max_iterations,GCR_diagnostics,GCR_precon_option,           &
     &  GCR_n_ADI_pseudo_timesteps,GCR_Restart_value,GCR_tol_res,       &
     &  GCR_tol_res2,GCR_tol_abs,GCR_tol_abs2,GCR_ADI_pseudo_timestep,  &

     &  G_term_tol, GCR_its_avg_step

! Semi-Lagrangian advection control options:


      INTEGER, PARAMETER :: Number_SL_choices = 4
      INTEGER, PARAMETER :: Theta_SL = 1
      INTEGER, PARAMETER :: moist_SL = 2
      INTEGER, PARAMETER :: Wind_SL = 3
      INTEGER, PARAMETER :: rho_SL = 4

      ! T if conservation to be enforced.
      LOGICAL:: L_conserv(Number_SL_choices)

      ! T if monotonicity to be enforced.
      LOGICAL:: L_mono(Number_SL_choices)

      ! T if high order interpolation scheme to be used.
      LOGICAL:: L_high(Number_SL_choices)

      ! T if high order scheme to be used in Ritchie routine.
      LOGICAL:: L_Ritchie_high

      ! T if monotone scheme to be used in Ritchie routine.
      LOGICAL:: L_Ritchie_mono

      ! T then only perform vector co-ordinate geometry in 2d.
      LOGICAL:: L_2d_sl_geometry


      ! T: sl code bit repoducible with any sensible halo size
      LOGICAL:: L_sl_halo_reprod

      ! improved nonhydrostatic conservation  of q, qcl and qcf
      LOGICAL:: L_moist_nonhydro_conserve

      ! Run tracer advection code with conservation on          
      LOGICAL :: L_conserve_tracers

      ! a code saying which high monotone scheme to use.
      ! 1 = tensor tri-cubic lagrange (j,i,k)
      ! no other options available at present.
      INTEGER :: high_order_scheme(Number_SL_choices)

      ! a code saying which  monotone scheme to use.
      ! 1 = tri-linear order (j,i,k)
      ! no other options available at present.
      INTEGER :: monotone_scheme(Number_SL_choices)
      REAL :: thmono_height       ! top for monotone fully-interp theta
      INTEGER :: thmono_levels ! levels for monotone fully-interp theta
!                 increments used to limit non-interp increments

      ! a code saying which high  order scheme to use in Ritchie routine
      INTEGER :: Ritchie_high_order_scheme

      ! a code saying which monotone scheme to use in Ritchie routine.
      INTEGER :: Ritchie_monotone_scheme

      INTEGER :: Depart_scheme     ! which departure point scheme to use

      ! for the chosen departure point scheme how many iterations/terms
      ! to use
      INTEGER :: Depart_order

      ! used in interpolation code:
      ! number of levels either side of default level to search.
      INTEGER :: interp_vertical_search_tol

      ! Set of look-up tables for searching on variable grids
      INTEGER :: look_lam(max_look)  ! lambda p search
      INTEGER :: look_phi(max_look)  ! phi p search

      ! Set of minumum search grid lengths for variable grids
      REAL ::  recip_dlam   ! smallest delta_lambda p
      REAL ::  recip_dphi   ! smallest delta_phi p

      INTEGER :: halo_lam  ! halo for lamp look-up table
      INTEGER :: halo_phi  ! halo for phip look-up table


      COMMON  /RUN_SL/L_conserv,L_mono,L_high,L_Ritchie_high,           &
     &  L_Ritchie_mono, thmono_height, thmono_levels, L_2d_sl_geometry, &
     &  L_sl_halo_reprod,high_order_scheme,monotone_scheme,             &
     &  Ritchie_high_order_scheme,Ritchie_monotone_scheme,              &
     &  Depart_scheme,Depart_order,interp_vertical_search_tol,          &
     &  recip_dlam, recip_dphi, halo_lam, halo_phi,                     &
     &  look_lam, look_phi,                                             &
     &  L_moist_nonhydro_conserve, L_conserve_tracers,                  &
     &  Instability_diagnostics

      NAMELIST/RUN_SL/L_conserv,L_mono,L_high,L_Ritchie_high,           &
     &  L_Ritchie_mono, thmono_height, L_2d_sl_geometry,                &
     &  L_sl_halo_reprod,high_order_scheme,monotone_scheme,             &
     &  Ritchie_high_order_scheme,Ritchie_monotone_scheme,              &
     &  Depart_scheme,Depart_order,interp_vertical_search_tol,          &
     &  L_moist_nonhydro_conserve, L_conserve_tracers

      ! Diffusion control options:
      LOGICAL :: L_diffusion
      LOGICAL :: L_upper_ramp   ! ramp upper-level diffusion
      LOGICAL :: L_adjust_theta ! activate convective adjustment
      LOGICAL :: L_vdiff_uv    ! activate targeted diffusion uv
      LOGICAL :: L_filter    ! activate polar filter or diffusion
      LOGICAL :: L_filter_incs    ! activate polar filter of incs
      LOGICAL :: L_pftheta   ! activate polar filter of theta
      LOGICAL :: L_pfuv      ! activate polar filter of u,v
      LOGICAL :: L_pfw       ! activate polar filter of w
      LOGICAL :: L_pofil_new  ! activate new polar filter
      LOGICAL :: L_pfexner  ! activate polar filter of Exner pressure
      LOGICAL :: L_pofil_hadgem2 ! run with HadGEM2 polar filter setting   
      LOGICAL :: L_diff_exner  ! activate diffusion of Exner pressure
      LOGICAL :: L_pfcomb   ! combined polar filter/diffusion active
      LOGICAL :: L_sponge   ! activate lateral boundariessponge zones
      LOGICAL :: L_pfincs    ! activate polar filter of incs
      LOGICAL :: L_diff_thermo ! horiz. diffusion of theta
      LOGICAL :: L_diff_wind   ! horiz. diffusion of u, v
      LOGICAL :: L_diff_w      ! horiz. diffusion of w
      LOGICAL :: L_diff_incs   ! horiz. diffusion of increments
      LOGICAL :: L_diff_auto   ! UM calculates diffusion parameters
      LOGICAL :: L_tardiff_q     ! activate targeted diffusion q
      LOGICAL :: L_diff_ctl      ! general diffusion control
! L_diff_ctl == L_diffusion .or. L_cdiffusion .or. L_vertical_diffusion
!               .or. L_divdamp .or.  L_tardiff_q .or. L_diag_print
      LOGICAL :: L_cdiffusion
      LOGICAL :: L_subfilter_horiz  ! subgrid turbulence in horizontal
      LOGICAL :: L_subfilter_vert   ! subgrid turbulence in vertical  
      LOGICAL :: L_subfilter_blend  ! blend diffusion coeffs

      ! value * del^2: order=2
      INTEGER  :: diffusion_order_thermo(max_model_levels)

      ! gives del^4 diffusion
      INTEGER  :: diffusion_order_wind(max_model_levels)
      INTEGER  :: diffusion_order_w(max_model_levels)
      INTEGER  :: diffusion_order_q(max_model_levels)
      INTEGER :: u_begin(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: u_end(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_begin(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_end(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: u_sweeps(max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_sweeps(max_121_rows) ! Sweep control on 121 filter
      INTEGER :: max_sweeps ! Max sweeps wanted for 121 filter
      INTEGER :: global_u_filter ! Sweep control on 121 filter
      INTEGER :: global_v_filter ! Sweep control on 121 filter
      INTEGER :: diff_order_thermo   ! diffusion order for theta
      INTEGER :: diff_order_wind ! diffusion order for winds
      INTEGER :: diff_timescale_thermo  ! diffusion timescale for theta
      INTEGER :: diff_timescale_wind    ! diffusion timescale for wind
      INTEGER :: vdiffuv_timescale ! diffusion e-folding timesteps
      INTEGER :: vdiffuv_start ! start level  targeted diffusion
      INTEGER :: vdiffuv_end ! end level targeted diffusion
      INTEGER :: adjust_theta_start ! start level convective adjustment
      INTEGER :: adjust_theta_end   ! end levelconvective adjustment
      INTEGER :: top_filt_start ! start level upper-level diffusion
      INTEGER :: top_filt_end   ! end level upper-level diffusion
      INTEGER :: sponge_ew   ! left/right boundaries sponge zone width
      INTEGER :: sponge_ns   ! north/south boundaries sponge zone width
      INTEGER :: sponge_power ! sponge zone weighting order
      INTEGER :: tardiffq_test ! test level test w targetted diffusion
      INTEGER :: tardiffq_start ! start level test w targetted diffusion
      INTEGER :: tardiffq_end ! end level test w targetted diffusion
      INTEGER :: turb_startlev_horiz   ! 1st lev for horiz subgrid turb
      INTEGER :: turb_endlev_horiz    ! last lev for horiz subgrid turb
      INTEGER :: turb_startlev_vert   ! 1st lev for vert subgrid turb
      INTEGER :: turb_endlev_vert     ! last lev for vert subgrid turb

      !  level - assume surfaces are horizontal
      INTEGER  :: horizontal_level
      INTEGER  ::  tar_horizontal ! steep slope test targeted diffusion

      REAL :: diffusion_coefficient_thermo(max_model_levels)
      REAL :: diffusion_coefficient_wind(max_model_levels)
      REAL :: diffusion_coefficient_w(max_model_levels)
      REAL :: diffusion_coefficient_q(max_model_levels)
      REAL :: vdiffuv_test ! test to activate shear diffusion of u,v
      REAL :: vdiffuv_factor ! vertical diffusion coeff for shear diff
      REAL :: scale_ratio ! Pass control on 121 filter
      REAL :: ref_lat_deg ! Reference lat for auto diffusion
      REAL :: top_diff    !  upper-level diffusion coefficient
      REAL :: up_diff_scale    !  upper-level diffusion ramping factor

      REAL :: up_diff(max_updiff_levels) ! upper-level diffusion coeff
      REAL :: sponge_wts_ew(max_sponge_width) ! sponge weights
      REAL :: sponge_wts_ns(max_sponge_width) ! sponge weights

      REAL :: diff_coeff_ref ! EW diffusion coefficient at polar cap
      REAL :: diff_coeff_thermo  ! NS theta diffusion coeff
      REAL :: diff_coeff_wind    ! NS u,v diffusion coeff
      REAL :: diff_coeff_phi    ! North-South diffusion coefficient
!   reference latitudes for filtering and diffusion
      REAL :: polar_cap  ! Apply 1-2-1 filter polewards
      REAL :: tardiffq_factor ! targeted diffusion coefficient
      REAL :: diff_factor
      REAL :: mix_factor

      ! Divergence damping control options:
      LOGICAL :: L_divdamp

      REAL :: div_damp_coefficient(max_model_levels)

      ! Polar filter control options:
      LOGICAL :: L_polar_filter
      ! T: use polar filter to filter increment
      LOGICAL :: L_polar_filter_incs

      REAL :: polar_filter_north_lat_limit
      REAL :: polar_filter_south_lat_limit
      REAL :: polar_filter_coefficient

      ! amount in radians to increment start latitude by per sweep
      REAL :: polar_filter_step_per_sweep

      ! max latitude at which filter can star
      REAL :: polar_filter_lat_limit

      ! number of sweeps of filter to do
      INTEGER :: polar_filter_n_sweeps

      ! Vertical Diffusion control options:
      LOGICAL :: L_vertical_diffusion
      LOGICAL :: L_ramp

      INTEGER :: level_start_wind
      INTEGER :: level_stop_wind
      INTEGER :: level_start_theta
      INTEGER :: level_stop_theta
      INTEGER :: level_start_q
      INTEGER :: level_stop_q

      REAL :: vert_diffusion_coeff_wind
      REAL :: vert_diffusion_coeff_theta
      REAL :: vert_diffusion_coeff_q
      REAL :: ramp_lat_radians

      ! Moisture resetting control options:
      LOGICAL :: l_qpos         ! logical to run qpos code

      ! true to do Method 2 false to do Method 1
      LOGICAL :: l_q_pos_local


      REAL :: qlimit         !    lowest allowed value of q

      COMMON  /RUN_Diffusion/L_diffusion,diffusion_order_thermo,        &
     &  diffusion_order_wind,diffusion_order_w,diffusion_order_q,       &
     &  diffusion_coefficient_thermo,diffusion_coefficient_wind,        &
     &  diffusion_coefficient_w,diffusion_coefficient_q,                &
     &  diff_order_thermo, diff_timescale_thermo, diff_coeff_thermo,    &
     &  diff_order_wind, diff_timescale_wind, diff_coeff_wind,          &
     &  horizontal_level, tar_horizontal, tropics_deg,                  &
     &  L_cdiffusion,L_ramp, ramp_lat_radians,                          &
     &  L_divdamp,div_damp_coefficient,                                 &
     &  L_polar_filter,L_polar_filter_incs,                             &
     &  polar_filter_north_lat_limit,polar_filter_south_lat_limit,      &
     &  polar_filter_coefficient,polar_filter_step_per_sweep,           &
     &  polar_filter_lat_limit,polar_filter_n_sweeps,                   &
     &  l_qpos,  l_q_pos_local, qlimit,                                 &
     &  L_diff_ctl, L_tardiff_q, w_conv_limit,                          &
     &  tardiffq_factor, tardiffq_test, tardiffq_start, tardiffq_end,   &
     &  L_subfilter_horiz, L_subfilter_vert, L_subfilter_blend,         &
     &  diff_factor, mix_factor, turb_startlev_horiz, turb_endlev_horiz,&
     &  turb_startlev_vert, turb_endlev_vert,                           &
     &  L_vertical_diffusion, L_pofil_hadgem2,                          &
     &  L_pofil_new, L_pfcomb, L_pftheta, L_pfuv, L_pfw, L_pfincs,      &
     &  L_pfexner, L_filter, L_filter_incs, L_diff_exner, L_diff_auto,  &
     &  L_diff_incs, L_diff_thermo, L_diff_wind, L_diff_w,              &
     &  L_vdiff_uv, L_adjust_theta,                                     &
     &  level_start_wind, level_stop_wind,                              &
     &  level_start_theta, level_stop_theta,                            &
     &  level_start_q, level_stop_q,                                    &
     &  polar_cap, scale_ratio, ref_lat_deg, max_sweeps,                &
     &  L_upper_ramp, top_filt_start, top_filt_end, top_diff,           &
     &  up_diff, up_diff_scale, sponge_wts_ew, sponge_wts_ns,           &
     &  u_begin, u_end, v_begin, v_end, u_sweeps, v_sweeps,             &
     &  global_u_filter, global_v_filter,                               &
     &  diff_coeff_ref, diff_coeff_phi, vdiffuv_timescale,              &
     &  vdiffuv_test, vdiffuv_factor, vdiffuv_start, vdiffuv_end,       &
     &  adjust_theta_start, adjust_theta_end,                           &
     &  L_sponge, sponge_power, sponge_ew, sponge_ns,                   &
     &  vert_diffusion_coeff_wind, vert_diffusion_coeff_theta,          &
     &  vert_diffusion_coeff_q

      NAMELIST/RUN_Diffusion/L_diffusion,diffusion_order_thermo,        &
     &  diffusion_order_wind,diffusion_order_w,diffusion_order_q,       &
     &  diffusion_coefficient_thermo,diffusion_coefficient_wind,        &
     &  diffusion_coefficient_w,diffusion_coefficient_q,                &
     &  diff_order_thermo, diff_timescale_thermo,                       &
     &  diff_order_wind, diff_timescale_wind,                           &
     &  horizontal_level, tar_horizontal,                               &
     &  L_cdiffusion,L_ramp, ramp_lat_radians,                          &
     &  L_divdamp,div_damp_coefficient,                                 &
     &  L_polar_filter,L_polar_filter_incs,                             &
     &  polar_filter_north_lat_limit,polar_filter_south_lat_limit,      &
     &  polar_filter_coefficient,polar_filter_step_per_sweep,           &
     &  polar_filter_lat_limit,polar_filter_n_sweeps,                   &
     &  l_qpos,  l_q_pos_local, qlimit,                                 &
     &  L_diff_ctl, L_tardiff_q, w_conv_limit,                          &
     &  tardiffq_factor, tardiffq_test, tardiffq_start, tardiffq_end,   &
     &  L_subfilter_horiz, L_subfilter_vert, L_subfilter_blend,         &
     &  diff_factor, mix_factor, turb_startlev_horiz, turb_endlev_horiz,&
     &  turb_startlev_vert, turb_endlev_vert,                           &
     &  L_vertical_diffusion,                                           &
     &  L_pftheta, L_pfuv, L_pfw, L_pfincs, L_pfexner, L_pofil_hadgem2, &
     &  L_diff_incs, L_diff_thermo, L_diff_wind, L_diff_w,              &
     &  level_start_wind, level_stop_wind,                              &
     &  level_start_theta, level_stop_theta,                            &
     &  level_start_q, level_stop_q,                                    &
     &  L_pofil_new, L_diff_auto,                                       &
     &  diff_coeff_ref, polar_cap, scale_ratio, ref_lat_deg,            &
     &  max_sweeps, L_upper_ramp, up_diff_scale, top_diff,              &
     &  top_filt_start, top_filt_end, L_vdiff_uv, vdiffuv_timescale,    &
     &  vdiffuv_test, vdiffuv_factor, vdiffuv_start, vdiffuv_end,       &
     &  L_adjust_theta, adjust_theta_start, adjust_theta_end,           &
     &  L_sponge, sponge_power, sponge_ew, sponge_ns,                   &
     &  vert_diffusion_coeff_wind, vert_diffusion_coeff_theta,          &
     &  vert_diffusion_coeff_q,                                         &
     &  L_diag_print, L_diag_print_ops, L_print_pe,                     &
     &  L_print_w, L_print_wmax, L_print_lapse, L_print_theta1,         &
     &  L_print_div, L_diag_wind, L_print_shear, L_print_max_wind,      &
     &  L_diag_noise, L_diag_L2norms, L_diag_L2helm,                    &
     &  norm_lev_start, norm_lev_end, first_norm_print,                 &
     &  print_step, diag_interval, w_print_limit, L_Flush6

!------------------   Dynamical core   -------------------------------
! Suarez-Held variables
      REAL :: SuHe_newtonian_timescale_ka
      REAL :: SuHe_newtonian_timescale_ks
      REAL :: SuHe_pole_equ_deltaT
      REAL :: SuHe_static_stab
      REAL :: base_frictional_timescale
      REAL :: SuHe_sigma_cutoff
      REAL :: SuHe_level_weight(max_model_levels)
      REAL :: friction_level(max_model_levels)

      INTEGER :: SuHe_relax
      INTEGER :: SuHe_fric

      LOGICAL :: L_SH_Williamson

      COMMON/Run_Dyncore/                                               &
     &  SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks,       &
     &  SuHe_pole_equ_deltaT, SuHe_static_stab,                         &
     &  base_frictional_timescale, SuHe_sigma_cutoff,                   &
     &  L_SH_Williamson, SuHe_relax, SuHe_fric,                         &
     &  SuHe_level_weight, frictional_timescale, friction_level

!------------------  Idealised model   ----------------------------

      INTEGER,PARAMETER:: max_num_profile_data = 100
      INTEGER,PARAMETER:: max_num_force_times = 100
      INTEGER,PARAMETER:: idl_max_num_bubbles = 3

! Idealised  variables
      REAL :: h_o
      REAL :: h_o_actual  ! height of growing mountain
      REAL :: h_o_per_step ! height change per step of growing mountain
      REAL :: lambda_fraction
      REAL :: phi_fraction
      REAL :: half_width_x
      REAL :: half_width_y
      REAL :: Witch_power
      REAL :: plat_size_x
      REAL :: plat_size_y
      REAL :: height_domain
      REAL :: delta_x
      REAL :: delta_y
      REAL :: big_factor
      REAL :: mag
      REAL :: first_theta_height
      REAL :: thin_theta_height
      REAL :: p_surface
      REAL :: theta_surface
      REAL :: dtheta_dz1(3)
      REAL :: height_dz1(3)
      REAL :: Brunt_Vaisala
      REAL :: u_in(4)
      REAL :: v_in(4)
      REAL :: height_u_in(3)
      REAL :: u_ramp_start
      REAL :: u_ramp_end
      REAL :: ujet_lat
      REAL :: ujet_width
      REAL :: t_horizfn_data(10)
      REAL :: q1
      REAL :: theta_ref(max_model_levels)
      REAL :: rho_ref(max_model_levels)
      REAL :: exner_ref(max_model_levels + 1)
      REAL :: q_ref(max_model_levels)
      REAL :: u_ref(max_model_levels)
      REAL :: v_ref(max_model_levels)
      REAL :: z_orog_print(0:max_model_levels)
      REAL :: f_plane
      REAL :: ff_plane
      REAL :: r_plane
      REAL :: zprofile_data(max_num_profile_data)
      REAL :: tprofile_data(max_num_profile_data)
      REAL :: qprofile_data(max_num_profile_data)
      REAL :: z_uvprofile_data(max_num_profile_data)
      REAL :: uprofile_data(max_num_profile_data)
      REAL :: vprofile_data(max_num_profile_data)
      REAL :: tforce_time_interval
      REAL :: qforce_time_interval
      REAL :: uvforce_time_interval
      REAL :: newtonian_timescale
      REAL :: z_tforce_data(max_num_profile_data)
      REAL :: z_qforce_data(max_num_profile_data)
      REAL :: z_uvforce_data(max_num_profile_data)
      REAL :: tforce_data(max_num_profile_data, max_num_force_times)
      REAL :: qforce_data(max_num_profile_data, max_num_force_times)
      REAL :: uforce_data(max_num_profile_data, max_num_force_times)
      REAL :: vforce_data(max_num_profile_data, max_num_force_times)
      REAL :: tforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: qforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: uforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: vforce_data_modlev(max_model_levels, max_num_force_times)
      REAL :: pforce_time_interval
      REAL :: p_surface_data(max_num_force_times)
      REAL :: perturb_factor
      REAL :: perturb_magnitude_t
      REAL :: perturb_magnitude_q
      REAL :: perturb_height(2)
      REAL :: orog_hgt_lbc
      REAL :: zprofile_orog
      REAL :: hf
      REAL :: cool_rate
      REAL :: IdlSurfFluxSeaParams(10) ! Idealised surface flux params
      REAL :: roughlen_z0m   
      REAL :: roughlen_z0h
      ! Idealised bubbles
      REAL :: idl_bubble_max(idl_max_num_bubbles) ! Bubble max amplitude
      REAL :: idl_bubble_height(idl_max_num_bubbles)  ! Bubble height
      REAL :: idl_bubble_width(idl_max_num_bubbles)   ! Bubble width
      REAL :: idl_bubble_depth(idl_max_num_bubbles)   ! Bubble depth
      ! Bubble x-offset, y-offset in normalised units (0:1)
      ! (0.5=domain centre)
      REAL :: idl_bubble_xoffset(idl_max_num_bubbles)
      REAL :: idl_bubble_yoffset(idl_max_num_bubbles)
      REAL :: DMPTIM, HDMP, ZDMP   ! Damping layer values
      REAL :: u_geo, v_geo         ! Geostrophic wind


      INTEGER :: n_rims_to_do   ! rim size for LAM
      INTEGER :: surface_type
      INTEGER :: grow_steps
      INTEGER :: grid_number
      INTEGER :: grid_flat
      INTEGER :: tprofile_number
      INTEGER :: qprofile_number
      INTEGER :: uvprofile_number
      INTEGER :: num_profile_data
      INTEGER :: num_uvprofile_data
      INTEGER :: t_horizfn_number
      INTEGER :: uv_horizfn_number
      INTEGER :: pforce_option
      INTEGER :: num_pforce_times
      INTEGER :: tforce_option
      INTEGER :: qforce_option
      INTEGER :: uvforce_option
      INTEGER :: num_tforce_levels
      INTEGER :: num_tforce_times
      INTEGER :: num_qforce_levels
      INTEGER :: num_qforce_times
      INTEGER :: num_uvforce_levels
      INTEGER :: num_uvforce_times
      INTEGER :: IdlSurfFluxSeaOption  ! Idealised surface flux option
      INTEGER :: first_constant_r_rho_level_new
      INTEGER :: big_layers
      INTEGER :: transit_layers
      INTEGER :: mod_layers
      INTEGER :: idl_bubble_option(idl_max_num_bubbles) ! Bubble option
      INTEGER :: idl_interp_option  ! Profile interpolation option
      INTEGER :: perturb_type

      LOGICAL :: L_initialise_data
      LOGICAL :: L_constant_dz
      LOGICAL :: L_trivial_trigs !.false. for Cartesian coords (lat=0.0)
      LOGICAL :: L_idl_bubble_saturate(idl_max_num_bubbles)
      LOGICAL :: L_fixed_lbcs
      LOGICAL :: L_fix_orog_hgt_lbc
      LOGICAL :: L_pressure_balance
      LOGICAL :: L_wind_balance
      LOGICAL :: L_rotate_winds
      LOGICAL :: L_polar_wind_zero
      LOGICAL :: L_vert_Coriolis
      LOGICAL :: L_rotating     ! .true. for Earth's rotation
      LOGICAL :: L_perturb      ! add random perturb. to surface theta
      LOGICAL :: L_code_test    ! User switch for testing code
      LOGICAL :: L_pforce
      LOGICAL :: L_baroclinic
      LOGICAL :: L_cyclone
      LOGICAL :: L_force
      LOGICAL :: L_force_lbc
      LOGICAL :: L_perturb_t
      LOGICAL :: L_perturb_q
      LOGICAL :: L_perturb_correlate_tq
      LOGICAL :: L_perturb_correlate_vert
      LOGICAL :: L_perturb_correlate_time
      LOGICAL :: L_damp      ! Logical for damping layer
      LOGICAL :: L_geo_for ! Logical for geostrophic wind forcing
      LOGICAL :: L_bomex     ! Logical for BOMEX set up
      LOGICAL :: L_spec_z0   ! specification of roughness length    

      COMMON  /RUN_Ideal/                                               &
     &  h_o, h_o_actual, h_o_per_step,                                  &
     &  lambda_fraction, phi_fraction, half_width_x, half_width_y,      &
     &  Witch_power, plat_size_x, plat_size_y,                          &
     &  height_domain, delta_x, delta_y, big_factor, mag,               &
     &  first_theta_height, thin_theta_height, p_surface,               &
     &  theta_surface, dtheta_dz1, height_dz1, Brunt_Vaisala,           &
     &  u_in, v_in, height_u_in, u_ramp_start, u_ramp_end, q1,          &
     &  ujet_lat, ujet_width,                                           &
     &  t_horizfn_number, t_horizfn_data, uv_horizfn_number,            &
     &  u_ref, v_ref, theta_ref, exner_ref, rho_ref, q_ref,             &
     &  z_orog_print, grow_steps,                                       &
     &  surface_type, grid_number, grid_flat,                           &
     &  tprofile_number, qprofile_number, uvprofile_number,             &
     &  num_profile_data, num_uvprofile_data,                           &
     &  tforce_option, qforce_option, uvforce_option,                   &
     &  num_tforce_levels, num_tforce_times,                            &
     &  num_qforce_levels, num_qforce_times,                            &
     &  num_uvforce_levels, num_uvforce_times,                          &
     &  L_pforce, pforce_option, num_pforce_times,                      &
     &  first_constant_r_rho_level_new,                                 &
     &  big_layers, transit_layers, mod_layers,                         &
     &  zprofile_data, tprofile_data, qprofile_data,                    &
     &  z_uvprofile_data, uprofile_data, vprofile_data,                 &
     &  tforce_time_interval, qforce_time_interval,                     &
     &  uvforce_time_interval, newtonian_timescale,                     &
     &  z_tforce_data, z_qforce_data, z_uvforce_data,                   &
     &  tforce_data, qforce_data, uforce_data, vforce_data,             &
     &  tforce_data_modlev, qforce_data_modlev,                         &
     &  uforce_data_modlev, vforce_data_modlev,                         &
     &  pforce_time_interval, p_surface_data,                           &
     &  L_initialise_data,                                              &
     &  L_perturb_t, perturb_magnitude_t,                               &
     &  L_perturb_q, perturb_magnitude_q,                               &
     &  L_perturb_correlate_tq,                                         &
     &  L_perturb_correlate_vert,                                       &
     &  L_perturb_correlate_time,                                       &
     &  perturb_type, perturb_height,                                   &
     &  L_constant_dz, L_polar_wind_zero,                               &
     &  L_wind_balance, L_rotate_winds,                                 &
     &  IdlSurfFluxSeaOption, IdlSurfFluxSeaParams,                     &
     &  L_spec_z0, roughlen_z0m, roughlen_z0h,                          &
     &  L_pressure_balance, L_vert_Coriolis,                            &
     &  cool_rate, L_force, L_force_lbc,                                &
     &  zprofile_orog, idl_interp_option, hf,                           &
     &  L_fix_orog_hgt_lbc, orog_hgt_lbc,                               &
     &  L_trivial_trigs, f_plane, ff_plane, r_plane,                    &
     &  idl_bubble_option, idl_bubble_max                               &
     &, idl_bubble_height, idl_bubble_width, idl_bubble_depth           &
     &, idl_bubble_xoffset,idl_bubble_yoffset                           &
     &, L_idl_bubble_saturate,                                          &
     &  L_rotating, L_fixed_lbcs, L_code_test,                          &
     &  L_baroclinic, L_cyclone,                                        &
     &  L_damp, L_geo_for, L_bomex,                                     &
     &  DMPTIM, HDMP, ZDMP,                                             &
     &  u_geo, v_geo
! CRUNTIMC end
! Description: COMDECK containing problem_number
!  for use in setting problem types
!
! Author : T. Davies
! History:
! Version  Date      Comment.
! 5.3      15/11/01  New code

      INTEGER, PARAMETER :: standard=0
      INTEGER, PARAMETER :: monsoon=1
      INTEGER, PARAMETER :: dynamical_core=2
      INTEGER, PARAMETER :: idealised_problem=3
      INTEGER, PARAMETER :: standard_namelist=4
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
! ----------------------- Comdeck: NATFORCE  ----------------------------
! Description: COMDECK containing a common block for natural forcing filenames
!
! Author : C.D.Jones
!
! History:
! Version  Date      Comment.
!  6.6.2  11/06/09  New code. C.D.Jones

      CHARACTER*(120) file_scvary
      CHARACTER*(120) file_volcts

      COMMON /FILENATFORCE/                                                 &
     &  file_scvary, file_volcts
      NAMELIST /FILENATFORCE/                                               &
     &  file_scvary, file_volcts

      INTEGER                                                           &
     &       LENRIMA(Nfld_max,NHalo_max,Nrima_max),                     &
                                                     ! IN LBC data len
     &       LBC_SIZEA(4,Nfld_max,NHalo_max,Nrima_max),                 &
                                                       ! IN LBC size
     &       LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max),                &
                                                         ! IN LBC start
     &       RIMWIDTHA(Nrima_max),                                      &
                                    ! IN RIM width
     &       global_row_length,                                         &
                                 ! IN total number of point in a row
     &       global_rows,                                               &
                           ! IN total number of rows in model
     &       MODEL_LEVELS,ROWS,N_ROWS,ROW_LENGTH,ICODE,                 &
     &       LAND_FIELD,                                                &
                           ! IN Number of land points in model from umui
     &       WET_LEVELS,                                                &
                           ! IN Number of wet levels in model, from umui
     &       first_constant_r_rho_level,                                &
                                         !(IN) 1st constant height level
     &       boundary_layer_levels,                                     &
                                    ! (IN) Num. of boundary layer levels
     &       height_gen_method,                                         &
                                  ! (IN) Method for generating heights
     &       CLOUD_LEVELS                                               &
                           ! IN No of cloudy levels in the model
     &      ,tr_levels, tr_vars, tr_ukca                                &
     &      ,yvolc, mvolc, ysol

      INTEGER RowDepCStart
                        ! IN Start of Row dependent constant
      CHARACTER*(80)                                                    &
     &       CMESSAGE              ! Error message if ICODE >0
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='SETCONA')

      REAL                                                              &
            ! Input arguments (grid constants)
     &       Delta_lambda_in                                            &
                             ! EW (x) grid spacing in degrees
     &      ,Delta_phi_in                                               &
                             ! NS (y) grid spacing in degrees
     &      ,Base_phi_in                                                &
                             ! Latitude of first theta point in degrees
     &      ,Base_lambda_in                                             &
                             ! Longitude of first theta point in degs
     &      ,lat_rot_NP_in                                              & 
                           ! Real latitude of 'pseudo' N pole in degs
     &      ,long_rot_NP_in                                             & 
                           ! Real longitude of 'pseudo' N pole in degs
     &      ,z_top_of_model                                             &
                            ! (IN) Height of top of model in metres
     &,      RIMWEIGHTSA(RIMWIDTHA(rima_type_norm))  ! IN RIM weights

      REAL                                                              &
            ! Input VarRes grid info in degrees
     &   Lambda_p_in(global_row_length)                                 &
                                               !IN EW and NS VarRes grid
     &  ,Lambda_u_in(global_row_length)                                 &
                                               !IN EW and NS u,v grid
     &  ,Phi_p_in(global_rows)                                          &
                                               !location in degrees
     &  ,Phi_v_in(global_rows)                 !location in degrees

      REAL                                                              &
            ! Array arguments with intent(IN):
     &     Smvcst(land_field)                                           &
                                      ! IN Volumetric saturation point
     &    ,eta_theta_in(0:model_levels)                                 &
                                        ! IN Eta values for theta levs
     &    ,eta_rho_in(model_levels)                                     &
                                        ! IN Eta values for rho levels
     &    ,Orog(row_length, rows)                                       &
                                        ! IN Orography (on all points)
     &    ,grad_x(row_length*rows)                                      &
                                        ! IN Orographic X-gradient
     &    ,grad_y(row_length*rows)                                      &
                                        ! IN Orographic Y-gradient
     &    ,rho(1-offx:row_length+offx, 1-offy:rows+offy,                &
     &       model_levels)                                              &
                           ! density*(r**2): used in call to Calc_P_star
     &    ,exner_rho_levels(1-offx:row_length+offx,                     &
                                                     ! Exner at p levels
     &                      1-offy:rows+offy, model_levels+1)           &
     &,      orog_lbc(LENRIMA(fld_type_p,halo_type_extended,            &
     &                        rima_type_orog))                          &
                                                 ! Orography LBC
     &,      exner_lbc(LENRIMA(fld_type_p,halo_type_extended,           &
                                                              !Exner LBC
     &                         rima_type_norm),MODEL_LEVELS+1)

      REAL                                                              &
            ! Output args (secondary arrays calculated from dump)
     &  p_theta_levels(1-offx:row_length+offx,                          &
                                               ! press on theta levs Pa
     &                 1-offy:rows+offy, model_levels)                  &
     &, p(1-offx:row_length+offx, 1-offy:rows+offy, model_levels+1)     &
                                                              ! in Pa
     &, p_star(row_length, rows)                                        &
                                 ! surface pressure (Pa)
     &, exner_theta_levels(1-offx:row_length+offx,                      &
                                                   ! Exner at theta levs
     &                     1-offy:rows+offy, model_levels)

      LOGICAL                                                           &
     &       Land_sea_mask(row_length,rows)

      REAL                                                              &
     &  orog_halos(1-halo_i:ROW_LENGTH+halo_i,                          &
     &           1-halo_j:ROWS+halo_j)

      LOGICAL       isSTASH  ! FLUME-STASH

!-----------------Start of TYPCONA------------------------------------
! History:
! Version  Date    Comment
!  3.4   18/05/94  Add new field sin_u_latitude. J F Thomson.
!  5.0   18/05/99  Removed B-grid dynamics terms, replaced by C-P C-grid
!                  dynamics terms. M.L.Gallani
!  5.1   25/02/00  Add minor trig variables to allow initialisation
!                  in Setcona instead of Atm_step. R Rawlins
!  5.3   01/10/01  Add fields for chequerboard radiation. S Cusack
! 6.1  04/08/04  Add arrays  diff_coeff_u, diff_coeff_v
!                                             Terry Davies
! 6.2  04/08/04  Add true_latitude. Yongming Tang
! 6.2  14/02/06  Add arrays for variable resolution  Terry Davies
! CMAXSIZE should be included first.

      ! Constants for ATMOSPHERE model.
      !  Constants for routines independent of resolution.
! CCONSTS start
! Description:
!    This file contains declarations for derived constants within
!   the atmospheric model. Where necessary PARAMETERS are defined to
!   dimension these constants. All constants are placed in the common
!   block CDERIVED, except hardwired constants, e.g. ETA_SPLIT and LENs.
!   file CMAXSIZE must be included first
!
!   The derived constants are calculated in the routine SETCONA1.
!
! PA, WI      <- programmer of some or all of previous code or changes
!
!  Model            Modification history from model version 3.0:
! version  Date
! 3.2   26/03/93  Remove resolution dependent variables for dynamic
!                 allocation. R.Rawlins
!   4.0   20/07/95  Sizes of tables expanded for version 3A
!                   of the radiation in sections 1 and 2. J.M. Edwards
!   5.0   21/06/99  Remove obsolete constants, for C-P C-grid dynamics.
!                   M L Gallani
!   5.1   25/02/00  Replace Data ETA_SPLIT by hard-wire h_split in
!                   Readlsta. Also remove obsolete references to
!                   radiation tables. R Rawlins
! Logical component: F011

      ! No of cloud types ie low/med/high
      INTEGER, PARAMETER :: NUM_CLOUD_TYPES = 3

      ! derived constants:
      INTEGER :: LOW_BOT_LEVEL      ! Bottom level of lowest cloud type
      INTEGER :: LOW_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: MED_BOT_LEVEL      ! Bottom   "    "  med      "    "
      INTEGER :: MED_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: HIGH_BOT_LEVEL     ! Bottom   "    "  top      "    "
      INTEGER :: HIGH_TOP_LEVEL     ! Top      "    "   "       "    "

      ! height values to split model levels into l/m/h cloud
      REAL ::    h_split(NUM_CLOUD_TYPES+1)

      LOGICAL :: ELF                ! T if atmosphere model on LAM grid

      ! Constants for dynamics output independent of resolution but
      ! dependent on choice of levels for output.
      REAL :: REQ_THETA_PV_LEVS(MAX_REQ_THPV_LEVS)

      COMMON /CDERIVED/                                                 &
     &  h_split,LOW_BOT_LEVEL,LOW_TOP_LEVEL,MED_BOT_LEVEL,MED_TOP_LEVEL,&
     &  HIGH_BOT_LEVEL, HIGH_TOP_LEVEL,ELF,REQ_THETA_PV_LEVS
! CCONSTS end

! typcona.h contained constants for the atmosphere.
! As of vn6.6 these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, AD_MASK_TROP_MOD, ROT_COEFF_MOD

! The following common block does not correpsond to the constants
! specified by (and subsequently commented out from) argcona.h, 
! and so it has remained here - 
! It could be moved to a more appropriate place in future

      ! Trigonometric co-ordinates in radians
      REAL:: Delta_lambda       ! EW (x) grid spacing in radians
      REAL:: Delta_phi          ! NS (y) grid spacing in radians
      REAL:: Base_phi           ! Lat of first theta point in radians
      REAL:: Base_lambda        ! Long of first theta point in radians
      REAL:: lat_rot_NP         ! Real lat of 'pseudo' N pole in radians
      REAL:: long_rot_NP        ! Real long of 'pseudo' N pole in radians

      COMMON/cderv_trig/                                                &
     &  Delta_lambda,Delta_phi,Base_phi,Base_lambda,                    &
     &  lat_rot_NP,long_rot_NP
!-------------End of TYPCONA---------------------------------------
! TYPLNDM
! Formerly integral part of TYPCONA, the variables below have been
! separated from the rest of TYPCONA as they are required by some
! of the Ocean routines in the Ocean-Atmosphere configuration of
! the UM whilest TYPCONA is not.

      ! Primary Arrays
      INTEGER::land_points     ! No. of land points  (can be 0)
      INTEGER::land_ice_points ! Number of land ice points
      INTEGER::soil_points     ! Number of soil points

!     INTEGER::land_index (land_field) ! set from land_sea_mask
!     INTEGER::land_ice_index (land_field)  ! Array of land ice points.
!     INTEGER::soil_index(land_field)       ! Array of soil points.
! sza fix conflict case when land_field=0
      INTEGER::land_index (max(1,land_field)) ! set from land_sea_mask
      INTEGER::land_ice_index (max(1,land_field))  ! Array of land ice points.
      INTEGER::soil_index(max(1,land_field))       ! Array of soil points.

      ! Gets some sizes transported around the model :
      COMMON /land_soil_dimensions/                                     &
     &  land_points , land_ice_points , soil_points

! TYPLNDA end
!*L -----------------COMDECK PHYSCONS----------------------------------
!
!  Purpose : contains physical constants required by the whole of the
!            model. It is made up of individual COMDECKS for sets of
!            of related constants, each routine can access one or
!            several of these COMDECKS seperately
!  System Component : F07
!  System task : Z
! END
!*----------------------------------------------------------------------
!*L------------------COMDECK C_A----------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Replace variable A by more meaningful name for
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Convert to Fixed/Free format. P. Selwood

      ! Mean radius of Earth in m.
      Real, Parameter  :: Earth_Radius = 6371000.

!*----------------------------------------------------------------------
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
! C_LHEAT start

! latent heat of condensation of water at 0degc
      REAL,PARAMETER:: LC=2.501E6

 ! latent heat of fusion at 0degc
      REAL,PARAMETER:: LF=0.334E6

! C_LHEAT end
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
! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------
!*L------------------COMDECK C_OMEGA------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable two_omega for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.3   12/10/01  two_omega initialised in SETCON     Terry Davies
!OMEGA IS MAGNITUDE OF EARTH'S ANGULAR VELOCITY
      REAL OMEGA                                                        &
     &,two_omega

       Common/Omega/Omega, two_omega
!  Angular speed of Earth's rotation Omega to be initialised in SETCON
!*----------------------------------------------------------------------
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
! C_KT_FT start

      REAL,PARAMETER:: KT2MS=1852.0/3600.0 ! Knots to m/s conversion
      REAL,PARAMETER:: FT2M =0.3048        ! Feet to meters conversion

! C_KT_FT end
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
!*L------------------COMDECK C_LAPSE ----------------------------------
      Real, Parameter :: lapse      = 0.0065  ! Near surface lapse rate
      Real, Parameter :: lapse_trop = 0.002   ! Tropopause lapse rate
!*---------------------------------------------------------------------
! C_ETA_PMSL start
! ETA_PMSL is the ETA value used to determine the model level
! used in PMSL reduction calculations.
      REAL,PARAMETER:: ETA_PMSL=0.795
      INTEGER LEVNO_PMSL_CALC   !  Model level for PMSL reduction calc.
      COMMON /PMSLCALC/ LEVNO_PMSL_CALC
! C_ETA_PMSL end
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
!
!  4.1      04/09/96 New comdeck for Data Assimilation on T3E
!                                   Deborah Salmond
!  4.3      18/4/97 Increase OBSNUMDIM  Stuart Bell
!  4.4      19/11/97 Increase OBSNUMDIM,OBSDIM,inxdim  Stuart Bell
!  5.2      12/12/00 Restore COUNTA,B,C dropped at 5.0   B Macpherson
!LL  5.0      24/06/99 Alter names for C-P C-grid dynamics. M.L.Gallani
!LL  5.3      08/06/01 Remove duplicate declarations.  A van der Wal
!    5.3      05/12/01 Re-size obsdim to smaller value (500,000 to
!                      250,000) to reduce memory usage.    S. Cusack
!    6.0     17/10/03 Increase obsdim and obsnumdim for SX6. Clive Jones
!    6.2     11/01/06 Remove hard-wired obsdim and obsnumdim.
!                       Camilla Mathison/R Barnes

      ! dimensions for obs allocated in subroutine AC
      ! dynamically allocated in AC and RDOBS

      ! dimension inxdim allocated in subroutine HORINF
      INTEGER,PARAMETER:: inxdim    = 15000

      ! common for Statistics Calcs in DIAGO ; Prints in RDOBS,GETOBS
      REAL :: R_STAT(MODEL_LEVELS_MAX,0:8)
      REAL :: S_STAT(MODEL_LEVELS_MAX,0:8)
      INTEGER :: COUNTA(NOBTYPMX)
      INTEGER :: COUNTB(NOBTYPMX)
      INTEGER :: COUNTC(NOBTYPMX)

      COMMON /mpp_ac/ R_STAT,S_STAT,COUNTA,COUNTB,COUNTC

      ! common to pass longitudes and latitudes for edges of local
      ! box from setcona to RDOBS and HINTCF

      REAL :: LONG_E
      REAL :: LONG_W
      REAL :: LAT_N,LAT_S
      REAL :: LONG_W_MODEL
      REAL :: LONG_E_MODEL
      COMMON/latlonmax/                                                 &
     &  LONG_E,LONG_W,LAT_N,LAT_S,LONG_W_MODEL,LONG_E_MODEL
! MPPAC end
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
!+ ---------------------------------------------------------------------
!  Description: Limits on the position of the tropopause used in
!               radiative calculations are defined.
!
!  Current Code Owner: J. M. Edwards
!
!  History:
!
!  Version  Date      Comment.
!  5.3      27/09/01  Original code.
!                     (J. M. Edwards)
!  6.2      27/01/06  Allows limits to be set in SCM (R. Wong)
!
!- ----------------------------------------------------------------------
!
!     Set limits on the height of the tropopause. Generous limits are
!     required as the tropopause can be very low in the Antarctic
!     winter. These limits should be reviewed for simulations of
!     climates very different from the present one, such as runs
!     with very high concentrations of CO2.
!
!     These limits are in part chosen to accord with the pressure
!     levels used in earlier configurations of the Unified Model.
!     In setting the upper limit consideration has been given to the
!     paper entitled "The tropical tropopause over the Western
!     Pacifiic: Wave driving, convection and the annual cycle,"
!     (J. Geophys. Res., 1996, Vol. 101 (D16), p. 21223) by
!     G. C. Reid and K. S. Gage.
      Real, parameter :: z_min_trop = 2.0e3 ! Lowest permitted height
      Real, parameter :: z_max_trop = 2.0e4 ! Maximum permitted height
!
! -----------------------------------------------------------------------

! HIGHOS starts:
! contains the allowed high order scheme options
!
! Author : Michail Diamantakis
! History:
! Version  Date      Comment.
! 5.3      25/10/01  New comdeck
!
      INTEGER,PARAMETER :: cubicLagrange         = 1
      INTEGER,PARAMETER :: quinticLagrange       = 2
      INTEGER,PARAMETER :: ECMWF_quasiCubic      = 3
      INTEGER,PARAMETER :: ECMWF_mono_quasiCubic = 4
      INTEGER,PARAMETER :: hCubic_vLin           = 5
      INTEGER,PARAMETER :: hQuasiCubic_vQuintic  = 6
      INTEGER,PARAMETER :: hCubic_vQuintic       = 7
! HIGHOS ends


! Local variables
      REAL                                                              &
     &   Delta_lambda_wk                                                &
                             ! EW (x) grid spacing in degrees
     &  ,Delta_phi_wk                                                   &
                             ! NS (y) grid spacing in degrees
     &  ,Base_phi_wk                                                    &
                             ! Latitude of first theta point in degrees
     &  ,Base_lambda_wk                                                 &
                             ! Longitude of first theta point in degrees
     &  ,latitude_rot(row_length, rows)                                 &
                                         ! rot. latit. in degrees (LAM)
     &  ,longitude_rot(row_length, rows)                                &
                                         ! rot. longit. in degrees (LAM)
     &  ,true_latitude_b(row_length, rows)                              &
                                           ! true latitude on B-grid
     &  ,true_longitude_b(row_length,rows)                              &
                                           ! true longitude on B-grid
     &  ,bear_rot_NP(row_length, rows)                                  & 
                                         ! Bearing of 'pseudo' N pole
     &  ,temp1                                                          &
                   ! Used in calculation of Coriolis terms (LAM only)
     &  ,temp2                                                          &
                   ! Used in calculation of Coriolis terms (LAM only)
     &  ,constant                                                       &
                   ! 1/kappa used to calc pressure on model levels
     &   , f_plane_rad                                                  &
                        ! f_plane latitude in radians
     &   , ff_plane_rad                                                 &
                         ! f_plane latitude in radians
     &   , f1_temp                                                      &
                    ! f1 term for f-plane or ff-plane
     &   , f2_temp                                                      &
                    ! f2 term for f-plane or ff-plane
     &   , f3_temp                                                      &
                    ! f3 term for f-plane or ff-plane
     &  , scale                                                         &
     &  ,CLOUD_BOUND(NUM_CLOUD_TYPES+1)                                 &
                                        ! boundaries of cloud types
     &  ,r_ref_theta(model_levels)                                      &
                                   ! Local dynamic array
     &  ,r_ref_rho(model_levels)   ! Local dynamic array

      LOGICAL                                                           &
     &   landmask(1-offx:row_length+offx, 1-offy:rows+offy)             &
     &  ,L_error                                                        &
     &  ,L_varres

      INTEGER                                                           &
              ! Mostly loop counters, but j also used for interp points
     &   I,KK,                                                          &
     &   J,GI,GJ,                                                       &
     &   j0,j1,k,                                                       &
     &   LEVEL                                                          &
                         ! Used to set up cloud type boundaries
     &  ,first_row                                                      &
     &  ,last_row                                                       &
     &  ,info                                                           &
     &  , active_levels

      INTEGER, PARAMETER :: height_gen_original = 1 ! methods for height
      INTEGER, PARAMETER :: height_gen_smooth   = 2 ! generation

!  External subroutines called :
      EXTERNAL                                                          &
     &   EQTOLL                                                         &
     &  ,W_Coeff                                                        &
     &  ,Ereport                                                        &
     &  ,Calc_Exner_at_theta                                            &
     &  ,Calc_P_from_Exner                                              &
     &  ,Calc_P_star                                                    &
     &  ,Swap_Bounds
      EXTERNAL pole_bearing,aspang

! Set flag for orography correction in module solinc_data
      L_orog = L_use_orog_corr .or. L_use_grad_corr
! Grid type in dump: L_varres = T / F for variable / unif grid
      L_varres = .FALSE.
      If (RowDepCStart >= 1) L_varres = .TRUE.
! ----------------------------------------------------------------------
! 0. Set-up vertical co-ordinate arrays
! ----------------------------------------------------------------------
! Set up orog_tmp - copy of orog but with wide halos
      DO j=1-halo_j,ROWS+halo_j
        DO i=1-halo_i,ROW_LENGTH+halo_i
          orog_halos(i,j)=0.0
        END DO
      END DO

      DO j=1,ROWS
        DO i=1,ROW_LENGTH
          orog_halos(i,j)=orog(i,j)
        END DO
      END DO

! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(orog_halos,row_length,rows,1,                    &
     &                 halo_i,halo_j,fld_type_p,.FALSE.)

! ----------------------------------------------------------------------
! 0.1 Write out some useful information at start of run
! ----------------------------------------------------------------------

      write(6,*) '  '
      write(6,*) '  ****  This run uses bit-reproducible code  ****'

      write(6,*) '  '
      write(6,*) '  ****  PE configuration for this run  ****   '
      if ( nproc_x > 1 .and. nproc_y > 1 ) then
        write(6,*) '      ', nproc_x, ' processors East-West ',         &
     &                     nproc_y, ' processors North South'
      else if ( nproc_x > 1 .and. nproc_y == 1 ) then
        write(6,*) '      ', nproc_x, ' processors East-West ',         &
     &                     nproc_y, ' processor North South'
      else if ( nproc_x == 1 .and. nproc_y > 1 ) then
        write(6,*) '      ', nproc_x, ' processor East-West ',          &
     &                     nproc_y, ' processors North South'
      else ! nproc_x = 1 .and. nproc_y =1 ) then
        write(6,*) '       Single processor 1x1 '
      end if ! nproc_x > 1 .and. nproc_y > 1 )

      if ( offx >= halo_i .or. offy >= halo_j ) then
        if (.not.isSTASH ) then                 ! FLUME-STASH
          write(6,*) '  '
          write(6,*) '  ****  DANGER  ****   '
          write(6,*) '  *** SMALL halo >= LARGE halo  *****   '
          write(6,*) '  This could result in overwriting  '
          ICODE    = 20
          CMESSAGE ='SETCONA: Large halo size is too small'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        end if
      end if ! offx >= halo_i .or. offy >= halo_j

      If (L_mix_ratio) Then
        write(6,*) '  '
        write(6,*) '  ***  This run uses mixing ratios for      ***'
        write(6,*) '  ***  the moist variables in the dynamics  ***'
      Else
        write(6,*) '  '
        write(6,*) '  ***  This run uses specific humidities for  ***'
        write(6,*) '  ***  the moist variables in the dynamics    ***'
      End If !  L_mix_ratio

      If (model_domain  ==  mt_LAM) Then
      
        write(6,*) '  '
        Write ( Unit = 6, fmt=*) ' ***   LAM set-up ***' 
        
        Write ( Unit = 6, fmt=*) 'LBC frame size for solver, '          &
     &,                          ' n_rims_to_do = ', n_rims_to_do
        If( L_LBC_balance ) then
          Write ( Unit = 6, fmt=*) 'L_LBC_balance = ', L_LBC_balance    &
     &,          ' Impose vertically balanced Exner pressures'          &
     &,          ' and rho and set w=0 in lbcs '
        Else !  L_LBC_balance = .false.
          Write ( Unit = 6, fmt=*) 'L_LBC_balance = ', L_LBC_balance    &
     &,                            ' No balancing of lbcs '
        EndIf ! L_LBC_balance      
        If( L_lbc_new ) then
          Write ( Unit = 6, fmt=*) 'L_lbc_new = ', L_lbc_new            &
     &,                              'Use new lbc algorithm '
        Else !  L_lbc_new = .false.
          Write ( Unit = 6, fmt=*) 'L_lbc_new = ', L_lbc_new            &
     &,                              'Use old lbc algorithm '
        EndIf ! L_lbc_new
        If( L_int_uvw_lbc) then
          Write ( Unit = 6, fmt=*) 'L_int_uvw_lbc =', L_int_uvw_lbc     &
     &,          ' Advecting winds interpolated in lateral boundaries '
        Else !  L_int_uvw_lbc = .false.
          Write ( Unit = 6, fmt=*) 'L_int_uvw_lbc =', L_int_uvw_lbc     &
     &,                ' Extrapolated advecting winds from lbc file '
        EndIf ! L_int_uvw_lbc
        write(6,*) '  '

! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
     &    ROW_LENGTH,ROWS,halo_i,halo_j,1,fld_type_p,orog_halos,        &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_orog),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_orog),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_orog),   &
     &    halo_i,halo_j,orog_lbc,RIMWIDTHA(rima_type_orog),             &
     &    RIMWIDTHA(rima_type_orog),RIMWEIGHTSA,                        &
     &    at_extremity,                                                 &
     &    .FALSE.,.TRUE.)

! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
     &    ROW_LENGTH,ROWS,Offx,Offy,MODEL_LEVELS+1,fld_type_p,          &
     &    exner_rho_levels,                                             &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm),   &
     &    halo_i,halo_j,exner_lbc,RIMWIDTHA(rima_type_norm),            &
     &    RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                        &
     &    at_extremity,                                                 &
     &    .FALSE.,.TRUE.)

      END IF ! IF (model_domain  ==  mt_LAM)

! Set reference height profile

!!! Allocate arrays for level_heights_mod module
      If (.not.allocated(eta_theta_levels)) Then
        Allocate (eta_theta_levels(0:model_levels))
      End If
      If (.not.allocated(eta_rho_levels)) Then
        Allocate (eta_rho_levels(model_levels))
      End If
      If (.not.allocated(r_theta_levels)) Then
        Allocate (r_theta_levels(1-halo_i:row_length+halo_i,           &
     &                           1-halo_j:rows+halo_j, 0:model_levels))
      End If
      If (.not.allocated(r_rho_levels)) Then
        Allocate (r_rho_levels(1-halo_i:row_length+halo_i,             &
     &                         1-halo_j:rows+halo_j, model_levels))
      End If

      eta_theta_levels(0) = eta_theta_in(0)
      Do k = 1, model_levels
        eta_theta_levels(k) = eta_theta_in(k)
        eta_rho_levels(k) = eta_rho_in(k)
        r_ref_theta(k) = eta_theta_levels(k) * z_top_of_model
      End Do
      Do k = 1, model_levels
        r_ref_rho(k) = eta_rho_levels(k) * z_top_of_model
      End Do
! set bottom level, ie: orography
      Do j = 1-halo_j, rows+halo_j
        Do i= 1-halo_i, row_length+halo_i
          r_theta_levels(i,j,0) = Orog_halos(i,j) + Earth_radius
        End Do
      End Do
! For constant levels set r to be a constant on the level
      Do k = first_constant_r_rho_level, model_levels
        Do j = 1-halo_j, rows+halo_j
          Do i= 1-halo_i, row_length+halo_i
            r_theta_levels(i,j,k) = Earth_radius + r_ref_theta(k)
            r_rho_levels(i,j,k) = Earth_radius + r_ref_rho(k)
          End Do
        End Do
      End Do

      Select Case( height_gen_method )
        Case( height_gen_original )
! The original version of height generation used in the SI dynamics
!
! For boundary layer levels set depth to be constant.
      Do k = 1, boundary_layer_levels
        Do j = 1-halo_j, rows+halo_j
          Do i= 1-halo_i, row_length+halo_i
            r_theta_levels(i,j,k) = r_theta_levels(i,j,0)               &
     &                                 + r_ref_theta(k)
            r_rho_levels(i,j,k) = r_theta_levels(i,j,0)                 &
     &                               + r_ref_rho(k)
          End Do
        End Do
      End Do
! For intermediate levels use linear relaxation to constant value.
! set orographic heights.
      Do k = boundary_layer_levels+1,                                   &
     &       first_constant_r_rho_level-1
        Do j = 1-halo_j, rows+halo_j
          Do i= 1-halo_i, row_length+halo_i
              r_rho_levels(i,j,k) =                                     &
     &          ( r_rho_levels(i,j,first_constant_r_rho_level) -        &
     &            r_theta_levels(i,j,boundary_layer_levels) ) *         &
     &          ( eta_rho_levels(k) -                                   &
     &            eta_theta_levels(boundary_layer_levels) ) /           &
     &          ( eta_rho_levels(first_constant_r_rho_level) -          &
     &            eta_theta_levels(boundary_layer_levels) )             &
     &           +  r_theta_levels(i,j,boundary_layer_levels)
            r_theta_levels(i,j,k) =                                     &
     &          ( r_rho_levels(i,j,first_constant_r_rho_level) -        &
     &            r_theta_levels(i,j,boundary_layer_levels) ) *         &
     &          ( eta_theta_levels(k) -                                 &
     &            eta_theta_levels(boundary_layer_levels) ) /           &
     &          ( eta_rho_levels(first_constant_r_rho_level) -          &
     &            eta_theta_levels(boundary_layer_levels) )             &
     &           +  r_theta_levels(i,j,boundary_layer_levels)
          End Do
        End Do
      End Do

        Case( height_gen_smooth )
! A smooth quadratic height generation
          Do k = 1, first_constant_r_rho_level-1
            Do j = 1-halo_j, rows+halo_j
              Do i= 1-halo_i, row_length+halo_i
              r_rho_levels(i,j,k) = eta_rho_levels(k) * z_top_of_model +&
     &         Earth_radius + Orog_halos(i,j) * (1.0 - eta_rho_levels(k)&
     &              /eta_rho_levels(first_constant_r_rho_level))**2
              r_theta_levels(i,j,k) = eta_theta_levels(k) *             &
     &             z_top_of_model + Earth_radius + Orog_halos(i,j) *    &
     &             (1.0 - eta_theta_levels(k) /                         &
     &              eta_rho_levels(first_constant_r_rho_level))**2
              End Do
            End Do
          End Do

        Case Default
          icode = 10
          Write (Cmessage,*) 'Unrecognised height generation method - ',&
     &                       'Dump needs to be reconfigured'
! DEPENDS ON: ereport
          Call Ereport( RoutineName, icode, Cmessage )
      End Select

! 0.1 Initialise secondary arrays.
! Exner at p (=rho) levels is obtained from dump.
! calculate p from exner_rho_levels
! [halos required for diagnostic calculations at T+0 in INITDIAG.]
      constant = 1./ kappa
      Do k = 1, model_levels+1
        Do j = 1-offy, rows+offy
          Do i = 1-offx, row_length+offx
            p(i,j,k)= (exner_rho_levels(i,j,k) ** constant)             &
     &                      * p_zero
          End Do
        End Do
      End Do

! calculate exner at theta_levels which is then used to get
! p at theta_levs
! DEPENDS ON: calc_exner_at_theta
      Call Calc_Exner_at_theta( r_theta_levels, r_rho_levels,           &
     &                exner_rho_levels,                                 &
     &                 row_length, rows, model_levels,                  &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 exner_theta_levels,.TRUE.)

! Calculate pressure from Exner at theta levels.

! DEPENDS ON: calc_p_from_exner
       call Calc_P_from_Exner(                                          &
     &                      p_theta_levels, kappa, p_zero,              &
     &                      row_length, rows, model_levels,             &
     &                      offx, offy,                                 &
     &                      exner_theta_levels,.TRUE.)

! calculate p_star using rho (from dump) and p on model levels
! DEPENDS ON: calc_p_star
      Call Calc_P_star (r_theta_levels, r_rho_levels, p, rho,           &
     &                  g, row_length, rows, model_levels,              &
     &                  offx, offy, halo_i, halo_j,                     &
     &                  p_star)

! ----------------------------------------------------------------------
! Initialise Constants and trigonometric functions.
! ----------------------------------------------------------------------

      If( L_rotating) Then
        OMEGA = 7.292116E-5 ! Angular speed of Earth's rotation
                            ! = 2*pi/siderial day (23h56m04s)
      Else
        OMEGA = 0.0         ! Zero angular speed of rotation for planet
      End If    ! L_rotating
      two_omega = 2. * OMEGA

! Convert grid-spacing and base values from degrees to radians
      If( L_regular ) then ! 
        base_lambda_wk = base_lambda_in
        base_phi_wk = base_phi_in
        delta_lambda_wk = delta_lambda_in
        delta_phi_wk = delta_phi_in
      Else
        base_lambda_wk = lambda_p_in(1)
        base_phi_wk = phi_p_in(1)
        delta_lambda_wk = (lambda_p_in(global_row_length) -             &
     &                    lambda_p_in(1))/(global_row_length - 1) 
        delta_phi_wk = (phi_p_in(global_rows) -                         &
     &                    phi_p_in(1))/(global_rows - 1) 
      End If 
      delta_lambda=delta_lambda_wk * Pi_over_180
      delta_phi=delta_phi_wk * Pi_over_180
      base_phi=base_phi_wk * Pi_over_180
      base_lambda=base_lambda_wk * Pi_over_180
      f_plane_rad = f_plane * Pi_over_180
      lat_rot_NP =lat_rot_NP_in  * Pi_over_180
      long_rot_NP=long_rot_NP_in * Pi_over_180

      If ( L_regular ) Then 

        Write ( Unit = 6, fmt=*) ' '
        Write ( Unit = 6, fmt=*) '*** The horizontal grid is regular'   &
     &                    , '  L_regular = ', L_regular, ' ***'
        Write ( Unit = 6, fmt=*) ' '

!   allocate small arrays to avoid out-of-bounds in atm_step 

        If (.not.allocated(glambda_p)) Then
          Allocate (glambda_p ( 1-halo_i : halo_i ))
        End If
        If (.not.allocated(glambda_u)) Then
          Allocate (glambda_u ( 1-halo_i : halo_i ))
        End If
        If (.not.allocated(gdlambda_p)) Then
          Allocate (gdlambda_p ( 1-halo_i : halo_i ))
        End If
        If (.not.allocated(gdlambda_u)) Then
          Allocate (gdlambda_u ( 1-halo_i : halo_i ))
        End If
        If (.not.allocated(grecip_dlamp)) Then
          Allocate (grecip_dlamp ( 1-halo_i : halo_i ))
        End If
        If (.not.allocated(grecip_dlamu)) Then
          Allocate (grecip_dlamu ( 1-halo_i : halo_i ))
        End If

      Else  !  Variable grid; set up parameters

        Write ( Unit = 6, fmt=*) ' '
        Write ( Unit = 6, fmt=*) '*** This run uses a variable '        &
     &           , ' horizontal grid L_regular = ', L_regular, ' ***'
        Write ( Unit = 6, fmt=*) ' '

! Allocate arrays for dyn_var_res_mod module
        If (.not.allocated(glambda_p)) Then
          Allocate (glambda_p ( 1-halo_i : global_row_length+halo_i ))
        End If
        If (.not.allocated(glambda_u)) Then
          Allocate (glambda_u ( 1-halo_i : global_row_length+halo_i ))
        End If
        If (.not.allocated(phi_p)) Then
          Allocate (phi_p ( 1-halo_i : row_length + halo_i,             &
     &                      1-halo_j : rows + halo_j ))
        End If
        If (.not.allocated(phi_v)) Then
          Allocate (phi_v ( 1-halo_i : row_length + halo_i,             &
     &                      1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(gdlambda_p)) Then
          Allocate (gdlambda_p ( 1-halo_i : global_row_length+halo_i ))
        End If
        If (.not.allocated(gdlambda_u)) Then
          Allocate (gdlambda_u ( 1-halo_i : global_row_length+halo_i ))
        End If
        If (.not.allocated(dphi_p)) Then
          Allocate (dphi_p ( 1-halo_i : row_length + halo_i,            &
     &                       1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(dphi_v)) Then
          Allocate (dphi_v ( 1-halo_i : row_length + halo_i,            &
     &                       1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(grecip_dlamp)) Then
          Allocate (grecip_dlamp ( 1-halo_i : global_row_length+halo_i))
        End If
        If (.not.allocated(grecip_dlamu)) Then
          Allocate (grecip_dlamu ( 1-halo_i : global_row_length+halo_i))
        End If
        If (.not.allocated(recip_dphip)) Then
          Allocate (recip_dphip ( 1-halo_i : row_length + halo_i,       &
     &                            1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(recip_dphiv)) Then
          Allocate (recip_dphiv ( 1-halo_i : row_length + halo_i,       &
     &                            1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(wt_lambda_p)) Then
          Allocate (wt_lambda_p ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(wt_lambda_u)) Then
          Allocate (wt_lambda_u ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(wt_phi_p)) Then
          Allocate (wt_phi_p ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(wt_phi_v)) Then
          Allocate (wt_phi_v ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(lambda_p_rm)) Then
          Allocate (lambda_p_rm ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(lambda_p_rp)) Then
          Allocate (lambda_p_rp ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(lambda_u_rm)) Then
          Allocate (lambda_u_rm ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(lambda_u_rp)) Then
          Allocate (lambda_u_rp ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(phi_p_rm)) Then
          Allocate (phi_p_rm ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(phi_p_rp)) Then
          Allocate (phi_p_rp ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(phi_v_rm)) Then
          Allocate (phi_v_rm ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(phi_v_rp)) Then
          Allocate (phi_v_rp ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(recip_lambda_p_m)) Then
          Allocate (recip_lambda_p_m ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_p_0)) Then
          Allocate (recip_lambda_p_0 ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_p_p)) Then
          Allocate (recip_lambda_p_p ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_p_p2)) Then
          Allocate (recip_lambda_p_p2 ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_u_m)) Then
          Allocate (recip_lambda_u_m ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_u_0)) Then
          Allocate (recip_lambda_u_0 ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_u_p)) Then
          Allocate (recip_lambda_u_p ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_u_p2)) Then
          Allocate (recip_lambda_u_p2 ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_phi_p_m)) Then
          Allocate (recip_phi_p_m ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_p_0)) Then
          Allocate (recip_phi_p_0 ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_p_p)) Then
          Allocate (recip_phi_p_p ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_p_p2)) Then
          Allocate (recip_phi_p_p2 ( 1-halo_i : row_length + halo_i,    &
     &                               1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_v_m)) Then
          Allocate (recip_phi_v_m ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_v_0)) Then
          Allocate (recip_phi_v_0 ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_v_p)) Then
          Allocate (recip_phi_v_p ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_v_p2)) Then
          Allocate (recip_phi_v_p2 ( 1-halo_i : row_length + halo_i,    &
     &                               1-halo_j : n_rows+halo_j ))
        End If

! DEPENDS ON: set_var_grid
      Call Set_var_grid (                                               &
     &                   Lambda_p_in, Lambda_u_in,                      &
     &                   Phi_p_in, Phi_v_in,                            &
     &                   global_row_length, global_rows,                &
     &                   row_length, rows,  n_rows,                     &
     &                   halo_i, halo_j, L_varres,                      &
     &                   delta_lambda_wk, delta_phi_wk,                 &
     &                   base_lambda_wk, base_phi_wk,                   &
     &                   lam_var, phi_var,                              &
     &                   var_ratio, lam_ratio, phi_ratio,               &
     &                   lam_frac, phi_frac, Pi_over_180,               &
     &                   glambda_p, glambda_u, phi_p, phi_v,            &
     &                   gdlambda_p, gdlambda_u, dphi_p, dphi_v,        &
     &                   grecip_dlamp, grecip_dlamu,                    &
     &                   recip_dphip, recip_dphiv,                      &
     &                   wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v,  &
     &                   lambda_p_rm, lambda_p_rp,                      &
     &                   lambda_u_rm, lambda_u_rp,                      &
     &                   phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,        &
     &                   lambda_p_end, lambda_u_end,                    &
     &                   phi_p_end, phi_v_end, dlambda_p_end,           &
     &                   dlambda_u_end, dphi_p_end, dphi_v_end,         &
     &                   recip_lambda_p_m, recip_lambda_p_0,            &
     &                   recip_lambda_p_p, recip_lambda_p_p2,           &
     &                   recip_lambda_u_m, recip_lambda_u_0,            &
     &                   recip_lambda_u_p, recip_lambda_u_p2,           &
     &                   recip_phi_p_m, recip_phi_p_0,                  &
     &                   recip_phi_p_p, recip_phi_p_p2,                 &
     &                   recip_phi_v_m, recip_phi_v_0,                  &
     &                   recip_phi_v_p, recip_phi_v_p2,                 &
     &                   max_look, recip_dlam, recip_dphi,              &
     &                   halo_lam, halo_phi, look_lam, look_phi,        &
     &                   model_domain, datastart )

      End If   !  L_regular

! 1. set trig fields and Coriolis components

!!! Allocate latitude arrays for trignometric_mod module
      If (.not.allocated(cos_theta_latitude)) Then
        Allocate (cos_theta_latitude (1-Offx:row_length+Offx,           &
     &                                1-Offy:rows+Offy))
      End If
      If (.not.allocated(sec_theta_latitude)) Then
        Allocate (sec_theta_latitude (1-Offx:row_length+Offx,           &
     &                                1-Offy:rows+Offy))
      End If
      If (.not.allocated(FV_cos_theta_latitude)) Then
        Allocate (FV_cos_theta_latitude (1-Offx:row_length+Offx,        &
     &                                   1-Offy:rows+Offy))
      End If
      If (.not.allocated(FV_sec_theta_latitude)) Then
        Allocate (FV_sec_theta_latitude (1-Offx:row_length+Offx,        &
     &                                   1-Offy:rows+Offy))
      End If
      If (.not.allocated(sin_theta_latitude)) Then
        Allocate (sin_theta_latitude (row_length, rows))
      End If
      If (.not.allocated(tan_theta_latitude)) Then
        Allocate (tan_theta_latitude (row_length, rows))
      End If
      If (.not.allocated(sin_v_latitude)) Then
        Allocate (sin_v_latitude (row_length, n_rows))
      End If
      If (.not.allocated(tan_v_latitude)) Then
        Allocate (tan_v_latitude (row_length, n_rows))
      End If
      If (.not.allocated(cos_v_latitude)) Then
        Allocate (cos_v_latitude (1-Offx:row_length+Offx,               &
     &                            1-Offy:n_rows+Offy))
      End If
      If (.not.allocated(sec_v_latitude)) Then
        Allocate (sec_v_latitude (1-Offx:row_length+Offx,               &
     &                            1-Offy:n_rows+Offy))
      End If

!!! Allocate longitude and 'true' arrays for trignometric_mod module
      If (.not.allocated(cos_theta_longitude)) Then
        Allocate (cos_theta_longitude (row_length, rows))
      End If
      If (.not.allocated(sin_theta_longitude)) Then
        Allocate (sin_theta_longitude (row_length, rows))
      End If
      If (.not.allocated(cos_u_longitude)) Then
        Allocate (cos_u_longitude (row_length, rows))
      End If
      If (.not.allocated(sin_u_longitude)) Then
        Allocate (sin_u_longitude (row_length, rows))
      End If
      If (.not.allocated(true_latitude)) Then
        Allocate (true_latitude (row_length, rows))
      End If
      If (.not.allocated(true_longitude)) Then
        Allocate (true_longitude (row_length, rows))
      End If

      If(model_domain  /=  mt_global) Then
!  For LAMs only, Trigs may be set to equatorial values
        If (L_trivial_trigs) then
        If(mype  ==  0) write(6,*)'WARNING: trig values set to trivial'
        Do j = 1, rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
            cos_theta_latitude(i, j) = 1.
            sec_theta_latitude(i, j) = 1.
        sin_theta_latitude(i, j) = 0.0
!       sin_theta_latitude(i, j) = sin(f_plane_rad)
            tan_theta_latitude(i, j) = 0.
            FV_cos_theta_latitude(i, j) = 1.
            FV_sec_theta_latitude(i, j) = 1.
          End Do
        End Do

        Do j = 1, n_rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
        sin_v_latitude(i, j) = 0.0
!       sin_v_latitude(i, j) = sin(f_plane_rad)
            cos_v_latitude(i, j) = 1.
            tan_v_latitude(i, j) = 0.
            sec_v_latitude(i, j) = 1.
          End Do
        End Do

        Do j = 1, rows
          Do i = 1, row_length
            gi = datastart(1) + i - 1
            sin_theta_longitude(i, j) = sin ((gi-1)*delta_lambda)
            cos_theta_longitude(i, j) = cos ((gi-1)*delta_lambda)
            sin_u_longitude(i, j) = sin ((gi-.5)*delta_lambda)
            cos_u_longitude(i, j) = cos ((gi-.5)*delta_lambda)
          End Do
        End Do

        End If        ! L_trivial_trigs = .true.
      End If        !model_domain  /=  mt_global
      If (model_domain  ==  mt_global .or. .not.L_trivial_trigs) Then

! non-trivial trigonometric info.
        Do j = 1, rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
            cos_theta_latitude(i, j) = cos(Base_phi+(gj-1)*delta_phi)
            sin_theta_latitude(i, j) = sin(Base_phi+(gj-1)*delta_phi)
            FV_cos_theta_latitude(i, j) = cos_theta_latitude(i, j)
          End Do
        End Do

        If (model_domain  ==  mt_global) Then

          j0 = 1
          j1 = rows
          If (at_extremity(PNorth)) Then
            Do i = 1, row_length
              cos_theta_latitude(i, rows) = 0.
              sin_theta_latitude(i, rows) = 1.
              FV_cos_theta_latitude(i, rows) = delta_phi/8.
              sec_theta_latitude(i, rows) = rmdi
              tan_theta_latitude(i, rows) = rmdi
            End Do

            j1 = rows - 1

          End If

          If (at_extremity(PSouth)) Then
            Do i = 1, row_length
              cos_theta_latitude(i, 1) = 0.
              sin_theta_latitude(i, 1) = -1.
              FV_cos_theta_latitude(i, 1) = delta_phi/8.
              sec_theta_latitude(i, 1) = rmdi
              tan_theta_latitude(i, 1) = rmdi
            End Do

            j0 = 2

          End If


          Do j = j0, j1
            Do i = 1, row_length
              sec_theta_latitude(i, j) = 1./cos_theta_latitude(i,j)
              tan_theta_latitude(i, j) = sin_theta_latitude(i,j)        &
     &                                   /cos_theta_latitude(i,j)
            End Do
          End Do


        Else
! Limited area model
          Do j = 1, rows
            Do i = 1, row_length
              sec_theta_latitude(i, j) = 1./cos_theta_latitude(i, j)
              tan_theta_latitude(i, j) = sin_theta_latitude(i, j)       &
     &                                   /cos_theta_latitude(i, j)
            End Do
          End Do
        End If ! model_domain

        Do j = 1, rows
          Do i = 1, row_length
            FV_sec_theta_latitude(i, j) = 1./FV_cos_theta_latitude(i,j)
          End Do
        End Do

        Do j = 1, n_rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
            sin_v_latitude(i, j) = sin(Base_phi + (gj-.5)*delta_phi)
            cos_v_latitude(i, j) = cos(Base_phi + (gj-.5)*delta_phi)
            sec_v_latitude(i, j) = 1./cos_v_latitude(i, j)
            tan_v_latitude(i, j) = sin_v_latitude(i, j)                 &
     &                             /cos_v_latitude(i, j)
          End Do
        End Do

        Do j = 1, rows
          Do i = 1, row_length
            gi = datastart(1) + i - 1
            sin_theta_longitude(i, j) =                                 &
     &                      sin (Base_lambda+(gi-1)*delta_lambda)
            cos_theta_longitude(i, j) =                                 &
     &                      cos (Base_lambda+(gi-1)*delta_lambda)
            sin_u_longitude(i, j) =                                     &
     &                      sin (Base_lambda+(gi-.5)*delta_lambda)
            cos_u_longitude(i, j) =                                     &
     &                      cos (Base_lambda+(gi-.5)*delta_lambda)
          End Do
        End Do

      End If !model_domain  ==  mt_global .or. .not.L_trivial_trigs

!!! Allocate arrays for dyn_coriolis_mod module
      If (.not.allocated(f1_at_v)) Then
        Allocate (f1_at_v (row_length, 0:n_rows+1))
      End If
      If (.not.allocated(f2_at_u)) Then
        Allocate (f2_at_u (0:row_length+1, rows))
      End If
      If (.not.allocated(f3_at_u)) Then
      Allocate (f3_at_u (1-offx:row_length+offx, 1-offy:rows+offy))
      End If
      If (.not.allocated(f3_at_v)) Then
      Allocate (f3_at_v (1-offx:row_length+offx, 1-offy:n_rows+offy))
      End If

      If (model_domain  ==  mt_global) Then
! global model
        Do j = 1, n_rows
          Do i = 1, row_length
            f1_at_v(i,j) = 0.
            f3_at_v(i,j) = two_omega * sin_v_latitude(i,j)
          End Do
        End Do
        Do j = 1, rows
          Do i = 1, row_length
            f2_at_u(i,j) = two_omega * cos_theta_latitude(i,j)
            f3_at_u(i,j) = two_omega * sin_theta_latitude(i,j)
          End Do
        End Do
        Do j = 1, rows
          Do i = 1, row_length
            gi = datastart(1) + i - 1
            true_longitude(i,j) = (gi-1)*delta_lambda
          End Do
        End Do
! set polar values of true longitude to be all the same

        If (at_extremity(PNorth)) then
          Do i = 1, row_length
            true_longitude(i,rows) = 0.
          End Do
        End If

        If (at_extremity(PSouth)) then
          Do i = 1, row_length
            true_longitude(i,1) = 0.
          End Do
        End If

! get parameters for common latlonmax within mppac
! these need to be in degrees
      lat_n = Base_phi_wk + (datastart(2)+rows-2) * Delta_phi_wk
      lat_s = Base_phi_wk

      long_w       = Base_lambda_wk
      long_w_model = long_w
      long_e       = Base_lambda_wk +                                   &
     &              (datastart(1)+row_length-2) * Delta_lambda_wk
      long_e_model = long_e

      if(long_w >  180.0) long_w = long_w-360.0
      if(long_e >  180.0) long_e = long_e-360.0

      Else
! limited area model
        temp1 = cos( lat_rot_NP)
        temp2 = sin( lat_rot_NP)

      If( L_trivial_trigs .or. f_plane  >   -89.0) Then
        f_plane_rad = f_plane * Pi_over_180
        ff_plane_rad = ff_plane * Pi_over_180
        f1_temp = - two_omega * temp1 * sin(ff_plane_rad)
        f2_temp = two_omega * (cos(f_plane_rad) * temp2                 &
     &               - sin(f_plane_rad) * temp1 * cos(ff_plane_rad))
        f3_temp = two_omega * ( sin(f_plane_rad) * temp2                &
     &               + cos(f_plane_rad) * temp1 * cos(ff_plane_rad))
        If (L_vert_Coriolis) Then
          f1_temp = 0.0
          f2_temp = 0.0
          f3_temp = two_omega * ( sin(f_plane_rad) * temp2              &
     &               + cos(f_plane_rad) * temp1 * cos(ff_plane_rad))
        End If  !  L_vert_Coriolis

        Do j = 1, n_rows
          Do i = 1, row_length
            f1_at_v(i,j) = f1_temp
            f3_at_v(i,j) = f3_temp
          End Do
        End Do
        Do j = 1, rows
          Do i = 1, row_length
            f2_at_u(i,j) = f2_temp
            f3_at_u(i,j) = f3_temp
          End Do
        End Do

      Else        !  L_trivial_trigs = .false.
        Do j = 1, n_rows
          Do i = 1, row_length
            f1_at_v(i,j) = - two_omega * temp1 *                        &
     &                       sin_theta_longitude(i,j)
            f3_at_v(i,j) = two_omega * ( sin_v_latitude(i,j) * temp2    &
     &                                  +cos_v_latitude(i,j) * temp1    &
     &                                  *cos_theta_longitude(i,j) )
          End Do
        End Do
        Do j = 1, rows
          Do i = 1, row_length
            f2_at_u(i,j) = two_omega * ( cos_theta_latitude(i,j) * temp2&
     &                                  -sin_theta_latitude(i,j) * temp1&
     &                                  *cos_u_longitude(i,j) )
            f3_at_u(i,j) = two_omega * ( sin_theta_latitude(i,j) * temp2&
     &                                  +cos_theta_latitude(i,j) * temp1&
     &                                  *cos_u_longitude(i,j) )
          End Do
        End Do
      End If   !  L_trivial_trigs

! calculate true longitude in radians
        Do j = 1, rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
            gi = datastart(1) + i - 1
            longitude_rot(i,j) = (base_lambda + (gi-1) * delta_lambda)  &
     &                         / Pi_over_180
            latitude_rot(i,j) = (base_phi + (gj-1) * delta_phi)         &
     &                         / Pi_over_180
          End Do
        End Do

! DEPENDS ON: eqtoll
        Call eqtoll(latitude_rot, longitude_rot                         &
     &,             true_latitude, true_longitude                       &
     &,             lat_rot_NP_in, long_rot_NP_in, rows*row_length)

        Do j = 1, rows
          Do i = 1, row_length
            true_longitude(i,j) = true_longitude(i,j) * Pi_over_180
            true_latitude(i,j) = true_latitude(i,j) * Pi_over_180
          End Do
        End Do

! get parameters for common latlonmax within mppac
! these need to be in degrees
      lat_n = latitude_rot(1,rows)
      lat_s = latitude_rot(1,1)

      long_w       = longitude_rot(1,1)
      long_w_model = long_w
      long_e       = longitude_rot(row_length,1)
      long_e_model = long_e

      if(long_w >  180.0) long_w = long_w-360.0
      if(long_e >  180.0) long_e = long_e-360.0
      if(long_w_model >= 360.0)long_w_model = long_w_model-360.
      if(long_e_model >= 360.0)long_e_model = long_e_model-360.

! calculate lat/longitude for points on equatorial grid for B grid
        Do j = 1, n_rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
            gi = datastart(1) + i - 1
            longitude_rot(i,j) = (base_lambda + (gi-.5) * delta_lambda) &
     &                         / Pi_over_180
            latitude_rot(i,j) = (base_phi + (gj-.5) * delta_phi)        &
     &                         / Pi_over_180
          End Do
        End Do

! DEPENDS ON: eqtoll
        Call eqtoll(latitude_rot, longitude_rot                         &
     &,             true_latitude_b, true_longitude_b                   &
     &,             lat_rot_NP_in, long_rot_NP_in                       &
     &,             n_rows*row_length)

! Calculate rotation coefficients for wind

!!! Allocate arrays for rot_coeff_mod module
      If (.not.allocated(rot_coeff1)) Then
        Allocate (rot_coeff1 ( row_length, n_rows ))
      End If
      If (.not.allocated(rot_coeff2)) Then
        Allocate (rot_coeff2 ( row_length, n_rows ))
      End If

! DEPENDS ON: w_coeff
        Call w_coeff(rot_coeff1, rot_coeff2,                            &
     &               true_longitude_b, longitude_rot,                   &
     &              lat_rot_NP_in, long_rot_NP_in,                      &
     &               n_rows*row_length )

      End If ! model_domain

! Call swap_bounds for those trig fields which require halos

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (FV_cos_theta_latitude, row_length, rows, 1,    &
     &                   Offx, Offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(FV_cos_theta_latitude,row_length,rows,   &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (FV_sec_theta_latitude, row_length, rows, 1,    &
     &                   Offx, Offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(FV_sec_theta_latitude,row_length,rows,   &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (cos_theta_latitude, row_length, rows, 1,       &
     &                   Offx, Offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(cos_theta_latitude,row_length,rows,      &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (sec_theta_latitude, row_length, rows, 1,       &
     &                   Offx, Offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(sec_theta_latitude,row_length,rows,      &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (cos_v_latitude, row_length, n_rows, 1,         &
     &                   Offx, Offy, fld_type_v, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(cos_v_latitude,row_length,n_rows,        &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (sec_v_latitude, row_length, n_rows, 1,         &
     &                   Offx, Offy, fld_type_v, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(sec_v_latitude,row_length,n_rows,        &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (f2_at_u, row_length, rows, 1,                  &
     &                   1, 0, fld_type_u, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(f2_at_u,row_length,rows,1,1,0)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (f1_at_v, row_length, n_rows, 1,                &
     &                   0, 1, fld_type_v, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(f1_at_v,row_length,n_rows,1,0,1)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (f3_at_u, row_length, rows, 1,                  &
     &                   Offx, Offy, fld_type_u, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(f3_at_u,row_length,rows,1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (f3_at_v, row_length, n_rows, 1,                &
     &                   Offx, Offy, fld_type_v, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(f3_at_v,row_length,n_rows,1,Offx,Offy)

! ----------------------------------------------------------------------
! 1. Set polar filtering and diffusion
! ----------------------------------------------------------------------
!!! Allocate arrays for diff_coeff_mod module
      If (.not.allocated(diff_coeff_u)) Then
        Allocate (diff_coeff_u ( 1-Offx : row_length+Offx,              &
     &                           1-Offy : rows+Offy ))
      End If
      If (.not.allocated(diff_coeff_v)) Then
        Allocate (diff_coeff_v ( 1-Offx : row_length+Offx,              &
     &                           1-Offy : n_rows+Offy ))
      End If

      If ( L_pofil_new ) Then

! DEPENDS ON: setdiff
        Call Setdiff(                                                   &
! Input size and control variables
     &      global_rows, row_length, rows, n_rows, model_levels,        &
     &      model_domain, at_extremity, datastart,                      &
     &      offx, offy, mype, nproc_y,                                  &
     &      max_model_levels, max_121_rows, max_sweeps,                 &
! other constants
     &      delta_lambda, delta_phi,                                    &
     &      pi, Pi_over_180, polar_cap, scale_ratio,                    &
     &      ref_lat_deg, diff_coeff_ref,                                &
     &      cos_theta_latitude, sin_theta_latitude,                     &
     &      cos_v_latitude, sin_v_latitude,                             &
! Output data
     &       global_u_filter, global_v_filter,                          &
     &       u_sweeps, v_sweeps,                                        &
     &       u_begin, u_end, v_begin, v_end,                            &
     &       diff_coeff_u, diff_coeff_v, diff_coeff_phi,                &
     &       diffusion_coefficient_thermo, diffusion_order_thermo,      &
     &       diffusion_coefficient_wind, diffusion_order_wind,          &
     &       diff_order_thermo, diff_order_wind,                        &
     &       diff_timescale_thermo, diff_timescale_wind,                &
     &       L_sponge, sponge_power, sponge_ew, sponge_ns,              &
     &       sponge_wts_ew, sponge_wts_ns,                              &
     &       L_diffusion, L_cdiffusion, L_filter, L_diff_auto,          &
     &       L_pfcomb, L_pfexner, L_pftheta, L_pfuv, L_pfw, L_pfincs,   &
     &       L_diff_thermo, L_diff_wind, L_diff_w, L_diff_incs)

      Else ! original combi/diffusion setting

! DEPENDS ON: setdiff_old
        Call Setdiff_old(                                               &
! Input size and control variables
     &      global_rows, row_length, rows, n_rows,                      &
     &      model_domain, at_extremity, datastart,                      &
     &      offx, offy, mype, nproc_y, max_121_rows,                    &
! other constants
     &      delta_lambda, delta_phi,                                    &
     &      pi, Pi_over_180, polar_cap, scale_ratio, diff_coeff_ref,    &
     &      cos_theta_latitude, sin_theta_latitude,                     &
     &      cos_v_latitude, sin_v_latitude,                             &
! Output data
     &       global_u_filter, global_v_filter,                          &
     &       u_sweeps, v_sweeps,                                        &
     &       u_begin, u_end, v_begin, v_end,                            &
     &       diff_coeff_u, diff_coeff_v,                                &
     &       diff_coeff_thermo, diff_coeff_wind,                        &
     &       diff_order_thermo, diff_order_wind,                        &
     &       diff_timescale_thermo, diff_timescale_wind,                &
     &       L_filter, L_pfcomb, L_pftheta, L_pfuv, L_pfw, L_pfincs,    &
     &       L_diff_thermo, L_diff_wind, L_diff_w, L_diff_incs)

      End If ! L_pofil_new

! ----------------------------------------------------------------------
! 1.5 Set filtering control variable
! ----------------------------------------------------------------------

      If ( L_polar_filter .or. L_pfcomb ) Then
        L_filter = .true.
        if ( model_domain  /=  mt_global) then
          write(6,*)' ******   Polar filter NOT ALLOWED in LAMs ****** '
          write(6,*)' ******       RESET UMUI or NAMELIST       ****** '
          ICODE    = 1
          CMESSAGE ='SETCONA:Polar filter NOT ALLOWED in LAMs'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        else if ( L_pofil_new ) then
          write(6,*)' ******  Newest filtering code being used ****** '
        end if ! model_domain  /=  mt_global
      Else
        write(6,*) ' ******   Polar filter is OFF   ****** '
      End If ! L_polar_filter .or. L_pfcomb
      If ( L_diff_thermo .or. L_diff_wind .or. L_diff_w ) Then
        if ( L_pofil_new ) then
          write(6,*) ' ******  Newest diffusion code being used ****** '
        endif !  L_pofil_new
        L_filter = .true.
      End If ! L_diff_thermo .or. L_diff_wind .or. L_diff_w
      If ( L_diff_incs .or. L_pfincs) Then
        L_filter_incs = .true.
      End If ! L_diff_incs .or. L_pfincs
      If ( L_pfexner  ) Then
        write(6,*) ' ** Polar filtering of Exner pressure is active ** '
      End If ! L_pfexner
      If ( L_diff_exner  ) Then
        write(6,*) ' ** Diffusion of Exner pressure is active ** '
      End If ! L_diff_exner

      If( L_tardiff_q ) Then
        write(6,*) ' '
        write(6,*) '  ***   Targeted diffusion of q is ACTIVE  *** '
        write(6, 913) w_conv_limit, tardiffq_factor
        write(6, 914) tardiffq_test
        write(6, 915) tardiffq_start, tardiffq_end
        write(6, 916) tar_horizontal
      Else
        write(6,*) ' '
        write(6,*) '  ***   Targeted diffusion of q is NOT ACTIVE  *** '
      End If !   L_tardiff_q

!L ----------------------------------------------------------
!  Print information to tell user what is active
!L ----------------------------------------------------------
!  vertical diffusion coefficients
      If ( L_vdiff_uv ) Then
        vdiffuv_factor = 0.25 *                                         &
     &                     (1.0 - EXP( -1.0 / vdiffuv_timescale) )
       Write ( Unit = 6, fmt=*) 'A factor delta_z**2/timestep is '      &
     & ,' allowed for in the application of the vertical diffusion.'
       Write ( Unit = 6, fmt=*) 'Vertical diffusion of order 1 with'    &
     &        , ' vertical diffusion coefficient = ', vdiffuv_factor
       Write ( Unit = 6, fmt=*) 'Vertical diffusion timescale is '      &
     &                           , vdiffuv_timescale,' steps '
      End If ! L_vdiff_uv

      if ( top_filt_start > model_levels ) then
        Write ( Unit = 6, fmt=*) 'No additional upper-level diffusion'
      else  ! top_filt_start <= model_levels
        active_levels = top_filt_end - top_filt_start + 1
        if ( active_levels > max_updiff_levels) then
          top_filt_start  = top_filt_end - max_updiff_levels + 1
          write(6,*) ' Max of uppermost ', max_updiff_levels            &
     &          , ' levels allowed for upper-level'                     &
     &          , ' diffusion - start level reset to ', top_filt_start
          active_levels = max_updiff_levels
        end if ! active_levels > max_updiff_levels
        Write ( Unit = 6, fmt=*) 'Extra upper-level diffusion applied'
        scale = 1.0
        if ( L_upper_ramp ) then
          scale = up_diff_scale
          write(6, 911) top_diff, top_filt_end
          write(6, 912) scale, top_filt_start
        else  !
          write(6, 910) top_diff, top_filt_start, top_filt_end
        end if  !  L_upper_ramp
        kk = active_levels
        up_diff(active_levels) = top_diff
        do k = top_filt_end - 1, top_filt_start, - 1
          kk = kk - 1
          up_diff(kk) = scale * up_diff(kk + 1)
        enddo !  k = top_filt_end - 1, top_filt_start, - 1
        Do k = 1, active_levels
          kk = k + top_filt_start - 1
          write(6,*)'Level', kk,' Diffusion factor = ',up_diff(k)
        End do !  k = 1, active_levels
      end if ! top_filt_start > model_levels

      if (L_adjust_theta) then
        if( adjust_theta_start < 2 ) then
          adjust_theta_start = 2
          Write ( Unit = 6, fmt=*) 'Start level for convective '        &
     & , ' adjustment reset to ', adjust_theta_start
        end if ! adjust_theta_start < 2
        Write ( Unit = 6, fmt=*) 'Convective adjustment applied when'   &
     &    ,' lapse rate is negative above level ',  adjust_theta_start  &
     &    ,' up to level ', adjust_theta_end
      end if  ! L_adjust_theta

! ----------------------------------------------------------------------
! Section 1.5a . Print FORMATTING
! ----------------------------------------------------------------------

 910  FORMAT(' Diffusion coefficient = ',F5.2                           &
     &      ,' from level ', I3,' to level ', I3)
 911  FORMAT(' Diffusion coefficient = ',F5.2,' at level ', I3)
 912  FORMAT(' decreasing by ', F6.3, ' at each level down to ', I3)
 913  FORMAT('   Vertical velocity test value is  ',F5.2,' m/s and ',   &
     &                                 'diffusion factor = ',F6.3)
 914  FORMAT('   Start testing at level  ',I3)
 915  FORMAT('   Apply from level  ',I3,' to level ',I3)
 916  FORMAT('   Slope test applied up to level ',I4)


! ----------------------------------------------------------------------
! 2. Semi-Lagrangian scheme options.
! ----------------------------------------------------------------------


! j holds number of points required for interpolation scheme
      j = 2
      Do i = 1, 3
        If (high_order_scheme(i)  ==  quinticLagrange) Then
          j = 3
          If (halo_j  <   5) Then
            write(6,*)' error halo_j too small ',halo_j
          End If
        Else
          If (halo_j  <   4) Then
            write(6,*)' error halo_j too small ',halo_j
          End If
        End If
      End Do

! Calculate max cfl possible in LAM
! Y-direction also applies in global model though variable is not
! used in code.

      LAM_max_cfl(1) = halo_i - j
      LAM_max_cfl(2) = halo_j - j

! ----------------------------------------------------------------------
! Set up data.
! ----------------------------------------------------------------------

! Calculate land_index
          land_points = 0
          Do j =1, rows
            Do i = 1, row_length
              If (land_sea_mask(i,j)) then
                land_points = land_points + 1
                land_index(land_points) = (j-1)*row_length + i
              End If
            End Do
          End Do

! set-up code for hydrology
          soil_points = 0
          land_ice_points = 0
          Do i= 1, land_points
! Test on soil moisture concentration at saturation
            If (smvcst(i) >   0.0) Then       ! Soil points
              soil_points = soil_points+1
              soil_index(soil_points)=i
            Else If (smvcst(i)  ==  0.0) Then   ! Land-ice points
              land_ice_points = land_ice_points+1
              land_ice_index(land_ice_points)=i
            End If
          End Do


! calculate r_at_u points on rho levels and
! r_at_v points on rho levels.

!!! Allocate r_at_u/v arrays for level_heights_mod module
      If (.not.allocated(r_at_u)) Then
        Allocate (r_at_u (1-halo_i:row_length+halo_i,                   &
     &                    1-halo_j:rows+halo_j, model_levels))
      End If
      If (.not.allocated(r_at_v)) Then
        Allocate (r_at_v (1-halo_i:row_length+halo_i,                   &
     &                    1-halo_j:n_rows+halo_j, model_levels))
      End If

        Do k = 1, model_levels
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i - 1
              r_at_u (i,j,k) = .5 * (r_rho_levels(i,j,k) +              &
     &                                 r_rho_levels(i+1,j,k) )
            End Do
          End Do
        End Do

        Do k = 1, model_levels
          Do j = 1-halo_j, n_rows+halo_j - 1
            Do i = 1-halo_i, row_length+halo_i
              r_at_v (i,j,k) = .5 * (r_rho_levels(i,j,k) +              &
     &                               r_rho_levels(i,j+1,k) )
            End Do
          End Do
        End Do

! call swap_bounds to set extra points
! DEPENDS ON: swap_bounds
       call Swap_Bounds(r_at_u, row_length, rows, model_levels,         &
     &                  halo_i, halo_j, fld_type_u, .false.)

! DEPENDS ON: swap_bounds
       call Swap_Bounds(r_at_v, row_length, n_rows, model_levels,       &
     &                  halo_i, halo_j, fld_type_v, .false.)

      If ( L_diag_L2norms ) then
        Write ( Unit = 6, fmt=*) 'Printing of norms during timestep'    &
     &, ' activated'
      endIf ! L_diag_L2norms
      If ( L_diag_L2helm ) then
        Write ( Unit = 6, fmt=*) 'Printing of coefficient norms  '      &
     &, ' in solver'
      endIf ! L_diag_L2helm
!  Initialise diagnostic printing items
! Sampling interval must not be larger than print interval
      if ( diag_interval > print_step) then
        diag_interval = print_step
      end if !  diag_interval > print_step
! Set array sizes needed for chosen diagnostic printing
! These are needed to do sums and max/mins across processors
      rpemax = 0
      rpemin = 0
      ipesum = 0
      rpesum = 0
      if (L_print_theta1 ) then
        rpemin = rpemin + 1
        ipesum = ipesum + 1
        rpesum = rpesum + 3
      end if !  L_print_theta1
      if (L_print_lapse ) then
        rpemin = rpemin + (model_levels - 1)
        ipesum = ipesum + 2 * (model_levels - 1)
        rpesum = rpesum + 2 * (model_levels - 1)
      end if !  L_print_lapse
      if (L_print_div) then
        rpemax = rpemax + model_levels
        rpemin = rpemin + model_levels
        ipesum = ipesum + 2 * model_levels
        rpesum = rpesum + 4 * model_levels
      end if !  L_print_div
      if (L_print_w .or. L_print_wmax) then
        rpemax = rpemax + model_levels - 1
        ipesum = ipesum + 2 * (model_levels - 1)
        rpesum = rpesum + 2 * (model_levels - 1)
      end if !  L_print_w .or. L_print_wmax
      if (L_print_shear) then
        rpemax = rpemax + model_levels - 1
        ipesum = ipesum + model_levels - 1
        rpesum = rpesum + 2 * (model_levels - 1)
      end if !  L_print_shear
      if (L_print_max_wind) then
        rpemax = rpemax + model_levels
        ipesum = ipesum + model_levels
        rpesum = rpesum + 2 * model_levels
      end if !  L_print_max_wind
      if ( L_diag_wind ) then
        rpesum = rpesum + 2 * model_levels
      end if  ! L_diag_wind
        min_theta1_run = 1000.0
        time_theta1_min = 0
      Do k = 1, model_levels
         max_w_run(k) =  0.0
         time_w_max(k) = 0
         max_div_run(k) =  0.0
         time_div_max(k) = 0
         min_div_run(k) =  0.0
         time_div_min(k) = 0
        min_lapse_run(k) =  1.0e6
         time_lapse_min(k) = 0
        max_shear_run(k) = 0.0
        time_max_shear(k) = 0
        max_wind_run(k) = 0.0
        time_max_wind(k) = 0
        max_KE_run(k) = 0.0
        min_KE_run(k) = 1.0e30
        time_KE_max(k) = 0
        time_KE_min(k) = 0
        time_noise_max(k) = 0
        max_noise_run(k) = 0.0
      End Do  ! k = 1, model_levels
      max_KE_run(model_levels + 1) = 0.0
      min_KE_run(model_levels + 1) = 1.0e30
      time_KE_max(model_levels + 1) = 0
      time_KE_min(model_levels + 1) = 0

!L ----------------------------------------------------------
!L 3. Set up cloud type boundaries for low/medium/high cloud.
!L ----------------------------------------------------------
       IF (NUM_CLOUD_TYPES  >   3) THEN
         ICODE    = 1
         CMESSAGE = 'SETCONA: Parameter NUM_CLOUD_TYPES exceeds 3'
         WRITE(6,*) 'NUM_CLOUD_TYPES=',NUM_CLOUD_TYPES
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ICODE,Cmessage)

       END IF

!  Diagnostics of cloud amount for low/medium/high cloud are calculated
!  by finding the max cloud cover over a range of model levels. The
!  ranges are determined by heights [ 1->2:2->3:3->4 = L:M:H ] held in
!  array h_split by comparison with heights of model levels (above
!  surface) in r_ref_theta.

      DO KK = 1, NUM_CLOUD_TYPES + 1
        LEVEL = 1
!
!       r_ref_theta turns out to be A_theta(k) - Earth_radius at every
!       model level because z_top_of_model = z_rho_top is chosen here.
!
        DO WHILE ((r_ref_theta(LEVEL)  <=  h_split(KK)) .AND.           &
     &                                       (LEVEL  <=  model_levels))
          LEVEL = LEVEL + 1

        END DO

        IF (LEVEL  >   model_levels) THEN
          ICODE    = 1
          CMESSAGE ='SETCONA:Error in locating levels for cloud layers'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)

        END IF
        CLOUD_BOUND(KK) = LEVEL

      END DO

      LOW_BOT_LEVEL  = CLOUD_BOUND(1)
      LOW_TOP_LEVEL  = CLOUD_BOUND(2) - 1
      MED_BOT_LEVEL  = CLOUD_BOUND(2)
      MED_TOP_LEVEL  = CLOUD_BOUND(3) - 1
      HIGH_BOT_LEVEL = CLOUD_BOUND(3)
      HIGH_TOP_LEVEL = CLOUD_BOUND(4) - 1

      IF (LOW_TOP_LEVEL  >   CLOUD_LEVELS) THEN
        ICODE    = 1
        CMESSAGE = 'SETCONA: No of cloud levels less than Top of Low'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF

      IF (MED_TOP_LEVEL >  CLOUD_LEVELS) THEN
        ICODE    = 1
        CMESSAGE = 'SETCONA:  No of cloud levels less than Top of Med'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF

      IF (HIGH_TOP_LEVEL  >   CLOUD_LEVELS) THEN
        ICODE    = 1
        CMESSAGE = 'SETCONA: No of cloud levels less than Top of High'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF

!L ----------------------------------------------------------
!L 4. Set up locally-derived parameters.
!L ----------------------------------------------------------

!     The tropopause diagnosed for radiative purposes
!     divides theta-levels considered to lie in the stratosphere
!     from those considered to lie in the troposphere: the
!     tropopause is therefore taken as a rho-level. This level
!     is constrained to lie between heights of z_min_trop
!     and z_max_trop. The level is used in the setting of the
!     background aerosol climatology and of ozone profiles,
!     subject to the choice of appropriate options; additionally
!     it is used in the calculation of diagnostics defined at
!     the tropopause.
!
!     Start at the second rho-level because the first is omitted
!     from the grid seen in the physics.
      min_trop_level=2
      Do ; If ( (r_ref_rho(min_trop_level) >= z_min_trop) .OR.          &
     &          (min_trop_level == model_levels) ) Exit
        min_trop_level = min_trop_level+1
      End Do
!
      max_trop_level=min_trop_level
      Do ; If ( (r_ref_rho(max_trop_level) > z_max_trop) .OR.           &
     &          (max_trop_level == model_levels) ) Exit
        max_trop_level = max_trop_level+1
      End Do
      max_trop_level = max_trop_level-1
!
! set up super tracer array size

          super_array_size=0

          if (l_CO2_interactive) super_array_size=super_array_size+1
          If (l_Soot) super_array_size=super_array_size+3
          If (l_Biomass) super_array_size=super_array_size+3
          If (l_Sulpc_so2)super_array_size=super_array_size+4
          if (L_sulpc_nh3)super_array_size=super_array_size+1
          if (L_sulpc_dms)super_array_size=super_array_size+1
          IF (L_DUST) super_array_size=super_array_size+6
          if (L_ocff) super_array_size=super_array_size+3
          IF (L_Murk_advect)super_array_size=super_array_size+1
          if (L_USE_CARIOLLE)super_array_size=super_array_size+1
          if (tr_levels==model_levels.and.tr_vars>0)                    &
     &            super_array_size=super_array_size+tr_vars
          if (tr_levels==model_levels.and.tr_ukca>0)                    &
     &            super_array_size=super_array_size+tr_ukca


          if(mype == 0)write(6,*)'super_array_size =  ',                &
     &                          super_array_size

! set up moisture array size

!        always do q,qcl,qcf
          moisture_array_size=3
          if (l_pc2) moisture_array_size=moisture_array_size+4
          if (l_mcr_qrain)moisture_array_size=moisture_array_size+1
          if (l_mcr_qcf2 )moisture_array_size=moisture_array_size+1
          if (l_mcr_qgraup)moisture_array_size=moisture_array_size+1
          if(mype == 0)write(6,*)'moisture_array_size =  ',             &
     &                          moisture_array_size

!      end of set up code

! ----------------------------------------------------------------
! 5. Calculate the coeffs for the interpolation of the radiation
!    quantities if spatial degradation of the radiation code is on.
!    Calculate the mask which ensures a chequer-board pattern over
!    the whole domain (and not just for a single PE).
! -----------------------------------------------------------------

! Create a land mask with a halo.
      IF (L_rad_deg) THEN

        landmask(1-offx:0, 1-offy: rows+offy) = .FALSE.
        landmask(row_length+1:row_length+offx, 1-offy:rows+offy)=.FALSE.
        landmask(1:row_length, 1-offy: 0) = .FALSE.
        landmask(1:row_length, rows+1: rows+offy) = .FALSE.

        DO j=1,ROWS
          DO i=1,ROW_LENGTH
            landmask(i,j)=Land_sea_mask(i,j)
          END DO
        END DO

! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(landmask,row_length,rows,1,                    &
     &                   offx,offy,fld_type_p,.FALSE.)

        first_row=1
        last_row=rows
        IF (model_domain  ==  mt_global) THEN
          IF (at_extremity(PNorth)) THEN
            last_row=rows-1
          END IF
          IF (at_extremity(PSouth)) THEN
            first_row=2
          END IF
        END IF

!!! Allocate arrays for rad_mask_trop_mod module
      If (.not.allocated(es_space_interp)) Then
        Allocate (es_space_interp(4, row_length, rows))
      End If
      If (.not.allocated(rad_mask)) Then
        Allocate (rad_mask(row_length, rows))
      End If

! DEPENDS ON: coeffs_degrade
        CALL COEFFS_DEGRADE(ES_SPACE_INTERP, landmask,                  &
     &                      first_row, last_row, model_domain,          &
     &                      mt_lam, at_extremity(PNorth),               &
     &                      at_extremity(PSouth), at_extremity(PWest),  &
     &                      at_extremity(PEast),                        &
     &                      row_length*rows,row_length, rows,offx, offy)
! DEPENDS ON: rad_degrade_mask
        CALL RAD_DEGRADE_MASK(RAD_MASK, datastart,                      &
     &                      Ndim_max, first_row, last_row, offx,        &
     &                      row_length, rows)

      ENDIF  !  If L_rad_deg loop

        If ( thmono_height >= 0.5 * r_ref_theta(1) ) then
          If ( thmono_height < r_ref_theta(1) ) then
            k = 1
          else if ( thmono_height > r_ref_theta(model_levels) ) then
            k = model_levels
          else
            k = 1
            do  !  cycle to find nearest level to  thmono_height
              k = k + 1
              if ( thmono_height < r_ref_theta(k) ) then
                if ( r_ref_theta(k) - thmono_height >                   &
     &               thmono_height - r_ref_theta(k-1) ) k = k - 1
                EXIT
              end if ! thmono_height < r_ref_theta(k)
              CYCLE
            end do !  cycle to find nearest level to  thmono_height
          End If ! thmono_height < r_ref_theta(1)
          thmono_levels = k
          Write ( Unit = 6, fmt=*) '3D Monotone limiter will be applied'&
     &,                            ' to advection of theta'
          Write ( Unit = 6, fmt=*) 'thmono_height was set to '          &
     &,                             thmono_height,' in the UMUI'
          Write ( Unit = 6, fmt=*) 'Limiter will be applied up to level'&
     &,   thmono_levels,', the nearest level to',thmono_height,'metres'
        Else   ! thmono_height < 0.5 * r_ref_theta(1)
          thmono_levels = 0
        End If ! thmono_height >= 0.5 * r_ref_theta(1)
!L ----------------------------------------------------------
!L 5. Set up dynamical core parameters
!L ----------------------------------------------------------
!    Set frictional_timescale = 0 everywhere since this value
!    is passed into PE_HELMHOLTZ from ATMSTEP
!    The value used in simple friction is calculated internally
!           ( stored in friction_level)
!     as is SuHe_level_weight (temp1) in the temperature relaxation

      Do k = 1, model_levels
        frictional_timescale(k) = 0.0
      End Do

      If( problem_number  ==  dynamical_core) Then

! DEPENDS ON: idl_set_suhe_params
        Call IDL_Set_SuHe_params                                        &
     &                      (row_length, rows, model_levels             &
     &,                     kappa, Cp, p_zero, Earth_radius, pi, R, g   &
     &,                     offx, offy, halo_i, halo_j                  &
     &,                     mype, nproc, at_extremity                   &
     &,                     datastart, gc_all_proc_group                &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     SuHe_pole_equ_deltaT, SuHe_static_stab      &
     &,                     p, p_theta_levels                           &
     &,                     base_frictional_timescale                   &
     &,                     friction_level, SuHe_sigma_cutoff           &
     &,                     SuHe_level_weight )

      else if( problem_number  ==  idealised_problem) Then
!     Clear arrays
        Do k = 1, model_levels
          SuHe_level_weight(k) = 0.0
          friction_level(k) = 0.0
        End Do    ! k = 1, model_levels

      end if    ! problem_number  ==  dynamical_core
! ----------------------------------------------------------------------
! 6. GCR diagnostics initialisation.
! ----------------------------------------------------------------------
      if (GCR_diagnostics == 3 )then
        L_error = .false.
        Do i = 1, 3
          if(GCR_its_avg_step(i) == 0 )L_error = .true.
        EndDo  !  i = 1, 3
        if ( L_error ) then
          Write ( Unit = 6, fmt=*) 'WARNING GCR iteration counting at'  &
     &,         ' timestep 0 or interval of 0 NOT PERMITTED '
! Following values will output iteration info at 6, 12 and
!  1440 timesteps (30 days for a 30 minute timestep)
          GCR_its_avg_step(1) = 6
          GCR_its_avg_step(2) = 12
          GCR_its_avg_step(3) = 1440
          write(6,*) ' Iteration count diagnostics reset for timesteps '&
     &          , GCR_its_avg_step(1), GCR_its_avg_step(2)              &
     &                               , GCR_its_avg_step(3)
          write(6,*)  '  and at intervals of ', GCR_its_avg_step(3),    &
     &            ' timesteps thereafter.'
          write(6,*)  ' Change in UMUI if you want different values '
        else ! L_error .false.
          write(6,*)  ' '
          write(6,*)  'Iteration count diagnostics at timesteps '       &
     &          , GCR_its_avg_step(1), GCR_its_avg_step(2)              &
     &          , GCR_its_avg_step(3)
          write(6,*)  ' and at intervals of ', GCR_its_avg_step(3),     &
     &            ' timesteps thereafter.'
        endif ! L_error
        GCR_its_switch =  1
        GCR_sum_its = 0
        GCR_min_its = 1000
        GCR_max_its = 0
        GCR_max_time = 0
        GCR_min_time = 0
      endif ! GCR_diagnostics == 3

! ----------------------------------------------------------------------
! 7. Error trapping for advection choices
! ----------------------------------------------------------------------

      do i=1,3
        if (.not. L_mono(i) .and. .not. L_high(i)) then
          if(i == 1)write ( Unit = 6, fmt=*)                            &
     &      'WARNING Advection choices incompatible for theta'
          if(i == 2)write ( Unit = 6, fmt=*)                            &
     &      'WARNING Advection choices incompatible for moisture'
          if(i == 3)write ( Unit = 6, fmt=*)                            &
     &      'WARNING Advection choices incompatible for winds'
          Write ( Unit = 6, fmt=*)                                      &
     &      'Both L_high and L_mono switches set to false '

        endif
      end do

! ----------------------------------------------------------------------
! Orography slope angle and aspect for the radiation code.
! ----------------------------------------------------------------------

      If (L_orog) Then

         If (model_domain == mt_global) Then
            bear_rot_NP = 0.0
         Else
! DEPENDS ON: pole_bearing
            Call pole_bearing(row_length, rows,                         &
     &            lat_rot_NP, long_rot_NP, true_longitude,              &
     &            f3_at_u(1:row_length,1:rows), two_Omega,              &
     &            bear_rot_NP)
         Endif

         If (L_use_grad_corr) Then
! DEPENDS ON: aspang_ancil
            Call aspang_ancil(row_length, rows, land_points,            &
     &            land_sea_mask, grad_x, grad_y, bear_rot_NP)
         Else
! DEPENDS ON: aspang
            Call aspang(row_length,rows, Delta_lambda, Delta_phi,       &
     &            Earth_radius, Orog_halos(0:row_length+1,0:rows+1),    &
     &            cos_theta_latitude(1:row_length,1:rows),              &
     &            bear_rot_NP)
         Endif

         If (model_domain == mt_global) Then
            If (at_extremity(PNorth)) slope_angle(:,rows) = 0.0
            If (at_extremity(PSouth)) slope_angle(:,1) = 0.0
         Endif

      Endif


! ----------------------------------------------------------------------
! Boundary layer solver options
! ----------------------------------------------------------------------
      If ( L_us_blsol ) Then
        If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
          write(6,*)
          write(6,*)  '*** New stable and non-oscillatory boundary-layer&
     & solver is ACTIVE ***'
          write(6,*)  '     It will run with Puns=',Puns,'Pstb=',Pstb
        End If 
      Else
        If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
          write(6,*)
          write(6,*)'*** Standard boundary-layer solver is ACTIVE ***'
        End If 
      End If
  
! ----------------------------------------------------------------------
! read in volcanic forcing, VOLCTS
! ----------------------------------------------------------------------
      if (L_VOLCTS) then
      write(6,*) 'SETCONA: Reading in volcanic forcing from file: ',    &
     &            FILE_VOLCTS
      OPEN(UNIT=56,file=FILE_VOLCTS,status='old') 
      do i=1850,2300 
        do j=1,12 
          read(56,*) yvolc,mvolc,volcts(:,j,i) 
        enddo 
      enddo 
      CLOSE(56)
      endif

! ----------------------------------------------------------------------
! read in varying Solar constant, SCVARY
! ----------------------------------------------------------------------
!
! Year-to-year variation of the solar "constant" - from 
! Solanki & Krivova (2003). + solar cycle into 21stC  
!
      if (L_SCVARY) then
      write(6,*) 'SETCONA: Reading in solar forcing from file: ',    &
     &            FILE_SCVARY
      OPEN(UNIT=56,file=FILE_SCVARY,status='old')
      do i=1,601
        read(56,*) ysol,scvary(i)
      enddo
      CLOSE(56)
      endif

  
      RETURN
      END SUBROUTINE Setcona

