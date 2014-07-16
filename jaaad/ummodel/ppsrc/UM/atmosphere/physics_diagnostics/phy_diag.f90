
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate various diagnostics related to dynamics variables.
!
! Subroutine Interface:
      SUBROUTINE Phy_diag(                                              &
! Primary data: in
     & p_star,p,rho,u,v,w,theta,q                                       &
     &,p_theta_levels,exner_rho_levels,exner_theta_levels               &
! Grid sizes and definition: in
     &,rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,global_row_length, global_rows                                   &
     &,eta_theta_levels,eta_rho_levels                                  &
     &,r_theta_levels,r_rho_levels                                      &
     &,delta_lambda,delta_phi,sec_theta_latitude                        &
! Control information: in
     &,Model_domain,npmsl_height                                        &
! Pressure levels for output arrays: in
     &,hts_p_press,t_p_press                                            &
     &,RHice_p_press,RHwat_p_press,wbpt_p_press                         &
! Flags to request each diagnostic output field: in
     &,qT_m                                                             &
     &,qhts_theta                                                       &
     &,qhts_p,qt_p,qRHice_p,qwbpt_p                                     &
     &,qp_MSL                                                           &
     &,qhts_rho                                                         &
     &,qRHwat_p                                                         &
! Diagnostics lengths: in
     &,hts_p_levs,t_p_levs                                              &
     &,RHice_p_levs,RHwat_p_levs,wbpt_p_levs                            &
! Diagnostic arrays: out
     &,T                                                                &
     &,hts_theta                                                        &
     &,hts_p,t_p,RHice_p,wbpt_p                                         &
     &,p_MSL                                                            &
     &,hts_rho                                                          &
     &,RHwat_p                                                          &
     &     )








      IMPLICIT NONE
!
! Description:
!   Calculate physics-related diagnostics - held in STASH section 16 -
!   which may include interpolation onto pressure surfaces. Diagnostics
!   currently supported:
!   STASH item
!     4 temperature on model levels
!   201 geopotential height on theta levels
!   202 geopotential height on pressure surfaces
!   203 temperature         on pressure surfaces
!   204 relative humidity wrt ice on pressure surfaces
!   205 wet bulb potential temperature on pressure surfaces
!   222 mean sea level pressure
!   255 geopotential height on rho levels
!   256 relative humidity wrt water on pressure surface
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Primary model data is input, and each diagnostic is calculated - if
!   its flag is set - in simple sequential order. Where the extraction
!   of the diagnostic quantity requires further calculation, a lower
!   level diagnostic-specific routine is called.
!
! Current Code Owner: <Rick Rawlins>
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.0  02/06/99   Extensive revision of PHYDIA1A deck for 'C-P C
!                  dynamics upgrade' project. Rick Rawlins.
!  5.1  14/02/00   Correct interface to routine vert_interp_mdi2 for
!                  RH diagnostic (16,204).
!                  Replace rmdi_pp by rmdi. R Rawlins
!  5.1  17/02/00   Add diagnostic for T on model levels. R Rawlins
!  5.2  19/01/01   Add diagnostic for wet bulb potential temperature
!                  (16,205). R Rawlins
!  5.2  19/03/01   Linear interpolate in exner instead of height.
!                  (T on p levels). C. Wilson
!  5.2  19/03/01   Linear interpolate in exner instead of height.
!                  C. Wilson
!  5.2  19/03/01   Change to use exner pressures instead of p in
!                  call to vert_interp_mdi2 (now vert_interp2)
! 6/6/01   5.3         Pass gc_all_proc_group to CALC_PMSL  P.Burton
!  5.3  07/11/01      Add Clive Wilson change to PMSL calc    A.Malcolm
!  5.3  19/11/01   Subroutine Thetaw has been generalized.
!                  Population of pressure array now done in PHYDIAG.
!                  D.M. Goddard
!  5.3  05/12/01   Add geopotential height on theta levels (16,201)
!                  and geopotential height on rho levels (16,255)
!                  D.M. Goddard
!  5.4  10/04/02   Add Relative Humidity wrt water on presure levels
!                  (16,256). D.M. Goddard
!  6.0  26/09/03   Fix reference to uninitialised height array
!                  when (16,255) is requested without (16,202).
!                  Restrict (16,201) to 1,model_levels.
!                  P. Selwood
!  6.1  08/07/04   send P_theta_levels into vert_h_onto_p
!                                                    Michael Hughes
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!   Documentation: UMDP 80
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):
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
! ----------------------- include file: INTERPOR -----------------------
! Description: Hold magic numbers for order of vertical interpolation,
!              where prime use is for atmosphere physics/dynamics
!              diagnostics.
!              interp_order should be set to a _default value where
!              used generally in the code, but can be set explicitly
!              locally where required.
!
! Current Code Owner: R. Rawlins
! History:
! Version  Date      Comment.
!  5.0 19/06/99  Original version. R. Rawlins.

      INTEGER,PARAMETER:: interp_order_linear  = 1
      INTEGER,PARAMETER:: interp_order_cubic   = 3
      INTEGER,PARAMETER:: interp_order_quintic = 5
      INTEGER,PARAMETER:: interp_order_default = interp_order_cubic

! INTERPOR end

! Subroutine arguments
!   Scalar arguments with intent(in):
! Grid sizes:
      INTEGER                                                           &
     & rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,global_row_length, global_rows                                   &
! Control information: in
     &,Model_domain                                                     &
                       ! Domain of atmosphere model:
!                                   global,LAM,cyclic LAM,single column

! Diagnostics lengths: IN
     &,hts_p_levs                                                       &
                      ! NO OF LEVS ON WHICH TO INTERP hts_p
     &,t_p_levs                                                         &
                      ! NO OF LEVS ON WHICH TO INTERP t_p
     &,RHice_p_levs                                                     &
                      ! NO OF LEVS ON WHICH TO INTERP RHice_p
     &,RHwat_p_levs                                                     &
                      ! NO OF LEVS ON WHICH TO INTERP RHwat_p
     &,wbpt_p_levs    ! NO OF LEVS ON WHICH TO INTERP wbpt_p

! Grid definition: IN
      REAL                                                              &
     & delta_lambda,delta_phi                                           &
! Orographic height threshold for new pmsl calculation: IN
     &,npmsl_height

! Flags to request each diagnostic output field: IN
      LOGICAL                                                           &
     & qhts_p                                                           &
                   ! Flag for geopotential heights   on pressure levels
     &,qt_p                                                             &
                   ! Flag for temperature            on pressure levels
     &,qRHice_p                                                         &
                   ! Flag for relatice humidity wrt ice
                   !   on pressure surfaces
     &,qRHwat_p                                                         &
                   ! Flag for relative humidity wrt water
                   !   on pressure surfaces
     &,qwbpt_p                                                          &
                   ! Flag for wet bulb potential T   on pressure levels
     &,qp_MSL                                                           &
                   ! Flag for mean sea level pressure
     &,qT_m                                                             &
                   ! Flag for temperature on model levels
     &,qhts_theta                                                       &
                   ! Flag for geopotential heights   on theta levels
     &,qhts_rho    ! Flag for geopotential heights   on rho levels
!   Array  arguments with intent(in):
! Primary data: IN
      REAL                                                              &
     & p_star(row_length, rows)                                         &
     &,p    (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,rho  (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,u    (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,v    (1-offx:row_length+offx, 1-offy:n_rows+offy,  model_levels) &
     &,w    (1-offx:row_length+offx, 1-offy:rows+offy,  0:model_levels) &
     &,theta(1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,q (1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,  wet_levels) &
     &,p_theta_levels                                                   &
     &      (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,exner_rho_levels                                                 &
     &      (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &
     &,exner_theta_levels                                               &
     &      (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels) &


! Vertical grid definition: IN
     &,eta_theta_levels(0:model_levels)                                 &
                                        ! vertical grid for theta vars
     &,eta_rho_levels    (model_levels)                                 &
                                        ! vertical grid for rho   vars
     &,r_theta_levels  (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j,       0:model_levels)     &
     &,r_rho_levels    (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j,         model_levels)     &
! Derived from horizontal grid: IN
     &,sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)      &

! Pressure levels for output arrays: IN
     &,hts_p_press(hts_p_levs)                                          &
                                      ! for heights      on p surfaces
     &,t_p_press  (  t_p_levs)                                          &
                                      ! for temperature  on p surfaces
     &,RHice_p_press ( RHice_p_levs)                                    &
                                      ! for r.h. wrt ice on p surf
     &,RHwat_p_press ( RHwat_p_levs)                                    &
                                      ! for r.h. wrt water on p surf
     &,wbpt_p_press(wbpt_p_levs)      ! for wet bulb T   on p surfaces

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
! Diagnostic arrays: OUT
      REAL                                                              &
     & hts_p(row_length,rows,hts_p_levs)                                &
                                         ! heights           at p levels
     &,t_p  (row_length,rows,t_p_levs)                                  &
                                         ! temperature       at p levels
     &,RHice_p (row_length,rows,RHice_p_levs)                           &
                                              ! rh wrt ice at p levels
     &,RHwat_p (row_length,rows,RHwat_p_levs)                           &
                                              ! rh wrt water at p-levels
     &,wbpt_p(row_length,rows,wbpt_p_levs)                              &
                                          !wet bulb pot temp at p levels
     &,p_MSL(row_length,rows)                                           &
                                         ! pressure at mean sea level
     &,T(row_length,rows,model_levels)                                  &
                                         ! temperature on model levels
     &,hts_theta(row_length,rows,model_levels)                          &
                                         ! heights on theta levels
     &,hts_rho(row_length,rows,model_levels)
                                         ! heights on rho levels
! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='Phy_diag')
      LOGICAL, PARAMETER ::                                             &
     & l_potential=.true.              ! Wet bulb potential temperature
                                       ! required from subroutine ThetaW
! Local scalars:
!   ErrorStatus
      INTEGER      ErrorStatus                                          &
                                        ! Error flag (0 = OK)
     &,i,j,k                                                            &
                                        ! loop counters
     &,interp_order    !  order of vertical interpolation

      CHARACTER*256                                                     &
     & CMessage         ! Error message if return code >0

      REAL                                                              &
     & dummy                                                            &
                                ! dummy argument - not referenced
     &,pressure_pa                                                      &
                                ! pressure in pascals
     &,pressure_ex              ! exner pressure

      ! Arguments for FlumeSendDiag
      INTEGER im ! model
      INTEGER is ! section
      INTEGER ie ! item
      INTEGER levels ! number of vertical levels

! Local dynamic arrays:
      REAL                                                              &
     & p_at_theta                                                       &
     &       (row_length,rows,model_levels)                             &
                                            ! pressure at theta points
     &,rh    (row_length,rows,model_levels)                             &
                                            ! workspace for RH
     &,height(row_length,rows,model_levels)                             &
                                            ! height at rho levels
     &,work1  (row_length,rows)                                         &
                                            ! temporary space
     &,work2  (row_length,rows)                                         &
                                            ! temporary space
     &,work3 (row_length,rows)              ! temporary space

! Function & Subroutine calls:
      External Calc_PMSL,Ereport,Qsat,Thetaw,T_vert_interp_to_p         &
     &,vert_interp2,vert_h_onto_p

!- End of header


! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

! set Error code to zero
      ErrorStatus = 0

! Set order of vertical interpolation
      interp_order = interp_order_linear

!  Calculate p at theta points. Store in p_at_theta

      DO k = 1, model_levels
        DO j = 1, rows
        DO i = 1, row_length
             p_at_theta(i,j,k) = p_theta_levels(i,j,k)
        ENDDO ! i
        ENDDO ! j
      ENDDO ! k

!   Remove radius of earth from rho levels to create geopotential height
!   of rho level.
!   Required for either 202 (geopotential height on pressure levels) or
!                       255 (geopotential height on rho levels)
      IF (qhts_p .OR. qhts_rho) THEN

      DO k = 1, model_levels
        DO j = 1, rows
        DO i = 1, row_length
             height(i,j,k) = r_rho_levels(i,j,k) - earth_radius
        ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      ENDIF ! on relevant STASHflags

!   Calculate temperature at theta points
      IF(qt_p .OR. qRHice_p .OR. qRHwat_p .OR. qT_m .OR. qwbpt_p) THEN

      DO k = 1, model_levels
        DO j = 1, rows
        DO i = 1, row_length
             T(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
        ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      ENDIF ! on relevant STASHflags

! ----------------------------------------------------------------------
! STASH item 201: geopotential height      on  theta levels
! ----------------------------------------------------------------------
      IF(qhts_theta) THEN

!   Remove radius of earth from height field

        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              hts_theta(i,j,k) = r_theta_levels(i,j,k) - earth_radius
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

      ENDIF ! on relevant STASHflags

      IF(qhts_rho) THEN
! ----------------------------------------------------------------------
! STASH item 255: geopotential height on  rho levels
! ----------------------------------------------------------------------

        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              hts_rho(i,j,k) = height(i,j,k)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

      ENDIF ! on relevant STASHflags
! ----------------------------------------------------------------------
! STASH item 202: geopotential height      on pressure surfaces
! ----------------------------------------------------------------------
      IF(qhts_p) THEN

        DO k = 1, hts_p_levs

         pressure_pa = hts_p_press(k)*100. ! convert to Pascals
! DEPENDS ON: vert_h_onto_p
         CALL vert_h_onto_p(                                            &
     &        height(1,1,1), row_length,rows, model_levels              &
     &       ,pressure_pa,r_rho_levels, r_theta_levels                  &
     &, p_theta_levels                                                  &
     &       ,theta, exner_theta_levels, exner_rho_levels               &
     &       ,R, g, Lapse,bl_levels                                     &
     &       ,offx, offy, halo_i, halo_j                                &
     &       ,p, interp_order,kappa, p_zero,cp, hts_p(1,1,k) )

        ENDDO ! over output STASH pressure levels


      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 203: temperature              on pressure surfaces
! ----------------------------------------------------------------------
      IF(qt_p) THEN

        DO k = 1, t_p_levs

         pressure_pa = t_p_press(k)*100. ! convert to Pascals
! DEPENDS ON: t_vert_interp_to_p
         CALL T_vert_interp_to_p(                                       &
     &        T, theta, row_length, rows                                &
     &        ,model_levels, pressure_pa, offx, offy,halo_i,halo_j      &
     &        ,p_theta_levels, Lapse, R, g                              &
     &        ,bl_levels                                                &
     &        ,exner_theta_levels                                       &
     &        ,r_theta_levels                                           &
     &        ,kappa, p_zero,  T_p(1,1,k) )

        ENDDO ! over output STASH pressure levels


      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 204, 256 : RH: relative humidity on pressure surfaces
! ----------------------------------------------------------------------

! STASH item 204: Relative humidity wrt ice
      IF(qRHice_p) THEN

        DO k = 1, wet_levels
!  Find humidity saturation at theta points - store in rh
! DEPENDS ON: qsat
          CALL QSAT(rh(1,1,k),T(1,1,k),p_at_theta(1,1,k),               &
     &              theta_field_size)
!  And convert to relative humidity
          DO j = 1, rows
            DO i = 1, row_length

              rh(i,j,k) = q(i,j,k)/rh(i,j,k)*100.
!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:

              IF(rh(i,j,k) <  0.0) THEN
                rh(i,j,k) = 0.
              END IF
            END DO ! i
          END DO ! j

        END DO ! k wet_levels

!  Interpolate
        DO k = 1, RHice_p_levs

          pressure_pa = RHice_p_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2(                                            &
     &          rh, row_length, rows, wet_levels                        &
     &        , pressure_ex                                             &
     &        , 0, 0, offx, offy                                        &
     &        , exner_theta_levels, interp_order                        &
     &        , RHice_p(1,1,k) )

        ENDDO ! k over output STASH pressure levels


      END IF  ! on STASHflag

! STASH item 256: Relative humidity wrt water
      IF(qRHwat_p) THEN

        DO k = 1, wet_levels
!  Find humidity saturation at theta points - store in rh
! DEPENDS ON: qsat_wat
          CALL QSAT_WAT(rh(1,1,k),T(1,1,k),p_at_theta(1,1,k),           &
     &              theta_field_size)
!  And convert to relative humidity
          DO j = 1, rows
            DO i = 1, row_length

              rh(i,j,k) = q(i,j,k)/rh(i,j,k)*100.
!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:

              IF(rh(i,j,k) <  0.0) THEN
                rh(i,j,k) = 0.
              END IF
            END DO ! i
          END DO ! j

        END DO ! k wet_levels

!  Interpolate
        DO k = 1, RHwat_p_levs

          pressure_pa = RHwat_p_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2(                                            &
     &          rh, row_length, rows, wet_levels                        &
     &        , pressure_ex                                             &
     &        , 0, 0, offx, offy                                        &
     &        , exner_theta_levels, interp_order                        &
     &        , RHwat_p(1,1,k) )

        ENDDO ! k over output STASH pressure levels


      END IF  ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 205: wet bulb potential temperature on pressure surfaces
! ----------------------------------------------------------------------
      IF(qwbpt_p) THEN

        DO k = 1, wbpt_p_levs

         pressure_pa = wbpt_p_press(k)*100. ! convert to Pascals

! Interpolate T onto required pressure level (work1)
! DEPENDS ON: t_vert_interp_to_p
         CALL T_vert_interp_to_p(                                       &
     &        T, theta, row_length, rows                                &
     &        ,model_levels, pressure_pa, offx, offy, halo_i, halo_j    &
     &        ,p_theta_levels, Lapse, R, g                              &
     &        ,bl_levels                                                &
     &        ,exner_theta_levels, r_theta_levels                       &
     &        ,kappa, p_zero, work1 )

! Interpolate q onto required pressure level (work2)
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2(                                            &
     &         q, row_length, rows, wet_levels                          &
     &        ,pressure_ex                                              &
     &        ,halo_i, halo_j, offx, offy                               &
     &        ,exner_theta_levels, interp_order                         &
     &        ,work2 )

! Generate pressure array for required pressure level (work3)
         DO j = 1, rows
           DO i = 1, row_length
             work3(i,j) = pressure_pa
           END DO
         END DO
! DEPENDS ON: thetaw
         CALL Thetaw(                                                   &
     &        theta_field_size,work1,work2,work3,l_potential,           &
                                                               ! in
     &        wbpt_p(1,1,k))                                  ! out

        ENDDO ! over output STASH pressure levels



      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 222: mean sea level pressure
! ----------------------------------------------------------------------
      IF(qp_MSL) THEN

! DEPENDS ON: calc_pmsl
         CALL Calc_PMSL(theta, exner_theta_levels, p                    &
     &                 ,r_theta_levels, r_rho_levels, eta_theta_levels  &
     &                 ,g, R, Lapse, earth_radius                       &
     &                 ,row_length, rows, model_levels, bl_levels       &
     &                 ,offx, offy, halo_i, halo_j                      &
     &                 , p_MSL, p_star                                  &
     &                 ,mype, nproc, nproc_x, nproc_y                   &
     &                 ,g_datastart                                     &
     &                 ,neighbour, at_extremity                         &
     &                 ,gc_all_proc_group, model_domain                 &
     &                 ,delta_lambda, delta_phi                         &
     &                 ,Cp, npmsl_height, sec_theta_latitude            &
     &                 ,global_row_length, global_rows)


      ENDIF ! on STASHflag

! ---------------------------------------------------------------------
! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      RETURN
      END SUBROUTINE Phy_diag
