
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate various diagnostics related to dynamics variables.
!
! Subroutine Interface:
      SUBROUTINE Dyn_diag(                                              &
! Primary data: in
     & exner_rho_levels,rho,u,v,w                                       &
     &,exner_theta_levels                                               &
     &,theta                                                            &
! Grid sizes and definition: in
     &,rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,global_rows,global_row_length                                    &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,eta_theta_levels,eta_rho_levels                                  &
     &,Model_domain                                                     &
! Grid coordinates: in
     &,delta_lambda,delta_phi                                           &
     &,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole      &
! Pre-calculated grid associated arrays: in
     &,r_at_u,r_at_v                                                    &
     &, r_theta_levels, r_rho_levels, sec_v_latitude                    &
     &, tan_v_latitude, sec_theta_latitude, f3_at_v                     &
     &,rot_coeff1,rot_coeff2                                            &
! Time information: in
     &,forecast_hrs                                                     &
! Theta levels for output arrays
     &, desired_theta                                                   &
! Pressure levels for output arrays: in
     &,ucomB_press,vcomB_press                                          &
     &,ucomp_press,vcomp_press,wcomp_press                              &
     &,testd_press                                                      &
     &,p_height,theta_height,rho_height                                 &
     &,w_height,u_height,v_height                                       &
     &,pv_press                                                         &
! Model levels for output arrays: in
     &,ucomB_model,vcomB_model                                          &
     &,testd_model                                                      &
     &,Htheta_model,Hrho_model                                          &
! Flags to request each diagnostic output field: in
! wind related diagnostics
     &,qucomB_m,qvcomB_m                                                &
     &,qucomB_p,qvcomB_p                                                &
     &,qucomp_p,qvcomp_p,qwcomp_p                                       &
     &,qu50mB_h,qv50mB_h                                                &
     &,qu50m_h,qv50m_h                                                  &
! PV related diagnostics
     &,qpotn_vort_theta,qtheta_potn_vort,qtheta_pv_points               &
     &,qpv_mod_levs,qpv_theta_mlev,qpotn_vort_press                     &
! test fields
     &,qdia1,qdia2,qdia3,qdia4                                          &
! flux diagnostics
     &,qrhow,qrhouw,qrhovw,qrhow_up                                     &
     &,qrhow_down,qrhowc_up,qrhowc_down                                 &
! height and height level diagnostics
     &,qHtheta_ml,qHrho_ml                                              &
     &,qpress_h,qtheta_h,qrho_h                                         &
     &,qu_h,qv_h,qw_h                                                   &
! other diagnostics
     &,spec_w,qtrue_density                                             &
! Flags for wind rotation (lam grid): in
     &,rot_uvcomB_p                                                     &
! Diagnostics lengths: in
     &,ucomB_m_levs,vcomB_m_levs                                        &
     &,ucomB_p_levs,vcomB_p_levs                                        &
     &,ucomp_p_levs,vcomp_p_levs,wcomp_p_levs                           &
     &,pv_theta_levs,pv_press_levs                                      &
     &,testd_p_levs,testd_m_levs                                        &
     &,Htheta_m_levs,Hrho_m_levs                                        &
     &,p_h_levs,theta_h_levs,rho_h_levs,w_h_levs,u_h_levs,v_h_levs      &
! Diagnostic arrays: out
! wind related diagnostics
     &,ucomB_m,vcomB_m                                                  &
     &,ucomB_p,vcomB_p                                                  &
     &,ucomp_p,vcomp_p,wcomp_p                                          &
     &,u50mB_h,v50mB_h                                                  &
     &,u50m_h,v50m_h                                                    &
! PV related diagnostics
     &,potn_vort_theta,theta_potn_vort,theta_pv_points                  &
     &,pv_mod_levs,pv_theta_mlev,potn_vort_press                        &
! test fields
     &,testdiag1,testdiag2,testdiag3,testdiag4                          &
! flux diagnostics
     &,rhow,rhouw,rhovw,rhow_up                                         &
     &,rhow_down,rhow_convup,rhow_convdown                              &
! height and height level diagnostics
     &,height_theta_ml,height_rho_ml                                    &
     &,press_h,theta_h,rho_h                                            &
     &,ucomp_h,vcomp_h,wcomp_h                                          &
! other diagnostics
     &,spec_3D,true_density)








      IMPLICIT NONE
!
! Description:
!   Calculate dynamics-related diagnostics - held in STASH section 15 -
!   which may include interpolation onto pressure surfaces. Diagnostics
!   currently supported:
!   [All diagnostics on native 'C' grid unless specified.]
!   STASH item
!     2 u component of wind on model levels      on 'B' grid
!     3 v component of wind on model levels      on 'B' grid
!   201 u component of wind on pressure surfaces on 'B' grid
!   202 v component of wind on pressure surfaces on 'B' grid
!   243 u component of wind on pressure surfaces
!   244 v component of wind on pressure surfaces
!   242 w component of wind on pressure surfaces
!   245 u component of wind at 50m height
!   246 v component of wind at 50m height
!   212 u component of wind at 50m height on 'B' grid
!   213 v component of wind at 50m height on 'B' grid
!   214 potential vorticity on theta levels
!   229 potential vorticity on pressure levels
!   215 theta on potential vorticity = +/-2 surface
!   216 theta at potential vorticity points
!   217 potential vorticity on model levels
!   218 potential vorticity on model theta grid and theta levels
!   231 test analytic field on v grid - single level
!   232 test analytic field on p grid - single level
!   233 test analytic field on p grid - pressure levels
!   234 test analytic field on p grid - model levels
!   260 mass flux (rhow) on model levels
!   261 momentum flux (rhouw) on model levels
!   262 momentum flux (rhovw) on model levels
!   263 upward mass flux (rhow, w> 0m/s) on model levels
!   264 downward mass flux (rhow, w<0m/s) on model levels
!   265 upward convective mass flux (rhow, w>1m/s) on model levels
!   266 downward convective mass flux (rhow, w<-1m/s) on model levels
!   101 Height from sea level on theta model levels
!   102 Height from sea level on rho model levels
!   108 Pressure on height (from sea) levels
!   119 Potential temperature on height (from sea) levels
!   127 Rho (density) on height (from sea) levels
!   142 W component of wind on height (from sea) levels
!   143 U component of wind on height (from sea) levels
!   144 V component of wind on height (from sea) levels
!   270 real part of spectra (model levels)
!   271 true unscaled density on model levels
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Primary model data is input, and each diagnostic is calculated - if
!   its flag is set - in simple sequential order. Where the extraction
!   of the diagnostic quantity requires further calculation, a lower
!   level diagnostic-specific routine is called.
!
! STASH items 002,003     : u,v wind components on model levs 'B' grid:
!    Simple horizontal interpolation from u,v to uv 'B' position.
!
! STASH items 201,202     : u,v wind components on p surfaces 'B' grid:
! 1. Simple horizontal interpolation of p onto u,v staggering.
! 2. Perform vertical interpolation of u,v onto each o/p p surface.
! 3. Simple horizontal interpolation from u,v to uv 'B' position.
! 4. (Lam only) Rotate winds from native eq lat-long to standard grid
!
! STASH item 214: pv on theta levels:
! 1. Interpolate theta onto PV points.
! 2. Calculate PV on model levels.
! 3. Perform vertical interpolation of PV onto theta levels.

! STASH items 243,244,242 : u,v,w wind components on pressure surfaces:
! 1. Simple horizontal interpolation of p onto u,v staggering.
! 2. Perform vertical interpolation of u,v onto each o/p p surface.
!
! STASH items 212,213 : u,v wind components at 50m height 'B' grid:
! STASH items 244,245 : u,v wind components at 50m height:
! 1. Re-use p_at_u/v arrays to get height surface 50m above orography
! 2. Perform vertical interpolation of u,v onto height surface
! 3. Simple horizontal interpolation from u,v to uv 'B' position.
!
! STASH items 231,232,233,234: Test diagnostics 1-4:
! Call TestDiag routine. See UM Doc Paper D7.
!
! STASH item 215: theta on pv=+/-2 surface:
! 1. Interpolate theta onto PV points.
! 2. Calculate PV on model levels.
! 3. Set mod_PV = |PV|
! 4. Perform vertical interpolation of theta onto PV=2 surface
!
! STASH item 216: theta at pv points:
! 1. Interpolate theta onto PV points.
!
! STASH item 217: pv on model levels:
! 1. Calculate PV on model levels.
!
! STASH item 218: pv on model theta levels:
! 1. Calculate PV on model levels, using 'calc_pv_at_theta'.
!
! Current Code Owner: <Rick Rawlins>
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.0  02/06/99   Extensive revision of DYNDIA1A deck for 'C-P C
!                  dynamics upgrade' project. Rick Rawlins.
!  5.1  25/01/00   Change u,v diagnostics from 'C' to 'B' grid. Retain
!                  old diagnostics but with new STASHcodes.
!                  Add u,v on model levels on 'B' grid.
!                  Replace rmdi_pp by rmdi. R Rawlins
!  5.2  14/09/00   Add PV on theta levels.  Z. Gardner
!  5.2  19/03/01   Change to use exner pressures instead of p in
!                  vertical interpolation and change
!                  call to vert_interp_mdi2 to vert_interp2
!  5.2  31/07/00   (1) Correction for (15,002-003): add
!                  swapbounds to update halos.
!                  (2) Introduce rotation of winds for a subset of
!                  diagnostics in lam. R Rawlins
!  5.2  15/11/00   Set up stash variables levels required for pv on
!                  theta levels diagnostic.  Zoe Gardner
!  5.3  27/04/01   Correct non-reprod u&v on pressure levels
!                  caused by acw4f502 changes.   S. Cusack
!  5.3  06/09/01   Set up stash variables levels required for pv on
!                  pressure levels diagnostic.  D.M. Goddard
!  5.4  28/05/02   Add theta on pv=2 or -2 surface also theta at pv poin
!                  and pv on model levels. T.J. Hinton
!  5.4  16/07/02   Introduce new diagnostics: u,v,w,theta,rho,p on
!                  geometric height levels and height (in meters) of
!                  theta, rho model levels from sea level.
!                               M. Diamantakis
!  5.4     28/08/02    Bug Fix (Bi-cyclic LAM)           Carol Roadnight
!  5.4  15/05/02   Correction to get bit comparison of potential
!                  vorticity diagnostics over different processor
!                  configurations. D.M. Goddard
!  5.4  11/03/02   Remove comment on same line as #include
!                                               S. Carroll
!  5.5  28/02/03   Add new mass flux diagnostics       Carol Roadnight
!  5.5  26/02/03   Get averaged spectra of vertical velocity field
!                                                       Carol Roadnight
!  6.0  18/07/03   Merge theta on PV = +/- 2 fields T.J.Hinton
!  6.1  17/05/04   Calc_PV arguments corrected. Adam Clayton
!  6.2  15/08/05   Free format fixes. P.Selwood
!  6.4  16/11/06   Get out true unscaled density          Andy Malcolm
!  6.4  16/11/06   Remove unnecessary swap_bounds calls   Andy Malcolm
!  6.4  16/11/06   Move & combine swap_bounds calls logic Andy Malcolm
!  6.4  16/11/06   Change code so only 1 call to CALC_PV  Andy Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
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

! Subroutine arguments
!   Scalar arguments with intent(in):
! Grid sizes:
      INTEGER                                                           &
     & rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &, global_rows,global_row_length                                   &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,Model_domain
! Grid coordinates: in
      REAL                                                              &
     & delta_lambda,delta_phi                                           &
     &,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole
! Time information: in
      REAL                                                              &
     & forecast_hrs       ! T+forecast_hrs: hours after analysis time
                          ! UM6.5 - MODEL_ANALYSIS_HRS changed to REAL -
                          ! requires FORECAST_HRS changed to REAL also  
! Diagnostics lengths: in
      INTEGER                                                           &
     & ucomB_m_levs                                                     & 
                          ! No of levels for output u on m 'B' grid
     &,vcomB_m_levs                                                     & 
                          ! No of levels for output v on m 'B' grid
     &,ucomB_p_levs                                                     & 
                          ! No of levels for output u on p 'B' grid
     &,vcomB_p_levs                                                     & 
                          ! No of levels for output v on p 'B' grid
     &,ucomp_p_levs                                                     &
                          ! No of levels for output of u on p
     &,vcomp_p_levs                                                     &
                          ! No of levels for output of v on p
     &,wcomp_p_levs                                                     &
                          ! No of levels for output of w on p
     &, pv_theta_levs                                                   &
                         !No of levels for output of pv on theta
     &, pv_press_levs                                                   &
                            !No of levels for output of pv on pressure
     &,testd_p_levs                                                     &
                          ! No of levels for output of testdiag3
     &,testd_m_levs                                                     &
                          ! No of levels for output of testdiag4
     &,Htheta_m_levs                                                    &
                              ! Num of levs for output H on theta lev
     &,Hrho_m_levs                                                      &
                              ! Num of levs for output H on rho lev
     &,p_h_levs                                                         &
                              ! Num of levs for output p on H levs
     &,theta_h_levs                                                     &
                              ! Num of levs for output theta on H levs
     &,rho_h_levs                                                       &
                              ! Num of levs for output rho on H levs
     &,w_h_levs                                                         &
                              ! Num of levs for output w on H levs
     &,u_h_levs                                                         &
                              ! Num of levs for output u on H levs
     &,v_h_levs               ! Num of levs for output v on H levs
! Flags to request each diagnostic output field: IN
      LOGICAL                                                           &
     & qucomB_m                                                         & 
                     ! Flag for U wind on model    levels  'B' grid
     &,qvcomB_m                                                         & 
                     ! Flag for V wind on model    levels  'B' grid
     &,qucomB_p                                                         & 
                     ! Flag for U wind on pressure levels  'B' grid
     &,qvcomB_p                                                         & 
                     ! Flag for V wind on pressure levels  'B' grid
     &,qucomp_p                                                         &
                     ! Flag for U wind component on pressure levels
     &,qvcomp_p                                                         &
                     ! Flag for V wind component on pressure levels
     &,qwcomp_p                                                         &
                     ! Flag for W wind component on pressure levels
     &,qu50mB_h                                                         & 
                     ! Flag for U wind comp at 50m height  'B' grid
     &,qv50mB_h                                                         & 
                     ! Flag for V wind comp at 50m height  'B' grid
     &,qu50m_h                                                          &
                     ! Flag for U wind component at 50m height
     &,qv50m_h                                                          &
                     ! Flag for V wind component at 50m height
     &, qpotn_vort_theta                                                &
                            !Flag for pv on theta levels
     &, qpotn_vort_press                                                &
                            !Flag for pv on pressure levels
     &,qtheta_potn_vort                                                 &
                           !Flag for theta on pv=+/-2 surface
     &,qtheta_pv_points                                                 &
                           !Flag for theta at pv points
     &,qpv_mod_levs                                                     &
                           !Flag for pv on model levels
     &,qpv_theta_mlev                                                   &
                           !Flag for pv on model theta points and levels
     &,qdia1                                                            &
                     ! Flag for test diagnostic 1
     &,qdia2                                                            &
                     ! Flag for test diagnostic 2
     &,qdia3                                                            &
                     ! Flag for test diagnostic 3
     &,qdia4                                                            &
                     ! Flag for test diagnostic 4
     &,qrhow                                                            &
                     ! Flag for rhow
     &,qrhouw                                                           &
                     ! Flag for rhouw
     &,qrhovw                                                           &
                     ! Flag for rhovw
     &,qrhow_up                                                         &
                     ! Flag for rhow_up
     &,qrhow_down                                                       &
                     ! Flag for rhow_down
     &,qrhowc_up                                                        &
                     ! Flag for rhow_convup
     &,qrhowc_down                                                      &
                     ! Flag for rhow_convdown
     &,qHtheta_ml                                                       &
                                 ! Flag for height on theta levs
     &,qHrho_ml                                                         &
                                 ! Flag for height on rho levs
     &,qpress_H                                                         &
                               ! Flag for press on height levs
     &,qtheta_H                                                         &
                               ! Flag for theta on height levs
     &,qrho_H                                                           &
                               ! Flag for rho on height levs
     &,qw_H                                                             &
                               ! Flag for w on height levs
     &,qu_H                                                             &
                               ! Flag for u on height levs
     &,qv_H                                                             &
                               ! Flag for v on height levs
     &,spec_w                                                           &
                               ! Flag for real part of spectra          
     &,qtrue_density           ! flag for true unscaled density
! Flags for wind rotation (lam grid)  (rotate if .T. ) : in
      LOGICAL, INTENT(IN) ::                                            &
     & rot_uvcomB_p          ! u,v (B grid) on p levels
!   Array  arguments with intent(in):
! Primary data: IN
      REAL                                                              &
     & u  (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels)   &
     &,v  (1-offx:row_length+offx, 1-offy:n_rows+offy,  model_levels)   &
     &,w  (1-offx:row_length+offx, 1-offy:rows+offy,  0:model_levels)   &
     &,rho(1-offx:row_length+offx, 1-offy:rows+offy,    model_levels)   &
     &,exner_rho_levels                                                 &
     &    (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)      &
     &,exner_theta_levels                                               &
     &    (1-offx:row_length+offx, 1-offy:rows+offy,    model_levels)   &
     &, theta(1-offx:row_length+offx, 1-offy:rows+offy,model_levels)    &
! Vertical grid definition: IN
     &,eta_theta_levels(0:model_levels)                                 &
                                        ! vertical grid for theta vars
     &,eta_rho_levels    (model_levels)                                 &
                                        ! vertical grid for rho   vars
! Pre-calculated grid associated arrays: IN
     &,r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:  rows+halo_j,      &
     &          model_levels)                                           &
     &,r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,      &
     &          model_levels)                                           &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, sec_v_latitude (1-offx:row_length+offx,1-offy:n_rows+offy)      &
     &, tan_v_latitude(row_length,n_rows)                               &
     &, sec_theta_latitude (1-offx:row_length+offx,1-offy:rows+offy)    &
     &, f3_at_v (1-offx:row_length+offx, 1-offy: n_rows+offy)           &
     &,rot_coeff1(row_length,n_rows)                                    &
                                      ! for lam wind rotations
     &,rot_coeff2(row_length,n_rows)                                    &
                                      !     (on B grid)
! Pressure levels (units mb) for output arrays: IN
     &,ucomB_press(ucomB_p_levs)                                        & 
                                        ! for u wind  'B' grid  press
     &,vcomB_press(vcomB_p_levs)                                        & 
                                        ! for v wind  'B' grid  press
     &,ucomp_press(ucomp_p_levs)                                        &
                                        ! for u wind
     &,vcomp_press(vcomp_p_levs)                                        &
                                        ! for v wind
     &,wcomp_press(wcomp_p_levs)                                        &
                                        ! for w wind
     &,testd_press(testd_p_levs)                                        &
                                        ! for test diagnostic3
     &,p_height(p_h_levs)                                               &
                                        ! for diagnostic 108 (press)
     &,theta_height(theta_h_levs)                                       &
                                        ! for diagnostic 119 (theta)
     &,rho_height(rho_h_levs)                                           &
                                        ! for diagnostic 127 (rho)
     &,w_height(w_h_levs)                                               &
                                        ! for diagnostic 142 (w)
     &,u_height(u_h_levs)                                               &
                                        ! for diagnostic 143 (u)
     &,v_height(v_h_levs)                                               &
                                        ! for diagnostic 144 (v)
     &, desired_theta(pv_theta_levs)                                    &
                                      ! for potential vorticity
                                      !  on theta levels
     &, pv_press(pv_press_levs)                                         &
                                      ! for potential vorticity
                                      !  on pressure levels
! Model    levels for output arrays: IN
     &,testd_model(testd_m_levs)        ! for test diagnostic4

      INTEGER                                                           &
     & ucomB_model(ucomB_m_levs)                                        & 
                                        ! for u wind  'B' grid  model
     &,vcomB_model(vcomB_m_levs)                                        & 
                                        ! for v wind  'B' grid  model
     &,Htheta_model(Htheta_m_levs)                                      &
                                        ! for diagnostic 101
     &,Hrho_model(Hrho_m_levs)          ! for diagnostic 102
!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
! Diagnostic arrays: OUT
      REAL                                                              &
     & ucomB_m(row_length,n_rows,ucomB_m_levs)                          & 
                                               ! u 'B' grid at mod levs
     &,vcomB_m(row_length,n_rows,vcomB_m_levs)                          & 
                                               ! v 'B' grid at mod levs
     &,ucomB_p(row_length,n_rows,ucomB_p_levs)                          & 
                                               ! u 'B' grid at pressures
     &,vcomB_p(row_length,n_rows,vcomB_p_levs)                          & 
                                               ! v 'B' grid at pressures
     &,ucomp_p(row_length,  rows,ucomp_p_levs)                          &
                                               ! u at selected pressures
     &,vcomp_p(row_length,n_rows,ucomp_p_levs)                          &
                                               ! v at selected pressures
     &,wcomp_p(row_length,rows,  wcomp_p_levs)                          &
                                               ! w at selected pressures
     &,u50mB_h(row_length,n_rows)                                       & 
                                               ! u at 50m ht 'B' grid
     &,v50mB_h(row_length,n_rows)                                       & 
                                               ! v at 50m ht 'B' grid
     &,u50m_h(row_length,  rows)                                        &
                                               ! u at 50m height
     &,v50m_h(row_length,n_rows)                                        &
                                               ! v at 50m height
     &, potn_vort_theta(row_length, n_rows, pv_theta_levs)              &
!                                  pv on theta levels
     &, potn_vort_press(row_length, n_rows, pv_press_levs)              &
                                               ! pv on pressure levels
     &,theta_potn_vort(row_length, n_rows)                              &
                                               ! theta on pv+/-2
     &,theta_pv_points(row_length,n_rows,model_levels)                  &
                                                       ! theta at pv poi
     &,pv_mod_levs(row_length,n_rows,model_levels)                      &
                                                   ! pv on model levels
     &,pv_theta_mlev(row_length,rows,model_levels)                      &
                                              ! pv on model theta levels
     &,testdiag1(row_length,n_rows)                                     &
                                               ! testdiag1 -single level
     &,testdiag2(row_length,  rows)                                     &
                                               ! testdiag2 -single level
     &,testdiag3(row_length,rows,testd_p_levs)                          &
                                               ! testdiag3 -press levels
     &,testdiag4(row_length,rows,testd_m_levs)                          &
                                               ! testdiag4 -model levels
     &,rhow(row_length,rows,model_levels)                               &
                                   ! mass flux on model levels
     &,rhow_up(row_length,rows,model_levels)                            &
                                   ! up mass flux on model levels
     &,rhow_down(row_length,rows,model_levels)                          &
                                   ! down mass flux on model levels
     &,rhow_convup(row_length,rows,model_levels)                        &
                                   ! up conv mass flux on model levels
     &,rhow_convdown(row_length,rows,model_levels)                      &
                                   ! down conv mass flux on model levels
     &,rhouw(row_length,rows,model_levels)                              &
                                             ! momentum flux
     &,rhovw(row_length,rows,model_levels)                              &
                                             ! momentum flux
! stash work arrays for height on theta and rho model levels
     &,height_theta_ml(row_length,rows,Htheta_m_levs)                   &
     &,height_rho_ml(row_length,rows,Hrho_m_levs)                       &
! stash work arrays for height level diagnostics
     &,press_h(row_length,rows,p_h_levs)                                &
     &,theta_h(row_length,rows,theta_h_levs)                            &
     &,rho_h(row_length,rows,rho_h_levs)                                &
     &,wcomp_h(row_length,rows,w_h_levs)                                &
     &,ucomp_h(row_length,rows,u_h_levs)                                &
     &,vcomp_h(row_length,n_rows,v_h_levs)                              &
     &,true_density(row_length,rows,model_levels) ! unscaled density


! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='Dyn_diag')
      REAL                                                              &
     & z_50m                    ! height=50m above earth's surface
      PARAMETER (                                                       &
     & z_50m = 50.                                                      &
     &)

! Local scalars:
!   ErrorStatus
      INTEGER      ErrorStatus                                          &
                                        ! Error flag (0 = OK)
     &,i,j,k                                                            &
                                        ! loop counters
     &,kk                                                               &
                                        ! k level index
     &,interp_order    !  order of vertical interpolation

      CHARACTER*256                                                     &
     & CMessage         ! Error message if return code >0

      LOGICAL                                                           &
     & LAM              ! T: limited area model

      REAL                                                              &
     & dummy                                                            &
                                ! dummy argument - not referenced
     &,desired_potn_vort                                                &
                                ! value of pv surface for theta
     &,pressure_pa                                                      &
                                ! pressure in pascals
     &,pressure_ex                                                      &
                                ! exner pressure
     &,weight                                                           &
                           ! weight for calculating w_on_rho
     &,desired_r              ! height to interpolate

! Local dynamic arrays:
      REAL                                                              &
     & exner_at_u(row_length,rows,  model_levels)                       &
                                                  !pressure at u points
     &,exner_at_v(row_length,n_rows,model_levels)                       &
                                                  !pressure at v points
     &,exner_at_w(row_length,rows,  model_levels)                       &
                                                  !pressure at w points
     &,exner_at_pv(row_length, n_rows, model_levels)                    &
                                               !pressure at pv points
![Note: exner_at_u,exner_at_v are optionally re-used as height points.]
     &,theta_at_PV(row_length, n_rows, model_levels)                    &
     &,T (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &,PV(row_length, n_rows, model_levels)                             &
     &,work_1(row_length,rows)                                          &
                                              ! workspace
     &,work_2(row_length,rows)                                          &
                                              ! workspace
     &,mod_PV(row_length, n_rows, model_levels)                         &
                                                 ! for pv=-2 diag
     &,w_on_rho(row_length,rows,model_levels)                           &
                                          ! vertical velocity on rho pts
     &,u_on_rho(row_length,rows,model_levels)                           &
                                              ! u velocity on rho pts
     &,v_on_rho(row_length,rows,model_levels)                           &
                                              ! v velocity on rho pts
     &,true_rho(row_length,rows,model_levels) ! rho/(r_rho_levels^2)

! Spectral variables

      Real                                                              &
     & spec_3D(row_length,rows,model_levels)                            &
     &,spec_2D(1-offx:row_length+offx,1-offy:rows+offy)

      Integer                                                           &
     & klev                                                             &
     &,gath_proc

! Function & Subroutine calls:
      External Ereport,vert_interp2,vert_interp_mdi,TestDiag            &
     &        ,uC_to_uB,vC_to_vB, calc_pv

!- End of header


! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

! Set Error code to zero
      ErrorStatus = 0

! Set order of vertical interpolation
      interp_order = interp_order_linear

! Determine whether limited area model
      IF (model_domain  ==  mt_lam .OR.                                 &
     &    model_domain  ==  mt_cyclic_lam .or.                          &
     &    model_domain  ==  mt_bi_cyclic_lam) THEN
        lam = .true.
      else
        lam = .false.
      ENDIF

!-----------------------------------------------
!    do all swapbounds calls needed
!-----------------------------------------------
      IF(qucomp_p .OR. qucomB_p .OR. qpotn_vort_press) THEN
! DEPENDS ON: swap_bounds
        call swap_bounds(exner_rho_levels,row_length,rows,              &
     &                 model_levels,offx,offy,fld_type_p,.false.)
      endif

      IF(qpotn_vort_theta .OR. qtheta_potn_vort .OR.                    &
     &   qtheta_pv_points .OR. qpotn_vort_press .OR. qpv_mod_levs .OR.  &
     &   qpv_theta_mlev) THEN  
! DEPENDS ON: swap_bounds
        CALL Swap_bounds(theta,row_length,rows,model_levels,            &
     &                   offx,offy,fld_type_p,.true.)
      endif

      IF(qucomB_m .OR. qpotn_vort_press .OR. qtheta_potn_vort .OR.      &
     &   qpv_mod_levs .OR. qrhouw .OR. qpotn_vort_theta .OR.            &
     &   qpv_theta_mlev) THEN  
! DEPENDS ON: swap_bounds
        CALL Swap_bounds(u,row_length,rows,model_levels,                &
     &                 offx,offy,fld_type_u,.true.)
      endif

      IF(qvcomB_m .OR. qpotn_vort_press .OR. qtheta_potn_vort .OR.      &
     &   qpv_mod_levs .OR. qrhovw .OR.qpotn_vort_theta .OR.             &
     &   qpv_theta_mlev) THEN  
! DEPENDS ON: swap_bounds
        CALL Swap_bounds(v,row_length,n_rows,model_levels,              &
     &                 offx,offy,fld_type_v,.true.)
      endif

      If (qpotn_vort_press .OR. qtheta_potn_vort .OR.                   &
     &    qpv_mod_levs .OR. qpotn_vort_theta .OR. qpv_theta_mlev) THEN
! DEPENDS ON: swap_bounds
        CALL Swap_bounds(rho,row_length,rows,model_levels,              &
     &                   offx,offy,fld_type_p,.true.)
      endif

      IF(qucomp_p.OR.qucomB_p) THEN
!  Calculate exner at u points. Store in exner_at_u
! Need call to swapbounds for exner_rho before interpolation

      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
             exner_at_u(i,j,k) = 0.5 *                                  &
     &             (exner_rho_levels(i,j,k) + exner_rho_levels(i+1,j,k))
           ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      ENDIF ! on STASHflag

      IF(qvcomp_p.OR.qvcomB_p) THEN
!  Calculate exner at v points. Store in exner_at_v

      DO k = 1, model_levels
        DO j = 1, n_rows
          DO i = 1, row_length
              exner_at_v(i,j,k) = 0.5 *                                 &
     &             (exner_rho_levels(i,j,k) + exner_rho_levels(i,j+1,k))
           ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      ENDIF ! on STASHflag

      IF(qwcomp_p) THEN
!  Calculate exner at w points. Store in exner_at_w

      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
             exner_at_w(i,j,k) = exner_theta_levels(i,j,k)
           ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      ENDIF ! on STASHflag

      IF(qpotn_vort_theta.OR.qtheta_potn_vort.OR.                       &
     &   qtheta_pv_points) THEN

! Calculate theta at PV points. Store in theta_at_pv
! first interpolate theta to rho levels. Use linear interpolation.
! Store in T as work space.

        Do j = 1-offy, rows+offy
          Do i = 1-offx, row_length+offx
             T(i,j,1) = theta(i,j,1)
          End Do
        End Do

        Do k = 2, model_levels
          Do j = 1-offy, rows+offy
            Do i = 1-offx, row_length+offx
                    T(i,j,k) = (theta(i,j,k) *                          &
     &                        (r_rho_levels(i,j,k) -                    &
     &                         r_theta_levels(i,j,k-1) ) +              &
     &                        theta(i,j,k-1) *                          &
     &                        (r_theta_levels(i,j,k) -                  &
     &                         r_rho_levels(i,j,k) ) ) /                &
     &                        (r_theta_levels(i,j,k) -                  &
     &                         r_theta_levels(i,j,k-1) )
            End Do
          End Do
        End Do

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
               theta_at_pv(i,j,k) = .25*(T(i+1,j,k) + T(i,j,k) +        &
     &                                   T(i+1,j+1,k) + T(i,j+1,k) )
            End Do
          End Do
        End Do
      ENDIF ! on STASHflag

      IF(qpotn_vort_press) THEN

! Interpolate exner onto PV points.

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
               exner_at_pv(i,j,k) = .25*(exner_rho_levels(i+1,j,k) +    &
     &                                   exner_rho_levels(i,j,k) +      &
     &                                   exner_rho_levels(i+1,j+1,k) +  &
     &                                   exner_rho_levels(i,j+1,k) )
            End Do
          End Do
        End Do
      ENDIF ! on STASHflag

! STASH items 002,003     : u,v wind components on model levs 'B' grid
! ----------------------------------------------------------------------

      IF(qucomB_m) THEN

        DO  k=1,ucomB_m_levs
! Perform simple horizontal interpolation from 'C' to 'B' grid
! Halos already populated, so interpolate directly:

          kk = ucomB_model(k) ! selected model level

          DO j=1,n_rows
          DO i=1,row_length
            ucomB_m(i,j,k) = (u(i,j,kk)+u(i,j+1,kk)) * 0.5
          ENDDO ! i
          ENDDO ! j

        ENDDO  ! k model levels loop


      ENDIF ! on STASHflag

      IF(qvcomB_m) THEN

        DO  k=1,vcomB_m_levs
! Perform simple horizontal interpolation from 'C' to 'B' grid
! Halos already populated, so interpolate directly:

          kk = vcomB_model(k) ! selected model level

          DO j=1,n_rows
          DO i=1,row_length
            vcomB_m(i,j,k) = (v(i,j,kk)+v(i+1,j,kk)) * 0.5
          ENDDO ! i
          ENDDO ! j

        ENDDO  ! k model levels loop


      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 201,202     : u,v wind components on p surfaces 'B' grid
! ----------------------------------------------------------------------

      IF(qucomB_p) THEN

        DO  k=1,ucomB_p_levs

          pressure_pa = ucomB_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2 (u, row_length, rows, model_levels          &
     &                          ,pressure_ex                            &
     &                          ,offx, offy, 0, 0                       &
     &                          ,exner_at_u, interp_order               &
     &                          ,work_1 )

! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: uc_to_ub
          CALL  uC_to_uB(work_1,                                        &
     &                   row_length,rows,n_rows,1,offx,offy,            &
     &                   ucomB_p(1,1,k))
        ENDDO  ! k pressure levels loop


      ENDIF ! on STASHflag

      IF(qvcomB_p) THEN
        DO  k=1,vcomB_p_levs

          pressure_pa = vcomB_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2 (v, row_length, n_rows, model_levels        &
     &                          ,pressure_ex                            &
     &                          ,offx, offy, 0, 0                       &
     &                          ,exner_at_v, interp_order               &
     &                          ,work_1 )

! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: vc_to_vb
          CALL  vC_to_vB(work_1,                                        &
     &                   row_length,n_rows,1,offx,offy,                 &
     &                   vcomB_p(1,1,k))

! Rotate winds from model to standard lat-long grid
          IF (lam .AND. rot_uvcomB_p) THEN
! First check valid requests: implicit assumption that u components
! and v components are requested for the same pressure levels
            IF(qucomB_p .AND. ucomB_press(k) == vcomB_press(k)) THEN

              DO j=1,n_rows
              DO i=1,row_length
                work_1(i,j) = ucomB_p(i,j,k)
                work_2(i,j) = vcomB_p(i,j,k)
              ENDDO
              ENDDO

! Rotation calculation on B grid
! DEPENDS ON: w_eqtoll
              CALL W_EqtoLL(rot_coeff1,rot_coeff2,work_1,work_2,        &
     &         ucomB_p(1,1,k),vcomB_p(1,1,k),v_field_size,v_field_size)

            ELSE

              ErrorStatus = -1        ! Warning
              Cmessage='wind diagnostics cannot be rotated: u and v '   &
     &           //'requested components on pressure levels must match'
! DEPENDS ON: ereport
              CALL Ereport(Routinename,ErrorStatus,Cmessage)
              exit                    ! jump out of levels loop

            ENDIF       ! Check valid request

          ENDIF  ! (lam wind rotation)

        ENDDO  ! k pressure levels loop


      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! Calculate PV for use with STASH items 229,214,215,217
! ----------------------------------------------------------------------

      If (qpotn_vort_press .OR. qtheta_potn_vort .OR.                   &
     &    qpv_mod_levs .or. qpotn_vort_theta) THEN

! Calculate PV on model levels.

! DEPENDS ON: calc_pv
        Call Calc_PV                                                    &
     &            (u, v, theta, rho,                                    &
     &             r_theta_levels, r_rho_levels,                        &
     &             r_at_u, r_at_v,                                      &
     &             sec_v_latitude, tan_v_latitude,                      &
     &             sec_theta_latitude, f3_at_v,                         &
     &             delta_lambda, delta_phi,                             &
     &             row_length, rows, n_rows, model_levels,              &
     &             offx, offy, halo_i, halo_j,                          &
     &             at_extremity,                                        &
     &             PV)

      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 229 : potential vorticity on pressure levels
! ----------------------------------------------------------------------

      If (qpotn_vort_press) THEN

        Do k = 1, pv_press_levs

          pressure_pa = pv_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          Call vert_interp2 (pv, row_length, n_rows,                    &
     &                       model_levels,                              &
     &                       pressure_ex,                               &
     &                       0, 0, 0, 0,                                &
     &                       exner_at_pv, interp_order_cubic,           &
     &                       potn_vort_press(1,1,k) )

        ENDDO  ! k pressure levels loop


      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 214 : potential vorticity on theta levels
! ----------------------------------------------------------------------

      If (qpotn_vort_theta) THEN

        Do k = 1, pv_theta_levs

! DEPENDS ON: vert_interp_mdi
          Call vert_interp_mdi (PV, row_length, n_rows,                 &
     &                          model_levels,                           &
     &                          desired_theta(k),                       &
     &                          0, 0, 0, 0,                             &
     &                          theta_at_pv, interp_order_cubic,        &
     &                          rmdi, potn_vort_theta(1,1,k) )

        ENDDO  ! k theta levels loop


      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 243,244,242 : u,v,w wind components on pressure surfaces
! ----------------------------------------------------------------------

      IF(qucomp_p) THEN

        DO  k=1,ucomp_p_levs

          pressure_pa = ucomp_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2 (u, row_length, rows, model_levels          &
     &                          ,pressure_ex                            &
     &                          ,offx, offy, 0, 0                       &
     &                          ,exner_at_u, interp_order               &
     &                          ,ucomp_p(1,1,k) )

        ENDDO  ! k pressure levels loop


      ENDIF ! on STASHflag

      IF(qvcomp_p) THEN
        DO  k=1,vcomp_p_levs

          pressure_pa = vcomp_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2 (v, row_length, n_rows, model_levels        &
     &                          ,pressure_ex                            &
     &                          ,offx, offy, 0, 0                       &
     &                          ,exner_at_v, interp_order               &
     &                          ,vcomp_p(1,1,k) )

        ENDDO  ! k pressure levels loop


      ENDIF ! on STASHflag

      IF(qwcomp_p) THEN
        DO  k=1,wcomp_p_levs

          pressure_pa = wcomp_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
! DEPENDS ON: vert_interp2
          CALL vert_interp2 (w, row_length, rows, model_levels          &
     &                          ,pressure_ex                            &
     &                          ,offx, offy, 0, 0                       &
     &                          ,exner_at_w, interp_order               &
     &                          ,wcomp_p(1,1,k) )

        ENDDO  ! k pressure levels loop


      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 212,213 : u,v wind components at 50m height 'B' grid
! STASH items 245,246 : u,v wind components at 50m height
! ----------------------------------------------------------------------

! Restrict interpolation to boundary levels, which should always
! greatly exceed levels in the vicinity of 50m.


      IF(qu50mB_h.OR.qu50m_h) THEN
! Generate height field above orography for lower (ie boundary) levels
! at u pts: (re-use array exner_at_u for workspace)
        DO k=1,bl_levels
          DO j=1,rows
            DO i=1,row_length
              exner_at_u(i,j,k)= r_at_u(i,j,k) - r_at_u(i,j,1)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

! DEPENDS ON: vert_interp_mdi
         CALL vert_interp_mdi (u, row_length, rows,                     &
     &                         bl_levels, z_50m,                        &
     &                         offx, offy,                              &
     &                         0, 0,                                    &
     &                         exner_at_u, interp_order,                &
     &                         rmdi, work_1 )
         IF(qu50m_h) THEN
           DO j=1,rows
             DO i=1,row_length
               u50m_h(i,j) = work_1(i,j)
             ENDDO ! i
           ENDDO ! j


         ENDIF ! on STASHflag qu50m_h

       IF(qu50mB_h) THEN
! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: uc_to_ub
            CALL  uC_to_uB(work_1,                                      &
     &                     row_length,rows,n_rows,1,offx,offy,          &
     &                     u50mB_h)

     
       ENDIF ! on STASHflag qu50mB_h

      ENDIF ! on STASHflags  qu50mB_h.OR.qu50m_h

      IF(qv50mB_h.OR.qv50m_h) THEN

! Generate height field above orography for lower (ie boundary) levels
! at v pts: (re-use array exner_at_v for workspace)
        DO k=1,bl_levels
          DO j=1,n_rows
            DO i=1,row_length
              exner_at_v(i,j,k)= r_at_v(i,j,k) - r_at_v(i,j,1)
            ENDDO ! i
          ENDDO ! j
        ENDDO ! k

! DEPENDS ON: vert_interp_mdi
         CALL vert_interp_mdi (v, row_length, n_rows,                   &
     &                         bl_levels, z_50m,                        &
     &                         offx, offy,                              &
     &                         0, 0,                                    &
     &                         exner_at_v, interp_order,                &
     &                         rmdi, work_1)
         IF(qv50m_h) THEN
           DO j=1,n_rows
             DO i=1,row_length
               v50m_h(i,j) = work_1(i,j)
             ENDDO ! i
           ENDDO ! j


         ENDIF ! on STASHflag qv50m_h

        IF(qv50mB_h) THEN
! Perform simple horizontal interpolation from 'C' to 'B' grid

! DEPENDS ON: vc_to_vb
            CALL  vC_to_vB(work_1,                                      &
     &                     row_length,n_rows,1,offx,offy,               &
     &                     v50mB_h)


         ENDIF ! on STASHflag qv50mB_h

      ENDIF ! on STASHflags qv50mB_h.OR.qv50m_h


! ----------------------------------------------------------------------
! STASH items 231,232,233,234: Test diagnostics 1-4
! ----------------------------------------------------------------------

      IF (qdia1.OR.qdia2.OR.qdia3.OR.qdia4) THEN

! DEPENDS ON: testdiag
        CALL TestDiag(                                                  &
     &  theta_field_size,v_field_size,rows,n_rows,row_length            &
     & ,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole     &
     & ,lam                                                             &
     & ,testd_press,testd_p_levs                                        &
     & ,testd_model,testd_m_levs,forecast_hrs                           &
     & ,testdiag1,testdiag2,testdiag3,testdiag4                         &
     & ,qdia1,qdia2,qdia3,qdia4)

      ENDIF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 215 : theta on potential vorticity = +/-2 surface
! ---------------------------------------------------------------------

      IF (qtheta_potn_vort) THEN
        desired_potn_vort = 0.0000020    ! set pv surface to pv=2x10^-6

! Take the absolute value of pv so that interpolation is to pv = +/- 2

        Do k = 1, model_levels
          Do j = 1, n_rows
            Do i = 1, row_length
              mod_PV(i,j,k) = ABS(PV(i,j,k))
            End Do
          End Do
        End Do

! Interpolate theta onto pv=2 surface

! DEPENDS ON: vert_interp_mdi
        Call vert_interp_mdi (theta_at_pv, row_length, n_rows,          &
     &                         model_levels,                            &
     &                         desired_potn_vort,                       &
     &                         0, 0, 0, 0,                              &
     &                         mod_PV, interp_order_linear,             &
     &                         rmdi, theta_potn_vort(1,1) )


      ENDIF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 216 : theta at potential vorticity points
! ---------------------------------------------------------------------

      IF (qtheta_pv_points) THEN

         Do k = 1, model_levels
            Do j = 1, n_rows
               Do i = 1, row_length
                  theta_pv_points(i,j,k) = theta_at_pv(i,j,k)
               End Do
            End Do
         End Do


      ENDIF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 217 : potential vorticity on model levels
! ---------------------------------------------------------------------

      IF (qpv_mod_levs) THEN

         Do k = 1, model_levels
            Do j = 1, n_rows
               Do i = 1, row_length
                  pv_mod_levs(i,j,k) = PV(i,j,k)
               End Do
            End Do
         End Do


      ENDIF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 218 : potential vorticity on model theta points and levels
! ---------------------------------------------------------------------

      IF (qpv_theta_mlev) THEN

! DEPENDS ON: calc_pv_at_theta
        Call Calc_PV_at_theta(u, v, theta, rho,                         &
     &                        r_theta_levels, r_rho_levels,             &
     &                        r_at_u, r_at_v,                           &
     &                        sec_v_latitude, tan_v_latitude,           &
     &                        sec_theta_latitude, f3_at_v,              &
     &                        delta_lambda, delta_phi,                  &
     &                        model_domain,                             &
     &                        pv_theta_mlev)


      ENDIF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 260-266 : mass and momentum fluxes
! 260 = mass flux = rhow
! 261 = rhouw
! 262 = rhovw
! 263 = upward mass flux = rhow with w >0m/s
! 264 = downward mass flux = rhow with w <0m/s
! 265 = upward convective mass flux = rhow with w> 1m/s
! 266 = upward convective mass flux = rhow with w < -1m/s
! ----------------------------------------------------------------------
      IF (qrhow .OR. qrhouw .OR. qrhovw .OR. qrhow_up .OR.              &
     &     qrhow_down .OR. qrhowc_up .OR. qrhowc_down) THEN
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              weight = (r_theta_levels(i,j,k) - r_rho_levels(i,j,k))/   &
     &                 (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
              w_on_rho(i,j,k)= weight        * w(i,j,k-1) +             &
     &                        (1.0 - weight) * w(i,j,k)
              true_rho(i,j,k)=rho(i,j,k)/                               &
     &                        (r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
              rhow(i,j,k) = true_rho(i,j,k)*w_on_rho(i,j,k)
              IF (qrhouw) THEN
                u_on_rho(i,j,k)=0.5*(u(i-1,j,k)+u(i,j,k))
                rhouw(i,j,k)=rhow(i,j,k)*u_on_rho(i,j,k)
              ENDIF
              IF (qrhovw) THEN
                v_on_rho(i,j,k)=0.5*(v(i,j-1,k)+v(i,j,k))
                rhovw(i,j,k)=rhow(i,j,k)*v_on_rho(i,j,k)
              ENDIF
              IF (qrhow_up) THEN
                IF (w_on_rho(i,j,k)  >   0.0) THEN
                  rhow_up(i,j,k)=rhow(i,j,k)
                ELSE
                  rhow_up(i,j,k)=0.0
                ENDIF
              ENDIF
              IF (qrhow_down) THEN
                IF (w_on_rho(i,j,k)  <   0.0) THEN
                  rhow_down(i,j,k) = rhow(i,j,k)
                ELSE
                  rhow_down(i,j,k)=0.0
                ENDIF
              ENDIF
              IF (qrhowc_up) THEN
                IF (w_on_rho(i,j,k)  >   1.0) THEN
                  rhow_convup(i,j,k)=rhow(i,j,k)
                ELSE
                  rhow_convup(i,j,k)=0.0
                ENDIF
              ENDIF
              IF (qrhowc_down) THEN
                IF (w_on_rho(i,j,k)  <   -1.0) THEN
                  rhow_convdown(i,j,k) = rhow(i,j,k)
                ELSE
                  rhow_convdown(i,j,k)=0.0
                ENDIF
              ENDIF
            End Do
          End Do
        End Do


      ENDIF ! on STASHflag
!-----------------------------------------------------------------------
! STASH items 101, 102, 108, 119, 127, 142, 143, 144
!-----------------------------------------------------------------------
! Height on theta-levels diagnostic

      If ( qHtheta_ml ) Then

         Do k=1, Htheta_m_levs
            Do j=1, rows
               Do i=1, row_length
                  kk = Htheta_model(k)
                  height_theta_ml(i,j,k)                                &
     &                        = r_theta_levels(i,j,kk) - Earth_Radius
               Enddo
            Enddo
         Enddo

      Endif

! Height on rho-levels diagnostic

      If ( qHrho_ml ) Then

         Do k=1, Hrho_m_levs
            Do j=1, rows
               Do i=1, row_length
                  kk = Hrho_model(k)
                  height_rho_ml(i,j,k)                                  &
     &                 = r_rho_levels(i,j,kk) - Earth_Radius
               Enddo
            Enddo
         Enddo

      End If

! pressure on height-levels diagnostic

      If ( qpress_H ) Then

         Do k=1, p_h_levs

            desired_r = p_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi (exner_rho_levels,row_length,rows,     &
     &           model_levels,desired_r,offx,offy,halo_i,halo_j,        &
     &           r_rho_levels,interp_order,rmdi,press_h(:,:,k) )

! Convert to standard pressure

            Do j=1, rows
               Do i=1, row_length
                  If ( press_h(i,j,k)  /=  rmdi ) Then
                     press_h(i,j,k) = p_zero *                          &
     &                    press_h(i,j,k)**(1./kappa)
                  End If
               End Do
            End Do

         End do  ! end k-loop

      End if

! potential temperature on height-levels diagnostic

      If ( qtheta_H ) Then

         Do k=1,theta_h_levs

            desired_r = theta_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi ( theta,row_length,rows,               &
     &           model_levels,desired_r,offx,offy,halo_i,halo_j,        &
     &           r_theta_levels(1-halo_i,1-halo_j,1),interp_order,      &
     &           rmdi,theta_h(:,:,k) )

         End do

      End if

! density (rho) on height-levels diagnostic

      If ( qrho_H ) Then

         Do k=1,rho_h_levs

            desired_r = rho_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi (rho,row_length,rows,model_levels,     &
     &           desired_r,offx,offy,halo_i,halo_j,r_rho_levels,        &
     &           interp_order,rmdi,rho_h(:,:,k))

! Convert to true density: divide by r_rho_levels**2

            Do j=1, rows
               Do i=1, row_length
                  If ( rho_h(i,j,k)  /=  rmdi ) Then
                     rho_h(i,j,k) = rho_h(i,j,k) /                      &
     &                    (r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
                  End if
               End do
            End do

         End do

      End if

! w on height-levels diagnostic

      If ( qw_H ) Then

         Do k=1, w_h_levs

            desired_r = w_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi (w,row_length,rows,model_levels+1,     &
     &           desired_r,offx,offy,halo_i,halo_j,r_theta_levels,      &
     &           interp_order,rmdi,wcomp_h(:,:,k) )

         End do

      End if

! u on height-levels diagnostic

      If ( qu_H ) Then

         Do k=1, u_h_levs

            desired_r = u_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi (u,row_length,rows,model_levels,       &
     &           desired_r,offx,offy,halo_i,halo_j,                     &
     &           r_at_u,interp_order,rmdi,ucomp_h(:,:,k))

         End do

      End if

! v on height-levels diagnostic

      If ( qv_H ) Then

         Do k=1, v_h_levs

            desired_r = v_height(k) + Earth_Radius

! DEPENDS ON: vert_interp_mdi
            Call vert_interp_mdi (v,row_length,n_rows,model_levels,     &
     &           desired_r,offx,offy,halo_i,halo_j,r_at_v,              &
     &           interp_order,rmdi,vcomp_h(:,:,k))

         End do

      End if

      IF (spec_w) then
        If (model_domain  ==  4) then
          gath_proc=0
          do klev=1,model_levels
! DEPENDS ON: calc_spectra
            CALL CALC_SPECTRA(w(1-offx,1-offy,klev),spec_2D,            &
     &             row_length+2*offx,rows+2*offy,                       &
     &             global_row_length,global_rows,                       &
     &             fld_type_p,halo_type_single,                         &
     &             gath_proc)
            do j = 1,rows
              do i = 1,row_length
                spec_3D(i,j,klev)=spec_2D(i,j)
              enddo
            enddo
          enddo
        Else
          write(6,*)'The spectra is not set up for this domain'
        Endif
      ENDIF

! -----------------------------------------------------
! stash item 271  : true density on model_levels
! -----------------------------------------------------
      if (qtrue_density) then
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              true_density(i,j,k)=rho(i,j,k)/                           &
     &                        (r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
            end do
          end do
        end do




      endif

! ---------------------------------------------------------------------
! Check error condition
      IF(ErrorStatus >  0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      ENDIF

      RETURN
      END SUBROUTINE Dyn_diag
