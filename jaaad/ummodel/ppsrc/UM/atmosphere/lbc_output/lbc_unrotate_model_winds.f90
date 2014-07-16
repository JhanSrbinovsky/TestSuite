
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Unrotates the model winds if on a rotated grid.
!
! Subroutine Interface:

      Subroutine LBC_Unrotate_Model_Winds (                             &
     &     u                                                            &
     &,    v                                                            &
     &,    u_size                                                       &
     &,    v_size                                                       &
     &,    row_length                                                   &
     &,    rows                                                         &
     &,    rows_v                                                       &
     &,    model_levels                                                 &
     &,    src_halo_type                                                &
     &,    src_pole_lat                                                 &
     &,    src_pole_long                                                &
     &,    src_first_lat                                                &
     &,    src_first_long                                               &
     &,    src_delta_lat                                                &
     &,    src_delta_long                                               &
     & )

      Implicit NONE
!
! Description:
!   Unrotates the model winds if on a rotated grid.
!
! Method:
!   1. Linear Interpolation of u-comp from u-grid to p-grid.
!   2. Linear Interpolation of v-comp from v-grid to p-grid.
!      (u and v now on same grid points)
!   3. Call W_EqToLL to unrotate the winds.
!   4. Linear Interpolation of u-comp from p-grid to u-grid.
!   5. Linear Interpolation of v-comp from p-grid to v-grid.
!      (u and v now back on u and v grids)
!
! Original Author : Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.5    16/10/02  Original code. Dave Robinson
!   6.1    07/12/04  Correct halosize used to calculate lambda_rot
!                    and phi_rot. Dave Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables :

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

! Subroutine arguments

      Integer :: u_size
      Integer :: v_size
      Integer :: row_length
      Integer :: rows
      Integer :: rows_v
      Integer :: model_levels
      Integer :: src_halo_type

      Real    :: src_delta_lat
      Real    :: src_delta_long
      Real    :: src_first_lat
      Real    :: src_first_long
      Real    :: src_pole_lat
      Real    :: src_pole_long

      Real :: u (u_size, model_levels)
      Real :: v (v_size, model_levels)

! Local parameters:

      Character (Len=*), Parameter ::                                   &
     &                   RoutineName = 'LBC_Unrotate_Model_Winds'

! Local scalars:

      Integer :: ipt, i, j   ! Loop indices
      Integer :: ij          ! Grid point number
      Integer :: points      ! No of grid points
      Integer :: halo_x      ! Halo size in EW direction
      Integer :: halo_y      ! Halo size in NS direction
      Integer :: level       ! Loop index
      Integer :: ErrorStatus ! Error Code

      Character (Len=80) :: CMessage

! Local dynamic arrays:

      Real, dimension (:), allocatable :: coeff3
      Real, dimension (:), allocatable :: coeff4
      Real, dimension (:), allocatable :: lambda_rot
      Real, dimension (:), allocatable :: phi_rot
      Real, dimension (:), allocatable :: lambda_true
      Real, dimension (:), allocatable :: phi_true
      Real, dimension (:), allocatable :: u_p
      Real, dimension (:), allocatable :: v_p
      Real, dimension (:), allocatable :: u_p_rot
      Real, dimension (:), allocatable :: v_p_rot
!
!- End of header

      ErrorStatus = 0
      CMessage = ' '

! --------------------
! Allocate work arrays
! --------------------

      points = g_lasize(1,fld_type_p,src_halo_type,mype) *              &
     &         g_lasize(2,fld_type_p,src_halo_type,mype)

      allocate ( coeff3      (points) )
      allocate ( coeff4      (points) )
      allocate ( lambda_rot  (points) )
      allocate ( phi_rot     (points) )
      allocate ( lambda_true (points) )
      allocate ( phi_true    (points) )

! -------------------
! Get halo dimensions
! -------------------

      halo_x = halosize (1, src_halo_type)
      halo_y = halosize (2, src_halo_type)

! -------------------------------------------
! Set up lat/longs for the rotated model grid
! -------------------------------------------

      ipt = 0
      Do j = 1,g_lasize(2,fld_type_p,src_halo_type,mype)
        Do i = 1,g_lasize(1,fld_type_p,src_halo_type,mype)

          ipt =ipt + 1

          Lambda_Rot(ipt) = Src_First_Long + Src_Delta_Long *           &
     &                     (I + g_datastart(1,mype) - halo_x - 2)

          Phi_Rot(ipt) = Src_First_Lat + Src_Delta_Lat *                &
     &                  (J + g_datastart(2,mype) - halo_y - 2)

        End Do
      End Do

! -------------------------------------
! Get true lat/longs for the model grid
! -------------------------------------

! DEPENDS ON: eqtoll
      Call EqToLL (Phi_Rot, Lambda_Rot, Phi_True, Lambda_True,          &
     &             Src_Pole_Lat, Src_Pole_Long, Points)

! ----------------------------------
! Compute wind rotation coefficients
! ----------------------------------

! DEPENDS ON: w_coeff
      Call W_Coeff (Coeff3, Coeff4, Lambda_True, Lambda_Rot,            &
     &              Src_Pole_Lat, Src_Pole_Long, Points)

      deallocate (lambda_rot)
      deallocate (phi_rot)
      deallocate (lambda_true)
      deallocate (phi_true)

! --------------------------------------
! Allocate work arrays for interpolation
! --------------------------------------

      allocate ( u_p     (points) )
      allocate ( v_p     (points) )
      allocate ( u_p_rot (points) )
      allocate ( v_p_rot (points) )

      Do level = 1, model_levels

! --------------------------------------------------
! Interpolate u from u-grid to p-grid ; u to u_p_rot
! --------------------------------------------------

! DEPENDS ON: lbc_u_to_p
        call lbc_u_to_p (u(1,level), u_p_rot,                           &
     &                   row_length, rows, halo_x, halo_y)

! --------------------------------------------------
! Interpolate v from v-grid to p-grid ; v to v_p_rot
! --------------------------------------------------

! DEPENDS ON: lbc_v_to_p
        call lbc_v_to_p (v(1,level), v_p_rot,                           &
     &                   row_length, rows, rows_v, halo_x, halo_y )

! -----------------------------------------------------
! Rotate winds to standard lat-lon ; u/v_p_rot to u/v_p
! -----------------------------------------------------

! DEPENDS ON: w_eqtoll
        call w_eqtoll(coeff3, coeff4, u_p_rot, v_p_rot, u_p, v_p,       &
     &                points, points)

! ---------------------------------------------------
! Interpolate u from p-grid back to u-grid ; u_p to u
! ---------------------------------------------------

! DEPENDS ON: lbc_p_to_u
        call lbc_p_to_u (u_p, u(1,level),                               &
     &                   row_length, rows, halo_x, halo_y)

! ---------------------------------------------------
! Interpolate v from p-grid back to v_grid ; v_p to v
! ---------------------------------------------------

! DEPENDS ON: lbc_p_to_v
        call lbc_p_to_v (v_p, v(1,level),                               &
     &                   row_length, rows, rows_v, halo_x, halo_y)

      End Do  !  Loop over levels

      deallocate ( u_p )
      deallocate ( v_p )
      deallocate ( u_p_rot )
      deallocate ( v_p_rot )
      deallocate ( coeff3 )
      deallocate ( coeff4 )

      Return
      END SUBROUTINE LBC_Unrotate_Model_Winds
