
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generate the LBC data from the model prognostics
!
! Subroutine Interface:

      SUBROUTINE MAKE_LBCS(                                             &
! Prognostic
     &  SOURCE_FIELD,                                                   &
     &  LOCAL_ROW_LENGTH,                                               &
     &  LOCAL_ROWS,                                                     &
     &  GLOBAL_ROW_LENGTH,                                              &
     &  GLOBAL_ROWS,                                                    &
     &  SOURCE_HALO_X,                                                  &
     &  SOURCE_HALO_Y,                                                  &
     &  SOURCE_LEVELS,                                                  &
     &  SOURCE_FLD_TYPE,                                                &
     &  SOURCE_HALO_TYPE,                                               &
     &  source_delta_lat,                                               &
     &  source_delta_long,                                              &
     &  source_first_lat,                                               &
     &  source_first_long,                                              &
     &  source_pole_lat,                                                &
     &  source_pole_long,                                               &
     &  source_cyclic,                                                  &
     &  source_rotated,                                                 &
! Orography Field
     &  OROG_FIELD,                                                     &
     &  OROG_LOCAL_ROW_LENGTH,                                          &
     &  OROG_LOCAL_ROWS,                                                &
     &  OROG_GLOBAL_ROW_LENGTH,                                         &
     &  OROG_GLOBAL_ROWS,                                               &
     &  OROG_FLD_TYPE,                                                  &
     &  OROG_HALO_TYPE,                                                 &
     &  OROG_HALO_X,                                                    &
     &  OROG_HALO_Y,                                                    &
     &  OROG_first_lat,                                                 &
     &  OROG_first_long,                                                &
! LBC field
     &  lbc_row_len,                                                    &
     &  lbc_rows,                                                       &
     &  lbc_levels,                                                     &
     &  lbc_delta_lat,                                                  &
     &  lbc_delta_long,                                                 &
     &  lbc_first_lat,                                                  &
     &  lbc_first_long,                                                 &
     &  lbc_pole_lat,                                                   &
     &  lbc_pole_long,                                                  &
     &  lbc_halo_x,                                                     &
     &  lbc_halo_y,                                                     &
     &  rimwidth,                                                       &
     &  LBC,                                                            &
     &  coeff1,                                                         &
     &  coeff2,                                                         &
     &  Lbc_size,                                                       &
     &  level_type,                                                     &
! VarRes grid info in degrees
     &  L_var_lbc,                                                      &
     &  Lambda_in,                                                      &
     &  Phi_in,                                                         &
! hi - indexes and weights
     &  l_calc_lbc_wts,                                                 &
     &  lbc_index_bl  , lbc_index_br  ,                                 &
     &  lbc_weights_tr, lbc_weights_br,                                 &
     &  lbc_weights_bl, lbc_weights_tl,                                 &
! vi
     &  max_seg_size,                                                   &
     &  N_SEGS,                                                         &
     &  MAX_LEVELS_PER_PE,                                              &
     &  GATHER_PE                                                       &
     &, L_VI                                                            &
! vi - src
     &,       src_model_levels                                          &
     &,       src_ht_gen_method                                         &
     &,       src_first_r_rho                                           &
     &,       src_z_top_model                                           &
     &,       src_eta_theta                                             &
     &,       src_eta_rho                                               &
     &,       src_first_level                                           &
     &,       src_last_level                                            &
! vi - lbc
     &,       lbc_model_levels                                          &
     &,       lbc_ht_gen_method                                         &
     &,       lbc_first_r_rho                                           &
     &,       lbc_z_top_model                                           &
     &,       lbc_eta_theta                                             &
     &,       lbc_eta_rho                                               &
     &,       lbc_first_level                                           &
     &,       lbc_last_level                                            &
     &,       i_uv                                                      &
     & )

      Implicit None

!
! Description:
!   Generates the Atmosphere LBCs from the model prognostics. Performs
!   horizontal and vertical interpolation as required.
!
! Method:
!   1. Horizontal Interpolation. (HI)
!      HI from source grid to lbc grid is done by level on a PE.
!      a) Determine no of levels to be interpolated on each PE.
!      b) Gather global source data on PEs from local data.
!      c) Determine the horizontal interpolation coefficients.
!      d) Do the HI through suboutine H_INT_BL.
!      e) If no vertical interpolation, gather all lbc data on PE 0.
!
!   2. Vertical interpolation (VI) (to follow)
!
! Current Code Owner: Dave Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables :

! Subroutine arguments

      INTEGER                                                           &
     &  LOCAL_ROW_LENGTH                                                &
                      ! IN  : East-West size of source field
     &, LOCAL_ROWS                                                      &
                      ! IN  : North-South size of source field
     &, GLOBAL_ROW_LENGTH                                               &
                      ! IN  : East-West Size of full gathered level
                      !     : glsize(1,fld_type)
     &, GLOBAL_ROWS                                                     &
                      ! IN  : North-South Size of full gathered level
                      !     : glsize(2,fld_type)
     &, SOURCE_HALO_X                                                   &
                      ! IN  : Size of East-West halo for source field
     &, SOURCE_HALO_Y                                                   &
                      ! IN  : Size of North-South halo for source field

     &, SOURCE_LEVELS                                                   &
                      ! IN  : Number of levels for source field
     &, SOURCE_FLD_TYPE                                                 &
                        ! IN  : Field type of source field
     &, SOURCE_HALO_TYPE ! IN  : Halo type of source field

      Logical                                                           &
     &  SOURCE_CYCLIC                                                   &
                      ! IN   : T if Source Grid is cyclic
     &, SOURCE_ROTATED! IN   : T if Source Grid is rotated

       Logical                                                          &
     &  L_CALC_LBC_WTS ! IN  : T if Interp coeffs to be calculated

      INTEGER                                                           &
     &  OROG_LOCAL_ROW_LENGTH                                           &
                      ! IN  : East-West size of source field
     &, OROG_LOCAL_ROWS                                                 &
                      ! IN  : North-South size of source field
     &, OROG_GLOBAL_ROW_LENGTH                                          &
                      ! IN  : East-West Size of full gathered level
                      !     : glsize(1,fld_type)
     &, OROG_GLOBAL_ROWS                                                &
                      ! IN  : North-South Size of full gathered level
                      !     : glsize(2,fld_type)
     &, OROG_HALO_X                                                     &
                    ! IN  : Size of East-West halo for source field
     &, OROG_HALO_Y                                                     &
                    ! IN  : Size of North-South halo for source field
     &, OROG_HALO_TYPE                                                  &
     &, OROG_FLD_TYPE

      Real                                                              &
     &  source_delta_lat                                                &
     &, source_delta_long                                               &
     &, source_first_lat                                                &
     &, source_first_long                                               &
     &, source_pole_lat                                                 &
     &, source_pole_long                                                &
     &, lbc_delta_lat                                                   &
     &, lbc_delta_long                                                  &
     &, lbc_first_lat                                                   &
     &, lbc_first_long                                                  &
     &, lbc_pole_lat                                                    &
     &, lbc_pole_long

! Input VarRes grid info in degrees      
      Logical L_var_lbc 
      REAL                                                              &
     &  Lambda_in ( 1-lbc_halo_x: lbc_row_len+lbc_halo_x )              & 
     &, Phi_in ( 1-lbc_halo_y: lbc_rows+lbc_halo_y )
      
      Real                                                              &
     &  orog_first_lat                                                  &
     &, orog_first_long

      Integer                                                           &
     &  rimwidth                                                        &
     &, lbc_halo_x                                                      &
     &, lbc_halo_y                                                      &
     &, LBC_SIZE                                                        &
                      ! IN  : Size of single level of LBC array
     &, LBC_LEVELS                                                      &
                      ! IN  : Number of levels of LBC data
     &, level_type                                                      &
     &, max_seg_size                                                    &
                      ! IN  : Size of data on each processor for vert.
                      !       int. step. Something like:
                      ! MAX(minimum_segment_size,((LBC_SIZE-1)/nproc)+1)
     &, N_SEGS                                                          &
                      ! IN  : Number of segments,
                      ! given by ((LBC_SIZE-1)/LBC_SEG_SIZE)+1
     &, MAX_LEVELS_PER_PE                                               &
                      ! IN  : Maximum no of levels per PE - basically:
                      !       MAX_LEVELS_PER_PE=(SOURCE_LEVELS-1)/NPROC
     &, GATHER_PE                                                       &
                      ! IN  : Processor that will contain the final
                      !       3D LBC field
                      !       (probably usually PE 0)
     &, i_uv          ! IN  : i_uv=1 => u or v field

      Real                                                              &
     & Source_Field(1-SOURCE_HALO_X:LOCAL_ROW_LENGTH+SOURCE_HALO_X,     &
     &              1-SOURCE_HALO_Y:LOCAL_ROWS+SOURCE_HALO_Y,           &
     &              SOURCE_LEVELS)

      Real                                                              &
     & Orog_Field(1-OROG_HALO_X:OROG_LOCAL_ROW_LENGTH+OROG_HALO_X,      &
     &            1-OROG_HALO_Y:OROG_LOCAL_ROWS+OROG_HALO_Y )

      Real                                                              &
     &  Coeff1 (lbc_size)                                               &
                            ! \ Coefficients for rotating winds
     &, Coeff2 (lbc_size)   ! / on lbc points


      Real lbc (lbc_size, lbc_levels)  !  OUT : LBC data

! vi
      LOGICAL                                                           &
     & L_VI           ! IN  : Perform vertical inerpolation?

! vi - src
      Integer   :: src_model_levels
      Integer   :: src_first_level
      Integer   :: src_last_level
      Integer   :: src_ht_gen_method
      Integer   :: src_first_r_rho
      Real      :: src_z_top_model
      Real      :: src_eta_theta (0:src_model_levels)
      Real      :: src_eta_rho   (  src_model_levels)

! vi - lbc
      Integer   :: lbc_model_levels
      Integer   :: lbc_first_level
      Integer   :: lbc_last_level
      Integer   :: lbc_ht_gen_method
      Integer   :: lbc_first_r_rho
      Real      :: lbc_z_top_model
      Real      :: lbc_eta_theta (0:lbc_model_levels)
      Real      :: lbc_eta_rho   (  lbc_model_levels)

! Local parameters:

      Character (Len=*), Parameter :: RoutineName= 'Make_LBCs'

! Local scalars:

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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! GC - General Communication primitives package. For use on
! multiprocessor shared memory and message passing systems.
!
!
! LICENSING TERMS
!
!  GC is provided free of charge. Unless otherwise agreed with SINTEF,
!  use and redistribution in source and binary forms are permitted
!  provided that
!
!      (1) source distributions retain all comments appearing within
!          this file header, and
!
!      (2) distributions including binaries display the following
!          acknowledgement:
!
!              "This product includes software developed by SINTEF.",
!
!          in the documentation or other materials provided with the
!          distribution and in all advertising materials mentioning
!          features or use of this software.
!
!  The name of SINTEF may not be used to endorse or promote products
!  derived from this software without specific prior written
!  permission.  SINTEF disclaims any warranty that this software will
!  be fit for any specific purposes. In no event shall SINTEF be liable
!  for any loss of performance or for indirect or consequential damage
!  or direct or indirect injury of any kind. In no case shall SINTEF
!  be liable for any representation or warranty make to any third party
!  by the users of this software.
!
!
! Fortran header file. PLEASE use the parameter variables in user
! routines calling GC and NOT the numeric values. The latter are
! subject to change without further notice.
!
!---------------------------------------------- ------------------------
! $Id: gps0h501,v 1.6 2000/04/17 10:05:47 t11ps Exp $
! (C) Jorn Amundsen, Roar Skaalin, SINTEF Industrial Mathematics.

!    4.4   30/09/97  Added code to permit the SHMEM/NAM timeout
!                    value to be set from a shell variable.
!                      Author: Bob Carruthers  Cray Research.
!    5.1   17/04/00  Fixed/Free format. P.Selwood.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!     GC general options
      INTEGER, PARAMETER :: GC_OK         =     0
      INTEGER, PARAMETER :: GC_FAIL       =    -1
      INTEGER, PARAMETER :: GC_NONE       =     0
      INTEGER, PARAMETER :: GC_ANY        =    -1
      INTEGER, PARAMETER :: GC_DONTCARE   =    -1
      INTEGER, PARAMETER :: GC_SHM_DIR    =     1
      INTEGER, PARAMETER :: GC_SHM_SAFE   =     2
      INTEGER, PARAMETER :: GC_NAM_TIMEOUT=     4
      INTEGER, PARAMETER :: GC_SHM_GET    = -9999
      INTEGER, PARAMETER :: GC_SHM_PUT    = -9998
      INTEGER, PARAMETER :: GC_USE_GET    = -9999
      INTEGER, PARAMETER :: GC_USE_PUT    = -9998

!     GC functions
      INTEGER GC_COMLEN, GC_ISIZE, GC_RSIZE, GC_ME, GC_NPROC

!     GC groups (GCG) support
      INTEGER, PARAMETER :: GC_ALLGROUP = 0
      INTEGER, PARAMETER :: GCG_ALL = GC_ALLGROUP

!     GC groups (GCG) functions
      INTEGER GCG_ME

!     GC reserved message tags
      INTEGER, PARAMETER :: GC_MTAG_LOW   = 999999901
      INTEGER, PARAMETER :: GC_MTAG_HIGH  = 999999999

!     GCG_RALLETOALLE index parameters
      INTEGER, PARAMETER :: S_DESTINATION_PE = 1
      INTEGER, PARAMETER :: S_BASE_ADDRESS_IN_SEND_ARRAY = 2
      INTEGER, PARAMETER :: S_NUMBER_OF_ELEMENTS_IN_ITEM = 3
      INTEGER, PARAMETER :: S_STRIDE_IN_SEND_ARRAY = 4
      INTEGER, PARAMETER :: S_ELEMENT_LENGTH = 5
      INTEGER, PARAMETER :: S_BASE_ADDRESS_IN_RECV_ARRAY = 6
      INTEGER, PARAMETER :: S_STRIDE_IN_RECV_ARRAY = 7

      INTEGER, PARAMETER :: R_SOURCE_PE = 1
      INTEGER, PARAMETER :: R_BASE_ADDRESS_IN_RECV_ARRAY = 2
      INTEGER, PARAMETER :: R_NUMBER_OF_ELEMENTS_IN_ITEM = 3
      INTEGER, PARAMETER :: R_STRIDE_IN_RECV_ARRAY = 4
      INTEGER, PARAMETER :: R_ELEMENT_LENGTH = 5
      INTEGER, PARAMETER :: R_BASE_ADDRESS_IN_SEND_ARRAY = 6
      INTEGER, PARAMETER :: R_STRIDE_IN_SEND_ARRAY = 7

! vi
      INTEGER                                                           &
     &  send_map(7,N_SEGS*MAX_LEVELS_PER_PE)                            &
     &, recv_map(7,N_SEGS*SOURCE_LEVELS)

      INTEGER                                                           &
     &  map(SOURCE_LEVELS)                                              &
                           ! mapping of processors holding which level
     &, LEVEL_ON_GAT_PE(SOURCE_LEVELS)                                  &
                                       ! what is the level number on the
                           !            gathering PE of this level
     &, LEVEL_COUNT(0:nproc-1)  ! Count of levels on each PE

      INTEGER                                                           &
     &  iproc                                                           &
                    ! loop counter over processors
     &, iseg                                                            &
                    ! loop counter over segments
     &, seg_start                                                       &
                    ! first point of segment
     &, seg_end                                                         &
                    ! last point of segment
     &, seg_size_pe                                                     &
                    ! local segment size
     &, seg_size(n_segs)                                                &
     &, n_send                                                          &
                    ! number of items of data to send
     &, n_recv                                                          &
                    ! number of items of data to receive
     &, k                                                               &
                    ! loop counter over levels
     &, info                                                            &
                    ! return code
     &, flag                                                            &
                    ! GCOM input code
     &, ErrorStatus ! Return Code

      Character (Len=80)  :: CMessage

! Local dynamic arrays:

! Whole levels of the source field gathered onto processors

      Real, dimension (:,:,:), allocatable :: Gathered_Source_Field
      Real, dimension (  :,:), allocatable :: Gathered_Orog_Field

! LBC data after horizontal interpolation

      Real, dimension (  :,:), allocatable :: lbc_data
      Real, dimension (    :), allocatable :: lbc_orog

! Indexes for bottom left & right points in bi-linear interpolation

      Integer, dimension (lbc_size) :: lbc_index_bl
      Integer, dimension (lbc_size) :: lbc_index_br

! Weights for 4 points in bi-linear interpolation

      Real,    dimension (lbc_size) :: lbc_weights_tr
      Real,    dimension (lbc_size) :: lbc_weights_br
      Real,    dimension (lbc_size) :: lbc_weights_tl
      Real,    dimension (lbc_size) :: lbc_weights_bl

! LBC data before/after vertical interpolation

      Real,    dimension (:,:), allocatable :: lbc_vi_data_in
      Real,    dimension (:,:), allocatable :: lbc_vi_data_out
      Real,    dimension (  :), allocatable :: lbc_vi_orog

! Temp

      Integer kk,ipt  ! temp for write statement.
      integer halo_x, halo_y, lbc_row_len, lbc_rows
      integer row,pt,level
      logical lprint

! Function & Subroutine calls:
      External Gather_Field, H_Int_BL, LBC_Interp_Coeffs
      EXTERNAL :: Gather_Field_ML

!- End of header
!--------------------------------------------------------------------

! 1.0 Set up arrays to map which levels are on which processors

      Do iproc=0,nproc-1
        level_count(iproc)=0
      Enddo

      Do k=1,source_levels
        map(k)              = MOD(k-1,nproc)
        level_count(map(k)) = level_count(map(k))+1
        level_on_gat_pe(k)  = level_count(map(k))
      Enddo

! 2.0 Gather levels to different PEs

      allocate ( Gathered_Source_Field ( global_row_length,             &
     &                                   global_rows,                   &
     &                                   max_levels_per_pe ) )

! DEPENDS ON: gather_field_ml
      Call Gather_Field_ML (                                            &
     &  Source_Field, Gathered_Source_Field,                            &
     &  Local_Row_Length + 2*Source_Halo_x,                             &
     &  Local_Rows + 2*Source_Halo_y,                                   &
     &  Source_Levels,                                                  &
     &  Global_Row_Length, Global_Rows, Max_Levels_per_pe,              &
     &  Map, Level_on_Gat_pe, Source_Fld_Type, Source_Halo_Type)


      If ( l_vi ) Then

! Set up a global field of orography on PE 0

        allocate ( Gathered_Orog_Field ( orog_global_row_length,        &
     &                                   orog_global_rows ) )

! DEPENDS ON: gather_field
        Call Gather_Field (Orog_Field,                                  &
     &                     Gathered_Orog_Field,                         &
     &                     Orog_Local_Row_Length, Orog_Local_Rows,      &
     &                     Orog_Global_Row_Length, Orog_Global_Rows,    &
     &                     Orog_Fld_Type, Orog_Halo_Type,               &
     &                     0, gc_all_proc_group, info, cmessage )

! 3.0 Perform horizontal interpolation.
!     This processor has LEVEL_COUNT(mype) levels to do

! 3.1 Allocate workspace for the interpolation coefficients

!      work space allocated in GEN_INTF_A

! 3.2 Allocate work space to horizontally interpolated orography

        allocate ( lbc_orog (lbc_size) )

        If (mype == 0) Then

! 3.3 Calculate the interpolation coefficients for orography field

! NB. LBC_Interp_Coeffs must be calculated for the orography before
!     the prognostic so that on exit COEFF1 and COEFF2 correspond to
!     the prognostic.

        If (l_calc_lbc_wts ) Then

! DEPENDS ON: lbc_interp_coeffs
          Call LBC_interp_coeffs (                                      &
     &         lbc_size                                                 &
     &,        orog_global_row_length                                   &
     &,        orog_global_rows                                         &
     &,        source_delta_lat                                         &
     &,        source_delta_long                                        &
     &,        orog_first_lat                                           &
     &,        orog_first_long                                          &
     &,        source_pole_lat                                          &
     &,        source_pole_long                                         &
     &,        source_cyclic                                            &
     &,        source_rotated                                           &
     &,        lbc_row_len                                              &
     &,        lbc_rows                                                 &
     &,        L_var_lbc                                                &
     &,        Lambda_in                                                &
     &,        Phi_in                                                   &
     &,        lbc_delta_lat                                            &
     &,        lbc_delta_long                                           &
     &,        lbc_first_lat                                            &
     &,        lbc_first_long                                           &
     &,        lbc_pole_lat                                             &
     &,        lbc_pole_long                                            &
     &,        rimwidth                                                 &
     &,        lbc_halo_x                                               &
     &,        lbc_halo_y                                               &
     &,        lbc_index_bl                                             &
     &,        lbc_index_br                                             &
     &,        lbc_weights_tr                                           &
     &,        lbc_weights_br                                           &
     &,        lbc_weights_bl                                           &
     &,        lbc_weights_tl                                           &
     &,        coeff1                                                   &
     &,        coeff2                                                   &
     &,        i_uv                                                     &
     & )

        End If  !  l_calc_lbc_wts

! 3.4 Do the bi-linear horizontal interpolation for orography

! DEPENDS ON: h_int_bl
          call h_int_bl (                                               &
     &         orog_global_rows                                         &
     &,        orog_global_row_length                                   &
     &,        lbc_size                                                 &
     &,        lbc_index_bl                                             &
     &,        lbc_index_br                                             &
     &,        gathered_orog_field                                      &
     &,        lbc_weights_bl                                           &
     &,        lbc_weights_br                                           &
     &,        lbc_weights_tl                                           &
     &,        lbc_weights_tr                                           &
     &,        lbc_orog                                                 &
     & )

          deallocate ( gathered_orog_field)

        End If  !  If mype=0

      End If

! 3.2 Calculate the interpolation coefficients for the prognostic

      If ( l_calc_lbc_wts) Then

! DEPENDS ON: lbc_interp_coeffs
      Call LBC_interp_coeffs (                                          &
     &     lbc_size                                                     &
     &,    global_row_length                                            &
     &,    global_rows                                                  &
     &,    source_delta_lat                                             &
     &,    source_delta_long                                            &
     &,    source_first_lat                                             &
     &,    source_first_long                                            &
     &,    source_pole_lat                                              &
     &,    source_pole_long                                             &
     &,    source_cyclic                                                &
     &,    source_rotated                                               &
     &,    lbc_row_len                                                  &
     &,    lbc_rows                                                     &
     &,    L_var_lbc                                                    &
     &,    Lambda_in                                                    &
     &,    Phi_in                                                       &
     &,    lbc_delta_lat                                                &
     &,    lbc_delta_long                                               &
     &,    lbc_first_lat                                                &
     &,    lbc_first_long                                               &
     &,    lbc_pole_lat                                                 &
     &,    lbc_pole_long                                                &
     &,    rimwidth                                                     &
     &,    lbc_halo_x                                                   &
     &,    lbc_halo_y                                                   &
     &,    lbc_index_bl                                                 &
     &,    lbc_index_br                                                 &
     &,    lbc_weights_tr                                               &
     &,    lbc_weights_br                                               &
     &,    lbc_weights_bl                                               &
     &,    lbc_weights_tl                                               &
     &,    coeff1                                                       &
     &,    coeff2                                                       &
     &,    i_uv                                                         &
     & )

      End If   !  l_calc_lbc_wts

! 3.3 Allocate work space for the horizontally interpolated data

      allocate ( lbc_data (lbc_size, max_levels_per_pe) )

! 3.4 Do the bi-linear horizontal interpolation

      Do level = 1, level_count(mype)  !  Levels to do on this pe.

! DEPENDS ON: h_int_bl
        call h_int_bl (                                                 &
     &       global_rows                                                &
                                   !  glsize(2)
     &,      global_row_length                                          &
                                   !  glsize(1)
     &,      lbc_size                                                   &
     &,      lbc_index_bl                                               &
     &,      lbc_index_br                                               &
     &,      gathered_source_field (1,1,level)                          &
     &,      lbc_weights_bl                                             &
     &,      lbc_weights_br                                             &
     &,      lbc_weights_tl                                             &
     &,      lbc_weights_tr                                             &
     &,      lbc_data (1,level)                                         &
     & )

      Enddo  !  Loop over levels

! 3.5 Deallocate space reserved for interpolation coefficients

!      work space deallocated in GEN_INTF_A

!     Deallocate space reserved for gathered source data

      deallocate ( gathered_source_field )

! We should now have an array LBC_DATA(:,:) where we have filled in
! 1:LEVEL_COUNT(mype) levels of data (this may be zero if there are
! more processors than levels !)

! 4.0 Perform Vertical Interpolation.

      If (L_VI) Then

!       Determine segment size for each segment.
!       Max_Seg_Size is maximum size of any segment.

        Do iseg = 1, n_segs
          seg_start=((iseg-1)*Max_Seg_Size)+1
          seg_end=MIN(LBC_SIZE,seg_start+Max_Seg_Size-1)
          seg_size(iseg) = seg_end - seg_start + 1
        End Do

!       Get segment size for this pe.
        If (mype < n_segs) Then
          seg_size_pe = seg_size(mype+1)
        Else
          seg_size_pe = 1
        End If

        ! Move data so that each processor has a segment of data, on
        ! all levels

        n_send=0    ! number of messages I send
        n_recv=0    ! number of messages I receive
        send_map (:,:) = 0
        recv_map (:,:) = 0

        DO iseg=1,N_SEGS

          iproc=iseg-1  ! This is the processor doing vertical
!                       ! interpolation on segment iseg

          seg_start=((iseg-1)*Max_Seg_Size)+1
          seg_end=MIN(LBC_SIZE,seg_start+Max_Seg_Size-1)

          DO k=1,SOURCE_LEVELS

            IF (MAP(k) == mype) THEN  ! I've got data to send
              n_send=n_send+1
              send_map(S_DESTINATION_PE,n_send)=iproc
              send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_send)=            &
     &          seg_start+((LEVEL_ON_GAT_PE(k)-1)*LBC_SIZE)
              send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_send)=1
              send_map(S_STRIDE_IN_SEND_ARRAY,n_send)=1
              send_map(S_ELEMENT_LENGTH,n_send)=seg_end-seg_start+1
              send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_send)=            &
     &          1+(k-1)*seg_size(iseg)
              send_map(S_STRIDE_IN_RECV_ARRAY,n_send)=1
            ENDIF

            IF (iproc == mype) Then  ! I am receiving data
              n_recv=n_recv+1
              recv_map(R_SOURCE_PE,n_recv)=MAP(k)
              recv_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recv)=            &
     &          1+(k-1)*seg_size(iseg)

              recv_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recv)=1
              recv_map(R_STRIDE_IN_RECV_ARRAY,n_recv)=1
              recv_map(R_ELEMENT_LENGTH,n_recv)=seg_end-seg_start+1
              recv_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recv)=            &
     &          seg_start+((LEVEL_ON_GAT_PE(k)-1)*LBC_SIZE)
              recv_map(R_STRIDE_IN_SEND_ARRAY,n_recv)=1
            ENDIF

          ENDDO ! k
        ENDDO ! iseg

        allocate ( lbc_vi_data_in (seg_size_pe, source_levels) )

        flag=GC_NONE
        info=GC_NONE

        CALL GCG_RALLTOALLE(LBC_DATA,send_map,n_send,                   &
     &                      LBC_SIZE*MAX_LEVELS_PER_PE,                 &
     &                      LBC_VI_DATA_IN,recv_map,n_recv,             &
     &                      SEG_SIZE_PE*SOURCE_LEVELS,                  &
     &                      gc_all_proc_group,flag,info)

        If (info /= GC_NONE) Then
          ErrorStatus = 10
          Write (cmessage,*) 'GCG_RALLTOALLE failed'
! DEPENDS ON: ereport
          Call Ereport ( RoutineName, ErrorStatus, cmessage)
        End If

        ! Now the prognostic is in LBC_VI_DATA_IN(:,SOURCE_LEVELS)

! We now need to scatter the orography to the same number of segments

        n_send=0    ! number of messages I send
        n_recv=0    ! number of messages I receive
        send_map (:,:) = 0
        recv_map (:,:) = 0

        Do iseg=1,N_SEGS

          iproc=iseg-1  ! This is the processor doing vertical
!                       ! interpolation on segment iseg

          seg_start=((iseg-1)*Max_Seg_Size)+1
          seg_end=MIN(LBC_SIZE,seg_start+Max_Seg_Size-1)

            IF ( mype == 0 ) Then  ! Only PE 0 has orog to send out
              n_send=n_send+1
              send_map(S_DESTINATION_PE,n_send)=iproc
              send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_send)= seg_start
              send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_send)=1
              send_map(S_STRIDE_IN_SEND_ARRAY,n_send)=1
              send_map(S_ELEMENT_LENGTH,n_send)=seg_end-seg_start+1
              send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_send)=1
              send_map(S_STRIDE_IN_RECV_ARRAY,n_send)=1
            ENDIF

            IF (iproc == mype) Then  ! This PE to receive data
              n_recv=n_recv+1
              recv_map(R_SOURCE_PE,n_recv)=0
              recv_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recv)=1
              recv_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recv)=1
              recv_map(R_STRIDE_IN_RECV_ARRAY,n_recv)=1
              recv_map(R_ELEMENT_LENGTH,n_recv)=seg_end-seg_start+1
              recv_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recv)=seg_start
              recv_map(R_STRIDE_IN_SEND_ARRAY,n_recv)=1
            ENDIF

        End Do ! iseg

        allocate ( lbc_vi_orog (seg_size_pe) )

        flag=GC_NONE
        info=GC_NONE

        CALL GCG_RALLTOALLE(LBC_OROG,send_map,n_send,                   &
     &                      LBC_SIZE,                                   &
     &                      LBC_VI_OROG,recv_map,n_recv,                &
     &                      seg_size_pe,                                &
     &                      gc_all_proc_group,flag,info)

        If (info /= GC_NONE) Then
          ErrorStatus = 10
          Write (cmessage,*) 'GCG_RALLTOALLE failed'
! DEPENDS ON: ereport
          Call Ereport ( RoutineName, ErrorStatus, cmessage)
        End If

        ! Now the orography is in LBC_VI_OROG(:)

        deallocate (lbc_orog)
        allocate ( lbc_vi_data_out (seg_size_pe, lbc_levels) )

        If (mype < N_SEGS) Then

          ! I have to do some vertical interpolation...

! DEPENDS ON: lbc_vert_interp
          Call LBC_Vert_Interp (                                        &
     &         LBC_VI_Data_in,                                          &
     &         LBC_VI_Orog,                                             &
     &         LBC_VI_Data_out,                                         &
     &         seg_size_pe,                                             &
     &         level_type                                               &
! src
     &,        src_model_levels                                         &
     &,        source_levels                                            &
     &,        src_first_level                                          &
     &,        src_last_level                                           &
     &,        src_ht_gen_method                                        &
     &,        src_first_r_rho                                          &
     &,        src_z_top_model                                          &
     &,        src_eta_theta                                            &
     &,        src_eta_rho                                              &
! lbc
     &,        lbc_model_levels                                         &
     &,        lbc_levels                                               &
     &,        lbc_first_level                                          &
     &,        lbc_last_level                                           &
     &,        lbc_ht_gen_method                                        &
     &,        lbc_first_r_rho                                          &
     &,        lbc_z_top_model                                          &
     &,        lbc_eta_theta                                            &
     &,        lbc_eta_rho                                              &
     & )

        ENDIF

        ! Now the data is in LBC_VI_DATA_OUT(:,LBC_LEVELS)
        ! This needs to be moved to LBC on processor GATHER_PE

        n_send=0
        n_recv=0
        send_map (:,:) = 0
        recv_map (:,:) = 0

        DO iseg=1,N_SEGS

          iproc=iseg-1  ! This is the processor doing vertical
!                       ! interpolation on segment iseg

          seg_start=((iseg-1)*Max_Seg_Size)+1
          seg_end=MIN(LBC_SIZE,seg_start+Max_Seg_Size-1)

          IF (mype == iproc) Then ! I've got data to send
            n_send=n_send+1
            send_map(S_DESTINATION_PE,n_send)=GATHER_PE
            send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_send)=1
            send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_send)=LBC_LEVELS
            send_map(S_STRIDE_IN_SEND_ARRAY,n_send)=SEG_SIZE(iseg)
            send_map(S_ELEMENT_LENGTH,n_send)=seg_end-seg_start+1
            send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_send)=seg_start
            send_map(S_STRIDE_IN_RECV_ARRAY,n_send)=LBC_SIZE
          ENDIF

          IF (mype == GATHER_PE) Then ! I've got data to receive
            n_recv=n_recv+1
            recv_map(R_SOURCE_PE,n_recv)=iproc
            recv_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recv)=seg_start
            recv_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recv)=LBC_LEVELS
            recv_map(R_STRIDE_IN_RECV_ARRAY,n_recv)=LBC_SIZE
            recv_map(R_ELEMENT_LENGTH,n_recv)=seg_end-seg_start+1
            recv_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recv)=1
            recv_map(R_STRIDE_IN_SEND_ARRAY,n_recv)=SEG_SIZE(iseg)
          ENDIF

        ENDDO ! iseg

        flag=GC_NONE
        info=GC_NONE

        CALL GCG_RALLTOALLE(LBC_VI_DATA_OUT,send_map,n_send,            &
     &                      SEG_SIZE_PE*LBC_LEVELS,                     &
     &                      LBC,recv_map,n_recv,                        &
     &                      LBC_SIZE*LBC_LEVELS,                        &
     &                      gc_all_proc_group,flag,info)

        If (info /= GC_NONE) Then
          ErrorStatus = 20
          Write (cmessage,*) 'GCG_RALLTOALLE failed'
! DEPENDS ON: ereport
          Call Ereport ( RoutineName, ErrorStatus, cmessage)
        End If

! deallocate arrays
        deallocate ( lbc_vi_data_in  )
        deallocate ( lbc_vi_data_out )
        deallocate ( lbc_vi_orog     )

      Else ! .NOT. L_VI - ie. no vertical interpolation

        ! Need to move the data from LBC_DATA with different levels
        ! on different PEs to LBC with all the data on GATHER_PE

        n_send=0
        n_recv=0
        send_map (:,:) = 0
        recv_map (:,:) = 0

        Do k=1,lbc_levels

          If (Map(k) == mype) Then ! I need to send data
            n_send=n_send+1
            send_map(S_DESTINATION_PE,n_send)=GATHER_PE
            send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_send)=              &
     &        1+((LEVEL_ON_GAT_PE(k)-1)*LBC_SIZE)
            send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_send)=1
            send_map(S_STRIDE_IN_SEND_ARRAY,n_send)=1
            send_map(S_ELEMENT_LENGTH,n_send)=LBC_SIZE
            send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_send)=              &
     &        1+(k-1)*LBC_SIZE
            send_map(S_STRIDE_IN_RECV_ARRAY,n_send)=1
          Endif

          If (mype == GATHER_PE) Then ! I need to receive data
            n_recv=n_recv+1
            recv_map(R_SOURCE_PE,n_recv)=MAP(k)
            recv_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recv)=              &
     &        1+(k-1)*LBC_SIZE
            recv_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recv)=1
            recv_map(R_STRIDE_IN_RECV_ARRAY,n_recv)=1
            recv_map(R_ELEMENT_LENGTH,n_recv)=LBC_SIZE
            recv_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recv)=              &
     &        1+((LEVEL_ON_GAT_PE(k)-1)*LBC_SIZE)
            recv_map(R_STRIDE_IN_SEND_ARRAY,n_recv)=1
          Endif

        Enddo ! k

        flag=GC_NONE
        info=GC_NONE

        Call GCG_RALLTOALLE(LBC_DATA,send_map,n_send,                   &
     &                      LBC_SIZE*MAX_LEVELS_PER_PE,                 &
     &                      LBC,recv_map,n_recv,                        &
     &                      LBC_SIZE*LBC_LEVELS,                        &
     &                      gc_all_proc_group,flag,info)

        If (info /= GC_NONE) THEN
          ErrorStatus = 30
          Write (cmessage,*) 'GCG_RALLTOALLE failed'
! DEPENDS ON: ereport
          Call Ereport ( RoutineName, ErrorStatus, cmessage)
        Endif

      Endif ! IF (L_VI)

      deallocate ( lbc_data )

      ! The array LBC on processor GATHER_PE should now contain the full
      ! 3D LBC array

      Return
      END SUBROUTINE MAKE_LBCS
