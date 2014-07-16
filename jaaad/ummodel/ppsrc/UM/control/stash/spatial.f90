
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: SPATIAL --------------------------------------------------
!LL
!LL  Purpose: Performs general spatial processing on an input field to
!LL           produce an output field or scalar within STASH.  Lower-
!LL           level routines are called to perform the various spatial
!LL           processing options.
!LL
!LL  Author:   T.Johns/S.Tett
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.3  30/03/94  Explicitly declare (sub-addressed) output field
!LL                  fieldout using 'lenout' dimension.  Tim Johns.
!LL   3.3  16/09/93  Pass LOGICAL lmasswt to processing routines to
!LL                  denote that level-by-level mass-weights exist.
!LL   4.3   9/12/96  Added mpp code.
!LL                  Moved calculation of weighting and masking terms
!LL                  up from processing routines.            P.Burton
!LL   4.4   13/06/97 mpp : Where reduction spatial meaning takes place
!LL                  processors not getting results should set
!LL                  their diagnostic space to zeros.          P.Burton
!LL   4.4   22/10/97 mpp : Prevent uninitialised points when
!LL                  pstar_weight on U or C grid S.D.Mullerworth
!LL   5.0   13/07/99 Changes for C-P C grid upgrade.
!LL                  R Rawlins.
!LL   5.0   22/06/99 Added halo_type argument            P.Burton
!LL   5.0   01/11/99 Added halo size arguments
!LL                  Add halos to weighting/mask arrays  P.Burton
!LL   5.0   17/11/99 Added no_rows as an argument, and use to
!LL                  dimension local weighting/masking arrays
!LL                                                      P.Burton
!LL   5.0    9/11/99 Correct to south-north order for ystart,yend
!LL                                                       P.Burton
!LL   5.1   15/05/0 S-N ordering consistency correction. R Rawlins
!     5.3   05/10/01 Correction to mes diagnostics for fields with halos
!                    ie prognostics u,v,etc. R Rawlins
!     5.4   02/09/02 Use sea mask for grid type=3. This is not the
!                    reverse of the land mask when coastal tiling
!                    is used.                            K.Williams
!     5.5   17/02/03 Amendment for Wave model. D.Holmes-Bell
!     6.0   06/10/03 Cater for river grid 23. C.Bunton
!     6.1   16/08/04 Increase first dimension in certain arrays
!                    that are passed to STEXTC. Jens-Olaf Beismann
!                    NEC.   Anthony A. Dickinson
!     6.2   18/10/05 Fix bugs introduced by gan1f601      Andy Malcolm
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D71
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!LL-----------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
      SUBROUTINE SPATIAL(fieldin,vx,vy,vz,GR,st_grid,                   &
     &                   fld_type,halo_type,                            &
     &                   halo_x,halo_y,                                 &
     &                   lcyclic,lmasswt,                               &
     &      n_cols_out,n_rows_out,                                      &
     &      this_index_lev,level_list,index_lev,no_of_levels,           &
     &      no_of_levels_masswt,                                        &
     &      p,pstar,                                                    &
     &      cos_v_latitude,cos_theta_latitude,land,sea,                 &
     &      row_length,rows,n_rows,no_rows,model_levels,                &
     &      fieldout,lenout,                                            &
     &      control,control_size,rmdi,                                  &
     &      icode,cmessage)
!
      IMPLICIT NONE
!
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

      INTEGER                                                           &
     &    vx,vy,vz,                                                     &
                                            ! IN size of fieldin
     &    lenout,                                                       &
                                            ! IN size of fieldout
     &    GR,                                                           &
                                            ! IN ppxref gridtype code
     &    st_grid,                                                      &
                                            ! IN STASH gridtype code
     &    fld_type,                                                     &
                                            ! IN field type (u/v/p)
     &    halo_type,                                                    &
                                            ! IN halo type
     &    halo_x,                                                       &
                                            ! IN EW halo of input
     &    halo_y,                                                       &
                                            ! IN NS halo of input
     &    n_rows_out,                                                   &
                                            ! OUT no. of output rows
     &    n_cols_out,                                                   &
                                            ! OUT no. of output cols
     &    this_index_lev,                                               &
                                            ! IN level index, this field
     &    row_length,rows,n_rows,                                       &
                                            ! IN horiz. sizes (C grid)
     &    no_rows,                                                      &
                                            ! IN number of rows used
     &    model_levels,                                                 &
                                            ! IN vertical size
     &    control_size,                                                 &
                                            ! IN size of control record
     &    control(control_size),                                        &
                                            ! IN control record
     &    icode,                                                        &
                                            ! OUT error code 0 if ok
     &    no_of_levels,                                                 &
                                            ! IN no of levels
     &    no_of_levels_masswt,                                          &
                                ! IN levels for mass weighting array
                                ! lmasswt F: =1; T: =no_of_levels
     &    index_lev(no_of_levels),                                      &
                                            ! IN index to levels
     &    level_list(no_of_levels)          ! IN model level list
      REAL                                                              &
     &    fieldin(vx,vy,vz),                                            &
                             ! IN fieldin which is to be acted on
     &    p(1-offx:row_length+offx,1-offy:rows+offy,model_levels),      &
                                            ! IN pressure (rho levels)
     &    pstar(row_length,rows),                                       &
                                            ! IN surface pressure
     &    cos_v_latitude(row_length,n_rows),                            &
                                            ! IN v-grid area fn
     &    cos_theta_latitude(row_length,rows),                          &
                                                ! IN T-grid area fn
     &    fieldout(lenout),                                             &
                                                ! OUT output field
     &    rmdi                                  ! IN  missing data indic
      LOGICAL                                                           &
     &    lcyclic,                                                      &
                                                ! IN .true. if cyclic EW
     &    lmasswt,                                                      &
                                                ! IN  TRUE if masswts OK
     &    land(row_length,rows),                                        &
                                                ! IN land mask
     &    sea(row_length,rows)                  ! IN sea mask
      CHARACTER*(*) cmessage                    ! OUT error message

!*----------------------------------------------------------------------
!
! STPARAM
!
!  Purpose: Meaningful PARAMETER names for STASH processing routines.
!           Both a long name and short name have been declared, to
!           reduce the existence of "magic" numbers in STASH.
!           Format is that first the address of the item is declare in
!           both long and short form. example is;
!             integer st_item_code,s_item  !Item number (declaration)
!             parameter(st_item_code=3,s_item=3)
!
!  Author:   S.Tett             Date:           22 January 1991
!
!  Model            Modification history from model version 3.0:
! version  Date
!   3.5    Mar. 95  Sub-models project.
!                   st_model_code=28 added to STLIST addresses
!                                   S.J.Swarbrick
!   4.2    27/11/96 mpp code: Added new stlist "magic numbers" :
!                   st_dump_output_length, st_dump_output_addr
!                                                       P.Burton
!   4.4    23/09/97 Add st_offset_code to the STASH list
!                   S.D. Mullerworth
!    4.4  02/12/96 Time mean timeseries added R A Stratton.
!    4.5  23/01/98 Added new stlist magic number
!                  st_dump_level_output_length
!    4.5  23/01/98 A
!    5.5  28/02/03 Original Modifications for WAM. M.Holt
!         06/08/00 Modification for parallelisation of WAM.
!                          Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!    6.0  08/09/03 Add st_riv_grid 23. C.Bunton
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!  Logical components covered: D70
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!--------------------------------------------------------------

      ! Internal model number address
      INTEGER,PARAMETER:: st_model_code = 28
      INTEGER,PARAMETER:: s_modl        = 28

      ! Section Number address
      INTEGER,PARAMETER:: st_sect_no_code = 2
      INTEGER,PARAMETER:: s_sect          = 2
      INTEGER,PARAMETER:: st_sect_code    = 2

      INTEGER,PARAMETER:: st_item_code=1,s_item=1 ! Item number address

      ! Processing Code address
      INTEGER,PARAMETER:: st_proc_no_code=3,s_proc=3

      ! subsidiary codes for st_proc_no_code now
      INTEGER,PARAMETER:: st_replace_code=1
      INTEGER,PARAMETER:: st_accum_code=2
      INTEGER,PARAMETER:: st_time_mean_code=3
      INTEGER,PARAMETER:: st_time_series_code=4
      INTEGER,PARAMETER:: st_max_code=5
      INTEGER,PARAMETER:: st_min_code=6
      INTEGER,PARAMETER:: st_append_traj_code=7
      INTEGER,PARAMETER:: st_time_series_mean=8
      INTEGER,PARAMETER:: st_variance_code=9

      ! Frequency (Input & output) addres
      INTEGER,PARAMETER:: st_freq_code=4,s_freq=4

      ! Offset for sampling
      INTEGER,PARAMETER:: st_offset_code=30,s_offs=30

      ! start timestep address
      INTEGER,PARAMETER:: st_start_time_code=5,s_times=5

      ! end timestep address
      INTEGER,PARAMETER:: st_end_time_code=6,s_timee=6

      ! period in timesteps address
      INTEGER,PARAMETER:: st_period_code=7,s_period=7

      ! infinite end/period value
      INTEGER,PARAMETER:: st_infinite_time=-1

      INTEGER,PARAMETER:: st_end_of_list=-1 !end-of-list marker in times

      ! grid point stuff
      ! gridpoint info address
      INTEGER,PARAMETER:: st_gridpoint_code=8,s_grid=8

      ! now subsid grid point stuff
      ! no masking done
      INTEGER,PARAMETER:: stash_null_mask_code=1,s_nomask=1

      ! land mask conds
      INTEGER,PARAMETER:: stash_land_mask_code=2,s_lndms=2

      ! sea mask code
      INTEGER,PARAMETER:: stash_sea_mask_code=3,s_seams =3

      ! processing options

      ! size of block for gridpoint code
      INTEGER,PARAMETER:: block_size=10

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: extract_base=block_size*0

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: extract_top=block_size*1

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_base=block_size*1

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_top=block_size*2

      ! max code for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_base=block_size*2

      ! base codes for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_top=block_size*3

      ! max code for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_base=block_size*3

      ! base codes for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_top=block_size*4

      ! max code for field mean subroutine
      INTEGER,PARAMETER:: field_mean_base=block_size*4

      ! base codes for field mean subroutine
      INTEGER,PARAMETER:: field_mean_top=block_size*5

      ! max code for global mean subroutine
      INTEGER,PARAMETER:: global_mean_base=block_size*5

      ! base codes for global mean subroutine
      INTEGER,PARAMETER:: global_mean_top=block_size*6

      ! Weighting

      ! weighting info address
      INTEGER,PARAMETER:: st_weight_code=9,s_weight=9

      INTEGER,PARAMETER:: stash_weight_null_code  =0,s_noweight  =0
      INTEGER,PARAMETER:: stash_weight_area_code  =1,s_areaweight=1
      INTEGER,PARAMETER:: stash_weight_volume_code=2,s_volweight =2
      INTEGER,PARAMETER:: stash_weight_mass_code  =3,s_massweight=3

      ! Domain definition

      ! row addresses
      INTEGER,PARAMETER:: st_north_code=12,s_north=12
      INTEGER,PARAMETER:: st_south_code=13,s_south=13
      INTEGER,PARAMETER:: st_west_code =14,s_west =14
      INTEGER,PARAMETER:: st_east_code =15,s_east =15

      ! Levels

      ! input bottom level address
      INTEGER,PARAMETER:: st_input_bottom=10,s_bottom =10

      ! special code
      INTEGER,PARAMETER:: st_special_code=100,s_special=100

      ! input top level address
      INTEGER,PARAMETER:: st_input_top=11,s_top=11

      ! output bottom level address
      INTEGER,PARAMETER:: st_output_bottom=21,s_outbot=21

      ! output top level address
      INTEGER,PARAMETER:: st_output_top=22,s_outtop=22

      INTEGER,PARAMETER:: st_model_level_code=1,s_model=1

      ! code for pressure leve
      INTEGER,PARAMETER:: st_pressure_level_code=2,s_press=2

      ! code for height levels
      INTEGER,PARAMETER:: st_height_level_code=3,s_height=3

      ! input code addres
      INTEGER,PARAMETER:: st_input_code=16,s_input=16

      ! input length of diagnostic address
      INTEGER,PARAMETER:: st_input_length=17,s_length=17

      ! output code address
      INTEGER,PARAMETER:: st_output_code=18,s_output=18

      ! Pointer to D1 addressing information
      ! Pos of item in D1 for relevant submodel
      INTEGER,PARAMETER:: st_position_in_d1=29,st_d1pos=29

      ! Output destination options

      INTEGER,PARAMETER:: st_dump=1
      INTEGER,PARAMETER:: st_secondary=2

      ! output length of diagnostic address
      INTEGER,PARAMETER:: st_output_length=19,s_outlen=19
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

      ! ptr to dump lookup header address
      INTEGER,PARAMETER:: st_lookup_ptr=23

      ! ptr into stash_series where control data address
      INTEGER,PARAMETER:: st_series_ptr=24

      ! subsid stuff for time series
      INTEGER,PARAMETER:: series_grid_type=1
      INTEGER,PARAMETER:: series_grid_code=0
      INTEGER,PARAMETER:: series_long_code=1
      INTEGER,PARAMETER:: series_size=2
      INTEGER,PARAMETER:: series_proc_code=3
      INTEGER,PARAMETER:: series_north=4
      INTEGER,PARAMETER:: series_south=5
      INTEGER,PARAMETER:: series_west=6
      INTEGER,PARAMETER:: series_east=7
      INTEGER,PARAMETER:: series_list_start=8
      INTEGER,PARAMETER:: series_list_end=9
      INTEGER,PARAMETER:: record_size=9

      ! Miscellaneous parameters

      ! system/user tag field in stlist address
      INTEGER,PARAMETER:: st_macrotag=25

      ! Pseudo-level list pointers

      ! pseudo-levels input list address
      INTEGER,PARAMETER:: st_pseudo_in=26

      ! pseudo-levels output list address
      INTEGER,PARAMETER:: st_pseudo_out=27

      ! Internal horizontal gridtype codes common to all diagnostics

      INTEGER,PARAMETER:: st_tp_grid =1 ! T-p grid
      INTEGER,PARAMETER:: st_uv_grid =2 ! u-v grid
      INTEGER,PARAMETER:: st_cu_grid =3 ! C-grid u point
      INTEGER,PARAMETER:: st_cv_grid =4 ! C-grid v point
      INTEGER,PARAMETER:: st_zt_grid =5 ! Zonal T-grid
      INTEGER,PARAMETER:: st_zu_grid =6 ! Zonal u-grid
      INTEGER,PARAMETER:: st_mt_grid =7 ! Meridional T-grid
      INTEGER,PARAMETER:: st_mu_grid =8 ! Meridional u-grid
      INTEGER,PARAMETER:: st_riv_grid= 23    ! river_routing grid
      INTEGER,PARAMETER:: st_scalar  =9 ! Scalar (ie. single value)
      INTEGER,PARAMETER:: st_wam_all= 60    ! Wam Field on Full Grid
      INTEGER,PARAMETER:: st_wam_sea= 62    ! Wam Field on Sea Points

! STPARAM end
!LL  Comdeck: STERR ----------------------------------------------------
!LL
!LL  Purpose: PARAMETER names for STASH processing error codes;
!LL           fatal errors have positive codes, warnings negative.
!LL
!LL  Author:   S.Tett
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.3  16/09/93  Add st_illegal_weight error code.
!LL                   Added st_no_data for mpp code
!LL                   (means a processor does not contain any data
!LL                    for a given subdomain)                 P.Burton
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D70
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!
! Warning codes
!
         integer st_upper_less_lower ! warning code for bad domain
         parameter(st_upper_less_lower=-1)

         integer st_not_supported ! warning code for unsupported routine
         parameter(st_not_supported=-2)
         integer st_no_data,st_nd ! indicates no data on a processor
         parameter(st_no_data=-3,st_nd=-3)
!
! Error codes
!
         integer st_bad_array_param ! error code for dodgy array params
         parameter(st_bad_array_param=1)

         integer st_bad_address     ! error code for address violation
         parameter(st_bad_address=2)

         integer st_unknown ! error code for unknown option
         parameter(st_unknown=3)

         integer st_bad_wraparound ! error code for illegal wraparound
         parameter(st_bad_wraparound=4)

         integer st_illegal_weight ! error code for illegal weighting
         parameter(st_illegal_weight=9)

         integer unknown_weight ! error code for an unknown weight
         parameter(unknown_weight=10)

         integer unknown_mask ! error code for an unknown mask
         parameter(unknown_mask=11)

         integer unknown_processing ! error code for unknown processing
         parameter(unknown_processing=12)

         integer nonsense ! error code for general nonsense request
         parameter(nonsense=13)
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
!L
!L external routines
!L
      EXTERNAL stextc ! extracts the field
      EXTERNAL stcolm ! computes the column mean
      EXTERNAL stzonm ! computes the zonal mean
      EXTERNAL stmerm ! computes the meridional mean
      EXTERNAL stglom ! computes the global mean
      EXTERNAL stfieldm ! computes the field mean
!L
!L local variables
!L
      LOGICAL lwrap                ! TRUE if output field wraparound EW
      LOGICAL lmasswt_strict       ! copy of lmasswt - but set to false
!                                  ! if mass weighting is not requested
      INTEGER xstart,ystart        ! lower left hand corner coords
      INTEGER xend,yend            ! upper right hand corner coords
      INTEGER processing_code      ! what kind of mean  will be done
      INTEGER what_level           ! what type of input level
      INTEGER what_mask            ! what mask is used
      INTEGER what_weight          ! what weighting is used

      INTEGER i,j,k                                                     &
                                   ! loop counters
     &,model_level                                                      &
                                   ! model level
     &,this_index_lev_masswt       ! level index for mass weighting
                                   ! (=1 if no mass weighting or no
                                   !  model level weights available)

      INTEGER                                                           &
! global versions of the extracted area domain limits
     &  global_xstart,global_xend,global_ystart,global_yend

! workspace arrays containining weighting factors and masks.
      REAL                                                              &
     &  area_weight(1-halo_x:row_length+halo_x+1,                       &
     &              1-halo_y:no_rows+halo_y)                            &
     &, pstar_weight(1-halo_x:row_length+halo_x+1,                      &
     &               1-halo_y:no_rows+halo_y,                           &
     &               no_of_levels_masswt)                               &
     &, pstar_interp(row_length,no_rows)
      LOGICAL                                                           &
     &  mask(1-halo_x:row_length+halo_x+1,                              &
     &       1-halo_y:no_rows+halo_y)


!L----------------------------------------------------------------------
!L 1. Set up local variables
!L
      xstart=control(st_west_code)
      xend=control(st_east_code)
      ystart=control(st_south_code)  ! NOTE: Grid is assumed to be
      yend=control(st_north_code)    !       oriented south-to-north

      global_xstart=xstart
      global_ystart=ystart
      global_xend=xend
      global_yend=yend

! and calculate what the local subdomain limits are:
! DEPENDS ON: global_to_local_subdomain
        CALL GLOBAL_TO_LOCAL_SUBDOMAIN( .TRUE.,.TRUE.,                  &
     &                                  GR,halo_type,mype,              &
     &                                  global_ystart,global_xend,      &
     &                                  global_yend,global_xstart,      &
     &                                  ystart,xend,yend,xstart)

! Check if wraparound field
      IF (xstart >  xend) THEN

        IF (lcyclic) THEN
          xend=xend+blsize(1,fld_type)
! subtract two halos as we don't wish to include halos at the end
! and start of the row within the wrap around domain
          lwrap=.TRUE.
        ELSE
          icode=st_bad_wraparound    ! wraparound illegal unless cyclic
          GOTO 999
        ENDIF

      ELSE
        lwrap=.FALSE.
      ENDIF

      IF (global_xstart  >   global_xend) THEN
        IF (lcyclic) THEN
          global_xend=global_xend+glsize(1,fld_type)
        ELSE
          icode=st_bad_wraparound  ! wraparound illegal unless cyclic
          GOTO 999
        ENDIF
      ENDIF

      processing_code=control(st_gridpoint_code)
      what_level=control(st_input_bottom)
      what_mask=mod(processing_code,block_size)
      what_weight=control(st_weight_code)
!L
!L 1.1 Prevent masking or weighting if input field is not 2D in extent
!L     - weighting and masking is assumed to have been done outside.
!L
      IF ( (.NOT.(st_grid == st_tp_grid.OR.st_grid == st_uv_grid.OR.    &
     &            st_grid == st_cu_grid.OR.st_grid == st_cv_grid.OR.    &
     &             st_grid == st_riv_grid))                             &
     &      .AND.(what_mask   /= stash_null_mask_code   .OR.            &
     &            what_weight /= stash_weight_null_code) ) THEN
       icode=st_not_supported
       cmessage='SPATIAL : Masking/weighting unsupported - non 2D field'
       GOTO 999
      ENDIF

! Check for supported weighting and masking options

      IF (.NOT. ((what_weight  ==  stash_weight_null_code) .OR.         &
     &           (what_weight  ==  stash_weight_area_code) .OR.         &
     &           (what_weight  ==  stash_weight_volume_code) .OR.       &
     &           (what_weight  ==  stash_weight_mass_code) ) ) THEN
        cmessage='SPATIAL : Unrecognized weighting option'
        icode=unknown_weight
        GOTO 999
      ENDIF

      IF (.NOT. ((what_mask  ==  stash_null_mask_code) .OR.             &
     &           (what_mask  ==  stash_land_mask_code) .OR.             &
     &           (what_mask  ==  stash_sea_mask_code ) ) ) THEN
        cmessage='SPATIAL : Unrecognized masking option'
        icode=unknown_mask
        GOTO 999
      ENDIF

      IF (what_weight  ==  stash_weight_volume_code) THEN
        cmessage='SPATIAL : Volume-weighting not supported'
        icode=st_illegal_weight
        GOTO 999
      ENDIF

! Set lmasswt_strict - copy of lmasswt, but set to false is mass
! weighting not requested

      lmasswt_strict=                                                   &
     &  (lmasswt .AND. (what_weight  ==  stash_weight_mass_code))

! Precalculate weighting and masking arrays
! I've used IF tests inside the loops, but since the logical
! expressions are invariant wrt i and j, the compiler will
! move them outside the DO loops. It makes the code a lot shorter!

! NOTE that neither area-weights or mass-weights are normalised so
! that the interpretation of weighted diagnostics is non-obvious. Also
! the PP lookup header has no switch for indicating whether or not such
! processing has taken place. More careful treatment of horizontal
! interpolation is required for stricter accuracy.


! area weighting
      IF (what_weight  ==  stash_weight_null_code) THEN
          ! Ensure initialisation of weight array including halos
          area_weight (:,:)   = 1.0
      ELSE ! some form of area weighting will be required
          DO j=1,no_rows
          DO i=1,row_length
             IF (st_grid  ==  st_cv_grid) THEN
                 area_weight(i,j)=cos_v_latitude(i,j)
             ELSE
! NOTE that this is will not be accurate for C-u grid variables for
! LAMs, since cos_theta_latitude will differ between theta,u positions
                 area_weight(i,j)=cos_theta_latitude(i,j)
             ENDIF
         ENDDO ! i
         ENDDO ! j

      ENDIF    ! what_weight


! mass weighting
      IF ((what_weight  ==  stash_weight_null_code) .OR.                &
     &    (what_weight  ==  stash_weight_area_code)) THEN
! No mass weighting is required
          this_index_lev_masswt = 1

          ! Ensure initialisation of weight array including halos
          pstar_weight(:,:,this_index_lev_masswt) = 1.0
      ELSE

! Mass weighting requested

! Ensure that halos are initialised
        DO j=no_rows-1,no_rows
        DO i=1,row_length
           pstar_interp(i,j) =1.0
        ENDDO !i
        ENDDO !j

! Interpolate pstar to correct horizontal staggering
        IF     (st_grid  ==  st_cu_grid) THEN
! NOT YET CORRECT: pstar has no halos! So set to nearby value
!          CALL P_TO_U(pstar,row_length,rows,1,0,0,pstar_interp)
           DO j=1,no_rows
           DO i=1,row_length
              pstar_interp(i,j) =pstar(i,j)
           ENDDO !i
           ENDDO !j

        ELSEIF (st_grid  ==  st_cv_grid) THEN
! NOT YET CORRECT: pstar has no halos! So set to nearby value
!          CALL P_TO_V(pstar,row_length,rows,1,0,0,pstar_interp)
           DO j=1,no_rows
           DO i=1,row_length
              pstar_interp(i,j) =pstar(i,j)
           ENDDO !i
           ENDDO !j

        ELSE
          DO j=1,no_rows
            DO i=1,row_length
              pstar_interp(i,j)=pstar(i,j)
            ENDDO
          ENDDO
        ENDIF

        IF(lmasswt) THEN  ! model level mass weights available

          this_index_lev_masswt = this_index_lev

          DO k=1,no_of_levels_masswt
            model_level = level_list(k)
            IF(model_level == model_levels) THEN  ! top level

              DO j=1,no_rows
              DO i=1,row_length
                pstar_weight(i,j,k) = p(i,j,model_level)
              ENDDO !i
              ENDDO !j

            ELSE      ! not top level

              DO j=1,no_rows
              DO i=1,row_length
! Only accurate for variables on theta levels
                pstar_weight(i,j,k) = p(i,j,model_level) -              &
     &                                  p(i,j,model_level+1)
              ENDDO !i
              ENDDO !j
            ENDIF

          ENDDO ! k

        ELSE              ! no model level mass weights available:
                          ! weight by pstar

            this_index_lev_masswt = 1

            DO j=1,no_rows
            DO i=1,row_length
               pstar_weight(i,j,this_index_lev_masswt) =                &
     &                                               pstar_interp(i,j)
            ENDDO !i
            ENDDO !j
! Horizontal interpolation may be required - this should be done here:

        ENDIF   ! lmasswt
      ENDIF

! masking

      DO j=1,no_rows
        DO i=1,row_length
          IF (what_mask  ==  stash_land_mask_code) THEN
            mask(i,j)=land(i,j)
          ELSEIF (what_mask  ==  stash_sea_mask_code) THEN
            mask(i,j)=sea(i,j)
          ELSE
            mask(i,j)=.TRUE.
          ENDIF
        ENDDO
      ENDDO

! Update the halo areas of the weighting/masking arrays

!  [Note that for lams at UM5.2 and before, valid values for rim
!   boundaries would be needed using FILL_EXTERNAL_HALOS calls for
!   area_weight and pstar_weight arrays, but this should now be
!   superseded by initialising full arrays]

      ! Update halos only if halos present (standard diagnostics have no
      ! halos, and if some weighting required (probably defunct
      ! functionality)
      IF( (halo_x  /=  0 .OR. halo_y  /=  0) .AND.                      &
     &    (what_weight  /=  stash_weight_null_code .OR.                 &
     &     what_mask    /=  stash_null_mask_code )     ) THEN

! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS( AREA_WEIGHT(1-halo_x:row_length+halo_x,:),      &
     &     ROW_LENGTH,NO_ROWS,1,halo_x,halo_y,                          &
     &                fld_type_p,.FALSE.)
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS( PSTAR_WEIGHT(1-halo_x:row_length+halo_x,:,:),   &
     &     ROW_LENGTH,NO_ROWS,                                          &
     &                no_of_levels_masswt,halo_x,halo_y,                &
     &                fld_type_p,.FALSE.)
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS( MASK(1-halo_x:row_length+halo_x,:),             &
     &     ROW_LENGTH,NO_ROWS,1,halo_x,halo_y,                          &
     &                fld_type_p,.FALSE.)
      ENDIF
!L----------------------------------------------------------------------
!L 2. Call service routine to perform required processing
!L
!L 2.1 Extract sub field (single level at a time)
!L
      IF (processing_code <  extract_top.and.                           &
     &    processing_code >= extract_base) THEN
        n_rows_out=(yend+1)-ystart
        n_cols_out=(xend+1)-xstart

        IF (                                                            &
     &   (xstart  /=  st_no_data) .AND. (xend  /=  st_no_data) .AND.    &
     &   (ystart  /=  st_no_data) .AND. (yend  /=  st_no_data)) THEN

! DEPENDS ON: stextc
        CALL STEXTC(fieldin,vx,vy,fld_type,halo_type,                   &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)

        ELSE  ! just set values to non NaN
          DO i=1,lenout
            fieldout(i)=0.0
          ENDDO
        ENDIF

!L
!L 2.2 Calculate column mean (over multiple levels indexed by index_lev)
!L
      ELSEIF (processing_code <  vert_mean_top.and.                     &
     &        processing_code >  vert_mean_base) THEN
        n_rows_out=yend+1-ystart
        n_cols_out=xend+1-xstart

        IF (                                                            &
     &   (xstart  /=  st_no_data) .AND. (xend  /=  st_no_data) .AND.    &
     &   (ystart  /=  st_no_data) .AND. (yend  /=  st_no_data)) THEN

! DEPENDS ON: stcolm
        CALL STCOLM(fieldin,vx,vy,vz,fld_type,halo_type,                &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              fieldout,index_lev,no_of_levels,                    &
     &              pstar_weight,                                       &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)

        ELSE  ! just set values to non NaN
          DO i=1,lenout
            fieldout(i)=0.0
          ENDDO
        ENDIF

!L
!L 2.3 Calculate zonal mean (single level at a time)
!L
      ELSEIF (processing_code <  zonal_mean_top.and.                    &
     &        processing_code >  zonal_mean_base) THEN
        n_rows_out=yend+1-ystart
        n_cols_out=1
! DEPENDS ON: stzonm
        CALL STZONM(fieldin,vx,vy,fld_type,gr,halo_type,                &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.4 Calculate meridional mean (single level at a time)
!L
      ELSEIF (processing_code <  merid_mean_top.and.                    &
     &        processing_code >  merid_mean_base) THEN
        n_rows_out=1
        n_cols_out=xend+1-xstart
! DEPENDS ON: stmerm
        CALL STMERM(fieldin,vx,vy,fld_type,gr,halo_type,                &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.5 Calculate field mean (single level at a time)
!L
      ELSEIF (processing_code <  field_mean_top.and.                    &
     &        processing_code >  field_mean_base) THEN
        n_rows_out=1
        n_cols_out=1
! DEPENDS ON: stfieldm
        CALL STFIELDM(fieldin,vx,vy,fld_type,gr,halo_type,              &
     &                lwrap,lmasswt_strict,                             &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,                                           &
     &              pstar_weight(1-halo_x,1-halo_y,                     &
     &                           this_index_lev_masswt),                &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.6 Calculate global mean (over multiple levels)
!L
      ELSEIF (processing_code <  global_mean_top.and.                   &
     &        processing_code >  global_mean_base) THEN
        n_rows_out=1
        n_cols_out=1
! DEPENDS ON: stglom
        CALL STGLOM(fieldin,vx,vy,vz,fld_type,gr,halo_type,             &
     &              lwrap,lmasswt_strict,                               &
     &              xstart,ystart,xend,yend,                            &
     &              global_xstart,global_ystart,global_xend,global_yend,&
     &              fieldout,index_lev,no_of_levels,                    &
     &              pstar_weight,                                       &
     &              area_weight,mask,                                   &
     &              what_level,what_mask,what_weight,rmdi,              &
     &              icode,cmessage)
!L
!L 2.7 Invalid processing option
!L
      ELSE
        icode=unknown_processing
        write(cmessage,111)'unknown processing option',                 &
     &    processing_code
      ENDIF
!L
  999 CONTINUE
111   format('SPATIAL : >>FATAL ERROR <<',a40,i5)
!
      RETURN
      END SUBROUTINE SPATIAL
