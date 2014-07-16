
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: STZONM ---------------------------------------------------
!LL
!LL  Purpose: Calculate weighted zonal mean within a region specified
!LL           by a lower left hand and upper right hand corner.
!LL           (STASH service routine).
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.Johns/S.Tett
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D713
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE STZONM(fieldin,vx,vy,fld_type,gr,halo_type,            &
     &                  lwrap,lmasswt,                                  &
     &                  xstart,ystart,xend,yend,                        &
     &                  global_xstart,global_ystart,                    &
     &                  global_xend,global_yend,                        &
     &                  fieldout,                                       &
     &                  pstar_weight,                                   &
     &                  area_weight,mask,                               &
     &                  level_code,mask_code,weight_code,rmdi,          &
     &                  icode,cmessage)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &    vx,vy,                                                        &
                                                ! IN  input field size
     &    fld_type,                                                     &
                                                ! IN field type (u/v/p)
     &    gr,                                                           &
                                                ! IN input fld grid
     &    halo_type,                                                    &
                                                ! IN halo type
     &    xstart,ystart,                                                &
                                                ! IN  lower LH corner
     &    xend,yend,                                                    &
                                                ! IN  upper RH corner
     &    global_xstart,global_ystart,                                  &
                                                ! IN global versions of
     &    global_xend,  global_yend,                                    &
                                                ! IN xstart etc.
     &    level_code,                                                   &
                                                ! IN  input level code
     &    mask_code,                                                    &
                                                ! IN  masking code
     &    weight_code,                                                  &
                                                ! IN  weighting code
     &    icode                                 ! OUT error return code
      CHARACTER*(*)                                                     &
     &    cmessage                              ! OUT error return msg
      LOGICAL                                                           &
     &    lwrap,                                                        &
                                                ! IN  TRUE if wraparound
     &    lmasswt,                                                      &
                                                ! IN  TRUE if masswts OK
     &    mask(vx+1,vy)                         ! IN  mask array
      REAL                                                              &
     &    fieldin(vx,vy),                                               &
                                                ! IN  input field
     &    fieldout(ystart:yend),                                        &
                                                ! OUT output field
     &    pstar_weight(vx+1,vy),                                        &
                                                ! IN  mass weight factor
     &    area_weight(vx+1,vy),                                         &
                                                ! IN  area weight factor
! (already interpolated to the correct grid and
!  set to 1.0 where no area weighting is required)
     &    rmdi                                  ! IN  missing data indic
!*----------------------------------------------------------------------
!
! External subroutines called
!
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
!
! Local variables
!
        INTEGER i,ii,j   ! ARRAY INDICES FOR VARIABLE


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
! Processor co-ordinates of processors at the corners of the
! processed subdomain
     &  proc_top_left_x,proc_top_left_y                                 &
     &, proc_bot_right_x,proc_bot_right_y                               &

! size of the full subarea in
     &, zonal_sum_global_len_x                                          &

! loop variables for loops over processors in subdomain
     &, proc_x,proc_y                                                   &

! real processor x co-ordinate - when proc_x > nproc_x is just
! proc_x-nproc_x
     &, eff_proc_x                                                      &

! processor id of processor (proc_x,proc_y)
     &, proc_id                                                         &

! definition of the extracted subarea array on processor proc_id
     &, local_array_top_left_x,local_array_top_left_y                   &
     &, local_array_bot_right_x,local_array_bot_right_y                 &

! definition of the real data contained within the extracted
! subarea on processor proc_id (ie. not including halos)
     &, local_data_top_left_x,local_data_top_left_y                     &
     &, local_data_bot_right_x,local_data_bot_right_y                   &

! size in the x dimension of the subarea array on proc_id
     &, local_array_size_x                                              &

! length of (partial) local zonal data to be sent
     &, local_send_size_x                                               &

! offset of data to be sent from start of local row
     &, local_send_off_x                                                &

! position in the full zonal row of (partial) data to be sent
     &, pos_in_zonal_array                                              &

! number of (partial) zonal mean rows on processor proc_id
     &, local_n_rows_to_send                                            &

! the first (partial) zonal mean row to be sent from proc_id
     &, local_send_off_y                                                &

! the global zonal mean number of the first local (partial) zonal
! mean row to be sent from proc_id
     &, global_zonal_row_number_start                                   &

! loop counter for loop over rows to send
     &, row                                                             &

! index of row in proc_id's array of local data
     &, local_row                                                       &

! global zonal row number of a row
     &, global_zonal_row_id                                             &

! processor which this zonal mean row will be sent to
     &, dest_proc                                                       &

! row number on the destination processor
     &, work_dest_row_id                                                &

! number of items of zonal data to send and receive
     &, n_send_data , n_rec_data                                        &

! number of final zonal means to send and receive
     &, n_send_means , n_rec_means                                      &

! number of rows of (full) zonal data on this processor
     &, n_rows_full_zonal_data                                          &

! size of local_sum_arrays and global_sum_arrays
     &, local_sum_array_len                                             &
     &, global_sum_array_len                                            &


! arguments for GCOM routines
     &, flag , info                                                     &

! dummy variables (unused return values from subroutine calls)
     &, dummy1,dummy2


      LOGICAL                                                           &

! indicates if the subarea requested for zonal meaning wraps over
! zero longitude
     &  lwrap_zonal_mean                                                &

! indicates if the subdomain contains processors which hold both
! the start and end of the subdomain, which wraps over zero
! longitude
     &, lwrap_proc                                                      &

! indicates that a full field is being zonal meaned
     &, fullfield

      REAL                                                              &
! temporary variables used in calculation  of zonal means
     &  zonal_sum_top,zonal_sum_bot

      INTEGER                                                           &

! Send/receive maps for zonal data arrays to be summed
     &  send_data_map(7,2*(yend-ystart+1))                              &
     &, rec_data_map(7,2*(global_yend-global_ystart+1)*nproc)           &

! send/receive maps for zonal means
     &, send_means_map(7,2*(global_yend-global_ystart+1))               &
     &, rec_means_map(7,2*(yend-ystart+1))

! Weighted version of fieldin
      REAL local_sum_array_top(xstart:xend,ystart:yend)
! Weights applied to fieldin
      REAL local_sum_array_bot(xstart:xend,ystart:yend)


      INTEGER                                                           &
! Sizes of the global_sum_arrays defined below
     &  global_sum_array_sizex,global_sum_array_sizey

      REAL                                                              &
! Collected versions of fieldin and the weights containing
! whole (subarea) rows of zonal data
     &  global_sum_array_top(global_xstart:global_xend,                 &
     &                       global_yend-global_ystart+1)               &
     &, global_sum_array_bot(global_xstart:global_xend,                 &
     &                       global_yend-global_ystart+1)               &

! Calculated zonal means on the calculating processor
     &, zonal_mean_array(global_yend-global_ystart+1)



! Integer function used for obtaining field type
      INTEGER GET_FLD_TYPE


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

        REAL SUMZTOP(ystart:yend)
        REAL SUMZBOT(ystart:yend)

!L----------------------------------------------------------------------
!L 0. Initialise sums
!L
!FPP$ NOINNER R
      DO j=ystart,yend
        SUMZTOP(j)=0.0
        SUMZBOT(j)=0.0
      ENDDO
!L----------------------------------------------------------------------

! pstar_weight and area_weight arrays contain appropriate
! weighting factors, interpolated to the correct grid, for
! mass weighting and area weighting respectively. If either type
! of weighting is not required, the relevant array is set to 1.0
! The mask array contains appropriate masking

! Create arrays of weighted data suitable to be summed

! Only do the calculations if some of the subarea is contained
! within this processor
      IF ((xstart  /=  st_no_data) .AND. (xend  /=  st_no_data) .AND.   &
     &    (ystart  /=  st_no_data) .AND. (yend  /=  st_no_data)) THEN

        DO i=xstart,xend
            IF ( lwrap .AND.                                            &
     &          (i  >   (lasize(1,fld_type,halo_type)-                  &
     &                   halosize(1,halo_type)))) THEN
! miss halos on wrap around
              ii=i-blsize(1,fld_type)
          ELSE
            ii=i
          ENDIF
          DO j=ystart,yend
            IF (mask(ii,j)) THEN
                local_sum_array_bot(i,j)=                               &
     &            pstar_weight(ii,j)
                local_sum_array_top(i,j)=                               &
     &            fieldin(ii,j)*pstar_weight(ii,j)
            ELSE
              local_sum_array_bot(i,j)=0.0
              local_sum_array_top(i,j)=0.0
            ENDIF
          ENDDO
        ENDDO

      ENDIF  ! if this processor contains any of the subarea

! Initialise fieldout array - so all PE's have valid data
! (Only PEs on left of subdomain get the zonal means)
      DO i=ystart,yend
        fieldout(i)=0.0
      ENDDO

! The local_sum_arrays must be distributed so that complete
! sub-area rows exist on processors, so that a reproducible sum
! can be carried out.
! The following code calculates where the local_sum_array data
! must be sent to, and where the final answers must be sent back to

! 0.0 : Initialise variables defining the size of the arrays
!       global_sum_arrays

      global_sum_array_sizex=global_xend-global_xstart+1
      global_sum_array_sizey=global_yend-global_ystart+1

      local_sum_array_len=((xend-xstart)+1)*((yend-ystart)+1)
      global_sum_array_len=global_sum_array_sizex*                      &
     &                     global_sum_array_sizey

! Set a logicial indicating if the area being meaned is the
! full field


      fullfield= ( (global_xstart  ==  1                 ) .AND.        &
     &             (global_xend    ==  glsize(1,fld_type)) .AND.        &
     &             (global_yend    ==  glsize(2,fld_type))      )

! Calculate the length of the full zonal subarea

      zonal_sum_global_len_x=global_xend-global_xstart+1

! 1.0 Find the set of processors covering the requested sub-area

! DEPENDS ON: global_to_local_rc
      CALL GLOBAL_TO_LOCAL_RC(gr,halo_type,                             &
     &   global_xstart , global_ystart,                                 &
     &   proc_top_left_x, proc_top_left_y,                              &
     &   dummy1,dummy2)

! DEPENDS ON: global_to_local_rc
      CALL GLOBAL_TO_LOCAL_RC(gr,halo_type,                             &
     &   global_xend,global_yend,                                       &
     &   proc_bot_right_x, proc_bot_right_y,                            &
     &   dummy1,dummy2)

! Set a logical to indicate if the zonal mean area required
! wraps over zero longitude

      lwrap_zonal_mean=                                                 &
     &  ((global_xend  >   glsize(1,fld_type)) .OR.                     &
     &   (global_xend  <   global_xstart))

! If there is a wrap around over 0 longitude, ensure that
! proc_bot_right_x > proc_top_left_x

      IF (lwrap_zonal_mean)                                             &
     &  proc_bot_right_x=proc_bot_right_x+nproc_x

! Set up a logical to indicate if a processor in the subdomain
! contains both the start and end of a zonal mean which wraps over
! zero longitude. If TRUE, some extra work is required at this
! processor as it contains data for two non-contiguous parts
! of the zonal mean

      lwrap_proc=(proc_bot_right_x  ==  proc_top_left_x+nproc_x)

! 2.0 Loop over all the processors in the subdomain, and set
!     up the send/receive maps defining the redistribution
!     of data

      n_send_data=0            ! number of items of data to send
      n_rec_data=0             ! number of items of data to receive
      n_send_means=0           ! number of zonal means I will send
      n_rec_means=0            ! number of zonal means I will receive
      n_rows_full_zonal_data=0 ! number of rows of data I will mean

      DO proc_y=proc_top_left_y , proc_bot_right_y

        DO proc_x=proc_top_left_x , proc_bot_right_x

          eff_proc_x=MOD(proc_x,nproc_x)
          proc_id=eff_proc_x+proc_y*nproc_x

! 2.1  Find the size of the array containing the zonal arrays on
!      processor proc_id

! DEPENDS ON: global_to_local_subdomain
          CALL GLOBAL_TO_LOCAL_SUBDOMAIN(.TRUE.,.TRUE.,                 &
     &      gr,halo_type,proc_id,                                       &
     &      global_ystart,global_xend,                                  &
     &      global_yend,global_xstart,                                  &
     &      local_array_top_left_y,local_array_bot_right_x,             &
     &      local_array_bot_right_y,local_array_top_left_x)

! 2.2 Using this information, calculate the size of this array in
!     the x dimension. If the data is wrapped round, the calculation
!     is done differently:

          IF (local_array_top_left_x  <=  local_array_bot_right_x)      &
     &    THEN
            local_array_size_x=                                         &
     &          local_array_bot_right_x-local_array_top_left_x+1
          ELSE ! wrap around
            local_array_size_x=                                         &
     &          local_array_bot_right_x-local_array_top_left_x+1+       &
     &          g_blsize(1,fld_type,proc_id)
          ENDIF

! 2.3 Find out the size of the actual zonal mean data within the
!     subarea array on processor proc_id

! DEPENDS ON: global_to_local_subdomain
          CALL GLOBAL_TO_LOCAL_SUBDOMAIN(.FALSE.,.FALSE.,               &
     &      gr,halo_type,proc_id,                                       &
     &      global_ystart,global_xend,                                  &
     &      global_yend,global_xstart,                                  &
     &      local_data_top_left_y,local_data_bot_right_x,               &
     &      local_data_bot_right_y,local_data_top_left_x)

! 2.4 Calculate various quantities, which require different
!     calculations depending on if  LWRAP_PROC is .TRUE.,
!     and if so, if this processor contains both start and
!     end of the zonal data
!     local_send_size_x  : the length of data to be sent
!     local_send_off_x   : the offset of this data from the
!                          start of row
!     pos_in_zonal_array : position of this data in the full
!                          zonal array

          IF ((LWRAP_PROC) .AND. (proc_x  ==  proc_top_left_x)) THEN
! Processor containing start and end of zonal mean - but here
! we're interested only in the start segment

            local_send_size_x=                                          &
     &        g_lasize(1,fld_type,halo_type,proc_id) -                  &
     &                             local_data_top_left_x -              &
     &                             halosize(1,halo_type)+1
            local_send_off_x=                                           &
     &        local_data_top_left_x-local_array_top_left_x
            pos_in_zonal_array=                                         &
     &        g_datastart(1,proc_id)+local_data_top_left_x              &
     &              - halosize(1,halo_type) - global_xstart

          ELSEIF ((LWRAP_PROC) .AND.                                    &
     &            (proc_x  ==  proc_bot_right_x)) THEN
! Processor containing start and end of zonal mean - but here
! we're interested only in the end segment

            local_send_size_x=local_data_bot_right_x                    &
     &                             - halosize(1,halo_type)
            local_send_off_x=local_array_size_x-local_send_size_x
            pos_in_zonal_array=                                         &
     &        zonal_sum_global_len_x-local_send_size_x+1

          ELSE
! all other processors

            local_send_size_x=                                          &
     &        local_data_bot_right_x-local_data_top_left_x+1
            local_send_off_x=                                           &
     &        local_data_top_left_x-local_array_top_left_x
            pos_in_zonal_array=                                         &
     &        g_datastart(1,proc_id)+local_data_top_left_x              &
     &                  - halosize(1,halo_type) - global_xstart

          ENDIF

          IF (pos_in_zonal_array  <   1) THEN
! This means the sub-area wraps over zero longitude - so to get
! the correct position in the array we add the global row length
            pos_in_zonal_array=pos_in_zonal_array+glsize(1,fld_type)
          ENDIF

! 2.5 Find the number of zonal mean rows to be sent from this
!     processor

          local_n_rows_to_send=                                         &
     &      local_data_bot_right_y-local_data_top_left_y+1

! 2.6 and the first row to be sent from this processor

          local_send_off_y=                                             &
     &      local_data_top_left_y-local_array_top_left_y

! 2.7 Calculate which global zonal mean is the first one to
!     send from this processor

          global_zonal_row_number_start=                                &
     &      g_datastart(2,proc_id)+local_data_top_left_y                &
     &            - halosize(2,halo_type) - global_ystart

! 2.8 Loop over rows and construct send/receive maps

          DO row=1,local_n_rows_to_send

! 2.8.1 Find the local_row index on proc_id, and the global zonal
!       row index of this row

            local_row=row+local_send_off_y
            global_zonal_row_id=global_zonal_row_number_start+row-1

! 2.8.2 and find the destination processor of this row, and
!       where on this processor it will be sent to

            dest_proc=MOD(global_zonal_row_id-1,nproc)
            work_dest_row_id=((global_zonal_row_id-1)/nproc)+1

! 2.8.3 If this processor is proc_id construct a send_data_map
!       entry for this row of data

            IF (mype  ==  proc_id) THEN

              n_send_data = n_send_data+1

              send_data_map(S_DESTINATION_PE,n_send_data)=              &
     &          dest_proc
              send_data_map(S_BASE_ADDRESS_IN_SEND_ARRAY,               &
     &                      n_send_data)=                               &
     &          (local_row-1)*local_array_size_x +                      &
     &          local_send_off_x+1
              send_data_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,               &
     &                      n_send_data)=1
              send_data_map(S_STRIDE_IN_SEND_ARRAY,n_send_data)=1
              send_data_map(S_ELEMENT_LENGTH,n_send_data)=              &
     &          local_send_size_x
              send_data_map(S_BASE_ADDRESS_IN_RECV_ARRAY,               &
     &                      n_send_data)=                               &
     &          (work_dest_row_id-1)*global_sum_array_sizex +           &
     &          pos_in_zonal_array
              send_data_map(S_STRIDE_IN_RECV_ARRAY,n_send_data)=1

! 2.8.3.1 If this processor is on LHS of the subarea, then it is
!         responsible for holding the final zonal mean values.
!         So we must set up a rec_means_map entry to allow the
!         zonal mean value for this row to be returned.

              IF (proc_x  ==  proc_top_left_x) THEN

                n_rec_means = n_rec_means+1

                rec_means_map(R_SOURCE_PE,n_rec_means)=dest_proc
                IF (fullfield) THEN ! We don't want halos
                  rec_means_map(R_BASE_ADDRESS_IN_RECV_ARRAY,           &
     &                          n_rec_means)=                           &
     &              local_row-halosize(2,halo_type)
                ELSE ! halos are automatically removed
                  rec_means_map(R_BASE_ADDRESS_IN_RECV_ARRAY,           &
     &                          n_rec_means)=                           &
     &              local_row
                ENDIF
                rec_means_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,             &
     &                        n_rec_means)= 1
                rec_means_map(R_STRIDE_IN_RECV_ARRAY,                   &
     &                        n_rec_means)= 1
                rec_means_map(R_ELEMENT_LENGTH,n_rec_means)=1
                rec_means_map(R_BASE_ADDRESS_IN_SEND_ARRAY,             &
     &                        n_rec_means)=                             &
     &            work_dest_row_id
                rec_means_map(R_STRIDE_IN_SEND_ARRAY,                   &
     &                        n_rec_means)=1
              ENDIF

            ENDIF

! 2.8.4 If this processor is dest_proc construct a rec_data_map
!       entry for this row of data

            IF (mype  ==  dest_proc) THEN

              IF (proc_x  ==  proc_top_left_x)                          &
! increment counter of full zonal rows on this processor
     &          n_rows_full_zonal_data=n_rows_full_zonal_data+1

              n_rec_data = n_rec_data+1

              rec_data_map(R_SOURCE_PE,n_rec_data)=                     &
     &          proc_id
              rec_data_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_rec_data)=    &
     &          (work_dest_row_id-1)*global_sum_array_sizex +           &
     &          pos_in_zonal_array
              rec_data_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_rec_data)=    &
     &          1
              rec_data_map(R_STRIDE_IN_RECV_ARRAY,n_rec_data)=1
              rec_data_map(R_ELEMENT_LENGTH,n_rec_data)=                &
     &          local_send_size_x
              rec_data_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_rec_data)=    &
     &          (local_row-1)*local_array_size_x +                      &
     &          local_send_off_x+1
              rec_data_map(R_STRIDE_IN_SEND_ARRAY,n_rec_data)=1

! 2.8.4.1 Set up the send_means_map entry for sending the
!         resulting zonal mean of this row back to
!         the processor at the LHS of the subarea.
!         We only need to do this once per row (not for
!         each value of proc_x).

              IF (proc_x  ==  proc_top_left_x) THEN

                n_send_means = n_send_means+1

                send_means_map(S_DESTINATION_PE,n_send_means)=          &
     &            proc_id
                send_means_map(S_BASE_ADDRESS_IN_SEND_ARRAY,            &
     &                         n_send_means)=work_dest_row_id
                send_means_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,            &
     &                         n_send_means)=1
                send_means_map(S_STRIDE_IN_SEND_ARRAY,                  &
     &                         n_send_means)=1
                send_means_map(S_ELEMENT_LENGTH,n_send_means)=1
                IF (fullfield) THEN ! we don't want halos
                  send_means_map(S_BASE_ADDRESS_IN_RECV_ARRAY,          &
     &                     n_send_means)=local_row-halosize(2,halo_type)
                ELSE ! halos are automatically removed
                  send_means_map(S_BASE_ADDRESS_IN_RECV_ARRAY,          &
     &                           n_send_means)=local_row
                ENDIF
                send_means_map(S_STRIDE_IN_RECV_ARRAY,                  &
     &                         n_send_means)=1

              ENDIF ! if at LHS of subarea

            ENDIF ! if mype  ==  dest_proc

          ENDDO ! row : loop over local rows on proc_id
        ENDDO ! proc_x : loop over processors in x dimension
      ENDDO ! proc_y : loop over processors in y dimension

! 3.0 Now the send and receive maps are set up, use
!     GCG_RALLTOALLE to redistribute the data
!     from the local_sum_arrays to the global_sum_arrays

      flag=GC_NONE ! flag argument is currently ignored by GCOM
      info=GC_NONE

      CALL GCG_RALLTOALLE(                                              &
     &  local_sum_array_top , send_data_map , n_send_data ,             &
     &  local_sum_array_len ,                                           &
     &  global_sum_array_top , rec_data_map , n_rec_data ,              &
     &  global_sum_array_len ,                                          &
     &  gc_all_proc_group , flag , info)

      info=GC_NONE

      CALL GCG_RALLTOALLE(                                              &
     &  local_sum_array_bot , send_data_map , n_send_data ,             &
     &  local_sum_array_len ,                                           &
     &  global_sum_array_bot , rec_data_map , n_rec_data ,              &
     &  global_sum_array_len ,                                          &
     &  gc_all_proc_group , flag , info)

! 4.0 Calculate mean of any zonal data on this processor

      DO j=1,n_rows_full_zonal_data

        zonal_sum_top=0.0
        zonal_sum_bot=0.0

        DO i=global_xstart,global_xend

          zonal_sum_top=zonal_sum_top+                                  &
     &                     global_sum_array_top(i,j)
          zonal_sum_bot=zonal_sum_bot+                                  &
     &                     global_sum_array_bot(i,j)
        ENDDO

        IF (zonal_sum_bot  ==  0.0) THEN
          zonal_mean_array(j)=rmdi
        ELSE
          zonal_mean_array(j)=zonal_sum_top/zonal_sum_bot
        ENDIF

      ENDDO

! 5.0 Send the calculated means back to the processors on the
!     LHS of the subarea, into the fieldout array

      flag=GC_NONE ! flag argument is currently ignored by GCOM
      info=GC_NONE

      CALL GCG_RALLTOALLE(                                              &
     &  zonal_mean_array , send_means_map , n_send_means ,              &
     &  global_sum_array_sizey,                                         &
     &  fieldout , rec_means_map , n_rec_means,                         &
     &  (yend-ystart)+1,                                                &
     &  gc_all_proc_group , flag , info)

!L
  999 CONTINUE
      RETURN
      END SUBROUTINE STZONM
