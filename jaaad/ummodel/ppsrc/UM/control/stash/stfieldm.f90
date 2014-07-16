
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: STFIELDM -------------------------------------------------
!LL
!LL  Purpose: Calculate weighted field mean within a region specified
!LL           by a lower left hand and upper right hand corner.
!LL           Single level fields only.
!LL           (STASH service routine).
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.Johns/S.Tett
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.3  16/09/93  Allow level-by-level mass-weighting if mass-weights
!LL                  are so defined, otherwise use P*.
!LL   4.3  28/01/97  Moved weighting and masking calculations up to
!LL                  SPATIAL.
!LL                  Significantly rewritten for mpp mode - data
!LL                  must be gathered to a processor for
!LL                  reproducible sums to be calculated.   P.Burton
!LL   4.4  13/06/97  mpp: Set fieldout to zero for processors in
!LL                  subdomain area which will not otherwise receive
!LL                  the result of the field mean.
!LL                  mpp: Correct bug in calculating SUMGBOT in non
!LL                       reproducible code                   P.Burton
!LL   4.5  12/01/98  Replaced usage of shmem common block by a
!LL                  dynamic array.                   P.Burton
!LL   4.5  09/01/98  Correct calculation of sum_pe      P.Burton
!LL   5.0  13/07/99  Changes for C-P C grid upgrade.
!LL                  R Rawlins.
!LL   5.0  17/11/99  Removed row_length,p_rows arguments
!LL                                             P.Burton
!LL   5.0  22/06/99  Added halo_type argument            P.Burton
!LL   5.0  15/9/99   Changes to South->North grid         P.Burton
!LL   5.1  15/05/00  S-N ordering consistency correction. R Rawlins
!     6.1  13/07/04  Add packing in 2 stage gather.
!                    P.Selwood/B. Carruthers (Cray)
!     6.2   18/10/05 Fix bugs introduced by gan1f601      Andy Malcolm
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D715
!LL
!LL  Project task: D7
!LL
!LL  External documentation:
!LL    Unified Model Doc Paper C4 - Storage handling and diagnostic
!LL                                 system (STASH)
!LL
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE STFIELDM(fieldin,vx,vy,fld_type,gr,halo_type,          &
     &                    lwrap,lmasswt,                                &
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
     &    fieldout,                                                     &
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
! Start i_stgfld

! Description:
!   This file contains an interface to STASH_GATHER_FIELD and
!   must be included whenever this routines is used so as to
!   get declarations of optional arguments correct.
!
! Current Code Owner: Paul Selwood
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   6.1  13/07/04   Original code. Paul Selwood.

      INTERFACE
        SUBROUTINE STASH_GATHER_FIELD (                                 &
     &    LOCAL_FIELD , GLOBAL_FIELD ,                                  &
     &    LOCAL_SIZE, GLOBAL_SIZE, LEVELS,                              &
     &    GLOBAL_NORTH , GLOBAL_EAST_IN , GLOBAL_SOUTH , GLOBAL_WEST,   &
     &    GRIDTYPE_CODE ,HALO_TYPE,                                     &
     &    GATHER_PE,                                                    &
     &    DATA_EXTRACTED,                                               &
     &    PACKING, IM_IDENT, LRLE, PACKING_TYPE,                        &
     &    NUM_OUT,                                                      &
     &    COMP_ACCRCY, loc_RMDI,                                        &
     &    ICODE, CMESSAGE)

        INTEGER, INTENT(IN) ::                                          &
     &    LOCAL_SIZE                                                    &
                          ! IN: size of level of LOCAL_FIELD
     &  , GLOBAL_SIZE                                                   &
                          ! IN: size of level of GLOBAL_FIELD
     &  , LEVELS                                                        &
                          ! IN: number of levels
     &  , GLOBAL_NORTH                                                  &
                          ! IN: specification of subdomain boundaries
     &  , GLOBAL_EAST_IN                                                &
                          ! IN: ""
     &  , GLOBAL_SOUTH                                                  &
                          ! IN: ""
     &  , GLOBAL_WEST                                                   &
                          ! IN: ""
     &  , GRIDTYPE_CODE                                                 &
                          ! IN: indicates the type of grid output
     &  , HALO_TYPE                                                     &
                          ! IN: type of halo on this field
     &  , GATHER_PE       ! IN: the PE to gather the global field to

        INTEGER, INTENT(OUT) ::                                         &
     &    ICODE           ! OUT: return code, 0=OK
!
! Optional Arguments to handle the COEX packing if necessary
!
        LOGICAL, INTENT(IN), OPTIONAL ::                                &
     &    PACKING                                                       &
                          ! IN: Set .true. if packing of the input
                          !     field is to be packed
     &  , LRLE            ! IN: True if Run Length Encoding is required

        INTEGER, INTENT(IN), OPTIONAL ::                                &
     &    IM_IDENT        ! IN: Internal model identifier

        INTEGER, INTENT(INOUT), OPTIONAL ::                             &
     &    PACKING_TYPE    ! IN/OUT: This flag is zero on input,
                          !         then stash packing is selected,
                          !         and the routine returns the
                          !         packing flag.
                          !
                          !         If the variable is set to 1 on input
                          !         then 32-bit packing for dumpfiles
                          !         is selected

        INTEGER, INTENT(OUT), OPTIONAL ::                               &
     &    NUM_OUT         ! OUT: Number of 32-bit IBM words in the
                          !      Packed field for WDGOS packing

        INTEGER, INTENT(IN), OPTIONAL ::                                &
     &    COMP_ACCRCY     ! IN: Packing Accuracy in Power of 2

        REAL, INTENT(IN), OPTIONAL ::                                   &
     &    loc_RMDI        ! IN: Missing data indicator
!
! Remaining Non-Optional Arguments
!
        LOGICAL, INTENT(IN) ::                                          &
     &    DATA_EXTRACTED  ! IN: TRUE if the data in LOCAL_FIELD has
                          !     already been extracted, or FALSE if
                          !     the extraction must be done here.

        REAL, INTENT(IN) ::                                             &
     &    LOCAL_FIELD(LOCAL_SIZE,LEVELS)
                          ! IN : local data

        REAL, INTENT(OUT) ::                                            &
     &    GLOBAL_FIELD(GLOBAL_SIZE,LEVELS)
                          ! OUT : (PE GATHER_PE only) - full gathered
                          !       field

        CHARACTER*(*), INTENT(OUT) ::                                   &
     &    CMESSAGE        ! OUT: Error message if ICODE .NE. 0

        END SUBROUTINE STASH_GATHER_FIELD
      END INTERFACE
! End i_stgfld
!
! Local variables
!
        INTEGER i,ii,j ! ARRAY INDICES FOR VARIABLE


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
! Co-ords to PE at top left of subarea
     &  proc_top_left_x , proc_top_left_y                               &

! unused return values from GLOBAL_TO_LOCAL_RC
     &, dummy1 , dummy2                                                 &

! PE number of PE at top left of subarea
     &, sum_pe                                                          &

! size of local and global arrays
     &, local_size,global_size

! Weighted version of fieldin
      REAL local_sum_array_top(xstart:xend,ystart:yend)
! Weights applied to fieldin
      REAL local_sum_array_bot(xstart:xend,ystart:yend)



      INTEGER                                                           &
! Sizes of the global_sum_arrays defined below
     &  global_sum_array_sizex,global_sum_array_sizey

      REAL                                                              &
! Collected versions of fieldin and the weights containing
! whole (subarea) columns of meridional data
     &  global_sum_array_top(global_xstart:global_xend,                 &
     &                       global_ystart:global_yend)                 &
     &, global_sum_array_bot(global_xstart:global_xend,                 &
     &                       global_ystart:global_yend)


        REAL SUMFTOP
        REAL SUMFBOT

!L----------------------------------------------------------------------
!L 0. Initialise sums
!L
      SUMFTOP=0.0
      SUMFBOT=0.0
!L----------------------------------------------------------------------

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
     &            pstar_weight(ii,j)*area_weight(ii,j)
                local_sum_array_top(i,j)=                               &
     &            fieldin(ii,j)*pstar_weight(ii,j)*area_weight(ii,j)
            ELSE
              local_sum_array_bot(i,j)=0.0
              local_sum_array_top(i,j)=0.0
            ENDIF
          ENDDO
        ENDDO

      ENDIF  ! if this processor contains any of the subarea

! Initialise fieldout - so all PE's have valid data
! (Only PEs on top left of subdomain get the field mean)

      fieldout=0.0


! The local_sum_arrays must be distributed so that the complete
! sub-area exists on a single processor, so that a reproducible sum
! can be carried out.

! 0.0 : Initialise variables defining the size of the arrays
!       global_sum_arrays

      global_sum_array_sizex=global_xend-global_xstart+1
      global_sum_array_sizey=global_yend-global_ystart+1

! 1.0 Gather the fields to a single processor

! DEPENDS ON: global_to_local_rc
      CALL GLOBAL_TO_LOCAL_RC(gr,halo_type,                             &
     &  global_xstart , global_ystart,                                  &
     &  proc_top_left_x, proc_top_left_y,                               &
     &  dummy1,dummy2)

      sum_pe=proc_top_left_x + nproc_x*proc_top_left_y

      local_size=(xend-xstart+1)*(yend-ystart+1)
      global_size=global_sum_array_sizex*global_sum_array_sizey

! DEPENDS ON: stash_gather_field
      CALL STASH_GATHER_FIELD (                                         &
     &  local_sum_array_top , global_sum_array_top ,                    &
     &  local_size          , global_size,                              &
     &  1,                                                              &
            ! 1 level
     &  global_yend, global_xend, global_ystart, global_xstart,         &
     &  gr , halo_type, sum_pe,                                         &
     &  .TRUE.,                                                         &
                ! data has been extracted
     &  ICODE=ICODE, CMESSAGE=CMESSAGE)

      IF (ICODE  /=  0) THEN
        WRITE(6,*) 'STFIELDM : MPP Error in STASH_GATHER_FIELD'
        WRITE(6,*) CMESSAGE
        GOTO 999
      ENDIF

! DEPENDS ON: stash_gather_field
      CALL STASH_GATHER_FIELD (                                         &
     &  local_sum_array_bot , global_sum_array_bot ,                    &
     &  local_size          , global_size,                              &
     &  1,                                                              &
            ! 1 level
     &  global_yend, global_xend, global_ystart, global_xstart,         &
     &  gr , halo_type, sum_pe,                                         &
     &  .TRUE.,                                                         &
                ! data has been extracted
     &  ICODE=ICODE, CMESSAGE=CMESSAGE)

      IF (ICODE  /=  0) THEN
        WRITE(6,*) 'STFIELDM : MPP Error in STASH_GATHER_FIELD'
        WRITE(6,*) CMESSAGE
        GOTO 999
      ENDIF

! 2.0 Calculate the sums

      IF (mype  ==  sum_pe) THEN

        DO i=global_xstart,global_xend
          DO j=global_ystart,global_yend
            SUMFTOP=SUMFTOP+global_sum_array_top(i,j)
            SUMFBOT=SUMFBOT+global_sum_array_bot(i,j)
          ENDDO
        ENDDO

        IF (SUMFBOT  ==  0.0) THEN
          fieldout=rmdi
        ELSE
         fieldout=SUMFTOP/SUMFBOT
        ENDIF

      ENDIF

  999 CONTINUE
      RETURN
      END SUBROUTINE STFIELDM
