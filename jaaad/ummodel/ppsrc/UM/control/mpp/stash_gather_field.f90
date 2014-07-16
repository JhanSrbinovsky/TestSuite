

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Gathers STASHed data from many processors to one processor
!
! Subroutine interface:
      SUBROUTINE STASH_GATHER_FIELD (                                   &
     &  LOCAL_FIELD , GLOBAL_FIELD ,                                    &
     &  LOCAL_SIZE, GLOBAL_SIZE, LEVELS,                                &
     &  GLOBAL_NORTH , GLOBAL_EAST_IN , GLOBAL_SOUTH , GLOBAL_WEST,     &
     &  GRIDTYPE_CODE ,HALO_TYPE,                                       &
     &  GATHER_PE,                                                      &
     &  DATA_EXTRACTED,                                                 &
     &  PACKING, IM_IDENT, LRLE, PACKING_TYPE,                          &
     &  NUM_OUT,                                                        &
     &  COMP_ACCRCY, loc_RMDI,                                          &
     &  ICODE, CMESSAGE)

      IMPLICIT NONE

! Description:
! Takes a decomposed, STASH processed field and gathers
! it to a single processor, ready for I/O,
!
! Method:
! See in-line documentation
!
! Current code owner : P.Selwood
!
! Subroutine arguments:


      INTEGER, INTENT(IN) ::                                            &
     &  LOCAL_SIZE                                                      &
                        ! IN: size of level of LOCAL_FIELD
     &, GLOBAL_SIZE                                                     &
                        ! IN: size of level of GLOBAL_FIELD
     &, LEVELS                                                          &
                        ! IN: number of levels
     &, GLOBAL_NORTH                                                    &
                        ! IN: specification of subdomain boundaries
     &, GLOBAL_EAST_IN                                                  &
                        ! IN: ""
     &, GLOBAL_SOUTH                                                    &
                        ! IN: ""
     &, GLOBAL_WEST                                                     &
                        ! IN: ""
     &, GRIDTYPE_CODE                                                   &
                        ! IN: indicates the type of grid output
     &, HALO_TYPE                                                       &
                        ! IN: type of halo on this field
     &, GATHER_PE       ! IN: the PE to gather the global field to

      INTEGER, INTENT(OUT) ::                                           &
     &  ICODE           ! OUT: return code, 0=OK
!
! Optional Arguments to handle the COEX packing if necessary
!
      LOGICAL, INTENT(IN), OPTIONAL ::                                  &
     &  PACKING                                                         &
                        ! IN: Set .true. if packing of the input
                        !     field is to be packed
     &, LRLE            ! IN: True if Run Length Encoding is required

      INTEGER, INTENT(IN), OPTIONAL ::                                  &
     &  IM_IDENT        ! IN: Internal model identifier

      INTEGER, INTENT(INOUT), OPTIONAL ::                               &
     &  PACKING_TYPE    ! IN/OUT: This flag is zero on input,
                        !         then stash packing is selected,
                        !         and the routine returns the
                        !         packing flag.
                        !
                        !         If the variable is set to 1 on input
                        !         then 32-bit packing for dumpfiles
                        !         is selected

      INTEGER, INTENT(OUT), OPTIONAL ::                                 &
     &  NUM_OUT         ! OUT: Number of 32-bit IBM words in the Packed
                        !      field for WDGOS packing

      INTEGER, INTENT(IN), OPTIONAL ::                                  &
     &  COMP_ACCRCY     ! IN: Packing Accuracy in Power of 2

      REAL, INTENT(IN), OPTIONAL ::                                     &
     &  loc_RMDI        ! IN: Missing data indicator
!
! Remaining Non-Optional Arguments
!
      LOGICAL, INTENT(IN) ::                                            &
     &  DATA_EXTRACTED  ! IN: TRUE if the data in LOCAL_FIELD has
                        !     already been extracted, or FALSE if
                        !     the extraction must be done here.

      REAL, INTENT(IN) ::                                               &
     &  LOCAL_FIELD(LOCAL_SIZE,LEVELS)
                        ! IN : local data

      REAL, INTENT(OUT) ::                                              &
     &  GLOBAL_FIELD(GLOBAL_SIZE,LEVELS)
                        ! OUT : (PE GATHER_PE only) - full gathered
                        !       field

      CHARACTER*(*), INTENT(OUT) ::                                     &
     &  CMESSAGE        ! OUT: Error message if ICODE  /=  0

! Parameters and common blocks
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

! Local variables

      INTEGER                                                           &
     &  GLOBAL_EAST                                                     &
                      ! copy of GLOBAL_EAST_IN with wrap around s.t.
!                     ! GLOBAL_EAST > GLOBAL_ROW_LEN
     &, global_x                                                        &
                      ! size of global data EW
     &, global_y                                                        &
                      ! size of global data NS
     &, fld_type                                                        &
                      ! indicates if field is on P or U grid
     &, level                                                           &
                      ! loop index for loop over levels
     &, i                                                               &
                      ! loop counter
     &, proc_topleft_x,proc_topleft_y                                   &
                                        ! processors at corners of
     &, proc_botright_x,proc_botright_y                                 &
                                        ! the subarea
     &, dummy1,dummy2                                                   &
                      ! ignored return arguments
     &, procx,procy                                                     &
                      ! loop indexes for loops over processors
     &, eff_procx                                                       &
                      ! real x co-ord of processor column procx
     &, procid                                                          &
                      ! processor id of (procx,procy)
     &, local_xstart,local_xend                                         &
                                 ! boundaries of subdomain for
     &, local_ystart,local_yend                                         &
                                 ! processor procid
     &, local_start_row                                                 &
                          ! first row to send from procid
     &, local_start_col                                                 &
                          ! first column to send from procid
     &, sendsize_x                                                      &
                          ! number of points on each row to send
     &, nrows_to_send                                                   &
                          ! number of rows to send from procid
     &, local_row_length                                                &
                          ! size of sending array EW
     &, global_start_row                                                &
                          ! first row to receive at on GATHER_PE
     &, global_start_col                                                &
                          ! first col. to recieve on GATHER_PE
     &, global_row_length                                               &
                          ! size of receiving array EW
     &, flag,info         ! GCOM arguments

! Copies of arguments / variables used to decide if we can use the
! send/receive maps used in the last call

      INTEGER                                                           &
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE                                &
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN                           &
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST                              &
     &, old_GRIDTYPE_CODE , old_GATHER_PE                               &
     &, old_current_decomp_type, old_HALO_TYPE

      INTEGER                                                           &
! variables defining send and receive maps to be passed to
! GC_RALL_TO_ALL, defining the data transposition
     &  send_map(7,2)                                                   &
     &, receive_map(7,2*MAXPROC)                                        &
     &, n_sends,n_recvs  ! number of sends and receives


      LOGICAL                                                           &
     &  wrap                                                            &
              ! if the subdomain wraps over 0 degree meridion
     &, wrap_special                                                    &
                     ! if there is a wrap around, which starts and
!                      ends on the same processor
     &, zonal_data                                                      &
                     ! if this is a zonal data grid
     &, fullfield                                                       &
                     ! if this is a full field - NOT a subarea
     &, l_packing    ! if packing of data is required

! Save all the variables that may be used in the next call
      SAVE                                                              &
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE                                &
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN                           &
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST                              &
     &, old_GRIDTYPE_CODE , old_GATHER_PE                               &
     &, old_current_decomp_type                                         &
     &, send_map,receive_map,n_sends,n_recvs                            &
     &, old_HALO_TYPE

! Set all the old_* variables to a number indicating they've
! not been used yet

      DATA                                                              &
     &  old_LOCAL_SIZE , old_GLOBAL_SIZE                                &
     &, old_GLOBAL_NORTH , old_GLOBAL_EAST_IN                           &
     &, old_GLOBAL_SOUTH , old_GLOBAL_WEST                              &
     &, old_GRIDTYPE_CODE , old_GATHER_PE                               &
     &, old_current_decomp_type, old_HALO_TYPE                          &
     &  / -1,-1,-1,-1,-1,-1,-1,-1,-1,-1 /

! Functions

      INTEGER GET_FLD_TYPE
! ------------------------------------------------------------------

      ICODE=0
      IF (PRESENT(PACKING)) THEN
        L_PACKING = PACKING
      ELSE
        L_PACKING = .FALSE.
      END IF

! DEPENDS ON: get_fld_type
       fld_type=GET_FLD_TYPE(GRIDTYPE_CODE)
! See if there is wrap around over meridion, and if so make
! sure that GLOBAL_EAST is > glsize(1)

      GLOBAL_EAST=GLOBAL_EAST_IN
       IF (GLOBAL_EAST  >   glsize(1,fld_type)) THEN
        wrap=.TRUE.
      ELSEIF (GLOBAL_EAST  <   GLOBAL_WEST) THEN
        wrap=.TRUE.
        GLOBAL_EAST=GLOBAL_EAST_IN+glsize(1,fld_type)
      ELSE
        wrap=.FALSE.
      ENDIF

      IF ((GRIDTYPE_CODE  ==  ppx_atm_tzonal) .OR.                      &
                                                   ! Atmos T zonal
     &   ( GRIDTYPE_CODE  ==  ppx_atm_uzonal) .OR.                      &
                                                   ! Atmos U zonal
     &   ( GRIDTYPE_CODE  ==  ppx_ocn_tzonal) .OR.                      &
                                                   ! Ocean T zonal
     &   ( GRIDTYPE_CODE  ==  ppx_ocn_uzonal))                          &
                                                   ! Atmos U zonal
     &  THEN

! This is a zonal field

        zonal_data=.TRUE.
        global_x=1

        IF ((GRIDTYPE_CODE  ==  ppx_atm_tzonal) .OR.                    &
                                                     ! Atmos T zonal
     &      ( GRIDTYPE_CODE  ==  ppx_ocn_tzonal))                       &
                                                     ! Ocean T zonal
     &  THEN
          fld_type=fld_type_p
        ELSE
          fld_type=fld_type_u
        ENDIF
      ELSE

! This is a normal field

        zonal_data=.FALSE.
        global_x=glsize(1,fld_type)
      ENDIF

      global_y=glsize(2,fld_type)

! Set up logical indicating if this is a full field, or just
! a subdomain

      IF (zonal_data) THEN

        fullfield= ( ( GLOBAL_SOUTH  ==  1) .AND.                       &
     &             ( GLOBAL_NORTH  ==  global_y))

      ELSE

        fullfield = (( GLOBAL_WEST  ==  1) .AND.                        &
     &               ( GLOBAL_EAST  ==  global_x) .AND.                 &
     &               ( GLOBAL_SOUTH  ==  1) .AND.                       &
     &               ( GLOBAL_NORTH  ==  global_y))

      ENDIF

! Dealing with fields not in model grid

      if((global_x == 0).or.(global_x == imdi)) then
        write(6,*)'local_size=',local_size
        write(6,*)'global_size=',global_size
        do level=1,levels
          do i=1,global_size
            global_field(i,level)=local_field(i,level)
          enddo
        enddo
      else

! If this a fullfield, we can simply use the standard
! GATHER_FIELD routine

      IF (fullfield) THEN

        IF (zonal_data) THEN

! DEPENDS ON: gather_zonal_field
          CALL GATHER_ZONAL_FIELD( LOCAL_FIELD,GLOBAL_FIELD,            &
     &                          lasize(2,fld_type,halo_type),global_y,  &
     &                            LEVELS,GRIDTYPE_CODE,fld_type,        &
     &                            halo_type,GATHER_PE)

        ELSE

          DO level=1,LEVELS

            IF (L_PACKING) THEN
! DEPENDS ON: gather_pack_field
              CALL GATHER_PACK_FIELD(                                   &
     &                        LOCAL_FIELD(1,level),                     &
     &                        GLOBAL_FIELD(1,level),                    &
     &                        lasize(1,fld_type,halo_type),             &
     &                        lasize(2,fld_type,halo_type),             &
     &                        global_x,global_y,                        &
     &                        fld_type,halo_type,                       &
     &                        GATHER_PE,GC_ALL_PROC_GROUP,              &
     &                        PACKING, IM_IDENT, LRLE, PACKING_TYPE,    &
     &                        NUM_OUT,                                  &
     &                        COMP_ACCRCY, loc_RMDI)
            ELSE
! DEPENDS ON: gather_field
              CALL GATHER_FIELD( LOCAL_FIELD(1,level),                  &
     &                           GLOBAL_FIELD(1,level),                 &
     &                           lasize(1,fld_type,halo_type),          &
     &                           lasize(2,fld_type,halo_type),          &
     &                           global_x, global_y,                    &
     &                           fld_type, halo_type,                   &
     &                           GATHER_PE, GC_ALL_PROC_GROUP,          &
     &                           ICODE, CMESSAGE)
            END IF

            IF (ICODE  /=  0) THEN
              WRITE(6,*)                                                &
     &         'STASH_GATHER_FIELD : Failed during GATHER_FIELD'
              WRITE(6,*) 'Error code : ',ICODE
              WRITE(6,*) 'Message : ',CMESSAGE

              ICODE=2
              CMESSAGE='Failed to gather field'
              GOTO 9999
            ENDIF

           ENDDO

         ENDIF
       ELSE
! for subdomains, life is not so easy - we must explicitly
! calculate our own send and receive maps, and use GCG_RALLTOALLE
! to shift the data around.

! If the same arguments are used as were used in the last call
! to this routine, we can just use the previously calculated
! send and receive maps, otherwise we need to calculate new maps

        IF (.NOT. (                                                     &
     &    (LOCAL_SIZE  ==  old_LOCAL_SIZE) .AND.                        &
     &    (GLOBAL_SIZE  ==  old_GLOBAL_SIZE) .AND.                      &
     &    (GLOBAL_NORTH  ==  old_GLOBAL_NORTH) .AND.                    &
     &    (GLOBAL_EAST_IN  ==  old_GLOBAL_EAST_IN) .AND.                &
     &    (GLOBAL_SOUTH  ==  old_GLOBAL_SOUTH) .AND.                    &
     &    (GLOBAL_WEST  ==  old_GLOBAL_WEST) .AND.                      &
     &    (GRIDTYPE_CODE  ==  old_GRIDTYPE_CODE) .AND.                  &
     &    (GATHER_PE  ==  old_GATHER_PE) .AND.                          &
     &    (HALO_TYPE  ==  old_HALO_TYPE) .AND.                          &
     &    (current_decomp_type  ==  old_current_decomp_type ))) THEN

          old_LOCAL_SIZE=LOCAL_SIZE
          old_GLOBAL_SIZE=GLOBAL_SIZE
          old_GLOBAL_NORTH=GLOBAL_NORTH
          old_GLOBAL_EAST_IN=GLOBAL_EAST_IN
          old_GLOBAL_SOUTH=GLOBAL_SOUTH
          old_GLOBAL_WEST=GLOBAL_WEST
          old_GRIDTYPE_CODE=GRIDTYPE_CODE
          old_GATHER_PE=GATHER_PE
          old_current_decomp_type=current_decomp_type
          old_HALO_TYPE=HALO_TYPE

! Find out what the boundaries of the subdomain area

! DEPENDS ON: global_to_local_rc
          CALL GLOBAL_TO_LOCAL_RC(GRIDTYPE_CODE,HALO_TYPE,              &
     &                            GLOBAL_WEST,GLOBAL_NORTH,             &
     &                            proc_topleft_x,proc_topleft_y,        &
     &                            dummy1,dummy2)
! DEPENDS ON: global_to_local_rc
          CALL GLOBAL_TO_LOCAL_RC(GRIDTYPE_CODE,HALO_TYPE,              &
     &                            GLOBAL_EAST,GLOBAL_SOUTH,             &
     &                            proc_botright_x,proc_botright_y,      &
     &                            dummy1,dummy2)

! Ensure that the processor x co-ords are such that the botright_x is
! always greater than (or equal to) top_left_x.
          IF (wrap) proc_botright_x=gridsize(1)+proc_botright_x

! wrap_special is set to true if there is a wrap around which starts
! and ends on the same processor. This case requires extra work as
! the processor in question
          IF (wrap .AND. (proc_topleft_x+gridsize(1)  ==                &
     &                    proc_botright_x)) THEN
            wrap_special=.TRUE.
          ELSE
            wrap_special=.FALSE.
          ENDIF

          n_sends=0
          n_recvs=0

          DO procy=proc_botright_y,proc_topleft_y
            DO procx=proc_topleft_x,proc_botright_x

              eff_procx=MOD(procx,gridsize(1))
              procid=eff_procx+procy*gridsize(1)

! DEPENDS ON: global_to_local_subdomain
              CALL GLOBAL_TO_LOCAL_SUBDOMAIN(                           &
     &          .TRUE.,.TRUE.,                                          &
     &          GRIDTYPE_CODE,HALO_TYPE,procid,                         &
     &          GLOBAL_SOUTH,GLOBAL_EAST,                               &
     &          GLOBAL_NORTH,GLOBAL_WEST,                               &
     &          local_ystart,local_xend,                                &
     &          local_yend  ,local_xstart)

! Calculate the shape of the arrays, and where to start sending/
! receiving data, and how many rows to send

              IF (DATA_EXTRACTED) THEN
                local_start_row=1
              ELSE
                local_start_row=local_ystart
              ENDIF
              nrows_to_send=local_yend-local_ystart+1

             global_start_row=g_datastart(2,procid)+local_ystart-       &
     &                         halosize(2,halo_type) - GLOBAL_SOUTH
              global_row_length=GLOBAL_EAST-GLOBAL_WEST+1

! Calculate the following variables:
! local_row_length : the X dimension size of the local array
! local_send_offx  : the offset into each row to start sending from
! sendsize_x       : the number of points on each row to send
! The calculation of these numbers is different for processors
! at the start and end of a wrap_special case
! Note that when DATA_EXTRACTED is true, then local_field has no
! halos.

              IF (wrap_special .AND. procx  ==  proc_topleft_x) THEN
                IF (DATA_EXTRACTED) THEN
      local_row_length = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
     &                   + local_xend - local_xstart + 1
      sendsize_x       = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
     &                   - local_xstart + 1
      local_start_col  = 1
                ELSE
      local_row_length = g_lasize(1,fld_type,halo_type,procid)
      sendsize_x       = g_lasize(1,fld_type,halo_type,procid)          &
     &                   - local_xstart
      local_start_col  = local_xstart

                ENDIF
                global_start_col=1

              ELSEIF (wrap_special .AND. procx  ==  proc_botright_x)    &
     &        THEN
                IF (DATA_EXTRACTED) THEN
      local_row_length = g_lasize(1,fld_type,halo_type_no_halo,procid)  &
     &                   + local_xend - local_xstart + 1
      local_start_col  = local_row_length - local_xend + 1
      sendsize_x       = local_xend
                ELSE
      local_row_length = g_lasize(1,fld_type,halo_type,procid)
      local_start_col  = Offx + 1
      sendsize_x       = local_xend - Offx
                ENDIF
                global_start_col=global_row_length-sendsize_x+1

              ELSE
                IF (DATA_EXTRACTED) THEN
                  local_row_length=local_xend-local_xstart+1
                  local_start_col=1
                ELSE
          local_row_length=g_lasize(1,fld_type,halo_type,procid)
                  local_start_col=local_xstart
                ENDIF
                sendsize_x=local_xend-local_xstart+1
                global_start_col=local_xstart-halosize(1,halo_type)+    &
     &                           g_datastart(1,procid)-GLOBAL_WEST
              ENDIF

              IF (global_start_col  <   0) THEN
! Wrapped around field, but this processor is not start or end
! processor
        global_start_col=global_start_col+glsize(1,fld_type)
              ENDIF

! Now we can set up the send and receive map entries for the data on
! this processor

              IF (mype  ==  procid) THEN  ! I need to send some data


                n_sends=n_sends+1

                send_map(S_DESTINATION_PE,n_sends) = GATHER_PE
                send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_sends) =        &
     &            (local_start_row-1)*local_row_length +                &
     &            local_start_col
                send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_sends) =        &
     &            nrows_to_send
                send_map(S_STRIDE_IN_SEND_ARRAY,n_sends) =              &
     &            local_row_length
                send_map(S_ELEMENT_LENGTH,n_sends) = sendsize_x
                send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_sends) =        &
     &            (global_start_row-1)*global_row_length +              &
     &            global_start_col
                send_map(S_STRIDE_IN_RECV_ARRAY,n_sends) =              &
     &            global_row_length

              ENDIF ! if I'm sending data

              IF (mype  ==  GATHER_PE) THEN  ! I need to receive data

                n_recvs=n_recvs+1

                receive_map(R_SOURCE_PE,n_recvs) = procid
                receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recvs) =     &
     &            (global_start_row-1)*global_row_length +              &
     &            global_start_col
                receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recvs) =     &
     &            nrows_to_send
                receive_map(R_STRIDE_IN_RECV_ARRAY,n_recvs) =           &
     &            global_row_length
                receive_map(R_ELEMENT_LENGTH,n_recvs) = sendsize_x
                receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recvs) =     &
     &            (local_start_row-1)*local_row_length +                &
     &            local_start_col
                receive_map(R_STRIDE_IN_SEND_ARRAY,n_recvs) =           &
     &            local_row_length

              ENDIF ! if I'm receiving data

            ENDDO ! procx : loop along processor row

          ENDDO ! procy : loop down processor column

        ENDIF ! if I need to recalculate my send/receive maps

! Send / receive the data using GCG_RALLTOALLE


        DO level=1,LEVELS

          flag=0  ! This is currently ignored at GCG v1.1
          info=GC_NONE

          CALL GCG_RALLTOALLE(                                          &
     &      LOCAL_FIELD(1,level)  ,                                     &
     &      send_map    , n_sends  ,LOCAL_SIZE  ,                       &
     &      GLOBAL_FIELD(1,level) ,                                     &
     &      receive_map , n_recvs , GLOBAL_SIZE ,                       &
     &      GC_ALL_PROC_GROUP , flag, info)

        ENDDO

        ENDIF ! if this is a full or extracted field

      endif


 9999 CONTINUE

      RETURN
      END SUBROUTINE STASH_GATHER_FIELD

