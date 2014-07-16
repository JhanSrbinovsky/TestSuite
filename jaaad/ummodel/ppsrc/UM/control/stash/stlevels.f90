
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine STLEVELS -----------------------------------------------
!LL
!LL  Purpose: Generate a level index from STASHrecord and level_lists
!LL           and number of levels tailored to a particular diagnostic.
!LL           Also set levels and pseudo-levels information for encoding
!LL           PPheader details.  (This subroutine based on a merger
!LL           between GEN_INDEX and PP_COMPUTE_LEVEL).
!LL                  New subroutine STLEVELS is based on GEN_INDEX and
!LL                  PP_COMPUTE_LEVEL with merged functionality.
!LL           A general note as levels list is an integer
!LL           real values are multiplied by a 1000.0.
!LL           When computing the real value of the level for the
!LL           pp header it is necessary to divide by a 1000.0.
!LL           Levels that are affected by this are theta, pressure and
!LL           height. S. Anderson.
!LL
!LL  Author:   T.Johns
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  14/01/93  Set pseudo_level to 0 not IMDI if no pseudo-level.
!LL        29/01/93  Correct FORMAT statements.
!LL  3.1     14/01/93 Include PV levels as levels divided by 1000.
!LL   3.2  19/04/93  Correct roundoff error for LEVEL type conversion.
!LL                  1.0E-10 is added after REAL divide by 1000 (TCJ).
!LL  4.0  14/12/95  Correct long-standing error in input levels range
!LL                 to output levels list conversion.  RTHBarnes.
!    4.4  02/12/96 Time mean timeseries added R A Stratton.
!    5.0  10/06/99 Remove references to ak,bk,ak_lev,bk_lev, which are
!                  not used in downstream processing. Rick Rawlins
!    5.3  23/07/01 Replace lbvc magic no.s by parameters. R Rawlins
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered : C4?
!LL
!LL  Project task: C4
!LL
!LL  External documentation : UMDP no C4
!LL
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE STLEVELS(stash_control,stash_control_size,             &
     &     stash_levels,num_stash_levels,num_level_lists,               &
     &     stash_pseudo_levels,num_stash_pseudo,num_pseudo_lists,       &
     &     max_stash_levs,num_levs_in,num_levs_out,index_size,          &
     &     index_lev,level_list,                                        &
     &     lbvcl,level,pseudo_level,                                    &
     &     icode,cmessage)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &       stash_control_size,                                        &
                                 ! IN size of stash control record
     &       stash_control(stash_control_size),                         &
                                               ! IN  stash control
     &       num_stash_levels,                                          &
                                 ! IN max. no of hts for a levels list
     &       num_level_lists,                                           &
                                 ! IN max. no of level lists
     &       stash_levels(num_stash_levels+1,num_level_lists),          &
                                                               ! IN
!                                !    lookup table for level lists
     &       num_stash_pseudo,num_pseudo_lists,                         &
                                               ! IN dims of pseudo_levs
     &       stash_pseudo_levels(num_stash_pseudo+1,num_pseudo_lists),  &
!                                ! IN lookup table for pseudo-lev lists
     &       max_stash_levs,                                            &
                                 ! IN max. no of output levels
     &       num_levs_in,                                               &
                                 ! OUT no of levels in input data
     &       num_levs_out,                                              &
                                 ! OUT no of levels in output data
     &       index_size,                                                &
                                 ! OUT no of levels in levels index
     &       index_lev(max_stash_levs),                                 &
                                        ! OUT index of output level
!                                               relative to input level
     &       level_list(max_stash_levs),                                &
                                         ! OUT value of model level
     &       pseudo_level(max_stash_levs),                              &
                                           ! OUT Value of pseudo levels
     &       lbvcl,                                                     &
                                 ! IN  vertical coordinate PP code
     &       icode               ! OUT error code
      REAL                                                              &
     &       level(max_stash_levs)  ! OUT Value of output levels (real)
      CHARACTER*(*)                                                     &
     &       cmessage            ! OUT error message
!*----------------------------------------------------------------------
! Parameters
!
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
!
! Local variables
!
      INTEGER                                                           &
     &       index_pseudo_lev(max_stash_levs),                          &
                                               ! Pseudo-level 1D index
     &       num_pseudo_in,num_pseudo_out,                              &
                                               ! Number of pseudo levs
     &       k2,ml,kl,                                                  &
                                       ! loop counts
     &       NI,NO,                                                     &
                                       ! Number In/Out
     &       indx1,                                                     &
                                       ! index count
     &       ilev,                                                      &
                                       ! Integer level/pseudo-level
     &       what_mean,what_proc       ! Meaning and processing code
!
! First compute the index for physical levels
!
      IF(STASH_CONTROL(st_input_bottom) <  0) THEN ! Input LEVELS list
        NI=-STASH_CONTROL(st_input_bottom)
        NUM_LEVS_IN=STASH_LEVELS(1,NI)
        IF(STASH_CONTROL(st_output_bottom) <  0) THEN ! LEVELS LIST out
          NO=-STASH_CONTROL(st_output_bottom)
          NUM_LEVS_OUT=STASH_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_LEVS_OUT
            ilev=STASH_LEVELS(ML+1,NO)    !  Level required
            DO KL=1,NUM_LEVS_IN
              IF(STASH_LEVELS(KL+1,NI) == ilev) THEN
                INDX1=INDX1+1
                INDEX_LEV(INDX1)=KL   ! Relative position of Input to Ou
                level_list(indx1)=ilev
                GOTO 400
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output level ',ilev,                   &
     &                          ' not found in input levels list'
            GOTO 999
 400        CONTINUE
            ENDDO
        ELSE           !  Output as a Level range
          NUM_LEVS_OUT=STASH_CONTROL(st_output_top)-                    &
     &                 STASH_CONTROL(st_output_bottom)+1
          ilev=STASH_CONTROL(st_output_bottom) !1st output model level
          DO KL=1,NUM_LEVS_IN
            IF(STASH_LEVELS(KL+1,NI) == ilev) THEN
              INDEX_LEV(1)=KL ! Relative posn of Input to the 1st level
              level_list(1)=ilev
              GOTO 401
            ENDIF
          ENDDO
          ICODE=nonsense
          WRITE(CMESSAGE,101) 'Output bottom model level ',ilev,        &
     &                        ' not found in input levels list'
          GOTO 999
 401      CONTINUE
          DO KL=2,NUM_LEVS_OUT
            INDEX_LEV(KL)=INDEX_LEV(KL-1)+1
            level_list(kl)=level_list(kl-1)+1
          ENDDO
        ENDIF
      ELSEIF(STASH_CONTROL(st_input_bottom) == 100) THEN !Special level
          NUM_LEVS_IN=1
          NUM_LEVS_OUT=1
          INDEX_LEV(1)=1
          level_list(1)=1 ! could be worth setting to some nonsense no.
      ELSE     !  Input is Model level range
        NUM_LEVS_IN=STASH_CONTROL(st_input_top)-                        &
     &              STASH_CONTROL(st_input_bottom)+1
        IF(STASH_CONTROL(st_output_bottom) <  0) THEN ! LEVELS LIST out
          NO=-STASH_CONTROL(st_output_bottom)
          NUM_LEVS_OUT=STASH_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_LEVS_OUT
            ilev=STASH_LEVELS(ML+1,NO)    ! Output level reqd
            DO KL=1,NUM_LEVS_IN
              IF((STASH_CONTROL(st_input_bottom)+KL-1) == ilev) THEN
                INDX1=INDX1+1
                INDEX_LEV(INDX1)=KL   ! Relative posn of output to inpt
                level_list(INDX1)=ilev
                GOTO 402
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output model level ',ilev,             &
     &                          ' not in input model level range'
            GOTO 999
 402        CONTINUE
          ENDDO
        ELSE     !   Output as model level range
! Do some consistency checks here to ensure valid processing request
! output bottom should be greater or equal to input bottom
          IF (stash_control(st_output_bottom) <                         &
     &       stash_control(st_input_bottom)) THEN
            icode=nonsense
            write(cmessage,103)'bad level spec, bot input>output',      &
     &       stash_control(st_input_bottom),                            &
     &       stash_control(st_output_bottom)
            goto 999 ! jump to error
          ELSEIF (stash_control(st_output_top) >                        &
     &         stash_control(st_input_top)) THEN
            icode=nonsense
            write(cmessage,103)'bad level spec, top input<output',      &
     &        stash_control(st_input_top),                              &
     &        stash_control(st_output_top)
              goto 999 ! jump to error
          ENDIF
          NUM_LEVS_OUT=STASH_CONTROL(st_output_top)-                    &
     &                 STASH_CONTROL(st_output_bottom)+1
          INDEX_LEV(1)=STASH_CONTROL(st_output_bottom)-                 &
     &                 STASH_CONTROL(st_input_bottom)+1
          level_list(1)=stash_control(st_output_bottom)
          DO kl=2,NUM_LEVS_OUT
            INDEX_LEV(kl)=INDEX_LEV(kl-1)+1
            level_list(kl)=level_list(kl-1)+1
          ENDDO
        ENDIF
      ENDIF
      index_size=num_levs_out
      IF (num_levs_out >  num_levs_in) THEN   ! things very badly wrong
        icode=nonsense
        write(cmessage,103)'asking for num_levs_out>num_levs_in',       &
     &   num_levs_out,num_levs_in
        goto 999 ! jump to return
      ENDIF
!
! Next, compute actual (physical) levels for encoding PPheaders
!
      IF (STASH_CONTROL(st_output_bottom) <  0) THEN ! Levels List ?
        NO=-STASH_CONTROL(st_output_bottom)     ! Index of Levels list

          ! Remove scaling (by factor 1000) of vertical level coord
          ! for certain types of STASH output [originally needed to
          ! store in an intermediary integer array]
          IF( LBVCL  ==  ppx_lbvc_height   .OR.                         &
                                                !  height levels
     &        LBVCL  ==  ppx_lbvc_pressure .OR.                         &
                                                ! pressure levels
     &        LBVCL  ==  ppx_lbvc_theta    .OR.                         &
                                                ! theta levels
     &        LBVCL  ==  ppx_lbvc_PV ) THEN     ! potential vorticity


          DO ML=1,NUM_LEVS_OUT
            LEVEL(ML)=REAL(STASH_LEVELS(ML+1,NO))*0.001+1.0E-10
          ENDDO
        ELSE
          DO ML=1,NUM_LEVS_OUT
            LEVEL(ML)=REAL(STASH_LEVELS(ML+1,NO))
          ENDDO
        ENDIF
      ELSEIF (STASH_CONTROL(st_output_bottom) == st_special_code) THEN
       ! Special level.
       ! The LEVEL array is not used by the model except to construct pp
       ! header items at output. The value of -1.0 is set as a flag for
       ! special levels so that routine PP_HEAD will insert the lbvc
       ! item in STASHmaster record.
        DO ML=1,NUM_LEVS_OUT
          LEVEL(ML)=-1.0
        ENDDO
      ELSE
        DO ML=1,NUM_LEVS_OUT
          LEVEL(ML)=REAL(STASH_CONTROL(st_output_bottom)+ML-1)
        ENDDO
      ENDIF
!
!
! Now reset the number of output levels to 1 if vertical compression is
! to be done in SPATIAL.  NB: index_lev and level_list need to be filled
! with values corresponding to the full range of levels processed.
!
      what_proc=STASH_CONTROL(st_proc_no_code)
      what_mean=(STASH_CONTROL(st_gridpoint_code)/block_size)*block_size
      IF(what_mean == vert_mean_base .OR. what_mean == global_mean_base &
     &   .OR. what_proc == st_time_series_code                          &
     &   .OR. what_proc == st_time_series_mean                          &
     &   .OR. what_proc == st_append_traj_code) num_levs_out=1
!
! Next compute the index for pseudo levels, if there are any
!
      IF(STASH_CONTROL(st_pseudo_in) >  0) THEN ! Input PSEUDO_LEVELS
        NI=STASH_CONTROL(st_pseudo_in)
        num_pseudo_in=STASH_PSEUDO_LEVELS(1,NI)
        IF(STASH_CONTROL(st_pseudo_out) >  0) THEN ! Output PSEUDO_LEVS
          NO=STASH_CONTROL(st_pseudo_out)
          num_pseudo_out=STASH_PSEUDO_LEVELS(1,NO)
          INDX1=0
          DO ML=1,NUM_PSEUDO_OUT
            ilev=STASH_PSEUDO_LEVELS(ML+1,NO)   !  Level required
            DO KL=1,NUM_PSEUDO_IN
              IF(STASH_PSEUDO_LEVELS(KL+1,NI) == ilev) THEN
                INDX1=INDX1+1
                INDEX_PSEUDO_LEV(INDX1)=KL
                pseudo_level(indx1)=ilev
                GOTO 500
              ENDIF
            ENDDO
            ICODE=nonsense
            WRITE(CMESSAGE,101) 'Output pseudo level ',ilev,            &
     &                          ' not found in input levels list'
            GOTO 999
 500        CONTINUE
          ENDDO
        ELSE  ! Illegal combination
          ICODE=nonsense
          WRITE(CMESSAGE,101) 'Input pseudo level list ',NI,            &
     &         ' has illegal output pseudo levels list'
          GOTO 999
        ENDIF
      ELSE  ! Only levels lists are supported for pseudo levels
        num_pseudo_out=0
      ENDIF
!
! Next expand the separate indexes and physical levels arrays into
! combined arrays if necessary, taking care not to overwrite earlier
! parts of the arrays.  If no pseudo-levels, set pseudo-level to 0.
!
      IF (num_pseudo_out >  0) THEN
        DO K2=num_pseudo_out,1,-1
          DO ML=1,num_levs_out
            INDEX_LEV(ML+(K2-1)*num_levs_out)=                          &
     &        (INDEX_PSEUDO_LEV(K2)-1)*num_levs_in+INDEX_LEV(ML)
            level(ML+(K2-1)*num_levs_out)=level(ML)
          ENDDO
          DO ML=num_levs_out,1,-1
            pseudo_level(ML+(K2-1)*num_levs_out)=pseudo_level(K2)
          ENDDO
        ENDDO
        num_levs_out=num_levs_out*num_pseudo_out
      ELSE
        DO ML=1,num_levs_out
          pseudo_level(ML)=0
        ENDDO
      ENDIF
!
999   CONTINUE ! jump here for error return
 101  FORMAT('STLEVELS : ',a,i6,a)
 103  FORMAT('STLEVELS : >> FATAL ERROR <<',a,2i5)
      RETURN
      END SUBROUTINE STLEVELS

