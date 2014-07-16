
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets the data_source array corresponding to the fields array.

Module Rcf_Set_Data_Source_Mod

!  Subroutine Rcf_Set_Data_Source - sets the data_source array.
!
! Description:
!   Allocates the data_source array and sets the elements according
!   to how the data in the the associated fields array should be
!   initialised.
!
! Method:
!   Uses data from the ITEMS namelist.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Set_data_source( data_source, fields_in, fields_out, &
                                field_count_in, field_count_out,    &
                                hdr_in, hdr_out )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Items_Mod             ! All of it

Use Rcf_Data_Source_Mod   ! Most of It

Use Rcf_Field_Type_Mod, Only : &
    Field_type

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Oper

Use Rcf_Submodel_Mod, Only : &
    Atmos_IM,            &
    Submodel_Ident

Use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Rcf_UMhead_Mod, Only : &
    UM_Header_type

Use Rcf_Recon_Mod, Only : &
    nice

Use Rcf_Stashcodes_Mod   ! So many of these are used the whole lot
                         ! have been included. Anything starting
                         ! with stashcode_ came from here.

Use Rcf_Recon_Mod, Only : &
    GRIB

Implicit None

! Arguments
Type( data_source_type ), Pointer  :: data_source(:)
Type( field_type), Pointer         :: fields_in(:)
Type( field_type ) , Pointer       :: fields_out(:)
Type( um_header_type ), Intent(In) :: hdr_in
Type( um_header_type ), Intent(In) :: hdr_out
Integer                            :: field_count_in
Integer                            :: field_count_out

! Comdecks
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
!*L------------------ COMDECK LOOKADD ----------------------------------
!LL
!LL Purpose : Contains information about the format
!LL           of the PP header
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   4.0  12/09/95   Change NPERIODS to LBUSER3, BRSVD1 to BULEV,
!LL                   BRSVD2 to BHULEV and definitions for BRLEV and
!LL                   BHRLEV. Corresponding changes made to STWORK1A
!LL                   and PPHEAD1A. (Andrew Brady)
!LL  4.0  12/10/95  Change item 45 from lbuser7 to model_code. RTHBarnes
!LL  5.1  17/04/00    Fixed/Free format. P.Selwood.
!LL  5.2  25/09/00    Add LBCC_xxxx variables for the compressed
!LL                   LBC LOOKUP array                  P.Burton
!LL
!LL Programming standard :
!LL
!LL Logical components covered : F092
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LLEND -----------------------------------------------------------------
!

! Validity time
      Integer, Parameter :: LBYR   =1   ! Year
      Integer, Parameter :: LBMON  =2   ! Month
      Integer, Parameter :: LBDAT  =3   ! Day of month
      Integer, Parameter :: LBHR   =4   ! Hour
      Integer, Parameter :: LBMIN  =5   ! Minute
      Integer, Parameter :: LBDAY  =6   ! Day number

! Data time
      Integer, Parameter :: LBYRD  =7   ! Year
      Integer, Parameter :: LBMOND =8   ! Month
      Integer, Parameter :: LBDATD =9   ! Day of month
      Integer, Parameter :: LBHRD  =10  ! Hour
      Integer, Parameter :: LBMIND =11  ! Minute
      Integer, Parameter :: LBDAYD =12  ! Day number

      Integer, Parameter :: LBTIM  =13  ! Time indicator
      Integer, Parameter :: LBFT   =14  ! Forcast period (hours)
      Integer, Parameter :: LBLREC =15  ! Length of data record
      Integer, Parameter :: LBCODE =16  ! Grid type code
      Integer, Parameter :: LBHEM  =17  ! Hemisphere indicator
      Integer, Parameter :: LBROW  =18  ! Number of rows in grid
      Integer, Parameter :: LBNPT  =19  ! Number of points per row
      Integer, Parameter :: LBEXT  =20  ! Length of extra data
      Integer, Parameter :: LBPACK =21  ! Packing method indicator
      Integer, Parameter :: LBREL  =22  ! Header release number
      Integer, Parameter :: LBFC   =23  ! Field code
      Integer, Parameter :: LBCFC  =24  ! Second field code
      Integer, Parameter :: LBPROC =25  ! Processing code
      Integer, Parameter :: LBVC   =26  ! Vertical coordinate type
      Integer, Parameter :: LBRVC  =27  ! Coordinate type for reference
                                        ! level

      Integer, Parameter :: LBEXP  =28  ! Experiment number
      Integer, Parameter :: LBEGIN =29  ! Start record
      Integer, Parameter :: LBNREC =30  ! No of records-Direct access
                                        ! only
      Integer, Parameter :: LBPROJ =31  ! Met-O-8 projection number
      Integer, Parameter :: LBTYP  =32  ! Met-O-8 field type
      Integer, Parameter :: LBLEV  =33  ! Met-O-8 level code
      Integer, Parameter :: LBRSVD1=34  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD2=35  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD3=36  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD4=37  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBSRCE =38  ! =1111 to indicate following
                                        ! apply to UM
      Integer, Parameter :: DATA_TYPE =39  ! Indicator for real/int or
                                           ! timeseries
      Integer, Parameter :: NADDR  =40  ! Start address in DATA_REAL or
                                        ! DATA_INT
      Integer, Parameter :: LBUSER3=41  ! Free for user-defined function
      Integer, Parameter :: ITEM_CODE =42  !Stash code
      Integer, Parameter :: LBPLEV =43  ! Pseudo-level indicator (if
                                        ! defined)
      Integer, Parameter :: LBUSER6=44  ! Free for user-defined function
      Integer, Parameter :: MODEL_CODE =45 ! internal model identifier

      Integer, Parameter :: BULEV  =46  ! Upper level boundary
      Integer, Parameter :: BHULEV =47  ! Upper level boundary
      Integer, Parameter :: BRSVD3 =48  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BRSVD4 =49  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BDATUM =50  ! Datum value
      Integer, Parameter :: BACC   =51  ! (Packed fields) Packing
                                        ! accuracy
      Integer, Parameter :: BLEV   =52  ! Level
      Integer, Parameter :: BRLEV  =53  ! Lower level boundary
      Integer, Parameter :: BHLEV  =54  ! (Hybrid levels) A-level of
                                        ! value
      Integer, Parameter :: BHRLEV =55  ! Lower level boundary
      Integer, Parameter :: BPLAT  =56  ! Real latitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BPLON  =57  ! Real longitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BGOR   =58  ! Grid orientation
      Integer, Parameter :: BZY    =59  ! Zeroth latitude
      Integer, Parameter :: BDY    =60  ! Latitude interval
      Integer, Parameter :: BZX    =61  ! Zeroth longitude
      Integer, Parameter :: BDX    =62  ! Longitude interval
      Integer, Parameter :: BMDI   =63  ! Missing data indicator
      Integer, Parameter :: BMKS   =64  ! M,K,S scaling factor

      Integer, Parameter :: LBCC_LBYR   = 1  ! Year
      Integer, Parameter :: LBCC_LBMON  = 2  ! Month
      Integer, Parameter :: LBCC_LBDAT  = 3  ! Day of the month
      Integer, Parameter :: LBCC_LBHR   = 4  ! Hour
      Integer, Parameter :: LBCC_LBMIN  = 5  ! Minute
      Integer, Parameter :: LBCC_LBDAY  = 6  ! Day number
      Integer, Parameter :: LBCC_LBEGIN = 7  ! Start record
      Integer, Parameter :: LBCC_NADDR  = 8  ! Start address of DATA
! Mapping of MPP_LOOKUP; analogous to mapping in PP header

      Integer, Parameter :: P_NADDR=1    ! Address on local PE
      Integer, Parameter :: P_LBLREC=2   ! Local length of record

!*----------------------------------------------------------------------
! NADDR IS LOCATION IN PP-HEADER (LOOKUP) FOR START POSN OF VARIABLE
! ITEM_CODE is the location in PP header for a code defined as
!           (section number)*1000+item number
! DATA_TYPE is the location in the PP header defining data as REAL or
!           INTEGER.
! LBNPT is the location defining the number of points per row
!

! Local variables
Integer                            :: item_index
Integer                            :: i
Integer                            :: j
Integer                            :: pos
Integer                            :: ErrorStatus
Integer                            :: pos_itc
Character (Len=*), Parameter       :: RoutineName='Rcf_Set_Data_Source'
Character (Len=80)                 :: Cmessage
Character (Len=*), Parameter       :: form="(2i5,4i4,' ',e9.3,' ',50a)"


!---------------------------------------------------------------
! check some of the arrays needed have been set up
!---------------------------------------------------------------

If ( Associated( data_source ) ) Then
  ErrorStatus = -10
  Cmessage = 'data_source is already set up!'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

If (.NOT. Associated( fields_out )  )Then
  ErrorSTatus = 20
  Cmessage = 'Fields have not yet been set up'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!---------------------------------------------------------------
! Allocate the space
!---------------------------------------------------------------
Allocate( data_source( field_count_out ) )

!---------------------------------------------------------------
! Find any approriate item from namelist and set values
!----------------------------------------------------------------
Do i = 1, field_count_out
  item_index = 0
  Do j = 1, num_items
    If ( fields_out( i ) % stashmaster % item == item_array( j ) .AND.     &
         fields_out( i ) % stashmaster % section == sctn_array( j ) ) Then
      item_index = j
      Exit
    End If
  End Do

  If ( item_index > 0 ) Then       ! Have found the item - fill in gaps

    data_source( i ) % source      = source_array( item_index )
    data_source( i ) % Domain      = area_array  ( item_index )
    data_source( i ) % Ancil_SctnC = upas_array  ( item_index )
    data_source( i ) % Ancil_ItemC = upaa_array  ( item_index )
    data_source( i ) % RConst      = uprc_array  ( item_index )
    data_source( i ) % Ancil_File  = upaf_array  ( item_index )

  Else                             ! No item found - use defaults

    data_source( i ) % source      = Input_Dump
    data_source( i ) % Domain      = Whole_Grid
    data_source( i ) % Ancil_SctnC = 0
    data_source( i ) % Ancil_ItemC = 0
    data_source( i ) % RConst      = 0.0
    data_source( i ) % Ancil_File  = ' '

  End If

 ! Setting Source for items required from input dump but not
  ! available (and known about with relevant code written)
  Call Rcf_Locate( fields_out( i ) % stashmaster % section,         &
                   fields_out( i ) % stashmaster % item,            &
                   fields_in, field_count_in, pos, .TRUE.)

  If ( pos == 0 .AND. data_source( i ) % source == Input_Dump ) Then

  ! ----------------------------------
  ! This item is not in the input dump
  ! ----------------------------------

    Select Case( fields_out( i ) % stashmaster % model )

    Case ( Atmos_IM )                 ! Atmosphere Sub-Model

     ! ----------------
     ! Atmosphere items
     ! ----------------

      Select Case( fields_out( i ) % stashmaster % section )

      Case( stashcode_prog_sec )

        Select Case( fields_out( i ) % stashmaster % item )

        Case( stashcode_3d_cca,          & ! 3D Convective cloud
              stashcode_cca )              ! 2D Convective cloud

          data_source( i ) % source = Field_Calcs

        Case( stashcode_npp_pft_acc,     & ! Carbon accumulation fields
              stashcode_g_lf_pft_acc,    &
              stashcode_g_ph_lf_pft_acc, &
              stashcode_rsp_w_pft_acc,   &
              stashcode_rsp_s_acc,       &
              stashcode_catch_snow,      & ! NLT canopy snow capacity
              stashcode_catch_tile,      & ! Tiled canopy capacity
              stashcode_z0_tile          & ! Tiled roughness length
                                   )
          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = -1.0

        Case( stashcode_can_water_tile,  & ! Vegetation fields
              stashcode_tstar_tile,      &
              stashcode_ice_temp_cat,    &
              stashcode_ice_conc_cat,    & ! ice conc (fraction) cats
              stashcode_ice_thick_cat,   & ! ice thickness categories
              stashcode_Gamtot,          & !\ large-scale hydrology
              stashcode_Zw,              & !/ fields
              stashcode_Sthzw,           & !/
              stashcode_Fsat,            & !\.
              stashcode_Fwetl,           & !/.
              stashcode_a_fsat,         & !\.
              stashcode_c_fsat,         & !/.
              stashcode_a_fwet,         & !\.
              stashcode_c_fwet,         & !/.
              stashcode_snow_tile)

          data_source( i ) % source = Field_Calcs

        Case( stashcode_infil_max_tile,  &
              stashcode_snow_grnd,       & ! Snow beneath NLT canopy
              stashcode_snow_on_ice,     &
              stashcode_ice_snow_depth,  &
              stashcode_sw_down_tile,    &
              stashcode_sw_down,         &
              stashcode_lw_up_diff,      &
              stashcode_mmr_smoke_fr,    &
              stashcode_mmr_smoke_ag,    &
              stashcode_mmr_smoke_cl,    &
              stashcode_biom_surf_em,    &
              stashcode_biom_elev_em,    &
              stashcode_dust1_mmr,       &
              stashcode_dust2_mmr,       &
              stashcode_dust3_mmr,       &
              stashcode_dust4_mmr,       &
              stashcode_dust5_mmr,       &
              stashcode_dust6_mmr,       &
              stashcode_albedo_sice,     &
              stashcode_albedo_land,     &
              stashcode_soilcarb_dpm,    & ! RothC soil C prognostic
              stashcode_soilcarb_rpm,    & ! RothC soil C prognostic
              stashcode_soilcarb_bio,    & ! RothC soil C prognostic
              stashcode_soilcarb_hum,    & ! RothC soil C prognostic
              stashcode_qcf2,            & ! second cloud ice prognostic
              stashcode_qrain,           & ! prognostic rain
              stashcode_qgraup,          & ! prognostic graupel
              stashcode_lcbase,          & ! lowest convective cloud base
              stashcode_3d_ccw           & ! CCW profile sent to radiation
                                        )

          data_source( i ) % source = set_to_zero

        Case( stashcode_rgrain )

          data_source( i ) % source = set_to_const
          data_source( i ) % rconst = 50.0

        Case( stashcode_tstar_land,     &
              stashcode_tstar_sea,      &
              stashcode_tstar_sice )

          data_source( i ) % source = Already_Processed

        End Select                       ! based on STASH item code

        !kdcorbin, 05/10 - added boundary layer section
        Case (stashcode_bl_sec)
              data_source(i) % source = set_to_zero
 
     End Select                         ! based on STASH section code

      !For files not found when reading from GRIB data
      If ( GRIB ) Then

        Select Case( fields_out( i ) % stashmaster % section )

        Case( stashcode_prog_sec )

          Select Case( fields_out( i ) % stashmaster % item )

          Case( stashcode_sea_ice_temp,     &
                stashcode_u_adv,            &
                stashcode_v_adv)

            data_source( i ) % source = Field_Calcs

          Case( stashcode_w_adv,            &
                stashcode_qcf,              &
                stashcode_cc_lwp,           &
                stashcode_unfrozen_soil,    &
                stashcode_frozen_soil,      &
                stashcode_qcl,              &
                stashcode_n_turb_mixlvs,    &
                stashcode_lvl_bse_dp_sc,    &
                stashcode_lvl_top_dp_sc,    &
                stashcode_bl_conv_flag,     &
                stashcode_turb_temp,        &
                stashcode_turb_humid,       &
                stashcode_area_cf,          &
                stashcode_bulk_cf,          &
                stashcode_liquid_cf,        &
                stashcode_frozen_cf,        &
                stashcode_sfc_zonal_cur,    &
                stashcode_sfc_merid_cur,    &
                stashcode_3d_cca,           & ! has to overide source=8
                stashcode_can_water_tile,   &
                stashcode_rho,              & ! rho calc'd after interp
                stashcode_exner,            &
                stashcode_can_conduct,      &
                stashcode_ccb,              &
                stashcode_cct,              &
                stashcode_mean_canopyw,     &
                stashcode_surf_z_curr,      &
                stashcode_surf_m_curr,      &
                stashcode_w)

            data_source( i ) % source = Set_To_Zero

          Case( stashcode_bl_depth )

            data_source( i ) % source = set_to_const
            data_source( i ) % rconst = 500.000

          Case( stashcode_z0 )

            data_source( i ) % source = set_to_const
            data_source( i ) % rconst = 0.500

          End Select                 ! select by item code
        End Select                   ! select by section code

      End If                       ! If GRIB section

    Case Default                 ! Couldn't find a proper Internal Model

      ! The Model type didn't match one of the three defined types.
      Write (Cmessage,*) " Couldn't find a Sub-Model ID type for : ", &
                  " Model ", fields_out( i ) % stashmaster % model,   &
                  " Section ", fields_out( i ) % stashmaster % model, &
                  " Item ",  fields_out( i ) % stashmaster % item
      ErrorStatus = 25
      Call Ereport( RoutineName, ErrorStatus, Cmessage )

    End Select                   ! select by Internal model

    ! Check that Source is now set correctly otherwise, fail
    If ( data_source( i ) % source == Input_Dump ) Then
      Write ( Cmessage, *) 'Section ',                              &
                           fields_out( i ) % stashmaster % section, &
                           'Item ',                                 &
                           fields_out( i ) % stashmaster % item ,   &
                           ' : Required field is not in input dump!'
      ErrorStatus = 30
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

  End If  ! If (pos == 0 ..)

  ! ---------------------------------------------------------------
  ! Some fields, even when found in the input dump, still need
  ! extra work. Mostly for fields that may have differing numbers
  ! of pseudo levels. But also fields needing to be recalculated to
  ! ensure consistancy with others.
  ! ---------------------------------------------------------------

  If ( pos /= 0 .AND. data_source( i ) % source == Input_Dump ) Then

  ! ----------------------------------
  ! This item is in the input dump
  ! ----------------------------------

    Select Case( fields_out( i ) % stashmaster % model )

    ! -----------
    ! Atmos items
    ! -----------

    Case ( Atmos_IM )
      
      Select Case( fields_out( i ) % stashmaster % section )

      Case( stashcode_prog_sec )

          Select Case( fields_out( i ) % stashmaster % item )

          Case( stashcode_icethick)

            If (h_int_active) Then
              data_source( i ) % source = Field_Calcs
            End If

          Case( stashcode_ice_temp_cat,       &
                stashcode_ice_conc_cat,       &
                stashcode_ice_thick_cat )

          !-------------------------------------------------------------
          ! Check number of pseudo levels in input dump is same as
          ! number of sea ice categories (nice) in output dump.
          ! If not -recalculate
          !-------------------------------------------------------------

            Call Rcf_Locate( fields_out( i ) % stashmaster % section, &
                             fields_out( i ) % stashmaster % item,    &
                             fields_in, field_count_in, pos_itc )

            If (fields_in(pos_itc) % levels /= nice) Then
               data_source( i ) % source = Field_Calcs
               If ( mype == 0 .AND.                                   &
                   Printstatus >= PrStatus_Normal ) Then
                Write (6,*) "Setting source of stashcode item ",      &
                  fields_out( i ) % stashmaster % item,               &
                  " to calculations due to difference in NICE"
              EndIf
            EndIf

          Case( stashcode_ice_snow_depth )

          !-------------------------------------------------------------
          ! Check number of pseudo levels in input dump is same as
          ! number of sea ice categories (nice) in output dump.
          ! If not -set to zero
          !-------------------------------------------------------------

            Call Rcf_Locate( fields_out( i ) % stashmaster % section, &
                             fields_out( i ) % stashmaster % item,    &
                             fields_in, field_count_in, pos_itc )

            If (fields_in(pos_itc) % levels /= nice) Then
               data_source( i ) % source = set_to_zero
               If ( mype == 0 .AND.                                   &
                   Printstatus >= PrStatus_Normal ) Then
                Write (6,*) "Setting source of stashcode item ",      &
                  fields_out( i ) % stashmaster % item,               &
                  " to -set-to-zero- due to difference in NICE"
              EndIf
            EndIf

          Case( stashcode_infil_max_tile, &
                stashcode_sw_down_tile )

          !-------------------------------------------------------------
          ! Check number of tiles in input dump is same as
          ! number of tiles in output dump.
          ! If not -set to zero
          !-------------------------------------------------------------

            Call Rcf_Locate( fields_out( i ) % stashmaster % section, &
                             fields_out( i ) % stashmaster % item,    &
                             fields_in, field_count_in, pos_itc )

            If (fields_in(pos_itc) % levels /=                        &
                fields_out( i )    % levels ) Then
               data_source( i ) % source = set_to_zero
               If ( mype == 0 .AND.                                   &
                   Printstatus >= PrStatus_Normal ) Then
                Write (6,*) "Setting source of stashcode item ",      &
                  fields_out( i ) % stashmaster % item,               &
                  " to -set-to-zero- due to difference in no. of tiles"
              EndIf
            EndIf

          Case( stashcode_catch_tile,     &
                stashcode_z0_tile )

          !-------------------------------------------------------------
          ! Check number of tiles in input dump is same as
          ! number of tiles in output dump.
          ! If not -set to user constant -1
          !-------------------------------------------------------------

            Call Rcf_Locate( fields_out( i ) % stashmaster % section, &
                             fields_out( i ) % stashmaster % item,    &
                             fields_in, field_count_in, pos_itc )

            If (fields_in(pos_itc) % levels /=                        &
                fields_out( i )    % levels ) Then
               data_source( i ) % source = set_to_const
               data_source( i ) % rconst = -1.000
               If ( mype == 0 .AND.                                   &
                   Printstatus >= PrStatus_Normal ) Then
                Write (6,*) "Setting source of stashcode item ",      &
                  fields_out( i ) % stashmaster % item,               &
                  " to -1 due to difference in no. of tiles"
              EndIf
            EndIf

          Case( stashcode_can_water_tile, &
                stashcode_tstar_tile,     &
                stashcode_snow_tile)

          !-------------------------------------------------------------
          ! Check number of tiles in input dump is same as
          ! number of tiles in output dump.
          ! If not -set to Field Calcs
          !-------------------------------------------------------------

            Call Rcf_Locate( fields_out( i ) % stashmaster % section, &
                             fields_out( i ) % stashmaster % item,    &
                             fields_in, field_count_in, pos_itc )

            If (fields_in(pos_itc) % levels /=                        &
                fields_out( i )    % levels ) Then
               data_source( i ) % source = Field_Calcs
               If ( mype == 0 .AND.                                   &
                   Printstatus >= PrStatus_Normal ) Then
                Write (6,*) "Setting source of stashcode item ",      &
                  fields_out( i ) % stashmaster % item,               &
                  " to -FielcCalcs- due to difference in no. of tiles"
              EndIf
            EndIf
        End Select                        ! Select on item number
      End Select                          ! Select on section number
    
    Case Default                 ! Couldn't find a proper Internal Model

      ! The Model type didn't match one of the three defined types.
      Write (Cmessage,*) " Couldn't find a Sub-Model ID type for : ",   &
                  " Model ", fields_out( i ) % stashmaster % model,     &
                  " Section ", fields_out( i ) % stashmaster % section, &
                  " Item ",  fields_out( i ) % stashmaster % item
      ErrorStatus = 70
      Call Ereport( RoutineName, ErrorStatus, Cmessage )

    End Select                            ! Select on Internal model

  End If                                  ! Item was in Dump

  ! Check on Rimwidtha if an atmos LBC file for copying
  If ( data_source( i ) % source == Input_Dump .AND.         &
     (fields_out( i ) % stashmaster % grid_type ==ppx_atm_lbc_theta.OR.&
      fields_out( i ) % stashmaster % grid_type == ppx_atm_lbc_u   .OR.&
      fields_out( i ) % stashmaster % grid_type == ppx_atm_lbc_v) ) Then

    If ( hdr_in  % Lookup( lbrow, fields_out( i ) % dump_pos ) /=  &
         hdr_out % Lookup( lbrow, fields_out( i ) % dump_pos) ) Then
      Cmessage = 'Rimwidth needs to be equal for input and output grids&
                 &if LBCs are copied'
      ErrorStatus = 40
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If
  End If

End Do

!------------------------------------------------------------------
! Print out some diagnostics if required....
!------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Oper ) Then
  Write (6,*)
  Write (6,*) 'Data source'

  Do i = 1, field_count_out
    Write (6,form) fields_out( i ) % stashmaster % section,            &
                   fields_out( i ) % stashmaster % item,               &
                   data_source( i ) % source,                          &
                   data_source( i ) % Domain,                          &
                   data_source( i ) % Ancil_SctnC,                     &
                   data_source( i ) % Ancil_ItemC,                     &
                   data_source( i ) % RConst,                          &
   data_source(i) % Ancil_File(1:Len_Trim(data_source(i) % Ancil_File))
  End Do
  Write (6,*)

End If


Return
End Subroutine Rcf_Set_Data_Source
End Module Rcf_Set_Data_Source_Mod
