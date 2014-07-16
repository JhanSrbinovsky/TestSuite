
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Ancillary processing for the Atmosphere model

Module Rcf_Ancil_Atmos_Mod

!  Subroutine Rcf_Ancil_Atmos  - Ancillary processing for Atmosphere
!
! Description:
!    Controls all ancillary processing for the atmosphere model
!
! Method:
!    1. Reads in the ancilmaster records to get information on the
!       ancillary fields and files.
!    2. Determines workspace required for ancillary lookups and data.
!    3. Calls INANCILA to read in ancillary lookups.
!    4. Calls REPLANCA to read in ancillary data.
!
! Current Code Owner: D. Robinson
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Ancil_Atmos ( Hdr_In, Hdr_Out,              &
                             Fields_In, Field_Count_In,    &
                             Fields_Out, Field_Count_Out )

Use Ancil_mod, Only : &
    anc_record,          anc_file,            &
    ancrecs,             max_ancrecs,         &
    ancFiles,            max_ancfiles,        &
    n_uancil,            uancfils,            &
    nlookup,             lookup_step,         &
    stashancil,          levels,              &
    ancil_add

Use Rcf_Ppx_Info_Mod, Only : &
    ppxRecs

Use calc_nlookups_Mod, Only : &
    calc_nlookups

Use Rcf_calc_len_ancil_Mod, Only : &
    Rcf_calc_len_ancil

Use getanc_fields_Mod, Only :  &
    getanc_fields

Use getanc_files_Mod, Only :  &
    getanc_files

Use hdancilm_Mod, Only : &
    hdancilm

Use rcf_interpolate_Mod, Only : &
    rcf_interpolate

Use Rcf_Items_Mod, Only :       &
    Num_Items,              &
    Sctn_Array,             &
    Item_Array,             &
    Source_Array

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

Use Rcf_Submodel_Mod, Only :    &
    N_Internal_Model,       &
    Internal_Model_List,    &
    Atmos_IM

Use Rcf_Grid_Type_Mod, Only :   &
    Input_Grid,             &
    Output_Grid

Use Rcf_HeadAddress_Mod, Only : &
    fh_vtyear,              &
    fh_vtmonth,             &
    fh_vtday,               &
    fh_vthour,              &
    fh_vtminute,            &
    fh_vtsecond,            &
    rc_latspacing,          &
    rc_firstlat

Use Rcf_UMhead_Mod, Only :  &
    Um_header_type,         &
    LenFixHd

Use Rcf_Lsm_Mod, Only :     &
    glob_lsm_out,           &
    local_lsm_out,          &
    glob_land_out,          &
    local_land_out

Use Rcf_DecompTP_Mod, Only :    &
    Decomp_rcf_input,       &
    Decomp_rcf_output

Use Rcf_Alloc_Field_mod, Only :  &
    Rcf_Alloc_Field,             &
    Rcf_Dealloc_Field

Use rcf_read_field_mod, Only : &
    Rcf_Read_Field

Use rcf_write_field_mod, Only : &
    Rcf_Write_Field

Use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Locate_Anc_mod, Only : &
    Locate_Anc_Field

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Parvars_Mod, Only : &
    mype

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   &
    interp_copy

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_tstar,           &
    stashcode_land_frac,       &
    stashcode_tstar_land,      &
    stashcode_tstar_sea,       &
    stashcode_tstar_sice,      &
    stashcode_icefrac,         &
    stashcode_tstar_anom,      &
    stashcode_prog_sec

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Rcf_CntlAtm_Mod, Only : &
    L_CTILE

Use Rcf_Recon_Mod, Only : &
    tpps_ozone_levels

Implicit None

! Arguments

Type (um_header_type),      Intent(In) :: hdr_in
Type (um_header_type),      Intent(In) :: hdr_out
Type (field_type), Pointer             :: fields_in (:)
Type (field_type), Pointer             :: fields_out (:)
Integer ,                   Intent(In) :: field_count_in
Integer ,                   Intent(In) :: Field_count_out

! ----------------------------------------------------------------
! Arrays to store headers from ancillary files
! (Could go into a module & add USE to INANCILA & REPLANCA ?)

! Length of fixed header is defined in umhead

Integer, Parameter :: LenInthd_anc  = 15
Integer, Parameter :: LenRealhd_anc = 6

! The following arrays are 2-dimensional, the second dimension
! being the number of atmosphere ancillary files.

Integer, dimension(:,:), allocatable :: fixhd_ancil
Integer, dimension(:,:), allocatable :: inthd_ancil
Real   , dimension(:,:), allocatable :: realhd_ancil
Integer, dimension(:,:), allocatable :: lookup_ancil

! ----------------------------------------------------------------
! Variables & Work arrays for ancillary processing

Integer                            :: nlookups
Integer                            :: len_ancil
Integer                            :: ipt
Real,    dimension(:), allocatable :: ancil_data

Integer, dimension(:), allocatable :: lookup_start

! ----------------------------------------------------------------
! Local Variables

Integer :: pos_IceFrac_In  ! Position of Ice Fraction   in Input  dump
Integer :: pos_IceFrac_Out ! Position of Ice Fraction   in Output dump
Integer :: pos_Tstar_In    ! Position of T Star         in Input  dump
Integer :: pos_Tstar_Out   ! Position of T Star         in Output dump
Integer :: pos_Land_Frac_In     ! Position of land frac in Input  dump
Integer :: pos_Land_Frac_Out    ! Position of land frac in Output dump
Integer :: pos_Tstar_Land_In    ! Position of T Star_Land in Input  dump
Integer :: pos_Tstar_Land_Out   ! Position of T Star_Land in Output
Integer :: pos_Tstar_Sea_In     ! Position of T Star_Sea in Input
Integer :: pos_Tstar_Sea_Out    ! Position of T Star_Sea in Output
Integer :: pos_Tstar_Sice_In    ! Position of T Star_Sice in Input
Integer :: pos_Tstar_Sice_Out   ! Position of T Star_Sice in Output
Integer :: ipos_111             ! Land fraction ancil position
Integer :: pos_Tstar_Anom  ! Position of T Star Anomoly in Input  dump
Integer :: pos

Type (field_type), Pointer   :: IceFrac_In
Type (field_type), Pointer   :: TStar_In
Type (field_type), Pointer   :: Land_Frac_In
Type (field_type), Pointer   :: TStar_Land_In
Type (field_type), Pointer   :: TStar_Sea_In
Type (field_type), Pointer   :: TStar_Sice_In
Type (field_type), Pointer   :: IceFrac_Out
Type (field_type), Pointer   :: TStar_Out
Type (field_type), Pointer   :: Land_Frac_Out
Type (field_type), Pointer   :: TStar_Land_Out
Type (field_type), Pointer   :: TStar_Sea_Out
Type (field_type), Pointer   :: TStar_Sice_Out
Type (field_type), Pointer   :: TStar_Anom

Type (field_type)            :: dummy
Type (field_type)            :: dummy_anc
Type (field_type), Target    :: dummy_IceFr
Type (field_type), Target    :: dummy_TStar
Type (field_type), Target    :: dummy_Land_Frac
Type (field_type), Target    :: dummy_TStar_Land
Type (field_type), Target    :: dummy_TStar_Sea
Type (field_type), Target    :: dummy_TStar_Sice
Type (field_type), Target    :: dummy_TStar_Anom

Integer      :: i,j,k          !  Loop indices
Integer      :: i_full
Integer      :: i_land
Integer      :: irec
Integer      :: n_anc_read
Integer      :: ianc_Mask
Integer      :: ianc_IceFrac
Integer      :: ianc_TStar
Integer      :: ianc_u_curr
Integer      :: ianc_v_curr
Integer      :: ErrorStatus

Logical            :: Found      ! T : ANCILmaster record found
Logical, Parameter :: UserAnc = .true.

Character (Len=13)           :: AM_filename ! ANCILmaster filename
Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Ancil_Atmos'

Integer             ::  ppxi  !  Dummy argument for INANCILA
Character (Len=1)   ::  ppxc  !  Dummy argument for INANCILA

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

! ----------------------------------------------------------------

! Make sure data portions are initially null
Nullify(dummy % Data)
Nullify(dummy % Data_Int)
Nullify(dummy % Data_Log)
Nullify(dummy_IceFr % Data)
Nullify(dummy_IceFr % Data_Int)
Nullify(dummy_IceFr % Data_Log)
Nullify(dummy_TStar % Data)
Nullify(dummy_TStar % Data_Int)
Nullify(dummy_TStar % Data_Log)
Nullify(dummy_Land_Frac % Data)
Nullify(dummy_Land_Frac % Data_Int)
Nullify(dummy_Land_Frac % Data_Log)
Nullify(dummy_TStar_Land % Data)
Nullify(dummy_TStar_Land % Data_Int)
Nullify(dummy_TStar_Land % Data_Log)
Nullify(dummy_TStar_Sea % Data)
Nullify(dummy_TStar_Sea % Data_Int)
Nullify(dummy_TStar_Sea % Data_Log)
Nullify(dummy_TStar_Sice % Data)
Nullify(dummy_TStar_Sice % Data_Int)
Nullify(dummy_TStar_Sice % Data_Log)

Nullify(dummy_TStar_Anom % Data)
Nullify(dummy_TStar_Anom % Data_Int)
Nullify(dummy_TStar_Anom % Data_Log)

Write (6,*)

! Read in ANCILmaster files and determine total number of records

Do i = 1, N_Internal_Model
  If ( Internal_Model_List(i) == Atmos_IM ) Then

    !  Ancillary Fields
    AM_filename = 'ANCILfields_A'
    Call hdancilm ( AM_filename, 'ANCILMSTR' )

    !  Ancillary Files
    AM_filename = 'ANCILfiles_A'
    Call hdancilm ( AM_filename, 'ANCILMSTR' )

  End If
End Do

! User ANCILmaster files - kdcorbin, 05/10
! Commented below:
!Do i = 1, n_uancil
!  AM_filename = uancfils(i)
!  Call hdancilm ( AM_filename, 'UANCLMSTR')
!End Do

!Added user information to ancillary fields and files lists:
AM_filename='UAFLDS_A'
Call hdancilm (AM_filename,'UANCLMSTR')

AM_filename='UAFILES_A'
Call hdancilm (AM_filename,'UANCLMSTR')

If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
  Write (6,*) 'Total No of ANCILmaster records (Fields) ',max_ancRecs
  Write (6,*) 'Total No of ANCILmaster records (Files ) ',max_ancFiles
End If

! ===============================================

Do i = 1, N_Internal_Model
  If ( Internal_Model_List(i) == Atmos_IM ) Then

    AM_filename = 'ANCILfields_A'
    Call getanc_fields ( AM_filename, 'ANCILMSTR' )

    AM_filename = 'ANCILfiles_A'
    Call getanc_files  ( AM_filename, 'ANCILMSTR' )

  End If
End Do

! User ANCILmaster files - kdcorbin, 05/10
! Commented below:
!Do i = 1, n_uancil
!  AM_filename = uancfils(i)
!  Call getanc_fields ( AM_filename, 'UANCLMSTR', .TRUE. )
!End Do

!Read in user information:
AM_filename = 'UAFLDS_A'
Call getanc_fields ( AM_filename,'UANCLMSTR')

AM_filename = 'UAFILES_A'
Call getanc_files( AM_filename, 'UANCLMSTR')

! ===============================================

! Not complete ; see old CONTROL
! Add check for L_SSTANOM & ITEMS=39/SOURCE=2

! Check that ANCILmaster & STASHmaster record exists for
! ITEMS namelists with SOURCE=2

do I=1,num_items
if (source_array(i) == 2) then

  found = .false.
  do irec=1,ancrecs
    if (item_array(i)   == anc_record(irec) % item_number .and. &
        sctn_array(i)   == anc_record(irec) % section_number) then
      found = .true.
      exit
    endif
  enddo

  if (.not.found) Then  !  No ANCILmaster record for this ITEM namelist
    ErrorStatus = 10
    write (CMessage,*) 'ITEMS namelist for Stash-Code', &
    sctn_array(i), item_array(i),' : No ANCILmaster record ?'
    call Ereport ( RoutineName, ErrorStatus, cmessage)
  endif

! Now check that a STASHmaster record also exists
  dummy_anc % stashmaster => rcf_exppx (atmos_im, sctn_Array(i), & 
                                                  item_Array(i))

endif   !  if source = 2
enddo

! Set up list of ancillary fields to be read in

anc_record (:) % anc_field_read = 0

do I=1,num_items
  do irec=1,ancrecs
    if (item_array(i)   == anc_record(irec) % item_number    .and. &
        sctn_array(i)   == anc_record(irec) % section_number .and. &
        source_array(i) == 2 ) then
      anc_record (irec) % anc_field_read = 1
    endif
  enddo
enddo

! Check if both surface currents are selected

call locate_anc_field (30, ianc_u_curr)
call locate_anc_field (31, ianc_v_curr)
if (anc_record (ianc_u_curr) % anc_field_read == 1 .or.    &
    anc_record (ianc_v_curr) % anc_field_read == 1) then
  anc_record (ianc_u_curr) % anc_field_read = 1
  anc_record (ianc_v_curr) % anc_field_read = 1
endif

! Sea surface temperature must be updated when sea ice is updated

call locate_anc_field (27, ianc_IceFrac)
call locate_anc_field (28, ianc_TStar)
if (anc_record (ianc_IceFrac) % anc_field_read == 1 .and.    &
    anc_record (ianc_TStar)   % anc_field_read <= 0) then
  anc_record (ianc_TStar) % anc_field_read = 1
endif

n_anc_read = 0
do irec=1,ancRecs
  if (anc_record (irec) % anc_field_read == 1) then
    n_anc_read = n_anc_read + 1
  endif
enddo

If ( PrintStatus >= PrStatus_Diag .and. mype == 0 ) Then
  If (n_anc_read > 0 ) then
    write (6,*)
    write (6,*) n_anc_read,' Ancillary fields to be read in :'
    Do i=1,ancRecs
      If (anc_record (i) % anc_field_read == 1) then
        write (6,*) anc_record (i) % ancil_ref_number,' ', &
                    anc_record (i) % anc_name
      End If
    End Do
  End If
Endif

! ===============================================

! Determine the no of lookup entries to be read in from ancillary files
Call calc_nlookups (nlookups)

If (PrintStatus >= PrStatus_Normal .and. mype == 0) Then
  write (6,*) ' rcf_ancil_atmos : nlookups ',nlookups
Endif

! Determine the workspace (len_ancil) required for the ancillaries.
Call rcf_calc_len_ancil (Output_Grid % loc_p_field,           &
                         Output_Grid % loc_r_field,           &
                         Output_Grid % loc_p_rows,  len_ancil )

If (PrintStatus >= PrStatus_Normal .AND. mype == 0 ) Then
  write (6,*) ' rcf_ancil_atmos : len_ancil ',len_ancil
Endif

! ===============================================

allocate ( fixhd_ancil(LenFixHd,ancFiles) )
allocate ( inthd_ancil(LenInthd_anc,ancFiles) )
allocate ( realhd_ancil(LenRealhd_anc,ancFiles) )

fixhd_ancil(:,:) = 0
inthd_ancil(:,:) = 0
realhd_ancil (:,:) = 0.0

! ===============================================

allocate ( lookup_ancil(Hdr_Out % Len1LookUp, nlookups) )
allocate ( nlookup(ancRecs) )
allocate ( lookup_step(ancRecs) )
allocate ( stashancil(ancRecs) )
allocate ( levels(ancRecs) )
allocate ( ancil_add(ancRecs) )
allocate ( lookup_start(ancFiles) )

ancil_add(:) = 0

! DEPENDS ON: inancila
CALL INANCILA(LenFixHd,                                    &
     LenInthd_anc,                                         &
     LenRealhd_anc,                                        &
     Hdr_Out % Len1LevDepC,                                &
     Hdr_Out % Len2LevDepC,                                &
     FIXHD_ANCIL,INTHD_ANCIL,REALHD_ANCIL,LOOKUP_ANCIL,    &
     Hdr_Out % RealC,                                      &
     Hdr_Out % LevDepC,                                    &
     NLOOKUPS,                                             &
     LOOKUP_START,                                         &
     Hdr_Out % Len1LookUp,                                 &
     Output_Grid % glob_p_row_length,                      &
     Output_Grid % loc_p_row_length,                       &
     Output_Grid % glob_p_rows,                            &
     Output_Grid % loc_p_rows,                             &
     Output_Grid % glob_u_rows,                            &
     Output_Grid % glob_r_row_length,                      &
     Output_Grid % glob_r_rows,                            &
     Output_Grid % loc_r_row_length,                       &
     Output_Grid % loc_r_rows,                             &
     Output_Grid % model_levels,                           &
     Output_Grid % tr_levels,                              &
     Output_Grid % st_levels,                              &
     Output_Grid % sm_levels,                              &
     Output_Grid % ozone_levels,                           &
     tpps_ozone_levels,                                    &
     ANCIL_ADD,                                            &
     PPXI,PPXC,ppxRecs,                                    &
     ErrorStatus,CMESSAGE)

if (ErrorStatus /= 0) then
  write (6,*) ' Error in INANCILA'
  write (6,*) ' CMESSAGE ',CMESSAGE
  write (6,*) ' ErrorStatus ',ErrorStatus
  Call Ereport ( RoutineName, ErrorStatus, cmessage)
endif

! ===============================================
! Locate Ice fraction in Output dump

call Rcf_Locate ( stashcode_prog_sec, stashcode_icefrac,             &
                  Fields_Out, Field_Count_Out, Pos_IceFrac_Out, .TRUE. )

If ( Pos_IceFrac_Out /= 0 ) Then  ! Ice Fraction in output dump

  IceFrac_Out => Fields_Out (Pos_IceFrac_Out)
  Call Rcf_Alloc_Field ( IceFrac_Out )

  ! Locate Ice fraction in input dump

  Call Rcf_Locate ( stashcode_prog_sec, stashcode_icefrac,           &
                    Fields_In, Field_Count_In, Pos_IceFrac_In, .TRUE. )

  If ( Pos_IceFrac_In /= 0) Then  ! Ice Fraction in input dump

    IceFrac_In => Fields_In (Pos_IceFrac_In)
    Call Rcf_Alloc_Field ( IceFrac_In )

! Read in Ice Fraction

    Call Rcf_Read_Field ( IceFrac_in, Hdr_In, Decomp_rcf_input )

    ! Set interpolation
    If (h_int_active) Then
      IceFrac_In % interp = interp_h_only
    Else
      IceFrac_In % interp = interp_copy
    End If

    Call rcf_interpolate( IceFrac_In, IceFrac_Out, Input_Grid,   &
                          Output_Grid, dummy, dummy )

    Call Rcf_Dealloc_Field ( IceFrac_In )

  Else ! Ice Fraction not in input dump

!  Ice_Fraction is not in the input dump. It is assumed here
!  that an ITEMS namelist with SOURCE=2 or similar exists for
!  the Ice Fraction Stash Code. The rcf will abort in Create_Dump
!  if the Output Ice Fraction field cannot be initialised.

  Endif

Else ! Ice Fraction not in output dump

! REPLANCA expects an Ice Fraction field of length p_field
! even if not required - set up and allocate.

  Dummy_IceFr % level_size = Output_Grid % loc_p_field
  Dummy_IceFr % levels     = 1
  Dummy_IceFr % stashmaster => Rcf_Exppx( 1, 0, stashcode_icefrac )

  IceFrac_Out => Dummy_IceFr
  Call Rcf_Alloc_Field (IceFrac_out)

Endif

! ===============================================

! Locate TStar in output dump.

call Rcf_Locate ( stashcode_prog_sec, stashcode_tstar,               &
                  Fields_Out, Field_Count_Out, Pos_TStar_Out, .TRUE. )

If ( Pos_TStar_Out /= 0 ) Then  ! TStar in output dump

  TStar_Out => Fields_Out (Pos_TStar_Out)
  Call Rcf_Alloc_Field ( TStar_Out )

  call Rcf_Locate ( stashcode_prog_sec, stashcode_tstar,             &
                    Fields_In, Field_Count_In, Pos_TStar_In, .TRUE. )

  If ( Pos_TStar_In /= 0) Then  ! TStar in input dump

    TStar_In => Fields_In (Pos_TStar_In)
    call Rcf_Alloc_Field ( TStar_In )

    call Rcf_Read_Field ( TStar_In, Hdr_In, Decomp_rcf_input )

   ! Set interpolation
    If (h_int_active) Then
      TStar_In % interp = interp_h_only
    Else
      TStar_In % interp = interp_copy
    End If

    Call rcf_interpolate( TStar_In, TStar_Out, Input_Grid,     &
                          Output_Grid, dummy, dummy )

    call Rcf_Dealloc_Field ( TStar_In )

  Else ! TStar not in input dump

!  TStar is not in the input dump. It is assumed here
!  that an ITEMS namelist with SOURCE=2 or similar exists for
!  the TStar Stash Code. The rcf will abort in Create_Dump
!  if the Output TStar field cannot be initialised.

  Endif

Else ! TStar not in output dump

! REPLANCA expects a TStar field of length p_field
! even if not required - set up and allocate.

  Dummy_TStar % level_size = Output_Grid % loc_p_field
  Dummy_TStar % levels     = 1
  Dummy_TStar % stashmaster => Rcf_Exppx( 1, 0, stashcode_tstar )

  TStar_Out => Dummy_TStar
  Call Rcf_Alloc_Field (TStar_out)

Endif

! ===============================================
! Locate Land_Frac in output dump.

call Rcf_Locate ( stashcode_prog_sec, stashcode_land_frac,           &
                  Fields_Out, Field_Count_Out, Pos_Land_Frac_Out, .TRUE. )

If ( Pos_Land_Frac_Out /= 0 ) Then  ! Land_Frac in output dump

  Land_Frac_Out => Fields_Out (Pos_Land_Frac_Out)
  Call Rcf_Alloc_Field ( Land_Frac_Out )

  Call Locate_Anc_Field (111, ipos_111)

! Check not reading land fraction from ancillary file:

  If(anc_record(ipos_111) % anc_field_read /= 1 )Then

    call Rcf_Locate ( stashcode_prog_sec, stashcode_land_frac,       &
                      Fields_In, Field_Count_In, Pos_Land_Frac_In, .TRUE. )

    If ( Pos_Land_Frac_In /= 0) Then  ! Land_Frac in input dump

      Land_Frac_In => Fields_In (Pos_Land_Frac_In)
      call Rcf_Alloc_Field ( Land_Frac_In )

      call Rcf_Read_Field ( Land_Frac_In, Hdr_In, Decomp_rcf_input )

! Dont allow horizonal interpolation for Land_Frac:
      If (h_int_active) Then
        ErrorStatus = 20
        write (CMessage,*) 'Horizontal interpolation of ' &
        ,'land fraction not allowed: Need ancillary'
        call Ereport ( RoutineName, ErrorStatus, cmessage)
      Else
        Land_Frac_In % interp = interp_copy
      End If

      Call rcf_interpolate( Land_Frac_In, Land_Frac_Out, Input_Grid, &
                          Output_Grid, dummy, dummy )


      call Rcf_Dealloc_Field ( Land_Frac_In )

    Else ! Land frac not in input dump

        ErrorStatus = 22
        write (CMessage,*) 'Land fraction is not in the input dump ' &
        ,'=> an ancillary must be provide.'
        call Ereport ( RoutineName, ErrorStatus, cmessage)

    Endif

  End If ! Land_Frac read from ancillary

Else ! Land_Frac not in output dump

    Write(6,*)'Land Frac is not in output dump => setting to dummy'

! REPLANCA expects a Land_Frac field of length p_field
! even if not required - set up and allocate.

  Dummy_Land_Frac % level_size = Output_Grid % loc_land_field
  Dummy_Land_Frac % levels     = 1
  Dummy_Land_Frac % stashmaster => Rcf_Exppx(1,0,stashcode_land_frac)

  Land_Frac_Out => Dummy_Land_Frac
  Call Rcf_Alloc_Field (Land_Frac_out)

Endif

! ===============================================
! Locate Tstar_Land in output dump.

call Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_land,          &
                  Fields_Out, Field_Count_Out, Pos_Tstar_Land_Out, .TRUE. )

If ( Pos_Tstar_Land_Out /= 0 ) Then  ! Tstar_Land in output dump

  Tstar_Land_Out => Fields_Out (Pos_Tstar_Land_Out)
  Call Rcf_Alloc_Field ( Tstar_Land_Out )

  call Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_land,        &
                    Fields_In, Field_Count_In, Pos_Tstar_Land_In, .TRUE. )

  If ( Pos_Tstar_Land_In /= 0) Then  ! Tstar_Land in input dump

    Tstar_Land_In => Fields_In (Pos_Tstar_Land_In)
    call Rcf_Alloc_Field ( Tstar_Land_In )

    call Rcf_Read_Field ( Tstar_Land_In, Hdr_In, Decomp_rcf_input )

   ! Set interpolation
    If (h_int_active) Then
      Tstar_Land_In % interp = interp_h_only
    Else
      Tstar_Land_In % interp = interp_copy
    End If

    Call rcf_interpolate( Tstar_Land_In, Tstar_Land_Out, Input_Grid, &
                          Output_Grid, dummy, dummy )

    call Rcf_Dealloc_Field ( Tstar_Land_In )

  Else ! Tstar_Land not in input dump

!  Tstar_Land is not in the input dump. It is assumed here
!  that an ITEMS namelist with SOURCE=2 or similar exists for
!  the Tstar_Land Stash Code. The rcf will abort in Create_Dump
!  if the Output Tstar_Land field cannot be initialised.

    Write(6,*)'Land T* is not in input dump => setting to gbm T*'
    Tstar_Land_out % Data (:,:) = Tstar_out % Data(:,:)

  Endif

Else ! Tstar_Land not in output dump

    Write(6,*)'Land T* is not in output dump => setting to dummy'

! REPLANCA expects a Tstar_Land field of length p_field
! even if not required - set up and allocate.

  Dummy_Tstar_Land % level_size = Output_Grid % loc_p_field
  Dummy_Tstar_Land % levels     = 1
  Dummy_Tstar_Land % stashmaster => Rcf_Exppx(1,0,stashcode_tstar_land)

  Tstar_Land_Out => Dummy_Tstar_Land
  Call Rcf_Alloc_Field (Tstar_Land_out)

Endif

! ===============================================
! Locate Tstar_Sea in output dump.

call Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_sea,           &
                  Fields_Out, Field_Count_Out, Pos_Tstar_Sea_Out, .TRUE. )

If ( Pos_Tstar_Sea_Out /= 0 ) Then  ! Tstar_Sea in output dump

  Tstar_Sea_Out => Fields_Out (Pos_Tstar_Sea_Out)
  Call Rcf_Alloc_Field ( Tstar_Sea_Out )

  call Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_sea,         &
                    Fields_In, Field_Count_In, Pos_Tstar_Sea_In, .TRUE. )

  If ( Pos_Tstar_Sea_In /= 0) Then  ! Tstar_Sea in input dump

    Tstar_Sea_In => Fields_In (Pos_Tstar_Sea_In)
    call Rcf_Alloc_Field ( Tstar_Sea_In )

    call Rcf_Read_Field ( Tstar_Sea_In, Hdr_In, Decomp_rcf_input )

   ! Set interpolation
    If (h_int_active) Then
      Tstar_Sea_In % interp = interp_h_only
    Else
      Tstar_Sea_In % interp = interp_copy
    End If

    Call rcf_interpolate( Tstar_Sea_In, Tstar_Sea_Out, Input_Grid, &
                          Output_Grid, dummy, dummy )

    call Rcf_Dealloc_Field ( Tstar_Sea_In )

  Else ! Tstar_Sea not in input dump

!  Tstar_Sea is not in the input dump. It is assumed here
!  that an ITEMS namelist with SOURCE=2 or similar exists for
!  the Tstar_Sea Stash Code. The rcf will abort in Create_Dump
!  if the Output Tstar_Sea field cannot be initialised.

    Write(6,*)'Open sea T* is not in input dump => setting to gbm T*'
    Tstar_Sea_out % Data (:,:) = Tstar_out % Data(:,:)

  Endif

Else ! Tstar_Sea not in output dump

! REPLANCA expects a Tstar_Sea field of length p_field
! even if not required - set up and allocate.

  Dummy_Tstar_Sea % level_size = Output_Grid % loc_p_field
  Dummy_Tstar_Sea % levels     = 1
  Dummy_Tstar_Sea % stashmaster => Rcf_Exppx(1,0,stashcode_tstar_sea)

  Tstar_Sea_Out => Dummy_Tstar_Sea
  Call Rcf_Alloc_Field (Tstar_Sea_out)

Endif

! ===============================================
! Locate Tstar_Sice in output dump.

call Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_sice,          &
                  Fields_Out, Field_Count_Out, Pos_Tstar_Sice_Out, .TRUE. )

If ( Pos_Tstar_Sice_Out /= 0 ) Then  ! Tstar_Sice in output dump

  Tstar_Sice_Out => Fields_Out (Pos_Tstar_Sice_Out)
  Call Rcf_Alloc_Field ( Tstar_Sice_Out )

  call Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_sice,        &
                    Fields_In, Field_Count_In, Pos_Tstar_Sice_In, .TRUE. )

  If ( Pos_Tstar_Sice_In /= 0) Then  ! Tstar_Sice in input dump

    Tstar_Sice_In => Fields_In (Pos_Tstar_Sice_In)
    call Rcf_Alloc_Field ( Tstar_Sice_In )

    call Rcf_Read_Field ( Tstar_Sice_In, Hdr_In, Decomp_rcf_input )

   ! Set interpolation
    If (h_int_active) Then
      Tstar_Sice_In % interp = interp_h_only
    Else
      Tstar_Sice_In % interp = interp_copy
    End If

    Call rcf_interpolate( Tstar_Sice_In, Tstar_Sice_Out, Input_Grid, &
                          Output_Grid, dummy, dummy )

    call Rcf_Dealloc_Field ( Tstar_Sice_In )

  Else ! Tstar_Sice not in input dump

!  Tstar_Sice is not in the input dump. It is assumed here
!  that an ITEMS namelist with SOURCE=2 or similar exists for
!  the Tstar_Sice Stash Code. The rcf will abort in Create_Dump
!  if the Output Tstar_Sice field cannot be initialised.

    Write(6,*)'Sea-ice T* is not in input dump => setting to gbm T*'
    Tstar_Sice_out % Data (:,:) = Tstar_out % Data(:,:)

  Endif

Else ! Tstar_Sice not in output dump

! REPLANCA expects a Tstar_Sice field of length p_field
! even if not required - set up and allocate.

  Dummy_Tstar_Sice % level_size = Output_Grid % loc_p_field
  Dummy_Tstar_Sice % levels     = 1
  Dummy_Tstar_Sice % stashmaster => Rcf_Exppx(1,0,stashcode_tstar_sice)

  Tstar_Sice_Out => Dummy_Tstar_Sice
  Call Rcf_Alloc_Field (Tstar_Sice_out)

Endif
! ===============================================
! T Anolomy Option to be sorted out yet.
! Locate TStar in output dump.

call Rcf_Locate ( stashcode_prog_sec, stashcode_tstar_anom,          &
                  Fields_Out, Field_Count_Out, Pos_TStar_Anom, .TRUE. )

If ( Pos_TStar_Anom /= 0 ) Then  ! TStar Anomoly in output dump

! Not catered for yet
  ErrorStatus = 111
  CMESSAGE = ' New Rcf doesnt handle the T* Anomoly field yet.'
  write (6,*) ' CMESSAGE ',CMESSAGE
  write (6,*) ' ErrorStatus ',ErrorStatus
  Call Ereport ( RoutineName, ErrorStatus, cmessage)

Else ! TStar Anomoly not in output dump

! REPLANCA expects a TStar Anomly field of length p_field
! even if not required - set up and allocate.

  Dummy_TStar_Anom % level_size = Output_Grid % loc_p_field
  Dummy_TStar_Anom % levels     = 1
  Dummy_TStar_Anom % stashmaster =>  &
                             Rcf_Exppx( 1, 0, stashcode_tstar_anom )

  TStar_Anom => Dummy_TStar_Anom
  Call Rcf_Alloc_Field (TStar_Anom)

Endif

! ===============================================
allocate ( ancil_data(len_ancil) )
ancil_data(:) = 0.0

! DEPENDS ON: replanca
CALL REPLANCA (                                           &
     Hdr_Out % FixHd (FH_VTYear),                         &
     Hdr_Out % FixHd (FH_VTMonth),                        &
     Hdr_Out % FixHd (FH_VTDay),                          &
     Hdr_Out % FixHd (FH_VTHour),                         &
     Hdr_Out % FixHd (FH_VTMinute),                       &
     Hdr_Out % FixHd (FH_VTSecond),                       &
     Output_Grid % loc_p_field,                           &
     Output_Grid % loc_p_rows,                            &
     Output_Grid % loc_u_field,                           &
     Output_Grid % loc_v_field,                           &
     Output_Grid % loc_r_field,                           &
     Output_Grid % loc_land_field,                        &
     ancil_data,                                          &
     local_lsm_out,                                       &
     IceFrac_Out % Data,                                  &
     TStar_Out % Data,                                    &
     Land_Frac_Out % Data,                                &
     TStar_Land_Out % Data,                               &
     TStar_Sea_Out % Data,                                &
     TStar_Sice_Out % Data,                               &
     TStar_Anom % Data,                                   &
     Hdr_Out % RealC (RC_LatSpacing),                     &
     Hdr_Out % RealC (RC_FirstLat),                       &
     Hdr_Out % Len1LookUp, LenFixHd,                      &
     LenInthd_anc,                                        &
     LEN_ANCIL,FIXHD_ANCIL,                               &
     INTHD_ANCIL,LOOKUP_ANCIL,LOOKUP_ANCIL,               &
     LOOKUP_START,NLOOKUPS,                               &
     ErrorStatus,CMESSAGE)

if (ErrorStatus /= 0) then
  write (6,*) ' Error in REPLANCA'
  write (6,*) ' CMESSAGE ',CMESSAGE
  write (6,*) ' ErrorStatus ',ErrorStatus
  Call Ereport ( RoutineName, ErrorStatus, cmessage)
endif

! ===============================================
! Write ancillary fields into output dump

call locate_anc_field (1, ianc_Mask)
call locate_anc_field (27,ianc_IceFrac)
call locate_anc_field (28,ianc_TStar)

do k=1,ancrecs

  if (anc_record(k) % anc_field_read == 1) then

    if (k == ianc_Mask .or. k == ianc_IceFrac .or. k == ianc_TStar) then

!     The ancillary fields : Land-Sea mask, Ice Fraction and TStar
!     are not written to the output dump in this loop.

      Cycle

    else ! write this ancillary to the output dump
!

! Locate position of ancillary field in output dump
      if (ancil_add(k).ne.0) then
        call Rcf_Locate ( anc_record(k) % section_number,             &
                          anc_record(k) % item_number,                &
                          Fields_Out, Field_Count_Out, Pos)
      endif

! For Land fields, compress field to land points first

      if (Fields_out(Pos) % stashmaster % grid_type ==                &
                                          ppx_atm_compressed) Then

        do j=1, Fields_out(pos) % levels

          i_full = ancil_add(k)+(j-1) * Output_Grid % loc_p_field
          i_land = ancil_add(k)+(j-1) * Fields_Out(pos) % level_size

! DEPENDS ON: to_land_points
          call to_land_points (                  &
               ancil_data( i_full ),             &
               ancil_data( i_land ),             &
               local_lsm_out,                    &
               Output_Grid % loc_p_field,        &
               local_land_out )

       enddo  !  Loop over j

      endif  !  If land-only field

! ===============================================

! Write the ancillary fields out to the dump

! DEPENDS ON: rcf_writflds
      call Rcf_writflds (Hdr_Out % UnitNum,                        &
                         Fields_out(pos) % levels,                 &
                         Fields_Out(pos) % dump_pos,               &
                         Hdr_Out % Lookup,                         &
                         Hdr_Out % Len1LookUp,                     &
                         ancil_data(ancil_add(k)),                 &
                         Fields_out(pos) % level_size,             &
                         Hdr_Out % FixHd,                          &
                         ErrorStatus,cmessage,.True.)

      if (ErrorStatus /= 0) then
        write (6,*) ' Error in WRITFLDS for Ancillary Field ',K
        write (6,*) ' CMESSAGE ',CMESSAGE
        write (6,*) ' ErrorStatus ',ErrorStatus
        Call Ereport ( RoutineName, ErrorStatus, cmessage)
      endif

    endif  ! if K etc

  endif  !  if ancillary field read in
enddo

deallocate ( lookup_ancil )
deallocate ( nlookup )
deallocate ( lookup_step )
deallocate ( stashancil )
deallocate ( levels )
deallocate ( ancil_add )

deallocate ( fixhd_ancil )
deallocate ( inthd_ancil )
deallocate ( realhd_ancil )

deallocate ( ancil_data )
deallocate ( lookup_start )

! ===============================================

! Write TStar to the output dump

If ( Pos_TStar_Out /= 0 ) Then

  if (mype == 0) then
    if (PrintStatus >= PrStatus_Diag ) Then
      write (6,*) 'Writing TStar to output dump.'
    endif
  endif

  Call Rcf_Write_Field ( TStar_Out, Hdr_Out, Decomp_rcf_output )

Endif

Call Rcf_Dealloc_Field ( TStar_Out )

! ===============================================

Call Rcf_Dealloc_Field ( Land_Frac_Out )

! ===============================================

! Write TStar_land to the output dump

If ( Pos_TStar_Land_Out /= 0 ) Then

  if (mype == 0) then
    if (PrintStatus >= PrStatus_Diag ) Then
      write (6,*) 'Writing TStar_Land to output dump.'
    endif
  endif

  Call Rcf_Write_Field ( TStar_Land_Out, Hdr_Out, Decomp_rcf_output )

Endif

Call Rcf_Dealloc_Field ( TStar_Land_Out )

! ===============================================

! Write TStar_sea to the output dump

If ( Pos_TStar_Sea_Out /= 0 ) Then

  if (mype == 0) then
    if (PrintStatus >= PrStatus_Diag ) Then
      write (6,*) 'Writing TStar_Sea to output dump.'
    endif
  endif

  Call Rcf_Write_Field ( TStar_Sea_Out, Hdr_Out, Decomp_rcf_output )

Endif

Call Rcf_Dealloc_Field ( TStar_Sea_Out )

! ===============================================

! Write TStar_sice to the output dump

If ( Pos_TStar_Sice_Out /= 0 ) Then

  if (mype == 0) then
    if (PrintStatus >= PrStatus_Diag ) Then
      write (6,*) 'Writing TStar_Sice to output dump.'
    endif
  endif

  Call Rcf_Write_Field ( TStar_Sice_Out, Hdr_Out, Decomp_rcf_output )

Endif

Call Rcf_Dealloc_Field ( TStar_Sice_Out )

! ===============================================
! ===============================================

! Write Ice Fraction to the dump

If ( Pos_IceFrac_Out /= 0 ) Then

  if (mype == 0) then
    if (PrintStatus >= PrStatus_Diag ) Then
      write (6,*) 'Writing Ice Fraction to output dump.'
    endif
  endif

  Call Rcf_Write_Field ( IceFrac_Out, Hdr_Out, Decomp_rcf_output )

Endif

call Rcf_Dealloc_Field ( IceFrac_Out )

! ===============================================

! Write T* Anomoly to the dump

If ( Pos_TStar_Anom /= 0 ) Then

  if (mype == 0) then
    if (PrintStatus >= PrStatus_Diag ) Then
      write (6,*) 'Writing TStar Anomoly to output dump.'
    endif
  endif

  Call Rcf_Write_Field ( TStar_Anom, Hdr_Out, Decomp_rcf_output )

Endif

call Rcf_Dealloc_Field ( TStar_Anom )

! ===============================================

! Deallocate workspace for ancilmaster records

Deallocate (anc_record)
Deallocate (anc_file)

Return
End Subroutine Rcf_Ancil_Atmos
End Module Rcf_Ancil_Atmos_Mod
