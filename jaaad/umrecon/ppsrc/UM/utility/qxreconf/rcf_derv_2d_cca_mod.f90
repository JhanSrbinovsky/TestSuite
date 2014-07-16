
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Deriving 2d convective cloud amounts from 2d.

Module Rcf_Derv_2D_CCA_Mod

!  Subroutine Rcf_Derv_2D_CCA
!
! Description:
!    Top level routine obtaining variables (and some conversions)
!    for the Rcf_Calc_2D_CCA routine.
!    To avoid a 3D interpolation, calculation is done with
!    input grid variables and the result is interpolated.
!    Some correction of ccb and cct is then required.
!
! Method:
!
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   30/08/01   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Derv_2D_CCA( fields_in, fields_out, field_count_in, &
                            field_count_out, hdr_in, hdr_out, cca_2d)

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,           &
    Output_Grid

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output,   &
    decomp_rcf_input

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   &
    interp_copy

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

Use Rcf_Interpolate_Mod, Only : &
    Rcf_Interpolate

Use Rcf_Calc_2d_CCA_Mod, Only : &
    Rcf_Calc_2d_CCA

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_3d_cca,          &
    stashcode_ccb,             &
    stashcode_cct,             &
    stashcode_prog_sec

Implicit None

! Arguments
Type( field_type ), Pointer        :: fields_in(:)
Type( field_type ), Pointer        :: fields_out(:)
Type( field_type ), Intent(InOut)  :: cca_2d
Type( um_header_type ), Intent(In) :: hdr_in
Type( um_header_type ), Intent(In) :: hdr_out
Integer, Intent(In)                :: field_count_in
Integer, Intent(In)                :: field_count_out

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

! Local variables
Integer                       :: i           ! looper
Integer                       :: pos         ! field position

Type( field_type ), Pointer   :: cca_3d
Type( field_type ), Pointer   :: ccb
Type( field_type ), Pointer   :: cct
Type( field_type ), Pointer   :: ccb_in
Type( field_type ), Pointer   :: cct_in
Type( field_type )            :: dummy
Type( field_type )            :: cca_2d_in

! Nullify field data
Nullify( dummy % Data )
Nullify( dummy % Data_Int )
Nullify( dummy % Data_Log )
Nullify( cca_2d_in % Data )
Nullify( cca_2d_in % Data_Int )
Nullify( cca_2d_in % Data_Log )

!--------------------------------------------------------------
! Write out out action if appropriate
!--------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Initialising 2D CCA'
End If

!--------------------------------------------------------------
! Find and setup 3D CCA from input
!--------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_3d_cca,               &
                 fields_in, field_count_in, pos )
cca_3d => fields_in(pos)
Call Rcf_Alloc_Field( cca_3d )
Call Rcf_Read_Field( cca_3d, hdr_in, decomp_rcf_input )

!---------------------------------------------------------------
! Set up the cca_2d_in field
!---------------------------------------------------------------
Call rcf_field_equals( cca_2d_in, cca_2d )
cca_2d_in % rows            = cca_3d % rows
cca_2d_in % row_len         = cca_3d % row_len
cca_2d_in % level_size      = cca_3d % level_size
cca_2d_in % glob_rows       = cca_3d % glob_rows
cca_2d_in % glob_row_len    = cca_3d % glob_row_len
cca_2d_in % glob_level_size = cca_3d % glob_level_size

Call Rcf_Alloc_Field( cca_2d_in )
!---------------------------------------------------------------
! Read conv cloud base and top as these are required
!---------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_ccb,                  &
                 fields_in, field_count_in, pos )
ccb_in => fields_in(pos)
Call Rcf_Alloc_Field( ccb_in )
Call Rcf_Read_Field( ccb_in, hdr_in, decomp_rcf_input )

Call Rcf_Locate( stashcode_prog_sec, stashcode_cct,                  &
                 fields_in, field_count_in, pos )
cct_in => fields_in(pos)
Call Rcf_Alloc_Field( cct_in )
Call Rcf_Read_Field( cct_in, hdr_in, decomp_rcf_input )

!-----------------------------------------------------------------
! Have all fields we require - call the calculation routine
!-----------------------------------------------------------------
Call Rcf_Calc_2D_CCA( cca_3d, ccb_in, cct_in, cca_2d_in )

!--------------------------------------------------------------
! Get rid of the fields we no longer require
!--------------------------------------------------------------
Call Rcf_Dealloc_Field( cca_3d )
Call Rcf_Dealloc_Field( ccb_in )
Call Rcf_Dealloc_Field( cct_in )


!--------------------------------------------------------------
! We now need to interpolate cca_2d_in to the output grid
!--------------------------------------------------------------
If (h_int_active) Then
  cca_2d_in % interp = interp_h_only
Else
  cca_2d_in % interp = interp_copy
End If

Call Rcf_Interpolate( cca_2d_in, cca_2d, input_grid, output_grid, &
                      dummy, dummy)

!----------------------------------------------------------------
! Interpolation can mean that ccb, cct and cca_2d are not
! consistent. A quick `clean up' should prevent later crashes.
!----------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_ccb,                  &
                 fields_out, field_count_out, pos )
ccb => fields_out(pos)
Call Rcf_Alloc_Field( ccb )
Call Rcf_Read_Field( ccb, hdr_out, decomp_rcf_output )

Call Rcf_Locate( stashcode_prog_sec, stashcode_cct,                  &
                 fields_out, field_count_out, pos )
cct => fields_out(pos)
Call Rcf_Alloc_Field( cct )
Call Rcf_Read_Field( cct, hdr_out, decomp_rcf_output )

Do i = 1, cca_2d % level_size
  ! Make sure cca is +ive
  If (cca_2d % Data(i,1) < 0. ) Then
    cca_2d % Data(i,1) = 0.
  End If

  ! Make sure cca only exists where ccb and cct are both non-zero
  If ( ccb % Data_Int(i,1) == 0 .AND.     &
       cct % Data_Int(i,1) == 0 .AND.     &
       cca_2d % Data(i,1) > 0. ) Then
    cca_2d % Data(i,1) = 0.0
  End If
End Do

!-------------------------------------------------------------
! Write out the required fields and tidy up memory
!-------------------------------------------------------------
Call Rcf_Write_Field( ccb, hdr_out, decomp_rcf_output )
Call Rcf_Write_Field( cct, hdr_out, decomp_rcf_output )

Call Rcf_Dealloc_Field( cca_2d_in )
Call Rcf_Dealloc_Field( ccb )
Call Rcf_Dealloc_Field( cct )

Return
End Subroutine Rcf_Derv_2D_CCA
End Module Rcf_Derv_2D_CCA_Mod
