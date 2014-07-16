
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ top level wrapper for interpolation

Module Rcf_interpolate_mod

!  Subroutine Rcf_Interpolate - top level interpolation wrapper
!
! Description:
! This module contains a top-level wrapper subroutine for
! interpolation (both horizontal and vertical). It handles
! conversion between data-types (log,int and real) .
! It is worth noting that the orography fields that are fed in
! need not exist if only horizontal interpolation is done.
!
! Method:
!  Intermediate temporary fields are set up - heights are generated,
!  horizontal and then vertical interpolation is done - and then
!  data is reconverted as required.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   05/10/01   Cater for extra rho level at top. D. Robinson
!   5.4   12/08/02   Fix small bug with height_gen_method for input
!                    dumps with 'Original' and Output 'Smooth' R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_interpolate( field_in, field_out, grid_in, grid_out, &
                            interp_orography_in, orography_out )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Grid_Type_Mod, Only : &
    Grid_type

Use Rcf_horizontal_mod, Only : &
    Rcf_horizontal

Use Rcf_vertical_mod, Only : &
    Rcf_vertical

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,             &
    Rcf_Dealloc_Field

Use Rcf_generate_heights_mod, Only : &
    Rcf_generate_heights

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output

Use Rcf_Vert_Cloud_Mod, Only : &
    Rcf_Vert_Cloud

Use Rcf_Set_Interp_Flags_Mod, Only :                         &
    interp_v_only,          interp_h_only,                   &
    interp_all,             interp_done,                     &
    interp_copied,          interp_copy

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_cct,             &
    stashcode_ccb,             &
    stashcode_prog_sec

Use Rcf_Change_Decomposition_Mod, Only : &
    Rcf_Change_Decomposition

Implicit None

! Arguments
Type (field_type), Intent(InOut) :: field_in
Type (field_type), Intent(InOut) :: field_out
Type (grid_type), Intent(In)     :: grid_in
Type (grid_type), Intent(In)     :: grid_out

! This argument (pos dummy) is for the *FINAL OUTPUT* orography
Type (field_type), Intent(In)    :: orography_out

! This argument (pos dummy) is for the *INTERPOLATED INPUT* orography
Type (field_type), Intent(In)    :: interp_orography_in

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

! Local variables/parameters
Type (grid_type)                 :: grid_middle
Type (field_type)                :: field_middle
Character (Len=*), Parameter     :: RoutineName = 'Rcf_Interpolate'
Character (Len=80)               :: Cmessage
Integer                          :: ErrorStatus
Integer                          :: i
Integer                          :: j
Integer                          :: log_to_int
Real, Allocatable                :: heights_middle( :, : )
Real, Allocatable                :: heights_out( :, : )
Logical                          :: vertical_required
Logical                          :: horizontal_required
Logical                          :: interp_required


Nullify( field_middle % Data )
Nullify( field_middle % Data_Int )
Nullify( field_middle % Data_Log )

vertical_required   = (field_in % interp == interp_v_only .OR. &
                       field_in % interp == interp_all )
horizontal_required = (field_in % interp == interp_h_only .OR. &
                       field_in % interp == interp_all )
interp_required     = (vertical_required .OR. horizontal_required)

!------------------------------------------------------------
! Tests to see if orography is set should it be required
!------------------------------------------------------------
If ( vertical_required ) Then
  If ( .NOT. Associated( interp_orography_in % Data ) .OR. &
            .NOT. Associated( orography_out % Data ) ) Then
    Cmessage = 'Orography data not present when vertical &
               &interpolation required'
    ErrorStatus = 10
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If
End If

!------------------------------------------------------------
! Test to see if Data is associated where necessary
!------------------------------------------------------------
If ( (.NOT. Associated( field_in % Data  )       .AND. &
      .NOT. Associated( field_in % Data_Int  )   .AND. &
      .NOT. Associated( field_in % Data_Log  ) ) .OR.  &
     (.NOT. Associated( field_out % Data )       .AND. &
      .NOT. Associated( field_out % Data_Int )   .AND. &
      .NOT. Associated( field_out % Data_Log ) ) ) Then
  Cmessage = 'An interpolation field has not had Data space allocated!'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If


!------------------------------------------------------------------
! Do the interpolation in the relevant way.
!------------------------------------------------------------------
! Force convective cloud base and convective cloud top to have
! lv_code = ppx_rho_level for the duration of this routine to
! generate correct height field
! Note that have to check section as well as item number
! THIS IS DANGEROUS AS IT ALTERS THE INTERNAL STASHMASTER!!!!
If ((field_in % stashmaster % item == stashcode_cct  .OR.    &
     field_in % stashmaster % item == stashcode_ccb) .AND.   &
     field_in % stashmaster % section == stashcode_prog_sec) Then
  field_in  % stashmaster % lv_code = ppx_rho_level
  field_out % stashmaster % lv_code = ppx_rho_level
End If

! middle grid/field have out horizont resolution and in vert.
grid_middle = grid_out     ! horizontal resolutions correct
grid_middle % model_levels = grid_in % model_levels
grid_middle % wet_levels   = grid_in % wet_levels
grid_middle % cloud_levels = grid_in % cloud_levels
grid_middle % st_levels    = grid_in % st_levels
grid_middle % sm_levels    = grid_in % sm_levels
grid_middle % bl_levels    = grid_in % bl_levels
grid_middle % ozone_levels = grid_in % ozone_levels
grid_middle % tr_levels    = grid_in % tr_levels
grid_middle % z_top_of_model = grid_in % z_top_of_model
grid_middle % first_constant_r_rho_level =                          &
                               grid_in % first_constant_r_rho_level
grid_middle % height_gen_method = grid_in % height_gen_method

grid_middle % eta_theta_levels => grid_in % eta_theta_levels
grid_middle % eta_rho_levels   => grid_in % eta_rho_levels
grid_middle % soil_depths      => grid_in % soil_depths

field_middle % levels          = field_in  % levels
field_middle % interp          = field_in  % interp
field_middle % rows            = field_out % rows
field_middle % row_len         = field_out % row_len
field_middle % level_size      = field_out % level_size
field_middle % glob_rows       = field_out % glob_rows
field_middle % glob_row_len    = field_out % glob_row_len
field_middle % glob_level_size = field_out % glob_level_size
field_middle % stashmaster     => field_out % stashmaster

! Allocate field_middle Data space
Call Rcf_Alloc_Field( field_middle )

! Allocate heights - need for vertical function call
Allocate( heights_middle( field_middle % level_size,            &
                      0 :  grid_middle % model_levels + 1) )


Allocate( heights_out( field_out % level_size,                  &
                        0 :  grid_out % model_levels + 1) )


! Only need to generate heights if vertical interpolation is actually
! done
If ( Associated( interp_orography_in % Data ) .AND.                   &
                           vertical_required) Then

  Call rcf_generate_heights( grid_middle, interp_orography_in % Data, &
                         field_middle % stashmaster % grid_type,      &
                         field_middle % stashmaster % lv_code,        &
                         heights_middle,                              &
                         field_middle % level_size )

  Call rcf_generate_heights( grid_out, orography_out % Data,      &
                         field_out % stashmaster % grid_type,     &
                         field_out % stashmaster % lv_code,       &
                         heights_out,                             &
                         field_out % level_size )
End If

!------------------------------------------------------------------
! If data is integer or logical, we will need to convert it to
! Real to do any interpolation on it. Only need to do the conversion
! if some interpolation is actually required!
!------------------------------------------------------------------
If ( interp_required ) Then

  If (field_in % stashmaster % data_type == ppx_type_int ) Then
    Allocate( field_in % Data( field_in % level_size,       &
                               field_in % levels ) )
    Allocate( field_out % Data( field_out % level_size,     &
                                field_out % levels ) )
    Allocate( field_middle % Data( field_middle % level_size,  &
                                   field_middle % levels ) )

    field_in % Data(:,:) = Real( field_in % Data_Int(:,:) )

  Else If (field_in % stashmaster % data_type == ppx_type_log ) Then

    Allocate( field_in % Data( field_in % level_size,       &
                               field_in % levels ) )
    Allocate( field_out % Data( field_out % level_size,     &
                                field_out % levels ) )
    Allocate( field_middle % Data( field_middle % level_size, &
                                   field_middle % levels ) )
    Do i = 1, field_in % level_size
      Do j = 1, field_in % levels
        If ( field_in % Data_Log(i,j) ) Then
          field_in % Data(i,j) = 1.0
        Else
          field_in % Data(i,j) = 0.0
        End If
      End Do
    End Do
  End If
End If

!-------------------------------------------------------------------
! Now can interpolate
!-------------------------------------------------------------------
Call Rcf_horizontal( field_in, field_middle, grid_in, grid_middle )

! Ensure we are on the correct (ie output) decomposition before
! we do vertical interpolation
Call Rcf_Change_Decomposition( decomp_rcf_output )

If( (field_out % stashmaster % item    == stashcode_cct  .OR.  &
     field_out % stashmaster % item    == stashcode_ccb) .AND. &
     field_out % stashmaster % section == stashcode_prog_sec ) Then

  Call Rcf_Vert_Cloud( field_middle, field_out, grid_middle, grid_out, &
                       heights_middle, heights_out)

  ! convert lv_code in stashmaster back to the correct value
  field_in  % stashmaster % lv_code = ppx_single_level
  field_out % stashmaster % lv_code = ppx_single_level
Else
  Call Rcf_vertical( field_middle, field_out, grid_middle, grid_out,&
                     heights_middle, heights_out)
End If

! Deallocate heights
Deallocate( heights_middle )
Deallocate( heights_out )

!-------------------------------------------------------------------
! Deallocate field_middle Data
!-------------------------------------------------------------------
Call Rcf_Dealloc_Field( field_middle )

!-------------------------------------------------------------------
! Need to convert back to integer or logical if field is such. Only
! need to do this if interpolation active!
!-------------------------------------------------------------------
If ( interp_required ) Then
  If (field_out % stashmaster % data_type == ppx_type_int ) Then
    field_out % Data_Int(:,:) = Nint( field_out % Data(:,:) )

    Deallocate( field_out % Data )
    Deallocate( field_in % Data )
    ! middle data already deallocated by dealloc_field above

  Else If (field_out % stashmaster % data_type == ppx_type_log ) Then
    Do i = 1, field_out % level_size
      Do j = 1, field_out % levels
        log_to_int = Int( field_out % Data(i,j) + 0.5 )
        If ( log_to_int == 1 ) Then
          field_out % Data_Log(i,j) =  .TRUE.
        Else If (log_to_int == 0 ) Then
          field_out % Data_Log(i,j) =  .FALSE.
        Else
          Cmessage = 'Conversion from Logical to Integer: Illegal value'
          ErrorStatus = 30
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
      End Do
    End Do

    Deallocate( field_out % Data )
    Deallocate( field_in % Data )
  End If
End If

!--------------------------------------------------------------------
! Set the interp flag for the output field
!--------------------------------------------------------------------
If (interp_required) Then
  field_out % interp = interp_done

Else If (field_in % interp == interp_copy) Then
  field_out % interp = interp_copied

End If

End Subroutine Rcf_interpolate

End Module Rcf_interpolate_mod
