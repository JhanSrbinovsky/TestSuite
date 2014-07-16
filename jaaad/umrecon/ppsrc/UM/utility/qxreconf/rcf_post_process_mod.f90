
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Post `main loop' processing of data

Module Rcf_Post_Process_Mod

!  Subroutine Rcf_Post_Process - process fields after `main loop'
!
! Description:
! This subroutine performs post-processing of data using transforms
! etc that require a number of output-grid fields. This involves
! some transforms to original data types. Performed after the main
! field data creation loop.
!
! Method:
!  The relevant processing is just run through - this may need to be
!  be made more modular in future.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Post_Process_Atmos( fields_in, field_count_in, orog_in,&
                                   fields_out, field_count_out, orog, &
                                   hdr_in, hdr_out, data_source)

Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Rcf_Submodel_Mod, Only : &
    Atmos_IM

Use Rcf_Calc_Rho_mod, Only : &
    Rcf_Calc_Rho

Use Rcf_Set_Interp_Flags_Mod, Only  : &
    interp_done,                      &
    interp_no_op

Use Rcf_Freeze_Soil_Mod, Only : &
    Rcf_Freeze_Soil

Use Rcf_Data_Source_Mod, Only : &
    Data_Source_Type,           &
    Ancillary_File,             &
    Input_Dump

Use Rcf_Conv_Cld_Chk_Mod, Only : &
    Rcf_Conv_Cld_Chk

Use Rcf_Cloud_Frac_Chk_Mod, Only : &
    Rcf_Cloud_Frac_Chk

Use Rcf_Snow_Amount_Chk_Mod, Only : &
    Rcf_Snow_Amount_Chk

Use Rcf_Sea_Ice_Frac_Chk_Mod, Only : &
    Rcf_Sea_Ice_Frac_Chk

Use Rcf_Soil_Moist_Chk_Mod, Only : &
    Rcf_Soil_Moist_Chk

Use Rcf_WriteUMhdr_Mod, Only : &
    Rcf_WriteUMhdr

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only :&
    um_header_type

Use Rcf_Grid_Type_Mod, Only : &
    output_grid

Use Rcf_DecompTP_Mod, Only : &
    Decomp_rcf_output

Use Rcf_HeadAddress_Mod, Only : &
    RC_PressureTop,         &
    RC_LongSpacing,         &
    RC_AtmMoist,            &
    RC_AtmMass,             &
    RC_AtmEnergy,           &
    RC_EnergyCorr

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_generate_heights_mod, Only : &
    Rcf_generate_heights

Use rcf_theta_t_convs_mod, Only : &
    Rcf_conv_theta_t,             &
    Rcf_conv_t_theta

Use Rcf_read_field_mod, Only : &
    Rcf_read_field

Use Rcf_Calc_Output_Exner_Mod, Only : &
    Rcf_Calc_Output_Exner

Use Rcf_Write_field_mod, Only : &
    Rcf_write_field

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_polar_u_Mod, Only : &
    Rcf_polar_u

Use Rcf_Exner_P_Convs_Mod, Only : &
    Rcf_Conv_Exner_P,             &
    Rcf_Conv_P_Exner

Use Rcf_PrintStatus_Mod, Only : &
    LTimer

Use Rcf_Stashcodes_Mod, Only :&
    stashcode_u,              stashcode_v,       &
    stashcode_theta,          stashcode_q,       &
    stashcode_rho,            stashcode_exner,        &
    stashcode_p,              stashcode_soil_moist,   &
    stashcode_soil_temp,      stashcode_orog,         &
    stashcode_prog_sec

Use Rcf_Soilstress_To_Soilmoist_Mod,only:&
    Rcf_Soilstress_to_soilmoist

Use Rcf_Recon_Mod, Only : &
    use_smc_stress

Implicit None

! Arguments
Type( field_type ), Pointer           :: fields_in(:)
Type( field_type ), Pointer           :: fields_out(:)
Type( field_type ), Intent(InOut)     :: orog_in
Type( field_type ), Intent(In)        :: orog
Type( data_source_type ), Pointer     :: data_source(:)
Type( um_header_type ), Intent(In)    :: hdr_in
Type( um_header_type ), Intent(InOut) :: hdr_out
Integer, Intent(In)                   :: field_count_in
Integer, Intent(In)                   :: field_count_out

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
Integer                       :: pos
Integer                       :: pos_st   !} positions in array of
Integer                       :: pos_sm   !} soil temp and moisture
Integer                       :: pos_orog !} and of orography
Integer                       :: ErrorStatus
Character (Len=*), Parameter  :: RoutineName = 'Rcf_Post_Process'
Character (Len=80)            :: Cmessage

Real                          :: theta_heights(                     &
                                    output_grid % loc_p_field,      &
                                    0 : output_grid % model_levels + 1)
                                    ! Height of theta levels
Real                          :: rho_heights(                       &
                                    output_grid % loc_p_field,      &
                                    0 : output_grid % model_levels + 1)
                                    ! Height of rho levels
Logical                       :: exner_out ! is exner in the output?
Logical                       :: l_soil_change   ! is soil changed

Type( field_type ), Pointer   :: u
Type( field_type ), Pointer   :: v
Type( field_type ), Pointer   :: exner
Type( field_type ), Pointer   :: theta
Type( field_type ), Pointer   :: q
Type( field_type ), Pointer   :: rho
Type( field_type ), Pointer   :: p
Type( field_type ), Target    :: pressure  ! exner OR p depending on
                                           ! circumstance

External Timer

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)

Nullify( pressure % Data     )
Nullify( pressure % Data_Int )
Nullify( pressure % Data_Log )

Call Rcf_Locate( stashcode_prog_sec, stashcode_soil_moist,           &
                 fields_out, field_count_out, pos_sm, .TRUE. )

If (use_smc_stress .and.                                             &
   Data_Source( pos_sm ) % Source == Input_Dump) Then
   Call Rcf_Soilstress_To_Soilmoist( fields_out, field_count_out,    &
                                     output_grid,decomp_rcf_output,  &
                                     hdr_out)
End If

!-------------------------------------------------------------------
! Find and setup Theta (will hold T if interpolated)
!-------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_theta,                &
                 fields_out, field_count_out, pos)
theta => fields_out( pos )
Call Rcf_Alloc_Field( theta )
Call Rcf_Read_Field( theta, hdr_out, decomp_rcf_output )

If ( h_int_active .OR. v_int_active ) Then
!-------------------------------------------------------------------
! Calculate required heights
!-------------------------------------------------------------------
  Call rcf_generate_heights( output_grid,orog % Data(1:,1:),   &
                             ppx_atm_tall, ppx_theta_level,    &
                             theta_heights, theta % level_size )

  Call rcf_generate_heights( output_grid,orog % Data(1:,1:),  &
                             ppx_atm_tall, ppx_rho_level,     &
                             rho_heights,  theta % level_size )
End If

!--------------------------------------------------------------------
! Find and read Q
!--------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_q,                    &
                 fields_out, field_count_out, pos )
q => fields_out( pos )
Call Rcf_Alloc_Field( q )
Call Rcf_Read_Field( q, hdr_out, decomp_rcf_output )

!-------------------------------------------------------------------
! Find exner and P (as required)
! Note assumption that will have exner *OR* P available
!-------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_exner,                &
                 fields_out, field_count_out, pos, .TRUE. )

If (pos /= 0) Then
  exner => fields_out( pos )
  Call Rcf_Alloc_Field( exner )
  exner_out = .TRUE.

  ! Also need space for P
  p => pressure
  Call Rcf_Field_Equals( p, exner )
  Call Rcf_Alloc_Field( p )
  p % stashmaster => Rcf_Exppx( Atmos_IM, 0, stashcode_p )

Else    ! exner not in output dump

  Call Rcf_Locate( stashcode_prog_sec, stashcode_p,                  &
                   fields_out, field_count_out, pos )
  p => fields_out( pos )
  Call Rcf_Alloc_Field( p )
  exner_out = .FALSE.

  ! Also need space for exner
  exner => pressure
  Call Rcf_Field_Equals( exner, p )
  Call Rcf_Alloc_Field( exner )
  exner % stashmaster => Rcf_Exppx( Atmos_IM, 0, stashcode_exner )

End If

!--------------------------------------------------------------------
! Either read or calculate exner (as appropriate)
!--------------------------------------------------------------------
If ( exner % interp /= interp_no_op ) Then
  If ( exner_out ) Then
    Call Rcf_Read_Field( exner, hdr_out, decomp_rcf_output )
  Else
    Call Rcf_Read_Field( p, hdr_out, decomp_rcf_output )
    exner % Data(:,:) = p % Data(:,:)
    exner % stashmaster => Rcf_Exppx( Atmos_IM, 0, stashcode_p )
    Call Rcf_Conv_P_Exner( exner )
  End If

Else
  If ( theta % interp /= interp_done ) Then
    ErrorStatus = 10
    Cmessage = 'Only have Theta available - need T to calculate exner'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  Call Rcf_Locate( stashcode_prog_sec, stashcode_orog,               &
                   fields_out, field_count_out, pos_orog )
  Call Rcf_Calc_Output_Exner( fields_in, field_count_in, orog_in,     &
                              hdr_in, orog, theta, q, exner,          &
                              Data_Source( pos_orog ) % source,       &
                              rho_heights, theta_heights )

  If (exner_out) Then
    Call Rcf_Write_Field( exner, hdr_out, decomp_rcf_output )
  End If

  ! Also calculate P - needed internally. Note need to temporarily
  ! fool with stashcode
  p % Data(:,:) = exner % Data(:,:)
  p % stashmaster => Rcf_Exppx( Atmos_IM, 0, stashcode_exner )
  Call Rcf_Conv_Exner_P( p )

  ! Write this newly calcuated field to dump
  If (.NOT. exner_out) Then
    Call Rcf_Write_Field( p, hdr_out, decomp_rcf_output )
  End If

End If

!********************************************************************
! Convert T (in dump) back to theta
!********************************************************************
If ( theta % interp == interp_done ) Then
  Call Rcf_Conv_T_Theta( theta, fields_out, field_count_out, hdr_out, &
                         decomp_rcf_output, rho_heights, theta_heights)

  Call Rcf_Write_Field( theta, hdr_out, decomp_rcf_output )
End If


!*******************************************************************
! Calculate rho
!*******************************************************************
Call Rcf_Locate( stashcode_prog_sec, stashcode_rho,                  &
                 fields_out, field_count_out, pos)
rho => fields_out(pos)
If ( rho % interp == interp_no_op ) Then

  ! No interpolation for rho need to calculate
  Call Rcf_Alloc_Field( rho )

  Call Rcf_Calc_Rho( theta, q, exner, p, theta_heights, rho_heights, &
                     rho)

  Call Rcf_Write_Field( rho, hdr_out, decomp_rcf_output )
  Call Rcf_DeAlloc_Field( rho )

End if

!--------------------------------------------------------------------
! Tidy up
!--------------------------------------------------------------------
Call Rcf_DeAlloc_Field( p )
Call Rcf_DeAlloc_Field( q )
Call Rcf_DeAlloc_Field( exner )
Call Rcf_DeAlloc_Field( theta )

!********************************************************************
! Reset Soil Moisture to match vegetation if interpolated
!********************************************************************
If (h_int_active) Then
  Call Rcf_Soil_Moist_Chk( fields_out, field_count_out, output_grid, &
                           decomp_rcf_output, hdr_out, l_soil_change) 
End If

!********************************************************************
! Adjust frozen/unfrozen soil if updated from ancillary
!********************************************************************
Call Rcf_Locate( stashcode_prog_sec, stashcode_soil_temp,            &
                 fields_out, field_count_out, pos_st, .TRUE. )


! Only need to correct if soil moisture/temperature are both in dump
If (pos_st /= 0 .AND. pos_sm /= 0) Then
  If ( Data_Source( pos_st ) % Source == Ancillary_File  .OR. &
       Data_Source( pos_sm ) % Source == Ancillary_File  .OR. &
       l_soil_change ) Then

    Call Rcf_Freeze_Soil( fields_out, field_count_out, hdr_out, &
                          output_grid )
  End If
End If

!********************************************************************
! Check convective cloud base and top are sensible
!********************************************************************
Call Rcf_Conv_Cld_Chk( fields_out, field_count_out, output_grid,  &
                       decomp_rcf_output, hdr_out )

!********************************************************************
! Check that cloud fraction fields are consistent with q fields.
!********************************************************************
Call Rcf_Cloud_Frac_Chk( fields_out, field_count_out, output_grid,  &
                         decomp_rcf_output, hdr_out )

!********************************************************************
! Check for negative snow amounts
!********************************************************************
If (h_int_active) Then
  Call Rcf_Snow_Amount_Chk( fields_out, field_count_out,             &
                            output_grid,  decomp_rcf_output, hdr_out )
End If

!********************************************************************
! Check for small sea-ice fractions if interpolated
!********************************************************************
If (h_int_active) Then
  Call Rcf_Sea_Ice_Frac_Chk( fields_out, field_count_out,            &
                             decomp_rcf_output, hdr_out, data_source)
End If

!********************************************************************
! If appropriate, set the polar rows for u (as a vector mean)
!********************************************************************
If ( output_grid % global ) Then
  Call Rcf_Locate( stashcode_prog_sec, stashcode_u,                  &
                   fields_out, field_count_out, pos)
  u => fields_out(pos)

  If ( u % interp == interp_done ) Then
    Call Rcf_Alloc_Field( u )
    Call Rcf_Read_Field( u, hdr_out, decomp_rcf_output )

    Call Rcf_Locate( stashcode_prog_sec, stashcode_v,                &
                     fields_out, field_count_out, pos )
    v => fields_out(pos)
    Call Rcf_Alloc_Field( v )
    Call Rcf_Read_Field( v, hdr_out, decomp_rcf_output )

    Call Rcf_Polar_U( u, v, hdr_out % RealC( RC_LongSpacing ) )

    Call Rcf_Write_Field( u, hdr_out, decomp_rcf_output )

    Call Rcf_DeAlloc_Field( u )
    Call Rcf_DeAlloc_Field( v )
  End If
End If

!------------------------------------------------------------
! Want to keep energy correction information if no interpolation 
! and is not an ocean.  This was placed here due to decision
! whether interpolation is on is also done in create_dump
! when orography ancillary can be different to input dump.
!------------------------------------------------------------
If ( .Not. ( V_Int_Active .Or. H_Int_Active ) ) Then
  Hdr_Out % RealC( RC_AtmMoist   ) = Hdr_In % RealC( RC_AtmMoist   )
  Hdr_Out % RealC( RC_AtmMass    ) = Hdr_In % RealC( RC_AtmMass    )
  Hdr_Out % RealC( RC_AtmEnergy  ) = Hdr_In % RealC( RC_AtmEnergy  )
  Hdr_Out % RealC( RC_EnergyCorr ) = Hdr_In % RealC( RC_EnergyCorr )

!---------------------------------------------------------------
! Write out the header here again due to change above.
!---------------------------------------------------------------
  Call Rcf_WriteUMhdr( Hdr_Out )

End If

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

Return
End Subroutine Rcf_Post_Process_Atmos
End Module Rcf_Post_Process_Mod
