
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Dump creation `main loop'

Module Rcf_Create_Dump_Mod

!  Subroutine Rcf_Create_Dump - main dump creation loop.
!
! Description:
!   Creates the output dump based on the input choices.
!
! Method:
!   Data sources are setup and then fields in the output dump are
!   looped over. Each one is set appropriately ( interpolated from
!   input dump, set to 0/MDI/constant/from file/etc). Pre and Post
!   processing is performed as required.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Create_Dump( Hdr_In, Hdr_Out, fields_in, fields_out, &
                            field_count_in, field_count_out )

Use Rcf_Set_Orography_Mod, Only : &
    Rcf_Set_Orography

Use Rcf_Pre_Interp_Transform_Mod, Only : &
    Rcf_Pre_Interp_Transform

Use Rcf_Post_Interp_Transform_Mod, Only : &
    Rcf_Post_Interp_Transform

Use Rcf_Post_Process_Mod, Only :&
    Rcf_Post_Process_Atmos

Use Rcf_Recon_Mod, Only : &
    Trans,            &
    Uars,             &
    GRIB

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Get_Unit,             &
    Rcf_Free_Unit

Use Rcf_FreeUMhdr_mod, Only : &
    Rcf_FreeUMhdr

Use Rcf_UMhead_Mod, Only : &
    UM_header_type

Use Ereport_mod, Only : &
    Ereport

Use Rcf_Data_Source_Mod    ! All of it...

Use rcf_read_field_mod, Only : &
    Rcf_Read_Field

Use rcf_write_field_mod, Only : &
    Rcf_Write_Field

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_DecompTP_Mod, Only : &
    Decomp_rcf_input,    &
    Decomp_rcf_output

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,           &
    Output_Grid

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_field_equals_mod, Only : &
    Rcf_field_equals

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active,         &
    v_int_active_soil

Use rcf_interpolate_mod, Only : &
    rcf_interpolate

Use rcf_set_data_source_mod, Only : &
    Rcf_Set_Data_Source

Use Rcf_Lsm_Mod, Only : &
    local_lsm_out

Use Rcf_aux_file_mod, Only : &
    tracers,                  uars_data,      &
    user_prog,                transplant,     &
    rcf_aux_file

Use Rcf_Field_Calcs_Mod, Only : &
    Rcf_Field_Calcs

Use Rcf_Set_Interp_Flags_Mod, Only :               &
    Rcf_Set_Interp_Flags,                          &
    interp_v_only,            interp_h_only,       &
    interp_all,               interp_no_op,        &
    interp_copy

Use Rcf_Rotate_Mod, Only : &
    Rcf_Rotate,            &
    ToStandard,            &
    FromStandard

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_lsm,             &
    stashcode_orog,            &
    stashcode_icefrac,         &
    stashcode_tstar,           &
    stashcode_tstar_land,      &
    stashcode_tstar_sea,       &
    stashcode_tstar_sice,      &
    stashcode_prog_sec

Use Rcf_GRIB_Interp_TnPstar_Mod, Only : &
    Rcf_GRIB_Interp_TnPstar

Implicit None

! Arguments
Integer, Intent(In)                       :: field_count_in
Integer, Intent(In)                       :: field_count_out
Type (UM_Header_Type), Intent(InOut)      :: Hdr_In
Type (UM_Header_Type), Intent(InOut)      :: Hdr_Out
Type (field_type), Pointer                :: fields_in(:)
Type (field_type), Pointer                :: fields_out(:)

! Comdecks
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

! Local vars
Integer                          :: pos      ! position in fields array
Integer                          :: i        ! looper
Integer                          :: err      ! for file open/close
Integer                          :: IOStatus
Integer                          :: ErrorStatus
Logical                          :: l_exist
Character (Len=80)               :: Cmessage
Character (Len=Max_Filename_Len) :: TracerFileName
Character (Len=*), Parameter     :: RoutineName='Rcf_Create_Dump'

Type (field_type), Target  :: interp_orog   ! interpolated input orog
Type (field_type), Pointer :: orog_out      ! ptr to output orog
Type (field_type), Pointer :: orog_in       ! ptr to intput orog
Type (field_type), Pointer :: input_field   ! ptr to input field
Type (field_type), Pointer :: output_field  ! ptr to output field
Type (UM_Header_Type)      :: Hdr_Aux       ! auxillary file header

Type (data_source_type), Pointer :: data_source(:)

! Formatting
Character (Len=*), Parameter :: form="(a20, i3, ' ( Section', i4, ' )', &
                                      &' ( Stashcode', i4, ' )', ' ', a36)"

External File_Open, File_Close, Fort_Get_Env
!---------------------------------------------------------------
! Initialise some data
!---------------------------------------------------------------
Nullify( interp_orog % data )
Nullify( interp_orog % data_int )
Nullify( interp_orog % data_log )
Nullify( data_source )

! set the source of the data for output dump
Call Rcf_Set_Data_Source( data_source, fields_in, fields_out,  &
                          field_count_in, field_count_out,     &
                          hdr_in, hdr_out )

!-----------------------------------------------------------
! Set up the input, output and interpolated orography fields
!-----------------------------------------------------------
Call Rcf_Set_Orography( fields_in, fields_out, field_count_in,      &
                        field_count_out, hdr_in, hdr_out,           &
                        data_source, orog_in, orog_out, interp_orog )

! If recon from GRIB T & Pstar need to interp'd horizontally
! for height generation
If ( GRIB ) Then
  Call Rcf_GRIB_Interp_TnPstar( fields_in, fields_out, Hdr_In,      &
                                  field_count_in, field_count_out)
End If

!-------------------------------------------------------------
! Land Sea Mask - need to deal with this a little differently
! as data already read in and processed.
!-------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_lsm,                  &
                 fields_out, field_count_out, pos, .TRUE.)
If ( pos /= 0 ) Then
  Call Rcf_Alloc_Field( fields_out ( pos ) )
  fields_out( pos ) % Data_Log(:,1) = local_lsm_out(:)

  Call Rcf_Write_Field( fields_out( pos ), Hdr_Out, decomp_rcf_output )
  Call Rcf_DeAlloc_Field( fields_out( pos ) )
End If

!-------------------------------------------------------------
! Setup interpolation flags
!-------------------------------------------------------------
Call Rcf_Set_Interp_Flags( fields_in, fields_out, field_count_in, &
                           field_count_out, data_source )

!------------------------------------------------------------
! Rotate input winds if so required
!------------------------------------------------------------
If ( Input_Grid % Rotated .AND. h_int_active ) Then
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
    Write (6,*) 'Rotating input winds'
  End If

  Call Rcf_Rotate( fields_in, field_count_in, Input_Grid, Hdr_In,  &
                   decomp_rcf_input, ToStandard)
End If

!-----------------------------------------------------
! Main loop
!-----------------------------------------------------
Do i = 1, field_count_out ! run through the output fields

  !--------------------------------------------
  ! Which fields have we already dealt with?
  ! Orography
  ! LSM
  ! T*             )  if source == 4
  ! Ice fraction   )  if source == 4
  !--------------------------------------------
  If ( fields_out( i ) % stashmaster % section == stashcode_prog_sec .AND. &
       ( fields_out( i ) % stashmaster % item == stashcode_orog .OR.       &
         fields_out( i ) % stashmaster % item == stashcode_lsm  .OR.       &
         (fields_out( i ) % stashmaster % item == stashcode_tstar .AND.    &
          data_source( i ) % source == Set_to_MDI )              .OR.      &
         (fields_out( i ) % stashmaster % item == stashcode_icefrac .AND.  &
          data_source( i ) % source == Set_to_MDI ) ) .OR.                 &
       data_source( i ) % source == Already_Processed ) Then

    If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
      Write (6,form) 'Already processed ', i,                       &
                      fields_out( i ) % stashmaster % section,      &
                      fields_out( i ) % stashmaster % item,         &
                      fields_out( i ) % stashmaster % name
    End If

    Cycle    ! just skip this iteration of the loop.

  End If

  ! Set the output data space - don't actually use this in the case
  ! of external files used for the data, but will do it anyway to
  ! simplify the algorithm

  output_field => fields_out( i )
  Call Rcf_Alloc_Field( output_field )

  Select Case( Data_Source( i ) % Source )

!---------------------------------------------------------------
! Data interpolated (or copied) from input dump
!---------------------------------------------------------------
    Case( Input_Dump )

      ! set up input fields
      Call Rcf_Locate( fields_out( i ) % stashmaster % section,    &
                       fields_out( i ) % stashmaster % item,       &
                       fields_in, field_count_in, pos )
      input_field => fields_in( pos )
      Call Rcf_Alloc_Field( input_field )

      Call Rcf_Read_Field( input_field, Hdr_In, decomp_rcf_input )

      ! Write out appropriate message
      Select Case( input_field % interp )

        Case( interp_h_only, interp_v_only, interp_all)
          If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
            Write (6,form) 'Interpolating Field', i,                &
                           input_field % stashmaster % section,     &
                           input_field % stashmaster % item,        &
                           input_field % stashmaster % name
          End If

        Case( interp_copy )
          If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
            Write (6,form) 'Copying Field', i,                      &
                           input_field % stashmaster % section,     &
                           input_field % stashmaster % item,        &
                           input_field % stashmaster % name
          End If

        Case( interp_no_op )
          If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
            Write (6,form) 'Skipping Field', i,                     &
                           input_field % stashmaster % section,     &
                           input_field % stashmaster % item,        &
                           input_field % stashmaster % name
          End If

          Call Rcf_DeAlloc_Field( input_field )
          Call Rcf_DeAlloc_Field( output_field )
          Cycle

       End Select

      ! convert fields approriately for interpolation
      Call Rcf_Pre_Interp_Transform( input_field, fields_in,     &
                                     field_count_in, hdr_in, orog_in )

      Call Rcf_Interpolate( input_field, output_field, Input_Grid, &
                              Output_Grid, interp_orog, orog_out )

      ! Convert fields back to original form (if possible) and
      ! perform any simple post-processing required
      Call Rcf_Post_Interp_Transform( output_field, fields_out,  &
                                      field_count_out )

      Call Rcf_DeAlloc_Field( input_field )

!------------------------------------------------------------------
! Data from Ancillary File
!------------------------------------------------------------------
    Case( Ancillary_File )

      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Ancillary field', i,                       &
                    fields_out( i ) % stashmaster % section,       &
                    fields_out( i ) % stashmaster % item,          &
                    fields_out( i ) % stashmaster % name
      End If

!------------------------------------------------------------------
! Data to be set to zero
!------------------------------------------------------------------
    Case( Set_To_Zero )

      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Set to Zero, field', i,                    &
                    output_field % stashmaster % section,          &
                    output_field % stashmaster % item,             &
                    output_field % stashmaster % name
      End If

      Select Case ( output_field % stashmaster % data_type )
        Case ( ppx_type_real )
          output_field % Data( :, : ) = 0.0

        Case ( ppx_type_int )
          output_field % Data_Int( : , : ) = 0

        Case ( ppx_type_log )
          output_field % Data_Log( : , : ) = .FALSE.
          If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
            Write (6,*) 'Assuming zero means FALSE for logicals'
          End If

        Case Default
          Write (6,*) 'Cannot set fields of this type to zero!'
          Write (6,*) 'ppx_type_? =',&
            output_field % stashmaster % data_type
      End Select

!------------------------------------------------------------------
! Data to be set as missing
!------------------------------------------------------------------
    Case( Set_To_MDI )

      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Set to MDI ', i,                          &
                    output_field % stashmaster % section,         &
                    output_field % stashmaster % item,            &
                    output_field % stashmaster % name
      End If

      Select Case ( output_field % stashmaster % data_type )
        Case ( ppx_type_real )
          output_field % Data( :, : ) = RMDI

        Case ( ppx_type_int )
          output_field % Data_Int( : , : ) = IMDI

        Case Default
          Write (6,*) 'Cannot set fields of this type to MDI!'
      End Select

!----------------------------------------------------------------
! Tracer data
     !!!kdcorbin, 03/10 - Replaced using ATRACER file with a 
     !!!   user prognostic ancillary file
!----------------------------------------------------------------
    Case( Tracer_File )

      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Tracer file data ', i,                   &
                    output_field % stashmaster % section,        &
                    output_field % stashmaster % item,           &
                    output_field % stashmaster % name
      End If

      ! Get Tracer File name from env var ATRACER
      !Call Fort_Get_Env( 'ATRACER', 7, TracerFileName, &
      !                   Max_Filename_Len, err )

      !If ( err /= 0 ) Then
      !  Write (6,*) 'Cannot get Tracer File name from Env Var ATRACER '
      !  cmessage = 'Cannot get Tracer File name from Env Var ATRACER'
      !  Call Ereport ( RoutineName, err, cmessage)
      !Endif

      TracerFileName=data_source(i)%Ancil_File

!     Check that Tracer File exists
      Inquire ( file=TracerFileName, exist=l_exist, iostat=IOStatus )

      If ( .not. l_exist ) then
        Write (6,*) 'Tracer File does not exist.'
        Write (6,*) 'File : ',TracerFileName
        ErrorStatus=10
        Cmessage = 'Tracer File does not exist.'
        Call Ereport ( RoutineName, ErrorStatus, Cmessage )
      End If

      Call Rcf_Get_Unit( Hdr_Aux % UnitNum )
! DEPENDS ON: file_open
      !Replaced ATRACER - kdcorbin, 03/10
      Call File_Open( Hdr_Aux % UnitNum, TracerFileName,Max_Filename_Len , &
             0 ,1 ,err )

      !Replaced IMDI with info in last two spots - kdcorbin, 08/10
      Call Rcf_Aux_File( Hdr_Aux, Hdr_Out, fields_out, field_count_out,&
                         tracers, output_field % stashmaster % section,&
                         output_field % stashmaster % item, &
                         output_field%stashmaster%section,  &
                         output_field%stashmaster%item)

! DEPENDS ON: file_close
      !Replace ATRACER with TracerFileName - kdcorbin, 03/10
      Call File_Close( Hdr_Aux % UnitNum,trim(TracerFileName), &
                   Max_Filename_Len ,1 ,0 ,err )
      Call Rcf_Free_Unit( Hdr_Aux % UnitNum)
      Call Rcf_FreeUMhdr( Hdr_Aux )

!----------------------------------------------------------------
! Data to be set to constant from the namelist
!----------------------------------------------------------------
    Case( Set_To_Const )
      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Set to user const', i,                    &
                    output_field % stashmaster % section,         &
                    output_field % stashmaster % item,            &
                    output_field % stashmaster % name
      End If

      Select Case ( fields_out( i ) % stashmaster % data_type )
        Case ( ppx_type_real )
          output_field % Data( :, : ) = data_source( i ) % RConst

        Case ( ppx_type_int )
          output_field % Data_Int( : , : ) =                         &
                                     Nint( data_source( i ) % RConst )
        Case ( ppx_type_log )
          If ( data_source( i ) % Rconst > 0.5 ) Then
            output_field % Data_Log( : ,: ) = .TRUE.
          Else
            output_field % Data_Log( :, : ) = .FALSE.
          End If

      End Select

!----------------------------------------------------------------
! User prognostics from external dump
!----------------------------------------------------------------
    Case( External_Dump )

      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'User Prognostic', i,                     &
                    output_field % stashmaster % section,        &
                    output_field % stashmaster % item,           &
                    output_field % stashmaster % name
      End If

!     Check that ancillary file exists
      Inquire ( file=data_source( i ) % Ancil_File, &
               exist=l_exist, iostat=IOStatus )

      If ( .not. l_exist ) then
        Write (6,*) 'User Prognostic File does not exist.'
        Write (6,*) 'File : ',data_source( i ) % Ancil_File
        ErrorStatus=10
        Cmessage = 'User Prognostic File does not exist.'
        Call Ereport ( RoutineName, ErrorStatus, Cmessage )
      End If

      Call Rcf_Get_Unit( Hdr_Aux % UnitNum )
! DEPENDS ON: file_open
      Call File_Open( Hdr_Aux % UnitNum, data_source( i ) % Ancil_File,&
                        Max_Filename_Len ,0 ,1 ,err )

      Call Rcf_Aux_File( Hdr_Aux, Hdr_Out, fields_out, field_count_out,&
                     user_prog, data_source( i ) % Ancil_SctnC,        &
                     data_source( i ) % Ancil_ItemC,                   &
                     output_field % stashmaster % section,             &
                     output_field % stashmaster % item)

! DEPENDS ON: file_close
      Call File_Close(Hdr_Aux % UnitNum, data_source( i ) % Ancil_File,&
                        Max_Filename_Len ,1 ,0 ,err )
      Call Rcf_Free_Unit( Hdr_Aux % UnitNum)
      Call Rcf_FreeUMhdr( Hdr_Aux )

!-----------------------------------------------------------------
! Calculations for fields that are missing in the input dump.
! These are skipped now and done in a later loop
!-----------------------------------------------------------------
    Case( Field_Calcs)
      If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Calcs. needed, skip', i,                 &
                    output_field % stashmaster % section,        &
                    output_field % stashmaster % item,           &
                    output_field % stashmaster % name
      End If

!-----------------------------------------------------------------
! Unrecognised source for data - set it to missing
!-----------------------------------------------------------------
    Case Default

      ErrorStatus = -10
      Cmessage = 'Source code not recognised - will set field to MDI'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )

      If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,form) 'Unknown source', i,                      &
                  & output_field % stashmaster % section,        &
                  & output_field % stashmaster % item, 'Set to MDI!'
      End If

      Select Case ( fields_out(i) % stashmaster % data_type )
        Case ( ppx_type_real )
          output_field % Data( :, : ) = RMDI

        Case ( ppx_type_int )
          output_field % Data_Int( : , : ) = IMDI

        Case ( ppx_type_log )
          output_field % Data_Log( :, : ) = .FALSE.
      End Select


  End Select

!-------------------------------------------------------------------
! Write out data if required
!-------------------------------------------------------------------
  If ( Data_Source( i ) % Source /= External_Dump  .AND.  &
       Data_Source( i ) % Source /= Tracer_File    .AND.  &
       Data_Source( i ) % Source /= Field_Calcs    .AND.  &
       Data_Source( i ) % Source /= Ancillary_File ) Then

    Call Rcf_Write_Field( output_field, Hdr_Out, decomp_rcf_output )

  End If

  Call Rcf_DeAlloc_Field( output_field )
End Do

!------------------------------------------------------------------
! Perform post-processing
!------------------------------------------------------------------
Call Rcf_Post_Process_Atmos( fields_in, field_count_in, orog_in,    &
                             fields_out, field_count_out, orog_out, &
                             hdr_in, hdr_out, data_source)

!------------------------------------------------------------------
! Tidy up orography fields
!------------------------------------------------------------------
Call Rcf_DeAlloc_Field( orog_in )
Call Rcf_DeAlloc_Field( orog_out )
Call Rcf_DeAlloc_Field( interp_orog )

!------------------------------------------------------------------
! Need to revisit the source = field_calcs (8) fields
!------------------------------------------------------------------

Call Rcf_Field_Calcs( fields_in, fields_out, field_count_in, &
                      field_count_out, data_source, hdr_in, hdr_out )

!-------------------------------------------------------------------
! Rotate the winds if so required
!-------------------------------------------------------------------
If (h_int_active .AND. Output_Grid % Rotated) Then
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
    Write (6,*) 'Rotating Output Winds'
  End If

  Call Rcf_Rotate( fields_out, field_count_out, Output_Grid, Hdr_Out, &
                   decomp_rcf_output, FromStandard)
End If

!-------------------------------------------------------------------
! Transplant Data
!-------------------------------------------------------------------
If ( TRANS ) Then

  Call Rcf_Get_Unit( Hdr_Aux % UnitNum )
! DEPENDS ON: file_open
  Call File_Open( Hdr_Aux % UnitNum, 'TRANSP',6 ,0 ,0 ,err )

  Call Rcf_Aux_File( Hdr_Aux, Hdr_Out, fields_out, field_count_out,&
                     transplant, IMDI, IMDI, IMDI, IMDI)

! DEPENDS ON: file_close
  Call File_Close( Hdr_Aux % UnitNum, 'TRANSP',6 ,0 ,0 ,err )
  Call Rcf_Free_Unit( Hdr_Aux % UnitNum)
  Call Rcf_FreeUMhdr( Hdr_Aux )

End If

!-------------------------------------------------------------------
! UARS Data
!-------------------------------------------------------------------
If ( UARS ) Then

  Call Rcf_Get_Unit( Hdr_Aux % UnitNum )
! DEPENDS ON: file_open
  Call File_Open( Hdr_Aux % UnitNum, 'SSU',3 ,0 ,0 ,err )

  Call Rcf_Aux_File( Hdr_Aux, Hdr_Out, fields_out, field_count_out,&
                     uars_data, IMDI, IMDI, IMDI, IMDI)

! DEPENDS ON: file_close
  Call File_Close( Hdr_Aux % UnitNum, 'SSU',3 ,0 ,0 ,err )
  Call Rcf_Free_Unit( Hdr_Aux % UnitNum)
  Call Rcf_FreeUMhdr( Hdr_Aux )

End If

End Subroutine Rcf_Create_Dump

End Module Rcf_Create_Dump_Mod
