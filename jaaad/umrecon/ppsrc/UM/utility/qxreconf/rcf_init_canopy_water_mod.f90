
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialisations for tiled canopy water field

Module Rcf_Init_Canopy_Water_Mod

! Subroutine Rcf_Init_Canopy_Water
!
! Description:
!   Initialises canopy water on tiles.
!
! Method:
!   Initialises canopy water on tiles from mean canopy water.  This
!   latter field is only found in the input dump, so needs to be
!   interpolated.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2   13/10/00   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Init_Canopy_Water( fields_in,  field_count_in,  hdr_in, &
                                  canopy_water)

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_mean_canopyw,    &
    stashcode_prog_sec

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Interpolate_Mod, Only : &
    Rcf_Interpolate

Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_input

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,               &
    Output_Grid

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   &
    interp_copy

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Implicit None

! Arguments
Type( field_type ), Pointer          :: fields_in(:)
Type( field_type ), Intent( InOut )  :: canopy_water
Type( um_header_type ), Intent( In ) :: hdr_in

Integer, Intent( In )                :: field_count_in

! Local variables
Type (field_type), Pointer           :: mean_canopyw_in
Type (field_type)                    :: mean_canopyw_out
Type (field_type)                    :: dummy  ! pretend orography
Integer                              :: pos    ! field position
Integer                              :: i      ! looper


! Nullify field data
Nullify(dummy % Data)
Nullify(dummy % Data_Int)
Nullify(dummy % Data_Log)
Nullify(mean_canopyw_out % Data)
Nullify(mean_canopyw_out % Data_Int)
Nullify(mean_canopyw_out % Data_Log)
!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Initialising tile canopy water from input mean &
              &canopy water'
End If

!----------------------------------------------------------------------
! Find mean canopy water in input fields and read it in
!----------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_mean_canopyw,         &
                 fields_in, field_count_in, pos)
mean_canopyw_in => fields_in(pos)
Call Rcf_Alloc_Field( mean_canopyw_in )
Call Rcf_Read_Field( mean_canopyw_in, hdr_in, decomp_rcf_input )

!----------------------------------------------------------------------
! Initialise mean_canopyw_out - and allocate space
!----------------------------------------------------------------------
Call Rcf_Field_Equals( mean_canopyw_out, mean_canopyw_in )
mean_canopyw_out % rows            = canopy_water % rows
mean_canopyw_out % row_len         = canopy_water % row_len
mean_canopyw_out % level_size      = canopy_water % level_size
mean_canopyw_out % glob_rows       = canopy_water % glob_rows
mean_canopyw_out % glob_row_len    = canopy_water % glob_row_len
mean_canopyw_out % glob_level_size = canopy_water % glob_level_size

If (h_int_active) Then
  mean_canopyw_in % interp = interp_h_only
Else
  mean_canopyw_in % interp = interp_copy
End If

Call Rcf_Alloc_Field( mean_canopyw_out )

!----------------------------------------------------------------------
! Interpolate mean_canopyw and copy result to canopy_water levels 1-6
! (levels 1-5 are vegetation, level 6 is urban).  On other levels
! (water, bare soil, ice) set canopy_water to 0.
!----------------------------------------------------------------------
Call Rcf_Interpolate( mean_canopyw_in, mean_canopyw_out, input_grid, &
                      output_grid, dummy, dummy )

!jhan: optio for n tles quai-hard-wired reflecting input dataset
If(canopy_water % levels == 1) then
  canopy_water % Data(:,1) = mean_canopyw_out % Data(:,1)
Else if(canopy_water % levels == 9) then
  Do i = 1, 6
    canopy_water % Data(:,i) = mean_canopyw_out % Data(:,1)
  End Do
  Do i = 7, 9
    canopy_water % Data(:,i) = 0.0
  End Do
!jhan{appropriate for  17 tiles testing
Else if(canopy_water % levels == 17) then
  Do i = 1, 14
    canopy_water % Data(:,i) = mean_canopyw_out % Data(:,1)
  End Do
  Do i = 15, 17
       canopy_water % Data(:,i) = 0.0
  End Do
!}jhan
Endif

!----------------------------------------------------------------------
! Tidy up
!----------------------------------------------------------------------
Call Rcf_Dealloc_Field( mean_canopyw_out )
Call Rcf_Dealloc_Field( mean_canopyw_in )

End Subroutine Rcf_Init_Canopy_Water

End Module Rcf_Init_Canopy_Water_Mod
