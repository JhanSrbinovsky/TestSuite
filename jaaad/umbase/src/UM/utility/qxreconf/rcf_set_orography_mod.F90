#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up the orography for the reconfiguration

Module Rcf_Set_Orography_Mod

!  Subroutine Rcf_Set_Orography - sets the orographies
!
! Description:
! This module sets up input, output and interpolated orography
! fields for interpolation (if required ) (particularly for height
! field generation.
!
! Method:
!  Read input and output orographies - if required interpolate the
!  input orography. Do orogrphic blending for LAM and check if
!  orographic changes force vertical interpolation.
!
!  Orographic blending blends from "outside to inside" of the domain
!  with the specified weight. The blending zone *includes* the
!  Rim - thus a weight of 1 should be specified to leave the Rim
!  untouched.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   16/10/01   Orographic blending for LAM. P.Selwood
!   5.3   08/05/01   Remove some excess vertical interpolation.
!                    P.Selwood
!   5.4   08/05/02   Bug fix for rare deadlock. P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Set_Orography( fields_in, fields_out, field_count_in, &
                              field_count_out, hdr_in, hdr_out,      &
                              data_source, orog_in, orog_out,        &
                              interp_orog )

Use Rcf_Gather_Field_Mod, Only : &
    Rcf_Gather_Field_Real

Use Rcf_Scatter_Field_Mod, Only : &
    Rcf_Scatter_Field_Real

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use Rcf_Parvars_Mod, Only : &
    mype,                   &
    nproc,                  &
    datastart,              &
    current_decomp_type,    &
    gc_all_proc_group

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,               &
    Output_Grid

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_input,        &
    decomp_rcf_output

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active

Use Rcf_field_equals_mod, Only  : &
    Rcf_field_equals

Use Rcf_Set_Interp_Flags_Mod, Only :  &
    interp_copy,                      &
    interp_h_only

Use Rcf_Data_Source_Mod, Only : &
    data_source_type,           &
    Ancillary_File

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Use Rcf_Interpolate_Mod, Only : &
    Rcf_Interpolate

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field

Use Rcf_Change_Decomposition_Mod, Only : &
    Rcf_Change_Decomposition

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_orog,            &
    stashcode_prog_sec

Use Rcf_ReadNL_Horizont_Mod, Only : &
    orog_blend_width,               &
    blend_weights

Implicit None

! Arguments
Type( field_type ), Pointer          :: fields_in(:)
Type( field_type ), Pointer          :: fields_out(:)
Type( field_type ), Pointer          :: orog_in
Type( field_type ), Pointer          :: orog_out
Type( field_type ), Target           :: interp_orog
Type( data_source_type ), Intent(In) :: data_source(:)
Type( um_header_type ), Intent(In)   :: hdr_in
Type( um_header_type ), Intent(In)   :: hdr_out
Integer, Intent(In)                  :: field_count_in
Integer, Intent(In)                  :: field_count_out

! Local variables
Integer                              :: i
Integer                              :: j
Integer                              :: k
Integer                              :: pos
Integer                              :: decomp_old   ! old decomposition
Integer                              :: stat         ! gcom status
Integer                              :: v_on = 0     ! vert interp flag
Logical                              :: ancil_orog
Real                                 :: weight       ! blending weight
Real, Allocatable                    :: orog_out_fullfield   ( :, : )
Real, Allocatable                    :: orog_interp_fullfield( :, : )
Type (field_type)                    :: dummy


Nullify( dummy % data )

If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Processing Orography (stashcode 33) '
End If

!---------------------------------------------------------------
! Find origin of orography
!----------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_orog,                 &
                 fields_out, field_count_out, pos )
If ( data_source( pos ) % source == Ancillary_File ) Then
  ancil_orog = .TRUE.
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
    Write (6,*) 'Using Ancillary Orography'
  End If
Else
  ancil_orog = .FALSE.
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
    If (h_int_active) Then
      Write (6,*) 'Using interpolated Orography'
    Else
      Write (6,*) 'Copying input Orography'
    End If
  End If
End If

!------------------------------------------------------------------
! find and read input orography
!------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_orog,                 &
                 fields_in, field_count_in, pos )
orog_in  => fields_in(pos)
Call Rcf_Alloc_Field( orog_in )
Call Rcf_Read_Field( orog_in, Hdr_In, decomp_rcf_input )

!------------------------------------------------------------------
! Find output orography
!------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_orog,                 &
                 fields_out, field_count_out, pos )
orog_out => fields_out(pos)
Call Rcf_Alloc_Field( orog_out )

! Set the sizes of interp_orog to be those of orog_out
Call rcf_field_equals(interp_orog, orog_out)

!------------------------------------------------------------------
! Setup inteperpolated orography
!------------------------------------------------------------------
If (h_int_active) Then
  orog_in % interp = interp_h_only
Else
  orog_in % interp = interp_copy
End If

Call Rcf_Alloc_Field( interp_orog )
Call rcf_interpolate( orog_in, interp_orog, Input_Grid,       &
                      Output_Grid, dummy, dummy )

!-------------------------------------------------------------------
! Check for which output orography is required
!-------------------------------------------------------------------
If (ancil_orog) Then
  ! read the ancillary back in from output dump
  Call Rcf_Read_Field( orog_out, Hdr_Out, decomp_rcf_output )

  !----------------------------------------------------------------
  ! Perform the Topog masking
  !----------------------------------------------------------------
  decomp_old = decomp_rcf_output
  If (current_decomp_type /= decomp_rcf_output) Then
    decomp_old = current_decomp_type
    Call Rcf_Change_Decomposition( decomp_rcf_output )
  End If

  If ( orog_blend_width > 0 ) Then

    Allocate( orog_out_fullfield(    orog_out % glob_row_len, &
                                     orog_out % glob_rows ) )
    Allocate( orog_interp_fullfield( orog_out % glob_row_len, &
                                     orog_out % glob_rows ) )

    ! Gather orogrophies on PE 0
    ! Cannot use generic routine as fullfield is 2D so doesn't match
    ! in rank!
    Call Rcf_Gather_Field_Real( orog_out % Data(:,1),                &
                                orog_out_fullfield,                  &
                                orog_out % row_len, orog_out % rows, &
                                orog_out % glob_row_len,             &
                                orog_out % glob_rows, 0,             &
                                gc_all_proc_group )

    Call Rcf_Gather_Field_Real( interp_orog % Data(:,1),             &
                                orog_interp_fullfield,               &
                                orog_out % row_len, orog_out % rows, &
                                orog_out % glob_row_len,             &
                                orog_out % glob_rows, 0,             &
                                gc_all_proc_group )

    ! Do the orography blending on PE 0
    If (mype == 0) Then

      ! Northern and Southern Strips (including corners)
      Do i = 1, orog_out % glob_row_len
        Do j = 1, orog_blend_width

          ! First determine which weight to use
          ! Western corners
          If ( i < orog_blend_width ) Then
            weight = blend_weights( min(i,j) )

          ! Eastern corners
          Else If ( i > orog_out % glob_row_len - orog_blend_width + 1)&
                                                                   Then
            weight = blend_weights(                                 &
                               min( orog_out % glob_row_len - i + 1, j))

          ! Middle section
          Else
            weight = blend_weights( j )

          End If

          ! Set the blended field for the Southern strip
          k = j
          orog_out_fullfield(i,k) =                                   &
                    orog_interp_fullfield(i,k) * weight +             &
                    orog_out_fullfield(i,k) * (1.0 - weight )

          ! Set the blended field for the Northern strip
          k = orog_out % glob_rows - j + 1
          orog_out_fullfield(i,k) =                                   &
                    orog_interp_fullfield(i,k) * weight +             &
                    orog_out_fullfield(i,k) * (1.0 - weight )
        End Do
      End Do

      ! Western and Eastern Strips (excluding corners)
      Do i = 1, orog_blend_width
        Do j = orog_blend_width + 1, orog_out % glob_rows -           &
                                     orog_blend_width

          ! Set the weight used
          weight = blend_weights( i )

          ! Set the blended field for the Western Strip
          k = i
          orog_out_fullfield(k,j) =                                   &
                    orog_interp_fullfield(k,j) * weight +             &
                    orog_out_fullfield(k,j) * (1.0 - weight )

          ! Set the blended field for the Eastern Strip
          k = orog_out % glob_row_len - i + 1
          orog_out_fullfield(k,j) =                                   &
                    orog_interp_fullfield(k,j) * weight +             &
                    orog_out_fullfield(k,j) * (1.0 - weight )

        End Do
      End Do

    End If

    ! Only need to scatter the orog_out_fullfield
    Call Rcf_Scatter_Field_Real( orog_out % Data, orog_out_fullfield, &
                                 orog_out % row_len, orog_out % rows, &
                                 orog_out % glob_row_len,             &
                                 orog_out % glob_rows, 0,             &
                                 gc_all_proc_group )

    Call Rcf_Write_Field( orog_out, Hdr_Out, decomp_rcf_output )

    Deallocate( orog_out_fullfield )
    Deallocate( orog_interp_fullfield )

  End If

  ! Change decomposition back
  If ( current_decomp_type /= decomp_old ) Then
    Call Rcf_Change_Decomposition( decomp_old )
  End If

Else
  ! Use interpolated orography for output
  orog_out => interp_orog
  Call Rcf_Write_Field( orog_out, Hdr_Out, decomp_rcf_output )
End If

!------------------------------------------------------------------
! If there is a difference between input and output orographies,
! we need to turn on vertical interpolation
!------------------------------------------------------------------
If ( .NOT. v_int_active ) Then
  If ( orog_in % glob_level_size /= orog_out % glob_level_size ) Then
    If (ancil_orog) Then    ! Only switch on interpolation if ancillary
      v_int_active = .TRUE.
      If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,*) 'Vertical interpolation has been switched on &
                    &due to a change in orography'
      End If
    End If
  Else
    Do i = 1, orog_in % level_size
      If ( Abs( orog_in % Data( i, 1) - orog_out % Data( i, 1 ) ) > &
           Epsilon( 1.0 ) ) Then
        v_on = 1
        Exit
      End If
    End Do

    ! Need to make sure all PEs turn on v interp if 1 does
    Call GC_IMAX( 1, nproc, stat, v_on )
    If ( v_on == 1 ) Then
      v_int_active = .TRUE.
      If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
        Write (6,*) 'Vertical interpolation has been switched on &
                    &due to a change in orography'
      End If
    End If
  End If
End If

Return
End Subroutine Rcf_Set_Orography
End Module Rcf_Set_Orography_Mod
#endif
