#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Horizontally interpolates T and Pstar for Height generation

Module Rcf_GRIB_Interp_TnPstar_Mod

!  Subroutine Rcf_GRIB_Interp_TnPstar
!
! Description: A routine to allocate space and call rcf_interpolate
!              to provide horizontally interpolated versions of
!              input T and Pstar on the output grid's horizontal
!              resolution. (copied and modified rcf_set_orography)
!
! Method: Allocate appropriate space
!         Set interp flag to horizontal only
!         Call rcf_interpolate
!
! Current Code Owner: R.Sharp
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.5   01/09/02   Original code.  R.Sharp
!   6.2   20/06/05   Made alterations for use of model level
!                    ECMWF GRIB data. Paul Earnshaw (frpe)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_GRIB_Interp_TnPstar( fields_in, fields_out, Hdr_In,    &
                                    field_count_in, field_count_out)

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,               &
    Output_Grid

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_input

Use Rcf_Set_Interp_Flags_Mod, Only :  &
    interp_copy,                      &
    interp_h_only

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_orog,          &
    stashcode_pstar,         &
    stashcode_theta,         &
    stashcode_prog_sec

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_field_equals_mod, Only  : &
    Rcf_field_equals

Use Rcf_Interpolate_Mod, Only : &
    Rcf_Interpolate

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_GRIB_T_n_Pstar_H_Interp_Mod, Only : &
    GRIB_T,            &
    GRIB_Pstar,        &
    GRIB_Levels

Implicit None

Type( field_type ), Pointer          :: Pstar_in

! Arguments
Type( field_type ), Pointer          :: fields_in(:)
Type( field_type ), Pointer          :: fields_out(:)
Type( um_header_type ), Intent(In)   :: hdr_in
Integer, Intent(In)                  :: field_count_in
Integer, Intent(In)                  :: field_count_out

! Local variables

Type( field_type ), Pointer          :: orog_out
Type( field_type ), Pointer          :: T_in
Type( field_type ), Pointer          :: T_out
Integer                              :: pos
Type (field_type)                    :: dummy

Nullify( dummy % data )

If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Horizontally Interpolating T and Pstar'
End If

!------------------------------------------------------------------
! find and read input Pstar
!------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_pstar,                &
                 fields_in, field_count_in, pos )
Pstar_in  => fields_in(pos)

Call Rcf_Alloc_Field( Pstar_in )
Call Rcf_Read_Field( Pstar_in, Hdr_In, decomp_rcf_input )

!------------------------------------------------------------------
! Setup Pstar field with output horizontal resolution
!------------------------------------------------------------------
! Pstar is (probably) not in the output dump, it does however need
! to be interpolated. Use the output orography for the 'donor' level
! descriptor.
Call Rcf_Locate( stashcode_prog_sec, stashcode_orog,                 &
                 fields_out, field_count_out, pos )
orog_out => fields_out(pos)
Call Rcf_Alloc_Field( orog_out )

! Set the sizes of GRIB_Pstar to be those of orog_out
Call rcf_field_equals(GRIB_Pstar, orog_out)

! Now need to reset some values which aren't correct
GRIB_Pstar % dump_pos        = Pstar_in % dump_pos
GRIB_Pstar % interp          = interp_h_only
GRIB_Pstar % stashmaster    => Pstar_in % stashmaster

!------------------------------------------------------------------
! Interpolate Pstar
!------------------------------------------------------------------
If (h_int_active) Then
  Pstar_in % interp = interp_h_only
Else
  Pstar_in % interp = interp_copy
End If


Call Rcf_Alloc_Field( GRIB_Pstar )
Call rcf_interpolate( Pstar_in, GRIB_Pstar, Input_Grid,       &
                      Output_Grid, dummy, dummy )

!------------------------------------------------------------------
! find and read input T
!------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_theta,                &
                 fields_in, field_count_in, pos )
T_in  => fields_in(pos)
Call Rcf_Alloc_Field( T_in )
Call Rcf_Read_Field( T_in, Hdr_In, decomp_rcf_input )

!------------------------------------------------------------------
! Setup T field with output horizontal resolution
!------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_theta,                &
                 fields_out, field_count_out, pos )
T_out => fields_out(pos)
Call Rcf_Alloc_Field( T_out )

! Set the sizes of GRIB_T to be those of T_out
Call rcf_field_equals(GRIB_T, T_out)

! Now need to reset some values to maintain same vertical levels
GRIB_T % levels          = T_in % levels

!------------------------------------------------------------------
! Interpolate T
!------------------------------------------------------------------
If (h_int_active) Then
  T_in % interp = interp_h_only
Else
  T_in % interp = interp_copy
End If

Call Rcf_Alloc_Field( GRIB_T )
Call rcf_interpolate( T_in, GRIB_T, Input_Grid,       &
                      Output_Grid, dummy, dummy )

!------------------------------------------------------------------
! Save level definitons
!------------------------------------------------------------------
Allocate ( GRIB_Levels(T_in % levels) )

GRIB_Levels(:) = Hdr_In % LevDepC (1:T_in % levels,1)

!------------------------------------------------------------------
! Clean up
!------------------------------------------------------------------
Call Rcf_DeAlloc_Field( Pstar_in )
Call Rcf_DeAlloc_Field( T_in )
Call Rcf_DeAlloc_Field( T_out )

Return
End Subroutine Rcf_GRIB_Interp_TnPstar
End Module Rcf_GRIB_Interp_TnPstar_Mod
#endif
