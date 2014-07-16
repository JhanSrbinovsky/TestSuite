#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Deriving 3d convective cloud amounts.

Module Rcf_Derv_3D_CCA_Mod

!  Subroutine Rcf_Derv_3D_CCA
!
! Description:
!    Top level routine obtaining variables (and some conversions)
!    for the Rcf_Calc_3D_CCA routine.
!
! Method:
!    Derived from UM4.5 code.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   25/10/01   Cater for extra rho level at top. D.Robinson
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Derv_3D_CCA( fields_in, fields_out, field_count_in, &
                            field_count_out, hdr_in, hdr_out, cca_3d)

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
    Interp_h_only

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

Use Rcf_Generate_Heights_Mod, Only : &
    Rcf_Generate_Heights

Use Rcf_Interpolate_Mod, Only : &
    Rcf_Interpolate

Use Rcf_Exner_P_Convs_Mod, Only : &
    Rcf_Conv_Exner_P

Use Rcf_Theta_T_Convs_Mod, Only : &
    Rcf_Conv_Theta_T

Use Rcf_Calc_3d_CCA_Mod, Only : &
    Rcf_Calc_3d_CCA

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_theta,              stashcode_cca, &
    stashcode_ccb,                stashcode_cct, &
    stashcode_orog,               stashcode_exner,&
    stashcode_prog_sec

Implicit None

! Arguments
Type( field_type ), Pointer        :: fields_in(:)
Type( field_type ), Pointer        :: fields_out(:)
Type( field_type ), Intent(InOut)  :: cca_3d
Type( um_header_type ), Intent(In) :: hdr_in
Type( um_header_type ), Intent(In) :: hdr_out
Integer, Intent(In)                :: field_count_in
Integer, Intent(In)                :: field_count_out

! Comdecks
#include "cppxref.h"

! Local variables
Integer                       :: i           ! looper
Integer                       :: pos         ! field position

Real                          :: theta_heights(                     &
                                    output_grid % loc_p_field,      &
                                    0 : output_grid % model_levels + 1)
                                    ! Height of theta levels
Real                          :: rho_heights(                       &
                                    output_grid % loc_p_field,      &
                                    0 : output_grid % model_levels + 1)
                                    ! Height of rho levels

Type( field_type ), Pointer   :: cca_in_2d
Type( field_type ), Pointer   :: theta
Type( field_type ), Pointer   :: exner
Type( field_type ), Pointer   :: ccb
Type( field_type ), Pointer   :: cct
Type( field_type ), Pointer   :: orog
Type( field_type )            :: cca_out_2d
Type( field_type )            :: dummy

! Nullify field data
Nullify(dummy % Data)
Nullify(dummy % Data_Int)
Nullify(dummy % Data_Log)
Nullify(cca_out_2d % Data )
Nullify(cca_out_2d % Data_Int )
Nullify(cca_out_2d % Data_Log )

!--------------------------------------------------------------
! Write out out action if appropriate
!--------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Initialising 3D CCA'
End If

!--------------------------------------------------------------
! Find and setup 2D CCA from input
!--------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_cca,                  &
                 fields_in, field_count_in, pos )
cca_in_2d => fields_in(pos)
cca_in_2d % interp = interp_h_only
Call Rcf_Alloc_Field( cca_in_2d )
Call Rcf_Read_Field( cca_in_2d, hdr_in, decomp_rcf_input )

!--------------------------------------------------------------
! cca_out_2d must be set up - uses a combination of data from
! cca_in_2d and cca_3d, then interpolate
!--------------------------------------------------------------
Call rcf_field_equals( cca_out_2d, cca_in_2d )
cca_out_2d % rows            = cca_3d % rows
cca_out_2d % row_len         = cca_3d % row_len
cca_out_2d % level_size      = cca_3d % level_size
cca_out_2d % glob_rows       = cca_3d % glob_rows
cca_out_2d % glob_row_len    = cca_3d % glob_row_len
cca_out_2d % glob_level_size = cca_3d % glob_level_size

Call Rcf_Alloc_Field( cca_out_2d )

Call rcf_interpolate( cca_in_2d, cca_out_2d, input_grid, output_grid, &
                      dummy, dummy )

Call Rcf_Dealloc_Field( cca_in_2d )

!---------------------------------------------------------------
! Read conv cloud base and top as these are required
!---------------------------------------------------------------
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

!----------------------------------------------------------------
! Interpolation can mean that ccb, cct and cca_out_2d are not
! consistent. A quick `clean up' should prevent crashes.
!----------------------------------------------------------------
Do i = 1, cca_out_2d % level_size
  ! Make sure cca is +ive
  If (cca_out_2d % Data(i,1) < 0. ) Then
    cca_out_2d % Data(i,1) = 0
  End If

  ! Make sure cca only exists where ccb and cct are both non-zero
  If ( ccb % Data_Int(i,1) == 0 .AND.     &
       cct % Data_Int(i,1) == 0 .AND.     &
       cca_out_2d % Data(i,1) > 0. ) Then
    cca_out_2d % Data(i,1) = 0.0
  End If
End Do

!----------------------------------------------------------------
! Need pressure and temperature on theta levels - read exner and
! theta and do the necessary calculations
!----------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_theta,                &
                 fields_out, field_count_out, pos )
theta => fields_out(pos)
Call Rcf_Alloc_Field( theta )
Call Rcf_Read_Field( theta, hdr_out, decomp_rcf_output)

! Need orography for heights
Call Rcf_Locate( stashcode_prog_sec, stashcode_orog,                 &
                 fields_out, field_count_out, pos )
orog => fields_out(pos)
Call Rcf_Alloc_Field( orog )
Call Rcf_Read_Field( orog, hdr_out, decomp_rcf_output )

! Generate heights
Call rcf_generate_heights( output_grid,orog % Data(1:,1:),ppx_atm_tall,&
                           ppx_theta_level, theta_heights,           &
                           theta % level_size )

Call rcf_generate_heights( output_grid,orog % Data(1:,1:),ppx_atm_tall,&
                           ppx_rho_level, rho_heights,               &
                           theta % level_size )


! Theta to T
Call Rcf_Conv_Theta_T( theta, fields_out, field_count_out, hdr_out, &
                       decomp_rcf_output, rho_heights, theta_heights )


! Read exner
Call Rcf_Locate( stashcode_prog_sec, stashcode_exner,                &
                 fields_out, field_count_out, pos )
exner => fields_out(pos)
Call Rcf_Alloc_Field( exner )
Call Rcf_Read_Field( exner, hdr_out, decomp_rcf_output)

! exner to P
Call Rcf_Conv_Exner_P( exner )

!-----------------------------------------------------------------
! Have all fields we require - call the calculation routine
!-----------------------------------------------------------------
Call Rcf_Calc_3D_CCA( cca_out_2d, ccb, cct, theta, exner, cca_3d )

!-----------------------------------------------------------------
! Write out 3D CCA and clean up
!-----------------------------------------------------------------
Call Rcf_Write_Field( cca_3d, hdr_out, decomp_rcf_output )

Call Rcf_Dealloc_Field( theta )
Call Rcf_Dealloc_Field( exner )
Call Rcf_Dealloc_Field( cca_out_2d )
Call Rcf_Dealloc_Field( ccb )
Call Rcf_Dealloc_Field( cct )
Call Rcf_Dealloc_Field( orog )

Return
End Subroutine Rcf_Derv_3D_CCA
End Module Rcf_Derv_3D_CCA_Mod
#endif
