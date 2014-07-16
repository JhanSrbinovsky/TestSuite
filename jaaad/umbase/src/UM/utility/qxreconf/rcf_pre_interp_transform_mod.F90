#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Perform transformations etc before a field is interpolated

Module Rcf_Pre_Interp_Transform_Mod

!  Subroutine Rcf_Pre_Interp_Transform
!
! Description:
!   Wrapper to perform pre-interpolation transforms to those fields
!   that require it.
!
! Method:
!   Choice of transform is based on stashcode.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Pre_Interp_Transform( input_field, fields_in, &
                                     field_count_in, hdr_in, orog_in )

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Exner_P_Convs_Mod, Only : &
    Rcf_Conv_Exner_P

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_input

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid

Use Rcf_theta_t_convs_mod, Only : &
    Rcf_conv_theta_t

Use Rcf_Generate_Heights_Mod, Only : &
    Rcf_Generate_Heights

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   &
    interp_v_only,                   &
    interp_all

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_theta,           &
    stashcode_exner,           &
    stashcode_prog_sec,        &
    stashcode_soil_moist

Use Rcf_Recon_Mod, Only : &
    GRIB,                 &
    use_smc_stress

Use Rcf_smc_stress_Mod,only: &
    Rcf_smc_stress

Implicit None

! Arguments
Type( field_type ), Intent( InOut )   :: input_field
Type( field_type ), Pointer           :: fields_in(:)
Type( field_type ), Intent(In)        :: orog_in
Type( um_header_type ), Intent(In)    :: hdr_in
Integer,            Intent(In)        :: field_count_in

! Comdecks
#include "cppxref.h"

! Local variables
Logical                       :: interpolate
Real                          :: theta_heights(                     &
                                    input_grid % loc_p_field,       &
                                    0 : input_grid % model_levels + 1)
                                    ! Height of theta levels
Real                          :: rho_heights(                       &
                                    input_grid % loc_p_field,       &
                                    0 : input_grid % model_levels + 1)
                                    ! Height of rho levels

interpolate = ( input_field % interp == interp_h_only .OR. &
                input_field % interp == interp_v_only .OR. &
                input_field % interp == interp_all )
!---------------------------------------------------------------
! Only do transforms if interpolation is switched on
!---------------------------------------------------------------
If ( interpolate ) Then

  Select Case( input_field % stashmaster % section )

  Case (stashcode_prog_sec)
    ! Which fields do we wish to apply transforms to?
    Select Case( input_field % stashmaster % item )

    Case( stashcode_exner )
      ! convert exner to P
      If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
        Write (6,'(a25)') 'Converting exner to P'
      End If
      Call Rcf_Conv_Exner_P( input_field )

    Case( stashcode_theta )
      ! convert theta to T
      ! Unless input data was in GRIB format when it is assumed that
      ! the field stored as Theta was actually read in as T already.
      If ( .NOT.GRIB ) Then
        If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
          Write (6,'(a25)') 'Converting Theta to T'
        End If

        Call rcf_generate_heights( input_grid, orog_in % Data(1:,1:), &
                                   ppx_atm_tall,ppx_theta_level,      &
                                   theta_heights,                     &
                                   input_field % level_size)

        Call rcf_generate_heights( input_grid, orog_in % Data(1:,1:), &
                                   ppx_atm_tall, ppx_rho_level,       &
                                   rho_heights,                       &
                                   input_field % level_size)

        Call Rcf_Conv_Theta_T( input_field, fields_in,field_count_in, &
                               hdr_in, decomp_rcf_input, rho_heights, &
                               theta_heights )
      Else
        If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
          Write (6,'(a41)') 'Assuming Theta is T from GRIB data source'
        End If

      End If
    Case( stashcode_soil_moist )
       If ( use_smc_stress ) Then
           Call Rcf_smc_stress( input_field, fields_in,               &
                                field_count_in, hdr_in )
       End if

    End Select
  End Select

End If

Return
End Subroutine Rcf_Pre_Interp_Transform
End Module Rcf_Pre_Interp_Transform_Mod
#endif
