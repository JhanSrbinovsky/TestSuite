#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  magic numbers and code for setting up field interpolation flags

Module Rcf_Set_Interp_Flags_Mod

!  Subroutine Rcf_Set_Interp_Flags - set field interpolation flags
!
! Description:
! This Module contains magic numbers and code for setting up
! field interpolation flags
!
! Method:
!  Based on the v_int_active and h_int_active logical switches.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.2   15/11/00   Allow dumps with P rather than exner. P.Selwood.
!   5.3   28/01/02   Disallow horizontal interpolation for
!                    land fraction if read from unput dump.
!                                                  Nic Gedney.
!   6.2   19/10/05   Allow vertical interpolation of ozone when levels
!                    not changing otherwise. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


!------------------------------------------------------------------
! Magic Numbers
!------------------------------------------------------------------
Integer, Parameter     :: interp_all    = 1    ! h and v to be done
Integer, Parameter     :: interp_h_only = 2    ! Do h - copy for v
Integer, Parameter     :: interp_v_only = 3    ! Do v - copy for h
Integer, Parameter     :: interp_copy   = 4    ! copy for h and v
Integer, Parameter     :: interp_no_op  = 5    ! No interp. operations
Integer, Parameter     :: interp_done   = 6    ! interp. completed
Integer, Parameter     :: interp_copied = 7    ! copy completed

Contains

!--------------------------------------------------------------------
! Subroutine to set interp flags for intput fields based on
! h_int_active, v_int_active and internal rules
!--------------------------------------------------------------------
Subroutine Rcf_Set_Interp_Flags( fields_in, fields_out, field_count_in,&
                                 field_count_out, data_source )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Data_Source_Mod, Only : &
    Data_Source_Type,           &
    Input_Dump

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active,         &
    v_int_active_soil

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_ccb,             &
    stashcode_cct,             &
    stashcode_rho,             &
    stashcode_exner,           &
    stashcode_land_frac,       &
    stashcode_p,               &
    stashcode_prog_sec

Implicit None

! Arguments
Type( field_type ), Pointer       :: fields_in(:)
Type( field_type ), Pointer       :: fields_out(:)
Type( data_source_type ), Pointer :: data_source(:)
Integer, Intent(In)               :: field_count_in
Integer, Intent(In)               :: field_count_out

! Comdecks
#include "cppxref.h"

! Local variables
Integer                           :: pos
Integer                           :: i
Logical                           :: vertical
Logical                           :: horizontal

!------------------------------------------------------------------
! Find the output fields that should be sourced from the input dump
!------------------------------------------------------------------
Do i = 1, field_count_out
  If ( data_source(i) % Source == Input_Dump ) Then
    Call Rcf_Locate( fields_out(i) % stashmaster % section,          &
                     fields_out(i) % stashmaster % item,             &
                     fields_in, field_count_in, pos )

    ! Interpolation decisions based on standard interpolation flags
    vertical   = .FALSE.
    horizontal = .FALSE.
    If (v_int_active) Then
      vertical = .TRUE.
    End If

    If (h_int_active) Then
      horizontal = .TRUE.
    End If

!-------------------------------------------------------------------
! Override standard decisions
!-------------------------------------------------------------------
    ! Only do h interpolation if single level (not ccb or cct)
    If ( fields_in(pos) % stashmaster % lv_code ==              &
                                        ppx_single_level ) Then
      If ( .NOT.                                                        &
          ( fields_in(pos) % stashmaster % section ==                   &
            stashcode_prog_sec .AND.                                    &
            ( fields_in(pos) % stashmaster % item == stashcode_ccb .OR. &
              fields_in(pos) % stashmaster % item == stashcode_cct ) ) ) Then

        vertical = .FALSE.

      End If
    End If

    ! Soil levels need special treatment
    If ( fields_in(pos) % stashmaster % lv_code == ppx_soil_level) Then
      If ( v_int_active_soil ) Then

        vertical = .TRUE.

      Else

        vertical = .FALSE.

      End If
    End If
    If ( fields_in(pos) % stashmaster % grid_type == ppx_atm_ozone .AND.&
         fields_in(pos) % levels /= fields_out(i) % levels )     Then

      ! If Ozone grid and levels differ then turn on vertical
      ! interpolation whatever the main interpolation flag
      vertical = .TRUE.
    End If


!-------------------------------------------------------------------
! Set the field interp flag
!-------------------------------------------------------------------
    If (vertical .AND. horizontal) Then
      fields_in(pos) % interp = interp_all

    Else If (vertical) Then
      fields_in(pos) % interp = interp_v_only

    Else If (horizontal) Then
      fields_in(pos) % interp = interp_h_only

    Else
      fields_in(pos) % interp = interp_copy

    End If

    ! Special treatment for rho, exner and p as need calculation rather
    ! than interpolation if would have been interpolated
    If ( fields_in(pos) % stashmaster % section == stashcode_prog_sec .AND. & 
         (fields_in(pos) % stashmaster % item == stashcode_rho        .OR.  &
          fields_in(pos) % stashmaster % item == stashcode_exner      .OR.  &
          fields_in(pos) % stashmaster % item == stashcode_p    )           &
         .AND. ( horizontal .OR. vertical )             ) Then

      fields_in(pos) % interp = interp_no_op

    End If

    If ( fields_in(pos) % stashmaster % section == stashcode_prog_sec .AND. &
         fields_in(pos) % stashmaster % item == stashcode_land_frac   .AND. &
         ( horizontal .OR. vertical )             ) Then
      If (vertical .AND. horizontal) Then
        fields_in(pos) % interp = interp_v_only
      Else If (vertical) Then
        fields_in(pos) % interp = interp_v_only
      Else
        fields_in(pos) % interp = interp_no_op
      End If
    End If


  End If
End Do

Return
End Subroutine Rcf_Set_Interp_Flags
End Module Rcf_Set_Interp_Flags_Mod
#endif
