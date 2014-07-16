#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculations between T and Theta

Module Rcf_theta_t_Convs_mod

!  Subroutine Rcf_Conv_Theta_T - converts theta to t
!  Subroutine Rcf_Conv_T_Theta - converts t to theta
!
! Description:
!   Performs calculations between T and Theta
!
! Method:
!   Derived from New Dynamics 2.7 code.
!
! Current Code Owner:
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.2   15/11/00   Allow use of P rather than Exner. P.Selwood.
!   5.3   27/07/01   Linearly interpolate exner directly, without
!                    using the intermediary pressure variable. S.Cusack
!   5.3   25/10/01   Cater for extra rho level at top. D. Robinson.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

!******************************************************************
! Routine to calculate T from Theta
!******************************************************************
Subroutine Rcf_Conv_theta_t( theta, fields, field_count, hdr, decomp, &
                             rho_heights, theta_heights)

Use Rcf_field_equals_mod, Only : &
    Rcf_Field_Equals

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Rcf_Submodel_Mod, Only : &
    Atmos_IM

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Calc_Exner_Theta_Mod, Only : &
    Rcf_Calc_Exner_Theta

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Rcf_Exner_P_Convs_Mod, Only : &
    Rcf_Conv_Exner_P,             &
    Rcf_Conv_P_Exner

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_exner,           &
    stashcode_p,               &
    stashcode_prog_sec

Implicit None

! Arguments
Type( field_type ), Pointer        :: fields(:)  ! array of fields
Type( field_type ), Intent(InOut)  :: theta      ! theta field - convert
Type( um_header_type ), Intent(In) :: hdr        ! header of dump
Integer, Intent(In)                :: decomp     ! working decomposition
Integer, Intent(In)                :: field_count ! # of fields in dump
Real, Intent(In)                   :: theta_heights(theta % level_size,&
                                                0 : theta % levels + 1)
Real, Intent(In)                   :: rho_heights( theta % level_size, &
                                                0 : theta % levels + 1)
! Comdecks
#include "cppxref.h"

! Local variables
Integer                            :: i
Integer                            :: j
Integer                            :: pos
Type( field_type )                 :: exner_theta
Type( field_type )                 :: exner

! Nullify data
Nullify( exner % Data )
Nullify( exner % Data_Int )
Nullify( exner % Data_Log )
Nullify( exner_theta % Data )
Nullify( exner_theta % Data_Int )
Nullify( exner_theta % Data_Log )

!------------------------------------------------------------------
! Get exner and convert it to P on rho levels
! (for VAR we will only have P - just read it in)
!------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_exner,                &
                 fields, field_count, pos, .TRUE. )

If (pos /= 0) Then
  Call Rcf_Field_Equals( exner, fields(pos) )
  Call Rcf_Alloc_Field( exner )
  Call Rcf_Read_Field( exner, hdr, decomp )

Else
  Call Rcf_Locate( stashcode_prog_sec, stashcode_p,                  &
                   fields, field_count, pos )
  Call Rcf_Field_Equals( exner, fields(pos) )
  Call Rcf_Alloc_Field( exner )
  Call Rcf_Read_Field( exner, hdr, decomp )

  Call Rcf_Conv_P_Exner( exner )
End If

!-------------------------------------------------------------------
! Convert exner from rho levels to theta levels
!-------------------------------------------------------------------
Call Rcf_Field_Equals( exner_theta, exner )
exner_theta % levels = theta % levels
Call Rcf_Alloc_Field( exner_theta )

Call Rcf_calc_exner_theta( exner % level_size, exner_theta % levels, &
                           theta_heights, rho_heights(1:,1:),        &
                           exner % Data, exner_theta % Data )

!-------------------------------------------------------------------
! Calculate T from theta and exner
!-------------------------------------------------------------------
Do j = 1, theta % levels
  Do i = 1, theta % level_size
    theta % Data(i,j) = theta % Data(i,j) * exner_theta % Data(i,j)
  End Do
End Do

!-------------------------------------------------------------------
! Tidy up
!-------------------------------------------------------------------
Call Rcf_Dealloc_Field( exner_theta )
Call Rcf_Dealloc_Field( exner )


Return
End Subroutine Rcf_Conv_theta_t





!******************************************************************
! Routine to calculate theta from T
!******************************************************************
Subroutine Rcf_Conv_t_theta( theta, fields, field_count, hdr, decomp, &
                             rho_heights, theta_heights )

Use Rcf_field_equals_mod, Only : &
    Rcf_Field_Equals

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Rcf_Submodel_Mod, Only : &
    Atmos_IM

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Calc_Exner_Theta_Mod, Only : &
    Rcf_Calc_exner_Theta

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Rcf_Exner_P_Convs_Mod, Only : &
    Rcf_Conv_Exner_P,             &
    Rcf_Conv_P_Exner

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_exner,           &
    stashcode_p,               &
    stashcode_prog_sec

Implicit None

! Arguments
Type( field_type ), Pointer        :: fields(:)  ! array of fields
Type( field_type ), Intent(InOut)  :: theta      ! theta field - convert
Type( um_header_type ), Intent(In) :: hdr        ! header of dump
Integer, Intent(In)                :: decomp     ! working decomposition
Integer, Intent(In)                :: field_count ! # of fields in dump
Real, Intent(In)                   :: theta_heights(theta % level_size,&
                                                   0 : theta % levels )
Real, Intent(In)                   :: rho_heights( theta % level_size, &
                                                0 : theta % levels + 1)

! Comdecks
#include "cppxref.h"

! Local variables
Integer                            :: i
Integer                            :: j
Integer                            :: pos
Type( field_type )                 :: exner_theta
Type( field_type )                 :: exner

! Nullify data of new fields
Nullify( exner % Data )
Nullify( exner % Data_Int )
Nullify( exner % Data_Log )
Nullify( exner_theta % Data )
Nullify( exner_theta % Data_Int )
Nullify( exner_theta % Data_Log )

!------------------------------------------------------------------
! Get exner and convert it to P on rho levels
! (for VAR we will only have P - just read it in)
!------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_exner, &
                 fields, field_count, pos, .TRUE. )

If (pos /= 0) Then
  Call Rcf_Field_Equals( exner, fields(pos) )
  Call Rcf_Alloc_Field( exner )
  Call Rcf_Read_Field( exner, hdr, decomp )

Else
  Call Rcf_Locate( stashcode_prog_sec, stashcode_p, &
                   fields, field_count, pos )
  Call Rcf_Field_Equals( exner, fields(pos) )
  Call Rcf_Alloc_Field( exner )
  Call Rcf_Read_Field( exner, hdr, decomp )

  Call Rcf_Conv_P_Exner( exner )
End If

!-------------------------------------------------------------------
! Convert exner from rho levels to theta levels
!-------------------------------------------------------------------
Call Rcf_Field_Equals( exner_theta, exner )
exner_theta % levels = theta % levels
Call Rcf_Alloc_Field( exner_theta )

Call Rcf_calc_exner_theta( exner % level_size, exner_theta % levels, &
                           theta_heights, rho_heights(1:,1:),        &
                           exner % Data, exner_theta % Data )

!-------------------------------------------------------------------
! Calculate theta from T and exner
!-------------------------------------------------------------------
Do j = 1, theta % levels
  Do i = 1, theta % level_size
    theta % Data(i,j) = theta % Data(i,j) / exner_theta % Data(i,j)
  End Do
End Do

!-------------------------------------------------------------------
! Tidy up
!-------------------------------------------------------------------
Call Rcf_Dealloc_Field( exner_theta )
Call Rcf_Dealloc_Field( exner )

Return
End Subroutine Rcf_Conv_t_theta

End Module Rcf_theta_t_Convs_mod

#endif
