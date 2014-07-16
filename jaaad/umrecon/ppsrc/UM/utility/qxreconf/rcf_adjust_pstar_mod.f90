
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Performs orographic adjustment of P*

Module Rcf_Adjust_Pstar_Mod

!  Subroutine Rcf_Adjust_Pstar    Adjusts P* when an ancillary
!                                 orography has been used with an
!                                 interpolated P*
!
! Description:
!   A similar method to that used in 4.5 is employed, with a
!   change to a height based vertical co-ordinate frame.
!
! Method:
!   Assuming a constant lapse rate, a surface temperature (Ts) is
!   calculated based on a temperature at a reference height (Tr)
!              Ts = Tr + lapse * (Zr - Z0_at_old_orography)
!
!   P* is then adjusted thus:-
!     P*' = P* [ Ts - lapse ( Zr - Z0_new_orography) ]  ** g/(R lapse)
!              [ ----------------------------------- ]
!              [                 Ts                  ]
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2   19/03/01   Original code.  P.Selwood
!   5.3   17/01/02   Change reference level to top of boundary layer
!                    to fix reproducability issues. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Adjust_Pstar( pstar, t, orog_in, orog_out, Input_Grid,  &
                             Output_Grid, r_theta_levels)

Use Rcf_Grid_Type_Mod, Only : &
    Grid_Type

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   &
    interp_copy,                     &
    interp_no_op

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Interpolate_Mod, Only : &
    Rcf_Interpolate

Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

Implicit None

! Arguments
Type( Grid_Type ), Intent(In)       :: Input_Grid
Type( Grid_Type ), Intent(In)       :: Output_Grid
Type( Field_Type ), Intent(InOut)   :: pstar
Type( Field_Type ), Intent(InOut)   :: orog_in
Type( Field_Type ), Intent(In)      :: orog_out
Type( Field_Type ), Intent(In)      :: t
Real, Intent(In)                    :: r_theta_levels(                 &
                                           output_grid % loc_p_field,  &
                                       0 : output_grid % model_levels )

! Local Variables/Parameters
Type( Field_Type )    :: orog_interp         ! interpolated input
                                             ! orography
Type( Field_Type )    :: dummy               ! dummy field
Real                  :: g_over_lapse_r
Real                  :: Tr                  ! T at ref. level.
Real                  :: Zr                  ! Height at ref. level.
Real                  :: Z0i                 ! Interpolated orog.
Real                  :: Z0o                 ! Output orog.
Integer               :: ref_lev             ! first theta level above
                                             ! ref_height
Integer               :: i                   ! Looper

! Comdecks
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_A----------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Replace variable A by more meaningful name for
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Convert to Fixed/Free format. P. Selwood

      ! Mean radius of Earth in m.
      Real, Parameter  :: Earth_Radius = 6371000.

!*----------------------------------------------------------------------
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
!*L------------------COMDECK C_LAPSE ----------------------------------
      Real, Parameter :: lapse      = 0.0065  ! Near surface lapse rate
      Real, Parameter :: lapse_trop = 0.002   ! Tropopause lapse rate
!*---------------------------------------------------------------------


Nullify( dummy % Data )

!--------------------------------------------------------------------
! Set up the orography fields to obtain interpolated orography
!--------------------------------------------------------------------
Call Rcf_Field_Equals( orog_interp, orog_out )

If (h_int_active) Then
  orog_in % interp   = interp_h_only
Else
  orog_in % interp   = interp_copy
End If

orog_interp % interp = interp_no_op

Call Rcf_Alloc_Field( orog_interp )
Call Rcf_Interpolate( orog_in, orog_interp, Input_Grid, Output_Grid, &
                      dummy, dummy )

! Reference level is last boundary layer level
ref_lev = Output_Grid % BL_levels

!--------------------------------------------------------------------
! Can now do the p* adjustment calculation as specified above
!--------------------------------------------------------------------
g_over_lapse_r = g/(lapse * R)
Do i = 1, pstar % level_size
  Tr  = t % Data( i, ref_lev )
  Zr  = r_theta_levels( i, ref_lev ) - Earth_Radius
  Z0i = orog_interp % Data( i, 1 )     ! interpolated orography
  Z0o = orog_out % Data( i, 1 )        ! output orography

  pstar % Data( i, 1 ) = pstar % Data( i, 1 ) *                       &
                      ( ( Tr + lapse * (Zr - Z0o) ) /                 &
                        ( Tr + lapse * (Zr - Z0i) ) ) ** g_over_lapse_r
End Do


!--------------------------------------------------------------------
! Tidy up
!--------------------------------------------------------------------
Call Rcf_Dealloc_Field( orog_interp )


Return
End Subroutine Rcf_Adjust_Pstar

End Module Rcf_Adjust_Pstar_Mod
