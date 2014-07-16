
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Rcf_Calc_Output_Exner_Mod - calculate output dump exner pressure

Module Rcf_Calc_Output_Exner_Mod

!  Subroutine Rcf_Calc_Output_Exner
!
! Description:
!    This module calculates exner pressure for the output dump
!
! Method:
!    This code is derived from that used in the New Dynamics
!    reconfiguration program. First P* is calculated on the input
!    grid and then interpolated horizontally onto the output grid.
!    This is then used to calculate the 1st level of exner which
!    in turn is used to calculate the higher levels.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   18/05/00   Original code.  P.Selwood
!   5.2   20/03/01   Adjust P* for ancillary orography. P.Selwood.
!   5.3   25/10/01   Cater for extra rho level at top. D.Robinson.
!   5.3   06/12/01   Correct inconsistency in T_v calculation A.Malcolm
!   5.5   21/01/03   Allow for read in of pstar if present. R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Calc_Output_Exner( fields_in, field_count_in, orog_in, &
                                  hdr_in, orog_out, T, Q, Exner,      &
                                  orog_source, r_rho_levels,          &
                                  r_theta_levels )
Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only : &
    UM_header_type

Use Rcf_Grid_Type_Mod, Only: &
    Input_Grid,          &
    Output_Grid

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_pstar,           &
    stashcode_prog_sec

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   &
    interp_copy,                     &
    interp_no_op

Use Rcf_Submodel_Mod, Only : &
    Atmos_IM

Use Rcf_Calc_P_Star_Mod, Only : &
    Rcf_Calc_P_Star

Use Rcf_DecompTP_Mod, Only : &
    Decomp_rcf_input

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Interpolate_Mod, Only : &
    Rcf_Interpolate

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Adjust_Pstar_Mod, Only : &
    Rcf_Adjust_Pstar

Use Rcf_Data_Source_Mod, Only : &
    Ancillary_File

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Implicit None

! Arguments
Integer, Intent(In)                :: field_count_in
Type( field_type ), Pointer        :: fields_in(:)
Type( field_type ), Intent(InOut)  :: orog_in
Type( field_type ), Intent(In)     :: orog_out
Type( field_type ), Intent(In)     :: T
Type( field_type ), Intent(In)     :: Q
Type( field_type ), Intent(InOut)  :: Exner
Type( UM_header_type ), Intent(In) :: hdr_in
Real, Intent(In)                   :: r_rho_levels(               &
                                      output_grid % loc_p_field,  &
                                      0 : output_grid % model_levels+1)
Real, Intent(In)                   :: r_theta_levels(             &
                                      output_grid % loc_p_field,  &
                                      0 : output_grid % model_levels+1)
Integer, Intent(In)                :: orog_source

! Local variables
Integer                            :: i          ! Looper
Integer                            :: k          ! Looper
Integer                            :: pos        ! position in fields_in
Real                               :: temp       ! temporary value
Real                               :: weight1    ! averaging weight
Real                               :: weight2    ! averaging weight
Real                               :: deltaz     ! vertical spacing
Type( field_type )                 :: pstar_in
Type( field_type )                 :: pstar_out
Type( field_type )                 :: dummy

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
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

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
!*L------------------COMDECK C_EPSLON-----------------------------------
! EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR

      Real, Parameter :: Epsilon   = 0.62198
      Real, Parameter :: C_Virtual = 1./Epsilon-1.

!*----------------------------------------------------------------------

Nullify( dummy % Data )

If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
  Write (6,*) 'Calculating Exner'
End If

!--------------------------------------------------------------------
! Setup p_star_in field - either read or calculate from input dump
!--------------------------------------------------------------------
! Check input dump for p_star
Call Rcf_Locate( stashcode_prog_sec, stashcode_Pstar,                &
                 fields_in, field_count_in, pos, .True. )

! If present, read in, else calculate.
If (pos /= 0) Then
  Pstar_in  = fields_in(pos)
  Call Rcf_Alloc_Field( Pstar_in )
  Call Rcf_Read_Field( Pstar_in, Hdr_In, decomp_rcf_input )

Else

  Nullify( pstar_in % Data )
  Nullify( pstar_in % Data_Int )
  Nullify( pstar_in % Data_Log )
  pstar_in % levels          = 1
  pstar_in % rows            = Input_Grid % loc_p_rows
  pstar_in % row_len         = Input_Grid % loc_p_row_length
  pstar_in % level_size      = Input_Grid % loc_p_field
  pstar_in % glob_rows       = Input_Grid % glob_p_rows
  pstar_in % glob_row_len    = Input_Grid % glob_p_row_length
  pstar_in % glob_level_size = Input_Grid % glob_p_field
  pstar_in % dump_pos        = IMDI

  pstar_in % stashmaster     => Rcf_Exppx(Atmos_IM, 0, stashcode_pstar)

  Call Rcf_Alloc_Field( pstar_in )
  Call Rcf_Calc_P_Star( Input_Grid, fields_in, field_count_in, hdr_in, &
                        decomp_rcf_input, orog_in, pstar_in )

End If

! Only need to do horizontal interpolation of P* if h_int_active
If (h_int_active) Then
  pstar_in % interp          = interp_h_only
Else
  pstar_in % interp          = interp_copy
End If

!--------------------------------------------------------------------
! Setup p_star_out field - interpolated from p_star_in
!--------------------------------------------------------------------
Nullify( pstar_out % Data )
Nullify( pstar_out % Data_Int )
Nullify( pstar_out % Data_Log )
pstar_out % levels          = 1
pstar_out % rows            = Output_Grid % loc_p_rows
pstar_out % row_len         = Output_Grid % loc_p_row_length
pstar_out % level_size      = Output_Grid % loc_p_field
pstar_out % glob_rows       = Output_Grid % glob_p_rows
pstar_out % glob_row_len    = Output_Grid % glob_p_row_length
pstar_out % glob_level_size = Output_Grid % glob_p_field
pstar_out % dump_pos        = IMDI
pstar_out % interp          = interp_no_op
pstar_out % stashmaster     => Rcf_Exppx( Atmos_IM, 0, stashcode_pstar )

Call Rcf_Alloc_Field( pstar_out )
Call Rcf_Interpolate( pstar_in, pstar_out, Input_Grid, Output_Grid, &
                      dummy, dummy)

!--------------------------------------------------------------------
! If required, adjust the p* field for orographic effects
!--------------------------------------------------------------------
If ( orog_source == Ancillary_File ) Then
  Call Rcf_Adjust_Pstar( pstar_out, t, orog_in, orog_out, Input_Grid, &
                         Output_Grid, r_theta_levels )
End If

!--------------------------------------------------------------------
! calculate exner
! equation is
!
! exner(k) (Cp T_v(k-1) + g weight1 (r(k)-r(k-1)) =
! exner(k-1) (Cp T_v(k-1) - (1-weight1) (r(k)-r(k-1))
!
! where the weights represent the interpolation coefficients
! required to calculate exner at a T point.
! Note: In the new integration scheme Theta below level 1 is assumed
!       to be equal to theta at level 1. The equation is thus
!
! exner(1) (Cp T_v(1) + g (r(k)-r(k-1)) =
! exner* (Cp T_v(1) )
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! Rho level 1
!--------------------------------------------------------------------

k = 1
Do i = 1, pstar_out % level_size

  temp = T % Data(i,k) * (1. + (1./epsilon -1.) * q % Data(i,k))
  exner % Data(i,k) = (CP * temp *                                &
                      (pstar_out % Data(i,1) / pref)**kappa ) /   &
                      (CP * temp + g * (r_rho_levels(i,k) -       &
                        r_theta_levels(i,k-1) ) )

End Do

!--------------------------------------------------------------------
! Rho levels 2 to model_levels + 1
!--------------------------------------------------------------------
Do k = 2, exner % levels
  If ( k <= Output_Grid % wet_levels + 1 ) Then  ! Wet Level
    Do i = 1, exner % level_size
      If (k > Output_Grid % model_levels) Then

        ! extra pressure levels is same height above top theta as
        ! pressure below is below top theta - hence weights are 0.5
        weight1 = 0.5
        deltaz  = 2.0 * (r_theta_levels(i,k-1) - r_rho_levels(i,k-1))

      Else

        deltaz  = r_rho_levels(i,k) - r_rho_levels(i,k-1)
        weight1 = (r_rho_levels(i,k) - r_theta_levels(i,k-1) ) / deltaz

      End If ! (K > model_levels)

      weight2 = 1. - weight1

      temp = T % Data(i,k-1) * (1. + (1./epsilon -1.) * q % Data(i,k-1))

      exner % Data(i,k) = (CP * temp - g * weight1 * deltaz) *         &
                          exner % Data(i,k-1) /                        &
                          (CP * temp + g * weight2 * deltaz)

    End Do
  Else    ! Dry Levels
    Do i = 1, exner % level_size
      If (k > Output_Grid % model_levels ) Then

        ! extra pressure levels is same height above top theta as
        ! pressure below is below top theta - hence weights are 0.5
        weight1 = 0.5
        deltaz = 2.0 * (r_theta_levels(i,k-1) - r_rho_levels(i,k-1))

      Else
        deltaz  = r_rho_levels(i,k) - r_rho_levels(i,k-1)
        weight1 = (r_rho_levels(i,k) - r_theta_levels(i,k-1) ) / deltaz

      End If ! (k > model_levels)

      weight2 = 1. - weight1

      temp = T % Data(i,k-1)
      exner % Data(i,k) = (CP * temp - g * weight1 * deltaz) *       &
                          exner % Data(i,k-1) /                      &
                          (CP * temp + g * weight2 * deltaz)

    End Do
  End If ! wet/dry layer
End Do ! k


!--------------------------------------------------------------------
! Clear up memory for the pstar fields
!--------------------------------------------------------------------
Call Rcf_Dealloc_Field( pstar_in )
Call Rcf_Dealloc_Field( pstar_out )

Return
End Subroutine Rcf_Calc_Output_Exner
End Module Rcf_Calc_Output_Exner_Mod
