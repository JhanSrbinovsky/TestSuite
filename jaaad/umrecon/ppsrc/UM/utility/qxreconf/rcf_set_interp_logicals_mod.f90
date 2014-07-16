
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initial setting of interpolation logical switches.

Module Rcf_Set_Interp_Logicals_Mod

!  Subroutine Rcf_Set_Interp_Logicals
!
! Description:
!    Sets the h_int_active, v_int_active and v_int_active_soil
!    logical switches that control horizontal and vertical
!    interpolation. This is just an initial setting and may be
!    changed at a later time.
!
! Method:
!    h_int_active - Checks for domain and row/row-length changes
!    v_int_active - Checks for # of model/wet level changes and
!                   changes to eta values or height field constants
!    v_int_active_soil - Checks for changes in number and values of
!                        soil levels.
!
! Current Code Owner: P.Selwood
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Set_Interp_Logicals( Input_Grid, Output_Grid, Hdr_In )

Use Rcf_Grid_Type_Mod, Only : &
    Grid_Type

Use Rcf_UMhead_Mod, Only : &
    UM_Header_Type

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active,         &
    v_int_active_soil

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Normal,            &
    PrStatus_Oper

Use Rcf_HeadAddress_Mod, Only :                                &
    RC_LongSpacing,        RC_LatSpacing,         RC_FirstLat, &
    RC_FirstLong,          RC_PoleLat,            RC_PoleLong

Use Rcf_readnl_horizont_mod, Only :                             &
    Delta_Lambda,          Delta_Phi,             Lambda_First, &
    Phi_First,             Lambda_Npole,          Phi_Npole

Use Rcf_Generate_Heights_Mod, Only : &
    height_gen_original

Use Rcf_Parvars_Mod, Only : &
    mype

Implicit None

! Arguments
Type( Grid_Type ), Intent(In)       :: Input_Grid
Type( Grid_Type ), Intent(In)       :: Output_Grid
Type( UM_Header_Type ), Intent(In)  :: Hdr_In

! Local variables
Integer                             :: i        ! Looper
Character (Len=40)                  :: reason   ! why is interp.
                                                ! switched on.

!--------------------------------------------------------------------
! Horizontal Interpolation
!--------------------------------------------------------------------
! Default
h_int_active = .FALSE.

If ( .NOT. Output_Grid % Global ) Then
  If ( abs( Hdr_In % RealC( RC_FirstLat )  - Phi_First ) >      &
                                            Epsilon( 1.0 ) ) Then
    h_int_active = .TRUE.
    reason = 'differing first latitude'
  End If
End If

If ( abs( Hdr_In % RealC( RC_FirstLong )  - Lambda_First )  >   &
                                            Epsilon( 1.0 ) .OR. &
     abs( Hdr_In % RealC( RC_PoleLat )    - Phi_Npole )     >   &
                                            Epsilon( 1.0 ) .OR. &
     abs( Hdr_In % RealC( RC_PoleLong)    - Lambda_Npole )   >  &
                                            Epsilon( 1.0 ) ) Then
  h_int_active = .TRUE.
  reason = 'changed domain'
End If

If ( Input_Grid % glob_p_rows /= Output_Grid % glob_p_rows .OR. &
     Input_Grid % glob_p_row_length /=                          &
                        Output_Grid % glob_p_row_length ) Then
  h_int_active = .TRUE.
  reason = 'changed horizontal resolution'
End If

If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*)
  If (h_int_active) Then
    Write (6,*) 'Horizontal interpolation is switched ON'
    If (PrintStatus >= PrStatus_Oper) Then
      Write (6,*) 'Because of ', reason
    End If
  Else
    Write (6,*) 'Horizontal interpolation is switched OFF'
  End If
End If

!---------------------------------------------------------------------
! Vertical interpolation
!---------------------------------------------------------------------
v_int_active = .FALSE.

! Turn on if number of levels varies
If (Input_Grid % model_levels /= Output_Grid % model_levels .OR.    &
    Input_Grid % wet_levels /= Output_Grid % wet_levels) Then

  v_int_active = .TRUE.
  reason = 'changed number of model or wet levels'

Else   ! Turn on if eta levels vary
  Do i = 1, Hdr_In % Len1LevDepC
    If ( Abs( Hdr_In % LevDepC( i, 1 ) -                            &
         Output_Grid % eta_theta_levels( i-1)) > Epsilon( 1.0 )) Then

      v_int_active = .TRUE.
      reason = 'change in eta_theta levels'
     End If
  End Do

  Do i = 1, Hdr_In % Len1LevDepC - 1
    If ( Abs( Hdr_In % LevDepC( i, 2 ) -                            &
         Output_Grid % eta_rho_levels( i ) ) > Epsilon( 1.0 ) ) Then

      v_int_active = .TRUE.
      reason = 'change in eta_rho levels'
    End If
  End Do
End If

! Turn on if height specifying constants vary
If ( (Input_Grid  % first_constant_r_rho_level /=     &
      Output_Grid % first_constant_r_rho_level ) .OR. &
     (Input_Grid % z_top_of_model - Output_Grid % z_top_of_model ) > &
                                    Epsilon( 1.0 ) ) Then
  v_int_active = .TRUE.
  reason = 'change in height defining constants'
End If

! Turn on if number of boundary levels change
If ( Input_Grid % Height_Gen_Method == height_gen_original .AND. &
     Input_Grid % BL_Levels /= Output_Grid % BL_Levels ) Then
  v_int_active = .TRUE.
  reason = 'change in number of boundary levels'
End If

! Turn on if change in way heights are calculated
If ( Input_Grid  % Height_Gen_Method /=     &
     Output_Grid % Height_Gen_Method ) Then
  v_int_active = .TRUE.
  reason = 'change in method of height generation'
End If

If (v_int_active) Then
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
    Write (6,*) 'Vertical interpolation is switched ON'
    If (PrintStatus >= PrStatus_Oper) Then
      Write (6,*) 'Because of ', reason
    End If
  End If
Else
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
    Write (6,*) 'Vertical interpolation is switched OFF'
  End If
End If

!---------------------------------------------------------------------
! Vertical interpolation for soil levels
!---------------------------------------------------------------------
v_int_active_soil = .FALSE.
If ( Input_Grid % sm_levels /= Output_Grid % sm_levels) Then
  v_int_active_soil = .TRUE.
  reason = 'changed number of soil depths'
Else
  Do i = 1, Input_Grid % sm_levels
    If ( Abs( Input_Grid % soil_depths(i) -                          &
              Output_Grid % soil_depths(i) ) > Epsilon(1.0) ) Then

      v_int_active_soil = .TRUE.
      reason = 'changed soil depth'
    End If
  End Do
End If

If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  If (v_int_active_soil) Then
    Write (6,*) 'Vertical interpolation between soil depths is &
                 &switched ON'
    If (PrintStatus >= PrStatus_Oper) Then
      Write (6,*) 'Because of ', reason
    End If
  Else
    Write (6,*) 'Vertical interpolation between soil depths is &
                  &switched OFF'
  End If
    Write (6,*)
End If

Return
End Subroutine Rcf_Set_Interp_Logicals
End Module Rcf_Set_Interp_Logicals_Mod
