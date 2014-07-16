
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  vertical interpolation of some cloud variables

Module Rcf_Vert_Cloud_Mod

!  Subroutine Rcf_Vert_Cloud
!
! Description:
!   Performs vertical interpolation for certain cloud variables
!   (convective cloud top and bottom for example)
!
! Method:
!    The variables store a model level - thus the nearest output level
!    to the input level is calculated as the interpolated level.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   08/05/01   Copy correct data field if only horizontal
!                    interpolation required.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Vert_Cloud( field_in, field_out, grid_in, grid_out, &
                           heights_in, heights_out )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_v_only,                   &
    interp_h_only,                   &
    interp_all,                      &
    interp_copy

Implicit None

! Arguments
Type( field_type ), Intent(In)    :: field_in
Type( field_type ), Intent(InOut) :: field_out
Type( grid_type ),  Intent(In)    :: grid_in
Type( grid_type ),  Intent(In)    :: grid_out
Real,               Intent(In)    :: heights_in(field_in % level_size,&
                                       0 : grid_in % model_levels)
Real,               Intent(In)    :: heights_out(field_out %level_size,&
                                       0 : grid_out % model_levels)
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

! Local variables
Integer                           :: ErrorStatus
Integer                           :: i
Integer                           :: j
Real                              :: level_height_in
Character (Len=*), Parameter      :: RoutineName = 'Rcf_Vert_Cloud'
Character (Len=80)                :: Cmessage

!------------------------------------------------------------------
! Initial checks
!------------------------------------------------------------------
If ( field_out % levels /= 1 .OR. field_in % levels /= 1 ) Then
  ErrorStatus = 10
  Cmessage = 'Cloud base/top interpolation requires single level field'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!------------------------------------------------------------------
! If no interpolation required, just copy
!------------------------------------------------------------------
Select Case( field_in % interp )
  Case( interp_h_only, interp_copy )
    If (field_in % level_size /= field_out % level_size ) Then
      ErrorStatus = 20
      Cmessage = 'Unable to copy data as initial and final grid sizes&
                 & are different'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    If ( field_in % interp == interp_copy) Then
      field_out % Data_Int(:,:) = field_in % Data_Int(:,:)
    Else
      field_out % Data(:,:) = field_in % Data(:,:)
    End If

  Case( interp_all, interp_v_only )

!-------------------------------------------------------------------
! Find output height nearest to input height
!-------------------------------------------------------------------
    Do i = 1, field_out % level_size
      ! Do nothing if level is 0
      If ( Nint( field_in % Data(i,1) ) == 0 ) Then
        field_out % Data(i,1) = 0.0
      Else
        ! find the input height
        level_height_in = heights_in( i, Nint( field_in % Data(i,1) ) )

        ! set output height to missing data
        field_out % Data(i,1) = RMDI

        Do j = 1, grid_out % wet_levels - 1
          If ( level_height_in >= heights_out(i, j) .AND.    &
               level_height_in < heights_out(i, j+1) ) Then

            field_out % Data(i,1) = Real(j)
            Exit
          End If

        End Do

        If (field_out % Data(i,1) == RMDI) Then
          field_out % Data(i,1) = Real( grid_out % model_levels )
        End If

     End If
   End Do

End Select

Return
End Subroutine Rcf_Vert_Cloud
End Module Rcf_Vert_Cloud_Mod
