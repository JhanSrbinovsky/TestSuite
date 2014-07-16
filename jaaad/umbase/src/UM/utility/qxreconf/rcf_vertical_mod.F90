#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Wrapper for vertical interpolation

Module Rcf_vertical_Mod

!  Subroutine Rcf_Vertical - wrapper for vertical interpolation
!
! Description:
! This module contains a wrapper subroutine for vertical
! interpolation. Data is left on the current decomposition
! as there should be no spatial dependencies.
!
! Method:
!   If no interp. required, data is copied. Otherwise, loop
!   through levels calling relevant interpolation routine for each.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   05/10/01   Cater for extra rho level at top. D.Robinson
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_vertical( field_in, field_out, grid_in, grid_out, &
                         heights_in, heights_out )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_order

Use Rcf_Submodel_Mod, Only : &
    submodel_ident

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_all,                      &
    interp_v_only,                   &
    interp_h_only,                   &
    interp_copy,                     &
    interp_no_op

Use Rcf_PrintStatus_Mod, Only : &
    LTimer

Use Rcf_Level_Code_Mod, Only : &
    Rcf_Level_Code

Implicit None

! Arguments
Type (field_type), Intent(InOut)  :: field_in
Type (field_type), Intent(InOut)  :: field_out
Type (grid_type), Intent(In)      :: grid_in
Type (grid_type), Intent(In)      :: grid_out
Real                              :: heights_in(field_in % level_size,&
                                       0 : grid_in % model_levels+1)
Real                              :: heights_out(field_out %level_size,&
                                       0 : grid_out % model_levels+1)

! Comdecks
#include "cppxref.h"

! Local Data
Character (Len=*), Parameter      :: RoutineName='Rcf_vertical'
Character (Len=80)                :: Cmessage
Integer                           :: ErrorStatus
Integer                           :: i
Integer                           :: j
Integer                           :: start_level
Integer                           :: end_level

External Timer

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3 )
!-----------------------------------------------------------------
! Is interpolation activated? If not, copy data across is all we
! will do.
!-----------------------------------------------------------------
Select Case( field_in % interp )
  Case( interp_copy, interp_h_only )         ! Do a copy for vertical

    ! sizes should be the same, but will check
    If ( field_in % level_size /= field_out % level_size .OR. &
         field_in % levels /= field_out % levels ) Then
      Cmessage = 'No interpolation, but data field sizes/levels are &
                 &different!'
      ErrorStatus = 10
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    Select Case( field_in % stashmaster % data_type )
      Case( ppx_type_real )
        field_out % Data(:,:) = field_in % Data(:,:)

      Case( ppx_type_int )
        If ( Associated( field_in %  Data ) ) Then
          field_out % Data(:,:) = field_in % Data(:,:)
        Else
          field_out % Data_Int(:,:) = field_in % Data_Int(:,:)
        End If

      Case ( ppx_type_log )
        If ( Associated( field_in % Data ) ) Then
          field_out % Data(:,:) = field_in % Data(:,:)
        Else
          field_out % Data_Log(:,:) = field_in % Data_Log(:,:)
        End If

      Case Default
        Cmessage = 'Unsupported Data-Type'
        ErrorStatus = -20
        Call Ereport( RoutineName, ErrorStatus, Cmessage )

    End Select

  Case( interp_all, interp_v_only )      ! Do vertical interpolation

    ! Check that level_size is the same for input and output
    If ( field_in % level_size /= field_out % level_size ) Then
      Cmessage = 'Input and output level_sizes do not match'
      ErrorStatus = 30
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    ! Find start and end levels of output grid through levcod
    Call Rcf_Level_Code( field_out % stashmaster % lb_code,       &
                         start_level,  grid_out )
    Call Rcf_Level_Code( field_out % stashmaster % lt_code,       &
                         end_level, grid_out )

    ! check that start/end level range matches number of levels
    If ( (end_level - start_level) + 1 /= field_out % levels ) Then
      Cmessage = 'Start and End levels cover a different range to field&
                & levels'
      ErrorStatus=40
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    ! Loop through levels of output grid
    j = 1
    Do i = start_level, end_level

! DEPENDS ON: vert_interp
      Call vert_interp( field_in % Data, field_in % level_size,      &
                        field_in % levels, heights_out(1,i),         &
                        heights_in(1,start_level), v_int_order,      &
                        field_out % Data(1,j) )

      j = j + 1

    End Do

  Case( interp_no_op )
  ! Do nothing

End Select

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4 )

Return
End Subroutine Rcf_vertical

End Module Rcf_vertical_Mod
#endif
