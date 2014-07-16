
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Interpolations between A and C grids

Module Rcf_AC_interp_Mod

!  Subroutine Rcf_U_to_P
!  Subroutine Rcf_V_to_P
!  Subroutine Rcf_P_to_U
!  Subroutine Rcf_P_to_V
!
! Description:
!   Interpolates full-level (ie no horizontal decomposition) fields
!   between A and C grids.
!
! Method:
!   Linear interpolation in relevant (x or y) direction.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.4   29/08/02   Enable use of bicyclic LAM. C. Roadnight
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

!------------------------------------------------------------------
! Routine to interpolate U field onto P points
!------------------------------------------------------------------
Subroutine Rcf_U_to_P( u_level, u_at_p, grid )

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Implicit None

! Arguments
Type( grid_type ), Intent(In)     :: grid
Real,              Intent(In)     :: u_level( grid % glob_u_row_length,&
                                              grid % glob_u_rows )
Real,              Intent(Out)    :: u_at_p(  grid % glob_p_row_length,&
                                              grid % glob_p_rows )

! Local variables
Integer     :: i
Integer     :: j

Do j = 1, grid % glob_u_rows
  Do i = 2, grid % glob_u_row_length
    u_at_p(i,j) = ( u_level(i-1,j) +  u_level(i,j) ) * 0.5
  End Do

!  If ( grid % glob_u_row_length == grid % glob_p_row_length ) Then
  ! Wrapping LAM

!   u_at_p(1,j) =(u_level(1,j) + u_level( grid % glob_u_row_length, j))&
!                 * 0.5

!  Else    ! non-wrapping

    ! Copy for end P points.
    u_at_p(1,j) = u_level(1,j)
    u_at_p(grid % glob_p_row_length,j) =                               &
                                    u_at_p(grid % glob_u_row_length, j)
!  End If
End Do

Return
End Subroutine Rcf_U_to_P

!------------------------------------------------------------------
! Routine to interpolate V field onto P points
!------------------------------------------------------------------
Subroutine Rcf_V_to_P( v_level, v_at_p, grid )

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_Parvars_Mod, Only : &
    bound,                  &
    bc_cyclic

Implicit None

! Arguments
Type( grid_type ), Intent(In)     :: grid
Real,              Intent(In)     :: v_level( grid % glob_v_row_length,&
                                              grid % glob_v_rows )
Real,              Intent(Out)    :: v_at_p(  grid % glob_p_row_length,&
                                              grid % glob_p_rows )

! Local variables
Integer                           :: i
Integer                           :: j

Do i = 1, grid % glob_v_row_length
  Do j = 2, grid % glob_v_rows
    v_at_p(i,j) = (v_level(i, j-1) + v_level(i,j) ) * 0.5
  End Do

  ! Copy for end points

  If (bound(2) == bc_cyclic)then
    v_at_p(i,1)=(v_level(i,1)+v_level(i,grid % glob_v_rows))* 0.5
  else
  v_at_p(i,1) = v_level(i,1)
  v_at_p(i,grid % glob_p_rows) = v_level(i,grid % glob_v_rows)
  Endif
End Do

Return
End Subroutine Rcf_V_to_P

!------------------------------------------------------------------
! Routine to interpolate P field onto U points
!------------------------------------------------------------------
Subroutine Rcf_P_to_U( u_at_p, u_level, grid )

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Implicit None

! Arguments
Type( grid_type ), Intent(In)     :: grid
Real,              Intent(In)     :: u_at_p(  grid % glob_p_row_length,&
                                              grid % glob_p_rows )
Real,              Intent(Out)    :: u_level( grid % glob_u_row_length,&
                                              grid % glob_u_rows )

! Local variables
Integer                           :: i
Integer                           :: j


Do j = 1, grid % glob_u_rows
  Do i = 1, grid % glob_p_row_length - 1
    u_level(i,j) = ( u_at_p(i,j) +  u_at_p(i+1,j) ) * 0.5
  End Do

!  If ( grid % glob_p_row_length == grid % glob_u_row_length ) Then
    ! Wrapping grid
!   u_level( grid % glob_u_row_length, j ) =                        &
!          (u_at_p( grid % glob_p_row_length, j ) + u_at_p(1, j) ) * 0.5
!  End If
End Do

Return
End Subroutine Rcf_P_to_U

!------------------------------------------------------------------
! Routine to interpolate P field onto V points
!------------------------------------------------------------------
Subroutine Rcf_P_to_V( v_at_p, v_level, grid )

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_Parvars_Mod, Only : &
    bound,                  &
    bc_cyclic


Implicit None

! Arguments
Type( grid_type ), Intent(In)     :: grid
Real,              Intent(In)     :: v_at_p(  grid % glob_p_row_length,&
                                              grid % glob_p_rows )
Real,              Intent(Out)    :: v_level( grid % glob_v_row_length,&
                                              grid % glob_v_rows )

! Local variables
Integer                           :: i
Integer                           :: j

Do i = 1, grid % glob_v_row_length
  Do j = 1, grid % glob_v_rows - 1
    v_level(i,j) = ( v_at_p(i,j) + v_at_p(i,j+1) ) * 0.5
  End Do
End Do
If(bound(2) == bc_cyclic)then
  Do i = 1, grid % glob_v_row_length
    v_level(i,grid % glob_v_rows)=(v_at_p(i,grid % glob_v_rows) + &
                             v_at_p(i,1) ) *0.5
  end do
else
  Do i = 1, grid % glob_v_row_length
    j = grid % glob_v_rows
    v_level(i,j) = ( v_at_p(i,j) + v_at_p(i,j+1) ) * 0.5
  End Do
endif

Return
End Subroutine Rcf_P_to_V


End Module Rcf_AC_interp_Mod
