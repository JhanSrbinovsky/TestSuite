
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interpolate wind components from P to U/V grids.
!
! This decks contains 4 subroutines :-
!
! lbc_v_a_to_c : controls interpolation of v-component
!                           from p-grid to v-grid.
! lbc_u_a_to_c : controls interpolation of u-component
!                           from p-grid to u-grid.
! v_a_to_c     : interpolates v-component from p-grid to v-grid
! u_a_to_c     : interpolates u-component from p-grid to u-grid
!
! Description:
!   Interpolates u-component and v-component from the p-grid to
!   u-grid and v-grid respectively.
!
! Method:
!   u-comp and v-comp are both at all p-grid points after being rotated.
!   The u_comp is interpolated in the x_direction to the u-grid.
!   The v_comp is interpolated in the y_direction to the v-grid.
!   Linear interpolation is used. The u-grid and v-grid are within
!   the p-grid so there is no approximation round the edges.
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!   5.3    08/06/01   Remove duplicate declarations.  A van der Wal
!   5.3    09/10/01   Correct bugs introduced at 5.2. D Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

! *****************************************************************
! Routine to control v-comp interpolation from p-grid to v-grid
! *****************************************************************


! *****************************************************************
! Routine to control u-comp interpolation from p-grid to u-grid
! *****************************************************************


! *****************************************************************
! Routine to interpolate v-comp from p-grid to v-grid
! *****************************************************************

      Subroutine v_a_to_c (                                             &
     &           v_p                                                    &
     &,          v                                                      &
     &,          lbc_row_len_p                                          &
     &,          lbc_rows_p                                             &
     &,          lbc_row_len_v                                          &
     &,          lbc_rows_v                                             &
     &,          offset                                                 &
     &,          l_var_lbc                                              &
     &,          offset_vr                                              &
     &,          lbc_rows                                               &
     &,          halo_y                                                 &
     &,          phi_p_in                                               &
     &,          phi_v_in                                               &
     & )

      Implicit None

      Integer ::  lbc_row_len_p
      Integer ::  lbc_row_len_v
      Integer ::  lbc_rows_p
      Integer ::  lbc_rows_v
      Integer ::  offset
      Integer ::  lbc_rows  
      Integer ::  halo_y      
      Integer ::  offset_vr 
      
      Real    ::  v_p (lbc_row_len_p, lbc_rows_p)
      Real    ::  v   (lbc_row_len_v, lbc_rows_v)

      Logical l_var_lbc
! Input VarRes grid info in degrees  
         
      Real    ::  phi_p_in ( 1-halo_y: lbc_rows + halo_y )   
      Real    ::  phi_v_in ( 1-halo_y: lbc_rows + halo_y ) 
     
      Integer :: i,j
      Real  weight1, weight2 

      If(l_var_lbc) Then
      
        Do j = 1, lbc_rows_v 
          weight1=(phi_v_in(offset_vr+j) - phi_p_in(offset_vr+j))/      &
     &           (phi_p_in(offset_vr+j+1) - phi_p_in(offset_vr+j))
          weight2= 1.0 - weight1
          Do i = 1, lbc_row_len_v
            v(i,j) = weight1 * v_p(offset+i,j)                          &
     &              + weight2 * v_p(offset+i,j+1) 
          End Do
        End Do

      Else      

        Do j = 1, lbc_rows_v
          Do i = 1, lbc_row_len_v
            v(i,j) = 0.5 * ( v_p(offset+i,j) + v_p(offset+i,j+1) )
          End Do
        End Do

      End if
      
      Return
      END SUBROUTINE v_a_to_c
