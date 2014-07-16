
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


! *****************************************************************
! Routine to interpolate u-comp from p-grid to u-grid
! *****************************************************************

      Subroutine u_a_to_c (                                             &
     &           u_p                                                    &
     &,          u                                                      &
     &,          lbc_row_len_p                                          &
     &,          lbc_rows_p                                             &
     &,          lbc_row_len_u                                          &
     &,          lbc_rows_u                                             &
     &,          offset                                                 &
     &,          l_var_lbc                                              &
     &,          offset_vr                                              &
     &,          lbc_row_len                                            &
     &,          halo_x                                                 &
     &,          lambda_p_in                                            &
     &,          lambda_u_in                                            &
     & )    

      Implicit None

      Integer ::  lbc_row_len_p
      Integer ::  lbc_row_len_u
      Integer ::  lbc_rows_p
      Integer ::  lbc_rows_u
      Integer ::  offset
      Integer ::  lbc_row_len
      Integer ::  halo_x      
      Integer ::  offset_vr  
      
      Real    ::  u_p (lbc_row_len_p, lbc_rows_p)
      Real    ::  u   (lbc_row_len_u, lbc_rows_u)
      
      Logical  l_var_lbc
! Input VarRes grid info in degrees      
      REAL  :: Lambda_p_in ( 1-halo_x: lbc_row_len + halo_x )  
      REAL  :: Lambda_u_in ( 1-halo_x: lbc_row_len + halo_x )
      
      Integer :: i,j
      Real  weight1, weight2 
      
      If(l_var_lbc) Then
                     
        Do i = 1, lbc_row_len_u
          weight1=(lambda_u_in(offset_vr+i) - lambda_p_in(offset_vr+i)) &
     &         /(lambda_p_in(offset_vr+i+1) - lambda_p_in(offset_vr+i))  
          weight2= 1.0 - weight1
         
          Do j = 1, lbc_rows_u        
            u(i,j) = weight1 *  u_p(i,offset+j)                         &
     &              + weight2 * u_p(i+1,offset+j)
          End Do
        End Do

      Else       

        Do j = 1, lbc_rows_u
          Do i = 1, lbc_row_len_u
            u(i,j) = 0.5 * ( u_p(i,offset+j) + u_p(i+1,offset+j) )
          End Do
        End Do

      End if 
      
      Return
      END SUBROUTINE u_a_to_c
