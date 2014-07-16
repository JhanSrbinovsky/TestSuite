
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

      Subroutine lbc_u_a_to_c (                                         &
     &           u_p                                                    &
     &,          u                                                      &
     &,          lbc_size_p                                             &
     &,          lbc_size_u                                             &
     &,          lbc_levels                                             &
     &,          lbc_row_len                                            &
     &,          lbc_rows                                               &
     &,          rimwidth                                               &
     &,          halo_x                                                 &
     &,          halo_y                                                 &
     &,          l_var_lbc                                              &
     &,          lambda_p_in                                            &
     &,          lambda_u_in                                            &
     & )

      Implicit None

      Integer ::  lbc_size_p
      Integer ::  lbc_size_u
      Integer ::  lbc_levels
      Integer ::  lbc_row_len
      Integer ::  lbc_rows
      Integer ::  rimwidth
      Integer ::  halo_x
      Integer ::  halo_y

      Real    ::  u_p (lbc_size_p, lbc_levels)
      Real    ::  u   (lbc_size_u, lbc_levels)

      Logical l_var_lbc
      
! Input VarRes grid info in degrees      
      Real :: lambda_p_in ( 1-halo_x: lbc_row_len + halo_x ) 
      Real :: lambda_u_in ( 1-halo_x: lbc_row_len + halo_x )
      
      Integer ::  ipt_p, ipt_u
      Integer ::  lbc_row_len_p
      Integer ::  lbc_row_len_u
      Integer ::  lbc_rows_p
      Integer ::  lbc_rows_u
      Integer ::  level
      Integer ::  offset
      Integer ::  offset_vr 
      
! North

      ipt_p = 1
      ipt_u = 1
      lbc_row_len_p = lbc_row_len + 2*halo_x
      lbc_rows_p    = rimwidth + halo_y + 1
      lbc_row_len_u = lbc_row_len_p - 1
      lbc_rows_u    = lbc_rows_p - 1
      offset        = 1
      offset_vr  = - halo_x
      
      do level =1, lbc_levels

! DEPENDS ON: u_a_to_c
       call u_a_to_c (                                                  &
            u_p(ipt_p,level),                                           &
            u(ipt_u,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_u,                                              &
            lbc_rows_u,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            offset_vr,                                                  &
            lbc_row_len,                                                &
            halo_x,                                                     &
            lambda_p_in,                                                &
            lambda_u_in )

      enddo

! East

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_u = ipt_u + lbc_row_len_u * lbc_rows_u
      lbc_row_len_p = rimwidth + halo_x + 1
      lbc_rows_p    = lbc_rows - 2 * rimwidth
      lbc_row_len_u = lbc_row_len_p - 1
      lbc_rows_u    = lbc_rows_p
      offset        = 0
      offset_vr = lbc_row_len - rimwidth  - 1
      
      do level =1, lbc_levels

! DEPENDS ON: u_a_to_c
       call u_a_to_c (                                                  &
            u_p(ipt_p,level),                                           &
            u(ipt_u,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_u,                                              &
            lbc_rows_u,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            offset_vr,                                                  &
            lbc_row_len,                                                &
            halo_x,                                                     &
            lambda_p_in,                                                &
            lambda_u_in )

      enddo

! South

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_u = ipt_u + lbc_row_len_u * lbc_rows_u
      lbc_row_len_p = lbc_row_len + 2*halo_x
      lbc_rows_p    = rimwidth + halo_y + 1
      lbc_row_len_u = lbc_row_len_p - 1
      lbc_rows_u    = lbc_rows_p - 1
      offset        = 0
      offset_vr = - halo_x
       
      do level =1, lbc_levels

! DEPENDS ON: u_a_to_c
       call u_a_to_c (                                                  &
            u_p(ipt_p,level),                                           &
            u(ipt_u,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_u,                                              &
            lbc_rows_u,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            offset_vr,                                                  &
            lbc_row_len,                                                &
            halo_x,                                                     &
            lambda_p_in,                                                &
            lambda_u_in )

      enddo

! West

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_u = ipt_u + lbc_row_len_u * lbc_rows_u
      lbc_row_len_p = rimwidth + halo_x + 1
      lbc_rows_p    = lbc_rows - 2 * rimwidth
      lbc_row_len_u = lbc_row_len_p - 1
      lbc_rows_u    = lbc_rows_p
      offset        = 0
      offset_vr = - halo_x 
      
      do level =1, lbc_levels

! DEPENDS ON: u_a_to_c
       call u_a_to_c (                                                  &
            u_p(ipt_p,level),                                           &
            u(ipt_u,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_u,                                              &
            lbc_rows_u,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            offset_vr,                                                  &
            lbc_row_len,                                                &
            halo_x,                                                     &
            lambda_p_in,                                                &
            lambda_u_in )

      enddo

      Return
      END SUBROUTINE lbc_u_a_to_c
