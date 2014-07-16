
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

      Subroutine lbc_v_a_to_c (                                         &
     &           v_p                                                    &
     &,          v                                                      &
     &,          lbc_size_p                                             &
     &,          lbc_size_v                                             &
     &,          lbc_levels                                             &
     &,          lbc_row_len                                            &
     &,          lbc_rows                                               &
     &,          rimwidth                                               &
     &,          halo_x                                                 &
     &,          halo_y                                                 &
     &,          l_var_lbc                                              &
     &,          phi_p_in                                               &
     &,          phi_v_in                                               &
     & )

      Implicit None

      Integer ::  lbc_size_p    ! grid size on p grid
      Integer ::  lbc_size_v    ! grid size on v grid
      Integer ::  lbc_levels    ! no of lbc levels
      Integer ::  lbc_row_len   ! row length of lbc grid
      Integer ::  lbc_rows      ! rows on lbc grid
      Integer ::  rimwidth      !
      Integer ::  halo_x
      Integer ::  halo_y
      
      Real    ::  v_p (lbc_size_p, lbc_levels)  ! lbc v on p grid
      Real    ::  v   (lbc_size_v, lbc_levels)  ! lbc v on v grid

      Logical l_var_lbc
      
! Input VarRes grid info in degrees
      Real    ::  phi_p_in ( 1-halo_y: lbc_rows + halo_y )
      Real    ::  phi_v_in ( 1-halo_y: lbc_rows + halo_y )

! Local variables
      Integer ::  ipt_p, ipt_v
      Integer ::  lbc_row_len_p
      Integer ::  lbc_row_len_v
      Integer ::  lbc_rows_p
      Integer ::  lbc_rows_v
      Integer ::  level
      Integer ::  offset
      Integer ::  offset_vr 
      
! To be replaced with loop over iside (5.3)
! Set up row length, rows and start addresses in parlbcs

! North

      ipt_p = 1
      ipt_v = 1
      lbc_row_len_p = lbc_row_len + 2*halo_x
      lbc_rows_p    = rimwidth + halo_y + 1
      lbc_row_len_v = lbc_row_len_p
      lbc_rows_v    = lbc_rows_p - 1
      offset        = 0
      offset_vr = lbc_rows - rimwidth - 1
       
      do level =1, lbc_levels

! DEPENDS ON: v_a_to_c
       call v_a_to_c (                                                  &
            v_p(ipt_p,level),                                           &
            v(ipt_v,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_v,                                              &
            lbc_rows_v,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            offset_vr,                                                  &
            lbc_rows,                                                   &
            halo_y,                                                     &
            phi_p_in,                                                   &
            phi_v_in )

      enddo

! East

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_v = ipt_v + lbc_row_len_v * lbc_rows_v
      lbc_row_len_p = rimwidth + halo_x + 1
      lbc_rows_p    = lbc_rows - 2 * rimwidth
      lbc_row_len_v = lbc_row_len_p - 1
      lbc_rows_v    = lbc_rows_p - 1
      offset        = 1
      offset_vr = rimwidth
      
      do level =1, lbc_levels

! DEPENDS ON: v_a_to_c
       call v_a_to_c (                                                  &
            v_p(ipt_p,level),                                           &
            v(ipt_v,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_v,                                              &
            lbc_rows_v,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            offset_vr,                                                  &
            lbc_rows,                                                   &
            halo_y,                                                     &
            phi_p_in,                                                   &
            phi_v_in )

      enddo

! South

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_v = ipt_v + lbc_row_len_v * lbc_rows_v
      lbc_row_len_p = lbc_row_len + 2*halo_x
      lbc_rows_p    = rimwidth + halo_y + 1
      lbc_row_len_v = lbc_row_len_p
      lbc_rows_v    = lbc_rows_p - 1
      offset        = 0
      offset_vr =  - halo_y
       
      do level =1, lbc_levels

! DEPENDS ON: v_a_to_c
       call v_a_to_c (                                                  &
            v_p(ipt_p,level),                                           &
            v(ipt_v,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_v,                                              &
            lbc_rows_v,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            offset_vr,                                                  &
            lbc_rows,                                                   &
            halo_y,                                                     &
            phi_p_in,                                                   &
            phi_v_in )

      enddo

! West

      ipt_p = ipt_p + lbc_row_len_p * lbc_rows_p
      ipt_v = ipt_v + lbc_row_len_v * lbc_rows_v
      lbc_row_len_p = rimwidth + halo_x + 1
      lbc_rows_p    = lbc_rows - 2 * rimwidth
      lbc_row_len_v = lbc_row_len_p - 1
      lbc_rows_v    = lbc_rows_p - 1
      offset        = 0
      offset_vr = rimwidth
       
      do level =1, lbc_levels

! DEPENDS ON: v_a_to_c
       call v_a_to_c (                                                  &
            v_p(ipt_p,level),                                           &
            v(ipt_v,level),                                             &
            lbc_row_len_p,                                              &
            lbc_rows_p,                                                 &
            lbc_row_len_v,                                              &
            lbc_rows_v,                                                 &
            offset,                                                     &
            l_var_lbc,                                                  &
            offset_vr,                                                  &
            lbc_rows,                                                   &
            halo_y,                                                     &
            phi_p_in,                                                   &
            phi_v_in )

      enddo

      Return
      END SUBROUTINE lbc_v_a_to_c
