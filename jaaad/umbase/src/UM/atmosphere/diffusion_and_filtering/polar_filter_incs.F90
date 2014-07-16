#if defined(A13_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Polar_filter_incs

      Subroutine Polar_filter_incs(                                     &
     &                      theta, theta_star, R_u, R_v, R_w,           &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      north_lat_limit, south_lat_limit,           &
     &                      filt_coeff,                                 &

     &                      n_sweeps, step_per_sweep,                   &
     &                      polar_start_lat_limit,                      &

     &                      off_x, off_y, halo_i, halo_j,               &
     &                      n_proc, n_procx, n_procy, l_datastart,      &
     &                      neighbour, at_extremity, proc_row_group)

! Purpose:
!          Filter fields on rows adjacent to pole and
!          set polar v components to wave 1.
!
! Method:
!          Using 1-2-1 Shapiro filter to zap 2-grid waves
!          T. Davies.
!
! Original Programmer: T. Davies
! Current code owner: T. Davies
!
! History:
! Date     Version     Comment
! ----     -------     -------
!   6.1   04/08/04  argument list change                 Terry  Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels.
     &, global_row_length                                               &
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, l_datastart(3)       ! First gridpoints held by this processor

      Real                                                              &
     &  cos_theta_longitude(row_length,rows)                            &
                                             ! cos of longitude
     &, sin_theta_longitude(row_length,rows)                            &
                                             ! sin of longitude
     &, sin_theta_latitude (row_length, rows)                           &
     &, sin_v_latitude (row_length, n_rows)

      Real                                                              &
           ! limits for filtering (in radians)
     &  north_lat_limit                                                 &
                         ! latitude above which 1-2-1 filter is applied
     &, south_lat_limit                                                 &
                         ! latitude below which 1-2-1 filter is applied
     &, filt_coeff                                                      &
                         !   filter coefficient
     &, step_per_sweep                                                  &
                         ! amount in radians to increment start latitude
                         ! by per sweep
     &, polar_start_lat_limit ! max latitude at which filter can start

      Integer                                                           &
     &  n_sweeps        ! number of sweeps of filter to do

      Integer                                                           &
     &  n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, proc_row_group                                                  &
     &, neighbour(4)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

! Arguments with Intent IN/OUT.

      Real                                                              &
     &  R_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)                                              &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &       model_levels)                                              &
     &, R_w(row_length, rows, model_levels)                             &
!ajm &, R_w(row_length, rows, 0:model_levels)
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
     &, theta_star(1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &             model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                     ! Loop indices
     &, j_start, j_end, i_sweep

      Real                                                              &
     &  u_filt(row_length,rows)                                         &
     &, v_filt(row_length,n_rows)                                       &
     &, w_filt(row_length,rows)                                         &
     &, theta_filt(row_length,rows)                                     &
     &, latitude_theta(rows)                                            &
     &, latitude_v(n_rows)                                              &
     &, R_w_halo(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels-1)                                        &
     &, theta_incs(1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &             model_levels)

      Real                                                              &
     &  scalar                                                          &
                         !   filter coefficient smoother
     &, local_north_lat_limit                                           &
     &, local_south_lat_limit

#include "fldtype.h"

! External Routines:
      External                                                          &
     &  swap_bounds

! ----------------------------------------------------------------------
! Section 0. Calulate latitude for each row.
!            Create theta_increments array and R_w array with halo
! ----------------------------------------------------------------------
! calculate latitude of each row

      Do j = 1, rows
        latitude_theta(j) = asin(sin_theta_latitude(1,j))
      End Do
      Do j = 1, n_rows
        latitude_v(j) = asin(sin_v_latitude(1,j))
      End Do

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_incs(i,j,k) = theta_star(i,j,k) - theta(i,j,k)
          End Do
         End Do
      End Do

      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
            R_w_halo(i,j,k) = R_w(i,j,k)
          End Do
         End Do
      End Do

!ajm    bug found by deborah Salmond, R_u out of bound
! DEPENDS ON: swap_bounds
      Call Swap_Bounds(                                                 &
     &                   R_u, row_length, rows, model_levels,           &
     &                   off_x, off_y, fld_type_u, .true.)
!ajm   R_v may need swapbounding, add in for safety
! DEPENDS ON: swap_bounds
      Call Swap_Bounds(                                                 &
     &                   R_v, row_length, n_rows, model_levels,         &
     &                   off_x, off_y, fld_type_v, .true.)
! DEPENDS ON: swap_bounds
      Call Swap_Bounds(                                                 &
     &                   theta_incs, row_length, rows, model_levels,    &
     &                   off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: swap_bounds
      Call Swap_Bounds(                                                 &
     &                   R_w_halo, row_length, rows, model_levels-1,    &
     &                   off_x, off_y, fld_type_p, .false.)

      Do i_sweep = 1, n_sweeps

        local_north_lat_limit = north_lat_limit + (i_sweep - 1) *       &
     &                          step_per_sweep
        If (local_north_lat_limit  >   polar_start_lat_limit ) Then
          local_north_lat_limit = polar_start_lat_limit
        End If
        local_south_lat_limit = south_lat_limit - (i_sweep - 1) *       &
     &                          step_per_sweep
        If (local_south_lat_limit  <   - polar_start_lat_limit ) Then
          local_south_lat_limit = - polar_start_lat_limit
        End If

! ----------------------------------------------------------------------
! Section 1. 1-2-1 filter desired area near north pole which lies on
!            this processor
! ----------------------------------------------------------------------

! calculate which rows to filter on this processor, default is none
      j_start = 0
      j_end = -1
      Do j = 1, rows
        If (latitude_theta(j)  >   local_north_lat_limit                &
     &       .and. j_start  ==  0) Then
          j_start = j
          j_end = rows
        End If
      End Do

! on north polar processor do not filter polar row
      If (at_extremity(PNorth) .and. j_end  ==  rows) Then
        j_end = rows - 1
      End If

      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 1, row_length
            u_filt(i,j) = filt_coeff * (R_u(i-1,j,k) - R_u(i,j,k) * 2.  &
     &                                    + R_u(i+1,j,k))
            theta_filt(i,j) = filt_coeff*(theta_incs(i-1,j,k)           &
     &                                    - theta_incs(i,j,k) * 2.      &
     &                                    + theta_incs(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!         scalar=1.0
          Do i = 1, row_length
            theta_incs(i,j,k) = theta_incs(i,j,k) + theta_filt(i,j)
            R_u(i,j,k) = R_u(i,j,k) + u_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do
      Do k = 1, model_levels-1
        Do j = j_start, j_end
          Do i = 1, row_length
            w_filt(i,j) = filt_coeff*(R_w_halo(i-1,j,k)                 &
     &                                - R_w_halo(i,j,k) * 2.            &
     &                                + R_w_halo(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!         scalar=1.0
          Do i = 1, row_length
            R_w_halo(i,j,k) = R_w_halo(i,j,k) + w_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do

! calculate which v rows to filter on this processor, default is none
      j_start = 0
      j_end = -1
      Do j = 1, n_rows
        If (latitude_v(j)  >   local_north_lat_limit                    &
     &      .and. j_start  ==  0) Then
          j_start = j
          j_end = n_rows
        End If
      End Do
      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 1, row_length
            v_filt(i,j) = filt_coeff * (R_v(i-1,j,k) - R_v(i,j,k) * 2.  &
     &                                  + R_v(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!            scalar=1.0
          Do i = 1, row_length
             R_v(i,j,k) = R_v(i,j,k) + v_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 2. 1-2-1 filter desired area near south pole which lies on
!            this processor
! ----------------------------------------------------------------------

! calculate which rows to filter on this processor
      If (latitude_theta(1)  <   local_south_lat_limit) Then
! at least one row to filter, default to all
        j_start = 1
        j_end = rows
        Do j = 1, rows
          If (latitude_theta(j)  >=  local_south_lat_limit              &
     &        .and. j_end  ==  rows) Then
! limit filtering to previous row
            j_end = j-1
          End If
        End Do
      Else
! no rows to filter
        j_start = 1
        j_end = 0
      End If

! on north polar processor do not filter polar row
      If (at_extremity(PSouth) .and. j_start  ==  1) Then
        j_start = 2
      End If

      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 1, row_length
            u_filt(i,j) = filt_coeff * (R_u(i-1,j,k) - R_u(i,j,k) * 2.  &
     &                                    + R_u(i+1,j,k))
            theta_filt(i,j) = filt_coeff*(theta_incs(i-1,j,k)           &
     &                                    - theta_incs(i,j,k) * 2.      &
     &                                    + theta_incs(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!         scalar=1.0
          Do i = 1, row_length
            theta_incs(i,j,k) = theta_incs(i,j,k) + theta_filt(i,j)
            R_u(i,j,k) = R_u(i,j,k) + u_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do

      Do k = 1, model_levels-1
        Do j = j_start, j_end
          Do i = 1, row_length
            w_filt(i,j) = filt_coeff*(R_w_halo(i-1,j,k)                 &
     &                                - R_w_halo(i,j,k) * 2.            &
     &                                + R_w_halo(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!         scalar=1.0
          Do i = 1, row_length
            R_w_halo(i,j,k) = R_w_halo(i,j,k) + w_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do

! calculate which rows to filter on this processor
      If (latitude_v(1)  <   local_south_lat_limit) Then
! at least one row to filter, default to all
        j_start = 1
        j_end = n_rows
        Do j = 1, n_rows
          If (latitude_v(j)  >=  local_south_lat_limit                  &
     &        .and. j_end  ==  n_rows) Then
! limit filtering to previous row
            j_end = j-1
          End If
        End Do
      Else
! no rows to filter
        j_start = 1
        j_end = 0
      End If

      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 1, row_length
            v_filt(i,j) = filt_coeff * (R_v(i-1,j,k) - R_v(i,j,k) * 2.  &
     &                                  + R_v(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!            scalar=1.0
          Do i = 1, row_length
             R_v(i,j,k) = R_v(i,j,k) + v_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 3. swop haloes for those fields that have been altered
! ----------------------------------------------------------------------

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   R_u, row_length, rows, model_levels,           &
     &                   off_x, off_y, fld_type_u, .true.)

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   R_v, row_length, n_rows, model_levels,         &
     &                   off_x, off_y, fld_type_v, .true.)

        If (i_sweep  /=  n_sweeps) Then
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   theta_incs, row_length, rows, model_levels,    &
     &                   off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   R_w_halo, row_length, rows, model_levels-1,    &
     &                   off_x, off_y, fld_type_p, .false.)
        End If

      End Do ! end loop over number of sweeps

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_star(i,j,k) = theta(i,j,k) + theta_incs(i,j,k)
          End Do
         End Do
      End Do

      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
            R_w(i,j,k) = R_w_halo(i,j,k)
          End Do
         End Do
      End Do

! DEPENDS ON: swap_bounds
      Call Swap_Bounds(                                                 &
     &                   theta_star, row_length, rows, model_levels,    &
     &                   off_x, off_y, fld_type_p, .false.)

! End of routine
      return
      END SUBROUTINE Polar_filter_incs

#endif
