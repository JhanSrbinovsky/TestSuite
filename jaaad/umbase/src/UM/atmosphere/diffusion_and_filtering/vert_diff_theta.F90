#if defined(A13_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine vert_diff_theta

      Subroutine vert_diff_theta                                        &
     &                     (theta, r_theta_levels, r_rho_levels,        &
     &                      sin_theta_latitude,pi,L_ramp,               &
     &                      ramp_lat_radians,                           &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      at_extremity,                               &
     &                      timestep, rows, row_length,                 &
     &                      model_levels, model_domain,                 &
     &                      level_start, level_stop,                    &
     &                      diffusion_coefficient,                      &
     &                      theta_star)

! Purpose:
!          Calculates vertical diffusion increment to theta.
!
! Method:
!          Is described in ;
!
!
! Original Programmer: Terry Davies
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date      Comment
! -------  ------     -------
!LL  5.2  24/10/00    Code introduced                    Andy Malcolm
!   5.4    11/07/02    Latitude dependent option added    Andy Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

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

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y         ! Size of small halo in j.

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, level_start                                                     &
     &, level_stop

      Logical :: L_ramp

      Real                                                              &
     &  timestep                                                        &
     &, diffusion_coefficient                                           &
     &, pi                                                              &
     &, ramp_lat_radians

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)

      Real                                                              &
     & theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &          model_levels)                                           &
     &, sin_theta_latitude(row_length, rows)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  theta_star(1-off_x:row_length+off_x,                            &
     &               1-off_y:rows+off_y, model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k, j0, j1      ! Loop indices

! Local arrays

      Real                                                              &
     &  flux (row_length, rows, 2)                                      &
     &, latitude_theta(rows)                                            &
     &, factor(rows)


      Integer                                                           &
     &  k_up, k_down, k_switch

      logical                                                           &
     &      do_row(rows)

! No External Routines:
! calculate latitude of each row
      Do j = 1, rows
        latitude_theta(j) = asin(sin_theta_latitude(1,j))
        do_row(j)=.false.
      End Do

      j0 = 1
      j1 = rows
      If (model_domain  ==  1) Then
        If (at_extremity(PSouth)) j0 = 2
        If (at_extremity(PNorth)) j1 = rows-1
      End If

      do j=j0,j1
        if(L_ramp)then
        if (latitude_theta(j) > -ramp_lat_radians .and.                 &
     &      latitude_theta(j) < ramp_lat_radians) then
           do_row(j)=.true.
           factor(j)=                                                   &
     &     (cos (latitude_theta(j) ) - cos(ramp_lat_radians) )          &
     &               / ( 1.0 -cos(ramp_lat_radians) )
        endif
        else
          do_row(j)=.true.
          factor(j)= 1.0
        endif
      end do

! ----------------------------------------------------------------------
! Section 1.   Vertical Diffusion increment
! ----------------------------------------------------------------------

      k_up = 2
      k_down = 1
      Do j = j0, j1
        Do i = 1, row_length
          flux(i,j,k_down) = 0.0
        End Do
      End Do

      Do k = level_start, level_stop
        Do j = j0, j1
      if (do_row(j)) then
          Do i = 1, row_length
            flux(i,j,k_up) = (theta(i,j,k+1) - theta(i,j,k))/           &
     &               (r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k))
            theta_star(i,j,k) = theta_star(i,j,k) +                     &
     &                          timestep * diffusion_coefficient *      &
     &    factor(j) *                                                   &
     &                          (flux(i,j,k_up) - flux(i,j,k_down))/    &
     &                  (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
          End Do
      endif
        End Do
        k_switch = k_up
        k_up     = k_down
        k_down   = k_switch
      End Do

! End of routine
      return
      END SUBROUTINE vert_diff_theta

#endif
