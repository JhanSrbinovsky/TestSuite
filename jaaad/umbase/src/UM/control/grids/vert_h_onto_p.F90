#if defined(C92_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine vert_h_onto_p

      Subroutine vert_h_onto_p(                                         &
     &                           data_in, row_length, data_rows,        &
     &                           data_levels, desired_p,                &
     &                           r_rho_levels, r_theta_levels,          &
     &  p_theta_levels,                                                 &
     &                           theta, exner_theta_levels,             &
     &                           exner_rho_levels,                      &
     &                           R, g, Lapse_rate,                      &
     &                           boundary_layer_levels,                 &
     &                           off_x, off_y, halo_i, halo_j,          &
     &                           p_at_data, interp_order,               &
     &                           kappa, p_zero,cp,data_out )

! Purpose:
!          Performs linear interpolation in exner of the height
!          field onto a pressure surface. Where surface is below
!          bottom of model, data isd extrapolated as in the UM.
!
! Method:
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Date     Version     Comment
! ----     -------     -------
!  5.2     Replaced log p by Exner in interpolation  C Wilson
!          Added hydrostatic extrapolation at top and removed mdi
!          Set height dependence(2000m) for extrapolation level at
!          bottom of model domain   C Wilson
!  6.1  15/06/04   Change calc of sub surface Geopotential heights
!                                                    Michael Hughes
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                        ! number of points on a row
     &, data_rows                                                       &
                        ! number of rows of data
     &, data_levels                                                     &
                        ! number of levels of data
     &, interp_order                                                    &
                        ! 1 = linear, 3 = cubic, 5=quintic
     &, boundary_layer_levels                                           &
     &, off_x, off_y                                                    &
                        ! small halo sizes
     &, halo_i, halo_j  ! large halo sizes

      Real                                                              &
     &  desired_p                                                       &
                        ! desired value to which data should be
                        ! interpolated to.
     &, desired_exner                                                   &
     &, mdi                                                             &
                        ! missing data indicator.
     &, R                                                               &
     &, g                                                               &
     &, Lapse_rate                                                      &
     &, p_zero                                                          &
     &, kappa                                                           &
     &, cp

      Real                                                              &
     &  data_in (row_length, data_rows, data_levels)                    &
     &, p_at_data (1-off_x:row_length+off_x, 1-off_y:data_rows+off_y,   &
     &             data_levels)                                         &
     &, p_theta_levels (1-off_x:row_length+off_x,                       &
     &                  1-off_y:data_rows+off_y,                        &
     &                  data_levels)                                    &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:data_rows+halo_j,0:data_levels)         &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:data_rows+halo_j, data_levels)            &
     &, theta(1-off_x:row_length+off_x, 1-off_y:data_rows+off_y,        &
     &          data_levels)                                            &
     &, exner_theta_levels(1-off_x:row_length+off_x,                    &
     &                     1-off_y:data_rows+off_y, data_levels)        &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                     1-off_y:data_rows+off_y, data_levels)

! Arguments with Intent OUT. ie: Output variables.
      Real                                                              &
     &  data_out (row_length, data_rows)

! Local variables

      Integer                                                           &
     & i,j,k                                                            &
     & , level_extrap (row_length, data_rows)

      Integer                                                           &
     &  level_below(row_length, data_rows)

      Real                                                              &
     &  power                                                           &
     &, T_ref_level_1                                                   &
     & , extrap_height

! ----------------------------------------------------------------------
! Section 1. Find level below which desired surface is.
!            Use first level above extrap_height if extrapolation needed
! ----------------------------------------------------------------------
      extrap_height = 2000 !m
      Do k = 1, data_levels
        Do j = 1, data_rows
          Do i = 1, row_length
            If ( (r_theta_levels(i,j,k)-r_theta_levels(i,j,0) )  <=     &
     &            extrap_height ) Then
              level_extrap(i,j) = k
            End If
          End Do
        End Do
      End Do

      Do j = 1, data_rows
        Do i = 1, row_length
          level_below(i,j) = 0
        End Do
      End Do

      Do k = 1, data_levels - 1
        Do j = 1, data_rows
          Do i = 1, row_length
            If ( p_at_data(i,j,k)  >=  desired_p ) Then
              level_below(i,j) = k
            End If
          End Do
        End Do
      End Do

! if requested level is above top of model, set to data_levels+1,
! which will be converted to missing data indicator.

      Do j = 1, data_rows
        Do i = 1, row_length
          If ( desired_p  <   p_at_data(i,j,data_levels) ) Then
            level_below(i,j) = data_levels+1
          End If
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 2. Vertical interpolation.
! ----------------------------------------------------------------------

! change desired pressure into exner equivalent.

      desired_exner = (desired_p/p_zero) ** kappa

! required to calculate values below bottom of model
      power = (R * Lapse_Rate) / g

      Do j = 1, data_rows
        Do i = 1, row_length

          If (level_below(i,j)  ==  0) Then

             T_ref_level_1= Lapse_Rate *                                &
     &           ( r_theta_levels(i,j,level_extrap(i,j)) -              &
     &             r_rho_levels(i,j,1)  ) /                             &
     &           (1.- ( p_theta_levels(i,j,level_extrap(i,j)) /         &
     &                  p_at_data(i,j,1) )**power )

            data_out(i,j) = data_in(i,j,1) +                            &
     &                     ( T_ref_level_1 / Lapse_Rate) *              &
     &                     (1. - (desired_p/p_at_data(i,j,1))**power)

          Else If (level_below(i,j)  ==  data_levels+1) Then
            data_out(i,j) = data_in(i,j,data_levels) -                  &
     &            cp*(desired_exner -exner_rho_levels(i,j,data_levels)) &
     &            *theta(i,j,data_levels)/g

          Else
! linearly interpolate
            data_out (i,j) = ( ( desired_exner -                        &
     &                       exner_rho_levels(i,j,level_below(i,j)) )   &
     &                          * data_in (i,j,level_below(i,j)+1)      &
     &                     -( desired_exner -                           &
     &                    exner_rho_levels(i,j,level_below(i,j)+1) ) *  &
     &                           data_in (i,j,level_below(i,j)) ) /     &
     &                   ( exner_rho_levels(i,j,level_below(i,j)+1) -   &
     &                       exner_rho_levels(i,j,level_below(i,j)) )


          End If

        End Do
      End Do

! end of routine

      Return
      END SUBROUTINE vert_h_onto_p

#endif
