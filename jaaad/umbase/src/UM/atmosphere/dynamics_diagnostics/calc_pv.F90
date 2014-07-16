#if defined(A15_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine Calc_PV

      Subroutine Calc_PV                                                &
     &                  (u, v, theta, rho,                              &
     &                   r_theta_levels, r_rho_levels,                  &
     &                   r_at_u, r_at_v,                                &
     &                   sec_v_latitude, tan_v_latitude,                &
     &                   sec_theta_latitude, f3_at_v,                   &
     &                   delta_lambda, delta_phi,                       &
     &                   row_length, rows, n_rows, model_levels,        &
     &                   off_x, off_y, halo_i, halo_j,                  &
     &                   at_extremity,                                  &
     &                   pv)

! Description:
!          Calculates potential vorticity at PV points; ie. midway
!          horizontally between v points on rho levels.
!
! Method:
!          Discretisation of equation 7 in UMDP 13.
!
! Original Progammer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Version Date     Comment
! ------- ----     -------
! 5.2     15/11/00 ND version imported into UM. Z. Gardner
!  5.4     28/08/02    Bug Fix (Bi-cyclic LAM)           Carol Roadnight
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      Implicit None

#include "cntlatm.h"
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
                        ! number of points on a row
     &, rows                                                            &
                        ! number of rows of data
     &, n_rows                                                          &
                        ! number of rows of data on a v row
     &, model_levels                                                    &
                        ! number of levels of data
     &, off_x                                                           &
     &, off_y                                                           &
     &, halo_i                                                          &
     &, halo_j

      Real                                                              &
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
     &  theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)

      Real                                                              &
     &  sec_v_latitude (1-off_x:row_length+off_x,                       &
     &                  1-off_y:n_rows+off_y)                           &
     &, tan_v_latitude(row_length,n_rows)                               &
     &, sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)

      Real                                                              &
           ! Coriolis term
     &  f3_at_v (1-off_x:row_length+off_x, 1-off_y: n_rows+off_y)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  pv (row_length, n_rows, model_levels)

! Local variables

      Integer                                                           &
     & i, j, k

      Real                                                              &
     &  recip_delta_lambda                                              &
     &, recip_delta_phi                                                 &
     &, weight1                                                         &
     &, weight2

      Real                                                              &
     &  r_at_pv(row_length,n_rows, model_levels)                        &
     &, r_at_pv_on_theta_levels(row_length,n_rows,model_levels)         &
     &, dtheta_dr(row_length+1,n_rows+1)                                &
     &, du_dr(row_length,n_rows+1)                                      &
     &, dv_dr(row_length+1,n_rows)                                      &
     &, dtheta_dx(row_length,n_rows+1)                                  &
     &, dtheta_dx_ave(row_length,n_rows,model_levels)                   &
     &, du_dr_ave(row_length,n_rows)                                    &
     &, dv_dr_ave(row_length,n_rows)                                    &
     &, dtheta_dy(row_length+1,n_rows)                                  &
     &, dtheta_dy_ave(row_length,n_rows,model_levels)                   &
     &, dtheta_dr_ave(row_length,n_rows)                                &
     &, x_term(row_length,n_rows)                                       &
     &, y_term(row_length,n_rows)                                       &
     &, z_term(row_length,n_rows)                                       &
     &, r_at_u_on_theta_levels(row_length,n_rows+1)                     &
     &, density(row_length+1,n_rows+1)                                  &
     &, density_ave(row_length,n_rows)

! ----------------------------------------------------------------------
! Section 0. Initialisation
! ----------------------------------------------------------------------

      recip_delta_lambda = 1./ delta_lambda
      recip_delta_phi = 1./ delta_phi

! ----------------------------------------------------------------------
! Section 1. Calculate horizontal part of terms
! ----------------------------------------------------------------------

! calculate dtheta/dx and dtheta/dy
      Do k = 1, model_levels
        Do j = 1, n_rows+1
          Do i = 1, row_length
            dtheta_dx(i,j) = (theta(i+1,j,k) - theta(i,j,k)) *          &
     &                          recip_delta_lambda *                    &
     &                          sec_theta_latitude(i,j) * 2. /          &
     &                         (r_theta_levels(i+1,j,k) +               &
     &                          r_theta_levels(i,j,k) )
          End Do
        End Do

        If (model_domain  /=  mt_bi_cyclic_LAM) then
        If (at_extremity(PSouth)) Then
          Do i = 1, row_length
            dtheta_dx(i,1) = 0.
          End Do
        End If
        If (at_extremity(PNorth)) Then
          Do i = 1, row_length
            dtheta_dx(i,rows) = 0.
          End Do
        End If
        Endif

        Do j = 1, n_rows
          Do i = 1, row_length + 1
            dtheta_dy(i,j) = (theta(i,j+1,k) - theta(i,j,k)) *          &
     &                          recip_delta_phi * 2. /                  &
     &                         (r_theta_levels(i,j+1,k) +               &
     &                          r_theta_levels(i,j,k) )
          End Do
        End Do

        Do j = 1, n_rows
          Do i = 1, row_length
            dtheta_dx_ave(i,j,k) = 0.5 * (dtheta_dx(i,j+1) +            &
     &                                    dtheta_dx(i,j) )
          End Do
        End Do

        Do j = 1, n_rows
          Do i = 1, row_length
            dtheta_dy_ave(i,j,k) = 0.5 * (dtheta_dy(i+1,j) +            &
     &                                    dtheta_dy(i,j) )
          End Do
        End Do

        Do j = 1, n_rows + 1
          Do i = 1, row_length
            r_at_u_on_theta_levels(i,j) =                               &
     &                      (r_theta_levels(i+1,j,k) +                  &
     &                       r_theta_levels(i,j,k) ) * 0.5
          End Do
        End Do

        Do j = 1, n_rows
          Do i = 1, row_length
            r_at_pv(i,j,k) = 0.5 * (r_at_u(i,j+1,k) + r_at_u(i,j,k))
            r_at_pv_on_theta_levels(i,j,k) = 0.5 *                      &
     &                             ( r_at_u_on_theta_levels(i,j+1) +    &
     &                               r_at_u_on_theta_levels(i,j) )
          End Do
        End Do

! end loop over levels
      End Do

! loop over model levels for rest of code.
      Do k = 1, model_levels

! ----------------------------------------------------------------------
! Section 1.1 Calculate x_term and y_term.
! ----------------------------------------------------------------------

        If (k  ==  1) Then
          Do j = 1, n_rows
            Do i = 1, row_length
              x_term(i,j) = - dtheta_dx_ave(i,j,k)
              y_term(i,j) = dtheta_dy_ave(i,j,k)
            End Do
          End Do

        Else

          Do j = 1, n_rows
            Do i = 1, row_length
              weight1 = (r_at_pv_on_theta_levels(i,j,k) -               &
     &                   r_at_pv(i,j,k) ) /                             &
     &                  (r_at_pv_on_theta_levels(i,j,k) -               &
     &                   r_at_pv_on_theta_levels(i,j,k-1) )
              weight2 = 1.0 - weight1

              x_term(i,j) = - ( weight2 * dtheta_dx_ave(i,j,k-1) +      &
     &                          weight1 * dtheta_dx_ave(i,j,k) )
              y_term(i,j) = ( weight2 * dtheta_dy_ave(i,j,k-1) +        &
     &                        weight1 * dtheta_dy_ave(i,j,k) )
            End Do
          End Do

        End If

! ----------------------------------------------------------------------
! Section 1.2 Calculate z_term.
! ----------------------------------------------------------------------

        Do j = 1, n_rows
          Do i = 1, row_length
            z_term(i,j) = f3_at_v(i,j) + (                              &
     &                  (v(i+1,j,k)-v(i,j,k)) * recip_delta_lambda *    &
     &                  sec_v_latitude(i,j) +                           &
     &                   0.5 * (u(i,j+1,k) + u(i,j,k)) *                &
     &                   tan_v_latitude(i,j) -                          &
     &                   (u(i,j+1,k) - u(i,j,k)) * recip_delta_phi      &
     &                   ) / r_at_pv(i,j,k)
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2. Multiply horizontal terms by vertical terms and form
!            full pv.
! ----------------------------------------------------------------------

        If (k  ==  1) Then
! Since theta gradient is un-defined use gradient from levels above.
          Do j = 1, n_rows + 1
            Do i = 1, row_length + 1
              dtheta_dr(i,j) = (theta(i,j,k+1) - theta(i,j,k)) /        &
     &                         (r_theta_levels(i,j,k+1) -               &
     &                          r_theta_levels(i,j,k) )
            End Do
          End Do

! Now average to required location.
          Do j = 1, n_rows
            Do i = 1, row_length
              dtheta_dr_ave(i,j) = 0.25 * (dtheta_dr(i,j) +             &
     &                                     dtheta_dr(i+1,j) +           &
     &                                     dtheta_dr(i,j+1) +           &
     &                                     dtheta_dr(i+1,j+1) )
            End Do
          End Do

        Else If (k  /=  2) Then
! no code required at level 2 since the value is the same as level 1.

! calculate theta gradient.
          Do j = 1, n_rows + 1
            Do i = 1, row_length + 1
              dtheta_dr(i,j) = (theta(i,j,k) - theta(i,j,k-1)) /        &
     &                         (r_theta_levels(i,j,k) -                 &
     &                          r_theta_levels(i,j,k-1) )
            End Do
          End Do

! Now average to required location.
          Do j = 1, n_rows
            Do i = 1, row_length
              dtheta_dr_ave(i,j) = 0.25 * (dtheta_dr(i,j) +             &
     &                                     dtheta_dr(i+1,j) +           &
     &                                     dtheta_dr(i,j+1) +           &
     &                                     dtheta_dr(i+1,j+1) )
            End Do
          End Do

        End If

        If (k  ==  1) Then
          Do j = 1, n_rows
            Do i = 1, row_length + 1
              dv_dr(i,j) = (v(i,j,k+1) - v(i,j,k)) /                    &
     &                     (r_at_v(i,j,k+1) -                           &
     &                      r_at_v(i,j,k) )
            End Do
          End Do
          Do j = 1, n_rows + 1
            Do i = 1, row_length
              du_dr(i,j) = (u(i,j,k+1) - u(i,j,k)) /                    &
     &                     (r_at_u(i,j,k+1) -                           &
     &                      r_at_u(i,j,k) )
            End Do
          End Do

        Else If (k  ==  model_levels) Then

          Do j = 1, n_rows
            Do i = 1, row_length + 1
              dv_dr(i,j) = (v(i,j,k) - v(i,j,k-1)) /                    &
     &                     (r_at_v(i,j,k) -                             &
     &                      r_at_v(i,j,k-1) )
            End Do
          End Do
          Do j = 1, n_rows + 1
            Do i = 1, row_length
              du_dr(i,j) = (u(i,j,k) - u(i,j,k-1)) /                    &
     &                     (r_at_u(i,j,k) -                             &
     &                      r_at_u(i,j,k-1) )
            End Do
          End Do

        Else

          Do j = 1, n_rows
            Do i = 1, row_length + 1
              dv_dr(i,j) = (v(i,j,k+1) - v(i,j,k-1)) /                  &
     &                     (r_at_v(i,j,k+1) -                           &
     &                      r_at_v(i,j,k-1) )
            End Do
          End Do
          Do j = 1, n_rows + 1
            Do i = 1, row_length
              du_dr(i,j) = (u(i,j,k+1) - u(i,j,k-1)) /                  &
     &                     (r_at_u(i,j,k+1) -                           &
     &                      r_at_u(i,j,k-1) )
            End Do
          End Do

        End If

! now average quantities
        Do j = 1, n_rows
          Do i = 1, row_length
            du_dr_ave(i,j) = 0.5 * (du_dr(i,j+1) + du_dr(i,j))
          End Do
        End Do

        Do j = 1, n_rows
          Do i = 1, row_length
            dv_dr_ave(i,j) = 0.5 * (dv_dr(i+1,j) + dv_dr(i,j))
          End Do
        End Do

! convert rho to true density by removing factor of r squared.
        Do j = 1, n_rows + 1
          Do i = 1, row_length + 1
            density(i,j) = rho(i,j,k) / (r_rho_levels(i,j,k) *          &
     &                                   r_rho_levels(i,j,k) )
          End Do
        End Do

        Do j = 1, n_rows
          Do i = 1, row_length
            density_ave(i,j) = 0.25 * (density(i,j+1) + density(i+1,j) +&
     &                                 density(i+1,j+1) + density(i,j) )
          End Do
        End Do

! Calculate full PV.
        Do j = 1, n_rows
          Do i = 1, row_length
            pv(i,j,k) = (x_term(i,j) * dv_dr_ave(i,j) +                 &
     &                   y_term(i,j) * du_dr_ave(i,j) +                 &
     &                   z_term(i,j) * dtheta_dr_ave(i,j) )             &
     &                  / density_ave(i,j)
          End Do
        End Do

! end loop over model levels
      End Do

! end of routine

      Return
      END SUBROUTINE Calc_PV

#endif
