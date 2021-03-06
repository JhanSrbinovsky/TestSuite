#if defined(A13_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine turb_diff_q

      Subroutine turb_diff_q(                                           &
     &                       moist,                                     &
     &                       r_theta_levels, r_rho_levels,              &
     &                       sec_theta_latitude,                        &
     &                       cos_v_latitude, sec_v_latitude,            &
     &                       off_x, off_y, halo_i, halo_j,              &
     &                       delta_lambda, delta_phi, timestep,         &
     &                       rows, n_rows, row_length,                  &
     &                       model_levels, wet_levels, levels,          &
     &                       diff_coeff_u, diff_coeff_v,                &
     &                       moist_star)

! Purpose:
!          Turbulent diffusion of a theta field
!          Based on conservative diffusion operator.
!
! Method:
!          Is described in ;
!
!
! Original Programmer:   T Davies
! Current code owner: Andrew J. Malcolm
!
! History:
! Version    Date      Comment
! ----     -------     -------
! 6.2     04/10/05    New code based on conservative diffusion. T Davies
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
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, wet_levels                                                      &
                       ! number of wet model levels.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, me                                                              &
                  ! processor id
     &, levels     ! number of levels to process

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi                                                       &
     &, timestep

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)

      Real                                                              &
           !  trigonometric functions
     &  sec_theta_latitude(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y)                          &
     &, cos_v_latitude(1-off_x:row_length+off_x,                        &
     &                     1-off_y:n_rows+off_y)                        &
     &, sec_v_latitude(1-off_x:row_length+off_x,                        &
     &                     1-off_y:n_rows+off_y)                        &
!  diffusion coefficient , only 1 side needs halo
     &, diff_coeff_u(1-off_x:row_length, rows, levels)                  &
     &, diff_coeff_v(row_length, 1-off_y:n_rows+off_y, levels)

      Real                                                              &
     &  moist (1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, wet_levels )

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  moist_star (1-off_x:row_length+off_x,                           &
     &              1-off_y:rows+off_y, wet_levels )

! Local Variables.

      Integer                                                           &
     &  i, j, k     ! Loop indices

      Real                                                              &
     &  recip_delta_lambda                                              &
     &, recip_delta_phi

! Local arrays

      Real                                                              &
     &  delta_z(1-off_x:row_length+off_x,                               &
     &          1-off_y:rows+off_y, levels )                            &
     &, dt_over_r_squared_delz(1-off_x:row_length+off_x,                &
     &                       1-off_y:rows+off_y, levels )               &
     &, temp(1-off_x:row_length+off_x, 1-off_y:rows)                    &
     &, lambda_term(row_length, rows)                                   &
     &, phi_term(row_length, rows)

#include "parparm.h"
#include "domtyp.h"

! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set values and calculate delta_z
! ----------------------------------------------------------------------

      recip_delta_lambda = 1.0 / delta_lambda
      recip_delta_phi = 1.0 / delta_phi

! calculate D(r)/D(eta) about theta levels
      Do k = 1, levels
        if( k < model_levels )then
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
!
! "conservative diffusion" diffuses dzQ rather than Q 
! delta_z is set to 1.0 to avoid instabilities ocurring in the presence
! of orography due to volume differences between neighbouring gridboxes.
!
              delta_z(i,j,k) = 1.0
!              delta_z(i,j,k) = r_rho_levels(i,j,k+1) -                 &
!     &                           r_rho_levels(i,j,k )
              dt_over_r_squared_delz(i,j,k) = timestep /                &
     &             ( r_theta_levels(i,j,k) * r_theta_levels(i,j,k) *    &
     &                                              delta_z(i,j,k) )
            End Do
          End Do
          else ! k = model_levels
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at top level and then no need to take special action in later loops
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              delta_z(i,j,k) = 1.0
              dt_over_r_squared_delz(i,j,k) = timestep /                &
     &             ( r_theta_levels(i,j,k) * r_theta_levels(i,j,k) )
            End Do
          End Do
        endif !  k < model_levels
      End Do   !  k = 1, levels

! ----------------------------------------------------------------------
! Section 2.0  Horizontal Diffusion
! ----------------------------------------------------------------------

      Do k = 1, levels

! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------
        Do j = 1, rows
          Do i = 1-off_x, row_length
            temp(i,j) = (moist(i+1,j,k) * delta_z(i+1,j,k) -            &
     &                   moist(i  ,j,k) * delta_z(i  ,j,k) ) *          &
     &                    recip_delta_lambda * diff_coeff_u(i,j,k) *    &
     &                   sec_v_latitude(i,j) * sec_v_latitude(i,j)
          End Do

          Do i = 1, row_length
            lambda_term(i,j) = recip_delta_lambda *                     &
     &                          ( temp(i,j) - temp(i-1,j))
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------

        Do j = 0, rows
          Do i = 1, row_length
            temp(i,j) = (moist(i,j+1,k) * delta_z(i,j+1,k) -            &
     &                   moist(i,j,  k) * delta_z(i,j,  k) ) *          &
     &                   recip_delta_phi * diff_coeff_v(i,j,k) *        &
     &                                  cos_v_latitude(i,j)
          End Do
        End Do

        Do j = 1, rows
          Do i = 1, row_length
            phi_term(i,j) = (temp(i,j) - temp(i,j-1)) *                 &
     &                        recip_delta_phi * sec_theta_latitude(i,j)
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.3   Calculate new variable.
! ----------------------------------------------------------------------

        Do j = 1, rows
          Do i = 1, row_length
              moist_star(i,j,k) = moist_star(i,j,k) +                   &
     &                            dt_over_r_squared_delz(i,j,k) *       &
     &                         ( lambda_term(i,j) + phi_term(i,j) )
          End Do
        End Do

      End Do ! k = 1, levels

! DEPENDS ON: swap_bounds
      call Swap_Bounds(                                                 &
     &                 moist_star, row_length, rows, levels,            &
     &                 off_x, off_y, fld_type_p, .false.)

      return
      END SUBROUTINE turb_diff_q

#endif
