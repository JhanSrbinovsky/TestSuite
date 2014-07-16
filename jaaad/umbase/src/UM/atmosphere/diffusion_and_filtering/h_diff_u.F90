#if defined(A13_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine H_diff_u

      Subroutine H_diff_u(                                              &
     &                      u, r_at_u, ground,                          &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      at_extremity, proc_row_group,               &
     &                      n_proc, n_procx, n_procy, neighbour,        &
     &                      delta_lambda, delta_phi,                    &
     &                      timestep, rows, row_length,                 &
     &                      model_levels, model_domain,                 &
     &                      diffusion_coefficient,                      &
     &                      diffusion_order,                            &
     &                      R_u,                                        &
     &                      horizontal_level)

! Purpose:
!          Calculates horizontal diffusion increment to u.
!
! Method:
!          Is described in ;
!
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Date     Version     Comment
! ----     -------     -------
!LL   5.1   11/02/00  Use DOMTYP parameters                    P.Burton
!LL  5.2 24/10/00  Add steep slope diffusion code          Andy Malcolm
!  5.2  27/09/00  Bug fixed                                   A. Malcolm
!  5.3  23/01/02  correct out of bounds memory call         A. Malcolm
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
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, proc_row_group                                                  &
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)         ! Array with the Ids of the four neighbours
                             ! in the horizontal plane

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, diffusion_order(model_levels)

      Real                                                              &
     &  timestep

      Real                                                              &
     &  diffusion_coefficient(model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_at_u (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:rows+halo_j, model_levels)                     &
     &, ground (1-halo_i:row_length+halo_i,                             &
     &                   1-halo_j:rows+halo_j)

      Real                                                              &
     &  u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &     model_levels)

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  R_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &        model_levels)


      Integer                                                           &
     &  horizontal_level  !  level at which steep slope test no
!                         !  longer operates

! Local Variables.

      Integer                                                           &
     &  i, j, k, ka, order, j0, j1      ! Loop indices


      Real                                                              &
     &  scalar1                                                         &
     &, scalar2                                                         &
     &, sign                                                            &
     &, lower_surface

! Local arrays

      Real                                                              &
     &  lambda_term(row_length, rows, model_levels)                     &
     &, phi_term(row_length, rows, model_levels)                        &
     &, timestep_over_r_squared(row_length, rows, model_levels)         &
     &, field(1-off_x:row_length+off_x,                                 &
     &        1-off_y:rows+off_y, model_levels)

      Integer                                                           &
     &  active_diffusion_order(model_levels)                            &
     &, n_active_levels                                                 &
     &, max_diffusion_order
      Integer                                                           &
     &  mask_i(1-off_x:row_length, rows, model_levels)                  &
     &, mask_j(row_length, 1-off_y:rows, model_levels)                  &
     &, level_base                                                      &
     &, level_flat


      Real                                                              &
     &  active_diff_coeff(model_levels)

#include "parparm.h"
#include "domtyp.h"

! No External Routines:

      scalar1 = 1. / (delta_lambda * delta_lambda)
      scalar2 = 1. / (delta_phi * delta_phi)
      j0 = 1
      j1 = rows
      If (model_domain  ==  mt_Global) Then
        If (at_extremity(PSouth)) j0 = 2
        If (at_extremity(PNorth)) j1 = rows-1
      End If

! ----------------------------------------------------------------------
! Section 1.   Copy u into field for those model levels with a
!              positive diffusion coefficient.
!              Copy coeffs and diffusion order for all active levels.
! ----------------------------------------------------------------------

      n_active_levels = 0
      max_diffusion_order = 0
      Do k = 1, model_levels
        If (diffusion_coefficient(k)  >   0.0 .and.                     &
     &      diffusion_order(k)  >   0 )Then
          n_active_levels = n_active_levels + 1
          active_diffusion_order(n_active_levels) = diffusion_order(k)
          active_diff_coeff(n_active_levels) = diffusion_coefficient(k)
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              field(i,j,n_active_levels) = u(i,j,k)
            End Do
          End Do
          Do j = 1, rows
            Do i = 1, row_length
              timestep_over_r_squared(i,j,n_active_levels) =            &
     &                        timestep /                                &
     &                     ( r_at_u(i,j,k) * r_at_u(i,j,k) )
            End Do
          End Do
          If (diffusion_order(k)  >   max_diffusion_order )             &
     &        max_diffusion_order = diffusion_order(k)
        End If
      End Do

        level_base = model_levels - n_active_levels
!  Last model-level upon which diffusion is inactive
! ----------------------------------------------------------------------
! Section 2.   Loop over the order to produce order * del squared
! ----------------------------------------------------------------------
!  Horizontal_level should be the first flat level and as such
!  does not need testing for a sloping surface.
!   i.e if level_flat < 2   then
!       diffusion behaves as original (along eta surfaces)
! To switch off slope test completetely then
! either hard-wire level_flat = 1 or ensure that the argument
!    horizontal level <= wet_model_levels - n_active_levels + 1
!      level_flat = 1
       level_flat = n_active_levels -                                   &
     &              model_levels + horizontal_level
       if (level_flat  <=  0) level_flat = 1

       Do k = 1, n_active_levels
        Do j = 1, rows
          Do i = 0, row_length
           mask_i(i,j,k)=1
          Enddo
        Enddo
        Do j = 0, rows
          Do i = 1, row_length
           mask_j(i,j,k)=1
          Enddo
        Enddo
      Enddo

      Do ka = 1, level_flat - 1
        k = ka + level_base
        Do j = 1, rows
          Do i = 0, row_length
          if(k == 1)then
            lower_surface=( ground(i,j)+ground(i+1,j) )/2.
          else
            lower_surface=r_at_u(i,j,k-1)
          endif
           if(r_at_u(i+1,j,k)  <   r_at_u(i,j,k+1)                      &
     &        .and.                                                     &
     &        r_at_u(i+1,j,k)  >   lower_surface)then
!ajm     &        r_at_u(i+1,j,k)  >   r_at_u(i,j,k-1))then
           mask_i(i,j,ka)=1
           else
           mask_i(i,j,ka)=0
           endif
          Enddo
        Enddo
        Do j = 0, rows
          Do i = 1, row_length
          if(k == 1)then
            lower_surface=( ground(i,j)+ground(i+1,j) )/2.
          else
            lower_surface=r_at_u(i,j,k-1)
          endif
           if(r_at_u(i,j+1,k)  <   r_at_u(i,j,k+1)                      &
     &        .and.                                                     &
     &        r_at_u(i,j+1,k)  >   lower_surface)then
!ajm     &        r_at_u(i,j+1,k)  >   r_at_u(i,j,k-1))then
           mask_j(i,j,ka)=1
           else
           mask_j(i,j,ka)=0
           endif
          Enddo
        Enddo
      Enddo

      Do order = 1, max_diffusion_order

! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------

        Do k = 1, n_active_levels
! Only do calculations for levels whose diffusion order is less than
! the current value
          If (order  <=  active_diffusion_order(k)) Then
            Do j = j0, j1
              Do i = 1, row_length
                lambda_term(i,j,k) = (field(i+1,j,k) - field(i,j,k)) *  &
     &                            mask_i(i,j,k)*active_diff_coeff(k)    &
     &                         - (field(i,j,k) - field(i-1,j,k)) *      &
     &                            mask_i(i-1,j,k)*active_diff_coeff(k)
              End Do
            End Do

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------

            Do j = j0, j1
              Do i = 1, row_length
                phi_term(i,j,k) = (field(i,j+1,k) - field(i,j,k)) *     &
     &                         mask_j(i,j,k)*active_diff_coeff(k)       &
     &                      - (field(i,j,k) - field(i,j-1,k)) *         &
     &                         mask_j(i,j-1,k)*active_diff_coeff(k)
              End Do
            End Do

          End If ! On if (order  <=  active_diffusion_order(k))
        End Do ! end loop over model levels

! ----------------------------------------------------------------------
! Section 2.3   Calculate new variable.
! ----------------------------------------------------------------------

        Do k = 1, n_active_levels
          If (order  <=  active_diffusion_order(k)) Then
            Do j = j0 ,j1
              Do i = 1, row_length
                field(i,j,k) = timestep_over_r_squared(i,j,k) *         &
     &                     (lambda_term(i,j,k) * scalar1 +              &
     &                      phi_term(i,j,k) * scalar2 )
              End Do
            End Do
          End If
        End Do

        If (model_domain  ==  mt_Global) Then
          If (at_extremity(PSouth)) Then
            Do k = 1, n_active_levels
              Do i = 1, row_length
                field(i,1,k) = 0.
              End Do
            End Do
          End If

          If (at_extremity(PNorth)) Then
            Do k = 1, n_active_levels
              Do i = 1, row_length
                field(i,rows,k) = 0.
              End Do
            End Do
          End If
        End If

        If (order  /=  max_diffusion_order ) then
! DEPENDS ON: swap_bounds
          call Swap_Bounds(                                             &
     &                   field, row_length, rows, n_active_levels,      &
     &                   off_x, off_y, fld_type_u, .true.)
        End If

! End loop over order
      End Do

! ----------------------------------------------------------------------
! Section 3.   Diffusion increment is field * appropriate sign.
! ----------------------------------------------------------------------

! use ka to calculate which level of field maps to which level of
! input data.
      ka = 0
      Do k = 1, model_levels
        If (diffusion_coefficient(k)  >   0.0 .and.                     &
     &      diffusion_order(k)  >   0 )Then
          ka = ka + 1
          sign = (-1) ** (diffusion_order(k)-1)
          Do j = j0, j1
            Do i = 1, row_length
              R_u(i,j,k) = R_u(i,j,k) + field(i,j,ka) * sign
            End Do
          End Do
        End If
      End Do

! End of routine
      return
      END SUBROUTINE H_diff_u

#endif
