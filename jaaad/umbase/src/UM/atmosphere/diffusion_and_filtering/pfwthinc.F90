#if defined(A13_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine pfwthinc

      Subroutine pfwthinc(                                              &
     &                    field, r_theta_levels, r_rho_levels,          &
     &                    off_x, off_y, halo_i, halo_j,                 &
     &                    sec_theta_latitude, cos_v_latitude,           &
     &                    me, n_procy, delta_lambda, delta_phi,         &
     &                    rows, n_rows, row_length, model_levels,       &
     &                    max_filter_rows, u_begin, u_end,              &
     &                    u_sweeps, global_u_filter,                    &
     &                    levels, horizontal_level,                     &
     &                    diff_coeff_phi, diff_coeff_u, L_diff )

! Purpose:
!          Filter/diffusion of a theta-type field ( theta or w )
!          Based on conservative diffusion operator.
!          Applied as multiple EW 1-2-1 filter near poles.
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
! 6.2     14/02/06    New code based on pofil_wth_incs + changes
!                                                         Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  L_diff   ! Horizontal diffusion active

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v-rows.
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
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, me        ! processor id

      Integer                                                           &
     &  max_filter_rows                                                 &
                        !  array dimension for u_begin etc.
     &, u_begin(0:max_filter_rows)                                      &
     &, u_end(0:max_filter_rows)                                        &
     &, u_sweeps(max_filter_rows)                                       &
     &, global_u_filter                                                 &
     &, levels                                                          &
                 ! levels active
     &, horizontal_level   ! level at which steep slope test no
!                          ! longer operates

      Real                                                              &
     &  delta_lambda                                                    &
     &, delta_phi                                                       &
     &, diff_coeff_phi   ! NS diffusion coefficient

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
!   EW  diffusion coefficient
     &, diff_coeff_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  field (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &               levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                    ! Loop indices
     &, j_begin, j_end                                                  &
                          ! Loop bounds
     &, j_store                                                         &
     &, level_flat                                                      &
     &, i_filter                                                        &
     &, i_sweep                                                         &
     &, times                                                           &
     &, num_pass

      Logical                                                           &
     &  L_cycle                                                         &
     &, L_combine

! Local arrays

      Real                                                              &
     &  delta_z(1-off_x:row_length+off_x,                               &
     &          1-off_y:rows+off_y, levels )                            &
     &, recip_r_squared_delz(1-off_x:row_length+off_x,                  &
     &                       1-off_y:rows+off_y, levels )               &
     &, temp(1-off_x:row_length, 1-off_y:rows)                          &
     &, mask_i(1-off_x:row_length, rows, levels)                        &
     &, mask_j(row_length, 1-off_y:rows, levels)

#include "parparm.h"
#include "domtyp.h"

! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set up loop bounds and pointers
!             u_begin and u_end for each potential sweep are set up
!             in SETCON (as u_begin and u_end since on same latitude)
! ----------------------------------------------------------------------
! To switch off slope test completetely then
! either hard-wire level_flat = 1 or ensure that the argument
!    horizontal level <= wet_model_levels - model_levels + 1
!      level_flat = 1
      level_flat =  horizontal_level
      if (level_flat  <=  0) level_flat = 1

! ----------------------------------------------------------------------
! Section 2.   Set field to be filtered and calculate delta_z
!              Diffusion coefficient = 0.25 set in mask array
! ----------------------------------------------------------------------

! calculate D(r)/D(eta) about theta levels
      Do k = 1, levels
        if( k < model_levels )then
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              delta_z(i,j,k) = r_rho_levels(i,j,k+1) -                  &
     &                           r_rho_levels(i,j,k )
              recip_r_squared_delz(i,j,k) = 1.0 /                       &
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
              recip_r_squared_delz(i,j,k) = 1.0 /                       &
     &             ( r_theta_levels(i,j,k) * r_theta_levels(i,j,k) )
            End Do
          End Do
        endif !  k < model_levels
      End Do   !  k = 1, levels

! ----------------------------------------------------------------------
! Section 3.   Switch off diffusion at steep slopes
! ----------------------------------------------------------------------
!  Horizontal_level should be the first flat level and as such
!  does not need testing for a sloping surface.
!   i.e if level_flat < 2   then
!       diffusion behaves as original (along eta surfaces)

      Do k = 1, levels
        Do j = 1, rows
          Do i = 0, row_length
           mask_i(i,j,k) = 1.0
          Enddo
        Enddo
        Do j = 0, rows
          Do i = 1, row_length
           mask_j(i,j,k) = 1.0
          Enddo
        Enddo
      Enddo ! k = 1, levels

      Do k = 1, level_flat - 1
        Do j = 1, rows
          Do i = 0, row_length
            if( r_theta_levels(i,j,k) < r_theta_levels(i+1,j,k-1) .or.  &
     &           r_theta_levels(i+1,j,k) < r_theta_levels(i,j,k-1) )    &
     &       mask_i(i,j,k) = 0.0
          Enddo
        Enddo
        Do j = 0, rows
          Do i = 1, row_length
            if( r_theta_levels(i,j,k) < r_theta_levels(i,j+1,k-1) .or.  &
     &           r_theta_levels(i,j+1,k) < r_theta_levels(i,j,k-1) )    &
     &      mask_j(i,j,k) = 0.0
          Enddo
        Enddo
      Enddo    ! k = 1, level_flat - 1

! ----------------------------------------------------------------------
! Section 4.0  Polar Filtering and diffusion
!              j loops will no-op for filtering
!              All processors active for swap bounds synchronization
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
!        write(6,*)'pofil_wth_incs global_u_filter = ', global_u_filter

      L_combine = .false.
      num_pass = global_u_filter
      j_begin = u_begin(1)
      j_end = u_end(1)
! Combine EW diffusion with first sweep of polar filter
      if ( L_diff ) then
        L_combine = .true.
        if( global_u_filter < 1 ) num_pass = 1
        if ( j_begin < 0 ) then
          j_begin = u_begin(0)
          j_end = u_end(0)
        elseif ( j_end < u_end(0) ) then
          j_end = u_end(0)
        elseif ( j_begin > u_begin(0) ) then
          j_begin = u_begin(0)
        endif !  j_begin < 0
      endif ! L_diff

      DO i_filter = 1, num_pass

        if( i_filter > 1 ) then
          L_combine = .false.
          j_begin = u_begin(i_filter)
          j_end = u_end(i_filter)
        endif !i_filter > 1

        L_cycle = .false.
!  if n_procy = 1, need to filter hemispheres separately
        if( n_procy == 1 )L_cycle = .true.

        i_sweep = 1
        Do  ! Sweeping loop  i_sweep = 1, u_sweeps(i_filter)
!          write(6,*)'  i_filter = ', i_filter ,' i_sweep = ', i_sweep
!     &          ,' j_begin = ', j_begin,' j_end = ', j_end

! ----------------------------------------------------------------------
! Section 4.1   EW filtering and diffusion
! ----------------------------------------------------------------------
          Do k = 1, levels

! There are 2 factors of 1/4 multiplied to give 1/16
            Do j = j_begin, j_end
              Do i = 1-off_x, row_length
                temp(i,j) = (field(i+1,j,k) * delta_z(i+1,j,k) -        &
     &                         field(i,j,k) * delta_z(i,j,k) ) *        &
     &                               mask_i(i,j,k) * diff_coeff_u(i,j) *&
     &        0.25 * (r_theta_levels(i,j,k) + r_theta_levels(i+1,j,k)) *&
     &               (r_theta_levels(i,j,k) + r_theta_levels(i+1,j,k))
              End Do

              Do i = 1, row_length
                field(i,j,k) = field(i,j,k) +                           &
     &                                  recip_r_squared_delz(i,j,k) *   &
     &                                    ( temp(i,j) - temp(i-1,j) )
              End Do
           End Do  ! j = j_begin, j_end

          End Do  ! k = 1, levels

! If n_procy=1, Northern hemisphere needs to be done after S. Hem
          If ( L_cycle ) then
            j_store = j_begin
            if ( i_sweep == 1 .and. L_combine ) then
              j_begin = j_end + 1
            else
              j_begin = rows - j_end + 1
            endif ! i_sweep == 1 .and. L_combine
            j_end = rows - j_store + 1
            L_cycle = .false.
        CYCLE
          endif ! L_cycle)

! Reset pointers for next sweep
! either because N Hem has just been done or
! 1st sweep was combined filter and diffusion
          j_begin = u_begin(i_filter)
          j_end = u_end(i_filter)
          if( n_procy == 1 ) L_cycle = .true.

! DEPENDS ON: swap_bounds
          call Swap_Bounds(                                             &
     &                     field, row_length, rows, levels,             &
     &                     off_x, off_y, fld_type_p, .false.)

          i_sweep = i_sweep + 1
          if ( i_sweep > u_sweeps(i_filter) ) EXIT
        CYCLE
        EndDo ! sweeping  loop  i_sweep = 1, u_sweeps(i_filter)

      EndDO   !  i_filter = 1, num_pass

! ----------------------------------------------------------------------
! Section 4.2  NS diffusion
! ----------------------------------------------------------------------

      if( L_diff ) then
        j_begin = u_begin(0)
        j_end = u_end(0)
        if( n_procy == 1 ) j_end = rows - j_begin + 1
!          write(6,*)' NS diffusion  diff_coeff_phi = ', diff_coeff_phi
!     &          ,' j_begin = ', j_begin,' j_end = ', j_end

        Do k = 1, levels

          Do j = j_begin - 1, j_end
            Do i = 1, row_length
              temp(i,j) = (field(i,j+1,k) * delta_z(i,j+1,k) -          &
     &                       field(i,j,  k) * delta_z(i,j,  k) ) *      &
     &                               mask_j(i,j,k) * diff_coeff_phi *   &
     &        0.25 * (r_theta_levels(i,j,k) + r_theta_levels(i,j+1,k)) *&
     &               (r_theta_levels(i,j,k) + r_theta_levels(i,j+1,k)) *&
     &                                        cos_v_latitude(i,j)
            End Do
          End Do

          Do j = j_begin, j_end
            Do i = 1, row_length
              field(i,j,k) = field(i,j,k) + sec_theta_latitude(i,j) *   &
     &                                   recip_r_squared_delz(i,j,k) *  &
     &                                      (temp(i,j) - temp(i,j-1))
            End Do
          End Do  ! j = j_begin, j_end

        End Do  ! k = 1, levels

! DEPENDS ON: swap_bounds
        call Swap_Bounds(                                               &
     &                   field, row_length, rows, levels,               &
     &                   off_x, off_y, fld_type_p, .false.)
      endif ! L_diff

      return
      END SUBROUTINE pfwthinc

#endif
