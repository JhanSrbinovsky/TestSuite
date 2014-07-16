#if defined(A13_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine NI_filter_Ctl

      Subroutine NI_filter_Ctl(                                         &
     &                      theta, u, v, w, Exner, rho,                 &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      r_theta_levels, r_rho_levels,               &
     &                      r_at_u, r_at_v, delta_lambda, delta_phi,    &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      cos_theta_latitude, sec_theta_latitude,     &
     &                      cos_v_latitude, sec_v_latitude,             &
     &                      north_lat_limit, south_lat_limit,           &
     &                      polar_filter_coefficient,                   &
     &                      polar_filter_n_sweeps, step_per_sweep,      &
     &                      polar_start_lat_limit,                      &
     &                      max_filter_rows, u_sweeps, v_sweeps,        &
     &                      global_u_filter, global_v_filter,           &
     &                      u_begin, u_end, v_begin, v_end,             &
     &                      diff_coeff_phi, diff_coeff_u, diff_coeff_v, &
     &                      diff_coeff_thermo, diff_coeff_wind,         &
     &                      diff_order_thermo, diff_order_wind,         &
     &                      first_constant_r_rho_level,                 &
     &                      first_constant_r_rho_level_m1,              &
     &                      top_filt_start, top_filt_end,               &
     &                      up_diff, max_upd_levels,                    &
     &                      horizontal_level, me,                       &
     &                      global_row_length, global_rows,             &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      n_proc, n_procx, n_procy, l_datastart,      &
     &                      neighbour, at_extremity, model_domain,      &
     &                      proc_row_group, proc_col_group,             &
     &                      L_polar_filter, L_pofil_new,                &
     &                      L_pfcomb, L_pftheta, L_pfuv,                &
     &                      L_pfw, L_pfexner, L_diff_exner,             &
     &                      L_diff_thermo, L_diff_wind, L_diff_w,       &
     &                      L_pofil_hadgem2, Ltimer,exner_theta_levels, &
#include "argsts.h"
     &                      STASHwork)

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
! Version    Date      Comment
! ----     -------     -------
!  6.1     04/08/04    Code Introduced            Terry Davies
!  6.2     25/12/05    Filter modifications           Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

!#include "parvars.h"
#include "csubmodl.h"
#include "typsts.h"

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  L_polar_filter                                                  &
                       ! switch for polar filter
     &, L_pofil_new                                                     &
                       ! switch for new polar filter
     &, Ltimer                                                          &
     &, L_pfcomb                                                        &
                       ! switch for combined polar filter/diffusion
     &, L_pftheta                                                       &
                       ! switch for polar filter for theta
     &, L_pfw                                                           &
                       ! switch for polar filter for w
     &, L_pfuv                                                          &
                       ! switch for polar filter for horizontal winds
     &, L_pfexner                                                       &
                       ! switch for polar filter for Exner pressure
     &, L_diff_Exner                                                    &
                       ! switch for diffusion of Exner pressure
     &, L_diff_thermo                                                   &
                       ! switch for horiz. diffusion of theta
     &, L_diff_wind                                                     &
                       ! switch for horiz. diffusion of u,v
     &, L_diff_w       ! switch for horiz. diffusion of w
      Logical :: L_pofil_hadgem2  ! use hadgem2 polar filtering settings

      Integer                                                           &
     &  max_filter_rows                                                 &
                           ! max array size for u_begin etc
     &, max_upd_levels                                                  &
                         ! max no. levels for upper diffusion
     &, global_u_filter                                                 &
                         ! number of filter sweeps; 0=diffusion only
     &, global_v_filter                                                 &
                         ! number of filter sweeps; 0=diffusion only
     &, row_length                                                      &
                         ! number of points on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels.
     &, global_row_length                                               &
     &, first_constant_r_rho_level                                      &
     &, first_constant_r_rho_level_m1                                   &
     &, global_rows                                                     &
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, l_datastart(3)                                                  &
                             ! First gridpoints held by this processor
     &, me                                                              &
                  !  this processor
     &, model_domain

      Real                                                              &
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_at_u (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:rows+halo_j, model_levels)                     &
     &, r_at_v (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:n_rows+halo_j, model_levels)
      Real                                                              &
     &  cos_theta_longitude(row_length,rows)                            &
                                             ! cos of longitude
     &, sin_theta_longitude(row_length,rows)                            &
                                             ! sin of longitude
     &, sin_theta_latitude (row_length, rows)                           &
     &, sin_v_latitude (row_length, n_rows)                             &
     &, cos_theta_latitude(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y )                         &
     &, sec_theta_latitude(1-off_x:row_length+off_x,                    &
     &                     1-off_y:rows+off_y )                         &
     &, cos_v_latitude(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y ) &
     &, sec_v_latitude(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y ) &
     &, up_diff(max_upd_levels)                                         &
                                 !  upper-level diffusion coefficients
!    diffusion coefficient for u/theta rows
     &, diff_coeff_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y )     &
!    diffusion coefficient for u/theta rows
     &, diff_coeff_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y )

      Real                                                              &
           ! limits for polar filtering (in radians)
     &  north_lat_limit                                                 &
     &, south_lat_limit                                                 &
     &, polar_filter_coefficient                                        &
                                       !   filter coefficient
     &, step_per_sweep                                                  &
                         ! amount in radians to increment start latitude
                         ! by per sweep
     &, polar_start_lat_limit                                           &
                              ! max latitude at which filter can start
     &, diff_coeff_phi                                                  &
                         ! NS diffusion coeff
     &, diff_coeff_thermo                                               &
                            ! NS diffusion coeff thermo
     &, diff_coeff_wind     ! NS diffusion coeff u,v



      Integer                                                           &
     &  n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, proc_row_group                                                  &
     &, diff_order_thermo                                               &
                                ! diffusion order
     &, diff_order_wind                                                 &
                              ! diffusion order
     &, proc_col_group                                                  &
     &, polar_filter_n_sweeps                                           &
                                     ! number of sweeps of filter to do
     &, horizontal_level                                                &
                                     ! steep slope control
     &, top_filt_start                                                  &
                             ! start level for upper-level diffusion
     &, top_filt_end         ! end level for upper-level diffusion

      Integer                                                           &
     &  neighbour(4)                                                    &
     &, u_sweeps(max_filter_rows)                                       &
                                   ! sweeps for 1-2-1 filter
     &, v_sweeps(max_filter_rows)                                       &
                                   ! sweeps for 1-2-1 filter
     &, u_begin(0:max_filter_rows)                                      &
                                    ! row pointers for 1-2-1 filter
     &, u_end(0:max_filter_rows)                                        &
                                    ! row pointers for 1-2-1 filter
     &, v_begin(0:max_filter_rows)                                      &
                                    ! row pointers for 1-2-1 filter
     &, v_end(0:max_filter_rows)    ! row pointers for 1-2-1 filter

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN/OUT.

      Real                                                              &
     &  u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &     model_levels)                                                &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      0:model_levels)                                             &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
     &, exner (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels + 1)                                        &
     &, rho (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)                                              &
     &, exner_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,&
     &         model_levels )

!   Array arguments with intent(out):
      REAL, INTENT(OUT) ::                                              &
     & STASHwork(*)   ! Output array holding diagnostic fields

#include "parparm.h"
#include "domtyp.h"

! Local variables
      Integer                                                           &
     &  i, j, k                                                         &
                          ! loop counters
     &, j_start, j_stop                                                 &
                          ! Loop bounds
     &, gi                                                              &
                       ! pointer relative to origin
     &, active_levels  ! number of levels for upper-level diffusion

      Real                                                              &
     &  pole_term     ! factor at pole calculated once

      Real                                                              &
     &  mag_vector_np (model_levels)                                    &
     &, dir_vector_np (model_levels)                                    &
     &, mag_vector_sp (model_levels)                                    &
     &, dir_vector_sp (model_levels)

      Real,DIMENSION(:,:,:),ALLOCATABLE ::                              &
     &  r_theta_at_u                                                    &
     &, r_theta_at_v                                                    &
     &, r_uv_b

      Real :: u_inc (row_length, rows, model_levels) 
      Real :: v_inc (row_length, n_rows, model_levels) 
      Real :: w_inc (row_length, rows, 0:model_levels) 
      Real :: T_inc (row_length, rows, model_levels) 
      Real :: exner_inc (row_length, rows, model_levels+1) 
      Integer :: sect
      Integer :: item
      Integer :: im_index      !  internal model index for STASH arrays
      INTEGER :: Errorstatus = 0  ! initial value for error code
      CHARACTER (LEN=80) :: CMessage !  Error message

! ----------------------------------------------------------------------
! 0.1  Section for setting up diagnostic arrays
! ----------------------------------------------------------------------
      sect=13
      im_index    = internal_model_index(atmos_im)
      Cmessage    = ''

      IF (sf(0,sect)) THEN
        IF(sf(381,sect))t_inc(:,:,:)= theta(1:row_length, 1:rows, :)
        IF(sf(385,sect))u_inc(:,:,:)=u(1:row_length, 1:rows, :)
        IF(sf(386,sect))v_inc(:,:,:)=v(1:row_length, 1:n_rows, :)
        IF(sf(387,sect))w_inc(:,:,:)=w(1:row_length, 1:rows, :)
        IF(sf(388,sect))exner_inc(:,:,:)= exner(1:row_length, 1:rows, :)
      ENDIF

! ----------------------------------------------------------------------
! 1.0  Section for original polar filter
! ----------------------------------------------------------------------

      If (L_polar_filter) Then

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('Polar_Filter',3)

! DEPENDS ON: polar_filter
        Call Polar_filter(                                              &
     &                      theta, u, v, w,                             &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      global_row_length, delta_lambda,            &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      north_lat_limit, south_lat_limit,           &
     &                      polar_filter_coefficient,                   &
     &                      polar_filter_n_sweeps, step_per_sweep,      &
     &                      polar_start_lat_limit,                      &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      n_proc, n_procx, n_procy, l_datastart,      &
     &                      neighbour, at_extremity, proc_row_group)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('Polar_Filter',4)

      End If    !  L_polar_filter

! ----------------------------------------------------------------------
! 2.0  Section for new polar filter
! ----------------------------------------------------------------------

      ALLOCATE ( r_theta_at_u(1-off_x:row_length+off_x,                 &
     &                        1-off_y:rows+off_y, 0:model_levels) )
      ALLOCATE ( r_theta_at_v(1-off_x:row_length+off_x,                 &
     &                        1-off_y:n_rows+off_y, 0:model_levels) )

      Do k = 0, model_levels
        Do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
! end loop bound OK since r_theta_levels is halo_i > off_x
            r_theta_at_u(i,j,k) = .5 * (r_theta_levels(i+1,j,k) +       &
     &                                  r_theta_levels(i  ,j,k) )
          End Do
        End Do
      End Do ! k = 0, model_levels
! No need for swap_bounds as halos have been calculated explicitly

      Do k = 0, model_levels
        Do j = 1-off_y, n_rows+off_y
! end loop bound OK since r_theta_levels is halo_j > off_y
          Do i = 1-off_x, row_length+off_x
            r_theta_at_v(i,j,k) = .5 * ( r_theta_levels(i,j+1,k) +      &
     &                                   r_theta_levels(i,j  ,k) )
          End Do
        End Do
      End Do ! k = 0, model_levels
! No need for swap_bounds as halos have been calculated explicitly

! ----------------------------------------------------------------------
! 2.1  polar filter/diffusion theta/w  polar filter Exner
! ----------------------------------------------------------------------

! Swap_bounds are done inside the polar filter sweeps
      if( L_pftheta .or. L_diff_thermo )then
      if ( L_pofil_new ) then
! DEPENDS ON: pofil_new
        Call pofil_new(                                                 &
     &                theta, fld_type_p, 0, 0, 0, 1, 0, 0,              &
     &                model_levels, model_levels, model_levels,         &
     &                first_constant_r_rho_level_m1,                    &
     &                rows, n_rows, rows, row_length,                   &
     &                r_theta_levels, r_rho_levels,                     &
     &                r_theta_at_u, r_theta_at_v,                       &
     &                off_x, off_y, off_x, off_y,                       &
     &                halo_i, halo_j, halo_i, halo_j,                   &
     &                off_x, off_y, off_x, off_y,                       &
     &                sec_theta_latitude, cos_v_latitude, 0.0,          &
     &                model_domain, at_extremity, n_procy,              &
     &                max_filter_rows, global_u_filter,                 &
     &                u_sweeps, u_begin, u_end,                         &
     &                horizontal_level, diff_coeff_phi,                 &
     &                diff_coeff_u, L_diff_thermo, .false.,             &
     &                L_pofil_hadgem2 )
      else
! DEPENDS ON: pofil_th_field
        Call pofil_th_field(                                            &
     &                  theta, r_theta_levels, r_rho_levels,            &
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  sec_theta_latitude, cos_v_latitude,             &
     &                  me, n_procy, delta_lambda, delta_phi,           &
     &                  rows, n_rows, row_length, model_levels,         &
     &                  max_filter_rows, u_begin, u_end, u_sweeps,      &
     &                  global_u_filter, model_levels, horizontal_level,&
     &                  diff_order_thermo, diff_coeff_thermo,           &
     &                  diff_coeff_u, L_diff_thermo )
      endif !  L_pofil_new
      endif !L_pftheta .or. L_diff_thermo

      if( L_pfw .or. L_diff_w )then
      if ( L_pofil_new ) then
! DEPENDS ON: pofil_new
        Call pofil_new(                                                 &
     &                w(1-off_x,1-off_y, 1), fld_type_p,                &
     &                0, 0, 0, 1, 0, 0,                                 &
     &                model_levels, model_levels, model_levels-1,       &
     &                first_constant_r_rho_level_m1,                    &
     &                rows, n_rows, rows, row_length,                   &
     &                r_theta_levels, r_rho_levels,                     &
     &                r_theta_at_u, r_theta_at_v,                       &
     &                off_x, off_y, off_x, off_y,                       &
     &                halo_i, halo_j, halo_i, halo_j,                   &
     &                off_x, off_y, off_x, off_y,                       &
     &                sec_theta_latitude, cos_v_latitude, 0.0,          &
     &                model_domain, at_extremity, n_procy,              &
     &                max_filter_rows, global_u_filter,                 &
     &                u_sweeps, u_begin, u_end,                         &
     &                horizontal_level, diff_coeff_phi,                 &
     &                diff_coeff_u, L_diff_w, .false.,                  &
     &                L_pofil_hadgem2 )
      else
! DEPENDS ON: pofil_th_field
        Call pofil_th_field(                                            &
     &                  w(1-off_x,1-off_y, 1),                          &
     &                  r_theta_levels, r_rho_levels,                   &
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  sec_theta_latitude, cos_v_latitude,             &
     &                  me, n_procy, delta_lambda, delta_phi,           &
     &                  rows, n_rows, row_length, model_levels,         &
     &                  max_filter_rows, u_begin, u_end, u_sweeps,      &
     &                  global_u_filter,                                &
     &                  model_levels - 1, horizontal_level,             &
     &                  diff_order_thermo, diff_coeff_thermo,           &
     &                  diff_coeff_u, L_diff_w )
      endif !  L_pofil_new
      endif ! L_pfw .or. L_diff_w

      if( L_pfexner .and. L_pofil_new )then
! DEPENDS ON: pofil_new
        Call pofil_new(                                                 &
     &                Exner, fld_type_p, 0, 0, 1, 0, 1, 1,              &
     &                model_levels+1, model_levels, model_levels+1,     &
     &                first_constant_r_rho_level,                       &
     &                rows, n_rows, rows, row_length,                   &
     &                r_rho_levels, r_theta_levels,                     &
     &                r_at_u, r_at_v,                                   &
     &                off_x, off_y, off_x, off_y,                       &
     &                halo_i, halo_j, halo_i, halo_j,                   &
     &                halo_i, halo_j, halo_i, halo_j,                   &
     &                sec_theta_latitude, cos_v_latitude, 0.0,          &
     &                model_domain, at_extremity, n_procy,              &
     &                max_filter_rows, global_u_filter,                 &
     &                u_sweeps, u_begin, u_end,                         &
     &                horizontal_level, diff_coeff_phi,                 &
     &                diff_coeff_u, L_diff_Exner, .false.,              &
     &                L_pofil_hadgem2 )
      endif !  L_pfexner .and. L_pofil_new

      if( L_pfuv .or. L_diff_wind )then

        ALLOCATE ( r_uv_b(1-off_x:row_length+off_x,                     &
     &                    1-off_y:n_rows+off_y, model_levels) )

!  r needed at centre of grid but no swap bound option for grid-type
!      so fill required halos explicitly
!      (fortunately have large halos to start)
        Do k = 1, model_levels
          Do j = 1-off_y, n_rows+off_y
! end loop bound OK since r_rho_levels is halo_j > off_y
            Do i = 1-off_x, row_length+off_x
              r_uv_b(i,j,k) = .5 * ( r_at_u(i,j,k) + r_at_u(i,j+1,k) )
            End Do
          End Do
        End Do ! k = 1, model_levels

! ----------------------------------------------------------------------
! 2.2  polar filter/diffusion u/v
! ----------------------------------------------------------------------
      if ( L_pofil_new ) then
! DEPENDS ON: pofil_new
        Call pofil_new(                                                 &
     &                 u, fld_type_u, 1, 0, 1, 0, 1, 1,                 &
     &                 model_levels, model_levels, model_levels,        &
     &                 first_constant_r_rho_level,                      &
     &                 rows, n_rows, rows, row_length,                  &
     &                 r_at_u, r_theta_at_u, r_rho_levels, r_uv_b,      &
     &                 off_x, off_y, off_x, off_y,                      &
     &                 halo_i, halo_j, off_x, off_y,                    &
     &                 halo_i, halo_j, off_x, off_y,                    &
     &                 sec_theta_latitude, cos_v_latitude, 0.0,         &
     &                 model_domain, at_extremity, n_procy,             &
     &                 max_filter_rows, global_u_filter,                &
     &                 u_sweeps, u_begin, u_end,                        &
     &                 horizontal_level, diff_coeff_phi,                &
     &                 diff_coeff_u, L_diff_wind, .true.,               &
     &                 L_pofil_hadgem2 )
      else
! DEPENDS ON: pofil_u
        Call pofil_u(                                                   &
     &              u, r_at_u, r_theta_levels, r_rho_levels,            &
     &              off_x, off_y, halo_i, halo_j,                       &
     &              sec_theta_latitude, cos_v_latitude,                 &
     &              me, n_procy, delta_lambda, delta_phi,               &
     &              rows, n_rows, row_length, model_levels,             &
     &              max_filter_rows, u_begin, u_end, u_sweeps,          &
     &              global_u_filter,                                    &
     &              horizontal_level, diff_order_wind,                  &
     &              diff_coeff_wind, diff_coeff_u, L_diff_wind )
      endif !  L_pofil_new

      if ( L_pofil_new ) then
! DEPENDS ON: pofil_new
      Call pofil_new(                                                   &
     &               v, fld_type_v, 0, 1, 1, 0, 1, 1,                   &
     &               model_levels, model_levels, model_levels,          &
     &               first_constant_r_rho_level,                        &
     &               n_rows, rows, rows, row_length,                    &
     &               r_at_v, r_theta_at_v, r_uv_b, r_rho_levels,        &
     &               off_x, off_y, off_x, off_y,                        &
     &               halo_i, halo_j, off_x, off_y,                      &
     &               off_x, off_y, halo_i, halo_j,                      &
     &               sec_v_latitude, cos_theta_latitude, 0.0,           &
     &               model_domain, at_extremity, n_procy,               &
     &               max_filter_rows, global_v_filter,                  &
     &               v_sweeps, v_begin, v_end,                          &
     &               horizontal_level, diff_coeff_phi,                  &
     &               diff_coeff_v, L_diff_wind, .true.,                 &
     &               L_pofil_hadgem2 )
      else
! DEPENDS ON: pofil_v
      Call pofil_v(                                                     &
     &              v, r_at_v, r_theta_levels, r_rho_levels,            &
     &              off_x, off_y, halo_i, halo_j,                       &
     &              sec_v_latitude, cos_theta_latitude,                 &
     &              l_datastart, global_row_length,                     &
     &              delta_lambda, delta_phi,                            &
     &              me, at_extremity, proc_row_group, n_procy,          &
     &              rows, n_rows, row_length, model_levels,             &
     &              max_filter_rows, v_begin, v_end, v_sweeps,          &
     &              global_v_filter,                                    &
     &              horizontal_level, diff_order_wind,                  &
     &              diff_coeff_wind, diff_coeff_v, L_diff_wind )
      endif !  L_pofil_new

      DEALLOCATE ( r_uv_b )
      DEALLOCATE ( r_theta_at_u )
      DEALLOCATE ( r_theta_at_v )

! ----------------------------------------------------------------------
! Section 2.3  Set u wind at poles since v has changed
!              only for global and if upper level diffusion is inactive
! ----------------------------------------------------------------------

      If ( model_domain == mt_global .and.                              &
     &     top_filt_start > model_levels ) Then

! DEPENDS ON: polar_vector_wind_n
        Call Polar_vector_wind_n(                                       &
     &                           v,                                     &
     &                           sin_theta_longitude,                   &
     &                           cos_theta_longitude, row_length,       &
     &                           n_rows, model_levels,                  &
     &                           mag_vector_np, dir_vector_np,          &
     &                           mag_vector_sp, dir_vector_sp,          &
     &                           off_x, off_y, global_row_length,       &
     &                           proc_row_group, at_extremity)

        If (at_extremity(PSouth) ) Then
          Do k = 1, model_levels
            Do i = 1, row_length
              gi = l_datastart(1) + i - 1
              u(i,1,k) = - mag_vector_sp(k) * sin ( (gi-.5)*            &
     &                     delta_lambda - dir_vector_sp(k))
            End Do
          End Do  ! k = 1, model_levels
        End If   !  at_extremity(PSouth)
        If (at_extremity(PNorth) ) Then
          Do k = 1, model_levels
            Do i = 1, row_length
              gi = l_datastart(1) + i - 1
              u(i,rows,k) = mag_vector_np(k) *                          &
     &                          sin ( (gi-.5)*delta_lambda -            &
     &                              dir_vector_np(k))
            End Do
          End Do ! k = 1, model_levels
        End If  ! at_extremity(PNorth

      endIf ! model_domain = mt_global and top_filt_start > model_levels


! u field has been changed at poles, therefore need swap_bounds call
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   u, row_length, rows, model_levels,             &
     &                   off_x, off_y, fld_type_u, .true.)

      endif ! L_pfuv .or. L_diff_wind


! ----------------------------------------------------------------------
! 3.0  Section for upper-level diffusion
! ----------------------------------------------------------------------

      If ( top_filt_start < model_levels + 1 ) then

        active_levels = top_filt_end - top_filt_start + 1
        pole_term = 8. / ( global_row_length * delta_phi )
        j_start = 1
        j_stop = rows
        If (model_domain == mt_Global) Then
          If (at_extremity(PSouth)) j_start = 2
          If (at_extremity(PNorth)) j_stop = rows - 1
        End If

! DEPENDS ON: diffupper
        call diffupper(                                                 &
     &                 theta, fld_type_p, 0,                            &
     &                 model_levels, active_levels,                     &
     &                 rows, n_rows, row_length,                        &
     &                 off_x, off_y, off_x, off_y,                      &
     &                 sec_theta_latitude, cos_v_latitude,              &
     &                 j_start, j_stop, pole_term,                      &
     &                 model_domain, at_extremity, proc_row_group,      &
     &                 top_filt_start, up_diff, max_upd_levels )

! DEPENDS ON: diffupper
        call diffupper(                                                 &
     &                 w(1-off_x, 1-off_y, 1), fld_type_p, 0,           &
     &                 model_levels - 1, active_levels - 1,             &
     &                 rows, n_rows, row_length,                        &
     &                 off_x, off_y, off_x, off_y,                      &
     &                 sec_theta_latitude, cos_v_latitude,              &
     &                 j_start, j_stop, pole_term,                      &
     &                 model_domain, at_extremity, proc_row_group,      &
     &                 top_filt_start, up_diff, max_upd_levels )

! DEPENDS ON: diffupper
        call diffupper(                                                 &
     &                 u, fld_type_u, 0,                                &
     &                 model_levels, active_levels,                     &
     &                 rows, n_rows, row_length,                        &
     &                 off_x, off_y, off_x, off_y,                      &
     &                 sec_theta_latitude, cos_v_latitude,              &
     &                 j_start, j_stop, pole_term,                      &
     &                 model_domain, at_extremity, proc_row_group,      &
     &                 top_filt_start, up_diff, max_upd_levels )

        j_start = 1
        j_stop = n_rows

! DEPENDS ON: diffupper
        call diffupper(                                                 &
     &                 v, fld_type_v, 1,                                &
     &                 model_levels, active_levels,                     &
     &                 n_rows, rows, row_length,                        &
     &                 off_x, off_y, off_x, off_y,                      &
     &                 sec_v_latitude, cos_theta_latitude,              &
     &                 j_start, j_stop, pole_term,                      &
     &                 model_domain, at_extremity, proc_row_group,      &
     &                 top_filt_start, up_diff, max_upd_levels )

! ----------------------------------------------------------------------
! Section 3.1  Set u wind at poles since v has changed
! ----------------------------------------------------------------------

        If (model_domain == mt_Global) Then

! DEPENDS ON: polar_vector_wind_n
          Call Polar_vector_wind_n(                                     &
     &                             v,                                   &
     &                             sin_theta_longitude,                 &
     &                             cos_theta_longitude, row_length,     &
     &                             n_rows, model_levels,                &
     &                             mag_vector_np, dir_vector_np,        &
     &                             mag_vector_sp, dir_vector_sp,        &
     &                             off_x, off_y, global_row_length,     &
     &                             proc_row_group, at_extremity)

          If (at_extremity(PSouth) ) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                gi = l_datastart(1) + i - 1
                u(i,1,k) = - mag_vector_sp(k) * sin ( (gi-.5)*          &
     &                       delta_lambda - dir_vector_sp(k))
              End Do
            End Do  ! k = 1, model_levels
          End If   !  at_extremity(PSouth)
          If (at_extremity(PNorth) ) Then
            Do k = 1, model_levels
              Do i = 1, row_length
                gi = l_datastart(1) + i - 1
                u(i,rows,k) = mag_vector_np(k) *                        &
     &                            sin ( (gi-.5)*delta_lambda -          &
     &                                dir_vector_np(k))
              End Do
            End Do ! k = 1, model_levels
          End If  ! at_extremity(PNorth

        endIf !model_domain == mt_Global

! DEPENDS ON: swap_bounds
        call Swap_Bounds(                                               &
     &                   theta, row_length, rows, model_levels,         &
     &                   off_x, off_y, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
        call Swap_Bounds(                                               &
     &                   w, row_length, rows, model_levels + 1,         &
     &                   off_x, off_y, fld_type_p, .false.)

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   u, row_length, rows, model_levels,             &
     &                   off_x, off_y, fld_type_u, .true.)

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &                   v, row_length, n_rows, model_levels,           &
     &                   off_x, off_y, fld_type_v, .true.)

      endif !  top_filt_start < model_levels + 1

!*******************  increments for STASH   **************************

      IF(sf(0,sect)) THEN

! T increment
        item = 381 
        IF(sf(item,sect)) THEN

          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                T_inc(i,j,k) = (theta(i,j,k) - T_inc(i,j,k))            &
     &                                  *exner_theta_levels(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        T_inc,                                                    &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

! u wind increment
        item=385
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                u_inc(i,j,k) = u(i,j,k) - u_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        u_inc,                                                    &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

! v wind increment
        item=386
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels
            DO j=1,n_rows
              DO i=1,row_length
                v_inc(i,j,k) = v(i,j,k) - v_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        v_inc,                                                    &
     &        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

! w wind increment
        item=387
        IF( sf(item,sect) ) THEN

          DO k=0,model_levels
            DO j=1,rows
              DO i=1,row_length
                w_inc(i,j,k) = w(i,j,k) - w_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        w_inc,                                                    &
     &        row_length,rows,model_levels+1,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

! exner increment
        item=388
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels+1
            DO j=1,rows
              DO i=1,row_length
                exner_inc(i,j,k) = exner(i,j,k) - exner_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        exner_inc,                                                &
     &        row_length,rows,model_levels+1,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

      ENDIF ! sf(0,sect) 

      return

      END SUBROUTINE NI_filter_Ctl

#endif
