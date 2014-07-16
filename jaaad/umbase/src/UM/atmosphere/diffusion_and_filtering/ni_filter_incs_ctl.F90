#if defined(A13_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine NI_filter_incs_Ctl

      Subroutine NI_filter_incs_Ctl(                                    &
     &                      theta, theta_star, R_u, R_v, R_w,           &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      r_theta_levels, r_rho_levels,               &
     &                      r_at_u, r_at_v,                             &
     &                      delta_lambda, delta_phi,                    &
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
     &                      horizontal_level, me,                       &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      n_proc, n_procx, n_procy, l_datastart,      &
     &                      neighbour, at_extremity, model_domain,      &
     &                      proc_row_group, proc_col_group,             &
     &                      L_polar_filter_incs, L_filter_incs,         &
     &                      L_pfcomb, L_pftheta, L_pfuv, L_pfw,         &
     &                      L_pofil_new, L_diff_incs, Ltimer,           &
     &                      exner_theta_levels,                         &
#include "argsts.h"
     &                      STASHwork)

! Purpose:
!          Call original polar filter code or new code
!          which also does horizontal diffusion as a filter
!
! Original Programmer: T. Davies
! Current code owner: T. Davies
!
! History:
! Version    Date      Comment
! ----     -------     -------
!  6.1     04/08/04    Code Introduced                     Terry Davies
!  6.2     25/12/05    Changes to filtering               Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

!include "parvars.h"
#include "csubmodl.h"
#include "typsts.h"

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  L_filter_incs                                                   &
                          ! switch for filtering of increments
     &, Ltimer                                                          &
     &, L_polar_filter_incs                                             &
                             ! switch for polar filtering of increments
     &, L_pofil_new                                                     &
                          ! new filter/diffusion switch
     &, L_pfcomb                                                        &
                          ! combined polar filter/diffusion active
     &, L_pftheta                                                       &
                          ! switch for polar filter for theta
     &, L_pfw                                                           &
                          ! switch for polar filter for w
     &, L_pfuv                                                          &
                          ! switch for polar filter for horizontal winds
     &, L_diff_incs       ! switch for horiz. diffusion of increments

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels.
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
     &, cos_v_latitude(1-off_x:row_length+off_x, 1-off_y:rows+off_y )   &
     &, sec_v_latitude(1-off_x:row_length+off_x, 1-off_y:rows+off_y )   &
!   diffusion coefficient for u/th rows
     &, diff_coeff_u(1-off_x:row_length+off_x,1-off_y:rows+off_y )      &
!   diffusion coefficient for v rows
     &, diff_coeff_v(1-off_x:row_length+off_x,1-off_y:n_rows+off_y )

      Real                                                              &
           ! limits for filtering (in radians)
     &  north_lat_limit                                                 &
                         ! latitude above which 1-2-1 filter is applied
     &, south_lat_limit                                                 &
                         ! latitude below which 1-2-1 filter is applied
     &, polar_filter_coefficient                                        &
                                     !   filter coefficient
     &, step_per_sweep                                                  &
                         ! amount in radians to increment start latitude
                         ! by per sweep
     &, polar_start_lat_limit                                           &
                              ! max latitude at which filter can start
     &, diff_coeff_wind                                                 &
                           ! NS diffusion coeff for winds
     &, diff_coeff_thermo                                               &
                           ! NS diffusion coeff for thermo
     &, diff_coeff_phi    ! NS diffusion coeff


      Integer                                                           &
     &  polar_filter_n_sweeps                                           &
                                     ! number of sweeps of filter to do
     &, max_filter_rows                                                 &
                           ! array size for u_begin etc
     &, global_u_filter                                                 &
                         ! sweep control for 1-2-1 filter
     &, global_v_filter                                                 &
                         ! sweep control for 1-2-1 filter
     &, u_begin(0:max_filter_rows)                                      &
                                    ! row pointers for 1-2-1 filter
     &, u_end(0:max_filter_rows)                                        &
                                    ! row pointers for 1-2-1 filter
     &, v_begin(0:max_filter_rows)                                      &
                                    ! row pointers for 1-2-1 filter
     &, v_end(0:max_filter_rows)                                        &
                                    ! row pointers for 1-2-1 filter
     &, u_sweeps(max_filter_rows)                                       &
                                   ! sweeps for 1-2-1 filter
     &, v_sweeps(max_filter_rows)                                       &
                                   ! sweeps for 1-2-1 filter
     &, horizontal_level                                                &
                               ! steep slope control
     &, diff_order_thermo                                               &
                           ! diffusion order for theta/w
     &, diff_order_wind                                                 &
                           ! diffusion order for u,v
     &, tar_horizontal                                                  &
                            ! steep slope control
     &, test_level                                                      &
                            ! test start level for targeted diffusion
     &, start_level                                                     &
                            ! start level for targeted diffusion
     &, end_level           ! end level for targeted diffusion

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

      Real :: exner_theta_levels                                        &
     &    (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)

      Integer                                                           &
     &  n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, proc_row_group                                                  &
     &, proc_col_group                                                  &
     &, neighbour(4)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN/OUT.
      Real                                                              &
     &  R_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)                                              &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &       model_levels)                                              &
     &, R_w(row_length, rows, model_levels)                             &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
     &, theta_star(1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &             model_levels)

!   Array arguments with intent(out):
      REAL, INTENT(OUT) ::                                              &
     & STASHwork(*)   ! Output array holding diagnostic fields

#include "fldtype.h"

! local  variables
      Integer                                                           &
     &  i, j, k     ! Loop indices

! Local arrays
      Real                                                              &
     &  field( 1-off_x:row_length+off_x, 1-off_y:rows+off_y             &
     &                   , model_levels)

      Real :: T_inc (row_length, rows, model_levels)
      Real :: u_inc (row_length, rows, model_levels)
      Real :: v_inc (row_length, n_rows, model_levels)
      Real :: w_inc (row_length, rows, model_levels)
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
        IF(sf(481,sect))t_inc(:,:,:)=theta_star(1:row_length, 1:rows, :)
        IF(sf(485,sect))u_inc(:,:,:)=R_u(1:row_length, 1:rows, :)
        IF(sf(486,sect))v_inc(:,:,:)=R_v(1:row_length,1:n_rows,:)
        IF(sf(487,sect))w_inc(:,:,:)=R_w(1:row_length, 1:rows, :)
      ENDIF

! ----------------------------------------------------------------------
! 1.0  Section for original polar filter
! ----------------------------------------------------------------------

      If (L_polar_filter_incs) Then

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('Polar_Filter_incs',3)

! DEPENDS ON: polar_filter_incs
        Call Polar_filter_incs(                                         &
     &                      theta, theta_star, R_u, R_v, R_w,           &
     &                      row_length, rows, n_rows, model_levels,     &
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
        IF (Ltimer) CALL TIMER('Polar_Filter_incs',4)

      End If   ! L_polar_filter_incs

! ----------------------------------------------------------------------
! 2.0  Section for new polar filter
! ----------------------------------------------------------------------

      if( L_pftheta .or. L_diff_incs )then

          Do k = 1, model_levels
            Do j = 1-off_y, rows+off_y
              Do i = 1-off_x, row_length+off_x
                field(i,j,k) = theta_star(i,j,k) - theta(i,j,k)
              End Do
            End Do
          End Do !  k = 1, model_levels


! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   field, row_length, rows, model_levels,         &
     &                   off_x, off_y, fld_type_p, .false.)

        if( L_pofil_new )then

! DEPENDS ON: pfwthinc
          Call pfwthinc(                                                &
     &                  field, r_theta_levels, r_rho_levels,            &
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  sec_theta_latitude, cos_v_latitude,             &
     &                  me, n_procy, delta_lambda, delta_phi,           &
     &                  rows, n_rows, row_length, model_levels,         &
     &                  max_filter_rows, u_begin, u_end,                &
     &                  u_sweeps, global_u_filter,                      &
     &                  model_levels, horizontal_level,                 &
     &                  diff_coeff_phi, diff_coeff_u, L_diff_incs )

        else  ! original combi filter

! DEPENDS ON: pofil_wth_incs
          Call pofil_wth_incs(                                          &
     &                      field, r_theta_levels, r_rho_levels,        &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      sec_theta_latitude, cos_v_latitude,         &
     &                      me, n_procy, delta_lambda, delta_phi,       &
     &                      rows, n_rows, row_length, model_levels,     &
     &                      max_filter_rows, u_begin, u_end,            &
     &                      u_sweeps, global_u_filter,                  &
     &                      model_levels, horizontal_level,             &
     &                      diff_order_thermo, diff_coeff_thermo,       &
     &                      diff_coeff_u, L_diff_incs )

        endif ! L_pofil_new

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                theta_star(i,j,k) = theta(i,j,k) + field(i,j,k)
              End Do
            End Do
          End Do

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     theta_star, row_length, rows, model_levels,  &
     &                     off_x, off_y, fld_type_p, .false.)

      endif ! L_pftheta .or. L_diff_incs

        if(L_pfw)then

          Do k = 1, model_levels - 1
            Do j = 1, rows
              Do i = 1, row_length
                field(i,j,k) = R_w(i,j,k)
              End Do
            End Do
          End Do  ! k = 1, model_levels - 1

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   field, row_length, rows, model_levels - 1,     &
     &                   off_x, off_y, fld_type_p, .false.)

        if( L_pofil_new )then

! DEPENDS ON: pfwthinc
          Call pfwthinc(                                                &
     &                  field, r_theta_levels, r_rho_levels,            &
     &                  off_x, off_y, halo_i, halo_j,                   &
     &                  sec_theta_latitude, cos_v_latitude,             &
     &                  me, n_procy, delta_lambda, delta_phi,           &
     &                  rows, n_rows, row_length, model_levels,         &
     &                  max_filter_rows, u_begin, u_end,                &
     &                  u_sweeps, global_u_filter,                      &
     &                  model_levels - 1, horizontal_level,             &
     &                  diff_coeff_phi, diff_coeff_u, L_diff_incs )

        else  ! original combi filter

! DEPENDS ON: pofil_wth_incs
          Call pofil_wth_incs(                                          &
     &                      field, r_theta_levels, r_rho_levels,        &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      sec_theta_latitude, cos_v_latitude,         &
     &                      me, n_procy, delta_lambda, delta_phi,       &
     &                      rows, n_rows, row_length, model_levels,     &
     &                      max_filter_rows, u_begin, u_end,            &
     &                      u_sweeps, global_u_filter,                  &
     &                      model_levels-1, horizontal_level,           &
     &                      diff_order_thermo, diff_coeff_thermo,       &
     &                      diff_coeff_u, L_diff_incs )

        endif ! L_pofil_new

          Do k = 1, model_levels - 1
            Do j = 1, rows
              Do i = 1, row_length
                R_w(i,j,k) = field(i,j,k)
              End Do
            End Do
          End Do  ! k = 1, model_levels - 1

        endif ! L_pfw

      if( L_pfuv .or. L_diff_incs) then

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   R_u, row_length, rows, model_levels,           &
     &                   off_x, off_y, fld_type_u, .true.)

        if( L_pofil_new )then

! DEPENDS ON: pfuinc
          Call pfuinc(                                                  &
     &                R_u, r_at_u,                                      &
     &                r_theta_levels, r_rho_levels,                     &
     &                off_x, off_y, halo_i, halo_j,                     &
     &                sec_theta_latitude, cos_v_latitude,               &
     &                me, n_procy, delta_lambda, delta_phi,             &
     &                rows, n_rows, row_length, model_levels,           &
     &                max_filter_rows, u_begin, u_end, u_sweeps,        &
     &                global_u_filter, horizontal_level,                &
     &                diff_coeff_phi, diff_coeff_u, L_diff_incs )

        else  ! original combi filter

! DEPENDS ON: pofil_u_incs
          Call pofil_u_incs(                                            &
     &                 R_u, r_at_u,                                     &
     &                 r_theta_levels, r_rho_levels,                    &
     &                 off_x, off_y, halo_i, halo_j,                    &
     &                 sec_theta_latitude, cos_v_latitude,              &
     &                 me, n_procy, delta_lambda, delta_phi,            &
     &                 rows, n_rows, row_length, model_levels,          &
     &                 max_filter_rows, u_begin, u_end, u_sweeps,       &
     &                 global_u_filter,                                 &
     &                 horizontal_level, diff_order_wind,               &
     &                 diff_coeff_wind, diff_coeff_u, L_diff_incs )

        endif ! L_pofil_new

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   R_v, row_length, n_rows, model_levels,         &
     &                   off_x, off_y, fld_type_v, .true.)

        if( L_pofil_new )then

! DEPENDS ON: pfvinc
          Call pfvinc(                                                  &
     &                 R_v, r_at_v, r_theta_levels, r_rho_levels,       &
     &                 off_x, off_y, halo_i, halo_j,                    &
     &                 sec_v_latitude, cos_theta_latitude,              &
     &                 me, n_procy, at_extremity,                       &
     &                 delta_lambda, delta_phi,                         &
     &                 rows, n_rows, row_length, model_levels,          &
     &                 max_filter_rows, v_begin, v_end, v_sweeps,       &
     &                 global_v_filter,                                 &
     &                 horizontal_level, model_domain,                  &
     &                 diff_coeff_phi, diff_coeff_v, L_diff_incs )

        else  ! original combi filter

! DEPENDS ON: pofil_v_incs
          Call pofil_v_incs(                                            &
     &                 R_v, r_at_v, r_theta_levels, r_rho_levels,       &
     &                 off_x, off_y, halo_i, halo_j,                    &
     &                 sec_v_latitude, cos_theta_latitude,              &
     &                 me, n_procy, delta_lambda, delta_phi,            &
     &                 rows, n_rows, row_length, model_levels,          &
     &                 max_filter_rows, v_begin, v_end, v_sweeps,       &
     &                 global_v_filter,                                 &
     &                 horizontal_level, diff_order_wind,               &
     &                 diff_coeff_wind, diff_coeff_v, L_diff_incs )

        endif ! L_pofil_new

      endif !  L_pfuv .or. L_diff_incs

!*******************  increments for STASH   **************************

      IF (sf(0,sect)) THEN

! T increment
      item = 481
        IF(sf(item,sect)) THEN

          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                T_inc(i,j,k) = (theta_star(i,j,k) - T_inc(i,j,k))       &
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
        item=485
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                u_inc(i,j,k) = R_u(i,j,k) - u_inc(i,j,k)
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
        item=486
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels
            DO j=1,n_rows
              DO i=1,row_length
                v_inc(i,j,k) = R_v(i,j,k) - v_inc(i,j,k)
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
        item=487
        IF( sf(item,sect) ) THEN

          DO k=1,model_levels
            DO j=1,rows
              DO i=1,row_length
                w_inc(i,j,k) = R_w(i,j,k) - w_inc(i,j,k)
              ENDDO  ! i
            ENDDO  ! j
          ENDDO  ! k

! DEPENDS ON: copydiag_3d
          CALL copydiag_3d(STASHwork(si(item,sect,im_index)),           &
     &        w_inc,                                                    &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

        ENDIF ! sf(item,sect)

      ENDIF ! sf(0,sect)

      return

      END SUBROUTINE NI_filter_incs_Ctl

#endif
