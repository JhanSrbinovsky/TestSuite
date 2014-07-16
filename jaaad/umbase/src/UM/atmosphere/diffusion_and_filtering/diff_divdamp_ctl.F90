#if defined(A13_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Diff_Divdamp_Ctl
          Subroutine Diff_Divdamp_Ctl(                                  &
     &                     L_diffusion, L_cdiffusion,L_subfilter_horiz, &
     &                     L_vertical_diffusion, L_divdamp,             &
     &                     L_ramp, ramp_lat_radians,                    &
     &                     L_Backwards, Ltimer,                         &
     &                     timestep, pos_timestep, neg_timestep,        &
     &                     theta, w, moist,                             &
     &                     u, v, rho, Exner, w_adv,                     &
     &                     r_theta_levels, r_rho_levels, r_at_u, r_at_v,&
     &                    eta_theta_levels, eta_rho_levels,             &
     &                    sec_theta_latitude,                           &
     &                    cos_theta_latitude, sec_v_latitude,           &
     &                     FV_sec_theta_latitude, cos_v_latitude,       &
     &                     sin_theta_latitude, sin_v_latitude, pi,      &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity, gc_proc_row_group,             &
     &                     mype, nproc, nproc_x, nproc_y, neighbour,    &
     &                     delta_lambda, delta_phi, L_regular,          &
     &                     lambda_p, phi_p, lambda_u, phi_v,            &
     &                     rows, row_length, n_rows,                    &
     &                     model_levels, wet_levels, model_domain,      &
     &                     global_row_length, global_rows,              &
     &                     diffusion_coefficient_thermo,                &
     &                     diffusion_coefficient_w,                     &
     &                     diffusion_coefficient_q,                     &
     &                     diffusion_coefficient_wind,                  &
     &                     diffusion_order_thermo,                      &
     &                     diffusion_order_w,                           &
     &                     diffusion_order_q,                           &
     &                     diffusion_order_wind,                        &
     &                     visc_m, visc_h,                              &
     &                     horizontal_level, tar_horizontal,            &
     &                     level_start_wind, level_stop_wind,           &
     &                     level_start_q, level_stop_q,                 &
     &                     level_start_theta, level_stop_theta,         &
     &                     L_tardiff_q, w_conv_limit, tardiffq_factor,  &
     &                     tardiffq_test, tardiffq_start, tardiffq_end, &
     &                     L_diag_w, w_local_mask,                      &
     &                     L_adjust_theta,                              &
     &                     adjust_theta_start, adjust_theta_end,        &
     &                     L_vdiff_uv, vdiffuv_test, epsilon,           &
     &                     vdiffuv_factor, vdiffuv_start, vdiffuv_end,  &
     &                     vert_diffusion_coeff_wind,                   &
     &                     vert_diffusion_coeff_q,                      &
     &                     vert_diffusion_coeff_theta,                  &
     &                     div_damp_coefficient,                        &
     &                     theta_star, R_w, moist_star, R_u, R_v        &
     &                     )

! Purpose:
!          Subroutine to interface to diffusion code
!
! Method:
!          Is described in ;
!
!
! Original Programmer: Andrew J. Malcolm
! Current code owner:  Andrew J. Malcolm
!
! History:
! Version   Date     Comment
! -------   ----     -------
!LL   5.2  24/10/00  original deck                     Andy Malcolm
!    5.3  14/09/01  correct bug with vertical diffusion    Andy Malcolm
!  5.4  11/07/02  add latitude dependent vertical option  Andy Malcolm
!   5.5  28/02/03  add targeted diffusion of q dependent on w
!                                                       Terry Davies
!   5.5  28/02/03  add diagnostic print routine         Terry Davies
!   5.5  28/02/03  add w diagnostics for targeted diffusion
!                                                       Terry Davies
! 6.0   18/08/03  Add steep slope test to targeted diffusion
!                                                         Terry Davies
!   6.1  03/08/04  remove call to w_div_diag        Terry Davies
!  6.2  25/12/05  Variable resolution changes            Yongming Tang
!   6.2  10/02/06  Add 3D subgrid turbulence scheme. Carol Halliwell
! 6.2   04/08/05  add calls to convective adjustment and
!                 new vertical diffusion of winds         Terry Davies
! 6.2   25/12/05  change name of vapour variables to moist to permit use
!               of specific quantities or mixing ratios.   Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

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
                         ! number of wet levels.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, offx                                                            &
                     ! Size of small halo in i
     &, offy                                                            &
                     ! Size of small halo in j.
     &, mype                                                            &
                     ! processor Id
     &, global_row_length                                               &
     &,  global_rows                                                    &
     &, gc_proc_row_group                                               &
     &, nproc                                                           &
                  ! Total number of processors
     &, nproc_x                                                         &
                   ! Number of processors in longitude
     &, nproc_y                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)         ! Array with the Ids of the four neighbours
                             ! in the horizontal plane

      REAL                                                              &
     & timestep                                                         &
                            ! atmosphere model timestep
     &,pos_timestep                                                     &
                    ! = +timestep.
     &,neg_timestep                                                     &
                    ! = -timestep.
     &,  w_conv_limit                                                   &
     &,  vdiffuv_test                                                   &
     &,  epsilon                                                        &
     &, pi

      Real                                                              &
           ! primary model variables
     &  u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)     &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)     &
     &, rho(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)     &
     &, Exner(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &        model_levels + 1)                                         &
     &, moist (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         wet_levels)                                              &
     &, theta (1-offx:row_length+offx, 1-offy:rows+offy,                &
     &         model_levels)                                            &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          0:model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)                                           &
     &, eta_theta_levels (0:model_levels)                               &
     &, eta_rho_levels (model_levels)                                   &
     &, cos_theta_latitude (1-offx:row_length+offx,                     &
     &                    1-offy:rows+offy)                             &
     &, sec_theta_latitude (1-offx:row_length+offx,                     &
     &                    1-offy:rows+offy)                             &
     &, sec_v_latitude (1-offx:row_length+offx, 1-offy:n_rows+offy)     &
     &, cos_v_latitude (1-offx:row_length+offx, 1-offy:n_rows+offy)     &
     &, FV_sec_theta_latitude (1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy)                        &
     &, sin_theta_latitude(row_length, rows)                            &
     &, sin_v_latitude(row_length, n_rows)

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, diffusion_order(model_levels)

      Real                                                              &
     &  diffusion_coefficient(model_levels)

      Logical                                                           &
     &  L_diffusion                                                     &
     &, L_cdiffusion                                                    &
     &, L_regular                                                       &
                        !  Variable resolution switch
     &, L_vertical_diffusion                                            &
     &, L_vdiff_uv                                                      &
     &, L_adjust_theta                                                  &
     &, L_tardiff_q                                                     &
     &, L_ramp                                                          &
     &, L_divdamp                                                       &
     &, L_Backwards                                                     &
     &, Ltimer                                                          &
     &, L_diag_w                                                        &
     &, L_subfilter_horiz   
                            ! subgrid turbulence scheme

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
            ! VarRes horizontal co-ordinate information
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, phi_p(1-halo_j:rows+halo_j)                                     &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_v(1-halo_j:n_rows+halo_j)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      INTEGER                                                           &
     & diffusion_order_thermo(model_levels)                             &
                                              ! value * del^2: order=2
     &,diffusion_order_wind(model_levels)                               &
                                              ! gives del^4 diffusion
     &,diffusion_order_w(model_levels-1)                                &
                                              !
     &,diffusion_order_q(wet_levels)          !

      REAL                                                              &
     & diffusion_coefficient_thermo(model_levels)                       &
     &,diffusion_coefficient_wind(model_levels)                         &
     &,diffusion_coefficient_w(model_levels-1)                          &
     &,diffusion_coefficient_q(wet_levels)                              &
     &, tardiffq_factor                                                 &
                              ! targeted diffusion coefficient
     &,visc_m(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j          &
     &                                             , model_levels)      &
!            ! diffusion coefficient for momentum from subgrid
!            ! turbulence scheme
     &,visc_h(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j          &
     &                                             , model_levels)      &
!            ! diffusion coefficient for heat and moisture from
!            ! subgrid turbulence scheme
     &,weight1, weight2, weight3 ! interpolation weights

      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     & diff_coeff_u                                                     &
                             ! visc_m interpolated onto
     &,diff_coeff_v                                                     &
                             ! appropriate points
     &,diff_coeff_th                                                    &
                             ! for use in the turb_diff_*
     &,diff_coeff_centre                                                &
                             ! diffusion routines.                      
     &,w_diff_coeff_u                                                   &
                             ! horizontal diffusion coefficients for    
     &,w_diff_coeff_v        ! diffusion of w

! local
       Integer                                                          &
     & i,j,k

      Integer                                                           &
     & horizontal_level                                                 &
     &, tar_horizontal

      Integer                                                           &
     & level_start_wind, level_stop_wind                                &
     &,level_start_theta, level_stop_theta                              &
     &,level_start_q, level_stop_q                                      &
     &, adjust_theta_start                                              &
     &, adjust_theta_end                                                &
     &, vdiffuv_start                                                   &
     &, vdiffuv_end                                                     &
     &, tardiffq_test                                                   &
     &, tardiffq_start                                                  &
     &, tardiffq_end

      Real                                                              &
     & vert_diffusion_coeff_wind                                        &
     &,vert_diffusion_coeff_theta                                       &
     &,vert_diffusion_coeff_q                                           &
     &, vdiffuv_factor                                                  &
     &, ramp_lat_radians

      Real                                                              &
     & div_damp_coefficient(model_levels)

!   variables with intent INOUT

      Real                                                              &
     &  R_u(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &        model_levels)                                             &
     &, R_v(1-offx:row_length+offx, 1-offy:n_rows+offy,                 &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels)                             &
     &, moist_star(1-offx:row_length+offx,                              &
     &             1-offy:rows+offy, wet_levels)                        &
     &, theta_star(1-offx:row_length+offx,                              &
     &               1-offy:rows+offy, model_levels)                    &
     &, w_local_mask(row_length,rows)

#include "domtyp.h"
! Local variables
      Integer                                                           &
     & levels

! ----------------------------------------------------------------------
! Section 1.0  Horizontal Diffusion
! ----------------------------------------------------------------------

      If ( L_subfilter_horiz) Then
        ALLOCATE (diff_coeff_u(1-offx:row_length, rows                  &
     &                                        , model_levels))
        ALLOCATE (diff_coeff_v(row_length, 1-offy:n_rows+offy           &
     &                                        , model_levels))
        ALLOCATE (diff_coeff_th(1-offx:row_length+offx, 1-offy:rows+offy&
     &                                        , model_levels))
        ALLOCATE (diff_coeff_centre(1-offx:row_length, 1-offy:rows,     &
     &                                          model_levels))
        ALLOCATE (w_diff_coeff_u(1-offx:row_length, rows                &
     &                                        , model_levels))
        ALLOCATE (w_diff_coeff_v(row_length, 1-offy:n_rows+offy         &
     &                                        , model_levels))  
      Else
        ALLOCATE (diff_coeff_u(1,1,1))
        ALLOCATE (diff_coeff_v(1,1,1))
        ALLOCATE (diff_coeff_th(1,1,1))
        ALLOCATE (diff_coeff_centre(1,1,1))
        ALLOCATE (w_diff_coeff_u(1,1,1))
        ALLOCATE (w_diff_coeff_v(1,1,1))
      End If

      If ( L_diffusion ) Then

!       Run with a positive timestep if integrating backwards.
        IF (L_Backwards) timestep = pos_timestep

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_diff_theta',3)
! DEPENDS ON: h_diff_theta
        Call h_diff_theta(                                              &
     &                     theta,                                       &
     &                     r_theta_levels,                              &
     &                     offx, offy, halo_i, halo_j, offx, offy,      &
     &                     at_extremity, gc_proc_row_group,             &
     &                     nproc, nproc_x, nproc_y, neighbour,          &
     &                     delta_lambda, delta_phi,                     &
     &                     timestep, rows, row_length,                  &
     &                     model_levels, model_domain,                  &
     &                     global_row_length,                           &
     &                     diffusion_coefficient_thermo,                &
     &                     diffusion_order_thermo, theta_star,          &
     &                     horizontal_level)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_diff_theta',4)

! diffusion of w not called, code left in case needs to be reinstated
!        IF (Ltimer) CALL TIMER('h_diff_w',3)
!        Call h_diff_theta(
!     &                     w(1-offx, 1-offy, 1),
!     &                     r_theta_levels,
!     &                     offx, offy, halo_i, halo_j, 0, 0,
!     &                     at_extremity, gc_proc_row_group,
!     &                     nproc, nproc_x, nproc_y, neighbour,
!     &                     delta_lambda, delta_phi,
!     &                     timestep, rows, row_length,
!     &                     model_levels-1, model_domain,
!     &                     global_row_length,
!     &                     diffusion_coefficient_w,
!     &                     diffusion_order_w, R_w,
!     &                     horizontal_level)
!        IF (Ltimer) CALL TIMER('h_diff_w',4)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_diff_q',3)
! DEPENDS ON: h_diff_q
        Call h_diff_q(                                                  &
     &                     moist,                                       &
     &                     r_theta_levels,                              &
     &                     offx, offy, halo_i, halo_j, offx, offy,      &
     &                     at_extremity, gc_proc_row_group,             &
     &                     nproc, nproc_x, nproc_y, neighbour,          &
     &                     delta_lambda, delta_phi,                     &
     &                     timestep, rows, row_length,                  &
     &                     model_levels,                                &
     &                     wet_levels, model_domain,                    &
     &                     global_row_length,                           &
     &                     diffusion_coefficient_q,                     &
     &                     diffusion_order_q, moist_star,               &
     &                     horizontal_level)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_diff_q',4)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_diff_u',3)
! DEPENDS ON: h_diff_u
        Call h_diff_u(                                                  &
     &                 u, r_at_u, r_theta_levels,                       &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 at_extremity, gc_proc_row_group,                 &
     &                 nproc, nproc_x, nproc_y, neighbour,              &
     &                 delta_lambda, delta_phi,                         &
     &                 timestep, rows, row_length,                      &
     &                 model_levels, model_domain,                      &
     &                 diffusion_coefficient_wind,                      &
     &                 diffusion_order_wind, R_u,                       &
     &                 horizontal_level)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_diff_u',4)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_diff_v',3)
! DEPENDS ON: h_diff_v
        Call h_diff_v(                                                  &
     &                 v, r_at_v, r_theta_levels, rows,                 &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 at_extremity, gc_proc_row_group,                 &
     &                 nproc, nproc_x, nproc_y, neighbour,              &
     &                 delta_lambda, delta_phi,                         &
     &                 timestep, n_rows, row_length,                    &
     &                 model_levels, model_domain,                      &
     &                 diffusion_coefficient_wind,                      &
     &                 diffusion_order_wind, R_v,                       &
     &                 horizontal_level)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_diff_v',4)

!       Go back to negative timestep if integrating backwards.
        IF (L_Backwards) timestep = neg_timestep

      End If        !  L_diffusion
! If L_cdiffusion .true. then call conservative diffusion
      If ( L_cdiffusion ) Then

!       Run with a positive timestep if integrating backwards.
        IF (L_Backwards) timestep = pos_timestep

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_cdiff_theta',3)
! DEPENDS ON: h_cdiff_theta
        Call h_cdiff_theta(                                             &
     &                     theta,                                       &
     &                     r_theta_levels, r_rho_levels,                &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     sec_theta_latitude, cos_v_latitude,          &
     &                     offx, offy, halo_i, halo_j, offx, offy,      &
     &                     at_extremity, gc_proc_row_group,             &
     &                     nproc, nproc_x, nproc_y, neighbour,          &
     &                     delta_lambda, delta_phi,                     &
     &                     timestep, rows, n_rows, row_length,          &
     &                     model_levels, model_domain,                  &
     &                     global_row_length,                           &
     &                     diffusion_coefficient_thermo,                &
     &                     diffusion_order_thermo, theta_star,          &
     &                     horizontal_level)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_cdiff_theta',4)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_cdiff_q',3)
! DEPENDS ON: h_cdiff_q
        Call h_cdiff_q(                                                 &
     &                     moist,                                       &
     &                     r_theta_levels, r_rho_levels,                &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     sec_theta_latitude, cos_v_latitude,          &
     &                     offx, offy, halo_i, halo_j, offx, offy,      &
     &                     at_extremity, gc_proc_row_group,             &
     &                     nproc, nproc_x, nproc_y, neighbour,          &
     &                     delta_lambda, delta_phi,                     &
     &                     timestep, rows, n_rows, row_length,          &
     &                     model_levels,                                &
     &                     wet_levels, model_domain,                    &
     &                     global_row_length,                           &
     &                     diffusion_coefficient_q,                     &
     &                     diffusion_order_q, moist_star,               &
     &                     horizontal_level)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_cdiff_q',4)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_cdiff_u',3)
! DEPENDS ON: h_cdiff_u
        Call h_cdiff_u(                                                 &
     &                 u, r_at_u, r_theta_levels,                       &
     &                 eta_theta_levels, eta_rho_levels,                &
     &                 sec_theta_latitude, cos_v_latitude,              &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 at_extremity, gc_proc_row_group,                 &
     &                 nproc, nproc_x, nproc_y, neighbour,              &
     &                 delta_lambda, delta_phi,                         &
     &                 timestep, rows, n_rows, row_length,              &
     &                 model_levels, model_domain,                      &
     &                 diffusion_coefficient_wind,                      &
     &                 diffusion_order_wind, R_u,                       &
     &                 horizontal_level)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_cdiff_u',4)

! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_cdiff_v',3)
! DEPENDS ON: h_cdiff_v
        Call h_cdiff_v(                                                 &
     &                 v, r_at_v, r_theta_levels, rows,                 &
     &                 eta_theta_levels, eta_rho_levels,                &
     &                 cos_theta_latitude, sec_v_latitude,              &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 at_extremity, gc_proc_row_group,                 &
     &                 nproc, nproc_x, nproc_y, neighbour,              &
     &                 delta_lambda, delta_phi,                         &
     &                 timestep, n_rows, row_length,                    &
     &                 model_levels, model_domain,                      &
     &                 diffusion_coefficient_wind,                      &
     &                 diffusion_order_wind, R_v,                       &
     &                 horizontal_level)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('h_cdiff_v',4)

!       Go back to negative timestep if integrating backwards.
        IF (L_Backwards) timestep = neg_timestep

      End If        !  L_cdiffusion


      if(L_adjust_theta) then
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('conv_adjust_theta',3)

        levels = adjust_theta_end - adjust_theta_start + 1
! DEPENDS ON: conv_adjust_theta
        Call conv_adjust_theta                                          &
     &                     (theta_star, Exner,                          &
     &                      offx, offy,                                 &
     &                      rows, row_length, model_levels,             &
     &                      adjust_theta_start, adjust_theta_end )
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('conv_adjust_theta',4)
      endif ! L_adjust_theta

        if(L_vdiff_uv) then
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('vert_diff_uv',3)

       levels = vdiffuv_end - vdiffuv_start + 1
! DEPENDS ON: vert_diff_uv
       Call vert_diff_uv(                                               &
     &                    u, v, w,                                      &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    r_at_u, r_at_v,                               &
     &                    offx, offy, halo_i, halo_j,                   &
     &                    mype, nproc, model_domain, at_extremity,      &
     &                    timestep, rows, n_rows, row_length,           &
     &                    model_levels, levels,                         &
     &                    vdiffuv_start, vdiffuv_end,                   &
     &                    vdiffuv_factor, vdiffuv_test,                 &
     &                    R_u, R_v)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('vert_diff_uv',4)
        endif ! L_vdiff_uv

! Test to apply targetted diffusion on q if w exceeds prescribed limit
        IF( L_tardiff_q) then
          if( tar_horizontal > 0 ) then
! DEPENDS ON: timer
            IF (Ltimer) CALL TIMER('tardiff_q_wss',3)
! DEPENDS ON: tardiff_q_wss
            Call tardiff_q_wss(                                         &
     &                         moist, w_adv,                            &
     &                         r_theta_levels, r_rho_levels,            &
     &                         sec_theta_latitude,                      &
     &                         cos_theta_latitude, cos_v_latitude,      &
     &                         offx, offy, halo_i, halo_j, offx, offy,  &
     &                         at_extremity, gc_proc_row_group,         &
     &                         delta_lambda, delta_phi,                 &
     &                         timestep, rows, n_rows, row_length,      &
     &                         model_levels, wet_levels,                &
     &                         model_domain, global_row_length,         &
     &                         moist_star, w_conv_limit, tar_horizontal,&
     &                         tardiffq_factor, tardiffq_test,          &
     &                         tardiffq_start, tardiffq_end,            &
     &                         L_diag_w, w_local_mask )
! DEPENDS ON: timer
            IF (Ltimer) CALL TIMER('tardiff_q_wss',4)
          else  !  tar_horizontal = 0
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('tardiff_q_w',3)
! DEPENDS ON: tardiff_q_w
          Call tardiff_q_w(                                             &
     &                    moist, w_adv,                                 &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_theta_latitude,                           &
     &                    cos_theta_latitude, cos_v_latitude,           &
     &                    offx, offy, halo_i, halo_j, offx, offy,       &
     &                    at_extremity, gc_proc_row_group,              &
     &                    delta_lambda, delta_phi,                      &
     &                    timestep, rows, n_rows, row_length,           &
     &                    model_levels, wet_levels,                     &
     &                    model_domain, global_row_length,              &
     &                    moist_star, w_conv_limit,                     &
     &                    tardiffq_factor, tardiffq_test,               &
     &                    tardiffq_start, tardiffq_end,                 &
     &                    L_diag_w, w_local_mask)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('tardiff_q_w',4)
          endif !  tar_horizontal > 0
        ENDIF !   L_tardiff_q


       If ( L_subfilter_horiz) Then
!        Run with a positive timestep if integrating backwards.
         If (L_Backwards) timestep = pos_timestep
!
!Interpolate visc_m onto appropriate points for the diffusion routines.
! At this point the diffusion coefficients are on w (theta) points.
!
         Do k = 2, model_levels

           Do j = 1, rows
             Do i = 1-offx,row_length
               diff_coeff_u(i,j,k) = 0.5*(                              &
     &                   visc_h(i+1,j,k) + visc_h(i,j,k) )
               w_diff_coeff_u(i,j,k) = 0.5*(                            &
     &                   visc_m(i+1,j,k) + visc_m(i,j,k) )     
             End Do
           End Do

           Do j = 1-offy, n_rows+offy
             Do i = 1, row_length
               diff_coeff_v(i,j,k) = 0.5*(                              &
     &                   visc_h(i,j+1,k) + visc_h(i,j,k) )
               w_diff_coeff_v(i,j,k) = 0.5*(                            &
     &                   visc_m(i,j+1,k) + visc_m(i,j,k) )  
             End Do
           End Do

           Do j = 1-offy, rows + offy
             Do i = 1-offx, row_length+offx
               weight1 =                                                &
     &           r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
               weight2 =                                                &
     &           r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
               weight3 =                                                &
     &           r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
               diff_coeff_th(i,j,k) = ( weight3*visc_m(i,j,k) +         &
     &               weight2*visc_m(i,j,k-1) )/weight1
             End Do
           End Do

           Do j = 1-offy, rows
             Do i = 1-offx, row_length
               diff_coeff_centre(i,j,k) = 0.25*(                         &
     &             diff_coeff_th(i,j,k) + diff_coeff_th(i+1,j,k) +       &
     &             diff_coeff_th(i,j+1,k) + diff_coeff_th(i+1,j+1,k)  )
             End Do
           End Do

         End Do

         Do j = 1, rows
           Do i = 1-offx, row_length
             diff_coeff_u(i,j,1) = 0.5*(                                &
     &                 visc_h(i+1,j,1) + visc_h(i,j,1) )
             w_diff_coeff_u(i,j,1) = 0.5*(                              &
     &                 visc_m(i+1,j,1) + visc_m(i,j,1) )
           End Do
         End Do

         Do j = 1-offy, n_rows+offy
           Do i = 1, row_length
             diff_coeff_v(i,j,1) = 0.5*(                                &
     &                   visc_h(i,j+1,1) + visc_h(i,j,1) )
             w_diff_coeff_v(i,j,1) = 0.5*(                              &
     &                   visc_m(i,j+1,1) + visc_m(i,j,1) ) 
           End Do
         End Do

         Do j = 1-offy, rows + offy
           Do i = 1-offx, row_length+offx
             diff_coeff_th(i,j,1) = visc_m(i,j,1)
           End Do
         End Do

         Do j = 1-offy, rows
           Do i = 1-offx, row_length
             diff_coeff_centre(i,j,1) = 0.25* (                         &
     &         visc_m(i,j,1) + visc_m(i,j+1,1) +                        &
     &         visc_m(i+1,j,1) + visc_m(i+1,j+1,1)  )
           Enddo
         Enddo
!
!  Call diffusion routines with (3D) coefficients from the subgrid
!  turbulence scheme.
!
! DEPENDS ON: timer
         IF (Ltimer) CALL TIMER('turb_diff_th',3)
! DEPENDS ON: turb_diff_th
         Call turb_diff_th(                                             &
     &                    theta,                                        &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_theta_latitude,                           &
     &                    cos_v_latitude, sec_v_latitude,               &
     &                    offx, offy, halo_i, halo_j,                   &
     &                    delta_lambda, delta_phi, timestep,            &
     &                    rows, n_rows, row_length, model_levels,       &
     &                    model_levels, diff_coeff_u, diff_coeff_v,     &
     &                    theta_star)
! DEPENDS ON: timer
         IF (Ltimer) CALL TIMER('turb_diff_th',4)

! DEPENDS ON: timer
         IF (Ltimer) CALL TIMER('turb_diff_q',3)
! DEPENDS ON: turb_diff_q
         Call turb_diff_q(                                              &
     &                    moist,                                        &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_theta_latitude,                           &
     &                    cos_v_latitude, sec_v_latitude,               &
     &                    offx, offy, halo_i, halo_j,                   &
     &                    delta_lambda, delta_phi, timestep,            &
     &                    rows, n_rows, row_length,                     &
     &                    model_levels, wet_levels, wet_levels,         &
     &                    diff_coeff_u, diff_coeff_v,                   &
     &                    moist_star)
! DEPENDS ON: timer
         IF (Ltimer) CALL TIMER('turb_diff_q',4)

! DEPENDS ON: timer
         IF (Ltimer) CALL TIMER('turb_diff_u',3)
! DEPENDS ON: turb_diff_u
         Call turb_diff_u(                                              &
     &                    u, r_at_u,                                    &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_theta_latitude,                           &
     &                    cos_v_latitude, sec_v_latitude,               &
     &                    offx, offy, halo_i, halo_j,                   &
     &                    delta_lambda, delta_phi, timestep,            &
     &                    rows, n_rows, row_length,                     &
     &                    model_levels, model_levels,                   &
     &                    diff_coeff_th, diff_coeff_centre,             &
     &                    R_u )
! DEPENDS ON: timer
         IF (Ltimer) CALL TIMER('turb_diff_u',4)

! DEPENDS ON: timer
         IF (Ltimer) CALL TIMER('turb_diff_v',3)
! DEPENDS ON: turb_diff_v
         Call turb_diff_v(                                              &
     &                    v, r_at_v,                                    &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_v_latitude,                               &
     &                    cos_theta_latitude, sec_theta_latitude,       &
     &                    offx, offy, halo_i, halo_j,                   &
     &                    delta_lambda, delta_phi, timestep,            &
     &                    rows, n_rows, row_length,                     &
     &                    model_levels, model_levels,                   &
     &                    diff_coeff_th, diff_coeff_centre,             &
     &                    R_v )
! DEPENDS ON: timer
         IF (Ltimer) CALL TIMER('turb_diff_v',4)

! DEPENDS ON: timer
         IF (Ltimer) CALL TIMER('turb_diff_w',3)
! DEPENDS ON: turb_diff_w
         Call turb_diff_w(                                              &
     &                   w,                                             &
     &                   r_theta_levels, r_rho_levels,                  &
     &                   sec_theta_latitude, cos_v_latitude,            &
     &                   sec_v_latitude,                                &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   delta_lambda, delta_phi,                       &
     &                   timestep, rows, n_rows, row_length,            &
     &                   model_levels, model_levels - 1,                &
     &                   w_diff_coeff_u, w_diff_coeff_v,                &
     &                   R_w )
! DEPENDS ON: timer
         IF (Ltimer) CALL TIMER('turb_diff_w',4)

!        Go back to negative timestep if integrating backwards.
         IF (L_Backwards) timestep = neg_timestep
       Endif   !L_subfilter_horiz

! ----------------------------------------------------------------------
! Section 2.0  Vertical Diffusion
! ----------------------------------------------------------------------

      IF (L_vertical_diffusion) then

!       Run with a positive timestep if integrating backwards.
        IF (L_Backwards) timestep = pos_timestep

        IF (vert_diffusion_coeff_wind  >   tiny(1.0)) then
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('vert_diff_u',3)
! DEPENDS ON: vert_diff_u
          Call vert_diff_u(                                             &
     &                     u, r_theta_levels, r_at_u,                   &
     &                     sin_theta_latitude,pi,L_ramp,                &
     &                     ramp_lat_radians,                            &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity,                                &
     &                     timestep, rows, row_length,                  &
     &                     model_levels, model_domain,                  &
     &                     level_start_wind, level_stop_wind,           &
     &                     vert_diffusion_coeff_wind,                   &
     &                     R_u )
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('vert_diff_u',4)

! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('vert_diff_v',3)
! DEPENDS ON: vert_diff_v
          Call vert_diff_v(                                             &
     &                     v, r_theta_levels, r_at_v,                   &
     &                     sin_v_latitude,pi,L_ramp,                    &
     &                     ramp_lat_radians,                            &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity,                                &
     &                     timestep, rows, n_rows, row_length,          &
     &                     model_levels, model_domain,                  &
     &                     level_start_wind, level_stop_wind,           &
     &                     vert_diffusion_coeff_wind,                   &
     &                     R_v )
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('vert_diff_v',4)
        Endif  ! vert_diffusion_coeff_wind

        IF (vert_diffusion_coeff_q  >   tiny(1.0)) then
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('vert_diff_q',3)
! DEPENDS ON: vert_diff_q
          Call vert_diff_q(                                             &
     &                     moist, r_theta_levels, r_rho_levels,         &
     &                     sin_theta_latitude,pi,L_ramp,                &
     &                     ramp_lat_radians,                            &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     offx, offy,                                  &
     &                     at_extremity,                                &
     &                     timestep, rows, row_length,                  &
     &                     model_levels, wet_levels, model_domain,      &
     &                     level_start_q, level_stop_q,                 &
     &                     vert_diffusion_coeff_q,                      &
     &                     moist_star )
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('vert_diff_q',4)
        Endif  ! vert_diffusion_coeff_q

        IF (vert_diffusion_coeff_theta  >   tiny(1.0)) then
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('vert_diff_theta',3)
! DEPENDS ON: vert_diff_theta
          Call vert_diff_theta(                                         &
     &                     theta, r_theta_levels, r_rho_levels,         &
     &                     sin_theta_latitude,pi,L_ramp,                &
     &                     ramp_lat_radians,                            &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity,                                &
     &                     timestep, rows, row_length,                  &
     &                     model_levels, model_domain,                  &
     &                     level_start_theta, level_stop_theta,         &
     &                     vert_diffusion_coeff_theta,                  &
     &                     theta_star)
! DEPENDS ON: timer
          IF (Ltimer) CALL TIMER('vert_diff_theta',4)
        Endif  ! vert_diffusion_coeff_q

!       Go back to negative timestep if integrating backwards.
        IF (L_Backwards) timestep = neg_timestep

      END IF     !  (L_vertical_diffusion)

! ----------------------------------------------------------------------
! Section 3.0  Divergence damping
! ----------------------------------------------------------------------

      If (L_divdamp) Then

!       Run with a positive timestep if integrating backwards.
        IF (L_Backwards) timestep = pos_timestep

! Call divergence damping code
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('div_damp',3)
! DEPENDS ON: div_damp
        Call div_damp(                                                  &
     &                u, v, rho,                                        &
     &                r_at_u, r_at_v, r_rho_levels,                     &
     &                FV_sec_theta_latitude, cos_v_latitude,            &
     &                offx, offy, halo_i, halo_j,                       &
     &                mype, at_extremity, gc_proc_row_group,            &
     &                nproc, nproc_x, nproc_y, neighbour,               &
     &                delta_lambda, delta_phi,                          &
     &                timestep, rows, n_rows, row_length,               &
     &                model_levels, model_domain,                       &
     &                global_row_length, div_damp_coefficient,          &
     &                R_u, R_v)
! DEPENDS ON: timer
        IF (Ltimer) CALL TIMER('div_damp',4)

!       Go back to negative timestep if integrating backwards.
        IF (L_Backwards) timestep = neg_timestep

      End If         !  L_divdamp

       If (L_subfilter_horiz) then
         DEALLOCATE (diff_coeff_u)
         DEALLOCATE (diff_coeff_v)
         DEALLOCATE (diff_coeff_th)
         DEALLOCATE (diff_coeff_centre)
         DEALLOCATE (w_diff_coeff_u)
         DEALLOCATE (w_diff_coeff_v) 
       End If
!   call to w_div_diag removed from diffctl (now called from atmstep)
      Return
      END SUBROUTINE Diff_Divdamp_Ctl
#endif
