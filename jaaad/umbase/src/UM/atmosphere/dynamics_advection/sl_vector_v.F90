#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SL_Vector_v
!

      Subroutine SL_Vector_v(                                           &
     &                         v, u_adv, v_adv, w_adv,                  &
     &                         Yu, Yv, Yw,                              &
     &                         eta_theta_levels, eta_rho_levels,        &
     &                         r_theta_levels, r_rho_levels,            &
     &                         r_at_u, r_at_v,                          &
     &                         depart_lambda,depart_phi,depart_r,       &
     &                         L_regular, L_interp_depart,              &
     &                         cos_theta_latitude, sec_theta_latitude,  &
     &                         sin_theta_latitude, cos_v_latitude,      &
     &                         sec_v_latitude, sin_v_latitude,          &
     &                         tan_theta_latitude, tan_v_latitude,      &
     &                         cos_theta_longitude,                     &
     &                         sin_theta_longitude,                     &
     &                         delta_lambda, delta_phi,                 &
     &                         glambda_p, phi_p, glambda_u, phi_v,      &
     &                         gdlambda_p, dphi_p, gdlambda_u, dphi_v,  &
     &                         grecip_dlamp, recip_dphip, grecip_dlamu, &
     &                         recip_dphiv, wt_lambda_p, wt_phi_p,      &
     &                         wt_lambda_u, wt_phi_v, lambda_p_rm,      &
     &                         lambda_p_rp, lambda_u_rm, lambda_u_rp,   &
     &                         phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,  &
     &                         recip_lambda_p_m, recip_lambda_p_0,      &
     &                         recip_lambda_p_p, recip_lambda_p_p2,     &
     &                         recip_lambda_u_m, recip_lambda_u_0,      &
     &                         recip_lambda_u_p, recip_lambda_u_p2,     &
     &                         recip_phi_p_m, recip_phi_p_0,            &
     &                         recip_phi_p_p, recip_phi_p_p2,           &
     &                         recip_phi_v_m, recip_phi_v_0,            &
     &                         recip_phi_v_p, recip_phi_v_p2,           &
     &                         Base_lambda, base_phi,                   &
     &                         lambda_p_end, phi_p_end, dlambda_p_end,  &
     &                         dphi_p_end, dphi_v_end,                  &
     &                         recip_dlam, recip_dphi, max_look,        &
     &                         look_lam, look_phi, halo_lam, halo_phi,  &
     &                         n_Y_arrays, n_Yd_arrays, n_Ydw_arrays,   &
     &                         timestep, Pi,                            &
     &                         alpha_3, alpha_4, LAM_max_cfl,           &
     &                         row_length, rows, n_rows, model_levels,  &
     &                         Depart_scheme, Depart_order,             &
     &                         high_order_scheme, monotone_scheme,      &
     &                         model_domain, L_2d_sl_geometry,          &
     &                         L_high, L_mono, L_conserv,               &
     &                         L_Robert_high, L_Robert_mono,            &
     &                         Robert_high_order_scheme,                &
     &                         Robert_monotone_scheme,                  &
     &                         check_bottom_levels,                     &
     &                         interp_vertical_search_tol,              &
     &                         first_constant_r_rho_level,              &
     &                         L_trivial_trigs,                         &
     &                         me, n_proc, n_procx, n_procy,            &
     &                         off_x, off_y, halo_i, halo_j,            &
     &                         global_row_length, global_rows,          &
     &                         l_datastart, at_extremity,               &
     &                         g_i_pe, proc_row_group, proc_col_group,  &
     &                         L_sl_halo_reprod,                        &
     &                         L_new_tdisc, CycleNo, R_v, Error_code)

! Purpose:
!          Performs vector semi-Lagrangian integration of values at
!          v points given the forcing functions for the first advection
!          step.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
!          and
!
!          A semi-Implicit scheme for the Unified Model.
!          F.R. Division working paper No 154.
!          M. J. P. Cullen and T. Davies.
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Date     Version     Comment
! ----     -------     -------
!LL   5.1   14/02/00  Use DOMTYP parameters                    P.Burton
!LL  5.2 26/09/00  interpolation in eta not z    Andy Malcolm
!LL   5.2   27/09/00  cyclic LAM fix                         A.Malcolm
!   5.3     15/09/01  add mt_bi_cyclic_LAM code            A. Malcolm
! 6.1   08/06/04  Add interpolation of departure points   A. Malcolm
! 6.2      04/10/05  Changes for cycling semi-Lagrangian scheme.
!                                                        M. Diamantakis
!  6.2   25/12/05  Variable resolution changes           Yongming Tang
!    6.4 07/02/07 Resolve multiple declaration of sin_lambda_ad
!                 S.D.Mullerworth
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
     &, me                                                              &
                         ! My processor number
     &, n_proc                                                          &
                         ! Total number of processors
     &, n_procx                                                         &
                         ! Number of processors in longitude
     &, n_procy                                                         &
                         ! Number of processors in latitude
     &, halo_i                                                          &
                         ! Size of halo in i.
     &, halo_j                                                          &
                         ! Size of halo in j.
     &, off_x                                                           &
                     ! Size of small halo in i
     &, off_y                                                           &
                     ! Size of small halo in j.
     &, l_datastart(3)                                                  &
                           ! First gridpoints held by this processor.
     &, proc_row_group                                                  &
                         ! Group id for processors on the same row
     &, proc_col_group                                                  &
                         ! Group id for processors on the same column
     &, max_look                                                        &
                         ! max size of look-up arrays for searches
     &, global_row_length                                               &
                              ! global number of points on a row
     &, global_rows                                                     &
                            ! global number of points in a column
     &, g_i_pe(1-halo_i:global_row_length+halo_i)                       &
                                                  ! processor on my
                    ! processor-row holding a given value in i direction
     &, n_Y_arrays                                                      &
                         ! = 1 for global, 3 for LAM
     &, n_Yd_arrays                                                     &
                         ! = 1 for global, 3 for LAM
     &, n_Ydw_arrays                                                    &
                         ! = 1 for global, 2 for LAM
     &, LAM_max_cfl(2)                                                  &
     &, CycleNo

      Logical                                                           &
     &  L_sl_halo_reprod                                                &
                         ! if true then sl code bit repoducible with
                         ! any sensible halo size
     &, L_new_tdisc      ! if true activate time discretization
                         ! for iterative scheme
      Real                                                              &
           ! physical constants
     &  Pi

      Integer                                                           &
     &  high_order_scheme                                               &
                           ! a code saying which high order scheme to
                           ! use. 1 = tensor tri-cubic lagrange order
                           ! (j,i,k) no other options available at
                           ! present.
     &, monotone_scheme                                                 &
                        ! a code saying which monotone scheme to use.
                        ! 1 = tri-linear order (j,i,k)
                        ! no other options available at present.
     &, Depart_scheme                                                   &
                        ! code saying which departure point scheme to
                        ! use.
     &, Depart_order                                                    &
                        ! for the chosen departure point scheme how
                        ! many iterations/terms to use.
     &, Robert_high_order_scheme                                        &
                                 ! code choosing high order
                          ! interpolation scheme used in Robert routine
     &, Robert_monotone_scheme                                          &
                                 ! code choosing monotone
                          ! interpolation scheme used in Robert routine
     &, interp_vertical_search_tol                                      &
                                   ! used in interpolation code.
     &, check_bottom_levels                                             &
                            ! used in interpolation code, and is
                            ! the number of levels to check to see
                            ! if the departure point lies inside the
                            ! orography.
     &, first_constant_r_rho_level

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Logical                                                           &
     &  L_2d_sl_geometry                                                &
                         ! True, then only perform vector co-ordinate
                         !       geometry in 2d.
     &, L_high                                                          &
                       ! True, if high order interpolation required.
     &, L_mono                                                          &
                       ! True, if interpolation required to be monotone.
     &, L_conserv                                                       &
                       ! True, if interpolation to be monotone and
                       !       conservative.
     &, L_interp_depart                                                 &
     &, L_regular

      Logical :: L_trivial_trigs ! True if trivial_trigs (Cartesian grid)

      Logical                                                           &
     &  L_Robert_high                                                   &
                      ! True if high order interpolation scheme to be
                      ! used in Robert scheme
     &, L_Robert_mono ! True if monotone interpolation scheme to be
                      ! used in Robert scheme

      Real                                                              &
     &  delta_lambda                                                    &
                         ! grid-length in lambda direction
     &, delta_phi                                                       &
                         ! grid-length in phi direction
     &, base_lambda                                                     &
     &, base_phi                                                        &
     &, lambda_p_end                                                    &
     &, phi_p_end                                                       &
     &, dlambda_p_end                                                   &
     &, dphi_p_end                                                      &
     &, dphi_v_end                                                      &
     &, recip_dlam                                                      &
     &, recip_dphi                                                      &
     &, alpha_3                                                         &
     &, alpha_4                                                         &
     &, timestep

! look-up table halos
      Integer                                                           &
     &  halo_lam                                                        &
     &, halo_phi

!VarRes horizontal co-ordinate look-up table
      Integer                                                           &
     &  look_lam(1-halo_lam : max_look-halo_lam)                        &
     &, look_phi(1-halo_phi : max_look-halo_phi)

!VarRes horizontal co-ordinate information etc.
      Real                                                              &
     &  glambda_p( 1-halo_i : global_row_length+halo_i)                 &
     &, glambda_u( 1-halo_i : global_row_length+halo_i)                 &
     &, phi_p    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : rows + halo_j )                           &
     &, phi_v    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : n_rows+halo_j )                           &
     &, gdlambda_p(1-halo_i : global_row_length+halo_i)                 &
     &, gdlambda_u(1-halo_i : global_row_length+halo_i)                 &
     &, dphi_p  ( 1-halo_i : row_length + halo_i                        &
     &,           1-halo_j : rows+halo_j )                              &
     &, dphi_v  ( 1-halo_i : row_length + halo_i                        &
     &,           1-halo_j : n_rows+halo_j )                            &
     &, grecip_dlamp(1-halo_i : global_row_length + halo_i)             &
     &, grecip_dlamu(1-halo_i : global_row_length + halo_i)             &
     &, recip_dphip( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, recip_dphiv( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : n_rows+halo_j )                         &
     &, wt_lambda_p(1-halo_i:row_length+halo_i)                         &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)      &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)    &
     &, lambda_p_rm (1-halo_i : row_length + halo_i)                    &
     &, lambda_p_rp (1-halo_i : row_length + halo_i)                    &
     &, lambda_u_rm (1-halo_i : row_length + halo_i)                    &
     &, lambda_u_rp (1-halo_i : row_length + halo_i)                    &
     &, phi_p_rm   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, phi_p_rp   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, phi_v_rm   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : n_rows+halo_j )                         &
     &, phi_v_rp   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : n_rows+halo_j )                         &
     &, recip_lambda_p_m(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_0(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_p(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_p2(1-halo_i : row_length+halo_i)                 &
     &, recip_lambda_u_m(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_u_0(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_u_p(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_u_p2(1-halo_i : row_length+halo_i)                 &
     &, recip_phi_p_m ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_0 ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_p ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_p2( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_v_m ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : n_rows+halo_j )                      &
     &, recip_phi_v_0 ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : n_rows+halo_j )                      &
     &, recip_phi_v_p ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : n_rows+halo_j )                      &
     &, recip_phi_v_p2( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : n_rows+halo_j )

      Real                                                              &
           ! trigonometric functions
     &  cos_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
     &, sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
     &, sin_theta_latitude (row_length, rows)                           &
     &, tan_theta_latitude (row_length, rows)                           &
     &, cos_v_latitude (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y) &
     &, sec_v_latitude (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y) &
     &, sin_v_latitude (row_length, n_rows)                             &
     &, tan_v_latitude (row_length, n_rows)                             &
     &, cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)

      Real                                                              &
           ! primary model variables
     &  v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y, model_levels) &
     &, u_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         model_levels)                                            &
     &, v_adv (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,      &
     &         model_levels)                                            &
     &, w_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         0:model_levels)

      Real                                                              &
     &  Yu (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      model_levels, n_Y_arrays)                                   &
     &, Yv (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,         &
     &      model_levels, n_Y_arrays)                                   &
     &, Yw (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      0:model_levels)

      Real                                                              &
           ! ! model level arrays
     &  eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent INOUT. ie: Input variables changed on output.

      Real                                                              &
           ! See WP154 and WP162 for definitions.
     &  R_v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,            &
     &       model_levels)

! Arguments with Intent OUT. ie: Output variables.

      Integer                                                           &
     &  Error_Code     ! Non-zero on exit if error detected.

! Local Variables.

      Integer                                                           &
     &  i, j, k, gi                                                     &
                      ! Loop indices
     &, imin, imax, jmin, jmax

      Integer                                                           &
     &  type                                                            &
                       ! a code saying what points the grid points are
                       ! on.
     &, number_of_inputs !the number of fields to interpolate in any one
                         ! call to the interpolation routine.

      Logical                                                           &
     &  L_vector       ! True, if data is a horizontal vector component,
                       ! False, then data is a scalar.

      Real                                                              &
     &  sin_lambda_ad                                                   &
                       ! sine(lambda_a - lambda_d)
     &, cos_lambda_ad                                                   &
                       ! cosine(lambda_a - lambda_d)
     &, cos_phi_d                                                       &
     &, sin_phi_d                                                       &
     &, trig_u                                                          &
     &, trig_v                                                          &
     &, trig_w                                                          &
     &, dummy                                                           &
     &, delta_r        ! dummy array

      Real                                                              &
     &  depart_lambda (row_length, n_rows, model_levels)                &
                                                        ! Lambda
                                                      ! co-ordinate of
                                                      ! departure point.
     &, depart_phi (row_length, n_rows, model_levels)                   &
                                                      ! Phi Co-ordinate
                                                     ! of co-ordinate of
                                                      ! departure point.
     &, depart_r (row_length, n_rows, model_levels)                     &
                                                      ! Vertical
                                                      ! co-ordinate of
                                                      ! departure point.
     &, lambda_a(row_length)


      Real                                                              &
     &  Yu_d (row_length, n_rows, model_levels, n_Yd_arrays)            &
     &, Yv_d (row_length, n_rows, model_levels, n_Yd_arrays)            &
     &, Yw_d (row_length, n_rows, model_levels, n_Ydw_arrays)

      Integer                                                           &
     &  i_out (row_length, n_rows, model_levels)                        &
     &, j_out (row_length, n_rows, model_levels)

      Real                                                              &
     &  weight_lambda (row_length, n_rows, model_levels)                &
     &, weight_phi    (row_length, n_rows, model_levels)

      Real                                                              &
     &  tmp1          ! temporary variable to handle int -> real

#include "parparm.h"
#include "domtyp.h"
! External Routines:
      External                                                          &
     &  Departure_Point, Interpolation

! ----------------------------------------------------------------------
! Section 1.   Calculate Departure_Point
! ----------------------------------------------------------------------

      type = 2  ! v points.
      Error_Code=0

      IF(.NOT. L_INTERP_DEPART)then
! DEPENDS ON: departure_point
      Call Departure_Point(                                             &
     &                     type, timestep, u_adv, v_adv, w_adv,         &
     &                     eta_rho_levels, eta_theta_levels,            &
     &                     r_rho_levels, r_theta_levels,                &
     &                     r_at_u, r_at_v,                              &
     &                     row_length, rows, n_rows, model_levels,      &
     &                     delta_lambda, delta_phi,                     &
     &                     glambda_p, phi_p, glambda_u, phi_v,          &
     &                     gdlambda_p, dphi_p, gdlambda_u, dphi_v,      &
     &                     grecip_dlamp, recip_dphip, grecip_dlamu,     &
     &                     recip_dphiv, wt_lambda_p, wt_phi_p,          &
     &                     wt_lambda_u, wt_phi_v, lambda_p_rm,          &
     &                     lambda_p_rp, lambda_u_rm, lambda_u_rp,       &
     &                     phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,      &
     &                     recip_lambda_p_m, recip_lambda_p_0,          &
     &                     recip_lambda_p_p, recip_lambda_p_p2,         &
     &                     recip_lambda_u_m, recip_lambda_u_0,          &
     &                     recip_lambda_u_p, recip_lambda_u_p2,         &
     &                     recip_phi_p_m, recip_phi_p_0,                &
     &                     recip_phi_p_p, recip_phi_p_p2,               &
     &                     recip_phi_v_m, recip_phi_v_0,                &
     &                     recip_phi_v_p, recip_phi_v_p2,               &
     &                     Base_lambda, base_phi,                       &
     &                     lambda_p_end, phi_p_end, dlambda_p_end,      &
     &                     dphi_p_end, dphi_v_end,                      &
     &                     recip_dlam, recip_dphi, max_look,            &
     &                     look_lam, look_phi, halo_lam, halo_phi,      &
     &                     cos_theta_latitude, sec_theta_latitude,      &
     &                     sin_theta_latitude, tan_theta_latitude,      &
     &                     cos_v_latitude, sec_v_latitude,              &
     &                     sin_v_latitude, tan_v_latitude,              &
     &                     sin_theta_longitude,                         &
     &                     cos_theta_longitude,                         &
     &                     model_domain, Pi, Depart_scheme,             &
     &                     Depart_order, L_Robert_high, L_Robert_mono,  &
     &                     Robert_high_order_scheme,                    &
     &                     Robert_monotone_scheme,                      &
     &                     check_bottom_levels,                         &
     &                     interp_vertical_search_tol,                  &
     &                     L_2d_sl_geometry, L_regular,                 &
     &                     first_constant_r_rho_level,                  &
     &                     LAM_max_cfl, n_rows,                         &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     off_x, off_y, halo_i, halo_j, l_datastart,   &
     &                     global_row_length, global_rows,              &
     &                     g_i_pe, proc_row_group,                      &
     &                     proc_col_group, at_extremity,                &
     &                     L_sl_halo_reprod,                            &
     &                     depart_lambda, depart_phi, depart_r)
      endif

! ----------------------------------------------------------------------
! Section 2.  Call interpolation routine for each term.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 2.1 Call interpolation routine for Yu and u.
! ----------------------------------------------------------------------
      If (Error_Code  ==  0 ) Then

! DEPENDS ON: calc_index
        Call Calc_Index(                                                &
     &                      row_length, n_rows, model_levels,           &
     &                      delta_lambda, delta_phi,                    &
     &                      base_lambda, base_phi,                      &
     &                      glambda_p, phi_p, grecip_dlamp, recip_dphip,&
     &                      recip_dlam, recip_dphi, max_look,           &
     &                      look_lam, look_phi, halo_lam, halo_phi,     &
     &                      L_regular, depart_lambda, depart_phi,       &
     &                      halo_i, halo_j,                             &
     &                      global_row_length,                          &
     &                      row_length, rows, l_datastart,              &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi)

        L_vector = .true.
        type = 1 ! data is at u points.
        number_of_inputs = 1

        If (model_domain  ==  mt_Global .or.                            &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) Then
          number_of_inputs = 1
! DEPENDS ON: interpolation
          Call Interpolation(                                           &
     &                     Yu, dummy, dummy,                            &
     &                     eta_rho_levels,                              &
     &                     r_at_u, delta_r, r_rho_levels, type,         &
     &                     number_of_inputs,                            &
     &                     check_bottom_levels,                         &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level,                  &
     &                     row_length, rows, model_levels,              &
     &                     rows,                                        &
     &                     row_length, n_rows, model_levels,            &
     &                     delta_lambda, delta_phi,                     &
     &                     base_lambda, base_phi,                       &
     &                     glambda_u, phi_v,                            &
     &                     gdlambda_p, dphi_p, gdlambda_u, dphi_v,      &
     &                     grecip_dlamu, recip_dphiv,                   &
     &                     lambda_u_rm, lambda_u_rp,                    &
     &                     phi_p_rm, phi_p_rp,                          &
     &                     recip_lambda_u_m, recip_lambda_u_0,          &
     &                     recip_lambda_u_p, recip_lambda_u_p2,         &
     &                     recip_phi_p_m, recip_phi_p_0,                &
     &                     recip_phi_p_p, recip_phi_p_p2,               &
     &                     i_out, j_out,                                &
     &                     weight_lambda, weight_phi,                   &
     &                     high_order_scheme, monotone_scheme,          &
     &                     cos_theta_latitude, L_regular,               &
     &                     L_vector, model_domain, L_high, L_mono,      &
     &                     .false.,                                     &
     &                     depart_r, depart_lambda, depart_phi,         &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j,                              &
     &                     global_row_length, global_rows,              &
     &                     row_length, rows, n_rows, rows,              &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     proc_row_group, proc_col_group, 2, 0, 0,     &
     &                     L_sl_halo_reprod, off_x,off_y,               &
     &                     Yu_d, dummy, dummy, Error_Code)

        Else if (model_domain  ==  mt_LAM .and.                         &
     &           L_2d_sl_geometry) Then

          number_of_inputs = 3
! DEPENDS ON: interpolation
          Call Interpolation(                                           &
     &                     Yu(1-halo_i, 1-halo_j, 1 ,3),                &
     &                     Yu(1-halo_i, 1-halo_j, 1 ,1),                &
     &                     Yu(1-halo_i, 1-halo_j, 1 ,2),                &
     &                     eta_rho_levels,                              &
     &                     r_at_u, delta_r, r_rho_levels, type,         &
     &                     number_of_inputs,                            &
     &                     check_bottom_levels,                         &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level,                  &
     &                     row_length, rows, model_levels,              &
     &                     rows,                                        &
     &                     row_length, n_rows, model_levels,            &
     &                     delta_lambda, delta_phi,                     &
     &                     base_lambda, base_phi,                       &
     &                     glambda_u, phi_v,                            &
     &                     gdlambda_p, dphi_p, gdlambda_u, dphi_v,      &
     &                     grecip_dlamu, recip_dphiv,                   &
     &                     lambda_u_rm, lambda_u_rp,                    &
     &                     phi_p_rm, phi_p_rp,                          &
     &                     recip_lambda_u_m, recip_lambda_u_0,          &
     &                     recip_lambda_u_p, recip_lambda_u_p2,         &
     &                     recip_phi_p_m, recip_phi_p_0,                &
     &                     recip_phi_p_p, recip_phi_p_p2,               &
     &                     i_out, j_out,                                &
     &                     weight_lambda, weight_phi,                   &
     &                     high_order_scheme, monotone_scheme,          &
     &                     cos_theta_latitude, L_regular,               &
     &                     L_vector, model_domain, L_high, L_mono,      &
     &                     .false.,                                     &
     &                     depart_r, depart_lambda, depart_phi,         &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j,                              &
     &                     global_row_length, global_rows,              &
     &                     row_length, rows, n_rows, rows,              &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     proc_row_group, proc_col_group, 2, 0, 0,     &
     &                     L_sl_halo_reprod,off_x,off_y,                &
     &                     Yu_d(1,1,1,1), Yu_d(1,1,1,2),                &
     &                     Yu_d(1,1,1,3), Error_Code)

        End If
      End If

! ----------------------------------------------------------------------
! Section 2.2 Call interpolation routine for Yv and v.
! ----------------------------------------------------------------------

      If (Error_Code  ==  0 ) Then

        type = 2 ! data at v points.
        L_vector = .true.

        If (model_domain  ==  mt_Global .or.                            &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) Then
          number_of_inputs = 1
! DEPENDS ON: interpolation
          Call Interpolation(                                           &
     &                     Yv, dummy, dummy,                            &
     &                     eta_rho_levels,                              &
     &                     r_at_v, delta_r, r_rho_levels, type,         &
     &                     number_of_inputs,                            &
     &                     check_bottom_levels,                         &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level,                  &
     &                     row_length, n_rows, model_levels,            &
     &                     rows,                                        &
     &                     row_length, n_rows, model_levels,            &
     &                     delta_lambda, delta_phi,                     &
     &                     base_lambda, base_phi,                       &
     &                     glambda_u, phi_v,                            &
     &                     gdlambda_p, dphi_p, gdlambda_u, dphi_v,      &
     &                     grecip_dlamu, recip_dphiv,                   &
     &                     lambda_p_rm, lambda_p_rp,                    &
     &                     phi_v_rm, phi_v_rp,                          &
     &                     recip_lambda_p_m, recip_lambda_p_0,          &
     &                     recip_lambda_p_p, recip_lambda_p_p2,         &
     &                     recip_phi_v_m, recip_phi_v_0,                &
     &                     recip_phi_v_p, recip_phi_v_p2,               &
     &                     i_out, j_out,                                &
     &                     weight_lambda, weight_phi,                   &
     &                     high_order_scheme, monotone_scheme,          &
     &                     cos_v_latitude, L_regular,                   &
     &                     L_vector, model_domain, L_high, L_mono,      &
     &                     .false.,                                     &
     &                     depart_r, depart_lambda, depart_phi,         &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j,                              &
     &                     global_row_length, global_rows,              &
     &                     row_length, rows, n_rows, n_rows,            &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     proc_row_group, proc_col_group, 2, 0, 0,     &
     &                     L_sl_halo_reprod,off_x,off_y,                &
     &                     Yv_d, dummy, dummy, Error_Code)

        Else

          number_of_inputs = 3

! DEPENDS ON: interpolation
          Call Interpolation(                                           &
     &                     Yv(1-halo_i, 1-halo_j, 1 ,3),                &
     &                     Yv(1-halo_i, 1-halo_j, 1 ,1),                &
     &                     Yv(1-halo_i, 1-halo_j, 1 ,2),                &
     &                     eta_rho_levels,                              &
     &                     r_at_v, delta_r, r_rho_levels, type,         &
     &                     number_of_inputs,                            &
     &                     check_bottom_levels,                         &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level,                  &
     &                     row_length, n_rows, model_levels,            &
     &                     rows,                                        &
     &                     row_length, n_rows, model_levels,            &
     &                     delta_lambda, delta_phi,                     &
     &                     base_lambda, base_phi,                       &
     &                     glambda_u, phi_v,                            &
     &                     gdlambda_p, dphi_p, gdlambda_u, dphi_v,      &
     &                     grecip_dlamu, recip_dphiv,                   &
     &                     lambda_p_rm, lambda_p_rp,                    &
     &                     phi_v_rm, phi_v_rp,                          &
     &                     recip_lambda_p_m, recip_lambda_p_0,          &
     &                     recip_lambda_p_p, recip_lambda_p_p2,         &
     &                     recip_phi_v_m, recip_phi_v_0,                &
     &                     recip_phi_v_p, recip_phi_v_p2,               &
     &                     i_out, j_out,                                &
     &                     weight_lambda, weight_phi,                   &
     &                     high_order_scheme, monotone_scheme,          &
     &                     cos_v_latitude, L_regular,                   &
     &                     L_vector, model_domain, L_high, L_mono,      &
     &                     .false.,                                     &
     &                     depart_r, depart_lambda, depart_phi,         &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j,                              &
     &                     global_row_length, global_rows,              &
     &                     row_length, rows, n_rows, n_rows,            &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     proc_row_group, proc_col_group, 2, 0, 0,     &
     &                     L_sl_halo_reprod,off_x,off_y,                &
     &                     Yv_d(1,1,1,1), Yv_d(1,1,1,2),                &
     &                     Yv_d(1,1,1,3), Error_Code)

        End If
      End If

! ----------------------------------------------------------------------
! Section 2.3 Call interpolation routine for Yw and w.
! ----------------------------------------------------------------------

      If (.not. L_2d_sl_geometry .and.                                  &
     &    (model_domain  ==  mt_Global .or.                             &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) ) Then
        If (Error_Code  ==  0 ) Then

          type = 3 ! data at w points.
          L_vector = .false.

          number_of_inputs = 1
! DEPENDS ON: interpolation
          Call Interpolation(                                           &
     &                     Yw, dummy, dummy,                            &
     &                     eta_theta_levels,                            &
     &                     r_theta_levels, delta_r, r_theta_levels,     &
     &                     type, number_of_inputs,                      &
     &                     check_bottom_levels,                         &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level+1,                &
     &                     row_length, rows, model_levels+1,            &
     &                     rows,                                        &
     &                     row_length, n_rows, model_levels,            &
     &                     delta_lambda, delta_phi,                     &
     &                     base_lambda, base_phi,                       &
     &                     glambda_u, phi_v,                            &
     &                     gdlambda_p, dphi_p, gdlambda_u, dphi_v,      &
     &                     grecip_dlamu, recip_dphiv,                   &
     &                     lambda_p_rm, lambda_p_rp,                    &
     &                     phi_p_rm, phi_p_rp,                          &
     &                     recip_lambda_p_m, recip_lambda_p_0,          &
     &                     recip_lambda_p_p, recip_lambda_p_p2,         &
     &                     recip_phi_p_m, recip_phi_p_0,                &
     &                     recip_phi_p_p, recip_phi_p_p2,               &
     &                     i_out, j_out,                                &
     &                     weight_lambda, weight_phi,                   &
     &                     high_order_scheme, monotone_scheme,          &
     &                     cos_theta_latitude, L_regular,               &
     &                     L_vector, model_domain, L_high, L_mono,      &
     &                     .false.,                                     &
     &                     depart_r, depart_lambda, depart_phi,         &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j,                              &
     &                     global_row_length, global_rows,              &
     &                     row_length, rows, n_rows, rows,              &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     proc_row_group, proc_col_group, 2, 0, 0,     &
     &                     L_sl_halo_reprod,off_x,off_y,                &
     &                     Yw_d, dummy, dummy, Error_Code)

        End If

      End If ! on 2d geometry

! ----------------------------------------------------------------------
! Section 3.  Form R_v from the interpolated terms using
!             expressions given in WP 162.
! ----------------------------------------------------------------------

      If ( L_regular ) then
        Do i = 1, row_length
          gi = l_datastart(1) + i - 1
          tmp1 = float(gi-1)
          lambda_a(i) = tmp1 * delta_lambda
        endDo
      else ! variable resolution
        Do i = 1, row_length
          gi = l_datastart(1) + i - 1
          lambda_a(i) = glambda_p(gi) - Base_lambda
        endDo
      end If ! L_regular
      If (Error_Code  ==  0 ) Then
        If (model_domain  ==  mt_Global .or.                            &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) Then
          If (L_2d_sl_geometry )Then
            Do k = 1, model_levels
              Do j = 1, n_rows
                Do i = 1, row_length
                  sin_lambda_ad= sin(lambda_a(i) - depart_lambda(i,j,k))
                  cos_lambda_ad= cos(lambda_a(i) - depart_lambda(i,j,k))
                  cos_phi_d = cos ( depart_phi(i,j,k))
                  sin_phi_d = sin ( depart_phi(i,j,k))

                  trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                  trig_v = cos_v_latitude(i,j) * cos_phi_d +            &
     &                     sin_v_latitude(i,j) * sin_phi_d *            &
     &                     cos_lambda_ad

                  R_v(i,j,k) = trig_u * Yu_d(i,j,k,1)                   &
     &                       + trig_v * Yv_d(i,j,k,1)                   &
     &                       + R_v(i,j,k)

                End Do
              End Do
            End Do

          Elseif(L_Trivial_trigs)Then
            Do k = 1, model_levels
              Do j = 1, n_rows
                Do i = 1, row_length
                  R_v(i,j,k) = Yv_d(i,j,k,1) + R_v(i,j,k)
                End Do
              End Do
            End Do

          Else
            Do k = 1, model_levels
              Do j = 1, n_rows
                Do i = 1, row_length
                  sin_lambda_ad= sin(lambda_a(i) - depart_lambda(i,j,k))
                  cos_lambda_ad= cos(lambda_a(i) - depart_lambda(i,j,k))
                  cos_phi_d = cos ( depart_phi(i,j,k))
                  sin_phi_d = sin ( depart_phi(i,j,k))

                  trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                  trig_v = cos_v_latitude(i,j) * cos_phi_d +            &
     &                     sin_v_latitude(i,j) * sin_phi_d *            &
     &                     cos_lambda_ad
                  trig_w = cos_v_latitude(i,j) * sin_phi_d -            &
     &                     sin_v_latitude(i,j) * cos_phi_d *            &
     &                     cos_lambda_ad

                  R_v(i,j,k) = trig_u * Yu_d(i,j,k,1)                   &
     &                       + trig_v * Yv_d(i,j,k,1)                   &
     &                       + trig_w * Yw_d(i,j,k,1)                   &
     &                       + R_v(i,j,k)

                End Do
              End Do
            End Do
          End If

        Else If (model_domain  ==  mt_LAM ) Then

          If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then

          If (L_2d_sl_geometry )Then

! internal area uses real alphas
            Do k = 1, model_levels
              Do j = 1, n_rows
                Do i = 1, row_length
                  sin_lambda_ad= sin(lambda_a(i) - depart_lambda(i,j,k))
                  cos_lambda_ad= cos(lambda_a(i) - depart_lambda(i,j,k))
                  cos_phi_d = cos ( depart_phi(i,j,k))
                  sin_phi_d = sin ( depart_phi(i,j,k))

                  trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                  trig_v = cos_v_latitude(i,j) * cos_phi_d +            &
     &                   sin_v_latitude(i,j) * sin_phi_d *              &
     &                   cos_lambda_ad

                  R_v(i,j,k) = trig_u * ( Yu_d(i,j,k,1) +               &
     &                       Yu_d(i,j,k,2) * (1. - alpha_3) -           &
     &                       Yu_d(i,j,k,3) * (1. - alpha_4) )           &
     &                     + trig_v *                                   &
     &                      (Yv_d(i,j,k,1) -                            &
     &                       Yv_d(i,j,k,2) * (1. - alpha_3) +           &
     &                       Yv_d(i,j,k,3) * (1. - alpha_4) )           &
     &                     - alpha_3 * Yv(i,j,k,1)                      &
     &                     + alpha_4 * Yv(i,j,k,2)                      &
     &                     - v(i,j,k)

                End Do
              End Do
            End Do

! boundary area uses alphas = 1
! southern area
            If (at_extremity(PSouth)) Then
              Do k = 1, model_levels
                Do j = 1, halo_j
                  Do i = 1, row_length
                  sin_lambda_ad= sin(lambda_a(i) - depart_lambda(i,j,k))
                  cos_lambda_ad= cos(lambda_a(i) - depart_lambda(i,j,k))
                    cos_phi_d = cos ( depart_phi(i,j,k))
                    sin_phi_d = sin ( depart_phi(i,j,k))

                    trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                    trig_v = cos_v_latitude(i,j) * cos_phi_d +          &
     &                       sin_v_latitude(i,j) * sin_phi_d *          &
     &                       cos_lambda_ad

                    R_v(i,j,k) = trig_u * Yu_d(i,j,k,1)                 &
     &                         + trig_v * Yv_d(i,j,k,1)                 &
     &                         - Yv(i,j,k,1)                            &
     &                         + Yv(i,j,k,2)                            &
     &                         - v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Western Area
            If (at_extremity(PWest)) Then
              Do k = 1, model_levels
                Do j = 1, n_rows
                  Do i = 1, halo_i
                  sin_lambda_ad= sin(lambda_a(i) - depart_lambda(i,j,k))
                  cos_lambda_ad= cos(lambda_a(i) - depart_lambda(i,j,k))
                    cos_phi_d = cos ( depart_phi(i,j,k))
                    sin_phi_d = sin ( depart_phi(i,j,k))

                    trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                    trig_v = cos_v_latitude(i,j) * cos_phi_d +          &
     &                       sin_v_latitude(i,j) * sin_phi_d *          &
     &                       cos_lambda_ad

                    R_v(i,j,k) = trig_u * Yu_d(i,j,k,1)                 &
     &                         + trig_v * Yv_d(i,j,k,1)                 &
     &                         - Yv(i,j,k,1)                            &
     &                         + Yv(i,j,k,2)                            &
     &                         - v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Eastern Area
            If (at_extremity(PEast)) Then
              Do k = 1, model_levels
                Do j = 1, n_rows
                  Do i = row_length-halo_i+1, row_length
                  sin_lambda_ad= sin(lambda_a(i) - depart_lambda(i,j,k))
                  cos_lambda_ad= cos(lambda_a(i) - depart_lambda(i,j,k))
                    cos_phi_d = cos ( depart_phi(i,j,k))
                    sin_phi_d = sin ( depart_phi(i,j,k))

                    trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                    trig_v = cos_v_latitude(i,j) * cos_phi_d +          &
     &                       sin_v_latitude(i,j) * sin_phi_d *          &
     &                       cos_lambda_ad

                    R_v(i,j,k) = trig_u * Yu_d(i,j,k,1)                 &
     &                         + trig_v * Yv_d(i,j,k,1)                 &
     &                         - Yv(i,j,k,1)                            &
     &                         + Yv(i,j,k,2)                            &
     &                         - v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Northern Area
            If (at_extremity(PNorth)) Then
              Do k = 1, model_levels
                Do j = n_rows-halo_j, n_rows
                  Do i = 1, row_length
                  sin_lambda_ad= sin(lambda_a(i) - depart_lambda(i,j,k))
                  cos_lambda_ad= cos(lambda_a(i) - depart_lambda(i,j,k))
                    cos_phi_d = cos ( depart_phi(i,j,k))
                    sin_phi_d = sin ( depart_phi(i,j,k))

                    trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                    trig_v = cos_v_latitude(i,j) * cos_phi_d +          &
     &                       sin_v_latitude(i,j) * sin_phi_d *          &
     &                       cos_lambda_ad

                    R_v(i,j,k) = trig_u * Yu_d(i,j,k,1)                 &
     &                         + trig_v * Yv_d(i,j,k,1)                 &
     &                         - Yv(i,j,k,1)                            &
     &                         + Yv(i,j,k,2)                            &
     &                         - v(i,j,k)
                  End Do
                End Do
              End Do
            End If

          Else
! No geometry
! internal area uses real alphas
            Do k = 1, model_levels
              Do j = 1, n_rows
                Do i = 1, row_length
                  R_v(i,j,k) =                                          &
     &                       Yv_d(i,j,k,1) -                            &
     &                       Yv_d(i,j,k,2) * (1. - alpha_3) +           &
     &                       Yv_d(i,j,k,3) * (1. - alpha_4)             &
     &                     - alpha_3 * Yv(i,j,k,1)                      &
     &                     + alpha_4 * Yv(i,j,k,2)                      &
     &                     - v(i,j,k)
                End Do
              End Do
            End Do

! boundary area uses alphas = 1
! southern area
            If (at_extremity(PSouth)) Then
              Do k = 1, model_levels
                Do j = 1, halo_j
                  Do i = 1, row_length

                    R_v(i,j,k) = Yv_d(i,j,k,1)                          &
     &                         - Yv(i,j,k,1)                            &
     &                         + Yv(i,j,k,2)                            &
     &                         - v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Western Area
            If (at_extremity(PWest)) Then
              Do k = 1, model_levels
                Do j = 1, n_rows
                  Do i = 1, halo_i

                    R_v(i,j,k) = Yv_d(i,j,k,1)                          &
     &                         - Yv(i,j,k,1)                            &
     &                         + Yv(i,j,k,2)                            &
     &                         - v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Eastern Area
            If (at_extremity(PEast)) Then
              Do k = 1, model_levels
                Do j = 1, n_rows
                  Do i = row_length-halo_i+1, row_length
                    R_v(i,j,k) = Yv_d(i,j,k,1)                          &
     &                         - Yv(i,j,k,1)                            &
     &                         + Yv(i,j,k,2)                            &
     &                         - v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Northern Area
            If (at_extremity(PNorth)) Then
              Do k = 1, model_levels
                Do j = n_rows-halo_j, n_rows
                  Do i = 1, row_length
                    R_v(i,j,k) = Yv_d(i,j,k,1)                          &
     &                         - Yv(i,j,k,1)                            &
     &                         + Yv(i,j,k,2)                            &
     &                         - v(i,j,k)
                  End Do
                End Do
              End Do
            End If

          End If ! on geometry

          Else ! If ( CycleNo > 1 and L_new_tdisc ) then

          If (L_2d_sl_geometry )Then
!
! First do the boundaries to set appropriate loop indices
! thus avoiding overwriting R_u or using temp storage
!
! Loop indices for R_u update when processor is away from boundaries
!
            imin = 1
            imax = row_length
            jmin = 1
            jmax = n_rows
!
! boundary area uses alphas = 1
!
! southern area
            If (at_extremity(PSouth)) Then
              jmin = halo_j+1
              Do k = 1, model_levels
                Do j = 1, halo_j
                  Do i = 1, row_length
                    sin_lambda_ad= sin(lambda_a(i)-depart_lambda(i,j,k))
                    cos_lambda_ad= cos(lambda_a(i)-depart_lambda(i,j,k))
                    cos_phi_d = cos ( depart_phi(i,j,k))
                    sin_phi_d = sin ( depart_phi(i,j,k))
                    trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                    trig_v = cos_v_latitude(i,j) * cos_phi_d +          &
     &                       sin_v_latitude(i,j) * sin_phi_d *          &
     &                       cos_lambda_ad
                    R_v(i,j,k) = trig_u * Yu_d(i,j,k,1)                 &
     &                         + trig_v * Yv_d(i,j,k,1) + R_v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Western Area
            If (at_extremity(PWest)) Then
              imin = halo_i+1
              Do k = 1, model_levels
! use jmin to avoid overwriting S boundary
                Do j = jmin, n_rows
                  Do i = 1, halo_i
                    sin_lambda_ad= sin(lambda_a(i)-depart_lambda(i,j,k))
                    cos_lambda_ad= cos(lambda_a(i)-depart_lambda(i,j,k))
                    cos_phi_d = cos ( depart_phi(i,j,k))
                    sin_phi_d = sin ( depart_phi(i,j,k))
                    trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                    trig_v = cos_v_latitude(i,j) * cos_phi_d +          &
     &                       sin_v_latitude(i,j) * sin_phi_d *          &
     &                       cos_lambda_ad
                    R_v(i,j,k) = trig_u * Yu_d(i,j,k,1)                 &
     &                         + trig_v * Yv_d(i,j,k,1) + R_v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Eastern Area
            If (at_extremity(PEast)) Then
              imax = row_length-halo_i
              Do k = 1, model_levels
! use jmin to avoid overwriting S boundary
                Do j = jmin, n_rows
                  Do i = row_length-halo_i+1, row_length
                    sin_lambda_ad= sin(lambda_a(i)-depart_lambda(i,j,k))
                    cos_lambda_ad= cos(lambda_a(i)-depart_lambda(i,j,k))
                    cos_phi_d = cos ( depart_phi(i,j,k))
                    sin_phi_d = sin ( depart_phi(i,j,k))
                    trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                    trig_v = cos_v_latitude(i,j) * cos_phi_d +          &
     &                       sin_v_latitude(i,j) * sin_phi_d *          &
     &                       cos_lambda_ad
                    R_v(i,j,k) = trig_u * Yu_d(i,j,k,1)                 &
     &                         + trig_v * Yv_d(i,j,k,1) + R_v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Northern Area
            If (at_extremity(PNorth)) Then
              jmax = n_rows-halo_j-1
              Do k = 1, model_levels
                Do j = n_rows-halo_j, n_rows
! use imin, imax to avoid overwriting E and W boundary
                  Do i = imin, imax
                    sin_lambda_ad= sin(lambda_a(i)-depart_lambda(i,j,k))
                    cos_lambda_ad= cos(lambda_a(i)-depart_lambda(i,j,k))
                    cos_phi_d = cos ( depart_phi(i,j,k))
                    sin_phi_d = sin ( depart_phi(i,j,k))
                    trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                    trig_v = cos_v_latitude(i,j) * cos_phi_d +          &
     &                       sin_v_latitude(i,j) * sin_phi_d *          &
     &                       cos_lambda_ad
                    R_v(i,j,k) = trig_u * Yu_d(i,j,k,1)                 &
     &                         + trig_v * Yv_d(i,j,k,1) + R_v(i,j,k)
                  End Do
                End Do
              End Do
            End If

! internal area uses real alphas
            Do k = 1, model_levels
              Do j = jmin, jmax
                Do i = imin, imax
                  sin_lambda_ad= sin(lambda_a(i)-depart_lambda(i,j,k))
                  cos_lambda_ad= cos(lambda_a(i)-depart_lambda(i,j,k))
                  cos_phi_d = cos ( depart_phi(i,j,k))
                  sin_phi_d = sin ( depart_phi(i,j,k))

                  trig_u = - sin_lambda_ad * sin_v_latitude(i,j)
                  trig_v = cos_v_latitude(i,j) * cos_phi_d +            &
     &                   sin_v_latitude(i,j) * sin_phi_d *              &
     &                   cos_lambda_ad

                  R_v(i,j,k) = trig_u * ( Yu_d(i,j,k,1) +               &
     &                       Yu_d(i,j,k,2) * (1. - alpha_3) -           &
     &                       Yu_d(i,j,k,3) * (1. - alpha_4) )           &
     &                     + trig_v *                                   &
     &                      (Yv_d(i,j,k,1) -                            &
     &                       Yv_d(i,j,k,2) * (1. - alpha_3) +           &
     &                       Yv_d(i,j,k,3) * (1. - alpha_4) )           &
     &                     + R_v(i,j,k)

                End Do
              End Do
            End Do

          Else
! No geometry
!
! First do the boundaries to set appropriate loop indices
! thus avoiding overwriting R_u or using temp storage
!
! Loop indices for R_u update when processor is away from boundaries
!
            imin = 1
            imax = row_length
            jmin = 1
            jmax = n_rows
!
! boundary area uses alphas = 1
!
! southern area
            If (at_extremity(PSouth)) Then
              jmin = halo_j+1
              Do k = 1, model_levels
                Do j = 1, halo_j
                  Do i = 1, row_length
                    R_v(i,j,k) = Yv_d(i,j,k,1) + R_v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Western Area
            If (at_extremity(PWest)) Then
              imin = halo_i+1
              Do k = 1, model_levels
! use jmin to avoid overwriting S boundary
                Do j = jmin, n_rows
                  Do i = 1, halo_i
                    R_v(i,j,k) = Yv_d(i,j,k,1) + R_v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Eastern Area
            If (at_extremity(PEast)) Then
              imax = row_length-halo_i
              Do k = 1, model_levels
! use jmin to avoid overwriting S boundary
                Do j = jmin, n_rows
                  Do i = row_length-halo_i+1, row_length
                    R_v(i,j,k) = Yv_d(i,j,k,1) + R_v(i,j,k)
                  End Do
                End Do
              End Do
            End If
! Northern Area
            If (at_extremity(PNorth)) Then
              jmax = n_rows-halo_j-1
              Do k = 1, model_levels
                Do j = n_rows-halo_j, n_rows
! use imin, imax to avoid overwriting E and W boundary
                  Do i = imin, imax
                    R_v(i,j,k) = Yv_d(i,j,k,1) + R_v(i,j,k)
                  End Do
                End Do
              End Do
            End If
!
! internal area uses real alphas
            Do k = 1, model_levels
              Do j = jmin, jmax
                Do i = imin, imax
                  R_v(i,j,k) = Yv_d(i,j,k,1)                            &
     &                     -  Yv_d(i,j,k,2) * (1. - alpha_3)            &
     &                     +  Yv_d(i,j,k,3) * (1. - alpha_4)            &
     &                     +  R_v(i,j,k)
                End Do
              End Do
            End Do

          End If ! on geometry

          End If ! CycleNo == 1
        End If ! on model domain

      End If ! on error code eq 0

! End of routine.
      return
      END SUBROUTINE SL_Vector_v

#endif
