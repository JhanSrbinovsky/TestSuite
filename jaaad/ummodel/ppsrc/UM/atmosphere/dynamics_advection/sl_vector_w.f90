
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SL_Vector_w
!

      Subroutine SL_Vector_w(                                           &
     &                         w, Yu, Yv, Yw,                           &
     &                         eta_theta_levels, eta_rho_levels,        &
     &                         r_theta_levels, r_rho_levels,            &
     &                         r_at_u, r_at_v,                          &
     &                         depart_lambda_w, depart_phi_w,           &
     &                         depart_r_w,                              &
     &                         cos_theta_latitude, sec_theta_latitude,  &
     &                         sin_theta_latitude,                      &
     &                         cos_theta_longitude,                     &
     &                         sin_theta_longitude,                     &
     &                         delta_lambda, delta_phi,                 &
     &                         glambda_p, phi_p, glambda_u, phi_v,      &
     &                         gdlambda_p, dphi_p, gdlambda_u, dphi_v,  &
     &                         grecip_dlamp, recip_dphip,               &
     &                         grecip_dlamu, recip_dphiv,               &
     &                         lambda_p_rm, lambda_p_rp,                &
     &                         lambda_u_rm, lambda_u_rp,                &
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
     &                         recip_dlam, recip_dphi, max_look,        &
     &                         look_lam, look_phi, halo_lam, halo_phi,  &
     &                         n_Y_arrays, n_Yw_arrays,                 &
     &                         n_Yd_arrays, n_Ydw_arrays,               &
     &                         timestep, alpha_3, alpha_4,              &
     &                         row_length, rows, n_rows, model_levels,  &
     &                         high_order_scheme, monotone_scheme,      &
     &                         model_domain, L_2d_sl_geometry,          &
     &                         L_high, L_mono, L_conserv,               &
     &                         check_bottom_levels,                     &
     &                         interp_vertical_search_tol,              &
     &                         first_constant_r_rho_level,              &
     &                         L_trivial_trigs,                         &
     &                         me, n_proc, n_procx, n_procy,            &
     &                         off_x, off_y, halo_i, halo_j,            &
     &                         global_row_length, global_rows,          &
     &                         l_datastart,  at_extremity,              &
     &                         g_i_pe, proc_row_group, proc_col_group,  &
     &                         L_regular, L_sl_halo_reprod,             &
     &                         L_new_tdisc, CycleNo, R_w, Error_code)

! Purpose:
!          Performs vector semi-Lagrangian integration of values at
!          w points given the forcing functions for the first advection
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
! 9/5/00   5.1         Fix code for North boundary of Mesoscale
!                                                      P.Burton
!LL   5.1   14/02/00  Use DOMTYP parameters                    P.Burton
!LL  5.2 26/09/00  interpolation in eta not z    Andy Malcolm
!LL   5.2   27/09/00  cyclic LAM fix                         A.Malcolm
!   5.3     15/09/01  add mt_bi_cyclic_LAM code            A. Malcolm
! 6.2      04/10/05  Changes for cycling semi-Lagrangian scheme.
!                                                        M. Diamantakis
!  6.2    25/12/05  Variable resolution changes          Yongming Tang
!    6.4 07/02/07 Resolve multiple declaration of cos_lambda_ad
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
                              ! global number of rows
     &, g_i_pe(1-halo_i:global_row_length+halo_i)                       &
                                                  ! processor on my
                    ! processor-row holding a given value in i direction
     &, n_Y_arrays                                                      &
                         ! = 1 for global, 3 for LAM
     &, n_Yw_arrays                                                     &
                         ! = 1 for global, 2 for LAM
     &, n_Yd_arrays                                                     &
                         ! = 1 for global, 3 for LAM
     &, n_Ydw_arrays                                                    &
                         ! = 1 for global, 2 for LAM
     &, CycleNo

      Logical                                                           &
     &  L_sl_halo_reprod                                                &
                         ! if true then sl code bit repoducible with
                         ! any sensible halo size
     &, L_new_tdisc

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
     &, interp_vertical_search_tol                                      &
                                   ! used in interpolation code.
     &, check_bottom_levels ! used in interpolation code, and is
                            ! the number of levels to check to see
                            ! if the departure point lies inside the
                            ! orography.

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, first_constant_r_rho_level

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
     &, L_regular

      Logical :: L_trivial_trigs ! True if trivial_trigs (Cartesian grid)

      Real                                                              &
     &  delta_lambda                                                    &
                         ! grid-length in lambda direction
     &, delta_phi                                                       &
                         ! grid-length in phi direction
     &, base_lambda                                                     &
     &, base_phi                                                        &
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
     &, cos_theta_longitude (row_length, rows)                          &
     &, sin_theta_longitude (row_length, rows)                          &
     &, cos_u_longitude (row_length, rows)                              &
     &, sin_u_longitude (row_length, rows)

      Real                                                              &
           ! primary model variables
     &  w(1-off_x:row_length+off_x, 1-off_y:rows+off_y, 0:model_levels)

      Real                                                              &
           ! co-ordinates of departure points for w.
     &  depart_lambda_w(row_length, rows, model_levels)                 &
     &, depart_phi_w(row_length, rows, model_levels)                    &
     &, depart_r_w(row_length, rows, model_levels)

      Real                                                              &
     &  Yu (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      model_levels, n_Y_arrays)                                   &
     &, Yv (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,         &
     &      model_levels, n_Y_arrays)                                   &
     &, Yw (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      0:model_levels, n_Yw_arrays)

      Real                                                              &
           ! model level arrays
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

      Integer                                                           &
     &  i_out (row_length, rows, model_levels)                          &
     &, j_out (row_length, rows, model_levels)

      Real                                                              &
     &  weight_lambda (row_length, rows, model_levels)                  &
     &, weight_phi    (row_length, rows, model_levels)

      Integer                                                           &
     &  neighbour(4)         ! Array with the Ids of the four neighbours
                             ! in the horizontal plane

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Parameters

! Arguments with Intent INOUT. ie: Input variables changed on output.

      Real                                                              &
           ! See WP154 and WP162 for definitions.
     &  R_w (row_length, rows, model_levels)

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
     & cos_lambda_ad                                                    &
                       ! cosine(lambda_a - lambda_d)
     &, sin_lambda_ad                                                   &
                       ! sine(lambda_a - lambda_d)
     &, cos_phi_d                                                       &
     &, sin_phi_d                                                       &
     &, trig_u                                                          &
     &, trig_v                                                          &
     &, trig_w                                                          &
     &, dummy                                                           &
     &, delta_r        ! dummy array

      Real                                                              &
     &  Yu_d (row_length, rows, model_levels-1, n_Yd_arrays)            &
     &, Yv_d (row_length, rows, model_levels-1, n_Yd_arrays)            &
     &, Yw_d (row_length, rows, model_levels-1, n_Ydw_arrays)           &
     &, lambda_a(row_length)

      Real                                                              &
     &  tmp1          ! temporary variable to handle int -> real

!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end
! External Routines:
      External                                                          &
     &  Interpolation

! ----------------------------------------------------------------------
! Section 1.  Call interpolation routine for each term.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 1.1 Call interpolation routine for Yu.
! ----------------------------------------------------------------------

! DEPENDS ON: calc_index
      Call Calc_Index(                                                  &
     &                      row_length, rows, model_levels-1,           &
     &                      delta_lambda, delta_phi,                    &
     &                      base_lambda, base_phi,                      &
     &                      glambda_p, phi_p, grecip_dlamp, recip_dphip,&
     &                      recip_dlam, recip_dphi, max_look,           &
     &                      look_lam, look_phi, halo_lam, halo_phi,     &
     &                      L_regular, depart_lambda_w, depart_phi_w,   &
     &                      halo_i, halo_j,                             &
     &                      global_row_length,                          &
     &                      row_length, rows, l_datastart,              &
     &                      i_out, j_out,                               &
     &                      weight_lambda, weight_phi)


      If (.not. L_2d_sl_geometry .and.                                  &
     &    (model_domain  ==  mt_Global .or.                             &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) ) Then
        If (Error_Code  ==  0 ) Then

          L_vector = .true.
          type = 1 ! data is at u points.

          number_of_inputs = 1
! Note: dummy arrays, not used, replaced by variable dummy.

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
     &                     row_length, rows, model_levels-1,            &
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
     &                     depart_r_w, depart_lambda_w, depart_phi_w,   &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j,                              &
     &                     global_row_length, global_rows,              &
     &                     row_length, rows, n_rows, rows,              &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     proc_row_group, proc_col_group, 1, 0, 0,     &
     &                     L_sl_halo_reprod, off_x,off_y,               &
     &                     Yu_d, dummy, dummy, Error_Code)

        End If

! ----------------------------------------------------------------------
! Section 2.2 Call interpolation routine for Yv.
! ----------------------------------------------------------------------

        If (Error_Code  ==  0 ) Then

          type = 2 ! data at v points.
          L_vector = .true.

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
     &                     row_length, rows, model_levels-1,            &
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
     &                     cos_theta_latitude, L_regular,               &
     &                     L_vector, model_domain, L_high, L_mono,      &
     &                     .false.,                                     &
     &                     depart_r_w, depart_lambda_w, depart_phi_w,   &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j,                              &
     &                     global_row_length, global_rows,              &
     &                     row_length, rows, n_rows, n_rows,            &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     proc_row_group, proc_col_group, 1, 0, 0,     &
     &                     L_sl_halo_reprod, off_x,off_y,               &
     &                     Yv_d, dummy, dummy, Error_Code)

        End If

      End If ! on 2d geometry

! ----------------------------------------------------------------------
! Section 2.3 Call interpolation routine for Yw.
! ----------------------------------------------------------------------

      If (Error_Code  ==  0 ) Then

        type = 3 ! data at w points.
        L_vector = .false.

        If (model_domain  ==  mt_Global .or.                            &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) Then
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
     &                     row_length, rows, model_levels-1,            &
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
     &                     L_conserv,                                   &
     &                     depart_r_w, depart_lambda_w, depart_phi_w,   &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j,                              &
     &                     global_row_length, global_rows,              &
     &                     row_length, rows, n_rows, rows,              &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     proc_row_group, proc_col_group, 1, 0, 0,     &
     &                     L_sl_halo_reprod, off_x,off_y,               &
     &                     Yw_d, dummy, dummy, Error_Code)

        Else

          number_of_inputs = 2
! DEPENDS ON: interpolation
          Call Interpolation(                                           &
     &                     Yw(1-halo_i,1-halo_j,0,2),                   &
     &                     Yw(1-halo_i,1-halo_j,0,1), dummy,            &
     &                     eta_theta_levels,                            &
     &                     r_theta_levels, delta_r, r_theta_levels,     &
     &                     type, number_of_inputs,                      &
     &                     check_bottom_levels,                         &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level+1,                &
     &                     row_length, rows, model_levels+1,            &
     &                     rows,                                        &
     &                     row_length, rows, model_levels-1,            &
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
     &                     L_conserv,                                   &
     &                     depart_r_w, depart_lambda_w, depart_phi_w,   &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j,                              &
     &                     global_row_length, global_rows,              &
     &                     row_length, rows, n_rows, rows,              &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     proc_row_group, proc_col_group, 1, 0, 0,     &
     &                     L_sl_halo_reprod, off_x,off_y,               &
     &                     Yw_d(1,1,1,1), Yw_d(1,1,1,2),                &
     &                     dummy, Error_Code)

        End If
      End If

! ----------------------------------------------------------------------
! Section 3.  Form R_w from the interpolated terms using
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

          If (L_2d_sl_geometry .or. L_trivial_trigs)Then

            Do k = 1, model_levels-1
              Do j = 1, rows
                Do i = 1, row_length
                  R_w(i,j,k) =  Yw_d(i,j,k,1) + R_w(i,j,k)
                End Do
              End Do
            End Do

          Else

            Do k = 1, model_levels-1
              Do j = 1, rows
                Do i = 1, row_length
                   cos_lambda_ad = cos( lambda_a(i) -                   &
     &                                  depart_lambda_w(i,j,k) )
                  sin_lambda_ad = sin( lambda_a(i) -                    &
     &                                  depart_lambda_w(i,j,k) )
                  cos_phi_d = cos ( depart_phi_w(i,j,k))
                  sin_phi_d = sin ( depart_phi_w(i,j,k))

                  trig_u = sin_lambda_ad * cos_theta_latitude(i,j)
                  trig_v = sin_theta_latitude(i,j) * cos_phi_d -        &
     &                     cos_theta_latitude(i,j) * sin_phi_d *        &
     &                     cos_lambda_ad
                  trig_w = sin_theta_latitude(i,j) * sin_phi_d +        &
     &                     cos_theta_latitude(i,j) * cos_phi_d *        &
     &                     cos_lambda_ad

                  R_w(i,j,k) = trig_u * Yu_d(i,j,k,1)                   &
     &                       + trig_v * Yv_d(i,j,k,1)                   &
     &                       + trig_w * Yw_d(i,j,k,1)                   &
     &                       + R_w(i,j,k)

                End Do
              End Do
            End Do

          End If

        Else If (model_domain  ==  mt_LAM ) Then
! NB: code is the same whether 2d geometry option on or off

          If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then

          Do k = 1, model_levels - 1
! internal area uses real alphas
            Do j = 1, rows
              Do i = 1, row_length
                R_w(i,j,k) = (Yw_d(i,j,k,1) +                           &
     &                       Yw_d(i,j,k,2) * (1. - alpha_4) )           &
     &                     + alpha_4 * Yw(i,j,k,1)                      &
     &                     - w(i,j,k)
              End Do
            End Do
          End Do

! boundary area uses alphas = 1
! southern area
          If (at_extremity(PSouth)) Then
            Do k = 1, model_levels - 1
              Do j = 1, halo_j
                Do i = 1, row_length
                 R_w(i,j,k) =  Yw_d(i,j,k,1)                            &
     &                       + Yw(i,j,k,1)                              &
     &                       - w(i,j,k)
                End Do
              End Do
            End Do
          End If
! Western Area
          If (at_extremity(PWest)) Then
            Do k = 1, model_levels - 1
              Do j = 1, rows
                Do i = 1, halo_i
                 R_w(i,j,k) =  Yw_d(i,j,k,1)                            &
     &                       + Yw(i,j,k,1)                              &
     &                       - w(i,j,k)
                End Do
              End Do
            End Do
          End If
! Eastern Area
          If (at_extremity(PEast)) Then
            Do k = 1, model_levels - 1
              Do j = 1, rows
                Do i = row_length-halo_i+1, row_length
                 R_w(i,j,k) =  Yw_d(i,j,k,1)                            &
     &                       + Yw(i,j,k,1)                              &
     &                       - w(i,j,k)
                End Do
              End Do
            End Do
          End If
! Northern Area
          If (at_extremity(PNorth)) Then
            Do k = 1, model_levels - 1
              Do j = rows-halo_j+1, rows
                Do i = 1, row_length
                 R_w(i,j,k) =  Yw_d(i,j,k,1)                            &
     &                       + Yw(i,j,k,1)                              &
     &                       - w(i,j,k)
                End Do
              End Do
            End Do
          End If

          Else ! if CycleNo > 1 and L_new_tdisc
!
! First do the boundaries to set appropriate loop indices
! thus avoiding overwriting R_u or using temp storage
!
! Loop indices for R_u update when processor is away from boundaries
!
            imin = 1
            imax = row_length
            jmin = 1
            jmax = row_length
!
! boundary area uses alphas = 1
!
! southern area
            If (at_extremity(PSouth)) Then
              Do k = 1, model_levels - 1
                Do j = 1, halo_j
                  Do i = 1, row_length
                   R_w(i,j,k) =  Yw_d(i,j,k,1) + R_w(i,j,k)
                End Do
              End Do
            End Do
          End If
! Western Area
          If (at_extremity(PWest)) Then
            Do k = 1, model_levels - 1
! use jmin to avoid overwriting S boundary
              Do j = jmin, rows
                Do i = 1, halo_i
                 R_w(i,j,k) =  Yw_d(i,j,k,1) + R_w(i,j,k)
                End Do
              End Do
            End Do
          End If
! Eastern Area
          If (at_extremity(PEast)) Then
            Do k = 1, model_levels - 1
! use jmin to avoid overwriting S boundary
              Do j = jmin, rows
                Do i = row_length-halo_i+1, row_length
                 R_w(i,j,k) =  Yw_d(i,j,k,1) + R_w(i,j,k)
                End Do
              End Do
            End Do
          End If
! Northern Area
          If (at_extremity(PNorth)) Then
            Do k = 1, model_levels - 1
              Do j = rows-halo_j+1, rows
! use imin, imax to avoid overwriting E and W boundary
                Do i = imin, imax
                  R_w(i,j,k) =  Yw_d(i,j,k,1) + R_w(i,j,k)
                End Do
              End Do
            End Do
          End If

          Do k = 1, model_levels - 1
! internal area uses real alphas
            Do j = 1, rows
              Do i = 1, row_length
                R_w(i,j,k) = (Yw_d(i,j,k,1) +                           &
     &                       Yw_d(i,j,k,2) * (1. - alpha_4) ) +         &
     &                       R_w(i,j,k)
              End Do
            End Do
          End Do

          End If ! CycleNo ==1 ...

        End If ! on model domain

        k = model_levels
          Do j = 1, rows
            Do i = 1, row_length
              R_w(i,j,k) = 0.
            End Do
          End Do

      End If

! End of routine.
      return
      END SUBROUTINE SL_Vector_w

