
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Calc_non_int_sl_theta
!

      Subroutine Calc_non_int_sl_theta(                                 &
     &                          w, theta, theta_lbcs,                   &
     &                          r_theta_levels, r_rho_levels,           &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          row_length, rows, model_levels,         &
     &                          rimwidth, rimweights, rims_to_do,       &
     &                          lenrim,lbc_size, lbc_start,             &
     &                          n_ext_fields,                           &
     &                          lambda_rm, lambda_rp, phi_rm, phi_rp,   &
     &                          recip_lambda_m, recip_lambda_0,         &
     &                          recip_lambda_p, recip_lambda_p2,        &
     &                          recip_phi_m, recip_phi_0,               &
     &                          recip_phi_p, recip_phi_p2,              &
     &                          i_out_in, j_out_in,                     &
     &                          weight_lambda, weight_phi,              &
     &                          model_domain, timestep, alpha_2,        &
     &                          high_order_scheme_theta,                &
     &                          monotone_scheme_theta,                  &
     &                          L_high_theta, L_mono_theta,             &
     &                          L_sl_halo_reprod, L_lbc_new,            &
     &                          L_regular, L_new_tdisc, CycleNo,        &
     &                          r_out, lambda_out, phi_out,             &
     &                          me, n_procx, n_procy,                   &
     &                          off_x, off_y, halo_i, halo_j,           &
     &                          datastart, g_i_pe, at_extremity,        &
     &                          g_row_length, proc_row_group,           &
     &                          g_rows, proc_col_group,                 &
     &                          theta_star, theta_np1,                  &
     &                          Error_Code)

! Purpose:
!          Performs non-interpolating in the vertical semi-Lagrangian
!          advection of theta.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date     Comment
! ----     -------     -------
! 09/03/00 5.1         Changed code dealing with LAM LBC updating
!                                                          P.Burton
! 24/03/00 5.1        Changes to bring up to v2p9.  Andy Malcolm
!LL   5.1   11/02/00  Use DOMTYP parameters                    P.Burton
! 5.1 14/03/00 Correct occasional failures from uninitialised data
!              accessed above top level when l_mono(sl_theta=1)
!              switched on. R Rawlins
!  5.2  26/09/00   make theta active at top level         Andy Malcolm
! 5.5      12/02/03     fix bug for L_mono_theta         Andy Malcolm
!  6.0  18/08/03  NEC SX-6 optimisation - R Barnes & J-C Rioual.
!  6.1  07/04/04   Reintroduce Correction for LAM unset variables
!                  Andy Malcolm/Paul Selwood
!  6.2  21/10/05   Remove commented out code. P.Selwood.
!  6.2  21/03/05   correct history header              Andy Malcolm
!  6.2  06/10/05   Fix non-reprod code when quintic chosen  Andy Malcolm
! 6.2      04/10/05  Changes for cycling semi-Lagrangian scheme.
!                                                        M. Diamantakis
!  6.2  25/12/05  Variable resolution changes            Yongming Tang
!  6.2  25/12/05  rims_to_do an argument for lbc updating
!                  needed for calc_non_int_sl_theta      Terry Davies
!  6.4  27/11/06  Fix to integer arithmetic bug          A.Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                      ! Dimension of Data_in in i direction.
     &, rows                                                            &
                      ! Dimension of Data_in in j direction.
     &, model_levels                                                    &
                      ! Dimension of Data_in in k direction.
     &, RIMWIDTH                                                        &
                      ! Width of boundaries in LBCs
     &, RIMWEIGHTS(RIMWIDTH)                                            &
                      ! Weights to apply to the LBCs
     &, LENRIM                                                          &
                      ! Size of single level of LBC data
     &, LBC_SIZE(4)                                                     &
                      ! Size of each side of LBC data
     &, LBC_START(4)                                                    &
                      ! Start of each side in LBC data
     &, rims_to_do                                                      &
                            ! rim zone where wt = 1
     &, me                                                              &
                      ! My processor number
     &, n_procx                                                         &
                      ! Number of processors in longitude
     &, n_procy                                                         &
                      ! Number of processors in latitude
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, datastart(3)                                                    &
                      ! First gridpoints held by this processor.
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, g_row_length                                                    &
                      ! global number of points on a row
     &, g_rows                                                          &
                      ! global number of rows
     &, g_i_pe(1-halo_i:g_row_length+halo_i)                            &
                                             ! processor on my procr-row
                             ! holding a given value in i direction
     &, CycleNo

      Logical                                                           &
     &  L_sl_halo_reprod                                                &
                         ! if true then sl code bit repoducible with
                         ! any sensible halo size
     &, L_regular                                                       &
                    ! false if variable resolution
     &, L_lbc_new                                                       &
                    !  true for new lbc treatment
     &, L_new_tdisc

      Integer                                                           &
     &  interp_vertical_search_tol                                      &
                                   !number of levels either side of
                                   ! default level to search.
     &, check_bottom_levels                                             &
                            ! used in interpolation code, and is
                            ! the number of levels to check to see
                            ! if the departure point lies inside the
                            ! orography.
     &, n_ext_fields    ! number of ext_data fields required, 1 for
                        ! global model, 2 for LAM.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Integer                                                           &
     &  high_order_scheme_theta                                         &
                                 ! a code saying which high order
                           ! scheme to use for theta.
     &, monotone_scheme_theta ! a code saying which monotone
                           ! scheme to use for theta.

      Logical                                                           &
     &  L_high_theta                                                    &
                       ! True, if high order interpolation required
                       !       for theta.
     &, L_mono_theta   ! True, if interpolation required to be monotone
                       !       for theta.

      Real                                                              &
     &  timestep                                                        &
     &, alpha_2

!VarRes horizontal co-ordinate information
      Real                                                              &
     &  lambda_rm   ( 1-halo_i : row_length+halo_i )                    &
     &, lambda_rp   ( 1-halo_i : row_length+halo_i )                    &
     &, phi_rm      ( 1-halo_i : row_length + halo_i                    &
     &,               1-halo_j : rows + halo_j )                        &
     &, phi_rp      ( 1-halo_i : row_length + halo_i                    &
     &,               1-halo_j : rows + halo_j )                        &
     &, recip_lambda_m (1-halo_i : row_length+halo_i)                   &
     &, recip_lambda_0 (1-halo_i : row_length+halo_i)                   &
     &, recip_lambda_p (1-halo_i : row_length+halo_i)                   &
     &, recip_lambda_p2(1-halo_i : row_length+halo_i)                   &
     &, recip_phi_m  ( 1-halo_i : row_length + halo_i                   &
     &,                1-halo_j : rows + halo_j )                       &
     &, recip_phi_0  ( 1-halo_i : row_length + halo_i                   &
     &,                1-halo_j : rows + halo_j )                       &
     &, recip_phi_p  ( 1-halo_i : row_length + halo_i                   &
     &,                1-halo_j : rows + halo_j )                       &
     &, recip_phi_p2 ( 1-halo_i : row_length + halo_i                   &
     &,                1-halo_j : rows + halo_j )

      Real                                                              &
     &  r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)

      Real                                                              &
     &  lambda_out (row_length, rows, model_levels)                     &
                                                      ! Lambda
                                                      ! co-ordinate of
                                                      ! output data on
                                                      ! input.
     &, phi_out (row_length, rows, model_levels)      ! Phi Co-ordinate
                                                      ! of output data
                                                      ! on input.
      Real                                                              &
     &  weight_lambda (row_length, rows, model_levels)                  &
     &, weight_phi (row_length, rows, model_levels)

      Integer                                                           &
     &  i_out_in (row_length, rows, model_levels)                       &
     &, j_out_in (row_length, rows, model_levels)

      Real                                                              &
     &  theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, theta_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &         model_levels)                                            &
     &, w (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &     0:model_levels)

      Real                                                              &
     &  THETA_LBCS(LENRIM, MODEL_LEVELS)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent IN/OUT. ie: Input and Output variables.
      Real                                                              &
     &  r_out (row_length, rows, model_levels)        ! Vertical
                                                      ! co-ordinate
                                                      ! of output data.

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  theta_star(1-off_x:row_length+off_x,                            &
     &             1-off_y:rows+off_y, model_levels)

      Integer                                                           &
     &  Error_Code     ! Non-zero on exit if error detected.

! Local Variables.

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
! HIGHOS starts:
! contains the allowed high order scheme options
!
! Author : Michail Diamantakis
! History:
! Version  Date      Comment.
! 5.3      25/10/01  New comdeck
!
      INTEGER,PARAMETER :: cubicLagrange         = 1
      INTEGER,PARAMETER :: quinticLagrange       = 2
      INTEGER,PARAMETER :: ECMWF_quasiCubic      = 3
      INTEGER,PARAMETER :: ECMWF_mono_quasiCubic = 4
      INTEGER,PARAMETER :: hCubic_vLin           = 5
      INTEGER,PARAMETER :: hQuasiCubic_vQuintic  = 6
      INTEGER,PARAMETER :: hCubic_vQuintic       = 7
! HIGHOS ends
! Description: COMDECK containing the allowed
!              monotone scheme options
!
! Author : Michail Diamantakis
! History:
! Version  Date      Comment.
! 5.3      29/10/01  New comdeck
!
      INTEGER                                                           &
     &     triLinear                                                    &
     &,    mono_quasiCubic

      PARAMETER(                                                        &
     &     triLinear       = 1                                          &
     &,    mono_quasiCubic = 2 )

      integer ixl(rows*row_length),jxl(rows*row_length),ijlc
      integer ixu(rows*row_length),jxu(rows*row_length),ijuc
      integer ij
! scalars

      Integer                                                           &
     &  i, j, k                                                         &
                       ! Loop indices
     &, index                                                           &
     &, lower_limit, upper_limit                                        &
     &, type                                                            &
     &, number_of_inputs

      Integer :: ErrorStatus

      Real                                                              &
     &  r_below                                                         &
     &, r_above                                                         &
     &, r_here                                                          &
     &, r_belowi(row_length,rows)                                       &
     &, r_abovei(row_length,rows)                                       &
     &, r_herei(row_length,rows)                                        &
     &, r_mid_above                                                     &
     &, r_mid_below                                                     &
     &, r_mid_abovei(row_length,rows)                                   &
     &, r_mid_belowi(row_length,rows)                                   &
     &, max_mono                                                        &
     &, min_mono

      Logical                                                           &
     &  L_vector                                                        &
     &, L_conserv

! arrays

      Logical                                                           &
     &  L_continue(row_length, rows, model_levels)

      Integer                                                           &
     &  i_out (row_length, rows, model_levels)                          &
     &, j_out (row_length, rows, model_levels)

      Integer                                                           &
     &  depart_level(row_length, rows, model_levels)                    &
                                                     ! closest model lev
                                                    ! to departure point
     &, level_below(row_length, rows, model_levels) ! model level just
                                                 ! below departure point

      Real                                                              &
     &  w_star(row_length, rows, model_levels) ! vertical velocity
                                               ! required to move parcel
                                              ! from closest model level
          ! to departure point to arrival point model level.

      Real                                                              &
     &  ext_r_rho_levels(1-halo_i:row_length+halo_i,                    &
     &                   1-halo_j:rows+halo_j, model_levels+1)

      Real                                                              &
     &  work (row_length, rows, model_levels)                           &
     &, work2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &         model_levels)                                            &
     &, w_minus_wstar (row_length, rows, model_levels)                  &
     &, ext_data (1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,      &
     &            -1:model_levels+2,n_ext_fields )                      &
     &, r_d_s (row_length, rows)                                        &
     &, a_coeff_a (row_length, rows)                                    &
     &, coeff_a                                                         &
     &, coeff_b                                                         &
     &, coeff_c                                                         &
     &, coeff_d                                                         &
     &, coeff_ai(row_length,rows)                                       &
     &, coeff_bi(row_length,rows)                                       &
     &, coeff_ci(row_length,rows)                                       &
     &, coeff_di(row_length,rows)                                       &
     &, theta_d_max(row_length, rows, model_levels)                     &
     &, theta_d_min(row_length, rows, model_levels)

      LOGICAL                                                           &
     &  L_do_halos                                                      &
                          ! update the halos?
     &, L_do_boundaries   ! update the boundaries?

! External Routines:

      External Bi_Linear_niv                                            &
     &,        Cubic_Lagrange_niv                                       &
     &,        Quintic_Lagrange_niv                                     &
     &,        ECMWF_quasi_cubic_niv                                    &
     &,        ECMWF_mono_quasi_cubic_niv                               &
     &,        Swap_Bounds

! Varibles applied in the "compute-on-demand" strategy

      Integer                                                           &
     &  ime, ibase, irecv, my_imin, my_imax, dim_e_out, h_factor        &
     &, nsend, nrecv, info, len, itmp, j0, j1                           &
     &, my_iminp, my_imaxp
      Logical L_continue_e
      Real r_d_s_e(model_levels*g_row_length)

      Integer sp_send(0:n_procx-1), sp_levels(0:n_procx-1,model_levels)
      Integer np_send(0:n_procx-1), np_levels(0:n_procx-1,model_levels)
      Integer kk, sender
      Real ctmp1(2,model_levels)

      Integer                                                           &
     &  n_sendto(0:n_procx-1), n_recvfrom(0:n_procx-1)                  &
     &, i_store(model_levels*g_row_length,0:n_procx-1)                  &
     &, j_store(model_levels*g_row_length,0:n_procx-1)                  &
     &, k_store(model_levels*g_row_length,0:n_procx-1)                  &
     &, i_out_e(model_levels*g_row_length)                              &
     &, j_out_e(model_levels*g_row_length)                              &
     &, depart_level_e(model_levels*g_row_length)                       &
     &, level_below_e(model_levels*g_row_length)                        &
     &, k_e(model_levels*g_row_length)                                  &
     &, isend_arr(3,model_levels*g_row_length,0:n_procx-1)              &
     &, irecv_arr(3,model_levels*g_row_length,0:n_procx-1)
      Real                                                              &
     &  rsend_arr(5,model_levels*g_row_length,0:n_procx-1)              &
     &, rrecv_arr(5,model_levels*g_row_length,0:n_procx-1)              &
     &, weight_lambda_e(model_levels*g_row_length)                      &
     &, weight_phi_e(model_levels*g_row_length)                         &
     &, coeff_a_e(model_levels*g_row_length)                            &
     &, coeff_b_e(model_levels*g_row_length)                            &
     &, coeff_c_e(model_levels*g_row_length)                            &
     &, coeff_d_e(model_levels*g_row_length)                            &
     &, r_out_e(model_levels*g_row_length)                              &
     &, w_e(model_levels*g_row_length)                                  &
     &, w_star_e(model_levels*g_row_length)                             &
     &, theta_d_max_e(model_levels*g_row_length)                        &
     &, theta_d_min_e(model_levels*g_row_length)                        &
     &, theta_star_e(model_levels*g_row_length)                         &
     &, r_theta_levels_e(model_levels*g_row_length)                     &
     &, send_data(model_levels*g_row_length,0:n_procx-1)                &
     &, recv_data(model_levels*g_row_length,0:n_procx-1)                &
     &, rsend_arr2(2,model_levels*g_row_length)                         &
     &, rrecv_arr2(2,model_levels*g_row_length,0:n_procx-1)

! Functions: None

! ----------------------------------------------------------------------
! Section 1.   Set some control variables.
! ----------------------------------------------------------------------

      type = 3 ! data at theta points
      L_vector = .false. ! field to be interpolated is not a horizontal
!                          vector component.
      L_conserv = .false. ! conservation not possible with this scheme
      number_of_inputs = 1 ! only one field to interpolate

! Execute rest of routine only if error code is still zero.

      If (Error_Code  ==  0 ) Then

! ----------------------------------------------------------------------
! Section 2.   Extend r array.
! ----------------------------------------------------------------------

      Do k = 1, model_levels

        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            Ext_r_rho_levels (i,j,k) = r_rho_levels(i,j,k)
          End Do
        End Do

      End Do

! add extra rho level just above top theta level.

        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            Ext_r_rho_levels(i,j,model_levels+1) =                      &
     &                  r_theta_levels(i,j,model_levels) + 5.0
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 3.   For each output point find i,j so that the point on the
!              output grid lies between i and i+1, j and j+1.
! ----------------------------------------------------------------------

! i_out and j_out must be copied because they change in this subroutine
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            i_out(i,j,k) = i_out_in(i,j,k)
            j_out(i,j,k) = j_out_in(i,j,k)
          End Do
        End Do
      End Do

        If (n_procx > 1 ) then
        If (model_domain == mt_Global ) then
! Send the points outside my region to the appropriate processor for
! interpolation. Only performed if the domain is decomposed in the
! i direction.

! The first and last point I can interpolate in, based on available
! data on this processor

          h_factor = 2
          If ( high_order_scheme_theta  ==  quinticLagrange .and.       &
     &         L_High_theta )                                           &
     &         h_factor = 3
          my_imin = datastart(1) - halo_i + h_factor
          my_imax = datastart(1) + row_length - 1 + halo_i - h_factor

! values for use in polar row to ensure pole is only calculated on one
! processor.
          my_iminp = datastart(1)
          my_imaxp = datastart(1)+row_length-1
          If(at_extremity(PWest) ) my_iminp = my_imin
          If(at_extremity(PEast) ) my_imaxp = my_imax

! The base processor on this row, and my address relative to that
! processor

          ibase = (me/n_procx) * n_procx
          ime = me - ibase

          Do i = 0, n_procx-1
            n_sendto(i) = 0
          End Do

          Do i = 0, n_procx-1
            sp_send(i) = 0
          End Do
          If (at_extremity(PSouth)) then
            Do k = 1, model_levels
                If (i_out(1,1,k)  >=  my_iminp .and.                    &
     &              i_out(1,1,k)  <=  my_imaxp)  then
                i_out(1,1,k) = i_out(1,1,k) - datastart(1) + 1
                sp_send(ime) = sp_send(ime) + 1
                sp_levels(ime,sp_send(ime)) = k
                Do i = 2, row_length
                  i_out(i,1,k) = i ! i_out(i,1,k) - datastart(1) + 1
                End Do
              Else
                sender = g_i_pe(i_out(1,1,k))
                sp_send(sender) = sp_send(sender) + 1
                sp_levels(sender,sp_send(sender)) = k
                Do i = 1, row_length
                  i_out(i,1,k) = i
                End Do
              End If
            End Do
          End If
          Do i = 0, n_procx-1
            np_send(i) = 0
          End Do
          If (at_extremity(PNorth)) then
            Do k = 1, model_levels
                If (i_out(1,rows,k)  >=  my_iminp .and.                 &
     &              i_out(1,rows,k)  <=  my_imaxp)  then
                np_send(ime) = np_send(ime) + 1
                np_levels(ime,np_send(ime)) = k
                i_out(1,rows,k) = i_out(1,rows,k)-datastart(1)+1
                Do i = 2, row_length
                  i_out(i,rows,k) = i ! i_out(i,rows,k)-datastart(1)+1
                End Do
              Else
                sender = g_i_pe(i_out(1,rows,k))
                np_send(sender) = np_send(sender) + 1
                np_levels(sender,np_send(sender)) = k
                Do i = 1, row_length
                  i_out(i,rows,k) = i
                End Do
              End If
            End Do
          End If

          j0 = 1
          j1 = rows
          If (at_extremity(PSouth)) j0 = 2
          If (at_extremity(PNorth)) j1 = rows-1

          If( L_sl_halo_reprod) Then

! On the global boundaries, use i_out < 1 or i_out > g_row_length
! if that makes local computation possible. Not required when
! L_sl_halo_reprod is false is other logic ensures this is done.

! This code unsafe if applied at poles, where it isn't required.

          If (at_extremity(PWest)) then
            Do k = 1, model_levels
              Do j = j0, j1
                Do i = 1, halo_i
                  If (i_out(i,j,k)  >   g_row_length-halo_i+h_factor)   &
     &                 i_out(i,j,k) = i_out(i,j,k) - g_row_length
                End Do
              End Do
            End Do
          End If
          If (at_extremity(PEast)) then
            Do k = 1, model_levels
              Do j = j0, j1
                Do i = row_length-halo_i+1, row_length
                  If (i_out(i,j,k)  <   halo_i-h_factor+1) then
                    i_out(i,j,k) = i_out(i,j,k) + g_row_length
                  End If
                End Do
              End Do
            End Do
          End If

          End If ! on L_sl_halo_reprod

          Do k = 1, model_levels
            Do j = j0, j1
              Do i = 1, row_length
                If (i_out(i,j,k)  >=  my_imin .and.                     &
     &               i_out(i,j,k)  <=  my_imax) then
! Process locally, so find the local destination
                  i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
                Else
!     CODE TO STOP BIT NON-REPRODUCIBILITY
                  if(i_out(i,j,k) > g_row_length+halo_i-h_factor)then
                    i_out(i,j,k)=i_out(i,j,k)-g_row_length
                  endif
                  if(i_out(i,j,k) < 1-halo_i+h_factor)then
                    i_out(i,j,k)=i_out(i,j,k)+g_row_length
                  endif
!     END CODE TO STOP BIT NON-REPRODUCIBILITY
! Send to a remote processor, given by the array g_i_pe
                  irecv = g_i_pe(i_out(i,j,k))
                  n_sendto(irecv) = n_sendto(irecv) + 1
                  itmp = n_sendto(irecv)
                  isend_arr(1,itmp,irecv) = i_out(i,j,k)
                  isend_arr(2,itmp,irecv) = j_out(i,j,k)
                  isend_arr(3,itmp,irecv) = k
                  rsend_arr(1,itmp,irecv) = weight_lambda(i,j,k)
                  rsend_arr(2,itmp,irecv) = weight_phi(i,j,k)
                  rsend_arr(3,itmp,irecv) = r_out(i,j,k)
                  rsend_arr(4,itmp,irecv) = r_theta_levels(i,j,k)
                  rsend_arr(5,itmp,irecv) = w(i,j,k)
                  i_store(itmp,irecv) = i
                  j_store(itmp,irecv) = j
                  k_store(itmp,irecv) = k
                  i_out(i,j,k) = i
                End If
              End Do
            End Do
          End Do

          nsend = 0
          Do i = 0,n_procx-1
            call gc_isend(10*(me+1)+ibase+i, 1, ibase+i, info,          &
     &           n_recvfrom(ime), n_sendto(i))
            nsend = nsend + n_sendto(i)
          End Do

          Call gcg_ssync(proc_row_group, info)

          nrecv = 0
          Do i = 0, n_procx-1
            call gc_irecv(10*(ibase+i+1)+me, 1, ibase+i, info,          &
     &           n_recvfrom(i), n_sendto(ime))
            nrecv = nrecv + n_recvfrom(i)
          End Do

!          Call gcg_ssync(proc_row_group, info)

          Do i = 0,n_procx-1
            len = 5*n_sendto(i)
            If (n_sendto(i)  >   0) then
              call gc_rsend(20*(me+1)+ibase+i, len, ibase+i, info,      &
     &             rrecv_arr(1,1,ime), rsend_arr(1,1,i))
            End If
          End Do

          Call gcg_ssync(proc_row_group, info)

          Do i = 0,n_procx-1
            len = 5*n_recvfrom(i)
            If (n_recvfrom(i)  >   0) then
              call gc_rrecv(20*(ibase+i+1)+me, len, ibase+i, info,      &
     &             rrecv_arr(1,1,i), rsend_arr(1,1,ime))
            End If
          End Do

!          Call gcg_ssync(proc_row_group, info)

          Do i = 0,n_procx-1
            len = 3*n_sendto(i)
            If (n_sendto(i)  >   0) then
              call gc_isend(30*(me+1)+ibase+i, len, ibase+i, info,      &
     &             irecv_arr(1,1,ime), isend_arr(1,1,i))
            End If
          End Do

          Call gcg_ssync(proc_row_group, info)

          Do i = 0,n_procx-1
            len = 3*n_recvfrom(i)
            If (n_recvfrom(i)  >   0) then
              call gc_irecv(30*(ibase+i+1)+me, len, ibase+i, info,      &
     &             irecv_arr(1,1,i), isend_arr(1,1,ime))
            End If
          End Do

        Else   ! model is a type of LAM

          k = 1
          j = 1
          Do i = 1, row_length * rows * model_levels
            i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
          End Do  ! i = row_length * rows * model_levels

        EndIf  !  model_domain == mt_Global

        End If  ! n_procx > 1

!     CODE TO STOP BIT NON-REPRODUCIBILITY
        if(model_domain == mt_global .and. n_procx == 1)then
          h_factor = 2
          If (high_order_scheme_theta == quinticLagrange .and.          &
     &          L_High_theta) Then
            h_factor = 3
          End If
          my_imin = datastart(1) - halo_i + h_factor - 1
          my_imax =datastart(1) + row_length - 1 + halo_i - h_factor + 1
          Do k = 1, model_levels
            Do j = 1,rows
              Do i = 1, row_length
                if(i_out(i,j,k) >= my_imax)then
                  i_out(i,j,k)=i_out(i,j,k) - g_row_length
                endif
                if(i_out(i,j,k) <= my_imin )then
                  i_out(i,j,k)=i_out(i,j,k)+g_row_length
                endif
              End Do
            End Do
          End Do
        endif   !model_domain == 1 .and. n_procx == 1
!     END CODE TO STOP BIT NON-REPRODUCIBILITY

! ----------------------------------------------------------------------
! Section 4.   Find closest model level to departure point and
!              calculate w_star.
!              Limit trajectories so that they do not go below bottom
!              data level.
!              Find max/min values of theta that surround departure
!              point. (only if L_theta_mono)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 4.1  W-star and trajectory limits for points on processor
! ----------------------------------------------------------------------

        Do k = 1, model_levels

! calculate horizontal interpolation weights
          Do j = 1, rows
            Do i = 1, row_length
              depart_level(i,j,k) = k
              w_star(i,j,k) = 0.
              L_continue(i,j,k) = .true.
            End Do
          End Do

          If ( k  <=  check_bottom_levels ) Then
! Perform check for below bottom data surface.

           Do j = 1, rows
              Do i = 1, row_length

              coeff_a = (1.-weight_lambda(i,j,k)) *                     &
     &                       (1.-weight_phi(i,j,k))
              coeff_b = weight_lambda(i,j,k) *                          &
     &                       (1.-weight_phi(i,j,k))
              coeff_c = (1.-weight_lambda(i,j,k)) *                     &
     &                       weight_phi(i,j,k)
              coeff_d = weight_lambda(i,j,k)*weight_phi(i,j,k)

                r_d_s(i,j) = coeff_a *                                  &
     &              r_theta_levels (i_out(i,j,k),j_out(i,j,k),1)        &
     &              + coeff_b                                           &
     &              * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k),1)    &
     &              + coeff_c                                           &
     &              * r_theta_levels (i_out(i,j,k),j_out(i,j,k)+1,1)    &
     &              + coeff_d *                                         &
     &              r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)+1,1)

                If (r_out(i,j,k) <   r_d_s(i,j) ) Then
! move trajectory up to lowest level of data.
                  r_out(i,j,k) = r_d_s(i,j)
! Set logical switch to say it was done
                  L_continue(i,j,k) = .false.
! set w_star = w to turn off vertical adjustment term
                  w_star(i,j,k) = w(i,j,k)
! depart_level = 1 to get bottom most value.
                  depart_level(i,j,k) = 1
                End If

              End Do
            End Do

          End If

! Find k point.
! use search over restricted levels.
! Find level which is just below r_out value, min possible is
! max(1,k- interp_vertical_search_tol), max possible is
! min(model_levels, k+interp_vertical_search_tol)-1

          If ( k  <=  check_bottom_levels ) Then
            lower_limit = max(1,k - check_bottom_levels)
            upper_limit = min(model_levels,                             &
     &                        k + check_bottom_levels)
          Else
            lower_limit = max(1,k - interp_vertical_search_tol)
            upper_limit = min(model_levels,                             &
     &                        k + interp_vertical_search_tol)
          End If

          If (k  ==  1) Then
! level 1 Only performs level check and upward search
            Do j = 1, rows
              Do i = 1, row_length
                coeff_a = (1.-weight_lambda(i,j,k)) *                   &
     &                       (1.-weight_phi(i,j,k))
                coeff_b = weight_lambda(i,j,k) *                        &
     &                       (1.-weight_phi(i,j,k))
                coeff_c = (1.-weight_lambda(i,j,k)) *                   &
     &                       weight_phi(i,j,k)
                coeff_d = weight_lambda(i,j,k)*weight_phi(i,j,k)

                If (L_continue(i,j,k)) Then
                  index = k
                  r_above =                                             &
     &   coeff_a * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)  ,index+1)&
     & + coeff_b * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)  ,index+1)&
     & + coeff_c * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)+1,index+1)&
     & + coeff_d * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)+1,index+1)

                  r_here = r_d_s(i,j)
                  r_mid_above= (r_above + r_here )/2.

                  If (r_out(i,j,k)  <   r_mid_above ) Then
! No need to check r_below as bottom adjust has done this.
! set w_star

                    w_star(i,j,k) = (r_theta_levels(i,j,k) - r_d_s(i,j))&
     &                               / timestep

                  Else
! upward search
                    Do while (index  <   upper_limit .and.              &
     &                      L_continue(i,j,k) )
                      index = index + 1

                      r_here = r_above
                      r_mid_below= r_mid_above
                      r_above =                                         &
     &   coeff_a * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)  ,index+1)&
     & + coeff_b * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)  ,index+1)&
     & + coeff_c * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)+1,index+1)&
     & + coeff_d * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)+1,index+1)

              r_mid_above= (r_here + r_above)/2.
                      If (r_out(i,j,k)  >=  r_mid_below .and.           &
     &                    r_out(i,j,k)  <   r_mid_above ) Then

                        depart_level(i,j,k) = index
                        w_star(i,j,k) = (r_theta_levels(i,j,k) - r_here)&
     &                                / timestep

                        L_continue(i,j,k) = .false.
                      End if
                    End Do

                  End If ! on level search

                End If ! on bottom adjustment
              End Do
            End Do

          Else If (k  ==  model_levels) Then
! AT upper limit data stays at top level so no work required.
! w_star and depart_level are already set
            Do j = 1, rows
              Do i = 1, row_length
                L_continue(i,j,k) = .false.
              End Do
            End Do

          Else
! Interior levels
! Calculate level k value
            ijlc=0
            ijuc=0
            Do j = 1, rows
              Do i = 1, row_length
                coeff_a = (1.-weight_lambda(i,j,k)) *                   &
     &                       (1.-weight_phi(i,j,k))
                coeff_b = weight_lambda(i,j,k) *                        &
     &                       (1.-weight_phi(i,j,k))
                coeff_c = (1.-weight_lambda(i,j,k)) *                   &
     &                       weight_phi(i,j,k)
                coeff_d = weight_lambda(i,j,k)*weight_phi(i,j,k)

                If (L_continue(i,j,k) ) Then
                  r_above =                                             &
     &   coeff_a * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)  ,k+1)   &
     & + coeff_b * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)  ,k+1)   &
     & + coeff_c * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)+1,k+1)   &
     & + coeff_d * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)+1,k+1)
                  r_here =                                              &
     &   coeff_a * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)  ,k)     &
     & + coeff_b * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)  ,k)     &
     & + coeff_c * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)+1,k)     &
     & + coeff_d * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)+1,k)
                  r_below =                                             &
     &   coeff_a * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)  ,k-1)   &
     & + coeff_b * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)  ,k-1)   &
     & + coeff_c * r_theta_levels (i_out(i,j,k)  ,j_out(i,j,k)+1,k-1)   &
     & + coeff_d * r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k)+1,k-1)

              coeff_ai(i,j) = coeff_a
              coeff_bi(i,j) = coeff_b
              coeff_ci(i,j) = coeff_c
              coeff_di(i,j) = coeff_d
              r_belowi(i,j) = r_below
              r_herei(i,j) = r_here
              r_abovei(i,j) = r_above

              r_mid_abovei(i,j) = (r_above+r_here)/2.
              r_mid_belowi(i,j) = (r_below+r_here)/2.
              endif
              enddo
            enddo

            Do j = 1, rows
              Do i = 1, row_length
                If (L_continue(i,j,k) ) Then
                  r_here = r_herei(i,j)
                  r_mid_above = r_mid_abovei(i,j)
                  r_mid_below = r_mid_belowi(i,j)

                  If (r_out(i,j,k)  >   r_mid_above ) Then
                  ijlc=ijlc+1
                  ixl(ijlc)=i
                  jxl(ijlc)=j

                  Else If (r_out(i,j,k)  <   r_mid_below ) Then
                  ijuc=ijuc+1
                  ixu(ijuc)=i
                  jxu(ijuc)=j


                  Else
! set w_star values as level index is correct
                    w_star(i,j,k) = (r_theta_levels(i,j,k) -            &
     &                           r_here                                 &
     &                              ) / timestep
                  End If

                End If ! end if on bottom adjustment
              End Do
            End Do

            do ij=1,ijlc
                  i=ixl(ij)
                  j=jxl(ij)
                  coeff_a = coeff_ai(i,j)
                  coeff_b = coeff_bi(i,j)
                  coeff_c = coeff_ci(i,j)
                  coeff_d = coeff_di(i,j)
                  r_below = r_belowi(i,j)
                  r_here = r_herei(i,j)
                  r_above = r_abovei(i,j)
                  r_mid_above = r_mid_abovei(i,j)
                  r_mid_below = r_mid_belowi(i,j)

                  index = k
! upward search
                    Do while (index  <   upper_limit .and.              &
     &                        L_continue(i,j,k) )
                      index = index + 1

                      r_below     = r_here
                      r_here      = r_above
                      r_mid_below = r_mid_above
                      r_above =                                         &
     &   coeff_a * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)  ,index+1)&
     & + coeff_b * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)  ,index+1)&
     & + coeff_c * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)+1,index+1)&
     & + coeff_d * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)+1,index+1)
                      r_mid_above=(r_above+r_here)/2.
                      If (r_out(i,j,k)  >=  r_mid_below .and.           &
     &                    r_out(i,j,k)  <   r_mid_above ) Then

                        depart_level(i,j,k) = index
                        w_star(i,j,k) = (r_theta_levels(i,j,k) - r_here)&
     &                               / timestep

                        L_continue(i,j,k) = .false.
                      End if
                    End Do
                 End Do

            do ij=1,ijuc
                  i=ixu(ij)
                  j=jxu(ij)
                  coeff_a = coeff_ai(i,j)
                  coeff_b = coeff_bi(i,j)
                  coeff_c = coeff_ci(i,j)
                  coeff_d = coeff_di(i,j)
                  r_below = r_belowi(i,j)
                  r_here = r_herei(i,j)
                  r_above = r_abovei(i,j)
                  r_mid_above = r_mid_abovei(i,j)
                  r_mid_below = r_mid_belowi(i,j)

                  index = k

                    Do while (index  >   lower_limit .and.              &
     &                        L_continue(i,j,k))
                      index = index - 1
                      r_above = r_here
                      r_here= r_below
                      r_mid_above=r_mid_below
                      r_below =                                         &
     &   coeff_a * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)  ,index-1)&
     & + coeff_b * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)  ,index-1)&
     & + coeff_c * r_theta_levels(i_out(i,j,k)  ,j_out(i,j,k)+1,index-1)&
     & + coeff_d * r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)+1,index-1)
                      r_mid_below=(r_below+r_here)/2.

                      If (r_out(i,j,k)  >=  r_mid_below .and.           &
     &                    r_out(i,j,k)  <   r_mid_above ) Then

                        depart_level(i,j,k) = index
                        w_star(i,j,k) = (r_theta_levels(i,j,k) - r_here)&
     &                                / timestep

                        L_continue(i,j,k) = .false.
                      End if
                    End Do
                 End Do

          End If ! on which model level

! ----------------------------------------------------------------------
! Section 4.2  max/min values for points on processor
! ----------------------------------------------------------------------

! if monotonocity required then find level just below departure point
! depart_level must be either just below or just above desired level
! so work out which it is.
          If (L_mono_theta ) Then
            Do j = 1, rows
              Do i = 1, row_length
                coeff_a = (1.-weight_lambda(i,j,k)) *                   &
     &                       (1.-weight_phi(i,j,k))
                coeff_b =      weight_lambda(i,j,k) *                   &
     &                       (1.-weight_phi(i,j,k))
                coeff_c = (1.-weight_lambda(i,j,k)) *                   &
     &                            weight_phi(i,j,k)
                coeff_d = weight_lambda(i,j,k)*weight_phi(i,j,k)

                r_below = coeff_a *                                     &
     & r_theta_levels(i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k))    &
     &                  + coeff_b *                                     &
     & r_theta_levels (i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k)) &
     &                + coeff_c *                                       &
     & r_theta_levels (i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k)) &
     &                  + coeff_d *                                     &
     & r_theta_levels(i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k))
                If(r_below  <   r_out(i,j,k) .and.                      &
     &             depart_level(i,j,k)  /=  1) Then
                  level_below(i,j,k) = depart_level(i,j,k)
                Else
                  level_below(i,j,k) = depart_level(i,j,k) -1
                End If
              End Do
            End Do
          End If  !  L_mono_theta

        End Do ! k = 1, model_levels


        If (L_mono_theta ) Then
! set up theta in ext_data
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                Ext_data (i,j,k,1) = theta(i,j,k)
              End Do
            End Do
          End Do

! Ensure data above the top level is initialised for later max/min
          Do j = 1,rows
            Do i = 1,row_length
              Ext_data (i,j,model_levels+1,1) = theta(i,j,model_levels)
            End Do
          End Do

! Ensure data below the bottom level is initialised for later max/min
          Do j = 1,rows
            Do i = 1,row_length
              Ext_data (i,j,0,1) = theta(i,j,1)
            End Do
          End Do
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   ext_data(1-halo_i,1-halo_j,0,1),               &
     &                   row_length, rows, model_levels+2,              &
     &                   halo_i, halo_j, fld_type_p, L_vector)

          If (model_domain  ==  mt_LAM) Then

! Set data on edge processor edge haloes from lateral boundary data
            L_do_halos=.TRUE.
            L_do_boundaries=.FALSE.

! DEPENDS ON: set_lateral_boundaries
            CALL SET_LATERAL_BOUNDARIES(                                &
     &        ROW_LENGTH,ROWS,halo_i,halo_j,                            &
     &        MODEL_LEVELS,fld_type_p,                                  &
     &        ext_data(1-halo_i,1-halo_j,1,1),                          &
     &        LENRIM,LBC_SIZE,LBC_START,HALO_I,HALO_J,THETA_LBCS,       &
     &        RIMWIDTH,RIMWIDTH,RIMWEIGHTS,AT_EXTREMITY,                &
     &        L_do_boundaries,L_do_halos)

! Ensure data below the bottom level is initialised for later max/min
! DEPENDS ON: set_lateral_boundaries
            CALL SET_LATERAL_BOUNDARIES(                                &
     &        ROW_LENGTH,ROWS,halo_i,halo_j,                            &
     &        1,fld_type_p,                                             &
     &        ext_data(1-halo_i,1-halo_j,0,1),                          &
     &        LENRIM,LBC_SIZE,LBC_START,HALO_I,HALO_J,THETA_LBCS,       &
     &        RIMWIDTH,RIMWIDTH,RIMWEIGHTS,AT_EXTREMITY,                &
     &        L_do_boundaries,L_do_halos)

! Ensure data above the top level is initialised for later max/min
! DEPENDS ON: set_lateral_boundaries
            CALL SET_LATERAL_BOUNDARIES(                                &
     &        ROW_LENGTH,ROWS,halo_i,halo_j,                            &
     &        1,fld_type_p,                                             &
     &        ext_data(1-halo_i,1-halo_j,model_levels+1,1),             &
     &        LENRIM,LBC_SIZE,LBC_START,HALO_I,HALO_J,                  &
     &        THETA_LBCS(1,model_levels),                               &
     &        RIMWIDTH,RIMWIDTH,RIMWEIGHTS,AT_EXTREMITY,                &
     &        L_do_boundaries,L_do_halos)

          End If ! model_domain == mt_LAM

! calculate max/min theta values at departure points
          Do k = 1, model_levels - 1

            Do j = 1, rows
              Do i = 1, row_length
                theta_d_max(i,j,k) = max (                              &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k),level_below(i,j,k),1),       &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),level_below(i,j,k),1),     &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,level_below(i,j,k),1),     &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,level_below(i,j,k),1),   &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k),level_below(i,j,k)+1,1),     &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),level_below(i,j,k)+1,1),   &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,level_below(i,j,k)+1,1),   &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,level_below(i,j,k)+1,1) )
                theta_d_min(i,j,k) = min (                              &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k),level_below(i,j,k),1),       &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),level_below(i,j,k),1),     &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,level_below(i,j,k),1),     &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,level_below(i,j,k),1),   &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k),level_below(i,j,k)+1,1),     &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),level_below(i,j,k)+1,1),   &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,level_below(i,j,k)+1,1),   &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,level_below(i,j,k)+1,1) )
              End Do
            End Do

          End Do  !  k = 1, model_levels - 1

        End If !  L_mono_theta

! ----------------------------------------------------------------------
! Section 4.3  W-star and trajectory limits for points requested from
!              other processors
! ----------------------------------------------------------------------

! receive extra points and process them

        If (n_procx  >   1 .and. model_domain  ==  mt_Global) then

          dim_e_out = 0
          Do i = 0, n_procx-1
            if (n_recvfrom(i)  >   0) then
              do j = 1, n_recvfrom(i)
                dim_e_out = dim_e_out + 1
                i_out_e(dim_e_out) =                                    &
     &                   irecv_arr(1,j,i) - datastart(1) + 1
                j_out_e(dim_e_out) = irecv_arr(2,j,i)
                k_e(dim_e_out) = irecv_arr(3,j,i)
                weight_lambda_e(dim_e_out) = rrecv_arr(1,j,i)
                weight_phi_e(dim_e_out) = rrecv_arr(2,j,i)
                r_out_e(dim_e_out) = rrecv_arr(3,j,i)
                r_theta_levels_e(dim_e_out) = rrecv_arr(4,j,i)
                w_e(dim_e_out) = rrecv_arr(5,j,i)
              enddo
            endif
          End Do

          If (dim_e_out  >  model_levels*g_row_length) Then
            ErrorStatus = 10
! DEPENDS ON: Ereport
            Call Ereport("calc_non_int_sl_theta", ErrorStatus,          &
     &           "over-writing due to dim_e_out size" )
          End If

          Do i = 1, dim_e_out

! calculate horizontal interpolation weights
            coeff_a_e(i) = (1.-weight_lambda_e(i)) *                    &
     &                    (1.-weight_phi_e(i))
            coeff_b_e(i) = weight_lambda_e(i) *                         &
     &                    (1.-weight_phi_e(i))
            coeff_c_e(i) = (1.-weight_lambda_e(i)) *                    &
     &                       weight_phi_e(i)
            coeff_d_e(i) = weight_lambda_e(i)*weight_phi_e(i)
            depart_level_e(i) = k_e(i)
            w_star_e(i) = 0.
            L_continue_e = .true.

            If ( k_e(i)  <=  check_bottom_levels ) Then
! Perform check for below bottom data surface.

              r_d_s_e(i) =                                              &
     &    coeff_a_e(i) * r_theta_levels(i_out_e(i)  ,j_out_e(i)  ,1)    &
     &  + coeff_b_e(i) * r_theta_levels(i_out_e(i)+1,j_out_e(i)  ,1)    &
     &  + coeff_c_e(i) * r_theta_levels(i_out_e(i)  ,j_out_e(i)+1,1)    &
     &  + coeff_d_e(i) * r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,1)

              If (r_out_e(i) <   r_d_s_e(i) ) Then
! move trajectory up to lowest level of data.
                r_out_e(i) = r_d_s_e(i)
! Set logical switch to say it was done
                L_continue_e = .false.
! set w_star = w to turn off vertical adjustment term
                w_star_e(i) = w_e(i)
! depart_level = 1 to get bottom most value.
                depart_level_e(i) = 1
              End If

            End If

! Find k point.
! use search over restricted levels.
! Find level which is just below r_out value, min possible is
! max(1,k- interp_vertical_search_tol), max possible is
! min(model_levels, k+interp_vertical_search_tol)-1

            If ( k_e(i)  <=  check_bottom_levels ) Then
              lower_limit = max(1,k_e(i) - check_bottom_levels)
              upper_limit = min(model_levels,                           &
     &                          k_e(i) + check_bottom_levels)
            Else
              lower_limit = max(1,k_e(i)-interp_vertical_search_tol)
              upper_limit = min(model_levels,                           &
     &                         k_e(i) + interp_vertical_search_tol)
            End If

            If (k_e(i)  ==  1) Then
! level 1 Only performs level check and upward search
              If (L_continue_e) Then
                index = k_e(i)
                r_above =                                               &
     &   coeff_a_e(i) *r_theta_levels(i_out_e(i)  ,j_out_e(i)  ,index+1)&
     & + coeff_b_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)  ,index+1)&
     & + coeff_c_e(i) *r_theta_levels(i_out_e(i)  ,j_out_e(i)+1,index+1)&
     & + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index+1)
                r_here=r_d_s_e(i)
                r_mid_above= (r_above +r_here)/2.

                If (r_out_e(i)  <   r_mid_above ) Then
! No need to check r_below as bottom adjust has done this.
! set w_star
                  w_star_e(i) = (r_theta_levels_e(i) - r_d_s_e(i) )     &
     &                              / timestep

                Else
! upward search
                  Do while (index  <   upper_limit .and.                &
     &                      L_continue_e )
                    index = index + 1

                    r_here=r_above
                    r_mid_below = r_mid_above
                    r_above =                                           &
     &  coeff_a_e(i) *r_theta_levels(i_out_e(i)  ,j_out_e(i)  ,index+1) &
     & +coeff_b_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)  ,index+1) &
     & +coeff_c_e(i) *r_theta_levels(i_out_e(i)  ,j_out_e(i)+1,index+1) &
     & +coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index+1)
                    r_mid_above=(r_here+r_above)/2.

                    If (r_out_e(i)  >=  r_mid_below .and.               &
     &                  r_out_e(i)  <   r_mid_above ) Then

                      depart_level_e(i) = index
                      w_star_e(i) = (r_theta_levels_e(i) - r_here )     &
     &                               / timestep

                      L_continue_e = .false.
                    End if
                  End Do

                End If ! on level search

              End If ! on bottom adjustment

            Else If (k_e(i)  ==  model_levels) Then
! AT upper limit data stays at top level so no work required.
! w_star and depart_level are already set
              L_continue_e = .false.
              depart_level_e(i) = k_e(i)
              w_star_e(i) = 0.0

            Else
! Interior levels
! Calculate level k value
              If (L_continue_e ) Then
                index = k_e(i)
                r_above =                                               &
     &   coeff_a_e(i) * r_theta_levels (i_out_e(i)  ,j_out_e(i),index+1)&
     & + coeff_b_e(i) * r_theta_levels (i_out_e(i)+1,j_out_e(i),index+1)&
     & + coeff_c_e(i) * r_theta_levels (i_out_e(i),j_out_e(i)+1,index+1)&
     & + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index+1)
                r_here =                                                &
     &   coeff_a_e(i) * r_theta_levels (i_out_e(i)  ,j_out_e(i),index)  &
     & + coeff_b_e(i) * r_theta_levels (i_out_e(i)+1,j_out_e(i),index)  &
     & + coeff_c_e(i) * r_theta_levels (i_out_e(i),j_out_e(i)+1,index)  &
     & + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index)
                r_below =                                               &
     &   coeff_a_e(i) * r_theta_levels (i_out_e(i)  ,j_out_e(i),index-1)&
     & + coeff_b_e(i) * r_theta_levels (i_out_e(i)+1,j_out_e(i),index-1)&
     & + coeff_c_e(i) * r_theta_levels (i_out_e(i),j_out_e(i)+1,index-1)&
     & + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index-1)
                r_mid_above=(r_here+r_above)/2.
                r_mid_below=(r_here+r_below)/2.

                If (r_out_e(i)  >   r_mid_above ) Then
! upward search
                  Do while (index  <   upper_limit .and.                &
     &                      L_continue_e )
                    index = index + 1

                    r_below = r_here
                    r_here  = r_above
                    r_mid_below = r_mid_above
                    r_above =                                           &
     &   coeff_a_e(i) * r_theta_levels (i_out_e(i)  ,j_out_e(i),index+1)&
     & + coeff_b_e(i) * r_theta_levels (i_out_e(i)+1,j_out_e(i),index+1)&
     & + coeff_c_e(i) * r_theta_levels (i_out_e(i),j_out_e(i)+1,index+1)&
     & + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index+1)
                    r_mid_above=(r_here+r_above)/2.

                    If (r_out_e(i)  >=  r_mid_below .and.               &
     &                  r_out_e(i)  <   r_mid_above ) Then

                      depart_level_e(i) = index
                      w_star_e(i) = (r_theta_levels_e(i) - r_here )     &
     &                               / timestep

                      L_continue_e = .false.
                    End if
                  End Do

                Else If (r_out_e(i)  <   r_mid_below ) Then

                  Do while (index  >   lower_limit .and.                &
     &                      L_continue_e)
                    index = index - 1
                    r_above = r_here
                    r_here  = r_below
                    r_mid_above = r_mid_below
                    r_below =                                           &
     &   coeff_a_e(i) * r_theta_levels (i_out_e(i)  ,j_out_e(i),index-1)&
     & + coeff_b_e(i) * r_theta_levels (i_out_e(i)+1,j_out_e(i),index-1)&
     & + coeff_c_e(i) * r_theta_levels (i_out_e(i),j_out_e(i)+1,index-1)&
     & + coeff_d_e(i) *r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,index-1)
                    r_mid_below=(r_here+r_below)/2.

                    If (r_out_e(i)  >=  r_mid_below .and.               &
     &                  r_out_e(i)  <   r_mid_above ) Then

                      depart_level_e(i) = index
                      w_star_e(i) = (r_theta_levels_e(i) - r_here)      &
     &                               / timestep

                      L_continue_e = .false.
                    End if
                  End Do

                Else
! set w_star values as level index is correct
                  w_star_e(i) = (r_theta_levels_e(i) - r_here)          &
     &                               / timestep
                End If

              End If ! end if on bottom adjustment

            End If ! on which model level

          End Do ! end loop over communication on demand points

! ----------------------------------------------------------------------
! Section 4.4  Max/min values for points requested by other processors
! ----------------------------------------------------------------------

! if monotonocity required then find level just below departure point
! depart_level must be either just below or just above desired level
! so work out which it is.
          If (L_mono_theta ) Then
            Do i = 1, dim_e_out
              r_below = coeff_a_e(i) *                                  &
     & r_theta_levels(i_out_e(i),j_out_e(i),depart_level_e(i))          &
     &                  + coeff_b_e(i) *                                &
     & r_theta_levels (i_out_e(i)+1,j_out_e(i),depart_level_e(i))       &
     &                + coeff_c_e(i) *                                  &
     & r_theta_levels (i_out_e(i),j_out_e(i)+1,depart_level_e(i))       &
     &                  + coeff_d_e(i) *                                &
     & r_theta_levels(i_out_e(i)+1,j_out_e(i)+1,depart_level_e(i))
              If(r_below  <   r_out_e(i) .and.                          &
     &             depart_level_e(i)  /=  1) Then
                level_below_e(i) = depart_level_e(i)
              Else
                level_below_e(i) = depart_level_e(i)  - 1
              End If

! calculate max/min theta values at departure points
                theta_d_max_e(i) = max (                                &
     &       Ext_Data(i_out_e(i),j_out_e(i),level_below_e(i),1),        &
     &       Ext_Data(i_out_e(i)+1,j_out_e(i),level_below_e(i),1),      &
     &       Ext_Data(i_out_e(i),j_out_e(i)+1,level_below_e(i),1),      &
     &       Ext_Data(i_out_e(i)+1,j_out_e(i)+1,level_below_e(i),1),    &
     &       Ext_Data(i_out_e(i),j_out_e(i),level_below_e(i)+1,1),      &
     &       Ext_Data(i_out_e(i)+1,j_out_e(i),level_below_e(i)+1,1),    &
     &       Ext_Data(i_out_e(i),j_out_e(i)+1,level_below_e(i)+1,1),    &
     &       Ext_Data(i_out_e(i)+1,j_out_e(i)+1,level_below_e(i)+1,1) )
                theta_d_min_e(i) = min (                                &
     &       Ext_Data(i_out_e(i),j_out_e(i),level_below_e(i),1),        &
     &       Ext_Data(i_out_e(i)+1,j_out_e(i),level_below_e(i),1),      &
     &       Ext_Data(i_out_e(i),j_out_e(i)+1,level_below_e(i),1),      &
     &       Ext_Data(i_out_e(i)+1,j_out_e(i)+1,level_below_e(i),1),    &
     &       Ext_Data(i_out_e(i),j_out_e(i),level_below_e(i)+1,1),      &
     &       Ext_Data(i_out_e(i)+1,j_out_e(i),level_below_e(i)+1,1),    &
     &       Ext_Data(i_out_e(i),j_out_e(i)+1,level_below_e(i)+1,1),    &
     &       Ext_Data(i_out_e(i)+1,j_out_e(i)+1,level_below_e(i)+1,1) )

            End Do ! end loop over communication on demand points
          End If

! ----------------------------------------------------------------------
! Section 4.5  send data requested back to other processors
! ----------------------------------------------------------------------

! send answers back to processor that asked for them
! 1. w_star
          nsend = 1
          Do i = 0, n_procx-1
            If (n_recvfrom(i)  >   0) then
              len = n_recvfrom(i)
              call gc_rsend(40*(me+1)+ibase+i, len, ibase+i, info,      &
     &                recv_data(1,ime), w_star_e(nsend))
              nsend = nsend + n_recvfrom(i)
            End If
          End Do

          Call gcg_ssync(proc_row_group, info)

          Do i = 0, n_procx-1
            If (n_sendto(i)  >   0) then
              len = n_sendto(i)
              call gc_rrecv(40*(ibase+i+1)+me, len, ibase+i, info,      &
     &                recv_data(1,i), w_star_e)
              Do j = 1, n_sendto(i)
                w_star(i_store(j,i), j_store(j,i),                      &
     &                   k_store(j,i)) = recv_data(j,i)
              End Do
            End If
          End Do

          If (at_extremity(PSouth)) then
            Do j = 0, n_procx-1
              If (sp_send(j)  >   0) then
                Do kk = 1, sp_send(j)
                  k = sp_levels(j,kk)
                  ctmp1(1,kk) = w(1,1,k)
                  ctmp1(2,kk) = w_star(1,1,k)
                End Do
                call gcg_rbcast(101, 2*sp_send(j), ibase + j,           &
     &                   proc_row_group, info, ctmp1)
                do kk = 1, sp_send(j)
                  k = sp_levels(j,kk)
                  If (ctmp1(1,kk)  ==  ctmp1(2,kk)) then
                    Do i = 1, row_length
                      w_star(i,1,k) = w(i,1,k)
                    End Do
                  Else
                    Do i = 1, row_length
                      w_star(i,1,k) = ctmp1(2,kk)
                    End Do
                  End If
                End Do
              End If
            End Do
          End If
          If (at_extremity(PNorth)) then
            Do j = 0, n_procx-1
              If (np_send(j)  >   0) then
                Do kk = 1, np_send(j)
                  k = np_levels(j,kk)
                  ctmp1(1,kk) = w(1,rows,k)
                  ctmp1(2,kk) = w_star(1,rows,k)
                End Do
                call gcg_rbcast(102, 2*np_send(j), ibase + j,           &
     &                   proc_row_group, info, ctmp1)
                do kk = 1, np_send(j)
                  k = np_levels(j,kk)
                  If (ctmp1(1,kk)  ==  ctmp1(2,kk)) then
                    Do i = 1, row_length
                      w_star(i,rows,k) = w(i,rows,k)
                    End Do
                  Else
                    Do i = 1, row_length
                      w_star(i,rows,k) = ctmp1(2,kk)
                    End Do
                  End If
                End Do
              End If
            End Do
          End If

! 2. max/min theta
          If (L_mono_theta) Then
            Do i = 1, dim_e_out
              rsend_arr2(1,i) = theta_d_max_e(i)
              rsend_arr2(2,i) = theta_d_min_e(i)
            End Do

          nsend = 1
          Do i = 0, n_procx-1
            If (n_recvfrom(i)  >   0) then
              len = 2 * n_recvfrom(i)
              call gc_rsend(50*(me+1)+ibase+i, len, ibase+i, info,      &
     &                rrecv_arr2(1,1,ime), rsend_arr2(1,nsend) )
              nsend = nsend + n_recvfrom(i)
            End If
          End Do

          Call gcg_ssync(proc_row_group, info)

          Do i = 0, n_procx-1
            If (n_sendto(i)  >   0) then
              len = n_sendto(i) * 2
              call gc_rrecv(50*(ibase+i+1)+me, len, ibase+i, info,      &
     &                rrecv_arr2(1,1,i), rsend_arr2 )
              Do j = 1, n_sendto(i)
                theta_d_max(i_store(j,i), j_store(j,i),                 &
     &                      k_store(j,i)) = rrecv_arr2(1,j,i)
                theta_d_min(i_store(j,i), j_store(j,i),                 &
     &                      k_store(j,i)) = rrecv_arr2(2,j,i)
              End Do
            End If
          End Do

          If (at_extremity(PSouth)) then
            Do j = 0, n_procx-1
              If (sp_send(j)  >   0) then
                Do kk = 1, sp_send(j)
                  k = sp_levels(j,kk)
                  ctmp1(1,kk) = theta_d_max(1,1,k)
                  ctmp1(2,kk) = theta_d_min(1,1,k)
                End Do
                call gcg_rbcast(101, 2*sp_send(j), ibase + j,           &
     &                   proc_row_group, info, ctmp1)
                do kk = 1, sp_send(j)
                  k = sp_levels(j,kk)
                  Do i = 1, row_length
                    theta_d_max(i,1,k) = ctmp1(1,kk)
                    theta_d_min(i,1,k) = ctmp1(2,kk)
                  End Do
                End Do
              End If
            End Do
          End If
          If (at_extremity(PNorth)) then
            Do j = 0, n_procx-1
              If (np_send(j)  >   0) then
                Do kk = 1, np_send(j)
                  k = np_levels(j,kk)
                  ctmp1(1,kk) = theta_d_max(1,rows,k)
                  ctmp1(2,kk) = theta_d_min(1,rows,k)
                End Do
                call gcg_rbcast(102, 2*np_send(j), ibase + j,           &
     &                   proc_row_group, info, ctmp1)
                do kk = 1, np_send(j)
                  k = np_levels(j,kk)
                  Do i = 1, row_length
                    theta_d_max(i,rows,k) = ctmp1(1,kk)
                    theta_d_min(i,rows,k) = ctmp1(2,kk)
                  End Do
                End Do
              End If
            End Do
          End If
          End If   ! on (L_mono_theta)

        End If  ! on n_procx  >   1 .and. model_domain  ==  mt_Global

! ----------------------------------------------------------------------
! Section 5.    Calculate Eulerian vertical advection term.
! ----------------------------------------------------------------------

        k = 1
          Do j = 1, rows
            Do i = 1, row_length
! CFL=(w-w*) * dt /dz
! CFL limit is 1
              w_minus_wstar(i,j,k) = timestep *                         &
     &                             ( w(i,j,k) - w_star(i,j,k) )         &
     &                    / (r_theta_levels(i,j,k+1) -                  &
     &                       r_theta_levels(i,j,k))
              If (w_minus_wstar(i,j,k)  >   1.0) Then
                w_minus_wstar(i,j,k) = 1.0
              Else If (w_minus_wstar(i,j,k)  <   -1.0) Then
                w_minus_wstar(i,j,k) = -1.0
              End If
              work(i,j,k) = w_minus_wstar(i,j,k) * (theta(i,j,k+1) -    &
     &                                              theta(i,j,k))
            End Do
          End Do
        Do k = 2, model_levels - 1
          Do j = 1, rows
            Do i = 1, row_length
! CFL=(w-w*) * dt / dz
! In this case it is 2 dz so CFL limit is 0.5 not 1
              w_minus_wstar(i,j,k)= timestep *                          &
     &                             ( w(i,j,k) - w_star(i,j,k) )         &
     &                    / (r_theta_levels(i,j,k+1) -                  &
     &                       r_theta_levels(i,j,k-1))
!             If (w_minus_wstar(i,j,k)  >   0.5) Then
!               w_minus_wstar(i,j,k) = 0.5
        w_minus_wstar(i,j,k) = max(-0.5,min(0.5,w_minus_wstar(i,j,k)))
!             Else If (w_minus_wstar(i,j,k)  <   -0.5) Then
!               w_minus_wstar(i,j,k) = -0.5
!             End If
              work(i,j,k) = w_minus_wstar(i,j,k) * (theta(i,j,k+1) -    &
     &                                              theta(i,j,k-1))
            End Do
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 6.    Perform Horizontal Interpolation
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 6.1   Extend data arrays.
! ----------------------------------------------------------------------

        If (model_domain  ==  mt_Global .or.                            &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) Then
          Do k = 1, model_levels - 1
            Do j = 1, rows
              Do i = 1, row_length
                Ext_data (i,j,k,1) = theta(i,j,k) -                     &
     &                               (1.-alpha_2) * work(i,j,k)
              End Do
            End Do
          End Do
          k = model_levels
            Do j = 1, rows
              Do i = 1, row_length
                Ext_data (i,j,k,1) = theta(i,j,k)
              End Do
            End Do

! DEPENDS ON: swap_bounds
            call Swap_Bounds(                                           &
     &                   ext_data(1-halo_i,1-halo_j,1,1),               &
     &                   row_length, rows, model_levels,                &
     &                   halo_i, halo_j, fld_type_p, L_vector)

        Else If (model_domain  ==  mt_LAM ) Then

          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                Ext_data (i,j,k,1) = theta(i,j,k)
              End Do
            End Do
          End Do

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   ext_data(1-halo_i,1-halo_j,1,1),               &
     &                   row_length, rows, model_levels,                &
     &                   halo_i, halo_j, fld_type_p, L_vector)

! Set data on edge processor edge haloes from lateral boundary data

            L_do_halos = .TRUE.
            L_do_boundaries = .FALSE.

! DEPENDS ON: set_lateral_boundaries
            CALL SET_LATERAL_BOUNDARIES(                                &
     &                                 ROW_LENGTH,ROWS, halo_i, halo_j, &
     &                                 MODEL_LEVELS, fld_type_p,        &
     &                                 ext_data(1-halo_i,1-halo_j,1,1), &
     &                                 LENRIM, LBC_SIZE, LBC_START,     &
     &                                 halo_i, halo_j, theta_lbcs,      &
     &                                 RIMWIDTH, RIMWIDTH, RIMWEIGHTS,  &
     &                                 AT_EXTREMITY,                    &
     &                                 L_do_boundaries, L_do_halos)

          Do k = 1, model_levels - 1

            Do j = 1, rows
              Do i = 1, row_length
                Ext_data (i,j,k,2) = work(i,j,k)
              End Do
            End Do
          End Do
            Do j = 1, rows
              Do i = 1, row_length
                Ext_data (i,j,model_levels,2) = 0.
              End Do
            End Do

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   ext_data(1-halo_i,1-halo_j,1,2),               &
     &                   row_length, rows, model_levels,                &
     &                   halo_i, halo_j, fld_type_p, L_vector)

! Set data on edge processor edge haloes to zero

            L_do_halos = .TRUE.
            L_do_boundaries = .FALSE.

! DEPENDS ON: zero_lateral_boundaries
            CALL ZERO_LATERAL_BOUNDARIES(                               &
     &                                 ROW_LENGTH, ROWS,                &
     &                                 halo_i, halo_j,                  &
     &                                 MODEL_LEVELS, fld_type_p,        &
     &                                 ext_data(1-halo_i,1-halo_j,1,2), &
     &                                 RIMWIDTH, AT_EXTREMITY,          &
     &                                 L_do_boundaries, L_do_halos)
        
        End If  !   model_domain

! ----------------------------------------------------------------------
! Section 6.2   Call interpolation code.
! ----------------------------------------------------------------------

! Call high order scheme if required.

        If (L_high_theta ) Then

         If (high_order_scheme_theta  ==  cubicLagrange  .or.           &
     &       high_order_scheme_theta  ==  hCubic_vLin .or.              &
     &       high_order_scheme_theta  ==  hCubic_vQuintic) Then

! DEPENDS ON: cubic_lagrange_niv
            Call Cubic_Lagrange_niv (                                   &
     &                             Ext_Data,                            &
     &                             row_length, rows, model_levels,      &
     &                             row_length, rows, model_levels,      &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, depart_level,          &
     &                             row_length, rows,                    &
     &                             lambda_rm, lambda_rp, phi_rm, phi_rp,&
     &                             recip_lambda_m, recip_lambda_0,      &
     &                             recip_lambda_p, recip_lambda_p2,     &
     &                             recip_phi_m, recip_phi_0,            &
     &                             recip_phi_p, recip_phi_p2,           &
     &                             L_regular,                           &
     &                             model_domain,                        &
     &                             at_extremity, n_procx, n_procy,      &
     &                             g_row_length, g_rows,                &
     &                             proc_col_group, proc_row_group,      &
     &                             datastart,                           &
     &                             off_x, off_y,                        &
     &                             theta_star)

            If (model_domain  ==  mt_LAM) Then

! DEPENDS ON: cubic_lagrange_niv
              Call Cubic_Lagrange_niv (                                 &
     &              Ext_Data(1-halo_i,1-halo_j,-1,2),                   &
     &                             row_length, rows, model_levels,      &
     &                             row_length, rows, model_levels,      &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, depart_level,          &
     &                             row_length, rows,                    &
     &                             lambda_rm, lambda_rp, phi_rm, phi_rp,&
     &                             recip_lambda_m, recip_lambda_0,      &
     &                             recip_lambda_p, recip_lambda_p2,     &
     &                             recip_phi_m, recip_phi_0,            &
     &                             recip_phi_p, recip_phi_p2,           &
     &                             L_regular,                           &
     &                             model_domain,                        &
     &                             at_extremity, n_procx, n_procy,      &
     &                             g_row_length, g_rows,                &
     &                             proc_col_group, proc_row_group,      &
     &                             datastart,                           &
     &                             halo_i, halo_j,                      &
     &                             work2)

            End If

          Else If (high_order_scheme_theta  ==  quinticLagrange) Then

! DEPENDS ON: quintic_lagrange_niv
            Call Quintic_Lagrange_niv (                                 &
     &                             Ext_Data,                            &
     &                             row_length, rows, model_levels,      &
     &                             row_length, rows, model_levels,      &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, depart_level,          &
     &                             off_x, off_y,                        &
     &                             theta_star)

            If (model_domain  ==  mt_LAM) Then

! DEPENDS ON: quintic_lagrange_niv
              Call Quintic_Lagrange_niv (                               &
     &              Ext_Data(1-halo_i,1-halo_j,-1,2),                   &
     &                             row_length, rows, model_levels,      &
     &                             row_length, rows, model_levels,      &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, depart_level,          &
     &                             halo_i, halo_j,                      &
     &                             work2)

            End If

          Else If (high_order_scheme_theta  ==  ECMWF_quasiCubic .or.   &
     &             high_order_scheme_theta  ==  hQuasiCubic_vQuintic)   &
     &       Then

! DEPENDS ON: ecmwf_quasi_cubic_niv
            Call ECMWF_quasi_cubic_niv (                                &
     &                             Ext_Data,                            &
     &                             row_length, rows, model_levels,      &
     &                             row_length, rows, model_levels,      &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, depart_level,          &
     &                             off_x, off_y,                        &
     &                             theta_star)

            If (model_domain  ==  mt_LAM) Then

! DEPENDS ON: ecmwf_quasi_cubic_niv
              Call ECMWF_quasi_cubic_niv (                              &
     &              Ext_Data(1-halo_i,1-halo_j,-1,2),                   &
     &                             row_length, rows, model_levels,      &
     &                             row_length, rows, model_levels,      &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, depart_level,          &
     &                             halo_i, halo_j,                      &
     &                             work2)

            End If

          End If

        End If  !  L_high_theta

        If (L_mono_theta .and. L_high_theta) Then
! Perform montonicity enforcement if high order scheme
! used and monotonicity required.

          Do k = 1, model_levels
            Do j = 1,rows
              Do i = 1, row_length

! Find max and min monotone values for the point concerned.

                max_mono = max (                                        &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k),1),      &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k),1),    &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k),1),    &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k),1) )

               min_mono = min (                                         &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k),depart_level(i,j,k),1),      &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k),depart_level(i,j,k),1),    &
     &  Ext_Data(i_out(i,j,k),j_out(i,j,k)+1,depart_level(i,j,k),1),    &
     &  Ext_Data(i_out(i,j,k)+1,j_out(i,j,k)+1,depart_level(i,j,k),1) )

                If (theta_star(i,j,k)  >   max_mono )                   &
     &              theta_star(i,j,k) = max_mono
                If (theta_star(i,j,k)  <   min_mono )                   &
     &              theta_star(i,j,k) = min_mono

              End Do
            End Do
          End Do

        Else If (L_mono_theta ) Then
! Call monotone scheme if required.

          If (monotone_scheme_theta  ==  triLinear ) Then

! DEPENDS ON: bi_linear_niv
            Call Bi_Linear_niv (                                        &
     &                          Ext_Data,                               &
     &                          row_length, rows, model_levels,         &
     &                          row_length, rows, model_levels,         &
     &                          halo_i, halo_j,                         &
     &                          weight_lambda, weight_phi,              &
     &                          i_out, j_out, depart_level,             &
     &                          off_x, off_y,                           &
     &                          theta_star)

            If (model_domain  ==  mt_LAM) Then

! DEPENDS ON: bi_linear_niv
              Call Bi_Linear_niv (                                      &
     &              Ext_Data(1-halo_i,1-halo_j,-1,2),                   &
     &                             row_length, rows, model_levels,      &
     &                             row_length, rows, model_levels,      &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, depart_level,          &
     &                             halo_i, halo_j,                      &
     &                             work2)

            End If

          Else If (monotone_scheme_theta  ==  mono_quasiCubic) Then

! DEPENDS ON: ecmwf_mono_quasi_cubic_niv
            Call ECMWF_mono_quasi_cubic_niv (                           &
     &                             Ext_Data,                            &
     &                             row_length, rows, model_levels,      &
     &                             row_length, rows, model_levels,      &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, depart_level,          &
     &                             off_x, off_y,                        &
     &                             theta_star)

            If (model_domain  ==  mt_LAM) Then

! DEPENDS ON: ecmwf_mono_quasi_cubic_niv
              Call ECMWF_mono_quasi_cubic_niv (                         &
     &              Ext_Data(1-halo_i,1-halo_j,-1,2),                   &
     &                             row_length, rows, model_levels,      &
     &                             row_length, rows, model_levels,      &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, depart_level,          &
     &                             halo_i, halo_j,                      &
     &                             work2)

            End If  !  model_domain  ==  mt_LAM

          End If  !  monotone_scheme_theta

        End If    !  L_mono_theta .and. L_high_theta

! And now do the extra points

        If (n_procx  >   1 .and. model_domain  ==  mt_Global) then

        If (dim_e_out  >   0) then

        If (L_high_theta ) Then

         If (high_order_scheme_theta  ==  cubicLagrange  .or.           &
     &       high_order_scheme_theta  ==  hCubic_vLin .or.              &
     &       high_order_scheme_theta  ==  hCubic_vQuintic) Then

! DEPENDS ON: cubic_lagrange_niv
            Call Cubic_Lagrange_niv (                                   &
     &                             Ext_Data,                            &
     &                             row_length, rows, model_levels,      &
     &                             dim_e_out, 1, 1,                     &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda_e, weight_phi_e,       &
     &                             i_out_e, j_out_e, depart_level_e,    &
     &                             row_length, rows,                    &
     &                             lambda_rm, lambda_rp, phi_rm, phi_rp,&
     &                             recip_lambda_m, recip_lambda_0,      &
     &                             recip_lambda_p, recip_lambda_p2,     &
     &                             recip_phi_m, recip_phi_0,            &
     &                             recip_phi_p, recip_phi_p2,           &
     &                             L_regular,                           &
     &                             model_domain,                        &
     &                             at_extremity, n_procx, n_procy,      &
     &                             g_row_length, g_rows,                &
     &                             proc_col_group, proc_row_group,      &
     &                             datastart,                           &
     &                             0, 0,                                &
     &                             theta_star_e)


          Else If (high_order_scheme_theta  ==  quinticLagrange) Then

! DEPENDS ON: quintic_lagrange_niv
            Call Quintic_Lagrange_niv (                                 &
     &                             Ext_Data,                            &
     &                             row_length, rows, model_levels,      &
     &                             dim_e_out, 1, 1,                     &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda_e, weight_phi_e,       &
     &                             i_out_e, j_out_e, depart_level_e,    &
     &                             0, 0,                                &
     &                             theta_star_e)


          Else If (high_order_scheme_theta  ==  ECMWF_quasiCubic .or.   &
     &             high_order_scheme_theta  ==  hQuasiCubic_vQuintic)   &
     &       Then

! DEPENDS ON: ecmwf_quasi_cubic_niv
            Call ECMWF_quasi_cubic_niv (                                &
     &                             Ext_Data,                            &
     &                             row_length, rows, model_levels,      &
     &                             dim_e_out, 1, 1,                     &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda_e, weight_phi_e,       &
     &                             i_out_e, j_out_e, depart_level_e,    &
     &                             0, 0,                                &
     &                             theta_star_e)


          End If

        End If

        If (L_mono_theta .and. L_high_theta) Then
! Perform montonicity enforcement if high order scheme
! used and monotonicity required.

          Do i = 1, dim_e_out

! Find max and min monotone values for the point concerned.

            max_mono = max (                                            &
     &         Ext_Data(i_out_e(i),j_out_e(i),depart_level_e(i),1),     &
     &         Ext_Data(i_out_e(i)+1,j_out_e(i),depart_level_e(i),1),   &
     &         Ext_Data(i_out_e(i),j_out_e(i)+1,depart_level_e(i),1),   &
     &         Ext_Data(i_out_e(i)+1,j_out_e(i)+1,depart_level_e(i),1) )

            min_mono = min (                                            &
     &         Ext_Data(i_out_e(i),j_out_e(i),depart_level_e(i),1),     &
     &         Ext_Data(i_out_e(i)+1,j_out_e(i),depart_level_e(i),1),   &
     &         Ext_Data(i_out_e(i),j_out_e(i)+1,depart_level_e(i),1),   &
     &         Ext_Data(i_out_e(i)+1,j_out_e(i)+1,depart_level_e(i),1) )

            If (theta_star_e(i)  >   max_mono )                         &
     &           theta_star_e(i) = max_mono
            If (theta_star_e(i)  <   min_mono )                         &
     &           theta_star_e(i) = min_mono

          End Do

        Else If (L_mono_theta ) Then
! Call monotone scheme if required.

          If (monotone_scheme_theta  ==  triLinear) Then

! DEPENDS ON: bi_linear_niv
            Call Bi_Linear_niv (                                        &
     &                          Ext_Data,                               &
     &                          row_length, rows, model_levels,         &
     &                          dim_e_out, 1, 1,                        &
     &                          halo_i, halo_j,                         &
     &                          weight_lambda_e, weight_phi_e,          &
     &                          i_out_e, j_out_e, depart_level_e,       &
     &                          0, 0,                                   &
     &                          theta_star_e)


          Else If (monotone_scheme_theta  ==  mono_quasiCubic) Then

! DEPENDS ON: ecmwf_mono_quasi_cubic_niv
            Call ECMWF_mono_quasi_cubic_niv (                           &
     &                             Ext_Data,                            &
     &                             row_length, rows, model_levels,      &
     &                             dim_e_out, 1, 1,                     &
     &                             halo_i, halo_j,                      &
     &                             weight_lambda_e, weight_phi_e,       &
     &                             i_out_e, j_out_e, depart_level_e,    &
     &                             0, 0,                                &
     &                             theta_star_e)


          End If

        End If

        End If ! (dim_e_out  >   0)

! NB: Message passing does not cope with model_domain eq 2 as it fails
!     to return work2_e to the originating processor.

        nsend = 0
        Do i = 0, n_procx-1
           If (n_recvfrom(i)  >   0) then
              len = n_recvfrom(i)
              Do j = 1, len
                 send_data(j,i) = theta_star_e(nsend+j)
                 send_data(len+j,i) = r_out_e(nsend+j)
              End Do
              call gc_rsend(60*(me+1)+ibase+i, 2*len, ibase+i, info,    &
     &             recv_data(1,ime), send_data(1,i))
              nsend = nsend + n_recvfrom(i)
           End If
        End Do

        call gcg_ssync(proc_row_group,info)

        Do i = 0, n_procx-1
           If (n_sendto(i)  >   0) then
              len = n_sendto(i)
              call gc_rrecv(60*(ibase+i+1)+me, 2*len, ibase+i, info,    &
     &             recv_data(1,i), send_data(1,ime))
              Do j = 1, len
                 theta_star(i_store(j,i), j_store(j,i),                 &
     &                k_store(j,i)) = recv_data(j,i)
                 r_out(i_store(j,i), j_store(j,i),                      &
     &                k_store(j,i)) = recv_data(len+j,i)
              End Do
           End If
        End Do

        If (at_extremity(PSouth)) then
           do j = 0, n_procx-1
              If (sp_send(j)  >   0) then
                 do kk = 1, sp_send(j)
                    k = sp_levels(j,kk)
                    ctmp1(1,kk) = theta_star(1,1,k)
                    ctmp1(2,kk) = r_out(1,1,k)
                 enddo
                 call gcg_rbcast(103, 2*sp_send(j), ibase+j,            &
     &             proc_row_group, info, ctmp1)
                 do kk = 1, sp_send(j)
                    k = sp_levels(j,kk)
                    Do i = 1, row_length
                       theta_star(i,1,k) = ctmp1(1,kk)
                       r_out(i,1,k) = ctmp1(2,kk)
                    End Do
                 End Do
              End If
           End Do
        End If
        If (at_extremity(PNorth)) then
           do j = 0, n_procx-1
              If (np_send(j)  >   0) then
                 do kk = 1, np_send(j)
                    k = np_levels(j,kk)
                    ctmp1(1,kk) = theta_star(1,rows,k)
                    ctmp1(2,kk) = r_out(1,rows,k)
                 enddo
                 call gcg_rbcast(104, 2*np_send(j), ibase+j,            &
     &                proc_row_group, info, ctmp1)
                 do kk = 1, np_send(j)
                    k = np_levels(j,kk)
                    Do i = 1, row_length
                       theta_star(i,rows,k) = ctmp1(1,kk)
                       r_out(i,rows,k) = ctmp1(2,kk)
                    End Do
                 End Do
              End If
           End Do
        End If

        End If  ! nprocx  >   1

! ----------------------------------------------------------------------
! Section 7.    Add on Eulerian vertical advection term.
! ----------------------------------------------------------------------

        If (model_domain  ==  mt_Global .or.                            &
     &       model_domain  ==  mt_cyclic_LAM .or.                       &
     &       model_domain  ==  mt_bi_cyclic_LAM) Then
          Do j=1,rows
            Do i = 1, row_length
              a_coeff_a(i,j) = alpha_2
            End Do
          End Do

          If (model_domain  ==  mt_cyclic_LAM) Then
            If (at_extremity(PSouth)) Then
              Do i = 1, row_length
                a_coeff_a(i,1) = 1.0
              End Do
            Else If (at_extremity(PNorth)) Then
              Do i = 1, row_length
                a_coeff_a(i,rows) = 1.0
              End Do
            End If
          End If

          If ( ( .NOT. L_new_tdisc ) .OR. CycleNo == 1 ) Then

! calculate first approximation to theta_star
            Do k = 1, model_levels - 1
              Do j = 1, rows
                Do i = 1, row_length
                  work2(i,j,k) = theta_star(i,j,k)                      &
     &                       - a_coeff_a(i,j) * work(i,j,k)
                End Do
              End Do
            End Do

! ----------------------------------------------------------------
!  Theta on top surface is advected only in the horizontal.
! ----------------------------------------------------------------
            k = model_levels
            Do j = 1, rows
              Do i = 1, row_length
                work2(i,j,k) = theta_star(i,j,k)
              End Do
            End Do

! calculate second approximation to theta_star
            k = 1
            Do j = 1, rows
              Do i = 1, row_length
                theta_star(i,j,k) = work2(i,j,k)                        &
     &              - a_coeff_a(i,j) * w_minus_wstar(i,j,k)             &
     &                            * (work2(i,j,k+1) - work2(i,j,k))     &
     &                            + a_coeff_a(i,j) * work(i,j,k)
              End Do
            End Do

            Do k = 2, model_levels - 1
              Do j = 1, rows
                Do i = 1, row_length
                  theta_star(i,j,k) = work2(i,j,k)                      &
     &                     - a_coeff_a(i,j) * w_minus_wstar(i,j,k)      &
     &                            * (work2(i,j,k+1) - work2(i,j,k-1))   &
     &                            + a_coeff_a(i,j) * work(i,j,k)
                End Do
              End Do
            End Do

          Else ! CycleNo >1 .and. L_new_tdisc=T

! calculate second approximation to theta_star
! Compute theta* = theta* - a2*dt*[w^n-w*]d2r(theta^(1))
! where theta^(1) the theta^n+1 estimate from the 1st sweep.

            k = 1
              Do j = 1, rows
                Do i = 1, row_length
                  theta_star(i,j,k) = theta_star(i,j,k)                 &
     &                       - a_coeff_a(i,j) * w_minus_wstar(i,j,k)    &
     &                       * (theta_np1(i,j,k+1) - theta_np1(i,j,k))
                End Do
              End Do

            Do k = 2, model_levels - 1
              Do j = 1, rows
                Do i = 1, row_length
                  theta_star(i,j,k) = theta_star(i,j,k)                 &
     &                      - a_coeff_a(i,j) * w_minus_wstar(i,j,k)     &
     &                      * (theta_np1(i,j,k+1) - theta_np1(i,j,k-1))
                End Do
              End Do
            End Do
          End If  ! .NOT. L_new_tdisc  .OR.  CycleNo == 1

        Else If (model_domain  ==  mt_LAM ) Then

! set up array for alpha_2 to take into account LAM boundaries
          Do j = 1, rows
            Do i = 1, row_length
              a_coeff_a(i,j) = alpha_2
            End Do
          End Do
! Northern Area
          If (at_extremity(PNorth)) Then
            Do j = rows-halo_j+1, rows
              Do i = 1, row_length
               a_coeff_a(i,j) = 1.0
              End Do
            End Do
          End If
! Southern Area
          If (at_extremity(PSouth)) Then
            Do j = 1, halo_j
              Do i = 1, row_length
               a_coeff_a(i,j) = 1.0
              End Do
            End Do
          End If
! Western Area
          If (at_extremity(PWest)) Then
            Do j = 1, rows
              Do i = 1, halo_i
                a_coeff_a(i,j) = 1.0
              End Do
            End Do
          End If
! Eastern Area
          If (at_extremity(PEast)) Then
            Do j = 1, rows
              Do i = row_length-halo_i+1, row_length
                a_coeff_a(i,j) = 1.0
              End Do
            End Do
          End If

          If ( ( .NOT. L_new_tdisc ) .OR. CycleNo == 1 ) Then

! calculate first approximation to theta_star
            Do k = 1, model_levels - 1
              Do j = 1, rows
                Do i = 1, row_length
                  work2(i,j,k) = theta_star(i,j,k)                      &
     &                         - a_coeff_a(i,j) * work(i,j,k)           &
     &                  - (1.0-a_coeff_a(i,j)) * work2(i,j,k)
                End Do
              End Do
            End Do

! ----------------------------------------------------------------
!  Theta on top surface is advected only in the horizontal.
! ----------------------------------------------------------------
            k = model_levels
            Do j = 1, rows
              Do i = 1, row_length
                work2(i,j,k) = theta_star(i,j,k)
              End Do
            End Do

! calculate second approximation to theta_star
            k = 1
              Do j = 1, rows
                Do i = 1, row_length
                  theta_star(i,j,k) = work2(i,j,k)                      &
     &                - a_coeff_a(i,j) * w_minus_wstar(i,j,k)           &
     &                              * (work2(i,j,k+1) - work2(i,j,k))   &
     &                            + a_coeff_a(i,j) * work(i,j,k)
                End Do
              End Do

            Do k = 2, model_levels - 1
              Do j = 1, rows
                Do i = 1, row_length
                  theta_star(i,j,k) = work2(i,j,k)                      &
     &                     - a_coeff_a(i,j) * w_minus_wstar(i,j,k)      &
     &                              * (work2(i,j,k+1) - work2(i,j,k-1)) &
     &                              + a_coeff_a(i,j) * work(i,j,k)
                End Do
              End Do
            End Do

          Else  ! L_new_tdisc .and. CycleNo > 1

! calculate work2=theta_dl - (1-a2)*dt*[(w^n-w^*)*d2r(theta)]_dl

            Do k = 1, model_levels - 1
              Do j = 1, rows
                Do i = 1, row_length
                  work2(i,j,k) = theta_star(i,j,k)                      &
     &                  - (1.0-a_coeff_a(i,j)) * work2(i,j,k)
                End Do
              End Do
            End Do

! calculate second approximation to theta_star:
! theta* = theta_dl - (1-a2)*dt*[(w^n-w^*)*d2r(theta)]_dl
!                   - a2*dt*(w^(1)_w*)*d2r(theta^(1))

            k = 1
              Do j = 1, rows
                Do i = 1, row_length
                  theta_star(i,j,k) = work2(i,j,k)                      &
     &                       - a_coeff_a(i,j) * w_minus_wstar(i,j,k)    &
     &                       * (theta_np1(i,j,k+1) - theta_np1(i,j,k))
                End Do
              End Do

            Do k = 2, model_levels - 1
              Do j = 1, rows
                Do i = 1, row_length
                  theta_star(i,j,k) = work2(i,j,k)                      &
     &                      - a_coeff_a(i,j) * w_minus_wstar(i,j,k)     &
     &                      * (theta_np1(i,j,k+1) - theta_np1(i,j,k-1))
                End Do
              End Do
            End Do

          End If !  .NOT. L_new_tdisc .OR. CycleNo == 1

        End If ! on model_domain


! end conditional on zero error code
      End If

! End of routine.
      return
      END SUBROUTINE Calc_non_int_sl_theta

