
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Interpolation_qcon.
!

      Subroutine Interpolation_qcon(                                    &
     &                          Data_in1, Data_in2, Data_in3,           &
     &                          eta_in,                                 &
     &                          r_in, delta_r_in, delta_r_np1,          &
     &                          r_in_w, type,                           &
     &                          number_of_inputs,                       &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          first_flat_level_in,                    &
     &                          dim_i_in, dim_j_in, dim_k_in,           &
     &                          dim_j_in_w,                             &
     &                          dim_i_out, dim_j_out, dim_k_out,        &
     &                          glambda_p, phi_p,                       &
     &                          lambda_rm, lambda_rp, phi_rm, phi_rp,   &
     &                          recip_lambda_m, recip_lambda_0,         &
     &                          recip_lambda_p, recip_lambda_p2,        &
     &                          recip_phi_m, recip_phi_0,               &
     &                          recip_phi_p, recip_phi_p2,              &
     &                          i_out_in, j_out_in,                     &
     &                          weight_lambda, weight_phi,              &
     &                          high_order_scheme, monotone_scheme,     &
     &                          cos_latitude, L_regular,                &
     &                          L_vector, model_domain, L_high, L_mono, &
     &                          L_conserv,                              &
     &                          r_out, lambda_out, phi_out,             &
     &                          me, n_proc, n_procx, n_procy,           &
     &                          halo_i, halo_j, g_row_length,           &
     &                          g_rows, row_length, rows,               &
     &                          datastart, at_extremity, g_i_pe,        &
     &                          proc_row_group, proc_col_group,         &
     &                          pole_handling_in,                       &
     &                          halo_data_out_i, halo_data_out_j,       &
     &                          L_sl_halo_reprod, off_x, off_y,         &
     &                          Data_out1, Data_out2, Data_out3,        &
     &                          Error_Code)

! Purpose:
!          Performs interpolation of a field or fields defined on one
!          grid to another grid. The current grids that can be handled
!          can be found in the documentation as can the current options
!          for interpolation schemes. The desired interpolation can be
!          monotone and conservative if desired. Input data can be on a
!          sphere or be a rectangular box. Requested output points must
!          lie inside the region defined for the input Data.
!          The number of fields interpolated is either 1,2 or 3 and is
!          controlled by the switch number_of_inputs. If only one field
!          is to be interpolated the fields Data_in2/3 and Data_out2/3
!          should be set to dummy arguments as they are not used.
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
! Version   Date      Comment
! ----     --------   -------
!LL5.2     02/11/00   New routine. INTERP2A + improved conservation
!                                  algorithm                 A. Malcolm
!  5.3    11/07/01   add interpolation options 6 & 7       Andy Malcolm
!  5.3    25/10/01   Change from magic numbers            M.Diamantakis
!  5.3    14/11/01   Dynamic allocation of memory introduced. S. Cusack
!  5.4     28/08/02    Bug Fix (Bi-cyclic LAM)           Carol Roadnight
!  6.0  18/08/03  NEC SX-6 optimisations - Allocatable replaced by
!                 Dynamic arrays and logical test loop by list.
!                 R Barnes & J-C Rioual.
!  6.0  18/08/03  Make GCOM tags fit within MPI limits. A. Malcolm
!  6.1  17/08/04  NEC optimised compile directives S.S.Wilson
!  6.2  25/12/05  Variable resolution changes            Yongming Tang
!  6.4  26/01/07  Send first_flat_level to eta_vert_weights_e  A.Malcolm
!  6.4  27/11/06  Fix to integer arithmetic bug          A.Malcolm

!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  dim_i_in                                                        &
                    ! Dimension of Data_in in i direction.
     &, dim_j_in                                                        &
                    ! Dimension of Data_in in j direction.
     &, dim_j_in_w                                                      &
                    ! Dimension of Data_in in j direction.
     &, dim_k_in                                                        &
                    ! Dimension of Data_in in k direction.
     &, dim_i_out                                                       &
                    ! Dimension of Data_out in i direction.
     &, dim_j_out                                                       &
                    ! Dimension of Data_out in j direction.
     &, dim_k_out                                                       &
                    ! Dimension of Data_out in k direction.
     &, me                                                              &
                    ! My processor number
     &, n_proc                                                          &
                    ! Total number of processors
     &, n_procx                                                         &
                    ! Number of processors in longitude
     &, n_procy                                                         &
                    ! Number of processors in latitude
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j                                                          &
                    ! Size of halo in j direction.
     &, off_x                                                           &
     &, off_y                                                           &
     &, halo_data_out_i                                                 &
                        ! size of data out halo in i direction
     &, halo_data_out_j                                                 &
                        ! size of data out halo in j direction
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, pole_handling_in                                                &
                          ! How to treat the poles:
                          !   0 - no calculations at the poles
                          !   1 - poles in one point
                          !   2 - poles in all points (default)
     &, g_row_length                                                    &
                             ! global number of points on a row
     &, g_rows                                                          &
                             ! global rows
     &, row_length                                                      &
                     ! row_length on this pe for dynamic arrays
     &, rows                                                            &
                     ! rows on this pe for dynamic arrays
     &, datastart(3)                                                    &
                     ! First gridpoints held by this processor.
     &, g_i_pe(1-halo_i:g_row_length+halo_i) ! processor on my procr-row
                             ! holding a given value in i direction

      Logical                                                           &
     &  L_sl_halo_reprod ! if true then sl code bit repoducible with
                         ! any sensible halo size

      Integer                                                           &
     &  type                                                            &
                    ! Defines via an integer code the nature of the
                    ! interpolation in terms of which grid the input
                    ! Data is on. The codes are given in
                    ! terms of a primary variable that would be held
                    ! at that point and are u=1, v=2, w=3.
     &, high_order_scheme                                               &
                           ! a code saying which high order scheme to
                           ! use.
     &, monotone_scheme                                                 &
                        ! a code saying which monotone scheme to use.
     &, number_of_inputs                                                &
                         !number of fields to interpolate.
     &, interp_vertical_search_tol                                      &
                                   !number of levels either side of
                                   ! default level to search.
     &, check_bottom_levels                                             &
                            ! used in interpolation code, and is
                            ! the number of levels to check to see
                            ! if the departure point lies inside the
                            ! orography.
     &, first_flat_level_in

      Logical                                                           &
     &  L_vector                                                        &
                       ! True, if data is a horizontal vector component,
                       ! False, then data is a scalar.
     &, L_high                                                          &
                       ! True, if high order interpolation required.
     &, L_mono                                                          &
                       ! True, if interpolation required to be monotone.
     &, L_conserv                                                       &
                       ! True, if interpolation to be monotone and
                       !       conservative.
     &, L_regular      ! False if variable resolution

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

!  VarRes horizontal co-ordinate spacing etc.
       Real                                                             &
     &  glambda_p(1-halo_i : g_row_length + halo_i)                     &
     &, phi_p( 1-halo_i : row_length + halo_i                           &
     &,        1-halo_j : rows + halo_j )                               &
     &, lambda_rm (1-halo_i : row_length + halo_i)                      &
     &, lambda_rp (1-halo_i : row_length + halo_i)                      &
     &, phi_rm    ( 1-halo_i : row_length + halo_i                      &
     &,             1-halo_j : rows + halo_j )                          &
     &, phi_rp    ( 1-halo_i : row_length + halo_i                      &
     &,             1-halo_j : rows + halo_j )                          &
     &, recip_lambda_m(1-halo_i : row_length + halo_i)                  &
     &, recip_lambda_0(1-halo_i : row_length + halo_i)                  &
     &, recip_lambda_p(1-halo_i : row_length + halo_i)                  &
     &, recip_lambda_p2(1-halo_i : row_length + halo_i)                 &
     &, recip_phi_m( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows + halo_j )                         &
     &, recip_phi_0( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows + halo_j )                         &
     &, recip_phi_p( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows + halo_j )                         &
     &, recip_phi_p2( 1-halo_i : row_length + halo_i                    &
     &,               1-halo_j : rows + halo_j )

      Real                                                              &
     &  eta_in(dim_k_in) ! eta coordinate levels.

      Real                                                              &
     &  Data_in1 (1-halo_i:dim_i_in+halo_i,                             &
     &            1-halo_j:dim_j_in+halo_j, dim_k_in)                   &
                                                      ! data to be
                                                      ! interpolated
     &, Data_in2 (1-halo_i:dim_i_in+halo_i,                             &
     &            1-halo_j:dim_j_in+halo_j, dim_k_in)                   &
                                                      ! optional second
                                                      ! field of data to
                                                      ! be interpolated
     &, Data_in3 (1-halo_i:dim_i_in+halo_i,                             &
     &            1-halo_j:dim_j_in+halo_j, dim_k_in)                   &
                                                      ! optional third
                                                      ! field of data to
                                                      ! be interpolated
     &, r_in (1-halo_i:dim_i_in+halo_i,                                 &
     &        1-halo_j:dim_j_in+halo_j, dim_k_in)                       &
                                                      ! Vertical
                                                      ! co-ordinate
                                                      ! of input data.
     &, r_in_w(1-halo_i:dim_i_in+halo_i,                                &
     &        1-halo_j:dim_j_in_w+halo_j, dim_k_in)                     &
                                                      ! Vertical
                                                      ! co-ordinate
                                                      ! of input data
                                                      ! on w grid
     &, delta_r_in (dim_i_in, dim_j_in, dim_k_in)                       &
                                                      ! Vertical
                                                      ! layer thickness
                                                      ! of input data.
     &, delta_r_np1 (dim_i_in, dim_j_in, dim_k_in)                      &
                                                      ! Vertical
                                                      ! layer thickness
                                                      ! of input data.
     &, cos_latitude (1-off_x:dim_i_in+off_x,                           &
     &                1-off_y:dim_j_in+off_y)
      Real                                                              &
     &  lambda_out (dim_i_out, dim_j_out, dim_k_out)                    &
                                                      ! Lambda
                                                      ! co-ordinate of
                                                      ! output data on
                                                      ! input.
     &, phi_out (dim_i_out, dim_j_out, dim_k_out)                       &
                                                      ! Phi Co-ordinate
                                                      ! of output data
                                                      ! on input.
     &, r_out (dim_i_out, dim_j_out, dim_k_out)       ! Vertical
                                                      ! co-ordinate
                                                      ! of output data.

      Real                                                              &
     &  weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)

      Integer                                                           &
     &  i_out_in (dim_i_out, dim_j_out, dim_k_out)                      &
     &, j_out_in (dim_i_out, dim_j_out, dim_k_out)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
           ! data interpolated to desired locations.
     &  Data_out1(1-halo_data_out_i:dim_i_out+halo_data_out_i,          &
     &            1-halo_data_out_j:dim_j_out+halo_data_out_j,          &
     &            dim_k_out)                                            &
     &, Data_out2(1-halo_data_out_i:dim_i_out+halo_data_out_i,          &
     &            1-halo_data_out_j:dim_j_out+halo_data_out_j,          &
     &            dim_k_out)                                            &
     &, Data_out3(1-halo_data_out_i:dim_i_out+halo_data_out_i,          &
     &            1-halo_data_out_j:dim_j_out+halo_data_out_j,          &
     &            dim_k_out)

      Integer                                                           &
     &  Error_Code     ! Non-zero on exit if error detected.

! Local Variables.

! scalars

      Logical                                                           &
     &  Conserv_fail

      Integer                                                           &
     &  i, j, k, n                                                      &
                   ! Loop indices
     &, lambda_start 
                     ! pointer for start of lambda_p/lambda_u on this pe

! arrays

      Integer                                                           &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
     &, j_out (dim_i_out, dim_j_out, dim_k_out)                         &
     &, k_out (dim_i_out, dim_j_out, dim_k_out)                         &
     &, i_out_w (dim_i_out, dim_j_out, dim_k_out)                       &
     &, j_out_w (dim_i_out, dim_j_out, dim_k_out)


      Real                                                              &
     &  Ext_Data(1-halo_i:dim_i_in+halo_i+1,1-halo_j:dim_j_in+halo_j,   &
     &           -1:dim_k_in+2,number_of_inputs)                        &
     &, Data_out_high (dim_i_out, dim_j_out, dim_k_out,                 &
     &                 number_of_inputs)                                &
     &, Data_out_mono (dim_i_out, dim_j_out, dim_k_out,                 &
     &                 number_of_inputs)                                &
     &, weight_lambda_w (dim_i_out, dim_j_out, dim_k_out)               &
     &, weight_phi_w (dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeff_z(dim_i_out, dim_j_out, dim_k_out, -2:3)                  &
     &, coeff_z_lin(dim_i_out, dim_j_out, dim_k_out, 0:1)


! Varibles applied in the "compute-on-demand" strategy

      Integer                                                           &
     &  ime, ibase, irecv, my_imin, my_imax, dim_e_out, h_factor        &
     &, nsend, nrecv, info, len, itmp, j0, j1, kk, sender               &
     &, pole_handling, my_iminp, my_imaxp

!
      Integer                                                           &
     &  n_sendto(0:n_procx-1), n_recvfrom(0:n_procx-1)                  &
     &, i_store(dim_k_out*g_row_length,0:n_procx-1)                     &
     &, j_store(dim_k_out*g_row_length,0:n_procx-1)                     &
     &, k_store(dim_k_out*g_row_length,0:n_procx-1)                     &
     &, i_out_e(dim_k_out*g_row_length)                                 &
     &, j_out_e(dim_k_out*g_row_length)                                 &
     &, k_out_e(dim_k_out*g_row_length)                                 &
     &, k_level_e(dim_k_out*g_row_length)                               &
                                          ! vertical levl of arrival pnt
     &, i_out_w_e(dim_k_out*g_row_length)                               &
     &, j_out_w_e(dim_k_out*g_row_length)                               &
     &, isend_arr(5,dim_k_out*g_row_length,0:n_procx-1)                 &
     &, irecv_arr(5,dim_k_out*g_row_length,0:n_procx-1)                 &
     &, sp_send(0:n_procx-1), sp_levels(0:n_procx-1,dim_k_out)          &
     &, np_send(0:n_procx-1), np_levels(0:n_procx-1,dim_k_out)
      Real                                                              &
     &  rsend_arr(5,dim_k_out*g_row_length,0:n_procx-1)                 &
     &, rrecv_arr(5,dim_k_out*g_row_length,0:n_procx-1)                 &
     &, weight_lambda_e(dim_k_out*g_row_length)                         &
     &, weight_phi_e(dim_k_out*g_row_length)                            &
     &, weight_lambda_w_e(dim_k_out*g_row_length)                       &
     &, weight_phi_w_e(dim_k_out*g_row_length)                          &
     &, coeff_z_e(dim_k_out*g_row_length, -2:3)                         &
     &, coeff_z_lin_e(dim_k_out*g_row_length, 0:1)                      &
     &, r_out_e(dim_k_out*g_row_length)                                 &
     &, Data_out_high_e(dim_k_out*g_row_length*number_of_inputs)        &
     &, Data_out_mono_e(dim_k_out*g_row_length*number_of_inputs)        &
     &, send_data(2*number_of_inputs*dim_k_out*g_row_length,0:n_procx-1)  &
     &, recv_data(2*number_of_inputs*dim_k_out*g_row_length,0:n_procx-1)  &
     &, bcast_data(4*number_of_inputs*dim_k_out)

      Integer :: ErrorStatus

! Variables used for index vectors

      Integer                                                           &
     &  ic,icnt,inx(dim_i_out)

      integer d_imin,d_imax
      Integer                                                           &
     &  lon_dem(3), gon_dem(3)
      Common / ondem / lon_dem, gon_dem

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

! External Routines:
! subroutines
      External Extend_data_linear, Extend_data_cubic                    &
     &,        Extend_data_quintic
      External Cubic_Lagrange, Quintic_Lagrange
      External ECMWF_quasi_cubic, ECMWF_mono_quasi_cubic
      External Tri_Linear, h_cubic_v_linear, Mono_Enforce
      External Mono_Conserv_qcon, eta_vert_weights, eta_vert_weights_e

! Functions: None

! ----------------------------------------------------------------------
!  Section 1.   Error trap Un-supported options.
! ----------------------------------------------------------------------

      Error_Code = 0
      pole_handling = pole_handling_in
      If (pole_handling  /=  0 .and. pole_handling  /=  1 .and.         &
     &     pole_handling  /=  2) pole_handling = 2

! Execute rest of routine only if error code is still zero.

      If (Error_Code  ==  0 ) Then

! ----------------------------------------------------------------------
! Section 2.   Call appropriate routine to extend data.
!              Minimum amount of extending done to cope with highest
!              order interpolation scheme requested. This can leave
!              unset values in Ext_Data and Ext_r_in, but
!              these values will not be used.
!
!              Parallel version: No extension in r required.
! ----------------------------------------------------------------------

        If ( (high_order_scheme  ==  quinticLagrange .or.               &
     &        high_order_scheme  ==  hQuasiCubic_vQuintic .or.          &
     &        high_order_scheme  ==  hCubic_vQuintic) .and. L_High )    &
     &        Then

! DEPENDS ON: extend_data_quintic
          Call Extend_data_quintic (Data_in1, Data_in2, Data_in3,       &
     &                              dim_i_in, dim_j_in, dim_k_in,       &
     &                              number_of_inputs, halo_i, halo_j,   &
     &                              Ext_data)

        Else If ( ( ( high_order_scheme  ==  cubicLagrange              &
     &           .or. high_order_scheme  ==  ECMWF_quasiCubic           &
     &           .or. high_order_scheme  ==  ECMWF_mono_quasiCubic      &
     &           .or. high_order_scheme  ==  hCubic_vLin )              &
     &          .and. L_High )                                          &
     &           .or. ( monotone_scheme  ==  mono_quasiCubic            &
     &          .and. L_Mono ) ) Then

! DEPENDS ON: extend_data_cubic
          Call Extend_data_cubic  (Data_in1, Data_in2, Data_in3,        &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             number_of_inputs, halo_i, halo_j,    &
     &                             Ext_data)

        Else

! DEPENDS ON: extend_data_linear
          Call Extend_data_linear (Data_in1, Data_in2, Data_in3,        &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             number_of_inputs, halo_i, halo_j,    &
     &                             Ext_data)

        End If

! ----------------------------------------------------------------------
! Section 3.   For each output point find i,j,k so that the point on the
!              output grid lies between i and i+1, j and j+1, k and k+1
! ----------------------------------------------------------------------

! i_out and j_out must be copied because they change in this subroutine
      Do k = 1, dim_k_out
        Do j = 1, dim_j_out
          Do i = 1, dim_i_out
            i_out(i,j,k) = i_out_in(i,j,k)
            j_out(i,j,k) = j_out_in(i,j,k)
            i_out_w(i,j,k) = i_out(i,j,k)
            weight_lambda_w(i,j,k) =  weight_lambda(i,j,k)
            j_out_w(i,j,k) = j_out(i,j,k)
            weight_phi_w(i,j,k) = weight_phi(i,j,k)
          End Do
        End Do
      End Do

       If ( n_procx > 1 ) then
       If ( model_domain == mt_Global ) then
! Send the points outside my region to the appropriate processor for
! interpolation. Only performed if the domain is decomposed in the
! i direction.

! The first and last point I can interpolate in, based on available
! data on this processor (minus/plus one to avoid use of ge/le)

          h_factor = 2
          If (high_order_scheme  ==  quinticLagrange .and. L_High) Then
            h_factor = 3
          End If
          my_imin = datastart(1) - halo_i + h_factor - 1
          my_imax = datastart(1) + dim_i_out - 1 + halo_i - h_factor    &
     &              + 1

! values for use in polar row to ensure pole is only calculated on one
! processor, (minus/plus one to avoid use of ge/le)
          my_iminp = datastart(1) - 1
          my_imaxp = datastart(1)+dim_i_out
          If(at_extremity(PWest) ) my_iminp = my_imin
          If(at_extremity(PEast) ) my_imaxp = my_imax

! The base processor on this row, and my address relative to that
! processor

          ibase = (me/n_procx) * n_procx
          ime = me - ibase

          Do i = 0, n_procx-1
            n_sendto(i) = 0
          End Do

          j0 = 1
          j1 = dim_j_out

! If the pole values not are going to be used, we we can save a
! large amount of communication by evaluating local dummys at the poles

          If (pole_handling  ==  0) then

            If (at_extremity(PSouth)) then
              j0 = 2
              Do k = 1, dim_k_out
                Do i = 1, dim_i_out
                  i_out(i,1,k) = i
                  i_out_w(i,1,k) = i
                End Do
              End Do
            End If

            If (at_extremity(PNorth)) then
              j1 = dim_j_out-1
              Do k = 1, dim_k_out
                Do i = 1, dim_i_out
                  i_out(i,dim_j_out,k) = i
                  i_out_w(i,dim_j_out,k) = i
                End Do
              End Do
            End If

          End If

! If all values along a polar row is evaluated in one point, we can
! save a large amount of communication by evaluating only one point
! correctly and then broadcast this value.

          If (pole_handling  ==  1) then

            Do i = 0, n_procx-1
              sp_send(i) = 0
            End Do
            If (at_extremity(PSouth)) then
              j0 = 2
              Do k = 1, dim_k_out
                If (i_out(1,1,k)  >   my_iminp .and.                    &
     &              i_out(1,1,k)  <   my_imaxp) then
                  i_out(1,1,k) = i_out(1,1,k) - datastart(1) + 1
                  i_out_w(1,1,k) = i_out_w(1,1,k) - datastart(1) + 1
                  sp_send(ime) = sp_send(ime) + 1
                  sp_levels(ime,sp_send(ime)) = k
                  Do i = 2, dim_i_out
                    i_out(i,1,k) = i ! i_out(i,1,k) - datastart(1) + 1
                    i_out_w(i,1,k) = i
                  End Do
                Else
                  sender = g_i_pe(i_out(1,1,k))
                  sp_send(sender) = sp_send(sender) + 1
                  sp_levels(sender,sp_send(sender)) = k
                  Do i = 1, dim_i_out
                    i_out(i,1,k) = i
                    i_out_w(i,1,k) = i
                  End Do
                End If
              End Do
            End If

            Do i = 0, n_procx-1
              np_send(i) = 0
            End Do
            If (at_extremity(PNorth)) then
              j1 = dim_j_out-1
              Do k = 1, dim_k_out
                If (i_out(1,dim_j_out,k)  >   my_iminp .and.            &
     &              i_out(1,dim_j_out,k)  <   my_imaxp)  then

                  np_send(ime) = np_send(ime) + 1
                  np_levels(ime,np_send(ime)) = k
                  i_out(1,dim_j_out,k) =                                &
     &                 i_out(1,dim_j_out,k) - datastart(1) + 1
                  i_out_w(1,dim_j_out,k) =                              &
     &                 i_out_w(1,dim_j_out,k) - datastart(1) + 1
                  Do i = 2, dim_i_out
                    i_out(i,dim_j_out,k) = i
                    i_out_w(i,dim_j_out,k) = i
                  End Do
                Else
                  sender = g_i_pe(i_out(1,dim_j_out,k))
                  np_send(sender) = np_send(sender) + 1
                  np_levels(sender,np_send(sender)) = k
                  Do i = 1, dim_i_out
                    i_out(i,dim_j_out,k) = i
                    i_out_w(i,dim_j_out,k) = i
                  End Do
                End If
              End Do
            End If

          End If ! on pole_handling

          If( L_sl_halo_reprod) Then

! On the global boundaries, use i_out < 1 or i_out > g_row_length
! if that makes local computation possible. Not required when
! L_sl_halo_reprod is false is other logic ensures this is done.

! This code unsafe if applied at poles, where it isn't required.

          If (at_extremity(PWest)) then
            Do k = 1, dim_k_out
              Do j = j0, j1
                Do i = 1, halo_i
                  If (i_out(i,j,k)  >   g_row_length-halo_i+h_factor)   &
     &              Then
                     i_out(i,j,k) = i_out(i,j,k) - g_row_length
                     i_out_w(i,j,k) = i_out_w(i,j,k) - g_row_length
                  End If
                End Do
              End Do
            End Do
          End If
          If (at_extremity(PEast)) then
            Do k = 1, dim_k_out
              Do j = j0, j1
                Do i = dim_i_out-halo_i+1, dim_i_out
                  If (i_out(i,j,k)  <   halo_i-h_factor+1)              &
     &              Then
                      i_out(i,j,k) = i_out(i,j,k) + g_row_length
                      i_out_w(i,j,k) = i_out_w(i,j,k) + g_row_length
                  End If
                End Do
              End Do
            End Do
          End If

          End If ! on L_sl_halo_reprod

! And now decide where a point should be evaluated

          d_imin=my_imin-datastart(1) + 1
          d_imax=my_imax-datastart(1) + 1

!               If (i_out(i,j,k)  >   my_imin .and.
!    &              i_out(i,j,k)  <   my_imax) then

          Do k = 1, dim_k_out
            Do j = j0, j1
              Do i = 1, dim_i_out
                  i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
                  i_out_w(i,j,k) = i_out_w(i,j,k) - datastart(1) + 1
              End Do
            End Do
          End Do

          itmp=0
          Do k = 1, dim_k_out
            Do j = j0, j1
              Do i = 1, dim_i_out
                If (i_out(i,j,k)  <=  d_imin .or.                       &
     &              i_out(i,j,k)  >=  d_imax) then
                  itmp=itmp+1
                End If
              End Do
            End Do
          End Do
          Do k = 1, dim_k_out
            Do j = j0, j1
              icnt=0
              Do i = 1, dim_i_out
                If (i_out(i,j,k)  <=  d_imin .or.                       &
     &              i_out(i,j,k)  >=  d_imax) then
                    icnt = icnt+1
                    inx(icnt) = i
                End If
              End Do

! Process locally, so find the local destination
              Do ic=1, icnt
                 i = inx(ic)
                  i_out(i,j,k) = i_out(i,j,k) + datastart(1) - 1
                  i_out_w(i,j,k) = i_out_w(i,j,k) + datastart(1) - 1
! Send to a remote processor, given by the array g_i_pe
!     CODE TO STOP BIT NON-REPRODUCIBILITY
      if(i_out(i,j,k) >  g_row_length+halo_i-h_factor)then
        i_out(i,j,k)=i_out(i,j,k)-g_row_length
        i_out_w(i,j,k)=i_out_w(i,j,k)-g_row_length
      endif
      if(i_out(i,j,k) <  1-halo_i+h_factor)then
        i_out(i,j,k)=i_out(i,j,k)+g_row_length
        i_out_w(i,j,k)=i_out_w(i,j,k)+g_row_length
      endif
!     END CODE TO STOP BIT NON-REPRODUCIBILITY
! sza
                  irecv = g_i_pe(i_out(i,j,k))
                  n_sendto(irecv) = n_sendto(irecv) + 1
                  itmp = n_sendto(irecv)
                  isend_arr(1,itmp,irecv) = i_out(i,j,k)
                  isend_arr(2,itmp,irecv) = j_out(i,j,k)
                  isend_arr(3,itmp,irecv) = i_out_w(i,j,k)
                  isend_arr(4,itmp,irecv) = j_out_w(i,j,k)
                  isend_arr(5,itmp,irecv) = k
                  rsend_arr(1,itmp,irecv) = weight_lambda(i,j,k)
                  rsend_arr(2,itmp,irecv) = weight_phi(i,j,k)
                  rsend_arr(3,itmp,irecv) = r_out(i,j,k)
                  rsend_arr(4,itmp,irecv) = weight_lambda_w(i,j,k)
                  rsend_arr(5,itmp,irecv) = weight_phi_w(i,j,k)
                  i_store(itmp,irecv) = i
                  j_store(itmp,irecv) = j
                  k_store(itmp,irecv) = k
                  i_out(i,j,k) = i
                  i_out_w(i,j,k) = i
!                End If
              End Do
            End Do
          End Do

! Send the points to be evaluated outside my region

          nsend = 0
          Do i = 0,n_procx-1
            call gc_isend(10*(me+1)+ibase+i, 1, ibase+i, info,          &
     &           n_recvfrom(ime), n_sendto(i))
            nsend = nsend + n_sendto(i)
          End Do
          lon_dem(1) = lon_dem(1) + nsend

          Call gcg_ssync(proc_row_group, info)

          nrecv = 0
          Do i = 0, n_procx-1
            call gc_irecv(10*(ibase+i+1)+me, 1, ibase+i, info,          &
     &           n_recvfrom(i), n_sendto(ime))
            nrecv = nrecv + n_recvfrom(i)
          End Do

!          Call gcg_ssync(proc_row_group, info)

          Do i = 0,n_procx-1
            If (n_sendto(i)  >   0) then
              len = 5*n_sendto(i)
              call gc_rsend(20*(me+1)+ibase+i, len, ibase+i, info,      &
     &             rrecv_arr(1,1,ime), rsend_arr(1,1,i))
            End If
          End Do

          Call gcg_ssync(proc_row_group, info)

          Do i = 0,n_procx-1
            If (n_recvfrom(i)  >   0) then
              len = 5*n_recvfrom(i)
              call gc_rrecv(20*(ibase+i+1)+me, len, ibase+i, info,      &
     &             rrecv_arr(1,1,i), rsend_arr(1,1,ime))
            End If
          End Do

!          Call gcg_ssync(proc_row_group, info)

          Do i = 0,n_procx-1
            If (n_sendto(i)  >   0) then
              len = 5*n_sendto(i)
              call gc_isend(30*(me+1)+ibase+i, len, ibase+i, info,      &
     &             irecv_arr(1,1,ime), isend_arr(1,1,i))
            End If
          End Do

          Call gcg_ssync(proc_row_group, info)

          Do i = 0,n_procx-1
            If (n_recvfrom(i)  >   0) then
              len = 5*n_recvfrom(i)
              call gc_irecv(30*(ibase+i+1)+me, len, ibase+i, info,      &
     &             irecv_arr(1,1,i), isend_arr(1,1,ime))
            End If
          End Do


          dim_e_out = 0
          Do i = 0, n_procx-1
            If (n_recvfrom(i)  >   0) then
              Do j = 1, n_recvfrom(i)
                dim_e_out = dim_e_out + 1
                i_out_e(dim_e_out) = irecv_arr(1,j,i) - datastart(1)+1
                j_out_e(dim_e_out) = irecv_arr(2,j,i)
                i_out_w_e(dim_e_out) = irecv_arr(3,j,i) - datastart(1)+1
                j_out_w_e(dim_e_out) = irecv_arr(4,j,i)
                k_level_e(dim_e_out) = irecv_arr(5,j,i)
                weight_lambda_e(dim_e_out) = rrecv_arr(1,j,i)
                weight_phi_e(dim_e_out) = rrecv_arr(2,j,i)
                r_out_e(dim_e_out) = rrecv_arr(3,j,i)
                weight_lambda_w_e(dim_e_out) = rrecv_arr(4,j,i)
                weight_phi_w_e(dim_e_out) = rrecv_arr(5,j,i)
              End Do
            End If
          End Do

          If (dim_e_out  >   dim_k_out*g_row_length) Then
            ErrorStatus = 10
! DEPENDS ON: Ereport
            Call Ereport("Interpolation_qcon", ErrorStatus,             &
     &           "over-writing due to dim_e_out size" )
          End If

          Call gc_ssync(n_proc, info)

          If (dim_e_out  >   0) Then
! calculate level just below departure point in vertical and
! vertical interpolation weights.
! DEPENDS ON: eta_vert_weights_e
            Call eta_vert_weights_e(                                    &
     &                          eta_in, r_in_w,                         &
     &                          check_bottom_levels,                    &
     &                          interp_vertical_search_tol,             &
     &                          first_flat_level_in,                    &
     &                          dim_i_in, dim_j_in_w, dim_k_in,         &
     &                          dim_k_out, dim_e_out,                   &
     &                          high_order_scheme, monotone_scheme,     &
     &                          model_domain, L_high, L_mono,           &
     &                          L_conserv,                              &
     &                          halo_i, halo_j,                         &
     &                          r_out_e, k_level_e,                     &
     &                          coeff_z_e, coeff_z_lin_e,               &
     &                          k_out_e, i_out_w_e, j_out_w_e,          &
     &                          weight_lambda_w_e, weight_phi_w_e )

          End If

       Else   !  model_domain is type of LAM

          k = 1
          j = 1
          Do i = 1, dim_i_out * dim_j_out * dim_k_out
            i_out(i,j,k) = i_out(i,j,k) - datastart(1) + 1
            i_out_w(i,j,k) = i_out_w(i,j,k) - datastart(1) + 1
          End Do ! i = 1, dim_i_out * dim_j_out * dim_k_out

       EndIf ! model_domain == mt_Global

       EndIf ! n_procx > 1


!     CODE TO STOP BIT NON-REPRODUCIBILITY
      if(model_domain == mt_global .and. n_procx == 1)then
        h_factor = 2
        If (high_order_scheme  ==  quinticLagrange .and. L_High) Then
          h_factor = 3
        End If
        my_imin = datastart(1) - halo_i + h_factor - 1
        my_imax = datastart(1) + dim_i_out - 1 + halo_i - h_factor + 1
        Do k = 1, dim_k_out
          Do j = 1,dim_j_out
            Do i = 1, dim_i_out
              if(i_out(i,j,k) >= my_imax)then
                i_out(i,j,k)=i_out(i,j,k) - g_row_length
                i_out_w(i,j,k)=i_out_w(i,j,k) - g_row_length
              endif
              if(i_out(i,j,k) <= my_imin )then
                i_out(i,j,k)=i_out(i,j,k)+g_row_length
                i_out_w(i,j,k)=i_out_w(i,j,k)+g_row_length
              endif
            End Do
          End Do
        End Do
      endif   !model_domain == 1 .and. n_procx == 1
!     END CODE TO STOP BIT NON-REPRODUCIBILITY

! calculate level just below departure point in vertical and
! vertical interpolation weights.
! DEPENDS ON: eta_vert_weights
        Call eta_vert_weights(                                          &
     &                        eta_in, r_in_w,                           &
     &                        check_bottom_levels,                      &
     &                        interp_vertical_search_tol,               &
     &                        first_flat_level_in,                      &
     &                        dim_i_in, dim_j_in_w, dim_k_in,           &
     &                        dim_i_out, dim_j_out, dim_k_out,          &
     &                        high_order_scheme, monotone_scheme,       &
     &                        model_domain, L_high, L_mono,             &
     &                        L_conserv,                                &
     &                        halo_i, halo_j,                           &
     &                        r_out,                                    &
     &                        coeff_z, coeff_z_lin,                     &
     &                        k_out, i_out_w, j_out_w,                  &
     &                        weight_lambda_w, weight_phi_w )

! ----------------------------------------------------------------------
! Section 4.   Perform required Interpolations.
! ----------------------------------------------------------------------


! Call high order scheme if required.

        If (L_high ) Then

          If (high_order_scheme  ==  cubicLagrange) Then

! DEPENDS ON: cubic_lagrange
            Call Cubic_Lagrange (                                       &
     &                           Ext_Data,                              &
     &                           dim_i_in, dim_j_in, dim_k_in,          &
     &                           dim_i_out, dim_j_out, dim_k_out,       &
     &                           halo_i, halo_j, number_of_inputs,      &
     &                           weight_lambda, weight_phi,             &
     &                           i_out, j_out, k_out,                   &
     &                           row_length, rows,                      &
     &                           lambda_rm, lambda_rp, phi_rm, phi_rp,  &
     &                           recip_lambda_m, recip_lambda_0,        &
     &                           recip_lambda_p, recip_lambda_p2,       &
     &                           recip_phi_m, recip_phi_0,              &
     &                           recip_phi_p, recip_phi_p2,             &
     &                           coeff_z, L_regular,                    &
     &                           model_domain,                          &
     &                           at_extremity, n_procx, n_procy,        &
     &                           g_row_length, g_rows,                  &
     &                           proc_col_group, proc_row_group,        &
     &                           datastart,                             &
     &                           Data_out_high)

          Else If (high_order_scheme  ==  quinticLagrange) Then

! DEPENDS ON: quintic_lagrange
            Call Quintic_Lagrange (                                     &
     &                               Ext_Data,                          &
     &                               dim_i_in, dim_j_in, dim_k_in,      &
     &                               dim_i_out, dim_j_out, dim_k_out,   &
     &                               halo_i, halo_j, number_of_inputs,  &
     &                               weight_lambda, weight_phi,         &
     &                               i_out, j_out, k_out,               &
     &                               coeff_z,                           &
     &                               Data_out_high )

          Else If (high_order_scheme  ==  ECMWF_quasiCubic) Then

! DEPENDS ON: ecmwf_quasi_cubic
            Call ECMWF_quasi_cubic (                                    &
     &                                Ext_Data,                         &
     &                                dim_i_in, dim_j_in, dim_k_in,     &
     &                                dim_i_out, dim_j_out, dim_k_out,  &
     &                                halo_i, halo_j, number_of_inputs, &
     &                                weight_lambda, weight_phi,        &
     &                                i_out, j_out, k_out,              &
     &                                coeff_z,                          &
     &                                Data_out_high )

          Else If (high_order_scheme  ==  ECMWF_mono_quasiCubic) Then

! DEPENDS ON: ecmwf_mono_quasi_cubic
            Call ECMWF_mono_quasi_cubic (                               &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_i_out, dim_j_out, dim_k_out,     &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, k_out,                 &
     &                             coeff_z,                             &
     &                             Data_out_high)

          Else If (high_order_scheme  ==  hCubic_vLin) Then

! DEPENDS ON: h_cubic_v_linear
            Call h_cubic_v_linear (                                     &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_i_out, dim_j_out, dim_k_out,     &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, k_out,                 &
     &                             coeff_z_lin,                         &
     &                             Data_out_high)

          Else If (high_order_scheme  ==  hQuasiCubic_vQuintic) Then

! DEPENDS ON: h_quasi_cubic_v_quintic
            Call h_quasi_cubic_v_quintic (                              &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_i_out, dim_j_out, dim_k_out,     &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, k_out,                 &
     &                             coeff_z,                             &
     &                             Data_out_high)

          Else If (high_order_scheme  ==  hCubic_vQuintic) Then

! DEPENDS ON: h_cubic_v_quintic
            Call h_cubic_v_quintic (                                    &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_i_out, dim_j_out, dim_k_out,     &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda, weight_phi,           &
     &                             i_out, j_out, k_out,                 &
     &                             coeff_z,                             &
     &                             Data_out_high)



          End If

        End If

        If (((L_High .and. L_mono) .or. L_Conserv)                      &
     &       .and. high_order_scheme  /=  ECMWF_mono_quasiCubic) Then

! DEPENDS ON: mono_enforce
          Call Mono_Enforce(                                            &
     &                      Ext_Data, number_of_inputs,                 &
     &                      dim_i_in, dim_j_in, dim_k_in,               &
     &                      dim_i_out, dim_j_out, dim_k_out,            &
     &                      halo_i, halo_j,                             &
     &                      i_out, j_out, k_out,                        &
     &                      Data_out_high)

        End If

! Call monotone scheme if required.

        If ( (L_mono .and. .not. L_high ) .or. L_Conserv ) Then

          If (monotone_scheme  ==  triLinear) Then

! DEPENDS ON: tri_linear
              Call Tri_Linear (                                         &
     &                         Ext_Data,                                &
     &                         dim_i_in, dim_j_in, dim_k_in,            &
     &                         dim_i_out, dim_j_out, dim_k_out,         &
     &                         halo_i, halo_j, number_of_inputs,        &
     &                         weight_lambda, weight_phi,               &
     &                         i_out, j_out, k_out,                     &
     &                         coeff_z_lin,                             &
     &                         Data_out_mono )

          Else If (monotone_scheme  ==  mono_quasiCubic) Then

! DEPENDS ON: ecmwf_mono_quasi_cubic
            Call ECMWF_mono_quasi_cubic (                               &
     &                       Ext_Data,                                  &
     &                       dim_i_in, dim_j_in, dim_k_in,              &
     &                       dim_i_out, dim_j_out, dim_k_out,           &
     &                       halo_i, halo_j, number_of_inputs,          &
     &                       weight_lambda, weight_phi,                 &
     &                       i_out, j_out, k_out,                       &
     &                       coeff_z,                                   &
     &                       Data_out_mono)

          End If

        End If

! Repeat the interpolation procedure for the "compute-on-demand" points

        If (n_procx  >   1 .and. model_domain  ==  mt_Global) then

        Call gc_ssync(n_proc, info)

! Call high order scheme if required.

        If (dim_e_out  >   0) then

        If (L_high ) Then

          If (high_order_scheme  ==  cubicLagrange) Then

! DEPENDS ON: cubic_lagrange
            Call Cubic_Lagrange (                                       &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_e_out, 1, 1,                     &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda_e,weight_phi_e,        &
     &                             i_out_e, j_out_e, k_out_e,           &
     &                             row_length, rows,                    &
     &                             lambda_rm, lambda_rp, phi_rm, phi_rp,&
     &                             recip_lambda_m, recip_lambda_0,      &
     &                             recip_lambda_p, recip_lambda_p2,     &
     &                             recip_phi_m, recip_phi_0,            &
     &                             recip_phi_p, recip_phi_p2,           &
     &                             coeff_z_e, L_regular,                &
     &                             model_domain,                        &
     &                             at_extremity, n_procx, n_procy,      &
     &                             g_row_length, g_rows,                &
     &                             proc_col_group, proc_row_group,      &
     &                             datastart,                           &
     &                             Data_out_high_e )

          Else If (high_order_scheme  ==  quinticLagrange) Then

! DEPENDS ON: quintic_lagrange
            Call Quintic_Lagrange (                                     &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_e_out, 1, 1,                     &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda_e,weight_phi_e,        &
     &                             i_out_e, j_out_e, k_out_e,           &
     &                             coeff_z_e,                           &
     &                             Data_out_high_e)

          Else If (high_order_scheme  ==  ECMWF_quasiCubic) Then

! DEPENDS ON: ecmwf_quasi_cubic
            Call ECMWF_quasi_cubic (                                    &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_e_out, 1, 1,                     &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda_e,weight_phi_e,        &
     &                             i_out_e, j_out_e, k_out_e,           &
     &                             coeff_z_e,                           &
     &                             Data_out_high_e)

          Else If (high_order_scheme  ==  ECMWF_mono_quasiCubic) Then

! DEPENDS ON: ecmwf_mono_quasi_cubic
            Call ECMWF_mono_quasi_cubic (                               &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_e_out, 1, 1,                     &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda_e,weight_phi_e,        &
     &                             i_out_e, j_out_e, k_out_e,           &
     &                             coeff_z_e,                           &
     &                             Data_out_high_e)

          Else If (high_order_scheme  ==  hCubic_vLin) Then


! DEPENDS ON: h_cubic_v_linear
              Call h_cubic_v_linear (                                   &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_e_out, 1, 1,                     &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda_e,weight_phi_e,        &
     &                             i_out_e, j_out_e, k_out_e,           &
     &                             coeff_z_lin_e,                       &
     &                             Data_out_high_e)

          Else If (high_order_scheme  ==  hQuasiCubic_vQuintic) Then

! DEPENDS ON: h_quasi_cubic_v_quintic
            Call h_quasi_cubic_v_quintic (                              &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_e_out, 1,1,                      &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda_e, weight_phi_e,       &
     &                             i_out_e, j_out_e, k_out_e,           &
     &                             coeff_z_e,                           &
     &                             Data_out_high_e)

          Else If (high_order_scheme  ==  hCubic_vQuintic) Then

! DEPENDS ON: h_cubic_v_quintic
            Call h_cubic_v_quintic (                                    &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_e_out, 1,1,                      &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda_e, weight_phi_e,       &
     &                             i_out_e, j_out_e, k_out_e,           &
     &                             coeff_z_e,                           &
     &                             Data_out_high_e)



          End If

        End If

        If (((L_High .and. L_mono) .or. L_Conserv)                      &
     &       .and. high_order_scheme  /=  ECMWF_mono_quasiCubic) Then

! DEPENDS ON: mono_enforce
          Call Mono_Enforce(                                            &
     &                      Ext_Data, number_of_inputs,                 &
     &                      dim_i_in, dim_j_in, dim_k_in,               &
     &                      dim_e_out, 1, 1,                            &
     &                      halo_i, halo_j,                             &
     &                      i_out_e, j_out_e, k_out_e,                  &
     &                      Data_out_high_e)

        End If

! Call monotone scheme if required.

        If ( (L_mono .and. .not. L_high ) .or. L_Conserv ) Then

          If (monotone_scheme  ==  triLinear) Then

! DEPENDS ON: tri_linear
              Call Tri_Linear (                                         &
     &                         Ext_Data,                                &
     &                         dim_i_in, dim_j_in, dim_k_in,            &
     &                         dim_e_out, 1, 1,                         &
     &                         halo_i, halo_j, number_of_inputs,        &
     &                         weight_lambda_e,weight_phi_e,            &
     &                         i_out_e, j_out_e, k_out_e,               &
     &                         coeff_z_lin_e,                           &
     &                         Data_out_mono_e)

          Else If (monotone_scheme  ==  mono_quasiCubic) Then

! DEPENDS ON: ecmwf_mono_quasi_cubic
            Call ECMWF_mono_quasi_cubic (                               &
     &                             Ext_Data,                            &
     &                             dim_i_in, dim_j_in, dim_k_in,        &
     &                             dim_e_out, 1, 1,                     &
     &                             halo_i, halo_j, number_of_inputs,    &
     &                             weight_lambda_e,weight_phi_e,        &
     &                             i_out_e, j_out_e, k_out_e,           &
     &                             coeff_z_e,                           &
     &                             Data_out_mono_e)

          End If

        End If

        End If ! (dim_e_out  >   0)

! Time to return the points I have computed on demand.

        nsend = 0

        Call gc_ssync(n_proc, info)

        Do i = 0, n_procx-1
          If (n_recvfrom(i)  >   0) then
            If (L_high .and. L_mono .and. L_conserv) then
              len = 2*n_recvfrom(i)*number_of_inputs
              Do k = 1, number_of_inputs
                Do j = 1, n_recvfrom(i)
                  send_data(2*(k-1)*n_recvfrom(i)+2*(j-1)+1,i)          &
     &                 = Data_out_mono_e(nsend+j+(k-1)*dim_e_out)
                  send_data(2*(k-1)*n_recvfrom(i)+2*(j-1)+2,i)          &
     &                 = Data_out_high_e(nsend+j+(k-1)*dim_e_out)
                End Do
              End Do
              call gc_rsend(40*(me+1)+ibase+i, len, ibase+i, info,      &
     &                      recv_data(1,ime), send_data(1,i))
              nsend = nsend + n_recvfrom(i)
            Else If (L_high) then
              len = n_recvfrom(i)*number_of_inputs
              Do k = 1, number_of_inputs
                Do j = 1, n_recvfrom(i)
                  send_data((k-1)*n_recvfrom(i)+j,i)                    &
     &                 = Data_out_high_e(nsend+j+(k-1)*dim_e_out)
                End Do
              End Do
              call gc_rsend(40*(me+1)+ibase+i, len, ibase+i, info,      &
     &             recv_data(1,ime), send_data(1,i))
              nsend = nsend + n_recvfrom(i)
            Else
              len = n_recvfrom(i)*number_of_inputs
              Do k = 1, number_of_inputs
                Do j = 1, n_recvfrom(i)
                  send_data((k-1)*n_recvfrom(i)+j,i)                    &
     &                 = Data_out_mono_e(nsend+j+(k-1)*dim_e_out)
                End Do
              End Do
              call gc_rsend(40*(me+1)+ibase+i, len, ibase+i, info,      &
     &             recv_data(1,ime), send_data(1,i))
              nsend = nsend + n_recvfrom(i)
            End If
          End If
        End Do

        Call gcg_ssync(proc_row_group, info)

        Do i = 0, n_procx-1
          If (n_sendto(i)  >   0) then
            If (L_high .and. L_mono .and. L_conserv) then
              len = 2*n_sendto(i)*number_of_inputs
            Else
              len = n_sendto(i)*number_of_inputs
            End If

            call gc_rrecv(40*(ibase+i+1)+me, len, ibase+i, info,        &
     &             recv_data(1,i), send_data(1,ime))

            If (L_high .and. L_mono .and. L_conserv) then
              Do n = 1, number_of_inputs
                Do j = 1, n_sendto(i)
                  Data_out_mono(i_store(j,i), j_store(j,i),             &
     &                   k_store(j,i), n)                               &
     &                   = recv_data(2*(n-1)*n_sendto(i)+2*(j-1)+1,i)
                  Data_out_high(i_store(j,i), j_store(j,i),             &
     &                   k_store(j,i), n)                               &
     &                   = recv_data(2*(n-1)*n_sendto(i)+2*(j-1)+2,i)
                End Do
              End Do
            Else If ( L_high ) then
              Do n = 1, number_of_inputs
                Do j = 1, n_sendto(i)
                  Data_out_high(i_store(j,i), j_store(j,i),             &
     &                 k_store(j,i), n)                                 &
     &                 = recv_data((n-1)*n_sendto(i)+j,i)
                End Do
              End Do
            Else
              Do n = 1, number_of_inputs
                Do j = 1, n_sendto(i)
                  Data_out_mono(i_store(j,i), j_store(j,i),             &
     &                 k_store(j,i), n)                                 &
     &                 = recv_data((n-1)*n_sendto(i)+j,i)
                End Do
              End Do
            End If
          End If
        End Do

!        Call gcg_ssync(proc_row_group, info)

! Now distribute the pole values.

        If (pole_handling  ==  1) then

          If (at_extremity(PSouth)) then
            If (L_high .and. L_mono .and. L_conserv) then
              Do j = 0, n_procx-1
                If (sp_send(j)  >   0) then
                  Do n = 1,number_of_inputs
                    Do kk = 1, sp_send(j)
                      k = sp_levels(j,kk)
                      bcast_data(2*(n-1)*sp_send(j)+2*(kk-1)+1) =       &
     &                     Data_out_mono(1,1,k,n)
                      bcast_data(2*(n-1)*sp_send(j)+2*(kk-1)+2) =       &
     &                     Data_out_high(1,1,k,n)
                    End Do
                  End Do
                  len = 2 * number_of_inputs * sp_send(j)
                  call gcg_rbcast(201, len, ibase+j,                    &
     &                 proc_row_group, info, bcast_data)
                  Do n = 1,number_of_inputs
                    Do kk = 1, sp_send(j)
                      k = sp_levels(j,kk)
                      Do i = 1, dim_i_out
                        Data_out_mono(i,1,k,n) =                        &
     &                      bcast_data(2*(n-1)*sp_send(j)+2*(kk-1)+1)
                        Data_out_high(i,1,k,n) =                        &
     &                      bcast_data(2*(n-1)*sp_send(j)+2*(kk-1)+2)
                      End Do
                    End Do
                  End Do
                End If
              End Do
            Else If (L_high) then
              Do j = 0, n_procx-1
                If (sp_send(j)  >   0) then
                  Do n = 1,number_of_inputs
                    Do kk = 1, sp_send(j)
                      k = sp_levels(j,kk)
                      bcast_data((n-1)*sp_send(j)+kk) =                 &
     &                     Data_out_high(1,1,k,n)
                    End Do
                  End Do
                  len = number_of_inputs * sp_send(j)
                  call gcg_rbcast(201, len, ibase+j,                    &
     &                 proc_row_group, info, bcast_data)
                  Do n = 1,number_of_inputs
                    Do kk = 1, sp_send(j)
                      k = sp_levels(j,kk)
                      Do i = 1, dim_i_out
                        Data_out_high(i,1,k,n) =                        &
     &                       bcast_data((n-1)*sp_send(j)+kk)
                      End Do
                    End Do
                  End Do
                End If
              End Do
            Else
              Do j = 0, n_procx-1
                If (sp_send(j)  >   0) then
                  Do n = 1,number_of_inputs
                    Do kk = 1, sp_send(j)
                      k = sp_levels(j,kk)
                      bcast_data((n-1)*sp_send(j)+kk) =                 &
     &                     Data_out_mono(1,1,k,n)
                    End Do
                  End Do
                  len = number_of_inputs * sp_send(j)
                  call gcg_rbcast(201, len, ibase+j,                    &
     &                 proc_row_group, info, bcast_data)
                  Do n = 1,number_of_inputs
                    Do kk = 1, sp_send(j)
                      k = sp_levels(j,kk)
                      Do i = 1, dim_i_out
                        Data_out_mono(i,1,k,n) =                        &
     &                       bcast_data((n-1)*sp_send(j)+kk)
                      End Do
                    End Do
                  End Do
                End If
              End Do
            End If
          End If

          If (at_extremity(PNorth)) then
            If (L_high .and. L_mono .and. L_conserv) then
              Do j = 0, n_procx-1
                If (np_send(j)  >   0) then
                  Do n = 1,number_of_inputs
                    Do kk = 1, np_send(j)
                      k = np_levels(j,kk)
                      bcast_data(2*(n-1)*np_send(j)+2*(kk-1)+1) =       &
     &                     Data_out_mono(1,dim_j_out,k,n)
                      bcast_data(2*(n-1)*np_send(j)+2*(kk-1)+2) =       &
     &                     Data_out_high(1,dim_j_out,k,n)
                    End Do
                  End Do
                  len = 2 * number_of_inputs * np_send(j)
                  call gcg_rbcast(202, len, ibase+j,                    &
     &                 proc_row_group, info, bcast_data)
                  Do n = 1,number_of_inputs
                    Do kk = 1, np_send(j)
                      k = np_levels(j,kk)
                      Do i = 1, dim_i_out
                        Data_out_mono(i,dim_j_out,k,n) =                &
     &                       bcast_data(2*(n-1)*np_send(j)+2*(kk-1)+1)
                        Data_out_high(i,dim_j_out,k,n) =                &
     &                       bcast_data(2*(n-1)*np_send(j)+2*(kk-1)+2)
                      End Do
                    End Do
                  End Do
                End If
              End Do
            Else If (L_high) then
              Do j = 0, n_procx-1
                If (np_send(j)  >   0) then
                  Do n = 1,number_of_inputs
                    Do kk = 1, np_send(j)
                      k = np_levels(j,kk)
                      bcast_data((n-1)*np_send(j)+kk) =                 &
     &                     Data_out_high(1,dim_j_out,k,n)
                    End Do
                  End Do
                  len = number_of_inputs * np_send(j)
                  call gcg_rbcast(202, len, ibase+j,                    &
     &                 proc_row_group, info, bcast_data)
                  Do n = 1,number_of_inputs
                    Do kk = 1, np_send(j)
                      k = np_levels(j,kk)
                      Do i = 1, dim_i_out
                        Data_out_high(i,dim_j_out,k,n) =                &
     &                       bcast_data((n-1)*np_send(j)+kk)
                      End Do
                    End Do
                  End Do
                End If
              End Do
            Else
              Do j = 0, n_procx-1
                If (np_send(j)  >   0) then
                  Do n = 1,number_of_inputs
                    Do kk = 1, np_send(j)
                      k = np_levels(j,kk)
                      bcast_data((n-1)*np_send(j)+kk) =                 &
     &                     Data_out_mono(1,dim_j_out,k,n)
                    End Do
                  End Do
                  len = number_of_inputs * np_send(j)
                  call gcg_rbcast(202, len, ibase+j,                    &
     &                 proc_row_group, info, bcast_data)
                  Do n = 1,number_of_inputs
                    Do kk = 1, np_send(j)
                      k = np_levels(j,kk)
                      Do i = 1, dim_i_out
                        Data_out_mono(i,dim_j_out,k,n) =                &
     &                       bcast_data((n-1)*np_send(j)+kk)
                      End Do
                    End Do
                  End Do
                End If
              End Do
            End If

          End If

        End If ! (pole_handling  ==  1)


        End If ! (n_procx  >   1 .and. model_domain  ==  mt_Global)

! ----------------------------------------------------------------------
! Section 5.   Perform conservation if desired.
!              Put interpolated field into Data_out.
! ----------------------------------------------------------------------

        If ( L_high .and. L_mono .and. L_conserv ) Then

          lambda_start = datastart(1) - halo_i

! DEPENDS ON: mono_conserv_qcon
          Call Mono_Conserv_qcon (                                      &
     &                       Data_out_high(1,1,1,1),                    &
     &                       Data_out_mono(1,1,1,1),                    &
     &                       Ext_Data(1-halo_i,1-halo_j,-1,1),          &
     &                       delta_r_in, delta_r_np1,                   &
     &                       dim_i_in, dim_j_in, dim_k_in,              &
     &                       dim_i_out, dim_j_out, dim_k_out,           &
     &                       row_length, rows,                          &
     &                       glambda_p(lambda_start), phi_p, L_regular, &
     &                       cos_latitude, off_x, off_y,                &
     &                       n_proc, halo_i, halo_j, model_domain,      &
     &                       me, proc_row_group, proc_col_group,        &
     &                       halo_data_out_i, halo_data_out_j,          &
     &                       Data_out1, conserv_fail)

            If ( conserv_fail .and. me  ==  0 ) Then
              print*,' non-conservation for field 1.'
            End If

          If (number_of_inputs  >=  2) Then

! DEPENDS ON: mono_conserv_qcon
            Call Mono_Conserv_qcon (                                    &
     &                       Data_out_high(1,1,1,2),                    &
     &                       Data_out_mono(1,1,1,2),                    &
     &                       Ext_Data(1-halo_i,1-halo_j,-1,2),          &
     &                       delta_r_in, delta_r_np1,                   &
     &                       dim_i_in, dim_j_in, dim_k_in,              &
     &                       dim_i_out, dim_j_out, dim_k_out,           &
     &                       row_length, rows,                          &
     &                       glambda_p(lambda_start), phi_p, L_regular, &
     &                       cos_latitude, off_x, off_y,                &
     &                       n_proc, halo_i, halo_j, model_domain,      &
     &                       me, proc_row_group, proc_col_group,        &
     &                       halo_data_out_i, halo_data_out_j,          &
     &                       Data_out2, conserv_fail)

            If ( conserv_fail .and. me  ==  0 ) Then
              print*,' non-conservation for field 2.'
            End If

          End If

          If (number_of_inputs  >=  3) Then

! DEPENDS ON: mono_conserv_qcon
            Call Mono_Conserv_qcon (                                    &
     &                       Data_out_high(1,1,1,3),                    &
     &                       Data_out_mono(1,1,1,3),                    &
     &                       Ext_Data(1-halo_i,1-halo_j,-1,3),          &
     &                       delta_r_in, delta_r_np1,                   &
     &                       dim_i_in, dim_j_in, dim_k_in,              &
     &                       dim_i_out, dim_j_out, dim_k_out,           &
     &                       row_length, rows,                          &
     &                       glambda_p(lambda_start), phi_p, L_regular, &
     &                       cos_latitude, off_x, off_y,                &
     &                       n_proc, halo_i, halo_j, model_domain,      &
     &                       me, proc_row_group, proc_col_group,        &
     &                       halo_data_out_i, halo_data_out_j,          &
     &                       Data_out3, conserv_fail)

            If ( conserv_fail .and. me  ==  0 ) Then
              print*,' non-conservation for field 3.'
            End If

          End If

        Else If ( L_high ) Then

          Do k = 1, dim_k_out
            Do j = 1, dim_j_out
              Do i = 1, dim_i_out
                Data_out1 (i,j,k) = Data_out_high (i,j,k,1)
              End Do
            End Do
          End Do

          If (number_of_inputs  >=  2) Then
            Do k = 1, dim_k_out
              Do j = 1, dim_j_out
                Do i = 1, dim_i_out
                  Data_out2 (i,j,k) = Data_out_high (i,j,k,2)
                End Do
              End Do
            End Do
          End If

          If (number_of_inputs  >=  3) Then
            Do k = 1, dim_k_out
              Do j = 1, dim_j_out
                Do i = 1, dim_i_out
                  Data_out3 (i,j,k) = Data_out_high (i,j,k,3)
                End Do
              End Do
            End Do
          End If

        Else

          Do k = 1, dim_k_out
            Do j = 1, dim_j_out
              Do i = 1, dim_i_out
                Data_out1 (i,j,k) = Data_out_mono (i,j,k,1)
              End Do
            End Do
          End Do

          If (number_of_inputs  >=  2) Then
            Do k = 1, dim_k_out
              Do j = 1, dim_j_out
                Do i = 1, dim_i_out
                  Data_out2 (i,j,k) = Data_out_mono (i,j,k,2)
                End Do
              End Do
            End Do
          End If

          If (number_of_inputs  >=  3) Then
            Do k = 1, dim_k_out
              Do j = 1, dim_j_out
                Do i = 1, dim_i_out
                  Data_out3 (i,j,k) = Data_out_mono (i,j,k,3)
                End Do
              End Do
            End Do
          End If

        End If


! End if statement on Error Code.
      End If

! End of routine.
      return
      END SUBROUTINE Interpolation_qcon

