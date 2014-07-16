#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Cubic_Lagrange
!

      Subroutine Cubic_Lagrange(                                        &
     &                          Ext_Data,                               &
     &                          dim_i_in, dim_j_in, dim_k_in,           &
     &                          dim_i_out, dim_j_out, dim_k_out,        &
     &                          halo_i, halo_j, number_of_inputs,       &
     &                          weight_lambda, weight_phi,              &
     &                          i_out, j_out, k_out,                    &
     &                          row_length_in, rows_in,                 &
     &                          lambda_rm, lambda_rp, phi_rm, phi_rp,   &
     &                          recip_lambda_m, recip_lambda_0,         &
     &                          recip_lambda_p, recip_lambda_p2,        &
     &                          recip_phi_m, recip_phi_0,               &
     &                          recip_phi_p, recip_phi_p2,              &
     &                          coeff_z, L_regular,                     &
     &                          model_domain,                           &
     &                          at_extremity, n_procx, n_procy,         &
     &                          global_row_length, global_rows,         &
     &                          proc_col_group, proc_row_group,         &
     &                          datastart,                              &
     &                          Data_out)

! Purpose:
!          Performs cubic Lagrange interpolation of the input field to a
!          set of points defined by i_out, j_out, k_out, and
!          weight_lambda,weight_phi,r_out.
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
! Date     Version     Comment
! ----     -------     -------
!  5.1  14/02/00   Upgrade from v2p7 to v2p8a           Andy Malcolm
!  6.0  18/08/03  NEC SX-6 optimisation - loop modifications.
!                 R Barnes & J-C Rioual.
!  6.1  20/07/04  Change code to be faster, especially for tracer
!                 super_arrays.                        Andy Malcolm
!  6.2  25/12/05  Variable resolution changes            Yongming Tang
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  L_regular                                                       &
                         ! false if variable resolution
     &, at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Integer                                                           &
     &  dim_i_in                                                        &
                    ! Dimension of Data_in in i direction.
     &, dim_j_in                                                        &
                    ! Dimension of Data_in in j direction.
     &, dim_k_in                                                        &
                    ! Dimension of Data_in in k direction.
     &, dim_i_out                                                       &
                    ! Dimension of Data_out in i direction.
     &, dim_j_out                                                       &
                    ! Dimension of Data_out in j direction.
     &, dim_k_out                                                       &
                    ! Dimension of Data_out in k direction.
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j                                                          &
                    ! Size of halo in j direction.
     &, number_of_inputs                                                &
                            ! number of fields to interpolate
     &, row_length_in                                                   &
                            ! for lambda dynamic arrays
     &, rows_in                                                         &
                            ! for phi dynamic arrays
     &, model_domain                                                    &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, proc_row_group                                                  &
     &, proc_col_group                                                  &
     &, global_row_length                                               &
                            ! number of points in domain row
     &, global_rows                                                     &
                            ! number of rows in domain
     &, datastart(3)

      Real                                                              &
     &  Ext_Data (1-halo_i:dim_i_in+halo_i+1,                           &
     &            1-halo_j:dim_j_in+halo_j,                             &
                                                   !data to be
     &           -1:dim_k_in+2, number_of_inputs)                       &
                                                   ! interpolated
     &, weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
                                                        ! a number
                                                      ! between 0 & 1
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)  ! a number between
                                                      ! 0 & 1

        Real                                                            &
     &  lambda_rm (1-halo_i : row_length_in+halo_i)                     &
     &, lambda_rp (1-halo_i : row_length_in+halo_i)                     &
     &, phi_rm    ( 1-halo_i : row_length_in + halo_i                   &
     &,             1-halo_j : rows_in + halo_j )                       &
     &, phi_rp    ( 1-halo_i : row_length_in + halo_i                   &
     &,             1-halo_j : rows_in + halo_j )                       &
     &, recip_lambda_m (1-halo_i : row_length_in+halo_i)                &
     &, recip_lambda_0 (1-halo_i : row_length_in+halo_i)                &
     &, recip_lambda_p (1-halo_i : row_length_in+halo_i)                &
     &, recip_lambda_p2(1-halo_i : row_length_in+halo_i)                &
     &, recip_phi_m ( 1-halo_i : row_length_in + halo_i                 &
     &,               1-halo_j : rows_in + halo_j )                     &
     &, recip_phi_0 ( 1-halo_i : row_length_in + halo_i                 &
     &,               1-halo_j : rows_in + halo_j )                     &
     &, recip_phi_p ( 1-halo_i : row_length_in + halo_i                 &
     &,               1-halo_j : rows_in + halo_j )                     &
     &, recip_phi_p2( 1-halo_i : row_length_in + halo_i                 &
     &,               1-halo_j : rows_in+ halo_j )

      Integer                                                           &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! point such that
     &, j_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! the desired
                                                      ! output point
     &, k_out (dim_i_out, dim_j_out, dim_k_out)       ! lies between it
                                                      ! and it+1

! Arguments with Intent IN/OUT. ie: Output variables.
      Real                                                              &
     &  coeff_z (dim_i_out, dim_j_out, dim_k_out,-2:3)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  Data_out (dim_i_out, dim_j_out, dim_k_out, number_of_inputs)
                      ! data interpolated to desired locations.

! Local Variables.

      Integer                                                           &
     &  i, j, k, index, n                                               &
                          ! Loop indices
     &, i_index, j_index


      Real                                                              &
     &  one_sixth                                                       &
     &, val_i_minus                                                     &
     &, val_i                                                           &
     &, val_i_plus                                                      &
     &, val_i_plus2                                                     &
     &, phi_cubed                                                       &
     &, phi_sq                                                          &
     &, lambda_cubed                                                    &
     &, lambda_sq                                                       &
                            ! x1-x6 are temporary variables
     &, x1, x3, x4, x5, x6                                              &
     &, coeff_i_minus(dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeff_i_zero(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_i_plus(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_i_plus2(dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeff_j_minus(dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeff_j_zero(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_j_plus(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_j_plus2(dim_i_out, dim_j_out, dim_k_out)

      Real                                                              &
     &  col_data (dim_i_out, dim_j_out, -1:2, number_of_inputs)
                                             ! Horizontally interpolated
                                             ! data ready for vertical
                                             ! interpolation.

! External Routines: None.

! ----------------------------------------------------------------------
! Section 0.   Calculate Horizontal interpolation weights
! ----------------------------------------------------------------------

      If ( L_regular ) Then

      one_sixth = 1./6.

! begin loop over levels

      Do k = 1, dim_k_out
        j=1
        Do i = 1, dim_i_out*dim_j_out

          phi_sq = weight_phi (i,j,k) * weight_phi (i,j,k)
          phi_cubed = weight_phi (i,j,k) * phi_sq

          coeff_j_plus2(i,j,k) = phi_cubed - weight_phi(i,j,k)
          coeff_j_plus(i,j,k)  =                                        &
     &             0.5 * (coeff_j_plus2(i,j,k) - phi_sq                 &
     &                     - weight_phi(i,j,k))
          coeff_j_zero(i,j,k)  = 0.5 * (coeff_j_plus2(i,j,k)            &
     &                             - 2. * phi_sq  + 2.)
          coeff_j_plus2(i,j,k) = one_sixth * coeff_j_plus2(i,j,k)
          coeff_j_minus(i,j,k) = one_sixth * (phi_cubed                 &
     &                                  - 3.*phi_sq +                   &
     &                                    2. * weight_phi(i,j,k))

          lambda_sq = weight_lambda (i,j,k) * weight_lambda (i,j,k)
          lambda_cubed = weight_lambda (i,j,k) * lambda_sq

          coeff_i_plus2(i,j,k)  = lambda_cubed - weight_lambda(i,j,k)
          coeff_i_plus(i,j,k)   =                                       &
     &       0.5 * (coeff_i_plus2(i,j,k) - lambda_sq                    &
     &                 - weight_lambda(i,j,k))
          coeff_i_zero(i,j,k)   = 0.5 * (coeff_i_plus2(i,j,k)           &
     &                        - 2.*lambda_sq + 2.)
          coeff_i_plus2(i,j,k)  = one_sixth * coeff_i_plus2(i,j,k)
          coeff_i_minus(i,j,k)  = one_sixth * (lambda_cubed             &
     &                                  - 3.*lambda_sq +                &
     &                                      2. * weight_lambda(i,j,k))

        End Do
      End Do

      else  ! variable resolution

      Do k = 1, dim_k_out

        j = 1
        Do i = 1, dim_i_out * dim_j_out

          x1 = weight_phi(i,j,k) + phi_rm(i_out(i,j,k),j_out(i,j,k))
          x3 = weight_phi(i,j,k) - 1.0
          x4 = x3 - phi_rp(i_out(i,j,k),j_out(i,j,k))
          x5 = x3 * x4
          x6 = x1 * weight_phi(i,j,k)

          coeff_j_plus2(i,j,k) = x6 * x3 *                              &
     &                           recip_phi_p2(i_out(i,j,k),j_out(i,j,k))
          coeff_j_plus(i,j,k)  = x6 * x4 *                              &
     &                           recip_phi_p(i_out(i,j,k),j_out(i,j,k))
          coeff_j_zero(i,j,k)  = x1 * x5 *                              &
     &                           recip_phi_0(i_out(i,j,k),j_out(i,j,k))
          coeff_j_minus(i,j,k) = weight_phi(i,j,k) * x5 *               &
     &                           recip_phi_m(i_out(i,j,k),j_out(i,j,k))

          x1 = weight_lambda(i,j,k) + lambda_rm(i_out(i,j,k))
          x3 = weight_lambda(i,j,k) - 1.0
          x4 = x3 - lambda_rp(i_out(i,j,k))
          x5 = x3 * x4
          x6 = x1 * weight_lambda(i,j,k)

          coeff_i_plus2(i,j,k) = x6 * x3 * recip_lambda_p2(i_out(i,j,k))
          coeff_i_plus(i,j,k)  = x6 * x4 * recip_lambda_p(i_out(i,j,k))
          coeff_i_zero(i,j,k)  = x1 * x5 * recip_lambda_0(i_out(i,j,k))
          coeff_i_minus(i,j,k) = weight_lambda(i,j,k) * x5 *            &
     &                                     recip_lambda_m(i_out(i,j,k))

        End Do  !  i = 1, dim_i_out * dim_j_out

      End Do  ! k = 1, dim_k_out

      end If ! L_regular
! loop over 4 levels required for cubic vertical interpolation
! ----------------------------------------------------------------------
! Section 1.   Perform cubic Lagrange Interpolation in j direction.
! ----------------------------------------------------------------------


! Interpolate data to get the data at the point for the k interpolation
      Do k = 1, dim_k_out
        Do n = 1, number_of_inputs
          j=1
          Do i = 1, dim_i_out*dim_j_out

!CDIR UNROLL=4
            Do index = -1, 2
                val_i_minus         = coeff_j_plus2(i,j,k) *            &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+2,                     &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_plus(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+1,                     &
     &               k_out(i,j,k)+index,n)                              &
     &                         + coeff_j_zero(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k),                       &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_minus(i,j,k) *                 &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)-1,                     &
     &               k_out(i,j,k)+index,n)

                val_i              = coeff_j_plus2(i,j,k) *             &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+2,                       &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_plus(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,                       &
     &               k_out(i,j,k)+index,n)                              &
     &                         + coeff_j_zero(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k),                         &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_minus(i,j,k) *                 &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)-1,                       &
     &               k_out(i,j,k)+index,n)

                val_i_plus         = coeff_j_plus2(i,j,k) *             &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+2,                     &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_plus(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,                     &
     &               k_out(i,j,k)+index,n)                              &
     &                         + coeff_j_zero(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k),                       &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_minus(i,j,k) *                 &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-1,                     &
     &               k_out(i,j,k)+index,n)

                val_i_plus2        = coeff_j_plus2(i,j,k) *             &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+2,                     &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_plus(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+1,                     &
     &               k_out(i,j,k)+index,n)                              &
     &                         + coeff_j_zero(i,j,k) *                  &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k),                       &
     &               k_out(i,j,k)+index,n)                              &
     &                         - coeff_j_minus(i,j,k) *                 &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)-1,                     &
     &               k_out(i,j,k)+index,n)

! Interpolate data to the 4 i points needed to do the i direction
! interpolation
! ----------------------------------------------------------------------
! Section 2.   Perform cubic Lagrange Interpolation in i direction.
! ----------------------------------------------------------------------

                col_data (i,j,index,n) = coeff_i_plus2(i,j,k)           &
     &                                   * val_i_plus2                  &
     &                                 - coeff_i_plus(i,j,k)            &
     &                                   * val_i_plus                   &
     &                                 + coeff_i_zero(i,j,k)            &
     &                                   * val_i                        &
     &                                 - coeff_i_minus(i,j,k)           &
     &                                   * val_i_minus
            End Do ! index = -1, 2
          End Do  !  i  = 1, dim_i_out*dim_j_out
        End Do      !n

! ----------------------------------------------------------------------
! Section 3.   Perform cubic Lagrange Interpolation in k direction.
! ----------------------------------------------------------------------

! Interpolate data
        Do n = 1, number_of_inputs
          j=1
            Do i = 1, dim_i_out*dim_j_out
              Data_out (i,j,k,n) = coeff_z(i,j,k,-1) *                  &
     &                             col_data (i,j,-1,n) +                &
     &                             coeff_z(i,j,k,0) *                   &
     &                             col_data (i,j,0,n) +                 &
     &                             coeff_z(i,j,k,1) *                   &
     &                             col_data (i,j,1,n) +                 &
     &                             coeff_z(i,j,k,2) *                   &
     &                             col_data (i,j,2,n)

            End Do !  i = 1, dim_i_out*dim_j_out
        End Do   !  n = number_of_inputs

      End Do   !  k = 1, dim_k_out

      return
      END SUBROUTINE Cubic_Lagrange

#endif
