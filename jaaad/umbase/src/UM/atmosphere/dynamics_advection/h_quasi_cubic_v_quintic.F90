#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!! Subroutine h_quasi_cubic_v_quintic
      Subroutine h_quasi_cubic_v_quintic                                &
     &                         (Ext_Data,                               &
     &                          dim_i_in, dim_j_in, dim_k_in,           &
     &                          dim_i_out, dim_j_out, dim_k_out,        &
     &                          halo_i, halo_j, number_of_inputs,       &
     &                          weight_lambda, weight_phi,              &
     &                          i_out, j_out, k_out,                    &
     &                          coeff_z,                                &
     &                          Data_out)

! Purpose:
!          Performs non-monotone quasi-cubic interpolation, as developed
!          by Mariano Hortal at ECMWF, of the input field to a
!          set of points defined by i_out, j_out, k_out, and
!          weight_lambda,weight_phi.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Original Progammer: Andrew J. Malcolm
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date     Comment
! -------  -------   -------
!  5.3     11/07/01  Original code                         Andy Malcolm
!  6.0  19/08/03  NEC SX-6 optimisation - loop collapse and
!                 unroll directives.  R Barnes & J-C Rioual.
!  6.0  18/08/03  NEC SX-6 optimisation - R Barnes & J-C Rioual.
!  6.1  07/07/04  Change code to be faster, especially for tracer
!                 super_arrays.                        Andy Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
! System component covered: ??
! System Task:              ??
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

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
     &, number_of_inputs ! number of fields to interpolate

      Real                                                              &
     &  Ext_Data (1-halo_i:dim_i_in+halo_i+1,                           &
     &            1-halo_j:dim_j_in+halo_j,                             &
                                                   !data to be
     &           -1:dim_k_in+2, number_of_inputs)                       &
                                                   ! interpolated
     &, weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
                                                        ! a number
                                                      ! between 0 & 1
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)  ! a number
                                                      ! between 0 & 1
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

! Local Parameters

      Real one_sixth
      Parameter ( one_sixth = 1./6. )

! Local Variables.

      Integer                                                           &
     &  i, j, k, index, n ! Loop indices

      Real                                                              &
     &  val_i_minus                                                     &
     &, val_i                                                           &
     &, val_i_plus                                                      &
     &, val_i_plus2                                                     &
     &, phi_cubed                                                       &
     &, phi_sq                                                          &
     &, lambda_cubed                                                    &
     &, lambda_sq                                                       &
     &, coeff_i_minus(dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeff_i_zero(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_i_plus(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_i_plus2(dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeff_j_minus(dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeff_j_zero(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_j_plus(dim_i_out, dim_j_out, dim_k_out)                   &
     &, coeff_j_plus2(dim_i_out, dim_j_out, dim_k_out)                  &
     &, coeff_j_lin1(dim_i_out, dim_j_out, dim_k_out)

      Real                                                              &
     &  col_data (dim_i_out, dim_j_out, -2:3, number_of_inputs)
                                           ! Horizontally interpolated
                                           ! data ready for vertical
                                           ! interpolation.

! External Routines: None.

! ----------------------------------------------------------------------
! Section 1.   Perform cubic Lagrange Interpolation in j direction and
!              then in i direction, for central lines, linear for
!              outside lines, and ensure all pieces are monotone.
! ----------------------------------------------------------------------

!     begin loop over levels
      Do k = 1, dim_k_out
!        Do j = 1, dim_j_out
        j=1
        Do i = 1, dim_i_out * dim_j_out
!         interpolate in j at each index level for the four points
!         needed for the i interpolation.

                phi_sq    = weight_phi (i,j,k) * weight_phi (i,j,k)
                phi_cubed = weight_phi (i,j,k) * phi_sq

                coeff_j_plus2(i,j,k) = phi_cubed - weight_phi(i,j,k)
                coeff_j_plus(i,j,k)  = 0.5 * (coeff_j_plus2(i,j,k) -    &
     &                             phi_sq  - weight_phi(i,j,k))
                coeff_j_zero(i,j,k)  = 0.5 *                            &
     &                      (coeff_j_plus2(i,j,k) - 2. * phi_sq + 2.)
                coeff_j_minus(i,j,k) = one_sixth * (phi_cubed           &
     &                                - 3.*phi_sq +                     &
     &                                    2. * weight_phi(i,j,k))
                coeff_j_plus2(i,j,k) = one_sixth * coeff_j_plus2(i,j,k)
                coeff_j_lin1(i,j,k)  = 1.- weight_phi(i,j,k)

! interpolate in i and store.
                lambda_sq = weight_lambda(i,j,k) * weight_lambda(i,j,k)
                lambda_cubed = weight_lambda (i,j,k) * lambda_sq

                coeff_i_plus2(i,j,k) = (lambda_cubed -                  &
     &                                  weight_lambda(i,j,k))
                coeff_i_plus(i,j,k)  = 0.5 * (coeff_i_plus2(i,j,k)      &
     &                         - lambda_sq  - weight_lambda(i,j,k))
                coeff_i_zero(i,j,k)  = 0.5 *                            &
     &                          (coeff_i_plus2(i,j,k) - 2.*lambda_sq    &
     &                                + 2.)
                coeff_i_plus2(i,j,k) = one_sixth * coeff_i_plus2(i,j,k)
                coeff_i_minus(i,j,k) = one_sixth *                      &
     &                          (lambda_cubed - 3.*lambda_sq +          &
     &                                      2. * weight_lambda(i,j,k))
        End Do
      End Do

      Do k = 1, dim_k_out
        Do n = 1, number_of_inputs
          Do index = -2, 3
!           Do j = 1, dim_j_out
            j=1
            Do i = 1, dim_i_out * dim_j_out

!CDIR UNROLL=6
                val_i_minus = weight_phi(i,j,k) *                       &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)+1,k_out(i,j,k)+index,n)&
     &                    + coeff_j_lin1(i,j,k) *                       &
     &     Ext_Data (i_out(i,j,k)-1,j_out(i,j,k)  ,k_out(i,j,k)+index,n)

                val_i =  coeff_j_plus2(i,j,k) *                         &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+2,k_out(i,j,k)+index,n)  &
     &                  - coeff_j_plus(i,j,k) *                         &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k)+index,n)  &
     &                  + coeff_j_zero(i,j,k) *                         &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)  ,k_out(i,j,k)+index,n)  &
     &                  - coeff_j_minus(i,j,k) *                        &
     &     Ext_Data (i_out(i,j,k),j_out(i,j,k)-1,k_out(i,j,k)+index,n)

                val_i_plus =   coeff_j_plus2(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+2,k_out(i,j,k)+index,n)&
     &                       - coeff_j_plus(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k)+index,n)&
     &                       + coeff_j_zero(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)  ,k_out(i,j,k)+index,n)&
     &                       - coeff_j_minus(i,j,k) *                   &
     &     Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)-1,k_out(i,j,k)+index,n)

                val_i_plus2 =  weight_phi(i,j,k) *                      &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)+1,k_out(i,j,k)+index,n)&
     &                       + coeff_j_lin1(i,j,k) *                    &
     &     Ext_Data (i_out(i,j,k)+2,j_out(i,j,k)  ,k_out(i,j,k)+index,n)


                col_data (i,j,index,n) =                                &
     &                               coeff_i_plus2(i,j,k) * val_i_plus2 &
     &                             - coeff_i_plus(i,j,k)  * val_i_plus  &
     &                             + coeff_i_zero(i,j,k)  * val_i       &
     &                             - coeff_i_minus(i,j,k) * val_i_minus
              End Do
!           End Do

          End Do ! end loop over index
        End Do ! end loop over number of inputs
! End loop over levels.

! ----------------------------------------------------------------------
! Section 2.   Perform quintic Lagrange Interpolation in k direction.
! ----------------------------------------------------------------------

! Interpolate data
        Do n = 1, number_of_inputs
!         Do j = 1, dim_j_out
          j=1
            Do i = 1, dim_i_out * dim_j_out
              Data_out(i,j,k,n) = coeff_z(i,j,k,-2) *col_data (i,j,-2,n)&
     &                          + coeff_z(i,j,k,-1) *col_data (i,j,-1,n)&
     &                          + coeff_z(i,j,k,0)  *col_data (i,j,0,n) &
     &                          + coeff_z(i,j,k,1)  *col_data (i,j,1,n) &
     &                          + coeff_z(i,j,k,2)  *col_data (i,j,2,n) &
     &                          + coeff_z(i,j,k,3)  *col_data (i,j,3,n)

            End Do
!         End Do
        End Do
! End loop over levels.
      End Do

! End of routine.
      return
      END SUBROUTINE h_quasi_cubic_v_quintic
#endif
