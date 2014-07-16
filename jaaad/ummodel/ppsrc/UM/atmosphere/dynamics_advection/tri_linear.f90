
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Tri_Linear

      Subroutine Tri_Linear(                                            &
     &                      Ext_Data,                                   &
     &                      dim_i_in, dim_j_in, dim_k_in,               &
     &                      dim_i_out, dim_j_out, dim_k_out,            &
     &                      halo_i, halo_j, number_of_inputs,           &
     &                      weight_lambda, weight_phi,                  &
     &                      i_out, j_out, k_out,                        &
     &                      coeff_z_lin,                                &
     &                      Data_out)

! Purpose:
!          Performs linear interpolation of the input field to a
!          set of points given by weight_lambda,weight_phi,r_out.
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
!  5.1  15/02/00   Upgrade from v2p7 to v2p8a           Andy Malcolm
!  6.0  18/08/03  NEC SX-6 optimisation - R Barnes & J-C Rioual.
!  6.1  17/05/04 change and optimise code for super-tracer arrays and
!                remove inconsistency if n_inputs was 1     Andy Malcolm
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
     &, coeff_z_lin(dim_i_out, dim_j_out, dim_k_out, 0:1)               &
     &, weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
                                                        ! a number
                                                      ! between 0 & 1
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)  ! a number between
                                                      ! 0 & 1

      Integer                                                           &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! point such that
     &, j_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! the desired
                                                      ! output point
     &, k_out (dim_i_out, dim_j_out, dim_k_out)       ! lies between it
                                                      ! and it+1

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  Data_out (dim_i_out, dim_j_out, dim_k_out, number_of_inputs)
                      ! data interpolated to desired locations.

! Local Variables.

      Integer                                                           &
     &  i, j, k, n   ! Loop indices

      Real                                                              &
     &  val_here                                                        &
     &, val_here_plus                                                   &
     &, a_coeff( dim_i_out, dim_j_out, dim_k_out )                      &
     &, b_coeff( dim_i_out, dim_j_out, dim_k_out )                      &
     &, c_coeff( dim_i_out, dim_j_out, dim_k_out )                      &
     &, d_coeff( dim_i_out, dim_j_out, dim_k_out )

! External Routines: None.

! ----------------------------------------------------------------------
! Section 1.   Perform Linear Interpolation horizontal directions
!              simultaneously.
! ----------------------------------------------------------------------

        Do k = 1, dim_k_out
          j=1
          Do i = 1, dim_i_out * dim_j_out
            d_coeff(i,j,k) = weight_lambda(i,j,k)*weight_phi(i,j,k)
            c_coeff(i,j,k) = weight_phi(i,j,k) - d_coeff(i,j,k)
            b_coeff(i,j,k) = weight_lambda(i,j,k) - d_coeff (i,j,k)
            a_coeff(i,j,k) = 1.-weight_lambda(i,j,k)-c_coeff(i,j,k)
          End Do
!         End Do
        End Do ! end loop over levels

        Do n = 1, number_of_inputs
          Do k = 1, dim_k_out
            j=1
            Do i = 1, dim_i_out * dim_j_out

              val_here = a_coeff(i,j,k) *                               &
     &        Ext_Data (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k),n)       &
     &                  +  b_coeff(i,j,k)                               &
     &      * Ext_Data (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k),n)     &
     &                  +  c_coeff(i,j,k)                               &
     &      * Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k),n)     &
     &                  +  d_coeff(i,j,k) *                             &
     &        Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k),n)

              val_here_plus = a_coeff(i,j,k)                            &
     &      * Ext_Data (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+1,n)     &
     &                  +  b_coeff(i,j,k)                               &
     &      * Ext_Data (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+1,n)   &
     &                  +  c_coeff(i,j,k)                               &
     &      * Ext_Data (i_out(i,j,k),j_out(i,j,k)+1,k_out(i,j,k)+1,n)   &
     &                  +  d_coeff(i,j,k) *                             &
     &        Ext_Data (i_out(i,j,k)+1,j_out(i,j,k)+1,k_out(i,j,k)+1,n)

! interpolate data vertically
                Data_out (i,j,k,n) = coeff_z_lin(i,j,k,1)               &
     &                         * val_here_plus                          &
     &                         + coeff_z_lin(i,j,k,0) * val_here

              End Do
!         End Do

          End Do ! end loop over levels
        End Do ! end loop over number of inputs

! End of routine.
      return
      END SUBROUTINE Tri_Linear

