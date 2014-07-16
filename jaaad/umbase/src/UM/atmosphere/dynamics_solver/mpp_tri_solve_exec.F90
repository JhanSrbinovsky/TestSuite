#if defined(A10_2A) || defined(A10_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine Mpp_tri_solve_exec

      Subroutine Mpp_tri_solve_exec(                                    &
     &                     row_length, rows, model_levels,              &
     &                     off_x, off_y, j_start, j_end,                &
     &                     n_proc, n_procx, me,                         &
     &                     proc_row_group,                              &
     &                     a_plus_x, a_minus_x,                         &
     &                     factor_forward, factor_backward,             &
     &                     recip_a_central_x,                           &
     &                     bv_a_matrix_0, bv_a_matrix_np1,              &
     &                     recip_bv_a_matrix_diag, bv_a_matrix_sup,     &
     &                     bv_factor_forward, bv_factor_backward,       &
     &                     bv_soln_n_term,                              &
     &                     bv_soln_1_term1,                             &
     &                     bv_soln_1_term2,                             &
     &                     rhs_xig, sss, soln)

! Purpose:
!          Calculates ADI pre-conditioning matrix coefficients.
!          This code only for global model at the moment.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Original Progammer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Version    Date   Comment
! ----     -------  -------
!  6.0    18/08/03  NEC SX-6 optimisation - change dimension to avoid
!                   bank conflict etc.  R Barnes & J-C Rioual.
!  6.2    03.11.04  Loop collasing to get sufficent vector length in case
!                   of high number of CPUs.   Klaus Ketelsen/R Barnes
!  6.2    10/03/06  Add section 10_2B           A. Malcolm
!  6.4    26/01/07  Make all loops in j have same range        A Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit none

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, me                                                              &
     &, proc_row_group                                                  &
     &, off_x                                                           &
     &, off_y                                                           &
     &, j_start, j_end

      real con1,con2
      Real                                                              &
     &  a_plus_x(row_length+1,rows,model_levels)                        &
     &, a_minus_x(row_length+1,rows,model_levels)                       &
     &, rhs_xig(row_length,rows,model_levels)

!kk   Dimension row_length+1 prevents copying to sss_x
      real,intent(in),                                                  &
     &     dimension(row_length+1,j_start:j_end,model_levels) :: sss

      Real                                                              &
     &  factor_forward(row_length+1,j_start:j_end,model_levels)         &
     &, factor_backward(row_length+1,j_start:j_end,model_levels)        &
     &, bv_a_matrix_0(2*n_procx,rows,model_levels)                      &
     &, bv_a_matrix_np1(2*n_procx,rows,model_levels)                    &
     &, bv_factor_forward(2,2*n_procx,rows,model_levels)                &
     &, bv_factor_backward(2*n_procx,rows,model_levels)                 &
     &, recip_bv_a_matrix_diag(2*n_procx,rows,model_levels)             &
     &, bv_a_matrix_sup(2*n_procx,rows,model_levels)                    &
     &, recip_a_central_x(row_length,rows,model_levels)

      Real                                                              &
     &  bv_soln_n_term(rows, model_levels)                              &
     &, bv_soln_1_term1(rows, model_levels)                             &
     &, bv_soln_1_term2(rows, model_levels)

! Arguments with Intent OUT. ie: Output variables.
      Real                                                              &
     &  soln(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)

! Local Variables.

      Integer                                                           &
     &  i,j,k,d                                                         &
     &, nbv_points                                                      &
     &, ime                                                             &
     &, ibase                                                           &
     &, info                                                            &
     &, len                                                             &
     &, start_point

! Local arrays.
      Real                                                              &
     &  bv_rhs_x(2*n_procx,rows,model_levels)                           &
     &, bv_soln(0:2*n_procx+1,rows,model_levels)                        &
     &, rbuf(2,rows,model_levels,0:n_procx-1)

      Real rhs_xx(row_length+1,j_start:j_end,model_levels)
      Real sss_x(row_length+1,j_start:j_end,model_levels)

!    External routines.

!-----------------------------------------------------------------------
!     Section 1. Forward Gaussian elimination sweep
!-----------------------------------------------------------------------
!     Create a local copy of array sss to avoid bank conflicts
!     in the subsequent  elimination sweep
!     The overhead of copying is much cheaper than the bankconflicts
!     However a better solution would be a thorough redesign of the
!     whole structure of the solver logic to avoid bad memory access
!     patterns by proper dimensioning of all arrays.


      Do k = 1, model_levels
        Do j = j_start, j_end
        rhs_xx(1,j,k)=sss(1,j,k)
        End Do
      End Do

!CDIR unroll=8
          Do i = 2, row_length
!CDIR collapse
      Do k = 1, model_levels
        Do j = j_start, j_end
! update right-hand-side term
            rhs_xx(i,j,k) = sss(i,j,k) + factor_forward(i,j,k) *        &
     &                                    rhs_xx(i-1,j,k)
          End Do
        End Do
      End Do

!CDIR unroll=8
          Do i = row_length - 1, 1, -1
!CDIR collapse
      Do k = 1, model_levels
        Do j = j_start, j_end
! update right-hand-side term
            rhs_xx(i,j,k) = rhs_xx(i,j,k) + factor_backward(i,j,k) *    &
     &                                    rhs_xx(i+1,j,k)
          End Do
        End Do
      End Do

!-----------------------------------------------------------------------
!     Section 3. Solve for boundary points
!-----------------------------------------------------------------------
! calculate position on row of processors, 0 = start, n_procx-1 = end
      ibase = (me/n_procx) * n_procx
      ime = me - ibase

! number of boundary data points is twice number of processors
      nbv_points = 2 * n_procx

! copy data to be broadcast from this processor into the buffer
      Do k = 1, model_levels
        Do j = j_start, j_end
          rbuf(1,j,k,ime) = rhs_xx(1,j,k)
          rbuf(2,j,k,ime) = rhs_xx(row_length,j,k)
        End Do
      End Do

! broadcast buffer data to all processors on this row.
      len = model_levels * rows * 2

      if(n_procx > 1)  then
        Call gcg_ssync(proc_row_group, info)
        Do k = 0, n_procx-1
          Call gcg_rbcast(200+k, len, ibase+k, proc_row_group,          &
     &                    info, rbuf(1,1,1,k) )
        End Do
        Call gcg_ssync(proc_row_group, info)
      end if

! Copy all boundary data values into new matrix
      Do d = 0, n_procx-1
        i = 2 * d
        Do k = 1, model_levels
          Do j = j_start, j_end
            bv_rhs_x(i+1,j,k) = rbuf(1,j,k,d)
            bv_rhs_x(i+2,j,k) = rbuf(2,j,k,d)
          End Do
        End Do
      End Do

!-----------------------------------------------------------------------
!     Section 3.1. Forward Gaussian elimination sweep on bv rhs
!-----------------------------------------------------------------------

      Do k = 1, model_levels
        Do i = j_start, j_end
        Do j = 2, nbv_points - 2, 2
! elimination from data on even rows
          bv_rhs_x(j+1,i,k) = bv_rhs_x(j+1,i,k)                         &
     &                      + bv_factor_forward(1,j,i,k)                &
     &                      * bv_rhs_x(j,i,k)
          bv_rhs_x(j+2,i,k) = bv_rhs_x(j+2,i,k)                         &
     &                      + bv_factor_forward(2,j,i,k)                &
     &                      * bv_rhs_x(j,i,k)

! elimination from data on odd rows
          bv_rhs_x(j+2,i,k) = bv_rhs_x(j+2,i,k)                         &
     &                      + bv_factor_forward(1,j+1,i,k)              &
     &                      * bv_rhs_x(j+1,i,k)
        End Do
        End Do
      End Do

!-----------------------------------------------------------------------
!     Section 3.2. Backward Gaussian elimination sweep on bv rhs
!-----------------------------------------------------------------------

      Do k = 1, model_levels
      Do i = j_start, j_end
        Do j = nbv_points - 1, 3, -2
           bv_rhs_x(j-2,i,k) = bv_rhs_x(j-2,i,k)                        &
     &                         + bv_factor_backward(j,i,k)              &
     &                         * bv_rhs_x(j,i,k)
        End Do
      End Do
      End Do

!-----------------------------------------------------------------------
!     Section 3.3. Solution for end points 1 and nbv_points
!-----------------------------------------------------------------------

      Do k = 1, model_levels
      Do i = j_start, j_end
        bv_soln(1,i,k) = ( bv_rhs_x(1,i,k)                              &
     &                     - bv_soln_1_term1(i,k) *                     &
     &                       bv_rhs_x(nbv_points,i,k) )                 &
     &                     * bv_soln_1_term2(i,k)
        bv_soln(nbv_points,i,k) = ( bv_rhs_x(nbv_points,i,k) -          &
     &                                bv_a_matrix_np1(nbv_points,i,k) * &
     &                                bv_soln(1,i,k) ) *                &
     &                              bv_soln_n_term(i,k)
      End Do
      End Do

! set 0, and nbv_points+1 solutions, ie make solution periodic
      Do k = 1, model_levels
        Do j = j_start, j_end
          bv_soln(0,j,k) = bv_soln(nbv_points,j,k)
          bv_soln(nbv_points+1,j,k) = bv_soln(1,j,k)
        End Do
      End Do

!-----------------------------------------------------------------------
!     Section 3.4. Find solution in inner of matrix
!-----------------------------------------------------------------------

      Do k = 1, model_levels
        Do j = j_start, j_end
! first all odd rows
        Do i = 3, nbv_points - 1, 2
          bv_soln(i,j,k) = ( bv_rhs_x(i,j,k)                            &
     &                     -   bv_a_matrix_0(i,j,k) * bv_soln(0,j,k)    &
     &                     -   bv_a_matrix_np1(i,j,k) *                 &
     &                         bv_soln(nbv_points+1,j,k) )              &
     &                     * recip_bv_a_matrix_diag(i,j,k)
        End Do
! now all even rows
        Do i = 2, nbv_points - 2, 2
          bv_soln(i,j,k) = ( bv_rhs_x(i,j,k)                            &
     &                     -   bv_a_matrix_sup(i,j,k) * bv_soln(i+1,j,k)&
     &                     -   bv_a_matrix_0(i,j,k) * bv_soln(0,j,k)    &
     &                     -   bv_a_matrix_np1(i,j,k) *                 &
     &                         bv_soln(nbv_points+1,j,k) )              &
     &                     * recip_bv_a_matrix_diag(i,j,k)
        End Do
        End Do
      End Do

!-----------------------------------------------------------------------
!     Section 3.5.  Copy bv_soln into main solution
!-----------------------------------------------------------------------

      start_point = 2 * ime
      Do k = 1, model_levels
        Do j = j_start, j_end
          soln(0,j,k) = bv_soln(start_point,j,k)
          soln(1,j,k) = bv_soln(start_point+1,j,k)
          soln(row_length,j,k) = bv_soln(start_point+2,j,k)
          soln(row_length+1,j,k) = bv_soln(start_point+3,j,k)
        End Do
      End Do

!-----------------------------------------------------------------------
!     Section 4. Solve for main solution not at boundaries
!-----------------------------------------------------------------------

      Do k = 1, model_levels
        Do j = j_start, j_end
        con1=soln(0,j,k)
        con2=soln(row_length+1,j,k)
!dir$ unroll
          Do i = 2, row_length-1
            soln(i,j,k) = ( rhs_xx(i,j,k)                               &
     &                         - a_minus_x(i,j,k) * con1                &
     &                         - a_plus_x(i,j,k) * con2 )               &
     &                     * recip_a_central_x(i,j,k)
          End Do
        End Do
      End Do

!     end of routine mpp_tri_solve_exec

      return
      END SUBROUTINE Mpp_tri_solve_exec

#endif
