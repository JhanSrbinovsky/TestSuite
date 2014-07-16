
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Polar_Reset_Mean

      Subroutine Polar_Reset_Mean(                                      &
     &                      p, rho, theta, w, q, qcl, qcf,              &
     &                      cf, cfl, cff,                               &
     &                      row_length, rows, model_levels,             &
     &                      wet_model_levels, global_row_length,        &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      n_proc, n_procy, proc_row_group,            &
     &                      at_extremity)

! Purpose:
!          Resets polar values to mean over row to remove any
!          deviations due to rounding error.
!
! Method:
!          Is described in ;
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
!
! 5.3     19/10/01  Use appropriate gcg routines.   S. Cusack
! 5.4     27/07/02  Include cloud fractions for PC2 scheme. D. Wilson
! 6.2     21/10/05  Replace GSYNC with SSYNC. P.Selwood
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
     &, model_levels                                                    &
                         ! number of model levels.
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held
     &, global_row_length                                               &
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j     ! Size of halo in j direction.

      Integer                                                           &
     &  n_proc                                                          &
     &, n_procy                                                         &
     &, proc_row_group

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

! Arguments with Intent IN/OUT.

      Real                                                              &
     &  p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      0:model_levels)                                             &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)                                           &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      wet_model_levels)                                           &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_model_levels)                                         &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_model_levels)                                         &
     &, cf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &      wet_model_levels)                                           &
     &, cfl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_model_levels)                                         &
     &, cff(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_model_levels)

! Local Variables.

      Integer                                                           &
     &  i, k                                                            &
                  ! Loop indices
     &, info

      Real                                                              &
     &  l_f_poles (row_length,4*model_levels+6*wet_model_levels)        &
     &, sum(4*model_levels+6*wet_model_levels)

! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Reset polar values to mean over the row.
! ----------------------------------------------------------------------

! variables on all model levels.

      If (n_procy  >   1) Then

        If(at_extremity(PNorth)) Then
          Do k = 1, model_levels
            Do i = 1, row_length
              l_f_poles(i,k) = p(i,rows,k)
              l_f_poles(i,k+model_levels) = rho(i,rows,k)
              l_f_poles(i,k+2*model_levels) = theta(i,rows,k)
              l_f_poles(i,k+3*model_levels) = w(i,rows,k)
            End Do
          End Do
          Do k = 1, wet_model_levels
            Do i = 1, row_length
              l_f_poles(i,k+4*model_levels) = q(i,rows,k)
              l_f_poles(i,k+4*model_levels+wet_model_levels)            &
     &          = qcl(i,rows,k)
              l_f_poles(i,k+4*model_levels+2*wet_model_levels)          &
     &          = qcf(i,rows,k)
              l_f_poles(i,k+4*model_levels+3*wet_model_levels)          &
     &          = cf(i,rows,k)
              l_f_poles(i,k+4*model_levels+4*wet_model_levels)          &
     &          = cfl(i,rows,k)
              l_f_poles(i,k+4*model_levels+5*wet_model_levels)          &
     &          = cff(i,rows,k)
            End Do
          End Do

        Else If(at_extremity(PSouth)) Then
          Do k = 1, model_levels
            Do i = 1, row_length
              l_f_poles(i,k) = p(i,1,k)
              l_f_poles(i,k+model_levels) = rho(i,1,k)
              l_f_poles(i,k+2*model_levels) = theta(i,1,k)
              l_f_poles(i,k+3*model_levels) = w(i,1,k)
            End Do
          End Do
          Do k = 1, wet_model_levels
            Do i = 1, row_length
              l_f_poles(i,k+4*model_levels) = q(i,1,k)
              l_f_poles(i,k+4*model_levels+wet_model_levels)            &
     &          = qcl(i,1,k)
              l_f_poles(i,k+4*model_levels+2*wet_model_levels)          &
     &          = qcf(i,1,k)
              l_f_poles(i,k+4*model_levels+3*wet_model_levels)          &
     &          = cf(i,1,k)
              l_f_poles(i,k+4*model_levels+4*wet_model_levels)          &
     &          = cfl(i,1,k)
              l_f_poles(i,k+4*model_levels+5*wet_model_levels)          &
     &          = cff(i,1,k)
            End Do
          End Do

        Else
          Do k = 1, model_levels
            Do i = 1, row_length
              l_f_poles(i,k) = 0.
              l_f_poles(i,k+model_levels) = 0.
              l_f_poles(i,k+2*model_levels) = 0.
              l_f_poles(i,k+3*model_levels) = 0.
            End Do
          End Do
          Do k = 1, wet_model_levels
            Do i = 1, row_length
              l_f_poles(i,k+4*model_levels) = 0.
              l_f_poles(i,k+4*model_levels+wet_model_levels)            &
     &          = 0.
              l_f_poles(i,k+4*model_levels+2*wet_model_levels)          &
     &          = 0.
              l_f_poles(i,k+4*model_levels+3*wet_model_levels)          &
     &          = 0.
              l_f_poles(i,k+4*model_levels+4*wet_model_levels)          &
     &          = 0.
              l_f_poles(i,k+4*model_levels+5*wet_model_levels)          &
     &          = 0.
            End Do
          End Do
        End If

        call gc_ssync(n_proc,info)

        call gcg_rvecsumr(row_length,row_length,1,                      &
     &                    4*model_levels+6*wet_model_levels,            &
     &                    l_f_poles,proc_row_group,info,sum)






        Do i = 1, 4*model_levels+6*wet_model_levels
          sum(i) = sum(i) / global_row_length
        End Do

        If(at_extremity(PNorth)) Then
          Do k = 1, model_levels
            Do i = 1, row_length
              p(i,rows,k) = sum(k)
              rho(i,rows,k) = sum(k+model_levels)
              theta(i,rows,k) = sum(k+2*model_levels)
              w(i,rows,k) = sum(k+3*model_levels)
            End Do
          End Do
          Do k = 1, wet_model_levels
            Do i = 1, row_length
              q(i,rows,k) = sum(k+4*model_levels)
              qcl(i,rows,k) = sum(k+4*model_levels+wet_model_levels)
              qcf(i,rows,k) = sum(k+4*model_levels+2*wet_model_levels)
              cf(i,rows,k)  = sum(k+4*model_levels+3*wet_model_levels)
              cfl(i,rows,k) = sum(k+4*model_levels+4*wet_model_levels)
              cff(i,rows,k) = sum(k+4*model_levels+5*wet_model_levels)
            End Do
          End Do

        Else If(at_extremity(PSouth)) Then

          Do k = 1, model_levels
            Do i = 1, row_length
              p(i,1,k) = sum(k)
              rho(i,1,k) = sum(k+model_levels)
              theta(i,1,k) = sum(k+2*model_levels)
              w(i,1,k) = sum(k+3*model_levels)
            End Do
          End Do
          Do k = 1, wet_model_levels
            Do i = 1, row_length
              q(i,1,k) = sum(k+4*model_levels)
              qcl(i,1,k) = sum(k+4*model_levels+wet_model_levels)
              qcf(i,1,k) = sum(k+4*model_levels+2*wet_model_levels)
              cf(i,1,k)  = sum(k+4*model_levels+3*wet_model_levels)
              cfl(i,1,k) = sum(k+4*model_levels+4*wet_model_levels)
              cff(i,1,k) = sum(k+4*model_levels+5*wet_model_levels)
            End Do
          End Do

        End If

      Else

! One processor so must have both North and South
! Do north pole
        Do k = 1, model_levels
          Do i = 1, row_length
            l_f_poles(i,k) = p(i,rows,k)
            l_f_poles(i,k+model_levels) = rho(i,rows,k)
            l_f_poles(i,k+2*model_levels) = theta(i,rows,k)
            l_f_poles(i,k+3*model_levels) = w(i,rows,k)
          End Do
        End Do
        Do k = 1, wet_model_levels
          Do i = 1, row_length
            l_f_poles(i,k+4*model_levels) = q(i,rows,k)
            l_f_poles(i,k+4*model_levels+wet_model_levels)              &
     &        = qcl(i,rows,k)
            l_f_poles(i,k+4*model_levels+2*wet_model_levels)            &
     &        = qcf(i,rows,k)
            l_f_poles(i,k+4*model_levels+3*wet_model_levels)            &
     &        = cf(i,rows,k)
            l_f_poles(i,k+4*model_levels+4*wet_model_levels)            &
     &        = cfl(i,rows,k)
            l_f_poles(i,k+4*model_levels+5*wet_model_levels)            &
     &        = cff(i,rows,k)
          End Do
        End Do

        call gc_ssync(n_proc,info)

        call gcg_rvecsumr(row_length,row_length,1,                      &
     &                    4*model_levels+6*wet_model_levels,            &
     &                    l_f_poles,proc_row_group,info,sum)






        Do i = 1, 4*model_levels+6*wet_model_levels
          sum(i) = sum(i) / global_row_length
        End Do

        Do k = 1, model_levels
          Do i = 1, row_length
            p(i,rows,k) = sum(k)
            rho(i,rows,k) = sum(k+model_levels)
            theta(i,rows,k) = sum(k+2*model_levels)
            w(i,rows,k) = sum(k+3*model_levels)
          End Do
        End Do
        Do k = 1, wet_model_levels
          Do i = 1, row_length
            q(i,rows,k) = sum(k+4*model_levels)
            qcl(i,rows,k) = sum(k+4*model_levels+wet_model_levels)
            qcf(i,rows,k) = sum(k+4*model_levels+2*wet_model_levels)
            cf(i,rows,k)  = sum(k+4*model_levels+3*wet_model_levels)
            cfl(i,rows,k) = sum(k+4*model_levels+4*wet_model_levels)
            cff(i,rows,k) = sum(k+4*model_levels+5*wet_model_levels)
          End Do
        End Do

! Do south pole
        Do k = 1, model_levels
          Do i = 1, row_length
            l_f_poles(i,k) = p(i,1,k)
            l_f_poles(i,k+model_levels) = rho(i,1,k)
            l_f_poles(i,k+2*model_levels) = theta(i,1,k)
            l_f_poles(i,k+3*model_levels) = w(i,1,k)
          End Do
        End Do
        Do k = 1, wet_model_levels
          Do i = 1, row_length
            l_f_poles(i,k+4*model_levels) = q(i,1,k)
            l_f_poles(i,k+4*model_levels+wet_model_levels)              &
     &          = qcl(i,1,k)
            l_f_poles(i,k+4*model_levels+2*wet_model_levels)            &
     &          = qcf(i,1,k)
            l_f_poles(i,k+4*model_levels+3*wet_model_levels)            &
     &          = cf(i,1,k)
            l_f_poles(i,k+4*model_levels+4*wet_model_levels)            &
     &          = cfl(i,1,k)
            l_f_poles(i,k+4*model_levels+5*wet_model_levels)            &
     &          = cff(i,1,k)
          End Do
        End Do

        call gc_ssync(n_proc,info)

        call gcg_rvecsumr(row_length,row_length,1,                      &
     &                    4*model_levels+6*wet_model_levels,            &
     &                    l_f_poles,proc_row_group,info,sum)






        Do i = 1, 4*model_levels+6*wet_model_levels
          sum(i) = sum(i) / global_row_length
        End Do

        Do k = 1, model_levels
          Do i = 1, row_length
            p(i,1,k) = sum(k)
            rho(i,1,k) = sum(k+model_levels)
            theta(i,1,k) = sum(k+2*model_levels)
            w(i,1,k) = sum(k+3*model_levels)
          End Do
        End Do
        Do k = 1, wet_model_levels
          Do i = 1, row_length
            q(i,1,k) = sum(k+4*model_levels)
            qcl(i,1,k) = sum(k+4*model_levels+wet_model_levels)
            qcf(i,1,k) = sum(k+4*model_levels+2*wet_model_levels)
            cf(i,1,k)  = sum(k+4*model_levels+3*wet_model_levels)
            cfl(i,1,k) = sum(k+4*model_levels+4*wet_model_levels)
            cff(i,1,k) = sum(k+4*model_levels+5*wet_model_levels)
          End Do
        End Do

      End If

! End of routine
      return
      END SUBROUTINE Polar_Reset_Mean

