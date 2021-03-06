#if defined(A10_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_error_print

      subroutine GCR_error_print(                                       &
     &                            Error, gc_proc_row_group,             &
     &                            global_row_length, global_rows,       &
     &                            g_rows, g_row_length, g_datastart,    &
     &                            n_proc, n_procx, n_procy, me,         &
     &                            row_length, rows, model_levels,       &
     &                            i_start, i_stop, j_start, j_stop,     &
     &                            model_domain, solver_row_length,      &
     &                            init_error_mean, switch )

! Purpose:
!          Diagnostic print of residual error norms
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
! Version   Date       Comment
! ----     -------     -------
!   6.2   25/12/05  Change loop bounds variable names    Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, switch                                                          &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, gc_proc_row_group                                               &
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, me                                                              &
     &, model_domain                                                    &
     &, solver_row_length                                                      &
     &, i_start                                                         &
                      ! loop bound set in PE_Helmholtz
     &, i_stop                                                          &
                      ! loop bound set in PE_Helmholtz
     &, j_start                                                         &
                      ! loop bound set in PE_Helmholtz
     &, j_stop        ! loop bound set in PE_Helmholtz

      Integer                                                           &
     &  g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        &
     &, g_datastart(3,0:n_proc-1)

      Real                                                              &
     &  Error (row_length, rows, model_levels)

! Local variables
      Real                                                              &
     &  init_Error_mean (global_rows)                                   &
     &, error_sq(row_length, rows)                                      &
     &, mean(global_rows)                                               &
     &, rbuf(global_rows)

      Real                                                              &
     &  non_zero

      Integer d, i, j, k, istat, info, gj

#include "parparm.h"
#include "domtyp.h"
!    External Routines:
      External                                                          &
     &  gcg_rvecsumr, gc_rsend, gc_rrecv, gc_ssync

! ----------------------------------------------------------------------
! Section 1.   Calculate error norms
! ----------------------------------------------------------------------

      non_zero = 1.e-10

! If switch = 0 then set up initial error mean
! else calculate row error means.

      Do j = 1, rows
        Do i = 1, row_length
          error_sq(i,j) = 0.0
        End Do
      End Do
      Do k = 1, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            error_sq(i,j) = error_sq(i,j) +                             &
     &                      error(i,j,k) * error(i,j,k)
          End Do
        End Do
      End Do

      call gcg_rvecsumr(row_length,row_length,1,                        &
     &   rows,error_sq,gc_proc_row_group,istat,mean)

      Do j = 1, rows
        mean(j) = SQRT( mean(j) ) / (model_levels * solver_row_length)
      End Do

! Gather all data to processor 0 from its column and form
! global rows fields

! step 1 send data to processor zero
      If (me  /=  0) then
        Do j = 1, rows
          rbuf(j) = mean(j)
        End Do

        call gc_rsend(100*me, g_rows(me),                               &
     &                  0, info, rbuf, rbuf)
      End If

      call gc_ssync(n_proc,info)

! step 2 processor zero receives data and puts in correct location

      If (me  ==  0) then
        Do d = 1, n_procx*n_procy-1
          call gc_rrecv(100*d, g_rows(d),                               &
     &                  d, info, rbuf, rbuf)
          Do j = 1, g_rows(d)
            gj = g_datastart(2,d) + j - 1
            mean(gj) = rbuf(j)
          End Do
        End Do

        If (switch  ==  0) Then
          Do j = 1, global_rows
            init_Error_mean(j) = max(mean(j),non_zero)
          End Do
        Else
          write(6,*) ' End of convergence tolerances achieved '
          Do j = 1, global_rows
            mean(j) = mean(j) / init_Error_mean(j)
            write(6,100) j, mean(j)
          End Do
        End If

      End If

      call gc_ssync(n_proc,info)

 100  Format(1x,I3,1x,E10.2)

      return
      END SUBROUTINE GCR_error_print

#endif
