#if defined(A10_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_abs_print

      subroutine GCR_abs_print(                                         &
     &                           Error, off_x, off_y,                   &
     &                            gc_proc_row_group,                    &
     &                            global_row_length, global_rows,       &
     &                            g_rows, g_row_length, g_datastart,    &
     &                            n_proc, n_procx, n_procy, me,         &
     &                            row_length, rows,                     &
     &                            model_levels, model_domain,           &
     &                            solver_row_length,                    &
     &                            i_start, i_stop, j_start, j_stop )

! Purpose:
!          Diagnostic print of absolute error norms
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
! Version  Date        Comment
! ----     -------     -------
!   6.2   20/10/05  Change loop bounds variable names    Terry Davies
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
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, gc_proc_row_group                                               &
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, me                                                              &
     &, off_x                                                           &
     &, off_y                                                           &
     &, model_domain                                                    &
     &, solver_row_length                                                      &
     &, i_start, i_stop                                                 &
                                   ! loop bounds set in PE_Helmholtz
     &, j_start, j_stop            ! loop bounds set in PE_Helmholtz

      Integer                                                           &
     &  g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        &
     &, g_datastart(3,0:n_proc-1)

      Real                                                              &
     &  Error (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)

! Local variables
      Real                                                              &
     &  error_sum(row_length, rows)                                     &
     &, mean(global_rows)                                               &
     &, rbuf(global_rows)                                               &
     &, max_val(model_levels)

      Integer d, i, j, k, istat, info, gj

#include "parparm.h"
#include "domtyp.h"
!    External Routines:
      External                                                          &
     &  gcg_rvecsumr, gc_imax, gc_rsend, gc_rrecv, gc_ssync

! ----------------------------------------------------------------------
! Section 1.   Calculate error norm over rows and levels.
! ----------------------------------------------------------------------


      Do k = 1, model_levels
        max_val(k) = 0.0
      End Do
      k = 1
      Do j = 1, rows
        Do i = 1, row_length
          error_sum(i,j) = 0.0
        End Do
      End Do
      Do k = 1, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            error_sum(i,j) = error_sum(i,j) + error(i,j,k)
            If (error(i,j,k)  >   max_val(k) )                          &
     &        max_val(k) = error(i,j,k)
          End Do
        End Do
      End Do

#if defined(REPROD)
      call gcg_rvecsumr(row_length, row_length, 1, rows,                &
     &                  error_sum, gc_proc_row_group, istat, mean)
#else
      call gcg_rvecsumf(row_length, row_length, 1, rows,                &
     &                  error_sum, gc_proc_row_group, istat, mean)
#endif

      Do j = 1, rows
        mean(j) = mean(j) / ( model_levels * solver_row_length )
      End Do

! Gather all data to processor 0 from its column and form
! global rows fields

      Call gc_imax(model_levels, n_proc, info, max_val)

! step 1 send data to processor zero
      If (me /= 0) then
        Do j = 1, rows
          rbuf(j) = mean(j)
        End Do

        call gc_rsend( 100 * me, g_rows(me), 0, info, rbuf, rbuf )
      End If ! me /= 0

      call gc_ssync(n_proc,info)

! step 2 processor zero receives data and puts in correct location

      If (me == 0) then
        Do d = 1, n_procx*n_procy-1
          call gc_rrecv(100*d, g_rows(d),                               &
     &                  d, info, rbuf, rbuf)
          Do j = 1, g_rows(d)
            gj = g_datastart(2,d) + j - 1
            mean(gj) = rbuf(j)
          End Do
        End Do

        write(6,*)' mean approximate absolute density change left '
        Do j = 1, global_rows
          write(6,100) j, mean(j)
        End Do
        write(6,*)' max value left per level '
        Do k = 1, model_levels
            write(6,100) k, max_val(k)
        End Do

      End If !  me == 0

      call gc_ssync(n_proc,info)

 100  Format(1x,I3,1x,E10.2)

      return
      END SUBROUTINE GCR_abs_print

#endif
