#if defined(A10_2B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine GCR_calc_abs_norm

      subroutine GCR_calc_abs_norm(                                     &
     &                            Error, off_x, off_y,                  &
     &                            i_start, i_stop, j_start, j_stop,     &
     &                            n_proc, row_length, rows,             &
     &                            model_levels, Abs_Norm)

! Purpose:
!          Calculate max value of error over all processors
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
!   6.2   25/12/05  Pass loop bounds variables as arguments
!                                                        Terry Davies
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
     &, n_proc                                                          &
     &, off_x                                                           &
     &, off_y                                                           &
     &, i_start, i_stop                                                 &
                                   ! loop bounds set in PE_Helmholtz
     &, j_start, j_stop            ! loop bounds set in PE_Helmholtz

      Real                                                              &
     &  Error (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)

! Arguments with Intent OUT. ie: variables Output only
      Real                                                              &
     &  Abs_norm

! Local variables

      Integer i, j, k, info

#include "parparm.h"
#include "domtyp.h"
!    External Routines:

      External                                                          &
     &  gc_rmax

! ----------------------------------------------------------------------
! Section 1.   Calculate maximum error over all processors
! ----------------------------------------------------------------------


      Abs_norm= 0.0
      Do k = 1, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            If (error(i,j,k)  >   abs_norm ) abs_norm = error(i,j,k)
          End Do
        End Do
      End Do

! Calculate max over all processors

      Call gc_rmax(1, n_proc, info, abs_norm)

      Return
      END SUBROUTINE GCR_calc_abs_norm

#endif
