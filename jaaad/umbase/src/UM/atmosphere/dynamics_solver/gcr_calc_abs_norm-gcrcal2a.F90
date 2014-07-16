#if defined(A10_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine GCR_calc_abs_norm

      subroutine GCR_calc_abs_norm(                                     &
     &                            Error, off_x, off_y,                  &
     &                            n_proc, model_domain,                 &
     &                            at_extremity,                         &
     &                            row_length, rows,                     &
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
! Date     Version     Comment
! ----     -------     -------
! 22/03/00 5.1        Changed j_start/stop for cyclic LAM
!                     Correct gc_imax to gc_rmax           Andy Malcolm
!LL   5.1   10/02/00  Use DOMTYP parameters                    P.Burton
!LL   5.2   27/09/00  tidy up code                           A.Malcolm
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
     &, model_domain

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


      Real                                                              &
     &  Error (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)

! Arguments with Intent OUT. ie: variables Output only
      Real                                                              &
     &  Abs_norm

! Local variables

      Integer i, j, k, info                                             &
     &, i_start, i_end                                                  &
     &, j_start, j_end

#include "parparm.h"
#include "domtyp.h"
!    External Routines:

      External                                                          &
     &  gc_rmax

! ----------------------------------------------------------------------
! Section 1.   Calculate maximum error over all processors
! ----------------------------------------------------------------------

      i_start = 1
      i_end = row_length
      j_start = 1
      j_end = rows
      If (model_domain  ==  mt_lam) Then
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_end = rows-1
        if(at_extremity(PEast)) i_end = row_length-1
        if(at_extremity(PWest)) i_start = 2
      Else If (model_domain  ==  mt_cyclic_lam) Then
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_end = rows-1
      End If

      Abs_norm= 0.0
      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = i_start, i_end
            If (error(i,j,k)  >   abs_norm ) abs_norm = error(i,j,k)
          End Do
        End Do
      End Do

! Calculate max over all processors

      Call gc_rmax(1, n_proc, info, abs_norm)

      Return
      END SUBROUTINE GCR_calc_abs_norm

#endif
