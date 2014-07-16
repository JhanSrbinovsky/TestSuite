
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE EXCFNL_CCI---------------------------------------------
!!!
!!!  Purpose: Compute Compressed Index
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  6.2    Jan 2006  New deck as part of optimisation of EXCF_NL
!!!                   by J-C Rioual.                 A P Lock
!!!
!!!  Code Description:
!!!    Language: FORTRAN 77 + common extensions.
!!!    This code is written to UMDP3 v6 programming standards.
!!!
!!!  Documentation: UMDP No.24
!!!
!!!--------------------------------------------------------------------
!
!!   Arguments :-
      SUBROUTINE EXCFNL_CCI (                                           &
     & row_length, rows, c_len, to_do, ind_todo                         &
     & )
!
      IMPLICIT NONE

        integer, intent(in) :: row_length, rows

        integer, intent(inout) :: c_len
        logical,dimension(row_length*rows), intent(inout) :: to_do
        integer,dimension(row_length*rows), intent(inout) :: ind_todo

! local variables
        integer                      :: n,m

        m = 0
!          Compress index for main loops
!CDIR Nodep
           do n=1,c_len
             if(to_do(n))   then
               m=m+1
               to_do(m)    = to_do(n)
               ind_todo(m) = ind_todo(n)
             end if
           end do
           c_len = m

        return

      END SUBROUTINE EXCFNL_CCI
