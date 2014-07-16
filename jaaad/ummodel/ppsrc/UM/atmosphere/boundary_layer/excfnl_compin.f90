
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE EXCFNL_COMPIN------------------------------------------
!!!
!!!  Purpose: Compute condition for inner interation and
!!!           compress index array
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  6.2    Jan 2006  New deck as part of optimisation of EXCF_NL
!!!                   by J-C Rioual.                 A P Lock
!!!  6.4    Nov 2006  Correct j1 indexing to avoid 
!!!                   out-of-bounds errors       Adrian Lock
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
      SUBROUTINE EXCFNL_COMPIN (                                        &
     & row_length, rows, c_len_i, ind_todo_i,                           &
     & todo_inner, UP, WB_RATIO, dec_thres, switch                      &
     & )
!
      IMPLICIT NONE

! Intent IN:

      integer, intent(in)           :: row_length, rows

      integer, intent(in), dimension(row_length*rows) :: UP
      real,    intent(in), dimension(row_length, rows)::                &
     & WB_RATIO                                                         &
                                    ! WBN_INT/WBP_INT
     &,DEC_THRES                    ! Local decoupling threshold

      integer, intent(in)           :: switch ! =1 for KSURF, 2 for KTOP

! Intent INOUT:

      integer, intent(inout)        :: c_len_i
      integer, intent(inout), dimension(row_length*rows) :: ind_todo_i
      logical, intent(inout), dimension(row_length*rows) :: todo_inner

! local variables
      integer                       :: n,m
      integer                       :: l, i1, j1

!     Check for active elements

        select case (switch)
          case(1)
!           ! For top of KSURF
            do n=1,c_len_i
              l = ind_todo_i(n)
              j1=(l-1)/row_length+1
              i1=l-(j1-1)*row_length

              todo_inner(n) = (                                         &
     &       (UP(l) == 1 .AND. WB_RATIO(i1,j1) <  DEC_THRES(i1,j1)) .OR.&
!               ! keep working up while wb_ratio lt thres
     &       (UP(l) == 0 .AND. WB_RATIO(i1,j1) >  DEC_THRES(i1,j1)) )
!               ! keep working down while wb_ratio gt thres

            end do

          case(2)
!           ! For base of KTOP
            do n=1,c_len_i
              l = ind_todo_i(n)
              j1=(l-1)/row_length+1
              i1=l-(j1-1)*row_length

              todo_inner(n) = (                                         &
     &       (UP(l) == 1 .AND. WB_RATIO(i1,j1) >  DEC_THRES(i1,j1)) .OR.&
!               ! keep working up while wb_ratio gt thres
     &       (UP(l) == 0 .AND. WB_RATIO(i1,j1) <  DEC_THRES(i1,j1)) )
!               ! keep working down while wb_ratio lt thres

            end do

          end select
!
!     Compress
!
        m = 0
!CDIR Nodep
        do n=1,c_len_i
          if(todo_inner(n))  then
            m = m+1
            todo_inner(m) = todo_inner(n)
            ind_todo_i(m) = ind_todo_i(n)
          end if
        end do
        c_len_i = m

        return

      END SUBROUTINE EXCFNL_COMPIN
