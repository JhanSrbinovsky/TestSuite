#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C)
#if defined(MPP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Parallel UM: Fills global halos with sensible numbers
!
!  Vn.    Date    Comment
!  6.0  17/09/03  Add def for new NEC opt section C96_1C. R Barnes
! Subroutine interface:
      SUBROUTINE SET_SIDES()

!     This routine temporarily deleted while compile errors
!     are sorted out

      WRITE(6,*) 'Error: SET_SIDES called : IMPLEMENT IT !!!'
!      CALL ABORT()
       RETURN
       END SUBROUTINE SET_SIDES

#endif
#endif
