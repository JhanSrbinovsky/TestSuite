#if defined(C84_0A)
!
! Subroutine COPYDIAG_3D -------------------------------------------------------
!
! Desciption: Dummy stub routine of COPYDIAG3D, called if umui stash
!             routines related to section 84 are not selected.
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
!
! Current code owner: Richard Barnes
!
!
! History
! Date       Version   Comment
! ------------------------------------------------------------------------------
! 24/04/07   6.5       New dummy stub added.                          Ricky Wong
! ------------------------------------------------------------------------------
!

SUBROUTINE copydiag_3d(diagout, diagin, row_length, rows, levels, offx_out,     &
                       offy_out, offx_in, offy_in, at_extremity, stlist,        &
                       len_stlist, stash_levels, len_stashlevels, im, is, ie,   &
                       ErrorStat, Message)

  IMPLICIT NONE

  ! Arguments:
  ! ----------
  INTEGER, INTENT(IN) :: row_length      ! Number of points in a row
  INTEGER, INTENT(IN) :: rows            ! Number of rows
  INTEGER, INTENT(IN) :: levels          ! Number of levels in input data
  INTEGER, INTENT(IN) :: offx_out        ! X-Offset dimensions of diagout
  INTEGER, INTENT(IN) :: offy_out        ! Y-Offset dimensions of diagout
  INTEGER, INTENT(IN) :: offx_in         ! X-Offset dimensions of diagin
  INTEGER, INTENT(IN) :: offy_in         ! Y-Offset dimensions of diagin

  REAL,    INTENT(INOUT)                                                        &
                      :: diagout(1-offx_out : row_length+offx_out,              &
                                 1-offy_out : rows      +offy_out,              &
                                 *)
                                         ! Output data

  REAL,    INTENT(IN) :: diagin (1-offx_in  : row_length+offx_in,               &
                                 1-offy_in  : rows      +offy_in,               &
                                 levels)                                                             
                                         ! Input data

  LOGICAL, INTENT(IN) :: at_extremity(4) ! Indicates if this processor is at
                                         ! north, south, east or west of the
                                         ! processor grid

  INTEGER, INTENT(IN) :: len_stlist
  INTEGER, INTENT(IN) :: stlist(len_stlist)              ! Stash list
  INTEGER, INTENT(IN) :: len_stashlevels
  INTEGER, INTENT(IN) :: stash_levels(len_stashlevels,*) ! Stash levels list
  INTEGER, INTENT(IN) :: im                              ! Stash Model number
  INTEGER, INTENT(IN) :: is                              ! Stash Section number
  INTEGER, INTENT(IN) :: ie                              ! Stash Item number


  CHARACTER(len=48), INTENT(INOUT) :: Message
  INTEGER,           INTENT(INOUT) :: ErrorStat          ! Return code:
                                                         !   0 = Normal exit
                                                         ! +ve = Fatal error
                                                         ! -ve = Warning

!-------------------------------------------------------------------------------
! Code Statements
!-------------------------------------------------------------------------------

  Message   = 'STASH service routines unavailable - see output.'
  ErrorStat = 1

  WRITE (6,*) '**ERROR**: COPYDIAG_3D called but is unavailable. '
  WRITE (6,*) '  COPYDIAG_3D should not be called if section     '
  WRITE (6,*) '  C84_1A is not selected.'

  RETURN

!-------------------------------------------------------------------------------
END SUBROUTINE copydiag_3d
#endif
