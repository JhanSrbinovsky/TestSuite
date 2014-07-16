#if defined(C84_0A)
!
! Subroutine COPYDIAG ----------------------------------------------------------
!
! Description: Dummy stub routine of COPYDIAG, called if umui stash
!              routines related to section 84 are not selected.
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
!
! Current code owner: Richard Barnes
!
! History
! Date       Version   Comment
! ------------------------------------------------------------------------------
! 24/04/07   6.5       New dummy stub added.                          Ricky Wong
! ------------------------------------------------------------------------------
!

SUBROUTINE copydiag(diagout, diagin, row_length, rows, offx_out, offy_out,      &
                    offx_in, offy_in, at_extremity, im, is, ie, errorstat,      &
                    message)

  IMPLICIT NONE

  ! Arguments:
  ! ----------
  INTEGER, INTENT(IN) :: row_length          ! Number of points in a row
  INTEGER, INTENT(IN) :: rows                ! Number of rows
  INTEGER, INTENT(IN) :: offx_out            ! X-Offset dimensions of diagout
  INTEGER, INTENT(IN) :: offy_out            ! Y-Offset dimensions of diagout               
  INTEGER, INTENT(IN) :: offx_in             ! X-Offset dimensions of diagin
  INTEGER, INTENT(IN) :: offy_in             ! Y-Offset dimensions of diagin

  REAL,    INTENT(INOUT)                                                        &
                      :: diagout(1-offx_out : row_length+offx_out,              &
                                 1-offy_out : rows      +offy_out)              
                                             ! Output field


  REAL,    INTENT(IN) :: diagin (1-offx_in  : row_length+offx_in,               &
                                 1-offy_in  : rows      +offy_in)                                  
                                             ! Input field

  LOGICAL, INTENT(IN) :: at_extremity(4)     ! Indicates if this processor is at
                                             ! north, south, east or west of the
                                             ! processor grid

  INTEGER, INTENT(IN) :: im                  ! Stash Model number
  INTEGER, INTENT(IN) :: is                  ! Stash Section number 
  INTEGER, INTENT(IN) :: ie                  ! Stash Item number

  CHARACTER(Len=48), INTENT(INOUT)   :: Message
  INTEGER,           INTENT(INOUT)   :: ErrorStat  ! Return code:
                                                   !   0 = Normal exit
                                                   ! +ve = Fatal Error
                                                   ! -ve = Warning


!-------------------------------------------------------------------------------
! Code Statements
!-------------------------------------------------------------------------------

  Message   = 'STASH service routines unavailable - see output.'
  ErrorStat = 1

  WRITE (6,*) '**ERROR**: COPYDIAG called but is unavailable.    '
  WRITE (6,*) '  COPYDIAG should not be called if section C84_1A '
  WRITE (6,*) '  is not selected.'

  RETURN

!-------------------------------------------------------------------------------

END SUBROUTINE copydiag
#endif
