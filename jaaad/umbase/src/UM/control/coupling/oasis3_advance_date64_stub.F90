#if !defined(OASIS3)
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
SUBROUTINE OASIS3_ADVANCE_DATE64(ierr64)

  ! Description: This subroutine is a stub routine for
  !              OASIS3_ADVANCE_DATE64. Run time logic
  !              should prevent it from actually being called.
  !
  ! Author: R. Hill
  ! Current Code Owner : R. Hill
  !
  !=============================================================

  IMPLICIT NONE

  INTEGER :: ierr64

  CHARACTER(Len=52)   :: Message
  CHARACTER(Len=20)   :: RoutineName
  INTEGER             :: ErrorStat        ! Return code:
  !   0 = Normal exit
  ! +ve = Fatal Error
  ! -ve = Warning

  ErrorStat   = 1
  RoutineName = 'oasis3_advance_date64_stub'
  Message     = 'OASIS3 Routines unavailable - see output.'

  WRITE (6,*) '**ERROR**: oasis3_advance_date unavailable.'
  WRITE (6,*) 'Check OASIS3 cpp key is set'

  ! DEPENDS ON: ereport
  CALL ereport(RoutineName, ErrorStat, Message)

  RETURN
END SUBROUTINE OASIS3_ADVANCE_DATE64
#endif
