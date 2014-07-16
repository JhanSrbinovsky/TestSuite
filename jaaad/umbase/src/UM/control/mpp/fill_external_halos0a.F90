#if  defined(C96_0A) && ! defined(UTILIO)
!
! Subroutine FILL_EXTERNAL_HALOS------------------------------------------------
!
! Description: 
!  Dummy stub version of routine FILL_EXTERNAL_HALOS to be used when 
!  MPP service routines are not selected.
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
!
!
! Author: Paul Burton
! Current code owner: Paul Selwood
!
! Date       Version   Comment
! ------------------------------------------------------------------------------
! 24/04/07   6.5       New dummy stub added.                          Ricky Wong
! ------------------------------------------------------------------------------
!

SUBROUTINE fill_external_halos(field, row_length, rows, levels, halo_x, halo_y)

  IMPLICIT NONE


! Arguments:
! ----------
  INTEGER, INTENT(IN) :: row_length ! Number of points on a row
                                    ! (excluding halos)
  INTEGER, INTENT(IN) :: rows       ! Number of rows in a theta field
                                    ! (not including halos)
  INTEGER, INTENT(IN) :: levels     ! Number of model levels
  INTEGER, INTENT(IN) :: halo_x     ! Size of halo in "i" direction
  INTEGER, INTENT(IN) :: halo_y     ! Size of halo in "j" direction


  REAL,    INTENT(INOUT) :: field(1-halo_x : row_length+halo_x,                 &
                                  1-halo_y : rows      +halo_y,                 &
                                  levels)
                                    ! Field to have its halos updated

! Local Parameters:
! ------------------

  INTEGER            :: ErrorStat
  CHARACTER (LEN=19) :: RoutineName
  CHARACTER (LEN=46) :: Message


#if ! defined(SCMA)

!-------------------------------------------------------------------------------
! Code Statements
!-------------------------------------------------------------------------------

  ErrorStat   =  1
  RoutineName = 'FILL_EXTERNAL_HALOS'
  Message    =  'MPP service routines unavailable - see output.'

  WRITE (6,*) '**ERROR**: FILL_EXTERNAL_HALOS called but is '
  WRITE (6,*) '  unavailable. Sections C96_1A, C96_1B or    '
  WRITE (6,*) '  C96_1C are required (MPP Service routines).'

! DEPENDS ON: ereport
  CALL ereport(RoutineName, ErrorStat, Message)

#endif

  RETURN

!-------------------------------------------------------------------------------

END SUBROUTINE fill_external_halos
#endif
