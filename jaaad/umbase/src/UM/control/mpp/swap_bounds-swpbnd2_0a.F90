#if defined(C96_0A) && ! defined(MAKEBC)

SUBROUTINE swap_bounds(field, row_length, rows, levels,                         &
                       halo_x, halo_y, field_type, l_vector)

! Description: 
!  Dummy stub version of routine SWAP_BOUNDS to be used when 
!  MPP service routines are not selected.
!
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
!
! Author: Paul Burton
! Current code owner: Paul Selwood
!
! Date       Version    Comment
! ------------------------------------------------------------------------------
! 24/04/07     6.5      New dummy stub routine added.                 Ricky Wong
! ------------------------------------------------------------------------------
!

  IMPLICIT NONE


! Arguments:
! ----------
  INTEGER, INTENT(IN) :: row_length ! Number of points on a row (excludes halos)
  INTEGER, INTENT(IN) :: rows       ! Number of rows in a theta field.
                                    ! (excludes halos)
  INTEGER, INTENT(IN) :: levels     ! Number of model levels
  INTEGER, INTENT(IN) :: halo_x     ! Size of halo in "i" direction
  INTEGER, INTENT(IN) :: halo_y     ! Size of halo in "j" direction

  REAL,    INTENT(INOUT)                                                        &
                      :: field(1-halo_x : row_length+halo_x,                    &
                               1-halo_y : rows      +halo_y,                    &
                               levels) 
                                    ! Field to have its halos updated

  INTEGER, INTENT(IN) :: field_type ! Defines the grid interpolation type of the
                                    ! input FIELD (u, v or w)
  LOGICAL, INTENT(IN) :: l_vector   ! TRUE : Data is a horizontal vector
                                    !        component
                                    ! FALSE: Data is a scalar

! Local Parameters:
! -----------------
  INTEGER             :: ErrorStat  
  CHARACTER(Len=11)   :: RoutineName
  CHARACTER(Len=46)   :: Message    


#if ! defined(SCMA) 

!-------------------------------------------------------------------------------
! Code Statements
!-------------------------------------------------------------------------------

  ErrorStat   = 1
  RoutineName = 'SWAP_BOUNDS'
  Message     = 'MPP service routines unavailable - see output.'

  WRITE (6,*) '**ERROR**: SWAP_BOUNDS called but is unavailable.'
  WRITE (6,*) '  Sections C96_1A, C96_1B or C96_1C are required '
  WRITE (6,*) '  (MPP Service routines)                         '

! DEPENDS ON: ereport
  CALL ereport(RoutineName, ErrorStat, Message)

#endif 

  RETURN

!-------------------------------------------------------------------------------

END SUBROUTINE SWAP_BOUNDS
#endif
