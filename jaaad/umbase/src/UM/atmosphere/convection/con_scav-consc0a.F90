#if defined(A05_0A)
!
! SUBROUTINE CONV_SCAV----------------------------------------------------------
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! ******************************COPYRIGHT***************************************
!

SUBROUTINE con_scav(timestep, row_length, rows, levels, cblevel, ctlevel, rain  &
         , snow, aerosol                                                        &
         )

!-------------------------------------------------------------------------------
! Description:
!   Dummy stub of CON_SCAV which is called if Convection Scheme is disabled.
!
!   Current code owner:  R.A. Stratton
!
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.
!
!-------------------------------------------------------------------------------

  IMPLICIT NONE

! Arguments
! ---------
  INTEGER, INTENT(IN) :: row_length                 ! Number of points in row
  INTEGER, INTENT(IN) :: rows                       ! Number of rows
  INTEGER, INTENT(IN) :: levels                     ! Number of levels
  INTEGER, INTENT(IN) :: cblevel (row_length, rows) ! Conv. cloud base level
  INTEGER, INTENT(IN) :: ctlevel (row_length, rows) ! Conv. cloud  top level


  Real, INTENT(IN) :: timestep                  ! Timestep (s).
  Real, INTENT(IN) :: rain   (row_length, rows) ! Rate of rainfall in this layer
                                                ! from above (kg per sq m per s)
  Real, INTENT(IN) :: snow   (row_length, rows) ! Rate of snowfall in this layer
                                                ! from above (kg per sq m per s)
  Real, INTENT(IN) :: aerosol(row_length, rows, levels) ! Aerosol mixing ratio

! Local Variables
! ---------------

  CHARACTER(Len=52)   :: Message
  CHARACTER(Len=8 )   :: RoutineName
  INTEGER             :: ErrorStat        ! Return code:
                                          !   0 = Normal exit
                                          ! +ve = Fatal Error
                                          ! -ve = Warning

!-------------------------------------------------------------------------------
! Code Statements
!-------------------------------------------------------------------------------
        
  ErrorStat   = 1
  RoutineName = 'CON_SCAV'
  Message     = 'Convection Scheme Routines unavailable - see output.'

  WRITE (6,*) '**ERROR**: CON_SCAV called but is unavailable.   '
  WRITE (6,*) '  Sections 5: Convection Scheme is required '


! DEPENDS ON: ereport
  CALL ereport(RoutineName, ErrorStat, Message)

  RETURN

!-------------------------------------------------------------------------------
END SUBROUTINE CON_SCAV
#endif
