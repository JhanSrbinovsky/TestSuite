#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      INTEGER FUNCTION IMET(pos1,lhalf)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Returns met grid indicies for half and full
!-                         grids
!-
!-   Outputs : IMET
!-   Inputs  : POS1,lhalf
!-   Inputs  : grid lengths: from module
!-   Outputs :
!-   Controls:
!-
!
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  5.0    31/08/01  Created.  C.E. Johnson
!  6.1    22/01/04  Commented out error message. M.G. Sanderson
!-
!----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!----------------------------------------------------------------------
      REAL,    INTENT(IN) :: pos1         ! position
      LOGICAL, INTENT(IN) :: lhalf        ! T for 1/2 grid points

      CHARACTER(LEN=72)   :: cmessage

      IF (lhalf) THEN                     ! 1/2 gridpoint
        imet = 1 + FLOOR(pos1/dlongm)
      ELSE                                ! Full gridpoint
        imet = FLOOR(pos1/dlongm + 0.5)
      END IF
      imet = imet - lnbound + 1

! Test to see whether point lies within this pe.
!     IF (imet < 1 .OR. imet > nlonpe) THEN
!       cmessage='Position out of bounds'
!       WRITE(6,*) cmessage,'PE: ',mype,' POS1: ',pos1
!       CALL EREPORT('IMET',1,cmessage)
!     END IF

      END FUNCTION IMET
#endif
