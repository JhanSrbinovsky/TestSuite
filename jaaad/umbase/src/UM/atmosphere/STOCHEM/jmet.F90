#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      INTEGER FUNCTION JMET(pos2,lhalf)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Returns met grid indicies for half and full
!-                         grids
!-
!-   Outputs : JMET
!-   Inputs  : POS2,lhalf
!-   Inputs  : LOBOUND,grid length: from module
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
      LOGICAL, INTENT(IN) :: lhalf       ! T for half grid points
      REAL,    INTENT(IN) :: pos2        ! position

!     CHARACTER(LEN=72)   :: cmessage

      IF (lhalf) THEN                    ! 1/2 gridpoint
        jmet = FLOOR(pos2/dlatm)+1-lobound+1
      ELSE                               ! Full gridpoint
        jmet = FLOOR(pos2/dlatm+0.5)-lobound+1
      END IF
! Test to see whether point lies within this pe.
!     IF (jmet < 1 .OR. jmet > nlatpe) THEN
!       cmessage = 'Position out of bounds'
!       WRITE(6,*) cmessage,'PE: ',mype,' POS2: ',pos2
!       CALL EREPORT('JMET',1,cmessage)
!     END IF

      END FUNCTION JMET
#endif
