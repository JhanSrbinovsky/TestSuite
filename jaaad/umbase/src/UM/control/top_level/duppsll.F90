#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+
! Subroutine Interface:


!- End of subroutine code -------------------------------------------


!+Determine whether a levels list is a duplicate of another levels list
! Subroutine Interface:


!- End of subroutine code ----------------------------------------


!+Determine whether a pseudo lev list is a duplicate of another one
! Subroutine Interface:

      SUBROUTINE DUPPSLL(LDUPLL,NDUPLL)
      IMPLICIT NONE
! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Apr. 95    Original code.  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:

#include "csubmodl.h"
#include "version.h"
#include "cstash.h"
#include "stextend.h"

! Subroutine arguments:

!   Scalar arguments with intent(out):

      LOGICAL LDUPLL
      LOGICAL LLOCAL
      INTEGER NDUPLL

! Local scalars:

      INTEGER I
      INTEGER J

!- End of Header -------------------------------------------------


      LDUPLL=.FALSE.
      NDUPLL=0
      DO 100 I=1,NPSLISTS-1
        IF(LENPLST(I) == LENPLST(NPSLISTS)) THEN
          LLOCAL=.TRUE.
        DO 200 J=1,LENPLST(NPSLISTS)
            IF(PSLIST_D(J,I) /= PSLIST_D(J,NPSLISTS)) THEN
              LLOCAL=.FALSE.
              GOTO 210
            END IF
 200      CONTINUE
 210      CONTINUE
          IF(LLOCAL) THEN
            LDUPLL=.TRUE.
            NDUPLL=I
            RETURN
          END IF
        END IF
 100  CONTINUE


      RETURN
      END SUBROUTINE DUPPSLL

!- End of subroutine code ---------------------------------------
#endif
