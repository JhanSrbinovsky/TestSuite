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

      SUBROUTINE DUPLEVL(NLEVELS,LDUPLL,NDUPLL)
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

!   Scalar arguments with intent(in):

      INTEGER NLEVELS

!   Scalar arguments with intent(out):

      LOGICAL LDUPLL
      LOGICAL LLOCAL
      INTEGER NDUPLL

! Local scalars:

      INTEGER I
      INTEGER J

!- End of Header -----------------------------------------------------

      LDUPLL=.FALSE.
      NDUPLL=0
      DO 100 I=1,NLEVELS-1
        IF((LEVLST_S(1,I) == LEVLST_S(1,NLEVELS)).AND.                  &
     &  (LLISTTY(I) == LLISTTY(NLEVELS))) THEN
          LLOCAL=.TRUE.
          DO 200 J=2,LEVLST_S(1,NLEVELS)+1
            IF(LLISTTY(NLEVELS) == 'I') THEN
              IF(LEVLST_S(J,I) /= LEVLST_S(J,NLEVELS)) THEN
                LLOCAL=.FALSE.
                GOTO 210
              END IF
            ELSE
              IF(RLEVLST_S(J,I) /= RLEVLST_S(J,NLEVELS)) THEN
                LLOCAL=.FALSE.
                GOTO 210
              END IF
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
      END SUBROUTINE DUPLEVL

!- End of subroutine code ----------------------------------------


!+Determine whether a pseudo lev list is a duplicate of another one
! Subroutine Interface:


!- End of subroutine code ---------------------------------------
#endif
