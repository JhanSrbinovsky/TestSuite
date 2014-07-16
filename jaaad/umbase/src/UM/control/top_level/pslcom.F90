#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Construct preliminary STASH list of user requests
!
! Subroutine Interface:


!- End of Subroutine code -------------------------------------------


!+ Compress out unused pseudo levels lists
      SUBROUTINE PSLCOM(NRECS)

!  Description:
!
!  Method:
!
!  Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
!  Global variables:

#include "csubmodl.h"
#include "stparam.h"
#include "version.h"
#include "cstash.h"
#include "stextend.h"

! Subroutine arguments:

!   Scalar arguments with intent(in):

      INTEGER NRECS

! Local Scalars:

      INTEGER ICOUNT

! Local arrays:

      INTEGER IPOS (NPSLISTP)  ! POSITION IN OLD LIST OF THE NEW
      INTEGER IPOS1(NPSLISTP)  ! POSITION IN THE NEW LIST OF THE OLD

!- End of Header ---------------------------------------------------


      ICOUNT=0

      DO I=1,NPSLISTS         ! LOOP DOWN THE LISTS
        IF(LENPLST(I) /= 0) THEN   ! USED LIST
          ICOUNT=ICOUNT+1
          IPOS(ICOUNT)=I
          IPOS1(I)=ICOUNT
        ELSE
          IPOS1(I)=0
        END IF
      END DO

      NPSLISTS=ICOUNT

      DO    I=1,NPSLISTS
        LENPLST(I)=LENPLST(IPOS(I))
        DO  J=1,NPSLEVP
          PSLIST_D(J,I)=PSLIST_D(J,IPOS(I))
        END DO
      END DO

      DO   I=NPSLISTS+1,NPSLISTP
        LENPLST(I)=0
        DO J=1,NPSLEVP
          PSLIST_D(J,I)=0
        END DO
      END DO

      DO I=1,NRECS
        IF(LIST_S(st_pseudo_out,I) /= 0) THEN
          LIST_S(st_pseudo_out,I)=IPOS1(LIST_S(st_pseudo_out,I))
        END IF
      END DO

      RETURN
      END SUBROUTINE PSLCOM

!- End of Subroutine code ------------------------------------------
#endif
