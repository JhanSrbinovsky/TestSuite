#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINES ADDINC  RELAXC ---------------------------------------
!LL
!LL  Purpose : Add analysis increments to model field.
!LL
!LL  For use on Cray Y-MP
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  S.Bell     <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  6.1  17/08/04: Use PE edge points in limited are case
!LL                 Adam Maycock (J. Bornemann lodged).
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!L*  ARGUMENTS
!LL  SUBROUTINE RELAXC -------------------------------------------------
!LL
!LL  For use on Cray
!LL
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  Version 1.0 : Written 20/4/90 by Dave Robinson.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LL  Purpose : Scale analysis increments by relaxation coefficients.
!LL
!LLEND------------------------------------------------------------------
      SUBROUTINE RELAXC (INCR,LEN,NPTS,NROWS,RELAX,ICODE,CMESSAGE)

!     SCALE AC INCREMENTS (INCR) BY RELAXATION COEFFICIENT (RELAX)
!     WHICH IS A FUNCTION OF ROW
!     DIMENSIONS GIVEN BY LEN,NPTS,NROWS

      IMPLICIT NONE

      EXTERNAL TIMER

#include "acparm.h"
#include "comacp.h"

      INTEGER NROWS,NPTS,LEN
      REAL                                                              &
     &               INCR(LEN),                                         &
     &               RELAX(NROWS)
!*
      INTEGER ICODE
      CHARACTER*256 CMESSAGE

      INTEGER IS,JROW,JPT

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RELAXC  ',3)

!     MULTIPLY BY RELAXATION COEFFICIENT
      IS=1
      DO 100 JROW=1,NROWS
       DO 110 JPT=0,NPTS-1
       INCR(IS+JPT) = INCR(IS+JPT)*RELAX(JROW)
110    CONTINUE
      IS=IS+NPTS
100   CONTINUE

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER ('RELAXC  ',4)
      RETURN
      END SUBROUTINE RELAXC
#endif
