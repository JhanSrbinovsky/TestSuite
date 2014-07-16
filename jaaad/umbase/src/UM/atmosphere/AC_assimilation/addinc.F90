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
      SUBROUTINE ADDINC (FIELD,INCR,LEN,NPTS,NROWS,NAVT,ICODE,CMESSAGE)

!  ADD AC INCREMENTS (INCR) TO MODEL PROGNOSTIC VARIABLES (FIELD)
!  WHICH IS A FUNCTION OF ROW
!  DIMENSIONS GIVEN BY LEN,NPTS,NROWS

!  AUTHOR : STUART BELL 17/01/90

      IMPLICIT NONE

      EXTERNAL TIMER

#include "acparm.h"
#include "comacp.h"

      INTEGER NPTS,NROWS,LEN
      REAL                                                              &
     &               FIELD(LEN),                                        &
     &               INCR(LEN)
!*
      INTEGER ICODE
      INTEGER NAVT
      CHARACTER*256 CMESSAGE

      INTEGER IS,IROW,IPT

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('ADDINC  ',3)

!   ADD INCREMENTS TO MODEL FIELD
#if defined(GLOBAL)
!   FOR GLOBAL DOMAIN

       DO 120 IPT=1,LEN
       FIELD(IPT) = FIELD(IPT)+INCR(IPT)
120    CONTINUE

#else
!   USE ONLY THOSE VALUES NOT NEAR EDGE OF LIMITED AREA

!FPP$ NODEPCHK
       DO 130 IROW=1,NROWS
       IS=(IROW-1)*NPTS
        DO 140 IPT=1,NPTS
        FIELD(IS+IPT) = FIELD(IS+IPT)+INCR(IS+IPT)
140     CONTINUE
130    CONTINUE
!
#endif

! ensure RH increments don't generate -ve moisture points
       IF(NAVT == 4)THEN
        DO 160 IPT=1,LEN
        IF(FIELD(IPT) <= 0.0)FIELD(IPT)=0.0
160     CONTINUE
       ENDIF

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('ADDINC  ',4)
      RETURN
      END SUBROUTINE ADDINC
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
#endif
