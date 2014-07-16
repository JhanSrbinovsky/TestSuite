#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: SET_HISTORY_VALUES----------------------------------------
!LL
!LL  Purpose: Updates certain history file values before any calls
!LL           to write out history file
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Author:   R A Stratton
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.3  08/02/94  Modify calls to TIME2SEC/SEC2TIME to output/input
!LL                  elapsed times in days & secs, for portability. TCJ
!LL   3.4  17/06/94  Argument LCAL360 added and passed to SEC2TIME,
!LL                   TIME2SEC.                        S.J.Swarbrick
!LL   3.4  27/10/94  Set year correctly for Gregorian calendar
!LL  4.1  18/04/96  Set RUN_ID from EXPT_ID & JOB_ID.  RTHBarnes.
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: H
!LL
!LL  Project task: H
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE SET_HISTORY_VALUES
!
      IMPLICIT NONE
!
!*----------------------------------------------------------------------
!  Common blocks
!
#include "cmaxsize.h"
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
#include "ctime.h"
#include "ccontrol.h"
!
!  Subroutines called
!
      EXTERNAL SEC2TIME,TIME2SEC,STP2TIME
!
!  Local variables
!
      INTEGER                                                           &
     &     ELAPSED_DAYS                                                 &
                                   ! Elapsed time since basis time(days)
     &,    ELAPSED_SECS                                                 &
                                   ! Elapsed time since basis time(secs)
     &,    ELAPSED_DAYA                                                 &
                                   ! time since basis time atmosphere
     &,    ELAPSED_SECA                                                 &
                                   ! time since basis time atmosphere
     &,    ELAPSED_DAYO                                                 &
                                   ! time since basis time Ocean
     &,    ELAPSED_SECO                                                 &
                                   ! time since basis time Ocean
     &,    FINAL_DAYS                                                   &
                                   ! target time (days)
     &,    FINAL_SECS                                                   &
                                   ! target time (seconds)
     &,    J_YEAR                                                       &
                                   ! Year
     &,    J_MONTH                                                      &
                                   ! month
     &,    J_DAY                                                        &
                                   ! day
     &,    J_HOUR                                                       &
                                   ! hour
     &,    J_MINUTE                                                     &
                                   ! minute
     &,    J_SECOND                                                     &
                                   ! second
     &,    J_DAY_NUMBER            ! day number
!L
!L----------------------------------------------------------------------
!L 1. Set Run target  reached so far ready for possible resubmission
!L
#if defined(ATMOS)
#if !defined(OCEAN)
! DEPENDS ON: stp2time
      CALL STP2TIME(STEPim(A_IM),                                       &
     &              STEPS_PER_PERIODim(A_IM),SECS_PER_PERIODim(A_IM),   &
     &              ELAPSED_DAYS,ELAPSED_SECS)
#else
! DEPENDS ON: stp2time
      CALL STP2TIME(STEPim(A_IM),                                       &
     &              STEPS_PER_PERIODim(A_IM),SECS_PER_PERIODim(A_IM),   &
     &              ELAPSED_DAYA,ELAPSED_SECA)
! DEPENDS ON: stp2time
      CALL STP2TIME(STEPim(O_IM),                                       &
     &              STEPS_PER_PERIODim(O_IM),SECS_PER_PERIODim(O_IM),   &
     &              ELAPSED_DAYO,ELAPSED_SECO)
      IF (ELAPSED_DAYA <  ELAPSED_DAYO) THEN
        ELAPSED_DAYS=ELAPSED_DAYA
        ELAPSED_SECS=ELAPSED_SECA
      ELSE
        ELAPSED_DAYS=ELAPSED_DAYO
        ELAPSED_SECS=ELAPSED_SECO
      ENDIF
#endif
#endif
#if defined(OCEAN)
! DEPENDS ON: stp2time
      CALL STP2TIME(STEPim(O_IM),                                       &
     &              STEPS_PER_PERIODim(O_IM),SECS_PER_PERIODim(O_IM),   &
     &               ELAPSED_DAYS,ELAPSED_SECS)
#endif
!L
!L  convert ELAPSED_DAYS/SECS back into a time relative to Basis time
!L
! DEPENDS ON: sec2time
      CALL SEC2TIME(ELAPSED_DAYS,ELAPSED_SECS,0,0,                      &
     &              J_YEAR,J_MONTH,J_DAY,J_HOUR,                        &
     &              J_MINUTE,J_SECOND,J_DAY_NUMBER,LCAL360)
      IF (LCAL360) THEN
      RUN_RESUBMIT_TARGET(1)=J_YEAR
      ELSE
        RUN_RESUBMIT_TARGET(1)=J_YEAR-1
      ENDIF
      RUN_RESUBMIT_TARGET(2)=J_MONTH-1
      RUN_RESUBMIT_TARGET(3)=J_DAY-1
      RUN_RESUBMIT_TARGET(4)=J_HOUR
      RUN_RESUBMIT_TARGET(5)=J_MINUTE
      RUN_RESUBMIT_TARGET(6)=J_SECOND

!L
!L check to see if final target date reached if resubmission is
!L  operating.
!L
      IF (RUN_RESUBMIT == 'Y') THEN
!L
!L  Work out final target date required  in seconds
!L
        J_YEAR   = RUN_TARGET_END(1)
        J_MONTH  = RUN_TARGET_END(2)
        J_DAY    = RUN_TARGET_END(3)
        J_HOUR   = RUN_TARGET_END(4)
        J_MINUTE = RUN_TARGET_END(5)
        J_SECOND = RUN_TARGET_END(6)
      IF (LCAL360) THEN
! DEPENDS ON: time2sec
        CALL TIME2SEC(J_YEAR,J_MONTH+1,J_DAY+1,J_HOUR,J_MINUTE,         &
     &                J_SECOND,0,0,FINAL_DAYS,FINAL_SECS,LCAL360)
      ELSE
! DEPENDS ON: time2sec
        CALL TIME2SEC(J_YEAR+1,J_MONTH+1,J_DAY+1,J_HOUR,J_MINUTE,       &
     &                J_SECOND,0,0,FINAL_DAYS,FINAL_SECS,LCAL360)
      END IF
        IF (ELAPSED_DAYS >= FINAL_DAYS.AND.ELAPSED_SECS >= FINAL_SECS)  &
     &  THEN
          RUN_RESUBMIT='N'
      WRITE(6,*)'SETHIST: final target date reached automatic resubmiss &
     &ion switched off'
        ENDIF
      ENDIF
!  Set RUN_IN for qxhistreport
      RUN_ID = EXPT_ID//JOB_ID

      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE SET_HISTORY_VALUES
!
#endif
