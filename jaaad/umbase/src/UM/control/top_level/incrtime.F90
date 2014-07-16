#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: INCRTIME -------------------------------------------------
!LL
!LL  Purpose: Increments the model time by one atmosphere or ocean
!LL           timestep according to which is currently being
!LL           integrated.  Also updates timestamps in dump LOOKUP
!LL           headers of PROGNOSTIC fields (diagnostic LOOKUP headers
!LL           are updated exclusively by STASH).
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE INCRTIME(                                              &
#include "argduma.h"
     &       internal_model,ICODE,CMESSAGE )
!
      IMPLICIT NONE
!
!*L Arguments
!L
#include "parparm.h"
#include "typsize.h"
#include "typduma.h"
#include "nstypes.h"
!
      INTEGER       internal_model    ! IN : internal model identifier
      INTEGER       ICODE             ! OUT: Error return code
      CHARACTER*80  CMESSAGE          ! OUT: Error return message
!
!*----------------------------------------------------------------------
!  Common blocks
!
#include "cmaxsize.h"
#include "csubmodl.h"
#include "chsunits.h"
#include "clookadd.h"
#include "ccontrol.h"
#include "chistory.h"
#include "ctime.h"
#include "cruntimc.h"
#include "ctfilt.h"
#include "lbc_coup.h"
! Include printstatus
#include "cprintst.h"
!
!  Local variables
!
      INTEGER                                                           &
     &     ELAPSED_DAYS,                                                &
                                   ! Elapsed days  since basis time
     &     ELAPSED_SECS,                                                &
                                   ! Elapsed secs  since basis time
     &     ELAPSED_DAYS_PREV,                                           &
                                   ! Elapsed days, end of previous step
     &     ELAPSED_SECS_PREV,                                           &
                                   ! Elapsed secs, end of previous step
     &     I                       ! Loop index
      INTEGER                                                           &
                                   ! Local scalars of internal model
     &  STEP                                                            &
                                   !  arrays.
     &, STEPS_PER_PERIOD                                                &
     &, SECS_PER_PERIOD
     
! UM6.5   MODEL_ANALYSIS_HRS replaced by MODEL_ANALYSIS_MINS in cntlall.h -
!         Requires ELAPSED_HRS changed to REAL 
      REAL ELAPSED_HRS             ! Elapsed hours since basis time

!L
!L----------------------------------------------------------------------
!L 1. General timestep, increment STEP by one and update atmos
!L    elapsed seconds (integer) relative to basis time
!L
        STEP            =STEPim(internal_model)
        STEPS_PER_PERIOD=STEPS_PER_PERIODim(internal_model)
        SECS_PER_PERIOD =SECS_PER_PERIODim(internal_model)

! DEPENDS ON: stp2time
        CALL STP2TIME(STEP,STEPS_PER_PERIOD,SECS_PER_PERIOD,            &
     &                ELAPSED_DAYS_PREV,ELAPSED_SECS_PREV)
        STEP = STEP+1
! DEPENDS ON: stp2time
        CALL STP2TIME(STEP,STEPS_PER_PERIOD,SECS_PER_PERIOD,            &
     &                ELAPSED_DAYS,ELAPSED_SECS)

        STEPim(internal_model) = STEP
        H_STEPim(internal_model) = STEP

#if defined(ATMOS)
!L
!L 1.1 If integrating backwards, negate elapsed times.
!L
      IF (L_Backwards .AND. internal_model  ==  ATMOS_SM) THEN
        ELAPSED_DAYS_PREV = -ELAPSED_DAYS_PREV
        ELAPSED_SECS_PREV = -ELAPSED_SECS_PREV
        ELAPSED_DAYS      = -ELAPSED_DAYS
        ELAPSED_SECS      = -ELAPSED_SECS
      END IF
#endif

!L
!L 1.2 Set FORECAST_HRS - the number of hours relative to the current
!L                        data time.
!L
      MODEL_ANALYSIS_HRS =  REAL(MODEL_ANALYSIS_MINS)/60.0
      ELAPSED_HRS = REAL(ELAPSED_SECS)/3600 + ELAPSED_DAYS*24

      IF ((RUN_ASSIM_MODE /= "None      " .OR. L_IAU) .AND.             &
     &     ELAPSED_HRS >= MODEL_ANALYSIS_HRS) THEN
        ! At or beyond data time reset step:
        FORECAST_HRS = ELAPSED_HRS - MODEL_ANALYSIS_HRS
      ELSE
        ! Data time of initial dump still valid:
        FORECAST_HRS = ELAPSED_HRS - DATA_MINUS_BASIS_HRS
      END IF

      IF ( PrintStatus  >  PrStatus_Oper ) THEN
        write(6,*) 
        write(6,*) 'incrtime: ELAPSED SECS       ', ELAPSED_SECS
        write(6,*) 'incrtime: ELAPSED HRS        ', ELAPSED_HRS
        write(6,*) 'incrtime: FORECAST HRS       ', FORECAST_HRS
        write(6,*) 'incrtime: MODEL_ANALYSIS_HRS ', MODEL_ANALYSIS_HRS
        write(6,*)
      END IF

!L
!L 1.3 Set LBC_FC_HRS - the number of hours relative to the last data
!L                      time that will apply during this run.
!L
      IF (L_LBC_COUP) THEN
        IF (RUN_ASSIM_MODE /= "None      " .OR. L_IAU) THEN
          ! Data time reset during run:
          LBC_FC_HRS = ELAPSED_HRS - MODEL_ANALYSIS_HRS
        ELSE
          ! No data time reset:
          LBC_FC_HRS = FORECAST_HRS
        END IF
      END IF

!L----------------------------------------------------------------------
!L 2. Convert elapsed seconds since basis time to calendar time/date
!L
! DEPENDS ON: sec2time
      CALL SEC2TIME(ELAPSED_DAYS_PREV,ELAPSED_SECS_PREV,                &
     &              BASIS_TIME_DAYS,BASIS_TIME_SECS,                    &
     &              PREVIOUS_TIME(1),PREVIOUS_TIME(2),PREVIOUS_TIME(3), &
     &              PREVIOUS_TIME(4),PREVIOUS_TIME(5),PREVIOUS_TIME(6), &
     &              PREVIOUS_TIME(7),LCAL360)
! DEPENDS ON: sec2time
      CALL SEC2TIME(ELAPSED_DAYS,ELAPSED_SECS,                          &
     &              BASIS_TIME_DAYS,BASIS_TIME_SECS,                    &
     &              I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,      &
     &              I_DAY_NUMBER,LCAL360)
!L----------------------------------------------------------------------
!L 3. Copy date/time information into the dump header and update
!L    VALIDITY TIME and FORECAST PERIOD in prognostic field LOOKUP
!L    headers.
!L
#if defined(ATMOS)
      IF (internal_model == atmos_im) THEN
        A_FIXHD(28) = I_YEAR
        A_FIXHD(29) = I_MONTH
        A_FIXHD(30) = I_DAY
        A_FIXHD(31) = I_HOUR
        A_FIXHD(32) = I_MINUTE
        A_FIXHD(33) = I_SECOND
        A_FIXHD(34) = I_DAY_NUMBER
        DO I=1,A_PROG_LOOKUP
          A_LOOKUP(LBYR  ,I)=I_YEAR
          A_LOOKUP(LBMON ,I)=I_MONTH
          A_LOOKUP(LBDAT ,I)=I_DAY
          A_LOOKUP(LBHR  ,I)=I_HOUR
          A_LOOKUP(LBMIN ,I)=I_MINUTE
          A_LOOKUP(LBDAY ,I)=I_DAY_NUMBER
          A_LOOKUP(LBFT,I)=FORECAST_HRS
        ENDDO
!L Update run date time
! These appear not to be really used.
        ACTUAL_ENDT(1,A_IM) = I_YEAR
        ACTUAL_ENDT(2,A_IM) = I_MONTH
        ACTUAL_ENDT(3,A_IM) = I_DAY
        ACTUAL_ENDT(4,A_IM) = I_HOUR
        ACTUAL_ENDT(5,A_IM) = I_MINUTE
        ACTUAL_ENDT(6,A_IM) = I_SECOND
!       LENGTH(A_IM) = STEP(A_IM)
!       A_LENGTH=STEP
      ENDIF
#endif
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE INCRTIME
!
#endif
