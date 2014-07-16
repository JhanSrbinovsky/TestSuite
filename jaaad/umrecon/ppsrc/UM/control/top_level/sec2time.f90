

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: SEC2TIM0 -------------------------------------------------
!LL
!LL  Purpose: Converts from an integer number of elapsed seconds since
!LL           the model basis time to a calendar date/time, using the
!LL           absolute calendar zero point as a reference.  30-day
!LL           month or standard calendar may be used.
!LL           NB: BASIS_TIME_SECS is the number of seconds from the
!LL           calendar zero point to the basis time for the run, and
!LL           is calculated in INITTIME.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL  3.1    12/01/93  Correct error on day 366 in leap year (-CAL360).
!LL   3.3  01/02/94  Modify calling arguments and internal calculations
!LL                  to split days/seconds, for 32-bit portability. TCJ
!LL   3.4  17/06/94  DEF CAL360 replaced by LOGICAL LCAL360
!LL                                                 S.J.Swarbrick
!LL   4.5  10/07/98  Correct possible out of range error in array
!LL                  DAYS_IN_MONTHS for leap years.  JC Thil
!LL   4.5  28/10/98  Changed def line to use superdef UTILS instead of
!LL                  DEF,UTILHIST,OR,DEF,FLDIO,OR,DEF,UTILIO which
!LL                  allowed addition of SCMA def.  K Rogers
!LL   5.1  13/04/00  Small mods to allow negative elapsed times.
!LL                  Adam Clayton
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: S62
!LL
!LL  Project task: S62
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE SEC2TIME(ELAPSED_DAYS,ELAPSED_SECS                     &
     &,                   BASIS_TIME_DAYS,BASIS_TIME_SECS               &
     &,                   I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND &
     &,                   I_DAY_NUMBER,LCAL360)
!
      IMPLICIT NONE
      LOGICAL LCAL360
!
      INTEGER                                                           &
     &     ELAPSED_DAYS,                                                &
                                   ! IN  - elapsed days since basis time
     &     ELAPSED_SECS,                                                &
                                   ! IN  - elapsed secs in part of day
     &     BASIS_TIME_DAYS,                                             &
                                   ! IN  - whole days to basis time
     &     BASIS_TIME_SECS,                                             &
                                   ! IN  - secs in day at basis time
!                                  !       relative to calendar zero
     &     I_SECOND,                                                    &
                                   ! OUT - model time (seconds)
     &     I_MINUTE,                                                    &
                                   ! OUT - model time (minutes)
     &     I_HOUR,                                                      &
                                   ! OUT - model time (hours)
     &     I_DAY,                                                       &
                                   ! OUT - model time (days)
     &     I_MONTH,                                                     &
                                   ! OUT - model time (months)
     &     I_YEAR,                                                      &
                                   ! OUT - model time (years)
     &     I_DAY_NUMBER            ! OUT - model time (day number)

!
!*----------------------------------------------------------------------
!
!  Common blocks
!
! CDATDATA start
!
! Constants needed by routines to calculate day/month/year from
! incremental day number since calendar zero point, and vice-versa,
! when using Gregorian calendar

      INTEGER,PARAMETER:: DAYS_PER_4C = 146097
      INTEGER,PARAMETER:: DAYS_PER_C  = 36524
      INTEGER,PARAMETER:: DAYS_PER_4Y = 1461
      INTEGER,PARAMETER:: DAYS_PER_Y  = 365

      INTEGER,PARAMETER:: DAYS_IN_MONTH(12) =                           &
     &  (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

      INTEGER,PARAMETER:: DAYS_TO_MONTH(12) =                           &
     &  (/0, 31, 59, 90,120,151,181,212,243,273,304,334/)

! CDAYDATA end
!
!  Local variables
!
      INTEGER                                                           &
     &       SECOND      ! number of seconds since calendar zero
!
!L----------------------------------------------------------------------
!L 1. Add elapsed time to basis time in days/seconds to get elapsed
!L     since calendar zero, and convert to hours, minutes, seconds and
!L     total days since calendar zero
!L
      SECOND   = BASIS_TIME_SECS+ELAPSED_SECS
      I_DAY    = BASIS_TIME_DAYS+ELAPSED_DAYS+SECOND/86400
      IF (SECOND  >=  0) THEN
        SECOND = MOD(SECOND,86400)
      ELSE
        SECOND = MOD(SECOND,86400)
       ! Correction for fact that rounding was upwards.
        IF (SECOND  /=  0) THEN
          I_DAY  = I_DAY - 1
          SECOND = SECOND + 86400
        END IF
      END IF
      I_HOUR   = MOD(SECOND/3600,24)
      I_MINUTE = MOD(SECOND/60  ,60)
      I_SECOND = MOD(SECOND,60)
!L----------------------------------------------------------------------
!L 2. Convert day number to date
!L
      IF (LCAL360) THEN
!L
!L 2.1 30-day month (360 day year) calendar
!L
      I_YEAR  = I_DAY/360
      I_MONTH = MOD(I_DAY/30,12)+1
      I_DAY   = MOD(I_DAY,30)+1
      I_DAY_NUMBER = I_DAY+30*(I_MONTH-1)

      ELSE
!L
!L 2.2 Gregorian calendar
!L
      I_YEAR = (I_DAY/DAYS_PER_4C)*400
      I_DAY = I_DAY-(I_DAY/DAYS_PER_4C)*DAYS_PER_4C
!L     Catch special case 31 Dec in leap years
      IF (I_DAY == 4*DAYS_PER_C) THEN
        I_YEAR = I_YEAR+400
        I_DAY = DAYS_PER_Y+1
      ELSE
        I_YEAR = I_YEAR+(I_DAY/DAYS_PER_C)*100
        I_DAY = I_DAY-(I_DAY/DAYS_PER_C)*DAYS_PER_C
        I_YEAR = I_YEAR+(I_DAY/DAYS_PER_4Y)*4
        I_DAY = I_DAY-(I_DAY/DAYS_PER_4Y)*DAYS_PER_4Y
        IF (I_DAY == 4*DAYS_PER_Y) THEN
          I_YEAR = I_YEAR+4
          I_DAY = DAYS_PER_Y+1
        ELSE
          I_YEAR = I_YEAR+(I_DAY/DAYS_PER_Y) + 1
          I_DAY = I_DAY-(I_DAY/DAYS_PER_Y)*DAYS_PER_Y + 1
        ENDIF
      ENDIF
      I_DAY_NUMBER = I_DAY
!L     Find month/day from day no in year
      I_MONTH = 1
      DO WHILE ((I_MONTH  <=  12) .and.                                 &
     &  (I_DAY  >   DAYS_IN_MONTH(I_MONTH)))
        I_DAY = I_DAY-DAYS_IN_MONTH(I_MONTH)
        I_MONTH = I_MONTH+1
      ENDDO
!L     Adjust if leap year and after February
      IF (I_MONTH >  2 .AND. MOD(I_YEAR,4) == 0 .AND.                   &
     &    (MOD(I_YEAR,400) == 0 .OR. MOD(I_YEAR,100) /= 0)) THEN
        I_DAY = I_DAY-1
        IF (I_DAY == 0) THEN
          I_MONTH = I_MONTH-1
          I_DAY = DAYS_IN_MONTH(I_MONTH)
          IF (I_MONTH == 2) I_DAY=29
        ENDIF
      ENDIF
      END IF  !  LCAL360
!L----------------------------------------------------------------------
      RETURN
      END SUBROUTINE SEC2TIME
!
