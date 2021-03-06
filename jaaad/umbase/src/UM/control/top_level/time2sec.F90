#if defined(CONTROL) || defined(RECON) || defined(UTILS)               \
 || defined(SCMA) || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: TIME2SEC -------------------------------------------------
!LL
!LL  Purpose: Converts from calendar date/time to an integer number
!LL           of elapsed seconds since the model basis time, using the
!LL           absolute calendar zero point as a reference.  30-day
!LL           month or standard calendar may be used.
!LL           NB: BASIS_TIME_SECS is the number of seconds from the
!LL           calendar zero point to the basis time for the model, and
!LL           is calculated in INITTIME.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.3  01/02/94  Modify calling arguments and internal calculations
!LL                  to split days/seconds, for 32-bit portability. TCJ
!LL   3.4  20/06/94  DEF CAL360 replaced by LOGICAL LCAL360
!LL                                                  S.J.Swarbrick
!LL   4.5  28/10/98  Changed def line to use superdef UTILS instead of
!LL                  DEF,UTILHIST,OR,DEF,FLDIO,OR,DEF,UTILIO which
!LL                  allowed addition of SCMA def.  K Rogers
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
      SUBROUTINE TIME2SEC (I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND&
     &,                    BASIS_TIME_DAYS,BASIS_TIME_SECS              &
     &,                    ELAPSED_DAYS,ELAPSED_SECS,LCAL360)
!
      IMPLICIT NONE
      LOGICAL LCAL360
!
      INTEGER                                                           &
     &     I_YEAR,                                                      &
                                   ! IN  - model time (years)
     &     I_MONTH,                                                     &
                                   ! IN  - model time (months)
     &     I_DAY,                                                       &
                                   ! IN  - model time (days)
     &     I_HOUR,                                                      &
                                   ! IN  - model time (hours)
     &     I_MINUTE,                                                    &
                                   ! IN  - model time (minutes)
     &     I_SECOND,                                                    &
                                   ! IN  - model time (seconds)
     &     BASIS_TIME_DAYS,                                             &
                                   ! IN  - whole days to basis time
     &     BASIS_TIME_SECS,                                             &
                                   ! IN  - secs in day at basis time
     &     ELAPSED_DAYS,                                                &
                                   ! OUT - elapsed days since basis time
     &     ELAPSED_SECS            ! OUT - elapsed secs in part of day
!                                  !       relative to basis time
!
!*----------------------------------------------------------------------
!
!  Common blocks
!
#include "cdaydata.h"
!
!  Local variables
!
      INTEGER                                                           &
     &       YEAR                                                       &
                         ! years
     &,      DAY         ! number of days since calendar zero
!
!L----------------------------------------------------------------------
!L 1. Add up days from time zero to specified time
!L
      IF (LCAL360) THEN
!L
!L 1.1 30-day month (360 day year) calendar
!L
      DAY = 360*I_YEAR+30*(I_MONTH-1)+I_DAY-1
!L
      ELSE
!L
!L 1.2 Gregorian calendar
!L
!L     If leap year and after 28 February, adjust day number by one
      IF (MOD(I_YEAR,4) == 0   .AND.                                    &
     &   (MOD(I_YEAR,400) == 0 .OR. MOD(I_YEAR,100) /= 0) .AND.         &
     &    I_MONTH >  2) THEN
        DAY = I_DAY+1
      ELSE
        DAY = I_DAY
      ENDIF
!L     Add on days in the preceding months in the year
      DAY = DAY + DAYS_TO_MONTH(I_MONTH) - 1
      YEAR = I_YEAR - 1
!L     Add on days up to the specified year
      DAY =  DAY+(YEAR/400)*DAYS_PER_4C
      YEAR = YEAR-(YEAR/400)*400
      DAY =  DAY+(YEAR/100)*DAYS_PER_C
      YEAR = YEAR-(YEAR/100)*100
      DAY =  DAY+(YEAR/4)*DAYS_PER_4Y
      YEAR = YEAR-(YEAR/4)*4
      DAY =  DAY+YEAR*DAYS_PER_Y
!L
      END IF       ! LCAL360
!L----------------------------------------------------------------------
!L 2. Convert days, hours and minutes to days/secs since calendar zero,
!L     and subtract basis time in days/secs to get elapsed time since
!L     basis, converted to whole days and +ve no of secs in partial day
!L
      ELAPSED_DAYS=DAY-BASIS_TIME_DAYS
      ELAPSED_SECS=3600*I_HOUR+60*I_MINUTE+I_SECOND-BASIS_TIME_SECS
      IF (ELAPSED_SECS >= 86400) THEN
        ELAPSED_DAYS=ELAPSED_DAYS+ELAPSED_SECS/86400
        ELAPSED_SECS=MOD(ELAPSED_SECS,86400)
      ELSEIF (ELAPSED_SECS <  0) THEN
        ELAPSED_DAYS=ELAPSED_DAYS+(ELAPSED_SECS-86399)/86400
        ELAPSED_SECS=MOD(ELAPSED_SECS,86400)+86400
      ENDIF
!
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE TIME2SEC
#endif
