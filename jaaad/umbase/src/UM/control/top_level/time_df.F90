#if defined(CONTROL) || defined(FLDOP) || defined(FLDC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: TIME_LT --------------------------------------------------
!LL
!LL  Purpose: Function to compare two model times in days/seconds, and
!LL           return a value .TRUE. or .FALSE. depending on whether the
!LL           first time is earlier than the second time.  NOTE: If the
!LL           days are more than 10 years apart, only days are compared.
!LL           Forms a service routine for model date/time and internal
!LL           clock purposes, written for 32-bit portability.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.0
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Model            Modification history from model version 3.2:
!LL version  date
!LL   3.3  05/05/94  Introduced as new deck in association with changes
!LL                  to internal clock for 32-bit portability. TCJ
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: S620
!LL
!LL  Project task: S62
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!

!LL  Routine: TIME_EQ --------------------------------------------------
!LL
!LL  Purpose: Function to compare two model times in days/seconds, and
!LL           return a value .TRUE. or .FALSE. depending on whether the
!LL           first time is equal to the second time.  NOTE: If the days
!LL           are more than 10 years apart, seconds are not checked.
!LL           Forms a service routine for model date/time and internal
!LL           clock purposes, written for 32-bit portability.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.0
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Model            Modification history from model version 3.2:
!LL version  date
!LL   3.3  05/05/94  Introduced as new deck in association with changes
!LL                  to internal clock for 32-bit portability. TCJ
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: S620
!LL
!LL  Project task: S62
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
#endif
#if defined(FLDC) || defined(CONTROL)
!LL  Routine: TIME_DF --------------------------------------------------
!LL
!LL  Purpose: Subroutine to obtain a new model time in days and seconds
!LL           from some reference point, given an increment in days and
!LL           seconds.  Note that the seconds and days increments are
!LL           treated independently so that -ve increments or seconds
!LL           increments larger than the no of seconds in a day are
!LL           handled correctly.
!LL           Forms a service routine for model date/time and internal
!LL           clock purposes, written for 32-bit portability.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.0
!LL
!LL  Author:   R.Rawlins
!LL
!LL  Model            Modification history from model version 3.2:
!LL version  date
!LL   3.3  05/05/94  Introduced as new deck in association with changes
!LL                  to internal clock for 32-bit portability. R.Rawlins
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE TIME_DF(DAYS1,SECS1,DEL_DAYS,DEL_SECS,DAYS2,SECS2)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &     DAYS1,SECS1,                                                 &
                                   ! IN  - days/seconds (input time)
     &     DEL_DAYS,DEL_SECS,                                           &
                                   ! IN  - days/seconds increments
     &     DAYS2,SECS2             ! OUT - days/seconds (output time)
!
      INTEGER                                                           &
     &     SECS_IN_DAY            ! No of seconds in a day
      PARAMETER                                                         &
     &    (SECS_IN_DAY=24*60*60)
!
      DAYS2 = DAYS1 + DEL_DAYS + DEL_SECS/SECS_IN_DAY
      SECS2 = SECS1 + MOD(DEL_SECS,SECS_IN_DAY)
!
      IF(SECS2 <  0) THEN
         SECS2 = SECS2 + SECS_IN_DAY
         DAYS2 = DAYS2 - 1
      ENDIF
!
      IF(SECS2 >= SECS_IN_DAY) THEN
         SECS2 = SECS2 - SECS_IN_DAY
         DAYS2 = DAYS2 + 1
      ENDIF
!
      RETURN
      END SUBROUTINE TIME_DF
#endif
