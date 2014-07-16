
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
      LOGICAL FUNCTION TIME_LT(DAYS1,SECS1,DAYS2,SECS2)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &     DAYS1,SECS1,                                                 &
                                   ! IN  - days/seconds (first time)
     &     DAYS2,SECS2             ! IN  - days/seconds (second time)
!
      IF (DAYS1 >= DAYS2 .AND. SECS1 >= SECS2) THEN
        TIME_LT = .FALSE.
      ELSEIF (DAYS1 <  DAYS2 .AND. SECS1 <  SECS2) THEN
        TIME_LT = .TRUE.
      ELSEIF (ABS(DAYS1-DAYS2) >  3600) THEN
        TIME_LT = DAYS1 <  DAYS2
      ELSE
        TIME_LT = ((DAYS1-DAYS2)*86400 + SECS1-SECS2) <  0
      ENDIF
!
      RETURN
      END FUNCTION TIME_LT

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
