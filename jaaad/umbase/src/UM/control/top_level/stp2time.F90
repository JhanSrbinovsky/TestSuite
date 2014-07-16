#if defined(CONTROL) || defined(UTILHIST)                              \
 || defined(FLDC) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: STP2TIME -------------------------------------------------
!LL
!LL  Purpose: Converts from an integer number of elapsed steps since
!LL           the model basis time to elapsed whole days and seconds,
!LL           (although the seconds may exceed a day in some cases).
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
!LL   3.3  08/02/94  Introduced as new deck in association with changes
!LL                  to internal clock for 32-bit portability. TCJ
!LL   3.5  07/04/95  Generalise code to allow SECS_IN_PERIOD to be a
!LL                  value corresponding to an exact multiple or exact
!LL                  numerator of seconds in a day. Previous code was
!LL                  only correct for a period of exactly 1 day.
!LL                  R.Rawlins
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
      SUBROUTINE STP2TIME(ELAPSED_STEPS,                                &
     &                    STEPS_IN_PERIOD,SECS_IN_PERIOD,               &
     &                    ELAPSED_DAYS,ELAPSED_SECS)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &     ELAPSED_STEPS,                                               &
                                   ! IN  - elapsed steps since ref time
     &     STEPS_IN_PERIOD,                                             &
                                   ! IN  - steps in period defining T/S
     &     SECS_IN_PERIOD,                                              &
                                   ! IN  - secs  in period defining T/S
     &     ELAPSED_DAYS,                                                &
                                   ! OUT - elapsed days since ref time
     &     ELAPSED_SECS            ! OUT - elapsed secs in part of day
!                                  !       or days since ref time
! Local Parameters
      INTEGER   SECS_IN_DAY        ! no. of seconds in 1 day
      PARAMETER(SECS_IN_DAY=3600*24)
! Local Scalars
      INTEGER   FACTOR             ! ratio of integers is exact factor
!L----------------------------------------------------------------------
!L 1. Perform integer arithmetic to compute elapsed days/seconds from
!L    elapsed steps and timestep definition.  The subroutine assumes
!L    that SECS_IN_PERIOD is a whole number multiple of a day, or a
!L    whole number divisor of a day, but does not check explicitly.
!L
      IF(SECS_IN_PERIOD >= SECS_IN_DAY) THEN
         FACTOR       = SECS_IN_PERIOD/SECS_IN_DAY ! no. days in period
         ELAPSED_DAYS = (ELAPSED_STEPS*FACTOR)/STEPS_IN_PERIOD
      ELSE
         FACTOR       = SECS_IN_DAY/SECS_IN_PERIOD
         ELAPSED_DAYS = (ELAPSED_STEPS/FACTOR)/STEPS_IN_PERIOD
      ENDIF
!
! Find number of steps per day
      FACTOR=(SECS_IN_DAY*STEPS_IN_PERIOD)/SECS_IN_PERIOD
      IF(FACTOR >  1) THEN    ! Timestep is less than 1 day
         ELAPSED_SECS = ( ELAPSED_STEPS -                               &
     &                   (ELAPSED_STEPS/FACTOR)*FACTOR  )*              &
     &                   (SECS_IN_PERIOD/STEPS_IN_PERIOD)
      ELSE                    ! Timestep is multiple of 1 day
         ELAPSED_SECS = 0
      ENDIF
!
      RETURN
      END SUBROUTINE STP2TIME
#endif
