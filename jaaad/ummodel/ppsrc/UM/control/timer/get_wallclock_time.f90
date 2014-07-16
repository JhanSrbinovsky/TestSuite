


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Gets the cpu time from the system

! Function Interface:

!----------------------------------------------------------------------
!
!+ Gets the elapsed time from the system

! Function Interface:
      REAL FUNCTION Get_Wallclock_Time()

      IMPLICIT NONE
!
! Description:
!   The system function SYSTEM_CLOCK is used to return the numbers of
!   seconds which have elapsed.
!
! Current Code Owner: Anette Van der Wal
!
! IE, PB, Richard Barnes  <- programmer of some or all of previous
!                            code or changes
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 5.3       24/04/01  Complete re-write to simplify.  A Van der Wal
!
! Code Description:
!   Language: FORTRAN 77 plus some Fortran 90 features
!   This code is written to UMDP3 v6 programming standards
!
! System component covered:
! System Task:
!
!- End of header

! Local variables
      INTEGER, SAVE :: Start_Count=-1
      INTEGER, SAVE :: Old_Count=0
      REAL, SAVE    :: Rollover=0.0
      INTEGER       :: Count, Count_Rate, Count_Max, Elapsed_Count
      REAL, SAVE    :: OneOver_Count_Rate=0.0

! Intrinsic procedures called:
      INTRINSIC SYSTEM_CLOCK

      CALL SYSTEM_CLOCK(Count=Count,Count_Rate=Count_Rate,              &
     &                  Count_Max=Count_Max)

      IF ((Old_Count  <   Start_Count) .AND.                            &
     &    ((Count  <   Old_Count) .OR. (Count  >   Start_Count))) THEN
        Rollover=Rollover+(REAL(Count_Max)/REAL(Count_Rate))
      ENDIF

      IF (Start_Count  ==  -1) THEN
        Start_Count=Count
        OneOver_Count_Rate=1.0/REAL(Count_Rate)
      ENDIF

      Elapsed_Count=Count-Start_Count
      IF (Elapsed_Count  <   0) Elapsed_Count=Elapsed_Count+Count_Max

      Get_Wallclock_Time = Rollover+                                    &
     &                    (REAL(Elapsed_Count)*OneOver_Count_Rate)
      Old_Count=Count

      RETURN
      END FUNCTION Get_Wallclock_Time
