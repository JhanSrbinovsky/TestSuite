


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Gets the cpu time from the system

! Function Interface:
      REAL FUNCTION Get_CPU_Time()

      IMPLICIT NONE
!
! Description:
!   There is no Fortran standard for calculating CPU time, so by
!   default we just call Get_Wallclock_Time.
!   Alternative "true" CPU times may be generated for other platforms
!   by inserting a call to the relevant routine.
!
! Current Code Owner: A Van der Wal
!
! IE, BC     <- programmer of some or all of previous code or changes
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 5.3       24/04/01  Complete re-write to simplify.  A Van der Wal
! 6.0       06/09/03  Added defs to allow use in makebc/fieldcalc.
!                     R.Sempers
! 6.1       25/05/04  Add proper CPU timing for NEC. P.Selwood
!
! Code Description:
!   Language: FORTRAN 77 plus some Fortran 90 features
!   This code is written to UMDP3 v6 programming standards
!
! System component covered:
! System Task:
!
!- End of header

! Local variables:

! External functions called:

      REAL, EXTERNAL :: Get_Wallclock_Time

! DEPENDS ON: get_wallclock_time
      Get_CPU_Time = Get_Wallclock_Time()

      RETURN
      END FUNCTION Get_CPU_Time

!----------------------------------------------------------------------
!
!+ Gets the elapsed time from the system

! Function Interface:
