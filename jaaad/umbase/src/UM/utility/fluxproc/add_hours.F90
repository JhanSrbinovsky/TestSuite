#if defined(FLUXPROC) || defined(FLXPLPR)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
! History:
! version  date         change
! 4.5      03/09/98     New code
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: add_hours,amend_times
!
! Purpose:  Flux processing routines.
!           add_hours adds an input number of hours
!           to an input date/time (called Clim1)
!           to produce an output date/time (called Valid).
!           amend_times amends date information in Int_Head.
!----------------------------------------------------------------------
      subroutine add_hours(                                             &
#include "aclm1tim.h"
#include "avaltim.h"
     &       IAddHrs )

      implicit none

! declaration of argument list

! input time: intent IN
#include "cclm1tim.h"
! output time: intent OUT
#include "cvaltim.h"

      integer IAddHrs  ! number of hours to add to input time

! no parameters, globals or local arrays used

! declaration of local scalars
      integer CDay  ! century day (of input)
      integer CHour ! century-hour
      integer New_CHour ! new century-hour (after adding input hours)
      integer New_CDay ! new century day

      external date31, date13

!----------------------------------------------------------------------

! 1. Convert the input date to century-day
! DEPENDS ON: date31
      call date31(Clim1Day, Clim1Month, Clim1Year,CDay)

! 2. Calculate century-hour
      CHour = (CDay-1)*24 + Clim1Hour

! 3. Calculate new century-hour
      New_CHour = CHour + IAddHrs

! 4. Calculate new century-day
      New_CDay = 1 + New_Chour/24

! 5. Convert new century-day to new date
! DEPENDS ON: date13
      call date13(New_CDay,ValidDay,ValidMonth,ValidYear)

! 6. Convert rest of values
      ValidHour = New_CHour - 24 * (New_CDay - 1)
      ValidMin = Clim1Min
      ValidSec = Clim1Sec

      return
      END SUBROUTINE add_hours
!----------------------------------------------------------------------
!----------------------------------------------------------------------
#endif
