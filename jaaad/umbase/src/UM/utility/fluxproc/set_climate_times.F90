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
! contains routines: set_climate_times
!
! Purpose: Flux processing routine.
!          Preliminaries for interpolating climate fields in time
!          sets date/times required to extract two climate fields
!          and calculates the weights to give to them
!
! WARNING: does not test ICODE properly yet
!----------------------------------------------------------------------
      subroutine set_climate_times ( stcode,                            &
#include "avaltim.h"
#include "aclm1tim.h"
#include "aclm2tim.h"
     &     weight1, weight2, icode  )

      implicit none

! declaration of argument list

      integer stcode ! IN stash code of field being accessed
! validity time of field (intent: IN)
#include "cvaltim.h"
! validity time (in lookup headers) of climate fields (intent: OUT)
#include "cclm1tim.h"
#include "cclm2tim.h"
      real weight1   ! OUT  weight to give to climate field 1
      real weight2   ! OUT  weight to give to climate field 2
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! no parameters

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"

! no local arrays

! declaration of local scalars
! mid_month_day_# is day number of middle of month determined from
! lookup tables of climate file
      integer mid_month_day_valid  ! for validity time
      integer mid_month_day_clim1  ! for climate field 1
      integer mid_month_day_clim2  ! for climate field 2
      integer Year1     ! year for climate field 1
      integer Year2     ! year for climate field 2
      integer CDay      ! century day
      integer C1Hour    ! century hour of climate field 1
      integer C2Hour    ! century hour of climate field 2
      integer ValHour   ! century hour of validity time

      external climate_month_date, date31
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'set_climate_times'  ! subroutine name for error messages

! 1. calculate the mid-month day of the validity time

! DEPENDS ON: climate_month_date
      call climate_month_date( stcode, ValidMonth,                      &
#include "aclm1tim.h"
     &       mid_month_day_valid, icode )


! 2. determine the months of the first and second climate fields to use

      if (  ValidDay  >   mid_month_day_valid ) then
        Clim1Month = ValidMonth
      else
        Clim1Month = ValidMonth - 1
      end if

      if ( Clim1Month  ==  0) then
        Clim1Month = 12
      end if

      Clim2Month = Clim1Month + 1

      if ( Clim2Month  ==  13) then
        Clim2Month = 1
      end if

! 3. find mid-month days of the first and second climate months
!    and the full dates in the lookup tables for these fields (one
!    of the main outputsfrom this routine)

! DEPENDS ON: climate_month_date
      call climate_month_date( stcode, Clim1Month,                      &
#include "aclm1tim.h"
     &       mid_month_day_clim1,icode )

! DEPENDS ON: climate_month_date
      call climate_month_date( stcode, Clim2Month,                      &
#include "aclm2tim.h"
     &       mid_month_day_clim2,icode )

! 4. find the weights to give two months when interpolating to
!    validity time

! 4.1 find the years for months 1 and 2

      if ( Clim1Month  ==  12 .and. ValidMonth  ==  1 ) then
        Year1 = ValidYear - 1
      else
        Year1 = ValidYear
      end if

      if ( Clim2Month  ==  1 .and. ValidMonth  ==  12 ) then
        Year2 = ValidYear + 1
      else
        Year2 = ValidYear
      end if

! 4.2 find the relative times (in hours) of the three dates

! DEPENDS ON: date31
      call date31(mid_month_day_clim1, Clim1Month, Year1,CDay)
      C1Hour = (CDay-1)*24

! DEPENDS ON: date31
      call date31(mid_month_day_clim2, Clim2Month, Year2,CDay)
      C2Hour = (CDay-1)*24

! DEPENDS ON: date31
      call date31(ValidDay, ValidMonth, ValidYear,CDay)
      ValHour = (CDay-1)*24 + ValidHour

! 4.3 calculate the weights
      weight1 = real( C2Hour - ValHour ) / real( C2Hour - C1Hour )
      weight2 = 1.0 - weight1

9999  continue
      return
      END SUBROUTINE set_climate_times
!----------------------------------------------------------------------
#endif
