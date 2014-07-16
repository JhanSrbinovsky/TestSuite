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
! contains routines: climate_month_date
!
! Purpose: Flux processing routine.
!          Calculates the full date of a climate field which matches
!          the input stash code and month number.
!          Also output the day of the middle of the month.
!
! WARNING: This routine contains mid_month_day_valid hard wired
!          as 15
!----------------------------------------------------------------------
      subroutine climate_month_date( stcode, ValidMonth,                &
#include "aclm1tim.h"
     &       mid_month_day_valid, icode)

      implicit none

! declaration of argument list
      integer stcode     ! IN stash code of field to look for
      integer ValidMonth ! IN month of field to look for
! validity time (in lookup header) of climate field (intent: OUT)
#include "cclm1tim.h"
      integer mid_month_day_valid ! OUT day number of middle of month
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of parameters
#include "clookadd.h"
#include "plookups.h"

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"
#include "clookups.h"

! no local arrays

! declaration of local scalars
      logical ItemFound     ! T => item has been found
      integer fld_no        ! number of lookup table of required field
      integer i             ! do loop index for lookup table number
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'climate_month_date'  ! subroutine name for error messages

! 1. Find a match in the climate lookup tables

      ItemFound = .false.
      do i = 1, Len2_ActualClimate
        if ( LookupClimate(LBMON,i)      ==  ValidMonth  .and.          &
     &       LookupClimate(ITEM_CODE,i)  ==  StCode ) then
          ItemFound = .True.
          fld_no = i
          go to 100
        end if
      end do

100   continue

      if ( .not. ItemFound ) then
        icode = 34
        write(UnWarn,*)CWarn,CSub,                                      &
     &  ' step 1. unable to find climate field with stcode ', stcode,   &
     &  ' for month ', ValidMonth
        go to 9999
      end if


! 2. Set the date in CCLM1TIM from the lookup table
      Clim1Year  = LookupClimate(LBYR, fld_no)
      Clim1Month = LookupClimate(LBMON, fld_no)
      Clim1Day   = LookupClimate(LBDAT, fld_no)
      Clim1Hour  = LookupClimate(LBHR, fld_no)
      Clim1Min   = LookupClimate(LBMIN, fld_no)
      Clim1Sec   = 0

! 3. Calculate the middle day in the month from the lookup table

!     if ( LookupClimate(LBDAT,  fld_no)  == 
!    #     LookupClimate(LBDATD, fld_no)       ) then
!        mid_month_day_valid = LookupClimate(LBDATD, fld_no)

!     else
!        mid_month_day_valid = 0.5 * (LookupClimate(LBDATD,fld_no) + 1)

!     end if

      mid_month_day_valid = 15
      if ( mid_month_day_valid  <   14 .or.                             &
     &     mid_month_day_valid  >   16       ) then
        icode = 35
        write(UnWarn,*)CErr,CSub,                                       &
     &  ' step 3. Lookup table times for climate fields are strange ',  &
     &  ' mid_month_day_valid = ', mid_month_day_valid
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE climate_month_date
!----------------------------------------------------------------------
#endif
