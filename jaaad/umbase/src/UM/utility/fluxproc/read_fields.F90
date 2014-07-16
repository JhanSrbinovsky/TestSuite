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
! contains routines: read_fields
!
! Purpose: Flux processing routine.
!          Finds a field according to user's search criteria and
!          returns it and its lookup table by the argument list
!
! Uses StCode to read NWP and climate files
!----------------------------------------------------------------------
      subroutine read_fields (StCode, IVTOffHr,                         &
     &               ldebug, Int_Head, Real_Head,                       &
     &               ncols, nrows, field,                               &
#include "argppx.h"
     &               icode)

      implicit none

! declaration of parameters
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"
#include "plookups.h"

! declaration of argument list

      integer StCode       ! IN   stash code of fields to search for

!       Reference date is used with IVTOffHr to define validity
!       time needed
      integer IVTOffHr     ! IN offset from validity time in hours

! debug control variable
      logical ldebug       ! IN T => output debugging info

! lookup tables
      integer Int_Head(Len_IntHd) ! OUT
      real Real_Head(Len_RealHd)  ! OUT

! output field
      integer ncols             ! IN  number of columns
      integer nrows             ! IN  number of rows
      real field(ncols,nrows)   ! OUT field values

! error code
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected


! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"

#include "clookups.h"

#include "creftim.h"
#include "cvaltim.h"

! declaration of local logical

      logical l_climate_field ! Set to false initially

! no local arrays

! declaration of externals
      external add_hours, read_one_field, read_climate_field

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_fields'  ! subroutine name for error messages
      l_climate_field = .false.
! 1. calculate validity time of NWP data required

! DEPENDS ON: add_hours
      call add_hours(                                                   &
#include "areftim.h"
#include "avaltim.h"
     &       IVTOffHr)

!----------------------------------------------------------------------
! 2. Extract field from preferred file if available
!----------------------------------------------------------------------
      if ( LPreferred ) then
! 2.1 try to read preferred field
! DEPENDS ON: read_one_field
        call read_one_field (UnitPreferred, ITEM_CODE, Stcode,          &
#include "avaltim.h"
     &       Len_FixHd, FixHdPreferred,Len1_Lookup,                     &
     &       Len2_ActualPreferred, LookupPreferred, LookFldNoPreferred, &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, field,                                       &
#include "argppx.h"
     &       icode)

        if ( icode  <=  0) then

! 2.2 if successful, issue standard message and exit routine
          write(UnStd,*)CStd//CSub//'NWP preferred field stash code ',  &
     &    StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
! 2.2.1  Write times to integer header
! DEPENDS ON: amend_times
          call amend_times (                                            &
#include "avaltim.h"
     &                   Int_Head,Len_IntHd )
          go to 9999
        else

! 2.3 else write warning message and reset icode
         write(UnWarn,*)CWarn//CSub//'NWP preferred field stash code ', &
     &    StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
        end if
        icode = 0     ! reset icode

      end if    ! LPreferred

!----------------------------------------------------------------------
! 3. Otherwise extract field from preferred file if available
!----------------------------------------------------------------------
      if ( LPrevious ) then

! 3.1 try to read previous field
! DEPENDS ON: read_one_field
       call read_one_field (UnitPrevious, ITEM_CODE, Stcode,            &
#include "avaltim.h"
     &       Len_FixHd, FixHdPrevious,Len1_Lookup,                      &
     &       Len2_ActualPrevious, LookupPrevious, LookFldNoPrevious,    &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, field,                                       &
#include "argppx.h"
     &       icode)


        if ( icode  <=  0) then

! 3.2 if successful, issue standard message and exit routine
          write(UnStd,*)CStd//CSub//'NWP previous field stash code ',   &
     &    Stcode, '; IVTOffHr = ', IVTOffHr, ' extracted'
! 3.2.1  Write times to integer header
! DEPENDS ON: amend_times
          call amend_times (                                            &
#include "avaltim.h"
     &                   Int_Head,Len_IntHd )
          go to 9999
        else

! 3.3 else write warning message and reset icode
          write(UnWarn,*)CWarn//CSub//'NWP previous field stash code ', &
     &    StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
        end if
        icode = 0     ! reset icode

      end if    ! LPrevious

!----------------------------------------------------------------------
! 4. Otherwise extract field from climate file if available
!----------------------------------------------------------------------
      if ( LClimate ) then
! DEPENDS ON: read_climate_field
        call read_climate_field(StCode, IVTOffHr,                       &
     &           ldebug, Int_Head, Real_Head,                           &
     &           ncols, nrows, field,                                   &
#include "argppx.h"
     &           icode)

        if ( icode  <=  0) then
          write(UnStd,*)CStd//CSub//'4. climate field extracted  ',     &
     &     ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
          go to 9999
        else

          write(UnWarn,*)CWarn//CSub//                                  &
     &     '4. failed to retrieve climate field ',                      &
     &     ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
          icode = 0
        endif     ! icode
      endif     ! LClimate

!----------------------------------------------------------------------
! 5. If no data has been successfully extracted return an error code
!----------------------------------------------------------------------
      icode = 5
      write(UnErr,*)CErr//CSub//'5. failed to extract any data',        &
     &    ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr

9999  continue
      return
      END SUBROUTINE read_fields
!----------------------------------------------------------------------
#endif
