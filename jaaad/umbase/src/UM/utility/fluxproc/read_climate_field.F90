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
! contains routines: read_climate_field
!
! Purpose: Flux processing routine.
!          Finds and interpolates in time a climate field specified by
!          user's search criteria and returns it and its lookup
!          table by the argument list
!----------------------------------------------------------------------
      subroutine read_climate_field(StCode, IVTOffHr,                   &
     &           ldebug, Int_Head, Real_Head,                           &
     &           ncols, nrows, field,                                   &
#include "argppx.h"
     &           icode)

      implicit none

! declaration of parameters
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"
#include "plookups.h"
#include "c_mdi.h"

! declaration of argument list

! user's search criteria
      integer StCode       ! IN
      integer IVTOffHr     ! IN offset from validity time in hours

! debug control variable
      logical ldebug          ! IN T => output debugging info
      logical l_climate_field ! Set to true if reading climate field

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
#include "cclm1tim.h"
#include "cclm2tim.h"

! declaration of local arrays
      real field1(ncols,nrows)  ! values from earlier climate field
      real field2(ncols,nrows)  ! values from later climate field

! declaration of local scalars
      real weight1   ! weight to give to 1st climate field
      real weight2   ! weight to give to 2nd climate field

! declaration of externals
      external add_hours, read_one_field, set_climate_times,            &
     &         interp_time
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_climate_field' ! subroutine name for error messages
      l_climate_field = .true.

      if (LClimate) then

! 1. calculate validity time of NWP data required

! DEPENDS ON: add_hours
        call add_hours(                                                 &
#include "areftim.h"
#include "avaltim.h"
     &       IVTOffHr)


! 2. set up times of fields to look for and
!    time interpolation coefficients

! DEPENDS ON: set_climate_times
        call set_climate_times ( StCode,                                &
#include "avaltim.h"
#include "aclm1tim.h"
#include "aclm2tim.h"
     &       weight1, weight2, icode )

        if ( icode  >   0) then
          write(UnWarn,*)CWarn//CSub//                                  &
     &    '2. failed setting climate times',                            &
     &    ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
          go to 9999
        end if

! 3. Extract climate field before validity time

! DEPENDS ON: read_one_field
        call read_one_field (UnitClimate, ITEM_CODE, Stcode,            &
#include "aclm1tim.h"
     &         Len_FixHd, FixHdClimate, Len1_Lookup,                    &
     &         Len2_ActualClimate, LookupClimate, LookFldNoClimate,     &
     &         ldebug, l_climate_field,                                 &
     &         Len_IntHd, Len_RealHd, Int_Head, Real_Head,              &
     &         ncols, nrows, field1,                                    &
#include "argppx.h"
     &         icode)


        if ( icode  >   0) then
          write(UnWarn,*)CWarn//CSub//                                  &
     &    ' 3. failed reading 1st climate field',                       &
     &    ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
          go to 9999
        end if

! 4. Extract climate field after validity time
! DEPENDS ON: read_one_field
        call read_one_field (UnitClimate, ITEM_CODE, Stcode,            &
#include "aclm2tim.h"
     &         Len_FixHd, FixHdClimate, Len1_Lookup,                    &
     &         Len2_ActualClimate, LookupClimate, LookFldNoClimate,     &
     &         ldebug, l_climate_field,                                 &
     &         Len_IntHd, Len_RealHd, Int_Head, Real_Head,              &
     &         ncols, nrows, field2,                                    &
#include "argppx.h"
     &         icode)

        if ( icode  >   0) then
          write(UnWarn,*)CWarn//CSub//                                  &
     &    '4. failed reading 2nd climate field',                        &
     &    'for stash code ', stcode, '; IVTOffHr = ', IVTOffHr
          go to 9999
        end if

! 5. If found: interpolate in time to validity time

! DEPENDS ON: interp_time
        call interp_time(Int_Head, ncols, nrows, rmdi,                  &
#include "avaltim.h"
     &         weight1, weight2, Field1, Field2, Field)

! 6.    Output standard message and exit routine

        write(UnStd,*)CStd//CSub//'climate field stcode ',              &
     &  stcode, '; IVTOffHr = ', IVTOffHr, ' extracted'

! 7.  Write times to integer headers
! DEPENDS ON: amend_times
        call amend_times (                                              &
#include "avaltim.h"
     &                   Int_Head,Len_IntHd )
         go to 9999


! 7. Else If there is no climate file return an error code

      else !  LClimate

        icode = 7
        write(UnWarn,*)CWarn//CSub//'7. Climate file is not open,',     &
     &    ' so no climate data can be extracted.'

      end if ! LClimate

9999  continue
      return
      END SUBROUTINE read_climate_field
!----------------------------------------------------------------------
#endif
