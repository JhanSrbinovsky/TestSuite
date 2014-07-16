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
! 6.1      30/07/04     Changes required for global and prelim merge.
!                       A. Hines.
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: read_accum_flds
!
! Purpose:  Flux processing routine.
!           Reads fields for validity time and validity time minus
!           six hours. These two fields are then manipulated to
!           obtain an accumulation for the six hour period.
!
! Uses:     StCode and to read NWP files;
!           xstcode to read climate fields
!----------------------------------------------------------------------
      subroutine read_accum_flds(StCode, IVTOffHr,                      &
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
#include "cvaloff.h"
      real time            ! timescale in seconds for division
!(SAS)**** parameter ( time = 21600 )

! declaration of argument list

! search criteria

!       Uses    StCode to read NWP files
!               stcode to read climate fields
      integer StCode       ! IN StCode value to test

!       Reference date is used with IVTOffHr to define validity
!       time needed
      integer IVTOffHr     ! IN offset from validity time in hours

! debug control variable
      logical ldebug          ! IN T => output debugging info
      logical l_climate_field ! Set to false initially

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
#include "c_mdi.h"
#include "clookups.h"

#include "creftim.h"
#include "cvaltim.h"


! declaration of local arrays
      real fieldVT(ncols,nrows)  ! field at validity time
      real fieldM6(ncols,nrows)  ! field at validity time minus 6 hours
      real fieldint(ncols,nrows) ! intermediate field for calculation

! declaration of local scalars
      real timediv        ! division scale for field (6x3600)
      integer IM6OffHr    ! secondary validity time offset for VT-6

! declaration of logicals
      logical l_preferred_VT   ! OUT test for preferred field at VT
      logical l_preferred_M6   ! OUT test for preferred field at VT-6
      logical l_previous_VT    ! OUT test for previous field at VT
      logical l_previous_M6    ! OUT test for previous field at VT-6

      character *256 cmessage  ! error message



! declaration of externals
      external add_hours, read_one_field, read_climate_field,           &
     &  FieldSub,ScalarMult,check_header


!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_accum_flds'  ! subroutine name for error messages
      l_preferred_M6 = .false.
      l_preferred_VT = .false.
      l_previous_M6 = .false.
      l_previous_VT = .false.
      l_climate_field = .false.

      time = ValidityPeriod * 3600

! 1. calculate validity time minus 6 hours of NWP data required
      IM6OffHr = IVTOffHr - ValidityPeriod
! DEPENDS ON: add_hours
      call add_hours(                                                   &
#include "areftim.h"
#include "avaltim.h"
     &       IM6OffHr)

!----------------------------------------------------------------------
! 2. Check headers for preferred and previous to see if they exist
!----------------------------------------------------------------------
      if ( LPreferred ) then
! DEPENDS ON: check_header
        call check_header (StCode,Len1_Lookup,                          &
     &                         Len2_ActualPreferred,                    &
     &                         LookupPreferred,                         &
#include "avaltim.h"
     &                         l_preferred_M6)
      endif
      if ( LPrevious ) then
! DEPENDS ON: check_header
        call check_header (StCode,Len1_Lookup,                          &
     &                         Len2_ActualPrevious,                     &
     &                         LookupPrevious,                          &
#include "avaltim.h"
     &                         l_previous_M6)
      endif

! 2.1 Calculate Validity Time and check if VT exists
! DEPENDS ON: add_hours
      call add_hours(                                                   &
#include "areftim.h"
#include "avaltim.h"
     &       IVTOffHr)
      if ( LPreferred) then
! DEPENDS ON: check_header
        call check_header (StCode,Len1_Lookup,                          &
     &                         Len2_ActualPreferred,                    &
     &                         LookupPreferred,                         &
#include "avaltim.h"
     &                         l_preferred_VT)
      endif
      if ( LPrevious ) then
! DEPENDS ON: check_header
        call check_header (StCode,Len1_Lookup,                          &
     &                         Len2_ActualPrevious,                     &
     &                         LookupPrevious,                          &
#include "avaltim.h"
     &                         l_previous_VT)
      endif

!----------------------------------------------------------------------
! 3. Read preferred VT&VT-6 if they exist else previous if they do
!----------------------------------------------------------------------
      if ( l_preferred_M6 .and. l_preferred_VT ) then
! DEPENDS ON: read_one_field
        call read_one_field (UnitPreferred, ITEM_CODE, StCode,          &
#include "avaltim.h"
     &       Len_FixHd, FixHdPreferred,Len1_Lookup,                     &
     &       Len2_ActualPreferred, LookupPreferred, LookFldNoPreferred, &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldVT,                                     &
#include "argppx.h"
     &       icode)

        if ( icode  <=  0) then
! 3.1 if successful, issue standard message and exit routine
          write(UnStd,*)CStd//CSub//                                    &
     &     'NWP preferred field (VT) StCode ',                          &
     &     StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
        else
! 3.2 else write warning message and reset icode
          write(UnWarn,*)CWarn//CSub//                                  &
     &     'NWP preferred field (VT) StCode ',                          &
     &     StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
          l_preferred_VT = .false.
        end if
        icode = 0     ! reset icode

! 3.3 If preferred VT has been read, then read preferred VT-6
        if ( l_preferred_VT ) then
! DEPENDS ON: add_hours
          call add_hours(                                               &
#include "areftim.h"
#include "avaltim.h"
     &          IM6OffHr)
! DEPENDS ON: read_one_field
          call read_one_field (UnitPreferred, ITEM_CODE, StCode,        &
#include "avaltim.h"
     &       Len_FixHd, FixHdPreferred,Len1_Lookup,                     &
     &       Len2_ActualPreferred, LookupPreferred, LookFldNoPreferred, &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldM6,                                     &
#include "argppx.h"
     &       icode)

          if ( icode  <=  0) then
! 3.4 if successful, issue standard message and exit routine
            write(UnStd,*)CStd//CSub//                                  &
     &       'NWP preferred field (VT-6) StCode ',                      &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
          else
! 3.5 else write warning message and reset icode
            write(UnWarn,*)CWarn//CSub//                                &
     &       'NWP preferred field (VT-6) StCode ',                      &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
            l_preferred_M6 = .false.
          end if
          icode = 0     ! reset icode
        endif    ! l_preferred_VT
      endif    ! l_preferred_VT / l_preferred_M6

! 3.6 If either preferred VT or preferred M6 has not been read
!     read previous VT and VT-6
      if ( (.not. l_preferred_M6 .or. .not. l_preferred_VT) .and.       &
     &     (      l_previous_M6 .and. l_previous_VT       ) .and.       &
     &            LPrevious ) then
! DEPENDS ON: add_hours
        call add_hours(                                                 &
#include "areftim.h"
#include "avaltim.h"
     &        IVTOffHr)
! DEPENDS ON: read_one_field
        call read_one_field (UnitPrevious, ITEM_CODE, StCode,           &
#include "avaltim.h"
     &       Len_FixHd, FixHdPrevious,Len1_Lookup,                      &
     &       Len2_ActualPrevious, LookupPrevious, LookFldNoPrevious,    &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldVT,                                     &
#include "argppx.h"
     &       icode)

        if ( icode  <=  0) then
! 3.7 if successful, issue standard message and exit routine
          write(UnStd,*)CStd//CSub//                                    &
     &     'NWP previous field (VT) StCode ',                           &
     &     StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
        else
! 3.8 else write warning message and reset icode
          write(UnWarn,*)CWarn//CSub//                                  &
     &     'NWP previous field (VT) StCode ',                           &
     &    StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
          l_previous_VT = .false.
        end if
        icode = 0     ! reset icode
        if ( l_previous_VT) then
! DEPENDS ON: add_hours
          call add_hours(                                               &
#include "areftim.h"
#include "avaltim.h"
     &          IM6OffHr)
! DEPENDS ON: read_one_field
          call read_one_field (UnitPrevious, ITEM_CODE, StCode,         &
#include "avaltim.h"
     &       Len_FixHd, FixHdPrevious,Len1_Lookup,                      &
     &       Len2_ActualPrevious, LookupPrevious, LookFldNoPrevious,    &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldM6,                                     &
#include "argppx.h"
     &       icode)

          if ( icode  <=  0) then
! 3.9 if successful, issue standard message and exit routine
            write(UnStd,*)CStd//CSub//                                  &
     &       'NWP previous field (VT-6) StCode ',                       &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
          else
! 3.10 else write warning message and reset icode
            write(UnWarn,*)CWarn//CSub//                                &
     &       'NWP previous field (VT-6) StCode ',                       &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
            l_previous_M6 = .false.
          end if
          icode = 0     ! reset icode
          endif   ! l_previous_VT
      endif

! Special case for first validity time in file: if accumulations to
! T+0 are not available, instead then first accumulation is from
! T+0, so only one field required.
!
      if ( IVTOffHr  ==  IValidOffHr(1) ) then
!
! 3.11 If no M6 fields have been read, but preferred VT is available
!      read preferred VT
        if ( (.not. l_preferred_M6 .and. .not. l_previous_M6) .and.     &
     &       l_preferred_VT ) then
! DEPENDS ON: add_hours
          call add_hours(                                               &
#include "areftim.h"
#include "avaltim.h"
     &        IVTOffHr)
! DEPENDS ON: read_one_field
          call read_one_field (UnitPreferred, ITEM_CODE, StCode,        &
#include "avaltim.h"
     &       Len_FixHd, FixHdPreferred,Len1_Lookup,                     &
     &       Len2_ActualPreferred, LookupPreferred, LookFldNoPreferred, &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldVT,                                     &
#include "argppx.h"
     &       icode)

          if ( icode  <=  0) then

! 3.12 if successful, issue standard message and exit routine
            write(UnStd,*)CStd//CSub//                                  &
     &       'NWP preferred field (VT) StCode ',                        &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
          else

! 3.13 else write warning message and reset icode and validity time
!      for read of previous file
            write(UnWarn,*)CWarn//CSub//                                &
     &       'NWP preferred field (VT) StCode ',                        &
     &      StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
            l_preferred_VT = .false.
          end if
          icode = 0     ! reset icode

        endif

! 3.14 If no M6 fields have been read, and preferred VT is not
!      available, read previous VT if available
        if ( (.not. l_preferred_M6 .and. .not. l_previous_M6 .and.      &
     &        .not. l_preferred_VT ) .and. l_previous_VT ) then
! DEPENDS ON: add_hours
          call add_hours(                                               &
#include "areftim.h"
#include "avaltim.h"
     &        IVTOffHr)
! DEPENDS ON: read_one_field
          call read_one_field (UnitPrevious, ITEM_CODE, StCode,         &
#include "avaltim.h"
     &       Len_FixHd, FixHdPrevious,Len1_Lookup,                      &
     &       Len2_ActualPrevious, LookupPrevious, LookFldNoPrevious,    &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldVT,                                     &
#include "argppx.h"
     &       icode)

          if ( icode  <=  0) then

! 3.15 if successful, issue standard message and exit routine
            write(UnStd,*)CStd//CSub//                                  &
     &       'NWP previous field (VT) StCode ',                         &
     &       StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
          else

! 3.16 else write warning message and reset icode
            write(UnWarn,*)CWarn//CSub//                                &
     &       'NWP previous field (VT) StCode ',                         &
     &      StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
            l_previous_VT = .false.
          end if
          icode = 0     ! reset icode

        endif

      endif  ! First validity time in file
!----------------------------------------------------------------------
! 4. If there is a preferred field for VT-6 and for VT
!    or previous for VT-6 and VT then do accumulation
!----------------------------------------------------------------------
      if ( (l_preferred_VT .and. l_preferred_M6) .or.                   &
     &     (l_previous_VT .and. l_previous_M6) ) then
! DEPENDS ON: fieldsub
         call FieldSub (ncols,nrows,rmdi,                               &
     &                   fieldVT,fieldM6,                               &
     &                   fieldint,                                      &
     &                   icode,cmessage)

! 4.1 Now divide the result by a scalar using ScalarMult
        timediv = 1.0 / time
! DEPENDS ON: scalarmult
        call ScalarMult (ncols,nrows,rmdi,timediv,                      &
     &                      fieldint,field,                             &
     &                      icode,cmessage)

! 4.2 Write times to integer header
! DEPENDS ON: add_hours
        call add_hours(                                                 &
#include "areftim.h"
#include "avaltim.h"
     &       IVTOffHr)
! DEPENDS ON: amend_times
        call amend_times (                                              &
#include "avaltim.h"
     &                   Int_Head,Len_IntHd )
        goto 9999

      elseif ( l_preferred_VT .OR. l_previous_VT ) then

        timediv = 1.0 / time
! DEPENDS ON: scalarmult
        call ScalarMult (ncols,nrows,rmdi,timediv,                      &
     &                      fieldVT,field,                              &
     &                      icode,cmessage)

! DEPENDS ON: add_hours
        call add_hours(                                                 &
#include "areftim.h"
#include "avaltim.h"
     &       IVTOffHr)

! DEPENDS ON: amend_times
        call amend_times (                                              &
#include "avaltim.h"
     &                   Int_Head,Len_IntHd )

        goto 9999

      endif   ! test for both fields

!----------------------------------------------------------------------
! 5. Otherwise extract field from climate file if available
!----------------------------------------------------------------------
      if (LClimate) then
! DEPENDS ON: read_climate_field
        call read_climate_field(StCode, IVTOffHr,                       &
     &           ldebug, Int_Head, Real_Head,                           &
     &           ncols, nrows, field,                                   &
#include "argppx.h"
     &           icode)

         if ( icode  <=  0) then
            write(UnStd,*)CStd//CSub//'5. climate field extracted  ',   &
     &       ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
            go to 9999
         else

            write(UnWarn,*)CWarn//CSub//                                &
     &       '5. failed to retrieve climate field ',                    &
     &       ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
            icode = 0
         end if   ! icode
      end if !  LClimate

!----------------------------------------------------------------------
! 6. If no data has been successfully extracted return an error code
!----------------------------------------------------------------------
      icode = 5
      write(UnErr,*)CErr//CSub//'6. failed to extract any data',        &
     &    ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr

9999  continue
      return
      END SUBROUTINE read_accum_flds
!----------------------------------------------------------------------
#endif
