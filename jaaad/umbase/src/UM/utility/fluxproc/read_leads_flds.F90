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
! contains routines: read_leads_flds
!
! Purpose: Flux processing routine.
!          Finds a field according to user's search criteria.
!          Each element in the ice fraction field is tested to see
!          if it is greater than a specified constant. If it is
!          not then climatology is used for the required field.
!          If l_leads is true, the constant is (1 - minleadsfrac),
!          else the constant is minicefrac.
!          Returns field and its lookup table by the argument list.
!
!  Uses:   StCode and StCAICE  to read NWP files;
!          stcode to read climate fields
!
!          WARNING: If StCode = 3231 (sublimation rate), the input
!                   NWP field must be divided by 1200 to get sublimation
!                   rate in kg/m^2/s. This is hard-wired.
!----------------------------------------------------------------------
      subroutine read_leads_flds(StCode,StCAICE,                        &
     &                    IVTOffHr, ldebug, l_leads,Int_Head,           &
     &                    Real_Head, ncols, nrows, field,               &
#include "argppx.h"
     &                    icode)

      implicit none

! declaration of parameters
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"
#include "plookups.h"


! declaration of argument list

! search criteria

!       Uses    StCode to read NWP files
!               stcode to read climate fields

      integer StCode     ! IN stash code value to test
      integer StCAICE    ! IN icefrac code to test

!       Reference date is used with IVTOffHr to define validity
!       time needed
      integer IVTOffHr     ! IN offset from validity time in hours

! declare logicals
! debug control variable
      logical ldebug       ! IN T => output debugging info
      logical l_leads      ! if T => then using minleadsfrac
                           ! if F => then using minicefrac

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

! declaration of logicals
      logical l_present_fieldNWP   ! test for NWP field
      logical l_present_icefrac    ! test for Icefrac field
      logical l_climate_field      ! set to false initially

! declaration of local arrays
      real fieldNWP (ncols,nrows)    ! nwp field
      real fieldClim (ncols,nrows)   ! Climate field
      real icefrac (ncols,nrows)     ! icefrac field
      real time                      ! division factor for sublimation
      parameter (time = 1200)
      real timediv                   ! 1 / time

! no local scalars
      integer i                         ! loop index for columns
      integer j                         ! loop index for rows
      character * 20 cmessage           ! error message for scalarmult

! declaration of externals
      external add_hours, read_one_field, read_climate_field,           &
     &         check_header,interleave

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_leads_flds'  ! subroutine name for error messages
      l_climate_field = .false.

! 1. calculate validity time of NWP data required

! DEPENDS ON: add_hours
      call add_hours(                                                   &
#include "areftim.h"
#include "avaltim.h"
     &       IVTOffHr)

!----------------------------------------------------------------------
! 2. Extract NWP field and icefrac field if available as preferred
!----------------------------------------------------------------------

      if ( LPreferred ) then
! DEPENDS ON: check_header
        call check_header (StCode,Len1_Lookup,                          &
     &                      Len2_ActualPreferred,                       &
     &                      LookupPreferred,                            &
#include "avaltim.h"
     &                       l_present_fieldNWP)
!
! DEPENDS ON: check_header
        call check_header (StCAICE,Len1_Lookup,                         &
     &                      Len2_ActualPreferred,                       &
     &                      LookupPreferred,                            &
#include "avaltim.h"
     &                       l_present_icefrac)
      endif    ! LPreferred

! 2.1 If both fields exist, read them both
      if ( l_present_fieldNWP .and. l_present_icefrac ) then
! DEPENDS ON: read_one_field
           call read_one_field (UnitPreferred,ITEM_CODE,StCode,         &
#include "avaltim.h"
     &     Len_FixHd, FixHdPreferred,Len1_Lookup,                       &
     &     Len2_ActualPreferred, LookupPreferred, LookFldNoPreferred,   &
     &     ldebug, l_climate_field,                                     &
     &     Len_IntHd, Len_RealHd, Int_Head, Real_Head,                  &
     &     ncols, nrows, fieldNWP,                                      &
#include "argppx.h"
     &     icode)

        if ( icode  <=  0) then

! 2.1.2 if NWP read successful, issue standard message
          write(UnStd,*)CStd//CSub//'NWP preferred field StCode ',      &
     &        StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'
        else
! 2.1.3 else write warning message and set l_present_preferred to false
          write(UnWarn,*)CWarn//CSub//'NWP preferred field StCode ',    &
     &        StCAICE, '; IVTOffHr = ', IVTOffHr, ' not found'
            l_present_fieldNWP = .false.
          icode = 0     ! reset icode
        endif

! 2.1.4 read icefrac field
! DEPENDS ON: read_one_field
        call read_one_field (UnitPreferred,ITEM_CODE,StCAICE,           &
#include "avaltim.h"
     &       Len_FixHd, FixHdPreferred,Len1_Lookup,                     &
     &       Len2_ActualPreferred, LookupPreferred, LookFldNoPreferred, &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, icefrac,                                     &
#include "argppx.h"
     &       icode)

        if ( icode  <=  0) then

! 2.1.5 if successful, issue standard message
          write(UnStd,*)CStd//CSub//'icefrac preferred field StCode ',  &
     &        StCAICE, '; IVTOffHr = ', IVTOffHr, ' extracted'
        else
! 2.1.6 else write warning message and set l_present_icefrac to false
          write(UnWarn,*)CWarn//CSub//                                  &
     &     'icefrac preferred field StCode ',                           &
     &     StCAICE, '; IVTOffHr = ', IVTOffHr, ' not found'
            l_present_icefrac = .false.
          icode = 0     ! reset icode
        endif
      endif    !  l_present_fieldNWP / l_present_icefrac

!----------------------------------------------------------------------
! 3. If either read fails extract previous fields if available
!----------------------------------------------------------------------
      if ( .not. l_present_fieldNWP .or.                                &
     &       .not. l_present_icefrac ) then
        if ( LPrevious ) then
! DEPENDS ON: check_header
          call check_header (StCode,Len1_Lookup,                        &
     &                      Len2_ActualPrevious,                        &
     &                      LookupPrevious,                             &
#include "avaltim.h"
     &                       l_present_fieldNWP)
!
! DEPENDS ON: check_header
          call check_header (StCAICE,Len1_Lookup,                       &
     &                      Len2_ActualPrevious,                        &
     &                      LookupPrevious,                             &
#include "avaltim.h"
     &                       l_present_icefrac)

! 3.1 If both are present, read previous fields
          if ( l_present_fieldNWP .and. l_present_icefrac ) then
! DEPENDS ON: read_one_field
            call read_one_field (UnitPrevious,                          &
     &       ITEM_CODE,StCode,                                          &
#include "avaltim.h"
     &       Len_FixHd, FixHdPrevious,Len1_Lookup,                      &
     &       Len2_ActualPrevious, LookupPrevious, LookFldNoPrevious,    &
     &       ldebug, l_climate_field,                                   &
     &       Len_IntHd, Len_RealHd, Int_Head, Real_Head,                &
     &       ncols, nrows, fieldNWP,                                    &
#include "argppx.h"
     &       icode)

            if ( icode  <=  0) then

! 3.1.1 if successful, issue standard message.

            write(UnStd,*)CStd//CSub//                                  &
     &             'NWP previous field StCode ',                        &
     &             StCode, '; IVTOffHr = ', IVTOffHr, ' extracted'

            else

! 3.1.2 else write warning message and reset icode
              write(UnWarn,*)CWarn//CSub//                              &
     &         'NWP previous field StCode ',                            &
     &         StCode, '; IVTOffHr = ', IVTOffHr, ' not found'
              l_present_fieldNWP = .false.
            end if        ! icode
            icode = 0     ! reset icode

! 3.2 Read previous icefrac field
! DEPENDS ON: read_one_field
            call read_one_field (UnitPrevious,                          &
     &        ITEM_CODE,StCAICE,                                        &
#include "avaltim.h"
     &        Len_FixHd, FixHdPrevious,Len1_Lookup,                     &
     &        Len2_ActualPrevious, LookupPrevious, LookFldNoPrevious,   &
     &        ldebug, l_climate_field,                                  &
     &        Len_IntHd, Len_RealHd, Int_Head, Real_Head,               &
     &        ncols, nrows, icefrac,                                    &
#include "argppx.h"
     &        icode)

              if ( icode  <=  0) then

! 3.2.1 if successful, issue standard message.

              write(UnStd,*)CStd//CSub//                                &
     &         'icefrac previous field StCode ',                        &
     &         StCAICE, '; IVTOffHr = ', IVTOffHr, ' extracted'

            else

! 3.2.2 else write warning message and reset icode
              write(UnWarn,*)CWarn//CSub//                              &
     &         'icefrac previous field StCode ',                        &
     &         StCAICE, '; IVTOffHr = ', IVTOffHr, ' not found'
               l_present_icefrac = .false.
            end if        ! icode
            icode = 0     ! reset icode

          endif         ! l_present_fieldNWP / l_present_icefrac
        endif       ! LPrevious
      endif     ! .not. l_present_fieldNWP .or. l_present_icefrac

!----------------------------------------------------------------------
! 4. If both fields exist, perform calculation
!----------------------------------------------------------------------
      if ( l_present_fieldNWP .and. l_present_icefrac) then

! 4.1.1 Convert units in fieldNWP if StCode is 3231
        if ( StCode  ==  3231 ) then
          timediv = 1.0 / time
! DEPENDS ON: scalarmult
          call ScalarMult (ncols, nrows, rmdi,                          &
     &            timediv, fieldNWP,                                    &
     &            fieldNWP,                                             &
     &            icode, cmessage)
        endif

      if (LClimate) then

! 4.1.2 Read climate field into fieldClim
! DEPENDS ON: read_climate_field
           call read_climate_field(StCode, IVTOffHr,                    &
     &           ldebug, Int_Head, Real_Head,                           &
     &           ncols, nrows, fieldClim,                               &
#include "argppx.h"
     &           icode)

! 4.1.2 If successful write out standard message
        if ( icode  <=  0) then
          write(UnStd,*)CStd//CSub//'Climate field extracted',          &
     &     ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
        else
! 4.1.2 If fails write out warning message
          write(UnErr,*)CErr//CSub//                                    &
     &     '4. failed to retrieve climate field ',                      &
     &     ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
          go to 9999
        end if

! 4.2 Perform calculation using interleave
! DEPENDS ON: interleave
        call interleave (ncols, nrows,                                  &
     &            fieldNWP, fieldClim,                                  &
     &            icefrac, rmdi,                                        &
     &            l_leads,field)

      else ! LClimate

      do i = 1, ncols
       do j = 1, nrows
         field(i,j) = fieldNWP(i,j)
       enddo
      enddo

      endif ! LClimate

! 4.3.  Write times to integer header
! DEPENDS ON: amend_times
        call amend_times (                                              &
#include "avaltim.h"
     &                   Int_Head,Len_IntHd )
        go to 9999
      endif     ! l_present_fieldNWP / l_present_icefrac

!----------------------------------------------------------------------
! 5. Otherwise extract field from climate file if available
!----------------------------------------------------------------------
      if (LClimate) then

! 5.1 read_climate field by calling read_climate field
! DEPENDS ON: read_climate_field
         call read_climate_field(StCode, IVTOffHr,                      &
     &           ldebug, Int_Head, Real_Head,                           &
     &           ncols, nrows, field,                                   &
#include "argppx.h"
     &           icode)

         if ( icode  <=  0) then
            write(UnStd,*)CStd//CSub//'4. climate field extracted  ',   &
     &      ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
            go to 9999
         else

            write(UnErr,*)CErr//CSub//                                  &
     &        '4. failed to retrieve climate field ',                   &
     &        ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr
            go to 9999
         end if

      end if !  LClimate/l_present_output

!----------------------------------------------------------------------
! 6. If no data has been successfully extracted return an error code
!----------------------------------------------------------------------
      icode = 5
      write(UnErr,*)CErr//CSub//'5. failed to extract any data',        &
     &    ' for stash code =', stcode, '; IVTOffHr = ', IVTOffHr

9999  continue
      return
      END SUBROUTINE read_leads_flds
!----------------------------------------------------------------------
#endif
