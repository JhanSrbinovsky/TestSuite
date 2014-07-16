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
! contains routines: read_control_files
!
! Purpose: Flux processing routine.
!          Reads all control files used by FOAM_Flux_Process
!          Units added for pressure and windspeed (S. Spall)
!----------------------------------------------------------------------
      subroutine read_control_files (icode)

      implicit none

! declaration of argument list
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"

#include "creftim.h"
#include "cvaloff.h"

#include "c_mdi.h"

! declaration of  local arrays
      integer iunit_base(7)   ! main unit numbers for output files

! declaration of local scalars
      integer ivt       ! loop index over validity times
      integer iunit     ! loop index over unit numbers
      integer IAdd      ! value to add to basic unit number
      integer IUnitOpen ! unit number to open

! namelist declaration
      NAMELIST / NamFluxSelect /                                        &
     &  ValidityPeriod,                                                 &
     &  NoValidTimes, IValidOffHr, IOutUnitOff,                         &
     &  NoAddTimesPreferred, ISrchOffHrPreferred, INewOffHrPreferred,   &
     &  NoAddTimesPrevious, ISrchOffHrPrevious, INewOffHrPrevious,      &
     &  NoAddTimesClimate, ISrchOffHrClimate, INewOffHrClimate,         &
     &  output_land_value

! declaration of external subroutines and functions
      external readhk, read_debug_cntl, open_file

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_control_files' ! subroutine name for error messages

! 1. Read house keeping file to set reference date
      RefSec = 0
      RefMin = 0
! DEPENDS ON: readhk
      call readhk(UnitHK, RefHour, RefDay, RefMonth, RefYear,icode)

      if ( icode  /=  0 ) then
        write (UnErr,*)CErr,CSub,                                       &
     &   '1. Failed to read housekeeping file'
        goto 9999
      endif
      write(UnStd,*) CStd,CSub,'reference time from housekeeping file:' &
     & , ' RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec = ',      &
     & RefYear, RefMonth, RefDay, RefHour, RefMin, RefSec

! 2. Read debug control file and open debug ouput file
! DEPENDS ON: read_debug_cntl
      call read_debug_cntl ( icode )
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 2.  Failed to read debug control file'
        go to 9999
      end if

! 2.1 Read select control file
! DEPENDS ON: read_select_cntl
      call read_select_cntl ( icode )
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 2.1  Failed to read select control file'
        go to 9999
      end if

! 3. Read validity times control file

! 3.0 Set defaults
      NoAddTimesPreferred = 0
      NoAddTimesPrevious  = 0
      NoAddTimesClimate   = 0

      do ivt = 1, MaxTimes
        IValidOffHr(ivt) = 0
        IOutUnitOff(ivt) = 0
        ISrchOffHrPreferred(ivt) = 0
        INewOffHrPreferred(ivt) = 0
        ISrchOffHrPrevious(ivt) = 0
        INewOffHrPrevious(ivt) = 0
        ISrchOffHrClimate(ivt) = 0
        INewOffHrClimate(ivt) = 0
      end do

      output_land_value = rmdi

! 3.1 read namelist

      read (UnitVT, NamFluxSelect,  iostat = icode)
      if ( icode  >   0) then
        write(UnErr,*)CErr,CSub,                                        &
     & ' step 3.1  Failed to read validity times control file'
        icode = 10
        go to 9999
      end if

      write(UnStd, NamFluxSelect)

! 4. set which units to open for output flux files

      do iunit = IUnOutLow, IUnOutHi
        LUnOutOpen(iunit) = .False.
      end do

      iunit_base(1) = UnitWindsOut
      iunit_base(2) = UnitHeatOut
      iunit_base(3) = UnitMoistureOut
      iunit_base(4) = UnitSeaIceOut
      iunit_base(5) = UnitReferencesOut
      iunit_base(6) = UnitPressureOut
      iunit_base(7) = UnitWindspdOut

      do ivt = 1, NoValidTimes
        IAdd = IOutUnitOff(ivt)
        do iunit = 1, 7
          IUnitOpen = iunit_base(iunit) + IAdd
          if ( IUnitOpen  <   IUnOutLow .or.                            &
     &        IUnitOpen  >   IUnOutHi ) then
            icode = 11
            write(UnErr,*)CErr,CSub,' step 4. Unit number chosen'       &
     &      ,' incorrectly; ivt,iunit =',ivt,iunit
            go to 9999
          else
            LUnOutOpen(IUnitOpen) = .True.
          end if
        end do ! iunit
      end do ! ivt

! 5. open output flux files
      do iunit = IUnOutLow, IUnOutHi
        if ( LUnOutOpen(iunit) ) then

! DEPENDS ON: open_file
          call open_file ( iunit, 'unformatted', 'unknown', icode )
          write(UnStd,*)CStd,CSub, ' step 5.  Opening file ', iunit

          if ( icode  >   0) then
            write(UnErr,*)CErr,CSub,                                    &
     &      ' step 5.  Failed to open output flux file ', iunit
            icode = 12
            go to 9999
          end if  ! icode

        end if ! LUnOutOpen(iunit)
      end do ! iunit

9999  continue
      return
      END SUBROUTINE read_control_files
!----------------------------------------------------------------------
#endif
