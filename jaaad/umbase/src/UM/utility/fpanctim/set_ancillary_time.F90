#if defined(FLUXPROC)
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
! contains progam:   set_ancillary_time
!
! Purpose: Flux processing routine.
!          Writes namelists which determine the validity times in the
!          fixed headers of the FOAM flux ancillary files
!----------------------------------------------------------------------
      Program set_ancillary_time

      implicit none

! declaration of globals

#include "cunitnos.h"
#include "creftim.h"

! declaration of local scalars

      integer fvhh,fvdd,fvmm,fvyy    ! first validity time
      integer lvhh,lvdd,lvmm,lvyy    ! last validity time
      integer ivhh,ivdd,ivmm,ivyy    ! interval between fields in file

      logical year360                ! true for 360-day calendar

      integer first_vt_offset        ! offset in hours of first validity
                                     ! time from reference time

      integer last_vt_offset         ! offset in hours of first validity
                                     ! time from reference time

      integer imins, isecs           ! minutes and seconds discarded

      integer icode                  ! error return code; > 0 => error


! declaration of namelists

      namelist /offset_vt/ first_vt_offset,last_vt_offset

      namelist /first_vt/ fvhh,fvdd,fvmm,fvyy
      namelist /last_vt/  lvhh,lvdd,lvmm,lvyy
      namelist /interval/ year360,ivhh,ivdd,ivmm,ivyy

! declaration of externals
      external open_file, readhk, add_hours
!----------------------------------------------------------------------

! 0. Preliminaries
      icode = 0

! 1. Open files for input and output

! 1.1 open housekeeping file
! DEPENDS ON: open_file
      call open_file(UnitHK, 'Formatted  ', 'Unknown', icode)
      if ( icode  /=  0 ) then
        write(6,*) ' ERROR: set_ancillary_time: 1.1 ',                  &
     &             ' unable to open housekeeping file'
        go to 9999
      end if


!  open Validity time selection file
! DEPENDS ON: open_file
      call open_file(UnitVT, 'Formatted  ', 'Unknown', icode)
      if ( icode  /=  0 ) then
        write(6,*) ' ERROR: set_ancillary_time: 1.2 ',                  &
     &             ' unable to open Validity time selection file'
        go to 9999
      end if

! open output file for ancillary file validity times namelist
! DEPENDS ON: open_file
      call open_file(UnitVTOut, 'Formatted  ', 'Unknown', icode)
      if ( icode  /=  0 ) then
        write(6,*) ' ERROR: set_ancillary_time: 1.3 ',                  &
     &             ' unable to open Validity times output file'
        go to 9999
      end if


! 1. Read house keeping file to set reference date
      RefSec = 0
      RefMin = 0
! DEPENDS ON: readhk
      call readhk(UnitHK, RefHour, RefDay, RefMonth, RefYear,icode)

! Read offset_vt  namelist
      first_vt_offset =  0
      last_vt_offset =  0

      read ( UnitVT, offset_vt )
      write ( 6, offset_vt )

! Read interval namelist

      year360=.false.
      ivhh=0
      ivdd=0
      ivmm=0
      ivyy=0

      read ( UnitVT, interval )
      write ( 6, interval )

! Set first validity time
! DEPENDS ON: add_hours
      call add_hours (                                                  &
#include "areftim.h"
     &  fvyy, fvmm, fvdd, fvhh, imins, isecs,                           &
     &  first_vt_offset )

      write (6, first_vt)

! Set last validity time
! DEPENDS ON: add_hours
      call add_hours (                                                  &
#include "areftim.h"
     &  lvyy, lvmm, lvdd, lvhh, imins, isecs,                           &
     &  last_vt_offset )

      write (6, last_vt)

! Write output namelists

      write ( UnitVTOut, first_vt )
      write ( UnitVTOut, interval )
      write ( UnitVTOut, last_vt )

! close files
      close (UnitHK)
      close (UnitVT)
      close (UnitVTOut)

9999  continue

      stop
      END PROGRAM set_ancillary_time
!----------------------------------------------------------------------
#endif
