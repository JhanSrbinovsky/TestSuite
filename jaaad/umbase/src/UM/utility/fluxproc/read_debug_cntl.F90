#if defined(FLUXPROC) || defined(FLXPLPR) || defined(FLXPLIN)
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
! contains routines: read_debug_cntl
!
! Purpose: Flux processing routine.
!          Reads files controlling diagnostic "debugging" output
!----------------------------------------------------------------------
      subroutine read_debug_cntl ( icode )

      implicit none

! declaration of argument list
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of parameters

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"
#include "cdebug.h"

! No local arrays

! declaration of local scalars
      integer i  ! loop index

! namelist declaration
      NAMELIST /NmLstDbg/                                               &
     &           NoDbgPts,                                              &
     &           IColDbg, JRowDbg,                                      &
     &           l_winds_dbg, l_heat_dbg, l_moisture_dbg,               &
     &           l_sea_ice_dbg, l_references_dbg, l_pressure_dbg,       &
     &           l_windspd_dbg

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_debug_cntl'  ! subroutine name for error messages

! 1. set default values for variables in NmLstDbg
      NoDbgPts = 0

      do i = 1, MaxNoDbgPts
        IColDbg(i) = 1
        JRowDbg(i) = 1
      end do

      l_winds_dbg = .false.
      l_heat_dbg = .false.
      l_moisture_dbg = .false.
      l_sea_ice_dbg = .false.
      l_references_dbg = .false.
      l_pressure_dbg = .false.
      l_windspd_dbg = .false.

! 2. read debug control namelist
      read (UnitDbg, NmLstDbg, iostat = icode)

      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &     ' step 2. unable to read debug control namelist'
        icode = 8
        go to 9999
      end if

! 3. open file for debugging output
! DEPENDS ON: open_file
      call Open_file(OutUnitDbg,'Formatted  ','Unknown',icode)
      if (icode  /=  0) then
        write(UnErr,*)CErr,CSub,                                        &
     &     ' step 3. unable to open file for debugging output'
        icode = 9
        go to 9999
      end if

! 4. read and write out contents of namelist
      write(OutUnitDbg, NmLstDbg)

      if ( NoDbgPts  <=  0) then
        write(OutUnitDbg, *) ' no points to output '
      else
        write(OutUnitDbg, *) ' columns of output: '
        write(OutUnitDbg, '(11I11)' ) (IColDbg(i), i=1,NoDbgPts)
        write(OutUnitDbg, *) ' rows of output: '
        write(OutUnitDbg, '(11I11)' ) (JRowDbg(i), i=1,NoDbgPts)
      end if

9999  continue
      return
      END SUBROUTINE read_debug_cntl
!----------------------------------------------------------------------
#endif
