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
! 5.3      10/10/00     Split to enable flux selection and grid
!                       interpolation to be split between 2 programs.
!                       M. J. Bell and A. Hines.
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines:  open_flux_control_files, open_grid_control_files,
!                     open_control_files
!
! Purpose: Flux processing routine.
!          Opens all control and log files used by FOAM_Flux_Process
!----------------------------------------------------------------------
#if !defined(FLXPLPR) && !defined(FLXPLIN)
!----------------------------------------------------------------------
#endif
#if defined(FLXPLPR)
      subroutine open_flux_control_files ( icode )

! Purpose: Opens all control and log files used by Flux_Process_Main

      implicit none

! declaration of argument list
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"
#include "cenviron.h"

! No local arrays

! declaration of local scalars
      integer i  ! do loop index
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'open_flux_control_files' ! subroutine name for errors

! 0.1 Before opening any files set values in CENVIRON
      do i = 1, NUNITS
        CEnv(i) = ' '
        LEnv(i) = 1
      end do

      CEnv(UnitPreferred) = 'FFPREFERRED'
      LEnv(UnitPreferred) = 11
      CEnv(UnitPrevious)  = 'FFPREVIOUS'
      LEnv(UnitPrevious)  = 10
      CEnv(UnitClimate)   = 'FFCLIMATE'
      LEnv(UnitClimate)   = 9

! 1. open log files

! 1.1 open error, warning and standard log files
! DEPENDS ON: open_file
      call open_file(UnErr, 'Formatted  ', 'Unknown', icode)
      if ( icode  >   0 ) then
        write(6,*)' ERROR opening error log file in FOAM flux '
        write(6,*)' processing. This job will have failed.  '
        write(6,*)CErr,CSub,' step 1.1 unable to open error log '
        icode = 1
        go to 9999
      end if

! 1.2 open warning log file
! DEPENDS ON: open_file
      call open_file(UnWarn, 'Formatted  ', 'Unknown', icode)
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1.2 unable to open warning log '
        icode = 2
        go to 9999
      end if

! 1.3 open standard log file
! DEPENDS ON: open_file
      call open_file(UnStd, 'Formatted  ', 'Unknown', icode)
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1.3 unable to open standard log '
        icode = 3
        go to 9999
      end if

! 2. open control files

! 2.1 open housekeeping file
! DEPENDS ON: open_file
      call open_file(UnitHK, 'Formatted  ', 'Unknown', icode)
      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.1 unable to open housekeeping file'
        icode = 4
        go to 9999
      end if

! 2.2 open Validity time selection file
! DEPENDS ON: open_file
      call open_file(UnitVT, 'Formatted  ', 'Unknown', icode)
      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &   ' step 2.2 unable to open Validity time selection file'
        icode = 5
        go to 9999
      end if

! 2.3 open debug control file
! DEPENDS ON: open_file
      call open_file(UnitDbg, 'Formatted  ', 'Unknown', icode)
      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.3 unable to open debug control file'
        icode = 6
        go to 9999
      end if

! 2.4 open flux selection control file
! DEPENDS ON: open_file
      call open_file(UnitSlt, 'Formatted  ', 'Unknown', icode)
      if ( icode  /=  0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2.4 unable to open select control file'
        icode = 7
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE open_flux_control_files
!----------------------------------------------------------------------
#endif
#if defined(FLXPLIN)
#endif

#endif
