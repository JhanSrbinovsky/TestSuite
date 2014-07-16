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
! 5.3      22/10/01     New routines to allow separation of
!                       flux processing and interpolation. A.Hines.
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines:  close_files, close_flux_files, close_grid_files
!
! Purpose: Flux processing routine.
!          Closes all files used by FOAM_Flux_Process
!----------------------------------------------------------------------
#if !defined(FLXPLPR) && !defined(FLXPLIN)
!----------------------------------------------------------------------
#endif
#if defined(FLXPLPR)
!----------------------------------------------------------------------
#endif
#if defined(FLXPLIN)
      subroutine close_grid_files
!----------------------------------------------------------------------

      implicit none

! no arguments

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"
#include "cdebug.h"

! No local arrays

! declaration of local scalars
      integer iun   ! loop index for unit number
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'close_grid_files'  ! subroutine name for error messages

! 1. Close input and output flux files ! will need changing !
      do iun = IUnOutLow, IUnOutHi
        if ( LUnOutOpen(iun) ) then
          close ( iun )
        end if ! LUnOutOpen(iun)
      end do ! iun


! 2. Close land/sea mask files
      close ( UnitNWPlsmt )
      close ( UnitFOAMlsmt )
      close ( UnitFOAMlsmu )

! 4. Close control files
      close ( UnitDbg )

! 5. Close log and debug output files
      close ( OutUnitDbg )
      close ( UnStd )
      close ( UnWarn )

      return
      END SUBROUTINE close_grid_files
!----------------------------------------------------------------------
#endif
#endif
