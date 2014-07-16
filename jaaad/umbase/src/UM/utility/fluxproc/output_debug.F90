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
! contains routines: output_debug
!
! Purpose: Flux processing routine.
!          Writes out field values at selected grid points
!----------------------------------------------------------------------
      subroutine output_debug(CMessage, nrows, ncols, Field)

      implicit none

! argument list: all intent IN
      character*(*) CMessage    ! message to output
      integer    nrows          ! input dimension; number of rows
      integer    ncols          ! input dimension; number of columns
      real       Field(ncols,nrows)   ! field to output

! globals
#include "cmess.h"
#include "cdebug.h"

! local scalars
      integer IDbg   ! loop index over selected points to debug

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'output_debug'  ! subroutine name for error messages


! 1. write out field descriptor
      if ( CMessage  /=  '         ') then
        write(OutUnitDbg, *) CMessage
      end if

! 2. convert values to output
      do IDbg = 1, NoDbgPts

        if ( IColDbg(IDbg)  <=  ncols .and.                             &
     &       JRowDbg(IDbg)  <=  nrows       ) then
           write(CValues(IDbg), '(G11.4)' )                             &
     &                          Field( IColDbg(IDbg), JRowDbg(IDbg) )
         else
           CValues(IDbg) = ' OOB '
         end if

      end do ! IDbg

! 3. output values
        write(OutUnitDbg, '(11A11)' ) (CValues(IDbg), IDbg=1,NoDbgPts)

      return
      END SUBROUTINE output_debug
!----------------------------------------------------------------------
#endif
