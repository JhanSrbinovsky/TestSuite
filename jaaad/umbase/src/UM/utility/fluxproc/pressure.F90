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
! 5.3      22/10/01     New if defined added. A. Hines
!
! Author:     S. A. Spall
!----------------------------------------------------------------------
! contains routines: pressure
!
! Purpose: Flux processing routine.
!          To produce a pp file containing:
!          Sea Surface pressure for the times required
!----------------------------------------------------------------------
      subroutine pressure(                                              &
#include "afields.h"
#include "argppx.h"
     &                 icode )

      implicit none

! declaration of argument list

! array dimensions, lsms, interpolation coeffs etc. : all intent IN
#include "cfields.h"

      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of parameters
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cfdcodes.h"
#include "plookups.h"

! declaration of globals used
#include "clookadd.h"
#include "cunitnos.h"
#include "cmess.h"
#include "c_mdi.h"
#include "cvaloff.h"
#include "cdebug.h"
#include "c_0_dg_c.h"

! declaration of local arrays
      integer Int_Head_SSP(Len_IntHd)  ! integer part of lookup table
      real Real_Head_SSP(Len_RealHd)   ! real part of lookup table
      real sea_surface_pressure(ncols, nrowst) ! ref SSP

! declaration of local scalars

      integer ivt           ! loop index over validity times
      integer IVTOffHr      ! offset of validity time from reference
      integer IOutUnit      ! output unit

      logical ldebug        ! T => output debugging info (set in 0.)

      character * 256 cmessage   ! error message


! declaration of externals
      external read_fields, write_one_field

!----------------------------------------------------------------------
! 0. Preliminaries
!----------------------------------------------------------------------
      CSub = 'pressure'  ! subroutine name for error messages

      ldebug = l_pressure_dbg    ! set by debug input control file

!----------------------------------------------------------------------
! 1. start loop over validity times
!----------------------------------------------------------------------
      do ivt = 1, NoValidTimes

        IVTOffHr = IValidOffHr(ivt)
        IOutUnit = IOutUnitOff(ivt) + UnitPressureOut

!----------------------------------------------------------------------
! 2. read in sea surface pressure
!----------------------------------------------------------------------
! DEPENDS ON: read_fields
        call read_fields(StCSSP, IVTOffHr,                              &
     &               ldebug, Int_Head_SSP, Real_Head_SSP,               &
     &               ncols, nrowst,                                     &
     &               sea_surface_pressure,                              &
#include "argppx.h"
     &               icode)

        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 2. unable to read sea surface pressure'
          icode = 1121
          go to 9999
        end if

!----------------------------------------------------------------------
! 3. Write out sea surface pressure
!----------------------------------------------------------------------
! DEPENDS ON: write_one_field
        call write_one_field (                                          &
#include "afields.h"
     &       OutStCSSP, FFSSP, PPSSP, IVTOffHr,                         &
     &       ITGrid, nrowst,                                            &
     &       Int_Head_SSP, Real_Head_SSP, IOutUnit, ldebug,             &
     &       sea_surface_pressure, icode)
        if ( icode  >   0 ) then
          write(UnErr,*)CErr,CSub,                                      &
     &       ' step 3. unable to write sea surface pressure'
          icode = 1122
          go to 9999
        end if

!----------------------------------------------------------------------
! 4. end loop over validity times
!----------------------------------------------------------------------
        enddo    !  ivt

9999  continue
      return
      END SUBROUTINE pressure
!----------------------------------------------------------------------
#endif
