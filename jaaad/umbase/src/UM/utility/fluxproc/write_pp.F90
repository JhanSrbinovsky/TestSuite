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
! 5.3      10/10/00     Changes to enable flux selection and grid
!                       interpolation to be split between 2 programs.
!                       M. J. Bell and A. Hines.
! 6.1      30/07/04     Correct call to spiral_s to use LCyclicO.
!                       A. Hines.
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: write_one_field,write_pp
!
! Purpose: Flux processing routine.
!          write_one_field:
!          This routine writes out one pp-field to the required output
!          pp file. The pp-header is atered for the correct stashcode
!          and the field written to the file using the routine write_pp.
!          The routine also sets land points to missing data and checks
!          for consistency in the field dimensions.
!          write_pp: Writes out a pp header then a pp field
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      subroutine write_pp ( IOutUnit, Int_Head, Real_Head,              &
     &                      ncolsOut, nrowsOut, field_out, icode)



      implicit none

! declaration of parameters
#include "plookups.h"

! declaration of argument list
      integer IOutUnit  ! IN output unit number
      integer Int_Head(Len_IntHd)  ! integer part of lookup table
      real Real_Head(Len_RealHd)   ! real part of lookup table
      integer ncolsOut             ! # of columns in output field
      integer nrowsOut             ! # of rows in output field
      real field_out( ncolsOut, nrowsOut ) ! field output
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"

! no local arrays

! declaration of local scalars
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'write_pp'  ! subroutine name for error messages

! 1. Write out header
      write (IOutUnit, IOStat = icode) Int_Head, Real_Head
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1. error writing out lookup table  '
        icode = 47
        go to 9999
      end if

! 2. Write out data
      write (IOutUnit, IOStat = icode) field_out
      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2. error writing out data field  '
        icode = 48
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE write_pp
!----------------------------------------------------------------------
#endif
