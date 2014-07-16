#if defined(FLUXPROC) || defined(FLXPLIN)
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
! contains routines: get_rows_cols
!
! Purpose: Flux processing routine.
!          Returns row and column dependent constants for one file
!----------------------------------------------------------------------
      subroutine get_rows_cols (UnitIn, icode,                          &
#include "dump_ar2.h"
#include "argppx.h"
     & Len_data, rowdepc, coldepc)

      implicit none

! declaration of argument list
      integer UnitIn  ! IN unit number of file to read
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected
! all arguments in DUMP_LEN are IN
#include "dump_len.h"

      integer Len_data ! IN length of data in file

! rowdepc, coldepc: intent is OUT; they are declared by DUMP_DIM

! no parameters

! declaration of globals used
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cmess.h"

! declaration of local arrays
! DUMP_DIM  declares local arrays for ancillary file header
#include "dump_dim.h"

! declaration of local scalars
      integer START_BLOCK       !  start of data block
      character*256 CMESSAGE    !  error message

      external READHEAD

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'get_rows_cols'  ! subroutine name for error messages

! 1. Read dump / ancillary file header

! DEPENDS ON: readhead
      call READHEAD(UnitIn,                                             &
#include "dump_ar1.h"
     & Len_data,                                                        &
#include "argppx.h"
     & START_BLOCK,ICODE,CMESSAGE)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &  ' step 1. failed to read ancillary file header; cmessage is ',  &
     &  cmessage
        icode = 22
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE get_rows_cols
!----------------------------------------------------------------------
#endif
