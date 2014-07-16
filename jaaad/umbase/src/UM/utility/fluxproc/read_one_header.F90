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
! contains routines: read_one_header
!
! Purpose: Flux processing routine.
!          Opens and reads fixed header and lookup table of one file
!----------------------------------------------------------------------
      subroutine read_one_header ( InUnit, ICode,                       &
     & Len_FixHd_P, Len1_Lookup_P, Len2_Lookup_P,                       &
     & Len2_Lookup_Actual, FixHd,                                       &
#include "argppx.h"
     & Lookup)

      implicit none

! declaration of argument list
      integer InUnit  ! IN     input unit number
      integer ICode   ! IN/OUT error code ; > 0 => fatal error detected

! dimensions used to declare arrays to be read in
      integer Len_FixHd_P   ! IN length of fixed header
      integer Len1_Lookup_P ! IN length of first dimension of Lookup
      integer Len2_Lookup_P ! IN max length of 2nd dimension of Lookup

! fixed header and lookup tables: intent OUT
      integer Len2_Lookup_Actual    ! OUT actual 2nd dimension of Lookup
      integer FixHd(Len_FixHd_P)                     ! OUT fixed header
      integer Lookup(Len1_Lookup_P, Len2_Lookup_P)   ! OUT lookup tables

!  declaration of globals
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cmess.h"
#include "cenviron.h"

! declaration of local arrays

! lengths of headers in fields file (local arrays)
! (this declares LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS)
#include "dump_len.h"

! declaration of local scalars
      integer Len_data  ! length of data in file
      character*256 CMessage ! error messages

      external READ_FLH, GET_DIM, setpos, get_lookup

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_one_header'  ! subroutine name for error messages
      CMessage = ' '

! 1 open  file
! DEPENDS ON: file_open
      call file_open(InUnit, CEnv(InUnit), LEnv(InUnit), 0, 0, icode)

      if (icode  >   0) then
        write(UnWarn,*)CWarn,CSub,                                      &
     &  ' step 1. failed to open file with environment name ',          &
     &  CEnv(InUnit)
        icode = 27
        go to 9999
      end if

! 2. Read fixed header
! DEPENDS ON: read_flh
      CALL READ_FLH(InUnit,FIXHD,LEN_FIXHD_P,ICODE,CMESSAGE)

      if ( icode  >   0) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 2. unable to read fixed header; cmessage is ',      &
     &       cmessage
        icode = 28
        go to 9999
      end if

! 3. get dimensions from lookup table, check them and set actual
!    2nd dimension of lookup table

! 3.0  get dimensions from lookup table
      LEN_FIXHD = LEN_FIXHD_P
! DEPENDS ON: get_dim
      CALL GET_DIM(FIXHD,                                               &
#include "dump_ar2.h"
     & Len_data)

! 3.1 Set to zero any dimensions of headers which are less than zero
!     (readhead etc. fail if this is not done)

      if ( LEN1_COLDEPC  <   0) LEN1_COLDEPC = 0
      if ( LEN2_COLDEPC  <   0) LEN2_COLDEPC = 0
      if ( LEN1_FLDDEPC  <   0) LEN1_FLDDEPC = 0
      if ( LEN2_FLDDEPC  <   0) LEN2_FLDDEPC = 0
      if ( LEN_EXTCNST  <   0)  LEN_EXTCNST  = 0
      if ( LEN_DUMPHIST  <   0) LEN_DUMPHIST = 0
      if ( LEN_CFI1  <   0)     LEN_CFI1 = 0
      if ( LEN_CFI2  <   0)     LEN_CFI2 = 0
      if ( LEN_CFI3  <   0)     LEN_CFI3 = 0

! 3.2 check lookup 2nd dimensions are not too large
      if ( Len2_Lookup_P  <   Len2_Lookup_Obs ) then
        write(UnErr,*)CErr,CSub,                                        &
     &  ' step 3.2 lookup table is not big enough; Len2_Lookup_P = ',   &
     &  Len2_Lookup_P,'; Len2_Lookup = ', Len2_Lookup_Obs
        icode = 29
        go to 9999
      end if

! 3.3 check first dimensions of lookup tables match
      if ( Len1_Lookup_P  /=  Len1_Lookup_Obs ) then
        write(UnErr,*)CErr,CSub,                                        &
     &  ' step 3.3 lookup first dimensions do not match ;',             &
     &  ' Len1_Lookup_P = ', Len1_Lookup_P,'; Len1_Lookup = ',          &
     &  Len1_Lookup_Obs
        icode = 30
        go to 9999
      end if

! 3.4 set actual 2nd dimension of lookup table
      Len2_Lookup_Actual = Len2_Lookup_Obs

! 4. set position to start of file to re-read the header
! DEPENDS ON: setpos
      call setpos(InUnit, 0, icode)

      if ( icode  >   0) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 4. setpos failed; icode = ', icode
        icode = 31
        go to 9999
      end if

! 5. get the lookup header

! DEPENDS ON: get_lookup
      call get_lookup(InUnit, icode,                                    &
#include "dump_ar2.h"
#include "argppx.h"
     & Len_data, LOOKUP)

      if ( icode  >   0) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 5. failed read lookup table '
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE read_one_header
!----------------------------------------------------------------------
#endif
