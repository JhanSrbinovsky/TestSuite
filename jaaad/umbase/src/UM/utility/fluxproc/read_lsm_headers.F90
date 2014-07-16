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
! Author:     M. J. Bell
! 5.3      14/08/01     Allow atmosphere fields on B or C grids MJ
! 6.1      30/07/04     Add setting of LCyclicO. A. Hines
!----------------------------------------------------------------------
! contains routines: read_lsm_headers
!
! Purpose: Flux processing routine.
!          Opens and reads lookup tables for land sea masks used by
!          FOAM_Flux_Process. Also sets LCyclic = T if atmosphere
!          grid has wrap-round points, and LCyclicO = T if
!          ocean grid has wrap-round points.
!----------------------------------------------------------------------
      subroutine read_lsm_headers (                                     &
#include "aflddims.h"
     &    ppxRecs,icode)

      implicit none

! declaration of argument list
! dimensions of ocean and atmosphere fields
#include "cflddims.h"
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of parameters
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "plookups.h"
#include "clookadd.h"

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"
#include "clookups.h"
#include "cselct.h"

! no local arrays

! declaration of local scalars
      integer Len2_Lookup_lsm     ! max 2nd dimension for lsms
      integer Len2_Lookup_Actual  ! actual 2nd dimension for lsms
      integer IROW_NUMBER
      character*80 cmessage

      external read_one_header
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_lsm_headers'! subroutine name for error messages
      Len2_Lookup_lsm = 1      ! all lsm ancillary files contain 1 field

! 0.1 Read StashMaster files
      IROW_NUMBER=0
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_A',IROW_NUMBER,                     &
#include "argppx.h"
     &  ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_O',IROW_NUMBER,                     &
#include "argppx.h"
     &  ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_S',IROW_NUMBER,                     &
#include "argppx.h"
     &  ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(22,2,'STASHmaster_W',IROW_NUMBER,                     &
#include "argppx.h"
     &  ICODE,CMESSAGE)


! 1. read atmosphere tracer land / sea mask fixed header and lookup
!    table from an an ancillary file
! DEPENDS ON: read_one_header
      call read_one_header(UnitNWPlsmt, icode,                          &
     &               Len_FixHd, Len1_Lookup, Len2_Lookup_lsm,           &
     &               Len2_Lookup_Actual, FixHdlsmt,                     &
#include "argppx.h"
     &               Lookuplsmt)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 1. unable to read NWP tracer land sea mask headers'
        go to 9999
      end if

! 1.1 extract the number of rows and columns from the lookup table
      ncols  = Lookuplsmt(LBNPT)
      nrowst = Lookuplsmt(LBROW)

! 2. Set dimensions for atmosphere velocity grids.
!    Calculations differ for B and C grids.
!    Note that Lookuplsmu and set_lookups_u are no longer used.
      if ( l_B_grid) then
        nrowsu = nrowst - 1
        nrowsv = nrowsu
        nrowsuv = nrowst - 1
      else
        nrowsu = nrowst
        nrowsv = nrowsu - 1
        nrowsuv = nrowst - 1
      end if

! 3. Set LCyclic (T if atmosphere grid has wrap points)
!    if fixhd(4)
      if ( MOD ( FixHdlsmt (4) , 100 )  /=  3 ) then
        LCyclic = .True.
      else
        LCyclic = .False.
      end if

! 4. read ocean tracer land / sea mask lookup table
! DEPENDS ON: read_one_header
      call read_one_header(UnitFOAMlsmt, icode,                         &
     &               Len_FixHd, Len1_Lookup, Len2_Lookup_lsm,           &
     &               Len2_Lookup_Actual, FixHdlsmtO,                    &
#include "argppx.h"
     &               LookuplsmtO)


      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &       ' step 4. unable to read ocean tracer land sea mask '
        go to 9999
      end if

! 4.1 extract the number of rows and columns from the lookup table
      ncolsO  = LookuplsmtO(LBNPT)
      nrowstO = LookuplsmtO(LBROW)

! 5. read ocean velocity land / sea mask lookup table
! DEPENDS ON: read_one_header
      call read_one_header(UnitFOAMlsmu, icode,                         &
     &               Len_FixHd, Len1_Lookup, Len2_Lookup_lsm,           &
     &               Len2_Lookup_Actual, FixHdlsmuO,                    &
#include "argppx.h"
     &               LookuplsmuO)

      if ( icode  >   0 ) then
        write(UnErr,*)CErr,CSub,                                        &
     &    ' step 5. unable to read ocean velocity land sea mask '
        go to 9999
      end if

! 5.1 extract the number of rows and columns from the lookup table
      nrowsuO = LookuplsmuO(LBROW)

! 6. Set LCyclicO (T if ocean grid has wrap points)
!    if fixhd(4)
      if ( MOD ( FixHdlsmtO (4) , 100 )  /=  3 ) then
        LCyclicO = .True.
      else
        LCyclicO = .False.
      end if

9999  continue
      return
      END SUBROUTINE read_lsm_headers
!----------------------------------------------------------------------
#endif
