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
! 5.3      10/10/00     Changes to enable flux selection and grid
!                       interpolation to be split between 2 programs.
!                       M. J. Bell and A. Hines.
! 6.1      29/07/04     Correction to setting of LCyclic. A. Hines.
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: read_field_headers, read_field_sizes
!
! Purpose: Flux processing routine.
!          Reads fixed headers and lookup tables of  flux files input
!          to FOAM_Flux_Process (i.e. Preferred or Previous fluxes and
!          climate fluxes)
!          Also works out dimensions ncols, nrowsu, nrowst
!----------------------------------------------------------------------
      subroutine read_field_headers (                                   &
#include "aflddima.h"
     &                                ppxRecs,icode )

      implicit none

! declaration of argument list
#include "cflddima.h"
#include "cselct.h"
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of parameters
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "plookups.h"

! declaration of globals used
#include "cunitnos.h"
#include "cmess.h"
#include "clookups.h"
#include "cvaloff.h"
! declaration of logicals
      logical l_climate_field       ! T => Climate Field being used
      integer IROW_NUMBER
      character*80 cmessage

      external read_one_header, add_lookups
!----------------------------------------------------------------------
! 0. Preliminaries
!----------------------------------------------------------------------
      CSub = 'read_field_headers' ! subroutine name for error messages

#if defined(FLXPLPR)
! 0.1 Set null values for row and column dimensions
      ncols=0
      nrowst=0
      nrowsu=0
      nrowsv=0
      nrowsuv=0
#endif

! 0.2 Read StashMaster files
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

!----------------------------------------------------------------------
! 1. Read and amend fixed header and lookups of preferred file
!----------------------------------------------------------------------


! 1.0 read headers

      LPreferred = .True.
! DEPENDS ON: read_one_header
      call read_one_header(UnitPreferred, icode,                        &
     &   Len_FixHd, Len1_Lookup, Len2_LookupPreferred,                  &
     &   Len2_ActualPreferred, FixHdPreferred,                          &
#include "argppx.h"
     &   LookupPreferred)

      if (icode  /=  0) then
        LPreferred = .False.
        write(UnWarn,*)CWarn,CSub,                                      &
     &  ' step 1.0 unable to open and read headers of' //               &
     &  ' preferred flux file '
        icode = 0
      end if

! 1.1 amend headers

      if ( LPreferred ) then

#if defined(FLXPLPR)
! 1.1.0 set LCyclic (T if atmosphere grid has wrap points)
      if ( MOD ( FixHdPreferred (4) , 100 )  /=  3 ) then
        LCyclic = .True.
      else
        LCyclic = .False.
      end if

! 1.1.1 attempt to set field dimensions
! DEPENDS ON: read_field_sizes
        call read_field_sizes(                                          &
     &      Len1_Lookup, Len2_ActualPreferred, LookupPreferred,         &
     &      ncols, nrowst, nrowsu, icode)
        if ( icode  >   0) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 1.1 read_field_sizes failed for preferred file'
          go to 9999
        end if

! 1.1.2 Set dimensions for atmosphere velocity grids.
!       Calculations differ for B and C grids.
      if ( l_B_grid) then
        nrowsv = nrowsu
        nrowsuv = nrowst - 1
      else
        nrowsv = nrowsu - 1
        nrowsuv = nrowst - 1
      end if

#endif

! 1.1.2 amend headers

        l_climate_field = .false.

! DEPENDS ON: add_lookups
        call add_lookups (                                              &
     &   NoAddTimesPreferred, ISrchOffHrPreferred, INewOffHrPreferred,  &
     &   l_climate_field, Len1_Lookup, Len2_LookupPreferred,            &
     &   Len2_ActualPreferred,                                          &
     &   LookupPreferred, LookFldNoPreferred, icode )

        if ( icode  >   0) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 1.1 add_lookups failed for preferred file'
          go to 9999
        end if

      end if   ! LPreferred


!----------------------------------------------------------------------
! 2. Read and amend fixed header and lookups of previous file
!----------------------------------------------------------------------

! 2.0 read headers

      LPrevious = .True.
! DEPENDS ON: read_one_header
      call read_one_header(UnitPrevious, icode,                         &
     &     Len_FixHd, Len1_Lookup, Len2_LookupPrevious,                 &
     &     Len2_ActualPrevious, FixHdPrevious,                          &
#include "argppx.h"
     &     LookupPrevious)

      if (icode  /=  0) then
        LPrevious = .False.
        write(UnWarn,*)CWarn,CSub,                                      &
     &  ' step 2.0 failed to open and read headers of' //               &
     &  ' previous flux file '
        icode = 0
      end if

! 2.1 amend headers

      if ( LPrevious ) then

#if defined(FLXPLPR)
! 2.1.0 set LCyclic (T if atmosphere grid has wrap points)
      if ( MOD ( FixHdPrevious (4) , 100 )  /=  3 ) then
        LCyclic = .True.
      else
        LCyclic = .False.
      end if

! 2.1.1 attempt to set field dimensions
! DEPENDS ON: read_field_sizes
        call read_field_sizes(                                          &
     &      Len1_Lookup, Len2_ActualPrevious, LookupPrevious,           &
     &      ncols, nrowst, nrowsu, icode)
        if ( icode  >   0) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 2.1 read_field_sizes failed for previous file'
          go to 9999
        end if

! 2.1.2 Set dimensions for atmosphere velocity grids.
!       Calculations differ for B and C grids.
      if ( l_B_grid) then
        nrowsv = nrowsu
        nrowsuv = nrowst - 1
      else
        nrowsv = nrowsu - 1
        nrowsuv = nrowst - 1
      end if

#endif

! 2.1.2 amend headers

        l_climate_field = .false.

! DEPENDS ON: add_lookups
        call add_lookups (                                              &
     &   NoAddTimesPrevious, ISrchOffHrPrevious, INewOffHrPrevious,     &
     &   l_climate_field, Len1_Lookup, Len2_LookupPrevious,             &
     &   Len2_ActualPrevious,                                           &
     &   LookupPrevious, LookFldNoPrevious, icode )

        if ( icode  >   0) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 2.1 add_lookups failed for Previous file'
          icode = icode + 2000
          go to 9999
        end if

      end if   ! LPrevious

!----------------------------------------------------------------------
! 3. Read and amend fixed header and lookups of climate file
!----------------------------------------------------------------------

! 3.0 read headers

      LClimate = .True.
! DEPENDS ON: read_one_header
      call read_one_header(UnitClimate, icode,                          &
     &   Len_FixHd, Len1_Lookup, Len2_LookupClimate,                    &
     &   Len2_ActualClimate, FixHdClimate,                              &
#include "argppx.h"
     &   LookupClimate)

      if (icode  >   0) then
        LClimate = .false.
        write(UnWarn,*)CWarn,CSub,                                      &
     &  ' step 3.0 failed to read headers of climate flux file '
        icode = 0
      end if

! 3.1 amend headers

      if ( LClimate ) then

#if defined(FLXPLPR)
! 3.1.0 set LCyclic (T if atmosphere grid has wrap points)
      if ( MOD ( FixHdClimate (4) , 100 )  /=  3 ) then
        LCyclic = .True.
      else
        LCyclic = .False.
      end if
! 3.1.1 attempt to set field dimensions
! DEPENDS ON: read_field_sizes
        call read_field_sizes(                                          &
     &      Len1_Lookup, Len2_ActualClimate, LookupClimate,             &
     &      ncols, nrowst, nrowsu, icode)
        if ( icode  >   0) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 2.2 read_field_sizes failed for climate file'
          go to 9999
        end if

! 3.1.2 Set dimensions for atmosphere velocity grids.
!       Calculations differ for B and C grids.
      if ( l_B_grid) then
        nrowsv = nrowsu
        nrowsuv = nrowst - 1
      else
        nrowsv = nrowsu - 1
        nrowsuv = nrowst - 1
      end if

#endif

! 3.1.2 amend headers

        l_climate_field = .true.

! DEPENDS ON: add_lookups
        call add_lookups (                                              &
     &   NoAddTimesClimate, ISrchOffHrClimate, INewOffHrClimate,        &
     &   l_climate_field, Len1_Lookup, Len2_LookupClimate,              &
     &   Len2_ActualClimate,                                            &
     &   LookupClimate, LookFldNoClimate, icode )

        if ( icode  >   0) then
          write(UnErr,*)CErr,CSub,                                      &
     &    ' step 3.1 add_lookups failed for Climate file'
          icode = icode + 2500
          go to 9999
        end if

      end if   ! LClimate

!----------------------------------------------------------------------
! 4. If no file headers have been read exit with a fatal error
!----------------------------------------------------------------------
      if ( .not. ( LPreferred .or. LPrevious .or. LClimate) ) then
        icode = 16
         write(UnErr,*)CErr,CSub,                                       &
     &  ' step 4. failed to read headers of any flux file '
        go to 9999
      end if

9999  continue
      return
      END SUBROUTINE read_field_headers
!----------------------------------------------------------------------
#if defined(FLXPLPR)
#endif
!----------------------------------------------------------------------
#endif
