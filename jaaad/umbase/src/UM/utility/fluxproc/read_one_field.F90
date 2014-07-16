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
! 6.1      30/07/04     Add call to set decomposition variables.
!                       A. Hines.
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! contains routines: read_one_field
!
! Purpose: Flux processing routine.
!          Selects a field from input data files and returns
!          data (in array field) and lookup table in PP_Int and PP_Real
!
! Uses:    readflds and coex
!
! also transfers bmks scale and offset from pp header to data field
!      outputs debugging information on field read in
!----------------------------------------------------------------------
      subroutine read_one_field ( InUnit, item, itemvalue,              &
#include "avaltim.h"
     &               Len_FixHd, FixHd, Len1_Lookup,                     &
     &               Len2_Lookup, Lookup, LookFldNo,                    &
     &               ldebug, l_climate_field,                           &
     &               Len_IntHd, Len_RealHd, IntHead, RealHead,          &
     &               ncols, nrows, field,                               &
#include "argppx.h"
     &               icode)

      implicit none

! declaration of parameters
#include "clookadd.h"
#include "c_mdi.h"

! Local parameters:
      integer len_full_word   ! The length of a FULL_WORD
      parameter ( len_full_word = 64 )

! declaration of argument list

! search conditions (all intent IN)
      integer InUnit    ! IN    unit number for input
      integer item      ! IN    lookup header item to test
      integer itemvalue ! IN    value to look for

! validity time to look for: intent IN
#include "cvaltim.h"

! fixed headers and lookup tables to use: all intent IN
      integer Len_FixHd                    ! length of fixed header
      integer FixHd(Len_FixHd)             ! fixed header
      integer Len1_Lookup, Len2_Lookup     ! true lengths of tables
      integer Lookup(Len1_Lookup, Len2_Lookup)  ! lookup tables
      integer LookFldNo(Len2_Lookup)       ! field nos for lookups

! control logical for debugging output
      logical ldebug          ! IN T => output debugging info
      logical l_climate_field ! IN T => trying to read climate field
                              !    F => trying to read NWP field

! lengths of Lookup table
      integer Len_IntHd         ! IN   length of integer part of lookup
      integer Len_RealHd        ! IN   length of real part of lookup

! lookup tables of field found
      integer IntHead(Len_IntHd) ! OUT integer part
      real RealHead(Len_RealHd)
! OUT real part

! output field
      integer ncols             ! IN  number of columns
      integer nrows             ! IN  number of rows
      real field(ncols,nrows)   ! OUT field values
      integer icode  ! IN/OUT error code ; > 0 => fatal error detected

! declaration of globals used
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cunitnos.h"
#include "cmess.h"
#include "cdebug.h"

! declaration of local arrays
      integer field_wgdos_packed(ncols*nrows)  ! wgdos packed field

! declaration of local scalars
      integer fld_no      ! number of field matching search conditions
      integer i           ! do loop index
      integer nvalues     ! # values in wgdos packed field
      integer idum        ! dummy integer
      integer ixx         ! # of columns in grid according to coex
      integer iyy         ! # of rows in grid according to coex
      integer gl_row_len  ! # of columns in global grid
      integer gl_n_rows   ! # of rows in global grid
      logical ItemFound   ! T => item to search for has been found
      logical l_data_time ! T => use data time; F => use validity time

      real offset         ! offset of data from zero
      real scale          ! MKS scaling factor
      real rmdit          ! rmdi read in from file

      character *256 cmessage  ! error message

      external time_to_use, readflds, copy_to_real, scalarmult,         &
     &         scalaradd,output_debug

!----------------------------------------------------------------------

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'read_one_field'  ! subroutine name for error messages
      idum = 0                 ! dummy integer used in call coex

! 1. Decide whether to search lookups using validity or data time

! DEPENDS ON: time_to_use
      call time_to_use ( itemvalue, l_climate_field, l_data_time)

! 2. Search lookup tables for match on date and stash item

! 2.1 search lookup tables using data time

      ItemFound = .false.

      if ( l_data_time ) then

        do i = 1, Len2_Lookup
          if ( Lookup(LBYRD,i)   ==  ValidYear   .and.                  &
     &       Lookup(LBMOND,i)  ==  ValidMonth  .and.                    &
     &       Lookup(LBDATD,i)  ==  ValidDay    .and.                    &
     &       Lookup(LBHRD,i)   ==  ValidHour   .and.                    &
     &       Lookup(LBMIND,i)  ==  ValidMin    .and.                    &
     &       Lookup(item,i)   ==  itemvalue ) then
            ItemFound = .true.
            fld_no = LookFldNo(i)
            go to 100
          end if
        end do

! 2.2 Search lookup tables using validity time
      else   !  .not. l_data_time

        do i = 1, Len2_Lookup
          if ( Lookup(LBYR,i)   ==  ValidYear   .and.                   &
     &       Lookup(LBMON,i)  ==  ValidMonth  .and.                     &
     &       Lookup(LBDAT,i)  ==  ValidDay    .and.                     &
     &       Lookup(LBHR,i)   ==  ValidHour   .and.                     &
     &       Lookup(LBMIN,i)  ==  ValidMin    .and.                     &
     &       Lookup(item,i)   ==  itemvalue ) then
            ItemFound = .true.
            fld_no = LookFldNo(i)
            go to 100
          end if
        end do

      endif   !  l_data_time

100   continue

! if item has not been found set icode > 0 and exit routine

      if ( .not. ItemFound ) then
        icode = 36
        write(UnWarn,*)CWarn,CSub,                                      &
     &       ' step 2. unable to find required field '
        go to 9999
      end if


! 3. check that nrows and ncols agree with those in lookup table

      if ( Lookup(LBNPT,fld_no)  /=  ncols  ) then
        icode = 37
        write(UnWarn,*)CWarn,CSub,'3.1 number of columns do ',          &
     &  'not agree: ' ,ncols, Lookup(LBNPT,fld_no)
        go to 9999
       end if

      if ( Lookup(LBROW,fld_no)  /=  nrows  ) then
        icode = 38
        write(UnWarn,*)CWarn,CSub,'3.2 number of rows do ',             &
     &  'not agree: ' ,nrows, Lookup(LBROW,fld_no)
        go to 9999
      end if

! 4. If found: Use READFIELDS to extract field


        gl_row_len=lookup(19,fld_no)
        gl_n_rows=lookup(18,fld_no)

! DEPENDS ON: decompose_smexe
        call decompose_smexe(gl_row_len, gl_n_rows,0,0,1)

! 4.1 extract field which is wgdos packed

      if ( MOD ( Lookup(LBPACK,fld_no) ,10)  ==  1) then

        nvalues = Lookup(LBLREC,fld_no)

! DEPENDS ON: readflds
        call readflds (InUnit , 1, fld_no, LOOKUP,                      &
     &      Len1_Lookup, field_wgdos_packed, nvalues, FIXHD,            &
#include "argppx.h"
     &      icode, cmessage)

        if ( icode  >   0 ) then
          write(UnWarn,*)CWarn,CSub,                                    &
     &       ' step 4.1 unable to read field: cmessage is ',            &
     &       cmessage
          icode = 39
          go to 9999
        end if

! DEPENDS ON: coex
        call coex(field,                                                &
                                       ! OUT unpacked field
     &            nrows*ncols,                                          &
                                       ! IN  size of unpacked field
     &            field_wgdos_packed,                                   &
                                       ! IN  packed field
     &            nvalues,                                              &
                                       ! IN  size of packed field
     &            ixx,iyy,                                              &
                                       ! OUT row and column sizes
     &            idum,idum,                                            &
                                       ! IN  not used
     &            .false.,                                              &
                                       ! IN  => expansion
     &            rmdi,                                                 &
                                       ! IN  real missing data value
     &            len_full_word,                                        &
                                       ! IN  length of a full word
     &            icode,                                                &
     &            cmessage)

        if ( ixx  /=  ncols  .or. iyy  /=  nrows ) then
          icode = 40
          write(UnWarn,*)CWarn,CSub,                                    &
     &       ' step 4.1 number of rows and columns garbled ?  ',        &
     &       ixx, ncols, iyy, nrows
          go to 9999
        end if

! 4.2 or extract field which is not packed

      else

        if ( Lookup(LBLREC,fld_no)  /=  ncols*nrows) then
           icode = 41
           write(UnWarn,*)CWarn,CSub,                                   &
     &     ' step 4.2 wrong number of data points in field ',           &
     &     Lookup(LBLREC,fld_no), ncols*nrows, ncols, nrows
           go to 9999
        end if

! DEPENDS ON: readflds
        call readflds (InUnit , 1, fld_no, LOOKUP,                      &
     &      Len1_Lookup, Field, ncols*nrows, FIXHD,                     &
#include "argppx.h"
     &      icode, cmessage)

        if ( icode  >   0 ) then
          write(UnWarn,*)CWarn,CSub,                                    &
     &       ' step 4.2 unable to read field: cmessage is ',            &
     &       cmessage
          icode = 42
          go to 9999
        end if

      end if ! Lookup(LBPACK,fld_no)

! 5.  convert lookup table to Int_Head and Real_Head
!     and field to a 2D field

      do i = 1, Len_IntHd
        IntHead(i) = Lookup(i,fld_no)
      end do

      do i = Len_IntHd+1, Len_IntHd+Len_RealHd
! DEPENDS ON: copy_to_real
        call copy_to_real( Lookup(i,fld_no),                            &
     &                        RealHead(i-Len_IntHd) )
      end do

! 6.  correct data offset and change to SI units
      rmdit = RealHead(BMDI   - Len_IntHd)

      if ( rmdit  /=  rmdi ) then
        icode = 43
        write(UnWarn,*)CWarn,CSub,                                      &
     &       ' step 6.1 real missing data indicators do not match: ',   &
     &       rmdit, rmdi
        go to 9999
      end if

      offset = RealHead(BDATUM - Len_IntHd)
      if ( offset  /=  rmdi .and. offset  /=  0.0 ) then
! DEPENDS ON: scalaradd
        call ScalarAdd(ncols, nrows, rmdi, offset,                      &
     &                 Field, Field, icode, cmessage)
        RealHead(BDATUM - Len_IntHd) = rmdi

        write(UnWarn,*)CWarn,CSub,                                      &
     &       ' step 6.2 adding offset factor  ', offset

      end if

      scale = RealHead(BMKS - Len_IntHd)
      if ( scale  /=  rmdi .and. scale  /=  1.0 ) then
! DEPENDS ON: scalarmult
        call ScalarMult(ncols, nrows, rmdi, scale,                      &
     &                 Field, Field, icode, cmessage)
        RealHead(BMKS - Len_IntHd) = rmdi

        write(UnWarn,*)CWarn,CSub,                                      &
     &       ' step 6.3 multiplying by factor  ', scale

      end if

! 7.  output debug info
      if (ldebug) then
        write(OutUnitDbg,*) ' read_data: unit ', InUnit,'; item ',      &
     &    '; itemvalue ', itemvalue
        CMessage = ' '
! DEPENDS ON: output_debug
        call  output_debug(CMessage, nrows, ncols, Field)
      end if

9999  continue
      return
      END SUBROUTINE read_one_field
!----------------------------------------------------------------------
#endif
