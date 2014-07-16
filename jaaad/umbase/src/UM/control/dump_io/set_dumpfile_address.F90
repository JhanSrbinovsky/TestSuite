#if defined(C80_1A) || defined(RECON) || defined(VAROPSVER) ||         \
defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: SDFADR1A ------------------------------------------------
!LL
!LL  Purpose: To set the LBEGIN and LBNREC fields in the LOOKUP Headers
!LL           for VN 16 Type Dumpfiles - addressed by location and
!LL           length which are rounded up the 'UM_SECTOR_SIZE' to make
!LL           Well-Formed I/O Requests.
!LL
!LL  Author:  Bob Carruthers, Cray Research.   Date: 20 May 1997
!LL
!LL  Modifications
!LL  V4.5     Check that there are data in file. C.P. Jones 12 Feb 1998
!    V5.3     Enable MPP as the only option for small executables
!             E.Leung  22 Nov 2001
!    6.2      Removed DIAG92 diagnositc define T. Edwards
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
      SUBROUTINE set_dumpfile_address(fixhd, len_fixhd,                 &
     &                                lookup, len1_lookup, len2_lookup, &
     &                                number_of_data_words_in_memory,   &
     &                                number_of_data_words_on_disk,     &
     &                                disk_address)
!
#if defined(RECON)
      Use Rcf_HeadAddress_Mod, Only :                                   &
     &    UM_Sector_Size
#endif
#if defined(RECON) || defined(VAROPSVER)
      Use Rcf_Parvars_Mod, Only :                                       &
     &    mype
#endif
      implicit none
#if defined(VAROPSVER)
      Integer   UM_Sector_Size
      parameter (UM_Sector_Size = 2048)
#endif

      integer                                                           &
     & len_fixhd                                                        &
                                       ! IN  Length of fixed length
                                       !     header
     &,len1_lookup                                                      &
                                       ! IN  1st dim of lookup
     &,len2_lookup                                                      &
                                       ! IN  2nd dim of lookup
     &,number_of_data_words_in_memory                                   &
                                       ! OUT Number of Data Words
                                       !     in memory
     &,number_of_data_words_on_disk                                     &
                                       ! OUT Number of data words
                                       !     on disk
     &,disk_address                    ! OUT Current rounded disk
                                       !     address and final data
                                       !     length

      integer                                                           &
     & fixhd(len_fixhd)                                                 &
                                       !IN Fixed length header
     &,lookup(len1_lookup,len2_lookup) !IN/OUT PP lookup tables

#include "clookadd.h"
#if !defined(RECON) && !defined(VAROPSVER)
#include "cntl_io.h"
#include "parvars.h"
#endif

      integer                                                           &
     & disk_length                                                      &
                                       ! current data length on disk
     &,i                                                                &
                                       ! Loop Index
     &,old_fixhd_160                   ! Original value of fixhd(160)
                                       ! checking as the new addresses
                                       ! are computed

! Including cprintst so that we can use the PrintStatus variable
#include "cprintst.h"

!
      if(fixhd(160) <  0) RETURN

!--check that the initial data address has been rounded up
!  to a sector boundary - REMEMBER all the code removes
!  one from this address because addresses start at zero.
      if((fixhd(160)-1)  /=                                             &
     & (((fixhd(160)-1)+um_sector_size-1)/um_sector_size)*              &
     & um_sector_size) then
!--save the current initial data address
        old_fixhd_160=fixhd(160)
!--round up the initial disk address
        fixhd(160)=(((fixhd(160)-1)+um_sector_size-1)/                  &
     &   um_sector_size)*um_sector_size+1
        if(mype == 0) then
         if(PrintStatus >= PrStatus_Diag) Then
          write(6,900) old_fixhd_160-1, fixhd(160)-1
900       format(/'SET_DUMPFILE_ADDRESS: Start of Data Address',        &
     &     ' on Disk reset from ',i10,' to ',i10)
#if defined(T3E)
          write(0,900) old_fixhd_160-1, fixhd(160)-1
#endif
         endif
        endif
      endif
!
!--adjust the Dumpfile version Number
!      if(fixhd(1) <  16) fixhd(1)=16
!
!--count the number of words on disk and in memory
      number_of_data_words_on_disk=0
      number_of_data_words_in_memory=0
!
!--find the initial data location on disk
      disk_address=fixhd(160)-1
!
!--loop over all the entries and alter the addresses and lengths
      do i=1, len2_lookup
!--check for a PP type file with an incomplete lookup table
        if(lookup(1, i) == -99) goto 200
!--check for packing to 32-bits
        if(lookup(lbpack,i)-                                            &
     &   ((lookup(lbpack,i)/10)*10) == 2) then
          disk_length=(lookup(lblrec,i)+1)/2
        else
          disk_length=lookup(lblrec,i)
        endif
!--count the number of words
        number_of_data_words_on_disk=                                   &
     &   number_of_data_words_on_disk+disk_length
        number_of_data_words_in_memory=                                 &
     &   number_of_data_words_in_memory+lookup(lblrec,i)
!--round up the length to a number of sectors
        disk_length=((disk_length+um_sector_size-1)/                    &
     &   um_sector_size)*um_sector_size
!--set the disk address
        lookup(lbegin,i)=disk_address
!--set the disk length
        lookup(lbnrec,i)=disk_length
!--increment the disk address
        disk_address=disk_address+lookup(lbnrec,i)
      end do
!--escape for PP type files
200   continue
!

      if(mype == 0) then
       if(PrintStatus >= PrStatus_Diag) Then
!--find the number of bytes in a word
        call word_length(i)
!--print the diagnostic message
        write(6,1000) fixhd(161), number_of_data_words_in_memory,       &
     &   number_of_data_words_on_disk, disk_address-fixhd(160),         &
     &   disk_address, disk_address*i
1000    format(/'SET_DUMPFILE_ADDRESS: Dumpfile LOOKUP Address',        &
     &   ' and Lengths Rewritten:'//                                    &
     &   i10,' Words Stored as Data Length in FIXHD(161)'/              &
     &   i10,' Words Used in Memory for Data'/                          &
     &   i10,' Words Used on Disk for Data'/                            &
     &   i10,' Words Used on Disk for Data after Rounding'/             &
     &   i10,' Words Used on Disk in Total for the File',               &
     &   '  (',i11,' Bytes)'/)
#if defined(T3E)
        write(0,1000) fixhd(161), number_of_data_words_in_memory,       &
     &   number_of_data_words_on_disk, disk_address-fixhd(160),         &
     &   disk_address, disk_address*i
#endif
       endif
      endif
      return
      END SUBROUTINE set_dumpfile_address
#endif
