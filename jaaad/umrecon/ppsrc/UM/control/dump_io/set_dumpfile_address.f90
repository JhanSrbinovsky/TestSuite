

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
!    V5.3     Enable mpp as the only option for small executables
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

      Use Rcf_HeadAddress_Mod, Only :                                   &
     &    UM_Sector_Size


      Use Rcf_Parvars_Mod, Only :                                       &
     &    mype

      implicit none





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

!*L------------------ COMDECK LOOKADD ----------------------------------
!LL
!LL Purpose : Contains information about the format
!LL           of the PP header
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   4.0  12/09/95   Change NPERIODS to LBUSER3, BRSVD1 to BULEV,
!LL                   BRSVD2 to BHULEV and definitions for BRLEV and
!LL                   BHRLEV. Corresponding changes made to STWORK1A
!LL                   and PPHEAD1A. (Andrew Brady)
!LL  4.0  12/10/95  Change item 45 from lbuser7 to model_code. RTHBarnes
!LL  5.1  17/04/00    Fixed/Free format. P.Selwood.
!LL  5.2  25/09/00    Add LBCC_xxxx variables for the compressed
!LL                   LBC LOOKUP array                  P.Burton
!LL
!LL Programming standard :
!LL
!LL Logical components covered : F092
!LL
!LL Project task :
!LL
!LL External documentation:
!LL
!LLEND -----------------------------------------------------------------
!

! Validity time
      Integer, Parameter :: LBYR   =1   ! Year
      Integer, Parameter :: LBMON  =2   ! Month
      Integer, Parameter :: LBDAT  =3   ! Day of month
      Integer, Parameter :: LBHR   =4   ! Hour
      Integer, Parameter :: LBMIN  =5   ! Minute
      Integer, Parameter :: LBDAY  =6   ! Day number

! Data time
      Integer, Parameter :: LBYRD  =7   ! Year
      Integer, Parameter :: LBMOND =8   ! Month
      Integer, Parameter :: LBDATD =9   ! Day of month
      Integer, Parameter :: LBHRD  =10  ! Hour
      Integer, Parameter :: LBMIND =11  ! Minute
      Integer, Parameter :: LBDAYD =12  ! Day number

      Integer, Parameter :: LBTIM  =13  ! Time indicator
      Integer, Parameter :: LBFT   =14  ! Forcast period (hours)
      Integer, Parameter :: LBLREC =15  ! Length of data record
      Integer, Parameter :: LBCODE =16  ! Grid type code
      Integer, Parameter :: LBHEM  =17  ! Hemisphere indicator
      Integer, Parameter :: LBROW  =18  ! Number of rows in grid
      Integer, Parameter :: LBNPT  =19  ! Number of points per row
      Integer, Parameter :: LBEXT  =20  ! Length of extra data
      Integer, Parameter :: LBPACK =21  ! Packing method indicator
      Integer, Parameter :: LBREL  =22  ! Header release number
      Integer, Parameter :: LBFC   =23  ! Field code
      Integer, Parameter :: LBCFC  =24  ! Second field code
      Integer, Parameter :: LBPROC =25  ! Processing code
      Integer, Parameter :: LBVC   =26  ! Vertical coordinate type
      Integer, Parameter :: LBRVC  =27  ! Coordinate type for reference
                                        ! level

      Integer, Parameter :: LBEXP  =28  ! Experiment number
      Integer, Parameter :: LBEGIN =29  ! Start record
      Integer, Parameter :: LBNREC =30  ! No of records-Direct access
                                        ! only
      Integer, Parameter :: LBPROJ =31  ! Met-O-8 projection number
      Integer, Parameter :: LBTYP  =32  ! Met-O-8 field type
      Integer, Parameter :: LBLEV  =33  ! Met-O-8 level code
      Integer, Parameter :: LBRSVD1=34  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD2=35  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD3=36  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBRSVD4=37  ! Reserved for future PP-package
                                        !  use
      Integer, Parameter :: LBSRCE =38  ! =1111 to indicate following
                                        ! apply to UM
      Integer, Parameter :: DATA_TYPE =39  ! Indicator for real/int or
                                           ! timeseries
      Integer, Parameter :: NADDR  =40  ! Start address in DATA_REAL or
                                        ! DATA_INT
      Integer, Parameter :: LBUSER3=41  ! Free for user-defined function
      Integer, Parameter :: ITEM_CODE =42  !Stash code
      Integer, Parameter :: LBPLEV =43  ! Pseudo-level indicator (if
                                        ! defined)
      Integer, Parameter :: LBUSER6=44  ! Free for user-defined function
      Integer, Parameter :: MODEL_CODE =45 ! internal model identifier

      Integer, Parameter :: BULEV  =46  ! Upper level boundary
      Integer, Parameter :: BHULEV =47  ! Upper level boundary
      Integer, Parameter :: BRSVD3 =48  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BRSVD4 =49  ! Reserved for future PP-package
                                        ! use
      Integer, Parameter :: BDATUM =50  ! Datum value
      Integer, Parameter :: BACC   =51  ! (Packed fields) Packing
                                        ! accuracy
      Integer, Parameter :: BLEV   =52  ! Level
      Integer, Parameter :: BRLEV  =53  ! Lower level boundary
      Integer, Parameter :: BHLEV  =54  ! (Hybrid levels) A-level of
                                        ! value
      Integer, Parameter :: BHRLEV =55  ! Lower level boundary
      Integer, Parameter :: BPLAT  =56  ! Real latitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BPLON  =57  ! Real longitude of 'pseudo'
                                        ! N Pole
      Integer, Parameter :: BGOR   =58  ! Grid orientation
      Integer, Parameter :: BZY    =59  ! Zeroth latitude
      Integer, Parameter :: BDY    =60  ! Latitude interval
      Integer, Parameter :: BZX    =61  ! Zeroth longitude
      Integer, Parameter :: BDX    =62  ! Longitude interval
      Integer, Parameter :: BMDI   =63  ! Missing data indicator
      Integer, Parameter :: BMKS   =64  ! M,K,S scaling factor

      Integer, Parameter :: LBCC_LBYR   = 1  ! Year
      Integer, Parameter :: LBCC_LBMON  = 2  ! Month
      Integer, Parameter :: LBCC_LBDAT  = 3  ! Day of the month
      Integer, Parameter :: LBCC_LBHR   = 4  ! Hour
      Integer, Parameter :: LBCC_LBMIN  = 5  ! Minute
      Integer, Parameter :: LBCC_LBDAY  = 6  ! Day number
      Integer, Parameter :: LBCC_LBEGIN = 7  ! Start record
      Integer, Parameter :: LBCC_NADDR  = 8  ! Start address of DATA
! Mapping of MPP_LOOKUP; analogous to mapping in PP header

      Integer, Parameter :: P_NADDR=1    ! Address on local PE
      Integer, Parameter :: P_LBLREC=2   ! Local length of record

!*----------------------------------------------------------------------
! NADDR IS LOCATION IN PP-HEADER (LOOKUP) FOR START POSN OF VARIABLE
! ITEM_CODE is the location in PP header for a code defined as
!           (section number)*1000+item number
! DATA_TYPE is the location in the PP header defining data as REAL or
!           INTEGER.
! LBNPT is the location defining the number of points per row
!

      integer                                                           &
     & disk_length                                                      &
                                       ! current data length on disk
     &,i                                                                &
                                       ! Loop Index
     &,old_fixhd_160                   ! Original value of fixhd(160)
                                       ! checking as the new addresses
                                       ! are computed

! Including cprintst so that we can use the PrintStatus variable
! CPRINTST defines print status for standard output messages

      ! Minimum output, only essential messages
      INTEGER,PARAMETER :: PrStatus_Min    = 1

      ! Normal informative messages + warnings
      INTEGER,PARAMETER :: PrStatus_Normal = 2

      ! Operational status, all informative messages
      INTEGER,PARAMETER :: PrStatus_Oper   = 3

      ! All informative + extra diagnostic messages
      INTEGER,PARAMETER :: PrStatus_Diag   = 4

      INTEGER PrintStatus ! Control volume of standard output messages
      COMMON/PrintSt/PrintStatus

! CPRINTST end

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
       endif
      endif
      return
      END SUBROUTINE set_dumpfile_address
