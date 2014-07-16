#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL
!LL    Subroutine:
!LL    ACUMPS
!LL
!LL    Purpose:
!LL    To accumulate partial sums of climate mean tagged diagnostics
!LL    and create dumps containing them. Also to overwrite the D1
!LL    diagnostice with the partial sum for use by MEANPS. This
!LL    saves MEANPS having to reread the partial sum dump.
!LL
!LL    Tested under compiler:
!LL    cft77
!LL
!LL    Tested under OS version:
!LL    UNICOS 5.1
!LL
!LL AD, DR      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  19/02/93  Use FIXHD(12) not FIXHD(1) as Version no in P21BITS
!LL   3.1   25/01/93 : Correct LBPACK for 32 bit dumps after changes.
!LL 3.4  16/6/94 : Change CHARACTER*(*) to CHARACTER*(80) N.Farnon
!     4.1  18/06/96  Changes to cope with changes in STASH addressing
!                    Author D.M. Goddard.
!LL   4.2  27/11/96  MPP changes for T3E.  Using READDUMP and
!LL                  WRITDUMP with partial mean files. K Rogers
!LL   4.3  22/01/97  Use MPP_LOOKUP to address D1 on MPP
!LL                  S.D.Mullerworth
!LL   4.3  10/04/97  Call READDUMP without reading header to avoid
!LL                  overwriting REALHD with previous values. K Rogers
!LL   4.4  22/09/97  Remove superfluous arrays. Add extra error trap
!LL                  S.D. Mullerworth
!LL   4.4  26/06/97  Changes to allow climate means with Gregorian
!LL                  calendar. Author: M. Gallani
!LL   4.4  16/06/97  Add Broadcast after the WRITDUMP, so
!LL                  that all the processors know the answer
!LL                    Author: Bob Carruthers, Cray Rsearch.
!LL   4.5  23/10/98  Remove unused arrays. S.D.Mullerworth
!LL   5.0  3/6/99    Remove ICODE and CMESSAGE from UM_READDUMP and
!LL                  UM_WRITDUMP.
!LL                  Change order of UM_WRITDUMP arguments
!LL                                                        P.Burton
!LL   5.1  24/01/00  Optimised to sum only tagged diagnostics.
!LL                  S.D.Mullerworth
!LL   5.1  22/02/00  Move PARVARS for TYPSIZE                P.Burton
!LL   5.2  24/07/00  Add new argument to GENERAL_GATHER/SCATTER_FIELD
!LL                                                          P.Burton
!LL   5.2  06/09/00  Adding checksum to spot any corruption of data
!LL                  during I/O operation.                 E.Leung
!LL   5.3  06/06/01  Put T3E CPP switch around code to check for
!LL                  alignment of WORK array.           P.Burton
!     5.5  05/02/03  Portability changes allowing for big_endian
!                    I/O on little_endian platforms.        P.Dando
!    6.1  06/08/04  Avoid out-of-bounds reference. B Ingleby
!LL
!LL    Programming standard:
!LL    UM Doc Paper 3
!LL
!LL    Logical system components covered:
!LL    C5
!LL
!LL    Project tasks:
!LL    C5,C51,C52
!LL
!LL    External documentation:
!LL    On-line UM document C5 - Control of means calculations
!LL
!*L    Interface and arguments:
      SUBROUTINE ACUMPS(                                                &
     &  N_OBJS_D1,D1_ADDR                                               &
     &  ,LEN_DATA,D1,LD1,ID1                                            &
     &  ,MAXSIZE,MEANS_TOTAL                                            &
     &  ,FLAG,NFTIN,NFTOUT,LCLIMREALYR,MEANLEV                          &
     &  ,I_MONTH,I_YEAR                                                 &
     &  ,HEAD_OUT,HEAD_LEN,HEAD_SIZE,                                   &
     &  TIMESTEP,CMITEMS,FIXHD12,                                       &
#include "argsts.h"
     &  ICODE,CMESSAGE)
!
      IMPLICIT NONE
!
#include "d1_addr.h"
      INTEGER                                                           &
     &  N_OBJS_D1               !IN No objects in D1 array

      INTEGER                                                           &
     &  D1_ADDR(D1_LIST_LEN,N_OBJS_D1) !IN Addressing of D1 array

      INTEGER                                                           &
     &  MAXSIZE,                                                        &
                                ! IN dimension of largest data field
     &  LEN_DATA,                                                       &
                                ! IN Length of model data
     &  FLAG,                                                           &
                                ! IN Flag for reading partial sum dump
     &  NFTIN,                                                          &
                                ! IN Unit no for reading partial sums
     &  NFTOUT,                                                         &
                                ! IN Unit no for writing partial sums
     &  ICODE,                                                          &
                                ! OUT Return code; successful=0
                                !                  error>0
     &  MEANLEV,                                                        &
                                ! IN level of climate meaning
     &  MEANS_TOTAL,                                                    &
                                ! IN Indicates a meaning period
     &  I_MONTH,                                                        &
                                ! IN Current model time (months)
     &  I_YEAR,                                                         &
                                ! IN Current model time (years)
     &  FIXHD12                                                         &
                                ! IN Version of model
     &  ,CMITEMS                                                        &
                                ! IN Number of items being meaned
     &  ,TIMESTEP               ! IN Submodel timestep
!
      CHARACTER *(80)                                                   &
     &       CMESSAGE             ! OUT Error message if ICODE>0
!
      INTEGER                                                           &
     &  ID1(LEN_DATA)           ! IN/OUT Integer equiv. of data block.
                                !        Overwritten with partial sums.
!
      REAL                                                              &
     &  D1(LEN_DATA)            ! IN/OUT Real equivalence of data block
!
      LOGICAL                                                           &
     &  LD1(LEN_DATA),                                                  &
                                ! IN/OUT Logical equiv. of data block
     &  LCLIMREALYR             ! IN Real-period climate meaning
!
!      Common blocks
!
#include "clookadd.h"
#include "c_mdi.h"
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "typsts.h"
#include "stparam.h"
      INTEGER                                                           &
     &  HEAD_LEN                                                        &
     &  ,HEAD_SIZE

      INTEGER                                                           &
     &  HEAD_OUT(HEAD_LEN,TOTITEMS)                                     &
                                    ! IN Header contains packing
                                !    info for output ps file
     &  ,HEAD_BUF(HEAD_SIZE)

! Header formatted as follows:
! HEAD_OUT(1,*): No of words per level in field
! HEAD_OUT(2,*): 2 for packed, 1 for unpacked
! HEAD_OUT(3,*): No of words per level on disk

! Align for well-formed io
!dir$ cache_align head_buf

!
!*L
!*L    External subroutines called:
!
!      Cray specific functions  UNIT,LENGTH
!
!      Local variables
!
      INTEGER                                                           &
     &  I,J,K                                                           &
                                ! Loop indices
     &  ,LEN_IO                                                         &
                                ! Actual IO length
     &  ,CITEMS                                                         &
                                ! Count variable
     &  ,PERIODLEN                                                      &
                                ! Current meaning period in days
     &  ,TAG                                                            &
                                ! Stash tag
     &  ,PTD1                                                           &
                                ! Pointer to D1_ADDR information
     &  ,address                                                        &
                                ! Address in local D1
     &  ,levels                                                         &
                                ! Number of levels per diagnostic
     &  ,length                                                         &
                                ! Length of each level in local D1
     &  ,global_length                                                  &
                                ! Length of global field
     &  ,offset                 ! Indexing offset for WORK array
!
      INTEGER                                                           &
     &  HEADER(2)                                                       &
                                ! Initial header
     &  ,HEAD_IN(HEAD_LEN,CMITEMS) ! Packing info for input ps file
                                ! Will differ from HEAD_OUT if packing
                                ! codes have changed mid-run

      REAL                                                              &
     &  IOSTAT,                                                         &
                                ! IO error code
     &  REALPERIODLEN,                                                  &
                                ! explicitly real equivalent
                                ! of PERIODLEN
     &  AWORK,                                                          &
                                ! Accumulative sum for WORK array
     &  CWORK,                                                          &
                                ! Accumulative sum for WORK array
     &  CKWORK,                                                         &
                                ! Checksum for WORK array
     &  CKWORKO                 ! Packed CKWORK
!
!      Local arrays
!
      REAL                                                              &
     &  D1_DATA(MAXSIZE)        ! Work area for fields
      REAL                                                              &
     &  WORK(MAXSIZE+4)        ! Work area and IO buffer

! Align for well-formed io
!dir$ cache_align work

!
      IF(ICODE /= 0)GOTO 999

! Arrays sent to BUFFIN/BUFFOUT need to be cache aligned for
! well-formed io to work, but the cache_align directive wasn't working
! correctly for WORK. The following adds an offset to the index of
! WORK which resolves the problem.
! This is an offset from 0, so an offset of 1 really means no offset!
      offset=1
#if defined(T3E)

      i=LOC(WORK)
      IF (MOD(i,32) /= 0)THEN
! Buffers must start on 32-byte boundary. If the WORK array does not
! start on a 32-byte boundary, then add an offset to the index so that
! the first element of WORK passed to the buffer routines *is* on a
! boundary
        offset=1+((i/32)*32+32-i)/8
      ENDIF
#endif

!   Set up variables needed for weighting accumulations if real-period
!   climate meaning is selected. Partial sums are normalised elsewhere.

      if (lclimrealyr) then
! DEPENDS ON: setperlen
        call setperlen(meanlev,i_month,i_year,periodlen)
        realperiodlen=real(periodlen)
      endif

! STEP 1: Read in headers of previous partial sum and write out
!         header of new.

      IF (FLAG /= 1) THEN       ! PS data exist on disk
! Read headers for input partial sum file
! DEPENDS ON: buffin
        CALL BUFFIN(NFTIN,HEAD_BUF,HEAD_SIZE,LEN_IO,IOSTAT)
        IF(IOSTAT /= -1.0.OR.LEN_IO /= HEAD_SIZE)THEN
          WRITE(6,*)'ACUMPS: Error reading header: IO code ',           &
     &      IOSTAT,' on unit ',NFTIN
          WRITE(6,*)'Words requested ',HEAD_SIZE,                       &
     &      ' Words read ',LEN_IO
          ICODE=1
          CMESSAGE='ACUMPS: BUFFIN error - see output'
          GOTO 999
        ENDIF
! Transfer header information from buffer to header arrays
        HEADER(1)=HEAD_BUF(1) ! Timestep of creation
        HEADER(2)=HEAD_BUF(2) ! Number of records
        K=3
        DO I=1,CMITEMS
          DO J=1,HEAD_LEN
            HEAD_IN(J,I)=HEAD_BUF(K)
            K=K+1
          ENDDO
        ENDDO

        IF (HEADER(1) >= TIMESTEP.OR.HEADER(2) /= CMITEMS)THEN
          WRITE(6,*)'ACUMPS1: Partial sum file inconsistent'
          WRITE(6,*)'PS file holds ',HEADER(2),' items and written      &
     &      at STEP ',HEADER(1)
          WRITE(6,*)'Expected timestep should be < ',TIMESTEP
          WRITE(6,*)'Expected number of items ',CMITEMS
          CMESSAGE='ACUMPS1: Partial sum file inconsistent. See Output'
          ICODE=2
          GOTO 999
        ENDIF
      ELSE
! No input sum, so initialise header array
        DO I=1,HEAD_SIZE
          HEAD_BUF(I)=0
        ENDDO
      ENDIF

! Write headers for new partial sum file
! Transfer information to io buffer
      HEAD_BUF(1)=TIMESTEP
      HEAD_BUF(2)=CMITEMS
      K=3
      DO I=1,CMITEMS
        DO J=1,HEAD_LEN
          HEAD_BUF(K)=HEAD_OUT(J,I)
          K=K+1
        ENDDO
      ENDDO
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,HEAD_BUF,HEAD_SIZE,LEN_IO,IOSTAT)
      IF(IOSTAT /= -1.0.OR.LEN_IO /= HEAD_SIZE)THEN
        WRITE(6,*)'ACUMPS: Error writing header: IO code ',             &
     &    IOSTAT,' on unit ',NFTOUT
        WRITE(6,*)'Words requested ',HEAD_SIZE,                         &
     &    ' Words written ',LEN_IO
        ICODE=4
        CMESSAGE='ACUMPS: BUFFOUT error - see output'
        GOTO 999
      ENDIF

! STEP 2 : Loop over all STASH items. For each tagged item, gather
!          current data to D1_DATA array, read partial sum into WORK
!          array (if there is a partial sum), sum the two and write
!          out to new partial sum file.
!           Also, if this is a meaning period, overwrite the field
!          in D1 with the complete sum, to be picked up by MEANPS.

!     Start of loop over STASH items
      CITEMS=0
      DO K=1,TOTITEMS
        TAG=STLIST(st_macrotag,K)/1000
        PTD1=STLIST(st_d1pos,K)
        IF(TAG/=0)THEN
          IF(STLIST(s_modl,k)==D1_ADDR(d1_imodl,PTD1))THEN
! Object tagged for climate meaning and in relevant internal model
          address=D1_ADDR(d1_address,PTD1)
          levels=D1_ADDR(d1_no_levels,PTD1)
          length=D1_ADDR(d1_length,PTD1)/levels
#if defined(MPP)
          global_length=STLIST(st_dump_level_output_length,K)
#else
          global_length=STLIST(st_output_length,K)/levels
#endif
          CITEMS=CITEMS+1
          DO J=1,levels
! Copy current field from D1 to D1_DATA
#if defined(MPP)
! by gathering full field to pe0
! DEPENDS ON: general_gather_field
            CALL GENERAL_GATHER_FIELD(                                  &
     &        D1(address),D1_DATA,length,                               &
     &        global_length,1,                                          &
     &        D1_ADDR(1,PTD1),0,                                        &
     &        ICODE,CMESSAGE)
            IF(ICODE /= 0)GOTO 999
#else
            DO I=1,global_length
              D1_DATA(I)=D1(address+I-1)
            ENDDO
#endif
            DO I=global_length+1,MAXSIZE
              D1_DATA(I)=0.
            ENDDO
! Set initial value for AWORK and CWORK
            AWORK=0.0
            CWORK=0.0
! If partial sum exists on disk, read it in and add to current field
            IF (FLAG /= 1) THEN ! PS data exist on disk
! Read in one level of partial sum field
              IF (HEAD_IN(2,CITEMS)  ==  2) THEN
! Data is packed using CRAY 32 bit method - note that we need to read
! in 2*HEAD_IN(3,CITEMS) 32 bit words using BUFFIN32
! DEPENDS ON: buffin32
                CALL BUFFIN32(NFTIN,WORK(offset),2*HEAD_IN(3,CITEMS)    &
     &            ,LEN_IO,IOSTAT)
! And then halve LEN_IO to satisfy tests against HEAD_IN(3,CITEMS)
                LEN_IO = LEN_IO/2
              ELSE
! For non-packed data
! DEPENDS ON: buffin
                CALL BUFFIN(NFTIN,WORK(offset),HEAD_IN(3,CITEMS)        &
     &            ,LEN_IO,IOSTAT)
              ENDIF
              IF(IOSTAT /= -1.0.OR.LEN_IO /= HEAD_IN(3,CITEMS))THEN
                WRITE(6,*)'ACUMPS: Error reading partial sum IO code ', &
     &            IOSTAT,' on unit ',NFTIN
                WRITE(6,*)'Words requested ',HEAD_IN(3,CITEMS),         &
     &            ' Words read ',LEN_IO
                ICODE=6
                CMESSAGE='ACUMPS: BUFFIN error - see output'
                GOTO 999
              ENDIF
#if defined(MPP)
              IF (mype == 0) THEN
! Valid data exists on pe0 only
#endif
! Unpack if data on disk was packed
                IF (HEAD_IN(2,CITEMS) == 2)THEN
! DEPENDS ON: expand32b
                  CALL EXPAND32B(GLOBAL_LENGTH+1,WORK(offset),FIXHD12)
! Calculate a checksum
                  DO I=1,GLOBAL_LENGTH
                    AWORK=AWORK+WORK(I+offset-1)
                  END DO
                  CKWORK=AWORK/INT(global_length)
! Pack and umpack checksum to force it losing precision in order to do
! the comparison
! DEPENDS ON: pack21
                  CALL PACK21(1,CKWORK,CKWORKO)
! DEPENDS ON: expand32b
                  CALL EXPAND32B(1,CKWORKO,FIXHD12)
                ELSE
                  DO I=1,GLOBAL_LENGTH
                    AWORK=AWORK+WORK(I+offset-1)
                  END DO
                  CKWORKO=AWORK/INT(global_length)
                ENDIF
                IF(CKWORKO /= WORK(global_length+offset)) THEN
                  WRITE(6,*)'WARNING: checksum detects a corruption'
                  ICODE=4
                  CMESSAGE='ACUMPS: Data corruption during I/O'
                  GOTO 999
                ENDIF
! Sum with field in D1 - Scale data if 365 day calendar
                IF (LCLIMREALYR)THEN
                  DO I=1,global_length
                    IF (WORK(I+offset-1) == RMDI)THEN
                      D1_DATA(I)=RMDI
                    ELSE
                      D1_DATA(I)=WORK(I+offset-1)+                      &
     &                  (realperiodlen*D1_DATA(I))
                    ENDIF
                  END DO
                ELSE
! 360 day calendar
                  DO I=1,global_length
                    IF (WORK(I+offset-1) == RMDI)THEN
                      D1_DATA(I)=RMDI
                    ELSE
                      D1_DATA(I)=WORK(I+offset-1)+D1_DATA(I)
                    ENDIF
                  END DO
                ENDIF
#if defined(MPP)
              ENDIF
#endif
            ELSE
! First data for this period - no partial sum to add
              IF (LCLIMREALYR)THEN
! Scale initial data if 365 day calendar
#if defined(MPP)
                IF (mype == 0) THEN
#endif
                  DO I=1,global_length
                    IF (D1_DATA(I) /= RMDI)THEN
                      D1_DATA(I)=realperiodlen*D1_DATA(I)
                    ENDIF
                  END DO
#if defined(MPP)
                ENDIF
#endif
              ENDIF
            ENDIF ! End of adding PS data

!         Write out sum to PS file
#if defined(MPP)
            IF (mype == 0) THEN
#endif
! Copy data to WORK array, packing if necessary
              IF (HEAD_OUT(2,CITEMS) == 2)THEN
                DO I=HEAD_OUT(1,CITEMS),MAXSIZE
                  WORK(I+offset-1)=0.
                ENDDO
! DEPENDS ON: pack21
                CALL PACK21(GLOBAL_LENGTH+1,D1_DATA,                    &
                            WORK(offset))
! DEPENDS ON: expand32b
                CALL EXPAND32B(GLOBAL_LENGTH+1,WORK(offset),FIXHD12)
                DO I=1,GLOBAL_LENGTH
                  CWORK=CWORK+WORK(I+offset-1)
                END DO
                WORK(global_length+offset)=CWORK/INT(global_length)
! DEPENDS ON: pack21
                CALL PACK21(GLOBAL_LENGTH+1,WORK(offset),               &
                            WORK(offset) )

! If data not packed, calculate checksum straight away
              ELSE
                DO I=1,GLOBAL_LENGTH
                  WORK(I+offset-1)=D1_DATA(I)
                  CWORK=CWORK+WORK(I+offset-1)
                ENDDO
                WORK(global_length+offset)=CWORK/INT(global_length)
              ENDIF
#if defined(MPP)
            ENDIF
#endif


! Output partial sum to file
            IF (HEAD_OUT(2,CITEMS)  ==  2) THEN
! Data is packed using CRAY 32 bit method - note that we need to write
! out 2*HEAD_OUT(3,CITEMS) 32 bit words using BUFFO32
! DEPENDS ON: buffo32
                CALL BUFFO32(NFTOUT,WORK(offset),                       &
     &            2*HEAD_OUT(3,CITEMS),LEN_IO,IOSTAT)
! And then halve LEN_IO to satisfy tests against HEAD_OUT(3,CITEMS)
                LEN_IO = LEN_IO/2
            ELSE
! For non-packed data
! DEPENDS ON: buffout
              CALL BUFFOUT(NFTOUT,WORK(offset),                         &
     &          HEAD_OUT(3,CITEMS),LEN_IO,IOSTAT)
            ENDIF
            IF(IOSTAT /= -1.0.OR.LEN_IO /= HEAD_OUT(3,CITEMS))THEN
              WRITE(6,*)'ACUMPS: Error writing partial sum. Code ',     &
     &          IOSTAT,' on unit ',NFTOUT
              WRITE(6,*)'Words requested ',HEAD_OUT(3,CITEMS),          &
     &          ' Words written ',LEN_IO
              ICODE=7
              CMESSAGE='ACUMPS: BUFFOUT error - see output'
              GOTO 999
            ENDIF
            IF (MEANS_TOTAL /= 0)THEN
! Overwrite field in D1 with partial sum for use by MEANPS
#if defined(MPP)
              IF (mype == 0)then
! Pack and unpack for bit comparison with old system
                IF (HEAD_OUT(2,CITEMS) == 2)THEN
                  DO I=1,HEAD_OUT(1,CITEMS)
                    D1_DATA(I)=WORK(I+offset-1)
                  ENDDO
! DEPENDS ON: expand32b
                  CALL EXPAND32B(GLOBAL_LENGTH,D1_DATA,FIXHD12)
                ENDIF
              ENDIF
! DEPENDS ON: general_scatter_field
              CALL GENERAL_SCATTER_FIELD(                               &
     &         D1(address),D1_DATA,LENGTH,global_length,1,              &
     &          D1_ADDR(1,PTD1),0,ICODE,CMESSAGE)
              IF(ICODE /= 0)GOTO 999
#else
              DO I=1,global_length
                D1(address+I-1)=D1_DATA(I)
              ENDDO
#endif
            ENDIF
            address=address + length ! Point to next level
          ENDDO                 ! End loop over levels
          ENDIF
        ENDIF                   ! End tagged for meaning
      END DO                    ! End of loop over STASH list

 999  CONTINUE
      RETURN
      END SUBROUTINE ACUMPS
#endif
