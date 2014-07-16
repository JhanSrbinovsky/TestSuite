#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL    Subroutine: TRANSOUT -----------------------------------------
!LL
!LL    Purpose:
!LL    To transfer dump data from memory to disk en masse
!LL
!LL    Tested under compiler: cft77
!LL    Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Model
!LL version  Date     Modification history:
!LL   3.3    07/10/93 Corrected order of *CALLs to comdecks TYPSIZE
!LL                   and TYPD1.  Tracey Smith
!LL
!LL    Programming standard:
!LL    UM Doc Paper 3
!LL
!LL    Logical system components covered: C2
!LL
!LL    Project tasks: C2
!LL
!LL    External documentation:
!LL    On-line UM document C5 - Control of means calculations
!LL  Model            Modification history from model version 3.4:
!LL version  Date
!LL 3.4  16/6/94 : Change CHARACTER*(*) to CHARACTER*(80) N.Farnon
!LL   3.5  24/03/95    Changed OPEN to FILE_OPEN and
!LL                    CLOSE to FILE_CLOSE    P.Burton
!LL   4.2  11/10/96    Enable atmos-ocean coupling for MPP.
!LL                    (2): Swap D1 memory.
!LL                    Image of D1 either copied directly from local
!LL                    memory or I/O from file (local to processor)
!LL                    under MPP. I/O kept for mean dumps. R.Rawlins
!LL   4.4  01/07/97    Make transfers to the input file
!LL                    well-formed.
!LL                      Author: Bob Carruthers, Cray Research.
!LL   4.4  11/10/97    Call CNTLALL for L_AO_D1_MEMORY. D. Robinson.
!LL   4.4  28/08/97    Minor tidy: replace SETPOS by SETPOS_SINGLE for
!LL                    MPP case. R.Rawlins
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
!LL
!LLEND---------------------------------------------------------------
!*L    Interface and arguments:
      SUBROUTINE TRANSOUT(                                              &
#include "argd1.h"
     &  LEN_DATA,NFTOUT,sm_ident                                        &
     & ,ICODE,CMESSAGE)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &       LEN_DATA,                                                  &
                                  ! IN Length of model data
     &       NFTOUT,                                                    &
                                  ! IN Unit no for data dump
     &       sm_ident,                                                  &
                                  ! IN submodel identifier
     &       ICODE                ! OUT Return code; successful=0
                                  !                  error>0
!
      CHARACTER*(80)                                                    &
     &       CMESSAGE             ! OUT Error message if ICODE>0
!
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#if defined(T3E)
#include "cntl_io.h"
!
      real local_buffer(um_sector_size)
!dir$ cache_align local_buffer
#endif
!
!      Cray specific functions  UNIT,LENGTH
!
!      External subroutines called
!
      EXTERNAL SETPOS,BUFFOUT
#if defined(MPP)
      EXTERNAL FORT_GET_ENV,OPEN_SINGLE,CLOSE_SINGLE
      EXTERNAL BUFFIN_SINGLE,SETPOS_SINGLE
#else
      EXTERNAL FILE_OPEN,FILE_CLOSE
#endif
!
!      Local variables
!
      INTEGER                                                           &
     &       LEN_IO                                                     &
                                  ! No of 64-bit words buffered in/out
     &      ,I                                                          &
                                  ! loop counter
     &      ,LEN_FILENAME                                               &
                                  ! Length of FILENAME variable
     &      ,LL                                                         &
                                  ! Character length of filename root
     &      ,disk_len_1                                                 &
                                  ! Input length for the first transfer
     &      ,disk_len_2           ! The remainder
!
      REAL                                                              &
     &       A                    ! Error code from UNIT

      LOGICAL                                                           &
     &       D1_COPY_IN_MEMORY    ! T or F: D1 copy in memory or disk

      CHARACTER                                                         &
     &       FILENAME*80          ! File name for copy of D1

#include "mpptrans.h"
#include "cenvir.h"
#include "csmid.h"
#include "chsunits.h"
#include "cntlall.h"

#if defined(MPP)
      D1_COPY_IN_MEMORY=L_AO_D1_MEMORY       ! from COMDECK CNTLALL
      IF(NFTOUT == FT_MEANDUMP_UNIT) THEN    ! Check for dump meaning
        D1_COPY_IN_MEMORY=.FALSE.
      ENDIF

      IF(D1_COPY_IN_MEMORY) THEN  ! Write to memory rather than disk
#if defined(ATMOS) && defined(OCEAN)
!L
!L     Copy D1 directly into memory for submodel
!L
      IF(sm_ident == atmos_sm) THEN
         DO I=1,LEN_DATA
          D1_A(I)=D1(I)
         ENDDO            ! I
      ELSEIF(sm_ident == ocean_sm) THEN
         DO I=1,LEN_DATA
          D1_O(I)=D1(I)
         ENDDO            ! I
      ELSE
         CMESSAGE='TRANSOUT: ERROR. Non-valid submodel identifier '
         write(6,*) CMESSAGE,sm_ident
         ICODE=1
         GO TO 999
      ENDIF
#endif
      write(6,*) 'TRANSOUT: Copied into memory LEN_DATA=',LEN_DATA,     &
     &                      'submodel=',sm_ident

      ELSE                       ! Write to disk rather than memory

      LEN_FILENAME=LEN(FILENAME)
      CALL FORT_GET_ENV(FT_ENVIRON(NFTOUT),LEN_FT_ENVIR(NFTOUT),        &
     &                  FILENAME,LEN_FILENAME,ICODE)

      IF(ICODE /= 0) THEN
         CMESSAGE='TRANSOUT: Environment variable not set '
         write(6,*) 'ERROR ',CMESSAGE,FT_ENVIRON(NFTOUT)
         GO TO 999
      ENDIF

!       Search for end of filename
      LL=0
      DO I=1,LEN_FILENAME
        IF(FILENAME(I:I) /= ' ') THEN
           LL=LL+1
        ENDIF
      ENDDO    ! I over characters

!       Construct filename with PE no. appended
      FILENAME(LL+1:LL+1)='.'
      WRITE(FILENAME(LL+2:LL+5),'(i4.4)') mype
!L
!L     Write out data
!L
        CALL OPEN_SINGLE(NFTOUT,FILENAME,LL+5,1,1,ICODE)
        CALL SETPOS_SINGLE(NFTOUT,0,ICODE)
#if defined(T3E)
!--compute the length of the first write
        disk_len_1=(len_data/um_sector_size)*um_sector_size
        call buffout_single(nftout, d1(1), disk_len_1, len_io, a)
        if(a /= -1. .or. len_io /= disk_len_1) then
          write(6,*) 'TRANSOUT: Error in data transfer to disk',        &
     &     '  A = ',a,'  LEN_IO = ',len_io,                             &
     &     '  Length Requested = ',disk_len_1
          icode=1
          cmessage='TRANSOUT: I/O write error'
          goto 999
        endif
!--now the remainder
        disk_len_2=len_data-disk_len_1
        if(disk_len_2 >  0) then
!--copy the rest of the data over
          do i=1, disk_len_2
            local_buffer(i)=d1(disk_len_1+i)
          end do
!--now output the remainder
          call buffout_single(nftout, local_buffer(1),                  &
     &     um_sector_size, len_io, a)
          if(a /= -1. .or. len_io /= um_sector_size) then
            write(6,*) 'TRANSOUT: Error in data transfer to disk',      &
     &       '  A = ',a,'  LEN_IO = ',len_io,                           &
     &       '  Length Requested = ',um_sector_size
            icode=1
            cmessage='TRANSOUT: I/O write error'
            goto 999
          endif
        endif
        call close_single(nftout, filename, ll+5, 1, 0, icode)
        write(6,*) 'TRANSOUT: Length transferred = ', len_data
#else
        CALL BUFFOUT_SINGLE(NFTOUT,D1(1),LEN_DATA,LEN_IO,A)
!L
!L----------------------------------------------------------------------
!L     Check for errors in data transfer to disk
!L----------------------------------------------------------------------
!L
        CALL CLOSE_SINGLE(NFTOUT,FILENAME,LL+5,1,0,ICODE)
          WRITE(6,*) 'TRANSOUT: Length transferred=',LEN_IO
          IF(A /= -1.0.OR.LEN_IO /= LEN_DATA)THEN
            WRITE(6,*) 'TRANSOUT: Error in data transfer to disk'
            ICODE=1
            CMESSAGE='TRANSOUT: I/O write error'
            GOTO 999
          ENDIF
#endif
      ENDIF                      ! End of disk/memory block

#else
!L
!L     Write out data (non-MPP)
!L
! DEPENDS ON: file_open
        CALL FILE_OPEN(NFTOUT,FT_ENVIRON(NFTOUT),                       &
     &            LEN_FT_ENVIR(NFTOUT),1,0,ICODE)
! DEPENDS ON: setpos
        CALL SETPOS(NFTOUT,0,ICODE)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,D1(1),LEN_DATA,LEN_IO,A)
!L
!L----------------------------------------------------------------------
!L     Check for errors in data transfer to disk
!L----------------------------------------------------------------------
!L
! DEPENDS ON: file_close
        CALL FILE_CLOSE(NFTOUT,FT_ENVIRON(NFTOUT),LEN_FT_ENVIR(NFTOUT)  &
     &                  ,0,0,ICODE)
          WRITE(6,*) 'TRANSOUT: Length transferred=',LEN_IO
          IF(A /= -1.0.OR.LEN_IO /= LEN_DATA)THEN
            WRITE(6,*) 'TRANSOUT: Error in data transfer to disk'
            ICODE=1
            CMESSAGE='TRANSOUT: I/O write error'
            GOTO 999
          ENDIF
#endif
!
 999  CONTINUE
      RETURN
      END SUBROUTINE TRANSOUT
#endif
