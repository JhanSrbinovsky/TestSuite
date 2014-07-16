#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL
!LL    Subroutine: TRANSIN  ---------------------------------------
!LL
!LL    Purpose:
!LL    To transfer dump data from disk to memory en masse
!LL
!LL    Tested under compiler: cft77
!LL
!LL    Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL 3.4  16/6/94 : Change CHARACTER*(*) to CHARACTER*(80) N.Farnon
!LL   3.3    07/10/93 Corrected order of *CALLs to comdecks TYPSIZE
!LL                   and TYPD1.  Tracey Smith
!LL   3.5  24/03/95    Changed OPEN to FILE_OPEN and
!LL                    CLOSE to FILE_CLOSE    P.Burton
!LL   4.2  11/10/96    Enable atmos-ocean coupling for MPP.
!LL                    (2): Swap D1 memory.
!LL                    Image of D1 either copied directly from local
!LL                    memory or I/O from file (local to processor)
!LL                    under MPP. I/O kept for mean dumps. R.Rawlins
!LL   4.3  30/01/97    Ensure that domain decomposition is consistent
!LL                    with submodel. R.Rawlins
!LL   4.4  01/07/97    Make transfers to the input file
!LL                    well-formed.
!LL                      Author: Bob Carruthers, Cray Research.
!LL   4.4  11/10/97    Call CNTLALL for L_AO_D1_MEMORY. D. Robinson.
!LL   4.4  28/08/97    Minor tidy: replace SETPOS by SETPOS_SINGLE for
!LL                    MPP case. R.Rawlins
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
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
!LL
!LLEND --------------------------------------------------------------
!*L    Interface and arguments:
      SUBROUTINE TRANSIN(                                               &
#include "argd1.h"
     &  LEN_DATA,NFTIN,sm_ident                                         &
     & ,ICODE,CMESSAGE)
!
      IMPLICIT NONE
!
      INTEGER                                                           &
     &       LEN_DATA,                                                  &
                                  ! IN Length of model data
     &       NFTIN,                                                     &
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
      EXTERNAL SETPOS
#if defined(MPP)
      EXTERNAL FORT_GET_ENV,OPEN_SINGLE,CLOSE_SINGLE
      EXTERNAL BUFFIN_SINGLE,SETPOS_SINGLE
      EXTERNAL CHANGE_DECOMPOSITION
#else
      EXTERNAL FILE_OPEN,FILE_CLOSE
      EXTERNAL BUFFIN
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
     &      ,decomp_standard                                            &
                                  ! MPP domain decomposition ident
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

#include "decomptp.h"
#include "mpptrans.h"
#include "cenvir.h"
#include "csmid.h"
#include "chsunits.h"
#include "cntlall.h"

#if defined(MPP)
      D1_COPY_IN_MEMORY=L_AO_D1_MEMORY       ! from COMDECK CNTLALL
      IF(NFTIN == FT_MEANDUMP_UNIT) THEN     ! Check for dump meaning
         D1_COPY_IN_MEMORY=.FALSE.
      ENDIF

      IF(D1_COPY_IN_MEMORY) THEN  ! Read from memory rather than disk
#if defined(ATMOS) && defined(OCEAN)
!L
!L     Copy D1 directly from memory for submodel
!L
      IF(sm_ident == atmos_sm) THEN
         DO I=1,LEN_DATA
          D1(I)=D1_A(I)
         ENDDO            ! I
      ELSEIF(sm_ident == ocean_sm) THEN
         DO I=1,LEN_DATA
          D1(I)=D1_O(I)
         ENDDO            ! I
      ELSE
         CMESSAGE='TRANSIN: ERROR. Non-valid submodel identifier '
         write(6,*) CMESSAGE,sm_ident
         ICODE=1
         GO TO 999
      ENDIF
#endif
      write(6,*) 'TRANSIN : Copied from memory LEN_DATA=',LEN_DATA,     &
     &                      'submodel=',sm_ident

      ELSE                       ! Read from disk rather than memory

!L
!L        Read from disk file rather than memory
!L
      LEN_FILENAME=LEN(FILENAME)
      CALL FORT_GET_ENV(FT_ENVIRON(NFTIN),LEN_FT_ENVIR(NFTIN),          &
     &                  FILENAME,LEN_FILENAME,ICODE)

      IF(ICODE /= 0) THEN
         CMESSAGE='TRANSIN : Environment variable not set '
         write(6,*) 'ERROR ',CMESSAGE,FT_ENVIRON(NFTIN)
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
!L     Read in data
!L
        CALL OPEN_SINGLE(NFTIN,FILENAME,LL+5,0,1,ICODE)
        CALL SETPOS_SINGLE(NFTIN,0,ICODE)
#if defined(T3E)
!--compute the length of the first read
        disk_len_1=(len_data/um_sector_size)*um_sector_size
        call buffin_single(nftin, d1(1), disk_len_1, len_io, a)
        if(a /= -1. .or. len_io /= disk_len_1) then
          write(6,*) 'TRANSIN: Error in data transfer from disk',       &
     &     '  A = ',a,'  LEN_IO = ',len_io,                             &
     &     '  Length Requested = ',disk_len_1
          icode=1
          cmessage='TRANSIN: I/O read error'
          goto 999
        endif
!--now the remainder
        disk_len_2=len_data-disk_len_1
        if(disk_len_2 >  0) then
          call buffin_single(nftin, local_buffer(1),                    &
     &     um_sector_size, len_io, a)
          if(a /= -1. .or. len_io /= um_sector_size) then
            write(6,*) 'TRANSIN: Error in data transfer from disk',     &
     &       '  A = ',a,'  LEN_IO = ',len_io,                           &
     &       '  Length Requested = ',um_sector_size
            icode=1
            cmessage='TRANSIN: I/O read error'
            goto 999
          endif
!--copy the rest of the data over
          do i=1, disk_len_2
            d1(disk_len_1+i)=local_buffer(i)
          end do
        endif
        call close_single(nftin, filename, ll+5, 1, 0, icode)
        write(6,*) 'TRANSIN: Length transferred = ', len_data
#else
        CALL BUFFIN_SINGLE(NFTIN,D1(1),LEN_DATA,LEN_IO,A)
!L
!L----------------------------------------------------------------------
!L     Check for errors in data transfer from disk
!L----------------------------------------------------------------------
!L
        CALL CLOSE_SINGLE(NFTIN,FILENAME,LL+5,1,0,ICODE)
          WRITE(6,*) 'TRANSIN: Length transferred=',LEN_IO
          IF(A /= -1.0.OR.LEN_IO /= LEN_DATA)THEN
            WRITE(6,*) 'TRANSIN: Error in data transfer from disk'
            ICODE=1
            CMESSAGE='TRANSIN: I/O read error'
            GOTO 999
          ENDIF
#endif
      ENDIF                       ! End of disk/memory block

#if defined(ATMOS) && defined(OCEAN)
!L
!L    Ensure that domain decomposition is consistent with submodel
!L
      IF(sm_ident == atmos_sm) THEN
         decomp_standard = decomp_standard_atmos
      ELSEIF(sm_ident == ocean_sm) THEN
         decomp_standard = decomp_standard_ocean
      ENDIF

! DEPENDS ON: change_decomposition
      CALL CHANGE_DECOMPOSITION(decomp_standard,ICODE)
#endif

#else

!L
!L     Read in data (non-MPP)
!L
! DEPENDS ON: file_open
        CALL FILE_OPEN(NFTIN,FT_ENVIRON(NFTIN),                         &
     &            LEN_FT_ENVIR(NFTIN),0,0,ICODE)
! DEPENDS ON: setpos
        CALL SETPOS(NFTIN,0,ICODE)
! DEPENDS ON: buffin
        CALL BUFFIN(NFTIN,D1(1),LEN_DATA,LEN_IO,A)
!L
!L----------------------------------------------------------------------
!L     Check for errors in data transfer from disk
!L----------------------------------------------------------------------
!L
! DEPENDS ON: file_close
        CALL FILE_CLOSE(NFTIN,FT_ENVIRON(NFTIN),LEN_FT_ENVIR(NFTIN),    &
     &                  0,0,ICODE)
          WRITE(6,*) 'TRANSIN: Length transferred=',LEN_IO
          IF(A /= -1.0.OR.LEN_IO /= LEN_DATA)THEN
            WRITE(6,*) 'TRANSIN: Error in data transfer from disk'
            ICODE=1
            CMESSAGE='TRANSIN: I/O read error'
            GOTO 999
          ENDIF
#endif
!
 999  CONTINUE
      RETURN
      END SUBROUTINE TRANSIN
#endif
