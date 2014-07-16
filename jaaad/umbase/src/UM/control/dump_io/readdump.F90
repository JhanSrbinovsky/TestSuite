#if defined(C80_1A) || defined(MAKEBC)
#if !defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE UM_READDUMP---------------------------------------
!LL
!LL  Purpose: Reads in model dump on unit NFTIN and checks model
!LL           and dump dimensions for consistency.
!LL
!LL Rewritten and tidied up for vn5.x
!LL
!LL Model   Date      Modification history from model version 5.0
!LL version
!LL 5.0     3/6/99    Rewritten UM_READDUMP for vn5.x       P.Burton
!LL 5.2     18/12/00  Add levels argument to read_multi  P.Burton
!LL 5.3     30/01/02  Remove MPP defs. Pass size of unpacked field
!LL                   to read_multi to match change in READFLDS.
!LL                   S.D.Mullerworth
!    5.5    02/08/00  Modification for parallelisation of WAM.
!                     Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
!    6.0    02/07/03  Remove MPP defs completely.   E.Leung
!    6.0    14/11/03  Fix to avoid problem with reading WF block
!                     into obs array - only select valid data,
!                     discarding packing area. M. Huddleston
!   6.1   22/10/04  Correct logic when printing file type. R. Hill
!    6.2    23/11/05  Removed all references to the wavemodel.
!                     T.Edwards
!    6.2    26/05/06  Correct External statements. P.Selwood
!    6.1    22/10/04  Correct logic when printing file type. R. Hill
!    6.2    06/12/05  Removed DIAG80 CPP define. T.Edwards
!LL
!LL
!LL

! Subroutine Interface

!LL  SUBROUTINE READDUMP---------------------------------------
!LL
!LL  Purpose: Reads in model obs dump on unit NFTIN and checks model
!LL           and dump dimensions for consistency.
!LL
!LL  Code mostly copied from original READDUMP
!LL
!LL  Model            Modification history from model version 4.3:
!LL version  Date
!LL   4.3  19/3/97   New deck introduced                    P.Burton
!     6.1  27/07/04  Correction to avoid problems when obs file only
!                    has a small number of obs. A. Hines.
!     6.2  01/06/06  Replace BUFFIN_ACOBS and BUFFIN_SHMEM with
!                    BUFFIN. P.Selwood
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL                        Version No 1 15/1/90
!LL
!LL  Logical component: R30
!LL
!LL  System task: F3
!LL
!LL  Documentation: Unified Model Documentation Paper No F3
!LL                 Version No 5 9/2/90
!LLEND---------------------------------------------------------
!
!*L Arguments:-------------------------------------------------
      SUBROUTINE READDUMP(NFTIN,FIXHD,LEN_FIXHD                         &
     & ,INTHD,LEN_INTHD                                                 &
     & ,REALHD,LEN_REALHD                                               &
     & ,LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC                               &
     & ,ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC                               &
     & ,COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC                               &
     & ,FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC                               &
     & ,EXTCNST,LEN_EXTCNST                                             &
     & ,DUMPHIST,LEN_DUMPHIST                                           &
     & ,CFI1,LEN_CFI1                                                   &
     & ,CFI2,LEN_CFI2                                                   &
     & ,CFI3,LEN_CFI3                                                   &
     & ,LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP                                  &
     &      ,LEN_DATA,D1,                                               &
#include "argppx.h"
     &  ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER                                                           &
     & NFTIN                                                            &
                     !IN Unit no of dump
     &,LEN_FIXHD                                                        &
                     !IN Length of fixed length header
     &,LEN_INTHD                                                        &
                     !IN Length of integer header
     &,LEN_REALHD                                                       &
                     !IN Length of real header
     &,LEN1_LEVDEPC                                                     &
                     !IN 1st dim of level dep consts
     &,LEN2_LEVDEPC                                                     &
                     !IN 2nd dim of level dep consts
     &,LEN1_ROWDEPC                                                     &
                     !IN 1st dim of row dep consts
     &,LEN2_ROWDEPC                                                     &
                     !IN 2nd dim of row dep consts
     &,LEN1_COLDEPC                                                     &
                     !IN 1st dim of column dep consts
     &,LEN2_COLDEPC                                                     &
                     !IN 2nd dim of column dep consts
     &,LEN1_FLDDEPC                                                     &
                     !IN 1st dim of field dep consts
     &,LEN2_FLDDEPC                                                     &
                     !IN 2nd dim of field dep consts
     &,LEN_EXTCNST                                                      &
                     !IN Length of extra constants
     &,LEN_DUMPHIST                                                     &
                     !IN Length of history block
     &,LEN_CFI1                                                         &
                     !IN Length of comp field index 1
     &,LEN_CFI2                                                         &
                     !IN Length of comp field index 2
     &,LEN_CFI3                                                         &
                     !IN Length of comp field index 3
     &,LEN1_LOOKUP                                                      &
                     !IN 1st dim of lookup
     &,LEN2_LOOKUP   !IN 2nd dim of lookup

      INTEGER                                                           &
     & LEN_DATA                                                         &
                      !IN Length of model data
     &,ICODE          !OUT Return code; successful=0
                      !                 error > 0

      CHARACTER*(80)                                                    &
     & CMESSAGE       !OUT Error message if ICODE > 0

      INTEGER                                                           &
     & FIXHD(LEN_FIXHD)                                                 &
                        !IN Fixed length header
     &,INTHD(LEN_INTHD)                                                 &
                        !IN Integer header
     &,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                                  &
                                       !IN PP lookup tables

     &,CFI1(LEN_CFI1+1)                                                 &
                        !IN Compressed field index no 1
     &,CFI2(LEN_CFI2+1)                                                 &
                        !IN Compressed field index no 2
     &,CFI3(LEN_CFI3+1) !IN Compressed field index no 3

      REAL                                                              &
     & REALHD(LEN_REALHD)                                               &
                          !IN Real header
     &,LEVDEPC(1+LEN1_LEVDEPC*LEN2_LEVDEPC)                             &
                                            !IN Lev dep consts
     &,ROWDEPC(1+LEN1_ROWDEPC*LEN2_ROWDEPC)                             &
                                            !IN Row dep consts
     &,COLDEPC(1+LEN1_COLDEPC*LEN2_COLDEPC)                             &
                                            !IN Col dep consts
     &,FLDDEPC(1+LEN1_FLDDEPC*LEN2_FLDDEPC)                             &
                                            !IN Field dep consts
     &,EXTCNST(LEN_EXTCNST+1)                                           &
                                !IN Extra constants
     &,DUMPHIST(LEN_DUMPHIST+1)                                         &
                                !IN History block
     &,D1(*)

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "clookadd.h"
#include "parvars.h"
#include "c_mdi.h"
#include "cntl_io.h"

! -------------------------------------------------------------
! Local arrays:------------------------------------------------

! -------------------------------------------------------------
!*L External subroutines called:-------------------------------
      EXTERNAL IOERROR,POSERROR,READHEAD,PR_LOOK                        &
     &,BUFFIN,EXPAND32B
! Cray specific functions  UNIT,LENGTH
!*-------------------------------------------------------------
! Local variables:---------------------------------------------
      INTEGER START_BLOCK                                               &
                           ! Pointer to current position in file
     &,LEN_IO                                                           &
                           ! No of 64-bit words buffered in
     &,K,I                                                              &
                           ! Loop counts
     &,IPTS                ! No of 64-bit words requested to be
                           ! buffered in
      REAL A               ! Error code returned by UNIT
!
      integer real_start_block                                          &
                               ! Real disk address
     & , l                                                              &
                               ! loop counter
     & , word_address                                                   &
                               ! word address on disk of the record
     & , um_sector_ipts                                                 &
                               ! number fo words to read, rounded up
                               ! to a sector size
     & , l_ipts                ! local value of ipts for address calc.
!--------------------------------------------------------------

      IF (mype  ==  0) THEN
      WRITE(6,'(/,'' READING UNIFIED MODEL DUMP ON UNIT'',I3)')NFTIN
      WRITE(6,'('' #####################################'',/)')
      ENDIF
      ICODE=0
      CMESSAGE=' '

!L 1. Read in all header records and check for consistency.
!     START_BLOCK points to position of model data block
!     on return

! DEPENDS ON: readhead
      CALL READHEAD(NFTIN,FIXHD,LEN_FIXHD,                              &
     &              INTHD,LEN_INTHD,                                    &
     &              REALHD,LEN_REALHD,                                  &
     &              LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                  &
     &              ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,                  &
     &              COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                  &
     &              FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,                  &
     &              EXTCNST,LEN_EXTCNST,                                &
     &              DUMPHIST,LEN_DUMPHIST,                              &
     &              CFI1,LEN_CFI1,                                      &
     &              CFI2,LEN_CFI2,                                      &
     &              CFI3,LEN_CFI3,                                      &
     &              LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,                     &
     &              LEN_DATA,                                           &
#include "argppx.h"
     &              START_BLOCK,ICODE,CMESSAGE)

      IF(ICODE >  0)RETURN


!L 2. Buffer in model data one field at a time for
!L    conversion from 32-bit to 64-bit numbers

      IF(FIXHD(160) >  0)THEN

! Check for error in file pointers
       real_start_block=start_block
       if(start_block /= fixhd(160)) then
! If new format Dumpfile, we must reset the start address
         if((lookup(lbnrec,1) == 0.and.lookup(lblrec,1) >  0) .or.      &
! Ocean ACOBS Files (?)
     &     ((lookup(lbnrec,1) == imdi) .or. (lookup(lbegin,1) == imdi)) &
     &     .or.                                                         &
! Prog lookups in dump before vn3.2:
     &     ((lookup(lbnrec,1) == imdi) .and. (fixhd(12) <= 301))) then
        CMESSAGE='READDUMP: Addressing conflict'
        ICODE=1
! DEPENDS ON: poserror
        CALL POSERROR('model data',                                     &
     &  START_BLOCK,160,FIXHD(160))
        RETURN
         else
           real_start_block=fixhd(160)
         endif
       ENDIF

!      Move to start of data.
! DEPENDS ON: setpos
       CALL SETPOS (NFTIN,FIXHD(160)-1,ICODE)

! Loop over number of fields in data block
       DO 200 K=1,FIXHD(152)

       IF (LOOKUP(LBLREC,K) >  0) THEN   !  Any data for this field ?

! Test whether data stored as 32-bit on disk
        IF (MOD((LOOKUP(LBPACK,K)),10) == 2) THEN
         IPTS=(LOOKUP(LBLREC,K)+1)/2
        ELSE
         IPTS=LOOKUP(LBLREC,K)
        ENDIF

! Compute word address in file from which to begin I/O

! Old Format dumpfiles
        if((lookup(lbnrec,k) == 0) .or.                                 &
! Ocean ACOBS Files (?)
     &    ((lookup(lbnrec,k) == imdi) .or. (lookup(lbegin,k) == imdi))  &
     &    .or.                                                          &
! Prog lookups in dump before vn3.2:
     &    ((lookup(lbnrec,k) == imdi) .and. (fixhd(12) <= 301))) then
! Dump and ancillary files
          word_address=1
          if(k >  1)then
            do l=2,k
              if(mod(lookup(lbpack,l-1),10) == 2) then
                l_ipts=(lookup(lblrec,l-1)+1)/2
              else
                l_ipts=(lookup(lblrec,l-1))
              endif
              word_address=word_address+l_ipts
            end do
          endif
          word_address=fixhd(160)+word_address-2
          um_sector_ipts=ipts

        else

! PP type files and new format Dumpfiles (vn4.4 onwards)
          word_address=lookup(lbegin,k)
! Use the stored round-up value
          um_sector_ipts=lookup(lbnrec,k)
       LEN_IO=lookup(lblrec,k)
        endif

! To avoid problems when reading acobs files, set the size of data
! to be read in to be the real size of the data, and not the
! size including the padding.

      UM_SECTOR_IPTS=IPTS

! Position file pointer
! DEPENDS ON: setpos
        call setpos(nftin,word_address,icode)


! Read data into final position
! Check that data_type is valid no: 1 to 3 or -1 to -3
        IF((LOOKUP(DATA_TYPE,K) >= 1.AND.LOOKUP(DATA_TYPE,K) <= 3) .OR. &
     &     (LOOKUP(DATA_TYPE,K) <= -1.AND.LOOKUP(DATA_TYPE,K) >= -3))   &
     &     THEN
        ipts=um_sector_ipts
! DEPENDS ON: buffin
       CALL BUFFIN(NFTIN,D1(LOOKUP(NADDR,K)),IPTS,LEN_IO,A)
       IF (MYPE == 0) WRITE(0,*) K,lookup(lblrec,k),                    &
     &                       'READ IN TO',LOOKUP(NADDR,K)
        IF ((A /= -1.0).OR.(LEN_IO /= IPTS)) THEN
          WRITE(6,*)'ERROR READING DUMP ON UNIT ',NFTIN
          ICODE=2
          CMESSAGE='READDUMP: BAD BUFFIN OF DATA'
! DEPENDS ON: ioerror
          CALL IOERROR('BUFFER IN FROM READDUMP',A,LEN_IO,IPTS)
          RETURN
        END IF
! Error in lookup(data_type,k)
      ELSE
        ICODE=3
        CMESSAGE='READDUMP:  Invalid code in LOOKUP(DATA_TYPE,K)'
      END IF

      ENDIF  !  Skip to here if no data for this field

       START_BLOCK=START_BLOCK+LOOKUP(LBLREC,K)
       real_start_block=real_start_block+um_sector_ipts

200   CONTINUE

      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       IF (FIXHD(5) >= 6 .AND. FIXHD(5) <= 8) THEN ! AC/Var Obs/ Cx file
         WRITE(6,'('' OBSERVATION DATA'')')
       ELSE
         WRITE(6,'('' MODEL DATA'')')
       ENDIF
       WRITE(6,'('' '',I8,'' words long'')')FIXHD(161)
      ENDIF ! mype  ==  0

      ENDIF

      IF (mype  ==  0) THEN
      WRITE(6,'('' '')')
      WRITE(6,'('' INITIAL DATA SUCCESSFULLY READ -'',I9,               &
     &'' WORDS FROM UNIT'',I3)')START_BLOCK,NFTIN
       if(real_start_block /= start_block) then
         write(6,'(/'' Number of Words Read from Disk was '',i9)')      &
     &    real_start_block
       endif
      ENDIF ! mype  ==  0

 9999 CONTINUE
      RETURN
      END SUBROUTINE READDUMP
#endif
#endif
