#if defined(C80_1A) || defined(UTILIO) || defined(RECON)               \
 || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE READHEAD---------------------------------------
!LL
!LL AD, DR      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL   3.2  06/05/93    Skip call to CHKLOOK if PP type file
!LL                    Author: A. Dickinson    Reviewer: D. Richardson
!LL
!LL  3.1   22/12/92     Allow use by ancillary field headers
!LL                     Author A. Dickinson    Reviewer C. Wilson
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL  3.2   12/05/93     Adapt to read prognostic fields only.
!LL                     Author D. Robinson     Reviewer A. Dickinson
!LL   3.5    28/03/95 MPP code: New code for parallel I/O
!LL                                              P.Burton
!     3.5  21/06/95  Set lookup(45) if initial dump pre version 3.5
!                    Author D.M.Goddard    Reviewer S Swarbrick
!     4.0  06/10/95  Set variable MODEL for all diagnostics in dump
!                    Author D.M. Goddard
!     4.1  23/05/96  Removed resetting of FIXHD(161) for MPP code
!                    P.Burton
!     4.1  18/06/96  Changes to cope with changes in STASH addressing
!                    Author D.M. Goddard.
!     5.1  10/04/00 Control stdout with PrintStatus and new
!                   reconfiguration support. P.Selwood.
!     5.3  22/11/01 Enable MPP as the only option for
!                   small executables         E.Leung
!     6.2  15/08/05 Free format fixes. P.Selwood

!LL  Programming standard: Unified Model Documentation Paper No 3
!LL                        Version No 1 15/1/90
!LL
!LL  Logical component: R30
!LL
!LL  System task: F3
!LL
!LL  Purpose: Reads in model dump header records on unit NFTIN and
!LL           checks model and dump dimensions for consistency.
!LL
!LL  Documentation: Unified Model Documentation Paper No F3
!LL                 Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE READHEAD(NFTIN,FIXHD,LEN_FIXHD,                        &
                                                       ! Intent (In)
     &                    INTHD,LEN_INTHD,                              &
     &                    REALHD,LEN_REALHD,                            &
     &                    LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,            &
     &                    ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,            &
     &                    COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,            &
     &                    FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,            &
     &                    EXTCNST,LEN_EXTCNST,                          &
     &                    DUMPHIST,LEN_DUMPHIST,                        &
     &                    CFI1,LEN_CFI1,                                &
     &                    CFI2,LEN_CFI2,                                &
     &                    CFI3,LEN_CFI3,                                &
     &                    LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,      &
#include "argppx.h"
     &                    START_BLOCK,ICODE,CMESSAGE)   ! Intent (Out)
#if defined(RECON) || defined(VAROPSVER)
      Use Rcf_Parvars_Mod, Only :                                       &
     &    mype

      Use Rcf_PrintStatus_Mod, Only :                                   &
     &    PrintStatus,                                                  &
     &    PrStatus_Oper,                                                &
     &    PrStatus_Normal
#endif

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
                     !IN 2ndt dim of level dep consts
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
     &,START_BLOCK                                                      &
                      !OUT Pointer to position of each block.
                      !Should point to start of model data block on exit
     &,ICODE          !OUT Return code; successful=0
                      !                 error > 0

#if defined(RECON) || defined(VAROPSVER)
! dummy declarations for reconfiguration
      Integer ppxi
      Integer ppxrecs
      Character ppxc
#endif
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
     &,DUMPHIST(LEN_DUMPHIST+1) !IN History block

! Local arrays:------------------------------------------------
! None
! -------------------------------------------------------------
! External subroutines called:---------------------------------
      EXTERNAL IOERROR,POSERROR,PR_FIXHD,CHK_LOOK,BUFFIN
!*-------------------------------------------------------------
! Comdecks:----------------------------------------------------------
#if !defined(RECON) && !defined(VAROPSVER)
#include "cprintst.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#endif
#include "c_mdi.h"
#include "clookadd.h"
#if !defined(RECON)
#include "parvars.h"
#endif
! Local variables:---------------------------------------------
      INTEGER K
      INTEGER LEN_IO
      INTEGER FIXHD_152    !  Original value of FIXHD(152)
      LOGICAL L_A_DUMP
      LOGICAL L_O_DUMP
      LOGICAL L_FF    ! FieldsFile Logical
      REAL A
! -------------------------------------------------------------

      ICODE=0
      CMESSAGE=' '

!L 1. Buffer in fixed length header record

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,FIXHD(1),LEN_FIXHD,LEN_IO,A)


! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= LEN_FIXHD)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header',A,LEN_IO        &
     &               ,LEN_FIXHD)
        CMESSAGE='READHEAD: I/O error'
        ICODE=1
        RETURN
      ENDIF

      START_BLOCK=LEN_FIXHD+1

      FIXHD_152 = FIXHD(152)    !  Store original value

!     Test if atmos dump read in
      L_A_DUMP = FIXHD(5) == 1 .AND. FIXHD(2) == 1                      &
     &           .AND. LEN_DATA /= IMDI

!     Test if ocean dump read in
      L_O_DUMP = FIXHD(5) == 1 .AND. FIXHD(2) == 2                      &
     &           .AND. LEN_DATA /= IMDI

!      Test if FieldsFile read in
      L_FF = FIXHD(5) == 3

#if defined(RECON) || defined(VAROPSVER)
      IF (L_A_DUMP .OR. L_O_DUMP .OR. L_FF) THEN
#else
      IF (L_A_DUMP .OR. L_O_DUMP) THEN
#endif
        IF (FIXHD(152) /= LEN2_LOOKUP) THEN
!XX       WRITE (6,*) 'FIXHD(152) being reset from ',FIXHD(152),' to ',
!XX  *    LEN2_LOOKUP
          FIXHD(152) = LEN2_LOOKUP
        ENDIF
#if !defined(MPP)
        IF (FIXHD(161) /= LEN_DATA) THEN
!XX       WRITE (6,*) 'FIXHD(161) being reset from ',FIXHD(161),' to ',
!XX  *    LEN_DATA
          FIXHD(161) = LEN_DATA
        ENDIF
#endif
      ENDIF

! Check validity of data and print out fixed header information

      IF (mype  ==  0) THEN
! DEPENDS ON: pr_fixhd
      CALL PR_FIXHD(FIXHD,LEN_FIXHD,LEN_INTHD,LEN_REALHD,LEN1_LEVDEPC   &
     &,LEN2_LEVDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC &
     &,LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,LEN_DUMPHIST,LEN_CFI1      &
     &,LEN_CFI2,LEN_CFI3,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA               &
     &,ICODE,CMESSAGE)

      IF(ICODE >  0)RETURN

      ENDIF
!L 2. Buffer in integer constants

      IF(FIXHD(100) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(100) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('integer constants',START_BLOCK,100,FIXHD(100))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=2
        RETURN
       ENDIF

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,INTHD(1),FIXHD(101),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(101))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of integer constants',A,LEN_IO,         &
     &               FIXHD(101))
        CMESSAGE='READHEAD: I/O error'
        ICODE=3
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(101)

      ENDIF

!L 3. Buffer in real constants

      IF(FIXHD(105) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(105) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('real constants',START_BLOCK,105,FIXHD(105))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=4
        RETURN
       ENDIF

! Check for I/O errors
! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,REALHD(1),FIXHD(106),LEN_IO,A)

       IF(A /= -1.0.OR.LEN_IO /= FIXHD(106))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of real constants',A,LEN_IO,            &
     &                FIXHD(106))
        CMESSAGE='READHEAD: I/O error'
        ICODE=5
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(106)


      ENDIF

!L 4. Buffer in level dependent constants

      IF(FIXHD(110) >  0.AND.LEN1_LEVDEPC /= 0)THEN

! Check for error in file pointers
       IF(FIXHD(110) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('level dependent constants',                      &
     &  START_BLOCK,110,FIXHD(110))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=6
        RETURN
       ENDIF

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,LEVDEPC(1),FIXHD(111)*FIXHD(112),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(111)*FIXHD(112))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of level dependent constants',A,LEN_IO, &
     &               FIXHD(111)*FIXHD(112))
        CMESSAGE='READHEAD: I/O error'
        ICODE=7
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(111)*FIXHD(112)

      If (PrintStatus >= PrStatus_Oper) Then
       IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' LEVEL DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(111)*FIXHD(112)

      ENDIF ! mype  ==  0
      End If
      ENDIF

!L 5. Buffer in row dependent constants

      IF(FIXHD(115) >  0.AND.LEN1_ROWDEPC /= 0)THEN

! Check for error in file pointers
       IF(FIXHD(115) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('row dependent constants',                        &
     &  START_BLOCK,115,FIXHD(115))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=8
        RETURN
       ENDIF

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,ROWDEPC(1),FIXHD(116)*FIXHD(117),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(116)*FIXHD(117))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of row dependent constants',A,LEN_IO,   &
     &                FIXHD(116)*FIXHD(117))
        CMESSAGE='READHEAD: I/O error'
        ICODE=9
        RETURN
      ENDIF


       START_BLOCK=START_BLOCK+FIXHD(116)*FIXHD(117)

      If (PrintStatus >= PrStatus_Oper) Then
       IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' ROW DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(116)*FIXHD(117)

      ENDIF ! mype  ==  0
      End If
      ENDIF

!L 6. Buffer in column dependent constants

      IF(FIXHD(120) >  0.AND.LEN1_COLDEPC /= 0)THEN

! Check for error in file pointers
       IF(FIXHD(120) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('column dependent constants',                     &
     &  START_BLOCK,120,FIXHD(120))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=10
        RETURN
       ENDIF

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,COLDEPC(1),FIXHD(121)*FIXHD(122),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(121)*FIXHD(122))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of column dependent constants',A,LEN_IO,&
     &               FIXHD(121)*FIXHD(122))
        CMESSAGE='READHEAD: I/O error'
        ICODE=11
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(121)*FIXHD(122)

      If (PrintStatus >= PrStatus_Oper) Then
       IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COLUMN DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(121)*FIXHD(122)

      ENDIF ! mype  ==  0
      End If
      ENDIF

!L 7. Buffer in constants stored as fields

      IF(FIXHD(125) >  0.AND.LEN1_FLDDEPC /= 0)THEN

! Check for error in file pointers
       IF(FIXHD(125) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('fields of constants',                            &
     &  START_BLOCK,125,FIXHD(125))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=12
        RETURN
       ENDIF

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,FLDDEPC(1),FIXHD(126)*FIXHD(127),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(126)*FIXHD(127))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of field dependent constants',A,LEN_IO, &
     &               FIXHD(126)*FIXHD(127))
        CMESSAGE='READHEAD: I/O error'
        ICODE=13
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(126)*FIXHD(127)

      If (PrintStatus >= PrStatus_Oper) Then
       IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' FIELD DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(126)*FIXHD(127)

      ENDIF ! mype  ==  0
      End If
      ENDIF

!L 8. Buffer in extra constants

      IF(FIXHD(130) >  0.AND.LEN_EXTCNST /= 0)THEN

! Check for error in file pointers
       IF(FIXHD(130) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('extra constants',                                &
     &  START_BLOCK,130,FIXHD(130))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=14
        RETURN
       ENDIF

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,EXTCNST(1),FIXHD(131),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(131))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in extra constants',A,LEN_IO,              &
     &               FIXHD(131))
        CMESSAGE='READHEAD: I/O error'
        ICODE=15
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(131)

      If (PrintStatus >= PrStatus_Oper) Then
       IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' EXTRA CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(131)

      ENDIF ! mype  ==  0
      End If
      ENDIF

!L 9. Buffer in temporary history block

      IF(FIXHD(135) >  0.AND.LEN_DUMPHIST /= 0)THEN

! Check for error in file pointers
       IF(FIXHD(135) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('history',                                        &
     &  START_BLOCK,136,FIXHD(136))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=16
        RETURN
       ENDIF

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,DUMPHIST(1),FIXHD(136),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(136))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of history file',A,LEN_IO,              &
     &               FIXHD(136))
        CMESSAGE='READHEAD: I/O error'
        ICODE=17
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(136)

      If (PrintStatus >= PrStatus_Oper) Then
       IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' TEMPORARY HISTORY BLOCK'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(136)

      ENDIF ! mype  ==  0
      End If
      ENDIF

!L 10. Buffer in compressed field index1

      IF(FIXHD(140) >  0.AND.LEN_CFI2 /= 0)THEN

! Check for error in file pointers

       IF(FIXHD(140) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('compressed field index1',                        &
     &  START_BLOCK,140,FIXHD(140))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=18
        RETURN
       ENDIF

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,CFI1(1),FIXHD(141),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(141))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of compressed index1',A,LEN_IO,         &
     &               FIXHD(141))
        CMESSAGE='READHEAD: I/O error'
        ICODE=19
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(141)

      If (PrintStatus >= PrStatus_Oper) Then
       IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 1'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(141)

      ENDIF ! mype  ==  0
      End If
      ENDIF

!L 11. Buffer in compressed field index2

      IF(FIXHD(142) >  0.AND.LEN_CFI2 /= 0)THEN

! Check for error in file pointers
       IF(FIXHD(142) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('compressed field index2',                        &
     &  START_BLOCK,142,FIXHD(142))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=20
        RETURN
       ENDIF

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,CFI2(1),FIXHD(143),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(143))THEN
! DEPENDS ON: ioerror
       CALL IOERROR('buffer in of compressed index2',A,LEN_IO,          &
     &               FIXHD(143))
        CMESSAGE='READHEAD: I/O error'
        ICODE=21
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(143)

      If (PrintStatus >= PrStatus_Oper) Then
       IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 2'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(143)

      ENDIF ! mype  ==  0
      End If
      ENDIF

!L 12. Buffer in compressed field index3

      IF(FIXHD(144) >  0.AND.LEN_CFI3 /= 0)THEN

! Check for error in file pointers
       IF(FIXHD(144) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('compressed field index3',                        &
     &  START_BLOCK,144,FIXHD(144))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=22
        RETURN
       ENDIF

! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,CFI3(1),FIXHD(145),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(145))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of compressed index3',A,LEN_IO,         &
     &               FIXHD(145))
        CMESSAGE='READHEAD: I/O error'
        ICODE=23
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(145)

      If (PrintStatus >= PrStatus_Oper) Then
       IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 3'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(145)

      ENDIF ! mype  ==  0
      End If
      ENDIF

!L 13. Buffer in lookup table

      IF(FIXHD(150) >  0)THEN

! Supress checking if not full dump
      IF(LEN_DUMPHIST /= 0)THEN
! Check for error in file pointers
       IF(FIXHD(150) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('lookup table',                                   &
     &  START_BLOCK,150,FIXHD(150))
        CMESSAGE='READHEAD: Addressing conflict'
        ICODE=24
        RETURN
       ENDIF
      ENDIF

! Move to start of Look Up Table
! DEPENDS ON: setpos
      CALL SETPOS(NFTIN,FIXHD(150)-1,ICODE)

! Read in fields from LOOKUP table
! DEPENDS ON: buffin
      CALL BUFFIN(NFTIN,LOOKUP(1,1),FIXHD(151)*FIXHD(152),LEN_IO,A)

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(151)*FIXHD(152))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of lookup table',A,LEN_IO,              &
     &               FIXHD(151)*FIXHD(152))
        CMESSAGE='READHEAD: I/O error'
        ICODE=25
        RETURN
       ENDIF

! Point to start of data section ( Use original FIXHD(152) )
       START_BLOCK=START_BLOCK+FIXHD(151)*FIXHD_152

      If (PrintStatus >= PrStatus_Oper) Then
       IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' LOOKUP TABLE'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(151)*FIXHD(152)

       IF (FIXHD(152) <  FIXHD_152) THEN
         WRITE(6,'('' '')')
         WRITE(6,'('' '',I6,'' Entries in Look Up Table.'')') FIXHD_152
         WRITE(6,'('' '',I6,'' Entries read in.'')') FIXHD(152)
       ENDIF

      ENDIF ! mype  ==  0
      End If
!---------------------------------------------------------------
! Reset LOOKUP(45) if not set
!---------------------------------------------------------------

        DO K=1,LEN2_LOOKUP
          IF(LOOKUP(45,K) == 0.OR.LOOKUP(45,K) == IMDI)THEN

!Section 0: Prognostic fields.
          IF(LOOKUP(42,K) <= 100.OR.                                    &
     &      (LOOKUP(42,K) >= 200.AND.LOOKUP(42,K) <= 205))THEN
            LOOKUP(45,K)=1

          ELSE IF((LOOKUP(42,K) >  100.AND.LOOKUP(42,K) <= 176).OR.     &
     &            (LOOKUP(42,K) >= 180.AND.LOOKUP(42,K) <  200))THEN
            LOOKUP(45,K)=2

          ELSE IF((LOOKUP(42,K) >= 177.AND.LOOKUP(42,K) <= 179).OR.     &
     &            (LOOKUP(42,K) >= 210.AND.LOOKUP(42,K) <= 212))THEN
            LOOKUP(45,K)=3

! Sections 1 - 99: Diagnostic fields
          ELSE IF(LOOKUP(42,K) >= 1000.AND.LOOKUP(42,K) <= 29999)THEN
            IF((LOOKUP(42,K) >= 21177.AND.LOOKUP(42,K) <= 21179).OR.    &
     &         (LOOKUP(42,K) >= 21225.AND.LOOKUP(42,K) <= 21227).OR.    &
     &         (LOOKUP(42,K) >= 22177.AND.LOOKUP(42,K) <= 22179).OR.    &
     &         (LOOKUP(42,K) >= 22225.AND.LOOKUP(42,K) <= 22227).OR.    &
     &         (LOOKUP(42,K) >= 23177.AND.LOOKUP(42,K) <= 23179).OR.    &
     &         (LOOKUP(42,K) >= 23225.AND.LOOKUP(42,K) <= 23227).OR.    &
     &         (LOOKUP(42,K) >= 24177.AND.LOOKUP(42,K) <= 24179).OR.    &
     &         (LOOKUP(42,K) >= 24225.AND.LOOKUP(42,K) <= 24227))THEN
              LOOKUP(45,K)=3        !Slab diagnostic

            ELSE
              LOOKUP(45,K)=1        !Atmosphere diagnostic

            END IF

          ELSE IF(LOOKUP(42,K) >= 30000.AND.LOOKUP(42,K) <= 99999)THEN
            IF(LOOKUP(42,K) >= 40000.AND.LOOKUP(42,K) <= 40999)THEN
              LOOKUP(45,K)=3        !Slab diagnostic

            ELSE
              LOOKUP(45,K)=2        !Ocean diagnostic

            END IF

          ELSE
       If (PrintStatus >= PrStatus_Normal) Then
            WRITE(6,*) 'WARNING: User defined field found - ',          &
     &                 'STASH code : ', LOOKUP(42,K)
            WRITE(6,*) ' Internal model number can not be defined.'
            WRITE(6,*) ' Setting internal model number to atmosphere.'
       End If
            LOOKUP(45,K)=1

          ENDIF

        ENDIF

      ENDDO
!---------------------------------------------------------------
!  Reset LOOKUP headers if dump created earlier than vn2.8
!---------------------------------------------------------------

      IF(FIXHD(12) <  208)THEN
! DEPENDS ON: newpack
        CALL NEWPACK(LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP)
      ENDIF

! Check LOOKUP for consistency with PARAMETER statements
      IF(LOOKUP(LBNREC,1) == 0 .OR.                                     &
!        Prog lookups in dump before vn3.2:
     &  (LOOKUP(LBNREC,1) == IMDI .AND. FIXHD(12) <= 301)) THEN
        IF(LEN_DATA /= IMDI)THEN
! DEPENDS ON: chk_look
      CALL CHK_LOOK(FIXHD,LOOKUP,LEN1_LOOKUP,LEN_DATA,                  &
#include "argppx.h"
     &              ICODE,CMESSAGE)
        ENDIF
      ENDIF

      ENDIF

      RETURN
      END SUBROUTINE READHEAD
#endif
