#if defined(C80_1A) || defined(UTILIO) || defined(RECON)               \
 || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE WRITHEAD---------------------------------------
!LL
!LL  Purpose: Writes out model dump header records on unit NFTOUT &
!LL           checks model and dump dimensions for consistency.
!LL           32-bit IEEE output option supported
!LL           64-bit IEEE output option supported
!LL
!LL  Written by A. Dickinson 31/01/90
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.5    28/03/95 MPP code: New code for parallel I/O
!LL                                              P.Burton
!     4.1    18/06/96 Changes to cope with changes in STASH addressing
!                     Author D.M. Goddard.
!    4.1  21/05/96  Correct conversion of LOOKUP(65-128) in obs
!                   files for IEEE/64 bits. D Robinson
!     4.4    17/07/97 Introduce conversion from ieee to Cray PVP
!                     numbers and reintroduce functionality for
!                     PVP machines
!                     Author: D.M. Goddard
!     4.5    25/08/98 Correct conversion of LOOKUP(65-128) for
!                     cx files and OPS obstores.
!                     Author D.M. Goddard
!     5.1    10/04/00 New reconfiguration and stdout control using
!                     PrintStatus. P.Selwood.
!     5.3    22/11/01 Enable MPP as the only option for
!                     small executables         E.Leung
!     6.0    11/09/03 Replaced ABORT with call to EREPORT.    P.Dando
!     6.0    20/08/03 Porting to NEC  E.Leung
!     6.2    28/10/05 Fixes for -64e option and remove obsolete code.
!                     P.Selwood
!     6.2    06/12/05 Removed the DIAG81 CPP DEF. T.Edwards
!     6.2    15/08/05 Minor reworkings for FCM. P.Selwood
!LL
!LL  Programming standard:
!LL           Unified Model Documentation Paper No 3
!LL           Version No 1 15/1/90
!LL
!LL  System component: W30
!LL
!LL  System task: F3
!LL
!LL  Documentation:
!LL           Unified Model Documentation Paper No F3
!LL           Version No 5 9/2/90
!LL
!LL------------------------------------------------------------
!*L Arguments:-------------------------------------------------
      SUBROUTINE WRITHEAD(NFTOUT,FIXHD,LEN_FIXHD,                       &
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
#if defined(IEEE)
     &                    IEEE_TYPE,                                    &
#endif
#include "argppx.h"
     &                    START_BLOCK,ICODE,CMESSAGE)   ! Intent (Out)

#if defined(RECON) || defined(VAROPSVER)
      Use Rcf_Parvars_Mod, Only :                                       &
     &    mype

      Use Rcf_PrintStatus_Mod, Only :                                   &
     &    PrintStatus,                                                  &
     &    PrStatus_Oper

      Use Ereport_Mod, Only :                                           &
     &    Ereport
#endif

      IMPLICIT NONE

      INTEGER                                                           &
     & NFTOUT                                                           &
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
     &,START_BLOCK    !OUT Pointer to position of each block.
                      !Should point to start of model data block on exit

      INTEGER ::  ICODE      !OUT Return code; successful=0
                             !                 error > 0

#if defined(IEEE)
      INTEGER ::  IEEE_TYPE  !IN IEEE precision
#endif

#if defined(RECON)
! Dummy variables from arguments
      Integer ppxi
      Integer ppxrecs
      Character ppxc
#endif

      CHARACTER*(80)                                                    &
     & CMESSAGE       !OUT Error message if ICODE > 0

      CHARACTER (Len=*), Parameter  ::  RoutineName = 'WRITHEAD'

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
#if defined(IEEE)
      INTEGER,EXTERNAL :: IEEE2IEG

      INTEGER                                                           &
     & FIXHD_OUT(LEN_FIXHD)                                             &
                            ! Fixed length header
     &,INTHD_OUT(LEN_INTHD)                                             &
                            ! Integer header
     &,LOOKUP_O32(LEN1_LOOKUP/2,LEN2_LOOKUP)                            &
                                             !PP lookup tables (32-bit)
     &,LOOKUP_O64(LEN1_LOOKUP,LEN2_LOOKUP)                              &
                                             !PP lookup tables (64-bit)
     &,CFI1_OUT(LEN_CFI1+1)                                             &
                            ! Compressed field index no 1
     &,CFI2_OUT(LEN_CFI2+1)                                             &
                            ! Compressed field index no 2
     &,CFI3_OUT(LEN_CFI3+1)                                             &
                            ! Compressed field index no 3
     &,I,J

      REAL                                                              &
     & REALHD_OUT(LEN_REALHD)                                           &
                              ! Real header
     &,LEVDEPC_OUT(1+LEN1_LEVDEPC*LEN2_LEVDEPC)                         &
                                                ! Lev dep consts
     &,ROWDEPC_OUT(1+LEN1_ROWDEPC*LEN2_ROWDEPC)                         &
                                                ! Row dep consts
     &,COLDEPC_OUT(1+LEN1_COLDEPC*LEN2_COLDEPC)                         &
                                                ! Col dep consts
     &,FLDDEPC_OUT(1+LEN1_FLDDEPC*LEN2_FLDDEPC)                         &
                                                ! Field dep consts
     &,EXTCNST_OUT(LEN_EXTCNST+1)                                       &
                                    ! Extra constants
     &,DUMPHIST_OUT(LEN_DUMPHIST+1) ! History block

#endif
#if !defined(RECON)  && !defined(VAROPSVER)
#include "parvars.h"
#endif
! -------------------------------------------------------------
! External subroutines called:---------------------------------
      EXTERNAL IOERROR,POSERROR,PR_FIXHD,BUFFOUT
!*-------------------------------------------------------------
! Local variables:---------------------------------------------
      INTEGER LEN_IO
      REAL A
! -------------------------------------------------------------
! Comdecks:----------------------------------------------------------
#if !defined(RECON)  && !defined(VAROPSVER)
#include "cprintst.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#endif
! -------------------------------------------------------------
      ICODE=0
      CMESSAGE=' '

!L 1. Buffer out fixed length header record

#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(2,LEN_FIXHD,FIXHD_OUT,0,FIXHD,1,64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,FIXHD_OUT,LEN_FIXHD,LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,FIXHD(1),LEN_FIXHD,LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,FIXHD(1),LEN_FIXHD,LEN_IO,A)
#endif


! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= LEN_FIXHD)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of fixed length header',A,LEN_IO       &
     &               ,LEN_FIXHD)
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=1
        RETURN
      ENDIF

      START_BLOCK=LEN_FIXHD+1

! Check validity of data and print out fixed header information
      IF(PrintStatus >= PrStatus_Oper) THEN

! DEPENDS ON: pr_fixhd
      CALL PR_FIXHD(FIXHD,LEN_FIXHD,LEN_INTHD,LEN_REALHD,LEN1_LEVDEPC   &
     &,LEN2_LEVDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC &
     &,LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,LEN_DUMPHIST,LEN_CFI1      &
     &,LEN_CFI2,LEN_CFI3,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA               &
     &,ICODE,CMESSAGE)

       IF(ICODE >  0) THEN
        CMESSAGE="WRITHEAD: Error returned from PR_FIXHD"
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       END IF
      END IF

!L 2. Buffer out integer constants

      IF(FIXHD(100) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(100) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('integer constants',START_BLOCK,100,FIXHD(100))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=2
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(2,FIXHD(101),INTHD_OUT,0,INTHD,1,64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,INTHD_OUT,FIXHD(101),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,INTHD(1),FIXHD(101),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,INTHD(1),FIXHD(101),LEN_IO,A)
#endif

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(101))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of integer constants',A,LEN_IO,        &
     &               FIXHD(101))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=3
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(101)

      ENDIF

!L 3. Buffer out real constants

      IF(FIXHD(105) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(105) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('real constants',START_BLOCK,105,FIXHD(105))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=4
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

! Check for I/O errors
#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(3,FIXHD(106),REALHD_OUT,0,REALHD,1,64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,REALHD_OUT,FIXHD(106),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
        REALHD_OUT(:) = REALHD(:)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,REALHD_OUT,FIXHD(106),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,REALHD(1),FIXHD(106),LEN_IO,A)
#endif

       IF(A /= -1.0.OR.LEN_IO /= FIXHD(106))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of real constants',A,LEN_IO,           &
     &                FIXHD(106))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=5
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(106)

      ENDIF

!L 4. Buffer out level dependent constants

      IF(FIXHD(110) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(110) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('level dependent constants',                      &
     &  START_BLOCK,110,FIXHD(110))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=6
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(3,FIXHD(111)*FIXHD(112),LEVDEPC_OUT,0,              &
     &             LEVDEPC,1,64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,LEVDEPC_OUT,FIXHD(111)*FIXHD(112),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
        LEVDEPC_OUT(:) = LEVDEPC(:)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,LEVDEPC_OUT,FIXHD(111)*FIXHD(112),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,LEVDEPC(1),FIXHD(111)*FIXHD(112),LEN_IO,A)
#endif

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(111)*FIXHD(112))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of level dependent constants',A,LEN_IO,&
     &               FIXHD(111)*FIXHD(112))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=7
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(111)*FIXHD(112)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' LEVEL DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(111)*FIXHD(112)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 5. Buffer out row dependent constants

      IF(FIXHD(115) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(115) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('row dependent constants',                        &
     &  START_BLOCK,115,FIXHD(115))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=8
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(3,FIXHD(116)*FIXHD(117),ROWDEPC_OUT,0,              &
     &             ROWDEPC,1,64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,ROWDEPC_OUT,FIXHD(116)*FIXHD(117),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
        ROWDEPC_OUT(:) = ROWDEPC(:)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,ROWDEPC_OUT,FIXHD(116)*FIXHD(117),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,ROWDEPC(1),FIXHD(116)*FIXHD(117),LEN_IO,A)
#endif

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(116)*FIXHD(117))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of row dependent constants',A,LEN_IO,  &
     &                FIXHD(116)*FIXHD(117))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=9
        RETURN
      ENDIF


       START_BLOCK=START_BLOCK+FIXHD(116)*FIXHD(117)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' ROW DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(116)*FIXHD(117)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 6. Buffer out column dependent constants

      IF(FIXHD(120) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(120) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('column dependent constants',                     &
     &  START_BLOCK,120,FIXHD(120))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=10
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(3,FIXHD(121)*FIXHD(122),COLDEPC_OUT,0,              &
     &             COLDEPC,1,64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,COLDEPC_OUT,FIXHD(121)*FIXHD(122),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
        COLDEPC_OUT(:) = COLDEPC(:)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,COLDEPC_OUT,FIXHD(121)*FIXHD(122),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,COLDEPC(1),FIXHD(121)*FIXHD(122),LEN_IO,A)
#endif

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(121)*FIXHD(122))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of column dependent constants',A,LEN_IO&
     &              ,FIXHD(121)*FIXHD(122))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=11
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(121)*FIXHD(122)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COLUMN DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(121)*FIXHD(122)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 7. Buffer out constants stored as fields

      IF(FIXHD(125) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(125) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('fields of constants',                            &
     &  START_BLOCK,125,FIXHD(125))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=12
        RETURN
       ENDIF

#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(3,FIXHD(126)*FIXHD(127),FLDDEPC_OUT,0,              &
     &             FLDDEPC,1,64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,FLDDEPC_OUT,FIXHD(126)*FIXHD(127),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
        FLDDEPC_OUT(:) = FLDDEPC(:)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,FLDDEPC_OUT,FIXHD(126)*FIXHD(127),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,FLDDEPC(1),FIXHD(126)*FIXHD(127),LEN_IO,A)
#endif

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(126)*FIXHD(127))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of field dependent constants',A,LEN_IO,&
     &               FIXHD(126)*FIXHD(127))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=13
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(126)*FIXHD(127)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' FIELD DEPENDENT CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(126)*FIXHD(127)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 8. Buffer out extra constants

      IF(FIXHD(130) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(130) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('extra constants',                                &
     &  START_BLOCK,130,FIXHD(130))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=14
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(3,FIXHD(131),EXTCNST_OUT,0,EXTCNST,1,               &
     &             64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,EXTCNST_OUT,FIXHD(131),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
        EXTCNST_OUT(:) = EXTCNST(:)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,EXTCNST_OUT,FIXHD(131),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,EXTCNST(1),FIXHD(131),LEN_IO,A)
#endif

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(131))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out extra constants',A,LEN_IO,             &
     &               FIXHD(131))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=15
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(131)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' EXTRA CONSTANTS'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(131)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 9. Buffer out temporary history block

      IF(FIXHD(135) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(135) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('history',                                        &
     &  START_BLOCK,136,FIXHD(136))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=16
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(3,FIXHD(136),DUMPHIST_OUT,0,DUMPHIST,1,             &
     &             64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,DUMPHIST_OUT,FIXHD(136),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
        DUMPHIST_OUT(:) = DUMPHIST(:)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,DUMPHIST_OUT,FIXHD(136),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,DUMPHIST(1),FIXHD(136),LEN_IO,A)
#endif

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(136))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of history file',A,LEN_IO,             &
     &               FIXHD(136))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=17
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(136)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' TEMPORARY HISTORY BLOCK'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(136)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 10. Buffer out compressed field index1

      IF(FIXHD(140) >  0)THEN

! Check for error in file pointers

       IF(FIXHD(140) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('compressed field index1',                        &
     &  START_BLOCK,140,FIXHD(140))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=18
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF


#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(2,FIXHD(141),CFI1_OUT,0,CFI1,1,64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,CFI1_OUT,FIXHD(141),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,CFI1(1),FIXHD(141),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,CFI1(1),FIXHD(141),LEN_IO,A)
#endif

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(141))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of compressed index1',A,LEN_IO,        &
     &               FIXHD(141))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=19
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(141)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 1'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(141)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 11. Buffer out compressed field index2

      IF(FIXHD(142) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(142) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('compressed field index2',                        &
     &  START_BLOCK,142,FIXHD(142))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=20
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(2,FIXHD(143),CFI2_OUT,0,CFI2,1,64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,CFI2_OUT,FIXHD(143),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,CFI2(1),FIXHD(143),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,CFI2(1),FIXHD(143),LEN_IO,A)
#endif

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(143))THEN
! DEPENDS ON: ioerror
       CALL IOERROR('buffer out of compressed index2',A,LEN_IO,         &
     &               FIXHD(143))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=21
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(143)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 2'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(143)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 12. Buffer out compressed field index3

      IF(FIXHD(144) >  0)THEN

! Check for error in file pointers
       IF(FIXHD(144) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('compressed field index3',                        &
     &  START_BLOCK,144,FIXHD(144))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=22
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        I= IEEE2IEG(2,FIXHD(145),CFI3_OUT,0,CFI3,1,64,IEEE_TYPE)
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,CFI3_OUT,FIXHD(145),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,CFI3(1),FIXHD(145),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,CFI3(1),FIXHD(145),LEN_IO,A)
#endif

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(145))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of compressed index3',A,LEN_IO,        &
     &               FIXHD(145))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=23
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(145)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' COMPRESSED FIELD INDEX NO 3'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(145)
      ENDIF ! if mype  ==  0
      End If

      ENDIF

!L 13. Buffer out lookup table

      IF(FIXHD(150) >  0)THEN
#if defined(IEEE)
!
      if(start_block /= fixhd(150)) then
        if(start_block >  fixhd(150)) then
          write(6,9975) start_block-1, fixhd(150)-1
9975      format(/                                                      &
     &     10(/'**** ERROR - Current Disk Address is greater than',     &
     &     ' the Address in the Fixed Length Header for the',           &
     &     ' Lookup Table *****'))
          CMessage = 'CONVIEEE: Fixed length Header Error'
          ICODE = 26
          call EReport(RoutineName,ICODE,CMessage)
        else
          write(6,9976) start_block-1, fixhd(150)-1
9976      format(                                                       &
     &     10(/'**** WARNING - Current Disk Address does not match',    &
     &     ' the Address in the Fixed Length Header for the',           &
     &     ' Lookup Table *****')//                                     &
     &     'Current Address altered from ',i10,' to ',i10,              &
     &     ' to match the Fixed Length Header'/)
          start_block=fixhd(150)
          if(ieee_type == 32) then
            call setpos32(nftout, start_block-1, j)
          else
! DEPENDS ON: setpos
            call setpos(nftout, start_block-1, j)
          endif
        endif
      endif
#endif

! Check for error in file pointers
       IF(FIXHD(150) /= START_BLOCK)THEN
! DEPENDS ON: poserror
        CALL POSERROR('lookup table',                                   &
     &  START_BLOCK,150,FIXHD(150))
        CMESSAGE='WRITHEAD: Addressing conflict'
        ICODE=24
        CALL EREPORT("WRITHEAD", ICODE, CMESSAGE)
       ENDIF

#if defined(IEEE)
      IF(IEEE_TYPE == 32)THEN
        DO I=1,FIXHD(152)
        J= IEEE2IEG(2,45,LOOKUP_O32(1,I),0,LOOKUP(1,I),1,               &
     &             64,IEEE_TYPE)
        J= IEEE2IEG(3,19,LOOKUP_O32(23,I),32,LOOKUP(46,I),1,            &
     &             64,IEEE_TYPE)
        IF (FIXHD(5) == 6.OR.FIXHD(5) == 7.OR.                          &
                                                 ! 6=ACOBS 7=VAROBS
     &      FIXHD(5) == 8.OR.FIXHD(5) == 10)THEN ! 8=CX   10=OBSTORE
          J= IEEE2IEG(2,64,LOOKUP_O32(33,I),0,LOOKUP(65,I),             &
     &              1,64,IEEE_TYPE)
        ENDIF
        ENDDO
! DEPENDS ON: buffo32
        CALL BUFFO32(NFTOUT,LOOKUP_O32,FIXHD(151)*FIXHD(152),LEN_IO,A)
      ELSEIF(IEEE_TYPE == 64)THEN
        LOOKUP_O64(:,:) = LOOKUP(:,:)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,LOOKUP_O64,FIXHD(151)*FIXHD(152),LEN_IO,A)
      ENDIF
#else
! DEPENDS ON: buffout
      CALL BUFFOUT(NFTOUT,LOOKUP(1,1),FIXHD(151)*FIXHD(152),LEN_IO,A)
#endif

! Check for I/O errors
       IF(A /= -1.0.OR.LEN_IO /= FIXHD(151)*FIXHD(152))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer out of lookup table',A,LEN_IO,             &
     &               FIXHD(151)*FIXHD(152))
        CMESSAGE='WRITHEAD: I/O error'
        ICODE=25
        RETURN
       ENDIF

       START_BLOCK=START_BLOCK+FIXHD(151)*FIXHD(152)

      If (PrintStatus >= PrStatus_Oper) Then
      IF (mype  ==  0) THEN
       WRITE(6,'('' '')')
       WRITE(6,'('' LOOKUP TABLE'')')
       WRITE(6,'('' '',I8,'' 64-bit words long'')')FIXHD(151)*FIXHD(152)
      ENDIF ! if mype  ==  0
      End If

! No consistency checks for parallel code. The LOOKUP headers don't
! match the data layout in memory within a PE.

      ENDIF

      RETURN
      END SUBROUTINE WRITHEAD
#endif
