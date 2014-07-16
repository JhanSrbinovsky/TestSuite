#if defined(FLDMOD)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: FLDMOD --------------------------------------------------
!LL
!LL  Purpose: To read a   direct access PP file  and convert it to a
!LL  sequential file read to be passed across to the IBM
!LL
!LL  Modification History:
!LL  Copy of FIELDMOD taken, named UMTHIN1 and code for thinning
!LL  fields added. I/O routines from portable model included so that
!LL  file names are passed through environmental variables instead
!LL  of assigns. Vic Blackman January-March 1995
!LL
!LL  Calls to OPEN replaced by FILE_OPEN due to change in UM v3.5
!LL  4.2   29/11/96 (1)Corrections to code to enable qxumthin1 to
!LL         produce a bit comparable dump to that produced by qxfieldmod
!LL                 (as used in operational suite) using the namelist
!LL                 /u/opfc/op/perm/in/qifieldmod as input.
!LL                 (2)Renamed FLDMOD
!LL                 Author I.Edmond
!LL  4.3   15/4/97 Get the current sector size for disk I/O from
!LL                environment variable UM_SECTOR_SIZE in calling
!LL                script, otherwise, use UM_SECTOR_SIZE=512.
!LL                Required by INITPP to make sure the data starts
!LL                on a sector bndry. IEdmond
!    4.5   05/06/98 Prevent failure if lookup table starts on a
!                   sector boundary
!                   Author D.M. Goddard
!LL  5.2  7/11/99  Enable run length encoding of ocean fieldsfiles to
!LL                compress the sequences of mdi values that represent
!LL                the land points. Ian Edmond
!    5.2   15/12/00 (1)Fixing problem with different array size in real
!                      and integer header array between vn4.x & vn5.1
!                   (2)Initialize icode in subroutine THIN_FIELD,
!                      SCALE_FIELD, WIND_10_M
!                    E.Leung
!    5.5   25/04/03 (a)Remove UNIT function for portability.
!                   (b)Rename FLDMOD to FLD_MOD.
!                   (c)Grib data format not support on non-CRAY platform
!                   (d)General tidy up.                      E.Leung
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Project task: ...
!LL
!LL  External documentation: On-line UM document ??? - ??????????
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!LL  Routine: CONTROL----------------------------------------------
!LL
!LL  Purpose: To control the calculation of the derived diagnostics
!LL  and output of the new LOOKUP table (called LOOKNEW)
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Modification History:
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Project task: ...
!LL
!LL  External documentation: On-line UM document ??? - ??????????
!LL
!LL  ---------
!*L  Interface and arguments: ------------------------------------------
!

!LL  Routine: DIMENS1--------------------------------------------
!LL
!LL  Purpose: To read a   direct access PP file  and convert it to a
!LL  sequential file read to be passed across to the IBM
!LL
!LL  Modification History:
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!

!LL  Routine: FIELDS ----------------------------------------------
!LL
!LL  Purpose: To calculate fields from the Fields File such as those
!LL  normaly derived in the Derived Printfile Program
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.190
!LL
!LL  Modification History:
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Project task: ...
!LL
!LL  External documentation: On-line UM document ??? - ??????????
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!


!LL  SUBROUTINE FLD_OUT------------------------------------------
!LL
!LL  REPLACES THE OUTPUT FROM STASH EITHER ON TO A PP FILE OR
!LL  BACK TO THE MAIN ARRAY D1
!LL
!LL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
!LL  VERSION 1, DATED 12/09/89
!LL
!LL  SYSTEM TASK: CONTROL PART OF C4
!LL
!LL  PURPOSE:   TO PROCESS DIAGNOSTICS CONTROLLED BY STASH
!LL  DOCUMENTATION:        ???
!LL
!LL
!LLEND-------------------------------------------------------------

!
!*L  ARGUMENTS:---------------------------------------------------
!LL  Routine: READPP--------------------------------------------
!LL
!LL  Purpose: To read a   direct access PP file  and convert it to a
!LL  sequential file read to be passed across to the IBM
!LL
!LL  Modification History:
!LL   5.1     31/03/00   Removed reference to READPP in the EXTERNAL
!LL                      statement (non-standard Fortran). D.P.Matthews
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE READPP(LEN_INTHD,LEN_REALHD,LEN1_LEVDPC,LEN2_LEVDPC,   &
     &   LEN1_ROWDPC,LEN2_ROWDPC,LEN1_COLDPC,LEN2_COLDPC,               &
     &   LEN1_LOOKUP,LEN2_LOOKUP,LEN_FIXHD,PP_FIXHD,LOOKUP,ROOKUP,      &
     &   PPUNIT1,PPUNIT2,ICODE,CMESSAGE)
      IMPLICIT NONE
      EXTERNAL POSERROR,IOERROR,GETPOS,SETPOS
      INTEGER                                                           &
     &     LEN_FIXHD                                                    &
     &    ,LEN_INTHD                                                    &
     &    ,LEN_REALHD                                                   &
     &    ,LEN_LEVDPC                                                   &
     &    ,LEN_ROWDPC                                                   &
     &    ,LEN_COLDPC                                                   &
     &    ,LEN1_LEVDPC                                                  &
     &    ,LEN2_LEVDPC                                                  &
     &    ,LEN1_ROWDPC                                                  &
     &    ,LEN2_ROWDPC                                                  &
     &    ,LEN1_COLDPC                                                  &
     &    ,LEN2_COLDPC                                                  &
     &    ,LEN_LOOKUP                                                   &
     &    ,LEN1_LOOKUP                                                  &
     &    ,LEN2_LOOKUP                                                  &
     &    ,LEN1_LOOKNEW                                                 &
     &    ,LEN2_LOOKNEW                                                 &
     &    ,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                              &
     &    ,PP_INTHD(LEN_INTHD)                                          &
     &    ,PP_FIXHD(LEN_FIXHD)                                          &
     &    ,LEN_IO                                                       &
     &    ,ICODE                                                        &
     &    ,PPUNIT1                                                      &
     &    ,PPUNIT2
      REAL                                                              &
     &     ROOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                              &
     &    ,PP_REALHD(LEN_REALHD)                                        &
     &    ,PP_LEVDPC(LEN1_LEVDPC*LEN2_LEVDPC+1)                         &
     &    ,PP_ROWDPC(LEN1_ROWDPC*LEN2_ROWDPC)                          &
     &    ,PP_COLDPC(LEN1_COLDPC*LEN2_COLDPC)                          &
     &    ,A_IO
      CHARACTER CMESSAGE*(*)
      CHARACTER OUTFILE*80
! Local variables
      INTEGER                                                           &
     &     START_BLOCK                                                  &
     &    ,NENT                                                         &
     &    ,K                                                            &
     &    ,Kk                                                           &
     &    ,iwa                                                          &
     &    ,RECL                                                         &
     &    ,IERR                                                         &
     &    ,PP_LEN_INTHD                                                 &
     &    ,PP_LEN_REALHD                                                 

!---------------------------------------------------------------------
      PP_LEN_REALHD=PP_FIXHD(106)
      PP_LEN_INTHD=PP_FIXHD(101)

!---------------------------------------------------------------------
      LEN_LEVDPC=LEN1_LEVDPC*LEN2_LEVDPC
      LEN_ROWDPC=LEN1_ROWDPC*LEN2_ROWDPC
      LEN_COLDPC=LEN1_COLDPC*LEN2_COLDPC
      LEN_LOOKUP=LEN1_LOOKUP*LEN2_LOOKUP
! The calculation of LEN_LEVDPC has PLUS 1 which is only true
! for PP headers and not model headers, hopefully the PLUS one will
! be removed as it is inconsistent)
      START_BLOCK=LEN_FIXHD+1
!L---------------------------------------------------------------
!L  Read in the integer constants
!L---------------------------------------------------------------
      IF(LEN_INTHD >  0) THEN  ! Integer constants to be read in
        IF(PP_FIXHD(100) /= START_BLOCK) THEN   ! Address incorrect
! DEPENDS ON: poserror
          CALL POSERROR('integer constants',START_BLOCK,100,            &
     &    PP_FIXHD(100))
          CMESSAGE=' READPP :  Adressing Conflict'
          ICODE=2
          RETURN
        ENDIF
! DEPENDS ON: buffin
        CALL BUFFIN(PPUNIT1,PP_INTHD,LEN_INTHD,LEN_IO,A_IO)
        WRITE(6,*)pp_inthd
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_INTHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR(' Buffer in of Integer constants',A_IO,LEN_IO    &
     &  ,  LEN_INTHD)
          CMESSAGE='READPP : I/O error'
          ICODE=3
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_INTHD
      ENDIF
!L---------------------------------------------------------------
!L  Read in the real constants
!L---------------------------------------------------------------
      IF(LEN_REALHD >  0) THEN  ! Real constants to be read in
        IF(PP_FIXHD(105) /= START_BLOCK) THEN   ! Address incorrect
! DEPENDS ON: poserror
          CALL POSERROR('Real constants',START_BLOCK,100,               &
     &    PP_FIXHD(105))
          CMESSAGE=' READPP :  Adressing Conflict'
          ICODE=4
          RETURN
        ENDIF
! DEPENDS ON: buffin
        CALL BUFFIN(PPUNIT1,PP_REALHD,LEN_REALHD,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_REALHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR(' Buffer in of Real constants',A_IO,LEN_IO       &
     &    ,LEN_REALHD)
          CMESSAGE='READPP : I/O error'
          ICODE=5
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_REALHD
      ENDIF
!L---------------------------------------------------------------
!L  Read in the level dependant constants
!L---------------------------------------------------------------
      IF(LEN_LEVDPC >  0) THEN  ! Level dep constants to be read in
        IF(PP_FIXHD(110) /= START_BLOCK) THEN   ! Address incorrect
! DEPENDS ON: poserror
          CALL POSERROR('Level depndt constants',START_BLOCK,100,       &
     &    PP_FIXHD(110))
          CMESSAGE=' READPP :  Adressing Conflict'
          ICODE=6
          RETURN
        ENDIF
! DEPENDS ON: buffin
        CALL BUFFIN(PPUNIT1,PP_LEVDPC,LEN_LEVDPC,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_LEVDPC) THEN
! DEPENDS ON: ioerror
          CALL IOERROR(' Buffer in of Level constants',A_IO,LEN_IO      &
     &    ,LEN_LEVDPC)
          CMESSAGE='READPP : I/O error'
          ICODE=7
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_LEVDPC
      ENDIF
!L---------------------------------------------------------------
!L  Read in the Row dependant constants
!L---------------------------------------------------------------
      IF(LEN_ROWDPC >  0) THEN  ! row dep constants to be read in
        IF(PP_FIXHD(115) /= START_BLOCK) THEN   ! Address incorrect
! DEPENDS ON: poserror
          CALL POSERROR('Row depndt constants',START_BLOCK,100,         &
     &                   PP_FIXHD(115))
          CMESSAGE=' READPP :  Adressing Conflict'
          ICODE=10
          RETURN
        ENDIF
! DEPENDS ON: buffin
        CALL BUFFIN(PPUNIT1,PP_ROWDPC,LEN_ROWDPC,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_ROWDPC) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('Buffer in of Row constants',A_IO,LEN_IO,        &
     &                  LEN_ROWDPC)
          CMESSAGE='READPP : I/O error'
          ICODE=11
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_ROWDPC
      ENDIF
!L---------------------------------------------------------------
!L  Read in the Col dependant constants
!L---------------------------------------------------------------
      IF(LEN_COLDPC >  0) THEN  ! col dep constants to be read in
        IF(PP_FIXHD(120) /= START_BLOCK) THEN   ! Address incorrect
! DEPENDS ON: poserror
          CALL POSERROR('Col depndt constants',START_BLOCK,100,         &
     &                   PP_FIXHD(120))
          CMESSAGE=' READPP :  Adressing Conflict'
          ICODE=20
          RETURN
        ENDIF
! DEPENDS ON: buffin
        CALL BUFFIN(PPUNIT1,PP_COLDPC,LEN_COLDPC,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_ROWDPC) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('Buffer in of Col constants',A_IO,LEN_IO,        &
     &                  LEN_COLDPC)
          CMESSAGE='READPP : I/O error'
          ICODE=21
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_COLDPC
      ENDIF      
!L---------------------------------------------------------------
!L  Read in the LOOKUP TABLE
!L---------------------------------------------------------------
      IF(LEN_LOOKUP >  0) THEN  ! Lookup Table to be read in
        IF(PP_FIXHD(150) /= START_BLOCK) THEN   ! Address incorrect
          WRITE(6,*) 'READPP : WARNING'
          WRITE(6,*) 'Conflict between start position of Lookup table'
          WRITE(6,*) 'block and pointer in fixed length header: ',      &
     &               'FIXHD(150) = ',PP_FIXHD(150)
          WRITE(6,*) 'Current position in file = ',START_BLOCK,         &
     &               ' words in'
          WRITE(6,*) 'Pointer moved to ',PP_FIXHD(150),' words in'
! DEPENDS ON: setpos
          CALL SETPOS(PPUNIT1,PP_FIXHD(150)-1,IERR)
        END IF
! DEPENDS ON: buffin
        CALL BUFFIN(PPUNIT1,LOOKUP,LEN_LOOKUP,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_LOOKUP) THEN
! DEPENDS ON: ioerror
          CALL IOERROR(' Buffer in of Lookup table   ',A_IO,LEN_IO      &
     &    ,LEN_LOOKUP)
          CMESSAGE='READPP : I/O error'
          ICODE=9
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_LOOKUP
      ENDIF
      WRITE(6,*)' ARRIVED HERE  ',START_BLOCK
      NENT=0
      DO 1 K=1,LEN2_LOOKUP
        IF(LOOKUP(1,K) >  0) THEN
          NENT=NENT+1
        ELSE
          GOTO 2
        ENDIF
    1 CONTINUE
    2 CONTINUE
      WRITE(6,*)' VALUE OF NENT   ',NENT
      do k=nent-2,nent+1
        WRITE(6,*)'k=',k
        WRITE(6,*) (lookup(kk,k),kk=1,44)
      enddo
!-----------------------------------------------------------------
!    OPEN NEW TARGET FIELDSFILE INITIALISING BY CALLING INITPP
!-----------------------------------------------------------------
!L
!L        Open named file on unit 60
!L
        WRITE(6,*)"*** Opening new file on unit ",pPUNIT2
        CALL GET_FILE(PPUNIT2,OUTFILE,80,ICODE)
! DEPENDS ON: file_open
        CALL FILE_OPEN(PPUNIT2,OUTFILE,80,1,1,ICODE)
!
! DEPENDS ON: init_pp
      CALL INIT_PP(PPUNIT2,'p',LEN1_LOOKUP,LEN2_LOOKUP,PP_FIXHD,        &
     &             PP_INTHD,PP_REALHD,PP_LEVDPC,PP_ROWDPC,PP_COLDPC,    &
     &             LEN_FIXHD,LEN_INTHD,LEN_REALHD,LEN1_LEVDPC,          &
     &             LEN2_LEVDPC,LEN1_ROWDPC,LEN2_ROWDPC,                 &
     &             LEN1_COLDPC,LEN2_COLDPC,PP_LEN_INTHD,                &
     &             PP_LEN_REALHD,ICODE,CMESSAGE)

      IF(ICODE /= 0) THEN
        WRITE(7,100) ICODE
        WRITE(7,110) CMESSAGE
        RETURN
 100  FORMAT(' ICODE EQUAL TO ',I2)
 110  FORMAT(A80)
      ENDIF
      LEN1_LOOKNEW=LEN1_LOOKUP
      LEN2_LOOKNEW=LEN2_LOOKUP
! DEPENDS ON: control
      CALL CONTROL(PPUNIT1,PPUNIT2,LEN1_LOOKNEW,LEN2_LOOKNEW,           &
     &             LOOKUP,PP_INTHD,LEN_INTHD,                           &
     &             PP_FIXHD,LEN_FIXHD,ICODE,CMESSAGE,NENT)
      IF(ICODE /= 0) THEN
        WRITE(7,120) ICODE
        WRITE(7,130) CMESSAGE
        RETURN
 120  FORMAT(' ICODE EQUAL TO ',I2)
 130  FORMAT(A80)
      ENDIF
      GOTO 901
  900 CONTINUE
      WRITE(6,*)' ERROR IN READPP OPENING THE PPUNIT2 FIELDS FILE'
  901 CONTINUE
 9999 CONTINUE
      RETURN
      END SUBROUTINE READPP
!LL  Routine: RE_PACK  -------------------------------------------------
!LL
!LL  Purpose: To repack data from the input array FIELD and return
!LL
!LL  Model            Modification history:
!LL version  Date
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
#endif
