#if defined(C98_1A) || defined(FLDMOD)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: FFREAD and others (see below) ----------------------------
!LL
!LL  Purpose: To read a   direct access PP file  and convert it to a
!LL  sequential file read to be passed across to the IBM
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  19/02/93  Use FIXHD(12) not FIXHD(1) as Version no in P21BITS
!LL  3.1     22/01/93 revert to using Standard CRAY BUFFER IN
!LL  3.2     08/07/93 code to use MDI from pp-header
!LL  3.4     09/09/94 add GRIB decode interface
!LL  4.1     12/12/96 Int lookup dimensioned by 45 (real part by 19)
!LL  5.1  7/11/99  Enable run length encoding of ocean fieldsfiles to
!LL                compress the sequences of mdi values that represent
!LL                the land points. Ian Edmond
!LL  5.5  25/04/03  Remove UNIT function for portability.
!LL                 Remove subroutine BUFFIN_207.  E.Leung
!LL  6.0  21/01/04  NECSX6 DEF added around call to DEGRIB as
!LL                 libgrib was then only available on NECSX6
!LL                 W Roseblade
!LL  6.2  15/08/05  Upped LEN1_RECORD to 30000 to handle larger files.
!LL                 Also added calls to EREPORT so that failures are
!LL                 recognized and registered. T. Edwards
!LL  6.2  16/08/05  Add a missing comma. P.Selwood
!    6.2  10/04/05  Removed calls to ABORT.  J. Gill
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!    IEXTRA(1) == 0  ! unpacking is required
!    IEXTRA(2) == 0  ! lookup table entry deleted after access.
!
      SUBROUTINE FFREAD(IPROJ,FCT,ITYPE,INT_LEVEL,PPUNIT,FIELD,IDIM,    &
     &ILABEL,RLABEL,IEXTRA,ICODE,CMESSAGE)
      IMPLICIT NONE
      EXTERNAL SETPOS,READ_REC,IOERROR,COEX

      CHARACTER CMESSAGE*(*)
      INTEGER                                                           &
     &     MAXFF                                                        &
                                  !OUT Max number of opened files
     &    ,LEN1_LOOKUP                                                  &
                                  !IN  first dimension of the lookup
     &    ,PP_LEN2_LOOKUP                                               &
                                  !OUT secnd dimension of the lookup
     &    ,LEN1_RECORD                                                  &
                                  !OUT First dimension of record
     &    ,RECLEN                                                       &
                                  !OUT Total length of record.
     &    ,IPROJ                                                        &
                                  !IN  map projection of data to read
     &    ,FCT                                                          &
                                  !IN  forecast period in hours
     &    ,ITYPE                                                        &
                                  !IN  M08 FIELD field type
     &    ,INT_LEVEL                                                    &
                                  !IN  LEVEL code (could be real)
     &    ,PPUNIT                                                       &
                                  !IN  unit no of required fieldsfile
     &    ,IDIM                                                         &
                                  !IN  dimension of FIELD
     &    ,ILABEL(45)                                                   &
                                  !OUT holds integer part of LOOKUP
     &    ,IEXTRA(10)                                                   &
                                  !INOUT Controls certain functions.
     &    ,ICODE                                                        &
                                  !OUT return code
     &    ,MAXPP                                                        &
                                  !    maximum number of unit number
     &    ,DATA_ADD                                                     &
                                  !OUT The word address of the data.
     &    ,PFNO                   !OUT No of fields files opened
      PARAMETER(MAXFF=10)
      PARAMETER(MAXPP=100)
      PARAMETER(LEN1_LOOKUP=64)
      PARAMETER(LEN1_RECORD=30000)!Max size of a lookup table allowed
      PARAMETER(RECLEN=LEN1_RECORD*MAXFF) ! Total length of RECORD
      REAL                                                              &
     &     FIELD(IDIM)                                                  &
                                  !OUT array holding final output data.
     &    ,RLABEL(19)                                                   &
                                  !OUT holds real part of LOOKUP
     &    ,REAL_LEVEL             !IN  LEVEL code (could be real)

!*-------------------------------------------------------------------
!     LOCAL VARIABLES
      INTEGER                                                           &
     &     TABLE(MAXPP)                                                 &
                                  ! associates unit no and file no
     &    ,PREV_PPUNIT(MAXFF)                                           &
                                  ! a record of unit nos already used
     &    ,I                                                            &
                                  ! local counter
     &    ,J                      ! local counter
      INTEGER                                                           &
     &     IX                                                           &
                                  ! used as a dummy variable in UNIT
     &    ,IWA                                                          &
                                  ! Word address in call SETPOS
     &    ,IK                                                           &
                                  ! Word address in call SETPOS
     &    ,ICOUNT                                                       &
                                  ! Counter
     &    ,LEN_IO                                                       &
                                  ! Length of data transferred from BUF
     &    ,LEN_FIXHD                                                    &
                                  ! Length of Fixed length header
     &    ,PP_FIXHD(256)                                                &
                                  ! Fixed length header
     &    ,IN_LBVC                ! Local copy of LBVC required to searc
      real                                                              &
     &     A_IO                   ! OUTPUT from UNIT command
      LOGICAL                                                           &
     &     TEST                                                         &
     &    ,RECORD(LEN1_RECORD,MAXFF)
      SAVE PREV_PPUNIT
      SAVE PFNO
      SAVE TABLE
      SAVE RECORD
!
      PARAMETER(LEN_FIXHD=256)
      DATA PREV_PPUNIT/MAXFF*0/
      DATA TABLE/MAXPP*0/
      DATA PFNO/0/
      DATA RECORD/RECLEN*.FALSE./
!    Remember that BUFFER OUT starts at address 0
      IF(PPUNIT >  MAXPP) THEN
        ICODE=1
        CMESSAGE=' FFREAD   the unit number is too large'
        RETURN
      ENDIF
      TEST=.TRUE.
!L Establish if a completely new FF is being read.
      DO 1 I=1,MAXFF
      IF(PPUNIT == PREV_PPUNIT(I)) THEN
        TEST=.FALSE.
      ENDIF
    1 CONTINUE
!L  A TABLE is set up associating FIELDS_FILE NO (1 to MAXFF) with a
!L  PP unit number. On succesive calls to the FF this table is used
!L  help record which LOOKUP table belongs to which FF (PPUNIT)
      IF(TEST) THEN   ! A FF never read in before.
        PFNO=PFNO+1
        PREV_PPUNIT(PFNO)=PPUNIT
        TABLE(PPUNIT)=PFNO
      ENDIF
!L  Read in the word address of the LOOKUP table (IWA), the length of
!L  the LOOKUP table (PP_LEN2_LOOKUP) and the start address of the data
!L  DATA_ADD
      PFNO=TABLE(PPUNIT)
      IWA=0
! DEPENDS ON: setpos
      CALL SETPOS(PPUNIT,IWA,ICODE) ! C coded routine
! DEPENDS ON: buffin
      CALL BUFFIN(PPUNIT,PP_FIXHD,LEN_FIXHD,LEN_IO,A_IO)
      IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
        CALL IOERROR('Buffer in fixed length header',A_IO,LEN_IO,       &
     &                LEN_FIXHD)
        CMESSAGE='  FFREAD : I/O error reading FIXED LENGTH HEADER'
        ICODE=2
        RETURN
      ENDIF
      IWA= PP_FIXHD(150)-1  ! NOTE for BUFFER IN the start address
!                             ! is zero for word 1
      DATA_ADD= PP_FIXHD(160)  ! The start address of the data
      PP_LEN2_LOOKUP=PP_FIXHD(152)
! DEPENDS ON: ffreada
      CALL FFREADA(IPROJ,FCT,ITYPE,INT_LEVEL,PPUNIT,FIELD,IDIM,         &
     &ILABEL,RLABEL,IEXTRA,TEST,PP_LEN2_LOOKUP,LEN1_LOOKUP,PP_FIXHD,    &
     &IWA,LEN1_RECORD,MAXFF,RECORD,PFNO,DATA_ADD,ICODE,CMESSAGE)

! DEPENDS ON: ereport
      IF(ICODE  /=  0) CALL EREPORT("FFREAD", ICODE, CMESSAGE)
!
 9999 CONTINUE
      RETURN
      END SUBROUTINE FFREAD
!LL  Routine: READ_REC--------------------------------------------------
!LL
!LL  Purpose: To read a data record from a  pp file
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!LL  Routine: LEVEL_RLEVEL ------------------------------------------
!LL
!LL  Purpose: To return a real value even though the routine is called
!LL  with integer arguments.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: ...
!LL
!LL  Project task: ...
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!LL  Routine: UN_PACK  -------------------------------------------------
!LL
!LL  Purpose: To unpack data from the input array FIELD and return
!LL  the data in FIELD.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: ...
!LL
!LL  Project task: ...
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!LL  Routine: LOGICAL_TO_REAL ------------------------------------------
!LL
!LL  Purpose: To convert logical data within FIELD to real data.
!LL  the data in FIELD.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: ...
!LL
!LL  Project task: ...
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!LL  Routine: INTEGER_TO_REAL ------------------------------------------
!LL
!LL  Purpose: To convert logical data within FIELD to real data.
!LL  the data in FIELD.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: ...
!LL
!LL  Project task: ...
!LL
!LL  External documentation:
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
#endif
