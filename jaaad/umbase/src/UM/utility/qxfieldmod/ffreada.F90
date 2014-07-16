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
      SUBROUTINE FFREADA(IPROJ,FCT,ITYPE,INT_LEVEL,PPUNIT,FIELD,IDIM,   &
     &ILABEL,RLABEL,IEXTRA,TEST,PP_LEN2_LOOKUP,LEN1_LOOKUP,PP_FIXHD,    &
     &IWA,LEN1_RECORD,MAXFF,RECORD,PFNO,DATA_ADD,ICODE,CMESSAGE)
      IMPLICIT NONE
      EXTERNAL SETPOS,READ_REC,IOERROR,COEX

      CHARACTER CMESSAGE*(*)
      INTEGER                                                           &
     &     MAXFF                                                        &
                                  !IN  Max number of opened files
     &    ,LEN1_LOOKUP                                                  &
                                  !IN  first dimension of the lookup
     &    ,LEN1_RECORD                                                  &
                                  !INOUT First dimension of record
     &    ,PP_LEN2_LOOKUP                                               &
                                  !IN  secnd dimension of the lookup
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
                                  !IN  spare for future use
     &    ,ICODE                                                        &
                                  !OUT return code
     &    ,MAXPP                                                        &
                                  !    maximum number of unit number
     &    ,DATA_ADD                                                     &
                                  !IN  The word address of the data.
     &    ,PFNO                   !INOUT No of fields files opened
      INTEGER                                                           &
     &     PP_FIXHD(*),                                                 &
                                              !IN PPfile fixed header
     &     LOOKUP(LEN1_LOOKUP,PP_LEN2_LOOKUP) !OUTinteger lookup
      REAL                                                              &
     &     FIELD(IDIM)                                                  &
                                  !OUT array holding final output data.
     &    ,RLABEL(19)                                                   &
                                  !OUT holds real part of LOOKUP
     &    ,REAL_LEVEL             !IN  LEVEL code (could be real)
      LOGICAL                                                           &
     &     RECORD(LEN1_RECORD,MAXFF) !INOUT Record of the field no read
!     LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
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
     &    ,LEN_IO_EXPECTED                                              &
                                  ! Length od data expected in transfer
     &    ,LENGTH_OF_DATA                                               &
                                  ! Length of a particular field
     &    ,ADDR                                                         &
                                  ! Address of a field in the data store
     &    ,IN_LBVC                ! Local copy of LBVC required to searc
      real                                                              &
     &     A_IO                   ! OUTPUT from UNIT command
      LOGICAL                                                           &
     &     TEST
!
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!     READ IN LOOKUP TABLE  IF FIRST TIME THRO
!
!     Read in the LOOKUP table.
!
! DEPENDS ON: setpos
        CALL SETPOS(PPUNIT,IWA,ICODE) ! C coded routine
        LEN_IO_EXPECTED=PP_LEN2_LOOKUP*LEN1_LOOKUP
! DEPENDS ON: buffin
        CALL BUFFIN(PPUNIT,LOOKUP,LEN_IO_EXPECTED,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_IO_EXPECTED) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('Buffer in lookup table   ',A_IO,LEN_IO,         &
     &    LEN_IO_EXPECTED)
          CMESSAGE='FFREADA: I/O error reading LOOKUP TABLE  '
          ICODE=3
          RETURN
        ENDIF
      IF(IEXTRA(2) == 0) THEN  ! Allows duplicate entries to be read
        IF(PP_LEN2_LOOKUP >  LEN1_RECORD) THEN
          CMESSAGE='FFREADA: LEN1_RECORD NOT LARGE ENOUGH    '
          ICODE=4
          RETURN
        ENDIF
         DO I=1,LEN1_RECORD
         IF(RECORD(I,PFNO)) THEN
            LOOKUP(14,I)=-99
          ENDIF
        ENDDO
      ENDIF
! DEPENDS ON: ffreadb
      CALL FFREADB      (IPROJ,FCT,ITYPE,INT_LEVEL,PPUNIT,FIELD,IDIM,   &
     &ILABEL,RLABEL,IEXTRA,PP_LEN2_LOOKUP,LEN1_LOOKUP,                  &
     &IWA,LEN1_RECORD,MAXFF,RECORD,PFNO,PP_FIXHD,LOOKUP,LOOKUP,DATA_ADD,&
     &ICODE,CMESSAGE)

! DEPENDS ON: ereport
      IF(ICODE  /=  0) CALL EREPORT("FFREADA", ICODE, CMESSAGE)

      RETURN
      END SUBROUTINE FFREADA
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
