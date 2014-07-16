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
      SUBROUTINE FFREADB(IPROJ,FCT,ITYPE,INT_LEVEL,PPUNIT,FIELD,IDIM,   &
     &ILABEL,RLABEL,IEXTRA,PP_LEN2_LOOKUP,LEN1_LOOKUP,                  &
     &IWA,LEN1_RECORD,MAXFF,RECORD,PFNO,PP_FIXHD,LOOKUP,ROOKUP,DATA_ADD,&
     &ICODE,CMESSAGE)
      IMPLICIT NONE
      EXTERNAL SETPOS,READ_REC,IOERROR,COEX
#include "clookadd.h"

      CHARACTER CMESSAGE*(*)
      INTEGER                                                           &
     &     MAXFF                                                        &
                                  !IN  Max number of opened files
     &    ,LEN1_LOOKUP                                                  &
                                  !IN  first dimension of the lookup
     &    ,LEN1_RECORD                                                  &
                                  !IN  First dimension of record
     &    ,PFNO                                                         &
                                  !IN  No of fields files opened
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
     &    ,LENBUF                                                       &
                                  !OUT input buffer length for data
     &    ,NUM_CRAY_WORDS                                               &
                                  !OUT no of values in an input field
     &    ,DATA_ADD                                                     &
                                  !IN  The word address of the data.
     &    ,NVALS                  !OUT The num of points in a data field
      INTEGER                                                           &
     &     PP_FIXHD(*),                                                 &
                                              !IN  PPfile fixed header
     &     LOOKUP(LEN1_LOOKUP,PP_LEN2_LOOKUP) !OUT integer lookup
      REAL                                                              &
     &     FIELD(IDIM)                                                  &
                                  !OUT array holding final output data.
     &    ,ROOKUP(LEN1_LOOKUP,PP_LEN2_LOOKUP)                           &
                                                !OUT real lookup
     &    ,RLABEL(19)                                                   &
                                  !OUT holds real part of LOOKUP
     &    ,REAL_LEVEL             !IN  LEVEL code (could be real)
      LOGICAL                                                           &
     &     RECORD(LEN1_RECORD,MAXFF) !IN Record of the field no read
!     LOCAL VARIABLES
      REAL                                                              &
     &     AMDI                   ! Missing data indicator
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
     &    ,LENGTH_OF_DATA                                               &
                                  ! Length of a particular field
     &    ,ADDR                                                         &
                                  ! Address of a field in the data store
     &    ,IN_LBVC                                                      &
                                  ! Local copy of LBVC required to searc
     &    ,PACK_TYPE                                                    &
                                  ! Packing type N1 of LBPACK
     &    ,PACK_TYPE_I                                                  &
                                  ! Packing type N1 of LBPACK in loop
     &    ,DATA_COMP                                                    &
                                  ! Data compression N2 of LBPACK
     &    ,DATA_COMP_DEF                                                &
                                  ! Compression definition N3 of LBPACK
     &    ,NUMBER_FORMAT          ! Number representation N4 of LBPACK
!
#include "c_mdi.h"
!
!     DO I=112,112
!       CALL PR_LOOK(LOOKUP(1,1),ROOKUP(1,1),I)
!     ENDDO
!
!----------------------------------------------------------------------
!L  Search for the required FIELD
!----------------------------------------------------------------------
      IF(IEXTRA(3) == 0) THEN ! Search on LBTYP/LBLEV/LBPROJ/LBFT
        DO  I=1,PP_LEN2_LOOKUP
          IF(ITYPE == LOOKUP(LBTYP,I)) THEN
            IF(INT_LEVEL == LOOKUP(LBLEV,I)) THEN
              IF(IPROJ == LOOKUP(LBPROJ,I)) THEN
                IF(FCT == LOOKUP(LBFT,I)) THEN
                  IK=I
                  GOTO 3
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ELSE
        IN_LBVC=IEXTRA(3)
        IF(IEXTRA(4) == 0) THEN  !  IEXTRA(3) HAS LBVC so search on LBVC
          IF(INT_LEVEL == 8888.OR.INT_LEVEL == 9999) THEN ! special lev
            DO  I=1,PP_LEN2_LOOKUP
              IF(ITYPE == LOOKUP(LBTYP,I)) THEN
                IF(INT_LEVEL == LOOKUP(LBLEV,I)) THEN
                  IF(IPROJ == LOOKUP(LBPROJ,I)) THEN
                    IF(FCT == LOOKUP(LBFT,I)) THEN
                      IK=I
                      GOTO 3
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ELSE  ! Not a special level so additional search on LBVC
            DO  I=1,PP_LEN2_LOOKUP
              IF(ITYPE == LOOKUP(LBTYP,I)) THEN
                IF(INT_LEVEL == LOOKUP(LBLEV,I)) THEN
                  IF(IPROJ == LOOKUP(LBPROJ,I)) THEN
                    IF(FCT == LOOKUP(LBFT,I)) THEN
                      IF(IN_LBVC == LOOKUP(LBVC,I)) THEN
                        IK=I
                        GOTO 3
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF   ! End of special level block
!----------------------------------------------------------------------
!     Search now on BLEV ie REAL_LEVEL except for special levels      C
!     and Data on model levels (BLEV would contain BK so would need   C
!     to search in this case on LBLEV).For special level search on    C
!     just LBVC LBFT LBPROJ and LBTYP.For model level convert the     C
!     real input value to integer and search as above plus LBLEV.     C
!     Note a special level cannot have an LBVC of 9                   C
!----------------------------------------------------------------------
        ELSE IF(IEXTRA(4) == 1) THEN !  IEXTRA(4) is not zero.
! DEPENDS ON: level_rlevel
          CALL LEVEL_RLEVEL(INT_LEVEL,INT_LEVEL,REAL_LEVEL)
          IF(REAL_LEVEL <= 0.0) THEN !  Special level indicated.
            DO  I=1,PP_LEN2_LOOKUP
              IF(ITYPE == LOOKUP(LBTYP,I)) THEN
                IF(IPROJ == LOOKUP(LBPROJ,I)) THEN
                  IF(FCT == LOOKUP(LBFT,I)) THEN
                    IF(IN_LBVC == LOOKUP(LBVC,I)) THEN
                      IK=I
                      GOTO 3
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ELSE IF(IN_LBVC == 9) THEN  ! model level data
            DO  I=1,PP_LEN2_LOOKUP
              IF(ITYPE == LOOKUP(LBTYP,I)) THEN
                INT_LEVEL=REAL_LEVEL+0.0000001  !
!               IF(REAL_LEVEL <= (ROOKUP(BLEV,I)+0.0001).AND.
!    *          REAL_LEVEL >= (ROOKUP(BLEV,I)-0.0001)) THEN
! That MOD is only for un-corrected model dumps
                IF(INT_LEVEL == LOOKUP(LBLEV,I)) THEN
                  IF(IPROJ == LOOKUP(LBPROJ,I)) THEN
                    IF(FCT == LOOKUP(LBFT,I)) THEN
                      IF(IN_LBVC == LOOKUP(LBVC,I)) THEN
                        IK=I
                        GOTO 3
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ELSE        ! not model level data or a special level
            DO  I=1,PP_LEN2_LOOKUP
              IF(ITYPE == LOOKUP(LBTYP,I)) THEN
                IF(REAL_LEVEL <= (ROOKUP(BLEV,I)+0.0001).AND.           &
     &          REAL_LEVEL >= (ROOKUP(BLEV,I)-0.0001)) THEN
                  IF(IPROJ == LOOKUP(LBPROJ,I)) THEN
                    IF(FCT == LOOKUP(LBFT,I)) THEN
                      IF(IN_LBVC == LOOKUP(LBVC,I)) THEN
                        IK=I
                        GOTO 3
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
!----------------------------------------------------------------------
!     Search now on BLEV ie REAL_LEVEL except for special levels      C
!     and Data on model levels (BLEV would contain BK so would need   C
!     to search in this case on LBLEV).For special level search on    C
!     just LBVC LBFT LBPROJ and LBTYP.For model level convert the     C
!     real input value to integer and search as above plus LBLEV.     C
!     Note a special level cannot have an LBVC of 9                   C
!----------------------------------------------------------------------
        ELSE IF(IEXTRA(4) == 2) THEN !  IEXTRA(4) is not zero.
! DEPENDS ON: level_rlevel
          CALL LEVEL_RLEVEL(INT_LEVEL,INT_LEVEL,REAL_LEVEL)
          IF(REAL_LEVEL <= 0.0) THEN !  Special level indicated.
            DO  I=1,PP_LEN2_LOOKUP
              IF(ITYPE == LOOKUP(LBFC,I)) THEN
                IF(IPROJ == LOOKUP(LBPROJ,I)) THEN
                  IF(FCT == LOOKUP(LBFT,I)) THEN
                    IF(IN_LBVC == LOOKUP(LBVC,I)) THEN
                      IK=I
                      GOTO 3
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ELSE IF(IN_LBVC == 9) THEN  ! model level data
            DO  I=1,PP_LEN2_LOOKUP
              IF(ITYPE == LOOKUP(LBFC,I)) THEN
                INT_LEVEL=REAL_LEVEL+0.0000001  !
                IF(INT_LEVEL == LOOKUP(LBLEV,I)) THEN
                  IF(IPROJ == LOOKUP(LBPROJ,I)) THEN
                    IF(FCT == LOOKUP(LBFT,I)) THEN
                      IF(IN_LBVC == LOOKUP(LBVC,I)) THEN
                        IK=I
                        GOTO 3
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ELSE        ! not model level data or a special level
            DO  I=1,PP_LEN2_LOOKUP
              IF(ITYPE == LOOKUP(LBFC,I)) THEN
                IF(REAL_LEVEL == ROOKUP(BLEV,I)) THEN
                  IF(IPROJ == LOOKUP(LBPROJ,I)) THEN
                    IF(FCT == LOOKUP(LBFT,I)) THEN
                      IF(IN_LBVC == LOOKUP(LBVC,I)) THEN
                        IK=I
                        GOTO 3
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDIF    ! IEXTRA(3) == 0 IF block
      ENDIF      ! IEXTRA(4) == 0 IF block

 108  FORMAT('  FIELD FOUND FOR ITYPE,LEVEL,IPROJ,FCT',4I5)
      IF (IEXTRA(4) >= 1.AND.IEXTRA(3) /= 0) THEN
        WRITE(6,112) ITYPE,REAL_LEVEL,IPROJ,FCT
      ELSE
        WRITE(6,104) ITYPE,INT_LEVEL,IPROJ,FCT
      END IF
 104  FORMAT('  FIELD NOT FOUND FOR ITYPE,INT_LEVEL,IPROJ,FCT',4I5)
 112  FORMAT('  FIELD NOT FOUND FOR ITYPE,REAL_LEVEL,IPROJ,FCT',I5,     &
     &      F7.1,2I5)
      ICODE=1
      CMESSAGE=' FFREAD  field not found'
      GOTO 9999
    3 CONTINUE
!=== Decode LBPACK code
      PACK_TYPE = MOD(LOOKUP(LBPACK,IK),10)
      DATA_COMP = MOD(LOOKUP(LBPACK,IK),100) - PACK_TYPE
      DATA_COMP_DEF = MOD(LOOKUP(LBPACK,IK),1000) -DATA_COMP -PACK_TYPE
      NUMBER_FORMAT = LOOKUP(LBPACK,IK)/1000
!=== Reading a model type dump =======================================
!    A model dump has no direct addressing only relative.
      IF(LOOKUP(LBNREC,IK) == 0) THEN ! A model dump
        IF(PACK_TYPE == 2) THEN          ! Is the field packed.
          NUM_CRAY_WORDS=LOOKUP(LBLREC,IK)/2
        ELSE
          NUM_CRAY_WORDS=LOOKUP(LBLREC,IK)
        ENDIF
        NVALS=LOOKUP(LBLREC,IK) ! No of data points
        ADDR=DATA_ADD
        DO I=1,IK-1
          PACK_TYPE_I = MOD(LOOKUP(LBPACK,I),10)
          IF(PACK_TYPE_I == 2) THEN ! 32 Bit packed
            LENGTH_OF_DATA=LOOKUP(LBLREC,I)/2
          ELSE
            LENGTH_OF_DATA=LOOKUP(LBLREC,I)
          ENDIF
          ADDR=ADDR+LENGTH_OF_DATA
        ENDDO
        IWA=ADDR-1
      ELSE
!=== Reading a PP type file.==========================================
        NUM_CRAY_WORDS=LOOKUP(LBLREC,IK) ! PP type file
        IWA=LOOKUP(29,IK)
        NVALS=LOOKUP(44,IK)
      ENDIF
      RECORD(IK,PFNO)=.TRUE.   ! Record which the no of the field read
      LENBUF=LOOKUP(LBNREC,IK) !
!==============================================================
      IF (IEXTRA(4) >= 1.AND.IEXTRA(3) /= 0) THEN
        WRITE(7,110) ITYPE,REAL_LEVEL,IPROJ,FCT,IK,NUM_CRAY_WORDS,NVALS
      ELSE
        WRITE(7,106) ITYPE,INT_LEVEL,IPROJ,FCT,IK,NUM_CRAY_WORDS,NVALS
      END IF
  106 FORMAT(' FIELD ','ITYPE=',I3,' LEVEL=',I5,' PROJ=',I4,' FCST=',   &
     &I5,' FIELD NO',I4,' NWORDS=',I5,' NVALS=',I5)
  110 FORMAT(' FIELD FOUND','ITYPE=',I4,'LEVEL=',F7.1,'PROJ=',I4,'FCST='&
     &,I5,'FIELD NO',I4,'NWORDS=',I5,'NVALS=',I5)
        IF(IDIM <  NUM_CRAY_WORDS) THEN
          ICODE=NUM_CRAY_WORDS
          CMESSAGE='FFREAD  Idim to small ICODE holds correct value'
          GOTO 9999
        ENDIF
      ICODE=0
!     RETURN
! DEPENDS ON: read_rec
      CALL READ_REC(FIELD,NUM_CRAY_WORDS,IWA,PPUNIT,ICODE,CMESSAGE)
 2212 FORMAT('  FIELDS FILE NUMBER ',I2,'  ON UNIT',I2,2X,'BEING READ')
!L    CLOSE(PPUNIT)
        IF(ICODE == 0) THEN
          DO 5 I=1,45
          ILABEL(I)=LOOKUP(I,IK)
    5     CONTINUE
          DO 6 I=1,19
           RLABEL(I)=ROOKUP(I+45,IK)
    6     CONTINUE
        ENDIF
!=======================================================================
! At this point FIELD holds the data either PACKED or UN-PACKED
! Is the packing indicator set and is un-packing required? If so then
! the data is temp un-packed into a work ARRAY of length IDIM
        IF(PACK_TYPE >  0) THEN               ! Is the field packed.
          IF(IEXTRA(1) == 0) THEN  ! unpacking is required
!           get missing data indicator from pp header
            AMDI = ROOKUP(BMDI,IK)
!           compare with MDI from COMDECK
      IF(AMDI /= RMDI) WRITE(6,*)' WARNING non-standard MDI in use'
! DEPENDS ON: un_pack
            CALL UN_PACK(PACK_TYPE,IDIM,FIELD,NUM_CRAY_WORDS            &
     &                  ,ILABEL,AMDI,PP_FIXHD,ICODE,CMESSAGE)
          ENDIF
        ELSEIF(LOOKUP(DATA_TYPE,IK) == 3) THEN   !Fld is logical
! DEPENDS ON: logical_to_real
          CALL LOGICAL_TO_REAL(IDIM,FIELD,FIELD,NVALS,                  &
     &                         ILABEL,ICODE,CMESSAGE)
        ELSEIF(LOOKUP(DATA_TYPE,IK) == 2) THEN   !Fld is integer
! DEPENDS ON: integer_to_real
          CALL INTEGER_TO_REAL(IDIM,FIELD,FIELD,NVALS,                  &
     &                         ILABEL,ICODE,CMESSAGE)
        ENDIF
!=======================================================================
 9999 CONTINUE
  100 FORMAT(//,32X,'   ARRAY        ',//,32(16F5.0/))
  101 FORMAT(//,32X,'   LOOKUP       ',//,32(16I5/))
  103 FORMAT('   LENIN  ',I12)
! DEPENDS ON: ereport
      IF(ICODE  /=  0) CALL EREPORT("FFREADB", ICODE, CMESSAGE)
      RETURN
      END SUBROUTINE FFREADB
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
