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
      SUBROUTINE UN_PACK(PACK_TYPE,IDIM,FIELD,NUM_CRAY_WORDS,           &
     &                   ILABEL,AMDI,PP_FIXHD,ICODE,CMESSAGE)
      INTEGER                                                           &
     &     PACK_TYPE                                                    &
                                !IN  The type of packing used
     &    ,IDIM                                                         &
                                !IN  The full unpacked size of a field
     &    ,ILABEL(45)                                                   &
                                !OUT holds integer part of LOOKUP
     &    ,ICODE                                                        &
                                !OUT Non zero for any error
     &    ,PP_FIXHD(*)          !IN  PPfile fixed length header
      REAL                                                              &
     &     FIELD(IDIM)                                                  &
                                !INOUT On Input contains data.On output
     &    ,AMDI                 !IN  Missing data indicator.
!                               ! contains the un-packed data.
      CHARACTER CMESSAGE*(*)    !OUT Will contain any error mesages.
!
!     LOCAL  VARIABLES
      REAL                                                              &
     &     WORK_ARRAY(IDIM)       !WORK array used for un_packing
      INTEGER                                                           &
     &     LEN_FULL_WORD                                                &
                                  ! The length of a FULL_WORD
     &    ,IXX                                                          &
                                  ! Returned X dimension from COEX
     &    ,IYY                                                          &
                                  ! Returned Y dimension from COEX
     &    ,IDUM                                                         &
                                  ! Dummy variable
     &    ,NUM_CRAY_WORDS                                               &
                                  ! IN no of values in an input field
     &    ,NUM_UNPACK_VALUES      ! Number of numbers originally packed
!
#include "clookadd.h"
!
      DATA LEN_FULL_WORD/64/
!
      IF(PACK_TYPE == 1) THEN     ! WGDOS packing
! DEPENDS ON: coex
        CALL COEX(WORK_ARRAY,IDIM,FIELD,NUM_CRAY_WORDS,IXX,IYY,         &
     &  IDUM,IDUM,.FALSE.,AMDI,LEN_FULL_WORD,ICODE,CMESSAGE)
        NUM_UNPACK_VALUES=IXX*IYY
      ELSEIF(PACK_TYPE == 3) THEN !  GRIB packing
#if defined(NECSX6)
! DEPENDS ON: degrib
        CALL DEGRIB(FIELD,WORK_ARRAY,IDIM,NUM_CRAY_WORDS,               &
     &              ILABEL,AMDI,NUM_UNPACK_VALUES,LEN_FULL_WORD)
#else
        WRITE(6,*) 'Grib unpacking only supported on NEC SX6.'
! DEPENDS ON: ereport
        CALL EREPORT('UN_PACK', 1000,                                   &
     &   'GRIB unpacking only support on NEC SX6/8')

#endif
      ELSEIF(PACK_TYPE == 4) THEN ! Run length encoded data
        NUM_UNPACK_VALUES = ILABEL(LBNPT) * ILABEL(LBROW)
! DEPENDS ON: runlen_decode
        CALL RUNLEN_DECODE(WORK_ARRAY,IDIM,FIELD,NUM_CRAY_WORDS,        &
     &                     AMDI,ICODE,CMESSAGE)
      ELSE
        ICODE=6
        CMESSAGE=' UNPACK - packing type not yet supported'
      ENDIF
      DO 8 I=1,NUM_UNPACK_VALUES
      FIELD(I)=WORK_ARRAY(I)
   8  CONTINUE
      ILABEL(DATA_TYPE)=1  ! The data type must now be real
      ILABEL(LBPACK)=ILABEL(LBPACK)-PACK_TYPE ! data no longer packed
! DEPENDS ON: ereport
      IF(ICODE  /=  0) CALL EREPORT("UN_PACK", ICODE, CMESSAGE)
      RETURN
      END SUBROUTINE UN_PACK
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
