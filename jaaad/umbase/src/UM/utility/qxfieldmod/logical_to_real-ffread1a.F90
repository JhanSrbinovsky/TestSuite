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
      SUBROUTINE LOGICAL_TO_REAL(IDIM,LOGICAL_FIELD,FIELD,NVALS,        &
     &                           ILABEL,ICODE,CMESSAGE)
      INTEGER                                                           &
     &     IDIM                                                         &
                                !IN  The full unpacked size of a field
     &    ,ILABEL(45)                                                   &
                                !OUT holds integer part of LOOKUP
     &    ,ICODE                !OUT Non zero for any error
      REAL                                                              &
     &     FIELD(IDIM)          !OUT On Input contains Real data.
      LOGICAL                                                           &
     &     LOGICAL_FIELD(IDIM)  !INOUT On Input contains logical data.
!                               ! contains the un-packed data.
      CHARACTER CMESSAGE*(*)    !OUT Will contain any error mesages.
!*
!     LOCAL  VARIABLES
      INTEGER                                                           &
     &     NVALS                  ! IN no of values in an input field
!
#include "clookadd.h"
!
      DO  I=1,NVALS
        IF(LOGICAL_FIELD(I))THEN
          FIELD(I)=1.0
        ELSE
          FIELD(I)=0.0
        ENDIF
      ENDDO
      ILABEL(DATA_TYPE)=1     ! The data type must now be real
      ICODE=0
      RETURN
      END SUBROUTINE LOGICAL_TO_REAL
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
