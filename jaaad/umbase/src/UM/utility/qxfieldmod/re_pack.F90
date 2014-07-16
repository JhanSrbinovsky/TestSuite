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
      SUBROUTINE RE_PACK(PACK_TYPE,IDIM,FIELD,NUM_CRAY_WORDS,           &
     &                   ILABEL,RLABEL,PP_FIXHD,ICODE,CMESSAGE)
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
     &    ,RLABEL(19)           !    holds real part of LOOKUP
      CHARACTER CMESSAGE*(*)    !OUT Will contain any error mesages.
!
!     LOCAL  VARIABLES
      REAL                                                              &
     &     WORK_ARRAY(IDIM)                                             &
                                  ! WORK array used for packing
     &    ,AMDI                   ! Missing data indicator.
      INTEGER                                                           &
     &     LEN_FULL_WORD                                                &
                                  ! The length of a FULL_WORD
     &    ,IXX                                                          &
                                  ! X dimension for COEX
     &    ,IYY                                                          &
                                  ! Y dimension for COEX
     &    ,ISC                                                          &
                                  ! Accuracy required for COEX
     &    ,IDUM                                                         &
                                  ! Dummy variable
     &    ,NUM_CRAY_WORDS                                               &
                                  ! IN no of values in an input field
     &    ,NUM_UNPACK_VALUES                                            &
                                  ! Number of numbers originally packed
     &    ,GRIB_PACKING           ! OUT - profile for packing
!
#include "clookadd.h"
!
      DATA LEN_FULL_WORD/64/
!
      AMDI=RLABEL(18)

      IF(PACK_TYPE == 1) THEN     ! WGDOS packing
        IXX=ILABEL(LBNPT)
        IYY=ILABEL(LBROW)
        ISC=NINT(RLABEL(6))
! DEPENDS ON: coex
        CALL COEX(FIELD,IDIM,WORK_ARRAY,IDIM,IXX,IYY,                   &
     &  NUM_CRAY_WORDS,ISC,.TRUE.,AMDI,LEN_FULL_WORD,                   &
     &  ICODE,CMESSAGE)
      ELSEIF(PACK_TYPE == 2) THEN !  32 Bit CRAY packing

#if defined(CRAY)
      ELSEIF(PACK_TYPE == 3) THEN !  GRIB PACKING
        GRIB_PACKING=1
!  RLABEL is returned from FFREAD and contains LOOKUP elements 45-64.
!  PP2GRIB requires this array to contain elements 46-64 from LOOKUP.
!  As a temporary measure the call to PP2GRIB has been amended to pass
!       CALL PP2GRIB(FIELD,WORK_ARRAY,IDIM,NUM_CRAY_WORDS,GRIB_PACKING,
!    &               ILABEL,RLABEL,ICODE,CMESSAGE)
! DEPENDS ON: pp2grib
        CALL PP2GRIB(FIELD,WORK_ARRAY,IDIM,NUM_CRAY_WORDS,GRIB_PACKING, &
     &               ILABEL,RLABEL(1),ICODE,CMESSAGE)
#endif
      ELSEIF(PACK_TYPE == 4) THEN ! Run length encoding
        IXX=ILABEL(LBROW)
        IYY=ILABEL(LBNPT)
! DEPENDS ON: runlen_encode
        CALL RUNLEN_ENCODE(FIELD,IXX*IYY,WORK_ARRAY,IXX*IYY,            &
     &                     NUM_CRAY_WORDS,AMDI,ICODE,CMESSAGE)
        ! Size of run length encoded data is greater than unpacked
        ! field therefore leave field unpacked.
          if (NUM_CRAY_WORDS  >=  IXX*IYY) then
            PACK_CODE = 0
            DO I=1,IXX*IYY
             WORK_ARRAY(I) = FIELD(I)
            END DO
            NUM_CRAY_WORDS = IXX*IYY
          endif

      ELSE
        ICODE=6
        CMESSAGE=' UNPACK - packing type not yet supported'
      ENDIF
      DO 8 I=1,NUM_cray_words
      FIELD(I)=WORK_ARRAY(I)
   8  CONTINUE
      ILABEL(DATA_TYPE)=1  ! The data type must now be real
      ILABEL(LBPACK)=ILABEL(LBPACK)+PACK_TYPE ! data now packed
      RETURN
      END SUBROUTINE RE_PACK
#endif
