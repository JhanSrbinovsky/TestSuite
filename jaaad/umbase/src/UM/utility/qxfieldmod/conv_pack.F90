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

      SUBROUTINE CONV_PACK(ILABEL,RLABEL,PACK_CODE,                     &
     &                     INPUT_PACK_TYPE,OUTPUT_PACK_TYPE,            &
     &                     FIELD,IDIM,LEN_FIELD,                        &
     &                     PP_FIXHD,ICODE,CMESSAGE)

      INTEGER                                                           &
     &     ILABEL(50)                                                   &
     &    ,PACK_CODE                                                    &
     &    ,IDIM                                                         &
     &    ,PP_FIXHD(*)                                                  &
     &    ,LEN_FIELD                                                    &
     &    ,ICODE
      REAL                                                              &
     &     RLABEL(19)                                                   &
     &    ,FIELD(IDIM)
      CHARACTER                                                         &
     &     INPUT_PACK_TYPE*6                                            &
     &    ,OUTPUT_PACK_TYPE*6                                           &
     &    ,CMESSAGE*(*)
      REAL                                                              &
     &     AMDI
#include "clookadd.h"

      AMDI=RLABEL(18)

! DEPENDS ON: un_pack
      CALL UN_PACK(PACK_CODE,IDIM,FIELD,LEN_FIELD,                      &
     &             ILABEL,AMDI,PP_FIXHD,ICODE,CMESSAGE)

      IF(ICODE /= 0) THEN
        WRITE(7,*) ICODE
        WRITE(7,*) CMESSAGE
      ENDIF

      LEN_FIELD = ILABEL(LBROW) * ILABEL(LBNPT)
      WRITE(6,*) INPUT_PACK_TYPE,' NOW UNPACKED'

      IF(OUTPUT_PACK_TYPE == 'NONE  ') THEN
        pack_code=0   ! no repacking needed
      ELSEIF(OUTPUT_PACK_TYPE == 'WGDOS ') THEN
        pack_code=1   ! repack using coex
! DEPENDS ON: re_pack
        CALL RE_PACK(PACK_CODE,IDIM,FIELD,LEN_FIELD,                    &
     &             ILABEL,RLABEL,PP_FIXHD,ICODE,CMESSAGE)
      ELSEIF(OUTPUT_PACK_TYPE == 'CRAY32') THEN
        pack_code=3   ! repack using cray 32
        WRITE(6,*) 'packing not supported'

      ELSEIF(OUTPUT_PACK_TYPE == 'GRIB  ') THEN
        pack_code=3   ! repack using grib
! DEPENDS ON: re_pack
        CALL RE_PACK(PACK_CODE,IDIM,FIELD,LEN_FIELD,                    &
     &             ILABEL,RLABEL,PP_FIXHD,ICODE,CMESSAGE)
      ELSEIF(OUTPUT_PACK_TYPE == 'RUNLEN') THEN
        pack_code=4   ! repack using run length encoding
! DEPENDS ON: re_pack
        CALL RE_PACK(PACK_CODE,IDIM,FIELD,LEN_FIELD,                    &
     &             ILABEL,RLABEL,PP_FIXHD,ICODE,CMESSAGE)
      ENDIF

      IF(ICODE /= 0) THEN
        WRITE(7,*) ICODE
        WRITE(7,*) CMESSAGE
      ENDIF
      ILABEL(LBLREC) = LEN_FIELD

      WRITE(6,*) 'NOW PACKED INTO ',OUTPUT_PACK_TYPE
!     WRITE(6,*)'len_field=',len_field

      RETURN
      END SUBROUTINE CONV_PACK
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
#endif
