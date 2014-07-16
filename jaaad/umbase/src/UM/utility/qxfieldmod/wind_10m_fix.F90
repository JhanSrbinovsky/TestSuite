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
      SUBROUTINE WIND_10M_FIX(PDATA,RDATA,PDATA_LEN,                    &
     &                        FCT,ITYPE,LEVEL,IPROJ,PPUNIT1,            &
     &                        WIND_10M_SCALE,WIND_10M_OROG,             &
     &                        MODEL_OROG,ILABEL_OROG,RLABEL_OROG,       &
     &                        IDIM,PACK_CODE,AMDI)
!LL
!LL    subroutine to unpack a 10m winds and replace if posible by
!LL    the level 1 wind scaled using wind_10m_scale
!LL
!LL    Author P J Smith      Date; 06 jan 95
!LL
!LL
      INTEGER IDIM,PDATA_LEN,PACK_CODE
      INTEGER ICODE
      CHARACTER CMESSAGE*80
      REAL RDATA(IDIM),FIELD(IDIM),FIELD1(IDIM),AMDI
      REAL MODEL_OROG(IDIM),RLABEL_OROG(19),RLABEL(19)
      REAL WIND_10M_OROG,WIND_10M_SCALE
      INTEGER PDATA(IDIM),NROW,NCOL,ISC,LWORD
      INTEGER ILABEL_OROG(45),ILABEL(45)
      INTEGER FCT,ITYPE,ITYPE1,LEVEL,LEVEL1,IPROJ,PPUNIT1
      INTEGER IEXTRA(10)
      LOGICAL OPACK
      DATA LWORD/64/

      DO I=1,10
        IEXTRA(I)=0
      ENDDO

      write(6,*) ' read level1 winds'
      ITYPE1 = 6
      IF(ITYPE == 75) ITYPE1 = 5
      LEVEL1 = 1
! DEPENDS ON: ffread
      CALL FFREAD(IPROJ,FCT,ITYPE1,LEVEL1,PPUNIT1,FIELD1,IDIM,          &
     &                ILABEL,RLABEL,IEXTRA,ICODE,CMESSAGE)

      write(6,*) 'icode=',icode
      IF(ICODE == 0) THEN
        IF(PACK_CODE == 1) THEN
          OPACK=.FALSE.
          WRITE(6,*)'call coex'
! DEPENDS ON: coex
          CALL COEX(FIELD,IDIM,PDATA,IDIM,NROW,NCOL,PDATA_LEN,          &
     &              ISC,OPACK,AMDI,LWORD,ICODE,CMESSAGE)

          WRITE(6,*)'loop field'
          DO I=1,NCOL*NROW
            IF(FIELD(I) /= AMDI) THEN
              IF(MODEL_OROG(I) >= WIND_10M_OROG) THEN
           WRITE(6,*)i,model_orog(i),field(i),field1(i),field1(i)*.8
                FIELD(I) = FIELD1(I) * WIND_10M_SCALE
              ENDIF
            ENDIF
          ENDDO

          OPACK=.TRUE.
! DEPENDS ON: coex
          CALL COEX(FIELD,IDIM,PDATA,IDIM,NROW,NCOL,PDATA_LEN,          &
     &              ISC,OPACK,AMDI,LWORD,ICODE,CMESSAGE)
        ELSEIF(PACK_CODE == 4) THEN
          OPACK=.FALSE.

! DEPENDS ON: runlen_decode
          CALL RUNLEN_DECODE(FIELD,IDIM,PDATA,IDIM,                     &
     &                       AMDI,ICODE,CMESSAGE )

          DO I=1,PDATA_LEN
            IF(FIELD(I) /= AMDI) THEN
              IF(MODEL_OROG(I) >= WIND_10M_OROG) THEN

                FIELD(I) = FIELD1(I) * WIND_10M_SCALE
              ENDIF
            ENDIF
          ENDDO
          IXX=ILABEL(18)
          IYY=ILABEL(19)
! DEPENDS ON: runlen_encode
          CALL RUNLEN_ENCODE(FIELD,IXX*IYY,PDATA,IXX*IYY,               &
     &                       PDATA_LEN,AMDI,ICODE,CMESSAGE)
        ELSEIF(PACK_CODE == 0) THEN
          DO I=1,PDATA_LEN
            IF(RDATA(I) /= AMDI) THEN
              IF(MODEL_OROG(I) >= WIND_10M_OROG) THEN
                RDATA(I) = FIELD1(I) * WIND_10M_SCALE
              ENDIF
            ENDIF
          ENDDO
        ELSE
          WRITE(6,*)pack_code,' not yet coded'
        ENDIF
      ELSE
        WRITE(7,*) ICODE
        WRITE(7,*) '10M WIND FAILED WITH THE FOLLOWING REASON:'
        WRITE(7,*) CMESSAGE
      ENDIF

      RETURN
      END SUBROUTINE WIND_10M_FIX
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
