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

      SUBROUTINE THIN_FIELD(PDATA,RDATA,PDATA_LEN,IXX,IYY,              &
     &                      IXXSTEP,IYYSTEP,IDIM,PACK_CODE,AMDI,        &
     &                      LREC, ICODE, CMESSAGE)
!
!    Subroutine to unpack a field, thin, then repack data.
!
!    Author V Blackman      Date; 12 JAN 95
!
!
      IMPLICIT NONE
      INTEGER IDIM,PDATA_LEN,PACK_CODE,IXXSTEP,IYYSTEP
      INTEGER LREC,ICODE
      INTEGER PDATA(IDIM),IXX,IYY,ISC,LWORD
      INTEGER i,j,k,kk,ix1,iy1
      INTEGER LREC
      INTEGER ICODE
      integer countx,county
      REAL RDATA(IDIM),FIELD(IDIM),AMDI
      LOGICAL OPACK
      CHARACTER*(80) CMESSAGE
      DATA LWORD/64/

      ICODE=0
      IF(PACK_CODE == 1) THEN
        OPACK=.FALSE.
! DEPENDS ON: coex
        CALL COEX(FIELD,IDIM,PDATA,IDIM,IXX,IYY,PDATA_LEN,              &
     &            ISC,OPACK,AMDI,LWORD,ICODE,CMESSAGE)

! If IXX and IYY are not decreased by 1 then GRDSET ( a PP routine)
! will fail and give the message 'BAD GRID DEFINITION'.
! Unfortunately the same failure occurs if IXX and IYY are decreased
! when a step size of 1 is specified so IXX and IYY will only be
! decreased for step sizes > 1 (in case anyone uses a step size of 1
! instead of using SELECT in the namelist)

        if(ixxstep >  1) then
          IX1 = IXX - 1
        else
          IX1 = IXX
        endif
        if(iyystep >  1) then
          IY1 = IYY - 1
        else
          IY1 = IYY
        endif

        county = 0

        K = 1
        DO J=1,IY1,IYYSTEP
          countx = 0
          DO I=1,IX1,IXXSTEP
            kk = (j-1) * ixx + i
            FIELD(K) = FIELD(KK)
            K = K + 1
            countx = countx + 1
          END DO
          county = county + 1
        END DO

        IXX = (IX1 + IXXSTEP - 1) / IXXSTEP
        IYY = (IY1 + IYYSTEP - 1) / IYYSTEP

        OPACK=.TRUE.
! DEPENDS ON: coex
        CALL COEX(FIELD,IDIM,PDATA,IDIM,IXX,IYY,PDATA_LEN,              &
     &            ISC,OPACK,AMDI,LWORD,ICODE,CMESSAGE)
        IF(ICODE /= 0) THEN
          WRITE(7,*) ICODE
          WRITE(7,*) 'THIN_FIELD FAILED WITH THE FOLLOWING REASON:'
          WRITE(7,*) CMESSAGE
        ENDIF

      ELSE IF(PACK_CODE == 4) THEN

! DEPENDS ON: runlen_decode
        CALL RUNLEN_DECODE(FIELD,IXX*IYY,PDATA,LREC,                    &
     &                     AMDI,ICODE,CMESSAGE )

        if(ixxstep >  1) then
          IX1 = IXX - 1
        else
          IX1 = IXX
        endif
        if(iyystep >  1) then
          IY1 = IYY - 1
        else
          IY1 = IYY
        endif

        county = 0

        K = 1
        DO J=1,IY1,IYYSTEP
          countx = 0
          DO I=1,IX1,IXXSTEP
            kk = (j-1) * ixx + i
            FIELD(K) = FIELD(KK)
            K = K + 1
            countx = countx + 1
          END DO
          county = county + 1
        END DO

        IXX = (IX1 + IXXSTEP - 1) / IXXSTEP
        IYY = (IY1 + IYYSTEP - 1) / IYYSTEP

! DEPENDS ON: runlen_encode
        CALL RUNLEN_ENCODE(FIELD,IXX*IYY,PDATA,IXX*IYY,                 &
     &                     PDATA_LEN,AMDI,ICODE,CMESSAGE)
      ELSE IF(PACK_CODE == 0) THEN

        if(ixxstep >  1) then
          IX1 = IXX - 1
        else
          IX1 = IXX
        endif
        if(iyystep >  1) then
          IY1 = IYY - 1
        else
          IY1 = IYY
        endif

        county = 0

        K = 1
        DO J=1,IY1,IYYSTEP
          countx = 0
          DO I=1,IX1,IXXSTEP
            kk = (j-1) * ixx + i
            RDATA(K) = RDATA(kk)
            K = K + 1
            countx = countx + 1
          END DO
          county = county + 1
        END DO

        IXX = (IX1 + IXXSTEP - 1) / IXXSTEP
        IYY = (IY1 + IYYSTEP - 1) / IYYSTEP
        PDATA_LEN = IXX * IYY

      ELSE
        WRITE(6,*)pack_code,' not yet coded'
      END IF

      RETURN
      END SUBROUTINE THIN_FIELD

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
