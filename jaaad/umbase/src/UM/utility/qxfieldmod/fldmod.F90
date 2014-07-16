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
      PROGRAM FLDMOD
      IMPLICIT NONE
      EXTERNAL DIMENS1
      CHARACTER CMESSAGE*80
      CHARACTER DIAGFILE*80
      CHARACTER INFILE*80
      CHARACTER*8 c_nproc            ! to get nproc_x and nproc_y from
!                                    ! environment variables.
!                         up to an EVEN no for conversion to IBM format
      INTEGER                                                           &
     &     LEN_FIXHD                                                    &
                                  !    Length of fixed length header
     &    ,LEN_INTHD                                                    &
     &    ,LEN_REALHD                                                   &
     &    ,LEN1_LEVDPC                                                  &
     &    ,LEN2_LEVDPC                                                  &
     &    ,LEN1_ROWDPC                                                  &
     &    ,LEN2_ROWDPC                                                  &     
     &    ,LEN1_COLDPC                                                  &
     &    ,LEN2_COLDPC                                                  &
     &    ,LEN1_LOOKUP                                                  &
     &    ,LEN2_LOOKUP                                                  &
     &    ,PPUNIT1                                                      &
                                  !OUT unit no of required fieldsfile 1
     &    ,PPUNIT2                                                      &
                                  !OUT unit no of required fieldsfile 2
     &    ,DIAG_UNIT                                                    &
                                  !unit number for diagnostics
     &    ,ICODE                                                        &
                                  !IN  return code
     &    ,ERR
!    LOCAL VARIABLES
      PARAMETER(LEN_FIXHD=256)
      INTEGER                                                           &
     &     I                                                            &
                                  ! local counter
     &    ,PP_FIXHD(LEN_FIXHD)                                          &
                                  !IN  Fixed length header
     &    ,IWA                                                          &
                                  !
     &    ,IX                                                           &
                                  !
     &    ,LEN_IO                 !
      REAL                                                              &
     &     A_IO                   !
!
#include "cntl_io.h"

!L
!L Get the current sector size for disk I/O
!L

      CALL FORT_GET_ENV('UM_SECTOR_SIZE',14,c_nproc,8,err)
      IF (err  /=  0) THEN
        WRITE(6,*) 'Warning: Environment variable UM_SECTOR_SIZE has ', &
     &             'not been set.'
        WRITE(6,*) 'Setting um_sector_size to 512'
        um_sector_size=512
      ELSE
        READ(c_nproc,'(I4)') um_sector_size
      ENDIF
!    OPEN DIAGNOSTIC FILE
      DIAG_UNIT = 7
      CALL GET_FILE(DIAG_UNIT,DIAGFILE,80,ICODE)
      OPEN(UNIT=DIAG_UNIT,FILE=DIAGFILE)

!*****************************************************************
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
      icode = 0
      cmessage= '                                         '
!
!     READ IN LOOKUP TABLE  IF FIRST TIME THRO
!*****************************************************************
      PPUNIT1=10
      PPUNIT2=11
!*****************************************************************
!     Buffer in the Fixed Length Header and obtain lengths
!*****************************************************************
      CALL GET_FILE(PPUNIT1,INFILE,80,ICODE)
! DEPENDS ON: file_open
      CALL FILE_OPEN(PPUNIT1,INFILE,80,0,1,ICODE)
! DEPENDS ON: buffin
      CALL BUFFIN(PPUNIT1,PP_FIXHD,LEN_FIXHD,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('Buffer in fixed length header',A_IO,LEN_IO,     &
     &                  LEN_FIXHD)
          CMESSAGE='FFREAD : I/O error reading FIXED LENGTH HEADER'
          ICODE=2
          WRITE(6,*)'umthin1 - I/O error reading FIXED LENGTH HEADER'
          STOP
        ENDIF
      LEN_INTHD=PP_FIXHD(101)
      LEN_REALHD=PP_FIXHD(106)
      LEN1_LEVDPC=PP_FIXHD(111)
      LEN2_LEVDPC=PP_FIXHD(112)
      LEN1_ROWDPC=PP_FIXHD(116)
      LEN2_ROWDPC=PP_FIXHD(117)
      LEN1_COLDPC=PP_FIXHD(121)
      LEN2_COLDPC=PP_FIXHD(122)
      LEN1_LOOKUP=PP_FIXHD(151)
      LEN2_LOOKUP=PP_FIXHD(152)
! DEPENDS ON: dimens1
      CALL DIMENS1(LEN_INTHD,LEN_REALHD,LEN1_LEVDPC,LEN2_LEVDPC,        &
     &   LEN1_ROWDPC,LEN2_ROWDPC,LEN1_COLDPC,LEN2_COLDPC,               &
     &   LEN1_LOOKUP,LEN2_LOOKUP,LEN_FIXHD,PP_FIXHD,PPUNIT1,PPUNIT2,    &
     &   ICODE,CMESSAGE)
      IF(ICODE /= 0) THEN
        WRITE(7,100) ICODE
        WRITE(7,110) CMESSAGE
      ENDIF
      STOP
 100  FORMAT(' ICODE EQUAL TO ',I2)
 110  FORMAT(A80)
      STOP
      END PROGRAM FLDMOD
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
#endif
