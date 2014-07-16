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
      SUBROUTINE CONTROL(PPUNIT1,PPUNIT2,LEN1_LOOKUP,LEN2_LOOKUP,       &
     &                   LOOKUP,PP_INTHD,LEN_INTHD,                     &
     &                   PP_FIXHD,LEN_FIXHD,ICODE,CMESSAGE,NENT)
      IMPLICIT NONE
      INTEGER                                                           &
     &     LEN_FIXHD                                                    &
     &    ,LEN_INTHD                                                    &
     &    ,LEN_LOOKUP                                                   &
     &    ,LEN1_LOOKUP                                                  &
     &    ,LEN2_LOOKUP                                                  &
     &    ,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                              &
     &    ,LOOKNEW(LEN1_LOOKUP,LEN2_LOOKUP)                             &
     &    ,PP_FIXHD(LEN_FIXHD)                                          &
     &    ,PP_INTHD(LEN_INTHD)                                          &
     &    ,LEN_IO                                                       &
     &    ,ICODE                                                        &
     &    ,PPUNIT1                                                      &
     &    ,PPUNIT2                                                      &
     &    ,NENT
      INTEGER                                                           &
     &     ROW_LENGTH                                                   &
     &    ,P_ROWS                                                       &
     &    ,P_FIELD                                                      &
     &    ,LENBUF                                                       &
     &    ,I                                                            &
     &    ,J
      REAL                                                              &
     &     A_IO

      INTEGER                                                           &
     &    STIME_MOD                                                     &
     &,   ETIME_MOD                                                     &
     &,   NFIELDS_MOD                                                   &
     &,   MTYPE_MOD(500)                                                &
     &,   MLEVS_MOD(500)                                                &
     &,   STIME_SEL                                                     &
     &,   ETIME_SEL                                                     &
     &,   NFIELDS_SEL                                                   &
     &,   MTYPE_SEL(500)                                                &
     &,   MLEVS_SEL(500)                                                &
     &,   STIME_REJ                                                     &
     &,   ETIME_REJ                                                     &
     &,   NFIELDS_REJ                                                   &
     &,   MTYPE_REJ(500)                                                &
     &,   MLEVS_REJ(500)                                                &
     &,   PPUNIT_OROG                                                   &
     &,   STIME_THI                                                     &
     &,   ETIME_THI                                                     &
     &,   NFIELDS_THI                                                   &
     &,   MTYPE_THI(500)                                                &
     &,   MLEVS_THI(500)                                                &
     &,   IXXSTEP_THI(500)                                              &
     &,   IYYSTEP_THI(500)

      CHARACTER                                                         &
     &    OUTPUT_PACK_TYPE*6

      REAL                                                              &
     &    AMULT(500)                                                    &
     &,   WIND_10M_OROG                                                 &
                                 !  LEVEL ABOVE WHICH 10M WIND FIXED
     &,   WIND_10M_SCALE         !  SCALE APPLIED TO LEVEL 1 WINDS

      LOGICAL                                                           &
     &    MODIFY                                                        &
     &,   REJECT                                                        &
     &,   SELECT                                                        &
     &,   WIND_10M                                                      &
     &,   THIN

      NAMELIST /MODS/                                                   &
     &  MODIFY,STIME_MOD,ETIME_MOD,NFIELDS_MOD,                         &
     &                                      MTYPE_MOD,MLEVS_MOD,AMULT,  &
     &  SELECT,STIME_SEL,ETIME_SEL,NFIELDS_SEL,MTYPE_SEL,MLEVS_SEL,     &
     &  REJECT,STIME_REJ,ETIME_REJ,NFIELDS_REJ,MTYPE_REJ,MLEVS_REJ,     &
     &  WIND_10M,WIND_10M_SCALE,WIND_10M_OROG,PPUNIT_OROG,              &
     &  THIN,STIME_THI,ETIME_THI,NFIELDS_THI,MTYPE_THI,MLEVS_THI,       &
     &                                        IXXSTEP_THI,IYYSTEP_THI,  &
     &  OUTPUT_PACK_TYPE

!-----------------------------------------------------------------------
      CHARACTER CMESSAGE*(*)
      EXTERNAL FIELDS
!
!L---------------------------------------------------------------
!L     init namelist
!L---------------------------------------------------------------
      MODIFY   = .FALSE.
      REJECT   = .FALSE.
      SELECT   = .FALSE.
      WIND_10M = .FALSE.
      THIN=.FALSE.
      STIME_MOD = -99
      ETIME_MOD = -99
      NFIELDS_MOD=0
      STIME_SEL = -99
      ETIME_SEL = -99
      NFIELDS_SEL=0
      STIME_REJ = -99
      ETIME_REJ = -99
      NFIELDS_REJ=0
      STIME_THI = -99
      ETIME_THI = -99
      NFIELDS_THI=0
      DO I=1,500
        MTYPE_MOD(I)=0
        MLEVS_MOD(I)=0
        AMULT(I)=1.0
        MTYPE_SEL(I)=0
        MLEVS_SEL(I)=0
        MTYPE_REJ(I)=0
        MLEVS_REJ(I)=0
        MTYPE_THI(I)=0
        MLEVS_THI(I)=0
        IXXSTEP_THI(I)=2
        IYYSTEP_THI(I)=2
      ENDDO
      WIND_10M_OROG  = -9999.
      WIND_10M_SCALE = .7
      PPUNIT_OROG    = 12
      OUTPUT_PACK_TYPE='WGDOS '
!
!L---------------------------------------------------------------
!L     read namelist
!L---------------------------------------------------------------
      READ(5,MODS)
      WRITE(7,MODS)

!L---------------------------------------------------------------
!L     Set up constants
!L---------------------------------------------------------------
      ROW_LENGTH=PP_INTHD(6)
      P_ROWS=PP_INTHD(7)
      P_FIELD=ROW_LENGTH*P_ROWS
      LENBUF=P_FIELD + 512

! DEPENDS ON: fields
      CALL FIELDS(PP_FIXHD,LEN_FIXHD,LENBUF,P_FIELD,                    &
     &             LOOKUP,LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,NENT,          &
     &             STIME_MOD,ETIME_MOD,NFIELDS_MOD,                     &
     &                                       MTYPE_MOD,MLEVS_MOD,AMULT, &
     &             STIME_SEL,ETIME_SEL,NFIELDS_SEL,MTYPE_SEL,MLEVS_SEL, &
     &             STIME_REJ,ETIME_REJ,NFIELDS_REJ,MTYPE_REJ,MLEVS_REJ, &
     &             STIME_THI,ETIME_THI,NFIELDS_THI,MTYPE_THI,MLEVS_THI, &
     &                                         IXXSTEP_THI,IYYSTEP_THI, &
     &             MODIFY,SELECT,REJECT,THIN,OUTPUT_PACK_TYPE,          &
     &             WIND_10M,WIND_10M_OROG,WIND_10M_SCALE,PPUNIT_OROG,   &
     &             PPUNIT1,PPUNIT2,ICODE,CMESSAGE)
 9999 CONTINUE
      RETURN
      END SUBROUTINE CONTROL

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
