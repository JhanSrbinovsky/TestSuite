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
      SUBROUTINE FLD_OUT                                                &
     &          (ICODE,CMESSAGE,BUFOUT,LENBUF,LEN_BUF_WORDS,NUM_WORDS,  &
     &           UNITPP,LEN1_LOOKUP,PP_LEN2_LOOKUP,IPPLOOK,RPPLOOK,     &
     &           ILABEL,RLABEL,IWL,DATA_ADDR)
      IMPLICIT NONE

      CHARACTER*(*) CMESSAGE !OUT OUT MESSAGE FROM ROUTINE
!

      INTEGER                                                           &
     &  ICODE                                                           &
                           !IN    RETURN CODE FROM ROUTINE
     &, LEN1_LOOKUP                                                     &
                           !IN    FIRST DIMENSION OF LOOKUP TABLE
     &, PP_LEN2_LOOKUP                                                  &
                           !IN    SECND DIMENSION OF LOOKUP TABLE
     &, LENBUF                                                          &
                           !IN     LENGTH OFF PP BUFFER
     &, UNITPP                                                          &
                           !IN     OUTPUT PP UNIT NUMBER
     &, LEN_BUF_WORDS                                                   &
                           !IN
     &, NUM_WORDS          !IN
!
      INTEGER                                                           &
     &  JJ            !IN    ITEM NUMBER
      INTEGER                                                           &
     &  IPPLOOK(LEN1_LOOKUP,PP_LEN2_LOOKUP)                             &
                                            !IN INTEGER LOOKUP TABLE
     &, ILABEL(45)                                                      &
                      ! INTEGER PART OF LOOKUP
     &, IWL                                                             &
                      !IN    Address of the PP LOOKUP Table
     &, DATA_ADDR     !IN    Address of start of data
!
      REAL                                                              &
     &  BUFOUT(LENBUF)                                                  &
                               !OUTPUT PP BUFFER (ROUNDED UP)
     &, RPPLOOK(LEN1_LOOKUP,PP_LEN2_LOOKUP)                             &
                                            !IN REAL LOOKUP TABLE
     &, RLABEL(19)    ! REAL PART OF LOOKUP

!*---------------------------------------------------------------------

!*L  WORKSPACE USAGE:-------------------------------------------------
!   DEFINE LOCAL WORKSPACE ARRAYS: 1 REAL ARRAY
!   AT FULL FIELD LENGTH
!
!*---------------------------------------------------------------------
!     EQUIVALENCE(IPPLOOK,RPPLOOK)
!
!*L EXTERNAL SUBROUTINES CALLED---------------------------------------
      EXTERNAL SETPOS
!*------------------------------------------------------------------
!L  MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS-1) * ROWLENGTH
!L---------------------------------------------------------------------
!----------------------------------------------------------------------
!    DEFINE LOCAL VARIABLES
      INTEGER                                                           &
     &  ADDR                                                            &
                      !
     &, IWA                                                             &
                      !     RECORD NUMBER
     &, IX                                                              &
                      !     RETURN VALUE FROM UNIT COMMAND
     &, LEN_IO                                                          &
                      !
     &, II                                                              &
                      !     COUNTER
     &, I                                                               &
                      !     COUNTER
     &,IERR           ! Error return from SETPOS

      real                                                              &
     &  A_IO          !

      INTEGER                                                           &
     &  LRESID                                                          &
                      !
     &, ICURRLL                                                         &
                      !
     &, IPAST                                                           &
                      !
     &, IPROJ         !     M08 PROJECTION NUMBER

      LOGICAL                                                           &
     &  FIRST              !
      DATA FIRST/.TRUE./
!
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS IPPLOOK GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!
      FIRST=.TRUE.
      ICURRLL=0

  103 FORMAT(//,32X,' ARRAY FROM START OF PPOUT  ',//,32(10F8.0/))
      LRESID=LEN_BUF_WORDS-NUM_WORDS
      DO 2 JJ=NUM_WORDS+1,LRESID
      BUFOUT(JJ)= 0.0
    2 CONTINUE
!
      IF(FIRST) THEN
        DO 3 JJ=1,PP_LEN2_LOOKUP
           IF(IPPLOOK(1,JJ) <  0) THEN  ! Search for last entry
             ICURRLL=JJ
               IF(JJ == 1) THEN
                 IWA=((IWL+511)/512)*512+PP_LEN2_LOOKUP*LEN1_LOOKUP
                 write(6,*) 'Start data',iwa,data_addr,iwa-1
                 IWA=DATA_ADDR
                 IWA=IWA-1
               ELSE
                IWA= IPPLOOK(29,JJ-1)+IPPLOOK(30,JJ-1) !ADDR+LGTH
               ENDIF
             GOTO 4
           ENDIF
    3   CONTINUE
          ICODE=1
          CMESSAGE="FROM PPOUT CANNOT FIND SUITABLE ENTRY IN LOOKUP"
          GOTO 999
    4     CONTINUE
      ELSE
          IPAST=ICURRLL-1
        WRITE(7,105) IPAST
 105    FORMAT('  FROM PPOUT AND FIRST IS FALSE IPAST=',I8)
          IWA=IPPLOOK(29,IPAST) + IPPLOOK(30,IPAST) ! ADDR + LENGTH
        WRITE(7,106) IWA
 106    FORMAT('  FROM PPOUT AND FIRST IS FALSE IWA=',I8)
      ENDIF
!
!     update lookup for this field
      DO I=1,45
        IPPLOOK(I,ICURRLL) = ILABEL(I)
      ENDDO
      DO I=1,19
        RPPLOOK(I+45,ICURRLL) = RLABEL(I)
      ENDDO
      IPPLOOK(29,ICURRLL)=IWA
      IPPLOOK(30,ICURRLL)=LEN_BUF_WORDS
      IPPLOOK(40,ICURRLL)=IWA
!

! DEPENDS ON: setpos
      CALL SETPOS(UNITPP,IWA,IERR)
! DEPENDS ON: buffout
      CALL BUFFOUT(unitpp,bufout,LEN_buf_words,LEN_IO,A_IO)
  100 FORMAT(//,32X,'   ARRAY BUFOUT AT END OF PPOUT ',//,32(10F8.0/))
  101 FORMAT(//,32X,'   IPPLOOK AT END OF  PPOUT   ',//,32(16I5/))
  102 FORMAT('     IWA  LEN_BUF_WORDS ',2I12)
  999 CONTINUE
      RETURN
      END SUBROUTINE FLD_OUT
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
