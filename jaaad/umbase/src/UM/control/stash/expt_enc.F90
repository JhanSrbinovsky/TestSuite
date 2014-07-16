#if defined(C84_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL SUBROUTINE EXPT_ENC--------------------------------------------
!LL
!LL     Given a valid five character experiment RUN_ID, an INTEGER*4
!LL   code number is generated. The valid experiment name characters
!LL   are A-Z (uppercase), 0-9. Each letter of the experiment name is
!LL   stored as a 6 bit number. The last letter is converted to the
!LL   6 least significant bits of the integer code, the next letter
!LL   the next 6 lsb's, etc. Hence 30 bits are used in all.
!LL     The lookup table is capable of holding 64 elements. This
!LL   number cannot be exceeded if the code number is to remain
!LL   INTEGER*4. Similarly, the experiment RUN_ID length (5 chars)
!LL   cannot be exceeded.
!LL     Subroutine called from PP_HEAD.
!LL
!LL   A. Brady <- programmer of some or all of previous code or changes
!LL
!LL    Model            Modification history from model version 3.0:
!LL   version  Date
!LL
!LL   Programming standard:
!LL
!LL   Logical components covered:
!LL
!LL   Project TASK:
!LL
!LL   External documentation:
!LL
!LLEND-------------------------------------------------------------

!*L  INTERFACE and ARGUMENTS:--------------------------------------

      SUBROUTINE EXPT_ENC(EXPTSTRG                                      &
     &  ,EXPTCODE                                                       &
     &  ,ICODE                                                          &
     &  ,CMESSAGE)
!*-----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER       ICODE       !OUT  Return code: successful=0
      CHARACTER*80  CMESSAGE    !OUT  Error message if ICODE > 0

      CHARACTER*5   EXPTSTRG    !IN   Experiment name string. Length
                                !     must equal parameter STRSIZE
      INTEGER       EXPTCODE    !OUT  Experiment code integer

!     Define local variables

      LOGICAL       TEST
      CHARACTER*1   LETTER
      INTEGER       I,J,NEWNUM,LETNUM,NSTRINGS,NBITS,STRSIZE

      PARAMETER(NSTRINGS=36,                                            &
     &  NBITS=6,                                                        &
     &  STRSIZE=5)

      CHARACTER*1   USTRINGS    ! Upper case strings
      CHARACTER*1   LSTRINGS    ! Lower case strings

      DIMENSION     USTRINGS(0:NSTRINGS-1)
      DIMENSION     LSTRINGS(0:NSTRINGS-1)

      DATA USTRINGS/'A','B','C','D','E','F','G','H','I','J',            &
     &  'K','L','M','N','O','P','Q','R','S','T',                        &
     &  'U','V','W','X','Y','Z','0','1','2','3',                        &
     &  '4','5','6','7','8','9'/

      DATA LSTRINGS/'a','b','c','d','e','f','g','h','i','j',            &
     &  'k','l','m','n','o','p','q','r','s','t',                        &
     &  'u','v','w','x','y','z','0','1','2','3',                        &
     &  '4','5','6','7','8','9'/

!     Begin main

      EXPTCODE=0
      LETNUM=STRSIZE

!     Loop over letters in EXPTSTRG
      DO 20 I=0,STRSIZE-1
        TEST=.FALSE.
        READ(EXPTSTRG(LETNUM:LETNUM),"(A1)")LETTER

!       Loop over letters in lookup table USTRINGS/LSTRINGS
        DO 10 J=0,NSTRINGS-1
          IF ((LETTER == USTRINGS(J)).OR.(LETTER == LSTRINGS(J))) THEN
            TEST=.TRUE.
            NEWNUM=J*(2**(I*NBITS))
            EXPTCODE=EXPTCODE+NEWNUM
!           Exit loop as we have found the code for this letter
            GOTO 15
          ENDIF
 10     CONTINUE
 15     CONTINUE

!       Check experiment name is valid
        IF (.NOT.TEST) THEN
          ICODE=99
          CMESSAGE='EXPT_ENC: Invalid letter in expt name (RUN_ID)'
          RETURN
        ENDIF
        LETNUM=LETNUM-1
 20   CONTINUE

      RETURN
      END SUBROUTINE EXPT_ENC
#endif
