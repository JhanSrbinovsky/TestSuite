#if defined(SETUP) || defined(COMB)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: WRITFTXX -------------------------------------------------
!LL
!LL  Purpose: To create file of user fortran unit assignment details
!LL           taken from common history block information.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.0
!LL
!LL  Author:   A.B.SANGSTER
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL  3.1    2/02/93 : Added comdeck CHSUNITS to define NUNITS for
!LL                   extra I/O
!LL
!LL   3.1  05/02/93    Portable Fortran unit no assigns
!LL                    Author: A. Dickinson    Reviewer: R. Stratton
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!*L  Interface and arguments:
!
      SUBROUTINE WRITFTXX                                               &
     &         ( UNITFTXX,ICODE,CMESSAGE )
!
      IMPLICIT NONE
!
      INTEGER      UNITFTXX  ! In  - FTXX file unit
      INTEGER       ICODE    ! Out - Return code from routine
      CHARACTER*(80) CMESSAGE ! Out - Return message if failure occured
       CHARACTER *80 FILENAME
!*
!
!L Common blocks
!
#include "chsunits.h"
#include "clfhist.h"
!
!*L EXTERNAL subroutines called
      EXTERNAL GET_FILE
!*
!L  local variables
!
      INTEGER        I  ! index for loop
!
!L
!L 1. Open, rewind and write a record
!L
      CALL GET_FILE(UNITFTXX,FILENAME,80,ICODE)
      OPEN(UNITFTXX,FILE=FILENAME,IOSTAT=ICODE)   ! OPEN THE OUTPUT FILE
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='WRITFTXX: Failed in OPEN of output unit'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
        WRITE(6,*)'WRITFTXX: Warning message on OPEN of output unit'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
      REWIND(UNITFTXX)
      DO I=1,NUNITS
        WRITE(UNITFTXX,100,IOSTAT=ICODE)MODEL_FT_UNIT(I)
      ENDDO
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='WRITFTXX: Failed in WRITE to output unit'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
        WRITE(6,*)'WRITFTXX: Warning message on WRITE to output unit'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
!
 100  FORMAT(A80)
 999  CONTINUE
!L
!L 2. Close and return
!L
      CLOSE(UNITFTXX)
      RETURN
      END SUBROUTINE WRITFTXX
#endif
