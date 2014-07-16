#if defined(HPRT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: WRITRSUB -------------------------------------------------
!LL
!LL  Purpose: To create file containing job details for qsub to resubmit
!LL           the climate model, using history block information.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Author:   A.B.SANGSTER
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL   3.1  05/02/93    Portable Fortran unit no assigns
!LL                    Author: A. Dickinson    Reviewer: R. Stratton
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL
!LL  Programming standard: UM Doc Paper 3, draft version
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!*L  Interface and arguments:
!
      SUBROUTINE WRITRSUB                                               &
     &         ( UNITRSUB,ICODE,CMESSAGE )
!
      IMPLICIT NONE
!
      INTEGER HIST_UNIT
      INTEGER      UNITRSUB  ! In  - RSUB file unit
      INTEGER       ICODE    ! Out - Return code from routine
      CHARACTER*(80) CMESSAGE ! Out - Return message if failure occured
       CHARACTER *80 FILENAME
      PARAMETER(HIST_UNIT=10)
!*
!
!L Common blocks
!
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
!
!*L EXTERNAL subroutines called
      EXTERNAL GET_FILE
!*
! 0. Read the top record of the permenant history
! file.
        CALL GET_FILE(HIST_UNIT,FILENAME,80,ICODE)
        OPEN(HIST_UNIT,FILE=FILENAME,FORM='FORMATTED',IOSTAT=ICODE)
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='HPRINT  : Failed in OPEN of history file'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
        WRITE(6,*)'HPRINT  : Warning message on OPEN of history file'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
      REWIND(HIST_UNIT)
      READ(HIST_UNIT,NLIHISTO)
      READ(HIST_UNIT,NLCHISTO)
      READ(HIST_UNIT,NLIHISTG)
      READ(HIST_UNIT,NLCHISTG)
      READ(HIST_UNIT,NLCFILES)
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='HPRINT  : Failed in READ of history file'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
        WRITE(6,*)'HPRINT  : Warning message on READ of history file'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
!L
!L 1. Open, rewind and write a record
!L
      CALL GET_FILE(UNITRSUB,FILENAME,80,ICODE)
      OPEN(UNITRSUB,FILE=FILENAME,IOSTAT=ICODE)
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='WRITRSUB: Failed in OPEN of output unit'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
        WRITE(6,*)'WRITRSUB: Warning message on OPEN of output unit'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
      REWIND(UNITRSUB)
      WRITE(UNITRSUB,100,IOSTAT=ICODE)                                  &
     &                   RUN_RESUBMIT,                                  &
     &                   RUN_RESUBMIT_Q,                                &
     &                   RUN_RESUBMIT_TIME,                             &
     &                   RUN_RESUBMIT_CPU,                              &
     &                   RUN_RESUBMIT_MEMORY,                           &
     &                   RUN_RESUBMIT_PRTY,                             &
     &                   RUN_RESUBMIT_JOBNAME,                          &
     &                   RUN_JOB_NAME
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='WRITRSUB: Failed in WRITE to output unit'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
        WRITE(6,*)'WRITRSUB: Warning message on WRITE to output unit'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
!
 100  FORMAT('FLAG   = ',A1/                                            &
     &       'QUEUE  = ',A12/                                           &
     &       'TIME   = ',A12/                                           &
     &       'CPU    = ',A12/                                           &
     &       'MEMORY = ',A12/                                           &
     &       'PRTY   = ',A12/                                           &
     &       'JOBNAME = ',A8/                                           &
     &       'THISJOB = ',A8)
 999  CONTINUE
!L
!L 2. Close and return
!L
      CLOSE(UNITRSUB)
      RETURN
      END SUBROUTINE WRITRSUB
#endif
