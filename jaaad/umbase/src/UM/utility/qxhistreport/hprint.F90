#if defined(HPRT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: HPRINT-------------------------------------------------
!LL
!LL  Purpose: Master routine for printing out summary reports from
!LL           History File records.
!LL
!LL           Also updates model resubmit job parameters from latest
!LL           history block information
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.0
!LL
!LL  Author:   A.Sangster
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL
!LL   3.1  05/02/93    Portable Fortran unit no assigns
!LL                    Author: A. Dickinson    Reviewer: R. Stratton
!LL  5.0  06/05/99  Pass RoutineName to EREPORT. M Gallani
!    6.2  10/04/05 Removed calls to ABORT.  J. Gill
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
!LL
!LL  Logical components covered: H5
!LL
!LL  Project task: H
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LL                  Version 5  18/6/90
!LL
!*L Interface and arguments
!
      PROGRAM HPRINT
!
      IMPLICIT NONE
!*
!
! Common blocks
!
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
!
!*L EXTERNAL subroutines called
      EXTERNAL INITCHST,PRINTHST,READHIST,WRITRSUB,EREPORT,ABORT
      EXTERNAL GET_FILE
!*
!
!  Local variables
!
      INTEGER  ICODE,                                                   &
                              !)Work- Return codes from called routines
     &         IABORT,                                                  &
                              !)
     &         ICOUNT,                                                  &
                              ! Work- History record counter
     &         HIST_UNIT,                                               &
                              ! Work- History file unit number
     &         RSUB_UNIT                                                &
                              ! Work- Resubmit parameter file unit no.
     &     ,NO_OF_RECORDS  !No of records user wants printed
!
      CHARACTER*256 CMESSAGE  ! Work- Return message if failure occured
      CHARACTER*80 FILENAME
      CHARACTER*(*) RoutineName
      PARAMETER (RoutineName='HPRINT')
!
      LOGICAL                                                           &
     &LONG,                                                             &
                    ! If true,  print out expanded history report
     &LAST_RECORD   ! If true,  process last history record only
!                   ! If false, process all history records
!
      PARAMETER(HIST_UNIT=10)
      PARAMETER(RSUB_UNIT=7)
!
      NAMELIST/PRINTOPT/                                                &
     &LONG,LAST_RECORD,NO_OF_RECORDS
!L
!L 0. Set default values and read namelist input
!L
      ICOUNT=0
      ICODE=0
      LONG          = .TRUE.
      LAST_RECORD   = .TRUE.
      NO_OF_RECORDS = 100
      CMESSAGE='HPRINT  : Problem reading namelist PRINTOPT'
      READ(5,PRINTOPT,END=50,ERR=999)
  50  CONTINUE
!L
!L 1. Set common block area to zero or blank
!L
! DEPENDS ON: initchst
      CALL INITCHST
!L
!L 2. Read History file and loop through records
!L
      IF(.NOT. LAST_RECORD) THEN
!
! Process each record in turn
!
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
!
 100  READ(HIST_UNIT,NLIHISTO,END=200,ERR=200)
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
!
      ICOUNT=ICOUNT+1
! DEPENDS ON: printhst
      CALL PRINTHST(ICOUNT,LONG)
      IF (ICOUNT  ==  (NO_OF_RECORDS - 1)) GOTO 200
!
      GOTO 100
 200  CONTINUE
      ELSE
!
! Process last record only
!
! DEPENDS ON: readhist
      CALL READHIST(HIST_UNIT,ICODE,CMESSAGE)
      IF(ICODE  >   0)GOTO 999
!
      LONG      = .TRUE.
!
! DEPENDS ON: printhst
      CALL PRINTHST(ICOUNT,LONG)
!
      ENDIF
!L
!L 3. Update resubmit file with current model resubmit information
!L    in history file
!L
! DEPENDS ON: writrsub
      CALL WRITRSUB(RSUB_UNIT,ICODE,CMESSAGE)
      IF(ICODE  >   0)GOTO 999
!L
!L 4.0 Normal exit
!L
      STOP
!L
!L 4.1 Output error message if problem
!L
 999  CONTINUE
      IF (ICODE == 0) ICODE=1
! DEPENDS ON: ereport
      CALL EREPORT(RoutineName,ICODE,CMESSAGE)

      STOP
      END PROGRAM HPRINT
#endif
