#if defined(HRES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: HRESET---------------------------------------------------
!LL
!LL  Purpose: To ensure only have current and previous history records
!LL           in permanent history file in operational runs.
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
!LL  3.4  22/12/94  Change H_SECT2 to length 3, like H_SECT. RTHBarnes.
!LL  3.5  16/05/95  Sub-models stage 1: namelist history file. RTHBarnes
!LL  4.0  03/11/95  Further mods for new history file structure. RTHB
!LL  5.0  06/05/99  Pass RoutineName to EREPORT. M Gallani
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
!LL
!LL  Logical components covered: F2
!LL
!LL  Project task: F
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LL
!LL
!
!*L Interface and arguments
!
      PROGRAM HRESET
!
      IMPLICIT NONE
!*
!
!L Common blocks
!
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
!
!*L EXTERNAL subroutines called
!
      EXTERNAL EREPORT,ABORT,GET_FILE
!*
!
!  Local variables
!
      INTEGER  ICODE,IABORT   ! Work- Return codes from called routines
      CHARACTER*256 CMESSAGE  ! Work- Return message if failure occured
      CHARACTER*80 FILENAME
      CHARACTER*(*) RoutineName
      PARAMETER (RoutineName='HRESET')
      INTEGER       I         ! Work- Loop counter
      INTEGER       ICOUNT    ! Work- History record counter
!L
!L 1. Set common block areas to zero or blank
!L
      DO 6 I=1,NUNITS
         MODEL_FT_UNIT(I)=' '
  6   CONTINUE
!L
!L 2. Count number of records in permanent history file
!L    and position to start of second last record.
!L    No action if only 2 or less records in file.
!L
      CALL GET_FILE(PHIST_UNIT,FILENAME,80,ICODE)
        OPEN(PHIST_UNIT,FILE=FILENAME,FORM='FORMATTED',IOSTAT=ICODE)
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='HRESET  : Failed in OPEN of permanent history file'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
      WRITE(6,*)'HRESET  : Warning message on OPEN of permanent history &
     & file'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
      REWIND(PHIST_UNIT)
      ICOUNT=0
      READ(PHIST_UNIT,NLIHISTO,END=100,ERR=200)
      READ(PHIST_UNIT,NLCHISTO,END=100,ERR=200)
      READ(PHIST_UNIT,NLIHISTG,END=100,ERR=200)
      READ(PHIST_UNIT,NLCHISTG,END=100,ERR=200)
      READ(PHIST_UNIT,NLCFILES,END=100,ERR=200)
      go to 300
!
! Check for error
!
  100 continue
        CMESSAGE='HRESET : Unexpected end in 1st READ of permanent histo&
     &ry file'
        ICODE = 100
        go to 999
  200 continue
        CMESSAGE='HRESET : Error in 1st READ of permanent history file'
        ICODE = 200
        go to 999
!
  300 continue
      ICOUNT=ICOUNT+1
!L
!L 3. Read next records
!L
      READ(PHIST_UNIT,NLIHISTO,END=400,ERR=500)
      READ(PHIST_UNIT,NLCHISTO,END=400,ERR=500)
      READ(PHIST_UNIT,NLIHISTG,END=400,ERR=500)
      READ(PHIST_UNIT,NLCHISTG,END=400,ERR=500)
      READ(PHIST_UNIT,NLCFILES,END=400,ERR=500)
      go to 600
!
! Check for error
!
  400 continue
      WRITE(6,*)' HRESET : End of file - only one set of history records&
     &present'
      icode = 0
        go to 999
  500 continue
        CMESSAGE='HRESET : Error in 2nd READ of permanent history file'
        ICODE = 500
        go to 999
!
  600 continue
!
!L
!L 4. Use ENDFILE to truncate phist file to 2 sets of records
!L
        ENDFILE (PHIST_UNIT,IOSTAT=ICODE)
        if (icode /= 0) go to 700
        go to 999
!
! Check for error
!
  700 continue
      CMESSAGE='HRESET : Error in ENDFILE trying to reduce history file &
     &to 2 sets of records'
        ICODE = 700
 999  CONTINUE
!L
!L 5. Output error message if problem
!L
      IABORT=ICODE
! DEPENDS ON: ereport
      IF(ICODE  /=  0) CALL EREPORT(RoutineName,ICODE,CMESSAGE)
!L
!L 6. Close and stop and abort if problem
!L
      CLOSE(PHIST_UNIT)
      IF(IABORT  >   0)CALL ABORT
      STOP
      END PROGRAM HRESET
#endif
