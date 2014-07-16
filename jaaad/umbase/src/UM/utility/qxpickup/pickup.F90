#if defined(PICK)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: PICKUP---------------------------------------------------
!LL
!LL  Purpose: To prepend the temporary history file record to the
!LL           beginning of the permanent history file
!LL
!LL           + other functions that will be added later
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.0
!LL
!LL  Author:   A.Sangster
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.1   1/02/93 : Added comdeck CHSUNITS to define NUNITS for
!LL                   extra i/o
!LL  3.5  03/05/95  Sub-models stage 1: History/control files. RTHBarnes
!LL  5.0  06/05/99  Pass RoutineName to EREPORT. M Gallani
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
!LL
!LL  Logical components covered:E2
!LL
!LL  Project task: E
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks
!LL                  Version 5  18/6/90
!LL
!
!*L Interface and arguments
!
      PROGRAM PICKUP
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
      EXTERNAL INITCHST,READHIST,WRITHIST,EREPORT,ABORT
!*
!
!  Local variables
!
      INTEGER  ICODE,IABORT   ! Work- Return codes from called routines
      INTEGER  CHIST_UNIT     ! Unit no. to which old history file is
!                             !  copied in script qspickup
      CHARACTER*256 CMESSAGE  ! Work- Return message if failure occured
      CHARACTER*(*) RoutineName
      PARAMETER (RoutineName='PICKUP')
!L
!L 1. Set common block area to zero or blank
!L
! DEPENDS ON: initchst
      CALL INITCHST
!L
!L 2. Read temporary history file records
!L
! DEPENDS ON: readhist
      CALL READHIST(THIST_UNIT,ICODE,CMESSAGE)
      IF(ICODE  >   0) GOTO 999
!L
!L 3. Prepend temporary history records to beginning of
!L     permanent history file
!L
      RUN_HIST_TYPE='Permanent'
      CHIST_UNIT = 9 ! Must check this against script qspickup
! DEPENDS ON: writhist
      CALL WRITHIST(PHIST_UNIT,CHIST_UNIT,ICODE,CMESSAGE)
      IF(ICODE  >   0) GOTO 999
 999  CONTINUE
!L
!L 4. Output error message if problem
!L
      IABORT=ICODE
! DEPENDS ON: ereport
      IF(ICODE  /=  0) CALL EREPORT(RoutineName,ICODE,CMESSAGE)
!L
!L 5. Stop and abort if problem
!L
      IF(IABORT  >   0)CALL ABORT
      STOP
      END PROGRAM PICKUP
#endif
