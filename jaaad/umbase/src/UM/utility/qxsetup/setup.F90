#if defined(SETUP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: SETUP--------------------------------------------------
!LL
!LL  Purpose: To create an Interim History File from the Master History
!LL           file namelist input and to create the associated
!LL           logical/physical filename lookup file for this run.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.0
!LL
!LL  Author:   A.Sangster
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL  3.1    29/01/93 : Added CHSUNITS to list of comdecks to define
!LL                    nunits for i/o arrays.
!LL  3.4    17/06/94   *CALL CCONTROL added: declares logical switches
!LL                   which replace *DEFs LCAL360 passed to READHK
!LL                                                  S.J.Swarbrick
!LL  3.5  30/04/95  Sub-models stage 1: History/control.  RTHBarnes.
!LL  4.5  10/11/98  Remove superfluous *CCONTROL and associated
!LL                 variables: this information is only available
!LL                 within the model and hence all such variables are
!LL                 uninitialised. R Rawlins
!LL  5.0  06/05/99  Pass RoutineName to EREPORT. M Gallani
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
!LL
!LL  Logical components covered: H81
!LL
!LL  Project task: H
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks Version 5  18/6/90
!LL
!*L Interface and arguments
!
      PROGRAM SETUP
!
      IMPLICIT NONE
!*
!
!L Common blocks
!
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
!*L EXTERNAL subroutines called
      EXTERNAL INITCHST,READMHIS,TEMPHIST,WRITFTXX,EREPORT,ABORT
!*
!
!  Local variables
!
      INTEGER IOERR
      INTEGER  ICODE,IABORT   ! Work- Return codes from called routines
      CHARACTER*80  CMESSAGE  ! Work- Return message if failure occured
      CHARACTER*(*) RoutineName
      PARAMETER (RoutineName='SETUP')
!L
!L 1. Set common block area to zero or blank
!L
! DEPENDS ON: initchst
      CALL INITCHST
!L
!L 2. Read Master History file namelist information
!L
! DEPENDS ON: readmhis
      CALL READMHIS(MCTL_UNIT,ICODE,CMESSAGE)
      IF(ICODE  >   0) GOTO 999
!
!L
!L 3. Read Operational housekeeping file
!L  Transfer to model
!
!L
!L 4. Write history common block data to Interim History File
!L
      RUN_HIST_TYPE='Interim'
! DEPENDS ON: temphist
      CALL TEMPHIST(IHIST_UNIT,ICODE,CMESSAGE)
      IF(ICODE  >   0) GOTO 999
!L
!L 5. Write logical/physical filename file
!L
! DEPENDS ON: writftxx
      CALL WRITFTXX(FTXX_UNIT,ICODE,CMESSAGE)
      IF(ICODE  >   0) GOTO 999
 999  CONTINUE
!L
!L 6. Output error message if problem
!L
      IABORT=ICODE
! DEPENDS ON: ereport
      IF(ICODE  /=  0) CALL EREPORT(RoutineName,ICODE,CMESSAGE)
!L
!L 7. Stop and abort if error occured
!L
      IF(IABORT  >   0)CALL ABORT
      STOP
      END PROGRAM SETUP
#endif
