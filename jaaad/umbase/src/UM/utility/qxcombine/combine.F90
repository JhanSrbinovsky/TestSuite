#if defined(COMB)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: COMBINE--------------------------------------------------
!LL
!LL  Purpose: To create an Interim History File by combining information
!LL           in the interim control file with that in the existing
!LL           permanent history file .If the run is operational, data
!LL           is also incorporated from the Operational Houskeeping File
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.0
!LL
!LL  Author:   A.Sangster
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.1  29/01/93 : added CHSUNITS to comdecks - defines nunits for i/o
!LL  3.4  17/06/94   *CALL CCONTROL added - declares logical switches
!LL                   which replace *DEFs             S.J.Swarbrick
!LL                  Argument LCAL360 passed to READHK
!LL  3.5  03/05/95  Sub-models stage 1: History/control files. RTHBarnes
!LL  4.5  10/11/98  Remove superfluous *CCONTROL and associated
!LL                 variables: this information is only available
!LL                 within the model and hence all such variables are
!LL                 uninitialised. R Rawlins
!LL  5.0  06/05/99  Pass RoutineName to EREPORT. M Gallani
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
!LL
!LL  Logical components covered: H82
!LL
!LL  Project task: H
!LL
!LL  Documentation:  Unified Model Documentation Paper
!LL                  H- History Bricks  Version 5  18/6/90
!LLEND --------------------------------------------------------------
!
!*L Interface and arguments
!
      PROGRAM COMBINE
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
      EXTERNAL INITCHST,READHIST,TEMPHIST,WRITFTXX,                     &
     &         EREPORT,ABORT
!*
!
!  Local variables
!
      INTEGER  ICODE,IABORT   ! Work- Return codes from called routines
      CHARACTER*256 CMESSAGE  ! Work- Return message if failure occured
      CHARACTER*(*) RoutineName
      INTEGER HIST_UNIT       ! Work- Input unit for history file
      PARAMETER(HIST_UNIT=10)
      PARAMETER (RoutineName='COMBINE')
!L
!L 1. Set common block area to zero or blank
!L
! DEPENDS ON: initchst
      CALL INITCHST
!L
!L 2. Read History file into common history block
!L
! DEPENDS ON: readhist
      CALL READHIST(HIST_UNIT,ICODE,CMESSAGE)
      IF(ICODE  >   0) GOTO 999
!L
!L 3. Read Interim Control file namelist information
!L
!!    CALL READINTC(ICTL_UNIT,ICODE,CMESSAGE)
!!    IF(ICODE  >   0) GOTO 999
!L
!L 4. Switch file assignments so old restart dump is new start dump
!L    ( Continuation runs only )
!L
      IF(H_STEPim(a_im)  >   0) ASTART(11:80) = ARESTART(11:80)
      IF(H_STEPim(o_im)  >   0) OSTART(11:80) = ORESTART(11:80)
!
!L
!L 5. Read Operational housekeeping file
!L
!  [Reading housekeeping moved from small execs into model at vn3.5.]
!
!L
!L
!L 6. Write history common block data to Interim History File
!L
      RUN_HIST_TYPE='Interim'
! DEPENDS ON: temphist
      CALL TEMPHIST(IHIST_UNIT,ICODE,CMESSAGE)
      IF(ICODE  >   0) GOTO 999
!L
!L 7. Re-write file of user assigned unit details
!L
! DEPENDS ON: writftxx
      CALL WRITFTXX(FTXX_UNIT,ICODE,CMESSAGE)
      IF(ICODE  >   0) GOTO 999
 999  CONTINUE
!L
!L 8. Output error message if problem
!L
      IABORT=ICODE
! DEPENDS ON: ereport
      IF(ICODE  /=  0) CALL EREPORT(RoutineName,ICODE,CMESSAGE)
!L
!L 9. Stop and abort if error has occurred
!L
      IF(IABORT  >   0)CALL ABORT
      STOP
      END PROGRAM COMBINE
#endif
