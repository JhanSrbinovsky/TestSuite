#if defined(C70_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: EXITPROC -------------------------------------------------
!LL
!LL  Purpose: Tidies up at the end of the run, and takes certain
!LL           actions in the case of model failure (as indicated by
!LL           ICODE input)
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL 3.4  16/6/94 : Change CHARACTER*(*) to CHARACTER*(80) N.Farnon
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C0 - The top-level
!LL                          control system
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE EXITPROC(ICODE,CMESSAGE)
      IMPLICIT NONE
      INTEGER ICODE            ! INOUT - Error code from model
      CHARACTER*(80) CMESSAGE   ! IN    - Error message from model
!
!*----------------------------------------------------------------------
!  Common blocks
!
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
!
!  Subroutines called
!
!L----------------------------------------------------------------------
!L 1. If fatal error occurred in main body of model, suppress resubmit
!L    switch to prevent model resubmit (if activated)
!L
      IF (ICODE >  0) THEN
        RUN_RESUBMIT="N"
      ENDIF
!L
!L 1.1  Reset error code
!L
      ICODE=0
 999  CONTINUE
!L----------------------------------------------------------------------
!L 2. Close named pipe unit used for communication with server
!L
      CLOSE(8)
!
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE EXITPROC
#endif
