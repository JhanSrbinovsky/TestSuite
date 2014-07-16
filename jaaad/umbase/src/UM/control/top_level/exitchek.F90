#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: EXITCHEK -------------------------------------------------
!LL
!LL  Purpose: Checks for end-of-run condition and returns a logical to
!LL           the top level.  There are three reasons for stopping,
!LL           namely:
!LL           (i)   Model has completed the required integration;
!LL           (ii)  Operator has requested model to stop;
!LL           (iii) CPU time remaining is insufficient to complete a
!LL                 further batch of timesteps.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Author:   T.C.Johns
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
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
      SUBROUTINE EXITCHEK                                               &
     &         ( internal_model, LEXITNOW )
!
      IMPLICIT NONE
      INTEGER internal_model ! In  - id of current internal model
      LOGICAL      LEXITNOW  ! Out - True/False flag for stopping
!
!*----------------------------------------------------------------------
!  Common blocks
!
#include "chsunits.h"
#include "cmaxsize.h"
#include "csubmodl.h"
#include "chistory.h"
#include "ccontrol.h"
#include "ctime.h"
#if defined(MPP)
#include "parvars.h"
#endif
!
! External Subroutines
!
      EXTERNAL GET_FILE
!
!  Local variables
!
      LOGICAL   LERRFLAG                                                &
                         ! error flag true/false indicates stop request
     &         ,LCHK14   ! true if to check unit 14
      CHARACTER*80 CERRMESS  ! explaination of stop request.
      CHARACTER*80 FILENAME
      INTEGER   ICODE        ! IO status
#if defined(MPP)
      INTEGER end_run   ! integer used in bcast of end_condition state
      INTEGER info
#endif
!
      DATA LCHK14/.TRUE./
!
!L----------------------------------------------------------------------
!L 1. Check for completed run (unless means remain to be completed)
!L
      IF (.NOT.RUN_MEANS_TO_DO == "Y") THEN

        LEXITNOW = ( STEPim(a_im) >= TARGET_END_STEPim(a_im) )

      ENDIF
!L----------------------------------------------------------------------
!L 2. Check for server stop request.
!L
      IF (MODEL_STATUS /= "Operational".AND.LCHK14) THEN
!     Ignore server failure until have reached end of ocean part of
!     dump period so correct restart dumps are available in coupled
!     models where the dump frequency is different to the coupling
!     period.  (May not be needed once a permanent fix is done.)
      IF (.NOT. (N_SUBMODEL_PARTITION  >  1 .AND.                       &
     &    SUBMODEL_FOR_SM(INTERNAL_MODEL)                               &
     &     /=  SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION))) then
#if defined(MPP)
      IF (mype  ==  0) THEN
#endif
        CALL GET_FILE(14,FILENAME,80,ICODE)
        OPEN(14,FILE=FILENAME,IOSTAT=ICODE)
        IF(ICODE  /=  0) THEN
          WRITE(6,*)'EXITCHK: Error trying to read unit 14 error flag'
          LCHK14=.FALSE.
          RETURN
        END IF
        READ(14,10,IOSTAT=ICODE) LERRFLAG,CERRMESS
        IF(ICODE  /=  0) THEN
          WRITE(6,*)'EXITCHK: Error trying to read unit 14 error flag'
          LCHK14=.FALSE.
          RETURN
        END IF
   10   FORMAT(L1,A80)
        CLOSE(14)
#if defined(MPP)
      IF (LERRFLAG) THEN
        end_run=1
      ELSE
        end_run=0
      ENDIF
      ENDIF

      CALL GC_IBCAST(1,1,0,nproc,info,end_run)

      IF (end_run  ==  1) THEN
        LERRFLAG=.TRUE.
        IF (mype  /=  0) THEN
          CERRMESS='PE 0 Signaled a sever stop request'
        ENDIF
      ELSE
        LERRFLAG=.FALSE.
      ENDIF
#endif
!
        IF (LERRFLAG) THEN
          WRITE(6,*)'EXITCHK: Request to stop model run received'
          WRITE(6,*)CERRMESS
          LEXITNOW=.TRUE.
        END IF
      END IF
      END IF
!L----------------------------------------------------------------------
!L 3. Check for insufficient time to complete a batch of timesteps
!L      .. not yet implemented
!L
      RETURN
      END SUBROUTINE EXITCHEK
#endif
