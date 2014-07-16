#if defined(HPRT)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: PRINTHST------------------------------------------------
!LL
!LL  Purpose: To list details of history file variables
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.0
!LL
!LL  R.Stratton <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.1  05/02/93 : Altered to make format unit name printout more
!LL                   flexible.
!LL
!LL  3.1  08/03/93  A_INTF_FREQ/END/START_HR removed from output as they
!LL                  are no longer kept in history block. Claire Douglas
!LL  3.3  09/08/93  Add A_CONV_STEP freq.of conv.scheme calls. R.Barnes.
!LL
!LL  Vn 3.3  22/11/93  Portable Unified Model needs a dummy declaration
!LL                    as it does not like the FORMAT statements. (N.Far
!LL  3.4  12/07/94  New FORMAT statement 410 inserted with corresponding
!LL                 WRITE statement - writes character variables which
!LL                 correspond to logical switches in history file.
!LL                 Also alter FORMATs 425,450,500 for H_SECT to A3.
!LL                                                S.J.Swarbrick
!LL  3.5  19/07/95  New FORMAT statements for reduced History file. RTHB
!LL  4.1  18/04/96  Add RUN_ID identifier and WAVE values.  RTHBarnes.
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
      SUBROUTINE PRINTHST                                               &
     &         (ICOUNT,LONG)
!
      IMPLICIT NONE
!
      INTEGER       ICOUNT   ! In  - Record ICOUNT being processed
      LOGICAL       LONG     ! In  - Print extra history details
!*
!
!L Common block
!
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
#include "c_global.h"
!
!*L EXTERNAL subroutines called
!     None
!*
!
!  Local variables
!
      INTEGER      I                                                    &
                        !Work - Loop index variable
     &,            J                                                    &
                        !Work - Loop index variable
     &,      DUMMY      ! Dummy variable to remove error for portability
!
      DUMMY=0
!
!
  5   FORMAT('1','**  HISTORY FILE PRINTOUT  **',                       &
     & 5X,'EXPERIMENT : ',A4,2X,'JOB : ',A,5X,'RECORD ',I4/)

!
 10   FORMAT('1','**  HISTORY FILE PRINTOUT  **',                       &
     & 5X,'EXPERIMENT : ',A4,2X,'JOB : ',A,5X,'LAST RECORD'/)

!
100   FORMAT(' '/                                                       &
     &'** Type of History File read **           : ',A10/)
!
!
!
150   FORMAT('0',                                                       &
     & 20X,'******** RUN PROGRESS/RESTART DETAILS **********'//         &
     &'Run Type                                   : ',A8/               &
     &'Run indicator     (Operational runs only)  : ',I1/               &
     &'Completion code this run    (NOT SET)      : ',A14/              &
     &'Internal models currently printed are ATMOS, OCEAN, SLAB, WAVE'/ &
     &'Total timesteps completed so far  : ',I8,3(4X,I8)/               &
     &'Timesteps per coupling period     : ',I8,3(4X,I8)/               &
     &'Model Data Time- YYYY:MM:DD:HH:MM:SS',                           &
     & 5X,I4,':',4(I2,':'),I2//                                         &
     &'Last Restart dump(s) written      : ',A14,3(4X,A14)/             &
     &'Current Restart dump(s) name      : ',A14,3(4X,A14)/             &
     &'Offset between mean reference time'/                             &
     &'           and model basis time   : ',I6,3(4X,I6)/               &
     &'Ocean to be first next run ?               : ',A1/)
!
!
!
200   FORMAT('0',                                                       &
     & 20X,'******** MEAN PROCESSING DETAILS **********'//              &
     &'Number of mean periods chosen per int.model: ',4I4/              &
     &'Last mean dump created by run              : ',A14/              &
     &'Next period to be processed by MEANCTL     : ',I1//              &
     &'Partial Sum unit number indicators ' /                           &
     & 30X,'    Period_1 : ',I1,3X,I1,3X,I1,3X,I1/                      &
     & 30X,'    Period_2 : ',I1,3X,I1,3X,I1,3X,I1/                      &
     & 30X,'    Period_3 : ',I1,3X,I1,3X,I1,3X,I1/                      &
     & 30X,'    Period_4 : ',I1,3X,I1,3X,I1,3X,I1//                     &
     &'No. of mean periods in 1st year (offset)   : ',4I4/)
!CC  *'Mean to do next run?  (NOT SET)            : ',A1/)
!
!
!
300   FORMAT('0',                                                       &
     & 20X,'********    JOB RESUBMIT DETAILS **********'//              &
     &'Jobname this run                           : ',A8/               &
     &'Automatic job resubmission on ?            : ',A1/               &
     &'Resubmit Target - YYYY:MM:DD:HH:MM:SS',                          &
     & 3X,I4,':',4(I2,':'),I2//                                         &
     &'Resubmit job queue                         : ',A12/              &
     &'Resubmit job time                          : ',A20/              &
     &'Resubmit job cpu time limit                : ',A6/               &
     &'Resubmit job memory requirement            : ',A6/               &
     &'Resubmit job intra queue priority          : ',A2)
!
!
 850  FORMAT('1',                                                       &
     & 20X,'*** LOGICAL/PHYSICAL FILE ASSOCIATIONS ****' /)
 851  FORMAT(1X,'UNIT ',I3,' : ',A80)
!
!L
!L 0. Title line
!L
      IF(ICOUNT  /=  0) THEN
        WRITE(6,5) RUN_ID(1:4),RUN_ID(5:5),ICOUNT
      ELSE
        WRITE(6,10)  RUN_ID(1:4),RUN_ID(5:5)
      ENDIF
!
!L
!L 1. History file type
!L
      WRITE(6,100)RUN_HIST_TYPE
!L
!L 2. Run progress details
!L
      WRITE(6,150)RUN_TYPE,                                             &
     &            RUN_INDIC_OP,                                         &
     &            RUN_COMPCODE,                                         &
     &            (H_STEPim(J),J=1,4),                                  &
     &            (H_GROUPim(J),J=1,4),                                 &
     &            (MODEL_DATA_TIME(J),J=1,6),                           &
     &            (END_DUMPim(J),J=1,4),                                &
     &            (RESTARTim(J),J=1,4),                                 &
     &            (OFFSET_DUMPSim(J),J=1,4),                            &
     &            RUN_OCEAN_FIRST
!L
!L 3. Mean Processing details
!L
      WRITE(6,200)(MEAN_NUMBERim(J),J=1,4),                             &
     &            RUN_LAST_MEAN,                                        &
     &            RUN_MEANCTL_RESTART,                                  &
     &            ((RUN_MEANCTL_INDICim(I,J),J=1,4),I=1,4),             &
     &            (MEAN_OFFSETim(J),J=1,4)
!CC  *           ,RUN_MEANS_TO_DO
!L
!L 4. Job resubmit details
!L
      WRITE(6,300)RUN_JOB_NAME,                                         &
     &            RUN_RESUBMIT,                                         &
     &            (RUN_RESUBMIT_TARGET(J),J=1,6),                       &
     &            RUN_RESUBMIT_Q,RUN_RESUBMIT_TIME,                     &
     &            RUN_RESUBMIT_CPU,RUN_RESUBMIT_MEMORY,                 &
     &            RUN_RESUBMIT_PRTY
!
      IF(LONG)THEN
!L
!L 5. Logical/physical file associations
!L
      WRITE(6,850)
      DO I=1,NUNITS
        WRITE(6,851)I,MODEL_FT_UNIT(I)
      ENDDO

      ENDIF ! LONG
!
!L 6. Return
!L
      RETURN
      END SUBROUTINE PRINTHST
#endif
