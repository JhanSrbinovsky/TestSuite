#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL
!LL    Subroutine:
!LL    INITMEAN
!LL
!LL    Purpose:
!LL    To set up and check the input parameters for the
!LL    meaning subroutines
!LL
!LL    Tested under compiler:           Tested under OS version:
!LL    cft77                            UNICOS 5.1
!LL
!LL    Programming standard:
!LL    UM Doc Paper 3
!LL
!LL    Logical system components covered:
!LL    C5
!LL
!LL    Project tasks:
!LL    C5
!LL
!LL    External documentation:
!LL    On-line UM document C5 - Control of means calculations
!LL
!*L    Interface and arguments:
      SUBROUTINE INITMEAN (                                             &
#include "argduma.h"
     &           ISUBMODL,ICODE,CMESSAGE)
!
      IMPLICIT NONE

!*L Arguments

#include "parparm.h"
#include "typsize.h"
#include "typduma.h"

      INTEGER                                                           &
     &       ISUBMODL                                                   &
                            ! IN  Submodel identifier (eg. atmos,ocean)
     &      ,ICODE          ! OUT return code; successful=0, error> 0
!
      CHARACTER*80                                                      &
     &       CMESSAGE       ! OUT Error message if ICODE > 0

!
!      Common blocks
!
#include "cmaxsize.h"
#include "csubmodl.h"
#include "chsunits.h"
#include "ctime.h"
#include "chistory.h"
#include "ccontrol.h"
#include "cmeanctl.h"
#include "clookadd.h"
!
!      Local variables and arrays
!
      INTEGER                                                           &
     &       IFIND,I,                                                   &
                                  ! Loop counts
     &       NMEANS,                                                    &
                                  ! No. of means chosen (fixed)
     &       MEAN_REFTIME_DAYS,                                         &
                                  ! Reference time for period means (D)
     &       MEAN_REFTIME_SECS,                                         &
                                  ! Reference time for period means (s)
     &       MEAN_START_DAYS,                                           &
                                  ! Start time for model run (days)
     &       MEAN_START_SECS,                                           &
                                  ! Start time for model run (s)
     &       MEAN_OFFSET_DAYS,                                          &
                                  ! Offset from mean ref. time (days)
     &       MEAN_OFFSET_SECS,                                          &
                                  ! Offset from mean ref. time (s)
     &       MEAN_OFFSET_STEPS,                                         &
                                  ! Offset from mean ref. time (steps)
     &       MEAN_FREQ_DUMPS,                                           &
                                  ! Mean frequency in dumps (updated)
     &       HOUR_SECS,                                                 &
                                  ! No. of secs in an hour
     &       YEAR,                                                      &
                                  ! Local variables for year,
     &       MONTH,                                                     &
                                  ! month etc
     &       DAY,                                                       &
     &       HOUR,                                                      &
     &       MINUTE,                                                    &
     &       SECOND
!
      CMESSAGE=' '
!
      ICODE=0
!
      HOUR_SECS=3600
!
!      Define mode of use for means program
!
#if defined(ATMOS)
      IF (ISUBMODL == atmos_sm) THEN
        WRITE(6,*)'INITMEAN: ***** Called in ATMOSPHERIC mode *****'
      END IF
#endif
!
!      Check input parameters
!
      NMEANS=0
!
      DO IFIND = 1,4
#if defined(ATMOS)
        IF (ISUBMODL == atmos_sm) THEN
          IF (MEANFREQim(IFIND,atmos_sm) >  0) THEN
            NMEANS = NMEANS+1
          END IF
          IF (MEANFREQim(IFIND,atmos_sm) == 1 .OR.                      &
     &        MEANFREQim(IFIND,atmos_sm) <  0) THEN
            ICODE = 1
            CMESSAGE='INITMEAN: Invalid atmos mean frequency'
            WRITE(6,*) 'INITMEAN: MEANFREQ(',IFIND,atmos_sm,') set to ',&
     &                          MEANFREQim(IFIND,atmos_sm)
            GO TO 9999
          END IF
        END IF
#endif
      END DO ! IFIND   
      IF (ISUBMODL == atmos_sm) THEN
        MEAN_NUMBERim(ISUBMODL) = NMEANS
        IF (MEAN_NUMBERim(ISUBMODL) == 0) THEN
          WRITE(6,*) 'INITMEAN: No means requested'
          GO TO 9999
        END IF
      END IF
!
!      If means are to be created:
!      Establish whether an offset exists between the
!      reference time for means creation and the start
!      of the integration
!      N.B. In the case of a restart, this check is ignored
!
      IF (ISUBMODL == atmos_sm) THEN
        WRITE(6,*)' H_STEP=',H_STEPim(ISUBMODL)
        IF (H_STEPim(ISUBMODL) == 0) THEN
!
          YEAR   = MEAN_REFTIMEim(1,ISUBMODL)
          MONTH  = MEAN_REFTIMEim(2,ISUBMODL)
          DAY    = MEAN_REFTIMEim(3,ISUBMODL)
          HOUR   = MEAN_REFTIMEim(4,ISUBMODL)
          MINUTE = MEAN_REFTIMEim(5,ISUBMODL)
          SECOND = MEAN_REFTIMEim(6,ISUBMODL)
!
          IF(YEAR   == 0.AND.                                           &
     &       MONTH  == 0.AND.                                           &
     &       DAY    == 0.AND.                                           &
     &       HOUR   == 0.AND.                                           &
     &       MINUTE == 0.AND.                                           &
     &       SECOND == 0)THEN
!
            MEAN_OFFSETim(ISUBMODL) = MEAN_NUMBERim(ISUBMODL)
            WRITE(6,*)'INITMEAN: No offset specified for means creation'
!
          ELSE
!
! DEPENDS ON: time2sec
            CALL TIME2SEC(YEAR,MONTH,DAY,HOUR,MINUTE,SECOND             &
     &                   ,0,0,MEAN_REFTIME_DAYS,MEAN_REFTIME_SECS,      &
     &                    LCAL360)
!
            YEAR   = MODEL_BASIS_TIME(1)
            MONTH  = MODEL_BASIS_TIME(2)
            DAY    = MODEL_BASIS_TIME(3)
            HOUR   = MODEL_BASIS_TIME(4)
            MINUTE = MODEL_BASIS_TIME(5)
            SECOND = MODEL_BASIS_TIME(6)
!
! DEPENDS ON: time2sec
            CALL TIME2SEC(YEAR,MONTH,DAY,HOUR,MINUTE,SECOND             &
     &                   ,0,0,MEAN_START_DAYS,MEAN_START_SECS,          &
     &                    LCAL360)
!
! DEPENDS ON: tim2step
            CALL TIM2STEP(MEAN_START_DAYS-MEAN_REFTIME_DAYS,            &
     &                    MEAN_START_SECS-MEAN_REFTIME_SECS,            &
     &        STEPS_PER_PERIODim(ISUBMODL),SECS_PER_PERIODim(ISUBMODL), &
     &                    MEAN_OFFSET_STEPS)
            OFFSET_DUMPSim(ISUBMODL) = MEAN_OFFSET_STEPS/               &
     &                                  DUMPFREQim(ISUBMODL)
            WRITE(6,*)'INITMEAN: Offset set up for means creation'
      WRITE(6,*)' OFFSET_DUMPSim(ISUBMODL)=',OFFSET_DUMPSim(ISUBMODL)
          MEAN_OFFSETim(ISUBMODL) = 0
            DO I = 1,MEAN_NUMBERim(ISUBMODL)
              IF (I == 1) THEN
                MEAN_FREQ_DUMPS = MEANFREQim(I,ISUBMODL)
              ELSE
                MEAN_FREQ_DUMPS = MEAN_FREQ_DUMPS*MEANFREQim(I,ISUBMODL)
              END IF
            IF (MOD(OFFSET_DUMPSim(ISUBMODL),MEAN_FREQ_DUMPS) == 0) THEN
                MEAN_OFFSETim(ISUBMODL) = MEAN_OFFSETim(ISUBMODL)+1
            END IF
            END DO ! I
          END IF
        END IF
      END IF
!L
!L Set IBUFLEN to the largest data field present in the dump
!L   for dimensioning of IO buffers in ACUMPS and MEANPS.
!L
#if defined(ATMOS)
      IF (ISUBMODL == atmos_sm) THEN
        IBUFLEN(atmos_sm) = 1
        DO I=1,A_LEN2_LOOKUP
          IF (A_LOOKUP(LBLREC,I) >  IBUFLEN(atmos_sm))                  &
     &      IBUFLEN(atmos_sm) = A_LOOKUP(LBLREC,I)
        END DO
      END IF
#endif
 9999 CONTINUE
      RETURN
      END SUBROUTINE INITMEAN
#endif
