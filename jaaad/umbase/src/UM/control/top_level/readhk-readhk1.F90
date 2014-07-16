#if ( defined(SETUP) || defined(COMB) || defined(CONTROL) ) \
 &&  !defined(FLUXPROC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: READHK------------------------------------------------
!LL
!LL  Purpose: Read the operational model houskeeping file and set up
!LL           operational model control variables
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.0
!LL
!LL  Author:   A.B.SANGSTER       Date:           22 January 1990
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  date
!LL   3.1   1/02/93 : Add comdeck CHSUNITS to def NUNITS for extra i/o
!LL
!LL   3.1  05/02/93    Portable Fortran unit no assigns
!LL                    Author: A. Dickinson    Reviewer: R. Stratton
!LL   3.3  08/02/94  Modify calls to TIME2SEC/SEC2TIME to output/input
!LL                  elapsed times in days & secs, for portability. TCJ
!LL   3.4  17/06/94  Argument LCAL360 added and passed to SEC2TIME,
!LL                                                      TIME2SEC
!LL                                                    S.J.Swarbrick
!LL   3.4  24/06/94  Use service routine TIME_DF correctly to increment
!LL                  HK by -analysis time to get model basis time. RTHB
!LL  3.5  11/05/95  Sub-models stage 1: History/control files. RTHBarnes
!LL  4.3  26/02/97  Correct error in argument list of SEC2TIME
!LL                 present since vn3.4.            RTHBarnes.
!LL  6.2  31/01/06  Added a check on the return code, icode. T.Edwards
!LL  6.4  09/01/07  Ensure not compiled for flux processing utilities
!LL                                                          P.Selwood
!LL
!LL  Programming standard: UM Doc Paper 3, draft version 3 (15/1/90)
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
      SUBROUTINE READHK                                                 &
     &         ( UNITHK,ICODE,CMESSAGE)
!
      IMPLICIT NONE
!
      INTEGER      UNITHK    ! In  - Operational Model housekeeping file
      INTEGER       ICODE    ! Out - Return code from routine
      CHARACTER*(*) CMESSAGE ! Out - Return message if failure occured
!*
!
!L Common blocks
!
#include "csubmodl.h"
#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"
!
!*L EXTERNAL subroutines called
      EXTERNAL TIME2SEC,GET_FILE,SEC2TIME,TIME_DF
!*
!
!  Local variables
!
      INTEGER YEAR,MONTH,DAY,HOUR ! Read but not used
      INTEGER     ITYPE           ! Indicator of operational run mode
!                                 ! 1 - Global main run
!                                 ! 2 - Global update run
!                                 ! 3 - Limited area run
      INTEGER     INDIC           ! Read but not used
      INTEGER     IDAY1,ISEC1     ! Days/seconds from calendar zero
      INTEGER     IDAY,ISEC       ! Days/seconds from calendar zero
      INTEGER     DEL_SEC         ! Increment in seconds
      CHARACTER*80 FILENAME
      CHARACTER*4 TYPE            ! Type of run
!                                 ! 'NRUN' -  New run
!                                 ! 'CRUN' -  Continuation
!                                 ! 'RRUN' -  Complete rerun
!L
!L 0. Initialise
!L
!       LCAL360 = .FALSE. ! already set in model; was not in qxcombine
!L
!L 1. Open rewind and read record
!L
      CALL GET_FILE(UNITHK,FILENAME,80,ICODE)
        OPEN(UNITHK,FILE=FILENAME,IOSTAT=ICODE)
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='READHK  : Failed in OPEN of input unit'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
        WRITE(6,*)'READHK  : Warning message on OPEN of input unit'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
      REWIND(UNITHK)
      READ(UNITHK,50,IOSTAT=ICODE)                                      &
     &                       YEAR,MONTH,DAY,HOUR,INDIC,ITYPE,TYPE
      IF(ICODE  <   0) CONTINUE
!
! Set MODEL_BASIS_TIME from HKfile time
!
! DEPENDS ON: time2sec
      CALL TIME2SEC(YEAR,MONTH,DAY,HOUR,0,0,                            &
     &              0,0,IDAY1,ISEC1,LCAL360)
     ! UM6.5 : MODEL_ANALYSIS_HRS replaced by MODEL_ANALYSIS_MINS
      DEL_SEC = -3600*REAL(MODEL_ANALYSIS_MINS)/60.0
! DEPENDS ON: time_df
      CALL TIME_DF(IDAY1,ISEC1,0,DEL_SEC,IDAY,ISEC)
! DEPENDS ON: sec2time
      CALL SEC2TIME(IDAY,ISEC,0,0,                                      &
     &              MODEL_BASIS_TIME(1),MODEL_BASIS_TIME(2),            &
     &              MODEL_BASIS_TIME(3),MODEL_BASIS_TIME(4),            &
     &              MODEL_BASIS_TIME(5),MODEL_BASIS_TIME(6),            &
     &              IDAY1,LCAL360)
!
! Check for error
!
      IF(ICODE  >  0)THEN
        CMESSAGE='READHK  : Failed in READ of housekeeping file'
        GOTO 999
      ELSEIF(ICODE  <   0)THEN
      WRITE(6,*)'READHK  : Warning message on READ of housekeeping file'
        WRITE(6,*)'IOSTAT= ',ICODE
      ENDIF
!
  50  FORMAT(I4,5I2,A4)
!L
!L 2. Set appropriate history variables
!L
      RUN_INDIC_OP=ITYPE
      IF(TYPE  ==  'NRUN')RUN_TYPE='Setup'
      IF(TYPE  ==  'CRUN')RUN_TYPE='Continue'
      IF(TYPE  ==  'RRUN')RUN_TYPE='Rerun'
 100  CONTINUE
 999  CONTINUE
!L
!L 3. Close and return
!L
      CLOSE(UNITHK)
      RETURN
      END SUBROUTINE READHK
#endif
