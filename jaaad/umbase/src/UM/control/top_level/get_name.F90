#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: GET_NAME -------------------------------------------------
!LL
!LL  Purpose: Generates an output file name of up to 14 characters using
!LL           the defined file naming convention, taking account of
!LL           file type, validity time, etc.
!LL           Obeys new filenaming convention introduced at version 2.7.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Author:   R A Stratton
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1  27/01/93 : correct error in 2 character months - change jan
!LL                   from jn to ja so that no clash with june.
!LL  3.1 2/02/93 : added comdeck CHSUNITS to define NUNITS for i/o
!LL
!LL   3.1  22/02/93 : Cater for filename changes for boundary datasets.
!LL                   FILETYPE=u-z for LAM areas 1-6. D. Robinson
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.2  15/04/93  Correct Y_HUNDS in Absolute_long convention (TCJ).
!LL   3.2  08/07/93  Correct 1 T/S offset in reinit'ed pp names  (TCJ).
!LL   3.3  08/02/94  Modify calls to TIME2SEC/SEC2TIME to output/input
!LL                  elapsed times in days & secs, for portability. TCJ
!LL   3.4  17/06/94  Argument LCAL360 added and passed to SEC2TIME
!LL                                                        S.J.Swarbrick
!LL  4.1  30/07/96  Introduce Wave sub-model.  M Holt
!LL  4.4  11/07/97  Allow character filenames for PP files
!LL                 reinitialised on real month boundaries.  M Gallani
!LL  4.5  29/07/98  New naming convention for reinitialised boundary
!LL                 files. D. Robinson.
!LL  5.2  13/10/00  Dumpname fixed when MODEL_STATUS=SCS. Stuart Bell
!LL  5.2  16/01/01  New naming convention for reinitialized Macro
!LL                 files. (R.Hatcher)
!LL  5.2  19/07/99  Correct new filenaming convention for reinitialised
!LL                 boundary files. D. Robinson.
!LL  6.2  27/10/05  Add "Sub-hourly" filename convention. R Barnes.
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: S51
!LL
!LL  Project task: S51
!LL
!LL  External documentation: UM documentation paper 7 - Filenaming
!LL                          conventions for the Unified Model
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE GET_NAME(EXPT_ID,JOB_ID,ISUBMODEL,MEANLEV,TOGGLE,      &
     &          REINIT_STEPS,FILETYPE,LETTER_3, MODEL_STATUS,           &
     &     TIME_CONVENTION,ANALYSIS_HRS,FILENAME,ICODE,CMESSAGE,        &
     &     LCAL360)
!
      IMPLICIT NONE
      LOGICAL LCAL360
!
      CHARACTER*4   EXPT_ID     ! IN  - Experiment ident or alias
      CHARACTER*1   JOB_ID      ! IN  - Job ident within experiment
      INTEGER       ISUBMODEL   ! IN  - Submodel indicator
      INTEGER       MEANLEV     ! IN  - Mean level indicator
      INTEGER       TOGGLE      ! IN  - Alternately 1/2 for partial sums
      REAL          ANALYSIS_HRS! IN  - Hrs from basis time to analysis
                                ! UM6.5 - Allow for fractional ANALYSIS_HRS - 
                                !          change from INTEGER to REAL
      INTEGER       REINIT_STEPS! IN  - timesteps between file reinit-
!                                      ialisation for non-mean pp files
!                                      or -ve for Gregorian reinit.

      CHARACTER*1   FILETYPE    ! IN  - Code for file type
      CHARACTER*1   LETTER_3    ! IN  - character for use in position 9
!                                       of non-mean pp files.
      CHARACTER*14  MODEL_STATUS! IN  - Operational/NonOperational
!
      CHARACTER*17  TIME_CONVENTION ! IN  - Relative/Timestep/
!                                       Absolute_standard/Absolute_long/
!                                       Absolute_short
      CHARACTER*14  FILENAME    ! OUT - Generated file name
      INTEGER ICODE             ! OUT - Error return code
      CHARACTER*80 CMESSAGE
!
!*----------------------------------------------------------------------
!  Common blocks
!
#include "cmaxsize.h"
#include "csubmodl.h"
#include "cntlgen.h"
#include "ctime.h"
!
! External subroutines called
!
      EXTERNAL SEC2TIME,DAY2CHAR,STP2TIME
!
!  Local variables
!
      INTEGER                                                           &
     & YYYY,MM,DD,HH,IMIN,ISEC   ! Current time values for filename
      INTEGER     COUNT,                                                &
                                 ! Counter for steps or hours
     &            DAYNO,                                                &
                                 ! day number
     &            DAYS,                                                 &
                                 ! Number of days for period
     &            HOURS,                                                &
                                 ! Number of hours for period
     &            I,                                                    &
                                 ! loop counter
     &            STEPS          ! number of steps
      INTEGER     END_DAYS       ! number of whole days from run start
      INTEGER     END_SECONDS    ! number of extra secs from run start
      INTEGER     MON            ! month for mean period
      INTEGER     A_STEPS_PER_HR    ! steps per hour for atmos sub-model
!
      CHARACTER*2 QW             ! Operational file prefix
!
      CHARACTER*1 SUBMODEL_ID    ! character for model a or o
      CHARACTER*1 FILETYPE_2     ! letter after FILETYPE in name
      CHARACTER*1 MEAN_PERIOD(4) ! default letter for mean period
!
      CHARACTER*1 Y_HUNDS        ! Character year identifier (hundreds)
      CHARACTER*1 Y_TENS         ! Character year identifier (tens)
      CHARACTER*1 Y_UNITS        ! Character year identifier (units)
      CHARACTER*1 M              ! Character month identifier
      CHARACTER*1 D              ! Character day-of-month identifier
      CHARACTER*1 H              ! Character hour identifier
      CHARACTER*1 HUNDREDS       ! Character for hundreds counter
      CHARACTER*1 TENS           ! Character for tens counter
      CHARACTER*1 UNITS          ! Character for units counter
      CHARACTER*1 DECI           ! Character for tenths (=mins)
      CHARACTER*1 CHAR_ID(36)    ! Valid characters for above (lookup)
      CHARACTER*1 SEPARATOR      ! character used as separator in name
      CHARACTER*1 STYLE          ! style of date in filename
      CHARACTER*3 CDAYNO         ! character day number
      CHARACTER*3 MONTH_3CHAR(12)! 3 character month identifier
      CHARACTER*2 MONTH_2CHAR(12)! 2 character month identifier
      CHARACTER*3 SEASON_3CHAR(12)! 3 character season identifier
      CHARACTER*2 SEASON_2CHAR(12)! 2 character season identifier
!
      DATA QW / 'qw'/
      DATA MEAN_PERIOD / '1', '2', '3', '4' /
      DATA CHAR_ID/ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',   &
     &              'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',   &
     &              'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',   &
     &              'u', 'v', 'w', 'x', 'y', 'z' /
      DATA MONTH_3CHAR/ 'jan', 'feb', 'mar', 'apr', 'may', 'jun',       &
     &                  'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/
      DATA MONTH_2CHAR/ 'ja', 'fb', 'mr', 'ar', 'my', 'jn',             &
     &                  'jl', 'ag', 'sp', 'ot', 'nv', 'dc'/
      DATA SEASON_3CHAR/ 'ndj', 'djf', 'jfm', 'fma', 'mam', 'amj',      &
     &                   'mjj', 'jja', 'jas', 'aso', 'son', 'ond'/
      DATA SEASON_2CHAR/ 'nj', 'df', 'jm', 'fa', 'mm', 'aj',            &
     &                   'mj', 'ja', 'js', 'ao', 'sn', 'od'/
!L
!L----------------------------------------------------------------------
!L 1. Determine submodel id - (used in filechar 6 or 7 if operational)
!L
      IF (ISUBMODEL == ATMOS_SM) THEN
       SUBMODEL_ID= 'a'
      ELSE IF (ISUBMODEL == OCEAN_SM) THEN
       SUBMODEL_ID= 'o'
      ELSE IF (ISUBMODEL == WAVE_SM) THEN
       SUBMODEL_ID= 'w'
      ELSE
       ICODE=2
       CMESSAGE='GET_NAME: Illegal sub-model specified'
       GOTO 999
      ENDIF
      IF (ISUBMODEL == ATMOS_SM) THEN
! 1.1 Compute steps per hour for atmosphere sub_model
      A_STEPS_PER_HR = 3600*STEPS_PER_PERIODim(a_im)/                   &
     &                       SECS_PER_PERIODim(a_im)
      ENDIF
!
!L----------------------------------------------------------------------
!L 2. Determine style filename and separator
!L
      IF (FILETYPE /= 's') THEN
!L
!L 2.1 Relative time convention
!L
        IF (TIME_CONVENTION == 'Relative        ') THEN
          SEPARATOR='_'
          STYLE='A'
          IF (ISUBMODEL == atmos_sm) THEN ! atmosphere
          COUNT = STEPim(a_im) / A_STEPS_PER_HR - ANALYSIS_HRS
          ELSE IF (ISUBMODEL == ocean_sm) THEN ! ocean
            COUNT = STEPim(o_im) * SECS_PER_PERIODim(o_im) /            &
     &                 STEPS_PER_PERIODim(o_im) / 3600                  &
     &             -ANALYSIS_HRS
          ELSE IF (ISUBMODEL == wave_sm) THEN ! WAVE
            COUNT = STEPim(w_im) * SECS_PER_PERIODim(w_im) /            &
     &                 STEPS_PER_PERIODim(w_im) / 3600                  &
     &             -ANALYSIS_HRS
          ENDIF
!L
!L 2.2 Step time convention
!L
        ELSE IF (TIME_CONVENTION == 'Timestep         ') THEN

          SEPARATOR='_'
          STYLE='A'
          IF (ISUBMODEL == atmos_sm) THEN ! Atmosphere
            COUNT = STEPim(a_im)
          ELSE if (ISUBMODEL == ocean_sm) then ! Ocean
            COUNT = STEPim(o_im)
          ELSE if (ISUBMODEL == wave_sm) then ! WAVE
            COUNT = STEPim(w_im)
          ENDIF
!L
!L 2.3 Absolute time convention -standard version
!L
        ELSE IF (TIME_CONVENTION == 'Absolute_standard') THEN
          SEPARATOR='.'
          STYLE='B'
!L
!L 2.4 Absolute time convention - short
!L
        ELSE IF (TIME_CONVENTION == 'Absolute_short   ') THEN
          SEPARATOR='-'
          STYLE='B'
!L
!L 2.5 Absolute time convention - long
!L
        ELSE IF (TIME_CONVENTION == 'Absolute_long    ') THEN
          SEPARATOR='@'
          STYLE='B'
!L
!L 2.6 Sub-hourly filenaming time convention
!L
        ELSE IF (TIME_CONVENTION == 'Sub-hourly      ') THEN
          SEPARATOR='_'
          STYLE='A'
          IF (ISUBMODEL == atmos_sm) THEN ! atmosphere
            COUNT = STEPim(a_im) * SECS_PER_STEPim(a_im)
      write(6,*)'GET_NAM - COUNT ',a_im,COUNT,STEPim(a_im),             &
     & SECS_PER_STEPim(a_im)
          ELSE IF (ISUBMODEL == ocean_sm) THEN ! ocean
            COUNT = STEPim(o_im) * SECS_PER_STEPim(o_im)
          ELSE IF (ISUBMODEL == wave_sm) THEN ! WAVE
            COUNT = STEPim(w_im) * SECS_PER_STEPim(w_im)
          ENDIF
!
!   UM6.5 - ANALYSIS_HRS changed to REAL - 
!                  requires recoding of sub-hourly file-naming.  
             COUNT = 100*(REAL(COUNT)/3600.0 - ANALYSIS_HRS)
             COUNT = (COUNT/100)*100 + (REAL(MOD(COUNT,100))/100.0)*60

      write(6,*)'GET_NAM - hhmm ',COUNT,ANALYSIS_HRS
        ELSE
            ICODE=1
            CMESSAGE='GET_NAME: Illegal TIME_CONVENTION specified'
            GOTO 999
        ENDIF
!L----------------------------------------------------------------------
!L
!L 3.0 work out encoding of date time and filetype_2
!L
        IF (STYLE == 'A') THEN
          IF (COUNT <  0) THEN
            COUNT=-COUNT
            FILETYPE_2='z'
          IF (TIME_CONVENTION == 'Sub-hourly      ') THEN
!!!            FILETYPE_2='h'
            IF (COUNT  >=  36000) THEN
              ICODE = COUNT
              CMESSAGE='GET_NAME: COUNT too big for '//                 &
     &                 'sub-hourly filenaming convention'
              GOTO 999
            END IF
            HUNDREDS = CHAR_ID(MOD(COUNT/1000,36)+1)
            TENS =     CHAR_ID(MOD(COUNT/100, 10)+1)
            UNITS =    CHAR_ID(MOD(COUNT/10,  10)+1)
            DECI =     CHAR_ID(MOD(COUNT,     10)+1)
          ELSE
            IF (COUNT  >=  3600) THEN
              ICODE = COUNT
              CMESSAGE='GET_NAME: COUNT too big for '//                 &
     &                 'hourly or timestep filenaming convention'
              GOTO 999
            END IF
            HUNDREDS = CHAR_ID(MOD(COUNT/100,36)+1)
            TENS =     CHAR_ID(MOD(COUNT/10 ,10)+1)
            UNITS =    CHAR_ID(MOD(COUNT,    10)+1)
          END IF
          ELSE
            IF (FILETYPE == 'p') THEN
              FILETYPE_2=LETTER_3
        ELSE IF (FILETYPE == 'b') THEN   !  Boundary File
          FILETYPE_2=LETTER_3
            ELSE IF (FILETYPE == 'c') THEN   ! Macro File
              FILETYPE_2=LETTER_3
            ELSE
              FILETYPE_2='a'
            ENDIF
          IF (TIME_CONVENTION == 'Sub-hourly      ') THEN
!!!            FILETYPE_2='h'
            IF (COUNT  >=  36000) THEN
              ICODE = COUNT
              CMESSAGE='GET_NAME: COUNT too big for '//                 &
     &                 'sub-hourly filenaming convention'
              GOTO 999
            END IF
            HUNDREDS = CHAR_ID(MOD(COUNT/1000,36)+1)
            TENS =     CHAR_ID(MOD(COUNT/100, 10)+1)
            UNITS =    CHAR_ID(MOD(COUNT/10,  10)+1)
            DECI =     CHAR_ID(MOD(COUNT,     10)+1)
          ELSE
            IF (COUNT  >=  3600) THEN
              ICODE = COUNT
              CMESSAGE='GET_NAME: COUNT too big for '//                 &
     &                 'hourly or timestep filenaming convention'
              GOTO 999
            END IF
            HUNDREDS = CHAR_ID(MOD(COUNT/100,36)+1)
            TENS =     CHAR_ID(MOD(COUNT/10 ,10)+1)
            UNITS =    CHAR_ID(MOD(COUNT,    10)+1)
          END IF
          ENDIF
        ELSE IF (STYLE == 'B') THEN   ! some sort of absolute time
!
! Current date time is
!
          YYYY=I_YEAR
          MM  =I_MONTH
          DD  =I_DAY
          HH  =I_HOUR
          DAYNO = I_DAY_NUMBER
!
! Instantaneous files
!
          IF (MEANLEV == 0) THEN
!  dumps
            IF (FILETYPE == 'd') THEN
              FILETYPE_2 = 'a'
            ENDIF
!
! Work out reintialisation period for pp and boundary files.
! Note assumes reinitialisation period is whole number of hours. This
! is not strictly true but is probably ok for this purpose.
!
         IF (  FILETYPE == 'p'                                          &
                                 !  PP File
     &    .or. FILETYPE == 'b'                                          &
                                 !  Boundary File
     &    .or. FILETYPE == 'c'                                          &
                                 !  Macro File
     & ) THEN
              IF (ISUBMODEL == atmos_sm) then
                HOURS = REINIT_STEPS/A_STEPS_PER_HR
              ELSE IF (ISUBMODEL == ocean_sm) THEN
                HOURS = REINIT_STEPS*SECS_PER_PERIODim(o_im)            &
     &                             /(STEPS_PER_PERIODim(o_im)*3600)
              ELSE IF (ISUBMODEL == wave_sm) THEN
                HOURS = REINIT_STEPS*SECS_PER_PERIODim(w_im)            &
     &                             /(STEPS_PER_PERIODim(w_im)*3600)
              ENDIF
              if (REINIT_STEPS <  0) then ! Gregorian reinitialisation
                HOURS=720 ! dummy: could be anything divisible by 24
              endif
!   do further checks if multiple of 1 day
              IF (MOD(HOURS,24) == 0) THEN  ! whole days in period, or
                DAYS=HOURS/24               ! Gregorian reinit.
                IF (FILETYPE == 'b'                                     &
                                      !  Boundary File
     &           .or. FILETYPE == 'c'                                   &
                                      !  Macro File
     &        ) THEN
                  FILETYPE_2=LETTER_3
                ENDIF
                IF (DD == 1 .and. TIME_CONVENTION /= 'Absolute_short   '&
     &              .and. (DAYS == 30 .or. REINIT_STEPS <  0)) then
!  Original code didn't allow style=C for 3-month reinit files but new
!  code (Gregorian reinit) does, at least in section 3.0.
                  STYLE='C'                ! month in characters
                ENDIF
              ELSE
                IF (FILETYPE == 'b'                                     &
                                      !  Boundary File
     &           .or. FILETYPE == 'c'                                   &
                                      !  Macro File
     &        ) THEN
                  FILETYPE_2=LETTER_3
                ENDIF
              ENDIF
!
! For instantaneous pp file need to work out end of period as call to
!   this routine occurrs on the first output timestep.
!
              IF (FILETYPE == 'p') THEN
                FILETYPE_2 = LETTER_3
                IF (STYLE /= 'C') THEN
                  IF (ISUBMODEL == atmos_sm) THEN
                    IF (STEPim(a_im) == 0) THEN
                     STEPS = REINIT_STEPS
                    ELSE
                     STEPS = STEPim(a_im) + REINIT_STEPS - 1
                    ENDIF
! DEPENDS ON: stp2time
                    CALL STP2TIME(STEPS,                                &
     &                            A_STEPS_PER_HR*24,86400,              &
     &                            END_DAYS,END_SECONDS)
                  ELSE IF (ISUBMODEL == ocean_sm) THEN
                    IF (STEPim(o_im) == 0) THEN
                     STEPS = REINIT_STEPS
                    ELSE
                     STEPS = STEPim(o_im) + REINIT_STEPS - 1
                    ENDIF
! DEPENDS ON: stp2time
                    CALL STP2TIME(STEPS,                                &
     &                STEPS_PER_PERIODim(o_im),SECS_PER_PERIODim(o_im), &
     &                            END_DAYS,END_SECONDS)
                  ELSE IF (ISUBMODEL == wave_sm) THEN
                    IF (STEPim(w_im) == 0) THEN
                     STEPS = REINIT_STEPS
                    ELSE
                     STEPS = STEPim(w_im) + REINIT_STEPS - 1
                    ENDIF
! DEPENDS ON: stp2time
                    CALL STP2TIME(STEPS,                                &
     &                STEPS_PER_PERIODim(w_im),SECS_PER_PERIODim(w_im), &
     &                            END_DAYS,END_SECONDS)
                  ENDIF
! DEPENDS ON: sec2time
                  CALL SEC2TIME(END_DAYS,END_SECONDS,                   &
     &                          BASIS_TIME_DAYS,BASIS_TIME_SECS,        &
     &                          YYYY,MM,DD,HH,IMIN,ISEC,DAYNO,LCAL360)
                ENDIF
              ENDIF
            ENDIF

          ELSE         !  MEANS
!
!L determine if special mean period
            IF (ISUBMODEL == ATMOS_SM) THEN
              HOURS=DUMPFREQim(a_im)/A_STEPS_PER_HR
              DO I=1,MEANLEV
                HOURS=HOURS*MEANFREQim(I,a_im) !hours per meaning period
              ENDDO
            ELSE IF (ISUBMODEL == OCEAN_SM) THEN
              HOURS=DUMPFREQim(o_im)*SECS_PER_PERIODim(o_im)            &
     &                        /(3600*STEPS_PER_PERIODim(o_im))
              DO I=1,MEANLEV
                HOURS=HOURS*MEANFREQim(I,o_im) !hours per meaning period
              ENDDO
            ELSE IF (ISUBMODEL == WAVE_SM) THEN
              HOURS=DUMPFREQim(w_im)*SECS_PER_PERIODim(w_im)            &
     &                        /(3600*STEPS_PER_PERIODim(w_im))
              DO I=1,MEANLEV
                HOURS=HOURS*MEANFREQim(I,w_im) !hours per meaning period
              ENDDO
            ENDIF
            IF (MOD(HOURS,24) == 0) THEN
              DAYS=HOURS/24
! DEPENDS ON: day2char
              CALL DAY2CHAR(DAYS,FILETYPE_2)
              IF (FILETYPE_2 == '0') THEN
                FILETYPE_2=MEAN_PERIOD(MEANLEV)
              ELSE if (FILETYPE_2 == 'm'.and.DD == 1                    &
     &              .AND.TIME_CONVENTION /= 'Absolute_short    ') THEN
                STYLE='C'      ! period starts at beginning of a month
                IF (MM == 1) THEN ! correct year if month december
                  YYYY=YYYY-1
                ENDIF
              ELSE if (FILETYPE_2 == 's'.and.DD == 1                    &
     &               .AND.TIME_CONVENTION /= 'Absolute_short    ') THEN
                STYLE='C'      ! period starts at beginning of aseason
              ENDIF
            ELSE
              FILETYPE_2=MEAN_PERIOD(MEANLEV)
            ENDIF
          ENDIF
!
          Y_UNITS = CHAR_ID(MOD(YYYY,10)+1)
          M = CHAR_ID(MM+1)
          D = CHAR_ID(DD+1)
          H = CHAR_ID(HH+1)
!
        ENDIF
!
      ELSE
! partial sum files - no date time information required
        SEPARATOR='_'
        FILETYPE_2 = MEAN_PERIOD(MEANLEV)
      ENDIF
!
!L----------------------------------------------------------------------
!L 3.1 Construct filename from the various components
!L
      FILENAME="              "
      IF (MODEL_STATUS == 'Operational'.OR.MODEL_STATUS == 'SCS') THEN
        FILENAME(1:2)  =QW
        FILENAME(3:6)  =EXPT_ID
        FILENAME(7:7)  =SUBMODEL_ID
        FILENAME(8:8)  =SEPARATOR
        FILENAME(9:9)  =FILETYPE
        FILENAME(10:10)=FILETYPE_2
        FILENAME(11:11)=HUNDREDS
        FILENAME(12:12)=TENS
        FILENAME(13:13)=UNITS
        IF (TIME_CONVENTION == 'Sub-hourly      ') THEN
          FILENAME(14:14)=DECI
        END IF
      ELSE
        FILENAME(1:4)  =EXPT_ID
        FILENAME(5:5)  =JOB_ID
        FILENAME(6:6)  =SUBMODEL_ID
        FILENAME(7:7)  =SEPARATOR
        FILENAME(8:8)  =FILETYPE
        FILENAME(9:9)  =FILETYPE_2
        IF (FILETYPE == 's') THEN
          IF (TOGGLE == 1) THEN
            FILENAME(10:10)='a'
          ELSE
            FILENAME(10:10)='b'
          ENDIF
        ELSE IF (STYLE == 'A') THEN
          FILENAME(10:10)=HUNDREDS
          FILENAME(11:11)=TENS
          FILENAME(12:12)=UNITS
          IF (TIME_CONVENTION == 'Sub-hourly      ') THEN
            FILENAME(13:13)=DECI
          END IF
        ELSE IF (STYLE == 'B') THEN
          IF (TIME_CONVENTION == 'Absolute_standard') THEN
!
! decades meansured relative to 1800
            Y_TENS  = CHAR_ID(MOD(YYYY/10,36)+1)
            FILENAME(10:10)=Y_TENS
            FILENAME(11:11)=Y_UNITS
            FILENAME(12:12)=M
            FILENAME(13:13)=D
            FILENAME(14:14)=H
          ELSE IF (TIME_CONVENTION == 'Absolute_long    ') THEN

! centuries  measured from 0 ie 1992  as j92, with wraparound at 3600
            Y_HUNDS = CHAR_ID(MOD(YYYY/100,36)+1)
            Y_TENS  = CHAR_ID((YYYY-(YYYY/100)*100)/10+1)
            FILENAME(10:10)=Y_HUNDS
            FILENAME(11:11)=Y_TENS
            FILENAME(12:12)=Y_UNITS
            FILENAME(13:13)=M
            FILENAME(14:14)=D
          ELSE IF (TIME_CONVENTION == 'Absolute_short   ') THEN
            FILENAME(10:10)=Y_UNITS
    1       FORMAT (I3)
    2       FORMAT ('0',I2)
    3       FORMAT ('00',I1)
            IF (DAYNO <  100) THEN
              IF (DAYNO <  10) THEN
                WRITE(CDAYNO,3) DAYNO
              ELSE
                WRITE(CDAYNO,2) DAYNO
              ENDIF
            ELSE
              WRITE(CDAYNO,1) DAYNO
            ENDIF
            FILENAME(11:13)= CDAYNO
            FILENAME(14:14)=H
          ENDIF
        ELSE  ! style C - Character date
          IF (TIME_CONVENTION == 'Absolute_standard') THEN
!
! decades meansured relative to 1800
            Y_TENS  = CHAR_ID(MOD(YYYY/10,36)+1)
            FILENAME(10:10)=Y_TENS
            FILENAME(11:11)=Y_UNITS
            IF (MEANLEV == 0) THEN
              FILENAME(12:14) = MONTH_3CHAR(MM)
            ELSE ! means date routine called is at beginning of next
!                   period
             MON=MM-1
             IF (MON == 0) THEN
               MON = 12
             ENDIF
             IF (FILETYPE_2 == 'm') THEN
               FILENAME(12:14) = MONTH_3CHAR(MON)
             ELSE
               FILENAME(12:14) = SEASON_3CHAR(MON)
             ENDIF
            ENDIF
          ELSE IF (TIME_CONVENTION == 'Absolute_long    ') THEN

! centuries  measured from 0 ie 1992  as j92, with wraparound at 3600
            Y_HUNDS = CHAR_ID(MOD(YYYY/100,36)+1)
            Y_TENS  = CHAR_ID((YYYY-(YYYY/100)*100)/10+1)
            FILENAME(10:10)=Y_HUNDS
            FILENAME(11:11)=Y_TENS
            FILENAME(12:12)=Y_UNITS
            IF (MEANLEV == 0) THEN
              FILENAME(13:14) = MONTH_2CHAR(MM)
            ELSE ! means date routine called is at beginning of next
!                   period
             MON=MM-1
             IF (MON == 0) THEN
               MON = 12
             ENDIF
             IF (FILETYPE_2 == 'm') THEN
               FILENAME(13:14) = MONTH_2CHAR(MON)
             ELSE
               FILENAME(13:14) = SEASON_2CHAR(MON)
             ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
!
 999  CONTINUE
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE GET_NAME
#endif
