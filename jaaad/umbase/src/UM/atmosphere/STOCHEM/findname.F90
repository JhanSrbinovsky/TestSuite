#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE FINDNAME(filetype,filetype_2,submodel_id,prev_h,       &
     &  prev_d,filename,month,year,hour,day)
! ----------------------------------------------------------------------
!-
!- Purpose and Methods: To construct UM style filenames.
!-   Simplified version of GET_NAM for use with REDIST_STOCHEM
!-   Allows hours and days to be subtracted from current time.
!-
!- Inputs: See below,
!-   PREV_H,PREV_D must be zero or negative.
!-   MONTH and YEAR are optional arguments.
!-
!- Outputs: FILENAME
!-
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  5.0    11/05/99  Created.  C.E. Johnson
!  5.1    26/10/00  month, year variables now optional
!  6.1    20/10/04  No change.
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
!
! ----------------------------------------------------------------------
      IMPLICIT NONE
! ----------------------------------------------------------------------
      INTEGER,       INTENT(IN)    :: prev_h      ! Hours to subtract
      INTEGER,       INTENT(IN)    :: prev_d      ! Days to subtract
      INTEGER, INTENT(IN),OPTIONAL :: month,day,hour ! month (optional)
      INTEGER, INTENT(IN),OPTIONAL :: year        ! year (optional)
      CHARACTER(1),  INTENT(IN)    :: filetype    ! Code for file type
      CHARACTER(1),  INTENT(IN)    :: filetype_2  ! Code for file type
      CHARACTER(1),  INTENT(IN)    :: submodel_id ! Code for submodel
      CHARACTER(14), INTENT(OUT)   :: filename    ! File name

! Global variables (*CALLed COMDECKs etc...):

#include "chsunits.h"
#include "cmaxsize.h"

! Common blocks

#include "csubmodl.h"
#include "ctime.h"
#include "cntlall.h"

! Local variables

      INTEGER :: yyyy,mm,dd,hh                       ! time variables

      CHARACTER(1)                :: y_tens          ! tens year id
      CHARACTER(1)                :: y_units         ! units year id
      CHARACTER(1)                :: m               ! month id
      CHARACTER(1)                :: d               ! day id
      CHARACTER(1)                :: h               ! hour id
      CHARACTER(1)                :: separator       ! separator
      CHARACTER(1), DIMENSION(36) :: char_id         ! lookup characters
      CHARACTER(3), DIMENSION(12) :: month_3char(12) ! ditto
      CHARACTER(72)               :: cmessage        ! Error message

      char_id = (/ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',    &
     &             'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',    &
     &             'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',    &
     &             'u', 'v', 'w', 'x', 'y', 'z' /)

      month_3char = (/ 'jan', 'feb', 'mar', 'apr', 'may', 'jun',        &
     &                 'jul', 'aug', 'sep', 'oct', 'nov', 'dec' /)

      separator = '.'

! Current time

      IF (present(year)) THEN
        yyyy = year
      ELSE
        yyyy = i_year
      END IF
      IF (present(month)) THEN
        mm = month
      ELSE
        mm = i_month
      END IF
      IF (present(day)) THEN
        dd = day
      ELSE
        dd = i_day
      END IF
      IF (present(hour)) THEN
        hh = hour
      ELSE
        hh = i_hour
      END IF

! Set prior time if required

      IF (prev_h < 0 .or. prev_d < 0) CALL PRIOR_DATE

      y_tens  = CHAR_ID(mod(yyyy/10,36)+1)
      y_units = CHAR_ID(mod(yyyy,10)+1)
      m       = CHAR_ID(mm+1)
      d       = CHAR_ID(dd+1)
      h       = CHAR_ID(hh+1)

! Construct the filename

      IF (filetype_2 == 'a') THEN
! hourly files
        filename(1:4)  =expt_id
        filename(5:5)  =job_id
        filename(6:6)  =submodel_id
        filename(7:7)  =separator
        filename(8:8)  =filetype
        filename(9:9)  =filetype_2
        filename(10:10)=y_tens
        filename(11:11)=y_units
        filename(12:12)=m
        filename(13:13)=d
        filename(14:14)=h
      ELSEIF (filetype_2 == 'm') THEN
! stochem formatted files or pp files
        filename(1:4)  =expt_id
        filename(5:5)  =job_id
        filename(6:6)  =submodel_id
        filename(7:7)  =separator
        filename(8:8)  =filetype
        filename(9:9)  =filetype_2
        filename(10:10)=y_tens
        filename(11:11)=y_units
        filename(12:14)=month_3char(mm)
      ELSE
        cmessage = 'filetype error in FINDNAME'
! DEPENDS ON: ereport
        CALL EREPORT('FINDNAME',1,cmessage)
      ENDIF

 999  CONTINUE
      RETURN

      CONTAINS

! #####################################################################
      SUBROUTINE PRIOR_DATE
! ----------------------------------------------------------------------
!- Purpose:
!-   To create a revised date with the option of subtraction of
!-   hours and days.
!-   PREV_H and PREV_D must be negative or zero.
!-   Internal procedure to FINDNAME.
!-
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  5.0    11/05/99  Created.  C.E. Johnson
!  5.1    17/10/00  360/365 day version. C.E. Johnson
!  6.1    20/10/04  No change.
!
! ----------------------------------------------------------------------
      IMPLICIT NONE

!      INTEGER, INTENT(INOUT) :: yyyy         ! Current year
!      INTEGER, INTENT(INOUT) :: mm           ! Current month
!      INTEGER, INTENT(INOUT) :: dd           ! Current day
!      INTEGER, INTENT(INOUT) :: hh           ! Current hour
!      INTEGER, INTENT(IN)    :: prev_h       ! No of hours previous
!      INTEGER, INTENT(IN)    :: prev_d       ! No of days previous

      INTEGER                :: julian
      INTEGER, DIMENSION(12) :: daym

! DEPENDS ON: set_daym
      CALL SET_DAYM(daym,yyyy)

! Subtract hours
      IF (prev_h < 0) THEN
        hh = hh + prev_h
        IF (hh < 0) THEN
          dd = dd - 1
          hh = 24 + hh
          IF (dd <= 0) THEN
            mm = mm - 1
            IF (mm == 0) THEN
              mm = 12
              yyyy = yyyy - 1
            END IF
            dd = daym(mm) + DD
          END IF
        END IF
      END IF

! Subtract days by conversion to julian numbers
      IF (prev_d < 0) THEN
! DEPENDS ON: datetojulian
        CALL DATETOJULIAN(yyyy,mm,dd,julian)

        julian = julian + prev_d
        IF (julian < 1 ) THEN
          yyyy = yyyy - 1
! DEPENDS ON: set_daym
          CALL SET_DAYM(daym,yyyy)
          julian = julian + SUM(daym)
        END IF

! DEPENDS ON: juliantodate
        CALL JULIANTODATE(yyyy,mm,dd,julian)
      END IF

      RETURN

      END SUBROUTINE PRIOR_DATE
! #####################################################################

      END SUBROUTINE FINDNAME
#endif
