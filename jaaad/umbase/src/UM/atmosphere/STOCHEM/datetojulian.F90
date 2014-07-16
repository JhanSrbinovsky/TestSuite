#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE DATETOJULIAN(YEAR,MONTH,DAY,JULIAN)
! ----------------------------------------------------------------------
!- Purpose:
!-   To convert date to Julian number
!-
!- Created:    C.E. Johnson 17/X/00
!- modified    C.E. Johnson 25/VII/01  - for vn5.2
! ----------------------------------------------------------------------

      INTEGER,  INTENT(IN) :: YEAR,MONTH,DAY
      INTEGER, INTENT(OUT) :: JULIAN

      INTEGER, DIMENSION(12) :: DAYM        ! Days of month
      CHARACTER*72           :: CMESSAGE    ! Error message

! DEPENDS ON: set_daym
      CALL SET_DAYM(DAYM,YEAR)
      IF(DAY > DAYM(MONTH)) THEN
        CMESSAGE='ERROR IN INPUT, Day > Days in month'
        WRITE(6,*) cmessage,' DMY: ',DAY,MONTH,YEAR
! DEPENDS ON: ereport
        CALL EREPORT('DatetoJulian',1,cmessage)
      ENDIF
      JULIAN=SUM(DAYM(1:MONTH-1))+DAY

      RETURN

      END SUBROUTINE DATETOJULIAN
#endif
