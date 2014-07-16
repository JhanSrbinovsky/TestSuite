#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE JULIANTODATE(YEAR,MONTH,DAY,JULIAN)
! ----------------------------------------------------------------------
!- Purpose:
!-   To convert Julian number to date
!-
!- Created:    C.E. Johnson 17/X/00
!- modified    C.E. Johnson 25/VII/01  - for vn5.2
! ----------------------------------------------------------------------

      INTEGER,       INTENT(IN)    :: JULIAN,YEAR
      INTEGER,       INTENT(OUT)   :: MONTH,DAY

      INTEGER, DIMENSION(12) :: DAYM        ! Days of month
      CHARACTER*72           :: CMESSAGE    ! Error message

! DEPENDS ON: set_daym
      CALL SET_DAYM(DAYM,YEAR)

      IF(JULIAN > SUM(DAYM) .OR. JULIAN < 1) THEN
        CMESSAGE='Error in julian number'
        WRITE(6,*) cmessage,' Julian: ',JULIAN,' Year: ',YEAR
! DEPENDS ON: ereport
        CALL EREPORT('JuliantoDate',1,cmessage)
      ENDIF
      DO MONTH=1,12
        IF (SUM(DAYM(1:MONTH)) >= JULIAN) EXIT
      ENDDO

      DO DAY=1,DAYM(MONTH)
        IF (SUM(DAYM(1:MONTH-1)) + DAY == JULIAN) EXIT
      ENDDO

      RETURN

      END SUBROUTINE JULIANTODATE
#endif
