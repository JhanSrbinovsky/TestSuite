#if defined(CONTROL) || defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Defines submodel and section/version configuration
!
! Subroutine Interface:

!- End of Subroutine code ----------------------------------------------

! Function Interface:

      INTEGER FUNCTION CTOITST(CHAR)
      CHARACTER*1  CHAR
!
      IF(CHAR == '0'.OR.CHAR == ' ') THEN
        CTOITST=0
      ELSE IF(CHAR == '1') THEN
        CTOITST=1
      ELSE IF(CHAR == '2') THEN
        CTOITST=2
      ELSE IF(CHAR == '3') THEN
        CTOITST=3
      ELSE IF(CHAR == '4') THEN
        CTOITST=4
      ELSE IF(CHAR == '5') THEN
        CTOITST=5
      ELSE IF(CHAR == '6') THEN
        CTOITST=6
      ELSE IF(CHAR == '7') THEN
        CTOITST=7
      ELSE IF(CHAR == '8') THEN
        CTOITST=8
      ELSE IF(CHAR == '9') THEN
        CTOITST=9
      ELSE IF(CHAR == 'A') THEN
        CTOITST=10
      ELSE IF(CHAR == 'B') THEN
        CTOITST=11
      ELSE IF(CHAR == 'C') THEN
        CTOITST=12
      ELSE IF(CHAR == 'D') THEN
        CTOITST=13
      ELSE IF(CHAR == 'E') THEN
        CTOITST=14
      ELSE IF(CHAR == 'F') THEN
        CTOITST=15
      ELSE IF(CHAR == 'G') THEN
        CTOITST=16
      ELSE IF(CHAR == 'H') THEN
        CTOITST=17
      ELSE IF(CHAR == 'I') THEN
        CTOITST=18
      ELSE IF(CHAR == 'J') THEN
        CTOITST=19
      ELSE IF(CHAR == 'K') THEN
        CTOITST=20
      ELSE IF(CHAR == '#') THEN
        CTOITST=0
      ELSE
        WRITE(6,*)'SETMODL: UNEXPECTED SECTION CHOICE  ',CHAR
      END IF

      END FUNCTION CTOITST
#endif
