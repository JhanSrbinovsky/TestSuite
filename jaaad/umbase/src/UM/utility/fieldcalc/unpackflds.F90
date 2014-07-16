#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines for packing and unpacking


!=======================================================================

!=======================================================================

SUBROUTINE UnPackFlds( PackType,        &  ! in
                       NumCols,         &  ! in
                       NumRows,         &  ! in
                       Field,           &  ! inout
                       ErrorStatus )       ! inout

! Description:
!
! Method:
!
! Owner: Dave Robinson
!
! History:
! Version Date     Comment
! ------- ----     -------
! 1.0     02/05/03 Original Code.  Sara James
! 6.0     12/09/03 Code implemented into UM. Dave Robinson
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  LenWord,                &
  PP_Field_type,          &
  PP_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: PackType      ! The type of packing used
INTEGER, INTENT(IN) :: NumCols
INTEGER, INTENT(IN) :: NumRows

TYPE(PP_Field_type), INTENT(INOUT) :: Field
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "UnPackFlds"
LOGICAL, PARAMETER :: Expand = .FALSE.

! Local Variables:
INTEGER :: FldSize
INTEGER :: idum                     ! Dummy integer
INTEGER :: NumWords                 ! Size of compressed field
CHARACTER(LEN=80) :: ErrMessage     ! Error message returned by COEX

! Local Arrays:
REAL :: WorkArray(NumCols, NumRows) ! array used for un_packing

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

NumWords = Field % Hdr % LBLREC
FldSize  = NumCols * NumRows

IF ( PackType == 1 ) THEN
 ! WGDOS packing
! DEPENDS ON: coex
  CALL COEX( WorkArray,  FldSize,  Field % RData, &
             NumWords,   NumCols,  NumRows,       &
             idum,       idum,     Expand,        &
             Field % Hdr % BMDI,   LenWord,       &
             ErrorStatus, ErrMessage)
  IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
    CALL EReport( RoutineName, StatusWarning, &
                  "COEX: " // ErrMessage )
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF

ELSE IF ( PackType == 3 ) THEN
 !  GRIB packing
  !CALL DEGRIB( Field % RData, WorkArray, FldSize, NumWords,    &
  !             Field % Hdr, Field % Hdr % BMDI, FldSize, LenWord)
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning, &
                "GRIB packing not supported." )
  ErrorStatus = StatusWarning
  GO TO 9999

ELSE IF ( PackType == 4 ) THEN
 ! Run length encoded data
  !CALL RUNLEN_DECODE( WorkArray, FldSize, Field % RData, NumWords,  &
  !                    Field % Hdr % BMDI, ErrorStatus, ErrMessage)
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning, &
                "RUNLEN packing not supported." )
  ErrorStatus = StatusWarning
  GO TO 9999

ELSE
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning, &
                "UNPACK - PackType not recognised." )
  ErrorStatus = StatusWarning
  GO TO 9999

ENDIF

Field % RData = WorkArray
Field % Hdr % LBPACK = Field % Hdr % LBPACK - PackType ! Data unpacked
Field % Hdr % LBLREC = FldSize
Field % Hdr % LBUser1 = 1  ! Data type is now real
Field % Hdr % BACC = 0.0

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE UnPackFlds

#endif
