#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines for packing and unpacking

SUBROUTINE Pack_Single( PackAcc,     &  ! in
                        PackType,    &  ! in
                        XpndField,   &  ! in
                        CompField,   &  ! inout
                        ErrorStatus )   ! inout

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
REAL, INTENT(IN)    :: PackAcc
CHARACTER(LEN=*), INTENT(IN) :: PackType
TYPE(PP_Field_type), INTENT(IN) :: XpndField

TYPE(PP_Field_type), INTENT(INOUT) :: CompField
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Pack_Single"
LOGICAL, PARAMETER :: Compress = .TRUE.

! Local variables:
INTEGER :: PackCode
INTEGER :: NumWords
INTEGER :: Num32BitWds
CHARACTER(LEN=80) :: ErrMessage

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

CompField % Hdr = XpndField % Hdr
IF ( ASSOCIATED( Compfield % RData ) ) THEN
  DEALLOCATE( CompField % RData )
  NULLIFY( CompField % RData )
END IF

IF( PackType == "WGDOS" ) THEN
  PackCode = 1
  CompField % Hdr % BACC = PackAcc
  ALLOCATE( CompField % RData( CompField % Hdr % LBLREC, 1 ) )
! DEPENDS ON: coex
  CALL COEX( XpndField % RData,            &
             XpndField % Hdr % LBLREC,     &
             CompField % RData,            &
             CompField % Hdr % LBLREC,     &
             XpndField % Hdr % NumCols,    &
             XpndField % Hdr % NumRows,    &
             Num32BitWds,                  &
             INT(PackAcc),                 &
             Compress,                     &
             XpndField % Hdr % BMDI,       &
             LenWord,                      &
             ErrorStatus, ErrMessage )
  IF( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
    CALL EReport( RoutineName, StatusWarning, ErrMessage )
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF
  NumWords = ( Num32BitWds-1 + (LenWord/32) ) * 32/LenWord
  CompField % Hdr % LBLREC = NumWords
  CompField % Hdr % LBPACK = XpndField % Hdr % LBPACK + 1
ELSE IF( PackType == "CRAY32" ) THEN
  WRITE(6,*) "CRAY32 Packing not supported - leaving unpacked"
ELSE IF( PackType == "GRIB" ) THEN
  WRITE(6,*) "GRIB Packing not supported - leaving unpacked"
ELSE IF( PackType == "RUNLEN" ) THEN
  WRITE(6,*) "RUNLEN Packing not supported - leaving unpacked"
ELSE
  WRITE(6,*) "PackType not recognised - leaving unpacked"
END IF

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE Pack_Single

!=======================================================================

!=======================================================================


#endif
