#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to set up a new UM header in Fieldcalc

SUBROUTINE New_UMHdr( UMHdr,        &  ! in
                      MaxFlds,      &  ! in
                      NewUMHdr,     &  ! inout
                      ErrorStatus )    ! inout

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

USE IO_Mod, ONLY: &
  LenFixHd,       &
  UM_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal
IMPLICIT None

! Subroutine Arguments:
TYPE(UM_Header_type), INTENT(IN)    :: UMHdr
INTEGER, INTENT(IN) :: MaxFlds

TYPE(UM_Header_type), INTENT(INOUT) :: NewUMHdr
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "New_UMHdr"
CHARACTER(LEN=*), PARAMETER :: PPFormat = "p"

! Local Variables:
INTEGER :: WordAddress
CHARACTER(LEN=80) :: ErrMessage

! End of header --------------------------------------------------------

IF( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! DEPENDS ON: file_open
CALL File_Open ( NewUMHdr % UnitNum,                     &
                 NewUMHdr % FileNameEnv,                 &
                 LEN_TRIM(NewUMHdr % FileNameEnv),       &
                 1, 0, ErrorStatus)
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, &
                "Failed to open output file" )
  ErrorStatus = StatusFatal
  GO TO 9999
END IF

!----------------------------------------------------
! Use header provided as start point for new header
!----------------------------------------------------
WordAddress = 0   ! Start of file
! DEPENDS ON: setpos
CALL SETPOS ( NewUMHdr % UnitNum, WordAddress, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! DEPENDS ON: init_pp
CALL INIT_PP ( NewUMHdr % UnitNum,     &
               PPFormat,               &
               UMHdr % FixHd(151),     &
               UMHdr % FixHd(152),     &
               UMHdr % FixHd,          &
               UMHdr % IntC,           &
               UMHdr % RealC,          &
               UMHdr % LevDepC,        &
               UMHdr % RowDepC,        &
               UMHdr % ColDepC,        &
               LenFixHd,               &
               UMHdr % LenIntC,        &
               UMHdr % LenRealC,       &
               UMHdr % Len1LevDepC,    &
               UMHdr % Len2LevDepC,    &
               UMHdr % Len1RowDepC,    &
               UMHdr % Len2RowDepC,    &
               UMHdr % Len1ColDepC,    &
               UMHdr % Len2ColDepC,    &
               UMHdr % LenIntC,        &
               UMHdr % LenRealC,       &
               ErrorStatus, ErrMessage )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, ErrMessage )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

!-------------------------------------
! Read the header written by INIT_PP
!-------------------------------------
ALLOCATE( NewUMHdr % FixHd(LenFixHd) )

WordAddress = 0   ! Start of file
! DEPENDS ON: setpos
CALL SETPOS ( NewUMHdr % UnitNum, WordAddress, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! DEPENDS ON: read_flh
CALL READ_FLH ( NewUMHdr % UnitNum, & ! in
                NewUMHdr % FixHd,   & ! out
                LenFixHd,           & ! in
                ErrorStatus,        & ! out
                ErrMessage )          ! out
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, ErrMessage )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! DEPENDS ON: setup_umhdr
CALL Setup_UMHdr( NewUMHdr )
NewUMHdr % Len2Lookup = MaxFlds
NewUMHdr % NumFlds = 0
NewUMHdr % IntC    = UMHdr % IntC
NewUMHdr % RealC   = UMHdr % RealC
NewUMHdr % LevDepC = UMHdr % LevDepC
NewUMHdr % RowDepC = UMHdr % RowDepC
NewUMHdr % ColDepC = UMHdr % ColDepC
ALLOCATE( NewUMHdr % Lookup(NewUMHdr % Len2Lookup) )

9999 CONTINUE

END SUBROUTINE New_UMHdr

#endif
