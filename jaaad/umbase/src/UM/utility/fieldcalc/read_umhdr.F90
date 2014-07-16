#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to read a UM Fieldsfile header

SUBROUTINE Read_UMHdr ( UMHdr,        &  ! inout
                        ErrorStatus )    ! inout

! Description:
!   This is a subroutine to open, read and store the contents of a UM
!   fieldsfile header.
!
! Method:
!   The environment variable containing the file name is held within the
!   derived type UM_Header_type.  This is used to open the file.  The
!   fixed length header is read, and UMHdr set up accordingly.  Space is
!   allocated for the lookup table, and the rest of the UM header read
!   in.  Memory is conserved by only allocating enough memory for the
!   exact number of lookup entries.
!
! Owner: Dave Robinson
!
! History:
! Version Date     Comment
! ------- ----     -------
! 1.0     02/05/03 Original Code.   Sara James
! 6.0     12/09/03 Code implemented into UM. Dave Robinson
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  LenFixHd,               &
  UM_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
TYPE(UM_Header_type), INTENT(INOUT) :: UMHdr
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Read_UMHdr"

! Local Variables:
INTEGER :: i
INTEGER :: WordAddress
INTEGER :: Len_IO
REAL :: Err_IO

INTEGER, POINTER :: Lookup(:,:)

CHARACTER(LEN=80) :: ErrMessage

! Dummy arguments for call to readhead.
INTEGER           :: ppxi
INTEGER           :: ppxrecs
CHARACTER         :: ppxc

! End of Header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
ENDIF

! DEPENDS ON: file_open
CALL File_Open ( UMHdr % UnitNum,                        &
                 UMHdr % FileNameEnv,                    &
                 LEN_TRIM(UMHdr % FileNameEnv),          &
                 0, 0, ErrorStatus)

IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, &
                "Failure to open input file" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

!-----------------------------------------------------------------------
! Allocate and read in Fixed_Length_Header
!-----------------------------------------------------------------------
ALLOCATE (UMHdr % FixHd(LenFixHd))

WordAddress = 0   ! Start of file
! DEPENDS ON: setpos
CALL SETPOS ( UMHdr % UnitNum, WordAddress, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! Is this just to check how much space should be allocated?
! DEPENDS ON: read_flh
CALL READ_FLH ( UMHdr % UnitNum, UMHdr % FixHd, LenFixHd,   &
                ErrorStatus, ErrMessage )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, ErrMessage )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! DEPENDS ON: setup_umhdr
CALL Setup_UMHdr( UMHdr )
! Need to change the data allocation
ALLOCATE ( Lookup (UMHdr % Len1Lookup, UMHdr % Len2Lookup) )

WordAddress = 0
! DEPENDS ON: setpos
CALL SETPOS ( UMHdr % UnitNum, WordAddress, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! DEPENDS ON: readhead
CALL READHEAD ( UMHdr % UnitNum, UMHdr % FixHd, LenFixHd,           &
       UMHdr % IntC,      UMhdr % LenIntC,                          &
       UMhdr % RealC,     UMhdr % LenRealC,                         &
       UMhdr % LevDepC,   UMhdr % Len1LevDepC, UMhdr % Len2LevDepC, &
       UMhdr % RowDepC,   UMhdr % Len1RowDepC, UMhdr % Len2RowDepC, &
       UMhdr % ColDepC,   UMhdr % Len1ColDepC, UMhdr % Len2ColDepC, &
       UMhdr % FldsOfC,   UMhdr % Len1FldsOfC, UMhdr % Len2FldsOfC, &
       UMhdr % ExtraC,    UMhdr % LenExtraC,                        &
       UMhdr % HistFile,  UMhdr % LenHistFile,                      &
       UMhdr % CompFldI1, UMhdr % LenCompFldI1,                     &
       UMhdr % CompFldI2, UMhdr % LenCompFldI2,                     &
       UMhdr % CompFldI3, UMhdr % LenCompFldI3,                     &
       Lookup, UMhdr % Len1Lookup, UMhdr % Len2Lookup,              &
       UMhdr % LenData,                                             &
#include "argppx.h"
       UMHdr % StartData,                                           &
       ErrorStatus, ErrMessage )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, ErrMessage )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

!-----------------------------------------------------------------------
! Deallocate blank lookup entries
!-----------------------------------------------------------------------
! Count Lookup entries
UMHdr % NumFlds = 0
DO i = 1,UMHdr % Len2Lookup
  IF ( Lookup(1, UMHdr%NumFlds+1) == -99 ) THEN
    UMHdr % Len2Lookup = UMHdr % NumFlds
    EXIT  ! Exit loop
  ELSE
    UMHdr % NumFlds = UMHdr % NumFlds+1
  END IF
END DO
DEALLOCATE( Lookup )

! DEPENDS ON: setpos
CALL SETPOS ( UMHdr % UnitNum, UMHdr % FixHd(150)-1, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

ALLOCATE( UMHdr % Lookup( UMHdr % Len2Lookup ) )
! DEPENDS ON: buffin
CALL BUFFIN( UMHdr % UnitNum, UMHdr % Lookup,         &
             UMHdr % Len1Lookup * UMHdr % Len2Lookup, &
             Len_IO, Err_IO )
IF ( Err_IO /= -1.0 ) THEN   ! Non-standard error code
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, "Failure in BUFFIN" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

9999 CONTINUE

END SUBROUTINE Read_UMHdr

#endif
