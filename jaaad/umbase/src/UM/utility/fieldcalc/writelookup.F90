#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to write lookup table to file

SUBROUTINE WriteLookup( UMHdr,        & ! in
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
  UM_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
TYPE(UM_Header_type), INTENT(IN) :: UMHdr

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "WriteLookup"

! Local variables:
INTEGER :: Len_IO       ! Length of buffer written by BUFFOUT
INTEGER :: Residue      ! Space for lookup not required

REAL    :: Err_IO       ! BUFFOUT error status

! End of header --------------------------------------------------------

! DEPENDS ON: setpos
CALL SETPOS( UMHdr % UnitNum, UMHdr % FixHd(150)-1, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning, "Failure in SETPOS" )
END IF

! DEPENDS ON: buffout
CALL BUFFOUT( UMHdr % UnitNum, UMHdr % Lookup,     &
              UMHdr%Len1Lookup * UMHdr%NumFlds,    &
              Len_IO, Err_IO )
IF ( Err_IO /= -1.0 ) THEN      ! BUFFOUT returns -1.0 when OK
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning, "Failure in BUFFOUT" )
END IF

! Should flush buffer here?

9999 CONTINUE

END SUBROUTINE WriteLookup

#endif
