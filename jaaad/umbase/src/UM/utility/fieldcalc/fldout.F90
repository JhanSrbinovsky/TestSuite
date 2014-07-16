#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to write fields to a UM fieldsfile


!=======================================================================

SUBROUTINE FldOut( Field,       &  ! inout
                   UMHdr,       &  ! inout
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
! 6.0     23/01/04 I/O optimisation. Dave Robinson.
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  PP_Field_type,          &
  PP_Header_type,         &
  UM_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(INOUT) :: Field
TYPE(UM_Header_type), INTENT(INOUT) :: UMHdr
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "FldOut"
#include "cntl_io.h"

! Local Variables:
INTEGER :: LookupPos
INTEGER :: DataAddress
INTEGER :: DataLen
INTEGER :: Len_IO
REAL :: Err_IO
REAL, Allocatable :: Field_Out (:,:)     !  Field written to disk

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Find and set next point to write data
IF( UMHdr % NumFlds == 0) THEN
  DataAddress = UMHdr % StartData - 1  ! First field
ELSE
  ! Start of prev fld + len of prev fld
  DataAddress = UMHdr % Lookup(UMHdr % NumFlds) % DataPos +  &
                UMHdr % Lookup(UMHdr % NumFlds) % LBNREC
END IF

! DEPENDS ON: setpos
CALL SETPOS ( UMHdr % UnitNum, DataAddress, ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus, "Failure in SETPOS" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! Set New Position in lookup - add to end of lookup
UMHdr % NumFlds = UMHdr % NumFlds + 1
LookupPos = UMHdr % NumFlds
DataLen = Field % Hdr % LBLREC

! update lookup for field - assume LBLREC is correct, but LBNREC isn't
! LBLREC is incorrect - it should contain the full grid size
Field % Hdr % DataPos = DataAddress
Field % Hdr % LBUser2 = DataAddress
Field % Hdr % LBNREC  = UM_Sector_Size *              &  ! round up
                      ( (DataLen + UM_Sector_Size-1) / UM_Sector_Size )
UMHdr % Lookup( LookupPos ) = Field % Hdr

! Field % RData contains the exact length of data (DataLen)
! Field_Out is rounded up to the next sector boundary

Allocate (Field_Out ( Field % Hdr % LBNREC, 1 ) )

! Copy data and set padding to next sector boundary to zero

Field_Out (1:DataLen, 1) = Field % RData (1:DataLen, 1)
Field_Out (DataLen+1:Field % Hdr % LBNREC, 1) = 0.0

! DEPENDS ON: buffout
CALL BUFFOUT( UMHdr % UnitNum, Field_Out,               &
              Field % Hdr % LBNREC, Len_IO, Err_IO )

! If checking the error code returned by BUFFOUT remember Err_IO == -1.0
! means write successful

Deallocate (Field_Out)

9999 CONTINUE

END SUBROUTINE FldOut

#endif
