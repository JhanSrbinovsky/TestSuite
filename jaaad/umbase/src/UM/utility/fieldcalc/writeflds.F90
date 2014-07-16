#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to write fields to a UM fieldsfile

SUBROUTINE WriteFlds ( NumFlds,     &  ! in
                       MaxFlds,     &  ! in
                       MaxFldsOut,  &  ! in
                       PackType,    &  ! in
                       Source,      &  ! in
                       PackAcc,     &  ! in
                       Fields,      &  ! inout
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
  StatusWarning,          &
  StatusFatal
USE FldCodes_mod, ONLY:   &
  ProcessNo
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds
INTEGER, INTENT(IN) :: MaxFlds
INTEGER, INTENT(IN) :: MaxFldsOut
CHARACTER(LEN=6), INTENT(IN) :: PackType
INTEGER, INTENT(IN) :: Source(NumFlds)
REAL,    INTENT(IN) :: PackAcc(NumFlds)

TYPE(PP_Field_type), INTENT(INOUT) :: Fields(MaxFlds)
TYPE(UM_Header_type), INTENT(INOUT) :: UMHdr
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "WriteFlds"
#include "cntl_io.h"

! Local Variables:
INTEGER :: i
INTEGER :: ifld
LOGICAL :: compressed = .FALSE.
TYPE(PP_Field_type) :: TempField

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

NULLIFY( TempField % RData )

! Check there is space in the lookup table
IF( NumFlds + UMHdr % NumFlds > UMHdr % Len2Lookup ) THEN
  WRITE(6,*) "Insufficient space allocated for Lookup Table"
  WRITE(6,*) "Increase the value of MaxFldsOut in main prog to correct"
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusFatal, &
               "Insufficient space allocated for Lookup Table." )
END IF

DO i = 1,NumFlds
  ifld = Source(i)
  compressed = .FALSE.
  Fields(ifld) % Hdr % LBUser7 = ProcessNo    ! Superstash code

  ! Packing if required
  IF ( PackType == "NONE" ) THEN
    TempField = Fields(ifld)
  ELSE
! DEPENDS ON: pack_single
    CALL Pack_Single( PackAcc(i), PackType, Fields(ifld), TempField, &
                      ErrorStatus )
    IF ( ErrorStatus == StatusOK ) THEN
      compressed = .TRUE.
    ELSE
      WRITE(6,*) "Error in Pack_Single - writing unpacked field."
      TempField = Fields(ifld)
    END IF
  END IF

! DEPENDS ON: fldout
  CALL FldOut( TempField, UMHdr, ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
! DEPENDS ON: ereport
    CALL EReport( RoutineName, StatusWarning, "Error writing field." )
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF

  IF ( compressed ) THEN
    DEALLOCATE( TempField % RData )  ! deallocate compressed field
    NULLIFY( TempField % RData )
  END IF

END DO

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE WriteFlds

!=======================================================================


#endif
