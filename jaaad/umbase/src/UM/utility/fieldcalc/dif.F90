#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generic routines for manipulating pp-fields within Fieldcalc


!=======================================================================
SUBROUTINE Dif( Field1,       &  ! in
                Field2,       &  ! in
                DifField,     &  ! inout
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

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: Field1
TYPE(PP_Field_type), INTENT(IN) :: Field2

TYPE(PP_Field_type), INTENT(INOUT) :: DifField ! Location of result
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Dif"
#include "c_mdi.h"

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( (Field1 % Hdr % NumCols /= Field2 % Hdr % NumCols) .OR. &
     (Field1 % Hdr % NumRows /= Field2 % Hdr % NumRows) ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning, &
                "Cannot subtract fields of different dimensions" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED( DifField % RData ) ) THEN
  DEALLOCATE( DifField % RData )
END IF
DifField % Hdr = Field1 % Hdr
ALLOCATE( DifField % RData(DifField % Hdr % NumCols, &
                           DifField % Hdr % NumRows) )
DifField % Hdr % BMDI = RMDI
DifField % Hdr % STCode = IMDI

WHERE ( (Field1 % RData /= Field1 % Hdr % BMDI) .AND.  &
        (Field2 % RData /= Field2 % Hdr % BMDI) )
  DifField % RData = Field1 % RData - Field2 % RData
ELSEWHERE
  DifField % RData = RMDI
END WHERE

9999 CONTINUE

END SUBROUTINE Dif

!=======================================================================

!=======================================================================

!=======================================================================

!=======================================================================
! Vector direction - NOT wind direction

#endif
