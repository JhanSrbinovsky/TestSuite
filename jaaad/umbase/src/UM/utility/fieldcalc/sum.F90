#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generic routines for manipulating pp-fields within Fieldcalc

SUBROUTINE Sum( Field1,       &  ! in
                Field2,       &  ! in
                SumField,     &  ! inout
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
TYPE(PP_Field_type), INTENT(IN) :: Field1       ! Fields to be added
TYPE(PP_Field_type), INTENT(IN) :: Field2

TYPE(PP_Field_type), INTENT(INOUT) :: SumField  ! Sum of two fields
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Sum"
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
                "Cannot add fields of different dimensions" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

! If there is data in SumField, and it is neither of the two input
! fields, get rid of it
IF ( ASSOCIATED( SumField % RData )            .AND.  &
    (SumField % ArrayPos /= Field1 % ArrayPos) .AND.  &
    (SumField % ArrayPos /= Field2 % ArrayPos) ) THEN
  DEALLOCATE( SumField % RData )
  NULLIFY( SumField % RData )
END IF
! IF SumField is now empty, allocate some space
IF ( .NOT. ASSOCIATED( SumField % RData ) ) THEN
  ALLOCATE( SumField % RData(Field1 % Hdr % NumCols, &
                             Field1 % Hdr % NumRows) )
  SumField % Hdr = Field1 % Hdr
END IF
SumField % Hdr % BMDI = RMDI
SumField % Hdr % STCode = IMDI

WHERE ( (Field1 % RData /= Field1 % Hdr % BMDI) .AND.  &
        (Field2 % RData /= Field2 % Hdr % BMDI) )
  SumField % RData = Field1 % RData + Field2 % RData
ELSEWHERE
  SumField % RData = RMDI
END WHERE

9999 CONTINUE

END SUBROUTINE Sum

!=======================================================================

!=======================================================================

!=======================================================================

!=======================================================================

!=======================================================================
! Vector direction - NOT wind direction

#endif
