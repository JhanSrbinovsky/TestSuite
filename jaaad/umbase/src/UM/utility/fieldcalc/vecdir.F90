#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generic routines for manipulating pp-fields within Fieldcalc


!=======================================================================

!=======================================================================

!=======================================================================

!=======================================================================

!=======================================================================
! Vector direction - NOT wind direction
SUBROUTINE VecDir( Field1,      &  ! in
                   Field2,      &  ! in
                   VDField,     &  ! inout
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
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN)    :: Field1
TYPE(PP_Field_type), INTENT(IN)    :: Field2
TYPE(PP_Field_type), INTENT(INOUT) :: VDField

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "VecDir"
#include "c_mdi.h"
#include "c_pi.h"

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ( Field1 % Hdr % NumCols /= Field2 % Hdr % NumCols )  .OR.   &
     ( Field1 % Hdr % NumRows /= Field2 % Hdr % NumRows ) )     THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning,   &
               "VecDir : Fields must have the same dimensions" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED(VDField % RData) ) THEN
  DEALLOCATE( VDField % RData )
END IF
VDField % Hdr = Field1 % Hdr
ALLOCATE( VDField % RData(VDField % Hdr % NumCols, &
                          VDField % Hdr % NumRows) )
VDField % Hdr % BMDI = RMDI
VDField % Hdr % STCode = IMDI

WHERE( (Field1 % RData == 0.0) .AND. (Field2 % RData == 0.0) )
  VDField % RData = 0.0
ELSEWHERE ( (Field1 % RData /= Field1 % Hdr % BMDI) .AND.  &
            (Field2 % RData /= Field2 % Hdr % BMDI) )
  VDField % RData = ATAN2( Field1 % RData, Field2 % RData )*(180.0/Pi)
ELSEWHERE
  VDField % RData = RMDI
END WHERE

9999 CONTINUE

END SUBROUTINE VecDir

#endif
