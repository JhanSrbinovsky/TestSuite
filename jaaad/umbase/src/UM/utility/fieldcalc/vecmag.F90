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
SUBROUTINE VecMag( Field1,      &  ! in
                   Field2,      &  ! in
                   VMField,     &  ! inout
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
! 6.2     22/06/06 Changed 2.0 to 2 - R. Sempers (frpz)
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

TYPE(PP_Field_type), INTENT(INOUT) :: VMField
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "VecMag"
#include "c_mdi.h"

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ( Field1 % Hdr % NumCols /= Field2 % Hdr % NumCols )  .OR. &
     ( Field1 % Hdr % NumRows /= Field2 % Hdr % NumRows ) ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning, &
               "Vector fields must have the same dimensions" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED( VMField % RData ) ) THEN
  DEALLOCATE( VMField % RData )
END IF
VMField % Hdr = Field1 % Hdr
VMField % Hdr % STCode  = IMDI
VMField % Hdr % BMDI    = RMDI
ALLOCATE( VMField % RData(VMField % Hdr % NumCols, &
                          VMField % Hdr % NumRows) )

WHERE ( (Field1 % RData /= Field1 % Hdr % BMDI) .AND.  &
        (Field2 % RData /= Field2 % Hdr % BMDI) )
  VMField % RData = SQRT( (Field1 % RData)**2 + &
                          (Field2 % RData)**2)
ELSEWHERE
  VMField % RData = RMDI
END WHERE

9999 CONTINUE

END SUBROUTINE VecMag

!=======================================================================
! Vector direction - NOT wind direction

#endif
