#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate divergence diagnostic

SUBROUTINE Diverg( UField,       &  ! in
                   VField,       &  ! in
                   DivField,     &  ! inout
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
USE FldCodes_Mod, ONLY:   &
  ST_Diverg, MO8_Diverg, PP_Diverg
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: UField       ! U Component of wind
TYPE(PP_Field_type), INTENT(IN) :: VField       ! V Component of wind

TYPE(PP_Field_type), INTENT(INOUT) :: DivField  ! Divergence field
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Diverg"
#include "c_mdi.h"

! Local Variables:
TYPE(PP_Field_type) :: DuDx      ! longitudinal derivative of u
TYPE(PP_Field_type) :: DvDy      ! latitudinal derivative of v

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

Nullify ( DuDx % RData )
Nullify ( DvDy % RData )

IF ( ( UField % Hdr % NumCols /= VField % Hdr % NumCols )  .OR.   &
     ( UField % Hdr % NumRows /= VField % Hdr % NumRows ) )     THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning, &
               "U&V Fields must have the same dimensions" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED( DivField % RData ) ) THEN
  DEALLOCATE( DivField % RData )
END IF
DivField % Hdr = UField % Hdr
DivField % Hdr % PPCode  =  PP_Diverg
DivField % Hdr % MO8Type = MO8_Diverg
DivField % Hdr % STCode  =  ST_Diverg
DivField % Hdr % BMDI    = RMDI
ALLOCATE( DivField % RData(DivField % Hdr % NumCols, &
                           DivField % Hdr % NumRows) )

! DEPENDS ON: diffx
CALL DiffX   ( UField, DuDx, ErrorStatus )  ! lat derivative
! DEPENDS ON: diffcosy
CALL DiffCOSY( VField, DvDy, ErrorStatus )  ! lon derivative
IF ( ErrorStatus /= StatusOK ) THEN
  GO TO 9999
END IF

WHERE( (DuDx % RData /= DuDx % Hdr % BMDI) .AND. &
       (DvDy % RData /= DvDy % Hdr % BMDI) )
  DivField % RData = DuDx % RData + DvDy % RData
ELSEWHERE
  DivField % RData = RMDI
END WHERE

DEALLOCATE( DuDx % RData )
DEALLOCATE( DvDy % RData )

9999 CONTINUE

END SUBROUTINE Diverg

#endif
