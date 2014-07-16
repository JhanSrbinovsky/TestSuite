#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate vorticity diagnostic

SUBROUTINE Vortic( UField,       &  ! in
                   VField,       &  ! in
                   VorField,     &  ! inout
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
  ST_Vortic, MO8_Vortic, PP_Vortic
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: UField    ! U Component of wind
TYPE(PP_Field_type), INTENT(IN) :: VField    ! V Component of wind

TYPE(PP_Field_type), INTENT(INOUT) :: VorField  ! Vorticity field
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Vortic"
#include "c_mdi.h"

! Local Variables:
TYPE(PP_Field_type) :: DuDy      ! latitudinal derivative of u
TYPE(PP_Field_type) :: DvDx      ! longitudinal derivative of v

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

Nullify ( DuDy % RData )
Nullify ( DvDx % RData )

IF ( ( UField % Hdr % NumCols /= VField % Hdr % NumCols )  .OR.   &
     ( UField % Hdr % NumRows /= VField % Hdr % NumRows ) )     THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning, &
               "U&V Fields must have the same dimensions" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED( VorField % RData ) ) THEN
  DEALLOCATE( VorField % RData )
END IF
VorField % Hdr = UField % Hdr
VorField % Hdr % PPcode  =  PP_Vortic
VorField % Hdr % MO8Type = MO8_Vortic
VorField % Hdr % STCode  =  ST_Vortic
VorField % Hdr % BMDI    = RMDI
ALLOCATE( VorField % RData(VorField % Hdr % NumCols, &
                           VorField % Hdr % NumRows) )

! DEPENDS ON: diffcosy
CALL DiffCOSY( UField, DuDy, ErrorStatus )  ! lon derivative
! DEPENDS ON: diffx
CALL DiffX   ( VField, DvDx, ErrorStatus )  ! lat derivative
IF ( ErrorStatus /= StatusOK ) THEN
  GO TO 9999
END IF

WHERE( (DvDx % RData /= DvDx % Hdr % BMDI) .AND. &
       (DuDy % RData /= DuDy % Hdr % BMDI) )
  VorField % RData = DvDx % RData - DuDy % RData
ELSEWHERE
  VorField % RData = RMDI
END WHERE

DEALLOCATE( DuDy % RData )
DEALLOCATE( DvDx % RData )

9999 CONTINUE

END SUBROUTINE Vortic

#endif
