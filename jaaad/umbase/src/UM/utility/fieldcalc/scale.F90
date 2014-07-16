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
SUBROUTINE Scale( Factor,       &  ! in
                  Field,        &  ! in
                  ScField,      &  ! inout
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
  StatusOK
USE FldCodes_mod, ONLY:   &
  UMSectnNo
IMPLICIT None

! Subroutine Arguments:
REAL, INTENT(IN) :: Factor
TYPE(PP_Field_type), INTENT(IN) :: Field

TYPE(PP_Field_type), INTENT(INOUT) :: ScField ! Location of result
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Scale"
#include "c_mdi.h"

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( ScField % RData )) THEN
  DEALLOCATE( ScField % RData )
END IF
ScField % Hdr = Field % Hdr
ALLOCATE( ScField % RData(ScField % Hdr % NumCols, &
                          ScField % Hdr % NumRows) )
ScField % Hdr % BMDI = RMDI
IF ( INT( Field % Hdr % STCode / 1000 ) /= UMSectnNo ) THEN
  ScField % Hdr % STCode = IMDI
ENDIF

WHERE ( Field % RData /= Field % Hdr % BMDI )
  ScField % RData = Factor * Field % RData
ELSEWHERE
  ScField % RData = RMDI
END WHERE

9999 CONTINUE

END SUBROUTINE Scale

!=======================================================================

!=======================================================================
! Vector direction - NOT wind direction

#endif
