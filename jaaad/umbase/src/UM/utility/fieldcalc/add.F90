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
SUBROUTINE Add( Factor,       &  ! in
                Field,        &  ! in
                AddField,     &  ! inout
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

TYPE(PP_Field_type), INTENT(INOUT) :: AddField ! Location of result
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Add"
#include "c_mdi.h"

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( AddField % RData )) THEN
  DEALLOCATE( AddField % RData )
END IF
AddField % Hdr = Field % Hdr
ALLOCATE( AddField % RData(AddField % Hdr % NumCols, &
                           AddField % Hdr % NumRows) )
AddField % Hdr % BMDI = RMDI
IF ( INT( Field % Hdr % STCode / 1000 ) /= UMSectnNo ) THEN
  AddField % Hdr % STCode = IMDI
ENDIF

WHERE ( Field % RData /= Field % Hdr % BMDI )
  AddField % RData = Factor + Field % RData
ELSEWHERE
  AddField % RData = RMDI
END WHERE

9999 CONTINUE

END SUBROUTINE Add

!=======================================================================

!=======================================================================

!=======================================================================
! Vector direction - NOT wind direction

#endif
