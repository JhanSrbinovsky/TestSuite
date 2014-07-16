#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to calculate horizontal derivatives.



!=======================================================================


!=======================================================================

SUBROUTINE DiffCOSY( FField,       &  ! in
                     dFdY,         &  ! inout
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
TYPE(PP_Field_type), INTENT(IN) :: FField   ! Input variable, F

TYPE(PP_Field_type), INTENT(INOUT) :: dFdY  ! Lat derivative of F
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "DiffCOSY"
#include "c_mdi.h"
#include "c_a.h"
#include "c_pi.h"

! Local Variables:
INTEGER :: i,j
REAL :: HdelY      ! temporary scalar
REAL :: alat(FField % Hdr % NumRows)

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Must have regular, non-rotated lat/lon grid
IF ( FField % Hdr % LBCODE /= 1 ) THEN
  ErrorStatus = FField % Hdr % LBCODE
! DEPENDS ON: ereport
  CALL EReport( RoutineName, ErrorStatus,                       &
                "Cannot calculate latitudinal derivative - " // &
                "incorrect grid type" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED( dFdY % RData ) ) THEN
  DEALLOCATE( dFdY % RData )
END IF
dFdY % Hdr = FField % Hdr
ALLOCATE( dFdY % RData(dFdY % Hdr % NumCols, &
                       dFdY % Hdr % NumRows) )
dFdY % Hdr % BMDI = RMDI
dFdY % Hdr % STcode = IMDI

DO j = 1, FField % Hdr % NumRows
  alat(j) = FField % Hdr % ZerothLat + &
            FField % Hdr % LatInt * FLOAT(j)
ENDDO
! Note latitudinal centred derivatives can only be found for j=(2:N-1).
DO j = 2, FField % Hdr % NumRows-1
  HdelY = 0.5/( Earth_Radius * Pi_Over_180 * FField % Hdr % LatInt *   &
                          COS(Pi_Over_180*alat(j)) )
  dFdY % RData(:,j) = HdelY *                                          &
                 (FField % RData(:,j+1) * COS(Pi_Over_180*alat(j+1)) - &
                  FField % RData(:,j-1) * COS(Pi_Over_180*alat(j-1))  )
ENDDO
! Set N & S boundaries to missing data
dFdY % RData(:,1) = RMDI
dFdY % RData(:,FField % Hdr % NumRows) = RMDI

9999 CONTINUE

END SUBROUTINE DiffCOSY

#endif
