#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to convert C-grid fields to B-grid U/V at UM5

SUBROUTINE CtoBgrid( CField,       &  ! in
                     B_Hdr,        &  ! in
                     BField,       &  ! inout
                     ErrorStatus )    ! inout

! Description:
!   Routine to interpolate ANY C-grid field to the B-grid UV points.
!
! Method:
!   Do x-direction, then y-direction.  Compare Zeroth Pt for C-grid and
!   B-grid.
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

! Description:  Routine to interpolate ANY C-grid field to the B-grid UV
!               points
! Method:       Do x-direction, then y-direction.  Compare Zeroth Pt for
!               C-grid and B-grid.

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: CField  ! C-grid field to interpolate
TYPE(PP_Header_type), INTENT(IN) :: B_Hdr  ! Template for B-grid header

TYPE(PP_Field_type), INTENT(INOUT) :: BField ! Resulting B-grid field
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "CtoBgrid"
#include "c_mdi.h"

! Local Variables:
INTEGER :: NumCols, NumRows
REAL :: XDiff, YDiff
REAL, ALLOCATABLE :: TempF(:,:)

! End of Header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! Check for similar grid spacing
IF ( (ABS( B_Hdr % LonInt - CField % Hdr % LonInt ) > 0.0001) .OR. &
     (ABS( B_Hdr % LatInt - CField % Hdr % LatInt ) > 0.0001) ) THEN
  WRITE(6,*) "CtoBgrid: Grid spacing incompatible"
  ErrorStatus = StatusWarning
  GO TO 9999
END IF
! Check how fields line up
XDiff = ( B_Hdr % ZerothLon - CField % Hdr % ZerothLon ) / B_Hdr%LonInt
YDiff = ( B_Hdr % ZerothLat - CField % Hdr % ZerothLat ) / B_Hdr%LatInt
IF ( (ABS(XDiff) > 1.0) .OR. (ABS(YDiff) > 1.0) ) THEN
  WRITE(6,*) "CtoBgrid: Grids more than 1 row / col out of alignment"
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

IF ( ASSOCIATED( BField % RData ) ) THEN
  DEALLOCATE( BField % RData )
END IF
BField % Hdr = CField % Hdr                ! Copy date and variable info
BField % Hdr % NumRows   = B_Hdr % NumRows ! Modify grid definition
BField % Hdr % NumCols   = B_Hdr % NumCols
BField % Hdr % LBLREC    = B_Hdr % NumCols * B_Hdr % NumRows
BField % Hdr % ZerothLat = B_Hdr % ZerothLat
BField % Hdr % ZerothLon = B_Hdr % ZerothLon
BField % Hdr % BMDI      = RMDI
ALLOCATE( BField % RData(BField % Hdr % NumCols, &
                         BField % Hdr % NumRows) )
! Insert an artificial STASH code to avoid redoing these interpolations
BField % Hdr % STCode = - CField % Hdr % STCode
BField % RData(:,:) = RMDI

NumCols = CField % Hdr % NumCols
NumRows = CField % Hdr % NumRows
ALLOCATE( TempF(NumCols, NumRows) )
TempF(:,:) = CField % RData(:,:)

IF ( XDiff > 0.1 ) THEN
  ! Do x-direction interpolations
  TempF(1:NumCols-1,1:NumRows) = 0.5 *                                 &
                               ( CField % RData(1:NumCols-1,1:NumRows) &
                               + CField % RData(2:NumCols  ,1:NumRows) )
  ! Check for wraparound - do edges of field
  IF ( (CField % Hdr % LBHEM == 0) .OR.  &
       (CField % Hdr % LBHEM == 4) ) THEN
    TempF(NumCols,1:NumRows)   = 0.5 *                                 &
                               ( CField % RData(NumCols,1:NumRows)     &
                               + CField % RData(1      ,1:NumRows) )
  ELSE
    TempF(NumCols,1:NumRows) = RMDI
  END IF

ELSE IF ( XDiff < -0.1 ) THEN
  TempF(2:NumCols,1:NumRows)   = 0.5 *                                 &
                               ( CField % RData(1:NumCols-1,1:NumRows) &
                               + CField % RData(2:NumCols  ,1:NumRows) )
  ! Check for wraparound - do edges of field
  IF ( (CField % Hdr % LBHEM == 0) .OR.  &
       (CField % Hdr % LBHEM == 4) ) THEN
    TempF(1,1:NumRows) = 0.5 * ( CField % RData(NumCols,1:NumRows)     &
                               + CField % RData(1      ,1:NumRows) )
  ELSE
    TempF(1,1:NumRows) = RMDI
  END IF

END IF

IF ( YDiff > 0.1 ) THEN
  ! Do x-direction interpolations
  TempF(1:NumCols,1:NumRows-1) = 0.5 * ( TempF(1:NumCols,1:NumRows-1)  &
                                       + TempF(1:NumCols,2:NumRows  ) )
  TempF(1:NumCols,NumRows) = RMDI

ELSE IF ( YDiff < -0.1 ) THEN
  TempF(1:NumCols,2:NumRows)   = 0.5 * ( TempF(1:NumCols,1:NumRows-1)  &
                                       + TempF(1:NumCols,2:NumRows  ) )
  TempF(1:NumCols,1) = RMDI

END IF

! Put TempF into correct position in BField % RData
NumCols = MIN( CField % Hdr % NumCols, BField % Hdr % NumCols )
NumRows = MIN( CField % Hdr % NumRows, BField % Hdr % NumRows )
BField % RData(1:NumCols,1:NumRows) = TempF(1:NumCols,1:NumRows)
DEALLOCATE( TempF )

9999 CONTINUE

END SUBROUTINE CtoBgrid

#endif
