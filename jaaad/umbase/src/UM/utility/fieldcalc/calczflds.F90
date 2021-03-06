#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate geopotential height from pressure at UM5

SUBROUTINE CalcZFlds( NumLevs,      &  ! in
                      Orog,         &  ! in
                      PFields,      &  ! in
                      ZFields,      &  ! inout
                      ErrorStatus )    ! inout

! Description:
!   This subroutine uses orography and header information to calculate
!   geopotential height fields from pressure fields.
!
! Method:
!   First an appropriate STASH code is ascertained.  The horizontal
!   grid of the output fields is given by that of the input orography
!   field, and the vertical levels are defined by the header information
!   of the input pressure fields.  The geoptl height is calculated using
!   the level definition at UM5:  Z = Zsea(k) + (C(k) * Orog)
!   (See UMDP F3, Appendix 1)
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

USE IO_mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal
USE FldCodes_mod, ONLY:   &
  ST_Prho,  ST_Ptheta,    &
  ST_Zrho,  ST_Ztheta,    &
  MO8_Z,    PP_Z,         &
  ST_Orog
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs                         ! No. levels
TYPE(PP_Field_type), INTENT(IN) :: Orog                ! Model Orography
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs)    ! Pressure fields

TYPE(PP_Field_type), INTENT(INOUT) :: ZFields(NumLevs) ! Geoptl height
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "CalcZFlds"
#include "c_mdi.h"

! Local Variables:
INTEGER :: i                         ! Loop counter
INTEGER :: STCode_Z                  ! STASH code of output fields

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( .NOT. ASSOCIATED( Orog % RData ) ) THEN
! DEPENDS ON: ereport
  CALL EReport( RoutineName, StatusWarning, &
            "Orography not available - cannot calculate height fields" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

STCode_Z = 0
! PFields define vertical co-ordinates, Orog defines horizontal grid
IF      ( ABS(PFields(1) % Hdr % STCode) == ST_Prho   ) THEN
  STCode_Z = ST_Zrho              ! Height of rho levels (either grid)
ELSE IF ( ABS(PFields(1) % Hdr % STCode) == ST_Ptheta ) THEN
  STCode_Z = ST_Ztheta            ! Height of theta levels (either grid)
END IF
IF ( Orog % Hdr % STCode == -ST_Orog ) THEN
  STCode_Z = -STCode_Z             ! B(UV)-grid
END IF
IF ( STCode_Z == 0 ) THEN
  STCode_Z = IMDI
END IF

DO i = 1,NumLevs
  IF ( (ZFields(i) % Hdr % STCode == STCode_Z) .AND. &
       (STCode_Z /= IMDI)                      .AND. &
       (ASSOCIATED( ZFields(i) % RData )) ) THEN
    CYCLE     ! Already exists - move on to next field
  END IF

  IF ( ASSOCIATED( ZFields(i) % RData ) ) THEN
    DEALLOCATE( ZFields(i) % RData )
  END IF
  ZFields(i) % Hdr = PFields(i) % Hdr
  ZFields(i) % Hdr % PPCode   =     PP_Z
  ZFields(i) % Hdr % MO8Type  =    MO8_Z
  ZFields(i) % Hdr % STCode   = STCode_Z
  ZFields(i) % Hdr % BMDI     = RMDI

  ALLOCATE( ZFields(i) % RData(ZFields(i) % Hdr % NumCols, &
                               ZFields(i) % Hdr % NumRows) )
  ! Z = Zsea(k) + C(k) * orog
  ZFields(i) % RData = PFields(i) % Hdr % RLevel + &
                       PFields(i) % Hdr % BHLEV  * Orog % RData

END DO

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE CalcZFlds

#endif
