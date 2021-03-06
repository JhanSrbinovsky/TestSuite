#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate the height of a temperature surface

SUBROUTINE Isotherm( NumLevs,      &  ! in
                     TempRef,      &  ! in
                     Orog,         &  ! in
                     PStar,        &  ! in
                     PFields,      &  ! in
                     TFields,      &  ! in
                     ZFields,      &  ! in
                     IsoThermZ,    &  ! inout
                     IsoThermP,    &  ! inout
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

USE IO_mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK
USE FldCodes_Mod, ONLY:   &
  ST_FreezZ, MO8_FreezZ,  &
  ST_FreezP, MO8_FreezP,  &
  ST_Iso20Z, MO8_Iso20Z,  &
  ST_Iso20P, MO8_Iso20P,  &
  PP_Z,      PP_P,        &
  VC_Freez,  VC_Iso20,    &
  LV_Special
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
REAL, INTENT(IN) :: TempRef
TYPE(PP_Field_type), INTENT(IN) :: Orog             ! Model Orography
TYPE(PP_Field_type), INTENT(IN) :: PStar            ! Surface Pressure
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs) ! Pressure, temp
TYPE(PP_Field_type), INTENT(IN) :: TFields(NumLevs) !  and height
TYPE(PP_Field_type), INTENT(IN) :: ZFields(NumLevs) !  on theta levels

TYPE(PP_Field_type), INTENT(INOUT) :: IsoThermZ   ! Height and pressure
TYPE(PP_Field_type), INTENT(INOUT) :: IsoThermP   !  of Isotherm
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Isotherm"
#include "c_g.h"
#include "c_r_cp.h"
#include "c_mdi.h"
INTEGER, PARAMETER :: IZeroC = 273
REAL, PARAMETER :: G_over_R = G / R

! Local Variables:
INTEGER :: i, j, k
INTEGER :: IsoLev(PStar%Hdr%NumCols,PStar%Hdr%NumRows)
REAL :: lapse_rate        ! lapse rate of layer containing TempRef

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( IsoThermZ % RData ) ) THEN
  DEALLOCATE( IsoThermZ % RData )
END IF
IF ( ASSOCIATED( IsoThermP % RData ) ) THEN
  DEALLOCATE( IsoThermP % RData )
END IF
IsoThermZ % Hdr = PStar % Hdr
IsoThermZ % Hdr % PPCode = PP_Z
IsoThermZ % Hdr % MO8Level = LV_Special
IsoThermZ % Hdr % STCode = IMDI
IsoThermZ % Hdr % BMDI = RMDI
IsoThermP % Hdr = IsoThermZ % Hdr
IsoThermP % Hdr % PPCode = PP_P
IF      ( INT( TempRef ) == IZeroC ) THEN     ! Freezing Level (to 0 dp)
  IsoThermZ % Hdr % LBVC    =  VC_Freez
  IsoThermZ % Hdr % MO8Type = MO8_FreezZ
  IsoThermZ % Hdr % STCode  =  ST_FreezZ
  IsoThermP % Hdr % LBVC    =  VC_Freez
  IsoThermP % Hdr % MO8Type = MO8_FreezP
  IsoThermP % Hdr % STcode  =  ST_FreezP
ELSE IF ( INT( TempRef ) == (IZeroC-20) ) THEN ! -20C Isotherm (to 0 dp)
  IsoThermZ % Hdr % LBVC    =  VC_Iso20
  IsoThermZ % Hdr % MO8Type = MO8_Iso20Z
  IsoThermZ % Hdr % STCode  =  ST_Iso20Z
  IsoThermP % Hdr % LBVC    =  VC_Iso20
  IsoThermP % Hdr % MO8Type = MO8_Iso20P
  IsoThermP % Hdr % STcode  =  ST_Iso20P
ELSE
  WRITE(6,*) "Warning: Isotherm value not recognised - STCode not set"
END IF
ALLOCATE( IsoThermZ % RData(IsoThermZ % Hdr % NumCols, &
                            IsoThermZ % Hdr % NumRows) )
ALLOCATE( IsoThermP % RData(IsoThermP % Hdr % NumCols, &
                            IsoThermP % Hdr % NumRows) )

IsoLev(:,:) = 0
WHERE( TFields(1) % RData(:,:) < TempRef )
  IsoLev(:,:) = 1
  IsoThermZ % RData(:,:) = Orog  % RData(:,:)
  IsoThermP % RData(:,:) = Pstar % RData(:,:)
END WHERE
DO k = 2,NumLevs
  WHERE( (TFields(k  ) % RData(:,:) < TempRef) .AND. &
         (TFields(k-1) % RData(:,:) > TempRef) )
    IsoLev(:,:) = k
  END WHERE
END DO

DO i = 1,IsoThermZ % Hdr % NumCols
  DO j = 1,IsoThermZ % Hdr % NumRows
    k = IsoLev(i,j)
    IF ( k >= 2 ) THEN

      lapse_rate = ( TFields(k-1)%RData(i,j) - TFields(k)%RData(i,j) ) &
                  / ( ZFields(k)%RData(i,j) - ZFields(k-1)%RData(i,j) )
      IsoThermZ % RData(i,j) = ZFields(k-1) % RData(i,j) +             &
                      (TFields(k-1) % RData(i,j)-TempRef)/lapse_rate
      IsoThermP % RData(i,j) = PFields(k-1) % RData(i,j) *             &
              (TempRef/TFields(k-1) % RData(i,j))**(G_over_R/lapse_rate)
    ENDIF
  END DO
END DO

WHERE( IsoLev(:,:) == 0 )
  IsoThermZ % RData(:,:) = -1
  IsoThermP % RData(:,:) = -1
END WHERE

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE Isotherm

#endif
