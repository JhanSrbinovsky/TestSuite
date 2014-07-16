#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate tropopause height

SUBROUTINE TropHeight( NumLevs,      &  ! in
                       PFields,      &  ! in
                       TFields,      &  ! in
                       ZFields,      &  ! in
                       TropP,        &  ! inout
                       TropT,        &  ! inout
                       TropZ,        &  ! inout
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
! 6.2     19/12/05 Optimisation mods. J-C Rioual (NEC)
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK
USE FldCodes_Mod, ONLY:         &
  ST_TropP,  MO8_TropP,  PP_P,  &
  ST_TropT,  MO8_TropT,  PP_T,  &
  ST_TropZ,  MO8_TropZ,  PP_Z,  &
  VC_Trop,   LV_Special
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs) ! Pressure, temp
TYPE(PP_Field_type), INTENT(IN) :: TFields(NumLevs) !   and height
TYPE(PP_Field_type), INTENT(IN) :: ZFields(NumLevs) !   on theta levels

TYPE(PP_Field_type), INTENT(INOUT) :: TropP         ! Pressure, temp
TYPE(PP_Field_type), INTENT(INOUT) :: TropT         !   and height
TYPE(PP_Field_type), INTENT(INOUT) :: TropZ         !   at Tropopause
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "TropHeight"
#include "c_g.h"
#include "c_r_cp.h"
#include "c_mdi.h"
#include "c_lapse.h"
REAL, PARAMETER :: G_over_R = G / R
REAL, PARAMETER :: VSmall = 1.e-6

! Local Variables:
INTEGER :: i,j,k           ! Loop counters
INTEGER :: TLev(PFields(1)%Hdr%NumCols,PFields(1)%Hdr%NumRows)
REAL :: PLwr, PUpr
REAL :: TLwr, TUpr
REAL :: ZLwr, ZUpr
REAL :: LapseUpr, LapseLwr
! Lapse rates below, at and above current layer
REAL :: Lapse_b(PFields(1)%Hdr%NumCols,PFields(1)%Hdr%NumRows)
REAL :: Lapse_ (PFields(1)%Hdr%NumCols,PFields(1)%Hdr%NumRows)
REAL :: Lapse_a(PFields(1)%Hdr%NumCols,PFields(1)%Hdr%NumRows)
REAL :: delta_lapse        ! Lapse_b-Lapse_a

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( TropP % RData ) ) THEN
  DEALLOCATE( TropP % RData )
END IF
IF ( ASSOCIATED( TropT % RData ) ) THEN
  DEALLOCATE( TropT % RData )
END IF
IF ( ASSOCIATED( TropZ % RData ) ) THEN
  DEALLOCATE( TropZ % RData )
END IF
TropP % Hdr = PFields(1) % Hdr
TropP % Hdr % LBVC     = VC_Trop
TropP % Hdr % MO8Level = LV_Special
TropP % Hdr % BULEV    = 0.0
TropP % Hdr % BHULEV   = 0.0
TropP % Hdr % RLevel   = 0.0
TropP % Hdr % RefLevel = 0.0
TropP % Hdr % BHLEV    = 0.0
TropP % Hdr % BHRLEV   = 0.0
TropP % Hdr % BMDI     = RMDI
TropT % Hdr = TropP % Hdr
TropZ % Hdr = TropP % Hdr
TropP % Hdr % PPCode   =  PP_P
TropP % Hdr % MO8Type  = MO8_TropP
TropP % Hdr % STCode   =  ST_TropP
TropT % Hdr % PPCode   =  PP_T
TropT % Hdr % MO8Type  = MO8_TropT
TropT % Hdr % STCode   =  ST_TropT
TropZ % Hdr % PPCode   =  PP_Z
TropZ % Hdr % MO8Type  = MO8_TropZ
TropZ % Hdr % STCode   =  ST_TropZ
ALLOCATE( TropP % RData(TropP % Hdr % NumCols, &
                        TropP % Hdr % NumRows) )
ALLOCATE( TropT % RData(TropT % Hdr % NumCols, &
                        TropT % Hdr % NumRows) )
ALLOCATE( TropZ % RData(TropZ % Hdr % NumCols, &
                        TropZ % Hdr % NumRows) )

TLev (:,:) = 0
Lapse_  = (TFields(1) % RData - TFields(2) % RData) / &
                 (ZFields(2) % RData - ZFields(1) % RData)
Lapse_a = (TFields(2) % RData - TFields(3) % RData) / &
                 (ZFields(3) % RData - ZFields(2) % RData)

DO k = 3,NumLevs
  Lapse_b = Lapse_
  Lapse_  = Lapse_a
  IF ( k < NumLevs ) THEN
    Lapse_a = (TFields(k) % RData - TFields(k+1) % RData) / &
                     (ZFields(k+1) % RData - ZFields(k) % RData)
  END IF

  WHERE( (PFields(k  ) % RData < 50000.0) .AND.  &
         (PFields(k-1) % RData > 5000.0 ) .AND.  &
         (Lapse_  < lapse_trop)           .AND.  &
         (Lapse_a < lapse_trop)           .AND.  &
         (Lapse_b > 0.0)                  .AND.  &
         (TLev == 0) )
    TLev = k
  END WHERE
END DO

! Optimisation introduced at vn6.2 results in differences in Global 
! as well as NAE and is reversed pending further investigation
k = NumLevs - 2 

! Optimisation added at vn6.2
! DO j = 1,TropP % Hdr % NumRows
!   DO i = 1,TropP % Hdr % NumCols

DO i = 1,TropP % Hdr % NumCols
  DO j = 1,TropP % Hdr % NumRows

! Optimisation added at vn6.2 - scalar variable must be initialised
! inside the loop to allow vectorisation
!    k = NumLevs - 2

    IF ( TLev(i,j) /= 0 ) THEN   ! This is a TEMPORARY fix to avoid
      k = TLev(i,j)              !  falling over when trop is not
    END IF                       !  found.  Will have to rewrite this.

    PLwr = PFields(k-1) % RData(i,j)
    PUpr = PFields(k  ) % RData(i,j)
    TLwr = TFields(k-1) % RData(i,j)
    TUpr = TFields(k  ) % RData(i,j)
    ZLwr = ZFields(k-1) % RData(i,j)
    ZUpr = ZFields(k  ) % RData(i,j)
    LapseUpr = (TUpr - TFields(k+1) % RData(i,j)) / &
                      (ZFields(k+1) % RData(i,j) - ZUpr)
    LapseLwr = (TFields(k-2) % RData(i,j) - TLwr) / &
                      (ZLwr - ZFields(k-2) % RData(i,j))

    ! Find Z at tropopause
    delta_lapse = LapseLwr - LapseUpr
    IF ( ABS(delta_lapse) < VSmall ) THEN
      IF ( delta_lapse >= 0 ) delta_lapse =  VSmall
      IF ( delta_lapse <  0 ) delta_lapse = -VSmall
    END IF
    TropZ % RData(i,j) = ( (TLwr+(LapseLwr*ZLwr)) -   &
                                (TUpr+(LapseUpr*ZUpr)) ) / delta_lapse

    IF ( TropZ % RData(i,j) < ZLwr ) THEN
      TropZ % RData(i,j) = ZLwr  ! ensure trop level doesn't undershoot
    END IF
    IF ( TropZ % RData(i,j) > ZUpr )  THEN
      TropZ % RData(i,j) = ZUpr  ! or overshoot
    END IF

    ! T at tropopause
    TropT % RData(i,j) = TLwr - LapseLwr * (TropZ % RData(i,j) - ZLwr)

    IF ( ABS(LapseLwr) < VSmall ) THEN
      IF ( LapseLwr >= 0 ) LapseLwr =  VSmall
      IF ( LapseLwr <  0 ) LapseLwr = -VSmall
    END IF

    ! P at tropopause
    TropP % RData(i,j) = PLwr *  &
         (TropT % RData(i,j)/TLwr)**(G_over_R/LapseLwr)

  END DO
END DO

WHERE (TropP % RData < 10100.0)      ! set arbitrary max tropopause
  TropP % RData = 10100.0
  TropT % RData = 199.0
  TropZ % RData = 16180.0
END WHERE

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE TropHeight

#endif
