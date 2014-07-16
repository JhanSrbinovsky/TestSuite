#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate Max Wind value and height


!=======================================================================

SUBROUTINE TopBase( NumLevs,      &  ! in
                    UFields,      &  ! in
                    VFields,      &  ! in
                    PFields,      &  ! in
                    MaxWindU,     &  ! in
                    MaxWindV,     &  ! in
                    MaxWindP,     &  ! in
                    MaxWindBase,  &  ! inout
                    MaxWindTop,   &  ! inout
                    ErrorStatus )    ! inout

! Description:
!
! Method:
!   Old version uses U&V on rho levels, and pressure on theta levels to
!   interpolate between.  To save space we will interpolate between rho
!   levels.  This may mean we use a different number of levels, and we
!   may want to widen the 100-700 hPa range.
!
! Owner: Dave Robinson
!
! History:
! Version Date     Comment
! ------- ----     -------
! 6.1     09/06/04 Original Code. Dave Robinson
! 6.2     19/12/05 Optimisation changes and call MaxWindSplineV.
!                  J-C Rioual (NEC)
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
  ST_MaxWU,  MO8_MaxWU,  PP_U,  &
  ST_MaxWV,  MO8_MaxWV,  PP_V,  &
  ST_MaxWP,  MO8_MaxWP,  PP_P,  &
  VC_Upper,  VC_Lower,          &
  LV_Special,                   &
  ST_MWBase,   ST_MWTop,        &
  MO8_MxWBase, MO8_MxWTop

IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
TYPE(PP_Field_type), INTENT(IN) :: UFields(NumLevs) ! U-wind on B-grid
TYPE(PP_Field_type), INTENT(IN) :: VFields(NumLevs) ! V-wind on B-grid
! Pressure on B-grid U/V points on rho levels:
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs)

TYPE(PP_Field_type), INTENT(IN) :: MaxWindU  ! U comp of max wind
TYPE(PP_Field_type), INTENT(IN) :: MaxWindV  ! V comp of max wind
TYPE(PP_Field_type), INTENT(IN) :: MaxWindP  ! P at max wind level

TYPE(PP_Field_type), INTENT(INOUT) :: MaxWindBase ! P at max wind base
TYPE(PP_Field_type), INTENT(INOUT) :: MaxWindTop  ! P at max wind top

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "TopBase"

#include "c_mdi.h"
#include "cprintst.h"
INTEGER, PARAMETER :: kmax=5            ! total levels needed for spline
INTEGER, PARAMETER :: khalf=(kmax+1)/2  ! number of levels each side
INTEGER, PARAMETER :: ninc=16           ! number of increments used

! Local variables:
INTEGER :: i, j, k
INTEGER :: LevBtm, LevTop
INTEGER :: Level_Below, Level_Above
INTEGER :: Lev_Below_MaxW (PFields(1)%Hdr%NumCols,  &
                           PFields(1)%Hdr%NumRows)
INTEGER :: Lev_Above_MaxW (PFields(1)%Hdr%NumCols,  &
                           PFields(1)%Hdr%NumRows)
INTEGER :: MaxWinc(PFields(1)%Hdr%NumCols, PFields(1)%Hdr%NumRows)
REAL    :: MaxW   (PFields(1)%Hdr%NumCols, PFields(1)%Hdr%NumRows)
REAL :: U, V, WS, WS_Diff
REAL :: Min_Diff
TYPE(PP_Field_Type) :: WindSpd

REAL, &
DIMENSION(PFields(1)%Hdr%NumCols, PFields(1)%Hdr%NumRows, 2*ninc):: &
     Uinc, Vinc, Pinc, Winc

! Threshold in knots for Max Wind Top and Base Fields
Real, Parameter :: Threshold_Knots = 80.0
Real, Parameter :: Conv_knots_ms = 0.5418
Real, Parameter :: Threshold_ms = Threshold_Knots * Conv_knots_ms

! End of Header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

! --------------
! Max Winds Base
! --------------

MaxWindBase % Hdr = MaxWindP % Hdr

MaxWindBase % Hdr % PPCode   = PP_P
MaxWindBase % Hdr % MO8Type  = MO8_MxWBase
MaxWindBase % Hdr % STcode   = ST_MWBase
MaxWindBase % Hdr % MO8Level = LV_Special
MaxWindBase % Hdr % LBVC     = VC_Lower

ALLOCATE( MaxWindBase % RData(MaxWindBase % Hdr % NumCols, &
                              MaxWindBase % Hdr % NumRows) )
!
! -------------
! Max Winds Top
! -------------

MaxWindTop % Hdr = MaxWindP % Hdr

MaxWindTop % Hdr % PPCode   = PP_P
MaxWindTop % Hdr % MO8Type  = MO8_MxWTop
MaxWindTop % Hdr % STcode   = ST_MWTop
MaxWindTop % Hdr % MO8Level = LV_Special
MaxWindTop % Hdr % LBVC     = VC_Upper

ALLOCATE( MaxWindTop % RData(MaxWindTop % Hdr % NumCols, &
                             MaxWindTop % Hdr % NumRows) )

LevBtm = khalf
LevTop = NumLevs+1 - khalf

NULLIFY( WindSpd % RData )

! Initialise MaxWindP fields for base/top to RMDI.
! Only points with MaxWindSpeed >= Threshold will have data

MaxWindBase % RData = RMDI
MaxWindTop  % RData = RMDI

! Calculate Windspeed from MaxWindU and MaxWindV

! DEPENDS ON: vecmag
Call VecMag (MaxWindU, MaxWindV, WindSpd, ErrorStatus)

Lev_Below_MaxW(:,:) = 0
Lev_Above_MaxW(:,:) = 0

! Now know MaxWindU, MaxWindV, WindSpeed, MaxWindP
!          U, V, P (all levels)

! If MaxWindSpeed < 80 knots, then P Level above/below remains RMDI

Do j=1, MaxWindU % Hdr % NumRows
  Do i=1, MaxWindU % Hdr % NumCols

    If ( WindSpd % RData (i,j) >= Threshold_ms ) Then

      ! --------------------
      ! Find P below & above
      ! --------------------

      Do k = NumLevs,1,-1

        If ( PFields(k) % RData (i,j) > MaxWindP % RData (i,j) ) Then
          Level_Below = k
          Level_Above = k+1
          Exit
        End If

      End Do

      ! Now know model level below and above MaxWindP
      ! Now search for model level nearest to Max Wind Base and Top

      ! -------------
      ! Max Wind Base
      ! -------------

      Min_Diff = 10000.0

      Do k = Level_Below, LevBtm, -1

        U  = UFields (k) % RData (i,j)
        V  = VFields (k) % RData (i,j)
        WS = SQRT ( U**2 + V**2 )

        WS_Diff = ABS ( WS - Threshold_ms )

        If (WS_Diff < Min_Diff) Then
          Min_Diff = WS_Diff
          Lev_Below_MaxW (i,j) = k
        End If

        If (WS < Threshold_ms) Then
          Exit
        End If

      End Do

      ! -------------
      ! Max Wind Top
      ! -------------

      Min_Diff = 10000.0

      Do k = Level_Above, LevTop

        U  = UFields (k) % RData (i,j)
        V  = VFields (k) % RData (i,j)
        WS = SQRT ( U**2 + V**2 )

        WS_Diff = ABS ( WS - Threshold_ms )

        If (WS_Diff < Min_Diff) Then
          Min_Diff = WS_Diff
          Lev_Above_MaxW (i,j) = k
        End If

        If (WS < Threshold_ms) Then
          Exit
        End If

      End Do

   End If   !  If WindSp >= Threshold

  End Do
End Do

! -------------
! Max Wind Base
! -------------

! DEPENDS ON: maxwindsplinev
Call MaxWindSplineV (MaxWindU % Hdr % NumCols,                  &
                     MaxWindU % Hdr % NumRows,                  &
                     NumLevs,                                   &
                     UFields, VFields, PFields, Lev_Below_MaxW, &
                     Uinc, Vinc, Pinc )

! Look along increments for new values for max wind and level
do k=1, 2*ninc
   do j=1, PFields(1)%Hdr%NumRows
      do i=1, PFields(1)%Hdr%NumCols
          if (Lev_Below_MaxW(i,j)/=0) then
             Winc(i,j,k) = SQRT(Uinc(i,j,k)**2 + Vinc(i,j,k)**2)
             Winc(i,j,k) = Winc(i,j,k) - Threshold_ms
             Winc(i,j,k) = ABS ( Winc(i,j,k) )
          end if
      end do
   end do
end do

MaxW   (:,:) = Winc(:,:,1)
MaxWinc(:,:) = 1
DO k=2,2*ninc
  WHERE(Winc(:,:,k) < MaxW )
    MaxW    = Winc(:,:,k)
    MaxWinc = k
  END WHERE
ENDDO

Do j=1, MaxWindU % Hdr % NumRows
  Do i=1, MaxWindU % Hdr % NumCols
    IF(Lev_Below_MaxW(i,j)>0) THEN
      MaxWindBase % RData(i,j) = Pinc(i,j,MaxWinc(i,j))
    ENDIF
  End Do
End Do

! ------------
! Max Wind Top
! ------------

! DEPENDS ON: maxwindsplinev
Call MaxWindSplineV (MaxWindU % Hdr % NumCols,                  &
                     MaxWindU % Hdr % NumRows,                  &
                     NumLevs,                                   &
                     UFields, VFields, PFields, Lev_Above_MaxW, &
                     Uinc, Vinc, Pinc )

! Look along increments for new values for max wind and level

do k=1, 2*ninc
   do j=1, PFields(1)%Hdr%NumRows
      do i=1, PFields(1)%Hdr%NumCols
         if (Lev_Below_MaxW(i,j)/=0) then
            Winc(i,j,k) = SQRT(Uinc(i,j,k)**2 + Vinc(i,j,k)**2)
            Winc(i,j,k) = Winc(i,j,k) - Threshold_ms
            Winc(i,j,k) = ABS ( Winc(i,j,k) )
         end if
      end do
   end do
end do

MaxW   (:,:) = Winc(:,:,1)
MaxWinc(:,:) = 1
DO k=2,2*ninc
  WHERE(Winc(:,:,k) < MaxW )
    MaxW    = Winc(:,:,k)
    MaxWinc = k
  END WHERE
ENDDO

Do j=1, MaxWindU % Hdr % NumRows
  Do i=1, MaxWindU % Hdr % NumCols
    IF(Lev_Above_MaxW(i,j)>0) THEN
      MaxWindTop % RData(i,j) = Pinc(i,j,MaxWinc(i,j))
    ENDIF
  End Do
End Do

DEALLOCATE ( WindSpd % RData )

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE TopBase

! ====================================================================



! ====================================================================



!=======================================================================


!=======================================================================


!=======================================================================

#endif
