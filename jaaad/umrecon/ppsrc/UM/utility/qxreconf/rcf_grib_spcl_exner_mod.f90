
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate Exner from pressure levels

Module Rcf_Grib_Spcl_Exner_Mod

! SUBROUTINE Rcf_Grib_Spcl_Exner
!
! Description: Calculate Exner from pressure level.
!
! Method: Using (P/P0)^k
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     21/08/02   Original code. Roddy Sharp (frtz)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains
Subroutine Rcf_Grib_Spcl_Exner(Current,FpData)

Use Rcf_GRIB_Block_Params_Mod, Only : &
    Grib_Record,     LenArrayMax,        &
    p_Lvl_Type,      p_Lvl_Desc_1,       &
    Tb3_Pressure

Use EReport_Mod, Only :     &
    EReport

Implicit None

! Declarations:
! Global variables (#include statements etc):
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------

! Subroutine arguments

!< Scalar arguments with intent(In):>

!< Array  arguments with intent(InOut):>
Type (Grib_Record),Pointer          :: Current
Real, Intent(InOut)                 :: FpData(LenArrayMax)


! Local constants
Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_Spcl_Exner'

! Local variables

Integer                          :: Pressure
Real                             :: Exner
Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport

!=======================================================================
!  Routine Code Start :
!=======================================================================

!=======================================================================
!  Data on a Pressure level
!=======================================================================

If (Current % Block_1 (p_Lvl_Type) == Tb3_Pressure ) Then

! Get pressure value for level
  Pressure = Current % Block_1 (p_Lvl_Desc_1) * 100 ! convert from HPa

! Calculate Exner for Level
  FpData(:) = ( Pressure / P_zero ) ** kappa

!=======================================================================
!  Data _NOT_ on a Pressure level
!=======================================================================

Else

  Cmessage(1) = 'Cant handle level type other than Pressure'
  ErrorStatus = 10
  Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
End If
!=======================================================================
!  Code End :
!=======================================================================

Return

End Subroutine Rcf_Grib_Spcl_Exner
End Module Rcf_Grib_Spcl_Exner_Mod
