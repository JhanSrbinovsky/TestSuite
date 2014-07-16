
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Converts between P and exner (on rho levels)

Module Rcf_Exner_P_Convs_Mod

!  Subroutine Rcf_Conv_Exner_P - converts exner to P
!  Subroutine Rcf_Conv_P_Exner - converts P to exner
!
! Description:
!   Performs conversions between exner and P on rho levels.
!
! Method:
!   Data is stored in the *original* field data pointer!
!   Derived from New Dynamics 2.7 code.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   23/10/01   Consistent stashmaster handling. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

!******************************************************************
! This routine converts Exner to P
!******************************************************************
Subroutine Rcf_Conv_Exner_P( exner_field )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_exner,           &
    stashcode_p

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Rcf_Submodel_Mod, Only : &
    Atmos_IM

Implicit none

! Arguments
Type( field_type), Intent( InOut )    :: exner_field

! Comdecks
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

!  Local variables
Integer                        :: i
Integer                        :: k
Integer                        :: ErrorStatus
Real, Parameter                :: recip_kappa=1./kappa
Character (Len=*), Parameter   :: RoutineName = 'Rcf_Conv_Exner_P'
Character (Len=80)             :: Cmessage

!-----------------------------------------------------------------
! Make sure field actually is exner...
!-----------------------------------------------------------------
If ( exner_field % stashmaster % item /= stashcode_exner ) Then
  ErrorStatus = 10
  Cmessage = 'Stashcode does not match expected (exner) data'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!------------------------------------------------------------------
! Convert exner to P
!------------------------------------------------------------------
Do k = 1, exner_field % levels
  Do i = 1, exner_field % level_size
    exner_field % Data(i,k) =                                       &
                  (exner_field % Data(i,k) ** recip_kappa)* pref
  End Do
End Do

!------------------------------------------------------------------
! Set field stashmaster to P
!------------------------------------------------------------------
exner_field % stashmaster => Rcf_Exppx( Atmos_IM, 0, stashcode_p )

Return
End Subroutine Rcf_Conv_Exner_P


!*******************************************************************
! Routine to convert P to exner
!*******************************************************************

Subroutine Rcf_Conv_P_Exner( exner_field )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_exner,           &
    stashcode_p

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Rcf_Submodel_Mod, Only : &
    Atmos_IM

Implicit none

! Arguments
Type( field_type), Intent( InOut )    :: exner_field

! Comdecks
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

!  Local variables
Integer                        :: i
Integer                        :: k
Integer                        :: ErrorStatus
Real, Parameter                :: recip_p_zero=1./pref
Character (Len=*), Parameter   :: RoutineName = 'Rcf_Conv_P_Exner'
Character (Len=80)             :: Cmessage

!-----------------------------------------------------------------
! Make sure field actually is P...
!-----------------------------------------------------------------
If ( exner_field % stashmaster % item /= stashcode_p ) Then
  ErrorStatus = 10
  Cmessage = 'Stashcode does not match expected (pressure) data'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!------------------------------------------------------------------
! Convert P to Exner
!------------------------------------------------------------------
Do k = 1, exner_field % levels
  Do i = 1, exner_field % level_size
    exner_field % Data(i,k) =                                       &
                  (exner_field % Data(i,k) * recip_p_zero) ** kappa
  End Do
End Do

!------------------------------------------------------------------
! Set field stashmaster to exner
!------------------------------------------------------------------
exner_field % stashmaster => Rcf_Exppx(Atmos_IM, 0, stashcode_exner)

Return
End Subroutine Rcf_Conv_P_Exner

End Module Rcf_Exner_P_Convs_Mod
