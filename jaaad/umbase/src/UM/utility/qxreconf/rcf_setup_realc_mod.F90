#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up the real constants in the output header

Module Rcf_setup_RealC_Mod

!  Subroutine Rcf_Setup_RealC - initialises real constants in header.
!
! Description:
!   Sets the real constants in the output dump header.
!
! Method:
!    Uses namelist information to set header.
!    UMDP F3 defines the real constants.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Setup_RealC( Hdr_In, Hdr_Out )

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_UMhead_Mod, Only : &
    Um_Header_type

Use rcf_headers_mod, Only : &
    RelHd

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_Readnl_Horizont_Mod, Only :             &
    Delta_Lambda,                 Delta_Phi,    &
    Phi_First,                    Lambda_First, &
    Phi_Npole,                    Lambda_Npole

Use Rcf_HeadAddress_Mod, Only :                        &
    RC_PoleLong,                  RC_LongSpacing,  &
    RC_LatSpacing,                RC_PoleLat,      &
    RC_FirstLong,                 RC_FirstLat,     &
    RC_ModelTop,                  RC_SWLDEG,       &
    RC_WEdgeDeg

Implicit None

! Arguments
Type (Um_Header_Type), Intent( In ) :: Hdr_In
Type (Um_Header_Type), Target       :: Hdr_Out

! Comdecks
#include "c_mdi.h"

! Local vars
Real, Pointer                  :: RealC(:)
Integer                        :: i

!---------------------------------------------------------------
! Clean start - RMDI for everything
!---------------------------------------------------------------
RealC => Hdr_Out % RealC
RealC(:) = RMDI

!--------------------------------------------------------------
! Set values we have numbers for
!--------------------------------------------------------------

RealC( RC_LongSpacing ) = Delta_Lambda
RealC( RC_LatSpacing  ) = Delta_Phi
RealC( RC_FirstLat    ) = Phi_First
RealC( RC_FirstLong   ) = Lambda_First
RealC( RC_PoleLat     ) = Phi_Npole
RealC( RC_PoleLong    ) = Lambda_Npole

! Submodel specifics
RealC( RC_ModelTop    ) = Output_Grid % z_top_of_model

!--------------------------------------------------------------
! Overrides from namelists
!--------------------------------------------------------------
Do i = 1, Hdr_Out % LenRealC
  If( RelHd(i) /= RMDI ) Then
    If (PrintStatus >= PrStatus_Oper) Then
      Write (6,*) 'RealC(',i,') has been reset from ', RealC(i), &
                  ' to ', RelHd(i)
    End If

    RealC(i) = RelHd(i)
  End If
End Do

Return
End Subroutine Rcf_Setup_RealC
End Module Rcf_Setup_RealC_Mod
#endif
