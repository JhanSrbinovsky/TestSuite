#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up Fitting parameters for large-scale hydrology

Module Rcf_Fit_Fsat_Mod

! Subroutine Rcf_Fit_Fsat
!
! Description:
!   Calls Calc_Fit_Fsat which calculates the fitting parameters for
!   LSH model and stores them.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Fit_Fsat( fields_out, field_count_out,        &
                         dum_out,a_fsat,c_fsat,              &
                         a_fwet,c_fwet,hdr_out )

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Stashcodes_Mod, Only :  &
    stashcode_sthzw,            &
    stashcode_vol_smc_sat,      &
    stashcode_unfrozen_soil,    &
    stashcode_ti_mean,          &
    stashcode_ti_sig,           &
    stashcode_gamtot,           &
    stashcode_fexp,             &
    stashcode_prog_sec

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Ereport_mod, Only : &
    Ereport

Use Rcf_CntlAtm_Mod, Only : &
    L_TOP


Implicit None

! Arguments
Type( field_type ), Pointer          :: fields_out(:)
Type( field_type ), Intent( InOut )  :: dum_out
Type( um_header_type ), Intent( In ) :: hdr_out

Real, Intent( InOut )  :: a_fsat &
                         (dum_out % levels,dum_out % level_size)
Real, Intent( InOut )  :: c_fsat &
                         (dum_out % levels,dum_out % level_size)
Real, Intent( InOut )  :: a_fwet &
                         (dum_out % levels,dum_out % level_size)
Real, Intent( InOut )  :: c_fwet &
                         (dum_out % levels,dum_out % level_size)

Integer, Intent( In )                :: field_count_out

! Local variables
Type( field_type ), Pointer          :: vol_smc_sat
Type( field_type ), Pointer          :: unfrozen_soil
Type( field_type ), Pointer          :: ti_mean
Type( field_type ), Pointer          :: ti_sig
Type( field_type ), Pointer          :: gamtot
Type( field_type ), Pointer          :: fexp

Integer                              :: soil_index  &
                                       (dum_out % level_size)
Integer                              :: soil_pts

Integer                              :: pos  ! field position
Integer                              :: i    ! looper
Integer                              :: n    ! looper
Integer                              :: j    ! looper

Integer      :: ErrorStatus

Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Fit_Fsat'

Real                                 :: soil_depth

! Comdecks:
#include "c_topog.h"
#include "soil_thick.h"

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Setting up call for calc_fit_fsat in rcf_fit_fsat'
End If

If (.NOT.L_TOP) Then
  ErrorStatus = 30
  Write(CMessage,*) 'Field only required in TOPMODEL-based hydrology'
  Call Ereport ( RoutineName, ErrorStatus, CMessage)
End If

!----------------------------------------------------------------------
! Find required fields in output dump and read them in:
!----------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_vol_smc_sat,            &
                 fields_out, field_count_out, pos)
vol_smc_sat => fields_out(pos)
Call Rcf_Alloc_Field( vol_smc_sat )
Call Rcf_Read_Field( vol_smc_sat, hdr_out, decomp_rcf_output )

Call Rcf_Locate( stashcode_prog_sec, stashcode_unfrozen_soil,          &
                 fields_out, field_count_out, pos)
unfrozen_soil => fields_out(pos)
Call Rcf_Alloc_Field( unfrozen_soil )
Call Rcf_Read_Field( unfrozen_soil, hdr_out, decomp_rcf_output )

Call Rcf_Locate( stashcode_prog_sec, stashcode_ti_mean,                &
                 fields_out, field_count_out, pos)
ti_mean => fields_out(pos)
Call Rcf_Alloc_Field( ti_mean )
Call Rcf_Read_Field( ti_mean, hdr_out, decomp_rcf_output )

Call Rcf_Locate( stashcode_prog_sec, stashcode_ti_sig,                 &
                 fields_out, field_count_out, pos)
ti_sig => fields_out(pos)
Call Rcf_Alloc_Field( ti_sig )
Call Rcf_Read_Field( ti_sig, hdr_out, decomp_rcf_output )

Call Rcf_Locate( stashcode_prog_sec, stashcode_gamtot,                 &
                 fields_out, field_count_out, pos)
gamtot => fields_out(pos)
Call Rcf_Alloc_Field( gamtot )
Call Rcf_Read_Field( gamtot, hdr_out, decomp_rcf_output )

Call Rcf_Locate( stashcode_prog_sec, stashcode_fexp,                   &
                 fields_out, field_count_out, pos)
fexp => fields_out(pos)
Call Rcf_Alloc_Field( fexp )
Call Rcf_Read_Field( fexp, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Initialise fitting variables to zero:
!----------------------------------------------------------------------
a_fsat(:,:) = 0.0
c_fsat(:,:) = 0.0
a_fwet(:,:) = 0.0
c_fwet(:,:) = 0.0

!----------------------------------------------------------------------
! Set up soil index:
!----------------------------------------------------------------------

soil_pts=0
soil_index(:)=0
Do i=1 , vol_smc_sat % level_size
  If (vol_smc_sat % Data(i,1) > 0.0) Then
    soil_pts = soil_pts + 1
    soil_index(soil_pts) = i
  End If
End Do

Soil_depth=0.0
Do n = 1,unfrozen_soil % levels
  soil_depth=soil_depth+dzsoil(n)
End Do

!----------------------------------------------------------------------
! Set up variables needed for call to calc_fit_fsat
!----------------------------------------------------------------------

Do j = 1, soil_pts
  i = soil_index(j)

  If (fexp % Data(i,1) <= 0.0) Then
    ErrorStatus = 20
    Write (CMessage,*)'Exp. decay in Vol_Smc_Sat wrongly/not set'
    Call Ereport ( RoutineName, ErrorStatus, CMessage)
  End If

End Do

!----------------------------------------------------------------------
! Get fitting parameters:
!----------------------------------------------------------------------

! DEPENDS ON: calc_fit_fsat
Call Calc_Fit_Fsat(                       &
                   soil_pts,              &
                   soil_index,            &
                   ti_mean % level_size,  &
                   fexp % Data,           &
                   ti_mean % Data,        &
                   ti_sig % Data,         &
                   gamtot % Data,         &
                   soil_depth,            &
                   a_fsat,                &
                   c_fsat,                &
                   a_fwet,                &
                   c_fwet)

!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------
Call Rcf_Dealloc_Field( vol_smc_sat )
Call Rcf_Dealloc_Field( unfrozen_soil )
Call Rcf_Dealloc_Field( ti_mean )
Call Rcf_Dealloc_Field( ti_sig )
Call Rcf_Dealloc_Field( gamtot )
Call Rcf_Dealloc_Field( fexp )

End Subroutine Rcf_Fit_Fsat

End Module Rcf_Fit_Fsat_Mod
#endif
