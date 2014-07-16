
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Set up Fsat/Fwetl for large-scale hydrology

Module Rcf_Calc_Fsat_Mod

! Subroutine Rcf_Calc_Fsat
!
! Description:
!   Initialises the surface saturation and/or wetland fraction
!
! Method:
!   The fields are calculated from the initial soil parameters and
!   soil moisture content and water table depth.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.5   17/01/03   Original code.  Nic Gedney.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Calc_Fsat( fields_out, field_count_out, &
                         stashcode_in, frac_field, hdr_out )

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_HeadAddress_Mod, Only : &
    IC_SoilMoistLevs

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_clapp_hb,        &
    stashcode_ksat,            &
    stashcode_vol_smc_sat,     &
    stashcode_unfrozen_soil,   &
    stashcode_frozen_soil,     &
    stashcode_ti_mean,         &
    stashcode_ti_sig,          &
    stashcode_fexp,            &
    stashcode_gamtot,          &
    stashcode_zw,              &
    stashcode_fsat,            &
    stashcode_fwetl,           &
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

Use Ereport_mod, Only : &
    Ereport

Use Rcf_CntlAtm_Mod, Only : &
    L_TOP

Implicit None

! Arguments
Type( field_type ), Pointer          :: fields_out(:)
Type( field_type ), Intent( InOut )  :: frac_field
Type( um_header_type ), Intent( In ) :: hdr_out

Integer, Intent( In )                :: field_count_out
Integer, Intent( In )                :: stashcode_in

! Local variables
Type( field_type ), Pointer          :: clapp_hb
Type( field_type ), Pointer          :: ksat
Type( field_type ), Pointer          :: vol_smc_sat
Type( field_type ), Pointer          :: unfrozen_soil
Type( field_type ), Pointer          :: frozen_soil
Type( field_type ), Pointer          :: ti_mean
Type( field_type ), Pointer          :: ti_sig
Type( field_type ), Pointer          :: fexp
Type( field_type ), Pointer          :: gamtot
Type( field_type ), Pointer          :: zw

Integer                              :: soil_index  &
                                       (frac_field % level_size)
Integer                              :: soil_pts
Integer                              :: size
Integer                              :: pos  ! field position
Integer                              :: i    ! looper
Integer                              :: j    ! looper
Integer                              :: n    ! looper

Real          :: wutot    (frac_field % level_size)
Real          :: top_crit (frac_field % level_size)
Real          :: fsat     (frac_field % level_size)
Real          :: fwetl    (frac_field % level_size)
Real          :: qbase    (frac_field % level_size)

Real          :: qbase_l (frac_field % level_size,                    &
                               hdr_out % IntC( IC_SoilMoistLevs ) +1 )
Real          :: Ksz     (frac_field % level_size,                    &
                               0: hdr_out % IntC( IC_SoilMoistLevs ) )
Real          :: zdepth  (0: hdr_out % IntC( IC_SoilMoistLevs) )

Integer       :: ErrorStatus

Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Calc_Fsat'

Logical :: L_Gamtot=.false.

! Comdecks
!---Soil layer thicknesses (m)
!---6 layers => CABLE else revert to MOSES
   REAL,PARAMETER:: DZSOIL(6) =(/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)


!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Setting up surface saturation/wetland fraction'
End If

If (.NOT.L_TOP) Then
   ErrorStatus = 10
   Write(CMessage,*) 'Field only required in TOPMODEL-based hydrology'
   call Ereport ( RoutineName, ErrorStatus, CMessage)
End If

!----------------------------------------------------------------------
! Find required fields in output dump and read them in:
!----------------------------------------------------------------------
  Call Rcf_Locate( stashcode_prog_sec, stashcode_clapp_hb,           &
                   fields_out, field_count_out, pos)
  clapp_hb => fields_out(pos)
  Call Rcf_Alloc_Field( clapp_hb )
  Call Rcf_Read_Field( clapp_hb, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_ksat,               &
                   fields_out, field_count_out, pos)
  ksat => fields_out(pos)
  Call Rcf_Alloc_Field( ksat )
  Call Rcf_Read_Field( ksat, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_vol_smc_sat,        &
                   fields_out, field_count_out, pos)
  vol_smc_sat => fields_out(pos)
  Call Rcf_Alloc_Field( vol_smc_sat )
  Call Rcf_Read_Field( vol_smc_sat, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_unfrozen_soil,      &
                   fields_out, field_count_out, pos)
  unfrozen_soil => fields_out(pos)
  Call Rcf_Alloc_Field( unfrozen_soil )
  Call Rcf_Read_Field( unfrozen_soil, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_frozen_soil,        &
                   fields_out, field_count_out, pos)
  frozen_soil => fields_out(pos)
  Call Rcf_Alloc_Field( frozen_soil )
  Call Rcf_Read_Field( frozen_soil, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_ti_mean,            &
                   fields_out, field_count_out, pos)
  ti_mean => fields_out(pos)
  Call Rcf_Alloc_Field( ti_mean )
  Call Rcf_Read_Field( ti_mean, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_ti_sig,             &
                   fields_out, field_count_out, pos)
  ti_sig => fields_out(pos)
  Call Rcf_Alloc_Field( ti_sig )
  Call Rcf_Read_Field( ti_sig, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_fexp,               &
                   fields_out, field_count_out, pos)
  fexp => fields_out(pos)
  Call Rcf_Alloc_Field( fexp )
  Call Rcf_Read_Field( fexp, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_gamtot,             &
                   fields_out, field_count_out, pos)
  gamtot => fields_out(pos)
  Call Rcf_Alloc_Field( gamtot )
  Call Rcf_Read_Field( gamtot, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_zw,                 &
                   fields_out, field_count_out, pos)
  zw => fields_out(pos)
  Call Rcf_Alloc_Field( zw )
  Call Rcf_Read_Field( zw, hdr_out, decomp_rcf_output )

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
  End do

  zdepth(0) = 0.0
  Do n = 1 , unfrozen_soil % levels
    zdepth(n) = zdepth(n-1) + dzsoil(n)
  End Do

!----------------------------------------------------------------------
! Set up saturated hydraulic conductivity for each soil layer:
!----------------------------------------------------------------------
  Do j = 1 , soil_pts
    i = soil_index(j)
    Do n = 0 , unfrozen_soil % levels
      Ksz(i,n) = ksat % Data(i,1)
    End do
  End do

!----------------------------------------------------------------------
! Initialise variables to zero:
!----------------------------------------------------------------------
  fsat(:) = 0
  fwetl(:) = 0
  qbase(:) = 0
  qbase_l(:,:) = 0
  frac_field % Data (:,:) = 0

!----------------------------------------------------------------------
! Calculate Baseflow:
!----------------------------------------------------------------------
! DEPENDS ON: calc_baseflow
      Call Calc_Baseflow(                    &
                      soil_pts,              &
                      soil_index,            &
                      ti_mean % level_size,  &
                      unfrozen_soil % levels,&
                      zdepth,                &
                      Ksz,                   &
                      clapp_hb % Data,       &
                      fexp % Data,           &
                      ti_mean % Data,        &
                      zw % Data,             &
                      frozen_soil % Data,    &
                      unfrozen_soil % Data,  &
                      wutot,                 &
                      top_crit,              &
                      qbase,                 &
                      qbase_l)

!----------------------------------------------------------------------
! Calculate surface saturation and wetland fraction:
!----------------------------------------------------------------------
! DEPENDS ON: calc_fsat
      Call Calc_Fsat(                        &
                      L_gamtot,              &
                      soil_pts,              &
                      soil_index,            &
                      ti_mean % level_size,  &
                      ti_mean % Data,        &
                      ti_sig % Data,         &
                      wutot,                 &
                      top_crit,              &
                      gamtot % Data,         &
                      fsat,                  &
                      fwetl)


!----------------------------------------------------------------------
! Set up output to surface saturation or wetland fraction:
!----------------------------------------------------------------------
        If (stashcode_in == stashcode_fsat) Then
          frac_field % Data(:,1)=fsat(:)
        End If
        If (stashcode_in == stashcode_fwetl) Then
          frac_field % Data(:,1)=fwetl(:)
        End If

!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------
  Call Rcf_Dealloc_Field( clapp_hb )
  Call Rcf_Dealloc_Field( ksat )
  Call Rcf_Dealloc_Field( vol_smc_sat )
  Call Rcf_Dealloc_Field( unfrozen_soil )
  Call Rcf_Dealloc_Field( frozen_soil )
  Call Rcf_Dealloc_Field( ti_mean )
  Call Rcf_Dealloc_Field( ti_sig )
  Call Rcf_Dealloc_Field( fexp )
  Call Rcf_Dealloc_Field( gamtot )
  Call Rcf_Dealloc_Field( zw )

End Subroutine Rcf_Calc_Fsat

End Module Rcf_Calc_Fsat_Mod
