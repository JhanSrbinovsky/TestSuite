
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Set up Gamtot for large-scale hydrology

Module Rcf_Calc_Gamtot_Mod

! Subroutine Rcf_Calc_Gamtot
!
! Description:
!   Gets the total gamma function for topographic index for use
!   as a scale factor in the large-scale hydrology scheme.
!
! Method:
!   Calls routine CALC_GAMMA
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

Subroutine Rcf_Calc_Gamtot( fields_out, field_count_out, &
                             gamtot, hdr_out )

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_gamtot,          &
    stashcode_vol_smc_sat,     &
    stashcode_ti_mean,         &
    stashcode_ti_sig,          &
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
Type( field_type ), Intent( InOut )  :: gamtot
Type( um_header_type ), Intent( In ) :: hdr_out

Integer, Intent( In )                :: field_count_out

! Local variables
Type( field_type ), Pointer          :: ti_mean
Type( field_type ), Pointer          :: ti_sig
Type( field_type ), Pointer          :: vol_smc_sat

Real                                 :: dummy  &
                                       (gamtot % level_size)


Integer                              :: soil_index  &
                                       (gamtot % level_size)
Integer                              :: soil_pts
Integer                              :: size
Integer                              :: pos  ! field position
Integer                              :: i    ! looper

Integer      :: ErrorStatus

Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Calc_Gamtot'

Logical :: L_Gamtot = .true.

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Setting up Gamtot from topographic variables'
End If

If (.NOT.L_TOP) Then
   ErrorStatus = 10
   Write(CMessage,*) 'Field only required in TOPMODEL-based hydrology'
   call Ereport ( RoutineName, ErrorStatus, CMessage)
End If

!----------------------------------------------------------------------
! Find required fields in output dump and read them in:
!----------------------------------------------------------------------
  Call Rcf_Locate( stashcode_prog_sec, stashcode_vol_smc_sat,        &
                   fields_out, field_count_out, pos)
  vol_smc_sat => fields_out(pos)
  Call Rcf_Alloc_Field( vol_smc_sat )
  Call Rcf_Read_Field( vol_smc_sat, hdr_out, decomp_rcf_output )

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

!----------------------------------------------------------------------
! Calculate gamtot:
!----------------------------------------------------------------------

  gamtot % Data(:,:) = 0.0

! DEPENDS ON: calc_fsat
  Call Calc_Fsat( L_gamtot,              &
                  soil_pts,              &
                  soil_index,            &
                  ti_mean % level_size,  &
                  ti_mean % Data,        &
                  ti_sig % Data,         &
                  dummy,                 &
                  dummy,                 &
                  gamtot % Data,         &
                  dummy,                 &
                  dummy)

!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------

  Call Rcf_Dealloc_Field( vol_smc_sat )
  Call Rcf_Dealloc_Field( ti_mean )
  Call Rcf_Dealloc_Field( ti_sig )

End Subroutine Rcf_Calc_Gamtot

End Module Rcf_Calc_Gamtot_Mod
