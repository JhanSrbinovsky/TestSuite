
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ read the recon namelist

Module Rcf_readnl_Recon_Mod

!  Subroutine Rcf_Readnl_Recon - read the recon namelist
!
! Description:
!   Read the recon namelist and set Output_Grid variables accordingly.
!
! Method:
!   Variables read into the Rcf_Recon_Mod module.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_readnl_recon( nft )

Use Rcf_Recon_Mod

Use Rcf_Grid_Type_Mod, Only : Output_Grid
Use Rcf_Lsm_Mod, Only       : glob_land_out
Use Rcf_Parvars_Mod, Only   : mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_Parvars_Mod, Only : &
    mype

Use Ereport_Mod, Only : &
    Ereport

Implicit None

! Arguments
Integer                      :: nft     ! Unit number

! Local variables
Integer                      :: ErrorStatus
Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Readnl_Recon'

! Initialisation
Grib             = .FALSE.
Var_Recon        = .FALSE.
Uars             = .FALSE.
L_IO_Timer       = .FALSE.
p_rows           = 0
row_length       = 0
model_levels     = 0
wet_levels       = 0
tr_vars          = 0
tr_ukca          = 0
ozone_levels     = 0
st_levels        = 0
sm_levels        = 0
bl_levels        = 0
tr_levels        = 0
conv_levels      = 0
w_zero_start     = -1
w_zero_end       = -1
rimwidtha        = 0
rimwidtho        = 0
q_min            = 0.0
river_row_length = 360    ! 1 degree default
river_rows       = 180    ! 1 degree default
use_smc_stress   = .FALSE.

! Read namelist
Read( Unit = nft, Nml = recon )
If (PrintStatus >= PrStatus_Oper .AND. mype == 0) Then
  Write( 6 , Nml = recon )
End If

! Set variables as required
Output_Grid % glob_p_rows       = p_rows
Output_Grid % glob_p_row_length = row_length
Output_Grid % glob_r_rows       = river_rows
Output_Grid % glob_r_row_length = river_row_length
Output_Grid % model_levels      = model_levels
Output_Grid % wet_levels        = wet_levels
Output_Grid % ozone_levels      = ozone_levels
Output_Grid % st_levels         = st_levels
Output_Grid % sm_levels         = sm_levels
Output_Grid % bl_levels         = bl_levels
Output_Grid % tr_levels         = tr_levels
Output_Grid % conv_levels       = conv_levels

! If either w_zero_start or w_zero_end are -1 then set to
! model_levels (from same namelist) as default value.
If (w_zero_start == -1) Then
  w_zero_start = model_levels
End If
If (w_zero_end == -1) Then
  w_zero_end = model_levels
End If
If (w_zero_start == 0) Then  !  W at surface is always set to zero.
  w_zero_start = 1
End If

! Check validity of w_zero_start and w_zero_end
If ( (w_zero_start > w_zero_end) .or.                             &
     (w_zero_start < 0 .or. w_zero_start > model_levels) .or.     &
     (w_zero_end   < 0 .or. w_zero_end   > model_levels) ) Then
  write (6,*) 'w_zero_start ',w_zero_start,' w_zero_end ',w_zero_end
  ErrorStatus = 10
  write (CMessage,*) 'w_zero_start and/or w_zero_end set incorrectly.'
  Call Ereport ( RoutineName, ErrorStatus, cmessage)
End If

! This is a default - not really needed in the dump
Output_Grid % cloud_levels      = IMDI

! Also set some land/sea mask values
glob_land_out                   = land_points
Output_Grid % glob_land_field   = land_points

! Set values for tropopause-based ozone.
! In future developments, obtain directly from namelist.
tpps_ozone_levels  = ozone_levels
l_tpps_ozone_zonal = lozone_zonal

Return
End Subroutine Rcf_readnl_recon
End Module Rcf_readnl_Recon_Mod
