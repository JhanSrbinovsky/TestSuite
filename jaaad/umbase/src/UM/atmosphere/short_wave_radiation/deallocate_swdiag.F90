#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Deallocate the SW diagnostics
!
! Current Code Owner: J.-C. Thelen
!
! History:
! Version     Date        Comment
! 6.2         19/01/06    Original Code (J.-C. Thelen
!
! Description of Code:
!    FORTRAN 77/90  with extensions listed in documentation.
!
!-----------------------------------------------------------------------
!
Subroutine deallocate_swdiag(j_sw)

  Use sw_diag_mod
  Implicit None

  Integer, Intent(in) :: j_sw              ! call to SW radiation

!
! Deallocate SW Diagnostics
!
  If (Associated(SW_diag(j_sw)%solar_out_toa)) &
    Deallocate(SW_diag(j_sw)%solar_out_toa)
  If (Associated(SW_diag(j_sw)%solar_out_clear)) &
    Deallocate(SW_diag(j_sw)%solar_out_clear)
  If (Associated(SW_diag(j_sw)%surface_down_flux)) &
    Deallocate(SW_diag(j_sw)%surface_down_flux)
  If (Associated(SW_diag(j_sw)%surf_down_clr)) &
    Deallocate(SW_diag(j_sw)%surf_down_clr)
  If (Associated(SW_diag(j_sw)%surf_up_clr)) &
    Deallocate(SW_diag(j_sw)%surf_up_clr)
  If (Associated(SW_diag(j_sw)%net_flux_trop)) &
    Deallocate(SW_diag(j_sw)%net_flux_trop)
  If (Associated(SW_diag(j_sw)%up_flux_trop)) &
    Deallocate(SW_diag(j_sw)%up_flux_trop)
  If (Associated(SW_diag(j_sw)%clear_hr)) &
    Deallocate(SW_diag(j_sw)%clear_hr)
! Radiances
  If (Associated(SW_diag(j_sw)%toa_radiance)) &
    Deallocate(SW_diag(j_sw)%toa_radiance)
! UV-Fluxes
  If (Associated(SW_diag(j_sw)%uvflux_direct)) &
    Deallocate(SW_diag(j_sw)%uvflux_direct)
  If (Associated(SW_diag(j_sw)%uvflux_up)) &
    Deallocate(SW_diag(j_sw)%uvflux_up)
  If (Associated(SW_diag(j_sw)%uvflux_net)) &
    Deallocate(SW_diag(j_sw)%uvflux_net)
! Direct and diffuse downward SW Fluxes
  If (Associated(SW_diag(j_sw)%flux_direct)) &
    Deallocate(SW_diag(j_sw)%flux_direct)
  If (Associated(SW_diag(j_sw)%flux_diffuse)) &
    Deallocate(SW_diag(j_sw)%flux_diffuse)
! Microphysical diagnostics
  If (Associated(SW_diag(j_sw)%re_strat)) &
    Deallocate(SW_diag(j_sw)%re_strat)
  If (Associated(SW_diag(j_sw)%wgt_strat)) &
    Deallocate(SW_diag(j_sw)%wgt_strat)
  If (Associated(SW_diag(j_sw)%lwp_strat)) &
    Deallocate(SW_diag(j_sw)%lwp_strat)
  If (Associated(SW_diag(j_sw)%re_conv)) &
    Deallocate(SW_diag(j_sw)%re_conv)
  If (Associated(SW_diag(j_sw)%wgt_conv)) &
    Deallocate(SW_diag(j_sw)%wgt_conv)
  If (Associated(SW_diag(j_sw)%ntot_diag)) &
    Deallocate(SW_diag(j_sw)%ntot_diag)
  If (Associated(SW_diag(j_sw)%strat_lwc_diag)) &
    Deallocate(SW_diag(j_sw)%strat_lwc_diag)
  If (Associated(SW_diag(j_sw)%so4_ccn_diag)) &
    Deallocate(SW_diag(j_sw)%so4_ccn_diag)
  If (Associated(SW_diag(j_sw)%cond_samp_wgt)) &
    Deallocate(SW_diag(j_sw)%cond_samp_wgt)
  If (Associated(SW_diag(j_sw)%weighted_re)) &
    Deallocate(SW_diag(j_sw)%weighted_re)
  If (Associated(SW_diag(j_sw)%sum_weight_re)) &
    Deallocate(SW_diag(j_sw)%sum_weight_re)
  If (Associated(SW_diag(j_sw)%weighted_warm_re)) &
    Deallocate(SW_diag(j_sw)%weighted_warm_re)
  If (Associated(SW_diag(j_sw)%sum_weight_warm_re)) &
    Deallocate(SW_diag(j_sw)%sum_weight_warm_re)
  If (Associated(SW_diag(j_sw)%Nc_diag)) &
    Deallocate(SW_diag(j_sw)%Nc_diag)
  If (Associated(SW_diag(j_sw)%Nc_weight)) &
    Deallocate(SW_diag(j_sw)%Nc_weight)
! Diagnostics for MOSES II
  If (Associated(SW_diag(j_sw)%FlxSolBelow690nmSurf)) &
    Deallocate(SW_diag(j_sw)%FlxSolBelow690nmSurf)
  If (Associated(SW_diag(j_sw)%FlxSeaBelow690nmSurf)) &
    Deallocate(SW_diag(j_sw)%FlxSeaBelow690nmSurf)
! STOCHEM Diagnostics
  If (Associated(SW_diag(j_sw)%flxdirparsurf)) &
    Deallocate(SW_diag(j_sw)%flxdirparsurf)
! Orography correction diagnostics:
  If (Associated(SW_diag(j_sw)%orog_corr)) &
    Deallocate(SW_diag(j_sw)%orog_corr)
  If (Associated(SW_diag(j_sw)%sol_bearing)) &
    Deallocate(SW_diag(j_sw)%sol_bearing)

End Subroutine deallocate_swdiag
#endif
