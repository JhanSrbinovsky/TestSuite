#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Allocate the required memory space for the radiation
!          diagnostics.
!
! Method: If a radiation diagnostic is required the correct amount
!         of memory space is allocated otherise a minimal amount
!         of space is allocated.
!         All radiation diagnostics are initialised to zero
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
Subroutine allocate_sw_diag(row_length, rows, model_levels,     &
                            cloud_levels,j_sw)

! Modules

   Use spec_sw_lw
   Use sw_diag_mod

   Implicit None

   Integer, Intent(in) :: row_length       ! Length of rows
   Integer, Intent(in) :: rows             ! Number of rows
   Integer, Intent(in) :: model_levels     ! Number of model levels
   Integer, Intent(in) :: cloud_levels     ! Number of cloud levels
   Integer, Intent(in) :: j_sw             ! call to SW radiation

   If (SW_diag(j_sw)%L_solar_out_toa) Then
      Allocate(SW_diag(j_sw)%solar_out_toa(row_length,rows))
   Else
      Allocate(SW_diag(j_sw)%solar_out_toa(1,1))
   End If
   SW_diag(j_sw)%solar_out_toa = 0.0
!
   If (SW_diag(j_sw)%L_solar_out_clear) Then
      Allocate(SW_diag(j_sw)%solar_out_clear(row_length,rows))
   Else
      Allocate(SW_diag(j_sw)%solar_out_clear(1,1))
   End If
   SW_diag(j_sw)%solar_out_clear = 0.0
!
   If (SW_diag(j_sw)%L_surface_down_flux) Then
      Allocate(SW_diag(j_sw)%surface_down_flux(row_length,rows))
   Else
      Allocate(SW_diag(j_sw)%surface_down_flux(1,1))
   End If
   SW_diag(j_sw)%surface_down_flux = 0.0
!
   If (SW_diag(j_sw)%L_surf_down_clr) Then
      Allocate(SW_diag(j_sw)%surf_down_clr(row_length,rows))
   Else
      Allocate(SW_diag(j_sw)%surf_down_clr(1,1))
   End If
   SW_diag(j_sw)%surf_down_clr = 0.0
!
   If (SW_diag(j_sw)%L_surf_up_clr) Then
      Allocate(SW_diag(j_sw)%surf_up_clr(row_length,rows))
   Else
      Allocate(SW_diag(j_sw)%surf_up_clr(1,1))
   End If
   SW_diag(j_sw)%surf_up_clr = 0.0
!
   If (SW_diag(j_sw)%L_net_flux_trop) Then
      Allocate(SW_diag(j_sw)%net_flux_trop(row_length,rows))
   Else
      Allocate(SW_diag(j_sw)%net_flux_trop(1,1))
   End If
   SW_diag(j_sw)%net_flux_trop = 0.0
!
   If (SW_diag(j_sw)%L_up_flux_trop) Then
      Allocate(SW_diag(j_sw)%up_flux_trop(row_length, rows))
   Else
      Allocate(SW_diag(j_sw)%up_flux_trop(1,1))
   End If
   SW_diag(j_sw)%up_flux_trop = 0.0
!
   If (SW_diag(j_sw)%L_clear_hr) Then
      Allocate(SW_diag(j_sw)%clear_hr(row_length, rows, model_levels))
   Else
      Allocate(SW_diag(j_sw)%clear_hr(1,1,1))
   End If
   SW_diag(j_sw)%clear_hr = 0.0

!
!  Radiance
!
   If (SW_diag(j_sw)%L_toa_radiance) Then
      Allocate(SW_diag(j_sw)%toa_radiance(row_length, rows,       &
                                      sw_spectrum(j_sw)%n_band))
   Else
      Allocate(SW_diag(j_sw)%toa_radiance(1,1,1))
   End If
   SW_diag(j_sw)%toa_radiance = 0.0
!
! Direct and Diffuse Downward Flux
!
   If (SW_diag(j_sw)%L_flux_direct) Then
      Allocate(SW_diag(j_sw)%flux_direct(row_length, rows,        &
                                         model_levels+1))
   Else
      Allocate(SW_diag(j_sw)%flux_direct(1,1,1))
   End If
   SW_diag(j_sw)%flux_direct = 0.0 
   
   If (SW_diag(j_sw)%L_flux_diffuse) Then
      Allocate(SW_diag(j_sw)%flux_diffuse(row_length, rows,       &
                                          model_levels+1))
   Else
      Allocate(SW_diag(j_sw)%flux_diffuse(1,1,1))
   End If
   SW_diag(j_sw)%flux_diffuse = 0.0  
!
! UV-Fluxes
!

   If (SW_diag(j_sw)%L_uvflux_direct) Then
      Allocate(SW_diag(j_sw)%uvflux_direct(row_length, rows,      &
                                          model_levels+1))
   Else
      Allocate(SW_diag(j_sw)%uvflux_direct(1,1,1))
   End If
   SW_diag(j_sw)%uvflux_direct = 0.0

   If (SW_diag(j_sw)%L_uvflux_up) Then
      Allocate(SW_diag(j_sw)%uvflux_up(row_length, rows,          &
                                          model_levels+1))
   Else
      Allocate(SW_diag(j_sw)%uvflux_up(1,1,1))
   End If
   SW_diag(j_sw)%uvflux_up = 0.0

   If (SW_diag(j_sw)%L_uvflux_net) Then
      Allocate(SW_diag(j_sw)%uvflux_net(row_length, rows,         &
                                          model_levels+1))
   Else
      Allocate(SW_diag(j_sw)%uvflux_net(1,1,1))
   End If
   SW_diag(j_sw)%uvflux_net = 0.0

!
! Microphysical diagnostics
!
   If (SW_diag(j_sw)%re_strat_flag) Then
     Allocate(SW_diag(j_sw)%re_strat(row_length, rows, cloud_levels))
   Else
     Allocate(SW_diag(j_sw)%re_strat(1,1,1))
   End If
   SW_diag(j_sw)%re_strat = 0.0

   If (SW_diag(j_sw)%wgt_strat_flag) Then
     Allocate(SW_diag(j_sw)%wgt_strat(row_length, rows, cloud_levels))
   Else
     Allocate(SW_diag(j_sw)%wgt_strat(1,1,1))
   End If
   SW_diag(j_sw)%wgt_strat = 0.0

   If (SW_diag(j_sw)%lwp_strat_flag) Then
     Allocate(SW_diag(j_sw)%lwp_strat(row_length, rows, cloud_levels))
   Else
     Allocate(SW_diag(j_sw)%lwp_strat(1,1,1))
   End If
   SW_diag(j_sw)%lwp_strat = 0.0

   If (SW_diag(j_sw)%re_conv_flag) Then
      Allocate(SW_diag(j_sw)%re_conv(row_length, rows, cloud_levels))
   Else
      Allocate(SW_diag(j_sw)%re_conv(1,1,1))
   End If
   SW_diag(j_sw)%re_conv = 0.0

   If (SW_diag(j_sw)%wgt_conv_flag) Then
      Allocate(SW_diag(j_sw)%wgt_conv(row_length, rows,cloud_levels))
   Else
      Allocate(SW_diag(j_sw)%wgt_conv(1,1,1))
   End If
   SW_diag(j_sw)%wgt_conv = 0.0

   If (SW_diag(j_sw)%ntot_diag_flag) Then
      Allocate(SW_diag(j_sw)%ntot_diag(row_length, rows, cloud_levels))
   Else
      Allocate(SW_diag(j_sw)%ntot_diag(1,1,1))
   End If
   SW_diag(j_sw)%ntot_diag = 0.0

   If (SW_diag(j_sw)%strat_lwc_diag_flag) Then
      Allocate(SW_diag(j_sw)%strat_lwc_diag(row_length, rows,      &
                                            cloud_levels))
   Else
      Allocate(SW_diag(j_sw)%strat_lwc_diag(1,1,1))
   End If
   SW_diag(j_sw)%strat_lwc_diag = 0.0

   If (SW_diag(j_sw)%so4_ccn_diag_flag) Then
      Allocate(SW_diag(j_sw)%so4_ccn_diag(row_length, rows,        &
                                          cloud_levels))
   Else
      Allocate(SW_diag(j_sw)%so4_ccn_diag(1,1,1))
   End If
   SW_diag(j_sw)%so4_ccn_diag = 0.0
!
   If (SW_diag(j_sw)%cond_samp_wgt_flag) Then
      Allocate(SW_diag(j_sw)%cond_samp_wgt(row_length, rows,       &
                                           cloud_levels))
   Else
      Allocate(SW_diag(j_sw)%cond_samp_wgt(1,1,1))
   End If
   SW_diag(j_sw)%cond_samp_wgt = 0.0
!
   If (SW_diag(j_sw)%weighted_re_flag) Then
      Allocate(SW_diag(j_sw)%weighted_re(row_length, rows))
   Else
      Allocate(SW_diag(j_sw)%weighted_re(1,1))
   End If
   SW_diag(j_sw)%weighted_re = 0.0
!
   If (SW_diag(j_sw)%sum_weight_re_flag) Then
      Allocate(SW_diag(j_sw)%sum_weight_re(row_length, rows))
   Else
      Allocate(SW_diag(j_sw)%sum_weight_re(1,1))
   End If
   SW_diag(j_sw)%sum_weight_re = 0.0
!
   If (SW_diag(j_sw)%wgtd_warm_re_flag) Then
      Allocate(SW_diag(j_sw)%weighted_warm_re(row_length, rows))
   Else
      Allocate(SW_diag(j_sw)%weighted_warm_re(1,1))
   End If
   SW_diag(j_sw)%weighted_warm_re = 0.0
!
   If (SW_diag(j_sw)%sum_wgt_warm_re_flag) Then
      Allocate(SW_diag(j_sw)%sum_weight_warm_re(row_length, rows))
   Else
      Allocate(SW_diag(j_sw)%sum_weight_warm_re(1,1))
   End If
   SW_diag(j_sw)%sum_weight_warm_re = 0.0
!
   If (SW_diag(j_sw)%Nc_diag_flag) Then
     Allocate(SW_diag(j_sw)%Nc_diag(row_length, rows))
   Else
      Allocate(SW_diag(j_sw)%Nc_diag(1,1))
   End If
   SW_diag(j_sw)%Nc_diag = 0.0
!
   If (SW_diag(j_sw)%Nc_weight_flag) Then
      Allocate(SW_diag(j_sw)%Nc_weight(row_length, rows))
   Else
      Allocate(SW_diag(j_sw)%Nc_weight(1,1))
   End If
   SW_diag(j_sw)%Nc_weight = 0.0
!
! Diagnostics for MOSES II:
!
   If (SW_diag(j_sw)%L_FlxSolBelow690nmSurf) Then
      Allocate(SW_diag(j_sw)%FlxSolBelow690nmSurf(row_length, rows))
   Else
      Allocate(SW_diag(j_sw)%FlxSolBelow690nmSurf(1,1))
   End If
   SW_diag(j_sw)%FlxSolBelow690nmSurf = 0.0
!
   If (SW_diag(j_sw)%L_FlxSeaBelow690nmSurf) Then
      Allocate(SW_diag(j_sw)%FlxSeaBelow690nmSurf(row_length, rows))
   Else
      Allocate(SW_diag(j_sw)%FlxSeaBelow690nmSurf(1,1))
   End If
   SW_diag(j_sw)%FlxSeaBelow690nmSurf = 0.0
!
!
! Direct component of PAR flux for STOCHEM
!
   If (SW_diag(j_sw)%L_direct_par) Then
      Allocate(SW_diag(j_sw)%flxdirparsurf(row_length, rows))
   Else
      Allocate(SW_diag(j_sw)%flxdirparsurf(1, 1))
   Endif
   SW_diag(j_sw)%flxdirparsurf = 0.0
!
! Orography correction diagnostics:
!
   If (SW_diag(j_sw)%L_orog_corr) Then
     Allocate(SW_diag(j_sw)%orog_corr(row_length, rows))
   Else
     Allocate(SW_diag(j_sw)%orog_corr(1,1))
   End If
   SW_diag(j_sw)%orog_corr = 0.0

   If (SW_diag(j_sw)%L_sol_bearing) Then
     Allocate(SW_diag(j_sw)%sol_bearing(row_length, rows))
   Else
     Allocate(SW_diag(j_sw)%sol_bearing(1,1))
   End If
   SW_diag(j_sw)%sol_bearing = 0.0
!
!
!        Extinction diagnostics:
!
   If (SW_diag(j_sw)%L_cloud_extinction) Then
      Allocate(SW_diag(j_sw)%cloud_extinction(row_length,rows,      &
                                              cloud_levels))
   Else
      Allocate(SW_diag(j_sw)%cloud_extinction(1,1,1))
   End If
   SW_diag(j_sw)%cloud_extinction = 0.0
!
   If (SW_diag(j_sw)%L_cloud_weight_extinction) Then
      Allocate(SW_diag(j_sw)%cloud_weight_extinction(row_length,    &
                                       rows,cloud_levels))
   Else
      Allocate(SW_diag(j_sw)%cloud_weight_extinction(1,1,1))
   End If
   SW_diag(j_sw)%cloud_weight_extinction = 0.0
!
   If (SW_diag(j_sw)%L_ls_cloud_extinction) Then
      Allocate(SW_diag(j_sw)%ls_cloud_extinction(row_length,        &
                                      rows, cloud_levels))
   Else
      Allocate(SW_diag(j_sw)%ls_cloud_extinction(1,1,1))
   End If
   SW_diag(j_sw)%ls_cloud_extinction = 0.0
!
   If (SW_diag(j_sw)%L_ls_cloud_weight_extinction) Then
       Allocate(SW_diag(j_sw)%ls_cloud_weight_extinction(           &
                                row_length,rows,cloud_levels))
   Else
       Allocate(SW_diag(j_sw)%ls_cloud_weight_extinction(1,1,1))
   End If
   SW_diag(j_sw)%ls_cloud_weight_extinction = 0.0
!
   If (SW_diag(j_sw)%L_cnv_cloud_extinction) Then
       Allocate(SW_diag(j_sw)%cnv_cloud_extinction(row_length,       &
                                           rows,cloud_levels))
   Else
       Allocate(SW_diag(j_sw)%cnv_cloud_extinction(1,1,1))
    End If
    SW_diag(j_sw)%cnv_cloud_extinction = 0.0
!
    If (SW_diag(j_sw)%L_cnv_cloud_weight_extinction) Then
       Allocate(SW_diag(j_sw)%cnv_cloud_weight_extinction(           &
                               row_length,rows, cloud_levels))
    Else
       Allocate(SW_diag(j_sw)%cnv_cloud_weight_extinction(1,1,1))
    End If
    SW_diag(j_sw)%cnv_cloud_weight_extinction = 0.0


End Subroutine allocate_sw_diag
#endif
