#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Initialise all the SW Diagnostics to False.
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
SUBROUTINE init_swdiag_logic(j_sw)

      Use sw_diag_mod
      Implicit None

! arguments with intent in

      integer, intent(in) :: j_sw      ! call to SW radiation

!
! Switches enabled as STASHflags: SW
!
! Fluxes and Heating Rates
!

        SW_diag(j_sw)%L_solar_out_toa     = .false.
        SW_diag(j_sw)%L_solar_out_clear   = .false.
        SW_diag(j_sw)%L_surface_down_flux = .false.
        SW_diag(j_sw)%L_surf_down_clr     = .false.
        SW_diag(j_sw)%L_surf_up_clr       = .false.
        SW_diag(j_sw)%L_clear_hr          = .false.
        SW_diag(j_sw)%L_net_flux_trop     = .false.
        SW_diag(j_sw)%L_up_flux_trop      = .false.
        SW_diag(j_sw)%L_flux_direct       = .false.
        SW_diag(j_sw)%L_flux_diffuse      = .false.

! UV-Fluxes

        SW_diag(j_sw)%L_uvflux_direct     = .false.
        SW_diag(j_sw)%L_uvflux_up         = .false.
        SW_diag(j_sw)%L_uvflux_net        = .false.

! Radiance

        SW_diag(j_sw)%L_toa_radiance      = .false.

! Microphysical diagnostics

        SW_diag(j_sw)%re_strat_flag       = .false.
        SW_diag(j_sw)%wgt_strat_flag      = .false.
        SW_diag(j_sw)%lwp_strat_flag      = .false.
        SW_diag(j_sw)%re_conv_flag        = .false.
        SW_diag(j_sw)%wgt_conv_flag       = .false.
        SW_diag(j_sw)%ntot_diag_flag      = .false.
        SW_diag(j_sw)%strat_lwc_diag_flag = .false.
        SW_diag(j_sw)%so4_ccn_diag_flag   = .false.
        SW_diag(j_sw)%cond_samp_wgt_flag  = .false.
        SW_diag(j_sw)%weighted_re_flag    = .false.
        SW_diag(j_sw)%sum_weight_re_flag  = .false.
        SW_diag(j_sw)%wgtd_warm_re_flag   = .false.
        SW_diag(j_sw)%sum_wgt_warm_re_flag= .false.
        SW_diag(j_sw)%seasalt_film_flag   = .false.
        SW_diag(j_sw)%seasalt_jet_flag    = .false.
        SW_diag(j_sw)%Nc_diag_flag        = .false.
        SW_diag(j_sw)%Nc_weight_flag      = .false.
        SW_diag(j_sw)%L_direct_par        = .false.

! Diagnostics for MOSES II

        SW_diag(j_sw)%L_FlxSolBelow690nmSurf = .false.
        SW_diag(j_sw)%L_FlxSeaBelow690nmSurf = .false.

! Diagnostics for orography correction

        SW_diag(j_sw)%L_sol_bearing = .false.
        SW_diag(j_sw)%L_orog_corr   = .false.

! Extinction and absorptivity diagnostics

        SW_diag(j_sw)%L_cloud_extinction           = .false.
        SW_diag(j_sw)%L_cloud_weight_extinction    = .false.
        SW_diag(j_sw)%L_ls_cloud_extinction        = .false.
        SW_diag(j_sw)%L_ls_cloud_weight_extinction = .false.
        SW_diag(j_sw)%L_cnv_cloud_extinction       = .false.
        SW_diag(j_sw)%L_cnv_cloud_weight_extinction= .false.

end subroutine init_swdiag_logic
#endif
