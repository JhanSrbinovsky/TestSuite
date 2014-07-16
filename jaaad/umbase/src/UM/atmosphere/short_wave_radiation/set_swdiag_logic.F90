#if defined(A70_1C) || defined(A70_1Z)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Set the SW diagnostic flags to true if necessary
!
! Method:
!
! If a radiation diagnostic has been chosen in STASH then the
! flag of the corresponding radiation diagnostic in the structure
! SW_diag is set to true.
!
! Current Code Owner: J.-C. Thelen
!
! History:
! Version     Date        Comment
! 6.2         19/01/06    Original Code (J.-C. Thelen)
!
! Description of Code:
!    FORTRAN 77/90  with extensions listed in documentation.
!
!-----------------------------------------------------------------------
!
SUBROUTINE set_swdiag_logic(sf,nitems,nsects,l_moses_II, l_ctile   &
           , L_flux_below_690nm_surf,j_sw, i_off)

      Use lw_diag_mod
      Use sw_diag_mod
      Implicit None

! arguments with intent in

      integer, intent(in) :: nitems    ! item number
      integer, intent(in) :: nsects    ! section number
      integer, intent(in) :: j_sw      ! call to SW radiation
      integer, intent(in) :: i_off     ! offset for diagnostics

      logical, intent(in) :: l_moses_II    ! Flag for Moses
      logical, intent(in) :: l_ctile       ! Flag for Tiling

      logical, intent(in) ::  L_flux_below_690nm_surf
!        flag for calculating surface flux in specified range.

      logical, intent(in) :: sf(0:nitems,0:nsects)
!        STASH Flags

#if !defined(SCMA)

! Switches enabled as STASHflags: SW

! Fluxes and Heating Rates

        SW_diag(j_sw)%L_solar_out_toa     = sf(208+i_off,1)
        SW_diag(j_sw)%L_solar_out_clear   = sf(209+i_off,1)
        SW_diag(j_sw)%L_surface_down_flux =(sf(235+i_off,1).OR.L_MOSES_II)
        SW_diag(j_sw)%L_surf_down_clr     = sf(210+i_off,1)
        SW_diag(j_sw)%L_surf_up_clr       = sf(211+i_off,1)
        SW_diag(j_sw)%L_clear_hr          = sf(233+i_off,1)
        SW_diag(j_sw)%L_net_flux_trop     = sf(237+i_off,1)
        SW_diag(j_sw)%L_up_flux_trop      = sf(238+i_off,1)
!
! Radiances
!
        SW_diag(j_sw)%L_toa_radiance      = sf(297+i_off,1)
!
! Diffuse & Direct Flux
!
        SW_diag(j_sw)%L_flux_direct       = sf(230+i_off,1)
        SW_diag(j_sw)%L_flux_diffuse      = sf(231+i_off,1)

! UV-Fluxes

        SW_diag(j_sw)%L_uvflux_direct     = sf(212+i_off,1)
        SW_diag(j_sw)%L_uvflux_up         = sf(213+i_off,1)
        SW_diag(j_sw)%L_uvflux_net        = sf(214+i_off,1)

! Microphysical diagnostics

        SW_diag(j_sw)%re_strat_flag       = sf(221+i_off,1)
        SW_diag(j_sw)%wgt_strat_flag      = sf(223+i_off,1)
        SW_diag(j_sw)%lwp_strat_flag      = sf(224+i_off,1)
        SW_diag(j_sw)%re_conv_flag        = sf(225+i_off,1)
        SW_diag(j_sw)%wgt_conv_flag       = sf(226+i_off,1)
        SW_diag(j_sw)%ntot_diag_flag      = sf(241+i_off,1)
        SW_diag(j_sw)%strat_lwc_diag_flag = sf(242+i_off,1)
        SW_diag(j_sw)%so4_ccn_diag_flag   = sf(243+i_off,1)
        SW_diag(j_sw)%cond_samp_wgt_flag  = sf(244+i_off,1)
        SW_diag(j_sw)%weighted_re_flag    = sf(245+i_off,1)
        SW_diag(j_sw)%sum_weight_re_flag  = sf(246+i_off,1)
        SW_diag(j_sw)%wgtd_warm_re_flag   = sf(254+i_off,1)
        SW_diag(j_sw)%sum_wgt_warm_re_flag= sf(255+i_off,1)
        SW_diag(j_sw)%seasalt_film_flag   = sf(247+i_off,1)
        SW_diag(j_sw)%seasalt_jet_flag    = sf(248+i_off,1)
        SW_diag(j_sw)%Nc_diag_flag        = sf(280+i_off,1)
        SW_diag(j_sw)%Nc_weight_flag      = sf(281+i_off,1)
        SW_diag(j_sw)%L_direct_par        = sf(291+i_off,1)

! Diagnostics for MOSES II

        SW_diag(j_sw)%L_FlxSolBelow690nmSurf = (sf(259+i_off,1) .OR.   &
         (L_flux_below_690nm_surf .AND. L_ctile) )
        SW_diag(j_sw)%L_FlxSeaBelow690nmSurf = (sf(260+i_off,1) .OR.   &
         (L_flux_below_690nm_surf .AND. L_ctile) )

! Diagnostics for orography correction

        SW_diag(j_sw)%L_sol_bearing = sf(292+i_off,1)
        SW_diag(j_sw)%L_orog_corr   = sf(295+i_off,1)

! Extinction and absorptivity diagnostics

        SW_diag(j_sw)%L_cloud_extinction           = sf(262+i_off,1)
        SW_diag(j_sw)%L_cloud_weight_extinction    = sf(263+i_off,1)
        SW_diag(j_sw)%L_ls_cloud_extinction        = (sf(264+i_off,1)  &
           .OR.(LW_diag(j_sw)%L_isccp_weights))
        SW_diag(j_sw)%L_ls_cloud_weight_extinction = (sf(265+i_off,1)  &
           .OR.(LW_diag(j_sw)%L_isccp_weights))
        SW_diag(j_sw)%L_cnv_cloud_extinction       = (sf(266+i_off,1)  &
           .OR.(LW_diag(j_sw)%L_isccp_weights))
        SW_diag(j_sw)%L_cnv_cloud_weight_extinction= (sf(267+i_off,1)  &
           .OR.(LW_diag(j_sw)%L_isccp_weights))
#else
!
!     SW diagnostics: Only the most common diagnostics are obtained
!     by default: this can be changed using a modset.
!
        SW_diag(j_sw)%L_solar_out_toa          = .true.
        SW_diag(j_sw)%L_solar_out_clear        = .true.
        SW_diag(j_sw)%L_surface_down_flux      = .true.
        SW_diag(j_sw)%L_surf_down_clr          = .true.
        SW_diag(j_sw)%L_surf_up_clr            = .true.
        SW_diag(j_sw)%L_clear_hr               = .true.
        SW_diag(j_sw)%re_strat_flag            = .true.
        SW_diag(j_sw)%wgt_strat_flag           = .true.
        SW_diag(j_sw)%lwp_strat_flag           = .true.
        SW_diag(j_sw)%re_conv_flag             = .true.
        SW_diag(j_sw)%wgt_conv_flag            = .true.
        SW_diag(j_sw)%L_FlxSolBelow690nmSurf   =  &
           (L_flux_below_690nm_surf .AND. L_ctile)
        SW_diag(j_sw)%L_FlxSeaBelow690nmSurf   =  &
           (L_flux_below_690nm_surf .AND. L_ctile)
#endif

end subroutine set_swdiag_logic
#endif
#endif
