#if defined(A70_1C) || defined(A70_1Z)
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Set the LW radiation diagnostic flags to true if necessary
!
! Method:
!
! If a radiation diagnostic has been chosen in STASH then the
! flag of the corresponding radiation diagnostic in the structure
! LW_diag is set to true.
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
SUBROUTINE set_lwdiag_logic(sf,nitems,nsects,l_moses_II, l_ctile   &
           ,j_lw, i_off)

      Use lw_diag_mod
      Implicit None

! arguments with intent in

      integer, intent(in) :: nitems       ! STASH item number
      integer, intent(in) :: nsects       ! STASH section number
      integer, intent(in) :: j_lw         ! call to LW radiation
      integer, intent(in) :: i_off        ! offset for diagnostics

      logical, intent(in) :: l_moses_II   ! Flag for MOSES
      logical, intent(in) :: l_ctile      ! Flag for tiling

      logical, intent(in) :: sf(0:nitems,0:nsects) ! STASH Flags

#if !defined(SCMA)

! Switches enabled as STASHflags: LW

! Fluxes and Heating Rate

        LW_diag(j_lw)%L_total_cloud_cover = sf(204+i_off,2)
        LW_diag(j_lw)%L_clear_olr         = sf(206+i_off,2)
        LW_diag(j_lw)%L_surface_down_flux =(sf(207+i_off,2).OR.L_MOSES_II)
        LW_diag(j_lw)%L_surf_down_clr     = sf(208+i_off,2)
        LW_diag(j_lw)%L_clear_hr          = sf(233+i_off,2)
        LW_diag(j_lw)%L_net_flux_trop     = sf(237+i_off,2)
        LW_diag(j_lw)%L_down_flux_trop    = sf(238+i_off,2)
        LW_diag(j_lw)%L_toa_radiance      = sf(297+i_off,2)

!  Isccp diagnostics

        LW_diag(j_lw)%L_isccp_weights                 = sf(269+i_off,2)
        LW_diag(j_lw)%L_isccp_cf                      = sf(270+i_off,2)
        LW_diag(j_lw)%L_isccp_cf_tau_0_to_p3          = sf(271+i_off,2)
        LW_diag(j_lw)%L_isccp_cf_tau_p3_to_1p3        = sf(272+i_off,2)
        LW_diag(j_lw)%L_isccp_cf_tau_1p3_to_3p6       = sf(273+i_off,2)
        LW_diag(j_lw)%L_isccp_cf_tau_3p6_to_9p4       = sf(274+i_off,2)
        LW_diag(j_lw)%L_isccp_cf_tau_9p4_to_23        = sf(275+i_off,2)
        LW_diag(j_lw)%L_isccp_cf_tau_23_to_60         = sf(276+i_off,2)
        LW_diag(j_lw)%L_isccp_cf_tau_ge_60            = sf(277+i_off,2)
        LW_diag(j_lw)%L_meanalbedocld                 = sf(290+i_off,2)
        LW_diag(j_lw)%L_meantaucld                    = sf(291+i_off,2)
        LW_diag(j_lw)%L_meanptop                      = sf(292+i_off,2)
        LW_diag(j_lw)%L_totalcldarea                  = sf(293+i_off,2)

! Extinction and absorptivity diagnostics

        LW_diag(j_lw)%L_cloud_absorptivity            = sf(262+i_off,2)
        LW_diag(j_lw)%L_cloud_weight_absorptivity     = sf(263+i_off,2)
        LW_diag(j_lw)%L_ls_cloud_absorptivity         =(sf(264+i_off,2) &
          .OR.(LW_diag(j_lw)%L_isccp_weights))
        LW_diag(j_lw)%L_ls_cloud_weight_absorptivity  =(sf(265+i_off,2) &
          .OR.(LW_diag(j_lw)%L_isccp_weights))
        LW_diag(j_lw)%L_cnv_cloud_absorptivity        =(sf(266+i_off,2) &
          .OR.(LW_diag(j_lw)%L_isccp_weights))
        LW_diag(j_lw)%L_cnv_cloud_weight_absorptivity =(sf(267+i_off,2) &
          .OR.(LW_diag(j_lw)%L_isccp_weights))
        LW_diag(j_lw)%L_total_cloud_on_levels         = sf(261+i_off,2)

! Grid-box mean cloud diagnostics as seen by radiation:

        LW_diag(j_lw)%L_ls_qcl_rad                    = sf(308+i_off,2)
        LW_diag(j_lw)%L_ls_qcf_rad                    = sf(309+i_off,2)
        LW_diag(j_lw)%L_cc_qcl_rad                    = sf(310+i_off,2)
        LW_diag(j_lw)%L_cc_qcf_rad                    = sf(311+i_off,2)

        if (312+i_off <= nitems) then
          LW_diag(j_lw)%L_ls_cl_rad                   = sf(312+i_off,2)
        end if
        if (313+i_off <= nitems) then
          LW_diag(j_lw)%L_ls_cf_rad                   = sf(313+i_off,2)
        end if
        if (314+i_off <= nitems) then
          LW_diag(j_lw)%L_cc_cl_rad                   = sf(314+i_off,2)
        end if
        if (315+i_off <= nitems) then
          LW_diag(j_lw)%L_cc_cf_rad                   = sf(315+i_off,2)
        end if

! Aerosol optical depth diagnostics
        LW_diag(j_lw)%L_aod_sulphate                  = sf(284+i_off,2)
        LW_diag(j_lw)%L_aod_dust                      = sf(285+i_off,2)
        LW_diag(j_lw)%L_aod_seasalt                   = sf(286+i_off,2)
        LW_diag(j_lw)%L_aod_soot                      = sf(287+i_off,2)
        LW_diag(j_lw)%L_aod_biomass                   = sf(288+i_off,2)
        LW_diag(j_lw)%L_aod_biogenic                  = sf(289+i_off,2)
        LW_diag(j_lw)%L_aod_ocff                      = sf(295+i_off,2)
        LW_diag(j_lw)%L_aod_delta                     = sf(296+i_off,2)

#else

        LW_diag(j_lw)%L_total_cloud_cover      = .true.
        LW_diag(j_lw)%L_clear_olr              = .true.
        LW_diag(j_lw)%L_surface_down_flux      = .true.
        LW_diag(j_lw)%L_surf_down_clr          = .true.
        LW_diag(j_lw)%L_clear_hr               = .true.
        LW_diag(j_lw)%L_ls_qcl_rad             = .true.
        LW_diag(j_lw)%L_ls_qcf_rad             = .true.
        LW_diag(j_lw)%L_cc_qcl_rad             = .true.
        LW_diag(j_lw)%L_cc_qcf_rad             = .true.
        LW_diag(j_lw)%L_ls_cl_rad              = .true.
        LW_diag(j_lw)%L_ls_cf_rad              = .true.
        LW_diag(j_lw)%L_cc_cl_rad              = .true.
        LW_diag(j_lw)%L_cc_cf_rad              = .true.

#endif

end subroutine set_lwdiag_logic
#endif
#endif
