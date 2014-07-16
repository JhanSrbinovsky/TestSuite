#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Deallocate the LW diagnostics.
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
SUBROUTINE deallocate_lwdiag(j_lw)

  Use lw_diag_mod
  Implicit None

  Integer, Intent(in) :: j_lw              ! call to LW radiation

!
! Deallocate LW Diagnostics
!
  If (Associated(LW_diag(j_lw)%total_cloud_cover)) &
    Deallocate(LW_diag(j_lw)%total_cloud_cover)
  If (Associated(LW_diag(j_lw)%clear_olr)) &
    Deallocate(LW_diag(j_lw)%clear_olr)
  If (Associated(LW_diag(j_lw)%surface_down_flux)) &
    Deallocate(LW_diag(j_lw)%surface_down_flux)
  If (Associated(LW_diag(j_lw)%surf_down_clr)) &
    Deallocate(LW_diag(j_lw)%surf_down_clr)
  If (Associated(LW_diag(j_lw)%clear_hr)) &
    Deallocate(LW_diag(j_lw)%clear_hr)
  If (Associated(LW_diag(j_lw)%net_flux_trop)) &
    Deallocate(LW_diag(j_lw)%net_flux_trop)
  If (Associated(LW_diag(j_lw)%down_flux_trop)) &
    Deallocate(LW_diag(j_lw)%down_flux_trop)
  If (Associated(LW_diag(j_lw)%total_cloud_on_levels)) &
    Deallocate(LW_diag(j_lw)%total_cloud_on_levels)
    
! Radiance
    
  If (Associated(LW_diag(j_lw)%toa_radiance)) &
    Deallocate(LW_diag(j_lw)%toa_radiance)
  If (Associated(LW_diag(j_lw)%cloud_absorptivity)) &
    Deallocate(LW_diag(j_lw)%cloud_absorptivity)
  If (Associated(LW_diag(j_lw)%cloud_weight_absorptivity)) &
    Deallocate(LW_diag(j_lw)%cloud_weight_absorptivity)
  If (Associated(LW_diag(j_lw)%ls_cloud_absorptivity)) &
    Deallocate(LW_diag(j_lw)%ls_cloud_absorptivity)
  If (Associated(LW_diag(j_lw)%ls_cloud_weight_absorptivity)) &
    Deallocate(LW_diag(j_lw)%ls_cloud_weight_absorptivity)
  If (Associated(LW_diag(j_lw)%cnv_cloud_absorptivity)) &
    Deallocate(LW_diag(j_lw)%cnv_cloud_absorptivity)
  If (Associated(LW_diag(j_lw)%cnv_cloud_weight_absorptivity)) &
    Deallocate(LW_diag(j_lw)%cnv_cloud_weight_absorptivity)
    
! Deallocate the aerosol optical depth diagnostics
    
  If (Associated(LW_diag(j_lw)%aod_sulphate)) &
    Deallocate(LW_diag(j_lw)%aod_sulphate)
  If (Associated(LW_diag(j_lw)%aod_dust)) &
    Deallocate(LW_diag(j_lw)%aod_dust)
  If (Associated(LW_diag(j_lw)%aod_seasalt)) &
    Deallocate(LW_diag(j_lw)%aod_seasalt)
  If (Associated(LW_diag(j_lw)%aod_soot)) &
    Deallocate(LW_diag(j_lw)%aod_soot)
  If (Associated(LW_diag(j_lw)%aod_biomass)) &
    Deallocate(LW_diag(j_lw)%aod_biomass)
  If (Associated(LW_diag(j_lw)%aod_biogenic)) &
    Deallocate(LW_diag(j_lw)%aod_biogenic)
  If (Associated(LW_diag(j_lw)%aod_ocff)) &
    Deallocate(LW_diag(j_lw)%aod_ocff)
  If (Associated(LW_diag(j_lw)%aod_delta)) &
    Deallocate(LW_diag(j_lw)%aod_delta)
    
! Deallocate the grid-box mean cloud diagnostics
    
  If (Associated(LW_diag(j_lw)%ls_qcl_rad)) &
    Deallocate(LW_diag(j_lw)%ls_qcl_rad)
  If (Associated(LW_diag(j_lw)%ls_qcf_rad)) &
    Deallocate(LW_diag(j_lw)%ls_qcf_rad)
  If (Associated(LW_diag(j_lw)%cc_qcl_rad)) &
    Deallocate(LW_diag(j_lw)%cc_qcl_rad)
  If (Associated(LW_diag(j_lw)%cc_qcf_rad)) &
    Deallocate(LW_diag(j_lw)%cc_qcf_rad)
  If (Associated(LW_diag(j_lw)%ls_cl_rad)) &
    Deallocate(LW_diag(j_lw)%ls_cl_rad)
  If (Associated(LW_diag(j_lw)%ls_cf_rad)) &
    Deallocate(LW_diag(j_lw)%ls_cf_rad)
  If (Associated(LW_diag(j_lw)%cc_cl_rad)) &
    Deallocate(LW_diag(j_lw)%cc_cl_rad)
  If (Associated(LW_diag(j_lw)%cc_cf_rad)) &
    Deallocate(LW_diag(j_lw)%cc_cf_rad)
    
! Deallocate the isccp diagnostics
    
  If (Associated(LW_diag(j_lw)%isccp_weights)) &
    Deallocate(LW_diag(j_lw)%isccp_weights)
  If (Associated(LW_diag(j_lw)%isccp_cf)) &
    Deallocate(LW_diag(j_lw)%isccp_cf)
  If (Associated(LW_diag(j_lw)%isccp_cf_tau_0_to_p3)) &
    Deallocate(LW_diag(j_lw)%isccp_cf_tau_0_to_p3)
  If (Associated(LW_diag(j_lw)%isccp_cf_tau_p3_to_1p3)) &
    Deallocate(LW_diag(j_lw)%isccp_cf_tau_p3_to_1p3)
  If (Associated(LW_diag(j_lw)%isccp_cf_tau_1p3_to_3p6)) &
    Deallocate(LW_diag(j_lw)%isccp_cf_tau_1p3_to_3p6)
  If (Associated(LW_diag(j_lw)%isccp_cf_tau_3p6_to_9p4)) &
    Deallocate(LW_diag(j_lw)%isccp_cf_tau_3p6_to_9p4)
  If (Associated(LW_diag(j_lw)%isccp_cf_tau_9p4_to_23)) &
    Deallocate(LW_diag(j_lw)%isccp_cf_tau_9p4_to_23)
  If (Associated(LW_diag(j_lw)%isccp_cf_tau_23_to_60)) &
    Deallocate(LW_diag(j_lw)%isccp_cf_tau_23_to_60)
  If (Associated(LW_diag(j_lw)%isccp_cf_tau_ge_60)) &
    Deallocate(LW_diag(j_lw)%isccp_cf_tau_ge_60)
  If (Associated(LW_diag(j_lw)%meanalbedocld)) &
    Deallocate(LW_diag(j_lw)%meanalbedocld)
  If (Associated(LW_diag(j_lw)%meantaucld)) &
    Deallocate(LW_diag(j_lw)%meantaucld)
  If (Associated(LW_diag(j_lw)%meanptop)) &
    Deallocate(LW_diag(j_lw)%meanptop)
  If (Associated(LW_diag(j_lw)%totalcldarea)) &
    Deallocate(LW_diag(j_lw)%totalcldarea)

end subroutine deallocate_lwdiag
#endif
