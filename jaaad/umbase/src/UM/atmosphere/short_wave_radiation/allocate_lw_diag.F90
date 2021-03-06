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
SUBROUTINE allocate_lw_diag(row_length, rows, model_levels,    &
                            cloud_levels, j_lw)


! Modules

   Use spec_sw_lw
   Use lw_diag_mod

   Implicit None

   Integer, Intent(in) :: row_length        ! Length of rows
   Integer, Intent(in) :: rows              ! Number of rows
   Integer, Intent(in) :: model_levels      ! Number of model levels
   Integer, Intent(in) :: cloud_levels      ! Number of cloud levels
   Integer, Intent(in) :: j_lw              ! call to LW radiation

   If (LW_diag(j_lw)%L_total_cloud_cover) then
      allocate(LW_diag(j_lw)%total_cloud_cover(row_length, rows))
      LW_diag(j_lw)%total_cloud_cover(:,:)=0.0
   else
      allocate(LW_diag(j_lw)%total_cloud_cover(1,1))
      LW_diag(j_lw)%total_cloud_cover(1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_clear_olr) then
      allocate(LW_diag(j_lw)%clear_olr(row_length, rows))
      LW_diag(j_lw)%clear_olr(:,:)=0.0
   else
      allocate(LW_diag(j_lw)%clear_olr(1,1))
      LW_diag(j_lw)%clear_olr(1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_surface_down_flux) then
      allocate(LW_diag(j_lw)%surface_down_flux(row_length, rows))
      LW_diag(j_lw)%surface_down_flux(:,:)=0.0
   else
      allocate(LW_diag(j_lw)%surface_down_flux(1,1))
      LW_diag(j_lw)%surface_down_flux(1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_surf_down_clr) then
      allocate(LW_diag(j_lw)%surf_down_clr(row_length, rows))
      LW_diag(j_lw)%surf_down_clr(:,:)=0.0
   else
      allocate(LW_diag(j_lw)%surf_down_clr(1,1))
      LW_diag(j_lw)%surf_down_clr(1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_clear_hr) then
      allocate(LW_diag(j_lw)%clear_hr(row_length, rows, model_levels))
      LW_diag(j_lw)%clear_hr(:,:,:)=0.0
   else
      allocate(LW_diag(j_lw)%clear_hr(1,1,1))
      LW_diag(j_lw)%clear_hr(1,1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_net_flux_trop) then
      allocate(LW_diag(j_lw)%net_flux_trop(row_length, rows))
      LW_diag(j_lw)%net_flux_trop(:,:)=0.0
   else
      allocate(LW_diag(j_lw)%net_flux_trop(1,1))
      LW_diag(j_lw)%net_flux_trop(1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_down_flux_trop) then
      allocate(LW_diag(j_lw)%down_flux_trop(row_length, rows))
      LW_diag(j_lw)%down_flux_trop(:,:)=0.0
   else
      allocate(LW_diag(j_lw)%down_flux_trop(1,1))
      LW_diag(j_lw)%down_flux_trop(1,1)=0.0
   endif
!
!  Radiance
!
   If (LW_diag(j_lw)%L_toa_radiance) then
      allocate(LW_diag(j_lw)%toa_radiance(row_length, rows,            &
                                         lw_spectrum(j_lw)%n_band))
      LW_diag(j_lw)%toa_radiance(:,:,:) = 0.0
   else
      allocate(LW_diag(j_lw)%toa_radiance(1,1,1))
      LW_diag(j_lw)%toa_radiance(1,1,1) = 0.0
   End if
!
! Absorptivity diagnostics:
!
   If (LW_diag(j_lw)%L_cloud_absorptivity) then
      allocate(LW_diag(j_lw)%cloud_absorptivity(row_length,            &
                                         rows, cloud_levels))
      LW_diag(j_lw)%cloud_absorptivity(:,:,:)=0.0
   else
      allocate(LW_diag(j_lw)%cloud_absorptivity(1,1,1))
      LW_diag(j_lw)%cloud_absorptivity(1,1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_cloud_weight_absorptivity) then
      allocate(LW_diag(j_lw)%cloud_weight_absorptivity(row_length,     &
                                         rows, cloud_levels))
      LW_diag(j_lw)%cloud_weight_absorptivity(:,:,:)=0.0
   else
      allocate(LW_diag(j_lw)%cloud_weight_absorptivity(1,1,1))
      LW_diag(j_lw)%cloud_weight_absorptivity(1,1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_ls_cloud_absorptivity) then
      allocate(LW_diag(j_lw)%ls_cloud_absorptivity(row_length,         &
                                         rows, cloud_levels))
      LW_diag(j_lw)%ls_cloud_absorptivity(:,:,:)=0.0
   else
      allocate(LW_diag(j_lw)%ls_cloud_absorptivity(1,1,1))
      LW_diag(j_lw)%ls_cloud_absorptivity(1,1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_ls_cloud_weight_absorptivity) then
      allocate(LW_diag(j_lw)%ls_cloud_weight_absorptivity(             &
                              row_length,rows, cloud_levels))
      LW_diag(j_lw)%ls_cloud_weight_absorptivity(:,:,:)=0.0
   else
      allocate(LW_diag(j_lw)%ls_cloud_weight_absorptivity(1,1,1))
      LW_diag(j_lw)%ls_cloud_weight_absorptivity(1,1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_cnv_cloud_absorptivity) then
      allocate(LW_diag(j_lw)%cnv_cloud_absorptivity(row_length,        &
                                            rows, cloud_levels))
      LW_diag(j_lw)%cnv_cloud_absorptivity(:,:,:)=0.0
   else
      allocate(LW_diag(j_lw)%cnv_cloud_absorptivity(1,1,1))
      LW_diag(j_lw)%cnv_cloud_absorptivity(1,1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_cnv_cloud_weight_absorptivity) then
      allocate(LW_diag(j_lw)%cnv_cloud_weight_absorptivity(            &
                                 row_length,rows, cloud_levels))
      LW_diag(j_lw)%cnv_cloud_weight_absorptivity(:,:,:)=0.0
   else
      allocate(LW_diag(j_lw)%cnv_cloud_weight_absorptivity(1,1,1))
      LW_diag(j_lw)%cnv_cloud_weight_absorptivity(1,1,1)=0.0
   endif

   If (LW_diag(j_lw)%L_total_cloud_on_levels) then
      allocate(LW_diag(j_lw)%total_cloud_on_levels(row_length,         &
                                            rows, cloud_levels))
      LW_diag(j_lw)%total_cloud_on_levels(:,:,:)=0.0
   else
      allocate(LW_diag(j_lw)%total_cloud_on_levels(1,1,1))
      LW_diag(j_lw)%total_cloud_on_levels(1,1,1)=0.0
   endif

! Aerosol optical depth diagnostics:

   If (LW_diag(j_lw)%L_aod_sulphate) then
      allocate(LW_diag(j_lw)%aod_sulphate(row_length, rows,            &
                                       lw_spectrum(j_lw)%N_AOD_WAVEL))
      LW_diag(j_lw)%aod_sulphate(:,:,:) = 0.0
   else
      allocate(LW_diag(j_lw)%aod_sulphate(1,1,1))
   Endif

   If (LW_diag(j_lw)%L_aod_dust) then
     allocate(LW_diag(j_lw)%aod_dust(row_length, rows,                 &
                                       lw_spectrum(j_lw)%N_AOD_WAVEL))
     LW_diag(j_lw)%aod_dust(:,:,:) = 0.0
   else
     allocate(LW_diag(j_lw)%aod_dust(1,1,1))
   Endif

   If (LW_diag(j_lw)%L_aod_seasalt) then
     allocate(LW_diag(j_lw)%aod_seasalt(row_length, rows,              &
                                       lw_spectrum(j_lw)%N_AOD_WAVEL))
     LW_diag(j_lw)%aod_seasalt(:,:,:) = 0.0
   else
     allocate(LW_diag(j_lw)%aod_seasalt(1,1,1))
   Endif

   If (LW_diag(j_lw)%L_aod_soot) then
     allocate(LW_diag(j_lw)%aod_soot(row_length, rows,                 &
                                       lw_spectrum(j_lw)%N_AOD_WAVEL))
     LW_diag(j_lw)%aod_soot(:,:,:) = 0.0
   else
     allocate(LW_diag(j_lw)%aod_soot(1,1,1))
   Endif

   If (LW_diag(j_lw)%L_aod_biomass) then
     allocate(LW_diag(j_lw)%aod_biomass(row_length, rows,              &
                                       lw_spectrum(j_lw)%N_AOD_WAVEL))
     LW_diag(j_lw)%aod_biomass(:,:,:) = 0.0
   else
     allocate(LW_diag(j_lw)%aod_biomass(1,1,1))
   Endif

   If (LW_diag(j_lw)%L_aod_biogenic) then
     allocate(LW_diag(j_lw)%aod_biogenic(row_length, rows,             &
                                       lw_spectrum(j_lw)%N_AOD_WAVEL))
     LW_diag(j_lw)%aod_biogenic(:,:,:) = 0.0
   else
     allocate(LW_diag(j_lw)%aod_biogenic(1,1,1))
   Endif

   If (LW_diag(j_lw)%L_aod_ocff) then
     allocate(LW_diag(j_lw)%aod_ocff(row_length, rows,                 &
                                       lw_spectrum(j_lw)%N_AOD_WAVEL))
     LW_diag(j_lw)%aod_ocff(:,:,:) = 0.0
   else
     allocate(LW_diag(j_lw)%aod_ocff(1,1,1))
   Endif

   If (LW_diag(j_lw)%L_aod_delta) then
     allocate(LW_diag(j_lw)%aod_delta(row_length, rows,                &
                                       lw_spectrum(j_lw)%N_AOD_WAVEL))
     LW_diag(j_lw)%aod_delta(:,:,:) = 0.0
   else
     allocate(LW_diag(j_lw)%aod_delta(1,1,1))
   Endif

!
! Grid-box mean cloud diagnostics as seen by radiation:
!
   If (LW_diag(j_lw)%L_ls_qcl_rad) then
     allocate(LW_diag(j_lw)%ls_qcl_rad(row_length, rows, model_levels))
   else
     allocate(LW_diag(j_lw)%ls_qcl_rad(1,1,1))
   endif
   LW_diag(j_lw)%ls_qcl_rad(:,:,:)=0.0

   If (LW_diag(j_lw)%L_ls_qcf_rad) then
     allocate(LW_diag(j_lw)%ls_qcf_rad(row_length, rows, model_levels))
   else
     allocate(LW_diag(j_lw)%ls_qcf_rad(1,1,1))
   endif
   LW_diag(j_lw)%ls_qcf_rad(:,:,:)=0.0

   If (LW_diag(j_lw)%L_cc_qcl_rad) then
     allocate(LW_diag(j_lw)%cc_qcl_rad(row_length, rows, model_levels))
   else
     allocate(LW_diag(j_lw)%cc_qcl_rad(1,1,1))
   endif
   LW_diag(j_lw)%cc_qcl_rad(:,:,:)=0.0

   If (LW_diag(j_lw)%L_cc_qcf_rad) then
     allocate(LW_diag(j_lw)%cc_qcf_rad(row_length, rows, model_levels))
   else
     allocate(LW_diag(j_lw)%cc_qcf_rad(1,1,1))
   endif
   LW_diag(j_lw)%cc_qcf_rad(:,:,:)=0.0

   If (LW_diag(j_lw)%L_ls_cl_rad) then
     allocate(LW_diag(j_lw)%ls_cl_rad(row_length, rows, model_levels))
   else
     allocate(LW_diag(j_lw)%ls_cl_rad(1,1,1))
   endif
   LW_diag(j_lw)%ls_cl_rad(:,:,:)=0.0

   If (LW_diag(j_lw)%L_ls_cf_rad) then
     allocate(LW_diag(j_lw)%ls_cf_rad(row_length, rows, model_levels))
   else
     allocate(LW_diag(j_lw)%ls_cf_rad(1,1,1))
   endif
   LW_diag(j_lw)%ls_cf_rad(:,:,:)=0.0

   If (LW_diag(j_lw)%L_cc_cl_rad) then
     allocate(LW_diag(j_lw)%cc_cl_rad(row_length, rows, model_levels))
   else
     allocate(LW_diag(j_lw)%cc_cl_rad(1,1,1))
   endif
   LW_diag(j_lw)%cc_cl_rad(:,:,:)=0.0

   If (LW_diag(j_lw)%L_cc_cf_rad) then
     allocate(LW_diag(j_lw)%cc_cf_rad(row_length, rows, model_levels))
   else
     allocate(LW_diag(j_lw)%cc_cf_rad(1,1,1))
   endif
   LW_diag(j_lw)%cc_cf_rad(:,:,:)=0.0


! ISCCP diagnostics

   If (LW_diag(j_lw)%L_isccp_weights) then
      allocate(LW_diag(j_lw)%isccp_weights(row_length,rows))
      LW_diag(j_lw)%isccp_weights(:,:) = 0.0
   else
      allocate(LW_diag(j_lw)%isccp_weights(1,1))
      LW_diag(j_lw)%isccp_weights(1,1) = 0.0
   End if

   If (LW_diag(j_lw)%L_isccp_cf) then
      allocate(LW_diag(j_lw)%isccp_cf(row_length,rows,7))
      LW_diag(j_lw)%isccp_cf(:,:,:) = 0.0
   else
      allocate(LW_diag(j_lw)%isccp_cf(1,1,1))
      LW_diag(j_lw)%isccp_cf(1,1,1) = 0.0
   End if

   If (LW_diag(j_lw)%L_isccp_cf_tau_0_to_p3) then
      allocate(LW_diag(j_lw)%isccp_cf_tau_0_to_p3(row_length,rows,7))
      LW_diag(j_lw)%isccp_cf_tau_0_to_p3(:,:,:) = 0.0
   else
      allocate(LW_diag(j_lw)%isccp_cf_tau_0_to_p3(1,1,1))
      LW_diag(j_lw)%isccp_cf_tau_0_to_p3(1,1,1) = 0.0
   End if

   If (LW_diag(j_lw)%L_isccp_cf_tau_p3_to_1p3) then
      allocate(LW_diag(j_lw)%isccp_cf_tau_p3_to_1p3(row_length,rows,7))
      LW_diag(j_lw)%isccp_cf_tau_p3_to_1p3(:,:,:) = 0.0
   else
      allocate(LW_diag(j_lw)%isccp_cf_tau_p3_to_1p3(1,1,1))
      LW_diag(j_lw)%isccp_cf_tau_p3_to_1p3(1,1,1) = 0.0
   End if

   If (LW_diag(j_lw)%L_isccp_cf_tau_1p3_to_3p6) then
      allocate(LW_diag(j_lw)%isccp_cf_tau_1p3_to_3p6(row_length,rows,7))
      LW_diag(j_lw)%isccp_cf_tau_1p3_to_3p6(:,:,:) = 0.0
   else
      allocate(LW_diag(j_lw)%isccp_cf_tau_1p3_to_3p6(1,1,1))
      LW_diag(j_lw)%isccp_cf_tau_1p3_to_3p6(1,1,1) = 0.0
   End if

  If (LW_diag(j_lw)%L_isccp_cf_tau_3p6_to_9p4) then
      allocate(LW_diag(j_lw)%isccp_cf_tau_3p6_to_9p4(row_length,rows,7))
      LW_diag(j_lw)%isccp_cf_tau_3p6_to_9p4(:,:,:) = 0.0
  else
      allocate(LW_diag(j_lw)%isccp_cf_tau_3p6_to_9p4(1,1,1))
      LW_diag(j_lw)%isccp_cf_tau_3p6_to_9p4(1,1,1) = 0.0
  End if

  If (LW_diag(j_lw)%L_isccp_cf_tau_9p4_to_23) then
     allocate(LW_diag(j_lw)%isccp_cf_tau_9p4_to_23(row_length,rows,7))
     LW_diag(j_lw)%isccp_cf_tau_9p4_to_23(:,:,:) = 0.0
  else
     allocate(LW_diag(j_lw)%isccp_cf_tau_9p4_to_23(1,1,1))
     LW_diag(j_lw)%isccp_cf_tau_9p4_to_23(1,1,1) = 0.0
  End if

  If (LW_diag(j_lw)%L_isccp_cf_tau_23_to_60) then
     allocate(LW_diag(j_lw)%isccp_cf_tau_23_to_60(row_length,rows,7))
     LW_diag(j_lw)%isccp_cf_tau_23_to_60(:,:,:) = 0.0
  else
     allocate(LW_diag(j_lw)%isccp_cf_tau_23_to_60(1,1,1))
     LW_diag(j_lw)%isccp_cf_tau_23_to_60(1,1,1) = 0.0
  End if

  If (LW_diag(j_lw)%L_isccp_cf_tau_ge_60) then
     allocate(LW_diag(j_lw)%isccp_cf_tau_ge_60(row_length,rows,7))
     LW_diag(j_lw)%isccp_cf_tau_ge_60(:,:,:) = 0.0
  else
     allocate(LW_diag(j_lw)%isccp_cf_tau_ge_60(1,1,1))
     LW_diag(j_lw)%isccp_cf_tau_ge_60(1,1,1) = 0.0
  End if
  If (LW_diag(j_lw)%L_meanalbedocld) then
       allocate(LW_diag(j_lw)%meanalbedocld(row_length,rows))
       LW_diag(j_lw)%meanalbedocld(:,:) = 0.0
  else
       allocate(LW_diag(j_lw)%meanalbedocld(1,1))
       LW_diag(j_lw)%meanalbedocld(1,1) = 0.0
  End if
  If (LW_diag(j_lw)%L_meantaucld) then
       allocate(LW_diag(j_lw)%meantaucld(row_length,rows))
       LW_diag(j_lw)%meantaucld(:,:) = 0.0
  else
       allocate(LW_diag(j_lw)%meantaucld(1,1))
       LW_diag(j_lw)%meantaucld(1,1) = 0.0
  End if
  If (LW_diag(j_lw)%L_meanptop) then
       allocate(LW_diag(j_lw)%meanptop(row_length,rows))
       LW_diag(j_lw)%meanptop(:,:) = 0.0
  else
       allocate(LW_diag(j_lw)%meanptop(1,1))
       LW_diag(j_lw)%meanptop(1,1) = 0.0
  End if
  If (LW_diag(j_lw)%L_totalcldarea) then
       allocate(LW_diag(j_lw)%totalcldarea(row_length,rows))
       LW_diag(j_lw)%totalcldarea(:,:) = 0.0
  else
       allocate(LW_diag(j_lw)%totalcldarea(1,1))
       LW_diag(j_lw)%totalcldarea(1,1)=0.0
  End if

end subroutine allocate_lw_diag
#endif
