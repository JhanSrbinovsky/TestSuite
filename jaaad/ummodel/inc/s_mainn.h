!include file: s_mainn.h
! Description:
!   Defines additional namelist data for the SCM.
!
! Current Code Owner: R Wong
!

#if defined(ATMOS)

!====================================================================
!     Specify namelists
!====================================================================
!     INDATA     Initial data required to run the model
!     INOBSFOR   Initial data required if observational forcing
!                is chosen
!     INPROF     Inital model profile of the primary variables (UMDP
!                    No1) plus T, but QCL and QCF initialised by
!                    subroutine INITQLCF
!     RUNDATA    Data required for run
!     INMOSES    Data for initialising soil moisture for MOSES code
!--------------------------------------------------------------------

      Namelist/INDATA/soil_type, tapeyear_init,                               &
        tapemonth_init, tapeday_init, tapehour_init, tapemin_init,            &
        tapesec_init, year_init, month_init, day_init, hour_init,             &
        min_init, sec_init, lat, long, gridbox_area,                          &
        salt_dim1, salt_dim2, salt_dim3, gather

      Namelist/INOBSFOR/l_windrlx, tau_rlx, l_vertadv, t_inc, q_star,         &
        u_inc, v_inc, w_inc, ichgf, flux_h, flux_e, tstar_forcing,            &
        tau_T_eurocs, tau_q_eurocs, tau_uv_eurocs

      Namelist/INPROF/ui, vi, wi, theta, qi, smci, canopy_gbi,                &
        tstari, t_deep_soili, snodepi, z0mseai, z0m_scm, z0h_scm,             &
        ccai, iccbi, iccti, tracer, sil_orog_land,                            &
        ho2r2_orog, ice_fract, di, u_0, v_0, p_in, w_advi

      Namelist/LOGIC/ancyc, altdat, tapein, tapeout, test,                    &
        obs, stats, prindump_step, prindump_day, prindump_days,               &
        prindump_obs, land_sea_mask, land_ice_mask, soil_mask,                &
        noforce, grafdump_day, grafdump_days, geoforce, geoinit,              &
        grafdump_step, cumulus, local_time,                                   &
        l_do_t, l_do_inc_vels, l_do_inc_q, prinstat,                          &
        L_spec_z0, radcloud_fixed,                                            &
        l_flux_bc, relaxT_eurocs, relaxq_eurocs, relaxuv_eurocs

      Namelist/RUNDATA/ndayin, nminin, nsecin, ntrad, timestep,               &
        ntrad1, resdump_days, exname_in, runno_in, exname_out,                &
        runno_out, dump_step, dump_days, min_trop_level,                      &
        max_trop_level, ntml, nbdsc, ntdsc,                                   &
        so2_em, nh3_em, dms_em, soot_em, co2start,                            &
        co2end, co2rate, co2, change_clim, cort, cord, corvn, corw,           &
        soot, cclwp, orog, ozone, zh, albsoil,                                &
        dolr_rts, aerosol_em, sum_eng_fluxes, sum_moist_flux,                 &
        aerosol, so2, so4_aitken, so4_accu, so4_diss,                         &
        dms, nh3, soot_new, soot_cld, soot_aged, soot_hilem,                  &
        fland_ctile, tstar_land, tstar_sea, tstar_sice, sice_alb,             &
        land_alb, co2_emits, co2flux,                                         &
    !kdcorbin, 05/10 - added tracer fluxes
        tracer_flux1,tracer_flux2,tracer_flux3,tracer_flux4,tracer_flux5,     &
        tracer_flux6,tracer_flux7,tracer_flux8,tracer_flux9,tracer_flux10,    &
        tracer_flux11,tracer_flux12,tracer_flux13,tracer_flux14,tracer_flux15,&
        tracer_flux16,tracer_flux17,tracer_flux18,tracer_flux19,tracer_flux20

      Namelist/INMOSES/smi_opt, smcli, fsmc, sth, canht, catch, snow_tile,    &
        lai, z0_tile, tstar_tile, canopy, frac_typ, frac_disturb, infil_tile, &
        rgrain, cs, gs, g_leaf_acc, g_leaf_phen_acc, npp_ft_acc,              &
        resp_w_ft_acc, resp_s_acc, lw_down

      Namelist/INGEOFOR/ug, vg

      Namelist/RADCLOUD/ cca_rad, iccb_rad, icct_rad,                         &
        layer_cloud_rad, qcl_rad, qcf_rad, ccwpin_rad

      Namelist/PHYSWITCH/ conv_mode

#endif
