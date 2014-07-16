!------------------------------------------------------------------------------
! Begin cv_run_nml.h
! Defines sub-namelist &Run_Convection read in from CNTLATM control file.
!
! Changes made to this list will affect both the Full UM and the SCM
!------------------------------------------------------------------------------

        Namelist/Run_Convection/                                              &

        ! Logical switches
        l_mom,       l_fix_udfactor,       l_scv,    l4a_kterm,    l_sdxs,    &
        l_xscomp,    l_convcld_hadgem1,    l_ccw,    l_eman_dd,    l_cape,    &
        l_conv4a,    l_cloud_deep,         l_rp,     l_rp2,    l_no_dcrit,    &

        ! General scheme options/variables
        n_conv_calls,         a_convect_segments,       a_convect_seg_size,   &
        convection_option,    icvdiag,                  sh_pert_opt,          &
        dd_opt,               deep_cmt_opt,             mid_cmt_opt,          &
        termconv,             adapt,                    r_det,                &               
        tice,                 qstice,                                         &
        entcoef_min,          entcoef_max,              bl_cnv_mix,           &
        mid_cnv_pmin,         ccw_for_precip_opt,                             &

        ! Cape related options/variables
        cape_opt,             cape_bottom,              cape_top,             &
        cape_timescale,       cape_timescale_min,       cape_timescale_max,   &
        cape_ts_w,            w_cape_limit,             cape_min,             &

        ! Convective cloud options/variables
        cld_life_opt,         rad_cloud_decay_opt,      cca_min,              &
        fixed_cld_life,       ud_factor,                mparwtr,              &

        ! CCRad options options/variables
        cca2d_sh_opt,         cca_sh_knob,          ccw_sh_knob,              &
        cca2d_md_opt,         cca_md_knob,          ccw_md_knob,              &
        cca2d_dp_opt,         cca_dp_knob,          ccw_dp_knob,              &

        ! Anvil scheme options
        anvil_factor,         tower_factor,         anv_opt,                  &
 
        ! Convection type options
        iconv_shallow,        iconv_mid,            iconv_deep,               &
        iconv_congestus

! End cv_run_nml.h
!------------------------------------------------------------------------------
