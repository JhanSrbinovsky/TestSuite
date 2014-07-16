!------------------------------------------------------------------------------
! Begin cv_run_nml_data.h
! Contains default for sub-namelist &Run_Convection read in from CNTLATM
! control file. 
!
! Changes made to this list will affect both the Full UM and the SCM
! Please see include file cv_run_nml.h for variable comments and possible 
! Options
!
!------------------------------------------------------------------------------

      ! Set &Run_convection defaults
      l_fix_udfactor      = .FALSE. 
      l_convcld_hadgem1   = .FALSE. 
      l_cloud_deep        = .TRUE.  
      l_no_dcrit          = .FALSE. 
      l_mom               = .FALSE. 
      l_scv               = .FALSE. 
      l_rp                = .FALSE. 
      l_rp2               = .FALSE. 
      l_conv4a            = .FALSE. 
      l4a_kterm           = .TRUE.  
      l_eman_dd           = .FALSE. 
      l_sdxs              = .FALSE. 
      l_xscomp            = .TRUE.  
      l_ccw               = .TRUE.  
      l_cape              = .TRUE.  
      n_conv_calls        = 1 
      a_convect_segments  = 1 
      a_convect_seg_size  = -99
      convection_option   = 3 
      cld_life_opt        = cld_life_constant
      rad_cloud_decay_opt = rad_decay_off
      anv_opt             = anv_model_levels
      iconv_shallow       = 0
      iconv_congestus     = 0
      iconv_mid           = 0
      iconv_deep          = 0
      deep_cmt_opt        = 0
      mid_cmt_opt         = 0
      ccw_for_precip_opt  = 0
      icvdiag             = 1
      adapt               = 0
      cape_opt            = 0
      cape_bottom         = IMDI
      cape_top            = IMDI
      sh_pert_opt         = 0
      dd_opt              = 0
      termconv            = 0
      bl_cnv_mix          = 0
      cca2d_sh_opt        = 0
      cca2d_md_opt        = 0
      cca2d_dp_opt        = 0
      cca_sh_knob         = 1.0
      cca_md_knob         = 1.0
      cca_dp_knob         = 1.0
      ccw_sh_knob         = 1.0
      ccw_md_knob         = 1.0
      ccw_dp_knob         = 1.0
      fixed_cld_life      = 7200.0
      cca_min             = 0.02
      r_det               = 0.75
      cape_ts_w           = RMDI
      cape_min            = RMDI
      w_cape_limit        = RMDI
      mparwtr             = RMDI
      anvil_factor        = RMDI
      tower_factor        = RMDI
      ud_factor           = RMDI
      tice                = 273.15
      qstice              = 3.5E-3
      entcoef_min         = 2.75
      entcoef_max         = 4.0
      mid_cnv_pmin        = 0.0
      cape_timescale      = RMDI
      cape_timescale_min  = 1800.0
      cape_timescale_max  = 3600.0

! End cv_run_nml_data.h
!------------------------------------------------------------------------------
