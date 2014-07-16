#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ List of stashcode magic numbers

MODULE Rcf_Stashcodes_Mod

! Description:
!    Magic numbers for stashcodes used in the reconfiguration
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

! Stashcode for standard sections useful in rcf_locate.
Integer, Parameter     :: stashcode_prog_sec          =   0
Integer, Parameter     :: stashcode_bl_sec            =   3 !kdcorbin, 05/10
Integer, Parameter     :: stashcode_tracer_sec        =  33
Integer, Parameter     :: stashcode_ukca_sec          =  34

! Stashcodes used in the reconfiguration

Integer, Parameter     :: stashcode_u                 =   2
Integer, Parameter     :: stashcode_v                 =   3
Integer, Parameter     :: stashcode_theta             =   4
Integer, Parameter     :: stashcode_soil_moist        =   9
Integer, Parameter     :: stashcode_q                 =  10
Integer, Parameter     :: stashcode_qcf               =  12
Integer, Parameter     :: stashcode_cca               =  13
Integer, Parameter     :: stashcode_ccb               =  14
Integer, Parameter     :: stashcode_cct               =  15
Integer, Parameter     :: stashcode_cc_lwp            =  16
Integer, Parameter     :: stashcode_soil_temp         =  20
Integer, Parameter     :: stashcode_lcbase            =  21
Integer, Parameter     :: stashcode_mean_canopyw      =  22
Integer, Parameter     :: stashcode_mean_snow         =  23
Integer, Parameter     :: stashcode_tstar             =  24
Integer, Parameter     :: stashcode_bl_depth          =  25
Integer, Parameter     :: stashcode_z0                =  26
Integer, Parameter     :: stashcode_surf_z_curr       =  28
Integer, Parameter     :: stashcode_surf_m_curr       =  29
Integer, Parameter     :: stashcode_lsm               =  30
Integer, Parameter     :: stashcode_icefrac           =  31
Integer, Parameter     :: stashcode_icethick          =  32
Integer, Parameter     :: stashcode_orog              =  33
Integer, Parameter     :: stashcode_tstar_anom        =  39
Integer, Parameter     :: stashcode_vol_smc_wilt      =  40
Integer, Parameter     :: stashcode_vol_smc_cri       =  41
Integer, Parameter     :: stashcode_vol_smc_sat       =  43
Integer, Parameter     :: stashcode_Ksat              =  44
Integer, Parameter     :: stashcode_soil_suction      =  48
Integer, Parameter     :: stashcode_sea_ice_temp      =  49
Integer, Parameter     :: stashcode_ozone             =  60
Integer, Parameter     :: stashcode_snow_on_ice       =  95
Integer, Parameter     :: stashcode_mmr_smoke_fr      = 111
Integer, Parameter     :: stashcode_mmr_smoke_ag      = 112
Integer, Parameter     :: stashcode_mmr_smoke_cl      = 113
Integer, Parameter     :: stashcode_3d_nat_so2_em     = 121
Integer, Parameter     :: stashcode_3d_oh_conc        = 122
Integer, Parameter     :: stashcode_biom_surf_em      = 130
Integer, Parameter     :: stashcode_biom_elev_em      = 131
Integer, Parameter     :: stashcode_w                 = 150
Integer, Parameter     :: stashcode_riv_sequence      = 151
Integer, Parameter     :: stashcode_riv_direction     = 152
Integer, Parameter     :: stashcode_riv_storage       = 153
Integer, Parameter     :: stashcode_u_compnt_pert     = 202
Integer, Parameter     :: stashcode_v_compnt_pert     = 203
Integer, Parameter     :: stashcode_clapp_hb          = 207
Integer, Parameter     :: stashcode_3d_cca            = 211
Integer, Parameter     :: stashcode_3d_ccw            = 212
Integer, Parameter     :: stashcode_can_conduct       = 213
Integer, Parameter     :: stashcode_unfrozen_soil     = 214
Integer, Parameter     :: stashcode_frozen_soil       = 215
Integer, Parameter     :: stashcode_npp_pft_acc       = 224
Integer, Parameter     :: stashcode_g_lf_pft_acc      = 225
Integer, Parameter     :: stashcode_g_ph_lf_pft_acc   = 226
Integer, Parameter     :: stashcode_rsp_w_pft_acc     = 227
Integer, Parameter     :: stashcode_rsp_s_acc         = 228
Integer, Parameter     :: stashcode_can_water_tile    = 229
Integer, Parameter     :: stashcode_catch_tile        = 230
Integer, Parameter     :: stashcode_rgrain            = 231
Integer, Parameter     :: stashcode_tstar_tile        = 233
Integer, Parameter     :: stashcode_z0_tile           = 234
Integer, Parameter     :: stashcode_infil_max_tile    = 236
Integer, Parameter     :: stashcode_sw_down_tile      = 237
Integer, Parameter     :: stashcode_sw_down           = 238
Integer, Parameter     :: stashcode_lw_up_diff        = 239
Integer, Parameter     :: stashcode_snow_tile         = 240
Integer, Parameter     :: stashcode_catch_snow        = 241
Integer, Parameter     :: stashcode_snow_grnd         = 242
Integer, Parameter     :: stashcode_rho               = 253
Integer, Parameter     :: stashcode_qcl               = 254
Integer, Parameter     :: stashcode_exner             = 255
Integer, Parameter     :: stashcode_u_adv             = 256
Integer, Parameter     :: stashcode_v_adv             = 257
Integer, Parameter     :: stashcode_w_adv             = 258
Integer, Parameter     :: stashcode_n_turb_mixlvs     = 259
Integer, Parameter     :: stashcode_lvl_bse_dp_sc     = 260
Integer, Parameter     :: stashcode_lvl_top_dp_sc     = 261
Integer, Parameter     :: stashcode_bl_conv_flag      = 262
Integer, Parameter     :: stashcode_turb_temp         = 263
Integer, Parameter     :: stashcode_turb_humid        = 264
Integer, Parameter     :: stashcode_area_cf           = 265
Integer, Parameter     :: stashcode_bulk_cf           = 266
Integer, Parameter     :: stashcode_liquid_cf         = 267
Integer, Parameter     :: stashcode_frozen_cf         = 268
Integer, Parameter     :: stashcode_sfc_zonal_cur     = 269
Integer, Parameter     :: stashcode_sfc_merid_cur     = 270
Integer, Parameter     :: stashcode_qcf2              = 271
Integer, Parameter     :: stashcode_qrain             = 272
Integer, Parameter     :: stashcode_qgraup            = 273
Integer, Parameter     :: stashcode_Ti_Mean           = 274
Integer, Parameter     :: stashcode_Ti_Sig            = 275
Integer, Parameter     :: stashcode_fexp              = 275
Integer, Parameter     :: stashcode_gamtot            = 277
Integer, Parameter     :: stashcode_zw                = 278
Integer, Parameter     :: stashcode_fsat              = 279
Integer, Parameter     :: stashcode_fwetl             = 280
Integer, Parameter     :: stashcode_sthzw             = 281
Integer, Parameter     :: stashcode_a_fsat            = 282
Integer, Parameter     :: stashcode_c_fsat            = 283
Integer, Parameter     :: stashcode_a_fwet            = 284
Integer, Parameter     :: stashcode_c_fwet            = 285
Integer, Parameter     :: stashcode_p                 = 407
Integer, Parameter     :: stashcode_pstar             = 409
Integer, Parameter     :: stashcode_ice_conc_cat      = 413
Integer, Parameter     :: stashcode_ice_thick_cat     = 414
Integer, Parameter     :: stashcode_ice_temp_cat      = 415
Integer, Parameter     :: stashcode_ice_snow_depth    = 416
Integer, Parameter     :: stashcode_dust1_mmr         = 431
Integer, Parameter     :: stashcode_dust2_mmr         = 432
Integer, Parameter     :: stashcode_dust3_mmr         = 433
Integer, Parameter     :: stashcode_dust4_mmr         = 434
Integer, Parameter     :: stashcode_dust5_mmr         = 435
Integer, Parameter     :: stashcode_dust6_mmr         = 436
Integer, Parameter     :: stashcode_soilcarb_dpm      = 466
Integer, Parameter     :: stashcode_soilcarb_rpm      = 467
Integer, Parameter     :: stashcode_soilcarb_bio      = 468
Integer, Parameter     :: stashcode_soilcarb_hum      = 469
Integer, Parameter     :: stashcode_ozone_tracer      = 480
Integer, Parameter     :: stashcode_o3_prod_loss      = 481
Integer, Parameter     :: stashcode_o3_p_l_vmr        = 482
Integer, Parameter     :: stashcode_o3_vmr            = 483
Integer, Parameter     :: stashcode_o3_p_l_temp       = 484
Integer, Parameter     :: stashcode_o3_temp           = 485
Integer, Parameter     :: stashcode_o3_p_l_colo3      = 486
Integer, Parameter     :: stashcode_o3_colo3          = 487
Integer, Parameter     :: stashcode_land_frac         = 505
Integer, Parameter     :: stashcode_tstar_land        = 506
Integer, Parameter     :: stashcode_tstar_sea         = 507
Integer, Parameter     :: stashcode_tstar_sice        = 508
Integer, Parameter     :: stashcode_albedo_sice       = 509
Integer, Parameter     :: stashcode_albedo_land       = 510

!----------------------------------------------------------
! Chemistry stashcodes - now section 34
!----------------------------------------------------------
Integer, Parameter     :: stashcode_NO2               = 34004
Integer, Parameter     :: stashcode_CH4               = 34009
Integer, Parameter     :: stashcode_CO                = 34010
Integer, Parameter     :: stashcode_HCHO              = 34011
Integer, Parameter     :: stashcode_O3                = 34001

END MODULE Rcf_Stashcodes_Mod

#endif
