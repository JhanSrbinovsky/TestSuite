#if defined(OCEAN)
! ARGOCTOP
! 4.5  14/08/97  Removed the start addresses to the boundary
!                 arrays (TBOUND_N,...) C.G. Jones
! 5.2  31/07/00  Included jocp_fkmz. M J Bell
! 5.4  31/08/02  Included jocp_cislbdy, jocp_iislby, jocp_jislbdy,
!                jocp_pe_islbdy. D.Storkey
! 5.4  29/08/02  Included jocp_t_strait, jocp_flux_strait,
!                jocp_t_strait_clm, jocp_fkmqx. D.Storkey
! 5.5  28/02/03  Included jocp_ainmin and jocp_hinmax. A.McLaren
!  6.2  Jun 2004  Include pointers to JP_MIK*, MAK* and MYSL for
!                 filtering master/slave indexing. R. Hill
!  6.1  Aug 2004   Reduce number of continuation lines to avoid 99
!                  line limit. R. Hill

      ! Start addresses for items in ARGOCFLW
     &  O_SPCON(jocp_kar), O_SPCON(jocp_isz), O_SPCON(jocp_iez),        &
     &  O_SPCON(jocp_ise), O_SPCON(jocp_iee), O_SPCON(jocp_isu),        &
     &  O_SPCON(jocp_ieu), O_SPCON(jocp_lse), O_SPCON(jocp_lsu),        &
     &  O_SPCON(jocp_iseg), O_SPCON(jocp_istf), O_SPCON(jocp_ietf),     &
     &  O_SPCON(jocp_isuf), O_SPCON(jocp_ieuf), O_SPCON(jocp_iszf),     &
     &  O_SPCON(jocp_iezf), O_SPCON(jocp_spsin), O_SPCON(jocp_spcos),   &
     &  O_SPCON(jocp_isis), O_SPCON(jocp_ieis), O_SPCON(jocp_jsis),     &
     & O_SPCON(jocp_jeis),O_SPCON(jocp_cislbdy),O_SPCON(jocp_iislbdy),  &
     & O_SPCON(jocp_jislbdy),O_SPCON(jocp_pe_islbdy),                   &
     &  O_SPCON(jocp_t_strait),O_SPCON(jocp_flux_strait),               &
     &  O_SPCON(jocp_t_strait_clm),                                     &
      ! Start addresses for items in ARGOCONE
     &  O_SPCON(jocp_dxt), O_SPCON(jocp_dxtr), O_SPCON(jocp_dxt2r),     &
     &  O_SPCON(jocp_dxu), O_SPCON(jocp_dxur), O_SPCON(jocp_dxu2r),     &
     &  O_SPCON(jocp_dxu4r), O_SPCON(jocp_dxt4r), O_SPCON(jocp_dyt),    &
     &  O_SPCON(jocp_dytr), O_SPCON(jocp_dyt2r), O_SPCON(jocp_dyu),     &
     &  O_SPCON(jocp_dyur), O_SPCON(jocp_dyu2r), O_SPCON(jocp_dyu2rj),  &
     &  O_SPCON(jocp_dyu4r), O_SPCON(jocp_dyt4r), O_SPCON(jocp_cs),     &
     &  O_SPCON(jocp_csr), O_SPCON(jocp_csrj), O_SPCON(jocp_cst),       &
     &  O_SPCON(jocp_cstr), O_SPCON(jocp_phi), O_SPCON(jocp_phit),      &
     &  O_SPCON(jocp_sine), O_SPCON(jocp_tng), O_SPCON(jocp_c2dz),      &
     &  O_SPCON(jocp_dz), O_SPCON(jocp_dz2r), O_SPCON(jocp_eeh),        &
     &  O_SPCON(jocp_eem), O_SPCON(jocp_ffh), O_SPCON(jocp_ffm),        &
     &  O_SPCON(jocp_zdz), O_SPCON(jocp_dzz), O_SPCON(jocp_dzz2r),      &
     &  O_SPCON(jocp_zdzz), O_SPCON(jocp_sol_pen), O_SPCON(jocp_dttsa), &
     &  O_SPCON(jocp_rz), O_SPCON(jocp_rzz),O_SPCON(jocp_rzz2r),        &
     &  O_SPCON(jocp_delpsl), O_SPCON(jocp_decay),O_SPCON(jocp_ahi),    &
     &  O_SPCON(jocp_amt), O_SPCON(jocp_amu),O_SPCON(jocp_kappabsi),    &
     &  O_SPCON(jocp_cosine),O_SPCON(jocp_rlambda),                     &
     &  O_SPCON(jocp_eddydiff),O_SPCON(jocp_amx),                       &
     &  O_SPCON(jocp_ainmin),O_SPCON(jocp_hinmax),O_SPCON(jocp_bbu),    &
     &  O_SPCON(jocp_ccu),O_SPCON(jocp_ddu),O_SPCON(jocp_ggu),          &
     &  O_SPCON(jocp_hhu),O_SPCON(jocp_athkdf),O_SPCON(jocp_kri),       &
     &  O_SPCON(JOCP_MAXLARGELEVELS),O_SPCON(JOCP_NOLEVSINLAYER),       &
     &  O_SPCON(jocp_csrjp),O_SPCON(jocp_dyu2rjp),O_SPCON(jocp_cstjp),  &
     &  O_SPCON(jocp_dytrjp),O_SPCON(jocp_csjm),O_SPCON(jocp_dyurjm),   &
     &  O_SPCON(jocp_dyujm),O_SPCON(jocp_dyt2rjm),                      &
     &  O_SPCON(jocp_cstrjp),O_SPCON(jocp_cstrjm),O_SPCON(jocp_cstjm),  &
      ! Start addresses for items in ARGOCFLD
     &  O_SPCON(jocp_hr), O_SPCON(jocp_hrj), O_SPCON(jocp_fkmp),        &
     & O_SPCON(jocp_fkmq_global),O_SPCON(jocp_fkmq),O_SPCON(jocp_fkmqx),&
     &  O_SPCON(jocp_fkmz),O_SPCON(jocp_coriolis),O_SPCON(jocp_em),     &
     &  O_SPCON(jocp_hrjp),O_SPCON(jocp_pjp),O_SPCON(jocp_pbjp),        &
     &  O_SPCON(jocp_fkmqjp),O_SPCON(JOCP_EMU),O_SPCON(JOCP_UBTBBCJP),  &
     &  O_SPCON(JOCP_VBTBBCJP),O_SPCON(JOCP_UBTJP),O_SPCON(JOCP_VBTJP), &
      ! Start addresses for items in ARGOCFIL
     &  O_SPCON(jocp_icbase),O_SPCON(jocp_idbase),O_SPCON(jocp_ind),    &
     &  O_SPCON(jocp_cossav),O_SPCON(jocp_denmsv),O_SPCON(jocp_cosnpi), &
     &  O_SPCON(JP_MCU),O_SPCON(JP_MCT),O_SPCON(JP_MCF),O_SPCON(JP_MPU),&
     &  O_SPCON(JP_MPT),O_SPCON(JP_MPF),O_SPCON(JP_MKU),O_SPCON(JP_MKT),&
     &  O_SPCON(JP_MSU),O_SPCON(JP_MST),O_SPCON(JP_MSF),O_SPCON(JP_MRF),&
     &O_SPCON(JP_SCU),O_SPCON(JP_SCT),O_SPCON(JP_SCF),O_SPCON(JP_MIKU), &
     &O_SPCON(JP_MIKT),O_SPCON(JP_MAKU), O_SPCON(JP_MAKT),              &
     &O_SPCON(JP_MYSL),O_SPCON(JP_MPZ),O_SPCON(JP_MJZ),O_SPCON(JP_MLZ), &
     &O_SPCON(JP_SPZ),O_SPCON(JP_SJZ),O_SPCON(JP_SLZ),                  &

      ! Start addresses for items in ARGOCMEA, ARGOCAC, ARGOCBIO
     &O_SPCON(jocp_isht),O_SPCON(jocp_ieht),O_SPCON(jocp_o_obs),O_SPCON(&
     &jocp_aice_obs),O_SPCON(jocp_ice_obs_dys),O_SPCON(jocp_ice_obs_scs)&
     &,O_SPCON(jocp_o_lon_m),O_SPCON(jocp_o_lat_m),O_SPCON(             &
     &jocp_o_dep_levs_m),O_SPCON(jocp_latbgc), O_SPCON(jocp_dlco),      &
! ARGOCTOP end
#endif
