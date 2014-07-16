!=========================== COMDECK COMOCDPT =========================
!LL  4.5  14/08/97  Removed the pointers, jocp_tbound_n,... to the
!LL                 boundary arrays (TBOUND_N,...) C.G. Jones
!LL  4.5  3/11/98  added pointers jocp_csrjp ... for remote arrays
!LL                from this PE
!LL  5.2  31/07/00 Included jocp_fkmz. M J Bell
!LL  5.4  31/08/02 Included jocp_cislbdy, jocp_iislbdy, jocp_jislbdy,
!LL                jocp_pe_islbdy. D.Storkey
!LL  5.4  29/08/02 Included jocp_t_strait, jocp_flux_strait,
!                  jocp_t_strait_clm, jocp_fkmqx. D. Storkey
!    5.5  28/02/03 Included jocp_ainmin and jocp_hinmax. A.McLaren
!  6.2  Jun 2004  Include JP*Z pointers for IFS. R. Hill
!  6.2  Jun 2004  Include JP_MIK*, MAK* and MYSL to hold F.
!                 filtering master/slave indexing. R. Hill
!    6.2  11/08/05 Fix continuations. P.Selwood

      ! Pointers for file TYPOCFLW
      COMMON /COMOCDPT/                                                 &
     &  jocp_kar,jocp_isz,jocp_iez,jocp_ise,jocp_iee,jocp_isu, jocp_ieu,&
     &  jocp_lse,jocp_lsu,jocp_iseg,jocp_istf,jocp_ietf,jocp_isuf,      &
     &  jocp_ieuf,jocp_iszf,jocp_iezf,jocp_spsin,jocp_spcos,jocp_isis,  &
     &  jocp_ieis,jocp_jsis,jocp_jeis,jocp_cislbdy,jocp_iislbdy,        &
     &  jocp_jislbdy,jocp_pe_islbdy,                                    &
     &  jocp_t_strait,jocp_flux_strait,jocp_t_strait_clm

      ! Pointers for file TYPOCONE


      COMMON /COMOCDPT/                                                 &
     &  jocp_dxt,jocp_dxtr,jocp_dxt2r,jocp_dxu,jocp_dxur,               &
     &  jocp_dxu2r,jocp_dxu4r,jocp_dxt4r,jocp_dyt,jocp_dytr,            &
     &  jocp_dyt2r,jocp_dyu,jocp_dyur,jocp_dyu2r,jocp_dyu2rj,           &
     &  jocp_dyu4r,jocp_dyt4r,jocp_cs,jocp_csr,jocp_csrj,               &
     &  jocp_cst,jocp_cstr,jocp_phi,jocp_phit,jocp_sine,                &
     &  jocp_tng,jocp_c2dz,jocp_dz,jocp_dz2r,jocp_eeh,                  &
     &  jocp_eem,jocp_ffh,jocp_ffm,jocp_zdz,jocp_dzz,                   &
     &  jocp_dzz2r,jocp_zdzz,jocp_sol_pen,jocp_dttsa,jocp_rz,           &
     &  jocp_rzz,jocp_rzz2r,jocp_delpsl,jocp_decay,                     &
     &  jocp_ahi,jocp_amt,jocp_amu,jocp_kappabsi,jocp_cosine,           &
     &  jocp_rlambda,jocp_eddydiff,jocp_amx,                            &
     &  jocp_ainmin,jocp_hinmax,jocp_bbu,jocp_ccu,                      &
     &  jocp_ddu,jocp_ggu,jocp_hhu,jocp_athkdf,jocp_kri,                &
     &  jocp_csrjp,jocp_dyu2rjp,jocp_cstjp,jocp_dytrjp,jocp_csjm,       &
     &  jocp_dyurjm,                                                    &
     &  jocp_dyujm,                                                     &
     &  jocp_dyt2rjm,                                                   &
     &  jocp_cstrjp,                                                    &
     &  jocp_cstrjm,                                                    &
     &  jocp_cstjm,                                                     &
     &  JOCP_MAXLARGELEVELS,JOCP_NOLEVSINLAYER

      ! Pointers for file TYPOCFLD

      COMMON /COMOCDPT/                                                 &
     &  jocp_hr,jocp_hrj,jocp_fkmp,jocp_fkmq_global,jocp_fkmq,          &
     &  jocp_fkmqx,jocp_fkmz,jocp_coriolis,jocp_em,jocp_hrjp,jocp_pjp,  &
     &  jocp_pbjp,jocp_fkmqjp,JOCP_EMU,                                 &
     &  JOCP_UBTBBCJP,JOCP_VBTBBCJP,JOCP_UBTJP,JOCP_VBTJP
      !  Pointers for file TYPOCFIL

      COMMON /COMOCDPT/                                                 &
     &  jocp_icbase,jocp_idbase,jocp_ind,jocp_cossav,jocp_denmsv,       &
     &  jocp_cosnpi,                                                    &
     &  JP_MCU,JP_MCT,JP_MCF,JP_MPU,JP_MPT,JP_MPF,JP_MKU,JP_MKT,        &
     &  JP_MSU,JP_MST,JP_MSF,JP_MRF,JP_SCU,JP_SCT,JP_SCF,               &
     & JP_MIKU,JP_MIKT,JP_MAKU,JP_MAKT,JP_MYSL,                         &
     &  JP_MPZ,JP_MJZ,JP_MLZ,JP_SPZ,JP_SJZ,JP_SLZ

      ! Pointers for file TYPOCMEA

      COMMON /COMOCDPT/jocp_isht,jocp_ieht

      ! Pointers for file TYPOCAC

      COMMON /COMOCDPT/                                                 &
     &  jocp_o_obs,jocp_aice_obs,jocp_ice_obs_dys,jocp_ice_obs_scs,     &
     &  jocp_o_lon_m,jocp_o_lat_m,jocp_o_dep_levs_m

      ! Pointers for file TYPOCBIO

      COMMON /COMOCDPT/jocp_latbgc,jocp_dlco

! COMOCDPT end
