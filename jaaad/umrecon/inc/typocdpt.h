! TYPOCDPT Pointers for dynamic allocation in ocean.
!  Pointer jocp_XXXX points to start of array XXX in super-array
!
! 4.5  14/08/97  Removed the pointers, jocp_tbound_n,... to the
!                boundary arrays (TBOUND_N,...) C.G. Jones
! 5.2  31/07/00  Included jocp_fkmz. M J Bell
! 5.4  31/08/02  Included jocp_cislbdy, jocp_iislby, jocp_jislbdy,
!                jocp_pe_islbdy. D.Storkey
! 5.4  29/08/02  Included jocp_t_strait, jocp_flux_strait,
!                jocp_t_strait_clm, jocp_fkmqx. D. Storkey
! 5.5  28/02/03  Included jocp_ainmin and jocp_hinmax.  A.McLaren
!  6.2  Jun 2004  Include JP*Z pointers for IFS. R. Hill
!  6.2  Jun 2004  Include JP_MIK*, MAK* and MYSL to hold F.
!                 filtering master/slave indexing. R. Hill

      ! Pointers for TYPOCFLW
      INTEGER :: jocp_kar
      INTEGER :: jocp_isz
      INTEGER :: jocp_iez
      INTEGER :: jocp_ise
      INTEGER :: jocp_iee
      INTEGER :: jocp_isu
      INTEGER :: jocp_ieu
      INTEGER :: jocp_lse
      INTEGER :: jocp_lsu
      INTEGER :: jocp_iseg
      INTEGER :: jocp_istf
      INTEGER :: jocp_ietf
      INTEGER :: jocp_isuf
      INTEGER :: jocp_ieuf
      INTEGER :: jocp_iszf
      INTEGER :: jocp_iezf
      INTEGER :: jocp_spsin
      INTEGER :: jocp_spcos
      INTEGER :: jocp_isis
      INTEGER :: jocp_ieis
      INTEGER :: jocp_jsis
      INTEGER :: jocp_jeis
      INTEGER :: jocp_cislbdy
      INTEGER :: jocp_iislbdy
      INTEGER :: jocp_jislbdy
      INTEGER :: jocp_pe_islbdy
      INTEGER :: jocp_t_strait
      INTEGER :: jocp_flux_strait
      INTEGER :: jocp_t_strait_clm
      ! Pointers for TYPOCONE
      INTEGER :: jocp_dxt
      INTEGER :: jocp_dxtr
      INTEGER :: jocp_dxt2r
      INTEGER :: jocp_dxu
      INTEGER :: jocp_dxur
      INTEGER :: jocp_dxu2r
      INTEGER :: jocp_dxu4r
      INTEGER :: jocp_dxt4r
      INTEGER :: jocp_dyt
      INTEGER :: jocp_dytr
      INTEGER :: jocp_dyt2r
      INTEGER :: jocp_dyu
      INTEGER :: jocp_dyur
      INTEGER :: jocp_dyu2r
      INTEGER :: jocp_dyu2rj
      INTEGER :: jocp_dyu4r
      INTEGER :: jocp_dyt4r
      INTEGER :: jocp_cs
      INTEGER :: jocp_csr
      INTEGER :: jocp_csrj
      INTEGER :: jocp_cst
      INTEGER :: jocp_cstr
      INTEGER :: jocp_phi
      INTEGER :: jocp_phit
      INTEGER :: jocp_sine
      INTEGER :: jocp_tng
      INTEGER :: jocp_c2dz
      INTEGER :: jocp_dz
      INTEGER :: jocp_dz2r
      INTEGER :: jocp_eeh
      INTEGER :: jocp_eem
      INTEGER :: jocp_ffh
      INTEGER :: jocp_ffm
      INTEGER :: jocp_zdz
      INTEGER :: jocp_dzz
      INTEGER :: jocp_dzz2r
      INTEGER :: jocp_zdzz
      INTEGER :: jocp_sol_pen
      INTEGER :: jocp_dttsa
      INTEGER :: jocp_rz
      INTEGER :: jocp_rzz
      INTEGER :: jocp_rzz2r
      INTEGER :: jocp_delpsl
      INTEGER :: jocp_decay
      INTEGER :: jocp_ahi
      INTEGER :: jocp_amt
      INTEGER :: jocp_amu
      INTEGER :: jocp_kappabsi
      INTEGER :: jocp_cosine
      INTEGER :: jocp_rlambda
      INTEGER :: jocp_eddydiff
      INTEGER :: jocp_amx
      INTEGER :: jocp_ainmin
      INTEGER :: jocp_hinmax
      INTEGER :: jocp_bbu
      INTEGER :: jocp_ccu
      INTEGER :: jocp_ddu
      INTEGER :: jocp_ggu
      INTEGER :: jocp_hhu
      INTEGER :: jocp_athkdf
      INTEGER :: jocp_kri
      INTEGER :: jocp_csrjp
      INTEGER :: jocp_dyu2rjp
      INTEGER :: jocp_cstjp
      INTEGER :: jocp_dytrjp
      INTEGER :: jocp_csjm
      INTEGER :: jocp_dyurjm
      INTEGER :: jocp_dyujm
      INTEGER :: jocp_dyt2rjm
      INTEGER :: jocp_cstrjp
      INTEGER :: jocp_cstrjm
      INTEGER :: jocp_cstjm
      INTEGER :: JOCP_MAXLARGELEVELS
      INTEGER :: JOCP_NOLEVSINLAYER
      ! Pointers for TYPOCFLD
      INTEGER :: jocp_hr
      INTEGER :: jocp_hrj
      INTEGER :: jocp_fkmp
      INTEGER :: jocp_fkmq_global
      INTEGER :: jocp_fkmq
      INTEGER :: jocp_fkmqx
      INTEGER :: jocp_fkmz
      INTEGER :: jocp_coriolis
      INTEGER :: jocp_em
      INTEGER :: jocp_hrjp
      INTEGER :: jocp_pjp
      INTEGER :: jocp_pbjp
      INTEGER :: jocp_fkmqjp
      INTEGER :: JOCP_EMU
      INTEGER :: JOCP_UBTBBCJP
      INTEGER :: JOCP_VBTBBCJP
      INTEGER :: JOCP_UBTJP
      INTEGER :: JOCP_VBTJP
      ! Pointers for TYPOCFIL
      INTEGER :: jocp_icbase
      INTEGER :: jocp_idbase
      INTEGER :: jocp_ind
      INTEGER :: jocp_cossav
      INTEGER :: jocp_denmsv
      INTEGER :: jocp_cosnpi
      INTEGER :: JP_MCU
      INTEGER :: JP_MCT
      INTEGER :: JP_MCF
      INTEGER :: JP_MPU
      INTEGER :: JP_MPT
      INTEGER :: JP_MPF
      INTEGER :: JP_MKU
      INTEGER :: JP_MKT
      INTEGER :: JP_MSU
      INTEGER :: JP_MST
      INTEGER :: JP_MSF
      INTEGER :: JP_MRF
      INTEGER :: JP_SCU
      INTEGER :: JP_SCT
      INTEGER :: JP_SCF
      INTEGER :: JP_MIKU,JP_MIKT,JP_MAKU,JP_MAKT,JP_MYSL
      INTEGER :: JP_MPZ
      INTEGER :: JP_MJZ
      INTEGER :: JP_MLZ
      INTEGER :: JP_SPZ, JP_SJZ, JP_SLZ
      ! Pointers for TYPOCMEA
      INTEGER :: jocp_isht
      INTEGER :: jocp_ieht
      ! Pointers for TYPOCAC
      INTEGER :: jocp_o_obs
      INTEGER :: jocp_aice_obs
      INTEGER :: jocp_ice_obs_dys
      INTEGER :: jocp_ice_obs_scs
      INTEGER :: jocp_o_lon_m
      INTEGER :: jocp_o_lat_m
      INTEGER :: jocp_o_dep_levs_m
      ! Pointers for TYPOCBIO
      INTEGER :: jocp_latbgc
      INTEGER :: jocp_dlco
      ! Pointers for TYPOASZ
      ! blank at moment
#include "comocdpt.h"
! TYPOCDPT end
