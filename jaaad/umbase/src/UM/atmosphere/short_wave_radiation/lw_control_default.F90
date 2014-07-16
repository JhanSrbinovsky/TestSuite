#if defined(A70_1C) || defined(A70_1Z) \
 || defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   This module defines the controlling structure for LW calculations
!   and perfroems suitable default initializations.
!
! Current Code Owner: J.-C. Thelen
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.2   20/08/04   Original code.
!                        J.-C. Thelen
!
! Code Description:
!   Language: FORTRAN 90
!
!- End of header
!
!
!   Subroutine to set the default values of the control structure.
!
SUBROUTINE lw_control_default(lw_control)
!
      USE control_struc

      TYPE (control_option) :: lw_control
!       The block of controlling options for the code
!
!
      lw_control%spectral_file=''

!     Default setting of range of bands
      lw_control%first_band=1
      lw_control%last_band=1
!
!     Miscellaneous options
      lw_control%isolir=2
      lw_control%i_scatter_method=1
      lw_control%l_ir_source_quad=.true.
      lw_control%l_extra_top=.false.
      lw_control%l_rad_deg=.false.
      lw_control%l_subsample=.false.
!
!     Properties of clouds
      lw_control%i_cloud=6
      lw_control%i_cloud_representation=4
      lw_control%i_st_water=0
      lw_control%i_cnv_water=0
      lw_control%i_st_ice=0
      lw_control%i_cnv_ice=0
      lw_control%l_local_cnv_partition=.false.
      lw_control%l_global_cloud_top=.true.
!
!     Physical processes
      lw_control%l_microphysics=.true.
      lw_control%l_gas=.true.
      lw_control%l_rayleigh=.false.
      lw_control%l_continuum=.true.
      lw_control%l_cloud=.true.
      lw_control%l_drop=.true.
      lw_control%l_ice=.true.
      lw_control%l_aerosol=.true.
      lw_control%l_aerosol_ccn=.true.
!
!     Gaseous absorption
      lw_control%i_gas_overlap=5
      lw_control%l_o2=.false.
      lw_control%l_n2o=.false.
      lw_control%l_ch4=.false.
      lw_control%l_cfc11=.false.
      lw_control%l_cfc12=.false.
      lw_control%l_cfc113=.false.
      lw_control%l_hcfc22=.false.
      lw_control%l_hfc125=.false.
      lw_control%l_hfc134a=.false.
!
!     Angular integration (including algorithmic options)
      lw_control%i_angular_integration=1
      lw_control%i_2stream=12
      lw_control%i_solver=15
      lw_control%n_order_gauss=0
      lw_control%i_truncation=3
      lw_control%i_sph_algorithm=1
      lw_control%n_order_phase_solar=1
      lw_control%ls_global_trunc=9
      lw_control%ms_min=0
      lw_control%ms_max=0
      lw_control%ls_brdf_trunc=0
      lw_control%accuracy_adaptive=1.0e-04
      lw_control%l_rescale=.true.
      lw_control%l_henyey_greenstein_pf=.true.
      lw_control%i_sph_mode=1
      lw_control%l_euler_trnf=.false.
!
!     Satellite data
      lw_control%l_geostationary=.false.
      lw_control%sat_desc="                    "  // &
                          "                    "  // &
                          "                    "  // &
                          "                    "
      lw_control%sat_hgt=0.0
      lw_control%sat_lon=0.0
      lw_control%sat_lat=0.0
      lw_control%max_view_lon=0.0
      lw_control%min_view_lon=0.0
      lw_control%max_view_lat=0.0
      lw_control%min_view_lat=0.0
!
      lw_control%l_layer=.true.
      lw_control%l_cloud_layer=.true.
      lw_control%l_2_stream_correct=.false.
!
!
END SUBROUTINE lw_control_default
!
#endif
