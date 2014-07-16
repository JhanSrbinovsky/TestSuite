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
!   This module declares the controlling structures for SW radiation
!   and performs suitable default initialization.
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
    SUBROUTINE sw_control_default(sw_control)
!
      use control_struc

      TYPE (control_option) :: sw_control
!       The block of controlling options for the code
!
!
      sw_control%spectral_file=''

!     Range of bands
      sw_control%first_band=1
      sw_control%last_band=1
!
!     Miscellaneous options
      sw_control%isolir=1
      sw_control%i_scatter_method=1
      sw_control%l_ir_source_quad=.false.
      sw_control%l_extra_top=.false.
      sw_control%l_rad_deg=.false.
      sw_control%l_subsample=.false.
!
!     Properties of clouds
      sw_control%i_cloud=6
      sw_control%i_cloud_representation=4
      sw_control%i_st_water=0
      sw_control%i_cnv_water=0
      sw_control%i_st_ice=0
      sw_control%i_cnv_ice=0
      sw_control%l_local_cnv_partition=.false.
      sw_control%l_global_cloud_top=.true.
!
!     Physical processes
      sw_control%l_microphysics=.true.
      sw_control%l_gas=.true.
      sw_control%l_rayleigh=.true.
      sw_control%l_continuum=.true.
      sw_control%l_cloud=.true.
      sw_control%l_drop=.true.
      sw_control%l_ice=.true.
      sw_control%l_aerosol=.true.
      sw_control%l_aerosol_ccn=.true.
!
!     Gaseous absorption
      sw_control%i_gas_overlap=5
      sw_control%l_o2=.false.
      sw_control%l_n2o=.false.
      sw_control%l_ch4=.false.
      sw_control%l_cfc11=.false.
      sw_control%l_cfc12=.false.
      sw_control%l_cfc113=.false.
      sw_control%l_hcfc22=.false.
      sw_control%l_hfc125=.false.
      sw_control%l_hfc134a=.false.
!
!     Angular integration (including algorithmic options)
      sw_control%i_angular_integration=1
      sw_control%i_2stream=16
      sw_control%i_solver=14
      sw_control%n_order_gauss=0
      sw_control%i_truncation=1
      sw_control%i_sph_algorithm=1
      sw_control%n_order_phase_solar=1
      sw_control%ls_global_trunc=9
      sw_control%ms_min=0
      sw_control%ms_max=0
      sw_control%ls_brdf_trunc=0
      sw_control%accuracy_adaptive=1.0e-04
      sw_control%l_rescale=.true.
      sw_control%l_henyey_greenstein_pf=.true.
      sw_control%i_sph_mode=1
      sw_control%l_euler_trnf=.false.
!
!     Satellite data
!
      sw_control%l_geostationary=.false.
      sw_control%sat_desc="                    "  // &
                          "                    "  // &
                          "                    "  // &
                          "                    "
      sw_control%sat_hgt=0.0
      sw_control%sat_lon=0.0
      sw_control%sat_lat=0.0
      sw_control%max_view_lon=0.0
      sw_control%min_view_lon=0.0
      sw_control%max_view_lat=0.0
      sw_control%min_view_lat=0.0
!
      sw_control%l_layer=.true.
      sw_control%L_cloud_layer=.true.
      sw_control%l_2_stream_correct=.false.
   END SUBROUTINE sw_control_default
#endif
