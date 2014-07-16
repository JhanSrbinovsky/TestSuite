#if defined(ATMOS)
#if defined(CONTROL) || defined(SCMA)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ Global data module for switches/options concerned with convection.

MODULE cv_run_mod

  ! Description:
  !   Module containing runtime logicals/options used by the convection code.
  !
  ! Method:
  !   All switches/options which are contained in the &Run_convection 
  !   sub-namelist in the CNTLATM control file are declared in this module.
  !   Default values have been declared where appropriate.
  !   
  !   Any routine wishing to use these options may do so with the 'Use' 
  !   statement.
  !
  ! Current Code Owner: Convection Scheme 
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 v7.4 programming standards.
  !
  ! Declarations:

  IMPLICIT NONE
  SAVE

  !===========================================================================
  ! Integer parameters to remove use of magic numbers
  !===========================================================================

  ! Basis for 2d convective cloud amount
  !   Total Condensed Water(TCW) original 4A
  !   LES mb/wsc scalings (Grant and Lock, 2004) (shallow)
  !   Surface rain from respective convective cloud type (Dp and Md-level)
  Integer, Parameter :: cca2d_total_condensed_water = 0
  Integer, Parameter :: cca2d_grant_lock            = 1 ! Shallow cnv (CCRad)
  Integer, Parameter :: cca2d_srf_precip_rate       = 1 ! mid/deep cnv (CCRad)


  ! Convective cloud decay
  Integer, Parameter :: rad_decay_off           = 0
  Integer, Parameter :: rad_decay_full_timestep = 1
  Integer, Parameter :: rad_decay_conv_substep  = 2 ! CCRad only


  ! Convective cloud decay timescale
  Integer, Parameter :: cld_life_constant = 0
  Integer, Parameter :: cld_life_func_hgt = 1 ! CCRad only


  ! Application of convective cloud anvils
  Integer, Parameter :: anv_pressure      = 0
  Integer, Parameter :: anv_height        = 1
  Integer, Parameter :: anv_model_levels  = 2


  !===========================================================================
  ! Logical switches set from UMUI
  !===========================================================================

  ! Comments for TRUE status

  Logical :: l_fix_udfactor    ! Fix application of UD_FACTOR (lock)

  Logical :: l_convcld_hadgem1 ! Use HadGEM1 cloud changes to deep convective
                               ! cclwp

  Logical :: l_cloud_deep      ! Use depth criterion for convective anvils

  Logical :: l_no_dcrit        ! Removes critical precip. depth criterion on
                               ! convective clouds (Ticket #996)

  Logical :: l_mom             ! Use Convective Momentum Transport

  Logical :: l_scv             ! Use STPH_SCV (Stochastic physics?)

  Logical :: l_rp              ! Use STPH_RP scheme

  Logical :: l_rp2             ! Use STPH_RP2 scheme

  Logical :: l_conv4a          ! Using 4a convection scheme

  Logical :: l4a_kterm         ! Use kterm instead for Deep CMT (4A)

  Logical :: l_eman_dd         ! Use Emanuel downdraught scheme (4A)

  Logical :: l_sdxs            ! Allow parcel excess to be set to s.d. of
                               ! turbulent fluctuations in lowest model layer

  Logical :: l_xscomp          ! Apply compensating cooling and drying of
                               ! environment in initiating model layer for
                               ! parcel excess

  Logical :: l_ccw             ! Do not includ rain in calculation of
                               ! convective cloud water path.

  Logical :: l_cape            ! Use Convective Available Potential Energy
                               ! (CAPE) closure


  !===========================================================================
  ! Integer options set from UMUI
  !===========================================================================

  ! Convection integer options set from UMUI Convection Scheme

  Integer :: n_conv_calls       ! Number of calls to convection
                                ! per physics timestep
  Integer :: a_convect_segments ! No of batches used in convection
  Integer :: a_convect_seg_size ! Size of convection segments. Can be 
                                ! specified as an alternative to  No. of
                                ! segments. 
                                ! -ive by default which means turned off.

  Integer :: convection_option  ! 0 Convection Scheme disabled
                                ! 1/2 Historical, no longer in use
                                ! 3 Convection Scheme enabled

  Integer :: cld_life_opt         ! Convective cloud decay time
  Integer :: rad_cloud_decay_opt  ! Convective cloud decay
  Integer :: anv_opt              ! Anvil cloud basis


  ! Convection Scheme Options (5A)
  ! NOTE: These options were valid at the time of writing. They are used in
  !       development code(5A) and so very likely to change.
  !       Users should consult the Convection Group for current available
  !       options.

  Integer :: iconv_shallow   ! Shallow
                             !   0: no scheme,
                             !   1: Gregory-Rowntree scheme
                             !   2: Turbulence scheme (non-precipitating)
                             !   3: Turbulence scheme (precipitating)

  Integer :: iconv_congestus ! Congestus
                             !   0: no scheme,
                             !   1: Gregory-Rowntree scheme
                             !   2: Future use
   
  Integer :: iconv_mid       ! Mid-level
                             !   0: no scheme,
                             !   1: Gregory-Rowntree scheme
                             !   2: Future use

  Integer :: iconv_deep      ! Deep
                             !   0: no scheme,
                             !   1: Gregory-Rowntree scheme
                             !   2: Future use

  Integer :: deep_cmt_opt    ! Deep CMT tunings
                             !   0: 4A Scheme (turbulence based) (Default)
                             !   1: Operational 70 level (modified 4A scheme)
                             !   2: Gregory-Kershaw scheme
                             !   3: New turbulence based scheme.
                             !   4: Future use

  Integer :: mid_cmt_opt     ! Mid CMT scheme to be used
                             !   0: Gregory-Kershaw scheme (Default)
                             !   1: Diffusive scheme
                             
  Integer :: icvdiag         ! Diagnosis calculation options
                             !   0: 4A Scheme (Default)
                             !   1: Improved 4A Scheme
                             !   2: Future/Experimental (dilute parcel)

  Integer :: adapt           ! Adaptive detrainment/entrainment options
                             !   0: 4A Scheme (Default)
                             !   1: Detrainment (mid and deep)
                             !   2: Future/Experimental (En/Detrainment)


  Integer :: cape_opt        ! CAPE closure options
                             !   0: RH based timescale
                             !   1: RH based timescale (timestep limited)
                             !   2: Fixed timescale
                             !   3: RH based timescale (timestep limited)
                             !      plus w test
                             !   4: Area scaled
                             !   5: Experimental

  Integer :: cape_bottom     ! Start level for w_max in column
  Integer :: cape_top        ! End   level for w_max in column

  Integer :: sh_pert_opt     ! Initial perturbation method for shallow
                             ! cumulus
                             ! 0 = Original code
                             ! 1 = Revised  code

  Integer :: dd_opt          ! Downdraught scheme options           
                             ! 0 = Original code
                             ! 1 = Revised  code

  Integer :: termconv        ! Not in umui, but available for future use
                             ! 0 for default
                             ! 1 for modified termination condition

  Integer :: bl_cnv_mix      ! Options for mixing convective increments in the BL
                             ! 0: original code (default)
                             ! 1: only mix the increments from the initial 
                             !    parcel perturbation

  Integer :: cca2d_sh_opt    ! Method to evaluate cca2d (Shallow)
  Integer :: cca2d_md_opt    ! Method to evaluate cca2d (Mid-level)
  Integer :: cca2d_dp_opt    ! Method to evaluate cca2d (Deep)


  Integer :: ccw_for_precip_opt   ! Option controlling critical cloud water 
                                  ! for the formation of precipitation
                                  ! 0 - original code
                                  ! 1 - no Dcrit option
                                  ! 2 - Manoj's first function 
                                  ! 3 - Manoj's congestus function


  !===========================================================================
  ! Real values set from UMUI
  !===========================================================================

  Real :: cca_sh_knob     ! Scales Shallow cloud fraction (CCRad)
  Real :: cca_md_knob     ! Scales Mid     cloud fraction (CCRad)
  Real :: cca_dp_knob     ! Scales Deep    cloud fraction (CCRad)

  Real :: ccw_sh_knob     ! Scales Shallow cloud water (CCRad)
  Real :: ccw_md_knob     ! Scales Mid     cloud water (CCRad)
  Real :: ccw_dp_knob     ! Scales Deep    cloud water (CCRad)


  Real :: fixed_cld_life  ! Fixed convective cloud lifetime decay value
                          ! (seconds)

  Real :: cca_min         ! Threshold value of convective cloud fraction
                          ! below which cca has neglible radiative impact and
                          ! is reset to zero

  Real :: r_det           ! Parameter controlling adaptive detrainment -
                          ! default = 0.75 (operational)
                          ! HadGEM1a recommended 0.5

  Real :: cape_ts_w       ! Cape timescale for scaling by w_max

  Real :: cape_min        ! Scale dependent min cape

  Real :: w_cape_limit    ! Test w for scale dependent cape timescale

  Real :: mparwtr         ! Reservoir of convective cloud water left in layer
                          ! whilst rest is rained out (kg/kg)

  Real :: mid_cnv_pmin    ! The minimum pressure (max height) at which mid 
                          ! level convection is allowed to initiate (Pa)

  ! Scales with cca2d to determine convective cloud amount with anvil
  Real :: anvil_factor    ! x cca2d = max cca, capped at 1.0
  Real :: tower_factor    ! x cca2d = min cca

  Real :: ud_factor       ! Updraught factor used in calculation of convective
                          ! water path

  Real :: tice            ! Phase change temperature in plume

  Real :: qstice          ! Qsat at phase change temperature
                          ! (freeze/melting temperature)

  Real :: entcoef_min         ! Minimum entrainment rate coefficient
  Real :: entcoef_max         ! Maximum entrainment rate coefficient
  Real :: cape_timescale      ! Timescale for CAPE closure.
  Real :: cape_timescale_min  ! Minimum CAPE closure timescale
  Real :: cape_timescale_max  ! Maximum CAPE closure timescale


END MODULE cv_run_mod
#endif
#endif
