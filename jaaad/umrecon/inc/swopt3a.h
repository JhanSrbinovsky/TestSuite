! SWOPT3A defines SW options for the radiation code.
!
!     NOTE: SWOPT3A and SWCAVR3A must be consistent.
!
!     Current Owner of Code: J. M. Edwards
!
!     History:
!     Version  Date      Comment.
!     5.2      05/12/00  Logical for extra top level added.
!                                J. M. Edwards
! 6.2  24/10/05 Functionality for radiative forcing, timestepping
!               and radiances under versions 3C and 3Z of radiation
!               code added                 (J.-C. Thelen)
!
!
#if defined(A01_3C) || defined(A01_3Z)
      CHARACTER  (LEN=132) :: spectral_file_sw
!       Spectral file
      INTEGER :: first_band_sw
!       First band in sw calculations
      INTEGER :: last_band_sw
!       Last band in sw calculations
!
#endif
      INTEGER :: I_2STREAM_SW     ! two_stream scheme selected
      INTEGER :: I_GAS_OVERLAP_SW ! treatment of gaseous overlaps
      INTEGER :: I_CLOUD_SW       ! treatment of cloud overlaps
      INTEGER :: I_CLOUD_REPRESENTATION_SW  ! representation of clouds
      INTEGER :: I_SOLVER_SW                ! solver selected

      LOGICAL :: L_O2_SW ! flag for oxygen
#if defined(A01_3C) || defined(A01_3Z)
      INTEGER :: i_angular_integration_sw
!       Method of angular integration
      INTEGER :: i_truncation_sw
!       Spherical truncation selected
      INTEGER :: ls_global_trunc_sw
!       Order of truncation selected
      INTEGER :: ms_min_sw
!       Minimum azimuthal order
      INTEGER :: ms_max_sw
!       Maximum azimuthal order
      INTEGER :: ls_brdf_trunc_sw
!       Order of truncation applied to brdfs
      INTEGER :: i_sph_algorithm_sw
!       Algorithm applied in the spherical harmonic solution
      INTEGER :: n_order_phase_solar_sw
!       Number of terms retained in the solar phase function
      REAL  :: accuracy_adaptive_sw
!       Accuracy of the adaptive truncation in the sw
#endif

#if defined(A01_3C) || defined(A01_3Z)
      LOGICAL :: l_henyey_greenstein_pf_sw
!       Flag for henyey_greenstein phase functions
      INTEGER :: i_sph_mode_sw
!       Mode of application of spherical harmonic code
      LOGICAL :: l_euler_trnf_sw
!       Flag to apply euler's transformation to sums of
!       amplitudes of harmonics
!
#endif
      ! types of droplets or ice crystals used for parametrizations
      INTEGER :: I_ST_WATER_SW  ! type for stratiform water
      INTEGER :: I_CNV_WATER_SW ! type for convective water
      INTEGER :: I_ST_ICE_SW    ! type for stratiform ice
      INTEGER :: I_CNV_ICE_SW   ! type for convective ice

      ! flag to partition convective clouds using the local temperature
      LOGICAL :: L_LOCAL_CNV_PARTITION_SW


      ! Flag to include an extra layer at the top of the atmosphere in
      ! radiative calculations
      LOGICAL :: L_EXTRA_TOP_SW

! SWOPT3A end
