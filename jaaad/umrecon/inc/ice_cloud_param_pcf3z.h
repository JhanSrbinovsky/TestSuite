!     ------------------------------------------------------------------
!     Module to set numbers for ice cloud schemes.
!
      INTEGER, PARAMETER :: NPD_ICE_CLOUD_FIT           = 10
!           Number of cloud fitting schemes
      INTEGER, PARAMETER :: IP_SLINGO_SCHRECKER_ICE     = 1
!           Parametrization of Slingo and Schrecker.
      INTEGER, PARAMETER :: IP_ICE_UNPARAMETRIZED       = 3
!           Unparametrized ice crystal data
      INTEGER, PARAMETER :: IP_SUN_SHINE_VN2_VIS        = 4
!           Sun and Shine's parametrization in the visible (version 2)
      INTEGER, PARAMETER :: IP_SUN_SHINE_VN2_IR         = 5
!           Sun and Shine's parametrization in the IR (version 2)
      INTEGER, PARAMETER :: IP_ICE_ADT                  = 6
!           ADT-based scheme for ice crystals
      INTEGER, PARAMETER :: IP_ICE_ADT_10               = 7
!           ADT-based scheme for ice crystals using 10th order
!           polynomials
#if ! defined(UM)
      INTEGER, PARAMETER :: IP_ICE_PARAMETRIZATION_TEST = 8
!           Test parametrization for ice crystals
#endif
      INTEGER, PARAMETER :: IP_ICE_FU_SOLAR             = 9
!           Fu's parametrization in the solar region of the spectrum
      INTEGER, PARAMETER :: IP_ICE_FU_IR                = 10
!           Fu's parametrization in the infra-red region of the spectrum
      INTEGER, PARAMETER :: IP_SLINGO_SCHR_ICE_PHF      = 11
!           Parametrization of Slingo and Schrecker
!           (Moments of phase function).
      INTEGER, PARAMETER :: IP_ICE_AGG_DE               = 12
!           Provisional agregate parametrization.
      INTEGER, PARAMETER :: IP_ICE_FU_PHF               = 13
!           Parametrization like Fu
!           (Moments of phase function).
!     ------------------------------------------------------------------
