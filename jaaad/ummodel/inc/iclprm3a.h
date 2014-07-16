! ICLPRM3A defines numbers for ice cloud schemes in two-stream radiation
! code.
!
!   5.5   Feb 2003     Addition of the aggregate parametrization.
!                                                 John Edwards
!   6.2   Jan 2006     Various options for radiation code
!                      3Z added.   (J.-C. Thelen)
!

      ! number of cloud fitting schemes
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      INTEGER,PARAMETER:: NPD_ICE_CLOUD_FIT=10
#else
      INTEGER,PARAMETER:: NPD_ICE_CLOUD_FIT=6
#endif

      ! parametrization of slingo and schrecker.
      INTEGER,PARAMETER:: IP_SLINGO_SCHRECKER_ICE=1

      ! unparametrized ice crystal data
       INTEGER,PARAMETER:: IP_ICE_UNPARAMETRIZED=3

      ! sun and shine's parametrization in the visible (version 2)
      INTEGER,PARAMETER:: IP_SUN_SHINE_VN2_VIS=4

      ! sun and shine's parametrization in the ir (version 2)
      INTEGER,PARAMETER:: IP_SUN_SHINE_VN2_IR=5

      ! scheme based on anomalous diffraction theory for ice crystals
      INTEGER,PARAMETER:: IP_ICE_ADT=6

#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
!           ADT-based scheme for ice crystals using 10th order
!           polynomials
      INTEGER,PARAMETER:: IP_ICE_ADT_10=7
#endif
      ! Provisional agregate parametrization.
      INTEGER,PARAMETER:: IP_ICE_AGG_DE=12
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
!           Fu's parametrization in the solar region of the spectrum
      INTEGER,PARAMETER:: IP_ICE_FU_SOLAR=9
!           Fu's parametrization in the infra-red region of the spectrum
      INTEGER,PARAMETER:: IP_ICE_FU_IR=10
!           Parametrization of Slingo and Schrecker
!           (Moments of phase function).
      INTEGER,PARAMETER:: IP_SLINGO_SCHR_ICE_PHF=11
!           Parametrization like Fu
!           (Moments of phase function).
      INTEGER,PARAMETER:: IP_ICE_FU_PHF=13
#endif

! ICLPRM3A end
