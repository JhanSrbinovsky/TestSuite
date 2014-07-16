!     ------------------------------------------------------------------
!     MODULE SETTING MAXIMUM DIMENSIONS OF ARRAYS IN THE RADIATION CODE.
!   4.5   Aug 1998     Increment by 2 the no. of aerosol species
!                      affecting the radiation.    Luke Robinson
!
!
!   5.2   Dec 2000     Increment by 2 the no. of aerosol species
!                      affecting the radiation.    Andy Jones
!   5.3   Oct 2001     Allow six more types of aerosol in preparation
!                      for runs with dust.      (J. M. Edwards)
!   5.5   Feb 2003     Increase aerosol species to include dust
!                      and biomass    (S. Woodward)
!   5.5   Feb 2003     Increase the number of types of ice crystal.
!                                                 John Edwards
!   6.2   Nov 2005     Increase number of types of data.
!                      Add number of wavelengths for aerosol optical
!                      depth.
!                                                 N. Bellouin
!   6.2   Jan 2006     Increase values of the cloud parameters and
!                      of the number of terms in the phase function
!                      under versions 3C and 3Z of the radiation
!                      code.                         (J.-C. Thelen)
!
      INTEGER,PARAMETER :: NPD_TYPE      = 16 ! number of types of data
      INTEGER,PARAMETER :: NPD_BAND      = 20 ! number of spectral bands
      INTEGER,PARAMETER :: NPD_EXCLUDE   = 2  ! nubmer of excluded bands
      INTEGER,PARAMETER :: NPD_SPECIES   = 11 ! number of gaseous species
      INTEGER,PARAMETER :: NPD_ESFT_TERM = 16 ! number of esft terms
      INTEGER,PARAMETER :: NPD_SCALE_FNC = 4  ! number of scaling funcs

      ! number of scaling variables
      INTEGER,PARAMETER :: NPD_SCALE_VARIABLE=4

      INTEGER,PARAMETER :: NPD_SURFACE     = 1 ! no of surface types
      INTEGER,PARAMETER :: NPD_ALBEDO_PARM = 4 ! no of albedo parameters
      INTEGER,PARAMETER :: NPD_CONTINUUM   = 2 ! no of continua
      INTEGER,PARAMETER :: NPD_DROP_TYPE   = 6 ! no of drop types
      INTEGER,PARAMETER :: NPD_ICE_TYPE    = 8 ! no of ice crystal types


      ! max no of cloud parameters
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
      INTEGER,PARAMETER :: NPD_CLOUD_PARAMETER=504
#else
      INTEGER,PARAMETER :: NPD_CLOUD_PARAMETER=30
#endif

      INTEGER, PARAMETER :: NPD_AEROSOL_SPECIES=24

      ! maximum no of humidities
      INTEGER,PARAMETER :: NPD_HUMIDITIES=21

      ! no of thermal coefficients
      INTEGER,PARAMETER :: NPD_THERMAL_COEFF=9

      ! number of wavelengths for aerosol optical depth
      INTEGER,PARAMETER :: NPD_AOD_WAVEL=6
#if defined(A01_3C) || defined(A02_3C) \
 || defined(A01_3Z) || defined(A02_3Z)
!       Size allocated for terms in phase functions
      INTEGER,PARAMETER :: NPD_PHASE_TERM=101
#endif
! MXSIZE3A end
