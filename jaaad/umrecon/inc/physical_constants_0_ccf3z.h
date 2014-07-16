!     ------------------------------------------------------------------
!     Module setting physical constants.
!
      REAL  (real64) ::                                                 &
     &    MOL_WEIGHT_AIR                                                &
!           Molar weight of dry air
#if defined(STANDARD)
     &  , SECONDS_PER_DAY                                               &
!           Number of seconds in a day
#endif
     &  , N2_MASS_FRAC
!           Mass fraction of nitrogen
!
      PARAMETER(                                                        &
     &    MOL_WEIGHT_AIR=28.966E-03_real64                              &
#if defined(STANDARD)
     &  , SECONDS_PER_DAY=8.6400E+04_real64                             &
#endif
     &  , N2_MASS_FRAC=0.781E+00_real64                                 &
     &  )
!
!     ------------------------------------------------------------------
