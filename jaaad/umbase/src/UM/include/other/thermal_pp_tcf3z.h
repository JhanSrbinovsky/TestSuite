!     ------------------------------------------------------------------
!     Module to set precision in calculating the thermal source.
!
      REAL  (real64) ::                                                 &
     &    P_MIN_THERMAL_RANGE                                           &
!           Minium range of temperature
     &  , TOL_THERMAL
!           Tolerance on integration
!
      PARAMETER(                                                        &
     &    P_MIN_THERMAL_RANGE=1.0E+00_real64                            &
     &  , TOL_THERMAL=1.0E-06_real64                                    &
     &  )
!
!     ------------------------------------------------------------------
