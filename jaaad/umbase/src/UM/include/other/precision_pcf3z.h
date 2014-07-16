!     ------------------------------------------------------------------
!     Module for setting algorithmic tolerances.
#if defined(UM)
!     (The comdeck prmch3a must always be included before this comdeck.)
#endif
#if defined(STANDARD)
!     ('prec_machine.ucf' must always be included with this file.)
#endif
!
      REAL  (real64) ::                                                 &
     &    TOL_DIV                                                       &
!           Tolerance for division
     &  , TOL_TEST
!           Tolerance for testing equality
!
      PARAMETER(                                                        &
     &    TOL_DIV=3.2E+01_real64*TOL_MACHINE                            &
     &  , TOL_TEST=1.6E+01_real64*TOL_MACHINE                           &
     &  )
!
!     ------------------------------------------------------------------
