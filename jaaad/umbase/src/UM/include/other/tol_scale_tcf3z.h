!     ------------------------------------------------------------------
!     Module to set tolerances for finding scaling function.
!
      REAL  (real64) ::                                                 &
     &    TOL_TRANSMISSION                                              &
!           Tolerance for linearization
     &  , TOL_OBJECTIVE
!           Tolerance for objective function
!
      PARAMETER(                                                        &
     &    TOL_TRANSMISSION=1.00E-08_real64                              &
     &  , TOL_OBJECTIVE=1.00E-08_real64                                 &
     &  )
!
!     ------------------------------------------------------------------
