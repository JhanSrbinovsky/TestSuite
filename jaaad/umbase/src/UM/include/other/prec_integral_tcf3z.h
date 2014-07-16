!     ------------------------------------------------------------------
!     Module to set tolerances for integration routines.
!
      REAL  (real64) ::                                                 &
     &    TOL_PANEL                                                     &
!           Threshold for neglect of panel
     &  , TOL_REFINEMENT
!           Tolerance on romberg step
!
      PARAMETER(                                                        &
     &    TOL_PANEL=1.0E-04_real64                                      &
     &  , TOL_REFINEMENT=1.0E-04_real64                                 &
     &  )
!
!     ------------------------------------------------------------------
