!     ------------------------------------------------------------------
!     Module to set tolerances for preprocessing.
!
      REAL  (real64) ::                                                 &
     &    TOL_PRESS_PP                                                  &
!           Tolerance for equality of p
     &  , TOL_TEMP_PP                                                   &
!           Tolerance for equality of t
     &  , TOL_AMOUNT_PP                                                 &
!           Tolerance for equality of amount
     &  , TOL_CUTOFF_PP
!           Tolerance on line cutoff
!
      PARAMETER(                                                        &
     &    TOL_PRESS_PP=1.0E-05_real64                                   &
     &  , TOL_TEMP_PP=1.0E-05_real64                                    &
     &  , TOL_AMOUNT_PP=1.0E-05_real64                                  &
     &  , TOL_CUTOFF_PP=1.0E-05_real64                                  &
     &  )
!
!     ------------------------------------------------------------------
