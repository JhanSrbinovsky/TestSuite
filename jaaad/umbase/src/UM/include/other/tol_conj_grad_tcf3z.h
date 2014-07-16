!     ------------------------------------------------------------------
!     Module to set tolerances for convergence of the
!     conjugate gradient algorithms
!
      REAL  (real64) ::                                                 &
     &    TOL_RESIDUAL                                                  &
!           Tolerance on residual
     &  , TOL_GRADIENT
!           Tolerance on gradient
!
      PARAMETER(                                                        &
     &    TOL_RESIDUAL=1.0E-10_real64                                   &
     &  , TOL_GRADIENT=1.0E-09_real64                                   &
     &  )
!
!     ------------------------------------------------------------------
