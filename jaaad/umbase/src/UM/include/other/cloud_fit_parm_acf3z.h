!     ------------------------------------------------------------------
!     Module to set parameters for cconjugate gradient cloud
!     parametrization.
!
      INTEGER                                                           &
     &    NP_MAX_ITERATION_CG                                           &
!           Maximum number of conjugate-gradient iterations
     &  , NP_MAX_LINE_SEARCH
!           Maximum number of line searches
      REAL  (real64) ::                                                 &
     &    S_MAX_RATIO
!           Maximum line search ratio
!
      PARAMETER(                                                        &
     &    NP_MAX_ITERATION_CG=1000                                      &
     &  , NP_MAX_LINE_SEARCH=100                                        &
     &  )
      PARAMETER(                                                        &
     &    S_MAX_RATIO=1.0E+03_real64                                    &
     &  )
!     ------------------------------------------------------------------
