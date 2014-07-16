!     ------------------------------------------------------------------
!     Module to set parameters for scaling function and
!     associated conjugate gradient algorithm.
!
      INTEGER                                                           &
     &    NP_MAX_ITERATION_CG                                           &
!           Maximum number of conjugate gradient cycles
     &  , NP_MAX_LINE_SEARCH
!           Maximum number of line searches
      REAL  (real64) ::                                                 &
     &    S_MAX_RATIO
!           Maximum line search parameter
!
      PARAMETER(                                                        &
     &    NP_MAX_ITERATION_CG=100                                       &
     &  , NP_MAX_LINE_SEARCH=100                                        &
     &  )
      PARAMETER(                                                        &
     &    S_MAX_RATIO=2.0E+04_real64                                    &
     &  )
!
!     ------------------------------------------------------------------
