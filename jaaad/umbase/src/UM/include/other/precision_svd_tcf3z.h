!     ------------------------------------------------------------------
!     Module to set tolerances for SVD algorithm.
!
      REAL  (real64) ::                                                 &
     &    TOL_SVD
!           Tolerance for SVD decomposition
!
      PARAMETER(                                                        &
     &    TOL_SVD=1.0E+03_real64*TOL_MACHINE                            &
     &  )
!
!     ------------------------------------------------------------------
