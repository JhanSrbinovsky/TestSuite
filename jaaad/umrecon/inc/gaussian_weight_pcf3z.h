!     ------------------------------------------------------------------
!     Module defining arrays of gaussian points and weights.
!
      INTEGER                                                           &
     &    NPD_GAUSS_ORD
!           Maximum order of gaussian quadrature
!
      PARAMETER(                                                        &
     &    NPD_GAUSS_ORD=10                                              &
     &  )
      REAL  (real64) ::                                                 &
     &    GAUSS_POINT(NPD_GAUSS_ORD, NPD_GAUSS_ORD)                     &
!           Points of gaussian integration
     &  , GAUSS_WEIGHT(NPD_GAUSS_ORD, NPD_GAUSS_ORD)
!           Weights of gaussian integration
!
!     ------------------------------------------------------------------
