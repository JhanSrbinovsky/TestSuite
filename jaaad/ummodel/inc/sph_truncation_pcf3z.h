!     ------------------------------------------------------------------
!     Module to set the types of spherical truncation.
!
      INTEGER                                                           &
     &     IP_TRUNC_TRIANGULAR                                          &
!             Trapezoidal truncation
     &   , IP_TRUNC_RHOMBOHEDRAL                                        &
!             Rhombohedral truncation
     &   , IP_TRUNC_AZIM_SYM                                            &
!             Truncation with azimuthal symmetry
     &   , IP_TRUNC_ADAPTIVE
!             Truncation set adaptively
!
      PARAMETER(                                                        &
     &     IP_TRUNC_TRIANGULAR=1                                        &
     &   , IP_TRUNC_RHOMBOHEDRAL=2                                      &
     &   , IP_TRUNC_AZIM_SYM=3                                          &
     &   , IP_TRUNC_ADAPTIVE=4                                          &
     &   )
!
!     ------------------------------------------------------------------
