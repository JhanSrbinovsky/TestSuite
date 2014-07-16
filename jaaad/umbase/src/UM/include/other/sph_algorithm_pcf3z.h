!     ------------------------------------------------------------------
!     Module to set the allowed spherical harmonic algorithms.
!
      INTEGER                                                           &
     &    IP_SPH_DIRECT                                                 &
!           Direct solution using spherical harmonics
     &  , IP_SPH_REDUCED_ITER
!           The spherical harmonic solution a reduced order of
!           truncation is used to define a source
!           term for integration along a line: this can be combined
!           with a higher order of truncation for the solar beam
!           to yield a solution almost identical to the TMS method
!           of Nakajima and Tanaka.
!
      PARAMETER(                                                        &
     &     IP_SPH_DIRECT=1                                              &
     &   , IP_SPH_REDUCED_ITER=2                                        &
     &   )
!
!     ------------------------------------------------------------------
