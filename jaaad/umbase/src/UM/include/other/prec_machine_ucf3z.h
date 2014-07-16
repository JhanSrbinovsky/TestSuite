!     ------------------------------------------------------------------
!     Module for setting machine precision.
!
      REAL  (real64) ::                                                 &
     &    TOL_MACHINE                                                   &
!           Machine tolerance
     &  , SQRT_TOL_MACHINE
!           Square root of machine tolerance
!
!
!     The precision should be about 2/2^(size of significand)
!
!     The IEEE-format uses 53 bits for the significand
!     in double precision
!
!     The CRAY format uses 47 bits in single precision.
!
      PARAMETER(                                                        &
     &    TOL_MACHINE=2.22E-16_real64                                   &
     &  , SQRT_TOL_MACHINE=1.5E-8                                       &
     &  )
!
!     ------------------------------------------------------------------
