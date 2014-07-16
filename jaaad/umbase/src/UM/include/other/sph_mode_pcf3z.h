!     ------------------------------------------------------------------
!     Module to set modes in which the spherical harmonic algorithm
!     can be used.
!
      INTEGER, Parameter :: IP_SPH_MODE_RAD  = 1
!             Spherical harmonics are used to calculate radiances
      INTEGER, Parameter :: IP_SPH_MODE_FLUX = 2
!             SPherical harmonics are sude to calculate fluxes
      INTEGER, Parameter :: IP_SPH_MODE_J    = 3
!             Spherical harmonics are used to calculate mean
!             radiances (actinic flux/4 pi)
!
!     ------------------------------------------------------------------
