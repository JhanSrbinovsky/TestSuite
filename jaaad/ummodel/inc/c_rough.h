!!----------------------------------------------------------------------
!!!-----------COMDECK C_ROUGH FOR SUBROUTINE SF_EXCH----------
! Sea ice parameters
! Z0FSEA = roughness length for free convective heat and moisture
!          transport over the sea (m).
!          DUMMY VARIABLE - Only used in 7A boundary layer scheme
! Z0HSEA = roughness length for heat and moisture transport
!          over the sea (m).
! Z0MIZ  = roughness length for heat, moisture and momentum over
!          the Marginal Ice Zone (m).
! Z0SICE = roughness length for heat, moisture and momentum over
!          sea-ice (m).
      REAL :: Z0HSEA
      REAL :: Z0MIZ
      REAL :: Z0SICE

      COMMON  /RUN_BLICE/Z0HSEA,Z0MIZ,Z0SICE

!!----------------------------------------------------------------------
