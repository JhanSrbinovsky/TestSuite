! LWFIX3A defining options to the edwards-slingo radiation code
! fixed in the unified model. options for longwave calculations.
!     Current Owner of Code: J. M. Edwards
!
!     History:
!     Version  Date      Comment.
!     5.3      04/10/01  Obsolete option for LW scattering removed.
!                                J. M. Edwards
!
      ! algorithmic options:

      ! spectral region
      INTEGER,PARAMETER:: ISOLIR_LW=IP_INFRA_RED

      ! method of angular integration
      INTEGER,PARAMETER:: I_ANGULAR_INTEGRATION_LW=IP_TWO_STREAM

      ! flag for properties in layers
      LOGICAL,PARAMETER:: L_LAYER_LW=.TRUE.

      ! flag for cloudy properties in layers
      LOGICAL,PARAMETER:: L_CLOUD_LAYER_LW=.TRUE.

      ! flag for corrections to 2-stream scheme
      LOGICAL,PARAMETER:: L_2_STREAM_CORRECT_LW=.FALSE.

      ! flag for rescaling of optical properties
      LOGICAL,PARAMETER:: L_RESCALE_LW=.TRUE.

      ! options invoking processes:

      LOGICAL,PARAMETER::  L_GAS_LW       =.TRUE. ! gaseous absorption
      LOGICAL,PARAMETER::  L_RAYLEIGH_LW  =.FALSE.! rayleigh scattering
      LOGICAL,PARAMETER::  L_CONTINUUM_LW =.TRUE. ! continuum absorption
      LOGICAL,PARAMETER::  L_CLOUD_LW     =.TRUE. ! clouds
      LOGICAL,PARAMETER::  L_DROP_LW      =.TRUE. ! droplets
      LOGICAL,PARAMETER::  L_ICE_LW       =.TRUE. ! ice crystals
      LOGICAL,PARAMETER::  L_AEROSOL_LW   =.TRUE. ! aerosols
      ! flag to use aerosols to determine ccn
      LOGICAL,PARAMETER::  L_AEROSOL_CCN_LW=.TRUE.

! LWFIX3A end
