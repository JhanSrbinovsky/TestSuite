! SWFIX3A defines options to the edwards-slingo radiation code
! fixed in the unified model. options for shortwave calculations.

      ! algorithmic options:

      INTEGER,PARAMETER:: ISOLIR_SW=IP_SOLAR ! spectral region

      ! method of angular integration
      INTEGER,PARAMETER:: I_ANGULAR_INTEGRATION_SW=IP_TWO_STREAM

      ! treatment of scattering
      INTEGER,PARAMETER:: I_SCATTER_METHOD_SW=IP_SCATTER_FULL

      ! flag for properties in layers
      LOGICAL,PARAMETER:: L_LAYER_SW=.TRUE.

      ! flag for cloudy properties in layers
      LOGICAL,PARAMETER:: L_CLOUD_LAYER_SW=.TRUE.

      ! flag for corrections to 2-stream scheme
      LOGICAL,PARAMETER:: L_2_STREAM_CORRECT_SW=.FALSE.

      ! flag for rescaling of optical properties
      LOGICAL,PARAMETER:: L_RESCALE_SW=.TRUE.

!     OPTIONS INVOKING PROCESSES:
!
!     flag for gaseous absorption
      LOGICAL,PARAMETER:: L_GAS_SW         =.TRUE. ! gaseous absorption

!     flag for rayleigh scattering
      LOGICAL,PARAMETER:: L_RAYLEIGH_SW    =.TRUE. ! rayleigh scattering

!     flag for continuum absorption
      LOGICAL,PARAMETER:: L_CONTINUUM_SW   =.TRUE. ! continuum absorption

!     flag for clouds
      LOGICAL,PARAMETER:: L_CLOUD_SW       =.TRUE. ! clouds

!     flag for droplets
      LOGICAL,PARAMETER:: L_DROP_SW        =.TRUE. ! droplets

!     flag for ice crystals
      LOGICAL,PARAMETER:: L_ICE_SW         =.TRUE. ! ice crystals

!     flag for aerosols
      LOGICAL,PARAMETER:: L_AEROSOL_SW     =.TRUE. ! aerosols

      ! flag aerosols to determine ccn
      LOGICAL,PARAMETER:: L_AEROSOL_CCN_SW =.TRUE.
! SWFIX3A end
