!     ------------------------------------------------------------------
!     Module to set components of clouds.
!
      INTEGER, PARAMETER :: IP_CLCMP_ST_WATER    = 1
!           Stratiform water droplets
      INTEGER, PARAMETER :: IP_CLCMP_ST_ICE      = 2
!           Stratiform ice crystals
      INTEGER, PARAMETER :: IP_CLCMP_CNV_WATER   = 3
!           Convective water droplets
      INTEGER, PARAMETER :: IP_CLCMP_CNV_ICE     = 4
!           Convective ice crystals

#if defined(XICE)
      INTEGER, PARAMETER :: IP_CLCMP_ST_ICE_COL  = 5
!           Columnar ice crystals in stratiform clouds
      INTEGER, PARAMETER :: IP_CLCMP_CNV_ICE_COL = 6
!           Columnar ice crystals in convective clouds
      INTEGER, PARAMETER :: IP_CLCMP_ST_ICE_ROS  = 7
!           Bullet rosette ice crystals in stratiform clouds
      INTEGER, PARAMETER :: IP_CLCMP_CNV_ICE_ROS = 8
!           Bullet rosette ice crystals in convective clouds
      INTEGER, PARAMETER :: IP_CLCMP_ST_ICE_PLT  = 9
!           Plate-like ice crystals in stratiform clouds
      INTEGER, PARAMETER :: IP_CLCMP_CNV_ICE_PLT = 10
!           Plate-like ice crystals in convective clouds
      INTEGER, PARAMETER :: IP_CLCMP_ST_ICE_PYC  = 11
!           Polycrystalline ice crystals in stratiform clouds
      INTEGER, PARAMETER :: IP_CLCMP_CNV_ICE_PYC = 12
!           Polycrystalline ice crystals in convective clouds
#endif
!     ------------------------------------------------------------------
