!MLO
      logical :: l_cable_lake
      real :: MLO_WATER_TEMP(land_points,ntiles,6)      ! Water temperature (K)
      real :: MLO_WATER_SALINITY(land_points,ntiles,6)  ! Water salinity (PSU)
      real :: MLO_WATER_UCUR(land_points,ntiles,6)      ! Water U current (m/s)
      real :: MLO_WATER_VCUR(land_points,ntiles,6)      ! Water V current (m/s)
      real :: MLO_SNOW_SURF_TEMP(land_points,ntiles)   ! Ice/snow surface surface (K)
      real :: MLO_SNOW_LEV1_TEMP(land_points,ntiles)   ! Snow Level 1 (if exists) temperature
      real :: MLO_ICE_LEV2_TEMP(land_points,ntiles)    ! Ice Level 2 (if exists) temperature
      real :: MLO_ICE_LEV3_TEMP(land_points,ntiles)    ! Ice level 3 (if exists) temperature
      real :: MLO_ICE_FRAC(land_points,ntiles)         ! Ice fraction (0-1)
      real :: MLO_ICE_DEPTH(land_points,ntiles)        ! Ice depth (m)
      real :: MLO_SNOW_DEPTH(land_points,ntiles)       ! Snow (on ice) depth (m)
      real :: MLO_ICE_ENERGY(land_points,ntiles)       ! Stored energy in ice (brine pockets) (Ws/m2)

