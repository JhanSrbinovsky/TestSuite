! CASPTR Holds address pointers for slab ocean model fields

      ! Pointers needed by SLABCNTL
      INTEGER:: JS_SOLARIN  ! Net downward shortwave heat flux
      INTEGER:: JS_EVAP     ! Surface evaporation
      INTEGER:: JS_LONGWAVE ! Net downward longwave heat flux
      INTEGER:: JS_SENSIBLE ! Sensible heat flux
      INTEGER:: JS_HEATCONV ! Heat convergence rate
      INTEGER:: JS_SNOWLS   ! Large scale snowfall rate
      INTEGER:: JS_SNOWCONV ! Convective snowfall rate
      INTEGER:: JS_SUBLIMZ  ! Rate of sublimation
      INTEGER:: JS_TOPMELTZ ! Rate of melting of snow (or ice)
      INTEGER:: JS_BOTMELTZ ! Diffusive heat flux through ice
      INTEGER:: JS_BLUEIN   ! Net surface SW over sea; Band 1
      INTEGER:: JS_USEA     ! Zonal component of surface current
      INTEGER:: JS_VSEA     ! Meridional component of surface current
      INTEGER:: JS_WSX      ! Zonal component of surface wind stress
      INTEGER:: JS_WSY      ! Meridional compt of surface wind stress

      COMMON /ASPTR/                                                    &
     &  JS_SOLARIN,JS_EVAP,JS_LONGWAVE,JS_SENSIBLE,JS_HEATCONV,         &
     &  JS_SNOWLS,JS_SNOWCONV,JS_SUBLIMZ,JS_TOPMELTZ,JS_BOTMELTZ,       &
     &  JS_BLUEIN,JS_USEA,JS_VSEA,JS_WSX,JS_WSY
! CASPTR end
