!     ------------------------------------------------------------------
!     Module to set values of astronomical constants.
!
      REAL  (real64) ::                                                 &
     &    T_EFFECTIVE_SOLAR                                             &
!           Effective solar temperature
     &  , D_EARTH_SUN                                                   &
!           Mean earth-sun distance
     &  , SUN_RADIUS                                                    &
!           Radius of sun
     &  , EARTH_RADIUS                                                  &
!           Radius of earth
     &  , ECCENTRICITY                                                  &
!           Eccentricity of earth's orbit
     &  , LENGTH_YEAR                                                   &
!           Number of days in year
     &  , DAY_PERIHELION
!           Time of annual perihelion
!
      PARAMETER(                                                        &
     &    T_EFFECTIVE_SOLAR=5.785E+03_real64                            &
     &  , SUN_RADIUS=6.96E+08_real64                                    &
     &  , D_EARTH_SUN=1.50E+11_real64                                   &
     &  , EARTH_RADIUS=6.371E+06_real64                                  &
     &  , ECCENTRICITY=1.67E-02_real64                                  &
     &  , LENGTH_YEAR=3.6525E+02_real64                                 &
     &  , DAY_PERIHELION=3.71E+00_real64                                &
     &  )
!
!     ------------------------------------------------------------------
