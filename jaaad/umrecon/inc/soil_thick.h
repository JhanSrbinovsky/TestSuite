!---Soil layer thicknesses (m)
!---6 layers => CABLE else revert to MOSES
#ifdef CABLE_SOIL_LAYERS
   REAL,PARAMETER:: DZSOIL(6) =(/0.022, 0.058, 0.154, 0.409, 1.085, 2.872/)
#else
   REAL,PARAMETER:: DZSOIL(4) = (/0.100, 0.250, 0.650, 2.000 /)
#endif

