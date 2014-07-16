!*L---------------COMDECK C_SLAB----------------------------------------
!PARAMETERS REQUIRED BY SLAB OCEAN MODEL AND NOT DEFINED ELSEWHERE
!
!CONRATIO IS THE RATIO OF THERMAL CONDUCTIVITIES OF ICE AND SNOW
!(DIMENSIONLESS)
!RHOCP IS THE VOLUMETRIC HEAT CAPACITY OF SEA WATER (J K-1 M-3)
!RHOICE IS THE DENSITY OF ICE (KG M-3)
!RHOSNOW IS THE DENSITY OF SNOW (KG M-3)
!NB ** RHOSNOW is also defined in the common deck C_SOILH, which
!cannot be used in the slab routines as it contains a duplicate
!definition of RHO_WATER, which is also found in C_DENSTY **
!** It should be noted that the value of RHOSNOW defined here matches
!   the value defined in C_SOIL_H, but differs from that currently
!   used in the ocean GCM (300 Kg m-3)
!
      REAL,PARAMETER:: CONRATIO=6.5656
      REAL,PARAMETER:: RHOCP=4.04E6
      REAL,PARAMETER:: RHOICE=900.0
      REAL,PARAMETER:: RHOSNOW=250.0
!
!*----------------------------------------------------------------------
