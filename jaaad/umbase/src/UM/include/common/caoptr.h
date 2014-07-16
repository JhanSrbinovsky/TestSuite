! CAOPTR Holds address pointers for atmosphere-to-ocean coupling fields
! required by SWAP_A2O and SWAP_O2A.
! 5.5 28/02/03 Add pointers for river outflow on ATMOS grid
! later remove runoff and ocenpts pointers  C.Bunton
! 6.2 24/02/06 Add pointers for 10m windspeed            J.Gunson
! 6.2 29/11/05 Add pointers for ice velocities on OCEAN grid  A.B.Keen
!
! 6.2 24/02/06 Add pointers for ocean DMS flux to atmosphere    J.Gunson
! Pointers needed by SWAP_A2O and SWAP_O2A

      INTEGER :: JO_SST    ! Sea-surface temperature on ocean grid
      INTEGER :: JO_UCURR  ! Surface zonal current on ocean grid
      INTEGER :: JA_AICE   ! Fractional ice conc. on atmos grid
      INTEGER :: JA_AICEN  ! Category frac ice conc. on atmos grid
      INTEGER :: JO_AICEN  ! Category frac ice conc. on ocean grid
      COMMON /AOPTR/JO_SST,JO_UCURR,JA_AICE,JA_AICEN,JO_AICEN

      ! Pointers needed by SWAP_A2O

      INTEGER:: JA_TAUX       ! Surface x-windstress on atmos grid
      INTEGER:: JO_TAUX       ! Surface x-windstress on ocean grid
      INTEGER:: JA_TAUY       ! Surface y-windstress on atmos grid
      INTEGER:: JO_TAUY       ! Surface y-windstress on ocean grid
      INTEGER:: JA_WINDMIX    ! Windmixing power on atmos grid
      INTEGER:: JO_WINDMIX    ! Windmixing power on ocean grid
      INTEGER:: JA_W10        ! 10m windspeed
      INTEGER:: JO_W10
      INTEGER:: JA_DDDL1_D1   ! Dust dry dep lyr 1
      INTEGER:: JA_DDDL1_D2
      INTEGER:: JA_DDDL1_D3
      INTEGER:: JA_DDDL1_D4
      INTEGER:: JA_DDDL1_D5
      INTEGER:: JA_DDDL1_D6
      INTEGER:: JA_DDDL2_D1   ! Dust dry dep lyr 2
      INTEGER:: JA_DDDL2_D2
      INTEGER:: JA_DDDL2_D3
      INTEGER:: JA_DDDL2_D4
      INTEGER:: JA_DDDL2_D5
      INTEGER:: JA_DDDL2_D6
      INTEGER:: JA_DWDLS_D1   ! Dust wet dep lrg scl rain
      INTEGER:: JA_DWDLS_D2
      INTEGER:: JA_DWDLS_D3
      INTEGER:: JA_DWDLS_D4
      INTEGER:: JA_DWDLS_D5
      INTEGER:: JA_DWDLS_D6
      INTEGER:: JA_DWDCN_D1   ! Dust wet dep convctv rain
      INTEGER:: JA_DWDCN_D2
      INTEGER:: JA_DWDCN_D3
      INTEGER:: JA_DWDCN_D4
      INTEGER:: JA_DWDCN_D5
      INTEGER:: JA_DWDCN_D6
      INTEGER:: JO_DDEPD1     ! Dust deposition : ocean
      INTEGER:: JO_DDEPD2
      INTEGER:: JO_DDEPD3
      INTEGER:: JO_DDEPD4
      INTEGER:: JO_DDEPD5
      INTEGER:: JO_DDEPD6
      INTEGER:: JA_SOLAR      ! Net downward SW at surf on atmos grid
      INTEGER:: JA_BLUE       ! Net blueband SW at surf on atmos grid
      INTEGER:: JO_BLUE       ! Net blueband SW at surf on ocean grid
      INTEGER:: JA_EVAP       ! Net evaporation over sea on atmos grid
      INTEGER:: JA_LONGWAVE   ! Net downward LW at surf on atmos grid
      INTEGER:: JA_SENSIBLE  ! Sensible heat flux over sea on atmos grid
      INTEGER:: JO_HEATFLUX   ! Non penetrative heatflux on ocean grid
      INTEGER:: JA_LSSNOW     ! Large-scale snowfall rate on atmos grid
      INTEGER:: JA_CVSNOW     ! Convective snowfall rate on atmos grid
      INTEGER:: JA_LSRAIN     ! Large-scale rainfall rate on atmos grid
      INTEGER:: JA_CVRAIN     ! Convective rainfall rate on atmos grid
      INTEGER:: JO_PMINUSE    ! Precipitation-evaporation on ocean grid
      INTEGER:: JA_SLOWRUNOFF ! Slow (sub-surface) runoff on atmos grid
      INTEGER:: JA_FASTRUNOFF ! Fast (surface) runoff on atmos grid
      INTEGER:: JA_RIVEROUT   ! Total river outflow on atmos grid
      INTEGER:: JA_OCENTPTS   ! Ocean entry point index to atmos landpts
      INTEGER:: JO_RIVEROUT   ! Total river outflow on ocean grid
      INTEGER:: JA_co2        ! atmos level 1 co2 conc.
      INTEGER:: JO_co2
      INTEGER:: JA_co2flux    ! ocean co2 flux.
      INTEGER:: JO_co2flux
      INTEGER:: JA_dmsflux    ! ocean DMS flux.
      INTEGER:: JO_dmsflux
      INTEGER:: JO_SNOWFALL   ! Snowfall rate on ocean grid
      INTEGER:: JA_SUBLIM     ! Sublimation on atmos grid
      INTEGER:: JO_SUBLIM     ! Sublimation on ocean grid
      INTEGER:: JA_BOTMELT    ! Diffusive heat thro ice on atmos grid
      INTEGER:: JA_BOTMELTN   ! Diff heat thro ice (ncat)
      INTEGER:: JO_BOTMELT    ! Diffusive heat thro ice on ocean grid
      INTEGER:: JA_TOPMELT    ! Seaice top melting flux on atmos grid
      INTEGER:: JA_TOPMELTN   ! Seaice top melting flux (ncat)
      INTEGER:: JO_TOPMELT    ! Seaice top melting flux on ocean grid
      INTEGER:: JO_AICE       ! Sea ice concentration on ocean grid
      INTEGER:: JA_PRESS      ! Surface pressure on atmos grid
      INTEGER :: JC_U10
      INTEGER :: JC_V10
      INTEGER :: JC_CO2
 
      COMMON /A2OPTR/                                                   &
     &  JA_TAUX,JO_TAUX,JA_TAUY,JO_TAUY,JA_WINDMIX,JO_WINDMIX,          &
     &  JA_W10, JO_W10,                                                 &
     &  JA_DDDL1_D1,JA_DDDL1_D2,JA_DDDL1_D3,JA_DDDL1_D4,JA_DDDL1_D5,    &
     &  JA_DDDL1_D6,JA_DDDL2_D1,JA_DDDL2_D2,JA_DDDL2_D3,JA_DDDL2_D4,    &
     &  JA_DDDL2_D5,JA_DDDL2_D6,JA_DWDLS_D1,JA_DWDLS_D2,JA_DWDLS_D3,    &
     &  JA_DWDLS_D4,JA_DWDLS_D5,JA_DWDLS_D6,JA_DWDCN_D1,JA_DWDCN_D2,    &
     &  JA_DWDCN_D3,JA_DWDCN_D4,JA_DWDCN_D5,JA_DWDCN_D6,                &
     &  JO_DDEPD1,JO_DDEPD2,JO_DDEPD3,JO_DDEPD4,JO_DDEPD5,JO_DDEPD6,    &
     &  JA_SOLAR,JA_BLUE,JO_BLUE,JA_EVAP,JA_LONGWAVE,JA_SENSIBLE,       &
     &  JO_HEATFLUX,JA_LSSNOW,JA_CVSNOW,JA_LSRAIN,JA_CVRAIN,JO_PMINUSE, &
     &  JA_SLOWRUNOFF,JA_FASTRUNOFF,JA_OCENTPTS,JA_RIVEROUT,            &
     &  JO_RIVEROUT,JA_co2, JO_co2, JA_co2flux, JO_co2flux,             &
     &  JA_dmsflux, JO_dmsflux,                                         &
     &  JO_SNOWFALL,JA_SUBLIM,JO_SUBLIM,JA_BOTMELT,JO_BOTMELT,          &
     &  JA_TOPMELT,JO_TOPMELT,JA_BOTMELTN,JA_TOPMELTN,JO_AICE,JA_PRESS, &
     &  JC_U10, JC_V10, JC_CO2 

      INTEGER:: JO_TSTAR           ! Surface temperature on ocean grid
      INTEGER:: JA_TSTAR           ! Surface temperature on atmos grid
      INTEGER:: JA_UCURR           ! Surface zonal current on atmos grid
      INTEGER:: JO_VCURR           ! Surface merid current on ocean grid
      INTEGER:: JA_VCURR           ! Surface merid current on atmos grid
      INTEGER:: JO_UICE            ! Zonal seaice vel on ocean grid
      INTEGER:: JO_VICE            ! Merid seaice vel on ocean grid
      INTEGER:: JO_ICEDEPTH        ! Ice depth on ocean grid
      INTEGER:: JA_ICEDEPTH        ! Ice depth on atmos grid
      INTEGER:: JA_ICEDEPTHN       ! Ice depth on atmos grid (ncat)
      INTEGER:: JO_SNOWDEPTH       ! Snow depth on ocean grid
      INTEGER:: JA_SNOWDEPTH       ! Snow depth on atmos grid
      INTEGER:: JA_SNOWDEPTHN      ! Snow depth on atmos grid (ncat)

      COMMON /O2APTR/                                                   &
     &  JO_TSTAR,JA_TSTAR,JA_UCURR,JO_VCURR,JA_VCURR,                   &
     &  JO_UICE,JO_VICE,                                                &
     &  JO_ICEDEPTH,JA_ICEDEPTH,JO_SNOWDEPTH,JA_SNOWDEPTH,              &
     &  JA_ICEDEPTHN,JA_SNOWDEPTHN
! CAOPTR end
