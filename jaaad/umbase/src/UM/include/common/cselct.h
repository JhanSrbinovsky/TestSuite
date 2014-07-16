!----------------------------------------------------------------------
! comdeck: CSELCT
! Purpose: declares and stores information needed for flux selectiong
! History:
! version  date         change
! 4.5      21/09/98     New code
! Author: S.Spall
! 5.3      14/08/01     Allow atmosphere fields on B or C grids MJB
!----------------------------------------------------------------------
! common block:
      common / CSelect /                                                &
     &   l_B_grid,                                                      &
     &   l_winds_slt,   l_heat_slt,       l_moisture_slt,               &
     &   l_sea_ice_slt, l_references_slt, l_pressure_slt,               &
     &   l_windspd_slt

! debug logical for each selected flux
      logical l_B_grid,                                                 &
     &        l_winds_slt,   l_heat_slt,      l_moisture_slt,           &
     &        l_sea_ice_slt, l_references_slt,l_pressure_slt,           &
     &        l_windspd_slt
!----------------------------------------------------------------------
