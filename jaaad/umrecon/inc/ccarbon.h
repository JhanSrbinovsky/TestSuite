!*L------------------COMDECK CCARBON------------------------------------
! Purpose: declares variables and parameters for the carbon cycle
! History:
! version  date         change
! 5.5      26/02/03     add M_CARBON. C Jones.
!----------------------------------------------------------------------
!carbon cycle and vegetation parameters
      REAL                                                              &
     & M_CO2                                                            &
                                  ! molecular weight of CO2
     &,M_AIR                                                            &
                                  ! molecular weight of dry air
     &,M_CARBON                                                         &
                                  ! molecular weight of carbon
     &,EPSILON                                                          &
                                  ! Ratio of molecular weights of water
!                                 !  and dry air.
     &,EPCO2                                                            &
                                  ! Ratio of molecular weights of CO2
!                                 !  and dry air.
     &,EPO2                                                             &
                                  ! Ratio of molecular weights of O2
!                                 !  and dry air.
     &,CO2CONV_A2O                                                      &
                                  ! conversion factor for atmos to
!                                 !  ocean passing of CO2 (mmr to ppmv)
     &,CO2CONV_O2A                ! conversion factor for ocean to
!                                 !  atmos passing of CO2 flux
!                                 !  (mol C/m2/yr to Kg CO2/m2/s)

      PARAMETER (M_AIR=28.966, EPCO2=1.5194, M_CO2=M_AIR*EPCO2,         &
     &           M_CARBON = 12.0, EPSILON = 0.62198, EPO2 = 1.106)

      PARAMETER (CO2CONV_A2O = M_AIR * 1E6 / M_CO2,                     &
     &           CO2CONV_O2A = M_CO2 * 1e-3 / (360.0 * 24.0 * 3600.0))
!*----------------------------------------------------------------------
