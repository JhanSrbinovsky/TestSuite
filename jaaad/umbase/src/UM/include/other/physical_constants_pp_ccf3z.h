!     ------------------------------------------------------------------
!     Module setting physical constants.
!
      REAL  (real64) ::                                                 &
     &    H_PLANCK                                                      &
!           Planck's constant
     &  , C_LIGHT                                                       &
!           Speed of light
     &  , K_BOLTZMANN                                                   &
!           Boltzmann's constant
     &  , RHO_N                                                         &
!           Depolarizing factor
     &  , N_AVOGADRO                                                    &
!           Avogadro's number
     &  , RHO_AIR_STP
!           Density of air at stp
!
      PARAMETER(                                                        &
     &    H_PLANCK=6.626176E-34_real64                                  &
     &  , C_LIGHT=2.9979245E+08_real64                                  &
     &  , K_BOLTZMANN=1.380662E-23_real64                               &
     &  , RHO_N=2.79E-02_real64                                         &
     &  , N_AVOGADRO=6.022045E+23_real64                                &
     &  , RHO_AIR_STP=1.293125E+00_real64                               &
     &  )
!
!     ------------------------------------------------------------------
