!     ------------------------------------------------------------------
!     Module to define reference numbers for solvers.
!
      INTEGER                                                           &
     &    IP_SOLVER_PENTADIAGONAL                                       &
!           Pentadiagonal scheme
     &  , IP_SOLVER_MIX_11                                              &
!           Mixed column scheme using full endekadiagonal matrix
     &  , IP_SOLVER_MIX_APP_SCAT                                        &
!           Mixed column scheme with approximate scattering
     &  , IP_SOLVER_MIX_NET_APP_SCAT                                    &
!           Mixed column scheme for net fluxes
!           With approximate scattering
     &  , IP_SOLVER_MIX_DIRECT                                          &
!           Direct mixed column scheme for full fluxes
     &  , IP_SOLVER_MIX_DIRECT_NET                                      &
!           Direct mixed column scheme for net fluxes
     &  , IP_SOLVER_HOMOGEN_DIRECT                                      &
!           Direct solver for a homogeneous column
     &  , IP_SOLVER_TRIPLE                                              &
!           Direct solver for triple column
     &  , IP_SOLVER_TRIPLE_APP_SCAT                                     &
!           Direct solver for triple column approximating scattering
     &  , IP_solver_mix_direct_hogan                                    &
!           Direct mixed column scheme for full fluxes (modified 
!           for correct treatment of shadowing by Robin Hogan)
     &  , IP_solver_triple_hogan
!           Direct solver for triple column (modified for correct
!           treatment of shadowing by Robin Hogan)
!
      PARAMETER(                                                        &
     &    IP_SOLVER_PENTADIAGONAL=1                                     &
     &  , IP_SOLVER_MIX_11=6                                            &
     &  , IP_SOLVER_MIX_APP_SCAT=9                                      &
     &  , IP_SOLVER_MIX_NET_APP_SCAT=10                                 &
     &  , IP_SOLVER_MIX_DIRECT=11                                       &
     &  , IP_SOLVER_MIX_DIRECT_NET=12                                   &
     &  , IP_SOLVER_HOMOGEN_DIRECT=13                                   &
     &  , IP_SOLVER_TRIPLE=14                                           &
     &  , IP_SOLVER_TRIPLE_APP_SCAT=15                                  &
     &  , IP_solver_mix_direct_hogan=16                                 &
     &  , IP_solver_triple_hogan=17                                     &
     &  )
!
!     ------------------------------------------------------------------
