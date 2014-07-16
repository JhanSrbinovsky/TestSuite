!     File to set heights for screen diagnostics and options for
!     diagnosis.
!
      REAL Z_OBS_TQ,Z_OBS_WIND
      PARAMETER (                                                       &
     & Z_OBS_TQ = 1.5                                                   &
                         ! Height of screen observations of temperature
!                        ! and humidity.
     &,Z_OBS_WIND = 10.0                                                &
                         ! Height of surface wind observations.
     &)
!
      INTEGER, Parameter :: IP_scrn_surf_sim = 0
!                           ! Diagnose the screen temperature using
!                           ! pure surface similarity theory
      INTEGER, Parameter :: IP_scrn_decpl_1=1
!                           ! Diagnose the screen temperature using 
!                           ! surface similarity theory, but allow 
!                           ! decoupling in very stable conditions
!                           ! based on the quasi-equilibrium radiative
!                           ! solution.
