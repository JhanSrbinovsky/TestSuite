! CSIZEOBS start
! The variables involved in dimensioning observational data arrays
! in the data assimilation section are stored in this comdeck.

      ! For ATMOSPHERE assimilations the values are computed in the
      ! initialisation routine INITAC and then passed to the main
      ! routine AC.

      INTEGER :: A_MAX_NO_OBS    !  No of observations in AC Obs files.
      INTEGER :: A_MAX_OBS_SIZE  !  No of obs values in AC Obs files.

      COMMON /CSIZEOBS/ A_MAX_NO_OBS, A_MAX_OBS_SIZE
! CSIZEOBS end
