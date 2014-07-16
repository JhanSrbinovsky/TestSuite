!     ------------------------------------------------------------------
!     Module to set data for scaling functions.
!
      INTEGER                                                           &
     &    N_SCALE_VARIABLE(0: 3)
!           Number of variables in scaling functions
!
      DATA                                                              &
     &    N_SCALE_VARIABLE(IP_SCALE_POWER_LAW)/2/                       &
     &  , N_SCALE_VARIABLE(IP_SCALE_FNC_NULL)/0/                        &
     &  , N_SCALE_VARIABLE(IP_SCALE_POWER_QUAD)/3/                      &
     &  , N_SCALE_VARIABLE(IP_SCALE_DOPPLER_QUAD)/4/
!
!     ------------------------------------------------------------------
