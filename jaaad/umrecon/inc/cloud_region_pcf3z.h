!     ------------------------------------------------------------------
!     Module to define reference numbers for regions of clouds.
!
      INTEGER                                                           &
     &    IP_REGION_CLEAR                                               &
!           Reference number for clear-sky region
     &  , IP_REGION_STRAT                                               &
!           Reference number for stratiform cloudy region
     &  , IP_REGION_CONV
!           Reference number for convective cloudy region
!
      PARAMETER(                                                        &
     &    IP_REGION_CLEAR=1                                             &
     &  , IP_REGION_STRAT=2                                             &
     &  , IP_REGION_CONV=3                                              &
     &  )
!
!     ------------------------------------------------------------------
