!     ------------------------------------------------------------------
!     Module to set scattering parameters:
!
      INTEGER                                                           &
     &    NPD_SCATTER_TYPE                                              &
!           Number of tyeps of scatterer
     &  , IP_TYPE_UNASSIGNED                                            &
!           Index of unassigned scatterer
     &  , IP_TYPE_AEROSOL                                               &
!           Index of aerosol
     &  , IP_TYPE_DROPLET                                               &
!           Index of droplets
     &  , IP_TYPE_ICE
!           Index of ice crystals
!
      PARAMETER(                                                        &
     &    NPD_SCATTER_TYPE=3                                            &
     &  , IP_TYPE_UNASSIGNED=0                                          &
     &  , IP_TYPE_AEROSOL=1                                             &
     &  , IP_TYPE_DROPLET=2                                             &
     &  , IP_TYPE_ICE=3                                                 &
     &  )
!
!     ------------------------------------------------------------------
