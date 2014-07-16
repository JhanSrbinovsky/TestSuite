!     ------------------------------------------------------------------
!     Module to set the methods of treating scattering.
!
      INTEGER                                                           &
     &    IP_SCATTER_FULL                                               &
!           Full treatment of scattering
     &  , IP_NO_SCATTER_ABS                                             &
!           Scattering ignored completely.
     &  , IP_NO_SCATTER_EXT
!           Scattering treated as absorption
!
      PARAMETER(                                                        &
     &    IP_SCATTER_FULL=1                                             &
     &  , IP_NO_SCATTER_ABS=2                                           &
     &  , IP_NO_SCATTER_EXT=3                                           &
     &  )
!
!     ------------------------------------------------------------------
