!     -------------------------------------------------------------------
!     Module to set permitted types of scattering algorithm.
!
      INTEGER                                                           &
     &    NPD_ALGORITHM_SCATTER                                         &
!           Number of algorithms
     &  , IP_ALGORITHM_MIE                                              &
!           Identification number for full mie calculation
     &  , IP_ALGORITHM_ADT                                              &
!           Identification number for adt
     &  , IP_TAKANO_LIOU_HEXCYL                                         &
!           Takano & liou's treatment of hexagonal cylinders
     &  , IP_ENHANCED_ADT
!           Takano & liou's treatment of hexagonal cylinders
!
      PARAMETER(                                                        &
     &    NPD_ALGORITHM_SCATTER=4                                       &
     &  , IP_ALGORITHM_MIE=1                                            &
     &  , IP_ALGORITHM_ADT=2                                            &
     &  , IP_TAKANO_LIOU_HEXCYL=3                                       &
     &  , IP_ENHANCED_ADT=4                                             &
     &  )
!
!     -------------------------------------------------------------------
