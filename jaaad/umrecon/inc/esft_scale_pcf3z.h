!     ------------------------------------------------------------------
!     Module setting types of ESFT scaling
!
      INTEGER                                                           &
     &    IP_SCALE_NULL                                                 &
!           No scaling at all
     &  , IP_SCALE_BAND                                                 &
!           Same scaling throughout band
     &  , IP_SCALE_TERM
!           Different scaling for each ESFT
!
      PARAMETER(                                                        &
     &    IP_SCALE_NULL=0                                               &
     &  , IP_SCALE_BAND=1                                               &
     &  , IP_SCALE_TERM=2                                               &
     &  )
!
!     ------------------------------------------------------------------
