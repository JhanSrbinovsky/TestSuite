!     ------------------------------------------------------------------
!     Module to set treatments of overlapping gaseous absorption.
!
      INTEGER                                                           &
     &    IP_OVERLAP_SINGLE                                             &
!           One species only
     &  , IP_OVERLAP_RANDOM                                             &
!           Random overlap
     &  , IP_OVERLAP_K_EQV                                              &
!           Equivalent extinction
     &  , IP_OVERLAP_SINGLE_INT
!           Interpolated treatment for principal species
!
      PARAMETER(                                                        &
     &    IP_OVERLAP_SINGLE=1                                           &
     &  , IP_OVERLAP_RANDOM=2                                           &
     &  , IP_OVERLAP_K_EQV=5                                            &
     &  , IP_OVERLAP_SINGLE_INT=6                                       &
     &  )
!
!     ------------------------------------------------------------------
