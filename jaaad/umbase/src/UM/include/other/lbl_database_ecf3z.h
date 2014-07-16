!     ------------------------------------------------------------------
!     Module to set identification numbers for the LBL databases.
!
      INTEGER                                                           &
     &    IP_LBL_DATA_UNASSIGNED                                        &
!           Unassigned line atlas
     &  , IP_LBL_DATA_HITRAN92                                          &
!           Line atlas HITRAN92
     &  , IP_LBL_DATA_HITRAN91                                          &
!           Line atlas HITRAN91
     &  , IP_LBL_DATA_HITRAN86                                          &
!           Line atlas HITRAN86
     &  , IP_LBL_DATA_HITRAN96
!           Line atlas HITRAN86
!
      PARAMETER(                                                        &
     &    IP_LBL_DATA_UNASSIGNED=0                                      &
     &  , IP_LBL_DATA_HITRAN92=1                                        &
     &  , IP_LBL_DATA_HITRAN91=2                                        &
     &  , IP_LBL_DATA_HITRAN86=3                                        &
     &  , IP_LBL_DATA_HITRAN96=4                                        &
     &  )
!
!     ------------------------------------------------------------------
