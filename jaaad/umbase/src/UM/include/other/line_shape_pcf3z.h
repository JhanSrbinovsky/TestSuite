!     ------------------------------------------------------------------
!     Module to set the various line shapes used.
!
      INTEGER                                                           &
     &     IP_LINE_SHAPE_UNASSIGNED                                     &
!           Unassigned line shape
     &  , IP_LINE_SHAPE_VOIGT                                           &
!           Voigt line shape
     &  , IP_LINE_SHAPE_LORENTZ                                         &
!           Lorentzian line shape
     &  , IP_LINE_SHAPE_VANVHUB                                         &
!           Van vleck & huber line shape
     &  , IP_LINE_SHAPE_NEWSH                                           &
!           New line shape
     &  , IP_LINE_SHAPE_DOPPLER                                         &
!           Doppler line shape
     &  , IP_LINE_SHAPE_VOIGTCO2                                        &
!           Sub-voigt co2 line shape
     &  , IP_LINE_SHAPE_XSECN
!           Cross-sectional line shape
!
      PARAMETER(                                                        &
     &    IP_LINE_SHAPE_UNASSIGNED=0                                    &
     &  , IP_LINE_SHAPE_VOIGT=1                                         &
     &  , IP_LINE_SHAPE_LORENTZ=2                                       &
     &  , IP_LINE_SHAPE_VANVHUB=3                                       &
     &  , IP_LINE_SHAPE_NEWSH=4                                         &
     &  , IP_LINE_SHAPE_DOPPLER=5                                       &
     &  , IP_LINE_SHAPE_VOIGTCO2=6                                      &
     &  , IP_LINE_SHAPE_XSECN=7                                         &
     &  )
!
!     ------------------------------------------------------------------
