!     ------------------------------------------------------------------
!     Module to set representation of optical properties of
!     layers in terms of levels.
!!    If weight_upper is set to 1 the optical properties for the layer
!!    are taken as those at the upper boundary: if it is 0.5, the
!!    optical properties of the layer are the mean of the boundary
!!    values. This freedom is permitted only for clouds and aerosols:
!!    it would not make sense for well-mixed gases.
!
      REAL  (real64) ::                                                 &
     &    WEIGHT_UPPER_CLOUD                                            &
!           Weighting for clouds
     &  , WEIGHT_UPPER_AEROSOL
!           Weighting for aerosols
!
      PARAMETER(                                                        &
     &    WEIGHT_UPPER_CLOUD=0.5E+00_real64                             &
     &  , WEIGHT_UPPER_AEROSOL=1.0E+00_real64                           &
     &  )
!
!     ------------------------------------------------------------------
