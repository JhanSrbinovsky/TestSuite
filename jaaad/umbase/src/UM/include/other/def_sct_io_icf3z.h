!     ------------------------------------------------------------------
!     Module to set unit numbers for scattering programs.
!
      INTEGER                                                           &
     &     IU_MIE                                                       &
!           Unit number for main output
     &  , IU_PHASE                                                      &
!           Unit number for phase function
     &  , IU_AVERAGE                                                    &
!           Unit number for averaged file
     &  , IU_CLOUD_FIT
!           Unit number for fitted data
!
      PARAMETER(                                                        &
     &    IU_MIE=61                                                     &
     &  , IU_PHASE=62                                                   &
     &  , IU_AVERAGE=63                                                 &
     &  , IU_CLOUD_FIT=64                                               &
     &  )
!
!     ------------------------------------------------------------------
