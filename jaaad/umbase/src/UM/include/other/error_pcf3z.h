!     ------------------------------------------------------------------
!     Module to set error flags in the radiation code.
!
      INTEGER                                                           &
     &    I_NORMAL                                                      &
!           Error free condition
     &  , I_ERR_FATAL                                                   &
!           Fatal error: immediate return
     &  , I_ABORT_CALCULATION                                           &
!           Calculation aborted
     &  , I_MISSING_DATA                                                &
!           Missing data error: conditional
     &  , I_ERR_IO                                                      &
!           I/O error
     &  , I_ERR_RANGE                                                   &
!           Interpolation range error
     &  , I_ERR_EXIST
!           Existence error
!
      PARAMETER(                                                        &
     &    I_NORMAL=0                                                    &
     &  , I_ERR_FATAL=1                                                 &
     &  , I_ABORT_CALCULATION=2                                         &
     &  , I_MISSING_DATA=3                                              &
     &  , I_ERR_IO=4                                                    &
     &  , I_ERR_RANGE=5                                                 &
     &  , I_ERR_EXIST=6                                                 &
     &  )
!
!     ------------------------------------------------------------------
