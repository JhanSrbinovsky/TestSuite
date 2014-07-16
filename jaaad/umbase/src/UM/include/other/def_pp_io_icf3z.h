!     ------------------------------------------------------------------
!     Module to set unit numbers for preprocessing.
!
      INTEGER                                                           &
     &    IU_WRSPC                                                      &
!           Unit number to write spectral file
     &  , IU_TRANSMISSION                                               &
!           Unit number for transmission file
     &  , IU_TRANSMISSION_SUMMARY                                       &
!           Unit numer for summary of transmission data
     &  , IU_C_TRANSMISSION                                             &
!           Unit number for continuum file
     &  , IU_C_TRANSMISSION_SUMMARY                                     &
!           Unit number for summary of continuum data
     &  , IU_TRANSMISSION_FORM                                          &
!           Unit number for formatted transmissions
     &  , IU_SOLAR                                                      &
!           Unit number for solar spectrum
     &  , IU_ESFT_OUT                                                   &
!           Unit number for ESFT output
     &  , IU_ESFT
!           Unit number for ESFT input
!
      PARAMETER(                                                        &
     &    IU_WRSPC=30                                                   &
     &  , IU_TRANSMISSION=31                                            &
     &  , IU_TRANSMISSION_SUMMARY=34                                    &
     &  , IU_TRANSMISSION_FORM=35                                       &
     &  , IU_C_TRANSMISSION=36                                          &
     &  , IU_C_TRANSMISSION_SUMMARY=37                                  &
     &  , IU_ESFT_OUT=32                                                &
     &  , IU_ESFT=33                                                    &
     &  , IU_SOLAR=40                                                   &
     &  )
!
!     ------------------------------------------------------------------
