!     ------------------------------------------------------------------
!     Module defining types of intermediate files.
!
      INTEGER                                                           &
     &    IT_FILE_LINE_TRANS_FORM                                       &
!           Formatted line transmission data
     &  , IT_FILE_CONT_TRANS_FORM                                       &
!           Formatted continuum transmission
     &  , IT_FILE_LINE_TRANS                                            &
!           File of line transmission data
     &  , IT_FILE_CONT_TRANS                                            &
!           File of continuum transmission data
     &  , IT_FILE_LINE_TRANS_SUMMARY                                    &
!           File of line transmission data
     &  , IT_FILE_CONT_TRANS_SUMMARY                                    &
!           File of continuum transmission data
     &  , IT_FILE_LINE_FIT                                              &
!           File of fits to line data
     &  , IT_FILE_CONT_FIT                                              &
!           File of fits to continuum data
     &  , IT_FILE_MIE_DRY                                               &
!           File of dry mie data
     &  , IT_FILE_MIE_HUMID                                             &
!           File of moist mie data
     &  , IT_FILE_AVE_MIE_DRY                                           &
!           File of dry averaged mie data
     &  , IT_FILE_AVE_MIE_HUMID                                         &
!           File of moist averaged mie data
     &  , IT_FILE_CLOUD_FIT                                             &
!           File of fits to cloud data
     &  , IT_FILE_CLOUD_OBS                                             &
!           File of observational cloud data
     &  , IT_FILE_ADT_DRY                                               &
!           File of dry scattering properties from ADT
     &  , IT_FILE_ADT_HUMID                                             &
!           File of dry scattering properties from ADT
     &  , IT_FILE_NSICE_AVE                                             &
!           File of averaged scattering properties for
!           non-spherical ice
     &  , IT_FILE_PHF_MIE_DRY                                           &
!           File of Mie scattering data with more than one moment
!           of the phase function
     &  , IT_FILE_PHF_MIE_HUMID                                         &
!           File of Mie scattering data with more than one moment
!           of the phase function
     &  , IT_FILE_AVE_PHF_MIE_DRY                                       &
!           File of averaged dry single scattering data including
!           higher moments of the phase function
     &  , IT_FILE_AVE_PHF_MIE_HUMID                                     &
!           File of averaged moist single scattering data including
!           higher moments of the phase function
     &  , IT_FILE_CLOUD_FIT_PHF
!           File of fitted single scattering data including
!           higher moments of the phase function
!
      PARAMETER(                                                        &
     &    IT_FILE_LINE_TRANS_FORM=2                                     &
     &  , IT_FILE_CONT_TRANS_FORM=3                                     &
     &  , IT_FILE_LINE_FIT=4                                            &
     &  , IT_FILE_CONT_FIT=5                                            &
     &  , IT_FILE_MIE_DRY=6                                             &
     &  , IT_FILE_MIE_HUMID=7                                           &
     &  , IT_FILE_AVE_MIE_DRY=8                                         &
     &  , IT_FILE_AVE_MIE_HUMID=9                                       &
     &  , IT_FILE_CLOUD_FIT=10                                          &
     &  , IT_FILE_CLOUD_OBS=11                                          &
     &  , IT_FILE_LINE_TRANS=12                                         &
     &  , IT_FILE_CONT_TRANS=13                                         &
     &  , IT_FILE_LINE_TRANS_SUMMARY=14                                 &
     &  , IT_FILE_CONT_TRANS_SUMMARY=15                                 &
     &  , IT_FILE_ADT_DRY=16                                            &
     &  , IT_FILE_NSICE_AVE=17                                          &
     &  , IT_FILE_PHF_MIE_DRY=18                                        &
     &  , IT_FILE_PHF_MIE_HUMID=19                                      &
     &  , IT_FILE_ADT_HUMID=20                                          &
     &  , IT_FILE_AVE_PHF_MIE_DRY=21                                    &
     &  , IT_FILE_AVE_PHF_MIE_HUMID=22                                  &
     &  , IT_FILE_CLOUD_FIT_PHF=24                                      &
     &  )
!
!     ------------------------------------------------------------------
