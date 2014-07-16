!     ------------------------------------------------------------------
!     Module to set internal dimensions tied to algorithms,
!     mostly for clouds.
!
      INTEGER                                                           &
     &    NPD_CLOUD_COMPONENT                                           &
!           Number of components of clouds
     &  , NPD_CLOUD_TYPE                                                &
!           Number of permitted types of clouds.
     &  , NPD_CLOUD_REPRESENTATION                                      &
!           Number of permitted representations of clouds.
     &  , NPD_OVERLAP_COEFF                                             &
!           Number of overlap coefficients for clouds
     &  , NPD_SOURCE_COEFF                                              &
!           Number of coefficients for two-stream sources
     &  , NPD_REGION                                                    &
!           Number of regions in a layer
#if defined(EXPCLOP)
     &  , NPD_CLOUD_ELEM
!           Number of elements of cloud allowed in a layer
#endif
!
      PARAMETER(                                                        &
#if ! defined(XICE)
     &     NPD_CLOUD_COMPONENT=4                                        &
#else
     &     NPD_CLOUD_COMPONENT=12                                       &
#endif
     &  , NPD_CLOUD_TYPE=4                                              &
#if ! defined(XICE)
     &  , NPD_CLOUD_REPRESENTATION=4                                    &
#else
     &  , NPD_CLOUD_REPRESENTATION=5                                    &
#endif
     &  , NPD_OVERLAP_COEFF=18                                          &
     &  , NPD_SOURCE_COEFF=2                                            &
     &  , NPD_REGION=3                                                  &
#if defined(EXPCLOP)
     &  , NPD_CLOUD_ELEM=2                                              &
#endif
     &  )
!
!     ------------------------------------------------------------------
