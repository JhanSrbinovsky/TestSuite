!     ------------------------------------------------------------------
!     Module to set dimensions for preprocessing routines.
!!
!!    Note on advised sizes. It may require too much space to set
!!    the sizes of all preprocessing arrays to  large values, and
!!    some variation in the sizes may be required between runs:
!!    in any event these criteria should be satisfied.
!!
!!       NPD_CASE > #bands considered * #p,t-pairs * #amounts
!!       NPD_GROUP > #bands considered * #p,t-pairs
!!       NPD_SET >  #bands considered * #p,t-pairs
!!                    * max(#species, #continua)
!!       NPD_SET_SCALE >  #p,t-pairs
!!       NPD_SET_GROUP >  max(#species, #continua)
!!
!!    The numbers used on the right-hand sides are those that will be
!!    used in the actual run, not the maxium values set below.
!
      INTEGER                                                           &
     &    NPD_SET                                                       &
!            Size allocated for sets of transmissions
     &  , NPD_SET_SCALE                                                 &
!           Size allocated for sets to scale at once
     &  , NPD_SET_GROUP                                                 &
!           Size allocated for sets in a group of transmission data
     &  , NPD_GROUP                                                     &
!           Size allocated for groups of data
     &  , NPD_PT                                                        &
!           Size allocated for p and T pairs
     &  , NPD_AMOUNT                                                    &
!           Size allocated for amounts
     &  , NPD_PARTIAL_P                                                 &
!           Size allocated for partial pressures
     &  , NPD_AMOUNT_CONTIN                                             &
!           Size allocated for continuum amounts
     &  , NPD_WAVENUMBER_TRANS                                          &
!           Size allocated for wavenumbers in a band
     &  , NPD_SOLAR_POINTS                                              &
!           Size allocated for the solar spectrum
     &  , NPD_SCATT_ANGLE                                               &
!           Size allocated for scattering angles
     &  , NPD_LOG_DERIV                                                 &
!           Size allocated for size of array of logarithmic
!           derivatives in the Mie scattering code
     &  , NPD_REFRACT                                                   &
!           Size allocated for (complex) refractive indices
     &  , NPD_WAVELENGTH_SCAT                                           &
!           Size allocated for scattering wavelengths
     &  , NPD_SIZE                                                      &
!           Size allocated for elements of a size distribution
     &  , NPD_MODE                                                      &
!           Size allocated for of modes in an analytical
!           size distribution
     &  , NPD_MIE_BLOCK                                                 &
!           Size allocated for blocks of Mie scattering data
     &  , NPD_THERMAL_ABSCISSA                                          &
!           Size allocated for number of points for integration with
!           respect to temperature
     &  , NPD_GENERAL_DATA                                              &
!           Size allocated for array allowed for general data
     &  , NPD_GAUSS_POINT                                               &
!           Size allocated for points of Gaussian integration
     &  , NPD_PHI_POINT                                                 &
!           Size allocated for integration in the phi-direction
     &  , NPD_BRDF_ORDER
!           Size allocated for orders of BRDFs
!
      PARAMETER(                                                        &
     &    NPD_SET=6900                                                  &
     &  , NPD_SET_SCALE=24                                              &
     &  , NPD_SET_GROUP=7                                               &
     &  , NPD_GROUP=240                                                 &
     &  , NPD_PT=24                                                     &
     &  , NPD_AMOUNT=330                                                &
     &  , NPD_PARTIAL_P=11                                              &
     &  , NPD_AMOUNT_CONTIN=NPD_AMOUNT/NPD_PARTIAL_P                    &
     &  , NPD_WAVENUMBER_TRANS=10000                                    &
     &  , NPD_SOLAR_POINTS=600                                          &
     &  , NPD_LOG_DERIV=250000                                          &
     &  , NPD_SCATT_ANGLE=200                                           &
     &  , NPD_REFRACT=600                                               &
     &  , NPD_WAVELENGTH_SCAT=500                                       &
     &  , NPD_SIZE=5000                                                 &
     &  , NPD_MODE=5                                                    &
     &  , NPD_MIE_BLOCK=400                                             &
     &  , NPD_THERMAL_ABSCISSA=30000                                    &
     &  , NPD_GENERAL_DATA=100                                          &
     &  , NPD_GAUSS_POINT=51                                            &
     &  , NPD_PHI_POINT=51                                              &
     &  , NPD_BRDF_ORDER=51                                             &
     &  )
!
!     ------------------------------------------------------------------
