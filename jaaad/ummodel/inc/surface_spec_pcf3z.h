!     ------------------------------------------------------------------
!     Module to set permitted methods of specifying the characteristics
!     of the surface.
!
!
!     Ways of specifiying the surface
!
      INTEGER                                                           &
     &    IP_SURFACE_SPECIFIED                                          &
!           Properties specified by surface type
     &  , IP_SURFACE_INTERNAL                                           &
!           Properties passed into code
     &  , IP_SURFACE_POLYNOMIAL                                         &
!           Direct albedo fitted as polynomial
     &  , IP_SURFACE_PAYNE                                              &
!           Fit in the functional form used by Payne
     &  , IP_SURFACE_LAMBERTIAN                                         &
!           BRDF represented by a Lambertian
     &  , IP_SURFACE_ROUJEAN                                            &
!           BRDF represented by a Roujean's basis
     &  , IP_SURFACE_LOMMEL_SEELIGER_AXI
!           BRDF represented by an axisymmetric Lommel-Seeliger
!           function
!
      PARAMETER(                                                        &
     &    IP_SURFACE_SPECIFIED=1                                        &
     &  , IP_SURFACE_INTERNAL=2                                         &
     &  , IP_SURFACE_POLYNOMIAL=3                                       &
     &  , IP_SURFACE_PAYNE=4                                            &
     &  , IP_SURFACE_LAMBERTIAN=5                                       &
     &  , IP_SURFACE_ROUJEAN=6                                          &
     &  , IP_SURFACE_LOMMEL_SEELIGER_AXI=7                              &
     &  )
!
!
!     Pointers to specific components of arrays
!
      INTEGER                                                           &
     &    IP_SURF_ALB_DIFF                                              &
!           Pointer to diffuse surface albedo
     &  , IP_SURF_ALB_DIR
!           Pointer to direct surface albedo
!
      PARAMETER(                                                        &
     &    IP_SURF_ALB_DIFF=1                                            &
     &  , IP_SURF_ALB_DIR=2                                             &
     &  )
!
!     ------------------------------------------------------------------
