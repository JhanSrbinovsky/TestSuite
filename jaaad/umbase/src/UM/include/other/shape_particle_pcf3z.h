!     ------------------------------------------------------------------
!     Module to set the shape of particles.
!
      INTEGER                                                           &
     &    IP_SHAPE_NULL                                                 &
!           Unassigned shape
     &  , IP_SHAPE_SPHERE                                               &
!           Spherical particles
     &  , IP_SHAPE_HEXCYL                                               &
!           Hexagonal cylinders
     &  , IP_SHAPE_POLYCRYSTAL                                          &
!           Polycrystals
     &  , IP_SHAPE_PLATE                                                &
!           Plates
     &  , IP_SHAPE_ROSETTE                                              &
!           Rosettes
     &  , IP_SHAPE_COLUMN
!           Hexagonal columns (kept distinct from hexagonal cylinders
!           For convenience: this shape uses more recent relations
!           For the aspect ratio)
!
      PARAMETER(                                                        &
     &    IP_SHAPE_NULL=0                                               &
     &  , IP_SHAPE_SPHERE=1                                             &
     &  , IP_SHAPE_HEXCYL=2                                             &
     &  , IP_SHAPE_POLYCRYSTAL=3                                        &
     &  , IP_SHAPE_PLATE=4                                              &
     &  , IP_SHAPE_ROSETTE=5                                            &
     &  , IP_SHAPE_COLUMN=6                                             &
     &  )
!
!     ------------------------------------------------------------------
