! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ data module for switches/options concerned with the land surface.

MODULE land_surf_mod

  ! Description:
  !   Module containing runtime options/data used by the land surface.
  !
  ! Method:
  !   Switches and associated data values used by the land-surface scheme
  !   are defined here and assigned default values. These may be overridden
  !   by namelist input.
  !
  !   Any routine wishing to use these options may do so with the 'USE' 
  !   statement.
  !
  ! Current Code Owner: Adrian Lock
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !
  ! Declarations:

  IMPLICIT NONE

!===========================================================================
! integer options set from RUN_LAND namelist
!===========================================================================
  INTEGER :: FRAC_SNOW_SUBL_MELT = 0
!                                     ! switch for use of snow-cover
!                                     ! fraction in the calculation of
!                                     ! sublimation and melting
!                                     ! 0 = off
!                                     ! 1 = on
  INTEGER :: SOILHC_METHOD       = 2
!                                     ! switch for the calculation method
!                                     ! of soil thermal conductivity
  INTEGER :: ALL_TILES           = 0
!                                     ! switch for doing calculations
!                                     ! of tile properties on all tiles
!                                     ! for all gridpoints even when the
!                                     ! tile fraction is zero
!                                     ! (except for land ice).
!===========================================================================
! logical options set from RUN_LAND namelist
!===========================================================================
  LOGICAL :: L_VG_SOIL           = .FALSE.
!                                       ! switch for the use of
!                                       ! Van Genuchten soil parametrization
!===========================================================================
! real values set from RUN_LAND namelist
!===========================================================================

  REAL :: MASKD                  = 0.2
!                   ! masking depth for fractional snow calculation


! Define the RUN_LAND namelist
!-------------------------------
      NAMELIST/RUN_LAND/                                                &
     &  FRAC_SNOW_SUBL_MELT                                             &
     &, MASKD                                                           &
     &, SOILHC_METHOD                                                   &
     &, L_VG_SOIL                                                       &
     &, ALL_TILES

END MODULE land_surf_mod

