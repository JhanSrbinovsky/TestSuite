#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to convert differential IR radiances to actual ones.
!
! Purpose:
!   This subroutine receives differntial IR radiances or fluxes
!   and returns actual values.
!
! Method:
!   Striaghtforward.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE ADJUST_IR_RADIANCE(N_PROFILE, N_LAYER, N_VIEWING_LEVEL &
     &  , N_DIRECTION, I_ANGULAR_INTEGRATION, I_SPH_MODE                &
     &  , PLANCK_FLUX, PLANCK_RADIANCE                                  &
     &  , FLUX_DOWN, FLUX_UP, RADIANCE                                  &
     &  , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE                          &
     &  , ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL                      &
     &  )
!
!
!
!
!
      IMPLICIT NONE
!
!
!     Dummy array sizes
      INTEGER, INTENT(IN) ::                                            &
     &    ND_FLUX_PROFILE                                               &
!           Size allocated for atmospheric profiles where
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for atmospheric profiles for
!           quantities used in calculations of radiances
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_VIEWING_LEVEL                                              &
!           Size allocated for levels for radiances
     &  , ND_DIRECTION
!           Size allocated for directions
!
!
!     Include header files.
#include "c_kinds.h"
#include "c_pi.h"
#include "angular_integration_pcf3z.h"
#include "spectral_region_pcf3z.h"
#include "sph_mode_pcf3z.h"
!
!     Dummy arguments
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of atmospheric profiles
     &  , N_LAYER                                                       &
!           Number of atmospheric layers
     &  , N_DIRECTION                                                   &
!           Number of directions
     &  , N_VIEWING_LEVEL
!           Number of levels at which to calculate radiances
      INTEGER, INTENT(IN) ::                                            &
     &    I_ANGULAR_INTEGRATION                                         &
!           Angular integration scheme
     &  , I_SPH_MODE
!           Mode in which the spherical solver is used
      REAL  (Real64), INTENT(IN) ::                                     &
     &    PLANCK_FLUX(ND_FLUX_PROFILE, 0: ND_LAYER)                     &
!           Planckian fluxes
     &  , PLANCK_RADIANCE(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL)
!           Planckian radiances
!
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    FLUX_DOWN(ND_FLUX_PROFILE, 0: ND_LAYER)                       &
!           Downward fluxes
     &  , FLUX_UP(ND_FLUX_PROFILE, 0: ND_LAYER)                         &
!           Upward fluxes
     &  , RADIANCE(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL, ND_DIRECTION)
!           Radiances in specified directions
!
!
!     Local arguments
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , ID                                                            &
!           Loop variable
     &  , L
!           Loop variable
!
!
!

      IF ( (I_ANGULAR_INTEGRATION == IP_TWO_STREAM).OR.                 &
     &     (I_ANGULAR_INTEGRATION == IP_IR_GAUSS) ) THEN

!
        DO I=0, N_LAYER
          DO L=1, N_PROFILE
            FLUX_UP(L, I)=FLUX_UP(L, I)+PLANCK_FLUX(L, I)
            FLUX_DOWN(L, I)=FLUX_DOWN(L, I)+PLANCK_FLUX(L, I)
          ENDDO
        ENDDO
!
      ELSE IF (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC) THEN
!
!       Planckian radiances are always used with spherical harmonics,
!       even when calculating fluxes. The number of levels should
!       be set appropriately above.
        IF (I_SPH_MODE == IP_SPH_MODE_FLUX) THEN
          DO I=0, N_LAYER
            DO L=1, N_PROFILE
              FLUX_UP(L, I)=FLUX_UP(L, I)+PI*PLANCK_RADIANCE(L, I+1)
              FLUX_DOWN(L, I)=FLUX_DOWN(L, I)                           &
     &          +PI*PLANCK_RADIANCE(L, I+1)
            ENDDO
          ENDDO
        ELSE IF (I_SPH_MODE == IP_SPH_MODE_RAD) THEN
          DO ID=1, N_DIRECTION
            DO I=1, N_VIEWING_LEVEL
              DO L=1, N_PROFILE
                RADIANCE(L, I, ID)=RADIANCE(L, I, ID)                   &
     &            +PLANCK_RADIANCE(L, I)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE ADJUST_IR_RADIANCE
#endif
