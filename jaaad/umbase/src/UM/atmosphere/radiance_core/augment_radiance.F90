#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to increment a radiances or fluxes.
!
! Method:
!       The arrays holding the summed fluxes or radiances are
!       incremented by a weighted sum of the variables suffixed
!       with _INCR. Arguments specify which arrays are to be
!       incremented.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE AUGMENT_RADIANCE(N_PROFILE, N_LAYER                    &
     &  , I_ANGULAR_INTEGRATION, I_SPH_MODE                             &
     &  , N_VIEWING_LEVEL, N_DIRECTION                                  &
     &  , ISOLIR, L_CLEAR                                               &
     &  , L_INITIAL, WEIGHT_INCR                                        &
     &  , L_BLUE_FLUX_SURF, WEIGHT_BLUE                                 &
!                     Actual radiances
     &  , FLUX_DIRECT, FLUX_DOWN, FLUX_UP                               &
     &  , FLUX_DIRECT_BLUE_SURF                                         &
     &  , FLUX_DOWN_BLUE_SURF, FLUX_UP_BLUE_SURF                        &
     &  , I_DIRECT, RADIANCE, PHOTOLYSIS                                &
     &  , FLUX_DIRECT_CLEAR, FLUX_DOWN_CLEAR, FLUX_UP_CLEAR             &
!                     Increments to radiances
     &  , FLUX_DIRECT_INCR, FLUX_TOTAL_INCR                             &
     &  , I_DIRECT_INCR, RADIANCE_INCR, PHOTOLYSIS_INCR                 &
     &  , FLUX_DIRECT_INCR_CLEAR, FLUX_TOTAL_INCR_CLEAR                 &
!                     Dimensions
     &  , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE            &
     &  , ND_LAYER, ND_VIEWING_LEVEL, ND_DIRECTION                      &
     &  )
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_FLUX_PROFILE                                               &
!           Size allocated for points where fluxes are calculated
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for points where radiances are calculated
     &  , ND_J_PROFILE                                                  &
!           Size allocated for points where photolysis is calculated
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ND_VIEWING_LEVEL                                              &
!           Size allocated for levels where radiances are calculated
     &  , ND_DIRECTION
!           Size allocated for viewing directions
!
!     Include header files.
#include "spectral_region_pcf3z.h"
#include "angular_integration_pcf3z.h"
#include "sph_mode_pcf3z.h"
#include "sph_algorithm_pcf3z.h"
#include "c_kinds.h"
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_VIEWING_LEVEL                                               &
!           Number of levels where the radiance is calculated
     &  , N_DIRECTION
!           Number of viewing directions
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR                                                        &
!           Spectral region
     &  , I_SPH_MODE                                                    &
!           Mode in which spherical harmonics are used
     &  , I_ANGULAR_INTEGRATION
!           Treatment of angular integration
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLEAR                                                       &
!           Clear fluxes calculated
     &  , L_INITIAL
!           Logical to perform initialization instead of incrementing
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    WEIGHT_INCR
!           Weight to apply to incrementing fluxes
!
!                     Increments to Fluxes
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FLUX_DIRECT_INCR(ND_FLUX_PROFILE, 0: ND_LAYER)                &
!           Increment to direct flux
     &  , FLUX_TOTAL_INCR(ND_FLUX_PROFILE, 2*ND_LAYER+2)                &
!           Increment to total flux
     &  , FLUX_DIRECT_INCR_CLEAR(ND_FLUX_PROFILE, 0: ND_LAYER)          &
!           Increment to clear direct flux
     &  , FLUX_TOTAL_INCR_CLEAR(ND_FLUX_PROFILE, 2*ND_LAYER+2)
!           Increment to clear total flux
!                     Increments to Radiances
      REAL  (Real64), INTENT(IN) ::                                     &
     &    I_DIRECT_INCR(ND_RADIANCE_PROFILE, 0: ND_LAYER)               &
!           Increments to the solar irradiance
     &  , RADIANCE_INCR(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL           &
     &      , ND_DIRECTION)
!           Increments to the radiance
!                     Increments to Rates of photolysis
      REAL  (Real64), INTENT(IN) ::                                     &
     &    PHOTOLYSIS_INCR(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Increments to the rates of photolysis
!
!                     Total Fluxes
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    FLUX_DIRECT(ND_FLUX_PROFILE, 0: ND_LAYER)                     &
!           Direct flux
     &  , FLUX_DOWN(ND_FLUX_PROFILE, 0: ND_LAYER)                       &
!           Total downward flux
     &  , FLUX_UP(ND_FLUX_PROFILE, 0: ND_LAYER)                         &
!           Upward flux
     &  , FLUX_DIRECT_CLEAR(ND_FLUX_PROFILE, 0: ND_LAYER)               &
!           Clear direct flux
     &  , FLUX_DOWN_CLEAR(ND_FLUX_PROFILE, 0: ND_LAYER)                 &
!           Clear total downward flux
     &  , FLUX_UP_CLEAR(ND_FLUX_PROFILE, 0: ND_LAYER)
!           Clear upward flux
!                     Total Radiances
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    I_DIRECT(ND_RADIANCE_PROFILE, 0: ND_LAYER)                    &
!           Solar irradiance
     &  , RADIANCE(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL                &
     &      , ND_DIRECTION)
!           Radiance
!                       Rates of photolysis
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    PHOTOLYSIS(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Rates of photolysis
!
!                        Special Diagnostics:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_BLUE_FLUX_SURF
!           Flag to calculate blue fluxes at the surface
      REAL  (Real64), INTENT(IN) ::                                     &
     &    WEIGHT_BLUE
!           Weights for blue fluxes in this band
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    FLUX_DIRECT_BLUE_SURF(ND_FLUX_PROFILE)                        &
!           Direct blue flux at the surface
     &  , FLUX_DOWN_BLUE_SURF(ND_FLUX_PROFILE)                          &
!           Total downward blue flux at the surface
     &  , FLUX_UP_BLUE_SURF(ND_FLUX_PROFILE)
!           Upward blue flux at the surface
!
!
!
!     Local arguments.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , K
!           Loop variable
!
!
!
      IF (.NOT.L_INITIAL) THEN
!
!       Most commonly, this routine will be called to increment
!       rather than to initialize fluxes.
!
        IF ( (I_ANGULAR_INTEGRATION == IP_TWO_STREAM).OR.               &
     &       (I_ANGULAR_INTEGRATION == IP_IR_GAUSS).OR.                 &
     &     ( (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC).AND.      &
     &       (I_SPH_MODE == IP_SPH_MODE_FLUX) ) ) THEN
!
!         Increment the actual fluxes.
          IF (ISOLIR == IP_SOLAR) THEN
            DO I=0, N_LAYER
              DO L=1, N_PROFILE
                FLUX_DIRECT(L, I)=FLUX_DIRECT(L, I)                     &
     &            +WEIGHT_INCR*FLUX_DIRECT_INCR(L, I)
              ENDDO
            ENDDO
            IF (L_BLUE_FLUX_SURF) THEN
              DO L=1, N_PROFILE
                FLUX_UP_BLUE_SURF(L)=FLUX_UP_BLUE_SURF(L)               &
     &            +WEIGHT_BLUE*FLUX_TOTAL_INCR(L, 2*N_LAYER+1)
                FLUX_DOWN_BLUE_SURF(L)=FLUX_DOWN_BLUE_SURF(L)           &
     &            +WEIGHT_BLUE*FLUX_TOTAL_INCR(L, 2*N_LAYER+2)
              ENDDO
              IF (ISOLIR == IP_SOLAR) THEN
                DO L=1, N_PROFILE
                  FLUX_DIRECT_BLUE_SURF(L)=FLUX_DIRECT_BLUE_SURF(L)     &
     &              +WEIGHT_BLUE*FLUX_DIRECT_INCR(L, N_LAYER)
                ENDDO
              ENDIF
            ENDIF
          ENDIF
          DO I=0, N_LAYER
            DO L=1, N_PROFILE
              FLUX_UP(L, I)=FLUX_UP(L, I)                               &
     &          +WEIGHT_INCR*FLUX_TOTAL_INCR(L, 2*I+1)
              FLUX_DOWN(L, I)=FLUX_DOWN(L, I)                           &
     &          +WEIGHT_INCR*FLUX_TOTAL_INCR(L, 2*I+2)
            ENDDO
          ENDDO
!
          IF (L_CLEAR) THEN
            IF (ISOLIR == IP_SOLAR) THEN
              DO I=0, N_LAYER
                DO L=1, N_PROFILE
                  FLUX_DIRECT_CLEAR(L, I)=FLUX_DIRECT_CLEAR(L, I)       &
     &              +WEIGHT_INCR*FLUX_DIRECT_INCR_CLEAR(L, I)
                ENDDO
              ENDDO
            ENDIF
            DO I=0, N_LAYER
              DO L=1, N_PROFILE
                FLUX_UP_CLEAR(L, I)=FLUX_UP_CLEAR(L, I)                 &
     &            +WEIGHT_INCR*FLUX_TOTAL_INCR_CLEAR(L, 2*I+1)
                FLUX_DOWN_CLEAR(L, I)=FLUX_DOWN_CLEAR(L, I)             &
     &            +WEIGHT_INCR*FLUX_TOTAL_INCR_CLEAR(L, 2*I+2)
              ENDDO
            ENDDO
          ENDIF
!
        ELSE IF ( (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC).AND. &
     &            (I_SPH_MODE == IP_SPH_MODE_RAD) ) THEN
!

          DO K=1, N_DIRECTION
            DO I=1, N_VIEWING_LEVEL
              DO L=1, N_PROFILE
                RADIANCE(L, I, K)=RADIANCE(L, I, K)                     &
     &            +WEIGHT_INCR*RADIANCE_INCR(L, I, K)
              ENDDO
            ENDDO
          ENDDO
!
          IF (ISOLIR == IP_SOLAR) THEN
            DO I=0, N_LAYER
              DO L=1, N_PROFILE
                I_DIRECT(L, I)=I_DIRECT(L, I)                           &
     &            +WEIGHT_INCR*I_DIRECT_INCR(L, I)
              ENDDO
            ENDDO
          ENDIF
!
        ELSE IF ( (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC).AND. &
     &            (I_SPH_MODE == IP_SPH_MODE_J) ) THEN
!
          DO I=1, N_VIEWING_LEVEL
            DO L=1, N_PROFILE
              PHOTOLYSIS(L, I)=PHOTOLYSIS(L, I)                         &
     &          +WEIGHT_INCR*PHOTOLYSIS_INCR(L, I)
            ENDDO
          ENDDO
!
        ENDIF
!
      ELSE
!
!       Initialization of the radiance field takes place here.
!
        IF ( (I_ANGULAR_INTEGRATION == IP_TWO_STREAM).OR.               &
     &       (I_ANGULAR_INTEGRATION == IP_IR_GAUSS).OR.                 &
     &     ( (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC).AND.      &
     &        (I_SPH_MODE == IP_SPH_MODE_FLUX) ) ) THEN
!
!         Increment the actual fluxes.
          IF (ISOLIR == IP_SOLAR) THEN
            DO I=0, N_LAYER
              DO L=1, N_PROFILE
                FLUX_DIRECT(L, I)=WEIGHT_INCR*FLUX_DIRECT_INCR(L, I)
              ENDDO
            ENDDO
            IF (L_BLUE_FLUX_SURF) THEN
              DO L=1, N_PROFILE
                FLUX_UP_BLUE_SURF(L)                                    &
     &            =WEIGHT_BLUE*FLUX_TOTAL_INCR(L, 2*N_LAYER+1)
                FLUX_DOWN_BLUE_SURF(L)                                  &
     &            =WEIGHT_BLUE*FLUX_TOTAL_INCR(L, 2*N_LAYER+2)
              ENDDO
              IF (ISOLIR == IP_SOLAR) THEN
                DO L=1, N_PROFILE
                  FLUX_DIRECT_BLUE_SURF(L)                              &
     &              =WEIGHT_BLUE*FLUX_DIRECT_INCR(L, N_LAYER)
                ENDDO
              ENDIF
            ENDIF
          ENDIF
          DO I=0, N_LAYER
            DO L=1, N_PROFILE
              FLUX_UP(L, I)=WEIGHT_INCR*FLUX_TOTAL_INCR(L, 2*I+1)
              FLUX_DOWN(L, I)=WEIGHT_INCR*FLUX_TOTAL_INCR(L, 2*I+2)
            ENDDO
          ENDDO
!
          IF (L_CLEAR) THEN
            IF (ISOLIR == IP_SOLAR) THEN
              DO I=0, N_LAYER
                DO L=1, N_PROFILE
                  FLUX_DIRECT_CLEAR(L, I)                               &
     &              =WEIGHT_INCR*FLUX_DIRECT_INCR_CLEAR(L, I)
                ENDDO
              ENDDO
            ENDIF
            DO I=0, N_LAYER
              DO L=1, N_PROFILE
                FLUX_UP_CLEAR(L, I)                                     &
     &            =WEIGHT_INCR*FLUX_TOTAL_INCR_CLEAR(L, 2*I+1)
                FLUX_DOWN_CLEAR(L, I)                                   &
     &            =WEIGHT_INCR*FLUX_TOTAL_INCR_CLEAR(L, 2*I+2)
              ENDDO
            ENDDO
          ENDIF
!
        ELSE IF ( (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC).AND. &
     &            (I_SPH_MODE == IP_SPH_MODE_RAD) ) THEN
!
!         Increment the radiances on levels where they are calculated.
          DO K=1, N_DIRECTION
            DO I=1, N_VIEWING_LEVEL
              DO L=1, N_PROFILE
                RADIANCE(L, I, K)=WEIGHT_INCR*RADIANCE_INCR(L, I, K)
              ENDDO
            ENDDO
          ENDDO
!
          IF (ISOLIR == IP_SOLAR) THEN
            DO I=0, N_LAYER
              DO L=1, N_PROFILE
                I_DIRECT(L, I)=WEIGHT_INCR*I_DIRECT_INCR(L, I)
              ENDDO
            ENDDO
          ENDIF
!
        ELSE IF ( (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC).AND. &
     &            (I_SPH_MODE == IP_SPH_MODE_J) ) THEN
!
          DO I=1, N_VIEWING_LEVEL
            DO L=1, N_PROFILE
              PHOTOLYSIS(L, I)=WEIGHT_INCR*PHOTOLYSIS_INCR(L, I)
            ENDDO
          ENDDO
!
        ENDIF
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE AUGMENT_RADIANCE
#endif
