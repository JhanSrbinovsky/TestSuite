#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to increment upward fluxes on a tiled surface.
!
! Method:
!       The arrays holding the local cumulative fluxes or radiances
!       on each tile are incremented by the variables suffixed
!       with _INCR, multiplied by appropriate weights. The routine
!       can be called to initialize fluxes.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE AUGMENT_TILED_RADIANCE(IERR                            &
     &  , N_POINT_TILE, N_TILE, LIST_TILE                               &
     &  , I_ANGULAR_INTEGRATION, ISOLIR, L_INITIAL                      &
     &  , WEIGHT_INCR, L_BLUE_FLUX_SURF, WEIGHT_BLUE_INCR               &
!                     Surface characteristics
     &  , RHO_ALB                                                       &
!                     Actual radiances
     &  , FLUX_UP_TILE, FLUX_UP_BLUE_TILE                               &
!                     Increments to radiances
     &  , FLUX_DIRECT_INCR, FLUX_DOWN_INCR                              &
     &  , PLANCK_FLUX_TILE, PLANCK_FLUX_AIR                             &
!                     Dimensions
     &  , ND_FLUX_PROFILE, ND_POINT_TILE, ND_TILE                       &
     &  , ND_BRDF_BASIS_FNC                                             &
     &  )
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_FLUX_PROFILE                                               &
!           Size allocated for points where fluxes are calculated
     &  , ND_POINT_TILE                                                 &
!           Size allocated for points where the surface is tiled
     &  , ND_TILE                                                       &
!           Size allocated for surface tiles
     &  , ND_BRDF_BASIS_FNC
!           Size allocated for BRDF basis functions
!
!     Include header files.
#include "c_kinds.h"
#include "spectral_region_pcf3z.h"
#include "angular_integration_pcf3z.h"
#include "error_pcf3z.h"
#include "def_std_io_icf3z.h"
#include "surface_spec_pcf3z.h"
!
!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_POINT_TILE                                                  &
!           Number of points where the surface is tiled
     &  , N_TILE                                                        &
!           Number of tiles used
     &  , LIST_TILE(ND_POINT_TILE)
!           List of tiled points
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR                                                        &
!           Spectral region
     &  , I_ANGULAR_INTEGRATION
!           Treatment of angular integration
      LOGICAL, INTENT(INOUT) ::                                         &
     &    L_INITIAL
!           Flag to call the routine to initialize the outputs
      LOGICAL, INTENT(IN) ::                                            &
     &    L_BLUE_FLUX_SURF
!           Flag to increment blue surface fluxes
      REAL  (Real64), INTENT(IN) ::                                     &
     &    WEIGHT_INCR                                                   &
!           Weight to apply to increments
     &  , WEIGHT_BLUE_INCR
!           Weight to apply to increments to blue fluxes
!
!                     Surface Characteristics
      REAL  (Real64), INTENT(IN) ::                                     &
     &    RHO_ALB(ND_POINT_TILE, ND_BRDF_BASIS_FNC, ND_TILE)
!           Weighting functions for BRDFs
!
!                     Increments to Fluxes
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FLUX_DIRECT_INCR(ND_FLUX_PROFILE)                             &
!           Increment to mean direct flux
     &  , FLUX_DOWN_INCR(ND_FLUX_PROFILE)
!           Increment to total downward flux
!
!                     Planckian Fluxes
      REAL  (Real64), INTENT(IN) ::                                     &
     &    PLANCK_FLUX_TILE(ND_POINT_TILE, ND_TILE)                      &
!           Local Planckian flux emitted from each tile
     &  , PLANCK_FLUX_AIR(ND_FLUX_PROFILE)
!           Hemispheric Planckian flux at the temperature of the air
!
!                     Total Fluxes
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    FLUX_UP_TILE(ND_POINT_TILE, ND_TILE)                          &
!           Local upward flux on each tile (not weighted by the
!           fractional coverage of the tile)
     &  , FLUX_UP_BLUE_TILE(ND_POINT_TILE, ND_TILE)
!           Local upward blue flux on each tile (not weighted by the
!           fractional coverage of the tile)
!
!
!     Local arguments.
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , LL                                                            &
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
     &       (I_ANGULAR_INTEGRATION == IP_IR_GAUSS) ) THEN
!
!         Increment the actual fluxes.
          IF (ISOLIR == IP_SOLAR) THEN
!
            DO K=1, N_TILE
              DO LL=1, N_POINT_TILE
                L=LIST_TILE(LL)
                FLUX_UP_TILE(LL, K)=FLUX_UP_TILE(LL, K)                 &
     &            +WEIGHT_INCR                                          &
     &            *(RHO_ALB(LL, IP_SURF_ALB_DIFF, K)*FLUX_DOWN_INCR(L)  &
     &            +(RHO_ALB(LL, IP_SURF_ALB_DIR, K)                     &
     &            -RHO_ALB(LL, IP_SURF_ALB_DIFF, K))                    &
     &            *FLUX_DIRECT_INCR(L))
              ENDDO
            ENDDO
!
            IF (L_BLUE_FLUX_SURF) THEN
              DO K=1, N_TILE
                DO LL=1, N_POINT_TILE
                  L=LIST_TILE(LL)
                  FLUX_UP_BLUE_TILE(LL, K)=FLUX_UP_BLUE_TILE(LL, K)     &
     &              +WEIGHT_BLUE_INCR                                   &
     &              *(RHO_ALB(LL, IP_SURF_ALB_DIFF, K)                  &
     &              *FLUX_DOWN_INCR(L)                                  &
     &              +(RHO_ALB(LL, IP_SURF_ALB_DIR, K)                   &
     &              -RHO_ALB(LL, IP_SURF_ALB_DIFF, K))                  &
     &              *FLUX_DIRECT_INCR(L))
                ENDDO
              ENDDO
            ENDIF
!
          ELSE IF (ISOLIR == IP_INFRA_RED) THEN
!
            DO K=1, N_TILE
              DO LL=1, N_POINT_TILE
                L=LIST_TILE(LL)
                FLUX_UP_TILE(LL, K)=FLUX_UP_TILE(LL, K)                 &
     &            +WEIGHT_INCR*(PLANCK_FLUX_TILE(LL, K)                 &
     &            +RHO_ALB(LL, IP_SURF_ALB_DIFF, K)                     &
     &            *(FLUX_DOWN_INCR(L)                                   &
     &            +PLANCK_FLUX_AIR(L)-PLANCK_FLUX_TILE(LL, K)))
              ENDDO
            ENDDO
          ENDIF
!
        ELSE IF (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC) THEN
!
          WRITE(IU_ERR, '(/A)')                                         &
     &      '*** Error: Tiled surfaces have not yet been '              &
     &      , 'implemented with the spherical harmonic solver.'
          IERR=I_ERR_FATAL
          RETURN
!
        ENDIF
!
      ELSE
!
!       Initialization of the radiance field takes place here.
!
        IF ( (I_ANGULAR_INTEGRATION == IP_TWO_STREAM).OR.               &
     &       (I_ANGULAR_INTEGRATION == IP_IR_GAUSS) ) THEN
!
!         Initialize the actual fluxes.
          IF (ISOLIR == IP_SOLAR) THEN
            DO K=1, N_TILE
              DO LL=1, N_POINT_TILE
                L=LIST_TILE(LL)
                FLUX_UP_TILE(LL, K)=WEIGHT_INCR                         &
     &            *(RHO_ALB(LL, IP_SURF_ALB_DIFF, K)*FLUX_DOWN_INCR(L)  &
     &            +(RHO_ALB(LL, IP_SURF_ALB_DIR, K)                     &
     &            -RHO_ALB(LL, IP_SURF_ALB_DIFF, K))                    &
     &            *FLUX_DIRECT_INCR(L))
              ENDDO
            ENDDO
!
            IF (L_BLUE_FLUX_SURF) THEN
              DO K=1, N_TILE
                DO LL=1, N_POINT_TILE
                  L=LIST_TILE(LL)
                  FLUX_UP_BLUE_TILE(LL, K)                              &
     &              =WEIGHT_BLUE_INCR*(RHO_ALB(LL, IP_SURF_ALB_DIFF, K) &
     &              *FLUX_DOWN_INCR(L)                                  &
     &              +(RHO_ALB(LL, IP_SURF_ALB_DIR, K)                   &
     &              -RHO_ALB(LL, IP_SURF_ALB_DIFF, K))                  &
     &              *FLUX_DIRECT_INCR(L))
                ENDDO
              ENDDO
            ENDIF
!
          ELSE IF (ISOLIR == IP_INFRA_RED) THEN
!
            DO K=1, N_TILE
              DO LL=1, N_POINT_TILE
                L=LIST_TILE(LL)
                FLUX_UP_TILE(LL, K)                                     &
     &            =WEIGHT_INCR*(PLANCK_FLUX_TILE(LL, K)                 &
     &            +RHO_ALB(LL, IP_SURF_ALB_DIFF, K)                     &
     &            *(FLUX_DOWN_INCR(L)                                   &
     &            +PLANCK_FLUX_AIR(L)-PLANCK_FLUX_TILE(LL, K)))
              ENDDO
            ENDDO
!
          ENDIF
!
        ELSE IF (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC) THEN
!
          WRITE(IU_ERR, '(/A)')                                         &
     &      '*** Error: Tiled surfaces have not yet been '              &
     &      , 'implemented with the spherical harmonic solver.'
          IERR=I_ERR_FATAL
          RETURN
!
        ENDIF
!
!       Now reset the initialization flag as the arrays have been set.
        L_INITIAL=.FALSE.
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE AUGMENT_TILED_RADIANCE
#endif
