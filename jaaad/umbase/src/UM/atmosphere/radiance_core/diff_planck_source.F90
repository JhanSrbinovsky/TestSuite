#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate differences in source functions.
!
! Method:
!       Using the polynomial fit to the Planck function, values
!       of this function at the boundaries of layers are found
!       and differences across layers are determined. If the
!       Planckian is being taken to vary quadratically across
!       the layer second differences are found.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE DIFF_PLANCK_SOURCE(N_PROFILE, N_LAYER                  &
     &  , N_DEG_FIT, THERMAL_COEFFICIENT                                &
     &  , T_REF_PLANCK, T_LEVEL, T_GROUND                               &
     &  , PLANCK_FLUX, DIFF_PLANCK, PLANCK_GROUND                       &
     &  , L_IR_SOURCE_QUAD, T, DIFF_PLANCK_2                            &
     &  , I_ANGULAR_INTEGRATION                                         &
     &  , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                  &
     &  , PLANCK_RADIANCE                                               &
     &  , L_TILE, N_POINT_TILE, N_TILE, LIST_TILE                       &
     &  , FRAC_TILE, T_TILE, PLANCK_FLUX_TILE                           &
     &  , ND_PROFILE, ND_LAYER, ND_THERMAL_COEFF                        &
     &  , ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL                         &
     &  , ND_POINT_TILE, ND_TILE                                        &
     &  )
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
!
!     Sizes of dummy arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_THERMAL_COEFF                                              &
!           Size allocated for thermal coefficients
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for profiles where radiances are calculated
     &  , ND_VIEWING_LEVEL                                              &
!           Size allocated for levels where radiances are calculated
     &  , ND_POINT_TILE                                                 &
!           Size allocated for points with surface tiling
     &  , ND_TILE
!           Size allocated for the number of tiles
!
!     Include header files.
#include "angular_integration_pcf3z.h"
#include "c_pi.h"
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_DEG_FIT
!           Degree of fitting function
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_VIEWING_LEVEL                                               &
!           Number of levels where radiances are calculated
     &  , I_RAD_LAYER(ND_VIEWING_LEVEL)
!           Layers in which to intercept radiances
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FRAC_RAD_LAYER(ND_VIEWING_LEVEL)
!           Fractions below the tops of the layers
!
      INTEGER, INTENT(IN) ::                                            &
     &    I_ANGULAR_INTEGRATION
!           Type of angular integration
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Flag for quadratic IR-source
      REAL  (Real64), INTENT(IN) ::                                     &
     &    THERMAL_COEFFICIENT(0: ND_THERMAL_COEFF-1)                    &
!           Coefficients of fit to the Planckian flux function
     &  , T_REF_PLANCK                                                  &
!           Planckian reference temperature
     &  , T_LEVEL(ND_PROFILE, 0: ND_LAYER)                              &
!           Temperatures on levels
     &  , T(ND_PROFILE, ND_LAYER)                                       &
!           Temperatures at centres of layers
     &  , T_GROUND(ND_PROFILE)
!           Temperatures at ground
!
!     Tiling of the surface:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_TILE
!           Local to allow tiling options
      INTEGER, INTENT(IN) ::                                            &
     &    N_POINT_TILE                                                  &
!           Number of points to tile
     &  , N_TILE                                                        &
!           Number of tiles used
     &  , LIST_TILE(ND_POINT_TILE)
!           List of points with surface tiling
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FRAC_TILE(ND_POINT_TILE, ND_TILE)                             &
!           Fraction of tiled grid-points occupied by each tile
     &  , T_TILE(ND_POINT_TILE, ND_TILE)
!           Local surface temperatures on individual tiles
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    PLANCK_FLUX(ND_PROFILE, 0: ND_LAYER)                          &
!           Planckian flux on levels
     &  , DIFF_PLANCK(ND_PROFILE, ND_LAYER)                             &
!           Differences in Planckian flux (bottom-top)
     &  , DIFF_PLANCK_2(ND_PROFILE, ND_LAYER)                           &
!           Twice 2nd differences in the Planckian flux
     &  , PLANCK_GROUND(ND_PROFILE)                                     &
!           Planckian flux at the surface temperature
     &  , PLANCK_RADIANCE(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL)        &
!           Planckian radiances at viewing levels
     &  , PLANCK_FLUX_TILE(ND_POINT_TILE, ND_TILE)
!           Local Planckian fluxes on surface tiles
!
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    T_RATIO(ND_PROFILE)
!           Temperature ratio
!
!
!
      IF (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC) THEN
!
!       Calculate the Planckian radiance on viewing levels.
        DO I=1, N_VIEWING_LEVEL
          DO L=1, N_PROFILE
!           Interpolate linearly in the temperature.
            T_RATIO(L)=(T_LEVEL(L, I_RAD_LAYER(I)-1)                    &
     &        +(T_LEVEL(L, I_RAD_LAYER(I))-T_LEVEL(L, I_RAD_LAYER(I)-1))&
     &        *FRAC_RAD_LAYER(I))/T_REF_PLANCK
!           Use the second differences of the Planckian as temporary
!           storage.
            PLANCK_RADIANCE(L, I)                                       &
     &        =THERMAL_COEFFICIENT(N_DEG_FIT)
          ENDDO
          DO J=N_DEG_FIT-1, 0, -1
            DO L=1, N_PROFILE
              PLANCK_RADIANCE(L, I)                                     &
     &          =PLANCK_RADIANCE(L, I)                                  &
     &          *T_RATIO(L)+THERMAL_COEFFICIENT(J)
            ENDDO
          ENDDO
!
          DO L=1, N_PROFILE
            PLANCK_RADIANCE(L, I)=PLANCK_RADIANCE(L, I)/PI
          ENDDO
!
        ENDDO
!
      ENDIF
!
!
!     Calculate the change in the Planckian flux across each layer.
      DO L=1, N_PROFILE
        T_RATIO(L)=T_LEVEL(L, 0)/T_REF_PLANCK
        PLANCK_FLUX(L, 0)                                               &
     &    =THERMAL_COEFFICIENT(N_DEG_FIT)
      ENDDO
      DO J=N_DEG_FIT-1, 0, -1
        DO L=1, N_PROFILE
          PLANCK_FLUX(L, 0)                                             &
     &      =PLANCK_FLUX(L, 0)                                          &
     &      *T_RATIO(L)+THERMAL_COEFFICIENT(J)
        ENDDO
      ENDDO
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
          T_RATIO(L)=T_LEVEL(L, I)/T_REF_PLANCK
          PLANCK_FLUX(L, I)                                             &
     &      =THERMAL_COEFFICIENT(N_DEG_FIT)
        ENDDO
        DO J=N_DEG_FIT-1, 0, -1
          DO L=1, N_PROFILE
            PLANCK_FLUX(L, I)                                           &
     &        =PLANCK_FLUX(L, I)                                        &
     &        *T_RATIO(L)+THERMAL_COEFFICIENT(J)
          ENDDO
        ENDDO
        DO L=1, N_PROFILE
          DIFF_PLANCK(L, I)=PLANCK_FLUX(L, I)                           &
     &      -PLANCK_FLUX(L, I-1)
        ENDDO
      ENDDO
!
!
!     Calculate the second difference if required.
      IF (L_IR_SOURCE_QUAD) THEN
        DO I=1, N_LAYER
!         Use the second difference for temporary storage.
!         of the Planckian at the middle of the layer.
          DO L=1, N_PROFILE
            T_RATIO(L)=T(L, I)/T_REF_PLANCK
            DIFF_PLANCK_2(L, I)                                         &
     &        =THERMAL_COEFFICIENT(N_DEG_FIT)
          ENDDO
          DO J=N_DEG_FIT-1, 0, -1
            DO L=1, N_PROFILE
              DIFF_PLANCK_2(L, I)                                       &
     &          =DIFF_PLANCK_2(L, I)                                    &
     &          *T_RATIO(L)+THERMAL_COEFFICIENT(J)
            ENDDO
          ENDDO
          DO L=1, N_PROFILE
            DIFF_PLANCK_2(L, I)=2.0E+00_Real64*(PLANCK_FLUX(L, I)       &
     &        +PLANCK_FLUX(L, I-1)-2.0E+00_Real64*DIFF_PLANCK_2(L, I))
          ENDDO
        ENDDO
      ENDIF
!
!
!     Planckian flux at the surface.
      DO L=1, N_PROFILE
        T_RATIO(L)=T_GROUND(L)/T_REF_PLANCK
        PLANCK_GROUND(L)=THERMAL_COEFFICIENT(N_DEG_FIT)
      ENDDO
      DO J=N_DEG_FIT-1, 0, -1
        DO L=1, N_PROFILE
          PLANCK_GROUND(L)=PLANCK_GROUND(L)*T_RATIO(L)                  &
     &      +THERMAL_COEFFICIENT(J)
        ENDDO
      ENDDO
!
!     Local Planckian fluxes will be required on tiled surfaces.
!     Furthermore, the overall Planckian will be calculated as a
!     weighted sum of the individual components: this allows for
!     variations in the Planckian between spectral bands more
!     satisfactorily than the use of an equivalent temperature
!     can.
      IF (L_TILE) THEN
!
        DO K=1, N_TILE
          DO L=1, N_POINT_TILE
            T_RATIO(L)=T_TILE(L, K)/T_REF_PLANCK
            PLANCK_FLUX_TILE(L, K)=THERMAL_COEFFICIENT(N_DEG_FIT)
          ENDDO
          DO J=N_DEG_FIT-1, 0, -1
            DO L=1, N_POINT_TILE
              PLANCK_FLUX_TILE(L, K)=PLANCK_FLUX_TILE(L, K)*T_RATIO(L)  &
     &          +THERMAL_COEFFICIENT(J)
            ENDDO
          ENDDO
        ENDDO
!
        DO L=1, N_POINT_TILE
          PLANCK_GROUND(LIST_TILE(L))                                   &
     &      =FRAC_TILE(L, 1)*PLANCK_FLUX_TILE(L, 1)
        ENDDO
        DO K=2, N_TILE
          DO L=1, N_POINT_TILE
            PLANCK_GROUND(LIST_TILE(L))=PLANCK_GROUND(LIST_TILE(L))     &
     &        +FRAC_TILE(L, K)*PLANCK_FLUX_TILE(L, K)
          ENDDO
        ENDDO
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE DIFF_PLANCK_SOURCE
#endif
