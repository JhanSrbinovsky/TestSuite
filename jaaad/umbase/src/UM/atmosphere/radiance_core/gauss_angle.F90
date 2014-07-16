#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate fluxes using gaussian quadrature.
!
! Method:
!       Fluxes are calculated by using gaussian quadrature for
!       the angular integration. This is not a full implementation
!       of gaussian quadrature for multiple scattering, but is
!       intended only for non-scattering calculations in the
!       infra-red. In this case, the fluxes can be calculated as
!       a weighted sum of two-stream fluxes where the diffusivity
!       factors for the two-stream approximations are determined
!       from the gaussian points.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE GAUSS_ANGLE(N_PROFILE, N_LAYER                         &
     &   , N_ORDER_GAUSS                                                &
     &   , TAU                                                          &
     &   , FLUX_INC_DOWN                                                &
     &   , DIFF_PLANCK, SOURCE_GROUND, ALBEDO_SURFACE_DIFF              &
     &   , FLUX_DIFFUSE                                                 &
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                              &
     &   , ND_PROFILE, ND_LAYER                                         &
     &   )
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ND_LAYER
!           Maximum number of layers
!
!     Include header files.
#include "c_kinds.h"
#include "c_pi.h"
#include "angular_integration_pcf3z.h"
#include "spectral_region_pcf3z.h"
#include "gaussian_weight_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_ORDER_GAUSS
!           Order of gaussian integration
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Use quadratic source term
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ND_LAYER)                                     &
!           Optical depth
     &  , ALBEDO_SURFACE_DIFF(ND_PROFILE)                               &
!           Diffuse albedo
     &  , FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident total flux
     &  , DIFF_PLANCK(ND_PROFILE, ND_LAYER)                             &
!           Difference in pi*Planckian function
     &  , SOURCE_GROUND(ND_PROFILE)                                     &
!           Ground source function
     &  , DIFF_PLANCK_2(ND_PROFILE, ND_LAYER)
!             2x2nd differences of Planckian
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIFFUSE(ND_PROFILE, 2*ND_LAYER+2)
!           Diffuse fluxes
!
!     Local variabales.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , I_DUM
!           Dummy variable
      REAL  (Real64) ::                                                 &
     &    FLUX_STREAM(ND_PROFILE, 2*ND_LAYER+2)                         &
!           Flux in stream
     &  , FLUX_NULL(ND_PROFILE, 2*ND_LAYER+2)                           &
!           Array of null fluxes
     &  , SECANT_RAY                                                    &
!           Secant of angle with vertical
     &  , DIFF_PLANCK_RAD(ND_PROFILE, ND_LAYER)                         &
!           Difference in pi*Planckian function
     &  , DIFF_PLANCK_RAD_2(ND_PROFILE, ND_LAYER)                       &
!             2x2nd differences of Planckian
     &  , SOURCE_GROUND_RAD(ND_PROFILE)                                 &
!           Ground source function
     &  , RADIANCE_INC(ND_PROFILE)                                      &
!           Incidnet radiance
     &  , WEIGHT_STREAM                                                 &
!           Weighting for stream
     &  , R_DUM
!           Dummy variable
!
!     Set the gaussian weights for integration.
#include "gaussian_weight_dcf3z.h"
!
!     Subroutines called:
      EXTERNAL                                                          &
     &     MONOCHROMATIC_IR_RADIANCE, AUGMENT_RADIANCE
!
!
!
!     Set the source function.
      DO L=1, N_PROFILE
        SOURCE_GROUND_RAD(L)=SOURCE_GROUND(L)/PI
        RADIANCE_INC(L)=FLUX_INC_DOWN(L)/PI
      ENDDO
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
          DIFF_PLANCK_RAD(L, I)=DIFF_PLANCK(L, I)/PI
        ENDDO
      ENDDO
      DO I=1, 2*N_LAYER+2
        DO L=1, N_PROFILE
          FLUX_DIFFUSE(L, I)=0.0
        ENDDO
      ENDDO
      IF (L_IR_SOURCE_QUAD) THEN
        DO I=1, N_LAYER
          DO L=1, N_PROFILE
            DIFF_PLANCK_RAD_2(L, I)=DIFF_PLANCK_2(L, I)/PI
          ENDDO
        ENDDO
      ENDIF
!
!     Calculate the fluxes with a number of diffusivity factors
!     and sum the results.
      DO K=1, N_ORDER_GAUSS
        SECANT_RAY=2.0E+00_Real64                                       &
     &    /(GAUSS_POINT(K, N_ORDER_GAUSS)+1.0E+00_Real64)
!
!       Calculate the radiance at this angle.
! DEPENDS ON: monochromatic_ir_radiance
        CALL MONOCHROMATIC_IR_RADIANCE(N_PROFILE, N_LAYER               &
     &    , TAU                                                         &
     &    , RADIANCE_INC                                                &
     &    , DIFF_PLANCK_RAD, SOURCE_GROUND_RAD, ALBEDO_SURFACE_DIFF     &
     &    , SECANT_RAY                                                  &
     &    , FLUX_STREAM                                                 &
     &    , ND_PROFILE, ND_LAYER                                        &
     &    )
!
!       Augment the flux by the amount in this stream.
        WEIGHT_STREAM=5.0E-01_Real64*PI*GAUSS_WEIGHT(K, N_ORDER_GAUSS)  &
     &    *(GAUSS_POINT(K, N_ORDER_GAUSS)+1.0E+00_Real64)
! DEPENDS ON: augment_radiance
        CALL AUGMENT_RADIANCE(N_PROFILE, N_LAYER                        &
     &    , IP_IR_GAUSS, I_DUM, I_DUM, I_DUM, I_DUM                     &
     &    , I_DUM, I_DUM, I_DUM                                         &
     &    , IP_INFRA_RED, .FALSE., .FALSE.                              &
     &    , WEIGHT_STREAM                                               &
!                     Actual radiances
     &    , FLUX_NULL, FLUX_DIFFUSE, FLUX_NULL                          &
     &    , R_DUM, R_DUM, R_DUM, R_DUM                                  &
     &    , FLUX_NULL, FLUX_NULL                                        &
!                       Increments to radiances
     &    , FLUX_NULL, FLUX_STREAM, FLUX_NULL                           &
     &    , R_DUM, R_DUM, R_DUM, R_DUM                                  &
     &    , FLUX_NULL, FLUX_NULL                                        &
!                     Dimensions
     &    , ND_PROFILE, 1, ND_LAYER                                     &
     &    , 1, 1, 1, 1                                                  &
     &    )
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE GAUSS_ANGLE
#endif
