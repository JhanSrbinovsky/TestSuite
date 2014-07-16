#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate fluxes using Gaussian quadrature.
!
! Method:
!       Fluxes are calculated by using Gaussian quadrature for
!       the angular integration. This is not a full implementation
!       of Gaussian quadrature for multiple scattering, but is
!       intended only for non-scattering calculations in the
!       infra-red. In this case, the fluxes can be calculated as
!       a weighted sum of two-stream fluxes where the diffusivity
!       factors for the two-stream approximations are determined
!       from the Gaussian points.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE GAUSS_ANGLE(N_PROFILE, N_LAYER, L_NET, N_AUGMENT       &
     &   , N_ORDER_GAUSS                                                &
     &   , TAU                                                          &
     &   , FLUX_INC_DOWN                                                &
     &   , DIFF_PLANCK, SOURCE_GROUND, ALBEDO_SURFACE_DIFF              &
     &   , FLUX_DIFFUSE                                                 &
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                              &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
#include "spcrg3a.h"
#include "gsswtp3a.h"
#include "c_pi.h"
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_AUGMENT                                                    &
!             SIZE OF ARRAY TO INCREMENT
     &   , N_ORDER_GAUSS
!             ORDER OF GAUSSIAN INTEGRATION
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_NET                                                        &
!             NET FLUXES REQUIRED
     &   , L_IR_SOURCE_QUAD
!             USE QUADRATIC SOURCE TERM
      REAL                                                              &
                !, INTENT(IN)
     &     TAU(NPD_PROFILE, NPD_LAYER)                                  &
!             OPTICAL DEPTH
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)                             &
!             DIFFUSE ALBEDO
     &   , FLUX_INC_DOWN(NPD_PROFILE)                                   &
!             INCIDENT TOTAL FLUX
     &   , DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)                          &
!             DIFFERENCE IN PI*PLANCK FUNCTION
     &   , SOURCE_GROUND(NPD_PROFILE)                                   &
!             GROUND SOURCE FUNCTION
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)
!             2x2ND DIFFERENCES OF PLANCKIAN
      REAL                                                              &
                !, INTENT(OUT)
     &     FLUX_DIFFUSE(NPD_PROFILE, 2*NPD_LAYER+2)
!             DIFFUSE FLUXES
!
!     LOCAL VARIABALES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , K
!             LOOP VARIABLE
      REAL                                                              &
     &     FLUX_STREAM(NPD_PROFILE, 2*NPD_LAYER+2)                      &
!             FLUX IN STREAM
     &   , FLUX_NULL(NPD_PROFILE, 2*NPD_LAYER+2)                        &
!             ARRAY OF NULL FLUXES
     &   , SECANT_RAY                                                   &
!             SECANT OF ANGLE WITH VERTICAL
     &   , DIFF_PLANCK_RAD(NPD_PROFILE, NPD_LAYER)                      &
!             DIFFERENCE IN PI*PLANCK FUNCTION
     &   , DIFF_PLANCK_RAD_2(NPD_PROFILE, NPD_LAYER)                    &
!             2x2ND DIFFERENCES OF PLANCKIAN
     &   , SOURCE_GROUND_RAD(NPD_PROFILE)                               &
!             GROUND SOURCE FUNCTION
     &   , RADIANCE_INC(NPD_PROFILE)                                    &
!             INCIDNET RADIANCE
     &   , WEIGHT_STREAM
!             WEIGHTING FOR STREAM
!
!     SET THE GAUSSIAN WEIGHTS FOR INTEGRATION.
#include "gsswtd3a.h"
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     MONOCHROMATIC_IR_RADIANCE, AUGMENT_FLUX
!
!
!
!     SET THE SOURCE FUNCTION.
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
!     CALCULATE THE FLUXES WITH A NUMBER OF DIFFUSIVITY FACTORS
!     AND SUM THE RESULTS.
      DO K=1, N_ORDER_GAUSS
         SECANT_RAY=2.0E+00/(GAUSS_POINT(K, N_ORDER_GAUSS)+1.0E+00)
!
!        CALCULATE THE RADIANCE AT THIS ANGLE.
! DEPENDS ON: monochromatic_ir_radiance
            CALL MONOCHROMATIC_IR_RADIANCE(N_PROFILE, N_LAYER           &
     &         , L_NET                                                  &
     &         , TAU                                                    &
     &         , RADIANCE_INC                                           &
     &         , DIFF_PLANCK_RAD, SOURCE_GROUND_RAD, ALBEDO_SURFACE_DIFF&
     &         , SECANT_RAY                                             &
     &         , FLUX_STREAM                                            &
     &         , NPD_PROFILE, NPD_LAYER                                 &
     &         )
!
!        AUGMENT THE FLUX BY THE AMOUNT IN THIS STREAM.
         WEIGHT_STREAM=5.0E-01*PI*GAUSS_WEIGHT(K, N_ORDER_GAUSS)        &
     &      *(GAUSS_POINT(K, N_ORDER_GAUSS)+1.0E+00)
! DEPENDS ON: augment_flux
         CALL AUGMENT_FLUX(N_PROFILE, N_LAYER, N_AUGMENT                &
     &      , IP_INFRA_RED, .FALSE.                                     &
     &      , WEIGHT_STREAM                                             &
     &      , FLUX_NULL, FLUX_DIFFUSE                                   &
     &      , FLUX_NULL, FLUX_STREAM                                    &
     &      , FLUX_NULL, FLUX_NULL                                      &
     &      , FLUX_NULL, FLUX_NULL                                      &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
!
      ENDDO
!
!
      RETURN
      END SUBROUTINE GAUSS_ANGLE
#endif
#endif
