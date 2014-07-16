#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to increment radiances at a given azimuthal order.
!
! Method:
!   The weights of the terms in the complementary function
!   of the direct solution by spherical harmonics, u_{imk}^+-
!   are now available. For each viewing level and direction
!   we multiply by the precalculated coefficients and the
!   factor representing the azimuthal dependence to complete the
!   calculation of the radiance.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE INCREMENT_RAD_CF(N_PROFILE                             &
     &  , N_DIRECTION, AZIM_FACTOR                                      &
     &  , N_VIEWING_LEVEL, I_RAD_LAYER                                  &
     &  , I_SPH_MODE, I_SPH_ALGORITHM, MS, LS_TRUNC, EULER_FACTOR       &
     &  , ISOLIR, MU_0, KAPPA, UP_LM                                    &
     &  , N_RED_EIGENSYSTEM, N_EQUATION, WEIGHT_U, UPM                  &
     &  , I_DIRECT, C_YLM, FLUX_DIRECT, FLUX_TOTAL                      &
     &  , RADIANCE, J_RADIANCE                                          &
     &  , ND_PROFILE, ND_FLUX_PROFILE                                   &
     &  , ND_RADIANCE_PROFILE, ND_J_PROFILE                             &
     &  , ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL                      &
     &  , ND_MAX_ORDER, ND_SPH_EQUATION, ND_SPH_CF_WEIGHT               &
     &  , ND_SPH_U_RANGE                                                &
     &  )
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for points where radiances are calculated
     &  , ND_FLUX_PROFILE                                               &
!           Size allocated for profiles where fluxes are calculated
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for profiles where radiances are calculated
     &  , ND_J_PROFILE                                                  &
!           Size allocated for profiles where mean radiances
!           are calculated
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_DIRECTION                                                  &
!           Size allocated for order of spherical calculation
     &  , ND_VIEWING_LEVEL                                              &
!           Size allocated for levels where radiances are calculated
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of direct calculation of
!           spherical harmonics
     &  , ND_SPH_EQUATION                                               &
!           Size allocated for spherical equations
     &  , ND_SPH_CF_WEIGHT                                              &
!           Size allocated for entities to be weighted by the C. F.
     &  , ND_SPH_U_RANGE
!           Size allowed for range of values of u^+|- contributing
!           on any viewing level
!
!     Include header files:
#include "c_kinds.h"
#include "c_pi.h"
#include "spectral_region_pcf3z.h"
#include "sph_mode_pcf3z.h"
#include "sph_algorithm_pcf3z.h"
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE
!           Number of profiles
!
!     Spectral decomposition
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Spectral region
!
!     Viewing geometry:
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION                                                   &
!           Number of directions in which radiances are calculated
     &  , N_VIEWING_LEVEL                                               &
!           Number of levels where the radiance is calculated
     &  , I_RAD_LAYER(ND_VIEWING_LEVEL)
!           Layers of the atmosphere in which viewing levels fall
      REAL  (Real64), INTENT(IN) ::                                     &
     &    AZIM_FACTOR(ND_PROFILE, ND_DIRECTION)                         &
!           Factors for representing the azimuthal dependence
     &  , MU_0(ND_PROFILE)                                              &
!           Cosines of the solar zenith angles
     &  , KAPPA(ND_MAX_ORDER/2, ND_MAX_ORDER/2)                         &
!           Integrals of Y_l^m*.Y_l^m over the downward hemisphere
     &  , UP_LM(ND_PROFILE, ND_MAX_ORDER+1, ND_DIRECTION)
!           Polar parts of spherical harmonics
!
!     Angular integration:
      INTEGER, INTENT(IN) ::                                            &
     &    I_SPH_MODE                                                    &
!           Mode in which the spherical harmonic code is called
     &  , I_SPH_ALGORITHM                                               &
!           Algorithm used to solve spherical harmonic problem
     &  , MS                                                            &
!           Azimuthal order of truncation
     &  , LS_TRUNC
!           Polar order of truncation
      REAL  (Real64), INTENT(IN) ::                                     &
     &    EULER_FACTOR
!           Factor weighting the last term of the series
!
!     Components of the solution of the linear system
      INTEGER, INTENT(IN) ::                                            &
     &    N_RED_EIGENSYSTEM                                             &
!           Size of the reduced eigensystem
     &  , N_EQUATION
!           Number of spherical equations
      REAL  (Real64), INTENT(IN) ::                                     &
     &    WEIGHT_U(ND_PROFILE, ND_VIEWING_LEVEL                         &
     &      , ND_SPH_CF_WEIGHT, ND_SPH_U_RANGE)                         &
!           Weights for coefficients in equations
     &  , UPM(ND_PROFILE, ND_SPH_EQUATION)
!           Variables u+|-
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    I_DIRECT(ND_PROFILE, 0: ND_LAYER)
!           Direct radiances
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    C_YLM(ND_PROFILE, ND_VIEWING_LEVEL, LS_TRUNC+1-MS)
!           Spherical harmonic coefficients for radiances
!
!     Calculated radiances
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    RADIANCE(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL, ND_DIRECTION)
!           Radiances
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_FLUX_PROFILE, 0: ND_LAYER)                     &
!           Direct fluxes
     &  , FLUX_TOTAL(ND_FLUX_PROFILE, 2*ND_LAYER+2)
!           Total fluxes
!
!           Mean radiances
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    J_RADIANCE(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Mean radiances
!
!
!     Local arguments.
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , ID                                                            &
!           Loop variable
     &  , IV                                                            &
!           Loop variable (viewing level)
     &  , IE                                                            &
!           Loop variable
     &  , LS                                                            &
!           Polar order
     &  , LSR                                                           &
!           Reduced polar order
     &  , OFFSET_U
!           Offset applied to the elements of u^+|- to move to
!           elements relevant to the current layer
      REAL  (Real64) ::                                                 &
     &    CONTRIBUTION                                                  &
!           Contribution of the current order to the flux
     &  , CNST_LS
!           Constant term involving the polar order
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    EVAL_UPLM
!
!
!
      IF (I_SPH_ALGORITHM == IP_SPH_DIRECT) THEN
!
!       Radiances or fluxes are calculated directly from the
!       spherical harmonics.
!
!       Determine the coefficients of the spherical harmonics
!       from the solution of the eigenproblem.
        DO IV=1, N_VIEWING_LEVEL
          OFFSET_U=2*N_RED_EIGENSYSTEM*(I_RAD_LAYER(IV)-1)
          DO K=1, 2*N_RED_EIGENSYSTEM
            DO LSR=1, LS_TRUNC+1-MS
              DO L=1, N_PROFILE
                C_YLM(L, IV, LSR)=C_YLM(L, IV, LSR)                     &
     &            +WEIGHT_U(L, IV, LSR, K)*UPM(L, K+OFFSET_U)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
!
        IF (I_SPH_MODE == IP_SPH_MODE_FLUX) THEN
!
!         Although this routine is called to increment radiances over
!         angular orders, when run to calculate fluxes it should
!         only be called once during each monochromatic calculation.
!
          DO IV=1, N_VIEWING_LEVEL
            DO L=1, N_PROFILE
              CONTRIBUTION=C_YLM(L, IV, 2)*SQRT(PI/3.0E+00_Real64)
!             Upward flux:
              FLUX_TOTAL(L, 2*IV-1)=CONTRIBUTION
!             Downward flux:
              FLUX_TOTAL(L, 2*IV)=-CONTRIBUTION
            ENDDO
          ENDDO
          DO LS=0, LS_TRUNC, 2
            CNST_LS=2.0E+00_Real64*KAPPA(1, (LS+2)/2)                   &
     &        *SQRT(PI/3.0E+00_Real64)
            DO IV=1, N_VIEWING_LEVEL
              DO L=1, N_PROFILE
                CONTRIBUTION=CNST_LS*C_YLM(L, IV, LS+1)
                FLUX_TOTAL(L, 2*IV-1)                                   &
     &            =FLUX_TOTAL(L, 2*IV-1)-CONTRIBUTION
                FLUX_TOTAL(L, 2*IV)                                     &
     &            =FLUX_TOTAL(L, 2*IV)-CONTRIBUTION
              ENDDO
            ENDDO
          ENDDO
!
          IF (ISOLIR == IP_SOLAR) THEN
            DO IV=1, N_VIEWING_LEVEL
              DO L=1, N_PROFILE
                FLUX_DIRECT(L, IV-1)=I_DIRECT(L, IV-1)*MU_0(L)
                FLUX_TOTAL(L, 2*IV)=FLUX_TOTAL(L, 2*IV)                 &
     &            +FLUX_DIRECT(L, IV-1)
              ENDDO
            ENDDO
          ENDIF
!
        ELSE IF (I_SPH_MODE == IP_SPH_MODE_J) THEN
!
!         Although this routine is called to increment radiances over
!         angular orders, when run to calculate mean radiances it should
!         be called only once during each monochromatic calculation.
!
          DO IV=1, N_VIEWING_LEVEL
            DO L=1, N_PROFILE
              J_RADIANCE(L, IV)=C_YLM(L, IV, 2)*SQRT(4.0E+00_Real64*PI)
            ENDDO
          ENDDO
!
          IF (ISOLIR == IP_SOLAR) THEN
            DO IV=1, N_VIEWING_LEVEL
              DO L=1, N_PROFILE
                J_RADIANCE(L, IV)=J_RADIANCE(L, IV)                     &
     &            +I_DIRECT(L, IV)
              ENDDO
            ENDDO
          ENDIF
!
        ELSE IF (I_SPH_MODE == IP_SPH_MODE_RAD) THEN
!
!         Determine the radiances directly from the amplitudes of
!         the harmonics.
          DO ID=1, N_DIRECTION
!
!           Add in the contributions on each viewing level. To improve
!           convergence of the alternating series the contribution
!           from the last term may be reduced in size.
            DO IV=1, N_VIEWING_LEVEL
              DO LSR=1, LS_TRUNC-MS
                DO L=1, N_PROFILE
                  RADIANCE(L, IV, ID)=RADIANCE(L, IV, ID)               &
     &              +AZIM_FACTOR(L, ID)*C_YLM(L, IV, LSR)               &
     &              *UP_LM(L, LSR, ID)
                ENDDO
              ENDDO
              DO L=1, N_PROFILE
                RADIANCE(L, IV, ID)=RADIANCE(L, IV, ID)+EULER_FACTOR    &
     &            *AZIM_FACTOR(L, ID)*C_YLM(L, IV, LS_TRUNC+1-MS)       &
     &            *UP_LM(L, LS_TRUNC+1-MS, ID)
              ENDDO
            ENDDO
!
          ENDDO
        ENDIF
!
      ELSE IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
!
        DO ID=1, N_DIRECTION
          DO IV=1, N_VIEWING_LEVEL
            DO IE=1, N_EQUATION
              DO L=1, N_PROFILE
                RADIANCE(L, IV, ID)=RADIANCE(L, IV, ID)                 &
     &            +AZIM_FACTOR(L, ID)                                   &
     &            *WEIGHT_U(L, IV, ID, IE)*UPM(L, IE)
              ENDDO
            ENDDO
          ENDDO
        ENDDO

!
      ENDIF
!
!
!
1000   format(4e14.3)
      RETURN
      END SUBROUTINE INCREMENT_RAD_CF
#endif
