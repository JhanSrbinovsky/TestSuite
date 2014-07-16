#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate fluxes including only gaseous absorption.
!
! Method:
!       Transmission coefficients for each layer are calculated
!       from the gaseous absorption alone. fluxes are propagated
!       upward or downward through the column using these
!       coefficients and source terms.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MONOCHROMATIC_GAS_FLUX(N_PROFILE, N_LAYER              &
     &  , TAU_GAS                                                       &
     &  , ISOLIR, SEC_0, FLUX_INC_DIRECT, FLUX_INC_DOWN                 &
     &  , DIFF_PLANCK, D_PLANCK_FLUX_SURFACE                            &
     &  , DIFFUSE_ALBEDO, DIRECT_ALBEDO                                 &
     &  , DIFFUSIVITY_FACTOR                                            &
     &  , FLUX_DIRECT, FLUX_DIFFUSE                                     &
     &  , ND_PROFILE, ND_LAYER                                          &
     &  )
!
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
!     Include header files
#include "c_kinds.h"
#include "spectral_region_pcf3z.h"
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , ISOLIR
!           Spectral region
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU_GAS(ND_PROFILE, ND_LAYER)                                 &
!           Gaseous optical depths
     &  , SEC_0(ND_PROFILE)                                             &
!           Secant of zenith angle
     &  , FLUX_INC_DIRECT(ND_PROFILE)                                   &
!           Incident direct flux
     &  , FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident diffuse flux
     &  , D_PLANCK_FLUX_SURFACE(ND_PROFILE)                             &
!           Difference in Planckian fluxes between the surface
!           and the overlying air
     &  , DIFF_PLANCK(ND_PROFILE, ND_LAYER)                             &
!           Difference in Planckian function
     &  , DIFFUSE_ALBEDO(ND_PROFILE)                                    &
!           Diffuse surface albedo
     &  , DIRECT_ALBEDO(ND_PROFILE)                                     &
!           Direct surface albedo
     &  , DIFFUSIVITY_FACTOR
!           Diffusivity factor
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_PROFILE, 0: ND_LAYER)                          &
!           Direct flux
     &  , FLUX_DIFFUSE(ND_PROFILE, 2*ND_LAYER+2)
!           Diffuse flux
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    TRANS(ND_PROFILE, ND_LAYER)                                   &
!           Transmissivities
     &  , SOURCE_UP(ND_PROFILE, ND_LAYER)                               &
!           Upward source function
     &  , SOURCE_DOWN(ND_PROFILE, ND_LAYER)
!           Downward source function
!
!     Variables related to the treatment of ill-conditioning
      REAL  (Real64) ::                                                 &
     &    EPS_R                                                         &
!           The smallest real number such that 1.0-EPS_R is not 1
!           to the computer's precision
     &  , SQ_EPS_R
!           The square root of the above
!
!
!
!     Set the tolerances used in avoiding ill-conditioning, testing
!     on any variable.
      EPS_R=EPSILON(TAU_GAS(1, 1))
      SQ_EPS_R=SQRT(EPS_R)
!
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
            TRANS(L, I)=EXP(-DIFFUSIVITY_FACTOR*TAU_GAS(L, I))
        ENDDO
      ENDDO
!
      IF (ISOLIR == IP_SOLAR) THEN
        DO I=1, N_LAYER
          DO L=1, N_PROFILE
            SOURCE_UP(L, I)=0.0E+00_Real64
            SOURCE_DOWN(L, I)=0.0E+00_Real64
          ENDDO
        ENDDO
      ELSE IF (ISOLIR == IP_INFRA_RED) THEN
        DO I=1, N_LAYER
          DO L=1, N_PROFILE
            SOURCE_UP(L, I)=(1.0E+00_Real64-TRANS(L, I)+SQ_EPS_R)       &
     &        *DIFF_PLANCK(L, I)                                        &
     &        /(DIFFUSIVITY_FACTOR*TAU_GAS(L, I)+SQ_EPS_R)
            SOURCE_DOWN(L, I)=-SOURCE_UP(L, I)
          ENDDO
        ENDDO
      ENDIF
!
!     The direct flux.
      IF (ISOLIR == IP_SOLAR) THEN
        DO L=1, N_PROFILE
          FLUX_DIRECT(L, 0)=FLUX_INC_DIRECT(L)
        ENDDO
        DO I=1, N_LAYER
          DO L=1, N_PROFILE
               FLUX_DIRECT(L, I)                                        &
     &            =FLUX_DIRECT(L, I-1)*EXP(-TAU_GAS(L, I)*SEC_0(L))
          ENDDO
        ENDDO
      ENDIF
!
!     Downward fluxes.
      DO L=1, N_PROFILE
        FLUX_DIFFUSE(L, 2)=FLUX_INC_DOWN(L)
      ENDDO
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
          FLUX_DIFFUSE(L, 2*I+2)=TRANS(L, I)*FLUX_DIFFUSE(L, 2*I)       &
     &      +SOURCE_DOWN(L, I)
        ENDDO
      ENDDO
!
!     Upward fluxes.
      IF (ISOLIR == IP_SOLAR) THEN
        DO L=1, N_PROFILE
          FLUX_DIFFUSE(L, 2*N_LAYER+1)=                                 &
     &      +DIFFUSE_ALBEDO(L)*FLUX_DIFFUSE(L, 2*N_LAYER+2)             &
     &      +DIRECT_ALBEDO(L)*FLUX_DIRECT(L, N_LAYER)
        ENDDO
      ELSE
        DO L=1, N_PROFILE
          FLUX_DIFFUSE(L, 2*N_LAYER+1)=D_PLANCK_FLUX_SURFACE(L)         &
     &      +DIFFUSE_ALBEDO(L)*FLUX_DIFFUSE(L, 2*N_LAYER+2)
        ENDDO
      ENDIF
      DO I=N_LAYER, 1, -1
        DO L=1, N_PROFILE
          FLUX_DIFFUSE(L, 2*I-1)=TRANS(L, I)*FLUX_DIFFUSE(L, 2*I+1)     &
     &      +SOURCE_UP(L, I)
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE MONOCHROMATIC_GAS_FLUX
#endif
