#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate fluxes in a homogeneous column directly.
!
! Method:
!       Straightforward.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVER_HOMOGEN_DIRECT(N_PROFILE, N_LAYER               &
     &  , TRANS, REFLECT                                                &
     &  , S_DOWN, S_UP                                                  &
     &  , ISOLIR, DIFFUSE_ALBEDO, DIRECT_ALBEDO                         &
     &  , FLUX_DIRECT_GROUND, FLUX_INC_DOWN                             &
     &  , D_PLANCK_FLUX_SURFACE                                         &
     &  , FLUX_TOTAL                                                    &
     &  , ND_PROFILE, ND_LAYER                                          &
     &  )
!
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
#include "spectral_region_pcf3z.h"
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER
!           Size allocated for atmospheric layers
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
     &    TRANS(ND_PROFILE, ND_LAYER)                                   &
!           Transmission coefficient
     &  , REFLECT(ND_PROFILE, ND_LAYER)                                 &
!           Reflection coefficient
     &  , S_DOWN(ND_PROFILE, ND_LAYER)                                  &
!           Downward diffuse source
     &  , S_UP(ND_PROFILE, ND_LAYER)                                    &
!           Upward diffuse source
     &  , DIFFUSE_ALBEDO(ND_PROFILE)                                    &
!           Diffuse surface albedo
     &  , DIRECT_ALBEDO(ND_PROFILE)                                     &
!           Direct surface albedo
     &  , D_PLANCK_FLUX_SURFACE(ND_PROFILE)                             &
!           Difference between the Planckian flux at the surface
!           temperature and that of the overlaying air
     &  , FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident total flux
     &  , FLUX_DIRECT_GROUND(ND_PROFILE)
!           Direct flux at ground level
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_TOTAL(ND_PROFILE, 2*ND_LAYER+2)
!           Total flux
!
!     Declaration of local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
      REAL  (Real64) ::                                                 &
     &    ALPHA(ND_PROFILE, ND_LAYER+1)                                 &
!           Combined albedo of lower layers
     &  , BETA(ND_PROFILE, ND_LAYER)                                    &
!           Working array
     &  , GAMMA(ND_PROFILE, ND_LAYER)                                   &
!           Working array
     &  , H(ND_PROFILE, ND_LAYER)                                       &
!           Working array
     &  , S_UP_PRIME(ND_PROFILE, ND_LAYER+1)
!           Modified upward source function
!
!
!
!     Initialization at the bottom for upward elimination:
      IF (ISOLIR == IP_SOLAR) THEN
        DO L=1, N_PROFILE
          ALPHA(L, N_LAYER+1)=DIFFUSE_ALBEDO(L)
          S_UP_PRIME(L, N_LAYER+1)                                      &
     &      =(DIRECT_ALBEDO(L)-DIFFUSE_ALBEDO(L))                       &
     &      *FLUX_DIRECT_GROUND(L)
        ENDDO
      ELSE IF (ISOLIR == IP_INFRA_RED) THEN
        DO L=1, N_PROFILE
          ALPHA(L, N_LAYER+1)=DIFFUSE_ALBEDO(L)
          S_UP_PRIME(L, N_LAYER+1)                                      &
     &      =(1.0E+00_Real64-DIFFUSE_ALBEDO(L))                         &
     &      *D_PLANCK_FLUX_SURFACE(L)
        ENDDO
      ENDIF
!
!     Eliminating loop:
      DO I=N_LAYER, 1, -1
        DO L=1, N_PROFILE
          BETA(L, I)=1.0E+00_Real64                                     &
     &      /(1.0E+00_Real64-ALPHA(L, I+1)*REFLECT(L, I))
          GAMMA(L, I)=ALPHA(L, I+1)*TRANS(L, I)
          H(L, I)=S_UP_PRIME(L, I+1)+ALPHA(L, I+1)*S_DOWN(L, I)
          ALPHA(L, I)=REFLECT(L, I)                                     &
     &      +BETA(L, I)*GAMMA(L, I)*TRANS(L, I)
          S_UP_PRIME(L, I)=S_UP(L, I)+BETA(L, I)*TRANS(L, I)*H(L, I)
        ENDDO
      ENDDO
!
!     Initialize for backward substitution.
      DO L=1, N_PROFILE
        FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
        FLUX_TOTAL(L, 1)=ALPHA(L, 1)*FLUX_TOTAL(L, 2)+S_UP_PRIME(L, 1)
      ENDDO
!
!     Backward substitution:
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
!         Upward flux
          FLUX_TOTAL(L, 2*I+1)                                          &
     &      =BETA(L, I)*(H(L, I)+GAMMA(L, I)*FLUX_TOTAL(L, 2*I))
!         Downward flux
          FLUX_TOTAL(L, 2*I+2)=S_DOWN(L, I)                             &
     &      +TRANS(L, I)*FLUX_TOTAL(L, 2*I)                             &
     &      +REFLECT(L, I)*FLUX_TOTAL(L, 2*I+1)
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SOLVER_HOMOGEN_DIRECT
#endif
