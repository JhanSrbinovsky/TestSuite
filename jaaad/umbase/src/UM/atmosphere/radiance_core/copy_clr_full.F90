#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set clear-sky optical properties.
!
! Method:
!       The arrays of clear-sky optical properties at the top
!       of the column and of total optical properties lower
!       down are combined to give a sinle array of clear-sky
!       optical properties.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE COPY_CLR_FULL(N_PROFILE, N_LAYER, N_CLOUD_TOP          &
     &  , N_ORDER_PHASE                                                 &
     &  , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR                             &
     &  , TAU, OMEGA, PHASE_FNC                                         &
     &  , TAU_CLR_F, OMEGA_CLR_F, PHASE_FNC_CLR_F                       &
!                       Sizes of arrays
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_MAX_ORDER       &
     &  )
!
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
!
!     Sizes of arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for totally clear layers
     &  , ID_CT                                                         &
!           Topmost declared layer for cloudy optical properties
     &  , ND_MAX_ORDER
!           Size allowed for orders of spherical harmonics
!
!                     Atmospheric properties
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_ORDER_PHASE
!           Number of terms in the phase function
!
!                     Optical properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU_CLR(ND_PROFILE, ND_LAYER_CLR)                             &
!           Optical depth in totally clear region
     &  , OMEGA_CLR(ND_PROFILE, ND_LAYER_CLR)                           &
!           Single scattering albedo in totally clear region
     &  , PHASE_FNC_CLR(ND_PROFILE, ND_LAYER_CLR, ND_MAX_ORDER)
!           Phase function in totally clear region
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ID_CT: ND_LAYER)                              &
!           Optical depth restricted to clear-sky regions
     &  , OMEGA(ND_PROFILE, ID_CT: ND_LAYER)                            &
!           ALbedo of single scattering restricted to clear-sky regions
     &  , PHASE_FNC(ND_PROFILE, ID_CT: ND_LAYER, ND_MAX_ORDER)
!           Phase function restricted to clear-sky regions
!
!                     Single scattering properties
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    TAU_CLR_F(ND_PROFILE, ND_LAYER)                               &
!           Optical depth
     &  , OMEGA_CLR_F(ND_PROFILE, ND_LAYER)                             &
!           Single scattering albedo
     &  , PHASE_FNC_CLR_F(ND_PROFILE, ND_LAYER, ND_MAX_ORDER)
!           Phase function
!
!
!
!     Local variables.
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , I                                                             &
!           Loop variable
     &  , K
!           Loop variable
!
!
!
!     Above cloud top.
      DO I=1, N_CLOUD_TOP-1
        DO L=1, N_PROFILE
          TAU_CLR_F(L, I)=TAU_CLR(L, I)
          OMEGA_CLR_F(L, I)=OMEGA_CLR(L, I)
          PHASE_FNC_CLR_F(L, I, 1)=PHASE_FNC_CLR(L, I, 1)
        ENDDO
        DO K=2, N_ORDER_PHASE
          DO L=1, N_PROFILE
            PHASE_FNC_CLR_F(L, I, K)=PHASE_FNC_CLR(L, I, K)
          ENDDO
        ENDDO
      ENDDO
!
!     Below cloud top.
      DO I=N_CLOUD_TOP, N_LAYER
        DO L=1, N_PROFILE
          TAU_CLR_F(L, I)=TAU(L, I)
          OMEGA_CLR_F(L, I)=OMEGA(L, I)
          PHASE_FNC_CLR_F(L, I, 1)=PHASE_FNC(L, I, 1)
        ENDDO
        DO K=2, N_ORDER_PHASE
          DO L=1, N_PROFILE
            PHASE_FNC_CLR_F(L, I, K)=PHASE_FNC(L, I, K)
          ENDDO
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE COPY_CLR_FULL
#endif
