#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve for fluxes treating scattering approximately.
!
! Method:
!       The routine is applicable in the infra-red. downward
!       differential fluxes are calculated first assuming that the
!       upward differential fluxes are 0. Upward fluxes are then
!       calculated using the previously calculated downward fluxes
!       in the reflected terms.
!
! Current owner of code: James Manners
!
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MIX_APP_SCAT(N_PROFILE, N_LAYER, N_CLOUD_TOP           &
     &   , T_FREE, R_FREE, S_DOWN_FREE, S_UP_FREE                       &
     &   , T_CLOUD, R_CLOUD, S_DOWN_CLOUD, S_UP_CLOUD                   &
     &   , G_FF, G_FC, G_CF, G_CC                                       &
     &   , B_FF, B_FC, B_CF, B_CC                                       &
     &   , FLUX_INC_DOWN                                                &
     &   , SOURCE_GROUND, ALBEDO_SURFACE_DIFF                           &
     &   , FLUX_DIFFUSE                                                 &
     &   , ND_PROFILE, ND_LAYER, ID_CT                                  &
     &   )
!
!
!
      IMPLICIT NONE
!
! Include Header Files
#include "c_kinds.h"
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ID_CT
!           Topmost declared cloudy layer
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP
!           Topmost cloudy layer
      REAL  (Real64), INTENT(IN) ::                                     &
     &    T_FREE(ND_PROFILE, ND_LAYER)                                  &
!           Free transmission
     &  , R_FREE(ND_PROFILE, ND_LAYER)                                  &
!           Free reflection
     &  , S_DOWN_FREE(ND_PROFILE, ND_LAYER)                             &
!           Free downward source function
     &  , S_UP_FREE(ND_PROFILE, ND_LAYER)                               &
!           Free upward source function
     &  , T_CLOUD(ND_PROFILE, ND_LAYER)                                 &
!           Cloudy transmission
     &  , R_CLOUD(ND_PROFILE, ND_LAYER)                                 &
!           Cloudy reflection
     &  , S_DOWN_CLOUD(ND_PROFILE, ND_LAYER)                            &
!           Downward cloudy source function
     &  , S_UP_CLOUD(ND_PROFILE, ND_LAYER)
!           Upward cloudy source function
      REAL  (Real64), INTENT(IN) ::                                     &
     &    G_FF(ND_PROFILE, ID_CT-1: ND_LAYER)                           &
!           Energy transfer coefficient
     &  , G_FC(ND_PROFILE, ID_CT-1: ND_LAYER)                           &
!           Energy transfer coefficient
     &  , G_CF(ND_PROFILE, ID_CT-1: ND_LAYER)                           &
!           Energy transfer coefficient
     &  , G_CC(ND_PROFILE, ID_CT-1: ND_LAYER)                           &
!           Energy transfer coefficient
     &  , B_FF(ND_PROFILE, ID_CT-1: ND_LAYER)                           &
!           Energy transfer coefficient
     &  , B_FC(ND_PROFILE, ID_CT-1: ND_LAYER)                           &
!           Energy transfer coefficient
     &  , B_CF(ND_PROFILE, ID_CT-1: ND_LAYER)                           &
!           Energy transfer coefficient
     &  , B_CC(ND_PROFILE, ID_CT-1: ND_LAYER)
!           Energy transfer coefficient
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident diffuse flux
     &  , SOURCE_GROUND(ND_PROFILE)                                     &
!           Source from ground
     &  , ALBEDO_SURFACE_DIFF(ND_PROFILE)
!           Diffuse albedo
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIFFUSE(ND_PROFILE, 2*ND_LAYER+2)
!           Diffuse flux
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
      REAL  (Real64) ::                                                 &
     &    FLUX_DOWN(ND_PROFILE, 0: ND_LAYER)                            &
!           Downward fluxes outside clouds just below i'th level
     &  , FLUX_DOWN_CLOUD(ND_PROFILE, 0: ND_LAYER)                      &
!           Downward fluxes inside clouds just below i'th level
     &  , FLUX_UP(ND_PROFILE, 0: ND_LAYER)                              &
!           Upward fluxes outside clouds just above i'th level
     &  , FLUX_UP_CLOUD(ND_PROFILE, 0: ND_LAYER)                        &
!           Upward fluxes inside clouds just above i'th level
     &  , FLUX_PROPAGATED                                               &
!           Temporary propagated flux outside cloud
     &  , FLUX_PROPAGATED_CLOUD                                         &
!           Temporary propagated flux inside cloud
     &  , FLUX_CLOUD_TOP(ND_PROFILE)
!           Total downward flux at top of cloud
!
!
!
!     The arrays flux_down and flux_up will eventually contain the total
!     fluxes, but initially they are used for the clear fluxes.
!     Note that downward fluxes refer to values just below the interface
!     and upward fluxes to values just above it.
!
!
!     Downward flux:
!
!     Region above clouds:
      DO L=1, N_PROFILE
        FLUX_DOWN(L, 0)=FLUX_INC_DOWN(L)
      ENDDO
      DO I=1, N_CLOUD_TOP-1
        DO L=1, N_PROFILE
          FLUX_DOWN(L, I)=T_FREE(L, I)*FLUX_DOWN(L, I-1)                &
     &      +S_DOWN_FREE(L, I)
        ENDDO
      ENDDO
      DO L=1, N_PROFILE
        FLUX_CLOUD_TOP(L)=FLUX_DOWN(L, N_CLOUD_TOP-1)
      ENDDO
!
!     Region of clouds:
      DO L=1, N_PROFILE
        FLUX_DOWN(L, N_CLOUD_TOP-1)                                     &
     &    =G_FF(L, N_CLOUD_TOP-1)*FLUX_CLOUD_TOP(L)
        FLUX_DOWN_CLOUD(L, N_CLOUD_TOP-1)                               &
     &    =G_FC(L, N_CLOUD_TOP-1)*FLUX_CLOUD_TOP(L)
      ENDDO
!
      DO I=N_CLOUD_TOP, N_LAYER-1
        DO L=1, N_PROFILE
!
!         Propagate downward fluxes through the layer.
          FLUX_PROPAGATED=T_FREE(L, I)*FLUX_DOWN(L, I-1)                &
     &      +S_DOWN_FREE(L, I)
          FLUX_PROPAGATED_CLOUD=T_CLOUD(L, I)*FLUX_DOWN_CLOUD(L, I-1)   &
     &      +S_DOWN_CLOUD(L, I)
!         Transfer downward fluxes across the interface.
          FLUX_DOWN(L, I)                                               &
     &      =G_FF(L, I)*FLUX_PROPAGATED                                 &
     &      +G_CF(L, I)*FLUX_PROPAGATED_CLOUD
          FLUX_DOWN_CLOUD(L, I)                                         &
     &      =G_CC(L, I)*FLUX_PROPAGATED_CLOUD                           &
     &      +G_FC(L, I)*FLUX_PROPAGATED
!
        ENDDO
      ENDDO
!
!     Propagate across the bottom layer, but without transferring
!     across the surface and form the reflected beams.
      DO L=1, N_PROFILE
!       Propagate downward fluxes through the layer.
        FLUX_DOWN(L, N_LAYER)                                           &
     &    =T_FREE(L, N_LAYER)*FLUX_DOWN(L, N_LAYER-1)                   &
     &    +S_DOWN_FREE(L, N_LAYER)
        FLUX_DOWN_CLOUD(L, N_LAYER)                                     &
     &    =T_CLOUD(L, N_LAYER)*FLUX_DOWN_CLOUD(L, N_LAYER-1)            &
     &    +S_DOWN_CLOUD(L, N_LAYER)
        FLUX_UP(L, N_LAYER)                                             &
     &    =ALBEDO_SURFACE_DIFF(L)*FLUX_DOWN(L, N_LAYER)                 &
     &    +B_FF(L, N_LAYER)*SOURCE_GROUND(L)
        FLUX_UP_CLOUD(L, N_LAYER)                                       &
     &    =ALBEDO_SURFACE_DIFF(L)*FLUX_DOWN_CLOUD(L, N_LAYER)           &
     &    +B_CF(L, N_LAYER)*SOURCE_GROUND(L)
      ENDDO
!
!
!     Calculate the upward fluxes using the previous downward fluxes
!     to approximate the scattering term.
      DO I=N_LAYER, N_CLOUD_TOP, -1
        DO L=1, N_PROFILE
!
!         Propagate upward fluxes through the layer.
          FLUX_PROPAGATED=T_FREE(L, I)*FLUX_UP(L, I)+S_UP_FREE(L, I)    &
     &      +R_FREE(L, I)*FLUX_DOWN(L, I-1)
          FLUX_PROPAGATED_CLOUD=T_CLOUD(L, I)*FLUX_UP_CLOUD(L, I)       &
     &      +S_UP_CLOUD(L, I)+R_CLOUD(L, I)*FLUX_DOWN_CLOUD(L, I-1)
!         Transfer upward fluxes across the interface.
          FLUX_UP(L, I-1)=B_FF(L, I-1)*FLUX_PROPAGATED                  &
     &      +B_FC(L, I-1)*FLUX_PROPAGATED_CLOUD
          FLUX_UP_CLOUD(L, I-1)=B_CC(L, I-1)*FLUX_PROPAGATED_CLOUD      &
     &      +B_CF(L, I-1)*FLUX_PROPAGATED
!
        ENDDO
      ENDDO
!
!     Continue through the region above clouds.
      DO I=N_CLOUD_TOP-1, 1, -1
        DO L=1, N_PROFILE
          FLUX_UP(L, I-1)=T_FREE(L, I)*FLUX_UP(L,I)+S_UP_FREE(L, I)     &
     &      +R_FREE(L, I)*FLUX_DOWN(L, I-1)
        ENDDO
      ENDDO
!
!
!
!     Calculate the overall flux.
      DO I=0, N_CLOUD_TOP-2
        DO L=1, N_PROFILE
          FLUX_DIFFUSE(L, 2*I+1)=FLUX_UP(L, I)
          FLUX_DIFFUSE(L, 2*I+2)=FLUX_DOWN(L, I)
        ENDDO
      ENDDO
      DO I=N_CLOUD_TOP-1, N_LAYER
        DO L=1, N_PROFILE
          FLUX_DIFFUSE(L, 2*I+1)=FLUX_UP(L, I)+FLUX_UP_CLOUD(L, I)
          FLUX_DIFFUSE(L, 2*I+2)=FLUX_DOWN(L, I)                        &
     &      +FLUX_DOWN_CLOUD(L, I)
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE MIX_APP_SCAT
#endif
