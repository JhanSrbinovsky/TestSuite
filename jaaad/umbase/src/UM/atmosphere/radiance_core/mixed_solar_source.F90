#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the solar source terms in a mixed column.
!
! Method:
!       The direct beam is calculated by propagating down through
!       the column. These direct fluxes are used to `define' the
!       source terms in each layer.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MIXED_SOLAR_SOURCE(N_PROFILE, N_LAYER, N_CLOUD_TOP     &
     &  , FLUX_INC_DIRECT                                               &
     &  , L_SCALE_SOLAR, ADJUST_SOLAR_KE                                &
     &  , TRANS_0_FREE, SOURCE_COEFF_FREE                               &
     &  , G_FF, G_FC, G_CF, G_CC                                        &
     &  , TRANS_0_CLOUD, SOURCE_COEFF_CLOUD                             &
     &  , FLUX_DIRECT                                                   &
     &  , FLUX_DIRECT_GROUND_CLOUD                                      &
     &  , S_UP_FREE, S_DOWN_FREE                                        &
     &  , S_UP_CLOUD, S_DOWN_CLOUD                                      &
     &  , ND_PROFILE, ND_LAYER, ID_CT, ND_SOURCE_COEFF                  &
     &  )
!
!
!
      USE solinc_data, ONLY: lg_orog_corr, L_orog
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_SOURCE_COEFF
!           Size allocated for coefficients in the source function
!
!     Include header files.
#include "c_kinds.h"
#include "source_coeff_pointer_pcf3z.h"
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP
!           Top cloudy layer
!
!     Special arrays for equivalent extinction:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_SCALE_SOLAR
!           Scaling applied to solar flux
      REAL  (Real64), INTENT(IN) ::                                     &
     &     ADJUST_SOLAR_KE(ND_PROFILE, ND_LAYER)
!           Adjustment to solar fluxes with equivalent extinction
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FLUX_INC_DIRECT(ND_PROFILE)
!           Incident direct solar flux
!
!     Clear-sky optical properties:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TRANS_0_FREE(ND_PROFILE, ND_LAYER)                            &
!           Free direct transmission
     &  , SOURCE_COEFF_FREE(ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF)
!           Clear-sky source coefficients
!
!     cloudy optical properties:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TRANS_0_CLOUD(ND_PROFILE, ND_LAYER)                           &
!           Cloudy transmission
     &  , SOURCE_COEFF_CLOUD(ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF)
!           Cloudy reflectance
!
!     Energy transfer coefficients:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    G_FF(ND_PROFILE, ID_CT-1: ND_LAYER)                           &
!           Energy transfer coefficient
     &  , G_FC(ND_PROFILE, ID_CT-1: ND_LAYER)                           &
!           Energy transfer coefficient
     &  , G_CF(ND_PROFILE, ID_CT-1: ND_LAYER)                           &
!           Energy transfer coefficient
     &  , G_CC(ND_PROFILE, ID_CT-1: ND_LAYER)
!           Energy transfer coefficient
!
!     Calculated direct flux and source terms:
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_PROFILE, 0: ND_LAYER)                          &
!           Direct flux
     &  , FLUX_DIRECT_GROUND_CLOUD(ND_PROFILE)                          &
!           Direct cloudy flux at ground
     &  , S_UP_FREE(ND_PROFILE, ND_LAYER)                               &
!           Free upward source function
     &  , S_DOWN_FREE(ND_PROFILE, ND_LAYER)                             &
!           Free downward source function
     &  , S_UP_CLOUD(ND_PROFILE, ND_LAYER)                              &
!           Cloudy upward source function
     &  , S_DOWN_CLOUD(ND_PROFILE, ND_LAYER)
!           Cloudy downward source function
!
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
      REAL  (Real64) ::                                                 &
     &    SOLAR_TOP_FREE(ND_PROFILE)                                    &
!           Free solar flux at top of layer
     &  , SOLAR_TOP_CLOUD(ND_PROFILE)                                   &
!           Cloudy solar flux at top of layer
     &  , SOLAR_BASE_FREE(ND_PROFILE)                                   &
!           Free solar flux at base of layer
     &  , SOLAR_BASE_CLOUD(ND_PROFILE)
!           Cloudy solar flux at base of layer
!
!
!
!     The clear and cloudy direct fluxes are calculated separately
!     and added together to form the total direct flux.
!
!     Set incident fluxes.
      DO L=1, N_PROFILE
        FLUX_DIRECT(L, 0)=FLUX_INC_DIRECT(L)
      ENDDO
!
!     With equivalent extinction the direct solar flux must be
!     corrected.
!
      IF (L_SCALE_SOLAR) THEN
!
        DO I=1, N_CLOUD_TOP-1
          DO L=1, N_PROFILE
            FLUX_DIRECT(L, I)                                           &
     &        =FLUX_DIRECT(L, I-1)*TRANS_0_FREE(L, I)                   &
     &        *ADJUST_SOLAR_KE(L, I)
            S_UP_FREE(L, I)=SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_UP)    &
     &        *FLUX_DIRECT(L, I-1)
            S_DOWN_FREE(L, I)                                           &
     &        =(SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_DOWN)              &
     &        -TRANS_0_FREE(L, I))*FLUX_DIRECT(L, I-1)                  &
     &        +FLUX_DIRECT(L, I)
          ENDDO
        ENDDO
!
      ELSE
!
        DO I=1, N_CLOUD_TOP-1
          DO L=1, N_PROFILE
            FLUX_DIRECT(L, I)                                           &
     &        =FLUX_DIRECT(L, I-1)*TRANS_0_FREE(L, I)
            S_UP_FREE(L, I)=SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_UP)    &
     &        *FLUX_DIRECT(L, I-1)
            S_DOWN_FREE(L, I)                                           &
     &        =SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_DOWN)               &
     &        *FLUX_DIRECT(L, I-1)
          ENDDO
        ENDDO
!
      ENDIF
!
!
!
!     Clear and cloudy region.
!     Initialize partial fluxes:
      DO L=1, N_PROFILE
        SOLAR_BASE_FREE(L)=FLUX_DIRECT(L, N_CLOUD_TOP-1)
        SOLAR_BASE_CLOUD(L)=0.0E+00_Real64
      ENDDO
!
!
      DO I=N_CLOUD_TOP, N_LAYER
!
!       Transfer fluxes across the interface. The use of only one
!       cloudy flux implicitly forces random overlap of different
!       subclouds within the cloudy parts of the layer.
!
        DO L=1, N_PROFILE
          SOLAR_TOP_CLOUD(L)=G_CC(L, I-1)*SOLAR_BASE_CLOUD(L)           &
     &      +G_FC(L, I-1)*SOLAR_BASE_FREE(L)
          SOLAR_TOP_FREE(L)=G_FF(L, I-1)*SOLAR_BASE_FREE(L)             &
     &      +G_CF(L, I-1)*SOLAR_BASE_CLOUD(L)
        ENDDO
!
!
!       Propagate the clear and cloudy fluxes through the layer:
!
        IF (L_SCALE_SOLAR) THEN
!
          DO L=1, N_PROFILE
            SOLAR_BASE_FREE(L)=SOLAR_TOP_FREE(L)                        &
     &        *TRANS_0_FREE(L, I)*ADJUST_SOLAR_KE(L, I)
            SOLAR_BASE_CLOUD(L)=SOLAR_TOP_CLOUD(L)                      &
     &        *TRANS_0_CLOUD(L, I)*ADJUST_SOLAR_KE(L, I)
            S_UP_FREE(L, I)=SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_UP)    &
     &        *SOLAR_TOP_FREE(L)
            S_DOWN_FREE(L, I)                                           &
     &        =(SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_DOWN)              &
     &        -TRANS_0_FREE(L, I))*SOLAR_TOP_FREE(L)                    &
     &        +SOLAR_BASE_FREE(L)
            S_UP_CLOUD(L, I)                                            &
     &        =SOURCE_COEFF_CLOUD(L, I, IP_SCF_SOLAR_UP)                &
     &        *SOLAR_TOP_CLOUD(L)
            S_DOWN_CLOUD(L, I)                                          &
     &        =(SOURCE_COEFF_CLOUD(L, I, IP_SCF_SOLAR_DOWN)             &
     &        -TRANS_0_CLOUD(L, I))*SOLAR_TOP_CLOUD(L)                  &
     &        +SOLAR_BASE_CLOUD(L)
          ENDDO
!
        ELSE
!
          DO L=1, N_PROFILE
            SOLAR_BASE_FREE(L)=SOLAR_TOP_FREE(L)                        &
     &        *TRANS_0_FREE(L, I)
            SOLAR_BASE_CLOUD(L)=SOLAR_TOP_CLOUD(L)                      &
     &        *TRANS_0_CLOUD(L, I)
            S_UP_FREE(L, I)=SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_UP)    &
     &        *SOLAR_TOP_FREE(L)
            S_DOWN_FREE(L, I)                                           &
     &        =SOURCE_COEFF_FREE(L, I, IP_SCF_SOLAR_DOWN)               &
     &        *SOLAR_TOP_FREE(L)
            S_UP_CLOUD(L, I)                                            &
     &        =SOURCE_COEFF_CLOUD(L, I, IP_SCF_SOLAR_UP)                &
     &        *SOLAR_TOP_CLOUD(L)
            S_DOWN_CLOUD(L, I)                                          &
     &        =SOURCE_COEFF_CLOUD(L, I, IP_SCF_SOLAR_DOWN)              &
     &        *SOLAR_TOP_CLOUD(L)
          ENDDO
!
        ENDIF
!
!
!       Calculate the total direct flux.
!
        DO L=1, N_PROFILE
          FLUX_DIRECT(L, I)=SOLAR_BASE_FREE(L)+SOLAR_BASE_CLOUD(L)
        ENDDO
!
      ENDDO
!
!     Pass the last value at the base of the cloud out.
      DO L=1, N_PROFILE
        FLUX_DIRECT_GROUND_CLOUD(L)=SOLAR_BASE_CLOUD(L)
      ENDDO
!
!
!     Correct the direct flux at the ground for sloping terrain

      IF (L_orog) THEN
         FLUX_DIRECT(1:N_PROFILE, N_LAYER) =                            &
     &      FLUX_DIRECT(1:N_PROFILE, N_LAYER) *                         &
     &      lg_orog_corr(1:N_PROFILE)

         FLUX_DIRECT_GROUND_CLOUD(1:N_PROFILE) =                        &
     &      FLUX_DIRECT_GROUND_CLOUD(1:N_PROFILE) *                     &
     &      lg_orog_corr(1:N_PROFILE)

         S_DOWN_FREE(1:N_PROFILE, N_LAYER) =                            &
     &         S_DOWN_FREE(1:N_PROFILE, N_LAYER) +                      &
     &         SOLAR_BASE_FREE(1:N_PROFILE) *                           &
     &         (lg_orog_corr(1:N_PROFILE) - 1.0)

         S_DOWN_CLOUD(1:N_PROFILE, N_LAYER) =                           &
     &         S_DOWN_CLOUD(1:N_PROFILE, N_LAYER) +                     &
     &         SOLAR_BASE_CLOUD(1:N_PROFILE) *                          &
     &         (lg_orog_corr(1:N_PROFILE) - 1.0)
      ENDIF

!
      RETURN
      END SUBROUTINE MIXED_SOLAR_SOURCE
#endif
