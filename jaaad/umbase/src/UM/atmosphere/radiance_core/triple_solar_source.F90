#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the solar solar terms in a triple column.
!
! Method:
!       The direct beam is calculated by propagating down through
!       the column. These direct fluxes are used to define the
!       source terms in each layer.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TRIPLE_SOLAR_SOURCE(N_PROFILE, N_LAYER, N_CLOUD_TOP    &
     &   , N_REGION, FLUX_INC_DIRECT                                    &
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE                               &
     &   , TRANS_0, SOURCE_COEFF                                        &
     &   , V11, V12, V13, V21, V22, V23, V31, V32, V33                  &
     &   , FLUX_DIRECT                                                  &
     &   , FLUX_DIRECT_GROUND                                           &
     &   , S_UP, S_DOWN                                                 &
     &   , ND_PROFILE, ND_LAYER, ID_CT, ND_SOURCE_COEFF, ND_REGION      &
     &   )
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
     &  , ND_SOURCE_COEFF                                               &
!           Size allocated for source coefficients
     &  , ND_REGION
!           Maximum number of cloudy regions
!
!     Include header files.
#include "c_kinds.h"
#include "source_coeff_pointer_pcf3z.h"
#include "cloud_region_pcf3z.h"
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP                                                   &
!           Top cloudy layer
     &  , N_REGION
!           Number of cloudy regions
!
!     Special arrays for equivalent extinction:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_SCALE_SOLAR
!           Scaling applied to solar flux
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ADJUST_SOLAR_KE(ND_PROFILE, ND_LAYER)
!           Adjustment to solar fluxes with equivalent extinction
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FLUX_INC_DIRECT(ND_PROFILE)
!           Incident direct solar flux
!
!     Optical properties:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TRANS_0(ND_PROFILE, ND_LAYER, ND_REGION)                      &
!           Direct transmission
     &  , SOURCE_COEFF(ND_PROFILE, ND_LAYER                             &
     &    , ND_SOURCE_COEFF, ND_REGION)
!           Source coefficients
!
!     Energy transfer coefficients:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    V11(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , V12(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , V13(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , V21(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , V22(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , V23(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , V31(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , V32(ND_PROFILE, ID_CT-1: ND_LAYER)                            &
!           Energy transfer coefficient
     &  , V33(ND_PROFILE, ID_CT-1: ND_LAYER)
!           Energy transfer coefficient
!
!     Calculated direct flux and source terms:
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_PROFILE, 0: ND_LAYER)                          &
!           Overall direct flux
     &  , FLUX_DIRECT_GROUND(ND_PROFILE, ND_REGION)                     &
!           Direct fluxes at ground beneath each region
     &  , S_UP(ND_PROFILE, ND_LAYER, ND_REGION)                         &
!           Upward source functions
     &  , S_DOWN(ND_PROFILE, ND_LAYER, ND_REGION)
!           Downward source functions
!
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , K
!           Loop variable
!
      REAL  (Real64) ::                                                 &
     &    SOLAR_TOP(ND_PROFILE, ND_REGION)                              &
!           Solar fluxes at top of layer
     &  , SOLAR_BASE(ND_PROFILE, ND_REGION)
!           Solar fluxes at base of layer
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
     &        =FLUX_DIRECT(L, I-1)*TRANS_0(L, I, IP_REGION_CLEAR)       &
     &        *ADJUST_SOLAR_KE(L, I)
            S_UP(L, I, IP_REGION_CLEAR)                                 &
     &        =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, IP_REGION_CLEAR)     &
     &        *FLUX_DIRECT(L, I-1)
            S_DOWN(L, I, IP_REGION_CLEAR)                               &
     &        =(SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN, IP_REGION_CLEAR)  &
     &        -TRANS_0(L, I, IP_REGION_CLEAR))*FLUX_DIRECT(L, I-1)      &
     &        +FLUX_DIRECT(L, I)
          ENDDO
        ENDDO
!
      ELSE
!
        DO I=1, N_CLOUD_TOP-1
          DO L=1, N_PROFILE
            FLUX_DIRECT(L, I)                                           &
     &        =FLUX_DIRECT(L, I-1)*TRANS_0(L, I, IP_REGION_CLEAR)
            S_UP(L, I, IP_REGION_CLEAR)                                 &
     &        =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, IP_REGION_CLEAR)     &
     &        *FLUX_DIRECT(L, I-1)
            S_DOWN(L, I, IP_REGION_CLEAR)                               &
     &        =SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN, IP_REGION_CLEAR)   &
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
        SOLAR_BASE(L, IP_REGION_CLEAR)=FLUX_DIRECT(L, N_CLOUD_TOP-1)
        SOLAR_BASE(L, IP_REGION_STRAT)=0.0E+00_Real64
        SOLAR_BASE(L, IP_REGION_CONV)=0.0E+00_Real64
      ENDDO
!
!
      DO I=N_CLOUD_TOP, N_LAYER
!
!       Transfer fluxes across the interface.
!
        DO L=1, N_PROFILE
          SOLAR_TOP(L, IP_REGION_CLEAR)                                 &
     &      =V11(L, I-1)*SOLAR_BASE(L, IP_REGION_CLEAR)                 &
     &      +V12(L, I-1)*SOLAR_BASE(L, IP_REGION_STRAT)                 &
     &      +V13(L, I-1)*SOLAR_BASE(L, IP_REGION_CONV)
          SOLAR_TOP(L, IP_REGION_STRAT)                                 &
     &      =V21(L, I-1)*SOLAR_BASE(L, IP_REGION_CLEAR)                 &
     &      +V22(L, I-1)*SOLAR_BASE(L, IP_REGION_STRAT)                 &
     &      +V23(L, I-1)*SOLAR_BASE(L, IP_REGION_CONV)
          SOLAR_TOP(L, IP_REGION_CONV)                                  &
     &      =V31(L, I-1)*SOLAR_BASE(L, IP_REGION_CLEAR)                 &
     &      +V32(L, I-1)*SOLAR_BASE(L, IP_REGION_STRAT)                 &
     &      +V33(L, I-1)*SOLAR_BASE(L, IP_REGION_CONV)
        ENDDO
!
!
!       Propagate the fluxes through the layer:
!
        IF (L_SCALE_SOLAR) THEN
!
          DO K=1, N_REGION
            DO L=1, N_PROFILE
              SOLAR_BASE(L, K)                                          &
     &          =SOLAR_TOP(L, K)                                        &
     &          *TRANS_0(L, I, K)*ADJUST_SOLAR_KE(L, I)
              S_UP(L, I, K)                                             &
     &          =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, K)                 &
     &          *SOLAR_TOP(L, K)
              S_DOWN(L, I, K)                                           &
     &          =(SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN, K)              &
     &          -TRANS_0(L, I, K))*SOLAR_TOP(L, K)                      &
     &          +SOLAR_BASE(L, K)
            ENDDO
          ENDDO
!
        ELSE
!
          DO K=1, N_REGION
            DO L=1, N_PROFILE
              SOLAR_BASE(L, K)=SOLAR_TOP(L, K)                          &
     &          *TRANS_0(L, I, K)
              S_UP(L, I, K)                                             &
     &          =SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP, K)                 &
     &          *SOLAR_TOP(L, K)
              S_DOWN(L, I, K)                                           &
     &          =SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN, K)               &
     &          *SOLAR_TOP(L, K)
            ENDDO
          ENDDO
!
        ENDIF
!
!
!       Calculate the total direct flux.
!
        DO L=1, N_PROFILE
          FLUX_DIRECT(L, I)=SOLAR_BASE(L, IP_REGION_CLEAR)              &
     &      +SOLAR_BASE(L, IP_REGION_STRAT)                             &
     &      +SOLAR_BASE(L, IP_REGION_CONV)
        ENDDO
!
      ENDDO
!
!     Pass the last value at the base of the cloud out.
      DO K=1, N_REGION
        DO L=1, N_PROFILE
          FLUX_DIRECT_GROUND(L, K)=SOLAR_BASE(L, K)
        ENDDO
      ENDDO
!
!
!     Correct the direct flux at the ground for sloping terrain

      IF (L_orog) THEN
         FLUX_DIRECT(1:N_PROFILE, N_LAYER) =                            &
     &      FLUX_DIRECT(1:N_PROFILE, N_LAYER) *                         &
     &      lg_orog_corr(1:N_PROFILE)

         DO K=1, N_REGION
            FLUX_DIRECT_GROUND(1:N_PROFILE, K) =                        &
     &         FLUX_DIRECT_GROUND(1:N_PROFILE, K) *                     &
     &         lg_orog_corr(1:N_PROFILE)

            S_DOWN(1:N_PROFILE, N_LAYER, K) =                           &
     &         S_DOWN(1:N_PROFILE, N_LAYER, K) +                        &
     &         SOLAR_BASE(1:N_PROFILE, K) *                             &
     &         (lg_orog_corr(1:N_PROFILE) - 1.0)

         ENDDO
      ENDIF

!
      RETURN
      END SUBROUTINE TRIPLE_SOLAR_SOURCE
#endif
