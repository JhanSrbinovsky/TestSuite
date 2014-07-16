#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve the two-stream equations in a triple column.
!
! Method:
!       The atmospheric column is divided into three regions
!       in each layer and the two-stream coefficients are determined
!       for each region. The equations are then solved using
!       appropriate coupling of the fluxes at the boundaries
!       of layers.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TRIPLE_COLUMN(IERR                                     &
!                     Atmospheric properties
     &   , N_PROFILE, N_LAYER                                           &
!                     Two-stream scheme
     &   , I_2STREAM                                                    &
!                     Options for solver
     &   , I_SOLVER, I_SCATTER_METHOD                                   &
!                     Options for equivalent extinction
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE                               &
!                     Spectral region
     &   , ISOLIR                                                       &
!                     Infra-red properties
     &   , DIFF_PLANCK                                                  &
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                              &
!                     Conditions at TOA
     &   , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                        &
!                     Conditions at surface
     &   , DIFFUSE_ALBEDO, DIRECT_ALBEDO, D_PLANCK_FLUX_SURFACE         &
!                     Single scattering properties
     &   , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR                            &
     &   , TAU, OMEGA, PHASE_FNC                                        &
!                     Cloud geometry
     &   , N_CLOUD_TOP                                                  &
     &   , N_CLOUD_TYPE, FRAC_CLOUD                                     &
     &   , N_REGION, I_REGION_CLOUD, FRAC_REGION                        &
     &   , W_FREE, W_CLOUD                                              &
     &   , CLOUD_OVERLAP                                                &
!                     fluxes calculated
     &   , FLUX_DIRECT, FLUX_TOTAL                                      &
!                     flags for clear-sky calculations
     &   , L_CLEAR, I_SOLVER_CLEAR                                      &
!                     Clear-sky fluxes calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                          &
!                     Dimensions of arrays
     &   , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT                    &
     &   , ND_MAX_ORDER, ND_SOURCE_COEFF                                &
     &   , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF                   &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for completely clear layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_SOURCE_COEFF                                               &
!           Size allocated for source coefficients
     &  , ND_CLOUD_TYPE                                                 &
!           Maximum number of types of cloud
     &  , ND_REGION                                                     &
!           Maximum number of cloudy regions
     &  , ND_OVERLAP_COEFF
!           Maximum number of overlap coefficients
!
!     Include header files.
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
#include "solver_pcf3z.h"
#include "spectral_region_pcf3z.h"
#include "scatter_method_pcf3z.h"
#include "cloud_region_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP                                                   &
!           Top cloudy layer
     &  , N_CLOUD_TYPE                                                  &
!           Number of types of clouds
     &  , ISOLIR                                                        &
!           Spectral region
     &  , I_2STREAM                                                     &
!           Two-stream scheme
     &  , I_SOLVER                                                      &
!           Solver used
     &  , I_SOLVER_CLEAR                                                &
!           Solver for clear-sky fluxes
     &  , I_SCATTER_METHOD
!           Method of treating scattering
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLEAR                                                       &
!           Calculate clear-sky fluxes
     &  , L_SCALE_SOLAR                                                 &
!           Flag to scale solar
     &  , L_IR_SOURCE_QUAD
!           Use quadratic source term
!
!     Optical properties:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU_CLR(ND_PROFILE, ND_LAYER_CLR)                             &
!           Clear-sky optical depth
     &  , OMEGA_CLR(ND_PROFILE, ND_LAYER_CLR)                           &
!           Clear-sky albedo of single scattering
     &  , PHASE_FNC_CLR(ND_PROFILE, ND_LAYER_CLR, ND_MAX_ORDER)         &
!           Clear-sky phase function
     &  , TAU(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)            &
!           Optical depth
     &  , OMEGA(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)          &
!           Albedo of single scattering
     &  , PHASE_FNC(ND_PROFILE, ID_CT: ND_LAYER                         &
     &      , ND_MAX_ORDER, 0: ND_CLOUD_TYPE)
!           Phase function
!
!     Cloud geometry:
      INTEGER, INTENT(IN) ::                                            &
     &    N_REGION                                                      &
!           Number of cloudy regions
     &  , I_REGION_CLOUD(ND_CLOUD_TYPE)
!           Regions in which types of clouds fall
      REAL  (Real64), INTENT(IN) ::                                     &
     &    W_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)                          &
!           Cloudy fractions in each layer
     &  , W_FREE(ND_PROFILE, ID_CT: ND_LAYER)                           &
!           Clear sky fractions in each layer
     &  , FRAC_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_CLOUD_TYPE)        &
!           Fractions of different types of cloud
     &  , CLOUD_OVERLAP(ND_PROFILE, ID_CT-1: ND_LAYER                   &
     &      , ND_OVERLAP_COEFF)                                         &
!           Energy transfer coefficients
     &  , FRAC_REGION(ND_PROFILE, ID_CT: ND_LAYER, ND_REGION)
!           Fractions of total cloud occupied by each region
      REAL  (Real64), INTENT(IN) ::                                     &
     &    SEC_0(ND_PROFILE)                                             &
!           Secant of solar zenith angle
     &  , DIFFUSE_ALBEDO(ND_PROFILE)                                    &
!           Diffuse albedo
     &  , DIRECT_ALBEDO(ND_PROFILE)                                     &
!           Direct albedo
     &  , FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident total flux
     &  , FLUX_INC_DIRECT(ND_PROFILE)                                   &
!           Incident direct flux
     &  , DIFF_PLANCK(ND_PROFILE, ND_LAYER)                             &
!           Change in Planckian function
     &  , D_PLANCK_FLUX_SURFACE(ND_PROFILE)                             &
!           Difference in Planckian fluxes between the surface
!           and the air above
     &  , ADJUST_SOLAR_KE(ND_PROFILE, ND_LAYER)                         &
!           Adjustment of solar beam with equivalent extinction
     &  , DIFF_PLANCK_2(ND_PROFILE, ND_LAYER)
!             2x2nd difference of Planckian
!
!     Fluxes calculated
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_PROFILE, 0: ND_LAYER)                          &
!           Direct flux
     &  , FLUX_TOTAL(ND_PROFILE, 2*ND_LAYER+2)                          &
!           Long flux vector
     &  , FLUX_DIRECT_CLEAR(ND_PROFILE, 0: ND_LAYER)                    &
!           Clear direct flux
     &  , FLUX_TOTAL_CLEAR(ND_PROFILE, 2*ND_LAYER+2)
!           Clear total flux
!
!
!
!     Local variabales.
      INTEGER                                                           &
     &    N_SOURCE_COEFF                                                &
!           Number of source coefficients
     &  , I                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , N_TOP
!           Top-most layer for calculation
!
!
!     Clear-sky coefficients:
      REAL  (Real64) ::                                                 &
     &    TRANS(ND_PROFILE, ND_LAYER, ND_REGION)                        &
!           Transmission coefficients
     &  , REFLECT(ND_PROFILE, ND_LAYER, ND_REGION)                      &
!           Reflection coefficients
     &  , TRANS_0(ND_PROFILE, ND_LAYER, ND_REGION)                      &
!           Direct transmission coefficients
     &  , SOURCE_COEFF(ND_PROFILE, ND_LAYER                             &
     &      , ND_SOURCE_COEFF, ND_REGION)                               &
!           Source coefficients
     &  , S_DOWN(ND_PROFILE, ND_LAYER, ND_REGION)                       &
!           Free downward source
     &  , S_UP(ND_PROFILE, ND_LAYER, ND_REGION)                         &
!           Free upward source
     &  , S_DOWN_CLEAR(ND_PROFILE, ND_LAYER)                            &
!           Clear downward source
     &  , S_UP_CLEAR(ND_PROFILE, ND_LAYER)
!           Clear upward source
!
!     Source functions at the ground
      REAL  (Real64) ::                                                 &
     &    SOURCE_FLUX_GROUND(ND_PROFILE, ND_REGION)                     &
!           Source of flux from ground
     &  , FLUX_DIRECT_GROUND(ND_PROFILE, ND_REGION)
!           Direct flux at ground in each region
!
!
!     Functions called:
      INTEGER                                                           &
     &    SET_N_SOURCE_COEFF
!           Function to set number of source coefficients
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    TWO_COEFF_REGION, TWO_COEFF_REGION_FAST_LW                    &
     &  , IR_SOURCE, TRIPLE_SOLAR_SOURCE                                &
     &  , SOLVER_TRIPLE, SOLVER_TRIPLE_APP_SCAT                         &
     &  , CLEAR_SUPPLEMENT
!
!
!
!     Set the number of source coefficients for the approximation
! DEPENDS ON: set_n_source_coeff
      N_SOURCE_COEFF=SET_N_SOURCE_COEFF(ISOLIR, L_IR_SOURCE_QUAD)
!
!
      IF (I_SCATTER_METHOD == IP_SCATTER_FULL) THEN
! DEPENDS ON: two_coeff_region
        CALL TWO_COEFF_REGION(IERR                                      &
     &    , N_PROFILE, N_LAYER, N_CLOUD_TOP                             &
     &    , I_2STREAM, L_IR_SOURCE_QUAD, N_SOURCE_COEFF                 &
     &    , N_CLOUD_TYPE, FRAC_CLOUD                                    &
     &    , N_REGION, I_REGION_CLOUD, FRAC_REGION                       &
     &    , PHASE_FNC_CLR, OMEGA_CLR, TAU_CLR                           &
     &    , PHASE_FNC, OMEGA, TAU                                       &
     &    , ISOLIR, SEC_0                                               &
     &    , TRANS, REFLECT, TRANS_0, SOURCE_COEFF                       &
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT                   &
     &    , ND_MAX_ORDER, ND_SOURCE_COEFF                               &
     &    , ND_CLOUD_TYPE, ND_REGION                                    &
     &    )
      ELSE IF (I_SCATTER_METHOD == IP_NO_SCATTER_ABS) THEN
! DEPENDS ON: two_coeff_region_fast_lw
        CALL TWO_COEFF_REGION_FAST_LW(IERR                              &
     &    , N_PROFILE, N_LAYER, N_CLOUD_TOP                             &
     &    , L_IR_SOURCE_QUAD, N_SOURCE_COEFF                            &
     &    , N_CLOUD_TYPE, FRAC_CLOUD                                    &
     &    , N_REGION, I_REGION_CLOUD, FRAC_REGION                       &
     &    , TAU_CLR, TAU                                                &
     &    , ISOLIR                                                      &
     &    , TRANS, REFLECT, SOURCE_COEFF                                &
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_SOURCE_COEFF  &
     &    , ND_CLOUD_TYPE, ND_REGION                                    &
     &    )
      ENDIF
      IF (IERR /= I_NORMAL) RETURN
!
!
      IF (ISOLIR == IP_INFRA_RED) THEN
!
        DO K=1, N_REGION
          IF (K == IP_REGION_CLEAR) THEN
            N_TOP=1
          ELSE
            N_TOP=N_CLOUD_TOP
          ENDIF
!
! DEPENDS ON: ir_source
          CALL IR_SOURCE(N_PROFILE, N_TOP, N_LAYER                      &
     &      , SOURCE_COEFF(1, 1, 1, K), DIFF_PLANCK                     &
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                           &
     &      , S_DOWN(1, 1, K), S_UP(1, 1, K)                            &
     &      , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                     &
     &      )
        ENDDO
!
!
!       Weight the source functions by the area fractions, but
!       save the clear-sky fractions for diagnostic use if
!       required.
        IF (L_CLEAR) THEN
          DO I=1, N_LAYER
            DO L=1, N_PROFILE
              S_DOWN_CLEAR(L, I)=S_DOWN(L, I, IP_REGION_CLEAR)
              S_UP_CLEAR(L, I)=S_UP(L, I, IP_REGION_CLEAR)
            ENDDO
          ENDDO
        ENDIF
        DO I=N_CLOUD_TOP, N_LAYER
          DO L=1, N_PROFILE
            S_DOWN(L, I, IP_REGION_CLEAR)                               &
     &        =W_FREE(L, I)*S_DOWN(L, I, IP_REGION_CLEAR)
            S_UP(L, I, IP_REGION_CLEAR)                                 &
     &        =W_FREE(L, I)*S_UP(L, I, IP_REGION_CLEAR)
            S_DOWN(L, I, IP_REGION_STRAT)                               &
     &        =W_CLOUD(L, I)                                            &
     &        *FRAC_REGION(L, I, IP_REGION_STRAT)                       &
     &        *S_DOWN(L, I, IP_REGION_STRAT)
            S_UP(L, I, IP_REGION_STRAT)                                 &
     &        =W_CLOUD(L, I)                                            &
     &        *FRAC_REGION(L, I, IP_REGION_STRAT)                       &
     &        *S_UP(L, I, IP_REGION_STRAT)
            S_DOWN(L, I, IP_REGION_CONV)                                &
     &        =W_CLOUD(L, I)                                            &
     &        *FRAC_REGION(L, I, IP_REGION_CONV)                        &
     &        *S_DOWN(L, I, IP_REGION_CONV)
            S_UP(L, I, IP_REGION_CONV)                                  &
     &        =W_CLOUD(L, I)                                            &
     &        *FRAC_REGION(L, I, IP_REGION_CONV)                        &
     &        *S_UP(L, I, IP_REGION_CONV)
          ENDDO
        ENDDO
!
      ENDIF
!
!
!     Calculate the appropriate source terms for the solar: cloudy
!     and clear properties are both needed here.
!
      IF (ISOLIR == IP_SOLAR) THEN
!
! DEPENDS ON: triple_solar_source
        CALL TRIPLE_SOLAR_SOURCE(N_PROFILE, N_LAYER, N_CLOUD_TOP        &
     &    , N_REGION, FLUX_INC_DIRECT                                   &
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
     &    , TRANS_0, SOURCE_COEFF                                       &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 1)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 2)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 3)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 4)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 5)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 6)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 7)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 8)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 9)                                &
     &    , FLUX_DIRECT, FLUX_DIRECT_GROUND                             &
     &    , S_UP, S_DOWN                                                &
     &    , ND_PROFILE, ND_LAYER, ID_CT, ND_SOURCE_COEFF, ND_REGION     &
     &    )
      ENDIF
!
!       Set the partitioned source functions at the ground.
        IF (ISOLIR == IP_SOLAR) THEN
          DO L=1, N_PROFILE
            SOURCE_FLUX_GROUND(L, IP_REGION_CLEAR)                      &
     &        =(DIRECT_ALBEDO(L)-DIFFUSE_ALBEDO(L))                     &
     &        *FLUX_DIRECT_GROUND(L, IP_REGION_CLEAR)
            SOURCE_FLUX_GROUND(L, IP_REGION_STRAT)                      &
     &        =(DIRECT_ALBEDO(L)-DIFFUSE_ALBEDO(L))                     &
     &        *FLUX_DIRECT_GROUND(L, IP_REGION_STRAT)
            SOURCE_FLUX_GROUND(L, IP_REGION_CONV)                       &
     &        =(DIRECT_ALBEDO(L)-DIFFUSE_ALBEDO(L))                     &
     &        *FLUX_DIRECT_GROUND(L, IP_REGION_CONV)
          ENDDO
        ELSE
          DO L=1, N_PROFILE
            SOURCE_FLUX_GROUND(L, IP_REGION_CLEAR)                      &
     &        =CLOUD_OVERLAP(L, N_LAYER, 10)                            &
     &        *(1.0E+00_Real64-DIFFUSE_ALBEDO(L))                       &
     &        *D_PLANCK_FLUX_SURFACE(L)
            SOURCE_FLUX_GROUND(L, IP_REGION_STRAT)                      &
     &        =CLOUD_OVERLAP(L, N_LAYER, 13)                            &
     &        *(1.0E+00_Real64-DIFFUSE_ALBEDO(L))                       &
     &        *D_PLANCK_FLUX_SURFACE(L)
            SOURCE_FLUX_GROUND(L, IP_REGION_CONV)                       &
     &        =CLOUD_OVERLAP(L, N_LAYER, 16)                            &
     &        *(1.0E+00_Real64-DIFFUSE_ALBEDO(L))                       &
     &        *D_PLANCK_FLUX_SURFACE(L)
          ENDDO
        ENDIF
!
!
!
      SELECT CASE (i_solver)

      CASE (IP_solver_triple)
!
! DEPENDS ON: solver_triple
        CALL SOLVER_TRIPLE(N_PROFILE, N_LAYER, N_CLOUD_TOP              &
     &    , TRANS(1, 1, IP_REGION_CLEAR)                                &
     &    , REFLECT(1, 1, IP_REGION_CLEAR)                              &
     &    , S_DOWN(1, 1, IP_REGION_CLEAR)                               &
     &    , S_UP(1, 1, IP_REGION_CLEAR)                                 &
     &    , TRANS(1, 1, IP_REGION_STRAT)                                &
     &    , REFLECT(1, 1, IP_REGION_STRAT)                              &
     &    , S_DOWN(1, 1, IP_REGION_STRAT)                               &
     &    , S_UP(1, 1, IP_REGION_STRAT)                                 &
     &    , TRANS(1, 1, IP_REGION_CONV)                                 &
     &    , REFLECT(1, 1, IP_REGION_CONV)                               &
     &    , S_DOWN(1, 1, IP_REGION_CONV)                                &
     &    , S_UP(1, 1, IP_REGION_CONV)                                  &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 1)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 2)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 3)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 4)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 5)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 6)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 7)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 8)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 9)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 10)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 11)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 12)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 13)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 14)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 15)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 16)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 17)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 18)                               &
     &    , FLUX_INC_DOWN                                               &
     &    , SOURCE_FLUX_GROUND(1, IP_REGION_CLEAR)                      &
     &    , SOURCE_FLUX_GROUND(1, IP_REGION_STRAT)                      &
     &    , SOURCE_FLUX_GROUND(1, IP_REGION_CONV)                       &
     &    , DIFFUSE_ALBEDO                                              &
     &    , FLUX_TOTAL                                                  &
     &    , ND_PROFILE, ND_LAYER, ID_CT                                 &
     &    )
!
      CASE (IP_solver_triple_hogan)
!
! DEPENDS ON: solver_triple_hogan
        CALL solver_triple_hogan(N_PROFILE, N_LAYER, N_CLOUD_TOP        &
     &    , TRANS(1, 1, IP_REGION_CLEAR)                                &
     &    , REFLECT(1, 1, IP_REGION_CLEAR)                              &
     &    , S_DOWN(1, 1, IP_REGION_CLEAR)                               &
     &    , S_UP(1, 1, IP_REGION_CLEAR)                                 &
     &    , TRANS(1, 1, IP_REGION_STRAT)                                &
     &    , REFLECT(1, 1, IP_REGION_STRAT)                              &
     &    , S_DOWN(1, 1, IP_REGION_STRAT)                               &
     &    , S_UP(1, 1, IP_REGION_STRAT)                                 &
     &    , TRANS(1, 1, IP_REGION_CONV)                                 &
     &    , REFLECT(1, 1, IP_REGION_CONV)                               &
     &    , S_DOWN(1, 1, IP_REGION_CONV)                                &
     &    , S_UP(1, 1, IP_REGION_CONV)                                  &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 1)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 2)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 3)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 4)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 5)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 6)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 7)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 8)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 9)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 10)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 11)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 12)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 13)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 14)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 15)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 16)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 17)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 18)                               &
     &    , FLUX_INC_DOWN                                               &
     &    , SOURCE_FLUX_GROUND(1, IP_REGION_CLEAR)                      &
     &    , SOURCE_FLUX_GROUND(1, IP_REGION_STRAT)                      &
     &    , SOURCE_FLUX_GROUND(1, IP_REGION_CONV)                       &
     &    , DIFFUSE_ALBEDO                                              &
     &    , FLUX_TOTAL                                                  &
     &    , ND_PROFILE, ND_LAYER, ID_CT                                 &
     &    )
!
      CASE (IP_solver_triple_app_scat)
!
! DEPENDS ON: solver_triple_app_scat
        CALL SOLVER_TRIPLE_APP_SCAT(N_PROFILE, N_LAYER, N_CLOUD_TOP     &
     &    , TRANS(1, 1, IP_REGION_CLEAR)                                &
     &    , REFLECT(1, 1, IP_REGION_CLEAR)                              &
     &    , S_DOWN(1, 1, IP_REGION_CLEAR)                               &
     &    , S_UP(1, 1, IP_REGION_CLEAR)                                 &
     &    , TRANS(1, 1, IP_REGION_STRAT)                                &
     &    , REFLECT(1, 1, IP_REGION_STRAT)                              &
     &    , S_DOWN(1, 1, IP_REGION_STRAT)                               &
     &    , S_UP(1, 1, IP_REGION_STRAT)                                 &
     &    , TRANS(1, 1, IP_REGION_CONV)                                 &
     &    , REFLECT(1, 1, IP_REGION_CONV)                               &
     &    , S_DOWN(1, 1, IP_REGION_CONV)                                &
     &    , S_UP(1, 1, IP_REGION_CONV)                                  &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 1)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 2)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 3)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 4)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 5)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 6)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 7)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 8)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 9)                                &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 10)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 11)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 12)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 13)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 14)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 15)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 16)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 17)                               &
     &    , CLOUD_OVERLAP(1, ID_CT-1, 18)                               &
     &    , FLUX_INC_DOWN                                               &
     &    , SOURCE_FLUX_GROUND(1, IP_REGION_CLEAR)                      &
     &    , SOURCE_FLUX_GROUND(1, IP_REGION_STRAT)                      &
     &    , SOURCE_FLUX_GROUND(1, IP_REGION_CONV)                       &
     &    , DIFFUSE_ALBEDO                                              &
     &    , FLUX_TOTAL                                                  &
     &    , ND_PROFILE, ND_LAYER, ID_CT                                 &
     &    )
!
      CASE DEFAULT
!
        WRITE(IU_ERR, '(/A)')                                           &
     &    '***Error: The solver specified is not valid here.'
        IERR=I_ERR_FATAL
        RETURN
!
      END SELECT
!
!
!
      IF (L_CLEAR) THEN
!
! DEPENDS ON: clear_supplement
        CALL CLEAR_SUPPLEMENT(IERR, N_PROFILE, N_LAYER, I_SOLVER_CLEAR  &
     &    , TRANS(1, 1, IP_REGION_CLEAR)                                &
     &    , REFLECT(1, 1, IP_REGION_CLEAR)                              &
     &    , TRANS_0(1, 1, IP_REGION_CLEAR)                              &
     &    , SOURCE_COEFF(1, 1, 1, IP_REGION_CLEAR)                      &
     &    , ISOLIR, FLUX_INC_DIRECT, FLUX_INC_DOWN                      &
     &    , S_DOWN_CLEAR, S_UP_CLEAR                                    &
     &    , DIFFUSE_ALBEDO, DIRECT_ALBEDO                               &
     &    , D_PLANCK_FLUX_SURFACE                                       &
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
     &    , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                         &
     &    , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                       &
     &    )
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE TRIPLE_COLUMN
#endif
