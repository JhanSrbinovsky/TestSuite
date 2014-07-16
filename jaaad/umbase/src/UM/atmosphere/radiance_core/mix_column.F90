#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve the two-stream equations in a mixed column.
!
! Method:
!       The two-stream coefficients are calculated in clear regions
!       and in stratiform and convective clouds. From these
!       coefficients transmission and reflection coefficients are
!       determined. The coefficients for convective and stratiform
!       clouds are appropriately mixed to form single cloudy values
!       and an appropriate solver is called.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MIX_COLUMN(IERR                                        &
!                     Atmospheric properties
     &   , N_PROFILE, N_LAYER, K_CLR                                    &
!                     Two-stream scheme
     &   , I_2STREAM                                                    &
!                     Options for solver
     &   , I_SOLVER                                                     &
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
!                     Clear-sky single scattering properties
     &   , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR                            &
     &   , TAU, OMEGA, PHASE_FNC                                        &
!                     Cloud geometry
     &   , N_CLOUD_TOP                                                  &
     &   , N_CLOUD_TYPE, FRAC_CLOUD                                     &
     &   , W_FREE, W_CLOUD                                              &
     &   , CLOUD_OVERLAP                                                &
!                     Calculated fluxes
     &   , FLUX_DIRECT, FLUX_TOTAL                                      &
!                     Flags for clear-sky calculations
     &   , L_CLEAR, I_SOLVER_CLEAR                                      &
!                     Calculated clear-sky fluxes
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                          &
!                     Dimensions of arrays
     &   , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT                    &
     &   , ND_MAX_ORDER, ND_SOURCE_COEFF                                &
     &   , ND_CLOUD_TYPE, ND_OVERLAP_COEFF                              &
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
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for completely clear layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_SOURCE_COEFF                                               &
!           Size allocated for coefficients in the source function
     &  , ND_CLOUD_TYPE                                                 &
!           Size allocated for types of clouds
     &  , ND_OVERLAP_COEFF
!           Size allocated for overlpa coefficients
!
!     Include header files.
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
#include "solver_pcf3z.h"
#include "spectral_region_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP                                                   &
!           Top cloudy layer
     &  , K_CLR                                                         &
!           Index of the clear-sky region
     &  , N_CLOUD_TYPE                                                  &
!           Number of types of clouds
     &  , ISOLIR                                                        &
!           Spectral region
     &  , I_2STREAM                                                     &
!           Two-stream scheme
     &  , I_SOLVER                                                      &
!           Solver used
     &  , I_SOLVER_CLEAR
!           Solver for clear-sky fluxes
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
      REAL  (Real64), INTENT(IN) ::                                     &
     &    W_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)                          &
!           Cloudy fractions in each layer
     &  , W_FREE(ND_PROFILE, ID_CT: ND_LAYER)                           &
!           Clear sky fractions in each layer
     &  , FRAC_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_CLOUD_TYPE)        &
!           Fractions of different types of cloud
     &  , CLOUD_OVERLAP(ND_PROFILE, ID_CT-1: ND_LAYER                   &
     &      , ND_OVERLAP_COEFF)
!           Energy transfer coefficients
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
!           Flux from surface
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
     &  , L
!           Loop variable
!
!     Pointers to sections of the array of overlap coefficients:
!     Here F denotes the clear-sky region and C the cloudy region;
!     the ordering of the final suffix is such that the suffix for
!     the area being left appears last: this is convenient to agree
!     with the documented notation, but is not directly consistent
!     with older versions of the code.
      INTEGER                                                           &
     &    I_OVP_DN_FF                                                   &
!           Pointer to section of the array of overlaps for downward
!           transmission from clear-sky to clear-sky
     &  , I_OVP_DN_FC                                                   &
!           Pointer to section of the array of overlaps for downward
!           transmission from cloud to clear-sky
     &  , I_OVP_DN_CF                                                   &
!           Pointer to section of the array of overlaps for downward
!           transmission from clear-sky to cloud
     &  , I_OVP_DN_CC                                                   &
!           Pointer to section of the array of overlaps for downward
!           transmission from cloud to cloud
     &  , I_OVP_UP_FF                                                   &
!           Pointer to section of the array of overlaps for upward
!           transmission from clear-sky to clear-sky
     &  , I_OVP_UP_FC                                                   &
!           Pointer to section of the array of overlaps for upward
!           transmission from cloud to clear-sky
     &  , I_OVP_UP_CF                                                   &
!           Pointer to section of the array of overlaps for upward
!           transmission from clear-sky to cloud
     &  , I_OVP_UP_CC
!           Pointer to section of the array of overlaps for upward
!           transmission from cloud to cloud
!
!
!     Clear-sky coefficients:
      REAL  (Real64) ::                                                 &
     &    TRANS_FREE(ND_PROFILE, ND_LAYER)                              &
!           Free transmission of layer
     &  , REFLECT_FREE(ND_PROFILE, ND_LAYER)                            &
!           Free reflectance of layer
     &  , TRANS_0_FREE(ND_PROFILE, ND_LAYER)                            &
!           Free direct transmission of layer
     &  , SOURCE_COEFF_FREE(ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF)      &
!           Free source coefficients
     &  , S_DOWN_FREE(ND_PROFILE, ND_LAYER)                             &
!           Free downward source
     &  , S_UP_FREE(ND_PROFILE, ND_LAYER)                               &
!           Free upward source
     &  , S_DOWN_CLEAR(ND_PROFILE, ND_LAYER)                            &
!           Clear downward source
     &  , S_UP_CLEAR(ND_PROFILE, ND_LAYER)
!           Clear upward source
!
!     Cloudy coefficients:
      REAL  (Real64) ::                                                 &
     &    TRANS_CLOUD(ND_PROFILE, ND_LAYER)                             &
!           Cloudy transmission of layer
     &  , REFLECT_CLOUD(ND_PROFILE, ND_LAYER)                           &
!           Cloudy reflectance of layer
     &  , TRANS_0_CLOUD(ND_PROFILE, ND_LAYER)                           &
!           Cloudy direct transmission of layer
     &  , SOURCE_COEFF_CLOUD(ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF)     &
!           Cloudy source coefficients
     &  , S_DOWN_CLOUD(ND_PROFILE, ND_LAYER)                            &
!           Cloudy downward source
     &  , S_UP_CLOUD(ND_PROFILE, ND_LAYER)
!           Cloudy upward source
!
!     Source functions at the surface
      REAL  (Real64) ::                                                 &
     &    SOURCE_GROUND_FREE(ND_PROFILE)                                &
!           Source from ground under clear skies
     &  , SOURCE_GROUND_CLOUD(ND_PROFILE)                               &
!           Source from ground under cloudy skies
     &  , FLUX_DIRECT_GROUND_CLOUD(ND_PROFILE)
!           Direct flux at ground under cloudy skies
!
!     Functions called:
      INTEGER                                                           &
     &    SET_N_SOURCE_COEFF
!           Function to set number of source coefficients
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    TWO_COEFF, TWO_COEFF_CLOUD, IR_SOURCE, MIXED_SOLAR_SOURCE     &
     &  , BAND_SOLVER , MIX_APP_SCAT, CLEAR_SUPPLEMENT
!
!
!
!     Set the pointers to the various types of transition.
      I_OVP_DN_FF=3*K_CLR-2
      I_OVP_DN_FC=K_CLR+1
      I_OVP_DN_CF=4-K_CLR
      I_OVP_DN_CC=7-3*K_CLR
      I_OVP_UP_FF=4+I_OVP_DN_FF
      I_OVP_UP_FC=4+I_OVP_DN_FC
      I_OVP_UP_CF=4+I_OVP_DN_CF
      I_OVP_UP_CC=4+I_OVP_DN_CC
!
!     Calculate the transmission and reflection coefficients and
!     source terms for the clear and cloudy parts of the column
!
!     Set the number of source coefficients for the approximation
! DEPENDS ON: set_n_source_coeff
      N_SOURCE_COEFF=SET_N_SOURCE_COEFF(ISOLIR, L_IR_SOURCE_QUAD)
!
! DEPENDS ON: two_coeff
      CALL TWO_COEFF(IERR                                               &
     &  , N_PROFILE, 1, N_CLOUD_TOP-1                                   &
     &  , I_2STREAM, L_IR_SOURCE_QUAD                                   &
     &  , PHASE_FNC_CLR(1, 1, 1), OMEGA_CLR, TAU_CLR                    &
     &  , ISOLIR, SEC_0                                                 &
     &  , TRANS_FREE, REFLECT_FREE, TRANS_0_FREE                        &
     &  , SOURCE_COEFF_FREE                                             &
     &  , ND_PROFILE, 1, ND_LAYER_CLR, 1, ND_LAYER, ND_SOURCE_COEFF     &
     &  )
! DEPENDS ON: two_coeff
      CALL TWO_COEFF(IERR                                               &
     &  , N_PROFILE, N_CLOUD_TOP, N_LAYER                               &
     &  , I_2STREAM, L_IR_SOURCE_QUAD                                   &
     &  , PHASE_FNC(1, 1, 1, 0), OMEGA(1, 1, 0), TAU(1, 1, 0)           &
     &  , ISOLIR, SEC_0                                                 &
     &  , TRANS_FREE, REFLECT_FREE, TRANS_0_FREE                        &
     &  , SOURCE_COEFF_FREE                                             &
     &  , ND_PROFILE, ID_CT, ND_LAYER, 1, ND_LAYER, ND_SOURCE_COEFF     &
     &  )
      IF (IERR /= I_NORMAL) RETURN
!
!
!     Infra-red source terms depend only on the layer and may be
!     calculated now. Solar terms depend on conditions in cloud
!     in overlying layers and must be calculated later.
!
      IF (ISOLIR == IP_INFRA_RED) THEN
!
! DEPENDS ON: ir_source
        CALL IR_SOURCE(N_PROFILE, 1, N_LAYER                            &
     &    , SOURCE_COEFF_FREE, DIFF_PLANCK                              &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                             &
     &    , S_DOWN_FREE, S_UP_FREE                                      &
     &    , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                       &
     &    )
!
!       If a clear-sky calculation is required these source terms must
!       be stored.
        IF (L_CLEAR) THEN
          DO I=1, N_LAYER
            DO L=1, N_PROFILE
              S_DOWN_CLEAR(L, I)=S_DOWN_FREE(L, I)
              S_UP_CLEAR(L, I)=S_UP_FREE(L, I)
            ENDDO
          ENDDO
        ENDIF
!
!       Scale the sources by the clear-sky fractions in the cloudy
!       layers. In higher layers the clear-sky fraction is 1.
        DO I=N_CLOUD_TOP, N_LAYER
          DO L=1, N_PROFILE
            S_DOWN_FREE(L, I)=W_FREE(L, I)*S_DOWN_FREE(L, I)
            S_UP_FREE(L, I)=W_FREE(L, I)*S_UP_FREE(L, I)
          ENDDO
        ENDDO
!
      ENDIF
!
!
!
!     Repeat the calculation for cloudy regions.
!
!     Clouds are indexed beginning with index 1 in the last
!     dimension of arrays of optical properties.
!
!
! DEPENDS ON: two_coeff_cloud
      CALL TWO_COEFF_CLOUD(IERR                                         &
     &  , N_PROFILE, N_CLOUD_TOP, N_LAYER                               &
     &  , I_2STREAM, L_IR_SOURCE_QUAD, N_SOURCE_COEFF                   &
     &  , N_CLOUD_TYPE, FRAC_CLOUD                                      &
     &  , PHASE_FNC(1, 1, 1, 1), OMEGA(1, 1, 1), TAU(1, 1, 1)           &
     &  , ISOLIR, SEC_0                                                 &
     &  , TRANS_CLOUD, REFLECT_CLOUD, TRANS_0_CLOUD                     &
     &  , SOURCE_COEFF_CLOUD                                            &
     &  , ND_PROFILE, ND_LAYER, ID_CT, ND_MAX_ORDER                     &
     &  , ND_SOURCE_COEFF, ND_CLOUD_TYPE                                &
     &  )
      IF (IERR /= I_NORMAL) RETURN
!
!
      IF (ISOLIR == IP_INFRA_RED) THEN
!
! DEPENDS ON: ir_source
        CALL IR_SOURCE(N_PROFILE, N_CLOUD_TOP, N_LAYER                  &
     &    , SOURCE_COEFF_CLOUD, DIFF_PLANCK                             &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                             &
     &    , S_DOWN_CLOUD, S_UP_CLOUD                                    &
     &    , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                       &
     &    )
!
        DO I=N_CLOUD_TOP, N_LAYER
          DO L=1, N_PROFILE
            S_DOWN_CLOUD(L, I)=W_CLOUD(L, I)*S_DOWN_CLOUD(L, I)
            S_UP_CLOUD(L, I)=W_CLOUD(L, I)*S_UP_CLOUD(L, I)
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
! DEPENDS ON: mixed_solar_source
        CALL MIXED_SOLAR_SOURCE(N_PROFILE, N_LAYER, N_CLOUD_TOP         &
     &    , FLUX_INC_DIRECT                                             &
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
     &    , TRANS_0_FREE, SOURCE_COEFF_FREE                             &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_FF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_CF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_FC)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_CC)                      &
     &    , TRANS_0_CLOUD, SOURCE_COEFF_CLOUD                           &
     &    , FLUX_DIRECT                                                 &
     &    , FLUX_DIRECT_GROUND_CLOUD                                    &
     &    , S_UP_FREE, S_DOWN_FREE                                      &
     &    , S_UP_CLOUD, S_DOWN_CLOUD                                    &
     &    , ND_PROFILE, ND_LAYER, ID_CT, ND_SOURCE_COEFF                &
     &  )
      ENDIF
!
!
!
!     Formulate the matrix equation for the fluxes.
!
      SELECT CASE (i_solver)

      CASE (IP_solver_mix_app_scat)
!
! DEPENDS ON: mix_app_scat
        CALL MIX_APP_SCAT(N_PROFILE, N_LAYER, N_CLOUD_TOP               &
     &    , TRANS_FREE, REFLECT_FREE, S_DOWN_FREE, S_UP_FREE            &
     &    , TRANS_CLOUD, REFLECT_CLOUD                                  &
     &    , S_DOWN_CLOUD, S_UP_CLOUD                                    &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_FF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_CF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_FC)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_CC)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_FF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_FC)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_CF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_CC)                      &
     &    , FLUX_INC_DOWN                                               &
     &    , D_PLANCK_FLUX_SURFACE, DIFFUSE_ALBEDO                       &
     &    , FLUX_TOTAL                                                  &
     &    , ND_PROFILE, ND_LAYER, ID_CT                                 &
     &    )
!
      CASE (IP_solver_mix_direct, IP_solver_mix_direct_hogan)
!
!       Set the partitioned source functions at the ground.
        IF (ISOLIR == IP_SOLAR) THEN
          DO L=1, N_PROFILE
            SOURCE_GROUND_FREE(L)=(DIRECT_ALBEDO(L)                     &
     &        -DIFFUSE_ALBEDO(L))                                       &
     &        *(FLUX_DIRECT(L, N_LAYER)                                 &
     &        -FLUX_DIRECT_GROUND_CLOUD(L))
            SOURCE_GROUND_CLOUD(L)=(DIRECT_ALBEDO(L)                    &
     &        -DIFFUSE_ALBEDO(L))                                       &
     &        *FLUX_DIRECT_GROUND_CLOUD(L)
          ENDDO
        ELSE
          DO L=1, N_PROFILE
            SOURCE_GROUND_FREE(L)                                       &
     &        =CLOUD_OVERLAP(L, N_LAYER, I_OVP_UP_FF)                   &
     &        *DIFFUSE_ALBEDO(L)*D_PLANCK_FLUX_SURFACE(L)
            SOURCE_GROUND_CLOUD(L)                                      &
     &        =CLOUD_OVERLAP(L, N_LAYER, I_OVP_UP_CF)                   &
     &        *DIFFUSE_ALBEDO(L)*D_PLANCK_FLUX_SURFACE(L)
          ENDDO
        ENDIF
!
        IF (i_solver == IP_solver_mix_direct) THEN
! DEPENDS ON: solver_mix_direct
          CALL SOLVER_MIX_DIRECT(N_PROFILE, N_LAYER, N_CLOUD_TOP        &
     &    , TRANS_FREE, REFLECT_FREE, S_DOWN_FREE, S_UP_FREE            &
     &    , TRANS_CLOUD, REFLECT_CLOUD                                  &
     &    , S_DOWN_CLOUD, S_UP_CLOUD                                    &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_FF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_CF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_FC)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_CC)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_FF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_FC)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_CF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_CC)                      &
     &    , FLUX_INC_DOWN                                               &
     &    , SOURCE_GROUND_FREE, SOURCE_GROUND_CLOUD                     &
     &    , DIFFUSE_ALBEDO                                              &
     &    , FLUX_TOTAL                                                  &
     &    , ND_PROFILE, ND_LAYER, ID_CT                                 &
     &    )
!
        ELSE IF (i_solver == IP_solver_mix_direct_hogan) THEN
! DEPENDS ON: solver_mix_direct_hogan
          CALL solver_mix_direct_hogan(N_PROFILE, N_LAYER, N_CLOUD_TOP  &
     &    , TRANS_FREE, REFLECT_FREE, S_DOWN_FREE, S_UP_FREE            &
     &    , TRANS_CLOUD, REFLECT_CLOUD                                  &
     &    , S_DOWN_CLOUD, S_UP_CLOUD                                    &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_FF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_CF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_FC)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_DN_CC)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_FF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_FC)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_CF)                      &
     &    , CLOUD_OVERLAP(1, ID_CT-1, I_OVP_UP_CC)                      &
     &    , FLUX_INC_DOWN                                               &
     &    , SOURCE_GROUND_FREE, SOURCE_GROUND_CLOUD                     &
     &    , DIFFUSE_ALBEDO                                              &
     &    , FLUX_TOTAL                                                  &
     &    , ND_PROFILE, ND_LAYER, ID_CT                                 &
     &    )

        ENDIF
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
     &    , TRANS_FREE, REFLECT_FREE, TRANS_0_FREE, SOURCE_COEFF_FREE   &
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
      END SUBROUTINE MIX_COLUMN
#endif
