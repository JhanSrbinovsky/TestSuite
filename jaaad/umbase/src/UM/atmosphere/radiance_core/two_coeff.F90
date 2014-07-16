#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate coefficients in the two-stream equations.
!
! Method:
!       The basic two-stream coefficients in the differential equations
!       are calculated. These are then used to determine the
!       transmission and reflection coefficients. Coefficients for
!       determining the solar or infra-red source terms are calculated.
!
! Current owner of code: James Manners
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TWO_COEFF(IERR                                         &
     &   , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                       &
     &   , I_2STREAM, L_IR_SOURCE_QUAD                                  &
     &   , ASYMMETRY, OMEGA, TAU                                        &
     &   , ISOLIR, SEC_0                                                &
     &   , TRANS, REFLECT, TRANS_0                                      &
     &   , SOURCE_COEFF                                                 &
     &   , ND_PROFILE                                                   &
     &   , ID_OP_LT, ID_OP_LB, ID_TRS_LT, ID_TRS_LB                     &
     &   , ND_SOURCE_COEFF                                              &
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
     &  , ID_OP_LT                                                      &
!           Topmost declared layer for optical properties
     &  , ID_OP_LB                                                      &
!           Bottom declared layer for optical properties
     &  , ID_TRS_LT                                                     &
!           Topmost declared layer for transmission coefficients
     &  , ID_TRS_LB                                                     &
!           Bottom declared layer for transmission coefficients
     &  , ND_SOURCE_COEFF
!           Size allocated for source coefficients
!
!     Include header files.
#include "c_kinds.h"
#include "spectral_region_pcf3z.h"
#include "error_pcf3z.h"
!
!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , I_LAYER_FIRST                                                 &
!           First layer to consider
     &  , I_LAYER_LAST                                                  &
!           Last layer to consider
     &  , ISOLIR                                                        &
!           Spectral region
     &  , I_2STREAM
!           Two stream scheme
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Use a quadratic source function
!
!     Optical properties of layer:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ASYMMETRY(ND_PROFILE, ID_OP_LT: ID_OP_LB)                     &
!           Asymmetry factor
     &  , OMEGA(ND_PROFILE, ID_OP_LT: ID_OP_LB)                         &
!           Albedo of single scattering
     &  , TAU(ND_PROFILE, ID_OP_LT: ID_OP_LB)
!           Optical depth
!
!     Solar beam
      REAL  (Real64), INTENT(IN) ::                                     &
     &    SEC_0(ND_PROFILE)
!           Secant of zenith angle
!
!
!     Coefficients in the two-stream equations:
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    TRANS(ND_PROFILE, ID_TRS_LT: ID_TRS_LB)                       &
!           Diffuse transmission coefficient
     &  , REFLECT(ND_PROFILE, ID_TRS_LT: ID_TRS_LB)                     &
!           Diffuse reflection coefficient
     &  , TRANS_0(ND_PROFILE, ID_TRS_LT: ID_TRS_LB)                     &
!           Direct transmission coefficient
     &  , SOURCE_COEFF(ND_PROFILE, ID_TRS_LT: ID_TRS_LB                 &
     &      , ND_SOURCE_COEFF)
!           Source coefficients in two-stream equations
!
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
!     Coefficients in the two-stream equations:
      REAL  (Real64) ::                                                 &
     &    LAMBDA(ND_PROFILE, ID_OP_LT: ID_OP_LB)                        &
!           Coefficients in two-stream equations
     &  , SUM(ND_PROFILE, ID_OP_LT: ID_OP_LB)                           &
!           Sum of alpha_1 and alpha_2
     &  , DIFF(ND_PROFILE, ID_OP_LT: ID_OP_LB)                          &
!           Difference of alpha_1 and alpha_2
     &  , GAMMA_UP(ND_PROFILE, ID_OP_LT: ID_OP_LB)                      &
!           Basic solar coefficient for upward radiation
     &  , GAMMA_DOWN(ND_PROFILE, ID_OP_LT: ID_OP_LB)
!           Basic solar coefficient for downward radiation
!
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    TWO_COEFF_BASIC, SOLAR_COEFFICIENT_BASIC                      &
     &  , TRANS_SOURCE_COEFF
!
!
!
!     Calculate the basic two-stream coefficients. (The single
!     scattering albedo has already been perturbed away from 1 in
!     SINGLE_SCATTERING.)
! DEPENDS ON: two_coeff_basic
      CALL TWO_COEFF_BASIC(IERR                                         &
     &  , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                        &
     &  , I_2STREAM                                                     &
     &  , ASYMMETRY, OMEGA                                              &
     &  , SUM, DIFF                                                     &
     &  , ND_PROFILE, ID_OP_LT, ID_OP_LB                                &
     &  )
      IF (IERR /= I_NORMAL) THEN
        RETURN
      ENDIF
!
!     LAMBDA is now calculated.
      DO I=I_LAYER_FIRST, I_LAYER_LAST
        DO L=1, N_PROFILE
          LAMBDA(L, I)=SQRT(SUM(L, I)*DIFF(L, I))
        ENDDO
      ENDDO
!
!
!     Calculate the basic coefficients for the solar source terms.
      IF (ISOLIR == IP_SOLAR) THEN
!       LAMBDA may be perturbed by this routine to avoid
!       ill-conditioning for the singular zenith angle.
! DEPENDS ON: solar_coefficient_basic
        CALL SOLAR_COEFFICIENT_BASIC(IERR                               &
     &    , N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST                      &
     &    , OMEGA, ASYMMETRY, SEC_0                                     &
     &    , I_2STREAM                                                   &
     &    , SUM, DIFF, LAMBDA                                           &
     &    , GAMMA_UP, GAMMA_DOWN                                        &
     &    , ND_PROFILE, ID_OP_LT, ID_OP_LB                              &
     &    )
        IF (IERR /= I_NORMAL) RETURN
      ENDIF
!
!
!     Determine the transmission and reflection coefficients.
! DEPENDS ON: trans_source_coeff
      CALL TRANS_SOURCE_COEFF(N_PROFILE, I_LAYER_FIRST, I_LAYER_LAST    &
     &  , ISOLIR, L_IR_SOURCE_QUAD                                      &
     &  , TAU, SUM, DIFF, LAMBDA, SEC_0                                 &
     &  , GAMMA_UP, GAMMA_DOWN                                          &
     &  , TRANS, REFLECT, TRANS_0, SOURCE_COEFF                         &
     &  , ND_PROFILE                                                    &
     &  , ID_OP_LT, ID_OP_LB, ID_TRS_LT, ID_TRS_LB                      &
     &  , ND_SOURCE_COEFF                                               &
     &  )
!
!
!
      RETURN
      END SUBROUTINE TWO_COEFF
#endif
