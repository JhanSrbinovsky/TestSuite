#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate transmission and reflection coefficients.
!
! Method:
!        Straightforward.
!
! Current owner of code: James Manners
!
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TRANS_SOURCE_COEFF(N_PROFILE                           &
     &   , I_LAYER_FIRST, I_LAYER_LAST                                  &
     &   , ISOLIR, L_IR_SOURCE_QUAD                                     &
     &   , TAU, SUM, DIFF, LAMBDA, SEC_0                                &
     &   , GAMMA_UP, GAMMA_DOWN                                         &
     &   , TRANS, REFLECT, TRANS_0, SOURCE_COEFF                        &
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
#include "source_coeff_pointer_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , I_LAYER_FIRST                                                 &
!           First layer to consider
     &  , I_LAYER_LAST
!           Last layer to consider
!
!     Algorithmic control
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Quadratic source in infra-red
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Spectral region
!
!     Optical properties of the layer
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ID_OP_LT: ID_OP_LB)                           &
!           Optical depths of layers
     &  , SUM(ND_PROFILE, ID_OP_LT: ID_OP_LB)                           &
!           Sum of alpha_1 and alpha_2
     &  , DIFF(ND_PROFILE, ID_OP_LT: ID_OP_LB)                          &
!           Difference of alpha_1 and alpha_2
     &  , LAMBDA(ND_PROFILE, ID_OP_LT: ID_OP_LB)                        &
!           Lambda
     &  , SEC_0(ND_PROFILE)                                             &
!           Secant of solar zenith angle
     &  , GAMMA_UP(ND_PROFILE, ID_OP_LT: ID_OP_LB)                      &
!           Basic solar coefficient for upward radiation
     &  , GAMMA_DOWN(ND_PROFILE, ID_OP_LT: ID_OP_LB)
!           Basic solar coefficient for downward radiation
!
!     Transmission and reflection coefficients and coefficients for
!     source terms.
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    TRANS(ND_PROFILE, ID_TRS_LT: ID_TRS_LB)                       &
!           Diffuse transmission coefficient
     &  , REFLECT(ND_PROFILE, ID_TRS_LT: ID_TRS_LB)                     &
!           Diffuse reflection coefficient
     &  , TRANS_0(ND_PROFILE, ID_TRS_LT: ID_TRS_LB)                     &
!           Direct transmission coefficient
     &  , SOURCE_COEFF(ND_PROFILE, ID_TRS_LT: ID_TRS_LB                 &
     &      , ND_SOURCE_COEFF)
!           Source coefficients
!
!
!     Local variables
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    GAMMA                                                         &
!           Gamma
     &  , EXPONENTIAL                                                   &
!           Exponential of scaled optical depth
     &  , GAMMA2                                                        &
!           Gamma squared
     &  , EXPONENTIAL2
!           Exponential squared
!
!     Variables related to the treatment of ill-conditioning
      REAL  (Real64) ::                                                 &
     &    EPS_R                                                         &
!           The smallest real number such that 1.0-EPS_R is not 1
!           to the computer's precision
     &  , SQ_EPS_R                                                      &
!           The square root of the above
     &  , TMP_INV
!           Temporary work variable
!
!
!
!
!     Set the tolerances used in avoiding ill-conditioning, testing
!     on any variable.
      EPS_R=EPSILON(EXPONENTIAL)
      SQ_EPS_R=SQRT(EPS_R)
!
!     Determine the diffuse transmission and reflection coefficients.
!
      DO I=I_LAYER_FIRST, I_LAYER_LAST
        DO L=1, N_PROFILE
          EXPONENTIAL=EXP(-LAMBDA(L, I)*TAU(L, I))
          EXPONENTIAL2=EXPONENTIAL*EXPONENTIAL
          GAMMA=(SUM(L, I)-LAMBDA(L, I)) / (SUM(L, I)+LAMBDA(L, I))
          GAMMA2=GAMMA*GAMMA
          TMP_INV=1.0_Real64 / ( 1.0_Real64 - EXPONENTIAL2*GAMMA2 )
          TRANS(L, I)=EXPONENTIAL*(1.0_Real64-GAMMA2)*TMP_INV
          REFLECT(L, I)=GAMMA*(1.0_Real64-EXPONENTIAL2)*TMP_INV
        ENDDO
      ENDDO
!
!
!
      IF (ISOLIR == IP_SOLAR) THEN
!
!       Calculate the direct transmission and the source coefficients
!       for the solar beam: in the solar case these are
!       the coefficients which will multiply the direct flux at the
!       top of the layer to give the source terms for the upward
!       diffuse flux and the total downward flux.
!
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
            TRANS_0(L, I)=EXP(-TAU(L, I)*SEC_0(L))
            SOURCE_COEFF(L, I, IP_SCF_SOLAR_UP)                         &
     &        =(GAMMA_UP(L, I)-REFLECT(L, I)                            &
     &        *(1.0_Real64+GAMMA_DOWN(L, I)))                           &
     &        -GAMMA_UP(L, I)*TRANS(L, I)*TRANS_0(L, I)
            SOURCE_COEFF(L, I, IP_SCF_SOLAR_DOWN)=TRANS_0(L, I)         &
     &        *(1.0_Real64+GAMMA_DOWN(L, I)                             &
     &        -GAMMA_UP(L, I)*REFLECT(L, I))                            &
     &        -(1.0_Real64+GAMMA_DOWN(L, I))*TRANS(L, I)
          ENDDO
        ENDDO
!
!
      ELSE IF (ISOLIR == IP_INFRA_RED) THEN
!
!       In the case of infra-red radiation, the first source
!       coefficient holds the multiplier for the first difference
!       of the Planckian function across the layer, and the second
!       that for the second difference.
!
        DO I=I_LAYER_FIRST, I_LAYER_LAST
          DO L=1, N_PROFILE
!
!           A tolerance is added to the numerator and the denomiator
!           to avoid ill-conditioning at small optical depths.
!
#if defined(UM)
            SOURCE_COEFF(L, I, IP_SCF_IR_1D)=(1.0_Real64-TRANS(L, I)    &
     &        +REFLECT(L, I)+SQ_EPS_R)                                  &
     &        /(SQ_EPS_R+TAU(L, I)*SUM(L, I))
#endif
#if defined(STANDARD)
            SOURCE_COEFF(L, I, IP_SCF_IR_1D)                            &
     &        =((1.0_Real64-TRANS(L, I)+REFLECT(L, I))                  &
     &        +(EPS_R/(SQ_EPS_R+TAU(L, I))))                            &
     &        /((EPS_R/(SQ_EPS_R+TAU(L, I)))                            &
     &        +TAU(L, I)*SUM(L, I))
#endif
!
          ENDDO
        ENDDO
!
!
        IF (L_IR_SOURCE_QUAD) THEN
!
!         Quadratic correction to source function.
!         This correction is very ill-conditioned for
!         small optical depths so the asymptotic form is then used.
!
          DO I=I_LAYER_FIRST, I_LAYER_LAST
            DO L=1, N_PROFILE
              IF (TAU(L, I) >  EXP(3.3E-01_Real64*LOG(EPS_R))) THEN
                SOURCE_COEFF(L, I, IP_SCF_IR_2D)                        &
     &            =-2.0_Real64                                          &
     &            *(1.0_Real64-TRANS(L, I)-REFLECT(L, I)                &
     &            +SQ_EPS_R)                                            &
     &            /(DIFF(L, I)*TAU(L, I)+SQ_EPS_R)
              ELSE
                SOURCE_COEFF(L, I, IP_SCF_IR_2D)                        &
     &            =-2.0_Real64+DIFF(L, I)*TAU(L, I)
              ENDIF
              SOURCE_COEFF(L, I, IP_SCF_IR_2D)                          &
     &          =-(1.0_Real64+REFLECT(L, I)+TRANS(L, I)                 &
     &          +SOURCE_COEFF(L, I, IP_SCF_IR_2D))                      &
     &          /(SUM(L, I)*TAU(L, I)+SQ_EPS_R)
            ENDDO
          ENDDO
!
        ENDIF
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE TRANS_SOURCE_COEFF
#endif
