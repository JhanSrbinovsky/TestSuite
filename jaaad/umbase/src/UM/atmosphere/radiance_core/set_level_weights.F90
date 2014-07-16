#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set level weights for calculating radiances.
!
! Purpose:
!   This routine yields the weights to be applied to the
!   solution of the equation for the complementary function.
!
! Method:
!
!   The particular integral is evaluated at each viewing level
!   within the current layer and weights are calculated for each
!   unknown.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_LEVEL_WEIGHTS(I                                    &
!                     Basic sizes
     &  , N_PROFILE, LS_TRUNC, MS, N_RED_EIGENSYSTEM                    &
!                     Numerical arrays of spherical terms
     &  , CG_COEFF, MU, EIG_VEC                                         &
!                     Solar variables
     &  , ISOLIR, Z_SOL, MU_0                                           &
!                     Infra-red variables
     &  , Q_0, L_IR_SOURCE_QUAD, Q_1                                    &
!                     Conditioning terms
     &  , UPM_C, K_SOL                                                  &
!                     Optical properies
     &  , TAU, SQS2                                                     &
!                     Levels where radiances are calculated
     &  , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                  &
     &  , L_ASSIGN, I_ASSIGN_LEVEL                                      &
!                     Output variables
     &  , C_YLM, WEIGHT_U                                               &
!                     Dimensions
     &  , ND_PROFILE, ND_VIEWING_LEVEL                                  &
     &  , ND_MAX_ORDER                                                  &
     &  , ND_RED_EIGENSYSTEM, ND_SPH_CF_WEIGHT                          &
     &  )
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_VIEWING_LEVEL                                              &
!           Allocated size for levels where radiances are calculated
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_RED_EIGENSYSTEM                                            &
!           Size allocated for the reduced eigensystem
     &  , ND_SPH_CF_WEIGHT
!           Size allocated for entities to be weighted by the C. F.
!
!     Include header files.
#include "c_kinds.h"
#include "c_pi.h"
#include "spectral_region_pcf3z.h"
!
!     Dummy arguments
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of atmospheric layers
     &  , N_RED_EIGENSYSTEM                                             &
!           Size of the reduced eigensystem
     &  , I
!           Current layer
      INTEGER, INTENT(IN) ::                                            &
     &    MS                                                            &
!           Azimuthal order
     &  , LS_TRUNC
!           The truncating order of the system of equations
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Flag for spectral region
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE)
!           Optical depths of the layers
      REAL  (Real64), INTENT(IN) ::                                     &
     &    MU_0(ND_PROFILE)                                              &
!           Cosine of solar zenith angle
     &  , Z_SOL(ND_PROFILE, LS_TRUNC+1-MS)
!           The direct solar radiance
      REAL  (Real64), INTENT(IN) ::                                     &
     &    Q_0(ND_PROFILE)                                               &
!           Term for thermal particular integral
     &  , Q_1(ND_PROFILE)
!           Term for thermal particular integral
!
      INTEGER, INTENT(IN) ::                                            &
     &    K_SOL(ND_PROFILE)
!           Index of eigenvalue closest to the cosine of the solar
!           zenith angle
      REAL  (Real64), INTENT(IN) ::                                     &
     &    UPM_C(ND_PROFILE, 2*ND_RED_EIGENSYSTEM)
!           Weights for exponentials in conditioning term
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Flag for quadratic source function in the IR
      LOGICAL, INTENT(INOUT) ::                                         &
     &    L_ASSIGN
!           Controlling logical for assigning levels
      INTEGER, INTENT(INOUT) ::                                         &
     &    I_ASSIGN_LEVEL
!           Current level where radiances are to be assigned
      REAL  (Real64), INTENT(IN) ::                                     &
     &    SQS2(ND_PROFILE, 0: ND_MAX_ORDER)                             &
!           S-coefficients
     &  , CG_COEFF(LS_TRUNC+1-MS)                                       &
!           Clebsch-Gordan coefficients
     &  , MU(ND_PROFILE, ND_RED_EIGENSYSTEM)                            &
!           Eigenvaluse of the reduced system
     &  , EIG_VEC(ND_PROFILE, 2*ND_RED_EIGENSYSTEM                      &
     &      , ND_RED_EIGENSYSTEM)
!           Eigenvectors of the full systems for positive eigenvalues
!           (these are scaled by the s-coefficients in the routine
!           EIG_SYS)
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_VIEWING_LEVEL                                               &
!           Number of levels where radiances are calculated
     &  , I_RAD_LAYER(ND_VIEWING_LEVEL)
!           Layers in which radiances are calculated
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FRAC_RAD_LAYER(ND_VIEWING_LEVEL)
!           Fractions below the tops of the layers
!
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    C_YLM(ND_PROFILE, ND_VIEWING_LEVEL, LS_TRUNC+1-MS)            &
!           Coefficients for radiances
     &  , WEIGHT_U(ND_PROFILE, ND_VIEWING_LEVEL, ND_SPH_CF_WEIGHT       &
     &      , 2*ND_RED_EIGENSYSTEM)
!           Weights to be applied to the vector U containing the
!           complementary functions
!
!
!     Local variables
      INTEGER                                                           &
     &    IVM                                                           &
!           Index for u^-
     &  , IVP
!           Index for u^+
      INTEGER                                                           &
     &    LSR                                                           &
!           Reduced polar order
     &  , M1LS                                                          &
!           -1^(l+m)
     &  , K                                                             &
!           Loop variable
     &  , L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    EXP_MINUS(ND_PROFILE, ND_RED_EIGENSYSTEM)                     &
!           Exponentials on viewing levels for negative terms
     &  , EXP_PLUS(ND_PROFILE, ND_RED_EIGENSYSTEM)                      &
!           Exponentials on viewing levels for positive terms
     &  , X_M(ND_PROFILE)                                               &
!           Work array connected with negative exponentials
     &  , X_P(ND_PROFILE)
!           Work array connected with positive exponentials
!
!
!
      DO WHILE (L_ASSIGN)
!
!       Calculate exponential terms for the whole routine for speed.
        DO K=1, N_RED_EIGENSYSTEM
          DO L=1, N_PROFILE
            EXP_MINUS(L, K)                                             &
     &        =EXP(-FRAC_RAD_LAYER(I_ASSIGN_LEVEL)*TAU(L)               &
     &        /MU(L, K))
            EXP_PLUS(L, K)                                              &
     &        =EXP((FRAC_RAD_LAYER(I_ASSIGN_LEVEL)-1.0E+00_Real64)      &
     &        *TAU(L)/MU(L, K))
          ENDDO
        ENDDO
!
!       Add on the particular integral.
        IF (ISOLIR == IP_SOLAR) THEN
          DO L=1, N_PROFILE
            X_M(L)                                                      &
     &        =EXP(-FRAC_RAD_LAYER(I_ASSIGN_LEVEL)*TAU(L)/MU_0(L))
          ENDDO
          DO LSR=1, LS_TRUNC-MS+1
            DO L=1, N_PROFILE
              C_YLM(L, I_ASSIGN_LEVEL, LSR)                             &
     &          =C_YLM(L, I_ASSIGN_LEVEL, LSR)+X_M(L)*Z_SOL(L, LSR)
            ENDDO
          ENDDO
!
!         Add on the homogeneous conditioning solution.
          DO L=1, N_PROFILE
            X_M(L)=UPM_C(L, K_SOL(L))*EXP_MINUS(L, K_SOL(L))
          ENDDO
          DO LSR=1, LS_TRUNC-MS+1
            M1LS=REAL(1-2*MOD((LSR-1),2), Real64)
            DO L=1, N_PROFILE
              C_YLM(L, I_ASSIGN_LEVEL, LSR)                             &
     &          =C_YLM(L, I_ASSIGN_LEVEL, LSR)+X_M(L)*M1LS              &
     &          *EIG_VEC(L, LSR, K_SOL(L))
            ENDDO
          ENDDO
!
        ELSE IF (ISOLIR == IP_INFRA_RED) THEN
!
          IF (MS == 0) THEN
!
            IF (L_IR_SOURCE_QUAD) THEN
!
              DO L=1, N_PROFILE
                C_YLM(L, I_ASSIGN_LEVEL, 1)                             &
     &            =C_YLM(L, I_ASSIGN_LEVEL, 1)                          &
     &            +CG_COEFF(1)*Q_1(L)/SQS2(L, 0)
                C_YLM(L, I_ASSIGN_LEVEL, 2)                             &
     &            =C_YLM(L, I_ASSIGN_LEVEL, 2)                          &
     &            +Q_0(L)+Q_1(L)                                        &
     &            *(FRAC_RAD_LAYER(I_ASSIGN_LEVEL)-0.5E+00_Real64)
              ENDDO
              IF (LS_TRUNC >  1) THEN
                DO L=1, N_PROFILE
                  C_YLM(L, I_ASSIGN_LEVEL, 3)                           &
     &              =C_YLM(L, I_ASSIGN_LEVEL, 3)                        &
     &              *CG_COEFF(2)*Q_1(L)/SQS2(L, 2)
                ENDDO
              ENDIF
!
            ELSE
!
              DO L=1, N_PROFILE
                C_YLM(L, I_ASSIGN_LEVEL, 2)                             &
     &            =C_YLM(L, I_ASSIGN_LEVEL, 2)+Q_0(L)
              ENDDO
!
!             Now add on the homogeneous conditioning solution.
              DO K=1, N_RED_EIGENSYSTEM
                DO L=1, N_PROFILE
                  X_M(L)=UPM_C(L, K)*EXP_MINUS(L, K)
                  X_P(L)=UPM_C(L, K+N_RED_EIGENSYSTEM)*EXP_PLUS(L, K)
                ENDDO
                DO LSR=1, LS_TRUNC+1-MS
                  M1LS=REAL(1-2*MOD(LSR-1, 2), Real64)
!                 Increment subsequent terms.
                  DO L=1, N_PROFILE
                    C_YLM(L, I_ASSIGN_LEVEL, LSR)                       &
     &                =C_YLM(L, I_ASSIGN_LEVEL, LSR)                    &
     &                +(X_M(L)*M1LS+X_P(L))*EIG_VEC(L, LSR, K)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDIF
!
        ENDIF
!
!       Calculate the appropriate weights.
        DO K=1, N_RED_EIGENSYSTEM
!         Variable numbers:
          IVM=K
          IVP=IVM+N_RED_EIGENSYSTEM
          DO LSR=1, LS_TRUNC+1-MS
            M1LS=REAL(1-2*MOD((LSR-1),2), Real64)
            DO L=1, N_PROFILE
              WEIGHT_U(L, I_ASSIGN_LEVEL, LSR, IVM)                     &
     &          =EIG_VEC(L, LSR, K)*EXP_MINUS(L, K)*M1LS
              WEIGHT_U(L, I_ASSIGN_LEVEL, LSR, IVP)                     &
     &          =EIG_VEC(L, LSR, K)*EXP_PLUS(L, K)
            ENDDO
          ENDDO
        ENDDO
!
!       Increment the level for assignments:
        I_ASSIGN_LEVEL=I_ASSIGN_LEVEL+1
        IF (I_ASSIGN_LEVEL <= N_VIEWING_LEVEL) THEN
          IF (I_RAD_LAYER(I_ASSIGN_LEVEL) >  I) L_ASSIGN=.FALSE.
        ELSE
          L_ASSIGN=.FALSE.
        ENDIF
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SET_LEVEL_WEIGHTS
#endif
