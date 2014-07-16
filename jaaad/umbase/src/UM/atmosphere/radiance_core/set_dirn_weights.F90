#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set weights for the C.F. along a direction.
!
! Purpose:
!   The complementary function for the radiation involves unknown
!   coefficients: we set the weights for these coefficients in the
!   current layer here.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_DIRN_WEIGHTS(N_PROFILE                             &
     &  , MS, LS_TRUNC, UP_LM                                           &
     &  , N_DIRECTION, MU_V, AZIM_FACTOR                                &
     &  , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER, I               &
     &  , N_RED_EIGENSYSTEM, MU, EIG_VEC                                &
     &  , ISOLIR, Z_SOL, MU_0                                           &
     &  , L_IR_SOURCE_QUAD, DIFF_PLANCK                                 &
     &  , UPM_C, K_SOL                                                  &
     &  , TAU, OMEGA, PHASE_FNC                                         &
     &  , WEIGHT_U, RADIANCE                                            &
     &  , ND_PROFILE, ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL          &
     &  , ND_RED_EIGENSYSTEM, ND_MAX_ORDER                              &
     &  )
!
!
!
      IMPLICIT NONE
!
!
!     Include header files
#include "c_kinds.h"
#include "c_pi.h"
#include "spectral_region_pcf3z.h"
!
!     Dummy arguments:
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_VIEWING_LEVEL                                              &
!           Size allocated for levels where radiances are calculated
     &  , ND_DIRECTION                                                  &
!           Size allocated for viewing directions
     &  , ND_RED_EIGENSYSTEM                                            &
!           Size allocated for the reduced eigensystem
     &  , ND_MAX_ORDER
!           Size allocated for orders of spherical harmonics
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of atmospheric profiles
     &  , N_DIRECTION                                                   &
!           Number of directions
     &  , N_VIEWING_LEVEL                                               &
!           Number of levels where the radiance is calculated
     &  , I_RAD_LAYER(ND_VIEWING_LEVEL)                                 &
!           Indices of layers containing viewing levels
     &  , N_RED_EIGENSYSTEM
!           Size of the reduced eigensystem
      INTEGER, INTENT(IN) ::                                            &
     &    I                                                             &
!           Current layer
     &  , MS                                                            &
!           Current azimuthal order
     &  , LS_TRUNC                                                      &
!           Order of polar truncation
     &  , ISOLIR
!           Index of spectral region
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    MU_V(ND_PROFILE, ND_DIRECTION)                                &
!           Cosines of polar viewing angles
     &  , AZIM_FACTOR(ND_PROFILE, ND_DIRECTION)                         &
!           Azimuthal factors
     &  , FRAC_RAD_LAYER(ND_VIEWING_LEVEL)                              &
!           Fraction optical depth into its layer of the
!           viewing level
     &  , MU_0(ND_PROFILE)                                              &
!           Cosines of solar zenith angle
     &  , TAU(ND_PROFILE, ND_LAYER)                                     &
!           Optical depths
     &  , OMEGA(ND_PROFILE, ND_LAYER)                                   &
!           Albedos of single scattering
     &  , PHASE_FNC(ND_PROFILE, ND_LAYER, ND_MAX_ORDER)                 &
!           Phase function
     &  , MU(ND_PROFILE, ND_RED_EIGENSYSTEM)                            &
!           Eigenvalues of the reduced eigensystem
     &  , EIG_VEC(ND_PROFILE, 2*ND_RED_EIGENSYSTEM                      &
     &      , ND_RED_EIGENSYSTEM)                                       &
!           Eigenvalues of the full eigensystem scaled by
!           the s-parameters
     &  , Z_SOL(ND_PROFILE, LS_TRUNC+1-MS)
!           Coefficient of the solar source function
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Flag for quadratic IR-sources
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DIFF_PLANCK(ND_PROFILE, ND_LAYER)
!           Differences in the hemispheric Planckian FLUX (bottom-top)
!           across the layer
      INTEGER, INTENT(IN) ::                                            &
     &    K_SOL(ND_PROFILE)
!           Indices of eigenvalues used to restore solar conditioning
      REAL  (Real64), INTENT(IN) ::                                     &
     &    UPM_C(ND_PROFILE, 2*ND_RED_EIGENSYSTEM)
!           Coefficients of homogeneous solution used to restore
!           conditioning
!
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    RADIANCE(ND_PROFILE, ND_VIEWING_LEVEL, ND_DIRECTION)
!           Radiances (to be incremented by the contribution of
!           the particular integral)
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    WEIGHT_U(ND_PROFILE, ND_VIEWING_LEVEL                         &
     &      , ND_DIRECTION, 2*ND_RED_EIGENSYSTEM)
!           Weights for the coefficients in the complementary
!           function
!
!
!     Local variables
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable (points)
     &  , IR                                                            &
!           Loop variable (radiative levels)
     &  , ID                                                            &
!           Loop variable (directions)
     &  , LS                                                            &
!           Loop variable (polar orders)
     &  , LL                                                            &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , II
!           Loop variable
      INTEGER                                                           &
     &    N_LIST_UP                                                     &
!           Numbers of points where the current viewing direction
!           points up
     &  , LIST_UP(ND_PROFILE)                                           &
!           List up points with where the current viewing direction
!           points up
     &  , N_LIST_DOWN                                                   &
!           Numbers of points where the current viewing direction
!           points down
     &  , LIST_DOWN(ND_PROFILE)
!           List up points with where the current viewing direction
!           points up
      REAL  (Real64) ::                                                 &
     &    GEOM_SOLAR(ND_PROFILE)                                        &
!           Geometrical factor for solar radiation
     &  , GEOM_INTEG_M(ND_PROFILE)                                      &
!           Geometrical factor for negative eigenvalues
     &  , GEOM_INTEG_P(ND_PROFILE)                                      &
!           Geometrical factor for positive eigenvalues
     &  , M_SLANT_DEPTH_NEAR(ND_PROFILE)                                &
!           Minus slantwise optical distance between the radiance
!           level and the nearer boundary of the current layer
     &  , M_SLANT_DEPTH_INC(ND_PROFILE)                                 &
!           Minus the increment in the slantwise optical distance
!           between the boundaries of the current layer or partial
!           layer when the viewing level lies within it
     &  , UP_LM(ND_PROFILE, ND_MAX_ORDER+1, ND_DIRECTION)               &
!           Spherical harmonics at a fixed azimuthal order
     &  , LS_SUM_S(ND_PROFILE)                                          &
!           Sum of terms over polar orders in the solar integral
     &  , LS_SUM_P(ND_PROFILE, ND_RED_EIGENSYSTEM)                      &
!           Sum of terms over polar orders in the integral over
!           eigenvalues
     &  , LS_SUM_M(ND_PROFILE, ND_RED_EIGENSYSTEM)                      &
!           Sum of terms over polar orders in the integral over
!           eigenvalues
     &  , TAU_I(ND_PROFILE)                                             &
!           Optical depth of the relevant part of the current layer
     &  , FRAC_TAU_I(ND_PROFILE)                                        &
!           Fractional of the optical depth of the current layer in
!           the relevant part
     &  , TRANS_TOP(ND_PROFILE)                                         &
!           Solar transmission from the top of the current layer to the
!           viewing level within the current layer
     &  , D_MU                                                          &
!           Difference in cosines of directions
     &  , X                                                             &
!           Temporary variable
     &  , M1LSMS
!           -1^(l+m)
!
!     Variables related to the treatment of ill-conditioning
      REAL  (Real64) ::                                                 &
     &    EPS_R                                                         &
!           The smallest real number such that 1.0-EPS_R is not 1
!           to the computer's precision
     &  , SQ_EPS_R                                                      &
!           The square root of the above
     &  , ETA                                                           &
!           The conditioning weight
     &  , ETA_NM
!           The conditioning multiplier applied in the numerator
!
!     Subroutines called
      EXTERNAL                                                          &
     &    EVAL_UPLM
!
!
!
!     Set the tolerances used in avoiding ill-conditioning, testing
!     on any variable.
      EPS_R=EPSILON(MU_0(1))
      SQ_EPS_R=SQRT(EPS_R)
!     Consider each level where radiances are required in turn
!     and calculate appropriate weightings. The expressions for
!     these weightings will look slightly different if the radiance
!     level lies within the current layer.
      DO ID=1, N_DIRECTION
!
!       Assemble the list of points where the direction is upward
!       or downard. Viewing along a horizontal direction is not
!       considered valid (and should be filtered out before this).
        N_LIST_UP=0
        N_LIST_DOWN=0
        DO L=1, N_PROFILE
          IF (MU_V(L, ID) >  0.0E+00_Real64) THEN
            N_LIST_UP=N_LIST_UP+1
            LIST_UP(N_LIST_UP)=L
          ELSE IF (MU_V(L, ID) <  0.0E+00_Real64) THEN
            N_LIST_DOWN=N_LIST_DOWN+1
            LIST_DOWN(N_LIST_DOWN)=L
          ENDIF
        ENDDO
!
!       Sum the terms which vary with the polar order, but
!       are independent of the viewing level.
!       First the contributions to the particular integral:
        IF (ISOLIR == IP_SOLAR) THEN
!
          IF (MS == 0) THEN
!           The zeroth moment of the phase function is not stored
!           because it is 1, but this means that we need some
!           special code to treat the exception.
            DO L=1, N_PROFILE
              LS_SUM_S(L)=UP_LM(L, 1, ID)*Z_SOL(L, 1)
            ENDDO
          ELSE
            DO L=1, N_PROFILE
              LS_SUM_S(L)=PHASE_FNC(L, I, MS)*UP_LM(L, 1, ID)           &
     &          *Z_SOL(L, 1)
            ENDDO
          ENDIF
          DO LS=MS+1, LS_TRUNC
            DO L=1, N_PROFILE
              LS_SUM_S(L)=LS_SUM_S(L)+Z_SOL(L, LS+1-MS)                 &
     &          *PHASE_FNC(L, I, LS)*UP_LM(L, LS+1-MS, ID)
            ENDDO
          ENDDO
!
        ENDIF
!
!
        DO K=1, N_RED_EIGENSYSTEM
          IF (MS == 0) THEN
            DO L=1, N_PROFILE
              LS_SUM_P(L, K)=UP_LM(L, 1, ID)*EIG_VEC(L, 1, K)
              LS_SUM_M(L, K)=LS_SUM_P(L, K)
            ENDDO
          ELSE
            DO L=1, N_PROFILE
              LS_SUM_P(L, K)=PHASE_FNC(L, I, MS)*UP_LM(L, 1, ID)        &
     &          *EIG_VEC(L, 1, K)
              LS_SUM_M(L, K)=LS_SUM_P(L, K)
            ENDDO
          ENDIF
          DO LS=MS+1, LS_TRUNC
            M1LSMS=REAL(1-2*MOD((LS+MS), 2), Real64)
            DO L=1, N_PROFILE
              X=PHASE_FNC(L, I, LS)*UP_LM(L, LS+1-MS, ID)               &
     &          *EIG_VEC(L, LS+1-MS, K)
              LS_SUM_P(L, K)=LS_SUM_P(L, K)+X
              LS_SUM_M(L, K)=LS_SUM_M(L, K)+X*M1LSMS
            ENDDO
          ENDDO
        ENDDO
!
!
        DO IR=1, N_VIEWING_LEVEL
!
!         Initialize the weights:
          DO K=1, 2*N_RED_EIGENSYSTEM
            DO L=1, N_PROFILE
              WEIGHT_U(L, IR, ID, K)=0.0E+00_Real64
            ENDDO
          ENDDO
!
!         Upward radiation:
!         No layers above the viewing level contribute here.
!
          IF (I >= I_RAD_LAYER(IR)) THEN
!
!           Calculate minus the slantwise optical depths to the
!           boundaries of the layer. If the level where the radiance
!           is required lies in the current layer we perform the
!           calculation for a temporary layer reaching from the
!           viewing level to the bottom of the current layer.
            IF (I >  I_RAD_LAYER(IR)) THEN
!             Full layers are required.
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                M_SLANT_DEPTH_NEAR(L)                                   &
     &            =(1.0E+00_Real64-FRAC_RAD_LAYER(IR))                  &
     &            *TAU(L, I_RAD_LAYER(IR))
              ENDDO
              DO II=I_RAD_LAYER(IR)+1, I-1
                DO LL=1, N_LIST_UP
                  L=LIST_UP(LL)
                  M_SLANT_DEPTH_NEAR(L)                                 &
     &              =M_SLANT_DEPTH_NEAR(L)+TAU(L, II)
                ENDDO
              ENDDO
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                M_SLANT_DEPTH_NEAR(L)                                   &
     &            =-M_SLANT_DEPTH_NEAR(L)/MU_V(L, ID)
                M_SLANT_DEPTH_INC(L)=-TAU(L, I)/MU_V(L, ID)
                TAU_I(L)=TAU(L, I)
                FRAC_TAU_I(L)=1.0E+00_Real64
              ENDDO
              IF (ISOLIR == IP_SOLAR) THEN
                DO LL=1, N_LIST_UP
                  L=LIST_UP(LL)
                  TRANS_TOP(L)=1.0E+00_Real64
                ENDDO
              ENDIF
            ELSE IF (I == I_RAD_LAYER(IR)) THEN
!             The viewing level lies in the current layer.
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                M_SLANT_DEPTH_NEAR(L)=0.0E+00_Real64
                FRAC_TAU_I(L)=1.0E+00_Real64-FRAC_RAD_LAYER(IR)
                TAU_I(L)=FRAC_TAU_I(L)*TAU(L, I)
                M_SLANT_DEPTH_INC(L)=-TAU_I(L)/MU_V(L, ID)
              ENDDO
              IF (ISOLIR == IP_SOLAR) THEN
                DO LL=1, N_LIST_UP
                  L=LIST_UP(LL)
                  TRANS_TOP(L)                                          &
     &              =EXP(-FRAC_RAD_LAYER(IR)*TAU(L, I)/MU_0(L))
                ENDDO
              ENDIF
            ENDIF
!
!
            IF (ISOLIR == IP_SOLAR) THEN
!             Set the geometrical terms for the solar integral.
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                GEOM_SOLAR(L)=(MU_0(L)/(MU_0(L)+MU_V(L, ID)))           &
     &            *EXP(M_SLANT_DEPTH_NEAR(L))*(1.0E+00_Real64           &
     &            -EXP(M_SLANT_DEPTH_INC(L)-TAU_I(L)/MU_0(L)))
              ENDDO
!             Add the contribution of the particular integral to the
!             radiance. TRANS_TOP is required to adjust the solar
!             particular integral from its value at the top of the
!             actual layer to its value at the top of the notional
!             layer when the viewing level lies within the current
!             layer.
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                RADIANCE(L, IR, ID)                                     &
     &            =RADIANCE(L, IR, ID)+AZIM_FACTOR(L, ID)               &
     &            *LS_SUM_S(L)*OMEGA(L, I)*GEOM_SOLAR(L)*TRANS_TOP(L)
              ENDDO
            ELSE IF (ISOLIR == IP_INFRA_RED) THEN
              IF (MS == 0) THEN
                IF (L_IR_SOURCE_QUAD) THEN
                  print*, 'Not done'
                ELSE
                  DO LL=1, N_LIST_UP
                    L=LIST_UP(LL)
!                   The azimuthal factor is omitted since it will be 1.
!                   Numerical ill-conditioning can arise when the
!                   optical depth is small, necessitating special
!                   treatment.
                    IF (M_SLANT_DEPTH_INC(L) <  -SQ_EPS_R) THEN
                      X=EXP(M_SLANT_DEPTH_NEAR(L))                      &
     &                *(1.0E+00_Real64-EXP(M_SLANT_DEPTH_INC(L)))       &
     &                /M_SLANT_DEPTH_INC(L)
                    ELSE
!                     Keep the first couple of terms from a power
!                     series.
                      X=-EXP(M_SLANT_DEPTH_NEAR(L))*(1.0E+00_Real64     &
     &                  +0.5E+00_Real64*M_SLANT_DEPTH_INC(L))
                    ENDIF
                    RADIANCE(L, IR, ID)                                 &
     &                =RADIANCE(L, IR, ID)                              &
     &                -(DIFF_PLANCK(L, I)/PI)*X*FRAC_TAU_I(L)           &
     &                /(1.0_Real64-OMEGA(L, I)*PHASE_FNC(L, I, 1))
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
!
!           Determine the contribution from each eigenvalue.
            DO K=1, N_RED_EIGENSYSTEM
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
!               The term for u^+:
!               This may exhibit ill-conditioning, so it is perturbed
!               using L'Hopital's rule (actually we keep two terms in
!               the expansion).
                D_MU=MU(L, K)-MU_V(L, ID)
                ETA=EPS_R/(D_MU+SIGN(SQ_EPS_R, D_MU))
                X=TAU_I(L)/(MU(L, K)*MU_V(L, ID))
                ETA_NM=1.0_Real64-ETA*X*(1.0_Real64+0.5*Real64*X*D_MU)
                GEOM_INTEG_P(L)                                         &
     &            =(MU(L, K)/(D_MU+ETA))*EXP(M_SLANT_DEPTH_NEAR(L))     &
     &            *(EXP(-TAU_I(L)/MU(L, K))*ETA_NM                      &
     &            -EXP(M_SLANT_DEPTH_INC(L)))
                GEOM_INTEG_M(L)                                         &
     &            =(MU(L, K)/(MU(L, K)+MU_V(L, ID)))                    &
     &            *EXP(M_SLANT_DEPTH_NEAR(L))                           &
     &            *(EXP(-(TAU(L, I)-TAU_I(L))/MU(L, K))                 &
     &            -EXP(M_SLANT_DEPTH_INC(L)-TAU(L, I)/MU(L, K)))
              ENDDO
!
!             Combine to form the weights for each element of the
!             solution. Only a poprtion of WEIGHT_U is passed to this
!             routine, so the offsetting over layers takes
!             care of itself.
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                WEIGHT_U(L, IR, ID, K)                                  &
     &            =OMEGA(L, I)*LS_SUM_M(L, K)*GEOM_INTEG_M(L)
                WEIGHT_U(L, IR, ID, K+N_RED_EIGENSYSTEM)                &
     &            =OMEGA(L, I)*LS_SUM_P(L, K)*GEOM_INTEG_P(L)
              ENDDO
            ENDDO
!
          ENDIF
!
!
!         Downward Radiation:
!         No layers below the current viewing level contribute here.
          IF (I <= I_RAD_LAYER(IR)) THEN
!
!           Calculate the slantwise optical depths to the
!           boundaries of the layer. If the observing level lies
!           within the current layer we perform the calculation for
!           a layer reaching from the top of the current layer to
!           the observing level.
            IF (I <  I_RAD_LAYER(IR)) THEN
              DO LL=1, N_LIST_DOWN
                L=LIST_DOWN(LL)
                M_SLANT_DEPTH_NEAR(L)                                   &
     &            =FRAC_RAD_LAYER(IR)*TAU(L, I_RAD_LAYER(IR))
              ENDDO
              DO II=I_RAD_LAYER(IR)-1, I+1, -1
                DO LL=1, N_LIST_DOWN
                  L=LIST_DOWN(LL)
                  M_SLANT_DEPTH_NEAR(L)                                 &
     &              =M_SLANT_DEPTH_NEAR(L)+TAU(L, II)
                ENDDO
              ENDDO
              DO LL=1, N_LIST_DOWN
                L=LIST_DOWN(LL)
                M_SLANT_DEPTH_NEAR(L)                                   &
     &            =M_SLANT_DEPTH_NEAR(L)/MU_V(L, ID)
                M_SLANT_DEPTH_INC(L)=TAU(L, I)/MU_V(L, ID)
                TAU_I(L)=TAU(L, I)
                FRAC_TAU_I(L)=1.0E+00_Real64
              ENDDO
            ELSE
!             The viewing level lies in the current layer.
              DO LL=1, N_LIST_DOWN
                L=LIST_DOWN(LL)
                TAU_I(L)=FRAC_RAD_LAYER(IR)*TAU(L, I)
                M_SLANT_DEPTH_NEAR(L)=0.0E+00_Real64
                M_SLANT_DEPTH_INC(L)=TAU_I(L)/MU_V(L, ID)
                FRAC_TAU_I(L)=FRAC_RAD_LAYER(IR)
              ENDDO
            ENDIF
!
!
            IF (ISOLIR == IP_SOLAR) THEN
!             Set the geometrical terms for the solar integral.
              DO LL=1, N_LIST_DOWN
                L=LIST_DOWN(LL)
!               This may exhibit ill-conditioning, so it is perturbed
!               using L'Hopital's rule (actually we keep two terms in
!               the expansion).
                D_MU=MU_0(L)+MU_V(L, ID)
                ETA=EPS_R/(D_MU+SIGN(SQ_EPS_R, D_MU))
                X=TAU_I(L)/(MU_0(L)*MU_V(L, ID))
                ETA_NM=1.0_Real64-ETA*X*(1.0_Real64+0.5_Real64*X*D_MU)
                GEOM_SOLAR(L)=(MU_0(L)/(D_MU+ETA))                      &
     &            *EXP(M_SLANT_DEPTH_NEAR(L))                           &
     &            *(EXP(-TAU_I(L)/MU_0(L))*ETA_NM                       &
     &            -EXP(M_SLANT_DEPTH_INC(L)))
              ENDDO
!             Add the contribution of the particular integral to the
!             radiance. In this case there is no factor representing
!             transmission from the top of the layer, since that is
!             intrinsically 1.
              DO LL=1, N_LIST_DOWN
                L=LIST_DOWN(LL)
                RADIANCE(L, IR, ID)                                     &
     &            =RADIANCE(L, IR, ID)+AZIM_FACTOR(L, ID)               &
     &            *LS_SUM_S(L)*OMEGA(L, I)*GEOM_SOLAR(L)
              ENDDO
            ELSE IF (ISOLIR == IP_INFRA_RED) THEN
              IF (MS == 0) THEN
                IF (L_IR_SOURCE_QUAD) THEN
                  print*, 'Not done'
                ELSE
                  DO LL=1, N_LIST_DOWN
                    L=LIST_DOWN(LL)
!                   The azimuthal factor is omitted since it will be 1.
                    IF (M_SLANT_DEPTH_INC(L) <  -SQ_EPS_R) THEN
                      X=EXP(M_SLANT_DEPTH_NEAR(L))                      &
     &                *(1.0E+00_Real64-EXP(M_SLANT_DEPTH_INC(L)))       &
     &                /M_SLANT_DEPTH_INC(L)
                    ELSE
!                     Keep the first couple of terms from a power
!                     series.
                      X=-EXP(M_SLANT_DEPTH_NEAR(L))*(1.0E+00_Real64     &
     &                  +0.5E+00_Real64*M_SLANT_DEPTH_INC(L))
                    ENDIF
                    RADIANCE(L, IR, ID)                                 &
     &                =RADIANCE(L, IR, ID)                              &
     &                +(DIFF_PLANCK(L, I)/PI)*X*FRAC_TAU_I(L)           &
     &                /(1.0_Real64-OMEGA(L, I)*PHASE_FNC(L, I, 1))
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
!
!           Determine the contribution from each eigenvalue.
            DO K=1, N_RED_EIGENSYSTEM
!             The term for u^+:
              DO LL=1, N_LIST_DOWN
                L=LIST_DOWN(LL)
                GEOM_INTEG_P(L)                                         &
     &            =(MU(L, K)/(MU(L, K)-MU_V(L, ID)))                    &
     &            *EXP(M_SLANT_DEPTH_NEAR(L))                           &
     &            *(EXP(-(TAU(L, I)-TAU_I(L))/MU(L, K))                 &
     &            -EXP(M_SLANT_DEPTH_INC(L)-TAU(L, I)/MU(L, K)))
!               The term for u^- may exhibit ill-conditioning,
!               so it is perturbed using L'Hopital's rule
!               (actually we keep two terms in the expansion).
                D_MU=MU(L, K)+MU_V(L, ID)
                ETA=EPS_R/(D_MU+SIGN(SQ_EPS_R, D_MU))
                X=TAU_I(L)/(MU(L, K)*MU_V(L, ID))
                ETA_NM=1.0_Real64-ETA*X*(1.0_Real64+0.5_Real64*X*D_MU)
                GEOM_INTEG_M(L)                                         &
     &            =(MU(L, K)/(D_MU+ETA))                                &
     &            *EXP(M_SLANT_DEPTH_NEAR(L))                           &
     &            *(EXP(-TAU_I(L)/MU(L, K))*ETA_NM                      &
     &            -EXP(M_SLANT_DEPTH_INC(L)))
              ENDDO
!
!             Combine to form the weights for each element of the
!             solution.
              DO LL=1, N_LIST_DOWN
                L=LIST_DOWN(LL)
                WEIGHT_U(L, IR, ID, K)                                  &
     &            =OMEGA(L, I)*LS_SUM_M(L, K)*GEOM_INTEG_M(L)
                WEIGHT_U(L, IR, ID, K+N_RED_EIGENSYSTEM)                &
     &            =OMEGA(L, I)*LS_SUM_P(L, K)*GEOM_INTEG_P(L)
              ENDDO
!
            ENDDO
!
          ENDIF
!
!         Add on the contribution from the conditioning homogeneous
!         solution. This includes some redundant calculation when
!         the weights will be zero for layers which cannot contribute.
!         Eventually, it may be better to tidy this up.
          IF (ISOLIR == IP_SOLAR) THEN
!
            DO L=1, N_PROFILE
              RADIANCE(L, IR, ID)=RADIANCE(L, IR, ID)                   &
     &          +AZIM_FACTOR(L, ID)                                     &
     &          *WEIGHT_U(L, IR, ID, K_SOL(L))*UPM_C(L, K_SOL(L))
            ENDDO
!
          ELSE IF (ISOLIR == IP_INFRA_RED) THEN
!
            DO K=1, N_RED_EIGENSYSTEM
              DO L=1, N_PROFILE
                RADIANCE(L, IR, ID)=RADIANCE(L, IR, ID)                 &
     &            +AZIM_FACTOR(L, ID)                                   &
     &            *(WEIGHT_U(L, IR, ID, K)*UPM_C(L, K)                  &
     &            +WEIGHT_U(L, IR, ID, K+N_RED_EIGENSYSTEM)             &
     &            *UPM_C(L, K+N_RED_EIGENSYSTEM))
              ENDDO
            ENDDO
!
          ENDIF
!
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE SET_DIRN_WEIGHTS
#endif
