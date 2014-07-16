#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set weights for the C.F. at the surface
!
! Purpose:
!   The contribution to the radiance of radiation reflected from the
!   surface is evaluated.
!
! Method:
!   The iterated expression for the radiance involves a contribution
!   from the radiance reflected from the surface. In principle, this
!   could be provided by the upward radiance, but in practice this
!   would be of low accuracy since the spherical harmonic series for
!   the radiance will be noisy at low orders of truncation. It is
!   better to evaluate the reflected radiance using the BRDFs, even
!   though it is more expensive to do so; this ensures that no
!   radiation will appear to be reflected from a non-reflecting
!   surface. Given these constraints the algorithm is essentially
!   straightforward.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_SURF_RAD(N_PROFILE, N_LAYER, TAU                  &
     &  , MS, LS_TRUNC, EULER_FACTOR                                    &
     &  , ISOLIR, I_DIRECT_SURF, MU_0, D_PLANCK_FLUX_SURFACE            &
     &  , N_BRDF_BASIS_FNC, LS_BRDF_TRUNC, F_BRDF                       &
     &  , RHO_ALB, BRDF_SOL, BRDF_HEMI, CGK                             &
     &  , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                  &
     &  , N_DIRECTION, MU_V, UP_LM, AZIM_FACTOR                         &
     &  , N_RED_EIGENSYSTEM, EIG_VEC, THETA, SOURCE_BASE                &
     &  , RADIANCE, WEIGHT_U                                            &
     &  , ND_PROFILE, ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL          &
     &  , ND_RED_EIGENSYSTEM, ND_MAX_ORDER, ND_BRDF_BASIS_FNC           &
     &  , ND_BRDF_TRUNC                                                 &
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
!     Dummy arguments.
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
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_BRDF_BASIS_FNC                                             &
!           Size allocated for basis functions of BRDFs
     &  , ND_BRDF_TRUNC
!           Size allocated for orders in basis functions of BRDFs
!
!
!     The atmosphere:
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of atmospheric profiles
     &  , N_LAYER
!           Number of atmospheric layers
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ND_LAYER)
!           Optical depths
!
!     Controlling spherical orders:
      INTEGER, INTENT(IN) ::                                            &
     &    MS                                                            &
!           Current azimuthal order
     &  , LS_TRUNC
!           Order of polar truncation
      REAL  (Real64), INTENT(IN) ::                                     &
     &    EULER_FACTOR
!           Factor applied to the last term of the series
!
!     Variables for solar or thermal sources
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Spectral region
      REAL  (Real64), INTENT(IN) ::                                     &
     &    I_DIRECT_SURF(ND_PROFILE)                                     &
!           The direct solar radiance at the surface
     &  , MU_0(ND_PROFILE)
!           Cosines of the solar zenith angle
      REAL  (Real64), INTENT(IN) ::                                     &
     &    D_PLANCK_FLUX_SURFACE(ND_PROFILE)
!           Differential Planckian flux at the surface
!
!     Variables related to the BRDFs
      INTEGER, INTENT(IN) ::                                            &
     &    N_BRDF_BASIS_FNC                                              &
!           Number of basis functions used in BRDFs
     &  , LS_BRDF_TRUNC
!           Order of polar truncation applied to BRDFs
      REAL  (Real64), INTENT(IN) ::                                     &
     &    F_BRDF(ND_BRDF_BASIS_FNC, 0: ND_BRDF_TRUNC/2                  &
     &      , 0: ND_BRDF_TRUNC/2, 0: ND_BRDF_TRUNC)                     &
!           BRDF basis functions
     &  , RHO_ALB(ND_PROFILE, ND_BRDF_BASIS_FNC)                        &
!           Weights of applied to the basis functions of the BRDF
     &  , BRDF_SOL(ND_PROFILE, ND_BRDF_BASIS_FNC, ND_DIRECTION)         &
!           The BRDF evaluated for scattering from the solar
!           beam into the viewing direction
     &  , BRDF_HEMI(ND_PROFILE, ND_BRDF_BASIS_FNC, ND_DIRECTION)
!           The BRDF evaluated for scattering from isotropic
!           radiation into the viewing direction
      REAL  (Real64), INTENT(IN) ::                                     &
     &    CGK(ND_BRDF_TRUNC/2+1, ND_MAX_ORDER)
!           Integrals of pairs of spherical harmonics over the downward
!           hemisphere

!
!     Viewing geometry:
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION                                                   &
!           Number of directions
     &  , N_VIEWING_LEVEL                                               &
!           Number of levels where the radiance is calculated
     &  , I_RAD_LAYER(ND_VIEWING_LEVEL)
!           Indices of layers containing viewing levels
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FRAC_RAD_LAYER(ND_VIEWING_LEVEL)                              &
!           Fraction optical depth into its layer of the
!           radiance level
     &  , MU_V(ND_PROFILE, ND_DIRECTION)                                &
!           Cosines of polar viewing angles
     &  , AZIM_FACTOR(ND_PROFILE, ND_DIRECTION)                         &
!           Azimuthal factors
     &  , UP_LM(ND_PROFILE, ND_MAX_ORDER+1, ND_DIRECTION)
!           Spherical harmonics at a fixed azimuthal order
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_RED_EIGENSYSTEM
!           Size of the reduced eigensystem
      REAL  (Real64), INTENT(IN) ::                                     &
     &    EIG_VEC(ND_PROFILE, 2*ND_RED_EIGENSYSTEM                      &
     &      , ND_RED_EIGENSYSTEM)                                       &
!           Eigenvalues of the full eigensystem scaled by
!           the s-parameters
     &  , THETA(ND_PROFILE, ND_RED_EIGENSYSTEM)                         &
!           Array of exponentials of optical depths along slant paths
     &  , SOURCE_BASE(ND_PROFILE, LS_TRUNC+1-MS)
!           Source function at the bottom of the layer
!
!
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    RADIANCE(ND_PROFILE, ND_VIEWING_LEVEL, ND_DIRECTION)
!           Radiances (to be incremented by the contribution of
!           the particular integral)
      REAL  (Real64), INTENT(INOUT) ::                                  &
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
!           Loop variable (viewing levels)
     &  , ID                                                            &
!           Loop variable (directions)
     &  , J                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , I                                                             &
!           Loop variable
     &  , LL
!           Loop variable
      INTEGER                                                           &
     &    LS                                                            &
!           Loop variable (polar orders)
     &  , LSR                                                           &
!           Loop variable (reduced polar orders)
     &  , LS_D                                                          &
!           Loop variable (polar orders)
     &  , LSR_D                                                         &
!           Loop variable (reduced polar orders)
     &  , LS_DD                                                         &
!           Loop variable (polar orders)
     &  , LSR_DD
!           Loop variable (reduced polar orders)
      INTEGER                                                           &
     &    N_LIST_UP                                                     &
!           Number of points where the viewing direction is upward
     &  , LIST_UP(ND_PROFILE)
!           List of points where the viewing direction is upward
      REAL  (Real64) ::                                                 &
     &    TRANS(ND_PROFILE)                                             &
!           Tranmission along the line of sight from the surface
!           to the viewing level
     &  , X
!           Temporary variable
!     Working arrays realated to the BRDF:
      REAL  (Real64) ::                                                 &
     &    XY(ND_PROFILE)                                                &
!           Product of (Clebsch-Gordan coefficient * kappa) and
!           the spherical harmonic at the current order in the
!           viewing direction
     &  , RYX(ND_PROFILE, LS_TRUNC-MS+1)                                &
!           Sum over basis functions of products of the above
!           and albedo weights
     &  , RVYX_M(ND_PROFILE, ND_RED_EIGENSYSTEM)                        &
!           Sum over polar orders of product of RYX and elements
!           of the eigenvalue for each eigenvalue for application
!           to terms in negative exponentials
     &  , RVYX_P(ND_PROFILE, ND_RED_EIGENSYSTEM)                        &
!           Sum over polar orders of product of RYX and elements
!           of the eigenvalue for each eigenvalue for application
!           to terms in positive exponentials
     &  , RSYX(ND_PROFILE)                                              &
!           Sum over polar orders of product of RYX and elements
!           of the source function
     &  , BRDF_FULL(ND_PROFILE)
!           Full BRDF weighted and summed over all basis functions
!
!
!
!     For each direction and observing level we calculate the
!     contribution of the particular integral to the radiance
!     from the surface and appropriate weightings for the
!     complementary function.

      DO ID=1, N_DIRECTION
!
!       Collect upward directions.
        N_LIST_UP=0
        DO L=1, N_PROFILE
          IF (MU_V(L, ID) >  0.0E+00_Real64) THEN
            N_LIST_UP=N_LIST_UP+1
            LIST_UP(N_LIST_UP)=L
          ENDIF
        ENDDO
!
!       Calculate the angular arrays related to the BRDF. At higher
!       azimuthal orders there will be contributions because all
!       terms of the BRDF that would contribute are beyond the
!       level of truncation and so are zero.
!
        IF (MS <= LS_BRDF_TRUNC-MOD(MS, 2)) THEN
!
          DO J=1, N_BRDF_BASIS_FNC
            DO LS=MS, LS_TRUNC
              LSR=LS-MS+1
              DO LS_D=MS, LS_BRDF_TRUNC-MOD(MS, 2), 2
                LSR_D=LS_D-MS+1
                X=0.0E+00_Real64
                DO LS_DD=MS, LS_BRDF_TRUNC-MOD(MS, 2), 2
                  LSR_DD=LS_DD-MS+1
                  X=X-CGK((LSR_DD+1)/2, LSR)*F_BRDF(J, LS_D, LS_DD, MS)
                ENDDO
                IF (LS_D == MS) THEN
!                 Initialize this time.
                  DO L=1, N_PROFILE
                    XY(L)=X*UP_LM(L, LSR_D, ID)
                  ENDDO
                ELSE
!                 Now add the increments.
                  DO L=1, N_PROFILE
                    XY(L)=XY(L)+X*UP_LM(L, LSR_D, ID)
                  ENDDO
                ENDIF
              ENDDO
              IF (J  == 1) THEN
!               Initialize this time.
                DO L=1, N_PROFILE
                  RYX(L, LSR)=RHO_ALB(L, 1)*XY(L)
                ENDDO
              ELSE
!               Increment for subsequent basis functions.
                DO L=1, N_PROFILE
                  RYX(L, LSR)=RYX(L, LSR)+RHO_ALB(L, J)*XY(L)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
!
          DO K=1, N_RED_EIGENSYSTEM
            DO L=1, N_PROFILE
              X=EULER_FACTOR*RYX(L, LS_TRUNC-MS+1)                      &
     &          *EIG_VEC(L, LS_TRUNC-MS+1, K)
              RVYX_M(L, K)=X
              RVYX_P(L, K)=X
            ENDDO
            DO LSR= LS_TRUNC-MS, 1, -1
              DO L=1, N_PROFILE
                X=RYX(L, LSR)*EIG_VEC(L, LSR, K)
                RVYX_M(L, K)=RVYX_M(L, K)                               &
     &            +X*REAL(1-2*MOD(LSR-1, 2), Real64)
                RVYX_P(L, K)=RVYX_P(L, K)+X
              ENDDO
            ENDDO
            DO L=1, N_PROFILE
              RVYX_M(L, K)=RVYX_M(L, K)*THETA(L, K)
            ENDDO
          ENDDO
!
          DO L=1, N_PROFILE
            RSYX(L)=EULER_FACTOR*RYX(L, LS_TRUNC-MS+1)                  &
     &        *SOURCE_BASE(L, LS_TRUNC-MS+1)
          ENDDO
          DO LSR= LS_TRUNC-MS, 1, -1
            DO L=1, N_PROFILE
              RSYX(L)=RSYX(L)+RYX(L, LSR)*SOURCE_BASE(L, LSR)
            ENDDO
          ENDDO
!
        ENDIF
!
!
        DO IR=1, N_VIEWING_LEVEL


!
!         Calculate minus the slantwise transmission from the
!         surface to the level in question. TRANS is used is
!         hold intermediate results.
          DO LL=1, N_LIST_UP
            L=LIST_UP(LL)
            TRANS(L)                                                    &
     &        =(1.0E+00_Real64                                          &
     &        -FRAC_RAD_LAYER(IR))*TAU(L, I_RAD_LAYER(IR))
          ENDDO
          DO I=I_RAD_LAYER(IR)+1, N_LAYER
            DO LL=1, N_LIST_UP
              L=LIST_UP(LL)
              TRANS(L)=TRANS(L)+TAU(L, I)
            ENDDO
          ENDDO
          DO LL=1, N_LIST_UP
            L=LIST_UP(LL)
            TRANS(L)                                                    &
     &        =EXP(-TRANS(L)/MU_V(L, ID))
          ENDDO
!
!         Add in the terms from the BRDF if in range.
          IF (MS <= LS_BRDF_TRUNC-MOD(MS,2)) THEN
!
            DO LL=1, N_LIST_UP
              L=LIST_UP(LL)
!             Add the contribution from the source function at the
!             base of the layer.
              RADIANCE(L, IR, ID)=RADIANCE(L, IR, ID)                   &
     &          +TRANS(L)*RSYX(L)*AZIM_FACTOR(L, ID)
            ENDDO
!
!           Increment the weights applied to the complementary function.
            DO K=1, N_RED_EIGENSYSTEM
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                WEIGHT_U(L, IR, ID, K)=WEIGHT_U(L, IR, ID, K)           &
     &            +TRANS(L)*RVYX_M(L, K)
                WEIGHT_U(L, IR, ID, K+N_RED_EIGENSYSTEM)                &
     &            =WEIGHT_U(L, IR, ID, K+N_RED_EIGENSYSTEM)             &
     &            +TRANS(L)*RVYX_P(L, K)
              ENDDO
            ENDDO
!
          ENDIF
!
!
!         Add the direct solar or thermal contributions to the radiance.
!         The azimuthal dependencies are included in the solar and
!         hemispheric parts of the BRDF, so the should be added in just
!         once, most naturally at the zeroth order.

          IF (MS == 0) THEN
            IF (ISOLIR == IP_SOLAR) THEN
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                BRDF_FULL(L)=RHO_ALB(L, 1)*BRDF_SOL(L, 1, ID)
              ENDDO

              DO J=2, N_BRDF_BASIS_FNC
                DO LL=1, N_LIST_UP
                  L=LIST_UP(LL)
                  BRDF_FULL(L)=BRDF_FULL(L)                             &
     &              +RHO_ALB(L, J)*BRDF_SOL(L, J, ID)
                ENDDO


              ENDDO
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                RADIANCE(L, IR, ID)=RADIANCE(L, IR, ID)                 &
     &            +TRANS(L)*I_DIRECT_SURF(L)*MU_0(L)                    &
     &            *BRDF_FULL(L)
              ENDDO

            ELSE IF (ISOLIR == IP_INFRA_RED) THEN
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                BRDF_FULL(L)=RHO_ALB(L, 1)*BRDF_HEMI(L, 1, ID)
              ENDDO
              DO J=2, N_BRDF_BASIS_FNC
                DO LL=1, N_LIST_UP
                  L=LIST_UP(LL)
                  BRDF_FULL(L)=BRDF_FULL(L)                             &
     &              +RHO_ALB(L, J)*BRDF_HEMI(L, J, ID)
                ENDDO
              ENDDO
              DO LL=1, N_LIST_UP
                L=LIST_UP(LL)
                RADIANCE(L, IR, ID)=RADIANCE(L, IR, ID)                 &
     &            +TRANS(L)                                             &
     &            *(1.0E+00_Real64-BRDF_FULL(L))                        &
     &            *D_PLANCK_FLUX_SURFACE(L)/PI
              ENDDO
            ENDIF
          ENDIF
!
        ENDDO
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE CALC_SURF_RAD
#endif
