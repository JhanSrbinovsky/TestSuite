#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve for radiances in harmonics.
!
! Method:
!       After setting the basic properties for the radiance solver
!       a matrix is built and solved.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SPH_SOLVER(IERR                                        &
!                     Atmospheric sizes
     &  , N_PROFILE, N_LAYER                                            &
!                     Angular integration
     &  , MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC                  &
     &  , CG_COEFF, UPLM_ZERO, IA_SPH_MM                                &
     &  , ACCURACY_ADAPTIVE, EULER_FACTOR                               &
     &  , I_SPH_ALGORITHM, I_SPH_MODE, L_RESCALE                        &
!                     Spectral Region
     &  , ISOLIR                                                        &
!                       Options for Equivalent Extinction
     &  , L_SCALE_SOLAR, ADJUST_SOLAR_KE                                &
!                     Solar Fields
     &  , I_DIRECT, MU_0, UPLM_SOL                                      &
!                     Infra-red Properties
     &  , DIFF_PLANCK, FLUX_INC_DOWN                                    &
     &  , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                               &
!                       Optical properies
     &  , TAU, OMEGA, PHASE_FNC, PHASE_FNC_SOLAR                        &
     &  , FORWARD_SCATTER                                               &
!                     Surface Conditions
     &  , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB                      &
     &  , F_BRDF, BRDF_SOL, BRDF_HEMI                                   &
     &  , D_PLANCK_FLUX_SURFACE                                         &
!                     Levels for calculating radiances
     &  , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                  &
!                       Viewing Geometry
     &  , N_DIRECTION, DIRECTION                                        &
!                     Calculated Radiances or Fluxes
     &  , FLUX_DIRECT, FLUX_TOTAL, RADIANCE_MONO, PHOTOLYSIS            &
!                     Dimensions of arrays
     &  , ND_PROFILE, ND_LAYER                                          &
     &  , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE            &
     &  , ND_MAX_ORDER, ND_SPH_COEFF                                    &
     &  , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC                              &
     &  , ND_RED_EIGENSYSTEM, ND_SPH_EQUATION, ND_SPH_DIAGONAL          &
     &  , ND_SPH_CF_WEIGHT, ND_SPH_U_RANGE                              &
     &  , ND_VIEWING_LEVEL, ND_DIRECTION                                &
     &  )
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ND_FLUX_PROFILE                                               &
!           Size allocated for profiles where fluxes are calculated
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for profiles where radiances are calculated
     &  , ND_J_PROFILE                                                  &
!           Size allocated for profiles where mean radiances
!           are calculated
     &  , ND_LAYER                                                      &
!           Maximum number of layers
     &  , ND_DIRECTION                                                  &
!           Size allowed for viewing directions
     &  , ND_VIEWING_LEVEL                                              &
!           Size allowed for levels where the radiance is calculated
     &  , ND_MAX_ORDER                                                  &
!           Size allowed for orders of spherical harmonics
     &  , ND_BRDF_BASIS_FNC                                             &
!           Size allowed for BRDF basis functions
     &  , ND_BRDF_TRUNC                                                 &
!           Size allowed for orders of BRDFs
     &  , ND_SPH_COEFF                                                  &
!           Size allowed for spherical harmonic coefficients
     &  , ND_RED_EIGENSYSTEM                                            &
!           Size allowed for the reduced eigensystem
     &  , ND_SPH_EQUATION                                               &
!           Size allowed for spherical harmonic equations
     &  , ND_SPH_DIAGONAL                                               &
!           Size allowed for diagonals of the spherical harmonic
!           matrix
     &  , ND_SPH_CF_WEIGHT                                              &
!           Size allowed for application of weights of the C. F.
     &  , ND_SPH_U_RANGE
!           Size allowed for the range of u^+|- contributing on any
!           viewing level
!
!     Include header files.
#include "c_kinds.h"
#include "c_pi.h"
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
#include "sph_mode_pcf3z.h"
#include "sph_algorithm_pcf3z.h"
#include "sph_truncation_pcf3z.h"
#include "spectral_region_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
!
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , ISOLIR
!           Spectral region
!                     Angular integration
      LOGICAL, INTENT(IN) ::                                            &
     &    L_RESCALE
!           Flag for rescaling of the optical properties
      INTEGER, INTENT(IN) ::                                            &
     &    MS_MIN                                                        &
!           Lowest azimuthal order calculated
     &  , MS_MAX                                                        &
!           Highest azimuthal order calculated
     &  , IA_SPH_MM(0: ND_MAX_ORDER)                                    &
!           Address of spherical coefficient for (m, m) for each m
     &  , LS_LOCAL_TRUNC(0: ND_MAX_ORDER)                               &
!           Orders of truncation at each azimuthal order
     &  , I_TRUNCATION                                                  &
!           Type of truncation used
     &  , I_SPH_MODE                                                    &
!           Mode in which the spherical solver is used
     &  , I_SPH_ALGORITHM
!           Algorithm used to solve the spherical system
      REAL  (Real64), INTENT(IN) ::                                     &
     &    CG_COEFF(ND_SPH_COEFF)                                        &
!           Clebsch-Gordan coefficients
     &  , UPLM_ZERO(ND_SPH_COEFF)                                       &
!           Values of spherical harmonics at polar angles pi/2
     &  , UPLM_SOL(ND_PROFILE, ND_SPH_COEFF)                            &
!           Values of spherical harmonics in the solar direction
     &  , ACCURACY_ADAPTIVE                                             &
!           Accuracy for adaptive truncation
     &  , EULER_FACTOR
!           Factor applied to the last polar order
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Use quadratic source term
!
!                       Variables for equivalent extinction
      LOGICAL, INTENT(IN) ::                                            &
     &    L_SCALE_SOLAR
!           Apply scaling to solar flux
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ADJUST_SOLAR_KE(ND_PROFILE, ND_LAYER)
!           Adjustment of solar beam with equivalent extinction
!
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ND_LAYER)                                     &
!           Optical depth
     &  , OMEGA(ND_PROFILE, ND_LAYER)                                   &
!           Albedos of single scattering
     &  , PHASE_FNC(ND_PROFILE, ND_LAYER, ND_MAX_ORDER)                 &
!           Moments of the phase function
     &  , PHASE_FNC_SOLAR(ND_RADIANCE_PROFILE, ND_LAYER, ND_DIRECTION)  &
!           The phase function evaluated for scattering from
!           the solar beam into the viewing directions
     &  , FORWARD_SCATTER(ND_PROFILE, ND_LAYER)
!           Forward scattering fractions

      REAL (Real64), INTENT(INOUT) ::                                   &
     &    MU_0(ND_PROFILE)
!           Cosines of solar zenith angles

      REAL (Real64), INTENT(IN) ::                                      &
     &    DIFF_PLANCK(ND_PROFILE, ND_LAYER)                             &
!           Difference in the FLUX Planckian function
     &  , DIFF_PLANCK_2(ND_PROFILE, ND_LAYER)                           &
!           2x2nd differences of Planckian
     &  , FLUX_INC_DOWN(ND_PROFILE)
!           Incident downward flux (in real calculations this is used
!           only in the IR where it is assumed to be Planckian, but
!           in may be used in the solar in idealized test cases)
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    I_DIRECT(ND_PROFILE, 0: ND_LAYER)
!           Direct solar radiance (only the first row is set on input)
      INTEGER, INTENT(IN) ::                                            &
     &    LS_BRDF_TRUNC                                                 &
!           Order of trunation of BRDFs
     &  , N_BRDF_BASIS_FNC
!           Number of BRDF basis functions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    D_PLANCK_FLUX_SURFACE(ND_PROFILE)                             &
!           Differential Planckian flux from the surface
     &  , RHO_ALB(ND_PROFILE, ND_BRDF_BASIS_FNC)                        &
!           Weights of the basis functions
     &  , F_BRDF(ND_BRDF_BASIS_FNC, 0: ND_BRDF_TRUNC/2                  &
     &      , 0: ND_BRDF_TRUNC/2, 0: ND_BRDF_TRUNC)                     &
!           Array of BRDF basis terms
     &  , BRDF_SOL(ND_PROFILE, ND_BRDF_BASIS_FNC, ND_DIRECTION)         &
!           The BRDF evaluated for scattering from the solar
!           beam into the viewing direction
     &  , BRDF_HEMI(ND_PROFILE, ND_BRDF_BASIS_FNC, ND_DIRECTION)
!           The BRDF evaluated for scattering from isotropic
!           radiation into the viewing direction
!
!                       Viewing Geometry
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION                                                   &
!           Number of viewing directions
     &  , N_VIEWING_LEVEL                                               &
!           Number of levels where radiances are calculated
     &  , I_RAD_LAYER(ND_VIEWING_LEVEL)
!           Layers in which radiances are calculated
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DIRECTION(ND_RADIANCE_PROFILE, ND_DIRECTION, 2)               &
!           Viewing directions
     &  , FRAC_RAD_LAYER(ND_VIEWING_LEVEL)
!           Fractions below the tops of the layers
!
!                       Radiances calculated
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_FLUX_PROFILE, 0: ND_LAYER)                     &
!           Direct Flux
     &  , FLUX_TOTAL(ND_FLUX_PROFILE, 2*ND_LAYER+2)                     &
!           Total Fluxes
     &  , RADIANCE_MONO(ND_RADIANCE_PROFILE                             &
     &      , ND_VIEWING_LEVEL, ND_DIRECTION)                           &
!           Radiances
     &  , PHOTOLYSIS(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Rates of photolysis
!
!
!     Local variabales.
      INTEGER                                                           &
     &    MS                                                            &
!           Azimuthal order
     &  , LSR                                                           &
!           Reduced polar order
     &  , N_RED_EIGENSYSTEM                                             &
!           Size of the reduced eigensystem
     &  , N_EQUATION                                                    &
!           Number of equations
     &  , I                                                             &
!           Loop variable
     &  , ID                                                            &
!           Loop variable (directions)
     &  , L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    KAPPA(ND_MAX_ORDER/2, ND_MAX_ORDER/2)                         &
!           Integrals of pairs of spherical harmonics over the downward
!           hemisphere
     &  , CGK(ND_BRDF_TRUNC/2+1, ND_MAX_ORDER)                          &
!           Products of the Clebsch-Gordan coefficients and the
!           hemispheric integrals
     &  , UP_LM(ND_PROFILE, ND_MAX_ORDER+1, ND_DIRECTION)               &
!           Polar parts of spherical harmonics
     &  , WEIGHT_U(ND_PROFILE, ND_VIEWING_LEVEL                         &
     &      , ND_SPH_CF_WEIGHT, ND_SPH_U_RANGE)                         &
!           Weights to be applied to the vector U containing the
!           complementary functions
     &  , C_YLM(ND_PROFILE, ND_VIEWING_LEVEL, ND_MAX_ORDER+1)           &
!           Coefficients in the expansion of the radiance
!           in spherical harmonics
     &  , A(ND_PROFILE, ND_SPH_EQUATION, ND_SPH_DIAGONAL)               &
!           Matrix on the LHS of the equation for spherical
!           harmonics
     &  , B(ND_PROFILE, ND_SPH_EQUATION)                                &
!           RHS of matrix equation
     &  , UPM(ND_PROFILE, ND_SPH_EQUATION)                              &
!           Variables u+|-
     &  , AZIM_FACTOR(ND_PROFILE, ND_DIRECTION)
!           Azimuthal factors
      INTEGER                                                           &
     &    LS_TRUNC_CALC                                                 &
!           Order of truncation required in calculations
     &  , LS_BRDF_TRUNC_CALC                                            &
!           Order of truncation of BRDFs required in calculations
     &  , LS_SIGNIFICANT
!           Maximum significant polar order
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    HEMI_SPH_INTEG, CG_KAPPA_MS                                   &
     &  , BUILD_SPH_MATRIX, SPH_MATRIX_SOLVER
!
!
!
!     Calculate the direct radiances which are independent of the
!     azimuthal order.
      IF (ISOLIR == IP_SOLAR) THEN
!
        IF (L_SCALE_SOLAR) THEN
          DO I=1, N_LAYER
            DO L=1, N_PROFILE
              I_DIRECT(L, I)=I_DIRECT(L, I-1)                           &
     &          *ADJUST_SOLAR_KE(L, I)*EXP(-TAU(L, I)/MU_0(L))
            ENDDO
          ENDDO
        ELSE
          DO I=1, N_LAYER
            DO L=1, N_PROFILE
              I_DIRECT(L, I)=I_DIRECT(L, I-1)*EXP(-TAU(L, I)/MU_0(L))
            ENDDO
          ENDDO
        ENDIF
!
      ENDIF
!
!     Initialize the monochromatic radiances if they are required.
      IF (I_SPH_MODE == IP_SPH_MODE_RAD) THEN
        DO ID=1, N_DIRECTION
          DO I=1, N_VIEWING_LEVEL
            DO L=1, N_PROFILE
              RADIANCE_MONO(L, I, ID)=0.0E+00_Real64
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
!
!
!     For each non-negative azimuthal order a matrix of coefficients
!     is built and solved.
      DO MS=MS_MIN, MS_MAX
!
!       Set the order of truncation required for the calculation.
!       If an adaptive truncation is used we may not calculate all
!       the declared polar orders, though space is reserved for them
!       and is initialized (this is necessary as a subsequent call to
!       this routine with different optical properties may require a
!       different number of terms). On entry we have no estimate of the
!       number of terms required, other than the assigned truncation.
        IF (I_TRUNCATION == IP_TRUNC_ADAPTIVE) THEN
          IF (MS == MS_MIN) THEN
            LS_TRUNC_CALC=LS_LOCAL_TRUNC(MS)
          ELSE
!           The polar order must exceed the azimuthal order by at least
!           one, and must give an even number of harmonics.
            LS_TRUNC_CALC=MAX((MS+1), LS_SIGNIFICANT)
            LS_TRUNC_CALC=LS_TRUNC_CALC+MOD((LS_TRUNC_CALC+MS+1), 2)
            LS_TRUNC_CALC=MIN(LS_TRUNC_CALC, LS_LOCAL_TRUNC(MS))
          ENDIF
        ELSE
          LS_TRUNC_CALC=LS_LOCAL_TRUNC(MS)
        ENDIF
!
!       Azimuthal factors are required if calculating radiances.
        IF (I_SPH_MODE == IP_SPH_MODE_RAD) THEN
!         The recalculation of azimuthal factors for each monochromatic
!         radiance is inefficient, but avoids the need for the storage
!         of arrays with different azimuthal orders. The equation of
!         transfer is solved only for positive azimuthal orders. Hence,
!         if the azimuthal order is 0 we add in the term, but if the
!         azimuthal order is non-zero it also represents the negative
!         order and acquires a weighting factor of 2.
          IF (MS == 0) THEN
            DO ID=1, N_DIRECTION
              DO L=1, N_PROFILE
                AZIM_FACTOR(L, ID)=1.0E+00_Real64
              ENDDO
            ENDDO
          ELSE
            DO ID=1, N_DIRECTION
              DO L=1, N_PROFILE
                AZIM_FACTOR(L, ID)                                      &
     &          =2.0E+00_Real64*COS(REAL(MS,Real64)*DIRECTION(L,ID,2))
              ENDDO
            ENDDO
          ENDIF
!         Calculate spherical harmonics in the viewing directions.
!         A judgement about storage has been made here. It is
!         deemed acceptable to store values for a fixed azimuthal
!         order, but deemed that too much storage would be required
!         for all azmuthal orders to be held at once.
          DO ID=1, N_DIRECTION
! DEPENDS ON: eval_uplm
            CALL EVAL_UPLM(MS, LS_TRUNC_CALC                            &
     &        , N_PROFILE, DIRECTION(1, ID, 1), UP_LM(1, 1, ID)         &
     &        , ND_PROFILE)
          ENDDO
        ENDIF
!
!       Calculate integrals of products of spherical harmonics for
!       use in Marshak's boundary conditions. These arrays are
!       recalculated each time to save storage.
! DEPENDS ON: hemi_sph_integ
        CALL HEMI_SPH_INTEG(LS_TRUNC_CALC, MS, UPLM_ZERO(IA_SPH_MM(MS)) &
     &    , KAPPA                                                       &
     &    , ND_MAX_ORDER                                                &
     &    )
!       Perform preliminary calculations for the BRDF: we need to check
!       that the order of truncation of BRDFs does not exceed the
!       order of calculation: it must also be even.
        LS_BRDF_TRUNC_CALC                                              &
     &    =MIN(LS_BRDF_TRUNC, LS_TRUNC_CALC-MOD(LS_TRUNC_CALC, 2))
! DEPENDS ON: cg_kappa_ms
        CALL CG_KAPPA_MS(MS, LS_TRUNC_CALC, LS_BRDF_TRUNC_CALC          &
     &    , CG_COEFF, KAPPA                                             &
     &    , CGK                                                         &
     &    , ND_MAX_ORDER, ND_BRDF_TRUNC                                 &
     &    )
!
!       Initialize the spherical harmonics at this azimuthal order.
        IF ( (I_SPH_MODE == IP_SPH_MODE_FLUX).OR.                       &
     &       (I_SPH_ALGORITHM == IP_SPH_DIRECT) ) THEN
          DO I=1, N_VIEWING_LEVEL
            DO LSR=1, LS_LOCAL_TRUNC(MS)+1-MS
              DO L=1, N_PROFILE
                C_YLM(L, I, LSR)=0.0E+00_Real64
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
!       In the infra-red region differential quantities are used. The
!       incident flux will then be determined from the Plankian flux
!       (which is used for consistency with two-stream calculations)
!       by dividing by pi: we must also multiply by sqrt(4.pi) to get
!       the weighting of the zeroth spherical harmonic. Other elements
!       of the array were zeroed before.
!
        N_RED_EIGENSYSTEM=(LS_TRUNC_CALC+1-MS)/2
!
! DEPENDS ON: build_sph_matrix
        CALL BUILD_SPH_MATRIX(I_SPH_ALGORITHM, EULER_FACTOR             &
!                       Basic sizes
     &    , N_PROFILE, N_LAYER, LS_TRUNC_CALC                           &
     &    , MS, N_RED_EIGENSYSTEM                                       &
!                       Numerical arrays of spherical terms
     &    , CG_COEFF(IA_SPH_MM(MS)), KAPPA, UP_LM                       &
!                       Solar variables
     &    , ISOLIR, I_DIRECT, MU_0, UPLM_SOL(1, IA_SPH_MM(MS))          &
     &    , AZIM_FACTOR                                                 &
!                       Infra-red variables
     &    , DIFF_PLANCK, L_IR_SOURCE_QUAD, DIFF_PLANCK_2                &
!                       Isotropic incident flux
     &    , FLUX_INC_DOWN                                               &
!                       Optical properties
     &    , TAU, OMEGA, PHASE_FNC                                       &
!                       Surface Fields
     &    , LS_BRDF_TRUNC_CALC, N_BRDF_BASIS_FNC, RHO_ALB               &
     &    , F_BRDF, BRDF_SOL, BRDF_HEMI, CGK                            &
     &    , D_PLANCK_FLUX_SURFACE                                       &
!                       Levels where radiances are calculated
     &    , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                &
!                       Viewing Geometry
     &    , N_DIRECTION, DIRECTION(1, 1, 1)                             &
!                       Output variables
     &    , A, B, C_YLM, WEIGHT_U, RADIANCE_MONO                        &
!                       Dimensions
     &    , ND_PROFILE, ND_RADIANCE_PROFILE                             &
     &    , ND_LAYER, ND_VIEWING_LEVEL, ND_DIRECTION                    &
     &    , ND_MAX_ORDER, ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC              &
     &    , ND_RED_EIGENSYSTEM, ND_SPH_EQUATION, ND_SPH_DIAGONAL        &
     &    , ND_SPH_CF_WEIGHT, ND_SPH_U_RANGE                            &
     &    )
        IF (IERR /= I_NORMAL) RETURN
!
!       Apply standard Gaussian elimination to obtain coefficients
!       u^+|-.
!
        N_EQUATION=2*N_LAYER*N_RED_EIGENSYSTEM
!
! DEPENDS ON: sph_matrix_solver
        CALL SPH_MATRIX_SOLVER(N_PROFILE, N_LAYER, N_RED_EIGENSYSTEM    &
     &    , A, B                                                        &
     &    , UPM                                                         &
     &    , ND_PROFILE, ND_SPH_EQUATION, ND_SPH_DIAGONAL                &
     &    )
!
!       Increment the radiances with the contributions from
!       the complementary function.
! DEPENDS ON: increment_rad_cf
        CALL INCREMENT_RAD_CF(N_PROFILE                                 &
     &    , N_DIRECTION, AZIM_FACTOR                                    &
     &    , N_VIEWING_LEVEL, I_RAD_LAYER                                &
     &    , I_SPH_MODE, I_SPH_ALGORITHM                                 &
     &    , MS, LS_TRUNC_CALC, EULER_FACTOR                             &
     &    , ISOLIR, MU_0, KAPPA, UP_LM                                  &
     &    , N_RED_EIGENSYSTEM, N_EQUATION, WEIGHT_U, UPM                &
     &    , I_DIRECT, C_YLM, FLUX_DIRECT, FLUX_TOTAL                    &
     &    , RADIANCE_MONO, PHOTOLYSIS                                   &
     &    , ND_PROFILE, ND_FLUX_PROFILE                                 &
     &    , ND_RADIANCE_PROFILE, ND_J_PROFILE                           &
     &    , ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL                    &
     &    , ND_MAX_ORDER, ND_SPH_EQUATION, ND_SPH_CF_WEIGHT             &
     &    , ND_SPH_U_RANGE                                              &
     &    )
!
        IF (I_TRUNCATION == IP_TRUNC_ADAPTIVE) THEN
!
!         Reduce the polar order of truncation if higher polar orders
!         make an insignificant contribution to the radiance field.
!         At least two polar orders are required.
          LS_SIGNIFICANT=MS+1
          DO LSR=MS+2, LS_TRUNC_CALC+1-MS
            DO I=1, N_VIEWING_LEVEL
              DO L=1, N_PROFILE
                IF (ABS(C_YLM(L, I, LSR)) >  ACCURACY_ADAPTIVE          &
     &            *ABS(C_YLM(L, I, 1))) LS_SIGNIFICANT=LSR+MS-1
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
      ENDDO
!
      IF ( (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER).AND.                &
     &     (ISOLIR == IP_SOLAR) ) THEN
!       Add in the singly scattered solar beam using the
!       potentially higher order of truncation.
! DEPENDS ON: single_scat_sol
        CALL SINGLE_SCAT_SOL(N_PROFILE, N_LAYER                         &
     &    , N_DIRECTION, DIRECTION                                      &
     &    , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                &
     &    , I_DIRECT, MU_0                                              &
     &    , TAU, OMEGA, PHASE_FNC_SOLAR                                 &
     &    , RADIANCE_MONO                                               &
     &    , ND_PROFILE, ND_RADIANCE_PROFILE                             &
     &    , ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL                    &
     &    )
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE SPH_SOLVER
#endif
