#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to build up the matrix for radiances.
!
! Purpose:
!   This routine assembles the stepped matrix to solve the equation
!   of transfer for a specific value of the azimuthal quantum
!   number.
!
! Method:
!
!   Labelling of variables and equations:
!     The variables are u_{ik}^{-|+} where i runs over all the layers
!   1,...,N_LAYER and k runs over positive eigenvalues of the reduced
!   eigensystem 1,...,N_RED_EIGENSYSTEM: there are thus 2n_e variables
!   describing the radiance in each layer, so the number of a variable
!   is
!      IV=2n_e(i-1)+k+n_e(1+|-1)/2
!   (Note that u_{ik}^- preceeds u_{ik}^+ by N_RED_EIGENSYSTEM).
!   At the top of the atmosphere (L'+1-m)/2 conditions are applied by
!   Marshak's conditions, where l'=m+1,...,L' in steps of 2, so for
!   this boundary
!      IE=(l'+1-m)/2
!   At the i'th interior boundary a condition of continuity is applied
!   to I_{lm}, where l=m,...,L'. To match the numbering of the equations
!   at the boundary values of l=m, m+2,...,L'-1 in steps of 2 are taken
!   first followed by those with l=m+1,...,L', so the number of the
!   equation is
!      IE=n_e(2i-1)+(l-m)/2+1,           l=m, m+2,...,L'-1
!      IE=n_e(2i-1)+(l+1-m)/2+n_e,       l=m+1,...,L',
!   allowing for n_e conditions at the top of the model and 2n_e
!   conditions at higher interfaces. At the bottom of the atmosphere
!   Marshak's condition is imposed using the harmonics l'=m+1,...,L'
!   in steps of 2, so the numbering of equations is
!      IE=n_e(2N_LAYER-1)+(l'+1-m)/2
!     Each of these equations couples together u_{ik}^{+|-} in the
!   layers above and below the interface. Hence, each equation
!   IE=(2i-1)n_e+1,...,(2i+1)n_e involves the variables IV=(2i-1)n_e+1,
!   ...,(2i+1)n_e, producing a stepped diagonal matrix which can be
!   encoded in an array of 4n_e columns with IE indexing the rows.
!   3n_e-1 sub-diagonals. The mapping is:
!     (IE, IV) --> (IE, IV-2*N_RED_EIGENSYSTEM*(I-1))
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE BUILD_SPH_MATRIX(I_SPH_ALGORITHM, EULER_FACTOR         &
!                     Basic sizes
     &  , N_PROFILE, N_LAYER, LS_TRUNC, MS, N_RED_EIGENSYSTEM           &
!                     Numerical arrays of spherical terms
     &  , CG_COEFF, KAPPA, UP_LM                                        &
!                     Solar variables
     &  , ISOLIR, I_DIRECT, MU_0, UPLM_SOL, AZIM_FACTOR                 &
!                     Infra-red variables
     &  , DIFF_PLANCK, L_IR_SOURCE_QUAD, DIFF_PLANCK_2                  &
!                     Diffuse incident field
     &  , FLUX_DOWN_INC                                                 &
!                     Optical properies
     &  , TAU, OMEGA, PHASE_FNC                                         &
!                     Surface Fields
     &  , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB                      &
     &  , F_BRDF, BRDF_SOL, BRDF_HEMI, CGK                              &
     &  , D_PLANCK_FLUX_SURFACE                                         &
!                     Levels where radiances are calculated
     &  , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                  &
!                       Viewing Geometry
     &  , N_DIRECTION, MU_V                                             &
!                     Output variables
     &  , A, B, C_YLM, WEIGHT_U, RADIANCE                               &
!                     Dimensions
     &  , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER                     &
     &  , ND_VIEWING_LEVEL, ND_DIRECTION                                &
     &  , ND_MAX_ORDER, ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC                &
     &  , ND_RED_EIGENSYSTEM, ND_SPH_EQUATION, ND_SPH_DIAGONAL          &
     &  , ND_SPH_CF_WEIGHT, ND_SPH_U_RANGE                              &
     &  )
!
!
!
!
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
!
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for atmospheric profiles where radiances
!           are calculated
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_VIEWING_LEVEL                                              &
!           Allocated size for levels where radiances are calculated
     &  , ND_DIRECTION                                                  &
!           Allocated size for viewing directions
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_BRDF_BASIS_FNC                                             &
!           Size allocated for BRDF basis functions
     &  , ND_BRDF_TRUNC                                                 &
!           Size allocated for orders in BRDFs
     &  , ND_RED_EIGENSYSTEM                                            &
!           Size allocated for the reduced eigensystem
     &  , ND_SPH_EQUATION                                               &
!           Size allocated for spherical harmonic equations
     &  , ND_SPH_DIAGONAL                                               &
!           Size allocated for diagonals in matrix for harmonics
     &  , ND_SPH_CF_WEIGHT                                              &
!           Size allocated for enetities to be incremented by the
!           complementary function
     &  , ND_SPH_U_RANGE
!           Range of values of u^+|- contributing on any viewing
!           level
!
!     Include header files.
#include "c_pi.h"
#include "sph_algorithm_pcf3z.h"
#include "spectral_region_pcf3z.h"
!
!     Dummy arguments
!     Atmospheric structrure:
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of atmospheric layers
     &  , N_LAYER
!           Number of atmospheric layers
!
!     Spherical harmonic structure:
      INTEGER, INTENT(IN) ::                                            &
     &    I_SPH_ALGORITHM                                               &
!           Algorithm for the spherical harmonic solution
     &  , LS_TRUNC                                                      &
!           The truncating order of the system of equations
     &  , MS                                                            &
!           Azimuthal order
     &  , N_RED_EIGENSYSTEM
!           Size of the reduced eigensystem
      REAL  (Real64), INTENT(IN) ::                                     &
     &    EULER_FACTOR
!           Factor applied to the last term of an alternating series
!
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Flag for spectral region
!
!     Optical properties:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ND_LAYER)                                     &
!           Optical depths of the layers
     &  , OMEGA(ND_PROFILE, ND_LAYER)                                   &
!           Albedos of single scattering of the layers
     &  , PHASE_FNC(ND_PROFILE, ND_LAYER, ND_MAX_ORDER)
!           Phase functions of the layers
!
!     Solar Fields:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    MU_0(ND_PROFILE)                                              &
!           Cosine of solar zenith angle
     &  , I_DIRECT(ND_PROFILE, 0: ND_LAYER)                             &
!           The direct solar radiance
     &  , UPLM_SOL(ND_PROFILE, LS_TRUNC+2-MS)
!           Spherical harmonics of the solar angle
!
!     Infra-red quantities:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Flag for quadratic source function in the IR
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DIFF_PLANCK(ND_PROFILE, ND_LAYER)                             &
!           Differences in the hemispheric Planckian FLUX (bottom-top)
!           across the layer
     &  , DIFF_PLANCK_2(ND_PROFILE, ND_LAYER)
!           Twice the second differences in the hemispheric Planckian
!           FLUX
      REAL  (Real64), INTENT(IN) ::                                     &
     &    CG_COEFF(LS_TRUNC+1-MS)                                       &
!           Clebsch-Gordan coefficients
     &  , KAPPA(ND_MAX_ORDER/2, ND_MAX_ORDER/2)                         &
!           Integrals of pairs of spherical harmonics over the downward
!           hemisphere
     &  , CGK(ND_BRDF_TRUNC/2+1, ND_MAX_ORDER)                          &
!           Products of the Clebsch-Gordan coefficients and the
!           hemispheric integrals
     &  , UP_LM(ND_PROFILE, ND_MAX_ORDER+1, ND_DIRECTION)               &
!           Polar parts of spherical harmonics in viewing directions
     &  , FLUX_DOWN_INC(ND_PROFILE)
!           Diffuse hemispherically isotropic incident flux
!
!     Surface Fields:
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
     &    MU_V(ND_PROFILE, ND_DIRECTION)                                &
!           Cosines of polar viewing directions
     &  , AZIM_FACTOR(ND_PROFILE, ND_DIRECTION)                         &
!           Azimuthal factors
     &  , FRAC_RAD_LAYER(ND_VIEWING_LEVEL)
!           Fractions below the tops of the layers
!
!
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    RADIANCE(ND_RADIANCE_PROFILE                                  &
     &      , ND_VIEWING_LEVEL, ND_DIRECTION)
!           Radiances to be incremented (note that at lower
!           levels this is declared withh ND_PROFILE, but this
!           is fine since those routines will be called only
!           when the two sizes are equal)
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    A(ND_PROFILE, ND_SPH_EQUATION, ND_SPH_DIAGONAL)               &
!           Matrix in the LHS of the equations
     &  , B(ND_PROFILE, ND_SPH_EQUATION)                                &
!           Vector of forcings for the matrix equation
     &  , WEIGHT_U(ND_PROFILE, ND_VIEWING_LEVEL, ND_SPH_CF_WEIGHT       &
     &      , ND_SPH_U_RANGE)
!           Weights to be applied to the vector U containing the
!           complementary functions
!
!                     Radiances
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    C_YLM(ND_PROFILE, ND_VIEWING_LEVEL, LS_TRUNC+1-MS)
!           Coefficients for radiances
!
!
!     Local variables
      INTEGER                                                           &
     &    IE                                                            &
!           Number of the equation
     &  , IVMA                                                          &
!           Index for u^- in the layer above the interface
!           This and the next three variables are used in two forms:
!           with an offset in indexing the matrix A and without
!           an offset in indexing the array WEIGHT_U.
     &  , IVMB                                                          &
!           Index for u^- in the layer below the interface
     &  , IVPA                                                          &
!           Index for u^+ in the layer above the interface
     &  , IVPB                                                          &
!           Index for u^+ in the layer below the interface
     &  , I_ABOVE                                                       &
!           Index for layer above in two-dimensional arrays
     &  , I_BELOW                                                       &
!           Index for layer below in two-dimensional arrays
     &  , I_ASSIGN_LEVEL
!           Level where a radiance is to be assigned
      INTEGER                                                           &
     &    LS_P                                                          &
!           Primed polar order
     &  , LS                                                            &
!           Polar order
     &  , LSR_P                                                         &
!           Reduced primed polar order (LSR_P is MS-1 less than
!           LS_P to facilitate addressing of arrays which do not
!           hold redundant space for m>l')
     &  , LSR                                                           &
!           Reduced polar order
     &  , LS_D                                                          &
!           Dummy polar order
     &  , LSR_D                                                         &
!           Reduced dummy polar order
     &  , LS_DD                                                         &
!           Dummy polar order
     &  , LSR_DD                                                        &
!           Reduced dummy polar order
     &  , I                                                             &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , L
!           Loop variable
      LOGICAL                                                           &
     &    L_ASSIGN
!           Controlling logical for assigning levels
      REAL  (Real64) ::                                                 &
     &    SS(ND_PROFILE, 0: ND_MAX_ORDER)                               &
!           S-coefficients for the current layer
     &  , SSRT(ND_PROFILE, 0: ND_MAX_ORDER)
!           Square roots of S-coefficients
      REAL  (Real64) ::                                                 &
     &    MU(ND_PROFILE, ND_RED_EIGENSYSTEM, 2)                         &
!           Eigenvaluse of the reduced system
     &  , EIG_VEC(ND_PROFILE, 2*ND_RED_EIGENSYSTEM                      &
     &      , ND_RED_EIGENSYSTEM, 2)                                    &
!           Eigenvectors of the full systems for positive eigenvalues
!           (these are scaled by the s-coefficients in the routine
!           EIG_SYS)
     &  , THETA(ND_PROFILE, ND_RED_EIGENSYSTEM, 2)                      &
!           Array of exponentials of optical depths along slant paths
     &  , SOURCE_TOP(ND_PROFILE, LS_TRUNC+1-MS, 2)                      &
!           Source function at the top of the layer
     &  , SOURCE_BOTTOM(ND_PROFILE, LS_TRUNC+1-MS, 2)
!           Source function at the bottom of the layer
      REAL  (Real64) ::                                                 &
     &    SURFACE_TERM(ND_PROFILE, LS_TRUNC+1-MS)                       &
!           Surface terms involving BRDFs
     &  , B_FACTOR(ND_PROFILE)                                          &
!           Contribution to the RHS of the equations
     &  , KSI                                                           &
!           Expression involving the BRDF
     &  , PHI                                                           &
!           Expression involving the BRDF
     &  , PHI_D                                                         &
!           Expression involving the BRDF
     &  , LAMBDA                                                        &
!           Expression involving the BRDF
     &  , LAMBDA_D
!           Expression involving the BRDF
      REAL  (Real64) ::                                                 &
     &    Z_SOL(ND_PROFILE, LS_TRUNC+1-MS)
!           Coefficient of the solar source function at the top of
!           the layer
      REAL  (Real64) ::                                                 &
     &    Q_0(ND_PROFILE)                                               &
!           Term for thermal particular integral
     &  , Q_1(ND_PROFILE)
!           Term for thermal particular integral
      INTEGER                                                           &
     &    K_SOL(ND_PROFILE)
!           Index of eigenvalue closest to the cosine of the solar
!           zenith angle
      REAL  (Real64) ::                                                 &
     &    UPM_C(ND_PROFILE, 2*ND_RED_EIGENSYSTEM)
!           Weights for exponentials in conditioning term
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    EIG_SYS, LAYER_PART_INTEG, SET_LEVEL_WEIGHTS                  &
     &  , SET_DIRN_WEIGHTS, CALC_SURF_RAD
!
!
!
!     Initialize the matrix.
      DO IE=1, 2*N_LAYER*N_RED_EIGENSYSTEM
        DO K=1, 6*N_RED_EIGENSYSTEM
          DO L=1, N_PROFILE
            A(L, IE, K)=0.0E+00_Real64
          ENDDO
        ENDDO
      ENDDO
!
!     To keep track of the layers in which radiances are required
!     I_ASSIGN_LEVEL is used: we search for the layer containing this
!     level, as indicated by I_RAD_LAYER and set the elements of
!     WEIGHT_U for later use with the vector giving the complementary
!     function. The terms of the particular integral are assigned to
!     C_YLM. Initialize to look for the first level.
      I_ASSIGN_LEVEL=1
!
!     I_BELOW and I_ABOVE hold variables for the layers below and
!     above the current interface. They are flipped to enable us
!     to use arrays with a dimension of 2, without the need to copy
!     lots of data.
      I_BELOW=1
!
!
!
!     Begin by determining the properties of the top layer.
      IF (MS == 0) THEN
        DO L=1, N_PROFILE
          SS(L, MS)=1.0E+00_Real64-OMEGA(L, 1)
          SSRT(L, MS)=SQRT(SS(L, MS))
        ENDDO
      ENDIF
      DO LS=MAX(1, MS), LS_TRUNC
        DO L=1, N_PROFILE
          SS(L, LS)=1.0E+00_Real64-OMEGA(L, 1)*PHASE_FNC(L, 1, LS)
          SSRT(L, LS)=SQRT(SS(L, LS))
        ENDDO
      ENDDO
!
!     Calculate the eigenvalues and eigenvectors for this layer.
! DEPENDS ON: eig_sys
      CALL EIG_SYS(N_PROFILE, LS_TRUNC, MS, N_RED_EIGENSYSTEM           &
     &  , CG_COEFF, SSRT(1, 0)                                          &
     &  , MU(1, 1, I_BELOW), EIG_VEC(1, 1, 1, I_BELOW)                  &
     &  , ND_PROFILE, ND_RED_EIGENSYSTEM, ND_MAX_ORDER                  &
     &  )
!
!     Calculate the exponential terms for this layer
      DO K=1, N_RED_EIGENSYSTEM
        DO L=1, N_PROFILE
          THETA(L, K, I_BELOW)=EXP(-TAU(L, 1)/MU(L, K, I_BELOW))
        ENDDO
      ENDDO
!
!     Find the particular integral in this layer.
! DEPENDS ON: layer_part_integ
      CALL LAYER_PART_INTEG(                                            &
     &    N_PROFILE, LS_TRUNC, MS, N_RED_EIGENSYSTEM                    &
     &  , CG_COEFF, MU(1, 1, I_BELOW)                                   &
     &  , EIG_VEC(1, 1, 1, I_BELOW), THETA(1, 1, I_BELOW)               &
     &  , ISOLIR, I_DIRECT(1, 0), I_DIRECT(1, 1), MU_0, UPLM_SOL        &
     &  , DIFF_PLANCK(1, 1), L_IR_SOURCE_QUAD, DIFF_PLANCK_2(1, 1)      &
     &  , TAU(1, 1), SS(1, 0)                                           &
     &  , SOURCE_TOP(1, 1, I_BELOW), SOURCE_BOTTOM(1, 1, I_BELOW)       &
     &  , UPM_C, K_SOL, Z_SOL, Q_0, Q_1                                 &
     &  , ND_PROFILE, ND_MAX_ORDER, ND_RED_EIGENSYSTEM                  &
     &  )
!
!
!     Impose Marshak's boundary conditions at the top of the atmosphere.
!     For each allowed order of l' (LS_P, or LSR_P in the reduced
!     notation), those with odd parity, the integral of Y_l'^m and the
!     boundary condition on the radiance is formed and integrated over
!     the downward hemisphere.
!
      DO LSR_P=2, LS_TRUNC+1-MS, 2
!
        IE=LSR_P/2
!
!       Begin with the exceptional case in which l=l' and KAPPA is 1/2.
        DO L=1, N_PROFILE
          B(L, IE)=-0.5E+00_Real64*SOURCE_TOP(L, LSR_P, I_BELOW)
        ENDDO
!       For other values of l, which must be odd when l' is even and
!       vice versa, the precalculated values are used. A hemispherically
!       isotropic incident radiance may exist if l=m=0, so we this
!       case exceptionally, adjusting the beginning of the loop.
        IF (MS == 0) THEN
          DO L=1, N_PROFILE
            B(L, IE)=B(L, IE)+KAPPA(LSR_P/2, 1)                         &
     &        *(2.0E+00_Real64*FLUX_DOWN_INC(L)/SQRT(PI)                &
     &        -SOURCE_TOP(L, 1, I_BELOW))
          ENDDO
        ENDIF
        DO LSR=MAX(3-2*MS, 1), LS_TRUNC-MS, 2
          DO L=1, N_PROFILE
            B(L, IE)=B(L, IE)                                           &
     &        -KAPPA(LSR_P/2, (LSR+1)/2)*SOURCE_TOP(L, LSR, I_BELOW)
          ENDDO
        ENDDO
!
!       Now calculate the coefficients of the matrix of unknowns,
!       u_{mik}^{+|-}.
        DO K=1, N_RED_EIGENSYSTEM
!         Variable numbers:
!         To accord with the general structure of the compressed matrix
!         the equations for the top boundary conditions are
!         right-justified by advancing the column by
!         2*N_RED_EIGENSYSTEM.
          IVMB=K+2*N_RED_EIGENSYSTEM
          IVPB=IVMB+N_RED_EIGENSYSTEM
!         In Marshak's procedure, l'+m will be odd, so apart from
!         the term where l'=l, l+m will be even in all the non-zero
!         terms of the sum over l, so it is easy to obtain
!         A(L, IE, IVMB) by a simple subtraction from A(L, IE, IVPB).
!         Begin with the term l=l'.
          DO L=1, N_PROFILE
            A(L, IE, IVPB)=0.5E+00_Real64*EIG_VEC(L, LSR_P, K, I_BELOW)
          ENDDO
          DO LS=MS, LS_TRUNC-1, 2
            LSR=LS+1-MS
            DO L=1, N_PROFILE
              A(L, IE, IVPB)                                            &
     &          =A(L, IE, IVPB)+KAPPA(LSR_P/2, (LSR+1)/2)               &
     &          *EIG_VEC(L, LSR, K, I_BELOW)
            ENDDO
          ENDDO
          DO L=1, N_PROFILE
            A(L, IE, IVMB)                                              &
     &        =A(L, IE, IVPB)-EIG_VEC(L, LSR_P, K, I_BELOW)
            A(L, IE, IVPB)=A(L, IE, IVPB)*THETA(L, K, I_BELOW)
          ENDDO
        ENDDO
      ENDDO
!
!     Set the weightings to be applied to the solution of the
!     linear system of equations.
      IF (I_SPH_ALGORITHM == IP_SPH_DIRECT) THEN
!       If we solve the problem directly the weightings will
!       apply to coefficients of the spherical harmonics.
!
!       The next test is done in two parts to ensure that it reamains
!       within bounds on I_RAD_LAYER.
        L_ASSIGN=(I_ASSIGN_LEVEL <= N_VIEWING_LEVEL)
        IF (L_ASSIGN) L_ASSIGN=(I_RAD_LAYER(I_ASSIGN_LEVEL) == 1)
!
! DEPENDS ON: set_level_weights
        CALL SET_LEVEL_WEIGHTS(1, N_PROFILE, LS_TRUNC                   &
     &    , MS, N_RED_EIGENSYSTEM                                       &
     &    , CG_COEFF, MU(1, 1, I_BELOW), EIG_VEC(1, 1, 1, I_BELOW)      &
     &    , ISOLIR, Z_SOL(1, 1), MU_0                                   &
     &    , Q_0, L_IR_SOURCE_QUAD, Q_1                                  &
     &    , UPM_C, K_SOL                                                &
     &    , TAU, SS                                                     &
     &    , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                &
     &    , L_ASSIGN, I_ASSIGN_LEVEL                                    &
     &    , C_YLM, WEIGHT_U(1, 1, 1, 1)                                 &
     &    , ND_PROFILE, ND_VIEWING_LEVEL                                &
     &    , ND_MAX_ORDER                                                &
     &    , ND_RED_EIGENSYSTEM, ND_SPH_CF_WEIGHT                        &
     &    )
      ELSE IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
!       Here the weights couple directly to radiances in
!       particular directions.
! DEPENDS ON: set_dirn_weights
        CALL SET_DIRN_WEIGHTS(N_PROFILE                                 &
     &    , MS, LS_TRUNC, UP_LM                                         &
     &    , N_DIRECTION, MU_V, AZIM_FACTOR                              &
     &    , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER, 1             &
     &    , N_RED_EIGENSYSTEM                                           &
     &    , MU(1, 1, I_BELOW), EIG_VEC(1, 1, 1, I_BELOW)                &
     &    , ISOLIR, Z_SOL(1, 1), MU_0                                   &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK                               &
     &    , UPM_C, K_SOL                                                &
     &    , TAU, OMEGA, PHASE_FNC                                       &
     &    , WEIGHT_U(1, 1, 1, 1), RADIANCE                              &
     &    , ND_PROFILE, ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL        &
     &    , ND_RED_EIGENSYSTEM, ND_MAX_ORDER                            &
     &    )
      ENDIF
!
!
!
!     For each interior level, 1,.., N_LAYER-1, continuity is imposed
!     on the (LS,MS)th component of the radiance field, for LS in the
!     range MS,..., LS_TRUNC (2*N_RED_EIGENSYSTEM orders). At each
!     stage we need information about the layers above and below the
!     layer. Arrays such as THETA therefore have an extra dimension of
!     size 2 to hold both values without the need to declare storage
!     for the whole column. To avoid copying values this last `index' is
!     accessed using the variables I_ABOVE and I_BELOW which are
!     flipped as we pass through each layer. To follow the indexing
!     note that when the loop variable is I we are looking at the
!     I-1st interface.
!
      DO I=2, N_LAYER
!
!       Flip the indices for the layer above and below.
        I_ABOVE=I_BELOW
        I_BELOW=3-I_BELOW
!
!       Calculate the condensed optical properties of the
!       current layer.
        IF (MS == 0) THEN
          DO L=1, N_PROFILE
            SS(L, MS)=1.0E+00_Real64-OMEGA(L, I)
            SSRT(L, MS)=SQRT(SS(L, MS))
          ENDDO
        ENDIF
        DO LS=MAX(1, MS), LS_TRUNC
          DO L=1, N_PROFILE
            SS(L, LS)=1.0E+00_Real64-OMEGA(L, I)*PHASE_FNC(L, I, LS)
            SSRT(L, LS)=SQRT(SS(L, LS))
          ENDDO
        ENDDO
!
!       Calculate the eigenvalues and eigenvectors for the current
!       layer which is that below the interface.
! DEPENDS ON: eig_sys
        CALL EIG_SYS(N_PROFILE, LS_TRUNC, MS, N_RED_EIGENSYSTEM         &
     &    , CG_COEFF, SSRT(1, 0)                                        &
     &    , MU(1, 1, I_BELOW), EIG_VEC(1, 1, 1, I_BELOW)                &
     &    , ND_PROFILE, ND_RED_EIGENSYSTEM, ND_MAX_ORDER                &
     &    )
!

!       Calculate the exponential terms for this layer
        DO K=1, N_RED_EIGENSYSTEM
          DO L=1, N_PROFILE
             THETA(L, K, I_BELOW)=EXP(-TAU(L, I)/MU(L, K, I_BELOW))
          ENDDO
        ENDDO
!
!       Find the particular integral in this layer.
! DEPENDS ON: layer_part_integ
        CALL LAYER_PART_INTEG(                                          &
     &      N_PROFILE, LS_TRUNC, MS, N_RED_EIGENSYSTEM                  &
     &    , CG_COEFF, MU(1, 1, I_BELOW)                                 &
     &    , EIG_VEC(1, 1, 1, I_BELOW), THETA(1, 1, I_BELOW)             &
     &    , ISOLIR, I_DIRECT(1, I-1), I_DIRECT(1, I), MU_0, UPLM_SOL    &
     &    , DIFF_PLANCK(1, I), L_IR_SOURCE_QUAD, DIFF_PLANCK_2(1, I)    &
     &    , TAU(1, I), SS(1, 0)                                         &
     &    , SOURCE_TOP(1, 1, I_BELOW), SOURCE_BOTTOM(1, 1, I_BELOW)     &
     &    , UPM_C, K_SOL, Z_SOL, Q_0, Q_1                               &
     &    , ND_PROFILE, ND_MAX_ORDER, ND_RED_EIGENSYSTEM                &
     &    )
!
!       Loop over the permitted orders of LS, compressing entries
!       into the matrix.
        DO LSR=1, 2*N_RED_EIGENSYSTEM
!
!         Number the equation:
          IF (MOD(LSR, 2) == 1) THEN
            IE=N_RED_EIGENSYSTEM*(2*I-3)+(LSR+1)/2
          ELSE IF (MOD(LSR, 2) == 0) THEN
            IE=N_RED_EIGENSYSTEM*(2*I-2)+LSR/2
          ENDIF
!
!         Loop over eigenvalues.
          DO K=1, N_RED_EIGENSYSTEM
!           Assign number to the variables in the equation
            IVMA=K
            IVPA=IVMA+N_RED_EIGENSYSTEM
            IVMB=IVPA+N_RED_EIGENSYSTEM
            IVPB=IVMB+N_RED_EIGENSYSTEM
            DO L=1, N_PROFILE
              A(L, IE, IVMA)=EIG_VEC(L, LSR, K, I_ABOVE)                &
     &          *THETA(L, K, I_ABOVE)*REAL(1-2*MOD(LSR-1, 2), Real64)
              A(L, IE, IVPA)=EIG_VEC(L, LSR, K, I_ABOVE)
              A(L, IE, IVMB)=-EIG_VEC(L, LSR, K, I_BELOW)               &
     &          *REAL(1-2*MOD(LSR-1, 2), Real64)
              A(L, IE, IVPB)=-EIG_VEC(L, LSR, K, I_BELOW)               &
     &          *THETA(L, K, I_BELOW)
            ENDDO
          ENDDO
!
          DO L=1, N_PROFILE
            B(L, IE)=SOURCE_TOP(L, LSR, I_BELOW)                        &
     &        -SOURCE_BOTTOM(L, LSR, I_ABOVE)
          ENDDO
!
        ENDDO
!
        IF (I_SPH_ALGORITHM == IP_SPH_DIRECT) THEN
!         If we solve the problem directly the weightings will
!         apply to coefficients of the spherical harmonics.
!         An assignment is required only if there are remaining
!         viewing levels and we are in the right layer.
!
!         The next test is done in two parts to ensure that it reamains
!         within bounds on I_RAD_LAYER.
          L_ASSIGN=(I_ASSIGN_LEVEL <= N_VIEWING_LEVEL)
          IF (L_ASSIGN) L_ASSIGN=(I_RAD_LAYER(I_ASSIGN_LEVEL) == I)
!
!         The different indexing of WEIGHT_U in the following two
!         calls is intentional. In the first case we interpolate
!         the radiance in one layer, so the final index runs only over
!         the eigensystem for that layer. In the second case, there
!         are contributions to the radiance at a particular level
!         from all layers, so the final index must be allow for
!         contributions from all layers.
!
! DEPENDS ON: set_level_weights
          CALL SET_LEVEL_WEIGHTS(I, N_PROFILE, LS_TRUNC                 &
     &      , MS, N_RED_EIGENSYSTEM                                     &
     &      , CG_COEFF, MU(1, 1, I_BELOW), EIG_VEC(1, 1, 1, I_BELOW)    &
     &      , ISOLIR, Z_SOL(1, 1), MU_0                                 &
     &      , Q_0, L_IR_SOURCE_QUAD, Q_1                                &
     &      , UPM_C, K_SOL                                              &
     &      , TAU(1, I), SS                                             &
     &      , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER              &
     &      , L_ASSIGN, I_ASSIGN_LEVEL                                  &
     &      , C_YLM, WEIGHT_U(1, 1, 1, 1)                               &
     &      , ND_PROFILE, ND_VIEWING_LEVEL                              &
     &      , ND_MAX_ORDER                                              &
     &      , ND_RED_EIGENSYSTEM, ND_SPH_CF_WEIGHT                      &
     &      )
        ELSE IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
!         Here the weights couple directly to radiances in
!         particular directions.
! DEPENDS ON: set_dirn_weights
          CALL SET_DIRN_WEIGHTS(N_PROFILE                               &
     &      , MS, LS_TRUNC, UP_LM                                       &
     &      , N_DIRECTION, MU_V, AZIM_FACTOR                            &
     &      , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER, I           &
     &      , N_RED_EIGENSYSTEM                                         &
     &      , MU(1, 1, I_BELOW), EIG_VEC(1, 1, 1, I_BELOW)              &
     &      , ISOLIR, Z_SOL(1, 1), MU_0                                 &
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK                             &
     &      , UPM_C, K_SOL                                              &
     &      , TAU, OMEGA, PHASE_FNC                                     &
     &      , WEIGHT_U(1, 1, 1, 1+2*N_RED_EIGENSYSTEM*(I-1)), RADIANCE  &
     &      , ND_PROFILE, ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL      &
     &      , ND_RED_EIGENSYSTEM, ND_MAX_ORDER                          &
     &      )
        ENDIF
!
!
      ENDDO
!
!
!
!
!     Impose the surface boundary condition using the appropriate
!     bidirectional reflection functions with Marshak's conditions.
!
!     Flip the `index' to the layer above the interface.
      I_ABOVE=I_BELOW
!
      DO LSR_P=2, LS_TRUNC+1-MS, 2
!
        LS_P=LSR_P+MS-1
!
        IE=N_RED_EIGENSYSTEM*(2*N_LAYER-1)+LSR_P/2
!
!       Initialize the RHS of the equations.
        DO L=1, N_PROFILE
          B(L, IE)=0.0E+00_Real64
        ENDDO
!
!       Terms in this equation fall into two groups: those which
!       involve the BRDF and those which do not. These latter
!       terms, which arise directly from Marshak's conditions
!       are treated first.
!
!       Begin with the exceptional case where l=l' so KAPPA is 1/2.
        DO L=1, N_PROFILE
          B(L,IE)=B(L,IE)-REAL(1-2*MOD(LS_P,2),Real64)*0.5E+00_Real64   &
     &      *SOURCE_BOTTOM(L, LSR_P, I_ABOVE)
        ENDDO
        IF ( (ISOLIR == IP_INFRA_RED).AND.(MS == 0).AND.                &
     &       (LSR_P == 1) ) THEN
!         The Planckian flux is used instead of the radiance for
!         consistency with the two-stream equations.
          DO L=1, N_PROFILE
            B(L, IE)=B(L, IE)+D_PLANCK_FLUX_SURFACE(L)/SQRT(PI)
          ENDDO
        ENDIF
!       For other values of l, which must be odd when l' is even and
!       vice versa, the precalculated values are used.
        DO LSR=1, LS_TRUNC-MS, 2
          DO L=1, N_PROFILE
            B(L, IE)=B(L, IE)                                           &
     &        -REAL(1-2*MOD(LSR+MS-1, 2), Real64)                       &
     &        *KAPPA(LSR_P/2, (LSR+1)/2)                                &
     &        *SOURCE_BOTTOM(L, LSR, I_ABOVE)
          ENDDO
          IF ( (ISOLIR == IP_INFRA_RED).AND.(MS == 0).AND.              &
     &         (LSR == 1) ) THEN
            DO L=1, N_PROFILE
              B(L, IE)=B(L, IE)+KAPPA(LSR_P/2, (LSR+1)/2)               &
     &          *D_PLANCK_FLUX_SURFACE(L)*2.0E+00_Real64/SQRT(PI)
            ENDDO
          ENDIF
        ENDDO
!
!       Now calculate the coefficients of the matrix of unknowns,
!       u_{mik}^{+|-}.
        DO K=1, N_RED_EIGENSYSTEM
!         Variable numbers:
          IVMA=K
          IVPA=IVMA+N_RED_EIGENSYSTEM
!         KAPPA has not been calculated for those values of l for
!         which it is 0: we therefore add the terms in two groups,
!         firstly those for l=l' and then those for other values of
!         l where KAPPA is non-zero. As at the top, of the atmosphere
!         it is possible to evaluate A(L, IE, IVMA) from
!         A(L, IE, IVPA)
          LS_P=LSR_P+MS-1
          DO L=1, N_PROFILE
            A(L, IE, IVPA)=A(L, IE, IVPA)+REAL(1-2*MOD(LS_P, 2), Real64)&
     &        *0.5E+00_Real64*EIG_VEC(L, LSR_P, K, I_ABOVE)
          ENDDO
          DO LS=MS, LS_TRUNC-1, 2
            LSR=LS+1-MS
            DO L=1, N_PROFILE
              A(L, IE, IVPA)=A(L, IE, IVPA)+REAL(1-2*MOD(LS, 2), Real64)&
     &          *KAPPA(LSR_P/2, (LSR+1)/2)*EIG_VEC(L, LSR, K, I_ABOVE)
            ENDDO
          ENDDO
!
          DO L=1, N_PROFILE
            A(L, IE, IVMA)=A(L, IE, IVPA)+REAL(1-2*MOD(MS, 2), Real64)  &
     &          *EIG_VEC(L, LSR_P, K, I_ABOVE)
          ENDDO
!
        ENDDO
!
!
!       The second group of terms involves the BRDF.
!       There will be no contribution from orders
!       above the order of trunction of the BRDF.
        IF (MS <= LS_BRDF_TRUNC) THEN
!         Add in the solar or infra-red contributions involving the
!         BRDF basis functions which do not involve terms in KAPPA.
          IF (ISOLIR == IP_SOLAR) THEN
!
            DO J=1, N_BRDF_BASIS_FNC
              DO L=1, N_PROFILE
                B_FACTOR(L)=0.0E+00_Real64
              ENDDO
              DO LS=MS, LS_BRDF_TRUNC-MOD(MS, 2), 2
                LSR=LS+1-MS
                KSI=KAPPA(LSR_P/2, 1)*F_BRDF(J, 0, LS/2, MS)
                DO LS_D=MS+2, LS_BRDF_TRUNC-MOD(MS, 2), 2
                  LSR_D=LS_D-MS+1
                  KSI=KSI+KAPPA(LSR_P/2, (LSR_D+1)/2)                   &
     &              *F_BRDF(J, LS_D/2, LS/2, MS)
                ENDDO
                DO L=1, N_PROFILE
                  B_FACTOR(L)=B_FACTOR(L)+KSI*UPLM_SOL(L, LSR)
                ENDDO
              ENDDO
              DO L=1, N_PROFILE
                B(L, IE)=B(L, IE)+I_DIRECT(L, N_LAYER)*MU_0(L)          &
     &            *REAL(1-2*MOD(MS, 2), Real64)                         &
     &            *RHO_ALB(L, J)*B_FACTOR(L)
              ENDDO
            ENDDO
!
          ELSE IF (ISOLIR == IP_INFRA_RED) THEN
            IF (MS == 0) THEN
              DO J=1, N_BRDF_BASIS_FNC
                LAMBDA=0.0E+00_Real64
                DO LS_D=0, LS_BRDF_TRUNC, 2
                  LSR_D=LS_D+1
                  LAMBDA_D=0.0E+00_Real64
                  DO LS_DD=0, LS_BRDF_TRUNC, 2
                    LSR_DD=LS_DD+1
                    LAMBDA_D=LAMBDA_D+KAPPA(LSR_P/2, (LSR_DD+1)/2)      &
     &                *F_BRDF(J, LS_DD/2, LS_D/2, MS)
                  ENDDO
                  LAMBDA=LAMBDA+KAPPA(1, (LSR_D+1)/2)*LAMBDA_D
                ENDDO
                DO L=1, N_PROFILE
                  B(L, IE)=B(L, IE)                                     &
     &              +RHO_ALB(L, J)*LAMBDA                               &
     &              *SQRT(4.0E+00_Real64*PI/3.0E+00_Real64)             &
     &              *D_PLANCK_FLUX_SURFACE(L)/PI
                ENDDO
              ENDDO
            ENDIF
          ENDIF
!
          DO LS=MS, LS_TRUNC
!
            LSR=LS+1-MS
!
            DO L=1, N_PROFILE
              SURFACE_TERM(L, LSR)=0.0E+00_Real64
            ENDDO
            DO J=1, N_BRDF_BASIS_FNC
              PHI=0.0E+00_Real64
              DO LS_D=MS, LS_BRDF_TRUNC-MOD(MS, 2), 2
                LSR_D=LS_D-MS+1
                PHI_D=0.0E+00_Real64
                DO LS_DD=MS, LS_BRDF_TRUNC-MOD(MS, 2), 2
                  LSR_DD=LS_DD-MS+1
                  PHI_D=PHI_D+CGK((LSR_DD+1)/2, LSR)                    &
     &              *F_BRDF(J, LS_D/2, LS_DD/2, MS)
                ENDDO
                PHI=PHI+KAPPA(LSR_P/2, (LSR_D+1)/2)*PHI_D
              ENDDO
              DO L=1, N_PROFILE
                SURFACE_TERM(L, LSR)=SURFACE_TERM(L, LSR)               &
     &            +RHO_ALB(L, J)*PHI*REAL(1-2*MOD(MS, 2), Real64)
              ENDDO
            ENDDO
!
!           Add on the contribution to the RHS.
            DO L=1, N_PROFILE
              B(L, IE)=B(L, IE)                                         &
     &          -SOURCE_BOTTOM(L, LSR, I_ABOVE)*SURFACE_TERM(L, LSR)
            ENDDO
!           Add in the contributions to the matrix on the LHS.
            DO K=1, N_RED_EIGENSYSTEM
!             Variable numbers:
              IVMA=K
              IVPA=IVMA+N_RED_EIGENSYSTEM
              DO L=1, N_PROFILE
                A(L, IE, IVMA)=A(L, IE, IVMA)                           &
     &            +SURFACE_TERM(L, LSR)*REAL(1-2*MOD(LSR-1, 2), Real64) &
     &            *EIG_VEC(L, LSR, K, I_ABOVE)
                A(L, IE, IVPA)=A(L, IE, IVPA)                           &
     &            +SURFACE_TERM(L, LSR)*EIG_VEC(L, LSR, K, I_ABOVE)
              ENDDO
!
            ENDDO
!
          ENDDO
!
        ENDIF
!
        DO K=1, N_RED_EIGENSYSTEM
          IVMA=K
          DO L=1, N_PROFILE
            A(L, IE, IVMA)=A(L, IE, IVMA)*THETA(L, K, I_ABOVE)
          ENDDO
        ENDDO
!
      ENDDO
!
!
      IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
!       Calculate the contribution of radiation reflected from the
!       surface.

! DEPENDS ON: calc_surf_rad
        CALL CALC_SURF_RAD(N_PROFILE, N_LAYER, TAU                      &
     &    , MS, LS_TRUNC, EULER_FACTOR                                  &
     &    , ISOLIR, I_DIRECT(1, N_LAYER), MU_0, D_PLANCK_FLUX_SURFACE   &
     &    , N_BRDF_BASIS_FNC, LS_BRDF_TRUNC, F_BRDF                     &
     &    , RHO_ALB, BRDF_SOL, BRDF_HEMI, CGK                           &
     &    , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                &
     &    , N_DIRECTION, MU_V, UP_LM, AZIM_FACTOR                       &
     &    , N_RED_EIGENSYSTEM, EIG_VEC(1, 1, 1, I_ABOVE)                &
     &    , THETA(1, 1, I_ABOVE), SOURCE_BOTTOM(1, 1, I_ABOVE)          &
     &    , RADIANCE                                                    &
     &    , WEIGHT_U(1, 1, 1, 1+2*N_RED_EIGENSYSTEM*(N_LAYER-1))        &
     &    , ND_PROFILE, ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL        &
     &    , ND_RED_EIGENSYSTEM, ND_MAX_ORDER, ND_BRDF_BASIS_FNC         &
     &    , ND_BRDF_TRUNC                                               &
     &    )

!       Isotropic incident fluxes are permitted (and required in the
!       differential formulation of the IR).
        IF (MS == 0) THEN
! DEPENDS ON: calc_top_rad
          CALL CALC_TOP_RAD(N_PROFILE, TAU                              &
     &      , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER              &
     &      , N_DIRECTION, MU_V                                         &
     &      , FLUX_DOWN_INC                                             &
     &      , RADIANCE                                                  &
     &      , ND_PROFILE, ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL      &
     &      )
        ENDIF
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE BUILD_SPH_MATRIX
#endif
