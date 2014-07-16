#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve for the monochromatic radiances.
!
! Method:
!       The final single scattering properties are calculated
!       and rescaled. An appropriate subroutine is called to
!       calculate the radiances depending on the treatment of
!       cloudiness.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MONOCHROMATIC_RADIANCE_SPH(IERR                        &
!                     Atmospheric Propertries
     &  , N_PROFILE, N_LAYER, D_MASS                                    &
!                     Angular Integration
     &  , N_ORDER_PHASE, MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC   &
     &  , ACCURACY_ADAPTIVE, EULER_FACTOR, I_SPH_ALGORITHM              &
     &  , I_SPH_MODE, L_RESCALE                                         &
!                       Precalculated angular arrays
     &  , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                      &
!                     Options for Equivalent Extinction
     &  , L_SCALE_SOLAR, ADJUST_SOLAR_KE                                &
!                     Spectral Region
     &  , ISOLIR                                                        &
!                     Infra-red Properties
     &  , DIFF_PLANCK, L_IR_SOURCE_QUAD, DIFF_PLANCK_2                  &
!                     Conditions at TOA
     &  , ZEN_0, FLUX_INC_DIRECT, FLUX_INC_DOWN                         &
     &  , I_DIRECT                                                      &
!                     Surface Properties
     &  , D_PLANCK_FLUX_SURFACE                                         &
     &  , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB                      &
     &  , F_BRDF, BRDF_SOL, BRDF_HEMI                                   &
!                     Optical Properties
     &  , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR                             &
     &  , FORWARD_SCATTER_CLR, PHASE_FNC_SOLAR_CLR                      &
     &  , TAU, OMEGA, PHASE_FNC                                         &
     &  , FORWARD_SCATTER, PHASE_FNC_SOLAR                              &
!                     Cloudy Properties
     &  , L_CLOUD, I_CLOUD                                              &
!                     Cloud Geometry
     &  , N_CLOUD_TOP                                                   &
     &  , N_CLOUD_TYPE, FRAC_CLOUD                                      &
     &  , N_REGION, K_CLR, I_REGION_CLOUD, FRAC_REGION                  &
     &  , W_FREE, W_CLOUD, CLOUD_OVERLAP                                &
     &  , N_COLUMN_SLV, LIST_COLUMN_SLV                                 &
     &  , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN                     &
!                       Levels for calculating radiances
     &  , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                  &
!                       Viewing geometry
     &  , N_DIRECTION, DIRECTION                                        &
!                     Calculated Fluxes
     &  , FLUX_DIRECT, FLUX_TOTAL                                       &
!                       Calculated radiances
     &  , RADIANCE                                                      &
!                       Calculated mean radiances
     &  , J_RADIANCE                                                    &
!                     Dimensions of Arrays
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN          &
     &  , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE            &
     &  , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF                    &
     &  , ND_MAX_ORDER, ND_SPH_COEFF                                    &
     &  , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC, ND_VIEWING_LEVEL            &
     &  , ND_DIRECTION, ND_SOURCE_COEFF                                 &
     &  )
!
!
      IMPLICIT NONE
!
! Include Header files
#include "c_kinds.h"
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ND_LAYER                                                      &
!           Maximum number of layers
     &  , ND_LAYER_CLR                                                  &
!           Maximum number of completely clear layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_FLUX_PROFILE                                               &
!           Maximum number of profiles in arrays of fluxes
     &  , ND_RADIANCE_PROFILE                                           &
!           Maximum number of profiles in arrays of radiances
     &  , ND_J_PROFILE                                                  &
!           Maximum number of profiles in arrays of mean radiances
     &  , ND_COLUMN                                                     &
!           Number of columns per point
     &  , ND_CLOUD_TYPE                                                 &
!           Maximum number of types of cloud
     &  , ND_REGION                                                     &
!           Maximum number of cloudy regions
     &  , ND_OVERLAP_COEFF                                              &
!           Maximum number of overlap coefficients
     &  , ND_MAX_ORDER                                                  &
!           Maximum order of spherical harmonics used
     &  , ND_SPH_COEFF                                                  &
!           Allocated size for spherical coefficients
     &  , ND_BRDF_BASIS_FNC                                             &
!           Size allowed for BRDF basis functions
     &  , ND_BRDF_TRUNC                                                 &
!           Size allowed for orders of BRDFs
     &  , ND_VIEWING_LEVEL                                              &
!           Allocated size for levels where radiances are calculated
     &  , ND_DIRECTION                                                  &
!           Allocated size for viewing directions
     &  , ND_SOURCE_COEFF
!           Size allocated for source coefficients
!
!     Include header files.
#include "def_std_io_icf3z.h"
#include "cloud_scheme_pcf3z.h"
#include "angular_integration_pcf3z.h"
#include "sph_algorithm_pcf3z.h"
#include "solver_pcf3z.h"
#include "surface_spec_pcf3z.h"
#include "spectral_region_pcf3z.h"
#include "error_pcf3z.h"
!
!
!
!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
!
!                     Atmospheric properties
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER
!           Number of layers
      REAL  (Real64), INTENT(IN) ::                                     &
     &    D_MASS(ND_PROFILE, ND_LAYER)
!           Mass thickness of each layer
!
!                     Angular integration
      INTEGER, INTENT(IN) ::                                            &
     &    N_ORDER_PHASE                                                 &
!           Highest order retained in the phase function
     &  , I_TRUNCATION                                                  &
!           Type of spherical truncation adopted
     &  , MS_MIN                                                        &
!           Lowest azimuthal order calculated
     &  , MS_MAX                                                        &
!           Highest azimuthal order calculated
     &  , IA_SPH_MM(0: ND_MAX_ORDER)                                    &
!           Address of spherical coefficient for (m, m) for each m
     &  , LS_LOCAL_TRUNC(0: ND_MAX_ORDER)                               &
!           Orders of truncation at each azimuthal order
     &  , I_SPH_MODE                                                    &
!           Mode in which the spherical harmonic solver is being used
     &  , I_SPH_ALGORITHM
!           Algorithm used for spherical harmonic calculation
      LOGICAL, INTENT(IN) ::                                            &
     &    L_RESCALE
!           Flag for rescaling of the optical properties
      REAL  (Real64) ::                                                 &
     &    CG_COEFF(ND_SPH_COEFF)                                        &
!           Clebsch-Gordan coefficients
     &  , UPLM_ZERO(ND_SPH_COEFF)                                       &
!           Values of spherical harmonics at polar angles pi/2
     &  , UPLM_SOL(ND_PROFILE, ND_SPH_COEFF)
!           Values of spherical harmonics in the solar direction
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ACCURACY_ADAPTIVE                                             &
!           Accuracy for adaptive truncation
     &  , EULER_FACTOR
!           Factor applied to the last term of an alternating series
!
!                     Variables for equivalent extinction
      LOGICAL, INTENT(IN) ::                                            &
     &    L_SCALE_SOLAR
!           Apply scaling to solar flux
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ADJUST_SOLAR_KE(ND_PROFILE, ND_LAYER)
!           Adjustment of solar beam with equivalent extinction
!
!                     Spectral region
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Visible or IR
!
!                     Infra-red properties
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Flag for quadratic IR-source
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DIFF_PLANCK(ND_PROFILE, ND_LAYER)                             &
!           DIfferences in the Planckian function across layers
     &  , DIFF_PLANCK_2(ND_PROFILE, ND_LAYER)
!           Twice the second differences of Planckian source function
!
!                     Conditions at TOA
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ZEN_0(ND_PROFILE)                                             &
!           Secants (two-stream) or cosines (spherical harmonics)
!           of the solar zenith angles
     &  , FLUX_INC_DIRECT(ND_PROFILE)                                   &
!           Incident direct flux
     &  , FLUX_INC_DOWN(ND_PROFILE)
!           Incident downward flux
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    I_DIRECT(ND_PROFILE, 0: ND_LAYER)
!           Direct radiance (the first row contains the incident
!           solar radiance: the other rows are calculated)
!
!                     Surface properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    D_PLANCK_FLUX_SURFACE(ND_PROFILE)
!           Differential Planckian flux from the surface
      INTEGER, INTENT(IN) ::                                            &
     &    LS_BRDF_TRUNC                                                 &
!           Order of trunation of BRDFs
     &  , N_BRDF_BASIS_FNC
!           Number of BRDF basis functions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    RHO_ALB(ND_PROFILE, ND_BRDF_BASIS_FNC)                        &
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
!                     Optical properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    PHASE_FNC_CLR(ND_PROFILE, ND_LAYER_CLR, ND_MAX_ORDER)         &
!           Clear-sky phase function
     &  , FORWARD_SCATTER_CLR(ND_PROFILE, ND_LAYER_CLR)                 &
!           Clear-sky forward scattering
     &  , PHASE_FNC_SOLAR_CLR(ND_RADIANCE_PROFILE, ND_LAYER_CLR         &
     &      , ND_DIRECTION)
!           Clear-sky solar phase fuction in viewing directions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    PHASE_FNC(ND_PROFILE, ID_CT: ND_LAYER                         &
     &      , ND_MAX_ORDER, 0: ND_CLOUD_TYPE)                           &
!           Phase functions
     &  , FORWARD_SCATTER(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)&
!           Forward scattering
     &  , PHASE_FNC_SOLAR(ND_RADIANCE_PROFILE, ID_CT: ND_LAYER          &
     &      , ND_DIRECTION, 0: ND_CLOUD_TYPE)
!           Phase function for the solar beam in viewing
!           directions
!
!                     Cloudy properties
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLOUD
!           Clouds required
      INTEGER, INTENT(IN) ::                                            &
     &    I_CLOUD
!           Cloud scheme used
!
!                     Cloud geometry
      INTEGER, INTENT(IN) ::                                            &
     &    N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_CLOUD_TYPE                                                  &
!           Number of types of clouds
     &  , N_REGION                                                      &
!           Number of cloudy regions
     &  , K_CLR                                                         &
!           Index of clear-sky region
     &  , I_REGION_CLOUD(ND_CLOUD_TYPE)
!           Regions in which types of clouds fall
!
!     Cloud geometry
      INTEGER, INTENT(IN) ::                                            &
     &    N_COLUMN_SLV(ND_PROFILE)                                      &
!           Number of columns to be solved in each profile
     &  , LIST_COLUMN_SLV(ND_PROFILE, ND_COLUMN)                        &
!           List of columns requiring an actual solution
     &  , I_CLM_LYR_CHN(ND_PROFILE, ND_COLUMN)                          &
!           Layer in the current column to change
     &  , I_CLM_CLD_TYP(ND_PROFILE, ND_COLUMN)
!           Type of cloud to introduce in the changed layer
      REAL  (Real64), INTENT(IN) ::                                     &
     &    W_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)                          &
!           Cloudy fraction
     &  , FRAC_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_CLOUD_TYPE)        &
!           Fractions of different types of cloud
     &  , W_FREE(ND_PROFILE, ID_CT: ND_LAYER)                           &
!           Clear-sky fraction
     &  , CLOUD_OVERLAP(ND_PROFILE, ID_CT-1: ND_LAYER                   &
     &      , ND_OVERLAP_COEFF)                                         &
!           Coefficients for energy transfer at interfaces
     &  , AREA_COLUMN(ND_PROFILE, ND_COLUMN)                            &
!           Areas of columns
     &  , FRAC_REGION(ND_PROFILE, ID_CT: ND_LAYER, ND_REGION)
!           Fractions of total cloud occupied by each region
!
!
!                     Levels where radiance are calculated
      INTEGER, INTENT(IN) ::                                            &
     &    N_VIEWING_LEVEL                                               &
!           Number of levels where radiances are calculated
     &  , I_RAD_LAYER(ND_VIEWING_LEVEL)
!           Layers in which radiances are calculated
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FRAC_RAD_LAYER(ND_VIEWING_LEVEL)
!           Fractions below the tops of the layers
!
!                       Viewing Geometry
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION
!           Number of viewing directions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DIRECTION(ND_RADIANCE_PROFILE, ND_DIRECTION, 2)
!           Viewing directions
!
!                     Calculated Fluxes
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_FLUX_PROFILE, 0: ND_LAYER)                     &
!           Direct flux
     &  , FLUX_TOTAL(ND_FLUX_PROFILE, 2*ND_LAYER+2)
!           Total flux
!
!                     Calculated radiances
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    RADIANCE(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL, ND_DIRECTION)
!           Radiances
!
!                     Calculated mean radiances
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    J_RADIANCE(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Mean radiances
!
!
!     Local variables.
      INTEGER                                                           &
     &    ND_RED_EIGENSYSTEM                                            &
!           Size allowed for the reduced eigensystem
     &  , ND_SPH_EQUATION                                               &
!           Size allowed for spherical harmonic equations
     &  , ND_SPH_DIAGONAL                                               &
!           Size allowed for diagonals of the spherical harmonic
!           matrix
     &  , ND_SPH_CF_WEIGHT                                              &
!           Size allowed for entities to be incremented by the
!           complementary function of the linear system
     &  , ND_SPH_U_RANGE                                                &
!           Size allowed for range of values of u^+|- contributing
!           on any viewing level
     &  , ND_PROFILE_COLUMN
!           Size allowed for profiles taken simultaneously in a
!           decomposition into columns
      INTEGER                                                           &
     &    L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    TAU_CLR(ND_PROFILE, ND_LAYER_CLR)                             &
!           Clear-sky optical depth
     &  , OMEGA_CLR(ND_PROFILE, ND_LAYER_CLR)                           &
!           Clear-sky albedo of single scattering
     &  , TAU(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)            &
!           Optical depth
     &  , OMEGA(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)
!           Single scattering albedo
!
      REAL  (Real64), ALLOCATABLE ::                                    &
     &    TAU_CLR_F(:, :)                                               &
!           Clear-sky optical depth for the whole column
     &  , OMEGA_CLR_F(:, :)                                             &
!           Clear-sky albedo of single scattering for the whole column
     &  , PHASE_FNC_CLR_F(:, :, :)                                      &
!           Clear-sky phase function for the whole column
     &  , FORWARD_SCATTER_CLR_F(:, :)                                   &
!           Clear-sky forward scattering for the whole column
     &  , PHASE_FNC_SOLAR_CLR_F(:, :, :)
!           Clear-sky solar phase function in viewing directions
!           for the whole column
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    COPY_CLR_FULL, COPY_CLR_SOL                                   &
     &  , SPH_SOLVER, CALC_RADIANCE_IPA
!
!
!
!     Split the method os folution according to the cloud scheme.
      IF (I_CLOUD == IP_CLOUD_CLEAR) THEN
!
!       Precalculate dimensions for the dynamically allocated
!       arrays.
        ND_RED_EIGENSYSTEM=(ND_MAX_ORDER+1)/2
        ND_SPH_EQUATION=2*ND_LAYER*ND_RED_EIGENSYSTEM
        ND_SPH_DIAGONAL=6*ND_RED_EIGENSYSTEM
        IF (I_SPH_ALGORITHM == IP_SPH_DIRECT) THEN
          ND_SPH_CF_WEIGHT=ND_MAX_ORDER+1
          ND_SPH_U_RANGE=2*ND_RED_EIGENSYSTEM
        ELSE IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
          ND_SPH_CF_WEIGHT=ND_DIRECTION
          ND_SPH_U_RANGE=ND_SPH_EQUATION
        ENDIF
!
!       Allocate and set dynamic arrays.
        ALLOCATE(TAU_CLR_F(ND_PROFILE, ND_LAYER))
        ALLOCATE(OMEGA_CLR_F(ND_PROFILE, ND_LAYER))
        ALLOCATE(PHASE_FNC_CLR_F(ND_PROFILE, ND_LAYER, ND_MAX_ORDER))
        ALLOCATE(FORWARD_SCATTER_CLR_F(ND_PROFILE, ND_LAYER))
        ALLOCATE(PHASE_FNC_SOLAR_CLR_F(ND_RADIANCE_PROFILE              &
     &    , ND_LAYER, ND_DIRECTION))
!
! DEPENDS ON: copy_clr_full
        CALL COPY_CLR_FULL(N_PROFILE, N_LAYER, N_CLOUD_TOP              &
     &    , N_ORDER_PHASE                                               &
     &    , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR                           &
     &    , TAU, OMEGA, PHASE_FNC                                       &
     &    , TAU_CLR_F, OMEGA_CLR_F, PHASE_FNC_CLR_F                     &
!                       Sizes of arrays
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_MAX_ORDER     &
     &    )
        IF ( (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER).AND.              &
     &       (ISOLIR == IP_SOLAR) ) THEN
! DEPENDS ON: copy_clr_sol
          CALL COPY_CLR_SOL(N_PROFILE, N_LAYER, N_CLOUD_TOP             &
     &      , N_DIRECTION, L_RESCALE                                    &
     &      , FORWARD_SCATTER_CLR                                       &
     &      , PHASE_FNC_SOLAR_CLR                                       &
     &      , FORWARD_SCATTER, PHASE_FNC_SOLAR                          &
     &      , FORWARD_SCATTER_CLR_F                                     &
     &      , PHASE_FNC_SOLAR_CLR_F                                     &
!                       Sizes of arrays
     &      , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_DIRECTION   &
     &      )
        ENDIF
!
! DEPENDS ON: sph_solver
        CALL SPH_SOLVER(IERR                                            &
!                       Atmospheric sizes
     &    , N_PROFILE, N_LAYER                                          &
!                       Angular integration
     &    , MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC                &
     &    , CG_COEFF, UPLM_ZERO, IA_SPH_MM                              &
     &    , ACCURACY_ADAPTIVE, EULER_FACTOR                             &
     &    , I_SPH_ALGORITHM, I_SPH_MODE, L_RESCALE                      &
!                       Spectral Region
     &    , ISOLIR                                                      &
!                     Options for Equivalent Extinction
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
!                       Solar Fields
     &    , I_DIRECT, ZEN_0, UPLM_SOL                                   &
!                       Infra-red Properties
     &    , DIFF_PLANCK, FLUX_INC_DOWN                                  &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                             &
!                       Optical properies
     &    , TAU_CLR_F, OMEGA_CLR_F, PHASE_FNC_CLR_F                     &
     &    , PHASE_FNC_SOLAR_CLR_F, FORWARD_SCATTER_CLR_F                &
!                       Surface Conditions
     &    , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB                    &
     &    , F_BRDF, BRDF_SOL, BRDF_HEMI                                 &
     &    , D_PLANCK_FLUX_SURFACE                                       &
!                       Levels for calculating radiances
     &    , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                &
!                       Viewing Geometry
     &    , N_DIRECTION, DIRECTION                                      &
!                       Radiances Calculated
     &    , FLUX_DIRECT, FLUX_TOTAL, RADIANCE, J_RADIANCE               &
!                       Dimensions of arrays
     &    , ND_PROFILE, ND_LAYER                                        &
     &    , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE          &
     &    , ND_MAX_ORDER, ND_SPH_COEFF                                  &
     &    , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC                            &
     &    , ND_RED_EIGENSYSTEM, ND_SPH_EQUATION, ND_SPH_DIAGONAL        &
     &    , ND_SPH_CF_WEIGHT, ND_SPH_U_RANGE                            &
     &    , ND_VIEWING_LEVEL, ND_DIRECTION                              &
     &    )
!
        IF (IERR /= I_NORMAL) RETURN
!
      ELSEIF ( (I_CLOUD == IP_CLOUD_MIX_MAX).OR.                        &
     &         (I_CLOUD == IP_CLOUD_MIX_RANDOM).OR.                     &
     &         (I_CLOUD == IP_CLOUD_TRIPLE).OR.                         &
     &         (I_CLOUD == IP_CLOUD_PART_CORR).OR.                      &
     &         (I_CLOUD == IP_CLOUD_PART_CORR_CNV) ) THEN
!
        WRITE(IU_ERR, '(/A)')                                           &
     &    '*** ERROR: Radiances cannot yet be computed using '          &
     &    //'coupled overlaps.'
        IERR=I_ERR_FATAL
        RETURN
!
      ELSEIF (I_CLOUD == IP_CLOUD_COLUMN_MAX) THEN
!
!       Clouds are treated using the independent pixel approximation,
!       as directed by the decompositional arrays.
!
!       Set a dimension to allow the subcolumns of several profiles
!       to be considered at once.
        ND_PROFILE_COLUMN=MAX(1, N_PROFILE)
        DO L=1, N_PROFILE
          ND_PROFILE_COLUMN=MAX(ND_PROFILE_COLUMN, N_COLUMN_SLV(L))
        ENDDO
!
!       Precalculate dimensions for the dynamically allocated
!       arrays.
        ND_RED_EIGENSYSTEM=(ND_MAX_ORDER+1)/2
        ND_SPH_EQUATION=2*ND_LAYER*ND_RED_EIGENSYSTEM
        ND_SPH_DIAGONAL=6*ND_RED_EIGENSYSTEM
        IF (I_SPH_ALGORITHM == IP_SPH_DIRECT) THEN
          ND_SPH_CF_WEIGHT=ND_MAX_ORDER+1
          ND_SPH_U_RANGE=2*ND_RED_EIGENSYSTEM
        ELSE IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
          ND_SPH_CF_WEIGHT=ND_DIRECTION
          ND_SPH_U_RANGE=ND_SPH_EQUATION
        ENDIF
!
!
! DEPENDS ON: calc_radiance_ipa
        CALL CALC_RADIANCE_IPA(IERR                                     &
!                     Atmospheric Properties
     &    , N_PROFILE, N_LAYER, N_CLOUD_TOP                             &
!                       Angular Integration
     &    , N_ORDER_PHASE, MS_MIN, MS_MAX, LS_LOCAL_TRUNC               &
     &    , I_TRUNCATION, ACCURACY_ADAPTIVE, EULER_FACTOR               &
     &    , I_SPH_ALGORITHM, I_SPH_MODE, L_RESCALE                      &
!                       Precalculated angular arrays
     &    , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                    &
!                     Options for Equivalent Extinction
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
!                     Spectral Region
     &    , ISOLIR                                                      &
!                     Infra-red Properties
     &    , DIFF_PLANCK                                                 &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                             &
!                     Conditions at TOA
     &    , FLUX_INC_DOWN, ZEN_0                                        &
!                     Conditions at Surface
     &    , D_PLANCK_FLUX_SURFACE                                       &
     &    , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB                    &
     &    , F_BRDF, BRDF_SOL, BRDF_HEMI                                 &
!                     Clear-sky Single Scattering Properties
     &    , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR, PHASE_FNC_SOLAR_CLR      &
     &    , FORWARD_SCATTER_CLR                                         &
     &    , TAU, OMEGA, PHASE_FNC, PHASE_FNC_SOLAR                      &
     &    , FORWARD_SCATTER                                             &
!                     Cloud Geometry
     &    , N_COLUMN_SLV, LIST_COLUMN_SLV                               &
     &    , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN                   &
!                       Levels for calculating radiances
     &    , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                &
!                       Viewing Geometry
     &    , N_DIRECTION, DIRECTION                                      &
!                       Calculated fluxes or radiances
     &    , FLUX_DIRECT, FLUX_TOTAL, I_DIRECT, RADIANCE, J_RADIANCE     &
!                     Dimensions of Arrays
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT                   &
     &    , ND_COLUMN, ND_CLOUD_TYPE                                    &
     &    , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE          &
     &    , ND_MAX_ORDER, ND_SPH_COEFF                                  &
     &    , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC                            &
     &    , ND_RED_EIGENSYSTEM, ND_SPH_EQUATION, ND_SPH_DIAGONAL        &
     &    , ND_SPH_CF_WEIGHT, ND_SPH_U_RANGE                            &
     &    , ND_VIEWING_LEVEL, ND_DIRECTION                              &
     &    , ND_PROFILE_COLUMN                                           &
     &    )
!
        IF (IERR /= I_NORMAL) RETURN
!
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE MONOCHROMATIC_RADIANCE_SPH
#endif
