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
      SUBROUTINE MONOCHROMATIC_RADIANCE(IERR                            &
!                     Atmospheric Propertries
     &  , N_PROFILE, N_LAYER, D_MASS                                    &
!                     Angular Integration
     &  , I_ANGULAR_INTEGRATION, I_2STREAM                              &
     &  , L_RESCALE, N_ORDER_GAUSS                                      &
     &  , N_ORDER_PHASE, MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC   &
     &  , ACCURACY_ADAPTIVE, EULER_FACTOR, I_SPH_ALGORITHM              &
     &  , I_SPH_MODE                                                    &
!                       Precalculated angular arrays
     &  , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                      &
!                     Treatment of Scattering
     &  , I_SCATTER_METHOD                                              &
!                     Options for Solver
     &  , I_SOLVER                                                      &
!                     Gaseous Properties
     &  , K_GAS_ABS                                                     &
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
     &  , K_GREY_TOT_CLR, K_EXT_SCAT_CLR, PHASE_FNC_CLR                 &
     &  , FORWARD_SCATTER_CLR, PHASE_FNC_SOLAR_CLR                      &
     &  , K_GREY_TOT, K_EXT_SCAT                                        &
     &  , PHASE_FNC, FORWARD_SCATTER, PHASE_FNC_SOLAR                   &
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
!                       Calculated Radiances
     &  , RADIANCE                                                      &
!                       Calculated mean radiances
     &  , J_RADIANCE                                                    &
!                     Flags for Clear-sky Calculation
     &  , L_CLEAR, I_SOLVER_CLEAR                                       &
!                     Clear-sky Fluxes Calculated
     &  , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                           &
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
! Include Header Files
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
#include "angular_integration_pcf3z.h"
#include "surface_spec_pcf3z.h"
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
     &    I_ANGULAR_INTEGRATION                                         &
!           Angular integration scheme
     &  , I_2STREAM                                                     &
!           Two-stream scheme
     &  , N_ORDER_GAUSS                                                 &
!           Order of Gaussian integration
     &  , N_ORDER_PHASE                                                 &
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
!           Mode in which teh spherical harmonic solver is being used
     &  , I_SPH_ALGORITHM
!           Algorithm used for spherical harmonic calculation
      LOGICAL, INTENT(IN) ::                                            &
     &    L_RESCALE
!           Rescale optical properties
      REAL  (Real64) ::                                                 &
     &    CG_COEFF(ND_SPH_COEFF)                                        &
!           Clebsch-Gordan coefficients
     &  , UPLM_ZERO(ND_SPH_COEFF)                                       &
!           Values of spherical harmonics at polar angles pi/2
     &  , UPLM_SOL(ND_RADIANCE_PROFILE, ND_SPH_COEFF)
!           Values of spherical harmonics in the solar direction
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ACCURACY_ADAPTIVE                                             &
!           Accuracy for adaptive truncation
     &  , EULER_FACTOR
!           Factor applied to the last term of an alternating series
!
!                     Treatment of scattering
      INTEGER, INTENT(IN) ::                                            &
     &    I_SCATTER_METHOD
!
!                     Options for solver
      INTEGER, INTENT(IN) ::                                            &
     &    I_SOLVER
!           Solver used
!
!                     Gaseous properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    K_GAS_ABS(ND_PROFILE, ND_LAYER)
!           Gaseous absorptive extinctions
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
     &    I_DIRECT(ND_RADIANCE_PROFILE, 0: ND_LAYER)
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
     &    K_GREY_TOT_CLR(ND_PROFILE, ND_LAYER_CLR)                      &
!           Clear-sky absorptive extinction
     &  , K_EXT_SCAT_CLR(ND_PROFILE, ND_LAYER_CLR)                      &
!           Clear-sky scattering extinction
     &  , PHASE_FNC_CLR(ND_PROFILE, ND_LAYER_CLR, ND_MAX_ORDER)         &
!           Clear-sky phase function
     &  , FORWARD_SCATTER_CLR(ND_PROFILE, ND_LAYER_CLR)                 &
!           Clear-sky forward scattering
     &  , PHASE_FNC_SOLAR_CLR(ND_RADIANCE_PROFILE, ND_LAYER_CLR         &
     &      , ND_DIRECTION)
!           Clear-sky solar phase fuction in viewing directions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    K_GREY_TOT(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)     &
!           Absorptive extinction
     &  , K_EXT_SCAT(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)     &
!           Scattering extinction
     &  , PHASE_FNC(ND_PROFILE, ID_CT: ND_LAYER                         &
     &      , ND_MAX_ORDER, 0: ND_CLOUD_TYPE)                           &
!           Phase functions
     &  , FORWARD_SCATTER(ND_PROFILE, ID_CT: ND_LAYER                   &
     &      , 0: ND_CLOUD_TYPE)                                         &
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
!                     Calculated mean radiances
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    J_RADIANCE(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Mean radiances
!
!                     Flags for clear-sky calculations
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLEAR
!           Calculate clear-sky properties
      INTEGER, INTENT(IN) ::                                            &
     &    I_SOLVER_CLEAR
!           Clear solver used
!
!                     Clear-sky fluxes calculated
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT_CLEAR(ND_FLUX_PROFILE, 0: ND_LAYER)               &
!           Clear-sky direct flux
     &  , FLUX_TOTAL_CLEAR(ND_FLUX_PROFILE, 2*ND_LAYER+2)
!           Clear-sky total flux
!
!
!
!     Local variables.
      INTEGER                                                           &
     &    K                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , I
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
     &    TAU_CLR_F(:, :)
!           Clear-sky optical depths for the whole column
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    SINGLE_SCATTERING_ALL, RESCALE_TAU_OMEGA                      &
     &  , MONOCHROMATIC_RADIANCE_TSEQ                                   &
     &  , MONOCHROMATIC_RADIANCE_SPH                                    &
     &  , GAUSS_ANGLE
!
!
!
!     Calculate the optical depths and albedos of single scattering.
!     The phase function is not involved here as that is constant
!     across the band, whereas these parameters vary with the gaseous
!     absorption.
!
! DEPENDS ON: single_scattering_all
      CALL SINGLE_SCATTERING_ALL(I_SCATTER_METHOD                       &
!                     Atmospheric properties
     &  , N_PROFILE, N_LAYER, D_MASS                                    &
!                     Cloudy properties
     &  , L_CLOUD, N_CLOUD_TOP, N_CLOUD_TYPE                            &
!                     Optical properties
     &  , K_GREY_TOT_CLR, K_EXT_SCAT_CLR                                &
     &  , K_GREY_TOT, K_EXT_SCAT                                        &
     &  , K_GAS_ABS                                                     &
!                     Single scattering properties
     &  , TAU_CLR, OMEGA_CLR                                            &
     &  , TAU, OMEGA                                                    &
!                     Dimensions of arrays
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_CLOUD_TYPE      &
     &  )
!
!
!
      IF ( (I_ANGULAR_INTEGRATION == IP_TWO_STREAM).OR.                 &
     &     (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC) ) THEN
!
!       Rescale the optical depth and albedo of single scattering.
!
        IF (L_RESCALE) THEN
!
! DEPENDS ON: rescale_tau_omega
          CALL RESCALE_TAU_OMEGA(N_PROFILE, 1, N_CLOUD_TOP-1            &
     &      , TAU_CLR, OMEGA_CLR, FORWARD_SCATTER_CLR                   &
     &      , ND_PROFILE, ND_LAYER_CLR, 1                               &
     &      )
! DEPENDS ON: rescale_tau_omega
          CALL RESCALE_TAU_OMEGA(N_PROFILE, N_CLOUD_TOP, N_LAYER        &
     &      , TAU, OMEGA, FORWARD_SCATTER                               &
     &      , ND_PROFILE, ND_LAYER, ID_CT                               &
     &      )
!
          IF (L_CLOUD) THEN
!
            DO K=1, N_CLOUD_TYPE
! DEPENDS ON: rescale_tau_omega
              CALL RESCALE_TAU_OMEGA(N_PROFILE, N_CLOUD_TOP             &
     &          , N_LAYER                                               &
     &          , TAU(1, ID_CT, K), OMEGA(1, ID_CT, K)                  &
     &          , FORWARD_SCATTER(1, ID_CT, K)                          &
     &          , ND_PROFILE, ND_LAYER, ID_CT                           &
     &          )
            ENDDO
!
          ENDIF
!
        ENDIF
!
      ENDIF
!
!
!     Now divide the algorithmic path depending on the option
!     for angular integration.
!
      IF (I_ANGULAR_INTEGRATION == IP_TWO_STREAM) THEN
!
!       The standard two-stream approximations.
! DEPENDS ON: monochromatic_radiance_tseq
        CALL MONOCHROMATIC_RADIANCE_TSEQ(IERR                           &
!                       Atmospheric Propertries
     &    , N_PROFILE, N_LAYER                                          &
!                       Options for Solver
     &    , I_2STREAM, I_SOLVER, I_SCATTER_METHOD                       &
!                       Optical Properties
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
!                       Spectral Region
     &    , ISOLIR                                                      &
!                       Infra-red Properties
     &    , DIFF_PLANCK, L_IR_SOURCE_QUAD, DIFF_PLANCK_2                &
!                       Conditions at TOA
     &    , ZEN_0, FLUX_INC_DIRECT, FLUX_INC_DOWN                       &
!                       Surface Properties
     &    , D_PLANCK_FLUX_SURFACE                                       &
     &    , RHO_ALB                                                     &
!                       Optical Properties
     &    , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR                           &
     &    , TAU, OMEGA, PHASE_FNC                                       &
!                       Cloudy Properties
     &    , I_CLOUD                                                     &
!                       Cloud Geometry
     &    , N_CLOUD_TOP                                                 &
     &    , N_CLOUD_TYPE, FRAC_CLOUD                                    &
     &    , N_REGION, K_CLR, I_REGION_CLOUD, FRAC_REGION                &
     &    , W_FREE, W_CLOUD, CLOUD_OVERLAP                              &
     &    , N_COLUMN_SLV, LIST_COLUMN_SLV                               &
     &    , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN                   &
!                       Fluxes Calculated
     &    , FLUX_DIRECT, FLUX_TOTAL                                     &
!                       Flags for Clear-sky Calculation
     &    , L_CLEAR, I_SOLVER_CLEAR                                     &
!                       Clear-sky Fluxes Calculated
     &    , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                         &
!                       Dimensions of Arrays
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN        &
     &    , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF                  &
     &    , ND_SOURCE_COEFF, ND_MAX_ORDER                               &
     &    )
!

      ELSE IF (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC) THEN
!
!       The spherical harmonic option:
! DEPENDS ON: monochromatic_radiance_sph
        CALL MONOCHROMATIC_RADIANCE_SPH(IERR                            &
!                       Atmospheric Propertries
     &    , N_PROFILE, N_LAYER, D_MASS                                  &
!                       Angular Integration
     &    , N_ORDER_PHASE, MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC &
     &    , ACCURACY_ADAPTIVE, EULER_FACTOR, I_SPH_ALGORITHM            &
     &    , I_SPH_MODE, L_RESCALE                                       &
!                       Precalculated angular arrays
     &    , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                    &
!                       Options for Equivalent Extinction
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
!                       Spectral Region
     &    , ISOLIR                                                      &
!                       Infra-red Properties
     &    , DIFF_PLANCK, L_IR_SOURCE_QUAD, DIFF_PLANCK_2                &
!                       Conditions at TOA
     &    , ZEN_0, FLUX_INC_DIRECT, FLUX_INC_DOWN                       &
     &    , I_DIRECT                                                    &
!                       Surface Properties
     &    , D_PLANCK_FLUX_SURFACE                                       &
     &    , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB                    &
     &    , F_BRDF, BRDF_SOL, BRDF_HEMI                                 &
!                       Clear-sky Optical Properties
     &    , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR, FORWARD_SCATTER_CLR      &
     &    , PHASE_FNC_SOLAR_CLR                                         &
     &    , TAU, OMEGA, PHASE_FNC, FORWARD_SCATTER                      &
     &    , PHASE_FNC_SOLAR                                             &
!                       Cloudy Properties
     &    , L_CLOUD, I_CLOUD                                            &
!                       Cloud Geometry
     &    , N_CLOUD_TOP                                                 &
     &    , N_CLOUD_TYPE, FRAC_CLOUD                                    &
     &    , N_REGION, K_CLR, I_REGION_CLOUD, FRAC_REGION                &
     &    , W_FREE, W_CLOUD, CLOUD_OVERLAP                              &
     &    , N_COLUMN_SLV, LIST_COLUMN_SLV                               &
     &    , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN                   &
!                       Levels for calculating radiances
     &    , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                &
!                       Viewing geometry
     &    , N_DIRECTION, DIRECTION                                      &
!                       Calculated Fluxes
     &    , FLUX_DIRECT, FLUX_TOTAL                                     &
!                       Calculated radiances
     &    , RADIANCE                                                    &
!                       Calculated mean radiances
     &    , J_RADIANCE                                                  &
!                       Dimensions of Arrays
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN        &
     &    , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE          &
     &    , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF                  &
     &    , ND_MAX_ORDER, ND_SPH_COEFF                                  &
     &    , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC, ND_VIEWING_LEVEL          &
     &    , ND_DIRECTION, ND_SOURCE_COEFF                               &
     &    )
!
      ELSE IF (I_ANGULAR_INTEGRATION == IP_IR_GAUSS) THEN
!
!       Full angular resolution using Gaussian integration.
!
        ALLOCATE(TAU_CLR_F(ND_PROFILE, ND_LAYER))
        DO I=1, N_CLOUD_TOP-1
          DO L=1, N_PROFILE
            TAU_CLR_F(L, I)=TAU_CLR(L, I)
          ENDDO
        ENDDO
        DO I=N_CLOUD_TOP, N_LAYER
          DO L=1, N_PROFILE
            TAU_CLR_F(L, I)=TAU_CLR(L, I)
          ENDDO
        ENDDO
!
! DEPENDS ON: gauss_angle
        CALL GAUSS_ANGLE(N_PROFILE, N_LAYER                             &
     &    , N_ORDER_GAUSS                                               &
     &    , TAU_CLR_F                                                   &
     &    , FLUX_INC_DOWN                                               &
     &    , DIFF_PLANCK, D_PLANCK_FLUX_SURFACE                          &
     &    , RHO_ALB(1, IP_SURF_ALB_DIFF)                                &
     &    , FLUX_TOTAL                                                  &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                             &
     &    , ND_PROFILE, ND_LAYER                                        &
     &    )
!
        DEALLOCATE(TAU_CLR_F)
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE MONOCHROMATIC_RADIANCE
#endif
