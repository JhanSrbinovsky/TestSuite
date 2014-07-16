#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the fluxes assuming random overlap.
!
! Method:
!       Monochromatic calculations are performed for each
!       combination of ESFT terms and the results are summed.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVE_BAND_RANDOM_OVERLAP(IERR                         &
!                     Atmospheric Column
     &  , N_PROFILE, N_LAYER, I_TOP, P, T, D_MASS                       &
!                     Angular Integration
     &  , I_ANGULAR_INTEGRATION, I_2STREAM                              &
     &  , N_ORDER_PHASE, L_RESCALE, N_ORDER_GAUSS                       &
     &  , MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC                  &
     &  , ACCURACY_ADAPTIVE, EULER_FACTOR                               &
     &  , I_SPH_ALGORITHM, I_SPH_MODE                                   &
!                     Precalculated angular arrays
     &  , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                      &
!                     Treatment of Scattering
     &  , I_SCATTER_METHOD                                              &
!                     Options for solver
     &  , I_SOLVER                                                      &
!                     Gaseous Properties
     &  , I_BAND, N_GAS                                                 &
     &  , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC          &
     &  , K_ESFT, W_ESFT, SCALE_VECTOR                                  &
     &  , P_REFERENCE, T_REFERENCE                                      &
     &  , GAS_MIX_RATIO, GAS_FRAC_RESCALED                              &
     &  , L_DOPPLER, DOPPLER_CORRECTION                                 &
!                     Spectral Region
     &  , ISOLIR                                                        &
!                     Solar Properties
     &  , ZEN_0, SOLAR_IRRAD                                            &
!                     Infra-red Properties
     &  , PLANCK_FLUX_TOP, PLANCK_FLUX_BOTTOM                           &
     &  , DIFF_PLANCK_BAND                                              &
     &  , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                          &
!                     Surface Properties
     &  , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB                      &
     &  , F_BRDF, BRDF_SOL, BRDF_HEMI                                   &
     &  , PLANCK_FLUX_GROUND                                            &
!                       Tiling of the surface
     &  , L_TILE, N_POINT_TILE, N_TILE, LIST_TILE, RHO_ALB_TILE         &
     &  , PLANCK_FLUX_TILE                                              &
!                     Optical Properties
     &  , K_GREY_TOT_CLR, K_EXT_SCAT_CLR, PHASE_FNC_CLR                 &
     &  , FORWARD_SCATTER_CLR, PHASE_FNC_SOLAR_CLR                      &
     &  , K_GREY_TOT, K_EXT_SCAT, PHASE_FNC                             &
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
!                       Viewing Geometry
     &  , N_DIRECTION, DIRECTION                                        &
!                       Weighting factor for the band
     &  , WEIGHT_BAND, L_INITIAL                                        &
!                     Fluxes Calculated
     &  , FLUX_DIRECT, FLUX_DOWN, FLUX_UP                               &
!                       Calculcated radiances
     &  , I_DIRECT, RADIANCE                                            &
!                       Calculcated rate of photolysis
     &  , PHOTOLYSIS                                                    &
!                     Flags for Clear-sky Fluxes
     &  , L_CLEAR, I_SOLVER_CLEAR                                       &
!                     Clear-sky Fluxes
     &  , FLUX_DIRECT_CLEAR, FLUX_DOWN_CLEAR, FLUX_UP_CLEAR             &
!                       Tiled Surface Fluxes
     &  , FLUX_UP_TILE, FLUX_UP_BLUE_TILE                               &
!                       Special Surface Fluxes
     &  , L_BLUE_FLUX_SURF, WEIGHT_BLUE                                 &
     &  , FLUX_DIRECT_BLUE_SURF                                         &
     &  , FLUX_DOWN_BLUE_SURF, FLUX_UP_BLUE_SURF                        &
!                       Dimensions
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN          &
     &  , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE            &
     &  , ND_BAND, ND_SPECIES                                           &
     &  , ND_ESFT_TERM, ND_SCALE_VARIABLE                               &
     &  , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF                    &
     &  , ND_MAX_ORDER, ND_SPH_COEFF                                    &
     &  , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC, ND_VIEWING_LEVEL            &
     &  , ND_DIRECTION, ND_SOURCE_COEFF                                 &
     &  , ND_POINT_TILE, ND_TILE                                        &
     &  )
!
!
      IMPLICIT NONE
!
! Include header files
#include "c_kinds.h"
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ND_LAYER                                                      &
!           Maximum number of layers
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for totally clear layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_FLUX_PROFILE                                               &
!           Size allocated for profiles in arrays of fluxes
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for profiles in arrays of radiances
     &  , ND_J_PROFILE                                                  &
!           Size allocated for profiles in arrays of mean radiances
     &  , ND_BAND                                                       &
!           Maximum number of spectral bands
     &  , ND_SPECIES                                                    &
!           Maximum number of species
     &  , ND_ESFT_TERM                                                  &
!           Maximum number of ESFT terms
     &  , ND_SCALE_VARIABLE                                             &
!           Maximum number of scale variables
     &  , ND_COLUMN                                                     &
!           Number of columns per point
     &  , ND_CLOUD_TYPE                                                 &
!           Size allocated for cloud types
     &  , ND_REGION                                                     &
!           Size allocated for cloudy regions
     &  , ND_OVERLAP_COEFF                                              &
!           Size allocated for cloudy overlap coefficients
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_SPH_COEFF                                                  &
!           Size allocated for spherical harmonic coefficients
     &  , ND_BRDF_BASIS_FNC                                             &
!           Size allowed for BRDF basis functions
     &  , ND_BRDF_TRUNC                                                 &
!           Size allowed for orders of BRDFs
     &  , ND_VIEWING_LEVEL                                              &
!           Size allocated for levels where radiances are calculated
     &  , ND_DIRECTION                                                  &
!           Size allocated for viewing directions
     &  , ND_SOURCE_COEFF                                               &
!           Size allocated for source coefficients
     &  , ND_POINT_TILE                                                 &
!           Size allocated for points where the surface is tiled
     &  , ND_TILE
!           Size allocated for surface tiles
!
!     Include header files.
#include "angular_integration_pcf3z.h"
#include "surface_spec_pcf3z.h"
#include "spectral_region_pcf3z.h"
#include "esft_scale_pcf3z.h"
#include "error_pcf3z.h"
!
!
!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
!
!                     Atmospheric column
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , I_TOP
!           Top of vertical grid
      REAL  (Real64), INTENT(IN) ::                                     &
     &    D_MASS(ND_PROFILE, ND_LAYER)                                  &
!           Mass thickness of each layer
     &  , P(ND_PROFILE, ND_LAYER)                                       &
!           Pressure
     &  , T(ND_PROFILE, ND_LAYER)
!           Temperature
!
!                     Angular integration
      INTEGER, INTENT(IN) ::                                            &
     &    I_ANGULAR_INTEGRATION                                         &
!           Angular integration scheme
     &  , I_2STREAM                                                     &
!           Two-stream scheme
     &  , N_ORDER_PHASE                                                 &
!           Maximum order of terms in the phase function used in
!           the direct calculation of spherical harmonics
     &  , N_ORDER_GAUSS                                                 &
!           Order of gaussian integration
     &  , MS_MIN                                                        &
!           Lowest azimuthal order used
     &  , MS_MAX                                                        &
!           Highest azimuthal order used
     &  , I_TRUNCATION                                                  &
!           Type of spherical truncation used
     &  , IA_SPH_MM(0: ND_MAX_ORDER)                                    &
!           Address of spherical coefficient of (m, m) for each m
     &  , LS_LOCAL_TRUNC(0: ND_MAX_ORDER)                               &
!           Orders of truncation at each azimuthal order
     &  , I_SPH_MODE                                                    &
!           Mode in which the spherical solver is to be used
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
     &  , UPLM_SOL(ND_RADIANCE_PROFILE, ND_SPH_COEFF)                   &
!           Values of spherical harmonics in the solar direction
     &  , ACCURACY_ADAPTIVE                                             &
!           Accuracy for adaptive truncation
     &  , EULER_FACTOR
!           Factor applied to the last term of an alternating series
      REAL  (Real64), INTENT(IN) ::                                     &
     &    WEIGHT_BAND
!           Weighting factor for the current band
      LOGICAL, INTENT(INOUT) ::                                         &
     &    L_INITIAL
!           Flag to initialize diagnostics
!
!                     Treatment of scattering
      INTEGER, INTENT(IN) ::                                            &
     &    I_SCATTER_METHOD
!           Method of treating scattering
!
!                     Options for solver
      INTEGER, INTENT(IN) ::                                            &
     &    I_SOLVER
!           Solver used
!
!                     Gaseous properties
      INTEGER, INTENT(IN) ::                                            &
     &    I_BAND                                                        &
!           Band being considered
     &  , N_GAS                                                         &
!           Number of gases in band
     &  , INDEX_ABSORB(ND_SPECIES, ND_BAND)                             &
!           List of absorbers in bands
     &  , I_BAND_ESFT(ND_BAND, ND_SPECIES)                              &
!           Number of terms in band
     &  , I_SCALE_ESFT(ND_BAND, ND_SPECIES)                             &
!           Type of ESFT scaling
     &  , I_SCALE_FNC(ND_BAND, ND_SPECIES)
!           Type of scaling function
      LOGICAL, INTENT(IN) ::                                            &
     &    L_DOPPLER(ND_SPECIES)
!           Doppler broadening included
      REAL  (Real64), INTENT(IN) ::                                     &
     &    K_ESFT(ND_ESFT_TERM, ND_BAND, ND_SPECIES)                     &
!           Exponential ESFT terms
     &  , W_ESFT(ND_ESFT_TERM, ND_BAND, ND_SPECIES)                     &
!           Weights for ESFT
     &  , SCALE_VECTOR(ND_SCALE_VARIABLE, ND_ESFT_TERM, ND_BAND         &
     &      , ND_SPECIES)                                               &
!           Absorber scaling parameters
     &  , P_REFERENCE(ND_SPECIES, ND_BAND)                              &
!           Reference scaling pressure
     &  , T_REFERENCE(ND_SPECIES, ND_BAND)                              &
!           Reference scaling temperature
     &  , GAS_MIX_RATIO(ND_PROFILE, ND_LAYER, ND_SPECIES)               &
!           Gas mass mixing ratios
     &  , DOPPLER_CORRECTION(ND_SPECIES)
!           Doppler broadening terms
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    GAS_FRAC_RESCALED(ND_PROFILE, ND_LAYER, ND_SPECIES)
!           Rescaled gas mass fractions
!
!                     Spectral region
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Spectral region
!
!                     Solar properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ZEN_0(ND_PROFILE)                                             &
!           Secants (two-stream) or cosines (spherical harmonics)
!           of the solar zenith angle
     &  , SOLAR_IRRAD(ND_PROFILE)
!           Incident solar irradiance in band
!
!                     Infra-red properties
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Use a quadratic source function
      REAL  (Real64), INTENT(IN) ::                                     &
     &    PLANCK_FLUX_TOP(ND_PROFILE)                                   &
!           Planckian flux at top
     &  , PLANCK_FLUX_BOTTOM(ND_PROFILE)                                &
!           Planckian flux at bottom
     &  , DIFF_PLANCK_BAND(ND_PROFILE, ND_LAYER)                        &
!           Thermal source function
     &  , DIFF_PLANCK_BAND_2(ND_PROFILE, ND_LAYER)
!           2x2nd difference of Planckian in band
!
!                     Surface properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    PLANCK_FLUX_GROUND(ND_PROFILE)
!           Planckian flux at the surface temperature
      INTEGER, INTENT(IN) ::                                            &
     &    LS_BRDF_TRUNC                                                 &
!           Order of truncation of BRDFs
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
!     Variables related to tiling of the surface
      LOGICAL, INTENT(IN) ::                                            &
     &    L_TILE
!           Logical to allow invoke options
      INTEGER, INTENT(IN) ::                                            &
     &    N_POINT_TILE                                                  &
!           Number of points to tile
     &  , N_TILE                                                        &
!           Number of tiles used
     &  , LIST_TILE(ND_POINT_TILE)
!           List of points with surface tiling
      REAL  (Real64), INTENT(IN) ::                                     &
     &    RHO_ALB_TILE(ND_POINT_TILE, ND_BRDF_BASIS_FNC, ND_TILE)       &
!           Weights for the basis functions of the BRDFs
!           at the tiled points
     &  , PLANCK_FLUX_TILE(ND_POINT_TILE, ND_TILE)
!           Local Planckian fluxes on surface tiles
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
!           Phase function
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
!           Cloud enabled
      INTEGER, INTENT(IN) ::                                            &
     &    I_CLOUD
!           Cloud scheme used
!
!                     Cloud geometry
      INTEGER, INTENT(IN) ::                                            &
     &    N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_CLOUD_TYPE                                                  &
!           Number of types of cloud
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
!           Fractions of types of clouds
     &  , W_FREE(ND_PROFILE, ID_CT: ND_LAYER)                           &
!           Clear-sky fraction
     &  , CLOUD_OVERLAP(ND_PROFILE, ID_CT-1: ND_LAYER                   &
     &      , ND_OVERLAP_COEFF)                                         &
!           Coefficients for transfer for energy at interfaces
     &  , AREA_COLUMN(ND_PROFILE, ND_COLUMN)                            &
!           Areas of columns
     &  , FRAC_REGION(ND_PROFILE, ID_CT: ND_LAYER, ND_REGION)
!           Fractions of total cloud occupied by each region
!
!
!
!                       Viewing Geometry
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION
!           Number of viewing directions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DIRECTION(ND_RADIANCE_PROFILE, ND_DIRECTION, 2)
!           Viewing directions
      INTEGER, INTENT(IN) ::                                            &
     &    N_VIEWING_LEVEL                                               &
!           Number of levels where radiances are calculated
     &  , I_RAD_LAYER(ND_VIEWING_LEVEL)
!           Layers in which radiances are calculated
      REAL  (Real64), INTENT(IN) ::                                     &
     &    FRAC_RAD_LAYER(ND_VIEWING_LEVEL)
!           Fractions below the tops of the layers
!
!                     Flags for clear-sky calculations
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLEAR
!           Calculate clear-sky properties
      INTEGER, INTENT(IN) ::                                            &
     &    I_SOLVER_CLEAR
!           Clear solver used
!
!                     Calculated Fluxes
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_FLUX_PROFILE, 0: ND_LAYER)                     &
!           Direct flux
     &  , FLUX_DOWN(ND_FLUX_PROFILE, 0: ND_LAYER)                       &
!           Total downward flux
     &  , FLUX_UP(ND_FLUX_PROFILE, 0: ND_LAYER)
!           Upward flux
!
!                       Calculated radiances
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    I_DIRECT(ND_RADIANCE_PROFILE, 0: ND_LAYER)                    &
!           Direct solar irradiance on levels
     &  , RADIANCE(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL                &
     &      , ND_DIRECTION)
!           Radiances
!
!                       Calculated mean radiances
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    PHOTOLYSIS(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Rates of photolysis
!
!                     Clear-sky fluxes calculated
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    FLUX_DIRECT_CLEAR(ND_FLUX_PROFILE, 0: ND_LAYER)               &
!           Clear-sky direct flux
     &  , FLUX_DOWN_CLEAR(ND_FLUX_PROFILE, 0: ND_LAYER)                 &
!           Clear-sky total downward flux in band
     &  , FLUX_UP_CLEAR(ND_FLUX_PROFILE, 0: ND_LAYER)                   &
!           Clear-sky upward flux
     &  , FLUX_UP_TILE(ND_POINT_TILE, ND_TILE)                          &
!           Upward fluxes at tiled surface points
     &  , FLUX_UP_BLUE_TILE(ND_POINT_TILE, ND_TILE)
!           Upward blue fluxes at tiled surface points
!
!                      Special Diagnostics:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_BLUE_FLUX_SURF
!           Flag to calculate the blue flux at the surface
      REAL  (Real64), INTENT(IN) ::                                     &
     &    WEIGHT_BLUE
!           Weights for blue fluxes in this band
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    FLUX_DIRECT_BLUE_SURF(ND_FLUX_PROFILE)                        &
!           Direct blue flux at the surface
     &  , FLUX_DOWN_BLUE_SURF(ND_FLUX_PROFILE)                          &
!           Total downward blue flux at the surface
     &  , FLUX_UP_BLUE_SURF(ND_FLUX_PROFILE)
!           Upward blue flux at the surface
!
!
!
!     Local variables.
      INTEGER                                                           &
     &     J                                                            &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , L
!           Loop variable
      INTEGER                                                           &
     &    I_GAS_BAND                                                    &
!           Index of active gas
     &  , I_GAS_POINTER(ND_SPECIES)                                     &
!           Pointer array for monochromatic ESFTs
     &  , I_ESFT_POINTER(ND_SPECIES)                                    &
!           Pointer to ESFT for gas
     &  , I_CHANGE                                                      &
!           Position of ESFT term to be altered
     &  , INDEX_CHANGE                                                  &
!           Index of term to be altered
     &  , INDEX_LAST                                                    &
!           Index of last gas in band
     &  , IEX
!           Index of ESFT term
      REAL  (Real64) ::                                                 &
     &    K_ESFT_MONO(ND_SPECIES)                                       &
!           ESFT monochromatic exponents
     &  , K_GAS_ABS(ND_PROFILE, ND_LAYER)                               &
!           Gaseous absorption
     &  , D_PLANCK_FLUX_SURFACE(ND_PROFILE)                             &
!           Difference in Planckian fluxes between the surface and
!           the air
     &  , FLUX_INC_DIRECT(ND_PROFILE)                                   &
!           Incident direct flux
     &  , FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident downward flux
     &  , PRODUCT_WEIGHT                                                &
!           Product of ESFT weights
     &  , DUMMY_KE(ND_PROFILE, ND_LAYER)
!           Dummy array (not used)
!
!     Monochromatic incrementing radiances:
      REAL  (Real64) ::                                                 &
     &    FLUX_DIRECT_PART(ND_FLUX_PROFILE, 0: ND_LAYER)                &
!           Partial direct flux
     &  , FLUX_TOTAL_PART(ND_FLUX_PROFILE, 2*ND_LAYER+2)                &
!           Partial total flux
     &  , FLUX_DIRECT_CLEAR_PART(ND_FLUX_PROFILE, 0: ND_LAYER)          &
!           Partial clear-sky direct flux
     &  , FLUX_TOTAL_CLEAR_PART(ND_FLUX_PROFILE, 2*ND_LAYER+2)
!           Partial clear-sky total flux
      REAL  (Real64) ::                                                 &
     &    I_DIRECT_PART(ND_RADIANCE_PROFILE, 0: ND_LAYER)               &
!           Partial solar irradiances
     &  , RADIANCE_PART(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL           &
     &      , ND_DIRECTION)
!           Partial radiances
      REAL  (Real64) ::                                                 &
     &    PHOTOLYSIS_PART(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Partial rates of photolysis
      REAL  (Real64) ::                                                 &
     &    WEIGHT_INCR                                                   &
!           Weight applied to increments
     &  , WEIGHT_BLUE_INCR
!           Weight applied to blue increments
!
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    SCALE_ABSORB, GAS_OPTICAL_PROPERTIES                          &
     &  , MONOCHROMATIC_RADIANCE, AUGMENT_RADIANCE
!
!
!
!     Set the number of active gases and initialize the pointers.
      DO K=1, N_GAS
        I_GAS_POINTER(K)=INDEX_ABSORB(K, I_BAND)
        I_ESFT_POINTER(INDEX_ABSORB(K, I_BAND))=1
      ENDDO
      INDEX_LAST=INDEX_ABSORB(N_GAS, I_BAND)
!
!     Perform the initial rescaling of the gases other than the last.
!     Note: we rescale amounts as required. It would be more
!     efficient to save the rescaled amounts, but the storage
!     needed would become excessive for a multicolumn code. In a
!     single code the overhead would be less significant.
      DO K=1, N_GAS-1
        I_GAS_BAND=I_GAS_POINTER(K)
!       Initialize the monochromatic absorption coefficients.
        K_ESFT_MONO(I_GAS_BAND)                                         &
     &    =K_ESFT(1, I_BAND, I_GAS_BAND)
        IF (I_SCALE_ESFT(I_BAND, I_GAS_BAND) == IP_SCALE_TERM) THEN
! DEPENDS ON: scale_absorb
          CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER                    &
     &      , GAS_MIX_RATIO(1, 1, I_GAS_BAND), P, T                     &
     &      , I_TOP                                                     &
     &      , GAS_FRAC_RESCALED(1, 1, I_GAS_BAND)                       &
     &      , I_SCALE_FNC(I_BAND, I_GAS_BAND)                           &
     &      , P_REFERENCE(I_GAS_BAND, I_BAND)                           &
     &      , T_REFERENCE(I_GAS_BAND, I_BAND)                           &
     &      , SCALE_VECTOR(1, 1, I_BAND, I_GAS_BAND)                    &
     &      , L_DOPPLER(I_GAS_BAND), DOPPLER_CORRECTION(I_GAS_BAND)     &
     &      , ND_PROFILE, ND_LAYER                                      &
     &      , ND_SCALE_VARIABLE                                         &
     &      )
          IF (IERR /= I_NORMAL) RETURN
        ENDIF
      ENDDO
!
!     Loop through the terms for the first absorber.
2000  I_ESFT_POINTER(INDEX_LAST)=0
      DO K=1, I_BAND_ESFT(I_BAND, INDEX_LAST)
        I_ESFT_POINTER(INDEX_LAST)                                      &
     &    =I_ESFT_POINTER(INDEX_LAST)+1
!
!       Set the ESFT coefficient and perform rescaling for the
!       last gas.
        IEX=I_ESFT_POINTER(INDEX_LAST)
        K_ESFT_MONO(INDEX_LAST)                                         &
     &    =K_ESFT(IEX, I_BAND, INDEX_LAST)
        IF (I_SCALE_ESFT(I_BAND, INDEX_LAST) == IP_SCALE_TERM) THEN
! DEPENDS ON: scale_absorb
          CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER                    &
     &      , GAS_MIX_RATIO(1, 1, INDEX_LAST), P, T                     &
     &      , I_TOP                                                     &
     &      , GAS_FRAC_RESCALED(1, 1, INDEX_LAST)                       &
     &      , I_SCALE_FNC(I_BAND, INDEX_LAST)                           &
     &      , P_REFERENCE(INDEX_LAST, I_BAND)                           &
     &      , T_REFERENCE(INDEX_LAST, I_BAND)                           &
     &      , SCALE_VECTOR(1, IEX, I_BAND, INDEX_LAST)                  &
     &      , L_DOPPLER(INDEX_LAST)                                     &
     &      , DOPPLER_CORRECTION(INDEX_LAST)                            &
     &      , ND_PROFILE, ND_LAYER                                      &
     &      , ND_SCALE_VARIABLE                                         &
     &      )
          IF (IERR /= I_NORMAL) RETURN
        ENDIF
!
!       Set the appropriate source terms for the two-stream
!       equations.
!       The product of the ESFT weights can be precalculated
!       for speed.
        PRODUCT_WEIGHT=1.0E+00_Real64
        DO J=1, N_GAS
          I_GAS_BAND=I_GAS_POINTER(J)
          IEX=I_ESFT_POINTER(I_GAS_BAND)
          PRODUCT_WEIGHT=PRODUCT_WEIGHT                                 &
     &      *W_ESFT(IEX, I_BAND, I_GAS_BAND)
        ENDDO
!
        IF ( (I_ANGULAR_INTEGRATION == IP_TWO_STREAM).OR.               &
     &       (I_ANGULAR_INTEGRATION == IP_IR_GAUSS) ) THEN
!
          IF (ISOLIR == IP_SOLAR) THEN
!
!           Solar region.
            DO L=1, N_PROFILE
              D_PLANCK_FLUX_SURFACE(L)=0.0E+00_Real64
              FLUX_INC_DOWN(L)=SOLAR_IRRAD(L)/ZEN_0(L)
              FLUX_INC_DIRECT(L)=SOLAR_IRRAD(L)/ZEN_0(L)
            ENDDO
!
          ELSEIF (ISOLIR == IP_INFRA_RED) THEN
!           Infra-red region.
!
            DO L=1, N_PROFILE
              FLUX_INC_DIRECT(L)=0.0E+00_Real64
              FLUX_DIRECT_PART(L, N_LAYER)=0.0E+00_Real64
              FLUX_INC_DOWN(L)=-PLANCK_FLUX_TOP(L)
              D_PLANCK_FLUX_SURFACE(L)                                  &
     &          =PLANCK_FLUX_GROUND(L)-PLANCK_FLUX_BOTTOM(L)
            ENDDO
            IF (L_CLEAR) THEN
              DO L=1, N_PROFILE
                FLUX_DIRECT_CLEAR_PART(L, N_LAYER)=0.0E+00_Real64
              ENDDO
            ENDIF
!
          ENDIF
!
        ELSE IF (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC) THEN
!
          IF (ISOLIR == IP_SOLAR) THEN
            DO L=1, N_PROFILE
              I_DIRECT_PART(L, 0)=SOLAR_IRRAD(L)
              FLUX_INC_DOWN(L)=0.0E+00_Real64
            ENDDO
          ELSE
            DO L=1, N_PROFILE
              FLUX_INC_DOWN(L)=-PLANCK_FLUX_TOP(L)
              D_PLANCK_FLUX_SURFACE(L)                                  &
     &          =PLANCK_FLUX_GROUND(L)-PLANCK_FLUX_BOTTOM(L)
            ENDDO
          ENDIF
!
        ENDIF
!
! DEPENDS ON: gas_optical_properties
        CALL GAS_OPTICAL_PROPERTIES(N_PROFILE, N_LAYER                  &
     &    , N_GAS, I_GAS_POINTER, K_ESFT_MONO                           &
     &    , GAS_FRAC_RESCALED                                           &
     &    , K_GAS_ABS                                                   &
     &    , ND_PROFILE, ND_LAYER, ND_SPECIES                            &
     &    )
!
!
! DEPENDS ON: monochromatic_radiance
        CALL MONOCHROMATIC_RADIANCE(IERR                                &
!                     Atmospheric properties
     &    , N_PROFILE, N_LAYER, D_MASS                                  &
!                     Angular integration
     &    , I_ANGULAR_INTEGRATION, I_2STREAM                            &
     &    , L_RESCALE, N_ORDER_GAUSS                                    &
     &    , N_ORDER_PHASE, MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC &
     &    , ACCURACY_ADAPTIVE, EULER_FACTOR                             &
     &    , I_SPH_ALGORITHM, I_SPH_MODE                                 &
!                       Precalculated angular arrays
     &    , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                    &
!                     Treatment of scattering
     &    , I_SCATTER_METHOD                                            &
!                     Options for solver
     &    , I_SOLVER                                                    &
!                     Gaseous propreties
     &    , K_GAS_ABS                                                   &
!                     Options for equivalent extinction
     &    , .FALSE., DUMMY_KE                                           &
!                     Spectral region
     &    , ISOLIR                                                      &
!                     Infra-red properties
     &    , DIFF_PLANCK_BAND                                            &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                        &
!                     Conditions at TOA
     &    , ZEN_0, FLUX_INC_DIRECT, FLUX_INC_DOWN                       &
     &    , I_DIRECT_PART                                               &
!                     Surface properties
     &    , D_PLANCK_FLUX_SURFACE                                       &
     &    , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB                    &
     &    , F_BRDF, BRDF_SOL, BRDF_HEMI                                 &
!                     Optical properties
     &    , K_GREY_TOT_CLR, K_EXT_SCAT_CLR                              &
     &    , PHASE_FNC_CLR, FORWARD_SCATTER_CLR                          &
     &    , PHASE_FNC_SOLAR_CLR                                         &
     &    , K_GREY_TOT, K_EXT_SCAT, PHASE_FNC, FORWARD_SCATTER          &
     &    , PHASE_FNC_SOLAR                                             &
!                     Cloudy properties
     &    , L_CLOUD, I_CLOUD                                            &
!                     Cloud geometry
     &    , N_CLOUD_TOP                                                 &
     &    , N_CLOUD_TYPE, FRAC_CLOUD                                    &
     &    , N_REGION, K_CLR, I_REGION_CLOUD, FRAC_REGION                &
     &    , W_FREE, W_CLOUD, CLOUD_OVERLAP                              &
     &    , N_COLUMN_SLV, LIST_COLUMN_SLV                               &
     &    , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN                   &
!                       Levels for calculating radiances
     &    , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                &
!                       Viewing Geometry
     &    , N_DIRECTION, DIRECTION                                      &
!                     Calculated fluxes
     &    , FLUX_DIRECT_PART, FLUX_TOTAL_PART                           &
!                       Calculated radiances
     &    , RADIANCE_PART                                               &
!                       Calculated rate of photolysis
     &    , PHOTOLYSIS_PART                                             &
!                     Flags for clear-sky calculations
     &    , L_CLEAR, I_SOLVER_CLEAR                                     &
!                     Clear-sky fluxes calculated
     &    , FLUX_DIRECT_CLEAR_PART, FLUX_TOTAL_CLEAR_PART               &
!                     Dimensions of arrays
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN        &
     &    , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE          &
     &    , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF                  &
     &    , ND_MAX_ORDER, ND_SPH_COEFF                                  &
     &    , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC, ND_VIEWING_LEVEL          &
     &    , ND_DIRECTION, ND_SOURCE_COEFF                               &
     &    )
        IF (IERR /= I_NORMAL) RETURN
!
!       Increment the fluxes within the band.
        WEIGHT_INCR=WEIGHT_BAND*PRODUCT_WEIGHT
        IF (L_BLUE_FLUX_SURF)                                           &
     &    WEIGHT_BLUE_INCR=WEIGHT_BLUE*PRODUCT_WEIGHT
! DEPENDS ON: augment_radiance
        CALL AUGMENT_RADIANCE(N_PROFILE, N_LAYER                        &
     &    , I_ANGULAR_INTEGRATION, I_SPH_MODE                           &
     &    , N_VIEWING_LEVEL, N_DIRECTION                                &
     &    , ISOLIR, L_CLEAR, L_INITIAL, WEIGHT_INCR                     &
     &    , L_BLUE_FLUX_SURF, WEIGHT_BLUE_INCR                          &
!                       Actual radiances
     &    , FLUX_DIRECT, FLUX_DOWN, FLUX_UP                             &
     &    , FLUX_DIRECT_BLUE_SURF                                       &
     &    , FLUX_DOWN_BLUE_SURF, FLUX_UP_BLUE_SURF                      &
     &    , I_DIRECT, RADIANCE, PHOTOLYSIS                              &
     &    , FLUX_DIRECT_CLEAR, FLUX_DOWN_CLEAR, FLUX_UP_CLEAR           &
!                       Increments to radiances
     &    , FLUX_DIRECT_PART, FLUX_TOTAL_PART                           &
     &    , I_DIRECT_PART, RADIANCE_PART, PHOTOLYSIS_PART               &
     &    , FLUX_DIRECT_CLEAR_PART, FLUX_TOTAL_CLEAR_PART               &
!                       Dimensions
     &    , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE          &
     &    , ND_LAYER, ND_VIEWING_LEVEL, ND_DIRECTION                    &
     &    )
!
!       Add in the increments from surface tiles
        IF (L_TILE) THEN
! DEPENDS ON: augment_tiled_radiance
          CALL AUGMENT_TILED_RADIANCE(IERR                              &
     &      , N_POINT_TILE, N_TILE, LIST_TILE                           &
     &      , I_ANGULAR_INTEGRATION, ISOLIR, L_INITIAL                  &
     &      , WEIGHT_INCR, L_BLUE_FLUX_SURF, WEIGHT_BLUE_INCR           &
!                       Surface characteristics
     &      , RHO_ALB_TILE                                              &
!                       Actual radiances
     &      , FLUX_UP_TILE, FLUX_UP_BLUE_TILE                           &
!                       Increments to radiances
     &      , FLUX_DIRECT_PART(1, N_LAYER)                              &
     &      , FLUX_TOTAL_PART(1, 2*N_LAYER+2)                           &
     &      , PLANCK_FLUX_TILE, PLANCK_FLUX_BOTTOM                      &
!                       Dimensions
     &      , ND_FLUX_PROFILE, ND_POINT_TILE, ND_TILE                   &
     &      , ND_BRDF_BASIS_FNC                                         &
     &      )
          IF (IERR /= I_NORMAL) RETURN
        ENDIF
!
!       After the first call to these routines quantities should be
!       incremented rather than initialized, until the flag is reset.
        L_INITIAL=.FALSE.
!
      ENDDO
!
      IF (N_GAS >  1) THEN
!       Increment the ESFT pointers for the next pass through
!       the loop above. I_CHANGE is the ordinal of the gas,
!       the pointer of which is to be changed.
        I_CHANGE=N_GAS-1
2001    INDEX_CHANGE=INDEX_ABSORB(I_CHANGE, I_BAND)
        IF (I_BAND_ESFT(I_BAND, INDEX_CHANGE)                           &
     &     >  I_ESFT_POINTER(INDEX_CHANGE)) THEN
          I_ESFT_POINTER(INDEX_CHANGE)                                  &
     &      =I_ESFT_POINTER(INDEX_CHANGE)+1
!         Rescale the amount of this gas and advance the ESFT term.
          K_ESFT_MONO(INDEX_CHANGE)                                     &
     &      =K_ESFT(I_ESFT_POINTER(INDEX_CHANGE)                        &
     &      , I_BAND, INDEX_CHANGE)
          IF (I_SCALE_ESFT(I_BAND, INDEX_CHANGE) == IP_SCALE_TERM)      &
     &         THEN
! DEPENDS ON: scale_absorb
            CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER                  &
     &        , GAS_MIX_RATIO(1, 1, INDEX_CHANGE), P, T                 &
     &        , I_TOP                                                   &
     &        , GAS_FRAC_RESCALED(1, 1, INDEX_CHANGE)                   &
     &        , I_SCALE_FNC(I_BAND, INDEX_CHANGE)                       &
     &        , P_REFERENCE(INDEX_CHANGE, I_BAND)                       &
     &        , T_REFERENCE(INDEX_CHANGE, I_BAND)                       &
     &        , SCALE_VECTOR(1, I_ESFT_POINTER(INDEX_CHANGE)            &
     &        , I_BAND, INDEX_CHANGE)                                   &
     &        , L_DOPPLER(INDEX_CHANGE)                                 &
     &        , DOPPLER_CORRECTION(INDEX_CHANGE)                        &
     &        , ND_PROFILE, ND_LAYER                                    &
     &        , ND_SCALE_VARIABLE                                       &
     &        )
            IF (IERR /= I_NORMAL) RETURN
          ENDIF
          GOTO 2000
        ELSE IF (I_CHANGE >  1) THEN
!         All terms for this absorber have been done:
!         reset its pointer to 1 and move to the next absorber.
          I_ESFT_POINTER(INDEX_CHANGE)=1
          K_ESFT_MONO(INDEX_CHANGE)=K_ESFT(1, I_BAND, INDEX_CHANGE)
          I_CHANGE=I_CHANGE-1
          GOTO 2001
        ENDIF
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE SOLVE_BAND_RANDOM_OVERLAP
#endif
