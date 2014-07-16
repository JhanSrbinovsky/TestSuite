#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the fluxes within the band with one gas.
!
! Method:
!       Monochromatic calculations are performed for each ESFT term
!       and the results are summed.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVE_BAND_ONE_GAS(IERR                                &
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
!                     Options for Solver
     &  , I_SOLVER                                                      &
!                     Gaseous Properties
     &  , I_BAND, I_GAS                                                 &
     &  , I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC                        &
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
!                     Calculated Fluxes
     &  , FLUX_DIRECT, FLUX_DOWN, FLUX_UP                               &
!                       Calculated radiances
     &  , I_DIRECT, RADIANCE                                            &
!                       Calculated rate of photolysis
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
!                     Dimensions of Arrays
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
!           Size allocated for profiles
     &  , ND_LAYER                                                      &
!           Size allocated for layers
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
     &  , ND_COLUMN                                                     &
!           Size allocated for sub-columns per point
     &  , ND_BAND                                                       &
!           Size allocated for bands
     &  , ND_SPECIES                                                    &
!           Size allocated for species
     &  , ND_ESFT_TERM                                                  &
!           Size allocated for ESFT variables
     &  , ND_SCALE_VARIABLE                                             &
!           Size allocated for scaling variables
     &  , ND_CLOUD_TYPE                                                 &
!           Size allocated for cloud types
     &  , ND_REGION                                                     &
!           Size allocated for cloudy regions
     &  , ND_OVERLAP_COEFF                                              &
!           Size allocated for cloudy overlap coefficients
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_SPH_COEFF                                                  &
!           Size allocated for coefficients of spherical harmonics
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
#include "spectral_region_pcf3z.h"
#include "surface_spec_pcf3z.h"
#include "angular_integration_pcf3z.h"
#include "esft_scale_pcf3z.h"
#include "error_pcf3z.h"
!
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
     &    P(ND_PROFILE, ND_LAYER)                                       &
!           Pressure
     &  , T(ND_PROFILE, ND_LAYER)                                       &
!           Temperature
     &  , D_MASS(ND_PROFILE, ND_LAYER)
!           Mass thickness of each layer
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
!           Type of truncation used
     &  , IA_SPH_MM(0: ND_MAX_ORDER)                                    &
!           Address of spherical coefficient of (m, m) for each m
     &  , LS_LOCAL_TRUNC(0: ND_MAX_ORDER)                               &
!           Orders of truncation at each azimuthal order
     &  , I_SPH_MODE                                                    &
!           Mode in which the spherical harmonic code is used
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
     &  , I_GAS                                                         &
!           Gas being considered
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
      REAL  (Real64), INTENT(INOUT) ::                                  &
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
     &     ZEN_0(ND_PROFILE)                                            &
!           Secant (two-stream) or cosine (spherical harmonics)
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
!           Planckian flux at the top of the layer
     &  , PLANCK_FLUX_BOTTOM(ND_PROFILE)                                &
!           Planckian source at the bottom of the layer
     &  , DIFF_PLANCK_BAND(ND_PROFILE, ND_LAYER)                        &
!           Thermal source function
     &  , DIFF_PLANCK_BAND_2(ND_PROFILE, ND_LAYER)
!           Twice second difference of Planckian in band
!
!                     Surface properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    PLANCK_FLUX_GROUND(ND_PROFILE)
!           Thermal source function at ground
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
!                     Clear-sky optical properties
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
     &  , PHASE_FNC(ND_PROFILE, ID_CT: ND_LAYER, ND_MAX_ORDER           &
     &      , 0: ND_CLOUD_TYPE)                                         &
!           Phase function
     &  , FORWARD_SCATTER(ND_PROFILE, ID_CT: ND_LAYER                   &
     &      , 0: ND_CLOUD_TYPE)                                         &
!           Forward scattering
     &  , PHASE_FNC_SOLAR(ND_RADIANCE_PROFILE, ID_CT: ND_LAYER          &
     &      , ND_DIRECTION, 0: ND_CLOUD_TYPE)
!           Phase function for the solar beam in viewing directions
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
!           Top cloudy layer
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
     &    W_FREE(ND_PROFILE, ID_CT: ND_LAYER)                           &
!           Clear-sky fraction
     &  , W_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)                          &
!           Cloudy fraction
     &  , FRAC_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_CLOUD_TYPE)        &
!           Fractions of types of clouds
     &  , CLOUD_OVERLAP(ND_PROFILE, ID_CT-1: ND_LAYER, ND_OVERLAP_COEFF)&
!           Coefficients for transfer for energy at interfaces
     &  , AREA_COLUMN(ND_PROFILE, ND_COLUMN)                            &
!           Areas of columns
     &  , FRAC_REGION(ND_PROFILE, ID_CT: ND_LAYER, ND_REGION)
!           Fractions of total cloud occupied by each region
!
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
!                       Viewing Geometry
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION
!           Number of viewing directions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DIRECTION(ND_RADIANCE_PROFILE, ND_DIRECTION, 2)
!           Viewing directions
!
!                     Calculated fluxes
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
     &  , RADIANCE(ND_RADIANCE_PROFILE,  ND_VIEWING_LEVEL               &
     &      , ND_DIRECTION)
!           Radiances
!
!                       Calculated mean radiances
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    PHOTOLYSIS(ND_J_PROFILE,  ND_VIEWING_LEVEL)
!           Mean rate of photolysis
!
!                     Flags for clear-sky calculations
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLEAR
!           Calculate net clear-sky properties
      INTEGER, INTENT(IN) ::                                            &
     &    I_SOLVER_CLEAR
!           Clear solver used
!
!                     Clear-sky fluxes calculated
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT_CLEAR(ND_PROFILE, 0: ND_LAYER)                    &
!           Clear-sky direct flux
     &  , FLUX_DOWN_CLEAR(ND_PROFILE, 0: ND_LAYER)                      &
!           Clear-sky total downward flux
     &  , FLUX_UP_CLEAR(ND_PROFILE, 0: ND_LAYER)                        &
!           Clear-sky upward flux
     &  , FLUX_UP_TILE(ND_POINT_TILE, ND_TILE)                          &
!           Upward fluxes at tiled surface points
     &  , FLUX_UP_BLUE_TILE(ND_POINT_TILE, ND_TILE)
!           Upward blue fluxes at tiled surface points
!
!                     Special Diagnostics:
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
     &    L
!           Loop variable
      INTEGER                                                           &
     &    I_GAS_POINTER(ND_SPECIES)                                     &
!           Pointer array for monochromatic ESFTs
     &  , IEX
!           Index of ESFT term
      REAL  (Real64) ::                                                 &
     &    K_ESFT_MONO(ND_SPECIES)                                       &
!           ESFT monochromatic exponents
     &  , K_GAS_ABS(ND_PROFILE, ND_LAYER)                               &
!           Gaseous absorptive extinction
     &  , D_PLANCK_FLUX_SURFACE(ND_PROFILE)                             &
!           Ground source function
     &  , FLUX_INC_DIRECT(ND_PROFILE)                                   &
!           Incident direct flux
     &  , FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident downward flux
     &  , DUMMY_KE(ND_PROFILE, ND_LAYER)
!           Dummy array (not used)
!
!     Monochromatic incrementing radiances:
      REAL  (Real64) ::                                                 &
     &    FLUX_DIRECT_PART(ND_PROFILE, 0: ND_LAYER)                     &
!           Partial direct flux
     &  , FLUX_TOTAL_PART(ND_PROFILE, 2*ND_LAYER+2)                     &
!           Partial total flux
     &  , FLUX_DIRECT_CLEAR_PART(ND_PROFILE, 0: ND_LAYER)               &
!           Partial clear-sky direct flux
     &  , FLUX_TOTAL_CLEAR_PART(ND_PROFILE, 2*ND_LAYER+2)
!           Partial clear-sky total flux
      REAL  (Real64) ::                                                 &
     &    I_DIRECT_PART(ND_RADIANCE_PROFILE, 0: ND_LAYER)               &
!           Partial solar irradiances
     &  , RADIANCE_PART(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL           &
     &      , ND_DIRECTION)
!           Partial radiances
      REAL  (Real64) ::                                                 &
     &    PHOTOLYSIS_PART(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Partial rate of photolysis
      REAL  (Real64) ::                                                 &
     &    WEIGHT_INCR                                                   &
!           Weight applied to increments
     &  , WEIGHT_BLUE_INCR
!           Weight applied to blue increments
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    SCALE_ABSORB, GAS_OPTICAL_PROPERTIES                          &
     &  , MONOCHROMATIC_RADIANCE, AUGMENT_RADIANCE
!
!
!
!     The ESFT terms for the first gas in the band alone are used.
      I_GAS_POINTER(1)=I_GAS
      DO IEX=1, I_BAND_ESFT(I_BAND, I_GAS)
!
!       Rescale for each ESFT term if that is required.
        IF (I_SCALE_ESFT(I_BAND, I_GAS) == IP_SCALE_TERM) THEN
! DEPENDS ON: scale_absorb
          CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER                    &
     &      , GAS_MIX_RATIO(1, 1, I_GAS), P, T                          &
     &      , I_TOP                                                     &
     &      , GAS_FRAC_RESCALED(1, 1, I_GAS)                            &
     &      , I_SCALE_FNC(I_BAND, I_GAS)                                &
     &      , P_REFERENCE(I_GAS, I_BAND)                                &
     &      , T_REFERENCE(I_GAS, I_BAND)                                &
     &      , SCALE_VECTOR(1, IEX, I_BAND, I_GAS)                       &
     &      , L_DOPPLER(I_GAS), DOPPLER_CORRECTION(I_GAS)               &
     &      , ND_PROFILE, ND_LAYER                                      &
     &      , ND_SCALE_VARIABLE                                         &
     &      )
          IF (IERR /= I_NORMAL) RETURN
        ENDIF
!
!       Set the appropriate boundary terms for the total
!       upward and downward fluxes.
!
        IF ( (I_ANGULAR_INTEGRATION == IP_TWO_STREAM).OR.               &
     &       (I_ANGULAR_INTEGRATION == IP_IR_GAUSS) ) THEN
!
          IF (ISOLIR == IP_SOLAR) THEN
!           Solar region.
            DO L=1, N_PROFILE
              D_PLANCK_FLUX_SURFACE(L)=0.0E+00_Real64
              FLUX_INC_DOWN(L)=SOLAR_IRRAD(L)/ZEN_0(L)
              FLUX_INC_DIRECT(L)=SOLAR_IRRAD(L)/ZEN_0(L)
            ENDDO
          ELSEIF (ISOLIR == IP_INFRA_RED) THEN
!           Infra-red region.
            DO L=1, N_PROFILE
              FLUX_INC_DIRECT(L)=0.0E+00_Real64
              FLUX_INC_DOWN(L)=-PLANCK_FLUX_TOP(L)
              D_PLANCK_FLUX_SURFACE(L)                                  &
     &          =(1.0E+00_Real64-RHO_ALB(L, IP_SURF_ALB_DIFF))          &
     &          *(PLANCK_FLUX_GROUND(L)-PLANCK_FLUX_BOTTOM(L))
            ENDDO
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
!       Assign the monochromatic absorption coefficient.
        K_ESFT_MONO(I_GAS)=K_ESFT(IEX, I_BAND, I_GAS)
!
! DEPENDS ON: gas_optical_properties
        CALL GAS_OPTICAL_PROPERTIES(N_PROFILE, N_LAYER                  &
     &    , 1, I_GAS_POINTER, K_ESFT_MONO                               &
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
     &    , PHASE_FNC_CLR, FORWARD_SCATTER_CLR, PHASE_FNC_SOLAR_CLR     &
     &    , K_GREY_TOT, K_EXT_SCAT, PHASE_FNC                           &
     &    , FORWARD_SCATTER, PHASE_FNC_SOLAR                            &
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
!                     Calculated flxues
     &    , FLUX_DIRECT_PART, FLUX_TOTAL_PART                           &
!                       Calculated Radiances
     &    , RADIANCE_PART                                               &
!                       Calculated rates of photolysis
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
!       Increment the radiances within the band. Each increment
!       represents a single k-term within a band weighted with
!       its own weighting factor, hence for each increment the
!       weighting is the product of these two factors: similarly
!       for the blue flux.
        WEIGHT_INCR=WEIGHT_BAND*W_ESFT(IEX, I_BAND,  I_GAS)
        IF (L_BLUE_FLUX_SURF)                                           &
     &    WEIGHT_BLUE_INCR=WEIGHT_BLUE*W_ESFT(IEX, I_BAND,  I_GAS)
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
!
!
      RETURN
      END SUBROUTINE SOLVE_BAND_ONE_GAS
#endif
