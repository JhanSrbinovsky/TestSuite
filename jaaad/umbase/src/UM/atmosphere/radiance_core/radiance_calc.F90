#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the radiance field.
!
! Method:
!       Properties independent of the spectral bands are set.
!       a loop over bands is then entered. Grey optical properties
!       are set and an appropriate subroutine is called to treat
!       the gaseous overlaps. The final radiances are assigned.
!
! Current owner of code: James Manners
!
! Description of code:
!   Fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE RADIANCE_CALC(IERR                                     &
!                     Logical flags for processes
     &  , L_RAYLEIGH, L_AEROSOL, L_GAS, L_CONTINUUM                     &
     &  , L_CLOUD, L_DROP, L_ICE                                        &
!                     Angular integration
     &  , I_ANGULAR_INTEGRATION, L_RESCALE, N_ORDER_FORWARD             &
     &  , I_2STREAM                                                     &
     &  , N_ORDER_GAUSS                                                 &
     &  , I_TRUNCATION, LS_GLOBAL_TRUNC, MS_MIN, MS_MAX                 &
     &  , ACCURACY_ADAPTIVE, EULER_FACTOR                               &
     &  , L_HENYEY_GREENSTEIN_PF, LS_BRDF_TRUNC                         &
     &  , I_SPH_ALGORITHM, N_ORDER_PHASE_SOLAR                          &
     &  , N_DIRECTION, DIRECTION                                        &
     &  , N_VIEWING_LEVEL, VIEWING_LEVEL, I_SPH_MODE                    &
!                     Treatment of scattering
     &  , I_SCATTER_METHOD                                              &
!                       Options for treating clouds
     &  , L_GLOBAL_CLOUD_TOP, N_GLOBAL_CLOUD_TOP                        &
     &  , L_INHOM_CLOUD, INHOM_CLOUD                                    &
!                     Options for solver
     &  , I_SOLVER                                                      &
!                     Properties of diagnostics
     &  , MAP_CHANNEL                                                   &
!                     General spectral properties
     &  , N_BAND, I_FIRST_BAND, I_LAST_BAND, WEIGHT_BAND                &
!                     General atmospheric properties
     &  , N_PROFILE, N_LAYER                                            &
     &  , P, T, T_GROUND, T_LEVEL, D_MASS                               &
!                     Spectral region
     &  , ISOLIR                                                        &
!                     Solar fields
     &  , ZEN_0, SOLAR_IRRAD, SOLAR_FLUX_BAND                           &
     &  , RAYLEIGH_COEFFICIENT                                          &
!                     Infra-red fields
     &  , N_DEG_FIT, THERMAL_COEFFICIENT, T_REF_PLANCK                  &
     &  , L_IR_SOURCE_QUAD                                              &
!                     Gaseous absorption
     &  , I_GAS_OVERLAP, I_GAS                                          &
     &  , GAS_MIX_RATIO, N_BAND_ABSORB, INDEX_ABSORB                    &
     &  , I_BAND_ESFT, W_ESFT, K_ESFT, I_SCALE_ESFT                     &
     &  , I_SCALE_FNC, SCALE_VECTOR                                     &
     &  , P_REFERENCE, T_REFERENCE, L_MOD_K_FLUX                        &
!                     Doppler broadening
     &  , L_DOPPLER, DOPPLER_CORRECTION                                 &
!                     Surface fields
     &  , N_BRDF_BASIS_FNC, RHO_ALB, F_BRDF                             &
!                     Tiling options for heterogeneous surfaces
     &  , L_TILE, N_POINT_TILE, N_TILE, LIST_TILE, RHO_ALB_TILE         &
     &  , FRAC_TILE, T_TILE                                             &
!                     Continuum absorption
     &  , N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER                &
     &  , K_CONTINUUM, I_SCALE_FNC_CONT, SCALE_CONTINUUM                &
     &  , P_REF_CONTINUUM, T_REF_CONTINUUM                              &
!                     Properties of aerosols
     &  , N_AEROSOL, AEROSOL_MIX_RATIO                                  &
     &  , AEROSOL_ABSORPTION, AEROSOL_SCATTERING                        &
     &  , N_AEROSOL_PHF_TERM, AEROSOL_PHASE_FNC                         &
     &  , I_AEROSOL_PARAMETRIZATION, NHUMIDITY, HUMIDITIES              &
#if defined(OBSERVED)
     &  , N_OPT_LEVEL_AEROSOL_PRSC, N_PHASE_TERM_AEROSOL_PRSC           &
     &  , AEROSOL_PRESSURE_PRSC, AEROSOL_ABSORPTION_PRSC                &
     &  , AEROSOL_SCATTERING_PRSC, AEROSOL_PHASE_FNC_PRSC               &
#endif
!                     Properties of clouds
     &  , N_CONDENSED, TYPE_CONDENSED                                   &
     &  , I_CLOUD, I_CLOUD_REPRESENTATION, W_CLOUD                      &
     &  , N_CLOUD_TYPE, FRAC_CLOUD                                      &
     &  , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                       &
     &  , I_CONDENSED_PARAM, CONDENSED_N_PHF, CONDENSED_PARAM_LIST      &
     &  , DP_CORR_STRAT, DP_CORR_CONV                                   &
#if defined(OBSERVED)
     &  , N_OPT_LEVEL_DROP_PRSC, N_PHASE_TERM_DROP_PRSC                 &
     &  , DROP_PRESSURE_PRSC, DROP_ABSORPTION_PRSC                      &
     &  , DROP_SCATTERING_PRSC, DROP_PHASE_FNC_PRSC                     &
     &  , N_OPT_LEVEL_ICE_PRSC, N_PHASE_TERM_ICE_PRSC                   &
     &  , ICE_PRESSURE_PRSC, ICE_ABSORPTION_PRSC                        &
     &  , ICE_SCATTERING_PRSC, ICE_PHASE_FNC_PRSC                       &
#endif
!                     Calculated Fluxes or Radiances
     &  , FLUX_DIRECT, FLUX_DIFFUSE, FLUX_DOWN, FLUX_UP                 &
     &  , UV_FLUX_DIRECT, UV_FLUX_DOWN, UV_FLUX_UP                      &
     &  , L_FLUX_DIFFUSE, L_UVFLUX_DIRECT, L_UVFLUX_DOWN, L_UVFLUX_UP   &
     &  , RADIANCE, PHOTOLYSIS                                          &
!                     Options for clear-sky fluxes
     &  , L_CLEAR, I_SOLVER_CLEAR                                       &
!                     Clear-sky fluxes calculated
     &  , FLUX_DIRECT_CLEAR, FLUX_DOWN_CLEAR, FLUX_UP_CLEAR             &
!                     Special Surface Fluxes
     &  , L_BLUE_FLUX_SURF, WEIGHT_BLUE,WEIGHT_UV                       &
     &  , FLUX_DIRECT_BLUE_SURF                                         &
     &  , FLUX_DOWN_BLUE_SURF, FLUX_UP_BLUE_SURF                        &
!                     Tiled Surface Fluxes
     &  , FLUX_UP_TILE, FLUX_UP_BLUE_TILE                               &
!                     Dimensions of arrays
     &  , ND_PROFILE, ND_LAYER, ND_COLUMN, ND_LAYER_CLR, ID_CT          &
     &  , ND_2SG_PROFILE, ND_FLUX_PROFILE, ND_RADIANCE_PROFILE          &
     &  , ND_J_PROFILE                                                  &
     &  , ND_CHANNEL, ND_BAND                                           &
     &  , ND_SPECIES, ND_ESFT_TERM, ND_SCALE_VARIABLE                   &
     &  , ND_CONTINUUM                                                  &
     &  , ND_AEROSOL_SPECIES, ND_HUMIDITIES                             &
     &  , ND_CLOUD_PARAMETER                                            &
     &  , ND_THERMAL_COEFF, ND_SOURCE_COEFF                             &
     &  , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC                              &
#if defined(OBSERVED)
     &  , ND_PROFILE_AEROSOL_PRSC, ND_PROFILE_CLOUD_PRSC                &
     &  , ND_OPT_LEVEL_AEROSOL_PRSC, ND_OPT_LEVEL_CLOUD_PRSC            &
#endif
     &  , ND_PHASE_TERM, ND_MAX_ORDER, ND_SPH_COEFF                     &
     &  , ND_DIRECTION, ND_VIEWING_LEVEL                                &
     &  , ND_REGION, ND_CLOUD_TYPE, ND_CLOUD_COMPONENT                  &
     &  , ND_OVERLAP_COEFF                                              &
     &  , ND_POINT_TILE, ND_TILE, I_CALL                                &
     &  )
!
!
!
      IMPLICIT NONE
!
! Include Header Files
#include "c_kinds.h"
!
!     Sizes of dummy arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for totally clear layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_2SG_PROFILE                                                &
!           Size allocated for profiles of fluxes
     &  , ND_FLUX_PROFILE                                               &
!           Size allocated for profiles of output fluxes
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for profiles of radiances
     &  , ND_J_PROFILE                                                  &
!           Size allocated for profiles of mean radiances
     &  , ND_CHANNEL                                                    &
!           Size allocated for channels of output
     &  , ND_BAND                                                       &
!           Size allocated for bands in spectral computation
     &  , ND_SPECIES                                                    &
!           Size allocated for gaseous species
     &  , ND_CONTINUUM                                                  &
!           Size allocated for types of continua
     &  , ND_AEROSOL_SPECIES                                            &
!           Size allocated for aerosol species
     &  , ND_HUMIDITIES                                                 &
!           Size allocated for humidities
     &  , ND_ESFT_TERM                                                  &
!           Size allocated for ESFT terms
     &  , ND_SCALE_VARIABLE                                             &
!           Size allocated for variables in scaling functions
     &  , ND_CLOUD_PARAMETER                                            &
!           Size allocated for cloud parameters
     &  , ND_THERMAL_COEFF                                              &
!           Size allocated for thermal coefficients
     &  , ND_SOURCE_COEFF                                               &
!           Size allocated for two-stream source coefficients
     &  , ND_BRDF_BASIS_FNC                                             &
!           Size allowed for BRDF basis functions
     &  , ND_BRDF_TRUNC                                                 &
!           Size allowed for truncation of BRDF basis functions
     &  , ND_COLUMN                                                     &
!           Size allocated for columns at each grid-point
     &  , ND_PHASE_TERM                                                 &
!           Size allocated for terms in the phase function
!           supplied to the routine
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for polar orders
     &  , ND_DIRECTION                                                  &
!           Size allocated for viewing directions at each point
     &  , ND_VIEWING_LEVEL                                              &
!           Size allocated for levels where the radiance
!           may be calculated
     &  , ND_REGION                                                     &
!           Size allocated for cloudy regions
     &  , ND_CLOUD_TYPE                                                 &
!           Size allocated for types of clouds
     &  , ND_CLOUD_COMPONENT                                            &
!           Size allocated for components in clouds
     &  , ND_OVERLAP_COEFF                                              &
!           Size allocated for overlap coefficients
     &  , ND_SPH_COEFF                                                  &
!           Size allocated for arrays of spherical coefficients
!           used in determining radiances
#if defined(OBSERVED)
     &  , ND_PROFILE_AEROSOL_PRSC                                       &
!           Size allocated for profiles of prescribed
!           aerosol optical properties
     &  , ND_PROFILE_CLOUD_PRSC                                         &
!           Size allocated for profiles of prescribed
!           cloudy optical properties
     &  , ND_OPT_LEVEL_AEROSOL_PRSC                                     &
!           Size allocated for levels of prescribed
!           aerosol optical properties
     &  , ND_OPT_LEVEL_CLOUD_PRSC                                       &
!           Size allocated for levels of prescribed
!           cloudy optical properties
#endif
     &  , ND_POINT_TILE                                                 &
!           Size allocated for points with surface tiling
     &  , ND_TILE
!           Size allocated for the number of tiles
!
      INTEGER, INTENT(IN) ::                                            &
     &    I_CALL
!           Number of radiation Call
!
!     Include header files.
#include "def_std_io_icf3z.h"
#include "gas_overlap_pcf3z.h"
#include "cloud_scheme_pcf3z.h"
#include "cloud_region_pcf3z.h"
#include "angular_integration_pcf3z.h"
#include "sph_algorithm_pcf3z.h"
#include "spectral_region_pcf3z.h"
#include "aerosol_parametrization_pcf3z.h"
#include "esft_scale_pcf3z.h"
#include "error_pcf3z.h"
!
!
!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
!
!     General logical switches:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLEAR                                                       &
!           Calculate clear-sky fluxes
     &  , L_IR_SOURCE_QUAD                                              &
!           Use a quadratic source function
     &  , L_RESCALE                                                     &
!           Flag for delta-rescaling
     &  , L_HENYEY_GREENSTEIN_PF
!           Use Henyey-Greenstein phase functions
!
!     Parameters controlling algorithms:
!     Representation of clouds:
      INTEGER, INTENT(IN) ::                                            &
     &    I_CLOUD
!           Cloud scheme used
!     Numerical algorithms:
      INTEGER, INTENT(IN) ::                                            &
     &    MAP_CHANNEL(ND_BAND)
!           Mapping of actual bands to the output channels
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR                                                        &
!           Visible or IR
     &  , I_SOLVER                                                      &
!           Solver used
     &  , I_SOLVER_CLEAR                                                &
!           Clear solver used
     &  , I_2STREAM                                                     &
!           Two-stream scheme
     &  , I_ANGULAR_INTEGRATION                                         &
!           Angular integration scheme
     &  , N_ORDER_GAUSS                                                 &
!           Order of Gaussian integration
     &  , I_TRUNCATION                                                  &
!           Type of spherical truncation
     &  , LS_GLOBAL_TRUNC                                               &
!           Truncating order of spherical harmonics
     &  , MS_MIN                                                        &
!           Lowest azimuthal order calculated
     &  , MS_MAX                                                        &
!           Highest azimuthal order calculated
     &  , N_ORDER_FORWARD                                               &
!           Order of the term used to `define' the forward scattering
!           fraction.
     &  , I_SPH_MODE                                                    &
!           Mode in which the spherical harmonic solver is being used
     &  , I_SPH_ALGORITHM
!           Algorithm used for spherical harmonic calculation
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ACCURACY_ADAPTIVE                                             &
!           Accuracy for adaptive truncation
     &  , EULER_FACTOR
!           Factor applied to the last term of an alternating series
      INTEGER, INTENT(INOUT) ::                                         &
     &    LS_BRDF_TRUNC
!           Order of truncation applied to BRDFs
!           (This will be reset to 0 if a Lambertian surface
!           is assumed)
!
!     Specification of the viewing geometry
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION                                                   &
!           Number of directions at which to calculate radiances
     &  , N_VIEWING_LEVEL
!           Number of levels where the radiance is required
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DIRECTION(ND_RADIANCE_PROFILE, ND_DIRECTION, 2)               &
!           Directions in which to calculate radiances
     &  , VIEWING_LEVEL(ND_VIEWING_LEVEL)
!           List of levels where the radiance is required
!
!     Range of spectral bands:
      INTEGER, INTENT(IN) ::                                            &
     &    I_FIRST_BAND                                                  &
!           First band
     &  , I_LAST_BAND
!           Last band
!
!     General properties of spectrum:
      INTEGER, INTENT(IN) ::                                            &
     &    N_BAND                                                        &
!           Number of spectral bands
     &  , N_AEROSOL
!           Number of aerosol species
!
!     Solar fields:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    SOLAR_IRRAD(ND_PROFILE)                                       &
!           Incident solar radiation
     &  , SOLAR_FLUX_BAND(ND_BAND)                                      &
!           Normalized flux in each spectral band
     &  , ZEN_0(ND_PROFILE)
!           Secant (two-stream) or cosine (spherical harmonics)
!           of solar zenith angle
!
!     Atmospheric profiles:
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER
!           Number of layers
      REAL  (Real64), INTENT(IN) ::                                     &
     &    P(ND_PROFILE, ND_LAYER)                                       &
!           Pressure
     &  , T(ND_PROFILE, ND_LAYER)                                       &
!           Temperature
     &  , T_GROUND(ND_PROFILE)                                          &
!           Temperature of ground
     &  , T_LEVEL(ND_PROFILE, 0: ND_LAYER)                              &
!           Temperature on levels
     &  , D_MASS(ND_PROFILE, ND_LAYER)                                  &
!           Mass thickness of each layer
     &  , GAS_MIX_RATIO(ND_PROFILE, ND_LAYER, ND_SPECIES)
!           Gaseous mass mixing ratios
!
!     Surface properties:
      INTEGER, INTENT(IN) ::                                            &
     &    N_BRDF_BASIS_FNC
!           Number of BRDF basis functions
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    RHO_ALB(ND_PROFILE, ND_BRDF_BASIS_FNC, ND_BAND)               &
!           Weights of the basis functions
     &  , F_BRDF(ND_BRDF_BASIS_FNC, 0: ND_BRDF_TRUNC/2                  &
     &      , 0: ND_BRDF_TRUNC/2, 0: ND_BRDF_TRUNC)
!           Array of BRDF basis terms
!
!     Arrays related to tiling of the surface
      LOGICAL, INTENT(IN) ::                                            &
     &    L_TILE
!           Local to allow tiling options
      INTEGER, INTENT(IN) ::                                            &
     &    N_POINT_TILE                                                  &
!           Number of points to tile
     &  , N_TILE                                                        &
!           Number of tiles used
     &  , LIST_TILE(ND_POINT_TILE)
!           List of points with surface tiling
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    RHO_ALB_TILE(ND_POINT_TILE, ND_BRDF_BASIS_FNC                 &
     &      , ND_TILE, ND_BAND)                                         &
!           Weights for the basis functions of the BRDFs
!           at the tiled points
     &  , FRAC_TILE(ND_POINT_TILE, ND_TILE)                             &
!           Fraction of tiled grid-points occupied by each tile
     &  , T_TILE(ND_POINT_TILE, ND_TILE)
!           Local surface temperatures on individual tiles
!
!     Rayleigh scattering:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_RAYLEIGH
!           Include rayleigh scattering in the calculation.
      REAL  (Real64), INTENT(IN) ::                                     &
     &     RAYLEIGH_COEFFICIENT(ND_BAND)
!           Rayleigh coefficients
!
!     fields for gaseous absorption:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_GAS
!           Include gas absorption in the calculation
!     gaseous overlaps:
      INTEGER, INTENT(IN) ::                                            &
     &    I_GAS_OVERLAP(ND_BAND)                                        &
!           Gas overlap assumption
     &  , I_GAS
!           Gas to be considered (one gas only)
!     ESFTs:
      INTEGER, INTENT(IN) ::                                            &
     &    N_BAND_ABSORB(ND_BAND)                                        &
!           Number of absorbers in band
     &  , INDEX_ABSORB(ND_SPECIES, ND_BAND)                             &
!           List of absorbers in bands
     &  , I_BAND_ESFT(ND_BAND, ND_SPECIES)                              &
!           Number of terms in band
     &  , I_SCALE_ESFT(ND_BAND, ND_SPECIES)                             &
!           Type of esft scaling
     &  , I_SCALE_FNC(ND_BAND, ND_SPECIES)
!           Type of scaling function
      REAL  (Real64), INTENT(IN) ::                                     &
     &    W_ESFT(ND_ESFT_TERM, ND_BAND, ND_SPECIES)                     &
!           Weights for ESFT
     &  , K_ESFT(ND_ESFT_TERM, ND_BAND, ND_SPECIES)                     &
!           Exponential ESFT terms
     &  , SCALE_VECTOR(ND_SCALE_VARIABLE, ND_ESFT_TERM, ND_BAND         &
     &        , ND_SPECIES)                                             &
!           Absorber scaling parameters
     &  , P_REFERENCE(ND_SPECIES, ND_BAND)                              &
!           Reference scaling pressure
     &  , T_REFERENCE(ND_SPECIES, ND_BAND)
!           Reference scaling temperature
!
!     Use modulus of fluxes to remove negative effective extinctions
      LOGICAL, INTENT(IN) :: L_MOD_K_FLUX

!     Spectral data for the continuum:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CONTINUUM
!           Include continuum absorption in the calculation
      INTEGER, INTENT(IN) ::                                            &
     &    N_BAND_CONTINUUM(ND_BAND)                                     &
!           Number of continua in bands
     &  , INDEX_CONTINUUM(ND_BAND, ND_CONTINUUM)                        &
!           Indices of continua
     &  , INDEX_WATER                                                   &
!           Index of water
     &  , I_SCALE_FNC_CONT(ND_BAND, ND_CONTINUUM)
!           Type of scaling function for continuum
      REAL  (Real64), INTENT(IN) ::                                     &
     &    K_CONTINUUM(ND_BAND, ND_CONTINUUM)                            &
!           Continuum extinction coefficients
     &  , SCALE_CONTINUUM(ND_SCALE_VARIABLE, ND_BAND, ND_CONTINUUM)     &
!           Continuum scaling parameters
     &  , P_REF_CONTINUUM(ND_CONTINUUM, ND_BAND)                        &
!           Continuum reference pressure
     &  , T_REF_CONTINUUM(ND_CONTINUUM, ND_BAND)
!           Continuum reference temperature
!
!
!     General cloud fields:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLOUD
!           Clouds are required in the calculation
      REAL  (Real64), INTENT(IN) ::                                     &
     &    W_CLOUD(ND_PROFILE, ID_CT: ND_LAYER)                          &
!           Amount of cloud
     &  , FRAC_CLOUD(ND_PROFILE, ID_CT: ND_LAYER, ND_CLOUD_TYPE)        &
!           Fractions of different types of cloud
     &  , DP_CORR_STRAT                                                 &
!           Decorrelation pressure scale for large scale cloud
     &  , DP_CORR_CONV
!           Decorrelation pressure scale for convective cloud
!
!     Fields for microphysical quantities:
!
      LOGICAL, INTENT(IN) ::                                            &
     &    L_DROP                                                        &
!           Include droplets in the calculation
     &  , L_ICE                                                         &
!           Include ice in the calculation
     &  , L_GLOBAL_CLOUD_TOP                                            &
!           Use a global value for the top of clouds
!           (This is for use in a GCM where the code must be
!           bit-reproducible across different configurations of PEs).
     &  , L_INHOM_CLOUD
!           Use scaling factors for inhomogeneous cloud
      INTEGER, INTENT(IN) ::                                            &
     &    N_CONDENSED                                                   &
!           Number of condensed components in clouds
     &  , TYPE_CONDENSED(ND_CLOUD_COMPONENT)                            &
!           Types of condensed components
     &  , I_CONDENSED_PARAM(ND_CLOUD_COMPONENT)                         &
!           Parametrization schemes for components
     &  , CONDENSED_N_PHF(ND_CLOUD_COMPONENT)                           &
!           Number of terms in the phase function
     &  , I_CLOUD_REPRESENTATION                                        &
!           Representation of mixing rule chosen
     &  , N_CLOUD_TYPE                                                  &
!           Number of types of cloud
     &  , N_GLOBAL_CLOUD_TOP
!           Global cloud top
!
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    CONDENSED_MIX_RATIO(ND_PROFILE, ID_CT: ND_LAYER               &
     &      , ND_CLOUD_COMPONENT)
!           Mixing ratios of condensed components

      REAL  (Real64), INTENT(IN) ::                                     &
     &    CONDENSED_DIM_CHAR(ND_PROFILE, ID_CT: ND_LAYER                &
     &      , ND_CLOUD_COMPONENT)                                       &
!           Characteristic dimensions of condensed components
     &  , CONDENSED_PARAM_LIST(ND_CLOUD_PARAMETER                       &
     &      , ND_CLOUD_COMPONENT, ND_BAND)                              &
!           Coefficients in parametrizations of condensed phases
     &  , INHOM_CLOUD(ND_CLOUD_COMPONENT)
!           Scaling factors for inhomogeneous cloud
!
#if defined(OBSERVED)
!
!     Fields for prescribed optical properties of droplets
      INTEGER, INTENT(IN) ::                                            &
     &    N_OPT_LEVEL_DROP_PRSC                                         &
!           Number of levels of prescribed
!           optical properties of droplets
     &  , N_PHASE_TERM_DROP_PRSC
!           Number of terms in the phase function for prescribed
!           water droplets
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DROP_PRESSURE_PRSC(ND_PROFILE_CLOUD_PRSC                      &
     &      , ND_OPT_LEVEL_CLOUD_PRSC)                                  &
!           Pressures at which optical properties of
!           droplets are prescribed
     &  , DROP_ABSORPTION_PRSC(ND_PROFILE_CLOUD_PRSC                    &
     &      , ND_OPT_LEVEL_CLOUD_PRSC, ND_BAND)                         &
!           Prescribed absorption by droplets
     &  , DROP_SCATTERING_PRSC(ND_PROFILE_CLOUD_PRSC                    &
     &      , ND_OPT_LEVEL_CLOUD_PRSC, ND_BAND)                         &
!           Prescribed scattering by droplets
     &  , DROP_PHASE_FNC_PRSC(ND_PROFILE_CLOUD_PRSC                     &
     &      , ND_OPT_LEVEL_CLOUD_PRSC, ND_PHASE_TERM, ND_BAND)
!           Prescribed phase function of droplets
!
!     Fields for prescribed optical properties of ice crystals
      INTEGER, INTENT(IN) ::                                            &
     &    N_OPT_LEVEL_ICE_PRSC                                          &
!           Number of levels of prescribed
!           optical properties of ice crystals
     &  , N_PHASE_TERM_ICE_PRSC
!           Number of terms in the phase function for prescribed
!           ice crystals
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ICE_PRESSURE_PRSC(ND_PROFILE_CLOUD_PRSC                       &
     &      , ND_OPT_LEVEL_CLOUD_PRSC)                                  &
!           Pressures at which optical properties of
!           ice crystals are prescribed
     &  , ICE_ABSORPTION_PRSC(ND_PROFILE_CLOUD_PRSC                     &
     &      , ND_OPT_LEVEL_CLOUD_PRSC, ND_BAND)                         &
!           Prescribed absorption by ice crystals
     &  , ICE_SCATTERING_PRSC(ND_PROFILE_CLOUD_PRSC                     &
     &      , ND_OPT_LEVEL_CLOUD_PRSC, ND_BAND)                         &
!           Prescribed scattering by ice crystals
     &  , ICE_PHASE_FNC_PRSC(ND_PROFILE_CLOUD_PRSC                      &
     &      , ND_OPT_LEVEL_CLOUD_PRSC, ND_PHASE_TERM, ND_BAND)
!           Prescribed phase functions of ice crystals
#endif
!
!
!     Fields for aerosols:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_AEROSOL
!           Include aerosols in the calculation
      INTEGER, INTENT(IN) ::                                            &
     &    I_AEROSOL_PARAMETRIZATION(ND_AEROSOL_SPECIES)                 &
!           Parametrization flags for aerosol
     &  , N_AEROSOL_PHF_TERM(ND_AEROSOL_SPECIES)
!           Number of terms in the phase function of aerosols
      INTEGER, INTENT(IN) ::                                            &
     &    NHUMIDITY(ND_AEROSOL_SPECIES)
!           Number of humidities
      REAL  (Real64), INTENT(IN) ::                                     &
     &    AEROSOL_MIX_RATIO(ND_PROFILE, ND_LAYER                        &
     &      , ND_AEROSOL_SPECIES)
!           Number density of aerosols
      REAL  (Real64), INTENT(IN) ::                                     &
     &    AEROSOL_ABSORPTION(ND_HUMIDITIES, ND_AEROSOL_SPECIES          &
     &      , ND_BAND)                                                  &
!           Absorption by aerosols
     &  , AEROSOL_SCATTERING(ND_HUMIDITIES, ND_AEROSOL_SPECIES          &
     &        , ND_BAND)                                                &
!           Scattering by aerosols
     &  , AEROSOL_PHASE_FNC(ND_HUMIDITIES, ND_PHASE_TERM                &
     &        , ND_AEROSOL_SPECIES, ND_BAND)                            &
!           Phase function of aerosols
     &  , HUMIDITIES(ND_HUMIDITIES, ND_AEROSOL_SPECIES)
!           Humidities for species
!
#if defined(OBSERVED)
      INTEGER, INTENT(IN) ::                                            &
     &    N_OPT_LEVEL_AEROSOL_PRSC(ND_AEROSOL_SPECIES)                  &
!           Number of levels of prescribed optical properties
!           of aerosols
     &  , N_PHASE_TERM_AEROSOL_PRSC(ND_AEROSOL_SPECIES)
!           Number of terms in the phase function for prescribed
!           aerosols
      REAL  (Real64), INTENT(IN) ::                                     &
     &    AEROSOL_PRESSURE_PRSC(ND_PROFILE_AEROSOL_PRSC                 &
     &      , ND_OPT_LEVEL_AEROSOL_PRSC, ND_AEROSOL_SPECIES)            &
!           Pressures at which optical properties of aerosols
!           are prescribed
     &  , AEROSOL_ABSORPTION_PRSC(ND_PROFILE_AEROSOL_PRSC               &
     &      , ND_OPT_LEVEL_AEROSOL_PRSC, ND_AEROSOL_SPECIES, ND_BAND)   &
!           Prescribed absorption by aerosols
     &  , AEROSOL_SCATTERING_PRSC(ND_PROFILE_AEROSOL_PRSC               &
     &      , ND_OPT_LEVEL_AEROSOL_PRSC, ND_AEROSOL_SPECIES, ND_BAND)   &
!           Prescribed scattering by aerosols
     &  , AEROSOL_PHASE_FNC_PRSC(ND_PROFILE_AEROSOL_PRSC                &
     &      , ND_OPT_LEVEL_AEROSOL_PRSC                                 &
     &      , ND_PHASE_TERM, ND_AEROSOL_SPECIES, ND_BAND)
!           Prescribed phase functions of aerosols
#endif
!
!     Fitting of the Planckian function:
      INTEGER, INTENT(IN) ::                                            &
     &    N_DEG_FIT
!           Degree of thermal fitting fnc.
      REAL  (Real64), INTENT(IN) ::                                     &
     &    THERMAL_COEFFICIENT(0: ND_THERMAL_COEFF-1, ND_BAND)           &
!           Coefficients of source terms
     &  , T_REF_PLANCK
!           Planckian reference temperature
!
!     Doppler broadening
      LOGICAL, INTENT(IN) ::                                            &
     &    L_DOPPLER(ND_SPECIES)
!           Flags to activate doppler corrections
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DOPPLER_CORRECTION(ND_SPECIES)
!           Doppler broadening term
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    WEIGHT_BAND(ND_BAND)
!           Weighting function for bands
!
!     Control of scattering:
      INTEGER, INTENT(IN) ::                                            &
     &    I_SCATTER_METHOD(ND_BAND)
!           Method of treating scattering in each band
!
!     Fluxes or radiances calculated:
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_FLUX_PROFILE, 0: ND_LAYER, ND_CHANNEL)         &
!           Direct flux
     &  , FLUX_DIFFUSE(ND_FLUX_PROFILE, 0: ND_LAYER, ND_CHANNEL)        &
!           Diffuse downward flux     
     &  , FLUX_DOWN(ND_FLUX_PROFILE, 0: ND_LAYER, ND_CHANNEL)           &
!           Downward flux
     &  , FLUX_UP(ND_FLUX_PROFILE, 0: ND_LAYER, ND_CHANNEL)             &
!           Upward flux
     &  , FLUX_DIRECT_CLEAR(ND_2SG_PROFILE, 0: ND_LAYER, ND_CHANNEL)    &
!           Clear direct flux
     &  , FLUX_DOWN_CLEAR(ND_2SG_PROFILE, 0: ND_LAYER, ND_CHANNEL)      &
!           Clear downward flux
     &  , FLUX_UP_CLEAR(ND_2SG_PROFILE, 0: ND_LAYER, ND_CHANNEL)        &
!           Clear upward flux
     &  , UV_FLUX_DIRECT(ND_FLUX_PROFILE, 0: ND_LAYER, ND_CHANNEL)      &
!           Direct UV-flux
     &  , UV_FLUX_UP(ND_FLUX_PROFILE, 0: ND_LAYER, ND_CHANNEL)          &
!           Upwards UV-flux
     &  , UV_FLUX_DOWN(ND_FLUX_PROFILE, 0: ND_LAYER, ND_CHANNEL)        &
!           Downward UV-Flux
     &  , WEIGHT_UV(ND_BAND)                                            &
!           Weights fo each band for the UV-interval
     &  , RADIANCE(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL                &
     &      , ND_DIRECTION, ND_CHANNEL)                                 &
!           Calculated radiances
     &  , PHOTOLYSIS(ND_J_PROFILE, ND_VIEWING_LEVEL, ND_CHANNEL)
!           Rate of photolysis
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_UP_TILE(ND_POINT_TILE, ND_TILE, ND_CHANNEL)              &
!           Upward fluxes at tiled surface points
     &  , FLUX_UP_BLUE_TILE(ND_POINT_TILE, ND_TILE, ND_CHANNEL)
!           Upward blue fluxes at tiled surface points

!
!     Special Diagnostics:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_BLUE_FLUX_SURF                                              &
!           Flag to calculate the blue surface fluxes
     &   , L_FLUX_DIFFUSE                                               &
!           Flag for the diffuse downward flux     
     &   , L_UVFLUX_DIRECT                                              &
!             FLAG FOR DIRECT ULTRAVIOLET FLUXES
     &   , L_UVFLUX_UP                                                  &
!             FLAF FOR UPWARD ULTRAVIOLET FLUXES
     &   , L_UVFLUX_DOWN
!             FLAG FOR DOWNWARD ULTRAVIOLET FLUXES
      REAL  (Real64), INTENT(IN) ::                                     &
     &    WEIGHT_BLUE(ND_BAND)
!           Weights for each band for blue fluxes
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT_BLUE_SURF(ND_FLUX_PROFILE)                        &
!           Direct blue flux at the surface
     &  , FLUX_DOWN_BLUE_SURF(ND_FLUX_PROFILE)                          &
!           Total downward blue flux at the surface
     &  , FLUX_UP_BLUE_SURF(ND_FLUX_PROFILE)
!           Upward blue flux at the surface
!
!
!
!     Local arguments.
!     General pointers:
      INTEGER                                                           &
     &    I_TOP                                                         &
!           Top level of profiles
     &  , I_BAND                                                        &
!           Spectral band
     &  , N_GAS                                                         &
!           Number of active gases
     &  , I_GAS_BAND                                                    &
!           Single variable for gas in band
     &  , N_CONTINUUM                                                   &
!           Number of continua in band
     &  , I_CONTINUUM                                                   &
!           Continuum number
     &  , I_CONTINUUM_POINTER(ND_CONTINUUM)
!           Pointers to continua
#if ! defined(UM)
      INTEGER                                                           &
     &    I_POINTER_WATER
!           Pointer to water vapour
#endif
!
!     Additional variables for angular integration:
      LOGICAL                                                           &
     &    L_SOLAR_PHF                                                   &
!           Logical to specify a separate treatment of the singly
!           scattered solar beam
     &  , L_RESCALE_SOLAR_PHF
!           Logical to apply rescaling to the singly scattered
!           solar phase function
      INTEGER                                                           &
     &    N_ORDER_PHASE                                                 &
!           Order of Legendre polynomials of the phase function
     &  , N_ORDER_PHASE_SOLAR
!           Order of Legendre polynomials of the phase function
!           retained for the singly scattered solar beam
!
!     Pointers to the contents of layers:
      INTEGER                                                           &
     &    N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_REGION                                                      &
!           Number of cloudy regions
     &  , N_CLOUD_PROFILE(ID_CT: ND_LAYER)                              &
!           Number of cloudy profiles
     &  , I_CLOUD_PROFILE(ND_PROFILE, ID_CT: ND_LAYER)
!           Profiles containing clouds
!
!     Pointers to types of clouds:
      LOGICAL                                                           &
     &    L_CLOUD_CMP(ND_CLOUD_COMPONENT)
!           Logical switches to `include' components
      INTEGER                                                           &
     &    I_PHASE_CMP(ND_CLOUD_COMPONENT)                               &
!           Phases of components
     &  , I_CLOUD_TYPE(ND_CLOUD_COMPONENT)                              &
!           Pypes of cloud to which each component contributes
     &  , TYPE_REGION(ND_REGION)                                        &
!           The types of the regions
     &  , K_CLR                                                         &
!           Index of clear-sky region
     &  , I_REGION_CLOUD(ND_CLOUD_TYPE)
!           Regions in which particular type of cloud fall
!
!     Fractional coverage of different regions:
      REAL  (Real64) ::                                                 &
     &    FRAC_REGION(ND_PROFILE, ID_CT: ND_LAYER, ND_REGION)
!           Fraction of total cloud occupied by specific regions
!
      REAL ::                                                           &
     &     TOT_CLOUD_COVER(ND_PROFILE)
!             Total cloud cover

!     Pointer to table of humidities:
      INTEGER                                                           &
     &    I_HUMIDITY_POINTER(ND_PROFILE, ND_LAYER)
!           Pointer to look-up table for aerosols
!
!     Controlling variables:
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , K                                                             &
!           Loop variable
     &  , L
!           Loop variable
!
!     Logical switches:
      LOGICAL                                                           &
     &    L_GAS_BAND                                                    &
!           Flag to `include' gaseous absorption in a particular band
     &  , L_MOIST_AEROSOL                                               &
!           Flag for moist aerosol
     &  , L_AEROSOL_DENSITY
!           Flag for calculation of atmospheric density for aerosols
!
      REAL  (Real64) ::                                                 &
     &    SOLAR_IRRAD_BAND(ND_PROFILE)
!           Solar irradiance in the band
      REAL  (Real64) ::                                                 &
     &    GAS_FRAC_RESCALED(ND_PROFILE, ND_LAYER, ND_SPECIES)           &
!           Rescaled gas mixing ratios
     &  , AMOUNT_CONTINUUM(ND_PROFILE, ND_LAYER, ND_CONTINUUM)          &
!           Amounts of continua
     &  , K_CONTINUUM_MONO(ND_CONTINUUM)
!           Monochromatic continuum components
!
!     Thermal arrays:
      REAL  (Real64) ::                                                 &
     &    PLANCK_FLUX_BAND(ND_PROFILE, 0: ND_LAYER)                     &
!           Planckian flux in band at edges of layers
     &  , DIFF_PLANCK_BAND(ND_PROFILE, ND_LAYER)                        &
!           Difference in the Planckian flux across layers
     &  , DIFF_PLANCK_BAND_2(ND_PROFILE, ND_LAYER)                      &
!           Twice the 2nd difference of in the Planckian flux across
!           layers
     &  , PLANCK_FLUX_GROUND(ND_PROFILE)
!           Planckian flux at the surface temperature
!
!     Surface BRDF terms
      LOGICAL                                                           &
     &    L_DIFF_ALB
!           Flag to calculate diffuse albedos
      REAL  (Real64) ::                                                 &
     &    BRDF_SOL(ND_PROFILE, ND_BRDF_BASIS_FNC, ND_DIRECTION)         &
!           The BRDF evaluated for scattering from the solar
!           beam into the viewing direction
     &  , BRDF_HEMI(ND_PROFILE, ND_BRDF_BASIS_FNC, ND_DIRECTION)        &
!           The BRDF evaluated for scattering from isotropic
!           radiation into the viewing direction
     &  , DIFFUSE_ALB_BASIS(ND_BRDF_BASIS_FNC)
!           The diffuse albedo of isotropic radiation for each
!           basis function
!
!     Atmospheric densities:
      REAL  (Real64) ::                                                 &
     &    DENSITY(ND_PROFILE, ND_LAYER)                                 &
!           Overall density
     &  , MOLAR_DENSITY_WATER(ND_PROFILE, ND_LAYER)                     &
!           Molar density of water
     &  , MOLAR_DENSITY_FRN(ND_PROFILE, ND_LAYER)
!           Molar density of foreign species
!
!     Fields for moist aerosols:
      REAL  (Real64) ::                                                 &
     &    DELTA_HUMIDITY                                                &
!           Increment in look-up table for hum.
     &  , MEAN_REL_HUMIDITY(ND_PROFILE, ND_LAYER)
!           Mean relative humidity of layers
!
!     Fundamental optical properties of layers:
!                       Clear-sky optical properties
      REAL  (Real64) ::                                                 &
     &    K_GREY_TOT_CLR(ND_PROFILE, ND_LAYER_CLR)                      &
!           Clear-sky total grey extinction above the layer
!           where clouds are allowed
     &  , K_EXT_SCAT_CLR(ND_PROFILE, ND_LAYER_CLR)                      &
!           Clear-sky scattering extinction
     &  , PHASE_FNC_CLR(ND_PROFILE, ND_LAYER_CLR, ND_MAX_ORDER)         &
!           Clear-sky phase function
     &  , PHASE_FNC_SOLAR_CLR(ND_RADIANCE_PROFILE, ND_LAYER_CLR         &
     &      , ND_DIRECTION)                                             &
!           Solar phase function for singly scattered radiation
     &  , FORWARD_SCATTER_CLR(ND_PROFILE, ND_LAYER_CLR)                 &
!           Clear-sky forward scattering
     &  , FORWARD_SOLAR_CLR(ND_PROFILE, ND_LAYER_CLR)
!           Clear-sky forward scattering for the solar beam
!                       All-sky optical properties
      REAL  (Real64) ::                                                 &
     &    K_GREY_TOT(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)     &
!           Total grey extinction
     &  , K_EXT_SCAT(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)     &
!           Scattering extinction
     &  , PHASE_FNC(ND_PROFILE, ID_CT: ND_LAYER                         &
     &        , ND_MAX_ORDER, 0: ND_CLOUD_TYPE)                         &
!           Phase function
     &  , PHASE_FNC_SOLAR(ND_RADIANCE_PROFILE, ID_CT: ND_LAYER          &
     &      , ND_DIRECTION, 0: ND_CLOUD_TYPE)                           &
!           Solar phase function for singly scattered radiation
     &  , FORWARD_SCATTER(ND_PROFILE, ID_CT: ND_LAYER                   &
     &      , 0: ND_CLOUD_TYPE)                                         &
!           Forward scattering
     &  , FORWARD_SOLAR(ND_PROFILE, ID_CT: ND_LAYER                     &
     &      , 0: ND_CLOUD_TYPE)
!           Forward scattering for the solar beam
!
!     Local variables for spherical harmonic integration
      INTEGER                                                           &
     &    LS_MAX_ORDER                                                  &
!           Maximum order of terms required
     &  , LS_LOCAL_TRUNC(0: ND_MAX_ORDER)                               &
!           Actual truncation for each particular value of m
     &  , MS_TRUNC(0: ND_MAX_ORDER)                                     &
!           Maximum azimuthal quantum number for each order
     &  , IA_SPH_MM(0: ND_MAX_ORDER)
!           Address of spherical coefficient of (m, m) for each m
      REAL  (Real64) ::                                                 &
     &    CG_COEFF(ND_SPH_COEFF)                                        &
!           Clebsch-gordon coefficients
     &  , UPLM_ZERO(ND_SPH_COEFF)                                       &
!           Upsilon terms
     &  , UPLM_SOL(ND_RADIANCE_PROFILE, ND_SPH_COEFF)                   &
!           Upsilon terms for solar radiation
     &  , COS_SOL_VIEW(ND_RADIANCE_PROFILE, ND_DIRECTION)
!           Cosines of the angles between the solar direction
!           and the viewing direction
!
!     Specification of the grid for radiances:
      INTEGER                                                           &
     &    I_RAD_LAYER(ND_VIEWING_LEVEL)
!           Layers in which to intercept radiances
      REAL  (Real64) ::                                                 &
     &    FRAC_RAD_LAYER(ND_VIEWING_LEVEL)
!           Fractions below the tops of the layers
!
      REAL  (Real64) ::                                                 &
     &    I_DIRECT(ND_RADIANCE_PROFILE, 0: ND_LAYER)
!           Direct solar irradiance on levels (not split by
!           diagnostic bands or returned, but retained for
!           future use)
      REAL  (Real64) ::                                                 &
     &    PLANCK_RADIANCE_BAND(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL)
!           Planckian radiance in the current band
      LOGICAL                                                           &
     &    L_INITIAL
!           Flag to run the routine incrementing diagnostics in
!           its initializing mode
!
!     Coefficients for the transfer of energy between
!     Partially cloudy layers:
      REAL  (Real64) ::                                                 &
     &    CLOUD_OVERLAP(ND_PROFILE, ID_CT-1: ND_LAYER                   &
     &      , ND_OVERLAP_COEFF)                                         &
!           Coefficients defining overlapping options for clouds:
!           these also depend on the solver selected.
     &  , W_FREE(ND_PROFILE, ND_LAYER)
!           Clear-sky fraction
!
!     Cloud geometry
      INTEGER                                                           &
     &    N_COLUMN_CLD(ND_PROFILE)                                      &
!           Number of columns in each profile (including those of
!           zero width)
     &  , N_COLUMN_SLV(ND_PROFILE)                                      &
!           Number of columns to be solved in each profile
     &  , LIST_COLUMN_SLV(ND_PROFILE, ND_COLUMN)                        &
!           List of columns requiring an actual solution
     &  , I_CLM_LYR_CHN(ND_PROFILE, ND_COLUMN)                          &
!           Layer in the current column to change
     &  , I_CLM_CLD_TYP(ND_PROFILE, ND_COLUMN)
!           Type of cloud to introduce in the changed layer
      REAL  (Real64) ::                                                 &
     &    AREA_COLUMN(ND_PROFILE, ND_COLUMN)
!           Areas of columns
#if defined(EXPCLOP)
      REAL  (Real64) ::                                                 &
     &   CLOUD_EDGE(ND_PROFILE, ID_CT: ND_LAYER, 2*ND_CLOUD_ELEM)
!           Array of edges of the cloud elements
#endif
!
!     Local variables for tiled fluxes:
      REAL  (Real64) ::                                                 &
     &    PLANCK_FLUX_TILE(ND_POINT_TILE, ND_TILE)
!           Local Planckian fluxes on surface tiles
!
!
!     Functions called:
      LOGICAL, EXTERNAL :: L_CLOUD_DENSITY
!           Flag for calculation of atmospheric densities for clouds
!
!
!
!
!     Initial determination of flags and switches:
!
      IF (I_ANGULAR_INTEGRATION == IP_TWO_STREAM) THEN
!
!       Only one term in the phase function is required.
        N_ORDER_PHASE=1
!
        L_SOLAR_PHF=.FALSE.
        L_RESCALE_SOLAR_PHF=.FALSE.
!
      ELSE IF (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC) THEN

!        RHO_ALB(:,:,:)=0.0
!        RHO_ALB_TILE(:,:,:,:)=0.0
!
!       Set limits on ranges of harmonics and set pointers to arrays.
! DEPENDS ON: set_truncation
        CALL SET_TRUNCATION(IERR                                        &
     &    , I_TRUNCATION, LS_GLOBAL_TRUNC                               &
     &    , LS_MAX_ORDER, LS_LOCAL_TRUNC                                &
     &    , MS_MIN, MS_MAX, MS_TRUNC                                    &
     &    , IA_SPH_MM, N_ORDER_PHASE                                    &
     &    , ND_MAX_ORDER                                                &
     &    )
!
!       Determine whether special treatment of the solar
!       beam is required.
        L_SOLAR_PHF=(ISOLIR == IP_SOLAR).AND.                           &
     &              (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER)
        L_RESCALE_SOLAR_PHF=L_RESCALE.AND.L_SOLAR_PHF
!       Calculate the solar scattering angles if treating the
!       solar beam separately.
        IF (L_SOLAR_PHF) THEN
! DEPENDS ON: sol_scat_cos
          CALL SOL_SCAT_COS(N_PROFILE, N_DIRECTION                      &
     &      , ZEN_0, DIRECTION, COS_SOL_VIEW                            &
     &      , ND_PROFILE, ND_DIRECTION)
        ENDIF
!
!       Calculate Clebsch-Gordan coefficients once and for all.
! DEPENDS ON: calc_cg_coeff
        CALL CALC_CG_COEFF(LS_MAX_ORDER                                 &
     &    , IA_SPH_MM, MS_MIN, MS_TRUNC                                 &
     &    , CG_COEFF                                                    &
     &    , ND_MAX_ORDER, ND_SPH_COEFF)
!
!       Calculate spherical harmonics at polar angles of pi/2 for
!       use in Marshak's boundary conditions.
! DEPENDS ON: calc_uplm_zero
        CALL CALC_UPLM_ZERO(MS_MIN, MS_MAX, IA_SPH_MM                   &
     &    , LS_LOCAL_TRUNC, UPLM_ZERO                                   &
     &    , ND_MAX_ORDER, ND_SPH_COEFF)
!
        IF (ISOLIR == IP_SOLAR) THEN
!         Calculate the spherical harmonics of the solar direction.
! DEPENDS ON: calc_uplm_sol
          CALL CALC_UPLM_SOL(N_PROFILE, MS_MIN, MS_MAX, IA_SPH_MM       &
     &      , LS_LOCAL_TRUNC, ZEN_0, UPLM_SOL                           &
     &      , ND_PROFILE, ND_MAX_ORDER, ND_SPH_COEFF)
        ENDIF
!
        IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
!         Calcuate some arrays of terms for the BRDF.
! DEPENDS ON: calc_brdf
          CALL CALC_BRDF(ISOLIR, MS_MIN, MS_MAX, IA_SPH_MM              &
     &      , UPLM_SOL, UPLM_ZERO                                       &
     &      , N_BRDF_BASIS_FNC, LS_BRDF_TRUNC, F_BRDF                   &
     &      , N_PROFILE, N_DIRECTION, DIRECTION                         &
     &      , BRDF_SOL, BRDF_HEMI                                       &
     &      , ND_PROFILE, ND_RADIANCE_PROFILE, ND_DIRECTION             &
     &      , ND_MAX_ORDER, ND_SPH_COEFF                                &
     &      , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC)
        ENDIF
!
!       For the calculation of equivalent extinction in the IR
!       we need the diffuse albedo for each basis function.
        L_DIFF_ALB=.FALSE.
        DO I_BAND=1, N_BAND
          L_DIFF_ALB=L_DIFF_ALB.OR.                                     &
     &      (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_K_EQV)
        ENDDO
        IF ( (ISOLIR == IP_INFRA_RED).AND.L_DIFF_ALB ) THEN
! DEPENDS ON: diff_albedo_basis
          CALL DIFF_ALBEDO_BASIS(N_BRDF_BASIS_FNC                       &
     &      , LS_BRDF_TRUNC, F_BRDF                                     &
     &      , UPLM_ZERO(IA_SPH_MM(0))                                   &
     &      , DIFFUSE_ALB_BASIS                                         &
     &      , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC, ND_SPH_COEFF            &
     &      )
        ENDIF
!
!       Determine which layers will be required to give radiances.
! DEPENDS ON: set_rad_layer
        CALL SET_RAD_LAYER(IERR                                         &
     &    , N_LAYER, N_VIEWING_LEVEL, VIEWING_LEVEL                     &
     &    , I_RAD_LAYER, FRAC_RAD_LAYER                                 &
     &    , ND_VIEWING_LEVEL                                            &
     &    )
!
!
      ENDIF
!
!     Set the top level of the profiles. This is currently reatined
!     for historical reasons.
      I_TOP=1
!
!
!
!     Initial calculations for aerosols:
!
!     Set the spectrally independent properties of moist aerosols.
      L_MOIST_AEROSOL=.FALSE.
! CHANGES INTRODUCED BY JCT 02/06/04
      DO J=1, N_AEROSOL
        L_MOIST_AEROSOL=L_MOIST_AEROSOL.OR.                             &
     &    (I_AEROSOL_PARAMETRIZATION(J)                                 &
     &     == IP_AEROSOL_PARAM_MOIST).OR.                               &
     &    (I_AEROSOL_PARAMETRIZATION(J)                                 &
     &     == IP_AEROSOL_PARAM_PHF_MOIST)
      ENDDO
!
      IF (L_MOIST_AEROSOL) THEN
! DEPENDS ON: set_moist_aerosol_properties
        CALL SET_MOIST_AEROSOL_PROPERTIES(IERR                          &
     &    , N_PROFILE, N_LAYER                                          &
     &    , N_AEROSOL, I_AEROSOL_PARAMETRIZATION, NHUMIDITY             &
     &    , GAS_MIX_RATIO(1, 1, INDEX_WATER), T, P, DELTA_HUMIDITY      &
     &    , MEAN_REL_HUMIDITY, I_HUMIDITY_POINTER                       &
     &    , ND_PROFILE, ND_LAYER, ND_AEROSOL_SPECIES                    &
     &    )
        IF (IERR /= I_NORMAL) RETURN
      ENDIF
!
!
!     Check whether the densities will be needed for
!     unparametrized aerosols.
      L_AEROSOL_DENSITY=.FALSE.
      IF (L_AEROSOL) THEN
        DO J=1, N_AEROSOL
          L_AEROSOL_DENSITY=L_AEROSOL_DENSITY.OR.                       &
     &      (I_AEROSOL_PARAMETRIZATION(J) ==                            &
     &       IP_AEROSOL_PARAM_MOIST)                                    &
#if defined(UM)
     &       .OR.(I_AEROSOL_PARAMETRIZATION(J) ==                       &
     &       IP_AEROSOL_UNPARAMETRIZED)
#endif
        ENDDO
      ENDIF
!
!
!
!     Initial calculations for clouds:
!
      IF (L_CLOUD) THEN
!
!       Set pointers to the types of cloud.
! DEPENDS ON: set_cloud_pointer
        CALL SET_CLOUD_POINTER(IERR                                     &
     &    , N_CONDENSED, TYPE_CONDENSED, I_CLOUD_REPRESENTATION         &
     &    , L_DROP, L_ICE                                               &
     &    , I_PHASE_CMP, I_CLOUD_TYPE, L_CLOUD_CMP                      &
     &    , ND_CLOUD_COMPONENT                                          &
     &    )
        IF (IERR /= I_NORMAL) RETURN
!
!
!       Set the geometry of the clouds.
! DEPENDS ON: set_cloud_geometry
        CALL SET_CLOUD_GEOMETRY(N_PROFILE, N_LAYER                      &
     &    , L_GLOBAL_CLOUD_TOP, N_GLOBAL_CLOUD_TOP, W_CLOUD             &
     &    , N_CLOUD_TOP, N_CLOUD_PROFILE, I_CLOUD_PROFILE               &
     &    , ND_PROFILE, ND_LAYER, ID_CT                                 &
     &    )

!       Scale the condensed water contents to simulate
!       inhomogeneities in the clouds.
        IF ( L_INHOM_CLOUD ) THEN
          DO I=1,ND_CLOUD_COMPONENT
            CONDENSED_MIX_RATIO(1:N_PROFILE,N_CLOUD_TOP:N_LAYER,I)      &
     &       =INHOM_CLOUD(I) *                                          &
     &        CONDENSED_MIX_RATIO(1:N_PROFILE,N_CLOUD_TOP:N_LAYER,I)
          ENDDO
        ENDIF
!
        K_CLR=1
        IF ( (I_CLOUD == IP_CLOUD_TRIPLE).OR.                           &
     &       (I_CLOUD == IP_CLOUD_PART_CORR_CNV) ) THEN
!         Aggregate clouds into regions for solving.
!         Three regions are used with this option. Additionally,
!         flag the clear-sky region.
          N_REGION=3
          TYPE_REGION(1)=IP_REGION_CLEAR
          TYPE_REGION(2)=IP_REGION_STRAT
          TYPE_REGION(3)=IP_REGION_CONV
! DEPENDS ON: aggregate_cloud
          CALL AGGREGATE_CLOUD(IERR                                     &
     &      , N_PROFILE, N_LAYER, N_CLOUD_TOP                           &
     &      , I_CLOUD, I_CLOUD_REPRESENTATION                           &
     &      , N_CLOUD_TYPE, FRAC_CLOUD                                  &
     &      , I_REGION_CLOUD, FRAC_REGION                               &
     &      , ND_PROFILE, ND_LAYER, ND_CLOUD_TYPE, ND_REGION            &
     &      , ID_CT                                                     &
     &      )
        ELSE IF ( (I_CLOUD == IP_CLOUD_MIX_MAX).OR.                     &
     &            (I_CLOUD == IP_CLOUD_MIX_RANDOM).OR.                  &
     &            (I_CLOUD == IP_CLOUD_PART_CORR) ) THEN
!         There will be only one cloudy region.
          N_REGION=2
          TYPE_REGION(1)=IP_REGION_CLEAR
          TYPE_REGION(2)=IP_REGION_STRAT
          DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
              FRAC_REGION(L, I, 2)=1.0E+00_Real64
            ENDDO
          ENDDO
        ENDIF
!
!       Calculate energy transfer coefficients in a mixed column,
!       or split the atmosphere into columns with a column model:
!
        IF ( (I_CLOUD == IP_CLOUD_MIX_MAX).OR.                          &
     &       (I_CLOUD == IP_CLOUD_MIX_RANDOM).OR.                       &
     &       (I_CLOUD == IP_CLOUD_TRIPLE).OR.                           &
     &       (I_CLOUD == IP_CLOUD_PART_CORR).OR.                        &
     &       (I_CLOUD == IP_CLOUD_PART_CORR_CNV) ) THEN
!
! DEPENDS ON: overlap_coupled
          CALL OVERLAP_COUPLED(N_PROFILE, N_LAYER, N_CLOUD_TOP          &
     &      , W_CLOUD, W_FREE, N_REGION, TYPE_REGION, FRAC_REGION, P    &
     &      , I_CLOUD                                                   &
     &      , CLOUD_OVERLAP                                             &
     &      , ND_PROFILE, ND_LAYER, ND_OVERLAP_COEFF, ND_REGION         &
     &      , ID_CT, DP_CORR_STRAT, DP_CORR_CONV, TOT_CLOUD_COVER       &
     &      )
!
        ELSE IF (I_CLOUD == IP_CLOUD_COLUMN_MAX) THEN
!
! DEPENDS ON: cloud_maxcs_split
            CALL CLOUD_MAXCS_SPLIT(IERR, N_PROFILE, N_LAYER, N_CLOUD_TOP&
     &        , W_CLOUD, FRAC_CLOUD                                     &
     &        , N_CLOUD_TYPE                                            &
     &        , N_COLUMN_CLD, N_COLUMN_SLV, LIST_COLUMN_SLV             &
     &        , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN               &
     &        , ND_PROFILE, ND_LAYER, ID_CT, ND_COLUMN, ND_CLOUD_TYPE   &
     &        )
!
        ENDIF
!
      ELSE
!
        N_CLOUD_TOP=N_LAYER+1
!
      ENDIF
!
!
!     Calculate the atmospheric densities:
!
      IF ( L_CONTINUUM                                                  &
     &      .OR.L_AEROSOL_DENSITY                                       &
     &      .OR.(L_CLOUD                                                &
! DEPENDS ON: l_cloud_density
     &      .AND.L_CLOUD_DENSITY(N_CONDENSED, I_PHASE_CMP, L_CLOUD_CMP  &
     &                        , I_CONDENSED_PARAM, ND_CLOUD_COMPONENT   &
     &                        ) ) ) THEN
#if defined(STANDARD)
!       Set the pointer for water vapour to a legal value: this must
!       be done for cases where water vapour is not included in the
!       spectral file, but densities are needed for aerosols.
        I_POINTER_WATER=MAX(INDEX_WATER, 1)
#endif
! DEPENDS ON: calculate_density
        CALL CALCULATE_DENSITY(N_PROFILE, N_LAYER, L_CONTINUUM          &
#if defined(UM)
     &    , GAS_MIX_RATIO(1, 1, INDEX_WATER)                            &
#endif
#if defined(STANDARD)
     &    , GAS_MIX_RATIO(1, 1, I_POINTER_WATER)                        &
#endif
     &    , P, T, I_TOP                                                 &
     &    , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN             &
     &    , ND_PROFILE, ND_LAYER                                        &
     &    )
      ENDIF
!
!
!     Check that there is enough information in the case of spherical
!     harmonics. This check is rather late in the logical order of
!     things, but we had to wait for certain other calculations to be
!     made.
      IF (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC) THEN
! DEPENDS ON: check_phf_term
        CALL CHECK_PHF_TERM(IERR                                        &
     &    , L_AEROSOL, N_AEROSOL, I_AEROSOL_PARAMETRIZATION             &
     &    , N_AEROSOL_PHF_TERM                                          &
#if defined(OBSERVED)
     &    , N_PHASE_TERM_AEROSOL_PRSC                                   &
#endif
     &    , L_CLOUD, N_CONDENSED, I_CONDENSED_PARAM, I_PHASE_CMP        &
     &    , CONDENSED_N_PHF                                             &
#if defined(OBSERVED)
     &    , N_PHASE_TERM_DROP_PRSC, N_PHASE_TERM_ICE_PRSC               &
#endif
     &    , N_ORDER_PHASE, L_HENYEY_GREENSTEIN_PF                       &
     &    , L_RESCALE, N_ORDER_FORWARD                                  &
     &    , L_SOLAR_PHF, N_ORDER_PHASE_SOLAR                            &
     &    , ND_AEROSOL_SPECIES, ND_CLOUD_COMPONENT                      &
     &    )
        IF (IERR /= I_NORMAL) RETURN
      ENDIF
!
!
!
!
!
!     Solve the equation of transfer in each band and
!     increment the fluxes.
!
      DO I_BAND=I_FIRST_BAND, I_LAST_BAND
!
!       Set the flag to initialize the diagnostic arrays.
        IF (I_BAND == I_FIRST_BAND) THEN
          L_INITIAL=.TRUE.
        ELSE
          L_INITIAL=(MAP_CHANNEL(I_BAND) >  MAP_CHANNEL(I_BAND-1))
        ENDIF
!
!
!       Determine whether gaseous absorption is included in this band.
        IF ( (L_GAS).AND.(N_BAND_ABSORB(I_BAND) >  0) ) THEN
!
!         Note: I_GAS_BAND is used extensively below since nested
!         array elements in a subroutine call (see later) can
!         confuse some compilers.
!
!         Normally the number of gases in the calculation will be
!         as in the spectral file, but particular options may result
!         in the omission of some gases.
!
          N_GAS=N_BAND_ABSORB(I_BAND)
!
          IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_SINGLE) THEN
!
!           There will be no gaseous absorption in this band
!           unless the selected gas appears.
            N_GAS=0
!
            DO I=1, N_BAND_ABSORB(I_BAND)
              IF (INDEX_ABSORB(I, I_BAND) == I_GAS) N_GAS=1
            ENDDO
!
          ENDIF
!
!
          IF (N_GAS >  0) THEN
!
!           Set the flag for gaseous absorption in the band.
            L_GAS_BAND=.TRUE.
!
            DO J=1, N_GAS
!
              I_GAS_BAND=INDEX_ABSORB(J, I_BAND)
!
!             Reset the pointer if there is just one gas.
!
              IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_SINGLE)           &
     &          THEN
!               Only the selected gas is active in the band.
                I_GAS_BAND=I_GAS
!
              ENDIF
!
              IF (I_SCALE_ESFT(I_BAND, I_GAS_BAND)                      &
     &            == IP_SCALE_BAND) THEN
!               Rescale the amount of gas for this band now.
! DEPENDS ON: scale_absorb
                CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER              &
     &            , GAS_MIX_RATIO(1, 1, I_GAS_BAND), P, T               &
     &            , I_TOP                                               &
     &            , GAS_FRAC_RESCALED(1, 1, I_GAS_BAND)                 &
     &            , I_SCALE_FNC(I_BAND, I_GAS_BAND)                     &
     &            , P_REFERENCE(I_GAS_BAND, I_BAND)                     &
     &            , T_REFERENCE(I_GAS_BAND, I_BAND)                     &
     &            , SCALE_VECTOR(1, 1, I_BAND, I_GAS_BAND)              &
     &            , L_DOPPLER(I_GAS_BAND)                               &
     &            , DOPPLER_CORRECTION(I_GAS_BAND)                      &
     &            , ND_PROFILE, ND_LAYER                                &
     &            , ND_SCALE_VARIABLE                                   &
     &            )
                IF (IERR /= I_NORMAL) RETURN
!
              ELSE IF (I_SCALE_ESFT(I_BAND, I_GAS_BAND)                 &
     &            == IP_SCALE_NULL) THEN
!               Copy across the unscaled array.
                DO I=I_TOP, N_LAYER
                  DO L=1, N_PROFILE
                    GAS_FRAC_RESCALED(L, I, I_GAS_BAND)                 &
     &                =GAS_MIX_RATIO(L, I, I_GAS_BAND)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ELSE
            L_GAS_BAND=.FALSE.
          ENDIF
!
        ELSE
          L_GAS_BAND=.FALSE.
        ENDIF
!
!
!
!       Rescale amounts of continua.
!
        IF (L_CONTINUUM) THEN
          N_CONTINUUM=N_BAND_CONTINUUM(I_BAND)
          DO I=1, N_CONTINUUM
            I_CONTINUUM_POINTER(I)=INDEX_CONTINUUM(I_BAND, I)
            I_CONTINUUM=I_CONTINUUM_POINTER(I)
            K_CONTINUUM_MONO(I_CONTINUUM)                               &
     &        =K_CONTINUUM(I_BAND, I_CONTINUUM)
! DEPENDS ON: rescale_continuum
            CALL RESCALE_CONTINUUM(N_PROFILE, N_LAYER, I_CONTINUUM      &
     &        , P, T, I_TOP                                             &
     &        , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN         &
     &        , GAS_MIX_RATIO(1, 1, INDEX_WATER)                        &
     &        , AMOUNT_CONTINUUM(1, 1, I_CONTINUUM)                     &
     &        , I_SCALE_FNC_CONT(I_BAND, I_CONTINUUM)                   &
     &        , P_REF_CONTINUUM(I_CONTINUUM, I_BAND)                    &
     &        , T_REF_CONTINUUM(I_CONTINUUM, I_BAND)                    &
     &        , SCALE_CONTINUUM(1, I_BAND, I_CONTINUUM)                 &
     &        , ND_PROFILE, ND_LAYER                                    &
     &        , ND_SCALE_VARIABLE                                       &
     &        )
          ENDDO
        ENDIF
!
!
!
!       Calculate the grey extinction within the band.
!
! DEPENDS ON: grey_opt_prop
        CALL GREY_OPT_PROP(IERR                                         &
     &    , N_PROFILE, N_LAYER, P, T, DENSITY                           &
     &    , N_ORDER_PHASE, L_RESCALE, N_ORDER_FORWARD                   &
     &    , L_HENYEY_GREENSTEIN_PF, L_SOLAR_PHF, N_ORDER_PHASE_SOLAR    &
     &    , N_DIRECTION, COS_SOL_VIEW                                   &
     &    , L_RAYLEIGH, RAYLEIGH_COEFFICIENT(I_BAND)                    &
     &    , L_CONTINUUM, N_CONTINUUM, I_CONTINUUM_POINTER               &
     &    , K_CONTINUUM_MONO, AMOUNT_CONTINUUM                          &
     &    , L_AEROSOL, N_AEROSOL, AEROSOL_MIX_RATIO                     &
     &    , I_AEROSOL_PARAMETRIZATION                                   &
     &    , I_HUMIDITY_POINTER, HUMIDITIES, DELTA_HUMIDITY              &
     &    , MEAN_REL_HUMIDITY                                           &
     &    , AEROSOL_ABSORPTION(1, 1, I_BAND)                            &
     &    , AEROSOL_SCATTERING(1, 1, I_BAND)                            &
     &    , AEROSOL_PHASE_FNC(1, 1, 1, I_BAND)                          &
#if defined(OBSERVED)
     &    , N_OPT_LEVEL_AEROSOL_PRSC, AEROSOL_PRESSURE_PRSC             &
     &    , AEROSOL_ABSORPTION_PRSC(1, 1, 1, I_BAND)                    &
     &    , AEROSOL_SCATTERING_PRSC(1, 1, 1, I_BAND)                    &
     &    , AEROSOL_PHASE_FNC_PRSC(1, 1, 1, 1, I_BAND)                  &
#endif
     &    , L_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE, N_CLOUD_TOP      &
     &    , N_CONDENSED, L_CLOUD_CMP, I_PHASE_CMP                       &
     &    , I_CONDENSED_PARAM, CONDENSED_N_PHF                          &
     &    , CONDENSED_PARAM_LIST(1, 1, I_BAND)                          &
     &    , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                     &
     &    , N_CLOUD_TYPE, I_CLOUD_TYPE                                  &
#if defined(OBSERVED)
     &    , N_OPT_LEVEL_DROP_PRSC, DROP_PRESSURE_PRSC                   &
     &    , DROP_ABSORPTION_PRSC(1, 1, I_BAND)                          &
     &    , DROP_SCATTERING_PRSC(1, 1, I_BAND)                          &
     &    , DROP_PHASE_FNC_PRSC(1, 1, 1, I_BAND)                        &
     &    , N_OPT_LEVEL_ICE_PRSC, ICE_PRESSURE_PRSC                     &
     &    , ICE_ABSORPTION_PRSC(1, 1, I_BAND)                           &
     &    , ICE_SCATTERING_PRSC(1, 1, I_BAND)                           &
     &    , ICE_PHASE_FNC_PRSC(1, 1, 1, I_BAND)                         &
#endif
     &    , K_GREY_TOT_CLR, K_EXT_SCAT_CLR, PHASE_FNC_CLR               &
     &    , FORWARD_SCATTER_CLR, FORWARD_SOLAR_CLR, PHASE_FNC_SOLAR_CLR &
     &    , K_GREY_TOT, K_EXT_SCAT, PHASE_FNC, FORWARD_SCATTER          &
     &    , FORWARD_SOLAR, PHASE_FNC_SOLAR                              &
     &    , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER                   &
     &    , ND_LAYER_CLR, ID_CT                                         &
     &    , ND_CONTINUUM, ND_AEROSOL_SPECIES, ND_HUMIDITIES             &
     &    , ND_CLOUD_PARAMETER, ND_CLOUD_COMPONENT, ND_CLOUD_TYPE       &
     &    , ND_PHASE_TERM, ND_MAX_ORDER, ND_DIRECTION                   &
#if defined(OBSERVED)
     &    , ND_PROFILE_AEROSOL_PRSC, ND_PROFILE_CLOUD_PRSC              &
     &    , ND_OPT_LEVEL_AEROSOL_PRSC, ND_OPT_LEVEL_CLOUD_PRSC          &
#endif
     &    )
        IF (IERR /= I_NORMAL) RETURN
!
!
        IF ( (I_ANGULAR_INTEGRATION == IP_TWO_STREAM).OR.               &
     &       (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC) ) THEN
!
!         Rescale the phase function and calculate the scattering
!         fractions. (These are grey and may be calculated outside
!         a loop over gases).
!
          IF (L_RESCALE) THEN
!
!           Rescale clear-sky phase function:
!
!           The section above clouds.
! DEPENDS ON: rescale_phase_fnc
            CALL RESCALE_PHASE_FNC(N_PROFILE, 1, N_CLOUD_TOP-1          &
     &        , N_DIRECTION, COS_SOL_VIEW                               &
     &        , N_ORDER_PHASE                                           &
     &        , PHASE_FNC_CLR, FORWARD_SCATTER_CLR                      &
     &        , FORWARD_SOLAR_CLR                                       &
     &        , L_RESCALE_SOLAR_PHF, N_ORDER_PHASE_SOLAR                &
     &        , PHASE_FNC_SOLAR_CLR                                     &
     &        , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER_CLR, 1        &
     &        , ND_DIRECTION, ND_MAX_ORDER                              &
     &        )
!           The section including clouds.
! DEPENDS ON: rescale_phase_fnc
            CALL RESCALE_PHASE_FNC(N_PROFILE, N_CLOUD_TOP               &
     &        , N_LAYER, N_DIRECTION, COS_SOL_VIEW                      &
     &        , N_ORDER_PHASE                                           &
     &        , PHASE_FNC(1, ID_CT, 1, 0)                               &
     &        , FORWARD_SCATTER(1, ID_CT, 0)                            &
     &        , FORWARD_SOLAR(1, ID_CT, 0)                              &
     &        , L_RESCALE_SOLAR_PHF, N_ORDER_PHASE_SOLAR                &
     &        , PHASE_FNC_SOLAR(1, ID_CT, 1, 0)                         &
     &        , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER, ID_CT        &
     &        , ND_DIRECTION, ND_MAX_ORDER                              &
     &        )
!
!
            IF (L_CLOUD) THEN
!
!             Rescale cloudy phase functions:
!
              DO K=1, N_CLOUD_TYPE
! DEPENDS ON: rescale_phase_fnc
                CALL RESCALE_PHASE_FNC(N_PROFILE, N_CLOUD_TOP           &
     &            , N_LAYER, N_DIRECTION, COS_SOL_VIEW                  &
     &            , N_ORDER_PHASE                                       &
     &            , PHASE_FNC(1, ID_CT, 1, K)                           &
     &            , FORWARD_SCATTER(1, ID_CT, K)                        &
     &            , FORWARD_SOLAR(1, ID_CT, K)                          &
     &            , L_RESCALE_SOLAR_PHF, N_ORDER_PHASE_SOLAR            &
     &            , PHASE_FNC_SOLAR(1, ID_CT, 1, K)                     &
     &            , ND_PROFILE, ND_RADIANCE_PROFILE, ND_LAYER, ID_CT    &
     &            , ND_DIRECTION, ND_MAX_ORDER                          &
     &            )
              ENDDO
!
            ENDIF
!
          ENDIF
!
        ENDIF
!
!
!
!
!       Preliminary calculations for source terms:
!
        IF (ISOLIR == IP_SOLAR) THEN
!         Convert normalized band fluxes to actual energy fluxes.
          DO L=1, N_PROFILE
            SOLAR_IRRAD_BAND(L)=SOLAR_IRRAD(L)                          &
     &        *SOLAR_FLUX_BAND(I_BAND)
          ENDDO
!
        ELSE IF (ISOLIR == IP_INFRA_RED) THEN
!
!         Calculate the change in the thermal source function
!         across each layer for the infra-red part of the spectrum.
!
! DEPENDS ON: diff_planck_source
          CALL DIFF_PLANCK_SOURCE(N_PROFILE, N_LAYER                    &
     &      , N_DEG_FIT, THERMAL_COEFFICIENT(0, I_BAND)                 &
     &      , T_REF_PLANCK, T_LEVEL, T_GROUND                           &
     &      , PLANCK_FLUX_BAND, DIFF_PLANCK_BAND                        &
     &      , PLANCK_FLUX_GROUND                                        &
     &      , L_IR_SOURCE_QUAD, T, DIFF_PLANCK_BAND_2                   &
     &      , I_ANGULAR_INTEGRATION                                     &
     &      , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER              &
     &      , PLANCK_RADIANCE_BAND                                      &
     &      , L_TILE, N_POINT_TILE, N_TILE, LIST_TILE                   &
     &      , FRAC_TILE, T_TILE, PLANCK_FLUX_TILE                       &
     &      , ND_PROFILE, ND_LAYER, ND_THERMAL_COEFF                    &
     &      , ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL                     &
     &      , ND_POINT_TILE, ND_TILE                                    &
     &      )
!
        ENDIF
!
!
!
!
!
!
!       Call a solver appropriate to the presence of gases and
!       the overlap assumed:
!
        IF (.NOT.L_GAS_BAND) THEN
!
!         There is no gaseous absorption. Solve for the
!         radiances directly.
!
! DEPENDS ON: solve_band_without_gas
          CALL SOLVE_BAND_WITHOUT_GAS(IERR                              &
!                     Atmospheric properties
     &      , N_PROFILE, N_LAYER, D_MASS                                &
!                     Angular integration
     &      , I_ANGULAR_INTEGRATION, I_2STREAM                          &
     &      , N_ORDER_PHASE, L_RESCALE, N_ORDER_GAUSS                   &
     &      , MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC              &
     &      , ACCURACY_ADAPTIVE, EULER_FACTOR, I_SPH_ALGORITHM          &
     &      , I_SPH_MODE                                                &
!                     Precalculated angular arrays
     &      , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                  &
!                     Treatment of scattering
     &      , I_SCATTER_METHOD(I_BAND)                                  &
!                     Options for solver
     &      , I_SOLVER                                                  &
!                     Spectral region
     &      , ISOLIR                                                    &
!                     Solar properties
     &      , ZEN_0, SOLAR_IRRAD_BAND                                   &
!                     Infra-red properties
     &      , PLANCK_FLUX_BAND(1, 0), PLANCK_FLUX_BAND(1, N_LAYER)      &
     &      , DIFF_PLANCK_BAND, L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2    &
!                     Surface properties
     &      , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB(1, 1, I_BAND)    &
     &      , F_BRDF, BRDF_SOL, BRDF_HEMI                               &
     &      , PLANCK_FLUX_GROUND                                        &
!                       Tiling of the surface
     &      , L_TILE, N_POINT_TILE, N_TILE, LIST_TILE                   &
     &      , RHO_ALB_TILE(1, 1, 1, I_BAND), PLANCK_FLUX_TILE           &
!                       Optical Properties
     &      , K_GREY_TOT_CLR, K_EXT_SCAT_CLR                            &
     &      , PHASE_FNC_CLR, FORWARD_SCATTER_CLR, PHASE_FNC_SOLAR_CLR   &
     &      , K_GREY_TOT, K_EXT_SCAT                                    &
     &      , PHASE_FNC, FORWARD_SCATTER, PHASE_FNC_SOLAR               &
!                     Cloudy properties
     &      , L_CLOUD, I_CLOUD                                          &
!                     Cloudy geometry
     &      , N_CLOUD_TOP                                               &
     &      , N_CLOUD_TYPE, FRAC_CLOUD                                  &
     &      , N_REGION, K_CLR, I_REGION_CLOUD, FRAC_REGION              &
     &      , W_FREE, W_CLOUD, CLOUD_OVERLAP                            &
     &      , N_COLUMN_SLV, LIST_COLUMN_SLV                             &
     &      , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN                 &
!                      Levels for calculating radiances
     &      , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER              &
!                       Viewing Geometry
     &      , N_DIRECTION, DIRECTION                                    &
!                     Weighting factor for the band
     &      , WEIGHT_BAND(I_BAND), L_INITIAL                            &
!                     Calculated fluxes
     &      , FLUX_DIRECT(1, 0, MAP_CHANNEL(I_BAND))                    &
     &      , FLUX_DOWN(1, 0, MAP_CHANNEL(I_BAND))                      &
     &      , FLUX_UP(1, 0, MAP_CHANNEL(I_BAND))                        &
!                     Radiances
     &      , I_DIRECT, RADIANCE(1, 1, 1, MAP_CHANNEL(I_BAND))          &
!                     Rate of photolysis
     &      , PHOTOLYSIS(1, 1, MAP_CHANNEL(I_BAND))                     &
!                     Flags for clear-sky fluxes
     &      , L_CLEAR, I_SOLVER_CLEAR                                   &
!                     Calculated clear-sky fluxes
     &      , FLUX_DIRECT_CLEAR(1, 0, MAP_CHANNEL(I_BAND))              &
     &      , FLUX_DOWN_CLEAR(1, 0, MAP_CHANNEL(I_BAND))                &
     &      , FLUX_UP_CLEAR(1, 0, MAP_CHANNEL(I_BAND))                  &
!                       Tiled Surface Fluxes
     &      , FLUX_UP_TILE(1, 1, MAP_CHANNEL(I_BAND))                   &
     &      , FLUX_UP_BLUE_TILE(1, 1, MAP_CHANNEL(I_BAND))              &
!                       Special Surface Fluxes
     &      , L_BLUE_FLUX_SURF, WEIGHT_BLUE(I_BAND)                     &
     &      , FLUX_DIRECT_BLUE_SURF                                     &
     &      , FLUX_DOWN_BLUE_SURF, FLUX_UP_BLUE_SURF                    &
!                     Dimensions of arrays
     &      , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN      &
     &      , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE        &
     &      , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF                &
     &      , ND_MAX_ORDER, ND_SPH_COEFF                                &
     &      , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC                          &
     &      , ND_VIEWING_LEVEL, ND_DIRECTION                            &
     &      , ND_SOURCE_COEFF, ND_POINT_TILE, ND_TILE                   &
     &      )
          IF (IERR /= I_NORMAL) RETURN
!
!
        ELSE
!
!         Gases are included.
!
!         Treat the gaseous overlaps as directed by
!         the overlap switch.
!
          IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_SINGLE) THEN
!
! DEPENDS ON: solve_band_one_gas
            CALL SOLVE_BAND_ONE_GAS(IERR                                &
!                     Atmospheric properties
     &        , N_PROFILE, N_LAYER, I_TOP, P, T, D_MASS                 &
!                     Angular integration
     &        , I_ANGULAR_INTEGRATION, I_2STREAM                        &
     &        , N_ORDER_PHASE, L_RESCALE, N_ORDER_GAUSS                 &
     &        , MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC            &
     &        , ACCURACY_ADAPTIVE, EULER_FACTOR                         &
     &        , I_SPH_ALGORITHM, I_SPH_MODE                             &
!                     Precalculated angular arrays
     &        , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                &
!                     Treatment of scattering
     &        , I_SCATTER_METHOD(I_BAND)                                &
!                     Options for solver
     &        , I_SOLVER                                                &
!                     Gaseous properties
     &        , I_BAND, I_GAS                                           &
     &        , I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC                  &
     &        , K_ESFT, W_ESFT, SCALE_VECTOR                            &
     &        , P_REFERENCE, T_REFERENCE                                &
     &        , GAS_MIX_RATIO, GAS_FRAC_RESCALED                        &
     &        , L_DOPPLER, DOPPLER_CORRECTION                           &
!                     Spectral region
     &        , ISOLIR                                                  &
!                     Solar properties
     &        , ZEN_0, SOLAR_IRRAD_BAND                                 &
!                     Infra-red properties
     &        , PLANCK_FLUX_BAND(1, 0)                                  &
     &        , PLANCK_FLUX_BAND(1, N_LAYER)                            &
     &        , DIFF_PLANCK_BAND                                        &
     &        , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                    &
!                     Surface properties
     &        , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB(1, 1, I_BAND)  &
     &        , F_BRDF, BRDF_SOL, BRDF_HEMI                             &
     &        , PLANCK_FLUX_GROUND                                      &
!                       Tiling of the surface
     &        , L_TILE, N_POINT_TILE, N_TILE, LIST_TILE                 &
     &        , RHO_ALB_TILE(1, 1, 1, I_BAND)                           &
     &        , PLANCK_FLUX_TILE                                        &
!                       Optical Properties
     &        , K_GREY_TOT_CLR, K_EXT_SCAT_CLR                          &
     &        , PHASE_FNC_CLR, FORWARD_SCATTER_CLR, PHASE_FNC_SOLAR_CLR &
     &        , K_GREY_TOT, K_EXT_SCAT, PHASE_FNC, FORWARD_SCATTER      &
     &        , PHASE_FNC_SOLAR                                         &
!                     Cloudy properties
     &        , L_CLOUD, I_CLOUD                                        &
!                     Cloud geometry
     &        , N_CLOUD_TOP                                             &
     &        , N_CLOUD_TYPE, FRAC_CLOUD                                &
     &        , N_REGION, K_CLR, I_REGION_CLOUD, FRAC_REGION            &
     &        , W_FREE, W_CLOUD, CLOUD_OVERLAP                          &
     &        , N_COLUMN_SLV, LIST_COLUMN_SLV                           &
     &        , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN               &
!                       Levels for calculating radiances
     &        , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER            &
!                       Viewing Geometry
     &        , N_DIRECTION, DIRECTION                                  &
!                     Weighting factor for the band
     &        , WEIGHT_BAND(I_BAND), L_INITIAL                          &
!                     Fluxes calculated
     &        , FLUX_DIRECT(1, 0, MAP_CHANNEL(I_BAND))                  &
     &        , FLUX_DOWN(1, 0, MAP_CHANNEL(I_BAND))                    &
     &        , FLUX_UP(1, 0, MAP_CHANNEL(I_BAND))                      &
!                     Radiances
     &        , I_DIRECT, RADIANCE(1, 1, 1, MAP_CHANNEL(I_BAND))        &
!                     Rate of photolysis
     &        , PHOTOLYSIS(1, 1, MAP_CHANNEL(I_BAND))                   &
!                     Flags for clear-sky calculations
     &        , L_CLEAR, I_SOLVER_CLEAR                                 &
!                     Clear-sky fluxes
     &        , FLUX_DIRECT_CLEAR(1, 0, MAP_CHANNEL(I_BAND))            &
     &        , FLUX_DOWN_CLEAR(1, 0, MAP_CHANNEL(I_BAND))              &
     &        , FLUX_UP_CLEAR(1, 0, MAP_CHANNEL(I_BAND))                &
!                       Tiled Surface Fluxes
     &        , FLUX_UP_TILE(1, 1, MAP_CHANNEL(I_BAND))                 &
     &        , FLUX_UP_BLUE_TILE(1, 1, MAP_CHANNEL(I_BAND))            &
!                       Special Surface Fluxes
     &        , L_BLUE_FLUX_SURF, WEIGHT_BLUE(I_BAND)                   &
     &        , FLUX_DIRECT_BLUE_SURF                                   &
     &        , FLUX_DOWN_BLUE_SURF, FLUX_UP_BLUE_SURF                  &
!                     Dimensions of arrays
     &        , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN    &
     &        , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE      &
     &        , ND_BAND, ND_SPECIES                                     &
     &        , ND_ESFT_TERM, ND_SCALE_VARIABLE                         &
     &        , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF              &
     &        , ND_MAX_ORDER, ND_SPH_COEFF                              &
     &        , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC                        &
     &        , ND_VIEWING_LEVEL, ND_DIRECTION                          &
     &        , ND_SOURCE_COEFF, ND_POINT_TILE, ND_TILE                 &
     &        )
!
          ELSE IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_RANDOM) THEN
!
! DEPENDS ON: solve_band_random_overlap
            CALL SOLVE_BAND_RANDOM_OVERLAP(IERR                         &
!                     Atmospheric properties
     &        , N_PROFILE, N_LAYER, I_TOP, P, T, D_MASS                 &
!                     Angular integration
     &        , I_ANGULAR_INTEGRATION, I_2STREAM                        &
     &        , N_ORDER_PHASE, L_RESCALE, N_ORDER_GAUSS                 &
     &        , MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC            &
     &        , ACCURACY_ADAPTIVE, EULER_FACTOR                         &
     &        , I_SPH_ALGORITHM, I_SPH_MODE                             &
!                     Precalculated angular arrays
     &        , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                &
!                     Treatment of scattering
     &        , I_SCATTER_METHOD(I_BAND)                                &
!                     Options for solver
     &        , I_SOLVER                                                &
!                     Gaseous properties
     &        , I_BAND, N_GAS                                           &
     &        , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC    &
     &        , K_ESFT, W_ESFT, SCALE_VECTOR                            &
     &        , P_REFERENCE, T_REFERENCE                                &
     &        , GAS_MIX_RATIO, GAS_FRAC_RESCALED                        &
     &        , L_DOPPLER, DOPPLER_CORRECTION                           &
!                     Spectral region
     &        , ISOLIR                                                  &
!                     Solar properties
     &        , ZEN_0, SOLAR_IRRAD_BAND                                 &
!                     Infra-red properties
     &        , PLANCK_FLUX_BAND(1, 0)                                  &
     &        , PLANCK_FLUX_BAND(1, N_LAYER)                            &
     &        , DIFF_PLANCK_BAND                                        &
     &        , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                    &
!                     Surface properties
     &        , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB(1, 1, I_BAND)  &
     &        , F_BRDF, BRDF_SOL, BRDF_HEMI                             &
     &        , PLANCK_FLUX_GROUND                                      &
!                       Tiling of the surface
     &        , L_TILE, N_POINT_TILE, N_TILE, LIST_TILE                 &
     &        , RHO_ALB_TILE(1, 1, 1, I_BAND)                           &
     &        , PLANCK_FLUX_TILE                                        &
!                       Optical Properties
     &        , K_GREY_TOT_CLR, K_EXT_SCAT_CLR, PHASE_FNC_CLR           &
     &        , FORWARD_SCATTER_CLR, PHASE_FNC_SOLAR_CLR                &
     &        , K_GREY_TOT, K_EXT_SCAT, PHASE_FNC, FORWARD_SCATTER      &
     &        , PHASE_FNC_SOLAR                                         &
!                     Cloudy properties
     &        , L_CLOUD, I_CLOUD                                        &
!                     Cloud geometry
     &        , N_CLOUD_TOP                                             &
     &        , N_CLOUD_TYPE, FRAC_CLOUD                                &
     &        , N_REGION, K_CLR, I_REGION_CLOUD, FRAC_REGION            &
     &        , W_FREE, W_CLOUD, CLOUD_OVERLAP                          &
     &        , N_COLUMN_SLV, LIST_COLUMN_SLV                           &
     &        , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN               &
!                       Levels for calculating radiances
     &        , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER            &
!                       Viewing Geometry
     &        , N_DIRECTION, DIRECTION                                  &
!                     Weighting factor for the band
     &        , WEIGHT_BAND(I_BAND), L_INITIAL                          &
!                     Fluxes calculated
     &        , FLUX_DIRECT(1, 0, MAP_CHANNEL(I_BAND))                  &
     &        , FLUX_DOWN(1, 0, MAP_CHANNEL(I_BAND))                    &
     &        , FLUX_UP(1, 0, MAP_CHANNEL(I_BAND))                      &
!                     Radiances
     &        , I_DIRECT, RADIANCE(1, 1, 1, MAP_CHANNEL(I_BAND))        &
!                     Rate of photolysis
     &        , PHOTOLYSIS(1, 1, MAP_CHANNEL(I_BAND))                   &
!                     Flags for clear-sky calculations
     &        , L_CLEAR, I_SOLVER_CLEAR                                 &
!                     Clear-sky fluxes
     &        , FLUX_DIRECT_CLEAR(1, 0, MAP_CHANNEL(I_BAND))            &
     &        , FLUX_DOWN_CLEAR(1, 0, MAP_CHANNEL(I_BAND))              &
     &        , FLUX_UP_CLEAR(1, 0, MAP_CHANNEL(I_BAND))                &
!                       Tiled Surface Fluxes
     &        , FLUX_UP_TILE(1, 1, MAP_CHANNEL(I_BAND))                 &
     &        , FLUX_UP_BLUE_TILE(1, 1, MAP_CHANNEL(I_BAND))            &
!                       Special Surface Fluxes
     &        , L_BLUE_FLUX_SURF, WEIGHT_BLUE(I_BAND)                   &
     &        , FLUX_DIRECT_BLUE_SURF                                   &
     &        , FLUX_DOWN_BLUE_SURF, FLUX_UP_BLUE_SURF                  &
!                     Dimensions of arrays
     &        , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN    &
     &        , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE      &
     &        , ND_BAND, ND_SPECIES                                     &
     &        , ND_ESFT_TERM, ND_SCALE_VARIABLE                         &
     &        , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF              &
     &        , ND_MAX_ORDER, ND_SPH_COEFF                              &
     &        , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC, ND_VIEWING_LEVEL      &
     &        , ND_DIRECTION, ND_SOURCE_COEFF, ND_POINT_TILE, ND_TILE   &
     &        )
!
          ELSE IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_K_EQV) THEN
!
! DEPENDS ON: solve_band_k_eqv
            CALL SOLVE_BAND_K_EQV(IERR                                  &
!                     Atmospheric properties
     &        , N_PROFILE, N_LAYER, I_TOP, P, T, D_MASS                 &
!                     Angular integration
     &        , I_ANGULAR_INTEGRATION, I_2STREAM                        &
     &        , N_ORDER_PHASE, L_RESCALE, N_ORDER_GAUSS                 &
     &        , MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC            &
     &        , ACCURACY_ADAPTIVE, EULER_FACTOR                         &
     &        , I_SPH_ALGORITHM, I_SPH_MODE                             &
!                     Precalculated angular arrays
     &        , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                &
!                     Treatment of scattering
     &        , I_SCATTER_METHOD(I_BAND)                                &
!                     Options for solver
     &        , I_SOLVER                                                &
!                     Gaseous properties
     &        , I_BAND, N_GAS                                           &
     &        , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC    &
     &        , K_ESFT, W_ESFT, SCALE_VECTOR                            &
     &        , P_REFERENCE, T_REFERENCE, L_MOD_K_FLUX                  &
     &        , GAS_MIX_RATIO, GAS_FRAC_RESCALED                        &
     &        , L_DOPPLER, DOPPLER_CORRECTION                           &
!                     Spectral region
     &        , ISOLIR                                                  &
!                     Solar properties
     &        , ZEN_0, SOLAR_IRRAD_BAND                                 &
!                     Infra-red properties
     &        , PLANCK_FLUX_BAND                                        &
     &        , DIFF_PLANCK_BAND                                        &
     &        , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                    &
!                     Surface properties
     &        , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB(1, 1, I_BAND)  &
     &        , F_BRDF, BRDF_SOL, BRDF_HEMI                             &
     &        , DIFFUSE_ALB_BASIS                                       &
     &        , PLANCK_FLUX_GROUND                                      &
!                       Tiling of the surface
     &        , L_TILE, N_POINT_TILE, N_TILE, LIST_TILE                 &
     &        , RHO_ALB_TILE(1, 1, 1, I_BAND)                           &
     &        , PLANCK_FLUX_TILE                                        &
!                       Optical Properties
     &        , K_GREY_TOT_CLR, K_EXT_SCAT_CLR, PHASE_FNC_CLR           &
     &        , FORWARD_SCATTER_CLR, PHASE_FNC_SOLAR_CLR                &
     &        , K_GREY_TOT, K_EXT_SCAT, PHASE_FNC, FORWARD_SCATTER      &
     &        , PHASE_FNC_SOLAR                                         &
!                     Cloudy properties
     &        , L_CLOUD, I_CLOUD                                        &
!                     Cloud geometry
     &        , N_CLOUD_TOP                                             &
     &        , N_CLOUD_TYPE, FRAC_CLOUD                                &
     &        , N_REGION, K_CLR, I_REGION_CLOUD, FRAC_REGION            &
     &        , W_FREE, W_CLOUD, CLOUD_OVERLAP                          &
     &        , N_COLUMN_SLV, LIST_COLUMN_SLV                           &
     &        , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN               &
!                       Levels for calculating radiances
     &        , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER            &
!                       Viewing Geometry
     &        , N_DIRECTION, DIRECTION                                  &
!                     Weighting factor for the band
     &        , WEIGHT_BAND(I_BAND), L_INITIAL                          &
!                     Fluxes calculated
     &        , FLUX_DIRECT(1, 0, MAP_CHANNEL(I_BAND))                  &
     &        , FLUX_DOWN(1, 0, MAP_CHANNEL(I_BAND))                    &
     &        , FLUX_UP(1, 0, MAP_CHANNEL(I_BAND))                      &
!                     Radiances
     &        , I_DIRECT, RADIANCE(1, 1, 1, MAP_CHANNEL(I_BAND))        &
!                     Rate of photolysis
     &        , PHOTOLYSIS(1, 1, MAP_CHANNEL(I_BAND))                   &
!                     Flags for clear-sky calculations
     &        , L_CLEAR, I_SOLVER_CLEAR                                 &
!                     Clear-sky fluxes calculated
     &        , FLUX_DIRECT_CLEAR(1, 0, MAP_CHANNEL(I_BAND))            &
     &        , FLUX_DOWN_CLEAR(1, 0, MAP_CHANNEL(I_BAND))              &
     &        , FLUX_UP_CLEAR(1, 0, MAP_CHANNEL(I_BAND))                &
!                       Tiled Surface Fluxes
     &        , FLUX_UP_TILE(1, 1, MAP_CHANNEL(I_BAND))                 &
     &        , FLUX_UP_BLUE_TILE(1, 1, MAP_CHANNEL(I_BAND))            &
!                       Special Surface Fluxes
     &        , L_BLUE_FLUX_SURF, WEIGHT_BLUE(I_BAND)                   &
     &        , FLUX_DIRECT_BLUE_SURF                                   &
     &        , FLUX_DOWN_BLUE_SURF, FLUX_UP_BLUE_SURF                  &
!                     Dimensions of arrays
     &        , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN    &
     &        , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE      &
     &        , ND_BAND, ND_SPECIES                                     &
     &        , ND_ESFT_TERM, ND_SCALE_VARIABLE                         &
     &        , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF              &
     &        , ND_MAX_ORDER, ND_SPH_COEFF                              &
     &        , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC                        &
     &        , ND_VIEWING_LEVEL, ND_DIRECTION                          &
     &        , ND_SOURCE_COEFF, ND_POINT_TILE, ND_TILE                 &
     &        )
!
          ELSE
            WRITE(IU_ERR, '(3(/A))')                                    &
     &        '*** Error: An appropriate gaseous overlap'               &
     &        , 'has not been specified, even though gaseous'           &
     &        , 'absorption is to be included.'
          ENDIF
        ENDIF

        IF (ISOLIR == IP_SOLAR) THEN
!
! Calculate the Diffuse Flux
!

          IF (L_FLUX_DIFFUSE) THEN
            DO I=0,N_LAYER
              DO L=1, N_PROFILE
                FLUX_DIFFUSE(L,I,1)=FLUX_DIFFUSE(L,I,1)                 &
     &                  +(FLUX_DOWN(L,I,1)-FLUX_DIRECT(L,I,1))
              ENDDO
            ENDDO
          ENDIF


!
! Calculate UV-Fluxes
!
          IF (L_UVFLUX_DIRECT) THEN
            DO I=0, N_LAYER
              DO L=1, N_PROFILE
                UV_FLUX_DIRECT(L,I,1)=UV_FLUX_DIRECT(L,I,1)             &
     &                + WEIGHT_UV(I_BAND)*FLUX_DIRECT(L,I,1)        
              ENDDO
            ENDDO       
          ENDIF
          IF (L_UVFLUX_DOWN) THEN 
            DO I=0, N_LAYER
              DO L=1, N_PROFILE
                UV_FLUX_UP(L,I,1)=UV_FLUX_UP(L,I,1)                     &
     &                + WEIGHT_UV(I_BAND)*FLUX_UP(L,I,1)            
              ENDDO
            ENDDO       
          ENDIF
          IF (L_UVFLUX_UP) THEN
            DO I=0, N_LAYER
              DO L=1, N_PROFILE
                UV_FLUX_DOWN(L,I,1)=UV_FLUX_DOWN(L,I,1)                 &
     &                + WEIGHT_UV(I_BAND)*FLUX_DOWN(L,I,1)            
              ENDDO
            ENDDO       
          ENDIF 

        ENDIF      
  
!
!       Make any adjustments to fluxes and radiances to convert
!       to actual values. This is done inside the loop over bands
!       to allow for division of the output fluxes between
!       separate diagnostic bands.


        IF (ISOLIR == IP_INFRA_RED) THEN
! DEPENDS ON: adjust_ir_radiance
          CALL ADJUST_IR_RADIANCE(N_PROFILE, N_LAYER, N_VIEWING_LEVEL   &
     &      , N_DIRECTION, I_ANGULAR_INTEGRATION, I_SPH_MODE            &
     &      , PLANCK_FLUX_BAND, PLANCK_RADIANCE_BAND                    &
     &      , FLUX_DOWN(1, 0, MAP_CHANNEL(I_BAND))                      &
     &      , FLUX_UP(1, 0, MAP_CHANNEL(I_BAND))                        &
     &      , RADIANCE(1, 1, 1, MAP_CHANNEL(I_BAND))                    &
     &      , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE                      &
     &      , ND_LAYER, ND_DIRECTION, ND_VIEWING_LEVEL                  &
     &      )
        ENDIF
!
!
      ENDDO
!
!
      RETURN
      END SUBROUTINE RADIANCE_CALC
#endif
