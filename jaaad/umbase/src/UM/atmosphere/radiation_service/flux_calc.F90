#if defined(A70_1B) || defined(A70_1C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate radiative fluxes.
!
! Method:
!       Properties independent of the spectral bands are set.
!       A loop over bands is then entered. Grey optical properties
!       are set and an appropriate subroutine is called to treat
!       the gaseous overlaps. The final fluxes are assigned.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE FLUX_CALC(IERR                                         &
!                       Logical Flags for Processes
     &   , L_RAYLEIGH, L_AEROSOL, L_GAS, L_CONTINUUM                    &
     &   , L_CLOUD, L_DROP, L_ICE, L_PC2                                &
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT         &
     &   , L_RESCALE, N_ORDER_GAUSS                                     &
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_BAND                                        &
!                       Options for treating clouds
     &   , L_GLOBAL_CLOUD_TOP, N_CLOUD_TOP_GLOBAL                       &
     &   , L_INHOM_CLOUD, INHOM_CLOUD                                   &
!    &   , L_INHOM_CLOUD, INHOM_CLOUD                                   &
!                       Options for Solver
     &   , I_SOLVER                                                     &
!                       General Spectral Properties
     &   , N_BAND, I_FIRST_BAND, I_LAST_BAND, WEIGHT_BAND               &
!                       General Atmospheric Properties
     &   , N_PROFILE, N_LAYER                                           &
     &   , L_LAYER, L_CLOUD_LAYER                                       &
     &   , P, T, T_GROUND, T_SOLID, T_SEA, T_LEVEL, D_MASS              &
!                       Spectral Region
     &   , ISOLIR                                                       &
!                       Solar Fields
     &   , SEC_0, SOLAR_TOA, SOLAR_FLUX_BAND, RAYLEIGH_COEFFICIENT      &
!                       Infra-red Fields
     &   , N_DEG_FIT, THERMAL_COEFFICIENT, T_REF_PLANCK                 &
     &   , L_IR_SOURCE_QUAD                                             &
!                       Gaseous Absorption
     &   , N_ABSORB, I_GAS_OVERLAP, I_GAS                               &
     &   , GAS_MIX_RATIO, N_BAND_ABSORB, INDEX_ABSORB                   &
     &   , I_BAND_ESFT, W_ESFT, K_ESFT, I_SCALE_ESFT                    &
     &   , I_SCALE_FNC, L_WENYI, SCALE_VECTOR                           &
     &   , P_REFERENCE, T_REFERENCE, L_MOD_K_FLUX                       &
!                       Doppler Broadening
     &   , L_DOPPLER, DOPPLER_CORRECTION                                &
!                       Surface Fields
     &   , L_SURFACE, I_SURFACE, I_SPEC_SURFACE, SURFACE_ALBEDO         &
     &   , ALBEDO_FIELD_DIFF, ALBEDO_FIELD_DIR                          &
     &   , N_DIR_ALBEDO_FIT, DIRECT_ALBEDO_PARM                         &
     &   , EMISSIVITY_GROUND, EMISSIVITY_FIELD                          &
!                       Continuum Absorption
     &   , N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER               &
     &   , K_CONTINUUM, I_SCALE_FNC_CONT, SCALE_CONTINUUM               &
     &   , P_REF_CONTINUUM, T_REF_CONTINUUM                             &
!                       Properties of Aerosols
     &   , N_AEROSOL, AEROSOL_MIX_RATIO                                 &
     &   , AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_ASYMMETRY    &
     &   , I_AEROSOL_PARAMETRIZATION, NHUMIDITY, HUMIDITIES             &
     &   , TYPE_AEROSOL, L_USE_CLEARRH                                  &
!                       Aerosol optical depth
     &   , N_AOD_WAVEL                                                  &
     &   , L_AOD_SULPHATE, AOD_SULPHATE, L_AOD_DUST,     AOD_DUST       &
     &   , L_AOD_SEASALT,  AOD_SEASALT,  L_AOD_SOOT,     AOD_SOOT       &
     &   , L_AOD_BIOMASS,  AOD_BIOMASS,  L_AOD_BIOGENIC, AOD_BIOGENIC   &
     &   , L_AOD_OCFF,     AOD_OCFF,     L_AOD_DELTA,    AOD_DELTA      &
     &   , AOD_ABSORPTION, AOD_SCATTERING, I_AOD_TYPE                   &
!                       Properties of Clouds
     &   , N_CONDENSED, TYPE_CONDENSED                                  &
     &   , I_CLOUD, I_CLOUD_REPRESENTATION, W_CLOUD, FRAC_CLOUD         &
     &   , TOT_CLOUD_COVER                                              &
     &   , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                      &
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST                      &
     &   , DP_CORR_STRAT, DP_CORR_CONV                                  &
     &   , LATITUDE                                                     &
!                       Fluxes Calculated
     &   , FLUX_DIRECT, FLUX_DIFFUSE, FLUX_DOWN, FLUX_UP                &
     &   , UV_FLUX_DIRECT, UV_FLUX_DOWN, UV_FLUX_UP                     &
     &   , L_FLUX_DIFFUSE, L_UVFLUX_DIRECT, L_UVFLUX_DOWN, L_UVFLUX_UP  &
!                       Options for Clear-sky Fluxes
     &   , L_CLEAR, I_SOLVER_CLEAR                                      &
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_DOWN_CLEAR, FLUX_UP_CLEAR            &
!                       Arrays specific to the UM
!                       Arrays for Coupling
     &   , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION             &
     &   , ALBEDO_SEA_DIFF, ALBEDO_SEA_DIR, FLANDG                      &
     &   , SEA_FLUX                                                     &
!                       Arrays for diagnostics specific to the UM
     &   , L_FLUX_BELOW_690NM_SURF, WEIGHT_690NM, WEIGHT_UV             &
     &   , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF               &
     &   , L_MOSES_II, L_CTILE                                          &
!     &   , L_MOSES_II, l_cable,L_CTILE                                  &
     &   , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF       &
     &   , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX                       &
     &   , L_SURF_DOWN_CLR, SURF_DOWN_CLR                               &
     &   , L_SURF_UP_CLR, SURF_UP_CLR                                   &
     &   , L_CLOUD_EXTINCTION, CLOUD_EXTINCTION                         &
     &   , CLOUD_WEIGHT_EXTINCTION                                      &
     &   , L_LS_CLOUD_EXTINCTION, LS_CLOUD_EXTINCTION                   &
     &   , LS_CLOUD_WEIGHT_EXTINCTION                                   &
     &   , L_CNV_CLOUD_EXTINCTION, CNV_CLOUD_EXTINCTION                 &
     &   , CNV_CLOUD_WEIGHT_EXTINCTION                                  &
     &   , L_CLOUD_ABSORPTIVITY, CLOUD_ABSORPTIVITY                     &
     &   , CLOUD_WEIGHT_ABSORPTIVITY                                    &
     &   , L_LS_CLOUD_ABSORPTIVITY, LS_CLOUD_ABSORPTIVITY               &
     &   , LS_CLOUD_WEIGHT_ABSORPTIVITY                                 &
     &   , L_CNV_CLOUD_ABSORPTIVITY, CNV_CLOUD_ABSORPTIVITY             &
     &   , CNV_CLOUD_WEIGHT_ABSORPTIVITY                                &
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN, NPD_FIELD                &
     &   , NPD_BAND                                                     &
     &   , NPD_SPECIES                                                  &
     &   , NPD_ESFT_TERM, NPD_SCALE_FNC, NPD_SCALE_VARIABLE             &
     &   , NPD_CONTINUUM                                                &
     &   , NPD_AEROSOL_SPECIES, NPD_HUMIDITIES                          &
     &   , NPD_CLOUD_PARAMETER                                          &
     &   , NPD_THERMAL_COEFF                                            &
     &   , NPD_SURFACE, NPD_ALBEDO_PARM, NPD_AOD_WAVEL                  &
     &   ,l_cable                                                       &
     &   )
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARRAY SIZES
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             FIELD SIZE IN CALLING PROGRAM
     &   , NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_BAND                                                     &
!             NUMBER OF BANDS
     &   , NPD_SPECIES                                                  &
!             NUMBER OF SPECIES
     &   , NPD_CONTINUUM                                                &
!             NUMBER OF CONTINUA
     &   , NPD_AEROSOL_SPECIES                                          &
!             NUMBER OF AEROSOL SPECIES
     &   , NPD_HUMIDITIES                                               &
!             MAXIMUM NUMBER OF HUMIDITIES
     &   , NPD_ESFT_TERM                                                &
!             MAXIMUM NUMBER OF ESFT TERMS
     &   , NPD_SCALE_FNC                                                &
!             NUMBER OF SCALING FUNCTIONS
     &   , NPD_SCALE_VARIABLE                                           &
!             NUMBER OF SCALING VARIABLES
     &   , NPD_CLOUD_PARAMETER                                          &
!             MAXIMUM NUMBER OF CLOUD PARAMETERS
     &   , NPD_THERMAL_COEFF                                            &
!             MAXIMUM NUMBER OF THERMAL COEFFICIENTS
     &   , NPD_SURFACE                                                  &
!             NUMBER OF SURFACE TYPES
     &   , NPD_ALBEDO_PARM                                              &
!             NUMBER OF PARAMETERS FOR DIRECT ALB.
     &   , NPD_COLUMN                                                   &
!             NUMBER OF COLUMNS PER POINT
     &   , NPD_AOD_WAVEL
!             NUMBER OF WAVELENGTHS FOR THE AEROSOL OPTICAL DEPTH
!
!     INCLUDE COMDECKS.
#include "dimfix3a.h"
#include "gasovl3a.h"
#include "clschm3a.h"
#include "clrepp3a.h"
#include "cldcmp3a.h"
#include "cldtyp3a.h"
#include "cldreg3a.h"
#include "angint3a.h"
#include "solver3a.h"
#include "twostr3a.h"
#include "spcrg3a.h"
#include "esftsc3a.h"
#include "aerprm3a.h"
#include "error3a.h"
#include "c_0_dg_c.h"
!
!
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!STR  GENERAL LOGICAL SWITCHES:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_LAYER                                                      &
!             VALUES GIVEN IN LAYERS
     &   , L_CLOUD_LAYER                                                &
!             CLOUD VALUES GIVEN IN LAYERS
     &   , L_CLEAR                                                      &
!             CALCULATE CLEAR-SKY FLUXES
     &   , L_IR_SOURCE_QUAD                                             &
!             USE A QUADRATIC SOURCE FUNCTION
     &   , L_RESCALE                                                    &
!             FLAG FOR DELTA-RESCALING
     &   , L_WENYI                                                      &
!             FLAG FOR WENYI P & T SCALING
     &   , L_2_STREAM_CORRECT                                           &
!             CORRECTION TO 2-STREAM SCHEME
     &   , L_pc2
!             Use PC2 cloud scheme
!
!STR  PARAMETERS CONTROLLING ALGORITHMS:
!     REPRESENTATION OF CLOUDS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_CLOUD
!             CLOUD SCHEME USED
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_GLOBAL_CLOUD_TOP                                           &
!             FLAG TO USE A GLOBAL VALUE FOR THE TOPS OF CLOUDS
     &   , L_INHOM_CLOUD
!             FLAG TO USE SCALING FACTORS FOR INHOMOGENEOUS CLOUD
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CLOUD_TOP_GLOBAL
!             GLOBAL TOPMOST CLOUDY LAYER
      REAL                                                              &
                !, INTENT(IN)
     &     INHOM_CLOUD(NPD_CLOUD_COMPONENT)
!             SCALING FACTORS FOR INHOMOGENEOUS CLOUD
!
!     NUMERICAL ALGORITHMS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     ISOLIR                                                       &
!             VISIBLE OR IR
     &   , I_SOLVER                                                     &
!             SOLVER USED
     &   , I_SOLVER_CLEAR                                               &
!             CLEAR SOLVER USED
     &   , I_2STREAM                                                    &
!             TWO-STREAM SCHEME
     &   , I_ANGULAR_INTEGRATION                                        &
!             ANGULAR INTEGRATION SCHEME
     &   , N_ORDER_GAUSS
!             ORDER OF GAUSSIAN INTEGRATION
!     RANGE OF SPECTRAL BANDS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_FIRST_BAND                                                 &
!             FIRST BAND
     &   , I_LAST_BAND
!             LAST BAND
!
!     GENERAL PROPERTIES OF SPECTRUM:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_BAND                                                       &
!             NUMBER OF SPECTRAL BANDS
     &   , N_ABSORB                                                     &
!             NUMBER OF ABSORBERS
     &   , N_AEROSOL
!             NUMBER OF AEROSOL SPECIES
!
!STR  SOLAR FIELDS:
      REAL                                                              &
                !, INTENT(IN)
     &     SOLAR_TOA(NPD_PROFILE)                                       &
!             INCIDENT SOLAR RADIATION
     &   , SOLAR_FLUX_BAND(NPD_BAND)                                    &
!             NORMALIZED FLUX IN BAND
     &   , SEC_0(NPD_PROFILE)
!             SECANT OF ZENITH ANGLE
!
!STR  ATMOSPHERIC PROFILES:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      REAL                                                              &
                !, INTENT(IN)
     &     P(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             PRESSURE
     &   , T(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             TEMPERATURE
     &   , T_GROUND(NPD_PROFILE)                                        &
!             TEMPERATURE OF GROUND
     &   , T_SOLID(NPD_PROFILE)                                         &
!             TEMPERATURE OF SOLID SURFACE
     &   , T_SEA(NPD_PROFILE)                                           &
!             SURFACE TEMPERATURE OVER OPEN SEA
     &   , T_LEVEL(NPD_PROFILE, 0: NPD_LAYER)                           &
!             TEMPERATURE ON LEVELS
     &   , D_MASS(NPD_PROFILE, NPD_LAYER)                               &
!             MASS THICKNESS OF EACH LAYER
     &   , GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES)
!             GASEOUS MASS MIXING RATIOS
!
!STR  SURFACE PROPERTIES:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_SURFACE(NPD_SURFACE)
!             TYPES OF SURFACES FOR WHICH DATA ARE PRESENT
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SURFACE(NPD_PROFILE)                                       &
!             TYPE OF SURFACE AT THE FOOT OF EACH PROFILE
     &   , I_SPEC_SURFACE(NPD_SURFACE)                                  &
!             METHOD OF SPECIFYING ALBEDO
     &   , N_DIR_ALBEDO_FIT(NPD_SURFACE)
!             NUMBER OF PARAMETERS IN FIT TO DIRECT ALBEDO
      REAL                                                              &
                !, INTENT(IN)
     &     SURFACE_ALBEDO(NPD_BAND, NPD_SURFACE)                        &
!             SURFACE ALBEDO
     &   , ALBEDO_FIELD_DIFF(NPD_PROFILE, NPD_BAND)                     &
!             SPECIFIED DIFFUSE ALBEDOS
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE, NPD_BAND)                      &
!             SPECIFIED DIRECT ALBEDOS
     &   , DIRECT_ALBEDO_PARM(0: NPD_ALBEDO_PARM, NPD_BAND, NPD_SURFACE)&
!             COEFFICIENTS FOR DIRECT ALBEDOS
     &   , EMISSIVITY_GROUND(NPD_BAND, NPD_SURFACE)                     &
!             SURFACE EMISSIVITIES
     &   , EMISSIVITY_FIELD(NPD_PROFILE, NPD_BAND)
!             SPECIFIED EMISSIVITIES
!
!STR  RAYLEIGH SCATTERING:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_RAYLEIGH
!             INCLUDE RAYLEIGH SCATTERING IN THE CALCULATION.
      REAL                                                              &
                !, INTENT(IN)
     &     RAYLEIGH_COEFFICIENT(NPD_BAND)
!             RAYLEIGH COEFFICIENTS
!
!STR  FIELDS FOR GASEOUS ABSORPTION:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_GAS
!             INCLUDE GAS ABSORPTION IN THE CALCULATION
!     GASEOUS OVERLAPS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_GAS_OVERLAP(NPD_BAND)                                      &
!             GAS OVERLAP ASSUMPTION
     &   , I_GAS
!             GAS TO BE CONSIDERED (ONE GAS ONLY)
!     ESFTS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_BAND_ABSORB(NPD_BAND)                                      &
!             NUMBER OF ABSORBERS IN BAND
     &   , INDEX_ABSORB(NPD_SPECIES, NPD_BAND)                          &
!             LIST OF ABSORBERS IN BANDS
     &   , I_BAND_ESFT(NPD_BAND, NPD_SPECIES)                           &
!             NUMBER OF TERMS IN BAND
     &   , I_SCALE_ESFT(NPD_BAND, NPD_SPECIES)                          &
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC(NPD_BAND, NPD_SPECIES)
!             TYPE OF SCALING FUNCTION
      REAL                                                              &
                !, INTENT(IN)
     &     W_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)                 &
!             WEIGHTS FOR ESFT
     &   , K_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)                 &
!             EXPONENTIAL ESFT TERMS
     &   , SCALE_VECTOR(NPD_SCALE_VARIABLE, NPD_ESFT_TERM, NPD_BAND     &
     &        , NPD_SPECIES)                                            &
!             ABSORBER SCALING PARAMETERS
     &   , P_REFERENCE(NPD_SPECIES, NPD_BAND)                           &
!             REFERENCE SCALING PRESSURE
     &   , T_REFERENCE(NPD_SPECIES, NPD_BAND)
!             REFERENCE SCALING TEMPERATURE
!
!     Use modulus of fluxes to remove negative effective extinctions
      LOGICAL, INTENT(IN) :: L_MOD_K_FLUX

!STR  SPECTRAL DATA FOR THE CONTINUUM:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CONTINUUM
!             INCLUDE CONTINUUM ABSORPTION IN THE CALCULATION
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_BAND_CONTINUUM(NPD_BAND)                                   &
!             NUMBER OF CONTINUA IN BANDS
     &   , INDEX_CONTINUUM(NPD_BAND, NPD_CONTINUUM)                     &
!             INDICES OF CONTINUA
     &   , INDEX_WATER                                                  &
!             INDEX OF WATER
     &   , I_SCALE_FNC_CONT(NPD_BAND, NPD_CONTINUUM)
!             TYPE OF SCALING FUNCTION FOR CONTINUUM
      REAL                                                              &
                !, INTENT(IN)
     &     K_CONTINUUM(NPD_BAND, NPD_CONTINUUM)                         &
!             CONTINUUM EXTINCTION COEFFICIENTS
     &   , SCALE_CONTINUUM(NPD_SCALE_VARIABLE, NPD_BAND, NPD_CONTINUUM) &
!             CONTINUUM SCALING PARAMETERS
     &   , P_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)                     &
!             CONTINUUM REFERENCE PRESSURE
     &   , T_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)
!             CONTINUUM REFERENCE TEMPERATURE
!
!
!STR  GENERAL CLOUD FIELDS:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLOUD
!             CLOUDS ARE REQUIRED IN THE CALCULATION
      REAL                                                              &
                !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             AMOUNT OF CLOUD
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)           &
!             FRACTIONS OF DIFFERENT TYPES OF CLOUD
     &   , DP_CORR_STRAT                                                &
!           Decorrelation pressure scale for large scale cloud
     &   , DP_CORR_CONV                                                 &
!           Decorrelation pressure scale for convective cloud
     &   , LATITUDE(NPD_FIELD)
!           Latitude field
      REAL                                                              &
                !, INTENT(INOUT)
     &     TOT_CLOUD_COVER(NPD_PROFILE)
!             Total cloud cover
!
!STR  FIELDS FOR MICROPHYSICAL QUANTITIES:
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_DROP                                                       &
!             INCLUDE DROPLETS IN THE CALCULATION
     &   , L_ICE
!             INCLUDE ICE IN THE CALCULATION
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CONDENSED                                                  &
!             NUMBER OF CONDENSED COMPONENTS IN CLOUDS
     &   , TYPE_CONDENSED(NPD_CLOUD_COMPONENT)                          &
!             TYPES OF CONDENSED COMPONENTS
     &   , I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)                       &
!             PARAMETRIZATION SCHEMES FOR COMPONENTS
     &   , I_CLOUD_REPRESENTATION
!             REPRESENTATION OF MIXING RULE CHOSEN
!
      REAL                                                              &
                !, INTENT(IN)
     &     CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER                &
     &        , NPD_CLOUD_COMPONENT)                                    &
!             MIXING RATIOS OF CONDENSED COMPONENTS
     &   , CONDENSED_DIM_CHAR(NPD_PROFILE, 0: NPD_LAYER                 &
     &        , NPD_CLOUD_COMPONENT)                                    &
!             EFFECTIVE RADII OF CONDENSED COMPONENTS
     &   , CONDENSED_PARAM_LIST(NPD_CLOUD_PARAMETER                     &
     &        , NPD_CLOUD_COMPONENT, NPD_BAND)
!             COEFFICIENTS IN PARAMETRIZATIONS OF CONDENSED PHASES
!
!
!
!STR  FIELDS FOR AEROSOLS:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_AEROSOL
!             INCLUDE AEROSOLS IN THE CALCULATION
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL_SPECIES)
!             PARAMETRIZATION FLAGS FOR AEROSOL
      INTEGER                                                           &
                !, INTENT(IN)
     &     NHUMIDITY(NPD_AEROSOL_SPECIES)
!             NUMBER OF HUMIDITIES
      REAL                                                              &
                !, INTENT(IN)
     &     AEROSOL_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER                  &
     &        , NPD_AEROSOL_SPECIES)
!             NUMBER DENSITY OF AEROSOLS
      REAL                                                              &
                !, INTENT(IN)
     &     AEROSOL_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES       &
     &        , NPD_BAND)                                               &
!             ABSORPTION BY AEROSOLS
     &   , AEROSOL_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES       &
     &        , NPD_BAND)                                               &
!             SCATTERING BY AEROSOLS
     &   , AEROSOL_ASYMMETRY(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES        &
     &        , NPD_BAND)                                               &
!             ASYMMETRY BY AEROSOLS
     &   , HUMIDITIES(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)
!             HUMIDITIES FOR SPECIES
      LOGICAL                                                           &
     &     L_USE_CLEARRH
!             TYPE OF REL HUMIDITY USED FOR HYGROSCOPIC AEROSOLS
!             (.TRUE. -> CLEAR-SKY MEAN, .FALSE. -> GRID-BOX MEAN)
!
!
!STR  FIELDS FOR AEROSOL OPTICAL DEPTH
#include "aodtype.h"

      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_AOD_SULPHATE,                                              &
     &     L_AOD_DUST,                                                  &
     &     L_AOD_SEASALT,                                               &
     &     L_AOD_SOOT,                                                  &
     &     L_AOD_BIOMASS,                                               &
     &     L_AOD_BIOGENIC,                                              &
     &     L_AOD_OCFF,                                                  &
     &     L_AOD_DELTA
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_AOD_WAVEL
      REAL                                                              &
                !, INTENT(IN)
     &     AOD_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES           &
     &                  , NPD_AOD_WAVEL)                                &
     &,    AOD_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES           &
     &                  , NPD_AOD_WAVEL)
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_AOD_TYPE(NPD_AEROSOL_SPECIES)
      INTEGER                                                           &
                !, INTENT(IN)
     &     TYPE_AEROSOL(NPD_AEROSOL_SPECIES)
      REAL                                                              &
                !, INTENT(OUT)
     &     AOD_SULPHATE(NPD_PROFILE, NPD_AOD_WAVEL),                    &
     &     AOD_DUST(NPD_PROFILE, NPD_AOD_WAVEL),                        &
     &     AOD_SEASALT(NPD_PROFILE, NPD_AOD_WAVEL),                     &
     &     AOD_SOOT(NPD_PROFILE, NPD_AOD_WAVEL),                        &
     &     AOD_BIOMASS(NPD_PROFILE, NPD_AOD_WAVEL),                     &
     &     AOD_BIOGENIC(NPD_PROFILE, NPD_AOD_WAVEL),                    &
     &     AOD_OCFF(NPD_PROFILE, NPD_AOD_WAVEL),                        &
     &     AOD_DELTA(NPD_PROFILE, NPD_AOD_WAVEL)
!
!
!STR  FITTING OF THE PLANCKIAN FUNCTION:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_DEG_FIT
!             DEGREE OF THERMAL FITTING FNC.
      REAL                                                              &
                !, INTENT(IN)
     &     THERMAL_COEFFICIENT(0: NPD_THERMAL_COEFF-1, NPD_BAND)        &
!             COEFFICIENTS OF SOURCE TERMS
     &   , T_REF_PLANCK
!             PLANCKIAN REFERENCE TEMPERATURE
!
!STR  DOPPLER BROADENING
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_DOPPLER(NPD_SPECIES)
!             FLAGS TO ACTIVATE DOPPLER CORRECTIONS
      REAL                                                              &
                !, INTENT(IN)
     &     DOPPLER_CORRECTION(NPD_SPECIES)
!             DOPPLER BROADENING TERM
      REAL                                                              &
                !, INTENT(OUT)
     &     WEIGHT_BAND(NPD_BAND)
!             WEIGHTING FUNCTION FOR BANDS
!
!STR  CONTROL OF SCATTERING:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SCATTER_METHOD_BAND(NPD_BAND)
!             METHOD OF TREATING SCATTERING IN EACH BAND
!
!
!STR  FLUXES CALCULATED:
      REAL                                                              &
                !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DIRECT FLUX
     &   , FLUX_DIFFUSE(NPD_PROFILE, 0: NPD_LAYER)                      &
!             DIFFUSE_FLUX     
     &   , FLUX_DOWN(NPD_PROFILE, 0: NPD_LAYER)                         &
!             DOWNWARD FLUX
     &   , FLUX_UP(NPD_PROFILE, 0: NPD_LAYER)                           &
!             UPWARD FLUX
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                 &
!             CLEAR DIRECT FLUX
     &   , FLUX_DOWN_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                   &
!             CLEAR DOWNWARD FLUX
     &   , FLUX_UP_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR UPWARD FLUX
!

!  UV FLUXES
      REAL                                                              &
                 !, INTENT(OUT)
     &     UV_FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)                    &
!             DIRECT UV-FLUX
     &   , UV_FLUX_UP(NPD_PROFILE, 0: NPD_LAYER)                        &
!             UPWARD UV-FLUX
     &   , UV_FLUX_DOWN(NPD_PROFILE, 0: NPD_LAYER)
!             DOWNWARD UV_FLUX

!STR  ARRAYS SPECIFIC TO THE UNIFIED MODEL
!
!     SWITCHES FOR DIAGNOSTICS:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_FLUX_BELOW_690NM_SURF                                      &
!             FLUX BELOW 690NM AT SURFACE TO BE CALCULATED
     &   , L_MOSES_II                                                   &
!             SURFACE SW REQUIRED FOR MOSES II
     &   , l_cable                                                   &
     &   , L_CTILE                                                      &
!             COASTAL TILING SWITCH
     &   , L_SURFACE_DOWN_FLUX                                          &
!             DOWNWARD SURFACE FLUX REQUIRED
     &   , L_SURF_DOWN_CLR                                              &
!             CALCULATE DOWNWARD CLEAR FLUX
     &   , L_SURF_UP_CLR                                                &
!             CALCULATE UPWARD CLEAR FLUX
     &   , L_CLOUD_ABSORPTIVITY                                         &
!             FLAG TO CALCULATE ABSORPTIVITY OF CLOUDS
!             (ONLY INFRA-RED)
     &   , L_CLOUD_EXTINCTION                                           &
!             FLAG TO CALCULATE EXTINCTION OF CLOUDS
!             (ONLY SOLAR)
     &   , L_LS_CLOUD_ABSORPTIVITY                                      &
!             FLAG TO CALCULATE ABSORPTIVITY OF LAYER CLOUDS
!             (ONLY INFRA-RED)
     &   , L_LS_CLOUD_EXTINCTION                                        &
!             FLAG TO CALCULATE EXTINCTION OF LAYER CLOUDS
!             (ONLY SOLAR)
     &   , L_CNV_CLOUD_ABSORPTIVITY                                     &
!             FLAG TO CALCULATE ABSORPTIVITY OF CONV.CLOUDS
!             (ONLY INFRA-RED)
     &   , L_CNV_CLOUD_EXTINCTION                                       &
!             FLAG TO CALCULATE EXTINCTION OF CONV.CLOUDS
!             (ONLY SOLAR)
     &   , L_FLUX_DIFFUSE                                               &
!             FLAG FOR THE DIFFUSE DOWNWORD FLUX     
     &   , L_UVFLUX_DIRECT                                              &
!             FLAG FOR DIRECT ULTRAVIOLET FLUXES
     &   , L_UVFLUX_UP                                                  &
!             FLAF FOR UPWARD ULTRAVIOLET FLUXES
     &   , L_UVFLUX_DOWN
!             FLAG FOR DOWNWARD ULTRAVIOLET FLUXES

!
!     ARRAYS FOR USE WITH COUPLING:

      REAL                                                              &
               !, INTENT(IN)
     &     FLANDG(NPD_PROFILE)
!            Land fraction
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_FRAC_SOL_POINT                                             &
!             NUMBER OF POINTS WITH FRACTIONAL ICE/LAND COVER
     &   , I_FRAC_SOL_POINT(NPD_PROFILE)
!             INDICES OF POINTS WITH FRACTIONAL ICE/LAND COVER
      REAL                                                              &
                !, INTENT(IN)
     &     ICE_FRACTION(NPD_PROFILE)
!             ICE FRACTION
      REAL                                                              &
                !, INTENT(IN)
     &     ALBEDO_SEA_DIFF(NPD_PROFILE, NPD_BAND)                       &
!             DIFFUSE ALBEDO FOR OPEN SEA
     &   , ALBEDO_SEA_DIR(NPD_PROFILE, NPD_BAND)
!             DIRECT ALBEDO FOR OPEN SEA
!
!     ARRAYS FOR USE WITH DIAGNOSTICS:
      REAL                                                              &
                !, INTENT(IN)
     &     WEIGHT_690NM(NPD_BAND)                                       &
!             WEIGHTS FOR EACH BAND FOR REGION BELOW 690 NM

     &   , WEIGHT_UV(NPD_BAND)
!             WEIGHTS FOR EACH BAND FOR THE UV INTERVAL
!
!     SURFACE FLUXES FOR COUPLING OR DIAGNOSTIC USE
      REAL                                                              &
                !, INTENT(OUT)
     &     SEA_FLUX(NPD_PROFILE)                                        &
!             NET DOWNWARD FLUX INTO SEA
     &   , SURFACE_DOWN_FLUX(NPD_PROFILE)                               &
!             DOWNWARD FLUX AT SURFACE
     &   , SURF_DOWN_CLR(NPD_PROFILE)                                   &
!             CLEAR-SKY DOWNWARD FLUX AT SURFACE
     &   , SURF_UP_CLR(NPD_PROFILE)                                     &
!             CLEAR-SKY UPWARD FLUX AT SURFACE
     &   , FLUX_BELOW_690NM_SURF(NPD_PROFILE)                           &
!          GRID BOX MEAN NET SURFACE FLUX BELOW 690NM
     &   , FL_SEA_BELOW_690NM_SURF(NPD_PROFILE)                         &
!          OPEN SEA NET SURFACE FLUX BELOW 690NM
     &   , SURF_VIS_DIR(NPD_PROFILE)                                    &
!             DOWNWARD SURFACE DIRECT BEAM VISIBLE FLUX
     &   , SURF_VIS_DIF(NPD_PROFILE)                                    &
!             DOWNWARD SURFACE DIFFUSE VISIBLE FLUX
     &   , SURF_NIR_DIR(NPD_PROFILE)                                    &
!             DOWNWARD SURFACE DIRECT BEAM NEAR-INFRARED FLUX
     &   , SURF_NIR_DIF(NPD_PROFILE)
!             DOWNWARD SURFACE DIFFUSE NEAR-INFRARED FLUX
!
!     DIAGNOSTICS FOR CLOUDS
      REAL                                                              &
                !, INTENT(OUT)
     &     CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)                   &
!             ABSORPTIVITY OF CLOUD WEIGHTED BY CLOUD FRACTION
!             AND UPWARD CLEAR-SKY INFRA-RED FLUX.
     &   , CLOUD_WEIGHT_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)            &
!             WEIGHTS TO BE APPLIED TO ABSORPTIVIES.
     &   , LS_CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)                &
!             ABSORPTIVITY OF LAYER CLOUD WEIGHTED BY CLOUD FRACTION
!             AND UPWARD CLEAR-SKY INFRA-RED FLUX.
     &   , LS_CLOUD_WEIGHT_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)         &
!             WEIGHTS TO BE APPLIED TO LAYER CLD. ABSORPTIVIES.
     &   , CNV_CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)               &
!             ABSORPTIVITY OF CONV.CLOUD WEIGHTED BY CLOUD FRACTION
!             AND UPWARD CLEAR-SKY INFRA-RED FLUX.
     &   , CNV_CLOUD_WEIGHT_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)        &
!             WEIGHTS TO BE APPLIED TO CONV.CLD ABSORPTIVIES.
     &   , CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                     &
!             ABSORPTIVITY OF CLOUD WEIGHTED BY CLOUD FRACTION
!             AND DOWNWARD CLEAR-SKY SOLAR FLUX.
     &   , CLOUD_WEIGHT_EXTINCTION(NPD_PROFILE, NPD_LAYER)              &
!             WEIGHTS TO BE APPLIED TO EXTINCTIONS.
     &   , LS_CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                  &
!             ABSORPTIVITY OF LAYER CLOUD WEIGHTED BY CLOUD FRACTION
!             AND DOWNWARD CLEAR-SKY SOLAR FLUX.
     &   , LS_CLOUD_WEIGHT_EXTINCTION(NPD_PROFILE, NPD_LAYER)           &
!             WEIGHTS TO BE APPLIED TO LAYER CLD. EXTINCTIONS.
     &   , CNV_CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                 &
!             ABSORPTIVITY OF CONV.CLOUD WEIGHTED BY CLOUD FRACTION
!             AND DOWNWARD CLEAR-SKY SOLAR FLUX.
     &   , CNV_CLOUD_WEIGHT_EXTINCTION(NPD_PROFILE, NPD_LAYER)
!             WEIGHTS TO BE APPLIED TO CONV. CLD. EXTINCTIONS.
!
!
!
!     LOCAL ARGUMENTS.
!     GENERAL POINTERS:
      INTEGER                                                           &
     &     I_TOP                                                        &
!             TOP LEVEL OF PROFILES
     &   , I_BAND                                                       &
!             SPECTRAL BAND
     &   , IEX                                                          &
!             INDEX OF ESFT TERM
     &   , N_AUGMENT                                                    &
!             LENGTH OF LONG FLUX VECTOR
     &   , N_GAS                                                        &
!             NUMBER OF ACTIVE GASES
     &   , I_GAS_BAND                                                   &
!             SINGLE VARIABLE FOR GAS IN BAND
     &   , N_CONTINUUM                                                  &
!             NUMBER OF CONTINUA IN BAND
     &   , I_CONTINUUM                                                  &
!             CONTINUUM NUMBER
     &   , I_CONTINUUM_POINTER(NPD_CONTINUUM)
!             POINTERS TO CONTINUA
!
!     VARIABLES FOR SURFACE PROPERTIES:
      INTEGER                                                           &
     &     N_POINT_TYPE(NPD_SURFACE)                                    &
!             NUMBER OF POINTS OF EACH TYPE
     &   , INDEX_SURFACE(NPD_PROFILE, NPD_SURFACE)
!             INDICES OF EACH SURFACE TYPE
!
!     POINTERS TO THE CONTENTS OF LAYERS:
      INTEGER                                                           &
     &     N_FREE_PROFILE(NPD_LAYER)                                    &
!             NUMBER OF FREE PROFILES
     &   , I_FREE_PROFILE(NPD_PROFILE, NPD_LAYER)                       &
!             COLUMNS CONTAINING FREE PROFILES
     &   , N_CLOUD_TOP                                                  &
!             TOPMOST CLOUDY LAYER
     &   , N_CLOUD_PROFILE(NPD_LAYER)                                   &
!             NUMBER OF CLOUDY PROFILES
     &   , I_CLOUD_PROFILE(NPD_PROFILE, NPD_LAYER)
!             PROFILES CONTAINING CLOUDS
!
!     POINTERS TO TYPES OF CLOUDS:
      LOGICAL                                                           &
     &     L_CLOUD_CMP(NPD_CLOUD_COMPONENT)
!             LOGICAL SWITCHES TO INCLUDE COMPONENTS
      INTEGER                                                           &
     &     I_PHASE_CMP(NPD_CLOUD_COMPONENT)                             &
!             PHASES OF COMPONENTS
     &   , I_CLOUD_TYPE(NPD_CLOUD_COMPONENT)                            &
!             TYPES OF CLOUD TO WHICH EACH COMPONENT CONTRIBUTES
     &   , TYPE_REGION(NPD_REGION)                                      &
!             The types of the regions
     &   , N_REGION                                                     &
!             Number of regions of aggregated cloud
     &   , I_REGION_CLOUD(NPD_CLOUD_TYPE)
!             REGIONS IN WHICH PARTICULAR TYPE OF CLOUD FALL
!
!     FRACTIONAL COVERAGE OF DIFFERENT REGIONS:
      REAL                                                              &
     &     FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTION OF TOTAL CLOUD OCCUPIED BY SPECIFIC REGIONS
!
!     POINTER TO TABLE OF HUMIDITY:
      INTEGER                                                           &
     &     I_HUMIDITY_POINTER(NPD_PROFILE, NPD_LAYER)
!             POINTER TO LOOK-UP TABLE FOR AEROSOLS
!
!     CONTROLLING VARIABLES:
      INTEGER                                                           &
     &     I                                                            &
!             LOOP INDEX
     &   , J                                                            &
!             LOOP INDEX
     &   , K                                                            &
!             LOOP INDEX
     &   , L
!             LOOP INDEX
!
!     LOGICAL SWITCHES:
      LOGICAL                                                           &
     &     L_GAS_BAND                                                   &
!             FLAG TO INCLUDE GASEOUS ABSORPTION IN A PARTICULAR BAND
     &   , L_MOIST_AEROSOL                                              &
!             FLAG FOR MOIST AEROSOL
     &   , L_AEROSOL_DENSITY
!             FLAG FOR CALCULATION OF ATMOSPHERIC DENSITY FOR AEROSOLS
!
!     SURFACE PROPERTIES:
      REAL                                                              &
     &     ALBEDO_SURFACE_DIFF(NPD_PROFILE)                             &
!             DIFFUSE SURFACE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT SURFACE ALBEDO
!
      REAL                                                              &
     &     INC_SOLAR_FLUX_BAND(NPD_PROFILE)
!             INCIDENT SOLAR FLUX IN BAND
      REAL                                                              &
     &     GAS_FRAC_RESCALED(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES)    &
!             RESCALED GAS MIXING RATIOS
     &   , AMOUNT_CONTINUUM(NPD_PROFILE, 0: NPD_LAYER, NPD_CONTINUUM)   &
!             AMOUNTS OF CONTINUA
     &   , K_CONTINUUM_MONO(NPD_CONTINUUM)
!             MONOCHROMATIC CONTINUUM COMPONENTS
!
!     THERMAL ARRAYS:
      REAL                                                              &
     &     PLANCK_SOURCE_BAND(NPD_PROFILE, 0: NPD_LAYER)                &
!             PLANCK FUNCTION IN BAND AT LEVELS
     &   , DIFF_PLANCK_BAND(NPD_PROFILE, NPD_LAYER)                     &
!             DIFFERENTIAL THERMAL SOURCE IN BAND
     &   , DIFF_PLANCK_BAND_2(NPD_PROFILE, NPD_LAYER)                   &
!             2 x 2ND DIFF. THERMAL SOURCE IN BAND
     &   , THERMAL_GROUND_BAND(NPD_PROFILE)
!             GROUND SOURCE FUNCTION IN BAND
!
!     ATMOSPHERIC DENSITIES:
      REAL                                                              &
     &     DENSITY(NPD_PROFILE, 0: NPD_LAYER)                           &
!             OVERALL DENSITY
     &   , MOLAR_DENSITY_WATER(NPD_PROFILE, 0: NPD_LAYER)               &
!             MOLAR DENSITY OF WATER
     &   , MOLAR_DENSITY_FRN(NPD_PROFILE, 0: NPD_LAYER)
!             MOLAR DENSITY OF FOREIGN SPECIES
!
!     FIELDS FOR MOIST AEROSOLS:
      REAL                                                              &
     &     DELTA_HUMIDITY                                               &
!             INCREMENT IN LOOK-UP TABLE FOR HUM.
     &   , MEAN_REL_HUMIDITY(NPD_PROFILE, NPD_LAYER)
!             MEAN RELATIVE HUMIDITY OF LAYERS
!
!     FUNDAMENTAL OPTICAL PROPERTIES OF LAYERS:
      REAL                                                              &
     &     K_GREY_TOT_FREE(NPD_PROFILE, NPD_LAYER)                      &
!             FREE TOTAL GREY EXTINCTION
     &   , K_EXT_SCAT_FREE(NPD_PROFILE, NPD_LAYER)                      &
!             FREE SCATTERING EXTINCTION
     &   , ASYMMETRY_FREE(NPD_PROFILE, NPD_LAYER)                       &
!             FREE ASYMMETRIES
     &   , FORWARD_SCATTER_FREE(NPD_PROFILE, NPD_LAYER)                 &
!             FREE FORWARD SCATTERING
     &   , K_GREY_TOT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)     &
!             TOTAL CLOUDY GREY EXTINCTION
     &   , K_EXT_SCAT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)     &
!             CLOUDY SCATTERING EXTINCTION
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)      &
!             CLOUDY ASYMMETRIES
     &   , FORWARD_SCATTER_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY FORWARD SCATTERING
!
!     LOCAL RADIATIVE FLUXES:
      REAL                                                              &
     &     FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)                       &
!             TOTAL FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)                 &
!             CLEAR TOTAL FLUX
     &   , FLUX_DIRECT_BAND(NPD_PROFILE, 0: NPD_LAYER)                  &
!             DIRECT FLUX IN BAND
     &   , FLUX_TOTAL_BAND(NPD_PROFILE, 2*NPD_LAYER+2)                  &
!             TOTAL FLUX IN BAND
     &   , FLUX_DIRECT_CLEAR_BAND(NPD_PROFILE, 0: NPD_LAYER)            &
!             DIRECT FLUX IN BAND
     &   , FLUX_TOTAL_CLEAR_BAND(NPD_PROFILE, 2*NPD_LAYER+2)            &
!             TOTAL FLUX IN BAND
     &   , PLANCK_FLUX(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCKIAN FLUX IN BAND
!
!     COEFFICIENTS FOR THE TRANSFER OF ENERGY BETWEEN
!     PARTIALLY CLOUDY LAYERS:
      REAL                                                              &
     &     CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)  &
!             COEFFICIENTS DEFINING OVERLAPPING OPTIONS FOR CLOUDS:
!             THESE ALSO DEPEND ON THE SOLVER SELECTED.
     &   , CLOUD_OVERLAP_TMP(NPD_PROFILE, 0: NPD_LAYER)                 &
     &   , W_FREE(NPD_PROFILE, NPD_LAYER)
!             CLEAR-SKY FRACTION
      INTEGER                                                           &
     &     N_COLUMN(NPD_PROFILE)
!             NUMBER OF COLUMNS REQUIRED
      LOGICAL                                                           &
     &     L_COLUMN(NPD_PROFILE, NPD_LAYER, NPD_COLUMN)
!             COLUMN FLAGS FOR COLUMNS
      REAL                                                              &
     &     AREA_COLUMN(NPD_PROFILE, NPD_COLUMN)
!             AREAS OF COLUMNS
!
!     LOCAL VARIABLES SPECIFIC TO THE UNIFIED MODEL.
      REAL                                                              &
     &     PLANCK_FREEZE_SEA                                            &
!             PLANCK FUNCTION OVER FREEZING SEA
     &   , PLANCK_LEADS_SEA(NPD_PROFILE)
!             PLANCK FUNCTION OVER SEA LEADS
!
!     SECONDARY ARRAYS FOR DIAGNOSTICS
      REAL                                                              &
     &     CLOUD_ABSORPTIVITY_BAND(NPD_PROFILE, NPD_LAYER)              &
!             ABSORPTIVITY OF CLOUD IN A PARTICULAR BAND
     &   , CLOUD_EXTINCTION_BAND(NPD_PROFILE, NPD_LAYER)                &
!             ABSORPTIVITY OF CLOUD IN A PARTICULAR BAND
     &   , LS_CLOUD_ABSORPTIVITY_BAND(NPD_PROFILE, NPD_LAYER)           &
!             ABSORPTIVITY OF CLOUD IN A PARTICULAR BAND
     &   , LS_CLOUD_EXTINCTION_BAND(NPD_PROFILE, NPD_LAYER)             &
!             ABSORPTIVITY OF CLOUD IN A PARTICULAR BAND
     &   , CNV_CLOUD_ABSORPTIVITY_BAND(NPD_PROFILE, NPD_LAYER)          &
!             ABSORPTIVITY OF CLOUD IN A PARTICULAR BAND
     &   , CNV_CLOUD_EXTINCTION_BAND(NPD_PROFILE, NPD_LAYER)
!             ABSORPTIVITY OF CLOUD IN A PARTICULAR BAND
!
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     SET_CLOUD_POINTER, SET_CLOUD_GEOMETRY                        &
     &   , OVERLAP_MIX_MAXIMUM, OVERLAP_MIX_RANDOM                      &
     &   , SPLIT_MAXIMUM, COLLECT_SURFACE, INITIALIZE_FLUX              &
     &   , CALCULATE_DENSITY, SET_MOIST_AEROSOL_PROPERTIES              &
     &   , SCALE_ABSORB, RESCALE_CONTINUUM, GREY_EXTINCTION             &
     &   , RESCALE_ASYMMETRY                                            &
     &   , DIFF_PLANCK_SOURCE, SET_SURFACE_PROPERTIES                   &
     &   , SOLVE_BAND_WITHOUT_GAS, SOLVE_BAND_ONE_GAS                   &
     &   , SOLVE_BAND_RANDOM_OVERLAP, SOLVE_BAND_FESFT                  &
     &   , SOLVE_BAND_CLR_FESFT, SOLVE_BAND_K_EQV                       &
     &   , AUGMENT_TOTAL_FLUX, ASSIGN_FLUX                              &
     &   , R2_INIT_COUPLE_DIAG, R2_COUPLE_DIAG                          &
     &   , OVERLAP_TRIPLE                                               &
     &   , COMPUTE_AOD
!     FUNCTIONS CALLED:
      LOGICAL                                                           &
     &     L_CLOUD_DENSITY
!             FLAG FOR CALCULATION OF ATMOSPHERIC DENSITIES FOR CLOUDS
      EXTERNAL                                                          &
     &     L_CLOUD_DENSITY
!
!
!     SETTING OF PROPERTIES OF ARRAYS:
#include "clrepd3a.h"
!
!
!
!
!     INITIAL DETERMINATION OF FLAGS AND SWITCHES:
!
!
!     Set the number of fluxes to be incremented at each grid-point.
!      print *,'flux_calc.F90'

      N_AUGMENT=2*(N_LAYER+1)
!
!     SET THE TOP LEVEL OF THE PROFILES.
      IF (L_LAYER) THEN
         I_TOP=1
      ELSE
         I_TOP=0
      ENDIF
!
!
!
!     INITIAL CALCULATIONS FOR SURFACE PROPERTIES:
!
!     COLLECT POINTS WITH THE SAME SURFACE SPECIFICATION.
! DEPENDS ON: collect_surface
      CALL COLLECT_SURFACE(N_PROFILE                                    &
     &   , I_SURFACE                                                    &
     &   , N_POINT_TYPE, INDEX_SURFACE                                  &
     &   , NPD_PROFILE, NPD_SURFACE                                     &
     &   )
!
!
!
!     INITIAL CALCULATIONS FOR AEROSOLS:
!
!     SET THE SPECTRALLY INDEPENDENT PROPERTIES OF MOIST AEROSOLS.
      L_MOIST_AEROSOL=.FALSE.
      DO J=1, N_AEROSOL
         L_MOIST_AEROSOL=L_MOIST_AEROSOL.OR.                            &
     &      (I_AEROSOL_PARAMETRIZATION(J)                               &
     &       == IP_AEROSOL_PARAM_MOIST)
      ENDDO
!
      IF (L_MOIST_AEROSOL) THEN
! DEPENDS ON: set_moist_aerosol_properties
         CALL SET_MOIST_AEROSOL_PROPERTIES(IERR                         &
     &      , N_PROFILE, N_LAYER, L_LAYER                               &
     &      , L_USE_CLEARRH                                             &
     &      , N_AEROSOL, I_AEROSOL_PARAMETRIZATION, NHUMIDITY           &
     &      , GAS_MIX_RATIO(1, 0, INDEX_WATER), T, P, W_CLOUD           &
     &      , DELTA_HUMIDITY, MEAN_REL_HUMIDITY, I_HUMIDITY_POINTER     &
     &      , NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES               &
     &      )
         IF (IERR /= I_NORMAL) RETURN

      ENDIF
!
!
!     CHECK WHETHER THE DENSITIES WILL BE NEEDED FOR
!     UNPARAMETRIZED AEROSOLS.
      L_AEROSOL_DENSITY=.FALSE.
      IF (L_AEROSOL) THEN
         DO J=1, N_AEROSOL
            L_AEROSOL_DENSITY=L_AEROSOL_DENSITY.OR.                     &
     &         (I_AEROSOL_PARAMETRIZATION(J) ==                         &
     &            IP_AEROSOL_PARAM_MOIST)                               &
     &         .OR.(I_AEROSOL_PARAMETRIZATION(J) ==                     &
     &            IP_AEROSOL_UNPARAMETRIZED)
         ENDDO
      ENDIF
!
!
!
!     INITIAL CALCULATIONS FOR CLOUDS:
!
!     SET POINTERS TO THE TYPES OF CLOUD.
! DEPENDS ON: set_cloud_pointer
      CALL SET_CLOUD_POINTER(IERR                                       &
     &   , N_CONDENSED, TYPE_CONDENSED, I_CLOUD_REPRESENTATION          &
     &   , L_DROP, L_ICE                                                &
     &   , I_PHASE_CMP, I_CLOUD_TYPE, L_CLOUD_CMP                       &
     &   )
      IF (IERR /= I_NORMAL) RETURN
!
!
!     SET THE GEOMETRY OF THE CLOUDS.
! DEPENDS ON: set_cloud_geometry
      CALL SET_CLOUD_GEOMETRY(N_PROFILE, N_LAYER                        &
     &   , L_GLOBAL_CLOUD_TOP, N_CLOUD_TOP_GLOBAL                       &
     &   , W_CLOUD                                                      &
     &   , N_CLOUD_PROFILE, I_CLOUD_PROFILE                             &
     &   , N_CLOUD_TOP                                                  &
     &   , N_FREE_PROFILE, I_FREE_PROFILE                               &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )

!     Scale the condensed water contents to simulate
!     inhomogeneities in the clouds.
      IF ( L_INHOM_CLOUD ) THEN
        DO I=1,NPD_CLOUD_COMPONENT
          CONDENSED_MIX_RATIO( 1:N_PROFILE, N_CLOUD_TOP:N_LAYER, I )    &
     &      = INHOM_CLOUD(I) *                                          &
     &        CONDENSED_MIX_RATIO(1:N_PROFILE,N_CLOUD_TOP:N_LAYER,I)
        ENDDO
      ENDIF
!
      IF ( (I_CLOUD == IP_CLOUD_TRIPLE).OR.                             &
     &     (I_CLOUD == IP_CLOUD_PART_CORR_CNV) ) THEN
!        There are three regions with this option
         N_REGION=3
         TYPE_REGION(1)=IP_REGION_CLEAR
         TYPE_REGION(2)=IP_REGION_STRAT
         TYPE_REGION(3)=IP_REGION_CONV
!        AGGREGATE CLOUDS INTO REGIONS FOR SOLVING.
! DEPENDS ON: aggregate_cloud
         CALL AGGREGATE_CLOUD(IERR                                      &
     &      , N_PROFILE, N_LAYER, N_CLOUD_TOP                           &
     &      , I_CLOUD, I_CLOUD_REPRESENTATION                           &
     &      , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION)                     &
     &      , FRAC_CLOUD                                                &
     &      , I_REGION_CLOUD, FRAC_REGION                               &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
      ELSE IF ( (I_CLOUD == IP_CLOUD_MIX_MAX).OR.                       &
     &          (I_CLOUD == IP_CLOUD_MIX_RANDOM).OR.                    &
     &          (I_CLOUD == IP_CLOUD_PART_CORR) ) THEN
!         There will be only one cloudy region.
          N_REGION=2
          TYPE_REGION(1)=IP_REGION_CLEAR
          TYPE_REGION(2)=IP_REGION_STRAT
          DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
              FRAC_REGION(L, I, 2)=1.0E+00
            ENDDO
          ENDDO
      ENDIF
!
!     CALCULATE ENERGY TRANSFER COEFFICIENTS IN A MIXED COLUMN,
!     OR SPLIT THE ATMOSPHERE INTO COLUMNS WITH A COLUMN MODEL:
!
      IF (I_CLOUD == IP_CLOUD_MIX_MAX) THEN
!
! DEPENDS ON: overlap_mix_maximum
         CALL OVERLAP_MIX_MAXIMUM(N_PROFILE, N_LAYER, N_CLOUD_TOP       &
     &      , ISOLIR, I_SOLVER                                          &
     &      , W_CLOUD, W_FREE                                           &
     &      , CLOUD_OVERLAP                                             &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
      ELSE IF (I_CLOUD == IP_CLOUD_MIX_RANDOM) THEN
! DEPENDS ON: overlap_mix_random
         CALL OVERLAP_MIX_RANDOM(N_PROFILE, N_LAYER, N_CLOUD_TOP        &
     &      , ISOLIR, I_SOLVER                                          &
     &      , W_CLOUD, W_FREE                                           &
     &      , CLOUD_OVERLAP                                             &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
!
      ELSE IF (I_CLOUD == IP_CLOUD_TRIPLE) THEN
!
!        CALCULATE OVERLAPS FOR THE TRIPLE DECOMPOSITION
!        OF THE COLUMN INTO STRATIFORM AND CONVECTIVE PARTS.
! DEPENDS ON: overlap_triple
         CALL OVERLAP_TRIPLE(N_PROFILE, N_LAYER, N_CLOUD_TOP            &
     &      , W_CLOUD, W_FREE, N_REGION, FRAC_REGION                    &
     &      , CLOUD_OVERLAP                                             &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
         IF (IERR /= I_NORMAL) RETURN
!
      ELSE IF ( (I_CLOUD == IP_CLOUD_PART_CORR).OR.                     &
     &          (I_CLOUD == IP_CLOUD_PART_CORR_CNV) ) THEN
!
! DEPENDS ON: overlap_coupled
          CALL OVERLAP_COUPLED(N_PROFILE, N_LAYER, N_CLOUD_TOP          &
     &      , W_CLOUD, W_FREE, N_REGION, TYPE_REGION, FRAC_REGION       &
     &      , P(:,1:), I_CLOUD                                          &
     &      , CLOUD_OVERLAP                                             &
     &      , NPD_PROFILE, NPD_LAYER, NPD_OVERLAP_COEFF, NPD_REGION     &
     &      , NPD_FIELD                                                 &
     &      , 1, DP_CORR_STRAT, DP_CORR_CONV, TOT_CLOUD_COVER           &
     &      , LATITUDE)

          IF (I_CLOUD == IP_CLOUD_PART_CORR) THEN
!           overlap_coupled uses a different order for the elements in
!           cloud_overlap when there is only one cloudy region:
            CLOUD_OVERLAP_TMP=CLOUD_OVERLAP(:,:,2)
            CLOUD_OVERLAP(:,:,2)=CLOUD_OVERLAP(:,:,3)
            CLOUD_OVERLAP(:,:,3)=CLOUD_OVERLAP_TMP
          ENDIF
!
      ELSE IF (I_CLOUD == IP_CLOUD_COLUMN_MAX) THEN
!
! DEPENDS ON: split_maximum
         CALL SPLIT_MAXIMUM(N_PROFILE, N_LAYER                          &
     &      , W_CLOUD                                                   &
     &      , N_COLUMN, AREA_COLUMN, L_COLUMN                           &
     &      , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                        &
     &      )
      ENDIF
!
!
!     CALCULATE THE ATMOSPHERIC DENSITIES:
!
      IF ( L_CONTINUUM                                                  &
     &      .OR.L_AEROSOL_DENSITY                                       &
     &      .OR.(L_CLOUD                                                &
! DEPENDS ON: l_cloud_density
     &      .AND.L_CLOUD_DENSITY(N_CONDENSED, I_PHASE_CMP, L_CLOUD_CMP  &
     &                        , I_CONDENSED_PARAM                       &
     &                        ) ) ) THEN
! DEPENDS ON: calculate_density
         CALL CALCULATE_DENSITY(N_PROFILE, N_LAYER, L_CONTINUUM         &
     &      , GAS_MIX_RATIO(1, 0, INDEX_WATER)                          &
     &      , P, T, I_TOP                                               &
     &      , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN           &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
      ENDIF
!
!
!
!     INITIALIZE THE TOTAL FLUXES.
!
! DEPENDS ON: initialize_flux
      CALL INITIALIZE_FLUX(N_PROFILE, N_LAYER, N_AUGMENT                &
     &   , ISOLIR                                                       &
     &   , FLUX_DIRECT, FLUX_TOTAL                                      &
     &   , L_CLEAR                                                      &
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                          &
     &   , 0.0E+00                                                      &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   , .FALSE.                                                      &
     &   )
!
!     INITIALIZE THE PLANCKIAN FLUXES IF REQUIRED.
      IF (ISOLIR == IP_INFRA_RED) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               PLANCK_FLUX(L, I)=0.0E+00
            ENDDO
         ENDDO
      ENDIF
!
!
!     INITIALIZATION OF DIAGNOSTICS AND COUPLING ARRAYS FOR THE
!     UNIFIED MODEL.
! DEPENDS ON: r2_init_couple_diag
      CALL R2_INIT_COUPLE_DIAG(N_PROFILE                                &
     &   , SEA_FLUX                                                     &
     &   , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX                       &
     &   , L_SURF_DOWN_CLR, SURF_DOWN_CLR                               &
     &   , L_SURF_UP_CLR, SURF_UP_CLR                                   &
     &   , L_FLUX_DIFFUSE, FLUX_DIFFUSE                                 &
     &   , L_UVFLUX_DIRECT, UV_FLUX_DIRECT                              &
     &   , L_UVFLUX_UP, UV_FLUX_UP                                      &
     &   , L_UVFLUX_DOWN, UV_FLUX_DOWN                                  &                                  
     &   , L_FLUX_BELOW_690NM_SURF                                      &
     &   , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF               &
     &   , L_MOSES_II                                                   &
     &   , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF       &
     &   , L_CLOUD_EXTINCTION, CLOUD_EXTINCTION                         &
     &   , CLOUD_WEIGHT_EXTINCTION                                      &
     &   , L_LS_CLOUD_EXTINCTION, LS_CLOUD_EXTINCTION                   &
     &   , LS_CLOUD_WEIGHT_EXTINCTION                                   &
     &   , L_CNV_CLOUD_EXTINCTION, CNV_CLOUD_EXTINCTION                 &
     &   , CNV_CLOUD_WEIGHT_EXTINCTION                                  &
     &   , L_CLOUD_ABSORPTIVITY, CLOUD_ABSORPTIVITY                     &
     &   , CLOUD_WEIGHT_ABSORPTIVITY                                    &
     &   , L_LS_CLOUD_ABSORPTIVITY, LS_CLOUD_ABSORPTIVITY               &
     &   , LS_CLOUD_WEIGHT_ABSORPTIVITY                                 &
     &   , L_CNV_CLOUD_ABSORPTIVITY, CNV_CLOUD_ABSORPTIVITY             &
     &   , CNV_CLOUD_WEIGHT_ABSORPTIVITY                                &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!     FOR EACH AEROSOL TYPE, COMPUTE THE OPTICAL DEPTH IF
!     IT WAS REQUESTED.
!
      IF (L_AOD_SULPHATE) THEN
! DEPENDS ON: compute_aod
        CALL COMPUTE_AOD(                                               &
     &    N_AEROSOL,      NPD_AEROSOL_SPECIES,                          &
     &    N_AOD_WAVEL,    NPD_AOD_WAVEL,                                &
     &    N_PROFILE,      NPD_PROFILE,                                  &
     &    N_LAYER,        NPD_LAYER,                                    &
     &    NPD_HUMIDITIES, NPD_HUMIDITIES,                               &
     &    TYPE_AEROSOL,                                                 &
     &    I_AOD_TYPE,     IP_TYPE_SULPHATE,                             &
     &    AEROSOL_MIX_RATIO, D_MASS,                                    &
     &    I_AEROSOL_PARAMETRIZATION,                                    &
     &    AOD_ABSORPTION, AOD_SCATTERING,                               &
     &    I_HUMIDITY_POINTER, MEAN_REL_HUMIDITY,                        &
     &    HUMIDITIES,     DELTA_HUMIDITY,                               &
     &    AOD_SULPHATE)
      ENDIF
      IF (L_AOD_DUST) THEN
! DEPENDS ON: compute_aod
        CALL COMPUTE_AOD(                                               &
     &    N_AEROSOL,      NPD_AEROSOL_SPECIES,                          &
     &    N_AOD_WAVEL,    NPD_AOD_WAVEL,                                &
     &    N_PROFILE,      NPD_PROFILE,                                  &
     &    N_LAYER,        NPD_LAYER,                                    &
     &    NPD_HUMIDITIES, NPD_HUMIDITIES,                               &
     &    TYPE_AEROSOL,                                                 &
     &    I_AOD_TYPE,     IP_TYPE_DUST,                                 &
     &    AEROSOL_MIX_RATIO, D_MASS,                                    &
     &    I_AEROSOL_PARAMETRIZATION,                                    &
     &    AOD_ABSORPTION, AOD_SCATTERING,                               &
     &    I_HUMIDITY_POINTER, MEAN_REL_HUMIDITY,                        &
     &    HUMIDITIES,     DELTA_HUMIDITY,                               &
     &    AOD_DUST)
      ENDIF
      IF (L_AOD_SEASALT) THEN
! DEPENDS ON: compute_aod
        CALL COMPUTE_AOD(                                               &
     &    N_AEROSOL,      NPD_AEROSOL_SPECIES,                          &
     &    N_AOD_WAVEL,    NPD_AOD_WAVEL,                                &
     &    N_PROFILE,      NPD_PROFILE,                                  &
     &    N_LAYER,        NPD_LAYER,                                    &
     &    NPD_HUMIDITIES, NPD_HUMIDITIES,                               &
     &    TYPE_AEROSOL,                                                 &
     &    I_AOD_TYPE,     IP_TYPE_SEASALT,                              &
     &    AEROSOL_MIX_RATIO, D_MASS,                                    &
     &    I_AEROSOL_PARAMETRIZATION,                                    &
     &    AOD_ABSORPTION, AOD_SCATTERING,                               &
     &    I_HUMIDITY_POINTER, MEAN_REL_HUMIDITY,                        &
     &    HUMIDITIES,     DELTA_HUMIDITY,                               &
     &    AOD_SEASALT)
      ENDIF
      IF (L_AOD_SOOT) THEN
! DEPENDS ON: compute_aod
        CALL COMPUTE_AOD(                                               &
     &    N_AEROSOL,      NPD_AEROSOL_SPECIES,                          &
     &    N_AOD_WAVEL,    NPD_AOD_WAVEL,                                &
     &    N_PROFILE,      NPD_PROFILE,                                  &
     &    N_LAYER,        NPD_LAYER,                                    &
     &    NPD_HUMIDITIES, NPD_HUMIDITIES,                               &
     &    TYPE_AEROSOL,                                                 &
     &    I_AOD_TYPE,     IP_TYPE_SOOT,                                 &
     &    AEROSOL_MIX_RATIO, D_MASS,                                    &
     &    I_AEROSOL_PARAMETRIZATION,                                    &
     &    AOD_ABSORPTION, AOD_SCATTERING,                               &
     &    I_HUMIDITY_POINTER, MEAN_REL_HUMIDITY,                        &
     &    HUMIDITIES,     DELTA_HUMIDITY,                               &
     &    AOD_SOOT)
      ENDIF
      IF (L_AOD_BIOMASS) THEN
! DEPENDS ON: compute_aod
        CALL COMPUTE_AOD(                                               &
     &    N_AEROSOL,      NPD_AEROSOL_SPECIES,                          &
     &    N_AOD_WAVEL,    NPD_AOD_WAVEL,                                &
     &    N_PROFILE,      NPD_PROFILE,                                  &
     &    N_LAYER,        NPD_LAYER,                                    &
     &    NPD_HUMIDITIES, NPD_HUMIDITIES,                               &
     &    TYPE_AEROSOL,                                                 &
     &    I_AOD_TYPE,     IP_TYPE_BIOMASS,                              &
     &    AEROSOL_MIX_RATIO, D_MASS,                                    &
     &    I_AEROSOL_PARAMETRIZATION,                                    &
     &    AOD_ABSORPTION, AOD_SCATTERING,                               &
     &    I_HUMIDITY_POINTER, MEAN_REL_HUMIDITY,                        &
     &    HUMIDITIES,     DELTA_HUMIDITY,                               &
     &    AOD_BIOMASS)
      ENDIF
      IF (L_AOD_BIOGENIC) THEN
! DEPENDS ON: compute_aod
        CALL COMPUTE_AOD(                                               &
     &    N_AEROSOL,      NPD_AEROSOL_SPECIES,                          &
     &    N_AOD_WAVEL,    NPD_AOD_WAVEL,                                &
     &    N_PROFILE,      NPD_PROFILE,                                  &
     &    N_LAYER,        NPD_LAYER,                                    &
     &    NPD_HUMIDITIES, NPD_HUMIDITIES,                               &
     &    TYPE_AEROSOL,                                                 &
     &    I_AOD_TYPE,     IP_TYPE_BIOGENIC,                             &
     &    AEROSOL_MIX_RATIO, D_MASS,                                    &
     &    I_AEROSOL_PARAMETRIZATION,                                    &
     &    AOD_ABSORPTION, AOD_SCATTERING,                               &
     &    I_HUMIDITY_POINTER, MEAN_REL_HUMIDITY,                        &
     &    HUMIDITIES,     DELTA_HUMIDITY,                               &
     &    AOD_BIOGENIC)
      ENDIF
      IF (L_AOD_OCFF) THEN
! DEPENDS ON: compute_aod
        CALL COMPUTE_AOD(                                               &
     &    N_AEROSOL,      NPD_AEROSOL_SPECIES,                          &
     &    N_AOD_WAVEL,    NPD_AOD_WAVEL,                                &
     &    N_PROFILE,      NPD_PROFILE,                                  &
     &    N_LAYER,        NPD_LAYER,                                    &
     &    NPD_HUMIDITIES, NPD_HUMIDITIES,                               &
     &    TYPE_AEROSOL,                                                 &
     &    I_AOD_TYPE,     IP_TYPE_OCFF,                                 &
     &    AEROSOL_MIX_RATIO, D_MASS,                                    &
     &    I_AEROSOL_PARAMETRIZATION,                                    &
     &    AOD_ABSORPTION, AOD_SCATTERING,                               &
     &    I_HUMIDITY_POINTER, MEAN_REL_HUMIDITY,                        &
     &    HUMIDITIES,     DELTA_HUMIDITY,                               &
     &    AOD_OCFF)
      ENDIF
      IF (L_AOD_DELTA) THEN
! DEPENDS ON: compute_aod
        CALL COMPUTE_AOD(                                               &
     &    N_AEROSOL,      NPD_AEROSOL_SPECIES,                          &
     &    N_AOD_WAVEL,    NPD_AOD_WAVEL,                                &
     &    N_PROFILE,      NPD_PROFILE,                                  &
     &    N_LAYER,        NPD_LAYER,                                    &
     &    NPD_HUMIDITIES, NPD_HUMIDITIES,                               &
     &    TYPE_AEROSOL,                                                 &
     &    I_AOD_TYPE,     IP_TYPE_DELTA,                                &
     &    AEROSOL_MIX_RATIO, D_MASS,                                    &
     &    I_AEROSOL_PARAMETRIZATION,                                    &
     &    AOD_ABSORPTION, AOD_SCATTERING,                               &
     &    I_HUMIDITY_POINTER, MEAN_REL_HUMIDITY,                        &
     &    HUMIDITIES,     DELTA_HUMIDITY,                               &
     &    AOD_DELTA)
      ENDIF
!
!
!     SOLVE THE EQUATION OF TRANSFER IN EACH BAND AND
!     INCREMENT THE FLUXES.
!
      DO I_BAND=I_FIRST_BAND, I_LAST_BAND

!
!
!        DETERMINE WHETHER GASEOUS ABSORPTION IS INCLUDED IN THIS BAND.
         IF ( (L_GAS).AND.(N_BAND_ABSORB(I_BAND) >  0) ) THEN
!
!           NOTE: I_GAS_BAND IS USED EXTENSIVELY BELOW SINCE NESTED
!           ARRAY ELEMENTS IN A SUBROUTINE CALL (SEE LATER) CAN
!           CONFUSE SOME COMPILERS.
!
!           NORMALLY THE NUMBER OF GASES IN THE CALCULATION WILL BE
!           AS IN THE SPECTRAL FILE, BUT PARTICULAR OPTIONS MAY RESULT
!           IN THE OMISSION OF SOME GASES.
!
            N_GAS=N_BAND_ABSORB(I_BAND)
!
            IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_SINGLE) THEN
!
!              THERE WILL BE NO GASEOUS ABSORPTION IN THIS BAND
!              UNLESS THE SELECTED GAS APPEARS.
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
!              SET THE FLAG FOR GASEOUS ABSORPTION IN THE BAND.
               L_GAS_BAND=.TRUE.
!
               DO J=1, N_GAS
!
                  I_GAS_BAND=INDEX_ABSORB(J, I_BAND)
!
!                 RESET THE POINTER IF THERE IS JUST ONE GAS.
!
                  IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_SINGLE)       &
     &               THEN
!                    ONLY THE SELECTED GAS IS ACTIVE IN THE BAND.
                     I_GAS_BAND=I_GAS
!
                  ENDIF
!
                  IF (I_SCALE_ESFT(I_BAND, I_GAS_BAND)                  &
     &                == IP_SCALE_BAND) THEN

                     IEX=1

!                    RESCALE THE AMOUNT OF GAS FOR THIS BAND NOW.
! DEPENDS ON: scale_absorb
                     CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER         &
     &                  , GAS_MIX_RATIO(1, 0, I_GAS_BAND), P, T         &
     &                  , L_LAYER, I_TOP                                &
     &                  , GAS_FRAC_RESCALED(1, 0, I_GAS_BAND)           &
     &                  , I_SCALE_FNC(I_BAND, I_GAS_BAND)               &
     &                  , P_REFERENCE(I_GAS_BAND, I_BAND)               &
     &                  , T_REFERENCE(I_GAS_BAND, I_BAND)               &
     &                  , SCALE_VECTOR(1, 1, I_BAND, I_GAS_BAND)        &
     &                  , L_WENYI, IEX, I_BAND                          & 
     &                  , L_DOPPLER(I_GAS_BAND)                         &
     &                  , DOPPLER_CORRECTION(I_GAS_BAND)                &
     &                  , NPD_PROFILE, NPD_LAYER, NPD_SCALE_FNC         &
     &                  , NPD_SCALE_VARIABLE                            &
     &                  )
                     IF (IERR /= I_NORMAL) RETURN
!
                  ELSE IF (I_SCALE_ESFT(I_BAND, I_GAS_BAND)             &
     &                == IP_SCALE_NULL) THEN
!                    COPY ACROSS THE UNSCALED ARRAY.
                        DO I=1, N_LAYER
                           DO L=1, N_PROFILE
                              GAS_FRAC_RESCALED(L, I, I_GAS_BAND)       &
     &                           =GAS_MIX_RATIO(L, I, I_GAS_BAND)
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
!        RESCALE AMOUNTS OF CONTINUA.
!
         IF (L_CONTINUUM) THEN
            N_CONTINUUM=N_BAND_CONTINUUM(I_BAND)
            DO I=1, N_CONTINUUM
               I_CONTINUUM_POINTER(I)=INDEX_CONTINUUM(I_BAND, I)
               I_CONTINUUM=I_CONTINUUM_POINTER(I)
               K_CONTINUUM_MONO(I_CONTINUUM)                            &
     &            =K_CONTINUUM(I_BAND, I_CONTINUUM)
! DEPENDS ON: rescale_continuum
               CALL RESCALE_CONTINUUM(N_PROFILE, N_LAYER, I_CONTINUUM   &
     &            , P, T, L_LAYER, I_TOP                                &
     &            , DENSITY, MOLAR_DENSITY_WATER, MOLAR_DENSITY_FRN     &
     &            , GAS_MIX_RATIO(1, 0, INDEX_WATER)                    &
     &            , AMOUNT_CONTINUUM(1, 0, I_CONTINUUM)                 &
     &            , I_SCALE_FNC_CONT(I_BAND, I_CONTINUUM)               &
     &            , P_REF_CONTINUUM(I_CONTINUUM, I_BAND)                &
     &            , T_REF_CONTINUUM(I_CONTINUUM, I_BAND)                &
     &            , SCALE_CONTINUUM(1, I_BAND, I_CONTINUUM)             &
     &            , NPD_PROFILE, NPD_LAYER, NPD_SCALE_FNC               &
     &            , NPD_SCALE_VARIABLE                                  &
     &            )
            ENDDO
         ENDIF
!
!
!
!        CALCULATE THE GREY EXTINCTION WITHIN THE BAND.
!
! DEPENDS ON: grey_extinction
         CALL GREY_EXTINCTION(IERR                                      &
     &      , N_PROFILE, N_LAYER, L_LAYER, P, T, DENSITY                &
     &      , L_RESCALE                                                 &
     &      , L_RAYLEIGH, RAYLEIGH_COEFFICIENT(I_BAND)                  &
     &      , L_CONTINUUM, N_CONTINUUM, I_CONTINUUM_POINTER             &
     &      , K_CONTINUUM_MONO, AMOUNT_CONTINUUM                        &
     &      , L_AEROSOL, N_AEROSOL, AEROSOL_MIX_RATIO                   &
     &      , I_AEROSOL_PARAMETRIZATION                                 &
     &      , I_HUMIDITY_POINTER, HUMIDITIES, DELTA_HUMIDITY            &
     &      , MEAN_REL_HUMIDITY                                         &
     &      , AEROSOL_ABSORPTION(1, 1, I_BAND)                          &
     &      , AEROSOL_SCATTERING(1, 1, I_BAND)                          &
     &      , AEROSOL_ASYMMETRY(1, 1, I_BAND)                           &
     &      , L_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE, N_CLOUD_TOP    &
     &      , L_CLOUD_LAYER, I_CLOUD                                    &
     &      , N_CONDENSED, L_CLOUD_CMP, I_PHASE_CMP                     &
     &      , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST(1, 1, I_BAND)     &
     &      , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                   &
     &      , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION)                     &
     &      , I_CLOUD_TYPE                                              &
     &      , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE          &
     &      , FORWARD_SCATTER_FREE                                      &
     &      , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                        &
     &      , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD                    &
     &      , FRAC_CLOUD, L_PC2                                         &
     &      , L_CLOUD_EXTINCTION, CLOUD_EXTINCTION_BAND                 &
     &      , L_CLOUD_ABSORPTIVITY, CLOUD_ABSORPTIVITY_BAND             &
     &      , L_LS_CLOUD_EXTINCTION, LS_CLOUD_EXTINCTION_BAND           &
     &      , L_LS_CLOUD_ABSORPTIVITY, LS_CLOUD_ABSORPTIVITY_BAND       &
     &      , L_CNV_CLOUD_EXTINCTION, CNV_CLOUD_EXTINCTION_BAND         &
     &      , L_CNV_CLOUD_ABSORPTIVITY, CNV_CLOUD_ABSORPTIVITY_BAND     &
     &      , NPD_PROFILE, NPD_LAYER, NPD_CONTINUUM                     &
     &      , NPD_AEROSOL_SPECIES, NPD_HUMIDITIES                       &
     &      , NPD_CLOUD_PARAMETER                                       &
     &      )
         IF (IERR /= I_NORMAL) RETURN
!
!
!
         IF (I_ANGULAR_INTEGRATION == IP_TWO_STREAM) THEN
!
!           RESCALE THE ASYMMETRY AND CALCULATE THE SCATTERING
!           FRACTIONS. (THESE ARE GREY AND MAY BE CALCULATED OUTSIDE
!           A LOOP OVER GASES).
!
            IF (L_RESCALE) THEN
!
!              RESCALE FREE ASYMMETRY:
!
! DEPENDS ON: rescale_asymmetry
               CALL RESCALE_ASYMMETRY(N_PROFILE, 1, N_LAYER             &
     &            , ASYMMETRY_FREE, FORWARD_SCATTER_FREE                &
     &            , NPD_PROFILE, NPD_LAYER, L_PC2                       &
     &            )
!
!
               IF (L_CLOUD) THEN
!
!                 RESCALE CLOUDY ASYMMETRY:
!
                  DO K=1, NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION)
! DEPENDS ON: rescale_asymmetry
                     CALL RESCALE_ASYMMETRY(N_PROFILE, N_CLOUD_TOP      &
     &                  , N_LAYER                                       &
     &                  , ASYMMETRY_CLOUD(1, 1, K)                      &
     &                  , FORWARD_SCATTER_CLOUD(1, 1, K)                &
     &                  , NPD_PROFILE, NPD_LAYER, L_PC2                 &
     &                  )
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
!        PRELIMINARY CALCULATIONS FOR SOURCE TERMS:
!
         IF (ISOLIR == IP_SOLAR) THEN
!           CONVERT NORMALIZED BAND FLUXES TO ACTUAL ENERGY FLUXES.
            DO L=1, N_PROFILE
               INC_SOLAR_FLUX_BAND(L)=SOLAR_TOA(L)                      &
     &            *SOLAR_FLUX_BAND(I_BAND)/SEC_0(L)
            ENDDO
!
         ELSE IF (ISOLIR == IP_INFRA_RED) THEN
!
!           CALCULATE THE CHANGE IN THE THERMAL SOURCE FUNCTION
!           ACROSS EACH LAYER FOR THE INFRA-RED PART OF THE SPECTRUM.
!
! DEPENDS ON: diff_planck_source
            CALL DIFF_PLANCK_SOURCE(N_PROFILE, N_LAYER                  &
     &         , N_DEG_FIT, THERMAL_COEFFICIENT(0, I_BAND)              &
     &         , T_REF_PLANCK, T_LEVEL, T_GROUND                        &
     &         , T_SOLID, T_SEA, L_CTILE                                &
     &         , PLANCK_SOURCE_BAND, DIFF_PLANCK_BAND                   &
     &         , THERMAL_GROUND_BAND                                    &
     &         , L_IR_SOURCE_QUAD, T, DIFF_PLANCK_BAND_2                &
     &         , N_FRAC_SOL_POINT, I_FRAC_SOL_POINT, ICE_FRACTION       &
     &         , FLANDG, PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA            &
     &         , NPD_PROFILE, NPD_LAYER, NPD_THERMAL_COEFF              &
     &         )
         ENDIF
!
!
!
!
!        SET THE SURFACE PROPERTIES:
!
! DEPENDS ON: set_surface_properties
         CALL SET_SURFACE_PROPERTIES(N_POINT_TYPE, INDEX_SURFACE        &
     &      , I_SPEC_SURFACE                                            &
     &      , ISOLIR, I_BAND                                            &
     &      , SURFACE_ALBEDO                                            &
     &      , ALBEDO_FIELD_DIFF(1, I_BAND), ALBEDO_FIELD_DIR(1, I_BAND) &
     &      , N_DIR_ALBEDO_FIT, DIRECT_ALBEDO_PARM, SEC_0               &
     &      , EMISSIVITY_GROUND, EMISSIVITY_FIELD(1, I_BAND)            &
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR                   &
     &      , THERMAL_GROUND_BAND                                       &
     &      , NPD_PROFILE, NPD_BAND, NPD_SURFACE, NPD_ALBEDO_PARM       &
     &      )
!
!
!
!
!
!        CALL A SOLVER APPROPRIATE TO THE PRESENCE OF GASES AND
!        THE OVERLAP ASSUMED:
!
         IF (.NOT.L_GAS_BAND) THEN
!
!           THERE IS NO GASEOUS ABSORPTION. SOLVE FOR THE
!           FLUXES DIRECTLY.
!
! DEPENDS ON: solve_band_without_gas
            CALL SOLVE_BAND_WITHOUT_GAS(IERR                            &
!                       Atmospheric Properties
     &         , N_PROFILE, N_LAYER, D_MASS                             &
!                       Angular integration
     &         , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT   &
     &         , L_RESCALE, N_ORDER_GAUSS                               &
!                       Treatment of scattering
     &         , I_SCATTER_METHOD_BAND(I_BAND)                          &
!                       Options for solver
     &         , I_SOLVER, .FALSE., N_AUGMENT                           &
!                       Spectral region
     &         , ISOLIR                                                 &
!                       Solar Properties
     &         , SEC_0, INC_SOLAR_FLUX_BAND                             &
!                       Infra-red Properties
     &         , PLANCK_SOURCE_BAND(1, 0)                               &
     &         , PLANCK_SOURCE_BAND(1, N_LAYER)                         &
     &         , DIFF_PLANCK_BAND, L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2 &
!                       Surface Properties
     &         , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR                &
     &         , THERMAL_GROUND_BAND                                    &
!                       Clear-sky optical properties
     &         , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE       &
     &         , FORWARD_SCATTER_FREE                                   &
!                       Cloudy properties
     &         , L_CLOUD, I_CLOUD                                       &
!                       Cloudy Geometry
     &         , N_CLOUD_TOP                                            &
     &         , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD      &
     &         , N_REGION, I_REGION_CLOUD, FRAC_REGION                  &
     &         , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE                 &
     &         , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE              &
     &         , CLOUD_OVERLAP                                          &
     &         , N_COLUMN, L_COLUMN, AREA_COLUMN                        &
!                       Cloudy optical properties
     &         , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                     &
     &         , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD                 &
!                       Calculated Fluxes
     &         , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                      &
!                       Flags for Clear-sky Fluxes
     &         , L_CLEAR, I_SOLVER_CLEAR                                &
!                       Calculated Clear-sky Fluxes
     &         , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND          &
!                       Planckian Function
     &         , PLANCK_SOURCE_BAND                                     &
!                       Dimensions of Arrays
     &         , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                     &
     &         )
            IF (IERR /= I_NORMAL) RETURN
!
!
         ELSE
!
!           GASES ARE INCLUDED.
!
!           INITIALIZE THE FLUX IN THE BAND TO ZERO.
! DEPENDS ON: initialize_flux
            CALL INITIALIZE_FLUX(N_PROFILE, N_LAYER, N_AUGMENT          &
     &         , ISOLIR                                                 &
     &         , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                      &
     &         , L_CLEAR                                                &
     &         , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND          &
     &         , 0.0E+00                                                &
     &         , NPD_PROFILE, NPD_LAYER                                 &
     &         , .FALSE.                                                &
     &         )
!
!           TREAT THE GASEOUS OVERLAPS AS DIRECTED BY
!           THE OVERLAP SWITCH.
!
            IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_SINGLE) THEN
!
! DEPENDS ON: solve_band_one_gas
               CALL SOLVE_BAND_ONE_GAS(IERR                             &
!                       Atmospheric Properties
     &            , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS    &
!                       Angular Integration
     &            , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT&
     &            , L_RESCALE, N_ORDER_GAUSS                            &
!                       Treatment of Scattering
     &            , I_SCATTER_METHOD_BAND(I_BAND)                       &
!                       Options for solver
     &            , I_SOLVER, .FALSE., N_AUGMENT                        &
!                       Gaseous Properties
     &            , I_BAND, I_GAS                                       &
     &            , I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC              &
     &            , K_ESFT, W_ESFT, L_WENYI, SCALE_VECTOR               &
     &            , P_REFERENCE, T_REFERENCE                            &
     &            , GAS_MIX_RATIO, GAS_FRAC_RESCALED                    &
     &            , L_DOPPLER, DOPPLER_CORRECTION                       &
!                       Spectral Region
     &            , ISOLIR                                              &
!                       Solar Properties
     &            , SEC_0, INC_SOLAR_FLUX_BAND                          &
!                       Infra-red Properties
     &            , PLANCK_SOURCE_BAND(1, 0)                            &
     &            , PLANCK_SOURCE_BAND(1, N_LAYER)                      &
     &            , DIFF_PLANCK_BAND                                    &
     &            , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                &
!                       Surface Properties
     &            , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR             &
     &            , THERMAL_GROUND_BAND                                 &
!                       Clear-sky optical Properties
     &            , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE    &
     &            , FORWARD_SCATTER_FREE                                &
!                       Cloudy properties
     &            , L_CLOUD, I_CLOUD                                    &
!                       Cloud Geometry
     &            , N_CLOUD_TOP                                         &
     &            , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD   &
     &            , N_REGION, I_REGION_CLOUD, FRAC_REGION               &
     &            , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE              &
     &            , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE           &
     &            , CLOUD_OVERLAP                                       &
     &            , N_COLUMN, L_COLUMN, AREA_COLUMN                     &
!                       Cloudy Optical Properties
     &            , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                  &
     &            , ASYMMETRY_CLOUD                                     &
     &            , FORWARD_SCATTER_CLOUD                               &
!                       Fluxes Calculated
     &            , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                   &
!                       Flags for clear-sky calculations
     &            , L_CLEAR, I_SOLVER_CLEAR                             &
!                       Clear-sky Fluxes
     &            , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND       &
!                       Planckian Function
     &            , PLANCK_SOURCE_BAND                                  &
!                       Dimensions of Arrays
     &            , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                  &
     &            , NPD_BAND, NPD_SPECIES                               &
     &            , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC    &
     &            )
!
            ELSE IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_RANDOM) THEN
!
! DEPENDS ON: solve_band_random_overlap
               CALL SOLVE_BAND_RANDOM_OVERLAP(IERR                      &
!                       Atmospheric Properties
     &            , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS    &
!                       Angular Integration
     &            , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT&
     &            , L_RESCALE, N_ORDER_GAUSS                            &
!                       Treatment of Scattering
     &            , I_SCATTER_METHOD_BAND(I_BAND)                       &
!                       Options for solver
     &            , I_SOLVER, .FALSE., N_AUGMENT                        &
!                       Gaseous Properties
     &            , I_BAND, N_GAS                                       &
     &            , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC&
     &            , K_ESFT, W_ESFT, L_WENYI, SCALE_VECTOR               &
     &            , P_REFERENCE, T_REFERENCE                            &
     &            , GAS_MIX_RATIO, GAS_FRAC_RESCALED                    &
     &            , L_DOPPLER, DOPPLER_CORRECTION                       &
!                       Spectral Region
     &            , ISOLIR                                              &
!                       Solar Properties
     &            , SEC_0, INC_SOLAR_FLUX_BAND                          &
!                       Infra-red Properties
     &            , PLANCK_SOURCE_BAND(1, 0)                            &
     &            , PLANCK_SOURCE_BAND(1, N_LAYER)                      &
     &            , DIFF_PLANCK_BAND                                    &
     &            , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                &
!                       Surface Properties
     &            , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR             &
     &            , THERMAL_GROUND_BAND                                 &
!                       Clear-sky optical Properties
     &            , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE    &
     &            , FORWARD_SCATTER_FREE                                &
!                       Cloudy Properties
     &            , L_CLOUD, I_CLOUD                                    &
!                       Cloud Geometry
     &            , N_CLOUD_TOP                                         &
     &            , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD   &
     &            , N_REGION, I_REGION_CLOUD, FRAC_REGION               &
     &            , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE              &
     &            , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE           &
     &            , CLOUD_OVERLAP                                       &
     &            , N_COLUMN, L_COLUMN, AREA_COLUMN                     &
!                       Cloudy optical Properties
     &            , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                  &
     &            , ASYMMETRY_CLOUD                                     &
     &            , FORWARD_SCATTER_CLOUD                               &
!                       Fluxes Calculated
     &            , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                   &
!                       Flags for Clear-sky Calculations
     &            , L_CLEAR, I_SOLVER_CLEAR                             &
!                       Clear-sky Fluxes
     &            , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND       &
!                       Planckian Function
     &            , PLANCK_SOURCE_BAND                                  &
!                       Dimensions of Arrays
     &            , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                  &
     &            , NPD_BAND, NPD_SPECIES                               &
     &            , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC    &
     &            )
!
            ELSE IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_FESFT) THEN
!
! DEPENDS ON: solve_band_fesft
               CALL SOLVE_BAND_FESFT(IERR                               &
!                       Atmospheric Properties
     &            , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS    &
!                       Angular Integration
     &            , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT&
     &            , L_RESCALE, N_ORDER_GAUSS                            &
!                       Treatment of Scattering
     &            , I_SCATTER_METHOD_BAND(I_BAND)                       &
!                       Options for solver
     &            , I_SOLVER, .FALSE., N_AUGMENT                        &
!                       Gaseous Properties
     &            , I_BAND, N_GAS                                       &
     &            , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC&
     &            , K_ESFT, W_ESFT, L_WENYI, SCALE_VECTOR               &
     &            , P_REFERENCE, T_REFERENCE                            &
     &            , GAS_MIX_RATIO, GAS_FRAC_RESCALED                    &
     &            , L_DOPPLER, DOPPLER_CORRECTION                       &
!                       Spectral Region
     &            , ISOLIR                                              &
!                       Solar Properties
     &            , SEC_0, INC_SOLAR_FLUX_BAND                          &
!                       Infra-red Properties
     &            , PLANCK_SOURCE_BAND(1, 0)                            &
     &            , PLANCK_SOURCE_BAND(1, N_LAYER)                      &
     &            , DIFF_PLANCK_BAND                                    &
     &            , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                &
!                       Surface Properties
     &            , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR             &
     &            , THERMAL_GROUND_BAND                                 &
!                       Clear-sky Optical Properties
     &            , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE    &
     &            , FORWARD_SCATTER_FREE                                &
!                       Cloudy Properties
     &            , L_CLOUD, I_CLOUD                                    &
!                       Cloud Geometry
     &            , N_CLOUD_TOP                                         &
     &            , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD   &
     &            , N_REGION, I_REGION_CLOUD, FRAC_REGION               &
     &            , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE              &
     &            , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE           &
     &            , CLOUD_OVERLAP                                       &
     &            , N_COLUMN, L_COLUMN, AREA_COLUMN                     &
!                       Cloudy Optical Properties
     &            , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                  &
     &            , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD              &
!                       Fluxes Calculated
     &            , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                   &
!                       Flags for Clear Fluxes
     &            , L_CLEAR, I_SOLVER_CLEAR                             &
!                       Clear-sky Fluxes Calculated
     &            , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND       &
!                       Planckian Source Function
     &            , PLANCK_SOURCE_BAND                                  &
!                       Sizes of Arrays
     &            , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                  &
     &            , NPD_BAND, NPD_SPECIES                               &
     &            , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC    &
     &            )
!
            ELSE IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_CLR_FESFT) THEN
!
! DEPENDS ON: solve_band_clr_fesft
               CALL SOLVE_BAND_CLR_FESFT(IERR                           &
!                       Atmospheric Properties
     &            , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS    &
!                       Angular Integration
     &            , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT&
     &            , L_RESCALE, N_ORDER_GAUSS                            &
!                       Treatment of Scattering
     &            , I_SCATTER_METHOD_BAND(I_BAND)                       &
!                       Options for Solver
     &            , I_SOLVER, .FALSE., N_AUGMENT                        &
!                       Gaseous Properties
     &            , I_BAND, N_GAS                                       &
     &            , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC&
     &            , K_ESFT, W_ESFT, L_WENYI, SCALE_VECTOR               &
     &            , P_REFERENCE, T_REFERENCE                            &
     &            , GAS_MIX_RATIO, GAS_FRAC_RESCALED                    &
     &            , L_DOPPLER, DOPPLER_CORRECTION                       &
!                       Spectral Region
     &            , ISOLIR                                              &
!                       Solar Properties
     &            , SEC_0, INC_SOLAR_FLUX_BAND                          &
!                       Infra-red Properties
     &            , PLANCK_SOURCE_BAND                                  &
     &            , DIFF_PLANCK_BAND                                    &
     &            , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                &
!                       Surface Properties
     &            , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR             &
     &            , THERMAL_GROUND_BAND                                 &
!                       Clear-sky Optical Properties
     &            , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE    &
     &            , FORWARD_SCATTER_FREE                                &
!                       Cloudy Properties
     &            , L_CLOUD, I_CLOUD                                    &
!                       Cloud Geometry
     &            , N_CLOUD_TOP                                         &
     &            , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD   &
     &            , N_REGION, I_REGION_CLOUD, FRAC_REGION               &
     &            , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE              &
     &            , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE           &
     &            , CLOUD_OVERLAP                                       &
     &            , N_COLUMN, L_COLUMN, AREA_COLUMN                     &
!                       Cloudy Optical Properties
     &            , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                  &
     &            , ASYMMETRY_CLOUD                                     &
     &            , FORWARD_SCATTER_CLOUD                               &
!                       Fluxes Calculated
     &            , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                   &
!                       Flags for Clear-sky Fluxes
     &            , L_CLEAR, I_SOLVER_CLEAR                             &
!                       Clear-sky Fluxes Calculated
     &            , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND       &
     &            , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                  &
     &            , NPD_BAND, NPD_SPECIES                               &
     &            , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC    &
     &            )
!
            ELSE IF (I_GAS_OVERLAP(I_BAND) == IP_OVERLAP_K_EQV) THEN

! DEPENDS ON: solve_band_k_eqv
               CALL SOLVE_BAND_K_EQV(IERR                               &
!                       Atmospheric Properties
     &            , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS    &
!                       Angular Integration
     &            , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT&
     &            , L_RESCALE, N_ORDER_GAUSS                            &
!                       Treatment of Scattering
     &            , I_SCATTER_METHOD_BAND(I_BAND)                       &
!                       Options for Solver
     &            , I_SOLVER, .FALSE., N_AUGMENT                        &
!                       Gaseous Properties
     &            , I_BAND, N_GAS                                       &
     &            , INDEX_ABSORB, I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC&
     &            , K_ESFT, W_ESFT, L_WENYI, SCALE_VECTOR               &
     &            , P_REFERENCE, T_REFERENCE, L_MOD_K_FLUX              &
     &            , GAS_MIX_RATIO, GAS_FRAC_RESCALED                    &
     &            , L_DOPPLER, DOPPLER_CORRECTION                       &
!                       Spectral Region
     &            , ISOLIR                                              &
!                       Solar Properties
     &            , SEC_0, INC_SOLAR_FLUX_BAND                          &
!                       Infra-red Properties
     &            , PLANCK_SOURCE_BAND                                  &
     &            , DIFF_PLANCK_BAND                                    &
     &            , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                &
!                       Surface Properties
     &            , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR             &
     &            , THERMAL_GROUND_BAND                                 &
!                       Clear-sky optical properties
     &            , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE    &
     &            , FORWARD_SCATTER_FREE                                &
!                       Cloudy Properties
     &            , L_CLOUD, I_CLOUD                                    &
!                       Cloud Geometry
     &            , N_CLOUD_TOP                                         &
     &            , NP_CLOUD_TYPE(I_CLOUD_REPRESENTATION), FRAC_CLOUD   &
     &            , N_REGION, I_REGION_CLOUD, FRAC_REGION               &
     &            , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE              &
     &            , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE           &
     &            , CLOUD_OVERLAP                                       &
     &            , N_COLUMN, L_COLUMN, AREA_COLUMN                     &
!                       Cloudy Optical Properties
     &            , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                  &
     &            , ASYMMETRY_CLOUD                                     &
     &            , FORWARD_SCATTER_CLOUD                               &
!                       Fluxes Calculated
     &            , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                   &
!                       Flags for Clear-sky Calculations
     &            , L_CLEAR, I_SOLVER_CLEAR                             &
!                       Clear-sky Fluxes Calculated
     &            , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND       &
!                       Dimensions of Arrays
     &            , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                  &
     &            , NPD_BAND, NPD_SPECIES                               &
     &            , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC    &
     &            )
!
            ENDIF
         ENDIF
!
!
!
!        INCREMENT THE TOTAL FLUXES.
!
! DEPENDS ON: augment_total_flux
         CALL AUGMENT_TOTAL_FLUX(N_PROFILE, N_LAYER, N_AUGMENT          &
     &      , ISOLIR, L_CLEAR, .FALSE.                                  &
     &      , WEIGHT_BAND(I_BAND), PLANCK_SOURCE_BAND                   &
     &      , FLUX_DIRECT, FLUX_TOTAL                                   &
     &      , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                         &
     &      , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                       &
     &      , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND             &
     &      , PLANCK_FLUX                                               &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      , ALBEDO_SURFACE_DIFF                                       &
     &      , ALBEDO_SURFACE_DIR                                        &
     &      )
!
!
!
!        INCREMENT THE BAND-DEPENDENT DIAGNOSTICS FOR THE
!        UNIFIED MODEL.
! DEPENDS ON: r2_couple_diag
         CALL R2_COUPLE_DIAG(N_PROFILE, N_LAYER, ISOLIR                 &
     &      , ALBEDO_FIELD_DIFF(1, I_BAND), ALBEDO_FIELD_DIR(1, I_BAND) &
     &      , ALBEDO_SEA_DIFF(1, I_BAND), ALBEDO_SEA_DIR(1, I_BAND)     &
     &      , FLANDG, ICE_FRACTION                                      &
     &      , PLANCK_FREEZE_SEA, PLANCK_LEADS_SEA                       &
     &      , PLANCK_SOURCE_BAND(1, N_LAYER), THERMAL_GROUND_BAND       &
     &      , FLUX_TOTAL_BAND(1, 2*N_LAYER+2)                           &
     &      , FLUX_TOTAL_BAND(1, 2*N_LAYER+1)                           &
     &      , FLUX_DIRECT_BAND(1, N_LAYER)                              &
     &      , FLUX_TOTAL_CLEAR_BAND(1, 2*N_LAYER+2)                     &
     &      , FLUX_TOTAL_CLEAR_BAND(1, 2*N_LAYER+1)                     &
     &      , FLUX_DIRECT_CLEAR_BAND(1, N_LAYER)                        &
     &      , WEIGHT_690NM(I_BAND), WEIGHT_UV(I_BAND)                   &
     &      , SEA_FLUX                                                  &
     &      , L_SURFACE_DOWN_FLUX, SURFACE_DOWN_FLUX                    &
     &      , L_SURF_DOWN_CLR, SURF_DOWN_CLR                            &
     &      , L_SURF_UP_CLR, SURF_UP_CLR                                &
     &      , L_FLUX_DIFFUSE, FLUX_DIFFUSE                              &
     &      , L_UVFLUX_DIRECT, UV_FLUX_DIRECT                           &
     &      , L_UVFLUX_UP, UV_FLUX_UP                                   &
     &      , L_UVFLUX_DOWN, UV_FLUX_DOWN                               &
     &      , L_FLUX_BELOW_690NM_SURF                                   &
     &      , FLUX_BELOW_690NM_SURF, FL_SEA_BELOW_690NM_SURF            &
     &      , L_MOSES_II, L_CTILE                                       &
     &      , SURF_VIS_DIR, SURF_VIS_DIF, SURF_NIR_DIR, SURF_NIR_DIF    &
     &      , L_CLOUD_EXTINCTION, CLOUD_EXTINCTION_BAND                 &
     &      , L_LS_CLOUD_EXTINCTION, LS_CLOUD_EXTINCTION_BAND           &
     &      , L_CNV_CLOUD_EXTINCTION, CNV_CLOUD_EXTINCTION_BAND         &
     &      , FLUX_DIRECT_BAND, FLUX_DIRECT_CLEAR_BAND                  &
     &      , N_CLOUD_TOP, N_CLOUD_PROFILE, I_CLOUD_PROFILE             &
     &      , W_CLOUD, FRAC_CLOUD                                       &
     &      , CLOUD_EXTINCTION, CLOUD_WEIGHT_EXTINCTION                 &
     &      , LS_CLOUD_EXTINCTION, LS_CLOUD_WEIGHT_EXTINCTION           &
     &      , CNV_CLOUD_EXTINCTION, CNV_CLOUD_WEIGHT_EXTINCTION         &
     &      , L_CLOUD_ABSORPTIVITY, CLOUD_ABSORPTIVITY_BAND             &
     &      , L_LS_CLOUD_ABSORPTIVITY, LS_CLOUD_ABSORPTIVITY_BAND       &
     &      , L_CNV_CLOUD_ABSORPTIVITY, CNV_CLOUD_ABSORPTIVITY_BAND     &
     &      , FLUX_TOTAL_BAND, FLUX_TOTAL_CLEAR_BAND                    &
     &      , CLOUD_ABSORPTIVITY, CLOUD_WEIGHT_ABSORPTIVITY             &
     &      , LS_CLOUD_ABSORPTIVITY, LS_CLOUD_WEIGHT_ABSORPTIVITY       &
     &      , CNV_CLOUD_ABSORPTIVITY, CNV_CLOUD_WEIGHT_ABSORPTIVITY     &
     &      , NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE                    &
     &      )
!
      ENDDO
!
!
!
!
!     PASS THE CALCULATED FLUXES INTO THE OUTPUT ARRAYS.
!
! DEPENDS ON: assign_flux
      CALL ASSIGN_FLUX(N_PROFILE, N_LAYER                               &
     &   , FLUX_TOTAL, FLUX_TOTAL_CLEAR                                 &
     &   , ISOLIR                                                       &
     &   , PLANCK_FLUX                                                  &
     &   , L_CLEAR, .FALSE.                                             &
     &   , FLUX_DOWN, FLUX_UP, FLUX_DOWN_CLEAR, FLUX_UP_CLEAR           &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
      RETURN
      END SUBROUTINE FLUX_CALC
#endif
