#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the fluxes within the band with no gases.
!
! Method:
!       Gaseous extinction is set to 0 and a monochromatic
!       calculation is performed.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVE_BAND_WITHOUT_GAS(IERR                            &
!                     Atmospheric column
     &  , N_PROFILE, N_LAYER, D_MASS                                    &
!                     Angular integration
     &  , I_ANGULAR_INTEGRATION, I_2STREAM                              &
     &  , N_ORDER_PHASE, L_RESCALE, N_ORDER_GAUSS                       &
     &  , MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC                  &
     &  , ACCURACY_ADAPTIVE, EULER_FACTOR, I_SPH_ALGORITHM              &
     &  , I_SPH_MODE                                                    &
!                     Precalculated angular arrays
     &  , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                      &
!                     Treatment of scattering
     &  , I_SCATTER_METHOD                                              &
!                     Options for solver
     &  , I_SOLVER                                                      &
!                     Spectral region
     &  , ISOLIR                                                        &
!                     Solar properties
     &  , ZEN_0, SOLAR_IRRAD                                            &
!                     Infra-red properties
     &  , PLANCK_FLUX_TOP, PLANCK_FLUX_BOTTOM                           &
     &  , DIFF_PLANCK_BAND, L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2        &
!                     Surface properties
     &  , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB                      &
     &  , F_BRDF, BRDF_SOL, BRDF_HEMI                                   &
     &  , PLANCK_FLUX_GROUND                                            &
!                       Tiling of the surface
     &  , L_TILE, N_POINT_TILE, N_TILE, LIST_TILE, RHO_ALB_TILE         &
     &  , PLANCK_FLUX_TILE                                              &
!                     Optical properties
     &  , K_GREY_TOT_CLR, K_EXT_SCAT_CLR, PHASE_FNC_CLR                 &
     &  , FORWARD_SCATTER_CLR, PHASE_FNC_SOLAR_CLR                      &
     &  , K_GREY_TOT, K_EXT_SCAT                                        &
     &  , PHASE_FNC, FORWARD_SCATTER, PHASE_FNC_SOLAR                   &
!                     Cloudy properties
     &  , L_CLOUD, I_CLOUD                                              &
!                     Cloud geometry
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
!                     Calculated fluxes
     &  , FLUX_DIRECT, FLUX_DOWN, FLUX_UP                               &
!                       Calculated radiances
     &  , I_DIRECT, RADIANCE                                            &
!                       Calculated rate of photolysis
     &  , PHOTOLYSIS                                                    &
!                       Flags for clear-sky fluxes
     &  , L_CLEAR, I_SOLVER_CLEAR                                       &
!                     Calculated clear-sky fluxes
     &  , FLUX_DIRECT_CLEAR, FLUX_DOWN_CLEAR, FLUX_UP_CLEAR             &
!                       Tiled Surface Fluxes
     &  , FLUX_UP_TILE, FLUX_UP_BLUE_TILE                               &
!                       Special Surface Fluxes
     &  , L_BLUE_FLUX_SURF, WEIGHT_BLUE                                 &
     &  , FLUX_DIRECT_BLUE_SURF                                         &
     &  , FLUX_DOWN_BLUE_SURF, FLUX_UP_BLUE_SURF                        &
!                     Dimensions of arrays
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN          &
     &  , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE            &
     &  , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF                    &
     &  , ND_MAX_ORDER, ND_SPH_COEFF                                    &
     &  , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC, ND_VIEWING_LEVEL            &
     &  , ND_DIRECTION, ND_SOURCE_COEFF                                 &
     &  , ND_POINT_TILE, ND_TILE                                        &
     &  )
!
      IMPLICIT NONE
!
! Include Header files
#include "c_kinds.h"
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
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
!           Size allocated for columns per point
     &  , ND_CLOUD_TYPE                                                 &
!           Size allocated for types of clouds
     &  , ND_REGION                                                     &
!           Size allocated for regions of clouds
     &  , ND_OVERLAP_COEFF                                              &
!           Size allocated for cloud overlap coefficients
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
#include "error_pcf3z.h"
#include "angular_integration_pcf3z.h"
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
!           Mode in which the spherical solver is being used
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
!
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
!           Two-stream solver used
!
!                     Spectral region
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Visible or IR
!
!                     Solar properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ZEN_0(ND_PROFILE)                                             &
!           Secants (two-stream) or cosines (spherical harmonics)
!           of the solar zenith angle
     &  , SOLAR_IRRAD(ND_PROFILE)
!           Incident solar irradiance in the band
!
!                     Infra-red properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    PLANCK_FLUX_TOP(ND_PROFILE)                                   &
!           Planck function at bottom of column
     &  , PLANCK_FLUX_BOTTOM(ND_PROFILE)                                &
!           Planck function at bottom of column
     &  , DIFF_PLANCK_BAND(ND_PROFILE, ND_LAYER)                        &
!           Differences in the Planckian function (bottom-top) across
!           layers
     &  , DIFF_PLANCK_BAND_2(ND_PROFILE, ND_LAYER)
!           Twice the second difference of Planckian in band
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Use a quadratic source function
!
!                     Surface properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    PLANCK_FLUX_GROUND(ND_PROFILE)
!           Thermal source at surface in band
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
!           Free absorptive extinction
     &  , K_EXT_SCAT_CLR(ND_PROFILE, ND_LAYER_CLR)                      &
!           Free scattering extinction
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
      REAL  (Real64), INTENT(INOUT) ::                                  &
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
!                       Calculated rates of photolysis
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    PHOTOLYSIS(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Rates of photolysis
!
!                     Flags for clear-sky fluxes
      LOGICAL, INTENT(IN) ::                                            &
     &    L_CLEAR
!           Calculate net clear-sky properties
      INTEGER, INTENT(IN) ::                                            &
     &    I_SOLVER_CLEAR
!           Clear solver used
!
!                     Calculated clear-sky fluxes
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT_CLEAR(ND_FLUX_PROFILE, 0: ND_LAYER)               &
!           Clear-sky direct flux
     &  , FLUX_DOWN_CLEAR(ND_FLUX_PROFILE, 0: ND_LAYER)                 &
!           Clear-sky total downward flux
     &  , FLUX_UP_CLEAR(ND_FLUX_PROFILE, 0: ND_LAYER)                   &
!           Clear-sky upward flux
     &  , FLUX_UP_TILE(ND_POINT_TILE, ND_TILE)                          &
!           Upward fluxes at tiled surface points
     &  , FLUX_UP_BLUE_TILE(ND_POINT_TILE, ND_TILE)
!           Upward blue fluxes at tiled surface points
!
!                        Special Diagnostics:
      LOGICAL, INTENT(IN) ::                                            &
     &    L_BLUE_FLUX_SURF
!           Flag to calculate blue fluxes at the surface
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
     &    I                                                             &
!           Loop variable
     &  , L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    FLUX_INC_DIRECT(ND_PROFILE)                                   &
!           Incident direct flux
     &  , FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident downward flux
     &  , D_PLANCK_FLUX_SURFACE(ND_PROFILE)                             &
!           Ground source function
     &  , K_NULL(ND_PROFILE, ND_LAYER)                                  &
!           Null vector for call to subroutine
     &  , DUMMY_KE(ND_PROFILE, ND_LAYER)
!           Dummy array (not used)
!
!     Monochromatic incrementing radiances:
      REAL  (Real64) ::                                                 &
     &    FLUX_DIRECT_BAND(ND_FLUX_PROFILE, 0: ND_LAYER)                &
!           Increment to direct flux
     &  , FLUX_TOTAL_BAND(ND_FLUX_PROFILE, 2*ND_LAYER+2)                &
!           Increment to total flux
     &  , FLUX_DIRECT_CLEAR_BAND(ND_FLUX_PROFILE, 0: ND_LAYER)          &
!           Increment to clear direct flux
     &  , FLUX_TOTAL_CLEAR_BAND(ND_FLUX_PROFILE, 2*ND_LAYER+2)
!           Increment to clear total flux
!                       Increments to Radiances
      REAL  (Real64) ::                                                 &
     &    I_DIRECT_BAND(ND_RADIANCE_PROFILE, 0: ND_LAYER)               &
!           Increments to the solar irradiance
     &  , RADIANCE_BAND(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL           &
     &      , ND_DIRECTION)
!           Increments to the radiance
      REAL  (Real64) ::                                                 &
     &    PHOTOLYSIS_BAND(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Increments to the rate of photolysis
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    MONOCHROMATIC_RADIANCE
!
!
!
!     Set the appropriate total upward and downward fluxes
!     at the boundaries.
!
      IF ( (I_ANGULAR_INTEGRATION == IP_TWO_STREAM).OR.                 &
     &     (I_ANGULAR_INTEGRATION == IP_IR_GAUSS) ) THEN
        IF (ISOLIR == IP_SOLAR) THEN
!         Visible region.
          DO L=1, N_PROFILE
            D_PLANCK_FLUX_SURFACE(L)=0.0E+00_Real64
            FLUX_INC_DOWN(L)=SOLAR_IRRAD(L)/ZEN_0(L)
            FLUX_INC_DIRECT(L)=SOLAR_IRRAD(L)/ZEN_0(L)
          ENDDO
        ELSEIF (ISOLIR == IP_INFRA_RED) THEN
!         Infra-red region.
          DO L=1, N_PROFILE
            FLUX_INC_DIRECT(L)=0.0E+00_Real64
            FLUX_DIRECT_BAND(L, N_LAYER)=0.0E+00_Real64
            FLUX_INC_DOWN(L)=-PLANCK_FLUX_TOP(L)
            D_PLANCK_FLUX_SURFACE(L)                                    &
     &        =PLANCK_FLUX_GROUND(L)-PLANCK_FLUX_BOTTOM(L)
          ENDDO
          IF (L_CLEAR) THEN
            DO L=1, N_PROFILE
              FLUX_DIRECT_CLEAR_BAND(L, N_LAYER)=0.0E+00_Real64
            ENDDO
          ENDIF
        ENDIF
!
      ELSE IF (I_ANGULAR_INTEGRATION == IP_SPHERICAL_HARMONIC) THEN
!
        IF (ISOLIR == IP_SOLAR) THEN
          DO L=1, N_PROFILE
            I_DIRECT_BAND(L, 0)=SOLAR_IRRAD(L)
            FLUX_INC_DOWN(L)=0.0E+00_Real64
          ENDDO
        ELSE
          DO L=1, N_PROFILE
            FLUX_INC_DOWN(L)=-PLANCK_FLUX_TOP(L)
            D_PLANCK_FLUX_SURFACE(L)                                    &
     &        =PLANCK_FLUX_GROUND(L)-PLANCK_FLUX_BOTTOM(L)
          ENDDO
        ENDIF
!
      ENDIF
!
      DO I=1, N_LAYER
        DO L=1, N_PROFILE
          K_NULL(L, I)=0.0E+00_Real64
        ENDDO
      ENDDO
!
!
! DEPENDS ON: monochromatic_radiance
      CALL MONOCHROMATIC_RADIANCE(IERR                                  &
!                     Atmospheric properties
     &  , N_PROFILE, N_LAYER, D_MASS                                    &
!                     Angular integration
     &  , I_ANGULAR_INTEGRATION, I_2STREAM                              &
     &  , L_RESCALE, N_ORDER_GAUSS                                      &
     &  , N_ORDER_PHASE, MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC   &
     &  , ACCURACY_ADAPTIVE, EULER_FACTOR, I_SPH_ALGORITHM              &
     &  , I_SPH_MODE                                                    &
!                     Precalculated angular arrays
     &  , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                      &
!                     Treatment of scattering
     &  , I_SCATTER_METHOD                                              &
!                     Options for solver
     &  , I_SOLVER                                                      &
!                     Gaseous propreties
     &  , K_NULL                                                        &
!                     Options for equivalent extinction
     &  , .FALSE., DUMMY_KE                                             &
!                     Spectral region
     &  , ISOLIR                                                        &
!                     Infra-red properties
     &  , DIFF_PLANCK_BAND, L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2        &
!                     Conditions at TOA
     &  , ZEN_0, FLUX_INC_DIRECT, FLUX_INC_DOWN                         &
     &  , I_DIRECT_BAND                                                 &
!                     Surface properties
     &  , D_PLANCK_FLUX_SURFACE                                         &
     &  , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB                      &
     &  , F_BRDF, BRDF_SOL, BRDF_HEMI                                   &
!                     Optical properties
     &  , K_GREY_TOT_CLR, K_EXT_SCAT_CLR                                &
     &  , PHASE_FNC_CLR, FORWARD_SCATTER_CLR, PHASE_FNC_SOLAR_CLR       &
     &  , K_GREY_TOT, K_EXT_SCAT, PHASE_FNC, FORWARD_SCATTER            &
     &  , PHASE_FNC_SOLAR                                               &
!                     Cloudy properties
     &  , L_CLOUD, I_CLOUD                                              &
!                     Cloud geometry
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
!                     Calculated Flxues
     &  , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                             &
!                       Calculated radiances
     &  , RADIANCE_BAND                                                 &
!                       Calculated rate of photolysis
     &  , PHOTOLYSIS_BAND                                               &
!                     Flags for clear-sky calculations
     &  , L_CLEAR, I_SOLVER_CLEAR                                       &
!                     Clear-sky fluxes calculated
     &  , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND                 &
!                     Dimensions of arrays
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN          &
     &  , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE            &
     &  , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF                    &
     &  , ND_MAX_ORDER, ND_SPH_COEFF                                    &
     &  , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC, ND_VIEWING_LEVEL            &
     &  , ND_DIRECTION, ND_SOURCE_COEFF                                 &
     &  )
!
!     Add the increments to the cumulative fluxes.
! DEPENDS ON: augment_radiance
      CALL AUGMENT_RADIANCE(N_PROFILE, N_LAYER                          &
     &  , I_ANGULAR_INTEGRATION, I_SPH_MODE                             &
     &  , N_VIEWING_LEVEL, N_DIRECTION                                  &
     &  , ISOLIR, L_CLEAR                                               &
     &  , L_INITIAL, WEIGHT_BAND                                        &
     &  , L_BLUE_FLUX_SURF, WEIGHT_BLUE                                 &
!                       Actual radiances
     &  , FLUX_DIRECT, FLUX_DOWN, FLUX_UP                               &
     &  , FLUX_DIRECT_BLUE_SURF                                         &
     &  , FLUX_DOWN_BLUE_SURF, FLUX_UP_BLUE_SURF                        &
     &  , I_DIRECT, RADIANCE, PHOTOLYSIS                                &
     &  , FLUX_DIRECT_CLEAR, FLUX_DOWN_CLEAR, FLUX_UP_CLEAR             &
!                       Increments to radiances
     &  , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                             &
     &  , I_DIRECT_BAND, RADIANCE_BAND, PHOTOLYSIS_BAND                 &
     &  , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND                 &
!                       Dimensions
     &  , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE            &
     &  , ND_LAYER, ND_VIEWING_LEVEL, ND_DIRECTION                      &
     &  )
!
!     Add in the increments from surface tiles
      IF (L_TILE) THEN
! DEPENDS ON: augment_tiled_radiance
        CALL AUGMENT_TILED_RADIANCE(IERR                                &
     &    , N_POINT_TILE, N_TILE, LIST_TILE                             &
     &    , I_ANGULAR_INTEGRATION, ISOLIR, L_INITIAL                    &
     &    , WEIGHT_BAND, L_BLUE_FLUX_SURF, WEIGHT_BLUE                  &
!                       Surface characteristics
     &    , RHO_ALB_TILE                                                &
!                       Actual radiances
     &    , FLUX_UP_TILE, FLUX_UP_BLUE_TILE                             &
!                       Increments to radiances
     &    , FLUX_DIRECT_BAND(1, N_LAYER)                                &
     &    , FLUX_TOTAL_BAND(1, 2*N_LAYER+2)                             &
     &    , PLANCK_FLUX_TILE, PLANCK_FLUX_BOTTOM                        &
!                       Dimensions
     &    , ND_FLUX_PROFILE, ND_POINT_TILE, ND_TILE                     &
     &    , ND_BRDF_BASIS_FNC                                           &
     &    )
        IF (IERR /= I_NORMAL) RETURN
      ENDIF
!
!     After the first call to these routines quantities should be
!     incremented rather than initialized, until the flag is reset.
      L_INITIAL=.FALSE.
!
!
!
      RETURN
      END SUBROUTINE SOLVE_BAND_WITHOUT_GAS
#endif
