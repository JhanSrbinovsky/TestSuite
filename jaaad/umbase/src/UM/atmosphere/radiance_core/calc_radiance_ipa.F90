#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate monochromatic radiances using IPA.
!
! Method:
!
!   In this subroutine a long vector for radiance calculations
!   is set up using the information on the types of cloud present.
!
! Current owner of code: James Manners
!
! Description of code:
!   Fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_RADIANCE_IPA(IERR                                 &
!                     Atmospheric Properties
     &  , N_PROFILE, N_LAYER, N_CLOUD_TOP                               &
!                       Angular Integration
     &  , N_ORDER_PHASE, MS_MIN, MS_MAX, LS_LOCAL_TRUNC                 &
     &  , I_TRUNCATION, ACCURACY_ADAPTIVE, EULER_FACTOR                 &
     &  , I_SPH_ALGORITHM, I_SPH_MODE, L_RESCALE                        &
!                       Precalculated angular arrays
     &  , IA_SPH_MM, CG_COEFF, UPLM_ZERO, UPLM_SOL                      &
!                     Options for Equivalent Extinction
     &  , L_SCALE_SOLAR, ADJUST_SOLAR_KE                                &
!                     Spectral Region
     &  , ISOLIR                                                        &
!                     Infra-red Properties
     &  , DIFF_PLANCK                                                   &
     &  , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                               &
!                     Conditions at TOA
     &  , FLUX_INC_DOWN, ZEN_0                                          &
!                     Conditions at Surface
     &  , D_PLANCK_FLUX_SURFACE                                         &
     &  , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB                      &
     &  , F_BRDF, BRDF_SOL, BRDF_HEMI                                   &
!                     Single Scattering Properties
     &  , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR, PHASE_FNC_SOLAR_CLR        &
     &  , FORWARD_SCATTER_CLR                                           &
     &  , TAU, OMEGA, PHASE_FNC, PHASE_FNC_SOLAR, FORWARD_SCATTER       &
!                     Cloud Geometry
     &  , N_COLUMN_SLV, LIST_COLUMN_SLV                                 &
     &  , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN                     &
!                       Levels for calculating radiances
     &  , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                  &
!                       Viewing Geometry
     &  , N_DIRECTION, DIRECTION                                        &
!                       Calculated fluxes or radiances
     &  , FLUX_DIRECT, FLUX_TOTAL, I_DIRECT, RADIANCE, J_RADIANCE       &
!                     Dimensions of Arrays
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT                     &
     &  , ND_COLUMN, ND_CLOUD_TYPE                                      &
     &  , ND_FLUX_PROFILE, ND_RADIANCE_PROFILE, ND_J_PROFILE            &
     &  , ND_MAX_ORDER, ND_SPH_COEFF                                    &
     &  , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC                              &
     &  , ND_RED_EIGENSYSTEM, ND_SPH_EQUATION, ND_SPH_DIAGONAL          &
     &  , ND_SPH_CF_WEIGHT, ND_SPH_U_RANGE                              &
     &  , ND_VIEWING_LEVEL, ND_DIRECTION                                &
     &  , ND_PROFILE_COLUMN                                             &
     &  )
!
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
!           Size allocated for atmospheric profiles
     &  , ND_FLUX_PROFILE                                               &
!           Size allocated for profiles of output fluxes
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for profiles of radiances
     &  , ND_J_PROFILE                                                  &
!           Size allocated for profiles of photolysis rates
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for completely clear layers
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_COLUMN                                                     &
!           Size allocated for columns at a grid-point
     &  , ND_CLOUD_TYPE                                                 &
!           Size allocated for types of clouds
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
!           Size allowed for the spherical harmonic eigensystem
     &  , ND_SPH_EQUATION                                               &
!           Size allowed for spherical harmonic equations
     &  , ND_SPH_DIAGONAL                                               &
!           Size allowed for diagonals of the spherical harmonic
!           matrix
     &  , ND_SPH_CF_WEIGHT                                              &
!           Size allowed for application of weights of the C. F.
     &  , ND_SPH_U_RANGE                                                &
!           Size allowed for range of values of u^+|- contributing
!           on any viewing level
     &  , ND_DIRECTION                                                  &
!           Size allocated for viewing dierctions
     &  , ND_PROFILE_COLUMN                                             &
!           Number of profiles of subcolumns considered at once
     &  , ID_CT
!           Topmost declared cloudy layer
!
!     Include header files.
#include "def_std_io_icf3z.h"
#include "spectral_region_pcf3z.h"
#include "sph_mode_pcf3z.h"
#include "sph_algorithm_pcf3z.h"
#include "solver_pcf3z.h"
#include "error_pcf3z.h"
!
!
!     Dummy arguments.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , N_CLOUD_TOP                                                   &
!           Topmost cloudy layer
     &  , N_ORDER_PHASE
!           Number of orders retained in the phase function
!
!                     Spherical arrays
      INTEGER, INTENT(IN) ::                                            &
     &    MS_MIN                                                        &
!           Lowest azimuthal order calculated
     &  , MS_MAX                                                        &
!           Highest azimuthal order calculated
     &  , I_TRUNCATION                                                  &
!           Type of speherical truncation
     &  , I_SPH_MODE                                                    &
!           Mode in which the spherical harmonic solver is used
     &  , I_SPH_ALGORITHM                                               &
!           Spherical harmonic algorithm
     &  , IA_SPH_MM(0: ND_MAX_ORDER)                                    &
!           Address of spherical coefficient for (m, m) for each m
     &  , LS_LOCAL_TRUNC(0: ND_MAX_ORDER)
!           Orders of truncation at each azimuthal order
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
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Spectral region
      LOGICAL, INTENT(IN) ::                                            &
     &    L_SCALE_SOLAR                                                 &
!           Scale solar beam
     &  , L_IR_SOURCE_QUAD                                              &
!           Use a quadratic source term
     &  , L_RESCALE
!           Flag for rescaling of the optical properties
!
!     Fields for equivalent extinction
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ADJUST_SOLAR_KE(ND_PROFILE, ND_LAYER)
!           Adjustment of solar beam with equivalent extinction
!
!     Clear-sky optical properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU_CLR(ND_PROFILE, ND_LAYER_CLR)                             &
!           Clear-sky optical depth
     &  , OMEGA_CLR(ND_PROFILE, ND_LAYER_CLR)                           &
!           Clear-sky albedos of single scattering
     &  , PHASE_FNC_CLR(ND_PROFILE, ND_LAYER_CLR, ND_MAX_ORDER)         &
!           Moments of the clear-sky phase function
     &  , PHASE_FNC_SOLAR_CLR(ND_RADIANCE_PROFILE                       &
     &      , ND_LAYER_CLR, ND_DIRECTION)                               &
!           Clear-sky solar phase functions
     &  , FORWARD_SCATTER_CLR(ND_PROFILE, ND_LAYER_CLR)
!           Clear-sky forward scattering fractions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)            &
!           Optical depth
     &  , OMEGA(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)          &
!           Single scattering albedos
     &  , PHASE_FNC(ND_PROFILE, ID_CT: ND_LAYER                         &
     &      , ND_MAX_ORDER, 0: ND_CLOUD_TYPE)                           &
!           Moments of the phase functions
     &  , PHASE_FNC_SOLAR(ND_RADIANCE_PROFILE, ID_CT: ND_LAYER          &
     &      , ND_DIRECTION, 0: ND_CLOUD_TYPE)                           &
!           Solar phase functions evaluated in the viewing directions
     &  , FORWARD_SCATTER(ND_PROFILE, ID_CT: ND_LAYER                   &
     &     , 0: ND_CLOUD_TYPE)
!           Forward scattering fractions in clouds
!
!     Planckian terms:
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DIFF_PLANCK(ND_PROFILE, ND_LAYER)                             &
!           Change in Planckian function
     &  , DIFF_PLANCK_2(ND_PROFILE, ND_LAYER)                           &
!           Twice 2nd differences in Planckian
     &  , D_PLANCK_FLUX_SURFACE(ND_PROFILE)
!           Differential Planckian flux from the surface
!
!     Conditions at TOA
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ZEN_0(ND_PROFILE)                                             &
!           Secant of zenith angle
     &  , FLUX_INC_DOWN(ND_PROFILE)
!           Incident total flux
!
!     Conditions at surface
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
     &    AREA_COLUMN(ND_PROFILE, ND_COLUMN)
!           Area of each column
!
!                     Levels at which radiances will be calculated
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
!                       Calculated Fluxes or Radiances
      REAL  (Real64), INTENT(INOUT) ::                                  &
     &    I_DIRECT(ND_PROFILE, 0: ND_LAYER)
!           Direct radiances
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_FLUX_PROFILE, 0: ND_LAYER)                     &
!           Direct Flux
     &  , FLUX_TOTAL(ND_FLUX_PROFILE, 2*ND_LAYER+2)                     &
!           Total Fluxes
     &  , RADIANCE(ND_RADIANCE_PROFILE, ND_VIEWING_LEVEL, ND_DIRECTION) &
!           Radiances
     &  , J_RADIANCE(ND_J_PROFILE, ND_VIEWING_LEVEL)
!           Photolysis rates
!
!
!
!     Local variables.
      INTEGER                                                           &
     &    LS                                                            &
!           Polar order of harmonic
     &  , MS                                                            &
!           Azimuthal order of harmonic
     &  , I                                                             &
!           Loop variable
     &  , J                                                             &
!           Loop variable
     &  , JS                                                            &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , ID                                                            &
!           Loop variable
     &  , LP                                                            &
!           Index of current real grid-point during assignments
     &  , LL                                                            &
!           Index in the long array of columns to be taken in one go
     &  , LL_COPY                                                       &
!           Index of column to be copied
     &  , ICL                                                           &
!           Index of notional sub-column
     &  , ICS                                                           &
!           Index of current sub-column where a solution is required
     &  , ICC                                                           &
!           Temporary variable listing the layer in the current column
!           where a change is required
     &  , ICT
!           Temporary variable listing the type of optical region moved
!           into be the current change
      INTEGER                                                           &
     &    N_LONG                                                        &
!           Length of long vector
     &  , TARGET(ND_PROFILE_COLUMN)
!           Actual target grid-point for point in the long array
      REAL  (Real64) ::                                                 &
     &    WEIGHT_COLUMN(ND_PROFILE_COLUMN)
!           Weight applied to each column in the sum
      LOGICAL                                                           &
     &    L_NEW
!           Flag to consider a new grid-point
!
!     Properties of vectors of subcolumns
      REAL  (Real64) ::                                                 &
     &    TAU_LONG(ND_PROFILE_COLUMN, ND_LAYER)                         &
!           Long vector of optical depth
     &  , OMEGA_LONG(ND_PROFILE_COLUMN, ND_LAYER)                       &
!           Long vector of albedo of single scattering
     &  , PHASE_FNC_LONG(ND_PROFILE_COLUMN, ND_LAYER, ND_MAX_ORDER)     &
!           Long vector of phase functions
     &  , PHASE_FNC_SOLAR_LONG(ND_PROFILE_COLUMN                        &
     &      , ND_LAYER, ND_DIRECTION)                                   &
!           Long vector of solar phase functions
     &  , FORWARD_SCATTER_LONG(ND_PROFILE_COLUMN, ND_LAYER)             &
!           Long vector of forward scattering fractions
     &  , ADJUST_SOLAR_KE_LONG(ND_PROFILE_COLUMN, ND_LAYER)             &
!           Long vector of solar scalings
     &  , ZEN_0_LONG(ND_PROFILE_COLUMN)                                 &
!           Long vector of cosines of the solar zenith angle
     &  , UPLM_SOL_LONG(ND_PROFILE_COLUMN, ND_SPH_COEFF)                &
!           Long vector of spherical harmonics at the solar angle
     &  , DIFF_PLANCK_LONG(ND_PROFILE_COLUMN, ND_LAYER)                 &
!           Long vector of differences in the Planckian
     &  , DIFF_PLANCK_2_LONG(ND_PROFILE_COLUMN, ND_LAYER)               &
!           Long vector of second differences in the Planckian
     &  , FLUX_INC_DOWN_LONG(ND_PROFILE_COLUMN)                         &
!           Long vector of incident downward fluxes
     &  , D_PLANCK_FLUX_SURFACE_LONG(ND_PROFILE_COLUMN)                 &
!           Long vector of differential Planckian fluxes
!           at the surface
     &  , RHO_ALB_LONG(ND_PROFILE_COLUMN, ND_BRDF_BASIS_FNC)            &
!           Long vector of weightings of BRDF basis functions
     &  , BRDF_SOL_LONG(ND_PROFILE_COLUMN, ND_BRDF_BASIS_FNC            &
     &      , ND_DIRECTION)                                             &
!           The BRDF evaluated for scattering from the solar
!           beam into the viewing direction
     &  , BRDF_HEMI_LONG(ND_PROFILE_COLUMN, ND_BRDF_BASIS_FNC           &
     &      , ND_DIRECTION)                                             &
!           The BRDF evaluated for scattering from isotropic
!           radiation into the viewing direction
     &  , DIRECTION_LONG(ND_PROFILE_COLUMN, ND_DIRECTION, 2)
!           Viewing directions
!
!     Calculated Fluxes or Radiances in subcolumns
      REAL  (Real64) ::                                                 &
     &    FLUX_DIRECT_COLUMN(ND_PROFILE_COLUMN, 0: ND_LAYER)            &
!           Direct Flux
     &  , FLUX_TOTAL_COLUMN(ND_PROFILE_COLUMN, 2*ND_LAYER+2)            &
!           Total Fluxes
     &  , I_DIRECT_COLUMN(ND_PROFILE_COLUMN, 0: ND_LAYER)               &
!           Direct radiances
     &  , RADIANCE_COLUMN(ND_PROFILE_COLUMN, ND_VIEWING_LEVEL           &
     &      , ND_DIRECTION)                                             &
!           Radiances
     &  , PHOTOLYSIS_COLUMN(ND_PROFILE_COLUMN, ND_VIEWING_LEVEL)
!           Photolysis rates
!
!
!
!
!     Functions called:
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    SPH_SOLVER
!
!
!
!     Zero the output arrays ready for incrementing.
      IF (I_SPH_MODE == IP_SPH_MODE_FLUX) THEN
!
        DO I=1, 2*N_LAYER+2
          DO L=1, N_PROFILE
            FLUX_TOTAL(L, I)=0.0E+00_Real64
          ENDDO
        ENDDO
!
        IF (ISOLIR == IP_SOLAR) THEN
          DO I=0, N_LAYER
            DO L=1, N_PROFILE
              FLUX_DIRECT(L, I)=0.0E+00_Real64
            ENDDO
          ENDDO
        ENDIF
!
      ELSE
!
        DO ID=1, N_DIRECTION
          DO I=1, N_VIEWING_LEVEL
            DO L=1, N_PROFILE
              RADIANCE(L, I, ID)=0.0E+00_Real64
            ENDDO
          ENDDO
        ENDDO
!
        IF (ISOLIR == IP_SOLAR) THEN
!         The top level contains the input: other values are zeroed
!         to allow incrementing.
          DO I=1, N_LAYER
            DO L=1, N_PROFILE
              I_DIRECT(L, I)=0.0E+00_Real64
            ENDDO
          ENDDO
        ENDIF
!
!
      ENDIF
!
!
!     Start feeding points into the long array. This is
!     not written to vectorize as that is quite complicated.
!
      LP=1
      L_NEW=.TRUE.
!
      DO WHILE (LP <= N_PROFILE)
!
        LL=0
!
        DO WHILE ( (LL <  ND_PROFILE_COLUMN).AND.(LP <= N_PROFILE) )
!
          LL=LL+1
          TARGET(LL)=LP
!
          IF (L_NEW) THEN
!
!           We consider a new grid-point and so must set the first
!           notional column which is contains no cloud.
            ICL=1
            ICS=1
            DO I=1, N_CLOUD_TOP-1
              TAU_LONG(LL, I)=TAU_CLR(LP, I)
              OMEGA_LONG(LL, I)=OMEGA_CLR(LP, I)
              DO LS=1, N_ORDER_PHASE
                PHASE_FNC_LONG(LL, I, LS)=PHASE_FNC_CLR(LP, I, LS)
              ENDDO
            ENDDO
            DO I=N_CLOUD_TOP, N_LAYER
              TAU_LONG(LL, I)=TAU(LP, I, 0)
              OMEGA_LONG(LL, I)=OMEGA(LP, I, 0)
              DO LS=1, N_ORDER_PHASE
                PHASE_FNC_LONG(LL, I, LS)=PHASE_FNC(LP, I, LS, 0)
              ENDDO
            ENDDO
            IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
              DO ID=1, N_DIRECTION
                DO I=1, N_CLOUD_TOP-1
                  PHASE_FNC_SOLAR_LONG(LL, I, ID)                       &
     &              =PHASE_FNC_SOLAR_CLR(LP, I, ID)
                ENDDO
                DO I=N_CLOUD_TOP, N_LAYER
                  PHASE_FNC_SOLAR_LONG(LL, I, ID)                       &
     &              =PHASE_FNC_SOLAR(LP, I, ID, 0)
                ENDDO
              ENDDO
              IF (L_RESCALE) THEN
                DO I=1, N_CLOUD_TOP-1
                  FORWARD_SCATTER_LONG(LL, I)                           &
     &              =FORWARD_SCATTER_CLR(LP, I)
                ENDDO
                DO I=N_CLOUD_TOP, N_LAYER
                  FORWARD_SCATTER_LONG(LL, I)                           &
     &              =FORWARD_SCATTER(LP, I, 0)
                ENDDO
              ENDIF
            ENDIF
!
!
            L_NEW=.FALSE.
!
!
          ELSE
!
!           Copy the previous column over. Normally this will be the
!           previous one, but if we are starting a new batch it will
!           be the one at the end of the previous batch.
            IF (LL >  1) THEN
              LL_COPY=LL-1
            ELSE
              LL_COPY=N_LONG
            ENDIF
!
            DO I=1, N_LAYER
              TAU_LONG(LL, I)=TAU_LONG(LL_COPY, I)
              OMEGA_LONG(LL, I)=OMEGA_LONG(LL_COPY, I)
              DO LS=1, N_ORDER_PHASE
                PHASE_FNC_LONG(LL, I, LS)                               &
     &            =PHASE_FNC_LONG(LL_COPY, I, LS)
              ENDDO
            ENDDO
            IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
              DO ID=1, N_DIRECTION
                DO I=1, N_LAYER
                  PHASE_FNC_SOLAR_LONG(LL, I, ID)                       &
     &              =PHASE_FNC_SOLAR_LONG(LL_COPY, I, ID)
                ENDDO
              ENDDO
              IF (L_RESCALE) THEN
                DO I=1, N_LAYER
                  FORWARD_SCATTER_LONG(LL, I)                           &
     &              =FORWARD_SCATTER_LONG(LL_COPY, I)
                ENDDO
              ENDIF
            ENDIF
!
          ENDIF
!
!         Move through the notional columns at this grid-point
!         adjusting individiual layers until we find one where the
!         equations are to be solved.
          DO WHILE (ICL <  LIST_COLUMN_SLV(LP, ICS))
            ICC=I_CLM_LYR_CHN(LP, ICL)
            ICT=I_CLM_CLD_TYP(LP, ICL)
!
            TAU_LONG(LL, ICC)=TAU(LP, ICC, ICT)
            OMEGA_LONG(LL, ICC)=OMEGA(LP, ICC, ICT)
            DO LS=1, N_ORDER_PHASE
              PHASE_FNC_LONG(LL, ICC, LS)                               &
     &          =PHASE_FNC(LP, ICC, LS, ICT)
            ENDDO
            IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
              DO ID=1, N_DIRECTION
                PHASE_FNC_SOLAR_LONG(LL, ICC, ID)                       &
     &            =PHASE_FNC_SOLAR(LP, ICC, ID, ICT)
              ENDDO
              IF (L_RESCALE) THEN
                FORWARD_SCATTER_LONG(LL, ICC)                           &
     &            =FORWARD_SCATTER(LP, ICC, ICT)
              ENDIF
            ENDIF
!
            ICL=ICL+1
          ENDDO
!
!
!         Set arrays which are independent of cloud changes.
          IF (ISOLIR == IP_SOLAR) THEN
!
            IF (L_SCALE_SOLAR) THEN
              DO I=1, N_LAYER
                ADJUST_SOLAR_KE_LONG(LL, I)=ADJUST_SOLAR_KE(LP, I)
              ENDDO
            ENDIF
!
            ZEN_0_LONG(LL)=ZEN_0(LP)
            I_DIRECT_COLUMN(LL, 0)=I_DIRECT(LP, 0)
            DO MS=MS_MIN, MS_MAX
              DO LS=MS, LS_LOCAL_TRUNC(MS)+1
                JS=IA_SPH_MM(MS)+LS-MS
                UPLM_SOL_LONG(LL, JS)=UPLM_SOL(LP, JS)
              ENDDO
            ENDDO
!
            IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
              DO ID=1, N_DIRECTION
                DO J=1, N_BRDF_BASIS_FNC
                  BRDF_SOL_LONG(LL, J, ID)=BRDF_SOL(LP, J, ID)
                ENDDO
              ENDDO
            ENDIF
!
          ELSE IF (ISOLIR == IP_INFRA_RED) THEN
!
            D_PLANCK_FLUX_SURFACE_LONG(LL)                              &
     &        =D_PLANCK_FLUX_SURFACE(LP)
            DO I=1, N_LAYER
              DIFF_PLANCK_LONG(LL, I)=DIFF_PLANCK(LP, I)
            ENDDO
            IF (L_IR_SOURCE_QUAD) THEN
              DO I=1, N_LAYER
                DIFF_PLANCK_2_LONG(LL, I)=DIFF_PLANCK_2(LP, I)
              ENDDO
            ENDIF
            IF (I_SPH_ALGORITHM == IP_SPH_REDUCED_ITER) THEN
              DO ID=1, N_DIRECTION
                DO J=1, N_BRDF_BASIS_FNC
                  BRDF_HEMI_LONG(LL, J, ID)=BRDF_HEMI(LP, J, ID)
                ENDDO
              ENDDO
            ENDIF
!
          ENDIF
!
!         Set the viewing directions.
          IF (I_SPH_MODE == IP_SPH_MODE_RAD) THEN
            DO ID=1, N_DIRECTION
              DIRECTION_LONG(LL, ID, 1)=DIRECTION(LP, ID, 1)
              DIRECTION_LONG(LL, ID, 2)=DIRECTION(LP, ID, 2)
            ENDDO
          ENDIF
!
          FLUX_INC_DOWN_LONG(LL)=FLUX_INC_DOWN(LP)
          DO JS=1, N_BRDF_BASIS_FNC
            RHO_ALB_LONG(LL, JS)=RHO_ALB(LP, JS)
          ENDDO
!
!         The curent notional column will contain the fraction of
!         the grid-box required for incrementing.
          WEIGHT_COLUMN(LL)=AREA_COLUMN(LP, ICL)
!
!
!         Prepare for the next column, moving on to the next grid-point
!         as required.
          ICS=ICS+1
          IF (ICS >  N_COLUMN_SLV(LP)) THEN
            LP=LP+1
            L_NEW=.TRUE.
          ENDIF
!
        ENDDO
!
!       Set N_LONG which will be required for the next batch after LL
!       has been reset.
        N_LONG=LL
!
!
! DEPENDS ON: sph_solver
        CALL SPH_SOLVER(IERR                                            &
!                       Atmospheric sizes
     &    , N_LONG, N_LAYER                                             &
!                       Angular integration
     &    , MS_MIN, MS_MAX, I_TRUNCATION, LS_LOCAL_TRUNC                &
     &    , CG_COEFF, UPLM_ZERO, IA_SPH_MM                              &
     &    , ACCURACY_ADAPTIVE, EULER_FACTOR                             &
     &    , I_SPH_ALGORITHM, I_SPH_MODE, L_RESCALE                      &
!                       Spectral Region
     &    , ISOLIR                                                      &
!                       Options for Equivalent Extinction
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE_LONG                         &
!                       Solar Fields
     &    , I_DIRECT_COLUMN, ZEN_0_LONG, UPLM_SOL_LONG                  &
!                       Infra-red Properties
     &    , DIFF_PLANCK_LONG, FLUX_INC_DOWN_LONG                        &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2_LONG                        &
!                       Optical properies
     &    , TAU_LONG, OMEGA_LONG, PHASE_FNC_LONG                        &
     &    , PHASE_FNC_SOLAR_LONG, FORWARD_SCATTER_LONG                  &
!                       Surface Conditions
     &    , LS_BRDF_TRUNC, N_BRDF_BASIS_FNC, RHO_ALB_LONG               &
     &    , F_BRDF, BRDF_SOL_LONG, BRDF_HEMI_LONG                       &
     &    , D_PLANCK_FLUX_SURFACE_LONG                                  &
!                       Levels for calculating radiances
     &    , N_VIEWING_LEVEL, I_RAD_LAYER, FRAC_RAD_LAYER                &
!                       Viewing Geometry
     &    , N_DIRECTION, DIRECTION_LONG                                 &
!                       Calculated Radiances or Fluxes
     &    , FLUX_DIRECT_COLUMN, FLUX_TOTAL_COLUMN, RADIANCE_COLUMN      &
     &    , PHOTOLYSIS_COLUMN                                           &
!                       Dimensions of arrays
     &    , ND_PROFILE_COLUMN, ND_LAYER                                 &
     &    , ND_PROFILE_COLUMN, ND_PROFILE_COLUMN, ND_PROFILE_COLUMN     &
     &    , ND_MAX_ORDER, ND_SPH_COEFF                                  &
     &    , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC                            &
     &    , ND_RED_EIGENSYSTEM, ND_SPH_EQUATION, ND_SPH_DIAGONAL        &
     &    , ND_SPH_CF_WEIGHT, ND_SPH_U_RANGE                            &
     &    , ND_VIEWING_LEVEL, ND_DIRECTION                              &
     &    )
!
!
!
!       Scatter the calculated fluxes or radiances back to their
!       appropriate grid-points.
        IF (I_SPH_MODE == IP_SPH_MODE_FLUX) THEN
!
          DO I=1, 2*N_LAYER+2
            DO LL=1, N_LONG
              L=TARGET(LL)
              FLUX_TOTAL(L, I)=FLUX_TOTAL(L, I)                         &
     &          +WEIGHT_COLUMN(LL)*FLUX_TOTAL_COLUMN(LL, I)
            ENDDO
          ENDDO
!
          IF (ISOLIR == IP_SOLAR) THEN
            DO I=0, N_LAYER
              DO LL=1, N_LONG
                L=TARGET(LL)
                FLUX_DIRECT(L, I)=FLUX_DIRECT(L, I)                     &
     &            +WEIGHT_COLUMN(LL)*FLUX_DIRECT_COLUMN(LL, I)
              ENDDO
            ENDDO
          ENDIF
!
        ELSE
!

          DO ID=1, N_DIRECTION
            DO I=1, N_VIEWING_LEVEL
              DO LL=1, N_LONG
                L=TARGET(LL)
                RADIANCE(L, I, ID)=RADIANCE(L, I, ID)                   &
     &            +WEIGHT_COLUMN(LL)*RADIANCE_COLUMN(LL, I, ID)

              ENDDO
            ENDDO
          ENDDO
!
          IF (ISOLIR == IP_SOLAR) THEN
            DO I=1, N_LAYER
              DO LL=1, N_LONG
                L=TARGET(LL)
                I_DIRECT(L, I)=I_DIRECT(L, I)                           &
     &            +WEIGHT_COLUMN(LL)*I_DIRECT_COLUMN(LL, I)
              ENDDO
            ENDDO
          ENDIF
!
        ENDIF
!
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE CALC_RADIANCE_IPA
#endif
