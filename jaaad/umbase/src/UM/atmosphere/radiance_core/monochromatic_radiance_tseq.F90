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
      SUBROUTINE MONOCHROMATIC_RADIANCE_TSEQ(IERR                       &
!                     Atmospheric Propertries
     &  , N_PROFILE, N_LAYER                                            &
!                     Options for Solver
     &  , I_2STREAM, I_SOLVER, I_SCATTER_METHOD                         &
!                     Optical Properties
     &  , L_SCALE_SOLAR, ADJUST_SOLAR_KE                                &
!                     Spectral Region
     &  , ISOLIR                                                        &
!                     Infra-red Properties
     &  , DIFF_PLANCK, L_IR_SOURCE_QUAD, DIFF_PLANCK_2                  &
!                     Conditions at TOA
     &  , SEC_0, FLUX_INC_DIRECT, FLUX_INC_DOWN                         &
!                     Surface Properties
     &  , D_PLANCK_FLUX_SURFACE                                         &
     &  , RHO_ALB                                                       &
!                     Optical Properties
     &  , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR                             &
     &  , TAU, OMEGA, PHASE_FNC                                         &
!                     Cloudy Properties
     &  , I_CLOUD                                                       &
!                     Cloud Geometry
     &  , N_CLOUD_TOP                                                   &
     &  , N_CLOUD_TYPE, FRAC_CLOUD                                      &
     &  , N_REGION, K_CLR, I_REGION_CLOUD, FRAC_REGION                  &
     &  , W_FREE, W_CLOUD, CLOUD_OVERLAP                                &
     &  , N_COLUMN_SLV, LIST_COLUMN_SLV                                 &
     &  , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN                     &
!                     Fluxes Calculated
     &  , FLUX_DIRECT, FLUX_TOTAL                                       &
!                     Flags for Clear-sky Calculation
     &  , L_CLEAR, I_SOLVER_CLEAR                                       &
!                     Clear-sky Fluxes Calculated
     &  , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                           &
!                     Dimensions of Arrays
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN          &
     &  , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF                    &
     &  , ND_SOURCE_COEFF, ND_MAX_ORDER                                 &
     &  )
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Maximum number of profiles
     &  , ND_LAYER                                                      &
!           Maximum number of layers
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for completely clear layers
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_COLUMN                                                     &
!           Number of columns per point
     &  , ND_CLOUD_TYPE                                                 &
!           Maximum number of types of cloud
     &  , ND_REGION                                                     &
!           Maximum number of cloudy regions
     &  , ND_OVERLAP_COEFF                                              &
!           Maximum number of overlap coefficients
     &  , ND_SOURCE_COEFF                                               &
!           Size allocated for source coefficients
     &  , ND_MAX_ORDER
!           Size allocated for spherical harmonics
!
!     Include header files.
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
#include "cloud_scheme_pcf3z.h"
#include "spectral_region_pcf3z.h"
#include "solver_pcf3z.h"
#include "surface_spec_pcf3z.h"
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
!
!                     Angular integration
      INTEGER, INTENT(IN) ::                                            &
     &    I_2STREAM                                                     &
!           Two-stream scheme
     &  , I_SCATTER_METHOD
!           Method of treating scattering
!
!                     Options for solver
      INTEGER, INTENT(IN) ::                                            &
     &    I_SOLVER
!           Solver used
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
     &    SEC_0(ND_PROFILE)                                             &
!           Secants (two-stream) or cosines (spherical harmonics)
!           of the solar zenith angles
     &  , FLUX_INC_DIRECT(ND_PROFILE)                                   &
!           Incident direct flux
     &  , FLUX_INC_DOWN(ND_PROFILE)
!           Incident downward flux
!
!                     Surface properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    D_PLANCK_FLUX_SURFACE(ND_PROFILE)
!           Differential Planckian flux from the surface
      REAL  (Real64), INTENT(IN) ::                                     &
     &    RHO_ALB(ND_PROFILE, 2)
!           Weights of the basis functions
!
!                     Optical properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU_CLR(ND_PROFILE, ND_LAYER_CLR)                             &
!           Clear-sky optical depth
     &  , OMEGA_CLR(ND_PROFILE, ND_LAYER_CLR)                           &
!           Clear-sky albedo of single scattering
     &  , PHASE_FNC_CLR(ND_PROFILE, ND_LAYER_CLR, ND_MAX_ORDER)
!           Clear-sky phase function
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)            &
!           Optical depth
     &  , OMEGA(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)          &
!           Albedo of single scattering
     &  , PHASE_FNC(ND_PROFILE, ID_CT: ND_LAYER                         &
     &      , ND_MAX_ORDER, 0: ND_CLOUD_TYPE)
!           Phase functions
!
!                     Cloudy properties
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
     &  , I_REGION_CLOUD(ND_CLOUD_TYPE)                                 &
!           Regions in which types of clouds fall
     &  , K_CLR
!           Index of clear-sky region
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
!                     Calculated Fluxes
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_PROFILE, 0: ND_LAYER)                          &
!           Direct flux
     &  , FLUX_TOTAL(ND_PROFILE, 2*ND_LAYER+2)
!           Total flux
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
     &    FLUX_DIRECT_CLEAR(ND_PROFILE, 0: ND_LAYER)                    &
!           Clear-sky direct flux
     &  , FLUX_TOTAL_CLEAR(ND_PROFILE, 2*ND_LAYER+2)
!           Clear-sky total flux
!
!
!
!     Local variables.
      INTEGER                                                           &
     &    ND_PROFILE_COLUMN
!           Size allowed for profiles taken simultaneously in a
!           decomposition into columns
      INTEGER                                                           &
     &    L                                                             &
!           Loop variable
     &  , I
!           Loop variable
!     Full clear-sky single-scattering properties
      REAL  (Real64), ALLOCATABLE ::                                    &
     &    TAU_CLR_F(:, :)                                               &
!           Clear-sky optical depth
     &  , OMEGA_CLR_F(:, :)                                             &
!           Clear-sky albedo of single scattering
     &  , PHASE_FNC_CLR_F(:, :, :)
!           Clear-sky phase function
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    TWO_STREAM, MIX_COLUMN, TRIPLE_COLUMN, CALC_FLUX_IPA          &
     &  , COPY_CLR_FULL
!
!
!
!     Choose an appropriate routine to calculate the fluxes as
!     determined by the cloud scheme selected.
!
      IF (I_CLOUD == IP_CLOUD_CLEAR) THEN
!       Allocate and set dynamic arrays.
        ALLOCATE(TAU_CLR_F(ND_PROFILE, ND_LAYER))
        ALLOCATE(OMEGA_CLR_F(ND_PROFILE, ND_LAYER))
        ALLOCATE(PHASE_FNC_CLR_F(ND_PROFILE, ND_LAYER, 1))
!
! DEPENDS ON: copy_clr_full
        CALL COPY_CLR_FULL(N_PROFILE, N_LAYER, N_CLOUD_TOP              &
     &    , 1                                                           &
     &    , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR(1, 1, 1)                  &
     &    , TAU, OMEGA, PHASE_FNC(1, 1, 1, 0)                           &
     &    , TAU_CLR_F, OMEGA_CLR_F, PHASE_FNC_CLR_F(1, 1, 1)            &
!                       Sizes of arrays
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, 1                &
     &    )
!
!       A two-stream scheme with no clouds.
! DEPENDS ON: two_stream
        CALL TWO_STREAM(IERR                                            &
!                     Atmospheric properties
     &    , N_PROFILE, N_LAYER                                          &
!                     Two-stream scheme
     &    , I_2STREAM                                                   &
!                     Options for solver
     &    , I_SOLVER                                                    &
!                     Options for equivalent extinction
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
!                     Spectral region
     &    , ISOLIR                                                      &
!                     Infra-red properties
     &    , DIFF_PLANCK                                                 &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                             &
!                     Conditions at TOA
     &    , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                       &
!                     Surface conditions
     &    , RHO_ALB(1, IP_SURF_ALB_DIFF)                                &
     &    , RHO_ALB(1, IP_SURF_ALB_DIR), D_PLANCK_FLUX_SURFACE          &
!                     Single scattering properties
     &    , TAU_CLR_F, OMEGA_CLR_F, PHASE_FNC_CLR_F(1, 1, 1)            &
!                     Fluxes calculated
     &    , FLUX_DIRECT, FLUX_TOTAL                                     &
!                     Sizes of arrays
     &    , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                       &
     &    )
!
!       Release temporary storage.
        DEALLOCATE(TAU_CLR_F)
        DEALLOCATE(OMEGA_CLR_F)
        DEALLOCATE(PHASE_FNC_CLR_F)
!
        IF (L_CLEAR) THEN
!         The clear fluxes here can be copied directly without
!         any further calculation.
          IF (ISOLIR == IP_SOLAR) THEN
            DO I=0, N_LAYER
              DO L=1, N_PROFILE
                FLUX_DIRECT_CLEAR(L, I)=FLUX_DIRECT(L, I)
              ENDDO
            ENDDO
          ENDIF
          DO I=1, 2*N_LAYER+2
            DO L=1, N_PROFILE
              FLUX_TOTAL_CLEAR(L, I)=FLUX_TOTAL(L, I)
            ENDDO
          ENDDO
        ENDIF
!
      ELSEIF ( (I_CLOUD == IP_CLOUD_MIX_MAX).OR.                        &
     &         (I_CLOUD == IP_CLOUD_MIX_RANDOM).OR.                     &
     &         ( (I_CLOUD == IP_CLOUD_PART_CORR).AND.                   &
     &           (N_REGION == 2) ) ) THEN
!
!       Clouds are treated using the coupled overlaps originally
!       introduced by Geleyn and Hollingsworth.
!
! DEPENDS ON: mix_column
        CALL MIX_COLUMN(IERR                                            &
!                     Atmospheric properties
     &    , N_PROFILE, N_LAYER, K_CLR                                   &
!                     Two-stream scheme
     &    , I_2STREAM                                                   &
!                     Options for solver
     &    , I_SOLVER                                                    &
!                     Options for equivalent extinction
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
!                     Spectral region
     &    , ISOLIR                                                      &
!                     Infra-red properties
     &    , DIFF_PLANCK                                                 &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                             &
!                     Conditions at TOA
     &    , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                       &
!                     Conditions at surface
     &    , RHO_ALB(1, IP_SURF_ALB_DIFF), RHO_ALB(1, IP_SURF_ALB_DIR)   &
     &    , D_PLANCK_FLUX_SURFACE                                       &
!                     Single scattering properties
     &    , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR(1, 1, 1)                  &
     &    , TAU, OMEGA, PHASE_FNC                                       &
!                     Cloud geometry
     &    , N_CLOUD_TOP                                                 &
     &    , N_CLOUD_TYPE, FRAC_CLOUD                                    &
     &    , W_FREE, W_CLOUD                                             &
     &    , CLOUD_OVERLAP                                               &
!                     Fluxes calculated
     &    , FLUX_DIRECT, FLUX_TOTAL                                     &
!                     Flags for clear-sky calculations
     &    , L_CLEAR, I_SOLVER_CLEAR                                     &
!                     Clear-sky fluxes calculated
     &    , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                         &
!                     Dimensions of arrays
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT                   &
     &    , ND_MAX_ORDER, ND_SOURCE_COEFF                               &
     &    , ND_CLOUD_TYPE, ND_OVERLAP_COEFF                             &
     &    )
        IF (IERR /= I_NORMAL) RETURN
!
      ELSEIF ( (I_CLOUD == IP_CLOUD_TRIPLE).OR.                         &
     &         ( (I_CLOUD == IP_CLOUD_PART_CORR_CNV).AND.               &
     &           (N_REGION == 3) ) ) THEN
!
!       Clouds are treated using a decomposition of the column
!       into clear-sky, stratiform and convective regions.
!
! DEPENDS ON: triple_column
        CALL TRIPLE_COLUMN(IERR                                         &
!                     Atmospheric properties
     &    , N_PROFILE, N_LAYER                                          &
!                     Two-stream scheme
     &    , I_2STREAM                                                   &
!                     Options for solver
     &    , I_SOLVER, I_SCATTER_METHOD                                  &
!                     Options for equivalent extinction
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
!                     Spectral region
     &    , ISOLIR                                                      &
!                     Infra-red properties
     &    , DIFF_PLANCK                                                 &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                             &
!                     Conditions at TOA
     &    , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                       &
!                     Conditions at surface
     &    , RHO_ALB(1, IP_SURF_ALB_DIFF), RHO_ALB(1, IP_SURF_ALB_DIR)   &
     &    , D_PLANCK_FLUX_SURFACE                                       &
!                     Single scattering properties
     &    , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR                           &
     &    , TAU, OMEGA, PHASE_FNC                                       &
!                     Cloud geometry
     &    , N_CLOUD_TOP                                                 &
     &    , N_CLOUD_TYPE, FRAC_CLOUD                                    &
     &    , N_REGION, I_REGION_CLOUD, FRAC_REGION                       &
     &    , W_FREE, W_CLOUD                                             &
     &    , CLOUD_OVERLAP                                               &
!                     Fluxes calculated
     &    , FLUX_DIRECT, FLUX_TOTAL                                     &
!                     Flags for clear-sky calculations
     &    , L_CLEAR, I_SOLVER_CLEAR                                     &
!                     Clear-sky fluxes calculated
     &    , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                         &
!                     Dimensions of arrays
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT                   &
     &    , ND_MAX_ORDER, ND_SOURCE_COEFF                               &
     &    , ND_CLOUD_TYPE, ND_REGION, ND_OVERLAP_COEFF                  &
     &    )
        IF (IERR /= I_NORMAL) RETURN
!
#if ! defined(EXCLOP)
      ELSEIF (I_CLOUD == IP_CLOUD_COLUMN_MAX) THEN
!       Clouds are treated on the assumption of maximum overlap
!       in a column model.
#else
      ELSEIF ( (I_CLOUD == IP_CLOUD_COLUMN_MAX).OR.                     &
     &         (I_CLOUD == IP_CLOUD_OVLPTEST) ) THEN
!       Clouds are treated using an explicit decomposition
!       into columns.
#endif
!
!       Set a dimension to allow the subcolumns of several profiles
!       to be considered at once.
        ND_PROFILE_COLUMN=MAX(1, N_PROFILE)
        DO L=1, N_PROFILE
          ND_PROFILE_COLUMN=MAX(ND_PROFILE_COLUMN, N_COLUMN_SLV(L))
        ENDDO
!
! DEPENDS ON: calc_flux_ipa
        CALL CALC_FLUX_IPA(IERR                                         &
!                     Atmospheric properties
     &    , N_PROFILE, N_LAYER, N_CLOUD_TOP                             &
!                     Options for equivalent extinction
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
!                     Algorithmic options
     &    , I_2STREAM, I_SOLVER                                         &
!                     Spectral region
     &    , ISOLIR                                                      &
!                     Infra-red properties
     &    , DIFF_PLANCK                                                 &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                             &
!                     Conditions at TOA
     &    , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                       &
!                     Conditions at surface
     &    , D_PLANCK_FLUX_SURFACE, RHO_ALB                              &
!                     Single scattering properties
     &    , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR                           &
     &    , TAU, OMEGA, PHASE_FNC                                       &
!                     Cloud geometry
     &    , N_COLUMN_SLV, LIST_COLUMN_SLV                               &
     &    , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN                   &
!                     Fluxes calculated
     &    , FLUX_DIRECT, FLUX_TOTAL                                     &
!                     Flags for clear-sky calculations
     &    , L_CLEAR, I_SOLVER_CLEAR                                     &
!                     Clear-sky fluxes calculated
     &    , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                         &
!                     Dimensions of arrays
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN        &
     &    , ND_MAX_ORDER, ND_CLOUD_TYPE                                 &
     &    , ND_PROFILE_COLUMN, ND_SOURCE_COEFF                          &
     &    )
        IF (IERR /= I_NORMAL) RETURN
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE MONOCHROMATIC_RADIANCE_TSEQ
#endif
