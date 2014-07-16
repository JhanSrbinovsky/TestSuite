#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate monochromatic fluxes using IPA.
!
! Method:
!
!   In this subroutine a long vector for two-stream flux calculations
!   is set up using the information on the types of cloud present.
!
! Current owner of code: James Manners
!
! Description of code:
!   Fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_FLUX_IPA(IERR                                     &
!                     Atmospheric Properties
     &  , N_PROFILE, N_LAYER, N_CLOUD_TOP                               &
!                     Options for Equivalent Extinction
     &  , L_SCALE_SOLAR, ADJUST_SOLAR_KE                                &
!                     Algorithmic options
     &  , I_2STREAM, I_SOLVER                                           &
!                     Spectral Region
     &  , ISOLIR                                                        &
!                     Infra-red Properties
     &  , DIFF_PLANCK                                                   &
     &  , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                               &
!                     Conditions at TOA
     &  , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                         &
!                     Conditions at Surface
     &  , D_PLANCK_FLUX_SURFACE, RHO_ALB                                &
!                     Single Scattering Properties
     &  , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR                             &
     &  , TAU, OMEGA, PHASE_FNC                                         &
!                     Cloud Geometry
     &  , N_COLUMN_SLV, LIST_COLUMN_SLV                                 &
     &  , I_CLM_LYR_CHN, I_CLM_CLD_TYP, AREA_COLUMN                     &
!                       Calculated fluxes
     &  , FLUX_DIRECT, FLUX_TOTAL                                       &
!                     Options for clear-sky fluxes
     &  , L_CLEAR, I_SOLVER_CLEAR                                       &
!                       Calculated fluxes
     &  , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                           &
!                     Dimensions of Arrays
     &  , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, ND_COLUMN          &
     &  , ND_MAX_ORDER, ND_CLOUD_TYPE                                   &
     &  , ND_PROFILE_COLUMN, ND_SOURCE_COEFF                            &
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
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_LAYER_CLR                                                  &
!           Size allocated for totally clear atmospheric layers
     &  , ND_COLUMN                                                     &
!           Size allocated for columns at a grid-point
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spectral calculations
!           (Here used only to ensure that dimensions are correct)
     &  , ND_CLOUD_TYPE                                                 &
!           Size allocated for types of clouds
     &  , ND_PROFILE_COLUMN                                             &
!           Number of profiles of subcolumns considered at once
     &  , ID_CT                                                         &
!           Topmost declared cloudy layer
     &  , ND_SOURCE_COEFF
!           Number of coefficients in the source function
!
!     Include header files.
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
#include "spectral_region_pcf3z.h"
#include "surface_spec_pcf3z.h"
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
     &  , N_CLOUD_TOP
!           Topmost cloudy layer
!
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR                                                        &
!           Spectral region
     &  , I_2STREAM                                                     &
!           Two-stream scheme selected
     &  , I_SOLVER
!           Solver selected
      LOGICAL, INTENT(IN) ::                                            &
     &    L_SCALE_SOLAR                                                 &
!           Scale solar beam
     &  , L_IR_SOURCE_QUAD
!           Use a quadratic source term
!           the singly scattered solar beam
!
!     Fields for equivalent extinction
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ADJUST_SOLAR_KE(ND_PROFILE, ND_LAYER)
!           Adjustment of solar beam with equivalent extinction
!
!     Optical properties
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU_CLR(ND_PROFILE, ND_LAYER_CLR)                             &
!           Clear-sky optical depth
     &  , OMEGA_CLR(ND_PROFILE, ND_LAYER_CLR)                           &
!           Clear-sky albedos of single scattering
     &  , PHASE_FNC_CLR(ND_PROFILE, ND_LAYER_CLR, ND_MAX_ORDER)
!           Moments of the clear-sky phase function
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)            &
!           Optical depth
     &  , OMEGA(ND_PROFILE, ID_CT: ND_LAYER, 0: ND_CLOUD_TYPE)          &
!           Single scattering albedos
     &  , PHASE_FNC(ND_PROFILE, ID_CT: ND_LAYER                         &
     &      , ND_MAX_ORDER, 0: ND_CLOUD_TYPE)
!           Moments of the phase functions
!           (Note that the third dimension is fixed at 1: the code
!           should not be called under other circumstances.)
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
     &    SEC_0(ND_PROFILE)                                             &
!           Secant of zenith angle
     &  , FLUX_INC_DIRECT(ND_PROFILE)                                   &
!           Incident direct flux
     &  , FLUX_INC_DOWN(ND_PROFILE)
!           Incident total flux
!
!     Conditions at surface
      REAL  (Real64), INTENT(IN) ::                                     &
     &    RHO_ALB(ND_PROFILE, 2)
!           Weights of the basis functions
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
!                       Calculated Fluxes
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_PROFILE, 0: ND_LAYER)                          &
!           Direct Flux
     &  , FLUX_TOTAL(ND_PROFILE, 2*ND_LAYER+2)
!           Total Fluxes
!
!                     Options for clear-sky fluxes
      LOGICAL                                                           &
     &    L_CLEAR
!           Flag for clear-sky fluxes
      INTEGER                                                           &
     &    I_SOLVER_CLEAR
!           Solver selected for clear-sky fluxes
!                       Calculated clear-sky fluxes
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT_CLEAR(ND_PROFILE, 0: ND_LAYER)                    &
!           Direct Clear-sky Flux
     &  , FLUX_TOTAL_CLEAR(ND_PROFILE, 2*ND_LAYER+2)
!           Total Clear-sky Flux
!
!
!
!     Local variables.
      INTEGER                                                           &
     &    I                                                             &
!           Loop variable
     &  , L                                                             &
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
     &    WEIGHT_LONG(ND_PROFILE_COLUMN)
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
     &  , ASYMMETRY_LONG(ND_PROFILE_COLUMN, ND_LAYER)                   &
!           Long vector of asymmetries
     &  , ADJUST_SOLAR_KE_LONG(ND_PROFILE_COLUMN, ND_LAYER)             &
!           Long vector of solar scalings
     &  , SEC_0_LONG(ND_PROFILE_COLUMN)                                 &
!           Long vector of cosines of the solar zenith angle
     &  , DIFF_PLANCK_LONG(ND_PROFILE_COLUMN, ND_LAYER)                 &
!           Long vector of differences in the Planckian
     &  , DIFF_PLANCK_2_LONG(ND_PROFILE_COLUMN, ND_LAYER)               &
!           Long vector of second differences in the Planckian
     &  , FLUX_INC_DIRECT_LONG(ND_PROFILE_COLUMN)                       &
!           Long vector of incident direct downward fluxes
     &  , FLUX_INC_DOWN_LONG(ND_PROFILE_COLUMN)                         &
!           Long vector of incident downward fluxes
     &  , D_PLANCK_FLUX_SURFACE_LONG(ND_PROFILE_COLUMN)                 &
!           Long vector of differential Planckian fluxes
!           at the surface
     &  , RHO_ALB_LONG(ND_PROFILE_COLUMN, 2)
!           Long vector of weightings of BRDF basis functions
!
!     Calculated Fluxes in subcolumns
      REAL  (Real64) ::                                                 &
     &    FLUX_DIRECT_LONG(ND_PROFILE_COLUMN, 0: ND_LAYER)              &
!           Direct Flux
     &  , FLUX_TOTAL_LONG(ND_PROFILE_COLUMN, 2*ND_LAYER+2)
!           Total Fluxes
!
!     Clear-sky optical properties of the whole column
      REAL  (Real64), ALLOCATABLE ::                                    &
     &    TAU_CLR_F(:, :)                                               &
!           Clear-sky optical depth for the whole column
     &  , OMEGA_CLR_F(:, :)                                             &
!           Clear-sky albedos of single scattering for the whole column
     &  , PHASE_FNC_CLR_F(:, :, :)
!           Moments of the clear-sky phase function for the whole column
!
!
!
!
!     Functions called:
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    TWO_STREAM
!
!
!
!     Zero the output arrays ready for incrementing.
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
              ASYMMETRY_LONG(LL, I)=PHASE_FNC_CLR(LP, I, 1)
            ENDDO
            DO I=N_CLOUD_TOP, N_LAYER
              TAU_LONG(LL, I)=TAU(LP, I, 0)
              OMEGA_LONG(LL, I)=OMEGA(LP, I, 0)
              ASYMMETRY_LONG(LL, I)=PHASE_FNC(LP, I, 1, 0)
            ENDDO
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
              ASYMMETRY_LONG(LL, I)                                     &
     &          =ASYMMETRY_LONG(LL_COPY, I)
            ENDDO
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
            ASYMMETRY_LONG(LL, ICC)                                     &
     &        =PHASE_FNC(LP, ICC, 1, ICT)
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
            SEC_0_LONG(LL)=SEC_0(LP)
            FLUX_INC_DIRECT_LONG(LL)=FLUX_INC_DIRECT(LP)
            D_PLANCK_FLUX_SURFACE_LONG(LL)=0.0E+00_Real64
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
!
          ENDIF
!
          FLUX_INC_DOWN_LONG(LL)=FLUX_INC_DOWN(LP)
          RHO_ALB_LONG(LL, IP_SURF_ALB_DIR)                             &
     &      =RHO_ALB(LP, IP_SURF_ALB_DIR)
          RHO_ALB_LONG(LL, IP_SURF_ALB_DIFF)                            &
     &      =RHO_ALB(LP, IP_SURF_ALB_DIFF)
!
!         The curent notional column will contain the fraction of
!         the grid-box required for incrementing.
          WEIGHT_LONG(LL)=AREA_COLUMN(LP, ICL)
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
!       N.B. The clear-sky option cannot be used here.
! DEPENDS ON: two_stream
        CALL TWO_STREAM(IERR                                            &
!                       Atmospheric properties
     &    , N_LONG, N_LAYER                                             &
!                       Two-stream scheme
     &    , I_2STREAM                                                   &
!                       Options for solver
     &    , I_SOLVER                                                    &
!                       Options for equivalent extinction
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE_LONG                         &
!                       Spectral region
     &    , ISOLIR                                                      &
!                       Infra-red properties
     &    , DIFF_PLANCK_LONG                                            &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2_LONG                        &
!                       Conditions at TOA
     &    , FLUX_INC_DOWN_LONG, FLUX_INC_DIRECT_LONG, SEC_0_LONG        &
!                       Surface conditions
     &    , RHO_ALB_LONG(1, IP_SURF_ALB_DIFF)                           &
     &    , RHO_ALB_LONG(1, IP_SURF_ALB_DIR), D_PLANCK_FLUX_SURFACE_LONG&
!                       Single scattering properties
     &    , TAU_LONG, OMEGA_LONG, ASYMMETRY_LONG(1, 1)                  &
!                       Fluxes calculated
     &    , FLUX_DIRECT_LONG, FLUX_TOTAL_LONG                           &
!                       Sizes of arrays
     &    , ND_PROFILE_COLUMN, ND_LAYER, ND_SOURCE_COEFF                &
     &    )
!
!
!
!       Scatter the calculated fluxes back to their
!       appropriate grid-points.
!
        DO I=1, 2*N_LAYER+2
          DO LL=1, N_LONG
            L=TARGET(LL)
            FLUX_TOTAL(L, I)=FLUX_TOTAL(L, I)                           &
     &        +WEIGHT_LONG(LL)*FLUX_TOTAL_LONG(LL, I)
          ENDDO
        ENDDO
!
        IF (ISOLIR == IP_SOLAR) THEN
          DO I=0, N_LAYER
            DO LL=1, N_LONG
              L=TARGET(LL)
              FLUX_DIRECT(L, I)=FLUX_DIRECT(L, I)                       &
     &          +WEIGHT_LONG(LL)*FLUX_DIRECT_LONG(LL, I)
            ENDDO
          ENDDO
        ENDIF
!
!
      ENDDO
!
!     Calculate the clear-sky fluxes if required.
      IF (L_CLEAR) THEN
!
!       Set aside space for the clear optical properties and copy
!       them across.
        ALLOCATE(TAU_CLR_F(ND_PROFILE, ND_LAYER))
        ALLOCATE(OMEGA_CLR_F(ND_PROFILE, ND_LAYER))
        ALLOCATE(PHASE_FNC_CLR_F(ND_PROFILE, ND_LAYER, 1))
!
! DEPENDS ON: copy_clr_full
        CALL COPY_CLR_FULL(N_PROFILE, N_LAYER, N_CLOUD_TOP, 1           &
     &    , TAU_CLR, OMEGA_CLR, PHASE_FNC_CLR                           &
     &    , TAU, OMEGA, PHASE_FNC                                       &
     &    , TAU_CLR_F, OMEGA_CLR_F, PHASE_FNC_CLR_F                     &
!                       Sizes of arrays
     &    , ND_PROFILE, ND_LAYER, ND_LAYER_CLR, ID_CT, 1                &
     &    )
!
! DEPENDS ON: two_stream
        CALL TWO_STREAM(IERR                                            &
!                       Atmospheric properties
     &    , N_PROFILE, N_LAYER                                          &
!                       Two-stream scheme
     &    , I_2STREAM                                                   &
!                       Options for solver
     &    , I_SOLVER_CLEAR                                              &
!                       Options for equivalent extinction
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
!                       Spectral region
     &    , ISOLIR                                                      &
!                       Infra-red properties
     &    , DIFF_PLANCK                                                 &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                             &
!                       Conditions at TOA
     &    , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                       &
!                       Surface conditions
     &    , RHO_ALB(1, IP_SURF_ALB_DIFF)                                &
     &    , RHO_ALB(1, IP_SURF_ALB_DIR), D_PLANCK_FLUX_SURFACE          &
!                       Single scattering properties
     &    , TAU_CLR_F, OMEGA_CLR_F, PHASE_FNC_CLR_F(1, 1, 1)            &
!                       Fluxes calculated
     &    , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                         &
!                       Sizes of arrays
     &    , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                       &
     &    )
!
!       Remove the arrays that are no longer required.
        DEALLOCATE(TAU_CLR_F)
        DEALLOCATE(OMEGA_CLR_F)
        DEALLOCATE(PHASE_FNC_CLR_F)
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE CALC_FLUX_IPA
#endif
