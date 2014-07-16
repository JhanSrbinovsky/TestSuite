#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve the two-stream equations in a mixed column.
!
! Method:
!       The two-stream coefficients are calculated in clear regions
!       and in stratiform and convective clouds. From these
!       coefficients transmission and reflection coefficients are
!       determined. The coefficients for convective and stratiform
!       clouds are appropriately mixed to form single cloudy values
!       and an appropriate solver is called.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-04-96                New solver added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Code for obsolete
!                                               solver removed.
!                                               (J. M. Edwards)
!       5.5             14-02-03                Obsolete header
!                                               files removed.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MIX_COLUMN(IERR                                        &
!                       Atmospheric Properties
     &   , N_PROFILE, N_LAYER                                           &
!                       Two-stream Scheme
     &   , I_2STREAM                                                    &
!                       Corrections to Two-stream Equations
     &   , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION           &
!                       Options for Solver
     &   , I_SOLVER, L_NET                                              &
!                       Options for Equivalent Extinction
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE                               &
!                       Spectral Region
     &   , ISOLIR                                                       &
!                       Infra-red Properties
     &   , DIFF_PLANCK                                                  &
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                              &
!                       Conditions at TOA
     &   , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                        &
!                       Conditions at Surface
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND       &
!                       Clear-sky Single Scattering Properties
     &   , TAU_FREE, OMEGA_FREE, ASYMMETRY_FREE                         &
!                       Cloud Geometry
     &   , N_CLOUD_TOP                                                  &
     &   , N_CLOUD_TYPE, FRAC_CLOUD                                     &
     &   , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE                       &
     &   , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE                    &
     &   , CLOUD_OVERLAP                                                &
!                       Cloudy Optical Properties
     &   , TAU_CLOUD, OMEGA_CLOUD, ASYMMETRY_CLOUD                      &
!                       Fluxes Calculated
     &   , FLUX_DIRECT, FLUX_TOTAL                                      &
!                       Flags for Clear-sky Calculations
     &   , L_CLEAR, I_SOLVER_CLEAR                                      &
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                          &
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
#include "stdio3a.h"
#include "dimfix3a.h"
#include "error3a.h"
#include "spcrg3a.h"
#include "solver3a.h"
#include "clcfpt3a.h"
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP                                                  &
!             TOP CLOUDY LAYER
     &   , N_CLOUD_TYPE                                                 &
!             NUMBER OF TYPES OF CLOUDS
     &   , N_FREE_PROFILE(NPD_LAYER)                                    &
!             NUMBER OF FREE PROFILES
     &   , I_FREE_PROFILE(NPD_PROFILE, NPD_LAYER)                       &
!             INDICES OF FREE PROFILES
     &   , N_CLOUD_PROFILE(NPD_LAYER)                                   &
!             NUMBER OF CLOUDY PROFILES
     &   , I_CLOUD_PROFILE(NPD_PROFILE, NPD_LAYER)                      &
!             INDICES OF CLOUDY PROFILES
     &   , ISOLIR                                                       &
!             SPECTRAL REGION
     &   , I_2STREAM                                                    &
!             TWO-STREAM SCHEME
     &   , I_SOLVER                                                     &
!             SOLVER USED
     &   , I_SOLVER_CLEAR
!             SOLVER FOR CLEAR-SKY FLUXES
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_NET                                                        &
!             CALCULATE NET FLUXES
     &   , L_CLEAR                                                      &
!             CALCULATE CLEAR-SKY FLUXES
     &   , L_SCALE_SOLAR                                                &
!             FLAG TO SCALE SOLAR
     &   , L_IR_SOURCE_QUAD                                             &
!             USE QUADRATIC SOURCE TERM
     &   , L_2_STREAM_CORRECT
!             EDGE CORRECTION TO 2-STREAM
!
!     OPTICAL PROPERTIES:
      REAL                                                              &
                !, INTENT(IN)
     &     TAU_FREE(NPD_PROFILE, NPD_LAYER)                             &
!             FREE OPTICAL DEPTH
     &   , OMEGA_FREE(NPD_PROFILE, NPD_LAYER)                           &
!             FREE ALBEDO OF SINGLE SCATTERING
     &   , ASYMMETRY_FREE(NPD_PROFILE, NPD_LAYER)                       &
!             CLEAR-SKY ASYMMETRY
     &   , TAU_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)            &
!             CLOUDY OPTICAL DEPTH
     &   , OMEGA_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)          &
!             CLOUDY ALBEDO OF SINGLE SCATTERING
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY ASYMMETRY
!
!     CLOUD GEOMETRY:
      REAL                                                              &
                !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             CLOUDY FRACTIONS IN EACH LAYER
     &   , W_FREE(NPD_PROFILE, NPD_LAYER)                               &
!             CLEAR SKY FRACTIONS IN EACH LAYER
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)           &
!             FRACTIONS OF DIFFERENT TYPES OF CLOUD
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)
!             ENERGY TRANSFER COEFFICIENTS
      REAL                                                              &
                !, INTENT(IN)
     &     SEC_0(NPD_PROFILE)                                           &
!             SECANT OF SOLAR ZENITH ANGLE
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)                             &
!             DIFFUSE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)                              &
!             DIRECT ALBEDO
     &   , FLUX_INC_DOWN(NPD_PROFILE)                                   &
!             INCIDENT TOTAL FLUX
     &   , FLUX_INC_DIRECT(NPD_PROFILE)                                 &
!             INCIDENT DIRECT FLUX
     &   , DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)                          &
!             CHANGE IN PLANCK FUNCTION
     &   , SOURCE_GROUND(NPD_PROFILE)                                   &
!             FLUX FROM SURFACE
     &   , ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)                      &
!             ADJUSTMENT OF SOLAR BEAM WITH EQUIVALENT EXTINCTION
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)                        &
!             2x2ND DIFFERENCE OF PLANCKIAN
     &   , PLANCK_SOURCE(NPD_PROFILE, 0: NPD_LAYER)                     &
!             PLANCKIAN SOURCE FUNCTION
     &   , GROUND_EMISSION(NPD_PROFILE)
!             TOTAL FLUX EMITTED FROM GROUND
!
!     FLUXES CALCULATED
      REAL                                                              &
                !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DIRECT FLUX
     &   , FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)                       &
!             LONG FLUX VECTOR
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                 &
!             CLEAR DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR TOTAL FLUX
!
!
!
!     LOCAL VARIABALES.
      INTEGER                                                           &
     &     N_SOURCE_COEFF                                               &
!             NUMBER OF SOURCE COEFFICIENTS
     &   , N_EQUATION                                                   &
!             NUMBER OF EQUATIONS
     &   , I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
!     CLEAR-SKY COEFFICIENTS:
      REAL                                                              &
     &     TRANS_FREE(NPD_PROFILE, NPD_LAYER)                           &
!             FREE TRANSMISSION OF LAYER
     &   , REFLECT_FREE(NPD_PROFILE, NPD_LAYER)                         &
!             FREE REFLECTANCE OF LAYER
     &   , TRANS_0_FREE(NPD_PROFILE, NPD_LAYER)                         &
!             FREE DIRECT TRANSMISSION OF LAYER
     &   , SOURCE_COEFF_FREE(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)  &
!             FREE SOURCE COEFFICIENTS
     &   , S_DOWN_FREE(NPD_PROFILE, NPD_LAYER)                          &
!             FREE DOWNWARD SOURCE
     &   , S_UP_FREE(NPD_PROFILE, NPD_LAYER)                            &
!             FREE UPWARD SOURCE
     &   , S_DOWN_CLEAR(NPD_PROFILE, NPD_LAYER)                         &
!             CLEAR DOWNWARD SOURCE
     &   , S_UP_CLEAR(NPD_PROFILE, NPD_LAYER)
!             CLEAR UPWARD SOURCE
!
!     CLOUDY COEFFICIENTS:
      REAL                                                              &
     &     TRANS_CLOUD(NPD_PROFILE, NPD_LAYER)                          &
!             CLOUDY TRANSMISSION OF LAYER
     &   , REFLECT_CLOUD(NPD_PROFILE, NPD_LAYER)                        &
!             CLOUDY REFLECTANCE OF LAYER
     &   , TRANS_0_CLOUD(NPD_PROFILE, NPD_LAYER)                        &
!             CLOUDY DIRECT TRANSMISSION OF LAYER
     &   , SOURCE_COEFF_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF) &
!             CLOUDY SOURCE COEFFICIENTS
     &   , S_DOWN_CLOUD(NPD_PROFILE, NPD_LAYER)                         &
!             CLOUDY DOWNWARD SOURCE
     &   , S_UP_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUDY UPWARD SOURCE
!
!     SOURCE FUNCTIONS AT THE CROUND
      REAL                                                              &
     &     SOURCE_GROUND_FREE(NPD_PROFILE)                              &
!             SOURCE FROM GROUND UNDER CLEAR SKIES
     &   , SOURCE_GROUND_CLOUD(NPD_PROFILE)                             &
!             SOURCE FROM GROUND UNDER CLOUDY SKIES
     &   , FLUX_DIRECT_GROUND_CLOUD(NPD_PROFILE)
!             DIRECT FLUX AT GROUND UNDER CLOUDY SKIES
!
!
!     NUMERICAL ARRAYS:
      REAL                                                              &
     &     A5(NPD_PROFILE, 5, 2*NPD_LAYER+2)                            &
!             PENTADIAGONAL MATRIX
     &   , B(NPD_PROFILE, 2*NPD_LAYER+2)                                &
!             RHS OF MATRIX EQUATION
     &   , WORK(NPD_PROFILE)
!             WORKING ARRAY FOR SOLVER
!
!     FUNCTIONS CALLED:
      INTEGER                                                           &
     &     SET_N_SOURCE_COEFF
!             FUNCTION TO SET NUMBER OF SOURCE COEFFICIENTS
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     TWO_COEFF, TWO_COEFF_CLOUD, IR_SOURCE, MIXED_SOLAR_SOURCE    &
     &   , BAND_SOLVER, MIX_COLUMN_FULL, MIX_APP_SCAT                   &
     &   , CLEAR_SUPPLEMENT
!
!
!
!     CALCULATE THE TRANSMISSION AND REFLECTION COEFFICIENTS AND
!     SOURCE TERMS FOR THE CLEAR AND CLOUDY PARTS OF THE COLUMN
!
!     SET THE NUMBER OF SOURCE COEFFICIENTS FOR THE APPROXIMATION
! DEPENDS ON: set_n_source_coeff
      N_SOURCE_COEFF=SET_N_SOURCE_COEFF(ISOLIR, L_IR_SOURCE_QUAD)
!
! DEPENDS ON: two_coeff
      CALL TWO_COEFF(IERR                                               &
     &   , N_PROFILE, 1, N_LAYER                                        &
     &   , I_2STREAM, L_IR_SOURCE_QUAD                                  &
     &   , ASYMMETRY_FREE, OMEGA_FREE, TAU_FREE                         &
     &   , ISOLIR, SEC_0                                                &
     &   , TRANS_FREE, REFLECT_FREE, TRANS_0_FREE                       &
     &   , SOURCE_COEFF_FREE                                            &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
      IF (IERR /= I_NORMAL) RETURN
!
!
!     INFRA-RED SOURCE TERMS DEPEND ONLY ON THE LAYER AND MAY BE
!     CALCULATED NOW. SOLAR TERMS DEPEND ON CONDITIONS IN CLOUD
!     IN OVERLYING LAYERS AND MUST BE CALCULATED LATER.
!
      IF (ISOLIR == IP_INFRA_RED) THEN
!
! DEPENDS ON: ir_source
         CALL IR_SOURCE(N_PROFILE, 1, N_LAYER                           &
     &      , SOURCE_COEFF_FREE, DIFF_PLANCK                            &
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                           &
     &      , L_2_STREAM_CORRECT, PLANCK_SOURCE                         &
     &      , GROUND_EMISSION, N_LAYER                                  &
     &      , TAU_FREE, TRANS_FREE                                      &
     &      , S_DOWN_FREE, S_UP_FREE                                    &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
!
!        IF A CLEAR-SKY CALCULATION IS REQUIRED THESE SOURCE TERMS MUST
!        BE STORED.
         IF (L_CLEAR) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  S_DOWN_CLEAR(L, I)=S_DOWN_FREE(L, I)
                  S_UP_CLEAR(L, I)=S_UP_FREE(L, I)
               ENDDO
            ENDDO
         ENDIF
!
!        SCALE THE SOURCES BY THE CLEAR-SKY FRACTIONS IN THE CLOUDY
!        LAYERS. IN HIGHER LAYERS THE CLEAR-SKY FRACTION IS 1.
         DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
               S_DOWN_FREE(L, I)=W_FREE(L, I)*S_DOWN_FREE(L, I)
               S_UP_FREE(L, I)=W_FREE(L, I)*S_UP_FREE(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
!     REPEAT THE CALCULATION FOR CLOUDY REGIONS.
!
!
! DEPENDS ON: two_coeff_cloud
      CALL TWO_COEFF_CLOUD(IERR                                         &
     &   , N_PROFILE, N_CLOUD_TOP, N_LAYER                              &
     &   , I_2STREAM, L_IR_SOURCE_QUAD, N_SOURCE_COEFF                  &
     &   , N_CLOUD_TYPE, FRAC_CLOUD                                     &
     &   , ASYMMETRY_CLOUD, OMEGA_CLOUD, TAU_CLOUD                      &
     &   , ISOLIR, SEC_0                                                &
     &   , TRANS_CLOUD, REFLECT_CLOUD, TRANS_0_CLOUD                    &
     &   , SOURCE_COEFF_CLOUD                                           &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
      IF (IERR /= I_NORMAL) RETURN
!
!
      IF (ISOLIR == IP_INFRA_RED) THEN
!
!        EDGE CORRECTIONS FOR THE TWO-STREAM EQUATIONS DO NOT
!        REALLY FIT WITH THIS METHOD OF TREATING CLOUDS. OPTICAL
!        DEPTHS AND TRANSMISSIONS MUST BE PASSED TO THE SUBROUTINE
!        TO FILL THE ARGUMENT LIST, BUT IT IS NOT INTENDED THAT
!        THESE ARRAYS WILL BE USED.
!
! DEPENDS ON: ir_source
         CALL IR_SOURCE(N_PROFILE, N_CLOUD_TOP, N_LAYER                 &
     &      , SOURCE_COEFF_CLOUD, DIFF_PLANCK                           &
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                           &
     &      , L_2_STREAM_CORRECT, PLANCK_SOURCE                         &
     &      , GROUND_EMISSION, N_LAYER                                  &
     &      , TAU_CLOUD, TRANS_CLOUD                                    &
     &      , S_DOWN_CLOUD, S_UP_CLOUD                                  &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
!
!
         DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
               S_DOWN_CLOUD(L, I)=W_CLOUD(L, I)*S_DOWN_CLOUD(L, I)
               S_UP_CLOUD(L, I)=W_CLOUD(L, I)*S_UP_CLOUD(L, I)
            ENDDO
         ENDDO
!
      ENDIF
!
!
!     CALCULATE THE APPROPRIATE SOURCE TERMS FOR THE SOLAR: CLOUDY
!     AND CLEAR PROPERTIES ARE BOTH NEEDED HERE.
!
      IF (ISOLIR == IP_SOLAR) THEN
!
! DEPENDS ON: mixed_solar_source
         CALL MIXED_SOLAR_SOURCE(N_PROFILE, N_LAYER, N_CLOUD_TOP        &
     &      , FLUX_INC_DIRECT                                           &
     &      , L_SCALE_SOLAR, ADJUST_SOLAR_KE                            &
     &      , TRANS_0_FREE, SOURCE_COEFF_FREE                           &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFC)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCC)                        &
     &      , TRANS_0_CLOUD, SOURCE_COEFF_CLOUD                         &
     &      , FLUX_DIRECT                                               &
     &      , FLUX_DIRECT_GROUND_CLOUD                                  &
     &      , S_UP_FREE, S_DOWN_FREE                                    &
     &      , S_UP_CLOUD, S_DOWN_CLOUD                                  &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &   )
      ENDIF
!
!
!
!     FORMULATE THE MATRIX EQUATIONS FOR THE FLUXES.
!
      SELECT CASE (i_solver)

      CASE (IP_solver_mix_app_scat)
!
! DEPENDS ON: mix_app_scat
         CALL MIX_APP_SCAT(N_PROFILE, N_LAYER, N_CLOUD_TOP              &
     &      , TRANS_FREE, REFLECT_FREE, S_DOWN_FREE, S_UP_FREE          &
     &      , TRANS_CLOUD, REFLECT_CLOUD                                &
     &      , S_DOWN_CLOUD, S_UP_CLOUD                                  &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFC)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCC)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BFF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BFC)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BCF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BCC)                        &
     &      , L_NET                                                     &
     &      , FLUX_INC_DOWN                                             &
     &      , SOURCE_GROUND, ALBEDO_SURFACE_DIFF                        &
     &      , FLUX_TOTAL                                                &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
!
      CASE (IP_solver_mix_direct, IP_solver_mix_direct_hogan)
!
!        SET THE PARTITIONED SOURCE FUNCTIONS AT THE GROUND.
         IF (ISOLIR == IP_SOLAR) THEN
            DO L=1, N_PROFILE
               SOURCE_GROUND_FREE(L)=(ALBEDO_SURFACE_DIR(L)             &
     &            -ALBEDO_SURFACE_DIFF(L))                              &
     &            *(FLUX_DIRECT(L, N_LAYER)                             &
     &            -FLUX_DIRECT_GROUND_CLOUD(L))
               SOURCE_GROUND_CLOUD(L)=(ALBEDO_SURFACE_DIR(L)            &
     &            -ALBEDO_SURFACE_DIFF(L))                              &
     &            *FLUX_DIRECT_GROUND_CLOUD(L)
            ENDDO
         ELSE
            DO L=1, N_PROFILE
               SOURCE_GROUND_FREE(L)                                    &
     &            =CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_BFF)             &
     &            *SOURCE_GROUND(L)
               SOURCE_GROUND_CLOUD(L)                                   &
     &            =CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_BCF)             &
     &            *SOURCE_GROUND(L)
            ENDDO
         ENDIF
!
         IF (i_solver == IP_solver_mix_direct) THEN
! DEPENDS ON: solver_mix_direct
           CALL SOLVER_MIX_DIRECT(N_PROFILE, N_LAYER, N_CLOUD_TOP       &
     &      , TRANS_FREE, REFLECT_FREE, S_DOWN_FREE, S_UP_FREE          &
     &      , TRANS_CLOUD, REFLECT_CLOUD                                &
     &      , S_DOWN_CLOUD, S_UP_CLOUD                                  &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFC)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCC)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BFF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BFC)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BCF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BCC)                        &
     &      , L_NET                                                     &
     &      , FLUX_INC_DOWN                                             &
     &      , SOURCE_GROUND_FREE, SOURCE_GROUND_CLOUD                   &
     &      , ALBEDO_SURFACE_DIFF                                       &
     &      , FLUX_TOTAL                                                &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )

         ELSE IF (i_solver == IP_solver_mix_direct_hogan) THEN
! DEPENDS ON: solver_mix_direct_hogan
           CALL solver_mix_direct_hogan(N_PROFILE, N_LAYER, N_CLOUD_TOP &
     &      , TRANS_FREE, REFLECT_FREE, S_DOWN_FREE, S_UP_FREE          &
     &      , TRANS_CLOUD, REFLECT_CLOUD                                &
     &      , S_DOWN_CLOUD, S_UP_CLOUD                                  &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GFC)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GCC)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BFF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BFC)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BCF)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BCC)                        &
     &      , L_NET                                                     &
     &      , FLUX_INC_DOWN                                             &
     &      , SOURCE_GROUND_FREE, SOURCE_GROUND_CLOUD                   &
     &      , ALBEDO_SURFACE_DIFF                                       &
     &      , FLUX_TOTAL                                                &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )

         ENDIF


      CASE (IP_solver_mix_11)
!
! DEPENDS ON: mix_column_full
         CALL MIX_COLUMN_FULL(N_PROFILE, N_LAYER, N_CLOUD_TOP           &
     &      , TRANS_FREE, REFLECT_FREE, S_DOWN_FREE, S_UP_FREE          &
     &      , TRANS_CLOUD, REFLECT_CLOUD                                &
     &      , S_DOWN_CLOUD, S_UP_CLOUD                                  &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GM)                         &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_GP)                         &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BM)                         &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_BP)                         &
     &      , FLUX_INC_DOWN                                             &
     &      , SOURCE_GROUND, ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR    &
     &      , FLUX_DIRECT(1, N_LAYER)                                   &
     &      , FLUX_TOTAL                                                &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
!
      CASE DEFAULT
!
         WRITE(IU_ERR, '(/A)')                                          &
     &      '***ERROR: THE SOLVER SPECIFIED IS NOT VALID HERE.'
         IERR=I_ERR_FATAL
         RETURN
!
      END SELECT
!
!
!
      IF (L_CLEAR) THEN
!
! DEPENDS ON: clear_supplement
         CALL CLEAR_SUPPLEMENT(IERR, N_PROFILE, N_LAYER, I_SOLVER_CLEAR &
     &      , TRANS_FREE, REFLECT_FREE, TRANS_0_FREE, SOURCE_COEFF_FREE &
     &      , ISOLIR, FLUX_INC_DIRECT, FLUX_INC_DOWN                    &
     &      , S_DOWN_CLEAR, S_UP_CLEAR                                  &
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR                   &
     &      , SOURCE_GROUND                                             &
     &      , L_SCALE_SOLAR, ADJUST_SOLAR_KE                            &
     &      , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                       &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE MIX_COLUMN
#endif
#endif
