


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve for the monochromatic fluxes.
!
! Method:
!       The final single scattering properties are calculated
!       and rescaled. An appropriate subroutine is called to
!       calculate the fluxes depending on the treatment of
!       cloudiness.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             08-08-96                Code for vertically
!                                               coherent cloud added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Variable for obsolete
!                                               solver removed.
!                                               Unused variables
!                                               removed from call
!                                               to TRPILE_COLUMN.
!                                               (J. M. Edwards)
!       5.3             04-10-01                Number of cloudy regions
!                                               passed explicitly. The
!                                               treatment of scattering
!                                               is passed to
!                                               TRIPLE_COLUMN.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE MONOCHROMATIC_FLUX(IERR                                &
!                       Atmospheric Propetries
     &   , N_PROFILE, N_LAYER, D_MASS                                   &
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT         &
     &   , L_RESCALE, N_ORDER_GAUSS                                     &
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_BAND                                        &
!                       Options for Solver
     &   , I_SOLVER, L_NET, N_AUGMENT                                   &
!                       Gaseous Propeties
     &   , K_GAS_ABS                                                    &
!                       Options for Equivalent Extinction
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE                               &
!                       Spectral Region
     &   , ISOLIR                                                       &
!                       Infra-red Properties
     &   , DIFF_PLANCK                                                  &
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                              &
!                       Conditions at TOA
     &   , SEC_0, FLUX_INC_DIRECT, FLUX_INC_DOWN                        &
!                       Surface Propeties
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND       &
     &   , GROUND_EMISSION                                              &
!                       Clear-sky Optical Propeties
     &   , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE             &
     &   , FORWARD_SCATTER_FREE                                         &
!                       Cloudy Properties
     &   , L_CLOUD, I_CLOUD                                             &
!                       Cloud Geometry
     &   , N_CLOUD_TOP                                                  &
     &   , N_CLOUD_TYPE, FRAC_CLOUD                                     &
     &   , N_REGION, I_REGION_CLOUD, FRAC_REGION                        &
     &   , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE                       &
     &   , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE                    &
     &   , CLOUD_OVERLAP                                                &
     &   , N_COLUMN, L_COLUMN, AREA_COLUMN                              &
!                       Cloudy Optical Propeties
     &   , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                           &
     &   , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD                       &
!                       Fluxes Calculated
     &   , FLUX_DIRECT, FLUX_TOTAL                                      &
!                       Flags for Clear-sky Calculation
     &   , L_CLEAR, I_SOLVER_CLEAR                                      &
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                          &
!                       Planckian Source
     &   , PLANCK_SOURCE                                                &
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                           &
     &   )
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
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_COLUMN
!             NUMBER OF COLUMNS PER POINT
!
!     INCLUDE COMDECKS.
! DIMFIX3A defines internal dimensions tied to algorithms for
! two-stream radiation code, mostly for clouds

      ! number of components of clouds
      INTEGER,PARAMETER:: NPD_CLOUD_COMPONENT=4

      ! number of permitted types of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_TYPE=4

      ! number of permitted representations of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_REPRESENTATION=4

      ! number of overlap coefficients for clouds
      INTEGER,PARAMETER:: NPD_OVERLAP_COEFF=18

      ! number of coefficients for two-stream sources
      INTEGER,PARAMETER:: NPD_SOURCE_COEFF=2

      ! number of regions in a layer
      INTEGER,PARAMETER:: NPD_REGION=3

! DIMFIX3A end
! CLSCHM3A defines reference numbers for cloud schemes in two-stream
! radiation code.

      ! maximum/random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_MAX=2

      ! random overlap in a mixed column
      INTEGER,PARAMETER:: IP_CLOUD_MIX_RANDOM=4

      ! maximum overlap in a column model
      INTEGER,PARAMETER:: IP_CLOUD_COLUMN_MAX=3

      ! clear column
      INTEGER,PARAMETER:: IP_CLOUD_CLEAR=5

      ! mixed column with split between  convective and layer cloud.
      INTEGER,PARAMETER:: IP_CLOUD_TRIPLE=6

      ! Coupled overlap with partial correlation of cloud
      INTEGER,Parameter:: IP_cloud_part_corr=7

      ! Coupled overlap with partial correlation of cloud
      ! with a separate treatment of convective cloud
      INTEGER,Parameter:: IP_cloud_part_corr_cnv=8

! CLSCHM3A end
! ANGINT3A defines types of angular integration for two-stream
! radiation code.

      INTEGER,PARAMETER:: IP_TWO_STREAM=1 ! two stream scheme

      ! gaussian integration in the IR
      INTEGER,PARAMETER:: IP_IR_GAUSS=2

! ANGINT3A end
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER                                                           &
     &     I_NORMAL                                                     &
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL                                                  &
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION                                          &
!             CALCULATION ABORTED
     &   , I_MISSING_DATA                                               &
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO                                                     &
!             I/O ERROR
     &   , I_ERR_RANGE                                                  &
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(                                                        &
     &     I_NORMAL=0                                                   &
     &   , I_ERR_FATAL=1                                                &
     &   , I_ABORT_CALCULATION=2                                        &
     &   , I_MISSING_DATA=3                                             &
     &   , I_ERR_IO=4                                                   &
     &   , I_ERR_RANGE=5                                                &
     &   , I_ERR_EXIST=6                                                &
     &   )
!
!     ------------------------------------------------------------------
!
!
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!                       Atmospheric Properties
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      REAL                                                              &
                !, INTENT(IN)
     &     D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESS OF EACH LAYER
!
!                       Angular Integration
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_ANGULAR_INTEGRATION                                        &
!             ANGULAR INTEGRATION SCHEME
     &   , I_2STREAM                                                    &
!             TWO-STREAM SCHEME
     &   , N_ORDER_GAUSS
!             ORDER OF GAUSSIAN INTEGRATION
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_2_STREAM_CORRECT                                           &
!             CORRECTION TO TWO-STREAM SCHEME
     &   , L_RESCALE
!             RESCALE OPTICAL PROPERTIES
!
!                       Treatment of Scattering
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SCATTER_METHOD_BAND
!
!                       Options for Solver
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SOLVER                                                     &
!             SOLVER USED
     &   , N_AUGMENT
!             LENGTH OF LONG FLUX VECTOR
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_NET
!             CALCULATE NET FLUXES
!
!                       Gaseous Properties
      REAL                                                              &
                !, INTENT(IN)
     &     K_GAS_ABS(NPD_PROFILE, NPD_LAYER)
!             GASEOUS ABSORPTIVE EXTINCTIONS
!
!                       Variables for Equivalent Extinction
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_SCALE_SOLAR
!             APPLY SCALING TO SOLAR FLUX
      REAL                                                              &
                !, INTENT(IN)
     &     ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             ADJUSTMENT OF SOLAR BEAM WITH EQUIVALENT EXTINCTION
!
!                       Spectral Region
      INTEGER                                                           &
                !, INTENT(IN)
     &     ISOLIR
!             VISIBLE OR IR
!
!                       Infra-red Properties
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             FLAG FOR QUADRATIC IR-SOURCE
      REAL                                                              &
                !, INTENT(IN)
     &     PLANCK_SOURCE(NPD_PROFILE, 0: NPD_LAYER)                     &
!             MONOCHROMATIC PLANCKIAN SOURCE
     &   , DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)                          &
!             THERMAL SOURCE FUNCTION
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)
!             2ND DIFF. OF THERMAL SOURCE FUNCTION
!
!                       Conditions at TOA
      REAL                                                              &
                !, INTENT(IN)
     &     SEC_0(NPD_PROFILE)                                           &
!             SECANT OF SOLAR ZENITH ANGLE
     &   , FLUX_INC_DIRECT(NPD_PROFILE)                                 &
!             INCIDENT DIRECT FLUX
     &   , FLUX_INC_DOWN(NPD_PROFILE)
!             INCIDENT DOWNWARD FLUX
!
!                       Surface Propeties
      REAL                                                              &
                !, INTENT(IN)
     &     ALBEDO_SURFACE_DIFF(NPD_PROFILE)                             &
!             DIFFUSE SURFACE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)                              &
!             DIRECT SURFACE ALBEDO
     &   , SOURCE_GROUND(NPD_PROFILE)
!             GROUND SOURCE FUNCTION
      REAL                                                              &
                !, INTENT(IN)
     &     GROUND_EMISSION(NPD_PROFILE)
!             TOTAL FLUX EMITTED FROM GROUND
!
!                       Optical Properties
      REAL                                                              &
                !, INTENT(IN)
     &     K_GREY_TOT_FREE(NPD_PROFILE, NPD_LAYER)                      &
!             FREE ABSORPTIVE EXTINCTION
     &   , K_EXT_SCAT_FREE(NPD_PROFILE, NPD_LAYER)                      &
!             FREE SCATTERING EXTINCTION
     &   , ASYMMETRY_FREE(NPD_PROFILE, NPD_LAYER)                       &
!             CLEAR-SKY ASYMMETRY
     &   , FORWARD_SCATTER_FREE(NPD_PROFILE, NPD_LAYER)
!             FREE FORWARD SCATTERING
!
!                       Cloudy Properties
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLOUD
!             CLOUDS REQUIRED
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_CLOUD
!             CLOUD SCHEME USED
!
!                       Cloud Geometry
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CLOUD_TOP                                                  &
!             TOPMOST CLOUDY LAYER
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
     &   , N_COLUMN(NPD_PROFILE)                                        &
!             NUMBER OF COLUMNS REQUIRED
     &   , N_REGION                                                     &
!             Number of cloudy regions
     &   , I_REGION_CLOUD(NPD_CLOUD_TYPE)
!             REGIONS IN WHICH TYPES OF CLOUDS FALL
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_COLUMN(NPD_PROFILE, NPD_LAYER, NPD_COLUMN)
!             FLAGS FOR CONTENTS OF COLUMNS
      REAL                                                              &
                !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             CLOUDY FRACTION
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)           &
!             FRACTIONS OF DIFFERENT TYPES OF CLOUD
     &   , W_FREE(NPD_PROFILE, NPD_LAYER)                               &
!             CLEAR-SKY FRACTION
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)  &
!             COEFFICIENTS FOR ENERGY TRANSFER AT INTERFACES
     &   , AREA_COLUMN(NPD_PROFILE, NPD_COLUMN)                         &
!             AREAS OF COLUMNS
     &   , FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTIONS OF TOTAL CLOUD OCCUPIED BY EACH REGION
!
!                       Cloudy Optical Properties
      REAL                                                              &
                !, INTENT(IN)
     &     K_GREY_TOT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)     &
!             CLOUDY ABSORPTIVE EXTINCTION
     &   , K_EXT_SCAT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)     &
!             CLOUDY SCATTERING EXTINCTION
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)      &
!             CLOUDY ASYMMETRY
     &   , FORWARD_SCATTER_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY FORWARD SCATTERING
!
!                       Fluxes Calculated
      REAL                                                              &
                !, INTENT(OUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DIRECT FLUX
     &   , FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX
!
!                       Flags for Clear-sky Calculations
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLEAR
!             CALCULATE CLEAR-SKY PROPERTIES
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SOLVER_CLEAR
!             CLEAR SOLVER USED
!
!                       Clear-sky Fluxes Calculated
      REAL                                                              &
                !, INTENT(OUT)
     &     FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                 &
!             CLEAR-SKY DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR-SKY TOTAL FLUX
!
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     K
!             LOOP VARIABLE
      REAL                                                              &
     &     TAU_FREE(NPD_PROFILE, NPD_LAYER)                             &
!             FREE OPTICAL DEPTH
     &   , OMEGA_FREE(NPD_PROFILE, NPD_LAYER)                           &
!             FREE ALBEDO OF S. S.
     &   , TAU_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)            &
!             CLOUDY OPTICAL DEPTH
     &   , OMEGA_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY SINGLE SCATTERING ALBEDO
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     SINGLE_SCATTERING_ALL, RESCALE_TAU_OMEGA                     &
     &   , TWO_STREAM                                                   &
     &   , MIX_COLUMN, CLOUD_COLUMN                                     &
     &   , GAUSS_ANGLE
!
!
!
!     CALCULATE SINGLE SCATTERING PROPERTIES FOR ALL ATMOSPHERIC
!     CONSTITUENTS.
!
! DEPENDS ON: single_scattering_all
      CALL SINGLE_SCATTERING_ALL(I_SCATTER_METHOD_BAND                  &
!                       Atmospheric Properties
     &   , N_PROFILE, N_LAYER, D_MASS                                   &
!                       Cloudy Properties
     &   , L_CLOUD, N_CLOUD_TOP, N_CLOUD_TYPE                           &
!                       Optical Properties
     &   , K_GREY_TOT_FREE, K_EXT_SCAT_FREE                             &
     &   , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                           &
     &   , K_GAS_ABS                                                    &
!                       Single Scattering Properties
     &   , TAU_FREE, OMEGA_FREE                                         &
     &   , TAU_CLOUD, OMEGA_CLOUD                                       &
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
!
      IF (I_ANGULAR_INTEGRATION == IP_TWO_STREAM) THEN
!
!        RESCALE TAU AND OMEGA. THE ASYMMETRY HAS ALREADY BEEN RESCALED.
!
         IF (L_RESCALE) THEN
!
! DEPENDS ON: rescale_tau_omega
            CALL RESCALE_TAU_OMEGA(N_PROFILE, 1, N_LAYER                &
     &         , TAU_FREE, OMEGA_FREE, FORWARD_SCATTER_FREE             &
     &         , NPD_PROFILE, NPD_LAYER                                 &
     &         )
!
            IF (L_CLOUD) THEN
!
               DO K=1, N_CLOUD_TYPE
! DEPENDS ON: rescale_tau_omega
                  CALL RESCALE_TAU_OMEGA(N_PROFILE, N_CLOUD_TOP         &
     &               , N_LAYER                                          &
     &               , TAU_CLOUD(1, 1, K), OMEGA_CLOUD(1, 1, K)         &
     &               , FORWARD_SCATTER_CLOUD(1, 1, K)                   &
     &               , NPD_PROFILE, NPD_LAYER                           &
     &               )
               ENDDO
!
            ENDIF
!
         ENDIF
!
!
!        SOLVE THE EQUATIONS USING THE SCHEME INDICATED BY THE VALUES
!        OF I_CLOUD AND I_SOLVER.
         IF (I_CLOUD == IP_CLOUD_CLEAR) THEN
!
!           A TWO-STREAM SCHEME WITH NO CLOUDS.
! DEPENDS ON: two_stream
            CALL TWO_STREAM(IERR                                        &
!                       Atmospheric Properties
     &         , N_PROFILE, N_LAYER                                     &
!                       Two-stream Scheme
     &         , I_2STREAM                                              &
!                       Corrections to Two-stream Equations
     &         , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION     &
!                       Options for Solver
     &         , L_NET, I_SOLVER                                        &
!                       Options for Equivalent Extinction
     &         , L_SCALE_SOLAR, ADJUST_SOLAR_KE                         &
!                       Spectral Region
     &         , ISOLIR                                                 &
!                       Infra-red Properties
     &         , DIFF_PLANCK                                            &
     &         , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                        &
!                       Conditions at TOA
     &         , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                  &
!                       Surface Conditions
     &         , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND &
!                       Single Scattering Propeties
     &         , TAU_FREE, OMEGA_FREE, ASYMMETRY_FREE                   &
!                       Fluxes Calculated
     &         , FLUX_DIRECT, FLUX_TOTAL                                &
!                       Flag for Clear-sky Fluxes
     &         , L_CLEAR                                                &
!                       Clear-sky Fluxes Calculated
     &         , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                    &
!                       Sizes of Arrays
     &         , NPD_PROFILE, NPD_LAYER                                 &
     &         )
               IF (IERR /= I_NORMAL) RETURN
!
      ELSEIF ( (I_CLOUD == IP_CLOUD_MIX_MAX).OR.                        &
     &         (I_CLOUD == IP_CLOUD_MIX_RANDOM).OR.                     &
     &         ( (I_CLOUD == IP_CLOUD_PART_CORR).AND.                   &
     &           (N_REGION == 2) ) ) THEN
!
!           CLOUDS ARE TREATED USING ZDUNKOWSKI'S MIXED-COLUMN SCHEME.
!           THE GEOMETRY HAS BEEN SET BEFORE.
!
! DEPENDS ON: mix_column
            CALL MIX_COLUMN(IERR                                        &
!                       Atmospheric Properties
     &         , N_PROFILE, N_LAYER                                     &
!                       Two-stream Scheme
     &         , I_2STREAM                                              &
!                       Corrections to Two-stream Equations
     &         , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION     &
!                       Options for Solver
     &         , I_SOLVER, L_NET                                        &
!                       Options for Equivalent Extinction
     &         , L_SCALE_SOLAR, ADJUST_SOLAR_KE                         &
!                       Spectral Region
     &         , ISOLIR                                                 &
!                       Infra-red Properties
     &         , DIFF_PLANCK                                            &
     &         , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                        &
!                       Conditions at TOA
     &         , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                  &
!                       Conditions at Surface
     &         , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND &
!                       Clear-sky Single Scattering Properties
     &         , TAU_FREE, OMEGA_FREE, ASYMMETRY_FREE                   &
!                       Cloud Geometry
     &         , N_CLOUD_TOP                                            &
     &         , N_CLOUD_TYPE, FRAC_CLOUD                               &
     &         , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE                 &
     &         , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE              &
     &         , CLOUD_OVERLAP                                          &
!                       Cloudy Optical Properties
     &         , TAU_CLOUD, OMEGA_CLOUD, ASYMMETRY_CLOUD                &
!                       Fluxes Calculated
     &         , FLUX_DIRECT, FLUX_TOTAL                                &
!                       Flags for Clear-sky Calculations
     &         , L_CLEAR, I_SOLVER_CLEAR                                &
!                       Clear-sky Fluxes Calculated
     &         , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                    &
!                       Dimensions of Arrays
     &         , NPD_PROFILE, NPD_LAYER                                 &
     &         )
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
            CALL TRIPLE_COLUMN(IERR                                     &
!                       Atmospheric Properties
     &         , N_PROFILE, N_LAYER                                     &
!                       Two-stream Scheme
     &         , I_2STREAM                                              &
!                       Corrections to Two-stream Equations
     &         , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION     &
!                       Options for Solver
     &         , I_SOLVER, L_NET, I_SCATTER_METHOD_BAND                 &
!                       Options for Equivalent Extinction
     &         , L_SCALE_SOLAR, ADJUST_SOLAR_KE                         &
!                       Spectral Region
     &         , ISOLIR                                                 &
!                       Infra-red Properties
     &         , DIFF_PLANCK                                            &
     &         , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                        &
!                       Conditions at TOA
     &         , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                  &
!                       Conditions at Surface
     &         , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND &
!                       Clear-sky Single Scattering Properties
     &         , TAU_FREE, OMEGA_FREE, ASYMMETRY_FREE                   &
!                       Cloud Geometry
     &         , N_CLOUD_TOP                                            &
     &         , N_CLOUD_TYPE, FRAC_CLOUD                               &
     &         , N_REGION, I_REGION_CLOUD, FRAC_REGION                  &
     &         , W_FREE, W_CLOUD                                        &
     &         , CLOUD_OVERLAP                                          &
!                       Cloudy Optical Properties
     &         , TAU_CLOUD, OMEGA_CLOUD, ASYMMETRY_CLOUD                &
!                       Fluxes Calculated
     &         , FLUX_DIRECT, FLUX_TOTAL                                &
!                       Flags for Clear-sky Calculations
     &         , L_CLEAR, I_SOLVER_CLEAR                                &
!                       Clear-sky Fluxes Calculated
     &         , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                    &
!                       Dimensions of Arrays
     &         , NPD_PROFILE, NPD_LAYER                                 &
     &         )
            IF (IERR /= I_NORMAL) RETURN
!
         ELSEIF (I_CLOUD == IP_CLOUD_COLUMN_MAX) THEN
!           CLOUDS ARE TREATED ON THE ASSUMPTION OF MAXIMUM OVERLAP
!           IN A COLUMN MODEL.
! DEPENDS ON: cloud_column
            CALL CLOUD_COLUMN(IERR                                      &
!                       Atmospheric Properties
     &         , N_PROFILE, N_LAYER                                     &
!                       Two-stream Scheme
     &         , I_2STREAM                                              &
!                       Corrections to Two-stream Equations
     &         , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION     &
!                       Options for Solver
     &         , I_SOLVER, N_AUGMENT                                    &
!                       Options for Equivalent Extinction
     &         , L_SCALE_SOLAR, ADJUST_SOLAR_KE                         &
!                       Spectral Region
     &         , ISOLIR                                                 &
!                       Infra-red Properties
     &         , DIFF_PLANCK                                            &
     &         , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                        &
!                       Conditions at TOA
     &         , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                  &
!                       Conditions at Surface
     &         , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND &
!                       Clear-sky Single Scattering Properties
     &         , TAU_FREE, OMEGA_FREE, ASYMMETRY_FREE                   &
!                       Cloud Geometry
     &         , N_CLOUD_TOP                                            &
     &         , N_CLOUD_TYPE, FRAC_CLOUD                               &
     &         , N_COLUMN, L_COLUMN, AREA_COLUMN                        &
!                       Cloudy Optical Properties
     &         , TAU_CLOUD, OMEGA_CLOUD, ASYMMETRY_CLOUD                &
!                       Fluxes Calculated
     &         , FLUX_DIRECT, FLUX_TOTAL                                &
!                       Flags for Clear-sky Calculations
     &         , L_CLEAR, I_SOLVER_CLEAR                                &
!                       Clear-sky Fluxes Calculated
     &         , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                    &
!                       Dimensions of Arrays
     &         , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                     &
     &         )
            IF (IERR /= I_NORMAL) RETURN
!
         ENDIF
!
      ELSE IF (I_ANGULAR_INTEGRATION == IP_IR_GAUSS) THEN
!
!        FULL ANGULAR RESOLUTION USING GASUSSIAN INTEGRATION.
! DEPENDS ON: gauss_angle
         CALL GAUSS_ANGLE(N_PROFILE, N_LAYER, L_NET, N_AUGMENT          &
     &      , N_ORDER_GAUSS                                             &
     &      , TAU_FREE                                                  &
     &      , FLUX_INC_DOWN                                             &
     &      , DIFF_PLANCK, SOURCE_GROUND, ALBEDO_SURFACE_DIFF           &
     &      , FLUX_TOTAL                                                &
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                           &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
         IF (IERR /= I_NORMAL) RETURN
!
      ENDIF
!
!
      RETURN
      END SUBROUTINE MONOCHROMATIC_FLUX
