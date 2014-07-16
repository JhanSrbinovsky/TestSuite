


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve the two-stream equations in a triple column.
!
! Method:
!       The atmospheric column is divided into three regions
!       in each layer and the two-stream coefficients are determined
!       for each region. The equations are then solved using
!       appropriate coupling of the fluxes at the boundaries
!       of layers.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.2             15-05-96                Original Code
!                                               (J. M. Edwards)
!       4.5             18-05-98                Variable for obsolete
!                                               solver removed. EXTERNAL
!                                               statement corrected.
!                                               Unused variables
!                                               removed.
!                                               (J. M. Edwards)
!       5.1             04-04-00                Obsolete header file
!                                               removed.
!                                               (J. M. Edwards)
!       5.3             04-10-01                Number of regions
!                                               passed explicitly.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TRIPLE_COLUMN(IERR                                     &
!                       Atmospheric Properties
     &   , N_PROFILE, N_LAYER                                           &
!                       Two-stream Scheme
     &   , I_2STREAM                                                    &
!                       Corrections to Two-stream Equations
     &   , L_2_STREAM_CORRECT, PLANCK_SOURCE, GROUND_EMISSION           &
!                       Options for Solver
     &   , I_SOLVER, L_NET, I_SCATTER_METHOD                            &
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
     &   , N_REGION, I_REGION_CLOUD, FRAC_REGION                        &
     &   , W_FREE, W_CLOUD                                              &
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
! STDIO3A defines unit numbers for standard i/o in two-stream radiation
! code.
      INTEGER,PARAMETER:: IU_STDIN=5
      INTEGER,PARAMETER:: IU_STDOUT=6
      INTEGER,PARAMETER:: IU_ERR=6
! STDIO3A end
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
! SCTMTH3A defines the methods of treating scattering for two-stream
! radiation code.
!     Current Owner of Code: J. M. Edwards
!
!     History:
!     Version  Date      Comment.
!     5.3      04/10/01  Obsolete option removed.
!                                J. M. Edwards
!

      ! full treatment of scattering
      INTEGER,PARAMETER:: IP_SCATTER_FULL=1

      ! scattering ignored completely.
      INTEGER,PARAMETER:: IP_NO_SCATTER_ABS=2

! SCTMTH3A end
! SPCRG3A defines flags for different portions of the spectrum in
! two-stream radiation code.
      INTEGER,PARAMETER:: IP_SOLAR=1
      INTEGER,PARAMETER:: IP_INFRA_RED=2
! SPCRG3A end
! SOLVER3A defines reference numbers for solvers for two-stream
! radiation code.

      ! pentadiagonal scheme
      INTEGER,PARAMETER:: IP_SOLVER_PENTADIAGONAL=1

      ! mixed column scheme using full endekadiagonal matrix
      INTEGER,PARAMETER:: IP_SOLVER_MIX_11=6

      ! mixed column scheme with approximate scattering
      INTEGER,PARAMETER:: IP_SOLVER_MIX_APP_SCAT=9

      ! direct mixed column scheme for full fluxes
      INTEGER,PARAMETER:: IP_SOLVER_MIX_DIRECT=11

      ! direct solver for a homogeneous column
      INTEGER,PARAMETER:: IP_SOLVER_HOMOGEN_DIRECT=13

      ! direct solver for triple column
      INTEGER,PARAMETER:: IP_SOLVER_TRIPLE=14

      ! direct solver for triple column approximating scattering
      INTEGER,PARAMETER:: IP_SOLVER_TRIPLE_APP_SCAT=15

      ! direct mixed column scheme for full fluxes (modified
      !   for correct treatment of shadowing by Robin Hogan)
      INTEGER,PARAMETER:: IP_SOLVER_MIX_DIRECT_HOGAN=16

      ! direct solver for triple column (modified for
      !   correct treatment of shadowing by Robin Hogan)
      INTEGER,PARAMETER:: IP_SOLVER_TRIPLE_HOGAN=17

! SOLVER3A end
! CLCFPT3A defines pointers in CLOUD_OVERLAP for two-stream radiation
! code.
!
! note that several pointers are identical since only certain
! groups of coefficients are relevant to a particular scheme.

      INTEGER,PARAMETER:: IP_CLOVLP_GFF=1 ! GAMMA-FREE-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_GFC=2 ! GAMMA-FREE-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_GCF=3 ! GAMMA-CLOUD-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_GCC=4 ! GAMMA-CLOUD-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_BFF=5 ! BETA-FREE-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_BFC=6 ! BETA-FREE-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_BCF=7 ! BETA-CLOUD-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_BCC=8 ! BETA-CLOUD-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_GFM=5 ! GAMMA_F-
      INTEGER,PARAMETER:: IP_CLOVLP_GFP=6 ! GAMMA_F+
      INTEGER,PARAMETER:: IP_CLOVLP_BFM=7 ! BETA_F-
      INTEGER,PARAMETER:: IP_CLOVLP_BFP=8 ! BETA_F+
      INTEGER,PARAMETER:: IP_CLOVLP_GM=5  ! GAMMA_-
      INTEGER,PARAMETER:: IP_CLOVLP_GP=6  ! GAMMA_+
      INTEGER,PARAMETER:: IP_CLOVLP_BM=7  ! BETA_-
      INTEGER,PARAMETER:: IP_CLOVLP_BP=8  ! BETA_+

      ! pointers for triple overlaps:

      INTEGER,PARAMETER:: IP_CLOVLP_V11=1
      INTEGER,PARAMETER:: IP_CLOVLP_V12=2
      INTEGER,PARAMETER:: IP_CLOVLP_V13=3
      INTEGER,PARAMETER:: IP_CLOVLP_V21=4
      INTEGER,PARAMETER:: IP_CLOVLP_V22=5
      INTEGER,PARAMETER:: IP_CLOVLP_V23=6
      INTEGER,PARAMETER:: IP_CLOVLP_V31=7
      INTEGER,PARAMETER:: IP_CLOVLP_V32=8
      INTEGER,PARAMETER:: IP_CLOVLP_V33=9
      INTEGER,PARAMETER:: IP_CLOVLP_U11=10
      INTEGER,PARAMETER:: IP_CLOVLP_U12=11
      INTEGER,PARAMETER:: IP_CLOVLP_U13=12
      INTEGER,PARAMETER:: IP_CLOVLP_U21=13
      INTEGER,PARAMETER:: IP_CLOVLP_U22=14
      INTEGER,PARAMETER:: IP_CLOVLP_U23=15
      INTEGER,PARAMETER:: IP_CLOVLP_U31=16
      INTEGER,PARAMETER:: IP_CLOVLP_U32=17
      INTEGER,PARAMETER:: IP_CLOVLP_U33=18

! CLCFPT3A end
! CLDREG3A defines reference numbers for regions of clouds.in two-stream
! radiation code.

      INTEGER,PARAMETER:: IP_REGION_CLEAR=1 ! clear-sky region
      INTEGER,PARAMETER:: IP_REGION_STRAT=2 ! stratiform cloudy region
      INTEGER,PARAMETER:: IP_REGION_CONV=3  ! convective cloudy region

! CLDREG3A end
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
     &   , ISOLIR                                                       &
!             SPECTRAL REGION
     &   , I_2STREAM                                                    &
!             TWO-STREAM SCHEME
     &   , I_SOLVER                                                     &
!             SOLVER USED
     &   , I_SOLVER_CLEAR                                               &
!             SOLVER FOR CLEAR-SKY FLUXES
     &   , I_SCATTER_METHOD
!             Method of treating scattering
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
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_REGION                                                     &
!             Number of cloudy regions
     &   , I_REGION_CLOUD(NPD_CLOUD_TYPE)
!             REGIONS IN WHICH TYPES OF CLOUDS FALL
      REAL                                                              &
                !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             CLOUDY FRACTIONS IN EACH LAYER
     &   , W_FREE(NPD_PROFILE, NPD_LAYER)                               &
!             CLEAR SKY FRACTIONS IN EACH LAYER
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)           &
!             FRACTIONS OF DIFFERENT TYPES OF CLOUD
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)  &
!             ENERGY TRANSFER COEFFICIENTS
     &   , FRAC_REGION(NPD_PROFILE, NPD_LAYER, NPD_REGION)
!             FRACTIONS OF TOTAL CLOUD OCCUPIED BY EACH REGION
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
     &   , I                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
!             LOOP VARIABLE
     &   , N_TOP
!             TOP-MOST LAYER FOR CALCULATION
!
!
!     CLEAR-SKY COEFFICIENTS:
      REAL                                                              &
     &     TRANS(NPD_PROFILE, NPD_LAYER, NPD_REGION)                    &
!             TRANSMISSION COEFFICIENTS
     &   , REFLECT(NPD_PROFILE, NPD_LAYER, NPD_REGION)                  &
!             REFLECTION COEFFICIENTS
     &   , TRANS_0(NPD_PROFILE, NPD_LAYER, NPD_REGION)                  &
!             DIRECT TRANSMISSION COEFFICIENTS
     &   , SOURCE_COEFF(NPD_PROFILE, NPD_LAYER                          &
     &      , NPD_SOURCE_COEFF, NPD_REGION)                             &
!             SOURCE COEFFICIENTS
     &   , S_DOWN(NPD_PROFILE, NPD_LAYER, NPD_REGION)                   &
!             FREE DOWNWARD SOURCE
     &   , S_UP(NPD_PROFILE, NPD_LAYER, NPD_REGION)                     &
!             FREE UPWARD SOURCE
     &   , S_DOWN_CLEAR(NPD_PROFILE, NPD_LAYER)                         &
!             CLEAR DOWNWARD SOURCE
     &   , S_UP_CLEAR(NPD_PROFILE, NPD_LAYER)
!             CLEAR UPWARD SOURCE
!
!     SOURCE FUNCTIONS AT THE CROUND
      REAL                                                              &
     &     SOURCE_FLUX_GROUND(NPD_PROFILE, NPD_REGION)                  &
!             SOURCE OF FLUX FROM GROUND
     &   , FLUX_DIRECT_GROUND(NPD_PROFILE, NPD_REGION)
!             DIRECT FLUX AT GROUND IN EACH REGION
!
!
!     FUNCTIONS CALLED:
      INTEGER                                                           &
     &     SET_N_SOURCE_COEFF
!             FUNCTION TO SET NUMBER OF SOURCE COEFFICIENTS
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     TWO_COEFF_REGION, IR_SOURCE, TRIPLE_SOLAR_SOURCE             &
     &   , SOLVER_TRIPLE, SOLVER_TRIPLE_APP_SCAT                        &
     &   , CLEAR_SUPPLEMENT, TWO_COEFF_REGION_FAST_LW
!
!
!
!     SET THE NUMBER OF SOURCE COEFFICIENTS FOR THE APPROXIMATION
! DEPENDS ON: set_n_source_coeff
      N_SOURCE_COEFF=SET_N_SOURCE_COEFF(ISOLIR, L_IR_SOURCE_QUAD)
!
!
      IF (I_SCATTER_METHOD == IP_SCATTER_FULL) THEN
! DEPENDS ON: two_coeff_region
         CALL TWO_COEFF_REGION(IERR                                     &
     &      , N_PROFILE, N_LAYER, N_CLOUD_TOP                           &
     &      , I_2STREAM, L_IR_SOURCE_QUAD, N_SOURCE_COEFF               &
     &      , N_CLOUD_TYPE, FRAC_CLOUD                                  &
     &      , N_REGION, I_REGION_CLOUD, FRAC_REGION                     &
     &      , ASYMMETRY_FREE, OMEGA_FREE, TAU_FREE                      &
     &      , ASYMMETRY_CLOUD, OMEGA_CLOUD, TAU_CLOUD                   &
     &      , ISOLIR, SEC_0                                             &
     &      , TRANS, REFLECT, TRANS_0, SOURCE_COEFF                     &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
      ELSE IF (I_SCATTER_METHOD == IP_NO_SCATTER_ABS) THEN
! DEPENDS ON: two_coeff_region_fast_lw
         CALL TWO_COEFF_REGION_FAST_LW(IERR                             &
     &      , N_PROFILE, N_LAYER, N_CLOUD_TOP                           &
     &      , I_2STREAM, L_IR_SOURCE_QUAD, N_SOURCE_COEFF               &
     &      , N_CLOUD_TYPE, FRAC_CLOUD                                  &
     &      , N_REGION, I_REGION_CLOUD, FRAC_REGION                     &
     &      , TAU_FREE, TAU_CLOUD                                       &
     &      , ISOLIR                                                    &
     &      , TRANS, REFLECT, SOURCE_COEFF                              &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
      ENDIF
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
         DO K=1, N_REGION
            IF (K == IP_REGION_CLEAR) THEN
               N_TOP=1
            ELSE
               N_TOP=N_CLOUD_TOP
            ENDIF
!
! DEPENDS ON: ir_source
            CALL IR_SOURCE(N_PROFILE, N_TOP, N_LAYER                    &
     &         , SOURCE_COEFF(1, 1, 1, K), DIFF_PLANCK                  &
     &         , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                        &
     &         , L_2_STREAM_CORRECT, PLANCK_SOURCE                      &
     &         , GROUND_EMISSION, N_LAYER                               &
     &         , TAU_FREE, TRANS                                        &
     &         , S_DOWN(1, 1, K), S_UP(1, 1, K)                         &
     &         , NPD_PROFILE, NPD_LAYER                                 &
     &         )
         ENDDO
!
!
!        WEIGHT THE SOURCE FUNCTIONS BY THE AREA FRACTIONS, BUT
!        SAVE THE CLEAR-SKY FRACTIONS FOR DIAGNOSTIC USE IF
!        REQUIRED.
         IF (L_CLEAR) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  S_DOWN_CLEAR(L, I)=S_DOWN(L, I, IP_REGION_CLEAR)
                  S_UP_CLEAR(L, I)=S_UP(L, I, IP_REGION_CLEAR)
               ENDDO
            ENDDO
         ENDIF
         DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
               S_DOWN(L, I, IP_REGION_CLEAR)                            &
     &            =W_FREE(L, I)*S_DOWN(L, I, IP_REGION_CLEAR)
               S_UP(L, I, IP_REGION_CLEAR)                              &
     &            =W_FREE(L, I)*S_UP(L, I, IP_REGION_CLEAR)
               S_DOWN(L, I, IP_REGION_STRAT)                            &
     &            =W_CLOUD(L, I)                                        &
     &            *FRAC_REGION(L, I, IP_REGION_STRAT)                   &
     &            *S_DOWN(L, I, IP_REGION_STRAT)
               S_UP(L, I, IP_REGION_STRAT)                              &
     &            =W_CLOUD(L, I)                                        &
     &            *FRAC_REGION(L, I, IP_REGION_STRAT)                   &
     &            *S_UP(L, I, IP_REGION_STRAT)
               S_DOWN(L, I, IP_REGION_CONV)                             &
     &            =W_CLOUD(L, I)                                        &
     &            *FRAC_REGION(L, I, IP_REGION_CONV)                    &
     &            *S_DOWN(L, I, IP_REGION_CONV)
               S_UP(L, I, IP_REGION_CONV)                               &
     &            =W_CLOUD(L, I)                                        &
     &            *FRAC_REGION(L, I, IP_REGION_CONV)                    &
     &            *S_UP(L, I, IP_REGION_CONV)
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
! DEPENDS ON: triple_solar_source
         CALL TRIPLE_SOLAR_SOURCE(N_PROFILE, N_LAYER, N_CLOUD_TOP       &
     &      , N_REGION, FLUX_INC_DIRECT                                 &
     &      , L_SCALE_SOLAR, ADJUST_SOLAR_KE                            &
     &      , TRANS_0, SOURCE_COEFF                                     &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V11)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V12)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V13)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V21)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V22)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V23)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V31)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V32)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V33)                        &
     &      , FLUX_DIRECT, FLUX_DIRECT_GROUND                           &
     &      , S_UP, S_DOWN                                              &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &   )
      ENDIF
!
!        SET THE PARTITIONED SOURCE FUNCTIONS AT THE GROUND.
         IF (ISOLIR == IP_SOLAR) THEN
            DO L=1, N_PROFILE
               SOURCE_FLUX_GROUND(L, IP_REGION_CLEAR)                   &
     &            =(ALBEDO_SURFACE_DIR(L)-ALBEDO_SURFACE_DIFF(L))       &
     &            *FLUX_DIRECT_GROUND(L, IP_REGION_CLEAR)
               SOURCE_FLUX_GROUND(L, IP_REGION_STRAT)                   &
     &            =(ALBEDO_SURFACE_DIR(L)-ALBEDO_SURFACE_DIFF(L))       &
     &            *FLUX_DIRECT_GROUND(L, IP_REGION_STRAT)
               SOURCE_FLUX_GROUND(L, IP_REGION_CONV)                    &
     &            =(ALBEDO_SURFACE_DIR(L)-ALBEDO_SURFACE_DIFF(L))       &
     &            *FLUX_DIRECT_GROUND(L, IP_REGION_CONV)
            ENDDO
         ELSE
            DO L=1, N_PROFILE
               SOURCE_FLUX_GROUND(L, IP_REGION_CLEAR)                   &
     &            =CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_U11)             &
     &            *SOURCE_GROUND(L)
               SOURCE_FLUX_GROUND(L, IP_REGION_STRAT)                   &
     &            =CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_U21)             &
     &            *SOURCE_GROUND(L)
               SOURCE_FLUX_GROUND(L, IP_REGION_CONV)                    &
     &            =CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_U31)             &
     &            *SOURCE_GROUND(L)
            ENDDO
         ENDIF
!
!
!
      SELECT CASE (i_solver)

      CASE (IP_solver_triple)
!
! DEPENDS ON: solver_triple
         CALL SOLVER_TRIPLE(N_PROFILE, N_LAYER, N_CLOUD_TOP             &
     &      , TRANS(1, 1, IP_REGION_CLEAR)                              &
     &      , REFLECT(1, 1, IP_REGION_CLEAR)                            &
     &      , S_DOWN(1, 1, IP_REGION_CLEAR)                             &
     &      , S_UP(1, 1, IP_REGION_CLEAR)                               &
     &      , TRANS(1, 1, IP_REGION_STRAT)                              &
     &      , REFLECT(1, 1, IP_REGION_STRAT)                            &
     &      , S_DOWN(1, 1, IP_REGION_STRAT)                             &
     &      , S_UP(1, 1, IP_REGION_STRAT)                               &
     &      , TRANS(1, 1, IP_REGION_CONV)                               &
     &      , REFLECT(1, 1, IP_REGION_CONV)                             &
     &      , S_DOWN(1, 1, IP_REGION_CONV)                              &
     &      , S_UP(1, 1, IP_REGION_CONV)                                &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V11)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V12)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V13)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V21)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V22)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V23)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V31)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V32)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V33)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U11)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U12)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U13)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U21)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U22)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U23)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U31)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U32)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U33)                        &
     &      , L_NET                                                     &
     &      , FLUX_INC_DOWN                                             &
     &      , SOURCE_FLUX_GROUND(1, IP_REGION_CLEAR)                    &
     &      , SOURCE_FLUX_GROUND(1, IP_REGION_STRAT)                    &
     &      , SOURCE_FLUX_GROUND(1, IP_REGION_CONV)                     &
     &      , ALBEDO_SURFACE_DIFF                                       &
     &      , FLUX_TOTAL                                                &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
!
      CASE (IP_solver_triple_hogan)
!
! DEPENDS ON: solver_triple_hogan
         CALL solver_triple_hogan(N_PROFILE, N_LAYER, N_CLOUD_TOP       &
     &      , TRANS(1, 1, IP_REGION_CLEAR)                              &
     &      , REFLECT(1, 1, IP_REGION_CLEAR)                            &
     &      , S_DOWN(1, 1, IP_REGION_CLEAR)                             &
     &      , S_UP(1, 1, IP_REGION_CLEAR)                               &
     &      , TRANS(1, 1, IP_REGION_STRAT)                              &
     &      , REFLECT(1, 1, IP_REGION_STRAT)                            &
     &      , S_DOWN(1, 1, IP_REGION_STRAT)                             &
     &      , S_UP(1, 1, IP_REGION_STRAT)                               &
     &      , TRANS(1, 1, IP_REGION_CONV)                               &
     &      , REFLECT(1, 1, IP_REGION_CONV)                             &
     &      , S_DOWN(1, 1, IP_REGION_CONV)                              &
     &      , S_UP(1, 1, IP_REGION_CONV)                                &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V11)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V12)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V13)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V21)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V22)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V23)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V31)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V32)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V33)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U11)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U12)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U13)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U21)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U22)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U23)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U31)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U32)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U33)                        &
     &      , L_NET                                                     &
     &      , FLUX_INC_DOWN                                             &
     &      , SOURCE_FLUX_GROUND(1, IP_REGION_CLEAR)                    &
     &      , SOURCE_FLUX_GROUND(1, IP_REGION_STRAT)                    &
     &      , SOURCE_FLUX_GROUND(1, IP_REGION_CONV)                     &
     &      , ALBEDO_SURFACE_DIFF                                       &
     &      , FLUX_TOTAL                                                &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
!
      CASE (IP_solver_triple_app_scat)
!
! DEPENDS ON: solver_triple_app_scat
         CALL SOLVER_TRIPLE_APP_SCAT(N_PROFILE, N_LAYER, N_CLOUD_TOP    &
     &      , TRANS(1, 1, IP_REGION_CLEAR)                              &
     &      , REFLECT(1, 1, IP_REGION_CLEAR)                            &
     &      , S_DOWN(1, 1, IP_REGION_CLEAR)                             &
     &      , S_UP(1, 1, IP_REGION_CLEAR)                               &
     &      , TRANS(1, 1, IP_REGION_STRAT)                              &
     &      , REFLECT(1, 1, IP_REGION_STRAT)                            &
     &      , S_DOWN(1, 1, IP_REGION_STRAT)                             &
     &      , S_UP(1, 1, IP_REGION_STRAT)                               &
     &      , TRANS(1, 1, IP_REGION_CONV)                               &
     &      , REFLECT(1, 1, IP_REGION_CONV)                             &
     &      , S_DOWN(1, 1, IP_REGION_CONV)                              &
     &      , S_UP(1, 1, IP_REGION_CONV)                                &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V11)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V12)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V13)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V21)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V22)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V23)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V31)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V32)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_V33)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U11)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U12)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U13)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U21)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U22)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U23)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U31)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U32)                        &
     &      , CLOUD_OVERLAP(1, 0, IP_CLOVLP_U33)                        &
     &      , L_NET                                                     &
     &      , FLUX_INC_DOWN                                             &
     &      , SOURCE_FLUX_GROUND(1, IP_REGION_CLEAR)                    &
     &      , SOURCE_FLUX_GROUND(1, IP_REGION_STRAT)                    &
     &      , SOURCE_FLUX_GROUND(1, IP_REGION_CONV)                     &
     &      , ALBEDO_SURFACE_DIFF                                       &
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
     &      , TRANS(1, 1, IP_REGION_CLEAR)                              &
     &      , REFLECT(1, 1, IP_REGION_CLEAR)                            &
     &      , TRANS_0(1, 1, IP_REGION_CLEAR)                            &
     &      , SOURCE_COEFF(1, 1, 1, IP_REGION_CLEAR)                    &
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
      END SUBROUTINE TRIPLE_COLUMN
