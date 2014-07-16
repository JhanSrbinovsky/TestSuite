


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
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             08-08-96                Code for vertically
!                                               coherent convective
!                                               cloud added.
!                                               (J. M. Edwards)
!       4.5             18-05-98                Variable for obsolete
!                                               solver removed.
!                                               (J. M. Edwards)
!       5.3             04-10-01                Number of cloudy regions
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
      SUBROUTINE SOLVE_BAND_WITHOUT_GAS(IERR                            &
!                       Atmospheric Column
     &   , N_PROFILE, N_LAYER, D_MASS                                   &
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT         &
     &   , L_RESCALE, N_ORDER_GAUSS                                     &
!                       Treatment of scattering
     &   , I_SCATTER_METHOD_BAND                                        &
!                       Options for solver
     &   , I_SOLVER, L_NET, N_AUGMENT                                   &
!                       Spectral region
     &   , ISOLIR                                                       &
!                       Solar properties
     &   , SEC_0, SOLAR_FLUX                                            &
!                       Infra-red properties
     &   , PLANCK_SOURCE_TOP, PLANCK_SOURCE_BOTTOM                      &
     &   , DIFF_PLANCK_BAND, L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2       &
!                       Surface properties
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, THERMAL_GROUND_BAND &
!                       Clear-sky optical properties
     &   , K_GREY_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE             &
     &   , FORWARD_SCATTER_FREE                                         &
!                       Cloudy properties
     &   , L_CLOUD, I_CLOUD                                             &
!                       Cloud Geometry
     &   , N_CLOUD_TOP                                                  &
     &   , N_CLOUD_TYPE, FRAC_CLOUD                                     &
     &   , N_REGION, I_REGION_CLOUD, FRAC_REGION                        &
     &   , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE                       &
     &   , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE                    &
     &   , CLOUD_OVERLAP                                                &
     &   , N_COLUMN, L_COLUMN, AREA_COLUMN                              &
!                       Cloudy optical properties
     &   , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                           &
     &   , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD                       &
!                       Calculated Fluxes
     &   , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                            &
!                       Flags For Clear-sky Fluxes
     &   , L_CLEAR, I_SOLVER_CLEAR                                      &
!                       Calculated Clear-sky Fluxes
     &   , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND                &
!                       Planckian Function
     &   , PLANCK_SOURCE_BAND                                           &
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                           &
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
! SPCRG3A defines flags for different portions of the spectrum in
! two-stream radiation code.
      INTEGER,PARAMETER:: IP_SOLAR=1
      INTEGER,PARAMETER:: IP_INFRA_RED=2
! SPCRG3A end
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!                       Atmospheric Column
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
!                       Angular integration
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
!             USE AN EDGE CORRECTION
     &   , L_RESCALE
!             RESCALE OPTICAL PROPERTIES
!
!                       Treatment of Scattering
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SCATTER_METHOD_BAND
!             METHOD OF TREATING SCATTERING
!
!                       Options for solver
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SOLVER                                                     &
!             SOLVER USED
     &   , N_AUGMENT
!             LENGTH OF LONG FLUX VECTOR
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_NET
!             SOLVE FOR NET FLUXES
!
!                       Spectral Region
      INTEGER                                                           &
                !, INTENT(IN)
     &     ISOLIR
!             VISIBLE OR IR
!
!                       Solar Properties
      REAL                                                              &
            !, INTENT(IN)
     &     SEC_0(NPD_PROFILE)                                           &
!             SECANT OF SOLAR ZENITH ANGLE
     &   , SOLAR_FLUX(NPD_PROFILE)
!             INCIDENT SOLAR FLUX
!
!                       Infra-red Properties
      REAL                                                              &
            !, INTENT(IN)
     &     PLANCK_SOURCE_TOP(NPD_PROFILE)                               &
!             PLANCK FUNCTION AT BOTTOM OF COLUMN
     &   , PLANCK_SOURCE_BOTTOM(NPD_PROFILE)                            &
!             PLANCK FUNCTION AT BOTTOM OF COLUMN
     &   , DIFF_PLANCK_BAND(NPD_PROFILE, NPD_LAYER)                     &
!             CHANGE IN PLANCK FUNCTION
     &   , DIFF_PLANCK_BAND_2(NPD_PROFILE, NPD_LAYER)
!             2x2ND DIFFERENCE OF PLANCKIAN IN BAND
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             USE A QUADRATIC SOURCE FUNCTION
!
!                       Surface Properties
      REAL                                                              &
            !, INTENT(IN)
     &     ALBEDO_SURFACE_DIFF(NPD_PROFILE)                             &
!             DIFFUSE SURFACE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)                              &
!             DIRECT SURFACE ALBEDO
     &   , THERMAL_GROUND_BAND(NPD_PROFILE)
!             THERMAL SOURCE AT SURFACE IN BAND
!
!                       Clear-sky Optical Properties
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
     &     W_FREE(NPD_PROFILE, NPD_LAYER)                               &
!             CLEAR-SKY FRACTION
     &   , W_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             CLOUDY FRACTION
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)           &
!             FRACTIONS OF TYPES OF CLOUDS
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)  &
!             COEFFICIENTS FOR TRANSFER FOR ENERGY AT INTERFACES
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
!                       Calculated Fluxes
      REAL                                                              &
            !, INTENT(OUT)
     &     FLUX_DIRECT_BAND(NPD_PROFILE, 0: NPD_LAYER)                  &
!             DIRECT FLUX
     &   , FLUX_TOTAL_BAND(NPD_PROFILE, 2*NPD_LAYER+2)
!
!                       Flags for Clear-sky Fluxes
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLEAR
!             CALCULATE NET CLEAR-SKY PROPERTIES
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SOLVER_CLEAR
!             CLEAR SOLVER USED
!
!                       Calculated Clear-sky Fluxes
      REAL                                                              &
            !, INTENT(OUT)
     &     FLUX_DIRECT_CLEAR_BAND(NPD_PROFILE, 0: NPD_LAYER)            &
!             CLEAR-SKY DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR_BAND(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR-SKY TOTAL FLUX
!
!                       Planckian Function
      REAL                                                              &
            !, INTENT(IN)
     &     PLANCK_SOURCE_BAND(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCKIAN SOURCE IN BAND
!
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      REAL                                                              &
     &     FLUX_INC_DIRECT(NPD_PROFILE)                                 &
!             INCIDENT DIRECT FLUX
     &   , FLUX_INC_DOWN(NPD_PROFILE)                                   &
!             INCIDENT DOWNWARD FLUX
     &   , SOURCE_GROUND(NPD_PROFILE)                                   &
!             GROUND SOURCE FUNCTION
     &   , K_NULL(NPD_PROFILE, NPD_LAYER)                               &
!             NULL VECTOR FOR CALL TO SUBROUTINE
     &   , DUMMY_KE(NPD_PROFILE, NPD_LAYER)
!             DUMMY ARRAY (NOT USED)
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     MONOCHROMATIC_FLUX
!
!
!
!     SET THE APPROPRIATE TOTAL UPWARD AND DOWNWARD FLUXES
!     AT THE BOUNDARIES.
!
      IF (ISOLIR == IP_SOLAR) THEN
!        VISIBLE REGION.
         DO L=1, N_PROFILE
            SOURCE_GROUND(L)=0.0E+00
            FLUX_INC_DOWN(L)=SOLAR_FLUX(L)
            FLUX_INC_DIRECT(L)=SOLAR_FLUX(L)
         ENDDO
      ELSEIF (ISOLIR == IP_INFRA_RED) THEN
!        INFRA-RED REGION.
         DO L=1, N_PROFILE
            FLUX_INC_DIRECT(L)=0.0E+00
            FLUX_DIRECT_BAND(L, N_LAYER)=0.0E+00
            FLUX_INC_DOWN(L)=-PLANCK_SOURCE_TOP(L)
            SOURCE_GROUND(L)=THERMAL_GROUND_BAND(L)                     &
     &         +(ALBEDO_SURFACE_DIFF(L)-1.0E+00)                        &
     &         *PLANCK_SOURCE_BOTTOM(L)
         ENDDO
         IF (L_CLEAR) THEN
            DO L=1, N_PROFILE
               FLUX_DIRECT_CLEAR_BAND(L, N_LAYER)=0.0E+00
            ENDDO
         ENDIF
      ENDIF
!
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            K_NULL(L, I)=0.0E+00
         ENDDO
      ENDDO
!
!
! DEPENDS ON: monochromatic_flux
      CALL MONOCHROMATIC_FLUX(IERR                                      &
!                       Atmospheric Properties
     &   , N_PROFILE, N_LAYER, D_MASS                                   &
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT         &
     &   , L_RESCALE, N_ORDER_GAUSS                                     &
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_BAND                                        &
!                       Options for Solver
     &   , I_SOLVER, L_NET, N_AUGMENT                                   &
!                       Gaseous Propreties
     &   , K_NULL                                                       &
!                       Options for Equivalent Extinction
     &   , .FALSE., DUMMY_KE                                            &
!                       Spectral Region
     &   , ISOLIR                                                       &
!                       Infra-red Properties
     &   , DIFF_PLANCK_BAND                                             &
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                         &
!                       Conditions at TOA
     &   , SEC_0, FLUX_INC_DIRECT, FLUX_INC_DOWN                        &
!                       Surface Properties
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND       &
     &   , THERMAL_GROUND_BAND                                          &
!                       Clear-sky Optical Properties
     &   , K_GREY_TOT_FREE, K_EXT_SCAT_FREE                             &
     &   , ASYMMETRY_FREE, FORWARD_SCATTER_FREE                         &
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
!                       Cloudy Optical Properties
     &   , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                           &
     &   , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD                       &
!                       Flxues Calculated
     &   , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                            &
!                       Flags for Clear-sky Calculations
     &   , L_CLEAR, I_SOLVER_CLEAR                                      &
!                       Clear-sky Fluxes Calculated
     &   , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND                &
!                       Planckian Function
     &   , PLANCK_SOURCE_BAND                                           &
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                           &
     &   )
!
!
!
      RETURN
      END SUBROUTINE SOLVE_BAND_WITHOUT_GAS
