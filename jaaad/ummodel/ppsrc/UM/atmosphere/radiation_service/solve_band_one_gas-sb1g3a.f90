


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the fluxes within the band with one gas.
!
! Method:
!       Monochromatic calculations are performed for each ESFT term
!       and the results are summed.
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
!       6.2             21/02/06                Updefs Added for version
!                                               control of radiation code
!                                               (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVE_BAND_ONE_GAS(IERR                                &
!                       Atmospheric Column
     &   , N_PROFILE, N_LAYER, L_LAYER, I_TOP, P, T, D_MASS             &
!                       Angular Integration
     &   , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT         &
     &   , L_RESCALE, N_ORDER_GAUSS                                     &
!                       Treatment of Scattering
     &   , I_SCATTER_METHOD_BAND                                        &
!                       Options for Solver
     &   , I_SOLVER, L_NET, N_AUGMENT                                   &
!                       Gaseous Properties
     &   , I_BAND, I_GAS                                                &
     &   , I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC                       &
     &   , K_ESFT, W_ESFT, L_WENYI, SCALE_VECTOR                        &
     &   , P_REFERENCE, T_REFERENCE                                     &
     &   , GAS_MIX_RATIO, GAS_FRAC_RESCALED                             &
     &   , L_DOPPLER, DOPPLER_CORRECTION                                &
!                       Spectral Region
     &   , ISOLIR                                                       &
!                       Solar Properties
     &   , SEC_0, SOLAR_FLUX                                            &
!                       Infra-red Properties
     &   , PLANCK_SOURCE_TOP, PLANCK_SOURCE_BOTTOM                      &
     &   , DIFF_PLANCK_BAND                                             &
     &   , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                         &
!                       Surface Properties
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, THERMAL_GROUND_BAND &
!                       Clear-sky Optical Properties
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
!                       Cloudy Optical Properties
     &   , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                           &
     &   , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD                       &
!                       Calculated Fluxes
     &   , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                            &
!                       Flags for Clear-sky Fluxes
     &   , L_CLEAR, I_SOLVER_CLEAR                                      &
!                       Clear-sky Fluxes
     &   , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND                &
!                       Planckian Function
     &   , PLANCK_SOURCE_BAND                                           &
!                       Dimensions of Arrays
     &   , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                           &
     &   , NPD_BAND, NPD_SPECIES                                        &
     &   , NPD_ESFT_TERM, NPD_SCALE_VARIABLE, NPD_SCALE_FNC             &
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
     &   , NPD_COLUMN                                                   &
!             NUMBER OF COLUMNS PER POINT
     &   , NPD_BAND                                                     &
!             MAXIMUM NUMBER OF BANDS
     &   , NPD_SPECIES                                                  &
!             MAXIMUM NUMBER OF SPECIES
     &   , NPD_ESFT_TERM                                                &
!             MAXIMUM NUMBER OF ESFT VARIABLES
     &   , NPD_SCALE_VARIABLE                                           &
!             MAXIMUM NUMBER OF SCALING VARIABLES
     &   , NPD_SCALE_FNC
!             MAXIMUM NUMBER OF SCALING FUNCTIONS
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
! ESFTSC3A defines types of esft scaling for two-stream radiation code.

      INTEGER,PARAMETER:: IP_SCALE_NULL=0 ! no scaling at all
      INTEGER,PARAMETER:: IP_SCALE_BAND=1 ! same scaling throughout band
      INTEGER,PARAMETER:: IP_SCALE_TERM=2 ! different for each esft

! ESFTSC3A end
! SPCRG3A defines flags for different portions of the spectrum in
! two-stream radiation code.
      INTEGER,PARAMETER:: IP_SOLAR=1
      INTEGER,PARAMETER:: IP_INFRA_RED=2
! SPCRG3A end
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
!                       Atmospheric Column
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , I_TOP
!             TOP OF VERTICAL GRID
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_LAYER
!             PROPERTIES GIVEN IN LAYERS
      REAL                                                              &
            !, INTENT(IN)
     &     P(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             PRESSURE
     &   , T(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             TEMPERATURE
     &   , D_MASS(NPD_PROFILE, NPD_LAYER)
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
!             USE AN EDGE CORRECTION
     &   , L_WENYI                                                      &
!             FLAG FOR WENYI PRESSURE & TEMPERATURE SCALING               
     &   , L_RESCALE
!             RESCALE OPTICAL PROPERTIES
!
!                       Treatment of Scattering
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SCATTER_METHOD_BAND
!             METHOD OF TREATING SCATTERING
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
!             SOLVE FOR NET FLUXES
!
!                       Gaseous Properties
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_BAND                                                       &
!             BAND BEING CONSIDERED
     &   , I_GAS                                                        &
!             GAS BEING CONSIDERED
     &   , I_BAND_ESFT(NPD_BAND, NPD_SPECIES)                           &
!             NUMBER OF TERMS IN BAND
     &   , I_SCALE_ESFT(NPD_BAND, NPD_SPECIES)                          &
!             TYPE OF ESFT SCALING
     &   , I_SCALE_FNC(NPD_BAND, NPD_SPECIES)
!             TYPE OF SCALING FUNCTION
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_DOPPLER(NPD_SPECIES)
!             DOPPLER BROADENING INCLUDED
      REAL                                                              &
            !, INTENT(IN)
     &     K_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)                 &
!             EXPONENTIAL ESFT TERMS
     &   , W_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)                 &
!             WEIGHTS FOR ESFT
     &   , SCALE_VECTOR(NPD_SCALE_VARIABLE, NPD_ESFT_TERM, NPD_BAND     &
     &        , NPD_SPECIES)                                            &
!             ABSORBER SCALING PARAMETERS
     &   , P_REFERENCE(NPD_SPECIES, NPD_BAND)                           &
!             REFERENCE SCALING PRESSURE
     &   , T_REFERENCE(NPD_SPECIES, NPD_BAND)                           &
!             REFERENCE SCALING TEMPERATURE
     &   , GAS_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES)        &
!             GAS MASS MIXING RATIOS
     &   , GAS_FRAC_RESCALED(NPD_PROFILE, 0: NPD_LAYER, NPD_SPECIES)    &
!             RESCALED GAS MASS FRACTIONS
     &   , DOPPLER_CORRECTION(NPD_SPECIES)
!             DOPPLER BROADENING TERMS
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
!             INCIDENT SOLAR FLUX IN BAND
!
!                       Infra-red Properties
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             USE A QUADRATIC SOURCE FUNCTION
      REAL                                                              &
            !, INTENT(IN)
     &     PLANCK_SOURCE_TOP(NPD_PROFILE)                               &
!             PLANCKIAN SOURCE AT TOP
     &   , PLANCK_SOURCE_BOTTOM(NPD_PROFILE)                            &
!             PLANCKIAN SOURCE AT BOTTOM
     &   , DIFF_PLANCK_BAND(NPD_PROFILE, NPD_LAYER)                     &
!             THERMAL SOURCE FUNCTION
     &   , DIFF_PLANCK_BAND_2(NPD_PROFILE, NPD_LAYER)
!             TWICE SECOND DIFFERENCE OF PLANCKIAN IN BAND
!
!                       Surface Properties
      REAL                                                              &
            !, INTENT(IN)
     &     ALBEDO_SURFACE_DIFF(NPD_PROFILE)                             &
!             DIFFUSE SURFACE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)                              &
!             DIRECT SURFACE ALBEDO
     &   , THERMAL_GROUND_BAND(NPD_PROFILE)
!             THERMAL SOURCE FUNCTION AT GROUND
!
!                       Clear-sky optical Properties
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
     &   , N_COLUMN(NPD_PROFILE)                                        &
!             NUMBER OF COLUMNS REQUIRED
     &   , N_REGION                                                     &
!             Number of cloudy regions
     &   , I_REGION_CLOUD(NPD_CLOUD_TYPE)
!             REGIONS IN WHICH TYPES OF CLOUDS FALL
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_COLUMN(NPD_PROFILE, NPD_LAYER, NPD_COLUMN)
!             FLAGS FOR CONTENT OF COLUMNS
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
!                       Fluxes Calculated
      REAL                                                              &
            !, INTENT(OUT)
     &     FLUX_DIRECT_BAND(NPD_PROFILE, 0: NPD_LAYER)                  &
!             DIRECT FLUX
     &   , FLUX_TOTAL_BAND(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX
!
!                       Flags for Clear-sky Calculations
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLEAR
!             CALCULATE NET CLEAR-SKY PROPERTIES
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SOLVER_CLEAR
!             CLEAR SOLVER USED
!
!                       Clear-sky Fluxes Calculated
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
     &     L
!             LOOP VARIABLE
      INTEGER                                                           &
     &     I_GAS_POINTER(NPD_SPECIES)                                   &
!             POINTER ARRAY FOR MONOCHROMATIC ESFTs
     &   , IEX
!             INDEX OF ESFT TERM
      REAL                                                              &
     &     K_ESFT_MONO(NPD_SPECIES)                                     &
!             ESFT MONOCHROMATIC EXPONENTS
     &   , K_GAS_ABS(NPD_PROFILE, NPD_LAYER)                            &
!             GASEOUS ABSORPTIVE EXTINCTION
     &   , SOURCE_GROUND(NPD_PROFILE)                                   &
!             GROUND SOURCE FUNCTION
     &   , FLUX_INC_DIRECT(NPD_PROFILE)                                 &
!             INCIDENT DIRECT FLUX
     &   , FLUX_INC_DOWN(NPD_PROFILE)                                   &
!             INCIDENT DOWNWARD FLUX
     &   , DUMMY_KE(NPD_PROFILE, NPD_LAYER)
!             DUMMY ARRAY (NOT USED)
      REAL                                                              &
     &     FLUX_DIRECT_PART(NPD_PROFILE, 0: NPD_LAYER)                  &
!             PARTIAL DIRECT FLUX
     &   , FLUX_TOTAL_PART(NPD_PROFILE, 2*NPD_LAYER+2)                  &
!             PARTIAL TOTAL FLUX
     &   , FLUX_DIRECT_CLEAR_PART(NPD_PROFILE, 0: NPD_LAYER)            &
!             PARTIAL CLEAR-SKY DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR_PART(NPD_PROFILE, 2*NPD_LAYER+2)
!             PARTIAL CLEAR-SKY TOTAL FLUX
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     SCALE_ABSORB, GAS_OPTICAL_PROPERTIES                         &
     &   , MONOCHROMATIC_FLUX, AUGMENT_FLUX
!
!
!

!     THE ESFT TERMS FOR THE FIRST GAS IN THE BAND ALONE ARE USED.
      I_GAS_POINTER(1)=I_GAS
      DO IEX=1, I_BAND_ESFT(I_BAND, I_GAS)
!
!        RESCALE FOR EACH ESFT TERM IF THAT IS REQUIRED.
         IF (I_SCALE_ESFT(I_BAND, I_GAS) == IP_SCALE_TERM) THEN
! DEPENDS ON: scale_absorb
            CALL SCALE_ABSORB(IERR, N_PROFILE, N_LAYER                  &
     &         , GAS_MIX_RATIO(1, 0, I_GAS), P, T                       &
     &         , L_LAYER, I_TOP                                         &
     &         , GAS_FRAC_RESCALED(1, 0, I_GAS)                         &
     &         , I_SCALE_FNC(I_BAND, I_GAS)                             &
     &         , P_REFERENCE(I_GAS, I_BAND)                             &
     &         , T_REFERENCE(I_GAS, I_BAND)                             &
     &         , SCALE_VECTOR(1, IEX, I_BAND, I_GAS)                    &
     &         , L_WENYI, IEX, I_BAND                                   &
     &         , L_DOPPLER(I_GAS), DOPPLER_CORRECTION(I_GAS)            &
     &         , NPD_PROFILE, NPD_LAYER, NPD_SCALE_FNC                  &
     &         , NPD_SCALE_VARIABLE                                     &
     &         )
            IF (IERR /= I_NORMAL) RETURN
         ENDIF
!
!        SET THE APPROPRIATE BOUNDARY TERMS FOR THE TOTAL
!        UPWARD AND DOWNWARD FLUXES.
!
         IF (ISOLIR == IP_SOLAR) THEN
!           VISIBLE REGION.
            DO L=1, N_PROFILE
               SOURCE_GROUND(L)=0.0E+00
               FLUX_INC_DOWN(L)=SOLAR_FLUX(L)
               FLUX_INC_DIRECT(L)=SOLAR_FLUX(L)
            ENDDO
         ELSEIF (ISOLIR == IP_INFRA_RED) THEN
!           INFRA-RED REGION.
            DO L=1, N_PROFILE
               FLUX_INC_DIRECT(L)=0.0E+00
               FLUX_INC_DOWN(L)=-PLANCK_SOURCE_TOP(L)
               SOURCE_GROUND(L)=THERMAL_GROUND_BAND(L)                  &
     &            -(1.0E+00-ALBEDO_SURFACE_DIFF(L))                     &
     &            *PLANCK_SOURCE_BOTTOM(L)
            ENDDO
         ENDIF
!
!        ASSIGN THE MONOCHROMATIC ABSORPTION COEFFICIENT.
         K_ESFT_MONO(I_GAS)=K_ESFT(IEX, I_BAND, I_GAS)
!
! DEPENDS ON: gas_optical_properties
         CALL GAS_OPTICAL_PROPERTIES(N_PROFILE, N_LAYER                 &
     &      , 1, I_GAS_POINTER, K_ESFT_MONO                             &
     &      , GAS_FRAC_RESCALED                                         &
     &      , K_GAS_ABS                                                 &
     &      , NPD_PROFILE, NPD_LAYER, NPD_SPECIES                       &
     &      )
!
!
! DEPENDS ON: monochromatic_flux
         CALL MONOCHROMATIC_FLUX(IERR                                   &
!                       Atmospheric Properties
     &      , N_PROFILE, N_LAYER, D_MASS                                &
!                       Angular Integration
     &      , I_ANGULAR_INTEGRATION, I_2STREAM, L_2_STREAM_CORRECT      &
     &      , L_RESCALE, N_ORDER_GAUSS                                  &
!                       Treatment of Scattering
     &      , I_SCATTER_METHOD_BAND                                     &
!                       Options for Solver
     &      , I_SOLVER, L_NET, N_AUGMENT                                &
!                       Gaseous Propreties
     &      , K_GAS_ABS                                                 &
!                       Options for Equivalent Extinction
     &      , .FALSE., DUMMY_KE                                         &
!                       Spectral Region
     &      , ISOLIR                                                    &
!                       Infra-red Properties
     &      , DIFF_PLANCK_BAND                                          &
     &      , L_IR_SOURCE_QUAD, DIFF_PLANCK_BAND_2                      &
!                       Conditions at TOA
     &      , SEC_0, FLUX_INC_DIRECT, FLUX_INC_DOWN                     &
!                       Surface Properties
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR, SOURCE_GROUND    &
     &      , THERMAL_GROUND_BAND                                       &
!                       Clear-sky Optical Properties
     &      , K_GREY_TOT_FREE, K_EXT_SCAT_FREE                          &
     &      , ASYMMETRY_FREE, FORWARD_SCATTER_FREE                      &
!                       Cloudy Properties
     &      , L_CLOUD, I_CLOUD                                          &
!                       Cloud Geometry
     &      , N_CLOUD_TOP                                               &
     &      , N_CLOUD_TYPE, FRAC_CLOUD                                  &
     &      , N_REGION, I_REGION_CLOUD, FRAC_REGION                     &
     &      , W_FREE, N_FREE_PROFILE, I_FREE_PROFILE                    &
     &      , W_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE                 &
     &      , CLOUD_OVERLAP                                             &
     &      , N_COLUMN, L_COLUMN, AREA_COLUMN                           &
!                       Cloudy Optical Properties
     &      , K_GREY_TOT_CLOUD, K_EXT_SCAT_CLOUD                        &
     &      , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD                    &
!                       Flxues Calculated
     &      , FLUX_DIRECT_PART, FLUX_TOTAL_PART                         &
!                       Flags for Clear-sky Calculations
     &      , L_CLEAR, I_SOLVER_CLEAR                                   &
!                       Clear-sky Fluxes Calculated
     &      , FLUX_DIRECT_CLEAR_PART, FLUX_TOTAL_CLEAR_PART             &
!                       Planckian Function
     &      , PLANCK_SOURCE_BAND                                        &
!                       Dimensions of Arrays
     &      , NPD_PROFILE, NPD_LAYER, NPD_COLUMN                        &
     &      )
         IF (IERR /= I_NORMAL) RETURN
!
!        INCREMENT THE FLUXES WITHIN THE BAND.
! DEPENDS ON: augment_flux
         CALL AUGMENT_FLUX(N_PROFILE, N_LAYER, N_AUGMENT                &
     &      , ISOLIR, L_CLEAR                                           &
     &      , W_ESFT(IEX, I_BAND,  I_GAS)                               &
     &      , FLUX_DIRECT_BAND, FLUX_TOTAL_BAND                         &
     &      , FLUX_DIRECT_PART, FLUX_TOTAL_PART                         &
     &      , FLUX_DIRECT_CLEAR_BAND, FLUX_TOTAL_CLEAR_BAND             &
     &      , FLUX_DIRECT_CLEAR_PART, FLUX_TOTAL_CLEAR_PART             &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
      ENDDO
!
!
      RETURN
      END SUBROUTINE SOLVE_BAND_ONE_GAS
