
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate grey extinctions.
!
! Method:
!       For each activated optical process, excluding gaseous
!       absorption, increments are calculated for the total and
!       scattering extinctions, and the products of the asymmetry
!       factor and the forward scattering factor in clear and
!       cloudy regions. These increments are summed, and the grey
!       total and scattering extinctions and the asymmetry and forward
!       scattering factors are thus calculated.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE GREY_EXTINCTION(IERR                                   &
     &   , N_PROFILE, N_LAYER, L_LAYER, P, T, DENSITY                   &
     &   , L_RESCALE                                                    &
     &   , L_RAYLEIGH, RAYLEIGH_COEFF                                   &
     &   , L_CONTINUUM, N_CONTINUUM, I_CONTINUUM_POINTER, K_CONTINUUM   &
     &   , AMOUNT_CONTINUUM                                             &
     &   , L_AEROSOL, N_AEROSOL, AEROSOL_MIX_RATIO                      &
     &   , I_AEROSOL_PARAMETRIZATION                                    &
     &   , I_HUMIDITY_POINTER, HUMIDITIES, DELTA_HUMIDITY               &
     &   , MEAN_REL_HUMIDITY                                            &
     &   , AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_ASYMMETRY    &
     &   , L_CLOUD, N_CLOUD_PROFILE, I_CLOUD_PROFILE, N_CLOUD_TOP       &
     &   , L_CLOUD_LAYER, I_CLOUD                                       &
     &   , N_CONDENSED, L_CLOUD_CMP, I_PHASE_CMP                        &
     &   , I_CONDENSED_PARAM, CONDENSED_PARAM_LIST                      &
     &   , CONDENSED_MIX_RATIO, CONDENSED_DIM_CHAR                      &
     &   , N_CLOUD_TYPE, I_CLOUD_TYPE                                   &
     &   , K_EXT_TOT_FREE, K_EXT_SCAT_FREE, ASYMMETRY_FREE              &
     &   , FORWARD_SCATTER_FREE                                         &
     &   , K_EXT_TOT_CLOUD, K_EXT_SCAT_CLOUD                            &
     &   , ASYMMETRY_CLOUD, FORWARD_SCATTER_CLOUD                       &
     &   , FRAC_CLOUD, l_pc2                                            &
     &   , L_CLOUD_EXTINCTION, CLOUD_EXTINCTION                         &
     &   , L_CLOUD_ABSORPTIVITY, CLOUD_ABSORPTIVITY                     &
     &   , L_LS_CLOUD_EXTINCTION, LS_CLOUD_EXTINCTION                   &
     &   , L_LS_CLOUD_ABSORPTIVITY, LS_CLOUD_ABSORPTIVITY               &
     &   , L_CNV_CLOUD_EXTINCTION, CNV_CLOUD_EXTINCTION                 &
     &   , L_CNV_CLOUD_ABSORPTIVITY, CNV_CLOUD_ABSORPTIVITY             &
     &   , NPD_PROFILE, NPD_LAYER, NPD_CONTINUUM                        &
     &   , NPD_AEROSOL_SPECIES, NPD_HUMIDITIES                          &
     &   , NPD_CLOUD_PARAMETER                                          &
     &   )
!
!
!
      IMPLICIT NONE
!
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_AEROSOL_SPECIES                                          &
!             MAXIMUM NUMBER OF AEROSOLS
     &   , NPD_HUMIDITIES                                               &
!             MAXIMUM NUMBER OF HUMIDITIES
     &   , NPD_CONTINUUM                                                &
!             MAXIMUM NUMBER OF CONTINUA
     &   , NPD_CLOUD_PARAMETER
!             MAXIMUM NUMBER OF CLOUD PARAMETERS
!
!     INCLUDE COMDECKS
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
! AERPRM3A aerosol parametrizations for two-stream radiation code.

      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_DRY=1
      INTEGER,PARAMETER:: IP_AEROSOL_PARAM_MOIST=2
      INTEGER,PARAMETER:: IP_AEROSOL_UNPARAMETRIZED=3 ! Observational

! AERPRM3A end
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
! PHASE3A defines indices for phases in two-stream radiation code.
      INTEGER,PARAMETER:: IP_PHASE_WATER = 1
      INTEGER,PARAMETER:: IP_PHASE_ICE   = 2
! PHASE3A end
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
! CLDTYP3A defines cloud types for TWO-STREAM RADIATION CODE.

      INTEGER,PARAMETER:: IP_CLOUD_TYPE_HOMOGEN=1 ! water and ice
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_WATER=1   ! Water only
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_ICE=2     ! Ice only

      ! mixed-phase stratiform cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_STRAT=1

      ! mixed-phase convective cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_CONV=2

      INTEGER,PARAMETER:: IP_CLOUD_TYPE_SW=1 ! stratiform water cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_SI=2 ! stratiform ice cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_CW=3 ! convective water cloud
      INTEGER,PARAMETER:: IP_CLOUD_TYPE_CI=4 ! convective ice cloud

! CLDTYP3A end
!
!
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
!
!
!     BASIC ATMOSPHERIC PROPERTIES:
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_LAYER
!             VARIABLES GIVEN IN LAYERS
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
!
      REAL                                                              &
                !, INTENT(IN)
     &     P(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             PRESSURE
     &   , T(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             TEMPERATURE
     &   , DENSITY(NPD_PROFILE, 0: NPD_LAYER)
!             DENSITY AT LEVELS
!
!
!     OPTICAL SWITCHES:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_RESCALE
!             DELTA-RESCALING REQUIRED
!
!
!     RAYLEIGH SCATTERING:
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_RAYLEIGH
!             RAYLEIGH SCATTERING ACTIVATED
!
      REAL                                                              &
                !, INTENT(IN)
     &     RAYLEIGH_COEFF
!             RAYLEIGH COEFFICIENT
!
!
!     CONTINUUM PROCESSES:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CONTINUUM
!             CONTINUUM ABSORPTION ACTIVATED
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CONTINUUM                                                  &
!             NUMBER OF CONTINUA
     &   , I_CONTINUUM_POINTER(NPD_CONTINUUM)
!             POINTERS TO ACTIVE CONTINUA
!
      REAL                                                              &
                !, INTENT(IN)
     &     K_CONTINUUM(NPD_CONTINUUM)                                   &
!             CONTINUUM EXTINCTION
     &   , AMOUNT_CONTINUUM(NPD_PROFILE, 0: NPD_LAYER, NPD_CONTINUUM)
!             AMOUNTS FOR CONTINUA
!
!
!     PROPERTIES OF AEROSOLS:
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_AEROSOL
!             AEROSOLS ACTIVATED
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_AEROSOL                                                    &
!             NUMBER OF AEROSOL SPECIES
     &   , I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL_SPECIES)               &
!             PARAMETRIZATIONS OF AEROSOLS
     &   , I_HUMIDITY_POINTER(NPD_PROFILE,  NPD_LAYER)
!             POINTER TO AEROSOL LOOK-UP TABLE
!
      REAL                                                              &
                !, INTENT(IN)
     &     AEROSOL_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER                  &
     &        , NPD_AEROSOL_SPECIES)                                    &
!             NUMBER DENSTY OF AEROSOLS
     &   , AEROSOL_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)      &
!             AEROSOL ABSORPTION IN BAND/MIX RT.
     &   , AEROSOL_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)      &
!             AEROSOL SCATTERING IN BAND/MIX RT.
     &   , AEROSOL_ASYMMETRY(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)       &
!             AEROSOL ASYMMETRY IN BAND
     &   , HUMIDITIES(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)              &
!             ARRAY OF HUMIDITIES
     &   , DELTA_HUMIDITY                                               &
!             INCREMENT IN HUMIDITY
     &   , MEAN_REL_HUMIDITY(NPD_PROFILE, NPD_LAYER)
!             RELATIVE HUMIDITY. MAY BE THE CLEAR-SKY MEAN OR
!             GRID-BOX MEAN, DEPENDING ON CALCULATIONS MADE IN
!             FLUX_CALC
!
!
!
!     PROPERTIES OF CLOUDS:
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLOUD                                                      &
!             CLOUDS ACTIVATED
     &   , L_pc2
!             Use PC2 cloud scheme
!
!     GEOMETRY OF CLOUDS:
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLOUD_LAYER
!             CLOUD VARIABLES GIVEN IN LAYERS
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CLOUD_TOP                                                  &
!             TOPMOST CLOUDY LAYER
     &   , I_CLOUD                                                      &
!             CLOUD SCHEME TO BE USED
     &   , N_CLOUD_TYPE                                                 &
!             NUMBER OF TYPES OF CLOUDS
     &   , N_CLOUD_PROFILE(NPD_LAYER)                                   &
!             NUMBER OF CLOUDY PROFILES IN EACH LAYER
     &   , I_CLOUD_PROFILE(NPD_PROFILE, NPD_LAYER)                      &
!             PROFILES CONTAINING CLOUDS
     &   , I_CLOUD_TYPE(NPD_CLOUD_COMPONENT)
!             TYPES OF CLOUD TO WHICH EACH COMPONENT CONTRIBUTES
!
!     MICROPHYSICAL QUANTITIES:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CONDENSED                                                  &
!             NUMBER OF CONDENSED COMPONENTS
     &   , I_PHASE_CMP(NPD_CLOUD_COMPONENT)                             &
!             PHASES OF CLOUDY COMPONENTS
     &   , I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             PARAMETRIZATION SCHEMES FOR CLOUDY COMPONENTS
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLOUD_CMP(NPD_CLOUD_COMPONENT)
!             FLAGS TO ACTIVATE CLOUDY COMPONENTS
!
      REAL                                                              &
                !, INTENT(IN)
     &     CONDENSED_PARAM_LIST(NPD_CLOUD_PARAMETER                     &
     &        , NPD_CLOUD_COMPONENT)                                    &
!             COEFFICIENTS IN PARAMETRIZATION SCHEMES
     &   , CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER                &
     &        , NPD_CLOUD_COMPONENT)                                    &
!             MIXING RATIOS OF CLOUDY COMPONENTS
     &   , CONDENSED_DIM_CHAR(NPD_PROFILE, 0: NPD_LAYER                 &
     &      , NPD_CLOUD_COMPONENT)
!             EFFECTIVE RADII OF CLOUDY COMPONENTS

!
!     VARIABLES REQUIRED FOR EXTRA DIAGNOSTIC CALCULATIONS.
!
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLOUD_EXTINCTION                                           &
!             FLAG FOR EXPLICIT CALCULATION OF EXTINCTION
     &   , L_CLOUD_ABSORPTIVITY                                         &
!             FLAG FOR EXPLICIT CALCULATION OF ABSORPTIVITY
     &   , L_LS_CLOUD_EXTINCTION                                        &
!             FLAG FOR EXPLICIT CALCULATION OF EXTINCTION
     &   , L_LS_CLOUD_ABSORPTIVITY                                      &
!             FLAG FOR EXPLICIT CALCULATION OF ABSORPTIVITY
     &   , L_CNV_CLOUD_EXTINCTION                                       &
!             FLAG FOR EXPLICIT CALCULATION OF EXTINCTION
     &   , L_CNV_CLOUD_ABSORPTIVITY
!             FLAG FOR EXPLICIT CALCULATION OF ABSORPTIVITY
!
      REAL                                                              &
                !, INTENT(IN)
     &     FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             FRACTIONS OF EACH TYPE OF CLOUD
!
!
!
!     CALCULATED OPTICAL PROPETIES:
!
      REAL                                                              &
                !, INTENT(OUT)
     &     K_EXT_SCAT_FREE(NPD_PROFILE, NPD_LAYER)                      &
!             FREE SCATTERING EXTINCTION
     &   , K_EXT_TOT_FREE(NPD_PROFILE, NPD_LAYER)                       &
!             TOTAL FREE EXTINCTION
     &   , ASYMMETRY_FREE(NPD_PROFILE, NPD_LAYER)                       &
!             FREE ASYMMETRIES
     &   , FORWARD_SCATTER_FREE(NPD_PROFILE, NPD_LAYER)                 &
!             FREE FORWARD SCATTERING
     &   , K_EXT_SCAT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)     &
!             CLOUDY SCATTERING EXTINCTION
     &   , K_EXT_TOT_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)      &
!             TOTAL CLOUDY EXTINCTION
     &   , ASYMMETRY_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)      &
!             CLOUDY ASYMMETRIES
     &   , FORWARD_SCATTER_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)
!             CLOUDY FORWARD SCATTERING
!
      REAL                                                              &
                !, INTENT(OUT)
     &     CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                     &
!             MEAN CLOUD EXTINCTION
     &   , CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)                   &
!             MEAN CLOUD ABSORPTIVITY
     &   , LS_CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                  &
!             MEAN CLOUD EXTINCTION
     &   , LS_CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)                &
!             MEAN CLOUD ABSORPTIVITY
     &   , CNV_CLOUD_EXTINCTION(NPD_PROFILE, NPD_LAYER)                 &
!             MEAN CLOUD EXTINCTION
     &   , CNV_CLOUD_ABSORPTIVITY(NPD_PROFILE, NPD_LAYER)
!             MEAN CLOUD ABSORPTIVITY
!
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     I_CONTINUUM                                                  &
!             TEMPORARY CONTINUUM INDEX
     &   , I_POINTER                                                    &
!             TEMPORARY POINTER
     &   , L                                                            &
!             LOOP VARIABLE
     &   , LL                                                           &
!             LOOP VARIABLE
     &   , I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
!             LOOP VARIABLE
     &   , N_INDEX                                                      &
!             NUMBER OF INDICES SATISFYING TEST
     &   , INDEX(NPD_PROFILE)
!             INDICES OF TESTED POINTS
!
!     TEMPORARY OPTICAL PROPERTIES:
!
      REAL                                                              &
     &     K_EXT_SCAT_CLOUD_COMP(NPD_PROFILE, NPD_LAYER)                &
!             SCATTERING EXTINCTION OF CLOUDY COMPONENT
     &   , K_EXT_TOT_CLOUD_COMP(NPD_PROFILE, NPD_LAYER)                 &
!             TOTAL EXTINCTION OF CLOUDY COMPONENT
     &   , ASYMMETRY_CLOUD_COMP(NPD_PROFILE, NPD_LAYER)                 &
!             ASYMMETRIES OF CLOUDY COMPONENT
     &   , FORWARD_SCATTER_CLOUD_COMP(NPD_PROFILE, NPD_LAYER)           &
!             FORWARD SCATTERING OF CLOUDY COMPONENT
     &   , K_SCATTER(NPD_PROFILE)                                       &
!             SCATTERING VARIABLE
     &   , ASYMMETRY_PROCESS(NPD_PROFILE)
!             ASYMMETRY FACTOR FOR CURRENT PROC.
!
!
      REAL                                                              &
     &     WEIGHT_UPPER                                                 &
!             UPPER WEIGHT FOR INTERPOLATION
     &   , WEIGHT_LOWER
!             LOWER WEIGHT FOR INTERPOLATION
!
!     Temporary variables for the divisions
      REAL                                                              &
     &  TMP_INV,TMP_INV1

!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     OPT_PROP_WATER_CLOUD, OPT_PROP_ICE_CLOUD
!
!     CRAY DIRECTIVES FOR THE WHOLE ROUTINE:
!     POINTS ARE NOT REPEATED IN THE INDEXING ARRAY, SO IT IS SAFE
!     TO VECTORIZE OVER INDIRECTLY ADDRESSED ARRAYS.
!fpp$ NODEPCHK R
!
!
!
!     INITIALIZE THE EXTINCTION COEFFICIENTS AND THE ASYMMETRY PRODUCT.

      IF(L_RESCALE) THEN

!     Forward scattering is required only in the visible where
!     delta-rescaling is performed.

         IF(.NOT.L_RAYLEIGH) THEN

            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  K_EXT_TOT_FREE(L, I)=0.0E+00
                  K_EXT_SCAT_FREE(L, I)=0.0E+00
                  ASYMMETRY_FREE(L, I)=0.0E+00
                  FORWARD_SCATTER_FREE(L, I)=0.0E+00
               ENDDO
            ENDDO

         ELSE IF(L_RAYLEIGH) THEN

!        INCLUDE RAYLEIGH SCATTERING.

            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  K_EXT_TOT_FREE(L, I)=0.0E+00
                  K_EXT_SCAT_FREE(L, I)=RAYLEIGH_COEFF
                  ASYMMETRY_FREE(L, I)=0.0E+00
                  FORWARD_SCATTER_FREE(L, I)=0.0E+00
               ENDDO
            ENDDO
         END IF

      ELSE IF(.NOT.L_RESCALE) THEN

         IF(.NOT.L_RAYLEIGH) THEN

            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  K_EXT_TOT_FREE(L, I)=0.0E+00
                  K_EXT_SCAT_FREE(L, I)=0.0E+00
                  ASYMMETRY_FREE(L, I)=0.0E+00
               ENDDO
            ENDDO

         ELSE IF(L_RAYLEIGH) THEN

!        INCLUDE RAYLEIGH SCATTERING.

            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  K_EXT_TOT_FREE(L, I)=0.0E+00
                  K_EXT_SCAT_FREE(L, I)=RAYLEIGH_COEFF
                  ASYMMETRY_FREE(L, I)=0.0E+00
               ENDDO
            ENDDO

         END IF

      END IF

!
      IF (L_AEROSOL) THEN
!        INCLUDE THE EFFECTS OF AEROSOL.
         DO J=1, N_AEROSOL
            IF (I_AEROSOL_PARAMETRIZATION(J)                            &
     &          == IP_AEROSOL_PARAM_DRY) THEN
               DO I=1, N_LAYER
                  DO L=1, N_PROFILE
                     K_EXT_TOT_FREE(L, I)=K_EXT_TOT_FREE(L, I)          &
     &                  +AEROSOL_MIX_RATIO(L, I, J)                     &
     &                  *AEROSOL_ABSORPTION(1, J)
                     K_SCATTER(L)=AEROSOL_MIX_RATIO(L, I, J)            &
     &                  *AEROSOL_SCATTERING(1, J)
                     K_EXT_SCAT_FREE(L, I)=K_EXT_SCAT_FREE(L, I)        &
     &                  +K_SCATTER(L)
                     ASYMMETRY_FREE(L, I)=ASYMMETRY_FREE(L, I)          &
     &                  +K_SCATTER(L)*AEROSOL_ASYMMETRY(1, J)
                  ENDDO
                  IF (L_RESCALE) THEN
!                    THIS BLOCK IS PLACED WITHIN THE LOOP OVER I TO SAVE
!                    STORAGE. THE COST OF RE-EXECUTING THE TEST IS QUITE
!                    SMALL.
                     DO L=1, N_PROFILE
                        FORWARD_SCATTER_FREE(L, I)                      &
     &                     =FORWARD_SCATTER_FREE(L, I)+K_SCATTER(L)     &
     &                     *(AEROSOL_ASYMMETRY(1, J))**2
                     ENDDO
                  ENDIF
               ENDDO
            ELSE IF (I_AEROSOL_PARAMETRIZATION(J)                       &
     &          == IP_AEROSOL_PARAM_MOIST) THEN
               DO I=1, N_LAYER
                  DO L=1, N_PROFILE
                     I_POINTER=I_HUMIDITY_POINTER(L, I)
                     WEIGHT_UPPER=(MEAN_REL_HUMIDITY(L, I)              &
     &                 -HUMIDITIES(I_POINTER, J))                       &
     &                 /DELTA_HUMIDITY
                     WEIGHT_LOWER=1.0E+00-WEIGHT_UPPER
                     K_EXT_TOT_FREE(L, I)=K_EXT_TOT_FREE(L, I)          &
     &                  +AEROSOL_MIX_RATIO(L, I, J)                     &
     &                  *(AEROSOL_ABSORPTION(I_POINTER, J)              &
     &                  *WEIGHT_LOWER+WEIGHT_UPPER                      &
     &                  *AEROSOL_ABSORPTION(I_POINTER+1, J))
                     K_SCATTER(L)=                                      &
     &                  AEROSOL_MIX_RATIO(L, I, J)                      &
     &                  *(AEROSOL_SCATTERING(I_POINTER, J)              &
     &                  *WEIGHT_LOWER+WEIGHT_UPPER                      &
     &                  *AEROSOL_SCATTERING(I_POINTER+1, J))
                     K_EXT_SCAT_FREE(L, I)=K_EXT_SCAT_FREE(L, I)        &
     &                  +K_SCATTER(L)
                     ASYMMETRY_PROCESS(L)=                              &
     &                  AEROSOL_ASYMMETRY(I_POINTER, J)                 &
     &                  *WEIGHT_LOWER+WEIGHT_UPPER                      &
     &                  *AEROSOL_ASYMMETRY(I_POINTER+1, J)
                     ASYMMETRY_FREE(L, I)=ASYMMETRY_FREE(L, I)          &
     &                  +K_SCATTER(L)*ASYMMETRY_PROCESS(L)
                  ENDDO
                  IF (L_RESCALE) THEN
                     DO L=1, N_PROFILE
                        FORWARD_SCATTER_FREE(L, I)                      &
     &                     =FORWARD_SCATTER_FREE(L, I)+K_SCATTER(L)     &
     &                     *(ASYMMETRY_PROCESS(L))**2
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               WRITE(IU_ERR, '(/A, I3, A)')                             &
     &            '*** ERROR : I_AEROSOL_PARAMETRIZATION FOR SPECIES '  &
     &            , J, ' HAS BEEN SET TO AN ILLEGAL VALUE.'
               IERR=I_ERR_FATAL
               RETURN

            ENDIF
         ENDDO
      ENDIF
!
      IF (L_CONTINUUM) THEN
!        INCLUDE CONTINUUM ABSORPTION.
         DO J=1, N_CONTINUUM
            I_CONTINUUM=I_CONTINUUM_POINTER(J)
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  K_EXT_TOT_FREE(L, I)=K_EXT_TOT_FREE(L, I)             &
     &               +K_CONTINUUM(I_CONTINUUM)                          &
     &               *AMOUNT_CONTINUUM(L, I, I_CONTINUUM)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
!
!
!     ADD THE SCATTERING ON TO THE TOTAL EXTINCTION. THE FINAL FREE
!     ASYMMETRY IS NOT CALCULATED HERE SINCE THE PRODUCT OF ASYMMETRY
!     AND SCATTERING IS ALSO NEEDED TO CALCULATE THE CLOUDY ASYMMETRY.
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            K_EXT_TOT_FREE(L, I)=K_EXT_TOT_FREE(L, I)                   &
     &         +K_EXT_SCAT_FREE(L, I)
         ENDDO
      ENDDO
!
!
!     IF THERE ARE NO CLOUDS CALCULATE THE FINAL OPTICAL PROPERTIES
!     AND RETURN TO THE CALLING ROUTINE.
!
      IF (.NOT.L_CLOUD) THEN

         IF(.NOT.L_RESCALE) THEN
            DO I=1, N_LAYER
               DO L=1, N_PROFILE
                  IF (K_EXT_SCAT_FREE(L, I) >                           &
     &               TINY(K_EXT_SCAT_FREE)) THEN
                     ASYMMETRY_FREE(L, I)=ASYMMETRY_FREE(L, I)          &
     &                    /K_EXT_SCAT_FREE(L, I)
                  ENDIF
               ENDDO
            ENDDO
         ELSE IF (L_RESCALE) THEN
               DO I=1, N_LAYER
                  DO L=1, N_PROFILE
                     IF (K_EXT_SCAT_FREE(L, I) >                        &
     &                  TINY(K_EXT_SCAT_FREE)) THEN
                        TMP_INV=1.0/K_EXT_SCAT_FREE(L,I)
                        FORWARD_SCATTER_FREE(L, I)                      &
     &                       =FORWARD_SCATTER_FREE(L, I)                &
     &                       *TMP_INV
                        ASYMMETRY_FREE(L, I)=ASYMMETRY_FREE(L, I)       &
     &                       *TMP_INV
                     ENDIF
                  ENDDO
               ENDDO

         ENDIF

         RETURN

      ENDIF


!
!     ADDITION OF CLOUDY PROPERTIES:
!
!
!     ADD IN BACKGROUND CONTIBUTIONS:
!
!
!     ALL THE PROCESSES OCCURRING OUTSIDE CLOUDS ALSO OCCUR
!     WITHIN THEM.
      DO K=1, N_CLOUD_TYPE
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               K_EXT_TOT_CLOUD(L, I, K)=K_EXT_TOT_FREE(L, I)
               K_EXT_SCAT_CLOUD(L, I, K)=K_EXT_SCAT_FREE(L, I)
               ASYMMETRY_CLOUD(L, I, K)=ASYMMETRY_FREE(L, I)
               FORWARD_SCATTER_CLOUD(L, I, K)                           &
     &            =FORWARD_SCATTER_FREE(L, I)
            ENDDO
         ENDDO
      ENDDO
!
!
!
!     INITIALIZE ARRAYS FOR DIAGNOSTIC USE.
!        N.B. This initialization was required at vn4.5 - it may not
!        be required at vn5.2 and later, but is included for now
!        for safety.
!
      IF (L_CLOUD_EXTINCTION) THEN
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               CLOUD_EXTINCTION(L, I)=0.0E+00
            ENDDO
         ENDDO
      ENDIF
!
      IF (L_CLOUD_ABSORPTIVITY) THEN
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               CLOUD_ABSORPTIVITY(L, I)=0.0E+00
            ENDDO
         ENDDO
      ENDIF
!
      IF (L_LS_CLOUD_EXTINCTION) THEN
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               LS_CLOUD_EXTINCTION(L, I)=0.0E+00
            ENDDO
         ENDDO
      ENDIF
!
      IF (L_LS_CLOUD_ABSORPTIVITY) THEN
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               LS_CLOUD_ABSORPTIVITY(L, I)=0.0E+00
            ENDDO
         ENDDO
      ENDIF
!
      IF (L_CNV_CLOUD_EXTINCTION) THEN
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               CNV_CLOUD_EXTINCTION(L, I)=0.0E+00
            ENDDO
         ENDDO
      ENDIF
!
      IF (L_CNV_CLOUD_ABSORPTIVITY) THEN
         DO I=1, N_LAYER
            DO L=1, N_PROFILE
               CNV_CLOUD_ABSORPTIVITY(L, I)=0.0E+00
            ENDDO
         ENDDO
      ENDIF
!
!
!     ADD ON THE TERMS REPRESENTING PROCESSES WITHIN CLOUDS.
!
!     LOOP OVER THE CONDENSED COMPONENTS, CALCULATING THEIR OPTICAL
!     PROPERTIES AND THEN ASSIGN THEM TO THE ARRAYS FOR THE TYPES OF
!     CLOUD.
!
      DO K=1, N_CONDENSED
!
!        FLAGS FOR DEALING WITH COMPONENTS WERE SET IN THE SUBROUTINE
!        SET_CLOUD_POINTER. WE NOW DETERMINE WHETHER THE COMPONENT IS
!        TO BE INCLUDED AND CALCULATE ITS OPTICAL PROPERTIES ACCORDING
!        TO THE PHASE OF THE COMPONENT. THESE CONTRIBUTIONS ARE ADDED
!        TO THE ARRAYS FOR THE SELECTED TYPE OF CLOUD.
!
         IF (L_CLOUD_CMP(K)) THEN
!
            IF (I_PHASE_CMP(K) == IP_PHASE_WATER) THEN
!
!              INCLUDE SCATTERING BY WATER DROPLETS.
!
! DEPENDS ON: opt_prop_water_cloud
               CALL OPT_PROP_WATER_CLOUD(IERR                           &
     &            , N_PROFILE, N_LAYER, N_CLOUD_TOP                     &
     &            , N_CLOUD_PROFILE, I_CLOUD_PROFILE                    &
     &            , L_RESCALE, L_LAYER, L_CLOUD_LAYER                   &
     &            , I_CONDENSED_PARAM(K), CONDENSED_PARAM_LIST(1, K)    &
     &            , CONDENSED_MIX_RATIO(1, 0, K)                        &
     &            , CONDENSED_DIM_CHAR(1, 0, K)                         &
     &            , K_EXT_TOT_CLOUD_COMP, K_EXT_SCAT_CLOUD_COMP         &
     &            , ASYMMETRY_CLOUD_COMP, FORWARD_SCATTER_CLOUD_COMP    &
     &            , NPD_PROFILE, NPD_LAYER                              &
     &            , NPD_CLOUD_PARAMETER                                 &
     &            )
!
            ELSE IF (I_PHASE_CMP(K) == IP_PHASE_ICE) THEN
!
!              INCLUDE SCATTERING BY ICE CRYSTALS.
!
! DEPENDS ON: opt_prop_ice_cloud
               CALL OPT_PROP_ICE_CLOUD(IERR                             &
     &            , N_PROFILE, N_LAYER, N_CLOUD_TOP                     &
     &            , N_CLOUD_PROFILE, I_CLOUD_PROFILE                    &
     &            , L_RESCALE, L_LAYER, L_CLOUD_LAYER                   &
     &            , I_CONDENSED_PARAM(K), CONDENSED_PARAM_LIST(1, K)    &
     &            , CONDENSED_MIX_RATIO(1, 0, K)                        &
     &            , CONDENSED_DIM_CHAR(1, 0, K)                         &
     &            , T, DENSITY                                          &
     &            , K_EXT_TOT_CLOUD_COMP, K_EXT_SCAT_CLOUD_COMP         &
     &            , ASYMMETRY_CLOUD_COMP, FORWARD_SCATTER_CLOUD_COMP    &
     &            , NPD_PROFILE, NPD_LAYER                              &
     &            , NPD_CLOUD_PARAMETER                                 &
     &            )
!
            ENDIF
!
!
!
!           INCREMENT THE ARRAYS OF OPTICAL PROPERTIES.
!
!
            DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
               DO LL=1, N_CLOUD_PROFILE(I)
                  L=I_CLOUD_PROFILE(LL, I)
                  K_EXT_TOT_CLOUD(L, I, I_CLOUD_TYPE(K))                &
     &               =K_EXT_TOT_CLOUD(L, I, I_CLOUD_TYPE(K))            &
     &               +K_EXT_TOT_CLOUD_COMP(L, I)
                  K_EXT_SCAT_CLOUD(L, I, I_CLOUD_TYPE(K))               &
     &               =K_EXT_SCAT_CLOUD(L, I, I_CLOUD_TYPE(K))           &
     &               +K_EXT_SCAT_CLOUD_COMP(L, I)
                  ASYMMETRY_CLOUD(L, I, I_CLOUD_TYPE(K))                &
     &               =ASYMMETRY_CLOUD(L, I, I_CLOUD_TYPE(K))            &
     &               +ASYMMETRY_CLOUD_COMP(L, I)
               ENDDO
            ENDDO
            IF (L_RESCALE) THEN
               DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
                  DO LL=1, N_CLOUD_PROFILE(I)
                     L=I_CLOUD_PROFILE(LL, I)
                     FORWARD_SCATTER_CLOUD(L, I, I_CLOUD_TYPE(K))       &
     &                  =FORWARD_SCATTER_CLOUD(L, I, I_CLOUD_TYPE(K))   &
     &                  +FORWARD_SCATTER_CLOUD_COMP(L, I)
                  ENDDO
               ENDDO
            ENDIF
!
!
!           EXTRA CALCULATIONS FOR DIAGNOSTICS.
!
            IF (L_CLOUD_EXTINCTION) THEN
               DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
                  DO LL=1, N_CLOUD_PROFILE(I)
                     L=I_CLOUD_PROFILE(LL, I)
                     CLOUD_EXTINCTION(L, I)                             &
     &                  =CLOUD_EXTINCTION(L, I)                         &
     &                  +K_EXT_TOT_CLOUD_COMP(L, I)                     &
     &                  *FRAC_CLOUD(L, I, I_CLOUD_TYPE(K))
                  ENDDO
               ENDDO
            ENDIF
!
!
            IF (L_CLOUD_ABSORPTIVITY) THEN
               DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
                  DO LL=1, N_CLOUD_PROFILE(I)
                     L=I_CLOUD_PROFILE(LL, I)
                     CLOUD_ABSORPTIVITY(L, I)                           &
     &                  =CLOUD_ABSORPTIVITY(L, I)                       &
     &                  +(K_EXT_TOT_CLOUD_COMP(L, I)                    &
     &                  -K_EXT_SCAT_CLOUD_COMP(L, I))                   &
     &                  *FRAC_CLOUD(L, I, I_CLOUD_TYPE(K))
                  ENDDO
               ENDDO
            ENDIF
!
            IF ((I_CLOUD_TYPE(K) == IP_CLOUD_TYPE_SW).OR.               &
     &                (I_CLOUD_TYPE(K) == IP_CLOUD_TYPE_SI)) THEN
!
            IF (L_LS_CLOUD_EXTINCTION) THEN
               DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
                  DO LL=1, N_CLOUD_PROFILE(I)
                     L=I_CLOUD_PROFILE(LL, I)
                     LS_CLOUD_EXTINCTION(L, I)                          &
     &                  =LS_CLOUD_EXTINCTION(L, I)                      &
     &                  +K_EXT_TOT_CLOUD_COMP(L, I)                     &
     &                  *FRAC_CLOUD(L, I, I_CLOUD_TYPE(K))
                  ENDDO
               ENDDO
            ENDIF
!
!
            IF (L_LS_CLOUD_ABSORPTIVITY) THEN
               DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
                  DO LL=1, N_CLOUD_PROFILE(I)
                     L=I_CLOUD_PROFILE(LL, I)
                     LS_CLOUD_ABSORPTIVITY(L, I)                        &
     &                  =LS_CLOUD_ABSORPTIVITY(L, I)                    &
     &                  +(K_EXT_TOT_CLOUD_COMP(L, I)                    &
     &                  -K_EXT_SCAT_CLOUD_COMP(L, I))                   &
     &                  *FRAC_CLOUD(L, I, I_CLOUD_TYPE(K))
                  ENDDO
               ENDDO
            ENDIF
!
            ELSE  !  Cloud is of convective type
!
            IF (L_CNV_CLOUD_EXTINCTION) THEN
               DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
                  DO LL=1, N_CLOUD_PROFILE(I)
                     L=I_CLOUD_PROFILE(LL, I)
                     CNV_CLOUD_EXTINCTION(L, I)                         &
     &                  =CNV_CLOUD_EXTINCTION(L, I)                     &
     &                  +K_EXT_TOT_CLOUD_COMP(L, I)                     &
     &                  *FRAC_CLOUD(L, I, I_CLOUD_TYPE(K))
                  ENDDO
               ENDDO
            ENDIF
!
!
            IF (L_CNV_CLOUD_ABSORPTIVITY) THEN
               DO I=N_CLOUD_TOP, N_LAYER
!CDIR NODEP
                  DO LL=1, N_CLOUD_PROFILE(I)
                     L=I_CLOUD_PROFILE(LL, I)
                     CNV_CLOUD_ABSORPTIVITY(L, I)                       &
     &                  =CNV_CLOUD_ABSORPTIVITY(L, I)                   &
     &                  +(K_EXT_TOT_CLOUD_COMP(L, I)                    &
     &                  -K_EXT_SCAT_CLOUD_COMP(L, I))                   &
     &                  *FRAC_CLOUD(L, I, I_CLOUD_TYPE(K))
                  ENDDO
               ENDDO
            ENDIF
            ENDIF
!
         ENDIF
!
      ENDDO
!
!
!
!
!     CALCULATE THE FINAL OPTICAL PROPERTIES.
!     THE SCATTERING WAS INCLUDED IN THE FREE TOTAL EXTINCTION EARLIER,
!     BUT WE HAVE YET TO DIVIDE THE PRODUCT OF THE ASYMMETRY AND THE
!     SCATTERING BY THE MEAN SCATTERING.
!
      DO I=1, N_LAYER
!
         N_INDEX=0
         DO L   =1,N_PROFILE
           IF (K_EXT_SCAT_FREE(L,I) >  TINY(K_EXT_SCAT_FREE)) THEN
             N_INDEX=N_INDEX+1
             INDEX(N_INDEX)=L
           END IF
         END DO
!
         IF(.NOT.L_RESCALE) THEN
!CDIR NODEP
            DO K=1, N_INDEX
               ASYMMETRY_FREE(INDEX(K), I)=ASYMMETRY_FREE(INDEX(K), I)  &
     &              /K_EXT_SCAT_FREE(INDEX(K), I)
            ENDDO
!
         ELSE IF (L_RESCALE.AND.N_INDEX >  0) THEN
            TMP_INV=1.0                                                 &
     &           /K_EXT_SCAT_FREE(INDEX(1), I)
            DO K=1, N_INDEX-1
               TMP_INV1=1.0                                             &
     &           /K_EXT_SCAT_FREE(INDEX(K+1), I)
               FORWARD_SCATTER_FREE(INDEX(K), I)                        &
     &            =FORWARD_SCATTER_FREE(INDEX(K), I)                    &
     &              *TMP_INV
               ASYMMETRY_FREE(INDEX(K), I)=ASYMMETRY_FREE(INDEX(K), I)  &
     &              *TMP_INV
               TMP_INV=TMP_INV1
            ENDDO
               FORWARD_SCATTER_FREE(INDEX(N_INDEX), I)                  &
     &            =FORWARD_SCATTER_FREE(INDEX(N_INDEX), I)              &
     &              *TMP_INV
               ASYMMETRY_FREE(INDEX(N_INDEX), I)=                       &
     &              ASYMMETRY_FREE(INDEX(N_INDEX), I)                   &
     &              *TMP_INV
         ENDIF
      ENDDO
!
!
!     REPEAT FOR CLOUDS.
      DO K=1, N_CLOUD_TYPE
         DO I=N_CLOUD_TOP, N_LAYER
!
            J      =1
            N_INDEX=0
            DO L   =1,N_PROFILE
              IF (K_EXT_SCAT_CLOUD(L,I,K) >                             &
     &          TINY(K_EXT_SCAT_CLOUD)) THEN
                INDEX(J)=L
                J       =J+1
                N_INDEX =N_INDEX+1
              END IF
            END DO

            IF(.NOT.L_RESCALE) THEN
!CDIR NODEP
              DO J=1, N_INDEX
                ASYMMETRY_CLOUD(INDEX(J), I, K)                         &
     &            =ASYMMETRY_CLOUD(INDEX(J), I, K)                      &
     &            /K_EXT_SCAT_CLOUD(INDEX(J), I, K)
              ENDDO

            ELSE IF (L_RESCALE.AND.N_INDEX >  0) THEN
               TMP_INV=1.0                                              &
     &              /K_EXT_SCAT_CLOUD(INDEX(1), I, K)
               DO J=1, N_INDEX-1
                  TMP_INV1=1.0                                          &
     &                 /K_EXT_SCAT_CLOUD(INDEX(J+1), I, K)
                  FORWARD_SCATTER_CLOUD(INDEX(J), I, K)                 &
     &               =FORWARD_SCATTER_CLOUD(INDEX(J), I, K)             &
     &                 *TMP_INV
                  ASYMMETRY_CLOUD(INDEX(J), I, K)                       &
     &                 =ASYMMETRY_CLOUD(INDEX(J), I, K)                 &
     &                 *TMP_INV
                  TMP_INV=TMP_INV1
               ENDDO
               FORWARD_SCATTER_CLOUD(INDEX(N_INDEX), I, K)              &
     &              =FORWARD_SCATTER_CLOUD(INDEX(N_INDEX), I, K)        &
     &              *TMP_INV
               ASYMMETRY_CLOUD(INDEX(N_INDEX), I, K)                    &
     &              =ASYMMETRY_CLOUD(INDEX(N_INDEX), I, K)              &
     &              *TMP_INV
            ENDIF

            IF (L_pc2) THEN
              DO J=1, N_INDEX
                ! Protect against asymmetry getting above 1
                IF (ASYMMETRY_CLOUD(INDEX(J), I, K) > 1.0) THEN
                  ASYMMETRY_CLOUD(INDEX(J), I, K) = 1.0
                ENDIF

                ! Protect against forward scatter getting above 1
                IF (FORWARD_SCATTER_CLOUD(INDEX(J),I,K) > 1.0) THEN
                  FORWARD_SCATTER_CLOUD(INDEX(J), I, K) = 1.0
                ENDIF
              ENDDO
            ENDIF  ! L_pc2

         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE GREY_EXTINCTION

