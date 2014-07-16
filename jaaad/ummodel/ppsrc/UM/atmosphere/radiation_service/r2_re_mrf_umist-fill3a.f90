

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to implement the MRF UMIST parametrization.
!
! Purpose:
!   Effective Radii are calculated in accordance with this
!   parametrization.
!
! Method:
!   The number density of CCN is found from the concentration
!   of aerosols, if available. This yields the number density of
!   droplets: if aerosols are not present, the number of droplets
!   is fixed. Effective radii are calculated from the number of
!   droplets and the LWC. Limits are applied to these values. In
!   deep convective clouds fixed values are assumed.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_RE_MRF_UMIST(N_PROFILE, N_LAYER, NCLDS              &
     &   , I_GATHER                                                     &
     &   , L_PC2, L_AEROSOL_CCN                                         &
     &   , L_BIOMASS_CCN, L_OCFF_CCN                                    &
     &   , SEA_SALT_FILM, SEA_SALT_JET                                  &
     &   , L_SEASALT_CCN, SALT_DIM_A, SALT_DIM_B                        &
     &   , L_USE_BIOGENIC, BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2       &
     &   , ACCUM_SULPHATE, DISS_SULPHATE, AITKEN_SULPHATE               &
     &   , AGED_BMASS, CLOUD_BMASS                                      &
     &   , AGED_OCFF, CLOUD_OCFF                                        &
     &   , LYING_SNOW_G                                                 &
     &   , I_CLOUD_REPRESENTATION                                       &
     &   , LAND_G, FLANDG_G                                             &
     &   , DENSITY_AIR, CONDENSED_MIX_RATIO, CC_DEPTH                   &
     &   , CONDENSED_RE                                                 &
     &   , D_MASS                                                       &
     &   , STRAT_LIQ_CLOUD_FRACTION                                     &
     &   , CONV_LIQ_CLOUD_FRACTION                                      &
     &   , NC_DIAG_FLAG                                                 &
     &   , NC_DIAG_G                                                    &
     &   , NC_WEIGHT_G                                                  &
     &   , NTOT_DIAG_G, Ntot_land, Ntot_sea                             &
     &   , STRAT_LWC_DIAG_G                                             &
     &   , SO4_CCN_DIAG_G                                               &
     &   , SULP_DIM1, SULP_DIM2                                         &
     &   , BMASS_DIM1, BMASS_DIM2                                       &
     &   , OCFF_DIM1, OCFF_DIM2                                         &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES       &
     &   )
!
!
!
      USE rad_switches_mod, ONLY:                                       &
          lrad_ccrad           
 
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED:
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------
! C_DENSTY for subroutine SF_EXCH
      REAL,PARAMETER:: RHOSEA = 1026.0 ! density of sea water (kg/m3)
      REAL,PARAMETER:: RHO_WATER = 1000.0! density of pure water (kg/m3)
! C_DENSTY end
!
! Description:
!
!  Contains various cloud droplet parameters, defined for
!  land and sea areas.
!
!  NTOT_* is the total number concentration (m-3) of cloud droplets;
!  KPARAM_* is the ratio of the cubes of the volume-mean radius
!                                           and the effective radius;
!  DCONRE_* is the effective radius (m) for deep convective clouds;
!  DEEP_CONVECTION_LIMIT_* is the threshold depth (m) bewteen shallow
!                                          and deep convective cloud.
!
! Current Code Owner: Andy Jones
!
! History:
!
! Version   Date     Comment
! -------   ----     -------
!    1     040894   Original code.    Andy Jones
!  5.2     111000   Updated in line with Bower et al. 1994 (J. Atmos.
!                   Sci., 51, 2722-2732) and subsequent pers. comms.
!                   Droplet concentrations now as used in HadAM4.
!                                     Andy Jones
!  5.4     02/09/02 Moved THOMO here from C_LSPMIC.      Damian Wilson
!  6.2     17/11/05 Remove variables that are now in UMUI. D. Wilson
!
!     REAL,PARAMETER:: NTOT_LAND is set in UMUI
!     REAL,PARAMETER:: NTOT_SEA is set in UMUI
      REAL,PARAMETER:: KPARAM_LAND = 0.67
      REAL,PARAMETER:: KPARAM_SEA = 0.80
      REAL,PARAMETER:: DCONRE_LAND = 9.5E-06
      REAL,PARAMETER:: DCONRE_SEA = 16.0E-06
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_LAND = 500.0
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_SEA = 1500.0
!
! Maximum Temp for homogenous nucleation (deg C)
      REAL,PARAMETER:: THOMO = -40.0
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
! CLDCMP3A sets components of clouds for two-stream radiation code.

      ! stratiform water droplets
      INTEGER,PARAMETER:: IP_CLCMP_ST_WATER=1

      ! stratiform ice crystals
      INTEGER,PARAMETER:: IP_CLCMP_ST_ICE=2

      ! convective water droplets
      INTEGER,PARAMETER:: IP_CLCMP_CNV_WATER=3

      ! convective ice crystals
      INTEGER,PARAMETER:: IP_CLCMP_CNV_ICE=4

! CLDCMP3A end
! CLREPP3A defines representations of clouds in two-stream radiation
! code.

      ! all components are mixed homogeneously
      INTEGER,PARAMETER:: IP_CLOUD_HOMOGEN     = 1

      ! ice and water clouds are treated separately
      INTEGER,PARAMETER:: IP_CLOUD_ICE_WATER   = 2

      ! clouds are divided into homogeneously mixed stratiform and
      ! convective parts
      INTEGER,PARAMETER:: IP_CLOUD_CONV_STRAT  = 3

      ! clouds divided into ice and water phases and into stratiform and
      ! convective components.
      INTEGER,PARAMETER:: IP_CLOUD_CSIW        = 4

! Types of clouds (values in CLREPD3A)

      ! number of type of clouds in representation
      INTEGER :: NP_CLOUD_TYPE(NPD_CLOUD_REPRESENTATION)

      ! map of components contributing to types of clouds
      INTEGER :: IP_CLOUD_TYPE_MAP(NPD_CLOUD_COMPONENT,                 &
     &  NPD_CLOUD_REPRESENTATION)

! CLREPP3A end
!
!
!     DUMMY ARGUMENTS:
!
!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             SIZE OF INPUT FIELDS TO THE RADIATION
     &   , NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_AEROSOL_SPECIES                                          &
!             MAXIMUM NUMBER OF AEROSOL SPECIES
     &   , SULP_DIM1                                                    &
!             1ST DIMENSION OF ARRAYS OF SULPHATE
     &   , SULP_DIM2                                                    &
!             2ND DIMENSION OF ARRAYS OF SULPHATE
     &   , BMASS_DIM1                                                   &
!             1ST DIMENSION OF ARRAYS OF BIOMASS SMOKE
     &   , BMASS_DIM2                                                   &
!             2ND DIMENSION OF ARRAYS OF BIOMASS SMOKE
     &   , OCFF_DIM1                                                    &
!             1ST DIMENSION OF ARRAYS OF FOSSIL-FUEL ORGANIC CARBON
     &   , OCFF_DIM2                                                    &
!             2ND DIMENSION OF ARRAYS OF FOSSIL-FUEL ORGANIC CARBON
     &   , SALT_DIM_A                                                   &
!             1ST DIMENSION OF ARRAYS OF SEA-SALT
     &   , SALT_DIM_B                                                   &
!             2ND DIMENSION OF ARRAYS OF SEA-SALT
     &   , BIOGENIC_DIM1                                                &
!             1ST DIMENSION OF BIOGENIC AEROSOL ARRAY
     &   , BIOGENIC_DIM2
!             2ND DIMENSION OF BIOGENIC AEROSOL ARRAY
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF ATMOSPHERIC PROFILES
     &   , N_LAYER                                                      &
!             Number of layers seen in radiation
     &   , NCLDS
!             NUMBER OF CLOUDY LEVELS
!
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO BE GATHERED
      LOGICAL                                                           &
                !, INTENT(IN)
     &     LAND_G(NPD_PROFILE)
!             GATHERED MASK FOR LAND POINTS
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_CLOUD_REPRESENTATION
!             REPRESENTATION OF CLOUDS
!
!     VARIABLES FOR PC2
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_PC2
!             PC2 cloud scheme is in use
!
!     VARIABLES FOR AEROSOLS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_AEROSOL_CCN                                                &
!             FLAG TO USE AEROSOLS TO FIND CCN.
     &   , L_SEASALT_CCN                                                &
!             FLAG TO USE SEA-SALT AEROSOL FOR CCN
     &   , L_USE_BIOGENIC                                               &
!             FLAG TO USE BIOGENIC AEROSOL FOR CCN
     &   , L_BIOMASS_CCN                                                &
!             FLAG TO USE BIOMASS SMOKE AEROSOL FOR CCN
     &   , L_OCFF_CCN                                                   &
!             FLAG TO USE FOSSIL-FUEL ORGANIC CARBON AEROSOL FOR CCN
     &   , NC_DIAG_FLAG
!             FLAG TO DIAGNOSE COLUMN-INTEGRATED DROPLET NUMBER
!
      REAL                                                              &
                !, INTENT(IN)
     &     ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)                         &
!             MIXING RATIOS OF ACCUMULATION MODE SULPHATE
     &   , AITKEN_SULPHATE(SULP_DIM1, SULP_DIM2)                        &
!             Mixing ratios of Aitken-mode sulphate
     &   , DISS_SULPHATE(SULP_DIM1, SULP_DIM2)                          &
!             MIXING RATIOS OF DISSOLVED SULPHATE
     &   , AGED_BMASS(BMASS_DIM1, BMASS_DIM2)                           &
!             MIXING RATIOS OF AGED BIOMASS SMOKE
     &   , CLOUD_BMASS(BMASS_DIM1, BMASS_DIM2)                          &
!             MIXING RATIOS OF IN-CLOUD BIOMASS SMOKE
     &   , AGED_OCFF(OCFF_DIM1, OCFF_DIM2)                              &
!             MIXING RATIOS OF AGED FOSSIL-FUEL ORGANIC CARBON
     &   , CLOUD_OCFF(OCFF_DIM1, OCFF_DIM2)                             &
!             MIXING RATIOS OF IN-CLOUD FOSSIL-FUEL ORGANIC CARBON
     &   , SEA_SALT_FILM(SALT_DIM_A, SALT_DIM_B)                        &
!             NUMBER CONCENTRATION OF FILM-MODE SEA-SALT AEROSOL
     &   , SEA_SALT_JET(SALT_DIM_A, SALT_DIM_B)                         &
!             NUMBER CONCENTRATION OF JET-MODE SEA-SALT AEROSOL
     &   , BIOGENIC(BIOGENIC_DIM1, BIOGENIC_DIM2)
!             M.M.R. OF BIOGENIC AEROSOL
!
      REAL                                                              &
                !, INTENT(IN)
     &     DENSITY_AIR(NPD_PROFILE, NPD_LAYER)
!             DENSITY OF AIR
!
      REAL                                                              &
                !, INTENT(IN)
     &     CONDENSED_MIX_RATIO(NPD_PROFILE, 0: NPD_LAYER                &
     &        , NPD_CLOUD_COMPONENT)                                    &
!             MIXING RATIOS OF CONDENSED SPECIES
     &   , CC_DEPTH(NPD_PROFILE)                                        &
!             DEPTH OF CONVECTIVE CLOUD
     &   , D_MASS(NPD_PROFILE, NPD_LAYER)                               &
!             MASS THICKNESS OF LAYER
     &   , STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)             &
!             STRATIFORM LIQUID CLOUD COVER IN LAYERS (T>273K)
     &   , CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)
!             CONVECTIVE LIQUID CLOUD COVER IN LAYERS (T>273K)



      Real                                                              &
                !, Intent(IN)
     &     Ntot_land                                                    &
                               ! Number of droplets over land / m-3
     &   , Ntot_sea            ! Number of droplets over sea / m-3

!
      REAL                                                              &
                !, INTENT(OUT)
     &     CONDENSED_RE(NPD_PROFILE, 0: NPD_LAYER, NPD_CLOUD_COMPONENT)
!             EFFECTIVE RADII OF CONDENSED COMPONENTS OF CLOUDS
!
      REAL                                                              &
                !, INTENT(OUT)
     &     NTOT_DIAG_G(NPD_PROFILE, NPD_LAYER)                          &
!             DIAGNOSTIC ARRAY FOR NTOT (GATHERED)
     &   , STRAT_LWC_DIAG_G(NPD_PROFILE, NPD_LAYER)                     &
!             DIAGNOSTIC ARRAY FOR STRATIFORM LWC (GATHERED)
     &   , SO4_CCN_DIAG_G(NPD_PROFILE, NPD_LAYER)                       &
!             DIAGNOSTIC ARRAY FOR SO4 CCN MASS CONC (GATHERED)
     &   , NC_DIAG_G(NPD_PROFILE)                                       &
!             DIAGNOSTIC ARRAY FOR INTEGRATED DROPLET NUMBER (GATHERED)
     &   , NC_WEIGHT_G(NPD_PROFILE)
!             DIAGNOSTIC ARRAY FOR INT DROP NO. SAMPLING WGT (GATHERED)
!
!
      REAL                                                              &
     &     LYING_SNOW_G(NPD_PROFILE)                                    &
!             GATHERED SNOW DEPTH (>5000m = LAND ICE SHEET)
     &   , FLANDG_G(NPD_PROFILE)
!             GATHERED global LAND FRACTION
!
!     LOCAL VARIABLES:
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , LG                                                           &
!             Index for gathering
     &   , IDUM(npd_field)                                              &
!             Temporary integer used to locate cloudy levels
     &   , SULPHATE_PTR_A                                               &
     &   , SULPHATE_PTR_B                                               &
!             POINTERS FOR SULPHATE ARRAYS
     &   , SEASALT_PTR_A                                                &
     &   , SEASALT_PTR_B                                                &
!             POINTERS FOR SEA-SALT ARRAYS
     &   , BIOMASS_PTR_A                                                &
     &   , BIOMASS_PTR_B                                                &
!             POINTERS FOR BIOMASS SMOKE ARRAYS
     &   , OCFF_PTR_A                                                   &
     &   , OCFF_PTR_B                                                   &
!             POINTERS FOR FOSSIL-FUEL ORGANIC CARBON ARRAYS
     &   , BIOGENIC_PTR_A                                               &
     &   , BIOGENIC_PTR_B
!             POINTERS FOR BIOGENIC ARRAY
!

      REAL                                                              &
     &     TOTAL_MIX_RATIO_ST(NPD_PROFILE)                              &
!             TOTAL MIXING RATIO OF WATER SUBSTANCE IN STRATIFORM CLOUD
     &   , TOTAL_MIX_RATIO_CNV(NPD_PROFILE)                             &
!             TOTAL MIXING RATIO OF WATER SUBSTANCE IN STRATIFORM CLOUD
     &   , TOTAL_STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE)                  &
!             TOTAL STRATIFORM LIQUID CLOUD COVER (T>273K)
     &   , TOTAL_CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE)
!             TOTAL CONVECTIVE LIQUID CLOUD COVER (T>273K)
!
      REAL                                                              &
     &     N_DROP(NPD_PROFILE, NPD_LAYER)                               &
!             NUMBER DENSITY OF DROPLETS
     &   , KPARAM                                                       &
!             RATIO OF CUBES OF VOLUME RADIUS TO EFFECTIVE RADIUS
     &   , TEMP
!             Temporary in calculation of effective radius
!
!     FIXED CONSTANTS OF THE PARAMETRIZATION:
      REAL                                                              &
     &     DEEP_CONVECTIVE_CLOUD
!             THRESHOLD VALUE FOR DEEP CONVECTIVE CLOUD
      PARAMETER(                                                        &
     &     DEEP_CONVECTIVE_CLOUD=5.0E+02                                &
     &   )
!
!
!     Subroutines called:
      EXTERNAL                                                          &
     &     R2_CALC_TOTAL_CLOUD_COVER                                    &
     &   , R2_COLUMN_DROPLET_CONC
!
!     Functions called:
      REAL                                                              &
     &     NUMBER_DROPLET
!             Function to calculate the number of clouds droplets
      EXTERNAL                                                          &
     &     NUMBER_DROPLET
!
!
!
!
!     CALCULATE THE NUMBER DENSITY OF DROPLETS
!
      IF (L_AEROSOL_CCN) THEN

         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE

               SULPHATE_PTR_A=I_GATHER(L)
               SULPHATE_PTR_B=N_LAYER+1-I

!  Diagnose SO4 aerosol concentrations. Mass mixing ratio of ammonium
!  sulphate is converted to microgrammes of the sulphate ion per m3
!  for diagnostic purposes.

               SO4_CCN_DIAG_G(L, I)=                                    &
                  (AITKEN_SULPHATE(SULPHATE_PTR_A, SULPHATE_PTR_B)      &
                   +ACCUM_SULPHATE(SULPHATE_PTR_A, SULPHATE_PTR_B)      &
                    +DISS_SULPHATE(SULPHATE_PTR_A, SULPHATE_PTR_B))     &
                    * DENSITY_AIR(L, I) * (96./132.) * 1.0E+09

            IF (L_SEASALT_CCN) THEN
               SEASALT_PTR_A=I_GATHER(L)
               SEASALT_PTR_B=N_LAYER+1-I
            ELSE
               SEASALT_PTR_A=1
               SEASALT_PTR_B=1
            ENDIF
            IF (L_BIOMASS_CCN) THEN
               BIOMASS_PTR_A=I_GATHER(L)
               BIOMASS_PTR_B=N_LAYER+1-I
            ELSE
               BIOMASS_PTR_A=1
               BIOMASS_PTR_B=1
            ENDIF
            IF (L_OCFF_CCN) THEN
               OCFF_PTR_A=I_GATHER(L)
               OCFF_PTR_B=N_LAYER+1-I
            ELSE
               OCFF_PTR_A=1
               OCFF_PTR_B=1
            ENDIF
            IF (L_USE_BIOGENIC) THEN
               BIOGENIC_PTR_A=I_GATHER(L)
               BIOGENIC_PTR_B=N_LAYER+1-I
            ELSE
               BIOGENIC_PTR_A=1
               BIOGENIC_PTR_B=1
            ENDIF
! DEPENDS ON: number_droplet
            N_DROP(L, I)=NUMBER_DROPLET(L_AEROSOL_CCN, .TRUE.           &
     &         , AITKEN_SULPHATE(SULPHATE_PTR_A, SULPHATE_PTR_B)        &
     &         , ACCUM_SULPHATE(SULPHATE_PTR_A, SULPHATE_PTR_B)         &
     &         , DISS_SULPHATE(SULPHATE_PTR_A, SULPHATE_PTR_B)          &
     &         , L_SEASALT_CCN                                          &
     &         , SEA_SALT_FILM(SEASALT_PTR_A, SEASALT_PTR_B)            &
     &         , SEA_SALT_JET(SEASALT_PTR_A, SEASALT_PTR_B)             &
     &         , L_USE_BIOGENIC                                         &
     &         , BIOGENIC(BIOGENIC_PTR_A, BIOGENIC_PTR_B)               &
     &         , L_BIOMASS_CCN                                          &
     &         , AGED_BMASS(BIOMASS_PTR_A, BIOMASS_PTR_B)               &
     &         , CLOUD_BMASS(BIOMASS_PTR_A, BIOMASS_PTR_B)              &
     &         , L_OCFF_CCN                                             &
     &         , AGED_OCFF(OCFF_PTR_A, OCFF_PTR_B)                      &
     &         , CLOUD_OCFF(OCFF_PTR_A, OCFF_PTR_B)                     &
     &         , DENSITY_AIR(L, I)                                      &
     &         , LYING_SNOW_G(L)                                        &
     &         , FLANDG_G(L)                                            &
     &         , Ntot_land, Ntot_sea                                    &
     &         )

            ENDDO
         ENDDO
      ELSE
!
!        Without aerosols, the number of droplets is fixed.
!
         DO L=1, N_PROFILE
            IF (FLANDG_G(L) >= 0.5) THEN
               N_DROP(L, N_LAYER+1-NCLDS:N_LAYER) = NTOT_LAND
            ELSE
               N_DROP(L, N_LAYER+1-NCLDS:N_LAYER) = NTOT_SEA
            ENDIF
         ENDDO
!
      ENDIF

!  Diagnose column-integrated cloud droplet number if required.

      IF (NC_DIAG_FLAG) THEN

! DEPENDS ON: r2_calc_total_cloud_cover
         CALL R2_CALC_TOTAL_CLOUD_COVER(N_PROFILE, N_LAYER, NCLDS       &
     &      , I_CLOUD_REPRESENTATION, STRAT_LIQ_CLOUD_FRACTION          &
     &      , TOTAL_STRAT_LIQ_CLOUD_FRACTION                            &
     &      , NPD_PROFILE, NPD_LAYER)

! DEPENDS ON: r2_calc_total_cloud_cover
         CALL R2_CALC_TOTAL_CLOUD_COVER(N_PROFILE, N_LAYER, NCLDS       &
     &      , I_CLOUD_REPRESENTATION, CONV_LIQ_CLOUD_FRACTION           &
     &      , TOTAL_CONV_LIQ_CLOUD_FRACTION                             &
     &      , NPD_PROFILE, NPD_LAYER)

! DEPENDS ON: r2_column_droplet_conc
         CALL R2_COLUMN_DROPLET_CONC(NPD_PROFILE, NPD_LAYER             &
     &      , N_PROFILE, N_LAYER, NCLDS                                 &
     &      , STRAT_LIQ_CLOUD_FRACTION                                  &
     &      , TOTAL_STRAT_LIQ_CLOUD_FRACTION                            &
     &      , CONV_LIQ_CLOUD_FRACTION                                   &
     &      , TOTAL_CONV_LIQ_CLOUD_FRACTION                             &
     &      , N_DROP, D_MASS, DENSITY_AIR                               &
     &      , NC_DIAG_G, NC_WEIGHT_G)

      ENDIF

!
      DO I=N_LAYER+1-NCLDS, N_LAYER
!
!        FIND THE TOTAL MIXING RATIO OF WATER SUBSTANCE IN THE CLOUD
!        AS IMPLIED BY THE REPRESENTATION.
         IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT) THEN
            DO L=1, N_PROFILE
               TOTAL_MIX_RATIO_ST(L)                                    &
     &            =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_ST_WATER)         &
     &            +CONDENSED_MIX_RATIO(L, I, IP_CLCMP_ST_ICE)
               TOTAL_MIX_RATIO_CNV(L)                                   &
     &            =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_CNV_WATER)        &
     &            +CONDENSED_MIX_RATIO(L, I, IP_CLCMP_CNV_ICE)
            ENDDO
         ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
            DO L=1, N_PROFILE
               TOTAL_MIX_RATIO_ST(L)                                    &
     &            =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_ST_WATER)
               TOTAL_MIX_RATIO_CNV(L)                                   &
     &            =CONDENSED_MIX_RATIO(L, I, IP_CLCMP_CNV_WATER)
            ENDDO
         ENDIF
         IF (L_PC2) THEN
           DO L=1, N_PROFILE
              IF (LAND_G(L)) THEN
                 KPARAM=KPARAM_LAND
              ELSE
                 KPARAM=KPARAM_SEA
              ENDIF
              TEMP                                                      &
     &         =(3.0E+00*TOTAL_MIX_RATIO_CNV(L)*DENSITY_AIR(L, I)       &
     &         /(4.0E+00*PI*RHO_WATER*KPARAM*N_DROP(L, I)))
              IF (TEMP  >=  0.0) THEN
                CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)                  &
     &          = TEMP**(1.0E+00/3.0E+00)
              ELSE
                CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER) = 0.0
              END IF
              TEMP                                                      &
     &         =(3.0E+00*TOTAL_MIX_RATIO_ST(L)*DENSITY_AIR(L, I)        &
     &         /(4.0E+00*PI*RHO_WATER*KPARAM*N_DROP(L, I)))
              IF (TEMP  >=  0.0) THEN
                CONDENSED_RE(L, I, IP_CLCMP_ST_WATER)                   &
     &          = TEMP**(1.0E+00/3.0E+00)
              ELSE
                CONDENSED_RE(L, I, IP_CLCMP_ST_WATER) = 0.0
              END IF
           ENDDO
         ELSE
           DO L=1, N_PROFILE
              IF (LAND_G(L)) THEN
                 KPARAM=KPARAM_LAND
              ELSE
                 KPARAM=KPARAM_SEA
              ENDIF
              CONDENSED_RE(L, I, IP_CLCMP_CNV_WATER)                    &
     &         =(3.0E+00*TOTAL_MIX_RATIO_CNV(L)*DENSITY_AIR(L, I)       &
     &         /(4.0E+00*PI*RHO_WATER*KPARAM*N_DROP(L, I)))             &
     &         **(1.0E+00/3.0E+00)
              CONDENSED_RE(L, I, IP_CLCMP_ST_WATER)                     &
     &         =(3.0E+00*TOTAL_MIX_RATIO_ST(L)*DENSITY_AIR(L, I)        &
     &         /(4.0E+00*PI*RHO_WATER*KPARAM*N_DROP(L, I)))             &
     &         **(1.0E+00/3.0E+00)
           ENDDO
         ENDIF
         DO L=1, N_PROFILE
            NTOT_DIAG_G(L, I)=N_DROP(L, I)*1.0E-06
            STRAT_LWC_DIAG_G(L, I)                                      &
     &         =TOTAL_MIX_RATIO_ST(L)*DENSITY_AIR(L, I)*1.0E03
         ENDDO
      ENDDO

      !-----------------------------------------------------------------------
      ! Reset the effective radii for deep convective clouds.
      !-----------------------------------------------------------------------

      IF (lrad_ccrad) THEN

        ! Loop over all gridpoints which are affected by radiation
        Do L=1, N_PROFILE 

          CC_DEPTH(L) = 0.0E+00    
          LG          = I_GATHER(L)
          idum        = 0            

          ! Loop from top layer down i.e. radiation indexes from the top of
          ! atmosphere down
          DO I=N_LAYER+1-NCLDS, N_LAYER


            IF (CONDENSED_MIX_RATIO(L,I,IP_CLCMP_CNV_WATER) > 0.0E+00) THEN

              IF (idum(LG) == 0) THEN
                ! Mark 1st cloud top (i.e highest cloud top in profile)
                idum(LG) = I
              END IF 

            ELSE

              IF (idum(LG) /= 0) THEN
                ! Cloud top already found so this is cloud base and have found
                ! a contigious cloud bank, begin to calculate depth

                DO J=idum(LG), I-1
                  CC_DEPTH(L) = CC_DEPTH(L)                                   &
                              + ( D_MASS(L,J)/DENSITY_AIR(L,J) )
                END DO


                ! Reset effective radii if clouds are physically deep
                DO J=idum(LG), I-1
                  IF (LAND_G(L)) THEN
               
                    ! LAND POINT
                    IF ((CC_DEPTH(L) > DEEP_CONVECTION_LIMIT_LAND) .OR.      &
                        (CONDENSED_RE(L,J,IP_CLCMP_CNV_WATER) > DCONRE_LAND))& 
                        THEN 

                      CONDENSED_RE(L,J,IP_CLCMP_CNV_WATER) = DCONRE_LAND
                    END IF

                  ELSE
    
                    ! SEA POINT
                    IF ((CC_DEPTH(L) > DEEP_CONVECTION_LIMIT_SEA) .OR.       &
                        (CONDENSED_RE(L,J,IP_CLCMP_CNV_WATER) > DCONRE_SEA)) & 
                        THEN 

                      CONDENSED_RE(L,J,IP_CLCMP_CNV_WATER) = DCONRE_SEA

                    END IF

                  END IF      ! (LAND_G)
                END DO      ! J

                ! Reset cloud top and carry on looking for cloud
                idum(LG) = 0 
           
              END IF      ! (idum)
            END IF      ! (CONDENSED_MIX_RATIO)
          END DO      ! I (N_LAYER)
        END DO      ! L (N_PROFILE)

      ELSE ! original code

        DO I=N_LAYER+1-NCLDS, N_LAYER
          DO L=1, N_PROFILE

            IF (LAND_G(L)) THEN

              IF (CC_DEPTH(L) > DEEP_CONVECTION_LIMIT_LAND) THEN
                CONDENSED_RE(L,I,IP_CLCMP_CNV_WATER) = DCONRE_LAND
              ELSE
                IF (CONDENSED_RE(L,I,IP_CLCMP_CNV_WATER) > DCONRE_LAND) THEN
                  CONDENSED_RE(L,I,IP_CLCMP_CNV_WATER) = DCONRE_LAND
                END IF
              END IF

            ELSE

              IF (CC_DEPTH(L) > DEEP_CONVECTION_LIMIT_SEA) THEN
                CONDENSED_RE(L,I,IP_CLCMP_CNV_WATER) = DCONRE_SEA
              ELSE
                IF (CONDENSED_RE(L,I,IP_CLCMP_CNV_WATER) > DCONRE_SEA) THEN 
                  CONDENSED_RE(L,I,IP_CLCMP_CNV_WATER) = DCONRE_SEA
                END IF
              END IF

            END IF      ! (LAND_G)

          END DO      ! L (N_PROFILE)
        END DO      ! I (N_LAYER)
      END IF      ! lrad_ccrad
!
!
!
      RETURN
      END SUBROUTINE R2_RE_MRF_UMIST
