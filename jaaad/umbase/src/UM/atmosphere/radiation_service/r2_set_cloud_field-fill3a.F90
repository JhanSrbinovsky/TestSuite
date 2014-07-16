#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the mixing ratios of gases.
!
! Purpose:
!   The full array of mass mixing ratios of gases is filled.
!
! Method:
!   The arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. For well-mixed
!   gases the constant mixing ratios are fed into this array.
!
! Current Owner of Code: J. M. Edwards
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set thermodynamic properties
!
! Purpose:
!   Pressures, temperatures at the centres and edges of layers
!   and the masses in layers are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to assign Properties of Clouds.
!
! Purpose:
!   The fractions of different types of clouds and their microphysical
!   preoperties are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_CLOUD_FIELD(N_PROFILE, NLEVS, N_LAYER, NCLDS    &
     &   , I_GATHER                                                     &
     &   , P, T, D_MASS                                                 &
     &   , CCB_IN, CCT_IN, CCA, CCCWP, CCW, LCBASE                      &
     &   , LCCWC1, LCCWC2, LCA_AREA, LCA_BULK                           &
     &   , L_PC2, L_MICROPHYSICS, L_AEROSOL_CCN                         &
     &   , SEA_SALT_FILM, SEA_SALT_JET                                  &
     &   , L_SEASALT_CCN, SALT_DIM_A, SALT_DIM_B                        &
     &   , L_USE_BIOGENIC, BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2       &
     &   , SULP_DIM1, SULP_DIM2, ACCUM_SULPHATE, DISS_SULPHATE          &
     &   , AITKEN_SULPHATE, L_BIOMASS_CCN                               &
     &   , BMASS_DIM1, BMASS_DIM2, AGED_BMASS, CLOUD_BMASS              &
     &   , L_OCFF_CCN, OCFF_DIM1, OCFF_DIM2, AGED_OCFF, CLOUD_OCFF      &
     &   , LYING_SNOW                                                   &
     &   , L_CLOUD_WATER_PARTITION, LAND_G, FLANDG_G                    &
     &   , I_CLOUD, I_CLOUD_REPRESENTATION, I_CONDENSED_PARAM           &
     &   , CONDENSED_MIN_DIM, CONDENSED_MAX_DIM                         &
     &   , N_CONDENSED, TYPE_CONDENSED                                  &
     &   , W_CLOUD, FRAC_CLOUD, L_LOCAL_CNV_PARTITION                   &
     &   , CONDENSED_MIX_RAT_AREA, CONDENSED_DIM_CHAR                   &
     &   , RE_CONV, RE_CONV_FLAG, RE_STRAT, RE_STRAT_FLAG               &
     &   , WGT_CONV, WGT_CONV_FLAG, WGT_STRAT, WGT_STRAT_FLAG           &
     &   , LWP_STRAT, LWP_STRAT_FLAG                                    &
     &   , NTOT_DIAG, NTOT_DIAG_FLAG                                    &
     &   , STRAT_LWC_DIAG, STRAT_LWC_DIAG_FLAG                          &
     &   , SO4_CCN_DIAG, SO4_CCN_DIAG_FLAG                              &
     &   , COND_SAMP_WGT, COND_SAMP_WGT_FLAG                            &
     &   , NC_DIAG, NC_DIAG_FLAG                                        &
     &   , NC_WEIGHT, NC_WEIGHT_FLAG                                    &
     &   , col_list, row_list, row_length, rows                         &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES       &
     &   , N_CCA_LEV, Ntot_land, Ntot_sea                               &
     &   )


      Use cv_cntl_mod, Only:                                            &
          lcv_3d_cca
      Use rad_switches_mod, ONLY:                                       &
          lrad_3d_ccw, lrad_ovrlap, lrad_ccw_scav, lrad_cldtop_t_fix,   &
          lrad_ccrad, lrad_tripleclouds

      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
#include "dimfix3a.h"
#include "cldcmp3a.h"
#include "cldtyp3a.h"
#include "clrepp3a.h"
#include "iclprm3a.h"
#include "c_0_dg_c.h"
#include "c_r_cp.h"
      Logical L_climat_aerosol, L_CLIM_AERO_HGT, L_HADGEM1_CLIM_AERO &
     &   , L_OVRLAP, L_CCW_SCAV, L_CLDTOP_T_FIX, L_SEC_VAR, L_EQT
      Real FW_STD, A_SW_SEGMENTS, A_LW_SEGMENTS, A_SW_SEG_SIZE, A_LW_SEG_SIZE &
     &   , CO2_MMR,INHOM_CLOUD_SW(NPD_CLOUD_COMPONENT)                  &
     &   , INHOM_CLOUD_LW(NPD_CLOUD_COMPONENT), DP_CORR_STRAT,DP_CORR_CONV

!     real FW_STD
      integer AERO_BL_LEVELS 
      COMMON  /RUN_Radiation/L_climat_aerosol, L_clim_aero_hgt,         &
     &  L_HadGEM1_Clim_Aero, FW_STD,                                    &
     &  A_SW_SEGMENTS,A_SW_SEG_SIZE,A_LW_SEGMENTS,A_LW_SEG_SIZE,CO2_MMR,&
     &  L_SEC_VAR,L_EqT,INHOM_CLOUD_SW,INHOM_CLOUD_LW,DP_CORR_STRAT,    &
     &  DP_CORR_CONV,AERO_BL_LEVELS, l_ovrlap, l_ccw_scav, l_cldtop_t_fix

!
#include "cloud_scheme_pcf3z.h"
!
!
!     DIMENSIONS OF ARRAYS:
      Integer, Intent(IN) :: row_length
!                              Number of grid-points in EW-direction
!                              in the local domain
      Integer, Intent(IN) :: rows
!                              Number of grid-points in NS-direction
!                              in the local domain
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             FIELD SIZE IN CALLING PROGRAM
     &   , NPD_PROFILE                                                  &
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER                                                    &
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_AEROSOL_SPECIES                                          &
!             MAXIMUM NUMBER OF AEROSOL_SPECIES
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
     &   , BIOGENIC_DIM2                                                &
!             2ND DIMENSION OF BIOGENIC AEROSOL ARRAY
     &   , N_CCA_LEV
!             NUMBER OF LEVELS FOR CONVECTIVE CLOUD AMOUNT
!
!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , NLEVS                                                        &
!             Number of layers used outside the radiation scheme
     &   , N_LAYER                                                      &
!             Number of layers seen by the radiation code
     &   , NCLDS
!             NUMBER OF CLOUDY LEVELS
!
!     GATHERING ARRAY:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO BE GATHERED
      Integer, Intent(IN) :: col_list(npd_field)
!                              EW indices of gathered points in the 2-D
!                              domain
      Integer, Intent(IN) :: row_list(npd_field)
!                              NS indices of gathered points in the 2-D
!                              domain
!
!     THERMODYNAMIC FIELDS:
      REAL                                                              &
                !, INTENT(IN)
     &     P(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             PRESSURES
     &   , T(NPD_PROFILE, 0: NPD_LAYER)                                 &
!             TEMPERATURES
     &   , D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESSES OF LAYERS
!
!     CONVECTIVE CLOUDS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     CCB_IN(NPD_FIELD)                                            &
!             BASE OF CONVECTIVE CLOUD
     &   , CCT_IN(NPD_FIELD)
!             TOP OF CONVECTIVE CLOUD

      INTEGER, INTENT(IN):: lcbase(npd_field) 
!             Lowest cloud base level in vertical profile.


      REAL                                                              &
                !, INTENT(IN)
     &     CCA(NPD_FIELD,N_CCA_LEV)                                     &
!             FRACTION OF CONVECTIVE CLOUD
     &   , CCCWP(NPD_FIELD)
!             WATER PATH OF CONVECTIVE CLOUD

      Real, INTENT(IN) :: ccw(npd_field, nlevs) 
!             Convective Cloud Water (kg/kg) (If CCRad =.true.)

      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_LOCAL_CNV_PARTITION                                        &
!             FLAG TO CARRY OUT THE PARTITIONING BETWEEN ICE
!             AND WATER IN CONVECTIVE CLOUDS AS A FUNCTION OF
!             THE LOCAL TEMPERATURE
     &   , L_SEASALT_CCN                                                &
!              FLAG FOR SEA-SALT PARAMETRIZATION FOR CCN
     &   , L_USE_BIOGENIC                                               &
!              FLAG TO USE BIOGENIC AEROSOLS AS CCN
     &   , L_BIOMASS_CCN                                                &
!              FLAG FOR BIOMASS PARAMETRIZATION FOR CCN
     &   , L_OCFF_CCN
!              FLAG FOR FOSSIL-FUEL ORG CARB PARAMETRIZATION FOR CCN
!
!     LAYER CLOUDS:
      REAL                                                              &
                !, INTENT(IN)
     &     LCCWC1(NPD_FIELD, NCLDS+1/(NCLDS+1))                         &
!             LIQUID WATER CONTENTS
     &   , LCCWC2(NPD_FIELD, NCLDS+1/(NCLDS+1))                         &
!             ICE WATER CONTENTS
     &   , LCA_AREA(NPD_FIELD, NCLDS+1/(NCLDS+1))                       &
!             AREA COVERAGE FRACTIONS OF LAYER CLOUDS
     &   , LCA_BULK(NPD_FIELD, NCLDS+1/(NCLDS+1))
!             BULK COVERAGE FRACTIONS OF LAYER CLOUDS
!
!     ARRAYS FOR MICROPHYSICS:
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_MICROPHYSICS                                               &
!             MICROPHYSICAL FLAG
     &   , L_PC2                                                        &
!             PC2 cloud scheme is in use
     &   , L_AEROSOL_CCN                                                &
!             FLAG TO USE AEROSOLS TO FIND CCN
     &   , L_CLOUD_WATER_PARTITION                                      &
!             FLAG TO USE PROGNOSTIC CLOUD ICE CONTENTS
     &   , LAND_G(NPD_PROFILE)
!             FLAG FOR LAND POINTS
      REAL                                                              &
                !, INTENT(IN)
     &     ACCUM_SULPHATE(SULP_DIM1, SULP_DIM2)                         &
!             MIXING RATIOS OF ACCUMULATION-MODE SULPHATE
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
!     REPRESENTATION OF CLOUDS
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_CLOUD                                                      &
!             CLOUD SCHEME USED
     &   , I_CLOUD_REPRESENTATION
!             REPRESENTATION OF CLOUDS
!
!     PARAMETRIZATIONS FOR CLOUDS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             TYPES OF PARAMETRIZATION USED FOR CONDENSED
!             COMPONENTS IN CLOUDS
!     LIMITS ON SIZES OF PARTICLES
      REAL                                                              &
                !, INTENT(IN)
     &     CONDENSED_MIN_DIM(NPD_CLOUD_COMPONENT)                       &
!             MINIMUM DIMENSION OF EACH CONDENSED COMPONENT
     &   , CONDENSED_MAX_DIM(NPD_CLOUD_COMPONENT)
!             MAXIMUM DIMENSION OF EACH CONDENSED COMPONENT
!
      Real                                                              &
                !, Intent(IN)
     &     Ntot_land                                                    &
                               ! Number of droplets over land / m-3
     &   , Ntot_sea            ! Number of droplets over sea / m-3
!
!     ASSIGNED CLOUD FIELDS:
      INTEGER                                                           &
                !, INTENT(OUT)
     &     N_CONDENSED                                                  &
!             NUMBER OF CONDENSED COMPONENTS
     &   , TYPE_CONDENSED(NPD_CLOUD_COMPONENT)
!             TYPES OF CONDENSED COMPONENTS
      REAL                                                              &
                !, INTENT(OUT)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)                              &
!             TOTAL AMOUNTS OF CLOUD
     &   , FRAC_CLOUD(NPD_PROFILE, NPD_LAYER, NPD_CLOUD_TYPE)           &
!             FRACTION OF EACH TYPE OF CLOUD
     &   , CONDENSED_DIM_CHAR(NPD_PROFILE, 0: NPD_LAYER                 &
     &      , NPD_CLOUD_COMPONENT)                                      &
!             CHARACTERISTIC DIMENSIONS OF CLOUDY COMPONENTS
     &   , CONDENSED_MIX_RAT_AREA(NPD_PROFILE, 0: NPD_LAYER             &
     &      , NPD_CLOUD_COMPONENT)                                      &
!             MASS MIXING RATIOS OF CONDENSED COMPONENTS USING AREA CLD
     &   , NTOT_DIAG_G(NPD_PROFILE, NPD_LAYER)                          &
!             DIAGNOSTIC ARRAY FOR NTOT (GATHERED)
     &   , STRAT_LWC_DIAG_G(NPD_PROFILE, NPD_LAYER)                     &
!             DIAGNOSTIC ARRAY FOR STRATIFORM LWC (GATHERED)
     &   , SO4_CCN_DIAG_G(NPD_PROFILE, NPD_LAYER)
!             DIAGNOSTIC ARRAY FOR SO4 CCN MASS CONC (GATHERED)
!
!
      REAL                                                              &
     &     LYING_SNOW(NPD_FIELD)                                        &
!            SNOW DEPTH (>5000m = LAND ICE SHEET)
     &   , LYING_SNOW_G(NPD_PROFILE)                                    &
!            GATHERED VERSION OF THE ABOVE
     &   , FLANDG_G(NPD_PROFILE)
!            GATHERED GLOBAL LAND FRACTION FIELD
!
!     MICROPHYSICAL DIAGNOSTICS:
      LOGICAL                                                           &
     &     RE_CONV_FLAG                                                 &
!             DIAGNOSE EFFECTIVE RADIUS*WEIGHT FOR CONVECTIVE CLOUD
     &   , RE_STRAT_FLAG                                                &
!             DIAGNOSE EFFECTIVE RADIUS*WEIGHT FOR STRATIFORM CLOUD
     &   , WGT_CONV_FLAG                                                &
!             DIAGNOSE WEIGHT FOR CONVECTIVE CLOUD
     &   , WGT_STRAT_FLAG                                               &
!             DIAGNOSE WEIGHT FOR STRATIFORM CLOUD
     &   , LWP_STRAT_FLAG                                               &
!             DIAGNOSE LIQUID WATER PATH*WEIGHT FOR STRATIFORM CLOUD
     &   , NTOT_DIAG_FLAG                                               &
!             DIAGNOSE DROPLET CONCENTRATION*WEIGHT
     &   , STRAT_LWC_DIAG_FLAG                                          &
!             DIAGNOSE STRATIFORM LWC*WEIGHT
     &   , SO4_CCN_DIAG_FLAG                                            &
!             DIAGNOSE SO4 CCN MASS CONC*COND. SAMP. WEIGHT
     &   , COND_SAMP_WGT_FLAG                                           &
!             DIAGNOSE CONDITIONAL SAMPLING WEIGHT
     &   , NC_DIAG_FLAG                                                 &
!             DIAGNOSE COLUMN DROPLET CONCENTRATION * SAMP. WEIGHT
     &   , NC_WEIGHT_FLAG
!             DIAGNOSE COLUMN DROPLET SAMPLING WEIGHT
!
      REAL                                                              &
     &     RE_CONV(row_length, rows, NCLDS)                             &
!             EFFECTIVE RADIUS*WEIGHT FOR CONVECTIVE CLOUD
     &   , RE_STRAT(row_length, rows, NCLDS)                            &
!             EFFECTIVE RADIUS*WEIGHT FOR STRATIFORM CLOUD
     &   , WGT_CONV(row_length, rows, NCLDS)                            &
!             WEIGHT FOR CONVECTIVE CLOUD
     &   , WGT_STRAT(row_length, rows, NCLDS)                           &
!             WEIGHT FOR STRATIFORM CLOUD
     &   , LWP_STRAT(row_length, rows, NCLDS)                           &
!             LIQUID WATER PATH*WEIGHT FOR STRATIFORM CLOUD
     &   , NTOT_DIAG(row_length, rows, NCLDS)                           &
!             DROPLET CONCENTRATION*WEIGHT
     &   , STRAT_LWC_DIAG(row_length, rows, NCLDS)                      &
!             STRATIFORM LWC*WEIGHT
     &   , SO4_CCN_DIAG(row_length, rows, NCLDS)                        &
!             SO4 CCN MASS CONC*COND. SAMP. WEIGHT
     &   , COND_SAMP_WGT(row_length, rows, NCLDS)                       &
!             CONDITIONAL SAMPLING WEIGHT
     &   , NC_DIAG(row_length, rows)                                    &
!             COLUMN DROPLET CONCENTRATION * SAMPLING WEIGHT
     &   , NC_WEIGHT(row_length, rows)
!             COLUMN DROPLET CONCENTRATION SAMPLING WEIGHT
!
!
!
!     LOCAL VARIABLES:
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , J                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , LG                                                           &
!             INDEX TO GATHER
     &   , MULT_CC_FLAG
!             FLAG FOR EXISTENCE OF MULTIPLE VERTICAL CONVECTIVE CLOUDS

      LOGICAL                                                           &
     &     L_GLACIATED_TOP(NPD_PROFILE)
!             LOGICAL FOR GLACIATED TOPS IN CONVECTIVE CLOUD.
!
      REAL                                                              &
     &     LIQ_FRAC(NPD_PROFILE)                                        &
!             FRACTION OF LIQUID CLOUD WATER
     &   , LIQ_FRAC_CONV(NPD_PROFILE)                                   &
!             FRACTION OF LIQUID WATER IN CONVECTIVE CLOUD
     &   , T_GATHER(NPD_PROFILE)                                        &
!             GATHERED TEMPERATURE FOR LSP_FOCWWIL
     &   , T_LIST(NPD_PROFILE)                                          &
!             LIST OF TEMPERATURES
     &   , TOTAL_MASS(NPD_PROFILE)                                      &
!             TOTAL MASS IN CONVECTIVE CLOUD
     &   , CC_DEPTH(NPD_PROFILE)                                        &
!             DEPTH OF CONVECTIVE CLOUD
     &   , CONDENSED_MIX_RAT_BULK(NPD_PROFILE, 0: NPD_LAYER             &
     &      , NPD_CLOUD_COMPONENT)                                      &
!             MASS MIXING RATIOS OF CONDENSED COMPONENTS USING BULK CLD
     &   , DENSITY_AIR(NPD_PROFILE, NPD_LAYER)                          &
!             DENSITY OF AIR
     &   , CONVECTIVE_CLOUD_LAYER(NPD_PROFILE)                          &
!             AMOUNT OF CONVECTIVE CLOUD IN TH CURRENT LAYER
     &   , STRAT_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)             &
!             STRATIFORM LIQUID CLOUD FRACTION (T>273K)
     &   , CONV_LIQ_CLOUD_FRACTION(NPD_PROFILE, NPD_LAYER)              &
!             CONVECTIVE LIQUID CLOUD FRACTION (T>273K)
     &   , NC_DIAG_G(NPD_PROFILE)                                       &
!             DIAGNOSTIC ARRAY FOR COLUMN DROPLET NUMBER (GATHERED)
     &   , NC_WEIGHT_G(NPD_PROFILE)
!             DIAGNOSTIC ARRAY FOR COL. DROP NO. SAMPLING WGT (GATHERED)
!
      REAL                                                              &
     &     FRAC_CLOUD_TMP(NPD_PROFILE, NPD_LAYER                        &
     &      , NPD_CLOUD_COMPONENT)                                      &
!             TEMPORARY STORAGE OF CLOUD FRACTION FIELD
     &   , CONDENSED_MIX_RAT_AREA_TMP(NPD_PROFILE, 0: NPD_LAYER         &
     &      , NPD_CLOUD_COMPONENT)                                      &
!             TEMPORARY STORAGE OF MASS MIXING RATIO
     &   , CONDENSED_MIX_RAT_BULK_TMP(NPD_PROFILE, 0: NPD_LAYER         &
     &      , NPD_CLOUD_COMPONENT)
!             TEMPORARY STORAGE OF MASS MIXING RATIO
! 
!

      !-----------------------------------------------------------------
      ! Convection variables - CCRad
      !-----------------------------------------------------------------
      INTEGER :: ccb(npd_field)
!             CONVECTIVE CLOUD BASE TO BE USED BY RADIATION.

      INTEGER :: cct(npd_field)
!             CONVECTIVE CLOUD TOP TO BE USED BY RADIATION.

      REAL    :: lca_of_grdbx (npd_field, MAX(nclds,1))
!             LARGE-SCALE CLOUD AREA OF GRIDBOX.
      REAL    :: IWC
!                Ice water content

!
!     Parameters for the aggregate parametrization.
      REAL, Parameter :: a0_agg_cold = 7.5094588E-04
      REAL, Parameter :: b0_agg_cold = 5.0830326E-07
      REAL, Parameter :: a0_agg_warm = 1.3505403E-04
      REAL, Parameter :: b0_agg_warm = 2.6517429E-05
      REAL, Parameter :: t_switch    = 216.208
      REAL, Parameter :: t0_agg      = 279.5
      REAL, Parameter :: s0_agg      = 0.05
!



!     SET THE COMPONENTS WITHIN THE CLOUDS. IN THE UNIFIED MODEL WE
!     HAVE FOUR COMPONENTS: STRATIFORM ICE AND WATER AND CONVECTIVE
!     ICE AND WATER.
      N_CONDENSED=4
      TYPE_CONDENSED(1)=IP_CLCMP_ST_WATER
      TYPE_CONDENSED(2)=IP_CLCMP_ST_ICE
      TYPE_CONDENSED(3)=IP_CLCMP_CNV_WATER
      TYPE_CONDENSED(4)=IP_CLCMP_CNV_ICE

      IF (lrad_cldtop_t_fix) THEN
        DO L=1, N_PROFILE 
          L_GLACIATED_TOP(L) = .FALSE. 
          T_GATHER(L)        =  0.0E+00 
        END DO 
      END IF      ! lrad_cldtop_t_fix

!     SET THE TOTAL AMOUNTS OF CLOUD AND THE FRACTIONS COMPRISED BY
!     CONVECTIVE AND STRATIFORM COMPONENTS.
!
!     ZERO THE AMOUNTS OF CLOUD IN THE UPPER LAYERS.
      DO I=1, N_LAYER-NCLDS
         DO L=1, N_PROFILE
            W_CLOUD(L, I)=0.0E+00
         ENDDO
      ENDDO
      

      IF (lrad_ccrad) THEN
        !---------------------------------------------------------------------
        ! Convection passes in ccb_in and cct_in on model layers, and these
        ! only to refer the top most contigiuos cloud bank.  Set local
        ! convective cloud base and top to:
        !
        !   ccb = LOWER BOUNDARY of LOWEST  LAYER with convective cloud
        !   cct = UPPER BOUNDARY of HIGHEST LAYER with convective cloud
        !
        ! as expected by the radiation scheme.
        !---------------------------------------------------------------------
        MULT_CC_FLAG = 0

        DO L=1, N_PROFILE
          LG=I_GATHER(L)
          
          CCB(LG) = LCBASE(LG)     ! Set ccb to lowest cloud boundary
                                   ! of lowest cloud 
          CCT(LG) = CCT_IN(LG) + 1 ! Set cct to upper boundary of highest
                                   ! Cloud

          mult_cc_flag = mult_cc_flag + ccb_in(lg) - lcbase(lg)
        END DO

        IF ((MULT_CC_FLAG /= 0) .AND. (.NOT. lrad_3d_ccw)) THEN
          WRITE(6,*) ''
          WRITE(6,*) ' Warning: Multiple convective clouds on gridboxs exists.'
          WRITE(6,*) '          L_3d_ccw is set to FALSE. Vertical cloud      '
          WRITE(6,*) '          structure seen by radiation will not be       '
          WRITE(6,*) '          consistent with that from convection scheme.  '
          WRITE(6,*) ''
        ENDIF

      ELSE ! Original 

        DO L=1, N_PROFILE
          LG=I_GATHER(L)
          CCB(LG) = CCB_IN(LG)
          CCT(LG) = CCT_IN(LG)
        END DO

      END IF      ! lrad_ccrad
      

!
      IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT .AND.           &
     &    .NOT. L_CLOUD_WATER_PARTITION) THEN
!  This cloud representation not available with new cloud microphysics
!
!        THE CLOUDS ARE DIVIDED INTO MIXED-PHASE STRATIFORM AND
!        CONVECTIVE CLOUDS: LSP_FOCWWIL GIVES THE PARTITIONING BETWEEN
!        ICE AND WATER IN STRATIFORM CLOUDS AND IN CONVECTIVE CLOUD,
!        UNLESS THE OPTION TO PARTITION AS A FUNCTION OF THE LOCAL
!        TEMPERATURE IS SELECTED. WITHIN CONVECTIVE CLOUD THE LIQUID
!        WATER CONTENT IS DISTRIBUTED UNIFORMLY THROUGHOUT THE CLOUD.
!
!        CONVECTIVE CLOUD:
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)=0.0E+00
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)=0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)=0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)=0.0E+00
            ENDDO
         ENDDO
!
!
         IF (L_LOCAL_CNV_PARTITION) THEN
!
!           PARTITION BETWEEN ICE AND WATER USING THE RELATIONSHIPS
!           GIVEN IN BOWER ET AL. (1996, Q.J. 122 p 1815-1844). ICE
!           IS ALLOWED IN A LAYER WARMER THAN THE FREEZING POINT
!           ONLY IF THE TOP OF THE CLOUD IS GLACIATED.
!

            IF (lrad_cldtop_t_fix) THEN
 
               ! USE MAX INSTEAD OF MIN, TO GET THE CORRECT CLOUD TOP
               ! TEMPERATURE AS INTENDED.

               DO L=1, N_PROFILE
                 IF ((N_LAYER+2-CCT(I_GATHER(L))) <= N_LAYER) THEN
                   L_GLACIATED_TOP(L) =                                 &
                     (T(L, MAX( N_LAYER+2-CCT(I_GATHER(L))              &
                              , N_LAYER-NCLDS+1)) <  TM)
                 END IF
               ENDDO

            ELSE 

               ! THE ORIGINAL (USING MIN) ALWAYS SELECTS THE TEMPERATURE
               ! AT THE TOP OF ATMOSPHERE

               DO L=1, N_PROFILE
                 L_GLACIATED_TOP(L) =                                   &
                   (T(L, MIN( N_LAYER+2-CCT(I_GATHER(L))                &
                            , N_LAYER-NCLDS+1)) <  TM)
               ENDDO

            END IF ! lrad_cldtop_fix


         ELSE
!
!           PARTITION BETWEEN ICE AND WATER AS DIRECTED BY THE
!           TEMPERATURE IN THE MIDDLE OF THE TOP LAYER OF THE CLOUD.
!           THE PARTITIONING MAY BE PRECALCULATED IN THIS CASE.
!

            IF (lrad_cldtop_t_fix) THEN

               ! USE MAX INSTEAD OF MIN, TO GET THE CORRECT CLOUD TOP
               ! TEMPERATURE AS INTENDED.

               DO L=1, N_PROFILE
                 IF ((N_LAYER+2-CCT(I_GATHER(L))) <= N_LAYER) THEN

                   T_GATHER(L) = T(L, MAX( N_LAYER+2-CCT(I_GATHER(L))   &
                                        , N_LAYER-NCLDS+1))
                 END IF
               ENDDO

            ELSE

               ! THE ORIGINAL (USING MIN) ALWAYS SELECTS THE TEMPERATURE
               ! AT THE TOP OF ATMOSPHERE

               DO L=1, N_PROFILE
                 T_GATHER(L) = T(L, MIN( N_LAYER+2-CCT(I_GATHER(L))     &
                                       , N_LAYER-NCLDS+1))
               ENDDO

            END IF ! (lrad_cld_top_t_fix)


! DEPENDS ON: lsp_focwwil
            CALL LSP_FOCWWIL(T_GATHER, N_PROFILE, LIQ_FRAC_CONV)
!
         ENDIF ! Test on local_cnv_partition
!
!
         DO L=1, N_PROFILE
            TOTAL_MASS(L)=0.0E+00
         ENDDO
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF ( (CCT(LG) >= N_LAYER+2-I).AND.                       &
     &              (CCB(LG) <= N_LAYER+1-I) ) THEN
                  TOTAL_MASS(L)=TOTAL_MASS(L)+D_MASS(L, I)
               ENDIF
            ENDDO
         ENDDO
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF ( (CCT(LG) >= N_LAYER+2-I).AND.                       &
     &              (CCB(LG) <= N_LAYER+1-I) ) THEN
                  IF (L_LOCAL_CNV_PARTITION) THEN
!                    THE PARTITIONING IS RECALCULATED FOR EACH LAYER
!                    OTHERWISE A GENERIC VALUE IS USED.
                     LIQ_FRAC_CONV(L)=MAX(0.0E+00, MIN(1.0E+00          &
     &                  , 1.61E-02*(T(L, I)-TM)+8.9E-01))
!                    Do not allow ice above 0 Celsius unless the top
!                    of the cloud is glaciated and force homogeneous
!                    nucleation at -40 Celsius.
                     IF ( (T(L, I) >  TM).AND.                          &
     &                    (.NOT.L_GLACIATED_TOP(L)) ) THEN
                       LIQ_FRAC_CONV(L)=1.0E+00
                     ELSE IF (T(L, I) <  TM-4.0E+01) THEN
                       LIQ_FRAC_CONV(L)=0.0E+00
                     ENDIF
                  ENDIF
                  CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)      &
     &               =CCCWP(LG)*LIQ_FRAC_CONV(L)                        &
     &               /(TOTAL_MASS(L)+TINY(CCCWP))
                  CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)        &
     &               =CCCWP(LG)*(1.0E+00-LIQ_FRAC_CONV(L))              &
     &               /(TOTAL_MASS(L)+TINY(CCCWP))
                  CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)      &
     &               =CCCWP(LG)*LIQ_FRAC_CONV(L)                        &
     &               /(TOTAL_MASS(L)+TINY(CCCWP))
                  CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)        &
     &               =CCCWP(LG)*(1.0E+00-LIQ_FRAC_CONV(L))              &
     &               /(TOTAL_MASS(L)+TINY(CCCWP))
               ENDIF
            ENDDO
         ENDDO
!
!
!        STRATIFORM CLOUDS:
!
!        PARTITION BETWEEN ICE AND WATER DEPENDING ON THE
!        LOCAL TEMPERATURE.
!
         DO I=1, NCLDS
! DEPENDS ON: lsp_focwwil
            CALL LSP_FOCWWIL(T(L, N_LAYER+1-I), N_PROFILE, LIQ_FRAC)
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF (LCA_AREA(LG, I) >  EPSILON(LCA_AREA)) THEN
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_WATER)                                &
     &               =(LCCWC1(LG, I)+LCCWC2(LG, I))                     &
     &               *LIQ_FRAC(L)/LCA_AREA(LG, I)
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_ICE)                                  &
     &               =(LCCWC1(LG, I)+LCCWC2(LG, I))                     &
     &               *(1.0E+00-LIQ_FRAC(L))/LCA_AREA(LG, I)
               ELSE
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_WATER)=0.0E+00
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_ICE)=0.0E+00
               ENDIF
!
               IF (LCA_BULK(LG, I) >  EPSILON(LCA_BULK)) THEN
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_WATER)                                &
     &               =(LCCWC1(LG, I)+LCCWC2(LG, I))                     &
     &               *LIQ_FRAC(L)/LCA_BULK(LG, I)
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_ICE)                                  &
     &               =(LCCWC1(LG, I)+LCCWC2(LG, I))                     &
     &               *(1.0E+00-LIQ_FRAC(L))/LCA_BULK(LG, I)
               ELSE
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_WATER)=0.0E+00
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_ICE)=0.0E+00
               ENDIF
            ENDDO
         ENDDO
!
!
!        CLOUD FRACTIONS:
!
       IF (LCV_3D_CCA) THEN
         DO I=1, NCLDS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               W_CLOUD(L, N_LAYER+1-I)                                  &
     &            =CCA(LG,I)+(1.0E+00-CCA(LG,I))*LCA_AREA(LG, I)
               FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_CONV)           &
     &            =CCA(LG,I)/(W_CLOUD(L, N_LAYER+1-I)+TINY(CCA))
               FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_STRAT)          &
     &            =1.0E+00-FRAC_CLOUD(L, N_LAYER+1-I                    &
     &            , IP_CLOUD_TYPE_CONV)
            ENDDO
         ENDDO
       ELSE
         DO I=1, NCLDS
            DO L=1, N_PROFILE
              LG=I_GATHER(L)
               IF ( (I <= CCT(LG)-1).AND.(I >= CCB(LG)) ) THEN
                  W_CLOUD(L, N_LAYER+1-I)                               &
     &               =CCA(LG,1)+(1.0E+00-CCA(LG,1))*LCA_AREA(LG, I)
                  FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_CONV)        &
     &               =CCA(LG,1)/(W_CLOUD(L, N_LAYER+1-I)+TINY(CCA))
               ELSE
                  W_CLOUD(L, N_LAYER+1-I)=LCA_AREA(LG, I)
                  FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_CONV)        &
     &               =0.0E+00
               ENDIF
               FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_STRAT)          &
     &            =1.0E+00-FRAC_CLOUD(L, N_LAYER+1-I                    &
     &            , IP_CLOUD_TYPE_CONV)
            ENDDO
         ENDDO
       ENDIF
!
!
!
!
      ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
!
!        HERE THE CLOUDS ARE SPLIT INTO FOUR SEPARATE TYPES.
!        THE PARTITIONING BETWEEN ICE AND WATER IS REGARDED AS
!        DETERMINING THE AREAS WITHIN THE GRID_BOX COVERED BY
!        ICE OR WATER CLOUD, RATHER THAN AS DETERMINING THE IN-CLOUD
!        MIXING RATIOS. THE GRID-BOX MEAN ICE WATER CONTENTS IN
!        STRATIFORM CLOUDS MAY BE PREDICTED BY THE ICE MICROPHYSICS
!        SCHEME OR MAY BE DETERMINED AS A FUNCTION OF THE TEMPERATURE
!        (LSP_FOCWWIL). IN CONVECTIVE CLOUDS THE PARTITIONING MAY BE
!        DONE USING THE SAME FUNCTION, LSP_FOCWWIL, BASED ON A SINGLE
!        TEMPERATURE, OR USING A PARTITION BASED ON THE LOCAL
!        TEMPERATURE.
!
!        CONVECTIVE CLOUD:
!

        IF (lrad_ccrad) THEN         

          ! Initialise Convection mixing ratio arrays on ALL levels, not just
          ! the levels that radiation is using (may affect diagnostics
          ! otherwise)

          DO I=0, npd_layer
            DO L=1, npd_profile  
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER) = 0.0E+00
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)   = 0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER) = 0.0E+00
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)   = 0.0E+00
            END DO      ! L
          END DO      ! I

        ELSE      ! Original Code

          DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
              CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)  = 0.0E+00
              CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)    = 0.0E+00
              CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)  = 0.0E+00
              CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)    = 0.0E+00
            END DO
          END DO
        END IF      ! lrad_ccrad



        IF (lrad_ccrad .AND. lrad_3d_ccw) THEN
          !-------------------------------------------------------------------
          ! USE CCW PROFILE WHICH WAS PASSED IN FROM CONVECTION
          !-------------------------------------------------------------------

          DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
              LG = I_GATHER(L)

              CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)                &
                                                         = CCW(LG,N_LAYER+1-I)
              CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)                &
                                                         = CCW(LG,N_LAYER+1-I)

              CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)                  &
                            = CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)
              CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)                  &
                            = CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)

            END DO      ! L (N_LAYER)
          END DO      ! I (N_PROFILE)

        ELSE ! Original Code, use cclwp spread from cloud bottom to top

          DO L=1, N_PROFILE
            TOTAL_MASS(L) = 0.0E+00
          END DO


          DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
              LG = I_GATHER(L)
              IF ( (CCT(LG) >= N_LAYER+2-I).AND.                              &
                   (CCB(LG) <= N_LAYER+1-I) ) THEN
                TOTAL_MASS(L) = TOTAL_MASS(L) + D_MASS(L, I)
              END IF
            END DO      ! L
          END DO      ! I


          DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
              LG = I_GATHER(L)

              IF ( (CCT(LG) >= N_LAYER+2-I).AND.                              &
                   (CCB(LG) <= N_LAYER+1-I) ) THEN

                CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)              &
                                       = CCCWP(LG)/(TOTAL_MASS(L)+TINY(CCCWP))
                CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)              &
                                       = CCCWP(LG)/(TOTAL_MASS(L)+TINY(CCCWP))
       
                CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)                &
                             =CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)
                CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)                &
                             =CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)

              END IF

            END DO      ! L
          END DO      ! I
        END IF      ! lrad_3d_ccw .AND. lrad_ccrad
!
!        STRATIFORM CLOUDS:
!
         DO I=1, NCLDS
            DO L=1, N_PROFILE
               LG=I_GATHER(L)
               IF (LCA_AREA(LG, I) >  EPSILON(LCA_AREA)) THEN
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_WATER)                                &
     &              =(LCCWC1(LG, I)+LCCWC2(LG, I))/LCA_AREA(LG, I)
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_ICE)                                  &
     &              =CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I              &
     &              , IP_CLCMP_ST_WATER)
               ELSE
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_WATER)=0.0E+00
                 CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_ICE)=0.0E+00
               ENDIF
!
               IF (LCA_BULK(LG, I) >  EPSILON(LCA_BULK)) THEN
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_WATER)                                &
     &              =(LCCWC1(LG, I)+LCCWC2(LG, I))/LCA_BULK(LG, I)
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_ICE)                                  &
     &              =CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I              &
     &              , IP_CLCMP_ST_WATER)
               ELSE
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_WATER)=0.0E+00
                 CONDENSED_MIX_RAT_BULK(L, N_LAYER+1-I                  &
     &              , IP_CLCMP_ST_ICE)=0.0E+00
               ENDIF
            ENDDO
         ENDDO


!
!        CLOUD FRACTIONS:
!
         IF (L_LOCAL_CNV_PARTITION) THEN
!
!           PARTITION BETWEEN ICE AND WATER USING THE RELATIONSHIPS
!           GIVEN IN BOWER ET AL. (1996, Q.J. 122 p 1815-1844). ICE
!           IS ALLOWED IN A LAYER WARMER THAN THE FREEZING POINT
!           ONLY IF THE TOP OF THE CLOUD IS GLACIATED.
!
            IF (lrad_cldtop_t_fix) THEN

               ! USE MAX INSTEAD OF MIN, TO GET THE CORRECT CLOUD
               ! TOP TEMPERATURE AS INTENDED.

               DO L=1, N_PROFILE
                 IF ((N_LAYER+2-CCT(I_GATHER(L))) <= N_LAYER) THEN
                   L_GLACIATED_TOP(L) =                                 &
                     (T(L, MAX(N_LAYER+2-CCT(I_GATHER(L))               &
                             , N_LAYER-NCLDS+1)) <  TM)
                 END IF
               END DO

            ELSE

               ! THE ORIGINAL (USING MIN) ALWAYS SELECTS THE TEMPERATURE
               ! AT THE TOP OF ATMOSPHERE.

               DO L=1, N_PROFILE
                 L_GLACIATED_TOP(L) =                                   &
                   (T(L, MIN(N_LAYER+2-CCT(I_GATHER(L))                 &
                           , N_LAYER-NCLDS+1)) <  TM)
               END DO
   
            END IF ! Test on lrad_cldtop_t_fix 

         ELSE
!
!           PARTITION BETWEEN ICE AND WATER AS DIRECTED BY THE
!           TEMPERATURE IN THE MIDDLE OF THE TOP LAYER OF THE CLOUD.
!           THE PARTITIONING MAY BE PRECALCULATED IN THIS CASE.
!
            IF (lrad_cldtop_t_fix) THEN

               ! USE MAX INSTEAD OF MIN, TO GET THE CORRECT CLOUD TOP
               ! TEMPERATURE AS INTENDED.

               DO L=1, N_PROFILE
                 IF ((N_LAYER+2-CCT(I_GATHER(L))) <= N_LAYER) THEN
                   T_GATHER(L) = T(L, MAX(N_LAYER+2-CCT(I_GATHER(L))    &
                                        , N_LAYER-NCLDS+1))
                 END IF
               END DO 

            ELSE

               ! THE ORIGINAL (USING MIN) ALWAYS SELECTS THE TEMPERATURE
               ! AT THE TOP OF ATMOSPHERE

               DO L=1, N_PROFILE
                 T_GATHER(L) = T(L, MIN(N_LAYER+2-CCT(I_GATHER(L))      &
                                      , N_LAYER-NCLDS+1))
               END DO 

            END IF ! Test on lrad_cldtop_t_fix

! DEPENDS ON: lsp_focwwil
            CALL LSP_FOCWWIL(T_GATHER, N_PROFILE, LIQ_FRAC_CONV)
!
         ENDIF



         IF (lrad_ccrad .AND. lrad_ovrlap) THEN

           DO I=N_LAYER+1-NCLDS, N_LAYER

             IF (.NOT. L_CLOUD_WATER_PARTITION) THEN
               ! Partition stratiform clouds using the local temperature.
! DEPENDS ON: lsp_focwwil
               CALL LSP_FOCWWIL(T(1, I), N_PROFILE, LIQ_FRAC)
             END IF


             !----------------------------------------------------------------
             ! Calculate cloud fractions of respective cloud types
             !----------------------------------------------------------------
             DO L=1, N_PROFILE

               ! Modify Convective/Large-scale cloud fractions to account. For
               ! overlap, CCA overlaps LCA and takes precedence, i.e. LCA must
               ! be > CCA in order to have a non-zero cloud fraction.

               LG                            = I_GATHER(L)
               W_CLOUD(L, I)                 = 0.0
               LCA_OF_GRDBX(LG, N_LAYER+1-I) = 0.0
               CONVECTIVE_CLOUD_LAYER(L)     = CCA(LG,N_LAYER+1-I)

               ! Calculate Total Cloud Cover in gridbox
               W_CLOUD(L, I) = MAX(  CONVECTIVE_CLOUD_LAYER(L)                &
                                   , LCA_AREA(LG, N_LAYER+1-I))


               !--------------------------------------------------------------
               ! Scavenge LS Cloud water by CCW
               !--------------------------------------------------------------
               IF (lrad_ccw_scav) THEN 
    
                 ! If Convective cloud and Large-scale cloud are non-zero and
                 ! occur on the same model level. Recalc. CCW to allow for the
                 ! additional condensate from the large-scale cloud
  
                 IF ((LCA_AREA(LG, N_LAYER+1-I) > 1.0E-10) .AND.              &
                     (CONVECTIVE_CLOUD_LAYER(L) > 1.0E-10)) THEN

                   CONDENSED_MIX_RAT_AREA(L,I,IP_CLCMP_CNV_WATER)             &
                     = CCW(LG,N_LAYER+1-I)                                    &
                     + (  (LCCWC1(LG, I)+LCCWC2(LG, I))                       &
                        * MIN(  LCA_AREA(LG, N_LAYER+1-I)                     &
                              , CONVECTIVE_CLOUD_LAYER(L))                    &
                        / CONVECTIVE_CLOUD_LAYER(L) )

                   CONDENSED_MIX_RAT_AREA(L,I,IP_CLCMP_CNV_ICE)               &
                     = CONDENSED_MIX_RAT_AREA(L,I,IP_CLCMP_CNV_WATER)

                 END IF
               END IF      ! lrad_ccw_scav


               !--------------------------------------------------------------
               ! Recalc LCA_OF_GRDBX if there is more large-scale cloud
               ! fraction than cca
               !--------------------------------------------------------------
               IF (  (LCA_AREA(LG, N_LAYER+1-I) - CONVECTIVE_CLOUD_LAYER(L))  &
                   > 1.0E-10) THEN

                 LCA_OF_GRDBX(LG, N_LAYER+1-I) = LCA_AREA(LG, N_LAYER+1-I)    &
                                               - CONVECTIVE_CLOUD_LAYER(L)
               END IF


               !--------------------------------------------------------------
               ! Partition stratiform clouds using the ratio of cloud water
               ! contents.
               !--------------------------------------------------------------
               IF (L_CLOUD_WATER_PARTITION) THEN
                 IF (( LCA_OF_GRDBX(LG, N_LAYER+1-I) > EPSILON(LCA_OF_GRDBX)) &
                    .AND.                                                     &
                    (( LCCWC1(LG,N_LAYER+1-I) + LCCWC2(LG,N_LAYER+1-I))       &
                      > 0.0)) THEN

                   LIQ_FRAC(L) = LCCWC1(LG, N_LAYER+1-I)                      &
                               / (  LCCWC1(LG, N_LAYER+1-I)                   &
                                  + LCCWC2(LG, N_LAYER+1-I))
                 ELSE
                   LIQ_FRAC(L) = 0.0E+00
                 END IF
               END IF      ! L_CLOUD_WATER_PARTITION


               !--------------------------------------------------------------
               ! Partition Convective liquid/ice cloud fraction based on local
               ! temperature.
               !--------------------------------------------------------------
               IF (L_LOCAL_CNV_PARTITION) THEN

                 ! NOTE: Ice is allowed above the freezing point ONLY if the
                 !       TOP of the cloud is glaciated.

                 LIQ_FRAC_CONV(L) = MAX(  0.0E+00                             &
                                        , MIN(  1.0E+00                       & 
                                              , 1.61E-02*(T(L, I)-TM)         &
                                                + 8.9E-01))
   
                 ! Do not allow ice above 0 Celsius unless the top of the
                 ! cloud is glaciated and force homogeneous nucleation at -40
                 ! Celsius

                 IF ( (T(L, I) > TM) .AND.                                    &
                      (.NOT.L_GLACIATED_TOP(L)) ) THEN
                   LIQ_FRAC_CONV(L) = 1.0E+00
                 ELSE IF (T(L, I) < TM-4.0E+01) THEN
                   LIQ_FRAC_CONV(L) = 0.0E+00
                 END IF

               END IF      ! L_LOCAL_CNV_PARTITION


               !--------------------------------------------------------------
               ! Split cloudly fraction of gridbox into ls/convective
               ! liquid/ice cloud fractions. i.e. Fraction of (fraction of
               ! gridbox).
               !--------------------------------------------------------------

               FRAC_CLOUD(L,I,IP_CLOUD_TYPE_SW)                               &
                                = (LIQ_FRAC(L) * LCA_OF_GRDBX(LG,N_LAYER+1-I))&
                                / (W_CLOUD(L,I) + TINY(W_CLOUD))

               FRAC_CLOUD(L,I,IP_CLOUD_TYPE_SI)                               &
                                = ((1.0E+00-LIQ_FRAC(L))                      &
                                * LCA_OF_GRDBX(LG, N_LAYER+1-I))              &
                                / (W_CLOUD(L, I) + TINY(W_CLOUD))

               IF (.NOT. L_LOCAL_CNV_PARTITION         &
                   .AND. lrad_cldtop_t_fix             &
                   .AND. (T(L,I) >= TM)) THEN
                   !-------------------------------------------------------- 
                   ! Convective cloud gridbox fraction (liquid) 
                   !-------------------------------------------------------- 
                   FRAC_CLOUD(L,I,IP_CLOUD_TYPE_CW)                         & 
                               = CONVECTIVE_CLOUD_LAYER(L)                  & 
                               / (W_CLOUD(L,I) + TINY(W_CLOUD)) 
 
                   !-------------------------------------------------------- 
                   ! Convective cloud gridbox fraction (ice) 
                   !-------------------------------------------------------- 
                   FRAC_CLOUD(L,I,IP_CLOUD_TYPE_CI) = 0.0 

               ELSE
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)                             &
                                = LIQ_FRAC_CONV(L) * CONVECTIVE_CLOUD_LAYER(L)&
                                / (W_CLOUD(L, I) + TINY(W_CLOUD))

               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CI)                             &
                                = (1.0E+00-LIQ_FRAC_CONV(L))                  &
                                * CONVECTIVE_CLOUD_LAYER(L)                   &
                                / (W_CLOUD(L, I) + TINY(W_CLOUD))
               END IF

             END DO      ! L (N_PROFILE)
           END DO      ! I (N_LAYER)

         ELSE     ! Original Code

           DO I=N_LAYER+1-NCLDS, N_LAYER

             IF (.NOT. L_CLOUD_WATER_PARTITION) THEN
               ! Partition stratiform clouds using the local temperature.

! DEPENDS ON: lsp_focwwil
               CALL LSP_FOCWWIL(T(1, I), N_PROFILE, LIQ_FRAC)
             END IF

             IF (LCV_3D_CCA) THEN
               DO L=1, N_PROFILE
                 LG = I_GATHER(L)
                 CONVECTIVE_CLOUD_LAYER(L)=CCA(LG,N_LAYER+1-I)
               END DO
             ELSE
               DO L=1, N_PROFILE
                 LG = I_GATHER(L)
                 IF ( (CCT(LG) >= N_LAYER+2-I).AND.                           &
                      (CCB(LG) <= N_LAYER+1-I) ) THEN
                   CONVECTIVE_CLOUD_LAYER(L)=CCA(LG,1)
                 ELSE
                   CONVECTIVE_CLOUD_LAYER(L)=0.0E+00
                 END IF
               END DO
             END IF


             DO L=1, N_PROFILE
               LG = I_GATHER(L)

               W_CLOUD(L, I)                                                  &
                 = CONVECTIVE_CLOUD_LAYER(L)                                  &
                 + (1.0E+00 - CONVECTIVE_CLOUD_LAYER(L))                      &
                 * LCA_AREA(LG, N_LAYER+1-I)



               IF (L_CLOUD_WATER_PARTITION) THEN

                 ! Partition stratiform clouds using the ratio of cloud water
                 ! contents.

                 IF (LCA_AREA(LG, N_LAYER+1-I) > EPSILON(LCA_AREA)            &
                    .AND.                                                     &
                    (LCCWC1(LG, N_LAYER+1-I) + LCCWC2(LG, N_LAYER+1-I))       &
                     > 0.0) THEN

                   LIQ_FRAC(L) =  LCCWC1(LG, N_LAYER+1-I)                     &
                               / (LCCWC1(LG, N_LAYER+1-I)                     &
                               +  LCCWC2(LG, N_LAYER+1-I))
                 ELSE
                   LIQ_FRAC(L) = 0.0E+00
                 END IF
               END IF



               IF (L_LOCAL_CNV_PARTITION) THEN

                 ! The partitioning between ice and water must be recalculated
                 ! for this layer as a function of the local temperature, but
                 ! ice is allowed above the freezing point only if the top of
                 ! the cloud is glaciated.

                 LIQ_FRAC_CONV(L)=MAX(  0.0E+00                               &
                                      , MIN(  1.0E+00                         &
                                            , 1.61E-02*(T(L, I)-TM)+8.9E-01))

                 ! Do not allow ice above 0 Celsius unless the top of the
                 ! cloud is glaciated and force homogeneous nucleation at -40
                 ! Celsius.

                 IF ( (T(L, I) > TM) .AND. (.NOT. L_GLACIATED_TOP(L)) ) THEN
                   LIQ_FRAC_CONV(L) = 1.0E+00
                 ELSE IF (T(L, I) <  TM-4.0E+01) THEN
                   LIQ_FRAC_CONV(L) = 0.0E+00
                 END IF

               END IF
             
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW)                             &
                         = LIQ_FRAC(L) * (1.0E+00-CONVECTIVE_CLOUD_LAYER(L))  &
                         * LCA_AREA(LG, N_LAYER+1-I)                          &
                         / (W_CLOUD(L, I) + TINY(W_CLOUD))

               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SI)                             &
                         = (1.0E+00 - LIQ_FRAC(L))                            &
                         * (1.0E+00 - CONVECTIVE_CLOUD_LAYER(L))              &
                         * LCA_AREA(LG, N_LAYER+1-I)                          &
                         / (W_CLOUD(L, I) + TINY(W_CLOUD))



               IF (.NOT. L_LOCAL_CNV_PARTITION                                &
                   .AND. lrad_cldtop_t_fix                                    &
                   .AND. (T(L,I) >= TM)) THEN
                   !-------------------------------------------------------- 
                   ! Convective cloud gridbox fraction (liquid) 
                   !-------------------------------------------------------- 
                   FRAC_CLOUD(L,I,IP_CLOUD_TYPE_CW)                           & 
                               = CONVECTIVE_CLOUD_LAYER(L)                    & 
                               / (W_CLOUD(L,I) + TINY(W_CLOUD)) 
 
                   !-------------------------------------------------------- 
                   ! Convective cloud gridbox fraction (ice) 
                   !-------------------------------------------------------- 
                   FRAC_CLOUD(L,I,IP_CLOUD_TYPE_CI) = 0.0 

               ELSE


               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)                             &
                         = LIQ_FRAC_CONV(L) * CONVECTIVE_CLOUD_LAYER(L)       &
                         / (W_CLOUD(L, I) + TINY(W_CLOUD)) 

               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CI)                             &
                         = (1.0E+00 - LIQ_FRAC_CONV(L))                       &
                         * CONVECTIVE_CLOUD_LAYER(L)                          &
                         / (W_CLOUD(L, I) + TINY(W_CLOUD))
               END IF
             END DO      ! L (N_PROFILE)
           END DO      ! I (N_LAYER)
         END IF      ! lrad_ccrad .AND. l_overlap



         !--------------------------------------------------------------------
         ! Cloud fraction/Mixing Ratio consistency check
         !--------------------------------------------------------------------
         IF (lrad_ccrad) THEN

           DO I=N_LAYER+1-NCLDS, N_LAYER
             DO L=1, N_PROFILE
               LG = I_GATHER(L)
               IF (FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW) <= 1.0E-10)             &
                        CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_WATER)  = 0.0
               IF (FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SI) <= 1.0E-10)             &
                        CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_ICE)    = 0.0
               IF (FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW) <= 1.0E-10)             &
                        CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER) = 0.0
               IF (FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CI) <= 1.0E-10)             &
                        CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)   = 0.0
             END DO      ! L (N_PROFILE)
           END DO      ! I (N_LAYER)


           ! For consistency with original code, bulk ratios
           ! are set to be the same as the area mixing ratios
           DO I=N_LAYER+1-NCLDS, N_LAYER
             DO L=1, N_PROFILE
               LG = I_GATHER(L)
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_ST_WATER)                &
                      = CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_WATER)
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_ST_ICE)                  &
                      = CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_ICE)
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)               &
                      = CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)                 &
                      = CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)
             END DO      ! L (N_PROFILE)
           END DO      ! I (N_LAYER)

         END IF      ! lrad_ccrad
!
! 
!
         IF ((LRAD_TRIPLECLOUDS) .AND. (I_CLOUD==IP_CLOUD_TRIPLE .OR.         &
                       I_CLOUD==IP_CLOUD_PART_CORR_CNV)) THEN
!
!          If the Tripleclouds flag is specified, we apply Tripleclouds.
!          First, copy the current versions of the fields FRAC_CLOUD, 
!          CONDENSED_MIX_RAT_AREA and CONDENSED_MIX_RAT_BULK to temporary 
!          arrays.
!
           FRAC_CLOUD_TMP = FRAC_CLOUD
           CONDENSED_MIX_RAT_AREA_TMP = CONDENSED_MIX_RAT_AREA
           CONDENSED_MIX_RAT_BULK_TMP = CONDENSED_MIX_RAT_BULK
!
!          Next, cycle through the profiles and layers and apply Tripleclouds 
!          to the mixing ratio quantities and cloud fractions.
!
           DO I=N_LAYER+1-NCLDS, N_LAYER
             DO L=1, N_PROFILE
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_WATER)                &
                      = CONDENSED_MIX_RAT_AREA_TMP(L, I, IP_CLCMP_ST_WATER)   &
                        *(1.00E+00-FW_STD)
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_WATER)               &
                      = CONDENSED_MIX_RAT_AREA_TMP(L, I, IP_CLCMP_ST_WATER)   &
                        *(1.00E+00+FW_STD)
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_ST_ICE)                  &
                      = CONDENSED_MIX_RAT_AREA_TMP(L, I, IP_CLCMP_ST_ICE)     &
                        *(1.00E+00-FW_STD)
               CONDENSED_MIX_RAT_AREA(L, I, IP_CLCMP_CNV_ICE)                 &
                      = CONDENSED_MIX_RAT_AREA_TMP(L, I, IP_CLCMP_ST_ICE)     &
                        *(1.00E+00+FW_STD)
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_ST_WATER)                &
                      = CONDENSED_MIX_RAT_BULK_TMP(L, I, IP_CLCMP_ST_WATER)   &
                        *(1.00E+00-FW_STD)
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_WATER)               &
                      = CONDENSED_MIX_RAT_BULK_TMP(L, I, IP_CLCMP_ST_WATER)   &
                        *(1.00E+00+FW_STD)
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_ST_ICE)                  &
                      = CONDENSED_MIX_RAT_BULK_TMP(L, I, IP_CLCMP_ST_ICE)     &
                        *(1.00E+00-FW_STD)
               CONDENSED_MIX_RAT_BULK(L, I, IP_CLCMP_CNV_ICE)                 &
                      = CONDENSED_MIX_RAT_BULK_TMP(L, I, IP_CLCMP_ST_ICE)     &
                        *(1.00E+00+FW_STD)
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW)                             &
                      = (FRAC_CLOUD_TMP(L, I, IP_CLOUD_TYPE_SW)               &
                         + FRAC_CLOUD_TMP(L, I, IP_CLOUD_TYPE_CW)) / 2.0 
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)                             &
                      = (FRAC_CLOUD_TMP(L, I, IP_CLOUD_TYPE_SW)               &
                         + FRAC_CLOUD_TMP(L, I, IP_CLOUD_TYPE_CW)) / 2.0  
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SI)                             &
                      = (FRAC_CLOUD_TMP(L, I, IP_CLOUD_TYPE_SI)               &
                         + FRAC_CLOUD_TMP(L, I, IP_CLOUD_TYPE_CI)) / 2.0
               FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CI)                             &
                      = (FRAC_CLOUD_TMP(L, I, IP_CLOUD_TYPE_SI)               &
                         + FRAC_CLOUD_TMP(L, I, IP_CLOUD_TYPE_CI)) / 2.0
            END DO 
           END DO
!
        END IF     ! LRAD_TRIPLECLOUDS
!
!
      END IF      ! I_CLOUD_REPRESENTATION


!
!
!
!     EFFECTIVE RADII OF WATER CLOUDS: A MICROPHYSICAL PARAMETRIZATION
!     IS AVAILABLE; OTHERWISE STANDARD VALUES ARE USED.
!
      IF (L_MICROPHYSICS) THEN
!
!        STANDARD VALUES ARE USED FOR ICE CRYSTALS, BUT
!        A PARAMETRIZATION PROVIDED BY UMIST AND MRF
!        IS USED FOR DROPLETS.
!
!        CALCULATE THE DENSITY OF AIR.
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               DENSITY_AIR(L, I)=P(L, I)/(R*T(L, I))
            ENDDO
         ENDDO


         IF (.NOT. lrad_ccrad) THEN
           ! Original code. With lrad_ccrad = T this has been moved into
           ! R2_RE_MRF_UMIST

           Do L=1, N_PROFILE
             CC_DEPTH(L) = 0.0E+00
           END DO

           DO L=1, N_PROFILE
             LG = I_GATHER(L)

             ! This loop should be safe even when convective cloud is not
             ! present, since CCB should not exceed CCT.
             DO I=N_LAYER+2-CCT(LG), N_LAYER+1-CCB(LG)
               CC_DEPTH(L) = CC_DEPTH(L) + (D_MASS(L, I)/DENSITY_AIR(L, I))
             END DO
           END DO
         END IF      ! lrad_ccrad


         DO L=1, N_PROFILE
            LYING_SNOW_G(L)=LYING_SNOW(I_GATHER(L))
         ENDDO
!
         IF (NC_DIAG_FLAG) THEN
            DO I=N_LAYER+1-NCLDS, N_LAYER
               DO L=1, N_PROFILE
                  IF (T(L, I)  >=  TM) THEN
                     STRAT_LIQ_CLOUD_FRACTION(L, I)=W_CLOUD(L, I)       &
     &                           *FRAC_CLOUD(L, I, IP_CLOUD_TYPE_SW)
                     CONV_LIQ_CLOUD_FRACTION(L, I)=W_CLOUD(L, I)        &
     &                           *FRAC_CLOUD(L, I, IP_CLOUD_TYPE_CW)
                  ELSE
                     STRAT_LIQ_CLOUD_FRACTION(L, I)=0.0
                     CONV_LIQ_CLOUD_FRACTION(L, I)=0.0
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
!
! DEPENDS ON: r2_re_mrf_umist
         CALL R2_RE_MRF_UMIST(N_PROFILE, N_LAYER, NCLDS                 &
     &      , I_GATHER                                                  &
     &      , L_PC2, L_AEROSOL_CCN                                      &
     &      , L_BIOMASS_CCN, L_OCFF_CCN                                 &
     &      , SEA_SALT_FILM, SEA_SALT_JET                               &
     &      , L_SEASALT_CCN, SALT_DIM_A, SALT_DIM_B                     &
     &      , L_USE_BIOGENIC, BIOGENIC, BIOGENIC_DIM1, BIOGENIC_DIM2    &
     &      , ACCUM_SULPHATE, DISS_SULPHATE, AITKEN_SULPHATE            &
     &      , AGED_BMASS, CLOUD_BMASS                                   &
     &      , AGED_OCFF, CLOUD_OCFF                                     &
     &      , LYING_SNOW_G, I_CLOUD_REPRESENTATION                      &
     &      , LAND_G, FLANDG_G, DENSITY_AIR                             &
     &      , CONDENSED_MIX_RAT_BULK, CC_DEPTH                          &
     &      , CONDENSED_DIM_CHAR                                        &
     &      , D_MASS                                                    &
     &      , STRAT_LIQ_CLOUD_FRACTION                                  &
     &      , CONV_LIQ_CLOUD_FRACTION                                   &
     &      , NC_DIAG_FLAG                                              &
     &      , NC_DIAG_G                                                 &
     &      , NC_WEIGHT_G                                               &
     &      , NTOT_DIAG_G, Ntot_land, Ntot_sea                          &
     &      , STRAT_LWC_DIAG_G                                          &
     &      , SO4_CCN_DIAG_G                                            &
     &      , SULP_DIM1, SULP_DIM2                                      &
     &      , BMASS_DIM1, BMASS_DIM2                                    &
     &      , OCFF_DIM1, OCFF_DIM2                                      &
     &      , NPD_FIELD, NPD_PROFILE, NPD_LAYER, NPD_AEROSOL_SPECIES    &
     &      )
!
!        CONSTRAIN THE SIZES OF DROPLETS TO LIE WITHIN THE RANGE OF
!        VALIDITY OF THE PARAMETRIZATION SCHEME.
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_WATER)              &
     &            =MAX(CONDENSED_MIN_DIM(IP_CLCMP_ST_WATER)             &
     &            , MIN(CONDENSED_MAX_DIM(IP_CLCMP_ST_WATER)            &
     &            , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_WATER)))
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_WATER)             &
     &            =MAX(CONDENSED_MIN_DIM(IP_CLCMP_CNV_WATER)            &
     &            , MIN(CONDENSED_MAX_DIM(IP_CLCMP_CNV_WATER)           &
     &            , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_WATER)))
            ENDDO
         ENDDO
!
!
!        SET MICROPHYSICAL DIAGNOSTICS. WEIGHTS FOR CLOUD CALCULATED
!        HERE ARE USED SOLELY FOR THE MICROPHYSICS AND DO NOT HAVE
!        AN INDEPENDENT MEANING.
!
         IF (WGT_CONV_FLAG) THEN
            IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     WGT_CONV(col_list(l), row_list(l), I)              &
     &                  =W_CLOUD(L, N_LAYER+1-I)                        &
     &                  *FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_CONV)
                  ENDDO
               ENDDO
            ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     WGT_CONV(col_list(l), row_list(l), I)              &
     &                  =W_CLOUD(L, N_LAYER+1-I)                        &
     &                  *FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_CW)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
!
         IF (RE_CONV_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
!                 EFFECTIVE RADII ARE GIVEN IN MICRONS.
                  RE_CONV(col_list(l), row_list(l), I)                  &
     &               =CONDENSED_DIM_CHAR(L, N_LAYER+1-I                 &
     &               , IP_CLCMP_CNV_WATER)                              &
     &               *WGT_CONV(col_list(l), row_list(l), I)*1.0E+06
               ENDDO
            ENDDO
         ENDIF
!
         IF (WGT_STRAT_FLAG) THEN
            IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CONV_STRAT) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     WGT_STRAT(col_list(l), row_list(l), I)             &
     &                  =W_CLOUD(L, N_LAYER+1-I)                        &
     &                  *FRAC_CLOUD(L, N_LAYER+1-I                      &
     &                  , IP_CLOUD_TYPE_STRAT)
                  ENDDO
               ENDDO
            ELSE IF (I_CLOUD_REPRESENTATION == IP_CLOUD_CSIW) THEN
               DO I=1, NCLDS
                  DO L=1, N_PROFILE
                     WGT_STRAT(col_list(l), row_list(l), I)             &
     &                  =W_CLOUD(L, N_LAYER+1-I)                        &
     &                  *FRAC_CLOUD(L, N_LAYER+1-I, IP_CLOUD_TYPE_SW)
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
!
         IF (RE_STRAT_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
!                 EFFECTIVE RADII ARE GIVEN IN MICRONS.
                  RE_STRAT(col_list(l), row_list(l), I)                 &
     &               =CONDENSED_DIM_CHAR(L, N_LAYER+1-I                 &
     &               , IP_CLCMP_ST_WATER)                               &
     &               *WGT_STRAT(col_list(l), row_list(l), I)*1.0E+06
               ENDDO
            ENDDO
         ENDIF

         IF (LWP_STRAT_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  LWP_STRAT(col_list(l), row_list(l), I)                &
     &               =CONDENSED_MIX_RAT_AREA(L, N_LAYER+1-I             &
     &               , IP_CLCMP_ST_WATER)*D_MASS(L, N_LAYER+1-I)        &
     &               *WGT_STRAT(col_list(l), row_list(l), I)
               ENDDO
            ENDDO
         ENDIF

         IF (NC_DIAG_FLAG .AND. NC_WEIGHT_FLAG) THEN
            DO L=1, N_PROFILE
               NC_DIAG(col_list(L), row_list(L))                        &
     &            =NC_DIAG_G(L)*NC_WEIGHT_G(L)
            ENDDO
         ENDIF

         IF (NC_WEIGHT_FLAG) THEN
            DO L=1, N_PROFILE
               NC_WEIGHT(col_list(L), row_list(L))=NC_WEIGHT_G(L)
            ENDDO
         ENDIF

         IF (NTOT_DIAG_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  NTOT_DIAG(col_list(l), row_list(l), I)                &
     &               =NTOT_DIAG_G(L, N_LAYER+1-I)                       &
     &               *WGT_STRAT(col_list(l), row_list(l), I)
               ENDDO
            ENDDO
         ENDIF

         IF (STRAT_LWC_DIAG_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  STRAT_LWC_DIAG(col_list(l), row_list(l), I)           &
     &               =STRAT_LWC_DIAG_G(L, N_LAYER+1-I)                  &
     &               *WGT_STRAT(col_list(l), row_list(l), I)
               ENDDO
            ENDDO
         ENDIF

! Non-cloud diagnostics are "weighted" by the conditional sampling
! weight COND_SAMP_WGT, but as this is 1.0 if the SW radiation is
! active, and 0.0 if it is not, there is no need to actually
! multiply by it.

         IF (COND_SAMP_WGT_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  COND_SAMP_WGT(col_list(l), row_list(l), I)=1.0
               ENDDO
            ENDDO
         ENDIF

         IF (SO4_CCN_DIAG_FLAG) THEN
            DO I=1, NCLDS
               DO L=1, N_PROFILE
                  SO4_CCN_DIAG(col_list(l), row_list(l), I)             &
     &                    =SO4_CCN_DIAG_G(L, N_LAYER+1-I)
               ENDDO
            ENDDO
         ENDIF
!
!
      ELSE
!
!        ALL EFFECTIVE RADII ARE SET TO STANDARD VALUES.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_WATER)=7.E-6
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_WATER)=7.E-6
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
!     SET THE CHARACTERISTIC DIMENSIONS OF ICE CRYSTALS:
!
!     ICE CRYSTALS IN STRATIFORM CLOUDS:
!
      IF (I_CONDENSED_PARAM(IP_CLCMP_ST_ICE) ==                         &
     &   IP_SLINGO_SCHRECKER_ICE) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE EFFECTIVE RADIUS
!        AND A STANDARD VALUE OF 30-MICRONS IS ASSUMED.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)=30.E-6
            ENDDO
         ENDDO
!
      ELSE IF (I_CONDENSED_PARAM(IP_CLCMP_ST_ICE) ==                    &
     &   IP_ICE_ADT) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE MEAN MAXIMUM
!        DIMENSION OF THE CRYSTAL, DETERMINED AS A FUNCTION OF
!        THE LOCAL TEMPERATURE. THE SIZE IS LIMITED TO ITS VALUE
!        AT THE FREEZING LEVEL.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)                &
     &            =MIN(7.198755E-04                                     &
     &            , EXP(5.522E-02*(T(L, I)-2.7965E+02))/9.702E+02)
            ENDDO
         ENDDO
!
      ELSE IF (I_CONDENSED_PARAM(IP_CLCMP_ST_ICE) ==                    &
     &   IP_ICE_AGG_DE) THEN
!
!      Aggregate parametrization based on effective dimension.
!      In the initial form, the same approach is used for stratiform
!      and convective cloud.
!
!      The fit provided here is based on Stephan Havemann's fit of
!      Dge with temperature, consistent with David Mitchell's treatment
!      of the variation of the size distribution with temperature. The
!      parametrization of the optical properties is based on De
!      (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
!      (=(2*SQRT(3)/3)*volume/projected area), which explains the
!      conversion factor. The fit to Dge is in two sections, because
!      Mitchell's relationship predicts a cusp at 216.208 K. Limits
!      of 8 and 124 microns are imposed on Dge: these are based on this
!      relationship and should be reviewed if it is changed. Note also
!      that the relationship given here is for polycrystals only.
       DO I=N_LAYER+1-NCLDS, N_LAYER
         DO L=1, N_PROFILE
!          Preliminary calculation of Dge.
           IF (T(L, I) < t_switch) THEN
             CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)                  &
     &         = a0_agg_cold*EXP(s0_agg*(T(L, I)-t0_agg))+b0_agg_cold
           ELSE
             CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)                  &
     &         = a0_agg_warm*EXP(s0_agg*(T(L, I)-t0_agg))+b0_agg_warm
           ENDIF
!          Limit and convert to De.
           CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)                    &
     &       = (3.0/2.0)*(3.0/(2.0*SQRT(3.0)))*                         &
     &         MIN(1.24E-04, MAX(8.0E-06,                               &
     &         CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)))
         ENDDO
       ENDDO
!
!
      ENDIF
!
!
!     ICE CRYSTALS IN CONVECTIVE CLOUDS:
!
      IF (I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE) ==                        &
     &   IP_SLINGO_SCHRECKER_ICE) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE EFFECTIVE RADIUS
!        AND A STANDARD VALUE OF 30-MICRONS IS ASSUMED.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)=30.E-6
            ENDDO
         ENDDO
!
      ELSE IF (I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE) ==                   &
     &   IP_ICE_ADT) THEN
!
!        THIS PARAMETRIZATION IS BASED ON THE MEAN MAXIMUM
!        DIMENSION OF THE CRYSTAL, DETERMINED AS A FUNCTION OF
!        THE LOCAL TEMPERATURE. THE SIZE IS LIMITED TO ITS VALUE
!        AT THE FREEZING LEVEL.
!
         DO I=N_LAYER+1-NCLDS, N_LAYER
            DO L=1, N_PROFILE
               CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)               &
     &            =MIN(7.198755E-04                                     &
     &            , EXP(5.522E-02*(T(L, I)-2.7965E+02))/9.702E+02)
            ENDDO
         ENDDO
      ELSE IF (I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE) ==                   &
     &   IP_ICE_AGG_DE) THEN
!
!      Aggregate parametrization based on effective dimension.
!      In the initial form, the same approach is used for stratiform
!      and convective cloud.
!
!      The fit provided here is based on Stephan Havemann's fit of
!      Dge with temperature, consistent with David Mitchell's treatment
!      of the variation of the size distribution with temperature. The
!      parametrization of the optical properties is based on De
!      (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
!      (=(2*SQRT(3)/3)*volume/projected area), which explains the
!      conversion factor. The fit to Dge is in two sections, because
!      Mitchell's relationship predicts a cusp at 216.208 K. Limits
!      of 8 and 124 microns are imposed on Dge: these are based on this
!      relationship and should be reviewed if it is changed. Note also
!      that the relationship given here is for polycrystals only.
       DO I=N_LAYER+1-NCLDS, N_LAYER
         DO L=1, N_PROFILE
!          Preliminary calculation of Dge.
           IF (T(L, I) < t_switch) THEN
             CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)                 &
     &         = a0_agg_cold*EXP(s0_agg*(T(L, I)-t0_agg))+b0_agg_cold
           ELSE
             CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)                 &
     &         = a0_agg_warm*EXP(s0_agg*(T(L, I)-t0_agg))+b0_agg_warm
           ENDIF
!          Limit and convert to De.
           CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)                   &
     &       = (3.0/2.0)*(3.0/(2.0*SQRT(3.0)))*                         &
     &         MIN(1.24E-04, MAX(8.0E-06,                               &
     &         CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)))
         ENDDO
       ENDDO
!
      ENDIF
!
!
!
!     CONSTRAIN THE SIZES OF ICE CRYSTALS TO LIE WITHIN THE RANGE
!     OF VALIDITY OF THE PARAMETRIZATION SCHEME.
      DO I=N_LAYER+1-NCLDS, N_LAYER
         DO L=1, N_PROFILE
            CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)                   &
     &         =MAX(CONDENSED_MIN_DIM(IP_CLCMP_ST_ICE)                  &
     &         , MIN(CONDENSED_MAX_DIM(IP_CLCMP_ST_ICE)                 &
     &         , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_ST_ICE)))
            CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)                  &
     &         =MAX(CONDENSED_MIN_DIM(IP_CLCMP_CNV_ICE)                 &
     &         , MIN(CONDENSED_MAX_DIM(IP_CLCMP_CNV_ICE)                &
     &         , CONDENSED_DIM_CHAR(L, I, IP_CLCMP_CNV_ICE)))
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE R2_SET_CLOUD_FIELD

!+ Subroutine to set the parametrization schemes for clouds.
!
! Purpose:
!   The parametrization schemes for each component within a cloud
!   are set.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of aerosols.
!
! Purpose:
!   The mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set fields of climatological aerosols in HADCM3.
!
! Purpose:
!   This routine sets the mixing ratios of climatological aerosols.
!   A separate subroutine is used to ensure that the mixing ratios
!   of these aerosols are bit-comparable with earlier versions of
!   the model where the choice of aerosols was more restricted:
!   keeping the code in its original form reduces the opportunity
!   for optimizations which compromise bit-reproducibilty.
!   The climatoogy used here is the one devised for HADCM3.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to calculate the total cloud cover.
!
! Purpose:
!   The total cloud cover at all grid-points is determined.
!
! Method:
!   A separate calculation is made for each different assumption about
!   the overlap.
!
! Current Owner of Code: J. M. Edwards
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
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
! Current Owner of Code: J. M. Edwards
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to calculate column-integrated cloud droplet number.
!
! Purpose:
!   To calculate a diagnostic of column-integrated cloud droplet
!   number which may be validated aginst satellite data.
!
! Method:
!   Column cloud droplet concentration (i.e. number of droplets per
!   unit area) is calculated as the vertically integrated droplet
!   number concentration averaged over the portion of the gridbox
!   covered by stratiform and convective liquid cloud with T>273K.
!
! Current Owner of Code: A. Jones
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ Subroutine to set the actual process options for the radiation code.
!
! Purpose:
!   To set a consistent set of process options for the radiation.
!
! Method:
!   The global options for the spectral region are compared with the
!   contents of the spectral file. The global options should be set
!   to reflect the capabilities of the code enabled in the model.
!
! Current Owner of Code: J. M. Edwards
!
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
#endif
#endif
